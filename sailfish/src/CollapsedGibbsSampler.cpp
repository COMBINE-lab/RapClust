#include <vector>
#include <unordered_map>
#include <atomic>
#include <random>

#include "tbb/task_scheduler_init.h"
#include "tbb/parallel_for.h"
#include "tbb/parallel_for_each.h"
#include "tbb/parallel_reduce.h"
#include "tbb/blocked_range.h"
#include "tbb/partitioner.h"

#include <boost/filesystem.hpp>

// C++ string formatting library
#include "spdlog/details/format.h"

#include "cuckoohash_map.hh"

#include "CollapsedGibbsSampler.hpp"
#include "Transcript.hpp"
#include "TranscriptGroup.hpp"
#include "SailfishMath.hpp"
#include "ReadExperiment.hpp"
#include "MultinomialSampler.hpp"
#include "BootstrapWriter.hpp"

using BlockedIndexRange =  tbb::blocked_range<size_t>;

// intelligently chosen value adopted from
// https://github.com/pachterlab/kallisto/blob/master/src/EMAlgorithm.h#L18
constexpr double minEQClassWeight = std::numeric_limits<double>::denorm_min();
constexpr double minWeight = std::numeric_limits<double>::denorm_min();

void initCountMap_(
        std::vector<std::pair<const TranscriptGroup, TGValue>>& eqVec,
	std::vector<Transcript>& transcriptsIn,
	double priorAlpha,
        MultinomialSampler& msamp,
        std::vector<uint64_t>& countMap,
        std::vector<double>& probMap,
	std::vector<int>& txpCounts) {

    size_t offset{0};
    for (auto& eqClass : eqVec) {
        uint64_t classCount = eqClass.second.count;

        // for each transcript in this class
        const TranscriptGroup& tgroup = eqClass.first;
        const size_t groupSize = tgroup.txps.size();
        if (tgroup.valid) {
            const std::vector<uint32_t>& txps = tgroup.txps;
            const std::vector<double>& auxs = eqClass.second.weights;

            double denom = 0.0;
            if (BOOST_LIKELY(groupSize > 1)) {

                for (size_t i = 0; i < groupSize; ++i) {
                    auto tid = txps[i];
                    auto aux = auxs[i];
                    denom += (priorAlpha + transcriptsIn[tid].mass(false)) * aux;
                    countMap[offset + i] = 0;
                }

		if (denom > ::minEQClassWeight) {
	   	   // Get the multinomial probabilities
		   double norm = 1.0 / denom;
		   for (size_t i = 0; i < groupSize; ++i) {
		     auto tid = txps[i];
		     auto aux = auxs[i];
		     probMap[offset + i] = norm *
                        ((priorAlpha + transcriptsIn[tid].mass(false)) * aux);
		    }

	   	    // re-sample
	            msamp(countMap.begin() + offset,
                      classCount,
                      groupSize,
                      probMap.begin() + offset);
		}
            } else {
                countMap[offset] = classCount;
            }


            for (size_t i = 0; i < groupSize; ++i) {
                auto tid = txps[i];
                txpCounts[tid] += countMap[offset + i];
            }

            offset += groupSize;
       } // valid group
    } // loop over all eq classes
}

void sampleRound_(
        std::vector<std::pair<const TranscriptGroup, TGValue>>& eqVec,
        std::vector<uint64_t>& countMap,
        std::vector<double>& probMap,
        double priorAlpha,
        std::vector<int>& txpCount,
        MultinomialSampler& msamp) {

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.25, 0.75);
    size_t offset{0};
    // Choose a fraction of this class to re-sample

    // The count substracted from each transcript
    std::vector<uint64_t> txpResamp;

    for (auto& eqClass : eqVec) {
        uint64_t classCount = eqClass.second.count;
        double sampleFrac = dis(gen);

        // for each transcript in this class
        const TranscriptGroup& tgroup = eqClass.first;
        const size_t groupSize = tgroup.txps.size();
        if (tgroup.valid) {
            const std::vector<uint32_t>& txps = tgroup.txps;
            const std::vector<double>& auxs = eqClass.second.weights;

            double denom = 0.0;
            // If this is a single-transcript group,
            // then it gets the full count --- otherwise,
	    // sample!
            if (BOOST_LIKELY(groupSize > 1)) {

                // Subtract some fraction of the current equivalence
                // class' contribution from each transcript.
		uint64_t numResampled{0};
		if (groupSize > txpResamp.size()) {
			txpResamp.resize(groupSize, 0);
		}

		// For each transcript in the group
                for (size_t i = 0; i < groupSize; ++i) {
                    auto tid = txps[i];
                    auto aux = auxs[i];
                    auto currCount = countMap[offset + i];
		    uint64_t currResamp = std::round(sampleFrac * currCount);
		    numResampled += currResamp;
		    txpResamp[i] = currResamp;
                    txpCount[tid] -= currResamp;
                    countMap[offset + i] -= currResamp;
                    denom += (priorAlpha + txpCount[tid]) * aux;
                }

		if (denom > ::minEQClassWeight) {
			// Get the multinomial probabilities
			double norm = 1.0 / denom;
			for (size_t i = 0; i < groupSize; ++i) {
			    auto tid = txps[i];
			    auto aux = auxs[i];
			    probMap[offset + i] = norm * ((priorAlpha + txpCount[tid]) * aux);
			}

			// re-sample
			msamp(txpResamp.begin(),        // count array to fill in
			      numResampled,		// multinomial n
			      groupSize,		// multinomial k
			      probMap.begin() + offset  // where to find multinomial probs
			      );

			for (size_t i = 0; i < groupSize; ++i) {
				auto tid = txps[i];
				countMap[offset + i] += txpResamp[i];
				txpCount[tid] += txpResamp[i];
			}

		} else { // We didn't sample
			// add back to txp-count!
			for (size_t i = 0; i < groupSize; ++i) {
			    auto tid = txps[i];
			    txpCount[tid] += txpResamp[i];
			    countMap[offset + i] += txpResamp[i];
			}
		}
            }

            offset += groupSize;
        } // valid group
    } // loop over all eq classes

}

CollapsedGibbsSampler::CollapsedGibbsSampler() {}

class DistStats {
	public:
	DistStats() : meanVal(0.0), minVal(std::numeric_limits<double>::max()), maxVal(0.0) {}
	double meanVal;
	double minVal;
	double maxVal;
};

template <typename ExpT>
bool CollapsedGibbsSampler::sample(ExpT& readExp,
        SailfishOpts& sopt,
        std::function<bool(const std::vector<int>&)>& writeSample, 
        uint32_t numSamples) {

    namespace bfs = boost::filesystem;
    tbb::task_scheduler_init tbbScheduler(sopt.numThreads);
    std::vector<Transcript>& transcripts = readExp.transcripts();

    std::vector<std::pair<const TranscriptGroup, TGValue>>& eqVec =
        readExp.equivalenceClassBuilder().eqVec();

    using VecT = CollapsedGibbsSampler::VecType;

    std::vector<std::vector<int>> allSamples(numSamples,
                                        std::vector<int>(transcripts.size(),0));
    double priorAlpha = 1e-8;
    auto numMappedFragments = readExp.numMappedFragments();


    for (auto& txp : transcripts) {
        txp.setMass(priorAlpha + (txp.mass(false) * numMappedFragments));
    }

    tbb::parallel_for(BlockedIndexRange(size_t(0), size_t(numSamples)),
                [&eqVec, &transcripts, priorAlpha, &writeSample,
                 &allSamples]( const BlockedIndexRange& range) -> void {

                std::random_device rd;
                MultinomialSampler ms(rd);

                size_t countMapSize{0};
                for (size_t i = 0; i < eqVec.size(); ++i) {
                    if (eqVec[i].first.valid) {
                        countMapSize += eqVec[i].first.txps.size();
                    }
                }

                size_t numTranscripts{transcripts.size()};

                // will hold estimated counts
                std::vector<double> alphas(numTranscripts, 0.0);
                std::vector<uint64_t> countMap(countMapSize, 0);
                std::vector<double> probMap(countMapSize, 0.0);

                initCountMap_(eqVec, transcripts, priorAlpha, ms, countMap, probMap, allSamples[range.begin()]);

                // For each sample this thread should generate
                bool isFirstSample{true};
                bool numInternalRounds = 10;
                for (auto sampleID : boost::irange(range.begin(), range.end())) {
                    if (sampleID % 100 == 0) {
                        std::cerr << "gibbs sampling " << sampleID << "\n";
                    }
                    if (!isFirstSample) {
                        // the counts start at what they were last round.
                        allSamples[sampleID] = allSamples[sampleID-1];
                    }
                    for (size_t i = 0; i < numInternalRounds; ++i){
                        sampleRound_(eqVec, countMap, probMap, priorAlpha,
                                allSamples[sampleID], ms);
                    }
		    /*
                    for (size_t tn = 0; tn < numTranscripts; ++tn) {
                        alphas[tn] = static_cast<double>(allSamples[sampleID][tn]);
                    }
		    */
                    //writeSample(alphas);//bootstrapWriter->writeBootstrap(alphas);
		    writeSample(allSamples[sampleID]);
                    isFirstSample = false;
                }
            });

    /*
    double cutoff = priorAlpha + 1e-8;
    // Truncate tiny expression values
    double txpSumTrunc = 0.0;
    for (size_t i = 0; i < transcripts.size(); ++i) {
	// maybe use the first decile instead of the mean for the cutoff;
	// this could let too much through
        if (transcripts[i].mass(false) <= cutoff) { transcripts[i].setMass(0.0); }
        txpSumTrunc += transcripts[i].mass(false);
    }

    for (size_t i = 0; i < transcripts.size(); ++i) {
        // Set the mass to the normalized (after truncation)
        // relative abundance
        transcripts[i].setMass(transcripts[i].mass(false) / txpSumTrunc);
    }
    */

    return true;
}

template
bool CollapsedGibbsSampler::sample<ReadExperiment>(ReadExperiment& readExp,
        SailfishOpts& sopt,
        std::function<bool(const std::vector<int>&)>& writeSample,
        uint32_t maxIter);
