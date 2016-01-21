#ifndef __SAMPLER__HPP__
#define __SAMPLER__HPP__

extern "C" {
#include "io_lib/scram.h"
#include "io_lib/os.h"
#undef max
#undef min
}

// for cpp-format
#include "format.h"

#include "spdlog/spdlog.h"

#include <tbb/atomic.h>
#include <iostream>
#include <fstream>
#include <atomic>
#include <vector>
#include <random>
#include <memory>
#include <exception>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <mutex>
#include <thread>
#include <condition_variable>

#include <tbb/concurrent_queue.h>

#include <boost/config.hpp>
#include <boost/timer/timer.hpp>
#include <boost/filesystem.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/program_options.hpp>
#include <boost/pending/disjoint_sets.hpp>
#include <boost/math/distributions/normal.hpp>

#include "ClusterForest.hpp"
#include "AlignmentLibrary.hpp"
#include "MiniBatchInfo.hpp"
#include "BAMQueue.hpp"
#include "SailfishMath.hpp"
#include "FASTAParser.hpp"
#include "LibraryFormat.hpp"
#include "Transcript.hpp"
#include "ReadPair.hpp"
#include "ErrorModel.hpp"
#include "AlignmentModel.hpp"
#include "FragmentLengthDistribution.hpp"
#include "TranscriptCluster.hpp"
#include "SailfishUtils.hpp"
#include "SalmonUtils.hpp"
#include "SalmonConfig.hpp"
#include "SalmonOpts.hpp"
#include "OutputUnmappedFilter.hpp"

namespace salmon {
    namespace sampler {

        namespace bfs = boost::filesystem;
        using sailfish::math::LOG_0;
        using sailfish::math::LOG_1;
        using sailfish::math::logAdd;
        using sailfish::math::logSub;

        template <typename FragT>
            using AlignmentBatch = std::vector<FragT>;

        template <typename FragT>
            using MiniBatchQueue = tbb::concurrent_queue<MiniBatchInfo<FragT>*>;

        template <typename FragT>
            using OutputQueue = tbb::concurrent_bounded_queue<FragT*>;

        template <typename FragT>
        void sampleMiniBatch(AlignmentLibrary<FragT>& alnLib,
                    MiniBatchQueue<AlignmentGroup<FragT*>>& workQueue,
                    std::condition_variable& workAvailable,
                    std::mutex& cvmutex,
                    volatile bool& doneParsing,
                    std::atomic<size_t>& activeBatches,
                    const SalmonOpts& salmonOpts,
                    bool& burnedIn,
                    std::atomic<size_t>& processedReads,
                    OutputQueue<FragT>& outputQueue) {

                // Seed with a real random value, if available
                std::random_device rd;
                auto log = spdlog::get("jointLog");

                // Create a random uniform distribution
                std::default_random_engine eng(rd());
                std::uniform_real_distribution<> uni(0.0, 1.0 + std::numeric_limits<double>::min());

                using sailfish::math::LOG_0;
                using sailfish::math::logAdd;
                using sailfish::math::logSub;

                auto& refs = alnLib.transcripts();
                auto& clusterForest = alnLib.clusterForest();
                auto& fragmentQueue = alnLib.fragmentQueue();
                auto& alignmentGroupQueue = alnLib.alignmentGroupQueue();
                auto& fragLengthDist = alnLib.fragmentLengthDistribution();
                auto& errMod = alnLib.alignmentModel();

                const auto expectedLibraryFormat = alnLib.format();

                std::chrono::microseconds sleepTime(1);
                MiniBatchInfo<AlignmentGroup<FragT*>>* miniBatch = nullptr;
                size_t numTranscripts = refs.size();

                while (!doneParsing or !workQueue.empty()) {
                    bool foundWork = workQueue.try_pop(miniBatch);
                    // If work wasn't immediately available, then wait for it.
                    if (!foundWork) {
                        std::unique_lock<std::mutex> l(cvmutex);
                        workAvailable.wait(l, [&miniBatch, &workQueue, &doneParsing]() { return workQueue.try_pop(miniBatch) or doneParsing; });
                    }
                    if (miniBatch != nullptr) {
                       ++activeBatches;
                        size_t batchReads{0};
                        std::vector<AlignmentGroup<FragT*>*>& alignmentGroups = *(miniBatch->alignments);

                        using TranscriptID = size_t;
                        using HitIDVector = std::vector<size_t>;
                        using HitProbVector = std::vector<double>;

                        std::unordered_map<TranscriptID, std::vector<FragT*>> hitList;
                        // Each alignment group corresponds to all of the potential
                        // mapping locations of a multi-mapping read
                        for (auto alnGroup : alignmentGroups) {
                            // Score all of these alignments and sample according to
                            // their probabilities
                            for (auto a : alnGroup->alignments()) {
                                auto transcriptID = a->transcriptID();
                                if (transcriptID < 0 or transcriptID >= refs.size()) {
                                    log->warn("Invalid Transcript ID: {}\n", transcriptID);
                                }
                                hitList[transcriptID].emplace_back(a);
                            }
                        }

                        {
                            // Iterate over each group of alignments (a group consists of all alignments reported
                            // for a single read).
                            for (auto alnGroup : alignmentGroups) {
                                double sumOfAlignProbs{LOG_0};
                                // update the cluster-level properties
                                bool transcriptUnique{true};
                                auto firstTranscriptID = alnGroup->alignments().front()->transcriptID();
                                for (auto& aln : alnGroup->alignments()) {
                                    auto transcriptID = aln->transcriptID();
                                    auto& transcript = refs[transcriptID];
                                    transcriptUnique = transcriptUnique and (transcriptID == firstTranscriptID);

                                    double refLength = transcript.RefLength > 0 ? transcript.RefLength : 1.0;

                                    double logFragProb = sailfish::math::LOG_1;

                                    if (!salmonOpts.noFragLengthDist) {
                                        if(aln->fragLen() == 0) {
                                            if (aln->isLeft() and transcript.RefLength - aln->left() < fragLengthDist.maxVal()) {
                                                logFragProb = fragLengthDist.cmf(transcript.RefLength - aln->left());
                                            } else if (aln->isRight() and aln->right() < fragLengthDist.maxVal()) {
                                                logFragProb = fragLengthDist.cmf(aln->right());
                                            }
                                        } else {
                                            logFragProb = fragLengthDist.pmf(static_cast<size_t>(aln->fragLen()));
                                        }
                                    }

                                    // The alignment probability is the product of a transcript-level term (based on abundance and) an alignment-level
                                    // term below which is P(Q_1) * P(Q_2) * P(F | T)
                                    double logRefLength = std::log(refLength);

                                    double logAlignCompatProb = (salmonOpts.useReadCompat) ?
                                        (salmon::utils::logAlignFormatProb(aln->libFormat(), expectedLibraryFormat, salmonOpts.incompatPrior)) :
                                        LOG_1;

                                    // P(Fn | Tn) = Probability of selecting a fragment of this length, given the transcript is t
                                    // d(Fn) / sum_x = 1^{lt} d(x)
                                    double qualProb = -logRefLength + logFragProb + aln->logQualProb() + logAlignCompatProb;
                                    double transcriptLogCount = transcript.mass(false);

                                    if ( transcriptLogCount != LOG_0 ) {
                                        double errLike = sailfish::math::LOG_1;

                                        if (burnedIn and salmonOpts.useErrorModel) {
                                            errLike = errMod.logLikelihood(*aln, transcript);
                                        }

                                        aln->logProb = transcriptLogCount + qualProb + errLike;

                                        sumOfAlignProbs = logAdd(sumOfAlignProbs, aln->logProb);
                                    } else {
                                        aln->logProb = LOG_0;
                                    }
                                }
                                // normalize the hits
                                if (sumOfAlignProbs == LOG_0) {
                                    auto aln = alnGroup->alignments().front();
                                    log->warn("0 probability fragment [{}] "
                                              "encountered\n", aln->getName());
                                    continue;
                                }


                                if (transcriptUnique) {
                                    // avoid r-value ref until we figure out what's
                                    // up with TBB 4.3
                                    auto* alnPtr = alnGroup->alignments().front()->clone();
                                    outputQueue.push(alnPtr);
                                } else { // read maps to multiple transcripts
                                    double r = uni(eng);
                                    double currentMass{0.0};
                                    double massInc{0.0};
                                    bool choseAlignment{false};
                                    for (auto& aln : alnGroup->alignments()) {
                                        if (aln->logProb == LOG_0) { continue; }
                                        aln->logProb -= sumOfAlignProbs;

                                        massInc = std::exp(aln->logProb);
                                        if (currentMass <= r and currentMass + massInc > r) {
                                            // Write out this read
                                            // avoid r-value ref until we figure out what's
                                            // up with TBB 4.3
                                           auto* alnPtr = aln->clone();
                                            outputQueue.push(alnPtr);
                                            currentMass += massInc;
                                            choseAlignment = true;
                                            break;
                                        }
                                        currentMass += massInc;
                                    } // end alignment group
                                    if (BOOST_UNLIKELY(!choseAlignment)) {
                                        log->warn("[Sampler.hpp]: Failed to sample an alignment for this read; "
                                                  "this shouldn't happen\n"
                                                  "currentMass = {}, r = {}\n", currentMass, r);
                                    }
                                } // non-unique read
                                ++batchReads;
                            } // end read group
                        }// end timer

                        miniBatch->release(fragmentQueue, alignmentGroupQueue);
                        delete miniBatch;
                        --activeBatches;
                        processedReads += batchReads;
                    }
                    miniBatch = nullptr;
                } // nothing left to process
            }

        /**
         *  Sample the alignments in the provided library given in current
         *  estimates of transcript abundance.
         */
        template <typename FragT>
        bool sampleLibrary(
                    AlignmentLibrary<FragT>& alnLib,
                    const SalmonOpts& salmonOpts,
                    bool burnedIn,
                    bfs::path& sampleFilePath,
                    bool sampleUnaligned) {

                fmt::MemoryWriter msgStr;
                auto log = spdlog::get("jointLog");

                msgStr << "Sampling alignments; outputting results to "
                       << sampleFilePath.string() << "\n";

                log->info() << msgStr.str();

                auto& refs = alnLib.transcripts();
                size_t numTranscripts = refs.size();
                size_t miniBatchSize{1000};
                size_t numObservedFragments{0};

                MiniBatchQueue<AlignmentGroup<FragT*>> workQueue;
                double logForgettingMass{std::log(1.0)};
                double forgettingFactor{0.60};
                size_t batchNum{0};

                /**
                * Output queue
                */
                volatile bool consumedAllInput{false};
                size_t defaultCapacity = 2000000;
                OutputQueue<FragT> outQueue;
                outQueue.set_capacity(defaultCapacity);

                std::unique_ptr<OutputUnmappedFilter<FragT>> outFilt = nullptr;

                if (sampleUnaligned) {
                    outFilt.reset(new OutputUnmappedFilter<FragT>(&outQueue));
                }

                // Reset our reader to the beginning
                if (!alnLib.reset(false, outFilt.get())) {
                    fmt::print(stderr,
                            "\n\n======== WARNING ========\n"
                            "A provided alignment file "
                            "is not a regular file and therefore can't be read from "
                            "more than once.\n\n"
                            "Therefore, we cannot provide sample output alignments "
                            "from this file.  No sampled BAM file will be generated. "
                            "Please consider re-running Salmon with these alignments "
                            "as a regular file!\n"
                            "==========================\n\n");

                    return false;
                }

                volatile bool doneParsing{false};
                std::condition_variable workAvailable;
                std::mutex cvmutex;
                std::vector<std::thread> workers;
                std::atomic<size_t> activeBatches{0};
                std::atomic<size_t> processedReads{0};

                size_t numProc{0};
                for (uint32_t i = 0; i < salmonOpts.numQuantThreads; ++i) {
                    workers.emplace_back(sampleMiniBatch<FragT>,
                            std::ref(alnLib),
                            std::ref(workQueue),
                            std::ref(workAvailable), std::ref(cvmutex),
                            std::ref(doneParsing), std::ref(activeBatches),
                            std::ref(salmonOpts),
                            std::ref(burnedIn),
                            std::ref(processedReads),
                            std::ref(outQueue));
                }

                std::thread outputThread(
                        [&consumedAllInput, &alnLib, &outQueue, &log, sampleFilePath] () -> void {

                            scram_fd* bf = scram_open(sampleFilePath.c_str(), "wb");
                            scram_set_option(bf, CRAM_OPT_NTHREADS, 3);
                            scram_set_header(bf, alnLib.header());
                            scram_write_header(bf);
                            if (bf == nullptr) {
                                fmt::MemoryWriter errstr;
                                errstr << ioutils::SET_RED << "ERROR: "
                                       << ioutils::RESET_COLOR
                                       << "Couldn't open output bam file "
                                       << sampleFilePath.string() << ". Exiting\n";
                                log->warn() << errstr.str();
                                std::exit(-1);
                            }

                            FragT* aln{nullptr};
                            while (!outQueue.empty() or !consumedAllInput) {
                                while (outQueue.try_pop(aln)) {
                                    if (aln != nullptr) {
                                        int ret = aln->writeToFile(bf);
                                        if (ret != 0) {
                                            std::cerr << "ret = " << ret << "\n";
                                            fmt::MemoryWriter errstr;
                                            errstr << ioutils::SET_RED << "ERROR:"
                                                << ioutils::RESET_COLOR << "Could not write "
                                                << "a sampled alignment to the output BAM "
                                                << "file. Please check that the file can "
                                                << "be created properly and that the disk "
                                                << "is not full.  Exiting.\n";
                                            log->warn() << errstr.str();
                                            std::exit(-1);
                                        }
                                        // Eventually, as we do in BAMQueue, we should
                                        // have queue of bam1_t structures that can be
                                        // re-used rather than continually calling
                                        // new and delete.
                                        delete aln;
                                        aln = nullptr;
                                    }
                                }
                            }

                            scram_close(bf); // will delete the header itself
                        });



                BAMQueue<FragT>& bq = alnLib.getAlignmentGroupQueue();
                std::vector<AlignmentGroup<FragT*>*>* alignments = new std::vector<AlignmentGroup<FragT*>*>;
                alignments->reserve(miniBatchSize);
                AlignmentGroup<FragT*>* ag;

                bool alignmentGroupsRemain = bq.getAlignmentGroup(ag);
                while (alignmentGroupsRemain or alignments->size() > 0) {
                    if (alignmentGroupsRemain) { alignments->push_back(ag); }
                    // If this minibatch has reached the size limit, or we have nothing
                    // left to fill it up with
                    if (alignments->size() >= miniBatchSize or !alignmentGroupsRemain) {
                        // Don't need to update the batch number or log forgetting mass in this phase
                        MiniBatchInfo<AlignmentGroup<FragT*>>* mbi =
                            new MiniBatchInfo<AlignmentGroup<FragT*>>(batchNum, alignments, logForgettingMass);
                        workQueue.push(mbi);
                        {
                            std::unique_lock<std::mutex> l(cvmutex);
                            workAvailable.notify_one();
                        }
                        alignments = new std::vector<AlignmentGroup<FragT*>*>;
                        alignments->reserve(miniBatchSize);
                    }
                    if (numProc % 1000000 == 0) {
                        const char RESET_COLOR[] = "\x1b[0m";
                        char green[] = "\x1b[30m";
                        green[3] = '0' + static_cast<char>(fmt::GREEN);
                        char red[] = "\x1b[30m";
                        red[3] = '0' + static_cast<char>(fmt::RED);
                        fmt::print(stderr, "\r\r{}processed{} {} {}reads{}", green, red, numProc, green, RESET_COLOR);
                    }
                    ++numProc;
                    alignmentGroupsRemain = bq.getAlignmentGroup(ag);

                }
                std::cerr << "\n";

                // Free the alignments and the vector holding them
                for (auto& aln : *alignments) {
                    aln->alignments().clear();
                    delete aln; aln = nullptr;
                }
                delete alignments;

                doneParsing = true;

                /**
                 * This could be a problem for small sets of alignments --- make sure the
                 * work queue is empty!!
                 * --- Thanks for finding a dataset that exposes this bug, Richard (Smith-Unna)!
                 */
                size_t t = 0;
                while (!workQueue.empty()) {
                    std::unique_lock<std::mutex> l(cvmutex);
                    workAvailable.notify_one();
                }

                size_t tnum{0};
                for (auto& t : workers) {
                    fmt::print(stderr, "killing thread {} . . . ", tnum++);
                    {
                        std::unique_lock<std::mutex> l(cvmutex);
                        workAvailable.notify_all();
                    }
                    t.join();
                    fmt::print(stderr, "done\r\r");
                }
                fmt::print(stderr, "\n");
                consumedAllInput = true;

                numObservedFragments += alnLib.numMappedReads();
                fmt::print(stderr, "# observed = {} mapped fragments.\033[F\033[F\033[F\033[F",
                        numObservedFragments);

                fmt::print(stderr, "Waiting on output thread\n");
                outputThread.join();
                fmt::print(stderr, "done\n");

                fmt::print(stderr, "\n\n\n\n");
                return true;
            }



    } // namespace sampler
} // namespace salmon

#endif // __SAMPLER__HPP__
