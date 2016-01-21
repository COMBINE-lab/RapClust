/**
>HEADER
    Copyright (c) 2013 Rob Patro robp@cs.cmu.edu

    This file is part of Sailfish.

    Sailfish is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Sailfish is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Sailfish.  If not, see <http://www.gnu.org/licenses/>.
<HEADER
**/


#ifndef COLLAPSED_ITERATIVE_OPTIMIZER_HPP
#define COLLAPSED_ITERATIVE_OPTIMIZER_HPP

#include <algorithm>
#include <cassert>
#include <cmath>
#include <unordered_map>
#include <map>
#include <vector>
#include <unordered_set>
#include <mutex>
#include <thread>
#include <sstream>
#include <exception>
#include <random>
#include <queue>
#include "btree_map.h"

/** Boost Includes */
#include <boost/container/flat_map.hpp>
#include <boost/dynamic_bitset/dynamic_bitset.hpp>
#include <boost/range/irange.hpp>
#include <boost/program_options.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/heap/fibonacci_heap.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/framework/accumulator_set.hpp>
#include <boost/accumulators/statistics/p_square_quantile.hpp>
#include <boost/accumulators/statistics/count.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/weighted_mean.hpp>
#include <boost/lockfree/queue.hpp>
#include <boost/thread/thread.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/distributions/beta.hpp>


//#include <Eigen/Core>

#include "tbb/concurrent_unordered_set.h"
#include "tbb/concurrent_vector.h"
#include "tbb/concurrent_unordered_map.h"
#include "tbb/concurrent_queue.h"
#include "tbb/parallel_for.h"
#include "tbb/parallel_for_each.h"
#include "tbb/parallel_reduce.h"
#include "tbb/blocked_range.h"
#include "tbb/task_scheduler_init.h"
#include "tbb/partitioner.h"

#include "BiasIndex.hpp"
#include "ezETAProgressBar.hpp"
#include "LookUpTableUtils.hpp"
#include "SailfishMath.hpp"

template <typename ReadHash>
class CollapsedIterativeOptimizer {

private:
    /**
    * Type aliases
    */
    using TranscriptID = uint32_t;
    using KmerID =  uint64_t;
    using Count = uint32_t;
    using KmerQuantity = double;
    using Promiscutity = double;
    using KmerMap = tbb::concurrent_unordered_map< uint64_t, tbb::concurrent_vector<uint32_t> >;

    struct TranscriptGeneVectors;
    using TranscriptIDVector = std::vector<TranscriptID>;
    using KmerIDMap = std::vector<TranscriptIDVector>;


    using TranscriptKmerSet = std::tuple<TranscriptID, std::vector<KmerID>>;
    using StringPtr = std::string*;
    using TranscriptScore = uint64_t;
   // using HashArray = jellyfish::invertible_hash::array<uint64_t, atomic::gcc, allocators::mmap>;
    using ReadLength = size_t;

    // Necessary forward declaration
    struct TranscriptData;
    using HeapPair = std::tuple<TranscriptScore, TranscriptID>;
    using Handle = typename boost::heap::fibonacci_heap<HeapPair>::handle_type;
    using BlockedIndexRange =  tbb::blocked_range<size_t>;

    struct TranscriptGeneVectors {
        tbb::concurrent_vector<uint32_t> transcripts;
        tbb::concurrent_vector<uint32_t> genes;
    };

    struct TranscriptData {
        TranscriptID id;
        StringPtr header;
        std::map<KmerID, KmerQuantity> binMers;
        KmerQuantity mean;
        size_t length;
        bool logSpace;
    };

    enum class CountSpace { LogSpace, NonLogSpace };

    class TranscriptInfo {
      public:

        TranscriptInfo() : binMers(boost::container::flat_map<KmerID, KmerQuantity>()),
                         //binMers(std::unordered_map<KmerID, KmerQuantity>()),
                         //logLikes(std::vector<KmerQuantity>()),
                         //weights(std::vector<KmerQuantity>()),
                           mean(0.0), fracLow(0.0), fracHigh(0.0), length(0), effectiveLength(0),
                           logInvEffectiveLength(sailfish::math::LOG_0){ updated.store(0); /* weightNum.store(0); totalWeight.store(0.0);*/ }

      TranscriptInfo(TranscriptInfo&& other) {
        std::swap(binMers, other.binMers);
        //std::swap(weights, other.weights);
        //std::swap(logLikes, other.logLikes);
        //totalWeight.store(other.totalWeight.load());
        //weightNum.store(other.weightNum.load());
        updated.store(other.updated.load());
        mean = other.mean;
        length = other.length;
        effectiveLength = other.effectiveLength;
      }
        void resetMass() {
            double retMass;
            switch (countSpace) {
            case CountSpace::LogSpace:
                retMass = sailfish::math::LOG_0;
                do {
                    retMass = totalMass.compare_and_swap(sailfish::math::LOG_0, totalMass);
                } while (retMass != sailfish::math::LOG_0);
                break;
            case CountSpace::NonLogSpace:
                retMass = 0.0;
                do {
                    retMass = totalMass.compare_and_swap(0.0, totalMass);
                } while (retMass != 0.0);
            }
        }

        void addMass(KmerQuantity mass) {
            double origMass, returnedMass, newMass;
            size_t tries{0};
            switch (countSpace) {
            case CountSpace::LogSpace:
                origMass = totalMass;
                returnedMass = totalMass;
                newMass = sailfish::math::LOG_0;
                do {
                    origMass = returnedMass;
                    newMass = sailfish::math::logAdd(origMass, mass);
                    returnedMass = totalMass.compare_and_swap(newMass, origMass);
                } while (returnedMass != origMass);
                break;
            case CountSpace::NonLogSpace:
                origMass = totalMass;
                returnedMass = totalMass;
                newMass = sailfish::math::LOG_0;
                do {
                    if (tries > 10) {
                      boost::this_thread::sleep_for(boost::chrono::microseconds(10));
                      std::cerr << "trying; expected: " << origMass << ", returned: " << returnedMass << "\n";
                      tries = 0;
                      std::cerr << "trying\n";
                    }
                    origMass = returnedMass;
                    newMass = origMass + mass;
                    returnedMass = totalMass.compare_and_swap(newMass, origMass);
                    ++tries;
                } while (returnedMass != origMass);
            }
        }

        //std::atomic<double> totalWeight;
        //btree::btree_map<KmerID, KmerQuantity> binMers;
        //std::vector<KmerQuantity> weights;
        //std::vector<KmerQuantity> logLikes;
        //std::atomic<uint32_t> weightNum;
        std::atomic<uint32_t> updated;
        //std::unordered_map<KmerID, KmerQuantity> binMers;
        boost::container::flat_map<KmerID, KmerQuantity> binMers;
        KmerQuantity mean;
        KmerQuantity fracLow, fracHigh;
        tbb::atomic<KmerQuantity> totalMass;
        ReadLength length;
        ReadLength effectiveLength;
        double logInvEffectiveLength;
        CountSpace countSpace;
        bool isAnchored;
    };

    // This struct represents a "job" (transcript) that needs to be processed
    struct TranscriptJob {
        StringPtr header;
        StringPtr seq;
        TranscriptID id;
    };

    struct TranscriptResult {
        TranscriptData *data;
        TranscriptKmerSet *ks;
    };

    struct BinmerUpdates {
        std::vector<KmerID> zeroedBinmers;
        std::vector<KmerID> updatedBinmers;
    };

    bool loggedCounts_;
    uint32_t numThreads_;
    size_t merLen_;
    ReadHash & readHash_;
    BiasIndex& biasIndex_;
    std::string kmerEquivClassFname_;

    // The number of occurences above whcih a kmer is considered promiscuous
    size_t promiscuousKmerCutoff_ {std::numeric_limits<size_t>::max()};

    // Map each kmer to the set of transcripts it occurs in
    KmerIDMap transcriptsForKmer_;

    // The actual data for each transcript
    std::vector<TranscriptInfo> transcripts_;

    TranscriptGeneMap& transcriptGeneMap_;

    tbb::concurrent_unordered_set<uint64_t> genePromiscuousKmers_;

    std::vector<Promiscutity> kmerGroupPromiscuities_;
    std::vector<Promiscutity> kmerGroupBiases_;
    std::vector<KmerQuantity> kmerGroupCounts_;
    std::vector<KmerQuantity> logKmerGroupCounts_;
    std::vector<Count> kmerGroupSizes_;


    /**
     * logs all counts in the kmerGroupCounts_ variable.  Groups with a
     * count of zero are assigned the special value LOG_0.
     */
    /*
    inline void logKmerGroupCounts_() {
        std::for_each(kmerGroupCounts_.begin(), kmerGroupCounts_.end(),
                      [](KmerQuantity& count) {
                          if (count > 0) {
                              count = std::log(count);
                          } else {
                              count = sailfish::math::LOG_0;
                          }
                      });
    }
    */

    /**
     * Compute the "Inverse Document Frequency" (IDF) of a kmer within a set of transcripts.
     * The inverse document frequency is the log of the number of documents (i.e. transcripts)
     * divded by the number of documents containing this term (i.e. kmer).
     * @param  k [kmer id for which the IDF should be computed]
     * @return   [IDF(k)]
     */
    inline double _idf( uint64_t k ) {
        double df = transcriptsForKmer_[k].size();
        return (df > 0.0) ? std::log(transcripts_.size() / df) : 0.0;
    }

    /**
     * Returns true if this kmer should be considered in our estimation, false
     * otherwise
     * @param  mer [kmer id to test for consideration]
     * @return     [true if we consider the kmer with this id, false otherwise]
     */
    inline bool _considered( uint64_t mer ) {
        // The kmer is only considered if it exists in the transcript set
        // (i.e. it's possible to cover) and it's less prmiscuous than the
        // cutoff.
        return true;
    }

    /**
     * The weight attributed to each appearence of the kmer with the given ID.
     * If the kmer with ID k occurs in m different transcripts, then
     * weight_(k) = 1 / m.
     * @param  k [The ID of the kmer whose weight is to be computed]
     * @return   [The weight of each appearance of the kmer with the given ID]
     */
    KmerQuantity weight_( KmerID k ) {
        return 1.0 / (kmerGroupPromiscuities_[k] );
    }

    KmerQuantity _computeMedian( const TranscriptInfo& ti ) {

      using namespace boost::accumulators;
      using Accumulator = accumulator_set<double, stats<tag::median(with_p_square_quantile)>>;

      Accumulator acc;
      for (auto binmer : ti.binMers) {
        acc(binmer.second);
      }

      return median(acc);
    }

    /**
     * Computes the sum of kmer counts within the transcript given by ti, but clamping
     * all non-zero counts to the given quantile.  For example, if quantile was 0.25, and
     * x and y represented the 1st and 3rd quantile of kmer counts, then every nonzero count c
     * would be transformed as c = max(x, min(y,c));
     *
     * @param  ti       [description]
     * @param  quantile [description]
     * @return          [description]
     */
    KmerQuantity _computeSumQuantile( const TranscriptInfo& ti, double quantile ) {
        using namespace boost::accumulators;
        using accumulator_t = accumulator_set<double, stats<tag::p_square_quantile> >;
        //using tail_t = accumulator_set<double, stats<tag::tail
        KmerQuantity sum = 0.0;


        accumulator_t accLow(quantile_probability = quantile);
        accumulator_t accHigh(quantile_probability = 1.0-quantile);
        for ( auto binmer : ti.binMers ) {
            accLow(binmer.second);
            accHigh(binmer.second);
        }

        auto cutLow = p_square_quantile(accLow);
        auto cutHigh = p_square_quantile(accHigh);

        for ( auto binmer : ti.binMers ) {
            KmerQuantity res = std::min( cutHigh, std::max(cutLow, binmer.second ) );
            if (res != binmer.second) {
              std::cerr << "cutLow = " << cutLow << ", cutHigh = " << cutHigh << " :: ";
              std::cerr << "res = " << res << ", binmer.second = " << binmer.second << "\n";
            }
            sum += std::min( cutHigh, std::max(cutLow, binmer.second ) );
        }
        return sum;
    }

    KmerQuantity _computeSum( const TranscriptInfo& ti ) {
        KmerQuantity sum = 0.0;
        for ( auto binmer : ti.binMers ) {
          sum += kmerGroupBiases_[binmer.first] * binmer.second;
        }
        return sum;
    }

    KmerQuantity _computeSumClamped( const TranscriptInfo& ti ) {
        if (ti.binMers.size() < 5) {return _computeSum(ti);}
        KmerQuantity sum = 0.0;

        auto maxQuant = std::numeric_limits<KmerQuantity>::max();
        auto minQuant = 0.0;

        auto startValue = ti.binMers.begin()->second;
        KmerQuantity lowestCount = startValue;
        KmerQuantity secondLowestCount = startValue;

        KmerQuantity highestCount = startValue;
        KmerQuantity secondHighestCount = startValue;

        for ( auto binmer : ti.binMers ) {
          auto prevLowest = lowestCount;
          lowestCount = std::min(binmer.second, lowestCount);
          if (lowestCount < prevLowest) { secondLowestCount = prevLowest; }
          auto prevHighest = highestCount;
          highestCount = std::max(binmer.second, highestCount);
          if (highestCount > prevHighest) { secondHighestCount = prevHighest; }
        }
        //std::cerr << lowestCount << ", " << secondLowestCount << ", " << highestCount << ", " << secondHighestCount << "\n";

        for ( auto binmer : ti.binMers ) {
          if (binmer.second > lowestCount and binmer.second < highestCount) {
            sum += kmerGroupBiases_[binmer.first] * binmer.second;
          } else if (binmer.second >= highestCount ) {
            sum += secondHighestCount;
          } else if (binmer.second <= lowestCount ) {
            sum += secondLowestCount;
          }
        }

        return sum;
    }


    bool _discard( const TranscriptInfo& ti) {
        if ( ti.mean == 0.0 ) {
            return false;
        } else {
            ti.mean = 0.0;
            ti.binMers.clear();
            return true;
        }
    }

    KmerQuantity _computeSumVec( const TranscriptInfo& ti ) {
        KmerQuantity sum = 0.0;
        for ( auto w : ti.weights ) {
          sum += w;
        }
        return sum;
    }

    KmerQuantity _computeClampedMean( const TranscriptInfo& ti ) {
        return (ti.effectiveLength > 0.0) ? (_computeSumClamped(ti) / ti.effectiveLength) : 0.0;
    }

    KmerQuantity _computeMean( const TranscriptInfo& ti ) {
        return (ti.effectiveLength > 0.0) ? (_computeSum(ti) / ti.effectiveLength) : 0.0;
        //return (ti.effectiveLength > 0.0) ? (ti.totalWeight.load() / ti.effectiveLength) : 0.0;
        //return (ti.effectiveLength > 0.0) ? (_computeSumVec(ti) / ti.effectiveLength) : 0.0;
    }

    KmerQuantity _computeWeightedMean( const TranscriptInfo& ti ) {
        using namespace boost::accumulators;
        accumulator_set<double, stats<tag::count, tag::weighted_mean>, double> acc;

        for ( auto binmer : ti.binMers ) {
          if ( this->genePromiscuousKmers_.find(binmer.first) == this->genePromiscuousKmers_.end() ){
            acc(binmer.second, weight=kmerGroupBiases_[binmer.first] * weight_(binmer.first));
          }
        }

        auto nnz = count(acc);

        if ( nnz < ti.effectiveLength ) {
            acc(0.0, weight=ti.effectiveLength-nnz);
        }

        auto sum = sum_of_weights(acc);
        return sum > 0.0 ? weighted_mean(acc) : 0.0;
    }

    double _effectiveLength( const TranscriptInfo& ts ) {
        double length = 0.0;
        for ( auto binmer : ts.binMers ) {
            length += weight_(binmer.first);
        }
        return length;
    }

    template <typename T>
    T dotProd_(std::vector<T>& u, std::vector<T>& v) {

      auto dot = tbb::parallel_reduce(
        BlockedIndexRange(size_t(0), v.size()),
            T(0.0),  // identity element for summation
            [&]( const BlockedIndexRange& r, T current_sum ) -> T {
             for (size_t i=r.begin(); i!=r.end(); ++i) {
               current_sum += (u[i]*v[i]);
             }
             return current_sum; // body returns updated value of the accumulator
             },
             []( double s1, double s2 ) -> double {
                return s1+s2;       // "joins" two accumulated values
      });

      return dot;
    }

    void normalizeTranscriptMeans_(){
        //auto sumMean = 0.0;
        //for ( auto ti : transcripts_ ) { sumMean += ti.mean; }

        auto sumMean = tbb::parallel_reduce(
            BlockedIndexRange(size_t(0), transcripts_.size()),
            double(0.0),  // identity element for summation
            [&, this]( const BlockedIndexRange& r, double current_sum ) -> double {
                 for (size_t i=r.begin(); i!=r.end(); ++i) {
                     double x = this->transcripts_[i].mean;
                     current_sum += x;
                 }
                 return current_sum; // body returns updated value of the accumulator
             },
             []( double s1, double s2 ) -> double {
                 return s1+s2;       // "joins" two accumulated values
             });

        // compute the new mean for each transcript
        tbb::parallel_for(BlockedIndexRange(size_t(0), size_t(transcripts_.size())),
            [this, sumMean](const BlockedIndexRange& range) -> void {
              for(size_t tid = range.begin(); tid != range.end(); ++tid) {
                this->transcripts_[tid].mean /= sumMean;
              }
            });

    }

    template <typename T>
    bool hasNan_(std::vector<T>& v) {
        std::atomic<bool> found{false};
        tbb::parallel_for(BlockedIndexRange(size_t(0), size_t(v.size())),
                          [this, &v, &found](const BlockedIndexRange& range) -> void {
                              for(size_t i = range.begin(); i != range.end(); ++i) {
                                  if (found) { return; }
                                  if (!std::isfinite(v[i])) { found = true; return; }
                              }
        });
        return found.load();
    }

    template <typename T>
    T psumAccumulate(const tbb::blocked_range<T*>& r, T value) {
            return std::accumulate(r.begin(),r.end(),value);
    }

    template <typename T>
    T psum_(std::vector<T>& v) {
      //auto func = std::bind( std::mem_fn(&CollapsedIterativeOptimizer<ReadHash>::psumAccumulate<T>),
      //                       this, std::placeholders::_1, std::placeholders::_2 );
      auto func = [this](const tbb::blocked_range<T*>& r, T value) -> T {
          return this->psumAccumulate(r, value);
      };
      auto sum = tbb::parallel_reduce(
        tbb::blocked_range<T*>(&v[0], &v[v.size()]),
          T{0},  // identity element for summation
          func,
          std::plus<T>()
        );
      return sum;
    }

    template <typename T>
    T pdiff_(std::vector<T>& v0, std::vector<T>& v1) {
        auto diff = tbb::parallel_reduce(
        BlockedIndexRange(size_t(0), v0.size()),
          double(0.0),  // identity element for difference
          [&](const BlockedIndexRange& r, double currentDiff ) -> double {
            for (size_t i=r.begin(); i!=r.end(); ++i) {
              currentDiff += v0[i] - v1[i];
            }
            return currentDiff; // body returns updated value of the accumulator
          },
          []( double s1, double s2 ) -> double {
               return s1+s2;       // "joins" two accumulated values
          }
        );

        return diff;
    }

    template <typename T>
    T pabsdiff_(std::vector<T>& v0, std::vector<T>& v1) {
        auto diff = tbb::parallel_reduce(
        BlockedIndexRange(size_t(0), v0.size()),
          double(0.0),  // identity element for difference
          [&]( const BlockedIndexRange& r, double currentDiff ) -> double {
            for (size_t i=r.begin(); i!=r.end(); ++i) {
              currentDiff = std::abs(v0[i] - v1[i]);
            }
            return currentDiff; // body returns updated value of the accumulator
          },
          []( double s1, double s2 ) -> double {
               return s1+s2;       // "joins" two accumulated values
          }
        );

        return diff;
    }

    template <typename T>
    std::vector<T> relAbsDiff_(std::vector<T>& v0, std::vector<T>& v1, T minVal) {
        std::vector<T> relDiff(v0.size(), T());
        // compute the new mean for each transcript
        tbb::parallel_for(BlockedIndexRange(size_t(0), size_t(v0.size())),
                          [&v0, &v1, &relDiff, minVal](const BlockedIndexRange& range ) -> void {
                              for (auto tid : boost::irange(range.begin(), range.end())) {
                                  T oldVal = v0[tid];
                                  T newVal = v1[tid];
                                  if (oldVal > minVal or newVal > minVal) {
                                      relDiff[tid] = std::abs(newVal - oldVal) / oldVal;
                                  }
                              }
                          });
        return relDiff;
    }

    void normalize_(std::vector<double>& means) {
        auto sumMean = psum_(means);
        auto invSumMean = (sumMean > 0.0) ? 1.0 / sumMean : 0.0;

        // compute the new mean for each transcript
        tbb::parallel_for(BlockedIndexRange(size_t(0), size_t(transcripts_.size())),
            [&means, invSumMean](const BlockedIndexRange& range ) -> void {
              for (auto tid : boost::irange(range.begin(), range.end())) {
                means[tid] *= invSumMean;
              }
            });
    }

    double averageCount(const TranscriptInfo& ts){
        if ( ts.binMers.size() == 0 ) { return 0.0; }
        double sum = 0.0;
        for ( auto binmer : ts.binMers ) {
            sum += kmerGroupBiases_[binmer.first] * binmer.second;
        }
        return sum / ts.binMers.size();

    }

    /**
     * This function should be run only after <b>after</b> an EM loop.
     * It estimates kmer specific biases based on how much a kmer's count
     * deviates from it's corresponding transcript's mean
     */
    void computeKmerFidelities_() {

      // For each transcript, compute it's overall fidelity.  This is related
      // to the variance of coverage across the transcript.  The more uniform
      // the coverage, the more we believe the transcript.
        std::vector<double> transcriptFidelities(transcripts_.size(), 0.0);
        tbb::parallel_for( BlockedIndexRange(size_t(0), transcripts_.size()),
            [this, &transcriptFidelities]( const BlockedIndexRange& range ) -> void {
               for (auto tid = range.begin(); tid != range.end(); ++tid) {
                double sumDiff = 0.0;
                //if (tid >= this->transcripts_.size()) { std::cerr << "attempting to access transcripts_ out of range\n";}
                auto& ts = this->transcripts_[tid];

                //std::cerr << "transcript " << tid << "\n";
                for ( auto& b : ts.binMers ) {
                    //if (b.first >= this->kmerGroupBiases_.size()) { std::cerr << "attempting to access kmerGroupBiases_ out of range\n";}
                    auto scaledMean = this->kmerGroupSizes_[b.first] * ts.mean;
                    auto diff = std::abs(b.second - scaledMean);
                    sumDiff += diff;//*diff;
                }
                // The rest of the positions have 0 coverage have an error
                // of |0 - \mu_t| = \mu_t.  There are l(t) - ts.binMers.size() of these.
                sumDiff += ts.mean * (ts.length - ts.binMers.size());
                auto fidelity = (ts.length > 0.0) ? sumDiff / ts.length : 0.0;
                fidelity = 1.0 / (fidelity + 1.0);
                //if (tid >= transcriptFidelities.size()) { std::cerr << "attempting to access transcriptFidelities out of range\n";}
                transcriptFidelities[tid] = fidelity;
                //std::cerr << "fidelity (" << tid << ") = " << fidelity << "\n";
               }
            });

        tbb::parallel_for(BlockedIndexRange(size_t(0), transcriptsForKmer_.size()),
            [this, &transcriptFidelities](const BlockedIndexRange& range) -> void {
              // Each transcript this kmer group appears in votes on the bias of this kmer.
              // Underrepresented kmers get bias values > 1.0 while overrepresented kmers get
              // bias values < 1.0.  The vote of each transcript is weigted by it's fidelity
                for (auto kid = range.begin(); kid != range.end(); ++kid) {
                  double totalBias = 0.0;
                  double totalFidelity = 0.0;
                  for( auto tid : this->transcriptsForKmer_[kid] ) {
                    auto& transcript = this->transcripts_[tid];
                    auto fidelity = transcriptFidelities[tid];
                    auto totalMean = transcript.mean * this->kmerGroupSizes_[kid];
                    auto curAlloc = transcript.binMers[kid];
                    totalBias += (curAlloc > 0.0) ? fidelity * (totalMean / curAlloc) : 0.0;
                    totalFidelity += fidelity;
                  }
                  double bias = totalBias / totalFidelity;
                  double alpha = 0.25; //std::min( 1.0, 10.0*averageFidelity / confidence);
                  double prevBias = this->kmerGroupBiases_[kid];
                  this->kmerGroupBiases_[kid] = alpha * bias + (1.0 - alpha) * prevBias;
              }
            }
        );


    }

    double logLikelihood_(std::vector<double>& means) {

      const auto numTranscripts = transcripts_.size();
      std::vector<double> likelihoods(numTranscripts, 0.0);

        // Compute the log-likelihood
        tbb::parallel_for(BlockedIndexRange(size_t(0), numTranscripts),
          // for each transcript
          [&likelihoods, &means, this](const BlockedIndexRange& range) ->void {
            auto epsilon = 1e-40;
            for (auto tid = range.begin(); tid != range.end(); ++tid) {
              auto& ti = transcripts_[tid];
              double relativeAbundance = means[tid];

              if (ti.binMers.size() > 0 ) { likelihoods[tid] = 1.0; }
              // For each kmer in this transcript
              for ( auto& binmer : ti.binMers ) {
                likelihoods[tid] *=  binmer.second /
                           (this->kmerGroupBiases_[binmer.first] * this->kmerGroupCounts_[binmer.first]);
              }
              likelihoods[tid] = (relativeAbundance > epsilon and likelihoods[tid] > epsilon) ?
                                 std::log(relativeAbundance * likelihoods[tid]) : 0.0;
            }
          });

      return psum_(likelihoods);

    }

    double expectedLogLikelihood_(std::vector<double>& means) {
      // alpha in arXiv:1104.3889v2
      std::vector<double> sampProbs(means.size(), 0.0);

      tbb::parallel_for(BlockedIndexRange(size_t(0), means.size()),
          // for each transcript
          [&sampProbs, &means, this](const BlockedIndexRange& range) ->void {
            for (auto tid = range.begin(); tid != range.end(); ++tid) {
              auto& ti = transcripts_[tid];
              double relativeAbundance = (ti.effectiveLength > 0.0) ? means[tid] / ti.effectiveLength : 0.0;
              //sampProbs[tid] = ti.length * relativeAbundance;
              sampProbs[tid] = relativeAbundance;
            }});

      normalize_(sampProbs);

      //return (means.size() > 0) ? logLikelihood2_(sampProbs) / means.size() : 0.0;
      return (means.size() > 0) ? logLikelihood3_(sampProbs) / means.size() : 0.0;
    }

    double logLikelihood3_(std::vector<double>& sampProbs) {

      std::vector<double> likelihoods(transcriptsForKmer_.size(), 0.0);

        // Compute the log-likelihood
        tbb::parallel_for(BlockedIndexRange(size_t(0), size_t(transcriptsForKmer_.size())),
          // for each transcript
          [&likelihoods, &sampProbs, this](const BlockedIndexRange& range) ->void {
            for (auto kid = range.begin(); kid != range.end(); ++kid) {
              double kmerLikelihood = 0.0;
              KmerQuantity totalKmerMass = kmerGroupCounts_[kid];
              for (auto& tid : this->transcriptsForKmer_[kid]) {
                  // double logProbSampleTID = (sampProbs[tid] > sailfish::math::EPSILON) ?
                  //     std::log(sampProbs[tid]) : sailfish::math::LOG_0;

                  double logProbSampleTID = (sampProbs[tid] > sailfish::math::EPSILON) ?
                      //std::log(sampProbs[tid]) : sailfish::math::LOG_0;
                  std::log(sampProbs[tid] / this->transcripts_[tid].effectiveLength) : sailfish::math::LOG_0;

                  if (logProbSampleTID != sailfish::math::LOG_0) {
                      kmerLikelihood += kmerGroupCounts_[kid] * logProbSampleTID;
                  }
              }
              // If this k-mer is present
              if (totalKmerMass > 0) {

                  if (kmerLikelihood == sailfish::math::LOG_0) {
                      std::cerr << "kmer group: " << kid << " has probability (" << kmerLikelihood << "); too low\n";
                  } else {
                      likelihoods[kid] = kmerLikelihood;
                  }

              }
            }
          });

      return psum_(likelihoods);

    }

    double logLikelihood2_(std::vector<double>& sampProbs) {

      std::vector<double> likelihoods(transcriptsForKmer_.size(), 0.0);

        // Compute the log-likelihood
        tbb::parallel_for(BlockedIndexRange(size_t(0), size_t(transcriptsForKmer_.size())),
          // for each transcript
          [&likelihoods, &sampProbs, this](const BlockedIndexRange& range) ->void {
            for (auto kid = range.begin(); kid != range.end(); ++kid) {
              double kmerLikelihood = 0.0;
              KmerQuantity totalKmerMass = 0.0;
              for (auto& tid : this->transcriptsForKmer_[kid]) {
                double kmerMass{this->transcripts_[tid].binMers[kid]};
                kmerLikelihood += kmerMass * (sampProbs[tid] / this->transcripts_[tid].length);
                totalKmerMass += kmerMass;
              }
              // If
              if (totalKmerMass > 0) {
                  if (kmerLikelihood < 1e-20) {
                      std::cerr << "kmer group: " << kid << " has probability (" << kmerLikelihood << "); too low\n";
                  } else {
                      likelihoods[kid] = std::log(kmerLikelihood);
                  }
              }
            }
          });

      return psum_(likelihoods);

    }


    // Since there's no built in hash code for vectors
  template <typename T>
  class my_hasher{
  public:
    size_t operator()(const T& x) const {
        if (x.size() == 0 ) { return 0; }
        size_t seed = x[0];
        for (auto i : boost::irange(size_t(1), x.size())) {
            boost::hash_combine(seed, static_cast<size_t>(x[i]));
        }
        return seed;
    }
  };

/**
 * Collapses all kmer which share the same transcript multiset.  Such kmers can be
 * treated as a "batch" with a count whose value is the sum of individual "batch"
 * members.
 * @param  isActiveKmer       [A bitvector which designates, for each kmer,
 *                             whether or not that kmer is active in the current
 *                             read set.]
 */
 void collapseKmers_( boost::dynamic_bitset<>& isActiveKmer ) {

    auto numTranscripts = transcriptGeneMap_.numTranscripts();

    /**
     * Map from a vector of transcript IDs to the list of kmers that have this
     * transcript list.  This allows us to collapse all kmers that exist in the
     * exact same set of transcripts into a single kmer group.
     */
    tbb::concurrent_unordered_map< TranscriptIDVector,
                                   tbb::concurrent_vector<KmerID>,
                                   my_hasher<std::vector<TranscriptID>> > m;

     // Asynchronously print out the progress of our hashing procedure
     std::atomic<size_t> prog{0};
     std::thread t([this, &prog]() -> void {
        ez::ezETAProgressBar pb(this->transcriptsForKmer_.size());
        pb.start();
        size_t prevProg{0};
        while ( prevProg < this->transcriptsForKmer_.size() ) {
            if (prog > prevProg) {
                auto diff = prog - prevProg;
                pb += diff;
                prevProg += diff;
            }
            boost::this_thread::sleep_for(boost::chrono::seconds(1));
        }
        if (!pb.isDone()) { pb.done(); }
     });

     //For every kmer, compute it's kmer group.
     tbb::parallel_for(BlockedIndexRange(size_t(0), transcriptsForKmer_.size()),
        [&](const BlockedIndexRange& range ) -> void {
          for (auto j = range.begin(); j != range.end(); ++j) {
            if (isActiveKmer[j]) {
              m[ transcriptsForKmer_[j]  ].push_back(j);
            }
            ++prog;
          }
     });

     // wait for the parallel hashing to finish
     t.join();

     // TESTING
     std::ofstream eqFile("KMER_EQUIV_CLASSES.txt");
     for (auto& kv : m ) {
       for (auto e : kv.second) { eqFile << e << "\t"; }
       eqFile << "\n";
     }
     eqFile.close();
     // END TESTING

     std::cerr << "Out of " << transcriptsForKmer_.size() << " potential kmers, "
               << "there were " << m.size() << " distinct groups\n";

     size_t totalKmers = 0;
     size_t index = 0;
     std::vector<KmerQuantity> kmerGroupCounts(m.size());
     std::vector<Promiscutity> kmerGroupPromiscuities(m.size());
     std::vector<TranscriptIDVector> transcriptsForKmer(m.size());
     kmerGroupSizes_.resize(m.size(), 0);

     using namespace boost::accumulators;
     std::cerr << "building collapsed transcript map\n";
     for ( auto& kv : m ) {

        // For each transcript covered by this kmer group, add this group to the set of kmer groups contained in
        // the transcript.  For efficiency, we also compute the kmer promiscuity values for each kmer
        // group here --- the promiscuity of a kmer group is simply the number of distinct transcripts in
        // which this group of kmers appears.
        auto prevTID = std::numeric_limits<TranscriptID>::max();
        KmerQuantity numDistinctTranscripts = 0.0;
        for ( auto& tid : kv.first ) {
          transcripts_[tid].binMers[index] += 1;
          // Since the transcript IDs are sorted we just have to check
          // if this id is different from the previous one
          if (tid != prevTID) { numDistinctTranscripts += 1.0; }
          prevTID = tid;
        }
        // Set the promiscuity and the set of transcripts for this kmer group
        kmerGroupPromiscuities[index] = numDistinctTranscripts;
        transcriptsForKmer[index] = kv.first;

        // Aggregate the counts attributable to each kmer into its repective
        // group's counts.
        for (auto kid : kv.second) {
            kmerGroupCounts[index] += readHash_.atIndex(kid);
        }
        kmerGroupSizes_[index] = kv.second.size();

        // Update the total number of kmers we're accounting for
        // and the index of the current kmer group.
        totalKmers += kv.second.size();
        ++index;
      }

      std::cerr << "Verifying that the unique set encodes " << totalKmers << " kmers\n";
      std::cerr << "collapsedCounts.size() = " << transcriptsForKmer.size() << "\n";

      // update the relevant structures holding info for the full kmer
      // set with those holding the info for our collapsed kmer sets
      std::swap(kmerGroupPromiscuities, kmerGroupPromiscuities_);
      std::swap(kmerGroupCounts, kmerGroupCounts_);
      std::swap(transcriptsForKmer, transcriptsForKmer_);

      /*
      uint64_t groupCounts = 0;
      for(auto c : kmerGroupCounts_) { groupCounts += c; }
      auto tmp = psum_(kmerGroupCounts_);
      std::cerr << "groupCount(" << groupCounts << ") - parallelCount(" << tmp << ") = " << groupCounts - tmp << "\n";
      std::atomic<uint64_t> individualCounts{0};
      tbb::parallel_for(size_t{0}, readHash_.size(),
        [&](size_t i) {
          if (isActiveKmer[i]) {
            individualCounts += readHash_.atIndex(i);
          } });
      auto diff = groupCounts - individualCounts;
      std::cerr << "groupTotal(" << groupCounts << ") - totalNumKmers(" << individualCounts << ") = " << diff << "\n";
      */
  }


  /**
   * NEW! Assuming that equivalence classes were computed in the index
   **/
  void prepareCollapsedMaps_(
                            const std::string& kmerEquivClassFname,
                            bool discardZeroCountKmers) {

    using std::cerr;

    cerr << "reading Kmer equivalence classes \n";

    auto memberships = LUTTools::readKmerEquivClasses(kmerEquivClassFname);
    size_t numKmers{readHash_.size()};
    size_t numKmerClasses{(*std::max_element(memberships.begin(), memberships.end())) + 1};
    boost::dynamic_bitset<> isActiveKmer(numKmers);

    kmerGroupCounts_.resize(numKmerClasses, 0.0);
    logKmerGroupCounts_.resize(numKmerClasses, sailfish::math::LOG_0);
    kmerGroupSizes_.resize(numKmerClasses, 0.0);
    kmerGroupPromiscuities_.resize(numKmerClasses, 0.0);

    cerr << "updating Kmer group counts\n";
    // Update the kmer group counts using the information from the read hash
    for (auto kid : boost::irange(size_t{0}, numKmers)) {

      size_t count = readHash_.atIndex(kid);
      if (!discardZeroCountKmers or count != 0) {
        auto kmerClassID = memberships[kid];
        isActiveKmer[kmerClassID] = true;
        kmerGroupCounts_[kmerClassID] += count;
        kmerGroupSizes_[kmerClassID]++;
      }

    }

    cerr << "updating transcript map\n";
    for (auto kmerClassID : boost::irange(size_t{0}, numKmerClasses)) {

      // For each transcript covered by this kmer group, add this group to the set of kmer groups contained in
      // the transcript.  For efficiency, we also compute the kmer promiscuity values for each kmer
      // group here --- the promiscuity of a kmer group is simply the number of distinct transcripts in
      // which this group of kmers appears.
      auto prevTID = std::numeric_limits<TranscriptID>::max();
      KmerQuantity numDistinctTranscripts = 0.0;

      //cerr << "numKmerClasses = " << numKmerClasses << ", kmerClassID = " << kmerClassID << ", transcriptsForKmer_.size() = " << transcriptsForKmer_.size() << "\n";
      for (auto& tid : transcriptsForKmer_[kmerClassID]) {
        transcripts_[tid].binMers[kmerClassID] += 1;
        // Since the transcript IDs are sorted we just have to check
        // if this id is different from the previous one
        if (tid != prevTID) { numDistinctTranscripts += 1.0; }
        prevTID = tid;
      }
      // Set the promiscuity and the set of transcripts for this kmer group
      kmerGroupPromiscuities_[kmerClassID] = numDistinctTranscripts;

      logKmerGroupCounts_[kmerClassID] = kmerGroupCounts_[kmerClassID] > 0 ? std::log(kmerGroupCounts_[kmerClassID]) : sailfish::math::LOG_0;
      //kmerGroupPromiscuities_[kmerClassID] = transcriptsForKmer_[kmerClassID].size();
    }
    cerr << "done\n";

    //promiscuousKmerCutoff_ = 50;
    for (auto kmerClassID : boost::irange(size_t{0}, numKmerClasses)) {
        if (kmerGroupPromiscuities_[kmerClassID] > promiscuousKmerCutoff_ ) {
            kmerGroupCounts_[kmerClassID] = 0;
        }
    }
    // By now we should have properly filled out the vairables
    // kmerGroupPromiscuities_
    // kmerGroupCounts_
    // transcripts_
  }


  /**
   * This function should be called before performing any optimization procedure.
   * It builds all of the necessary data-structures which are used during the transcript
   * quantification procedure.
   * @param  klutfname [The name of the file containing the kmer lookup table.]
   * @param  tlutfname [The name of the file containing the transcript lookup table.]
   */
    void initialize_(
        const std::string& klutfname,
        const std::string& tlutfname,
        const std::string& kmerEquivClassFname,
        const bool discardZeroCountKmers) {

        // So we can concisely identify each transcript
        TranscriptID transcriptIndex {0};

        size_t numTranscripts = transcriptGeneMap_.numTranscripts();
        size_t numKmers = readHash_.size();
        auto merSize = readHash_.kmerLength();

        size_t numActors = numThreads_;
        std::vector<std::thread> threads;

        transcripts_.resize(transcriptGeneMap_.numTranscripts());

        // Get the kmer look-up-table from file
        LUTTools::readKmerLUT(klutfname, transcriptsForKmer_);

/*        size_t numContainingTranscripts = transcriptsForKmer_[kid].size();
          assert(numContainingTranscripts > 0);
          if (numContainingTranscripts == 1 or
              std::all_of(transcriptsForKmer_[kid].begin(), transcriptsForKmer_[kid].end(), [this, kid](TranscriptID tid) -> bool {
                return tid == this->transcriptsForKmer_[kid][0]; } )) {
              uniquelyAnchoredTranscripts.insert(transcriptsForKmer_[kid][0]);
          }

        }
*/
        //std::cerr << "\nIn the original set, there were " << uniquelyAnchoredTranscripts.size() << " transcripts with unique nonzero anchors\n";

        std::cerr << "\n";
        //  collapseKmers_(isActiveKmer); // equiv-classes
        prepareCollapsedMaps_(kmerEquivClassFname, discardZeroCountKmers);


        // we have no k-mer-specific biases currently
        kmerGroupBiases_.resize(transcriptsForKmer_.size(), 1.0);

        // Get transcript lengths
        std::ifstream ifile(tlutfname, std::ios::binary);
        size_t numRecords {0};
        ifile.read(reinterpret_cast<char *>(&numRecords), sizeof(numRecords));

        std::cerr << "Transcript LUT contained " << numRecords << " records\n";
        for (auto i : boost::irange(size_t(0), numRecords)) {
            auto ti = LUTTools::readTranscriptInfo(ifile);
            // copy over the length, then we're done.
            auto& ts = transcripts_[ti->transcriptID];
            ts.length = ti->length;
            // would be length - k + 1, but we could have other Ns in the transcript
            ts.effectiveLength = ti->kmers.size();
            ts.isAnchored = false;
            ts.logInvEffectiveLength = (ts.effectiveLength > 0) ? std::log(1.0 / ts.effectiveLength) : sailfish::math::LOG_0;
        }
        ifile.close();
        // --- done ---

       // tbb::parallel_for( size_t(0), size_t(transcripts_.size()),
       //     [this]( size_t idx ) {
       //         auto& transcript = this->transcripts_[idx];
       //         transcript.effectiveLength = transcript.effectiveLength - transcript.binMers.size();
       //         for (auto binmer : transcript.binMers) {
       //          transcript.effectiveLength += this->_weight(binmer.first);
       //         }
       // });

        size_t numRes = 0;
        std::cerr << "\n\nRemoving duplicates from kmer transcript lists ... ";
        tbb::parallel_for(BlockedIndexRange(size_t(0), size_t(transcriptsForKmer_.size())),
            [&numRes, this](const BlockedIndexRange& range) -> void {
                for (auto idx = range.begin(); idx != range.end(); ++idx) {
                  auto& transcripts = this->transcriptsForKmer_[idx];
                  // should already be sorted -- extra check can be removed eventually
                  std::is_sorted(transcripts.begin(), transcripts.end());
                  // Uniqify the transcripts
                  auto it = std::unique(transcripts.begin(), transcripts.end());
                  transcripts.resize(std::distance(transcripts.begin(), it));
                  ++numRes;
                }
         });

         std::cerr << "done\n";

         std::cerr << "Computing kmer group promiscuity rates\n";
         /* -- done
         kmerGroupPromiscuities_.resize(transcriptsForKmer_.size());
         tbb::parallel_for( size_t{0}, kmerGroupPromiscuities_.size(),
            [this]( KmerID kid ) -> void { this->kmerGroupPromiscuities_[kid] = this->_weight(kid); }
         );
         */

        tbb::parallel_for(BlockedIndexRange(size_t{0}, transcripts_.size()),
          [&, this](const BlockedIndexRange& range) -> void {
            for (auto tid = range.begin(); tid != range.end(); ++tid) {
              auto& ti = this->transcripts_[tid];
              for (auto& binmer : ti.binMers) {
                if (binmer.second > promiscuousKmerCutoff_) {
                  ti.effectiveLength -= 1.0;
                }
              }
            }
        });

        /**
         * gene-promiscuous kmers can never be added to a transcript's counts, so
         * it's unfair to consider them in the transcripts effective length.
         */
        /*
        std::for_each( genePromiscuousKmers_.begin(), genePromiscuousKmers_.end(),
            [this]( KmerID kmerId ) -> void {
                for ( auto tid : transcriptsForKmer_[kmerId] ) {
                    transcripts_[tid].effectiveLength -= 1.0;
                }
            });
        */
        std::cerr << "done\n";

        //return mappedReads;
    }

    void _dumpCoverage( const boost::filesystem::path &cfname) {
        auto memberships = LUTTools::readKmerEquivClasses(kmerEquivClassFname_);

        size_t numTrans = transcripts_.size();
        size_t numProc = 0;
        std::ifstream kmerStructFile("zm_transcript.map", std::ios::in | std::ios::binary);
        std::unordered_map<std::string, std::vector<KmerID>> kstruct;
        uint64_t tlen{0};
        uint32_t nameLen{0};
        uint64_t tid{0};
        for (size_t i = 0; i < numTrans; ++i) {
            kmerStructFile.read(reinterpret_cast<char*>(&nameLen), sizeof(nameLen));
            std::string name(nameLen, ' ');
            kmerStructFile.read(reinterpret_cast<char*>(&name[0]), nameLen * sizeof(name[0]));
            kmerStructFile.read(reinterpret_cast<char*>(&tlen), sizeof(tlen));
            kstruct[name].resize(tlen);
            auto& tvec = kstruct[name];
            kmerStructFile.read(reinterpret_cast<char*>(&tvec[0]), tlen * sizeof(kstruct[name][0]));
        }

        kmerStructFile.close();

       std::ofstream ofile(cfname.native());

        ofile << "# numtranscripts_\n";
        ofile << "# transcript_name_{1} num_kmer_classes{1} class_1_count class_2_count ... class_{num_eq_clases}_count\n";
        ofile << "# ... \n";

        ofile << transcripts_.size() << "\n";

        std::cerr << "Dumping coverage statistics to " << cfname << "\n";


        tbb::concurrent_queue<StringPtr> covQueue;

        tbb::parallel_for(BlockedIndexRange(size_t{0}, transcripts_.size()),
             [this, &covQueue, &memberships, &kstruct] (const BlockedIndexRange& range) -> void {
                for (auto index = range.begin(); index != range.end(); ++index) {
                  auto& td = this->transcripts_[index];
                  const auto& name = this->transcriptGeneMap_.transcriptName(index);

                  std::stringstream ostream;
                  ostream << this->transcriptGeneMap_.transcriptName(index) << " " << kstruct[name].size();

                  std::map<uint64_t, double> kmerClassRelativeMass;

                  for (auto bm : kstruct[name]) {
                      auto kclass = memberships[bm];
                      double totalMass = 0.0;
                      for (auto tid : this->transcriptsForKmer_[kclass]) {
                          totalMass += this->transcripts_[tid].mean;
                      }
                      kmerClassRelativeMass[kclass] = (totalMass > 0.0) ? td.mean / totalMass : 0.0;
                  }

                  for (auto bm : kstruct[name]) {
                      auto kclass = memberships[bm];
                      auto count = readHash_.atIndex(bm);
                      ostream << " " << count * kmerClassRelativeMass[kclass];
                  }
                  ostream << "\n";
                  std::string* ostr = new std::string(ostream.str());
                  covQueue.push(ostr);
              }
            }
        );

        ez::ezETAProgressBar pb(transcripts_.size());
        pb.start();

        std::string* sptr = nullptr;
        while ( numProc < numTrans ) {
            while( covQueue.try_pop(sptr) ) {
                ofile << (*sptr);
                ++pb;
                ++numProc;
                delete sptr;
            }
        }

        ofile.close();

    }

    void filterByCoverage_( double coverageCutoff, std::vector<double>& means  ) {

        tbb::parallel_for(BlockedIndexRange(size_t{0}, transcripts_.size()),
                          [this, coverageCutoff, &means] (const BlockedIndexRange& range) -> void {
                              for (auto index = range.begin(); index != range.end(); ++index) {
                                  double numCovered{1.0};
                                  double totalNumKmers{1.0};
                                  bool unAnchored{true};
                                  const auto& td = this->transcripts_[index];

                                  for ( auto bm : td.binMers ) {
                                      double numKmersToCover = this->kmerGroupSizes_[bm.first];
                                      numCovered += (bm.second > 0.0) ? numKmersToCover : 0.0;
                                      totalNumKmers += numKmersToCover;
                                  }

                              if (unAnchored and ((numCovered / totalNumKmers) < coverageCutoff)) {
                                  means[index] = this->transcripts_[index].mean = 0.0;
                              }
                              }
                          }
                          );

    }

    void EMUpdate_( const std::vector<double>& meansIn, std::vector<double>& meansOut, bool accel) {
      assert(meansIn.size() == meansOut.size());

      auto reqNumJobs = transcriptsForKmer_.size();

      std::atomic<size_t> numJobs{0};
      std::atomic<size_t> completedJobs{0};

      // Print out our progress
      auto pbthread = std::thread(
        [&completedJobs, reqNumJobs]() -> bool {
          auto prevNumJobs = 0;
          ez::ezETAProgressBar show_progress(reqNumJobs);
          show_progress.start();
          while ( prevNumJobs < reqNumJobs ) {
            if ( prevNumJobs < completedJobs ) {
              show_progress += completedJobs - prevNumJobs;
            }
            prevNumJobs = completedJobs.load();
            boost::this_thread::sleep_for(boost::chrono::milliseconds(1));
          }
          if (!show_progress.isDone()) { show_progress.done(); }
          return true;
        });

      std::vector<double> rho(transcripts_.size(), 0.0);

      size_t numTranscripts = transcripts_.size();
      double priorAlpha = 0.01;
      double totalKmerCount = (priorAlpha * numTranscripts) + psum_(kmerGroupCounts_);
      double logAlpha0 = boost::math::digamma(totalKmerCount);
      tbb::parallel_for(BlockedIndexRange(size_t(0), size_t(transcripts_.size())),
          // for each kmer group
          [&completedJobs, &rho, logAlpha0, totalKmerCount, &meansIn, &meansOut, accel, reqNumJobs, this](const BlockedIndexRange& range) -> void {
            for (auto tid : boost::irange(range.begin(), range.end())) {
                auto& t = this->transcripts_[tid];
                double tMass = t.totalMass;
                double currCount = (accel) ? meansIn[tid] * totalKmerCount : tMass;
                if (currCount >= 1.0 and t.effectiveLength > 0) {
                    rho[tid] = boost::math::digamma(currCount) - logAlpha0 + t.logInvEffectiveLength;
                } else {
                    rho[tid] = sailfish::math::LOG_0;
                }
            }
      });

      tbb::parallel_for(BlockedIndexRange(size_t{0}, transcripts_.size()), [this](const BlockedIndexRange& range) -> void {
              for (auto tid : boost::irange(range.begin(), range.end())) { this->transcripts_[tid].resetMass();}}
      );

      //  E-Step : reassign the kmer group counts proportionally to each transcript
      tbb::parallel_for(BlockedIndexRange(size_t(0), size_t(transcriptsForKmer_.size())),
          // for each kmer group
          [&completedJobs, &meansIn, &rho, &meansOut, priorAlpha, reqNumJobs, this](const BlockedIndexRange& range) -> void {
            for (auto kid : boost::irange(range.begin(), range.end())) {
                // for each transcript containing this kmer group
                auto kmer = kid;
                auto& transcripts = this->transcriptsForKmer_[kmer];

                double totalMass = 0.0;
                //double logTotalMass = sailfish::math::LOG_0;
                /**
                 * Compute the total mass of all transcripts containing this k-mer
                 */
                for ( auto tid : transcripts ) {
                    auto& trans = this->transcripts_[tid];
                    //totalMass += meansIn[tid] / trans.effectiveLength;
                    if (rho[tid] != sailfish::math::LOG_0) {
                        totalMass += std::exp(rho[tid]);
                        //logTotalMass = sailfish::math::logAdd(logTotalMass, rho[tid]);
                    }
                }

                //double localEpsilon = 1e-15;
                //double norm = (totalMass >  localEpsilon) ? 1.0 / totalMass : 0.0;
                double norm = (totalMass >  sailfish::math::EPSILON) ? 1.0 / totalMass : 0.0;

                for ( auto tid : transcripts ) {
                  auto& trans = this->transcripts_[tid];
                  auto lastIndex = trans.binMers.size()  - 1;

                  if (trans.effectiveLength > 0) {
                      if (rho[tid] != sailfish::math::LOG_0) {
                          auto newMass = std::exp(rho[tid]) * norm * this->kmerGroupCounts_[kmer];
                          trans.addMass(newMass);
                      }  else {
                          trans.binMers[kmer] = 0;
                      }

                  }

                  // M-Step
                  // If we've seen all of the k-mers that appear in this transcript,
                  // then we can compute it's new estimated abundance
                  if (trans.updated++ == lastIndex) {
                      //trans.mean = meansOut[tid] = this->_computeMean(trans);
                      //meansOut[tid] = trans.totalMass;// / this->transcripts_[tid].effectiveLength;
                      meansOut[tid] = priorAlpha + trans.totalMass;
                      trans.updated.store(0);
                  }

                }

              ++completedJobs;
            } // for kid in range
          });

          // wait for all kmer groups to be processed
          pbthread.join();

          /*
          tbb::parallel_for(BlockedIndexRange(size_t(0), size_t(meansIn.size())),
                  // for each kmer group
                  [&completedJobs, &meansIn, &rho, &meansOut, priorAlpha, reqNumJobs, this](const BlockedIndexRange& range) -> void {
                      // M-Step
                      for (auto tid : boost::irange(range.begin(), range.end())) {
                            auto& trans = transcripts_[tid];
                            meansOut[tid] = priorAlpha + trans.totalMass;
                       }
                  });

          */
          // Make the output a proper probability vector
          normalize_(meansOut);
    }


public:
    /**
     * Construct the solver with the read and transcript hashes
     */
    CollapsedIterativeOptimizer( ReadHash &readHash, TranscriptGeneMap& transcriptGeneMap,
                                 BiasIndex& biasIndex, uint32_t numThreads ) :
                                 readHash_(readHash), merLen_(readHash.kmerLength()),
                                 transcriptGeneMap_(transcriptGeneMap), biasIndex_(biasIndex),
                                 numThreads_(numThreads) {}


    KmerQuantity optimize(const std::string& klutfname,
                          const std::string& tlutfname,
                          const std::string& kmerEquivClassFname,
                          size_t numIt,
                          double minMean,
                          double maxDelta) {


        using sailfish::math::LOG_1;
        using sailfish::math::LOG_0;
        using sailfish::math::logAdd;
        using sailfish::math::logSub;

        kmerEquivClassFname_ = kmerEquivClassFname;
        const bool discardZeroCountKmers = true;
        initialize_(klutfname, tlutfname, kmerEquivClassFname, discardZeroCountKmers);

        KmerQuantity globalError {0.0};
        bool done {false};
        std::atomic<size_t> numJobs {0};
        std::atomic<size_t> completedJobs {0};
        std::vector<KmerID> kmerList( transcriptsForKmer_.size(), 0 );
        size_t idx = 0;

        tbb::task_scheduler_init tbb_init(numThreads_);

        std::cerr << "Computing initial coverage estimates ... ";

        std::vector<double> means0(transcripts_.size(), 0.0);
        std::vector<double> means1(transcripts_.size(), 0.0);
        std::vector<double> means2(transcripts_.size(), 0.0);
        std::vector<double> meansPrime(transcripts_.size(), 0.0);

        std::vector<double> r(transcripts_.size(), 0.0);
        std::vector<double> sv2(transcripts_.size(), 0.0);
        std::vector<double> v(transcripts_.size(), 0.0);
        std::atomic<size_t> uniquelyAnchoredTranscripts{0};
        std::atomic<size_t> nonZeroTranscripts{0};

        loggedCounts_ = false;
        //logKmerGroupCounts_();

        // Compute the initial mean for each transcript
        tbb::parallel_for(BlockedIndexRange(size_t(0), size_t(transcriptGeneMap_.numTranscripts())),
        [this, &uniquelyAnchoredTranscripts, &nonZeroTranscripts, &means0](const BlockedIndexRange& range) -> void {
            for (auto tid = range.begin(); tid != range.end(); ++tid) {
              auto& transcriptData = this->transcripts_[tid];
              KmerQuantity total = 0.0;
              bool notTotallyPromiscuous{false};
              bool nonZero{false};

              transcriptData.countSpace = CountSpace::NonLogSpace;
              transcriptData.resetMass();

              for ( auto & kv : transcriptData.binMers ) {
                auto kmer = kv.first;
                if ( this->genePromiscuousKmers_.find(kmer) == this->genePromiscuousKmers_.end() ){
                    // count is the number of times kmer appears in transcript (tid)
                    auto count = kv.second;
                    kv.second = count * this->kmerGroupCounts_[kmer] * this->weight_(kmer);
                    transcriptData.addMass(kv.second);
                    if (kv.second > 0 and kmerGroupPromiscuities_[kmer] == 1) {
                        notTotallyPromiscuous = true;
                        transcriptData.isAnchored = true;
                    }
                    if (kv.second > 0) { nonZero = true; }
                }
              }
              means0[tid] = transcriptData.totalMass;
              if (notTotallyPromiscuous) { ++uniquelyAnchoredTranscripts; }
              if (nonZero) { ++nonZeroTranscripts; }
            }
        }
        );
        normalize_(means0);

        std::cerr << "\nThere were " << uniquelyAnchoredTranscripts.load() << " uniquely anchored transcripts\n";
        std::cerr << "There were " << nonZeroTranscripts.load() << " transcripts with at least one overlapping k-mer\n";
        std::cerr << "done\n";
        size_t outerIterations = 1;

        /**
         * Defaults for these values taken from the R implementation of
         * [SQUAREM](http://cran.r-project.org/web/packages/SQUAREM/index.html).
         */
        double minStep0, minStep, maxStep0, maxStep, mStep, nonMonotonicity;
        minStep0 = 1.0; minStep = 1.0;
        maxStep0 = 1.0; maxStep = 1.0;
        mStep = 4.0;
        nonMonotonicity = 1.0;

        double negLogLikelihoodOld = std::numeric_limits<double>::infinity();
        double negLogLikelihoodNew = std::numeric_limits<double>::infinity();
        bool accel{false};

        std::function<bool(std::vector<double>&, std::vector<double>&)> hasConverged;
        if (std::isfinite(maxDelta)) {
            hasConverged = [maxDelta, &accel, this] (std::vector<double>& v0, std::vector<double>& v1) -> bool {
                double minVal = 1e-7;
                auto relDiff = this->relAbsDiff_(v0, v1, minVal);
                double maxRelativeChange = *std::max_element( relDiff.begin(), relDiff.end() );
                std::cerr << "max relative change: " << maxRelativeChange << "\n";
                if (maxRelativeChange < 10.0*maxDelta) { accel = true; }// else { accel = true; }
                //accel = false;
                return maxRelativeChange < maxDelta;
            };
        } else {
            hasConverged = [&accel] (std::vector<double>& v0, std::vector<double>& v1) -> bool {
                std::cerr << "no data-driven convergence criterion specified\n";
                accel = true;
                return false;
            };
        }

        std::string clearline = "                                                                                \r\r";

        // Until we've reached the specified maximum number of iterations, or hit ourt
        // tolerance threshold
        for ( size_t iter = 0; iter < numIt; ++iter ) {
            std::string jumpBack = "\x1b[A";
            std::cerr << clearline << "SQUAREM iteraton [" << iter << " / max(" << numIt << ")]\n";
            jumpBack += "\x1b[A";

          // Theta_1 = EMUpdate(Theta_0)
          std::cerr << clearline << "1/3\n";
          EMUpdate_(means0, means1, accel);
          jumpBack += "\x1b[A\x1b[A";

          // Check for data-driven convergence criteria
          if (hasConverged(means0, means1)) {
              std::cerr << "convergence criteria met; terminating SQUAREM\n";
              break;
          } else if (!accel) {
              std::cerr << jumpBack;
              continue;
          }


          if (!std::isfinite(negLogLikelihoodOld)) {
            negLogLikelihoodOld = -expectedLogLikelihood_(means0);
          }

          // Theta_2 = EMUpdate(Theta_1)
          std::cerr << clearline << "2/3\n";
          EMUpdate_(means1, means2, accel);
          jumpBack += "\x1b[A\x1b[A";

          double delta = pabsdiff_(means1, means2);
          std::cerr << clearline << "delta = " << std::setw(6) << std::setfill(' ') << delta << "\n";
          jumpBack += "\x1b[A";

          // r = Theta_1 - Theta_0
          // v = (Theta_2 - Theta_1) - r
          tbb::parallel_for(BlockedIndexRange(size_t(0), transcripts_.size()),
              [&means0, &means1, &means2, &r, &v](const BlockedIndexRange& range) -> void {
              for (auto tid = range.begin(); tid != range.end(); ++tid) {
                r[tid] = means1[tid] - means0[tid];
                v[tid] = (means2[tid] - means1[tid]) - r[tid];
              }
            }
          );

          double rNorm = std::sqrt(dotProd_(r,r));
          double vNorm = std::sqrt(dotProd_(v,v));
          double alphaS = rNorm / vNorm;

          alphaS = std::max(minStep, std::min(maxStep, alphaS));

          tbb::parallel_for(BlockedIndexRange(size_t(0), transcripts_.size()),
            [&r, &v, alphaS, &means0, &meansPrime](const BlockedIndexRange& range) -> void {
              for (auto tid = range.begin(); tid != range.end(); ++tid) {
               // Looking into Nick Bray's problem
               meansPrime[tid] = std::max(means0[tid] * 0.1, means0[tid] + 2*alphaS*r[tid] + (alphaS*alphaS)*v[tid]);
               //meansPrime[tid] = std::max(0.0, means0[tid] + 2*alphaS*r[tid] + (alphaS*alphaS)*v[tid]);
              }
            }
          );

          // Stabilization step
          if (std::abs(alphaS - 1.0) > 0.01) {
              std::cerr << clearline << "alpha = " <<  std::setw(6) << std::setfill(' ') << alphaS << ". ";
            std::cerr << "Performing a stabilization step.\n";
            EMUpdate_(meansPrime, meansPrime, accel);
            jumpBack += "\x1b[A\x1b[A";
          }


          /** Check for an error in meansPrime **/
          /*if (hasNan_(meansPrime)) {
              std::swap(meansPrime, means2);
              negLogLikelihoodNew = -expectedLogLikelihood_(meansPrime);
              if (alphaS == maxStep) { maxStep = std::max(maxStep0, maxStep / mStep); }
              alphaS = 1.0;
          } else {
          */
          /** If there is **/
          if (std::isfinite(nonMonotonicity)) {
            negLogLikelihoodNew = -expectedLogLikelihood_(meansPrime);
          } else {
            negLogLikelihoodNew = negLogLikelihoodOld;
          }

          if (negLogLikelihoodNew > negLogLikelihoodOld + nonMonotonicity) {
            std::swap(meansPrime, means2);
            negLogLikelihoodNew = -expectedLogLikelihood_(meansPrime);
            if (alphaS == maxStep) { maxStep = std::max(maxStep0, maxStep/mStep); }
            alphaS = 1.0;
          }
          std::cerr << clearline << "alpha = " << std::setw(6) << std::setfill(' ') << alphaS << ", ";
          //}

          if (alphaS == maxStep) { maxStep = mStep * maxStep; }
          if (minStep < 0 and alphaS == minStep) { minStep = mStep * minStep; }
          std::swap(meansPrime, means0);

          if (!std::isnan(negLogLikelihoodNew)) {
            std::cerr << "negLogLikelihood = " << std::setw(6) << std::setfill(' ') << negLogLikelihoodNew << "\n";
            jumpBack += "\x1b[A";
            negLogLikelihoodOld = negLogLikelihoodNew;

          }

          // If it's not the final iteration, then jump back and clear the console.
          if (iter < numIt - 1) { std::cerr << jumpBack; }

        }
        KmerQuantity q{0.0};
        return q;
    }


    void writeAbundances(const boost::filesystem::path& outputFilePath,
                         const std::string& headerLines,
                         double minAbundance,
                         bool haveCI) {

        haveCI = false;

        std::cerr << "Writing output\n";
        ez::ezETAProgressBar pb(transcripts_.size());
        pb.start();

        auto estimatedGroupTotal = psum_(kmerGroupCounts_);
        auto totalNumKmers = readHash_.totalLength() * (std::ceil(readHash_.averageLength()) - readHash_.kmerLength() + 1);
        std::vector<KmerQuantity> fracTran(transcripts_.size(), 0.0);
        std::vector<KmerQuantity> fracTranLow(transcripts_.size(), 0.0);
        std::vector<KmerQuantity> fracTranHigh(transcripts_.size(), 0.0);

        std::vector<KmerQuantity> numKmersFromTranscript(transcripts_.size(), 0.0);

        // Compute nucleotide fraction (\nu_i in RSEM)
        tbb::parallel_for(BlockedIndexRange(size_t(0), size_t(transcripts_.size())),
           [this, &numKmersFromTranscript, &fracTran, totalNumKmers, estimatedGroupTotal](const BlockedIndexRange& range) -> void {
            for (auto tid = range.begin(); tid != range.end(); ++tid) {
              auto& ts = this->transcripts_[tid];
              //fracTran[tid] = this->_computeMean(ts);

              // EM
              //numKmersFromTranscript[tid] = this->_computeSum(ts);
              if (ts.totalMass == sailfish::math::LOG_0) {
                  ts.totalMass = 0.0;
              }
              numKmersFromTranscript[tid] = ts.totalMass;

              // VB
              //numKmersFromTranscript[tid] = ts.mean * estimatedGroupTotal;//this->_computeSum(ts);
              //fracTran[tid] = this->_computeClampedMean( ts );
            }
        });
        // noise transcript
        // fracTran[transcripts_.size()] = totalNumKmers - estimatedGroupTotal;

        // Normalize to obtain the fraction of kmers that originated from each transcript
        std::vector<KmerQuantity> fracNuc(numKmersFromTranscript);
        normalize_(fracNuc);

        // Compute transcript fraction (\tau_i in RSEM)
        tbb::parallel_for(BlockedIndexRange(size_t(0), transcripts_.size()),
          [&, this](const BlockedIndexRange& range) -> void {
            for (auto i = range.begin(); i != range.end(); ++i) {
              auto& ts = transcripts_[i];

              //fracTran[i] = this->_computeMean(ts);
              fracTran[i] = (ts.effectiveLength > 0.0) ? fracNuc[i] / ts.effectiveLength : 0.0;
              if (!std::isfinite(fracTran[i])) {
                  std::cerr << "Transcript had non-finite transcript fraction " << fracTran[i] << ",";
                  std::cerr << " fracNuc = " << fracNuc[i] << ", effLength = " << ts.effectiveLength << "\n";
              }

              ts.mean = numKmersFromTranscript[i] / ts.effectiveLength;
              //fracTran[i] = (ts.effectiveLength > 0) ? ts.mean / ts.effectiveLength : 0.0;
              fracTranLow[i] = (ts.effectiveLength > 0) ? ts.fracLow / ts.effectiveLength : 0.0;
              fracTranHigh[i] = (ts.effectiveLength > 0) ? ts.fracHigh / ts.effectiveLength : 0.0;

            }
          });
        //normalize_(fracTran);

        auto fracNorm = 1.0 / psum_(fracTran);
        tbb::parallel_for(BlockedIndexRange(size_t(0), transcripts_.size()),
                          [&, fracNorm](const BlockedIndexRange& range) -> void {
                              for (auto i = range.begin(); i != range.end(); ++i) {
                                  auto& ts = transcripts_[i];
                                  fracTran[i] *= fracNorm;
                                  fracTranLow[i] *= fracNorm;
                                  fracTranHigh[i] *= fracNorm;
                              }
        });

        auto writeCoverageInfo = false;
        if ( writeCoverageInfo ) {
            boost::filesystem::path covPath = outputFilePath.parent_path();
            covPath /= "equivClassCoverage.txt";
            _dumpCoverage( covPath );
        }

        std::ofstream ofile( outputFilePath.string() );
        size_t index = 0;
        double million = std::pow(10.0, 6);
        double billion = std::pow(10.0, 9);
        double estimatedReadLength = readHash_.averageLength();
        double kmersPerRead = ((estimatedReadLength - merLen_) + 1);
        double totalMappedKmers = std::accumulate(numKmersFromTranscript.begin(),
                                                  numKmersFromTranscript.end(),
                                                  0.0);
        double totalMappedReads = totalMappedKmers / kmersPerRead;

        ofile << headerLines;

        ofile << "# " <<
            "Transcript" << '\t' <<
            "Length" << '\t' <<
            "TPM" << '\t' <<
            "RPKM" << '\t' <<
            "KPKM" << '\t' <<
            "EstimatedNumKmers" << '\t' <<
            "EstimatedNumReads";

        if (haveCI) {
            ofile << '\t' << "TPM_LOW" << '\t' << "TPM_HIGH";
        }

        ofile << '\n';

        for ( auto i : boost::irange(size_t{0}, transcripts_.size()) ) {
          auto& ts = transcripts_[i];
          // expected # of kmers coming from transcript i

          // auto ci = estimatedGroupTotal * fracNuc[i];
          auto ci = numKmersFromTranscript[i];

          // expected # of reads coming from transcript i
          auto ri = ci / kmersPerRead;
          double effectiveLength = ts.effectiveLength + (merLen_ - 1) - std::floor(estimatedReadLength) + 1;//ts.length - std::floor(estimatedReadLength) + 1;
          double effectiveLengthKmer = ts.effectiveLength;//ts.length - merLen_ + 1;

          auto kpkm = (effectiveLengthKmer > 0) ?
            (ci * billion) / (effectiveLengthKmer * totalMappedKmers) : 0.0;
          kpkm = (kpkm < minAbundance) ? 0.0 : kpkm;

          auto rpkm = (effectiveLength > 0) ?
            (ri * billion) / (effectiveLength * totalMappedReads) : 0.0;
          rpkm = (kpkm < minAbundance) ? 0.0 : rpkm;

          // PREVIOUS
          // auto ci = estimatedGroupTotal * fracNuc[i];
          /*
          auto kpkm = (effectiveLengthKmer > 0) ?
            (billion * (fracNuc[i] / ts.effectiveLength)) : 0.0;
          kpkm = (kpkm < minAbundance) ? 0.0 : kpkm;
          auto rpkm = (effectiveLength > 0) ?
          (billion * (fracNuc[i] / ts.effectiveLength)) : 0.0;
          rpkm = (kpkm < minAbundance) ? 0.0 : rpkm;

          */


          auto tpm = fracTran[i] * million;
          tpm = (kpkm < minAbundance) ? 0.0 : tpm;

          ofile << transcriptGeneMap_.transcriptName(index) <<
                   '\t' << ts.length << '\t' <<
                   tpm << '\t' <<
                   rpkm << '\t' <<
                   kpkm << '\t' <<
                   ci << '\t' <<
                   ri;
          if (haveCI) {
              ofile << '\t' <<
                  fracTranLow[i] * million << '\t' <<
                  fracTranHigh[i] * million;
          }
          ofile << '\n';

          ++index;
          ++pb;
        }
        ofile.close();
    }


    /**
     *  Instead of using the EM (SQUAREM) algorithm, infer the posterior distribution
     *  using Streaming, Distributed, Asynchronous (SDA) variational Bayes
     */
    KmerQuantity optimizeVB(const std::string& klutfname,
                          const std::string& tlutfname,
                          const std::string& kmerEquivClassFname,
                          size_t numIt,
                          double minMean,
                          double maxDelta) {


         kmerEquivClassFname_ = kmerEquivClassFname;
        // Prepare the necessary structures
        const bool discardZeroCountKmers = false;
        initialize_(klutfname, tlutfname, kmerEquivClassFname, discardZeroCountKmers);

        const size_t numTranscripts = transcripts_.size();
        const size_t numKmers = transcriptsForKmer_.size();

        // Set the appropriate, user-specified, convergence criteria
        std::function<bool(std::vector<double>&, std::vector<double>&)> hasConverged;
        if (std::isfinite(maxDelta)) {
            hasConverged = [maxDelta, this] (std::vector<double>& v0, std::vector<double>& v1) -> bool {
                double maxVal = *std::max_element(v1.begin(), v1.end());
                std::cerr << "maxVal: " << maxVal << "\n";
                double minVal = 1e-7;
                auto relDiff = this->relAbsDiff_(v0, v1, minVal);
                double maxRelativeChange = *std::max_element( relDiff.begin(), relDiff.end() );
                std::cerr << "max relative change: " << maxRelativeChange << "\n";
                return maxRelativeChange < maxDelta;
            };
        } else {
            hasConverged = [] (std::vector<double>& v0, std::vector<double>& v1) -> bool {
                std::cerr << "no data-driven convergence criterion specified\n";
                return false;
            };
        }


        // We'll use these to check convergence
        std::vector<double> meansOld(transcripts_.size(), 0.0);
        std::vector<double> meansNew(transcripts_.size(), 0.0);


        std::atomic<size_t> uniquelyAnchoredTranscripts{0};
        std::atomic<size_t> nonZeroTranscripts{0};

        const double DirichletPriorAlpha = 0.1;
        std::vector<double> priorAlphas(transcripts_.size(), 0.0);
        std::vector<double> posteriorAlphas(transcripts_.size(), 0.0);
        /*
        std::vector<double> kmerWeights(numKmers, 0.0);
        for (auto tid : boost::irange(size_t{0}, numTranscripts)) {
            auto& transcriptData = transcripts_[tid];
            for ( auto & kv : transcriptData.binMers ) {
                kmerWeights[kv.first] += kv.second * weight_(kv.first);
            }
        }
        for (auto kw : kmerWeights) {
            if (std::abs(kw - 1.0) > 1e-3) {
                std::cerr << "Kmer had weight " << kw << "\n";
            }
        }
        */

        // Compute the initial mean for each transcript
        tbb::parallel_for(BlockedIndexRange(size_t(0), numTranscripts),
                          [this, DirichletPriorAlpha, &uniquelyAnchoredTranscripts, &nonZeroTranscripts, &meansOld, &posteriorAlphas](const BlockedIndexRange& range) -> void {
            for (auto tid = range.begin(); tid != range.end(); ++tid) {
              auto& transcriptData = this->transcripts_[tid];
              KmerQuantity total = 0.0;

              bool notTotallyPromiscuous{false};
              bool nonZero{false};

              auto effLength = transcriptData.effectiveLength;
              for ( auto & kv : transcriptData.binMers ) {
                auto kmer = kv.first;

                // count is the number of times occurrences of the k-mer in this transcript
                auto count = kv.second;

                // weight(kmer) = 1 / total # occurences
                // count = # occurences in this transcript
                // kmerGroupCounts(kmer) = # of observations of kmer in the read set
                // The weight attributed to this transcript (kv.second) is:
                // count * weight(kmer) * kmerGroupCounts(kmer) = (# occurrences in ts / total # occurrences) * total count
                kv.second = (effLength  > 0.0) ?
                    (count * this->kmerGroupCounts_[kmer] * this->weight_(kmer)) :
                    0.0;

                if (kv.second > 0.0 and kmerGroupPromiscuities_[kmer] == 1) { notTotallyPromiscuous = true; }
                if (kv.second > 0.0) { nonZero = true; }
              }

              //transcriptData.mean = this->_computeMean(transcriptData);
              auto weightedKmerSum = this->_computeSum(transcriptData);
              // Set \alpha_t = \alpha_0 + \mu_t
              posteriorAlphas[tid] = DirichletPriorAlpha ;//+ weightedKmerSum;
              if (notTotallyPromiscuous) { ++uniquelyAnchoredTranscripts; }
              if (nonZero) { ++nonZeroTranscripts; }
            }
        }
        );

        auto posteriorAlphaSum = psum_(posteriorAlphas);
        std::cerr << "posteriorAlphaSum : " << posteriorAlphaSum << " [before iterations]\n";
        for (size_t i : boost::irange({0}, numTranscripts)) {
            meansOld[i] = posteriorAlphas[i] / posteriorAlphaSum;
        }

        std::vector<double> expctedLogThetas(transcripts_.size(), 0.0);
        std::vector<double> logRho(transcripts_.size(), -std::numeric_limits<double>::infinity());
        std::string clearline = "                                                                                \r\r";


        for (size_t currIt : boost::irange({0}, numIt)) {
            std::string jumpBack = "\x1b[A";
            std::cerr << clearline << "iteration : " << currIt << "\n";
            // log rho_{ntsoa} = E_{theta}[log theta_t] + log P(S_n | T_n) + [other terms sum to 0]
            // E_{theta}[log theta_t] = digamma(alpha_t) - digamma(\sum_{t'} alpha_{t'})
            double sumAlpha = psum_(posteriorAlphas);
            double digammaSumAlpha = boost::math::digamma(sumAlpha);
            for (size_t i : boost::irange({0}, numTranscripts)) {
                auto& ts = transcripts_[i];
                //double digammaAlphaT = boost::math::digamma(posteriorAlphas[i]);

                logRho[i] = (ts.effectiveLength > 0 and posteriorAlphas[i] > 0.0) ?
                    (boost::math::digamma(posteriorAlphas[i]) - digammaSumAlpha) + std::log(1.0 / ts.effectiveLength) :
                    -std::numeric_limits<double>::infinity();
            }

            std::atomic<size_t> numUpdated{0};

            //  E-Step : reassign the kmer group counts proportionally to each transcript
            tbb::parallel_for(BlockedIndexRange(size_t(0), numKmers),
                 // for each kmer group
                 [&meansOld, &meansNew, &logRho, &posteriorAlphas, &numUpdated, DirichletPriorAlpha, this](const BlockedIndexRange& range) -> void {
                     for (auto kid : boost::irange(range.begin(), range.end())) {
                         auto kmer = kid;
                         // for each transcript containing this kmer group
                         auto& transcripts = this->transcriptsForKmer_[kmer];

                         double totalMass = 0.0;
                         for ( auto tid : transcripts ) {
                             auto el = this->transcripts_[tid].effectiveLength;
                             double rho = (el > 0) ? std::exp(logRho[tid]) : 0.0;
                             totalMass += rho;
                         }

                         double norm = (totalMass > 0.0) ? (1.0 / totalMass) : 0.0;
                         for ( auto tid : transcripts ) {
                             auto& trans = this->transcripts_[tid];
                             auto lastIndex = trans.binMers.size()  - 1;
                             auto el = trans.effectiveLength;
                             double rho = (el > 0) ? std::exp(logRho[tid]) : 0.0;
                             trans.binMers[kmer] = rho * norm *
                                 kmerGroupBiases_[kmer] * this->kmerGroupCounts_[kmer];

                             // If we've seen all of the k-mers that appear in this transcript,
                             // then we can compute it's new estimated abundance
                             if (trans.updated++ == lastIndex) {
                                 ++numUpdated;

                                 // We're folding the length term into the binMer weights now
                                 // should this computeSum be a computeMean?
                                 //trans.mean = this->_computeSum(trans);//this->_computeMean(trans);
                                 auto sumWeightedReadMass = this->_computeSum(trans);

                                 // with filter
                                 // posteriorAlphas[tid] = DirichletPriorAlpha + sumWeightedReadMass;
                                 posteriorAlphas[tid] = (posteriorAlphas[tid] > 0.0) ? DirichletPriorAlpha + sumWeightedReadMass : 0.0;

                                 // without filter
                                 //posteriorAlphas[tid] = DirichletPriorAlpha + trans.mean;
                                 if (std::isnan(posteriorAlphas[tid])) {
                                     std::cerr << "posterior is nan at " << tid << "\n";
                                     std::cerr << "mean: " << trans.mean << "\n";
                                     std::cerr << "norm: " << norm << "\n";
                                 }
                                 trans.updated.store(0);
                             }
                         }
                     } // for kid in range
           });

            std::cerr << clearline << "numUpdated = " << numUpdated << " : numTranscripts = " << numTranscripts << "\n";
            jumpBack += "\x1b[A";

            auto posteriorAlphaSum = psum_(posteriorAlphas);
            std::cerr << clearline << "posterior alpha sum: " << posteriorAlphaSum << "\n";
            jumpBack += "\x1b[A";
            for (size_t i : boost::irange({0}, numTranscripts)) {
                meansNew[i] = posteriorAlphas[i] / posteriorAlphaSum;
                transcripts_[i].mean = meansNew[i];
            }

            // Check for data-driven convergence criteria
            if (hasConverged(meansOld, meansNew)) {
              std::cerr << clearline << "convergence criteria met; terminating VB\n";
              break;
            }
            jumpBack += "\x1b[A";

            std::swap(meansNew, meansOld);

            bool doFilter{false};
            double itFrac = currIt / static_cast<double>(numIt);
            double cutoff = std::min(1.0, std::exp(itFrac - 0.3));
            if (doFilter and currIt > 0 and currIt % 100 == 0) {
                // cutoff =  e^{x^2 - 0.5}
                filterByCoverage_(cutoff, posteriorAlphas);
            }

            std::cerr << clearline << "cutoff = " << cutoff << "\n";
            jumpBack += "\x1b[A";
            std::cerr << jumpBack;
        }

        // Compute the confidence intervals for the expression values
        // The marginal for each transcript fraction takes the form of a Beta distribution
        // we compute the confidence intervals by the appropriate quantiles of the distribution
        double sumOfAlphas = psum_(posteriorAlphas);
        tbb::parallel_for(BlockedIndexRange(size_t(0), size_t(transcriptGeneMap_.numTranscripts())),
             [this, &posteriorAlphas, sumOfAlphas](const BlockedIndexRange& range) -> void {
                    for (auto tid = range.begin(); tid != range.end(); ++tid) {
                        auto& transcriptData = this->transcripts_[tid];
                        auto alphaT = posteriorAlphas[tid];
                        if (alphaT > 0) {
                            boost::math::beta_distribution<double> marginalDist(alphaT, sumOfAlphas - alphaT);
                            transcriptData.mean = boost::math::mean(marginalDist);
                            transcriptData.fracLow = boost::math::quantile(marginalDist, 0.025);
                            transcriptData.fracHigh = boost::math::quantile(marginalDist, 0.975);
                        } else {
                            transcriptData.mean = transcriptData.fracLow = transcriptData.fracHigh = 0.0;
                        }
                    }
        });
        KmerQuantity q{0.0};
        return q;
    }

    void applyCoverageFilter_(const boost::filesystem::path& estIn,
                              const boost::filesystem::path& transcriptKmerMap,
                              const boost::filesystem::path& filteredEstOut,
                              double minAbundance) {

        using boost::tokenizer;

        auto memberships = LUTTools::readKmerEquivClasses(kmerEquivClassFname_);
        size_t numTrans = transcripts_.size();
        size_t numProc = 0;
        std::ifstream kmerStructFile(transcriptKmerMap.native(), std::ios::in | std::ios::binary);
        std::unordered_map<std::string, std::vector<KmerID>> kstruct;
        uint64_t tlen{0};
        uint32_t nameLen{0};
        uint64_t tid{0};
        for (size_t i = 0; i < numTrans; ++i) {
            kmerStructFile.read(reinterpret_cast<char*>(&nameLen), sizeof(nameLen));
            std::string name(nameLen, ' ');
            kmerStructFile.read(reinterpret_cast<char*>(&name[0]), nameLen * sizeof(name[0]));
            kmerStructFile.read(reinterpret_cast<char*>(&tlen), sizeof(tlen));
            kstruct[name].resize(tlen);
            auto& tvec = kstruct[name];
            kmerStructFile.read(reinterpret_cast<char*>(&tvec[0]), tlen * sizeof(kstruct[name][0]));
        }

        kmerStructFile.close();

        // Make a map of names to indices
        std::unordered_map<std::string, size_t> nameToID;
        auto indices = boost::irange({0}, transcripts_.size());
        std::for_each(indices.begin(), indices.end(),
                      [this, &nameToID](size_t index) -> void {
                          const auto& name = this->transcriptGeneMap_.transcriptName(index);
                          nameToID[name] = index;
                      });
        /*
        std::unordered_map<size_t, size_t> uniqueMap;
        for (auto& kv : kstruct) {
            auto& name = kv.first;
            auto id = nameToID[name];
            auto& kmers = kv.second;
            for (auto k : kmers) {
                auto r = uniqueMap.find(k);
                if (r == uniqueMap.end()) {
                    uniqueMap[k] = id;
                } else {
                    if (r->second != id) {
                        r->second = std::numeric_limits<size_t>::max();
                    }
                }
            }
        }
        */
        std::ifstream ifile(estIn.native());
        std::ofstream ofile(filteredEstOut.native());


        //double thresh{1.14};
        auto numReads = readHash_.numLengths();
        auto readLength = readHash_.averageLength();
        auto kmerLength = readHash_.kmerLength();
        double kmersPerRead = static_cast<double>(readLength - kmerLength) + 1;

        double numMappedKmers = psum_(kmerGroupCounts_);
        double mappingRate = numMappedKmers / static_cast<double>(kmersPerRead * numReads);
        double numMillionReads = numReads / 1000000.0;
        double perKilobase = 1.0 / 1000.0;
        //double thresh = (kmersPerRead * mappingRate * numMillionReads) * 0.1 * perKilobase;
        double numMillionKmers = kmersPerRead * numReads / 1000000.0;
        double thresh = (mappingRate * numMillionKmers) * 10.0 * minAbundance * perKilobase;

        std::cerr << "Mapping rate = " << mappingRate << "\n";
        std::cerr << "Threshold = " << thresh << "\n";

        auto sufficientCoverage = [&, this](const std::string& transcriptName, double expression, double threshold) -> bool {
            size_t numProperlyCovered{0};
            size_t index = nameToID[transcriptName];
            auto& td = this->transcripts_[index];

            td.isAnchored = false;
            /*
            size_t trivCount = 9;
            */
            std::vector<double> masses(kstruct[transcriptName].size(), 0.0);
            size_t ind{0};
            for (auto bm : kstruct[transcriptName]) {
                auto kclass = memberships[bm];
                double totalMass = 0.0;
                for (auto tid : this->transcriptsForKmer_[kclass]) {
                    auto& ts = this->transcripts_[tid];
                    totalMass += ts.totalMass;
                }
                auto count = readHash_.atIndex(bm);
                auto relMass = (totalMass > 0.0) ? td.totalMass / totalMass : 0.0;

                // if (uniqueMap[bm] == index and count > 0) { //}and this->kmerGroupCounts_[kclass] > trivCount ) {
                //     //if (this->kmerGroupPromiscuities_[kclass] == 1 and count > 1000 ) { //}and this->kmerGroupCounts_[kclass] > trivCount ) {
                //     //std::cerr << "Anchored: " << transcriptName << ", expression =" << expression << ", count = " << count << "\n";
                //     td.isAnchored = true;
                // }
                // if (this->kmerGroupPromiscuities_[kclass] == 1 and count > 0)  {
                //     if (uniqueMap[bm] != index) { //}and this->kmerGroupCounts_[kclass] > trivCount ) {
                //         std::cerr << "unique map says that kmer is not unique to " << transcriptName << " but kmerGroupPromiscuities says it is!\n";
                //     }
                // }
                masses[ind] = relMass * count;
                numProperlyCovered += (masses[ind] > thresh) ? 1 : 0;
                ++ind;
            }

            if (td.isAnchored) { return true; }
            if (masses.size() == 0) { return false; }
            size_t n = (masses.size() / 2);
            std::nth_element(masses.begin(), masses.begin()+n, masses.end());
            double rmedian = (masses.size() % 2 == 0) ? (masses[n-1] + masses[n]) * 0.5 : masses[n];
            return rmedian > threshold;
            //return (numProperlyCovered > (mappingRate * masses.size()));
        };

        size_t lineNum{0};
        boost::char_separator<char> sep("\t ");
        for (std::string line; std::getline(ifile, line); ) {
            // If this is a comment line, then pass it through directly
            // to the output file.
            if (line.front() == '#') {
                ofile << line << "\n";
                continue;
            }

            std::vector<std::string> toks;
            tokenizer<boost::char_separator<char>> tok(line, sep);
            std::for_each(tok.begin(), tok.end(), [&toks](const std::string& s) -> void { toks.push_back(s); });
            // VERSION DEPENDENT
            double expression = boost::lexical_cast<double>(toks[4]);
            if (expression <= minAbundance or sufficientCoverage(toks[0], expression, thresh)) {
                ofile << line << "\n";
            } else {
                ofile << toks[0] << '\t'  << toks[1];
                for (auto i : boost::irange({2}, toks.size())) { ofile << "\t0.0"; }
                ofile << '\n';
           }
           if (lineNum % 10000 == 0) {
               std::cerr << "Applied coverage filter to: " << lineNum << " transcripts\r\r";
           }
           ++lineNum;
        }
        std::cerr << "\n";

        ofile.close();

    }
};

#endif // ITERATIVE_OPTIMIZER_HPP
