#include <vector>
#include <unordered_map>
#include <atomic>

#include "tbb/task_scheduler_init.h"
#include "tbb/parallel_for.h"
#include "tbb/parallel_for_each.h"
#include "tbb/parallel_reduce.h"
#include "tbb/blocked_range.h"
#include "tbb/partitioner.h"
#include "concurrentqueue.h"

#include <boost/math/special_functions/digamma.hpp>

// C++ string formatting library
#include "spdlog/details/format.h"

#include "cuckoohash_map.hh"
#include "Eigen/Dense"
#include "Eigen/Sparse"

#include "CollapsedEMOptimizer.hpp"
#include "Transcript.hpp"
#include "TranscriptGroup.hpp"
#include "SailfishMath.hpp"
#include "ReadExperiment.hpp"
#include "BootstrapWriter.hpp"
#include "MultinomialSampler.hpp"


typedef Eigen::SparseMatrix<uint64_t> eqMat ;

void clusterTranscripts(std::vector<std::pair<const TranscriptGroup, TGValue>>& eqVec,
        std::vector<Transcript>& transcripts){


    // make a matrix out of it
    // iterate over equivalence classes make matrix
    size_t numEqClasses = eqVec.size();
    size_t numTranscripts = transcripts.size();
    //eqMat.resize(numEqClasses,numTranscripts);

    for (size_t eqID = 0; eqID < numEqClasses; ++eqID) {
        auto& kv = eqVec[eqID];
        uint64_t count = kv.second.count;
        const TranscriptGroup& tgroup = kv.first;
        //for each transcript in this class
        if(tgroup.valid){
            const std::vector<uint32_t>& txps = tgroup.txps;
            size_t groupSize = txps.size();
            //iterate
            if (BOOST_LIKELY(groupSize > 1)) {
                for (size_t i = 0; i < groupSize; ++i) {
                    // column = tid
                    // row = eqID
                    auto tid = txps[i] ;
                    eqMat.insert(eqID,tid) = count ;

                }
            }else{
                eqMat.insert(eqID,tid) = count ;
            }
        }
    }

    cout<<"Number of nonzero elelments: "<<sm1.nonZeros();

}
