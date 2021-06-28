#pragma once

/*
 * CollapsedHaps.hpp
 *
 *  Created on: Jun 12, 2021
 *      Author: nick
 */


#include "elucidator/common.h"
#include "elucidator/objects/BioDataObject.h"

#include <njhseq/concurrency/pools/AlignerPool.hpp>

#include <njhseq/objects/Meta/MultipleGroupMetaData.hpp>

namespace njhseq {



class CollapsedHaps{
public:

	std::vector<std::shared_ptr<seqInfo>> seqs_;
	std::vector<std::unordered_set<std::string>> names_;

	std::unordered_map<std::string, uint32_t> subNamesToMainSeqPos_;

	bool verbose_{false};

	void setSubNamesToMainSeqPos();

	void revCompSeqs();

	void setFrequencies(uint32_t total);

	void setFrequencies();

	uint32_t getTotalHapCount() const; /**< The total number of input haplotypes */
	uint32_t getTotalUniqueHapCount() const; /**< the total number of unique haplotypes */
	std::vector<uint32_t> getReadLenVec() const;

	std::unordered_map<uint32_t, uint32_t> getReadLenMap() const;

	static CollapsedHaps readInReads(const SeqIOOptions & inOpts,
			std::unique_ptr<MultipleGroupMetaData> meta = nullptr,
			std::unordered_map<std::string, std::string> metaValuesToAvoid = std::unordered_map<std::string, std::string>{});

	std::vector<comparison> getCompsAgainstRef(const seqInfo & refSeq, aligner & alignerObj, uint32_t numThreads = 1) const;
	std::vector<std::vector<comparison>> getPairwiseComps(aligner & alignerObj, uint32_t numThreads = 1) const;

	struct AvgPairwiseMeasures{
		double avgPercentId {0};
		double avgNumOfDiffs {0};
	};

	AvgPairwiseMeasures getAvgPairwiseMeasures(const std::vector<std::vector<comparison>> & allComps) const;

};



}  // namespace njhseq




