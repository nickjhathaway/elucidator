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


namespace njhseq {



class CollapsedHaps{
public:

	std::vector<std::shared_ptr<seqInfo>> seqs_;
	std::vector<std::unordered_set<std::string>> names_;

	std::unordered_map<std::string, uint32_t> subNamesToMainSeqPos_;

	void setSubNamesToMainSeqPos();

	void revCompSeqs();

	void setFrequencies(uint32_t total);

	void setFrequencies();

	uint32_t getTotalHapCount() const;

	std::vector<uint32_t> getReadLenVec() const;

	std::unordered_map<uint32_t, uint32_t> getReadLenMap() const;

	static CollapsedHaps readInReads(const SeqIOOptions & inOpts);

	std::vector<comparison> getCompsAgainstRef(const seqInfo & refSeq, aligner & alignerObj, uint32_t numThreads = 1) const;
};




}  // namespace njhseq




