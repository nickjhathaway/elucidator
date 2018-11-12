#pragma once
/*
 * BamCountingUtils.hpp
 *
 *  Created on: Mar 2, 2016
 *      Author: nick
 */

#include <njhseq/BamToolsUtils/BamCountExtractStats.hpp>

#include "elucidator/common.h"
#include "elucidator/objects/counters/RefCounting.h"


namespace njhseq {



void increaseCountsForAln(const BamTools::BamAlignment & bAln,
		const std::string & rName, RefCounter & counter, aligner & alignerObj,
		uint32_t qualCutOff, bool ionTorrent, bool reAlignTorrent, uint32_t hRunSize,
		uint64_t & highQualityBases, uint64_t & lowQualityBases, bool debug, uint32_t readPosCutOff,
		size_t start = 0, size_t end = std::numeric_limits<size_t>::max());

bool determinQualityIonTorrent(aligner & alignerObj, uint32_t hRunSize,
		const seqInfo & query, uint64_t pos, uint32_t qualCutOff, bool debug);


void extractRef(
		concurrent::BamReaderPool & bamPool, concurrent::AlignerPool & alnPool,
		njh::concurrent::LockableQueue<std::string>& chroms,
		MultiRefCounter& masterCounter,std::unordered_map<std::string, BamCountExtractStats>& refExtracts,
		std::mutex & refExtractsMut,
		 uint32_t mappingQual, size_t minLen, uint32_t qualCutOff,
		bool ionTorrent, bool reAlignTorrent, uint32_t hRunSize, bool verbose, uint32_t readPosCutOff);

}  // namespace njhseq



