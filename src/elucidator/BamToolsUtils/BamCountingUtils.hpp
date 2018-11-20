#pragma once
/*
 * BamCountingUtils.hpp
 *
 *  Created on: Mar 2, 2016
 *      Author: nick
 */

// elucidator - A library for analyzing sequence data
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
//
// This file is part of elucidator.
//
// elucidator is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// elucidator is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with elucidator.  If not, see <http://www.gnu.org/licenses/>.
//

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



