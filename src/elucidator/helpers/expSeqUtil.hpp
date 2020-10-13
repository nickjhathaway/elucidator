#pragma once
/*
 * expSeqUtil.hpp
 *
 *  Created on: Dec 9, 2016
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

#include "elucidator/common.h"

#include <njhseq/objects/kmer/kmer.hpp>


namespace njhseq {

class expSeqUtil {
public:

	static std::unordered_map<std::string, uint32_t> getFuzzyKmerCount(
			const std::string &seq, uint32_t kLength, uint32_t allowableMutations,
			bool checkComplement);

	static VecStr getFuzzySharedMotif(const VecStr &strs, uint32_t kLength,
			uint32_t allowableMutations, bool checkComplement);

	static std::map<std::string, kmer> adjustKmerCountsForMismatches(
			const std::map<kmer, int> &kmers, int allowableMismatches);
	static std::map<std::string, kmer> adjustKmerCountsForMismatches(
			const std::map<std::string, kmer> &kmers, int allowableMismatches);

	static int getCyclopeptideLengthFromSprectumLength(uint64_t length);
	static std::vector<std::vector<char>> getPossibleCyclopeptideFromSpretrum(
			const std::vector<int> &spectrum);

	static int getNumberOfPossibleLinearPeptides(uint64_t lengthOfProtein);

};


}  // namespace njhseq




