#pragma once
/*
 * expSeqUtil.hpp
 *
 *  Created on: Dec 9, 2016
 *      Author: nick
 */

#include "elucidator/common.h"

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




