#pragma once

/*
 * RoughIlluminaProfiler.hpp
 *
 *  Created on: May 2, 2019
 *      Author: nicholashathaway
 */

#include "elucidator/BamToolsUtils.h"

namespace njhseq {

class RoughIlluminaProfiler{
	/**@todo add counts for overlapping errors, this could help estimate PCR error
	 *
	 */
public:
	struct Counts{
		std::unordered_map<uint32_t, uint32_t> positionErrorCounts;
		std::unordered_map<uint32_t, uint32_t> positionTotalCounts;

		std::unordered_map<char, std::unordered_map<char, uint32_t>> baseChangeCounts;

		std::unordered_map<uint32_t, std::unordered_map<char, std::unordered_map<char, uint32_t>>> baseChangeCountsPerPosition;

		std::unordered_map<uint32_t, std::unordered_map<uint32_t, uint32_t>> qualCounts;
		std::unordered_map<uint32_t, std::unordered_map<uint32_t, uint32_t>> qualErrorsCounts;

		void addOtherCounts(const Counts & other);

		void increaseCounts(const AlignmentResults & res);

		void writeProfiles(const std::string & prefix, bool overWrite);
	};

	Counts r1_counts;
	Counts r2_counts;

	void addOther(const RoughIlluminaProfiler & other);

	void increaseCounts(
			const BamTools::BamAlignment & bAln,
			const BamTools::RefVector & refData,
			TwoBit::TwoBitFile & tReader) ;

};


}  // namespace njhseq



