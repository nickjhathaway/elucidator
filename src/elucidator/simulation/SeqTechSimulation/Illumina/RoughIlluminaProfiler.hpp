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

		struct Indel {
			Indel(const gap & g, std::string refHomopolymer) :
					gapinfo_(g), refHomopolymer_(refHomopolymer) {
			}
			gap gapinfo_;
			std::string refHomopolymer_;
		};

		std::unordered_map<uint32_t, uint32_t> positionErrorCounts;
		std::unordered_map<uint32_t, uint32_t> positionTotalCounts;

		std::unordered_map<char, std::unordered_map<char, uint32_t>> baseChangeCounts;
		std::unordered_map<char, std::unordered_map<char, uint32_t>> allBaseCounts;

		std::unordered_map<uint32_t, std::unordered_map<char, std::unordered_map<char, uint32_t>>> baseChangeCountsPerPosition;

		std::unordered_map<uint32_t, std::unordered_map<uint32_t, uint32_t>> qualCounts;
		std::unordered_map<uint32_t, std::unordered_map<uint32_t, uint32_t>> qualErrorsCounts;

		std::unordered_map<uint32_t, std::vector<Indel>> deletions_;
		std::unordered_map<uint32_t, std::vector<Indel>> insertions_;;

		void addOtherCounts(const Counts & other);

		void increaseCounts(
				const seqInfo & refAln,
				const seqInfo & queryAln,
				const comparison & comp);

		void increaseCounts(const AlignmentResults & res);
		void increaseCounts(const ReAlignedSeq & res);

		void writeProfiles(const std::string & prefix, bool overWrite);
		void writeIndels(const std::string & prefix, bool overWrite);
		uint32_t softClipCutOff_{10};

		std::vector<double> percentIds_;
		uint32_t perfectHits_{0};
	};

	Counts r1_counts;
	Counts r2_counts;

	void addOther(const RoughIlluminaProfiler & other);

	void increaseCounts(
			const BamTools::BamAlignment & bAln,
			const BamTools::RefVector & refData,
			TwoBit::TwoBitFile & tReader) ;

	void increaseCounts(const ReAlignedSeq & res);


};


}  // namespace njhseq



