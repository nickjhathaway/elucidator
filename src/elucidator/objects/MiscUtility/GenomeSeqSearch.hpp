/*
 * GenomeSeqSearch.hpp
 *
 *  Created on: Aug 26, 2017
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

#include <njhseq/objects/BioDataObject/GenomicRegion.hpp>
#include <njhseq/objects/kmer/kmerInfo.hpp>



namespace njhseq {

class GenomeSeqSearch {
public:
	class SeqStreak {
	public:
		SeqStreak(size_t seqPos, const GenomicRegion & region);

		size_t seqPos_;
		GenomicRegion region_;
		bool growing_ = true;

		bool isNextPos(const std::string & chrom, const size_t gPos,
				const bool isRevCompl, const size_t seqPos, const uint32_t kLen);

		void increaseEndPos(uint32_t size = 1);


	};

	GenomeSeqSearch(const OutOptions & outOpts);

	std::vector<SeqStreak> growingStreaks_;

	std::mutex mut_;
	std::shared_ptr<std::ofstream> out_;

	uint32_t countCapForGrowing_ = 10000;

	void growStreaks(const std::string & kmer, size_t seqPos, bool isRevComp,
			const std::unordered_map<std::string, std::unique_ptr<kmerInfo>> & infos);

	void pruneStreaks(uint32_t minSize);

	void purgeAllStreaks(uint32_t minSize);


};



}  // namespace njhseq
