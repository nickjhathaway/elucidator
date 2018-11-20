#pragma once
//
//  baseContainer.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 3/7/14.
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


namespace njhseq {
class variant {
public:
	variant(const seqInfo & seqBase);
	seqInfo seqBase_;
	//maps keys are position in reference
	std::unordered_map<uint32_t, mismatch> mismatches_;
	std::unordered_map<uint32_t, gap> insertions_;
	std::unordered_map<uint32_t, gap> deletions_;

	void outputInfo(std::ofstream & out, const seqInfo & ref)const;
	/*
	std::string cigarString()const{

		 * R'(\d+)(\w)'

	}


	 * M	BAM_CMATCH	0
I	BAM_CINS	1
D	BAM_CDEL	2
N	BAM_CREF_SKIP	3
S	BAM_CSOFT_CLIP	4
H	BAM_CHARD_CLIP	5
P	BAM_CPAD	6
=	BAM_CEQUAL	7
X	BAM_CDIFF	8

	std::vector<std::pair<uint32_t, uint32_t>> cigarRepresentation()const{

	}*/

};

class refVariants {
public:

	refVariants(const seqInfo & seqBase);
	seqInfo seqBase_;
	std::vector<variant> variants_;

	void addVariant(const seqInfo & var, aligner & alignerObj,
			bool weighHomopolymer);
	std::vector<uint32_t> getVariantSnpLoci() const;
	std::map<uint32_t, std::vector<char>> getVariantSnpLociMap() const;

	std::vector<uint32_t> getVariantSnpLoci(VecStr names,
			uint32_t expand = 0) const;

	std::map<uint32_t, std::vector<char>> getVariantSnpLociMap(VecStr names,
			uint32_t expand = 0) const;

	void outPut(std::ofstream & out) const;

};

}  // namespace njh


