#pragma once
/*
 * BedUtility.hpp
 *
 *  Created on: Jan 17, 2018
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


#include <njhseq/objects/BioDataObject/BedRecordCore.hpp>

namespace njhseq {


class BedUtility {
public:

	static Bed3RecordCore & extendLeft(Bed3RecordCore & b, uint32_t extendLeft);
	static Bed3RecordCore & extendRight(Bed3RecordCore & b, uint32_t extendRight, uint32_t chromLength);
	static Bed3RecordCore & extendRight(Bed3RecordCore & b, uint32_t extendRight);
	static Bed3RecordCore & extendLeftRight(Bed3RecordCore & b, uint32_t extendLeft, uint32_t extendRight, uint32_t chromLength);
	static Bed3RecordCore & extendLeftRight(Bed3RecordCore & b, uint32_t extendLeft, uint32_t extendRight);

	template<typename T>
	static void coordSort(std::vector<T> & beds, bool decending = false){
		auto bedCoordSorterFunc =
				[](const T & reg1In, const T & reg2In) {
					const auto & reg1 = getRef(reg1In);
					const auto & reg2 = getRef(reg2In);
					if(reg1.chrom_ == reg2.chrom_) {
						if(reg1.chromStart_ == reg2.chromStart_) {
							return reg1.chromEnd_ < reg2.chromEnd_;
						} else {
							return reg1.chromStart_ < reg2.chromStart_;
						}
					} else {
						return reg1.chrom_ < reg2.chrom_;
					}
				};
		if(decending){
			std::sort(beds.rbegin(), beds.rend(), bedCoordSorterFunc);
		}else{
			njh::sort(beds, bedCoordSorterFunc);
		}
	}

	template<typename BED>
	static uint32_t getPlotIDForBed(const BED & record,
			std::unordered_map<std::string,
					std::unordered_map<uint32_t, std::vector<uint32_t>>> & alreadyTakenIds) {
		uint32_t id = 0;
		std::set<uint32_t> alreadyTaken;
		for (const auto pos : iter::range(getRef(record).chromStart_, getRef(record).chromEnd_)) {
			for (const auto & otherId : alreadyTakenIds[getRef(record).chrom_][pos]) {
				alreadyTaken.emplace(otherId);
			}
		}
		while (njh::in(id, alreadyTaken)) {
			++id;
		}
		for (const auto pos : iter::range(getRef(record).chromStart_, getRef(record).chromEnd_)) {
			alreadyTakenIds[getRef(record).chrom_][pos].emplace_back(id);
		}
		return id;
	}





};

} /* namespace njhseq */

