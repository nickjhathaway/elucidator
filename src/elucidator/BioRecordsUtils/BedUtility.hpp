#pragma once
/*
 * BedUtility.hpp
 *
 *  Created on: Jan 17, 2018
 *      Author: nick
 */

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



};

} /* namespace njhseq */

