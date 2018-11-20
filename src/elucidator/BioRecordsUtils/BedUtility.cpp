/*
 * BedUtility.cpp
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


#include "BedUtility.hpp"

namespace njhseq {

Bed3RecordCore & BedUtility::extendLeft(Bed3RecordCore & b, uint32_t extendLeftLen){
	b.chromStart_ = (b.chromStart_ > extendLeftLen ? b.chromStart_ - extendLeftLen : 0);
	return b;
}

Bed3RecordCore & BedUtility::extendRight(Bed3RecordCore & b, uint32_t extendRightLen, uint32_t chromLength){
	b.chromEnd_ = (b.chromEnd_ + extendRightLen > chromLength ? chromLength : b.chromEnd_ + extendRightLen);
	return b;
}

Bed3RecordCore & BedUtility::extendRight(Bed3RecordCore & b, uint32_t extendRightLen){
	extendRight(b,extendRightLen, std::numeric_limits<uint32_t>::max());
	return b;
}

Bed3RecordCore & BedUtility::extendLeftRight(Bed3RecordCore & b, uint32_t extendLeftLen, uint32_t extendRightLen, uint32_t chromLength){
	extendRight(extendLeft(b, extendLeftLen), extendRightLen, chromLength);
	return b;
}
Bed3RecordCore & BedUtility::extendLeftRight(Bed3RecordCore & b, uint32_t extendLeftLen, uint32_t extendRightLen){
	extendLeftRight(b, extendLeftLen, extendRightLen, std::numeric_limits<uint32_t>::max());
	return b;
}

} /* namespace njhseq */
