/*
 * BedUtility.cpp
 *
 *  Created on: Jan 17, 2018
 *      Author: nick
 */

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
