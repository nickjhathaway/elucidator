#pragma once
/*
 * SlimCounterCommon.hpp
 *
 *  Created on: Jun 16, 2016
 *      Author: nick
 */




#include "elucidator/common.h"

namespace njhseq {

namespace SlimCounter {
static std::array<uint32_t, 127> colIndex { 4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U,
		4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U,
		4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U,
		4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 0U,
		4U, 1U, 4U, 4U, 4U, 2U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 3U,
		4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U,
		4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U, 4U,
		4U, 4U, 4U, 4U };


std::array<uint32_t, 127> genAColIndex(bool countLowerCase);

void genColIndex(bool countLowerCase = false);

} // namespace  SlimCounter

}  // namespace njhseq




