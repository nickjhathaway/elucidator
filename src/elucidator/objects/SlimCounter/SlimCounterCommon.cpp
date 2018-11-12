/*
 * SlimCounterCommon.cpp
 *
 *  Created on: Jun 16, 2016
 *      Author: nick
 */


#include "SlimCounterCommon.hpp"

namespace njhseq {
namespace SlimCounter {

std::array<uint32_t, 127> genAColIndex(bool countLowerCase){
	std::array<uint32_t, 127> ret;
	ret.fill(4);
	ret['A'] = 0;
	ret['C'] = 1;
	ret['G'] = 2;
	ret['T'] = 3;
	if(countLowerCase){
		ret['a'] = 0;
		ret['c'] = 1;
		ret['g'] = 2;
		ret['t'] = 3;
	}
	return ret;
}

void genColIndex(bool countLowerCase){
	colIndex = genAColIndex(countLowerCase);
}

}  // namespace SlimCounter


}  // namespace njhseq
