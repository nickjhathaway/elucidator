/*
 * randStrGen.cpp
 *
 *  Created on: Jul 27, 2014
 *      Author: nickhathaway
 */
//
// njhseq - A library for analyzing sequence data
// Copyright (C) 2012-2016 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
// Jeffrey Bailey <Jeffrey.Bailey@umassmed.edu>
//
// This file is part of njhseq.
//
// njhseq is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// njhseq is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with njhseq.  If not, see <http://www.gnu.org/licenses/>.
//
#include "randStrGen.hpp"

namespace njhseq {
std::string randStrGen::rStr(uint64_t size){
  std::string ans;
  ans.reserve(size);
  while (ans.size() < size) {
  	ans.push_back(charGen_.genObj());
  }
  return ans;
}

VecStr randStrGen::rStrs(uint64_t size, uint32_t strNum){
  std::vector<std::string> ans(strNum);
  std::generate_n(ans.begin(), strNum,
                  [&]() { return rStr(size) ; });
  return ans;
}

std::string randStrGen::rStr(uint64_t minSize, uint64_t maxSize){
	//
	return rStr(rGen_.unifRand(minSize, maxSize + 1));
}

VecStr randStrGen::rStrs(uint64_t minSize, uint64_t maxSize, uint32_t num){
  std::vector<std::string> ans(num);
  std::generate_n(ans.begin(), num,
                  [&]() { return rStr(rGen_.unifRand(minSize, maxSize + 1)) ; });
  return ans;
}
} /* namespace njh */
