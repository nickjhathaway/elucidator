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
#include "simulationCommon.hpp"


namespace njhseq {

std::map<int32_t, int32_t, std::greater<int>> generate(randomGenerator& gen,
                                               uint32_t start, uint32_t stop,
                                               uint32_t num, bool verbose) {
  std::map<int32_t, int32_t, std::greater<int32_t>> counts;
  for (auto i : iter::range(start, stop)) {
    counts[i] = 0;
  }
  auto randNums = gen.unifRandVector(start, stop, num);
  for (auto i : randNums) {
    ++counts[i];
  }
  if (verbose) {
    for (const auto& c : counts) {
      std::cout << c.first << " : " << c.second << std::endl;
    }
    std::cout << std::endl;
  }
  return counts;
}

namespace simulation {

std::array<double, 100> makeQualErrorArr() {
  std::array<double, 100> arr;
  for (auto i : iter::range(100)) {
    arr[i] = std::pow(10.0, (-i / 10.0));
  }
  return arr;
}
//randomGenerator.hpp
const std::array<double, 100> constants::QualErrorArr = makeQualErrorArr();

}  // simulation
}  // njh
