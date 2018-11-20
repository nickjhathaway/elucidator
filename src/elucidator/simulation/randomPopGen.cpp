/*
 * randomPopGen.cpp
 *
 *  Created on: Mar 23, 2014
 *      Author: nickhathaway
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
#include "randomPopGen.hpp"

namespace njhseq {

void randomPopGen::runPcr(std::map<std::string, uint32_t>& startingReads,
                          double basalErrorRate, uint32_t rounds) {
  for (uint32_t round = 0; round < rounds; ++round) {
    njh::scopedStopWatch roundTime("round_" + estd::to_string(round),true);
    std::cout << "starting round: " << round << std::endl;
    runOnePcr(startingReads, basalErrorRate);
  }
}
void randomPopGen::runOnePcr(std::map<std::string, uint32_t>& startingReads,
                             double basalErrorRate) {
  std::map<std::string, uint32_t> currentRound;
  for (const auto& reads : startingReads) {
    // now attempt to mutate for the current round of pcr, do twice the ammount
    // due to two strands being duplicated
  	for (uint32_t mut = 0; mut <reads.second * 2 ; ++mut){
    //for (const auto& mut : iter::range(reads.second * 2)) {
      ++currentRound[eProfile_.mutateSeqSameErrorRate(
            reads.first, gen_, eProfile_.alphabet_, basalErrorRate)];
    }
  }
  // now add the newly mutated reads to the current reads
  for (const auto& newReads : currentRound) {
    startingReads[newReads.first] += newReads.second;
  }
}

} /* namespace njh */
