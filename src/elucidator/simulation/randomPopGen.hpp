
#pragma once
/*
 * randomPopGen.hpp
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
#include "elucidator/common.h"
#include "elucidator/simulation/errorProfile.hpp"

namespace njhseq {

class randomPopGen {

 public:
  // Constructors
  randomPopGen(const simulation::mismatchProfile& eProfile,
               const randomGenerator& gen)
      : eProfile_(eProfile), gen_(gen) {}

  // members
  simulation::mismatchProfile eProfile_;
  randomGenerator gen_;

  // functions
  void runPcr(std::map<std::string, uint32_t>& startingReads,
              double basalErrorRate, uint32_t rounds);
  void runOnePcr(std::map<std::string, uint32_t>& startingReads,
                 double basalErrorRate);
};

} /* namespace njhseq */


