#pragma once
//
//  errorProfile.h
//  sequenceTools
//
//  Created by Nicholas Hathaway on 1/2/14.
//  Copyright (c) 2014 Nicholas Hathaway. All rights reserved.
//
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






#include "elucidator/common.h"
#include "elucidator/simulation/simulationCommon.hpp"


namespace njhseq {
namespace simulation {

class mismatchProfile {

 public:
  // constructor
  /*! \brief Default Construct DNA Alphabet
   *
   */
  mismatchProfile() : alphabet_({'A', 'C', 'G', 'T'}) { reset(); }
  /*! \brief Construct with Given Alphabet
   *
   */
  mismatchProfile(const std::vector<char>& alphabet) : alphabet_(alphabet) {
    reset();
  }
  // Members, only upper case letters allowed
  /**
   * @todo increase size of array to include lower case letters
   */
  std::array<std::array<uint32_t, 26>, 26> counts_;
  std::array<std::array<double, 26>, 26> fractions_;
  std::vector<char> alphabet_;
  std::unordered_map<char, std::multimap<double, char> > probsByLet_;

  // functions
  void reset();
  void setEqualProb();
  // set the fractions
  void setFractions();
  // increase the counts of mismatches
  void increaseCountAmount(const std::string& ref, const std::string& compare,
                           uint32_t amount, const char& ignore = '-');

  void increaseCountAmount(char firstBase, char secondBase, uint32_t amount);
  // mutators based on errors
  char mutate(char firstBase, randomGenerator& gen,
              const std::vector<char>& mutateTo) const;
  void mutateInPlace(char& firstBase, randomGenerator& gen,
                     const std::vector<char>& mutateTo) const;

  void mutateSeqInPlace(std::string& seq, randomGenerator& gen,
                        const std::vector<char>& mutateTo,
                        const std::vector<double>& likelihood);

  std::string mutateSeq(const std::string& seq, randomGenerator& gen,
                        const std::vector<char>& mutateTo,
                        const std::vector<double>& likelihood);

  void mutateSeqInPlaceSameErrorRate(std::string& seq, randomGenerator& gen,
                                     const std::vector<char>& mutateTo,
                                     double errorRate);
  std::string mutateSeqSameErrorRate(const std::string& seq,
                                     randomGenerator& gen,
                                     const std::vector<char>& mutateTo,
                                     double errorRate);
  // printing
  void quickPrintProfile(std::ostream& out);
  void quickPrintCounts(std::ostream& out);
  void quickPrintProbs(std::ostream& out);
};

}  // simulation
}  // njhseq


