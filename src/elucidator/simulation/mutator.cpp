////  mutator.cpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 1/1/14.
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
#include "mutator.hpp"

#include <njhseq/utils.h>


namespace njhseq {

VecStr mutator::getSingleMutations(const std::string &originalSeq, bool sort) {
  VecStr ans;
  for (uint64_t i : iter::range(originalSeq.size())) {
    addOtherVec(ans, mutateAtPosition(originalSeq, i));
  }
  if (sort) std::sort(ans.begin(), ans.end());
  return ans;
}
VecStr mutator::getDoubleMutations(const std::string &originalSeq, bool sort) {
  VecStr ans;
  for (uint64_t i : iter::range(originalSeq.size())) {
    for (uint64_t j : iter::range(originalSeq.size())) {
      if (i == j) {
        continue;
      }
      addOtherVec(ans, mutateAtTwoPositions(originalSeq, i, j));
    }
  }
  removeDuplicates(ans);
  if (sort) std::sort(ans.begin(), ans.end());
  return ans;
}
VecStr mutator::getUpToDoubleMutations(const std::string &originalSeq,
                                       bool sort) {
  VecStr ans = getSingleMutations(originalSeq, false);
  addOtherVec(ans, getDoubleMutations(originalSeq, false));
  if (sort) std::sort(ans.begin(), ans.end());
  return ans;
}

VecStr mutator::getTripleMutations(const std::string &originalSeq, bool sort) {
  VecStr ans;
  for (uint64_t i : iter::range(originalSeq.size())) {
    for (uint64_t j : iter::range(originalSeq.size())) {
      if (i == j) {
        continue;
      }
      for (uint64_t k : iter::range(originalSeq.size())) {
        if (k == i || k == j) {
          continue;
        }
        addOtherVec(ans, mutateAtThreePositions(originalSeq, i, j, k));
      }
    }
  }
  removeDuplicates(ans);
  if (sort) std::sort(ans.begin(), ans.end());
  return ans;
  // return ans;
}
VecStr mutator::getUpToTripleMutations(const std::string &originalSeq,
                                       bool sort) {
  VecStr ans = getUpToDoubleMutations(originalSeq, false);
  addOtherVec(ans, getTripleMutations(originalSeq, false));
  if (sort) std::sort(ans.begin(), ans.end());
  return ans;
}

const VecStr mutator::mutateAtPosition(const std::string &originalSeq,
                                       uint64_t pos) {
  VecStr ans;
  std::string copy = originalSeq;
  if (originalSeq[pos] != 'C') {
    copy[pos] = 'C';
    ans.push_back(copy);
  }
  if (originalSeq[pos] != 'A') {
    copy[pos] = 'A';
    ans.push_back(copy);
  }
  if (originalSeq[pos] != 'T') {
    copy[pos] = 'T';
    ans.push_back(copy);
  }
  if (originalSeq[pos] != 'G') {
    copy[pos] = 'G';
    ans.push_back(copy);
  }
  return ans;
}
const VecStr mutator::mutateAtTwoPositions(const std::string &originalSeq,
                                           uint64_t pos1, uint64_t pos2) {
  VecStr ans;
  std::string copy = originalSeq;
  if (originalSeq[pos1] != 'C') {
    copy[pos1] = 'C';
    addOtherVec(ans, mutateAtPosition(copy, pos2));
  }
  if (originalSeq[pos1] != 'A') {
    copy[pos1] = 'A';
    addOtherVec(ans, mutateAtPosition(copy, pos2));
  }
  if (originalSeq[pos1] != 'T') {
    copy[pos1] = 'T';
    addOtherVec(ans, mutateAtPosition(copy, pos2));
  }
  if (originalSeq[pos1] != 'G') {
    copy[pos1] = 'G';
    addOtherVec(ans, mutateAtPosition(copy, pos2));
  }
  return ans;
}
const VecStr mutator::mutateAtThreePositions(const std::string &originalSeq,
                                             uint64_t pos1, uint64_t pos2,
                                             uint64_t pos3) {
  VecStr ans;
  std::string copy = originalSeq;
  if (originalSeq[pos1] != 'C') {
    copy[pos1] = 'C';
    addOtherVec(ans, mutateAtTwoPositions(copy, pos2, pos3));
  }
  if (originalSeq[pos1] != 'A') {
    copy[pos1] = 'A';
    addOtherVec(ans, mutateAtTwoPositions(copy, pos2, pos3));
  }
  if (originalSeq[pos1] != 'T') {
    copy[pos1] = 'T';
    addOtherVec(ans, mutateAtTwoPositions(copy, pos2, pos3));
  }
  if (originalSeq[pos1] != 'G') {
    copy[pos1] = 'G';
    addOtherVec(ans, mutateAtTwoPositions(copy, pos2, pos3));
  }
  return ans;
}
std::string mutator::mutateString(std::string seq,
                                  const std::vector<uint32_t> &qual,
                                  const simulation::mismatchProfile &profile,
                                  randomGenerator &gen,
                                  const std::vector<char> &mutateTo,
                                  uint32_t &mutateCount,
                                  const std::array<double, 100> &errorLookUp) {
  for (auto i : iter::range(seq.size())) {
    double rando = gen.unifRand();
    if (rando <= errorLookUp[qual[i]]) {
      ++mutateCount;
      profile.mutateInPlace(seq[i], gen, mutateTo);
    }
  }
  return seq;
}
}  // namespace njh
