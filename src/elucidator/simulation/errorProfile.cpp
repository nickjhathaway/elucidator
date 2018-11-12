//
//  errorProfile.cpp
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
#include "errorProfile.hpp"

namespace njhseq {
namespace simulation {

void mismatchProfile::increaseCountAmount(const std::string &ref,
                                       const std::string &compare,
                                       uint32_t amount, const char &ignore) {
  auto mis = std::make_pair(ref.begin(), compare.begin());
  while (mis.first != ref.end()) {
    mis = std::mismatch(mis.first, ref.end(), mis.second);
    if (mis.first != ref.end()) {
      if (*mis.first != ignore && *mis.second != ignore) {
        increaseCountAmount(*mis.first, *mis.second, amount);
      }
      ++mis.first;
      ++mis.second;
    }
  }
}

void mismatchProfile::increaseCountAmount(char firstBase, char secondBase,
                                       uint32_t amount) {
  counts_[firstBase - 'A'][secondBase - 'A'] += amount;
}

void mismatchProfile::setFractions() {
  // sum across letter counts to get full amount of changes and then set the
  // fractions
  for (auto i : iter::range(counts_.size())) {
    uint32_t currentSum = 0;
    for (auto j : iter::range(counts_[i].size())) {
      currentSum += counts_[i][j];
    }
    if (currentSum != 0) {
      for (auto j : iter::range(counts_[i].size())) {
        fractions_[i][j] = counts_[i][j] / static_cast<double>(currentSum);
      }
    }
  }
  // now calc the prob of change to each letter in set alphabet
  for (const auto &first : alphabet_) {
    std::multimap<double, char> currentProbsByChar;
    for (const auto &second : alphabet_) {
      currentProbsByChar.emplace(fractions_[first - 'A'][second - 'A'], second);
    }
    probsByLet_[first] = currentProbsByChar;
  }
}
void mismatchProfile::reset() {
  for (auto i : iter::range(fractions_.size())) {
    for (auto j : iter::range(fractions_[i].size())) {
      fractions_[i][j] = 0;
      counts_[i][j] = 0;
    }
  }
  probsByLet_.clear();
}

// mutate and return a char
char mismatchProfile::mutate(char firstBase, njh::randomGenerator &gen,
                          const std::vector<char> &mutateTo) const {
  mutateInPlace(firstBase, gen, mutateTo);
  return firstBase;
}

// mutate a char to another char
void mismatchProfile::mutateInPlace(char &firstBase, njh::randomGenerator &gen,
                                 const std::vector<char> &mutateTo) const {
  double sum = 0;
  double rando = gen.unifRand();
  auto search = probsByLet_.find(firstBase);
  if(search != probsByLet_.end()){
    for (const auto &prob : search->second) {
      sum += prob.first;
      if (sum > rando) {
        firstBase = prob.second;
        return;
      }
    }
  }
}
// mutate the given sequence
void mismatchProfile::mutateSeqInPlace(std::string &seq, njh::randomGenerator &gen,
                                    const std::vector<char> &mutateTo,
                                    const std::vector<double> &likelihood) {
  std::vector<double> rands = gen.unifRandVector(likelihood.size());
  for (const auto &i : iter::range(rands.size())) {
    if (rands[i] <= likelihood[i]) {
      mutateInPlace(seq[i], gen, mutateTo);
    }
  }
}
// mutate and return a new sequence
std::string mismatchProfile::mutateSeq(const std::string &seq,
                                    njh::randomGenerator &gen,
                                    const std::vector<char> &mutateTo,
                                    const std::vector<double> &likelihood) {
  std::string ans = seq;
  mutateSeqInPlace(ans, gen, mutateTo, likelihood);
  return ans;
}
void mismatchProfile::mutateSeqInPlaceSameErrorRate(
    std::string &seq, njh::randomGenerator &gen, const std::vector<char> &mutateTo,
    double errorRate) {
  std::vector<double> rands = gen.unifRandVector(seq.size());
  for (const auto &i : iter::range(rands.size())) {
    if (rands[i] <= errorRate) {
      mutateInPlace(seq[i], gen, mutateTo);
    }
  }
}
std::string mismatchProfile::mutateSeqSameErrorRate(
    const std::string &seq, njh::randomGenerator &gen,
    const std::vector<char> &mutateTo, double errorRate) {
  std::string ans = seq;
  mutateSeqInPlaceSameErrorRate(ans, gen, mutateTo, errorRate);
  return ans;
}
void mismatchProfile::quickPrintProfile(std::ostream &out) {
  std::vector<VecStr> outTable;
  VecStr header = concatVecs({"let"}, numVecToVecStr(alphabet_));
  for (auto first : alphabet_) {
    VecStr currentRow;
    currentRow.emplace_back(estd::to_string(first));
    for (auto second : alphabet_) {
      currentRow.emplace_back(estd::to_string(fractions_[first - 'A'][second - 'A']));
    }
    outTable.emplace_back(currentRow);
  }
  printTableOrganized(outTable, header, out);
}

void mismatchProfile::quickPrintCounts(std::ostream &out) {
  std::vector<VecStr> outTable;
  VecStr header = concatVecs({"let"}, numVecToVecStr(alphabet_));
  for (auto first : alphabet_) {
    VecStr currentRow;
    currentRow.emplace_back(estd::to_string(first));
    for (auto second : alphabet_) {
      currentRow.emplace_back(estd::to_string(counts_[first - 'A'][second - 'A']));
    }
    outTable.emplace_back(currentRow);
  }
  printTableOrganized(outTable, header, out);
}
void mismatchProfile::quickPrintProbs(std::ostream &out) {
  setFractions();

  std::vector<VecStr> outTable;
  VecStr header = {"let", "secondLet", "prob"};

  for (const auto &let : alphabet_) {
    for (const auto &probs : probsByLet_[let]) {
      outTable.emplace_back(VecStr{estd::to_string(let), estd::to_string(probs.second),
                                   estd::to_string(probs.first)});
      // out << let <<"\t" << probs.first << "\t" << probs.second << std::endl;
    }
  }
  printTableOrganized(outTable, header, out);
}

void mismatchProfile::setEqualProb() {
  reset();
  for (const auto &let : alphabet_) {
    for (const auto &letSecond : alphabet_) {
      if (let == letSecond) {
        continue;
      }
      counts_[let][letSecond] = 1;
    }
  }
  setFractions();
}
}  // sim
}  // njh
