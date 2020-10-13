/*
 * expSeqUtil.cpp
 *
 *  Created on: Dec 9, 2016
 *      Author: nick
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



#include "expSeqUtil.hpp"
#include "elucidator/simulation/mutator.hpp"
#include <njhseq/helpers/seqUtil.hpp>

namespace njhseq {

std::unordered_map<std::string, uint32_t> expSeqUtil::getFuzzyKmerCount(const std::string &seq,
                                                      uint32_t kLength,
                                                      uint32_t allowableMutations,
                                                      bool checkComplement) {
	std::unordered_map<std::string, uint32_t> ans;
  std::unordered_map<std::string, VecStr> alreadyMutated;
  for (auto i : iter::range(seq.size() - kLength + 1)) {
    std::string currentKmer = seq.substr(i, kLength);
    std::string currentKmerComplement = seqUtil::reverseComplement(currentKmer, "DNA");
    if (alreadyMutated.find(currentKmer) == alreadyMutated.end()) {
      if (allowableMutations == 1) {
        alreadyMutated[currentKmer] =
            mutator::getSingleMutations(currentKmer, true);
        if (checkComplement) {
          addOtherVec(alreadyMutated[currentKmer],
                      mutator::getSingleMutations(currentKmerComplement, true));
        }
      } else if (allowableMutations == 2) {
        alreadyMutated[currentKmer] =
            mutator::getUpToDoubleMutations(currentKmer, true);
        if (checkComplement) {
          addOtherVec(
              alreadyMutated[currentKmer],
              mutator::getUpToDoubleMutations(currentKmerComplement, true));
        }
      } else if (allowableMutations == 3) {
        alreadyMutated[currentKmer] =
            mutator::getUpToTripleMutations(currentKmer, true);
        if (checkComplement) {
          addOtherVec(
              alreadyMutated[currentKmer],
              mutator::getUpToTripleMutations(currentKmerComplement, true));
        }
      } else if (allowableMutations == 0) {
        // no mutations allowed don't add any mutated strings
        alreadyMutated[currentKmer] = {};
      } else {
      	std::stringstream ss;
        ss << "Only 1,2, or 3 mutation(s) supported, can't do "
                  << allowableMutations << std::endl;
        throw std::runtime_error{ss.str()};
      }
    }
    ++ans[currentKmer];
    if (checkComplement) {
      ++ans[currentKmerComplement];
    }
    for (const auto &mutant : alreadyMutated[currentKmer]) {
      ++ans[mutant];
    }
  }
  return ans;
}






VecStr expSeqUtil::getFuzzySharedMotif(const VecStr &strs, uint32_t kLength,
		uint32_t allowableMutations,
                                    bool checkComplement) {
  std::unordered_map<uint32_t, std::unordered_map<std::string, uint32_t>> kmers;
  VecStr ans;
  uint32_t count = 0;
  for (const auto &str : strs) {
    kmers[count] =
        getFuzzyKmerCount(str, kLength, allowableMutations, checkComplement);
    ++count;
  }
  for (const auto &firstMers : kmers[0]) {
    bool eachContains = true;
    for (auto i : iter::range(1, (int)kmers.size())) {
      if (kmers[i].find(firstMers.first) == kmers[i].end()) {
        eachContains = false;
        break;
      }
    }
    if (eachContains) {
      ans.push_back(firstMers.first);
    }
  }
  return ans;
}



std::map<std::string, kmer> expSeqUtil::adjustKmerCountsForMismatches(
    const std::map<kmer, int> &kmers, int allowableMismatches) {
  std::map<std::string, kmer> ans;
  for (const auto &k : kmers) {
    ans[k.first.k_] = k.first;
  }
  for (auto firstIter = kmers.begin(); firstIter != kmers.end(); ++firstIter) {
    for (auto secondIter = firstIter; secondIter != kmers.end(); ++secondIter) {
      if (secondIter == firstIter) {
        continue;
      }
      if (seqUtil::checkTwoEqualSeqs(firstIter->first.k_, secondIter->first.k_,
                            allowableMismatches)) {
        ans[firstIter->first.k_].count_ += secondIter->first.count_;
        addOtherVec(ans[firstIter->first.k_].positions_,
                    secondIter->first.positions_);
        ans[secondIter->first.k_].count_ += firstIter->first.count_;
        addOtherVec(ans[secondIter->first.k_].positions_,
                    firstIter->first.positions_);
      }
    }
  }
  return ans;
}
std::map<std::string, kmer> expSeqUtil::adjustKmerCountsForMismatches(
    const std::map<std::string, kmer> &kmers, int allowableMismatches) {
  std::map<std::string, kmer> ans = kmers;
  for (auto firstIter = kmers.begin(); firstIter != kmers.end(); ++firstIter) {
    for (auto secondIter = firstIter; secondIter != kmers.end(); ++secondIter) {
      if (secondIter == firstIter) {
        continue;
      }
      if (secondIter->second.count_ == 0 && firstIter->second.count_ == 0) {
        continue;
      }
      if (seqUtil::checkTwoEqualSeqs(firstIter->second.k_, secondIter->second.k_,
                            allowableMismatches)) {
        ans[firstIter->second.k_].count_ += secondIter->second.count_;
        addOtherVec(ans[firstIter->second.k_].positions_,
                    secondIter->second.positions_);
        ans[secondIter->second.k_].count_ += firstIter->second.count_;
        addOtherVec(ans[secondIter->second.k_].positions_,
                    firstIter->second.positions_);
      }
    }
  }
  return ans;
}

int expSeqUtil::getCyclopeptideLengthFromSprectumLength(uint64_t length) {
  return std::ceil(std::sqrt(length));
}
std::vector<std::vector<char>> expSeqUtil::getPossibleCyclopeptideFromSpretrum(
    const std::vector<int> &spectrum) {
  std::vector<std::vector<char>> ans;
  int lengthOfPeptide = 0;
  uint64_t lengthOfSpectrum = 0;
  bool containsZero = false;
  if (spectrum[0] == 0) {
    containsZero = true;
    lengthOfSpectrum = spectrum.size() - 2;
  } else {
    lengthOfSpectrum = spectrum.size() - 1;
  }
  lengthOfPeptide = getCyclopeptideLengthFromSprectumLength(lengthOfSpectrum);
  if (containsZero) {
    for (auto i : iter::range(lengthOfPeptide + 1)) {
      auto &spec = spectrum[i];
      if (spec == 0) {
        continue;
      }
      ans.push_back(aminoAcidInfo::infos::weightIntToAminoAcid.at(spec));
    }
  } else {
    for (auto i : iter::range(lengthOfPeptide)) {
      auto &spec = spectrum[i];
      ans.push_back(aminoAcidInfo::infos::weightIntToAminoAcid.at(spec));
    }
  }
  return ans;
}

int expSeqUtil::getNumberOfPossibleLinearPeptides(uint64_t lengthOfProtein) {
  int num = 0;
  for (auto i : iter::range(lengthOfProtein + 1)) {
    num += lengthOfProtein - i;
  }
  // for zero pepitide ""
  ++num;
  return num;
}



}  // namespace njhseq


