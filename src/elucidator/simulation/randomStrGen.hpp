#pragma once
//
//  randomStrGen.h
//  sequenceTools
//
//  Created by Nicholas Hathaway on 1/27/14.
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
#include "elucidator/simulation/simulationCommon.hpp"
namespace njhseq {


namespace simulation {
//even likelihood for all chars generators
std::multimap<double, std::string, std::less<double>> getEvenLikelihood(
    const std::vector<char> &letters);
std::multimap<double, std::string, std::less<double>> createLikelihood(
    const std::vector<char> &letters, const std::vector<uint32_t> & counts);
std::string evenRandStr(uint32_t size, const std::vector<char> &letters,
                        randomGenerator &gen);
VecStr evenRandStrs(uint32_t size, const std::vector<char> &letters,
                    randomGenerator &gen, uint32_t strNum) ;
VecStr evenRandStrsRandLen(uint32_t minLen, uint32_t maxLen,
                           const std::vector<char> &letters,
                           randomGenerator &gen, uint32_t strNum) ;

template<typename TYPE>
std::string randStrMap(
    uint32_t size,
    const std::multimap<double, TYPE, std::less<double>> &likelihoods,
    randomGenerator &gen) {
  std::string ans;
  ans.reserve(size);
  //std::vector<double> randos = gen.unifRandVector(size);
  while(ans.size() < size){
  	double sum = 0;
  	double rando = gen.unifRand();
		for (const auto &likelihood : likelihoods) {
			sum += likelihood.first;
			if (sum > rando) {
				addToStr(ans, likelihood.second);
				//ans.append(likelihood.second);
				break;
			}
		}
  }
  /*
  for (uint32_t i = 0; i < size; ++i) {
    double sum = 0;
    for (const auto &likelihood : likelihoods) {
      sum += likelihood.first;
      if (sum > randos[i]) {
      	addToStr(ans, likelihood.second);
        //ans.append(likelihood.second);
        break;
      }
    }
  }*/
  return ans;
}
template<typename TYPE>
std::string randStrMap(
    const std::vector<std::multimap<double, TYPE, std::less<double>>> &
        allLikelihoods,
    randomGenerator &gen) {
  std::string ans = "";
  ans.reserve(allLikelihoods.size());
  std::vector<double> randos = gen.unifRandVector(allLikelihoods.size());
  for (const auto &i : iter::range(allLikelihoods.size())) {
    double sum = 0;
    for (const auto &likelihood : allLikelihoods[i]) {
      sum += likelihood.first;
      if (sum > randos[i]) {
      	addToStr(ans, likelihood.second);
        //ans.append(likelihood.second);
        break;
      }
    }
  }
  return ans;
}

template<typename COUNTER>
std::string randStr(uint32_t size, COUNTER count, randomGenerator &gen) {
  auto likelihoods = count.createLikelihoodMaps(true);
  return randStrMap(size, likelihoods, gen);
}


template<typename COUNTER, typename TYPE>
std::string randStr(std::vector<COUNTER> counts, randomGenerator &gen) {
  std::vector<std::multimap<double, TYPE, std::less<double>>>
      allLikelihoods;
  for (auto &count : counts) {
    allLikelihoods.emplace_back(count.createLikelihoodMaps(true));
  }
  return randStrMap(allLikelihoods, gen);
}
//std::string randStr(std::vector<letterCounter> counts, randomGenerator &gen);


template<typename TYPE>
VecStr randStrsMap(
    uint32_t size,
    const std::multimap<double, TYPE, std::less<double>> &likelihoods,
    randomGenerator &gen, uint32_t strNum) {
  VecStr ans(strNum);
  std::generate_n(ans.begin(), strNum,
                  [&]() { return randStrMap(size, likelihoods, gen); });
  return ans;
}

template<typename COUNTER>
VecStr randStrs(uint32_t size, COUNTER count, randomGenerator &gen,
                uint32_t strNum) {
  auto likelihoods = count.createLikelihoodMaps(true);
  return randStrsMap(size, likelihoods, gen, strNum);
}

template<typename COUNTER, typename TYPE>
VecStr randStrs(std::vector<COUNTER> counts, randomGenerator &gen,
                uint32_t strNum) {
  std::vector<std::multimap<double, TYPE, std::less<double>>>
      allLikelihoods;
  for (auto &count : counts) {
    allLikelihoods.emplace_back(count.createLikelihoodMaps(true));
  }
  return randStrsMap(allLikelihoods, gen, strNum);
}

template<typename TYPE>
VecStr randStrsMap(
    const std::vector<std::multimap<double, TYPE, std::less<double>>> &
        allLikelihoods,
    randomGenerator &gen, uint32_t strNum) {
  VecStr ans(strNum);
  std::generate_n(ans.begin(), strNum,
                  [&]() { return randStr(allLikelihoods, gen); });
  return ans;
}

template<typename TYPE>
VecStr randStrsRandLenMap(
    uint32_t minLen, uint32_t maxLen,
    const std::multimap<double, TYPE, std::less<double>> &likelihoods,
    randomGenerator &gen, uint32_t strNum) {
  VecStr ans;
  ans.reserve(strNum);
  std::vector<uint32_t> randomLengths =
      gen.unifRandVector(minLen, maxLen + 1, strNum);
  for (const auto &len : randomLengths) {
    ans.emplace_back(randStrMap(len, likelihoods, gen));
  }
  return ans;
}

// multiple string generation at different lengths
template<typename COUNTER>
VecStr randStrsRandLen(uint32_t minLen, uint32_t maxLen, COUNTER count,
                       randomGenerator &gen, uint32_t strNum) {
  auto likelihoods = count.createLikelihoodMaps(true);
  return randStrsRandLenMap(minLen, maxLen, likelihoods, gen, strNum);
}



// from beginning
template<typename COUNTER, typename TYPE>
VecStr randStrsRandLen(uint32_t minLen, std::vector<COUNTER> counts,
                       randomGenerator &gen, uint32_t strNum) {
  std::vector<std::multimap<double, TYPE, std::less<double>>>
      allLikelihoods;
  for (auto &count : counts) {
    allLikelihoods.emplace_back(count.createLikelihoodMaps(true));
  }
  return randStrsRandLenMap(minLen, allLikelihoods, gen, strNum);
}
template<typename TYPE>
VecStr randStrsRandLenMap(
    uint32_t minLen,
    const std::vector<std::multimap<double, TYPE, std::less<double>>> &
        allLikelihoods,
    randomGenerator &gen, uint32_t strNum) {
  VecStr ans;
  ans.reserve(strNum);
  std::vector<uint32_t> randomLengths =
      gen.unifRandVector<uint32_t>(minLen, len(allLikelihoods) + 1, strNum);
  std::string randomString = "";
  ans.reserve(allLikelihoods.size());
  for (const auto &len : randomLengths) {
    randomString.clear();
    std::vector<double> randos = gen.unifRandVector(len);
    for (const auto &i : iter::range<uint32_t>(0, len)) {
      double sum = 0;
      for (const auto &likelihood : allLikelihoods[i]) {
        sum += likelihood.first;
        if (sum > randos[i]) {
        	addToStr(randomString, likelihood.second);
          //randomString.append(likelihood.second);
          break;
        }
      }
    }
    ans.emplace_back(randomString);
  }
  return ans;
}

// random start positions
template<typename COUNTER, typename TYPE>
VecStr randStrsRandLenRandPos(uint32_t minLen,
                              std::vector<COUNTER> counts,
                              randomGenerator &gen, uint32_t strNum) {
  std::vector<std::multimap<double, TYPE, std::less<double>>>
      allLikelihoods;
  for (auto &count : counts) {
    allLikelihoods.emplace_back(count.createLikelihoodMaps(true));
  }
  return randStrsRandLenRandPosMap(minLen, allLikelihoods, gen, strNum);
}
template<typename TYPE>
VecStr randStrsRandLenRandPosMap(
    uint32_t minLen,
    const std::vector<std::multimap<double, TYPE, std::less<double>>> &
        allLikelihoods,
    randomGenerator &gen, uint32_t strNum) {
  VecStr ans;
  ans.reserve(strNum);
  std::vector<uint32_t> randomLengths =
      gen.unifRandVector<uint32_t>(minLen, len(allLikelihoods) + 1, strNum);
  std::string randomString = "";
  ans.reserve(allLikelihoods.size());
  for (const auto &leng : randomLengths) {
    uint32_t start = gen.unifRand<uint32_t>(0, len(allLikelihoods) - leng);
    randomString.clear();
    std::vector<double> randos = gen.unifRandVector(leng);
    uint32_t pos = 0;
    for (const auto &i : iter::range<uint32_t>(start, leng + start)) {
      double sum = 0;
      for (const auto &likelihood : allLikelihoods[i]) {
        sum += likelihood.first;
        if (sum > randos[pos]) {
        	addToStr(randomString, likelihood.second);
        	//randomString.append(likelihood.second);
          break;
        }
      }
      ++pos;
    }
    ans.emplace_back(randomString);
  }
  return ans;
}

}  // simulation
}  // njhseq


