//
//  randomStrGen.cpp
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
#include "randomStrGen.hpp"
#include "elucidator/simulation/simulationCommon.hpp"
#include "elucidator/simulation/simulationCommon.hpp"
namespace njhseq {
namespace simulation {

//even likelihood for all chars generators
std::multimap<double, std::string, std::less<double>> getEvenLikelihood(
    const std::vector<char> &letters) {
  std::multimap<double, std::string, std::less<double>> likelihoods;
  for (const auto &let : letters) {
    likelihoods.emplace(1.0 / letters.size(), std::string(1, let));
  }
  return likelihoods;
}
std::multimap<double, std::string, std::less<double>> createLikelihood(
    const std::vector<char> &letters, const std::vector<uint32_t> & counts){
	if(counts.size() != letters.size()){
		std::stringstream ss;
		ss << "Error in createLikelihood(const std::vector<char> &letters,"
				" const std::vector<uint32_t> & counts)" << std::endl;
		ss << "Size of counts differs from size of letters" << std::endl;
		ss << "Size of counts: " << counts.size() << std::endl;
		ss << "Counts: "; printVector(counts, ", ");
		ss << "Size of letters: " << letters.size() << std::endl;
		ss << "Letters: "; printVector(letters, ", ");
		throw std::runtime_error{ss.str()};
	}
	double countsSum = vectorSum(counts);

  std::multimap<double, std::string, std::less<double>> likelihoods;
  for (const auto pos : iter::range(letters.size())) {
    likelihoods.emplace(counts[pos]/countsSum, std::string(1, letters[pos]));
  }
  return likelihoods;
}

std::string evenRandStr(uint32_t size, const std::vector<char> &letters,
                        randomGenerator &gen) {
  return randStrMap(size, getEvenLikelihood(letters), gen);
}
VecStr evenRandStrs(uint32_t size, const std::vector<char> &letters,
                    randomGenerator &gen, uint32_t strNum) {
  return randStrsMap(size, getEvenLikelihood(letters), gen, strNum);
}
VecStr evenRandStrsRandLen(uint32_t minLen, uint32_t maxLen,
                           const std::vector<char> &letters,
                           randomGenerator &gen, uint32_t strNum) {
  return randStrsRandLenMap(minLen, maxLen, getEvenLikelihood(letters), gen,
                         strNum);
}
/*
std::string randStr(std::vector<letterCounter> counts, randomGenerator &gen){
	return randStr<letterCounter, char>(counts, gen);
}*/

}  // sim
}  // njh
