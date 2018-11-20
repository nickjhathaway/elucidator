#pragma once
//

//  rabbitPopulation.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 11/18/13.
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
namespace njhseq {
class rabitPopulationImmortal {

 public:
  rabitPopulationImmortal(int numOfOffspringProdcued)
      : numberOfMatureRabitPairs(0),
        numberOfNewBornRabitsPairs(1),
        numberOfOffspringProduced(numOfOffspringProdcued) {}
  uint64_t numberOfMatureRabitPairs;
  uint64_t numberOfNewBornRabitsPairs;
  uint64_t numberOfOffspringProduced;
  void oneMonthPasses() {
    uint64_t newBorns = numberOfMatureRabitPairs * numberOfOffspringProduced;
    numberOfMatureRabitPairs += numberOfNewBornRabitsPairs;
    numberOfNewBornRabitsPairs = newBorns;
  }
  void multipleMonthsPass(int numberOfMonths, std::ostream& out) {
    out << "Month: " << 1 << ", #RabitPairs: " << getNumberOfRabitPairs()
        << ", #newBorns: " << numberOfNewBornRabitsPairs
        << ", #mature: " << numberOfMatureRabitPairs << std::endl;
    for (const auto& i : iter::range(1, numberOfMonths)) {
      oneMonthPasses();
      out << "Month: " << i + 1 << ", #RabitPairs: " << getNumberOfRabitPairs()
          << ", #newBorns: " << numberOfNewBornRabitsPairs
          << ", #mature: " << numberOfMatureRabitPairs << std::endl;
    }
    std::cout << std::endl;
  }
  uint64_t getNumberOfRabitPairs() {
    return numberOfMatureRabitPairs + numberOfNewBornRabitsPairs;
  }
};

class rabitPopulationMortal {

 public:
  rabitPopulationMortal(int numOfOffspringProdcued, int lifeSpan)
      : _numberOfOffspringProduced(numOfOffspringProdcued),
        _lifeSpan(lifeSpan) {
    population[0] = 1;
    for (auto i : iter::range(1, lifeSpan)) {
      population[i] = 0;
    }
  }
  std::unordered_map<uint, uint64_t> population;
  int _numberOfOffspringProduced;
  int _lifeSpan;
  void oneMonthPasses() {
    uint64_t newBorns = 0;
    // get newborns
    for (auto i : iter::range(1, _lifeSpan)) {
      newBorns += population[i] * _numberOfOffspringProduced;
    }
    // age population
    for (auto i : iter::range(_lifeSpan - 1, 0, -1)) {
      population[i] = population[i - 1];
    }
    // add newBorns
    population[0] = newBorns;
  }
  void multipleMonthsPass(int numberOfMonths, std::ostream& out) {

    out << "Month: " << 1 << ", #RabitPairs: " << getNumberOfRabitPairs()
        << ", #newBorns: " << population[0] << std::endl;
    for (const auto& i : iter::range(1, numberOfMonths)) {
      oneMonthPasses();
      out << "Month: " << i + 1 << ", #RabitPairs: " << getNumberOfRabitPairs()
          << ", #newBorns: " << population[0] << std::endl;
    }
    std::cout << std::endl;
  }
  uint64_t getNumberOfRabitPairs() {
    uint64_t sizeOfPopulation = 0;
    for (auto pop : population) {
      sizeOfPopulation += pop.second;
    }
    return sizeOfPopulation;
  }
};
}  // namesapce njh

