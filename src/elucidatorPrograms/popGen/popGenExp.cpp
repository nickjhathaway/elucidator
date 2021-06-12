/*
 * popGenExp.cpp
 *
 *  Created on: Jun 11, 2021
 *      Author: nicholas hathaway
 */

// elucidator - A library for analyzing sequence data
// Copyright (C) 2012-2021 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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

#include "popGenExp.hpp"
#include "elucidator/objects/BioDataObject.h"
#include "elucidator/objects/counters/DNABaseCounter.hpp"
#include "elucidator/PopulationGenetics.h"

namespace njhseq {



popGenExpRunner::popGenExpRunner()
    : njh::progutils::ProgramRunner(
          {
					 addFunc("tajimatest", tajimatest, false),
					 addFunc("tajimatest_testingExample", tajimatest_testingExample, false),
					 addFunc("doPairwiseComparisonOnHapsSharing", doPairwiseComparisonOnHapsSharing, false),
					 addFunc("doPairwiseComparisonOnHapsSharingDev", doPairwiseComparisonOnHapsSharing, true),

					 addFunc("callVariantsAgainstRefSeq", callVariantsAgainstRefSeq, false),
					 addFunc("callVariantsAgainstRefSeqIndividual", callVariantsAgainstRefSeqIndividual, false),
					 addFunc("quickHaplotypeInformation", quickHaplotypeInformation, false),
					 addFunc("getHapPopDifAndVariantsInfo", getHapPopDifAndVariantsInfo, true),
					 addFunc("oldQuickHaplotypeInformationAndVariants", oldQuickHaplotypeInformationAndVariants, false),
					 addFunc("quickHaplotypeVariantsWithRegion", quickHaplotypeVariantsWithRegion, false),
					 //
           },
          "popGenExp") {}




int popGenExpRunner::tajimatest(
		const njh::progutils::CmdArgs & inputCommands) {
	double n = std::numeric_limits<uint32_t>::max();
	uint32_t S = std::numeric_limits<uint32_t>::max();
	double khat = std::numeric_limits<uint32_t>::max();
//
//  khat <- mean(dist.dna(x, "N"))
//  S <- length(seg.sites(x))
	seqSetUp setUp(inputCommands);
	setUp.setOption(n, "--numOfSeqs", "Number of input seqs", true);
	setUp.setOption(S, "--segSites", "Number of segregating sites", true);
	setUp.setOption(khat, "--meanPairDiff", "The mean number of pairwise differences", true);
	setUp.finishSetUp(std::cout);
  auto res = PopGenCalculator::calcTajimaTest(n,S, khat);
  std::cout << "D = " << res.d_ << std::endl;
  std::cout << "Pval_normal = " << res.pval_normal_ << std::endl;
  std::cout << "Pval_beta = " << res.pval_beta_ << std::endl;

	return 0;
}

int popGenExpRunner::tajimatest_testingExample(
		const njh::progutils::CmdArgs & inputCommands) {
	double n = 995;
	uint32_t S = 77;
	double khat = 25.64669;
//
//  khat <- mean(dist.dna(x, "N"))
//  S <- length(seg.sites(x))
	seqSetUp setUp(inputCommands);
	setUp.setOption(n, "--numOfSeqs", "Number of input seqs");
	setUp.setOption(S, "--segSites", "Number of segregating sites");
	setUp.setOption(khat, "--meanPairDiff", "The mean number of pairwise differences");
	setUp.finishSetUp(std::cout);


  if (n < 4) {
  	std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "Tajima test requires at least 4 sequences"<< "\n";
		throw std::runtime_error{ss.str()};
  }

  if (S < 1) {
  	std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "No segrating sites"<< "\n";
		throw std::runtime_error{ss.str()};
  }
  std::vector<uint32_t> tmp(n-1);njh::iota<uint32_t>(tmp, 1);
  double a1 = 0;
  double a2 = 0;
  for(const auto & t : tmp){
  	a1 += 1/static_cast<double>(t);
  	a2 += 1/std::pow(static_cast<double>(t), 2);
  }
  double b1 = (n + 1)/(3 * (n - 1));
  double b2 = 2 * (std::pow(n,2) + n + 3)/(9 * n * (n - 1));
  double c1 = b1 - 1/a1;
  double c2 = b2 - (n + 2)/(a1 * n) + a2/std::pow(a1,2);
  double e1 = c1/a1;
  double e2 = c2/(std::pow(a1,2) + a2);
  double D = (khat - S/a1)/std::sqrt(e1 * S + e2 * S * (S - 1));
  double Dmin = (2/n - 1/a1)/std::sqrt(e2);
  double Dmax = ((n + 1)/(2 * n) - 1/a1)/std::sqrt(e2);
  double tmp1 = 1 + Dmin * Dmax;
  double tmp2 = Dmax - Dmin;
  double a = -tmp1 * Dmax/tmp2;
  double b = tmp1 * Dmin/tmp2;
  boost::math::beta_distribution<double> betadist(b, a);
  double p = boost::math::cdf(betadist,(D - Dmin)/tmp2);
  if (p < 0.5){
  	p =  2 * p;
  }else{
  	p =  2 * (1 - p);
  }
  boost::math::normal ndist;
  double Pval_normal = 2 * boost::math::cdf(ndist, -std::abs(D));
  double Pval_beta = p;

  std::cout << "D = " << D << std::endl;
  std::cout << "Pval_normal = " << Pval_normal << std::endl;
  std::cout << "Pval_beta = " << Pval_beta << std::endl;

  auto res = PopGenCalculator::calcTajimaTest(n,S, khat);
  std::cout << "D = " << res.d_ << std::endl;
  std::cout << "Pval_normal = " << res.pval_normal_ << std::endl;
  std::cout << "Pval_beta = " << res.pval_beta_ << std::endl;

	return 0;
}





}  // namespace njhseq



