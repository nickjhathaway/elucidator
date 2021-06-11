/*
 * statsExp.cpp
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

#include "statsExp.hpp"
#include "elucidator/objects/BioDataObject.h"
#include "elucidator/objects/counters/DNABaseCounter.hpp"
#include <random>
#include <boost/math/distributions/normal.hpp>

namespace njhseq {




statsExpRunner::statsExpRunner()
    : njh::progutils::ProgramRunner(
          {
					 addFunc("pnorm", pnorm, false),
					 addFunc("tajimatest", tajimatest, false),
					 addFunc("boost_normal_dist_example", boost_normal_dist_example, false),
					 //
           },
          "statsExp") {}



int statsExpRunner::tajimatest(
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
//  double Dmin = (2/n - 1/a1)/std::sqrt(e2);
//  double Dmax = ((n + 1)/(2 * n) - 1/a1)/std::sqrt(e2);
//  double tmp1 = 1 + Dmin * Dmax;
//  double tmp2 = Dmax - Dmin;
//  double a = -tmp1 * Dmax/tmp2;
//  double b = tmp1 * Dmin/tmp2;
//  double p = pbeta((D - Dmin)/tmp2, b, a);
//  if (p < 0.5){
//  	p =  2 * p;
//  }else{
//  	p =  2 * (1 - p);
//  }
//  double Pval_normal = 2 * pnorm(-std::abs(D));
//  double Pval_beta = p;
  std::cout << "D = " << D << std::endl;
//  std::cout << "Pval_normal = " << Pval_normal << std::endl;
//  std::cout << "Pval_beta = " << Pval_beta << std::endl;

	return 0;
}





int statsExpRunner::boost_normal_dist_example(
		const njh::progutils::CmdArgs &inputCommands) {
	seqSetUp setUp(inputCommands);

	setUp.finishSetUp(std::cout);
	double step = 1.; // in z
	double range = 4; // min and max z = -range to +range.
	int precision = 17; // traditional tables are only computed to much lower precision.

	// Construct a standard normal distribution s
	boost::math::normal s; // (default mean = zero, and standard deviation = unity)
	std::cout  << "Standard normal distribution, mean = "<< s.mean()
	    << ", standard deviation = " << s.standard_deviation() << std::endl;

	std::cout << "Probability distribution function values" << std::endl;
	std::cout << "  z " "      pdf " << std::endl;
	std::cout.precision(5);
	for (double z = -range; z < range + step; z += step)
	{
	  std::cout << std::left << std::setprecision(3) << std::setw(6) << z << " "
	    << std::setprecision(precision) << std::setw(12) << boost::math::pdf(s, z) << std::endl;
	}
	std::cout.precision(6); // default


	// For a standard normal distribution
	std::cout << "Standard normal mean = "<< s.mean()
	  << ", standard deviation = " << s.standard_deviation() << std::endl;
	std::cout << "Integral (area under the curve) from - infinity up to z " << std::endl;
	std::cout << "  z " "      cdf " << std::endl;
	for (double z = -range; z < range + step; z += step)
	{
	  std::cout << std::left << std::setprecision(3) << std::setw(6) << z << " "
	    << std::setprecision(precision) << std::setw(12) << cdf(s, z) << std::endl;
	}
	std::cout.precision(6); // default

	double z = 2.;
	std::cout << "Area for z = " << z << " is " << cdf(s, z) << std::endl; // to get the area for z.

  std::cout << "95% of area has a z below " << quantile(s, 0.95) << std::endl;
// 95% of area has a z below 1.64485

  std::cout << "95% of area has a z between " << quantile(s, 0.975)
    << " and " << -quantile(s, 0.975) << std::endl;
// 95% of area has a z between 1.95996 and -1.95996
  double alpha1 = cdf(s, -1) * 2; // 0.3173105078629142
  std::cout << std::setprecision(17) << "Significance level for z == 1 is " << alpha1 << std::endl;

  double alpha[] = {0.3173105078629142, // z for 1 standard deviation.
    0.20, 0.1, 0.05, 0.01, 0.001, 0.0001, 0.00001 };

  std::cout << "level of significance (alpha)" << std::setprecision(4) << std::endl;
  std::cout << "2-sided       1 -sided          z(alpha) " << std::endl;
  for (int i = 0; i < sizeof(alpha)/sizeof(alpha[0]); ++i)
  {
    std::cout << std::setw(15) << alpha[i] << std::setw(15) << alpha[i] /2 << std::setw(10) << quantile(complement(s,  alpha[i]/2)) << std::endl;
    // Use quantile(complement(s, alpha[i]/2)) to avoid potential loss of accuracy from quantile(s,  1 - alpha[i]/2)
  }
  std::cout << std::endl;

  std::cout.precision(3);
  std::cout << std::showpoint << "cdf(s, s.standard_deviation()) = "
    << cdf(s, s.standard_deviation()) << std::endl;  // from -infinity to 1 sd
  std::cout << "cdf(complement(s, s.standard_deviation())) = "
    << cdf(complement(s, s.standard_deviation())) << std::endl;
  std::cout << "Fraction 1 standard deviation within either side of mean is "
    << 1 -  cdf(complement(s, s.standard_deviation())) * 2 << std::endl;
  std::cout << "Fraction 2 standard deviations within either side of mean is "
    << 1 -  cdf(complement(s, 2 * s.standard_deviation())) * 2 << std::endl;
  std::cout << "Fraction 3 standard deviations within either side of mean is "
    << 1 -  cdf(complement(s, 3 * s.standard_deviation())) * 2 << std::endl;

	return 0;
}


int statsExpRunner::pnorm(
		const njh::progutils::CmdArgs & inputCommands) {
	double mean = 0;
	double sd = 1;
	double val = std::numeric_limits<double>::max();
	seqSetUp setUp(inputCommands);
	setUp.setOption(mean, "--mean", "mean");
	setUp.setOption(sd, "--sd", "standard deviation");
	setUp.setOption(val, "--val", "val", true);

	setUp.finishSetUp(std::cout);

	boost::math::normal_distribution<double> ndist(mean, sd);

	std::cout << "val: " << val << ", cdf:" << boost::math::cdf(ndist,val) << std::endl;

	return 0;
}



}  // namespace njhseq



