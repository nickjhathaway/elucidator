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
#include <njhseq/objects/counters/DNABaseCounter.hpp>
#include <njhseq/PopulationGenetics/PopGenCalcs.hpp>
#include <njhseq/objects/dataContainers/tables/TableReader.hpp>

#include <boost/math/distributions.hpp>

namespace njhseq {




statsExpRunner::statsExpRunner()
    : njh::progutils::ProgramRunner(
          {
					 addFunc("stats_pnorm", stats_pnorm, false),
					 addFunc("stats_pbeta", stats_pbeta, false),
					 addFunc("boost_normal_dist_example", boost_normal_dist_example, false),
					 addFunc("fisher_exact", fisher_exact, false),
					 addFunc("fisher_exact_test", fisher_exact_test, false),
					 addFunc("fisher_exact_tableInput", fisher_exact_tableInput, false),
					 //
           },
          "statsExp") {}






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




int statsExpRunner::stats_pnorm(
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

int statsExpRunner::stats_pbeta(
		const njh::progutils::CmdArgs & inputCommands) {
	double shape1 = 0;
	double shape2 = 1;
	double val = std::numeric_limits<double>::max();
	seqSetUp setUp(inputCommands);
	setUp.setOption(shape1, "--shape1", "shape 1");
	setUp.setOption(shape2, "--shape2", "shape 2");
	setUp.setOption(val, "--val", "val", true);

	setUp.finishSetUp(std::cout);

	boost::math::beta_distribution<double> betadist(shape1, shape2);

	std::cout << "val: " << val << ", cdf:" << boost::math::cdf(betadist,val) << std::endl;

	return 0;
}



int statsExpRunner::fisher_exact_test(
				const njh::progutils::CmdArgs & inputCommands) {

	seqSetUp setUp(inputCommands);


	setUp.finishSetUp(std::cout);


	//int a = atoi(argv[1]), b = atoi(argv[2]), c = atoi(argv[3]), d = atoi(argv[4]);

	//int32_t a = 45, b =12, c = 18, d = 39;
	//     [,1] [,2]
	// [1,] TP   FN
	// [2,] FP   TN


	int32_t TP = 45, FP =12, FN = 18, TN = 39;
	int32_t k = TP + FN;
	int32_t m = TP + FP;
	int32_t total = TP + FP + FN + TN;
	int32_t n = FN + TN;
	int32_t lo = std::max(0, k - n);
	int32_t hi = std::min(k, m);
	//x <- x[1L, 1L]
	int32_t x = TP;
	std::vector<int32_t> support(hi + 1 - lo);
	njh::iota(support, lo);
	//m, k, total
	boost::math::hypergeometric_distribution<double> hg_dist(m, k, total);
	std::vector<double> logdc;
	logdc.reserve(support.size())	;
	for(const auto val : support){
		logdc.emplace_back(log(boost::math::pdf(hg_dist, val)));
		std::cout << "val: " << val << ", pdf: " << boost::math::pdf(hg_dist, val) << ", pdflog: " << log(boost::math::pdf(hg_dist, val)) << std::endl;
	}
	std::cout << "exp(1): " << exp(1) << std::endl;
	auto dnhyper = [&logdc,&support](double ncp){
		std::vector<double> ret(logdc.size());
		for(const auto idx : iter::range(logdc.size())){
			ret[idx] = logdc[idx] + log(ncp) * support[idx];
		}
		auto maxret = vectorMaximum(ret);
		for(const auto idx : iter::range(logdc.size())){
			ret[idx] = exp(ret[idx] - maxret);
		}
		auto sumret = vectorSum(ret);
		for(const auto idx : iter::range(logdc.size())){
			ret[idx] = ret[idx]/sumret;
		}
		return ret;
	};
	auto test_dnhyper = dnhyper(1);
	std::cout << "test_dnhyper: "  << std::endl;
	for(const auto idx : iter::range(test_dnhyper.size())){
		std::cout << test_dnhyper[idx] << std::endl;
	}

	auto mnhyper = [&lo,&hi,&dnhyper,&support](double ncp){
		if (ncp == 0){
			return static_cast<double>(lo);
		}
		if (ncp == std::numeric_limits<double>::infinity()){
			return static_cast<double>(hi);
		}
		double sum = 0;
		auto dnhyper_res = dnhyper(ncp);
		for(const auto idx : iter::range(support.size())){
			sum += dnhyper_res[idx] * support[idx];
		}
		return sum;
	};


	auto pnhyper = [&lo,&hi,&dnhyper,&support,&hg_dist,&x](double q, double ncp =1, bool upperTail = false)  {
		if (ncp == 1) {
			if(upperTail){
				//1 - phyper(..., lower.tail = TRUE) is the same as phyper(..., lower.tail = FALSE)
				return 1 - boost::math::cdf(hg_dist, x -1);
			} else {
				return boost::math::cdf(hg_dist, x);
			}
		}
		if (ncp == 0) {
			if(upperTail){
				return q <= lo ? 1. : 0.;
			}else{
				return q >= lo ? 1. : 0.;
			}
		}
		if (ncp == std::numeric_limits<double>::infinity()){
			if(upperTail){
				return q <= hi ? 1. : 0.;
			}else{
				return q >= hi ? 1. : 0.;
			}
		}
		double sum = 0;
		auto dnhyper_res = dnhyper(ncp);
		if(upperTail){
			for(const auto idx : iter::range(support.size())){
				if(support[idx] >= q){
					sum += dnhyper_res[idx];
				}
			}
		}else{
			for(const auto idx : iter::range(support.size())){
				if(support[idx] <= q){
					sum += dnhyper_res[idx];
				}
			}
		}
		return sum;
	};

	auto pnhyper_twotailed = [&lo,&hi,&dnhyper,&support,&x](const double test_odds)  {
		if (test_odds == 0){
			return x == lo ? 1. : 0.;
		} else if (test_odds == std::numeric_limits<double>::infinity()){
			return x == hi ? 1. : 0.;
		} else {
			double relErr = 1 + 1e-07;
			auto dnhyper_res = dnhyper(test_odds);
			double sum = 0;
			double testAgainst = dnhyper_res[x - lo] * relErr;
//			std::cout << "x - lo + 1: " << x - lo + 1 << std::endl;
//			std::cout << "x - lo: " << x - lo<< std::endl;
//			std::cout << "dnhyper_res[x - lo]: " << dnhyper_res[x - lo] << std::endl;
//			std::cout << "testAgainst: " << testAgainst  << std::endl;
//			std::cout << "relErr: " << relErr << std::endl;
			for(const auto idx : iter::range(support.size())){
//				std::cout << "dnhyper_res[" << idx << "]: " << dnhyper_res[idx] << ", " << njh::colorBool(dnhyper_res[idx] <= testAgainst)<< std::endl;
				if(dnhyper_res[idx] <= testAgainst){
					sum += dnhyper_res[idx];
				}
			}
			return sum;
		}
	};
	auto mle = [&lo,&hi,&mnhyper](const double x) {
		if (x == lo) {
			return 0.0;
		}
		if (x == hi){
			return std::numeric_limits<double>::infinity();
		}
	boost::uintmax_t max_iter=500; // Set max iterations.

	boost::math::tools::eps_tolerance<double> tol(30); //Set the eps tolerance.

//



		double mu = mnhyper(1);
		if (mu > x) {
			auto funcToRoot = [&mnhyper,&x](double input){
				return mnhyper(input) - x;
			};
			auto r1 = boost::math::tools::toms748_solve(funcToRoot, 0.0, 1.0, tol, max_iter); // use the toms solve algorithm.
//			std::cout << r1.first << std::endl;
//			std::cout << r1.second << std::endl;
			return r1.first;
			//uniroot(function(t) mnhyper(t) - x, c(0, 1))$root
		} else if (mu < x) {
			auto funcToRoot = [&mnhyper,&x](double input){
				return mnhyper(1/input) - x;
			};
			auto r1 = boost::math::tools::toms748_solve(funcToRoot, std::numeric_limits<double>::epsilon(), 1.0, tol, max_iter); // use the toms solve algorithm.
//			std::cout << r1.first << std::endl;
//			std::cout << r1.second << std::endl;
//			std::cout << 1/static_cast<double>(r1.first) << std::endl;
//			std::cout << 1/static_cast<double>(r1.second) << std::endl;
			return 1/r1.first;
			//1 / uniroot(function(t) mnhyper(1 / t) - x, c(.Machine$double.eps, 1))$root
		}
		return 1.0;
	};

	auto upperConf = [&pnhyper,&hi](const double x, const double alpha){
		if (x == hi) {
			return std::numeric_limits<double>::infinity();
		}
		boost::uintmax_t max_iter=500; // Set max iterations.

		boost::math::tools::eps_tolerance<double> tol(30); //Set the eps tolerance.

		double p = pnhyper(x, 1);
		if (p < alpha) {
			auto funcToRoot = [&pnhyper,&x,&alpha](double t){
				return pnhyper(x, t) - alpha;
			};
			auto r1 = boost::math::tools::toms748_solve(funcToRoot, 0.0, 1.0, tol, max_iter); // use the toms solve algorithm.
//			std::cout << r1.first << std::endl;
//			std::cout << r1.second << std::endl;
			return r1.first;
			//uniroot(function(t) pnhyper(x, t) - alpha, c(0, 1))$root
		} else if (p > alpha) {
			auto funcToRoot = [&pnhyper,&x,&alpha](double t){
				return pnhyper(x, 1/t) - alpha;
			};
			auto r1 = boost::math::tools::toms748_solve(funcToRoot, std::numeric_limits<double>::epsilon(), 1.0, tol, max_iter); // use the toms solve algorithm.
//			std::cout << r1.first << std::endl;
//			std::cout << r1.second << std::endl;
			return 1/r1.first;
			//1 / uniroot(function(t) pnhyper(x, 1 / t) - alpha, c(.Machine$double.eps, 1))$root
		} else{
			return 1.0;
		}
	};

	auto lowerConf = [&pnhyper,&lo](const double x, const double alpha){
		if (x == lo) {
			return std::numeric_limits<double>::infinity();
		}
		boost::uintmax_t max_iter=500; // Set max iterations.

		boost::math::tools::eps_tolerance<double> tol(30); //Set the eps tolerance.

		double p = pnhyper(x, 1, true);
//		std::cout << "p: " << p << std::endl;
		if (p > alpha) {
			auto funcToRoot = [&pnhyper,&x,&alpha](double t){
				return pnhyper(x, t, true) - alpha;
			};
			auto r1 = boost::math::tools::toms748_solve(funcToRoot, 0.0, 1.0, tol, max_iter); // use the toms solve algorithm.
//			std::cout << r1.first << std::endl;
//			std::cout << r1.second << std::endl;
			return r1.first;
			//uniroot(function(t) pnhyper(x, t) - alpha, c(0, 1))$root
		} else if (p < alpha) {
			auto funcToRoot = [&pnhyper,&x,&alpha](double t){
				return pnhyper(x, 1/t, true) - alpha;
			};
			auto r1 = boost::math::tools::toms748_solve(funcToRoot, std::numeric_limits<double>::epsilon(), 1.0, tol, max_iter); // use the toms solve algorithm.
//			std::cout << r1.first << std::endl;
//			std::cout << r1.second << std::endl;
			return 1/r1.first;
			//1 / uniroot(function(t) pnhyper(x, 1 / t) - alpha, c(.Machine$double.eps, 1))$root
		} else{
			return 1.0;
		}
	};
//	ncp.U <- function(x, alpha) {
//		if (x == hi)
//			return(Inf)
//		p <- pnhyper(x, 1)
//		if (p < alpha)
//			uniroot(function(t) pnhyper(x, t) - alpha,
//						c(0, 1))$root
//		else if (p > alpha)
//			1/uniroot(function(t) pnhyper(x, 1/t) - alpha,
//						c(.Machine$double.eps, 1))$root
//		else 1
//	}
//	ncp.L <- function(x, alpha) {
//		if (x == lo)
//			return(0)
//		p <- pnhyper(x, 1, upper.tail = TRUE)
//		if (p > alpha)
//			uniroot(function(t) pnhyper(x, t, upper.tail = TRUE) -
//													alpha, c(0, 1))$root
//		else if (p < alpha)
//			1/uniroot(function(t) pnhyper(x, 1/t, upper.tail = TRUE) -
//														alpha, c(.Machine$double.eps, 1))$root
//		else 1
//	}

//	mle <- function(x) {
//		if (x == lo)
//			return(0)
//		if (x == hi)
//			return(Inf)
//										mu <- mnhyper(1)
//		if (mu > x)
//			uniroot(function(t) mnhyper(t) - x, c(0, 1))$root
//		else if (mu < x)
//			1/uniroot(function(t) mnhyper(1/t) - x, c(.Machine$double.eps,
//						1))$root
//		else 1
//	}

	double confLevel = 0.95;
	std::cout << "test_mnhyper: "  << mnhyper(1) << std::endl;

	std::cout << "test_pnhyper less: "  << pnhyper(45, 1, false) << std::endl;
	std::cout << "test_pnhyper greater: "  << pnhyper(45, 1, true) << std::endl;

	std::cout << "pnhyper_twotailed: " << pnhyper_twotailed(1) << std::endl;
	std::cout << "std::numeric_limits<double>::epsilon(): " << std::numeric_limits<double>::epsilon() << std::endl;
	std::cout << "mle_test: " << mle(x) << std::endl;

	std::cout << "lowerConf(45, 1 - confLevel)" << lowerConf(45, 1 - confLevel) << std::endl;
	std::cout << "upperConf(45, 1 - confLevel): " <<  upperConf(45, 1 - confLevel) << std::endl;

	std::cout << "lowerConf(45, (1 - confLevel)/2): " << lowerConf(45, (1 - confLevel)/2) << std::endl;
	std::cout << "upperConf(45, (1 - confLevel)/2): " << upperConf(45, (1 - confLevel)/2) << std::endl;


//	std::cout << "boost::math::pdf(hg_dist, x): " << boost::math::pdf(hg_dist, x) << std::endl;
//	std::cout << "1 - boost::math::pdf(hg_dist, x): " << 1 - boost::math::pdf(hg_dist, x) << std::endl;
//
//	std::cout << "boost::math::pdf(hg_dist, x - 1): " << boost::math::pdf(hg_dist, x - 1) << std::endl;
//	std::cout << "1 - boost::math::pdf(hg_dist, x - 1): " << 1 - boost::math::pdf(hg_dist, x - 1) << std::endl;
//
//	std::cout << "boost::math::cdf(hg_dist, x): " << boost::math::cdf(hg_dist, x) << std::endl;
//	std::cout << "1 - boost::math::cdf(hg_dist, x): " << 1 - boost::math::cdf(hg_dist, x) << std::endl;
//
//	std::cout << "boost::math::cdf(hg_dist, x - 1): " << boost::math::cdf(hg_dist, x - 1) << std::endl;
//	std::cout << "1 - boost::math::cdf(hg_dist, x - 1): " << 1 - boost::math::cdf(hg_dist, x - 1) << std::endl;

//	dnhyper <- function(ncp) {
//		d <- logdc + log(ncp) * support
//		d <- exp(d - max(d))
//		d/sum(d)
//	}
//	mnhyper <- function(ncp) {
//		if (ncp == 0)
//			return(lo)
//		if (ncp == Inf)
//			return(hi)
//		sum(support * dnhyper(ncp))
//	}
//	pnhyper <- function(q, ncp = 1, upper.tail = FALSE) {
//		if (ncp == 1) {
//			return(if (upper.tail) phyper(x - 1, m, n, k,
//																		lower.tail = FALSE) else phyper(x, m, n, k))
//		}
//		if (ncp == 0) {
//			return(as.numeric(if (upper.tail) q <= lo else q >=
//																										 lo))
//		}
//		if (ncp == Inf) {
//			return(as.numeric(if (upper.tail) q <= hi else q >=
//																										 hi))
//		}
//		sum(dnhyper(ncp)[if (upper.tail) support >= q else support <=
//																											 q])
//	}
//	mle <- function(x) {
//		if (x == lo)
//			return(0)
//		if (x == hi)
//			return(Inf)
//										mu <- mnhyper(1)
//		if (mu > x)
//			uniroot(function(t) mnhyper(t) - x, c(0, 1))$root
//		else if (mu < x)
//			1/uniroot(function(t) mnhyper(1/t) - x, c(.Machine$double.eps,
//						1))$root
//		else 1
//	}


//	boost::uintmax_t max_iter=500; // Set max iterations.
//
//	boost::math::tools::eps_tolerance<double> tol(30); //Set the eps tolerance.
//
//	Result r1 = boost::math::tools::toms748_solve(myBinom95, 0, 1, tol, max_iter); // use the toms solve algorithm.

/*	dnhyper <- function(ncp) {
		d <- logdc + log(ncp) * support
		d <- exp(d - max(d))
		d/sum(d)
	}*/


/*
	int n = a + b + c + d;
	double* logFacs = new double[n+1]; // *** dynamically allocate memory logFacs[0..n] ***
	initLogFacs(logFacs, n); // *** initialize logFacs array ***
	double logpCutoff = logHypergeometricProb(logFacs,a,b,c,d); // *** logFacs added
	double pFraction = 0;
	for(int x=0; x <= n; ++x) {
		if ( a+b-x >= 0 && a+c-x >= 0 && d-a+x >=0 ) {
			double l = logHypergeometricProb(logFacs,x,a+b-x,a+c-x,d-a+x);
			if ( l <= logpCutoff ) pFraction += exp(l - logpCutoff);
		}
	}
	double logpValue = logpCutoff + log(pFraction);
	std::cout << "Two-sided log10-p-value is " << logpValue/log(10.) << std::endl;
	std::cout << "Two-sided p-value is " << exp(logpValue) << std::endl;
	delete [] logFacs;

*/


	return 0;


}




int statsExpRunner::fisher_exact(
				const njh::progutils::CmdArgs & inputCommands) {
	//     [,1] [,2]
	// [1,] TP   FN
	// [2,] FP   TN

	PopGenCalculator::FisherExactFor2x2::FisherExactFor2x2Input input;
	input.TP = 45;
	input.FP = 12;
	input.FN = 18;
	input.TN = 39;
	input.confInterval = 0.95;
	seqSetUp setUp(inputCommands);

	setUp.setOption(input.TP, "--truePositives,--TP", "Number of true positives");
	setUp.setOption(input.FP, "--falsePositives,--FP", "Number of false positives");
	setUp.setOption(input.FN, "--falseNegatives,--FN", "Number of false negatives");
	setUp.setOption(input.TN, "--trueNegatives,--TN", "Number of true negatives");
	setUp.setOption(input.confInterval, "--confInterval", "Confidence interval");

	setUp.finishSetUp(std::cout);

	auto res = PopGenCalculator::FisherExactFor2x2::runFisherExactOn2x2(input);
	std::cout << "odds ratio: " << res.oddsRatio_ << std::endl;
	std::cout << "lower " << roundDecPlaces(100*input.confInterval, 2)<< "% confidence: " << res.lowerConfInterval_ << std::endl;
	std::cout << "upper " << roundDecPlaces(100*input.confInterval, 2)<< "% confidence: " << res.upperConfInterval_ << std::endl;
	std::cout << "p-value: " << res.pValue_ << std::endl;

	return 0;


}

int statsExpRunner::fisher_exact_tableInput(
				const njh::progutils::CmdArgs & inputCommands) {
	//     [,1] [,2]
	// [1,] TP   FN
	// [2,] FP   TN
	double maxPvalue = 1;
	double minLowerConfInterval = 0;
	bool addOneToAll = false;
	VecStr IDColumnNames{ "ID"};
	bfs::path inputTable;
	OutOptions outOpts(bfs::path(""), ".tsv");
	seqSetUp setUp(inputCommands);
	setUp.setOption(addOneToAll, "--addOneToAll", "add One To All");
	setUp.setOption(maxPvalue, "--maxPvalue", "max P-value");
	setUp.setOption(minLowerConfInterval, "--minLowerConfInterval", "min Lower Conf Interval");

	setUp.setOption(IDColumnNames, "--IDColumnNames", "ID Column Names");
	setUp.setOption(inputTable, "--inputTable", "inputTable, needs 5 columns; IDColumnName,TP, FP, FN, TN", true);
	setUp.processWritingOptions(outOpts);

	setUp.finishSetUp(std::cout);
	TableReader tabReader(TableIOOpts::genTabFileIn(inputTable));

	tabReader.header_.checkForColumnsThrow(VecStr{"TP", "FP", "FN", "TN"}, __PRETTY_FUNCTION__ );
	tabReader.header_.checkForColumnsThrow(IDColumnNames, __PRETTY_FUNCTION__ );


	OutputStream out(outOpts);
	out << njh::conToStr(IDColumnNames, "\t") << "\t"<< "oddsRatio\tlowerConf\tupperConf\tp-value" << std::endl;
	VecStr row;
	while(tabReader.getNextRow(row)){
		PopGenCalculator::FisherExactFor2x2::FisherExactFor2x2Input fisherInput;
		fisherInput.TP = njh::StrToNumConverter::stoToNum<uint32_t>(row[tabReader.header_.getColPos("TP")]);
		fisherInput.FP = njh::StrToNumConverter::stoToNum<uint32_t>(row[tabReader.header_.getColPos("FP")]);
		fisherInput.FN = njh::StrToNumConverter::stoToNum<uint32_t>(row[tabReader.header_.getColPos("FN")]);
		fisherInput.TN = njh::StrToNumConverter::stoToNum<uint32_t>(row[tabReader.header_.getColPos("TN")]);
		if(addOneToAll){
			++fisherInput.TP;
			++fisherInput.FP;
			++fisherInput.FN;
			++fisherInput.TN;
		}
		auto res = PopGenCalculator::FisherExactFor2x2::runFisherExactOn2x2(fisherInput);
		if(res.pValue_ <= maxPvalue && res.lowerConfInterval_ >= minLowerConfInterval){
			for(const auto & colname : IDColumnNames){
				out << row[tabReader.header_.getColPos(colname)] << "\t";
			}
			out
							<< res.oddsRatio_
							<< "\t" << res.lowerConfInterval_
							<< "\t" << res.upperConfInterval_
							<< "\t" << res.pValue_ << std::endl;
		}
	}


	return 0;


}



}  // namespace njhseq



