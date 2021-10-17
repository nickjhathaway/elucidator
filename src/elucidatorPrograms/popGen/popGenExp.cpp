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
#include <njhseq/PopulationGenetics.h>
#include <njhseq/objects/seqContainers/CollapsedHaps.hpp>

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
					 addFunc("callVariantsAgainstRefGenome", callVariantsAgainstRefGenome, false),
					 addFunc("randomSamplingPloidyTest", randomSamplingPloidyTest, false),
					 addFunc("randomSamplingPloidyTest2", randomSamplingPloidyTest2, false),
					 //
           },
          "popGenExp") {}






int popGenExpRunner::randomSamplingPloidyTest(
		const njh::progutils::CmdArgs & inputCommands) {
	OutOptions outOpts("", ".tab.txt");
	bfs::path countsTableFnp;
	uint32_t ploidyTest = 2;
	double runs = 100;
	seqSetUp setUp(inputCommands);
	setUp.setOption(countsTableFnp, "--countsTableFnp", "counts Table Fnp, 2 columns, 1) name, 2) count", true);
	setUp.setOption(ploidyTest, "--ploidyTest", "ploidy to test");
	setUp.setOption(runs, "--runs", "sims to run");

	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	table countsTab(countsTableFnp, "\t", true);
	countsTab.changeHeaderToLowerCase();
	countsTab.checkForColumnsThrow(VecStr{"name", "count"}, __PRETTY_FUNCTION__);

	OutputStream out(outOpts);


	VecStr namesStr = countsTab.getColumn("name");
	VecStr countsStr = countsTab.getColumn("count");


	std::vector<double> counts;
	for (const auto &s : countsStr) {
		counts.emplace_back(njh::StrToNumConverter::stoToNum<double>(s));
	}

	std::vector<long double> freqs;
	long double total = vectorSum(counts);
	for(const auto count : counts){
		freqs.emplace_back(count/total);
	}



	std::vector<uint32_t> names(namesStr.size());
	njh::iota(names, 0U);

	njh::randObjectGen rObjGen(names, counts);
	std::unordered_map<uint32_t, uint32_t> ploidyCounts;
	for(uint32_t run = 0; run < runs; ++run){
		std::unordered_set<uint32_t> generated;
		for(uint32_t ploidy= 0; ploidy < ploidyTest; ++ploidy){
			generated.emplace(rObjGen.genObj());
		}
		++ploidyCounts[generated.size()];
	}
	out << "ploidyObserved\tcount\tfreqObs\tfreqExpected" << std::endl;
	std::unordered_map<uint32_t, std::string> expectedFreq;
	for(uint32_t ploidy= 1; ploidy < ploidyTest + 1; ++ploidy){
		expectedFreq[ploidy] = "NotImplemented";
	}
	//monoclonal
	long double sumOfMonoSquareFreqs = 0;
	for(const auto freq : freqs){
		sumOfMonoSquareFreqs += std::pow(freq, ploidyTest);
	}
	expectedFreq[1] = estd::to_string(sumOfMonoSquareFreqs);

	if (freqs.size() < ploidyTest){
		for(uint32_t coi = freqs.size() + 1 ; coi < ploidyTest; ++coi){
			expectedFreq[coi] = "0";
		}
	}
	if (2 == ploidyTest) {
		//for 2 clones it's just expected heterozygosity
		expectedFreq[2] = estd::to_string(1 - sumOfMonoSquareFreqs);
	} else if(3 == ploidyTest ){
		long double sumOfNotTroidy = 0;
		long double sumTwoClones = 0;
		for(const auto freq : freqs){
			long double monoClonal = std::pow(freq, 3);
			long double twoClones = 3 * std::pow(freq, 2) * (1 - freq);
			sumTwoClones += twoClones;
			sumOfNotTroidy += monoClonal + twoClones;
		}
		if (freqs.size() > 1){
			expectedFreq[2] = estd::to_string(sumTwoClones);
		}
		if (freqs.size() > 2) {
			expectedFreq[3] = estd::to_string(1 - sumOfNotTroidy);
		}
	} else if (ploidyTest == 4) {
		long double sumOf2Clones = 0;
		for(const auto freqPos : iter::range(freqs.size())){
			sumOf2Clones += std::pow(freqs[freqPos], 4 - 1) * (1 - freqs[freqPos]) * 4;
			for(const auto otherFreqPos : iter::range(freqs.size())){
				//plus all the times that you pick this hap twice and then pick the other hap twice as well, this occurs 3 times with a ploidy of 4
				if(otherFreqPos != freqPos){
					sumOf2Clones += std::pow(freqs[freqPos], 4 - 2) * std::pow(freqs[otherFreqPos], 4-2)  * (4 - 1);
				}
			}
		}
		long double sumOf3Clones = 0;
		for (const auto freqPos : iter::range(freqs.size())) {
			long double square = std::pow(freqs[freqPos], 2.0);
			for (const auto otherFreqPos : iter::range(freqs.size())) {
				//plus all the times that you pick this hap twice and then pick two haps that arent this hap or the hap you just choose, with a ploidy of 4 this happens 6 times
				if (otherFreqPos != freqPos) {
					long double allOtherFreqs = 1 - freqs[freqPos] - freqs[otherFreqPos];
					sumOf3Clones += square * freqs[otherFreqPos] * allOtherFreqs * 6;
				}
			}
		}
		if(freqs.size() > 1){
			expectedFreq[2] = estd::to_string(sumOf2Clones);
		}
		if(freqs.size() > 2){
			expectedFreq[3] = estd::to_string(sumOf3Clones);
		}
		if(freqs.size() > 3){
			expectedFreq[4] = estd::to_string(1 - sumOf3Clones - sumOf2Clones - sumOfMonoSquareFreqs);
		}
	} else if (5 == ploidyTest){
		long double sumOf2Clones = 0;
		for(const auto freqPos : iter::range(freqs.size())){

			for(const auto otherFreqPos : iter::range(freqs.size())){
				//plus all the times that you pick this hap twice and then pick the other hap twice as well, this occurs 3 times with a ploidy of 4
				if(otherFreqPos != freqPos){

//					sumOf2Clones += std::pow(freqs[freqPos], 4.0) * freqs[otherFreqPos] * 5;
//					sumOf2Clones += std::pow(freqs[freqPos], 3.0) * std::pow(freqs[otherFreqPos], 2.0) * 9;
					//sumOf2Clones += std::pow(freqs[freqPos], 2.0) * std::pow(freqs[otherFreqPos], 3.0) * 9;
					sumOf2Clones += std::pow(freqs[freqPos], 4.0) * freqs[otherFreqPos] * 5;
					sumOf2Clones += std::pow(freqs[freqPos], 3.0) * std::pow(freqs[otherFreqPos], 2.0) * 10;
				}
			}
		}
//		long double test = 0.50156;
//
//		long double sumOf3Clones = 0;
//		{
//
//			struct Ans {
//				Ans(uint32_t coef1, uint32_t coef2, double testVal, double testValDiff) :
//						coeff1_(coef1), coeff2_(coef2), testVal_(testVal), testValDiff_(
//								testValDiff) {
//
//				}
//				uint32_t coeff1_;
//				uint32_t coeff2_;
//				double testVal_;
//				double testValDiff_;
//			};
//			std::vector<Ans> answers;
//			for(uint32_t coeff1 = 5; coeff1 < 30; ++coeff1){
//				for(uint32_t coeff2 = 5; coeff2 < 30; ++coeff2){
//
//					long double currentTesting = 0;
//					//long double
//					for (const auto freqPos : iter::range(freqs.size())) {
//						for (const auto qFreqPos : iter::range(freqs.size())) {
//							if(qFreqPos == freqPos){
//								continue;
//							}
//							double zFreq = 1 - freqs[freqPos] - freqs[qFreqPos];
//							currentTesting += std::pow(freqs[freqPos], 3.0) * std::pow(freqs[qFreqPos], 1.0) * zFreq * coeff1;
//							currentTesting += std::pow(freqs[freqPos], 2.0) * std::pow(freqs[qFreqPos], 2.0) * zFreq * coeff2;
//							//sumOf3Clones += std::pow(freqs[freqPos], 3.0) * std::pow(freqs[qFreqPos], 1.0) * zFreq * 20;
//							//sumOf3Clones += std::pow(freqs[freqPos], 2.0) * std::pow(freqs[qFreqPos], 2.0) * zFreq * 30;
//							//sumOf3Clones += std::pow(freqs[freqPos], 1.0) * std::pow(freqs[qFreqPos], 3.0) * zFreq * 20;
//						}
//					}
//					answers.emplace_back(Ans(coeff1, coeff2, currentTesting, std::abs(currentTesting - test)));
//				}
//			}
//			njh::sort(answers, [](const Ans & a1, const Ans & a2){
//				return a1.testValDiff_ < a2.testValDiff_;
//			});
//			for(uint32_t top = 0; top < 10; ++top){
//				std::cout << "index: " << top  << std::endl;
//				std::cout << "answers[top].coeff1_     : " << answers[top].coeff1_ << std::endl;
//				std::cout << "answers[top].coeff2_     : " << answers[top].coeff2_ << std::endl;
//				std::cout << "answers[top].testVal_    : " << answers[top].testVal_ << std::endl;
//				std::cout << "answers[top].testValDiff_: " << answers[top].testValDiff_ << std::endl;
//
//				std::cout << std::endl;
//			}
//		}
		long double sumOf3Clones = 0;
		for (const auto freqPos : iter::range(freqs.size())) {
			for (const auto qFreqPos : iter::range(freqs.size())) {
				if(qFreqPos == freqPos){
					continue;
				}
				double zFreq = 1 - freqs[freqPos] - freqs[qFreqPos];
				sumOf3Clones += std::pow(freqs[freqPos], 3.0) * std::pow(freqs[qFreqPos], 1.0) * zFreq * 10;
				sumOf3Clones += std::pow(freqs[freqPos], 2.0) * std::pow(freqs[qFreqPos], 2.0) * zFreq * 15;
				//sumOf3Clones += std::pow(freqs[freqPos], 3.0) * std::pow(freqs[qFreqPos], 1.0) * zFreq * 20;
				//sumOf3Clones += std::pow(freqs[freqPos], 2.0) * std::pow(freqs[qFreqPos], 2.0) * zFreq * 30;
				//sumOf3Clones += std::pow(freqs[freqPos], 1.0) * std::pow(freqs[qFreqPos], 3.0) * zFreq * 20;
			}
		}
//				if (otherFreqPos != freqPos) {
//					long double allOtherFreqs = 1 - freqs[freqPos] - freqs[otherFreqPos];
//					sumOf3Clones += square * freqs[otherFreqPos] * allOtherFreqs * 6;
//				}
//			}
//		}
		//something isn't quite right here
	//	sumOf3Clones = sumOf3Clones * 2.0/3.0;
		//sumOf3Clones = sumOf3Clones * 2.0/3.0;
		//sumOf3Clones = sumOf3Clones * 0.67845;


		long double sumOf4Clones = 0;
		for (const auto freqPos : iter::range(freqs.size())) {
			for (const auto qFreqPos : iter::range(freqs.size())) {
				if(qFreqPos == freqPos){
					continue;
				}
				for (const auto rFreqPos : iter::range(freqs.size())) {
					if(rFreqPos == freqPos || rFreqPos == qFreqPos){
						continue;
					}
					double zFreq = 1 - freqs[freqPos] - freqs[qFreqPos] - freqs[rFreqPos];
					//sumOf4Clones += std::pow(freqs[freqPos], 2.0) * freqs[qFreqPos] * freqs[rFreqPos] * zFreq * 36;
					sumOf4Clones += std::pow(freqs[freqPos], 2.0) * freqs[qFreqPos] * freqs[rFreqPos] * zFreq * 10;
				}
			}
		}

		if(freqs.size() > 1){
			expectedFreq[2] = estd::to_string(sumOf2Clones);
		}
		if(freqs.size() > 2){
			expectedFreq[3] = estd::to_string(sumOf3Clones);
		}
		if(freqs.size() > 3){
			expectedFreq[4] = estd::to_string(sumOf4Clones);
		}
		if(freqs.size() > 4){
			expectedFreq[5] = estd::to_string(1 - sumOf4Clones - sumOf3Clones - sumOf2Clones - sumOfMonoSquareFreqs);
		}
	}


	for(uint32_t ploidy= 1; ploidy < ploidyTest + 1; ++ploidy){
		out << ploidy
				<< "\t" << ploidyCounts[ploidy]
				<< "\t" << ploidyCounts[ploidy]/static_cast<double>(runs)
				<< "\t" << expectedFreq[ploidy]
				<< std::endl;
	}


	return 0;
}



int popGenExpRunner::randomSamplingPloidyTest2(
		const njh::progutils::CmdArgs & inputCommands) {
	OutOptions outOpts("", ".tab.txt");
	bfs::path countsTableFnp;
	double runs = 100;
	seqSetUp setUp(inputCommands);
	setUp.setOption(countsTableFnp, "--countsTableFnp", "counts Table Fnp, 2 columns, 1) name, 2) count", true);
	setUp.setOption(runs, "--runs", "sims to run");

	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	table countsTab(countsTableFnp, "\t", true);
	countsTab.changeHeaderToLowerCase();
	countsTab.checkForColumnsThrow(VecStr{"name", "count"}, __PRETTY_FUNCTION__);

	OutputStream out(outOpts);
	out << "ploidyTest\tCOI\tobserved\tfreqObs\tfreqExpected\texpPoly" << std::endl;


	VecStr namesStr = countsTab.getColumn("name");
	VecStr countsStr = countsTab.getColumn("count");


	std::vector<double> counts;
	for (const auto &s : countsStr) {
		counts.emplace_back(njh::StrToNumConverter::stoToNum<double>(s));
	}

	std::vector<long double> freqs;
	long double total = vectorSum(counts);
	for(const auto count : counts){
		freqs.emplace_back(count/total);
	}

	std::vector<uint32_t> names(namesStr.size());
	njh::iota(names, 0U);

	njh::randObjectGen rObjGen(names, counts);
	for(uint32_t ploidyTest = 2; ploidyTest <=5; ++ploidyTest){
		std::unordered_map<uint32_t, uint32_t> ploidyCounts;
		for(uint32_t run = 0; run < runs; ++run){
			std::unordered_set<uint32_t> generated;
			for(uint32_t ploidy= 0; ploidy < ploidyTest; ++ploidy){
				generated.emplace(rObjGen.genObj());
			}
			++ploidyCounts[generated.size()];
		}
		auto infoForPloidy = PopGenCalculator::ExpectedPloidyInfo::genPloidyInfo(ploidyTest, freqs);
		for(uint32_t ploidy= 1; ploidy <=ploidyTest; ++ploidy){
			out << ploidyTest
					<< "\t" << ploidy
					<< "\t" << ploidyCounts[ploidy]
					<< "\t" << ploidyCounts[ploidy]/static_cast<double>(runs)
					<< "\t" << infoForPloidy.expectedCOIForPloidy_[ploidy]
				  << "\t" << infoForPloidy.expectedPolyClonal_
				  << std::endl;

		}
	}
	return 0;
}


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



