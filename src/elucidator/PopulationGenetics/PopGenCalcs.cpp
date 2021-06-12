/*
 * PopGenCalcs.cpp
 *
 *  Created on: Mar 15, 2018
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

#include "PopGenCalcs.hpp"

namespace njhseq {

PopGenCalculator::TajimaTestRes PopGenCalculator::calcTajimaTest(uint32_t nInputSeqs, uint32_t nSegragtingSites, double meanPairwiseDifferences){
	double n = nInputSeqs;
  if (n < 4) {
  	std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "Tajima test requires at least 4 sequences"<< "\n";
		throw std::runtime_error{ss.str()};
  }

  if (nSegragtingSites < 1) {
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
  double D = (meanPairwiseDifferences - nSegragtingSites/a1)/std::sqrt(e1 * nSegragtingSites + e2 * nSegragtingSites * (nSegragtingSites - 1));
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
  return TajimaTestRes(D, Pval_normal, Pval_beta);
}






PopGenCalculator::PopDifferentiationMeasures PopGenCalculator::getOverallPopDiff(std::unordered_map<std::string, std::vector<PopHapInfo> > hapsForPopulations){

	if(hapsForPopulations.size() < 2){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << " error, popSeqs should at least be size 2 not " << hapsForPopulations.size() << "\n";
		throw std::runtime_error{ss.str()};
	}
	PopDifferentiationMeasures ret;

	uint32_t popK = hapsForPopulations.size();

	std::unordered_map<std::string, uint32_t> subPopSizes;
	std::unordered_set<uint32_t> allHapSeqs;
	std::unordered_map<std::string, std::unordered_map<uint32_t, double>> freqsForPopForHapSeq;
	std::unordered_map<std::string, std::unordered_map<uint32_t, uint32_t>> countsForPopForHapSeq;


	for(auto & hapsForPopulation : hapsForPopulations){
		auto totalForPop = PopHapInfo::getTotalPopCount(hapsForPopulation.second);
		subPopSizes[hapsForPopulation.first] = totalForPop;
		PopHapInfo::setProb(hapsForPopulation.second, totalForPop);
		double sumOfSquares = 0;
		for(const auto & hap : hapsForPopulation.second){

			allHapSeqs.emplace(hap.popUid_);
			freqsForPopForHapSeq[hapsForPopulation.first][hap.popUid_] = hap.prob_;
			countsForPopForHapSeq[hapsForPopulation.first][hap.popUid_] = hap.count_;
			sumOfSquares += std::pow(hap.prob_, 2.0);
		}
		ret.hjsSample_[hapsForPopulation.first] = 1 - sumOfSquares;
	}
	uint32_t maxHapId = *std::max_element(allHapSeqs.begin(), allHapSeqs.end());
	//this only works if population identifier is set in a logical way e.g. ids with ranging from 0 to haplotype count
	std::vector<PopHapInfo> acrossPopulationsHaps;
	for(const auto pos : iter::range(maxHapId + 1)){
		acrossPopulationsHaps.emplace_back(PopHapInfo(pos, 0));
	}
	for(auto & hapsForPopulation : hapsForPopulations){
		for(const auto & hap : hapsForPopulation.second){
			acrossPopulationsHaps[hap.popUid_].count_ += hap.count_;
		}
	}
	PopHapInfo::setProb(acrossPopulationsHaps);

	std::unordered_map<uint32_t, double> parametricProbsPerHap;
	for (const auto &popUid : allHapSeqs) {
		for (auto &freqWithPop : freqsForPopForHapSeq) {
			parametricProbsPerHap[popUid] += freqWithPop.second[popUid]/popK;
		}
	}
	//informativeness for assignment
	{

		for(const auto & popUid : allHapSeqs){

			double firstTerm = -parametricProbsPerHap[popUid] * log(parametricProbsPerHap[popUid]);
			double secondTerm = 0;
			for( auto & freqWithPop : freqsForPopForHapSeq){
				if(freqWithPop.second[popUid] > 0){
					secondTerm += freqWithPop.second[popUid]/hapsForPopulations.size() * log(freqWithPop.second[popUid]);
				}
			}
			ret.informativenessForAssignPerHap_[popUid] = firstTerm + secondTerm;
			ret.informativenessForAssign_ += firstTerm + secondTerm;
		}


		for(const auto & pop : freqsForPopForHapSeq){
			for(const auto & hap : pop.second){
				if(hap.second > 0){
					ret.informativenessForAssignPerPopulation_[pop.first] += hap.second/hapsForPopulations.size() * log(hap.second/parametricProbsPerHap[hap.first]);
				}
			}
		}
	}

	double sumOfHjsSample = 0;
	for(const auto & subPop : ret.hjsSample_){
		sumOfHjsSample += subPop.second;
	}
	ret.hsSample_ = sumOfHjsSample/ret.hjsSample_.size();
	double jtSample = 0;
	for(const auto & hapSeq : allHapSeqs){
		double freqSum = 0;
		for( auto & subPop : freqsForPopForHapSeq){
			freqSum += subPop.second[hapSeq];
		}
		jtSample += std::pow(freqSum/hapsForPopulations.size(), 2.0);
	}
	ret.htSample_ = 1 - jtSample;

	double harmonicMean = 0;
	double sumOfInverses = 0;
	uint32_t totalHaps = 0;
	for(const auto & popSize : subPopSizes){
		totalHaps += popSize.second;
		sumOfInverses += 1.0/popSize.second;
	}
	harmonicMean = subPopSizes.size()/sumOfInverses;

	ret.hsEst_ = (harmonicMean/(harmonicMean - 1)) * ret.hsSample_;
	ret.htEst_ = ret.htSample_ + (ret.hsEst_)/(harmonicMean * hapsForPopulations.size());

	ret.gst_ = (ret.htSample_ - ret.hsSample_)/ret.htSample_;
	ret.jostD_ = ((ret.htSample_ - ret.hsSample_)/(1 - ret.hsSample_)) * (hapsForPopulations.size()/(hapsForPopulations.size() - 1));

	ret.gstEst_ =  (ret.htEst_ - ret.hsEst_)/ret.htEst_;
	ret.jostDEst_ = ((ret.htEst_ - ret.hsEst_)/(1 - ret.hsEst_)) * (hapsForPopulations.size()/(hapsForPopulations.size() - 1));


	double a = 0;
	for(const auto & hapSeq : allHapSeqs){
		double sumOfFreqs = 0;
		double sumOfSqaureFreqs = 0;
		for( auto & subPop : freqsForPopForHapSeq){
			sumOfFreqs += subPop.second[hapSeq];
			sumOfSqaureFreqs += std::pow(subPop.second[hapSeq], 2.0);
		}
		a += (std::pow(sumOfFreqs, 2.0) - sumOfSqaureFreqs)/(subPopSizes.size() - 1);
	}
	ret.chaoA_ = a;
	double b = 0;
	for(const auto & hapSeq : allHapSeqs){
		for(auto & subPop : countsForPopForHapSeq){
			if(subPop.second[hapSeq] > 0){
				b += (subPop.second[hapSeq] *(subPop.second[hapSeq] - 1) )/static_cast<double>(subPopSizes[subPop.first] * (subPopSizes[subPop.first] - 1));
			}
		}
	}
	ret.chaoB_ = b;
	ret.jostDChaoEst_ = 1 - (a/b);
	return ret;
}




}  // namespace njhseq

