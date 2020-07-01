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

PopGenCalculator::PopDifferentiationMeasures PopGenCalculator::getOverallPopDiff(std::unordered_map<std::string, std::vector<PopHapInfo> > hapsForPopulations){

	if(hapsForPopulations.size() < 2){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << " error, popSeqs should at least be size 2 not " << hapsForPopulations.size() << "\n";
		throw std::runtime_error{ss.str()};
	}
	PopDifferentiationMeasures ret;


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
	std::vector<PopHapInfo> totalPopulationHaps;
	for(const auto pos : iter::range(maxHapId + 1)){
		totalPopulationHaps.emplace_back(PopHapInfo(pos, 0));
	}
	for(auto & hapsForPopulation : hapsForPopulations){
		for(const auto & hap : hapsForPopulation.second){
			totalPopulationHaps[hap.popUid_].count_ += hap.count_;
		}
	}
	PopHapInfo::setProb(totalPopulationHaps);


	//informativeness for assignment
	{

		for(const auto & popUid : allHapSeqs){

			double firstTerm = -totalPopulationHaps[popUid].prob_ * log(totalPopulationHaps[popUid].prob_);
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
					ret.informativenessForAssignPerPopulation_[pop.first] += hap.second/hapsForPopulations.size() * log(hap.second/totalPopulationHaps[hap.first].prob_);
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

