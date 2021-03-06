#pragma once
/*
 * PopGenCalcs.hpp
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
#include "elucidator/common.h"

namespace njhseq {

struct DiversityMeasures{

	uint32_t alleleNumber_ = 0; //number of unique alleles
	uint32_t doublets_ = 0; //number of haplotypes found twice
	uint32_t singlets_ = 0; //number of haplotypes found only once
	double expShannonEntropy_ = std::numeric_limits<double>::max(); //exp of shannon entropy base e
	double ShannonEntropyE_ = std::numeric_limits<double>::max(); //shannon entropy base e
	double effectiveNumOfAlleles_ = std::numeric_limits<double>::max();
	double heterozygostiy_  = std::numeric_limits<double>::max();


};

/**@brief Get several general measures of diversity, assumes haps are already collapsed to unique haplotypes and have frequencies set
 *
 * @param haps a vector of unique haplotypes
 * @return a struc wtih several diversity measurements
 */
template<typename T>
DiversityMeasures getGeneralMeasuresOfDiversity(const std::vector<T> & haps){
	DiversityMeasures res;

	res.alleleNumber_ = haps.size();
	double sumOfSquares = 0;
	double sumOfLogFreqTimesFreq = 0;

	for (const auto & hap : haps) {
		const seqInfo & seqRef = getSeqBase(hap);
		if (1 == seqRef.cnt_) {
			++res.singlets_;
		} else if (2 == seqRef.cnt_) {
			++res.doublets_;
		}
		sumOfSquares += std::pow(seqRef.frac_, 2.0);
		sumOfLogFreqTimesFreq += seqRef.frac_ * std::log(seqRef.frac_);
	}
	res.heterozygostiy_ = 1 - sumOfSquares;
	res.effectiveNumOfAlleles_ = std::pow(sumOfSquares, -1);
	res.ShannonEntropyE_ = -sumOfLogFreqTimesFreq;
	res.expShannonEntropy_ = std::exp(-sumOfLogFreqTimesFreq);

	return res;
}

struct PopDifferentiationMeasures{

	std::unordered_map<std::string, double> hjsSample_;
	double hsSample_ = std::numeric_limits<double>::max();
	double htSample_ = std::numeric_limits<double>::max();

	double hsEst_ = std::numeric_limits<double>::max();
	double htEst_ = std::numeric_limits<double>::max();

	double gst_ = std::numeric_limits<double>::max();
	double jostD_ = std::numeric_limits<double>::max();

	double gstEst_ = std::numeric_limits<double>::max();
	double jostDEst_ = std::numeric_limits<double>::max();
	double chaoA_ = std::numeric_limits<double>::max();
	double chaoB_ = std::numeric_limits<double>::max();
	double jostDChaoEst_ = std::numeric_limits<double>::max();

};



template<typename T>
PopDifferentiationMeasures getOverallPopDiff(const std::unordered_map<std::string, std::shared_ptr<std::vector<T>>> & popSeqs){
	if(popSeqs.size() < 2){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << " error, popSeqs should at least be size 2 not " << popSeqs.size() << "\n";
		throw std::runtime_error{ss.str()};
	}
	PopDifferentiationMeasures ret;

	std::unordered_map<std::string, uint32_t> subPopSizes;
	std::unordered_set<std::string> allHapSeqs;
	std::unordered_map<std::string, std::unordered_map<std::string, double>> freqsForPopForHapSeq;
	std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>> countsForPopForHapSeq;

	for(const auto & subPop : popSeqs){
		subPopSizes[subPop.first] =
				std::accumulate(subPop.second->begin(), subPop.second->end(), 0,
						[](uint32_t total, const T & seq) {
							return total + std::round(getSeqBase(seq).cnt_);
						});
		double sumOfSquares = 0;
		for(const auto & seq : *subPop.second){
			allHapSeqs.emplace(getSeqBase(seq).seq_);
			freqsForPopForHapSeq[subPop.first][getSeqBase(seq).seq_] = getSeqBase(seq).frac_;
			countsForPopForHapSeq[subPop.first][getSeqBase(seq).seq_] = std::round(getSeqBase(seq).cnt_);
			sumOfSquares += std::pow(getSeqBase(seq).frac_, 2.0);
		}
		ret.hjsSample_[subPop.first] = 1 - sumOfSquares;
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
		jtSample += std::pow(freqSum/popSeqs.size(), 2.0);
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
	ret.htEst_ = ret.htSample_ + (ret.hsEst_)/(harmonicMean * popSeqs.size());

	ret.gst_ = (ret.htSample_ - ret.hsSample_)/ret.htSample_;
	ret.jostD_ = ((ret.htSample_ - ret.hsSample_)/(1 - ret.hsSample_)) * (popSeqs.size()/(popSeqs.size() - 1));

	ret.gstEst_ =  (ret.htEst_ - ret.hsEst_)/ret.htEst_;
	ret.jostDEst_ = ((ret.htEst_ - ret.hsEst_)/(1 - ret.hsEst_)) * (popSeqs.size()/(popSeqs.size() - 1));


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

template<typename T>
std::unordered_map<std::string,
		std::unordered_map<std::string, PopDifferentiationMeasures>> getPairwisePopDiff(
		const std::unordered_map<std::string, std::shared_ptr<std::vector<T>>> & popSeqs) {
	if(popSeqs.size() < 2){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << " error, popSeqs should at least be size 2 not " << popSeqs.size() << "\n";
		throw std::runtime_error{ss.str()};
	}
	std::unordered_map<std::string,
			std::unordered_map<std::string, PopDifferentiationMeasures>> ret;

	auto keys = getVectorOfMapKeys(popSeqs);
	njh::sort(keys);

	for(const auto & keyPos : iter::range(keys.size())){
		for(const auto & secondKeyPos : iter::range(keyPos)){
			std::unordered_map<std::string, std::shared_ptr<std::vector<T>>> currentPair;
			currentPair[keys[keyPos]] = popSeqs.at(keys[keyPos]);
			currentPair[keys[secondKeyPos]] = popSeqs.at(keys[secondKeyPos]);
			auto popMeasures = getOverallPopDiff(currentPair);
			ret[keys[keyPos]][keys[secondKeyPos]] = popMeasures;
		}
	}
	return ret;
}


}  // namespace njhseq





