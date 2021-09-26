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



PopGenCalculator::ExpectedPloidyInfo PopGenCalculator::ExpectedPloidyInfo::genPloidyInfo(uint32_t ploidy, const std::vector<long double> & freqs){

	PopGenCalculator::ExpectedPloidyInfo ret;
	ret.ploidy_ = ploidy;
	if(0 == ploidy || 1 == ploidy){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "ploidy can't be " << ploidy << "\n";
		throw std::runtime_error{ss.str()};
	}
	if(ploidy > 5){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "currently calculations for ploidy greater than 5 not implemented"<< "\n";
		throw std::runtime_error{ss.str()};
	}
	//monoclonal
	long double sumOfMonoSquareFreqs = 0;
	for(const auto freq : freqs){
		sumOfMonoSquareFreqs += std::pow(freq, ploidy);
	}
	ret.expectedPolyClonal_ = 1 - sumOfMonoSquareFreqs;
	ret.expectedCOIForPloidy_[1] = sumOfMonoSquareFreqs;

	if (freqs.size() < ploidy){
		for(uint32_t coi = freqs.size() + 1 ; coi < ploidy + 1; ++coi){
			ret.expectedCOIForPloidy_[coi] = 0;
		}
	}

	if (2 == ploidy) {
		//for 2 clones it's just expected heterozygosity
		ret.expectedCOIForPloidy_[2] = 1 - sumOfMonoSquareFreqs;
	} else if (3 == ploidy) {
		long double sumOfNotTroidy = 0;
		long double sumTwoClones = 0;
		for(const auto freq : freqs){
			long double monoClonal = std::pow(freq, 3);
			long double twoClones = 3 * std::pow(freq, 2) * (1 - freq);
			sumTwoClones += twoClones;
			sumOfNotTroidy += monoClonal + twoClones;
		}
		if (freqs.size() > 1){
			ret.expectedCOIForPloidy_[2] = sumTwoClones;
		}
		if (freqs.size() > 2) {
			ret.expectedCOIForPloidy_[3] = 1 - sumOfNotTroidy;
		}
	} else if (4 == ploidy) {

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
			ret.expectedCOIForPloidy_[2] = sumOf2Clones;
		}
		if(freqs.size() > 2){
			ret.expectedCOIForPloidy_[3] = sumOf3Clones;
		}
		if(freqs.size() > 3){
			ret.expectedCOIForPloidy_[4] = 1 - sumOf3Clones - sumOf2Clones - sumOfMonoSquareFreqs;
		}
	} else if (5 == ploidy) {
		long double sumOf2Clones = 0;
		for(const auto freqPos : iter::range(freqs.size())){
			for(const auto otherFreqPos : iter::range(freqs.size())){
				//plus all the times that you pick this hap twice and then pick the other hap twice as well, this occurs 3 times with a ploidy of 4
				if(otherFreqPos != freqPos){
					sumOf2Clones += std::pow(freqs[freqPos], 4.0) * freqs[otherFreqPos] * 5;
					sumOf2Clones += std::pow(freqs[freqPos], 3.0) * std::pow(freqs[otherFreqPos], 2.0) * 10;
				}
			}
		}
		long double sumOf3Clones = 0;
		for (const auto freqPos : iter::range(freqs.size())) {
			for (const auto qFreqPos : iter::range(freqs.size())) {
				if(qFreqPos == freqPos){
					continue;
				}
				double zFreq = 1 - freqs[freqPos] - freqs[qFreqPos];
				sumOf3Clones += std::pow(freqs[freqPos], 3.0) * std::pow(freqs[qFreqPos], 1.0) * zFreq * 10;
				sumOf3Clones += std::pow(freqs[freqPos], 2.0) * std::pow(freqs[qFreqPos], 2.0) * zFreq * 15;
			}
		}
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
					sumOf4Clones += std::pow(freqs[freqPos], 2.0) * freqs[qFreqPos] * freqs[rFreqPos] * zFreq * 10;
				}
			}
		}
		if(freqs.size() > 1){
			ret.expectedCOIForPloidy_[2] = sumOf2Clones;
		}
		if(freqs.size() > 2){
			ret.expectedCOIForPloidy_[3] = sumOf3Clones;
		}
		if(freqs.size() > 3){
			ret.expectedCOIForPloidy_[4] = sumOf4Clones;
		}
		if(freqs.size() > 4){
			ret.expectedCOIForPloidy_[5] = 1 - sumOf4Clones - sumOf3Clones - sumOf2Clones - sumOfMonoSquareFreqs;
		}
	}

	return ret;
}



PopGenCalculator::DiversityMeasures PopGenCalculator::getGeneralMeasuresOfDiversity(const std::vector<PopHapInfo> & haps){
	DiversityMeasures res;

	res.alleleNumber_ = haps.size();
	double sumOfSquares = 0;
	double sumOfLogFreqTimesFreq = 0;
	double totalHaps = PopHapInfo::getTotalPopCount(haps);
	double sumTopOfSimpson = 0;
	std::vector<long double> freqs;
	for (const auto & hap : haps) {
		freqs.emplace_back(hap.count_/totalHaps);
		if (1 == hap.count_) {
			++res.singlets_;
		} else if (2 == hap.count_) {
			++res.doublets_;
		}
		sumTopOfSimpson += hap.count_ * (hap.count_ - 1);
		double prob = hap.count_/totalHaps;
		sumOfSquares += std::pow(prob, 2.0);
		sumOfLogFreqTimesFreq += prob * std::log(prob);
	}
	res.heterozygostiy_ = 1 - sumOfSquares;
	res.effectiveNumOfAlleles_ = std::pow(sumOfSquares, -1);
	res.ShannonEntropyE_ = -sumOfLogFreqTimesFreq;
	res.expShannonEntropy_ = std::exp(-sumOfLogFreqTimesFreq);
	if(totalHaps > 1){
		res.simpsonIndex_ = 1 - sumTopOfSimpson/(totalHaps * (totalHaps - 1));
	}else{
		res.simpsonIndex_ = 0;
	}


	//ploidy 2
	res.ploidy2_ = ExpectedPloidyInfo::genPloidyInfo(2, freqs);
	//ploidy of 2 is kind of unnecssary since it's just expected hetereozygosity

	//ploidy 3
	res.ploidy3_ = ExpectedPloidyInfo::genPloidyInfo(3, freqs);

	//ploidy 4
	res.ploidy4_ = ExpectedPloidyInfo::genPloidyInfo(4, freqs);

	//ploidy 5
	res.ploidy5_ = ExpectedPloidyInfo::genPloidyInfo(5, freqs);


	return res;
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


PopGenCalculator::PopDifferentiationMeasuresPairWise PopGenCalculator::getPopDiff(
		const std::string & pop1, const std::vector<PopHapInfo> & pop1Haps,
		const std::string & pop2, const std::vector<PopHapInfo> & pop2Haps,
		const std::unordered_set<uint32_t> & allPossibleHaps,
		const std::unordered_map<uint32_t, std::unordered_map<uint32_t, double>> & pairwiseDistances){
	PopDifferentiationMeasuresPairWise ret(getOverallPopDiff(std::unordered_map<std::string, std::vector<PopHapInfo>>{{pop1, pop1Haps}, {pop2, pop2Haps}}),
			pop1, pop2);


	std::unordered_map<uint32_t, uint32_t> pop1HapsCounts;
	std::unordered_map<uint32_t, uint32_t> pop2HapsCounts;

	std::unordered_map<uint32_t, double> pop1HapsFreqs;
	std::unordered_map<uint32_t, double> pop2HapsFreqs;

	std::unordered_set<uint32_t> allHaps;

	uint32_t pop1Total = 0;
	uint32_t pop2Total = 0;


	for(const auto & pop1Hap : pop1Haps){
		allHaps.emplace(pop1Hap.popUid_);
		pop1HapsCounts[pop1Hap.popUid_] = pop1Hap.count_;
		pop1Total += pop1Hap.count_;
	}
	for(auto & pop1Hap : pop1HapsCounts){
		pop1HapsFreqs[pop1Hap.first] = pop1Hap.second/static_cast<double>(pop1Total);
	}

	for(const auto & pop2Hap : pop2Haps){
		allHaps.emplace(pop2Hap.popUid_);
		pop2HapsCounts[pop2Hap.popUid_] = pop2Hap.count_;
		pop2Total += pop2Hap.count_;
	}
	for(auto & pop2Hap : pop2HapsCounts){
		pop2HapsFreqs[pop2Hap.first] = pop2Hap.second/static_cast<double>(pop2Total);
	}

	for(const auto & pop1Hap : pop1Haps){
		if(!njh::in(pop1Hap.popUid_, pop2HapsCounts)){
			++ret.uniqueHapsInPop1_;
			ret.uniqueHapsInPop1CumFreq_ += pop2HapsFreqs[pop1Hap.popUid_];
		}
	}
	for(const auto & pop2Hap : pop2Haps){
		if(!njh::in(pop2Hap.popUid_, pop1HapsCounts)){
			++ret.uniqueHapsInPop2_;
			ret.uniqueHapsInPop2CumFreq_ += pop2HapsFreqs[pop2Hap.popUid_];
		}
		if(njh::in(pop2Hap.popUid_, pop1HapsCounts)){
			++ret.uniqueHapsShared_;
		}
	}
	ret.uniqueHapsAll_ = allHaps.size();

	//sorensen
	{
		ret.sorensenDistance_ = (ret.uniqueHapsInPop1_ + ret.uniqueHapsInPop2_)/static_cast<double>(ret.uniqueHapsInPop1_ + ret.uniqueHapsInPop2_ + ret.uniqueHapsShared_ * 2.0);
	}

	//jaccard
	{
		ret.jaccardIndexDissim_ = (ret.uniqueHapsInPop1_ + ret.uniqueHapsInPop2_)/static_cast<double>(ret.uniqueHapsInPop1_ + ret.uniqueHapsInPop2_ + ret.uniqueHapsShared_);
	}

	//brayCurtis
	{
		double sumOfAbsDiffs = 0;
		double total = 0;
		for(const auto & hap : allHaps){
			sumOfAbsDiffs += uAbsdiff(pop1HapsCounts[hap], pop2HapsCounts[hap]);
			total += pop1HapsCounts[hap] + pop2HapsCounts[hap];
		}
		ret.brayCurtisDissim_ = sumOfAbsDiffs/total;
	}

	//brayCurtis relative
	{
		double sumOfAbsDiffs = 0;
		double total = 0;
		for(const auto & hap : allHaps){
			sumOfAbsDiffs += std::abs(pop1HapsFreqs[hap] - pop2HapsFreqs[hap]);
			total += pop1HapsFreqs[hap] + pop2HapsFreqs[hap];
		}
		ret.brayCurtisRelativeDissim_ = sumOfAbsDiffs/total;
	}

	//RMSE
	{
		double sumOfSqaures = 0;
		for(const auto & hap : allHaps){
			sumOfSqaures += std::pow(pop1HapsFreqs[hap] - pop2HapsFreqs[hap], 2.0);
		}
		ret.RMSE_ = std::sqrt(sumOfSqaures/allHaps.size());
	}

	//half R and matching coefficien
	{
		//if lots and lots of rare haps this number might be skewed but all pops should be effected
		double a = 0; //found in both
		double b = 0; //found only in pop 1
		double c = 0; //found only in pop 2
		double d = 0; //not found in either
		for(const auto & hap : allPossibleHaps){
			//have to do > 0 cause index in above will create empty records
			bool inPop1 = pop1HapsCounts[hap] > 0;
			bool inPop2 = pop2HapsCounts[hap] > 0;
			if(inPop1 && inPop2){
				++a;
			}else if(inPop1 && !inPop2){
				++b;
			}else if(!inPop1 && inPop2){
				++c;
			}else{
				++d;
			}
		}

		//matching coefficient
		ret.matchingCoefficientDistance_ = (b + c)/(a + b + c + d);

		//half r
		//I feel like we have to add 1 to all
		++a;
		++b;
		++c;
		++d;
		ret.halfR_ = 0.5 * (1 - (((a*d) - (b *c))/std::sqrt((a + b) * (c + d) * (a + c) * (b + d))));

	}



	// discriminatingAvalance_
	if(pairwiseDistances.empty()){
		//no distances supplied
		ret.discriminatingAvalance_  = ret.plainAvalance_ ;
	}else{
		double sum = 0;
		for(const auto & hap1 : allHaps){
			for(const auto & hap2 : allHaps){
				sum += pairwiseDistances.at(hap1).at(hap2) * std::abs(pop1HapsFreqs[hap1] - pop2HapsFreqs[hap1]) * std::abs(pop1HapsFreqs[hap2] - pop2HapsFreqs[hap2]);
			}
		}
		ret.plainAvalance_ = sum * 0.5;
	}

	// plainAvalance_
	{
		double sum = 0;
		for(const auto & hap1 : allHaps){
			for(const auto & hap2 : allHaps){
				sum += std::abs(pop1HapsFreqs[hap1] - pop2HapsFreqs[hap1]) * std::abs(pop1HapsFreqs[hap2] - pop2HapsFreqs[hap2]);
			}
		}
		ret.plainAvalance_ = sum * 0.5;
	}


	return ret;
}



std::unordered_map<std::string,
		std::unordered_map<std::string, PopGenCalculator::PopDifferentiationMeasuresPairWise>> PopGenCalculator::getPairwisePopDiff(
		const std::unordered_map<std::string, std::vector<PopHapInfo>> & hapsForPopulations,
		const std::unordered_map<uint32_t, std::unordered_map<uint32_t, double>> & pairwiseDists) {
	if(hapsForPopulations.size() < 2){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << " error, popSeqs should at least be size 2 not " << hapsForPopulations.size() << "\n";
		throw std::runtime_error{ss.str()};
	}
	std::unordered_map<std::string, std::unordered_map<std::string, PopDifferentiationMeasuresPairWise>> ret;
	auto keys = njh::getVecOfMapKeys(hapsForPopulations);
	njh::sort(keys);
	std::unordered_set<uint32_t> allHaps;
	for(const auto & hapsForPopulation : hapsForPopulations){
		for(const auto & hap : hapsForPopulation.second){
			allHaps.emplace(hap.popUid_);
		}
	}
	for(const auto keyPos : iter::range(keys.size())){
		for(const auto secondKeyPos : iter::range(keyPos)){
			auto popMeasures = getPopDiff(
					keys[keyPos],       hapsForPopulations.at(keys[keyPos]),
					keys[secondKeyPos], hapsForPopulations.at(keys[secondKeyPos]),
					allHaps, pairwiseDists
					);
			ret[keys[keyPos]][keys[secondKeyPos]] = popMeasures;
		}
	}
	return ret;
}


}  // namespace njhseq

