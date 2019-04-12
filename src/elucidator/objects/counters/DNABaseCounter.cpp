/*
 * DNABaseCounter.cpp
 *
 *  Created on: Feb 26, 2016
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
#include "DNABaseCounter.hpp"



namespace njhseq {

DNABaseCounter::DNABaseCounter() {
	bases_.fill(0);
	baseFracs_.fill(0);
}

DNABaseCounter::DNABaseCounter(std::set<char> alphabet) :
		alphabet_(alphabet), originalAlphabet_(alphabet) {
	bases_.fill(0);
	baseFracs_.fill(0);
}


Json::Value DNABaseCounter::toJson() const{
	Json::Value ret;
	ret["class"] = njh::json::toJson("njhseq::charCounter");
	auto & bases = ret["bases_"];
	for (auto let : alphabet_) {
		bases[std::string(1, let)] = bases_[let];
	}
	auto & fracs = ret["fractions_"];
	for (auto let : alphabet_) {
		fracs[std::string(1, let)] = baseFracs_[let];
	}

	ret["alphabet_"] = njh::json::toJson(alphabet_);
	ret["originalAlphabet_"] = njh::json::toJson(originalAlphabet_);

	ret["gcContent"] = njh::json::toJson(calcGcContent());
	ret["totalCount"] = njh::json::toJson(getTotalCount());

	return ret;
}

void DNABaseCounter::resetAlphabet(bool keepOld) {
	std::set<char> determineAlph;
	for (const auto base : bases_.bases()) {
		if (bases_[base] > 0) {
			determineAlph.insert(static_cast<char>(base));
		}
	}
	if (keepOld) {
		alphabet_.clear();
		std::set_union(originalAlphabet_.begin(), originalAlphabet_.end(),
				determineAlph.begin(), determineAlph.end(),
				std::inserter(alphabet_, alphabet_.begin()));
	} else {
		alphabet_ = determineAlph;
	}
}

void DNABaseCounter::reset() {
	for (const auto base : bases_.bases()) {
		bases_[base] = 0;
		baseFracs_[base] = 0;
	}
}

//sequences without qualities
void DNABaseCounter::increase(char base, uint32_t cnt) {
	bases_[base] += cnt;
}
void DNABaseCounter::increase(const std::string &seq, uint32_t cnt) {
	for (const auto c : seq) {
		increase(c, cnt);
	}
}

void DNABaseCounter::increaseWithSubStr(const std::string & str, size_t pos,
		size_t len, uint32_t cnt) {
	increase(str.substr(pos, len), cnt);
}

void DNABaseCounter::setFractions(){
	setFractions(alphabet_);
}
void DNABaseCounter::setFractions(const std::set<char> & alphabet){
	double total = getTotalCount(alphabet);
	if(total > 0){
		for(const auto c : alphabet){
			baseFracs_[c] = bases_[c]/total;
		}
	}else{
		for(const auto c : alphabet){
			baseFracs_[c] = 0;
		}
	}
}

void DNABaseCounter::addOtherCounts(const DNABaseCounter & otherCounter, bool setFractionsAfterwards){
	for(const auto c : alphabet_){
		bases_[c] += otherCounter.bases_[c];
	}
	if(setFractionsAfterwards){
		setFractions();
	}
}

uint32_t DNABaseCounter::getTotalCount() const{
	return getTotalCount(alphabet_);
}
uint32_t DNABaseCounter::getTotalCount(const std::set<char> & alphabet) const{
	uint32_t total = 0;
	for(const auto base : alphabet){
		total += bases_[base];
	}
	return total;
}

std::multimap<double, char, std::less<double>> DNABaseCounter::createLikelihoodMaps(
    bool setFractionFirst){
	if (setFractionFirst) {
		setFractions();
	}
	std::multimap<double, char, std::less<double>> likelihoods;
	for (const auto c : alphabet_) {
		likelihoods.emplace(baseFracs_[c], c);
	}
	return likelihoods;
}

uint32_t DNABaseCounter::getCountForBases(const std::vector<char> & bases) const{
	uint32_t total = 0;
	for(const auto base : bases){
		total += bases_[base];
	}
	return total;
}

//GC
uint32_t DNABaseCounter::getGcCount() const{
	return getCountForBases(std::vector<char>{'G', 'C', 'g' ,'c'});
}

double DNABaseCounter::calcGcContent()const{
	return (bases_['G'] + bases_['C'] + bases_['g'] + bases_['c'])/static_cast<double>(getTotalCount());
}
int DNABaseCounter::getGcDifference()const{
	return (bases_['G'] + bases_['g']) - (bases_['C']+ bases_['c']);
}
// compute entropy
double DNABaseCounter::computeEntrophy()const{
	double total = getTotalCount();
	double sum = 0;
	for (const auto c : alphabet_) {
		if (0 != bases_[c]/total) {
			if (1 == bases_[c]/total) {
				return 0;
			}
			sum += bases_[c]/total * std::log2(bases_[c]/total);
		}
	}
	return (-1 * sum);
}

char DNABaseCounter::getDegenativeBase() const{
	if (bases_['A'] == 0 && bases_['C'] == 0 && bases_['G'] == 0
			&& bases_['T'] > 0) {
		return 'T';
	} else if (bases_['A'] == 0 && bases_['C'] == 0 && bases_['G'] > 0
			&& bases_['T'] == 0) {
		return 'G';
	} else if (bases_['A'] == 0 && bases_['C'] > 0 && bases_['G'] == 0
			&& bases_['T'] == 0) {
		return 'C';
	} else if (bases_['A'] > 0 && bases_['C'] == 0 && bases_['G'] == 0
			&& bases_['T'] == 0) {
		return 'A';
	} else if (bases_['A'] == 0 && bases_['C'] == 0 && bases_['G'] > 0
			&& bases_['T'] > 0) {
		return 'K';
	} else if (bases_['A'] == 0 && bases_['C'] > 0 && bases_['G'] == 0
			&& bases_['T'] > 0) {
		return 'Y';
	} else if (bases_['A'] > 0 && bases_['C'] == 0 && bases_['G'] == 0
			&& bases_['T'] > 0) {
		return 'W';
	} else if (bases_['A'] == 0 && bases_['C'] > 0 && bases_['G'] > 0
			&& bases_['T'] == 0) {
		return 'S';
	} else if (bases_['A'] > 0 && bases_['C'] == 0 && bases_['G'] > 0
			&& bases_['T'] == 0) {
		return 'R';
	} else if (bases_['A'] > 0 && bases_['C'] > 0 && bases_['G'] == 0
			&& bases_['T'] == 0) {
		return 'M';
	} else if (bases_['A'] == 0 && bases_['C'] > 0 && bases_['G'] > 0
			&& bases_['T'] > 0) {
		return 'B';
	} else if (bases_['A'] > 0 && bases_['C'] == 0 && bases_['G'] > 0
			&& bases_['T'] > 0) {
		return 'D';
	} else if (bases_['A'] > 0 && bases_['C'] > 0 && bases_['G'] == 0
			&& bases_['T'] > 0) {
		return 'H';
	} else if (bases_['A'] > 0 && bases_['C'] > 0 && bases_['G'] > 0
			&& bases_['T'] == 0) {
		return 'V';
	} else {
		return 'N';
	}
}

std::string DNABaseCounter::outPutInfoHeader(){
	return njh::conToStr(toVecStr("base", "count", "fraction"), "\t");
}
// output data
void DNABaseCounter::outPutInfo(std::ostream &out, bool header) const{
	if(header){
		out << outPutInfoHeader() << std::endl;
	}
	for(auto c : alphabet_ ){
		out << njh::conToStr(toVecStr(c, bases_[c], baseFracs_[c]), "\t") << std::endl;
	}
}

double DNABaseCounter::getFracDifference(const DNABaseCounter & otherCounter, const std::vector<char> & alph)const{
	double sum = 0;
	for (const auto & let : alph) {
		sum += std::abs(baseFracs_[let] - otherCounter.baseFracs_[let]);
	}
	return sum;
}

void DNABaseCounter::increaseRates(substituteMatrix & mat, char refBase) const{
	for(char let : alphabet_){
		mat.mat_[refBase][let] += bases_[let];
	}
}



}  // namespace njhseq


