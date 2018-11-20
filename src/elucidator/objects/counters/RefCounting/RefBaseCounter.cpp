/*
 * RefBaseCounter.cpp
 *
 *  Created on: Dec 24, 2015
 *      Author: nick
 */
//
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
#include "RefBaseCounter.hpp"

namespace njhseq {

void RefBaseCounter::BaseCount::increaseRates(substituteMatrix & mat, char base)const{
	counts_.at(true).at(true).increaseRates(mat, base);
	counts_.at(true).at(false).increaseRates(mat, base);
	counts_.at(false).at(true).increaseRates(mat, base);
	counts_.at(false).at(false).increaseRates(mat, base);
}

void RefBaseCounter::BaseCount::increaseRatesHq(substituteMatrix & mat, char base)const{
	counts_.at(true).at(true).increaseRates(mat, base);
	counts_.at(false).at(true).increaseRates(mat, base);
}

void RefBaseCounter::BaseCount::increaseRatesLq(substituteMatrix & mat, char base)const{
	counts_.at(true).at(false).increaseRates(mat, base);
	counts_.at(false).at(false).increaseRates(mat, base);
}

void RefBaseCounter::BaseCount::increaseRatesForward(substituteMatrix & mat, char base)const{
	counts_.at(true).at(true).increaseRates(mat, base);
	counts_.at(true).at(false).increaseRates(mat, base);
}

void RefBaseCounter::BaseCount::increaseRatesReverse(substituteMatrix & mat, char base)const{
	counts_.at(false).at(true).increaseRates(mat, base);
	counts_.at(false).at(false).increaseRates(mat, base);
}

void RefBaseCounter::BaseCount::increaseRatesHqForward(substituteMatrix & mat, char base)const{
	counts_.at(true).at(true).increaseRates(mat, base);
}

void RefBaseCounter::BaseCount::increaseRatesLqForward(substituteMatrix & mat, char base)const{
	counts_.at(true).at(false).increaseRates(mat, base);
}

void RefBaseCounter::BaseCount::increaseRatesHqReverse(substituteMatrix & mat, char base)const{
	counts_.at(false).at(true).increaseRates(mat, base);
}

void RefBaseCounter::BaseCount::increaseRatesLqReverse(substituteMatrix & mat, char base)const{
	counts_.at(false).at(false).increaseRates(mat, base);
}

void RefBaseCounter::BaseCount::merge(const BaseCount & otherCounter){
	counts_[true][true].addOtherCounts(otherCounter.counts_.at(true).at(true), false);
	counts_[true][false].addOtherCounts(otherCounter.counts_.at(true).at(false), false);
	counts_[false][true].addOtherCounts(otherCounter.counts_.at(false).at(true), false);
	counts_[false][false].addOtherCounts(otherCounter.counts_.at(false).at(false), false);
}

void RefBaseCounter::merge(const RefBaseCounter & otherCounter){
	for(auto & ref : otherCounter.counts_){
		for(auto & pos : ref.second){
			counts_[ref.first][pos.first].merge(pos.second);
		}
	}
}

RefBaseCounter::BaseCount::BaseCount() {
	std::set<char> alphabet { 'A', 'C', 'G', 'T' };
	counts_[true].emplace(true, DNABaseCounter(alphabet));
	counts_[true].emplace(false, DNABaseCounter(alphabet));
	counts_[false].emplace(true, DNABaseCounter(alphabet));
	counts_[false].emplace(false, DNABaseCounter(alphabet));
}

void RefBaseCounter::BaseCount::increaseCount(bool plusStrand, bool highQuality,
		char base, uint32_t count) {
	counts_[plusStrand][highQuality].increase(base, count);
}

void RefBaseCounter::BaseCount::increaseCount(bool plusStrand, char base,
		uint32_t count) {
	counts_[plusStrand][true].increase(base, count);
}

void RefBaseCounter::BaseCount::resetAlphabets() {
	counts_[true][true].resetAlphabet(true);
	counts_[true][false].resetAlphabet(true);
	counts_[false][true].resetAlphabet(true);
	counts_[false][false].resetAlphabet(true);
}

void RefBaseCounter::BaseCount::setFractions() {
	counts_[true][true].setFractions();
	counts_[true][false].setFractions();
	counts_[false][true].setFractions();
	counts_[false][false].setFractions();
}

std::vector<char> RefBaseCounter::BaseCount::getAlphabet() const {
	std::set<char> allLetsSet;
	for (const auto & let : counts_.at(true).at(true).alphabet_) {
		allLetsSet.emplace(let);
	}
	for (const auto & let : counts_.at(true).at(false).alphabet_) {
		allLetsSet.emplace(let);
	}
	for (const auto & let : counts_.at(false).at(true).alphabet_) {
		allLetsSet.emplace(let);
	}
	for (const auto & let : counts_.at(false).at(false).alphabet_) {
		allLetsSet.emplace(let);
	}
	return std::vector<char>(allLetsSet.begin(), allLetsSet.end());
}

uint64_t RefBaseCounter::BaseCount::getTotalCount() const {
	return counts_.at(true).at(true).getTotalCount()
			+ counts_.at(false).at(true).getTotalCount()
			+ counts_.at(true).at(false).getTotalCount()
			+ counts_.at(false).at(false).getTotalCount();
}

uint64_t RefBaseCounter::BaseCount::getHighQTotalCount() const {
	return counts_.at(true).at(true).getTotalCount()
			+ counts_.at(false).at(true).getTotalCount();
}

double RefBaseCounter::BaseCount::getHighQFrac() const {
	return getHighQTotalCount() / static_cast<double>(getTotalCount());
}

uint64_t RefBaseCounter::BaseCount::getLowQTotalCount() const {
	return counts_.at(true).at(false).getTotalCount()
			+ counts_.at(false).at(false).getTotalCount();
}

double RefBaseCounter::BaseCount::getLowQFrac() const {
	return getLowQTotalCount() / static_cast<double>(getTotalCount());
}

uint64_t RefBaseCounter::BaseCount::getPlusStrandTotalCount() const {
	return counts_.at(true).at(true).getTotalCount()
			+ counts_.at(true).at(false).getTotalCount();
}

double RefBaseCounter::BaseCount::getPlusStrandFrac() const {
	return getPlusStrandTotalCount() / static_cast<double>(getTotalCount());
}

uint64_t RefBaseCounter::BaseCount::getNegStrandTotalCount() const {
	return counts_.at(false).at(true).getTotalCount()
			+ counts_.at(false).at(false).getTotalCount();
}

double RefBaseCounter::BaseCount::getNegStrandFrac() const {
	return getNegStrandTotalCount() / static_cast<double>(getTotalCount());
}

uint64_t RefBaseCounter::BaseCount::getHighQPlusStrandTotalCount() const {
	return counts_.at(true).at(true).getTotalCount();
}

double RefBaseCounter::BaseCount::getHighQFracPlusStrand() const {
	return getHighQPlusStrandTotalCount()
			/ static_cast<double>(getPlusStrandTotalCount());
}

uint64_t RefBaseCounter::BaseCount::getHighQNegStrandTotalCount() const {
	return counts_.at(false).at(true).getTotalCount();
}

uint64_t RefBaseCounter::BaseCount::getLowQPlusStrandTotalCount() const {
	return counts_.at(true).at(false).getTotalCount();
}

double RefBaseCounter::BaseCount::getHighQFracNegStrand() const {
	return getHighQNegStrandTotalCount()
			/ static_cast<double>(getNegStrandTotalCount());
}

uint64_t RefBaseCounter::BaseCount::getLowQNegStrandTotalCount() const {
	return counts_.at(false).at(false).getTotalCount();
}

uint64_t RefBaseCounter::BaseCount::getCount(char base) const {
	return counts_.at(true).at(true).bases_[base]
			+ counts_.at(false).at(true).bases_[base]
			+ counts_.at(true).at(false).bases_[base]
			+ counts_.at(false).at(false).bases_[base];
}

double RefBaseCounter::BaseCount::getFrac(char base) const {
	return getCount(base) / static_cast<double>(getTotalCount());
}

uint64_t RefBaseCounter::BaseCount::getHighQCount(char base) const {
	return counts_.at(true).at(true).bases_[base]
			+ counts_.at(false).at(true).bases_[base];
}

uint64_t RefBaseCounter::BaseCount::getLowQCount(char base) const {
	return counts_.at(true).at(false).bases_[base]
			+ counts_.at(false).at(false).bases_[base];
}

double RefBaseCounter::BaseCount::getHighQFrac(char base) const {
	return getHighQCount(base)
			/ static_cast<double>(getHighQCount(base) + getLowQCount(base));
}

double RefBaseCounter::BaseCount::getBaseFracHighQ(char base) const {
	return getHighQCount(base) / static_cast<double>(getHighQTotalCount());
}

double RefBaseCounter::BaseCount::getBaseFracLowQ(char base) const {
	return getLowQCount(base) / static_cast<double>(getLowQTotalCount());
}

uint64_t RefBaseCounter::BaseCount::getPlusStrandCount(char base) const {
	return counts_.at(true).at(true).bases_[base]
			+ counts_.at(true).at(false).bases_[base];
}

uint64_t RefBaseCounter::BaseCount::getNegStrandCount(char base) const {
	return counts_.at(false).at(true).bases_[base]
			+ counts_.at(false).at(false).bases_[base];
}

double RefBaseCounter::BaseCount::getPlusStrandFrac(char base) const {
	return getPlusStrandCount(base)
			/ static_cast<double>(getPlusStrandCount(base) + getNegStrandCount(base));
}

uint64_t RefBaseCounter::BaseCount::getHighQPlusStrandCount(char base) const {
	return counts_.at(true).at(true).bases_[base];
}

uint64_t RefBaseCounter::BaseCount::getHighQNegStrandCount(char base) const {
	return counts_.at(false).at(true).bases_[base];
}

double RefBaseCounter::BaseCount::getPlusStrandFracHighQ(char base) const {
	return getHighQPlusStrandCount(base)
			/ static_cast<double>(getHighQPlusStrandCount(base)
					+ getHighQNegStrandCount(base));
}

uint64_t RefBaseCounter::BaseCount::getLowQPlusStrandCount(char base) const {
	return counts_.at(true).at(false).bases_[base];
}

uint64_t RefBaseCounter::BaseCount::getLowQNegStrandCount(char base) const {
	return counts_.at(false).at(false).bases_[base];
}

double RefBaseCounter::BaseCount::getPlusStrandFracLowQ(char base) const {
	return getLowQPlusStrandCount(base)
			/ static_cast<double>(getLowQPlusStrandCount(base)
					+ getLowQNegStrandCount(base));
}

std::string RefBaseCounter::BaseCount::infoStratifiedByQualityHeader(
		const std::string & delim) {
	return vectorToString(VecStr { "baseCountHighQ", "fStrandHighQ",
			"rStrandHighQ", "percForwardHighQ", "baseFractionHighQ", "totalCountHighQ",
			"highQualFraction", "baseCountLowQ", "fStrandLowQ",
			"rStrandLowQ", "percForwardLowQ", "baseFractionLowQ",
			"totalCountLowQ", "totalBaseCount" }, delim);
}

void RefBaseCounter::BaseCount::infoStratifiedByQuality(std::ostream & out, char base, const std::string delim)const{
	out << getHighQCount(base)
							<< delim << getHighQPlusStrandCount(base)
							<< delim << getHighQNegStrandCount(base)
							<< delim << getPlusStrandFracHighQ(base)
							<< delim << getBaseFracHighQ(base)
							<< delim << getHighQTotalCount()
							<< delim << getHighQFrac()
							<< delim << getLowQCount(base)
							<< delim << getLowQPlusStrandCount(base)
							<< delim << getLowQNegStrandCount(base)
							<< delim << getPlusStrandFracLowQ(base)
							<< delim << getBaseFracLowQ(base)
							<< delim << getLowQTotalCount()
							<< delim << getTotalCount();
}

std::string RefBaseCounter::BaseCount::infoHeader(const std::string & delim){
	return vectorToString(VecStr { "baseCount", "fStrand",
			"rStrand", "percForward", "baseFraction", "totalCount",
			"highQualFraction"}, delim);
}

void RefBaseCounter::BaseCount::info(std::ostream & out, char base, const std::string delim)const{
	out << getCount(base)
							<< delim << getPlusStrandCount(base)
							<< delim << getNegStrandCount(base)
							<< delim << getPlusStrandFrac(base)
							<< delim << getFrac(base)
							<< delim << getTotalCount()
							<< delim << getHighQFrac();
}

VecStr RefBaseCounter::getRefNames() const {
	auto refKeys = getVectorOfMapKeys(counts_);
	njh::sort(refKeys);
	return refKeys;
}

void RefBaseCounter::setAlphabetAndFractions() {
	if (needsToSet) {
		for (auto & ref : counts_) {
			for (auto & pos : ref.second) {
				pos.second.resetAlphabets();
				pos.second.setFractions();
			}
		}
		needsToSet = false;
	}
}

void RefBaseCounter::increaseCount(const std::string & refName, size_t refPos,
		bool plusStrand, bool highQuality, char base, uint32_t count) {
	counts_[refName][refPos].increaseCount(plusStrand, highQuality, base, count);
	needsToSet = true;
}

void RefBaseCounter::increaseCount(const std::string & refName, size_t refPos,
		bool plusStrand, char base, uint32_t count) {
	counts_[refName][refPos].increaseCount(plusStrand, base, count);
	needsToSet = true;
}

std::vector<size_t> RefBaseCounter::getPositionsForRef(
		const std::string & refName) const {
	auto search = counts_.find(refName);
	std::vector<size_t> ret;
	if (search != counts_.end()) {
		ret = getVectorOfMapKeys(search->second);
		njh::sort(ret);
	}
	return ret;
}



}  // namespace njhseq
