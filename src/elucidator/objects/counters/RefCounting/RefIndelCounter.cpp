/*
 * RefIndelCounter.cpp
 *
 *  Created on: Dec 24, 2015
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


#include "RefIndelCounter.hpp"


namespace njhseq {
void RefIndelCounter::IndelCounter::merge(const IndelCounter & otherCounter){
	counts_[true][true].addOtherCounts(otherCounter.counts_.at(true).at(true));
	counts_[true][false].addOtherCounts(otherCounter.counts_.at(true).at(false));
	counts_[false][true].addOtherCounts(otherCounter.counts_.at(false).at(true));
	counts_[false][false].addOtherCounts(otherCounter.counts_.at(false).at(false));
}

void RefIndelCounter::merge(const RefIndelCounter & otherCounter){
	for(auto & ref : otherCounter.counts_){
		for(auto & pos : ref.second){
			counts_[ref.first][pos.first].merge(pos.second);
		}
	}
}



std::string RefIndelCounter::IndelCounter::infoHeader(
		const std::string & delim) {
	return vectorToString(VecStr { "count", "fStrand", "rStrand", "percForward" },
			delim);
}

std::string RefIndelCounter::IndelCounter::sizeInfoHeader(
		const std::string & delim) {
	return vectorToString(VecStr { "refName", "refPos", "refBase", "base", "size",
			"count", "fCount", "rCount", "percForward", "fraction" }, delim);
}

std::string RefIndelCounter::IndelCounter::seqInfoHeader(
		const std::string & delim) {
	return vectorToString(VecStr { "refName", "refPos", "refBase", "base", "size",
			"count", "fCount", "rCount", "percForward", "fraction" }, delim);
}



void RefIndelCounter::IndelCounter::insertionSizeInfo(std::ostream & out,
		const std::string refName, size_t refPos, char refBase,
		const std::string delim) const {
	//strand (true = + strand), size, count
	std::unordered_map<size_t, std::unordered_map<bool, uint32_t>> sizeCounts;
	uint32_t total = 0;
	for (const auto & strand : counts_) {
		for (const auto & count : strand.second.at(true).counts_) {
			sizeCounts[count.first.size()][strand.first] += count.second;
			sizeCounts[count.first.size()][!strand.first] += 0;
			total += count.second;
		}
	}
	for(const auto & size : sizeCounts){
		auto currentTotal = size.second.at(true) + size.second.at(false);
		out << refName
				<< delim << refPos
				<< delim << refBase
				<< delim << "INS"
				<< delim << size.first
				<< delim << currentTotal
				<< delim << size.second.at(true)
				<< delim << size.second.at(false)
				<< delim << size.second.at(true)/static_cast<double>(currentTotal)
				<< delim << currentTotal/static_cast<double>(total)
				<< std::endl;;
	}
	if (sizeCounts.empty()) {
		out << refName << delim << refPos << delim << refBase << delim << "INS"
				<< delim << 0 << delim << 0 << delim << 0 << delim << 0 << delim << 0
				<< delim << 0 << std::endl;
	}
}


void RefIndelCounter::IndelCounter::insertionSeqInfo(std::ostream & out,
		const std::string refName, size_t refPos, char refBase,
		const std::string delim) const {
	//strand (true = + strand), size, count
	std::unordered_map<std::string, std::unordered_map<bool, uint32_t>> seqCounts;
	uint32_t total = 0;
	for (const auto & strand : counts_) {
		for (const auto & count : strand.second.at(true).counts_) {
			seqCounts[count.first][strand.first] += count.second;
			seqCounts[count.first][!strand.first] += 0;
			total += count.second;
		}
	}
	for(const auto & seq : seqCounts){
		auto currentTotal = seq.second.at(true) + seq.second.at(false);
		out << refName
				<< delim << refPos
				<< delim << refBase
				<< delim << "INS"
				<< delim << seq.first
				<< delim << currentTotal
				<< delim << seq.second.at(true)
				<< delim << seq.second.at(false)
				<< delim << seq.second.at(true)/static_cast<double>(currentTotal)
				<< delim << currentTotal/static_cast<double>(total)
				<< std::endl;
	}
	if (seqCounts.empty()) {
		out << refName << delim << refPos << delim << refBase << delim << "INS"
				<< delim << "" << delim << 0 << delim << 0 << delim << 0 << delim << 0
				<< delim << 0 << std::endl;
	}
}

void RefIndelCounter::IndelCounter::deletionSizeInfo(std::ostream & out,
		const std::string refName, size_t refPos, char refBase,
		const std::string delim) const {
	//strand (true = + strand), size, count
	std::unordered_map<size_t, std::unordered_map<bool, uint32_t>> sizeCounts;
	uint32_t total = 0;
	for (const auto & strand : counts_) {
		for (const auto & count : strand.second.at(false).counts_) {
			sizeCounts[count.first.size()][strand.first] += count.second;
			sizeCounts[count.first.size()][!strand.first] += 0;
			total += count.second;
		}
	}

	for(const auto & size : sizeCounts){
		auto currentTotal = size.second.at(true) + size.second.at(false);
		out << refName
				<< delim << refPos
				<< delim << refBase
				<< delim << "DEL"
				<< delim << size.first
				<< delim << currentTotal
				<< delim << size.second.at(true)
				<< delim << size.second.at(false)
				<< delim << size.second.at(true)/static_cast<double>(currentTotal)
				<< delim << currentTotal/static_cast<double>(total)
				<< std::endl;
	}
	if (sizeCounts.empty()) {
		out << refName << delim << refPos << delim << refBase << delim << "DEL"
				<< delim << "" << delim << 0 << delim << 0 << delim << 0 << delim << 0
				<< delim << 0 << std::endl;
	}
}

void RefIndelCounter::IndelCounter::deletionSeqInfo(std::ostream & out,
		const std::string refName, size_t refPos, char refBase,
		const std::string delim) const {
	//strand (true = + strand), size, count
	std::unordered_map<std::string, std::unordered_map<bool, uint32_t>> seqCounts;
	uint32_t total = 0;
	for (const auto & strand : counts_) {
		for (const auto & count : strand.second.at(false).counts_) {
			seqCounts[count.first][strand.first] += count.second;
			seqCounts[count.first][!strand.first] += 0;
			total += count.second;
		}
	}
	for(const auto & seq : seqCounts){
		auto currentTotal = seq.second.at(true) + seq.second.at(false);
		out << refName
				<< delim << refPos
				<< delim << refBase
				<< delim << "DEL"
				<< delim << seq.first
				<< delim << currentTotal
				<< delim << seq.second.at(true)
				<< delim << seq.second.at(false)
				<< delim << seq.second.at(true)/static_cast<double>(currentTotal)
				<< delim << currentTotal/static_cast<double>(total)
				<< std::endl;;
	}
	if (seqCounts.empty()) {
		out << refName << delim << refPos << delim << refBase << delim << "DEL"
				<< delim << 0 << delim << 0 << delim << 0 << delim << 0 << delim << 0
				<< delim << 0 << std::endl;
	}
}

RefIndelCounter::IndelCounter::IndelCounter() {
	counts_[true][true] = strCounter();
	counts_[true][false] = strCounter();
	counts_[false][true] = strCounter();
	counts_[false][false] = strCounter();
}

void RefIndelCounter::IndelCounter::increaseCount(bool plusStrand,
		bool insertion, const std::string & gapStr, uint32_t count ) {
	counts_[plusStrand][insertion].increaseCountByString(gapStr, count);
}

uint32_t RefIndelCounter::IndelCounter::getInsTotal() const {
	return counts_.at(true).at(true).getTotalCount()
			+ counts_.at(false).at(true).getTotalCount();
}

uint32_t RefIndelCounter::IndelCounter::getInsTotalPlusStrand() const {
	return counts_.at(true).at(true).getTotalCount();
}

uint32_t RefIndelCounter::IndelCounter::getInsTotalNegStrand() const {
	return counts_.at(false).at(true).getTotalCount();
}

double RefIndelCounter::IndelCounter::getInsFracPlusStrand() const {
	return getInsTotalPlusStrand() / static_cast<double>(getInsTotal());
}

uint32_t RefIndelCounter::IndelCounter::getDelTotal() const {
	return counts_.at(true).at(false).getTotalCount()
			+ counts_.at(false).at(false).getTotalCount();
}

uint32_t RefIndelCounter::IndelCounter::getDelTotalPlusStrand() const {
	return counts_.at(true).at(false).getTotalCount();
}

uint32_t RefIndelCounter::IndelCounter::getDelTotalNegStrand() const {
	return counts_.at(false).at(false).getTotalCount();
}

double RefIndelCounter::IndelCounter::getDelFracPlusStrand() const {
	return getDelTotalPlusStrand() / static_cast<double>(getDelTotal());
}



void RefIndelCounter::IndelCounter::insertionInfo(std::ostream & out,
		const std::string delim) const {
	out << getInsTotal() << delim << getInsTotalPlusStrand() << delim
			<< getInsTotalNegStrand() << delim << getInsFracPlusStrand();
}



void RefIndelCounter::IndelCounter::deletionInfo(std::ostream & out,
		const std::string delim) const {
	out << getDelTotal() << delim << getDelTotalPlusStrand() << delim
			<< getDelTotalNegStrand() << delim << getDelFracPlusStrand();
}

void RefIndelCounter::increaseCount(const std::string & refName, size_t refPos,
		bool plusStrand, bool insertion, const std::string & gapStr,
		uint32_t count ) {
	counts_[refName][refPos].increaseCount(plusStrand, insertion, gapStr,
			count);
}

VecStr RefIndelCounter::getRefNames() const {
	auto refKeys = getVectorOfMapKeys(counts_);
	njh::sort(refKeys);
	return refKeys;
}

std::vector<size_t> RefIndelCounter::getPositionsForRef(const std::string & refName)const {
	auto search = counts_.find(refName);
	std::vector<size_t> ret;
	if(search != counts_.end()) {
		getVectorOfMapKeys(search->second);
		njh::sort(ret);
	}
	return ret;
}


}  // namespace njhseq
