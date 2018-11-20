/*
 * RefCounter.cpp
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
#include "RefCounter.hpp"

namespace njhseq {

RefCounter::RefCounter(const std::string & twoBitFilename) :
		reader_(twoBitFilename) {
}
RefCounter::RefCounter(const RefCounter & that) :
		reader_(that.reader_.getFilename()), hardCutOff_(that.hardCutOff_), baseCounter_(
				that.baseCounter_), indelCounter_(that.indelCounter_) {

}


void RefCounter::setHardCutOff(uint32_t hardCutOff){
	std::lock_guard<std::mutex> lock(mut_);
	hardCutOff_ = hardCutOff;
}

void RefCounter::increaseBaseCount(const std::string & refName, size_t refPos,
		bool plusStrand, bool highQuality, char base, uint32_t count) {
	std::lock_guard<std::mutex> lock(mut_);
	baseCounter_.increaseCount(refName, refPos, plusStrand, highQuality, base,
			count);
}

void RefCounter::increaseBaseCount(const std::string & refName, size_t refPos,
		bool plusStrand, char base, uint32_t count) {
	std::lock_guard<std::mutex> lock(mut_);
	baseCounter_.increaseCount(refName, refPos, plusStrand, base, count);
}

void RefCounter::increaseIndelCount(const std::string & refName, size_t refPos,
		bool plusStrand, bool insertion, const std::string & gapStr,
		uint32_t count) {
	std::lock_guard<std::mutex>lock(mut_);
	indelCounter_.increaseCount(refName, refPos, plusStrand, insertion, gapStr,
			count);
}

void RefCounter::merge(const RefCounter & otherCounter){
	indelCounter_.merge(otherCounter.indelCounter_);
	baseCounter_.merge(otherCounter.baseCounter_);
}

VecStr RefCounter::getRefNames() const {
	auto refKeys = getVectorOfMapKeys(baseCounter_.counts_);
	addOtherVec(refKeys, getVectorOfMapKeys(indelCounter_.counts_));
	removeDuplicates(refKeys);
	return refKeys;
}

std::vector<size_t> RefCounter::getPositionsForRef(const std::string & refName) const {
	std::vector<size_t> ret;
	auto baseSearch = baseCounter_.counts_.find(refName);
	if (baseSearch != baseCounter_.counts_.end()) {
		ret = getVectorOfMapKeys(baseSearch->second);
	}
	auto indelSearch = indelCounter_.counts_.find(refName);
	if (indelSearch != indelCounter_.counts_.end()) {
		addOtherVec(ret, getVectorOfMapKeys(indelSearch->second));
	}
	removeDuplicates(ret);
	return ret;
}

std::string RefCounter::printInfoHeader() {
	return "refName\trefPos\trefBase\tbase\t"
			+ RefBaseCounter::BaseCount::infoHeader("\t");
}

std::string RefCounter::printInfoStratifiedByQualHeader() {
	return "refName\trefPos\trefBase\tbase\t"
			+ RefBaseCounter::BaseCount::infoStratifiedByQualityHeader("\t");
}

void RefCounter::printIndelSeqPosInfo(std::ostream & out,
		const std::string & refKey, size_t posKey) {
	auto & baseRef = baseCounter_.counts_[refKey];
	auto & basePosCounts = baseRef[posKey];
	auto & indelRef = indelCounter_.counts_[refKey];
	auto & indelPosCounts = indelRef[posKey];
	if (basePosCounts.getTotalCount() + indelPosCounts.getInsTotal()
			+ indelPosCounts.getDelTotal() <= hardCutOff_) {
		return;
	}
	indelPosCounts.insertionSeqInfo(out, refKey, posKey,
			referencesCache_[refKey][posKey], "\t");
	indelPosCounts.deletionSeqInfo(out, refKey, posKey,
			referencesCache_[refKey][posKey], "\t");
}

void RefCounter::printIndelSizePosInfo(std::ostream & out,
		const std::string & refKey, size_t posKey) {
	auto & baseRef = baseCounter_.counts_[refKey];
	auto & basePosCounts = baseRef[posKey];
	auto & indelRef = indelCounter_.counts_[refKey];
	auto & indelPosCounts = indelRef[posKey];
	if (basePosCounts.getTotalCount() + indelPosCounts.getInsTotal()
			+ indelPosCounts.getDelTotal() <= hardCutOff_) {
		return;
	}
	indelPosCounts.insertionSizeInfo(out, refKey, posKey,
			referencesCache_[refKey][posKey], "\t");
	indelPosCounts.deletionSizeInfo(out, refKey, posKey,
			referencesCache_[refKey][posKey], "\t");
}

void RefCounter::printPosInfo(std::ostream & out, const std::string & refKey,
		size_t posKey) {

	auto & baseRef = baseCounter_.counts_[refKey];
	auto & indelRef = indelCounter_.counts_[refKey];
	auto & basePosCounts = baseRef[posKey];
	if (basePosCounts.getTotalCount() <= hardCutOff_) {
		return;
	}
	auto & indelPosCounts = indelRef[posKey];
	std::vector<char> allLets = basePosCounts.getAlphabet();
	for(auto & let : allLets){
		out << refKey
				<< "\t" << posKey
				<< "\t" << referencesCache_[refKey][posKey]
				<< "\t" << let
				<< "\t";
		basePosCounts.info(out, let, "\t");
		out << std::endl;
	}
	out << refKey
			<< "\t" << posKey
			<< "\t" << referencesCache_[refKey][posKey]
			<< "\t" << "INS"
			<< "\t";
	indelPosCounts.insertionInfo(out, "\t");
	out << "\t" << 0
			<< "\t" << basePosCounts.getTotalCount()
			<< "\t" << basePosCounts.getHighQFrac()
			<< std::endl;
	out << refKey
			<< "\t" << posKey
			<< "\t" << referencesCache_[refKey][posKey]
			<< "\t" << "DEL"
			<< "\t";
	indelPosCounts.deletionInfo(out, "\t");
	out << "\t" << 0
			<< "\t" << basePosCounts.getTotalCount()
			<< "\t" << basePosCounts.getHighQFrac()
			<< std::endl;
}

void RefCounter::printPosInfoStratifiedByQuality(std::ostream & out, const std::string & refKey,
		size_t posKey) {
	auto & baseRef = baseCounter_.counts_[refKey];
	auto & indelRef = indelCounter_.counts_[refKey];
	auto & basePosCounts = baseRef[posKey];
	auto & indelPosCounts = indelRef[posKey];
	if (basePosCounts.getTotalCount() <= hardCutOff_) {
		return;
	}
	std::vector<char> allLets = basePosCounts.getAlphabet();
	for(auto & let : allLets){
		out << refKey
				<< "\t" << posKey
				<< "\t" << referencesCache_[refKey][posKey]
				<< "\t" << let
				<< "\t";
		basePosCounts.infoStratifiedByQuality(out, let, "\t");
		out << std::endl;
	}
	out << refKey
			<< "\t" << posKey
			<< "\t" << referencesCache_[refKey][posKey]
			<< "\t" << "INS"
			<< "\t";
	indelPosCounts.insertionInfo(out, "\t");
	out << "\t" << 0
			<< "\t" << basePosCounts.getHighQTotalCount()
			<< "\t" << basePosCounts.getHighQFrac()
			<< "\t" << "\t" << "\t" << "\t" << "\t"
			<< "\t" << "\t" << std::endl;
	out << refKey
			<< "\t" << posKey
			<< "\t" << referencesCache_[refKey][posKey]
			<< "\t" << "DEL"
			<< "\t";
	indelPosCounts.deletionInfo(out, "\t");
	out << "\t" << 0
			<< "\t" << basePosCounts.getHighQTotalCount()
			<< "\t" << basePosCounts.getHighQFrac()
			<< "\t" << "\t" << "\t" << "\t" << "\t"
			<< "\t" << "\t" << std::endl;
}

void RefCounter::printPosInfo(std::ostream & out, const std::string & refKey,
		size_t posKey, size_t depthCutOff, double freqCutOff) {

	auto & baseRef = baseCounter_.counts_[refKey];
	auto & basePosCounts = baseRef[posKey];
	if (basePosCounts.getTotalCount() >= depthCutOff) {
		uint32_t aboveCutOff = 0;
		std::vector<char> allLets = basePosCounts.getAlphabet();
		for (auto & let : allLets) {
			if (basePosCounts.getFrac(let) > freqCutOff) {
				++aboveCutOff;
			}
		}
		if (aboveCutOff > 1) {
			printPosInfo(out, refKey, posKey);
		}
	}
}

void RefCounter::printPosInfoStratifiedByQuality(std::ostream & out,
		const std::string & refKey, size_t posKey, size_t depthCutOff,
		double freqCutOff) {
	auto & baseRef = baseCounter_.counts_[refKey];
	auto & basePosCounts = baseRef[posKey];
	if (basePosCounts.getHighQTotalCount() >= depthCutOff) {
		uint32_t aboveCutOff = 0;
		std::vector<char> allLets = basePosCounts.getAlphabet();
		for (auto & let : allLets) {
			if (basePosCounts.getBaseFracHighQ(let) > freqCutOff) {
				++aboveCutOff;
			}
		}
		if (aboveCutOff > 1) {
			printPosInfoStratifiedByQuality(out, refKey, posKey);
		}
	}
}

void RefCounter::printInfo(std::ostream & out) {
	std::lock_guard<std::mutex> lock(mut_);
	out << printInfoHeader() << std::endl;
	baseCounter_.setAlphabetAndFractions();
	auto refKeys = getRefNames();
	for (auto & refKey : refKeys) {
		//get reference, load sequence if necessary
		if (referencesCache_.find(refKey) == referencesCache_.end()) {
			std::string temp;
			reader_[refKey]->getSequence(temp);
			referencesCache_[refKey] = temp;
		}
		//get positions
		auto posKeys = getPositionsForRef(refKey);
		for (auto & posKey : posKeys) {
			printPosInfo(out, refKey, posKey);
		}
	}
}

void RefCounter::printInfoStratifiedByQual(std::ostream & out) {
	std::lock_guard<std::mutex> lock(mut_);
	out << printInfoStratifiedByQualHeader() << std::endl;
	baseCounter_.setAlphabetAndFractions();
	auto refKeys = getRefNames();
	for (auto & refKey : refKeys) {
		//get reference, load sequence if necessary
		if (referencesCache_.find(refKey) == referencesCache_.end()) {
			std::string temp;
			reader_[refKey]->getSequence(temp);
			referencesCache_[refKey] = temp;
		}
		//get positions
		auto posKeys = getPositionsForRef(refKey);
		for (auto & posKey : posKeys) {
			printPosInfoStratifiedByQuality(out, refKey, posKey);
		}
	}
}

void RefCounter::printInfo(std::ostream & out, const std::string & refKey,
		size_t start, size_t stop) {
	std::lock_guard<std::mutex> lock(mut_);
	out << printInfoHeader() << std::endl;
	baseCounter_.setAlphabetAndFractions();
	auto refKeys = getRefNames();
	//get reference, load sequence if necessary
	if (referencesCache_.find(refKey) == referencesCache_.end()) {
		std::string temp;
		reader_[refKey]->getSequence(temp);
		referencesCache_[refKey] = temp;
	}
	//get positions
	auto posKeys = getPositionsForRef(refKey);
	for (auto & posKey : posKeys) {
		if (posKey >= start && posKey < stop) {
			printPosInfo(out, refKey, posKey);
		}
	}
}

void RefCounter::printInfoStratifiedByQual(std::ostream & out,
		const std::string & refKey, size_t start, size_t stop) {
	std::lock_guard<std::mutex> lock(mut_);
	out << printInfoStratifiedByQualHeader() << std::endl;
	baseCounter_.setAlphabetAndFractions();
	auto refKeys = getRefNames();
	//get reference, load sequence if necessary
	if (referencesCache_.find(refKey) == referencesCache_.end()) {
		std::string temp;
		reader_[refKey]->getSequence(temp);
		referencesCache_[refKey] = temp;
	}
	//get positions
	auto posKeys = getPositionsForRef(refKey);
	for (auto & posKey : posKeys) {
		if (posKey >= start && posKey < stop) {
			printPosInfoStratifiedByQuality(out, refKey, posKey);
		}
	}
}

void RefCounter::printInfo(std::ostream & out, size_t depthCutOff,
		double freqCutOff) {
	std::lock_guard<std::mutex> lock(mut_);
	out << printInfoHeader() << std::endl;
	baseCounter_.setAlphabetAndFractions();
	auto refKeys = getRefNames();
	for (auto & refKey : refKeys) {
		//get reference, load sequence if necessary
		if (referencesCache_.find(refKey) == referencesCache_.end()) {
			std::string temp;
			reader_[refKey]->getSequence(temp);
			referencesCache_[refKey] = temp;
		}
		//get positions
		auto posKeys = getPositionsForRef(refKey);
		for (auto & posKey : posKeys) {
			printPosInfo(out, refKey, posKey, depthCutOff, freqCutOff);
		}
	}
}

void RefCounter::printInfoStratifiedByQual(std::ostream & out,
		size_t depthCutOff, double freqCutOff) {
	std::lock_guard<std::mutex> lock(mut_);
	out << printInfoStratifiedByQualHeader() << std::endl;
	baseCounter_.setAlphabetAndFractions();
	auto refKeys = getRefNames();
	for (auto & refKey : refKeys) {
		//get reference, load sequence if necessary
		if (referencesCache_.find(refKey) == referencesCache_.end()) {
			std::string temp;
			reader_[refKey]->getSequence(temp);
			referencesCache_[refKey] = temp;
		}
		//get positions
		auto posKeys = getPositionsForRef(refKey);
		for (auto & posKey : posKeys) {
			printPosInfoStratifiedByQuality(out, refKey, posKey, depthCutOff,
					freqCutOff);
		}
	}
}

void RefCounter::printInfo(std::ostream & out, const std::string & refKey,
		size_t start, size_t stop, size_t depthCutOff, double freqCutOff) {
	std::lock_guard<std::mutex> lock(mut_);
	out << printInfoHeader() << std::endl;
	baseCounter_.setAlphabetAndFractions();
	auto refKeys = getRefNames();
	//get reference, load sequence if necessary
	if (referencesCache_.find(refKey) == referencesCache_.end()) {
		std::string temp;
		reader_[refKey]->getSequence(temp);
		referencesCache_[refKey] = temp;
	}
	//get positions
	auto posKeys = getPositionsForRef(refKey);
	for (auto & posKey : posKeys) {
		if (posKey >= start && posKey < stop) {
			printPosInfo(out, refKey, posKey, depthCutOff, freqCutOff);
		}
	}
}

void RefCounter::printInfoStratifiedByQual(std::ostream & out,
		const std::string & refKey, size_t start, size_t stop, size_t depthCutOff,
		double freqCutOff) {
	std::lock_guard<std::mutex> lock(mut_);
	out << printInfoStratifiedByQualHeader() << std::endl;
	baseCounter_.setAlphabetAndFractions();
	auto refKeys = getRefNames();
	//get reference, load sequence if necessary
	if (referencesCache_.find(refKey) == referencesCache_.end()) {
		std::string temp;
		reader_[refKey]->getSequence(temp);
		referencesCache_[refKey] = temp;
	}
	//get positions
	auto posKeys = getPositionsForRef(refKey);
	for (auto & posKey : posKeys) {
		if (posKey >= start && posKey < stop) {
			printPosInfoStratifiedByQuality(out, refKey, posKey, depthCutOff,
					freqCutOff);
		}
	}
}

std::string RefCounter::printIndelSeqInfoHeader(){
	return RefIndelCounter::IndelCounter::seqInfoHeader("\t");
}

std::string RefCounter::printIndelSizeInfoHeader(){
	return RefIndelCounter::IndelCounter::sizeInfoHeader("\t");
}


void RefCounter::printIndelSeqInfo(std::ostream & out){
	std::lock_guard<std::mutex> lock(mut_);
	out << printIndelSeqInfoHeader() << std::endl;
	auto refKeys = getRefNames();
	for (auto & refKey : refKeys) {
		//get reference, load sequence if necessary
		if (referencesCache_.find(refKey) == referencesCache_.end()) {
			std::string temp;
			reader_[refKey]->getSequence(temp);
			referencesCache_[refKey] = temp;
		}
		//get positions
		auto posKeys = getPositionsForRef(refKey);
		for (auto & posKey : posKeys) {
			printIndelSeqPosInfo(out, refKey, posKey);
		}
	}
}

void RefCounter::printIndelSeqInfo(std::ostream & out,
		const std::string & refKey, size_t start, size_t stop) {
	std::lock_guard<std::mutex> lock(mut_);
	out << printIndelSeqInfoHeader() << std::endl;
	auto refKeys = getRefNames();
	//get reference, load sequence if necessary
	if (referencesCache_.find(refKey) == referencesCache_.end()) {
		std::string temp;
		reader_[refKey]->getSequence(temp);
		referencesCache_[refKey] = temp;
	}
	//get positions
	auto posKeys = getPositionsForRef(refKey);
	for (auto & posKey : posKeys) {
		if (posKey >= start && posKey < stop) {
			printIndelSeqPosInfo(out, refKey, posKey);
		}
	}
}

void RefCounter::printIndelSizeInfo(std::ostream & out){
	std::lock_guard<std::mutex> lock(mut_);
	out << printIndelSizeInfoHeader() << std::endl;
	auto refKeys = getRefNames();
	for (auto & refKey : refKeys) {
		//get reference, load sequence if necessary
		if (referencesCache_.find(refKey) == referencesCache_.end()) {
			std::string temp;
			reader_[refKey]->getSequence(temp);
			referencesCache_[refKey] = temp;
		}
		//get positions
		auto posKeys = getPositionsForRef(refKey);
		for (auto & posKey : posKeys) {
			printIndelSizePosInfo(out, refKey, posKey);
		}
	}
}

void RefCounter::printIndelSizeInfo(std::ostream & out, const std::string & refKey,
		size_t start, size_t stop){
	std::lock_guard<std::mutex> lock(mut_);
	out << printIndelSizeInfoHeader() << std::endl;
	auto refKeys = getRefNames();
	//get reference, load sequence if necessary
	if (referencesCache_.find(refKey) == referencesCache_.end()) {
		std::string temp;
		reader_[refKey]->getSequence(temp);
		referencesCache_[refKey] = temp;
	}
	//get positions
	auto posKeys = getPositionsForRef(refKey);
	for (auto & posKey : posKeys) {
		if (posKey >= start && posKey < stop) {
			printIndelSizePosInfo(out, refKey, posKey);
		}
	}
}

substituteMatrix RefCounter::getRatesAgainstRef() {
	std::lock_guard<std::mutex> lock(mut_);
	substituteMatrix ret;
	ret.setWithZeros();
	/*
	for(size_t i = 0; i < ret.mat_.size(); ++i){
		for(size_t j = 0; i < ret.mat_[i].size(); ++j){
			ret.mat_[i][j] = 0;
		}
	}*/
	auto refKeys = getRefNames();
	for (auto & refKey : refKeys) {
		//get reference, load sequence if necessary
		if (referencesCache_.find(refKey) == referencesCache_.end()) {
			std::string temp;
			reader_[refKey]->getSequence(temp);
			referencesCache_[refKey] = temp;
		}
		//get positions
		auto posKeys = getPositionsForRef(refKey);
		for (auto & posKey : posKeys) {
			auto & baseRef = baseCounter_.counts_[refKey];
			auto & basePosCounts = baseRef[posKey];
			if (basePosCounts.getTotalCount() > hardCutOff_) {
				baseCounter_.counts_.at(refKey).at(posKey).increaseRates(ret, referencesCache_[refKey][posKey]);
			}
		}
	}
	return ret;
}

substituteMatrix RefCounter::getRatesAgainstRef(const std::string & refKey,
		size_t start, size_t stop)  {
	std::lock_guard<std::mutex> lock(mut_);
	substituteMatrix ret;
	auto refKeys = getRefNames();
	//get reference, load sequence if necessary
	if (referencesCache_.find(refKey) == referencesCache_.end()) {
		std::string temp;
		reader_[refKey]->getSequence(temp);
		referencesCache_[refKey] = temp;
	}
	//get positions
	auto posKeys = getPositionsForRef(refKey);
	for (auto & posKey : posKeys) {
		if (posKey >= start && posKey < stop) {
			auto & baseRef = baseCounter_.counts_[refKey];
			auto & basePosCounts = baseRef[posKey];
			if (basePosCounts.getTotalCount() > hardCutOff_) {
				baseCounter_.counts_.at(refKey).at(posKey).increaseRates(ret, referencesCache_[refKey][posKey]);
			}
		}
	}
	return ret;
}

substituteMatrix RefCounter::getRatesAgainstRefHq(){
	std::lock_guard<std::mutex> lock(mut_);
	substituteMatrix ret;
	ret.setWithZeros();
	auto refKeys = getRefNames();
	for (auto & refKey : refKeys) {
		//get reference, load sequence if necessary
		if (referencesCache_.find(refKey) == referencesCache_.end()) {
			std::string temp;
			reader_[refKey]->getSequence(temp);
			referencesCache_[refKey] = temp;
		}
		//get positions
		auto posKeys = getPositionsForRef(refKey);
		for (auto & posKey : posKeys) {
			auto & baseRef = baseCounter_.counts_[refKey];
			auto & basePosCounts = baseRef[posKey];
			if (basePosCounts.getTotalCount() > hardCutOff_) {
				baseCounter_.counts_.at(refKey).at(posKey).increaseRatesHq(ret, referencesCache_[refKey][posKey]);
			}
		}
	}
	return ret;
}

substituteMatrix RefCounter::getRatesAgainstRefHq(const std::string & refKey,
		size_t start, size_t stop){
	std::lock_guard<std::mutex> lock(mut_);
	substituteMatrix ret;
	auto refKeys = getRefNames();
	//get reference, load sequence if necessary
	if (referencesCache_.find(refKey) == referencesCache_.end()) {
		std::string temp;
		reader_[refKey]->getSequence(temp);
		referencesCache_[refKey] = temp;
	}
	//get positions
	auto posKeys = getPositionsForRef(refKey);
	for (auto & posKey : posKeys) {
		if (posKey >= start && posKey < stop) {
			auto & baseRef = baseCounter_.counts_[refKey];
			auto & basePosCounts = baseRef[posKey];
			if (basePosCounts.getTotalCount() > hardCutOff_) {
				baseCounter_.counts_.at(refKey).at(posKey).increaseRatesHq(ret, referencesCache_[refKey][posKey]);
			}
		}
	}
	return ret;
}

void RefCounter::resetCounts(){
	baseCounter_ = RefBaseCounter();
	indelCounter_ = RefIndelCounter();
}

}  // namespace njhseq
