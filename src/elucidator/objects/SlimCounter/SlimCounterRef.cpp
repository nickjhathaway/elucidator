/*
 * SlimCounterRef.cpp
 *
 *  Created on: Jun 16, 2016
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

#include "SlimCounterRef.hpp"

#include <njhseq/utils/utils.hpp>


namespace njhseq {

SlimCounterRef::SlimCounterRef(const uint32_t & refId, bfs::path outputFile) :
		refId_(refId), outputFile_(njh::appendAsNeededRet(outputFile.string(), ".bgz")) {
}


void SlimCounterRef::initiateResetPositions(uint32_t start, uint32_t end){
	std::unique_lock<std::shared_timed_mutex> lock(mut_);
	if(end < start){
		std::stringstream ss;
		ss << __FILE__ << ": " << __LINE__ << " : " << __PRETTY_FUNCTION__ << " Error, end shouldn't be less than start" << "\n";
		throw std::runtime_error{ss.str()};
	}
	for(const auto pos : iter::range(start, end)){
		counts_[pos] = SlimCounterPos();
	}
}

void SlimCounterRef::increaseCountLockFree(uint32_t pos, char base, bool reverseStand, bool highQuality,
			uint32_t count){
	counts_[pos].increaseCount(base, reverseStand, highQuality, count);
}

void SlimCounterRef::increaseInsertionCountLockFree(uint32_t pos, bool reverseStand, uint32_t count){
	counts_[pos].increaseInsertionCount(reverseStand, count);
}

void SlimCounterRef::increaseDeletionCountLockFree(uint32_t pos, bool reverseStand, uint32_t count){
	counts_[pos].increaseDeletionCount(reverseStand, count);
}

uint32_t SlimCounterRef::maxPos() {
	std::shared_lock<std::shared_timed_mutex> lock(mut_);
	auto keys = getVectorOfMapKeys(counts_);
	return *std::max_element(keys.begin(), keys.end());
}

uint32_t SlimCounterRef::minPos() {
	std::shared_lock<std::shared_timed_mutex> lock(mut_);
	auto keys = getVectorOfMapKeys(counts_);
	return *std::min_element(keys.begin(), keys.end());
}
uint32_t SlimCounterRef::minPosLockFree()const {
	auto keys = getVectorOfMapKeys(counts_);
	return *std::min_element(keys.begin(), keys.end());
}



uint32_t SlimCounterRef::numOfCountsLockFree() const{
	return counts_.size();
}

uint32_t SlimCounterRef::numOfCounts() {
	std::shared_lock<std::shared_timed_mutex> lock(mut_);
	return counts_.size();
}

bool SlimCounterRef::flushCountsToFile(uint32_t fromPositionIncluding, uint32_t upToPosNotIncluding){

	if(upToPosNotIncluding < fromPositionIncluding ){
		std::stringstream ss;
		ss << __FILE__ << " : " << __LINE__ << " : " << __PRETTY_FUNCTION__ << ", error, end position shouldn't be less than starting position" << "\n";
		ss << "fromPositionIncluding: " << fromPositionIncluding << ", upToPosNotIncluding: " << upToPosNotIncluding << "\n";
		throw std::runtime_error{ss.str()};
	}

	std::unique_lock<std::shared_timed_mutex> lock(mut_);
	if(0 == numOfCountsLockFree()){
		return false;
	}
	//if attempting to flush counting in a range when there are still counts less than the input range then return nothing

	if(fromPositionIncluding > minPosLockFree() ){
		return false;
	}

	//open the file for writing if it isn't already
	if(nullptr == bgzOutFile_){
		bgzOutFile_ = std::make_unique<BGZFCPP>(outputFile_.string().c_str(), "w");
	}

	auto keys = getVectorOfMapKeys(counts_);
	njh::sort(keys);
	std::vector<uint32_t> keysWritten;
	bool wrote = false;
	for (const auto key : keys) {
		if(key > upToPosNotIncluding){
			break;
		}
		auto outVec = counts_[key].genOutVec();
		outVec[0] = refId_;
		outVec[1] = key;
		bgzOutFile_->writeVec(outVec);
		keysWritten.push_back(key);
		wrote = true;
	}
	for(const auto key : keysWritten){
		counts_.erase(key);
	}
	return wrote;
}

void SlimCounterRef::indexFile(){
	std::unique_lock<std::shared_timed_mutex> lock(mut_);
	if(nullptr != bgzOutFile_){
		/**@todo closing file by setting bgzOutFile_ to nullptr which should call the destructor
		 * of bgzOutFile_ which should close the bgzf file,
		 * need to make sure this will have the file closed before
		 *  reoping of reading for the creating of the index */
		bgzOutFile_ = nullptr;
	}

	bgzOutFile_ = std::make_unique<BGZFCPP>(outputFile_, "r");
	uint32_t count = 0;
	auto outFile = njh::appendAsNeededRet(outputFile_.string(), ".idx.gz");
	if (bfs::exists(outFile)) {
		bfs::remove(outFile);
	}
	//std::cout << "IN:  " << inFile  << std::endl;
	//std::cout << "OUT: " << outFile << std::endl;
	auto gzFileOut = gzopen(outFile.c_str(), "w");
	std::vector<uint32_t> d(SlimCounterPos::NumOfElements);
	while (true) {
		//get file loc
		auto fileLoc = bgzOutFile_->tell();
		if(!bgzOutFile_->readVec(d)){
			break;
		}
		SlimCounterPos pos(d);
		long long int refPos = d[1];
		long long int cov = pos.getHqBaseTotal();
		std::vector<long long int> outVec { refPos, cov, fileLoc };
		gzwrite(gzFileOut, outVec.data(), sizeof(long long int) * IndexNumOfElements);
		++count;
	}
	gzclose(gzFileOut);
}



}  // namespace njhseq

