/*
 * BamCountingUtils.cpp
 *
 *  Created on: Mar 2, 2016
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

#include "BamCountingUtils.hpp"
#include <njhseq/BamToolsUtils/BamToolsUtils.hpp>
#include <njhseq/seqToolsUtils/seqToolsUtils.hpp>


namespace njhseq {




bool determinQualityIonTorrent(aligner & alignerObj, uint32_t hRunSize,
		const seqInfo & query, uint64_t pos, uint32_t qualCutOff, bool debug) {
	auto seqPos = getRealPosForAlnPos(alignerObj.alignObjectB_.seqBase_.seq_,
			pos);
	//determine if seq base is the last or first base in a homopolymer of 3 bases or larger
	//first base
	bool firstBase = false;
	bool passHRunSizeTest = false;
	if (seqPos < query.seq_.size() - (hRunSize - 1) && seqPos > 0) {
		if (query.seq_[seqPos - 1] != query.seq_[seqPos]) {
			firstBase = true;
			bool passSize = true;
			for (auto hPos : iter::range<size_t>(1, hRunSize)) {
				if (query.seq_[seqPos + hPos] != query.seq_[seqPos]) {
					passSize = false;
					break;
				}
			}
			passHRunSizeTest = passSize;
		}
	}
	//if it wasn't the first base in a run, check if it is the last base
	bool lastBase = false;
	if (!firstBase) {
		if (seqPos > (hRunSize - 2) && seqPos + 1 != query.seq_.size()) {
			if (query.seq_[seqPos + 1] != query.seq_[seqPos]) {
				lastBase = true;
				bool passSize = true;
				for (auto hPos : iter::range<size_t>(1, hRunSize)) {
					if (query.seq_[seqPos - hPos] != query.seq_[seqPos]) {
						passSize = false;
						break;
					}
				}
				passHRunSizeTest = passSize;
			}
		}
	}
	if (passHRunSizeTest) {
		size_t otherEndPos = std::numeric_limits<size_t>::max();
		if (lastBase) {
			otherEndPos = seqPos - 2;
			while (otherEndPos > 0
					&& query.seq_[otherEndPos - 1] == query.seq_[seqPos]) {
				--otherEndPos;
			}
		} else {
			otherEndPos = seqPos + 2;
			while (otherEndPos + 1 < len(query)
					&& query.seq_[otherEndPos + 1] == query.seq_[seqPos]) {
				++otherEndPos;
			}
		}
		if (debug) {
			if (lastBase) {
				std::cout << "lastBase:" << seqPos << ":" << otherEndPos << std::endl;
				for (auto pos : iter::range(otherEndPos, seqPos + 1)) {
					std::cout << pos << ":" << query.seq_[pos] << ",";
				}
				std::cout << std::endl;
			} else {
				std::cout << "firstBase:" << seqPos << ":" << otherEndPos << std::endl;
				for (auto pos : iter::range(seqPos, otherEndPos + 1)) {
					std::cout << pos << ":" << query.seq_[pos] << ",";
				}
				std::cout << std::endl;
			}
		}
		return query.qual_[seqPos] > qualCutOff
				&& query.qual_[otherEndPos] > qualCutOff;
	}else{
		return alignerObj.alignObjectB_.seqBase_.qual_[pos] > qualCutOff;
	}
}

void increaseCountsForAln(const BamTools::BamAlignment & bAln,
		const std::string & rName, RefCounter & counter, aligner & alignerObj,
		uint32_t qualCutOff, bool ionTorrent, bool reAlignTorrent, uint32_t hRunSize,
		uint64_t & highQualityBases, uint64_t & lowQualityBases, bool debug, uint32_t readPosCutOff,
		size_t start, size_t end) {

	//convert to alnInfo
	if(debug){
		std::cout <<__PRETTY_FUNCTION__ << ": starting" << std::endl;
		std::cout <<"[" << std::this_thread::get_id() << "]: Starting Counting" << std::endl;
		std::cout <<"[" << std::this_thread::get_id() << "]: Converting to aln info" << std::endl;
	}
	auto alnInfos = bamAlnToAlnInfoLocal(bAln);
	if(debug){
		std::cout <<"[" << std::this_thread::get_id() << "]: Converting bam aln to seqInfo" << std::endl;
	}
	auto query = seqInfo(bAln.Name, bAln.QueryBases, bAln.Qualities,
			SangerQualOffset);
	if(debug){
		std::cout <<"[" << std::this_thread::get_id() << "]: Starting on alnInfos" << std::endl;
	}
	for(const auto & info : alnInfos){
		//get the sequence for the reference at this alignment location
		std::string seq = "";
		auto ref = seqInfo("ref",
				counter.reader_[rName]->getSequence(seq, info.first,
						info.first + info.second.localASize_));
		//convert bam alignment to njhseq::seqInfo object
		alignerObj.alignObjectA_ = baseReadObject(ref);
		alignerObj.alignObjectB_ = baseReadObject(query);
		//set the alignment info and and rearrange the sequences so they can be profiled with gaps
		alignerObj.parts_.lHolder_ = info.second;
		alignerObj.rearrangeObjsLocal(alignerObj.alignObjectA_,
				alignerObj.alignObjectB_);
		//a quick check to see if the rearrangement went well
		/*if(alignerObj.alignObjectB_.seqBase_.seq_ != bAln.AlignedBases){
			alignerObj.alignObjectA_.seqBase_.outPutSeqAnsi(std::cerr);
			auto queryAln = seqInfo(bAln.Name, bAln.AlignedBases);
			queryAln.outPutSeqAnsi(std::cerr);
			exit(1);
		}*/
		uint64_t currentMaxLen = alignerObj.parts_.maxSize_ - 1;
		if(ionTorrent && reAlignTorrent){
			readVec::getMaxLength(ref, currentMaxLen);
			readVec::getMaxLength(query, currentMaxLen);
			if(currentMaxLen > alignerObj.parts_.maxSize_){
				alignerObj.parts_.setMaxSize(currentMaxLen);
			}
			alignerObj.alignCacheGlobal(ref, query);
		}

		uint32_t aOffSet = 0;
		uint32_t bOffSet = 0;
		if(debug){
			query.outPutSeqAnsi(std::cout);
			auto seqTab = getSeqPosTab(query.seq_);
			seqTab.outPutContentOrganized(std::cout);
		}
		alignerObj.weighHomopolymers_ = false;
		alignerObj.profilePrimerAlignment(ref, query);
		if(debug){
			std::cout << "Now counting base information" << std::endl;
		}
		for (auto pos : iter::range(len(alignerObj.alignObjectB_))) {
			if ('-' == alignerObj.alignObjectA_.seqBase_.seq_[pos]) {
				++aOffSet;
			} else if ('-' == alignerObj.alignObjectB_.seqBase_.seq_[pos]) {
				++bOffSet;
			} else {
				if(info.first - aOffSet + pos >= start
						&& info.first - aOffSet + pos < end
						&& pos - bOffSet >= readPosCutOff
						&& pos - bOffSet < len(query) - readPosCutOff){
					bool highQuality = alignerObj.alignObjectB_.seqBase_.qual_[pos] > qualCutOff;
					if(ionTorrent){
						highQuality = determinQualityIonTorrent(alignerObj, hRunSize, query, pos, qualCutOff, debug);
					}
					if (highQuality) {
						++highQualityBases;
					} else {
						++lowQualityBases;
					}
					counter.increaseBaseCount(rName, info.first - aOffSet + pos,
							!bAln.IsReverseStrand(), highQuality,
							alignerObj.alignObjectB_.seqBase_.seq_[pos], 1);
				}
			}
		}
		if(debug){
			std::cout << "Now counting gap information" << std::endl;
		}
		for (const auto & alnGap : alignerObj.comp_.distances_.alignmentGaps_) {
			if (info.first + alnGap.second.refPos_ >= start
					&& info.first + alnGap.second.refPos_ < end
					&& alnGap.second.seqPos_ > readPosCutOff
					&& alnGap.second.seqPos_ < len(query) - readPosCutOff) {
				counter.increaseIndelCount(rName, info.first + alnGap.second.refPos_,
						!bAln.IsReverseStrand(), alnGap.second.ref_,
						alnGap.second.gapedSequence_, 1);
			}
		}
		if(debug){
			alignerObj.alignObjectA_.seqBase_.outPutSeqAnsi(std::cout);
			alignerObj.alignObjectB_.seqBase_.outPutSeqAnsi(std::cout);
		}
	}
}

void extractRef(
		concurrent::BamReaderPool & bamPool, concurrent::AlignerPool & alnPool,
		njh::concurrent::LockableQueue<std::string>& chroms,
		MultiRefCounter& masterCounter,std::unordered_map<std::string, BamCountExtractStats>& refExtracts,
		std::mutex & refExtractsMut,
		 uint32_t mappingQual, size_t minLen, uint32_t qualCutOff,
		bool ionTorrent, bool reAlignTorrent, uint32_t hRunSize,
		bool verbose, uint32_t readPosCutOff){
	if(verbose){
		std::cout <<"[" << std::this_thread::get_id() << "]: Starting Thread" << std::endl;
	}
	std::string refName = "";
	if(verbose){
		std::cout <<"[" << std::this_thread::get_id() << "]: Getting Reader" << std::endl;
	}
	auto reader = bamPool.popReader();
	if(verbose){
		std::cout <<"[" << std::this_thread::get_id() << "]: Getting Aligner" << std::endl;
	}
	auto alignerObj = alnPool.popAligner();
	if(verbose){
		std::cout <<"[" << std::this_thread::get_id() << "]: Starting extraction" << std::endl;
	}
	while(chroms.getVal(refName)){

		if(verbose){
			std::cout <<"[" << std::this_thread::get_id() << "]: Extracting on " << refName << std::endl;
		}
		BamCountExtractStats refStats;
		RefCounter counter(masterCounter.masterCounter_.reader_.getFilename().string());
		counter.setHardCutOff(masterCounter.masterCounter_.hardCutOff_);
		if(verbose){
			std::cout <<"[" << std::this_thread::get_id() << "]: Getting ref data " << std::endl;
		}
		auto refId = reader->GetReferenceID(refName);
		auto refData = reader->GetReferenceData();
		BamTools::BamAlignment bAln;
		if(verbose){
			std::cout <<"[" << std::this_thread::get_id() << "]: Setting region" << std::endl;
		}
		if(!reader->SetRegion(refId, 0, refId, refData[refId].RefLength)){
			std::stringstream ss;
			ss << reader->GetErrorString() << std::endl;
			ss << "Failed to get region in " << std::this_thread::get_id() << std::endl;
			ss << "Region: " << refId << std::endl;
			throw std::runtime_error{ss.str()	};
		}

		if(verbose){
			std::cout <<"[" << std::this_thread::get_id() << "]: Reading alignments" << std::endl;
		}
		while (reader->GetNextAlignment(bAln)) {
			++refStats.totalReads;
			if(verbose){
				std::cout <<"[" << std::this_thread::get_id() << "]: On " << refStats.totalReads << std::endl;
			}
			if (bAln.IsMapped() && bAln.MapQuality > mappingQual) {
				++refStats.readsMapped;
				//convert cigar data into alnInfo
				auto alnInfos = bamAlnToAlnInfoLocal(bAln);
				uint32_t totalSize = 0;
				for(const auto & info : alnInfos){
					totalSize += info.second.localBSize_;
				}
				if(totalSize < minLen){
					++refStats.readsBellowMinLen;
					continue;
				}
				++refStats.readsUsed;
				if(verbose){
					std::cout <<"[" << std::this_thread::get_id() << "]: Increasing counts for " << refStats.totalReads << std::endl;
				}

				increaseCountsForAln(bAln, refName, counter, *alignerObj, qualCutOff,
						ionTorrent, reAlignTorrent, hRunSize, refStats.highQualityBases,
						refStats.lowQaulityBases, verbose, readPosCutOff);
				if(verbose){
					std::cout <<"[" << std::this_thread::get_id() << "]: Finished increasing counts for " << refStats.totalReads << std::endl;
				}
			} else {
				if (!bAln.IsMapped()) {
					++refStats.readsNotMapped;
				} else if (bAln.MapQuality <= mappingQual) {
					++refStats.readsBellowMappingQuality;
				}
			}
		}

		{
			if(verbose){
				std::cout <<"[" << std::this_thread::get_id() << "]: Adding extraction stats for " << refName << std::endl;
			}
			std::lock_guard<std::mutex> statsLock(refExtractsMut);
			refExtracts[refName] = refStats;
		}
		if(verbose){
			std::cout <<"[" << std::this_thread::get_id() << "]: Adding to counter for " << refName << std::endl;
		}
		masterCounter.addCounter(refName, counter);
		if(verbose){
			std::cout <<"[" << std::this_thread::get_id() << "]: Finished Extracting for " << refName << std::endl;
		}
	}
}

}  // namespace njhse
