/*
 * bamCountingExpRunner.cpp
 *
 *  Created on: Jan 25, 2015
 *      Author: nickhathaway
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

#include "bamCountingExp.hpp"

#include <TwoBit.h>

#include "elucidator/BamToolsUtils.h"
#include "elucidator/objects/SlimCounter.h"
#include "elucidator/filesystem/GzSimpleBinFile.hpp"
#include "elucidator/objects/Gatherers.h"

namespace njhseq {
bamCountingExpRunner::bamCountingExpRunner()
    : njh::progutils::ProgramRunner(
          {
					 addFunc("bamBaseCountingSlim", bamBaseCountingSlim, false),
					 addFunc("slimCounterBinToTxt", slimCounterBinToTxt, false),
					 addFunc("readGzIndex", readGzIndex, false),
					 addFunc("coverageInfo", coverageInfo, false),
					 addFunc("gatherVariantsOnChrom", gatherVariantsOnChrom, false),
					 addFunc("gatherVariantsOnSample", gatherVariantsOnSample, false),
					 addFunc("gatherLocOnSample", gatherLocOnSample, false),
					 addFunc("gatherVariantsOnSampleOnLoc", gatherVariantsOnSampleOnLoc, false),
					 addFunc("gatherBaseCoverageOnSampleOnLoc", gatherBaseCoverageOnSampleOnLoc, false),
           },//,
          "bamCountingExp") {}


void increaseCountsForAlnSlim(const BamTools::BamAlignment & bAln,
		SlimCounterRef & counter, aligner & alignerObj,
		uint32_t qualCutOff, bool ionTorrent, bool reAlignTorrent, uint32_t hRunSize,
		uint64_t & highQualityBases, uint64_t & lowQualityBases, bool debug, uint32_t readPosCutOff,
		size_t start, size_t end) {

	//convert to alnInfo
	if(debug){
		std::cout <<__PRETTY_FUNCTION__ << ": starting" << std::endl;
		std::cout <<"[TID=" << std::this_thread::get_id() << "]: Starting Counting" << std::endl;
		std::cout <<"[TID=" << std::this_thread::get_id() << "]: Converting to aln info" << std::endl;
	}
	auto alnInfos = bamAlnToAlnInfoLocal(bAln);
	if(debug){
		std::cout <<"[TID=" << std::this_thread::get_id() << "]: Converting bam aln to seqInfo" << std::endl;
	}
	auto query = seqInfo(bAln.Name, bAln.QueryBases, bAln.Qualities,
			SangerQualOffset);
	if(debug){
		std::cout <<"[TID=" << std::this_thread::get_id() << "]: Starting on alnInfos" << std::endl;
	}

	std::shared_lock<std::shared_timed_mutex> lock(counter.mut_);
	for(const auto & info : alnInfos){
		//get the sequence for the reference at this alignment location
		std::string seq = std::string(info.second.localASize_, 'X');
		auto ref = seqInfo("ref", seq);
		//convert bam alignment to njhseq::seqInfo object
		alignerObj.alignObjectA_ = baseReadObject(ref);
		alignerObj.alignObjectB_ = baseReadObject(query);
		//set the alignment info and and rearrange the sequences so they can be profiled with gaps
		alignerObj.parts_.lHolder_ = info.second;
		alignerObj.rearrangeObjsLocal(alignerObj.alignObjectA_,
				alignerObj.alignObjectB_);
		uint64_t currentMaxLen = alignerObj.parts_.maxSize_;
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
					bool highQuality = alignerObj.alignObjectB_.seqBase_.qual_[pos]
							> qualCutOff;
					if (ionTorrent) {
						highQuality = determinQualityIonTorrent(alignerObj, hRunSize, query,
								pos, qualCutOff, debug);
					}
					if (highQuality) {
						++highQualityBases;
					} else {
						++lowQualityBases;
					}
					counter.increaseCountLockFree(info.first - aOffSet + pos,
							alignerObj.alignObjectB_.seqBase_.seq_[pos],
							bAln.IsReverseStrand(), highQuality, 1);
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
				if(alnGap.second.ref_){
					counter.increaseInsertionCountLockFree(info.first + alnGap.second.refPos_,bAln.IsReverseStrand(),1);
				}else{
					counter.increaseDeletionCountLockFree(info.first + alnGap.second.refPos_,bAln.IsReverseStrand(),1);
				}
			}
		}
		if(debug){
			alignerObj.alignObjectA_.seqBase_.outPutSeqAnsi(std::cout);
			alignerObj.alignObjectB_.seqBase_.outPutSeqAnsi(std::cout);
		}
	}
}

struct PosInfoForCount{
	enum class Type{
		INS, DEL, BASE
	};
	size_t pos_;
	char base_;
	uint32_t qual_;
	bool highQuality_;
	bool isReverseStrand_;
	Type type_;
};

void increaseCountsForAlnSlimForOverlapPairs(const BamTools::BamAlignment & firstPosMate,
		const BamTools::BamAlignment & secondPosMate,
		SlimCounterRef & counter, aligner & alignerObj,
		uint32_t qualCutOff,
		uint64_t & highQualityBases, uint64_t & lowQualityBases, bool debug, uint32_t readPosCutOff,
		size_t start, size_t end) {

	//convert to alnInfo
	if(debug){
		std::cout <<__PRETTY_FUNCTION__ << ": starting" << std::endl;
		std::cout <<"[TID=" << std::this_thread::get_id() << "]: Starting Counting" << std::endl;
		std::cout <<"[TID=" << std::this_thread::get_id() << "]: Converting to aln info" << std::endl;
	}
	auto firstAlnInfos = bamAlnToAlnInfoLocal(firstPosMate);
	auto secondAlnInfos = bamAlnToAlnInfoLocal(secondPosMate);
	if(debug){
		std::cout <<"[TID=" << std::this_thread::get_id() << "]: Converting bam aln to seqInfo" << std::endl;
	}

	auto firstQuery = seqInfo(firstPosMate.Name, firstPosMate.QueryBases, firstPosMate.Qualities,
			SangerQualOffset);
	auto secondQuery =  seqInfo(secondPosMate.Name, secondPosMate.QueryBases, secondPosMate.Qualities,
			SangerQualOffset);
	if(debug){
		std::cout <<"[TID=" << std::this_thread::get_id() << "]: Starting on alnInfos" << std::endl;
	}

	std::unordered_map<size_t, PosInfoForCount> firstOverLapBaseCounts;


	std::unordered_map<size_t, PosInfoForCount> secondOverLapBaseCounts;


	size_t minOverlapPos = std::numeric_limits<size_t>::max();
	size_t maxOverlapPos = 0;
	std::shared_lock<std::shared_timed_mutex> lock(counter.mut_);
	for(const auto & info : firstAlnInfos){
		//get the sequence for the reference at this alignment location
		std::string seq = std::string(info.second.localASize_, 'X');
		auto ref = seqInfo("ref", seq);
		//convert bam alignment to njhseq::seqInfo object
		alignerObj.alignObjectA_ = baseReadObject(ref);
		alignerObj.alignObjectB_ = baseReadObject(firstQuery);
		//set the alignment info and and rearrange the sequences so they can be profiled with gaps
		alignerObj.parts_.lHolder_ = info.second;
		alignerObj.rearrangeObjsLocal(alignerObj.alignObjectA_,
				alignerObj.alignObjectB_);


		uint32_t aOffSet = 0;
		uint32_t bOffSet = 0;
		if(debug){
			firstQuery.outPutSeqAnsi(std::cout);
			auto seqTab = getSeqPosTab(firstQuery.seq_);
			seqTab.outPutContentOrganized(std::cout);
		}
		alignerObj.weighHomopolymers_ = false;
		alignerObj.profilePrimerAlignment(ref, firstQuery);
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
						&& pos - bOffSet < len(firstQuery) - readPosCutOff){
					bool highQuality = alignerObj.alignObjectB_.seqBase_.qual_[pos]
							> qualCutOff;
					if (highQuality) {
						++highQualityBases;
					} else {
						++lowQualityBases;
					}
					if(static_cast<int64_t>(info.first - aOffSet + pos) >= secondPosMate.Position){
						if(info.first - aOffSet + pos < minOverlapPos){
							minOverlapPos = info.first - aOffSet + pos;
						}
						if (info.first - aOffSet + pos > maxOverlapPos) {
							maxOverlapPos = info.first - aOffSet + pos;
						}
						firstOverLapBaseCounts[info.first - aOffSet + pos] = {info.first - aOffSet + pos,
							alignerObj.alignObjectB_.seqBase_.seq_[pos],
							alignerObj.alignObjectB_.seqBase_.qual_[pos],
							highQuality,
							firstPosMate.IsReverseStrand(),
							PosInfoForCount::Type::BASE};
					}else{
						counter.increaseCountLockFree(info.first - aOffSet + pos,
													alignerObj.alignObjectB_.seqBase_.seq_[pos],
													firstPosMate.IsReverseStrand(), highQuality, 1);
					}
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
					&& alnGap.second.seqPos_ < len(firstQuery) - readPosCutOff) {
				if(alnGap.second.ref_){
					if(static_cast<int64_t>(info.first + alnGap.second.refPos_ )>= secondPosMate.Position){
						if(info.first + alnGap.second.refPos_ < minOverlapPos){
							minOverlapPos = info.first + alnGap.second.refPos_;
						}
						if (info.first + alnGap.second.refPos_ > maxOverlapPos) {
							maxOverlapPos = info.first + alnGap.second.refPos_;
						}
						auto meanQual = vectorMean(alnGap.second.qualities_);
						firstOverLapBaseCounts[info.first + alnGap.second.refPos_] = {info.first + alnGap.second.refPos_,
							'-',
							static_cast<uint32_t>(std::round(meanQual)),
							std::round(meanQual)>qualCutOff,
							firstPosMate.IsReverseStrand(),
							PosInfoForCount::Type::INS};
					}else{
						counter.increaseInsertionCountLockFree(info.first + alnGap.second.refPos_, firstPosMate.IsReverseStrand(),1);
					}
				}else{
					if(static_cast<int64_t>(info.first + alnGap.second.refPos_) >= secondPosMate.Position){
						if(info.first + alnGap.second.refPos_ < minOverlapPos){
							minOverlapPos = info.first + alnGap.second.refPos_;
						}
						if (info.first + alnGap.second.refPos_ > maxOverlapPos) {
							maxOverlapPos = info.first + alnGap.second.refPos_;
						}
						firstOverLapBaseCounts[info.first + alnGap.second.refPos_] = {info.first + alnGap.second.refPos_,
							'-',
							40,
							true,
							firstPosMate.IsReverseStrand(),
							PosInfoForCount::Type::DEL};
					}else{
						counter.increaseDeletionCountLockFree(info.first + alnGap.second.refPos_, firstPosMate.IsReverseStrand(),1);
					}
				}
			}
		}
		if(debug){
			alignerObj.alignObjectA_.seqBase_.outPutSeqAnsi(std::cout);
			alignerObj.alignObjectB_.seqBase_.outPutSeqAnsi(std::cout);
		}
	}

	for(const auto & info : secondAlnInfos){
		//get the sequence for the reference at this alignment location
		std::string seq = std::string(info.second.localASize_, 'X');
		auto ref = seqInfo("ref", seq);
		//convert bam alignment to njhseq::seqInfo object
		alignerObj.alignObjectA_ = baseReadObject(ref);
		alignerObj.alignObjectB_ = baseReadObject(secondQuery);
		//set the alignment info and and rearrange the sequences so they can be profiled with gaps
		alignerObj.parts_.lHolder_ = info.second;
		alignerObj.rearrangeObjsLocal(alignerObj.alignObjectA_,
				alignerObj.alignObjectB_);


		uint32_t aOffSet = 0;
		uint32_t bOffSet = 0;
		if(debug){
			secondQuery.outPutSeqAnsi(std::cout);
			auto seqTab = getSeqPosTab(secondQuery.seq_);
			seqTab.outPutContentOrganized(std::cout);
		}
		alignerObj.weighHomopolymers_ = false;
		alignerObj.profilePrimerAlignment(ref, secondQuery);
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
						&& pos - bOffSet < len(secondQuery) - readPosCutOff){
					bool highQuality = alignerObj.alignObjectB_.seqBase_.qual_[pos]
							> qualCutOff;
					if (highQuality) {
						++highQualityBases;
					} else {
						++lowQualityBases;
					}
					if(static_cast<int64_t>(info.first - aOffSet + pos )<= firstPosMate.GetEndPosition()){
						if(info.first - aOffSet + pos < minOverlapPos){
							minOverlapPos = info.first - aOffSet + pos;
						}
						if (info.first - aOffSet + pos > maxOverlapPos) {
							maxOverlapPos = info.first - aOffSet + pos;
						}
						secondOverLapBaseCounts[info.first - aOffSet + pos] = {info.first - aOffSet + pos,
							alignerObj.alignObjectB_.seqBase_.seq_[pos],
							alignerObj.alignObjectB_.seqBase_.qual_[pos],
							highQuality,
							secondPosMate.IsReverseStrand(),
							PosInfoForCount::Type::BASE};
					}else{
						counter.increaseCountLockFree(info.first - aOffSet + pos,
													alignerObj.alignObjectB_.seqBase_.seq_[pos],
													secondPosMate.IsReverseStrand(), highQuality, 1);
					}
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
					&& alnGap.second.seqPos_ < len(secondQuery) - readPosCutOff) {
				if(alnGap.second.ref_){
					if(static_cast<int64_t>(info.first + alnGap.second.refPos_) <= firstPosMate.GetEndPosition()){
						if(info.first + alnGap.second.refPos_ < minOverlapPos){
							minOverlapPos = info.first + alnGap.second.refPos_;
						}
						if (info.first + alnGap.second.refPos_ > maxOverlapPos) {
							maxOverlapPos = info.first + alnGap.second.refPos_;
						}
						auto meanQual = vectorMean(alnGap.second.qualities_);
						secondOverLapBaseCounts[info.first + alnGap.second.refPos_] = {info.first + alnGap.second.refPos_,
							'-',
							static_cast<uint32_t>(std::round(meanQual)),
							std::round(meanQual)>qualCutOff,
							secondPosMate.IsReverseStrand(),
							PosInfoForCount::Type::INS};
					}else{
						counter.increaseInsertionCountLockFree(info.first + alnGap.second.refPos_, secondPosMate.IsReverseStrand(),1);
					}
				}else{
					if(static_cast<int64_t>(info.first + alnGap.second.refPos_) <= firstPosMate.GetEndPosition()){
						if(info.first + alnGap.second.refPos_ < minOverlapPos){
							minOverlapPos = info.first + alnGap.second.refPos_;
						}
						if (info.first + alnGap.second.refPos_ > maxOverlapPos) {
							maxOverlapPos = info.first + alnGap.second.refPos_;
						}
						secondOverLapBaseCounts[info.first + alnGap.second.refPos_] = {info.first + alnGap.second.refPos_,
							'-',
							40,
							true,
							secondPosMate.IsReverseStrand(),
							PosInfoForCount::Type::DEL};

					}else{
						counter.increaseDeletionCountLockFree(info.first + alnGap.second.refPos_, secondPosMate.IsReverseStrand(),1);
					}
				}
			}
		}
		if(debug){
			alignerObj.alignObjectA_.seqBase_.outPutSeqAnsi(std::cout);
			alignerObj.alignObjectB_.seqBase_.outPutSeqAnsi(std::cout);
		}
	}
	//now handle overlap counts

	for(const auto & pos : iter::range(minOverlapPos, maxOverlapPos + 1)){
		if (njh::in(pos, firstOverLapBaseCounts)
				&& njh::in(pos, secondOverLapBaseCounts)) {
			switch (firstOverLapBaseCounts[pos].type_) {
			case PosInfoForCount::Type::BASE:
				switch (secondOverLapBaseCounts[pos].type_) {
				case PosInfoForCount::Type::BASE:
					if (firstOverLapBaseCounts[pos].qual_
							== secondOverLapBaseCounts[pos].qual_) {
						if (firstPosMate.MapQuality >= secondPosMate.MapQuality) {
							counter.increaseCountLockFree(pos, firstOverLapBaseCounts[pos].base_,
									firstOverLapBaseCounts[pos].isReverseStrand_,
									firstOverLapBaseCounts[pos].highQuality_, 1);
						} else {
							counter.increaseCountLockFree(pos, secondOverLapBaseCounts[pos].base_,
									secondOverLapBaseCounts[pos].isReverseStrand_,
									secondOverLapBaseCounts[pos].highQuality_, 1);
						}
					} else if (firstOverLapBaseCounts[pos].qual_
							>= secondOverLapBaseCounts[pos].qual_) {
						counter.increaseCountLockFree(pos, firstOverLapBaseCounts[pos].base_,
								firstOverLapBaseCounts[pos].isReverseStrand_,
								firstOverLapBaseCounts[pos].highQuality_, 1);
					} else {
						counter.increaseCountLockFree(pos, secondOverLapBaseCounts[pos].base_,
								secondOverLapBaseCounts[pos].isReverseStrand_,
								secondOverLapBaseCounts[pos].highQuality_, 1);
					}
					break;
				case PosInfoForCount::Type::INS:
					if(firstOverLapBaseCounts[pos].highQuality_){
						counter.increaseCountLockFree(pos, firstOverLapBaseCounts[pos].base_,
								firstOverLapBaseCounts[pos].isReverseStrand_,
								firstOverLapBaseCounts[pos].highQuality_, 1);
					}else{
						counter.increaseInsertionCountLockFree(pos,
								secondOverLapBaseCounts[pos].isReverseStrand_, 1);
					}
					break;
				case PosInfoForCount::Type::DEL:
					if(firstOverLapBaseCounts[pos].highQuality_){
						counter.increaseCountLockFree(pos, firstOverLapBaseCounts[pos].base_,
								firstOverLapBaseCounts[pos].isReverseStrand_,
								firstOverLapBaseCounts[pos].highQuality_, 1);
					}else{
						counter.increaseDeletionCountLockFree(pos,
								secondOverLapBaseCounts[pos].isReverseStrand_, 1);
					}
					break;
				default:
					std::stringstream ss;
					ss << "Error in " << __PRETTY_FUNCTION__ << " unhandled case\n";
					throw std::runtime_error { ss.str() };
					break;
				}
				break;
			case PosInfoForCount::Type::INS:
				switch (secondOverLapBaseCounts[pos].type_) {
				case PosInfoForCount::Type::BASE:
					if(secondOverLapBaseCounts[pos].highQuality_){
						counter.increaseCountLockFree(pos, secondOverLapBaseCounts[pos].base_,
								secondOverLapBaseCounts[pos].isReverseStrand_,
								secondOverLapBaseCounts[pos].highQuality_, 1);
					}else{
						counter.increaseInsertionCountLockFree(pos,
								firstOverLapBaseCounts[pos].isReverseStrand_, 1);
					}
					break;
				case PosInfoForCount::Type::INS:
					counter.increaseInsertionCountLockFree(pos,
							firstOverLapBaseCounts[pos].isReverseStrand_, 1);
					break;
				case PosInfoForCount::Type::DEL:
					if(firstOverLapBaseCounts[pos].highQuality_){
						counter.increaseInsertionCountLockFree(pos,
								firstOverLapBaseCounts[pos].isReverseStrand_, 1);
					}else{
						counter.increaseDeletionCountLockFree(pos,
								secondOverLapBaseCounts[pos].isReverseStrand_, 1);
					}
					break;
				default:
					std::stringstream ss;
					ss << "Error in " << __PRETTY_FUNCTION__ << " unhandled case\n";
					throw std::runtime_error { ss.str() };
					break;
				}
				break;
			case PosInfoForCount::Type::DEL:
				switch (secondOverLapBaseCounts[pos].type_) {
				case PosInfoForCount::Type::BASE:
					counter.increaseCountLockFree(pos, secondOverLapBaseCounts[pos].base_,
							secondOverLapBaseCounts[pos].isReverseStrand_,
							secondOverLapBaseCounts[pos].highQuality_, 1);
					break;
				case PosInfoForCount::Type::INS:
					if(secondOverLapBaseCounts[pos].highQuality_){
						counter.increaseInsertionCountLockFree(pos,
								secondOverLapBaseCounts[pos].isReverseStrand_, 1);
					}else{
						counter.increaseInsertionCountLockFree(pos,
								firstOverLapBaseCounts[pos].isReverseStrand_, 1);
					}
					break;
				case PosInfoForCount::Type::DEL:
					counter.increaseDeletionCountLockFree(pos,
							firstOverLapBaseCounts[pos].isReverseStrand_, 1);
					break;
				default:
					std::stringstream ss;
					ss << "Error in " << __PRETTY_FUNCTION__ << " unhandled case\n";
					throw std::runtime_error { ss.str() };
					break;
				}
				break;
			default:
				std::stringstream ss;
				ss << "Error in " << __PRETTY_FUNCTION__ << " unhandled case\n";
				throw std::runtime_error { ss.str() };
				break;
			}
		}else if(njh::in(pos, firstOverLapBaseCounts)){
			switch (firstOverLapBaseCounts[pos].type_) {
			case PosInfoForCount::Type::BASE:
				counter.increaseCountLockFree(pos, firstOverLapBaseCounts[pos].base_,
						firstOverLapBaseCounts[pos].isReverseStrand_,
						firstOverLapBaseCounts[pos].highQuality_, 1);
				break;
			case PosInfoForCount::Type::INS:
				counter.increaseInsertionCountLockFree(pos,
						firstOverLapBaseCounts[pos].isReverseStrand_, 1);
				break;
			case PosInfoForCount::Type::DEL:
				counter.increaseDeletionCountLockFree(pos,
						firstOverLapBaseCounts[pos].isReverseStrand_, 1);
				break;
			default:
				std::stringstream ss;
				ss << "Error in " << __PRETTY_FUNCTION__ << " unhandled case\n";
				throw std::runtime_error { ss.str() };
				break;
			}
		}else if(njh::in(pos, secondOverLapBaseCounts)){
			switch (secondOverLapBaseCounts[pos].type_) {
			case PosInfoForCount::Type::BASE:
				counter.increaseCountLockFree(pos, secondOverLapBaseCounts[pos].base_,
						secondOverLapBaseCounts[pos].isReverseStrand_,
						secondOverLapBaseCounts[pos].highQuality_, 1);
				break;
			case PosInfoForCount::Type::INS:
				counter.increaseInsertionCountLockFree(pos,
						secondOverLapBaseCounts[pos].isReverseStrand_, 1);
				break;
			case PosInfoForCount::Type::DEL:
				counter.increaseDeletionCountLockFree(pos,
						secondOverLapBaseCounts[pos].isReverseStrand_, 1);
				break;
			default:
				std::stringstream ss;
				ss << "Error in " << __PRETTY_FUNCTION__ << " unhandled case\n";
				throw std::runtime_error { ss.str() };
				break;
			}
		}
	}
}

//
//
//hRunSize
//setUp.pars_.debug_, readPosCutOff


void BamBaseCountingSlim(BamTools::BamReader & bReader, aligner & alignerObj,
		uint32_t mappingQual, uint32_t minLen, uint32_t insertLengthCutOff,
		uint32_t qualCutOff, bool ionTorrent, bool reAlignTorrent,
		uint32_t hRunSize, uint32_t readPosCutOff,
		SlimCounterMaster & masterCounter,
		std::shared_ptr<SlimCounterRef> & currentRef,
		BamCountExtractStats & extractStats, BamAlnsCache & alnCache,
		uint32_t testNumber, bool countAlones, bool debug, bool verbose,
		bool noOverlapHandling) {

	BamTools::BamAlignment bAln;
	while (bReader.GetNextAlignment(bAln)) {
		if(!bAln.IsPrimaryAlignment()){
			//if not primary alignment continue on
			continue;
		}
		if(currentRef == nullptr || static_cast<int32_t>(currentRef->refId_) != bAln.RefID){
			if(nullptr != currentRef){
				if(len(alnCache) > 0){
					auto currentId = bAln.RefID;
					extractStats.failedToFindMate += len(alnCache);
					auto cachedNames = alnCache.getNames();
					if(countAlones){
						for(const auto & name : cachedNames){
							auto aloneAln = alnCache.get(name);
							if(currentRef == nullptr || static_cast<int32_t>(currentRef->refId_) != aloneAln->RefID){
								currentRef = masterCounter.getSlimCounterRef(aloneAln->RefID);
							}
							++extractStats.readsUsed;
							increaseCountsForAlnSlim(*aloneAln, *currentRef, alignerObj, qualCutOff,
												ionTorrent, reAlignTorrent, hRunSize, extractStats.highQualityBases,
												extractStats.lowQaulityBases, debug, readPosCutOff,
												aloneAln->Position, aloneAln->GetEndPosition());
							alnCache.remove(name);
						}
					}else{
						for(const auto & name : cachedNames){
							alnCache.remove(name);
						}
					}
					currentRef = masterCounter.getSlimCounterRef(currentId);
					if(len(alnCache) > 0){
						std::cerr << "Cache size " << len(alnCache) << std::endl;
					}
				}
				currentRef->flushCountsToFile(0, std::numeric_limits<uint32_t>::max());
			}
			currentRef = masterCounter.getSlimCounterRef(bAln.RefID);
		}
		++extractStats.totalReads;
		//flushing reads in chunks of 100,000
		if(0 == extractStats.totalReads % 100000){
			auto mins = alnCache.getMinPos();
			for(const auto & min : mins){
				if(min.first == static_cast<int32_t>(currentRef->refId_)){
					currentRef->flushCountsToFile(0, min.second);
					break;
				}
			}
		}

		if(verbose && 0 == extractStats.totalReads % 10000){
			std::cout << "\r" << extractStats.totalReads;
			std::cout.flush();
		}
		if (bAln.IsPaired()) {
			++extractStats.pairedReadsTotal;
			if (bAln.IsMapped() && bAln.IsMateMapped()) {
				++extractStats.readsMapped;
			} else {
				//one or the other is not mapped
				++extractStats.readsNotMapped;
				continue;
			}

			if (bAln.RefID != bAln.MateRefID) {
				//mapped to different chromosome
				++extractStats.disCordantMapping;
				continue;
			}

			if (std::abs(bAln.InsertSize) > insertLengthCutOff) {
				//insert length is too long
				++extractStats.largeInsertSize;
				continue;
			}
			if (bAln.MatePosition == bAln.Position) {
				if (!alnCache.has(bAln.Name)) {
					//enter into cache for until mate is encountered
					alnCache.add(bAln);
					continue;
				}
			}
			if (bAln.MatePosition <= bAln.Position) {
				//need to check to see if there is an overlap
				if (!alnCache.has(bAln.Name)) {
					//can't find mate
					++extractStats.failedToFindMate;
					if (countAlones) {
						//single read count, baln
						increaseCountsForAlnSlim(bAln, *currentRef, alignerObj, qualCutOff,
								ionTorrent, reAlignTorrent, hRunSize,
								extractStats.highQualityBases, extractStats.lowQaulityBases,
								debug, readPosCutOff, bAln.Position, bAln.GetEndPosition());
					}
					continue;
					//std::stringstream ss;
					//ss << "Couldn't find mate for " << bAln.Name << std::endl;
					//throw std::runtime_error{ss.str()};
					//std::cerr << ss.str() << std::endl;
				} else {
					auto search = alnCache.get(bAln.Name);
					if (bAln.MapQuality <= mappingQual
							|| search->MapQuality <= mappingQual) {
						//bad mapping quality
						extractStats.readsBellowMappingQuality+=2;
						//remove from cache
						alnCache.remove(search->Name);
						continue;
					}
					if (bAln.GetEndPosition() - bAln.Position <= static_cast<int64_t>(minLen)
							|| search->GetEndPosition() - search->Position <= static_cast<int64_t>(minLen)) {
						//bad length
						extractStats.readsBellowMinLen+=2;
						//remove from cache
						alnCache.remove(search->Name);
						continue;
					}
					if(search->GetEndPosition() > bAln.Position){
						extractStats.readsUsed+=2;
						extractStats.pairedReadsUsed+=2;
						extractStats.overLappingPairedReads+=2;
						//overlap, overlap count together
						if(noOverlapHandling){
							//single read count, baln
								increaseCountsForAlnSlim(bAln, *currentRef, alignerObj, qualCutOff,
										ionTorrent, reAlignTorrent, hRunSize, extractStats.highQualityBases,
										extractStats.lowQaulityBases, debug, readPosCutOff,
										bAln.Position, bAln.GetEndPosition());
								//single read count, search
								increaseCountsForAlnSlim(*search, *currentRef, alignerObj, qualCutOff,
										ionTorrent, reAlignTorrent, hRunSize, extractStats.highQualityBases,
										extractStats.lowQaulityBases, debug, readPosCutOff,
										search->Position, search->GetEndPosition());
						}else{
							increaseCountsForAlnSlimForOverlapPairs(*search, bAln, *currentRef, alignerObj, qualCutOff, extractStats.highQualityBases,
									extractStats.lowQaulityBases, debug, readPosCutOff, search->Position, bAln.GetEndPosition());
						}
					}else{
						//no overlap, single read count both
						extractStats.readsUsed+=2;
						extractStats.pairedReadsUsed+=2;
						//single read count, baln
						increaseCountsForAlnSlim(bAln, *currentRef, alignerObj, qualCutOff,
								ionTorrent, reAlignTorrent, hRunSize, extractStats.highQualityBases,
								extractStats.lowQaulityBases, debug, readPosCutOff,
								bAln.Position, bAln.GetEndPosition());
						//single read count, search
						increaseCountsForAlnSlim(*search, *currentRef, alignerObj, qualCutOff,
								ionTorrent, reAlignTorrent, hRunSize, extractStats.highQualityBases,
								extractStats.lowQaulityBases, debug, readPosCutOff,
								search->Position, search->GetEndPosition());
					}
					//remove from cache
					alnCache.remove(search->Name);
				}
			} else {
				//enter into cache for until mate is encountered
				alnCache.add(bAln);
			}
		} else {
			if (bAln.IsMapped()) {
				++extractStats.readsMapped;
			} else {
				//is not mapped
				++extractStats.readsNotMapped;
				continue;
			}
			if (bAln.MapQuality <= mappingQual) {
				//mapping quality fails
				++extractStats.readsBellowMappingQuality;
				continue;
			}
			if (bAln.GetEndPosition() - bAln.Position <= static_cast<int64_t>(minLen)) {
				//bad length
				++extractStats.readsBellowMinLen;
				continue;
			}
			++extractStats.readsUsed;
			//single end count

			increaseCountsForAlnSlim(bAln, *currentRef, alignerObj, qualCutOff,
					ionTorrent, reAlignTorrent, hRunSize, extractStats.highQualityBases,
					extractStats.lowQaulityBases, debug, readPosCutOff,
					bAln.Position, bAln.GetEndPosition());
		}
		if(extractStats.totalReads > testNumber){
			break;
		}
	}

	if(verbose){
		std::cout << std::endl;
	}

	if(len(alnCache) > 0){
		extractStats.failedToFindMate += len(alnCache);
		auto cachedNames = alnCache.getNames();
		if(countAlones){
			for(const auto & name : cachedNames){
				auto aloneAln = alnCache.get(name);
				if(currentRef == nullptr || static_cast<int32_t>(currentRef->refId_) != aloneAln->RefID){
					currentRef = masterCounter.getSlimCounterRef(aloneAln->RefID);
				}
				++extractStats.readsUsed;
				increaseCountsForAlnSlim(*aloneAln, *currentRef, alignerObj, qualCutOff,
									ionTorrent, reAlignTorrent, hRunSize, extractStats.highQualityBases,
									extractStats.lowQaulityBases, debug, readPosCutOff,
									aloneAln->Position, aloneAln->GetEndPosition());
				alnCache.remove(name);
			}
		}else{
			for(const auto & name : cachedNames){
				alnCache.remove(name);
			}
		}
		if(len(alnCache) > 0){
			std::cerr << "Cache size " << len(alnCache) << std::endl;
		}
	}
}


int bamCountingExpRunner::slimCounterBinToTxt(const njh::progutils::CmdArgs & inputCommands){
	seqSetUp setUp(inputCommands);
	uint32_t testNum = 10;
	std::string filename = "";
	std::string twoBitFilename = "";
	std::string refIdLookupFile = "";
	setUp.setOption(filename, "--file", "Filename to index", true);
	setUp.setOption(testNum, "--testNum", "testNum");
	setUp.setOption(refIdLookupFile, "--refIdLookupFile", "Ref Id Lookup File", true);
	setUp.setOption(twoBitFilename, "--twoBitFilename", "TwoBit Filename", true);
	setUp.finishSetUp(std::cout);
	SlimCounter::genColIndex(false);
	uint32_t numBytes = SlimCounterPos::NumOfElements * sizeof(uint32_t);
	uint32_t count = 0;
	TwoBit::TwoBitFile twoBitfile(twoBitFilename);
	std::unordered_map<std::string, std::string> gSeqs;
	auto refInfos = tabToRefDataVec(table(TableIOOpts::genTabFileIn(refIdLookupFile, true)));
	std::cout << vectorToString(SlimCounterPos::fullDetailHeader(), "\t") << std::endl;
	auto gzInFile = gzopen(filename.c_str(), "r");
#if ZLIB_VERNUM >= 0x1280
	gzbuffer(gzInFile, 128 * 1024);
#endif
	std::string genomicSeq = "";
	uint32_t currentId = std::numeric_limits<uint32_t>::max();

	while (!gzeof(gzInFile)) {
		std::vector<uint32_t> d(SlimCounterPos::NumOfElements);
		int bytes_read = gzread(gzInFile, d.data(), numBytes);
		if (0 == bytes_read && gzeof(gzInFile)) {
			break;
		}
		if (bytes_read != static_cast<int64_t>(numBytes)) {
			std::stringstream ss;
			ss << "Error in reading: " << filename << std::endl;
			ss << "Read in " << bytes_read << " when expected " << numBytes
					<< std::endl;
			throw std::runtime_error { ss.str() };
		}
		SlimCounterPos pos(d);
		uint32_t refPos = d[1];
		auto rName = refInfos[d[0]].RefName;
		if(currentId != d[0]){
			twoBitfile[rName]->getSequence(genomicSeq);
			currentId = d[0];
		}

		pos.fullDetailOut(rName, refPos, genomicSeq[refPos], std::cout);
		++count;
		if(count > testNum){
			break;
		}
	}
	gzclose(gzInFile);

	return 0;
}

int bamCountingExpRunner::readGzIndex(const njh::progutils::CmdArgs & inputCommands){
	seqSetUp setUp(inputCommands);
	std::string filename = "";
	TableIOOpts outOpts = TableIOOpts::genTabFileOut("out.tab.txt", false);
	setUp.setOption(filename, "--file", "Chrom index file name", true);
	setUp.processWritingOptions(outOpts.out_);
	setUp.finishSetUp(std::cout);

	GzSimpleBinFile<uint64_t, SlimCounterRef::IndexNumOfElements> gzReader(filename);
	auto datcon = gzReader.genDataContainer();
	uint32_t count = 0;
	std::ofstream outfile;
	openTextFile(outfile, outOpts.out_);
	while(gzReader.read(datcon) ){
		auto refPos = datcon.data_[0];
		auto cov = datcon.data_[1];
		auto fileLoc = datcon.data_[2];
		outfile << refPos << "\t" << cov << '\t' << fileLoc << std::endl;
		++count;
	}

	return 0;
}

int bamCountingExpRunner::coverageInfo(const njh::progutils::CmdArgs & inputCommands){
	seqSetUp setUp(inputCommands);
	std::string filename = "";
	TableIOOpts outOpts = TableIOOpts::genTabFileOut("", false);
	setUp.setOption(filename, "--file", "Chrom index file name", true);
	setUp.processWritingOptions(outOpts.out_);
	setUp.finishSetUp(std::cout);

	GzSimpleBinFile<uint64_t, SlimCounterRef::IndexNumOfElements> gzReader(filename);
	auto datcon = gzReader.genDataContainer();
	uint32_t count = 0;
	std::vector<uint32_t> coverage;
	while(gzReader.read(datcon) ){
		coverage.push_back(datcon.data_[1]);
		++count;
	}

	auto stats = getStatsOnVec(coverage);
	table outTab(stats,VecStr{"stat", "value"});
	outTab.outPutContents(outOpts);
	return 0;
}


int bamCountingExpRunner::gatherVariantsOnChrom(const njh::progutils::CmdArgs & inputCommands){
	seqSetUp setUp(inputCommands);
	std::string filename = "";
	std::string outFilename = "";
	TableIOOpts opts = TableIOOpts::genTabFileOut(outFilename);
	double fracCutOff = 0.01;
	double depthCutOff = 10;
	double strandBiasCutOff = 0.3;
	std::string sampleName = "";
	std::string twoBitFilename = "";
	std::string refIdLookupFile = "";
	//bool addIndels = true;
	VariantGatherer::GathVarPars varPars;
	bool noIndels = false;
	setUp.processWritingOptions(opts.out_);
	setUp.setOption(noIndels, "--noIndels", "No Indels");
	varPars.addIndels_ = !noIndels;
	setUp.setOption(varPars.onlyBiallelicVariants_, "--onlyBiallelic", "Only report biallelic variants");
	setUp.setOption(varPars.addNonRefSnps_, "--addNonRefSnps", "Add info on locations that don't have different variants but have 100% different base from ref");
	setUp.setOption(refIdLookupFile, "--refIdLookupFile", "Ref Id Lookup File", true);
	setUp.setOption(twoBitFilename, "--twoBitFilename", "TwoBit Filename", true);
	setUp.setOption(filename, "--file", "Filename", true);
	setUp.setOption(depthCutOff, "--depthCutOff", "Depth in order to be considered");
	setUp.setOption(fracCutOff, "--fracCutOff", "Variant Frequency Cut Off");
	setUp.setOption(sampleName, "--sampleName", "Sample Name", true);
	setUp.setOption(strandBiasCutOff, "--strandBiasCutOff", "Strand Bias Cut off, if the percentage of the coverage from one strand direction of a snp is less than this cut off, then exclude it");
	setUp.finishSetUp(std::cout);


	TwoBit::TwoBitFile twoBitfile(twoBitFilename);
	std::unordered_map<std::string, std::string> gSeqs;
	auto refInfos = tabToRefDataVec(table(TableIOOpts::genTabFileIn(refIdLookupFile, true)));

	VariantGatherer varGetter(refInfos,fracCutOff, depthCutOff, strandBiasCutOff);
	auto outTab = varGetter.gatherVariantsFromChromFile(filename, sampleName, varPars, twoBitfile);
	opts.outOrganized_ = true;
	outTab.outPutContents(opts);
	return 0;
}





int bamCountingExpRunner::gatherVariantsOnSample(const njh::progutils::CmdArgs & inputCommands){
	seqSetUp setUp(inputCommands);
	std::string dirname = "";
	TableIOOpts opts = TableIOOpts::genTabFileOut("");
	double fracCutOff = 0.01;
	double depthCutOff = 10;
	double strandBiasCutOff = 0.3;
	std::string sampleName = "";
	std::string twoBitFilename = "";
	bool noIndels = false;
	VariantGatherer::GathVarPars varPars;
	setUp.processWritingOptions(opts.out_);
	setUp.setOption(twoBitFilename, "--twoBitFilename", "TwoBit Filename", true);
	setUp.setOption(dirname, "--dirname", "Directory Output of SlimBamBaseCounter", true);
	setUp.setOption(sampleName, "--sampleName", "Sample Name", true);
	setUp.setOption(depthCutOff, "--depthCutOff", "Depth in order to be considered");
	setUp.setOption(fracCutOff, "--fracCutOff", "Variant Fraction Cut Off");
	setUp.setOption(noIndels, "--noIndels", "No Indels");
	varPars.addIndels_ = !noIndels;
	setUp.setOption(varPars.onlyBiallelicVariants_, "--onlyBiallelic", "Only report biallelic variants");
	setUp.setOption(varPars.addNonRefSnps_, "--addNonRefSnps", "Add info on locations that don't have different variants but have different 100% different base from ref");
	setUp.setOption(strandBiasCutOff, "--strandBiasCutOff", "Strand Bias Cut off, if the percentage of the coverage from one strand direction of a snp is less than this cut off, then exclude it");
	setUp.processVerbose();
	setUp.finishSetUp(std::cout);

	SlimCounterMaster::checkDirectoryStructure(dirname);
	auto refIdLookupFile = njh::files::make_path(dirname, "data", "refIdLookup.tab.txt");
	TwoBit::TwoBitFile twoBitfile(twoBitFilename);
	auto refInfos = tabToRefDataVec(table(TableIOOpts::genTabFileIn(refIdLookupFile.string(), true)));
	VariantGatherer varGetter(refInfos, fracCutOff, depthCutOff, strandBiasCutOff);
	table outTab;

	auto countsDir = njh::files::make_path(dirname, "data", "counts");
	auto countsFiles = njh::files::listAllFiles(countsDir.string(), false, {std::regex{R"(.*\.bgz$)"}});

	for (const auto & file : countsFiles) {
		if (setUp.pars_.verbose_) {
			std::cout << "Currently on " << file.first.string() << std::endl;
		}
		auto chromTab = varGetter.gatherVariantsFromChromFile(file.first.string(),
				sampleName, varPars, twoBitfile);
		if (outTab.empty()) {
			outTab = chromTab;
		} else {
			outTab.rbind(chromTab, false);
		}
	}
	opts.outOrganized_ = true;
	outTab.outPutContents(opts);
	return 0;
}


int bamCountingExpRunner::gatherBaseCoverageOnSampleOnLoc(const njh::progutils::CmdArgs & inputCommands){

	std::string dirname = "";
	TableIOOpts opts = TableIOOpts::genTabFileOut("");
	std::string sampleName = "";
	std::string bedFile = "";
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processWritingOptions(opts.out_);
	setUp.setOption(dirname, "--dirname", "Directory Output of SlimBamBaseCounter", true);
	setUp.setOption(sampleName, "--sampleName", "Sample Name", true);
	setUp.setOption(bedFile, "--bedFile", "bed file where the chromosome name and start will be used to report counts on", true);
	setUp.finishSetUp(std::cout);

	SlimCounterMaster::checkDirectoryStructure(dirname);
	auto refIdLookupFile = njh::files::make_path(dirname, "data", "refIdLookup.tab.txt");
	auto refInfos = tabToRefDataVec(table(TableIOOpts::genTabFileIn(refIdLookupFile.string(), true)));
	std::map<std::string, uint32_t> refLensByName;
	for(const auto & refInfo :refInfos){
		refLensByName[refInfo.RefName] = refInfo.RefLength;
	}
	std::map<std::string, uint32_t> refIdByName;
	for(const auto & refInfoPos : iter::range(refInfos.size())){
		refIdByName[refInfos[refInfoPos].RefName] = refInfoPos;
	}

	auto countsDir = njh::files::make_path(dirname, "data", "counts");
	auto countsFiles = njh::files::listAllFiles(countsDir.string(), false, {std::regex{R"(.*\.bgz$)"}});

	//organize the counts file by chromosome
	std::unordered_map<std::string, bfs::path> filesByChrom;
	for(const auto & file : countsFiles){
		filesByChrom[bfs::basename(file.first.string())] = file.first;
	}

	//organize the regions by chromosome
	auto gRegions = gatherRegions(bedFile, "", setUp.pars_.verbose_);
	std::map<std::string, std::vector<GenomicRegion>> regionsByChrom;
	for(const auto & reg : gRegions){
		regionsByChrom[reg.chrom_].emplace_back(reg);
	}

	for( auto & reg : regionsByChrom){
		njh::sort(reg.second, [](const GenomicRegion & reg1, const GenomicRegion & reg2){
			return reg1.start_ < reg2.start_;
		});
	}
	OutputStream out(opts.out_);
	out << "#chrom\tstart\tend\tname\tscore\tstrand\tsample\tminCov\tmaxCov\tmeanCov\tmedianCov" << std::endl;
	for(const auto & chrom : regionsByChrom){
		if(njh::in(chrom.first, filesByChrom)){
			auto index = SlimCounterMaster::readGzChromIndex(filesByChrom[chrom.first].string() + ".idx.gz");
			BGZFCPP bgzReader(filesByChrom[chrom.first], "r");
			//GzSimpleBinFile<uint32_t, SlimCounterPos::NumOfElements> inFile(filesByChrom[chrom.first]);
			std::vector<uint32_t> d(SlimCounterPos::NumOfElements);
			for(const auto & reg : chrom.second){
				if(setUp.pars_.debug_){
					std::cout << reg.uid_ << std::endl;
					std::cout << reg.start_ << std::endl;
					std::cout << reg.end_ << std::endl;
				}
				auto rPosLast = reg.start_;
				bool haveSought = false;
				std::vector<uint32_t> baseCoverages;
				for(const auto & rPos : iter::range(reg.start_, reg.end_)){
					if (njh::in(rPos, index)) {
						if(rPosLast + 1 != rPos || !haveSought) {
							bgzReader.seek(index[rPos]);
							haveSought = true;
						}
						bgzReader.readVec(d);
						SlimCounterPos pos(d);
						if(setUp.pars_.debug_){
							std::cout << "rPos: " << rPos << std::endl;
							std::cout << "refPos: " << d[1] << std::endl;
						}
						//uint32_t refPos = d[1];
						auto rName = refInfos[d[0]].RefName;
						baseCoverages.emplace_back(pos.getBaseTotal());
					} else {
						if(refLensByName[reg.chrom_] > rPos){
							rPosLast = rPos;
							continue;
						}
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ": Error position: " << rPos
								<< " is greater than the length of the chromosome: "
								<< reg.chrom_ << ", len:  " << refLensByName[reg.chrom_]
								<< "\n";
						throw std::runtime_error { ss.str() };
					}
					rPosLast = rPos;
				}
				out << reg.genBedRecordCore().toDelimStr()
						<< "\t" << sampleName
						<< "\t" << vectorMinimum(baseCoverages)
						<< "\t" << vectorMaximum(baseCoverages)
						<< "\t" << vectorMean(baseCoverages)
						<< "\t" << vectorMedianRef(baseCoverages)
						<< std::endl;
			}
		}else{
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ": No count file found for " << chrom.first << "\n"
					<< "options are: " << njh::conToStr(getVectorOfMapKeys(filesByChrom), ", ") << "\n";
			throw std::runtime_error { ss.str() };
		}
	}
	return 0;
}

int bamCountingExpRunner::gatherVariantsOnSampleOnLoc(const njh::progutils::CmdArgs & inputCommands){

	std::string dirname = "";
	TableIOOpts opts = TableIOOpts::genTabFileOut("");
	double fracCutOff = 0.01;
	double depthCutOff = 10;
	double strandBiasCutOff = 0.3;
	std::string sampleName = "";
	std::string twoBitFilename = "";
	std::string bedFile = "";
	std::string gffFile = "";
	bool noIndels = false;
	VariantGatherer::GathVarPars varPars;
	bool addFiller = false;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processWritingOptions(opts.out_);
	setUp.setOption(twoBitFilename, "--twoBitFilename", "TwoBit Filename", true);
	setUp.setOption(dirname, "--dirname", "Directory Output of SlimBamBaseCounter", true);
	setUp.setOption(sampleName, "--sampleName", "Sample Name", true);
	setUp.setOption(depthCutOff, "--depthCutOff", "Depth in order to be considered");
	setUp.setOption(fracCutOff, "--fracCutOff", "Variant Fraction Cut Off");
	setUp.setOption(noIndels, "--noIndels", "No Indels");
	setUp.setOption(addFiller, "--addFiller", "If a location doesn't have reads mapped to it output empty values so it's still in the output file");

	varPars.addIndels_ = !noIndels;
	setUp.setOption(varPars.onlyBiallelicVariants_, "--onlyBiallelic", "Only report biallelic variants");
	setUp.setOption(varPars.addNonRefSnps_, "--addNonRefSnps", "Add info on locations that don't have different variants but have different 100% different base from ref");
	setUp.setOption(strandBiasCutOff, "--strandBiasCutOff", "Strand Bias Cut off, if the percentage of the coverage from one strand direction of a snp is less than this cut off, then exclude it");

	setUp.setOption(bedFile, "--bedFile", "bed file where the chromosome name and start will be used to report counts on");
	setUp.setOption(gffFile, "--gffFile", "gff file where the chromosome name and start will be used to report counts on");
	if("" == bedFile && "" == gffFile){
		setUp.failed_ = true;
		setUp.addWarning("Need to supply either a gff (--gffFile) or bed (--bedFile) file");
	}
	setUp.finishSetUp(std::cout);

	SlimCounterMaster::checkDirectoryStructure(dirname);
	auto refIdLookupFile = njh::files::make_path(dirname, "data", "refIdLookup.tab.txt");
	TwoBit::TwoBitFile twoBitfile(twoBitFilename);
	auto refInfos = tabToRefDataVec(table(TableIOOpts::genTabFileIn(refIdLookupFile.string(), true)));
	std::map<std::string, uint32_t> refLensByName;
	for(const auto & refInfo :refInfos){
		refLensByName[refInfo.RefName] = refInfo.RefLength;
	}
	std::map<std::string, uint32_t> refIdByName;
	for(const auto & refInfoPos : iter::range(refInfos.size())){
		refIdByName[refInfos[refInfoPos].RefName] = refInfoPos;
	}
	VariantGatherer varGetter(refInfos, fracCutOff, depthCutOff, strandBiasCutOff);

	auto countsDir = njh::files::make_path(dirname, "data", "counts");
	auto countsFiles = njh::files::listAllFiles(countsDir.string(), false, {std::regex{R"(.*\.bgz$)"}});

	//organize the counts file by chromosome
	std::unordered_map<std::string, bfs::path> filesByChrom;
	for(const auto & file : countsFiles){
		filesByChrom[bfs::basename(file.first.string())] = file.first;
	}

	//organize the regions by chromosome
	auto gRegions = gatherRegions(bedFile, gffFile, setUp.pars_.verbose_);
	std::map<std::string, std::vector<GenomicRegion>> regionsByChrom;
	for(const auto & reg : gRegions){
		regionsByChrom[reg.chrom_].emplace_back(reg);
	}

	for( auto & reg : regionsByChrom){
		njh::sort(reg.second, [](const GenomicRegion & reg1, const GenomicRegion & reg2){
			return reg1.start_ < reg2.start_;
		});
	}
	std::string genomicSeq;
	table outTab(
				concatVecs(VecStr { "sample", "alleleCount" }, SlimCounterPos::hqDetailHeader()));
	for(const auto & chrom : regionsByChrom){

		if(njh::in(chrom.first, filesByChrom)){
			twoBitfile[chrom.first]->getSequence(genomicSeq);
			auto index = SlimCounterMaster::readGzChromIndex(filesByChrom[chrom.first].string() + ".idx.gz");
			if(setUp.pars_.debug_){
				std::ofstream outFile("temp" + chrom.first + "_debug.tab.txt");
				auto keys = njh::getVecOfMapKeys(index);
				njh::sort(keys);
				for(const auto key : keys){
					outFile << key << "\t" << index[key] << std::endl;
				}
			}
			BGZFCPP bgzReader(filesByChrom[chrom.first], "r");
			//GzSimpleBinFile<uint32_t, SlimCounterPos::NumOfElements> inFile(filesByChrom[chrom.first]);
			std::vector<uint32_t> d(SlimCounterPos::NumOfElements);
			for(const auto & reg : chrom.second){
				if(setUp.pars_.debug_){
					std::cout << reg.uid_ << std::endl;
					std::cout << reg.start_ << std::endl;
					std::cout << reg.end_ << std::endl;
				}
				auto rPosLast = reg.start_;
				bool haveSought = false;
				for(const auto & rPos : iter::range(reg.start_, reg.end_)){
					if (njh::in(rPos, index)) {
						if(rPosLast + 1 != rPos || !haveSought) {
							bgzReader.seek(index[rPos]);
							haveSought = true;
						}
						bgzReader.readVec(d);
						SlimCounterPos pos(d);
						if(setUp.pars_.debug_){
							std::cout << "rPos: " << rPos << std::endl;
							std::cout << "refPos: " << d[1] << std::endl;
						}
						uint32_t refPos = d[1];
						auto rName = refInfos[d[0]].RefName;
						const char refBase = genomicSeq[refPos];
						auto varResults = varGetter.determineIfPosVariant(pos, refBase, varPars);
						if(varResults.variant_){
							for (const auto & row : pos.hqDetailOutVec(rName, refPos,
									refBase)) {
								outTab.addRow(
										concatVecs(toVecStr(sampleName, varResults.numberOfAlleles_), row));
							}
						}
					} else {
						if(refLensByName[reg.chrom_] > rPos){
							rPosLast = rPos;
							if(addFiller){
								std::fill(d.begin() + 2, d.end(), 0);
								d[0] = refIdByName[reg.chrom_];
								d[1] = rPos;
								SlimCounterPos pos(d);
								uint32_t refPos = d[1];
								auto rName = refInfos[d[0]].RefName;
								for (const auto & row : pos.hqDetailOutVec(rName, refPos,
										genomicSeq[refPos])) {
									outTab.addRow(
											concatVecs(VecStr { sampleName, reg.uid_ }, row));
								}
							}
							continue;
						}
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ": Error position: " << rPos
								<< " is greater than the length of the chromosome: "
								<< reg.chrom_ << ", len:  " << refLensByName[reg.chrom_]
								<< "\n";
						throw std::runtime_error { ss.str() };
					}
					rPosLast = rPos;
				}
			}
		}else{
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ": No count file found for " << chrom.first << "\n"
					<< "options are: " << njh::conToStr(getVectorOfMapKeys(filesByChrom), ", ") << "\n";
			throw std::runtime_error { ss.str() };
		}
	}
	opts.outOrganized_ = true;
	outTab.outPutContents(opts);
	return 0;
}

int bamCountingExpRunner::gatherLocOnSample(const njh::progutils::CmdArgs & inputCommands){

	std::string dirname = "";
	std::string outFilename = "";
	TableIOOpts opts = TableIOOpts::genTabFileOut(outFilename);
	std::string sampleName = "";
	std::string twoBitFilename = "";
	std::string bedFile = "";
	std::string gffFile = "";
	bool addFiller = false;
	bool onlyRef = false;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processWritingOptions(opts.out_);
	setUp.setOption(twoBitFilename, "--twoBitFilename", "TwoBit Filename", true);
	setUp.setOption(dirname, "--dirname", "Directory Output of SlimBamBaseCounter", true);
	setUp.setOption(sampleName, "--sampleName", "Sample Name", true);
	setUp.setOption(addFiller, "--addFiller", "If a location doesn't have reads mapped to it output empty values so it's still in the output file");
	setUp.setOption(onlyRef, "--onlyRef", "Only report on reference position base");
	setUp.setOption(bedFile, "--bedFile", "bed file to report counts on");
	setUp.setOption(gffFile, "--gffFile", "gff file to report counts on");
	if("" == bedFile && "" == gffFile){
		setUp.failed_ = true;
		setUp.addWarning("Need to supply either a gff (--gffFile) or bed (--bedFile) file");
	}
	setUp.finishSetUp(std::cout);

	SlimCounterMaster::checkDirectoryStructure(dirname);
	auto refIdLookupFile = njh::files::make_path(dirname, "data", "refIdLookup.tab.txt");

	TwoBit::TwoBitFile twoBitfile(twoBitFilename);
	std::unordered_map<std::string, std::string> gSeqs;

	auto gRegions = gatherRegions(bedFile, gffFile, setUp.pars_.verbose_);
	auto refInfos = tabToRefDataVec(table(TableIOOpts::genTabFileIn(refIdLookupFile.string(), true)));
	std::map<std::string, uint32_t> refLensByName;
	for(const auto & refInfo :refInfos){
		refLensByName[refInfo.RefName] = refInfo.RefLength;
	}
	std::map<std::string, uint32_t> refIdByName;
	for(const auto & refInfoPos : iter::range(refInfos.size())){
		refIdByName[refInfos[refInfoPos].RefName] = refInfoPos;
	}

	auto countsDir = njh::files::make_path(dirname, "data", "counts");
	auto countsFiles = njh::files::listAllFiles(countsDir.string(), false, {std::regex{R"(.*\.bgz$)"}});

	//organize the counts file by chromosome
	std::unordered_map<std::string, bfs::path> filesByChrom;
	for(const auto & file : countsFiles){
		filesByChrom[bfs::basename(file.first.string())] = file.first;
	}

	//organize the regions by chromosome
	std::map<std::string, std::vector<GenomicRegion>> regionsByChrom;
	for(const auto & reg : gRegions){
		regionsByChrom[reg.chrom_].emplace_back(reg);
	}

	for( auto & reg : regionsByChrom){
		njh::sort(reg.second, [](const GenomicRegion & reg1, const GenomicRegion & reg2){
			return reg1.start_ < reg2.start_;
		});
	}

	std::string genomicSeq;
	table outTab(concatVecs(VecStr{"sample", "locName"},SlimCounterPos::hqDetailHeader()));
	for(const auto & chrom : regionsByChrom){
		if(njh::in(chrom.first, filesByChrom)){
			twoBitfile[chrom.first]->getSequence(genomicSeq);
			auto index = SlimCounterMaster::readGzChromIndex(filesByChrom[chrom.first].string() + ".idx.gz");
			if(setUp.pars_.debug_){
				std::ofstream outFile("temp" + chrom.first + "_debug.tab.txt");
				auto keys = njh::getVecOfMapKeys(index);
				njh::sort(keys);
				for(const auto key : keys){
					outFile << key << "\t" << index[key] << std::endl;
				}
			}
			BGZFCPP bgzReader(filesByChrom[chrom.first], "r");
			//GzSimpleBinFile<uint32_t, SlimCounterPos::NumOfElements> inFile(filesByChrom[chrom.first]);
			std::vector<uint32_t> d(SlimCounterPos::NumOfElements);
			for(const auto & reg : chrom.second){
				if(setUp.pars_.debug_){
					std::cout << reg.uid_ << std::endl;
					std::cout << reg.start_ << std::endl;
					std::cout << reg.end_ << std::endl;
				}
				auto rPosLast = reg.start_;
				bool haveSought = false;
				for(const auto & rPos : iter::range(reg.start_, reg.end_)){
					if (njh::in(rPos, index)) {
						if(rPosLast + 1 != rPos || !haveSought) {
							bgzReader.seek(index[rPos]);
							haveSought = true;
						}
						bgzReader.readVec(d);
						SlimCounterPos pos(d);
						if(setUp.pars_.debug_){
							std::cout << "rPos: " << rPos << std::endl;
							std::cout << "refPos: " << d[1] << std::endl;
						}
						uint32_t refPos = d[1];
						auto rName = refInfos[d[0]].RefName;
						for (const auto & row : pos.hqDetailOutVec(rName, refPos,
								genomicSeq[refPos])) {
							outTab.addRow(
									concatVecs(VecStr { sampleName, reg.uid_ }, row));
						}
					} else {
						if(refLensByName[reg.chrom_] > rPos){
							rPosLast = rPos;
							if(addFiller){
								std::fill(d.begin() + 2, d.end(), 0);
								d[0] = refIdByName[reg.chrom_];
								d[1] = rPos;
								SlimCounterPos pos(d);
								uint32_t refPos = d[1];
								auto rName = refInfos[d[0]].RefName;
								for (const auto & row : pos.hqDetailOutVec(rName, refPos,
										genomicSeq[refPos])) {
									outTab.addRow(
											concatVecs(VecStr { sampleName, reg.uid_ }, row));
								}
							}
							continue;
						}
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ": Error position: " << rPos
								<< " is greater than the length of the chromosome: "
								<< reg.chrom_ << ", len:  " << refLensByName[reg.chrom_]
								<< "\n";
						throw std::runtime_error { ss.str() };
					}
					rPosLast = rPos;
				}
			}
		}else{
			if(addFiller){
				twoBitfile[chrom.first]->getSequence(genomicSeq);
				for(const auto & reg : chrom.second){
					for(const auto & rPos : iter::range(reg.start_, reg.end_)){
						std::vector<uint32_t> d(SlimCounterPos::NumOfElements);
						std::fill(d.begin() + 2, d.end(), 0);
						d[0] = refIdByName[reg.chrom_];
						d[1] = rPos;
						SlimCounterPos pos(d);
						uint32_t refPos = d[1];
						auto rName = refInfos[d[0]].RefName;
						for (const auto & row : pos.hqDetailOutVec(rName, refPos,
								genomicSeq[refPos])) {
							outTab.addRow(
									concatVecs(VecStr { sampleName, reg.uid_ }, row));
						}
					}
				}
			}else{
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ": No count file found for " << chrom.first << "\n"
						<< "options are: " << njh::conToStr(getVectorOfMapKeys(filesByChrom), ", ") << "\n";
				throw std::runtime_error { ss.str() };
			}
		}
	}

	opts.outOrganized_ = true;
	//outTab.sortTable("refId", "pos", "base", false);
	if(onlyRef){
		table selectedTab(outTab.columnNames_);
		for(const auto & row : outTab.content_){
			if(row[outTab.getColPos("refBase")] == row[outTab.getColPos("base")] ){
				selectedTab.addRow(row);
			}
		}
		selectedTab.outPutContents(opts);
	}else{
		outTab.outPutContents(opts);
	}

	return 0;
}



//

int bamCountingExpRunner::bamBaseCountingSlim(const njh::progutils::CmdArgs & inputCommands){
	uint32_t testNumber = std::numeric_limits<uint32_t>::max();
	uint32_t minLen = 20;
	uint32_t insertLengthCutOff = 2000;
	uint32_t mappingQual = 10;
	uint32_t qualCutOff = 20;
	uint32_t readPosCutOff = 0;
	uint32_t numThreads = 1;
	bool ionTorrent = false;
	bool reAlignTorrent = false;
	size_t hRunSize = 3;
	bool countAlones = false;
	bool noOverlapHandling = false;
	seqSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processVerbose();
	setUp.processAlnInfoInput();
	setUp.setOption(noOverlapHandling, "--noOverlapHandling", "Don't handle overlap for counting");
	setUp.setOption(countAlones, "--countAlones",
			"Count reads where their mates couldn't be found");
	setUp.setOption(numThreads, "--numThreads", "Number of Threads to Use");
	setUp.setOption(readPosCutOff, "--readPosCutOff",
			"Don't count bases appear this many bases from the front of the end of the mapping reads");
	setUp.processDefaultReader(VecStr { "--bam" }, true);
	setUp.setOption(ionTorrent, "--ionTorrent",
			"Whether or not the sequencing technology was Ion Torrent");
	setUp.setOption(reAlignTorrent, "--reAlignTorrent",
			"Whether the sequences should be realigned when sequencing technology was Ion Torrent");
	setUp.setOption(hRunSize, "--hRunSize",
			"Set the size of the homopolymer run when the --ionTorrent flag is used to test both the first and last base for quality, must be at least be 3");
	if (hRunSize < 3) {
		setUp.failed_ = true;
		setUp.addWarning(
				"Error, --hRunSize must be at least 3, was set to "
						+ estd::to_string(hRunSize));
	}
	setUp.setOption(mappingQual, "--mappingQual",
			"Mapping quality score cut off for a read");
	setUp.setOption(qualCutOff, "--qualCutOff",
			"Base Quality score cut off for a snp");
	setUp.setOption(minLen, "--minLen", "Minimum Length of read");
	setUp.setOption(testNumber, "--testNumber",
			"Testing Number of reads to test out the running of the program");
	setUp.pars_.gap_ = "5,1";
	setUp.pars_.gapLeft_ = "5,1";
	setUp.pars_.gapRight_ = "5,1";
	setUp.processGap();
	setUp.processScoringPars();
	setUp.processDirectoryOutputName(true);
	setUp.finishSetUp(std::cout);

	setUp.startARunLog(setUp.pars_.directoryName_);
	BamTools::BamReader bReader;
	BamTools::BamAlignment bAln;
	bReader.Open(setUp.pars_.ioOptions_.firstName_.string());
	checkBamOpenThrow(bReader, setUp.pars_.ioOptions_.firstName_.string());
	loadBamIndexThrow(bReader);
	auto refs = bReader.GetReferenceData();

	if(setUp.pars_.verbose_){
		std::cout << "Creating BamReader Pool" << std::endl;
	}
	concurrent::BamReaderPool bamPool(setUp.pars_.ioOptions_.firstName_, numThreads);
	if(setUp.pars_.verbose_){
		std::cout << "Initializing BamReader Pool" << std::endl;
	}
	bamPool.openBamFile();

	if(setUp.pars_.verbose_){
		std::cout << "Creating Aligner Pool" << std::endl;
	}
	uint64_t maxLen = 800;
	concurrent::AlignerPool alnPool(maxLen, setUp.pars_.gapInfo_, setUp.pars_.scoring_, numThreads);
	if(setUp.pars_.verbose_){
		std::cout << "Initializing Aligner Pool" << std::endl;
	}
	alnPool.initAligners();

	if(setUp.pars_.verbose_){
		std::cout << "Creating chroms queue" << std::endl;
	}

	njh::concurrent::LockableQueue<BamTools::RefData> chroms(refs);


	std::unordered_map<std::string, BamCountExtractStats> refExtracts;
	std::mutex refExtractsMut;

	SlimCounterMaster masterCounter(setUp.pars_.directoryName_ + "data", refs);

	auto countFunc =
			[mappingQual, minLen,
			insertLengthCutOff, qualCutOff, ionTorrent, reAlignTorrent, hRunSize,
			readPosCutOff, testNumber, countAlones, noOverlapHandling](
					concurrent::BamReaderPool & bamPool, concurrent::AlignerPool & alnPool,
					std::unordered_map<std::string, BamCountExtractStats>& refExtracts,
					std::mutex & refExtractsMut,
					SlimCounterMaster & masterCounter,
					njh::concurrent::LockableQueue<BamTools::RefData> & chroms,
					bool debug, bool verbose) {

				if(verbose) {
					std::lock_guard<std::mutex> statsLock(refExtractsMut);
					std::cout <<"[TID=" << std::this_thread::get_id() << "]: Starting Thread" << std::endl;
				}
				BamTools::RefData refId;
				if(verbose) {
					std::lock_guard<std::mutex> statsLock(refExtractsMut);
					std::cout <<"[TID=" << std::this_thread::get_id() << "]: Getting Reader" << std::endl;
				}
				auto reader = bamPool.popReader();
				if(verbose) {
					std::lock_guard<std::mutex> statsLock(refExtractsMut);
					std::cout <<"[TID=" << std::this_thread::get_id() << "]: Getting Aligner" << std::endl;
				}
				auto alignerObj = alnPool.popAligner();
				if(verbose) {
					std::lock_guard<std::mutex> statsLock(refExtractsMut);
					std::cout <<"[TID=" << std::this_thread::get_id() << "]: Starting extraction" << std::endl;
				}
				std::shared_ptr<SlimCounterRef> currentRef = nullptr;
				BamAlnsCache alnCache;
				while(chroms.getVal(refId)) {
					if(verbose) {
						std::lock_guard<std::mutex> statsLock(refExtractsMut);
						std::cout <<"[TID=" << std::this_thread::get_id() << "]: Extracting on " << refId.RefName << std::endl;
					}
					reader->SetRegion(reader->GetReferenceID(refId.RefName), 0, reader->GetReferenceID(refId.RefName), refId.RefLength);
					BamCountExtractStats extractStats;
					BamBaseCountingSlim(*reader, *alignerObj, mappingQual, minLen,
							insertLengthCutOff, qualCutOff, ionTorrent, reAlignTorrent, hRunSize,
							readPosCutOff, masterCounter, currentRef, extractStats, alnCache,
							testNumber,countAlones, debug, verbose, noOverlapHandling);
					{
						std::lock_guard<std::mutex> statsLock(refExtractsMut);
						if(verbose){
							std::cout <<"[TID=" << std::this_thread::get_id() << "]: Adding extraction stats for " << refId.RefName << std::endl;
						}
						refExtracts[refId.RefName] = extractStats;
					}
					if(nullptr != currentRef){
						currentRef->flushCountsToFile(0, std::numeric_limits<uint32_t>::max());
					}
					if(verbose){
						std::lock_guard<std::mutex> statsLock(refExtractsMut);
						std::cout <<"[TID=" << std::this_thread::get_id() << "]: Done Extracting on " << refId.RefName << std::endl;
					}
				}
			};
	std::vector<std::thread> threads;

	for (uint32_t threadNum = 0; threadNum < numThreads; ++threadNum) {
		threads.emplace_back(
						std::thread(countFunc,std::ref(bamPool),std::ref(alnPool),
								std::ref(refExtracts),std::ref(refExtractsMut),
								std::ref(masterCounter),std::ref(chroms),
								setUp.pars_.debug_,setUp.pars_.verbose_));
	}

	for(auto & thread : threads){
		thread.join();
	}
	if(setUp.pars_.debug_){
		std::cout << "Finished Extraction, setting master counter" << std::endl;
	}


	if(setUp.pars_.verbose_){
		std::cout << std::endl;
	}



	bfs::path statsDir = njh::files::makeDir(setUp.pars_.directoryName_,
			njh::files::MkdirPar("stats"));

	BamCountExtractStats::createPairedReadInfoFileMultiple(refExtracts,
			OutOptions(bfs::path(statsDir.string() + "pairedReadsInfo.tab.txt")));

	BamCountExtractStats::createBaseFilterFileMultiple(refExtracts,
			OutOptions(bfs::path(statsDir.string() + "baseFilterInfo.tab.txt")),
			qualCutOff);

	BamCountExtractStats::createReadFilterFileMultiple(refExtracts,
			OutOptions(bfs::path(statsDir.string() + "readFilterInfo.tab.txt")), minLen,
			mappingQual, insertLengthCutOff);

	for(const auto & ref : masterCounter.counts_){
		ref.second->flushCountsToFile(0, std::numeric_limits<uint32_t>::max());
		if(bfs::exists(ref.second->outputFile_)){
			ref.second->indexFile();
		}
	}

	/*alignerObj.processAlnInfoOutput(setUp.pars_.outAlnInfoDirName_,
			setUp.pars_.verbose_);*/
	return 0;
}

} /* namespace njhseq */
