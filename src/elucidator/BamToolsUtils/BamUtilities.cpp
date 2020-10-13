/*
 * BamUtilities.cpp
 *
 *  Created on: Jun 10, 2018
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


#include "BamUtilities.hpp"
#include <njhseq/concurrency/pools/BamReaderPool.hpp>
#include <njhseq/objects/BioDataObject.h>
#include <njhseq/readVectorManipulation/readVectorHelpers/readVecTrimmer.hpp>

namespace njhseq {


void RunCoverageFinderMulti(const RunCoverageFinderMultiPars & pars){
	VecStr header = {"#chrom", "start", "end", "name", "score", "strand"};

	std::vector<bfs::path> bamFnps;
	if ("" == pars.bams || std::string::npos != pars.bams.find(",")) {
		bamFnps = njh::files::gatherFilesByPatOrNames(std::regex { pars.pat },
				pars.bams);
	} else {
		bamFnps = {pars.bams};
	}

	checkBamFilesForIndexesAndAbilityToOpen(bamFnps);
	OutputStream out(pars.outOpts);
	std::mutex outMut;
	out << njh::conToStr(header, "\t");
	for(const auto & bamFnp : bamFnps){
		out << "\t" << bamFnp.filename().string();
	}
	out << std::endl;



	class GenomeRegionSlider {
	public:
		GenomeRegionSlider(const BamTools::RefVector & refData, uint32_t windowSize,
				uint32_t step) :
				refLengths_(refData), windowSize_(windowSize), step_(step) {
		}
		BamTools::RefVector refLengths_;
		uint32_t windowSize_;
		uint32_t step_;
	private:
		uint32_t currentPos_ = 0;
		uint32_t currentRefId_ = 0;
		std::mutex mut_;
		bool getRegionNoLock(GenomicRegion & region) {
			region.chrom_ = "*";
			if (currentRefId_ < refLengths_.size()) {
				if (currentPos_ + windowSize_ < static_cast<uint64_t>(refLengths_[currentRefId_].RefLength)) {
					region = GenomicRegion("", refLengths_[currentRefId_].RefName, currentPos_, currentPos_ + windowSize_, false);
					region.setUidWtihCoords();
					currentPos_ += step_;
					return true;
				} else {
					++currentRefId_;
					currentPos_ = 0;
					if (currentRefId_ < refLengths_.size()) {
						if (currentPos_ + windowSize_
								< static_cast<uint64_t>(refLengths_[currentRefId_].RefLength)) {
							region = GenomicRegion("", refLengths_[currentRefId_].RefName, currentPos_, currentPos_ + windowSize_, false);
							region.setUidWtihCoords();
							currentPos_ += step_;
							return true;
						}
					}
				}
			}
			return false;
		}

	public:

		bool getRegion(GenomicRegion & region) {
			std::lock_guard<std::mutex> lock(mut_);
			return getRegionNoLock(region);
		}

		bool getRegions(std::vector<GenomicRegion> & regions, uint32_t batchSize){
			std::lock_guard<std::mutex> lock(mut_);
			regions.clear();
			bool added = false;
			for(uint32_t batch = 0; batch < batchSize; ++batch){
				GenomicRegion region;
				if(getRegionNoLock(region)){
					regions.emplace_back(region);
					added = true;
				}else{
					break;
				}
			}
			return added;
		}

	};

	struct BamFnpRegionPair {
		BamFnpRegionPair(const bfs::path & bamFnp, const GenomicRegion & region) :
				bamFnp_(bamFnp), region_(region) {
			regUid_ = region_.createUidFromCoords();
			bamFname_ = bamFnp_.filename().string();
		}
		bfs::path bamFnp_;
		GenomicRegion region_;
		uint32_t coverage_ = 0;

		std::string regUid_;
		std::string bamFname_;
	};

	BamTools::BamReader bReaderForRefData;
	bReaderForRefData.Open(bamFnps.front().string());
	/**@todo should do a check that all bam fnps have the same ref data*/
	GenomeRegionSlider slider(bReaderForRefData.GetReferenceData(), pars.window, pars.step);

	auto refDataRaw = bReaderForRefData.GetReferenceData();
	BamTools::RefVector refData;
	for(const auto & refDatum : refDataRaw){
		if(!njh::in(refDatum.RefName, pars.chromsToSkip)){
			refData.emplace_back(refDatum);
		}
	}
	std::function<void()> getCov = [&slider,&bamFnps,&outMut,&out,&refData, &pars](){
		std::vector<GenomicRegion> regions;
		while(slider.getRegions(regions, pars.regionBatchSize)){
			std::vector<std::vector<BamFnpRegionPair>> allCoverages;
			for(const auto & region : regions){
				std::vector<BamFnpRegionPair> coverages;
				for(const auto & bamFnp : bamFnps){
					BamFnpRegionPair coverage(bamFnp, region);
					BamTools::BamReader bReader;
					bReader.Open(bamFnp.string());
					bReader.LocateIndex();
					setBamFileRegionThrow(bReader, region);
					BamTools::BamAlignment bAln;
					while(bReader.GetNextAlignmentCore(bAln) && bAln.IsPrimaryAlignment()){
						if(pars.byBases){
							if(bAln.IsMapped()){
								GenomicRegion alnRegion(bAln, refData);
								coverage.coverage_ += region.getOverlapLen(alnRegion);
							}
						}else{
							if(bAln.IsMapped() && bAln.IsPrimaryAlignment()){
								coverage.coverage_ += 1;
							}
						}
					}
					bool pass = false;
					if(pars.byBases){
						if(coverage.coverage_/static_cast<double>(coverage.region_.getLen()) >= pars.coverageCutOff){
							coverages.emplace_back(coverage);
							pass = true;
						}
					}else{
						if(coverage.coverage_ >= pars.coverageCutOff){
							coverages.emplace_back(coverage);
							pass = true;
						}
					}
					if(!pass){
						break;
					}
				}
				if(coverages.size() == bamFnps.size()){
					allCoverages.emplace_back(coverages);
				}
			}
			if(!allCoverages.empty()){
				std::lock_guard<std::mutex> lock(outMut);
				for(const auto & allCoverage : allCoverages){
					out << allCoverage.front().region_.genBedRecordCore().toDelimStr() ;
					for(const auto & cov : allCoverage){
						out << "\t" << cov.coverage_;
					}
					out << std::endl;
				}
			}
		}
	};
	{
		njh::concurrent::runVoidFunctionThreaded(getCov, pars.numThreads);
	}
}


void RunCoverageFinderSingle(const RunCoverageFinderSinglePars & pars){
	VecStr header = {"#chrom", "start", "end", "name", "score", "strand"};


	std::vector<bfs::path> bamFnps{pars.bamFnp};
	OutputStream out(pars.outOpts);
	std::mutex outMut;
	out << njh::conToStr(header, "\t");
	for(const auto & bamFnp : bamFnps){
		out << "\t" << bamFnp.filename().string();
	}
	out << std::endl;



	class GenomeRegionSlider {
	public:
		GenomeRegionSlider(const BamTools::RefVector & refData, uint32_t windowSize,
				uint32_t step) :
				refLengths_(refData), windowSize_(windowSize), step_(step) {
		}
		BamTools::RefVector refLengths_;
		uint32_t windowSize_;
		uint32_t step_;
	private:
		uint32_t currentPos_ = 0;
		uint32_t currentRefId_ = 0;
		std::mutex mut_;
		bool getRegionNoLock(GenomicRegion & region) {
			region.chrom_ = "*";
			if (currentRefId_ < refLengths_.size()) {
				if (currentPos_ + windowSize_ < static_cast<uint64_t>(refLengths_[currentRefId_].RefLength)) {
					region = GenomicRegion("", refLengths_[currentRefId_].RefName, currentPos_, currentPos_ + windowSize_, false);
					region.setUidWtihCoords();
					currentPos_ += step_;
					return true;
				} else {
					++currentRefId_;
					currentPos_ = 0;
					if (currentRefId_ < refLengths_.size()) {
						if (currentPos_ + windowSize_
								< static_cast<uint64_t>(refLengths_[currentRefId_].RefLength)) {
							region = GenomicRegion("", refLengths_[currentRefId_].RefName, currentPos_, currentPos_ + windowSize_, false);
							region.setUidWtihCoords();
							currentPos_ += step_;
							return true;
						}
					}
				}
			}
			return false;
		}

	public:

		bool getRegion(GenomicRegion & region) {
			std::lock_guard<std::mutex> lock(mut_);
			return getRegionNoLock(region);
		}

		bool getRegions(std::vector<GenomicRegion> & regions, uint32_t batchSize){
			std::lock_guard<std::mutex> lock(mut_);
			regions.clear();
			bool added = false;
			for(uint32_t batch = 0; batch < batchSize; ++batch){
				GenomicRegion region;
				if(getRegionNoLock(region)){
					regions.emplace_back(region);
					added = true;
				}else{
					break;
				}
			}
			return added;
		}

	};

	struct BamFnpRegionPair {
		BamFnpRegionPair(const bfs::path & bamFnp, const GenomicRegion & region) :
				bamFnp_(bamFnp), region_(region) {
			regUid_ = region_.createUidFromCoords();
			bamFname_ = bamFnp_.filename().string();
		}
		bfs::path bamFnp_;
		GenomicRegion region_;
		uint32_t coverage_ = 0;

		std::string regUid_;
		std::string bamFname_;
	};


	auto bamFnp = pars.bamFnp;
	njhseq::concurrent::BamReaderPool bPool(pars.bamFnp, pars.numThreads);
	bPool.openBamFile();
	BamTools::RefVector refData;
	{
		auto bReader = bPool.popReader();
		BamTools::RefVector refDataRaw = bReader->GetReferenceData();
		;
		for(const auto & refDatum : refDataRaw){
			if(!njh::in(refDatum.RefName, pars.chromsToSkip)){
				refData.emplace_back(refDatum);
			}
		}
	}
	GenomeRegionSlider slider(refData, pars.window, pars.step);
	std::function<void()> getCovForBam = [&slider,&bamFnp,
																				&outMut,&out,
																				&refData,&pars,
																				&bPool](){
		auto bReader = bPool.popReader();
		std::vector<GenomicRegion> regions;
		while(slider.getRegions(regions, pars.regionBatchSize)){

			std::vector<BamFnpRegionPair> allCoverages;
			for(const auto & region : regions){
				BamFnpRegionPair coverage(bamFnp, region);
				setBamFileRegionThrow(*bReader, region);
				BamTools::BamAlignment bAln;
				while(bReader->GetNextAlignmentCore(bAln) && bAln.IsPrimaryAlignment()){
					if(pars.byBases){
						if(bAln.IsMapped()){
							GenomicRegion alnRegion(bAln, refData);
							coverage.coverage_ += region.getOverlapLen(alnRegion);
						}
					}else{
						if(bAln.IsMapped() && bAln.IsPrimaryAlignment()){
							coverage.coverage_ += 1;
						}
					}
				}
				if(pars.byBases){
					if(coverage.coverage_/static_cast<double>(coverage.region_.getLen()) >= pars.coverageCutOff){
						allCoverages.emplace_back(coverage);
					}
				}else{
					if(coverage.coverage_ >= pars.coverageCutOff){
						allCoverages.emplace_back(coverage);
					}
				}
			}
			if(!allCoverages.empty()){
				std::lock_guard<std::mutex> lock(outMut);
				for(const auto & coverage : allCoverages){
					out << coverage.region_.genBedRecordCore().toDelimStr() ;
					out << "\t" << coverage.coverage_;
					out << "\n";
				}
			}
		}
	};
	{
		njh::concurrent::runVoidFunctionThreaded(getCovForBam, pars.numThreads);
	}
}
std::vector<std::shared_ptr<Bed6RecordCore>> RunRegionRefinement(const RegionRefinementPars & pars){
	concurrent::BamReaderPool bamPool(pars.bamFnp, pars.numThreads);
	bamPool.openBamFile();
	OutputStream out(pars.outOpts);
	auto beds = getBeds(pars.bedFnp);
	njh::concurrent::LockableQueue<std::shared_ptr<Bed6RecordCore>> bedQueue(beds);
	std::function<void()> refineRegions =[&bamPool,&bedQueue,&pars](){
		std::shared_ptr<Bed6RecordCore> region;
		BamTools::BamAlignment bAln;
		while(bedQueue.getVal(region)){
			auto bamReader = bamPool.popReader();
			setBamFileRegionThrow(*bamReader, *region);
			uint32_t minStart = std::numeric_limits<uint32_t>::max();
			uint32_t maxEnd = 0;
			uint32_t revCount = 0;
			uint32_t forCount = 0;
			while(bamReader->GetNextAlignmentCore(bAln)){
				if(bAln.Position < static_cast<int64_t>(minStart)){
					minStart = bAln.Position;
				}
				if(bAln.GetEndPosition() > static_cast<int64_t>(maxEnd)){
					maxEnd = bAln.GetEndPosition();
				}
				if(!bAln.IsPaired()){
					if(bAln.IsReverseStrand() ){
						++revCount;
					}
					if(!bAln.IsReverseStrand() ){
						++forCount;
					}
				}
			}
			bool changeScore = region->score_ == region->length();

			if(minStart > region->chromStart_){
				region->chromStart_ = minStart;
			}
			if(maxEnd < region->chromEnd_){
				region->chromEnd_ = maxEnd;
			}
			if(changeScore){
				region->score_ = region->length();
			}
			if(pars.reOrient){
				if((revCount + forCount >0 ) && static_cast<double>(revCount)/(forCount + revCount) > 0.50){
					region->strand_ = '-';
				}else{
					region->strand_ = '+';
				}
			}
		}
	};

	njh::concurrent::runVoidFunctionThreaded(refineRegions, pars.numThreads);

	for(auto & bed : beds){
		bed->name_ = bed->genUIDFromCoordsWithStrand();
		out << bed->toDelimStrWithExtra() << std::endl;
	}
	return beds;
}


BamCountSpecficRegionsPars::BamCountSpecficRegionsPars(){
	pairPars.minOverlap_ = 8;
}

void BamCountSpecficRegionsPars::setDefaults(seqSetUp & setUp){
	setUp.setOption(forcePlusStrand, "--forcePlusStrand", "force Plus Strand orientation, otherwise use the strand orientation in bed file");
	setUp.setOption(baseQuality, "--baseQuality", "Base Quality Cut Off");
	setUp.setOption(mappingQuality, "--mappingQuality", "Mapping Quality Cut Off");
	setUp.setOption(matchIDCutOff, "--matchIDCutOff", "Match ID Cut Off");
	setUp.setOption(numThreads, "--numThreads", "Number Threads");
	setUp.setOption(pairPars.minOverlap_, "--minOverlap", "Min Overlap in overlap when stitching pairs");
	setUp.setOption(pairPars.errorAllowed_, "--errorAllowed", "Error allowed in overlap when stitching pairs");
	setUp.setOption(pairPars.hardMismatchCutOff_, "--hardMismatchCutOff", "Hard Mismatch Cut Off allowed in overlaop when stitching pairs");
	setUp.setOption(pairPars.lqMismatchCutOff, "--lqMismatchCutOff", "low quality Mismatch Cut Off allowed in overlaop when stitching pairs");
	setUp.setOption(pairPars.hqMismatchCutOff, "--hqMismatchCutOff", "hq Mismatch Cut Off allowed in overlaop when stitching pairs");
	setUp.setOption(extendAmount, "--extendAmount", "extend this amount around the region to get better trimming");
	setUp.setOption(countDuplicates, "--countDuplicates", "Skip reads marked as duplicates");
	setUp.setOption(totalCountCutOff, "--totalCountCutOff", "Total Count Cut Off");
	setUp.setOption(perBaseCountCutOff, "--perBaseCountCutOff", "Per Base Count Cut Off");

}

std::unordered_map<std::string, std::unordered_map<std::string, uint64_t>> BamCountSpecficRegions(
		const std::vector<GenomicRegion> & inputRegions,
		const std::unordered_map<std::string, seqInfo> & regionSeqs,
		const bfs::path bamFnp,
		const BamCountSpecficRegionsPars & pars){

	njh::concurrent::LockableQueue<GenomicRegion> regionsQueue(inputRegions);
	njhseq::concurrent::BamReaderPool bamPool(bamFnp, pars.numThreads);
	bamPool.openBamFile();
	std::map<std::string, uint32_t> spanningReadCounts;
	std::map<std::string, uint32_t> totalReadCounts;
	std::mutex spanRCountsMut;

	uint64_t maxlenForRegions =0 ;
	for(const auto & regSeq : regionSeqs){
		readVec::getMaxLength(regSeq.second, maxlenForRegions);
	}
	std::unordered_map<std::string, std::unordered_map<std::string, uint64_t>> allCounts;

	std::function<void()> extractReadsForRegion = [&regionsQueue,
																&bamPool,&pars,
																&spanningReadCounts,&totalReadCounts,
																&spanRCountsMut,
																&regionSeqs,
																&maxlenForRegions,
																&allCounts](){

		aligner alignerObj(std::max<uint64_t>(maxlenForRegions, 500), gapScoringParameters(10,1,0,0,0,0));
		GenomicRegion currentRegion;
		auto currentBReader = bamPool.popReader();
		auto refData = currentBReader->GetReferenceData();
		BamTools::BamAlignment bAln;
		std::unordered_map<std::string, uint32_t> currentSpanningReadCounts;
		std::unordered_map<std::string, uint32_t> currentTotalReadCounts;
		std::unordered_map<std::string, std::unordered_map<std::string, uint64_t>> currentCounts;

		uint64_t maxlen = 0;
		PairedReadProcessor pProcess(pars.pairPars);
		PairedReadProcessor::ProcessedResultsCounts processCounts;
		while(regionsQueue.getVal(currentRegion)){
			const auto & refSeq = regionSeqs.at(currentRegion.uid_);
			readVecTrimmer::GlobalAlnTrimPars trimPars;
			trimPars.needJustOneEnd_ = false;
			trimPars.startInclusive_ = pars.extendAmount;
			trimPars.endInclusive_ = len(refSeq) - 1 - pars.extendAmount;

			BamAlnsCache cache;
			setBamFileRegionThrow(*currentBReader, currentRegion);
			while(currentBReader->GetNextAlignment(bAln)){
				if(bAln.IsPrimaryAlignment() && bAln.IsMapped()){
					if(bAln.MapQuality < pars.mappingQuality){
						continue;
					}
					if(bAln.IsDuplicate() && !pars.countDuplicates){
						continue;
					}
					if(bAln.IsPaired() && bAln.IsMateMapped() && bAln.MatePosition < currentRegion.end_){
						++currentTotalReadCounts[currentRegion.uid_];
						if(bAln.Position < bAln.MatePosition){
							cache.add(bAln);
						}else {
							if(cache.has(bAln.Name)) {
								auto mate = cache.get(bAln.Name);
								//see if stitching is even plausible
								if(mate->GetEndPosition() >  bAln.Position){
									if(mate->Position <= currentRegion.start_ &&
										 bAln.GetEndPosition() >= currentRegion.end_ &&
										 mate->IsReverseStrand() != bAln.IsReverseStrand()){
										//stitching plausible and spanning
										//spanning read
										seqInfo firstMate;
										seqInfo secondMate;
										bool reverseCompFirstMate = false;
										if(bAln.IsFirstMate()){
											firstMate = bamAlnToSeqInfo(bAln, false);
											secondMate = bamAlnToSeqInfo(*mate, false);
											secondMate.reverseComplementRead(false, true);
											reverseCompFirstMate = bAln.IsReverseStrand();
										}else{
											firstMate = bamAlnToSeqInfo(*mate, false);
											secondMate = bamAlnToSeqInfo(bAln, false);
											secondMate.reverseComplementRead(false, true);
											reverseCompFirstMate = mate->IsReverseStrand();
										}
										readVec::getMaxLength(firstMate, maxlen);
										readVec::getMaxLength(secondMate, maxlen);
										PairedRead pseq(firstMate, secondMate);
										pseq.mateRComplemented_ = false;
										auto pairRes = pProcess.processPairedEnd(pseq, processCounts, alignerObj);
										if(nullptr != pairRes.combinedSeq_){
											readVec::getMaxLength(*pairRes.combinedSeq_, maxlen);
											alignerObj.parts_.setMaxSize(maxlen);
											if(reverseCompFirstMate != currentRegion.reverseSrand_){
												pairRes.combinedSeq_->reverseComplementRead(false, true);
											}
											seqInfo queryCopy = *pairRes.combinedSeq_;
											readVecTrimmer::trimSeqToRefByGlobalAln(*pairRes.combinedSeq_, refSeq, trimPars, alignerObj);
											alignerObj.profilePrimerAlignment(refSeq, queryCopy);
											if(pairRes.combinedSeq_->on_ && alignerObj.comp_.distances_.eventBasedIdentity_ > pars.matchIDCutOff){
												if(vectorMean(pairRes.combinedSeq_->qual_) > pars.baseQuality){
													//spanning read
													++currentSpanningReadCounts[currentRegion.uid_];
													++currentCounts[currentRegion.uid_][pairRes.combinedSeq_->seq_];
												}
											}
										}else if(bAln.Position <= currentRegion.start_ && bAln.GetEndPosition() >= currentRegion.end_){
											//spanning read
											seqInfo querySeq = bamAlnToSeqInfo(bAln, false);
											if(bAln.IsReverseStrand() != currentRegion.reverseSrand_){
												querySeq.reverseComplementRead(false, true);
											}
											seqInfo querySeqCopy = querySeq;

											readVec::getMaxLength(querySeq, maxlen);
											alignerObj.parts_.setMaxSize(maxlen);
											readVecTrimmer::trimSeqToRefByGlobalAln(querySeq, refSeq,trimPars, alignerObj);
											alignerObj.profilePrimerAlignment(refSeq, querySeqCopy);
											if(querySeq.on_ && alignerObj.comp_.distances_.eventBasedIdentity_ > pars.matchIDCutOff){
												if(vectorMean(querySeq.qual_) > pars.baseQuality){
													//spanning read
													++currentSpanningReadCounts[currentRegion.uid_];
													++currentCounts[currentRegion.uid_][querySeq.seq_];
												}
											}
										}else if(mate->Position <= currentRegion.start_ && mate->GetEndPosition() >= currentRegion.end_){
											//spanning read
											seqInfo querySeq = bamAlnToSeqInfo(*mate, false);
											if(mate->IsReverseStrand() != currentRegion.reverseSrand_){
												querySeq.reverseComplementRead(false, true);
											}
											seqInfo querySeqCopy = querySeq;

											readVec::getMaxLength(querySeq, maxlen);
											alignerObj.parts_.setMaxSize(maxlen);
											readVecTrimmer::trimSeqToRefByGlobalAln(querySeq, refSeq,trimPars, alignerObj);
											alignerObj.profilePrimerAlignment(refSeq, querySeqCopy);

											if(querySeq.on_ && alignerObj.comp_.distances_.eventBasedIdentity_ > pars.matchIDCutOff){
												if(vectorMean(querySeq.qual_) > pars.baseQuality){
													//spanning read
													++currentSpanningReadCounts[currentRegion.uid_];
													++currentCounts[currentRegion.uid_][querySeq.seq_];
												}
											}
										}
									}
								}
								cache.remove(bAln.Name);
							} else {
								//doesn't have mate, it's mate probably doesn't map to this region
								if(bAln.Position <= currentRegion.start_ && bAln.GetEndPosition() >= currentRegion.end_){
									//spanning read
									seqInfo querySeq = bamAlnToSeqInfo(bAln, false);
									if(bAln.IsReverseStrand() != currentRegion.reverseSrand_){
										querySeq.reverseComplementRead(false, true);
									}
									seqInfo querySeqCopy = querySeq;

									readVec::getMaxLength(querySeq, maxlen);
									alignerObj.parts_.setMaxSize(maxlen);
									readVecTrimmer::trimSeqToRefByGlobalAln(querySeq, refSeq,trimPars, alignerObj);
									alignerObj.profilePrimerAlignment(refSeq, querySeqCopy);

									if(querySeq.on_ && alignerObj.comp_.distances_.eventBasedIdentity_ > pars.matchIDCutOff){
										if(vectorMean(querySeq.qual_) > pars.baseQuality){
											//spanning read
											++currentSpanningReadCounts[currentRegion.uid_];
											++currentCounts[currentRegion.uid_][querySeq.seq_];
										}
									}
								}
							}
						}
					}else{
						++currentTotalReadCounts[currentRegion.uid_];
						if(bAln.Position <= currentRegion.start_ && bAln.GetEndPosition() >= currentRegion.end_){
							//spanning read
							seqInfo querySeq = bamAlnToSeqInfo(bAln, false);
							if(bAln.IsReverseStrand() != currentRegion.reverseSrand_){
								querySeq.reverseComplementRead(false, true);
							}
							seqInfo querySeqCopy = querySeq;
							readVec::getMaxLength(querySeq, maxlen);
							alignerObj.parts_.setMaxSize(maxlen);
							readVecTrimmer::trimSeqToRefByGlobalAln(querySeq, refSeq, trimPars, alignerObj);
							alignerObj.profilePrimerAlignment(refSeq, querySeqCopy);
							if(querySeq.on_ && alignerObj.comp_.distances_.eventBasedIdentity_ > pars.matchIDCutOff){
								if(vectorMean(querySeq.qual_) > pars.baseQuality){
									//spanning read
									++currentSpanningReadCounts[currentRegion.uid_];
									++currentCounts[currentRegion.uid_][querySeq.seq_];
								}
							}
						}
					}
				}
			}
		}
		{
			std::lock_guard<std::mutex> lock(spanRCountsMut);
			for(const auto & spanCount : currentSpanningReadCounts){
				spanningReadCounts[spanCount.first] = spanCount.second;
			}

			for(const auto & totalCount : currentTotalReadCounts){
				totalReadCounts[totalCount.first] = totalCount.second;
			}

			for(const auto & regionCount : currentCounts){
				allCounts[regionCount.first] = regionCount.second;
			}

		}
	};

	njh::concurrent::runVoidFunctionThreaded(extractReadsForRegion, pars.numThreads);

	//zero out counts;
	for(const auto & region : inputRegions){
		spanningReadCounts[region.uid_] += 0;
	}

	std::unordered_map<std::string, std::unordered_map<std::string, uint64_t>> ret;
	for(const auto & count : allCounts){
		uint64_t total = 0;
		for(const auto & subCount : count.second){
			if(subCount.second > pars.perBaseCountCutOff){
				total += subCount.second;
			}
		}
		if(total > pars.totalCountCutOff){
			for(const auto & subCount : count.second){
				if(subCount.second > pars.perBaseCountCutOff){
					ret[count.first][subCount.first] = subCount.second;
				}
			}
		}
	}
	return ret;
}
}  // namespace njhseq
