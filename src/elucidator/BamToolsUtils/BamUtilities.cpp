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


namespace njhseq {


void RunCoverageFinder(const CoverageFinderPars & pars){
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

	auto refData = bReaderForRefData.GetReferenceData();
	auto getCov = [&slider,&bamFnps,&outMut,&out,&refData, &pars](){
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
					if(coverage.coverage_ >= pars.coverageCutOff){
						coverages.emplace_back(coverage);
					}else{
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
		std::vector<std::thread> threads;
		for(uint32_t t = 0; t < pars.numThreads; ++t){
			threads.emplace_back(std::thread(getCov));
		}
		for(auto & t : threads){
			t.join();
		}
	}
}
void RunRegionRefinement(const RegionRefinementPars & pars){
	concurrent::BamReaderPool bamPool(pars.bamFnp, pars.numThreads);
	bamPool.openBamFile();
	OutputStream out(pars.outOpts);
	auto beds = getBeds(pars.bedFnp);
	njh::concurrent::LockableQueue<std::shared_ptr<Bed6RecordCore>> bedQueue(beds);
	auto refineRegions =[&bamPool,&bedQueue,&pars](){
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

	std::vector<std::thread> threads;
	for(uint32_t t = 0; t < pars.numThreads; ++t){
		threads.emplace_back(std::thread(refineRegions));
	}
	njh::concurrent::joinAllThreads(threads);
	for(const auto & bed : beds){
		out << bed->toDelimStrWithExtra() << std::endl;
	}
}

}  // namespace njhseq
