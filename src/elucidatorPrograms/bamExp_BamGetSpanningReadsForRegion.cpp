/*
 * bamExp_BamGetSpanningReadsForRegion.cpp
 *
 *  Created on: Aug 28, 2019
 *      Author: nicholashathaway
 */







#include "bamExp.hpp"
#include "elucidator/simulation.h"
#include <SeekDeep/objects/PairedReadProcessor.hpp>

namespace njhseq {


int bamExpRunner::BamGetSpanningReadsForRegion(
		const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bedFnp = "";
	bfs::path twoBitFnp = "";
	uint32_t extendAmount = 25;
	uint32_t numThreads = 1;
	PairedReadProcessor::ProcessParams pairPars;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processReadInNames(VecStr{"--bam"});
	setUp.processDirectoryOutputName(true);
	setUp.setOption(numThreads, "--numThreads", "Number Threads");
	setUp.setOption(bedFnp, "--bed", "Bed file", true);
	setUp.setOption(twoBitFnp, "--twoBitFnp", "Two Bit file to genome that sequences were aligned to", true);
	setUp.setOption(pairPars.minOverlap_, "--minOverlap", "Min Overlap in overlap when stitching pairs");
	setUp.setOption(pairPars.errorAllowed_, "--errorAllowed", "Error allowed in overlap when stitching pairs");
	setUp.setOption(pairPars.hardMismatchCutOff_, "--hardMismatchCutOff", "Hard Mismatch Cut Off allowed in overlaop when stitching pairs");
	setUp.setOption(pairPars.lqMismatchCutOff, "--lqMismatchCutOff", "low quality Mismatch Cut Off allowed in overlaop when stitching pairs");
	setUp.setOption(pairPars.hqMismatchCutOff, "--hqMismatchCutOff", "hq Mismatch Cut Off allowed in overlaop when stitching pairs");
	setUp.setOption(extendAmount, "--extendAmount", "extend this amount around the region to get better trimming");

	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);

	auto inputRegions = bedPtrsToGenomicRegs(getBeds(bedFnp));
	std::vector<std::string> regionNames;
	std::set<std::string> repeatNames;
	for(const auto & region : inputRegions){
		if(njh::in(region.uid_, regionNames)){
			repeatNames.emplace(region.uid_);
		}
		regionNames.emplace_back(region.uid_);
	}
	if(!repeatNames.empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error "
				<< "can't have regions with the same names, following names were repeated:"
				<< njh::conToStr(repeatNames) << "\n";
		throw std::runtime_error{ss.str()};
	}

	njh::concurrent::LockableQueue<GenomicRegion> regionsQueue(inputRegions);
	njhseq::concurrent::BamReaderPool bamPool(setUp.pars_.ioOptions_.firstName_, numThreads);
	bamPool.openBamFile();

	MultiSeqOutCache<seqInfo> seqIOOutCache;
	for(const auto & reg : inputRegions){
		seqIOOutCache.addReader(reg.uid_, SeqIOOptions::genFastqOutGz(njh::files::make_path(setUp.pars_.directoryName_, reg.uid_)));
	}
	std::map<std::string, uint32_t> spanningReadCounts;
	std::map<std::string, uint32_t> totalReadCounts;

	std::mutex spanRCountsMut;
	TwoBit::TwoBitFile topRReader(twoBitFnp);
	std::unordered_map<std::string, seqInfo> regionSeqs;
	uint64_t maxlenForRegions = 0;
	auto chromLengths = topRReader.getSeqLens();
	for(const auto & reg : inputRegions){
		auto regCopy = reg;
		BedUtility::extendLeftRight(regCopy, extendAmount, extendAmount, njh::mapAt(chromLengths, regCopy.chrom_));
		regionSeqs[reg.uid_] = regCopy.extractSeq(topRReader);
		readVec::getMaxLength(regionSeqs[reg.uid_], maxlenForRegions);
	}

	std::function<void()> extractReadsForRegion = [&seqIOOutCache,&regionsQueue,
																&bamPool,&pairPars,
																&spanningReadCounts,&totalReadCounts,
																&spanRCountsMut,
																&regionSeqs,&maxlenForRegions,
																&extendAmount](){
		aligner alignerObj(std::max<uint64_t>(maxlenForRegions, 500), gapScoringParameters(5,1,0,0,0,0));
		GenomicRegion currentRegion;
		auto currentBReader = bamPool.popReader();
		auto refData = currentBReader->GetReferenceData();
		BamTools::BamAlignment bAln;
		std::unordered_map<std::string, uint32_t> currentSpanningReadCounts;
		std::unordered_map<std::string, uint32_t> currentTotalReadCounts;

		uint64_t maxlen = 0;
		PairedReadProcessor pProcess(pairPars);
		PairedReadProcessor::ProcessedResultsCounts processCounts;
		while(regionsQueue.getVal(currentRegion)){

			const auto & refSeq = regionSeqs.at(currentRegion.uid_);

			readVecTrimmer::GlobalAlnTrimPars trimPars;
			trimPars.needJustOneEnd_ = false;
			trimPars.startInclusive_ = extendAmount;
			trimPars.endInclusive_ = len(refSeq) - 1 - extendAmount;

			BamAlnsCache cache;
			setBamFileRegionThrow(*currentBReader, currentRegion);
			while(currentBReader->GetNextAlignment(bAln)){
				if(bAln.IsPrimaryAlignment() && bAln.IsMapped()){
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
										if(bAln.IsFirstMate()){
											firstMate = bamAlnToSeqInfo(bAln, false);
											secondMate = bamAlnToSeqInfo(*mate, false);
											secondMate.reverseComplementRead(false, true);
										}else{
											firstMate = bamAlnToSeqInfo(*mate, false);
											secondMate = bamAlnToSeqInfo(bAln, false);
											secondMate.reverseComplementRead(false, true);
										}
										readVec::getMaxLength(firstMate, maxlen);
										readVec::getMaxLength(secondMate, maxlen);
										PairedRead pseq(firstMate, secondMate);
										pseq.mateRComplemented_ = false;
										auto pairRes = pProcess.processPairedEnd(pseq, processCounts, alignerObj);
										if(nullptr != pairRes.combinedSeq_){
											readVec::getMaxLength(*pairRes.combinedSeq_, maxlen);
											alignerObj.parts_.setMaxSize(maxlen);
											readVecTrimmer::trimSeqToRefByGlobalAln(*pairRes.combinedSeq_, refSeq,trimPars, alignerObj);
											if(pairRes.combinedSeq_->on_){
												//spanning read
												++currentSpanningReadCounts[currentRegion.uid_];
												seqIOOutCache.add(currentRegion.uid_, *pairRes.combinedSeq_);
											}
										}else if(bAln.Position <= currentRegion.start_ && bAln.GetEndPosition() >= currentRegion.end_){
											//spanning read
											seqInfo querySeq = bamAlnToSeqInfo(bAln, false);
											readVec::getMaxLength(querySeq, maxlen);
											alignerObj.parts_.setMaxSize(maxlen);
											readVecTrimmer::trimSeqToRefByGlobalAln(querySeq, refSeq,trimPars, alignerObj);
											if(querySeq.on_){
												//spanning read
												++currentSpanningReadCounts[currentRegion.uid_];
												seqIOOutCache.add(currentRegion.uid_, querySeq);
											}
										}else if(mate->Position <= currentRegion.start_ && mate->GetEndPosition() >= currentRegion.end_){
											//spanning read
											seqInfo querySeq = bamAlnToSeqInfo(*mate, false);
											readVec::getMaxLength(querySeq, maxlen);
											alignerObj.parts_.setMaxSize(maxlen);
											readVecTrimmer::trimSeqToRefByGlobalAln(querySeq, refSeq,trimPars, alignerObj);
											if(querySeq.on_){
												//spanning read
												++currentSpanningReadCounts[currentRegion.uid_];
												seqIOOutCache.add(currentRegion.uid_, querySeq);
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
									readVec::getMaxLength(querySeq, maxlen);
									alignerObj.parts_.setMaxSize(maxlen);
									readVecTrimmer::trimSeqToRefByGlobalAln(querySeq, refSeq,trimPars, alignerObj);
									if(querySeq.on_){
										//spanning read
										++currentSpanningReadCounts[currentRegion.uid_];
										seqIOOutCache.add(currentRegion.uid_, querySeq);
									}
								}
							}
						}
					}else{
						++currentTotalReadCounts[currentRegion.uid_];
						if(bAln.Position <= currentRegion.start_ && bAln.GetEndPosition() >= currentRegion.end_){
							//spanning read
							seqInfo querySeq = bamAlnToSeqInfo(bAln, false);
							readVec::getMaxLength(querySeq, maxlen);
							alignerObj.parts_.setMaxSize(maxlen);
							readVecTrimmer::trimSeqToRefByGlobalAln(querySeq, refSeq, trimPars, alignerObj);
							if(querySeq.on_){
								//spanning read
								++currentSpanningReadCounts[currentRegion.uid_];
								seqIOOutCache.add(currentRegion.uid_, querySeq);
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
		}
	};

	njh::concurrent::runVoidFunctionThreaded(extractReadsForRegion, numThreads);

	//zero out counts;
	for(const auto & region : inputRegions){
		spanningReadCounts[region.uid_] += 0;
	}
	OutputStream outCounts(njh::files::make_path(setUp.pars_.directoryName_, "spanningReadCounts.tab.txt"));
	outCounts << "region\tspanningReads\ttotalReads\tregionLength" << std::endl;
	for(const auto & spanCount : spanningReadCounts){
		outCounts << spanCount.first
				<< "\t" << spanCount.second
				<< "\t" << totalReadCounts[spanCount.first]
				<< "\t" << len(regionSeqs[spanCount.first]) - extendAmount * 2<< std::endl;
	}
	return 0;
}



}  // namespace njhseq
