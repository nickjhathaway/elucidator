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

struct BamCountSpecficRegionsPars{
	BamCountSpecficRegionsPars(){
		pairPars.minOverlap_ = 8;
	}
	uint32_t extendAmount = 30;
	uint32_t numThreads = 1;
	uint32_t mappingQuality = 20;
	uint32_t baseQuality = 25;
	double matchIDCutOff = 0.70;

	uint32_t totalCountCutOff = 5;
	uint32_t perBaseCountCutOff = 3;
	bool forcePlusStrand = false;
	bool countDuplicates = false;
	PairedReadProcessor::ProcessParams pairPars;

	void setDefaults(seqSetUp & setUp){
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
};

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
	std::unordered_map<std::string, std::unordered_map<std::string, uint64_t>> counts;

	std::function<void()> extractReadsForRegion = [&regionsQueue,
																&bamPool,&pars,
																&spanningReadCounts,&totalReadCounts,
																&spanRCountsMut,
																&regionSeqs,
																&maxlenForRegions,
																&counts](){

		aligner alignerObj(std::max<uint64_t>(maxlenForRegions, 500), gapScoringParameters(10,1,0,0,0,0));
		GenomicRegion currentRegion;
		auto currentBReader = bamPool.popReader();
		auto refData = currentBReader->GetReferenceData();
		BamTools::BamAlignment bAln;
		std::unordered_map<std::string, uint32_t> currentSpanningReadCounts;
		std::unordered_map<std::string, uint32_t> currentTotalReadCounts;

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
													++counts[currentRegion.uid_][pairRes.combinedSeq_->seq_];
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
													++counts[currentRegion.uid_][querySeq.seq_];
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
													++counts[currentRegion.uid_][querySeq.seq_];
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
											++counts[currentRegion.uid_][querySeq.seq_];
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
									++counts[currentRegion.uid_][querySeq.seq_];
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
		}
	};

	njh::concurrent::runVoidFunctionThreaded(extractReadsForRegion, pars.numThreads);

	//zero out counts;
	for(const auto & region : inputRegions){
		spanningReadCounts[region.uid_] += 0;
	}

	std::unordered_map<std::string, std::unordered_map<std::string, uint64_t>> ret;
	for(const auto & count : counts){
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



int bamExpRunner::MultipleBamGetPileupForRegion(
		const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bedFnp = "";
	bfs::path twoBitFnp = "";
	std::string bams = "";
	std::string pat = ".*.bam$";
	bfs::path dir = "./";


	BamCountSpecficRegionsPars countPars;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(bedFnp, "--bed", "Bed file", true);
	setUp.setOption(twoBitFnp, "--twoBitFnp", "Two Bit file to genome that sequences were aligned to", true);
	setUp.setOption(dir, "--dir", "Directory to search in");
	setUp.setOption(pat, "--pat", "Pattern in current directory to get coverage for");
	setUp.setOption(bams, "--bams", "Either a file with the name of a bam file on each line or a comma separated value of bam file paths");
	setUp.processDirectoryOutputName(njh::pasteAsStr("MultipleBamGetPileupForRegion_", basename(bedFnp), "_", "TODAY"), true);
	countPars.setDefaults(setUp);
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
	if(countPars.forcePlusStrand){
		for(auto & region : inputRegions){
			region.reverseSrand_ = false;
		}
	}

	TwoBit::TwoBitFile topRReader(twoBitFnp);
	std::unordered_map<std::string, seqInfo> regionSeqs;
	std::unordered_map<std::string, seqInfo> inputRegionSeqs;

	uint64_t maxlenForRegions = 0;
	auto chromLengths = topRReader.getSeqLens();
	for(const auto & reg : inputRegions){
		auto regCopy = reg;
		BedUtility::extendLeftRight(regCopy, countPars.extendAmount, countPars.extendAmount, njh::mapAt(chromLengths, regCopy.chrom_));
		regionSeqs[reg.uid_] = regCopy.extractSeq(topRReader);
		readVec::getMaxLength(regionSeqs[reg.uid_], maxlenForRegions);
		inputRegionSeqs[reg.uid_] = reg.extractSeq(topRReader);
	}
	auto bamFnps = njh::files::gatherFilesByPatOrNames(dir, std::regex{pat}, bams);
	checkBamFilesForIndexesAndAbilityToOpen(bamFnps);
	OutputStream outCounts(njh::files::make_path(setUp.pars_.directoryName_, "seqCounts.tab.txt.gz"));
	outCounts << "bamFile\tregion\trefSeq\tseq\tcount" << std::endl;

	for(const auto & bam : bamFnps){
		auto counts = BamCountSpecficRegions(inputRegions, regionSeqs, bam, countPars);
		for(const auto & region : inputRegions){
			std::set<std::string> subCounts;
			njh::addVecToSet(getVectorOfMapKeys(counts[region.uid_]), subCounts);
			subCounts.emplace(inputRegionSeqs[region.uid_].seq_);
			uint64_t total = 0;
			for(const auto & seq : subCounts){
				total += counts[region.uid_][seq];
			}

			for(const auto & seq : subCounts){
				outCounts << bam.filename().string()
						<< "\t" << region.uid_
						<< "\t" << inputRegionSeqs[region.uid_].seq_
						<< "\t" << seq
						<< "\t" << counts[region.uid_][seq]
						<< std::endl;
			}
		}
	}


	return 0;
}


int bamExpRunner::BamGetPileupForRegion(
		const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bedFnp = "";
	bfs::path twoBitFnp = "";
	BamCountSpecficRegionsPars countPars;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processReadInNames(VecStr{"--bam"});
	setUp.processDirectoryOutputName(true);
	setUp.setOption(bedFnp, "--bed", "Bed file", true);
	setUp.setOption(twoBitFnp, "--twoBitFnp", "Two Bit file to genome that sequences were aligned to", true);
	countPars.setDefaults(setUp);
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
	if(countPars.forcePlusStrand){
		for(auto & region : inputRegions){
			region.reverseSrand_ = false;
		}
	}

	TwoBit::TwoBitFile topRReader(twoBitFnp);
	std::unordered_map<std::string, seqInfo> regionSeqs;
	std::unordered_map<std::string, seqInfo> inputRegionSeqs;

	uint64_t maxlenForRegions = 0;
	auto chromLengths = topRReader.getSeqLens();
	for(const auto & reg : inputRegions){
		auto regCopy = reg;
		BedUtility::extendLeftRight(regCopy, countPars.extendAmount, countPars.extendAmount, njh::mapAt(chromLengths, regCopy.chrom_));
		regionSeqs[reg.uid_] = regCopy.extractSeq(topRReader);
		readVec::getMaxLength(regionSeqs[reg.uid_], maxlenForRegions);
		inputRegionSeqs[reg.uid_] = reg.extractSeq(topRReader);
	}

	auto counts = BamCountSpecficRegions(inputRegions, regionSeqs, setUp.pars_.ioOptions_.firstName_, countPars);
	OutputStream outCounts(njh::files::make_path(setUp.pars_.directoryName_, "seqCounts.tab.txt.gz"));
	outCounts << "region\trefSeq\tseq\tcount" << std::endl;
	for(const auto & region : inputRegions){
		std::set<std::string> subCounts;
		njh::addVecToSet(getVectorOfMapKeys(counts[region.uid_]), subCounts);
		subCounts.emplace(inputRegionSeqs[region.uid_].seq_);
		uint64_t total = 0;
		for(const auto & seq : subCounts){
			total += counts[region.uid_][seq];
		}
		for(const auto & seq : subCounts){
			outCounts << region.uid_
					<< "\t" << inputRegionSeqs[region.uid_].seq_
					<< "\t" << seq
					<< "\t" << counts[region.uid_][seq]
					<< std::endl;
		}
	}
	return 0;
}


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
											readVecTrimmer::trimSeqToRefByGlobalAln(*pairRes.combinedSeq_, refSeq,trimPars, alignerObj);
											if(pairRes.combinedSeq_->on_){
												//spanning read
												++currentSpanningReadCounts[currentRegion.uid_];
												seqIOOutCache.add(currentRegion.uid_, *pairRes.combinedSeq_);
											}
										}else if(bAln.Position <= currentRegion.start_ && bAln.GetEndPosition() >= currentRegion.end_){
											//spanning read
											seqInfo querySeq = bamAlnToSeqInfo(bAln, false);
											if(bAln.IsReverseStrand() != currentRegion.reverseSrand_){
												querySeq.reverseComplementRead(false, true);
											}
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
											if(mate->IsReverseStrand() != currentRegion.reverseSrand_){
												querySeq.reverseComplementRead(false, true);
											}
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
									if(bAln.IsReverseStrand() != currentRegion.reverseSrand_){
										querySeq.reverseComplementRead(false, true);
									}
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
							if(bAln.IsReverseStrand() != currentRegion.reverseSrand_){
								querySeq.reverseComplementRead(false, true);
							}
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

