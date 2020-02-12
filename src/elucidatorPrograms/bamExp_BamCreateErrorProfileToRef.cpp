/*
 * bamExp_BamCreateErrorProfileToRef.cpp
 *
 *  Created on: Aug 27, 2019
 *      Author: nicholashathaway
 */



#include "bamExp.hpp"
#include "elucidator/simulation.h"


namespace njhseq {


int bamExpRunner::BamCreateErrorProfileToRef(
		const njh::progutils::CmdArgs & inputCommands) {
	ReAlignedSeq::genRealignmentPars ralnPars;

	bfs::path bedFnp = "";
	bfs::path twoBitFnp = "";
	uint32_t numThreads = 1;
	uint32_t lengthLimit = 0;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processReadInNames(VecStr{"--bam"});
	setUp.processDirectoryOutputName(true);
	setUp.setOption(lengthLimit, "--lengthLimit", "Length Limit");
	setUp.setOption(numThreads, "--numThreads", "Number Threads");
	setUp.setOption(bedFnp, "--bed", "Bed file");
	setUp.setOption(twoBitFnp, "--twoBitFnp", "2bit file for genome aligned to", true);
	setUp.setOption(ralnPars.extendAmount, "--extendAmount", "Extend Amount");

	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	BamTools::BamReader bReader;
	bReader.Open(setUp.pars_.ioOptions_.firstName_.string());
	checkBamOpenThrow(bReader, setUp.pars_.ioOptions_.firstName_);
	loadBamIndexThrow(bReader);
	auto refs = bReader.GetReferenceData();
	BamAlnsCache alnCache;
	auto refData = bReader.GetReferenceData();
	std::unordered_map<std::string, uint32_t> refNameToId;
	for (auto pos : iter::range(refData.size())) {
		refNameToId[refData[pos].RefName] = pos;
	}

	TwoBit::TwoBitFile tReader(twoBitFnp);
	auto chromLengths = tReader.getSeqLens();
	RoughIlluminaProfiler profiler;

	if("" != bedFnp){
		auto regions = bedPtrsToGenomicRegs(getBeds(bedFnp));
		njh::concurrent::LockableQueue<GenomicRegion> regionsQueue(regions);
		concurrent::BamReaderPool bamPools(setUp.pars_.ioOptions_.firstName_, numThreads);
		bamPools.openBamFile();
		std::mutex mut;
		std::function<void()> increaseCount = [&regionsQueue,&bamPools,&mut,&twoBitFnp,&refData,&setUp,&lengthLimit,&chromLengths,&ralnPars,&profiler](){
			aligner alignerObj(700, gapScoringParameters(5,1,0,0,0,0));
			GenomicRegion region;
			TwoBit::TwoBitFile tReader(twoBitFnp);
			auto bReader = bamPools.popReader();
			BamTools::BamAlignment bAln;
			RoughIlluminaProfiler currentProfiler;

			while(regionsQueue.getVal(region)){
				if(setUp.pars_.verbose_){
					std::lock_guard<std::mutex> lock(mut);
					std::cout << region.uid_ << std::endl;
				}
				setBamFileRegionThrow(*bReader, region);
				while (bReader->GetNextAlignment(bAln)) {
					if(bAln.IsMapped() && bAln.IsPrimaryAlignment()){
						if(bAln.GetEndPosition() - bAln.Position > lengthLimit){
							auto realignment = ReAlignedSeq::genRealignment(bAln, refData, alignerObj, chromLengths, tReader, ralnPars);
							currentProfiler.increaseCounts(realignment);
						}
					}
				}
			}
			{
				std::lock_guard<std::mutex> lock(mut);
				profiler.addOther(currentProfiler);
			}
		};
		njh::concurrent::runVoidFunctionThreaded(increaseCount, numThreads);
	}else{
		aligner alignerObj(700, gapScoringParameters(5,1,0,0,0,0));
		BamTools::BamAlignment bAln;
		while (bReader.GetNextAlignment(bAln)) {
			if(bAln.IsMapped() && bAln.IsPrimaryAlignment()){
				if(bAln.GetEndPosition() - bAln.Position > lengthLimit){
					auto realignment = ReAlignedSeq::genRealignment(bAln, refData, alignerObj, chromLengths, tReader, ralnPars);
					profiler.increaseCounts(realignment);
				}
			}
		}
	}
	profiler.r1_counts.writeProfiles(njh::files::make_path(setUp.pars_.directoryName_, "r1").string(), false);
	profiler.r1_counts.writeIndels(njh::files::make_path(setUp.pars_.directoryName_, "r1").string(), false);

	profiler.r2_counts.writeProfiles(njh::files::make_path(setUp.pars_.directoryName_, "r2").string(), false);
	profiler.r2_counts.writeIndels(njh::files::make_path(setUp.pars_.directoryName_, "r2").string(), false);

	return 0;
}

}  // namespace njhseq
