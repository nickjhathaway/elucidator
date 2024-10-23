/*
 * genExp_profileErrorsReadsToReferenceSeq.cpp
 *
 *  Created on: Jul 12, 2020
 *      Author: nicholashathaway
 */


#include <njhseq/seqToolsUtils/tandemRepeatUtils.hpp>

#include "genExp.hpp"
#include "elucidator/simulation/SeqTechSimulation/Illumina/RoughIlluminaProfiler.hpp"


namespace njhseq {


int genExpRunner::profileErrorsReadsToReferenceSeq(const njh::progutils::CmdArgs & inputCommands) {
	seqInfo refInfo("ref");
	bool countHomopolymerErrorRates = false;
	uint32_t minHomopolymerSize = 4;
	bool countSimpleTandemRepeatSize = false;
	uint32_t minSimpleTandemRepeatSize = 2;
	uint32_t maxSimpleTandemRepeatSize = 3;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();

	setUp.setOption(countHomopolymerErrorRates, "--countHomopolymerErrorRates", "count Homopolymer Error Rates");
	setUp.setOption(minHomopolymerSize, "--minHomopolymerSize", "min Homopolymer Size");
	setUp.setOption(countSimpleTandemRepeatSize, "--countSimpleTandemRepeatSize", "count Simple Tandem Repeat Size");
	setUp.setOption(minSimpleTandemRepeatSize, "--minSimpleTandemRepeatSize", "min Simple Tandem Repeat Size");
	setUp.setOption(maxSimpleTandemRepeatSize, "--maxSimpleTandemRepeatSize", "max Simple Tandem Repeat Size");


	setUp.processSeq(refInfo, "--ref", "Reference sequence to compare to", true);
	setUp.processReadInNames();
	setUp.processDirectoryOutputName(true);
	setUp.processAlnInfoInput();
	setUp.finishSetUp(std::cout);

	setUp.startARunLog(setUp.pars_.directoryName_);


	uint64_t maxLen = 0;
	uint32_t total = 0;
	{
		seqInfo seq;

		SeqInput reader(setUp.pars_.ioOptions_);
		reader.openIn();
		while(reader.readNextRead(seq)){
			++total;
			readVec::getMaxLength(seq, maxLen);
		}
	}
	readVec::getMaxLength(refInfo, maxLen);

	aligner alignerObj(maxLen, gapScoringParameters(5,1,0,0,0,0));
	alignerObj.processAlnInfoInputNoCheck(setUp.pars_.alnInfoDirName_, setUp.pars_.verbose_);

	RoughIlluminaProfiler::Counts errorCounter;
	{
		njh::ProgressBar pBar(total);
		seqInfo seq;
		SeqInput reader(setUp.pars_.ioOptions_);
		reader.openIn();
		while(reader.readNextRead(seq)){
			if(setUp.pars_.verbose_){
				pBar.outputProgAdd(std::cout, 1, true);
			}
			alignerObj.alignCacheGlobal(refInfo, seq);
			alignerObj.profileAlignment(refInfo, seq, false, true, false);
			errorCounter.increaseCounts(alignerObj.alignObjectA_.seqBase_, alignerObj.alignObjectB_.seqBase_, alignerObj.comp_);
		}
	}
	errorCounter.writeProfiles(njh::files::make_path(setUp.pars_.directoryName_, "all").string(), true);
	errorCounter.writeIndels(njh::files::make_path(setUp.pars_.directoryName_, "all").string(), true);

	if(countHomopolymerErrorRates || countSimpleTandemRepeatSize) {
		SimpleTandemRepeatFinder::SimpleTRFinderLocsPars pars;
		if(countHomopolymerErrorRates) {
			pars.minRepeatUnitSize = 1;
		} else {
			pars.minRepeatUnitSize = minSimpleTandemRepeatSize;
			pars.maxRepeatUnitSize = maxSimpleTandemRepeatSize;
		}

	}

	alignerObj.processAlnInfoOutputNoCheck(setUp.pars_.outAlnInfoDirName_, setUp.pars_.verbose_);

	return 0;

}



}  // namespace njhseq


