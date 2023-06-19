/*
 * genExp_profileErrorsReadsToReferenceSeq.cpp
 *
 *  Created on: Jul 12, 2020
 *      Author: nicholashathaway
 */




#include "genExp.hpp"
#include "elucidator/simulation/SeqTechSimulation/Illumina/RoughIlluminaProfiler.hpp"


namespace njhseq {


int genExpRunner::profileErrorsReadsToReferenceSeq(const njh::progutils::CmdArgs & inputCommands){

	seqInfo refInfo("ref");

	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
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


	alignerObj.processAlnInfoOutputNoCheck(setUp.pars_.outAlnInfoDirName_, setUp.pars_.verbose_);

	return 0;

}



}  // namespace njhseq


