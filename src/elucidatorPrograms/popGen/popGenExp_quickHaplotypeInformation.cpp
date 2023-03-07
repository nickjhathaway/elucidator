/*
 * popGenExp_quickHaplotypeInformation.cpp
 *
 *  Created on: Jun 23, 2021
 *      Author: nick
 */
#include "popGenExp.hpp"
#include "elucidator/objects/BioDataObject.h"
#include <njhseq/PopulationGenetics.h>
#include <njhseq/objects/seqContainers/CollapsedHaps.hpp>
#include <njhseq/concurrency/PairwisePairFactory.hpp>

namespace njhseq {

int popGenExpRunner::quickHaplotypeInformation(const njh::progutils::CmdArgs & inputCommands) {
	std::string identifier = "region";
	uint32_t numThreads = 1;
  CollapsedHaps::GenPopMeasuresPar calcPopMeasuresPars;
	bool noDiagAlnPairwiseComps = false;


	seqSetUp setUp(inputCommands);
	setUp.pars_.ioOptions_.out_.outFilename_ = "";
	setUp.pars_.ioOptions_.out_.outFileFormat_ = "tab";
	setUp.pars_.ioOptions_.out_.outExtention_ = ".tab.txt";
	setUp.processReadInNames(true);
	setUp.setOption(numThreads, "--numThreads", "number of threads");
	calcPopMeasuresPars.numThreads = numThreads;
	setUp.setOption(noDiagAlnPairwiseComps, "--noDiagAlnPairwiseComps", "Use diagonal Alignment for Pairwise Comparisons");
	calcPopMeasuresPars.diagAlnPairwiseComps = !noDiagAlnPairwiseComps;

	setUp.setOption(calcPopMeasuresPars.numSegSites_, "--numSegSites", "If known, the number of segrating sites than tajima's d can be calculated");
	setUp.setOption(calcPopMeasuresPars.lowVarFreq, "--lowVariantCutOff", "Used for length variation, if more than this frequency differ from the most common length than it reports as length variation");
	setUp.setOption(calcPopMeasuresPars.getPairwiseComps, "--getPairwiseComps", "Get Pairwise Comps");
	setUp.setOption(identifier, "--identifier", "Give a identifier name for info");
	setUp.processWritingOptions(setUp.pars_.ioOptions_.out_);
	setUp.processVerbose();
	setUp.processAlignerDefualts();
	setUp.finishSetUp(std::cout);
	bool writeHeader = true;
	if(bfs::exists(setUp.pars_.ioOptions_.out_.outName()) && setUp.pars_.ioOptions_.out_.append_){
		writeHeader = false;
	}
	auto inputSeqs = CollapsedHaps::readInReads(setUp.pars_.ioOptions_);
	inputSeqs.verbose_ = setUp.pars_.verbose_;
	OutputStream out(setUp.pars_.ioOptions_.out_);
	std::shared_ptr<aligner> alignerObj;
	if(calcPopMeasuresPars.getPairwiseComps){
		uint64_t maxLen = readVec::getMaxLength(inputSeqs.seqs_);
		alignerObj = std::make_shared<aligner>(maxLen, setUp.pars_.gapInfo_, setUp.pars_.scoring_);
		alignerObj->processAlnInfoInput(setUp.pars_.alnInfoDirName_, setUp.pars_.verbose_);
	}

	auto divMeasures = inputSeqs.getGeneralMeasuresOfDiversity(calcPopMeasuresPars, alignerObj);
	if(writeHeader){
		out << njh::conToStr(calcPopMeasuresPars.genHeader(), "\t") << std::endl;
	}
	out << njh::conToStr(divMeasures.getOut(inputSeqs, identifier, calcPopMeasuresPars), "\t") << std::endl;
	if(calcPopMeasuresPars.getPairwiseComps){
		alignerObj->processAlnInfoOutput(setUp.pars_.outAlnInfoDirName_, setUp.pars_.verbose_);
	}
	return 0;
}


}  // namespace njhseq

