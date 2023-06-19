/*
 * bamExp_getUnmappedAndShortAlnsFromBam.cpp
 *
 *  Created on: Feb 29, 2020
 *      Author: nicholashathaway
 */



#include "bamExp.hpp"

#include "elucidator/BamToolsUtils.h"
#include "elucidator/objects/BioDataObject.h"



namespace njhseq {


int bamExpRunner::getUnmappedAndShortAlnsFromBam(const njh::progutils::CmdArgs & inputCommands){
	BamExtractor::writeUnMappedSeqsAndSmallAlnsWithFiltersPars filtPars;
	uint32_t qualCheck = 30;
	double qualCheckCutOff = 0.50;
	bool noFilterForEntropy = false;

	bfs::path bedFnp = "";
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processDefaultReader({"--bam"}, true);
	setUp.setOption(noFilterForEntropy, "--noFilterForEntropy", "No Filter For Entropy");
	filtPars.filterOffLowEntropyShortAlnsRecruits_ = !noFilterForEntropy;
	setUp.setOption(filtPars.filterOffLowEntropyShortAlnsRecruitsCutOff_, "--filterOffLowEntropyShortAlnsRecruitsCutOff", "Filter Off Low Entropy Short Alns Recruits Cut Off");
	setUp.setOption(filtPars.maxAlnSize_, "--maxAlnSize", "Max Aln Size");
	setUp.setOption(filtPars.minQuerySize_, "--minQuerySize", "min Query Size");
	setUp.setOption(filtPars.minSotClip_, "--minSotClip", "min Sot Clip");
	setUp.setOption(filtPars.tryToFindOrphanMate_, "--tryToFindOrphanMate", "try To Find Orphan Mate");
	setUp.setOption(filtPars.entropyKlen_, "--entropyKlen", "K-mer length for entropy filter");


	setUp.setOption(qualCheck, "--qualCheck", "qual Check");
	setUp.setOption(qualCheckCutOff, "--qualCheckCutOff", "qual Check Cut Off");

	setUp.setOption(bedFnp, "--bed", "regions to skip");

	setUp.finishSetUp(std::cout);
	std::vector<GenomicRegion> regions;
	if("" != bedFnp){
		regions = bedPtrsToGenomicRegs(getBeds(bedFnp));
	}
	ReadCheckerQualCheck qualChecker(qualCheck, qualCheckCutOff, false);
	BamExtractor bExtractor(setUp.pars_.verbose_);
	auto outFiles = bExtractor.writeUnMappedSeqsAndSmallAlnsWithFilters(setUp.pars_.ioOptions_,qualChecker, filtPars, regions);

	return 0;
}



} // namespace njhseq
