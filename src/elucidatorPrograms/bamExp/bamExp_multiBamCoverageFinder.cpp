/*
 * bapExp_multiBamCoverageFinder.cpp
 *
 *  Created on: May 22, 2020
 *      Author: nick
 */


#include "bamExp.hpp"

#include "elucidator/BamToolsUtils.h"
#include "elucidator/objects/BioDataObject.h"


namespace njhseq {


int bamExpRunner::multiBamCoverageFinder(const njh::progutils::CmdArgs & inputCommands){
	RunCoverageFinderMultiPars covPars;

	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(covPars.byBases, "--byBases", "Output coverage in bases rather than reads");
	setUp.setOption(covPars.regionBatchSize, "--regionBatchSize", "Region Batch Size");
	setUp.setOption(covPars.coverageCutOff, "--coverageCutOff", "Coverage Cut Off to include region in output");
	setUp.setOption(covPars.window, "--window", "window size for sliding window");
	setUp.setOption(covPars.step, "--step", "step size for advancing the window");
	setUp.setOption(covPars.numThreads, "--numThreads", "Number of threads to use");
	setUp.processWritingOptions(covPars.outOpts);
	setUp.setOption(covPars.pat, "--pat", "Pattern in current directory to get coverage for");
	setUp.setOption(covPars.bams, "--bams", "Either a file with the name of a bam file on each line or a comma separated value of bam file paths");
	setUp.setOption(covPars.chromsToSkip, "--chromsToSkip", "Skip these chromosomes");

	setUp.finishSetUp(std::cout);


	RunCoverageFinderMulti(covPars);

	return 0;
}

int bamExpRunner::multiBamCoverageFinderBases(const njh::progutils::CmdArgs & inputCommands){
	RunCoverageFinderMultiPars covPars;
	covPars.byBases = true;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(covPars.regionBatchSize, "--regionBatchSize", "Region Batch Size");
	setUp.setOption(covPars.coverageCutOff, "--coverageCutOff", "Coverage Cut Off to include region in output");
	setUp.setOption(covPars.window, "--window", "window size for sliding window");
	setUp.setOption(covPars.step, "--step", "step size for advancing the window");
	setUp.setOption(covPars.numThreads, "--numThreads", "Number of threads to use");
	setUp.processWritingOptions(covPars.outOpts);
	setUp.setOption(covPars.pat, "--pat", "Pattern in current directory to get coverage for");
	setUp.setOption(covPars.bams, "--bams", "Either a file with the name of a bam file on each line or a comma separated value of bam file paths");
	setUp.setOption(covPars.chromsToSkip, "--chromsToSkip", "Skip these chromosomes");

	setUp.finishSetUp(std::cout);



	RunCoverageFinderMulti(covPars);

	return 0;
}


}  // namespace njhseq



