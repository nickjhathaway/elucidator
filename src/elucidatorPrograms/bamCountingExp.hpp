#pragma once
/*
 * bamCountingExpRunner.h
 *
 *  Created on: Feb 4, 2015
 *      Author: nickhathaway
 */

#include <njhcpp/progutils.h>
#include <njhseq.h>

namespace njhseq {

class bamCountingExpRunner : public njh::progutils::ProgramRunner {
 public:
	bamCountingExpRunner();

	static int bamBaseCountingSlim(const njh::progutils::CmdArgs & inputCommands);

	static int slimCounterBinToTxt(const njh::progutils::CmdArgs & inputCommands);
	static int readGzIndex(const njh::progutils::CmdArgs & inputCommands);
	static int coverageInfo(const njh::progutils::CmdArgs & inputCommands);

	static int gatherVariantsOnChrom(const njh::progutils::CmdArgs & inputCommands);
	static int gatherVariantsOnSample(const njh::progutils::CmdArgs & inputCommands);
	static int gatherVariantsOnSampleOnLoc(const njh::progutils::CmdArgs & inputCommands);
	static int gatherBaseCoverageOnSampleOnLoc(const njh::progutils::CmdArgs & inputCommands);

	static int gatherLocOnSample(const njh::progutils::CmdArgs & inputCommands);


};

} /* namespace njhseq */


