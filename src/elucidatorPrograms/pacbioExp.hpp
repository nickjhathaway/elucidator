#pragma once
/*
 * pacbioExpRunner.h
 *
 *  Created on: Feb 4, 2015
 *      Author: nickhathaway
 */

#include <njhcpp/progutils.h>
#include <njhseq.h>


namespace njhseq {

class pacbioExpRunner : public njh::progutils::ProgramRunner {
 public:
	pacbioExpRunner();
	static int profilePacbioReads(const njh::progutils::CmdArgs & inputCommands);
	static int createCssFromRawPacbio(const njh::progutils::CmdArgs & inputCommands);
	static int testCssConsensusBuilding(const njh::progutils::CmdArgs & inputCommands);
	static int checkCssReadsToRaw(const njh::progutils::CmdArgs & inputCommands);
	static int indexErrorsRawToCcsProcessed(const njh::progutils::CmdArgs & inputCommands);
	static int indexErrorsRawToCcs(const njh::progutils::CmdArgs & inputCommands);
	static int processSnpTable(const njh::progutils::CmdArgs & inputCommands);

	static int testVariantCalling(const njh::progutils::CmdArgs & inputCommands);

	static int printPacbioIDs(const njh::progutils::CmdArgs & inputCommands);
	static int getRoundsForPacbioIDs(const njh::progutils::CmdArgs & inputCommands);

	//take each read and simulate like it made it through pacbio ccs creation
	static int simPacbioPerRead(const njh::progutils::CmdArgs & inputCommands);

};
} /* namespace njhseq */


