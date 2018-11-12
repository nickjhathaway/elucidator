#pragma once
/*
 * programWrappers.hpp
 *
 *  Created on: Feb 2, 2017
 *      Author: nick
 */



#include <njhcpp/progutils.h>
#include <njhseq.h>

namespace njhseq {

class programWrapperRunner : public njh::progutils::ProgramRunner {
 public:
	programWrapperRunner();

	static int runTrinityOnRegion(const njh::progutils::CmdArgs & inputCommands);

	static int runSSAKEOnRegion(const njh::progutils::CmdArgs & inputCommands);

	static int runLastz(const njh::progutils::CmdArgs & inputCommands);


	static int sraIsPairedEnd(const njh::progutils::CmdArgs & inputCommands);

	static int sraFastqDump(const njh::progutils::CmdArgs & inputCommands);

	static int runTrimmomatic(const njh::progutils::CmdArgs & inputCommands);
	static int runBwaOnTrimmomaticOutputPE(const njh::progutils::CmdArgs & inputCommands);


	static int parsePrimer3OutputToJson(const njh::progutils::CmdArgs & inputCommands);
	static int parsePrimer3OutputToBed(const njh::progutils::CmdArgs & inputCommands);

	static int parsePrimer3OutputToPossibleMipArms(const njh::progutils::CmdArgs & inputCommands);


	static int findNonUniquePrimerArms(const njh::progutils::CmdArgs & inputCommands);

	static int convertShorahSupportSeqs(const njh::progutils::CmdArgs & inputCommands);

	static int runShorahAmplian(const njh::progutils::CmdArgs & inputCommands);

	static int runSamtoolsFlagStat(const njh::progutils::CmdArgs & inputCommands);

	static int runAdapterRemoval(const njh::progutils::CmdArgs & inputCommands);
	static int runAdapterRemovalSE(const njh::progutils::CmdArgs & inputCommands);
	static int setUpRunAdapterRemoval(const njh::progutils::CmdArgs & inputCommands);
	static int runBwaOnAdapterReomvalOutputSinglesCombined(const njh::progutils::CmdArgs & inputCommands);
	static int runBwaOnAdapterReomvalOutputSE(const njh::progutils::CmdArgs & inputCommands);


	static int runBwa(const njh::progutils::CmdArgs & inputCommands);

	static int testHasProgram(const njh::progutils::CmdArgs & inputCommands);
	static int generatingPrime3TemplatesBasedOnMALN(const njh::progutils::CmdArgs & inputCommands);



};

} /* namespace njhseq */





