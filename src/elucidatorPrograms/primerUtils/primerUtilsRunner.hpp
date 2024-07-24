#pragma once
//
// Created by Nicholas Hathaway on 5/22/23.
//


#include <njhcpp/progutils.h>
#include <njhseq/common.h>
#include <njhseq/programUtils/seqSetUp.hpp>

namespace njhseq {

class primerUtilsRunner : public njh::progutils::ProgramRunner {
public:
	primerUtilsRunner();

	static int computeDimerizationScore(const njh::progutils::CmdArgs & inputCommands);
	static int creatingMultiplexAmpliconPools(const njh::progutils::CmdArgs & inputCommands);
	static int testForPrimerDimers(const njh::progutils::CmdArgs & inputCommands);

	static int testWithBlastForUnspecificAmplification(const njh::progutils::CmdArgs & inputCommands);

};

} //  namespace njhseq

