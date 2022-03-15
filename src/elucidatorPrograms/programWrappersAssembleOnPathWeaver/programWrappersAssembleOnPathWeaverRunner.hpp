#pragma once

/*
 * programWrappersAssembleOnPathWeaverRunner.hpp
 *
 *  Created on: Sep 8, 2021
 *      Author: nick
 */


#include <njhcpp/progutils.h>
#include <njhseq/common.h>
#include <njhseq/programUtils/seqSetUp.hpp>

namespace njhseq {


class programWrappersAssembleOnPathWeaverRunner : public njh::progutils::ProgramRunner {
 public:
	programWrappersAssembleOnPathWeaverRunner();

	static int runUnicyclerOnPathWeaverRegions(const njh::progutils::CmdArgs & inputCommands);
	static int runSpadesOnPathWeaverRegions(const njh::progutils::CmdArgs & inputCommands);
	static int runMegahitOnPathWeaverRegions(const njh::progutils::CmdArgs & inputCommands);
	static int runSavageOnPathWeaverRegions(const njh::progutils::CmdArgs & inputCommands);
	static int runVelvetOptimizerAndMetaVelvetOnPathWeaverRegions(const njh::progutils::CmdArgs & inputCommands);
	static int runPRICEOnPathWeaverRegions(const njh::progutils::CmdArgs & inputCommands);
	static int runTrinityOnPathWeaverRegions(const njh::progutils::CmdArgs & inputCommands);
	static int runRayOnPathWeaverRegions(const njh::progutils::CmdArgs & inputCommands);
	static int runIDBAUDOnPathWeaverRegions(const njh::progutils::CmdArgs & inputCommands);
	static int runMIRAOnPathWeaverRegions(const njh::progutils::CmdArgs & inputCommands);

};



}  // namespace njhseq






