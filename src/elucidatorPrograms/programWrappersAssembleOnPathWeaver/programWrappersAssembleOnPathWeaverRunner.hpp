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

	static int runSpadesOnPathWeaverRegions(const njh::progutils::CmdArgs & inputCommands);
	static int runMegahitOnPathWeaverRegions(const njh::progutils::CmdArgs & inputCommands);




};



}  // namespace njhseq






