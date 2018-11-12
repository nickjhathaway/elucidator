#pragma once
/*
 * jsonExpRunner.h
 *
 *  Created on: Dec 19, 2017
 *      Author: nickhathaway
 */

#include <njhcpp/progutils.h>
#include <njhseq.h>

namespace njhseq {

class jsonExpRunner : public njh::progutils::ProgramRunner {
 public:
	jsonExpRunner();
	static int jsonPrintFieldInArray(const njh::progutils::CmdArgs & inputCommands);
	static int jsonExtractInArrayMatchingField(const njh::progutils::CmdArgs & inputCommands);

};
} /* namespace njhseq */


