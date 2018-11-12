#pragma once
/*
 * fileFormatExpRunner.h
 *
 *  Created on: May 18, 2015
 *      Author: nickhathaway
 */

#include <njhcpp/progutils.h>
#include <njhseq.h>

namespace njhseq {

class fileFormatExpRunner : public njh::progutils::ProgramRunner {
 public:
	fileFormatExpRunner();
	static int extractRefSeqRecords(const njh::progutils::CmdArgs & inputCommands);

	static int parsePf3kEmblFilesToGff3(const njh::progutils::CmdArgs & inputCommands);
};
} /* namespace njhseq */


