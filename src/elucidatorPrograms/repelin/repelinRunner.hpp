#pragma once
/*
 * repelinRunner.hpp
 *
 *  Created on: Jan 11, 2015
 *      Author: nickhathaway
 */

#include <njhseq.h>
#include <njhcpp/progutils/programRunner.hpp>


namespace njhseq {

class repelinRunner : public njh::progutils::ProgramRunner {
public:
	repelinRunner();
	static int parseRepeatMaskerOutputUnitTest(const njh::progutils::CmdArgs & inputCommands);

	static int parseRepeatMaskerOutputForElement(const njh::progutils::CmdArgs & inputCommands);

	static int extractElementSequences(const njh::progutils::CmdArgs & inputCommands);
	static int extractFullElementSequences(const njh::progutils::CmdArgs & inputCommands);

	static int parseRMForElementWithoutIntervening(const njh::progutils::CmdArgs & inputCommands);


	static int parseTandemRepeatFinderOutputUnitTest(const njh::progutils::CmdArgs & inputCommands);

	static int TandemRepeatFinderOutputToBed(const njh::progutils::CmdArgs & inputCommands);

};

} /* namespace njhseq */


