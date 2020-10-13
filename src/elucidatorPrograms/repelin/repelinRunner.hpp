#pragma once
/*
 * repelinRunner.hpp
 *
 *  Created on: Jan 11, 2015
 *      Author: nickhathaway
 */
// elucidator - A library for analyzing sequence data
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
//
// This file is part of elucidator.
//
// elucidator is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// elucidator is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with elucidator.  If not, see <http://www.gnu.org/licenses/>.
//
#include <njhseq/common.h>
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

	static int runTRF(const njh::progutils::CmdArgs & inputCommands);

};

} /* namespace njhseq */


