#pragma once
/*
 * bamCountingExpRunner.h
 *
 *  Created on: Feb 4, 2015
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


