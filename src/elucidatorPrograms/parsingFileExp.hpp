#pragma once
/*
 * parsingFileExpRunner.h
 *
 *  Created on: July 4, 2017
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
#include <njhseq/common.h>
#include <njhseq/programUtils/seqSetUp.hpp>


namespace njhseq {

class parsingFileExpRunner : public njh::progutils::ProgramRunner {
 public:
	parsingFileExpRunner();

	static int parsingBioSampleSetNCBIJson(const njh::progutils::CmdArgs & inputCommands);
	static int getStrainInfoTableBioSampleSetNCBIJson(const njh::progutils::CmdArgs & inputCommands);

	static int getAttributeLevelsBioSampleSetNCBIJson(const njh::progutils::CmdArgs & inputCommands);
	static int parseBlastpHitsTab(const njh::progutils::CmdArgs & inputCommands);
	static int BlastpHitsTabToBed(const njh::progutils::CmdArgs & inputCommands);

	static int parsePrimerFastaToPrimerTxt(const njh::progutils::CmdArgs & inputCommands);


	static int parseSTOCKHOLM(const njh::progutils::CmdArgs & inputCommands);
	static int parseSTOCKHOLMToFasta(const njh::progutils::CmdArgs & inputCommands);

	static int parsehmmerDomainHitTab(const njh::progutils::CmdArgs & inputCommands);


	static int quickCountFasta(const njh::progutils::CmdArgs & inputCommands);
	static int quickCountFastq(const njh::progutils::CmdArgs & inputCommands);

	static int quickCountDirectory(const njh::progutils::CmdArgs & inputCommands);

	static int parseNucmerResultsToBed(const njh::progutils::CmdArgs & inputCommands);

};
} /* namespace njhseq */


