#pragma once
/*
 * programWrappers.hpp
 *
 *  Created on: Feb 2, 2017
 *      Author: nick
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
	static int runBowtieOnAdapterReomvalOutputSinglesCombined(const njh::progutils::CmdArgs & inputCommands);
	static int processAdaptorRemovalLog(const njh::progutils::CmdArgs & inputCommands);


	static int runPicardMarkDups(const njh::progutils::CmdArgs & inputCommands);

	static int runBwaOnAdapterReomvalOutputSE(const njh::progutils::CmdArgs & inputCommands);


	static int runBwa(const njh::progutils::CmdArgs & inputCommands);

	static int testHasProgram(const njh::progutils::CmdArgs & inputCommands);
	static int generatingPrime3TemplatesBasedOnMALN(const njh::progutils::CmdArgs & inputCommands);

	static int runDada2(const njh::progutils::CmdArgs & inputCommands);
	static int runDada2SingleSamplePaired(const njh::progutils::CmdArgs & inputCommands);


};

} /* namespace njhseq */





