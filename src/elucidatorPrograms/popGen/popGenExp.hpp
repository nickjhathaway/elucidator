#pragma once
/*
 * popGenExp.hpp
 *
 *  Created on: Jun 11, 2021
 *      Author: nicholas hathaway
 */

// elucidator - A library for analyzing sequence data
// Copyright (C) 2012-2021 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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

class popGenExpRunner : public njh::progutils::ProgramRunner {
 public:
	popGenExpRunner();

	static int doPairwiseComparisonOnHapsSharing(const njh::progutils::CmdArgs & inputCommands);
	static int tajimatest_testingExample(const njh::progutils::CmdArgs & inputCommands);
	static int tajimatest(const njh::progutils::CmdArgs & inputCommands);


  static int callVariantsAgainstRefSeq(const njh::progutils::CmdArgs & inputCommands);
  static int callVariantsAgainstRefSeqIndividual(const njh::progutils::CmdArgs & inputCommands);

  static int quickHaplotypeInformation(const njh::progutils::CmdArgs & inputCommands);
  static int getHapPopDifAndVariantsInfo(const njh::progutils::CmdArgs & inputCommands);
  static int oldQuickHaplotypeInformationAndVariants(const njh::progutils::CmdArgs & inputCommands);
  static int quickHaplotypeVariantsWithRegion(const njh::progutils::CmdArgs & inputCommands);

};

} //  namespace njhseq







