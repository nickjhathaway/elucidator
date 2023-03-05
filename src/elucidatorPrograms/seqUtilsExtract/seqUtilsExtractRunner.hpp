#pragma once
//

//  seqUtilsExtractRunner.hpp
//
//  Created by Nicholas Hathaway on 2015/05/28.
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

#include "seqUtilsExtractSetUp.hpp"

namespace njhseq {

class seqUtilsExtractRunner : public njh::progutils::ProgramRunner {
 public:
  seqUtilsExtractRunner();
  
  static int extractSeqsBeginsWith(const njh::progutils::CmdArgs & inputCommands);
  static int extractSeqsEndsWith(const njh::progutils::CmdArgs & inputCommands);
  static int extractSeqsBeginsWithEndsWith(const njh::progutils::CmdArgs & inputCommands);

  static int extractByMIDs(const njh::progutils::CmdArgs & inputCommands);

  //extractByName
	static int extractByName(const njh::progutils::CmdArgs & inputCommands);

  static int extractSameSeqs(const njh::progutils::CmdArgs & inputCommands);
  static int extractBySeq(const njh::progutils::CmdArgs & inputCommands);


  static int binOnNucComp(const njh::progutils::CmdArgs & inputCommands);
  static int binOnNucCompFaster(const njh::progutils::CmdArgs & inputCommands);
  static int greedyKmerCluster(const njh::progutils::CmdArgs & inputCommands);
	static int clusterByKmerSim(const njh::progutils::CmdArgs & inputCommands);

  static int countIlluminaAaptors(const njh::progutils::CmdArgs & inputCommands);
  static int extractByIlluminaAaptors(const njh::progutils::CmdArgs & inputCommands);
};
} // namespace njhseq
