#pragma once
//

//  seqUtilsSplitRunner.hpp
//
//  Created by Nicholas Hathaway on 2018/11/18.
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
#include <njhseq/programUtils/seqSetUp.hpp>


namespace njhseq {

class seqUtilsSplitRunner : public njh::progutils::ProgramRunner {
 public:
  seqUtilsSplitRunner();

  static int SeqSplitOnLenBelow(const njh::progutils::CmdArgs & inputCommands);
  static int SeqSplitOnLenWithin(const njh::progutils::CmdArgs & inputCommands);
  static int SeqSplitOnNameContains(const njh::progutils::CmdArgs & inputCommands);
  static int SeqSplitOnNameContainsPattern(const njh::progutils::CmdArgs & inputCommands);

  static int SeqSplitOnSeqContains(const njh::progutils::CmdArgs & inputCommands);
  static int SeqSplitOnLenAbove(const njh::progutils::CmdArgs & inputCommands);
  static int SeqSplitOnLenBetween(const njh::progutils::CmdArgs & inputCommands);
  static int SeqSplitOnQualityWindow(const njh::progutils::CmdArgs & inputCommands);
  static int SeqSplitOnQualityCheck(const njh::progutils::CmdArgs & inputCommands);

  static int SeqSplitOnNucelotideComp(const njh::progutils::CmdArgs & inputCommands);
  static int SeqSplitOnCount(const njh::progutils::CmdArgs & inputCommands);


  static int getSimilarSequences(const njh::progutils::CmdArgs & inputCommands);
  static int getSimilarSequencesByKDist(const njh::progutils::CmdArgs & inputCommands);

};
} // namespace njhseq
