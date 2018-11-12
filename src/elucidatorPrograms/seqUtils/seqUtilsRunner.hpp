#pragma once
//
//  seqUtilsRunner.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 10/17/13.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
//

#include "seqUtilsSetUp.hpp"

namespace njhseq {

class seqUtilsRunner : public njh::progutils::ProgramRunner {
 public:
  seqUtilsRunner();

  static int createDegenerativeStr(const njh::progutils::CmdArgs & inputCommands);
  static int createConsensus(const njh::progutils::CmdArgs & inputCommands);

  static int compareAllByAll(const njh::progutils::CmdArgs & inputCommands);
  static int compareToRef(const njh::progutils::CmdArgs & inputCommands);
  static int mapCount(const njh::progutils::CmdArgs & inputCommands);

  static int alignToSequence(const njh::progutils::CmdArgs & inputCommands);
  static int checkTwoReadFiles(const njh::progutils::CmdArgs & inputCommands);
};






}  // namespace njhseq
