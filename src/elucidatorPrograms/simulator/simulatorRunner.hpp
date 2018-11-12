#pragma once
//
//  simulator.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 01/26/2013.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
//

#include "simulatorSetUp.hpp"

namespace njhseq {

class simulatorRunner : public njh::progutils::ProgramRunner {
 public:
  simulatorRunner();

  static int randomStrsWithKmers(const njh::progutils::CmdArgs & inputCommands);
  static int randomStrings(const njh::progutils::CmdArgs & inputCommands);
  static int randomSeqFile(const njh::progutils::CmdArgs & inputCommands);
  static int randomSampleFile(const njh::progutils::CmdArgs & inputCommands);
  static int randomSampleFast(const njh::progutils::CmdArgs & inputCommands);

};
}  // namespace njhseq
