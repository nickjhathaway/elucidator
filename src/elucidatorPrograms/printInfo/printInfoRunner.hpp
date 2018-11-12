#pragma once
//

//  printInfoRunner.hpp
//
//  Created by Nicholas Hathaway on 2015/05/29.
//  Copyright (c) 2015 Nicholas Hathaway. All rights reserved.
//

#include "printInfoSetUp.hpp"

namespace njhseq {

class printInfoRunner : public njh::progutils::ProgramRunner {
 public:
  printInfoRunner();
  
  static int printDegen(const njh::progutils::CmdArgs & inputCommands);
  static int printAminoAcidInfo(const njh::progutils::CmdArgs & inputCommands);
  static int printFastqAscII(const njh::progutils::CmdArgs & inputCommands);

};
} // namespace njhseq
