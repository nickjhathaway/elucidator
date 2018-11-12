#pragma once
//
//  graphicsUtils.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 05/30/2013.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
//

#include "graphicsUtilsSetUp.hpp"

namespace njhseq {

class graphicsUtilsRunner : public njh::progutils::ProgramRunner {
 public:
  graphicsUtilsRunner();

  static int printAnsiColors(const njh::progutils::CmdArgs & inputCommands);
  static int getColors(const njh::progutils::CmdArgs & inputCommands);
  static int colorInfo(const njh::progutils::CmdArgs & inputCommands);
  static int multipleColorsInfo(const njh::progutils::CmdArgs & inputCommands);

};
}  // namespace njhseq

