#pragma once
//

//  seqUtilsConverterRunner.hpp
//
//  Created by Nicholas Hathaway on 2015/05/28.
//  Copyright (c) 2015 Nicholas Hathaway. All rights reserved.
//

#include "seqUtilsConverterSetUp.hpp"

namespace njhseq {

class seqUtilsConverterRunner : public njh::progutils::ProgramRunner {
 public:
  seqUtilsConverterRunner();
  
  static int convertFiles(const njh::progutils::CmdArgs & inputCommands);
  static int convertTxtToFasta(const njh::progutils::CmdArgs & inputCommands);
  static int convertFastaToTxt(const njh::progutils::CmdArgs & inputCommands);
  static int convertTxtOneLineToFasta(const njh::progutils::CmdArgs & inputCommands);
  static int sffInfo(const njh::progutils::CmdArgs & inputCommands);
  static int plasmoDBTxtToFasta(const njh::progutils::CmdArgs & inputCommands);

};
} // namespace njhseq
