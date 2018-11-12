#pragma once
//

//  seqUtilsModRunner.hpp
//
//  Created by Nicholas Hathaway on 2015/05/28.
//  Copyright (c) 2015 Nicholas Hathaway. All rights reserved.
//

#include "seqUtilsModSetUp.hpp"

namespace njhseq {

class seqUtilsModRunner : public njh::progutils::ProgramRunner {
 public:
  seqUtilsModRunner();
  
  static int prependReads(const njh::progutils::CmdArgs & inputCommands);
  static int appendReads(const njh::progutils::CmdArgs & inputCommands);

  static int prependNames(const njh::progutils::CmdArgs & inputCommands);
  static int appendNames(const njh::progutils::CmdArgs & inputCommands);


  static int sortReads(const njh::progutils::CmdArgs & inputCommands);
  static int split(const njh::progutils::CmdArgs & inputCommands);
  static int renameIDs(const njh::progutils::CmdArgs & inputCommands);

  static int removeLowQualityBases(const njh::progutils::CmdArgs & inputCommands);
  static int translate(const njh::progutils::CmdArgs & inputCommands);
  static int revCompSeq(const njh::progutils::CmdArgs & inputCommands);
  static int reOrientReads(const njh::progutils::CmdArgs & inputCommands);
  static int collapseToUnique(const njh::progutils::CmdArgs & inputCommands);
  static int collapseToUniqueWithInMetaField(const njh::progutils::CmdArgs & inputCommands);

  static int dereplicate(const njh::progutils::CmdArgs & inputCommands);


  static int inverseLetterCase(const njh::progutils::CmdArgs & inputCommands);
  static int changeLetterCase(const njh::progutils::CmdArgs & inputCommands);

};
} // namespace njhseq
