#pragma once
//

//  seqUtilsSplitRunner.hpp
//
//  Created by Nicholas Hathaway on 2018/11/18.
//  Copyright (c) 2016 Nicholas Hathaway. All rights reserved.
//

#include <njhseq.h>


namespace njhseq {

class seqUtilsSplitRunner : public njh::progutils::ProgramRunner {
 public:
  seqUtilsSplitRunner();

  static int SeqSplitOnLenBelow(const njh::progutils::CmdArgs & inputCommands);
  static int SeqSplitOnLenWithin(const njh::progutils::CmdArgs & inputCommands);
  static int SeqSplitOnNameContains(const njh::progutils::CmdArgs & inputCommands);
  static int SeqSplitOnSeqContains(const njh::progutils::CmdArgs & inputCommands);
  static int SeqSplitOnLenAbove(const njh::progutils::CmdArgs & inputCommands);
  static int SeqSplitOnLenBetween(const njh::progutils::CmdArgs & inputCommands);
  static int SeqSplitOnQualityWindow(const njh::progutils::CmdArgs & inputCommands);
  static int SeqSplitOnNucelotideComp(const njh::progutils::CmdArgs & inputCommands);


  static int getSimilarSequences(const njh::progutils::CmdArgs & inputCommands);
  static int getSimilarSequencesByKDist(const njh::progutils::CmdArgs & inputCommands);

};
} // namespace njhseq
