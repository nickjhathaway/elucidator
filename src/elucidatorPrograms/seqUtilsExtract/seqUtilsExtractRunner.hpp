#pragma once
//

//  seqUtilsExtractRunner.hpp
//
//  Created by Nicholas Hathaway on 2015/05/28.
//  Copyright (c) 2015 Nicholas Hathaway. All rights reserved.
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


  static int extractSameSeqs(const njh::progutils::CmdArgs & inputCommands);
  static int extractBySeq(const njh::progutils::CmdArgs & inputCommands);
  static int getSimilarSequences(const njh::progutils::CmdArgs & inputCommands);
  static int getSimilarSequencesByKDist(const njh::progutils::CmdArgs & inputCommands);

  static int binOnNucComp(const njh::progutils::CmdArgs & inputCommands);
  static int binOnNucCompFaster(const njh::progutils::CmdArgs & inputCommands);
  static int greedyKmerCluster(const njh::progutils::CmdArgs & inputCommands);
};
} // namespace njhseq
