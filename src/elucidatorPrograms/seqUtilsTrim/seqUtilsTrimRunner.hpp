#pragma once
//

//  seqUtilsTrimRunner.hpp
//
//  Created by Nicholas Hathaway on 2016/10/05.
//  Copyright (c) 2016 Nicholas Hathaway. All rights reserved.
//

#include "seqUtilsTrimSetUp.hpp"

namespace njhseq {

class seqUtilsTrimRunner : public njh::progutils::ProgramRunner {
 public:
  seqUtilsTrimRunner();


  static int trimFront(const njh::progutils::CmdArgs & inputCommands);
  static int trimEnd(const njh::progutils::CmdArgs & inputCommands);
  static int trimEdges(const njh::progutils::CmdArgs & inputCommands);
  static int trimBeforeSeq(const njh::progutils::CmdArgs & inputCommands);
  static int trimFromSeq(const njh::progutils::CmdArgs & inputCommands);
  static int trimBetweenSeqs(const njh::progutils::CmdArgs & inputCommands);
  static int trimToLen(const njh::progutils::CmdArgs & inputCommands);
  static int trimToSimilarSeq(const njh::progutils::CmdArgs & inputCommands);

  static int trimAtFirstQual(const njh::progutils::CmdArgs & inputCommands);
  static int trimAtFirstBase(const njh::progutils::CmdArgs & inputCommands);
  static int trimAtLastBase(const njh::progutils::CmdArgs & inputCommands);

  static int trimFromMostProbableSharedKmer(const njh::progutils::CmdArgs & inputCommands);
  static int trimToMostProbableSharedKmer(const njh::progutils::CmdArgs & inputCommands);

  static int trimFromMostCommonKmer(const njh::progutils::CmdArgs & inputCommands);
  static int trimToMostCommonKmer(const njh::progutils::CmdArgs & inputCommands);
  static int trimBetweenMostCommonKmers(const njh::progutils::CmdArgs & inputCommands);

  static int trimWithMuscle(const njh::progutils::CmdArgs & inputCommands);

  static int trimWithMuscleMaxStartMinEnd(const njh::progutils::CmdArgs & inputCommands);

  static int trimWithMuscleToRef(const njh::progutils::CmdArgs & inputCommands);
  static int trimWithMuscleToRefInStreaks(const njh::progutils::CmdArgs & inputCommands);


  static int trimLstripBase(const njh::progutils::CmdArgs & inputCommands);
  static int trimRstripBase(const njh::progutils::CmdArgs & inputCommands);
  static int trimStripBase(const njh::progutils::CmdArgs & inputCommands);
  static int trimLstripQual(const njh::progutils::CmdArgs & inputCommands);
  static int trimRstripQual(const njh::progutils::CmdArgs & inputCommands);
  static int trimStripQual(const njh::progutils::CmdArgs & inputCommands);

  static int trimToPositions(const njh::progutils::CmdArgs & inputCommands);
  static int trimToPositionsForEachName(const njh::progutils::CmdArgs & inputCommands);

  static int trimToRefWithGlobalAlignment(const njh::progutils::CmdArgs & inputCommands);
  static int trimToRefWithGlobalAlignmentToRefPositions(const njh::progutils::CmdArgs & inputCommands);

  static int trimToRefWithGlobalAlignmentToRefMultiplePositions(const njh::progutils::CmdArgs & inputCommands);


  static int trimWithSlidingQualityAvgWindow(const njh::progutils::CmdArgs & inputCommands);


  static int breakupAtRegexPat(const njh::progutils::CmdArgs & inputCommands);

  //
  //
};
} // namespace njhseq
