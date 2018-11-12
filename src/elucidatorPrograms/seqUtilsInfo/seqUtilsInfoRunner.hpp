#pragma once
//

//  seqUtilsInfoRunner.hpp
//
//  Created by Nicholas Hathaway on 2015/05/28.
//  Copyright (c) 2015 Nicholas Hathaway. All rights reserved.
//

#include "seqUtilsInfoSetUp.hpp"

namespace njhseq {

class seqUtilsInfoRunner : public njh::progutils::ProgramRunner {
 public:
  seqUtilsInfoRunner();
  
  static int findSeq(const njh::progutils::CmdArgs & inputCommands);
  static int profileReadsToReference(const njh::progutils::CmdArgs & inputCommands);
  static int countSeqsFile(const njh::progutils::CmdArgs & inputCommands);
  static int fastaIdenticalInfo(const njh::progutils::CmdArgs & inputCommands);
  static int getHpProfile(const njh::progutils::CmdArgs & inputCommands);
  static int countHPRuns(const njh::progutils::CmdArgs & inputCommands);
  static int countOtus(const njh::progutils::CmdArgs & inputCommands);
  static int profileQualityScores(const njh::progutils::CmdArgs & inputCommands);
  static int quickLenInfo(const njh::progutils::CmdArgs & inputCommands);
  static int countLetters(const njh::progutils::CmdArgs & inputCommands);
  static int countSeqPortion(const njh::progutils::CmdArgs & inputCommands);
  static int printTandems(const njh::progutils::CmdArgs & inputCommands);
  static int quickMismatchDist(const njh::progutils::CmdArgs & inputCommands);
  static int countKmers(const njh::progutils::CmdArgs & inputCommands);
  static int countKmersPlusStats(const njh::progutils::CmdArgs & inputCommands);
  static int profileErrors(const njh::progutils::CmdArgs & inputCommands);
  static int countAllSeqs(const njh::progutils::CmdArgs & inputCommands);
  static int printNames(const njh::progutils::CmdArgs & inputCommands);
  static int printSeqs(const njh::progutils::CmdArgs & inputCommands);
  static int getGCContent(const njh::progutils::CmdArgs & inputCommands);

  static int fracInfo(const njh::progutils::CmdArgs & inputCommands);
  static int genPsuedoMismatchMinTree(const njh::progutils::CmdArgs & inputCommands);
  static int genPsuedoAllMinTree(const njh::progutils::CmdArgs & inputCommands);


  static int qualCounts(const njh::progutils::CmdArgs & inputCommands);

  static int quickHaplotypeInformation(const njh::progutils::CmdArgs & inputCommands);
  static int quickHaplotypeInformationAndVariants(const njh::progutils::CmdArgs & inputCommands);
  static int quickHaplotypeVariantsWithRegion(const njh::progutils::CmdArgs & inputCommands);

  static int multipleAlnProteinToPcaInput(const njh::progutils::CmdArgs & inputCommands);

  static int getReadLens(const njh::progutils::CmdArgs & inputCommands);
  static int readLengthDistribution(const njh::progutils::CmdArgs & inputCommands);


};
} // namespace njhseq
