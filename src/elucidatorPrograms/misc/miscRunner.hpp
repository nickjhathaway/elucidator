#pragma once
//
//  miscRunner.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 11/3/13.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
//
#include "miscSetUp.hpp"

namespace njhseq {

class miscRunner : public njh::progutils::ProgramRunner {
 public:
  miscRunner();
  // misc useful
  static int listAllFiles(const njh::progutils::CmdArgs & inputCommands);

  // specific programs

  static int hangman(const njh::progutils::CmdArgs & inputCommands);


  static int getSharedLines(const njh::progutils::CmdArgs & inputCommands);
  static int catColor(const njh::progutils::CmdArgs & inputCommands);

  static int aminoAcidQuiz(const njh::progutils::CmdArgs & inputCommands);

  static int genInputForEstimateS(const njh::progutils::CmdArgs & inputCommands);
  static int getHighestHapFrac(const njh::progutils::CmdArgs & inputCommands);

  static int increaseQualityScores(const njh::progutils::CmdArgs & inputCommands);

  static int findMotifLocations(const njh::progutils::CmdArgs & inputCommands);
  static int findTandemMotifLocations(const njh::progutils::CmdArgs & inputCommands);

  static int parseBamForForwardPerfectHits(const njh::progutils::CmdArgs & inputCommands);

  static int expandOutCollapsedToUnique(const njh::progutils::CmdArgs & inputCommands);

  static int guessAProteinFromSeq(const njh::progutils::CmdArgs & inputCommands);


  static int expandTableBySeparatingColumn(const njh::progutils::CmdArgs & inputCommands);

  static int codeComparison(const njh::progutils::CmdArgs & inputCommands);

  static int getSlidingQualityWindowMeans(const njh::progutils::CmdArgs & inputCommands);

  static int isFileEmpty(const njh::progutils::CmdArgs & inputCommands);

  static int createConnectedHaplotypeNetwork(const njh::progutils::CmdArgs & inputCommands);

  static int getLinkedInfoFromAdjList(const njh::progutils::CmdArgs & inputCommands);
  static int getAlnPosToRealPosTable(const njh::progutils::CmdArgs & inputCommands);
  static int calculateShannonEntropySimple(const njh::progutils::CmdArgs & inputCommands);




};
}  // namespace njhseq


