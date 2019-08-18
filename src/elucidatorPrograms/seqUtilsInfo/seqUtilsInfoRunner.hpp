#pragma once
//

//  seqUtilsInfoRunner.hpp
//
//  Created by Nicholas Hathaway on 2015/05/28.
// elucidator - A library for analyzing sequence data
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
//
// This file is part of elucidator.
//
// elucidator is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// elucidator is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with elucidator.  If not, see <http://www.gnu.org/licenses/>.
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
  static int getHapPopDifAndVariantsInfo(const njh::progutils::CmdArgs & inputCommands);
  static int oldQuickHaplotypeInformationAndVariants(const njh::progutils::CmdArgs & inputCommands);
  static int quickHaplotypeVariantsWithRegion(const njh::progutils::CmdArgs & inputCommands);

  static int multipleAlnProteinToPcaInput(const njh::progutils::CmdArgs & inputCommands);

  static int getReadLens(const njh::progutils::CmdArgs & inputCommands);
  static int readLengthDistribution(const njh::progutils::CmdArgs & inputCommands);


};
} // namespace njhseq
