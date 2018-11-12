#pragma once
//
//  readSimulator.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 12/23/15.
//  Copyright (c) 2015 Nicholas Hathaway. All rights reserved.
//

#include "readSimulatorSetUp.hpp"

namespace njhseq {

class readSimulatorRunner : public njh::progutils::ProgramRunner {
 public:
  readSimulatorRunner();

  static int chimeraSim(const njh::progutils::CmdArgs & inputCommands);
  static int createFractionAbundanceFile(const njh::progutils::CmdArgs & inputCommands);

  static int createRandomSequenceMixtures(const njh::progutils::CmdArgs & inputCommands);
  static int createMinorVariantsFromKnown(const njh::progutils::CmdArgs & inputCommands);

  static int simSingleSequencePCR(const njh::progutils::CmdArgs & inputCommands);
  static int simPcrShotgunSequences(const njh::progutils::CmdArgs & inputCommands);
  static int simMultipleMixturePCR(const njh::progutils::CmdArgs & inputCommands);

  static int createIlluminaErrorProfile(const njh::progutils::CmdArgs & inputCommands);
  static int simMultipleMixture(const njh::progutils::CmdArgs & inputCommands);
  static int createLibrarySimMultipleMixture(const njh::progutils::CmdArgs & inputCommands);
  static int createLibrarySimMultipleMixtureDrugResistant(const njh::progutils::CmdArgs & inputCommands);


  static int effectsOfPcrErrorRatePerRounds(const njh::progutils::CmdArgs & inputCommands);
};
}  // namespace njhseq

