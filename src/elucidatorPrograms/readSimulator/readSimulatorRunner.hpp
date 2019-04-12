#pragma once
//
//  readSimulator.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 12/23/15.
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
  //static int createLibrarySimMultipleMixture(const njh::progutils::CmdArgs & inputCommands);
  static int createLibrarySimMultipleMixtureDrugResistant(const njh::progutils::CmdArgs & inputCommands);


  static int effectsOfPcrErrorRatePerRounds(const njh::progutils::CmdArgs & inputCommands);
};
}  // namespace njhseq

