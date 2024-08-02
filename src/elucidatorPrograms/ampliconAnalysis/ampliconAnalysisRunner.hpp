#pragma once
//
//  ampliconAnalysisRunner.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 10/24/13.
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

#include "ampliconAnalysisSetUp.hpp"

namespace njhseq {

class ampliconAnalysisRunner : public njh::progutils::ProgramRunner {

 public:
  ampliconAnalysisRunner();

  static int collapseTandems(const njh::progutils::CmdArgs & inputCommands);
  static int markChimeras(const njh::progutils::CmdArgs & inputCommands);
  static int greedyCluster(const njh::progutils::CmdArgs & inputCommands);
  static int singleLinkageClusteringOnPerId(const njh::progutils::CmdArgs & inputCommands);

	static int processRawExtractByKmerPathWeaverResults(const njh::progutils::CmdArgs & inputCommands);

	static int determinePossibleMaskFromSeqs(const njh::progutils::CmdArgs & inputCommands);
	static int maskRegionBasedOnRefSubRegions(const njh::progutils::CmdArgs & inputCommands);


	// PMO
  static int extractedTarAmpInfoFileToJson(const njh::progutils::CmdArgs & inputCommands);
  static int finalClustersFileToJson(const njh::progutils::CmdArgs & inputCommands);
  static int specimenInfoFileToJson(const njh::progutils::CmdArgs & inputCommands);
	static int experimentInfoFileToJson(const njh::progutils::CmdArgs & inputCommands);
	static int demultiplexedExperimentSampleFileToJson(const njh::progutils::CmdArgs & inputCommands);

	static int combingAllIntoPMOJson(const njh::progutils::CmdArgs & inputCommands);

};
}  // namespace njhseq


