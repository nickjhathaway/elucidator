#pragma once
//
//  ampliconAnalysisRunner.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 10/24/13.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
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

};
}  // namespace njhseq


