#pragma once
//
// Created by Nicholas Hathaway on 8/20/25.
//


#include <njhcpp/progutils.h>
#include <njhseq/common.h>
#include <njhseq/programUtils/seqSetUp.hpp>

namespace njhseq {

class PMOUtilsRunner : public njh::progutils::ProgramRunner {
public:
  PMOUtilsRunner();

  static int read_pmo(const njh::progutils::CmdArgs & inputCommands);


};

} //  namespace njhseq

