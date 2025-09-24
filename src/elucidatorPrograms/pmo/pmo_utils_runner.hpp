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
  static int get_overlap_between_panels_in_pmos(const njh::progutils::CmdArgs & inputCommands);

  static int add_protein_variant_info_to_pmo(const njh::progutils::CmdArgs & inputCommands);
  static int count_protein_variant_info_to_pmo(const njh::progutils::CmdArgs & inputCommands);

};

} //  namespace njhseq

