#pragma once
//
// Created by Nicholas Hathaway on 8/20/25.
//


#include <njhcpp/progutils.h>
#include <njhseq/common.h>
#include <njhseq/programUtils/seqSetUp.hpp>

namespace njhseq {

class PhipSeqUtilsRunner : public njh::progutils::ProgramRunner {
public:
  PhipSeqUtilsRunner();

  static int fragmentSequencesForPhipSeq(const njh::progutils::CmdArgs & inputCommands);
  static int generateAllPossibleNucleotidePossibleFromProtein(const njh::progutils::CmdArgs & inputCommands);
  static int generateNucleotidePossibleFromProteins(const njh::progutils::CmdArgs & inputCommands);

  static int markGroupsByHammingDistanceCutOff(const njh::progutils::CmdArgs & inputCommands);

  static int appendRandomBarcode(const njh::progutils::CmdArgs & inputCommands);
  static int countPossiblePhipSeqRandomBarcodes(const njh::progutils::CmdArgs & inputCommands);

};

} //  namespace njhseq



