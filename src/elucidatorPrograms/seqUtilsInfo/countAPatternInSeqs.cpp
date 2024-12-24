//
// Created by Nicholas Hathaway on 12/24/24.
//



#include "seqUtilsInfoRunner.hpp"
#include <njhseq/IO/SeqIO/SeqIO.hpp>


namespace njhseq {

int seqUtilsInfoRunner::countAPatternInSeqs(const njh::progutils::CmdArgs & inputCommands) {

  OutOptions outOpts(bfs::path(""), ".tsv");
  std::string patternStr = "G{10,}$";
  seqSetUp setUp(inputCommands);
  setUp.description_ = "count the presence of a pattern within input sequences";

  setUp.processVerbose();
  setUp.processDebug();
  setUp.processWritingOptions(outOpts);
  setUp.setOption(patternStr, "--pattern", "pattern to search for", true);
  setUp.processReadInNames(true);
  setUp.finishSetUp(std::cout);

  SeqIO reader(setUp.pars_.ioOptions_);
  reader.openIn();
  OutputStream out(outOpts);

  std::regex pattern{patternStr};
  if (setUp.pars_.ioOptions_.isPairedIn()) {
    out << "inputFnp\tpattern\tr1_foundCount\tr2_coundCount\tinBoth_foundCount\ttotalInputPairedCount" << std::endl;

    PairedRead pseq;
    uint32_t r1_foundCount = 0;
    uint32_t r2_coundCount = 0;
    uint32_t inBoth_foundCount = 0;
    uint32_t totalInputPairedCount = 0;
    while(reader.readNextRead(pseq)){
      ++totalInputPairedCount;
      std::smatch firstMate_match;
      if (std::regex_search(pseq.seqBase_.seq_, firstMate_match, pattern)) {
        ++r1_foundCount;
      }
      std::smatch secondMate_match;
      if (std::regex_search(pseq.mateSeqBase_.seq_, secondMate_match, pattern)) {
        ++r2_coundCount;
      }
      if (firstMate_match.size() == 1 && secondMate_match.size() == 1) {
        ++inBoth_foundCount;
      }
    }
    out << setUp.pars_.ioOptions_.firstName_
        << "\t" << patternStr
        << "\t" << r1_foundCount
        << "\t" << r2_coundCount
        << "\t" << inBoth_foundCount
        << "\t" << totalInputPairedCount << std::endl;
  } else {
    seqInfo seq;
    out << "inputFnp\tpattern\tfoundCount\ttotalInputCount" << std::endl;
    uint32_t totalInputCount = 0;
    uint32_t foundCount = 0;
    while(reader.readNextRead(seq)){
      ++totalInputCount;
      std::smatch match;
      if (std::regex_search(seq.seq_, match, pattern)) {
        ++foundCount;
      }
    }
    out << setUp.pars_.ioOptions_.firstName_
        << "\t" << patternStr
        << "\t" << foundCount
        << "\t" << totalInputCount << std::endl;
  }
  if(setUp.pars_.verbose_){
    setUp.logRunTime(std::cout);
  }
  return 0;
}


} //namespace njhseq

