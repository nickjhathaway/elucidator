//
// Created by Nicholas Hathaway on 7/1/24.
//
#include "kmerSetExp.hpp"

#include <njhseq/IO/SeqIO.h>
#include <njhseq/objects/kmer/SimpleKmerHash.hpp>
#include <SeekDeep/parameters/setUpPars.hpp>


namespace njhseq {

int kmerSetExpRunner::getBasicKmerSetCountInfo(const njh::progutils::CmdArgs & inputCommands) {
  uint32_t kmerLength = 7;

  OutOptions outOpts(bfs::path(""), ".tab.txt.gz");
  seqSetUp setUp(inputCommands);
  setUp.description_ = "Get info on how many kmers can be found and at what counts";
  setUp.processVerbose();
  setUp.processDebug();
  setUp.setOption(kmerLength, "--kmerLength", "kmer Length", true);
  setUp.processReadInNames(true);
  setUp.processWritingOptions(outOpts);
  setUp.finishSetUp(std::cout);

  SeqInput reader(setUp.pars_.ioOptions_);
  reader.openIn();

  OutputStream out(outOpts);

  seqInfo seq;
  std::unordered_map<std::string, uint32_t> kmerCounts;
  std::unordered_map<std::string, std::unordered_set<uint32_t>> kmerSampleCounts;
  uint32_t seqCounts = 0;
  while(reader.readNextRead(seq)) {
    if(len(getSeqBase(seq).seq_) >=kmerLength) {
      for (uint32_t pos = 0; pos < len(getSeqBase(seq).seq_) - kmerLength + 1; ++pos) {
        auto k = getSeqBase(seq).seq_.substr(pos, kmerLength);
        ++kmerCounts[k];
        kmerSampleCounts[k].emplace(seqCounts);
      }
      ++seqCounts;
    }
  }
  out << "totalSeqs\tkmerLength\tkmerSeqOccurrences\tfreq\tnumberOfKmers" << std::endl;
  std::map<uint32_t, uint32_t> countOfKmersWithSampleCounts;
  for(const auto & kCount : kmerSampleCounts) {
    ++countOfKmersWithSampleCounts[kCount.second.size()];
  }

  for (const auto &count: countOfKmersWithSampleCounts) {
    out << seqCounts << "\t" << kmerLength
        << "\t" << count.first
        << "\t" << count.first / static_cast<double>(seqCounts)
        << "\t" << count.second
        << std::endl;
  }
  return 0;
}

} //namespace njhseq