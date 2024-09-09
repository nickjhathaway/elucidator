//
// Created by Nicholas Hathaway on 9/4/24.
//

#include "kmerSetExp.hpp"

#include <njhseq/IO/SeqIO.h>
#include <njhseq/objects/kmer/KmerGatherer.hpp>
#include <njhseq/objects/kmer/SimpleKmerHash.hpp>
#include <SeekDeep/parameters/setUpPars.hpp>


namespace njhseq {


int kmerSetExpRunner::createKmerPresenceMatrixFromSets(const njh::progutils::CmdArgs & inputCommands) {
  OutOptions outOpts("", ".tsv.gz");
  KmerGatherer::KmerGathererPars countPars;
  countPars.noRevComp_ = true;
  countPars.kmerLength_ = 7;
  countPars.entropyFilter_ = 0;
  std::vector<bfs::path> fastaFiles;
  seqSetUp setUp(inputCommands);
  setUp.description_ = "Get info on how many kmers can be found and at what counts";
  setUp.processVerbose();
  setUp.processDebug();
  setUp.setOption(countPars.kmerLength_, "--kmerLength", "kmer Length", true);
  setUp.setOption(fastaFiles, "--fastaFiles", "kmer files", true);
  setUp.setOption(countPars.numThreads_, "--numThreads", "num Threads");
  setUp.processWritingOptions(outOpts);
  setUp.finishSetUp(std::cout);

  KmerGatherer kGather(countPars);
  OutputStream out(outOpts);

  auto allKmers = kGather.getUniqueKmersSetFromFastas(fastaFiles);
  std::set<std::string> allUniqueKmers;
  for(const auto & set : allKmers) {
    allUniqueKmers.insert(set.second.begin(), set.second.end());
  }
  out << "kmer";
  std::vector<std::string> fastaFilesStrs;
  fastaFilesStrs.reserve(fastaFiles.size());
  for(const auto & fnp : fastaFiles) {
    fastaFilesStrs.emplace_back(fnp.string());
  }
  njh::naturalSortNameSet(fastaFilesStrs);
  for(const auto & fnp : fastaFilesStrs) {
    out << "\t" << bfs::basename(fnp);
  }
  out << std::endl;

  for(const auto & k : allUniqueKmers) {
    out << k;
    for(const auto & fnp : fastaFilesStrs) {
      out << "\t" << (njh::in(k, allKmers[fnp]) ? 1 : 0);
    }
    out << std::endl;
  }
  return 0;
}


} //namespace njhseq

