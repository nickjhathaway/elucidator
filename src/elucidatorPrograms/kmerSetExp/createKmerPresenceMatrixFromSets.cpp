//
// Created by Nicholas Hathaway on 9/4/24.
//

#include <elucidatorPrograms/programWrappersAssembleOnPathWeaver/otherAssemblersUtils.hpp>

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
	// naturalSortNameSet(fastaFilesStrs);
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





int kmerSetExpRunner::getUniqueKmersFromRandomSubsamples(const njh::progutils::CmdArgs & inputCommands) {
  OutOptions outOpts("", ".tsv.gz");
  uint32_t kmerLength = 16;
  uint32_t subsampleStart = 10;
  uint32_t subsampleEnd = 100;
  uint32_t subsampleStep = 10;
  uint32_t subsampleRuns = 10;
  std::string id = "id";
  seqSetUp setUp(inputCommands);
  setUp.description_ = "Get info on how many kmers can be found and at what counts";
  setUp.processVerbose();
  setUp.processDebug();
  setUp.setOption(id, "--id", "id for file", true);

  setUp.setOption(kmerLength, "--kmerLength", "kmer Length", true);
  setUp.setOption(subsampleStart, "--subsampleStart", "subsample Start", true);
  setUp.setOption(subsampleEnd, "--subsampleEnd", "subsample End", true);
  setUp.setOption(subsampleStep, "--subsampleStep", "subsample Step", true);
  setUp.setOption(subsampleRuns, "--subsampleRuns", "subsample Runs", true);

  setUp.processReadInNames(true);
  setUp.processWritingOptions(outOpts);
  setUp.finishSetUp(std::cout);

  auto input = createKmerReadVec(SeqInput::getSeqVec<readObject>(setUp.pars_.ioOptions_), kmerLength, false);

  OutputStream out(outOpts);
  out << "id\tkmerLength\tsubSampleAmount\tsubSampleRunID\tuniqueKmerCount" << std::endl;
  njh::randomGenerator rGen;
  std::vector<uint32_t> positions(input.size());
  njh::iota(positions,0U);
  for(uint32_t sub = subsampleStart; sub < subsampleEnd && sub <= input.size(); sub += subsampleStep) {
    for(uint32_t run = 0; run < subsampleRuns; run++) {
      auto sampledPositions = rGen.unifRandSelectionVec(positions, sub, false);
      std::unordered_set<std::string> kmers;
      for(const auto & position : sampledPositions) {
        njh::addVecToUOSet(njh::getVecOfMapKeys(input[position]->kInfo_.kmers_), kmers);
      }
      out << id
          << "\t" << kmerLength
          << "\t" << sub
          << "\t" << run
          << "\t" << kmers.size() << std::endl;
    }
  }
  return 0;
}



} //namespace njhseq

