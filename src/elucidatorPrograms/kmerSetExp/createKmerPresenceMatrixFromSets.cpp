//
// Created by Nicholas Hathaway on 9/4/24.
//

#include <elucidatorPrograms/programWrappersAssembleOnPathWeaver/otherAssemblersUtils.hpp>

#include "kmerSetExp.hpp"

#include <njhseq/IO/SeqIO.h>
#include <njhseq/objects/kmer/KmerGatherer.hpp>
#include <njhseq/objects/kmer/SimpleKmerHash.hpp>
#include <SeekDeep/parameters/setUpPars.hpp>
#include <boost/math/special_functions/binomial.hpp>

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
    out << "\t" << bfs::basename(njh::files::replaceExtension(fnp, ""));
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


long double lchoose_aprox(uint64_t n, uint64_t k) {
 // https://math.stackexchange.com/questions/64716/approximating-the-logarithm-of-the-binomial-coefficient?newreg=2fc6dc6d632741b88f67bb756f5354e3
  return n * logl(n) - k * logl(k) - (n - k) * logl(n - k) + 0.5 * (logl(n) - logl(k) - logl(n - k) - logl(2 * M_PI));
}


int kmerSetExpRunner::estimateKmerSubSamples(const njh::progutils::CmdArgs & inputCommands) {
  OutOptions outOpts("", ".tsv.gz");
  uint32_t kmerLength = 16;
  uint32_t subsampleStart = 10;
  uint32_t subsampleEnd = 100;
  uint32_t subsampleStep = 10;
  bool byReadLength = false;
  std::string id = "id";
  seqSetUp setUp(inputCommands);
  setUp.description_ = "Get info on how many kmers can be found and at what counts";
  setUp.processVerbose();
  setUp.processDebug();

  setUp.setOption(id, "--id", "id for file", true);
  setUp.setOption(byReadLength, "--byReadLength", "aproximate the sub-sampling k-mer count by stepping by aproximate read addition, e.g. if mean read length is 100 instead of stepping by step do (read_len - kmerlen +1) steps");


  setUp.setOption(kmerLength, "--kmerLength", "kmer Length", true);
  setUp.setOption(subsampleStart, "--subsampleStart", "subsample Start", true);
  setUp.setOption(subsampleEnd, "--subsampleEnd", "subsample End", true);
  setUp.setOption(subsampleStep, "--subsampleStep", "subsample Step", true);


  setUp.processReadInNames(true);
  setUp.processWritingOptions(outOpts);
  setUp.finishSetUp(std::cout);

  auto input = createKmerReadVec(SeqInput::getSeqVec<readObject>(setUp.pars_.ioOptions_), kmerLength, false);
  std::vector<uint64_t> r_lengths;
  r_lengths.reserve(input.size());
  for (const auto & seq : input) {
    r_lengths.emplace_back(len(seq->seqBase_));
  }
  uint64_t meanReadLen =static_cast<uint64_t>(std::round(vectorMedianCopy(r_lengths)));

  std::unordered_map<std::string, uint64_t> allKmerCounts;
  uint64_t totalKmers = 0;
  for (const auto & seq : input) {
    for (const auto & kmer : seq->kInfo_.kmers_) {
      allKmerCounts[kmer.first] += kmer.second.count_;
      totalKmers += kmer.second.count_;
    }
  }


  if (setUp.pars_.verbose_) {
    std::cout << "nput.size(): " << input.size() << std::endl;
    std::cout << "meanReadLen: " << meanReadLen << std::endl;
    std::cout << "totalKmers: " << totalKmers << std::endl;
    std::cout << "allKmerCounts.size(): " << allKmerCounts.size() << std::endl;
    if (byReadLength) {
      std::cout << "lchoose_aprox(" << totalKmers << ", " << subsampleStart * (meanReadLen - kmerLength + 1) << ") " << ": " << lchoose_aprox(totalKmers, subsampleStart * (meanReadLen - kmerLength + 1)) << std::endl;
      std::cout << "lchoose_aprox(" << totalKmers << ", " << subsampleEnd * (meanReadLen - kmerLength + 1)<< ") " << ": " << lchoose_aprox(totalKmers, subsampleEnd * (meanReadLen - kmerLength + 1)) << std::endl;
    } else {
      std::cout << "lchoose_aprox(" << totalKmers << ", " << subsampleStart << ") " << ": " << lchoose_aprox(totalKmers, subsampleStart) << std::endl;
      std::cout << "lchoose_aprox(" << totalKmers << ", " << subsampleEnd << ") " << ": " << lchoose_aprox(totalKmers, subsampleEnd) << std::endl;
    }
  }

  OutputStream out(outOpts);
  if (byReadLength) {
    out << "id\tkmerLength\tseqs\tsubSampleAmount\trarefiedUniqueKmerCount\ttotalKmers\ttotalUnique" << std::endl;
    for(uint64_t read_sub = subsampleStart; read_sub < subsampleEnd && read_sub <= input.size() && read_sub * (meanReadLen - kmerLength + 1 ) < totalKmers; read_sub += subsampleStep) {
      uint64_t sub = read_sub * (meanReadLen - kmerLength + 1);
      long double bottom = lchoose_aprox(totalKmers, sub);
      long double sum = 0;
      for (const auto & k : allKmerCounts) {
        uint64_t bigN_minus_species_n = totalKmers - k.second;
        long double top = lchoose_aprox(bigN_minus_species_n, sub);
        sum += 1 - exp(top - bottom);
      }
      out << id
      << "\t" << kmerLength
      << "\t" << read_sub
      << "\t" << sub
      << "\t" << sum
      << "\t" << totalKmers
      << "\t" << allKmerCounts.size() << std::endl;
    }
  } else {
    out << "id\tkmerLength\tsubSampleAmount\trarefiedUniqueKmerCount\ttotalKmers\ttotalUnique" << std::endl;
    for(uint64_t sub = subsampleStart; sub < subsampleEnd && sub <= totalKmers; sub += subsampleStep) {
      long double bottom = lchoose_aprox(totalKmers, sub);
      long double sum = 0;
      for (const auto & k : allKmerCounts) {
        uint64_t bigN_minus_species_n = totalKmers - k.second;
        long double top = lchoose_aprox(bigN_minus_species_n, sub);
        sum += 1 - exp(top - bottom);
      }
      out << id
      << "\t" << kmerLength
      << "\t" << sub
      << "\t" << sum
      << "\t" << totalKmers
      << "\t" << allKmerCounts.size() << std::endl;
    }
  }



  return 0;
}

int kmerSetExpRunner::rarifySequencesSubSamples(const njh::progutils::CmdArgs & inputCommands) {
  OutOptions outOpts("", ".tsv.gz");
  uint32_t subsampleStart = 10;
  uint32_t subsampleEnd = 100;
  uint32_t subsampleStep = 10;
  bool byReadLength = false;
  std::string id = "id";
  seqSetUp setUp(inputCommands);
  setUp.description_ = "Get info on how many kmers can be found and at what counts";
  setUp.processVerbose();
  setUp.processDebug();

  setUp.setOption(id, "--id", "id for file", true);
  setUp.setOption(subsampleStart, "--subsampleStart", "subsample Start", true);
  setUp.setOption(subsampleEnd, "--subsampleEnd", "subsample End", true);
  setUp.setOption(subsampleStep, "--subsampleStep", "subsample Step", true);


  setUp.processReadInNames(true);
  setUp.processWritingOptions(outOpts);
  setUp.finishSetUp(std::cout);

  auto input = SeqInput::getSeqVec<seqInfo>(setUp.pars_.ioOptions_);

  std::unordered_map<std::string, uint64_t> allSeqsCounts;
  uint64_t totalSeqs = 0;
  for (const auto & seq : input) {
    allSeqsCounts[seq.seq_] += seq.cnt_;
    totalSeqs += seq.cnt_;
  }


  if (setUp.pars_.verbose_) {
    std::cout << "input.size(): " << input.size() << std::endl;
    std::cout << "totalSeqs: " << totalSeqs << std::endl;
    std::cout << "allSeqsCounts.size(): " << allSeqsCounts.size() << std::endl;
  }

  OutputStream out(outOpts);
  out << "id\tsubSampleAmount\trarefiedUniqueSeqCount\ttotalSeqs\ttotalUniqueSeqs" << std::endl;
  for(uint64_t sub = subsampleStart; sub < subsampleEnd && sub <= totalSeqs; sub += subsampleStep) {
    long double bottom = lchoose_aprox(totalSeqs, sub);
    long double sum = 0;
    for (const auto & seq : allSeqsCounts) {
      uint64_t bigN_minus_species_n = totalSeqs - seq.second;
      long double top = lchoose_aprox(bigN_minus_species_n, sub);
      sum += 1 - exp(top - bottom);
    }
    out << id
    << "\t" << sub
    << "\t" << sum
    << "\t" << totalSeqs
    << "\t" << allSeqsCounts.size() << std::endl;
  }
  return 0;
}



} //namespace njhseq

