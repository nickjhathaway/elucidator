//
// Created by Nicholas Hathaway on 2/19/24.
//

#include "kmerSetExp.hpp"

#include <njhseq/IO/SeqIO.h>
#include <njhseq/objects/kmer/SimpleKmerHash.hpp>
#include <njhseq/objects/kmer/KmerSetInvertedIndex.hpp>
#include <njhseq/objects/kmer/KmerCollectionReadMatcher.hpp>
#include "elucidator/helpers/UniqueKmerSetHelper.hpp"


namespace njhseq {

int kmerSetExpRunner::extractByKmerUniqueSets(const njh::progutils::CmdArgs &inputCommands) {
  bfs::path kmersPerTargetFnp;
  uint32_t minLenCutOff = 100;
  uint32_t maxLenCutOff = 10000;
  uint32_t hardKmersPerTarget = 5;
  double kmerFracPerTarget = 0.10;
  uint32_t numThreads = 1;
  std::string sampleName = "sample";
  bool rename = false;

  seqSetUp setUp(inputCommands);
  setUp.description_ = "Assign reads to the target whose k-mers they best match, the k-mer step of SeekDeep's extractorByKmerMatching without the primer, MID, length and quality checks";
  setUp.setOption(kmersPerTargetFnp, "--kmersPerTarget,--uniqueKmersPerTarget", "Kmers Per Target, 1)set,2)kmer, all kmers per target (targets can overlap) or only unique kmers", true, "IDs");

  //in and out
  setUp.processDefaultReader(VecStr{"--fastq", "--fastqgz", "--fasta", "--fastagz"});
  setUp.processDirectoryOutputName(true);
  if("sample" == sampleName){
    sampleName = bfs::basename(setUp.pars_.ioOptions_.firstName_);
    if(njh::endsWith(setUp.pars_.ioOptions_.firstName_.string(), ".gz")){
      sampleName = bfs::path(sampleName).replace_extension("").string();
    }
  }
  setUp.setOption(sampleName, "--sampleName", "sample name");
  setUp.setOption(rename, "--rename", "rename input extracted reads");

  // filtering
  setUp.setOption(minLenCutOff, "--minLenCutOff", "Hard cut off min length", false, "filtering");
  setUp.setOption(maxLenCutOff, "--maxLenCutOff", "Hard cut off max length", false, "filtering");
  setUp.setOption(hardKmersPerTarget, "--hardKmersPerTarget", "hard occurrence cut off Kmers Per Target", false, "filtering");
  setUp.setOption(kmerFracPerTarget, "--kmerFracPerTarget", "at minimum have to have this fraction of kmers in a Target", false, "filtering");

  //running
  setUp.setOption(numThreads, "--numThreads", "number of threads");
  bool noInvertedKmerIndex = false;
  setUp.setOption(noInvertedKmerIndex, "--noInvertedKmerIndex",
      "Probe every target's k-mer set for every read k-mer instead of using one k-mer -> target(s) index, the counts are identical");

  setUp.finishSetUp(std::cout);
  setUp.startARunLog(setUp.pars_.directoryName_);
  setUp.writeParametersFile(setUp.pars_.directoryName_ + "parametersUsed.txt", false, false);

  uint32_t extractionKmer = UniqueKmerSetHelper::getKmerLenFromUniqueKmerTable(kmersPerTargetFnp);
  auto allKmersPerTarget = UniqueKmerSetHelper::readInUniqueKmerTablePerSet(kmersPerTargetFnp);
  if(allKmersPerTarget.empty()){
    std::stringstream ss;
    ss << __PRETTY_FUNCTION__ << ", error " << "no targets read in from " << kmersPerTargetFnp << "\n";
    throw std::runtime_error { ss.str() };
  }
  //all kmers pick the winning target, kmers unique to a target count how many targets a read hits
  const KmerCollectionReadMatcher kmerMatcher(allKmersPerTarget, extractionKmer, !noInvertedKmerIndex);
  if(!noInvertedKmerIndex){
    auto invIndex = KmerSetInvertedIndex::get(allKmersPerTarget);
    std::cout << "kmerSetExclusivity\tsets\t" << invIndex->nSets_
        << "\ttotalKmers\t" << invIndex->totalKmers_
        << "\tsharedKmers\t" << invIndex->nSharedKmers_
        << "\texclusive\t" << njh::boolToStr(0 == invIndex->nSharedKmers_) << std::endl;
  }

  SeqInput reader(setUp.pars_.ioOptions_);
  reader.openIn();

  VecStr targets = njh::getVecOfMapKeys(allKmersPerTarget);
  njh::sort(targets);
  MultiSeqIO seqOut;
  auto extractionDir = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar{"extractedReads"});
  for(const auto & name : targets){
    seqOut.addReader(name, SeqIOOptions::genFastqOutGz(njh::files::make_path(extractionDir, name)));
  }
  seqOut.addReader("undetermined", SeqIOOptions::genFastqOutGz(njh::files::make_path(setUp.pars_.directoryName_, "undetermined")));
  seqOut.addReader("multihit", SeqIOOptions::genFastqOutGz(njh::files::make_path(setUp.pars_.directoryName_, "multihit")));

  uint64_t totalReads = 0;
  uint64_t smallFrags = 0;
  uint64_t largeFrags = 0;
  std::unordered_map<std::string, uint64_t> readsPerSet;
  std::unordered_map<std::string, uint64_t> readsPerSetRevComp;
  //reads with enough of the winning target's k-mers in both orientations
  std::unordered_map<std::string, uint64_t> inverseHitsPerSet;

  std::mutex mut;
  std::function<void()> extractReads = [&]() {
    auto kmerMatch = kmerMatcher.genMatch();
    seqInfo seq;
    uint64_t totalReadsCurrent = 0;
    uint64_t smallFragsCurrent = 0;
    uint64_t largeFragsCurrent = 0;
    std::unordered_map<std::string, uint64_t> readsPerSetCurrent;
    std::unordered_map<std::string, uint64_t> readsPerSetRevCompCurrent;
    std::unordered_map<std::string, uint64_t> inverseHitsPerSetCurrent;

    while(reader.readNextReadLock(seq)){
      ++totalReadsCurrent;
      if(len(seq) < minLenCutOff){
        ++smallFragsCurrent;
        continue;
      }
      if(len(seq) > maxLenCutOff){
        ++largeFragsCurrent;
        continue;
      }
      //all kmers pick the winner, kmers unique to a target count the targets hit
      kmerMatch.reset();
      kmerMatcher.addKmers(seq.seq_, true, kmerMatch);
      kmerMatcher.countFound(kmerMatch);

      //zero out targets that don't pass the hard count and fraction cut offs
      double nInputKmers = static_cast<double>(kmerMatch.kmers.size());
      for(auto * found : {&kmerMatch.found, &kmerMatch.foundRevComp, &kmerMatch.uniqueFound, &kmerMatch.uniqueFoundRevComp}){
        for(auto & count : *found){
          if(count < hardKmersPerTarget || count / nInputKmers <= kmerFracPerTarget){
            count = 0;
          }
        }
      }
      kmerMatcher.determineWinner(kmerMatch);
      std::string winnerSet = kmerMatch.hasWinner() ? kmerMatcher.names()[kmerMatch.winnerIdx] : "undetermined";
      const bool winnerRevComp = kmerMatch.winnerRevComp;
      if(KmerCollectionReadMatcher::countCollectionsHit(kmerMatch, hardKmersPerTarget) > 1){
        winnerSet = "multihit";
      }
      if(winnerRevComp){
        ++readsPerSetRevCompCurrent[winnerSet];
        seq.reverseComplementRead(true, true);
      } else {
        ++readsPerSetCurrent[winnerSet];
      }
      if("undetermined" != winnerSet && "multihit" != winnerSet &&
         kmerMatch.found[kmerMatch.winnerIdx] >= hardKmersPerTarget && kmerMatch.foundRevComp[kmerMatch.winnerIdx] >= hardKmersPerTarget){
        ++inverseHitsPerSetCurrent[winnerSet];
      }
      if(rename){
        auto threadId = estd::to_string(std::this_thread::get_id());
        seq.name_ = njh::pasteAsStr(sampleName, ".", winnerSet, ".", readsPerSetCurrent[winnerSet] + readsPerSetRevCompCurrent[winnerSet], ".", threadId);
        if(winnerRevComp){
          seq.name_.append("_Comp");
        }
      }
      seqOut.openWrite(winnerSet, seq);
    }
    {
      std::lock_guard<std::mutex> lockGuard(mut);
      totalReads += totalReadsCurrent;
      smallFrags += smallFragsCurrent;
      largeFrags += largeFragsCurrent;
      for(const auto & count : readsPerSetCurrent){
        readsPerSet[count.first] += count.second;
      }
      for(const auto & count : readsPerSetRevCompCurrent){
        readsPerSetRevComp[count.first] += count.second;
      }
      for(const auto & count : inverseHitsPerSetCurrent){
        inverseHitsPerSet[count.first] += count.second;
      }
    }
  };
  njh::concurrent::runVoidFunctionThreaded(extractReads, numThreads);

  auto frac = [](uint64_t num, uint64_t denom){
    return 0 == denom ? 0.0 : static_cast<double>(num) / static_cast<double>(denom);
  };
  uint64_t totalReadsProcessed = totalReads - smallFrags - largeFrags;

  OutputStream outCounts(njh::files::make_path(setUp.pars_.directoryName_, "extractionProfile.tab.txt"));
  outCounts << "sample\ttotalReadsProcessed\ttarget\tcount\tfrac\tforwardCount\tfracForward\tinverseHits\tinverseHitsFrac" << std::endl;
  uint64_t totalExtracted = 0;
  uint64_t totalExtractedForward = 0;
  auto outNames = targets;
  outNames.emplace_back("multihit");
  outNames.emplace_back("undetermined");
  for(const auto & name : outNames){
    uint64_t extracted = readsPerSet[name] + readsPerSetRevComp[name];
    if("multihit" != name && "undetermined" != name){
      totalExtracted += extracted;
      totalExtractedForward += readsPerSet[name];
    }
    outCounts << sampleName
              << "\t" << totalReadsProcessed
              << "\t" << name
              << "\t" << extracted
              << "\t" << frac(extracted, totalReadsProcessed)
              << "\t" << readsPerSet[name]
              << "\t" << frac(readsPerSet[name], extracted)
              << "\t" << inverseHitsPerSet[name]
              << "\t" << frac(inverseHitsPerSet[name], extracted)
              << std::endl;
  }

  OutputStream outStats(njh::files::make_path(setUp.pars_.directoryName_, "extractionStats.tab.txt"));
  outStats << "sampleName\ttotalReads\tfailedMinLen_" << minLenCutOff << "\tfailedMinLenFrac"
           << "\tfailedMaxLen_" << maxLenCutOff << "\tfailedMaxLenFrac"
           << "\tundetermined\tundeterminedFrac\tmultihit\tmultihitFrac"
           << "\textracted\textractedFrac\textractedForward\textractedForwardFrac" << std::endl;
  uint64_t undetermined = readsPerSet["undetermined"] + readsPerSetRevComp["undetermined"];
  uint64_t multihit = readsPerSet["multihit"] + readsPerSetRevComp["multihit"];
  outStats << sampleName
           << "\t" << totalReads
           << "\t" << smallFrags
           << "\t" << frac(smallFrags, totalReads)
           << "\t" << largeFrags
           << "\t" << frac(largeFrags, totalReads)
           << "\t" << undetermined
           << "\t" << frac(undetermined, totalReads)
           << "\t" << multihit
           << "\t" << frac(multihit, totalReads)
           << "\t" << totalExtracted
           << "\t" << frac(totalExtracted, totalReads)
           << "\t" << totalExtractedForward
           << "\t" << frac(totalExtractedForward, totalExtracted)
           << std::endl;
  return 0;
}


} //namespace njhseq
