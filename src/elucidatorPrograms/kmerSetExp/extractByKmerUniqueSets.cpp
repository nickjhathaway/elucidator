//
// Created by Nicholas Hathaway on 2/19/24.
//

#include "kmerSetExp.hpp"

#include <njhseq/IO/SeqIO.h>
#include <njhseq/objects/kmer/SimpleKmerHash.hpp>
#include <SeekDeep/parameters/setUpPars.hpp>


namespace njhseq {

int kmerSetExpRunner::extractByKmerUniqueSets(const njh::progutils::CmdArgs &inputCommands) {
  bfs::path idFnp;
  bfs::path lenCutOffsFnp;
  bfs::path uniqueKmersPerTargetFnp;
  CoreExtractorPars corePars;
  corePars.smallFragmentCutoff = 100;


  uint32_t numThreads = 1;

  std::string sampleName = "sample";
  bool rename = false;

  seqSetUp setUp(inputCommands);
  //id
  setUp.setOption(idFnp, "--ids,--id", "Primers file", true, "IDs");

  setUp.setOption(lenCutOffsFnp, "--lenCutOffsFnp,--lenCutOffs", "length Cut Offs per target", true, "IDs");
  setUp.setOption(uniqueKmersPerTargetFnp, "--uniqueKmersPerTarget", "unique Kmers Per Target", true, "IDs");

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
  setUp.setOption(corePars.smallFragmentCutoff, "--minLenCutOff", "Hard cut off min length", false, "filtering");

  //running
  setUp.setOption(numThreads, "--numThreads", "number of threads");



  setUp.finishSetUp(std::cout);
  // run log
  setUp.startARunLog(setUp.pars_.directoryName_);
  // parameter file
  setUp.writeParametersFile(setUp.pars_.directoryName_ + "parametersUsed.txt", false, false);

  // create Primers and MIDs
  PrimersAndMids ids(idFnp);

  if(ids.getTargets().empty()){
    std::stringstream ss;
    ss << __PRETTY_FUNCTION__ << ", error " << "no targets read in from " << idFnp << "\n";
    throw std::runtime_error { ss.str() };
  }


  //init primer determinator
  ids.initPrimerDeterminator();

  //read in extra info
  ids.addLenCutOffs(lenCutOffsFnp);
  uint32_t extractionKmer = ids.addUniqKmerCounts(uniqueKmersPerTargetFnp);


  // set up input
  seqInfo seq;
  SeqInput reader(setUp.pars_.ioOptions_);
  reader.openIn();



  //qual check off;
  std::unique_ptr<ReadChecker> qualChecker = std::make_unique<ReadCheckerQualCheck>(corePars.qPars_.qualCheck_,
                                                                                    corePars.qPars_.qualCheckCutOff_, true);

  //set up extraction counts outputs
  ExtractionStator masterCounts;
  OutOptions outCountsOpts(njh::files::make_path(setUp.pars_.directoryName_, "extractionProfile.tab.txt"));
  OutputStream outCounts(outCountsOpts);

  OutOptions outStatsOpts(njh::files::make_path(setUp.pars_.directoryName_, "extractionStats.tab.txt"));
  OutputStream outStats(outStatsOpts);

  std::unordered_map<std::string, uint32_t> readsPerSet;
  std::unordered_map<std::string, uint32_t> readsPerSetRevComp;

  VecStr names = ids.getTargets();
  MultiSeqIO seqOut;
  auto initialExtractionDir = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar{"extractedReads"});
  auto failedExtractionDir = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar{"failedExtractedReads"});
  for(const auto & name : names){
    auto seqOutOpts = SeqIOOptions::genFastqOutGz(njh::files::make_path(initialExtractionDir, name) );
    seqOut.addReader(name, seqOutOpts);
  }

  for(const auto & name : names){
    auto seqOutOpts = SeqIOOptions::genFastqOutGz(njh::files::make_path(setUp.pars_.directoryName_, njh::pasteAsStr(name, sampleName) ) );
    seqOut.addReader(njh::pasteAsStr(name, "-passed"), seqOutOpts);
  }

  for(const auto & name : names){
    auto seqOutOpts = SeqIOOptions::genFastqOutGz(njh::files::make_path(failedExtractionDir, name) );
    seqOut.addReader(njh::pasteAsStr(name, "-failed"), seqOutOpts);
  }

  seqOut.addReader("undetermined", SeqIOOptions::genFastqOutGz(njh::files::make_path(setUp.pars_.directoryName_, "undetermined")));


  std::mutex mut;
  std::function<void()> readInComp = [&reader, &ids, &readsPerSet,&readsPerSetRevComp,&mut,
                                      &extractionKmer,&seqOut,&sampleName,&rename,
                                      &corePars,
                                      &masterCounts]() {
    SimpleKmerHash hasher;
    seqInfo seq;
    std::unordered_map<std::string, uint32_t> readsPerSetCurrent;
    std::unordered_map<std::string, uint32_t> readsPerSetRevCompCurrent;


    ExtractionStator masterCountsCurrent;

    while(reader.readNextReadLock(seq)){
      ++masterCountsCurrent.totalReadCount_;
      if(len(seq) < corePars.smallFragmentCutoff){
        ++masterCountsCurrent.smallFrags_;
        continue;
      }
      std::unordered_map<uint64_t, uint64_t> hashedInputKmers;
      std::unordered_map<uint64_t, uint64_t> hashedInputKmersRevComp;
      if(len(seq.seq_) > extractionKmer){
        for(uint32_t pos = 0; pos < len(seq.seq_) - extractionKmer + 1; ++pos){
          auto hash = hasher.hash(seq.seq_.substr(pos, extractionKmer));
          ++hashedInputKmers[hash];
        }
      }
      if(len(seq.seq_) > extractionKmer){
        for(uint32_t pos = 0; pos < len(seq.seq_) - extractionKmer + 1; ++pos){
          auto hash = hasher.revCompHash(seq.seq_.substr(pos, extractionKmer));
          ++hashedInputKmersRevComp[hash];
        }
      }

      std::unordered_map<std::string, uint32_t> foundPerSet;
      std::unordered_map<std::string, uint32_t> foundPerSetRevComp;
      for(const auto & setName  : njh::getVecOfMapKeys(ids.uniqueKmersPerTarget_)){
        foundPerSet[setName] = 0;
        foundPerSetRevComp[setName] = 0;
      }
      for(const auto & hashedKmer : hashedInputKmers){
        for(const auto & uniqueKmers : ids.uniqueKmersPerTarget_){
          if(njh::in(hashedKmer.first, uniqueKmers.second)){
            ++foundPerSet[uniqueKmers.first];
          }
        }
      }
      for(const auto & hashedKmer : hashedInputKmersRevComp){
        for(const auto & uniqueKmers : ids.uniqueKmersPerTarget_){
          if(njh::in(hashedKmer.first, uniqueKmers.second)){
            ++foundPerSetRevComp[uniqueKmers.first];
          }
        }
      }
      std::string winnerSet = "undetermined";
      double bestFrac = 0;
      bool winnerRevComp = false;

      for(const auto & setName  : njh::getVecOfMapKeys(ids.uniqueKmersPerTarget_)){
        if(static_cast<double>(foundPerSet[setName])/static_cast<double>(hashedInputKmers.size()) > bestFrac){
          bestFrac = static_cast<double>(foundPerSet[setName])/static_cast<double>(hashedInputKmers.size());
          winnerSet = setName;
        }
        if(static_cast<double>(foundPerSetRevComp[setName])/static_cast<double>(hashedInputKmers.size()) > bestFrac){
          bestFrac = static_cast<double>(foundPerSetRevComp[setName])/static_cast<double>(hashedInputKmers.size());
          winnerSet = setName;
          winnerRevComp = true;
        }
      }
      if(winnerRevComp){
        ++readsPerSetRevCompCurrent[winnerSet];
        seq.reverseComplementRead(true, true);
      }else{
        ++readsPerSetCurrent[winnerSet];
      }
      if(rename){
        auto threadId = estd::to_string(std::this_thread::get_id());
        seq.name_ = njh::pasteAsStr(sampleName, ".",winnerSet, ".", readsPerSetCurrent[winnerSet] + readsPerSetRevCompCurrent[winnerSet], ".", threadId);
        if(winnerRevComp){
          seq.name_.append("_Comp");
        }
      }
      seqOut.openWrite(winnerSet, seq);
      ++masterCountsCurrent.readsUnrecBarcode_;
    }
    {
      std::lock_guard<std::mutex> lockGuard(mut);

      masterCounts.addOtherExtractorCounts(masterCountsCurrent);

      for(const auto & readsPerSetCount : readsPerSetCurrent){
        readsPerSet[readsPerSetCount.first] += readsPerSetCount.second;
      }
      for(const auto & readsPerSetRevCompCount : readsPerSetRevCompCurrent){
        readsPerSetRevComp[readsPerSetRevCompCount.first] += readsPerSetRevCompCount.second;
      }
    }
  };

  njh::concurrent::runVoidFunctionThreaded(readInComp, numThreads);


  uint64_t totalReadsProcessed = 0;
  for(const auto & readsPerSetCount : readsPerSet){
    totalReadsProcessed += readsPerSetCount.second;
  }
  for(const auto & readsPerSetRevCompCount : readsPerSetRevComp){
    totalReadsProcessed += readsPerSetRevCompCount.second;
  }

  uint64_t totalReadsPassedAllFilters = 0;


  outCounts << "sample\ttotalReadsProcessed\ttarget\tcount\tfrac\tforwardCount\tfracForward";
  outCounts << "\tminLenFailed\tminLenFailedFrac";
  outCounts << "\tmaxLenFailed\tmaxLenFailedFrac";
  outCounts << "\tqualityFailed\tqualityFailedFrac";
  outCounts << "\tforwardPrimerFailed\tforwardPrimerFailedFrac";
  outCounts << "\treversePrimerFailed\treversePrimerFailedFrac";
  outCounts << "\tbothForRevPrimerFailed\tbothForRevPrimerFailedFrac";
  outCounts << "\ttotalFailed\ttotalFailedFrac";
  outCounts << "\tpassed\tpassedFrac";
  outCounts << std::endl;


  uint64_t totalExtractedAllTargets = 0;
  uint64_t totalExtractedAllTargetsForward = 0;
  uint64_t totalExtractedUndetermined = readsPerSet["undetermined"] + readsPerSetRevComp["undetermined"];


  for(const auto & setName : ids.getTargets()){
    uint64_t totalExtracted = readsPerSet[setName] + readsPerSetRevComp[setName];
    totalExtractedAllTargets += totalExtracted;
		totalExtractedAllTargetsForward += readsPerSet[setName];

    uint64_t totalBad = masterCounts.counts_[setName][false].minLenBad_ + masterCounts.counts_[setName][true].minLenBad_ +
        masterCounts.counts_[setName][false].maxLenBad_ + masterCounts.counts_[setName][true].maxLenBad_+
        masterCounts.counts_[setName][false].qualityFailed_ + masterCounts.counts_[setName][true].qualityFailed_ +
        masterCounts.counts_[setName][false].badForward_ + masterCounts.counts_[setName][true].badForward_ +
        masterCounts.counts_[setName][false].badReverse_ + masterCounts.counts_[setName][true].badReverse_ +
        masterCounts.counts_[setName][false].failedBothPrimers_ + masterCounts.counts_[setName][true].failedBothPrimers_;

    outCounts << sampleName
              << "\t" << totalReadsProcessed
              << "\t" << setName
              << "\t" << totalExtracted
              << "\t" << static_cast<double>(totalExtracted) / static_cast<double>(totalReadsProcessed)
              << "\t" << readsPerSet[setName]
              << "\t" << (totalExtracted == 0 ? 0 : static_cast<double>(readsPerSet[setName]) / static_cast<double>(totalExtracted));
    //minlen
    outCounts<< "\t" << masterCounts.counts_[setName][false].minLenBad_ + masterCounts.counts_[setName][true].minLenBad_
    << "\t" << (totalBad == 0 ? 0 : static_cast<double>(masterCounts.counts_[setName][false].minLenBad_ + masterCounts.counts_[setName][true].minLenBad_)/static_cast<double>(totalBad));
    //maxlen
    outCounts << "\t" << masterCounts.counts_[setName][false].maxLenBad_ + masterCounts.counts_[setName][true].maxLenBad_
    << "\t" << (totalBad == 0 ? 0 : static_cast<double>(masterCounts.counts_[setName][false].maxLenBad_ + masterCounts.counts_[setName][true].maxLenBad_)/static_cast<double>(totalBad));
    //quality
    outCounts << "\t" << masterCounts.counts_[setName][false].qualityFailed_ + masterCounts.counts_[setName][true].qualityFailed_
              << "\t" << (totalBad == 0 ? 0 : static_cast<double>(masterCounts.counts_[setName][false].qualityFailed_ + masterCounts.counts_[setName][true].qualityFailed_)/static_cast<double>(totalBad));
    //failed forward
    outCounts << "\t" << masterCounts.counts_[setName][false].badForward_ + masterCounts.counts_[setName][true].badForward_
              << "\t" << (totalBad == 0 ? 0 : static_cast<double>(masterCounts.counts_[setName][false].badForward_ + masterCounts.counts_[setName][true].badForward_)/static_cast<double>(totalBad));
    //bad reverse
    outCounts << "\t" << masterCounts.counts_[setName][false].badReverse_ + masterCounts.counts_[setName][true].badReverse_
              << "\t" << (totalBad == 0 ? 0 : static_cast<double>(masterCounts.counts_[setName][false].badReverse_ + masterCounts.counts_[setName][true].badReverse_)/static_cast<double>(totalBad));
    //both bad reverse and failed forward
    outCounts << "\t" << masterCounts.counts_[setName][false].failedBothPrimers_ + masterCounts.counts_[setName][true].failedBothPrimers_
              << "\t" << (totalBad == 0 ? 0 : static_cast<double>(masterCounts.counts_[setName][false].failedBothPrimers_ + masterCounts.counts_[setName][true].failedBothPrimers_)/static_cast<double>(totalBad));
    //bad
    outCounts << "\t" << totalBad
              << "\t" << (totalExtracted == 0 ? 0 : static_cast<double>(totalBad) / static_cast<double>(totalExtracted));
    //good
    outCounts << "\t" << masterCounts.counts_[setName][false].good_ + masterCounts.counts_[setName][true].good_
              << "\t" << (totalExtracted == 0 ? 0 : static_cast<double>(masterCounts.counts_[setName][false].good_ + masterCounts.counts_[setName][true].good_) / static_cast<double>(totalExtracted))
              << std::endl;
    totalReadsPassedAllFilters += masterCounts.counts_[setName][false].good_ + masterCounts.counts_[setName][true].good_;
  }
  {
    uint64_t totalExtracted = readsPerSet["undetermined"] + readsPerSetRevComp["undetermined"];
    outCounts << sampleName
              << "\t" << totalReadsProcessed
              << "\t" << "undetermined"
              << "\t" << totalExtracted
              << "\t" << (totalReadsProcessed == 0 ? 0 : static_cast<double>(totalExtracted) / static_cast<double>(totalReadsProcessed))
              << "\t" << readsPerSet["undetermined"]
              << "\t" << (totalExtracted == 0 ? 0 : static_cast<double>(readsPerSet["undetermined"]) / static_cast<double>(totalExtracted))
        << "\t" << 0
        << "\t" << 0
        << "\t" << 0
        << "\t" << 0
        << "\t" << 0
        << "\t" << 0
        << "\t" << 0
        << "\t" << 0
        << "\t" << 0
        << "\t" << 0
        << "\t" << 0
        << "\t" << 0
        << "\t" << 0
        << "\t" << 0
        << "\t" << 0
        << "\t" << 0
              << std::endl;
  }

  outStats << "sampleName\ttotalReadsProcessed\tfailedMinLen_" << corePars.smallFragmentCutoff << "\tfailedMinLenFrac\tundetermined\tundeterminedFrac\textracted\textractedFrac\textractedForward\textractedForwardFrac\tpassed\tpassedFrac" << std::endl;
  outStats << sampleName
					 << "\t" << totalReadsProcessed + masterCounts.smallFrags_
					 << "\t" << masterCounts.smallFrags_
					 << "\t" << static_cast<double>(masterCounts.smallFrags_) / static_cast<double>(totalReadsProcessed + masterCounts.smallFrags_)
					 << "\t" << totalExtractedUndetermined
					 << "\t" << static_cast<double>(totalExtractedUndetermined) / static_cast<double>(totalExtractedUndetermined + totalExtractedAllTargets)
					 << "\t" << totalExtractedAllTargets
					 << "\t" << static_cast<double>(totalExtractedAllTargets) / static_cast<double>(totalReadsProcessed + masterCounts.smallFrags_)
					 << "\t" << totalExtractedAllTargetsForward
					 << "\t" << static_cast<double>(totalExtractedAllTargetsForward) / static_cast<double>(totalExtractedAllTargets)
           << "\t" << totalReadsPassedAllFilters
           << "\t" << static_cast<double>(totalReadsPassedAllFilters) / static_cast<double>(totalReadsProcessed + masterCounts.smallFrags_)
           << std::endl;

  return 0;
}



} //namespace njhseq

