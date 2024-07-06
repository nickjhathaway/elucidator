//
// Created by Nicholas Hathaway on 7/1/24.
//
#include "kmerSetExp.hpp"

#include <njhseq/IO/SeqIO.h>
#include <njhseq/objects/seqObjects/seqKmers/KmerVecUtils.hpp>
#include <njhseq/objects/kmer/KmerUtils.hpp>
#include <SeekDeep/parameters/setUpPars.hpp>


namespace njhseq {

int kmerSetExpRunner::getBasicKmerSetCountInfo(const njh::progutils::CmdArgs & inputCommands) {
  uint32_t kmerLength = 7;


  seqSetUp setUp(inputCommands);
  setUp.description_ = "Get info on how many kmers can be found and at what counts";
  setUp.processVerbose();
  setUp.processDebug();
  setUp.setOption(kmerLength, "--kmerLength", "kmer Length", true);
  setUp.processReadInNames(true);
  setUp.processDirectoryOutputName(true);
  setUp.finishSetUp(std::cout);
  setUp.startARunLog(setUp.pars_.directoryName_);
  SeqInput reader(setUp.pars_.ioOptions_);
  reader.openIn();

  OutOptions outKmerInfoOpts(njh::files::make_path(setUp.pars_.directoryName_, "kmerCountsBasics.tsv"));
  OutputStream outKmerInfo(outKmerInfoOpts);

  OutOptions outKmerInfoPerSeqOpts(njh::files::make_path(setUp.pars_.directoryName_, "kmerCountsPerSeqBasics.tsv"));
  OutputStream outKmerInfoPerSeq(outKmerInfoPerSeqOpts);

  outKmerInfo << "totalSeqs\tkmerLength\tkmerSeqOccurrences\tfreq\tnumberOfKmers" << std::endl;
  outKmerInfoPerSeq << "seqName\tlength\tuniqueKmers\tuniqueKmerFrac\tavgNonUniqueKmerOccurrences\tmedianNonUniqueKmerOccurrences\tavgKmerOccurrences\tmedianKmerOccurrences" << std::endl;


  std::unordered_map<std::string, uint32_t> kmerCounts;
  std::unordered_map<std::string, std::unordered_set<uint32_t>> kmerSampleCounts;

  auto seqs = createKmerReadVec(setUp.pars_.ioOptions_, kmerLength, false);

  {
    for(const auto & seq : iter::enumerate(seqs)) {
      for(const auto & k : seq.element->kInfo_.kmers_) {
        kmerCounts[k.second.k_] += k.second.positions_.size();
        kmerSampleCounts[k.second.k_].emplace(seq.index);
      }
    }
  }

  std::unordered_map<std::string, uint32_t> kmerCountsFinal;
  std::unordered_map<std::string, std::unordered_set<uint32_t>> kmerSampleCountsFinal;

  for(const auto & seqEnum : iter::enumerate(seqs)) {
    const auto & seq = seqEnum.element;
    std::vector<uint32_t> kmerOccurrences;
    uint32_t uniqueKmers = 0;
    for(const auto & k : seq->kInfo_.kmers_) {
      if(kmerSampleCounts[k.second.k_].size() == 1) {
        ++uniqueKmers;
      } else {
        kmerOccurrences.emplace_back(kmerSampleCounts[k.second.k_].size());
      }
    }
    outKmerInfoPerSeq << seq->seqBase_.name_
        << "\t" << len(seq->seqBase_)
        << "\t" << uniqueKmers
        << "\t" << uniqueKmers / static_cast<double>(len(getSeqBase(seq).seq_) - kmerLength + 1)
        << "\t" << vectorMean(kmerOccurrences)
        << "\t" << vectorMedianRef(kmerOccurrences);
    addOtherVec(kmerOccurrences, std::vector<uint32_t>(uniqueKmers, 1));
    outKmerInfoPerSeq << "\t" << vectorMean(kmerOccurrences)
        << "\t" << vectorMedianRef(kmerOccurrences);
    outKmerInfoPerSeq << std::endl;

  }

  std::map<uint32_t, uint32_t> countOfKmersWithSampleCounts;
  for(const auto & kCount : kmerSampleCounts) {
    ++countOfKmersWithSampleCounts[kCount.second.size()];
  }

  for (const auto &count: countOfKmersWithSampleCounts) {
    outKmerInfo << seqs.size() << "\t" << kmerLength
        << "\t" << count.first
        << "\t" << count.first / static_cast<double>(seqs.size())
        << "\t" << count.second
        << std::endl;
  }
  return 0;
}

int kmerSetExpRunner::createMinimallyNonRedundantDownSampledSet(const njh::progutils::CmdArgs & inputCommands) {
  uint32_t kmerLength = 15;
  uint32_t kmerMinOccurrenceCoverage = 3;
  double fractionUniqueToInclude = std::numeric_limits<double>::max();
  // OutOptions outOpts(bfs::path(""), ".tab.txt.gz");
  SeqIOOptions initialSeqOpts;
  seqSetUp setUp(inputCommands);
  setUp.description_ = "Down sample a set of sequences to optimize a set with minimally non redundant kmers";
  setUp.processVerbose();
  setUp.processDebug();
  setUp.setOption(fractionUniqueToInclude, "--fractionUniqueToInclude", "If a seq has this fraction of its kmers unique then include automatically");
  setUp.setOption(kmerMinOccurrenceCoverage, "--kmerMinOccurrenceCoverage", "kmer Min Occurrence Coverage");
  setUp.setOption(kmerLength, "--kmerLength", "kmer Length", true);
  setUp.processDefaultReader(true);
  bool initialSeqOptsSet = setUp.processSeqIoFilename(initialSeqOpts,"initialSeqs", false);
  // setUp.processWritingOptions(outOpts);
  setUp.finishSetUp(std::cout);

  std::vector<std::unique_ptr<seqWithKmerInfo> > initialSeqs;
  if (initialSeqOptsSet) {
    initialSeqs = createKmerReadVec(initialSeqOpts, kmerLength, false);
  }

  SeqOutput writer(setUp.pars_.ioOptions_);
  writer.openOut();

  // std::unordered_map<std::string, uint32_t> kmerCounts;
  std::unordered_map<std::string, std::unordered_set<uint32_t>> kmerSampleCounts;
  // std::unordered_map<char, std::unordered_map<std::string, std::unordered_set<uint32_t>>> kmerSampleCounts;

  auto seqs = createKmerReadVec(setUp.pars_.ioOptions_, kmerLength, false);

  {
    for(const auto & seq : iter::enumerate(seqs)) {
      for(const auto & k : seq.element->kInfo_.kmers_) {
        // kmerCounts[k.second.k_] += k.second.positions_.size();
        kmerSampleCounts[k.second.k_].emplace(seq.index);
      }
    }
  }


  // std::unordered_map<std::string, uint32_t> kmerCountsFinal;
  std::unordered_map<std::string, std::unordered_set<uint32_t>> kmerSampleCountsFinal;
  // std::unordered_map<char, std::unordered_map<std::string, std::unordered_set<uint32_t>>> kmerSampleCountsFinal;
  std::vector<size_t> outputIndexs;

  if(!initialSeqs.empty()) {
    for(const auto & seq : iter::enumerate(initialSeqs)) {
      for( auto & k : seq.element->kInfo_.kmers_) {
        // kmerCounts[k.second.k_] += k.second.positions_.size();
        kmerSampleCounts[k.second.k_].emplace(seq.index);
        k.second.readCnt_ = kmerSampleCounts[k.second.k_].size();
        // kmerCountsFinal[k.second.k_] += k.second.positions_.size();
        kmerSampleCountsFinal[k.second.k_].emplace(seq.index);
      }
    }
  }
  for(auto & seq : seqs) {
    for( auto & k : seq->kInfo_.kmers_) {
      k.second.readCnt_ = kmerSampleCounts[k.second.k_].size();
    }
  }


  double kmersAboveCutOffCoveredFrac = 0;
  uint32_t totalNumberOfKmersAboveCutOff = 0;
  for(const auto & kCounts : kmerSampleCounts) {
    if(kCounts.second.size() >= kmerMinOccurrenceCoverage) {
      ++totalNumberOfKmersAboveCutOff;
    }
  }

  if(fractionUniqueToInclude <= 1) {
    for(const auto & seqEnum : iter::enumerate(seqs)) {
      if(njh::notIn(seqEnum.index, outputIndexs)) {
        const auto & seq = seqEnum.element;
        uint32_t uniqueKmers = 0;
        for(const auto & k : seq->kInfo_.kmers_) {
          if(kmerSampleCounts[k.second.k_].size() == 1) {
            ++uniqueKmers;
          }
        }
        if( uniqueKmers / static_cast<double>(len(getSeqBase(seq).seq_) - kmerLength + 1) >= fractionUniqueToInclude) {
          outputIndexs.emplace_back(seqEnum.index);
          for(const auto & k : seqs[seqEnum.index]->kInfo_.kmers_) {
            // kmerCountsFinal[k.second.k_] += k.second.positions_.size();
            kmerSampleCountsFinal[k.second.k_].emplace(seqEnum.index);
            if(k.second.readCnt_ >= kmerMinOccurrenceCoverage) {
              for(const auto & idx : kmerSampleCounts[k.second.k_]) {
                seqs[idx]->kInfo_.kmers_[k.second.k_].on_ = false;
              }
            }
          }
        }
      }
    }
    uint32_t numberOfKmersAboveCutOffCovered = 0;
    for(const auto & kCounts : kmerSampleCounts) {
      if(kCounts.second.size() >= kmerMinOccurrenceCoverage) {
        if(njh::in(kCounts.first, kmerSampleCountsFinal)) {
          ++numberOfKmersAboveCutOffCovered;
        }
      }
    }
    kmersAboveCutOffCoveredFrac = static_cast<double>(numberOfKmersAboveCutOffCovered)/totalNumberOfKmersAboveCutOff;
    if(setUp.pars_.debug_) {
      std::cout << "after doing initial uniqueness cut off" << std::endl;
      std::cout <<"number of output seqs: " << outputIndexs.size() + initialSeqs.size() << "/" << seqs.size()+ initialSeqs.size() << " : "
      << static_cast<double>(outputIndexs.size() + initialSeqs.size())/static_cast<double>(seqs.size()+ initialSeqs.size()) * 100 << "%"<< std::endl;
      std::cout << "\tnumberOfKmersAboveCutOffCovered: " << numberOfKmersAboveCutOffCovered << "/" << totalNumberOfKmersAboveCutOff << std::endl;
      std::cout << "\tkmersAboveCutOffCoveredFrac: " << kmersAboveCutOffCoveredFrac << std::endl;
      std::cout << std::endl;
    }
  }


  /*
  *            if(k.second.readCnt_ >= kmerMinOccurrenceCoverage) {
              for(const auto & allSeqEnum : iter::enumerate(seqs)) {
                if(njh::notIn(allSeqEnum.index, outputIndexs)) {
                  // njh::stopWatch subwatch;
                  // subwatch.setLapName("get occurrences");
                  const auto & allSeq = allSeqEnum.element;
                  if(njh::in(k.second.k_, allSeq->kInfo_.kmers_)) {
                    allSeq->kInfo_.kmers_[k.second.k_].on_ = false;
                  }
                }
              }
            }
   **/

  while(kmersAboveCutOffCoveredFrac < 1) {
    njh::stopWatch watch;
    watch.setLapName("start");
    uint32_t bestIndex = std::numeric_limits<uint32_t>::max();
    // double highestMedianOccurrence = 0;
    // double highestMaxOccurrence = 0;
    uint32_t highestSumOccurrence = 0;
    watch.startNewLap("determine kmer occurrences");
    for(const auto & seqEnum : iter::enumerate(seqs)) {
      if(njh::notIn(seqEnum.index, outputIndexs)) {
        // njh::stopWatch subwatch;
        // subwatch.setLapName("get occurrences");
        const auto & seq = seqEnum.element;
        // std::vector<uint32_t> kmerOccurrences;
        uint32_t kmerOccurrencesSum = 0;
        // uint32_t uniqueKmers = 0;
        for(const auto & k : seq->kInfo_.kmers_) {
          if(k.second.on_) {
          // if(njh::notIn(k.second.k_, kmerSampleCountsFinal)) {
            // if(kmerSampleCounts[k.second.k_].size() == 1) {
            //   ++uniqueKmers;
            // } else {
            //   kmerOccurrences.emplace_back(kmerSampleCounts[k.second.k_].size());
            // }
            // auto kmerSampleCount = kmerSampleCounts[k.second.k_].size();
            // if(kmerSampleCount >= kmerMinOccurrenceCoverage) {
            //   kmerOccurrences.emplace_back(kmerSampleCount);
            // }
            if(k.second.readCnt_ >= kmerMinOccurrenceCoverage) {
              // kmerOccurrences.emplace_back(k.second.readCnt_);
              kmerOccurrencesSum += k.second.readCnt_;
            }
          }
        }
        // addOtherVec(kmerOccurrences, std::vector<uint32_t>(uniqueKmers, 1));
        // if(vectorMaximum(kmerOccurrences) > highestMaxOccurrence) {
        //   highestMaxOccurrence = vectorMaximum(kmerOccurrences);
        //   bestIndex = seqEnum.index;
        // }
        // std::cout << seq->seqBase_.name_ << "\n\t" << njh::conToStr(kmerOccurrences, ",") << std::endl;
        // subwatch.startNewLap("compute median");
        // if(vectorMedianRef(kmerOccurrences) > highestMedianOccurrence) {
        //   highestMedianOccurrence = vectorMedianRef(kmerOccurrences);
        //   bestIndex = seqEnum.index;
        //   if(njh::in(seqEnum.element->seqBase_.name_, VecStr{"PfCD01_040018100.1_protein_t1", "PfKE01_040019700.1_protein_t1"})) {
        //     std::cout << "adding: " << seqEnum.element->seqBase_.name_ << std::endl;
        //     std::cout << "kmerOccurrences.size() : " << kmerOccurrences.size() << std::endl;
        //     std::cout << "kmerOccurrences: " << njh::conToStr(kmerOccurrences, ",") << std::endl;
        //   } //
        // }
        if(kmerOccurrencesSum > highestSumOccurrence) {
          highestSumOccurrence = kmerOccurrencesSum;
          bestIndex = seqEnum.index;
        }
        // if(vectorMean(kmerOccurrences) > highestMedianOccurrence) {
        //   highestMedianOccurrence = vectorMean(kmerOccurrences);
        //   bestIndex = seqEnum.index;
        // }
        // if(setUp.pars_.debug_) {
        //   subwatch.logLapTimesOrganized(std::cout);
        // }
      }
    }
    watch.startNewLap("update final kmer maps");
    outputIndexs.emplace_back(bestIndex);
    for(const auto & k : seqs[bestIndex]->kInfo_.kmers_) {
      // kmerCountsFinal[k.second.k_] += k.second.positions_.size();
      kmerSampleCountsFinal[k.second.k_].emplace(bestIndex);
      if(k.second.readCnt_ >= kmerMinOccurrenceCoverage) {
        for(const auto & idx : kmerSampleCounts[k.second.k_]) {
          seqs[idx]->kInfo_.kmers_[k.second.k_].on_ = false;
        }
      }
    }
    watch.startNewLap("get number of kmers covered");
    uint32_t numberOfKmersAboveCutOffCovered = 0;
    for(const auto & kCounts : kmerSampleCounts) {
      if(kCounts.second.size() >= kmerMinOccurrenceCoverage) {
        if(njh::in(kCounts.first, kmerSampleCountsFinal)) {
          ++numberOfKmersAboveCutOffCovered;
        }
      }
    }
    kmersAboveCutOffCoveredFrac = static_cast<double>(numberOfKmersAboveCutOffCovered)/totalNumberOfKmersAboveCutOff;
    if(setUp.pars_.debug_) {
      std::cout <<"number of output seqs: " << outputIndexs.size() + initialSeqs.size() << "/" << seqs.size()+ initialSeqs.size() << " : "
      << static_cast<double>(outputIndexs.size() + initialSeqs.size())/static_cast<double>(seqs.size()+ initialSeqs.size()) * 100 << "%"<< std::endl;
      std::cout << "\tnumberOfKmersAboveCutOffCovered: " << numberOfKmersAboveCutOffCovered << "/" << totalNumberOfKmersAboveCutOff << std::endl;
      std::cout << "\tkmersAboveCutOffCoveredFrac: " << kmersAboveCutOffCoveredFrac << std::endl;
      watch.logLapTimesOrganized(std::cout);
      std::cout << std::endl;
    }
  }

  for(const auto & seq : initialSeqs) {
    writer.write(seq);
  }
  for(const auto & idx : outputIndexs) {
    writer.write(seqs[idx]);
  }

  return 0;
}




} //namespace njhseq