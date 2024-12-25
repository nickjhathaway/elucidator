//
// Created by Nicholas Hathaway on 12/24/24.
//
#include "seqUtilsExtractRunner.hpp"

#include <njhseq/IO/SeqIO.h>
#include <njhseq/objects/helperObjects/motif.hpp>


namespace njhseq {

struct MotifError {
  MotifError(const std::string & motifStr, uint32_t allowableError ): motifStr_(motifStr), allowableError_(allowableError) {}
  std::string motifStr_;
  uint32_t allowableError_{0};

  uint32_t passingScore_ = std::numeric_limits<uint32_t>::max();
  std::shared_ptr<motif> motif_;

  void setMotif() {
    motif_  = std::make_shared<motif>(motifStr_);
    if (allowableError_ >= motif_->size()) {
      std::stringstream ss;
      ss << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << ", error for motif: " << motifStr_ << "motifAllowableError: " << allowableError_ << " can't be equal or greater than filter_motif.size(): " << motif_->size() << "\n";
      throw std::runtime_error{ss.str()};
    }
    passingScore_ = motif_->size() - allowableError_;
  }
};

int seqUtilsExtractRunner::filterOffSeqsWithEndingMotif(const njh::progutils::CmdArgs & inputCommands) {


  OutOptions outOpts("", ".tsv");

  std::vector<MotifError> motifErrors{
    MotifError("GGGGGGGGGGGGGG", 1),
    MotifError("GGGGGGGGGGGGGGGGGGGG", 3),
    MotifError("CCCCCCCCCCCCCC", 1),
    MotifError("CCCCCCCCCCCCCCCCCCCC", 3),
  };

  bool writeOutFiltered = false;

  bool needBothPairs = false;

  uint32_t numOfThreads = 1;
  seqUtilsExtractSetUp setUp(inputCommands);
  setUp.processVerbose();
  setUp.processDebug();

  bfs::path fastq1Fnp = "";
  bfs::path fastq2Fnp = "";
  bfs::path fastqFnp = "";
  bool setFastq1 = setUp.setOption(fastq1Fnp, "--fastq1,--fastq1gz", "Fastq first mate File");
  setUp.setOption(fastq2Fnp, "--fastq2,--fastq2gz", "Fastq second mate File", setFastq1);
  bool revCompMate = false;
  setUp.setOption(revCompMate, "--revCompMate", "Reverse Complement Sequences in mate file");
  setUp.setOption(fastqFnp, "--fastq,--fastqgz", "Fastq File", !setFastq1);
  bool overWrite = false;
  setUp.setOption(overWrite, "--overWrite", "Overwrite output files");
  if (setFastq1) {
    setUp.pars_.ioOptions_.firstName_ = fastq1Fnp;
    setUp.pars_.ioOptions_.revComplMate_ = revCompMate;
    setUp.pars_.ioOptions_.inFormat_ = SeqIOOptions::inFormats::FASTQPAIRED;
    setUp.pars_.ioOptions_.outFormat_ = SeqIOOptions::outFormats::FASTQPAIREDGZ;
  } else {
    setUp.pars_.ioOptions_.firstName_ = fastqFnp;
    setUp.pars_.ioOptions_.inFormat_ = SeqIOOptions::inFormats::FASTQ;
    setUp.pars_.ioOptions_.outFormat_ = SeqIOOptions::outFormats::FASTQGZ;
  }
  setUp.pars_.ioOptions_.out_.overWriteFile_ = overWrite;

  auto pairedIn = SeqIOOptions::genPairedIn(fastq1Fnp, fastq2Fnp);
  pairedIn.revComplMate_ = revCompMate;
  pairedIn.out_.overWriteFile_ = overWrite;

  auto singleIn = SeqIOOptions::genFastqIn(fastqFnp);
  singleIn.out_.overWriteFile_ = overWrite;

  if (!setUp.pars_.ioOptions_.firstName_.empty()) {
    outOpts.outFilename_ = njh::files::prependFileBasename(njh::files::removeExtension(setUp.pars_.ioOptions_.firstName_), "filteredCounts_");
  }

  std::string outputStub = "out";

  bfs::path motifTableFnp;
  setUp.setOption(motifTableFnp, "--motifTable", "a table with first column motif, second column amount of allowable error, this will replace the default motifs of Gs and Cs");

  setUp.setOption(outputStub, "--outputStub", "output stub for writing out the kept reads");
  setUp.setOption(outOpts.outFilename_, "--filteredCountsOutputFnp", "filtered Counts Output file");

  setUp.setOption(writeOutFiltered, "--writeOutFiltered", "write Out Filtered");
  setUp.setOption(needBothPairs, "--needBothPairs", "need both pairs to have pattern, by default requires that just one has the pattern");
  setUp.setOption(numOfThreads, "--numOfThreads", "number of threads");

  setUp.finishSetUp(std::cout);


  if ("" != motifTableFnp) {
    table motifTab(motifTableFnp, "\t", false);
    if (motifTab.nCol() != 2) {
      std::stringstream ss;
      ss << __PRETTY_FUNCTION__ << ", error " << "error, " << motifTableFnp << " should have 2 columns" << "\n";
      throw std::runtime_error{ss.str()};
    }
    motifErrors.clear();
    for (const auto & row : motifTab) {
      motifErrors.emplace_back(row[0], njh::StrToNumConverter::stoToNum<uint32_t>(row[1]));
    }
  }
  for (auto & m : motifErrors) {
    m.setMotif();

  }

  outOpts.transferOverwriteOpts(setUp.pars_.ioOptions_.out_);

  OutputStream out(outOpts);

  MultiSeqIO multi_seq_io;
  multi_seq_io.addReader("single", SeqIOOptions::genFastqOutGz(outputStub));
  multi_seq_io.addReader("pairs", SeqIOOptions::genPairedOutGz(outputStub));
  if (writeOutFiltered) {
    multi_seq_io.addReader("single-filtered", SeqIOOptions::genFastqOutGz(njh::files::prependFileBasename(outputStub, "filtered_")));
    multi_seq_io.addReader("pairs-filtered", SeqIOOptions::genPairedOutGz(njh::files::prependFileBasename(outputStub, "filtered_")));
  }
  if  (overWrite) {
    multi_seq_io.setAllReaderToOverwrite();
  }

  out << "inputFnp\ttotalReads\ttotalFiltered\tmotif\tmotifAllowError\tfilteredForMotif\n";

  if (!fastq1Fnp.empty() && bfs::exists(fastq1Fnp)) {


    uint32_t totalReads = 0;
    uint32_t totalFiltered = 0;
    std::unordered_map<uint32_t, uint32_t> filteredCounts;
    for (const auto & e : iter::enumerate(motifErrors)) {
      filteredCounts[e.index] = 0;
    }
    std::mutex countMut;
    SeqInput reader(pairedIn);
    reader.openIn();

    std::function<void()> filteredSeqs = [&reader,&multi_seq_io,
          &filteredCounts,
          &motifErrors,
          &writeOutFiltered,
          &needBothPairs,
          &totalReads, &totalFiltered, &countMut]() {
      uint32_t current_totalReads = 0;
      uint32_t current_totalFiltered = 0;
      std::unordered_map<uint32_t, uint32_t> current_filteredCounts;
      PairedRead pseq;
      while (reader.readNextReadLock(pseq)) {
        ++current_totalReads;
        bool pass = true;
        for (const auto &e: iter::enumerate(motifErrors)) {
          bool firstMateHasMotif = pseq.seqBase_.seq_.size() > e.element.motif_->size() &&
                                   e.element.motif_->passMotifParameter(
                                     pseq.seqBase_.seq_.begin() + (pseq.seqBase_.seq_.size() - e.element.motif_->size()),
                                     pseq.seqBase_.seq_.end(), e.element.passingScore_ );
          bool secondMateHasMotif = pseq.mateSeqBase_.seq_.size() > e.element.motif_->size() &&
                                    e.element.motif_->passMotifParameter(
                                      pseq.mateSeqBase_.seq_.begin() + (
                                        pseq.mateSeqBase_.seq_.size() - e.element.motif_->size()),
                                      pseq.mateSeqBase_.seq_.end(), e.element.passingScore_ );
          if ((needBothPairs && firstMateHasMotif && secondMateHasMotif) || (
                !needBothPairs && (firstMateHasMotif || secondMateHasMotif))) {
            ++current_totalFiltered;
            current_filteredCounts[e.index] += 1;
            pass = false;
            break;
          }
        }
        if (pass) {
          multi_seq_io.openWrite("pairs",pseq);
        } else if (writeOutFiltered) {
          multi_seq_io.openWrite("pairs-filtered", pseq);
        }
      }
      {
        std::lock_guard<std::mutex> lock(countMut);
        totalFiltered += current_totalFiltered;
        totalReads += current_totalReads;
        for (const auto & count : current_filteredCounts) {
          filteredCounts[count.first] += count.second;
        }
      }
    };
    njh::concurrent::runVoidFunctionThreaded(filteredSeqs, numOfThreads);
    for (const auto & motifError : iter::enumerate(motifErrors)) {
      out << fastq1Fnp
      << "\t" << totalReads
      << "\t" << totalFiltered
      << "\t" << motifError.element.motifStr_
      << "\t" << motifError.element.allowableError_
      << "\t" << filteredCounts[motifError.index] << std::endl;
    }
  }
  if (!fastqFnp.empty() && bfs::exists(fastqFnp)) {

    uint32_t totalReads = 0;
    uint32_t totalFiltered = 0;
    std::unordered_map<uint32_t, uint32_t> filteredCounts;
    for (const auto & e : iter::enumerate(motifErrors)) {
      filteredCounts[e.index] = 0;
    }
    std::mutex countMut;
    SeqInput reader(singleIn);
    reader.openIn();

    std::function<void()> filteredSeqs = [&reader,&multi_seq_io,
          &filteredCounts,
          &motifErrors,
          &writeOutFiltered,
          &totalReads, &totalFiltered, &countMut]() {
      uint32_t current_totalReads = 0;
      uint32_t current_totalFiltered = 0;
      std::unordered_map<uint32_t, uint32_t> current_filteredCounts;
      seqInfo seq;
      while (reader.readNextReadLock(seq)) {
        ++current_totalReads;
        bool pass = true;
        for (const auto &e: iter::enumerate(motifErrors)) {
          if (seq.seq_.size() > e.element.motif_->size() &&
                                   e.element.motif_->passMotifParameter(
                                     seq.seq_.begin() + (seq.seq_.size() - e.element.motif_->size()),
                                     seq.seq_.end(), e.element.passingScore_ )) {
            ++current_totalFiltered;
            current_filteredCounts[e.index] += 1;
            pass = false;
            break;
          }
        }
        if (pass) {
          multi_seq_io.openWrite("single",seq);
        } else if (writeOutFiltered) {
          multi_seq_io.openWrite("single-filtered", seq);
        }
      }
      {
        std::lock_guard<std::mutex> lock(countMut);
        totalFiltered += current_totalFiltered;
        totalReads += current_totalReads;
        for (const auto & count : current_filteredCounts) {
          filteredCounts[count.first] += count.second;
        }
      }
    };
    njh::concurrent::runVoidFunctionThreaded(filteredSeqs, numOfThreads);
    for (const auto & motifError : iter::enumerate(motifErrors)) {
      out << fastqFnp
      << "\t" << totalReads
      << "\t" << totalFiltered
      << "\t" << motifError.element.motifStr_
      << "\t" << motifError.element.allowableError_
      << "\t" << filteredCounts[motifError.index] << std::endl;
    }
  }
  return 0;
}


int seqUtilsExtractRunner::filterOffSeqsWithBeginningMotif(const njh::progutils::CmdArgs & inputCommands) {



  OutOptions outOpts("", ".tsv");

  std::vector<MotifError> motifErrors{
    MotifError("GGGGGGGGGGGGGG", 1),
    MotifError("GGGGGGGGGGGGGGGGGGGG", 3),
    MotifError("CCCCCCCCCCCCCC", 1),
    MotifError("CCCCCCCCCCCCCCCCCCCC", 3),
  };

  bool writeOutFiltered = false;

  bool needBothPairs = false;

  uint32_t numOfThreads = 1;
  seqUtilsExtractSetUp setUp(inputCommands);
  setUp.processVerbose();
  setUp.processDebug();

  bfs::path fastq1Fnp = "";
  bfs::path fastq2Fnp = "";
  bfs::path fastqFnp = "";
  bool setFastq1 = setUp.setOption(fastq1Fnp, "--fastq1,--fastq1gz", "Fastq first mate File");
  setUp.setOption(fastq2Fnp, "--fastq2,--fastq2gz", "Fastq second mate File", setFastq1);
  bool revCompMate = false;
  setUp.setOption(revCompMate, "--revCompMate", "Reverse Complement Sequences in mate file");
  setUp.setOption(fastqFnp, "--fastq,--fastqgz", "Fastq File", !setFastq1);
  bool overWrite = false;
  setUp.setOption(overWrite, "--overWrite", "Overwrite output files");
  if (setFastq1) {
    setUp.pars_.ioOptions_.firstName_ = fastq1Fnp;
    setUp.pars_.ioOptions_.revComplMate_ = revCompMate;
    setUp.pars_.ioOptions_.inFormat_ = SeqIOOptions::inFormats::FASTQPAIRED;
    setUp.pars_.ioOptions_.outFormat_ = SeqIOOptions::outFormats::FASTQPAIREDGZ;
  } else {
    setUp.pars_.ioOptions_.firstName_ = fastqFnp;
    setUp.pars_.ioOptions_.inFormat_ = SeqIOOptions::inFormats::FASTQ;
    setUp.pars_.ioOptions_.outFormat_ = SeqIOOptions::outFormats::FASTQGZ;
  }
  setUp.pars_.ioOptions_.out_.overWriteFile_ = overWrite;

  auto pairedIn = SeqIOOptions::genPairedIn(fastq1Fnp, fastq2Fnp);
  pairedIn.revComplMate_ = revCompMate;
  pairedIn.out_.overWriteFile_ = overWrite;

  auto singleIn = SeqIOOptions::genFastqIn(fastqFnp);
  singleIn.out_.overWriteFile_ = overWrite;

  if (!setUp.pars_.ioOptions_.firstName_.empty()) {
    outOpts.outFilename_ = njh::files::prependFileBasename(njh::files::removeExtension(setUp.pars_.ioOptions_.firstName_), "filteredCounts_");
  }

  std::string outputStub = "out";

  bfs::path motifTableFnp;
  setUp.setOption(motifTableFnp, "--motifTable", "a table with first column motif, second column amount of allowable error, this will replace the default motifs of Gs and Cs");

  setUp.setOption(outputStub, "--outputStub", "output stub for writing out the kept reads");
  setUp.setOption(outOpts.outFilename_, "--filteredCountsOutputFnp", "filtered Counts Output file");

  setUp.setOption(writeOutFiltered, "--writeOutFiltered", "write Out Filtered");
  setUp.setOption(needBothPairs, "--needBothPairs", "need both pairs to have pattern, by default requires that just one has the pattern");
  setUp.setOption(numOfThreads, "--numOfThreads", "number of threads");

  setUp.finishSetUp(std::cout);


  if ("" != motifTableFnp) {
    table motifTab(motifTableFnp, "\t", false);
    if (motifTab.nCol() != 2) {
      std::stringstream ss;
      ss << __PRETTY_FUNCTION__ << ", error " << "error, " << motifTableFnp << " should have 2 columns" << "\n";
      throw std::runtime_error{ss.str()};
    }
    motifErrors.clear();
    for (const auto & row : motifTab) {
      motifErrors.emplace_back(row[0], njh::StrToNumConverter::stoToNum<uint32_t>(row[1]));
    }
  }
  for (auto & m : motifErrors) {
    m.setMotif();

  }

  outOpts.transferOverwriteOpts(setUp.pars_.ioOptions_.out_);

  OutputStream out(outOpts);

  MultiSeqIO multi_seq_io;
  multi_seq_io.addReader("single", SeqIOOptions::genFastqOutGz(outputStub));
  multi_seq_io.addReader("pairs", SeqIOOptions::genPairedOutGz(outputStub));
  if (writeOutFiltered) {
    multi_seq_io.addReader("single-filtered", SeqIOOptions::genFastqOutGz(njh::files::prependFileBasename(outputStub, "filtered_")));
    multi_seq_io.addReader("pairs-filtered", SeqIOOptions::genPairedOutGz(njh::files::prependFileBasename(outputStub, "filtered_")));
  }
  if  (overWrite) {
    multi_seq_io.setAllReaderToOverwrite();
  }

  out << "inputFnp\ttotalReads\ttotalFiltered\tmotif\tmotifAllowError\tfilteredForMotif\n";

  if (!fastq1Fnp.empty() && bfs::exists(fastq1Fnp)) {


    uint32_t totalReads = 0;
    uint32_t totalFiltered = 0;
    std::unordered_map<uint32_t, uint32_t> filteredCounts;
    for (const auto & e : iter::enumerate(motifErrors)) {
      filteredCounts[e.index] = 0;
    }
    std::mutex countMut;
    SeqInput reader(pairedIn);
    reader.openIn();

    std::function<void()> filteredSeqs = [&reader,&multi_seq_io,
          &filteredCounts,
          &motifErrors,
          &writeOutFiltered,
          &needBothPairs,
          &totalReads, &totalFiltered, &countMut]() {
      uint32_t current_totalReads = 0;
      uint32_t current_totalFiltered = 0;
      std::unordered_map<uint32_t, uint32_t> current_filteredCounts;
      PairedRead pseq;
      while (reader.readNextReadLock(pseq)) {
        ++current_totalReads;
        bool pass = true;
        for (const auto &e: iter::enumerate(motifErrors)) {
          bool firstMateHasMotif = pseq.seqBase_.seq_.size() > e.element.motif_->size() &&
                                   e.element.motif_->passMotifParameter(
                                     pseq.seqBase_.seq_.begin(),
                                     pseq.seqBase_.seq_.begin() + e.element.motif_->size(),
                                     e.element.passingScore_);
          bool secondMateHasMotif = pseq.mateSeqBase_.seq_.size() > e.element.motif_->size() &&
                                    e.element.motif_->passMotifParameter(
                                      pseq.mateSeqBase_.seq_.begin(),
                                      pseq.mateSeqBase_.seq_.begin() + e.element.motif_->size(),
                                      e.element.passingScore_);
          if ((needBothPairs && firstMateHasMotif && secondMateHasMotif) || (
                !needBothPairs && (firstMateHasMotif || secondMateHasMotif))) {
            ++current_totalFiltered;
            current_filteredCounts[e.index] += 1;
            pass = false;
            break;
          }
        }
        if (pass) {
          multi_seq_io.openWrite("pairs",pseq);
        } else if (writeOutFiltered) {
          multi_seq_io.openWrite("pairs-filtered", pseq);
        }
      }
      {
        std::lock_guard<std::mutex> lock(countMut);
        totalFiltered += current_totalFiltered;
        totalReads += current_totalReads;
        for (const auto & count : current_filteredCounts) {
          filteredCounts[count.first] += count.second;
        }
      }
    };
    njh::concurrent::runVoidFunctionThreaded(filteredSeqs, numOfThreads);
    for (const auto & motifError : iter::enumerate(motifErrors)) {
      out << fastq1Fnp
      << "\t" << totalReads
      << "\t" << totalFiltered
      << "\t" << motifError.element.motifStr_
      << "\t" << motifError.element.allowableError_
      << "\t" << filteredCounts[motifError.index] << std::endl;
    }
  }
  if (!fastqFnp.empty() && bfs::exists(fastqFnp)) {

    uint32_t totalReads = 0;
    uint32_t totalFiltered = 0;
    std::unordered_map<uint32_t, uint32_t> filteredCounts;
    for (const auto & e : iter::enumerate(motifErrors)) {
      filteredCounts[e.index] = 0;
    }
    std::mutex countMut;
    SeqInput reader(singleIn);
    reader.openIn();

    std::function<void()> filteredSeqs = [&reader,&multi_seq_io,
          &filteredCounts,
          &motifErrors,
          &writeOutFiltered,
          &totalReads, &totalFiltered, &countMut]() {
      uint32_t current_totalReads = 0;
      uint32_t current_totalFiltered = 0;
      std::unordered_map<uint32_t, uint32_t> current_filteredCounts;
      seqInfo seq;
      while (reader.readNextReadLock(seq)) {
        ++current_totalReads;
        bool pass = true;
        for (const auto &e: iter::enumerate(motifErrors)) {
          if (seq.seq_.size() > e.element.motif_->size() &&
              e.element.motif_->passMotifParameter(seq.seq_.begin(),
                                                   seq.seq_.begin() + e.element.motif_->size(),
                                                   e.element.passingScore_)) {
            ++current_totalFiltered;
            current_filteredCounts[e.index] += 1;
            pass = false;
            break;
          }
        }
        if (pass) {
          multi_seq_io.openWrite("single",seq);
        } else if (writeOutFiltered) {
          multi_seq_io.openWrite("single-filtered", seq);
        }
      }
      {
        std::lock_guard<std::mutex> lock(countMut);
        totalFiltered += current_totalFiltered;
        totalReads += current_totalReads;
        for (const auto & count : current_filteredCounts) {
          filteredCounts[count.first] += count.second;
        }
      }
    };
    njh::concurrent::runVoidFunctionThreaded(filteredSeqs, numOfThreads);
    for (const auto & motifError : iter::enumerate(motifErrors)) {
      out << fastqFnp
      << "\t" << totalReads
      << "\t" << totalFiltered
      << "\t" << motifError.element.motifStr_
      << "\t" << motifError.element.allowableError_
      << "\t" << filteredCounts[motifError.index] << std::endl;
    }
  }
  return 0;
}

int seqUtilsExtractRunner::filterOffSeqsWithEdgeMotif(const njh::progutils::CmdArgs & inputCommands) {

  OutOptions outOpts("", ".tsv");

  std::vector<MotifError> motifErrors{
    MotifError("GGGGGGGGGGGGGG", 1),
    MotifError("GGGGGGGGGGGGGGGGGGGG", 3),
    MotifError("CCCCCCCCCCCCCC", 1),
    MotifError("CCCCCCCCCCCCCCCCCCCC", 3),
  };

  bool writeOutFiltered = false;

  bool needBothPairs = false;

  uint32_t numOfThreads = 1;
  seqUtilsExtractSetUp setUp(inputCommands);
  setUp.processVerbose();
  setUp.processDebug();

  bfs::path fastq1Fnp = "";
  bfs::path fastq2Fnp = "";
  bfs::path fastqFnp = "";
  bool setFastq1 = setUp.setOption(fastq1Fnp, "--fastq1,--fastq1gz", "Fastq first mate File");
  setUp.setOption(fastq2Fnp, "--fastq2,--fastq2gz", "Fastq second mate File", setFastq1);
  bool revCompMate = false;
  setUp.setOption(revCompMate, "--revCompMate", "Reverse Complement Sequences in mate file");
  setUp.setOption(fastqFnp, "--fastq,--fastqgz", "Fastq File", !setFastq1);
  bool overWrite = false;
  setUp.setOption(overWrite, "--overWrite", "Overwrite output files");
  if (setFastq1) {
    setUp.pars_.ioOptions_.firstName_ = fastq1Fnp;
    setUp.pars_.ioOptions_.revComplMate_ = revCompMate;
    setUp.pars_.ioOptions_.inFormat_ = SeqIOOptions::inFormats::FASTQPAIRED;
    setUp.pars_.ioOptions_.outFormat_ = SeqIOOptions::outFormats::FASTQPAIREDGZ;
  } else {
    setUp.pars_.ioOptions_.firstName_ = fastqFnp;
    setUp.pars_.ioOptions_.inFormat_ = SeqIOOptions::inFormats::FASTQ;
    setUp.pars_.ioOptions_.outFormat_ = SeqIOOptions::outFormats::FASTQGZ;
  }
  setUp.pars_.ioOptions_.out_.overWriteFile_ = overWrite;

  auto pairedIn = SeqIOOptions::genPairedIn(fastq1Fnp, fastq2Fnp);
  pairedIn.revComplMate_ = revCompMate;
  pairedIn.out_.overWriteFile_ = overWrite;

  auto singleIn = SeqIOOptions::genFastqIn(fastqFnp);
  singleIn.out_.overWriteFile_ = overWrite;

  if (!setUp.pars_.ioOptions_.firstName_.empty()) {
    outOpts.outFilename_ = njh::files::prependFileBasename(njh::files::removeExtension(setUp.pars_.ioOptions_.firstName_), "filteredCounts_");
  }

  std::string outputStub = "out";

  bfs::path motifTableFnp;
  setUp.setOption(motifTableFnp, "--motifTable", "a table with first column motif, second column amount of allowable error, this will replace the default motifs of Gs and Cs");

  setUp.setOption(outputStub, "--outputStub", "output stub for writing out the kept reads");
  setUp.setOption(outOpts.outFilename_, "--filteredCountsOutputFnp", "filtered Counts Output file");

  setUp.setOption(writeOutFiltered, "--writeOutFiltered", "write Out Filtered");
  setUp.setOption(needBothPairs, "--needBothPairs", "need both pairs to have pattern, by default requires that just one has the pattern");
  setUp.setOption(numOfThreads, "--numOfThreads", "number of threads");

  setUp.finishSetUp(std::cout);


  if ("" != motifTableFnp) {
    table motifTab(motifTableFnp, "\t", false);
    if (motifTab.nCol() != 2) {
      std::stringstream ss;
      ss << __PRETTY_FUNCTION__ << ", error " << "error, " << motifTableFnp << " should have 2 columns" << "\n";
      throw std::runtime_error{ss.str()};
    }
    motifErrors.clear();
    for (const auto & row : motifTab) {
      motifErrors.emplace_back(row[0], njh::StrToNumConverter::stoToNum<uint32_t>(row[1]));
    }
  }
  for (auto & m : motifErrors) {
    m.setMotif();

  }

  outOpts.transferOverwriteOpts(setUp.pars_.ioOptions_.out_);

  OutputStream out(outOpts);

  MultiSeqIO multi_seq_io;
  multi_seq_io.addReader("single", SeqIOOptions::genFastqOutGz(outputStub));
  multi_seq_io.addReader("pairs", SeqIOOptions::genPairedOutGz(outputStub));
  if (writeOutFiltered) {
    multi_seq_io.addReader("single-filtered", SeqIOOptions::genFastqOutGz(njh::files::prependFileBasename(outputStub, "filtered_")));
    multi_seq_io.addReader("pairs-filtered", SeqIOOptions::genPairedOutGz(njh::files::prependFileBasename(outputStub, "filtered_")));
  }
  if  (overWrite) {
    multi_seq_io.setAllReaderToOverwrite();
  }

  out << "inputFnp\ttotalReads\ttotalFiltered\tmotif\tmotifAllowError\tedge\tfilteredForMotif\n";

  if (!fastq1Fnp.empty() && bfs::exists(fastq1Fnp)) {


    uint32_t totalReads = 0;
    uint32_t totalFiltered = 0;
    std::unordered_map<uint32_t, uint32_t> filteredCounts_end;
    std::unordered_map<uint32_t, uint32_t> filteredCounts_beg;
    for (const auto & e : iter::enumerate(motifErrors)) {
      filteredCounts_end[e.index] = 0;
      filteredCounts_beg[e.index] = 0;
    }
    std::mutex countMut;
    SeqInput reader(pairedIn);
    reader.openIn();

    std::function<void()> filteredSeqs = [&reader,&multi_seq_io,
          &filteredCounts_end, &filteredCounts_beg,
          &motifErrors,
          &writeOutFiltered,
          &needBothPairs,
          &totalReads, &totalFiltered, &countMut]() {
      uint32_t current_totalReads = 0;
      uint32_t current_totalFiltered = 0;
      std::unordered_map<uint32_t, uint32_t> current_filteredCounts_end;
      std::unordered_map<uint32_t, uint32_t> current_filteredCounts_beg;
      PairedRead pseq;
      while (reader.readNextReadLock(pseq)) {
        ++current_totalReads;
        bool pass = true;
        for (const auto &e: iter::enumerate(motifErrors)) {
          bool firstMateEndHasMotif = pseq.seqBase_.seq_.size() > e.element.motif_->size() &&
                                   e.element.motif_->passMotifParameter(
                                     pseq.seqBase_.seq_.begin() + (pseq.seqBase_.seq_.size() - e.element.motif_->size()),
                                     pseq.seqBase_.seq_.end(), e.element.passingScore_ );
          bool secondMateEndHasMotif = pseq.mateSeqBase_.seq_.size() > e.element.motif_->size() &&
                                    e.element.motif_->passMotifParameter(
                                      pseq.mateSeqBase_.seq_.begin() + (
                                        pseq.mateSeqBase_.seq_.size() - e.element.motif_->size()),
                                      pseq.mateSeqBase_.seq_.end(), e.element.passingScore_ );
          if ((needBothPairs && firstMateEndHasMotif && secondMateEndHasMotif) || (
                !needBothPairs && (firstMateEndHasMotif || secondMateEndHasMotif))) {
            ++current_totalFiltered;
            current_filteredCounts_end[e.index] += 1;
            pass = false;
            break;
          } else {
            bool firstMateBegHasMotif = pseq.seqBase_.seq_.size() > e.element.motif_->size() &&
                                        e.element.motif_->passMotifParameter(
                                          pseq.seqBase_.seq_.begin(),
                                          pseq.seqBase_.seq_.begin() + e.element.motif_->size(),
                                          e.element.passingScore_);
            bool secondMateBegHasMotif = pseq.mateSeqBase_.seq_.size() > e.element.motif_->size() &&
                                         e.element.motif_->passMotifParameter(
                                           pseq.mateSeqBase_.seq_.begin(),
                                           pseq.mateSeqBase_.seq_.begin() + e.element.motif_->size(),
                                           e.element.passingScore_);
            if ((needBothPairs && firstMateBegHasMotif && secondMateBegHasMotif) || (
                  !needBothPairs && (firstMateBegHasMotif || secondMateBegHasMotif))) {
              ++current_totalFiltered;
              current_filteredCounts_beg[e.index] += 1;
              pass = false;
              break;
            }
          }
        }
        if (pass) {
          multi_seq_io.openWrite("pairs",pseq);
        } else if (writeOutFiltered) {
          multi_seq_io.openWrite("pairs-filtered", pseq);
        }
      }
      {
        std::lock_guard<std::mutex> lock(countMut);
        totalFiltered += current_totalFiltered;
        totalReads += current_totalReads;
        for (const auto & count : current_filteredCounts_end) {
          filteredCounts_end[count.first] += count.second;
        }
        for (const auto & count : current_filteredCounts_beg) {
          filteredCounts_beg[count.first] += count.second;
        }
      }
    };
    njh::concurrent::runVoidFunctionThreaded(filteredSeqs, numOfThreads);
    for (const auto &motifError: iter::enumerate(motifErrors)) {
      out << fastq1Fnp
          << "\t" << totalReads
          << "\t" << totalFiltered
          << "\t" << motifError.element.motifStr_
          << "\t" << motifError.element.allowableError_
          << "\t" << "end"
          << "\t" << filteredCounts_end[motifError.index] << std::endl;
      out << fastq1Fnp
          << "\t" << totalReads
          << "\t" << totalFiltered
          << "\t" << motifError.element.motifStr_
          << "\t" << motifError.element.allowableError_
          << "\t" << "beginning"
          << "\t" << filteredCounts_beg[motifError.index] << std::endl;
    }
  }
  if (!fastqFnp.empty() && bfs::exists(fastqFnp)) {

    uint32_t totalReads = 0;
    uint32_t totalFiltered = 0;
    std::unordered_map<uint32_t, uint32_t> filteredCounts_end;
    std::unordered_map<uint32_t, uint32_t> filteredCounts_beg;
    for (const auto & e : iter::enumerate(motifErrors)) {
      filteredCounts_end[e.index] = 0;
      filteredCounts_beg[e.index] = 0;
    }
    std::mutex countMut;
    SeqInput reader(singleIn);
    reader.openIn();

    std::function<void()> filteredSeqs = [&reader,&multi_seq_io,
          &filteredCounts_end, &filteredCounts_beg,
          &motifErrors,
          &writeOutFiltered,
          &totalReads, &totalFiltered, &countMut]() {
      uint32_t current_totalReads = 0;
      uint32_t current_totalFiltered = 0;
      std::unordered_map<uint32_t, uint32_t> current_filteredCounts_end;
      std::unordered_map<uint32_t, uint32_t> current_filteredCounts_beg;
      seqInfo seq;
      while (reader.readNextReadLock(seq)) {
        ++current_totalReads;
        bool pass = true;
        for (const auto &e: iter::enumerate(motifErrors)) {
          if (seq.seq_.size() > e.element.motif_->size() &&
                                   e.element.motif_->passMotifParameter(
                                     seq.seq_.begin() + (seq.seq_.size() - e.element.motif_->size()),
                                     seq.seq_.end(), e.element.passingScore_ )) {
            ++current_totalFiltered;
            current_filteredCounts_end[e.index] += 1;
            pass = false;
            break;
          } else if (seq.seq_.size() > e.element.motif_->size() &&
                     e.element.motif_->passMotifParameter(seq.seq_.begin(),
                                                          seq.seq_.begin() + e.element.motif_->size(),
                                                          e.element.passingScore_)) {
            ++current_totalFiltered;
            current_filteredCounts_beg[e.index] += 1;
            pass = false;
            break;
          }
        }
        if (pass) {
          multi_seq_io.openWrite("single",seq);
        } else if (writeOutFiltered) {
          multi_seq_io.openWrite("single-filtered", seq);
        }
      }
      {
        std::lock_guard<std::mutex> lock(countMut);
        totalFiltered += current_totalFiltered;
        totalReads += current_totalReads;
        for (const auto & count : current_filteredCounts_end) {
          filteredCounts_end[count.first] += count.second;
        }
        for (const auto & count : current_filteredCounts_beg) {
          filteredCounts_beg[count.first] += count.second;
        }
      }
    };
    njh::concurrent::runVoidFunctionThreaded(filteredSeqs, numOfThreads);
    for (const auto &motifError: iter::enumerate(motifErrors)) {
      out << fastqFnp
          << "\t" << totalReads
          << "\t" << totalFiltered
          << "\t" << motifError.element.motifStr_
          << "\t" << motifError.element.allowableError_
          << "\t" << "end"
          << "\t" << filteredCounts_end[motifError.index] << std::endl;
      out << fastqFnp
          << "\t" << totalReads
          << "\t" << totalFiltered
          << "\t" << motifError.element.motifStr_
          << "\t" << motifError.element.allowableError_
          << "\t" << "beginning"
          << "\t" << filteredCounts_beg[motifError.index] << std::endl;
    }
  }
  return 0;
}






} //namespace njhseq

