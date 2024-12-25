//
// Created by Nicholas Hathaway on 12/24/24.
//
#include "seqUtilsExtractRunner.hpp"

#include <njhseq/IO/SeqIO.h>
#include <njhseq/objects/helperObjects/motif.hpp>


namespace njhseq {

int seqUtilsExtractRunner::filterOffSeqsWithEndingMotif(const njh::progutils::CmdArgs & inputCommands) {
  OutOptions outOpts("", ".tsv");
  std::string motifStr = "GGGGGGGGGGGGGG";
  uint32_t motifAllowableError = 1;

  std::string secondaryMotifStr = "GGGGGGGGGGGGGGGGGGGG";
  uint32_t secondaryMotifAllowableError = 3;
  bool noSecondMotif = false;
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


  setUp.setOption(outputStub, "--outputStub", "output stub for writing out the kept reads");
  setUp.setOption(outOpts.outFilename_, "--filteredCountsOutputFnp", "filtered Counts Output file");
  setUp.setOption(motifStr, "--motif", "motif to search for");
  setUp.setOption(motifAllowableError, "--motifAllowableError", "motif allowable Error");

  setUp.setOption(noSecondMotif, "--noSecondMotif", "no Second Motif");
  if (noSecondMotif) {
    secondaryMotifStr = "";
  }
  setUp.setOption(secondaryMotifStr, "--secondaryMotif", "secondary Motif");
  setUp.setOption(secondaryMotifAllowableError, "--secondaryMotifAllowableError", "secondary motif allowable error");
  setUp.setOption(writeOutFiltered, "--writeOutFiltered", "write Out Filtered");
  setUp.setOption(needBothPairs, "--needBothPairs", "need both pairs to have pattern, by default requires that just one has the pattern");
  setUp.setOption(numOfThreads, "--numOfThreads", "number of threads");

  setUp.finishSetUp(std::cout);

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
  motif filter_motif(motifStr);
  if (motifAllowableError >= filter_motif.size()) {
    std::stringstream ss;
    ss << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << ", error " << "motifAllowableError: " << motifAllowableError << " can't be equal or greater than filter_motif.size(): " << filter_motif.size() << "\n";
    throw std::runtime_error{ss.str()};
  }
  uint32_t motif_passing_score = filter_motif.size() - motifAllowableError;
  std::unique_ptr<motif> secondary_filter_motif;
  if (!secondaryMotifStr.empty()) {
    secondary_filter_motif = std::make_unique<motif>(secondaryMotifStr);
    if (secondaryMotifAllowableError >= secondary_filter_motif->size()) {
      std::stringstream ss;
      ss << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << ", error " << "secondaryMotifAllowableError: " << secondaryMotifAllowableError << " can't be equal or greater than secondary_filter_motif->size(): " << secondary_filter_motif->size() << "\n";
      throw std::runtime_error{ss.str()};
    }
  }
  uint32_t secondary_motif_passing_score = 0;
  if (!secondaryMotifStr.empty()) {
    secondary_motif_passing_score = secondary_filter_motif->size() - secondaryMotifAllowableError;
  }
  out << "inputFnp\ttotalReads\ttotalFiltered\ttotalSecondaryFiltered\n";

  if (!fastq1Fnp.empty() && bfs::exists(fastq1Fnp)) {

    std::function<void()> filteredSeqs;
    uint32_t totalReads = 0;
    uint32_t totalFiltered = 0;
    uint32_t totalSecondaryFiltered = 0;
    std::mutex countMut;
    SeqInput reader(pairedIn);
    reader.openIn();
    if (secondaryMotifStr.empty()) {
      filteredSeqs = [&reader,&multi_seq_io,
            &filter_motif, &motif_passing_score,
            &writeOutFiltered,
            &needBothPairs,
            &totalReads, &totalFiltered, &countMut]() {
            uint32_t current_totalReads = 0;
            uint32_t current_totalFiltered = 0;
            PairedRead pseq;
            while (reader.readNextReadLock(pseq)) {
              ++current_totalReads;
              bool firstMateHasMotif = pseq.seqBase_.seq_.size() > filter_motif.size() &&
                                       filter_motif.passMotifParameter(
                                         pseq.seqBase_.seq_.begin() + (pseq.seqBase_.seq_.size() - filter_motif.size()),
                                         pseq.seqBase_.seq_.end(), motif_passing_score);
              bool secondMateHasMotif = pseq.mateSeqBase_.seq_.size() > filter_motif.size() &&
                                        filter_motif.passMotifParameter(
                                          pseq.mateSeqBase_.seq_.begin() + (
                                            pseq.mateSeqBase_.seq_.size() - filter_motif.size()),
                                          pseq.mateSeqBase_.seq_.end(), motif_passing_score);
              if ((needBothPairs && firstMateHasMotif && secondMateHasMotif) || (
                    !needBothPairs && (firstMateHasMotif || secondMateHasMotif))) {
                if (writeOutFiltered) {
                  multi_seq_io.openWrite("pairs-filtered", pseq);
                }
                ++current_totalFiltered;
              } else {
                multi_seq_io.openWrite("pairs",pseq);
              }
            }
            {
              std::lock_guard<std::mutex> lock(countMut);
              totalFiltered += current_totalFiltered;
              totalReads += current_totalReads;
            }
          };

    } else {
      filteredSeqs = [&reader, &multi_seq_io,
            &filter_motif, &motif_passing_score,
            &secondary_filter_motif, &secondary_motif_passing_score,
            &writeOutFiltered,
            &needBothPairs,
            &totalReads, &totalFiltered, &totalSecondaryFiltered,&countMut]() {
            uint32_t current_totalReads = 0;
            uint32_t current_totalFiltered = 0;
            uint32_t current_totalSecondaryFiltered = 0;
            PairedRead pseq;
            while (reader.readNextReadLock(pseq)) {
              ++current_totalReads;

              bool firstMateHasMotif = pseq.seqBase_.seq_.size() > filter_motif.size() &&
                                       filter_motif.passMotifParameter(
                                         pseq.seqBase_.seq_.begin() + (pseq.seqBase_.seq_.size() - filter_motif.size()),
                                         pseq.seqBase_.seq_.end(), motif_passing_score);
              bool secondMateHasMotif = pseq.mateSeqBase_.seq_.size() > filter_motif.size() &&
                                        filter_motif.passMotifParameter(
                                          pseq.mateSeqBase_.seq_.begin() + (
                                            pseq.mateSeqBase_.seq_.size() - filter_motif.size()),
                                          pseq.mateSeqBase_.seq_.end(), motif_passing_score);
              if ((needBothPairs && firstMateHasMotif && secondMateHasMotif) || (
                    !needBothPairs && (firstMateHasMotif || secondMateHasMotif))) {
                if (writeOutFiltered) {
                  multi_seq_io.openWrite("pairs-filtered", pseq);
                }
                ++current_totalFiltered;
              } else {
                bool firstMateHasSecondaryMotif = pseq.seqBase_.seq_.size() > secondary_filter_motif->size() &&
                                                  secondary_filter_motif->passMotifParameter(
                                                    pseq.seqBase_.seq_.begin() + (
                                                      pseq.seqBase_.seq_.size() - secondary_filter_motif->size()),
                                                    pseq.seqBase_.seq_.end(), secondary_motif_passing_score);
                bool secondMateHasSecondaryMotif = pseq.mateSeqBase_.seq_.size() > secondary_filter_motif->size() &&
                                                   secondary_filter_motif->passMotifParameter(
                                                     pseq.mateSeqBase_.seq_.begin() + (
                                                       pseq.mateSeqBase_.seq_.size() - secondary_filter_motif->size()),
                                                     pseq.mateSeqBase_.seq_.end(), secondary_motif_passing_score);
                if ((needBothPairs && firstMateHasSecondaryMotif && secondMateHasSecondaryMotif) || (
                      !needBothPairs && (firstMateHasSecondaryMotif || secondMateHasSecondaryMotif))) {
                  if (writeOutFiltered) {
                    multi_seq_io.openWrite("pairs-filtered", pseq);
                  }
                  ++current_totalFiltered;
                  ++current_totalSecondaryFiltered;
                } else {
                  multi_seq_io.openWrite("pairs",pseq);
                }
              }
            }
            {
              std::lock_guard<std::mutex> lock(countMut);
              totalFiltered += current_totalFiltered;
              totalSecondaryFiltered += current_totalSecondaryFiltered;
              totalReads += current_totalReads;
            }
          };
    }
    njh::concurrent::runVoidFunctionThreaded(filteredSeqs, numOfThreads);
    out << fastq1Fnp << "\t" << totalReads << "\t" << totalFiltered << "\t" << totalSecondaryFiltered << std::endl;
  }
  if (!fastqFnp.empty() && bfs::exists(fastqFnp)) {

    std::function<void()> filteredSeqs;
    uint32_t totalReads = 0;
    uint32_t totalFiltered = 0;
    uint32_t totalSecondaryFiltered = 0;
    std::mutex countMut;
    SeqInput reader(singleIn);
    reader.openIn();
    if (secondaryMotifStr.empty()) {
      filteredSeqs = [&reader,&multi_seq_io,
            &filter_motif, &motif_passing_score,
            &writeOutFiltered,
            &totalReads, &totalFiltered,&countMut]() {
            uint32_t current_totalReads = 0;
            uint32_t current_totalFiltered = 0;
            seqInfo seq;
            while (reader.readNextReadLock(seq)) {
              ++current_totalReads;
              if (seq.seq_.size() > filter_motif.size() &&
                                       filter_motif.passMotifParameter(
                                         seq.seq_.begin() + (seq.seq_.size() - filter_motif.size()),
                                         seq.seq_.end(), motif_passing_score)) {
                if (writeOutFiltered) {
                  multi_seq_io.openWrite("single-filtered", seq);
                }
                ++current_totalFiltered;
              } else {
                multi_seq_io.openWrite("single", seq);
              }
            }
            {
              std::lock_guard<std::mutex> lock(countMut);
              totalFiltered += current_totalFiltered;
              totalReads += current_totalReads;
            }
          };

    } else {
      filteredSeqs = [&reader, &multi_seq_io,
            &filter_motif, &motif_passing_score,
            &secondary_filter_motif, &secondary_motif_passing_score,
            &writeOutFiltered,
            &totalReads, &totalFiltered, &totalSecondaryFiltered,&countMut]() {
            uint32_t current_totalReads = 0;
            uint32_t current_totalFiltered = 0;
            uint32_t current_totalSecondaryFiltered = 0;

            seqInfo seq;
            while (reader.readNextReadLock(seq)) {
              ++current_totalReads;
              if (seq.seq_.size() > filter_motif.size() &&
                                       filter_motif.passMotifParameter(
                                         seq.seq_.begin() + (seq.seq_.size() - filter_motif.size()),
                                         seq.seq_.end(), motif_passing_score)) {
                if (writeOutFiltered) {
                  multi_seq_io.openWrite("single-filtered", seq);
                }
                ++current_totalFiltered;
              } else {
                if (seq.seq_.size() > secondary_filter_motif->size() &&
                                                  secondary_filter_motif->passMotifParameter(
                                                    seq.seq_.begin() + (
                                                      seq.seq_.size() - secondary_filter_motif->size()),
                                                    seq.seq_.end(), secondary_motif_passing_score)) {
                  if (writeOutFiltered) {
                    multi_seq_io.openWrite("single-filtered", seq);
                  }
                  ++current_totalFiltered;
                  ++current_totalSecondaryFiltered;
                } else {
                  multi_seq_io.openWrite("single", seq);
                }
              }
            }
            {
              std::lock_guard<std::mutex> lock(countMut);
              totalFiltered += current_totalFiltered;
              totalSecondaryFiltered += current_totalSecondaryFiltered;
              totalReads += current_totalReads;
            }
          };
    }
    njh::concurrent::runVoidFunctionThreaded(filteredSeqs, numOfThreads);
    out << fastqFnp << "\t" << totalReads << "\t" << totalFiltered << "\t" << totalSecondaryFiltered<< std::endl;
  }


  return 0;
}


} //namespace njhseq

