#include "ampliconAnalysisSetUp.hpp"
#include <njhcpp/bashUtils.h>

// elucidator - A library for analyzing sequence data
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
//
// This file is part of elucidator.
//
// elucidator is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// elucidator is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with elucidator.  If not, see <http://www.gnu.org/licenses/>.
//
namespace njhseq {

void ampliconAnalysisSetUp::setUpCollapseTandems(
    double& freqCutoff, bool& extra, bool& additionalOut,
    std::string& additionalOutLocationFile) {
  if (needsHelp()) {
    std::cout << "collapseTandems version 1" << std::endl;
    std::cout << "Commands, order not necessary" << std::endl;
    std::cout << "Required commands" << std::endl;
    std::cout
        << "1a) -stub [option]: Stub name for the .fasta and .fasta.qual files"
        << std::endl;
    std::cout << "1b) -fasta [option]: Full name for the fasta file"
              << std::endl;
    std::cout << "1c) -qual [option]: Full name for the qual file" << std::endl;
    std::cout << "1d) -fastq [option]: full name for the fastq file"
              << std::endl;
    std::cout << "Optional commands" << std::endl;
    std::cout << "1) -out [option]: Name of an output directory, will default "
                 "to the stub name plus the date" << std::endl;
    std::cout << "2) -qualThres [option]: Quality threshold for high qual "
                 "mismatch, given in the format of 20,15, where 20 is the "
                 "primary quality and 15 is the second quality" << std::endl;
    std::cout << "3) -gap [option]: Gap penalty, given in the format 7,0.5 "
                 "where 5 is the gap open penalty and 0.5 is the gap extension"
              << std::endl;
    std::cout << "4) -kLength [options]: The length of the k mer check"
              << std::endl;
    std::cout << "5) -kAnywhere : Check k mers anywhere" << std::endl;
    std::cout << "6) -freqCutoff [options] : The frequency at which to "
                 "collapse the sequences into" << std::endl;
    std::cout << "example, collapseTandems -stub MID -out outPut, "
                 "collapseTandems -fasta MID.fasta -v" << std::endl;
    exit(1);
  }
  pars_.ioOptions_.out_.outFilename_ = "tandemsCollapsed";
  pars_.ioOptions_.processed_ = true;
  processDefaultReader(true);
  bool mustMakeDirectory = false;
  processDirectoryOutputName(mustMakeDirectory);
  processVerbose();
  setOption(freqCutoff, "-freqcutoff", "Frequency_multipler_cutoff");
  processRefFilename();
  processKmerProfilingOptions();
  // get the qualities
  processQualThres();
  std::cout << "p: " << pars_.qScorePars_.primaryQual_ << std::endl;
  std::cout << "s: " << pars_.qScorePars_.secondaryQual_ << std::endl;
  // get the gap penatliy
  processGap();
  std::cout << "go: " << pars_.gapInfo_.gapOpen_ << std::endl;
  std::cout << "ge: " << pars_.gapInfo_.gapExtend_ << std::endl;
  setOption(extra, "-extra", "Extra");
  setOption(additionalOutLocationFile, "-additionalOut",
            "AdditionalOutFilename");
  setOption(additionalOut, "-additionalOut", "AdditionalOutputing");
  finishSetUp(std::cout);
}

void ampliconAnalysisSetUp::setUpMarkChimeras() {
  if (needsHelp()) {
    std::cout << "markChimeras version 1" << std::endl;
    std::cout << "Commands, order not necessary" << std::endl;
    std::cout << "Required commands" << std::endl;
    std::cout << "Input, just requires either -stub, -fasta (optional plus "
                 "-qual), or -fastq, this will also determine the output format"
              << std::endl;
    std::cout
        << "1a) -stub [option]: Stub name for the .fasta and .fasta.qual files"
        << std::endl;
    std::cout << "1b) -fasta [option]: Full name of the fasta file"
              << std::endl;
    std::cout << "1c) -qual [option]: Full name of the quality file"
              << std::endl;
    std::cout << "1d) -fastq [option]: Full name of the fastq file"
              << std::endl;
    std::cout << "Optional commands" << std::endl;
    std::cout << "1) -parFreq [option]: Frequency multiplier to be considered "
                 "a parent of the read, defaults to 2.0" << std::endl;
    // std::cout<<"2) -scoreMatrix [option]: Name of an external scoring matrix
    // file for alignment, will default to simple scoring"<<std::endl;
    std::cout << "2) -dout [option]: Name of an output directory, will default "
                 "to outing in the current directory" << std::endl;
    std::cout << "3) -out [option]: Name of an output filename, will default "
                 "to markedChimeras" << std::endl;
    // std::cout<<"4) -qualThres [option]: Quality threshold for high qual
    // mismatch, given in the format of 20,15, where 20 is the primary quality
    // and 15 is the second quality, 20 and 15 are the default"<<std::endl;
    std::cout << "4) -gap [option]: Gap penalty, given in the format 7,0.5 "
                 "where 7 is the gap open penalty and 0.5 is the gap "
                 "extension, 7 and 0.5 are default" << std::endl;
    std::cout << "5) -shorah : add if the names of the reads are in the format "
                 "of NAME_fractionNumber, example FirstRead_0.59" << std::endl;
    // std::cout<<"6) -kAnywhere : Check k mers anywhere, defaults to checking k
    // mers at position found"<<std::endl;
    std::cout << "Examples, sequenceTools markChimeras -fasta Reads.fasta, "
                 "sequenceTools markChimeras -fasta Reads.fasta -parFreq 3"
              << std::endl;
    exit(1);
  }
  pars_.ioOptions_.out_.outFilename_ = "markedChimeras";
  pars_.ioOptions_.processed_ = true;
  // input files
  processDefaultReader(true);
  // out
  processDirectoryOutputName(false);
  // frequency multiplier
  setOption(pars_.chiOpts_.parentFreqs_, "-parfreq", "Parent_frequence_multipler_Cutoff");
  processKmerProfilingOptions();
  processRefFilename();
  // get the qualities
  processQualThres();
  processGap();
  finishSetUp(std::cout);
}

}  // namespace njhseq
