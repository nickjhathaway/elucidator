
//  seqUtilsConverterSetUp.cpp
//
//  Created by Nicholas Hathaway on 2015/05/28.
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
#include "seqUtilsConverterSetUp.hpp"
    
    
namespace njhseq {


void seqUtilsConverterSetUp::setUpConvertFiles() {
  if (needsHelp()) {

    std::cout << "convertFiles" << std::endl;
    std::cout << "Commands, order not necessary" << std::endl;
    std::cout << "Required commands" << std::endl;
    std::cout
        << "1a) -stub [option]: Stub name for the .fasta and .fasta.qual files"
        << std::endl;
    std::cout << "1b) -fasta [option]: Full name of the fasta file"
              << std::endl;
    std::cout << "1c) -qual [option]: Full name of the quality file"
              << std::endl;
    std::cout << "1d) -fastq [option]: Full name of the fastq file"
              << std::endl;
    std::cout << "1e) -sff [option]: Full name of a sff.txt file (converted "
                 "into a text file), output will be fastq" << std::endl;
    std::cout << "Optional commands" << std::endl;
    std::cout << "examples, convertFiles -stub MID, convertFiles -fastq "
                 "MID.fastq, convertFiles -fasta MID.fasta -qual MID.fasta.qual"
              << std::endl;
    exit(1);
  }
  // input file info
  processDefaultReader(true);
  if (pars_.ioOptions_.out_.outFilename_ == "out") {
  	pars_.ioOptions_.out_.outFilename_ = njh::files::removeExtension(pars_.ioOptions_.firstName_);
  }
  if (pars_.ioOptions_.inFormat_ == SeqIOOptions::inFormats::FASTQ) {
  	pars_.ioOptions_.outFormat_ = SeqIOOptions::outFormats::FASTAQUAL;
  	pars_.ioOptions_.out_.outExtention_ = SeqIOOptions::getOutExtension(pars_.ioOptions_.outFormat_);
  } else if (pars_.ioOptions_.inFormat_ == SeqIOOptions::inFormats::SFFBIN ||
  		pars_.ioOptions_.inFormat_ == SeqIOOptions::inFormats::SFFTXT) {
  	pars_.ioOptions_.outFormat_ = SeqIOOptions::outFormats::FASTQ;
  	pars_.ioOptions_.out_.outExtention_ = SeqIOOptions::getOutExtension(pars_.ioOptions_.outFormat_);
  } else if (pars_.ioOptions_.inFormat_ == SeqIOOptions::inFormats::FASTAQUAL) {
  	pars_.ioOptions_.outFormat_ = SeqIOOptions::outFormats::FASTQ;
  	pars_.ioOptions_.out_.outExtention_ = SeqIOOptions::getOutExtension(pars_.ioOptions_.outFormat_);
  } else if (pars_.ioOptions_.inFormat_ == SeqIOOptions::inFormats::BAM) {
  	pars_.ioOptions_.outFormat_ = SeqIOOptions::outFormats::FASTQ;
  	pars_.ioOptions_.out_.outExtention_ = SeqIOOptions::getOutExtension(pars_.ioOptions_.outFormat_);
  } else if (pars_.ioOptions_.inFormat_ == SeqIOOptions::inFormats::FASTA) {
  	pars_.ioOptions_.outFormat_ = SeqIOOptions::outFormats::FASTQ;
  	pars_.ioOptions_.out_.outExtention_ = SeqIOOptions::getOutExtension(pars_.ioOptions_.outFormat_);
  }

  /**@todo make it so you can set outFormat
  if(commands_.hasFlagCaseInsen("-outFormat")){
  	pars_.ioOptions_.outFormat_ = SeqIOOptions::getOutFormat()
  }*/

  finishSetUp(std::cout);
}
} // namespace njhseq
