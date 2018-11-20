
//  seqUtilsExtractSetUp.cpp
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
    
#include "seqUtilsExtractSetUp.hpp"
    
    
namespace njhseq {




//
//void seqUtilsExtractSetUp::setUpExtractByMID(ExtractByMIDPars & pars) {
//  if (needsHelp()) {
//    std::cout << "extractByMID" << std::endl;
//    std::cout << "Commands, order not necessary" << std::endl;
//    std::cout << "Required commands" << std::endl;
//    std::cout
//        << "1a) -stub [option]: Stub name for the .fasta and .fasta.qual files"
//        << std::endl;
//    std::cout << "1b) -fasta [option]: Full name of the fasta file"
//              << std::endl;
//    std::cout << "1c) -qual [option]: Full name of the quality file"
//              << std::endl;
//    std::cout << "1d) -fastq [option]: Full name of the fastq file"
//              << std::endl;
//    std::cout << "2) -mids [option]: Name of the mids file" << std::endl;
//    std::cout << "Optional commands" << std::endl;
//    std::cout << "1) -midFileDelim [option]: Delimiter of mids file, defaults "
//                 "to whitespace delimited" << std::endl;
//    std::cout << "2) -checkComplement: add if the complement sequences should "
//                 "also be checked" << std::endl;
//    std::cout << "example, programName -stub MID2 -mids midFilename.tab.txt"
//              << std::endl;
//    exit(1);
//  }
//  processVerbose();
//  processDebug();
//  pars_.ioOptions_.lowerCaseBases_ = "remove";
//  processDefaultReader(true);
//  setOption(pars.idFilename, "-mids", "MID_Filename", true);
//  setOption(pars.idFileDelim, "-midFileDelim", "File Delimiter For ID File");
//  setOption(pars.mDetPars.checkComplement_, "-checkComplement", "Check reverse complement of the reads as well");
//  processDirectoryOutputName(true);
//  setOption(pars.mDetPars.variableStop_, "-variableStart", "Up to where to search for barcodes in the sequence");
//  setOption(pars.mDetPars.barcodesBothEnds_, "-barcodesBothEnds", "If the reads have barcodes on both ends");
//	setOption(pars.mDetPars.checkForShorten_, "-checkShortenBars",
//			"Check for shorten Barcodes if the first base may have been trimmed off");
//
//  setOption(pars.barcodeErrors, "-barcodeErrors", "How many errors to allow in the barcode");
//  setOption(pars.smallFragmentCutoff, "-smallFragmentCutoff", "Size of the sequences to be considered small and not searched for barcodes");
//  finishSetUp(std::cout);
//}



void seqUtilsExtractSetUp::setUpExtractSameSeqs(std::string& compareFileName) {
  if (needsHelp()) {
    std::cout << "extraSameSeqs" << std::endl;
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
    std::cout << "2) -compareFileName [option]: name of the file to compare to"
              << std::endl;
    exit(1);
  }
  // input file info
  processDefaultReader(true);
  setOption(compareFileName, "-compare", "ComprisonFile", true);
  if (pars_.ioOptions_.out_.outFilename_ == "out") {
  	pars_.ioOptions_.out_.outFilename_ =
        "extracted_" + njh::files::removeExtension(pars_.ioOptions_.firstName_);
  }
  finishSetUp(std::cout);
}

} // namespace njhseq
