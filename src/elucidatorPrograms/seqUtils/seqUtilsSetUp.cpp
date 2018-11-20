#include "seqUtilsSetUp.hpp"


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





void seqUtilsSetUp::setUpCompareToRef() {


  // input file info
  processReadInNames(true);
  processRefFilename(true);
  processVerbose();
  pars_.gapLeft_ = "0,0";
  pars_.colOpts_.iTOpts_.weighHomopolyer_ = false;
  processAlignerDefualts();
  processDirectoryOutputName(false);
  finishSetUp(std::cout);
  if(pars_.verbose_){
    std::cout << "Gap open: " << pars_.gapInfo_.gapOpen_ << std::endl;
    std::cout << "Gap extend: " << pars_.gapInfo_.gapExtend_ << std::endl;
    std::cout << "Gap Left Query open: " << pars_.gapInfo_.gapLeftQueryOpen_ << std::endl;
    std::cout << "Gap Left Query extend: " << pars_.gapInfo_.gapLeftQueryExtend_ << std::endl;
    std::cout << "Gap Right Query open: " << pars_.gapInfo_.gapRightQueryOpen_ << std::endl;
    std::cout << "Gap Right Query extend: " << pars_.gapInfo_.gapRightQueryExtend_ << std::endl;
    std::cout << "Gap Left Ref open: " << pars_.gapInfo_.gapLeftRefOpen_ << std::endl;
    std::cout << "Gap Left Ref extend: " << pars_.gapInfo_.gapLeftRefExtend_ << std::endl;
    std::cout << "Gap Right Ref open: " << pars_.gapInfo_.gapRightRefOpen_ << std::endl;
    std::cout << "Gap Right Ref extend: " << pars_.gapInfo_.gapRightRefExtend_ << std::endl;
  }
}


void seqUtilsSetUp::setUpMapToReferenceCount(bool& extra) {
	pars_.ioOptions_.lowerCaseBases_ = "upper";
  if (needsHelp()) {
    std::stringstream tempOut;
    tempOut << "mapCount version 1" << std::endl;
    tempOut
        << "Takes a sequence file and maps to several references sequences and "
           " counts the best matches to each reference seq, will split the "
           "count on ties" << std::endl;
    tempOut << "Commands, order not necessary" << std::endl;
    tempOut << njh::bashCT::bold <<
    		njh::bashCT::black <<"Required commands" << njh::bashCT::reset << std::endl;
    //printInputUsage(tempOut);
    //printReferenceComparisonUsage(tempOut);
    tempOut << njh::bashCT::bold <<
    		njh::bashCT::black <<"Optional commands" << njh::bashCT::reset << std::endl;
    tempOut << "1) -dout [option]: Name of an output directory, will default "
               "to InFileName_NameOfProgram_CURRENT_DATE" << std::endl;
    //printAlignmentUsage(tempOut);
    //printAlnInfoDirUsage(tempOut);
    std::cout << cleanOut(tempOut.str(), width_, indent_);
    tempOut.str(std::string());
    std::cout << "example, mapCount -fasta seqs.fasta -ref refSeqs.fasta"
              << std::endl;
    tempOut << njh::bashCT::bold <<
    		njh::bashCT::black <<"Output Files: " << njh::bashCT::reset << std::endl;
    tempOut << "1) mapFreqInfo.tab.txt : File containg the frequncy information"
            << std::endl;
    tempOut << "2) mismatchInfo.tab.txt : File containing infomation about the "
               "mistmatches to reference the sequences had" << std::endl;
    tempOut << "3) clusters : Directory with all the reads split into fasta "
               "files by their best reference match" << std::endl;
    tempOut << "4) runLog.txt : File with information about the running of the "
               "program" << std::endl;
    std::cout << cleanOut(tempOut.str(), width_, indent_);
    tempOut.str(std::string());
    exit(1);
  }

  processDefaultReader(true);
  processRefFilename(true);
  processDirectoryOutputName(true);
  processVerbose();
  // get the gap penatliy
  pars_.gapLeft_ = "0,0";
  pars_.gap_ = "5,1";
  processAlignerDefualts();
  if (pars_.verbose_){
    std::cout << "Gap open: " << pars_.gapInfo_.gapOpen_ << std::endl;
    std::cout << "Gap extend: " << pars_.gapInfo_.gapExtend_ << std::endl;
    std::cout << "Gap Left Query open: " << pars_.gapInfo_.gapLeftQueryOpen_ << std::endl;
    std::cout << "Gap Left Query extend: " << pars_.gapInfo_.gapLeftQueryExtend_ << std::endl;
    std::cout << "Gap Right Query open: " << pars_.gapInfo_.gapRightQueryOpen_ << std::endl;
    std::cout << "Gap Right Query extend: " << pars_.gapInfo_.gapRightQueryExtend_ << std::endl;
    std::cout << "Gap Left Ref open: " << pars_.gapInfo_.gapLeftRefOpen_ << std::endl;
    std::cout << "Gap Left Ref extend: " << pars_.gapInfo_.gapLeftRefExtend_ << std::endl;
    std::cout << "Gap Right Ref open: " << pars_.gapInfo_.gapRightRefOpen_ << std::endl;
    std::cout << "Gap Right Ref extend: " << pars_.gapInfo_.gapRightRefExtend_ << std::endl;
  }

  setOption(extra, "-extra", "Extra");
  finishSetUp(std::cout);
}







}  // namespace njh
