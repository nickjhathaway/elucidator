
//  seqUtilsInfoSetUp.cpp
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
    
#include "seqUtilsInfoSetUp.hpp"
    
    
namespace njhseq {


void seqUtilsInfoSetUp::setUpFastaIdenticalInfo(std::string& qualRep,
                                            int& qualCheck) {
  if (needsHelp()) {
    std::cout << commands_.getProgramName() + " version 1" << std::endl;
    std::cout << "Commands, order not necessary" << std::endl;
    std::cout << "Required command" << std::endl;
    std::cout
        << "1a) -stub [option]: Stub name for the .fasta and .fasta.qual files"
        << std::endl;
    std::cout << "1b) -fasta [option]: Full name for the fasta file"
              << std::endl;
    std::cout << "1c) -qual [option]: Full name for the qual file" << std::endl;
    std::cout << "1d) -fastq [option]: full name for the fastq file"
              << std::endl;
    std::cout << "Optional command" << std::endl;
    std::cout << "1) -out [option]: Name of an output directory, will default "
                 "to the stub name plus the date" << std::endl;
    std::cout << "2) -lower [option]: Give remove, upper, or nothing to say "
                 "what to do about lower case bases" << std::endl;
    std::cout << "3) -qualRep [option]: Given the options of worst, "
                 "median(default), bestSeq, bestQual, or average to set the "
                 "setting of the identical qual rep if using quality"
              << std::endl;
    std::cout << "4) -qualCheck [option]: Give a single integer to count bases "
                 "above that quality " << std::endl;
    std::cout << "example, " + commands_.getProgramName() + " -fastq MID.fastq -out MID_info" << std::endl;
    exit(1);
  }

  // input file info
  processDefaultReader(true);
  processDirectoryOutputName(true);
  processRefFilename();
  setOption(qualRep, "-qualRep", "QualityRepresentativeForClusters");
  setOption(qualCheck, "-qualCheck", "QualityCheckCutoff");
  // check for unused commands
  finishSetUp(std::cout);
}

} // namespace njhseq
