
//  seqUtilsInfoSetUp.cpp
//
//  Created by Nicholas Hathaway on 2015/05/28.
//  Copyright (c) 2015 Nicholas Hathaway. All rights reserved.
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
