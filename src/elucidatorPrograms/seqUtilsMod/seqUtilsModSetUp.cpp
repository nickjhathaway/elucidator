
//  seqUtilsModSetUp.cpp
//
//  Created by Nicholas Hathaway on 2015/05/28.
//  Copyright (c) 2015 Nicholas Hathaway. All rights reserved.
//

    
#include "seqUtilsModSetUp.hpp"
    
    
namespace njhseq {

void seqUtilsModSetUp::setUpSortReads(std::string& sortBy, bool& decending) {
  if (needsHelp()) {

    std::cout << "sortReads" << std::endl;
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
    std::cout << "SortOptions: " << std::endl;
    std::cout
        << "-ascending : The default sorting is decending but this can be";
    std::cout << "used to change it to ascending" << std::endl;
    std::cout << "example, sequenceTools sortReads -stub MID2 -sortBy seq"
              << std::endl;
    exit(1);
  }
  // input file info
  processDefaultReader(true);
  if (pars_.ioOptions_.out_.outFilename_ == "out") {
  	pars_.ioOptions_.out_.outFilename_ = "sorted_" + njh::files::removeExtension(pars_.ioOptions_.firstName_);
  }
  setOption(sortBy, "--sortBy", "Sort by Option");
  bool ascending = false;
  setOption(ascending, "--ascending", "Ascending Sort");
  decending = !ascending;
  finishSetUp(std::cout);
}



void seqUtilsModSetUp::setUpRenameIDs(std::string& stub, std::string& sortBy,
                                   bool& keepChimeraFlag) {
  if (needsHelp()) {

    std::cout << "renameIDs" << std::endl;
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
    std::cout << "2) -name [option]: new stub name for the ids" << std::endl;
    std::cout << "Optional commands" << std::endl;
    std::cout << "1) -sortBy [option]: some to sort the reads by first, "
                 "defaults to the order they came in" << std::endl;
    std::cout << "example, renameIDs -stub MID2 -name newName " << std::endl;
    exit(1);
  }

  // input file info
  processDefaultReader(true);
	if (pars_.ioOptions_.out_.outFilename_ == "out") {
		pars_.ioOptions_.out_.outFilename_ = njh::files::prependFileBasename(pars_.ioOptions_.firstName_, "renamed_").string();
	}
  setOption(stub, "-name", "NewStubName");
  setOption(sortBy, "-sortBy", "SortOption");
  setOption(keepChimeraFlag, "-keepChimeraFlag", "KeepChimeraFlag");
  finishSetUp(std::cout);
}




void seqUtilsModSetUp::setUpComplementSeq(std::string &seqType) {
  if (needsHelp()) {
    std::cout << "complementSeq" << std::endl;
    std::cout << "Required commands" << std::endl;
    std::cout << "1) -seq [option]: The sequence to be counted" << std::endl;
    std::cout << "Optional Commands" << std::endl;
    std::cout << "2) -seqType [options]: The type of sequence, choices are DNA "
                 "(default) or RNA" << std::endl;
    exit(1);
  }
  if (!processDefaultReader(false)) {
    	processSeq(true);
  } else {
    if (pars_.ioOptions_.out_.outFilename_ == "out") {
    	pars_.ioOptions_.out_.outFilename_ = "C_" + njh::files::removeExtension(pars_.ioOptions_.firstName_);
    }
  }

  setOption(seqType, "-seqType", "SequenceType");
  finishSetUp(std::cout);
}

void seqUtilsModSetUp::setUpTranslate(uint64_t &start, bool &complement,
                                   bool &reverse) {
  if (needsHelp()) {
    std::cout << "translate" << std::endl;
    std::cout << "Required commands" << std::endl;
    std::cout << "1) -seq [option]: The sequence to be counted" << std::endl;
    std::cout << "Optional Commands" << std::endl;
    std::cout << "2) -start [options]: The start of the reading frame"
              << std::endl;
    std::cout << "2) -complement : Whether to complement sequence first" << std::endl;
    std::cout << "3) -reverse : Whether to reverse the sequence first" << std::endl;
    exit(1);
  }
	processSeq(!processDefaultReader(false));
	if (pars_.ioOptions_.out_.outFilename_ == "out" || "" == pars_.ioOptions_.out_.outFilename_) {
		auto noExt = pars_.ioOptions_.firstName_;
		noExt.replace_extension("");
		pars_.ioOptions_.out_.outFilename_ = njh::files::prependFileBasename(noExt, "translated_");
		pars_.ioOptions_.outFormat_ = SeqIOOptions::outFormats::FASTA;
		pars_.ioOptions_.out_.outExtention_ = ".fasta";
	}
	setOption(start, "-start", "Reading Frame Start");
	setOption(complement, "-complement", "Complement sequence");
	setOption(reverse, "-reverse", "Reverse sequence");
  finishSetUp(std::cout);
}


} // namespace njhseq
