
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




void seqUtilsModSetUp::setUpSplit(std::string& splitOption, uint32_t& minLen,
		uint32_t& within, std::string& runCutoffString,
                               std::string& nameContains,
                               std::string& seqContains, uint32_t& occurences,
															 uint32_t& maxLength, uint32_t& qualWindowSize,
															 uint32_t& qualWindowStep, uint32_t& qualWindowThres) {
  if (needsHelp()) {
    std::cout << "split" << std::endl;
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
    std::cout << "SplitOptions" << std::endl;
    std::cout << "-splitOption [option]: The spliting option" << std::endl;
    std::cout << "onReadLength, below, withinGiven, withinMean, nameContains, "
                 "seqContains, clusterSize" << std::endl;
    std::cout << "split -splitOption onReadLength" << std::endl;
    std::cout << "\tThis will exclude any reads that have a read length "
                 "greater than 2 standard deviations away from the mean read "
                 "length" << std::endl;
    std::cout << "split -splitOption below -minLen [option]" << std::endl;
    std::cout << "\tThis will exclude any reads below the given -minLen"
              << std::endl;
    std::cout
        << "split -splitOption withinGiven -given [option] -within [option]"
        << std::endl;
    std::cout << "\tExclude reads with length outside the -given length option "
                 "plus/minus the -within option" << std::endl;
    std::cout << "split -splitOption withMean -withIn [option]" << std::endl;
    std::cout
        << "\tExclude reads outside the mean read length by the -within option"
        << std::endl;
    std::cout << "split -splitOption nameContains -contains [option] "
              << std::endl;
    std::cout << "\tExclude reads that names contains the -contains option"
              << std::endl;
    std::cout << "split -splitOption seqContains -contains [option]"
              << std::endl;
    std::cout
        << "\tExclude reads that sequence contains the given -contains option"
        << std::endl;
    std::cout << "split -splitOption clusterSize -sizeCutoff [option]"
              << std::endl;
    std::cout << "\tExclude reads that have a count equal to or less than the "
                 "given -sizeCutoff option, size is donated by the name "
                 "extension _XX, where XX is a frequency" << std::endl;
    std::cout << "example, split -stub MID2 -minLen 200" << std::endl;
    exit(1);
  }

  // input file info
  processDefaultReader(true);

  if (setOption(splitOption, "-splitOption", "SplitOption", true)) {
    splitOption = stringToLowerReturn(splitOption);
  }
  if (splitOption == "below") {
    setOption(minLen, "-minLen", "MinimumLength");
  } else if (splitOption == "onreadlength") {

  } else if (splitOption == "withingiven") {
    setOption(within, "-within", "Within");
    setOption(minLen, "-given", "MinimumLength");
  } else if (splitOption == "withinmean") {
    setOption(within, "-within", "Within");
  } else if (splitOption == "namecontains") {
    setOption(nameContains, "-namecontains,-contains", "NameContains", true);
  } else if (splitOption == "clustersize") {
    setOption(runCutoffString, "-runCutoff,-sizeCutoff", "ClusterSizeCutoff");
  } else if (splitOption == "seqcontains") {
    setOption(occurences, "-occurences", "Occurences");
    setOption(seqContains, "-seqContains,-contains", "SeqContains", true);
  } else if (splitOption == "above") {
    setOption(maxLength, "-maxLen,-maxLength", "maximumLength");
  } else if (splitOption == "between") {
    setOption(maxLength, "-maxLen,-maxLength", "maximumLength");
    setOption(minLen, "-given", "MinimumLength");
  } else if (splitOption == "nucleotidecomp" ||
             splitOption == "nucleotidecomposition") {
    splitOption = "nucleotidecomp";
  } else if (splitOption == "qualitywindow") {
    std::string qualWindowString = "";
    setOption(qualWindowString, "-qualWindow", "QualityWindowParams");
    if (setOption(qualWindowString, "-qualWindow", "Sliding Quality Window, goes WindowSize,WindowStep,Thresdhold")) {
    	uint32_t qualityWindowLengthLocal = 50;
    	uint32_t qualityWindowStepLocal = 5;
    	uint32_t qualityWindowThresLocal = 25;
      seqUtil::processQualityWindowString(qualWindowString, qualityWindowLengthLocal,
      		qualityWindowStepLocal, qualityWindowThresLocal);
      qualWindowSize = qualityWindowLengthLocal;
      qualWindowStep = qualityWindowStepLocal;
      qualWindowThres = qualityWindowThresLocal;
    } else {
    	qualWindowSize = 50;
      qualWindowStep = 5;
      qualWindowThres = 25;
    }
  }
  // get the out directory
  processDirectoryOutputName(false);

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
