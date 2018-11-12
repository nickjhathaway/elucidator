
//  seqUtilsExtractSetUp.cpp
//
//  Created by Nicholas Hathaway on 2015/05/28.
//  Copyright (c) 2015 Nicholas Hathaway. All rights reserved.
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


void seqUtilsExtractSetUp::setUpGetSimilarSequences(bool& useNucComp, double& maxNucCompDiff, double& idCutOff,
		double& gapCutoff, double& queryCutoff){
	if (needsHelp()) {
		std::stringstream tempOut;
		tempOut << njh::bashCT::bold
				<< njh::centerText("getSimilarSequences", width_)
		<< njh::bashCT::reset << std::endl;
		tempOut << "Commands, order not necessary" << std::endl;
		tempOut << "Required commands" << std::endl;
		tempOut
				<< "1a) -stub [option]: Stub name for the .fasta and .fasta.qual files"
				<< std::endl;
		tempOut << "1b) -fasta [option]: Full name of the fasta file"
							<< std::endl;
		tempOut << "1c) -qual [option]: Full name of the quality file"
							<< std::endl;
		tempOut << "1d) -fastq [option]: Full name of the fastq file"
							<< std::endl;
		tempOut << "2) -seq [option] : can be the actual sequence or in txt file where the first line is the sequence "
				"or fasta file and the first read will be used.  This is the sequence to compare against" << std::endl;
		tempOut << "Optional command options, order not necessary" << std::endl;
		tempOut <<"1) -id [option] : percent identify cut off should be between 0 and 1, this is the percent identity"
				" of the bases that are overlapping, not all the bases in the sequence, defaults to 0.97" << std::endl;
		tempOut <<"2) -gapCutoff [option] : the percentage of gaps in the overlapping bases allowed, should be between 0 and 1,"
				" the higher the number the more gaps allowed,defaults to 0.1" << std::endl;
		tempOut <<"3) -coverage [option] : The amount of sequence's bases that have to be overlapping the target/ref sequence,"
				" should be between 0 and 1, defaults to 0.75" << std::endl;
		tempOut <<"4) -useNucComp : this means it will count sequences as dissimilar if the sum of differences in the nucleotide "
				" compositions of the two sequences is greater than -maxNucCompDiff, saves on time since the sequences don't have to be"
				" aligned, the default is to not use this option" << std::endl;
		tempOut <<"5) -maxNucCompDiff [option] : the max sum of difference in nucleotide compositions allowed if -useNucComp is used,"
				" should be between 0 and 1, default is 0.1" << std::endl;
		//printAlignmentUsage(tempOut);
		tempOut << "examples "  << std::endl;
		tempOut << "sequenceTools getSimilarSequences -fastq seqs.fastq -seq ref.fasta -id 0.9"
							<< std::endl;
		tempOut << "sequenceTools getSimilarSequences -fastq seqs.fastq -seq ref.fasta -id 0.9"
									<< std::endl;
		std::cout << cleanOut(tempOut.str(), width_, indent_);
		exit(1);
	}

  setOption(useNucComp, "-useNucComp", "useNucComp");
  setOption(maxNucCompDiff, "-maxNucCompDiff", "maxNucCompDiff");
  setOption(idCutOff, "-idCutOff,-id", "idCutOff");
  setOption(gapCutoff, "-gapCutoff", "gapCutoff");
  setOption(queryCutoff, "-queryCutoff,-cover,-coverage", "queryCutoff");
  processDefaultReader(true);
  processDirectoryOutputName(true);
  processSeq(true);
  processAlignerDefualts();
  processVerbose();
  finishSetUp(std::cout);

}


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
