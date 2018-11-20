
//  seqUtilsTrimSetUp.cpp
//
//
//  Created by Nicholas Hathaway on 2016/10/05.
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
    
#include "seqUtilsTrimSetUp.hpp"
    
    
namespace njhseq {


void seqUtilsTrimSetUp::setUpTrimFront(FullTrimReadsPars & pars) {

  processVerbose();
  processDebug();
	setOption(pars.numberOfFowardBases, "--forwardBases", "Number of bases to trim off front", true);;
	processIoOptions();
	finishSetUp(std::cout);
}

void seqUtilsTrimSetUp::setUpTrimAtFirstQual(FullTrimReadsPars & pars) {

  processVerbose();
  processDebug();
	setOption(pars.qual, "--qScore", "Quality to trim at", true);;
	processIoOptions();
	finishSetUp(std::cout);
}




void seqUtilsTrimSetUp::setUpTrimAtLastBase(FullTrimReadsPars & pars) {

  processVerbose();
  processDebug();
  std::string baseStr = std::string(1, pars.base);
	setOption(baseStr, "--base", "Base to trim at", true);
	pars.base = baseStr.front();
	processIoOptions();
	finishSetUp(std::cout);
}


void seqUtilsTrimSetUp::setUpTrimAtFirstBase(FullTrimReadsPars & pars) {

  processVerbose();
  processDebug();
  std::string baseStr = std::string(1, pars.base);
	setOption(baseStr, "--base", "Base to trim at", true);
	pars.base = baseStr.front();
	processIoOptions();
	finishSetUp(std::cout);
}

void seqUtilsTrimSetUp::setUpTrimEnd(FullTrimReadsPars & pars) {

  processVerbose();
  processDebug();
	setOption(pars.numberOfEndBases, "--endBases", "Number of bases to trim off end", true);
	processIoOptions();
	finishSetUp(std::cout);
}

void seqUtilsTrimSetUp::setUpTrimEnds(FullTrimReadsPars & pars) {

  processVerbose();
  processDebug();
	setOption(pars.numberOfEndBases, "--endBases", "Number of bases to trim off end", true);
	setOption(pars.numberOfFowardBases, "--forwardBases", "Number of bases to trim off front", true);;
	processIoOptions();
	finishSetUp(std::cout);
}

void seqUtilsTrimSetUp::setUpTrimBeforeSeq(FullTrimReadsPars & pars) {

  processVerbose();
  processDebug();
	setOption(pars.forwardSeq, "--forwardSeq", "Trim to this sequence", true);

	processIoOptions();
	procesingTrimmingWithSeqsOpts(pars);
	finishSetUp(std::cout);
}

void seqUtilsTrimSetUp::setUpTrimFromSeq(FullTrimReadsPars & pars) {

  processVerbose();
  processDebug();

  setOption(pars.backSeq, "--backSeq", "Trim from this sequence on", true);

	processIoOptions();
	procesingTrimmingWithSeqsOpts(pars);
	finishSetUp(std::cout);
}


void seqUtilsTrimSetUp::setUpTrimBetweenSeqs(FullTrimReadsPars & pars) {

  processVerbose();
  processDebug();

  setOption(pars.backSeq, "--backSeq", "Trim from this sequence on", true);
	setOption(pars.forwardSeq, "--forwardSeq", "Trim to this sequence", true);

	processIoOptions();
	procesingTrimmingWithSeqsOpts(pars);
	finishSetUp(std::cout);
}

void seqUtilsTrimSetUp::setUpTrimToLen(FullTrimReadsPars & pars) {

  processVerbose();
  processDebug();
	setOption(pars.maxLength, "--length,-l", "Length to trim to (only reads of this len or greater be written)", true);
	setOption(pars.keepOnlyOn, "--keepOnlyOn", "Keep Only the Reads that are still on");

	processIoOptions();
	finishSetUp(std::cout);
}

void seqUtilsTrimSetUp::setUpTrimToSimilarSeq(FullTrimReadsPars & pars) {

  processVerbose();
  processDebug();
  /**@todo add a way to trim to similar condensed sequence */
  //bool useCondensed = false;
  //setUp.setOption(useCondensed, "--useCondensed", "Whether to use the condensed sequence instead");
  setOption(pars.windowLength, "--windowLen,-w", "Window length to expand the kmer search at which to trim");
  setOption(pars.maxLength,   "--trimLen,-t",   "Length to trim at approximately", true);
  setOption(pars.kmerLength,   "--kmerLen,-k",   "Kmer Length");
  if(pars.windowLength > pars.maxLength){
  	std::stringstream ss;
  	ss << __PRETTY_FUNCTION__ << ": Error, Window length shouldn't be longer than trim length" << std::endl;
  	ss <<"Window length: " << pars.windowLength  << ", trim length: " << pars.maxLength << std::endl;
  	addWarning(ss.str());
  	failed_ = true;
  }

	processIoOptions();
	finishSetUp(std::cout);
}




void seqUtilsTrimSetUp::procesingTrimmingWithSeqsOpts(FullTrimReadsPars & pars){

	bool global = false;
	setOption(global, "--global", "Match with Global alignment");
	pars.tSeqPars_.local_ = !global;
	if(global){
		//if left make it semi-global by defautl;
		pars_.gapRight_ = "0,0";
		pars_.gapLeft_ = "0,0";
	}
  pars_.scoring_ = substituteMatrix::createDegenScoreMatrix(pars_.generalMatch_, pars_.generalMismatch_);
  processAlignerDefualts();
	setOption(pars.tSeqPars_.includeSequence_, "--keepSeq", "Keep Sequence used for trimming");
	setOption(pars.tSeqPars_.sequenceToLowerCase_, "--seqToLower", "Make sequence used for trimming lowercase");
  setOption(pars.tSeqPars_.removePreviousSameBases_, "--removePreviousBases", "Trim any bases match the bases at the end of the sequence used for trimming");
	setOption(pars.keepOnlyOn, "--keepOnlyOn", "Keep Only the Reads that are still on");
	setOption(pars.tSeqPars_.alwaysTrim, "--alwaysTrim", "Always Trim even when the sequence doesn't match well");

	processComparison(pars.allowableErrors);

	setOption(pars.allowableErrors.distances_.query_.coverage_, "--coverage",
			"Amount of query sequence that needs to be found to trim to");
	setOption(pars.tSeqPars_.within_, "--within","Look for sequence only within this many bases of the front or back");
}

void seqUtilsTrimSetUp::processIoOptions(){
	processIoOptions(readInFormatsAvailable_);
}
void seqUtilsTrimSetUp::processIoOptions(const VecStr & formats){
	processDefaultReader(formats, true);
	if (pars_.ioOptions_.out_.outFilename_ == "out") {
		auto inputPath = njh::files::bfs::path(pars_.ioOptions_.firstName_);
		inputPath.filename().replace_extension("");
		pars_.ioOptions_.out_.outFilename_ = njh::files::prependFileBasename(inputPath, "trimmed_");
	}
}



} // namespace njhseq
