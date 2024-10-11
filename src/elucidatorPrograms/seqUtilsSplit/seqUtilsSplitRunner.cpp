
//  seqUtilsSplitRunner.cpp
//
//
//  Created by Nicholas Hathaway on 2018/11/18.
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
    
#include "seqUtilsSplitRunner.hpp"
#include <njhseq.h>
    
namespace njhseq {




seqUtilsSplitRunner::seqUtilsSplitRunner()
    : njh::progutils::ProgramRunner({
  		addFunc("SeqSplitOnLenBelow", SeqSplitOnLenBelow, false),
  		addFunc("SeqSplitOnLenWithin", SeqSplitOnLenWithin, false),
  		addFunc("SeqSplitOnLenWithinMedianLen", SeqSplitOnLenWithinMedianLen, false),


  		addFunc("SeqSplitOnNameContains", SeqSplitOnNameContains, false),
  		addFunc("SeqSplitOnSeqContains", SeqSplitOnSeqContains, false),
  		addFunc("SeqSplitOnLenAbove", SeqSplitOnLenAbove, false),
  		addFunc("SeqSplitOnLenBetween", SeqSplitOnLenBetween, false),
  		addFunc("SeqSplitOnQualityWindow", SeqSplitOnQualityWindow, false),
  		addFunc("SeqSplitOnQualityCheck", SeqSplitOnQualityCheck, false),
  		addFunc("SeqSplitOnNucelotideComp", SeqSplitOnNucelotideComp, false),
			addFunc("getSimilarSequences", getSimilarSequences, false),
			addFunc("getSimilarSequencesByKDist", getSimilarSequencesByKDist, false),
			addFunc("SeqSplitOnCount", SeqSplitOnCount, false),
			addFunc("SeqSplitOnNameContainsPattern", SeqSplitOnNameContainsPattern, false),
			addFunc("SeqSplitOnQualityCheck", SeqSplitOnQualityCheck, false),


},
                    "seqUtilsSplit") {}

//




int seqUtilsSplitRunner::getSimilarSequences(const njh::progutils::CmdArgs & inputCommands) {
  // remove sequences very dissimilar to input compare sequence
  seqSetUp setUp(inputCommands);
  double idCutOff = 0.97;
  double gapCutoff = 0.03;
  double queryCutoff = 0.75;
  bool useNucComp = false;
  bool useScore = false;
  double maxNucCompDiff = 0.1;
  double scoreCutOff = 0;
  bool trimSimilar = false;
  double percentile = 0.95;
  bool checkComplement = false;
	setUp.setOption(percentile, "-percentile", "percentile");
	setUp.setOption(trimSimilar, "-trimSimilar,-trim", "trimSimilar");
  if(setUp.setOption(useScore, "-useAlnScore,-useScore,-useAlignmentScore", "useAlignmentScore")){
  	setUp.setOption(scoreCutOff, "-scoreCutOff,-cutOff", "scoreCutOff");
  }
  setUp.setOption(checkComplement, "-checkComplement,-complement", "checkComplement");
  setUp.setOption(useNucComp, "-useNucComp", "useNucComp");
  setUp.setOption(maxNucCompDiff, "-maxNucCompDiff", "maxNucCompDiff");
  setUp.setOption(idCutOff, "-idCutOff,-id", "idCutOff");
  setUp.setOption(gapCutoff, "-gapCutoff", "gapCutoff");
  setUp.setOption(queryCutoff, "-queryCutoff,-cover,-coverage", "queryCutoff");
  setUp.processDefaultReader(true);
  setUp.processDirectoryOutputName(true);
  setUp.processSeq(true);
  setUp.processAlignerDefualts();
  setUp.processVerbose();
  setUp.startARunLog(setUp.pars_.directoryName_);
  SeqInput reader(setUp.pars_.ioOptions_);
  	reader.openIn();
	auto inReads = reader.readAllReads<readObject>();
  uint64_t maxReadLength = 0;
  readVec::getMaxLength(inReads, maxReadLength);
  if (setUp.pars_.seqObj_.seqBase_.seq_.size() > maxReadLength) {
    maxReadLength = len(setUp.pars_.seqObj_.seqBase_.seq_);
  }

  // aligner object
  auto scoringMatrixMap = setUp.pars_.scoring_;
  KmerMaps emptyMaps;
  gapScoringParameters gapPars(setUp.pars_.gapInfo_);
  aligner alignerObj = aligner(maxReadLength, gapPars,
  		 scoringMatrixMap, false);

  std::vector<readObject> similar;
  std::vector<readObject> dissimilar;
  uint32_t counter = 0;
  if(setUp.pars_.verbose_){
  	setUp.pars_.seqObj_.seqBase_.outPutSeq(std::cout);
  }
  readVec::allSetLetterCount(inReads);
  setUp.pars_.seqObj_.setLetterCount();

  for (const auto readPos : iter::range(len(inReads))) {
    ++counter;
    if (counter % 50 == 0) {
      std::cout << "On " << counter << " of " << inReads.size() << std::endl;
    }
    if(useNucComp){
    	double sum = inReads[readPos].counter_.getFracDifference(setUp.pars_.seqObj_.counter_,
    			inReads[readPos].counter_.alphabet_);
    	if(sum > maxNucCompDiff){
    		dissimilar.emplace_back(inReads[readPos]);
    		continue;
    	}
    }
    alignerObj.alignCache(setUp.pars_.seqObj_, inReads[readPos], false);
    alignerObj.profilePrimerAlignment(setUp.pars_.seqObj_, inReads[readPos]);
    if(useScore){
      if(alignerObj.parts_.score_ > scoreCutOff){
      	similar.emplace_back(inReads[readPos]);
      }else{
      	dissimilar.emplace_back(inReads[readPos]);
      }
    }else{
      if(alignerObj.comp_.distances_.query_.coverage_ > queryCutoff &&
      		alignerObj.comp_.distances_.percentMatch_ > idCutOff &&
      		alignerObj.comp_.distances_.percentGaps_ < gapCutoff){
      	similar.emplace_back(inReads[readPos]);
      }else{
      	dissimilar.emplace_back(inReads[readPos]);
      }
    }
  }
  counter = 0;
  if(checkComplement){
    for (const auto readPos : iter::range(len(dissimilar))) {
    	auto & read = dissimilar[readPos];
    	read.seqBase_.reverseComplementRead(true);
    	read.setLetterCount();
      ++counter;
      if (counter % 50 == 0) {
        std::cout << "On " << counter << " of " << dissimilar.size() << std::endl;
      }
      if(useNucComp){
      	double sum = read.counter_.getFracDifference(setUp.pars_.seqObj_.counter_,
      			read.counter_.alphabet_);
      	if(sum > maxNucCompDiff){
      		continue;
      	}
      }
      alignerObj.alignCache(setUp.pars_.seqObj_, read, false);
      alignerObj.profilePrimerAlignment(setUp.pars_.seqObj_, read);
      if(useScore){
        if(alignerObj.parts_.score_ > scoreCutOff){
        	similar.emplace_back(read);
        	read.remove = true;
        }
      }else{
        if(alignerObj.comp_.distances_.query_.coverage_ > queryCutoff &&
        		alignerObj.comp_.distances_.percentMatch_ > idCutOff &&
        		alignerObj.comp_.distances_.percentGaps_ < gapCutoff){
        	similar.emplace_back(read);
        	read.remove = true;
        }
      }
    }
    dissimilar = readVecSplitter::splitVectorOnRemove(dissimilar).first;
  }

  SeqOutput::write(similar, setUp.pars_.directoryName_ + "similar", setUp.pars_.ioOptions_);
  SeqOutput::write(dissimilar, setUp.pars_.directoryName_ + "dissimilar", setUp.pars_.ioOptions_);
  std::ofstream outInfo;
  openTextFile(outInfo, setUp.pars_.directoryName_ + "info", ".tab.txt", false, false);
  outInfo << "readsType\treadNum" << std::endl;
  outInfo << "similar\t" << getPercentageString(readVec::getTotalReadCount(similar),
  		readVec::getTotalReadCount(inReads)) << std::endl;
  outInfo << "dissimilar\t" << getPercentageString(readVec::getTotalReadCount(dissimilar),
  		readVec::getTotalReadCount(inReads)) << std::endl;
  if(trimSimilar){
  	std::vector<size_t> startSites;
  	std::vector<size_t> stopSites;
  	std::vector<readObject> dis;
  	for(auto & read : similar){
  		alignerObj.alignCache(setUp.pars_.seqObj_.seqBase_, read.seqBase_, false);
  		size_t firstBase = alignerObj.alignObjectB_.seqBase_.seq_.find_first_not_of("-");
  		startSites.emplace_back(firstBase);
  		size_t lastBase = alignerObj.alignObjectB_.seqBase_.seq_.find_last_not_of("-");
  		uint32_t gapCount = std::count(alignerObj.alignObjectA_.seqBase_.seq_.begin(),
  				alignerObj.alignObjectA_.seqBase_.seq_.begin() + lastBase, '-');
  		stopSites.emplace_back(lastBase - gapCount);
  	}
  	double currentPercentile = 0;
  	size_t minStart = vectorMinimum(startSites);
  	while(currentPercentile < percentile){
  		double currentCount = std::count_if(startSites.begin(), startSites.end(),[&](size_t num){ return num <=minStart;});
  		currentPercentile = currentCount/startSites.size();
  		++minStart;
  	}
  	currentPercentile = 0;
  	size_t maxStop = vectorMaximum(stopSites);
  	while(currentPercentile < percentile){
  		double currentCount = std::count_if(stopSites.begin(), stopSites.end(),[&](size_t num){ return num >=maxStop;});
  		currentPercentile = currentCount/stopSites.size();
  		--maxStop;
  	}
  	for(auto & read : similar){
  		alignerObj.alignCache(setUp.pars_.seqObj_.seqBase_, read.seqBase_, false);
  		size_t firstBase = alignerObj.alignObjectB_.seqBase_.seq_.find_first_not_of("-");
  		size_t lastBase = alignerObj.alignObjectB_.seqBase_.seq_.find_last_not_of("-");
  		uint32_t gapCount = std::count(alignerObj.alignObjectA_.seqBase_.seq_.begin(),
  				alignerObj.alignObjectA_.seqBase_.seq_.begin() + lastBase, '-');
  		lastBase = lastBase - gapCount;
  		if(firstBase > minStart || lastBase < maxStop){
  			dis.emplace_back(read);
  			read.remove = true;
  		}else{
  			uint32_t startPos = alignerObj.getAlignPosForSeqAPos(minStart);
  			uint32_t gapCount = std::count(alignerObj.alignObjectB_.seqBase_.seq_.begin(),
  					alignerObj.alignObjectB_.seqBase_.seq_.begin() + startPos + 1, '-' );
  			startPos -= gapCount;
  			uint32_t stopPos = alignerObj.getAlignPosForSeqAPos(maxStop);
  			uint32_t gapCountStop = std::count(alignerObj.alignObjectB_.seqBase_.seq_.begin(),
  					alignerObj.alignObjectB_.seqBase_.seq_.begin() + stopPos + 1, '-' );
  			stopPos -= gapCountStop;
  			read.trimBack(stopPos + 1);
  			read.trimFront(startPos);
  		}
  	}
  	std::vector<readObject> trimmed = readVecSplitter::splitVectorOnRemove(similar).first;
  	SeqOutput::write(trimmed, setUp.pars_.directoryName_ + "trimmed_similar", setUp.pars_.ioOptions_);
  	SeqOutput::write(dis, setUp.pars_.directoryName_ + "discarded_similar", setUp.pars_.ioOptions_);
    outInfo << "trimmed\t" << getPercentageString(readVec::getTotalReadCount(trimmed),
    		readVec::getTotalReadCount(similar, true)) << std::endl;
    outInfo << "discarded\t" << getPercentageString(readVec::getTotalReadCount(dis),
    		readVec::getTotalReadCount(similar, true)) << std::endl;
  }
  setUp.logRunTime(std::cout);
  return 0;
}


int seqUtilsSplitRunner::getSimilarSequencesByKDist(const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	bool checkComplement = false;
	double cutOff = 0.70;
	uint32_t kmerLen = 7;

	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(cutOff, "--kmerDistCutOff", "Kmer Similarity Score Cut off");
	setUp.setOption(checkComplement, "--checkComplement", "checkComplement");
	setUp.setOption(kmerLen, "--kLen", "Kmer Length");
	setUp.processKmerLenOptions();
	setUp.processDefaultReader(true);
	setUp.processRefFilename(true);
	setUp.finishSetUp(std::cout);


	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto similarOutOpts = setUp.pars_.ioOptions_;
	if("out" == setUp.pars_.ioOptions_.out_.outFilename_ ){
		similarOutOpts.out_.outFilename_ = njh::files::prependFileBasename(similarOutOpts.firstName_, "similar_");
	}else{
		similarOutOpts.out_.outFilename_ = njh::files::prependFileBasename(similarOutOpts.out_.outFilename_, "similar_");
	}
	auto disSimilarOutOpts = setUp.pars_.ioOptions_;
	if("out" == setUp.pars_.ioOptions_.out_.outFilename_ ){
		disSimilarOutOpts.out_.outFilename_ = njh::files::prependFileBasename(disSimilarOutOpts.firstName_, "dissimilar_");
	}else{
		disSimilarOutOpts.out_.outFilename_ = njh::files::prependFileBasename(disSimilarOutOpts.out_.outFilename_, "dissimilar_");
	}
	SeqOutput simWriter(similarOutOpts);
	SeqOutput disWriter(disSimilarOutOpts);

	simWriter.openOut();
	disWriter.openOut();

	SeqInput refReader(setUp.pars_.refIoOptions_);
	auto refSeqs = refReader.readAllReads<seqInfo>();
	std::vector<std::unique_ptr<seqWithKmerInfo>> refKmerReads;
	for (const auto & seq : refSeqs) {
		refKmerReads.emplace_back(std::make_unique<seqWithKmerInfo>(seq, kmerLen, false));
	}

	seqInfo seq;
	if(checkComplement){
		while (reader.readNextRead(seq)) {
			readVec::handelLowerCaseBases(seq, setUp.pars_.ioOptions_.lowerCaseBases_);
			auto seqKmer = std::make_unique<seqWithKmerInfo>(seq, kmerLen, checkComplement);
			uint32_t forwardWinners = 0;
			uint32_t revWinners = 0;
			for (const auto & refSeq : refKmerReads) {
				auto forDist = refSeq->compareKmers(*seqKmer);
				auto revDist = refSeq->compareKmersRevComp(*seqKmer);
				if(forDist.second > cutOff && revDist.second > cutOff){
					if (forDist.first < revDist.first) {
						++revWinners;
					} else {
						++forwardWinners;
					}
				}else if(forDist.second > cutOff){
					++forwardWinners;
				}else if(revDist.second > cutOff){
					++revWinners;
				}
			}
			if(revWinners + forwardWinners > 0){
				if (revWinners > forwardWinners) {
					seq.reverseComplementRead(true, true);
				}
				simWriter.write(seq);
			}else{
				disWriter.write(seq);
			}
		}
	}else{
		while (reader.readNextRead(seq)) {
			readVec::handelLowerCaseBases(seq, setUp.pars_.ioOptions_.lowerCaseBases_);
			auto seqKmer = std::make_unique<seqWithKmerInfo>(seq, kmerLen, checkComplement);
			uint32_t forwardWinners = 0;
			for (const auto & refSeq : refKmerReads) {
				auto forDist = refSeq->compareKmers(*seqKmer);
				if(forDist.second >= cutOff){
					++forwardWinners;
					break;
				}
			}
			if (forwardWinners >= 1) {
				simWriter.write(seq);
			}else{
				disWriter.write(seq);
			}
		}
	}
	return 0;
}


struct defaultSplitPars{
  bool mark_ = false;
	bool include_ = false;
  SeqIOOptions incOpts_;
  SeqIOOptions excOpts_;

};
void defaultSplitSetUpOptions(seqSetUp & setUp, defaultSplitPars & pars){
	setUp.setOption(pars.include_, "--include", "Switch the exclusion and inclusion, by default most splitters exclude the search criteria");
	setUp.setOption(pars.mark_, "--mark", "Mark the sequence names");
  setUp.processVerbose();
  // input file info
  setUp.processDefaultReader(true);

  pars.excOpts_ = setUp.pars_.ioOptions_;
  pars.incOpts_ = setUp.pars_.ioOptions_;
  bfs::path bName = pars.incOpts_.out_.outFilename_;
  if(pars.incOpts_.out_.outFilename_ == "out"){
  	bName = njh::files::bfs::path(setUp.pars_.ioOptions_.firstName_);
  }
  // bName.replace_extension("");
  bName = njh::files::removeExtension(bName);
  pars.incOpts_.out_.outFilename_ = bName.string() + "_included";
	pars.excOpts_.out_.outFilename_ = bName.string() + "_excluded";


}






int seqUtilsSplitRunner::SeqSplitOnLenBelow(const njh::progutils::CmdArgs & inputCommands) {
	defaultSplitPars dSplitPars;
	uint32_t minLen = 0;
  seqSetUp setUp(inputCommands);
  defaultSplitSetUpOptions(setUp, dSplitPars);

	setUp.setOption(minLen, "--minLen", "Exclude sequence below this Minimum Length", true);


  setUp.finishSetUp(std::cout);

  MultiSeqOutCache<seqInfo> seqOuts;
  seqOuts.addReader("include", dSplitPars.incOpts_);
  seqOuts.addReader("exclude", dSplitPars.excOpts_);

  auto checker = std::make_unique<const ReadCheckerLenAbove>( minLen, dSplitPars.mark_);

  SeqIO reader(setUp.pars_.ioOptions_);
  reader.openIn();
  seqInfo seq;
  while(reader.readNextRead(seq)){
  	checker->checkRead(seq);
  	if(dSplitPars.include_) {
  		seq.on_ = !seq.on_;
  	}
    std::string condition;
    if(seq.on_){
    	condition = "include";
    }else{
    	condition = "exclude";
    }
    seqOuts.add(condition, seq);
  }

  if(setUp.pars_.verbose_){
  	setUp.logRunTime(std::cout);
  }
  return 0;
}

int seqUtilsSplitRunner::SeqSplitOnCount(const njh::progutils::CmdArgs & inputCommands) {
	defaultSplitPars dSplitPars;
	uint32_t count = 0;
  seqSetUp setUp(inputCommands);
  defaultSplitSetUpOptions(setUp, dSplitPars);

	setUp.setOption(count, "--count", "Exclude seqs with counts equal to and below this number", true);


  setUp.finishSetUp(std::cout);

  MultiSeqOutCache<seqInfo> seqOuts;
  seqOuts.addReader("include", dSplitPars.incOpts_);
  seqOuts.addReader("exclude", dSplitPars.excOpts_);


  SeqIO reader(setUp.pars_.ioOptions_);
  reader.openIn();
  seqInfo seq;
  while(reader.readNextRead(seq)){
  	if(seq.cnt_ <= count){
  		seq.on_ = false;
  	}
    std::string condition = "";
    if(seq.on_){
    	condition = "include";
    }else{
    	condition = "exclude";
    }
    seqOuts.add(condition, seq);
  }

  if(setUp.pars_.verbose_){
  	setUp.logRunTime(std::cout);
  }
  return 0;
}


int seqUtilsSplitRunner::SeqSplitOnLenWithinMedianLen(const njh::progutils::CmdArgs & inputCommands) {
	defaultSplitPars dSplitPars;
	uint32_t within = 0;
	double withinMultiplier = 0.15;
  seqSetUp setUp(inputCommands);
  defaultSplitSetUpOptions(setUp, dSplitPars);

	setUp.setOption(within, "--within", "Within");
	setUp.setOption(withinMultiplier, "--withinMultiplier", "withinMultiplier");

  setUp.finishSetUp(std::cout);

  MultiSeqOutCache<seqInfo> seqOuts;
  seqOuts.addReader("include", dSplitPars.incOpts_);
  seqOuts.addReader("exclude", dSplitPars.excOpts_);
	uint32_t length = 0;
	{
		std::vector<uint32_t> readLengths;
	  SeqIO reader(setUp.pars_.ioOptions_);
	  reader.openIn();
	  seqInfo seq;
	  while(reader.readNextRead(seq)){
	  	readLengths.emplace_back(len(seq));
	  }
	  length = std::round(vectorMedianRef(readLengths));
	}

	if(0 == within){
		within = withinMultiplier * length;
	}

	auto checker = std::make_unique<const ReadCheckerLenWithin>(within, length, dSplitPars.mark_);


  SeqIO reader(setUp.pars_.ioOptions_);
  reader.openIn();
  seqInfo seq;
  while(reader.readNextRead(seq)){
  	checker->checkRead(seq);
  	if(dSplitPars.include_) {
  		seq.on_ = !seq.on_;
  	}
    std::string condition;
    if(seq.on_){
    	condition = "include";
    }else{
    	condition = "exclude";
    }
    seqOuts.add(condition, seq);
  }

  if(setUp.pars_.verbose_){
  	setUp.logRunTime(std::cout);
  }
  return 0;
}
int seqUtilsSplitRunner::SeqSplitOnLenWithin(const njh::progutils::CmdArgs & inputCommands) {
	defaultSplitPars dSplitPars;
	uint32_t length = 0;
	uint32_t within = 0;

  seqSetUp setUp(inputCommands);
  defaultSplitSetUpOptions(setUp, dSplitPars);

	setUp.setOption(within, "--within", "Within");
	setUp.setOption(length, "--length", "Target Length", true);

  setUp.finishSetUp(std::cout);

  MultiSeqOutCache<seqInfo> seqOuts;
  seqOuts.addReader("include", dSplitPars.incOpts_);
  seqOuts.addReader("exclude", dSplitPars.excOpts_);

	auto checker = std::make_unique<const ReadCheckerLenWithin>(within, length, dSplitPars.mark_);


  SeqIO reader(setUp.pars_.ioOptions_);
  reader.openIn();
  seqInfo seq;
  while(reader.readNextRead(seq)){
  	checker->checkRead(seq);
  	if(dSplitPars.include_) {
  		seq.on_ = !seq.on_;
  	}
    std::string condition;
    if(seq.on_){
    	condition = "include";
    }else{
    	condition = "exclude";
    }
    seqOuts.add(condition, seq);
  }

  if(setUp.pars_.verbose_){
  	setUp.logRunTime(std::cout);
  }
  return 0;
}




int seqUtilsSplitRunner::SeqSplitOnNameContainsPattern(const njh::progutils::CmdArgs & inputCommands) {
	defaultSplitPars dSplitPars;
	std::string nameContainsPattern;
  seqSetUp setUp(inputCommands);
  setUp.processVerbose();
  // input file info
  setUp.processDefaultReader(true);
	setUp.setOption(dSplitPars.include_, "--include", "Switch the exclusion and inclusion, by default most splitters exclude the search criteria");

  dSplitPars.excOpts_ = setUp.pars_.ioOptions_;
  dSplitPars.incOpts_ = setUp.pars_.ioOptions_;
  bfs::path bName = dSplitPars.incOpts_.out_.outFilename_;
  if(dSplitPars.incOpts_.out_.outFilename_ == "out"){
  	bName = njh::files::bfs::path(setUp.pars_.ioOptions_.firstName_);
  }
  bName.replace_extension("");
  dSplitPars.incOpts_.out_.outFilename_ = bName.string() + "_included";
  dSplitPars.excOpts_.out_.outFilename_ = bName.string() + "_excluded";

	setUp.setOption(nameContainsPattern, "--nameContains", "Exclude if Name Contains this regex pattern", true);

  setUp.finishSetUp(std::cout);

  MultiSeqOutCache<seqInfo> seqOuts;
  seqOuts.addReader("include", dSplitPars.incOpts_);
  seqOuts.addReader("exclude", dSplitPars.excOpts_);

  std::regex pat{nameContainsPattern};

  SeqIO reader(setUp.pars_.ioOptions_);
  reader.openIn();
	seqInfo seq;
	while (reader.readNextRead(seq)) {
		auto searchResults = std::regex_match(seq.name_, pat);
		std::string condition;
		if(dSplitPars.include_) {
			searchResults = !searchResults;
		}
		if (searchResults) {
			condition = "exclude";
		} else {
			condition = "include";
		}

		seqOuts.add(condition, seq);
	}
  if(setUp.pars_.verbose_){
  	setUp.logRunTime(std::cout);
  }
  return 0;
}


int seqUtilsSplitRunner::SeqSplitOnNameContains(const njh::progutils::CmdArgs & inputCommands) {
	defaultSplitPars dSplitPars;
	std::string nameContains;
  seqSetUp setUp(inputCommands);
  defaultSplitSetUpOptions(setUp, dSplitPars);
	setUp.setOption(nameContains, "--nameContains", "Exclude if Name Contains", true);

  setUp.finishSetUp(std::cout);

  MultiSeqOutCache<seqInfo> seqOuts;
  seqOuts.addReader("include", dSplitPars.incOpts_);
  seqOuts.addReader("exclude", dSplitPars.excOpts_);

  auto checker = std::make_unique<const ReadCheckerOnNameContaining>( nameContains, dSplitPars.mark_);

  SeqIO reader(setUp.pars_.ioOptions_);
  reader.openIn();
  seqInfo seq;
  while(reader.readNextRead(seq)){
  	checker->checkRead(seq);
  	if(dSplitPars.include_) {
  		seq.on_ = !seq.on_;
  	}
    std::string condition;
    if(seq.on_){
    	condition = "include";
    }else{
    	condition = "exclude";
    }
    seqOuts.add(condition, seq);
  }

  if(setUp.pars_.verbose_){
  	setUp.logRunTime(std::cout);
  }
  return 0;
}

int seqUtilsSplitRunner::SeqSplitOnSeqContains(const njh::progutils::CmdArgs & inputCommands) {
	defaultSplitPars dSplitPars;
	std::string seqContains;
  uint32_t occurences = 1;

  seqSetUp setUp(inputCommands);
  defaultSplitSetUpOptions(setUp, dSplitPars);
	setUp.setOption(occurences, "--occurences", "Minimum number of times a sequence must occur");
	setUp.setOption(seqContains, "--seqContains", "Exclude if Seq Contains this", true);

  setUp.finishSetUp(std::cout);

  MultiSeqOutCache<seqInfo> seqOuts;
  seqOuts.addReader("include", dSplitPars.incOpts_);
  seqOuts.addReader("exclude", dSplitPars.excOpts_);
	auto checker = std::make_unique<const ReadCheckerOnSeqContaining>( seqContains, occurences, dSplitPars.mark_);

  SeqIO reader(setUp.pars_.ioOptions_);
  reader.openIn();
  seqInfo seq;
  while(reader.readNextRead(seq)){
  	checker->checkRead(seq);
  	if(dSplitPars.include_) {
  		seq.on_ = !seq.on_;
  	}
    std::string condition;
    if(seq.on_){
    	condition = "include";
    }else{
    	condition = "exclude";
    }
    seqOuts.add(condition, seq);
  }

  if(setUp.pars_.verbose_){
  	setUp.logRunTime(std::cout);
  }
  return 0;
}

int seqUtilsSplitRunner::SeqSplitOnLenAbove(const njh::progutils::CmdArgs & inputCommands) {
	defaultSplitPars dSplitPars;
	uint32_t maxLength = 0;
  seqSetUp setUp(inputCommands);
	setUp.setOption(maxLength, "--maxLen", "Exclude sequence with lengths above this maximum Length", true);

  defaultSplitSetUpOptions(setUp, dSplitPars);
  setUp.finishSetUp(std::cout);

  MultiSeqOutCache<seqInfo> seqOuts;
  seqOuts.addReader("include", dSplitPars.incOpts_);
  seqOuts.addReader("exclude", dSplitPars.excOpts_);

  auto   	checker = std::make_unique<const ReadCheckerLenBelow>( maxLength, dSplitPars.mark_);

  SeqIO reader(setUp.pars_.ioOptions_);
  reader.openIn();
  seqInfo seq;
  while(reader.readNextRead(seq)){
  	checker->checkRead(seq);
  	if(dSplitPars.include_) {
  		seq.on_ = !seq.on_;
  	}
    std::string condition;
    if(seq.on_){
    	condition = "include";
    }else{
    	condition = "exclude";
    }
    seqOuts.add(condition, seq);
  }

  if(setUp.pars_.verbose_){
  	setUp.logRunTime(std::cout);
  }
  return 0;
}

int seqUtilsSplitRunner::SeqSplitOnLenBetween(const njh::progutils::CmdArgs & inputCommands) {
	defaultSplitPars dSplitPars;
	uint32_t minLen = 0;
	uint32_t maxLength = 0;
  seqSetUp setUp(inputCommands);
  defaultSplitSetUpOptions(setUp, dSplitPars);
	setUp.setOption(maxLength, "--maxLen", "Exclude sequence with lengths above this maximum Length", true);
	setUp.setOption(minLen, "--minLen", "Exclude sequence below this Minimum Length", true);
	if(minLen >maxLength){
		setUp.failed_ = true;
		setUp.addWarning(njh::pasteAsStr("Minimum length must be less than or equal to max length"));
	}
  setUp.finishSetUp(std::cout);

  MultiSeqOutCache<seqInfo> seqOuts;
  seqOuts.addReader("include", dSplitPars.incOpts_);
  seqOuts.addReader("exclude", dSplitPars.excOpts_);

  auto checker = std::make_unique<const ReadCheckerLenBetween>( maxLength, minLen, dSplitPars.mark_);

  SeqIO reader(setUp.pars_.ioOptions_);
  reader.openIn();
  seqInfo seq;
  while(reader.readNextRead(seq)){
  	checker->checkRead(seq);
  	if(dSplitPars.include_) {
  		seq.on_ = !seq.on_;
  	}
    std::string condition;
    if(seq.on_){
    	condition = "include";
    }else{
    	condition = "exclude";
    }
    seqOuts.add(condition, seq);
  }

  if(setUp.pars_.verbose_){
  	setUp.logRunTime(std::cout);
  }
  return 0;
}

int seqUtilsSplitRunner::SeqSplitOnQualityWindow(const njh::progutils::CmdArgs & inputCommands) {
	defaultSplitPars dSplitPars;
  std::string qualWindowString = "50,5,25";

  seqSetUp setUp(inputCommands);
  defaultSplitSetUpOptions(setUp, dSplitPars);
  setUp.setOption(qualWindowString, "-qualWindow", "Sliding Quality Window, goes WindowSize,WindowStep,Thresdhold", true);

	uint32_t qualWindowLength = 50;
	uint32_t qualWindowStep = 5;
	uint8_t qualWindowThres = 25;
  seqUtil::processQualityWindowString(qualWindowString, qualWindowLength,
  		qualWindowStep, qualWindowThres);

  setUp.finishSetUp(std::cout);

  MultiSeqOutCache<seqInfo> seqOuts;
  seqOuts.addReader("include", dSplitPars.incOpts_);
  seqOuts.addReader("exclude", dSplitPars.excOpts_);
  auto checker = std::make_unique<const ReadCheckerOnQualityWindow>( qualWindowLength, qualWindowStep, qualWindowThres, dSplitPars.mark_);

  SeqIO reader(setUp.pars_.ioOptions_);
  reader.openIn();
  seqInfo seq;
  while(reader.readNextRead(seq)){
  	checker->checkRead(seq);
  	if(dSplitPars.include_) {
  		seq.on_ = !seq.on_;
  	}
    std::string condition;
    if(seq.on_){
    	condition = "include";
    }else{
    	condition = "exclude";
    }
    seqOuts.add(condition, seq);
  }

  if(setUp.pars_.verbose_){
  	setUp.logRunTime(std::cout);
  }
  return 0;
}



int seqUtilsSplitRunner::SeqSplitOnQualityCheck(const njh::progutils::CmdArgs & inputCommands) {
	defaultSplitPars dSplitPars;

  seqSetUp setUp(inputCommands);
  defaultSplitSetUpOptions(setUp, dSplitPars);
  double qualCheckCutOff = .90;
  uint32_t qualCheck = 30;
  setUp.setOption(qualCheckCutOff, "--qualCheckCutOff", "the fraction of the bases that have to be above the --qualCheck flag, ranges 0-1", true );
  setUp.setOption(qualCheck, "--qualCheck", "Quality score Check to count", true);

  setUp.finishSetUp(std::cout);

  MultiSeqOutCache<seqInfo> seqOuts;
  seqOuts.addReader("include", dSplitPars.incOpts_);
  seqOuts.addReader("exclude", dSplitPars.excOpts_);
  auto checker = std::make_unique<const ReadCheckerQualCheck>( qualCheck,qualCheckCutOff , dSplitPars.mark_);

  SeqIO reader(setUp.pars_.ioOptions_);
  reader.openIn();
  seqInfo seq;
  while(reader.readNextRead(seq)){
  	checker->checkRead(seq);
  	if(dSplitPars.include_) {
  		seq.on_ = !seq.on_;
  	}
    std::string condition;
    if(seq.on_){
    	condition = "include";
    }else{
    	condition = "exclude";
    }
    seqOuts.add(condition, seq);
  }

  if(setUp.pars_.verbose_){
  	setUp.logRunTime(std::cout);
  }
  return 0;
}


int seqUtilsSplitRunner::SeqSplitOnNucelotideComp(const njh::progutils::CmdArgs & inputCommands) {
	defaultSplitPars dSplitPars;

  seqSetUp setUp(inputCommands);
  defaultSplitSetUpOptions(setUp, dSplitPars);
  setUp.finishSetUp(std::cout);

  MultiSeqOutCache<seqInfo> seqOuts;
  seqOuts.addReader("include", dSplitPars.incOpts_);
  seqOuts.addReader("exclude", dSplitPars.excOpts_);


	SeqIO nucReader(setUp.pars_.ioOptions_);
	nucReader.openIn();
	seqInfo nucSeq;
	charCounter counter;
	while(nucReader.readNextRead(nucSeq)){
		counter.increaseCountByString(nucSeq.seq_, nucSeq.cnt_);
	}
	counter.resetAlphabet(true);
	counter.setFractions();
	nucReader.closeIn();
	nucReader.openIn();
	std::vector<double> differences;
	while(nucReader.readNextRead(nucSeq)){
		charCounter currentCounter(counter.alphabet_);
		currentCounter.increaseCountByString(nucSeq.seq_, nucSeq.cnt_);
		currentCounter.setFractions(counter.alphabet_);
		differences.push_back(counter.getFracDifference(currentCounter, counter.alphabet_));
	}
	double stdCalc = vectorStandardDeviationSamp(differences);
	double meanCalc = vectorMean(differences);


	auto checker = std::make_unique<const ReadCheckerOnNucComp>( counter, meanCalc + 2 * stdCalc , dSplitPars.mark_);

  SeqIO reader(setUp.pars_.ioOptions_);
  reader.openIn();
  seqInfo seq;
  while(reader.readNextRead(seq)){
  	checker->checkRead(seq);
  	if(dSplitPars.include_) {
  		seq.on_ = !seq.on_;
  	}
    std::string condition;
    if(seq.on_){
    	condition = "include";
    }else{
    	condition = "exclude";
    }
    seqOuts.add(condition, seq);
  }

  if(setUp.pars_.verbose_){
  	setUp.logRunTime(std::cout);
  }
  return 0;
}




} // namespace njhseq
