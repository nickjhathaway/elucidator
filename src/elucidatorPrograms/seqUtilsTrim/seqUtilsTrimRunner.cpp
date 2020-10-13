
//  seqUtilsTrimRunner.cpp
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

    
#include "seqUtilsTrimRunner.hpp"
#include <njhseq/IO/SeqIO/SeqIO.hpp>
#include <njhseq/readVectorManipulation/readVectorHelpers/readVecTrimmer.hpp>


namespace njhseq {

seqUtilsTrimRunner::seqUtilsTrimRunner()
    : njh::progutils::ProgramRunner({
  		addFunc("trimFront", trimFront, false),
			addFunc("trimEnd", trimEnd, false),
			addFunc("trimBeforeSeq", trimBeforeSeq, false),
			addFunc("trimFromSeq", trimFromSeq, false),
			addFunc("trimBetweenSeqs", trimBetweenSeqs, false),
			addFunc("trimToLen", trimToLen, false),
			addFunc("trimEdges", trimEdges, false),
			addFunc("trimToSimilarSeq", trimToSimilarSeq, false),
			addFunc("trimAtFirstBase", trimAtFirstBase, false),
			addFunc("trimAtLastBase", trimAtLastBase, false),
			addFunc("trimAtFirstQual", trimAtFirstQual, false),
			addFunc("trimFromMostProbableSharedKmer", trimFromMostProbableSharedKmer, false),
			addFunc("trimToMostProbableSharedKmer", trimToMostProbableSharedKmer, false),
			addFunc("trimFromMostCommonKmer", trimFromMostCommonKmer, false),
			addFunc("trimToMostCommonKmer", trimToMostCommonKmer, false),
			addFunc("trimBetweenMostCommonKmers", trimBetweenMostCommonKmers, false),

			addFunc("trimWithMuscle", trimWithMuscle, false),
			addFunc("trimWithMuscleToRef", trimWithMuscleToRef, false),

			addFunc("trimLstripBase", trimLstripBase, false),
			addFunc("trimRstripBase", trimRstripBase, false),
			addFunc("trimStripBase", trimStripBase, false),
			addFunc("trimLstripQual", trimLstripQual, false),
			addFunc("trimRstripQual", trimRstripQual, false),
			addFunc("trimStripQual", trimStripQual, false),
			addFunc("trimWithMuscleToRefInStreaks", trimWithMuscleToRefInStreaks, false),
			addFunc("trimWithMuscleMaxStartMinEnd", trimWithMuscleMaxStartMinEnd, false),
			addFunc("trimToPositions", trimToPositions, false),
			addFunc("trimToRefWithGlobalAlignment", trimToRefWithGlobalAlignment, false),
			addFunc("trimToRefWithGlobalAlignmentToRefPositions", trimToRefWithGlobalAlignmentToRefPositions, false),
			addFunc("trimToRefWithGlobalAlignmentToRefMultiplePositions", trimToRefWithGlobalAlignmentToRefMultiplePositions, false),
			addFunc("trimWithSlidingQualityAvgWindow", trimWithSlidingQualityAvgWindow, false),
			addFunc("trimToPositionsForEachName", trimToPositionsForEachName, false),
			addFunc("breakupAtRegexPat", breakupAtRegexPat, false),
			addFunc("trimFrontForTandemRepeat", trimFrontForTandemRepeat, false),
			addFunc("trimEdgesForLowEntropy", trimEdgesForLowEntropy, false),
			addFunc("leftTrimToMeanAlignSite", leftTrimToMeanAlignSite, false),
},
                    "seqUtilsTrim") {}
//



int seqUtilsTrimRunner::trimWithSlidingQualityAvgWindow(const njh::progutils::CmdArgs & inputCommands){
	uint32_t windowStep = 1;
	uint32_t windowSize = 10;
	double threshold = 20;
	uint32_t thresholdStreakNumber = 3;

	seqUtilsTrimSetUp setUp(inputCommands);
	setUp.description_ = "Trim sequences where sliding quality windows fall below a certain thresdhold";

	FullTrimReadsPars pars;
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(windowStep, "--windowStep", "Window Step");
	setUp.setOption(windowSize, "--windowSize", "Window Size");
	setUp.setOption(threshold, "--threshold",    "Threshold for quality average window to fall below");
	setUp.setOption(thresholdStreakNumber, "--thresholdStreakNumber", "The number of adjacent windows to fall below the given threshold to trigger trimming");
	setUp.processIoOptions(VecStr{"--fastq", "--fastqgz"});
	setUp.finishSetUp(std::cout);

	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	seqInfo seq;

	while(reader.readNextRead(seq)){
		uint32_t streakNumber = 0;
		for (auto pos : iter::range<uint32_t>(0, len(seq) - windowSize + 1,
				windowStep)) {
			auto qualMean = std::accumulate(seq.qual_.begin() + pos, seq.qual_.begin() + pos + windowSize, 0.0)/windowSize;
			if(qualMean <=threshold){
				++streakNumber;
			}else{
				streakNumber = 0;
			}
			if(streakNumber == thresholdStreakNumber){
				readVecTrimmer::trimToMaxLength(seq, pos);
				break;
			}
		}
		reader.openWrite(seq);
	}

	if(setUp.pars_.verbose_){
		setUp.logRunTime(std::cout);
	}
	return 0;
}

int seqUtilsTrimRunner::trimToSimilarSeq(const njh::progutils::CmdArgs & inputCommands) {
	seqUtilsTrimSetUp setUp(inputCommands);
	FullTrimReadsPars pars;

	setUp.setUpTrimToSimilarSeq(pars);

  SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	reader.openOut();
	seqInfo seq;
  probabilityProfile profile(pars.kmerLength);

  while(reader.readNextRead(seq)){
  	readVec::handelLowerCaseBases(seq, setUp.pars_.ioOptions_.lowerCaseBases_);
  	//first turn off short sequences if they are too short and then continue on
  	if(len(seq) < pars.maxLength ){
  		seq.on_ = false;
  		continue;
  	}
  	uint32_t startPos = pars.maxLength - pars.windowLength;
  	uint32_t stopPos = std::min<uint32_t>(len(seq) - pars.kmerLength + 1, pars.maxLength + pars.windowLength + 1);
  	for(const auto pos : iter::range(startPos, stopPos)){
  		profile.add(seq.seq_.substr(pos, pars.kmerLength), false);
  	}
  }

  profile.updateProfile();
  uint32_t count = 0;
  reader.closeIn();
  reader.openIn();
  while(reader.readNextRead(seq)){
  	readVec::handelLowerCaseBases(seq, setUp.pars_.ioOptions_.lowerCaseBases_);
  	//first turn off short sequences if they are too short and then continue on
  	if(len(seq) < pars.maxLength ){
  		seq.on_ = false;
  		continue;
  	}
  	++count;
  	uint32_t startPos = pars.maxLength - pars.windowLength;
  	uint32_t stopPos = std::min<uint32_t>(len(seq) - pars.kmerLength + 1, pars.maxLength + pars.windowLength + 1);
  	double bestProb = 0;
  	std::vector<uint32_t> bestPos;
  	for(const auto pos : iter::range(startPos, stopPos)){
  		auto currentProb = roundDecPlaces(profile.getProbabilityOfKmer(seq.seq_.substr(pos, pars.kmerLength)), 10);
  		if(currentProb == bestProb){
  			bestProb = currentProb;
  			bestPos.emplace_back(pos);
  		}else if (currentProb > bestProb){
  			bestProb = currentProb;
  			bestPos.clear();
  			bestPos.emplace_back(pos);
  		}
  	}
  	auto maxLen = std::max_element(bestPos.begin(), bestPos.end());
  	if(setUp.pars_.debug_){
  		reader.out_.writeNoCheck(seq);
  	}
  	seq.trimBack(*maxLen);
  	reader.out_.writeNoCheck(seq);
  	if(setUp.pars_.debug_){
    	std::cout << seq.name_ << std::endl;
    	std::cout << "bestProb: " << bestProb << std::endl;
    	std::cout << "bestPos: " << vectorToString(bestPos, ",") << std::endl;
  	}
  }
  return 0;
}
//
int seqUtilsTrimRunner::trimToPositionsForEachName(const njh::progutils::CmdArgs & inputCommands){
	bfs::path trimPositionsFnp;
	seqUtilsTrimSetUp setUp(inputCommands);
	FullTrimReadsPars pars;
  setUp.processVerbose();
  setUp.processDebug();
  setUp.setOption(trimPositionsFnp, "--trimPositionsFnp", "trim Positions file, have at least 3 columns, 1)name,2)start,3)end, positions are 0-based, each additional column is added as meta to the trimmed output", true);
  setUp.processIoOptions();
  setUp.setOption(pars.keepOnlyOn, "--onlyIfInFile", "Trim Only sequence if in it's in the trim file, otherwise will error out");
  setUp.finishSetUp(std::cout);
	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	njh::files::checkExistenceThrow(trimPositionsFnp);
	reader.openOut();
	table trimPositionsTab(trimPositionsFnp, "\t", true);
	trimPositionsTab.checkForColumnsThrow(VecStr{"name", "start", "end"}, __PRETTY_FUNCTION__);
	struct TrimPosInfo {
		TrimPosInfo(uint32_t start, uint32_t end) :
				start_(start), end_(end) {

		}
		uint32_t start_;
		uint32_t end_;
		MetaDataInName meta_;
	};
	VecStr metaCols;
	for(const auto & col : trimPositionsTab.columnNames_){
		if(!njh::in(col, VecStr{"name", "start", "end"})){
			metaCols.emplace_back(col);
		}
	}
	uint32_t nameColPos = trimPositionsTab.getColPos("name");
	uint32_t endColPos = trimPositionsTab.getColPos("end");
	uint32_t startColPos = trimPositionsTab.getColPos("start");
	std::unordered_map<std::string, std::vector<TrimPosInfo>> trimPositions;

	for(const auto & row : trimPositionsTab.content_){
		TrimPosInfo tInfo(
				njh::StrToNumConverter::stoToNum<uint32_t>(row[startColPos]),
				njh::StrToNumConverter::stoToNum<uint32_t>(row[endColPos]));
		std::string name = row[nameColPos];
		if (tInfo.start_ >= tInfo.end_) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error in processing row "
					<< njh::conToStr(row, ", ") << " start, " << tInfo.start_
					<< " can't be equal to or less than end " << tInfo.end_ << "\n";
			throw std::runtime_error { ss.str() };
		}
		tInfo.meta_.addMeta("start", tInfo.start_);
		tInfo.meta_.addMeta("end", tInfo.end_);

		for(const auto & col : metaCols){
			tInfo.meta_.addMeta(col, row[trimPositionsTab.getColPos(col)]);
		}
		trimPositions[name].emplace_back(tInfo);
	}
	seqInfo seq;

	while (reader.readNextRead(seq)) {
		if (pars.keepOnlyOn && !njh::in(seq.name_, trimPositions)) {
			continue;
		}
		for (const auto & tInfo : njh::mapAt(trimPositions, seq.name_)) {
			if (tInfo.start_ >= len(seq) || tInfo.end_ > len(seq)) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error, can't trim seq " << seq.name_
						<< "with length " << len(seq) << ", to positions: " << tInfo.start_
						<< ":" << tInfo.end_ << "\n";
				throw std::runtime_error { ss.str() };
			}
			seqInfo trimmedSeq = seq.getSubRead(tInfo.start_,
					tInfo.end_ - tInfo.start_);
			trimmedSeq.name_.append(tInfo.meta_.createMetaName());
			reader.write(trimmedSeq);
		}

	}

	if(setUp.pars_.verbose_){
		setUp.logRunTime(std::cout);
	}
	return 0;
}

int seqUtilsTrimRunner::trimToPositions(const njh::progutils::CmdArgs & inputCommands){
	uint32_t start = std::numeric_limits<uint32_t>::max();
	uint32_t end  = std::numeric_limits<uint32_t>::max();
	seqUtilsTrimSetUp setUp(inputCommands);
	FullTrimReadsPars pars;
  setUp.processVerbose();
  setUp.processDebug();
  setUp.setOption(start, "--start", "Start position (zero-based, inclusive)", true);
  setUp.setOption(end, "--end", "End position (zero-based, exclusive)", true);;
  setUp.processIoOptions();
  setUp.finishSetUp(std::cout);
	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	seqInfo seq;

	while(reader.readNextRead(seq)){
		seq = seq.getSubRead(start, end - start);
		if (seq.on_ || !pars.keepOnlyOn) {
			reader.openWrite(seq);
		}
	}

	if(setUp.pars_.verbose_){
		setUp.logRunTime(std::cout);
	}
	return 0;
}

int seqUtilsTrimRunner::trimFront(const njh::progutils::CmdArgs & inputCommands){
	seqUtilsTrimSetUp setUp(inputCommands);
	FullTrimReadsPars pars;
	/**@todo add trimmer specific trimmer */
	setUp.setUpTrimFront(pars);
	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	seqInfo seq;

	while(reader.readNextRead(seq)){
		readVecTrimmer::trimOffForwardBases(seq, pars.numberOfFowardBases);
		if (seq.on_ || !pars.keepOnlyOn) {
			reader.openWrite(seq);
		}
	}

	if(setUp.pars_.verbose_){
		setUp.logRunTime(std::cout);
	}
	return 0;
}

int seqUtilsTrimRunner::trimAtFirstQual(const njh::progutils::CmdArgs & inputCommands){
	seqUtilsTrimSetUp setUp(inputCommands);
	FullTrimReadsPars pars;
	setUp.setUpTrimAtFirstQual(pars);
	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	seqInfo seq;

	while(reader.readNextRead(seq)){
		readVecTrimmer::trimAtFirstQualScore(seq, pars.qual);
		if (seq.on_ || !pars.keepOnlyOn) {
			reader.openWrite(seq);
		}
	}

	if(setUp.pars_.verbose_){
		setUp.logRunTime(std::cout);
	}
	return 0;
}

int seqUtilsTrimRunner::trimAtFirstBase(const njh::progutils::CmdArgs & inputCommands){
	seqUtilsTrimSetUp setUp(inputCommands);
	FullTrimReadsPars pars;
	/**@todo add trimmer specific trimmer */
	setUp.setUpTrimAtFirstBase(pars);
	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	seqInfo seq;

	while(reader.readNextRead(seq)){
		readVecTrimmer::trimAtFirstBase(seq, pars.base);
		if (seq.on_ || !pars.keepOnlyOn) {
			reader.openWrite(seq);
		}
	}

	if(setUp.pars_.verbose_){
		setUp.logRunTime(std::cout);
	}
	return 0;
}


int seqUtilsTrimRunner::trimAtLastBase(const njh::progutils::CmdArgs & inputCommands){
	seqUtilsTrimSetUp setUp(inputCommands);
	FullTrimReadsPars pars;
	/**@todo add trimmer specific trimmer */
	setUp.setUpTrimAtLastBase(pars);
	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	seqInfo seq;

	while(reader.readNextRead(seq)){
		readVecTrimmer::trimAtLastBase(seq, pars.base);
		if (seq.on_ || !pars.keepOnlyOn) {
			reader.openWrite(seq);
		}
	}

	if(setUp.pars_.verbose_){
		setUp.logRunTime(std::cout);
	}
	return 0;
}



int seqUtilsTrimRunner::trimEnd(const njh::progutils::CmdArgs & inputCommands){
	seqUtilsTrimSetUp setUp(inputCommands);
	FullTrimReadsPars pars;
	/**@todo add trimmer specific trimmer */
	setUp.setUpTrimEnd(pars);
	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	seqInfo seq;

	while(reader.readNextRead(seq)){
		readVecTrimmer::trimOffEndBases(seq, pars.numberOfEndBases);
		if (seq.on_ || !pars.keepOnlyOn) {
			reader.openWrite(seq);
		}
	}

	if(setUp.pars_.verbose_){
		setUp.logRunTime(std::cout);
	}

	return 0;
}

int seqUtilsTrimRunner::trimEdges(const njh::progutils::CmdArgs & inputCommands){
	seqUtilsTrimSetUp setUp(inputCommands);
	FullTrimReadsPars pars;
	/**@todo add trimmer specific trimmer */
	setUp.setUpTrimEnds(pars);
	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	seqInfo seq;

	while(reader.readNextRead(seq)){
		readVecTrimmer::trimOffForwardBases(seq, pars.numberOfFowardBases);
		readVecTrimmer::trimOffEndBases(seq, pars.numberOfEndBases);
		if (seq.on_ || !pars.keepOnlyOn) {
			reader.openWrite(seq);
		}
	}

	if(setUp.pars_.verbose_){
		setUp.logRunTime(std::cout);
	}

	return 0;
}

int seqUtilsTrimRunner::trimBeforeSeq(const njh::progutils::CmdArgs & inputCommands){
	FullTrimReadsPars pars;
	seqUtilsTrimSetUp setUp(inputCommands);
	/**@todo add trimmer specific trimmer */
	setUp.setUpTrimBeforeSeq(pars);
	uint64_t maxReadSize = 0;
	seqInfo seq;
	{
		SeqIO reader(setUp.pars_.ioOptions_);
		reader.openIn();
		while(reader.readNextRead(seq)){
			readVec::getMaxLength(seq, maxReadSize);
		}
	}

	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	seqInfo forwardObject("forwardSeq", pars.forwardSeq);
	readVec::getMaxLength(forwardObject, maxReadSize);
	aligner alignObj(maxReadSize, setUp.pars_.gapInfo_,
			setUp.pars_.scoring_, setUp.pars_.colOpts_.alignOpts_.countEndGaps_);
	while(reader.readNextRead(seq)){
		readVecTrimmer::trimBeforeSequence(seq, forwardObject, alignObj,
				pars.allowableErrors, pars.tSeqPars_);
		if (seq.on_ || !pars.keepOnlyOn) {
			reader.openWrite(seq);
		}
	}

	if(setUp.pars_.verbose_){
		setUp.logRunTime(std::cout);
	}
	return 0;
}

int seqUtilsTrimRunner::trimFromSeq(const njh::progutils::CmdArgs & inputCommands){
	FullTrimReadsPars pars;
	seqUtilsTrimSetUp setUp(inputCommands);
	/**@todo add trimmer specific trimmer */
	setUp.setUpTrimFromSeq(pars);
	uint64_t maxReadSize = 0;
	seqInfo seq;
	{
		SeqIO reader(setUp.pars_.ioOptions_);
		reader.openIn();
		while(reader.readNextRead(seq)){
			readVec::getMaxLength(seq, maxReadSize);
		}
	}

	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	seqInfo backObject("backSeq", pars.backSeq);
	readVec::getMaxLength(backObject, maxReadSize);
	aligner alignObj(maxReadSize, setUp.pars_.gapInfo_, setUp.pars_.scoring_,
			setUp.pars_.colOpts_.alignOpts_.countEndGaps_);

	while (reader.readNextRead(seq)) {
		readVecTrimmer::trimAtSequence(seq, backObject, alignObj,
				pars.allowableErrors, pars.tSeqPars_);
		if (seq.on_ || !pars.keepOnlyOn) {
			reader.openWrite(seq);
		}
	}

	if (setUp.pars_.verbose_) {
		setUp.logRunTime(std::cout);
	}

	return 0;
}

int seqUtilsTrimRunner::trimBetweenSeqs(
		const njh::progutils::CmdArgs & inputCommands) {

	FullTrimReadsPars pars;
	seqUtilsTrimSetUp setUp(inputCommands);
	/**@todo add trimmer specific trimmer */
	setUp.setUpTrimBetweenSeqs(pars);
	uint64_t maxReadSize = 0;
	seqInfo seq;
	{
		SeqIO reader(setUp.pars_.ioOptions_);
		reader.openIn();
		while (reader.readNextRead(seq)) {
			readVec::getMaxLength(seq, maxReadSize);
		}
	}

	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	seqInfo forwardObject("forwardSeq", pars.forwardSeq);
	readVec::getMaxLength(forwardObject, maxReadSize);
	seqInfo backObject("backSeq", pars.backSeq);
	readVec::getMaxLength(backObject, maxReadSize);
	aligner alignObj(maxReadSize, setUp.pars_.gapInfo_, setUp.pars_.scoring_,
			setUp.pars_.colOpts_.alignOpts_.countEndGaps_);

	while (reader.readNextRead(seq)) {
		readVecTrimmer::trimBetweenSequences(seq, forwardObject, backObject,
				alignObj, pars.allowableErrors, pars.tSeqPars_);

		if (seq.on_ || !pars.keepOnlyOn) {
			reader.openWrite(seq);
		}
	}

	if (setUp.pars_.verbose_) {
		setUp.logRunTime(std::cout);
	}
	return 0;
}

int seqUtilsTrimRunner::trimToLen(const njh::progutils::CmdArgs & inputCommands){
	seqUtilsTrimSetUp setUp(inputCommands);
	FullTrimReadsPars pars;
	/**@todo add trimmer specific trimmer */
	setUp.setUpTrimToLen(pars);
	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	seqInfo seq;

	while(reader.readNextRead(seq)){
		readVecTrimmer::trimToMaxLength(seq, pars.maxLength);
		if (seq.on_ || !pars.keepOnlyOn) {
			reader.openWrite(seq);
		}
	}

	if(setUp.pars_.verbose_){
		setUp.logRunTime(std::cout);
	}

	return 0;
}



int seqUtilsTrimRunner::trimToMostProbableSharedKmer(const njh::progutils::CmdArgs & inputCommands) {

	FullTrimReadsPars pars;
	pars.initForKSharedTrim();

  seqUtilsTrimSetUp setUp(inputCommands);
  setUp.procesingTrimmingWithSeqsOpts(pars);
  setUp.processDebug();
  setUp.processVerbose();
  setUp.processDefaultReader(true);
  setUp.setOption(pars.windowLength, "--windowLen", "Window length to expand the kmer search at which to trim");
  setUp.setOption(pars.kmerLength,   "--kmerLen",   "Kmer Length");
  setUp.finishSetUp(std::cout);
  SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<readObject>();
  readVec::handelLowerCaseBases(inReads, setUp.pars_.ioOptions_.lowerCaseBases_);
  probabilityProfile profile(pars.kmerLength);

  for(const auto readPos : iter::range(len(inReads))){
  	auto & read = inReads[readPos];
  	//first turn off short sequences if they are too short and then continue on
  	if(len(read) < pars.windowLength ){
  		read.seqBase_.on_ = false;
  		continue;
  	}
  	uint32_t startPos = len(read) - pars.windowLength;
  	uint32_t stopPos = len(read) - pars.kmerLength + 1;
  	for(const auto pos : iter::range(startPos, stopPos)){
  		profile.add(read.seqBase_.seq_.substr(pos, pars.kmerLength), false);
  	}
  }

  profile.updateProfile();
  if(setUp.pars_.debug_){
  	profile.printProfile(std::cout, "\t");
  }
  uint32_t count = 0;

  std::unordered_map<std::string, std::vector<double>> kBestCount;
  for(auto & read : inReads){
  	//skip if the length was too small
  	if(!read.seqBase_.on_){
  		continue;
  	}
		++count;
		uint32_t startPos = len(read) - pars.windowLength;
		uint32_t stopPos = len(read) - pars.kmerLength + 1;
  	double bestProb = 0;
  	std::vector<uint32_t> bestPos;
  	for(const auto pos : iter::range(startPos, stopPos)){
  		auto currentProb = roundDecPlaces(profile.getProbabilityOfKmer(read.seqBase_.seq_.substr(pos, pars.kmerLength)), pars.precision);
  		if(currentProb == bestProb){
  			bestProb = currentProb;
  			bestPos.emplace_back(pos);
  		}else if (currentProb > bestProb){
  			bestProb = currentProb;
  			bestPos.clear();
  			bestPos.emplace_back(pos);
  		}
  	}
		for (const auto & pos : bestPos) {
			kBestCount[read.seqBase_.seq_.substr(pos, pars.kmerLength)].emplace_back(
					profile.getProbabilityOfKmer(
							read.seqBase_.seq_.substr(pos, pars.kmerLength)));
		}
  	//auto maxLen = std::max_element(bestPos.begin(), bestPos.end());
  	//kBestCount[read.seqBase_.seq_.substr(*maxLen, pars.kmerLength)].emplace_back(profile.getProbabilityOfKmer(read.seqBase_.seq_.substr(*maxLen, pars.kmerLength)));
  	if(setUp.pars_.debug_){
    	std::cout << read.seqBase_.name_ << std::endl;
    	std::cout << "bestProb: " << bestProb << std::endl;
    	std::cout << "bestPos: " << vectorToString(bestPos, ",") << std::endl;
  	}
  }
  std::string bestK = "";
  double bestProb = 0;
  uint32_t bestCount = 0;
  if(setUp.pars_.debug_){
  	for (const auto & kBest : kBestCount) {
  		std::cout << kBest.first << "\t" <<kBest.second.size()  << std::endl;
  	}
  }
	for (const auto & kBest : kBestCount) {
		if(kBest.second.size() > bestCount){
			bestK = kBest.first;
			bestCount = kBest.second.size();
			bestProb = profile.getProbabilityOfKmer(kBest.first);
		}else if(kBest.second.size() == bestCount){
			auto prob = profile.getProbabilityOfKmer(kBest.first);
			if(prob > bestProb){
				bestK = kBest.first;
				bestProb = prob;
			}
		}
	}
	uint64_t maxReadSize = 0;
	readVec::getMaxLength(inReads, maxReadSize);
	aligner alignObj(maxReadSize, setUp.pars_.gapInfo_,
			setUp.pars_.scoring_, setUp.pars_.colOpts_.alignOpts_.countEndGaps_);

	seqInfo bestSeq(bestK + ":" + estd::to_string(bestProb), bestK);
  SeqOutput writer(setUp.pars_.ioOptions_);
  writer.openOut();
  pars.tSeqPars_.within_ = pars.windowLength + 5;
  for(auto & read : inReads){
  	//skip if the length was too small
  	if(!read.seqBase_.on_){
  		continue;
  	}
  	if(setUp.pars_.debug_){
  		writer.writeNoCheck(read);
		}
		readVecTrimmer::trimAtSequence(read, bestSeq, alignObj,
				pars.allowableErrors, pars.tSeqPars_);
  	writer.writeNoCheck(read);
  }
  return 0;
}

int seqUtilsTrimRunner::trimFromMostProbableSharedKmer(const njh::progutils::CmdArgs & inputCommands) {

  FullTrimReadsPars pars;
  pars.initForKSharedTrim();
  seqUtilsTrimSetUp setUp(inputCommands);
  setUp.procesingTrimmingWithSeqsOpts(pars);
  setUp.processDebug();
  setUp.processVerbose();
  setUp.processDefaultReader(true);
  setUp.setOption(pars.windowLength, "--windowLen", "Window length to expand the kmer search at which to trim");
  setUp.setOption(pars.kmerLength,   "--kmerLen",   "Kmer Length");
  setUp.finishSetUp(std::cout);
  SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<readObject>();
  readVec::handelLowerCaseBases(inReads, setUp.pars_.ioOptions_.lowerCaseBases_);
  probabilityProfile profile(pars.kmerLength);
  for(const auto readPos : iter::range(len(inReads))){
  	auto & read = inReads[readPos];
  	//first turn off short sequences if they are too short and then continue on
  	if(len(read) < pars.windowLength ){
  		read.seqBase_.on_ = false;
  		continue;
  	}
  	uint32_t startPos = 0;
  	uint32_t stopPos = pars.windowLength - pars.kmerLength + 1;
  	for(const auto pos : iter::range(startPos, stopPos)){
  		profile.add(read.seqBase_.seq_.substr(pos, pars.kmerLength), false);
  	}
  }

  profile.updateProfile();
  if(setUp.pars_.debug_){
  	profile.printProfile(std::cout, "\t");
  }
  uint32_t count = 0;

  std::unordered_map<std::string, std::vector<double>> kBestCount;
  for(auto & read : inReads){
  	//skip if the length was too small
  	if(!read.seqBase_.on_){
  		continue;
  	}
		++count;
  	uint32_t startPos = 0;
  	uint32_t stopPos = pars.windowLength - pars.kmerLength + 1;
  	double bestProb = 0;
  	std::vector<uint32_t> bestPos;
  	for(const auto pos : iter::range(startPos, stopPos)){
  		auto currentProb = roundDecPlaces(profile.getProbabilityOfKmer(read.seqBase_.seq_.substr(pos, pars.kmerLength)), pars.precision);
  		//auto currentProb = profile.getProbabilityOfKmer(read.seqBase_.seq_.substr(pos, pars.kmerLength));
  		if(setUp.pars_.debug_){
  			std::cout << read.seqBase_.seq_.substr(pos, pars.kmerLength) << "\t" << currentProb << std::endl;
  		}
  		if(currentProb == bestProb){
  			bestProb = currentProb;
  			bestPos.emplace_back(pos);
  		}else if (currentProb > bestProb){
  			bestProb = currentProb;
  			bestPos.clear();
  			bestPos.emplace_back(pos);
  		}
  	}
		for (const auto & pos : bestPos) {
			kBestCount[read.seqBase_.seq_.substr(pos, pars.kmerLength)].emplace_back(
					profile.getProbabilityOfKmer(
							read.seqBase_.seq_.substr(pos, pars.kmerLength)));
		}
  	//auto maxLen = std::max_element(bestPos.begin(), bestPos.end());
  	//kBestCount[read.seqBase_.seq_.substr(*maxLen, pars.kmerLength)].emplace_back(profile.getProbabilityOfKmer(read.seqBase_.seq_.substr(*maxLen, pars.kmerLength)));
  	if(setUp.pars_.debug_){
    	std::cout << read.seqBase_.name_ << std::endl;
    	std::cout << "bestProb: " << bestProb << std::endl;
    	std::cout << "bestPos: " << vectorToString(bestPos, ",") << std::endl;
  	}
  }
  std::string bestK = "";
  double bestProb = 0;
  uint32_t bestCount = 0;
  if(setUp.pars_.debug_){
  	for (const auto & kBest : kBestCount) {
  		std::cout << kBest.first << "\t" <<kBest.second.size()  << std::endl;
  	}
  }
	for (const auto & kBest : kBestCount) {
		if(kBest.second.size() > bestCount){
			bestK = kBest.first;
			bestCount = kBest.second.size();
			bestProb = profile.getProbabilityOfKmer(kBest.first);
		}else if(kBest.second.size() == bestCount){
			auto prob = profile.getProbabilityOfKmer(kBest.first);
			if(prob > bestProb){
				bestK = kBest.first;
				bestProb = prob;
			}
		}
	}
	uint64_t maxReadSize = 0;
	readVec::getMaxLength(inReads, maxReadSize);
	aligner alignObj(maxReadSize, setUp.pars_.gapInfo_,
			setUp.pars_.scoring_, setUp.pars_.colOpts_.alignOpts_.countEndGaps_);

	seqInfo bestSeq(bestK + ":" + estd::to_string(bestProb), bestK);
  SeqOutput writer(setUp.pars_.ioOptions_);
  writer.openOut();
  pars.tSeqPars_.within_ = pars.windowLength + 5;
  for(auto & read : inReads){
  	//skip if the length was too small
  	if(!read.seqBase_.on_){
  		continue;
  	}
  	if(setUp.pars_.debug_){
  		writer.writeNoCheck(read);
		}
		readVecTrimmer::trimAtSequence(read, bestSeq, alignObj,
				pars.allowableErrors, pars.tSeqPars_);
  	writer.writeNoCheck(read);
  }
  return 0;
}


int seqUtilsTrimRunner::trimToMostCommonKmer(const njh::progutils::CmdArgs & inputCommands) {

  FullTrimReadsPars pars;
  pars.initForKSharedTrim();
  seqUtilsTrimSetUp setUp(inputCommands);
  setUp.procesingTrimmingWithSeqsOpts(pars);
  setUp.processDebug();
  setUp.processVerbose();
  setUp.processDefaultReader(true);
	if (setUp.pars_.ioOptions_.out_.outFilename_ == "out") {
		auto inputPath = njh::files::bfs::path(setUp.pars_.ioOptions_.firstName_);
		inputPath.filename().replace_extension("");
		setUp.pars_.ioOptions_.out_.outFilename_ = njh::files::prependFileBasename(inputPath, "trimmed_");
	}
  setUp.setOption(pars.windowLength, "--windowLen", "Window length to expand the kmer search at which to trim");
  setUp.setOption(pars.kmerLength,   "--kmerLen",   "Kmer Length");
  setUp.finishSetUp(std::cout);
  SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<seqInfo>();
  readVec::handelLowerCaseBases(inReads,
  		setUp.pars_.ioOptions_.lowerCaseBases_);

	uint64_t maxReadSize = 0;
	readVec::getMaxLength(inReads, maxReadSize);

	aligner alignObj(maxReadSize, setUp.pars_.gapInfo_, setUp.pars_.scoring_,
			setUp.pars_.colOpts_.alignOpts_.countEndGaps_);


	readVecTrimmer::trimToMostCommonKmer(inReads, pars, alignObj);

  SeqOutput writer(setUp.pars_.ioOptions_);
  writer.openOut();
  for(auto & seq : inReads){
  	writer.writeNoCheck(seq);
  }

  return 0;
}



int seqUtilsTrimRunner::trimFromMostCommonKmer(const njh::progutils::CmdArgs & inputCommands) {

  FullTrimReadsPars pars;
  pars.initForKSharedTrim();
  seqUtilsTrimSetUp setUp(inputCommands);
  setUp.procesingTrimmingWithSeqsOpts(pars);
  setUp.processDebug();
  setUp.processVerbose();
  setUp.processDefaultReader(true);
  setUp.setOption(pars.windowLength, "--windowLen", "Window length to expand the kmer search at which to trim");
  setUp.setOption(pars.kmerLength,   "--kmerLen",   "Kmer Length");
  setUp.finishSetUp(std::cout);
  SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<seqInfo>();
  readVec::handelLowerCaseBases(inReads,
  		setUp.pars_.ioOptions_.lowerCaseBases_);

	uint64_t maxReadSize = 0;
	readVec::getMaxLength(inReads, maxReadSize);

	aligner alignObj(maxReadSize, setUp.pars_.gapInfo_,
			setUp.pars_.scoring_, setUp.pars_.colOpts_.alignOpts_.countEndGaps_);


	readVecTrimmer::trimFromMostCommonKmer(inReads, pars, alignObj);

  SeqOutput writer(setUp.pars_.ioOptions_);
  writer.openOut();
  for(auto & seq : inReads){
  	writer.writeNoCheck(seq);
  }
  return 0;
}

//
int seqUtilsTrimRunner::trimBetweenMostCommonKmers(const njh::progutils::CmdArgs & inputCommands) {

  FullTrimReadsPars pars;
  pars.initForKSharedTrim();
  seqUtilsTrimSetUp setUp(inputCommands);
  setUp.procesingTrimmingWithSeqsOpts(pars);
  setUp.processDebug();
  setUp.processVerbose();
  setUp.processDefaultReader(true);
  setUp.setOption(pars.windowLength, "--windowLen", "Window length to expand the kmer search at which to trim");
  setUp.setOption(pars.kmerLength,   "--kmerLen",   "Kmer Length");
  setUp.finishSetUp(std::cout);
  SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<seqInfo>();
  readVec::handelLowerCaseBases(inReads,
  		setUp.pars_.ioOptions_.lowerCaseBases_);

	uint64_t maxReadSize = 0;
	readVec::getMaxLength(inReads, maxReadSize);

	aligner alignObj(maxReadSize, setUp.pars_.gapInfo_,
			setUp.pars_.scoring_, setUp.pars_.colOpts_.alignOpts_.countEndGaps_);


	readVecTrimmer::trimBetweenMostCommonKmers(inReads, pars, alignObj);

  SeqOutput writer(setUp.pars_.ioOptions_);
  writer.openOut();
  for(auto & seq : inReads){
  	writer.writeNoCheck(seq);
  }
  return 0;
}



} // namespace njhseq
