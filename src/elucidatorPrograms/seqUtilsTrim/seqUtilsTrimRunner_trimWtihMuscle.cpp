/*
 * seqUtilsTrimRunner_trimWtihMuscle.cpp
 *
 *  Created on: Jun 5, 2017
 *      Author: nick
 */

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


namespace njhseq {

int seqUtilsTrimRunner::trimWithMuscleToRefInStreaks(const njh::progutils::CmdArgs & inputCommands) {
	seqUtilsTrimSetUp setUp(inputCommands);
	FullTrimReadsPars pars;
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processIoOptions();
	setUp.processRefFilename(true);
	setUp.setOption(pars.keepOnlyOn, "--keepOnlyOn",
			"Keep Only the Reads that are still on");
	setUp.finishSetUp(std::cout);

	njh::sys::requireExternalProgramThrow("muscle");

	SeqIO reader(setUp.pars_.ioOptions_);
	auto seqs = reader.in_.readAllReads<seqInfo>();
	reader.openOut();

	SeqInput refReader(setUp.pars_.refIoOptions_);
	auto refSeqs = refReader.readAllReads<seqInfo>();

	if(refSeqs.empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error no references sequences read in from " << setUp.pars_.refIoOptions_.firstName_ << "\n";
		throw std::runtime_error{ss.str()};
	}

	if(seqs.empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error no references sequences read in from " << setUp.pars_.refIoOptions_.firstName_ << "\n";
		throw std::runtime_error{ss.str()};
	}

	Muscler muscleRunner;

	bool fail = false;
	std::stringstream ss;
	if (refSeqs.empty()) {
		fail = true;
		ss << __PRETTY_FUNCTION__ << ", error refSeqs is empty " << "\n";
	}

	if (seqs.empty()) {
		fail = true;
		ss << __PRETTY_FUNCTION__ << ", error seqs is emtpy " << "\n";
	}
	if (fail) {
		throw std::runtime_error { ss.str() };
	}
	std::vector<seqInfo> allSeqs;
	for(const auto & refSeq : refSeqs){
		allSeqs.emplace_back(getSeqBase(refSeq));
	}
	for(const auto & inputSeq : seqs){
		allSeqs.emplace_back(getSeqBase(inputSeq));
	}
	muscleRunner.muscleSeqs(allSeqs);

	std::vector<seqInfo> alignedRefs = std::vector<seqInfo>(allSeqs.begin(), allSeqs.begin() + refSeqs.size());


	//at least no gaps for at least 10 bases
	uint32_t streakLenCutOff = 10;
	double baseCutOff = 1;


	//auto totalInputSeqs = len(alignedRefs);
	std::function<bool(const std::shared_ptr<Muscler::AlnPosScore> &)> scorePred = [&baseCutOff](const std::shared_ptr<Muscler::AlnPosScore> & score){
		if(static_cast<double>(score->baseCount_)/score->getSpanningCount() >= baseCutOff){
			return true;
		}else{
			return false;
		}
	};

	auto refStartsStop = muscleRunner.getMAlnStartsAndStops(alignedRefs);
	auto allStartsStop = muscleRunner.getMAlnStartsAndStops(allSeqs);
	//count up each location
	auto scores = Muscler::getPileupCounts(alignedRefs);

	auto streaks = Muscler::getAlignmentStreaksPositions(alignedRefs, scorePred,
			streakLenCutOff);
	auto allScores = Muscler::getPileupCounts(allSeqs);

	if (streaks.empty()) {
		std::cerr << "No streaks passed current filters, try being less stringent"
				<< std::endl;
	} else {
		uint32_t streakCount = 0;
		std::cout << "StreakNumber\tstart\tend\tlength\tseqSpanningToNextStreak\tseqSpanningToNextStreakPerc" << "\n";
		if(streaks.size() > 1){
			for(const auto streakPos : iter::range(streaks.size() - 1)){
				const auto & streak = streaks[streakPos];
				uint32_t toNextStreakCount = 0;
				for (const auto pos : iter::range(refSeqs.size(), allSeqs.size())) {
					if (allStartsStop[pos].start_ < streak.end_
							&& allStartsStop[pos].stop_ >= streaks[streakPos + 1].start_) {
						++toNextStreakCount;
					}
				}
				std::cout << streakCount
						<< "\t" << streak.start_
						<< "\t" << streak.end_
						<< "\t" << streak.getLen()
						<< "\t" << toNextStreakCount
						<< "\t" << getPercentageString(toNextStreakCount, seqs.size())
						<< std::endl;
				++streakCount;
			}
		}
	}

	return 0;
}

template<typename INPUTSEQ, typename REF>
void trimSeqsToMultiAlnRef(std::vector<INPUTSEQ> & inputSeqs,
		const std::vector<REF> & refSeqs, const Muscler::TrimWithMusclePars & pars,
		const std::function<bool(const std::shared_ptr<Muscler::AlnPosScore> &)> scorePred) {
	bool fail = false;
	std::stringstream ss;
	if (refSeqs.empty()) {
		fail = true;
		ss << __PRETTY_FUNCTION__ << ", error refSeqs is empty " << "\n";
	}

	if (inputSeqs.empty()) {
		fail = true;
		ss << __PRETTY_FUNCTION__ << ", error inputSeqs is emtpy " << "\n";
	}
	if (fail) {
		throw std::runtime_error { ss.str() };
	}
	std::vector<seqInfo> allSeqs;
	for (const auto & refSeq : refSeqs) {
		allSeqs.emplace_back(getSeqBase(refSeq));
	}
	for (const auto & inputSeq : inputSeqs) {
		allSeqs.emplace_back(getSeqBase(inputSeq));
	}
	Muscler mt;
	mt.muscleSeqs(allSeqs);

	std::vector<seqInfo> alignedRefs = std::vector<seqInfo>(allSeqs.begin(),
			allSeqs.begin() + refSeqs.size());
	auto tempOpts = SeqIOOptions::genFastaOut(bfs::path("tempAlignmentFile.fasta"));
	tempOpts.out_.overWriteFile_ = true;
	SeqOutput::write(allSeqs,tempOpts );

	//[^\|]\|[^|]
//		uint32_t streakLenCutOff = 3; // at least 3 positions in a row must pass the threshold below
//		double spanningCutOff =    .50; // at least 50% of the ref seqs must start here
//		double baseCutOff =        .50; // at least 50% of the bases at this location must have a base
//		uint32_t hardGapCutOff =   2; //there can't be two gaps here

	auto refStartsStop = Muscler::getMAlnStartsAndStops(alignedRefs);
	//count up each location
	auto scores = Muscler::getPileupCounts(alignedRefs);
	auto streaks = Muscler::getAlignmentStreaksPositions(alignedRefs, scorePred,
			pars.streakLenCutOff);

	std::cout << "ref\tstart\tstop" << std::endl;
	for(const auto refPositionsPos : iter::range(refStartsStop.size())){
		std::cout << getSeqBase(refSeqs[refPositionsPos]).name_
				<< "\t" << refStartsStop[refPositionsPos].start_
				<< "\t" << refStartsStop[refPositionsPos].stop_
				<< std::endl;
	}
	std::cout << std::endl;

	std::cout << "streaks" << std::endl;
	for (const auto & streak : streaks) {
		std::cout << streak.start_ << "\t" << streak.end_ << std::endl;
	}

	if (streaks.empty()) {
		std::cerr << "No streaks passed current filters, try being less stringent"
				<< std::endl;
	} else {
		Muscler::trimAlnSeqsToFirstAndLastStreak(allSeqs, streaks);
		readVec::removeGapsFromReads(allSeqs);
		for (const auto pos : iter::range(refSeqs.size(), allSeqs.size())) {
			const auto inputSeqPos = pos - refSeqs.size();
			getSeqBase(inputSeqs[inputSeqPos]) = allSeqs[pos];
		}
	}
}

template<typename INPUTSEQ, typename REF>
void trimSeqsToMultiAlnRef(std::vector<INPUTSEQ> & inputSeqs,
		const std::vector<REF> & refSeqs,
		const Muscler::TrimWithMusclePars & pars) {
	auto totalInputSeqs = len(refSeqs);
	std::function<bool(const std::shared_ptr<Muscler::AlnPosScore> &)> scorePred = [&pars,&totalInputSeqs](const std::shared_ptr<Muscler::AlnPosScore> & score){
		if(score->getBaseSpannedPerc()  >= pars.baseCutOff &&
				score->gapCount_ <= pars.hardGapCutOff &&
				score->getPercentOfSequencesSpanningPosition(totalInputSeqs) >= pars.spanningCutOff){
			return true;
		}else{
			return false;
		}//" | "
	};
	trimSeqsToMultiAlnRef(inputSeqs, refSeqs,pars, scorePred);
}

int seqUtilsTrimRunner::trimWithMuscleToRef(const njh::progutils::CmdArgs & inputCommands) {
	Muscler::TrimWithMusclePars mtPars;
	seqUtilsTrimSetUp setUp(inputCommands);
	FullTrimReadsPars pars;
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processIoOptions();
	setUp.processRefFilename(true);
	setUp.setOption(mtPars.spanningCutOff, "--spanningCutOff", "Spanning Cut Off");
	setUp.setOption(mtPars.baseCutOff, "--gapCutOff", "gap Cut Off");
	setUp.setOption(mtPars.hardGapCutOff, "--hardGapCutOff", "Hard Gap Cut Off");
	setUp.setOption(mtPars.streakLenCutOff, "--streakLenCutOff", "streak Len Cut Off");
	setUp.setOption(pars.keepOnlyOn, "--keepOnlyOn",
			"Keep Only the Reads that are still on");
	setUp.finishSetUp(std::cout);

	njh::sys::requireExternalProgramThrow("muscle");

	SeqIO reader(setUp.pars_.ioOptions_);
	auto seqs = reader.in_.readAllReads<seqInfo>();
	reader.openOut();

	SeqInput refReader(setUp.pars_.refIoOptions_);
	auto refSeqs = refReader.readAllReads<seqInfo>();

	if(refSeqs.empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error no references sequences read in from " << setUp.pars_.refIoOptions_.firstName_ << "\n";
		throw std::runtime_error{ss.str()};
	}

	if(seqs.empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error no references sequences read in from " << setUp.pars_.refIoOptions_.firstName_ << "\n";
		throw std::runtime_error{ss.str()};
	}

	Muscler muscleRunner;

	//muscleRunner.trimSeqsToMultiAlnRef(seqs, refSeqs, mtPars);
	trimSeqsToMultiAlnRef(seqs, refSeqs, mtPars);
	for(const auto & seq : seqs){
		if(!pars.keepOnlyOn || seq.on_){
			reader.write(seq);
		}
	}

	return 0;
}







int seqUtilsTrimRunner::trimWithMuscleMaxStartMinEnd(const njh::progutils::CmdArgs & inputCommands) {
	seqUtilsTrimSetUp setUp(inputCommands);
	FullTrimReadsPars pars;
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processIoOptions();
	setUp.setOption(pars.keepOnlyOn, "--keepOnlyOn",
			"Keep Only the Reads that are still on");
	setUp.finishSetUp(std::cout);

	njh::sys::requireExternalProgramThrow("muscle");




	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	reader.openOut();
	auto seqs = reader.in_.readAllReads<seqInfo>();

	Muscler muscleRunner;
	auto alnSeqs = muscleRunner.muscleSeqsRet(seqs);

	//first check to make sure all seqs are the same size;
	for(const auto & seq : alnSeqs){
		if(len(seq) != len(alnSeqs.front())){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error, all sequences should be the same size: " << len(seq) << " vs " << len(alnSeqs.front()) << "\n";
		}
	}
	auto startsAndStops = Muscler::getMAlnStartsAndStops(alnSeqs);
	//find max start and min end
	auto maxStart = std::max_element(startsAndStops.begin(), startsAndStops.end(), [](const Muscler::StartStopMALNPos & pos1,
			const Muscler::StartStopMALNPos & pos2){
		if(pos1.start_ > pos2.start_){
			return true;
		}else{
			return false;
		}
	})->start_;
	auto minEnd = std::min_element(startsAndStops.begin(), startsAndStops.end(), [](const Muscler::StartStopMALNPos & pos1,
			const Muscler::StartStopMALNPos & pos2){
		if(pos1.stop_ < pos2.stop_){
			return true;
		}else{
			return false;
		}
	})->stop_;

	if(minEnd <= maxStart){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "Min end: " << minEnd << " is less than or equal to maxStart: " << maxStart << std::endl;
		throw std::runtime_error{ss.str()};
	}

	readVecTrimmer::trimEnds(alnSeqs, maxStart, len(alnSeqs.front()) - minEnd + 1);
	readVec::removeGapsFromReads(alnSeqs);
	for(const auto & seq : seqs){
		reader.write(seq);
	}
	return 0;
}




int seqUtilsTrimRunner::trimWithMuscle(const njh::progutils::CmdArgs & inputCommands) {

	Muscler::TrimWithMusclePars mtPars;
	mtPars.streakLenCutOff = 11;
	mtPars.spanningCutOff = .90;
	mtPars.baseCutOff = .95;
	mtPars.hardGapCutOff = 2;
	seqUtilsTrimSetUp setUp(inputCommands);
	FullTrimReadsPars pars;
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(mtPars.spanningCutOff, "--spanningCutOff", "Spanning Cut Off");
	setUp.setOption(mtPars.baseCutOff, "--baseCutOff", "Base Cut Off");
	setUp.setOption(mtPars.hardGapCutOff, "--hardGapCutOff", "Hard Gap Cut Off");
	setUp.setOption(mtPars.streakLenCutOff, "--streakLenCutOff", "streak Len Cut Off");
	setUp.processIoOptions();
	setUp.setOption(pars.keepOnlyOn, "--keepOnlyOn",
			"Keep Only the Reads that are still on");
	setUp.finishSetUp(std::cout);

	njh::sys::requireExternalProgramThrow("muscle");




	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	reader.openOut();
	auto seqs = reader.in_.readAllReads<seqInfo>();
	auto totalInputSeqs = len(seqs);
	std::function<bool(const std::shared_ptr<Muscler::AlnPosScore> &)> scorePred = [&mtPars,&totalInputSeqs](const std::shared_ptr<Muscler::AlnPosScore> & score){
		if(score->getBaseSpannedPerc() >= mtPars.baseCutOff &&
				score->gapCount_ <= mtPars.hardGapCutOff &&
				score->getPercentOfSequencesSpanningPosition(totalInputSeqs) >= mtPars.spanningCutOff){
			return true;
		}else{
			return false;
		}
	};
//expandOutCollapsedToUnique
	Muscler muscleRunner;
	muscleRunner.trimSeqsByMultipleAlignment(seqs, mtPars, scorePred);

	for(const auto & seq : seqs){
		if(!pars.keepOnlyOn || seq.on_){
			reader.write(seq);
		}
	}
	return 0;
}

} // namespace njhseq
