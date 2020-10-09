/*
 * genExp_findOutliersWithMuscle.cpp
 *
 *  Created on: Jul 16, 2017
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

#include <TwoBit.h>
#include "genExp.hpp"
#include "elucidator/BamToolsUtils.h"
#include "elucidator/objects/BioDataObject.h"
#include "elucidator/objects/dataContainers/graphs.h"
#include <njhseq/GenomeUtils.h>


namespace njhseq {


int genExpRunner::createAgreementSegmentsWithMuscleRef(const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	OutOptions outOpts;
	outOpts.outFilename_ = "";
	outOpts.outExtention_ = ".tsv";
	uint32_t minSegmentLength = 50;
	uint32_t allowableErrors = 2;

	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(minSegmentLength, "--minSegmentLength", "Minimum segment length agreement");
	setUp.setOption(allowableErrors, "--allowableErrors", "allowable Errors to occur in a segment");
	setUp.processReadInNames(true);
	setUp.processWritingOptions(outOpts);
	setUp.processRefFilename(true);
	setUp.finishSetUp(std::cout);

	njh::sys::requireExternalProgramThrow("muscle");

	SeqInput reader(setUp.pars_.ioOptions_);
	auto seqs = reader.readAllReads<seqInfo>();

	OutputStream out(outOpts);

	SeqInput refReader(setUp.pars_.refIoOptions_);
	auto refSeqs = refReader.readAllReads<seqInfo>();
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

	Muscler muscleRunner;



	std::vector<seqInfo> allSeqs;
	for(const auto & refSeq : refSeqs){
		allSeqs.emplace_back(getSeqBase(refSeq));
	}
	for(const auto & inputSeq : seqs){
		allSeqs.emplace_back(getSeqBase(inputSeq));
	}
	muscleRunner.muscleSeqs(allSeqs);

	std::vector<seqInfo> alignedRefs = std::vector<seqInfo>(allSeqs.begin(), allSeqs.begin() + refSeqs.size());

	out << "seqName\trefName\tstart\tsize";

	out << std::endl;

	for(const auto pos : iter::range(refSeqs.size(), allSeqs.size())){
		for(const auto refPos : iter::range(refSeqs.size())){
			uint32_t currentSegLength = 0;
			uint32_t segstart = std::numeric_limits<uint32_t>::max();
			uint32_t errorsEncountered = 0;
			uint32_t basesSinceLastError = 0;
			for(const auto seqPos : iter::range(allSeqs[pos].seq_.size())){
				if(allSeqs[pos].seq_[seqPos] == allSeqs[refPos].seq_[seqPos]){
					++currentSegLength;
					if(std::numeric_limits<uint32_t>::max() == segstart ){
						segstart = seqPos;
					}
					++basesSinceLastError;
					if(basesSinceLastError > minSegmentLength){
						errorsEncountered = 0;
					}
				}else{
					basesSinceLastError = 0;
					++errorsEncountered;
					if(errorsEncountered > allowableErrors){
						if(currentSegLength >= minSegmentLength){
							out << allSeqs[pos].name_
									<< "\t" << allSeqs[refPos].name_
									<< "\t" << segstart
									<< "\t" << currentSegLength
							<< std::endl;
						}
						//reset counts
						currentSegLength = 0;
						segstart = std::numeric_limits<uint32_t>::max();
					}else{
						++currentSegLength;
						if(std::numeric_limits<uint32_t>::max() == segstart ){
							segstart = seqPos;
						}
					}
				}
			}
			if(currentSegLength >= minSegmentLength){
				out << allSeqs[pos].name_
						<< "\t" << allSeqs[refPos].name_
						<< "\t" << segstart
						<< "\t" << currentSegLength
				<< std::endl;
			}
		}
	}
	return 0;
}

int genExpRunner::createAgreementMatrixWithMuscleRef(const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	OutOptions outOpts;
	outOpts.outFilename_ = "";
	outOpts.outExtention_ = ".tsv";
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processReadInNames(true);
	setUp.processWritingOptions(outOpts);
	setUp.processRefFilename(true);
	setUp.finishSetUp(std::cout);

	njh::sys::requireExternalProgramThrow("muscle");

	SeqInput reader(setUp.pars_.ioOptions_);
	auto seqs = reader.readAllReads<seqInfo>();

	OutputStream out(outOpts);

	SeqInput refReader(setUp.pars_.refIoOptions_);
	auto refSeqs = refReader.readAllReads<seqInfo>();
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

	Muscler muscleRunner;



	std::vector<seqInfo> allSeqs;
	for(const auto & refSeq : refSeqs){
		allSeqs.emplace_back(getSeqBase(refSeq));
	}
	for(const auto & inputSeq : seqs){
		allSeqs.emplace_back(getSeqBase(inputSeq));
	}
	muscleRunner.muscleSeqs(allSeqs);

	std::vector<seqInfo> alignedRefs = std::vector<seqInfo>(allSeqs.begin(), allSeqs.begin() + refSeqs.size());

	out << "seqName\trefName";
	for(const auto pos : iter::range(allSeqs.front().seq_.size())){
		out << "\t" << pos;
	}
	out << std::endl;

	for(const auto pos : iter::range(refSeqs.size(), allSeqs.size())){
		for(const auto refPos : iter::range(refSeqs.size())){
			std::vector<uint32_t> refAgreement;
			for(const auto seqPos : iter::range(allSeqs[pos].seq_.size())){
				if(allSeqs[pos].seq_[seqPos] == allSeqs[refPos].seq_[seqPos]){
					refAgreement.emplace_back(1);
				}else{
					refAgreement.emplace_back(0);
				}
			}
			out << allSeqs[pos].name_
					<< "\t" << allSeqs[refPos].name_
					<< "\t" << njh::conToStr(refAgreement, "\t")
			<< std::endl;
		}
	}
	return 0;
}

int genExpRunner::printMlnScores(const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	OutOptions outOpts;
	outOpts.outFilename_ = "";
	outOpts.outExtention_ = ".tsv";
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processReadInNames(true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	njh::sys::requireExternalProgramThrow("muscle");

	SeqInput reader(setUp.pars_.ioOptions_);
	auto seqs = reader.readAllReads<seqInfo>();
	Muscler::checkAlignSeqsLensThrow(seqs, __PRETTY_FUNCTION__);
	OutputStream out(outOpts);

	auto scores = Muscler::getPileupCounts(seqs);
	Muscler::writeScores(out,seqs, scores);

	return 0;
}

int genExpRunner::printMlnStreaks(const njh::progutils::CmdArgs & inputCommands) {
	OutOptions outOpts;
	outOpts.outExtention_ = ".bed";
	uint32_t streakLenCutOff = 10;
	double baseCutOff = 1;
	bool muscleFirst = false;
	bfs::path muscleWorkingDir = "/tmp/";
	bool keepMuscleTempFiles = false;
	seqSetUp setUp(inputCommands);
	setUp.setOption(keepMuscleTempFiles, "--keepMuscleTempFiles", "Keep Muscle Temp Files");
	setUp.setOption(muscleWorkingDir, "--muscleWorkingDir", "muscle Working Dir");

	setUp.setOption(muscleFirst, "--muscleFirst", "Muscle seqs First if they aren't already aligned");
	setUp.setOption(baseCutOff, "--baseCutOff", "The fraction of seqs with bases here");
	setUp.setOption(streakLenCutOff, "--streakLenCutOff", "Streak Length Cut Off");

	outOpts.outFilename_ = "";
	outOpts.outExtention_ = ".tsv";
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processReadInNames(true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	njh::sys::requireExternalProgramThrow("muscle");

	SeqInput reader(setUp.pars_.ioOptions_);
	auto seqs = reader.readAllReads<seqInfo>();
	Muscler mRunner;
	mRunner.keepTemp_ = keepMuscleTempFiles;
	mRunner.workingPath_ = muscleWorkingDir;
	if(muscleFirst){
		mRunner.muscleSeqs(seqs);
	}
	Muscler::checkAlignSeqsLensThrow(seqs, __PRETTY_FUNCTION__);
	OutputStream out(outOpts);

	//at least no gaps for at least 10 bases


	std::function<bool(const std::shared_ptr<Muscler::AlnPosScore> &)> scorePred = [&baseCutOff](const std::shared_ptr<Muscler::AlnPosScore> & score){
		if(static_cast<double>(score->baseCount_)/score->getSpanningCount() >= baseCutOff){
			return true;
		}else{
			return false;
		}
	};

	auto alnStartsStop = mRunner.getMAlnStartsAndStops(seqs);
	//count up each location
	auto scores = Muscler::getPileupCounts(seqs);

	auto streaks = Muscler::getAlignmentStreaksPositions(seqs, scorePred,
			streakLenCutOff);

	auto allScores = Muscler::getPileupCounts(seqs);

	if (streaks.empty()) {
		std::cerr << "No streaks passed current filters, try being less stringent"
				<< std::endl;
	} else {
		auto names = readVec::getNames(seqs);
		for(const auto streakPos : iter::range(streaks.size())){
			const auto & streak = streaks[streakPos];
			for(const auto namePos : iter::range(names.size())){
				auto seqStart =  getRealPosForAlnPos(seqs[namePos].seq_, streak.start_);
				auto seqStop = getRealPosForAlnPos(seqs[namePos].seq_, streak.end_ - 1) + 1;
				out << names[namePos]
						<< "\t" << seqStart
						<< "\t" << seqStop
						<< "\t" << "streak" << streakPos << "-" << streak.start_ << "-" << streak.end_
						<< "\t" << seqStop - seqStart
						<< "\t" << '+'
						<< std::endl;
			}
		}
	}
	return 0;
}



int genExpRunner::findOutliersWithMuscleToRefs(const njh::progutils::CmdArgs & inputCommands) {
	FullTrimReadsPars pars;
	OutOptions outOpts;
	double cutOff = 0.80;
	bfs::path muscleWorkingDir = "/tmp/";
	bool keepMuscleTempFiles = false;

	seqSetUp setUp(inputCommands);
	setUp.setOption(keepMuscleTempFiles, "--keepMuscleTempFiles", "Keep Muscle Temp Files");
	setUp.setOption(muscleWorkingDir, "--muscleWorkingDir", "muscle Working Dir");
	setUp.setOption(cutOff, "--cutOff", "Write out seqs above this cut off");
	outOpts.outFilename_ = "";
	outOpts.outExtention_ = ".tsv";
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processReadInNames(true);
	setUp.processWritingOptions(outOpts);
	setUp.processRefFilename(true);
	setUp.setOption(pars.keepOnlyOn, "--keepOnlyOn",
			"Keep Only the Reads that are still on");
	setUp.finishSetUp(std::cout);

	njh::sys::requireExternalProgramThrow("muscle");

	SeqInput reader(setUp.pars_.ioOptions_);
	auto seqs = reader.readAllReads<seqInfo>();
	auto outSeqOpts = SeqIOOptions::genFastaOut(outOpts.outFilename_.string() + "_malnSeqs");
	outSeqOpts.out_.transferOverwriteOpts(outOpts);
	SeqOutput writer(outSeqOpts);

	writer.openOut();

	auto filteredSeqOpts = SeqIOOptions::genFastaOut(outOpts.outFilename_);
	filteredSeqOpts.out_.transferOverwriteOpts(outOpts);
	SeqOutput filtered_writer(filteredSeqOpts);
	if(pars.keepOnlyOn){
		filtered_writer.openOut();
	}

	auto filteredOffSeqOpts = SeqIOOptions::genFastaOut(njh::files::nameAppendBeforeExt(outOpts.outFilename_, "_filteredOff"));
	filteredOffSeqOpts.out_.transferOverwriteOpts(outOpts);
	SeqOutput filteredOff_writer(filteredOffSeqOpts);
	if(pars.keepOnlyOn){
		filteredOff_writer.openOut();
	}


	OutputStream out(outOpts);

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
	muscleRunner.keepTemp_ = keepMuscleTempFiles;
	muscleRunner.workingPath_ = muscleWorkingDir;

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


	uint32_t streakLenCutOff = 3;
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
		std::vector<seqInfo> aboveCutOff;
		std::vector<seqInfo> belowCutOff;
		out << "seqName\tscore\tspanningSpots\tadjustedScore\tnonZeroSpanningSpots\tfractionNonZeroSpanning\tspanningSpotsNoGaps\tadjustedScoreNoGaps\tnonZeroSpanningSpotsNoGaps\tfractionNonZeroSpanningNoGaps" << "\n";
		auto alpha = determineAlphabet(allSeqs);
		std::unordered_map<uint32_t, std::unique_ptr<charCounter>> profiles;
		for (const auto & streak : streaks) {
			for(uint32_t pos : iter::range(streak.start_, streak.end_)){
				profiles[pos] = std::make_unique<charCounter>(alpha);
				for(const char c : alpha){
					profiles[pos]->increaseCountOfBase(c);
				}
				for(const auto & ref : alignedRefs){
					profiles[pos]->increaseCountOfBase(ref.seq_[pos]);
				}
				profiles[pos]->setFractions();
			}
		}

		for (const auto alnSeqPos : iter::range(refSeqs.size(), allSeqs.size())) {
			uint32_t spanningSpots = 0;
			uint32_t nonZeroSpanningSpots = 0;
			double score = 0;
			uint32_t spanningSpotsNoGaps = 0;
			uint32_t nonZeroSpanningSpotsNoGaps = 0;
			double scoreNoGaps = 0;
			for (const auto & streak : streaks) {
				for(uint32_t pos : iter::range(streak.start_, streak.end_)){
					if(pos>= allStartsStop[alnSeqPos].start_ && pos <= allStartsStop[alnSeqPos].stop_){
						if(0 == spanningSpots ){
							score = 1;
						}
						if(1 != profiles[pos]->chars_[allSeqs[alnSeqPos].seq_[pos]]){
							++nonZeroSpanningSpots;
						}
						score *= profiles[pos]->fractions_[allSeqs[alnSeqPos].seq_[pos]];
						++spanningSpots;
						if('-' != allSeqs[alnSeqPos].seq_[pos]){
							if(0 == spanningSpotsNoGaps ){
								scoreNoGaps = 1;
							}
							if(1 != profiles[pos]->chars_[allSeqs[alnSeqPos].seq_[pos]]){
								++nonZeroSpanningSpotsNoGaps;
							}
							score *= profiles[pos]->fractions_[allSeqs[alnSeqPos].seq_[pos]];
							++spanningSpotsNoGaps;
						}
					}
				}
			}
			out << allSeqs[alnSeqPos].name_
					<< "\t" << score
					<< "\t" << spanningSpots
					<< "\t" << (0 == spanningSpots ? 0 :score/spanningSpots)
					<< "\t" << nonZeroSpanningSpots
					<< "\t" << (0 == spanningSpots ? 0 :static_cast<double>(nonZeroSpanningSpots)/spanningSpots)
					<< "\t" << spanningSpotsNoGaps
					<< "\t" << (0 == spanningSpotsNoGaps ? 0 :scoreNoGaps/spanningSpotsNoGaps)
					<< "\t" << nonZeroSpanningSpotsNoGaps
					<< "\t" << (0 == spanningSpotsNoGaps ? 0 :static_cast<double>(nonZeroSpanningSpotsNoGaps)/spanningSpotsNoGaps)
					<< std::endl;
			if((0 == spanningSpotsNoGaps ? 0 :static_cast<double>(nonZeroSpanningSpotsNoGaps)/spanningSpotsNoGaps) >=cutOff){
				aboveCutOff.emplace_back(seqs[alnSeqPos - refSeqs.size()]);
			}else{
				belowCutOff.emplace_back(seqs[alnSeqPos - refSeqs.size()]);
			}
		}
		if(pars.keepOnlyOn){
			filtered_writer.write(aboveCutOff);
			filteredOff_writer.write(belowCutOff);
		}
	}

	writer.write(allSeqs);

	return 0;
}

}  // namespace njhseq

