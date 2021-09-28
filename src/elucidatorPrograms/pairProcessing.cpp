/*
 * pairProcessing.cpp
 *
 *  Created on: Nov 12, 2018
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

#include "pairProcessing.hpp"
#include "elucidator/objects/BioDataObject.h"
#include "elucidator/objects/counters/DNABaseCounter.hpp"
#include <SeekDeep/objects/IlluminaUtils/PairedReadProcessor.hpp>
#include <njhseq/objects/helperObjects/motif.hpp>
#include <njhseq/readVectorManipulation/readVectorHelpers/readVecTrimmer.hpp>


namespace njhseq {




pairProcessingRunner::pairProcessingRunner()
    : njh::progutils::ProgramRunner(
          {
					 addFunc("detectPossiblePrimers", detectPossiblePrimers, false),
					 addFunc("StitchPairedReads", StitchPairedReads, false),
           },
          "pairProcessing") {}



int pairProcessingRunner::StitchPairedReads(
		const njh::progutils::CmdArgs & inputCommands) {
	PairedReadProcessor::ProcessParams params;

	TableIOOpts tabOpts = TableIOOpts::genTabFileOut(bfs::path(""), true);
	OutOptions outOpts(bfs::path("out"));
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	params.verbose_ = setUp.pars_.verbose_;
	setUp.processDebug();
	params.debug_ = setUp.pars_.debug_;
	setUp.pars_.gap_ = "10,1";
	setUp.pars_.gapInfo_.gapOpen_ = 10;
	setUp.pars_.gapInfo_.gapExtend_ = 1;
	setUp.processGap();
	setUp.pars_.caseInsensitiveScoring_ = true;
	setUp.pars_.lessNScoring_ = true;
	setUp.processScoringPars();
	setUp.setOption(params.testNumber_, "--testNumber", "Test Number");
	setUp.setOption(params.hardMismatchCutOff_, "--hardMismatchCutOff", "Hard Mismatch Cut Off, also don't allow this many mismatches");
	setUp.setOption(params.lqMismatchCutOff,    "--lqMismatchCutOff",   "Low qaulity Mismatch Cut Off, also don't allow this many mismatches");
	setUp.setOption(params.hqMismatchCutOff,    "--hqMismatchCutOff",   "High quality Mismatch Cut Off, also don't allow this many mismatches");

	setUp.setOption(params.minOverlap_, "--minOverlap", "Minimum overlap");
	setUp.setOption(params.writeOverHangs_, "--writeOverHangs", "Write Over Hangs");

	setUp.setOption(params.errorAllowed_, "--errorAllowed", "Percent Error Allowed, between 0 and 1");
	if(params.errorAllowed_ > 1 || params.errorAllowed_ < 0){
		setUp.failed_ = true;
		std::stringstream ss;
		ss << "Error in setting errorAllowed, error should be between 0 and 1, not " << params.errorAllowed_ << "\n";
		setUp.addWarning(ss.str());
	}

	setUp.setOption(params.qualWindowPar_.avgQualCutOff_, "--qWindowTrimAvgQualCutOff", "Quality Window Trim Avg Qual Cut Off");
	setUp.setOption(params.qualWindowPar_.windowSize_, "--qWindowSize", "Quality Window Trim Size");
	setUp.setOption(params.qualWindowPar_.windowStep_, "--qWindowStep", "Quality Window Trim Step");
	bool noTrimLowQualWindows = false;
	setUp.setOption(noTrimLowQualWindows, "--noTrimLowQualWindows", "Don't Trim Low Qual Windows");
	params.trimLowQaulWindows_ = !noTrimLowQualWindows;


	setUp.processReadInNames(VecStr { "--fastq1", "--fastq2", "--fastq1gz", "--fastq2gz" }, true);
	setUp.processWritingOptions(outOpts);
	setUp.setOption(setUp.pars_.ioOptions_.revComplMate_, "--revCompMate", "revCompMate");
	setUp.finishSetUp(std::cout);

	PairedReadProcessor processor(params);

	uint64_t maxsize = processor.guessMaxReadLenFromFile(setUp.pars_.ioOptions_);

	auto alnGapPars = gapScoringParameters(
			setUp.pars_.gapInfo_.gapOpen_,
			setUp.pars_.gapInfo_.gapExtend_,
			0,0,
			0,0);
	//std::cout << alnGapPars.toJson() << std::endl;

	aligner alignerObj(maxsize, alnGapPars, setUp.pars_.scoring_, false);
	alignerObj.qScorePars_.qualThresWindow_ = 0;

	PairedReadProcessor::ProcessorOutWriters writers(outOpts);

	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();

	auto processResults = processor.processPairedEnd(reader, writers, alignerObj);

	if(setUp.pars_.verbose_){
		std::cout << processResults.toJsonCounts() << std::endl;
	}
	return 0;
}

int pairProcessingRunner::detectPossiblePrimers(
		const njh::progutils::CmdArgs & inputCommands) {
	bool writeOutBaseQualCompositionProfile = false;
	bool keepOverHangs = false;
	uint32_t consensusCountCutOff = 5;
	uint32_t frontCheck = 7;
	uint32_t minOverlap = 10;

	uint32_t minOverhangLen = 5;

	double entropyCutOff = 0.80;

	double errorAllowed = 0.01;
	uint32_t hardMismatchCutOff = 10;
	uint32_t checkAmount = 10000;
	double fracCutOff = 0.01;
	uint32_t qualCutOff = 12;
	TableIOOpts tabOpts = TableIOOpts::genTabFileOut(bfs::path(""), true);
	OutOptions outOpts(bfs::path("out"));
	uint32_t testNumber = 20000000;//std::numeric_limits<uint32_t>::max();
	uint32_t gatheredNumber = 40000; //std::numeric_limits<uint32_t>::max();
	uint32_t outPrimerLengthMax = 70;
	uint32_t outPrimerLengthMin = 7;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.pars_.gap_ = "10,1";
	setUp.pars_.gapInfo_.gapOpen_ = 10;
	setUp.pars_.gapInfo_.gapExtend_ = 1;
	setUp.processGap();
	setUp.pars_.caseInsensitiveScoring_ = true;
	setUp.pars_.lessNScoring_ = true;
	setUp.processScoringPars();//
	setUp.setOption(writeOutBaseQualCompositionProfile, "--writeOutBaseQualCompositionProfile", "Write Out Base Qual Composition Profile");
	setUp.setOption(keepOverHangs, "--keepOverHangs", "Keep the overhang files");
	setUp.setOption(qualCutOff, "--qualCutOff", "Quality Cut Off for base to be considered for consensus");
	setUp.setOption(consensusCountCutOff, "--consensusCountCutOff", "Consensus Count Cut Off to be considered for the consensus building");
	setUp.setOption(fracCutOff, "--fracCutOff", "Frac Cut off to be considered for possible primer investigation");
	setUp.setOption(frontCheck, "--frontCheck", "How much of the front of the possible primers to use for grouping");
	setUp.setOption(testNumber, "--testNumber", "The total number of sequences to try to get overhangs from, once this number is hit no matter how many overhangs there are no more will be collected");
	setUp.setOption(gatheredNumber, "--gatheredNumber", "The number of gathered overhangs needed, once this number is hit no more overhangs will be collected");
	setUp.setOption(outPrimerLengthMax, "--outPrimerLengthMax","Output Primer Maximum length to output");
	setUp.setOption(outPrimerLengthMin, "--outPrimerLengthMin","Output Primer Minimum length to output");

	setUp.setOption(hardMismatchCutOff, "--hardMismatchCutOff", "Hard Mismatch Cut Off, also don't allow this many mismatches");
	setUp.setOption(minOverlap, "--minOverlap", "Minimum overlap");
	setUp.setOption(minOverhangLen, "--minOverhangLen", "min Overhang Len");
	setUp.setOption(entropyCutOff, "--entropyCutOff", "Entropy of the overhang Cut Off");



	setUp.setOption(errorAllowed, "--errorAllowed", "Percent Error Allowed, between 0 and 1");
	if(errorAllowed > 1 || errorAllowed < 0){
		setUp.failed_ = true;
		std::stringstream ss;
		ss << "Error in setting errorAllowed, error should be between 0 and 1, not " << errorAllowed << "\n";
		setUp.addWarning(ss.str());
	}
	double percentId = 1 - errorAllowed;
	setUp.processReadInNames(VecStr { "--fastq1", "--fastq2", "--fastq1gz", "--fastq2gz" }, true);
	setUp.processWritingOptions(outOpts);
	bool noRevCompMate = false;;
	setUp.setOption(noRevCompMate, "--noRevCompMate", "Don't Reverse Complement Mate");
	setUp.pars_.ioOptions_.revComplMate_ = !noRevCompMate;
	//setUp.setOption(setUp.pars_.ioOptions_.revComplMate_, "--revCompMate", "revCompMate");
	setUp.processDirectoryOutputName(true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);

	PairedRead seq;
	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	uint64_t maxsize = 0;
	uint32_t seqCount = 0;
	while(reader.readNextRead(seq)){
		readVec::getMaxLength(seq.seqBase_, maxsize);
		readVec::getMaxLength(seq.mateSeqBase_, maxsize);
		++seqCount;
		if(seqCount > checkAmount){
			break;
		}
	}
	maxsize = maxsize * 2;
	auto alnGapPars = gapScoringParameters(
			setUp.pars_.gapInfo_.gapOpen_,
			setUp.pars_.gapInfo_.gapExtend_,
			0,0,
			0,0);
	aligner alignerObj(maxsize, alnGapPars, setUp.pars_.scoring_, false);
	alignerObj.qScorePars_.qualThresWindow_ = 0;
	reader.reOpenIn();


	uint32_t total = 0;
	auto overhangsOpts =   SeqIOOptions::genPairedOut(njh::files::make_path(setUp.pars_.directoryName_, "overhangs"));
	overhangsOpts.out_.transferOverwriteOpts(outOpts);
	SeqOutput overhangsWriter(overhangsOpts);
	uint32_t overlapFail = 0;
	uint32_t overhangFail = 0;
	uint32_t perfectOverlapCombined = 0;
	uint32_t r1EndOverR2BegCombined = 0;
	uint32_t r1BegOverR2EndCombined = 0;

	std::unique_ptr<SeqOutput> debugWriter;
	if(setUp.pars_.debug_){
		auto debugOpts =   SeqIOOptions::genFastqOut(njh::files::make_path(setUp.pars_.directoryName_, "seqs_overhangs_aligned"));
		debugWriter = std::make_unique<SeqOutput>(debugOpts);
	}

	auto isHomopolymer =
			[](const std::string & k) {
				return std::all_of(k.begin(), k.end(),[&k](const char c) {return k.front() == c;});
			};
	auto isDiNucRepeat =
			[](const std::string & k) {
				bool oddsEquals = true;
				bool evenEquals = true;
				for(const auto pos : iter::range<size_t>(3, k.size(), 2)){
					if(k[pos] != k[1]){
						oddsEquals = false;
						break;
					}
				}
				for(const auto pos : iter::range<size_t>(2, k.size(), 2)) {
					if(k[pos] != k[0]){
						evenEquals = false;
						break;
					}
				}
				return oddsEquals && evenEquals;
			};
	double complexityCutOff = 0.55;
	auto lowComplexity =
			[&complexityCutOff](const std::string & k) {
				DNABaseCounter complexityCounter;
				complexityCounter.increase(k);
				complexityCounter.setFractions();

				for(const auto & a : complexityCounter.alphabet_){
					if(complexityCounter.baseFracs_[a] >= complexityCutOff){
						return true;

					}
				}
				return false;
			};



	uint32_t totalAligned = 0;
	while(reader.readNextRead(seq)){
		++total;
		if(total % 25000 == 0 && setUp.pars_.verbose_){
			//std::cout << total << std::endl;
			std::cout << total << " " << len(seq.seqBase_) << " : " <<  len(seq.mateSeqBase_)<< std::endl;
		}

		readVec::getMaxLength(seq.seqBase_, maxsize);
		readVec::getMaxLength(seq.mateSeqBase_, maxsize);
		if(maxsize >= alignerObj.parts_.maxSize_){
			alignerObj.parts_.setMaxSize(maxsize);
		}
		motif frontMot(seq.seqBase_.seq_.substr(0, minOverlap));
		if(0 == frontMot.findPositionsFull(seq.mateSeqBase_.seq_, std::max<uint32_t>(1, round(log10(frontMot.motifOriginal_.size())))).size()){
			++overlapFail;
			continue;
		}
		/**@todo check to see if the sequence is mostly must tandem repeat */
		++totalAligned;
		alignerObj.alignRegGlobalNoInternalGaps(seq.seqBase_, seq.mateSeqBase_);
		alignerObj.profileAlignment(seq.seqBase_, seq.mateSeqBase_, false, true, true);

		if( alignerObj.comp_.distances_.eventBasedIdentityHq_ >= percentId &&
				alignerObj.comp_.distances_.basesInAln_ >= minOverlap &&
				alignerObj.comp_.hqMismatches_ + alignerObj.comp_.lqMismatches_ <= hardMismatchCutOff){

			PairedReadProcessor::AlignOverlapEnd frontCase = PairedReadProcessor::AlignOverlapEnd::UNHANDLEED;
			PairedReadProcessor::AlignOverlapEnd backCase  = PairedReadProcessor::AlignOverlapEnd::UNHANDLEED;
			if( '-' != alignerObj.alignObjectA_.seqBase_.seq_.front() &&
					'-' != alignerObj.alignObjectB_.seqBase_.seq_.front()){
				frontCase = PairedReadProcessor::AlignOverlapEnd::NOOVERHANG;
			}else if('-' != alignerObj.alignObjectA_.seqBase_.seq_.front() &&
					     '-' == alignerObj.alignObjectB_.seqBase_.seq_.front()){
				frontCase = PairedReadProcessor::AlignOverlapEnd::R1OVERHANG;
			}else if('-' == alignerObj.alignObjectA_.seqBase_.seq_.front() &&
			         '-' != alignerObj.alignObjectB_.seqBase_.seq_.front()){
				frontCase = PairedReadProcessor::AlignOverlapEnd::R2OVERHANG;
			}else{
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error not handled case for " << seq.seqBase_.name_ << "\n";
				ss << "R1.front(): " << alignerObj.alignObjectA_.seqBase_.seq_.front() << ", R2.front(): " << alignerObj.alignObjectB_.seqBase_.seq_.front();
				throw std::runtime_error{ss.str()};
			}
			if( '-' != alignerObj.alignObjectA_.seqBase_.seq_.back() &&
					'-' != alignerObj.alignObjectB_.seqBase_.seq_.back()){
				backCase = PairedReadProcessor::AlignOverlapEnd::NOOVERHANG;
			}else if('-' != alignerObj.alignObjectA_.seqBase_.seq_.back() &&
					     '-' == alignerObj.alignObjectB_.seqBase_.seq_.back()){
				backCase = PairedReadProcessor::AlignOverlapEnd::R1OVERHANG;
			}else if('-' == alignerObj.alignObjectA_.seqBase_.seq_.back() &&
			         '-' != alignerObj.alignObjectB_.seqBase_.seq_.back()){
				backCase = PairedReadProcessor::AlignOverlapEnd::R2OVERHANG;
			}else{
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error not handled case for " << seq.seqBase_.name_ << "\n";
				ss << "R1.front(): " << alignerObj.alignObjectA_.seqBase_.seq_.front() << ", R2.front(): " << alignerObj.alignObjectB_.seqBase_.seq_.front();
				throw std::runtime_error{ss.str()};
			}

			if(PairedReadProcessor::AlignOverlapEnd::NOOVERHANG == frontCase && PairedReadProcessor::AlignOverlapEnd::NOOVERHANG == backCase){
				//no over hangs, perfect overlap
				++perfectOverlapCombined;
			}else if((PairedReadProcessor::AlignOverlapEnd::NOOVERHANG == frontCase || PairedReadProcessor::AlignOverlapEnd::R1OVERHANG == frontCase) &&
					     (PairedReadProcessor::AlignOverlapEnd::NOOVERHANG == backCase  || PairedReadProcessor::AlignOverlapEnd::R2OVERHANG == backCase)){
				//ideal situation, R1 end overlaps R2 beg
				++r1EndOverR2BegCombined;
			}else if((PairedReadProcessor::AlignOverlapEnd::NOOVERHANG == frontCase || PairedReadProcessor::AlignOverlapEnd::R2OVERHANG == frontCase) &&
					     (PairedReadProcessor::AlignOverlapEnd::NOOVERHANG == backCase  || PairedReadProcessor::AlignOverlapEnd::R1OVERHANG == backCase)){
				//read through situation, R2 end overlaps R1 beg, overhang is likely illumina adaptor/primer
				uint32_t r1Start =alignerObj.alignObjectA_.seqBase_.seq_.find_first_not_of('-');
				uint32_t r2End =alignerObj.alignObjectB_.seqBase_.seq_.find_last_not_of('-') + 1;
				auto firstBase = *(alignerObj.alignObjectA_.seqBase_.seq_.begin() + r1Start);
				DNABaseCounter overlapCounter;
				overlapCounter.increaseWithRange(alignerObj.alignObjectA_.seqBase_.seq_.begin() + r1Start, alignerObj.alignObjectA_.seqBase_.seq_.begin() + r2End);
				//skip if the overlap alignment is low complexity
				if(std::all_of(alignerObj.alignObjectA_.seqBase_.seq_.begin() + r1Start,alignerObj.alignObjectA_.seqBase_.seq_.begin() + r2End,[&firstBase](const char c) {return firstBase == c;}) ||
						(1.0 - overlapCounter.calcGcContent()) <= 0.05 ||
						overlapCounter.calcGcContent() <= 0.05) {
					++overlapFail;
					continue;
				}
				//write out overhangs
				seqInfo back = alignerObj.alignObjectB_.seqBase_.getSubRead(0, r1Start);
				seqInfo front = alignerObj.alignObjectA_.seqBase_.getSubRead(r2End);
				//skip
				if(len(back) <minOverhangLen || len(front) < minOverhangLen){
					++overhangFail;
					continue;
				}
				DNABaseCounter backCounter;
				backCounter.increase(back.seq_);
				DNABaseCounter frontCounter;
				frontCounter.increase(front.seq_);
				if(backCounter.computeEntrophyBasedOffAlph(4) < entropyCutOff || frontCounter.computeEntrophyBasedOffAlph(4) < entropyCutOff){
					++overhangFail;
					continue;
				}
				back.reverseComplementRead(false, true);
				std::shared_ptr<PairedRead> overhang = std::make_shared<PairedRead>(front, back);
				overhangsWriter.openWrite(overhang);
				++r1BegOverR2EndCombined;
				if(setUp.pars_.debug_){
					debugWriter->openWrite(alignerObj.alignObjectA_);
					debugWriter->openWrite(alignerObj.alignObjectB_);
				}
			}else{
				//failure
				++overhangFail;
			}
		}else{
			++overlapFail;
		}
		if(total >= testNumber){
			break;
		}
		if(r1BegOverR2EndCombined >= gatheredNumber){
			break;
		}
	}
//	std::cout << __FILE__ << " " << __LINE__ << std::endl;
	Json::Value outVal;
	outVal["overlapFail"] = overlapFail;
	outVal["overlapFailPerc"] = (overlapFail/static_cast<double>(total)) * 100;
	outVal["overhangFail"] = overhangFail;
	outVal["overhangFailPerc"] = (overhangFail/static_cast<double>(total)) * 100;
	outVal["perfectOverlapCombined"] = perfectOverlapCombined;
	outVal["perfectOverlapCombinedPerc"] = (perfectOverlapCombined/static_cast<double>(total)) * 100;
	outVal["r1EndOverR2BegCombined"] = r1EndOverR2BegCombined;
	outVal["r1EndOverR2BegCombinedPerc"] = (r1EndOverR2BegCombined/static_cast<double>(total)) *100;
	outVal["r1BegOverR2EndCombined"] = r1BegOverR2EndCombined;
	outVal["r1BegOverR2EndCombinedPerc"] = (r1BegOverR2EndCombined/static_cast<double>(total)) *100;
	outVal["total"] = total;

	if (setUp.pars_.verbose_) {
		std::cout << outVal << std::endl;
	}
	OutputStream jsonCountsOut(
			OutOptions(
					njh::files::make_path(setUp.pars_.directoryName_, "counts.json")));
	jsonCountsOut << outVal << std::endl;
	if (overhangsWriter.outOpen()) {
		std::vector<bfs::path> rawOverHangsFnps;
		bfs::path r1Fnp = overhangsWriter.getPrimaryOutFnp();
		bfs::path r2Fnp = overhangsWriter.getSecondaryOutFnp();
		if (bfs::exists(r1Fnp)) {
			//		uint32_t windowSize = 5;
			//		uint32_t thresholdStreakNumber = 3;
			//		uint32_t threshold = 20;
			std::string polyA = "AAAA";
			std::string polyG = "GGGG";
			std::string polyT = "TTTT";
			std::string polyC = "CCCC";
			std::string polyN = "NNNN";
			auto trimAtPolys =
					[&polyA,&polyG,&polyT,&polyC,&polyN](seqInfo & seqToTrim) {
						auto polyAPos = seqToTrim.seq_.find(polyA);
						if(std::string::npos != polyAPos) {
							readVecTrimmer::trimToMaxLength(seqToTrim, polyAPos);
						}
						auto polyGPos = seqToTrim.seq_.find(polyG);
						if (std::string::npos != polyGPos) {
							readVecTrimmer::trimToMaxLength(seqToTrim, polyGPos);
						}
						auto polyTPos = seqToTrim.seq_.find(polyT);
						if(std::string::npos != polyTPos) {
							readVecTrimmer::trimToMaxLength(seqToTrim, polyTPos);
						}
						auto polyCPos = seqToTrim.seq_.find(polyC);
						if (std::string::npos != polyCPos) {
							readVecTrimmer::trimToMaxLength(seqToTrim, polyCPos);
						}
						auto polyNPos = seqToTrim.seq_.find(polyN);
						if (std::string::npos != polyNPos) {
							readVecTrimmer::trimToMaxLength(seqToTrim, polyNPos);
						}
						if(seqToTrim.seq_.size() >= 58) {
							readVecTrimmer::trimAtRstripBase(seqToTrim, 'A');
						}
					};

			overhangsWriter.closeOut();
			{
				//r1 process
				SeqInput r1Reader(SeqIOOptions::genFastqIn(r1Fnp));
				r1Reader.openIn();
				seqInfo seq;
				std::unordered_map<std::string, uint32_t> frontCounts;
				uint32_t total = 0;
				while (r1Reader.readNextRead(seq)) {
					if (len(seq) >= frontCheck) {
						++frontCounts[seq.seq_.substr(0, frontCheck)];

					}
				}
				std::unordered_map<std::string, double> frontFracs;
				for (const auto & count : frontCounts) {
					//get rid of low complexity
					if (!lowComplexity(count.first)
							&& !isHomopolymer(count.first.substr(0, frontCheck / 2))
							&& !isHomopolymer(count.first.substr(frontCheck / 2))
							&& !isDiNucRepeat(count.first)
							&& std::string::npos == count.first.find("N")) {
						total += count.second;
					}
				}
				for (const auto & count : frontCounts) {
					//get rid of low complexity
					if (!lowComplexity(count.first)
							&& !isHomopolymer(count.first.substr(0, frontCheck / 2))
							&& !isHomopolymer(count.first.substr(frontCheck / 2))
							&& !isDiNucRepeat(count.first)
							&& std::string::npos == count.first.find("N")) {
						frontFracs[count.first] = count.second / static_cast<double>(total);
					}
				}
				auto frontFracsKeys = njh::getVecOfMapKeys(frontFracs);
				njh::sort(frontFracsKeys,
						[&frontFracs](const std::string & key1, const std::string & key2) {
							return frontFracs[key1] > frontFracs[key2];
						});
				OutOptions top5FrontsR1Opts(
						njh::files::make_path(setUp.pars_.directoryName_,
								"topFronts_r1.txt"));
				OutputStream top5FrontsR1Out(top5FrontsR1Opts);
				top5FrontsR1Out << "seq\tcount\tfrac" << std::endl;
				uint32_t count = 0;
				VecStr aboveCutOffKeys;
				for (const auto & key : frontFracsKeys) {
					top5FrontsR1Out << key << "\t" << frontCounts[key] << "\t"
							<< frontFracs[key] << std::endl;
					bool oneOff = false;
					for (const auto & otherKey : aboveCutOffKeys) {
						if (numberOfMismatches(otherKey, key) <= 1) {
							oneOff = true;
						}
					}
					if (!oneOff && frontFracs[key] >= fracCutOff
							&& frontCounts[key] >= consensusCountCutOff) {
						aboveCutOffKeys.emplace_back(key);
					}
					++count;
					if (count >= 5) {
						break;
					}
				}
				if (aboveCutOffKeys.empty() && !frontFracsKeys.empty()
						&& frontCounts[frontFracsKeys.front()] >= consensusCountCutOff) {
					aboveCutOffKeys.emplace_back(frontFracsKeys.front());
				}
				if (!aboveCutOffKeys.empty()) {
					MultiSeqIO writers;
					for (const auto & aboveCutOffKey : aboveCutOffKeys) {
						writers.addReader(aboveCutOffKey,
								SeqIOOptions::genFastqOut(
										njh::files::make_path(setUp.pars_.directoryName_,
												"r1_" + aboveCutOffKey)));
						writers.addReader(aboveCutOffKey + "-raw",
								SeqIOOptions::genFastqOut(
										njh::files::make_path(setUp.pars_.directoryName_, "raw_r1_" + aboveCutOffKey))) ;
						rawOverHangsFnps.emplace_back(njh::files::make_path(setUp.pars_.directoryName_, "raw_r1_" + aboveCutOffKey + ".fastq"));
					}
					r1Reader.reOpenIn();
					std::unordered_map<std::string, std::map<uint32_t, DNABaseCounter>> consensusCounts;
					while (r1Reader.readNextRead(seq)) {
						if (len(seq.seq_) < frontCheck) {
							continue;
						}
						for (const auto & aboveCutOffKey : aboveCutOffKeys) {
							if (numberOfMismatches(seq.seq_.substr(0, frontCheck),
									aboveCutOffKey) <= 1) {
								//							uint32_t streakNumber = 0;
								//							for (auto pos : iter::range<uint32_t>(0, len(seq) - windowSize + 1, 1)) {
								//								auto qualMean = std::accumulate(seq.qual_.begin() + pos, seq.qual_.begin() + pos + windowSize, 0.0)/windowSize;
								//								if(qualMean <=threshold){
								//									++streakNumber;
								//								}else{
								//									streakNumber = 0;
								//								}
								//								if(streakNumber == thresholdStreakNumber){
								//									readVecTrimmer::trimToMaxLength(seq, pos);
								//									break;
								//								}
								//							}
								if(writeOutBaseQualCompositionProfile){
									writers.openWrite(aboveCutOffKey + "-raw", seq);
								}
								trimAtPolys(seq);

								//readVecTrimmer::trimAtRstripBase(seq, 'A');
								if (len(seq) >= frontCheck) {
									for (const auto pos : iter::range(len(seq))) {
										if (seq.qual_[pos] > qualCutOff) {
											consensusCounts[aboveCutOffKey][pos].increase(
													seq.seq_[pos]);
										}
									}
									if (setUp.pars_.debug_) {
										writers.openWrite(aboveCutOffKey, seq);
									}
								}
								break;
							}
						}
					}
					SeqOutput conWriter(
							SeqIOOptions::genFastaOut(
									njh::files::make_path(setUp.pars_.directoryName_,
											"r1_consensus")));
					std::vector<seqInfo> cons;
					for (const auto & conCounts : consensusCounts) {
						std::string conSeq;
						for (const auto & count : conCounts.second) {
							if (count.second.getTotalCount() >= consensusCountCutOff) {
								uint32_t bestCount = 0;
								std::vector<char> bestBases;
								for (const auto & b : count.second.alphabet_) {
									if (count.second.bases_[b] > bestCount) {
										bestBases.clear();
										bestBases.push_back(b);
										bestCount = count.second.bases_[b];
									} else if (count.second.bases_[b] == bestCount) {
										bestBases.push_back(b);
									}
								}
								if (bestBases.size() == 1) {
									conSeq.push_back(bestBases.front());
								} else if (bestBases.size() > 1) {
									//multiple bases tied
									conSeq.push_back('N');
								}
							} else {
								break;
							}
						}
						seqInfo conSeqInfo(conCounts.first, conSeq);
						readVecTrimmer::trimToMaxLength(conSeqInfo, outPrimerLengthMax);
						trimAtPolys(conSeqInfo);
						readVecTrimmer::trimAtRstripBase(conSeqInfo, 'N');
						if ("" != conSeqInfo.seq_ && len(conSeqInfo) > outPrimerLengthMin) {
							bool alreadyAdded = false;
							for (const auto & otherCon : cons) {
								if (otherCon.seq_ == conSeqInfo.seq_) {
									alreadyAdded = true;
									break;
								}
							}
							if (!alreadyAdded) {
								cons.emplace_back(conSeqInfo);
							}
						}
					}
					if (!cons.empty()) {
						conWriter.openWrite(cons);
					}
				}
			}

			{
				//r2 process
				SeqInput r2Reader(SeqIOOptions::genFastqIn(r2Fnp));
				r2Reader.openIn();
				seqInfo seq;
				std::unordered_map<std::string, uint32_t> frontCounts;
				uint32_t total = 0;
				while (r2Reader.readNextRead(seq)) {
					if (len(seq) >= frontCheck) {
						++frontCounts[seq.seq_.substr(0, frontCheck)];

					}
				}
				std::unordered_map<std::string, double> frontFracs;
				for (const auto & count : frontCounts) {
					//get rid of low complexity
					if (!lowComplexity(count.first)
							&& !isHomopolymer(count.first.substr(0, frontCheck / 2))
							&& !isHomopolymer(count.first.substr(frontCheck / 2))
							&& !isDiNucRepeat(count.first)
							&& std::string::npos == count.first.find("N")) {
						total += count.second;
					}
				}
				for (const auto & count : frontCounts) {
					//get rid of low complexity
					if (!lowComplexity(count.first)
							&& !isHomopolymer(count.first.substr(0, frontCheck / 2))
							&& !isHomopolymer(count.first.substr(frontCheck / 2))
							&& !isDiNucRepeat(count.first)
							&& std::string::npos == count.first.find("N")) {
						frontFracs[count.first] = count.second / static_cast<double>(total);
					}
				}
				auto frontFracsKeys = njh::getVecOfMapKeys(frontFracs);
				njh::sort(frontFracsKeys,
						[&frontFracs](const std::string & key1, const std::string & key2) {
							return frontFracs[key1] > frontFracs[key2];
						});
				OutOptions top5Frontsr2Opts(
						njh::files::make_path(setUp.pars_.directoryName_,
								"topFronts_r2.txt"));
				OutputStream top5Frontsr2Out(top5Frontsr2Opts);
				top5Frontsr2Out << "seq\tcount\tfrac" << std::endl;
				uint32_t count = 0;
				VecStr aboveCutOffKeys;
				for (const auto & key : frontFracsKeys) {
					top5Frontsr2Out << key << "\t" << frontCounts[key] << "\t"
							<< frontFracs[key] << std::endl;
					bool oneOff = false;
					for (const auto & otherKey : aboveCutOffKeys) {
						if (numberOfMismatches(otherKey, key) <= 1) {
							oneOff = true;
						}
					}
					if (!oneOff && frontFracs[key] >= fracCutOff
							&& frontCounts[key] >= consensusCountCutOff) {
						aboveCutOffKeys.emplace_back(key);
					}
					++count;
					if (count >= 5) {
						break;
					}
				}
				if (aboveCutOffKeys.empty() && !frontFracsKeys.empty()
						&& frontCounts[frontFracsKeys.front()] >= consensusCountCutOff) {
					aboveCutOffKeys.emplace_back(frontFracsKeys.front());
				}
				if (!aboveCutOffKeys.empty()) {
					MultiSeqIO writers;
					for (const auto & aboveCutOffKey : aboveCutOffKeys) {
						writers.addReader(aboveCutOffKey,
								SeqIOOptions::genFastqOut(
										njh::files::make_path(setUp.pars_.directoryName_,
												"r2_" + aboveCutOffKey)));
						writers.addReader(aboveCutOffKey + "-raw",
								SeqIOOptions::genFastqOut(
										njh::files::make_path(setUp.pars_.directoryName_, "raw_r2_" + aboveCutOffKey)));
						rawOverHangsFnps.emplace_back(njh::files::make_path(setUp.pars_.directoryName_, "raw_r2_" + aboveCutOffKey + ".fastq"));

					}
					r2Reader.reOpenIn();
					std::unordered_map<std::string, std::map<uint32_t, DNABaseCounter>> consensusCounts;
					while (r2Reader.readNextRead(seq)) {
						if (len(seq.seq_) < frontCheck) {
							continue;
						}
						for (const auto & aboveCutOffKey : aboveCutOffKeys) {
							if (numberOfMismatches(seq.seq_.substr(0, frontCheck),
									aboveCutOffKey) <= 1) {
								//							uint32_t streakNumber = 0;
								//							for (auto pos : iter::range<uint32_t>(0, len(seq) - windowSize + 1, 1)) {
								//								auto qualMean = std::accumulate(seq.qual_.begin() + pos, seq.qual_.begin() + pos + windowSize, 0.0)/windowSize;
								//								if(qualMean <=threshold){
								//									++streakNumber;
								//								}else{
								//									streakNumber = 0;
								//								}
								//								if(streakNumber == thresholdStreakNumber){
								//									readVecTrimmer::trimToMaxLength(seq, pos);
								//									break;
								//								}
								//							}
								if(writeOutBaseQualCompositionProfile){
									writers.openWrite(aboveCutOffKey + "-raw", seq);
								}
								trimAtPolys(seq);

								//readVecTrimmer::trimAtRstripBase(seq, 'A');
								if (len(seq) >= frontCheck) {
									for (const auto pos : iter::range(len(seq))) {
										if (seq.qual_[pos] > qualCutOff) {
											consensusCounts[aboveCutOffKey][pos].increase(
													seq.seq_[pos]);
										}
									}
									if (setUp.pars_.debug_) {
										writers.openWrite(aboveCutOffKey, seq);
									}
								}
								break;
							}
						}
					}
					SeqOutput conWriter(
							SeqIOOptions::genFastaOut(
									njh::files::make_path(setUp.pars_.directoryName_,
											"r2_consensus")));
					std::vector<seqInfo> cons;
					for (const auto & conCounts : consensusCounts) {
						std::string conSeq;
						for (const auto & count : conCounts.second) {
							if (count.second.getTotalCount() >= consensusCountCutOff) {
								uint32_t bestCount = 0;
								std::vector<char> bestBases;
								for (const auto & b : count.second.alphabet_) {
									if (count.second.bases_[b] > bestCount) {
										bestBases.clear();
										bestBases.push_back(b);
										bestCount = count.second.bases_[b];
									} else if (count.second.bases_[b] == bestCount) {
										bestBases.push_back(b);
									}
								}
								if (bestBases.size() == 1) {
									conSeq.push_back(bestBases.front());
								} else if (bestBases.size() > 1) {
									//multiple bases tied
									conSeq.push_back('N');
								}
							} else {
								break;
							}
						}
						seqInfo conSeqInfo(conCounts.first, conSeq);
						readVecTrimmer::trimToMaxLength(conSeqInfo, outPrimerLengthMax);
						trimAtPolys(conSeqInfo);
						readVecTrimmer::trimAtRstripBase(conSeqInfo, 'N');
						if ("" != conSeqInfo.seq_ && len(conSeqInfo) > outPrimerLengthMin) {
							bool alreadyAdded = false;
							for (const auto & otherCon : cons) {
								if (otherCon.seq_ == conSeqInfo.seq_) {
									alreadyAdded = true;
									break;
								}
							}
							if (!alreadyAdded) {
								cons.emplace_back(conSeqInfo);
							}
						}
					}
					if (!cons.empty()) {
						conWriter.openWrite(cons);
					}
				}
			}
			if (!keepOverHangs) {
				bfs::remove(r1Fnp);
				bfs::remove(r2Fnp);
			}
		}
		if(writeOutBaseQualCompositionProfile){
			for(const auto & fnp : rawOverHangsFnps){
				if(bfs::exists(fnp)){
					std::unordered_map<uint32_t, std::unordered_map<char, std::unordered_map<uint32_t, uint32_t>>> counts;
					auto fnpOpts = SeqIOOptions::genFastqIn(fnp, false);
					seqInfo seq;
					SeqInput reader(fnpOpts);
					reader.openIn();
					while(reader.readNextRead(seq)){
						for(uint32_t pos = 0; pos < len(seq); ++pos){
							++counts[pos][seq.seq_[pos]][seq.qual_[pos]];
						}
					}
					OutOptions qualCountsFnpOpt(bfs::path(fnp).replace_extension(".tab.txt"));
					OutputStream qualCountsFnpOut(qualCountsFnpOpt);
					qualCountsFnpOut << "pos\tchar\tqual\tcount" << std::endl;
					auto posKeys = getVectorOfMapKeys(counts);
					njh::sort(posKeys);
					for(const auto & posKey : posKeys){
						auto baseKeys = getVectorOfMapKeys(counts[posKey]);
						njh::sort(baseKeys);
						for(const auto & baseKey : baseKeys){
							auto qualKeys = getVectorOfMapKeys(counts[posKey][baseKey]);
							njh::sort(qualKeys);
							for(const auto & qualKey : qualKeys){
								qualCountsFnpOut << posKey
										<< "\t" << baseKey
										<< "\t" << qualKey
										<< "\t" << counts[posKey][baseKey][qualKey] << std::endl;
							}
						}
					}
				}
			}
		}
	} else {
		std::cerr << "No overhangs" << std::endl;
	}

	return 0;
}



}  // namespace njhseq



