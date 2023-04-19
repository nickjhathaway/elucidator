/*
 * bedExpRunner.cpp
 *
 *  Created on: May 18, 2015
 *      Author: nickhathaway
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


#include "bedExp.hpp"
#include <njhseq/objects/BioDataObject.h>
#include "elucidator/BioRecordsUtils/BedUtility.hpp"
#include <njhseq/objects/counters/DNABaseCounter.hpp>
#include <njhseq/programUtils/seqSetUp.hpp>
#include <njhseq/concurrency/pools.h>
#include <njhseq/BamToolsUtils.h>


#include <TwoBit.h>


namespace njhseq {
bedExpRunner::bedExpRunner()
    : njh::progutils::ProgramRunner(
          {
					 addFunc("getFastaWithBed", getFastaWithBed, false),
					 addFunc("getSeqFromTwoBit", getSeqFromTwoBit, false),
					 addFunc("mergeOverlappingRegionsSameName", mergeOverlappingRegionsSameName, false),


					 addFunc("splitBedFile", splitBedFile, false),
					 addFunc("separateOutRecordsInBedFile", separateOutRecordsInBedFile, false),
					 addFunc("bedUnqiue", bedUnqiue, false),
					 addFunc("getCloseBedRegions", getCloseBedRegions, false),
					 addFunc("extendBedRegions", extendBedRegions, false),
					 addFunc("bed3ToBed6", bed3ToBed6, false),
					 addFunc("filterBedRecordsByLength", filterBedRecordsByLength, false),
					 addFunc("getOverlappingBedRegions", getOverlappingBedRegions, false),
					 addFunc("getNonOverlappingBedRegions", getNonOverlappingBedRegions, false),
					 addFunc("bedRenameRepeatUids", bedRenameRepeatUids, false),
					 addFunc("bedRenameWithCoords", bedRenameWithCoords, false),
					 addFunc("getUpstreamRegion", getUpstreamRegion, false),
					 addFunc("getDownstreamRegion", getDownstreamRegion, false),
					 addFunc("extendUpstreamRegion", extendUpstreamRegion, false),
					 addFunc("extendDownstreamRegion", extendDownstreamRegion, false),
					 addFunc("bedToggleStrand", bedToggleStrand, false),
					 addFunc("bedCreateSpanningRegions", bedCreateSpanningRegions, false),
					 addFunc("getBestNonOverlapingRegions", getBestNonOverlapingRegions, false),
					 addFunc("bedAddSmartIDForPlotting", bedAddSmartIDForPlotting, false),
					 addFunc("extendToEndOfChrom", extendToEndOfChrom, false),
					 addFunc("extendToStartOfChrom", extendToStartOfChrom, false),
					 addFunc("bedRenameChromosomes", bedRenameChromosomes, false),
					 addFunc("bedChangeScoreToLength", bedChangeScoreToLength, false),
					 addFunc("bedCoordSort", bedCoordSort, false),
					 addFunc("extractBedRecordsWithName", extractBedRecordsWithName, false),
					 addFunc("reorientBasedOnSingleReadsOrientationCounts", reorientBasedOnSingleReadsOrientationCounts, false),
					 addFunc("getFirstRegionPerChrom", getFirstRegionPerChrom, false),
					 addFunc("getLastRegionPerChrom", getLastRegionPerChrom, false),
					 addFunc("getGCContentOfRegions", getGCContentOfRegion, false),
					 addFunc("getEntropyOfRegion", getEntropyOfRegion, false),
					 addFunc("getNumOfNsOfRegion", getNumOfNsOfRegion, false),
					 addFunc("getNumOfNsEntropyGCContentOfRegion", getNumOfNsEntropyGCContentOfRegion, false),
					 addFunc("bedKeepRegionsCompletelyInOther", bedKeepRegionsCompletelyInOther, false),
					 addFunc("bedGetRegionsEncompassingOthers", bedGetRegionsEncompassingOthers, false),
					 addFunc("bedGetOverlappinBetweenPanels", bedGetOverlappinBetweenPanels, false),
					 addFunc("bedGetOverlappinCoverageBetweenPanels", bedGetOverlappinCoverageBetweenPanels, false),
					 addFunc("bedGetOverlappinPositionsBetweenPanels", bedGetOverlappinPositionsBetweenPanels, false),
//

					 addFunc("getLongestHomopolymerLengthInRegion", getLongestHomopolymerLengthInRegion, false),
					 addFunc("roughSmoothingForBedCoverage", roughSmoothingForBedCoverage, false),
					 addFunc("fastaToBed", fastaToBed, false),
					 addFunc("getInterveningRegions", getInterveningRegions, false),
					 addFunc("fillInRegionsByBestScore", fillInRegionsByBestScore, false),
					 addFunc("centerBedRegionWithFixSize", centerBedRegionWithFixSize, false),

					 addFunc("trimBedRegions", trimBedRegions, false),
					 addFunc("trimUpstreamRegion", trimUpstreamRegion, false),
					 addFunc("trimDownstreamRegion", trimDownstreamRegion, false),
					 addFunc("createWindowsInbetweenRegions", createWindowsInbetweenRegions, false),
					 addFunc("createWindowsInRegions", createWindowsInRegions, false),
					 addFunc("reverseComplementRegion", reverseComplementRegion, false),
					 addFunc("bedRemoveOveringLappingRegions", bedRemoveOveringLappingRegions, false),
					 addFunc("bedFilterRegionsCompletelyInOther", bedFilterRegionsCompletelyInOther, false),
					 addFunc("differentSubRegionCombos", differentSubRegionCombos, false),
					 addFunc("bedFilterRegionsStartingOrEndingInOther", bedFilterRegionsStartingOrEndingInOther, false),
					 addFunc("bedKeepRegionsStartingOrEndingInOther", bedKeepRegionsStartingOrEndingInOther, false),
					 addFunc("getAverageDistanceToOtherRegions", getAverageDistanceToOtherRegions, false),
					 addFunc("getDistanceToClostestRegion", getDistanceToClostestRegion, false),
					 addFunc("getIntersectionBetweenTwoBedFiles", getIntersectionBetweenTwoBedFiles, false),
					 addFunc("removeSubRegionsFromBedFile", removeSubRegionsFromBedFile, false),
					 addFunc("getDegreeOfOverlappingBedRegions", getDegreeOfOverlappingBedRegions, false),
					 addFunc("bedGetOverlappingRegionsFromPanels", bedGetOverlappingRegionsFromPanels, false),
					 addFunc("bedAddRegionsNonCompletelyInOther", bedAddRegionsNonCompletelyInOther, false),
					 addFunc("bedBinCloseRegions", bedBinCloseRegions, false),
           addFunc("createBedRegionFromName", createBedRegionFromName, false),
           },//
          "bedExp") {}





int bedExpRunner::bedRemoveOveringLappingRegions(const njh::progutils::CmdArgs & inputCommands) {

	bfs::path bedFile;
	bfs::path intersectWithBed;
	OutOptions outOpts;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.setOption(bedFile, "--bed", "Bed file to parse", true);
	setUp.setOption(intersectWithBed, "--intersectWithBed", "Bed file to intersect with", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	OutputStream out(outOpts);
	auto rawIntersectingBeds = getBeds(intersectWithBed);
	auto regions = getBeds(bedFile);
	BedUtility::coordSort(regions, false);
	BedUtility::coordSort(rawIntersectingBeds, false);
	std::vector<std::shared_ptr<Bed6RecordCore>> intersectingBeds;
	intersectingBeds.emplace_back(
			std::make_shared<Bed6RecordCore>(*rawIntersectingBeds.front()));

	for (const auto regPos : iter::range<uint32_t>(1, rawIntersectingBeds.size())) {
		if (intersectingBeds.back()->overlaps(*rawIntersectingBeds[regPos], 1) || (intersectingBeds.back()->chrom_ == rawIntersectingBeds[regPos]->chrom_ && intersectingBeds.back()->chromEnd_ == rawIntersectingBeds[regPos]->chromStart_) ) {
			intersectingBeds.back()->chromEnd_ = std::max(
					intersectingBeds.back()->chromEnd_,
					rawIntersectingBeds[regPos]->chromEnd_);
		} else {
			intersectingBeds.emplace_back(
					std::make_shared<Bed6RecordCore>(*rawIntersectingBeds[regPos]));
		}
	}
	if(setUp.pars_.verbose_){
		std::cout << "intersectingBeds.size(): " << intersectingBeds.size() << std::endl;
	}


	for(const auto & reg : regions){
		std::vector<std::shared_ptr<Bed6RecordCore>> overlappingRegions;
		for(const auto & inter : intersectingBeds){
			if(inter->chrom_ < reg->chrom_){
				//bother regions are sorted if we haven't reached this region's chromosome yet
				continue;
			}
			if(inter->chrom_ > reg->chrom_){
				//both regions are sorted so if we run into this we can break
				break;
			}
			if(inter->chrom_ == reg->chrom_ && inter->chromStart_ >= reg->chromEnd_){
				//both regions are sorted so if we run into this we can break
				break;
			}
			if(inter->chrom_ == reg->chrom_ && inter->chromEnd_ < reg->chromStart_){
				continue;
			}
			if(reg->overlaps(*inter, 1)){
				overlappingRegions.emplace_back(inter);
			}
		}
		if(overlappingRegions.empty()){
			out << reg->toDelimStr() << "\n";
		}else{
			if(overlappingRegions.front()->chromStart_ > reg->chromStart_){
				Bed6RecordCore outReg(reg->chrom_, reg->chromStart_, overlappingRegions.front()->chromStart_, "",0, reg->strand_);
				outReg.score_ = outReg.length();
				outReg.name_ = njh::pasteAsStr(outReg.chrom_, "-", outReg.chromStart_, "-", outReg.chromEnd_, "-", outReg.strand_ == '+' ? "for" : "rev");
				out << outReg.toDelimStrWithExtra() << "\n";
			}
			if(overlappingRegions.back()->chromEnd_ < reg->chromEnd_){
				Bed6RecordCore outReg(reg->chrom_, overlappingRegions.back()->chromEnd_, reg->chromEnd_, "",0, reg->strand_);
				outReg.score_ = outReg.length();
				outReg.name_ = njh::pasteAsStr(outReg.chrom_, "-", outReg.chromStart_, "-", outReg.chromEnd_, "-", outReg.strand_ == '+' ? "for" : "rev");
				out << outReg.toDelimStrWithExtra() << "\n";
			}
			if(overlappingRegions.size() > 1){
				for(uint32_t pos = 0; pos < overlappingRegions.size() -1; ++ pos){
					if(overlappingRegions[pos]->chromEnd_ != overlappingRegions[pos + 1]->chromStart_){
						Bed6RecordCore outReg(reg->chrom_, overlappingRegions[pos]->chromEnd_, overlappingRegions[pos + 1]->chromStart_, "",0, reg->strand_);
						outReg.score_ = outReg.length();
						outReg.name_ = njh::pasteAsStr(outReg.chrom_, "-", outReg.chromStart_, "-", outReg.chromEnd_, "-", outReg.strand_ == '+' ? "for" : "rev");
						out << outReg.toDelimStrWithExtra() << "\n";
					}
				}
			}
		}
	}
	return 0;
}

int bedExpRunner::reverseComplementRegion(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path chromLengthsTable = "";

	bfs::path bedFile;
	OutOptions outOpts;
	outOpts.outExtention_ = ".bed";
	seqSetUp setUp(inputCommands);
	setUp.setOption(bedFile, "--bed", "Bed file", true);
	setUp.setOption(chromLengthsTable, "--chromLengthsTable", "Chromosome Lengths Table, two columns, first chrom name, second length", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	table chromTab(chromLengthsTable, "\t", false);

	if(2 != chromTab.columnNames_.size()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << chromLengthsTable << " should be a table with two columns, first chrom name, second length" << "\n";
		throw std::runtime_error{ss.str()};
	}
	std::unordered_map<std::string, uint64_t> lengths;
	for(const auto & row : chromTab.content_){
		lengths[row[0]] = njh::StrToNumConverter::stoToNum<uint64_t>(row[1]);
	}
	BioDataFileIO<Bed3RecordCore> reader{IoOptions(InOptions(bedFile), outOpts)};
	reader.openIn();
	reader.openOut();
	std::unordered_map<std::string, std::unordered_map<uint32_t, std::vector<uint32_t>>> alreadyTakenIds;
	std::shared_ptr<Bed3RecordCore> b = reader.readNextRecord();
	while(nullptr != b){
		const auto & chromLen = njh::mapAt(lengths, b->chrom_);
		if(njh::mapAt(lengths, b->chrom_) < b->chromStart_ || njh::mapAt(lengths, b->chrom_) < b->chromEnd_){
			std::stringstream ss;
			std::string name = njh::pasteAsStr(b->chrom_, "-", b->chromStart_, "-", b->chromEnd_);
			if(!b->extraFields_.empty()){
				name = b->extraFields_[0];
			}
			ss << __PRETTY_FUNCTION__ << ", error for region: " <<  name << "start: " << b->chromStart_ << " or end: " << b->chromEnd_ << " is greater than given chromosome length: " << chromLen<< "\n";
			throw std::runtime_error{ss.str()};
		}
		auto newEnd = chromLen - b->chromStart_;
		auto newStart = chromLen - b->chromEnd_;
		b->chromStart_ = newStart;
		b->chromEnd_ = newEnd;
		if(b->extraFields_.size() >= 3){
			auto possibleStrandField = b->extraFields_[2];
			if("-" == possibleStrandField){
				possibleStrandField = "+";
			}else if("+" == possibleStrandField){
				possibleStrandField = "-";
			}
			b->extraFields_[2] = possibleStrandField;
		}
		reader.write(*b, [](const Bed3RecordCore & bed, std::ostream & out){
			out << bed.toDelimStrWithExtra() << std::endl;
		});
		b = reader.readNextRecord();
	}

	return 0;
}


int bedExpRunner::centerBedRegionWithFixSize(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bedFile;
	bool skipOutOfRangeRegions = false;
	OutOptions outOpts(bfs::path(""), ".bed");
	uint32_t windowSize = 25;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.setOption(skipOutOfRangeRegions, "--skipOutOfRangeRegions", "skip Out Of Range Regions (e.g. start=5 with windowsize = 100");

	setUp.setOption(windowSize, "--windowSize", "Window Size", true);
	setUp.setOption(bedFile, "--bed", "Bed file to parse", true);

	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);


	BioDataFileIO<Bed3RecordCore> reader{IoOptions(InOptions(bedFile))};
	reader.openIn();

	std::ofstream outFile;
	std::ostream out(determineOutBuf(outFile, outOpts));
	Bed3RecordCore reg;
	while(reader.readNextRecord(reg)){
		if(reg.length() > windowSize - 2){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << reg.toDelimStr() << " size is longer than desired window size"<< "\n";
			if(skipOutOfRangeRegions){
				if(setUp.pars_.verbose_){
					std::cerr << ss.str();
				}
				out << reg.toDelimStrWithExtra() << std::endl;
			}else{
				throw std::runtime_error{ss.str()};
			}
		}else{
			uint32_t left = (windowSize - reg.length())/2;
			uint32_t right = windowSize - reg.length() - left;
			if(left > reg.chromStart_){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "can't extend region " <<reg.toDelimStr() << " to left by " << left << "\n";
				if(skipOutOfRangeRegions){
					if(setUp.pars_.verbose_){
						std::cerr << ss.str();
					}
					continue;
				}else{
					throw std::runtime_error{ss.str()};
				}
			}

			auto oldLength = reg.length();
			reg.chromStart_ = reg.chromStart_ - left;
			reg.chromEnd_ += right;
			//reg.score_ = reg.length();
			if(reg.extraFields_.size() >=3){
				bool scoreIsLen = estd::to_string(oldLength) == reg.extraFields_[1];
				if(scoreIsLen){
					reg.extraFields_[1] = estd::to_string(reg.length());
				}
			}
			out << reg.toDelimStrWithExtra() << std::endl;
		}
	}
	return 0;
}

int bedExpRunner::fastaToBed(const njh::progutils::CmdArgs & inputCommands) {
	OutOptions outOpts(bfs::path(""), ".bed");
	seqSetUp setUp(inputCommands);
	setUp.processReadInNames(true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	seqInfo seq;
	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	OutputStream out(outOpts);
	while(reader.readNextRead(seq)){
		out << seq.name_
				<< "\t" << 0
				<< "\t" << len(seq) << std::endl;
	}



	return 0;
}


int bedExpRunner::getFirstRegionPerChrom(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bedFile = "";
	OutOptions outOpts;
	outOpts.outExtention_ = ".bed";
	seqSetUp setUp(inputCommands);
	setUp.setOption(bedFile, "--bed", "Bed file", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	BioDataFileIO<Bed3RecordCore> reader{IoOptions(InOptions(bedFile), outOpts)};
	reader.openIn();
	reader.openOut();
	Bed3RecordCore record;
	std::unordered_map<std::string, Bed3RecordCore> records;
	while(reader.readNextRecord(record)){
		if(!njh::in(record.chrom_, records)){
			records[record.chrom_] = record;
		}else{
			if(record.chromStart_ < records[record.chrom_].chromStart_ ){
				records[record.chrom_] = record;
			}else if(record.chromStart_ == records[record.chrom_].chromStart_ && record.chromEnd_ > records[record.chrom_].chromEnd_){
				records[record.chrom_] = record;
			}
		}
	}

	auto chromNames = njh::getVecOfMapKeys(records);
	njh::sort(chromNames);
	for(const auto & chromName : chromNames){
		reader.write(records[chromName], [](const Bed3RecordCore & reg, std::ostream & out){
			out << reg.toDelimStrWithExtra() << std::endl;
		});
	}

	return 0;
}

int bedExpRunner::getLastRegionPerChrom(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bedFile = "";
	OutOptions outOpts;
	outOpts.outExtention_ = ".bed";
	seqSetUp setUp(inputCommands);
	setUp.setOption(bedFile, "--bed", "Bed file", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	BioDataFileIO<Bed3RecordCore> reader{IoOptions(InOptions(bedFile), outOpts)};
	reader.openIn();
	reader.openOut();
	Bed3RecordCore record;
	std::unordered_map<std::string, Bed3RecordCore> records;
	while(reader.readNextRecord(record)){
		if(!njh::in(record.chrom_, records)){
			records[record.chrom_] = record;
		}else{
			if(record.chromStart_ > records[record.chrom_].chromStart_ ){
				records[record.chrom_] = record;
			}else if(record.chromStart_ == records[record.chrom_].chromStart_ && record.chromEnd_ > records[record.chrom_].chromEnd_){
				records[record.chrom_] = record;
			}
		}
	}

	auto chromNames = njh::getVecOfMapKeys(records);
	njh::sort(chromNames);
	for(const auto & chromName : chromNames){
		reader.write(records[chromName], [](const Bed3RecordCore & reg, std::ostream & out){
			out << reg.toDelimStrWithExtra() << std::endl;
		});
	}


	return 0;
}



int bedExpRunner::bedCoordSort(const njh::progutils::CmdArgs & inputCommands) {
	bool decending = false;
	bool useStrandInfo = false;
	bfs::path bedFile = "";
	OutOptions outOpts;
	outOpts.outExtention_ = ".bed";
	seqSetUp setUp(inputCommands);
	setUp.setOption(useStrandInfo, "--useStrandInfo", "Use Strand Info");
	setUp.setOption(decending, "--decending", "decending");
	setUp.setOption(bedFile, "--bed", "Bed file", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	if(useStrandInfo){
		BioDataFileIO<Bed6RecordCore> reader{IoOptions(InOptions(bedFile), outOpts)};
		reader.openIn();
		reader.openOut();
		std::vector<Bed6RecordCore> allBeds;
		Bed6RecordCore record;
		while(reader.readNextRecord(record)){
			allBeds.emplace_back(record);
		}

		std::function<bool(const Bed6RecordCore &, const Bed6RecordCore &)> bedCoordSorterFunc =
				[](const Bed6RecordCore & reg1, const Bed6RecordCore & reg2) {
					if(reg1.chrom_ == reg2.chrom_) {
						uint32_t reg1Start = ('+' == reg1.strand_ ? reg1.chromStart_ : reg1.chromEnd_);
						uint32_t reg2Start = ('+' == reg2.strand_ ? reg2.chromStart_ : reg2.chromEnd_);

						uint32_t reg1End = ('+' == reg1.strand_ ? reg1.chromEnd_ : reg1.chromStart_);
						uint32_t reg2End = ('+' == reg2.strand_ ? reg2.chromEnd_ : reg2.chromStart_);

						if(reg1Start == reg2Start) {
							return reg1End < reg2End;
						} else {
							return reg1Start < reg2Start;
						}
					} else {
						return reg1.chrom_ < reg2.chrom_;
					}
				};
		if(decending){
			std::sort(allBeds.rbegin(), allBeds.rend(), bedCoordSorterFunc);
		}else{
			njh::sort(allBeds, bedCoordSorterFunc);
		}

		reader.write(allBeds, [](const Bed6RecordCore & bed, std::ostream & out){
			out << bed.toDelimStrWithExtra() << std::endl;
		});
	}else{
		BioDataFileIO<Bed3RecordCore> reader{IoOptions(InOptions(bedFile), outOpts)};
		reader.openIn();
		reader.openOut();
		std::vector<Bed3RecordCore> allBeds;
		Bed3RecordCore record;
		while(reader.readNextRecord(record)){
			allBeds.emplace_back(record);
		}

		std::function<bool(const Bed3RecordCore &, const Bed3RecordCore &)> bedCoordSorterFunc=
				[](const Bed3RecordCore & reg1, const Bed3RecordCore & reg2) {
					if(reg1.chrom_ == reg2.chrom_) {
						if(reg1.chromStart_ == reg2.chromStart_) {
							return reg1.chromEnd_ < reg2.chromEnd_;
						} else {
							return reg1.chromStart_ < reg2.chromStart_;
						}
					} else {
						return reg1.chrom_ < reg2.chrom_;
					}
				};
		if(decending){
			std::sort(allBeds.rbegin(), allBeds.rend(), bedCoordSorterFunc);
		}else{
			njh::sort(allBeds, bedCoordSorterFunc);
		}

		reader.write(allBeds, [](const Bed3RecordCore & bed, std::ostream & out){
			out << bed.toDelimStrWithExtra() << std::endl;
		});
	}
	return 0;
}


int bedExpRunner::bedChangeScoreToLength(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bedFile = "";
	OutOptions outOpts;
	outOpts.outExtention_ = ".bed";
	seqSetUp setUp(inputCommands);
	setUp.setOption(bedFile, "--bed", "Bed file", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	BioDataFileIO<Bed6RecordCore> reader{IoOptions(InOptions(bedFile), outOpts)};
	reader.openIn();
	reader.openOut();
	std::shared_ptr<Bed6RecordCore> b = reader.readNextRecord();
	while(nullptr != b){
		b->score_ = b->length();
		reader.write(*b, [](const Bed6RecordCore & bed, std::ostream & out){
			out << bed.toDelimStrWithExtra() << std::endl;
		});
		b = reader.readNextRecord();
	}
	return 0;
}

int bedExpRunner::bedRenameChromosomes(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path chromKeyTableFnp = "";

	bfs::path bedFile = "";
	OutOptions outOpts;
	outOpts.outExtention_ = ".bed";
	seqSetUp setUp(inputCommands);
	setUp.setOption(bedFile, "--bed", "Bed file", true);
	setUp.setOption(chromKeyTableFnp, "--chromKeyTableFnp", "A key table, no header, 1) current chromosome names in bed file, 2) name to be renamed to", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	table chromKeyTab(chromKeyTableFnp, "\t", false);

	if(2 != chromKeyTab.columnNames_.size()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << chromKeyTableFnp << " should be a table with two columns, first chrom name, second new name" << "\n";
		throw std::runtime_error{ss.str()};
	}

	std::unordered_map<std::string, std::string> chromKey;
	for(const auto & row : chromKeyTab.content_){
		chromKey[row[0]] = row[1];
	}
	BioDataFileIO<Bed3RecordCore> reader{IoOptions(InOptions(bedFile), outOpts)};
	reader.openIn();
	reader.openOut();
	std::shared_ptr<Bed3RecordCore> b = reader.readNextRecord();
	while(nullptr != b){
		if(!njh::in(b->chrom_, chromKey)){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error couldn't find " << b->chrom_ << " in table, options are:" << "\n";
			ss << njh::conToStr(njh::getVecOfMapKeys(chromKey) , ", ") << '\n';
			throw std::runtime_error{ss.str()};
		}
		b->chrom_ = chromKey[b->chrom_];
		reader.write(*b, [](const Bed3RecordCore & bed, std::ostream & out){
			out << bed.toDelimStrWithExtra() << std::endl;
		});
		b = reader.readNextRecord();
	}
	return 0;
}


int bedExpRunner::extendToEndOfChrom(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path chromLengthsTable = "";

	bfs::path bedFile;
	OutOptions outOpts;
	outOpts.outExtention_ = ".bed";
	seqSetUp setUp(inputCommands);
	setUp.setOption(bedFile, "--bed", "Bed file", true);
	setUp.setOption(chromLengthsTable, "--chromLengthsTable", "Chromosome Lengths Table, two columns, first chrom name, second length", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	table chromTab(chromLengthsTable, "\t", false);

	if(2 != chromTab.columnNames_.size()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << chromLengthsTable << " should be a table with two columns, first chrom name, second length" << "\n";
		throw std::runtime_error{ss.str()};
	}
	std::unordered_map<std::string, uint64_t> lengths;
	for(const auto & row : chromTab.content_){
		lengths[row[0]] = njh::StrToNumConverter::stoToNum<uint64_t>(row[1]);
	}
	BioDataFileIO<Bed3RecordCore> reader{IoOptions(InOptions(bedFile), outOpts)};
	reader.openIn();
	reader.openOut();
	std::unordered_map<std::string, std::unordered_map<uint32_t, std::vector<uint32_t>>> alreadyTakenIds;
	std::shared_ptr<Bed3RecordCore> b = reader.readNextRecord();
	while(nullptr != b){
		bool scoreIsLen = false;
		if(b->extraFields_.size() >=2 && njh::strAllDigits(b->extraFields_[1])){
			scoreIsLen = njh::StrToNumConverter::stoToNum<uint32_t>(b->extraFields_[1]) == b->length();
		}
		b->chromEnd_ = njh::mapAt(lengths, b->chrom_);
		if(scoreIsLen){
			b->extraFields_[1] = estd::to_string(b->length());
		}
		reader.write(*b, [](const Bed3RecordCore & bed, std::ostream & out){
			out << bed.toDelimStrWithExtra() << std::endl;
		});
		b = reader.readNextRecord();
	}
	return 0;
}

int bedExpRunner::extendToStartOfChrom(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bedFile;
	OutOptions outOpts;
	outOpts.outExtention_ = ".bed";
	seqSetUp setUp(inputCommands);
	setUp.setOption(bedFile, "--bed", "Bed file", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	BioDataFileIO<Bed3RecordCore> reader{IoOptions(InOptions(bedFile), outOpts)};
	reader.openIn();
	reader.openOut();
	std::unordered_map<std::string, std::unordered_map<uint32_t, std::vector<uint32_t>>> alreadyTakenIds;
	std::shared_ptr<Bed3RecordCore> b = reader.readNextRecord();
	while(nullptr != b){
		bool scoreIsLen = false;
		if(b->extraFields_.size() >=2 && njh::strAllDigits(b->extraFields_[1])){
			scoreIsLen = njh::StrToNumConverter::stoToNum<uint32_t>(b->extraFields_[1]) == b->length();
		}
		b->chromStart_ = 0;
		if(scoreIsLen){
			b->extraFields_[1] = estd::to_string(b->length());
		}
		reader.write(*b, [](const Bed3RecordCore & bed, std::ostream & out){
			out << bed.toDelimStrWithExtra() << std::endl;
		});
		b = reader.readNextRecord();
	}
	return 0;
}




int bedExpRunner::bedAddSmartIDForPlotting(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bedFile;
	OutOptions outOpts;
	outOpts.outExtention_ = ".bed";
	seqSetUp setUp(inputCommands);
	setUp.setOption(bedFile, "--bed", "Bed file", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	BioDataFileIO<Bed3RecordCore> reader{IoOptions(InOptions(bedFile), outOpts)};
	reader.openIn();
	reader.openOut();
	std::unordered_map<std::string, std::unordered_map<uint32_t, std::vector<uint32_t>>> alreadyTakenIds;
	std::shared_ptr<Bed3RecordCore> b = reader.readNextRecord();
	while(nullptr != b){
		auto id = BedUtility::getPlotIDForBed(b, alreadyTakenIds);
		b->extraFields_.emplace_back(estd::to_string(id));
		reader.write(*b, [](const Bed3RecordCore & bed, std::ostream & out){
			out << bed.toDelimStrWithExtra() << std::endl;
		});
		b = reader.readNextRecord();
	}

	return 0;

}

int bedExpRunner::getBestNonOverlapingRegions(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bed1File;
	uint32_t minAllowableOverlap = 0;
	bool bottom = false;
	uint32_t selectionAmount = 10;
	OutOptions outOpts;
	seqSetUp setUp(inputCommands);
	setUp.description_ = "Output minimally overlapping regions with either the top or bottom scores";
	setUp.examples_.emplace_back("MASTERPROGRAM SUBPROGRAM --bed input.bed");
	setUp.examples_.emplace_back("MASTERPROGRAM SUBPROGRAM --bed input.bed --selectionAmount 20 #get top 20 scores that don't overlap");
	setUp.examples_.emplace_back("MASTERPROGRAM SUBPROGRAM --bed input.bed --selectionAmount 20 --minAllowableOverlap 5 #get top 20 scores that overlap by a minimal of 5 bases");
	setUp.examples_.emplace_back("MASTERPROGRAM SUBPROGRAM --bed input.bed --selectionAmount 20 --bottom #get bottom 20 scores that don't overlap");
	setUp.processVerbose();
	setUp.setOption(bed1File, "--bed", "Bed6 (at least first 6 columns) file, regions start from this bed", true);
	setUp.setOption(minAllowableOverlap, "--minAllowableOverlap", "Minimal amount of overlap allowed between regions");
	setUp.setOption(bottom, "--bottom", "Get the lowest scores instead of the highest scores");
	setUp.setOption(selectionAmount, "--selectionAmount", "The number of regions to output");
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	BioDataFileIO<Bed6RecordCore> reader{IoOptions(InOptions(bed1File), outOpts)};
	std::shared_ptr<Bed6RecordCore> reg;
	std::vector<std::shared_ptr<Bed6RecordCore>> startRegs;
	reader.openIn();
	reader.openOut();
	uint32_t amountWritten = 0;
	reg = reader.readNextRecord();
	while(nullptr != reg){
		startRegs.emplace_back(reg);
		reg = reader.readNextRecord();
	}
	if(bottom){
		njh::sort(startRegs, [](const std::shared_ptr<Bed6RecordCore> & b1,const std::shared_ptr<Bed6RecordCore> & b2 ){
			return b1->score_ < b2->score_;
		});
	}else{
		njh::sort(startRegs, [](const std::shared_ptr<Bed6RecordCore> & b1,const std::shared_ptr<Bed6RecordCore> & b2 ){
			return b1->score_ > b2->score_;
		});
	}

	std::vector<std::shared_ptr<GenomicRegion>> writtenRegs;

	for(const auto & reg : startRegs){
		bool foundOverLap = false;
		std::shared_ptr<GenomicRegion> currentReg = std::make_shared<GenomicRegion>(*reg);

		for(const auto & writtenReg : writtenRegs){
			if(writtenReg->getOverlapLen(*currentReg) > minAllowableOverlap){
				foundOverLap = true;
				break;
			}
		}
		if(!foundOverLap){
			amountWritten++;
			reader.write(*reg, [](const Bed6RecordCore & bRecord, std::ostream & out){
				out << bRecord.toDelimStrWithExtra() << "\n";
			});
			writtenRegs.emplace_back(currentReg);
			if(amountWritten >= selectionAmount){
				break;
			}
		}
	}



	return 0;
}

int bedExpRunner::bedCreateSpanningRegions(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bed1File;
	bfs::path bed2File;
	uint32_t minLength = 1;
	uint32_t maxLength = std::numeric_limits<uint32_t>::max();
	bool ignoreStrands = false;
	bool sameStrand = false;
	OutOptions outOpts;
	seqSetUp setUp(inputCommands);
	setUp.description_ = "Create spanning regions between two different beds starting from one set of beds to another set";
	setUp.examples_.emplace_back("MASTERPROGRAM SUBPROGRAM --bed1 starts.bed --bed2 ends.bed");
	setUp.examples_.emplace_back("MASTERPROGRAM SUBPROGRAM --bed1 starts.bed --bed2 ends.bed --ignoreStrands");
	setUp.processVerbose();
	setUp.setOption(bed1File, "--bed1", "Bed6 (at least first 6 columns) file, regions start from this bed", true);
	setUp.setOption(bed2File, "--bed2", "Bed6 (at least first 6 columns) file, regions end   from this bed", true);
	setUp.setOption(ignoreStrands, "--ignoreStrands", "Ignore the strands the two regions are coming from, strand will default to bed1 in this case");
	setUp.setOption(sameStrand, "--sameStrand", "Create spanning regions only if the regions are on the same strand, default is opposite strands");
	setUp.setOption(minLength, "--minLength", "Minimum length of the spanning region");
	setUp.setOption(maxLength, "--maxLength", "Maximum length of the spanning region");
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	BioDataFileIO<Bed6RecordCore> reader{IoOptions(InOptions(bed1File), outOpts)};
	std::shared_ptr<Bed6RecordCore> reg;
	std::vector<std::shared_ptr<Bed6RecordCore>> startRegs;
	reader.openIn();
	reader.openOut();

	reg = reader.readNextRecord();
	while(nullptr != reg){
		startRegs.emplace_back(reg);
		reg = reader.readNextRecord();

	}

	std::vector<std::shared_ptr<Bed6RecordCore>> endRegions;
	if(bed1File == bed2File){
		endRegions = startRegs;
	}else{
		BioDataFileIO<Bed6RecordCore> reader2{IoOptions(InOptions(bed2File))};
		reader2.openIn();
		std::shared_ptr<Bed6RecordCore> reg2;

		reg2 = reader2.readNextRecord();
		while(nullptr != reg2){
			endRegions.emplace_back(reg2);
			reg2 = reader2.readNextRecord();
		}
	}

	njh::sort(startRegs, [](const std::shared_ptr<Bed6RecordCore> & b1,const std::shared_ptr<Bed6RecordCore> & b2 ){
		return b1->chrom_ == b2->chrom_ ? b1->chromStart_ < b2->chromStart_ : b1->chrom_ < b2->chrom_;
	});
	njh::sort(endRegions, [](const std::shared_ptr<Bed6RecordCore> & b1,const std::shared_ptr<Bed6RecordCore> & b2 ){
		return b1->chrom_ == b2->chrom_ ? b1->chromStart_ < b2->chromStart_ : b1->chrom_ < b2->chrom_;
	});


	for(const auto & start : startRegs){
		for(const auto & end : endRegions){
			if(end->chrom_ != start->chrom_){
				continue;
			}
			if(end->chromStart_ > start->chromStart_ &&
					end->chromEnd_ - start->chromStart_ > maxLength){
				break;
			}else if(end->chromStart_ < start->chromStart_ ){
				continue;
			}else{
				uint32_t length = end->chromEnd_ - start->chromStart_ ;
				bool write = false;
				if(ignoreStrands){
					if(length >=minLength && length <= maxLength){
						write = true;
					}
				}else{
					if((sameStrand && end->strand_ == start->strand_) ||
							(!sameStrand && end->strand_ != start->strand_)){
						if(length >=minLength && length <= maxLength){
							write = true;
						}
					}
				}
				if(write){
					Bed6RecordCore out(start->chrom_,
							start->chromStart_,
							end->chromEnd_,
							njh::pasteAsStr(start->name_, "-", end->name_),
							length,
							start->strand_);
					reader.write(out, [](const Bed6RecordCore & bRecord, std::ostream & out){
						out << bRecord.toDelimStrWithExtra() << "\n";
					});
				}
			}
		}
	}

	return 0;
}


int bedExpRunner::bedToggleStrand(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bedFile;
	OutOptions outOpts;
	seqSetUp setUp(inputCommands);
	setUp.description_ = "Simply switch the strand of the region in the bed file";
	setUp.examples_.emplace_back("MASTERPROGRAM SUBPROGRAM --bed inputFile.bed");
	setUp.processVerbose();
	setUp.setOption(bedFile, "--bed", "Bed6 (at least first 6 columns) file", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	BioDataFileIO<Bed6RecordCore> reader{IoOptions(InOptions(bedFile), outOpts)};
	std::shared_ptr<Bed6RecordCore> reg;
	reader.openIn();
	reader.openOut();

	reg = reader.readNextRecord();
	while(nullptr != reg){
		if ('-' == reg->strand_) {
			reg->strand_ = '+';
		} else {
			reg->strand_ = '-';
		}
		reader.write(*reg, [](const Bed6RecordCore & bRecord, std::ostream & out){
			out << bRecord.toDelimStrWithExtra() << "\n";
		});
		reg = reader.readNextRecord();
	}

	return 0;
}

int bedExpRunner::extendUpstreamRegion(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bedFile;
	OutOptions outOpts;
	uint32_t length = 0;
	seqSetUp setUp(inputCommands);
	setUp.description_ = "Extend the upstream region from the bed file depending on the strand orientation so it's truly upstream";
	setUp.examples_.emplace_back("MASTERPROGRAM SUBPROGRAM --bed inputFile.bed --length");
	setUp.examples_.emplace_back("MASTERPROGRAM SUBPROGRAM --bed inputFile.bed --length");
	setUp.processVerbose();
	setUp.setOption(bedFile, "--bed", "Bed6 file extract upstream region from", true);
	setUp.setOption(length, "--length", "how much upstream to output", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	BioDataFileIO<Bed6RecordCore> reader{IoOptions(InOptions(bedFile), outOpts)};
	std::shared_ptr<Bed6RecordCore> reg;
	reader.openIn();
	reader.openOut();
	std::unordered_map<std::string, std::vector<std::shared_ptr<Bed6RecordCore>>> beds;
	reg = reader.readNextRecord();
	std::vector<std::shared_ptr<Bed6RecordCore>> bedVec;
	while(nullptr != reg){
		GenomicRegion inputRegion(*reg);
		GenomicRegion upstreamRegion = inputRegion;
		if (upstreamRegion.reverseSrand_) {
			upstreamRegion.end_ = inputRegion.end_ + length;
		} else {
			upstreamRegion.start_ = length < inputRegion.start_ ? inputRegion.start_ - length : 0;
		}
		reader.write(upstreamRegion.genBedRecordCore(), [](const Bed6RecordCore & bRecord, std::ostream & out){
			out << bRecord.toDelimStrWithExtra() << "\n";
		});
		reg = reader.readNextRecord();
	}
	return 0;
}

int bedExpRunner::trimUpstreamRegion(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bedFile;
	OutOptions outOpts;
	uint32_t length = 0;
	seqSetUp setUp(inputCommands);
	setUp.description_ = "Trim the upstream region from the bed file depending on the strand orientation so it's truly upstream";
	setUp.examples_.emplace_back("MASTERPROGRAM SUBPROGRAM --bed inputFile.bed --length");
	setUp.examples_.emplace_back("MASTERPROGRAM SUBPROGRAM --bed inputFile.bed --length");
	setUp.processVerbose();
	setUp.setOption(bedFile, "--bed", "Bed6 file extract upstream region from", true);
	setUp.setOption(length, "--length", "how much upstream to output", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	BioDataFileIO<Bed6RecordCore> reader{IoOptions(InOptions(bedFile), outOpts)};
	std::shared_ptr<Bed6RecordCore> reg;
	reader.openIn();
	reader.openOut();
	std::unordered_map<std::string, std::vector<std::shared_ptr<Bed6RecordCore>>> beds;
	reg = reader.readNextRecord();
	std::vector<std::shared_ptr<Bed6RecordCore>> bedVec;
	while(nullptr != reg){
		GenomicRegion inputRegion(*reg);
		GenomicRegion upstreamRegion = inputRegion;
		if (upstreamRegion.reverseSrand_) {
			upstreamRegion.end_ = length < inputRegion.end_ ? inputRegion.end_ - length : 0;
		} else {
			upstreamRegion.start_ = inputRegion.start_ + length;
		}
		if(upstreamRegion.end_ <= upstreamRegion.start_){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << inputRegion.createUidFromCoordsStrand() << " now has an end that is before or equal to it's start: " << upstreamRegion.createUidFromCoordsStrand()<< "\n";
			throw std::runtime_error{ss.str()};
		}
		reader.write(upstreamRegion.genBedRecordCore(), [](const Bed6RecordCore & bRecord, std::ostream & out){
			out << bRecord.toDelimStrWithExtra() << "\n";
		});
		reg = reader.readNextRecord();
	}
	return 0;
}




int bedExpRunner::extendDownstreamRegion(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bedFile;
	OutOptions outOpts;
	uint32_t length = 0;
	seqSetUp setUp(inputCommands);
	setUp.description_ = "Extend the downstream region from the bed file depending on the strand orientation so it's truly downstream";
	setUp.examples_.emplace_back("MASTERPROGRAM SUBPROGRAM --bed inputFile.bed --length");
	setUp.examples_.emplace_back("MASTERPROGRAM SUBPROGRAM --bed inputFile.bed --length --include 100 #include 100 bases of the region from which you are extracting");
	setUp.processVerbose();
	setUp.setOption(bedFile, "--bed", "Bed6 file extract downstream region from", true);
	setUp.setOption(length, "--length", "how much downstream to output", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	BioDataFileIO<Bed6RecordCore> reader{IoOptions(InOptions(bedFile), outOpts)};
	std::shared_ptr<Bed6RecordCore> reg;
	reader.openIn();
	reader.openOut();
	std::unordered_map<std::string, std::vector<std::shared_ptr<Bed6RecordCore>>> beds;
	reg = reader.readNextRecord();
	std::vector<std::shared_ptr<Bed6RecordCore>> bedVec;
	while(nullptr != reg){
		GenomicRegion inputRegion(*reg);
		GenomicRegion downstreamRegion = inputRegion;
		if (downstreamRegion.reverseSrand_) {
			downstreamRegion.start_ = length < inputRegion.start_ ? inputRegion.start_ - length : 0;
		} else {
			downstreamRegion.end_ = inputRegion.end_ + length;
		}


		reader.write(downstreamRegion.genBedRecordCore(), [](const Bed6RecordCore & bRecord, std::ostream & out){
			out << bRecord.toDelimStrWithExtra() << "\n";
		});
		reg = reader.readNextRecord();
	}

	return 0;
}


int bedExpRunner::trimDownstreamRegion(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bedFile;
	OutOptions outOpts;
	uint32_t length = 0;
	seqSetUp setUp(inputCommands);
	setUp.description_ = "Extend the downstream region from the bed file depending on the strand orientation so it's truly downstream";
	setUp.examples_.emplace_back("MASTERPROGRAM SUBPROGRAM --bed inputFile.bed --length");
	setUp.examples_.emplace_back("MASTERPROGRAM SUBPROGRAM --bed inputFile.bed --length --include 100 #include 100 bases of the region from which you are extracting");
	setUp.processVerbose();
	setUp.setOption(bedFile, "--bed", "Bed6 file extract downstream region from", true);
	setUp.setOption(length, "--length", "how much downstream to output", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	BioDataFileIO<Bed6RecordCore> reader{IoOptions(InOptions(bedFile), outOpts)};
	std::shared_ptr<Bed6RecordCore> reg;
	reader.openIn();
	reader.openOut();
	std::unordered_map<std::string, std::vector<std::shared_ptr<Bed6RecordCore>>> beds;
	reg = reader.readNextRecord();
	std::vector<std::shared_ptr<Bed6RecordCore>> bedVec;
	while(nullptr != reg){
		GenomicRegion inputRegion(*reg);
		GenomicRegion downstreamRegion = inputRegion;
		if (downstreamRegion.reverseSrand_) {
			downstreamRegion.start_ = inputRegion.start_ + length;
		} else {
			downstreamRegion.end_ = length < inputRegion.end_ ? inputRegion.end_ - length : 0;
		}

		if(downstreamRegion.end_ <= downstreamRegion.start_){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << inputRegion.createUidFromCoordsStrand() << " now has an end that is before or equal to it's start: " << downstreamRegion.createUidFromCoordsStrand()<< "\n";
			throw std::runtime_error{ss.str()};
		}

		reader.write(downstreamRegion.genBedRecordCore(), [](const Bed6RecordCore & bRecord, std::ostream & out){
			out << bRecord.toDelimStrWithExtra() << "\n";
		});
		reg = reader.readNextRecord();
	}

	return 0;
}



int bedExpRunner::getUpstreamRegion(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bedFile;
	OutOptions outOpts;
	uint32_t include = 0;
	uint32_t length = 0;
	seqSetUp setUp(inputCommands);
	setUp.description_ = "Extract the upstream region from the bed file depending on the strand orientation so it's truly upstream";
	setUp.examples_.emplace_back("MASTERPROGRAM SUBPROGRAM --bed inputFile.bed --length");
	setUp.examples_.emplace_back("MASTERPROGRAM SUBPROGRAM --bed inputFile.bed --length --include 100 #include 100 bases of the region from which you are extracting");
	setUp.processVerbose();
	setUp.setOption(bedFile, "--bed", "Bed6 file extract upstream region from", true);
	setUp.setOption(length, "--length", "how much upstream to output", true);
	setUp.setOption(include, "--include", "include this much of the region to extract the upstream region from");
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	BioDataFileIO<Bed6RecordCore> reader{IoOptions(InOptions(bedFile), outOpts)};
	std::shared_ptr<Bed6RecordCore> reg;
	reader.openIn();
	reader.openOut();
	std::unordered_map<std::string, std::vector<std::shared_ptr<Bed6RecordCore>>> beds;
	reg = reader.readNextRecord();
	std::vector<std::shared_ptr<Bed6RecordCore>> bedVec;
	while(nullptr != reg){
		GenomicRegion inputRegion(*reg);
		GenomicRegion upstreamRegion = inputRegion;
		upstreamRegion.uid_ = "upstream-" + upstreamRegion.uid_;
		if (upstreamRegion.reverseSrand_) {
			upstreamRegion.start_ = include < inputRegion.end_ ? inputRegion.end_ - include : 0;
			upstreamRegion.end_ = inputRegion.end_ + length;
		} else {
			upstreamRegion.start_ = length < inputRegion.start_ ? inputRegion.start_ - length : 0;
			upstreamRegion.end_ = inputRegion.start_ + include;
		}
		reader.write(upstreamRegion.genBedRecordCore(), [](const Bed6RecordCore & bRecord, std::ostream & out){
			out << bRecord.toDelimStrWithExtra() << "\n";
		});
		reg = reader.readNextRecord();
	}

	return 0;
}

int bedExpRunner::getDownstreamRegion(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bedFile;
	OutOptions outOpts;
	uint32_t include = 0;
	uint32_t length = 0;
	seqSetUp setUp(inputCommands);
	setUp.description_ = "Extract the downstream region from the bed file depending on the strand orientation so it's truly downstream";
	setUp.examples_.emplace_back("MASTERPROGRAM SUBPROGRAM --bed inputFile.bed --length");
	setUp.examples_.emplace_back("MASTERPROGRAM SUBPROGRAM --bed inputFile.bed --length --include 100 #include 100 bases of the region from which you are extracting");
	setUp.processVerbose();
	setUp.setOption(bedFile, "--bed", "Bed6 file extract downstream region from", true);
	setUp.setOption(length, "--length", "how much downstream to output", true);
	setUp.setOption(include, "--include", "include this much of the region to extract the downstream region from");
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	BioDataFileIO<Bed6RecordCore> reader{IoOptions(InOptions(bedFile), outOpts)};
	std::shared_ptr<Bed6RecordCore> reg;
	reader.openIn();
	reader.openOut();
	std::unordered_map<std::string, std::vector<std::shared_ptr<Bed6RecordCore>>> beds;
	reg = reader.readNextRecord();
	std::vector<std::shared_ptr<Bed6RecordCore>> bedVec;
	while(nullptr != reg){
		GenomicRegion inputRegion(*reg);
		GenomicRegion downstreamRegion = inputRegion;
		downstreamRegion.uid_ = "downstream-" + downstreamRegion.uid_;
		if (downstreamRegion.reverseSrand_) {
			downstreamRegion.start_ = length < inputRegion.start_ ? inputRegion.start_ - length : 0;
			downstreamRegion.end_ = inputRegion.start_ + include;
		} else {
			downstreamRegion.start_ = include < inputRegion.end_ ? inputRegion.end_ - include : 0;
			downstreamRegion.end_ = inputRegion.end_ + length;
		}
		reader.write(downstreamRegion.genBedRecordCore(), [](const Bed6RecordCore & bRecord, std::ostream & out){
			out << bRecord.toDelimStrWithExtra() << "\n";
		});
		reg = reader.readNextRecord();
	}

	return 0;
}

int bedExpRunner::bedRenameWithCoords(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bedFile;
	OutOptions outOpts;
	bfs::path keyOutFnp = "";

	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.setOption(keyOutFnp, "--key", "Name of a file for key for rename");
	setUp.setOption(bedFile, "--bed", "Bed6 file to rename", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);
	if ("STDIN" == bedFile && "" == outOpts.outFilename_) {
		outOpts.outFilename_ = "STDOUT";
	}
	if ("" == outOpts.outFilename_) {
		outOpts.outFilename_ = njh::files::prependFileBasename(bedFile, "renamed_");
	}
	OutOptions keyOpts(keyOutFnp);
	keyOpts.transferOverwriteOpts(outOpts);
	std::unique_ptr<OutputStream> keyOut;
	if("" != keyOutFnp){
		keyOut = std::make_unique<OutputStream>(keyOpts);
		(*keyOut) << "oldName\tnewName" << std::endl;
	}

	BioDataFileIO<Bed6RecordCore> reader{IoOptions(InOptions(bedFile), outOpts)};
	std::shared_ptr<Bed6RecordCore> reg;

	reader.openIn();
	reader.openOut();
	std::unordered_map<std::string, std::vector<std::shared_ptr<Bed6RecordCore>>> beds;
	reg = reader.readNextRecord();
	std::vector<std::shared_ptr<Bed6RecordCore>> bedVec;
	while(nullptr != reg){
		GenomicRegion gRegion(*reg);
		gRegion.setUidWtihCoordsStrand();
		if("" != keyOutFnp){
			(*keyOut) << reg->name_ << "\t" << gRegion.uid_ << std::endl;
		}
		reg->name_ = gRegion.uid_;
		reader.write(*reg, [](const Bed6RecordCore & bRecord, std::ostream & out){
			out << bRecord.toDelimStrWithExtra() << "\n";
		});
		reg = reader.readNextRecord();
	}

	return 0;
}

int bedExpRunner::bedRenameRepeatUids(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bedFile;
	OutOptions outOpts;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.setOption(bedFile, "--bed", "Bed file to parse", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	BioDataFileIO<Bed6RecordCore> reader{IoOptions(InOptions(bedFile), outOpts)};
	std::shared_ptr<Bed6RecordCore> reg;

	reader.openIn();
	reader.openOut();
	std::unordered_map<std::string, std::vector<std::shared_ptr<Bed6RecordCore>>> beds;
	reg = reader.readNextRecord();
	std::vector<std::shared_ptr<Bed6RecordCore>> bedVec;
	while(nullptr != reg){
		beds[reg->name_].emplace_back(reg);
		bedVec.emplace_back(reg);
		reg = reader.readNextRecord();
	}
	for(auto & bed : beds	){
		if(bed.second.size() > 1){
			uint32_t count = 0;
			for(auto & subBed : bed.second){
				if(count > 0){
					subBed->name_ = subBed->name_ + "." + leftPadNumStr<uint32_t>(count, bed.second.size());
				}
				++count;
			}
		}
	}
	for(const auto & bed : bedVec){
		reader.write(*bed, [](const Bed6RecordCore & bRecord, std::ostream & out){
			out << bRecord.toDelimStrWithExtra() << "\n";
		});
	}
	return 0;
}




int bedExpRunner::getDegreeOfOverlappingBedRegions(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bedFile;
	bfs::path intersectWithBed;
	OutOptions outOpts;
	uint32_t overlapLen = 1;
	bool ignoreSameID = false;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.setOption(overlapLen, "--overlapLen", "The minimum overlap length to be considered overlap, 0 means one region must be completely found in another");
	setUp.setOption(bedFile, "--bed", "Bed file to parse", true);
	setUp.setOption(ignoreSameID, "--ignoreSameID", "Don't compare regions with the same id if ids present");
	setUp.setOption(intersectWithBed, "--intersectWithBed", "Bed file to intersect with", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	BioDataFileIO<Bed6RecordCore> reader{IoOptions(InOptions(bedFile), outOpts)};
	auto intersectingBed3s = getBed3s(intersectWithBed);
	auto intersectingBed6s = convertBed3ToBed6(intersectingBed3s);
	Bed6RecordCore reg;

	reader.openIn();
	reader.openOut();
	(*reader.out_) << "#chrom\tstart\tend\treg1name\tscore\tstrand\toverlap\treg1OverlapPerc\treg2Name\treg2OverlapPerc" << std::endl;
	while(reader.readNextRecord(reg)){
		std::vector<Bed6RecordCore> intersectingRegions;
		for(const auto & bed : intersectingBed6s){
			if(ignoreSameID && !reg.extraFields_.empty() &&  bed->name_ == reg.extraFields_[0]){
				continue;
			}
			//a overlapLen of 0 indicates to find only regions that completely overlap with each other
			if(overlapLen == 0){
				auto overlapAmount = reg.getOverlapLen(*bed);
				if(bed->length() == overlapAmount ||reg.length() == overlapAmount){
					intersectingRegions.emplace_back(*bed);
				}
			}else{
				if(reg.overlaps(*bed, overlapLen)){
					intersectingRegions.emplace_back(*bed);
				}
			}
		}
		if(!intersectingRegions.empty()){
			reader.write(reg, [&intersectingRegions](const Bed6RecordCore & bedReg, std::ostream & out){
				for(const auto & ireg : intersectingRegions){
					auto overlapAmount = ireg.getOverlapLen(bedReg);
					out << bedReg.toDelimStr()
							<< "\t" << overlapAmount
							<< "\t" << static_cast<double>(overlapAmount)/bedReg.length()
							<< "\t" << ireg.name_
							<< "\t" << static_cast<double>(overlapAmount)/ireg.length()
							<< std::endl;
				}
			});
		}
	}
	return 0;
}

int bedExpRunner::getOverlappingBedRegions(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bedFile;
	bfs::path intersectWithBed;
	OutOptions outOpts;
	uint32_t overlapLen = 1;
	bool ignoreSameID = false;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.setOption(overlapLen, "--overlapLen", "The minimum overlap length to be considered overlap, 0 means one region must be completely found in another");
	setUp.setOption(bedFile, "--bed", "Bed file to parse", true);
	setUp.setOption(ignoreSameID, "--ignoreSameID", "Don't compare regions with the same id if ids present");
	setUp.setOption(intersectWithBed, "--intersectWithBed", "Bed file to intersect with", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	BioDataFileIO<Bed3RecordCore> reader{IoOptions(InOptions(bedFile), outOpts)};
	auto intersectingBed3s = getBed3s(intersectWithBed);
	auto intersectingBed6s = convertBed3ToBed6(intersectingBed3s);
	Bed3RecordCore reg;

	reader.openIn();
	reader.openOut();
	while(reader.readNextRecord(reg)){
		VecStr intersectingRegions;
		for(const auto & bed : intersectingBed6s){
			if(ignoreSameID && !reg.extraFields_.empty() &&  bed->name_ == reg.extraFields_[0]){
				continue;
			}
			//a overlapLen of 0 indicates to find only regions that completely overlap with each other
			if(overlapLen == 0){
				auto overlapAmount = reg.getOverlapLen(*bed);
				if(bed->length() == overlapAmount ||reg.length() == overlapAmount){
					intersectingRegions.emplace_back(bed->name_);
				}
			}else{
				if(reg.overlaps(*bed, overlapLen)){
					intersectingRegions.emplace_back(bed->name_);
				}
			}
		}
		if(!intersectingRegions.empty()){
			reader.write(reg, [&intersectingRegions](const Bed3RecordCore & bedReg, std::ostream & out){
				out << bedReg.toDelimStrWithExtra()
						<< "\t" << njh::conToStr(intersectingRegions, ",")<< std::endl;
			});
		}
	}
	return 0;
}

int bedExpRunner::getNonOverlappingBedRegions(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bedFile;
	bfs::path intersectWithBed;
	OutOptions outOpts;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.setOption(bedFile, "--bed", "Bed file to parse", true);
	setUp.setOption(intersectWithBed, "--intersectWithBed", "Bed file to intersect with", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);


	BioDataFileIO<Bed3RecordCore> reader{IoOptions(InOptions(bedFile), outOpts)};
	auto intersectingBed3s = getBed3s(intersectWithBed);
	auto intersectingBed6s = convertBed3ToBed6(intersectingBed3s);

	Bed3RecordCore reg;


	reader.openIn();
	reader.openOut();

	while(reader.readNextRecord(reg)){
		bool foundIntersection = false;
		for(const auto & bed : intersectingBed6s){
			if(reg.overlaps(*bed, 1)){
				foundIntersection = true;
				break;
			}
		}
		if(!foundIntersection){
			reader.write(reg, [](const Bed3RecordCore & bedReg, std::ostream & out){
				out << bedReg.toDelimStrWithExtra() << std::endl;
			});
		}
	}

	return 0;
}

int bedExpRunner::filterBedRecordsByLength(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bedFile;
	OutOptions outOpts;
	uint32_t minLen = 0;
	uint32_t maxLen = std::numeric_limits<uint32_t>::max();
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.setOption(minLen,  "--minLen", "Min length to allow (inclusive)");
	setUp.setOption(maxLen,  "--maxLen", "Max length to allow (inclusive)");
	setUp.setOption(bedFile, "--bed", "Bed file to parse", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);


	BioDataFileIO<Bed3RecordCore> reader{IoOptions(InOptions(bedFile), outOpts)};
	Bed3RecordCore reg;


	reader.openIn();
	reader.openOut();

	while(reader.readNextRecord(reg)){
		if(reg.length()>= minLen && reg.length() <= maxLen){
			reader.write(reg, [](const Bed3RecordCore & bedReg, std::ostream & out){
				out << bedReg.toDelimStrWithExtra() << std::endl;
			});
		}
	}

	return 0;
}

int bedExpRunner::bed3ToBed6(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bedFile;
	OutOptions outOpts;
	bool reverseStrand = false;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.setOption(bedFile, "--bed", "Bed file to parse", true);
	setUp.setOption(reverseStrand, "--reverseStrand", "Output in Reverse Strand");
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);


	BioDataFileIO<Bed3RecordCore> reader{IoOptions(InOptions(bedFile))};
	Bed3RecordCore reg;
	BioDataFileIO<Bed6RecordCore> writer{IoOptions(outOpts)};

	reader.openIn();
	writer.openOut();

	while(reader.readNextRecord(reg)){

		auto outBedRec = GenomicRegion(reg).genBedRecordCore();
		if(reverseStrand){
			outBedRec.strand_ = '-';
		}
		writer.write(outBedRec, []( const Bed6RecordCore & bedReg, std::ostream & out){
			out << bedReg.toDelimStr() << std::endl;
		});
	}

	return 0;
}



int bedExpRunner::extendBedRegions(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bedFile;
	OutOptions outOpts;
	uint32_t left = 0;
	uint32_t right = 0;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();

	setUp.setOption(left, "--left", "Expand regions to the left this much");
	setUp.setOption(right, "--right", "Expand regions to the right this much");
	setUp.setOption(bedFile, "--bed", "Bed file to parse", true);

	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);


	BioDataFileIO<Bed3RecordCore> reader{IoOptions(InOptions(bedFile))};
	reader.openIn();

	std::ofstream outFile;
	std::ostream out(determineOutBuf(outFile, outOpts));
	Bed3RecordCore reg;
	while(reader.readNextRecord(reg)){
		bool scoreIsLen = false;
		if(reg.extraFields_.size() >=3){
			if(njh::strAllDigits(reg.extraFields_[1])){
				auto score = njh::StrToNumConverter::stoToNum<double>(reg.extraFields_[1]);
				if(reg.length() == score){
					scoreIsLen = true;
				}
			}
		}
		reg.chromStart_ = (reg.chromStart_ > left ? reg.chromStart_ - left : 0);
		reg.chromEnd_ += right;
		//reg.score_ = reg.length();
		if(reg.extraFields_.size() >=3){
			if(scoreIsLen){
				reg.extraFields_[1] = estd::to_string(reg.length());
			}
		}
		out << reg.toDelimStrWithExtra() << std::endl;
	}
	return 0;
}

int bedExpRunner::trimBedRegions(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bedFile;
	OutOptions outOpts;
	uint32_t left = 0;
	uint32_t right = 0;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();

	setUp.setOption(left, "--left", "Trim regions to the left this much");
	setUp.setOption(right, "--right", "Trim regions to the right this much");
	setUp.setOption(bedFile, "--bed", "Bed file to parse", true);

	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);


	BioDataFileIO<Bed3RecordCore> reader{IoOptions(InOptions(bedFile))};
	reader.openIn();

	OutputStream outFile(outOpts);
	Bed3RecordCore reg;
	while(reader.readNextRecord(reg)){
		bool scoreIsLen = false;
		if(reg.extraFields_.size() >=3){
			if(njh::strAllDigits(reg.extraFields_[1])){
				auto score = njh::StrToNumConverter::stoToNum<double>(reg.extraFields_[1]);
				if(reg.length() == score){
					scoreIsLen = true;
				}
			}
		}
		auto uid = njh::pasteAsStr(reg.chrom_, "-", reg.chromStart_, "-", reg.chromEnd_);
		reg.chromStart_ += left;
		if(right > reg.chromEnd_ ){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "right: " << right << " larger than end position: " << reg.chromEnd_<< "\n";
			throw std::runtime_error{ss.str()};
		}
		reg.chromEnd_ -= right;
		if(reg.chromEnd_ <= reg.chromStart_){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "region: " << uid << " will now have a end that is before the start: " <<  njh::pasteAsStr(reg.chrom_, "-", reg.chromStart_, "-", reg.chromEnd_)<< "\n";
			throw std::runtime_error{ss.str()};
		}
		if(reg.extraFields_.size() >=3){
			if(scoreIsLen){
				reg.extraFields_[1] = estd::to_string(reg.length());
			}
		}
		outFile << reg.toDelimStrWithExtra() << std::endl;
	}
	return 0;
}



int bedExpRunner::getCloseBedRegions(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bedFile;
	bfs::path targetRegionFnp;
	OutOptions outOpts;
	uint32_t maxDistance = 1000;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();

	setUp.setOption(maxDistance, "--maxDistance", "Max distance away from target region to output region", true);
	setUp.setOption(bedFile, "--bed", "Bed file to parse", true);
	setUp.setOption(targetRegionFnp, "--targetRegion", "Region to get other regions close to", true);

	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);


	auto regions =	bedPtrsToGenomicRegs(getBeds(bedFile.string()));
	auto targetRegionInput =	bedPtrsToGenomicRegs(getBeds(targetRegionFnp.string()));
	if(targetRegionInput.empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", target region file " << targetRegionFnp << " was empty" << "\n";
		throw std::runtime_error{ss.str()};
	}
	const GenomicRegion targetRegion = targetRegionInput.front();
	std::vector<std::shared_ptr<GenomicRegion>> outRegions;
	for(const auto & reg : regions){
		if(reg.chrom_ == targetRegion.chrom_){
			std::vector<uint32_t> distances{uAbsdiff(reg.start_, targetRegion.start_),
																			uAbsdiff(reg.end_, targetRegion.start_),
																			uAbsdiff(reg.start_, targetRegion.end_),
																			uAbsdiff(reg.end_, targetRegion.end_)
			};
			if(vectorMinimum(distances) <= maxDistance){
				outRegions.emplace_back(std::make_shared<GenomicRegion>(reg));
			}
		}
	}

	std::ofstream outFile;
	std::ostream out(determineOutBuf(outFile, outOpts));

	for(const auto & reg : outRegions){
		out << reg->genBedRecordCore().toDelimStr() << std::endl;
	}
	return 0;
}

int bedExpRunner::bedUnqiue(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bedFile;
	bfs::path outputDir;
	OutOptions outOpts;
	bool ignoreStrand = false;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.setOption(bedFile, "--bed", "Bed file to collapse", true);
	setUp.setOption(ignoreStrand, "--ignoreStrand", "Ignore strand compare");

	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);


	auto regions = njh::convert<std::shared_ptr<Bed6RecordCore>, std::shared_ptr<GenomicRegion>>(getBeds(bedFile.string()),
				[](const std::shared_ptr<Bed6RecordCore> & bed) {return std::make_shared<GenomicRegion>(*bed);});
	njh::sort(regions, [](const std::shared_ptr<GenomicRegion> & reg1, const std::shared_ptr<GenomicRegion> & reg2){
		return reg1->createUidFromCoords() < reg2->createUidFromCoords();
	});
	std::vector<std::shared_ptr<GenomicRegion>> outRegions;
	for(const auto & reg : regions){
		if(!outRegions.empty() &&
				outRegions.back()->sameRegion(*reg) &&
				(ignoreStrand || outRegions.back()->reverseSrand_ == reg->reverseSrand_)){
			outRegions.back()->uid_.append("-" + reg->uid_);
		}else{
			outRegions.emplace_back(reg);
		}
	}

	std::ofstream outFile;
	std::ostream out(determineOutBuf(outFile, outOpts));

	for(const auto & reg : outRegions){
		out << reg->genBedRecordCore().toDelimStr() << std::endl;
	}
	return 0;
}



int bedExpRunner::getLongestHomopolymerLengthInRegion(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path filename = "";
	OutOptions outOpts(bfs::path(""), ".bed");
	bfs::path twoBitFilename = "";

	seqSetUp setUp(inputCommands);
	setUp.setOption(twoBitFilename, "--twoBit", "File path of the 2bit file", true);
	setUp.setOption(filename, "--bed", "BED6 file", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);
	OutputStream out(outOpts);
	TwoBit::TwoBitFile twoBitFile(twoBitFilename);
	auto seqNames = twoBitFile.sequenceNames();
	BioDataFileIO<Bed6RecordCore> bedReader{IoOptions(InOptions(filename))};
	bedReader.openIn();
	Bed6RecordCore record;
	out << "#chrom\tstart\tend\tname\tscore\tstrand\tHPbase\tHPposition\tHPlength" << std::endl;
	while (bedReader.readNextRecord(record)) {
		if (!njh::in(record.chrom_, seqNames)) {
			std::cerr << "chromosome name not found in seq names, skipping"
					<< std::endl;
			std::cerr << "chr: " << record.chrom_ << std::endl;
			std::cerr << "possibleNames: " << vectorToString(seqNames, ",")
					<< std::endl;
		} else {
			std::string seq = "";
			twoBitFile[record.chrom_]->getSequence(seq, record.chromStart_,
					record.chromEnd_);
			if (record.reverseStrand()) {
				seq = seqUtil::reverseComplement(seq, "DNA");
			}
			readObject seqObj{seqInfo{record.name_, seq}};
			seqObj.createCondensedSeq();
			uint32_t longestHomopolymerLen = 0;
			uint32_t longestHomopolymerPos = 0;
			char longestHomopolymerBase = ' ';
			uint32_t seqPos = 0;
			for(const auto homPos : iter::range(seqObj.condensedSeqCount.size())){
				if(seqObj.condensedSeqCount[homPos] > longestHomopolymerLen){
					longestHomopolymerLen = seqObj.condensedSeqCount[homPos];
					longestHomopolymerBase = seqObj.condensedSeq[homPos];
					longestHomopolymerPos = seqPos;
				}
				seqPos += seqObj.condensedSeqCount[homPos];
			}
			out << record.toDelimStr()
					<< "\t" << longestHomopolymerBase
					<< "\t" << GenomicRegion(record).getRelativePositionFromStartStrandAware(longestHomopolymerPos)
					<< "\t" << longestHomopolymerLen << std::endl;
		}
	}
	return 0;
}


int bedExpRunner::getNumOfNsEntropyGCContentOfRegion(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path filename = "";
	OutOptions outOpts(bfs::path(""), ".bed");
	bfs::path twoBitFilename = "";
	bool raw = false;
	seqSetUp setUp(inputCommands);
	setUp.setOption(twoBitFilename, "--twoBit", "File path of the 2bit file", true);
	setUp.setOption(filename, "--bed", "BED6 file", true);
	setUp.setOption(raw, "--raw", "raw");

	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);
	OutputStream out(outOpts);
	TwoBit::TwoBitFile twoBitFile(twoBitFilename);
	auto seqNames = twoBitFile.sequenceNames();
	BioDataFileIO<Bed6RecordCore> bedReader{IoOptions(InOptions(filename))};
	bedReader.openIn();
	Bed6RecordCore record;
	out << "#chrom\tstart\tend\tname\tscore\tstrand\tns\tnFrac\tentropy\tgcBases\tgcContent" << std::endl;
	while (bedReader.readNextRecord(record)) {
		if (!njh::in(record.chrom_, seqNames)) {
			std::cerr << "chromosome name not found in seq names, skipping"
					<< std::endl;
			std::cerr << "chr: " << record.chrom_ << std::endl;
			std::cerr << "possibleNames: " << vectorToString(seqNames, ",")
					<< std::endl;
		} else {
			std::string seq = "";
			twoBitFile[record.chrom_]->getSequence(seq, record.chromStart_,
					record.chromEnd_);
			if (record.reverseStrand()) {
				seq = seqUtil::reverseComplement(seq, "DNA");
			}
			if(!raw){
				njh::for_each(seq,[](char & c){c = toupper(c);});
			}
			DNABaseCounter counter;
			counter.increase(seq);
			counter.resetAlphabet(false);
			out << record.toDelimStr()
					    << "\t" << counter.bases_['N'] + counter.bases_['n']
							<< "\t" << (counter.bases_['N'] + counter.bases_['n'])/static_cast<double>(seq.size())
							<< "\t" << counter.computeEntrophyBasedOffAlph(4)
							<< "\t" << counter.getGcCount()
							<< "\t" << counter.getGcCount()/static_cast<double>(seq.size())
							<< std::endl;
		}
	}
	return 0;
}


int bedExpRunner::getNumOfNsOfRegion(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path filename = "";
	OutOptions outOpts(bfs::path(""), ".bed");
	bfs::path twoBitFilename = "";

	seqSetUp setUp(inputCommands);
	setUp.setOption(twoBitFilename, "--twoBit", "File path of the 2bit file", true);
	setUp.setOption(filename, "--bed", "BED6 file", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);
	OutputStream out(outOpts);
	TwoBit::TwoBitFile twoBitFile(twoBitFilename);
	auto seqNames = twoBitFile.sequenceNames();
	BioDataFileIO<Bed6RecordCore> bedReader{IoOptions(InOptions(filename))};
	bedReader.openIn();
	Bed6RecordCore record;
	out << "#chrom\tstart\tend\tname\tscore\tstrand\tns\tnFrac" << std::endl;
	while (bedReader.readNextRecord(record)) {
		if (!njh::in(record.chrom_, seqNames)) {
			std::cerr << "chromosome name not found in seq names, skipping"
					<< std::endl;
			std::cerr << "chr: " << record.chrom_ << std::endl;
			std::cerr << "possibleNames: " << vectorToString(seqNames, ",")
					<< std::endl;
		} else {
			std::string seq = "";
			twoBitFile[record.chrom_]->getSequence(seq, record.chromStart_,
					record.chromEnd_);
			if (record.reverseStrand()) {
				seq = seqUtil::reverseComplement(seq, "DNA");
			}
			DNABaseCounter counter;
			counter.increase(seq);
			counter.resetAlphabet(false);
			out << record.toDelimStr()
					<< "\t" << counter.bases_['N'] + counter.bases_['n']
					<< "\t" << (counter.bases_['N'] + counter.bases_['n'])/static_cast<double>(seq.size()) << std::endl;
		}
	}
	return 0;
}

int bedExpRunner::getGCContentOfRegion(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path filename = "";
	OutOptions outOpts(bfs::path(""), ".bed");
	bfs::path twoBitFilename = "";

	seqSetUp setUp(inputCommands);
	setUp.setOption(twoBitFilename, "--twoBit", "File path of the 2bit file", true);
	setUp.setOption(filename, "--bed", "BED6 file", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);
	OutputStream out(outOpts);
	TwoBit::TwoBitFile twoBitFile(twoBitFilename);
	auto seqNames = twoBitFile.sequenceNames();
	BioDataFileIO<Bed6RecordCore> bedReader{IoOptions(InOptions(filename))};
	bedReader.openIn();
	Bed6RecordCore record;
	out << "#chrom\tstart\tend\tname\tscore\tstrand\tgcBases\tgcContent" << std::endl;
	while (bedReader.readNextRecord(record)) {
		if (!njh::in(record.chrom_, seqNames)) {
			std::cerr << "chromosome name not found in seq names, skipping"
					<< std::endl;
			std::cerr << "chr: " << record.chrom_ << std::endl;
			std::cerr << "possibleNames: " << vectorToString(seqNames, ",")
					<< std::endl;
		} else {
			std::string seq = "";
			twoBitFile[record.chrom_]->getSequence(seq, record.chromStart_,
					record.chromEnd_);
			if (record.reverseStrand()) {
				seq = seqUtil::reverseComplement(seq, "DNA");
			}
			DNABaseCounter counter;
			counter.increase(seq);
			out << record.toDelimStr()
					<< "\t" << counter.getGcCount()
					<< "\t" << counter.getGcCount()/static_cast<double>(seq.size()) << std::endl;
		}
	}
	return 0;
}

int bedExpRunner::getEntropyOfRegion(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path filename = "";
	OutOptions outOpts(bfs::path(""), ".bed");
	bfs::path twoBitFilename = "";
	bool raw = false;
	seqSetUp setUp(inputCommands);
	setUp.setOption(twoBitFilename, "--twoBit", "File path of the 2bit file", true);
	setUp.setOption(filename, "--bed", "BED6 file", true);
	setUp.setOption(raw, "--raw", "raw");

	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);
	OutputStream out(outOpts);
	TwoBit::TwoBitFile twoBitFile(twoBitFilename);
	auto seqNames = twoBitFile.sequenceNames();
	BioDataFileIO<Bed6RecordCore> bedReader{IoOptions(InOptions(filename))};
	bedReader.openIn();
	Bed6RecordCore record;
	out << "#chrom\tstart\tend\tname\tscore\tstrand\tentropy" << std::endl;
	while (bedReader.readNextRecord(record)) {
		if (!njh::in(record.chrom_, seqNames)) {
			std::cerr << "chromosome name not found in seq names, skipping"
					<< std::endl;
			std::cerr << "chr: " << record.chrom_ << std::endl;
			std::cerr << "possibleNames: " << vectorToString(seqNames, ",")
					<< std::endl;
		} else {
			std::string seq = "";
			twoBitFile[record.chrom_]->getSequence(seq, record.chromStart_,
					record.chromEnd_);
			if (record.reverseStrand()) {
				seq = seqUtil::reverseComplement(seq, "DNA");
			}
			if(!raw){
				njh::for_each(seq,[](char & c){c = toupper(c);});
			}
			DNABaseCounter counter;
			counter.increase(seq);
			counter.resetAlphabet(false);
			out << record.toDelimStr()
					<< "\t" << counter.computeEntrophyBasedOffAlph(4) << std::endl;
		}
	}
	return 0;
}



int bedExpRunner::getSeqFromTwoBit(const njh::progutils::CmdArgs & inputCommands) {
	OutOptions outOpts(bfs::path(""));
	outOpts.outExtention_ = ".fasta";
	bfs::path twoBitFilename = "";
	std::string chrom = "";
	uint32_t start = std::numeric_limits<uint32_t>::max();
	uint32_t end = std::numeric_limits<uint32_t>::max();
	bool reverseStrand = false;
	seqSetUp setUp(inputCommands);
	setUp.setOption(twoBitFilename, "--twoBit,--2bit", "File path of the 2bit file", true);
	setUp.setOption(chrom, "--chrom", "chromosome", true);
	setUp.setOption(start, "--start", "start", true);
	end = start + 1;
	setUp.setOption(end, "--end", "end");
	setUp.setOption(reverseStrand, "--reverseStrand", "Get Reverse Strand");

	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);
	OutputStream out(outOpts);
	TwoBit::TwoBitFile twoBitFile(twoBitFilename);
	auto seqNames = twoBitFile.sequenceNames();

	Bed6RecordCore record(chrom, start, end, njh::pasteAsStr(chrom, "-", start, "-", end), end - start, (reverseStrand? '-' : '+'));

	if (!njh::in(record.chrom_, seqNames)) {
		std::cerr << "chromosome name not found in seq names, skipping"
				<< std::endl;
		std::cerr << "chr: " << record.chrom_ << std::endl;
		std::cerr << "possibleNames: " << vectorToString(seqNames, ",")
				<< std::endl;
	} else {
		std::string seq = "";
		twoBitFile[record.chrom_]->getSequence(seq, record.chromStart_,
				record.chromEnd_);
		if (record.reverseStrand()) {
			seq = seqUtil::reverseComplement(seq, "DNA");
		}
		out << ">" << record.name_ << std::endl;
		out << seq << std::endl;
	}

	return 0;
}

int bedExpRunner::getFastaWithBed(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path filename = "";
	OutOptions outOpts(bfs::path(""));
	outOpts.outExtention_ = ".fasta";
	bfs::path twoBitFilename = "";

	seqSetUp setUp(inputCommands);
	setUp.setOption(twoBitFilename, "--twoBit", "File path of the 2bit file", true);
	setUp.setOption(filename, "--bed", "BED6 file", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);
	OutputStream out(outOpts);
	TwoBit::TwoBitFile twoBitFile(twoBitFilename);
	auto seqNames = twoBitFile.sequenceNames();
	BioDataFileIO<Bed6RecordCore> bedReader{IoOptions(InOptions(filename))};
	bedReader.openIn();
	Bed6RecordCore record;
	while (bedReader.readNextRecord(record)) {
		if (!njh::in(record.chrom_, seqNames)) {
			std::cerr << "chromosome name not found in seq names, skipping"
					<< std::endl;
			std::cerr << "chr: " << record.chrom_ << std::endl;
			std::cerr << "possibleNames: " << vectorToString(seqNames, ",")
					<< std::endl;
		} else {
			std::string seq = "";
			twoBitFile[record.chrom_]->getSequence(seq, record.chromStart_,
					record.chromEnd_);
			if (record.reverseStrand()) {
				seq = seqUtil::reverseComplement(seq, "DNA");
			}
			out << ">" << record.name_ << std::endl;
			out << seq << std::endl;
		}
	}
	return 0;
}


int bedExpRunner::separateOutRecordsInBedFile(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bedFile;
	bfs::path outputDir;
	OutOptions outOpts;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.setOption(bedFile, "--bed", "Bed file to split", true);
	setUp.setOption(outOpts.overWriteFile_, "--overWrite", "Overwrite output bed files if they exist");
	outputDir = bedFile.parent_path();
	setUp.setOption(outputDir, "--outputDir", "Output directory");
	setUp.finishSetUp(std::cout);
	auto beds = getBeds(bedFile.string());
	auto regions = bedPtrsToGenomicRegs(beds);
	std::unique_ptr<OutputStream> keyOut;
	for(const auto regPos : iter::range(regions.size())){
		const auto & reg = regions[regPos];
		auto bRec = reg.genBedRecordCore();
		OutOptions currentOpts(njh::files::make_path(outputDir,bRec.name_ + ".bed"));
		currentOpts.transferOverwriteOpts(outOpts);
		OutputStream out(currentOpts);
		out << bRec.toDelimStr() << "\n";
	}
	return 0;
}




int bedExpRunner::extractBedRecordsWithName(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path filename = "";
	OutOptions outOpts(bfs::path(""), ".bed");
	std::string names = "";

	seqSetUp setUp(inputCommands);
	setUp.setOption(names, "--names", "Names", true);
	setUp.setOption(filename, "--bed", "BED6 file", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);
	OutputStream out(outOpts);
	auto namesVec = getInputValues(names, ",");

	BioDataFileIO<Bed6RecordCore> bedReader{IoOptions(InOptions(filename))};
	bedReader.openIn();
	Bed6RecordCore record;
	while (bedReader.readNextRecord(record)) {
		if (njh::in(record.name_, namesVec)) {
			out << record.toDelimStrWithExtra() << std::endl;
		}
	}
	return 0;
}


int bedExpRunner::splitBedFile(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bedFile;
	bfs::path outputDir;
	OutOptions outOpts;
	uint32_t step = 100;
	uint32_t windowSize = 300;
	bool separateFiles = false;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.setOption(bedFile, "--bed", "Bed file to split", true);
	outputDir = bedFile.parent_path();
	setUp.setOption(separateFiles, "--separateFiles", "Create a separate bed file for each bed region");
	setUp.setOption(outputDir, "--outputDir", "Output directory if outputting separate bed files");
	setUp.processWritingOptions(outOpts);
	setUp.setOption(step, "--step", "Step size");
	setUp.setOption(windowSize, "--windowSize", "Window size");
	setUp.finishSetUp(std::cout);

	if(!separateFiles){
		outOpts.throwIfOutExistsNoOverWrite(__PRETTY_FUNCTION__);
	}
	auto regions = gatherRegions(bedFile.string(), "", setUp.pars_.verbose_);
	if(separateFiles){
		table splitRegions(VecStr{"chrom", "start", "end", "name", "score", "strand"});
		for(const auto & region : regions){
			if(region.getLen() >= windowSize){
				splitRegions.addRow(region.chrom_, region.start_, region.end_, njh::pasteAsStr(region.chrom_, "_", region.start_, "_", region.end_), region.end_ - region.start_, region.reverseSrand_ ? '-':'+');
			}else{
				for(const auto pos : iter::range<size_t>(region.start_, region.end_ - windowSize, step)){
					auto actualEnd = std::min(pos + windowSize, region.end_);
					auto basesLeft = region.end_ - actualEnd;
					if(basesLeft > 0 && basesLeft < step){
						if(basesLeft > (windowSize/2.0)){
							splitRegions.addRow(region.chrom_, pos, actualEnd, njh::pasteAsStr(region.chrom_, "_", pos, "_", actualEnd), actualEnd - pos, region.reverseSrand_ ? '-':'+');
							splitRegions.addRow(region.chrom_, pos + step, region.end_, njh::pasteAsStr(region.chrom_, "_", pos + step, "_", region.end_), region.end_ - (pos + step), region.reverseSrand_ ? '-':'+');
						}else{
							actualEnd += basesLeft;
							splitRegions.addRow(region.chrom_, pos, actualEnd, njh::pasteAsStr(region.chrom_, "_", pos, "_", actualEnd), actualEnd - pos, region.reverseSrand_ ? '-':'+');
						}
					}else{
						splitRegions.addRow(region.chrom_, pos, actualEnd, njh::pasteAsStr(region.chrom_, "_", pos, "_", actualEnd), actualEnd - pos, region.reverseSrand_ ? '-':'+');
					}
				}
			}
		}
		for(const auto & row : splitRegions.content_){
			OutOptions currentOpts(njh::files::make_path(outputDir, row[splitRegions.getColPos("name")] + ".bed"));
			currentOpts.overWriteFile_ = outOpts.overWriteFile_;
			std::ofstream out;
			currentOpts.openFile(out);
			out << njh::conToStr(row, "\t") << "\n";
		}
	}else{
		OutputStream out(outOpts);
		for(const auto & region : regions){
			if (region.getLen() < windowSize) {
				out
						<< GenomicRegion(region.uid_, region.chrom_, region.start_, region.end_,
								region.reverseSrand_).genBedRecordCore().toDelimStr()
						<< std::endl;
			} else {
				for(const auto pos : iter::range<size_t>(region.start_, region.end_ - windowSize, step)){

					auto actualEnd = std::min(pos + windowSize, region.end_);
					auto basesLeft = region.end_ - actualEnd;
					if (basesLeft > 0 && basesLeft < step) {
						if (basesLeft > (windowSize / 2.0)) {
							out
									<< GenomicRegion("", region.chrom_, pos, actualEnd,
											region.reverseSrand_).genBedRecordCore().toDelimStr()
									<< std::endl;
							out
									<< GenomicRegion("", region.chrom_, pos + step, region.end_,
											region.reverseSrand_).genBedRecordCore().toDelimStr()
									<< std::endl;
						} else {
							actualEnd += basesLeft;
							out
									<< GenomicRegion("", region.chrom_, pos, actualEnd,
											region.reverseSrand_).genBedRecordCore().toDelimStr()
									<< std::endl;
						}
					} else {
						out
								<< GenomicRegion("", region.chrom_, pos, actualEnd,
										region.reverseSrand_).genBedRecordCore().toDelimStr()
								<< std::endl;
					}
				}
			}
		}
	}

	return 0;
}


int bedExpRunner::reorientBasedOnSingleReadsOrientationCounts(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bedFnp;
	bfs::path bamFnp;
	uint32_t numThreads = 1;
	OutOptions outOpts(bfs::path(""), ".bed");

	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(numThreads, "--numThreads", "numThreads");

	setUp.setOption(bedFnp, "--bedFnp", "bed file", true);
	setUp.setOption(bamFnp, "--bamFnp", "bam file", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);

	njh::files::checkExistenceThrow(bedFnp, __PRETTY_FUNCTION__);
	njh::files::checkExistenceThrow(bamFnp, __PRETTY_FUNCTION__);
	bool debug = setUp.pars_.debug_	;
	concurrent::BamReaderPool bamPool(bamFnp, numThreads);
	bamPool.openBamFile();
	OutputStream out(outOpts);
	auto beds = getBeds(bedFnp);
	njh::concurrent::LockableQueue<std::shared_ptr<Bed6RecordCore>> bedQueue(beds);
	std::function<void()> refineRegions =[&bamPool,&bedQueue,&debug](){
		std::shared_ptr<Bed6RecordCore> region;
		BamTools::BamAlignment bAln;
		while(bedQueue.getVal(region)){
			auto bamReader = bamPool.popReader();
			setBamFileRegionThrow(*bamReader, *region);
			uint32_t forCount = 0;
			uint32_t revCount  = 0;
			while(bamReader->GetNextAlignmentCore(bAln)){
				if(!bAln.IsPaired()){
					if(bAln.IsReverseStrand() ){
						++revCount;
					}
					if(!bAln.IsReverseStrand() ){
						++forCount;
					}
				}
			}
			if((revCount + forCount >0 ) && static_cast<double>(revCount)/(forCount + revCount) > 0.50){
				region->strand_ = '-';
			}else{
				region->strand_ = '+';
			}
			if(debug){
				MetaDataInName meta;
				meta.addMeta("forCount", forCount);
				meta.addMeta("revCount", revCount);
				region->extraFields_.emplace_back(meta.createMetaName());
			}
		}
	};

	njh::concurrent::runVoidFunctionThreaded(refineRegions, numThreads);
	for(const auto & bed : beds){
		out << bed->toDelimStrWithExtra() << std::endl;
	}


	return 0;
}


} /* namespace njhseq */
