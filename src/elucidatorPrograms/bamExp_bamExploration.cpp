/*
 * bamExp_bamExploration.cpp
 *
 *  Created on: May 16, 2017
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


#include "bamExp.hpp"
#include <TwoBit.h>

#include "elucidator/BamToolsUtils.h"
#include <njhseq/objects/BioDataObject.h>
#include <njhseq/GenomeUtils/GenomeExtraction/ParsingAlignmentInfo/AlignmentResults.hpp>

namespace njhseq {





int bamExpRunner::BamGetImproperPairCounts(const njh::progutils::CmdArgs & inputCommands){
	bfs::path bedFile = "";
	OutOptions outOpts(bfs::path(""), ".tab.txt");
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.setOption(bedFile, "--bed", "Bed file");
	setUp.processReadInNames({"--bam"}, true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	BamTools::BamReader bReader;
	bReader.Open(setUp.pars_.ioOptions_.firstName_.string());
	checkBamOpenThrow(bReader, setUp.pars_.ioOptions_.firstName_.string());
	loadBamIndexThrow(bReader);
	BamTools::BamAlignment bAln;

	if("" != bedFile){
		auto regions = bedPtrsToGenomicRegs(getBeds(bedFile));
		setBamFileRegionThrow(bReader, regions.front());
	}
	auto refData = bReader.GetReferenceData();
	OutputStream out(outOpts);

	//key1=isPaired, key2=isMapped, key3=isMateMapped, key4=ProperPair, value = count
	std::unordered_map<bool, std::unordered_map<bool, std::unordered_map<bool, std::unordered_map<bool, uint64_t>>>> pairedMappedMateMappedProperCounts;
	while(bReader.GetNextAlignment(bAln)){
		if(!bAln.IsPrimaryAlignment()){
			continue;
		}
		++pairedMappedMateMappedProperCounts[bAln.IsPaired()][bAln.IsMapped()][bAln.IsMateMapped()][bAln.IsProperPair()];
	}

	out << "isPaired\tisMapped\tisMateMapped\tisProperPair\tcount" << std::endl;
	for(const auto & paired : pairedMappedMateMappedProperCounts){
		for(const auto & mapped : paired.second){
			for(const auto & mateMapped : mapped.second){
				for(const auto & properPair : mateMapped.second){
					out << njh::boolToStr(paired.first)
					<< "\t" << njh::boolToStr(mapped.first)
					<< "\t" << njh::boolToStr(mateMapped.first)
					<< "\t" << njh::boolToStr(properPair.first)
					<< "\t" << properPair.second << std::endl;
				}
			}
		}
	}


	return 0;
}



int bamExpRunner::bamToFastqAlns(const njh::progutils::CmdArgs & inputCommands){
	bfs::path bedFile = "";
	OutOptions outOpts(bfs::path(""));
	bfs::path twoBitFnp = "";
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.setOption(bedFile, "--bedFile", "Bed file");
	setUp.setOption(twoBitFnp, "--2bit", "Two bit file of the genome aligned to", true);
	setUp.processReadInNames({"--bam"}, true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	BamTools::BamReader bReader;
	bReader.Open(setUp.pars_.ioOptions_.firstName_.string());
	checkBamOpenThrow(bReader, setUp.pars_.ioOptions_.firstName_.string());
	loadBamIndexThrow(bReader);
	BamTools::BamAlignment bAln;

	if("" != bedFile){
		auto regions = bedPtrsToGenomicRegs(getBeds(bedFile));
		setBamFileRegionThrow(bReader, regions.front());
	}
	auto refData = bReader.GetReferenceData();
	TwoBit::TwoBitFile tReader(twoBitFnp);
	auto seqOutOpts = SeqIOOptions::genFastqOut(outOpts.outFilename_);
	seqOutOpts.out_.transferOverwriteOpts(outOpts);
	SeqOutput seqOut(seqOutOpts);
	while(bReader.GetNextAlignment(bAln)){
		if(!bAln.IsPrimaryAlignment() || !bAln.IsMapped()){
			continue;
		}
		AlignmentResults results(bAln, refData);
		results.setRefSeq(tReader);
		results.setAlignedObjects();
		seqOut.openWrite(results.refSeqAligned_);
		seqOut.openWrite(results.alnSeqAligned_);
	}
	if(setUp.pars_.verbose_){
		std::cout << std::endl;
	}

	return 0;
}


int bamExpRunner::countUnMappedMateStatus(const njh::progutils::CmdArgs & inputCommands){
	bfs::path bedFile = "";
	OutOptions outOpts(bfs::path(""));
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.setOption(bedFile, "--bedFile", "Bed file");
	setUp.processReadInNames({"--bam"}, true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	BamTools::BamReader bReader;
	bReader.Open(setUp.pars_.ioOptions_.firstName_.string());
	checkBamOpenThrow(bReader, setUp.pars_.ioOptions_.firstName_.string());
	loadBamIndexThrow(bReader);
	BamTools::BamAlignment bAln;

	if("" != bedFile){
		auto regions = bedPtrsToGenomicRegs(getBeds(bedFile));
		setBamFileRegionThrow(bReader, regions.front());
	}
	std::unordered_map<bool, std::unordered_map<bool, std::unordered_map<bool, std::unordered_map<bool, std::unordered_map<bool, uint32_t>>>>> mateOrientationCount;
	uint32_t count = 0;
	while(bReader.GetNextAlignmentCore(bAln)){
		++count;
		if(setUp.pars_.verbose_){
			std::cout << "\r" << count;
			std::cout.flush();
		}
		if(!bAln.IsPrimaryAlignment()){
			continue;
		}
		if(bAln.IsPaired()){
			bool concordant = false;
			if(bAln.RefID == bAln.MateRefID && std::abs(bAln.InsertSize) < 1000){
				concordant = true;
			}
			++mateOrientationCount[concordant][bAln.IsMapped()][bAln.IsMateMapped()][bAln.IsReverseStrand()][bAln.IsMateReverseStrand()];
		}
	}
	if(setUp.pars_.verbose_){
		std::cout << std::endl;
	}
	table outTab(VecStr { "Concordant", "IsMapped", "IsMateMapped",
			"IsReverseStrand", "IsMateReverseStrand", "count" });
	for (const auto & concordant : mateOrientationCount) {
		for (const auto & mapped : concordant.second) {
			for (const auto & mateMapped : mapped.second) {
				for (const auto & revStrand : mateMapped.second) {
					for (const auto & mateRevStrand : revStrand.second) {
						outTab.addRow(concordant.first, mapped.first, mateMapped.first,
								revStrand.first, mateRevStrand.first, mateRevStrand.second);
					}
				}
			}
		}
	}

	std::ofstream outFile;
	std::ostream out(determineOutBuf(outFile, outOpts));
	outTab.outPutContents(out, "\t");


	return 0;
}


int bamExpRunner::BamFindDifferenceInUnmappedFileIndexPosition(const njh::progutils::CmdArgs & inputCommands){
//	bfs::path bedFile = "";
	std::string name = "";
	OutOptions outOpts(bfs::path(""));
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
//	setUp.setOption(bedFile, "--bedFile", "Bed file", true);
	setUp.processReadInNames({"--bam"}, true);
	setUp.processWritingOptions(outOpts);
	setUp.setOption(name, "--name", "Name");
	setUp.finishSetUp(std::cout);


	OutputStream out(outOpts);
	struct MatePositions{
		int64_t mappedMatePos_ = std::numeric_limits<int64_t>::max();
		int64_t unmappedMatePos_ = std::numeric_limits<int64_t>::max();
		int64_t mappedMateRefId = std::numeric_limits<int64_t>::max();
		int64_t unmappedMateRefId = std::numeric_limits<int64_t>::max();
		int64_t mappedMateMappedPosition = std::numeric_limits<int64_t>::max();
		int64_t unmappedMateMapppedPosition = std::numeric_limits<int64_t>::max();
	};
	std::unordered_map<std::string, MatePositions> matePositionsInfo;
	{//first pass
		BamTools::BamReader bReader;
		bReader.Open(setUp.pars_.ioOptions_.firstName_.string());
		checkBamOpenThrow(bReader, setUp.pars_.ioOptions_.firstName_.string());
		loadBamIndexThrow(bReader);
		BamTools::BamAlignment bAln;
		int64_t count = 0;
		while(bReader.GetNextAlignment(bAln)){
			if(!bAln.IsPrimaryAlignment()){
				continue;
			}
			if(setUp.pars_.verbose_){
				std::cout << "\r" << count;
				std::cout.flush();
			}
			if(bAln.IsPaired()){
				if(bAln.IsMapped() && !bAln.IsMateMapped()){
					matePositionsInfo[bAln.Name].mappedMatePos_ = count;
					matePositionsInfo[bAln.Name].mappedMateRefId = bAln.RefID;
					matePositionsInfo[bAln.Name].mappedMateMappedPosition = bAln.Position;
				}else if(!bAln.IsMapped() && bAln.IsMateMapped()){
					matePositionsInfo[bAln.Name].unmappedMatePos_ = count;
					matePositionsInfo[bAln.Name].unmappedMateRefId = bAln.RefID;
					matePositionsInfo[bAln.Name].unmappedMateMapppedPosition = bAln.Position;
				}
			}
			++count;
		}
	}
	if(setUp.pars_.verbose_){
		std::cout << std::endl;
	}

	std::map<int32_t, uint32_t> countsOfPositionDifference;
	std::map<int32_t, uint32_t> countsOfMappPositionDifference;

	uint32_t unmappedMateNotFound = 0;
	uint32_t mappedMateNotFound = 0;
	uint32_t bothNotFound = 0;
	VecStr mappedMateNotFoundNames;
	VecStr unmappedMateNotFoundNames;
	uint32_t refIdDifferent = 0;
	for (const auto & info : matePositionsInfo) {
		if (info.second.unmappedMatePos_ == std::numeric_limits<int64_t>::max()
				&& info.second.mappedMatePos_ == std::numeric_limits<int64_t>::max()) {
			++bothNotFound;
		} else if (info.second.unmappedMatePos_
				== std::numeric_limits<int64_t>::max()) {
			unmappedMateNotFoundNames.emplace_back(info.first);
			++unmappedMateNotFound;
		} else if (info.second.mappedMatePos_
				== std::numeric_limits<int64_t>::max()) {
			mappedMateNotFoundNames.emplace_back(info.first);
			++mappedMateNotFound;
		} else {
			++countsOfPositionDifference[info.second.unmappedMatePos_
					- info.second.mappedMatePos_];
			if(info.second.mappedMateRefId == info.second.unmappedMateRefId){
				++countsOfMappPositionDifference[info.second.unmappedMateMapppedPosition - info.second.mappedMateMappedPosition];
			}else{
				++refIdDifferent;
			}
		}
	}
	out << "FileIndexDifference\tcount" << "\n";
	for(const auto & diff : countsOfPositionDifference){
		out << diff.first << "\t" << diff.second << std::endl;
	}
	out << std::endl;
	out << "MappedDifference\tcount" << "\n";
	for(const auto & diff : countsOfMappPositionDifference){
		out << diff.first << "\t" << diff.second << std::endl;
	}
	out << std::endl;
	out << "bothNotFound" << "\t" << bothNotFound << std::endl;
	out << "unmappedMateNotFound" << "\t" << unmappedMateNotFound << std::endl;
	out << "mappedMateNotFound" << "\t" << mappedMateNotFound << std::endl << std::endl;

	out << "refIdDifferent" << "\t" << refIdDifferent << std::endl;
	out << "mappedMateNotFoundNames" << "\t" << njh::conToStr(mappedMateNotFoundNames, ",") << std::endl;
	out << "unmappedMateNotFoundNames" << "\t" << njh::conToStr(unmappedMateNotFoundNames, ",") << std::endl;


	return 0;
}

int bamExpRunner::BamGetFileIndexPositionOfName(const njh::progutils::CmdArgs & inputCommands){
	OutOptions outOpts(bfs::path(""), ".tab.txt");
	VecStr names;
	seqSetUp setUp(inputCommands);
	setUp.processReadInNames({"--bam"}, true);
	setUp.setOption(names, "--names", "names", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	BamTools::BamReader bReader;
	bReader.Open(setUp.pars_.ioOptions_.firstName_.string());
	checkBamOpenThrow(bReader, setUp.pars_.ioOptions_.firstName_.string());
	//loadBamIndexThrow(bReader);
	auto refData = bReader.GetReferenceData();
	BamTools::BamAlignment bAln;
	uint32_t count = 0;
	OutputStream out(outOpts);
	out << "FileIndexPosition\tName\tPosition\tEndPosition\tQuerySize\tAlnSize\tRefId\tRefName\tMateRefID\tMatePosition\tIsMateMapped\tIsMapped\tIsPaired\tIsPrimaryAlignment\tIsFirstMate\tDuplicate\tRevComp\tProperPair\tMapQ" << "\n";
	while(bReader.GetNextAlignment(bAln)){
		if (njh::in(bAln.Name, names)) {
			out << count
					<< "\t" << bAln.Name
					<< "\t" << bAln.Position
					<< "\t" << bAln.GetEndPosition()
					<< "\t" << bAln.QueryBases.size()
					<< "\t" << bAln.GetEndPosition() - bAln.Position
					<< "\t" << bAln.RefID
					<< "\t" << (bAln.RefID < 0 || bAln.RefID >= refData.size() ? "*" : refData[bAln.RefID].RefName)
					<< "\t" << (bAln.MateRefID < 0 || bAln.MateRefID >= refData.size() ? "*" : refData[bAln.MateRefID].RefName)
					<< "\t" << bAln.MatePosition
					<< "\t" << njh::boolToStr(bAln.IsMateMapped())
					<< "\t" << njh::boolToStr(bAln.IsMapped())
					<< "\t" << njh::boolToStr(bAln.IsPaired())
					<< "\t" << njh::boolToStr(bAln.IsPrimaryAlignment())
					<< "\t" << njh::boolToStr(bAln.IsFirstMate())
				  << "\t" << njh::boolToStr(bAln.IsDuplicate())
					<< "\t" << njh::boolToStr(bAln.IsReverseStrand())
					<< "\t" << njh::boolToStr(bAln.IsProperPair())
			    << "\t" << static_cast<uint32_t>(bAln.MapQuality) << std::endl;
		}
		++count;
	}
	return 0;
}

int bamExpRunner::getMateMapStatus(const njh::progutils::CmdArgs & inputCommands){
	bfs::path bedFile = "";
	OutOptions outOpts(bfs::path("mate_map_status.tab.txt"));
	seqSetUp setUp(inputCommands);
	setUp.processWritingOptions(outOpts);
	setUp.setOption(bedFile, "--bedFile", "Bed file", true);
	setUp.processReadInNames({"--bam"}, true);
	setUp.finishSetUp(std::cout);

	std::ofstream outFile;
	outOpts.openFile(outFile);
	BamTools::BamReader bReader;
	bReader.Open(setUp.pars_.ioOptions_.firstName_.string());
	checkBamOpenThrow(bReader, setUp.pars_.ioOptions_.firstName_.string());
	loadBamIndexThrow(bReader);
	BamTools::BamAlignment bAln;
	auto regions = bedPtrsToGenomicRegs(getBeds(bedFile));

	setBamFileRegionThrow(bReader, regions.front());
	struct MateMapStatusCount{
		uint32_t concordantMapping = 0;
		uint32_t crossChromeMapping = 0;
		uint32_t mateNotMapped = 0;

		uint32_t total() const{
			return concordantMapping + crossChromeMapping + mateNotMapped;
		}
	};
	std::unordered_map<uint32_t, MateMapStatusCount> pairCounts;
	while(bReader.GetNextAlignment(bAln)){
		if(bAln.Position > regions.front().start_ ){
			if(bAln.IsPaired()){
				if(bAln.IsMapped()){
					if(bAln.IsMateMapped()){
						if(bAln.RefID != bAln.MateRefID){
							//cross chrom mapping
							++pairCounts[bAln.Position].crossChromeMapping;
						}else{
							//same chrome mapping
							++pairCounts[bAln.Position].concordantMapping;
						}
					}else{
						//mate not mapping
						++pairCounts[bAln.Position].mateNotMapped;
					}
				}
			}
		}
	}
	auto keys = getVectorOfMapKeys(pairCounts);

	njh::sort(keys);
	std::vector<uint32_t> allPositions(regions.front().getLen());
	njh::iota<uint32_t>(allPositions, regions.front().start_);
	outFile << njh::conToStr(VecStr{"position", "concordantMapping", "crossChromeMapping", "mateNotMapped", "total"}, "\t")<< std::endl;
	for(const auto & key : allPositions){
		MateMapStatusCount current;
		if(njh::in(key, keys)){
			current = pairCounts[key];
		}
		outFile << njh::conToStr(toVecStr(key, current.concordantMapping, current.crossChromeMapping, current.mateNotMapped, current.total()) , "\t")<< std::endl;
	}
	return 0;
}

//
int bamExpRunner::getInsertSizeChanges(const njh::progutils::CmdArgs & inputCommands){
	bfs::path bedFile = "";
	OutOptions outOpts(bfs::path("insertSizes.tab.txt"));
	seqSetUp setUp(inputCommands);
	setUp.processWritingOptions(outOpts);
	setUp.setOption(bedFile, "--bedFile", "Bed file", true);
	setUp.processReadInNames({"--bam"}, true);
	setUp.finishSetUp(std::cout);

	std::ofstream outFile;
	outOpts.openFile(outFile);


	BamTools::BamReader bReader;
	bReader.Open(setUp.pars_.ioOptions_.firstName_.string());
	checkBamOpenThrow(bReader, setUp.pars_.ioOptions_.firstName_.string());
	loadBamIndexThrow(bReader);
	BamTools::BamAlignment bAln;

	auto regions = bedPtrsToGenomicRegs(getBeds(bedFile));

	setBamFileRegionThrow(bReader, regions.front());
	std::unordered_map<uint32_t, std::unordered_map<uint32_t, uint32_t>> pairCounts;
	while(bReader.GetNextAlignment(bAln)){
		if(bAln.Position > regions.front().start_ && bAln.GetEndPosition() < regions.front().end_){
			if(bAln.IsPaired()){
				if(bAln.IsMapped() && bAln.IsMateMapped() && bAln.Position < bAln.MatePosition){
					++pairCounts[bAln.Position][std::abs(bAln.InsertSize)];
				}
			}
		}
	}

	auto keys = getVectorOfMapKeys(pairCounts);

	njh::sort(keys);
	std::vector<uint32_t> allPositions(regions.front().getLen());
	njh::iota<uint32_t>(allPositions, regions.front().start_);
	outFile << njh::conToStr(VecStr{"position", "insertSize", "counts"}, "\t")<< std::endl;

	for (const auto & key : allPositions) {
		if (njh::in(key, keys)) {
			for (const auto & inSize : pairCounts[key]) {
				outFile << key << "\t" << inSize.first << "\t" << inSize.second
						<< std::endl;
			}
		} else {
			outFile << key << "\t" << "NA" << "\t" << "NA"
					<< std::endl;
		}
	}


	return 0;
}

int bamExpRunner::getMapQualityCounts(const njh::progutils::CmdArgs & inputCommands){
	bfs::path bedFile = "";
	OutOptions outOpts(bfs::path("mapQualCounts.tab.txt"));
	seqSetUp setUp(inputCommands);
	setUp.processWritingOptions(outOpts);
	setUp.setOption(bedFile, "--bedFile", "Bed file", true);
	setUp.processReadInNames({"--bam"}, true);
	setUp.finishSetUp(std::cout);

	std::ofstream outFile;
	outOpts.openFile(outFile);


	BamTools::BamReader bReader;
	bReader.Open(setUp.pars_.ioOptions_.firstName_.string());
	checkBamOpenThrow(bReader, setUp.pars_.ioOptions_.firstName_.string());
	loadBamIndexThrow(bReader);
	BamTools::BamAlignment bAln;

	auto regions = bedPtrsToGenomicRegs(getBeds(bedFile));

	setBamFileRegionThrow(bReader, regions.front());
	std::unordered_map<uint32_t, std::unordered_map<int32_t, uint32_t>> pairCounts;
	while(bReader.GetNextAlignment(bAln)){
		if(bAln.Position > regions.front().start_ && bAln.GetEndPosition() < regions.front().end_){
			if(bAln.IsMapped()){
				++pairCounts[bAln.Position][bAln.MapQuality];
			}
		}
	}

	auto keys = getVectorOfMapKeys(pairCounts);

	njh::sort(keys);
	std::vector<uint32_t> allPositions(regions.front().getLen());
	njh::iota<uint32_t>(allPositions, regions.front().start_);
	outFile << njh::conToStr(VecStr{"position", "insertSize", "counts"}, "\t")<< std::endl;

	for (const auto & key : allPositions) {
		if (njh::in(key, keys)) {
			for (const auto & inSize : pairCounts[key]) {
				outFile << key << "\t" << inSize.first << "\t" << inSize.second
						<< std::endl;
			}
		} else {
			outFile << key << "\t" << "NA" << "\t" << "NA"
					<< std::endl;
		}
	}
	return 0;
}



int bamExpRunner::testingBamToFastqToContigs(const njh::progutils::CmdArgs & inputCommands){
	SeqIOOptions reextractedSeqsOptsSingles;
	BamExtractor::extractReadsFromBamToSameOrientationContigsPars bamContigsExtractPars;
	bool keepThrownAwayRemappedUnMappedMate = false;
	reextractedSeqsOptsSingles.out_.outFilename_  = "out";
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processReadInNames({"--bam"}, true);
	setUp.processWritingOptions(reextractedSeqsOptsSingles.out_);
	setUp.setOption(keepThrownAwayRemappedUnMappedMate, "--keepThrownAwayRemappedUnMappedMate", "Keep Thrown Away Re-mapped Unmapped Mate");
	bamContigsExtractPars.throwAwayUnmmpaedMates = !keepThrownAwayRemappedUnMappedMate;
	setUp.setOption(bamContigsExtractPars.centerClipCutOff, "--centerClipCutOff", "Center Clip Cut Off");
	setUp.setOption(bamContigsExtractPars.forSoftClipFilterDistanceToEdges, "--forSoftClipFilterDistanceToEdges", "For Soft Clip Filter Distance To Edges");
	setUp.finishSetUp(std::cout);

	BamExtractor bExtractor(setUp.pars_.verbose_);
	bExtractor.debug_ = setUp.pars_.debug_;

	BamExtractor::ExtractedFilesOpts rextractedSeqsSingles;
	reextractedSeqsOptsSingles.firstName_ = setUp.pars_.ioOptions_.firstName_;
	rextractedSeqsSingles = bExtractor.extractReadsFromBamToSameOrientationContigs(reextractedSeqsOptsSingles, bamContigsExtractPars);
	rextractedSeqsSingles.log(std::cout, setUp.pars_.ioOptions_.firstName_);
	return 0;
}



int bamExpRunner::testingBamToFastq(const njh::progutils::CmdArgs & inputCommands){
	SeqIOOptions reextractedSeqsOptsSingles;
	BamExtractor::extractReadsFromBamToSameOrientationContigsPars bamContigsExtractPars;
	bool keepThrownAwayRemappedUnMappedMate = false;
	reextractedSeqsOptsSingles.out_.outFilename_  = "out";
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processReadInNames({"--bam"}, true);
	setUp.processWritingOptions(reextractedSeqsOptsSingles.out_);
	setUp.setOption(keepThrownAwayRemappedUnMappedMate, "--keepThrownAwayRemappedUnMappedMate", "Keep Thrown Away Re-mapped Unmapped Mate");
	bamContigsExtractPars.throwAwayUnmmpaedMates = !keepThrownAwayRemappedUnMappedMate;
	setUp.setOption(bamContigsExtractPars.centerClipCutOff, "--centerClipCutOff", "Center Clip Cut Off");
	setUp.setOption(bamContigsExtractPars.forSoftClipFilterDistanceToEdges, "--forSoftClipFilterDistanceToEdges", "For Soft Clip Filter Distance To Edges");
	setUp.finishSetUp(std::cout);

	BamExtractor bExtractor(setUp.pars_.verbose_);
	bExtractor.debug_ = setUp.pars_.debug_;

	BamExtractor::ExtractedFilesOpts rextractedSeqsSingles;
	reextractedSeqsOptsSingles.firstName_ = setUp.pars_.ioOptions_.firstName_;
	rextractedSeqsSingles = bExtractor.extractReadsFromBamToSameOrientationContigs(reextractedSeqsOptsSingles, bamContigsExtractPars);
	rextractedSeqsSingles.log(std::cout, setUp.pars_.ioOptions_.firstName_);

	return 0;
}



}



