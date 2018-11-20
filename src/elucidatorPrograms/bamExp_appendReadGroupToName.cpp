/*
 * bamExp_appendReadGroupToName.cpp
 *
 *  Created on: Nov 10, 2017
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


namespace njhseq {


int bamExpRunner::BamRenameRefHeader(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path otherBam;
	OutOptions outOpts(bfs::path(""));
	outOpts.outExtention_ = ".bam";
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processReadInNames( { "--bam" }, true);
	setUp.setOption(otherBam, "--otherBam", "The other bam to rename the reference to", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	BamTools::BamReader bReader;
	bReader.Open(setUp.pars_.ioOptions_.firstName_.string());
	checkBamOpenThrow(bReader, setUp.pars_.ioOptions_.firstName_.string());

	BamTools::BamReader bReaderOther;
	bReaderOther.Open(otherBam.string());
	checkBamOpenThrow(bReaderOther, otherBam.string());

	BamTools::BamAlignment bAln;
	if("" == outOpts.outFilename_){
		outOpts.outFilename_ = njh::files::prependFileBasename(setUp.pars_.ioOptions_.firstName_, "renamed_");
	}
	if(outOpts.outExists() && !outOpts.overWriteFile_){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << outOpts.outName() << " already exists, use --overWrite to over write it" << "\n";
		throw std::runtime_error{ss.str()};
	}

	//check to see the two contain the same reference data, just in different order, if same order do nothing

	auto otherRefData = bReaderOther.GetReferenceData();
	auto currentRefData = bReader.GetReferenceData();

	std::unordered_map<std::string, int32_t> otherRefIdKey;
	std::unordered_map<int32_t, std::string> currentRefIdReverseKey;
	if (otherRefData.size() == currentRefData.size()) {
		//check to see if same header
		int32_t numberOfHits = 0;
		for (auto sIter = bReaderOther.GetConstSamHeader().Sequences.Begin();
				sIter != bReaderOther.GetConstSamHeader().Sequences.End(); ++sIter) {
			for (auto currentSIter = bReader.GetConstSamHeader().Sequences.Begin();
					currentSIter != bReader.GetConstSamHeader().Sequences.End(); ++currentSIter) {
				if((*sIter) == (*currentSIter)){
					++numberOfHits;
					break;
				}
			}
		}
		if(numberOfHits != bReaderOther.GetConstSamHeader().Sequences.Size()){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error was expecting to get " << bReaderOther.GetConstSamHeader().Sequences.Size() << " hits while checking references data match up but got " << numberOfHits << " instead" << "\n";
			throw std::runtime_error{ss.str()};
		}
		bool sameRefDataOrder = true;
		for(const auto & pos : iter::range(otherRefData.size())){
			//just need to check names as the above check should make sure the rest (length etc is the same)
			if(currentRefData[pos].RefName != otherRefData[pos].RefName){
				sameRefDataOrder = false;
				break;
			}
		}
		if (sameRefDataOrder) {
			std::cout
					<< "References are already in the same order, no need to re organize "
					<< setUp.pars_.ioOptions_.firstName_ << " to match " << otherBam
					<< std::endl;
			return 0;
		}
		for(const auto pos : iter::range(otherRefData.size())){
			otherRefIdKey[otherRefData[pos].RefName] = pos;
		}
		otherRefIdKey["*"] = -1;

		for(const auto pos : iter::range(currentRefData.size())){
			currentRefIdReverseKey[pos] = currentRefData[pos].RefName;
		}
		currentRefIdReverseKey[-1] = "*";
	} else {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__
				<< ", error reference data sizes are different, can't be same reference sequences"
				<< "\n";
		ss << otherBam << " refSize: " << otherRefData.size() << "\n";
		ss << setUp.pars_.ioOptions_.firstName_ << " refSize: " << currentRefData.size() << '\n';
		throw std::runtime_error { ss.str() };
	}

	BamTools::BamWriter bWriter;
	auto outHeader = bReader.GetHeader();
	outHeader.Sequences = bReaderOther.GetHeader().Sequences;
	bWriter.Open(outOpts.outName().string(), outHeader, otherRefData);

	while(bReader.GetNextAlignment(bAln)){
		bAln.RefID = otherRefIdKey[currentRefIdReverseKey[bAln.RefID]];
		bAln.MateRefID = otherRefIdKey[currentRefIdReverseKey[bAln.MateRefID]];
		bWriter.SaveAlignment(bAln);
	}

	return 0;
}


int bamExpRunner::appendReadGroupToName(
		const njh::progutils::CmdArgs & inputCommands) {
	OutOptions outOpts(bfs::path(""));

	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processReadInNames( { "--bam" }, true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	BamTools::BamReader bReader;
	bReader.Open(setUp.pars_.ioOptions_.firstName_.string());
	checkBamOpenThrow(bReader, setUp.pars_.ioOptions_.firstName_.string());
	BamTools::BamAlignment bAln;
	if("" == outOpts.outFilename_){
		outOpts.outFilename_ = njh::files::prependFileBasename(setUp.pars_.ioOptions_.firstName_, "renamed_");
	}
	if(outOpts.outExists() && !outOpts.overWriteFile_){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << outOpts.outName() << " already exists, use --overWrite to over write it" << "\n";
		throw std::runtime_error{ss.str()};
	}
	BamTools::BamWriter bWriter;
	bWriter.Open(outOpts.outName().string(), bReader.GetHeader(), bReader.GetReferenceData());
	auto header = bReader.GetHeader();
	std::unordered_map<std::string,std::string> readGroupKeys;
	if(header.ReadGroups.IsEmpty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << outOpts.outName() << " header doesn't contain any read group(RG) information" << "\n";
		throw std::runtime_error{ss.str()};
	}
	for(const auto & it : header.ReadGroups){
		readGroupKeys[it.ID] = it.Sample;
	}
	while(bReader.GetNextAlignment(bAln)){
		std::string rg = "";
		bool success = bAln.GetTag("RG", rg);
		if(!success){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error to get read group(RG) information for " << bAln.Name << "\n";
			throw std::runtime_error{ss.str()};
		}
		if(!njh::in(rg, readGroupKeys)){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error read group(RG) " << rg << " not found in header" << "\n";
			ss << "Options are: " << njh::conToStr(njh::getVecOfMapKeys(readGroupKeys), ", ") << "\n";
			throw std::runtime_error{ss.str()};
		}
		bAln.Name.append("[sample=" + readGroupKeys[rg]+"]");
		bWriter.SaveAlignment(bAln);
	}

	return 0;
}


} //namespace njhseq
