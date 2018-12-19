/*
 * geneExpRunner.cpp
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

#include "geneExp.hpp"
#include <njhseq/objects/Gene/GeneFromGffs.hpp>
#include "elucidator/objects/BioDataObject.h"
#include "elucidator/seqToolsUtils/seqToolsUtils.hpp"

#include <TwoBit.h>

namespace njhseq {
geneExpRunner::geneExpRunner()
    : njh::progutils::ProgramRunner(
          {
					 addFunc("cDNAPosTogDNAPos", cDNAPosTogDNAPos, false),
					 addFunc("getBedOfAminoAcidPositions", getBedOfAminoAcidPositions, false),
					 addFunc("getBedOfCDnaPositions", getBedOfCDnaPositions, false),
					 addFunc("gffRecordIDToGeneInfo", gffRecordIDToGeneInfo, false),
           },
          "geneExp") {}
//




int geneExpRunner::gffRecordIDToGeneInfo(const njh::progutils::CmdArgs & inputCommands){
	GeneFromGffs::gffRecordIDsToGeneInfoPars pars;
	std::string idInput = "";
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.setOption(pars.inputFile, "--gff", "Input gff file", true);
	setUp.setOption(pars.twoBitFnp, "--2bit", "Two bit file", true);
	setUp.setOption(idInput, "--id", "Id(s) to extract", true);
	setUp.processWritingOptions(pars.outOpts);
	setUp.processDirectoryOutputName("./", false);
	setUp.finishSetUp(std::cout);

	pars.outOpts.outFilename_ = njh::files::make_path(setUp.pars_.directoryName_, pars.outOpts.outFilename_);

	auto idsVec = getInputValues(idInput, ",");
	pars.ids = std::set<std::string>(idsVec.begin(), idsVec.end());

	if(setUp.pars_.verbose_){
		std::cout << "Read in " << pars.ids.size() << " ids, IDs: " << njh::conToStr(pars.ids, ", ") << std::endl;
	}

	GeneFromGffs::gffRecordIDsToGeneInfo(pars);

	return 0;
}

int geneExpRunner::cDNAPosTogDNAPos(const njh::progutils::CmdArgs & inputCommands) {
	seqInfo cDna;
	seqInfo gDna;
	bfs::path bedFile = "";
	GeneSeqInfo::GeneSeqInfoPars pars;
	auto outOpts = TableIOOpts::genTabFileOut("out.tab.txt", true);
	seqSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processVerbose();
	setUp.processWritingOptions(outOpts.out_);
	setUp.processSeq(cDna, "--cDNA",
			"cDNA fasta file where the first record is the cDNA", true);
	setUp.setOption(bedFile, "--bedFile",
				"A bed file with the location of the gene, zero based", true);
	setUp.setOption(pars.oneBasedPos_, "--oneBasedPos",
				"Output 1 based positioning");
	setUp.processSeq(gDna, "--gDna",
			"gDNA fasta file where the first record is the gDNA (introns should be in lower case)", true);
	setUp.processGap();
	setUp.processScoringPars();
	setUp.finishSetUp(std::cout);
	uint64_t maxLen = 0;
	readVec::getMaxLength(cDna, maxLen);
	readVec::getMaxLength(gDna, maxLen);
	aligner alignerObj(maxLen, setUp.pars_.gapInfo_, setUp.pars_.scoring_,
			setUp.pars_.colOpts_.alignOpts_.countEndGaps_);

	auto regions = gatherRegions(bedFile.string(), "", setUp.pars_.verbose_);
	pars.region_ = regions.front();
	GeneSeqInfo info(cDna, gDna, pars);
	info.setCDnaAlnByAligning(alignerObj);
	info.setTable();
	if(setUp.pars_.debug_){
		alignerObj.alignObjectA_.seqBase_.outPutSeq(std::cout);
		alignerObj.alignObjectB_.seqBase_.outPutSeq(std::cout);
	}
	info.infoTab_.outPutContents(outOpts);
	return 0;
}

int geneExpRunner::getBedOfAminoAcidPositions(const njh::progutils::CmdArgs & inputCommands) {
	seqInfo cDna;
	seqInfo gDna;
	bfs::path bedFile = "";
	uint32_t aaStart = 0;
	uint32_t aaStop = 0;
	GeneSeqInfo::GeneSeqInfoPars pars;
	OutOptions outOpts (bfs::path("out.bed"));
	seqSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processVerbose();
	setUp.processWritingOptions(outOpts);
	setUp.setOption(aaStart, "--aaStart",
			"aa Start position", true);
	setUp.setOption(aaStop, "--aaStop",
			"aa Stop position", true);
	setUp.processSeq(cDna, "--cDNA",
			"cDNA fasta file where the first record is the cDNA", true);
	setUp.setOption(bedFile, "--bedFile",
				"A bed file with the location of the gene, zero based", true);
	setUp.setOption(pars.oneBasedPos_, "--oneBasedPos",
				"Output 1 based positioning");
	setUp.processSeq(gDna, "--gDna",
			"gDNA fasta file where the first record is the gDNA (introns should be in lower case)", true);
	setUp.processGap();
	setUp.processScoringPars();
	setUp.finishSetUp(std::cout);
	uint64_t maxLen = 0;
	readVec::getMaxLength(cDna, maxLen);
	readVec::getMaxLength(gDna, maxLen);
	aligner alignerObj(maxLen, setUp.pars_.gapInfo_, setUp.pars_.scoring_,
			setUp.pars_.colOpts_.alignOpts_.countEndGaps_);

	auto regions = gatherRegions(bedFile.string(), "", setUp.pars_.verbose_);
	pars.region_ = regions.front();
	GeneSeqInfo info(cDna, gDna, pars);
	info.setCDnaAlnByAligning(alignerObj);
	info.setTable();
	if(setUp.pars_.debug_){
		alignerObj.alignObjectA_.seqBase_.outPutSeq(std::cout);
		alignerObj.alignObjectB_.seqBase_.outPutSeq(std::cout);
	}

	std::ofstream outBed;
	std::ostream out(determineOutBuf(outBed, outOpts));

	out << info.genBedFromAAPositions(aaStart, aaStop).toDelimStr() << std::endl;

	return 0;
}

int geneExpRunner::getBedOfCDnaPositions(const njh::progutils::CmdArgs & inputCommands) {
	seqInfo cDna;
	seqInfo gDna;
	bfs::path bedFile = "";
	uint32_t start = 0;
	uint32_t stop = 0;
	GeneSeqInfo::GeneSeqInfoPars pars;
	OutOptions outOpts (bfs::path("out.bed"));
	seqSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processVerbose();
	setUp.processWritingOptions(outOpts);
	setUp.setOption(start, "--start",
			"aa Start position", true);
	setUp.setOption(stop, "--stop",
			"aa Stop position", true);
	setUp.processSeq(cDna, "--cDNA",
			"cDNA fasta file where the first record is the cDNA", true);
	setUp.setOption(bedFile, "--bedFile",
				"A bed file with the location of the gene, zero based", true);
	setUp.setOption(pars.oneBasedPos_, "--oneBasedPos",
				"Output 1 based positioning");
	setUp.processSeq(gDna, "--gDna",
			"gDNA fasta file where the first record is the gDNA (introns should be in lower case)", true);
	setUp.processGap();
	setUp.processScoringPars();
	setUp.finishSetUp(std::cout);
	uint64_t maxLen = 0;
	readVec::getMaxLength(cDna, maxLen);
	readVec::getMaxLength(gDna, maxLen);
	aligner alignerObj(maxLen, setUp.pars_.gapInfo_, setUp.pars_.scoring_,
			setUp.pars_.colOpts_.alignOpts_.countEndGaps_);

	auto regions = gatherRegions(bedFile.string(), "", setUp.pars_.verbose_);
	pars.region_ = regions.front();
	GeneSeqInfo info(cDna, gDna, pars);
	info.setCDnaAlnByAligning(alignerObj);
	info.setTable();
	if(setUp.pars_.debug_){
		alignerObj.alignObjectA_.seqBase_.outPutSeq(std::cout);
		alignerObj.alignObjectB_.seqBase_.outPutSeq(std::cout);
	}

	std::ofstream outBed;
	std::ostream out(determineOutBuf(outBed, outOpts));

	out << info.genBedFromCDNAPositions(start, stop).toDelimStr() << std::endl;

	return 0;
}



} /* namespace njhseq */
