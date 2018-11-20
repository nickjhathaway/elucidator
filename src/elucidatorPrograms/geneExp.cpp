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
	bfs::path inputFile = "";
	bfs::path twoBitFnp = "";
	OutOptions outOpts(bfs::path("out"));
	outOpts.outExtention_ = ".tab.txt";
	std::string idInput = "";
	seqSetUp setUp(inputCommands);
	setUp.setOption(inputFile, "--gff", "Input gff file", true);
	setUp.setOption(twoBitFnp, "--2bit", "Two bit file", true);
	setUp.setOption(idInput, "--id", "Id(s) to extract", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	auto idsVec = getInputValues(idInput, ",");
	std::set<std::string> ids(idsVec.begin(), idsVec.end());
	auto genes = GeneFromGffs::getGenesFromGffForIds(inputFile, ids);
	TwoBit::TwoBitFile tReader(twoBitFnp);

	std::string gffHeader = "";
	{
		std::stringstream gffHeaderStream;
		//write header
		std::ifstream infile(inputFile.string());
		std::string line = "";
		while('#' == infile.peek()){
			njh::files::crossPlatGetline(infile, line);
			gffHeaderStream << line << std::endl;
		}
		gffHeader = gffHeaderStream.str();
	}

	for(const auto & gene : genes){
		//std::cout << gene.first << "\t" << gene.second->getOneGeneDetailedName() << std::endl;
//		auto names = gene.second->getGeneDetailedName();
//		for(const auto & name : names){
//			std::cout << name.first << "\t" << name.second << std::endl;
//		}
		auto gsInfos = gene.second->generateGeneSeqInfo(tReader, false);
		//gff
		auto gffOpts = OutOptions(bfs::path(outOpts.outFilename_.string() + "_" + gene.second->gene_->getAttr("ID") + ".gff"));
		gffOpts.transferOverwriteOpts(outOpts);
		OutputStream gffOut(gffOpts);
		gffOut << gffHeader;
		gene.second->writeGffRecords(gffOut);

		for(const auto & transcript : gene.second->mRNAs_){

			GenomicRegion mRnaRegion(*transcript);
			OutOptions transcriptOut(bfs::path(outOpts.outFilename_.string() + "_" + transcript->getAttr("ID") + "_basePositions"), ".tab.txt");
			transcriptOut.transferOverwriteOpts(outOpts);
			auto tOutFile = transcriptOut.openFile();
			auto gsInfo = gsInfos[transcript->getIDAttr()];
			gsInfo->infoTab_.addColumn(VecStr{transcript->getIDAttr()}, "transcript");
			gsInfo->infoTab_.addColumn(VecStr{gene.second->gene_->getIDAttr()}, "GeneID");
			gsInfo->infoTab_.addColumn(VecStr{std::string(1, transcript->strand_)}, "strand");
			gsInfo->infoTab_.outPutContents(*tOutFile, "\t");
			auto gDNAOpts = SeqIOOptions::genFastaOut(    outOpts.outFilename_.string() + "_" + transcript->getAttr("ID") + "_gDNA");
			auto cDNAOpts = SeqIOOptions::genFastaOut(    outOpts.outFilename_.string() + "_" + transcript->getAttr("ID") + "_cDNA");
			auto proteinOpts = SeqIOOptions::genFastaOut( outOpts.outFilename_.string() + "_" + transcript->getAttr("ID") + "_protein");
			auto tableOpts = TableIOOpts::genTabFileOut(  outOpts.outFilename_.string() + "_" + transcript->getAttr("ID") + "_exonIntronPositions", true);
			auto bedOpts = OutOptions(bfs::path(          outOpts.outFilename_.string() + "_" + transcript->getAttr("ID") + "_exonIntronPositions.bed"));
			auto transcriptBedOpts = OutOptions(bfs::path(outOpts.outFilename_.string() + "_" + transcript->getAttr("ID") + ".bed"));

			bedOpts.transferOverwriteOpts(outOpts);
			transcriptBedOpts.transferOverwriteOpts(outOpts);
			tableOpts.out_.transferOverwriteOpts(outOpts);
			OutputStream transcriptBedOut(transcriptBedOpts);
			transcriptBedOut << GenomicRegion(*transcript).genBedRecordCore().toDelimStrWithExtra() <<std::endl;
			auto exonIntronPositions = gene.second->getIntronExonTables();
			auto exonIntronBeds = gene.second->getIntronExonBedLocs();
			exonIntronPositions[transcript->getIDAttr()].outPutContents(tableOpts);
			BioDataFileIO<Bed6RecordCore> reader{IoOptions(bedOpts)};
			reader.openWrite(exonIntronBeds[transcript->getIDAttr()], [](const Bed6RecordCore & record, std::ostream & out){
				out << record.toDelimStrWithExtra() << std::endl;
			});
			gDNAOpts.out_.transferOverwriteOpts(outOpts);
			cDNAOpts.out_.transferOverwriteOpts(outOpts);
			proteinOpts.out_.transferOverwriteOpts(outOpts);
			gsInfo->cDna_.name_ = transcript->getAttr("ID") + "_CodingDNA";
			gsInfo->gDna_.name_ = transcript->getAttr("ID") + "_GenomicDNA";
			SeqOutput::write(std::vector<seqInfo>{gsInfo->gDna_}, gDNAOpts);
			SeqOutput::write(std::vector<seqInfo>{gsInfo->cDna_}, cDNAOpts);
			gsInfo->cDna_.name_ = transcript->getAttr("ID") + "_protein";
			gsInfo->cDna_.translate(false, false);
			if('*' == gsInfo->cDna_.seq_.back()){
				gsInfo->cDna_.trimBack(gsInfo->cDna_.seq_.size() - 1);
			}
			SeqOutput::write(std::vector<seqInfo>{gsInfo->cDna_}, proteinOpts);


		}
	}
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
