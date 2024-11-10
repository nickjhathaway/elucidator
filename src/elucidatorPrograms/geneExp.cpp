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
#include <njhseq/GenomeUtils/GenomeMapping/MultiGenomeMapper.hpp>
#include <njhseq/objects/helperObjects/AminoAcidPositionInfo.hpp>

namespace njhseq {
geneExpRunner::geneExpRunner()
    : njh::progutils::ProgramRunner(
          {
					 addFunc("cDNAPosTogDNAPos", cDNAPosTogDNAPos, false),
					 addFunc("getBedOfAminoAcidPositions", getBedOfAminoAcidPositions, false),
					 addFunc("getBedOfCDnaPositions", getBedOfCDnaPositions, false),
					 addFunc("gffRecordIDToGeneInfo", gffRecordIDToGeneInfo, false),
          	addFunc("getBedOfAminoAcidPositionsFromGff", getBedOfAminoAcidPositionsFromGff, false),
          	addFunc("multiGenomeExtractGenesWithDescription", multiGenomeExtractGenesWithDescription, false),
          	addFunc("bedGetOverlappingAminoAcidPositions", bedGetOverlappingAminoAcidPositions, false),
           },
          "geneExp") {}
//


template <typename BEDREC>
std::unordered_map<uint32_t, std::vector<GFFCore>> getGffRecordsIntersectWithBeds(
		std::vector<BEDREC> & beds,
		const bfs::path & gffFnp,
		const VecStr & selectFeatures = VecStr{"gene", "protein_coding_gene"}){
	std::unordered_map<uint32_t, std::vector<GFFCore>> ret;

	std::unordered_map<std::string, std::vector<uint32_t>> bedsByChrome;

	BioDataFileIO<GFFCore> reader { IoOptions(InOptions(gffFnp)) };
	reader.openIn();
	// uint32_t count = 0;
	std::string line;
	std::shared_ptr<GFFCore> gRecord = reader.readNextRecord();
	for (const auto bPos : iter::range(beds.size())) {
		bedsByChrome[getRef(beds[bPos]).chrom_].emplace_back(bPos);
	}

	while (nullptr != gRecord) {
		if (selectFeatures.empty() || njh::in(gRecord->type_, selectFeatures)) {
			for (auto & inputRegionPos : bedsByChrome[gRecord->seqid_]) {
				if (Bed3RecordCore::getOverlapLen(getRef(beds[inputRegionPos]).chrom_,
					getRef(beds[inputRegionPos]).chromStart_,
					getRef(beds[inputRegionPos]).chromEnd_,
					gRecord->seqid_,
					gRecord->start_ - 1,
					gRecord->end_) >=1){
					ret[inputRegionPos].emplace_back(*gRecord);
				}
			}
		}
		bool end = false;
		while ('#' == reader.inFile_->peek()) {
			if (njh::files::nextLineBeginsWith(*reader.inFile_, "##FASTA")) {
				end = true;
				break;
			}
			njh::files::crossPlatGetline(*reader.inFile_, line);
		}
		if (end) {
			break;
		}
		gRecord = reader.readNextRecord();
		// ++count;
	}
	return ret;
}


template<typename BEDREC>
std::unordered_map<uint32_t, std::vector<MultiGenomeMapper::IntersectedProteinInfo>>
addIntersectingGeneInfosToLocs(
	std::vector<BEDREC>&regions,
	const intersectBedLocsWtihGffRecordsPars & pars,
	const bfs::path & twoBitFnp) {
	std::unordered_map<uint32_t, std::vector<MultiGenomeMapper::IntersectedProteinInfo>> ret;
	auto gffRecords = getGffRecordsIntersectWithBeds(regions, pars.gffFnp_, pars.selectFeatures_);
	std::set<std::string> rawReneIdsSet;
	for (const auto & gRecForBed : gffRecords) {
		for (const auto & g : gRecForBed.second) {
			rawReneIdsSet.insert(g.getIDAttr());
		}
	}

	auto rawGenes = GeneFromGffs::getGenesFromGffForIds(pars.gffFnp_, rawReneIdsSet);
	std::unordered_map<std::string, std::shared_ptr<GeneFromGffs>> genes;
	for (const auto&gene: rawGenes) {
		bool failFilter = false;
		for (const auto&transcript: gene.second->mRNAs_) {
			if (njh::in(njh::strToLowerRet(transcript->type_), VecStr{"rrna", "trna", "snorna", "snrna", "ncrna"})) {
				failFilter = true;
				break;
			}
		}
		if (!failFilter) {
			genes[gene.first] = gene.second;
		}
	}
	for (const auto&regPos: iter::range(regions.size())) {
		if (njh::notIn(regPos, gffRecords)) {
			continue;
		}

		for (const auto & gffRec: gffRecords[regPos]) {
			if (njh::in(gffRec.getIDAttr(), genes)) {
				TwoBit::TwoBitFile tReader(twoBitFnp);
				auto infos = njh::mapAt(genes, gffRec.getIDAttr())->generateGeneSeqInfo(tReader, false);
				auto detailedName = njh::mapAt(genes, gffRec.getIDAttr())->getGeneDetailedName();
				for (const auto&info: infos) {

					auto posInfos = info.second->getInfosByGDNAPos();
					auto minPos = vectorMinimum(getVectorOfMapKeys(posInfos));
					auto maxPos = vectorMaximum(getVectorOfMapKeys(posInfos));
					auto startPos = std::max(minPos, getRef(regions[regPos]).chromStart_);
					auto stopPos = std::min(maxPos, getRef(regions[regPos]).chromEnd_);
					if (stopPos == getRef(regions[regPos]).chromEnd_) {
						stopPos = getRef(regions[regPos]).chromEnd_ - 1;
					}
					//stop position is inclusive
					uint32_t aaStartPos = std::numeric_limits<uint32_t>::max();
					uint32_t aaStopPos = 0;

					for (const auto posInGene: iter::range(startPos, stopPos + 1)) {
						if (std::numeric_limits<uint32_t>::max() != posInfos[posInGene].aaPos_) {
							if (posInfos[posInGene].aaPos_ < aaStartPos) {
								aaStartPos = posInfos[posInGene].aaPos_;
							}
							if (posInfos[posInGene].aaPos_ > aaStopPos) {
								aaStopPos = posInfos[posInGene].aaPos_;
							}
						}
					}
					//if the region is completely within the intron the final stop and stops will be max() and 0;
					if (std::numeric_limits<uint32_t>::max() == aaStartPos) {
						aaStopPos = std::numeric_limits<uint32_t>::max();
					} else {
						//make 1 based
						++aaStopPos;
						++aaStartPos; //
					}
					ret[regPos].emplace_back(info.second->transcriptID_, aaStartPos, aaStopPos, detailedName[info.first]);
					ret[regPos].back().allMeta_.meta_ = gffRec.attributes_;
				}
			}
		}
	}
	return ret;
}

int geneExpRunner::bedGetOverlappingAminoAcidPositions(const njh::progutils::CmdArgs & inputCommands){
	intersectBedLocsWtihGffRecordsPars pars;
	pars.selectFeatures_ = VecStr{"gene", "protein_coding_gene"};
	pars.extraAttributes_ = VecStr{"ID"};
	bfs::path genomeTwoBit;
	bfs::path input;
	OutOptions outOpts(bfs::path(""), ".bed");
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.setOption(pars.gffFnp_, "--gff", "Input gff file", true);
	setUp.setOption(genomeTwoBit, "--genomeTwoBit", "genome Two Bit", true);
	setUp.setOption(pars.selectFeatures_, "--selectFeatures", "select gff Features");
	setUp.setOption(pars.extraAttributes_, "--extraAttributes", "extra gff Attributes to add to output");


	setUp.setOption(input, "--bed", "input bed", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	auto beds = getBeds(input);
	uint32_t extraFieldsMaxCount = 0;
	for (const auto & b : beds) {
		if (b->extraFields_.size() > extraFieldsMaxCount) {
			extraFieldsMaxCount = b->extraFields_.size();
		}
	}

	OutputStream out(outOpts);
	out << "#chrom\tstart\tend\tname\tscore\tstrand";
	for (auto pos = 0U; pos < extraFieldsMaxCount; ++pos) {
		out << "\t" << "extraField." << njh::leftPadNumStr(pos, extraFieldsMaxCount) ;
	}
	out << "\ttranscriptID\tAAStart-1based\tAAEnd-1based\tGeneDescription";
	out << "\t" << njh::conToStr(pars.extraAttributes_, "\t");
	out << std::endl;

	auto results = addIntersectingGeneInfosToLocs(beds, pars, genomeTwoBit);
	for (const auto bedPos: iter::range(beds.size())) {
		if (njh::in(bedPos, results)) {
			for (const auto & result: results[bedPos]) {
				out << beds[bedPos]->toDelimStr();
				for (const auto& extra: beds[bedPos]->extraFields_) {
					out << "\t" << extra;
				}
				for (auto pos = 0U; pos < extraFieldsMaxCount - beds[bedPos]->extraFields_.size(); ++pos) {
					out << "\tNA";
				}
				out << "\t" << result.id_
				<< "\t" << result.aaStart_
				<< "\t" << result.aaStop_
				<< "\t" << result.description_;
				for (const auto & field : pars.extraAttributes_) {
					if (result.allMeta_.containsMeta(field)) {
						out << "\t" << result.allMeta_.getMeta(field);
					} else {
						out << "\tNA";
					}
				}
				out << std::endl;
			}
		} else {
			out << beds[bedPos]->toDelimStr();
			for (const auto& extra: beds[bedPos]->extraFields_) {
				out << "\t" << extra;
			}
			for (auto pos = 0U; pos < extraFieldsMaxCount - beds[bedPos]->extraFields_.size(); ++pos) {
				out << "\tNA";
			}
			for (auto emptyPos = 0; emptyPos < 4 + pars.extraAttributes_.size(); ++emptyPos) {
				out << "\tNA";
			}
			out << std::endl;
		}
	}
	return 0;
}

int geneExpRunner::gffRecordIDToGeneInfo(const njh::progutils::CmdArgs & inputCommands){
	GeneFromGffs::gffRecordIDsToGeneInfoPars pars;
	std::string idInput;
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
	if("./" != setUp.pars_.directoryName_){
		setUp.startARunLog(setUp.pars_.directoryName_);
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

int geneExpRunner::getBedOfAminoAcidPositionsFromGff(const njh::progutils::CmdArgs & inputCommands) {

	std::string id;
	bfs::path gffFnp;
	bfs::path twoBitFnp;
	std::string addMetaField;
	bfs::path genomeFnp;
	uint32_t input_aaStart = 0;
	uint32_t input_aaStop = 0;
	std::set<uint32_t> input_aaPositions;
	bool aaPositionsZeroBased = false;
	GeneSeqInfo::GeneSeqInfoPars pars;
	bool collapsePerId = false;
	OutOptions outOpts (bfs::path(""), ".bed");
	seqSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processVerbose();
	setUp.processWritingOptions(outOpts);
	setUp.setOption(addMetaField, "--addMetaField", "Add Meta Field");

	setUp.setOption(id, "--geneID", "gene or transcript ID", true);
	setUp.setOption(gffFnp, "--gff", "gff", true);
	setUp.setOption(twoBitFnp, "--2bit", "2Bit Fnp", true);
	bool aaPositionsIsSet = setUp.setOption(input_aaPositions, "--aaPositions", "aa Positions", false);
	setUp.setOption(input_aaStart, "--aaStart", "aa Start position", !aaPositionsIsSet);
	setUp.setOption(input_aaStop, "--aaStop","aa Stop position", !aaPositionsIsSet);
	setUp.setOption(collapsePerId, "--collapsePerId","collapse Per Id");
	setUp.setOption(pars.oneBasedPos_, "--oneBasedPos","Output 1 based positioning");
	setUp.setOption(aaPositionsZeroBased, "--aaPositionsZeroBased","Input amino acid positions zero based");
	setUp.finishSetUp(std::cout);

	OutputStream out(outOpts);
	table aminoAcidPositionsTable(input_aaPositions.empty()? VecStr{"id", "aaStart", "aaStop"}: VecStr{"id", "aaPosition"});
	if(input_aaPositions.empty()) {
		aminoAcidPositionsTable.addRow(id, input_aaStart, input_aaStop);
	} else {
		for(const auto & pos : input_aaPositions) {
			aminoAcidPositionsTable.addRow(id, pos);
		}
	}
	AminoAcidPositionInfo aaInfos(aminoAcidPositionsTable, aaPositionsZeroBased);

	auto genes = GeneFromGffs::getGenesFromGffForGuessedTranscriptOrGeneIds(gffFnp, aaInfos.ids_);

	TwoBit::TwoBitFile tReader(twoBitFnp);

	if(aaInfos.byRange()){

		for(const auto & row : aaInfos.infoTab_){
			std::string idColName = "id";
			if(njh::in(std::string("transcriptid"), aaInfos.infoTab_.columnNames_)){
				idColName = "transcriptid";
			}
			bool byTranscript = false;
			auto idName = row[aaInfos.infoTab_.getColPos(idColName)];
			std::string geneID = idName;
			for(const auto & gene : genes){
				for(const auto & mRNA : gene.second->mRNAs_){
					if(mRNA->getIDAttr() == idName){
						byTranscript = true;
						geneID = gene.first;
						break;
					}
				}
			}

			auto gsInfos = njh::mapAt(genes, geneID)->generateGeneSeqInfo(tReader, false);
			MetaDataInName metaForCollapse;
			auto aaStart =
							aaInfos.zeroBased_ ?
							njh::StrToNumConverter::stoToNum<uint32_t>(
											row[aaInfos.infoTab_.getColPos("aastart")]) :
							njh::StrToNumConverter::stoToNum<uint32_t>(
											row[aaInfos.infoTab_.getColPos("aastart")]) - 1;
			auto aastop =
							aaInfos.zeroBased_ ?
							njh::StrToNumConverter::stoToNum<uint32_t>(
											row[aaInfos.infoTab_.getColPos("aastop")]) :
							njh::StrToNumConverter::stoToNum<uint32_t>(
											row[aaInfos.infoTab_.getColPos("aastop")]);

			std::vector<uint32_t> aaPositions(aastop - aaStart);
			njh::iota(aaPositions, aaStart);

			std::unordered_map<std::string, VecStr> combinedMeta;
			std::unordered_map<std::string, std::unordered_set<std::string>> combinedMetaCollapsed;
			for(const auto & pos : aaPositions){
				for(const auto & metaFields : aaInfos.metaDataForAAPos_[idName][pos].meta_){
					combinedMetaCollapsed[metaFields.first].emplace(metaFields.second);
					combinedMeta[metaFields.first].emplace_back(njh::pasteAsStr(pos, "-", metaFields.second));
				}
			}
			for(const auto & field : combinedMeta){
				if(combinedMetaCollapsed[field.first].size() == 1){
					std::string metaField = *combinedMetaCollapsed[field.first].begin();
					metaForCollapse.addMeta(field.first,metaField);
				}else{
					metaForCollapse.addMeta(field.first, njh::conToStr(field.second, ","));
				}
			}
			metaForCollapse.addMeta("GeneID", geneID);
			if(byTranscript){
				auto gsInfo = njh::mapAt(gsInfos, idName);
				metaForCollapse.addMeta("transcript", idName);
				std::vector<uint32_t> posVec(aaPositions.begin(),
																		 aaPositions.end());
				auto minAAPos = vectorMinimum(posVec);
				auto maxAAPos = vectorMaximum(posVec);
				auto posBed = gsInfo->genBedFromAAPositions(minAAPos, maxAAPos + 1);
				posBed.extraFields_.emplace_back(metaForCollapse.createMetaName());
				if (aaPositionsZeroBased) {
					posBed.name_ = njh::pasteAsStr(idName, "-", "[AA", minAAPos,
																				 "-", maxAAPos + 1, ")");
				} else {
					posBed.name_ = njh::pasteAsStr(idName, "-", "[AA", minAAPos + 1, "-", maxAAPos + 1, "]");
				}
				out << posBed.toDelimStrWithExtra() << std::endl;
			}else{
				for(const auto & gsInfo : gsInfos){
					std::vector<uint32_t> posVec(aaPositions.begin(),
																			 aaPositions.end());
					auto minAAPos = vectorMinimum(posVec);
					auto maxAAPos = vectorMaximum(posVec);
					metaForCollapse.addMeta("transcript", gsInfo.first);
					auto posBed = gsInfo.second->genBedFromAAPositions(minAAPos, maxAAPos + 1);
					posBed.extraFields_.emplace_back(metaForCollapse.createMetaName());

					if (aaPositionsZeroBased) {
						posBed.name_ = njh::pasteAsStr(gsInfo.first, "-", "[AA", minAAPos,
																					 "-", maxAAPos + 1, ")");
					} else {
						posBed.name_ = njh::pasteAsStr(gsInfo.first, "-", "[AA", minAAPos + 1, "-", maxAAPos + 1, "]");
					}
					out << posBed.toDelimStrWithExtra() << std::endl;
				}
			}
		}

		for (const auto & positions : aaInfos.aminoPositionsPerId_) {
			// bool byTranscript = false;
			std::string geneID = positions.first;
			for(const auto & gene : genes){
				for(const auto & mRNA : gene.second->mRNAs_){
					if(mRNA->getIDAttr() == positions.first){
						// byTranscript = true;
						geneID = gene.first;
						break;
					}
				}
			}
		}
	} else {
		for (const auto & positions : aaInfos.aminoPositionsPerId_) {
			bool byTranscript = false;
			std::string geneID = positions.first;
			for(const auto & gene : genes){
				for(const auto & mRNA : gene.second->mRNAs_){
					if(mRNA->getIDAttr() == positions.first){
						byTranscript = true;
						geneID = gene.first;
						break;
					}
				}
			}

			auto gsInfos = njh::mapAt(genes, geneID)->generateGeneSeqInfo(tReader, false);
			MetaDataInName metaForCollapse;
			if (collapsePerId && positions.second.size() > 1) {
				std::unordered_map<std::string, VecStr> combinedMeta;
				std::unordered_map<std::string, std::unordered_set<std::string>> combinedMetaCollapsed;
				for(const auto & pos : positions.second){
					for(const auto & metaFields : aaInfos.metaDataForAAPos_[positions.first][pos].meta_){
						combinedMetaCollapsed[metaFields.first].emplace(metaFields.second);
						combinedMeta[metaFields.first].emplace_back(njh::pasteAsStr(pos, "-", metaFields.second));
					}
				}
				for(const auto & field : combinedMeta){
					if(combinedMetaCollapsed[field.first].size() == 1){
						std::string metaField = *combinedMetaCollapsed[field.first].begin();
						metaForCollapse.addMeta(field.first,metaField);
					}else{
						metaForCollapse.addMeta(field.first, njh::conToStr(field.second, ","));
					}
				}
				metaForCollapse.addMeta("GeneID", geneID);
			}
			if(byTranscript){
				auto gsInfo = njh::mapAt(gsInfos, positions.first);
				if (collapsePerId && positions.second.size() > 1) {
					metaForCollapse.addMeta("transcript", positions.first);
					std::vector<uint32_t> posVec(positions.second.begin(),
																			 positions.second.end());
					auto minAAPos = vectorMinimum(posVec);
					auto maxAAPos = vectorMaximum(posVec);
					auto posBed = gsInfo->genBedFromAAPositions(minAAPos, maxAAPos + 1);
					posBed.extraFields_.emplace_back(metaForCollapse.createMetaName());
					if (aaPositionsZeroBased) {
						posBed.name_ = njh::pasteAsStr(positions.first, "-", "[AA", minAAPos,
																					 "-", maxAAPos + 1, ")");
					} else {
						posBed.name_ = njh::pasteAsStr(positions.first, "-", "[AA", minAAPos + 1, "-", maxAAPos + 1, "]");
					}
					out << posBed.toDelimStrWithExtra() << std::endl;
				} else {
					for (const auto & pos : positions.second) {
						MetaDataInName meta = aaInfos.metaDataForAAPos_[positions.first][pos];
						meta.addMeta("transcript", positions.first);
						meta.addMeta("GeneID", geneID);

						auto posBed = gsInfo->genBedFromAAPositions(pos, pos + 1);
						posBed.extraFields_.emplace_back(meta.createMetaName());
						std::string add;
						if(!addMetaField.empty() && meta.containsMeta(addMetaField)){
							add = "-" + meta.getMeta(addMetaField);
						}
						if (aaPositionsZeroBased) {
							posBed.name_ = njh::pasteAsStr(positions.first, add, "-", "AA", pos);
						} else {
							posBed.name_ = njh::pasteAsStr(positions.first, add, "-", "AA", pos + 1);
						}
						out << posBed.toDelimStrWithExtra() << std::endl;
					}
				}
			}else{
				for(const auto & gsInfo : gsInfos){
					if (collapsePerId && positions.second.size() > 1) {
						std::vector<uint32_t> posVec(positions.second.begin(),
																				 positions.second.end());
						auto minAAPos = vectorMinimum(posVec);
						auto maxAAPos = vectorMaximum(posVec);
						metaForCollapse.addMeta("transcript", gsInfo.first);
						auto posBed = gsInfo.second->genBedFromAAPositions(minAAPos, maxAAPos + 1);
						posBed.extraFields_.emplace_back(metaForCollapse.createMetaName());

						if (aaPositionsZeroBased) {
							posBed.name_ = njh::pasteAsStr(gsInfo.first, "-", "[AA", minAAPos,
																						 "-", maxAAPos + 1, ")");
						} else {
							posBed.name_ = njh::pasteAsStr(gsInfo.first, "-", "[AA", minAAPos + 1, "-", maxAAPos + 1, "]");
						}
						out << posBed.toDelimStrWithExtra() << std::endl;
					} else {
						for (const auto & pos : positions.second) {
							MetaDataInName meta = aaInfos.metaDataForAAPos_[positions.first][pos];
							meta.addMeta("transcript", gsInfo.first);
							meta.addMeta("GeneID", geneID);

							auto posBed = gsInfo.second->genBedFromAAPositions(pos, pos + 1);
							posBed.extraFields_.emplace_back(meta.createMetaName());
							std::string add;
							if(!addMetaField.empty() && meta.containsMeta(addMetaField)){
								add = "-" + meta.getMeta(addMetaField);
							}
							if (aaPositionsZeroBased) {
								posBed.name_ = njh::pasteAsStr(gsInfo.first, add, "-", "AA", pos);
							} else {
								posBed.name_ = njh::pasteAsStr(gsInfo.first, add, "-", "AA", pos + 1);
							}
							out << posBed.toDelimStrWithExtra() << std::endl;
						}
					}
				}
			}
		}
	}
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
