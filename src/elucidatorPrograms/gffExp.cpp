/*
 * gffExpRunner.cpp
 *
 *  Created on: May 18, 2015
 *      Author: nick hathaway
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


#include "gffExp.hpp"
#include "elucidator/objects/BioDataObject.h"
#include "elucidator/seqToolsUtils/seqToolsUtils.hpp"

#include <njhseq/objects/Gene.h>
#include <TwoBit.h>

namespace njhseq {
gffExpRunner::gffExpRunner()
    : njh::progutils::ProgramRunner(
          {
					 addFunc("testingGffReading", testingGffReading, false),
					 addFunc("gffFeatureCount", gffFeatureCount, false),
					 addFunc("gffDescriptionsCount", gffDescriptionsCount, false),
					 addFunc("extractGffFeature", extractGffFeature, false),
					 addFunc("gffToBed", gffToBed, false),
					 addFunc("gffToBedByDescription", gffToBedByDescription, false),
					 addFunc("gffToBedByChrom", gffToBedByChrom, false),
					 addFunc("gffToBedByBedLoc", gffToBedByBedLoc, false),
					 addFunc("gffAttributesCount", gffAttributesCount, false),
					 addFunc("gffCountAttribute", gffCountAttribute, false),
					 addFunc("gffToBedByAttribute", gffToBedByAttribute, false),
					 addFunc("gffToBedByName", gffToBedByName, false),
					 addFunc("gffToJsonByID", gffToJsonByID, false),
					 addFunc("extractGffRecordWithChildren", extractGffRecordWithChildren, false),
					 addFunc("gffToBedByAttributeIncludeExonInfo", gffToBedByAttributeIncludeExonInfo, false),
					 addFunc("bedGetIntersectingGenesInGff", bedGetIntersectingGenesInGff, false),
					 addFunc("bedGetIntersectingRecordsInGff", bedGetIntersectingRecordsInGff, false),
					 addFunc("gffPrintIds", gffPrintIds, false),
					 addFunc("reorientBedToIntersectingGeneInGff", reorientBedToIntersectingGeneInGff, false),
					 addFunc("removeFastaFromGffFile", removeFastaFromGffFile, false),
					 addFunc("setBedPositionsToIntersectingGeneInGff", setBedPositionsToIntersectingGeneInGff, false),
					 addFunc("gffGetNumOfTranscriptsForGenes", gffGetNumOfTranscriptsForGenes, false),
					 addFunc("aaPositionsToBed", aaPositionsToBed, false),
					 addFunc("bedGetRegionsCompletelyInGenesInGff", bedGetRegionsCompletelyInGenesInGff, false),
           },//
          "gffExp") {}
class AmionoAcidPositionInfo {
public:
	AmionoAcidPositionInfo(const bfs::path & inputInfoFnp,
			bool zeroBased) :
			inputInfoFnp_(inputInfoFnp),
			infoTab_(inputInfoFnp, "\t", true),
			zeroBased_(zeroBased){
		infoTab_.checkForColumnsThrow(VecStr { "ID", "AAPosition" },
				__PRETTY_FUNCTION__);
		/**@todo add readding in the rest of the columns as possible meta data to be associated with amino acid position*/
		auto idsVec = infoTab_.getColumnLevels("ID");
		ids_ = std::set<std::string>(idsVec.begin(), idsVec.end());
		for (const auto & row : infoTab_.content_) {
			aminoPositionsPerId_[row[infoTab_.getColPos("ID")]].emplace(
					zeroBased ?
							njh::StrToNumConverter::stoToNum<uint32_t>(
									row[infoTab_.getColPos("AAPosition")]) :
							njh::StrToNumConverter::stoToNum<uint32_t>(
									row[infoTab_.getColPos("AAPosition")]) - 1);
		}
	}
	const bfs::path inputInfoFnp_;
	table infoTab_;
	bool zeroBased_ = false;
	std::set<std::string> ids_;
	std::map<std::string, std::set<uint32_t>> aminoPositionsPerId_;

};//
int gffExpRunner::aaPositionsToBed(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path proteinMutantTypingFnp = "";
	bfs::path gffFnp = "";
	bfs::path twoBitFnp = "";
	bool zeroBased = false;
	bool collapsePerId = false;
	OutOptions outOpts(bfs::path("out.bed"));
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.setOption(proteinMutantTypingFnp, "--aaPositionsFile", "Amino Acid positions file, must have at least two columns; ID and aaPosition", true);
	setUp.setOption(gffFnp, "--gff", "GFF (gene feature format) file", true);
	setUp.setOption(twoBitFnp, "--2bit", "2bit file of genome", true);
	setUp.setOption(zeroBased, "--zeroBased", "Zero based");
	setUp.setOption(collapsePerId, "--collapsePerId", "Create a bed record that collapses all the amino acids from id rather than each position individually");
	setUp.processWritingOptions(outOpts);

	setUp.finishSetUp(std::cout);
	OutputStream out(outOpts);
	AmionoAcidPositionInfo aaInfos(proteinMutantTypingFnp, zeroBased);
	auto genes = GeneFromGffs::getGenesFromGffForIds(gffFnp, aaInfos.ids_);

	TwoBit::TwoBitFile tReader(twoBitFnp);
	for (const auto & positions : aaInfos.aminoPositionsPerId_) {
		auto gsInfos = njh::mapAt(genes, positions.first)->generateGeneSeqInfo(tReader, false);
		/**@todo need to fix, just taking the first transcipt so this compiles for the time*/
		auto gsInfo = gsInfos.begin()->second;
		if (collapsePerId && positions.second.size() > 1) {
			std::vector<uint32_t> posVec(positions.second.begin(),
					positions.second.end());
			auto minAAPos = vectorMinimum(posVec);
			auto maxAAPos = vectorMaximum(posVec);
			auto posBed = gsInfo->genBedFromAAPositions(minAAPos, maxAAPos + 1);
			if (zeroBased) {
				posBed.name_ = njh::pasteAsStr(positions.first, "-", "[AA", minAAPos,
						"-", maxAAPos + 1, ")");
			} else {
				posBed.name_ = njh::pasteAsStr(positions.first, "-", "[AA",
						minAAPos + 1, "-", maxAAPos + 1, "]");
			}
			out << posBed.toDelimStr() << std::endl;
		} else {
			for (const auto & pos : positions.second) {
				auto posBed = gsInfo->genBedFromAAPositions(pos, pos + 1);
				if (zeroBased) {
					posBed.name_ = njh::pasteAsStr(positions.first, "-", "AA", pos);
				} else {
					posBed.name_ = njh::pasteAsStr(positions.first, "-", "AA", pos + 1);
				}
				out << posBed.toDelimStr() << std::endl;
			}
		}
	}
	return 0;
}



int gffExpRunner::gffGetNumOfTranscriptsForGenes(const njh::progutils::CmdArgs & inputCommands){
	bfs::path inputFile;
	OutOptions outOpts(bfs::path(""), ".tab.txt");
	seqSetUp setUp(inputCommands);
	setUp.description_ = "Get the number of transcripts per gene record, for annotations like the one for PlasmoDB this won't work since separate transcript are organized under their own gene record";
	setUp.setOption(inputFile, "--gff", "Input gff file", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	BioDataFileIO<GFFCore> reader{(IoOptions(InOptions(inputFile)))};
	reader.openIn();
	OutputStream out(outOpts);
	out << "GeneID\tTranscript" << std::endl;
	uint32_t count = 0;
	std::string line = "";
	std::shared_ptr<GFFCore> gRecord = reader.readNextRecord();

	std::set<std::string> parents;
	std::vector<std::shared_ptr<GFFCore>> gffRecs;
	std::shared_ptr<GeneFromGffs> genes;
	while (nullptr != gRecord) {
		if ((!gRecord->hasAttr("Parent")
				|| !njh::in(gRecord->getAttr("Parent"), parents)) && !gffRecs.empty()) {
			auto currentGene = std::make_shared<GeneFromGffs>(gffRecs);
			out << currentGene->gene_->getIDAttr()
					<< "\t" << currentGene->mRNAs_.size() << std::endl;
			gffRecs.clear();
			parents.clear();
		} else if (gRecord->hasAttr("Parent")
				&& njh::in(gRecord->getAttr("Parent"), parents)) {
			parents.insert(gRecord->getAttr("ID"));
			gffRecs.emplace_back(std::make_shared<GFFCore>(*gRecord));
		}
		if ("gene" == gRecord->type_) {
			parents.insert(gRecord->getAttr("ID"));
			gffRecs.emplace_back(std::make_shared<GFFCore>(*gRecord));
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
		++count;
	}

	return 0;
}


/**
 *

 * @param inputCommands
 * @return
 */


int gffExpRunner::removeFastaFromGffFile(const njh::progutils::CmdArgs & inputCommands){
	bfs::path inputFile;
	OutOptions outOpts(bfs::path("out.gff"));
	seqSetUp setUp(inputCommands);
	setUp.setOption(inputFile, "--gff", "Input gff file", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	BioDataFileIO<GFFCore> reader{(IoOptions(InOptions(inputFile)))};
	reader.openIn();
	uint32_t count = 0;
	std::string line = "";
	std::shared_ptr<GFFCore> gRecord = reader.readNextRecord();
	std::ofstream outFile;
	outOpts.openFile(outFile);
	{
		//write header
		std::ifstream infile(inputFile.string());
		while('#' == infile.peek()){
			njh::files::crossPlatGetline(infile, line);
			outFile << line << std::endl;
		}
	}
	std::set<std::string> parents;
	while (nullptr != gRecord) {
		gRecord->writeGffRecord(outFile);
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
		++count;
	}
	return 0;
}


int gffExpRunner::gffToBedByName(const njh::progutils::CmdArgs & inputCommands) {

	bfs::path filename = "";
	std::string name = "";
	OutOptions outOpts(bfs::path(""), ".bed");

	seqSetUp setUp(inputCommands);
	setUp.processWritingOptions(outOpts);
	setUp.setOption(filename, "--gff", "Gff File", true);
	setUp.setOption(name, "--name", "Name", true);
	setUp.finishSetUp(std::cout);


	BioDataFileIO<GFFCore> gffIo(IoOptions(InOptions(filename), outOpts));
	gffIo.openIn();
	gffIo.openOut();
	GFFCore gff;
	std::string line = "";
	while(gffIo.readNextRecord(gff)){
		if(gff.hasAttr("Name") && gff.getAttr("Name") == name){
			gffIo.write(gff, [](const GFFCore & g, std::ostream & out){
					out << GenomicRegion(g).genBedRecordCore().toDelimStr() << std::endl;
				;});
		}
		bool end = false;
		while('#' == gffIo.inFile_->peek()){
			if (njh::files::nextLineBeginsWith(*gffIo.inFile_, "##FASTA")) {
				end = true;
				break;
			}
			njh::files::crossPlatGetline(*gffIo.inFile_, line);
		}
		if(end){
			break;
		}
	}
	return 0;
}

int gffExpRunner::gffToJsonByID(const njh::progutils::CmdArgs & inputCommands) {

	bfs::path filename = "";
	std::string ID = "";
	OutOptions outOpts(bfs::path(""), ".json");

	seqSetUp setUp(inputCommands);
	setUp.processWritingOptions(outOpts);
	setUp.setOption(filename, "--gff", "Gff File", true);
	setUp.setOption(ID, "--ID", "Gff recrod ID", true);
	setUp.finishSetUp(std::cout);


	BioDataFileIO<GFFCore> gffIo{IoOptions(InOptions(filename))};
	gffIo.openIn();
	OutputStream out(outOpts);
	GFFCore gff;
	std::string line = "";
	while(gffIo.readNextRecord(gff)){
		if(gff.hasAttr("ID") && gff.getAttr("ID") == ID){
			out << gff.toJson() << std::endl;
		}
		bool end = false;
		while('#' == gffIo.inFile_->peek()){
			if (njh::files::nextLineBeginsWith(*gffIo.inFile_, "##FASTA")) {
				end = true;
				break;
			}
			njh::files::crossPlatGetline(*gffIo.inFile_, line);
		}
		if(end){
			break;
		}
	}
	return 0;
}

int gffExpRunner::testingGffReading(const njh::progutils::CmdArgs & inputCommands){
	bfs::path inputFile;
	seqSetUp setUp(inputCommands);
	setUp.setOption(inputFile, "--gff", "Input gff file", true);
	setUp.finishSetUp(std::cout);


	std::vector<std::shared_ptr<GFFCore>> ret;
	BioDataFileIO<GFFCore> reader((IoOptions(InOptions(inputFile))));
	reader.openIn();
	uint32_t count = 0;
	std::string line = "";
	std::shared_ptr<GFFCore> gRecord = reader.readNextRecord();
	while (nullptr != gRecord) {
		ret.emplace_back(std::move(gRecord));
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
		++count;

	}
	for(const auto & g : ret){
		if(g->getAttr("ID") == "PF3D7_0100100" || g->getAttr("Parent") == "PF3D7_0100100"){
			std::cout << g->toJson() << std::endl;
		}
		if(std::numeric_limits<uint16_t>::max() != g->phase_){
			std::cout << g->toJson() << std::endl;
			break;
		}
	}

	return 0;
}

int gffExpRunner::gffFeatureCount(const njh::progutils::CmdArgs & inputCommands){
	bfs::path inputFile;
	seqSetUp setUp(inputCommands);
	setUp.setOption(inputFile, "--gff", "Input gff file", true);
	setUp.finishSetUp(std::cout);

	BioDataFileIO<GFFCore> reader((IoOptions(InOptions(inputFile))));
	reader.openIn();
	uint32_t count = 0;
	std::string line = "";
	std::shared_ptr<GFFCore> gRecord = reader.readNextRecord();
	std::unordered_map<std::string, uint32_t> featureCounts;
	while (nullptr != gRecord) {
		++featureCounts[gRecord->type_];
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
		++count;

	}

	table out(featureCounts, VecStr{"feature", "count"});

	out.sortTable("feature", false);
	out.outPutContents(std::cout, "\t");
	return 0;
}

int gffExpRunner::gffDescriptionsCount(const njh::progutils::CmdArgs & inputCommands){
	bfs::path inputFile;
	seqSetUp setUp(inputCommands);
	setUp.setOption(inputFile, "--gff", "Input gff file", true);
	setUp.finishSetUp(std::cout);

	BioDataFileIO<GFFCore> reader((IoOptions(InOptions(inputFile))));
	reader.openIn();
	uint32_t count = 0;
	std::string line = "";
	std::shared_ptr<GFFCore> gRecord = reader.readNextRecord();
	std::unordered_map<std::string, uint32_t> descriptionCounts;
	while (nullptr != gRecord) {
		if(gRecord->hasAttr("description")){
			++descriptionCounts[gRecord->getAttr("description")];
		}else{
			++descriptionCounts["none"];
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
		++count;

	}

	table out(descriptionCounts, VecStr{"description", "count"});

	out.sortTable("count", false);
	out.outPutContents(std::cout, "\t");
	return 0;
}

int gffExpRunner::gffCountAttribute(const njh::progutils::CmdArgs & inputCommands){
	bfs::path inputFile = "";
	std::string attribute = "";
	seqSetUp setUp(inputCommands);
	setUp.setOption(inputFile, "--gff", "Input gff file", true);
	setUp.setOption(attribute, "--attribute", "Attribute to count", true);
	setUp.finishSetUp(std::cout);

	BioDataFileIO<GFFCore> reader((IoOptions(InOptions(inputFile))));
	reader.openIn();
	uint32_t count = 0;
	std::string line = "";
	std::shared_ptr<GFFCore> gRecord = reader.readNextRecord();
	std::unordered_map<std::string, uint32_t> descriptionCounts;
	while (nullptr != gRecord) {
		if(gRecord->hasAttr(attribute)){
			++descriptionCounts[gRecord->getAttr(attribute)];
		}else{
			++descriptionCounts["none"];
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
		++count;

	}

	table out(descriptionCounts, VecStr{attribute, "count"});

	out.sortTable("count", false);
	out.outPutContents(std::cout, "\t");
	return 0;
}

int gffExpRunner::gffAttributesCount(const njh::progutils::CmdArgs & inputCommands){
	bfs::path inputFile;
	seqSetUp setUp(inputCommands);
	setUp.setOption(inputFile, "--gff", "Input gff file", true);
	setUp.finishSetUp(std::cout);

	BioDataFileIO<GFFCore> reader((IoOptions(InOptions(inputFile))));
	reader.openIn();
	uint32_t count = 0;
	std::string line = "";
	std::shared_ptr<GFFCore> gRecord = reader.readNextRecord();
	std::unordered_map<std::string, uint32_t> descriptionCounts;
	while (nullptr != gRecord) {
		for(const auto & att : gRecord->attributes_){
			++descriptionCounts[att.first];
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
		++count;

	}

	table out(descriptionCounts, VecStr{"attribute", "count"});

	out.sortTable("count", false);
	out.outPutContents(std::cout, "\t");
	return 0;
}







int gffExpRunner::gffToBedByAttributeIncludeExonInfo(const njh::progutils::CmdArgs & inputCommands){
	bfs::path inputFile = "";
	OutOptions outOpts(bfs::path(""));
	outOpts.outExtention_ = ".bed";
	std::string attribute = "";
	std::string attributeLevel = "";

	seqSetUp setUp(inputCommands);
	setUp.setOption(inputFile, "--gff", "Input gff file", true);
	setUp.setOption(attribute, "--attribute", "attribute to check of gff to extract", true);
	setUp.setOption(attributeLevel, "--attributeLevel", "attribute level to match in attribute of gff to extract", true);

	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	BioDataFileIO<GFFCore> reader{IoOptions(InOptions(inputFile))};
	reader.openIn();
	uint32_t count = 0;
	std::string line = "";
	std::shared_ptr<GFFCore> gRecord = reader.readNextRecord();
	OutOptions outOptsJson(outOpts.outFilename_);
	outOptsJson.outExtention_ = ".json";
	outOptsJson.transferOverwriteOpts(outOpts);
	std::ofstream outJsonFile;
	outOptsJson.openFile(outJsonFile);

	std::ofstream outputFile;
	std::ostream out(outOpts.determineOutBuf(outputFile));

	std::unordered_map<std::string, std::set<std::string>> parents;
	std::unordered_map<std::string, std::vector<std::shared_ptr<GFFCore>>> gffRecs;
	while (nullptr != gRecord) {
		if(gRecord->hasAttr(attribute) && gRecord->getAttr(attribute) == attributeLevel){
			parents[gRecord->getAttr("ID")].insert(gRecord->getAttr("ID"));
			gffRecs[gRecord->getAttr("ID")].emplace_back(std::make_shared<GFFCore>(*gRecord));
		}else if(gRecord->hasAttr("Parent")){
			for(const auto & pars : parents){
				if(njh::in(gRecord->getAttr("Parent"), pars.second)){
					parents[pars.first].insert(gRecord->getAttr("ID"));
					gffRecs[pars.first].emplace_back(std::make_shared<GFFCore>(*gRecord));
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
		++count;
	}
	std::unordered_map<std::string, std::shared_ptr<GeneFromGffs>> gffGenes;
	for(const auto & gffRec : gffRecs){
		gffGenes.emplace(gffRec.first, std::make_shared<GeneFromGffs>(gffRec.second));
	}
	Json::Value outJson;
	for(const auto & gffGene : gffGenes){
		GenomicRegion geneRegion (*gffGene.second->gene_);
		outJson[gffGene.second->gene_->getAttr("Name")] = gffGene.second->gene_->toJson();
		out << GenomicRegion(*gffGene.second->gene_).genBedRecordCore().toDelimStr() << std::endl;
		VecStr exonsAlreadyAdded;
		uint32_t exonCount = 0;
		for(const auto & exonByTranscript : gffGene.second->exons_){
			for(const auto & exon : exonByTranscript.second){
				GenomicRegion exonRegion(*exon);

				if(!njh::in(exonRegion.createUidFromCoords(), exonsAlreadyAdded)){
					++exonCount;
					exonRegion.uid_ = geneRegion.uid_ + "_exon_" + estd::to_string(exonCount);
					exonsAlreadyAdded.emplace_back(exonRegion.createUidFromCoords());
					outJson[exon->getAttr("Name")] = exon->toJson();
					out << exonRegion.genBedRecordCore().toDelimStr() << std::endl;
				}
			}
		}
	}


	outJsonFile << outJson << std::endl;
	return 0;
}

int gffExpRunner::gffPrintIds(const njh::progutils::CmdArgs & inputCommands){
	bfs::path inputFile = "";
	OutOptions outOpts(bfs::path(""));
	outOpts.outExtention_ = ".txt";
	seqSetUp setUp(inputCommands);
	setUp.setOption(inputFile, "--gff", "Input gff file", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	BioDataFileIO<GFFCore> reader{IoOptions(InOptions(inputFile))};
	reader.openIn();
	uint32_t count = 0;
	std::string line = "";
	std::shared_ptr<GFFCore> gRecord = reader.readNextRecord();
	OutputStream out(outOpts);

	while (nullptr != gRecord) {
		if(gRecord->hasAttr("ID")){
			out << gRecord->getAttr("ID") << std::endl;
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
		++count;
	}

	return 0;
}





int gffExpRunner::extractGffRecordWithChildren(const njh::progutils::CmdArgs & inputCommands){
	bfs::path inputFile;
	OutOptions outOpts(bfs::path("out.gff"));
	std::string id = "";
	seqSetUp setUp(inputCommands);
	setUp.setOption(inputFile, "--gff", "Input gff file", true);
	setUp.setOption(id, "--id", "Id to extract", true);

	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	BioDataFileIO<GFFCore> reader((IoOptions(InOptions(inputFile))));
	reader.openIn();
	uint32_t count = 0;
	std::string line = "";
	std::shared_ptr<GFFCore> gRecord = reader.readNextRecord();
	std::ofstream outFile;
	outOpts.openFile(outFile);
	{
		//write header
		std::ifstream infile(inputFile.string());
		while('#' == infile.peek()){
			njh::files::crossPlatGetline(infile, line);
			outFile << line << std::endl;
		}
	}
	std::set<std::string> parents;
	while (nullptr != gRecord) {
		if(gRecord->hasAttr("ID") && gRecord->getAttr("ID") == id){
			parents.insert(gRecord->getAttr("ID"));
			gRecord->writeGffRecord(outFile);
		}else if(gRecord->hasAttr("Parent") && njh::in(gRecord->getAttr("Parent"), parents)){
			parents.insert(gRecord->getAttr("ID"));
			gRecord->writeGffRecord(outFile);
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
		++count;
	}
	return 0;
}

int gffExpRunner::extractGffFeature(const njh::progutils::CmdArgs & inputCommands){
	bfs::path inputFile;
	OutOptions outOpts(bfs::path("out.gff"));
	std::string feature = "";
	seqSetUp setUp(inputCommands);
	setUp.setOption(inputFile, "--gff", "Input gff file", true);
	setUp.setOption(feature, "--feature", "Feature to extract", true);

	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	BioDataFileIO<GFFCore> reader((IoOptions(InOptions(inputFile))));
	reader.openIn();
	uint32_t count = 0;
	std::string line = "";
	std::shared_ptr<GFFCore> gRecord = reader.readNextRecord();
	std::unordered_map<std::string, uint32_t> featureCounts;
	std::ofstream outFile;
	outOpts.openFile(outFile);
	{
		//write header
		std::ifstream infile(inputFile.string());
		while('#' == infile.peek()){
			njh::files::crossPlatGetline(infile, line);
			outFile << line << std::endl;
		}
	}
	while (nullptr != gRecord) {
		++featureCounts[gRecord->type_];
		if(feature == gRecord->type_){
			gRecord->writeGffRecord(outFile);
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
		++count;
	}

	if (setUp.pars_.verbose_) {
		table out(featureCounts, VecStr { "feature", "count" });

		out.sortTable("feature", false);
		out.outPutContents(std::cout, "\t");
	}
	return 0;
}

int gffExpRunner::gffToBed(const njh::progutils::CmdArgs & inputCommands){
	bfs::path inputFile;
	OutOptions outOpts(bfs::path("out"));
	outOpts.outExtention_ = ".bed";

	seqSetUp setUp(inputCommands);
	setUp.setOption(inputFile, "--gff", "Input gff file", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	BioDataFileIO<GFFCore> reader((IoOptions(InOptions(inputFile))));
	reader.openIn();
	uint32_t count = 0;
	std::string line = "";
	std::shared_ptr<GFFCore> gRecord = reader.readNextRecord();
	std::ofstream outFile;
	outOpts.openFile(outFile);

	OutOptions outOptsJson(outOpts.outFilename_);
	outOptsJson.outExtention_ = ".json";
	outOptsJson.transferOverwriteOpts(outOpts);
	std::ofstream outJsonFile;
	outOptsJson.openFile(outJsonFile);

	Json::Value outJson;
	while (nullptr != gRecord) {

		outJson[gRecord->getAttr("Name")] = gRecord->toJson();
		outFile << GenomicRegion(*gRecord).genBedRecordCore().toDelimStr() << std::endl;
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
		++count;
	}
	outJsonFile << outJson << std::endl;
	return 0;
}

int gffExpRunner::gffToBedByDescription(const njh::progutils::CmdArgs & inputCommands){
	bfs::path inputFile;
	OutOptions outOpts(bfs::path("out"));
	outOpts.outExtention_ = ".bed";
	VecStr description;
	seqSetUp setUp(inputCommands);
	setUp.setOption(inputFile, "--gff", "Input gff file", true);
	setUp.setOption(description, "--description", "description of gff to extract", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	BioDataFileIO<GFFCore> reader((IoOptions(InOptions(inputFile))));
	reader.openIn();
	uint32_t count = 0;
	std::string line = "";
	std::shared_ptr<GFFCore> gRecord = reader.readNextRecord();
	std::ofstream outFile;
	outOpts.openFile(outFile);

	OutOptions outOptsJson(outOpts.outFilename_);
	outOptsJson.outExtention_ = ".json";
	outOptsJson.transferOverwriteOpts(outOpts);
	std::ofstream outJsonFile;
	outOptsJson.openFile(outJsonFile);

	Json::Value outJson;
	while (nullptr != gRecord) {
		if(njh::in(gRecord->getAttr("description"), description)){
			outJson[gRecord->getAttr("Name")] = gRecord->toJson();
			outFile << GenomicRegion(*gRecord).genBedRecordCore().toDelimStr() << std::endl;
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
		++count;
	}
	outJsonFile << outJson << std::endl;
	return 0;
}


int gffExpRunner::gffToBedByAttribute(const njh::progutils::CmdArgs & inputCommands){
	bfs::path inputFile;
	OutOptions outOpts(bfs::path("out"));
	outOpts.outExtention_ = ".bed";
	std::string attribute = "";
	std::string attributeLevel = "";

	seqSetUp setUp(inputCommands);
	setUp.setOption(inputFile, "--gff", "Input gff file", true);
	setUp.setOption(attribute, "--attribute", "attribute to check of gff to extract", true);
	setUp.setOption(attributeLevel, "--attributeLevel", "attribute level to match in attribute of gff to extract", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	BioDataFileIO<GFFCore> reader((IoOptions(InOptions(inputFile))));
	reader.openIn();
	uint32_t count = 0;
	std::string line = "";
	std::shared_ptr<GFFCore> gRecord = reader.readNextRecord();
	std::ofstream outFile;
	outOpts.openFile(outFile);

	OutOptions outOptsJson(outOpts.outFilename_);
	outOptsJson.outExtention_ = ".json";
	outOptsJson.transferOverwriteOpts(outOpts);
	std::ofstream outJsonFile;
	outOptsJson.openFile(outJsonFile);

	Json::Value outJson;
	while (nullptr != gRecord) {
		if(gRecord->hasAttr(attribute) && gRecord->getAttr(attribute) == attributeLevel){
			outJson[gRecord->getAttr("Name")] = gRecord->toJson();
			outFile << GenomicRegion(*gRecord).genBedRecordCore().toDelimStr() << std::endl;
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
		++count;
	}
	outJsonFile << outJson << std::endl;
	return 0;
}



int gffExpRunner::gffToBedByChrom(const njh::progutils::CmdArgs & inputCommands){
	bfs::path inputFile;
	OutOptions outOpts(bfs::path("out"));
	outOpts.outExtention_ = ".bed";
	std::string chrom = "";
	seqSetUp setUp(inputCommands);
	setUp.setOption(inputFile, "--gff", "Input gff file", true);
	setUp.setOption(chrom, "--chrom", "chrom in gff to extract", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	BioDataFileIO<GFFCore> reader((IoOptions(InOptions(inputFile))));
	reader.openIn();
	uint32_t count = 0;
	std::string line = "";
	std::shared_ptr<GFFCore> gRecord = reader.readNextRecord();
	std::ofstream outFile;
	outOpts.openFile(outFile);

	OutOptions outOptsJson(outOpts.outFilename_);
	outOptsJson.outExtention_ = ".json";
	outOptsJson.transferOverwriteOpts(outOpts);
	std::ofstream outJsonFile;
	outOptsJson.openFile(outJsonFile);

	Json::Value outJson;
	while (nullptr != gRecord) {
		if(gRecord->seqid_ == chrom){
			outJson[gRecord->getAttr("Name")] = gRecord->toJson();
			outFile << GenomicRegion(*gRecord).genBedRecordCore().toDelimStr() << std::endl;
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
		++count;
	}
	outJsonFile << outJson << std::endl;
	return 0;
}

int gffExpRunner::gffToBedByBedLoc(const njh::progutils::CmdArgs & inputCommands){
	bfs::path inputFile;
	OutOptions outOpts(bfs::path("out"));
	outOpts.outExtention_ = ".bed";
	bfs::path bedFnp = "";
	size_t overlapMin = 1;
	std::string feature = "";
	std::string extraAttributesStr = "";

	seqSetUp setUp(inputCommands);
	setUp.setOption(extraAttributesStr, "--extraAttributes", "Extra Attributes to output");
	setUp.setOption(overlapMin, "--overlapMin", "overlap minimum");
	setUp.setOption(feature, "--feature", "Only extract records of that match this feature");
	setUp.setOption(inputFile, "--gff", "Input gff file", true);
	setUp.setOption(bedFnp, "--bed", "Bed regions to extract", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	auto inputRegions = gatherRegions(bedFnp.string(), "", setUp.pars_.verbose_);
	BioDataFileIO<GFFCore> reader((IoOptions(InOptions(inputFile))));
	reader.openIn();
	uint32_t count = 0;
	std::string line = "";
	std::shared_ptr<GFFCore> gRecord = reader.readNextRecord();
	std::ofstream outFile;
	outOpts.openFile(outFile);

	OutOptions outOptsJson(outOpts.outFilename_);
	outOptsJson.outExtention_ = ".json";
	outOptsJson.transferOverwriteOpts(outOpts);
	std::ofstream outJsonFile;
	outOptsJson.openFile(outJsonFile);
	GenomicRegion gTest;
	auto extraAttributes = tokenizeString(extraAttributesStr, ",");
	if(!extraAttributes.empty() && !njh::in(std::string("ID"), extraAttributes)){
		extraAttributes.emplace_back("ID");
	}

	Json::Value outJson;
	while (nullptr != gRecord) {

		auto gRegion = GenomicRegion(*gRecord);
		for(const auto & inputRegion : inputRegions){

			if(inputRegion.overlaps(gRegion)){
				if("" != feature && gRecord->type_ != feature){
					continue;
				}
				outJson[gRecord->getAttr("Name")] = gRecord->toJson();
				auto outRegion = GenomicRegion(*gRecord).genBedRecordCore();

				if("" != extraAttributesStr){
					std::string extra = "";
					extra.append("[");
					for(const auto & attr : extraAttributes){
						if(gRecord->hasAttr(attr)){
							extra.append(attr + "=" + gRecord->getAttr(attr) + ";");
						}else{
							extra.append(attr + "=" + "NA" + ";");
						}
					}
					extra.append("]");
					outRegion.extraFields_.emplace_back(extra);
				}

				outFile << outRegion.toDelimStrWithExtra() << std::endl;
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
		++count;
	}
	outJsonFile << outJson << std::endl;
	return 0;
}



int gffExpRunner::setBedPositionsToIntersectingGeneInGff(const njh::progutils::CmdArgs & inputCommands){
	bfs::path inputFile;
	OutOptions outOpts(bfs::path(""));
	outOpts.outExtention_ = ".bed";
	bfs::path bedFnp = "";
	size_t overlapMin = 1;
	bool reOrient = false;
	seqSetUp setUp(inputCommands);
	setUp.description_ = "Will set the positions of a bed record to fall only within the region of the first intersected gene";
	setUp.setOption(overlapMin, "--overlapMin", "overlap minimum");
	setUp.setOption(inputFile, "--gff", "Input gff file", true);
	setUp.setOption(bedFnp, "--bed", "Bed regions to extract", true);
	setUp.setOption(reOrient, "--reOrient", "Will also set strand to the gene as well");
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	auto beds = getBeds(bedFnp);

	BioDataFileIO<GFFCore> reader{IoOptions(InOptions(inputFile))};
	reader.openIn();
	uint32_t count = 0;
	std::string line = "";
	std::shared_ptr<GFFCore> gRecord = reader.readNextRecord();
	OutputStream out(outOpts);

//	OutOptions outOptsJson(outOpts.outFilename_);
//	outOptsJson.outExtention_ = ".json";
//	outOptsJson.transferOverwriteOpts(outOpts);
//	std::ofstream outJsonFile;
//	outOptsJson.openFile(outJsonFile);
//	Json::Value outputInfo;

//	for(const auto & inputRegion : beds){
//		inputRegion->extraFields_.emplace_back("");
//	}
	Json::Value outJson;

	while (nullptr != gRecord) {
		if("gene" == gRecord->type_){
			auto gRegion = GenomicRegion(*gRecord);
			for(auto & inputRegion : beds){
				if(GenomicRegion(*inputRegion).overlaps(gRegion)){
					if(reOrient){
						inputRegion->strand_ = gRegion.reverseSrand_ ? '-' : '+';
					}
					bool resetScore = inputRegion->score_ == inputRegion->length();
					if(inputRegion->chromStart_ < gRegion.start_){
						inputRegion->chromStart_ = gRegion.start_;
					}
					if (inputRegion->chromEnd_ > gRegion.end_) {
						inputRegion->chromEnd_ = gRegion.end_;
					}
					if (resetScore) {
						inputRegion->score_ = inputRegion->length();
					}
					break;
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
		++count;
	}
	for(const auto & bed : beds){
		out << bed->toDelimStrWithExtra() << std::endl;
	}
//	outJsonFile << outJson << std::endl;
	return 0;
}


int gffExpRunner::reorientBedToIntersectingGeneInGff(const njh::progutils::CmdArgs & inputCommands){
	bfs::path inputFile;
	OutOptions outOpts(bfs::path(""));
	outOpts.outExtention_ = ".bed";
	bfs::path bedFnp = "";
	size_t overlapMin = 1;
	VecStr features{"gene"};
	seqSetUp setUp(inputCommands);
	setUp.description_ = "Will set strand to first gene intersected with in the gff strand";
	setUp.setOption(features, "--features", "features to use");

	setUp.setOption(overlapMin, "--overlapMin", "overlap minimum");
	setUp.setOption(inputFile, "--gff", "Input gff file", true);
	setUp.setOption(bedFnp, "--bed", "Bed regions to extract", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	auto beds = getBeds(bedFnp);

	BioDataFileIO<GFFCore> reader{IoOptions(InOptions(inputFile))};
	reader.openIn();
	uint32_t count = 0;
	std::string line = "";
	std::shared_ptr<GFFCore> gRecord = reader.readNextRecord();
	OutputStream out(outOpts);

	OutOptions outOptsJson(outOpts.outFilename_);
	outOptsJson.outExtention_ = ".json";
	outOptsJson.transferOverwriteOpts(outOpts);

	while (nullptr != gRecord) {

		if(njh::in(gRecord->type_, features)){
			auto gRegion = GenomicRegion(*gRecord);
			for(auto & inputRegion : beds){
				if(gRegion.overlaps(*inputRegion)){
					inputRegion->strand_ = gRegion.reverseSrand_ ? '-' : '+';
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
		++count;
	}
	for(const auto & bed : beds){
		out << bed->toDelimStrWithExtra() << std::endl;
	}
//	outJsonFile << outJson << std::endl;
	return 0;
}

template <typename BEDREC>
std::vector<BEDREC> getBedsCompletelyInGenesInGFf(
		const std::vector<BEDREC> & beds,
		const intersectBedLocsWtihGffRecordsPars & pars){
	std::vector<BEDREC>  ret;
	std::unordered_map<std::string, std::vector<uint32_t>> bedsByChrome;

	BioDataFileIO<GFFCore> reader { IoOptions(InOptions(pars.gffFnp_)) };
	reader.openIn();
	uint32_t count = 0;
	std::string line = "";
	std::shared_ptr<GFFCore> gRecord = reader.readNextRecord();
	for (const auto & bPos : iter::range(beds.size())) {
		bedsByChrome[getRef(beds[bPos]).chrom_].emplace_back(bPos);
	}
	while (nullptr != gRecord) {
		if (pars.selectFeatures_.empty() || njh::in(gRecord->type_, pars.selectFeatures_)) {
			auto gRegion = GenomicRegion(*gRecord);
			for (auto & inputRegionPos : bedsByChrome[gRegion.chrom_]) {
				if (gRegion.getOverlapLen(getRef(beds[inputRegionPos])) == getRef(beds[inputRegionPos]).length()) {
					ret.push_back(beds[inputRegionPos]);
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
		++count;
	}
	return ret;
}




int gffExpRunner::bedGetRegionsCompletelyInGenesInGff(const njh::progutils::CmdArgs & inputCommands){
	intersectBedLocsWtihGffRecordsPars pars;
	OutOptions outOpts("out", ".bed");
	bfs::path bedFnp = "";
	pars.selectFeatures_ = VecStr{"gene"};

	seqSetUp setUp(inputCommands);
	setUp.setOption(pars.selectFeatures_, "--selectFeatures", "Gff Features to consider");

	setUp.setOption(pars.gffFnp_, "--gff", "Input gff file", true);
	setUp.setOption(bedFnp, "--bed", "Bed regions to extract", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);
	auto beds = getBed3s(bedFnp);
	njh::files::checkExistenceThrow(pars.gffFnp_, __PRETTY_FUNCTION__);
	OutputStream out(outOpts);
	auto intersectedBeds = getBedsCompletelyInGenesInGFf(beds, pars);
	for(const auto & bed : intersectedBeds){
		out << bed->toDelimStrWithExtra() << std::endl;
	}
	return 0;
}

int gffExpRunner::bedGetIntersectingGenesInGff(const njh::progutils::CmdArgs & inputCommands){
	intersectBedLocsWtihGffRecordsPars pars;
	OutOptions outOpts("out", ".bed");
	bfs::path bedFnp = "";
	pars.selectFeatures_ = VecStr{"gene"};
	seqSetUp setUp(inputCommands);
	setUp.setOption(pars.extraAttributes_, "--extraAttributes", "Extra Attributes to output");
	setUp.setOption(pars.selectFeatures_, "--selectFeatures", "Gff Features to consider");
	setUp.setOption(pars.gffFnp_, "--gff", "Input gff file", true);
	setUp.setOption(bedFnp, "--bed", "Bed regions to extract", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);
	auto beds = getBed3s(bedFnp);
	njh::files::checkExistenceThrow(pars.gffFnp_, __PRETTY_FUNCTION__);

	OutputStream out(outOpts);
	OutOptions outOptsJson(outOpts.outFilename_);
	outOptsJson.outExtention_ = ".json";
	outOptsJson.transferOverwriteOpts(outOpts);
	OutputStream outJsonOut(outOptsJson);
	Json::Value outJson = intersectBedLocsWtihGffRecords(beds, pars);
	for(const auto & bed : beds){
		out << bed->toDelimStrWithExtra() << std::endl;
	}
	outJsonOut << outJson << std::endl;
	return 0;
}

int gffExpRunner::bedGetIntersectingRecordsInGff(const njh::progutils::CmdArgs & inputCommands){
	intersectBedLocsWtihGffRecordsPars pars;
	OutOptions outOpts(bfs::path("out"));
	outOpts.outExtention_ = ".bed";
	bfs::path bedFnp = "";
	VecStr extraAttributes;
	size_t overlapMin = 1;
	VecStr features;
	seqSetUp setUp(inputCommands);
	setUp.setOption(extraAttributes, "--extraAttributes", "Extra Attributes to output, comma separated");
	setUp.setOption(features, "--features", "Feature to extract, if left blank all features are extracted, or can be comma separated for multiple features");
	setUp.setOption(overlapMin, "--overlapMin", "overlap minimum");
	setUp.setOption(pars.gffFnp_, "--gff", "Input gff file", true);
	setUp.setOption(bedFnp, "--bed", "Bed regions to extract", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);


	pars.extraAttributes_ = extraAttributes;
	pars.selectFeatures_ = features;
	njh::files::checkExistenceThrow(pars.gffFnp_, __PRETTY_FUNCTION__);

	auto beds = getBed3s(bedFnp);
	OutputStream out(outOpts);
	OutOptions outOptsJson(outOpts.outFilename_);
	outOptsJson.outExtention_ = ".json";
	outOptsJson.transferOverwriteOpts(outOpts);
	OutputStream outJsonOut(outOptsJson);
	Json::Value outJson = intersectBedLocsWtihGffRecords(beds, pars);
	for(const auto & bed : beds){
		out << bed->toDelimStrWithExtra() << std::endl;
	}
	outJsonOut << outJson << std::endl;


	return 0;
}


} /* namespace njhseq */
