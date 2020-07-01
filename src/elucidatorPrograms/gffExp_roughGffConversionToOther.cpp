/*
 * gffExp_roughGffConversionToOther.cpp
 *
 *  Created on: Jun 29, 2020
 *      Author: nick
 */




#include "gffExp.hpp"
#include "elucidator/objects/BioDataObject.h"
#include "elucidator/seqToolsUtils/seqToolsUtils.hpp"

#include <njhseq/objects/Gene.h>
#include <TwoBit.h>




namespace njhseq {

int gffExpRunner::roughGffConversionToOther(const njh::progutils::CmdArgs & inputCommands){
	bfs::path inputFile;
	std::string identifier;
	VecStr toreplace;
	VecStr removeAttributes{"Dbxref", "translation", "polypeptide_domain", "orthologous_to"};
	VecStr removeTypes{"region"};
	OutOptions outOpts(bfs::path("out.gff"));
	seqInfo refChrom("ref");
	seqInfo seqChrom("seq");
	seqSetUp setUp(inputCommands);

	setUp.setOption(inputFile, "--gff", "Input gff file", true);
	setUp.processSeq(refChrom, "--refChrom", "Reference chromosome to convert gff to input seq", true);
	setUp.processSeq(seqChrom, "--seqChrom", "Seq chromosome to convert gff to input seq", true);
	setUp.setOption(identifier, "--identifier", "Identifier to use in the ID etc gff fields", true);
	setUp.setOption(toreplace, "--toreplace", "A series of text to replace in ID gff fields", true);
	setUp.setOption(removeAttributes, "--removeAttributes", "Remove Attributes");
	setUp.setOption(removeTypes, "--removeTypes", "Remove feature types gff records");

	setUp.processWritingOptions(outOpts);

	setUp.pars_.gapLeft_ = "0,0";
	setUp.pars_.gapRight_ = "0,0";
	setUp.pars_.gapInfo_ = gapScoringParameters(5,1,0,0,0,0);
	setUp.pars_.colOpts_.iTOpts_.weighHomopolyer_ = false;
	setUp.processAlignerDefualts();
  if(setUp.pars_.verbose_){
    std::cout << "Gap open: " << setUp.pars_.gapInfo_.gapOpen_ << std::endl;
    std::cout << "Gap extend: " << setUp.pars_.gapInfo_.gapExtend_ << std::endl;
    std::cout << "Gap Left Query open: " << setUp.pars_.gapInfo_.gapLeftQueryOpen_ << std::endl;
    std::cout << "Gap Left Query extend: " << setUp.pars_.gapInfo_.gapLeftQueryExtend_ << std::endl;
    std::cout << "Gap Right Query open: " << setUp.pars_.gapInfo_.gapRightQueryOpen_ << std::endl;
    std::cout << "Gap Right Query extend: " << setUp.pars_.gapInfo_.gapRightQueryExtend_ << std::endl;
    std::cout << "Gap Left Ref open: " << setUp.pars_.gapInfo_.gapLeftRefOpen_ << std::endl;
    std::cout << "Gap Left Ref extend: " << setUp.pars_.gapInfo_.gapLeftRefExtend_ << std::endl;
    std::cout << "Gap Right Ref open: " << setUp.pars_.gapInfo_.gapRightRefOpen_ << std::endl;
    std::cout << "Gap Right Ref extend: " << setUp.pars_.gapInfo_.gapRightRefExtend_ << std::endl;
  }
	setUp.finishSetUp(std::cout);

	trimAtFirstWhitespace(refChrom.name_);
	trimAtFirstWhitespace(seqChrom.name_);


	BioDataFileIO<GFFCore> reader{(IoOptions(InOptions(inputFile)))};
	reader.openIn();
	uint32_t count = 0;
	std::string line = "";

	OutputStream outFile(outOpts);
	outFile << "##gff-version 3" << std::endl;
	outFile << "##sequence-region   " << seqChrom.name_<< " 1 " << seqChrom.seq_.size() << std::endl;
	std::vector<std::shared_ptr<GFFCore>> allRecords;
	std::shared_ptr<GFFCore> gRecord = reader.readNextRecord();
	std::stringstream endOfFile;
	while(nullptr != gRecord) {
		if(refChrom.name_ == gRecord->seqid_ && !njh::in(gRecord->type_,removeTypes)){
			allRecords.emplace_back(gRecord);
		}
		bool end = false;
		while ('#' == reader.inFile_->peek()) {
			if (njh::files::nextLineBeginsWith(*reader.inFile_, "##FASTA")) {
				//write out the fasta if there
				while('#' == reader.inFile_->peek()){
					njh::files::crossPlatGetline(*reader.inFile_, line);
					endOfFile << line << "\n";
				}
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


	uint64_t maxLen = 0;
	readVec::getMaxLength(refChrom, maxLen);
	readVec::getMaxLength(seqChrom, maxLen);

	KmerMaps emptyMaps;
	aligner alignerObj(maxLen, setUp.pars_.gapInfo_,
			setUp.pars_.scoring_, emptyMaps, setUp.pars_.qScorePars_,
			setUp.pars_.colOpts_.alignOpts_.countEndGaps_,
			setUp.pars_.colOpts_.iTOpts_.weighHomopolyer_);


	alignerObj.alignCacheGlobal(refChrom, seqChrom);
	alignerObj.rearrangeObjsGlobal(refChrom, seqChrom);

	std::unordered_map<uint32_t, uint32_t> refToSeqPositionKey;
	std::unordered_map<uint32_t, uint32_t> seqToRefPositionKey;

	uint32_t refPos = 0;
	uint32_t seqPos = 0;

	for(const auto pos : iter::range(alignerObj.alignObjectA_.seqBase_.seq_.size())){
		if('-' != alignerObj.alignObjectA_.seqBase_.seq_[pos]){
			refToSeqPositionKey[refPos] = seqPos;
			++refPos;
		}
		if('-' != alignerObj.alignObjectB_.seqBase_.seq_[pos]){
			seqToRefPositionKey[seqPos] = refPos;
			++seqPos;
		}
	}

	for(const auto & rec : allRecords){
		rec->seqid_ = seqChrom.name_;
		for(const auto & attr : removeAttributes){
			if(rec->hasAttr(attr)){
				rec->attributes_.erase(attr);
			}
		}
		if("gene" == rec->type_ || "pseudogene" == rec->type_){
			if(rec->hasAttr("orthologous_to")){
				rec->attributes_["orthologous_to"].append(";");
			}
			rec->attributes_["orthologous_to"].append(rec->attributes_["ID"]);
		}
		if(rec->hasAttr("ID")){
			for(const auto & rep : toreplace){
				rec->attributes_["ID"] = njh::replaceString(rec->attributes_["ID"], rep, identifier);
			}
		}
		if(rec->hasAttr("Parent")){
			for(const auto & rep : toreplace){
				rec->attributes_["Parent"] = njh::replaceString(rec->attributes_["Parent"], rep, identifier);
			}
		}
		rec->source_ = njh::pasteAsStr("pseudo", rec->source_);
		rec->start_ = refToSeqPositionKey[rec->start_ - 1] + 1;
		rec->end_ = refToSeqPositionKey[rec->end_ - 1] + 1;
		rec->writeGffRecord(outFile);
	}


	return 0;
}
}  // namespace njhseq
