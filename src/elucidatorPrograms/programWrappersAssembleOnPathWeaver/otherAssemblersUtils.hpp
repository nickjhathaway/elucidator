#pragma once

// Created by Nicholas Hathaway on 3/19/22.
/* 
    
*/

#include "elucidator/objects/seqContainers.h"
#include "elucidator/objects/BioDataObject.h"
#include "elucidator/seqToolsUtils/seqToolsUtils.hpp"


#include <njhseq/objects/seqObjects/seqKmers.h>
#include <njhseq/system.h>
#include <njhseq/BamToolsUtils.h>
#include <PathWeaver/objects/bam/RegionInvestigatorInBam.hpp>

namespace njhseq {



struct RefSeqsWithKmers{
	RefSeqsWithKmers(const bfs::path & refFnp, uint32_t reOrientingKmerLength): refFnp_(refFnp){

		SeqInput refReader(SeqIOOptions::genFastaIn(refFnp_));
		refSeqs_ = refReader.readAllReads<seqInfo>();
		refKmerReads_.reserve(refSeqs_.size());
		for (const auto & seq : refSeqs_) {
			refKmerReads_.emplace_back(std::make_shared<seqWithKmerInfo>(seq));
		}
		allSetKmers(refKmerReads_, reOrientingKmerLength, true);
		revComp_refSeqs_.reserve(refSeqs_.size());
		revComp_refSeqsKInfos_.reserve(refSeqs_.size());
		for(const auto & rSeq : refSeqs_){
			revComp_refSeqs_.emplace_back(rSeq);
			revComp_refSeqs_.back().reverseComplementRead(false, true);
			revComp_refSeqsKInfos_.emplace_back(revComp_refSeqs_.back().seq_, 7, false);
		}
		refSeqsKmerInfos_.reserve(refSeqs_.size());
		for(const auto & input : refSeqs_){
			refSeqsKmerInfos_.emplace_back(input.seq_, 7, false);
		}
	}
	bfs::path refFnp_;
	std::vector<seqInfo> refSeqs_;
	std::vector<std::shared_ptr<seqWithKmerInfo>> refKmerReads_;
	std::vector<kmerInfo> refSeqsKmerInfos_;
	std::vector<seqInfo> revComp_refSeqs_;
	std::vector<kmerInfo> revComp_refSeqsKInfos_;
};


inline std::vector<std::shared_ptr<seqWithKmerInfo>> trimToFinalSeqs(const std::vector<std::shared_ptr<seqWithKmerInfo>> & contigsKmerReads, const RefSeqsWithKmers & refSeqs){
	uint64_t maxLen = 0;
	readVec::getMaxLength(refSeqs.refSeqs_, maxLen);
	readVec::getMaxLength(contigsKmerReads, maxLen);
	aligner alignerObj(maxLen, gapScoringParameters(5,1,0,0,0,0), substituteMatrix(2,-2), false);
	std::vector<std::shared_ptr<seqWithKmerInfo>> finalSeqs;
	for(const auto & seq : contigsKmerReads){
		auto trimmed = readVecTrimmer::trimSeqToRefByGlobalAlnBestNoOverlapIncludeRevComp(seq, refSeqs.refSeqs_, refSeqs.revComp_refSeqs_, refSeqs.refSeqsKmerInfos_, refSeqs.revComp_refSeqsKInfos_, alignerObj, false);
		for(auto & trimmedSeq : trimmed){
			bool found = false;
			for(const auto & finalSeq : finalSeqs){
				if(finalSeq->seqBase_.seq_ == trimmedSeq.seq_){
					found = true;
					break;
				}
			}
			if(!found){
				finalSeqs.emplace_back(std::make_shared<seqWithKmerInfo>(trimmedSeq, 7, false));
			}
		}
	}
	return finalSeqs;
}


/*
struct FinalResultFiles{

	bfs::path finalDirectory = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar("final"));
	bfs::path partialDirectory = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar("partial"));
	auto allFinalSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(finalDirectory, "allFinal.fasta"));
	auto allPartialSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(partialDirectory, "allPartial.fasta"));
	SeqOutput allFinalWriter(allFinalSeqOpts);
	SeqOutput allPartialWriter(allPartialSeqOpts);
	allFinalWriter.openOut();
	allPartialWriter.openOut();
	std::mutex allFinalWriterMut;
	std::mutex allPartialWriterMut;
};
*/


struct ExtractCountsRaw : BamExtractor::ExtractCounts{

	double medianReadLength_ = 0;
};

inline ExtractCountsRaw rawWriteExtractReadsFromBamOnlyMapped(const bfs::path &bamFnp,
																											 const OutOptions &outOpts,
																											 bool keepRefStrainDir = false) {
	ExtractCountsRaw ret;
	BamTools::BamReader bReader;
	BamTools::BamAlignment bAln;
	bReader.Open(bamFnp.string());
	checkBamOpenThrow(bReader, bamFnp);
	//loadBamIndexThrow(bReader);
	auto refs = bReader.GetReferenceData();
	BamAlnsCache alnCache;
	auto refData = bReader.GetReferenceData();
	std::unordered_map<std::string, uint32_t> refNameToId;
	for (auto pos : iter::range(refData.size())) {
		refNameToId[refData[pos].RefName] = pos;
	}
	//pair writer
	SeqIOOptions outOptsPaired(outOpts.outFilename_,
														 SeqIOOptions::outFormats::FASTQPAIRED, outOpts);
	SeqOutput pairWriter(outOptsPaired);

	//non paired writer
	SeqIOOptions outOptsSingle(outOpts.outFilename_,
														 SeqIOOptions::outFormats::FASTQ, outOpts);
	SeqOutput writer(outOptsSingle);
	std::vector<uint64_t> readLens;
	while (bReader.GetNextAlignment(bAln)) {
		//skip secondary alignments
		if (!bAln.IsPrimaryAlignment()) {
			continue;
		}
		//won't handled sequences with pairs that have a unmmaped mate, will make them orphan reads
		//should improve upon
		if(!bAln.IsMapped()){
			continue;
		}
		/**@todo consider skipping duplicates */
		if (bAln.IsPaired()) {
			if (!alnCache.has(bAln.Name)) {
				//pair hasn't been added to cache yet so add to cache
				//this only works if mate and first mate have the same name
				alnCache.add(bAln);
				continue;
			} else {
				auto search = alnCache.get(bAln.Name);
				readLens.emplace_back(bAln.QueryBases.size());
				readLens.emplace_back(search->QueryBases.size());
				if (bAln.IsFirstMate()) {
					pairWriter.openWrite(
									PairedRead(bamAlnToSeqInfo(bAln, keepRefStrainDir), bamAlnToSeqInfo(*search, keepRefStrainDir),
														 false));
				} else {
					pairWriter.openWrite(
									PairedRead(bamAlnToSeqInfo(*search, keepRefStrainDir), bamAlnToSeqInfo(bAln, keepRefStrainDir),
														 false));
				}
				++ret.pairedReads_;
				// now that operations have been computed, remove first mate found from cache
				alnCache.remove(search->Name);
				continue;
			}
		} else {
			//unpaired read
			++ret.unpaiedReads_;
			if (bAln.IsMapped()) {
				// do unpaired read operation
				readLens.emplace_back(bAln.QueryBases.size());
				writer.openWrite(bamAlnToSeqInfo(bAln, keepRefStrainDir));
			}
		}
	}
	//save the orphans;
	if (len(alnCache) > 0) {
		auto names = alnCache.getNames();
		for (const auto & name : names) {
			++ret.orphans_;
			auto search = alnCache.get(name);
			readLens.emplace_back(search->QueryBases.size());
			writer.openWrite(bamAlnToSeqInfo(*search, keepRefStrainDir));
			alnCache.remove(name);
		}
	}
	ret.medianReadLength_ = vectorMedianRef(readLens);
	return ret;
}


struct FermiLiteNameParse{

	explicit FermiLiteNameParse(std::string fullname):fullname_(std::move(fullname)){
		std::smatch match;
		if(!std::regex_match(fullname_, match, fermiLiteNamePat_)){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error in processing " << fullname_ << " for basic assembly info" << "\n";
			throw std::runtime_error{ss.str()};
		}
		shortName_ = match[1];
		coverage_ =  njh::StrToNumConverter::stoToNum<uint32_t>(match[2]);
		modFullname_ = njh::pasteAsStr(match[1], " ", match[2], " ", match[3], " ", match[4]);
	}
	std::string fullname_;
	std::string modFullname_;
	double coverage_{0};
	std::string shortName_;
	std::regex fermiLiteNamePat_ = std::regex(R"((\d+:\d+)\s+(\d+)\s+([\d\.;,]+)\s+([\d\.;,]+)(.*))");
};

struct DefaultAssembleNameInfo{

//	DefaultAssembleNameInfo(const std::string & fullname):fullname_(fullname){
//		setInfoFromName();
//	}

	explicit DefaultAssembleNameInfo(std::string fullname, bool megahit = false):fullname_(std::move(fullname)){
		if(megahit){
			setInfoFromNameMegahit();
		}else{
			setInfoFromName();
		}
	}

	/**@brief set node name, length and coverage if all exists within name, must be in that group order in the regex pattern
	 *
	 * @param fullname
	 * @param pat
	 */
	DefaultAssembleNameInfo(std::string fullname, const std::regex & pat):fullname_(std::move(fullname)){
		std::smatch match;
		if(!std::regex_match(fullname_, match, pat)){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error in processing " << fullname_ << " for basic assembly info" << "\n";
			throw std::runtime_error{ss.str()};
		}
		nodeName_ = match[1];
		len_ =  njh::StrToNumConverter::stoToNum<uint32_t>(match[2]);
		coverage_ =  njh::StrToNumConverter::stoToNum<double>(match[3]);
	}

	std::string fullname_;

	void setInfoFromName(){
		std::smatch match;
		std::regex pat{R"((NODE_\d+)_length_(\d+)_cov_([0-9.]+).*)"};
		if(!std::regex_match(fullname_, match, pat)){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error in processing " << fullname_ << " for basic assembly info" << "\n";
			throw std::runtime_error{ss.str()};
		}
		nodeName_ = match[1];
		len_ =  njh::StrToNumConverter::stoToNum<uint32_t>(match[2]);
		coverage_ =  njh::StrToNumConverter::stoToNum<double>(match[3]);
	}
	//k99_0 flag=1 multi=94.5593 len=1027
	void setInfoFromNameMegahit(){
		std::smatch match;
		std::regex pat{R"((k[0-9]+_\d+) flag=\d+ multi=([0-9.]+) len=(\d+).*)"};
		if(!std::regex_match(fullname_, match, pat)){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error in processing " << fullname_ << " for basic assembly info" << "\n";
			throw std::runtime_error{ss.str()};
		}
		nodeName_ = match[1];

		coverage_ =  njh::StrToNumConverter::stoToNum<double>(match[2]);

		len_ =  njh::StrToNumConverter::stoToNum<uint32_t>(match[3]);
	}

	[[maybe_unused]] std::string nodeName_;
	uint32_t len_ {0};
	double coverage_ {0.0};

};

struct MIRAAssembleNameInfo{

//	DefaultAssembleNameInfo(const std::string & fullname):fullname_(fullname){
//		setInfoFromName();
//	}

	explicit MIRAAssembleNameInfo(std::string fullname):fullname_(std::move(fullname)){
		setInfoFromName();
	}

	std::string fullname_;

	void setInfoFromName(){
		std::smatch match;
		std::regex pat{R"(([A-z1-9_]+\d+)\s+cov=([0-9.]+)\s+len=(\d+)\s+gc=([0-9.]+)\s+nseq=(\d+).*)"};
		if(!std::regex_match(fullname_, match, pat)){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error in processing " << fullname_ << " for basic assembly info" << "\n";
			throw std::runtime_error{ss.str()};
		}
		nodeName_ = match[1];
		coverage_ =   njh::StrToNumConverter::stoToNum<double>(match[2]);
		len_ =        njh::StrToNumConverter::stoToNum<uint32_t>(match[3]);
		gcContent_ =  njh::StrToNumConverter::stoToNum<double>(match[4]);
		seqNumber_ =  njh::StrToNumConverter::stoToNum<uint32_t>(match[5]);

	}

	[[maybe_unused]] std::string nodeName_;
	[[maybe_unused]] uint32_t len_ {std::numeric_limits<uint32_t>::max()};
	double coverage_{std::numeric_limits<double>::max()};
	[[maybe_unused]] double gcContent_{std::numeric_limits<double>::max()};
	uint32_t seqNumber_ {std::numeric_limits<uint32_t>::max()};

};



}  // namespace njhseq

