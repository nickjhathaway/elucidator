/*
 * seqUtilsInfoRunner_getHapInfoRegion.cpp
 *
 *  Created on: Apr 25, 2019
 *      Author: nicholashathaway
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

#include "seqUtilsInfoRunner.hpp"
#include "elucidator/seqToolsUtils/seqToolsUtils.hpp"
#include "elucidator/objects/BioDataObject.h"
#include <njhseq/GenomeUtils.h>
#include <njhseq/objects/Gene/GeneFromGffs.hpp>

#include "elucidator/BioRecordsUtils/BedUtility.hpp"
#include "elucidator/PopulationGenetics.h"
#include <njhseq/objects/Gene/GenomicAminoAcidPositionTyper.hpp>
#include <njhseq.h>

namespace njhseq {



int seqUtilsInfoRunner::quickHaplotypeVariantsWithRegion(const njh::progutils::CmdArgs & inputCommands) {
	std::string identifier = "";
	std::string fstMeta = "";
	double lowVariantCutOff = 0.005;
	uint32_t occurrenceCutOff = 5;
	uint32_t lengthDiffForLengthPoly = 6;
	bfs::path bedFnp = "";
	bfs::path intersectBedFnp = "";

	bfs::path genomefnp = "";
	uint32_t outwardsExpand = 5;
	OutOptions outOpts(bfs::path(""));
	outOpts.outExtention_ = ".tab.txt";
	MultiGenomeMapper::inputParameters gPars;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processReadInNames(true);
	setUp.processWritingOptions(outOpts);
	setUp.processAlnInfoInput();
	setUp.setOption(lengthDiffForLengthPoly, "--lengthDiffForLengthPoly", "length Diff For Length Poly");
	setUp.setOption(occurrenceCutOff, "--occurrenceCutOff", "Occurrence Cut Off, don't report variants above this percentage");
	setUp.setOption(outwardsExpand, "--outwardsExpand", "The amount to expand outwards from given region when determining variants positions with extracted ref seq");
	setUp.setOption(intersectBedFnp,    "--intersectBed",    "A bed to get intersecting variants for", true);
	setUp.setOption(bedFnp,    "--bed",    "A bed file of the location for the extraction", true);
	setUp.setOption(genomefnp, "--genome", "A reference genome to compare against", true);
	gPars.genomeDir_ = genomefnp.parent_path();
	gPars.primaryGenome_ = bfs::basename(genomefnp);
	gPars.primaryGenome_ = gPars.primaryGenome_.substr(0, gPars.primaryGenome_.rfind("."));

	setUp.setOption(identifier, "--identifier", "Give a identifier name for info");
	setUp.setOption(fstMeta,    "--fstMeta",    "Meta field to calculate fsts by", true);
	setUp.finishSetUp(std::cout);


	njh::files::checkExistenceThrow(intersectBedFnp, __PRETTY_FUNCTION__);
	njh::files::checkExistenceThrow(bedFnp,          __PRETTY_FUNCTION__);
	njh::files::checkExistenceThrow(genomefnp,       __PRETTY_FUNCTION__);

	uint64_t maxLen = 0;
	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();

	std::unordered_map<std::string, std::shared_ptr<std::vector<identicalCluster>>> uniqueSeqsByMeta;
	std::vector<identicalCluster> clusters;
	seqInfo seq;
	uint32_t totalInput = 0;
	bool calculatingFst = "" != fstMeta;
	std::set<std::string> samples;
	std::set<std::string> allMetaKeys;
	std::unordered_map<std::string, std::set<std::string>> samplesByMeta;
	// read in reads and collapse to unique
	while(reader.readNextRead(seq)) {
		//get meta keys if available
		if(MetaDataInName::nameHasMetaData(getSeqBase(seq).name_)){
			MetaDataInName metaData(getSeqBase(seq).name_);
			for(const auto & meta : metaData.meta_){
				allMetaKeys.emplace(meta.first);
			}
		}

		readVec::handelLowerCaseBases(seq, setUp.pars_.ioOptions_.lowerCaseBases_);
		if(setUp.pars_.ioOptions_.removeGaps_){
			seq.removeGaps();
		}
		readVec::getMaxLength(seq, maxLen);
		if(seq.nameHasMetaData()){
			MetaDataInName seqMeta(seq.name_);
			if(seqMeta.containsMeta("sample")){
				samples.emplace(seqMeta.getMeta("sample"));
			}
		}
		++totalInput;
		bool found = false;
		for (auto &cIter : clusters) {
			if (cIter.seqBase_.seq_ == seq.seq_) {
				cIter.addRead(seq);
				found = true;
				break;
			}
		}
		if (!found) {
			clusters.emplace_back(seq);
		}
		//if calculating sub population differences collect unique sequences for the given field
		if(calculatingFst){
			MetaDataInName seqMeta(seq.name_);
			bool found = false;
			if(!njh::in(seqMeta.getMeta(fstMeta),  uniqueSeqsByMeta)){
				uniqueSeqsByMeta[seqMeta.getMeta(fstMeta)] = std::make_shared<std::vector<identicalCluster>>();
			}
			if(seqMeta.containsMeta("sample")){
				samplesByMeta[seqMeta.getMeta(fstMeta)].emplace(seqMeta.getMeta("sample"));
			}

			for (auto &cIter : *uniqueSeqsByMeta[seqMeta.getMeta(fstMeta)]) {
				if (cIter.seqBase_.seq_ == seq.seq_) {
					cIter.addRead(seq);
					found = true;
					break;
				}
			}
			if (!found) {
				uniqueSeqsByMeta[seqMeta.getMeta(fstMeta)]->emplace_back(seq);
			}
		}
	}

	VecStr allMetaKeysVec(allMetaKeys.begin(), allMetaKeys.end());

	uint32_t samplesCalled = totalInput;
	if(!samples.empty()){
		samplesCalled = samples.size();
	}
	auto beds = getBeds(bedFnp);
	if(beds.empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error no records found in " << bedFnp << "\n";
		throw std::runtime_error{ss.str()};
	}
	auto intersectingBeds = getBeds(intersectBedFnp);
	if(intersectingBeds.empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error no records found in " << intersectBedFnp << "\n";
		throw std::runtime_error{ss.str()};
	}

	MultiGenomeMapper gMapper(gPars);
	gMapper.loadInGenomes();
	gMapper.setUpGenomes();
	auto gPos = beds.front();

	//get identifer for sequences
	if("" == identifier){
		identifier = gPos->name_;
	}
	njh::sort(clusters);
	uint32_t seqId = 0;
	std::unordered_map<std::string, std::string> nameLookUp;
	for (auto &cIter : clusters) {
		MetaDataInName popMeta;
		popMeta.addMeta("HapPopUIDCount", static_cast<uint32_t>(std::round(cIter.seqBase_.cnt_)));
		cIter.seqBase_.name_ = njh::pasteAsStr(identifier, ".", leftPadNumStr<uint32_t>(seqId, clusters.size()),popMeta.createMetaName());
		nameLookUp[cIter.seqBase_.seq_] = cIter.seqBase_.name_;
		++seqId;
	}
	if(calculatingFst){
		std::unordered_map<std::string, uint32_t> seqIdForField;
		for(auto & field : uniqueSeqsByMeta){
			readVec::allSetFractionByTotalCount(*field.second);
			for (auto &cIter : *field.second) {
				MetaDataInName popMeta;
				std::string popName = nameLookUp[cIter.seqBase_.seq_];
				MetaDataInName popNameMeta(popName);
				MetaDataInName::removeMetaDataInName(popName);
				popMeta.addMeta("HapPopUID", popName);
				popMeta.addMeta(popNameMeta, true);
				popMeta.addMeta("HapUniqueCount", static_cast<uint32_t>(std::round(cIter.seqBase_.cnt_)));
				cIter.seqBase_.name_ = njh::pasteAsStr(field.first, ".", leftPadNumStr<uint32_t>(seqIdForField[field.first], field.second->size()), popMeta.createMetaName());
				++seqIdForField[field.first];
			}
		}
	}


	//get region and ref seq for mapping of variants
	BedUtility::extendLeftRight(*gPos, outwardsExpand, outwardsExpand,
			gMapper.genomes_.at(gPars.primaryGenome_)->chromosomeLengths_.at(
					gPos->chrom_));
	GenomicRegion refRegion(*gPos);

	TwoBit::TwoBitFile tReader(
			gMapper.genomes_.at(gPars.primaryGenome_)->fnpTwoBit_);
	auto refSeq = refRegion.extractSeq(tReader);

	readVec::getMaxLength(refSeq, maxLen);
	aligner alignerObj(maxLen, gapScoringParameters(5,1,0,0,0,0), substituteMatrix(2,-2), false);
	alignerObj.weighHomopolymers_ =false;
	alignerObj.processAlnInfoInput(setUp.pars_.alnInfoDirName_, setUp.pars_.verbose_);

	table outputTab(VecStr{"CHROM", "POS", "ID", "REF", "ALT", fstMeta, "count", "total"});
	if(calculatingFst){
		for(auto & field : uniqueSeqsByMeta){
			njh::sort(*field.second);
			double seqCount = readVec::getTotalReadCount(*field.second);

			std::unordered_map<uint32_t, std::unordered_map<char,        uint32_t>> snps;
			std::unordered_map<uint32_t, std::unordered_map<std::string, uint32_t>> insertions;
			std::unordered_map<uint32_t, std::unordered_map<std::string, uint32_t>> deletions;

			std::unordered_map<uint32_t, std::unordered_map<char,        uint32_t>> snpsFinal;
			std::unordered_map<uint32_t, std::unordered_map<std::string, uint32_t>> insertionsFinal;
			std::unordered_map<uint32_t, std::unordered_map<std::string, uint32_t>> deletionsFinal;

			for(const auto & seq : *field.second){
				alignerObj.alignCacheGlobal(refSeq, seq);
				alignerObj.profileAlignment(refSeq, seq, false, false, false);
				for(const auto & m : alignerObj.comp_.distances_.mismatches_){
					snps[m.second.refBasePos][m.second.seqBase]+= seq.seqBase_.cnt_;
				}
				for(const auto & g : alignerObj.comp_.distances_.alignmentGaps_){
					if(g.second.ref_){
						//insertion
						insertions[g.second.refPos_ - 1][g.second.gapedSequence_]++;
					}else{
						//deletion
						deletions[g.second.refPos_ - 1][g.second.gapedSequence_]++;
					}
				}
			}
			std::vector<uint32_t> variablePositions;
			//filter snps and indels by occurrence cut off
			for(const auto & snp : snps){
				for(const auto & b : snp.second){
					if(b.second < occurrenceCutOff){
						continue;
					}
					snpsFinal[snp.first][b.first] = b.second;
					if(b.second/static_cast<double>(seqCount) > lowVariantCutOff){
						variablePositions.emplace_back(snp.first);
					}
				}
			}
			for(const auto & del : deletions){
				for(const auto & d : del.second){
					if(d.second < occurrenceCutOff){
						continue;
					}
					deletionsFinal[del.first][d.first] = d.second;
					if(d.second/static_cast<double>(seqCount) > lowVariantCutOff){
						variablePositions.emplace_back(del.first);
					}
				}
			}
			for(const auto & ins : insertions){
				for(const auto & i : ins.second){
					if(i.second < occurrenceCutOff){
						continue;
					}
					insertionsFinal[ins.first][i.first] = i.second;
					if(i.second/static_cast<double>(seqCount) > lowVariantCutOff){
						variablePositions.emplace_back(ins.first);
					}
				}
			}
			if(!snpsFinal.empty() || ! deletionsFinal.empty() || !insertionsFinal.empty()){

				std::unordered_set<uint32_t> positionsSet;
				for(const auto & snps : snpsFinal){
					positionsSet.emplace(snps.first);
				}
				for(const auto & ins : insertionsFinal){
					positionsSet.emplace(ins.first);
				}
				for(const auto & del : deletionsFinal){
					positionsSet.emplace(del.first);
				}
				std::vector<uint32_t> positions(positionsSet.begin(), positionsSet.end());
				uint32_t samplesCalledForField = seqCount;
				if(!samplesByMeta[field.first].empty()){
					samplesCalledForField = samplesByMeta[field.first].size();
				}
				njh::sort(positions);
				if(refRegion.reverseSrand_){
					njh::reverse(positions);
				}
				for(const auto & pos : positions){
					auto relativePos = refRegion.getRelativePositionFromStartStrandAware(pos);
					if (intersectingBeds.front()->chrom_ == refRegion.chrom_
							&& relativePos < intersectingBeds.front()->chromEnd_
							&& relativePos >= intersectingBeds.front()->chromStart_) {
						//	table outputTab(VecStr{"CHROM", "POS", "ID", "REF", "ALT", fstMeta, "count", "total"});
						//intersectingBeds.front()->chrom_,relativePos,identifier,
						std::string ref = njh::pasteAsStr(refSeq.seq_[pos]);
						if(refRegion.reverseSrand_){
							ref = seqUtil::reverseComplement(std::string(1, refSeq.seq_[pos]), "DNA");
						}
						if(njh::in(pos, snpsFinal)){
							for(const auto & b : snpsFinal[pos]){
								outputTab.addRow(intersectingBeds.front()->chrom_,relativePos,identifier, ref, b.first, field.first, b.second, seqCount);
							}
						}
					}
				}
			}
		}
	}
	if(outputTab.nRow() > 0){
		OutputStream out(outOpts);
		outputTab.outPutContents(out, "\t");
	}
	return 0;
}








}

