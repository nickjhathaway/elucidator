/*
 * genExp_bowtie2ExtractAndCompare.cpp
 *
 *  Created on: Feb 3, 2018
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

#include <TwoBit.h>
#include "genExp.hpp"
#include "elucidator/BamToolsUtils.h"
#include "elucidator/objects/BioDataObject.h"
#include "elucidator/objects/dataContainers/graphs.h"
#include <njhseq/GenomeUtils.h>
#include "elucidator/concurrency/LockableJsonLog.hpp"



namespace njhseq {

int genExpRunner::bowtie2ExtractAndCompareMultiple(const njh::progutils::CmdArgs & inputCommands){
	std::string genomes = "";
	bfs::path genomesDir = "";
	uint32_t numThreads = 1;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(numThreads, "--numThreads", "number of cpus to use");
  setUp.setOption(genomes, "--genomes", "Names of the genomes to extract from (should not have extension, e.g. Pf3d7 for Pf3d7.fasta", true);
  setUp.setOption(genomesDir, "--genomeDir", "Names of the genome directory where the genomes are stored", true);
  setUp.processReadInNames();
  setUp.processDirectoryOutputName(true);
	setUp.finishSetUp(std::cout);

	setUp.startARunLog(setUp.pars_.directoryName_);
	BioCmdsUtils bRunner(setUp.pars_.verbose_);

	MultiGenomeMapper gMapper(genomesDir, genomes.substr(0, genomes.find(",")));
	gMapper.pars_.numThreads_ = numThreads;
	gMapper.setSelectedGenomes(njh::strToSet<std::string>(genomes, ","));
	gMapper.loadInGenomes();
	gMapper.setUpGenomes();

	std::map<std::string, std::unordered_map<std::string, std::vector<std::shared_ptr<AlignmentResults>>>> allAlnResults;
	std::unordered_map<std::string, uint32_t> unmappedCounts;
	for(const auto & genome : gMapper.genomes_){
		bfs::path genomeDir = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar{genome.first});
		//align seqs
		auto seqOpts = setUp.pars_.ioOptions_;
		seqOpts.out_.outFilename_ = njh::files::make_path(genomeDir, "alignedSeqs.sorted.bam");
		auto runOut = bRunner.bowtie2Align(seqOpts, genome.second->fnp_, "-a");
		OutOptions bowtie2AlignLogOpts(njh::files::make_path(genomeDir, "bowtie2Log.json"));
		OutputStream bowtie2AlignLogOut(bowtie2AlignLogOpts);
		bowtie2AlignLogOut << njh::json::toJson(runOut) << std::endl;
		//extract locations and mapping stats

		BamTools::BamAlignment bAln;
		std::map<std::string, std::vector<BamTools::BamAlignment>> bamAligns;
		std::vector<BamTools::BamAlignment> unmapped;
		BamTools::BamReader bReader;
		bReader.Open(seqOpts.out_.outFilename_.string());
		checkBamOpenThrow(bReader, seqOpts.out_.outFilename_.string());
		auto refData = bReader.GetReferenceData();

		while (bReader.GetNextAlignment(bAln)) {
			if (bAln.IsMapped()) {
				bamAligns[bAln.Name].emplace_back(bAln);
			} else {
				unmapped.emplace_back(bAln);
				++unmappedCounts[bAln.Name];
			}
		}
		TwoBit::TwoBitFile twobitReader(genome.second->fnpTwoBit_);
		std::vector<GenomicRegion> alignedRegions;
		std::map<uint32_t, uint32_t> mapCounts;
		mapCounts[0] = unmapped.size();

		auto regionsExtractedOpts = SeqIOOptions::genFastaOut(njh::files::make_path(genomeDir, "regions"));
		SeqOutput extractedWriter(regionsExtractedOpts);
		extractedWriter.openOut();

		//alignment comparisons
		OutOptions comparisonOpts(njh::files::make_path(genomeDir, "refComparisonInfo.tab.txt"));
		OutputStream comparisonOut(comparisonOpts);

		comparisonOut << "ReadNumber\tReadId\tBestRef\tscore\talnScore\thqScore\tkDist-5\t1bIndel\t2bIndel\t>2bIndel\tlqMismatch\thqMismatch" << std::endl;
		uint32_t readNumber = 0;
		std::vector<std::shared_ptr<seqInfo>> refSeqs;
		std::unordered_map<std::string, VecStr> readNamesToRefSeqs;
		for (const auto & alnForRead : bamAligns) {
			++mapCounts[alnForRead.second.size()];
			uint32_t extractionCount = 0;
			for (const auto & aln : alnForRead.second) {
				auto results = std::make_shared<AlignmentResults>(aln, refData);
				results->setRefSeq(twobitReader);
				MetaDataInName refMeta;
				refMeta.addMeta("genome", genome.first);
				refMeta.addMeta("chrom", results->gRegion_.chrom_);
				refMeta.addMeta("start", results->gRegion_.start_);
				refMeta.addMeta("end", results->gRegion_.end_);
				results->refSeq_->name_.append(refMeta.createMetaName());
				if(!njh::in(results->refSeq_->name_, readNamesToRefSeqs)){
					refSeqs.emplace_back(std::make_shared<seqInfo>(*(results->refSeq_)));
				}
				readNamesToRefSeqs[results->refSeq_->name_].emplace_back(aln.Name);
				kmerInfo refInfo(results->refSeq_->seq_, 5, false);
				kmerInfo seqKInfo(results->alnSeq_->seq_, 5, false);
				results->setComparison(false);
				alignedRegions.emplace_back(results->gRegion_);
				comparisonOut << readNumber
						<< '\t' << aln.Name
						<< '\t' << results->gRegion_.createUidFromCoordsStrand()
						<< '\t' << results->comp_.distances_.eventBasedIdentityHq_
						<< '\t' << results->comp_.alnScore_
						<< '\t' << refInfo.compareKmers(seqKInfo).second
						<< '\t' << results->comp_.oneBaseIndel_
						<< '\t' << results->comp_.twoBaseIndel_
						<< '\t' << results->comp_.largeBaseIndel_
						<< '\t' << results->comp_.lqMismatches_
						<< '\t' << results->comp_.hqMismatches_ << std::endl;
				allAlnResults[aln.Name][genome.first].emplace_back(results);
				++extractionCount;
			}
			++readNumber;
		}
		//write out ref seqs;
		//read names for ref seqs
		OutOptions readNamesOpts(njh::files::make_path(genomeDir, "readNamesForRefSeqs.tab.txt"));
		OutputStream readNamesOut(readNamesOpts);
		readNamesOut << "refName\treadNames" << std::endl;
		VecStr refNames = njh::getVecOfMapKeys(readNamesToRefSeqs);
		njh::sort(refNames, [&readNamesToRefSeqs](const std::string & ref1, const std::string & ref2){
			return readNamesToRefSeqs[ref1].size() == readNamesToRefSeqs[ref2].size() ? ref1 < ref2 :  readNamesToRefSeqs[ref1].size() > readNamesToRefSeqs[ref2].size();
		});
		for(const auto & refName : refNames){
			readNamesOut << refName
					<< "\t" << njh::conToStr(readNamesToRefSeqs[refName], ",") << std::endl;
		}
		//collapse refseqs
		for(auto & refSeq : refSeqs){
			refSeq->cnt_ = readNamesToRefSeqs[refSeq->name_].size();
			MetaDataInName refMeta(refSeq->name_);
			refMeta.addMeta("extractCount", refSeq->cnt_);
			refMeta.resetMetaInName(refSeq->name_);
		}
		njh::sort(refSeqs, [](const std::shared_ptr<seqInfo> & ref1, const std::shared_ptr<seqInfo> & ref2){
			return ref1->cnt_ == ref2->cnt_ ? ref1->name_ < ref2->name_ : ref1->cnt_ > ref2->cnt_;
		});
		extractedWriter.write(refSeqs);


		//genomic regions hit
		OutOptions regionsOpts(njh::files::make_path(genomeDir, "regions.bed"));
		OutputStream regionsOut(regionsOpts);
		for(const auto & reg : alignedRegions){
			regionsOut << reg.genBedRecordCore().toDelimStr() << std::endl;
		}

		//map counts
		OutOptions mapCountsOpts(njh::files::make_path(genomeDir, "mapCounts.tab.txt"));
		OutputStream mapCountsOut(mapCountsOpts);
		table mapCountsTab(mapCounts, VecStr{"hits", "total"});
		mapCountsTab.outPutContents(mapCountsOut, "\t");

		//names of unmapped sequences
		OutOptions unmmapedOpts(njh::files::make_path(genomeDir, "unmappedReads.txt"));
		OutputStream unmappedOut(unmmapedOpts);
		for(const auto & unmappedAln : unmapped){
			unmappedOut << unmappedAln.Name << std::endl;
		}
	}


	//get best hits only
	auto regionsExtractedOpts = SeqIOOptions::genFastaOut(njh::files::make_path(setUp.pars_.directoryName_, "bestRegions"));
	SeqOutput extractedWriter(regionsExtractedOpts);
	extractedWriter.openOut();
	std::vector<std::shared_ptr<seqInfo>> refSeqs;
	std::unordered_map<std::string, VecStr> readNamesToRefSeqs;

	//alignment comparisons
	OutOptions comparisonOpts(njh::files::make_path(setUp.pars_.directoryName_, "refComparisonInfo.tab.txt"));
	OutputStream comparisonOut(comparisonOpts);

	comparisonOut << "ReadNumber\tReadId\tBestRef\tscore\talnScore\tkDist-5\t1bIndel\t2bIndel\t>2bIndel\tlqMismatch\thqMismatch" << std::endl;
	uint32_t readNumber = 0;
	for(const auto & alnResults : allAlnResults){
		double bestScore = std::numeric_limits<double>::lowest();
		std::vector<std::shared_ptr<AlignmentResults>> bestResults;
		for(const auto & genomeRes : alnResults.second){
			for(const auto & res : genomeRes.second){
				if(res->comp_.alnScore_ > bestScore){
					bestScore = res->comp_.alnScore_;
					bestResults.clear();
					bestResults.emplace_back(res);
				}else if(res->comp_.alnScore_  == bestScore){
					bestResults.emplace_back(res);
				}
			}
		}
		for(const auto & results : bestResults){
			if(!njh::in(results->refSeq_->name_, readNamesToRefSeqs)){
				refSeqs.emplace_back(results->refSeq_);
			}
			readNamesToRefSeqs[results->refSeq_->name_].emplace_back(results->bAln_.Name);

			kmerInfo refInfo(results->refSeq_->seq_, 5, false);
			kmerInfo seqKInfo(results->alnSeq_->seq_, 5, false);
			comparisonOut << readNumber
					<< '\t' << results->bAln_.Name
					<< '\t' << results->gRegion_.createUidFromCoordsStrand()
					<< '\t' << results->comp_.distances_.eventBasedIdentityHq_
					<< '\t' << results->comp_.alnScore_
					<< '\t' << refInfo.compareKmers(seqKInfo).second
					<< '\t' << results->comp_.oneBaseIndel_
					<< '\t' << results->comp_.twoBaseIndel_
					<< '\t' << results->comp_.largeBaseIndel_
					<< '\t' << results->comp_.lqMismatches_
					<< '\t' << results->comp_.hqMismatches_ << std::endl;
		}
		++readNumber;
	}

	//write out ref seqs;
	//read names for ref seqs
	OutOptions readNamesOpts(
			njh::files::make_path(setUp.pars_.directoryName_, "readNamesForRefSeqs.tab.txt"));
	OutputStream readNamesOut(readNamesOpts);
	readNamesOut << "refName\treadNames" << std::endl;
	VecStr refNames = njh::getVecOfMapKeys(readNamesToRefSeqs);
	njh::sort(refNames,
			[&readNamesToRefSeqs](const std::string & ref1, const std::string & ref2) {
				return readNamesToRefSeqs[ref1].size() == readNamesToRefSeqs[ref2].size() ? ref1 < ref2 : readNamesToRefSeqs[ref1].size() > readNamesToRefSeqs[ref2].size();
			});
	for (const auto & refName : refNames) {
		readNamesOut << refName << "\t"
				<< njh::conToStr(readNamesToRefSeqs[refName], ",") << std::endl;
	}
	//collapse refseqs
	for (auto & refSeq : refSeqs) {
		refSeq->cnt_ = readNamesToRefSeqs[refSeq->name_].size();
		MetaDataInName refMeta(refSeq->name_);
		refMeta.addMeta("extractCount", refSeq->cnt_, true);
		refMeta.resetMetaInName(refSeq->name_);
	}
	njh::sort(refSeqs,
			[](const std::shared_ptr<seqInfo> & ref1, const std::shared_ptr<seqInfo> & ref2) {
				return ref1->cnt_ == ref2->cnt_ ? ref1->name_ < ref2->name_ : ref1->cnt_ > ref2->cnt_;
			});
	extractedWriter.write(refSeqs);

	//names of unmapped sequences
	VecStr unmappedToAllGenomes;
	for(const auto & unmappedCount : unmappedCounts){
		if(gMapper.genomes_.size() == unmappedCount.second ){
			unmappedToAllGenomes.emplace_back(unmappedCount.first);
			comparisonOut << readNumber
					<< '\t' << unmappedCount.first
					<< '\t' << "*"
					<< '\t' << "*"
					<< '\t' << "*"
					<< '\t' << "*"
					<< '\t' << "*"
					<< '\t' << "*"
					<< '\t' << "*"
					<< '\t' << "*"
					<< '\t' << "*" << std::endl;
			++readNumber;
		}
	}
	if(!unmappedToAllGenomes.empty()){
		OutOptions unmmapedOpts(njh::files::make_path(setUp.pars_.directoryName_, "unmappedReads.txt"));
		OutputStream unmappedOut(unmmapedOpts);
		for(const auto & unmappedAln : unmappedToAllGenomes){
			unmappedOut << unmappedAln << std::endl;
		}
	}




	return 0;
}


int genExpRunner::bowtie2ExtractAndCompare(const njh::progutils::CmdArgs & inputCommands){
	bfs::path genomeFnp = "";
	uint32_t numThreads = 1;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(numThreads, "--numThreads", "number of cpus to use");
  setUp.setOption(genomeFnp, "--genomeFnp", "Name of the genome to extract from", true);
  setUp.processReadInNames();
  setUp.processDirectoryOutputName(true);
	setUp.finishSetUp(std::cout);

	setUp.startARunLog(setUp.pars_.directoryName_);
	BioCmdsUtils bRunner(setUp.pars_.verbose_);
	// index for alignment and for genome extract
	bRunner.RunBowtie2Index(genomeFnp);
	bRunner.RunFaToTwoBit(genomeFnp);

	//align seqs
	auto seqOpts = setUp.pars_.ioOptions_;
	seqOpts.out_.outFilename_ = njh::files::make_path(setUp.pars_.directoryName_, "alignedSeqs.sorted.bam");
	auto runOut = bRunner.bowtie2Align(seqOpts, genomeFnp);
	OutOptions bowtie2AlignLogOpts(njh::files::make_path(setUp.pars_.directoryName_, "bowtie2Log.json"));
	OutputStream bowtie2AlignLogOut(bowtie2AlignLogOpts);
	bowtie2AlignLogOut << njh::json::toJson(runOut) << std::endl;
	//extract locations and mapping stats

	BamTools::BamAlignment bAln;

	std::unordered_map<std::string, std::vector<BamTools::BamAlignment>> bamAligns;
	std::vector<BamTools::BamAlignment> unmapped;
	BamTools::BamReader bReader;
	bReader.Open(seqOpts.out_.outFilename_.string());
	checkBamOpenThrow(bReader, seqOpts.out_.outFilename_.string());
	auto refData = bReader.GetReferenceData();

	while (bReader.GetNextAlignment(bAln)) {
		if (bAln.IsMapped()) {
			bamAligns[bAln.Name].emplace_back(bAln);
		} else {
			unmapped.emplace_back(bAln);
		}
	}

	auto twoBitfnp = bfs::path(genomeFnp).replace_extension("").string() + ".2bit";
	TwoBit::TwoBitFile twobitReader(twoBitfnp);
	std::vector<GenomicRegion> alignedRegions;
	std::map<uint32_t, uint32_t> mapCounts;
	mapCounts[0] = unmapped.size();

	auto regionsExtractedOpts = SeqIOOptions::genFastaOut(njh::files::make_path(setUp.pars_.directoryName_, "regions"));
	SeqOutput extractedWriter(regionsExtractedOpts);
	extractedWriter.openOut();

	//alignment comparisons
	OutOptions comparisonOpts(njh::files::make_path(setUp.pars_.directoryName_, "refComparisonInfo.tab.txt"));
	OutputStream comparisonOut(comparisonOpts);

	comparisonOut << "ReadNumber\tReadId\tBestRef\tscore\talnScore\tkDist-5\t1bIndel\t2bIndel\t>2bIndel\tlqMismatch\thqMismatch" << std::endl;
	uint32_t readNumber = 0;
	std::vector<std::shared_ptr<seqInfo>> refSeqs;
	std::unordered_map<std::string, VecStr> readNamesToRefSeqs;
	for (const auto & alnForRead : bamAligns) {
		++mapCounts[alnForRead.second.size()];
		uint32_t extractionCount = 0;
		for (const auto & aln : alnForRead.second) {
			auto results = std::make_shared<AlignmentResults>(aln, refData);
			results->setRefSeq(twobitReader);
			MetaDataInName refMeta;
			refMeta.addMeta("genome", bfs::path(genomeFnp).filename().replace_extension("").string());
			refMeta.addMeta("chrom", results->gRegion_.chrom_);
			refMeta.addMeta("start", results->gRegion_.start_);
			refMeta.addMeta("end", results->gRegion_.end_);
			results->refSeq_->name_.append(refMeta.createMetaName());
			if(!njh::in(results->refSeq_->name_, readNamesToRefSeqs)){
				refSeqs.emplace_back(results->refSeq_);
			}
			readNamesToRefSeqs[results->refSeq_->name_].emplace_back(aln.Name);
			kmerInfo refInfo(results->refSeq_->seq_, 5, false);
			kmerInfo seqKInfo(results->alnSeq_->seq_, 5, false);
			results->setComparison(false);
			alignedRegions.emplace_back(results->gRegion_);
			comparisonOut << readNumber
					<< '\t' << aln.Name
					<< '\t' << results->gRegion_.createUidFromCoordsStrand()
					<< '\t' << results->comp_.distances_.eventBasedIdentityHq_
					<< '\t' << results->comp_.alnScore_
					<< '\t' << refInfo.compareKmers(seqKInfo).second
					<< '\t' << results->comp_.oneBaseIndel_
					<< '\t' << results->comp_.twoBaseIndel_
					<< '\t' << results->comp_.largeBaseIndel_
					<< '\t' << results->comp_.lqMismatches_
					<< '\t' << results->comp_.hqMismatches_ << std::endl;
			++extractionCount;
		}
		++readNumber;
	}
	//write out ref seqs;
	//read names for ref seqs
	OutOptions readNamesOpts(njh::files::make_path(setUp.pars_.directoryName_, "readNamesForRefSeqs.tab.txt"));
	OutputStream readNamesOut(readNamesOpts);
	readNamesOut << "refName\treadNames" << std::endl;
	VecStr refNames = njh::getVecOfMapKeys(readNamesToRefSeqs);
	njh::sort(refNames, [&readNamesToRefSeqs](const std::string & ref1, const std::string & ref2){
		return readNamesToRefSeqs[ref1].size() == readNamesToRefSeqs[ref2].size() ? ref1 < ref2 :  readNamesToRefSeqs[ref1].size() > readNamesToRefSeqs[ref2].size();
	});
	for(const auto & refName : refNames){
		readNamesOut << refName
				<< "\t" << njh::conToStr(readNamesToRefSeqs[refName], ",") << std::endl;
	}
	//collapse refseqs
	for(auto & refSeq : refSeqs){
		refSeq->cnt_ = readNamesToRefSeqs[refSeq->name_].size();
		MetaDataInName refMeta(refSeq->name_);
		refMeta.addMeta("extractCount", refSeq->cnt_);
		refMeta.resetMetaInName(refSeq->name_);
	}
	njh::sort(refSeqs, [](const std::shared_ptr<seqInfo> & ref1, const std::shared_ptr<seqInfo> & ref2){
		return ref1->cnt_ == ref2->cnt_ ? ref1->name_ < ref2->name_ : ref1->cnt_ > ref2->cnt_;
	});
	extractedWriter.write(refSeqs);

	//genomic regions hit
	OutOptions regionsOpts(njh::files::make_path(setUp.pars_.directoryName_, "regions.bed"));
	OutputStream regionsOut(regionsOpts);
	for(const auto & reg : alignedRegions){
		regionsOut << reg.genBedRecordCore().toDelimStr() << std::endl;
	}

	//map counts
	OutOptions mapCountsOpts(njh::files::make_path(setUp.pars_.directoryName_, "mapCounts.tab.txt"));
	OutputStream mapCountsOut(mapCountsOpts);
	table mapCountsTab(mapCounts, VecStr{"hits", "total"});
	mapCountsTab.outPutContents(mapCountsOut, "\t");

	//names of unmapped sequences
	OutOptions unmmapedOpts(njh::files::make_path(setUp.pars_.directoryName_, "unmappedReads.txt"));
	OutputStream unmappedOut(unmmapedOpts);
	for(const auto & unmappedAln : unmapped){
		unmappedOut << unmappedAln.Name << std::endl;
	}


	return 0;
}

}  // namespace njhseq
