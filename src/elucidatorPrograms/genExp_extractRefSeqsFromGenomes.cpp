/*
 * genExp_extractRefSeqsFromGenomes.cpp
 *
 *  Created on: Jul 8, 2017
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


#include <SeekDeep/utils.h>

namespace njhseq {

int genExpRunner::bioIndexGenomes(const njh::progutils::CmdArgs & inputCommands){
	bfs::path genomeDir = "";
	std::string primaryGenome = "";
	std::string selectedGenomes = "";
	uint32_t numThreads = 1;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(numThreads, "--numThreads", "number of cpus to use");
  setUp.setOption(genomeDir, "--genomeDir", "Name of the genome file fnp", true);
	setUp.setOption(primaryGenome, "--primaryGenome", "Name of the primary genome, should be in --genomeDir");
	bfs::path genomeFnp = njh::files::make_path(genomeDir, njh::appendAsNeededRet(primaryGenome, ".fasta"));
	setUp.setOption(selectedGenomes, "--selectedGenomes", "Name of the other genomes in --genomeDir to be read in, leave blank to just do all fastas");
	setUp.finishSetUp(std::cout);
	auto fastaFiles = njh::files::gatherFiles(genomeDir, ".fasta", false);
	if(fastaFiles.empty()){
		std::cerr << "Found no .fasta files in " << genomeDir << std::endl;
		return 1;
	}
	std::unique_ptr<MultiGenomeMapper> gMapper;
	//set up genome mapper;
	gMapper = std::make_unique<MultiGenomeMapper>(genomeDir, bfs::basename(fastaFiles.front()));

	//set threads;
	gMapper->pars_.numThreads_ = numThreads;
	gMapper->pars_.verbose_ = setUp.pars_.verbose_;
	//set up selected genomes
	gMapper->setSelectedGenomes(selectedGenomes);
	gMapper->loadInGenomes();
	gMapper->bioIndexAllGenomes();

	return 0;
}

int genExpRunner::bioIndexGenome(const njh::progutils::CmdArgs & inputCommands){
	bfs::path genomeFnp = "";
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
  setUp.setOption(genomeFnp, "--genomeFnp", "the file name path to the genome", true);
  setUp.finishSetUp(std::cout);

  BioCmdsUtils bioRunner(setUp.pars_.verbose_);

  bioRunner.runAllPossibleIndexes(genomeFnp);

	return 0;
}

int genExpRunner::extractRefSeqsFromGenomesWithPrimers(
		const njh::progutils::CmdArgs & inputCommands) {
	extractBetweenSeqsPars pars;
	std::set<std::string> forwardPrimers;
	std::set<std::string> reversePrimers;
	std::string targetName = "";
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	pars.setUpCoreOptions(setUp, false);
	bool setPrimersFile = setUp.setOption(pars.primersFile, "--primers", "A file that contains three columns, targetName,forwardPrimer,reversePrimer 5` to 3` directions");
	if(!setPrimersFile){
		setUp.setOption(forwardPrimers, "--forwardPrimer", "Forward Primer (5`to 3` direction)", true);
		setUp.setOption(reversePrimers, "--reversePrimer", "Reverse Primer (5`to 3` direction)", true);
		setUp.setOption(targetName, "--targetName", "Name of the target associated with forward and reverse primers", true);
	}
	//setUp.processDirectoryOutputName("extractedRegions_TODAY", true);
	setUp.finishSetUp(std::cout);

	table primerTable;

	if ("" != pars.primersFile) {
		primerTable = seqUtil::readPrimers(pars.primersFile.string(), "whitespace",
				false);
	} else {
		primerTable = table(
				VecStr { "targetName", "forwardPrimer", "reversePrimer" });
		primerTable.addRow(targetName, njh::conToStr(forwardPrimers, ","), njh::conToStr(reversePrimers, ","));
	}

	std::unordered_map<std::string, PrimersAndMids::Target> targets;

	for(const auto & row : primerTable){
		if(njh::in(row[0], targets)){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error already have target: " << row[0] << "\n";
			throw std::runtime_error{ss.str()};
		}
		targets.emplace(row[0], PrimersAndMids::Target(row[0], row[1], row[2]));
	}

	PrimersAndMids ids(targets);
	if(0 == ids.getTargets().size() ){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error in reading in target primers file " << pars.primersFile << "\n";
		ss << "Make sure there is a line that starts with target in file" << "\n";
		throw std::runtime_error{ss.str()};
	}
	ids.initPrimerDeterminator();

	extractBetweenSeqs(ids, pars);

	setUp.startARunLog(pars.outputDirPars.dirName_.string());

	return 0;
}




std::vector<GenomicRegion> mergeRegionsStrandAware(const std::vector<GenomicRegion> & regions){
	std::vector<GenomicRegion> ret;
	std::unordered_map<bool, std::vector<GenomicRegion>> regionsbyStrand;
	for(const auto & region : regions){
		regionsbyStrand[region.reverseSrand_].emplace_back(region);
	}
	std::vector<GenomicRegion> fwdRegions;
	sortGRegionsByStart(regionsbyStrand[false]);
	for(const auto & reg : regionsbyStrand[false]){
		if(fwdRegions.empty()){
			fwdRegions.emplace_back(reg);
		}else{
			if(reg.start_ >= fwdRegions.back().start_ && reg.start_ <= fwdRegions.back().end_){
				fwdRegions.back().end_ = reg.end_;
			}else{
				fwdRegions.emplace_back(reg);
			}
		}
	}
	std::vector<GenomicRegion> revRegions;
	sortGRegionsByStart(regionsbyStrand[true]);
	for(const auto & reg : regionsbyStrand[true]){
		if(revRegions.empty()){
			revRegions.emplace_back(reg);
		}else{
			if(reg.start_ >= revRegions.back().start_ && reg.start_ <= revRegions.back().end_){
				revRegions.back().end_ = reg.end_;
			}else{
				revRegions.emplace_back(reg);
			}
		}
	}
	addOtherVec(ret, fwdRegions);
	addOtherVec(ret, revRegions);
	sortGRegionsByStart(ret);
	return ret;
}


int genExpRunner::extractRefSeqsFromGenomes(
		const njh::progutils::CmdArgs & inputCommands) {

	MultiGenomeMapper::inputParameters genomeMappingPars;
	bfs::path bedFile = "";
	bfs::path outputDir = "";
	bool keepRefAlignments = false;
	bool overWriteDirs = false;
	bool combineByGenome = false;
	MultiGenomeMapper::getRefSeqsWithPrimaryGenomePars extractPars;
	extractPars.lzPars.identity = 80;
	extractPars.lzPars.coverage = 90;

	bool writeOutAllSeqsFile = false;

	bfs::path gffDir = "";
	std::string gffExtraAttributesStr = "description";
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(extractPars.shortNames, "--shortNames", "Short Names");
	setUp.setOption(writeOutAllSeqsFile, "--writeOutAllSeqsFile", "Write Out All Seqs File");

	setUp.setOption(combineByGenome, "--combineByGenome", "Combine and merge regions extracted and create a bed file for each genome");
	setUp.setOption(extractPars.extendAndTrim, "--extendAndTrim", "Extend the determine region and then trim back, can be helpful for when variation falls at the very ends of the sequence");
	setUp.setOption(extractPars.extendAndTrimLen, "--extendAndTrimLen", "When extending and trimming, use this length");
	setUp.setOption(extractPars.lzPars.identity, "--identity", "Identity to use for when searching with lastz");
	setUp.setOption(extractPars.lzPars.coverage, "--coverage", "Coverage to use for when searching with lastz");
	setUp.setOption(genomeMappingPars.numThreads_, "--numThreads", "Number of CPUs to utilize");
  setUp.setOption(gffDir, "--gffDir", "A directory with a gff for the genomes in --genomeDir, should be named GENOME.gff (for GENOME.fasta)");
  setUp.setOption(gffExtraAttributesStr, "--gffExtraAttributes", "Extra attributes to add to genome that has an accompany gff");
	setUp.setOption(extractPars.keepBestOnly, "--keepBestOnly", "Keep best hits only");
	setUp.setOption(keepRefAlignments, "--keepRefAlignments", "Keep Ref Alignments");

	setUp.setOption(bedFile, "--bed", "Bed file, first entry is used", true);
  setUp.setOption(genomeMappingPars.genomeDir_, "--genomeDir", "Name of the genome file fnp", true);
	setUp.setOption(genomeMappingPars.primaryGenome_, "--primaryGenome", "Name of the primary genome, should be in --genomeDir", true);
	//genomeMappingPars.setGenomeFnp();
	std::string selectedGenomesStr = "";
	setUp.setOption(selectedGenomesStr, "--selectedGenomes", "Name of the other genomes in --genomeDir to be read in, leave blank to just do all fastas");

	std::string extraGffAttributes = "";
	setUp.setOption(extraGffAttributes, "--extraGffAttributes", "Extra Gff Attributes");
	if("" != extraGffAttributes){
		genomeMappingPars.gffIntersectPars_.extraAttributes_ = tokenizeString(extraGffAttributes, ",");
	}
	//genomeMappingPars.numThreads_
	setUp.setOption(overWriteDirs, "--overWriteDirs","Over write region dirs");
	setUp.setOption(outputDir, "--outputDir","Output Directory to write files, will be created if it doesn't already exist", true);
	setUp.finishSetUp(std::cout);

	njh::files::makeDirP(njh::files::MkdirPar(outputDir));
	setUp.startARunLog(outputDir.string());
	auto regions = gatherRegions(bedFile.string(), "", setUp.pars_.verbose_);
	uint32_t maxlen = 0;
	for(const auto & region : regions){
		if(region.getLen() > maxlen){
			maxlen = region.getLen();
		}
	}
	njh::concurrent::LockableQueue<GenomicRegion> regionsQueue(regions);

	std::unique_ptr<MultiGenomeMapper> gMapper;
	//set up genome mapper;
	gMapper = std::make_unique<MultiGenomeMapper>(genomeMappingPars);
	//set threads;
	uint32_t originalNumThreads = genomeMappingPars.numThreads_;
	if(genomeMappingPars.numThreads_ >= 4 && regions.size() >=4){
		gMapper->pars_.numThreads_ = 2;
		genomeMappingPars.numThreads_ = genomeMappingPars.numThreads_/2;
	}
	//set gff dir
	gMapper->pars_.gffDir_ = gffDir;

	//set up selected genomes
	gMapper->setSelectedGenomes(selectedGenomesStr);
	gMapper->init();


	auto gapPars = gapScoringParameters(5, 1, 0, 0, 0, 0);
	auto scoring = substituteMatrix(2, -2);

	uint64_t maxAlignLen = maxlen + extractPars.extendAndTrimLen * 2 + 20;

	std::vector<aligner> aligners;
	for(uint32_t t = 0; t < genomeMappingPars.numThreads_; ++t){
		aligners.emplace_back(aligner(maxAlignLen, gapPars, scoring));
	}
//	std::cout << aligners[0].parts_.gapScores_.getIdentifer() << std::endl;
//	aligners[0].parts_.scoring_.printScores(std::cout);
//  auto refSeqs = gMapper->getRefSeqsWithPrimaryGenome(region, refAlignsDir, lzPars, keepBestOnly);
	auto extractPathway =
			[&regionsQueue,&keepRefAlignments,&outputDir,&extractPars,&overWriteDirs,&aligners,&writeOutAllSeqsFile
			 ](
				const std::unique_ptr<MultiGenomeMapper> & gMapper, uint32_t threadNumber) {
				GenomicRegion region;
				while(regionsQueue.getVal(region)) {
					std::string rUid = region.createUidFromCoordsStrand();
					region.setUidWtihCoordsStrand();
					auto regionDir = njh::files::make_path(outputDir,rUid);

					bool needToRun = true;
					if(bfs::exists(regionDir) && !overWriteDirs){
						needToRun = false;
					}
					if(needToRun){
						auto regDirPars = njh::files::MkdirPar(regionDir);
						regDirPars.overWriteDir_ = true;
						njh::files::makeDir(regDirPars);
						auto refAlignsDir = njh::files::makeDir(regionDir,
								njh::files::MkdirPar("refAlignments"));
						auto bedsDir = njh::files::makeDir(regionDir,njh::files::MkdirPar("beds"));
						auto refSeqs = gMapper->getRefSeqsWithPrimaryGenome(region,
						refAlignsDir, extractPars, aligners[threadNumber]);
						if(writeOutAllSeqsFile){
							auto separatedSeqFileOpts = SeqIOOptions::genFastaOut(njh::files::make_path(regionDir, "allRefsSeparated.fasta"));
							SeqOutput separatedWriter(separatedSeqFileOpts);
							separatedWriter.openOut();
							for(const auto & refSeq : refSeqs){
								auto toks = tokenizeString(refSeq.name_, "-");//unfortuantley there's not a better way to separate so hopefully the ref genomes don't have - in their names
								for(const auto & tok : toks){
									separatedWriter.write(seqInfo{tok, refSeq.seq_});
								}
							}
						}
						SeqOutput::write(refSeqs,
								SeqIOOptions::genFastaOut(
										njh::files::make_path(regionDir, "allRefs.fasta")));

						bfs::path extractionCountsFnp = njh::files::make_path(refAlignsDir, "extractionCounts.tab.txt");
						bfs::copy_file(extractionCountsFnp,
																njh::files::make_path(bedsDir , "extractionCounts.tab.txt" ));
						for(const auto & genome : getVectorOfMapKeys(gMapper->genomes_)) {
							auto gRegionPath = njh::files::make_path(refAlignsDir, genome + "_regions.bed" );
							if(bfs::exists(gRegionPath)) {
								bfs::rename(gRegionPath,
										njh::files::make_path(bedsDir , genome + "_region.bed" ));
							}
						}
						if(!keepRefAlignments){
							njh::files::rmDirForce(refAlignsDir);
						}
					}
				}
			};
	std::vector<std::thread> threads;
	for(uint32_t t = 0; t < genomeMappingPars.numThreads_; ++t){
		threads.emplace_back(std::thread(extractPathway,
				std::cref(gMapper), t));
	}
	njh::concurrent::joinAllThreads(threads);


	if(combineByGenome){
		njh::files::MkdirPar combinedByGenomeMkPars(njh::files::make_path(outputDir, "combinedByGenome"));
		combinedByGenomeMkPars.overWriteDir_ = overWriteDirs;
		njh::files::makeDir(combinedByGenomeMkPars);
		auto genomeNames = getVectorOfMapKeys(gMapper->genomes_);
		njh::sort(genomeNames);
		njh::concurrent::LockableQueue<std::string> genomeQueue(genomeNames);

		std::function<void()> combineForGenome = [&genomeQueue,&outputDir,&regions,&combinedByGenomeMkPars](){
			std::string genome = "";
			while(genomeQueue.getVal(genome)){
				std::vector<GenomicRegion> rawRegionsForGenome;
				for(const auto & region : regions){
					auto bedFnp = njh::files::make_path(outputDir, region.createUidFromCoordsStrand(), "beds", genome + "_region.bed");
					if(bfs::exists(bedFnp)){
						addOtherVec(rawRegionsForGenome, bedPtrsToGenomicRegs(getBeds(bedFnp)));
					}
				}
				std::vector<GenomicRegion> regionsForGenome = mergeRegionsStrandAware(rawRegionsForGenome);
				std::vector<Bed6RecordCore> bedsOuts;
				for(const auto & reg : regionsForGenome){
					bedsOuts.emplace_back(reg.genBedRecordCore());
				}
//				if(bfs::exists(gMapper->genomes_.at(genome)->gffFnp_)){
//					intersectBedLocsWtihGffRecordsPars interPars = gMapper->pars_.gffIntersectPars_;
//					interPars.gffFnp_ = gMapper->genomes_.at(genome)->gffFnp_;
//					intersectBedLocsWtihGffRecords(bedsOuts, interPars);
//				}
				OutputStream bedOut(njh::files::make_path(combinedByGenomeMkPars.dirName_, genome + ".bed"));
				for(const auto & bed : bedsOuts){
					bedOut << bed.toDelimStrWithExtra() << std::endl;
				}
			}
		};
		njh::concurrent::runVoidFunctionThreaded(combineForGenome, originalNumThreads);
	}


	return 0;
}





int genExpRunner::extractExpectedRefSeqsFromGenomes(
		const njh::progutils::CmdArgs & inputCommands) {

	MultiGenomeMapper::inputParameters genomeMappingPars;
	bfs::path bedFile = "";
	bfs::path outputDir = "";
	bool keepRefAlignments = false;
	bool overWriteDirs = false;
	MultiGenomeMapper::getRefSeqsWithPrimaryGenomePars extractPars;
	extractPars.extendAndTrimLen = 200;
	extractPars.lzPars.identity = 80;
	extractPars.lzPars.coverage = 80;

	bfs::path gffDir = "";
	std::string gffExtraAttributesStr = "description";

	VecStr descriptions{};
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(extractPars.extendAndTrim, "--extendAndTrim", "Extend the determine region and then trim back, can be helpful for when variation falls at the very ends of the sequence");
	setUp.setOption(extractPars.extendAndTrimLen, "--extendAndTrimLen", "When extending and trimming, use this length");
	setUp.setOption(extractPars.lzPars.identity, "--identity", "Identity to use for when searching with lastz");
	setUp.setOption(extractPars.lzPars.coverage, "--coverage", "Coverage to use for when searching with lastz");
	setUp.setOption(genomeMappingPars.numThreads_, "--numThreads", "Number of CPUs to utilize");
  setUp.setOption(gffDir, "--gffDir", "A directory with a gff for the genomes in --genomeDir, should be named GENOME.gff (for GENOME.fasta)", true);
  setUp.setOption(descriptions, "--descriptions", "Descriptions of genes to add to expected sequences", true);
  setUp.setOption(gffExtraAttributesStr, "--gffExtraAttributes", "Extra attributes to add to genome that has an accompany gff");
	setUp.setOption(extractPars.keepBestOnly, "--keepBestOnly", "Keep best hits only");
	setUp.setOption(keepRefAlignments, "--keepRefAlignments", "Keep Ref Alignments");

	setUp.setOption(bedFile, "--bed", "Bed file, first entry is used", true);
  setUp.setOption(genomeMappingPars.genomeDir_, "--genomeDir", "Name of the genome file fnp", true);
	setUp.setOption(genomeMappingPars.primaryGenome_, "--primaryGenome", "Name of the primary genome, should be in --genomeDir", true);
	//genomeMappingPars.setGenomeFnp();
	std::string selectedGenomesStr = "";
	setUp.setOption(selectedGenomesStr, "--selectedGenomes", "Name of the other genomes in --genomeDir to be read in, leave blank to just do all fastas");
	setUp.setOption(genomeMappingPars.gffIntersectPars_.extraAttributes_, "--extraGffAttributes", "Extra Gff Attributes");
	setUp.setOption(genomeMappingPars.gffIntersectPars_.selectFeatures_, "--features", "Features to extract that match given descriptions to add to expected");
	setUp.setOption(overWriteDirs, "--overWriteDir","Over write region dirs");
	setUp.setOption(outputDir, "--outputDir","Output Directory to write files, will be created if it doesn't already exist", true);
	setUp.finishSetUp(std::cout);


	njh::files::makeDir(njh::files::MkdirPar(outputDir,overWriteDirs));
	setUp.startARunLog(outputDir.string());
	auto regions = gatherRegions(bedFile.string(), "", setUp.pars_.verbose_);
	uint32_t maxlen = 0;
	for(const auto & region : regions){
		if(region.getLen() > maxlen){
			maxlen = region.getLen();
		}
	}
	njh::concurrent::LockableQueue<GenomicRegion> regionsQueue(regions);

	std::unique_ptr<MultiGenomeMapper> gMapper;
	//set up genome mapper;
	gMapper = std::make_unique<MultiGenomeMapper>(genomeMappingPars);
	//set threads;
	uint32_t originalNumThreads = genomeMappingPars.numThreads_;
	if(genomeMappingPars.numThreads_ >= 4 && regions.size() >=4){
		gMapper->pars_.numThreads_ = 2;
		genomeMappingPars.numThreads_ = genomeMappingPars.numThreads_/2;
	}
	//set gff dir
	gMapper->pars_.gffDir_ = gffDir;

	//set up selected genomes
	gMapper->setSelectedGenomes(selectedGenomesStr);
	gMapper->init();


	auto gapPars = gapScoringParameters(5, 1, 0, 0, 0, 0);
	auto scoring = substituteMatrix(2, -2);

	uint64_t maxAlignLen = maxlen + extractPars.extendAndTrimLen * 2 + 20;

	std::vector<aligner> aligners;
	for(uint32_t t = 0; t < genomeMappingPars.numThreads_; ++t){
		aligners.emplace_back(aligner(maxAlignLen, gapPars, scoring));
	}
//	std::cout << aligners[0].parts_.gapScores_.getIdentifer() << std::endl;
//	aligners[0].parts_.scoring_.printScores(std::cout);
//  auto refSeqs = gMapper->getRefSeqsWithPrimaryGenome(region, refAlignsDir, lzPars, keepBestOnly);
	auto extractPathway =
			[&regionsQueue,&keepRefAlignments,&outputDir,&extractPars,&overWriteDirs,&aligners](
					const std::unique_ptr<MultiGenomeMapper> & gMapper, uint32_t threadNumber) {
				GenomicRegion region;
				while(regionsQueue.getVal(region)) {
					std::string rUid = region.createUidFromCoordsStrand();
					region.setUidWtihCoordsStrand();
					auto regionDir = njh::files::make_path(outputDir,rUid);

					bool needToRun = true;
					if(bfs::exists(regionDir) && !overWriteDirs){
						needToRun = false;
					}
					if(needToRun){
						auto regDirPars = njh::files::MkdirPar(regionDir);
						regDirPars.overWriteDir_ = true;
						njh::files::makeDir(regDirPars);
						auto refAlignsDir = njh::files::makeDir(regionDir,
								njh::files::MkdirPar("refAlignments"));
						auto bedsDir = njh::files::makeDir(regionDir,njh::files::MkdirPar("beds"));
						auto refSeqs = gMapper->getRefSeqsWithPrimaryGenome(region, refAlignsDir, extractPars, aligners[threadNumber]);
						SeqOutput::write(refSeqs,
								SeqIOOptions::genFastaOut(
										njh::files::make_path(regionDir, "allRefs.fasta")));
						bfs::path extractionCountsFnp = njh::files::make_path(refAlignsDir, "extractionCounts.tab.txt");
						bfs::copy_file(extractionCountsFnp,
																njh::files::make_path(bedsDir , "extractionCounts.tab.txt" ));
						for(const auto & genome : getVectorOfMapKeys(gMapper->genomes_)) {
							auto gRegionPath = njh::files::make_path(refAlignsDir, genome + "_regions.bed" );
							if(bfs::exists(gRegionPath)) {
								bfs::rename(gRegionPath,
										njh::files::make_path(bedsDir , genome + "_region.bed" ));
							}
						}
						if(!keepRefAlignments){
							njh::files::rmDirForce(refAlignsDir);
						}
					}
				}
			};
	std::vector<std::thread> threads;
	for(uint32_t t = 0; t < genomeMappingPars.numThreads_; ++t){
		threads.emplace_back(std::thread(extractPathway,
				std::cref(gMapper), t));
	}
	njh::concurrent::joinAllThreads(threads);



	njh::files::MkdirPar combinedByGenomeMkPars(njh::files::make_path(outputDir, "combinedByGenome"));
	combinedByGenomeMkPars.overWriteDir_ = overWriteDirs;
	njh::files::makeDir(combinedByGenomeMkPars);
	auto genomeNames = getVectorOfMapKeys(gMapper->genomes_);
	njh::sort(genomeNames);
	njh::concurrent::LockableQueue<std::string> genomeQueue(genomeNames);

	std::function<void()> combineForGenome = [&gMapper, &genomeQueue,&outputDir,&regions,&combinedByGenomeMkPars](){
		std::string genome = "";
		while(genomeQueue.getVal(genome)){
			std::vector<GenomicRegion> rawRegionsForGenome;
			for(const auto & region : regions){
				auto bedFnp = njh::files::make_path(outputDir, region.createUidFromCoordsStrand(), "beds", genome + "_region.bed");
				if(bfs::exists(bedFnp)){
					addOtherVec(rawRegionsForGenome, bedPtrsToGenomicRegs(getBeds(bedFnp)));
				}
			}
			std::vector<GenomicRegion> regionsForGenome = mergeRegionsStrandAware(rawRegionsForGenome);
			std::vector<Bed6RecordCore> bedsOuts;
			for(const auto & reg : regionsForGenome){
				auto toAdd = reg.genBedRecordCore();
				toAdd.extraFields_.clear();
				bedsOuts.emplace_back(toAdd);
			}
			if(bfs::exists(gMapper->genomes_.at(genome)->gffFnp_)){
				intersectBedLocsWtihGffRecordsPars interPars = gMapper->pars_.gffIntersectPars_;
				interPars.gffFnp_ = gMapper->genomes_.at(genome)->gffFnp_;
				intersectBedLocsWtihGffRecords(bedsOuts, interPars);
			}
			OutputStream bedOut(njh::files::make_path(combinedByGenomeMkPars.dirName_, genome + ".bed"));
			for(const auto & bed : bedsOuts){
				bedOut << bed.toDelimStrWithExtra() << std::endl;
			}
		}
	};
	njh::concurrent::runVoidFunctionThreaded(combineForGenome, originalNumThreads);


	njh::files::MkdirPar expectedByGenomeMkPars(njh::files::make_path(outputDir, "expectedRegionsByGenome"));
	expectedByGenomeMkPars.overWriteDir_ = overWriteDirs;
	njh::files::makeDir(expectedByGenomeMkPars);



	njh::concurrent::LockableQueue<std::string> againGenomeQueue(genomeNames);

	std::function<void()> combineForExpectedRegions = [&gMapper,&againGenomeQueue,
																										 &expectedByGenomeMkPars,&combinedByGenomeMkPars,
																										 &descriptions](){
		std::string genome = "";
		while(againGenomeQueue.getVal(genome)){


			std::vector<GenomicRegion> expRegionsByDescription;
			{
				BioDataFileIO<GFFCore> reader((IoOptions(InOptions(gMapper->genomes_.at(genome)->gffFnp_))));
				reader.openIn();
				uint32_t count = 0;
				std::string line = "";
				std::shared_ptr<GFFCore> gRecord = reader.readNextRecord();
				while (nullptr != gRecord) {
					if(gMapper->pars_.gffIntersectPars_.selectFeatures_.empty() || njh::in(gRecord->type_, gMapper->pars_.gffIntersectPars_.selectFeatures_)){
						if(njh::in(gRecord->getAttr("description"), descriptions)){
							auto bedOut = GenomicRegion(*gRecord).genBedRecordCore();
							bedOut.extraFields_.emplace_back(njh::pasteAsStr("[",
									"id=",gRecord->getAttr("ID"), ";",
									"description=",gRecord->getAttr("description"), ";",
									"]"));
							expRegionsByDescription.emplace_back(bedOut);
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
			}

			auto extractBedFnp = njh::files::make_path(combinedByGenomeMkPars.dirName_, genome + ".bed");
			std::vector<GenomicRegion> extractedExpRegionsNoID;
			std::vector<GenomicRegion> extractedExpRegionsByID;

			if (bfs::exists(extractBedFnp)) {
				auto extractedRegions = bedPtrsToGenomicRegs(getBeds(extractBedFnp));
				std::set<std::string> idsNotPreviouslyFound;
				for (const auto & extractedReg : extractedRegions) {
					bool noOverlap = true;
					for (const auto & expRegByDes : expRegionsByDescription) {
						if (expRegByDes.overlaps(extractedReg)) {
							noOverlap = false;
							break;
						}
					}
					if (noOverlap) {
						if (extractedReg.meta_.containsMeta("ID")) {
							idsNotPreviouslyFound.emplace(extractedReg.meta_.getMeta("ID"));
						} else {
							extractedExpRegionsNoID.emplace_back(extractedReg);
						}
					}
				}
				if(!idsNotPreviouslyFound.empty()){
					BioDataFileIO<GFFCore> reader((IoOptions(InOptions(gMapper->genomes_.at(genome)->gffFnp_))));
					reader.openIn();
					uint32_t count = 0;
					std::string line = "";
					std::shared_ptr<GFFCore> gRecord = reader.readNextRecord();
					while (nullptr != gRecord) {
						if(gMapper->pars_.gffIntersectPars_.selectFeatures_.empty() || njh::in(gRecord->type_, gMapper->pars_.gffIntersectPars_.selectFeatures_)){
							if(njh::in(gRecord->getAttr("ID"), idsNotPreviouslyFound)){
								auto bedOut = GenomicRegion(*gRecord).genBedRecordCore();
								bedOut.extraFields_.emplace_back(njh::pasteAsStr("[",
										"id=",gRecord->getAttr("ID"), ";",
										"description=",gRecord->getAttr("description"), ";",
										"]"));
								extractedExpRegionsByID.emplace_back(bedOut);
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
				}
			}
			OutputStream bedOut(njh::files::make_path(expectedByGenomeMkPars.dirName_, genome + ".bed"));
			std::vector<GenomicRegion> allRegions;
			addOtherVec(allRegions, extractedExpRegionsNoID);
			addOtherVec(allRegions, extractedExpRegionsByID);
			addOtherVec(allRegions, expRegionsByDescription);
			sortGRegionsByStart(allRegions);
			allRegions = mergeRegionsStrandAware(allRegions);
			for(const auto & outReg : allRegions){
				bedOut << outReg.genBedRecordCore().toDelimStrWithExtra() << std::endl;
			}
		}
	};

	njh::concurrent::runVoidFunctionThreaded(combineForExpectedRegions, originalNumThreads);



	return 0;
}


} // namespace njhseq


