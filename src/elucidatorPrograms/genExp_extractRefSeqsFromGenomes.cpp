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
	std::string forwardPrimer = "";
	std::string reversePrimer = "";
	std::string targetName = "";
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	pars.setUpCoreOptions(setUp, false);
	bool setPrimersFile = setUp.setOption(pars.primersFile, "--primers", "A file that contains three columns, targetName,forwardPrimer,reversePrimer 5` to 3` directions");
	if(!setPrimersFile){
		setUp.processSeq(forwardPrimer, "--forwardPrimer", "Forward Primer (5`to 3` direction)", true);
		setUp.processSeq(reversePrimer, "--reversePrimer", "Reverse Primer (5`to 3` direction)", true);
		setUp.processSeq(targetName, "--targetName", "Name of the target associated with forward and reverse primers", true);

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
		primerTable.addRow(targetName, forwardPrimer, reversePrimer);
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

int genExpRunner::extractRefSeqsFromGenomes(
		const njh::progutils::CmdArgs & inputCommands) {
	BioCmdsUtils::LastZPars lzPars;
	lzPars.identity = 80;
	lzPars.coverage = 90;
	MultiGenomeMapper::inputParameters genomeMappingPars;
	bfs::path bedFile = "";
	bfs::path outputDir = "";
	bool keepRefAlignments = false;
	bool overWriteDirs = false;
	bool keepBestOnly = false;

	bool extendAndTrim = false;
	uint32_t extendAndTrimLen = 10;
	bfs::path gffDir = "";
	std::string gffExtraAttributesStr = "description";
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(extendAndTrim, "--extendAndTrim", "Extend the determine region and then trim back, can be helpful for when variation falls at the very ends of the sequence");
	setUp.setOption(extendAndTrimLen, "--extendAndTrimLen", "When extending and trimming, use this length");
	setUp.setOption(lzPars.identity, "--identity", "Identity to use for when searching with lastz");
	setUp.setOption(lzPars.coverage, "--coverage", "Coverage to use for when searching with lastz");
	setUp.setOption(genomeMappingPars.numThreads_, "--numThreads", "Number of CPUs to utilize");
  setUp.setOption(gffDir, "--gffDir", "A directory with a gff for the genomes in --genomeDir, should be named GENOME.gff (for GENOME.fasta)");
  setUp.setOption(gffExtraAttributesStr, "--gffExtraAttributes", "Extra attributes to add to genome that has an accompany gff");
	setUp.setOption(keepBestOnly, "--keepBestOnly", "Keep best hits only");
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
	if(genomeMappingPars.numThreads_ >= 4){
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

	uint64_t maxAlignLen = maxlen + extendAndTrimLen * 2 + 20;

	std::vector<aligner> aligners;
	for(uint32_t t = 0; t < genomeMappingPars.numThreads_; ++t){
		aligners.emplace_back(aligner(maxAlignLen, gapPars, scoring));
	}
//	std::cout << aligners[0].parts_.gapScores_.getIdentifer() << std::endl;
//	aligners[0].parts_.scoring_.printScores(std::cout);
//  auto refSeqs = gMapper->getRefSeqsWithPrimaryGenome(region, refAlignsDir, lzPars, keepBestOnly);
	auto extractPathway =
			[&regionsQueue,&keepRefAlignments,&outputDir,&lzPars,&overWriteDirs,&keepBestOnly,&aligners,
			 &extendAndTrim,&extendAndTrimLen](
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
						auto refSeqs = gMapper->getRefSeqsWithPrimaryGenome(region, refAlignsDir, lzPars, keepBestOnly, extendAndTrim, extendAndTrimLen, aligners[threadNumber]);
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

	return 0;
}

} // namespace njhseq


