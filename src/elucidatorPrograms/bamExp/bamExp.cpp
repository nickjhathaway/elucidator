/*
 * bamExpRunner.cpp
 *
 *  Created on: Jan 25, 2015
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


#include "bamExp.hpp"
#include <TwoBit.h>

#include "elucidator/BamToolsUtils.h"
#include <njhseq/GenomeUtils.h>
#include <njhseq/objects/BioDataObject.h>
#include <njhseq/programUtils/seqSetUp.hpp>




namespace njhseq {
bamExpRunner::bamExpRunner()
    : njh::progutils::ProgramRunner(
          {
					 addFunc("bamReadProfiling", bamReadProfiling, false),
					 addFunc("bamToFastq", bamToFastq, false),
					 addFunc("getInsertSizesStats", getInsertSizesStats, false),
					 addFunc("BamExtractReadsFromRegeion", BamExtractReadsFromRegion, false),
					 addFunc("printBamRefIds", printBamRefIds, false),
					 addFunc("bamToBed", bamToBed, false),
					 addFunc("getUnmappedAlnsFromBam", getUnmappedAlnsFromBam, false),
					 addFunc("getMateMapStatus", getMateMapStatus, false),
					 addFunc("getInsertSizeChanges", getInsertSizeChanges, false),
					 addFunc("getMapQualityCounts", getMapQualityCounts, false),
					 addFunc("bamMulticov", bamMulticov, false),
					 addFunc("bamMulticovBases", bamMulticovBases, false),
					 addFunc("determineRegion", determineRegion, false),
					 addFunc("determineRegionLastz", determineRegionLastz, false),
					 addFunc("getBestBedRegionFromBam", getBestBedRegionFromBam, false),
					 addFunc("mergeBedRegionFromBam", mergeBedRegionFromBam, false),
					 addFunc("getUnmappedAlnsFromBam", getUnmappedAlnsFromBam, false),
					 addFunc("countUnMappedMateStatus", countUnMappedMateStatus, false),
					 addFunc("BamGetFileIndexPositionOfName", BamGetFileIndexPositionOfName, false),
					 addFunc("BamFindDifferenceInUnmappedFileIndexPosition", BamFindDifferenceInUnmappedFileIndexPosition, false),
					 addFunc("isBamSorted", isBamSorted, false),
					 addFunc("multiBamCoverageFinder", multiBamCoverageFinder, false),
					 addFunc("genBedFromMappedGenome", genBedFromMappedGenome, false),
					 addFunc("appendReadGroupToName", appendReadGroupToName, false),
					 addFunc("GetBamLocationStats", GetBamLocationStats, false),
					 addFunc("outputSoftClipCounts", outputSoftClipCounts, false),
					 addFunc("bamMultiPairStats", bamMultiPairStats, false),
					 addFunc("GetKmerCoverageForRegions", GetKmerCoverageForRegions, false),
					 addFunc("BamRenameRefHeader", BamRenameRefHeader, false),
					 addFunc("testingBamToFastq", testingBamToFastq, false),
					 addFunc("printHeaderRefIndexes", printHeaderRefIndexes, false),
					 addFunc("BamFilterByChroms", BamFilterByChroms, false),
					 addFunc("refineBedRegionFromBam", refineBedRegionFromBam, false),
					 addFunc("bamMulticovBasesRough", bamMulticovBasesRough, false),
					 addFunc("bamDupCounts", bamDupCounts, false),
					 addFunc("BamRefIdsToBed", BamRefIdsToBed, false),
					 addFunc("bamToFastqAlns", bamToFastqAlns, false),
					 addFunc("BamExtractReAlignToRef", BamExtractReAlignToRef, false),
					 addFunc("BamCreateErrorProfileToRef", BamCreateErrorProfileToRef, false),
					 addFunc("BamGetSpanningReadsForRegion", BamGetSpanningReadsForRegion, false),
					 addFunc("BamGetPileupForRegion", BamGetPileupForRegion, false),
					 addFunc("MultipleBamGetPileupForRegion", MultipleBamGetPileupForRegion, false),
					 addFunc("BamFilterAlnsOnLenAndQual", BamFilterAlnsOnLenAndQual, false),
					 addFunc("testingBamToFastqToContigs", testingBamToFastqToContigs, false),
					 addFunc("getUnmappedAndShortAlnsFromBam", getUnmappedAndShortAlnsFromBam, false),
					 addFunc("multiBamCoverageFinderBases", multiBamCoverageFinderBases, false),
					 addFunc("BamGetImproperPairsOnChroms", BamGetImproperPairsOnChroms, false),
					 addFunc("BamGetPairedReadInfoForRegions", BamGetPairedReadInfoForRegions, false),
					 addFunc("BamGetImproperPairCounts", BamGetImproperPairCounts, false),
					 addFunc("bamCovBasesRough", bamCovBasesRough, false),
					 addFunc("bamCov", bamCov, false),
					 addFunc("BamFilterByChromsToBam", BamFilterByChromsToBam, false),
					 addFunc("BamGetSpanningReadsForRegionLongReads", BamGetSpanningReadsForRegionLongReads, false),
           addFunc("getMateMapLocationForRegion", getMateMapLocationForRegion, false),
          	 addFunc("mergeImproperMatePositions", mergeImproperMatePositions, false),
          	addFunc("connectFaceAwayMatesRegions", connectFaceAwayMatesRegions, false),
          },//

				"bamExp") {
}

int bamExpRunner::refineBedRegionFromBam(const njh::progutils::CmdArgs & inputCommands) {
	RegionRefinementPars refinePars;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(refinePars.numThreads, "--numThreads", "Number of threads to use");
	setUp.setOption(refinePars.bedFnp, "--bedFnp", "Bed Regions", true);
	setUp.setOption(refinePars.bamFnp, "--bam", "Bam file", true);
	setUp.processWritingOptions(refinePars.outOpts);
	setUp.finishSetUp(std::cout);

	RunRegionRefinement(refinePars);

	return 0;
}

int bamExpRunner::printHeaderRefIndexes(
		const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bamFnp = "";
	OutOptions outOpts(bfs::path(""), ".txt");

	bool header = false;

	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(bamFnp, "--bamFnp", "Bam file", true);
	setUp.setOption(header, "--header", "Print a column header");
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	njh::files::checkExistenceThrow(bamFnp, __PRETTY_FUNCTION__);
	BamTools::BamReader bReader;
	bReader.Open(bamFnp.string());
	OutputStream out(outOpts);
	auto refData = bReader.GetReferenceData();
	if (header) {
		out << "Index\tRefName\tRefLength" << std::endl;
	}
	for (const auto pos : iter::range(refData.size())) {
		out << pos << "\t" << refData[pos].RefName << "\t" << refData[pos].RefLength
				<< std::endl;
	}
	return 0;
}

int bamExpRunner::getUnmappedAlnsFromBam(const njh::progutils::CmdArgs & inputCommands){
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processDefaultReader({"--bam"}, true);
	setUp.finishSetUp(std::cout);

	BamExtractor bExtractor(setUp.pars_.verbose_);
	auto outFiles = bExtractor.writeUnMappedSeqs(setUp.pars_.ioOptions_);

	return 0;
}


int bamExpRunner::genBedFromMappedGenome(
		const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bamFnp = "";
	OutOptions outOpts;
	outOpts.outExtention_ = ".bed";
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(bamFnp, "--bam", "Bam file to check", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	BamTools::BamReader bReader;
	bReader.Open(bamFnp.string());
	checkBamOpenThrow(bReader, bamFnp);
	auto refData = bReader.GetReferenceData();
	auto regions = genGenRegionsFromRefData(refData);

	OutputStream out(outOpts);
	for(const auto & reg : regions){
		out << reg.genBedRecordCore().toDelimStrWithExtra() << std::endl;
	}
	return 0;
}


int bamExpRunner::isBamSorted(
		const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bamFnp = "";
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(bamFnp, "--bam", "Bam file to check", true);
	setUp.finishSetUp(std::cout);

	auto res = IsBamSorted(bamFnp.string(), setUp.pars_.verbose_);
	std::cout << njh::boolToStr(res) << std::endl;

	return 0;
}



int bamExpRunner::determineRegion(
		const njh::progutils::CmdArgs & inputCommands) {
	bfs::path genomeFnp = "";
	std::string name;
	bool keepIntermediateFiles = false;
	bool individual = false;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();

	setUp.pars_.ioOptions_.out_.outFilename_ = "";
	setUp.processDefaultReader(VecStr { "--fasta", "--fastq"});
	setUp.pars_.ioOptions_.out_.outExtention_ = ".bed";
	setUp.setOption(genomeFnp, "--genome",
			"Path to a genome to determine the location in", true);
	setUp.setOption(name, "--name",
				"Name to give the determined region");
	setUp.setOption(individual, "--individual",
					"Report a region for each input sequence rather than just the best overall regions for all seqs");
	setUp.setOption(keepIntermediateFiles, "--keepIntermediateFiles",
									"keep Intermediate Files");
	setUp.finishSetUp(std::cout);
	auto genomeName = bfs::basename(genomeFnp);
	if (std::string::npos != genomeName.rfind('.')) {
		genomeName = genomeName.substr(0, genomeName.rfind('.'));
	}

	MultiGenomeMapper genomeMapper(genomeFnp.parent_path(),
			genomeName);
	genomeMapper.setSelectedGenomes(VecStr { genomeName });
	genomeMapper.loadInGenomes();
	//open out file
	OutputStream out(setUp.pars_.ioOptions_.out_);
	//map reads
	auto outputs = genomeMapper.alignToGenomes(setUp.pars_.ioOptions_, "./");
	std::unordered_map<std::string, bfs::path> bamFnps;
	for (const auto & output : outputs) {
		bamFnps[output.first] = output.second.alignedFnp_;
	}

	//get regions
	if(individual){
		auto bamFnp = bamFnps.begin()->second;
		BamTools::BamReader bReader;
		BamTools::BamAlignment bAln;
		bReader.Open(bamFnp.string());
		checkBamOpenThrow(bReader, bamFnp.string());
		auto refDAta = bReader.GetReferenceData();
		while(bReader.GetNextAlignment(bAln)){
			GenomicRegion region{bAln, refDAta};
			out << region.genBedRecordCore().toDelimStr() << std::endl;
		}
	}else{
		auto determinedRegions = genomeMapper.getRegionsFromBams(bamFnps);
		if(!determinedRegions.at(genomeName).empty()){
			auto region = determinedRegions.at(genomeName).front();
			if(!name.empty()){
				region.uid_ = name;
			}
			//write out all the regions determined
			out << region.genBedRecordCore().toDelimStr() << std::endl;
		}
	}
	if(!keepIntermediateFiles){
		for(const auto & bamFnp : bamFnps){
			bfs::remove(bamFnp.second);
			bfs::remove(bamFnp.second.string() + ".bai");
		}
	}
	return 0;
}

int bamExpRunner::determineRegionLastz(
		const njh::progutils::CmdArgs & inputCommands) {
	bfs::path genomeFnp = "";
	std::string name;
	bool individual = false;
	BioCmdsUtils::LastZPars lzPars;
	bool keepIntermediateFiles = false;


	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(lzPars.coverage, "--coverage", "Coverage for lastz");
	setUp.setOption(lzPars.identity, "--identity", "Identity for lastz");
	setUp.setOption(lzPars.extraLastzArgs, "--extraLastzArgs", "Extra Lastz Arguments");
	setUp.pars_.ioOptions_.out_.outFilename_ = "";
	setUp.processDefaultReader(VecStr { "--fasta", "--fastq"});
	setUp.pars_.ioOptions_.out_.outExtention_ = ".bed";
	setUp.setOption(genomeFnp, "--genome",
			"Path to a genome to determine the location in", true);
	setUp.setOption(name, "--name",
				"Name to give the determined region");
	setUp.setOption(individual, "--individual",
					"Report a region for each input sequence rather than just the best overall regions for all seqs");
	setUp.setOption(keepIntermediateFiles, "--keepIntermediateFiles",
									"keep Intermediate Files");
	setUp.finishSetUp(std::cout);

	njh::files::checkExistenceThrow(genomeFnp, __PRETTY_FUNCTION__);

	auto genomeName = bfs::basename(genomeFnp);
	if (std::string::npos != genomeName.rfind(".")) {
		genomeName = genomeName.substr(0, genomeName.rfind("."));
	}
	MultiGenomeMapper genomeMapper(genomeFnp.parent_path(),
			genomeName);
	genomeMapper.setSelectedGenomes(VecStr { genomeName });
	genomeMapper.loadInGenomes();
	if(genomeMapper.genomes_.empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", " << genomeFnp << " not found" << "\n";
		throw std::runtime_error{ss.str()};
	}

	//work around for read names that are too long
	auto tempSeqOpts = SeqIOOptions::genFastaOut(njh::files::findNonexitantFile(njh::files::prependFileBasename(setUp.pars_.ioOptions_.firstName_, "temp_")));
	std::unordered_map<std::string, std::string> nameKey;
	{
		SeqOutput tempWriter(tempSeqOpts);
		tempWriter.openOut();
		seqInfo seq;
		SeqInput reader(setUp.pars_.ioOptions_);
		reader.openIn();
		uint32_t count = 0;
		while(reader.readNextRead(seq)){
			nameKey[estd::to_string(count)] = seq.name_;
			seq.name_ = estd::to_string(count);
			tempWriter.write(seq);
			++count;

		}
	}

	//open out file
	std::ofstream outFile;
	OutputStream out(setUp.pars_.ioOptions_.out_);


	//map reads
	auto tempInOpts = SeqIOOptions::genFastaIn(tempSeqOpts.out_.outName());
	tempInOpts.out_.transferOverwriteOpts(setUp.pars_.ioOptions_.out_);
	auto outputs = genomeMapper.alignToGenomesLastz(tempInOpts, bfs::basename(setUp.pars_.ioOptions_.out_.outFilename_) + "_", lzPars);
	std::unordered_map<std::string, bfs::path> bamFnps;
	for (const auto & output : outputs) {
		bamFnps[output.first] = output.second.alignedFnp_;
	}

	//get regions
	if(individual){

		VecStr added;
		auto bamFnp = bamFnps.begin()->second;
		BamTools::BamReader bReader;
		BamTools::BamAlignment bAln;
		bReader.Open(bamFnp.string());
		checkBamOpenThrow(bReader, bamFnp.string());
		auto refDAta = bReader.GetReferenceData();
		while(bReader.GetNextAlignment(bAln)){
			GenomicRegion region{bAln, refDAta};
			region.uid_ = nameKey[region.uid_];
			out << region.genBedRecordCore().toDelimStr() << std::endl;
			added.emplace_back(region.uid_);
		}
		seqInfo seq;
		SeqInput reader(setUp.pars_.ioOptions_);
		reader.openIn();
		while(reader.readNextRead(seq)){
			if(!njh::in(seq.name_, added)){
				GenomicRegion region{seq.name_, "*", std::numeric_limits<uint32_t>::max()-1, std::numeric_limits<uint32_t>::max(), false};
				out << region.genBedRecordCore().toDelimStr() << std::endl;
			}
		}
	} else {
		auto determinedRegions = genomeMapper.getRegionsFromBams(bamFnps);
		if(!determinedRegions.at(genomeName).empty()){
			auto region = determinedRegions.at(genomeName).front();
			if(!name.empty()){
				region.uid_ = name;
			}
			//write out all the regions determined
			out << region.genBedRecordCore().toDelimStr() << std::endl;
		}
	}



	if(!keepIntermediateFiles){
		for(const auto & bamFnp : bamFnps){
			bfs::remove(bamFnp.second);
			bfs::remove(bamFnp.second.string() + ".bai");
			bfs::remove(tempSeqOpts.out_.outName());
		}
	}
	return 0;
}








int bamExpRunner::outputSoftClipCounts(const njh::progutils::CmdArgs & inputCommands){
	OutOptions outOpts(bfs::path(""));
	//uint32_t numThreads = 1;
	bfs::path bedFnp = "";
	std::string bam;
	outOpts.outExtention_ = ".bed";
	bool noHeader = false;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	//setUp.setOption(numThreads, "--numThreads", "Number of threads to use");
	setUp.processWritingOptions(outOpts);
	setUp.setOption(bedFnp, "--bedFnp", "Bed file of regions to get soft clip information for");
	setUp.setOption(bam, "--bam", "Bam file", true);
	setUp.setOption(noHeader, "--noHeader", "Don't output a header so output can be treated like a real bed file");
	setUp.finishSetUp(std::cout);


	checkBamFilesForIndexesAndAbilityToOpen({bam});
	BamTools::BamReader bReader;
	BamTools::BamAlignment bAln;
	bReader.Open(bam);
	bReader.LocateIndex();
	auto refData = bReader.GetReferenceData();
	OutputStream out(outOpts);
	if(!noHeader){
		out << "#chrom\tstart\tend\tname\tscore\tstrand\tforwardClip\tbackClip" << "\n";
	}
	if(!bedFnp.empty()){
		auto regions = bed3PtrsToGenomicRegs(getBed3s(bedFnp));
		for(const auto & region : regions){
			setBamFileRegionThrow(bReader, region);
			while (bReader.GetNextAlignmentCore(bAln)) {
				if (bAln.IsMapped() && bAln.IsPrimaryAlignment()) {
					if(bAln.CigarData.front().Type == 'S'|| bAln.CigarData.back().Type == 'S'){
						GenomicRegion reg(bAln, refData);
						uint32_t forwardClip = 0;
						uint32_t backClip = 0;
						if(bAln.CigarData.front().Type == 'S'){
							forwardClip = bAln.CigarData.front().Length;
						}
						if(bAln.CigarData.back().Type == 'S'){
							backClip = bAln.CigarData.back().Length;
						}
						auto bedReg = reg.genBedRecordCore();
						bedReg.score_ = bAln.Length;
						bedReg.name_ = region.uid_;
						out << bedReg.toDelimStr()
								<< "\t" << forwardClip
								<< "\t" << backClip << std::endl;;
					}
				}
			}
		}
	} else {
		while (bReader.GetNextAlignmentCore(bAln)) {
			if (bAln.IsMapped()) {
				if(bAln.CigarData.front().Type == 'S'|| bAln.CigarData.back().Type == 'S'){
					GenomicRegion reg(bAln, refData);
					uint32_t forwardClip = 0;
					uint32_t backClip = 0;
					if(bAln.CigarData.front().Type == 'S'){
						forwardClip = bAln.CigarData.front().Length;
					}
					if(bAln.CigarData.back().Type == 'S'){
						backClip = bAln.CigarData.back().Length;
					}
					auto bedReg = reg.genBedRecordCore();
					bedReg.score_ = bAln.Length;
					out << bedReg.toDelimStr()
							<< "\t" << forwardClip
							<< "\t" << backClip << std::endl;
				}
			}
		}
	}

	return 0;
}



int bamExpRunner::bamMulticov(const njh::progutils::CmdArgs & inputCommands){
	OutOptions outOpts(bfs::path(""));
	uint32_t numThreads = 1;
	bfs::path bedFnp = "";
	std::string bams;
	std::string pat = ".*.bam$";
	outOpts.outExtention_ = ".tab.txt";
	bool noHeader = false;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(numThreads, "--numThreads", "Number of threads to use");
	setUp.processWritingOptions(outOpts);
	setUp.setOption(bedFnp, "--bedFnp", "Bed file of regions to get coverage for", true);
	setUp.setOption(pat, "--pat", "Pattern in current directory to get coverage for");
	setUp.setOption(bams, "--bams", "Either a file with the name of a bam file on each line or a comma separated value of bam file paths");
	setUp.setOption(noHeader, "--noHeader", "Don't output a header so output can be treated like a bed file");
	setUp.finishSetUp(std::cout);


	std::ofstream outFile;
	std::ostream out(determineOutBuf(outFile, outOpts));


	auto bamFnps = njh::files::gatherFilesByPatOrNames(std::regex(pat), bams);
	checkBamFilesForIndexesAndAbilityToOpen(bamFnps);
	if(setUp.pars_.verbose_){
		printVector(bamFnps, "\t", std::cout);
	}
	auto bed3s = getBed3s(bedFnp);
	std::vector<GenomicRegion> inputRegions;
	bool allAbove3 = true;
	for(const auto & b : bed3s){
		if(b->extraFields_.size() < 3){
			allAbove3 = false;
			break;
		}
	}
	if (allAbove3) {
		inputRegions = bedPtrsToGenomicRegs(getBeds(bedFnp));
	} else {
		inputRegions = bed3PtrsToGenomicRegs(bed3s);
	}

	//collapse identical regions
	std::vector<GenomicRegion> regions;
	for(const auto & inputRegion : inputRegions){
		if(regions.empty()){
			regions.emplace_back(inputRegion);
		}else{
			if(regions.back().sameRegion(inputRegion)){
				regions.back().uid_ += "," + inputRegion.uid_;
			}else{
				regions.emplace_back(inputRegion);
			}
		}
	}


	struct BamFnpRegionPair {
		BamFnpRegionPair(const bfs::path & bamFnp, const GenomicRegion & region) :
				bamFnp_(bamFnp), region_(region) {
			regUid_ = region_.createUidFromCoords();
			bamFname_ = bamFnp_.filename().string();
		}
		bfs::path bamFnp_;
		GenomicRegion region_;
		uint32_t coverage_ = 0;

		std::string regUid_;
		std::string bamFname_;
	};

	std::vector<std::shared_ptr<BamFnpRegionPair>> pairs;
	for(const auto & bamFnp : bamFnps){
		for(const auto & region : regions){
			pairs.emplace_back(std::make_shared<BamFnpRegionPair>(bamFnp, region));
		}
	}

	njh::concurrent::LockableVec<std::shared_ptr<BamFnpRegionPair>> pairsList(pairs);

	std::function<void()> getCov = std::function<void()>([&pairsList](){
		std::shared_ptr<BamFnpRegionPair> val;
		while(pairsList.getVal(val)){
			BamTools::BamReader bReader;
			bReader.Open(val->bamFnp_.string());
			bReader.LocateIndex();
			setBamFileRegionThrow(bReader, val->region_);
			BamTools::BamAlignment bAln;
			while(bReader.GetNextAlignmentCore(bAln)){
				if(bAln.IsMapped()){
					val->coverage_ += 1;
				}
			}
		}
	});

	{
		njh::concurrent::runVoidFunctionThreaded(getCov, numThreads);
	}

	VecStr header = {"chrom", "start", "end", "name", "score", "strand"};
	std::unordered_map<std::string, uint32_t> bamFnpToCol;
	uint32_t bIndex = 6;
	for(const auto & b : bamFnps){
		header.emplace_back(b.filename().string());
		bamFnpToCol[b.filename().string()] = bIndex;
		++bIndex;
	}

	njh::sort(regions,
			[](const GenomicRegion & reg1,
					const GenomicRegion & reg2) {
				return reg1.createUidFromCoords() < reg2.createUidFromCoords();
			});
	table output(header);
	output.content_ = std::vector<std::vector<std::string>>{regions.size(), std::vector<std::string>{header.size()}};
	uint32_t regRowCount = 0;
	std::unordered_map<std::string, uint32_t> regUidToRowPos;
	for(const auto & reg : regions){
		auto regAsBed = reg.genBedRecordCore();
		auto toks = tokenizeString(regAsBed.toDelimStr(), "\t");
		for(const auto col : iter::range(toks.size())){
			output.content_[regRowCount][col] = toks[col];
		}
		regUidToRowPos[reg.createUidFromCoords()] = regRowCount;
		++regRowCount;
	}
	pairsList.reset();
	std::function<void()> fillTable = [&output, &pairsList, &regUidToRowPos, &bamFnpToCol](){
		std::shared_ptr<BamFnpRegionPair> val;
		while(pairsList.getVal(val)){
			output.content_[regUidToRowPos[val->regUid_]][bamFnpToCol[val->bamFname_]] = estd::to_string(val->coverage_);
		}
	};

	{
		njh::concurrent::runVoidFunctionThreaded(fillTable, numThreads);
	}

	if(noHeader){
		output.hasHeader_ = false;
	}
	output.outPutContents(out, "\t");

	return 0;
}

int bamExpRunner::bamMultiPairStats(const njh::progutils::CmdArgs & inputCommands){
	OutOptions outOpts(bfs::path("out"));
	uint32_t numThreads = 1;
	bfs::path bedFnp = "";
	std::string bams;
	std::string pat = ".*.bam$";
	uint32_t insertSizeCutOff = 5000;
	double percentInRegion = 0.5;
	//outOpts.outExtention_ = ".tab.txt";
	bool noHeader = false;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(insertSizeCutOff, "--insertSizeCutOff", "Insert Size Cut Off to be considered a proper pair");
	setUp.setOption(numThreads, "--numThreads", "Number of threads to use");
	setUp.processWritingOptions(outOpts);
	setUp.setOption(bedFnp, "--bedFnp", "Bed file of regions to get coverage for", true);
	setUp.setOption(pat, "--pat", "Pattern in current directory to get coverage for");
	setUp.setOption(bams, "--bams", "Either a file with the name of a bam file on each line or a comma separated value of bam file paths");
	setUp.setOption(noHeader, "--noHeader", "Don't output a header so output can be treated like a bed file");
	setUp.setOption(percentInRegion, "--percentInRegion", "Percent In Region");
	setUp.finishSetUp(std::cout);


	OutOptions pairCountsOpts(njh::files::nameAppendBeforeExt(outOpts.outName(), "_totalPairCounts"), ".tab.txt");
	OutOptions properPairCountsOpts(njh::files::nameAppendBeforeExt(outOpts.outName(), "_properPairCounts"), ".tab.txt");
	OutOptions mateUnMappedPairCountsOpts(njh::files::nameAppendBeforeExt(outOpts.outName(), "_mateUnmappedPairCounts"), ".tab.txt");
	OutOptions discordantCountsOpts(njh::files::nameAppendBeforeExt(outOpts.outName(), "_discordantPairCounts"), ".tab.txt");

	pairCountsOpts.transferOverwriteOpts(outOpts);
	properPairCountsOpts.transferOverwriteOpts(outOpts);
	mateUnMappedPairCountsOpts.transferOverwriteOpts(outOpts);
	discordantCountsOpts.transferOverwriteOpts(outOpts);

	OutputStream pairCountsOut(pairCountsOpts);
	OutputStream properPairCountsOut(properPairCountsOpts);
	OutputStream mateUnMappedPairCountsOut(mateUnMappedPairCountsOpts);
	OutputStream discordantCountsOut(discordantCountsOpts);


	auto bamFnps = njh::files::gatherFilesByPatOrNames(std::regex{pat}, bams);
	checkBamFilesForIndexesAndAbilityToOpen(bamFnps);
	if(setUp.pars_.verbose_){
		printVector(bamFnps, "\t", std::cout);
	}
	auto bed3s = getBed3s(bedFnp);
	std::vector<GenomicRegion> inputRegions;
	bool allAbove3 = true;
	for(const auto & b : bed3s){
		if(b->extraFields_.size() < 3){
			allAbove3 = false;
			break;
		}
	}
	if (allAbove3) {
		inputRegions = bedPtrsToGenomicRegs(getBeds(bedFnp));
	} else {
		inputRegions = bed3PtrsToGenomicRegs(bed3s);
	}


	njh::sort(inputRegions,
			[](const GenomicRegion & reg1,
					const GenomicRegion & reg2) {
				return reg1.createUidFromCoords() < reg2.createUidFromCoords();
			});

	//collapse identical regions
	std::vector<GenomicRegion> regions;
	for(const auto & inputRegion : inputRegions){
		if(regions.empty()){
			regions.emplace_back(inputRegion);
		}else{
			if(regions.back().sameRegion(inputRegion)){
				regions.back().uid_ += "," + inputRegion.uid_;
			}else{
				regions.emplace_back(inputRegion);
			}
		}
	}


	struct BamFnpRegionPair {
		BamFnpRegionPair(const bfs::path & bamFnp, const GenomicRegion & region) :
				bamFnp_(bamFnp), region_(region) {
			regUid_ = region_.createUidFromCoords();
			bamFname_ = bamFnp_.filename().string();
		}
		bfs::path bamFnp_;
		GenomicRegion region_;
		uint32_t totalPairCnt_ = 0;
		uint32_t properPairCnt_ = 0;
		uint32_t mateUnmappedCnt_ = 0;
		uint32_t discordantCnt_ = 0;

		std::string regUid_;
		std::string bamFname_;
	};

	std::vector<std::shared_ptr<BamFnpRegionPair>> pairs;
	for(const auto & bamFnp : bamFnps){
		for(const auto & region : regions){
			pairs.emplace_back(std::make_shared<BamFnpRegionPair>(bamFnp, region));
		}
	}

	njh::concurrent::LockableVec<std::shared_ptr<BamFnpRegionPair>> pairsList(pairs);

	std::function<void()> getCov = [&pairsList,&insertSizeCutOff,&percentInRegion](){
		std::shared_ptr<BamFnpRegionPair> val;
		while(pairsList.getVal(val)){
			BamTools::BamReader bReader;

			bReader.Open(val->bamFnp_.string());
			bReader.LocateIndex();
			auto refData = bReader.GetReferenceData();
			setBamFileRegionThrow(bReader, val->region_);
			BamTools::BamAlignment bAln;
			while(bReader.GetNextAlignmentCore(bAln)){
				if(!bAln.IsPrimaryAlignment()){
					continue;
				}
				if(!bAln.IsMapped()){
					continue;
				}
				GenomicRegion bAlnReg(bAln, refData);
				double total = std::min(val->region_.getLen(), bAlnReg.getLen());
				if(percentInRegion > 0 && (val->region_.getOverlapLen(bAlnReg)/total) < percentInRegion){
					continue;
				}
				//just counting pairs
				if(bAln.IsPaired()){
					++(val->totalPairCnt_);
					if(bAln.IsMateMapped() && bAln.MateRefID == bAln.RefID && std::abs(bAln.InsertSize) <= insertSizeCutOff){
						++(val->properPairCnt_);
					}else if(!bAln.IsMateMapped()){
						++(val->mateUnmappedCnt_);
					}else if(bAln.IsMateMapped() && (bAln.MateRefID != bAln.RefID || std::abs(bAln.InsertSize) > insertSizeCutOff) ){
						++(val->discordantCnt_);
					}
				}
			}
		}
	};

	{
		njh::concurrent::runVoidFunctionThreaded(getCov, numThreads);
	}

	VecStr header = {"#chrom", "start", "end", "name", "score", "strand"};
	std::unordered_map<std::string, uint32_t> bamFnpToCol;
	uint32_t bIndex = 6;
	for(const auto & b : bamFnps){
		header.emplace_back(b.filename().string());
		bamFnpToCol[b.filename().string()] = bIndex;
		++bIndex;
	}


	table outputProperPairs(header);
	outputProperPairs.content_ = std::vector<std::vector<std::string>>(regions.size(), VecStr(header.size()));
	uint32_t regRowCount = 0;
	std::unordered_map<std::string, uint32_t> regUidToRowPos;
	for(const auto & reg : regions){
		auto regAsBed = reg.genBedRecordCore();
		auto toks = tokenizeString(regAsBed.toDelimStr(), "\t");
		for(const auto col : iter::range(toks.size())){
			outputProperPairs.content_[regRowCount][col] = toks[col];
		}
		regUidToRowPos[reg.createUidFromCoords()] = regRowCount;
		++regRowCount;
	}

	table outputTotalPairs(header);
	outputTotalPairs.content_ = std::vector<std::vector<std::string> >(regions.size(), std::vector<std::string>(header.size()));
	regRowCount = 0;
	for(const auto & reg : regions){
		auto regAsBed = reg.genBedRecordCore();
		auto toks = tokenizeString(regAsBed.toDelimStr(), "\t");
		for(const auto col : iter::range(toks.size())){
			outputTotalPairs.content_[regRowCount][col] = toks[col];
		}
		++regRowCount;
	}

	table outputMateUnmappedPairs(header);
	outputMateUnmappedPairs.content_ = std::vector<std::vector<std::string> >{regions.size(), std::vector<std::string>{header.size()}};
	regRowCount = 0;
	for(const auto & reg : regions){
		auto regAsBed = reg.genBedRecordCore();
		auto toks = tokenizeString(regAsBed.toDelimStr(), "\t");
		for(const auto col : iter::range(toks.size())){
			outputMateUnmappedPairs.content_[regRowCount][col] = toks[col];
		}
		++regRowCount;
	}

	table outputDiscordantPairs(header);
	outputDiscordantPairs.content_ = std::vector<std::vector<std::string> >{regions.size(), std::vector<std::string>{header.size()}};
	regRowCount = 0;
	for(const auto & reg : regions){
		auto regAsBed = reg.genBedRecordCore();
		auto toks = tokenizeString(regAsBed.toDelimStr(), "\t");
		for(const auto col : iter::range(toks.size())){
			outputDiscordantPairs.content_[regRowCount][col] = toks[col];
		}
		++regRowCount;
	}

	pairsList.reset();
	std::function<void()> fillTables = [&outputTotalPairs,&outputProperPairs,&outputMateUnmappedPairs,&outputDiscordantPairs, &pairsList, &regUidToRowPos, &bamFnpToCol](){
		std::shared_ptr<BamFnpRegionPair> val;
		while(pairsList.getVal(val)){
			outputTotalPairs.content_[regUidToRowPos[val->regUid_]][bamFnpToCol[val->bamFname_]] = estd::to_string(val->totalPairCnt_);
			outputProperPairs.content_[regUidToRowPos[val->regUid_]][bamFnpToCol[val->bamFname_]] = estd::to_string(val->properPairCnt_);
			outputMateUnmappedPairs.content_[regUidToRowPos[val->regUid_]][bamFnpToCol[val->bamFname_]] = estd::to_string(val->mateUnmappedCnt_);
			outputDiscordantPairs.content_[regUidToRowPos[val->regUid_]][bamFnpToCol[val->bamFname_]] = estd::to_string(val->discordantCnt_);
		}
	};

	{
		njh::concurrent::runVoidFunctionThreaded(fillTables, numThreads);
	}

	if(noHeader){
		outputTotalPairs.hasHeader_ = false;
	}
	outputTotalPairs.outPutContents(pairCountsOut, "\t");

	if(noHeader){
		outputProperPairs.hasHeader_ = false;
	}
	outputProperPairs.outPutContents(properPairCountsOut, "\t");

	if(noHeader){
		outputMateUnmappedPairs.hasHeader_ = false;
	}
	outputMateUnmappedPairs.outPutContents(mateUnMappedPairCountsOut, "\t");

	if(noHeader){
		outputDiscordantPairs.hasHeader_ = false;
	}
	outputDiscordantPairs.outPutContents(discordantCountsOut, "\t");

	return 0;
}



int bamExpRunner::bamCov(const njh::progutils::CmdArgs & inputCommands){
	OutOptions outOpts(bfs::path(""));
	uint32_t numThreads = 1;
	bfs::path bedFnp = "";
	bfs::path bamFnp = "";
	bool noHeader = false;
	outOpts.outExtention_ = ".tab.txt";
	bool countDups = false;
	uint32_t mapQualityCutOff = 20;
	seqSetUp setUp(inputCommands);
	setUp.description_ = "Get the coverage in base count for bam files for certain regions";
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(mapQualityCutOff, "--mapQualityCutOff",   "Only reads that are this mapping quality and above (inclusive)");
	setUp.setOption(countDups, "--countDups",   "Count records marked duplicate");
	setUp.setOption(numThreads, "--numThreads", "Number of threads to use");
	setUp.processWritingOptions(outOpts);
	setUp.setOption(bedFnp, "--bedFnp", "Bed file of regions to get coverage for", true, "Input");
	setUp.setOption(bamFnp, "--bamFnp", "Bam Input file", false, "Input");
	setUp.setOption(noHeader, "--noHeader", "Don't output a header so output can be treated like a bed file", false, "Writing Output");
	setUp.finishSetUp(std::cout);


	OutputStream out(outOpts);
	checkBamFilesForIndexesAndAbilityToOpen({bamFnp});
	auto bed3s = getBed3s(bedFnp);
	std::vector<GenomicRegion> inputRegions;
	bool allAbove3 = true;
	for(const auto & b : bed3s){
		if(b->extraFields_.size() < 3){
			allAbove3 = false;
			break;
		}
	}
	if (allAbove3) {
		inputRegions = bedPtrsToGenomicRegs(getBeds(bedFnp));
	} else {
		inputRegions = bed3PtrsToGenomicRegs(bed3s);
	}

	//collapse identical regions
	std::vector<GenomicRegion> regions;
	for(const auto & inputRegion : inputRegions){
		if(regions.empty()){
			regions.emplace_back(inputRegion);
		}else{
			if(regions.back().sameRegion(inputRegion)){
				regions.back().uid_ += "," + inputRegion.uid_;
			}else{
				regions.emplace_back(inputRegion);
			}
		}
	}



	out << "#chrom\tstart\tend\tname\tscore\tstrand\t"<< bamFnp << std::endl;

	njh::concurrent::LockableVec<GenomicRegion> regionsQueue(regions);
	std::mutex outMut;
	std::function<void()> getCov = [&regionsQueue,&countDups,&mapQualityCutOff,&bamFnp,&outMut,&out](){
		GenomicRegion val;
		while(regionsQueue.getVal(val)){

			BamTools::BamReader bReader;
			bReader.Open(bamFnp.string());
			bReader.LocateIndex();
			setBamFileRegionThrow(bReader, val);
			auto refData = bReader.GetReferenceData();
			BamTools::BamAlignment bAln;
			uint64_t coverage = 0;
			while(bReader.GetNextAlignmentCore(bAln)){
				if(bAln.IsMapped() && bAln.IsPrimaryAlignment()){
					if(bAln.IsDuplicate() && !countDups){
						continue;
					}
					if(bAln.MapQuality < mapQualityCutOff){
						continue;
					}
					/**@todo this doesn't take into account gaps, so base coverage isn't precise right here, more of an approximation, should improve */
					coverage += 1;
				}
			}
			{
				std::lock_guard<std::mutex> lockQuard(outMut);
				out << val.genBedRecordCore().toDelimStr() << "\t" << coverage << std::endl;
			}
		}
	};

	{
		njh::concurrent::runVoidFunctionThreaded(getCov, numThreads);
	}
	return 0;
}


int bamExpRunner::bamCovBasesRough(const njh::progutils::CmdArgs & inputCommands){
	OutOptions outOpts(bfs::path(""));
	uint32_t numThreads = 1;
	bfs::path bedFnp = "";
	bfs::path bamFnp = "";
	bool noHeader = false;
	outOpts.outExtention_ = ".tab.txt";
	bool countDups = false;
	uint32_t mapQualityCutOff = 20;
	seqSetUp setUp(inputCommands);
	setUp.description_ = "Get the coverage in base count for bam files for certain regions";
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(mapQualityCutOff, "--mapQualityCutOff",   "Only reads that are this mapping quality and above (inclusive)");
	setUp.setOption(countDups, "--countDups",   "Count records marked duplicate");
	setUp.setOption(numThreads, "--numThreads", "Number of threads to use");
	setUp.processWritingOptions(outOpts);
	setUp.setOption(bedFnp, "--bedFnp", "Bed file of regions to get coverage for", true, "Input");
	setUp.setOption(bamFnp, "--bamFnp", "Bam Input file", false, "Input");
	setUp.setOption(noHeader, "--noHeader", "Don't output a header so output can be treated like a bed file", false, "Writing Output");
	setUp.finishSetUp(std::cout);


	OutputStream out(outOpts);
	checkBamFilesForIndexesAndAbilityToOpen({bamFnp});
	auto bed3s = getBed3s(bedFnp);
	std::vector<GenomicRegion> inputRegions;
	bool allAbove3 = true;
	for(const auto & b : bed3s){
		if(b->extraFields_.size() < 3){
			allAbove3 = false;
			break;
		}
	}
	if (allAbove3) {
		inputRegions = bedPtrsToGenomicRegs(getBeds(bedFnp));
	} else {
		inputRegions = bed3PtrsToGenomicRegs(bed3s);
	}

	//collapse identical regions
	std::vector<GenomicRegion> regions;
	for(const auto & inputRegion : inputRegions){
		if(regions.empty()){
			regions.emplace_back(inputRegion);
		}else{
			if(regions.back().sameRegion(inputRegion)){
				regions.back().uid_ += "," + inputRegion.uid_;
			}else{
				regions.emplace_back(inputRegion);
			}
		}
	}



	out << "#chrom\tstart\tend\tname\tscore\tstrand\t"<< bamFnp << std::endl;

	njh::concurrent::LockableVec<GenomicRegion> regionsQueue(regions);
	std::mutex outMut;
	std::function<void()> getCov = [&regionsQueue,&countDups,&mapQualityCutOff,&bamFnp,&outMut,&out](){
		GenomicRegion val;
		while(regionsQueue.getVal(val)){

			BamTools::BamReader bReader;
			bReader.Open(bamFnp.string());
			bReader.LocateIndex();
			setBamFileRegionThrow(bReader, val);
			auto refData = bReader.GetReferenceData();
			BamTools::BamAlignment bAln;
			uint64_t coverage = 0;
			while(bReader.GetNextAlignmentCore(bAln)){
				if(bAln.IsMapped() && bAln.IsPrimaryAlignment()){
					if(bAln.IsDuplicate() && !countDups){
						continue;
					}
					if(bAln.MapQuality < mapQualityCutOff){
						continue;
					}
					/**@todo this doesn't take into account gaps, so base coverage isn't precise right here, more of an approximation, should improve */
					coverage += val.getOverlapLen(refData[bAln.RefID].RefName, bAln.Position, bAln.GetEndPosition());
				}
			}
			{
				std::lock_guard<std::mutex> lockQuard(outMut);
				out << val.genBedRecordCore().toDelimStr() << "\t" << coverage << std::endl;
			}
		}
	};

	{
		njh::concurrent::runVoidFunctionThreaded(getCov, numThreads);
	}
	return 0;
}




int bamExpRunner::bamMulticovBasesRough(const njh::progutils::CmdArgs & inputCommands){
	OutOptions outOpts(bfs::path(""));
	uint32_t numThreads = 1;
	bfs::path bedFnp = "";
	std::string bams;
	std::string pat = ".*.bam$";
	bool noHeader = false;
	outOpts.outExtention_ = ".tab.txt";
	bool countDups = false;
	uint32_t mapQualityCutOff = 20;
	seqSetUp setUp(inputCommands);
	setUp.description_ = "Get the coverage in base count for bam files for certain regions";
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(mapQualityCutOff, "--mapQualityCutOff",   "Only reads that are this mapping quality and above (inclusive)");
	setUp.setOption(countDups, "--countDups",   "Count records marked duplicate");
	setUp.setOption(numThreads, "--numThreads", "Number of threads to use");
	setUp.processWritingOptions(outOpts);
	setUp.setOption(bedFnp, "--bedFnp", "Bed file of regions to get coverage for", true, "Input");
	setUp.setOption(pat, "--pat", "Pattern in current directory to get coverage for", false, "Input");
	setUp.setOption(bams, "--bams", "Either a file with the name of a bam file on each line or a comma separated value of bam file paths", false, "Input");
	setUp.setOption(noHeader, "--noHeader", "Don't output a header so output can be treated like a bed file", false, "Writing Output");
	setUp.finishSetUp(std::cout);


	OutputStream outFile(outOpts);

	auto bamFnps = njh::files::gatherFilesByPatOrNames(std::regex{pat}, bams);
	checkBamFilesForIndexesAndAbilityToOpen(bamFnps);
	if(setUp.pars_.verbose_){
		printVector(bamFnps, "\t", std::cout);
	}
	auto bed3s = getBed3s(bedFnp);
	std::vector<GenomicRegion> inputRegions;
	bool allAbove3 = true;
	for(const auto & b : bed3s){
		if(b->extraFields_.size() < 3){
			allAbove3 = false;
			break;
		}
	}
	if (allAbove3) {
		inputRegions = bedPtrsToGenomicRegs(getBeds(bedFnp));
	} else {
		inputRegions = bed3PtrsToGenomicRegs(bed3s);
	}

	//collapse identical regions
	std::vector<GenomicRegion> regions;
	for(const auto & inputRegion : inputRegions){
		if(regions.empty()){
			regions.emplace_back(inputRegion);
		}else{
			if(regions.back().sameRegion(inputRegion)){
				regions.back().uid_ += "," + inputRegion.uid_;
			}else{
				regions.emplace_back(inputRegion);
			}
		}
	}

	struct BamFnpRegionPair {
		BamFnpRegionPair(const bfs::path & bamFnp, const GenomicRegion & region) :
				bamFnp_(bamFnp), region_(region) {
			regUid_ = region_.createUidFromCoords();
			bamFname_ = bamFnp_.filename().string();
		}
		bfs::path bamFnp_;
		GenomicRegion region_;
		uint32_t coverage_ = 0;

		std::string regUid_;
		std::string bamFname_;
	};

	std::vector<std::shared_ptr<BamFnpRegionPair>> pairs;
	for(const auto & bamFnp : bamFnps){
		for(const auto & region : regions){
			pairs.emplace_back(std::make_shared<BamFnpRegionPair>(bamFnp, region));
		}
	}

	njh::concurrent::LockableVec<std::shared_ptr<BamFnpRegionPair>> pairsList(pairs);

	std::function<void()> getCov = [&pairsList,&countDups,&mapQualityCutOff](){
		std::shared_ptr<BamFnpRegionPair> val;
		while(pairsList.getVal(val)){
			BamTools::BamReader bReader;
			bReader.Open(val->bamFnp_.string());
			bReader.LocateIndex();
			setBamFileRegionThrow(bReader, val->region_);
			auto refData = bReader.GetReferenceData();
			BamTools::BamAlignment bAln;
			while(bReader.GetNextAlignmentCore(bAln)){
				if(bAln.IsMapped() && bAln.IsPrimaryAlignment()){
					if(bAln.IsDuplicate() && !countDups){
						continue;
					}
					if(bAln.MapQuality < mapQualityCutOff){
						continue;
					}
					/**@todo this doesn't take into account gaps, so base coverage isn't precise right here, more of an approximation, should improve */
					GenomicRegion alnRegion(bAln, refData);
					val->coverage_ += val->region_.getOverlapLen(alnRegion);
				}
			}
		}
	};

	{
		njh::concurrent::runVoidFunctionThreaded(getCov, numThreads);
	}

	VecStr header = {"chrom", "start", "end", "name", "score", "strand"};
	std::unordered_map<std::string, uint32_t> bamFnpToCol;
	uint32_t bIndex = 6;
	for(const auto & b : bamFnps){
		header.emplace_back(b.filename().string());
		bamFnpToCol[b.filename().string()] = bIndex;
		++bIndex;
	}

	njh::sort(regions,
			[](const GenomicRegion & reg1,
					const GenomicRegion & reg2) {
				return reg1.createUidFromCoords() < reg2.createUidFromCoords();
			});
	table output(header);
	output.content_ = std::vector<std::vector<std::string>>{regions.size(), std::vector<std::string>{header.size()}};
	uint32_t regRowCount = 0;
	std::unordered_map<std::string, uint32_t> regUidToRowPos;
	for(const auto & reg : regions){
		auto regAsBed = reg.genBedRecordCore();
		auto toks = tokenizeString(regAsBed.toDelimStr(), "\t");
		for(const auto col : iter::range(toks.size())){
			output.content_[regRowCount][col] = toks[col];
		}
		regUidToRowPos[reg.createUidFromCoords()] = regRowCount;
		++regRowCount;
	}
	pairsList.reset();
	std::function<void()> fillTable = [&output, &pairsList, &regUidToRowPos, &bamFnpToCol](){
		std::shared_ptr<BamFnpRegionPair> val;
		while(pairsList.getVal(val)){
			output.content_[regUidToRowPos[val->regUid_]][bamFnpToCol[val->bamFname_]] = estd::to_string(val->coverage_);
		}
	};

	{
		njh::concurrent::runVoidFunctionThreaded(fillTable, numThreads);
	}
	if(noHeader){
		output.hasHeader_ = false;
	}
	output.outPutContents(outFile, "\t");

	return 0;
}



int bamExpRunner::bamToBed(const njh::progutils::CmdArgs & inputCommands){
	bool includeUnmapped = false;
	bool onlyPrimary = false;
	bool appendCigarString = false;
	uint32_t maxAlnSize = std::numeric_limits<uint32_t>::max();
	uint32_t minAlnSize = 0;

	seqSetUp setUp(inputCommands);
	setUp.pars_.ioOptions_.out_.outFilename_ = "";
	setUp.processDefaultReader(VecStr{"--bam"});
	setUp.pars_.ioOptions_.out_.outExtention_ = ".bed";
	setUp.pars_.ioOptions_.out_.outFileFormat_ = "bed";
	setUp.setOption(includeUnmapped, "--includeUnmapped", "Include Unmapped");
	setUp.setOption(onlyPrimary, "--onlyPrimary", "Only Primary");
	setUp.setOption(appendCigarString, "--appendCigarString", "Append Cigar String");
	setUp.setOption(maxAlnSize, "--maxAlnSize", "maximum Alignment Size");
	setUp.setOption(minAlnSize, "--minAlnSize", "minimum Alignment Size");

	setUp.finishSetUp(std::cout);

	BamTools::BamReader bReader;
	bReader.Open(setUp.pars_.ioOptions_.firstName_.string());
	checkBamOpenThrow(bReader, setUp.pars_.ioOptions_.firstName_.string());

	OutputStream out(setUp.pars_.ioOptions_.out_);

	BamTools::BamAlignment bAln;
	auto refData = bReader.GetReferenceData();
	while(bReader.GetNextAlignment(bAln)){
		if(bAln.IsMapped()){
			if(!bAln.IsPrimaryAlignment() && onlyPrimary){
				continue;
			}
			std::string app = "";
			if(appendCigarString){
				app += "-";
				for(const auto & cigar : bAln.CigarData){
					app += njh::pasteAsStr(cigar.Type, cigar.Length);
				}
			}
			auto alnSize = bAln.GetEndPosition() - bAln.Position ;
			if(alnSize < maxAlnSize && alnSize > minAlnSize){
				out << refData[bAln.RefID].RefName
						<< "\t" << bAln.Position
						<< "\t" << bAln.GetEndPosition()
						<< "\t" << bAln.Name << app
						<< "\t" << bAln.GetEndPosition() - bAln.Position
						<< "\t" << (bAln.IsReverseStrand() ? '-' : '+') << std::endl;
			}
		} else if (includeUnmapped) {
			out << "*"
					<< "\t" << "*"
					<< "\t" << "*"
					<< "\t" << bAln.Name
					<< "\t" << bAln.QueryBases.size()
					<< "\t" << (bAln.IsReverseStrand() ? '-' : '+') << std::endl;
		}
	}

	return 0;
}

int bamExpRunner::mergeBedRegionFromBam(const njh::progutils::CmdArgs & inputCommands){
	//only meant for very small bam files
	seqSetUp setUp(inputCommands);
	setUp.pars_.ioOptions_.out_.outFilename_ = "";
	setUp.processDefaultReader(VecStr{"--bam"});
	setUp.pars_.ioOptions_.out_.outExtention_ = ".bed";
	setUp.pars_.ioOptions_.out_.outFileFormat_ = "bed";
	setUp.finishSetUp(std::cout);

	BamTools::BamReader bReader;
	bReader.Open(setUp.pars_.ioOptions_.firstName_.string());
	checkBamOpenThrow(bReader, setUp.pars_.ioOptions_.firstName_.string());

	std::ofstream outBed;
	std::ostream out(determineOutBuf(outBed, setUp.pars_.ioOptions_.out_));

	BamTools::BamAlignment bAln;
	auto refData = bReader.GetReferenceData();
	std::unordered_map<std::string, std::unordered_map<uint32_t, uint32_t>> coverage;
	std::vector<GenomicRegion> regionsRev;
	std::vector<GenomicRegion> regionsFor;
	while(bReader.GetNextAlignment(bAln)){
		if(bAln.IsMapped()){
			if(bAln.IsReverseStrand()){
				regionsRev.emplace_back(bAln, refData);
			}else{
				regionsFor.emplace_back(bAln, refData);
			}
		}
	}

	njh::sort(regionsRev, [](const GenomicRegion & reg1, const GenomicRegion & reg2){
		if(reg1.chrom_ == reg2.chrom_){
			return reg1.start_ < reg2.start_;
		}else{
			return reg1.chrom_ < reg2.chrom_;
		}
	});

	njh::sort(regionsFor, [](const GenomicRegion & reg1, const GenomicRegion & reg2){
		if(reg1.chrom_ == reg2.chrom_){
			return reg1.start_ < reg2.start_;
		}else{
			return reg1.chrom_ < reg2.chrom_;
		}
	});

	if(regionsFor.size() == 1){
		out << regionsFor.front().genBedRecordCore().toDelimStr() << std::endl;
	}else if(regionsFor.size() > 1){
		GenomicRegion currentRegion = regionsFor.front();
		for(const auto pos : iter::range<uint32_t>(1, regionsFor.size()) ){
			if(currentRegion.overlaps(regionsFor[pos])){
				currentRegion.start_ = std::min(currentRegion.start_, regionsFor[pos].start_);
				currentRegion.end_ = std::max(currentRegion.end_, regionsFor[pos].end_);
				currentRegion.alternateUids_.emplace_back(regionsFor[pos].uid_);
			}else{
				out << currentRegion.genBedRecordCore().toDelimStr() << std::endl;
				currentRegion = regionsFor[pos];
			}
		}
		out << currentRegion.genBedRecordCore().toDelimStr() << std::endl;
	}

	if(regionsRev.size() == 1){
		out << regionsRev.front().genBedRecordCore().toDelimStr() << std::endl;
	}else if(regionsRev.size() > 1){
		GenomicRegion currentRegion = regionsRev.front();
		for(const auto pos : iter::range<uint32_t>(1, regionsRev.size()) ){
			if(currentRegion.overlaps(regionsRev[pos])){
				currentRegion.start_ = std::min(currentRegion.start_, regionsRev[pos].start_);
				currentRegion.end_ = std::max(currentRegion.end_, regionsRev[pos].end_);
				currentRegion.alternateUids_.emplace_back(regionsRev[pos].uid_);
			}else{
				out << currentRegion.genBedRecordCore().toDelimStr() << std::endl;
				currentRegion = regionsRev[pos];
			}
		}
		out << currentRegion.genBedRecordCore().toDelimStr() << std::endl;
	}


	return 0;
}


int bamExpRunner::getBestBedRegionFromBam(const njh::progutils::CmdArgs & inputCommands){
	//only meant for very small bam files
	seqSetUp setUp(inputCommands);
	setUp.pars_.ioOptions_.out_.outFilename_ = "";
	setUp.processDefaultReader(VecStr{"--bam"});
	setUp.pars_.ioOptions_.out_.outExtention_ = ".bed";
	setUp.pars_.ioOptions_.out_.outFileFormat_ = "bed";
	setUp.finishSetUp(std::cout);

	BamTools::BamReader bReader;
	bReader.Open(setUp.pars_.ioOptions_.firstName_.string());
	checkBamOpenThrow(bReader, setUp.pars_.ioOptions_.firstName_.string());

	std::ofstream outBed;
	std::ostream out(determineOutBuf(outBed, setUp.pars_.ioOptions_.out_));

	BamTools::BamAlignment bAln;
	auto refData = bReader.GetReferenceData();
	std::unordered_map<std::string, std::unordered_map<uint32_t, uint32_t>> coverage;
	while(bReader.GetNextAlignment(bAln)){

		if(bAln.IsMapped()){
			for(const auto pos : iter::range(bAln.Position, bAln.GetEndPosition())){
				++coverage[refData[bAln.RefID].RefName][pos];
			}
		}
	}

	for(const auto & chrom : coverage){
		uint32_t bestCoverage = 0;
		for(const auto & cov : chrom.second){
			if(cov.second > bestCoverage){
				bestCoverage = cov.second;
			}
		}
		std::vector<uint32_t> positionsWithBestCoverage;
		for(const auto & cov : chrom.second){
			if(cov.second  == bestCoverage){
				positionsWithBestCoverage.emplace_back(cov.first);
			}
		}
		auto start = *std::min_element(positionsWithBestCoverage.begin(), positionsWithBestCoverage.end());
		auto end = (*std::max_element(positionsWithBestCoverage.begin(), positionsWithBestCoverage.end())) + 1;
		out << chrom.first
				<< "\t" << start
				<< "\t" << end
				<< "\t" << chrom.first << "-" << start << "-" << end
				<< "\t" << end - start
				<< "\t" << "+" << std::endl;
	}
	return 0;
}





int bamExpRunner::BamExtractReadsFromRegion(
		const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bedFile = "";
	double percInRegion = .50;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(bedFile, "--bed", "Bed file, first entry is used", true);
	setUp.processDefaultReader( { "-bam" }, true);
	setUp.setOption(percInRegion, "--percInRegion", "Percent of bases n Region");

	setUp.finishSetUp(std::cout);
	BamExtractor bExtractor(setUp.pars_.verbose_);

	//auto regions = bedPtrsToGenomicRegs(getBeds(bedFile));
	auto regions = gatherRegions(bedFile.string(), "", setUp.pars_.verbose_);
	if(1 == regions.size()){
		bExtractor.writeExtractReadsFromBamRegion(setUp.pars_.ioOptions_.firstName_,
				regions.front(), percInRegion, setUp.pars_.ioOptions_.out_);
	}else if(regions.empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error no regions read in from "  << bedFile << std::endl;
	}else{
		std::unordered_map<std::string, uint32_t> regionNameCounts;
		for(const auto & region : regions){
			if(setUp.pars_.verbose_){
				std::cout << "Extracting: " << region.uid_ << std::endl;
			}
			std::string name = region.uid_;
			if(njh::in(region.uid_, regionNameCounts)){
				name += "_" + estd::to_string(regionNameCounts[region.uid_]);
			}
			++regionNameCounts[region.uid_];
			OutOptions regOutOpts(bfs::path(setUp.pars_.ioOptions_.out_.outFilename_.string() + "_" + name));
			regOutOpts.transferOverwriteOpts(setUp.pars_.ioOptions_.out_);
			bExtractor.writeExtractReadsFromBamRegion(setUp.pars_.ioOptions_.firstName_,
							region, percInRegion, regOutOpts);
		}
	}


	return 0;
}


int bamExpRunner::getInsertSizesStats(const njh::progutils::CmdArgs & inputCommands) {
	uint32_t testNum = std::numeric_limits<uint32_t>::max();
	uint32_t numThreads = 1;
	OutOptions outOpts(bfs::path(""), ".tab.txt");
	bfs::path bedFnp = "";
	uint32_t hardInsertSizeCutOff = 10000;
	uint32_t mapQualityCutOff = 20;
	GenomeRegionsGenerator::Params genPars;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.setOption(hardInsertSizeCutOff, "--hardInsertSizeCutOff", "Hard Insert Size Cut Off");
	setUp.setOption(testNum, "--testNum", "Test Number");
	setUp.processReadInNames(VecStr { "--bam" }, true);
	setUp.processWritingOptions(outOpts);
	setUp.setOption(bedFnp, "--bed", "Check just for these regions");
	setUp.setOption(mapQualityCutOff, "--mapQualityCutOff", "Mapping Quality Cut Off");
	setUp.setOption(genPars.step_, "--windowStep", "Window Step", njh::progutils::ProgramSetUp::CheckCase::NONZERO);
	setUp.setOption(genPars.window_, "--windowSize", "Window Size", njh::progutils::ProgramSetUp::CheckCase::NONZERO);
	setUp.setOption(numThreads, "--numThreads", "Number of Threads", njh::progutils::ProgramSetUp::CheckCase::NONZERO);


	setUp.finishSetUp(std::cout);


	std::unordered_map<std::string, std::vector<uint32_t>> insertSizes;
	std::mutex insertSizesMut;
	OutputStream out(outOpts);
	concurrent::BamReaderPool bPool(setUp.pars_.ioOptions_.firstName_, numThreads);
	bPool.openBamFile();
	std::vector<GenomicRegion> regions;
	if("" != bedFnp){
		regions =	bed3PtrsToGenomicRegs(getBed3s(bedFnp));
	}else{
		auto bReader = bPool.popReader();
		regions = genChromosomeGenRegions(bReader->GetReferenceData());
	}

	GenomeRegionsGenerator regionFactory(regions, genPars);
	std::atomic_uint32_t totalCount{0};

	std::function<void()> getInsertSizes = [&insertSizesMut,&insertSizes,&bPool,&regionFactory,&hardInsertSizeCutOff,&mapQualityCutOff,&totalCount,&testNum](){
		auto bReader = bPool.popReader();
		GenomicRegion reg;
		std::unordered_map<std::string, std::vector<uint32_t>> currentInsertSizes;
		BamTools::BamAlignment aln;
		while(regionFactory.genRegion(reg)){
			setBamFileRegionThrow(*bReader, reg);
			while (bReader->GetNextAlignmentCore(aln)) {
				//skip alignments that don't start in this region
				//this way if regions are close to each other it will avoid counting twice
				if(aln.Position <reg.start_ || aln.Position >= reg.end_){
					continue;
				}
				if (aln.IsMapped() &&
						aln.IsMateMapped() &&
						aln.IsPrimaryAlignment() &&
						aln.IsReverseStrand() != aln.IsMateReverseStrand()) {
					if (aln.RefID == aln.MateRefID) {
						if(std::abs(aln.InsertSize) > hardInsertSizeCutOff || aln.MapQuality < mapQualityCutOff){
							continue;
						}
						++totalCount;
						std::string readGroupName = "none";
						aln.BuildCharData();
						aln.GetTag("RG", readGroupName);
						currentInsertSizes[readGroupName].emplace_back(std::abs(aln.InsertSize));
					}
				}
			}
			if(totalCount.load() > testNum){
				break;
			}
		}
		{
			std::lock_guard<std::mutex> lock(insertSizesMut);
			for(const auto & rg : currentInsertSizes){
				addOtherVec(insertSizes[rg.first], rg.second);
			}
		}
	};
	njh::concurrent::runVoidFunctionThreaded(getInsertSizes, numThreads);



	table outStatTable(VecStr {"ReadGroup", "stat", "value" });
	for( auto & rg : insertSizes){
		if(rg.second.empty()){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "read group: " << rg.first << " had no insert sizes" << "\n";
			throw std::runtime_error{ss.str()};
		}
		njh::sort(rg.second);
		auto std = vectorStandardDeviationPop(rg.second);
		auto max = rg.second.back();
		auto min = rg.second.front();
		double sum = 0;
		for(const auto val : rg.second){
			sum += val;
		}
		auto mean = sum/rg.second.size();
		double median = rg.second[1 + (rg.second.size()/2)];
		if(0 == rg.second.size() % 2){
			median += rg.second[rg.second.size()/2];
			median /=2;
		}

		table rgStatTable(VecStr {"ReadGroup", "stat", "value" });
		rgStatTable.addRow(rg.first,"max", max);
		rgStatTable.addRow(rg.first,"mean", mean);
		rgStatTable.addRow(rg.first,"median", median);
		rgStatTable.addRow(rg.first,"std", std);
		rgStatTable.addRow(rg.first,"min", min);
		outStatTable.rbind(rgStatTable, false);
	}
	outStatTable.sortTable("ReadGroup", "stat", false);
	outStatTable.outPutContents(out, "\t");

	return 0;
}


int bamExpRunner::printBamRefIds(const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processVerbose();
	setUp.processDefaultReader( { "--bam" }, true);

	setUp.finishSetUp(std::cout);

	BamTools::BamReader bReader;
	bReader.Open(setUp.pars_.ioOptions_.firstName_.string());
	checkBamOpenThrow(bReader,setUp.pars_.ioOptions_.firstName_.string());
	auto refData = bReader.GetReferenceData();

	for(const auto & ref : refData){
		std::cout << ref.RefName << "\t"<< ref.RefLength << std::endl;
	}

	return 0;
}

int bamExpRunner::BamRefIdsToBed(const njh::progutils::CmdArgs & inputCommands) {
	OutOptions outOpts(bfs::path(""), ".bed");
	bfs::path bamFnp = "";
	seqSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processVerbose();
	setUp.setOption(bamFnp, "--bam", "bam file", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	BamTools::BamReader bReader;
	bReader.Open(bamFnp.string());
	checkBamOpenThrow(bReader,bamFnp.string());
	auto refData = bReader.GetReferenceData();
	OutputStream out(outOpts);
	for(const auto & ref : refData){
		out << ref.RefName << "\t" << "0" << "\t" << ref.RefLength << std::endl;
	}
	return 0;
}

int bamExpRunner::bamToFastq(const njh::progutils::CmdArgs & inputCommands) {
	bool onlyMapped = false;

	std::unordered_set<std::string> names;
	seqSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processVerbose();
	setUp.setOption(onlyMapped, "--onlyMapped", "Extract only mapped alignments");
	setUp.processDefaultReader( { "--bam" }, true);
	setUp.setOption(names, "--names", "Extract only alignments with these names");

	setUp.finishSetUp(std::cout);
	BamExtractor bExtractor(setUp.pars_.verbose_);
	if (onlyMapped) {
		if(!names.empty()){
			bExtractor.writeExtractReadsFromBamOnlyMapped(setUp.pars_.ioOptions_.firstName_, setUp.pars_.ioOptions_.out_, names);
		}else{
			bExtractor.writeExtractReadsFromBamOnlyMapped(setUp.pars_.ioOptions_.firstName_, setUp.pars_.ioOptions_.out_);
		}
	} else {
		if(!names.empty()){
			bExtractor.writeExtractReadsFromBam(setUp.pars_.ioOptions_.firstName_, setUp.pars_.ioOptions_.out_, names);
		}else{
			bExtractor.writeExtractReadsFromBam(setUp.pars_.ioOptions_.firstName_, setUp.pars_.ioOptions_.out_);
		}
	}
	return 0;
}

namespace bfs = boost::filesystem;

bool checkLine(std::istream & inFile, size_t pos, std::unordered_map<std::string, uint32_t>& colNamePos){
	auto nextLine = njh::files::peekLine(inFile);
	if("" != nextLine){
		auto nextLineToks = tokenizeString(nextLine, "\t");
		if(njh::lexical_cast<size_t>(nextLineToks[colNamePos["refPos"]]) == pos){
			return true;
		}
	}
	return false;
}

table processBaseCounts(const std::string & filename,
		const std::string & sampleName, const std::string fracColname,
		double depthCutOff, double fracCutOff,
		double strandBias){

	std::string percentStrandColname = "percForward";
	std::string totalColname = "totalCount";
	if(njh::endsWith(fracColname, "HighQ")){
		percentStrandColname.append("HighQ");
		totalColname.append("HighQ");
	}

	std::string year = "NA";
	std::string grouping = "NA";
	if(sampleName.size() > 2){
		if(isalpha(sampleName[0]) && isalpha(sampleName[1])){
			grouping = sampleName.substr(0,2);
		}
		if(sampleName.size() == 6){
			if(isdigit(sampleName[2])
					&& isdigit(sampleName[3])
					&& isdigit(sampleName[4])
					&& isdigit(sampleName[5])){
				year = sampleName.substr(2,4);
			}
		}
	}
	auto firstLine = njh::files::getFirstLine(filename);
	auto toks = tokenizeString(firstLine, "\t");
	std::unordered_map<std::string, uint32_t> colNamePos;
	for(const auto tokPos : iter::range(toks.size())){
		colNamePos[toks[tokPos]] = tokPos;
	}
	table outTab(toks);
	std::string line;

	std::ifstream inFile(filename);
	std::map<size_t, VecStr> lines;
	size_t currentPos = std::numeric_limits<size_t>::max();
	while(njh::files::crossPlatGetline(inFile, line)){
		if(njh::beginsWith(line, "refName")){
			continue;
		}
		if(lines.empty()){
			auto lineToks = tokenizeString(line, "\t", true);
			size_t depth = njh::lexical_cast<size_t>(lineToks[colNamePos[totalColname]]);
			if(depth < depthCutOff){
				continue;
			}
			currentPos = njh::lexical_cast<size_t>(lineToks[colNamePos["refPos"]]);
			lines[lines.size()] = lineToks;
			while(checkLine(inFile, currentPos, colNamePos)){
				njh::files::crossPlatGetline(inFile, line);
				auto lineToks = tokenizeString(line, "\t", true);
				lines[lines.size()] = lineToks;
			}
			std::vector<double> fracs;
			for (const auto & l : lines) {
				std::string percentForwardStr =
						l.second[colNamePos[percentStrandColname]];
				double percentForward = 1;
				if (!njh::containsSubString(stringToLowerReturn(percentForwardStr),
						"nan")) {
					percentForward = njh::lexical_cast<double>(
							l.second[colNamePos[percentStrandColname]]);
				}
				if (percentForward < (1 - strandBias) && percentForward > strandBias) {
					fracs.emplace_back(
							njh::lexical_cast<double>(l.second[colNamePos[fracColname]]));
				}
			}
			uint32_t numAboveCutOff = 0;
			for(const auto & frac : fracs){
				if(frac >= fracCutOff){
					++numAboveCutOff;
					if(numAboveCutOff > 1){
						break;
					}
				}
			}
			if(numAboveCutOff > 1){
				for(const auto & l : lines){
					outTab.content_.emplace_back(l.second);
				}
			}
			lines.clear();
		}
	}
	outTab.addColumn(VecStr{sampleName}, "sampleName");
	outTab.addColumn(VecStr{grouping}, "grouping");
	outTab.addColumn(VecStr{year}, "year");
	outTab = outTab.getColumns(concatVecs({"sampleName", "year", "grouping"}, toks));
	return outTab;
}



int bamExpRunner::bamReadProfiling(const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	//std::string bedFile = "";
	uint32_t mappingQual = 10;
	//setUp.setOption(bedFile, "--bedFile", "Bed file to print stats on these regions", true);
	setUp.processDefaultReader(VecStr{"-bam"}, true);
	setUp.processDirectoryOutputName(true);
	setUp.setOption(mappingQual, "--mappingQual", "Mapping quality score cut off for a read");
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);

	BamTools::BamReader bReader;
	BamTools::BamAlignment bAln;
	bReader.Open(setUp.pars_.ioOptions_.firstName_.string());
	if (!bReader.IsOpen()) {
		std::cerr << "Error in opening " << setUp.pars_.ioOptions_.firstName_
				<< std::endl;
		exit(1);
	}
	auto refs = bReader.GetReferenceData();
	uint64_t totalReads = 0;
	uint64_t readCount = 0;
	uint32_t readsNotMapped = 0;
	uint32_t readsBellowMappingQuality = 0;
	std::unordered_map<std::string,
			std::unordered_map<size_t,
					std::unordered_map<size_t, std::unordered_map<bool, uint64_t>>> >refPosSizeStrandCounts;
	std::unordered_map<std::string,
			std::unordered_map<size_t, std::unordered_map<bool, uint64_t>>>refSizeStrandCounts;
	std::unordered_map<std::string,
			std::unordered_map<size_t, std::unordered_map<bool, uint64_t>>>refPosStrandCounts;
	std::unordered_map<std::string, std::unordered_map<bool, uint64_t>> refStrandCounts;
	std::unordered_map<std::string, std::unordered_map<uint32_t,std::unordered_map<bool, uint64_t>>> queryCoverageCount;

	while (bReader.GetNextAlignment(bAln)) {
		++totalReads;
		if (bAln.IsMapped() && bAln.MapQuality > mappingQual) {
			++readCount;
			if (setUp.pars_.verbose_) {
				std::cout << '\r' << readCount;
				std::cout.flush();
			}
			//get reference name
			auto rName = refs[bAln.RefID].RefName;
			auto seq = bAln.AlignedBases;
			removeChar(seq, '-');
			double coverage = roundDecPlaces(
					seq.size() / static_cast<double>(bAln.QueryBases.size()), 3);

			++queryCoverageCount[rName][static_cast<uint32_t>(coverage * 1000)][!bAln.IsReverseStrand()];;

			++refPosSizeStrandCounts[rName][bAln.Position][seq.size()][!bAln.IsReverseStrand()];
			++refPosStrandCounts[rName][bAln.Position][!bAln.IsReverseStrand()];
			++refSizeStrandCounts[rName][seq.size()][!bAln.IsReverseStrand()];
			++refStrandCounts[rName][!bAln.IsReverseStrand()];

			refPosSizeStrandCounts[rName][bAln.Position][seq.size()][bAln.IsReverseStrand()] += 0;
			refPosStrandCounts[rName][bAln.Position][bAln.IsReverseStrand()] += 0;
			refSizeStrandCounts[rName][seq.size()][bAln.IsReverseStrand()] += 0;
			refStrandCounts[rName][bAln.IsReverseStrand()] += 0;
			queryCoverageCount[rName][static_cast<uint32_t>(coverage * 1000)][bAln.IsReverseStrand()] += 0;

		} else {
			if (!bAln.IsMapped()) {
				++readsNotMapped;
			} else if (bAln.MapQuality <= mappingQual) {
				++readsBellowMappingQuality;
			}
		}
	}

	table sizePosCountsTab(VecStr { "ref", "pos", "size", "totalCount", "fCount",
			"rCount", "percForward" });
	table posCountsTab(VecStr { "ref", "pos", "totalCount", "fCount", "rCount",
			"percForward" });
	table sizeCountsTab(VecStr { "ref", "size", "totalCount", "fCount", "rCount",
			"percForward" });
	table refCountsTab(VecStr { "ref", "totalCount", "fCount", "rCount",
			"percForward" });
	table refQueryCoverageTab(VecStr { "ref", "coverage", "totalCount", "fCount",
			"rCount", "percForward" });

	{
		auto refKeys = getVectorOfMapKeys(refStrandCounts);
		njh::sort(refKeys);
		for (const auto & refKey : refKeys) {
			const auto & ref = refPosSizeStrandCounts[refKey];
			auto posKeys = getVectorOfMapKeys(ref);
			njh::sort(posKeys);
			for (const auto & posKey : posKeys) {
				const auto & pos = ref.at(posKey);
				auto sizeKeys = getVectorOfMapKeys(pos);
				njh::sort(sizeKeys);
				for (const auto & sizeKey : sizeKeys) {
					const auto & size = pos.at(sizeKey);
					uint64_t forward = size.at(true);
					uint64_t reverse = size.at(false);
					uint64_t total = forward + reverse;
					sizePosCountsTab.content_.emplace_back(
							toVecStr(refKey, posKey, sizeKey, total, forward, reverse,
									forward / static_cast<double>(total)));
				}
			}
		}
	}

	{
		auto refKeys = getVectorOfMapKeys(refStrandCounts);
		njh::sort(refKeys);
		for (const auto & refKey : refKeys) {
			const auto & ref = refPosStrandCounts[refKey];
			auto posKeys = getVectorOfMapKeys(ref);
			njh::sort(posKeys);
			for (const auto & posKey : posKeys) {
				const auto & pos = ref.at(posKey);
				uint64_t forward = pos.at(true);
				uint64_t reverse = pos.at(false);
				uint64_t total = forward + reverse;
				posCountsTab.content_.emplace_back(
						toVecStr(refKey, posKey, total, forward, reverse,
								forward / static_cast<double>(total)));
			}
		}
	}

	{
		auto refKeys = getVectorOfMapKeys(refStrandCounts);
		njh::sort(refKeys);
		for (const auto & refKey : refKeys) {
			const auto & ref = refSizeStrandCounts[refKey];
			auto sizeKeys = getVectorOfMapKeys(ref);
			njh::sort(sizeKeys);
			for (const auto & sizeKey : sizeKeys) {
				const auto & size = ref.at(sizeKey);
				uint64_t forward = size.at(true);
				uint64_t reverse = size.at(false);
				uint64_t total = forward + reverse;
				sizeCountsTab.content_.emplace_back(
						toVecStr(refKey, sizeKey, total, forward, reverse,
								forward / static_cast<double>(total)));
			}
		}
	}

	{
		auto refKeys = getVectorOfMapKeys(refStrandCounts);
		njh::sort(refKeys);
		for (const auto & refKey : refKeys) {
			const auto & ref = refStrandCounts[refKey];
			uint64_t forward = ref.at(true);
			uint64_t reverse = ref.at(false);
			uint64_t total = forward + reverse;
			refCountsTab.content_.emplace_back(
					toVecStr(refKey, total, forward, reverse,
							forward / static_cast<double>(total)));
		}
	}

	{
		auto refKeys = getVectorOfMapKeys(queryCoverageCount);
		njh::sort(refKeys);
		for (const auto & refKey : refKeys) {
			const auto & ref = queryCoverageCount[refKey];
			auto covKeys = getVectorOfMapKeys(ref);
			njh::sort(covKeys);
			for (const auto & covKey : covKeys) {

				const auto & size = ref.at(covKey);
				uint64_t forward = size.at(true);
				uint64_t reverse = size.at(false);
				uint64_t total = forward + reverse;
				refQueryCoverageTab.content_.emplace_back(
						toVecStr(refKey, covKey/1000.0, total, forward, reverse,
								forward / static_cast<double>(total)));
			}
		}
	}


	sizePosCountsTab.outPutContents(
			TableIOOpts(OutOptions(setUp.pars_.directoryName_ + "sizePosCounts", ".tab.txt")));
	posCountsTab.outPutContents(
			TableIOOpts(OutOptions(setUp.pars_.directoryName_ + "posCounts", ".tab.txt")));
	sizeCountsTab.outPutContents(
			TableIOOpts(OutOptions(setUp.pars_.directoryName_ + "sizeCounts", ".tab.txt")));
	refCountsTab.outPutContents(
			TableIOOpts(OutOptions(setUp.pars_.directoryName_ + "refCounts", ".tab.txt")));
	refQueryCoverageTab.outPutContents(
			TableIOOpts(OutOptions(setUp.pars_.directoryName_ + "qCoveageCounts", ".tab.txt")));

	std::ofstream filterInfoFile;
	openTextFile(filterInfoFile,
			OutOptions(bfs::path(setUp.pars_.directoryName_ + "readFilterInfo.tab.txt")));
	filterInfoFile << "totalReads\treadsCounted\tMappingQuality<" << mappingQual
			<< "\treadsNotMapped" << std::endl;
	filterInfoFile << totalReads << "\t"
			<< getPercentageString(readCount, totalReads) << "\t"
			<< getPercentageString(readsBellowMappingQuality, totalReads) << "\t"
			<< getPercentageString(readsNotMapped, totalReads) << std::endl;

	return 0;
}







void outputAllCountsInfo(RefCounter & counter, std::string dirName,
		uint32_t depthCutOff, double fraqCutOff, std::vector<char> baseAlphabet) {
	njh::appendAsNeeded(dirName, "/");
	std::ofstream outLowFreqFile;
	std::ofstream outAllFile;
	openTextFile(outAllFile,
			OutOptions(bfs::path(dirName + "baseCountsAll.tab.csv")));
	openTextFile(outLowFreqFile,
			OutOptions(bfs::path(dirName + "baseCounts_lowFreq.tab.csv")));
	counter.printInfo(outAllFile);
	counter.printInfo(outLowFreqFile, depthCutOff, fraqCutOff);
	std::ofstream outLowFreqByQualFile;
	std::ofstream outAllByQualFile;
	openTextFile(outAllByQualFile,
			OutOptions(bfs::path(dirName + "baseCountsAllByQual.tab.csv")));
	openTextFile(outLowFreqByQualFile,
			OutOptions(bfs::path(dirName + "baseCounts_lowFreqByQual.tab.csv")));
	counter.printInfoStratifiedByQual(outAllByQualFile);
	counter.printInfoStratifiedByQual(outLowFreqByQualFile, depthCutOff,
			fraqCutOff);

	{
		auto baseChanges = counter.getRatesAgainstRef();
		auto rates = baseChanges.getRates(baseAlphabet);
		table ratesTab(VecStr { "refBase", "seqBase", "freq", "rate" });
		for (auto refBase : baseAlphabet) {
			for (auto seqBase : baseAlphabet) {
				ratesTab.content_.emplace_back(
						toVecStr(refBase, seqBase, baseChanges.mat_[refBase][seqBase],
								rates[refBase][seqBase]));
			}
		}
		ratesTab.outPutContents(
				TableIOOpts(OutOptions(njh::files::join(dirName, "baseRates.tab.csv")),
						"\t", ratesTab.hasHeader_));
	}
	{
		auto baseChanges = counter.getRatesAgainstRefHq();
		auto rates = baseChanges.getRates(baseAlphabet);
		table ratesTab(VecStr { "refBase", "seqBase", "freq", "rate" });
		for (auto refBase : baseAlphabet) {
			for (auto seqBase : baseAlphabet) {
				ratesTab.content_.emplace_back(
						toVecStr(refBase, seqBase, baseChanges.mat_[refBase][seqBase],
								rates[refBase][seqBase]));
			}
		}
		ratesTab.outPutContents(
				TableIOOpts(OutOptions(njh::files::join(dirName, "baseRatesHq.tab.csv")),
						"\t", ratesTab.hasHeader_));
	}

	bfs::path indelInfoDir = njh::files::makeDir(dirName, njh::files::MkdirPar("indelSpecificInfo", false));
	std::ofstream outIndelSeqFile;
	std::ofstream outIndelSizeFile;
	openTextFile(outIndelSeqFile,
			OutOptions(bfs::path(indelInfoDir.string() + "indelSeq.tab.csv")));
	openTextFile(outIndelSizeFile,
			OutOptions(bfs::path(indelInfoDir.string() + "indelSize.tab.csv")));
	counter.printIndelSeqInfo(outIndelSeqFile);
	counter.printIndelSizeInfo(outIndelSizeFile);
}

void outputCountsInfo(RefCounter & counter, std::string dirName,
		uint32_t depthCutOff, double fraqCutOff, std::vector<char> baseAlphabet,
		const std::string & refName, size_t start, size_t stop) {
	njh::appendAsNeeded(dirName, "/");
	std::ofstream outLowFreqFile;
	std::ofstream outAllFile;
	openTextFile(outAllFile,
			OutOptions(bfs::path(dirName + "baseCountsAll.tab.csv")));
	openTextFile(outLowFreqFile,
			OutOptions(bfs::path(dirName + "baseCounts_lowFreq.tab.csv")));
	counter.printInfo(outAllFile, refName, start,
			stop);
	counter.printInfo(outLowFreqFile, refName, start,
			stop, depthCutOff, fraqCutOff);
	std::ofstream outLowFreqByQualFile;
	std::ofstream outAllByQualFile;
	openTextFile(outAllByQualFile,
			OutOptions(bfs::path(dirName + "baseCountsAllByQual.tab.csv")));
	openTextFile(outLowFreqByQualFile,
			OutOptions(
					bfs::path(dirName + "baseCounts_lowFreqByQual.tab.csv")));
	counter.printInfoStratifiedByQual(outAllByQualFile, refName,
			start, stop);
	counter.printInfoStratifiedByQual(outLowFreqByQualFile, refName,
			start, stop, depthCutOff, fraqCutOff);

	bfs::path indelInfoDir = njh::files::makeDir(dirName, njh::files::MkdirPar("indelSpecificInfo", false));
	std::ofstream outIndelSeqFile;
	std::ofstream outIndelSizeFile;
	openTextFile(outIndelSeqFile,
			OutOptions(bfs::path(indelInfoDir.string() + "indelSeq.tab.csv")));
	openTextFile(outIndelSizeFile,
			OutOptions(bfs::path(indelInfoDir.string() + "indelSize.tab.csv")));
	counter.printIndelSeqInfo( outIndelSeqFile,refName, start, stop);
	counter.printIndelSizeInfo(outIndelSizeFile,refName, start, stop);

	{
		auto baseChanges = counter.getRatesAgainstRef(refName, start, stop);
		auto rates = baseChanges.getRates(baseAlphabet);
		table ratesTab(VecStr { "refBase", "seqBase", "freq", "rate" });
		for (auto refBase : baseAlphabet) {
			for (auto seqBase : baseAlphabet) {
				ratesTab.content_.emplace_back(
						toVecStr(refBase, seqBase, baseChanges.mat_[refBase][seqBase],
								rates[refBase][seqBase]));
			}
		}
		ratesTab.outPutContents(
						TableIOOpts(OutOptions(njh::files::join(dirName, "baseRates.tab.csv")),
								"\t", ratesTab.hasHeader_));
	}
	{
		auto baseChangesHq = counter.getRatesAgainstRefHq(refName, start, stop);
		auto ratesHq = baseChangesHq.getRates(baseAlphabet);
		table ratesTabHq(VecStr { "refBase", "seqBase", "freq", "rate" });
		for (auto refBase : baseAlphabet) {
			for (auto seqBase : baseAlphabet) {
				ratesTabHq.content_.emplace_back(
						toVecStr(refBase, seqBase, baseChangesHq.mat_[refBase][seqBase],
								ratesHq[refBase][seqBase]));
			}
		}
		ratesTabHq.outPutContents(
				TableIOOpts(OutOptions(njh::files::join(dirName, "baseRatesHq.tab.csv")),
						"\t", ratesTabHq.hasHeader_));
	}
}



} /* namespace njhseq */
