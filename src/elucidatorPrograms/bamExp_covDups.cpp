/*
 * bamExp_covDups.cpp
 *
 *  Created on: May 26, 2019
 *      Author: nicholashathaway
 */




#include "bamExp.hpp"
#include <TwoBit.h>

#include "elucidator/BamToolsUtils.h"
#include "elucidator/objects/BioDataObject.h"
#include "elucidator/objects/dataContainers/graphs.h"
#include <njhseq/GenomeUtils.h>

#include "elucidator/concurrency/LockableJsonLog.hpp"


namespace njhseq {



int bamExpRunner::bamMulticovBases(const njh::progutils::CmdArgs & inputCommands){
	OutOptions outOpts(bfs::path(""));
	uint32_t numThreads = 1;
	bfs::path bedFnp = "";
	std::string bams = "";
	std::string pat = ".*.bam$";
	bool noHeader = false;
	outOpts.outExtention_ = ".tab.txt";
	seqSetUp setUp(inputCommands);
	setUp.description_ = "Get the coverage in base count for bam files for certain regions";
	setUp.processVerbose();
	setUp.processDebug();
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

	auto getCov = [&pairsList](){
		std::shared_ptr<BamFnpRegionPair> val;
		while(pairsList.getVal(val)){
			BamTools::BamReader bReader;

			bReader.Open(val->bamFnp_.string());
			bReader.LocateIndex();
			setBamFileRegionThrow(bReader, val->region_);
			auto refData = bReader.GetReferenceData();
			BamTools::BamAlignment bAln;
			while(bReader.GetNextAlignmentCore(bAln)){
				if(bAln.IsMapped()){
					/**@todo this doesn't take into account gaps, so base coverage isn't precise right here, more of an approximation, should improve */
					GenomicRegion alnRegion(bAln, refData);
					val->coverage_ += val->region_.getOverlapLen(alnRegion);
				}
			}
		}
	};

	{
		std::vector<std::thread> threads;
		for(uint32_t t = 0; t < numThreads; ++t){
			threads.emplace_back(std::thread(getCov));
		}

		for(auto & t : threads){
			t.join();
		}
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
		for(const auto & col : iter::range(toks.size())){
			output.content_[regRowCount][col] = toks[col];
		}
		regUidToRowPos[reg.createUidFromCoords()] = regRowCount;
		++regRowCount;
	}
	pairsList.reset();
	auto fillTable = [&output, &pairsList, &regUidToRowPos, &bamFnpToCol](){
		std::shared_ptr<BamFnpRegionPair> val;
		while(pairsList.getVal(val)){
			output.content_[regUidToRowPos[val->regUid_]][bamFnpToCol[val->bamFname_]] = estd::to_string(val->coverage_);
		}
	};

	{
		std::vector<std::thread> threads;
		for(uint32_t t = 0; t < numThreads; ++t){
			threads.emplace_back(std::thread(fillTable));
		}

		for(auto & t : threads){
			t.join();
		}
	}
	if(noHeader){
		output.hasHeader_ = false;
	}
	output.outPutContents(outFile, "\t");

	return 0;
}



int bamExpRunner::bamDupCounts(const njh::progutils::CmdArgs & inputCommands){
	OutOptions outOpts(bfs::path(""));
	uint32_t numThreads = 1;
	bfs::path bedFnp = "";
	std::string bams = "";
	std::string pat = ".*.bam$";
	bool noHeader = false;
	outOpts.outExtention_ = ".tab.txt";
	seqSetUp setUp(inputCommands);
	setUp.description_ = "Get the coverage in base count for bam files for certain regions";
	setUp.processVerbose();
	setUp.processDebug();
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

	auto getCov = [&pairsList](){
		std::shared_ptr<BamFnpRegionPair> val;
		while(pairsList.getVal(val)){
			BamTools::BamReader bReader;

			bReader.Open(val->bamFnp_.string());
			bReader.LocateIndex();
			setBamFileRegionThrow(bReader, val->region_);
			auto refData = bReader.GetReferenceData();
			BamTools::BamAlignment bAln;
			while(bReader.GetNextAlignmentCore(bAln)){
				if(bAln.IsMapped()){
					/**@todo this doesn't take into account gaps, so base coverage isn't precise right here, more of an approximation, should improve */
					GenomicRegion alnRegion(bAln, refData);
					val->coverage_ += val->region_.getOverlapLen(alnRegion);
				}
			}
		}
	};

	{
		std::vector<std::thread> threads;
		for(uint32_t t = 0; t < numThreads; ++t){
			threads.emplace_back(std::thread(getCov));
		}

		for(auto & t : threads){
			t.join();
		}
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
		for(const auto & col : iter::range(toks.size())){
			output.content_[regRowCount][col] = toks[col];
		}
		regUidToRowPos[reg.createUidFromCoords()] = regRowCount;
		++regRowCount;
	}
	pairsList.reset();
	auto fillTable = [&output, &pairsList, &regUidToRowPos, &bamFnpToCol](){
		std::shared_ptr<BamFnpRegionPair> val;
		while(pairsList.getVal(val)){
			output.content_[regUidToRowPos[val->regUid_]][bamFnpToCol[val->bamFname_]] = estd::to_string(val->coverage_);
		}
	};

	{
		std::vector<std::thread> threads;
		for(uint32_t t = 0; t < numThreads; ++t){
			threads.emplace_back(std::thread(fillTable));
		}

		for(auto & t : threads){
			t.join();
		}
	}
	if(noHeader){
		output.hasHeader_ = false;
	}
	output.outPutContents(outFile, "\t");

	return 0;
}



} // namespace njhseq

