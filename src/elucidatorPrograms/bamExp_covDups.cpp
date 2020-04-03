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



namespace njhseq {





int bamExpRunner::bamMulticovBases(const njh::progutils::CmdArgs & inputCommands){
	bool writeOutTemporaryCoverageFiles = false;
	bool keepTempCovFiles = false;
	OutOptions outOpts(bfs::path(""), ".tab.txt");
	uint32_t numThreads = 1;
	bfs::path bedFnp = "";
	std::string bams = "";
	std::string pat = ".*.bam$";
	bool noHeader = false;
	bool countDups = false;
	uint32_t mapQualityCutOff = 20;
	bool dontHandlePairs = false;
	bfs::path directory = "./";
	seqSetUp setUp(inputCommands);
	setUp.description_ = "Get the coverage in base count for bam files for certain regions";
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(writeOutTemporaryCoverageFiles, "--writeOutTemporaryCoverageFiles", "Write Out Temporary Coverage Files");
	setUp.setOption(directory, "--directory", "Directory to search for bam files");
	setUp.setOption(dontHandlePairs, "--dontHandlePairs",   "Don't Handle Paired reads, this means pairs covering the same region will count twice");
	setUp.setOption(countDups, "--countDups",   "Count records marked duplicate");
	setUp.setOption(mapQualityCutOff, "--mapQualityCutOff",   "Only reads that are this mapping quality and above (inclusive)");
	setUp.setOption(keepTempCovFiles, "--keepTempCovFiles", "Temporary Coverage Files");
	setUp.setOption(numThreads, "--numThreads", "Number of threads to use");
	setUp.processWritingOptions(outOpts);
	setUp.setOption(bedFnp, "--bedFnp", "Bed file of regions to get coverage for", true, "Input");
	setUp.setOption(pat, "--pat", "Pattern in current directory to get coverage for", false, "Input");
	setUp.setOption(bams, "--bams", "Either a file with the name of a bam file on each line or a comma separated value of bam file paths", false, "Input");
	setUp.setOption(noHeader, "--noHeader", "Don't output a header so output can be treated like a bed file", false, "Writing Output");
	setUp.finishSetUp(std::cout);

	njh::stopWatch watch;
	watch.setLapName("Initial Set Up");

	OutputStream outFile(outOpts);

	auto bamFnps = njh::files::gatherFilesByPatOrNames(directory, std::regex{pat}, bams);
	if(bamFnps.empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "no bam files found with pattern: " << pat << " in " << directory << "\n";
		throw std::runtime_error{ss.str()};
	}
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
	std::mutex pairsMut;
	watch.startNewLap("Getting coverage");
	std::vector<bfs::path> tempCovFiles;

	for(const auto & bamFnp : bamFnps){
		njhseq::concurrent::BamReaderPool bamPool(bamFnp, numThreads);
		bamPool.openBamFile();
		njh::concurrent::LockableQueue<GenomicRegion> regionsQueue(regions);
		uint32_t pairsStart = pairs.size();
		std::function<void()> getCov = [&bamFnp,&bamPool,&pairs, &pairsMut,&regionsQueue,&countDups,
																		&mapQualityCutOff,&dontHandlePairs](){
			std::vector<std::shared_ptr<BamFnpRegionPair>> currentBamRegionsPairs;
			GenomicRegion region;
			auto bReader = bamPool.popReader();
			BamTools::BamAlignment bAln;
			auto refData = bReader->GetReferenceData();

			while(regionsQueue.getVal(region)){
				auto val = std::make_shared<BamFnpRegionPair>(bamFnp,region);
				setBamFileRegionThrow(*bReader, region);
				BamAlnsCache bCache;
				while(bReader->GetNextAlignmentCore(bAln)){
					if(bAln.IsMapped() && bAln.IsPrimaryAlignment()){
						if(bAln.IsDuplicate() && !countDups){
							continue;
						}
						if(bAln.MapQuality <  mapQualityCutOff){
							continue;
						}
						//only try to find and adjust for the mate's overlap if is mapped and could possibly fall within the region
						if(!dontHandlePairs && bAln.IsPaired() && bAln.IsMateMapped() && bAln.MatePosition < val->region_.end_){
							bAln.BuildCharData();
							if(bCache.has(bAln.Name)){
								//get and adjust current read region
								GenomicRegion alnRegion(bAln, refData);
								alnRegion.start_ = std::max(alnRegion.start_, val->region_.start_);
								alnRegion.end_ = std::min(alnRegion.end_, val->region_.end_);
								//get and adjust the mate's read region
								auto mate = bCache.get(bAln.Name);
								GenomicRegion mateRegion(*mate, refData);
								mateRegion.start_ = std::max(mateRegion.start_, val->region_.start_);
								mateRegion.end_ = std::min(mateRegion.end_, val->region_.end_);
								//now that the two regions have been adjusted to be just what overlaps the current region
								//coverage will be the length of the two regions minus the overlap between the two
								val->coverage_ +=alnRegion.getLen() + mateRegion.getLen() - alnRegion.getOverlapLen(mateRegion);
								bCache.remove(bAln.Name);
							}else{
								bCache.add(bAln);
							}
						}else{
							/**@todo this doesn't take into account gaps, so base coverage isn't precise right here, more of an approximation, should improve */
							GenomicRegion alnRegion(bAln, refData);
							val->coverage_ += val->region_.getOverlapLen(alnRegion);
						}
					}
				}
				//save alignments where mate couldn't be found, this could be for several reasons like mate was dup or mate's mapping quality or other filtering reasons
				for(const auto & name : bCache.getNames()){
					auto aln = bCache.get(name);
					GenomicRegion alnRegion(*aln, refData);
					val->coverage_ += val->region_.getOverlapLen(alnRegion);
				}
				currentBamRegionsPairs.emplace_back(val);
			}
			{
				std::lock_guard<std::mutex> lock(pairsMut);
				addOtherVec(pairs, currentBamRegionsPairs);
			}
		};
		//get coverage for this bam file threaded over the regions
		njh::concurrent::runVoidFunctionThreaded(getCov, numThreads);

		if(writeOutTemporaryCoverageFiles){
			auto sampOutName = njh::files::prependFileBasename(outOpts.outName(), bamFnp.filename().string() + "_");
			tempCovFiles.emplace_back(sampOutName);
			OutOptions sampOutOpts(sampOutName);
			sampOutOpts.transferOverwriteOpts(outOpts);
			OutputStream sampOut(sampOutOpts);
			VecStr header = {"chrom", "start", "end", "name", "score", "strand", bamFnp.filename().string()};
			sampOut << njh::conToStr(header, "\t") << std::endl;
			for(uint32_t pos = pairsStart; pos < pairs.size(); ++pos){
				sampOut << pairs[pos]->region_.genBedRecordCore().toDelimStr() << "\t" << pairs[pos]->coverage_ << "\n";
			}
		}
	}

	watch.startNewLap("Filling table");
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
	njh::concurrent::LockableVec<std::shared_ptr<BamFnpRegionPair>> pairsList(pairs);

	//pairsList.reset();
	std::function<void()> fillTable = [&output, &pairsList, &regUidToRowPos, &bamFnpToCol](){
		std::shared_ptr<BamFnpRegionPair> val;
		while(pairsList.getVal(val)){
			output.content_[regUidToRowPos[val->regUid_]][bamFnpToCol[val->bamFname_]] = estd::to_string(val->coverage_);
		}
	};

	njh::concurrent::runVoidFunctionThreaded(fillTable, numThreads);
	if(noHeader){
		output.hasHeader_ = false;
	}

	watch.startNewLap("Writing table");
	output.outPutContents(outFile, "\t");

	if(!keepTempCovFiles){
		for(const auto & tFile : tempCovFiles){
			if(bfs::exists(tFile)){
				bfs::remove(tFile);
			}
		}
	}

	if(setUp.pars_.debug_){
		watch.logLapTimes(std::cout, true, 6, true);
	}
	return 0;
}





int bamExpRunner::bamDupCounts(const njh::progutils::CmdArgs & inputCommands){
	OutOptions outStubOpts(bfs::path(""), ".tab.txt");
	uint32_t numThreads = 1;
	bfs::path bedFnp = "";
	std::string bams = "";
	std::string pat = ".*.bam$";
	bool noHeader = false;
	seqSetUp setUp(inputCommands);
	setUp.description_ = "Get the count of duplicates for bam files for certain regions";
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(numThreads, "--numThreads", "Number of threads to use");
	setUp.processWritingOptions(outStubOpts);
	setUp.setOption(bedFnp, "--bedFnp", "Bed file of regions to get coverage for", true, "Input");
	setUp.setOption(pat, "--pat", "Pattern in current directory to get coverage for", false, "Input");
	setUp.setOption(bams, "--bams", "Either a file with the name of a bam file on each line or a comma separated value of bam file paths", false, "Input");
	setUp.setOption(noHeader, "--noHeader", "Don't output a header so output can be treated like a bed file", false, "Writing Output");
	setUp.finishSetUp(std::cout);



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

	if("" == outStubOpts.outFilename_){
		outStubOpts.outFilename_ = bedFnp.filename().replace_extension("");
	}
	auto outOptsDubs = outStubOpts;
	outOptsDubs.outFilename_ = outOptsDubs.outFilename_.string() + "_duplicateCounts";
	auto outOptsTotal = outStubOpts;
	outOptsTotal.outFilename_ = outOptsDubs.outFilename_.string() + "_totalCounts";
	OutputStream outFileDubCounts(outOptsDubs);
	OutputStream outFileTotalCounts(outOptsTotal);

	struct BamFnpRegionPair {
		BamFnpRegionPair(const bfs::path & bamFnp, const GenomicRegion & region) :
				bamFnp_(bamFnp), region_(region) {
			regUid_ = region_.createUidFromCoords();
			bamFname_ = bamFnp_.filename().string();
		}
		bfs::path bamFnp_;
		GenomicRegion region_;
		uint32_t totalDubs_ = 0;
		uint32_t totalReads_ = 0;

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

	std::function<void()> getCov = [&pairsList](){
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
					++val->totalReads_;
					if(bAln.IsDuplicate()){
						++val->totalDubs_;
					}
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
	table outputDubs(header);
	outputDubs.content_ = std::vector<std::vector<std::string>>{regions.size(), std::vector<std::string>{header.size()}};
	table outputTotal(header);
	outputTotal.content_ = std::vector<std::vector<std::string>>{regions.size(), std::vector<std::string>{header.size()}};
	uint32_t regRowCount = 0;
	std::unordered_map<std::string, uint32_t> regUidToRowPos;
	for(const auto & reg : regions){
		auto regAsBed = reg.genBedRecordCore();
		auto toks = tokenizeString(regAsBed.toDelimStr(), "\t");
		for(const auto & col : iter::range(toks.size())){
			outputDubs.content_[regRowCount][col] = toks[col];
			outputTotal.content_[regRowCount][col] = toks[col];
		}
		regUidToRowPos[reg.createUidFromCoords()] = regRowCount;
		++regRowCount;
	}
	pairsList.reset();
	std::function<void()> fillTable = [&outputDubs,&outputTotal, &pairsList, &regUidToRowPos, &bamFnpToCol](){
		std::shared_ptr<BamFnpRegionPair> val;
		while(pairsList.getVal(val)){
			outputDubs.content_[regUidToRowPos[val->regUid_]][bamFnpToCol[val->bamFname_]] = estd::to_string(val->totalDubs_);
			outputTotal.content_[regUidToRowPos[val->regUid_]][bamFnpToCol[val->bamFname_]] = estd::to_string(val->totalReads_);
		}
	};
	njh::concurrent::runVoidFunctionThreaded(fillTable, numThreads);
	if(noHeader){
		outputDubs.hasHeader_ = false;
		outputTotal.hasHeader_ = false;

	}
	outputDubs.outPutContents(outFileDubCounts, "\t");
	outputTotal.outPutContents(outFileTotalCounts, "\t");

	return 0;
}



} // namespace njhseq

