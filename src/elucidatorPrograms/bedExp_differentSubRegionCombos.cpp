/*
 * bedExp_differentSubRegionCombos.cpp
 *
 *  Created on: Jun 2, 2020
 *      Author: nick
 */



#include "bedExp.hpp"
#include <njhseq/objects/BioDataObject.h>

#include "elucidator/objects/BioDataObject.h"
#include "elucidator/BioRecordsUtils/BedUtility.hpp"
#include <njhseq/objects/counters/DNABaseCounter.hpp>



namespace njhseq {




int bedExpRunner::differentSubRegionCombos(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bedFile;
	bfs::path genomeFnp;
	uint32_t minLen = 100;
	uint32_t maxLen = std::numeric_limits<uint32_t>::max();
	uint32_t subRegionSize = 20;
	bool includeFromFront  = false;
	bool includeFromBack   = false;
	bool justToNextRegion = false;

	OutOptions outOpts;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.setOption(bedFile, "--bed", "Bed file to parse", true);
	setUp.setOption(genomeFnp, "--genomeFnp", "genome fnp", true);

	setUp.setOption(minLen, "--minLen", "minLen");
	setUp.setOption(maxLen, "--maxLen", "maxLen");


	setUp.setOption(subRegionSize,    "--subRegionSize", "Sub Region Size");
	setUp.setOption(includeFromFront, "--includeFromFront", "include From Front");
	setUp.setOption(includeFromBack,  "--includeFromBack", "include From Back");
	setUp.setOption(justToNextRegion,  "--onlyAdjacentRegions", "Instead of making every possible regions, do just between adjacent regions");




	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	OutputStream out(outOpts);
	auto regionsRaw = getBeds(bedFile);

	std::unordered_map<std::string, uint32_t> genomeLen;
	auto genomeOpt = SeqIOOptions::genFastaIn(genomeFnp);
	genomeOpt.includeWhiteSpaceInName_ = false;
	{
		seqInfo seq;
		SeqInput reader(genomeOpt);
		reader.openIn();
		while(reader.readNextRead(seq)){
			genomeLen[seq.name_] = len(seq);
		}
	}


	BedUtility::coordSort(regionsRaw, false);
	std::vector<std::shared_ptr<Bed6RecordCore>> regions;
	regions.emplace_back(std::make_shared<Bed6RecordCore>(*regionsRaw.front()));

	for (const auto regPos : iter::range<uint32_t>(1, regionsRaw.size())) {
		if (regions.back()->overlaps(*regionsRaw[regPos], 1) || (regions.back()->chrom_ == regionsRaw[regPos]->chrom_ && regions.back()->chromEnd_ == regionsRaw[regPos]->chromStart_) ) {
			regions.back()->chromEnd_ = std::max(
					regions.back()->chromEnd_,
					regionsRaw[regPos]->chromEnd_);
		} else {
			regions.emplace_back(
					std::make_shared<Bed6RecordCore>(*regionsRaw[regPos]));
		}
	}


	std::unordered_map<std::string, std::vector<std::shared_ptr<Bed6RecordCore>>> regionsByChrom;
	for(const auto & region : regions){
		regionsByChrom[region->chrom_].emplace_back(region);
	}

	if(includeFromFront){
		for(auto & byChrom : regionsByChrom){
			if(byChrom.second.front()->chromStart_ > 10){
				byChrom.second.emplace_back(std::make_shared<Bed6RecordCore>(byChrom.first, 0, 1, "start", 1, '+'));
				BedUtility::coordSort(byChrom.second, false);
			}
		}
	}
	if(includeFromBack){
		for(auto & byChrom : regionsByChrom){
			if(!njh::in(byChrom.first, genomeLen)){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "don't have chromosome " << byChrom.first << " in genome file " << genomeOpt.firstName_<< "\n";
				ss << genomeOpt.firstName_ << "has: " << njh::conToStr(njh::getVecOfMapKeys(genomeLen), ",") << "\n";
				throw std::runtime_error{ss.str()};
			}
			if(genomeLen[byChrom.first] - byChrom.second.back()->chromEnd_ > 10){
				byChrom.second.emplace_back(std::make_shared<Bed6RecordCore>(byChrom.first, genomeLen[byChrom.first] -1, genomeLen[byChrom.first], "back", 1, '+'));
			}
		}
	}


	for(auto & byChrom : regionsByChrom){
		if(byChrom.second.size() > 1){
			for(const auto  pos : iter::range(byChrom.second.size() -1 )){
				for(const auto  downStreamPos : iter::range(pos + 1, byChrom.second.size())){
					Bed6RecordCore startReg = *byChrom.second[pos];
					Bed6RecordCore endReg   = *byChrom.second[downStreamPos];
					if(startReg.length() > subRegionSize){
						startReg.chromStart_ = startReg.chromEnd_ - subRegionSize;
					}
					if(endReg.length() > subRegionSize){
						endReg.chromEnd_ = endReg.chromStart_ + subRegionSize;
					}
					Bed6RecordCore outRegion = startReg;
					outRegion.chromEnd_ = endReg.chromEnd_;
					outRegion.name_ = njh::pasteAsStr(startReg.name_, "--", endReg.name_);
					outRegion.score_ = outRegion.length();
					if(outRegion.length() >= minLen && outRegion.length() <= maxLen){
						out << outRegion.toDelimStr() << std::endl;
					}
					if(justToNextRegion && outRegion.length() >= minLen && outRegion.length() <= maxLen){
						break;
					}
				}
			}
		}
	}



	return 0;
}





}  // namespace njhseq

