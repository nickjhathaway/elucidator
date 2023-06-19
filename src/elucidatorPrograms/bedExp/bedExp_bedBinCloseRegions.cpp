/*
 * bedExp_bedBinCloseRegions.cpp
 *
 *  Created on: Oct 2, 2021
 *      Author: nicholashathaway
 */



#include "bedExp.hpp"
#include <njhseq/objects/BioDataObject.h>

#include "elucidator/objects/BioDataObject.h"
#include "elucidator/BioRecordsUtils/BedUtility.hpp"
#include <njhseq/objects/counters/DNABaseCounter.hpp>




namespace njhseq {




int bedExpRunner::bedBinCloseRegions(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bedFile;
	uint32_t lengthToNextRegion = 10000;
	OutOptions outOpts(bfs::path(""), ".bed");
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.setOption(bedFile, "--bed", "Bed file to parse", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	OutputStream out(outOpts);
	auto regions = getBed3s(bedFile);
	BedUtility::coordSort(regions, false);
	if(regions.size() > 1){
		uint32_t startRegionPos = 0;
		uint32_t currentRegionPos = 0;
		for(const auto pos : iter::range(1UL,regions.size())){
			if(regions[currentRegionPos]->overlaps(*regions[pos], 1) || (regions[currentRegionPos]->chrom_ == regions[pos]->chrom_ && uAbsdiff(regions[currentRegionPos]->chromEnd_, regions[pos]->chromStart_) < lengthToNextRegion)){
				currentRegionPos = pos;
			}else{
				std::string regionName = njh::pasteAsStr(regions[startRegionPos]->genUIDFromCoords(), "--", regions[currentRegionPos]->genUIDFromCoords());
				for(const auto groupPos : iter::range(startRegionPos, currentRegionPos + 1)){
					regions[groupPos]->extraFields_.emplace_back(regionName);
				}
				startRegionPos = pos;
				currentRegionPos = pos;
			}
		}
		std::string regionName = njh::pasteAsStr(regions[startRegionPos]->genUIDFromCoords(), "--", regions[currentRegionPos]->genUIDFromCoords());
		for(const auto groupPos : iter::range(startRegionPos, currentRegionPos + 1)){
			regions[groupPos]->extraFields_.emplace_back(regionName);
		}
	}
	for(const auto & reg : regions){
		out << reg->toDelimStrWithExtra() << std::endl;
	}


	return 0;
}



}  // namespace njhseq
