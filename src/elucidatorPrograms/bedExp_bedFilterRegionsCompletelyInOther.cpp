/*
 * bedExp_bedFilterRegionsCompletelyInOther.cpp
 *
 *  Created on: May 31, 2020
 *      Author: nick
 */


#include "bedExp.hpp"
#include <njhseq/objects/BioDataObject.h>

#include "elucidator/objects/BioDataObject.h"
#include "elucidator/BioRecordsUtils/BedUtility.hpp"
#include "elucidator/objects/counters/DNABaseCounter.hpp"


namespace njhseq {

int bedExpRunner::bedFilterRegionsCompletelyInOther(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bedFile;
	bfs::path intersectWithBed;
	OutOptions outOpts;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.setOption(bedFile, "--bed", "Bed file to parse", true);
	setUp.setOption(intersectWithBed, "--intersectWithBed", "Bed file to intersect with", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	OutputStream out(outOpts);
	auto regions = getBeds(bedFile);

	if(bedFile == intersectWithBed){
		//internal filtration
		auto bedCoordSorterFunc =
				[](const std::shared_ptr<Bed6RecordCore> & reg1In, const std::shared_ptr<Bed6RecordCore> & reg2In) {
					const auto & reg1 = getRef(reg1In);
					const auto & reg2 = getRef(reg2In);
					if(reg1.chrom_ == reg2.chrom_) {
						if(reg1.chromStart_ == reg2.chromStart_) {
							if(reg1.chromEnd_ == reg2.chromEnd_){
								return reg1.name_ < reg2.name_;
							}else{
								return reg1.chromEnd_ > reg2.chromEnd_;
							}
						} else {
							return reg1.chromStart_ < reg2.chromStart_;
						}
					} else {
						return reg1.chrom_ < reg2.chrom_;
					}
		};
		njh::sort(regions, bedCoordSorterFunc);
		std::vector<std::shared_ptr<Bed6RecordCore>> intersectingBeds;
		for(const auto & reg : regions){
			std::vector<std::shared_ptr<Bed6RecordCore>> overlappingRegions;
			for(const auto & inter : intersectingBeds){
				if(inter->chrom_ < reg->chrom_){
					//bother regions are sorted if we haven't reached this region's chromosome yet
					continue;
				}
				if(inter->chrom_ > reg->chrom_){
					//both regions are sorted so if we run into this we can break
					break;
				}
				if(inter->chrom_ == reg->chrom_ && inter->chromStart_ >= reg->chromEnd_){
					//both regions are sorted so if we run into this we can break
					break;
				}
				if(inter->chrom_ == reg->chrom_ && inter->chromEnd_ < reg->chromStart_){
					continue;
				}
				if(reg->chromStart_ >= inter->chromStart_ && reg->chromEnd_ <= inter->chromEnd_){
					overlappingRegions.emplace_back(inter);
				}
			}
			if(overlappingRegions.empty()){
				//add region to intersecting regions to filter further regions;
				intersectingBeds.emplace_back(reg);
				out << reg->toDelimStr() << "\n";
			}
		}
	}else{
		auto intersectingBeds = getBeds(intersectWithBed);
		BedUtility::coordSort(regions, false);
		BedUtility::coordSort(intersectingBeds, false);

		for(const auto & reg : regions){
			std::vector<std::shared_ptr<Bed6RecordCore>> overlappingRegions;
			for(const auto & inter : intersectingBeds){
				if(inter->chrom_ < reg->chrom_){
					//bother regions are sorted if we haven't reached this region's chromosome yet
					continue;
				}
				if(inter->chrom_ > reg->chrom_){
					//both regions are sorted so if we run into this we can break
					break;
				}
				if(inter->chrom_ == reg->chrom_ && inter->chromStart_ >= reg->chromEnd_){
					//both regions are sorted so if we run into this we can break
					break;
				}
				if(inter->chrom_ == reg->chrom_ && inter->chromEnd_ < reg->chromStart_){
					continue;
				}
				if(reg->chromStart_ >= inter->chromStart_ && reg->chromEnd_ <= inter->chromEnd_){
					overlappingRegions.emplace_back(inter);
				}
			}
			if(overlappingRegions.empty()){
				out << reg->toDelimStr() << "\n";
			}
		}
	}

	return 0;
}

}  // namespace njhseq

