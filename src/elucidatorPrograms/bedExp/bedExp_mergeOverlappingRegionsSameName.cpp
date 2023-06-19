//
// Created by Nicholas Hathaway on 8/16/22.
//
#include "bedExp.hpp"

#include "elucidator/BamToolsUtils.h"
#include "elucidator/objects/BioDataObject.h"


namespace njhseq {



int bedExpRunner::mergeOverlappingRegionsSameName(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bedFnp;
	OutOptions outOpts(bfs::path(""), ".bed");
	seqSetUp setUp(inputCommands);
	setUp.description_ = "Some base coverage smoothing for output by bamMulticovBases";
	setUp.setOption(bedFnp, "--bedFnp", "bed file to merge", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	auto regions = getBeds(bedFnp);
	OutputStream out(outOpts);
	auto bedCoordSorterFunc =
					[](const std::shared_ptr<Bed6RecordCore> & reg1In, const  std::shared_ptr<Bed6RecordCore> & reg2In) {
						const auto & reg1 = getRef(reg1In);
						const auto & reg2 = getRef(reg2In);
						if(reg1.name_ == reg2.name_){
							if(reg1.chrom_ == reg2.chrom_) {
								if(reg1.chromStart_ == reg2.chromStart_) {
									return reg1.chromEnd_ < reg2.chromEnd_;
								} else {
									return reg1.chromStart_ < reg2.chromStart_;
								}
							} else {
								return reg1.chrom_ < reg2.chrom_;
							}
						}else{
							return reg1.name_ < reg2.name_;
						}
					};

	njh::sort(regions, bedCoordSorterFunc);

	std::vector<Bed6RecordCore> outRegions;
	for (const auto &region: regions) {
		if (outRegions.empty()) {
			outRegions.emplace_back(*region);
		} else {
			if (outRegions.back().name_ == region->name_ &&
					(region->overlaps(outRegions.back(), 1) ||
					 region->chromStart_ == outRegions.back().chromEnd_)) {
				bool scoreIsLen = outRegions.back().length() == outRegions.back().score_;
				outRegions.back().chromEnd_ = std::max(outRegions.back().chromEnd_, region->chromEnd_);
				outRegions.back().chromStart_ = std::max(outRegions.back().chromStart_, region->chromStart_);
				if (scoreIsLen) {
					outRegions.back().score_ = outRegions.back().length();
				}
			} else {
				outRegions.emplace_back(*region);
			}
		}
	}
	BedUtility::coordSort(outRegions);
	for(const auto & reg : outRegions){
		out << reg.toDelimStrWithExtra() << std::endl;
	}
	return 0;
}


}  // namespace njhseq

