//
// Created by Nicholas Hathaway on 12/3/23.
//



#include "bamExp.hpp"
#include <TwoBit.h>

#include "elucidator/BamToolsUtils.h"
#include <njhseq/objects/BioDataObject.h>


namespace njhseq {

int bamExpRunner::mergeImproperMatePositions(
		const njh::progutils::CmdArgs & inputCommands) {
	OutOptions outOpts(bfs::path(""), ".bed");
	bfs::path bedFnp;
	uint32_t minReadAmount = 10;
	uint32_t minInsertSize = 5000;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.setOption(bedFnp, "--bedFnp", "Bed file to get the merged improper mates", true);
	setUp.setOption(minReadAmount, "--minReadAmount", "min Read Amount");
	setUp.setOption(minInsertSize, "--minInsertSize", "minimum Insert Size to include in merge");


	setUp.processReadInNames( { "--bam" }, true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	BamTools::BamReader bReader;
	bReader.Open(setUp.pars_.ioOptions_.firstName_.string());
	checkBamOpenThrow(bReader, setUp.pars_.ioOptions_.firstName_.string());
	loadBamIndexThrow(bReader);
	auto refData = bReader.GetReferenceData();

	BamTools::BamAlignment bAln;

	auto regions = bedPtrsToGenomicRegs(getBeds(bedFnp));
	OutputStream out(outOpts);
	for (const auto& reg: regions) {
		setBamFileRegionThrow(bReader, reg);
		std::vector<Bed6RecordCore> otherRegions;

		while (bReader.GetNextAlignment(bAln)) {
			if (bAln.IsMapped() && bAln.IsMateMapped() && !bAln.IsProperPair() && std::abs(bAln.InsertSize) >=minInsertSize) {
				// std::cout << bAln.InsertSize << std::endl;
				otherRegions.emplace_back(
				                          refData[bAln.MateRefID].RefName,
				                          bAln.MatePosition,
				                          bAln.MatePosition + bAln.QueryBases.length(),
				                          bAln.Name,
				                          1,
				                          bAln.IsMateReverseStrand() ? '-': '+');
			}
		}
		BedUtility::coordSort(otherRegions);
		//merge regions
		if(otherRegions.size() > 1) {
			std::vector<Bed6RecordCore> mergedOtherRegions;
			mergedOtherRegions.emplace_back(otherRegions[0]);
			for(const auto & regPos : iter::range(1UL, otherRegions.size())) {
				if(mergedOtherRegions.back().overlaps(otherRegions[regPos],1)) {
					mergedOtherRegions.back().chromEnd_ = std::max(mergedOtherRegions.back().chromEnd_, otherRegions[regPos].chromEnd_);
					mergedOtherRegions.back().score_ += otherRegions[regPos].score_;
				} else {
					mergedOtherRegions.emplace_back(otherRegions[regPos]);
				}
			}
			for(const auto & merged : mergedOtherRegions) {
				if(merged.score_ >= minReadAmount) {
					auto outRegion = merged;
					outRegion.name_ = reg.uid_ + "__" + merged.genUIDFromCoords();
					out << outRegion.toDelimStrWithExtra() << "\t" << outRegion.length()  << "\t" << reg.uid_ << std::endl;
				}
			}
		}
	}


	return 0;
}


} //namespace njhseq


