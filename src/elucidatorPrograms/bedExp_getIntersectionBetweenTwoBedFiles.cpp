/*
 * bedExp_getIntersectionBetweenTwoBedFiles.cpp
 *
 *  Created on: Oct 20, 2020
 *      Author: nick
 */





#include "bedExp.hpp"

#include "elucidator/BamToolsUtils.h"
#include "elucidator/objects/BioDataObject.h"


namespace njhseq {



int bedExpRunner::getIntersectionBetweenTwoBedFiles(const njh::progutils::CmdArgs & inputCommands) {

	OutOptions outOpts(bfs::path(""), ".bed");
	bfs::path bed1 = "";
	bfs::path bed2 = "";
	uint32_t minOverlap = 1;
	bool doNotSkipSameName = false;
	seqSetUp setUp(inputCommands);
	setUp.description_ = "Get regions that overlap between bed regions and output sub regions that only appear in both";
	setUp.setOption(bed1, "--bed1", "bed1", true);
	setUp.setOption(bed2, "--bed2", "bed2", true);
	setUp.setOption(minOverlap, "--minOverlap", "minOverlap", njh::progutils::ProgramSetUp::CheckCase::NONZERO);
	setUp.setOption(doNotSkipSameName, "--doNotSkipSameName", "do Not Skip Same Name regions, regions with same name are skipped by default so same bed file can be given to get distance from");

	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);



	auto inputRegions = getBeds(bed1);
	auto compareRegions = getBeds(bed2);

	OutputStream out(outOpts);
	BedUtility::coordSort(inputRegions, false);
	BedUtility::coordSort(compareRegions, false);



	std::vector<uint32_t> allDistances;

	for(const auto & inputRegion : inputRegions){
		std::vector<std::shared_ptr<Bed6RecordCore>> overlappingRegions;
		for(const auto & compRegion : compareRegions){
			if(compRegion->chrom_ < inputRegion->chrom_){
				continue;
			}
			if(compRegion->chrom_ > inputRegion->chrom_){
				break;
			}
			if(compRegion->name_ == inputRegion->name_ && !doNotSkipSameName){
				continue;
			}
			if(compRegion->strand_ != inputRegion->strand_){
				continue;
			}
			if(inputRegion->overlaps(*compRegion, minOverlap)){
				overlappingRegions.emplace_back(compRegion);
			}
		}
		if(overlappingRegions.size() > 0){
			for(const auto & overReg : overlappingRegions){
				bool scoreIsLength = inputRegion->score_ == inputRegion->length();
				Bed6RecordCore outRegion = *inputRegion;
				outRegion.chromStart_ = std::max(overReg->chromStart_, inputRegion->chromStart_);
				outRegion.chromEnd_ = std::min(overReg->chromEnd_, inputRegion->chromEnd_);
				outRegion.name_ = outRegion.genUIDFromCoordsWithStrand();
				if(scoreIsLength){
					outRegion.score_ = outRegion.length();
				}
				out << outRegion.toDelimStr() << std::endl;
			}
		}
	}

	return 0;
}





int bedExpRunner::removeSubRegionsFromBedFile(const njh::progutils::CmdArgs & inputCommands) {

	OutOptions outOpts(bfs::path(""), ".bed");
	bfs::path bed1 = "";
	bfs::path bed2 = "";
	uint32_t minOverlap = 1;
	bool doNotSkipSameName = false;
	seqSetUp setUp(inputCommands);
	setUp.description_ = "Intersect Bed regions and remove the regions from the interesting bed";
	setUp.setOption(bed1, "--bed", "bed", true);
	setUp.setOption(bed2, "--removeBed", "regions to remove bed file", true);
	setUp.setOption(minOverlap, "--minOverlap", "minOverlap", njh::progutils::ProgramSetUp::CheckCase::NONZERO);
	setUp.setOption(doNotSkipSameName, "--doNotSkipSameName", "do Not Skip Same Name regions, regions with same name are skipped by default so same bed file can be given to get distance from");

	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);



	auto inputRegions = getBeds(bed1);
	auto compareRegions = getBeds(bed2);

	OutputStream out(outOpts);
	BedUtility::coordSort(inputRegions, false);
	BedUtility::coordSort(compareRegions, false);



	std::vector<uint32_t> allDistances;

	for(const auto & inputRegion : inputRegions){
		std::vector<std::shared_ptr<Bed6RecordCore>> overlappingRegions;
		for(const auto & compRegion : compareRegions){
			if(compRegion->chrom_ < inputRegion->chrom_){
				continue;
			}
			if(compRegion->chrom_ > inputRegion->chrom_){
				break;
			}
			if(compRegion->name_ == inputRegion->name_ && !doNotSkipSameName){
				continue;
			}
			if(compRegion->strand_ != inputRegion->strand_){
				continue;
			}
			if(inputRegion->overlaps(*compRegion, minOverlap)){
				overlappingRegions.emplace_back(compRegion);
			}
		}
		if(overlappingRegions.size() > 0){
			std::map<uint32_t, uint32_t> positionsCovered;
			for(const auto pos : iter::range(inputRegion->length())){
				positionsCovered[pos] = 0;
			}
			for(const auto & overReg : overlappingRegions){
				for(const auto pos : iter::range(overReg->chromStart_, overReg->chromEnd_)){
					positionsCovered[pos] += 1;
				}
			}
			struct StartEnd{
				StartEnd(uint32_t start, uint32_t end):start_(start), end_(end){

				}
				StartEnd(uint32_t start):start_(start), end_(start + 1){

				}
				StartEnd():start_(std::numeric_limits<uint32_t>::max() -1 ), end_(std::numeric_limits<uint32_t>::max()){

				}
				uint32_t start_;
				uint32_t end_;
			};

			std::vector<StartEnd> regionsNotCovered;
			for(const auto & posCov : positionsCovered){
				if(0 == posCov.second ){
					if(regionsNotCovered.empty() || regionsNotCovered.back().end_ != posCov.first){
						regionsNotCovered.emplace_back(posCov.first);
					}else{
						++regionsNotCovered.back().end_;
					}
				}
			}
			for(const auto & regionNotCovered : regionsNotCovered){
				bool scoreIsLength = inputRegion->score_ == inputRegion->length();
				Bed6RecordCore outRegion = *inputRegion;
				outRegion.chromStart_ = regionNotCovered.start_;
				outRegion.chromEnd_ =regionNotCovered.end_;
				outRegion.name_ = outRegion.genUIDFromCoordsWithStrand();
				if(scoreIsLength){
					outRegion.score_ = outRegion.length();
				}
				out << outRegion.toDelimStr() << std::endl;
			}
		}else{
			out << inputRegion->toDelimStr() << std::endl;
		}
	}

	return 0;
}

}  // namespace njhseq

