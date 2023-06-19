/*
 * bedExp_getAverageDistanceToOtherRegions.cpp
 *
 *  Created on: Jul 13, 2020
 *      Author: nick
 */

#include "bedExp.hpp"

#include "elucidator/BamToolsUtils.h"
#include "elucidator/objects/BioDataObject.h"
#include "elucidator/BioRecordsUtils/BedUtility.hpp"

namespace njhseq {



int bedExpRunner::getAverageDistanceToOtherRegions(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bedFile = "";
	bfs::path otherRegions = "";
	uint32_t buffer = 0;
	bool doNotSkipSameName = false;
	OutOptions outOpts(bfs::path(""), ".tab.txt");
	uint32_t centimorgans = 0;

	seqSetUp setUp(inputCommands);
	setUp.setOption(bedFile, "--bed", "Bed file", true);
	setUp.setOption(otherRegions, "--compareBed", "regions to get distances from", true);
	setUp.setOption(centimorgans, "--centimorgans", "centimorgans", true);

	setUp.setOption(buffer, "--buffer", "Skip regions that are this distance or below from the calculations to avoid calculating from related targets");


	setUp.setOption(doNotSkipSameName, "--doNotSkipSameName", "do Not Skip Same Name regions, regions with same name are skipped by default so same bed file can be given to get distance from");

	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);


	auto inputRegions = getBeds(bedFile);
	auto compareRegions = getBeds(otherRegions);

	OutputStream out(outOpts);
	out << "name\tmeanBPDist\tmedianBPDist\tmeanCentimorgans\tmedianCentimorgans" << std::endl;
	BedUtility::coordSort(inputRegions, false);
	BedUtility::coordSort(compareRegions, false);



	std::vector<uint32_t> allDistances;

	for(const auto & inputRegion : inputRegions){
		std::vector<uint32_t> distsForInputRegion;
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
			auto dist = inputRegion->getDistanceBetween(*compRegion);
			if(dist >= buffer){
				distsForInputRegion.emplace_back(dist);
				allDistances.emplace_back(dist);
			}
		}
		if(distsForInputRegion.size() > 0){
			std::vector<double> centimorgansDists;
			for(const auto & dist : distsForInputRegion){
				centimorgansDists.emplace_back(dist/static_cast<double>(centimorgans));
			}
			out << inputRegion->name_
					<< "\t" << vectorMean(distsForInputRegion)
					<< "\t" << vectorMedianRef(distsForInputRegion)
					<< "\t" << vectorMean(centimorgansDists)
					<< "\t" << vectorMedianRef(centimorgansDists) << std::endl;
		}else{
			out << inputRegion->name_
					<< "\t" << "NA"
					<< "\t" << "NA"
					<< "\t" << "NA"
					<< "\t" << "NA"  << std::endl;
		}
	}

	if(allDistances.size() > 0){
		std::vector<double> centimorgansDists;
		for(const auto & dist : allDistances){
			centimorgansDists.emplace_back(dist/static_cast<double>(centimorgans));
		}
		out << "all"
				<< "\t" << vectorMean(allDistances)
				<< "\t" << vectorMedianRef(allDistances)
				<< "\t" << vectorMean(centimorgansDists)
				<< "\t" << vectorMedianRef(centimorgansDists) << std::endl;
	}else{
		out << "all"
				<< "\t" << "NA"
				<< "\t" << "NA"
				<< "\t" << "NA"
				<< "\t" << "NA"  << std::endl;
	}

	return 0;
}






int bedExpRunner::getDistanceToClostestRegion(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bedFile = "";
	bfs::path otherRegions = "";
	uint32_t buffer = 0;
	bool doNotSkipSameName = false;
	OutOptions outOpts(bfs::path(""), ".tab.txt");
	uint32_t centimorgans = 0;

	seqSetUp setUp(inputCommands);
	setUp.setOption(bedFile, "--bed", "Bed file", true);
	setUp.setOption(otherRegions, "--compareBed", "regions to get distances from", true);
	setUp.setOption(centimorgans, "--centimorgans", "centimorgans", true);

	setUp.setOption(buffer, "--buffer", "Skip regions that are this distance or below from the calculations to avoid calculating from related targets");


	setUp.setOption(doNotSkipSameName, "--doNotSkipSameName", "do Not Skip Same Name regions, regions with same name are skipped by default so same bed file can be given to get distance from");

	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);


	auto inputRegions = getBeds(bedFile);
	auto compareRegions = getBeds(otherRegions);

	OutputStream out(outOpts);
	out << "name\tclosestRegion\tdistanceInBases\tdistanceInCentimorgans" << std::endl;
	BedUtility::coordSort(inputRegions, false);
	BedUtility::coordSort(compareRegions, false);

	for(const auto & inputRegion : inputRegions){
		uint32_t smallestDistance = std::numeric_limits<uint32_t>::max();
		std::shared_ptr<Bed6RecordCore> closestRegion;
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

			auto dist = inputRegion->getDistanceBetween(*compRegion);
			if(dist >= buffer){
				if(dist < smallestDistance){
					smallestDistance = dist;
					closestRegion = compRegion;
				}
			}
		}
		if(closestRegion){
			out << inputRegion->name_
					<< "\t" << closestRegion->name_
					<< "\t" << smallestDistance
					<< "\t" << smallestDistance/static_cast<double>(centimorgans)
					<< std::endl;
		}else{
			out << inputRegion->name_
					<< "\t" << "NA"
					<< "\t" << "NA"
					<< "\t" << "NA" << std::endl;
		}
	}



	return 0;
}


}  // namespace njhseq



