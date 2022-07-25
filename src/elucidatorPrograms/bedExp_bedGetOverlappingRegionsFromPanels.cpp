/*
 * bedExp_bedGetOverlappingRegionsFromPanels.cpp
 *
 *  Created on: Jul 16, 2021
 *      Author: nick
 */





#include "bedExp.hpp"
#include <njhseq/objects/BioDataObject.h>

#include "elucidator/objects/BioDataObject.h"
#include "elucidator/BioRecordsUtils/BedUtility.hpp"
#include "elucidator/objects/counters/DNABaseCounter.hpp"


namespace njhseq {




int bedExpRunner::bedGetOverlappinBetweenPanels(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path intersectWithBeds;
	bool doNotAdjustForSinglePosIntersection = false;
	OutOptions outOpts(bfs::path(""), ".tab.txt");
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.setOption(doNotAdjustForSinglePosIntersection, "--doNotAdjustForSinglePosIntersection", "do Not Adjust For Single Pos Intersection");
	setUp.setOption(intersectWithBeds, "--intersectWithBeds", "A table with bed files from previous panels to intersect with,two columns 1)panel,2)bed", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	OutputStream out(outOpts);
	// read in beds to intercept with
	table otherBedsTab(intersectWithBeds, "\t", true);
	otherBedsTab.checkForColumnsThrow(VecStr{"panel", "bed"}, __PRETTY_FUNCTION__);
	std::map<std::string, bfs::path> bedsByPanel;
	uint32_t panelPos = otherBedsTab.getColPos("panel");
	uint32_t bedPos = otherBedsTab.getColPos("bed");
	std::unordered_map<std::string, std::string> panelByBeds;
	for(const auto & row : otherBedsTab){
		auto panel = row[panelPos];
		auto bed = row[bedPos];
		if(njh::in(panel, bedsByPanel)){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "already have panel: " << panel<< "\n";
			throw std::runtime_error{ss.str()};
		}
		if(njh::in(bed, panelByBeds)){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "already have bed file: " << bed << " for panel: " << panelByBeds[bed] << "\n";
			throw std::runtime_error{ss.str()};
		}
		bedsByPanel[panel] = bed;
		panelByBeds[bed] = panel;
		if(!bfs::exists(bed)){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << bed << " doesn't exist" << "\n";
			throw std::runtime_error{ss.str()};
		}
	}
	std::vector<std::string> panelNames = getVectorOfMapKeys(bedsByPanel);
	njh::sort(panelNames);

	out << "#chrom\tstart\tend\tname\tlen\tstrand";
	out << "\t" << "panel";
	out << "\t" << njh::conToStr(panelNames, "\t") << std::endl;

	for(const auto & currentPanel : panelNames){

		auto currenntPanelBed = bedsByPanel[currentPanel];
		auto regions = getBeds(currenntPanelBed);

		std::vector<std::map<std::string, double>> coverageForEachPanel;
		std::map<std::string, double> panelCovTemplate;
		for(const auto & panel : panelNames){
			panelCovTemplate[panel] = 0;
		}
		for(uint32_t pos = 0; pos < regions.size(); ++pos){
			coverageForEachPanel.emplace_back(panelCovTemplate);
		}

		for(const auto & bedByPanel : bedsByPanel){
			auto panelRegions = getBeds(bedByPanel.second);

			for(const auto bedPos : iter::range(regions.size())){
				const auto & bed = regions[bedPos];
				std::vector<Bed6RecordCore> overlappingRegions;
				for(const auto & panelBed : panelRegions){
					if(bed->overlaps(*panelBed, 1)){
						overlappingRegions.emplace_back(*panelBed);
					}
				}
				double maxCov = 0;

				for(const auto & overlap : overlappingRegions){
					double cov = bed->getOverlapLen(overlap)/static_cast<double>(bed->length());
					if(!doNotAdjustForSinglePosIntersection){
						if(overlap.length() <=3){
							cov =  bed->getOverlapLen(overlap)/static_cast<double>(overlap.length());
						}
					}
					if(cov > maxCov){
						maxCov = cov;
					}
				}
				coverageForEachPanel[bedPos][bedByPanel.first] = maxCov;
			}
		}

		for(const auto bedPos : iter::range(regions.size())){
			out << regions[bedPos]->toDelimStr();
			out << "\t" << currentPanel;
			for(const auto & panel : panelNames){
				out << "\t" << coverageForEachPanel[bedPos][panel];
			}
			out << std::endl;
		}
	}


	return 0;
}


int bedExpRunner::bedGetOverlappingRegionsFromPanels(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bedFile;
	bfs::path intersectWithBeds;
	bool doNotAdjustForSinglePosIntersection = false;
	OutOptions outOpts(bfs::path(""), ".tab.txt");
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.setOption(doNotAdjustForSinglePosIntersection, "--doNotAdjustForSinglePosIntersection", "do Not Adjust For Single Pos Intersection");

	setUp.setOption(bedFile, "--bed", "Bed file to parse", true);
	setUp.setOption(intersectWithBeds, "--intersectWithBeds", "A table with bed files from previous panels to intersect with,two columns 1)panel,2)bed", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	OutputStream out(outOpts);
	auto regions = getBeds(bedFile);

	uint32_t maxExtraFields = 0;
	for(const auto & region : regions){
		if(region->extraFields_.size() > maxExtraFields){
			maxExtraFields = region->extraFields_.size();
		}
	}
	// read in beds to intercept with
	table otherBedsTab(intersectWithBeds, "\t", true);
	otherBedsTab.checkForColumnsThrow(VecStr{"panel", "bed"}, __PRETTY_FUNCTION__);
	std::map<std::string, bfs::path> bedsByPanel;
	uint32_t panelPos = otherBedsTab.getColPos("panel");
	uint32_t bedPos = otherBedsTab.getColPos("bed");
	std::unordered_map<std::string, std::string> panelByBeds;
	for(const auto & row : otherBedsTab){
		auto panel = row[panelPos];
		auto bed = row[bedPos];
		if(njh::in(panel, bedsByPanel)){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "already have panel: " << panel<< "\n";
			throw std::runtime_error{ss.str()};
		}
		if(njh::in(bed, panelByBeds)){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "already have bed file: " << bed << " for panel: " << panelByBeds[bed] << "\n";
			throw std::runtime_error{ss.str()};
		}
		bedsByPanel[panel] = bed;
		panelByBeds[bed] = panel;
		if(!bfs::exists(bed)){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << bed << " doesn't exist" << "\n";
			throw std::runtime_error{ss.str()};
		}
	}
	std::vector<std::string> panelNames = getVectorOfMapKeys(bedsByPanel);
	njh::sort(panelNames);
	std::vector<std::map<std::string, double>> coverageForEachPanel;
	std::map<std::string, double> panelCovTemplate;
	for(const auto & panel : panelNames){
		panelCovTemplate[panel] = 0;
	}
	for(uint32_t pos = 0; pos < regions.size(); ++pos){
		coverageForEachPanel.emplace_back(panelCovTemplate);
	}

	for(const auto & bedByPanel : bedsByPanel){
		auto panelRegions = getBeds(bedByPanel.second);

		for(const auto bedPos : iter::range(regions.size())){
			const auto & bed = regions[bedPos];
			std::vector<Bed6RecordCore> overlappingRegions;
			for(const auto & panelBed : panelRegions){
				if(bed->overlaps(*panelBed, 1)){
					overlappingRegions.emplace_back(*panelBed);
				}
			}
			double maxCov = 0;

			for(const auto & overlap : overlappingRegions){
				double cov = bed->getOverlapLen(overlap)/static_cast<double>(bed->length());
				if(!doNotAdjustForSinglePosIntersection){
					if(overlap.length() <= 3){
						cov =  bed->getOverlapLen(overlap)/static_cast<double>(overlap.length());
					}
				}
				if(cov > maxCov){
					maxCov = cov;
				}
			}
			coverageForEachPanel[bedPos][bedByPanel.first] = maxCov;
		}
	}
	out << "#chrom\tstart\tend\tname\tlen\tstrand";
	for(uint32_t field = 0; field < maxExtraFields; ++field){
		out << "\t" << "extrafield" << field;
	}
	out << "\t" << njh::conToStr(panelNames, "\t") << std::endl;
	for(const auto bedPos : iter::range(regions.size())){
		out << regions[bedPos]->toDelimStr();
		for(uint32_t field = 0; field < maxExtraFields; ++field){
			out << "\t";
			if(field < regions[bedPos]->extraFields_.size()){
				out << regions[bedPos]->extraFields_[field];
			}else{
				out << "NA";
			}
		}
		for(const auto & panel : panelNames){
			out << "\t" << coverageForEachPanel[bedPos][panel];
		}
		out << std::endl;
	}

	return 0;
}

} // namespace njheq


