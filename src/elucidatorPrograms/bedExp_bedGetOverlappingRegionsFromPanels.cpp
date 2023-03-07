/*
 * bedExp_bedGetOverlappingRegionsFromPanels.cpp
 *
 *  Created on: Jul 16, 2021
 *      Author: nick
 */





#include "bedExp.hpp"
#include <njhseq/objects/BioDataObject.h>

#include "elucidator/objects/BioDataObject.h"

#include <utility>
#include "elucidator/BioRecordsUtils/BedUtility.hpp"
#include <njhseq/objects/counters/DNABaseCounter.hpp>


namespace njhseq {


int bedExpRunner::bedGetOverlappinPositionsBetweenPanels(const njh::progutils::CmdArgs & inputCommands) {
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
	uint32_t bedColPos = otherBedsTab.getColPos("bed");
	std::unordered_map<std::string, std::string> panelByBeds;
	for(const auto & row : otherBedsTab){
		auto panel = row[panelPos];
		auto bed = row[bedColPos];
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
	out << "panel1";
	out << "\t" << "panel1_chrom\tpanel1_start\tpanel1_end\tpanel1_name\tpanel1_len\tpanel1_strand";
	out << "\t" << "panel2";
	out << "\t" << "panel2_chrom\tpanel2_start\tpanel2_end\tpanel2_name\tpanel2_len\tpanel2_strand";
	out << "\t" << "overlap_chrom\toverlap_start\toverlap_end\toverlap_name\toverlap_len\toverlap_strand";
	out << "\t" << "panel1_coverage" << "\t" << "panel2_coverage" << std::endl;
	for(const auto & currentPanel : panelNames){
		auto currenntPanelBed = bedsByPanel[currentPanel];
		auto regions = getBeds(currenntPanelBed);
		for(const auto & bedByPanel : bedsByPanel){
			if(currentPanel == bedByPanel.first){
				continue;
			}
			auto panelRegions = getBeds(bedByPanel.second);
			for(const auto bedPos : iter::range(regions.size())){
				const auto & bed = regions[bedPos];
				std::vector<std::shared_ptr<Bed6RecordCore>> overlapingRegions;
				for(const auto & panelBed : panelRegions){
					if(bed->overlaps(*panelBed, 1)){
						overlapingRegions.emplace_back(panelBed);
					}
				}
				for(const auto & overlap : overlapingRegions){
					out << currentPanel;
					out << "\t" << bed->toDelimStr();
					out << "\t" << bedByPanel.first;
					out << "\t" << overlap->toDelimStr();
					uint32_t overlapStart = std::max(overlap->chromStart_, bed->chromStart_);
					uint32_t overlapStop = std::min(overlap->chromEnd_, bed->chromEnd_);
					uint32_t overlapLen = overlapStop - overlapStart;
					out << "\t" << overlap->chrom_ << "\t" << overlapStart << "\t" << overlapStop << "\t"
							<< njh::pasteAsStr(overlap->chrom_, "-", overlapStart, "-", overlapStop) << "\t" << overlapLen << "\t"
							<< overlap->strand_;
					out << "\t" << static_cast<double>(overlapLen)/bed->length() << "\t" << static_cast<double>(overlapLen)/overlap->length() << std::endl;
				}
			}
		}
	}
	return 0;
}

int bedExpRunner::bedGetOverlappinCoverageBetweenPanels(const njh::progutils::CmdArgs & inputCommands) {
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
	uint32_t bedColPos = otherBedsTab.getColPos("bed");
	std::unordered_map<std::string, std::string> panelByBeds;
	for(const auto & row : otherBedsTab){
		auto panel = row[panelPos];
		auto bed = row[bedColPos];
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

		struct RegionWithOtherRegionsCoverage{
			explicit RegionWithOtherRegionsCoverage(Bed6RecordCore region): region_(std::move(region)){

			}
			Bed6RecordCore region_;
			std::vector<Bed6RecordCore> overlappingRegions_;

			[[nodiscard]] uint32_t getNumberOfBasesCovered() const {
				uint32_t ret = 0;
				std::unordered_map<uint32_t, uint32_t> coveragePerPosition;
				for(const auto pos : iter::range(region_.chromStart_, region_.chromEnd_)){
					coveragePerPosition[pos] = 0;
				}
				for(const auto & overlap : overlappingRegions_){
					for(const auto pos : iter::range(std::max(overlap.chromStart_, region_.chromStart_), std::min(overlap.chromEnd_, region_.chromEnd_))){
						coveragePerPosition[pos] += 1;
					}
				}
				for(const auto & cov : coveragePerPosition){
					if(cov.second > 0){
						++ret;
					}
				}
				return ret;
			}
			double getFractionCovered() const{
				return static_cast<double>(getNumberOfBasesCovered())/region_.length();
			}
			uint32_t getNumOfOverlaps() const{
				return overlappingRegions_.size();
			}
		};
		std::unordered_map<std::string, std::vector<RegionWithOtherRegionsCoverage>> overlappingRegionsPerPanel;

		for(const auto & bedByPanel : bedsByPanel){
			auto panelRegions = getBeds(bedByPanel.second);
			for(const auto bedPos : iter::range(regions.size())){
				overlappingRegionsPerPanel[bedByPanel.first].emplace_back(*regions[bedPos]);
			}
			for(const auto bedPos : iter::range(regions.size())){
				const auto & bed = regions[bedPos];
				for(const auto & panelBed : panelRegions){
					if(bed->overlaps(*panelBed, 1)){
						overlappingRegionsPerPanel[bedByPanel.first][bedPos].overlappingRegions_.emplace_back(*panelBed);
					}
				}
			}
		}

		for(const auto bedPos : iter::range(regions.size())){
			out << regions[bedPos]->toDelimStr();
			out << "\t" << currentPanel;
			for(const auto & panel : panelNames){
				out << "\t" << overlappingRegionsPerPanel[panel][bedPos].getFractionCovered();
			}
			out << std::endl;
		}
	}
	return 0;
}

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


