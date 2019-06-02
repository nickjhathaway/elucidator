/*
 * bedExp_fillInRegionsByBestScore.cpp
 *
 *  Created on: May 29, 2019
 *      Author: nicholashathaway
 */



#include "bedExp.hpp"

#include "elucidator/BamToolsUtils.h"
#include "elucidator/objects/BioDataObject.h"
#include "elucidator/BioRecordsUtils/BedUtility.hpp"


namespace njhseq {


template<typename T>
std::vector<GenomicRegion> mergeAndSort(std::vector<T> & beds) {
	BedUtility::coordSort(beds, false);
	std::vector<GenomicRegion> ret;
	if (beds.size() > 1) {
		ret.emplace_back(GenomicRegion(getRef(beds.front())));
		for (const auto & regPos : iter::range<uint32_t>(1, beds.size())) {
			//merge if they overlap or start where the last region starts
			if (ret.back().overlaps(getRef(beds[regPos])) || ret.back().end_ == getRef(beds[regPos]).chromStart_) {
				ret.back().end_ = getRef(beds[regPos]).chromEnd_;
			} else {
				ret.emplace_back(GenomicRegion(getRef(beds[regPos])));
			}
		}
	} else if (1 == beds.size()) {
		ret.emplace_back(GenomicRegion(getRef(beds.front())));
	} else {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "no regions read in input" << "\n";
		throw std::runtime_error { ss.str() };
	}
	return ret;
}


std::vector<GenomicRegion> readInBedsMergeAndSort(const bfs::path & bedFile) {
	auto inRegions = getBed3s(bedFile);
	return mergeAndSort(inRegions);
}



int bedExpRunner::getInterveningRegions(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bedFile = "";
	bfs::path genomeFnp = "";
	OutOptions outOpts(bfs::path(""), ".bed");
	uint32_t padding = 0;
	seqSetUp setUp(inputCommands);
	setUp.setOption(bedFile, "--bed", "Bed file", true);
	setUp.setOption(padding, "--padding", "Padding to add to the end and beginning of regions");
	setUp.setOption(genomeFnp, "--genome", "Genome file, this is used to determine the length of the chromosome", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);
	//get chrom lengths
	std::unordered_map<std::string, uint32_t> chromosomeLengths;
	{
		seqInfo seq;
		auto genomeOpts = SeqIOOptions::genFastaIn(genomeFnp);
		genomeOpts.includeWhiteSpaceInName_ = false;
		SeqInput reader(genomeOpts);
		reader.openIn();
		while(reader.readNextRead(seq)){
			chromosomeLengths[seq.name_] = len(seq);
		}
	}
	OutputStream out(outOpts);
	auto inRegions = getBed3s(bedFile);
	for(auto & reg : inRegions){
		BedUtility::extendLeftRight(*reg, padding, padding, njh::mapAt(chromosomeLengths,reg->chrom_));
	}
	auto mergedRegions = mergeAndSort(inRegions);
	//do the first region from beginning of chromosome
	if(0 != mergedRegions.front().start_){
		GenomicRegion front("", mergedRegions.front().chrom_, 0, mergedRegions.front().start_, false);
		out << front.genBed3RecordCore().toDelimStr() << std::endl;
	}
	if(mergedRegions.size() > 1){
		for(const auto regPos : iter::range<uint32_t>(1, mergedRegions.size())){
			if (mergedRegions[regPos - 1].chrom_ == mergedRegions[regPos].chrom_) {
				//on same chromosome, just do the region inbetween the two regions
				GenomicRegion inbetween("", mergedRegions[regPos].chrom_, mergedRegions[regPos-1].end_, mergedRegions[regPos].start_, false);
				out << inbetween.genBed3RecordCore().toDelimStr() << std::endl;
			} else {
				//start of a new chromsome need to add the back of the last chromosome and the front of current chromosome
				if(mergedRegions[regPos-1].end_ != njh::mapAt(chromosomeLengths,mergedRegions[regPos-1].chrom_)){
					GenomicRegion back("", mergedRegions[regPos-1].chrom_,mergedRegions[regPos-1].end_, njh::mapAt(chromosomeLengths,mergedRegions[regPos-1].chrom_), false);
					out << back.genBed3RecordCore().toDelimStr() << std::endl;
				}
				if(0 != mergedRegions[regPos].start_){
					GenomicRegion front("", mergedRegions[regPos].chrom_, 0, mergedRegions[regPos].start_, false);
					out << front.genBed3RecordCore().toDelimStr() << std::endl;
				}
			}
		}
	}
	//do the last region to the end of the chromosome
	if(mergedRegions.back().end_ != njh::mapAt(chromosomeLengths,mergedRegions.back().chrom_)){
		GenomicRegion back("", mergedRegions.back().chrom_,mergedRegions.back().end_, njh::mapAt(chromosomeLengths,mergedRegions.back().chrom_), false);
		out << back.genBed3RecordCore().toDelimStr() << std::endl;
	}


	return 0;
}






int bedExpRunner::fillInRegionsByBestScore(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path withinBed = "";
	bfs::path bedFnp = "";
	uint32_t within = 100000;
	OutOptions outOpts(bfs::path(""), ".bed");
	seqSetUp setUp(inputCommands);
	setUp.description_ = "";
	setUp.setOption(withinBed, "--withinBed", "Pick Regions that intersect with these regions", true);
	setUp.setOption(bedFnp, "--bedFnp", "Bed file to choose from", true);
	setUp.setOption(within, "--within", "within", false);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	auto bedsToSelect = getBeds(bedFnp);
	auto withinRegions = bed3PtrsToGenomicRegs(getBed3s(withinBed));

	BedUtility::coordSort(bedsToSelect);
	OutputStream out(outOpts);

	while(!withinRegions.empty()){

		std::unordered_map<std::string, std::vector<std::shared_ptr<Bed6RecordCore>>> byRegions;
		std::unordered_map<std::string, GenomicRegion> regionKey;

		for(const auto & reg : withinRegions){
			regionKey[reg.createUidFromCoords()] = reg;
			for(const auto & selectBed : bedsToSelect){
				if(selectBed->chrom_ == reg.chrom_){
					if(reg.overlaps(*selectBed)){
						byRegions[reg.createUidFromCoords()].emplace_back(selectBed);
					}
				}
			}
		}
		withinRegions.clear();
		//if byRegions is empty there were no overlapping regions
		if(byRegions.empty()){
			break;
		}

		for (auto & byRegion : byRegions) {
			if (byRegion.second.size() > 1) {
				njh::sort(byRegion.second,[](const std::shared_ptr<Bed6RecordCore> & reg1, const std::shared_ptr<Bed6RecordCore> & reg2){
					return reg1->score_ > reg2->score_;
				});
				auto outBed = *byRegion.second.front();
				outBed.extraFields_.emplace_back(byRegion.first);
				out << outBed.toDelimStrWithExtra() << std::endl;
				//front
				if(regionKey[byRegion.first].start_ + within < outBed.chromStart_){
					GenomicRegion front = regionKey[byRegion.first];
					front.end_ = outBed.chromStart_ - within;
					withinRegions.emplace_back(front);
				}
				if(regionKey[byRegion.first].end_ > within && regionKey[byRegion.first].end_ - within > outBed.chromEnd_){
					GenomicRegion back = regionKey[byRegion.first];
					back.start_ = outBed.chromEnd_ + within;
					withinRegions.emplace_back(back);
				}
			} else {
				auto outBed = *byRegion.second.front();
				outBed.extraFields_.emplace_back(byRegion.first);
				out << outBed.toDelimStrWithExtra() << std::endl;
			}
		}
	}


	return 0;
}





}  // namespace njhseq


