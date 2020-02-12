/*
 * bedExp_roughSmoothingForBedCoverage.cpp
 *
 *  Created on: May 27, 2019
 *      Author: nicholashathaway
 */

#include "bedExp.hpp"

#include "elucidator/BamToolsUtils.h"
#include "elucidator/objects/BioDataObject.h"


namespace njhseq {



int bedExpRunner::roughSmoothingForBedCoverage(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path coverageFile = "";
	OutOptions outOpts(bfs::path(""), ".bed");
	uint32_t within = 10000;
	uint32_t numThreads = 1;
	seqSetUp setUp(inputCommands);
	setUp.description_ = "Some base coverage smoothing for output by bamMulticovBases";
	setUp.setOption(coverageFile, "--coverage", "bed file with coverage information", true);
	setUp.setOption(numThreads, "--numThreads", "number Threads", njh::progutils::ProgramSetUp::CheckCase::NONZERO);
	setUp.setOption(within, "--within", "Use windows within this many bp to smooth coverage", njh::progutils::ProgramSetUp::CheckCase::NONZERO);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	table inTable(coverageFile, "\t", true);
	if(inTable.nCol() < 7){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "in table should at least be 7 not " << inTable.nCol()<< "\n";
		throw std::runtime_error{ss.str()};
	}

	OutputStream out(outOpts);
	VecStr samples = getSubVector(inTable.columnNames_, 6);

	std::vector<GenomicRegion> regions;
	table outTable(inTable.columnNames_);
	std::unordered_map<std::string, std::vector<double>> coverageBySample;
	for(const auto & row : inTable){
		outTable.addRow(row);
		regions.emplace_back(GenomicRegion(row[3], row[0],
				njh::StrToNumConverter::stoToNum<uint32_t>(row[1]),
				njh::StrToNumConverter::stoToNum<uint32_t>(row[2]), "-" == row[5]));
	}
	std::unordered_map<std::string, uint32_t> sampToColPos;
	for(const auto & colPos : iter::range<uint32_t>(6,inTable.nCol())){
		coverageBySample[inTable.columnNames_[colPos]] = vecStrToVecNum<double>(inTable.getColumn(colPos));
		sampToColPos[inTable.columnNames_[colPos]] = colPos;
	}

//	uint32_t testNum = 10;
//	uint32_t num = 0;

	struct WindowCovStat {
		WindowCovStat(double meanCov, double stdCov) :
				meanCov_(meanCov), stdCov_(stdCov) {

		}
		double meanCov_;
		double stdCov_;
	};


	std::vector<uint32_t> regionsPos(regions.size());
	njh::iota<uint32_t>(regionsPos, 0);
	njh::concurrent::LockableVec<uint32_t> regionsPosQueue(regionsPos);

	std::function<void()> roundRegions = [&regionsPosQueue,&regions,
											 &outTable,&sampToColPos,&coverageBySample,
											 &samples,&within](){
		uint32_t pos = std::numeric_limits<uint32_t>::max();
		while(regionsPosQueue.getVal(pos)){
			//		std::cout << "\r" << pos << "/" << regions.size();
			//		std::cout.flush();
					std::vector<uint32_t> otherRegions;
					for(const auto & otherPos : iter::range(regions.size())){
						if(regions[otherPos].chrom_ == regions[pos].chrom_){
							if(regions[otherPos].distBetweenRegions(regions[pos]) <= within){
								otherRegions.emplace_back(otherPos);
							}
							if(regions[otherPos].start_ > regions[pos].end_ && regions[otherPos].start_ - regions[pos].end_ > within){
								break;
							}
						}
					}
			//		std::cout << "pos: " << pos << std::endl;
			//		std::cout << regions[pos].genBedRecordCore().toDelimStr() << std::endl;
			//		for(const auto & otherPos : otherRegions){
			//			std::cout << "otherPos: " << otherPos << std::endl;
			//			std::cout << "\t" << regions[otherPos].genBedRecordCore().toDelimStr() << std::endl;
			//			std::cout << "\tdist: " << distBetweenRegions(regions[otherPos], regions[pos]) << std::endl;
			//		}
					std::vector<std::vector<uint32_t>> windows;
					for(const auto & posInOthers : iter::range(otherRegions.size())){
						uint32_t positionInRegions = otherRegions[posInOthers];
						if(positionInRegions > pos){
							break;
						}
						std::vector<uint32_t> window{positionInRegions};
						for(const auto & selectPos : iter::range(posInOthers + 1, otherRegions.size())){
							uint32_t nextRegion_positionInRegions = otherRegions[selectPos];
							auto dist = regions[positionInRegions].distBetweenRegions(regions[nextRegion_positionInRegions]);
							if(dist <=within){
								window.emplace_back(nextRegion_positionInRegions);
							}
							if(dist > within){
								break;
							}
						}
						windows.emplace_back(window);
					}
					std::unordered_map<std::string, std::vector<WindowCovStat>> statsPerWindowPerSamp;
					for(const auto & windowNumber : iter::range(windows.size())){
			//			std::cout << "window: " << windowNumber << ", size: " << windows[windowNumber].size()
			//					<< ", diff between front and back: "
			//					<< distBetweenRegions(regions[windows[windowNumber].front()], regions[windows[windowNumber].back()])<< std::endl;
			//			for(const auto & windowPos : windows[windowNumber]){
			//				std::cout << "\t";
			//				if(windowPos == pos){
			//					std::cout << njh::bashCT::green;
			//				}
			//				std::cout << regions[windowPos].genBedRecordCore().toDelimStr() ;
			//				if(windowPos == pos){
			//					std::cout << njh::bashCT::reset;
			//				}
			//				std::cout << std::endl;
			//			}
						//coverageBySample

						for (const auto samp : samples) {
							std::vector<double> covs;
			//				std::cout << samp << ":";
							for (const auto & windowPos : windows[windowNumber]) {
								covs.emplace_back(coverageBySample[samp][windowPos]/regions[windowPos].getLen());
			//					std::cout << " " << coverageBySample[samp][windowPos];
							}
							statsPerWindowPerSamp[samp].emplace_back(vectorMean(covs), vectorStandardDeviationPop(covs));
			//				std::cout << std::endl;
						}
					}
					for(const auto & samp : statsPerWindowPerSamp){
						double minStd = std::numeric_limits<double>::max();
						double avgPerBaseCovForMinStd = 0;
						for(const auto & windowStat : samp.second){
			//				std::cout << samp.first << "\t" << windowStat.meanCov_ << "\t" << windowStat.stdCov_ << std::endl;
							if(windowStat.stdCov_ < minStd){
								minStd = windowStat.stdCov_;
								avgPerBaseCovForMinStd = windowStat.meanCov_;
							}
						}
			//			std::cout << samp.first << " winner: " << avgPerBaseCovForMinStd << " " << regions[pos].getLen() * avgPerBaseCovForMinStd << std::endl;
			//			std::cout << "new: " << regions[pos].getLen() * avgPerBaseCovForMinStd << std::endl;
			//			std::cout << "was: " << coverageBySample[samp.first][pos] << std::endl;

						if(coverageBySample[samp.first][pos] > 0){
							outTable.content_[pos][sampToColPos[samp.first]] = estd::to_string(regions[pos].getLen() * avgPerBaseCovForMinStd);
						}
					}
			//		++num;
			//		if(num > testNum){
			//			break;
			//		}
		}
	};

	njh::concurrent::runVoidFunctionThreaded(roundRegions, numThreads);
	outTable.outPutContents(out, "\t");
	return 0;
}


}  // namespace njhseq

