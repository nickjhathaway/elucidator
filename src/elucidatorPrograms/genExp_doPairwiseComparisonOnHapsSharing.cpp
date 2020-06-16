/*
 * genExp_doPairwiseComparisonOnHapsSharing.cpp
 *
 *  Created on: Jun 5, 2020
 *      Author: nick
 */


// elucidator - A library for analyzing sequence data
// Copyright (C) 2012-20200 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
//
// This file is part of elucidator.
//
// elucidator is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// elucidator is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with elucidator.  If not, see <http://www.gnu.org/licenses/>.
//

#include <TwoBit.h>
#include "genExp.hpp"
#include "elucidator/BamToolsUtils.h"
#include "elucidator/objects/BioDataObject.h"
#include "elucidator/objects/dataContainers/graphs.h"
#include <njhseq/GenomeUtils.h>



namespace njhseq {

//s_Sample s_studyid p_name                 h_popUID                  c_relAbund

int genExpRunner::doPairwiseComparisonOnHapsSharing(const njh::progutils::CmdArgs & inputCommands){

	bfs::path tableFnp = "";

	std::string sampleCol = "";
	std::string targetNameCol = "";
	std::string popIDCol = "";
	std::string relAbundCol = "";

	uint32_t numThreads = 1;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();

	setUp.setOption(numThreads, "--numThreads", "number of cpus to use");
  setUp.setOption(tableFnp, "--tableFnp", "Table to read in (should be tab delimited)", true);

  setUp.setOption(sampleCol, "--sampleCol", "sampleCol", true);
  setUp.setOption(targetNameCol, "--targetNameCol", "targetNameCol", true);
  setUp.setOption(popIDCol, "--popIDCol", "popIDCol", true);
  setUp.setOption(relAbundCol, "--relAbundCol", "relAbundCol", true);

  setUp.processDirectoryOutputName(bfs::path(bfs::basename(tableFnp)).string() + "_doPairwiseComparisonOnHapsSharing_TODAY", true);
	setUp.finishSetUp(std::cout);

	setUp.startARunLog(setUp.pars_.directoryName_);

	table hapTab(tableFnp, "\t", true);
	hapTab.checkForColumnsThrow(VecStr{sampleCol, targetNameCol, popIDCol, relAbundCol}, __PRETTY_FUNCTION__);

	struct SampPopulation{

		struct Tar{
			Tar(){

			}
			Tar(const std::unordered_map<uint32_t, double> & haps): haps_(haps){

			}
			std::unordered_map<uint32_t, double> haps_;

			uint32_t hapCount()const{
				return haps_.size();
			}



			uint32_t hapsShared(const std::unordered_map<uint32_t, double> & otherHaps) const{
				uint32_t ret = 0;
				for(const auto & otherHap : otherHaps){
					if(njh::in(otherHap.first, haps_)){
						++ret;
					}
				}
				return ret;
			}
		};

		struct Samp{
			Samp(){

			}
			Samp(const std::unordered_map<uint32_t, Tar> & targetsWithHaps): targetsWithHaps_(targetsWithHaps){

			}
			std::unordered_map<uint32_t, Tar> targetsWithHaps_;

			struct SampComp {
				uint32_t targetsShared_{0};
				uint32_t targetsPossibleToShare_{0};
				uint32_t hapsShared_{0};
				uint32_t hapsPossibleToShareInSharedTargets_{0};

				std::vector<double> jacardsByTargets_;

				double percTargetsWithAtLeastOneHap() const {
					return
							targetsPossibleToShare_ > 0 ?
									targetsShared_
											/ static_cast<double>(targetsPossibleToShare_) :
									0.0;
				}
				double jacardByHaps() const {
					return
							hapsShared_ > 0 ?
									hapsPossibleToShareInSharedTargets_
											/ static_cast<double>(hapsPossibleToShareInSharedTargets_) :
									0.0;
				}

				double avgJacard() const {
					return vectorMean(jacardsByTargets_);
				}
			};
			SampComp compToOtherSamp(
					const std::unordered_map<uint32_t, Tar> &otherTargets) const {
				SampComp ret;
				for (const auto &other : otherTargets) {
//					std::cout << other.first << std::endl;
					if (njh::in(other.first, targetsWithHaps_)) {
						++ret.targetsPossibleToShare_;
						auto hapsShared = targetsWithHaps_.at(other.first).hapsShared(
								other.second.haps_);
						ret.hapsShared_ += hapsShared;
						if (hapsShared > 0) {
							++ret.targetsShared_;
						}
						std::unordered_set<uint32_t> haps;
						std::vector<uint32_t> currentHaps = getVectorOfMapKeys(targetsWithHaps_);
						haps.insert(currentHaps.begin(), currentHaps.end());
						std::vector<uint32_t> currentOtherHaps = getVectorOfMapKeys(otherTargets);
						haps.insert(currentOtherHaps.begin(), currentOtherHaps.end());
						ret.hapsPossibleToShareInSharedTargets_ = haps.size();
						ret.jacardsByTargets_.emplace_back(static_cast<double>(hapsShared)/haps.size());
					}
				}
				return ret;
			}
			SampComp compToOtherSamp(const Samp &otherSamp) const {
				return compToOtherSamp(otherSamp.targetsWithHaps_);
			}

		};
		VecStr getSampNames()const{
			return getVectorOfMapKeys(samples_);
		}

		std::unordered_map<std::string, Samp> samples_;
	};



	SampPopulation pop;

	std::unordered_set<std::string> tarNames;
	std::unordered_map<std::string, std::unordered_set<std::string>> hapNamesForTars;

	std::unordered_map<std::string, uint32_t> tarNameKey;
	std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>> hapNamesKey;

	uint64_t totalHaps{0};
	{
		for(const auto & row : hapTab){
			auto tar = row[hapTab.getColPos(targetNameCol)];
			auto hapName = row[hapTab.getColPos(popIDCol)];
			tarNames.emplace(tar);
			hapNamesForTars[tar].emplace(hapName);
		}
		std::vector<std::string> tarNamesVec{tarNames.begin(), tarNames.end()};
		for(const auto pos : iter::range(tarNames.size())){
			tarNameKey[tarNamesVec[pos]] = pos;
		}

		for(const auto & hapName : hapNamesForTars){
			std::vector<std::string> hapNamesVec{hapName.second.begin(), hapName.second.end()};
			for(const auto pos : iter::range(hapNamesVec.size())){
				hapNamesKey[hapName.first][hapNamesVec[pos]] = pos;
				++totalHaps;
			}
		}
	}

	std::cout << totalHaps << std::endl;

	for(const auto & row : hapTab){
		auto samp = row[hapTab.getColPos(sampleCol)];
		auto tar = row[hapTab.getColPos(targetNameCol)];
		auto hapName = row[hapTab.getColPos(popIDCol)];
		auto rBund = row[hapTab.getColPos(relAbundCol)];
		pop.samples_[samp].targetsWithHaps_[tarNameKey[tar]].haps_[hapNamesKey[tar][hapName]] += njh::StrToNumConverter::stoToNum<double>(rBund);
	}




	auto sampNames = pop.getSampNames();


	std::vector<std::vector<double>> byTarget{sampNames.size(), std::vector<double>(sampNames.size())};
	std::vector<std::vector<double>> byHap{sampNames.size(), std::vector<double>(sampNames.size())};
	std::vector<std::vector<double>> avgJacard{sampNames.size(), std::vector<double>(sampNames.size())};

	for(uint32_t pos = 0; pos < sampNames.size(); ++pos){
		byTarget[pos][pos] = 1;
		byHap[pos][pos] = 1;
		avgJacard[pos][pos] = 1;
	}
	PairwisePairFactory pFactor(sampNames.size());


	njh::ProgressBar progpar(pFactor.totalCompares_);

	std::function<void()> compSamps = [&pFactor,&byTarget,&byHap,&pop,&sampNames,&avgJacard](){
		PairwisePairFactory::PairwisePairVec pairVec;

		while(pFactor.setNextPairs(pairVec, 100)){
			for(const auto & pairPos : iter::range(pairVec.pairs_.size())){
				const auto & pair = pairVec.pairs_[pairPos];
//				if(setUp.pars_.verbose_){
//					progpar.outputProgAdd(std::cout, 1, true);
//					std::cout << sampNames[pair.row_] << " vs " << sampNames[pair.col_] << std::endl;
//				}
				auto compRes = pop.samples_[sampNames[pair.row_]].compToOtherSamp(pop.samples_[sampNames[pair.col_]]);
				byTarget[pair.row_][pair.col_] = compRes.percTargetsWithAtLeastOneHap();
				byHap[pair.row_][pair.col_] = compRes.jacardByHaps();
				avgJacard[pair.row_][pair.col_] = compRes.avgJacard();

				byTarget[pair.col_][pair.row_] = compRes.percTargetsWithAtLeastOneHap();
				byHap[pair.col_][pair.row_] = compRes.jacardByHaps();
				avgJacard[pair.col_][pair.row_] = compRes.avgJacard();
			}
		}
	};


	njh::concurrent::runVoidFunctionThreaded(compSamps, numThreads);



	OutputStream outSampNamesOut(njh::files::make_path(setUp.pars_.directoryName_, "sampleNames.tab.txt"));
	OutputStream byTargetOut(njh::files::make_path(setUp.pars_.directoryName_, "percOfTarSharingAtLeastOneHap.tab.txt.gz"));
	OutputStream byHapOut(njh::files::make_path(setUp.pars_.directoryName_, "jacardByAllHap.tab.txt.gz"));
	OutputStream avgHapOut(njh::files::make_path(setUp.pars_.directoryName_, "avgJacardPerTarget.tab.txt.gz"));

	outSampNamesOut << njh::conToStr(sampNames, "\n") << std::endl;
	for(const auto & it : byTarget){
		byTargetOut << njh::conToStr(it, "\t") << std::endl;
	}
	for(const auto & ih : byHap){
		byHapOut << njh::conToStr(ih, "\t") << std::endl;
	}
	for(const auto & ih : avgJacard){
		avgHapOut << njh::conToStr(ih, "\t") << std::endl;
	}


	return 0;
}


int genExpRunner::doPairwiseComparisonOnHapsSharingDev(const njh::progutils::CmdArgs & inputCommands){

	bfs::path tableFnp = "";

	std::string sampleCol = "";
	std::string targetNameCol = "";
	std::string popIDCol = "";
	std::string relAbundCol = "";

	uint32_t numThreads = 1;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();

	setUp.setOption(numThreads, "--numThreads", "number of cpus to use");
  setUp.setOption(tableFnp, "--tableFnp", "Table to read in (should be tab delimited)", true);

  setUp.setOption(sampleCol, "--sampleCol", "sampleCol", true);
  setUp.setOption(targetNameCol, "--targetNameCol", "targetNameCol", true);
  setUp.setOption(popIDCol, "--popIDCol", "popIDCol", true);
  setUp.setOption(relAbundCol, "--relAbundCol", "relAbundCol", true);

  setUp.processDirectoryOutputName(bfs::path(bfs::basename(tableFnp)).string() + "_doPairwiseComparisonOnHapsSharing_TODAY", true);
	setUp.finishSetUp(std::cout);

	setUp.startARunLog(setUp.pars_.directoryName_);

	table hapTab(tableFnp, "\t", true);
	hapTab.checkForColumnsThrow(VecStr{sampleCol, targetNameCol, popIDCol, relAbundCol}, __PRETTY_FUNCTION__);


	std::unordered_set<std::string> sampNames;
	std::unordered_map<std::string, uint32_t> sampNamesKey;

	std::unordered_set<std::string> tarNames;
	std::unordered_map<std::string, uint32_t> tarNameKey;

	std::unordered_map<std::string, std::unordered_set<std::string>> hapNamesForTars;
	std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>> hapNamesKey;



	uint64_t totalHaps{0};
	for(const auto & row : hapTab){
		auto samp = row[hapTab.getColPos(sampleCol)];
		auto tar = row[hapTab.getColPos(targetNameCol)];
		auto hapName = row[hapTab.getColPos(popIDCol)];

		sampNames.emplace(samp);
		tarNames.emplace(tar);
		hapNamesForTars[tar].emplace(hapName);
	}

	std::vector<std::string> sampNamesVec{sampNames.begin(), sampNames.end()};
	std::vector<std::string> tarNamesVec{tarNames.begin(), tarNames.end()};
	njh::sort(sampNamesVec);
	njh::sort(tarNamesVec);



	for(const auto pos : iter::range(tarNamesVec.size())){
		tarNameKey[tarNamesVec[pos]] = pos;
	}

	for(const auto pos : iter::range(sampNamesVec.size())){
		sampNamesKey[sampNamesVec[pos]] = pos;
	}

	for(const auto & hapName : hapNamesForTars){
		std::vector<std::string> hapNamesVec{hapName.second.begin(), hapName.second.end()};
		for(const auto pos : iter::range(hapNamesVec.size())){
			hapNamesKey[hapName.first][hapNamesVec[pos]] = pos;
			++totalHaps;
		}
	}
	std::vector<uint32_t> numberOfHapsPerTarget(tarNamesVec.size());
	for(const auto & tarHaps : hapNamesForTars){
		numberOfHapsPerTarget[tarNameKey[tarHaps.first]] = tarHaps.second.size();
	}
	std::vector<uint32_t> tarStart(tarNamesVec.size());

	{

		uint32_t pos = 0;
		for(const auto & tarPos : iter::range(tarNamesVec.size())){
			tarStart[tarPos] = pos;
			pos += numberOfHapsPerTarget[tarPos];
		}
	}

	std::cout << "totalHaps: " << totalHaps << std::endl;
//	std::cout << "tarStart: " << njh::conToStr(tarStart, ",") << std::endl;

	std::vector<std::vector<uint8_t>> hapsEncodeBySamp(sampNames.size());
	for(const auto & samp : sampNames){
		hapsEncodeBySamp[sampNamesKey[samp]] = std::vector<uint8_t>(totalHaps, 0);
	}

	std::vector<std::vector<uint8_t>> targetsEncodeBySamp(sampNames.size());
	for(const auto & samp : sampNames){
		targetsEncodeBySamp[sampNamesKey[samp]] = std::vector<uint8_t>(tarNamesVec.size(), 0);
	}

	for(const auto & row : hapTab){
		auto samp = row[hapTab.getColPos(sampleCol)];
		auto tar = row[hapTab.getColPos(targetNameCol)];
		auto hapName = row[hapTab.getColPos(popIDCol)];
		auto rBund = row[hapTab.getColPos(relAbundCol)];

		auto tKey = tarNameKey[tar];
		auto hKey = hapNamesKey[tar][hapName];
//		std::cout << "sampNamesKey[samp]                        : " << sampNamesKey[samp] << std::endl;
//		std::cout << "tKey                                      : " << tKey << std::endl;
//		std::cout << "[tarStart[tKey]                           : " << tarStart[tKey] << std::endl;
//		std::cout << "hKey                                      : " << hKey << std::endl;
//		std::cout << "hKey                                      : " << hKey << std::endl;
//		std::cout << "hapsEncodeBySamp[sampNamesKey[samp]]size(): " << hapsEncodeBySamp[sampNamesKey[samp]].size() << std::endl;
		hapsEncodeBySamp[sampNamesKey[samp]][tarStart[tKey] + hKey] = 1;
		targetsEncodeBySamp[sampNamesKey[samp]][tKey] = 1;
		//pop.samples_[samp].targetsWithHaps_[tarNameKey[tar]].haps_[hapNamesKey[tar][hapName]] += njh::StrToNumConverter::stoToNum<double>(rBund);
	}


	std::vector<std::vector<double>> byAllHaps{sampNames.size(), std::vector<double>(sampNames.size())};
	std::vector<std::vector<double>> byHapsTarShared{sampNames.size(), std::vector<double>(sampNames.size())};

	std::vector<std::vector<double>> avgJacard{sampNames.size(), std::vector<double>(sampNames.size())};
	std::vector<std::vector<double>> byTarget{sampNames.size(), std::vector<double>(sampNames.size())};

	for(uint32_t pos = 0; pos < sampNames.size(); ++pos){
		byHapsTarShared[pos][pos] = 1;
		byAllHaps[pos][pos] = 1;

		avgJacard[pos][pos] = 1;
		byTarget[pos][pos] = 1;
	}
	PairwisePairFactory pFactor(sampNames.size());


	njh::ProgressBar progpar(pFactor.totalCompares_);
	if(setUp.pars_.verbose_){
		std::cout << "totalComps: " << pFactor.totalCompares_ << std::endl;
	}

	std::function<void()> compSamps = [&pFactor,&byHapsTarShared,&byAllHaps,&targetsEncodeBySamp,&progpar,
																		 &totalHaps,&hapsEncodeBySamp,&setUp,&tarNamesVec,&byTarget,
																		 &numberOfHapsPerTarget,&tarStart,
																		 &avgJacard](){
		PairwisePairFactory::PairwisePairVec pairVec;
		while(pFactor.setNextPairs(pairVec, 100)){
			for(const auto & pairPos : iter::range(pairVec.pairs_.size())){
				const auto & pair = pairVec.pairs_[pairPos];
				if(setUp.pars_.verbose_){
					progpar.outputProgAdd(std::cout, 1, true);
				}


				//for shared targets
				{
					uint32_t totalSet = 0;
					uint32_t totalShared = 0;
					uint32_t totalTarsWithDataForBoth = 0;
					uint32_t totalTarsWithAtLeastOneHapShared = 0;
					std::vector<double> jacardsByTarget;
					for(const auto tpos : iter::range(tarNamesVec.size())){
						uint8_t tarRes = targetsEncodeBySamp[pair.col_][tpos] + targetsEncodeBySamp[pair.row_][tpos];
						//should be 2 if both samples have this target
						if(2 == tarRes){
							++totalTarsWithDataForBoth;
							uint32_t totalSetForTar = 0;
							uint32_t totalSharedForTar = 0;
							for(const auto & hapPos : iter::range(numberOfHapsPerTarget[tpos])){
								//position in the encoded vector should be the target start ranged over the possible haplotypes for that target
								uint8_t res = hapsEncodeBySamp[pair.col_][tarStart[tpos] + hapPos] + hapsEncodeBySamp[pair.row_][tarStart[tpos] + hapPos];
								//if res is 2 then haps are shared
								if(2 == res) {
									++totalShared;
									++totalSharedForTar;
								}
								//if results is either 1 or 2 then at least one of them has this hap
								if(res > 0) {
									++totalSet;
									++totalSetForTar;
								}
							}
							if(totalSharedForTar > 0){
								//at least one haplotype is shared between samps
								++totalTarsWithAtLeastOneHapShared;
							}
							jacardsByTarget.emplace_back(totalSharedForTar/static_cast<double>(totalSetForTar));
						}
					}
					byHapsTarShared[pair.row_][pair.col_] = totalShared/static_cast<double>(totalSet);
					byHapsTarShared[pair.col_][pair.row_] = totalShared/static_cast<double>(totalSet);
					avgJacard[pair.row_][pair.col_] = vectorMean(jacardsByTarget);
					avgJacard[pair.col_][pair.row_] = vectorMean(jacardsByTarget);
					byTarget[pair.row_][pair.col_] = totalTarsWithAtLeastOneHapShared/static_cast<double>(totalTarsWithDataForBoth);
					byTarget[pair.col_][pair.row_] = totalTarsWithAtLeastOneHapShared/static_cast<double>(totalTarsWithDataForBoth);
				}
				//for all
				{
					uint32_t totalSet = 0;
					uint32_t totalShared = 0;
					for(const auto pos : iter::range(totalHaps)){
						uint8_t res = hapsEncodeBySamp[pair.col_][pos] + hapsEncodeBySamp[pair.row_][pos];
						if(2 == res){
							++totalShared;
						}
						if(res > 0){
							++totalSet;
						}
					}
					byAllHaps[pair.row_][pair.col_] = totalShared/static_cast<double>(totalSet);
					byAllHaps[pair.col_][pair.row_] = totalShared/static_cast<double>(totalSet);
				}
			}
		}
	};


	njh::concurrent::runVoidFunctionThreaded(compSamps, numThreads);



	OutputStream outSampNamesOut(njh::files::make_path(setUp.pars_.directoryName_, "sampleNames.tab.txt"));
	OutputStream byTargetOut(njh::files::make_path(setUp.pars_.directoryName_, "percOfTarSharingAtLeastOneHap.tab.txt.gz"));
	OutputStream byHapOut(njh::files::make_path(setUp.pars_.directoryName_, "jacardByAllHap.tab.txt.gz"));
	OutputStream byHapTarSharedOut(njh::files::make_path(setUp.pars_.directoryName_, "jacardByHapsTarShared.tab.txt.gz"));

	OutputStream avgHapOut(njh::files::make_path(setUp.pars_.directoryName_, "avgJacardPerTarget.tab.txt.gz"));

	outSampNamesOut << njh::conToStr(sampNames, "\n") << std::endl;
	for(const auto & it : byTarget){
		byTargetOut << njh::conToStr(it, "\t") << std::endl;
	}
	for(const auto & ih : byHapsTarShared){
		byHapTarSharedOut << njh::conToStr(ih, "\t") << std::endl;
	}
	for(const auto & ih : byAllHaps){
		byHapOut << njh::conToStr(ih, "\t") << std::endl;
	}
	for(const auto & ih : avgJacard){
		avgHapOut << njh::conToStr(ih, "\t") << std::endl;
	}




	return 0;
}

}  //namespace njhseq



