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
#include "elucidator/PopulationGenetics.h"



namespace njhseq {

//s_Sample s_studyid p_name                 h_popUID                  c_relAbund

//int genExpRunner::doPairwiseComparisonOnHapsSharing(const njh::progutils::CmdArgs & inputCommands){
//
//	bfs::path tableFnp = "";
//
//	std::string sampleCol = "";
//	std::string targetNameCol = "";
//	std::string popIDCol = "";
//	std::string relAbundCol = "";
//
//	uint32_t numThreads = 1;
//	seqSetUp setUp(inputCommands);
//	setUp.processVerbose();
//	setUp.processDebug();
//
//	setUp.setOption(numThreads, "--numThreads", "number of cpus to use");
//  setUp.setOption(tableFnp, "--tableFnp", "Table to read in (should be tab delimited)", true);
//
//  setUp.setOption(sampleCol, "--sampleCol", "sampleCol", true);
//  setUp.setOption(targetNameCol, "--targetNameCol", "targetNameCol", true);
//  setUp.setOption(popIDCol, "--popIDCol", "popIDCol", true);
//  setUp.setOption(relAbundCol, "--relAbundCol", "relAbundCol", true);
//
//  setUp.processDirectoryOutputName(bfs::path(bfs::basename(tableFnp)).string() + "_doPairwiseComparisonOnHapsSharing_TODAY", true);
//	setUp.finishSetUp(std::cout);
//
//	setUp.startARunLog(setUp.pars_.directoryName_);
//
//	table hapTab(tableFnp, "\t", true);
//	hapTab.checkForColumnsThrow(VecStr{sampleCol, targetNameCol, popIDCol, relAbundCol}, __PRETTY_FUNCTION__);
//
//	struct SampPopulation{
//
//		struct Tar{
//			Tar(){
//
//			}
//			Tar(const std::unordered_map<uint32_t, double> & haps): haps_(haps){
//
//			}
//			std::unordered_map<uint32_t, double> haps_;
//
//			uint32_t hapCount()const{
//				return haps_.size();
//			}
//
//
//
//			uint32_t hapsShared(const std::unordered_map<uint32_t, double> & otherHaps) const{
//				uint32_t ret = 0;
//				for(const auto & otherHap : otherHaps){
//					if(njh::in(otherHap.first, haps_)){
//						++ret;
//					}
//				}
//				return ret;
//			}
//		};
//
//		struct Samp{
//			Samp(){
//
//			}
//			Samp(const std::unordered_map<uint32_t, Tar> & targetsWithHaps): targetsWithHaps_(targetsWithHaps){
//
//			}
//			std::unordered_map<uint32_t, Tar> targetsWithHaps_;
//
//			struct SampComp {
//				uint32_t targetsShared_{0};
//				uint32_t targetsPossibleToShare_{0};
//				uint32_t hapsShared_{0};
//				uint32_t hapsPossibleToShareInSharedTargets_{0};
//
//				std::vector<double> jacardsByTargets_;
//
//				double percTargetsWithAtLeastOneHap() const {
//					return
//							targetsPossibleToShare_ > 0 ?
//									targetsShared_
//											/ static_cast<double>(targetsPossibleToShare_) :
//									0.0;
//				}
//				double jacardByHaps() const {
//					return
//							hapsShared_ > 0 ?
//									hapsPossibleToShareInSharedTargets_
//											/ static_cast<double>(hapsPossibleToShareInSharedTargets_) :
//									0.0;
//				}
//
//				double avgJacard() const {
//					return vectorMean(jacardsByTargets_);
//				}
//			};
//			SampComp compToOtherSamp(
//					const std::unordered_map<uint32_t, Tar> &otherTargets) const {
//				SampComp ret;
//				for (const auto &other : otherTargets) {
////					std::cout << other.first << std::endl;
//					if (njh::in(other.first, targetsWithHaps_)) {
//						++ret.targetsPossibleToShare_;
//						auto hapsShared = targetsWithHaps_.at(other.first).hapsShared(
//								other.second.haps_);
//						ret.hapsShared_ += hapsShared;
//						if (hapsShared > 0) {
//							++ret.targetsShared_;
//						}
//						std::unordered_set<uint32_t> haps;
//						std::vector<uint32_t> currentHaps = getVectorOfMapKeys(targetsWithHaps_);
//						haps.insert(currentHaps.begin(), currentHaps.end());
//						std::vector<uint32_t> currentOtherHaps = getVectorOfMapKeys(otherTargets);
//						haps.insert(currentOtherHaps.begin(), currentOtherHaps.end());
//						ret.hapsPossibleToShareInSharedTargets_ = haps.size();
//						ret.jacardsByTargets_.emplace_back(static_cast<double>(hapsShared)/haps.size());
//					}
//				}
//				return ret;
//			}
//			SampComp compToOtherSamp(const Samp &otherSamp) const {
//				return compToOtherSamp(otherSamp.targetsWithHaps_);
//			}
//
//		};
//		VecStr getSampNames()const{
//			return getVectorOfMapKeys(samples_);
//		}
//
//		std::unordered_map<std::string, Samp> samples_;
//	};
//
//
//
//	SampPopulation pop;
//
//	std::unordered_set<std::string> tarNames;
//	std::unordered_map<std::string, std::unordered_set<std::string>> hapNamesForTars;
//
//	std::unordered_map<std::string, uint32_t> tarNameKey;
//	std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>> hapNamesKey;
//
//	uint64_t totalHaps{0};
//	{
//		for(const auto & row : hapTab){
//			auto tar = row[hapTab.getColPos(targetNameCol)];
//			auto hapName = row[hapTab.getColPos(popIDCol)];
//			tarNames.emplace(tar);
//			hapNamesForTars[tar].emplace(hapName);
//		}
//		std::vector<std::string> tarNamesVec{tarNames.begin(), tarNames.end()};
//		for(const auto pos : iter::range(tarNames.size())){
//			tarNameKey[tarNamesVec[pos]] = pos;
//		}
//
//		for(const auto & hapName : hapNamesForTars){
//			std::vector<std::string> hapNamesVec{hapName.second.begin(), hapName.second.end()};
//			for(const auto pos : iter::range(hapNamesVec.size())){
//				hapNamesKey[hapName.first][hapNamesVec[pos]] = pos;
//				++totalHaps;
//			}
//		}
//	}
//
//	if(setUp.pars_.verbose_){
//		std::cout << totalHaps << std::endl;
//	}
//
//	for(const auto & row : hapTab){
//		auto samp = row[hapTab.getColPos(sampleCol)];
//		auto tar = row[hapTab.getColPos(targetNameCol)];
//		auto hapName = row[hapTab.getColPos(popIDCol)];
//		auto rBund = row[hapTab.getColPos(relAbundCol)];
//		pop.samples_[samp].targetsWithHaps_[tarNameKey[tar]].haps_[hapNamesKey[tar][hapName]] += njh::StrToNumConverter::stoToNum<double>(rBund);
//	}
//
//
//
//
//	auto sampNames = pop.getSampNames();
//
//
//	std::vector<std::vector<double>> byTarget{sampNames.size(), std::vector<double>(sampNames.size())};
//	std::vector<std::vector<double>> byHap{sampNames.size(), std::vector<double>(sampNames.size())};
//	std::vector<std::vector<double>> avgJacard{sampNames.size(), std::vector<double>(sampNames.size())};
//
//	for(uint32_t pos = 0; pos < sampNames.size(); ++pos){
//		byTarget[pos][pos] = 1;
//		byHap[pos][pos] = 1;
//		avgJacard[pos][pos] = 1;
//	}
//	PairwisePairFactory pFactor(sampNames.size());
//
//
//	njh::ProgressBar progpar(pFactor.totalCompares_);
//
//	std::function<void()> compSamps = [&pFactor,&byTarget,&byHap,&pop,&sampNames,&avgJacard](){
//		PairwisePairFactory::PairwisePairVec pairVec;
//
//		while(pFactor.setNextPairs(pairVec, 100)){
//			for(const auto & pairPos : iter::range(pairVec.pairs_.size())){
//				const auto & pair = pairVec.pairs_[pairPos];
////				if(setUp.pars_.verbose_){
////					progpar.outputProgAdd(std::cout, 1, true);
////					std::cout << sampNames[pair.row_] << " vs " << sampNames[pair.col_] << std::endl;
////				}
//				auto compRes = pop.samples_[sampNames[pair.row_]].compToOtherSamp(pop.samples_[sampNames[pair.col_]]);
//				byTarget[pair.row_][pair.col_] = compRes.percTargetsWithAtLeastOneHap();
//				byHap[pair.row_][pair.col_] = compRes.jacardByHaps();
//				avgJacard[pair.row_][pair.col_] = compRes.avgJacard();
//
//				byTarget[pair.col_][pair.row_] = compRes.percTargetsWithAtLeastOneHap();
//				byHap[pair.col_][pair.row_] = compRes.jacardByHaps();
//				avgJacard[pair.col_][pair.row_] = compRes.avgJacard();
//			}
//		}
//	};
//
//
//	njh::concurrent::runVoidFunctionThreaded(compSamps, numThreads);
//
//
//
//	OutputStream outSampNamesOut(njh::files::make_path(setUp.pars_.directoryName_, "sampleNames.tab.txt"));
//	OutputStream byTargetOut(njh::files::make_path(setUp.pars_.directoryName_, "percOfTarSharingAtLeastOneHap.tab.txt.gz"));
//	OutputStream byHapOut(njh::files::make_path(setUp.pars_.directoryName_, "jacardByAllHap.tab.txt.gz"));
//	OutputStream avgHapOut(njh::files::make_path(setUp.pars_.directoryName_, "avgJacardPerTarget.tab.txt.gz"));
//
//	outSampNamesOut << njh::conToStr(sampNames, "\n") << std::endl;
//	for(const auto & it : byTarget){
//		byTargetOut << njh::conToStr(it, "\t") << std::endl;
//	}
//	for(const auto & ih : byHap){
//		byHapOut << njh::conToStr(ih, "\t") << std::endl;
//	}
//	for(const auto & ih : avgJacard){
//		avgHapOut << njh::conToStr(ih, "\t") << std::endl;
//	}
//
//
//	return 0;
//}



class HapsEncodedMatrix{

public:


	struct SetWithExternalPars{
		bfs::path tableFnp = "";

		std::string sampleCol = "";
		std::string targetNameCol = "";
		std::string popIDCol = "";
		std::string relAbundCol = "";

		void setDefaults(seqSetUp & setUp){
		  setUp.setOption(tableFnp, "--tableFnp", "Table to read in (should be tab delimited)", true);

		  setUp.setOption(sampleCol, "--sampleCol", "sampleCol", true);
		  setUp.setOption(targetNameCol, "--targetNameCol", "targetNameCol", true);
		  setUp.setOption(popIDCol, "--popIDCol", "popIDCol", true);
		  setUp.setOption(relAbundCol, "--relAbundCol", "relAbundCol", true);
		}
	};

	HapsEncodedMatrix(const SetWithExternalPars & pars){
		TableReader hapTab(TableIOOpts(InOptions(pars.tableFnp), "\t", true));
		hapTab.header_.checkForColumnsThrow(VecStr{pars.sampleCol, pars.targetNameCol, pars.popIDCol, pars.relAbundCol}, __PRETTY_FUNCTION__);

		//set all info
		{
			//read in first to gather information on the table
			VecStr row;
			while(hapTab.getNextRow(row)){
				auto samp = row[hapTab.header_.getColPos(pars.sampleCol)];
				auto tar = row[hapTab.header_.getColPos(pars.targetNameCol)];
				auto hapName = row[hapTab.header_.getColPos(pars.popIDCol)];
				addSampTarHapForEncoding(samp, tar, hapName);
			}
		}
		//set encoding
		setEncodeKeys();
		//re-read and encode
		{
			TableReader reReadHapTab(TableIOOpts(InOptions(pars.tableFnp), "\t", true));
			VecStr row;
			while(reReadHapTab.getNextRow(row)){
				auto samp = row[reReadHapTab.header_.getColPos(pars.sampleCol)];
				auto tar = row[reReadHapTab.header_.getColPos(pars.targetNameCol)];
				auto hapName = row[reReadHapTab.header_.getColPos(pars.popIDCol)];
				auto rBund = row[reReadHapTab.header_.getColPos(pars.relAbundCol)]; //doing nothing right now with this

				encodeSampTarHap(samp, tar, hapName);
			}
		}
	}


	std::unordered_set<std::string> sampNames_; //! used during the encoding phase
	std::vector<std::string> sampNamesVec_; //! sample names, index is the key to the sample
	std::unordered_map<std::string, uint32_t> sampNamesKey_;//! key to sample name to sample index

	std::unordered_set<std::string> tarNames_; //! used during the encoding phase
	std::vector<std::string> tarNamesVec_; //! target names, index is the key to the targets
	std::unordered_map<std::string, uint32_t> tarNameKey_; //! key to target name to target index

	std::unordered_map<std::string, std::unordered_set<std::string>> hapNamesForTars_; //! all haplotype names for each target
	std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>> hapNamesKey_; //! key for each haplotype within each target, will range from 0 to haplotype count for target


	std::vector<uint32_t> numberOfHapsPerTarget_; //! number of haplotypes per target, index is the target as denoted from tarNameKey_
	std::vector<uint32_t> tarStart_; //! the starting position of target in targetsEncodeBySamp_ columns

	std::vector<std::vector<uint8_t>> hapsEncodeBySamp_; //! each row is a sample, each column is hap, 0 for not present, 1 for present
	std::vector<std::vector<uint8_t>> targetsEncodeBySamp_; //! each row is a sample, each column is a target, 0 if sample has no data for target, 1 for has data

	std::vector<double> hapsProbs_;

	uint64_t totalHaps_{0};

	bool encodeKeysSet_{false};

	std::shared_ptr<MultipleGroupMetaData> meta_;


	struct SampTarHap {
		SampTarHap(const std::string &samp, const std::string &tar,
				const std::string &hap) :
				samp_(samp), tar_(tar), hap_(hap) {

		}
		std::string samp_;
		std::string tar_;
		std::string hap_;

	};
	void addSampTarHapForEncoding(const std::string &samp, const std::string &tar,
			const std::string &hap){
		sampNames_.emplace(samp);
		tarNames_.emplace(tar);
		hapNamesForTars_[tar].emplace(hap);
	}

	void addSampTarHapForEncoding(const SampTarHap & adding){
		addSampTarHapForEncoding(adding.samp_, adding.tar_, adding.hap_);
	}

	void encodeSampTarHap(const std::string &samp, const std::string &tar,
			const std::string &hap){
		auto tKey = tarNameKey_[tar];
		auto hKey = hapNamesKey_[tar][hap];
//		std::cout << "sampNamesKey[samp]                        : " << sampNamesKey[samp] << std::endl;
//		std::cout << "tKey                                      : " << tKey << std::endl;
//		std::cout << "[tarStart[tKey]                           : " << tarStart[tKey] << std::endl;
//		std::cout << "hKey                                      : " << hKey << std::endl;
//		std::cout << "hKey                                      : " << hKey << std::endl;
//		std::cout << "hapsEncodeBySamp[sampNamesKey[samp]]size(): " << hapsEncodeBySamp[sampNamesKey[samp]].size() << std::endl;
		hapsEncodeBySamp_[sampNamesKey_[samp]][tarStart_[tKey] + hKey] = 1;
		targetsEncodeBySamp_[sampNamesKey_[samp]][tKey] = 1;
	}

	void encodeSampTarHap(const SampTarHap & adding){
		encodeSampTarHap(adding.samp_, adding.tar_, adding.hap_);
	}


	void calcHapProbs(){
		hapsProbs_ = std::vector<double>(totalHaps_, 0);
		for(uint32_t tarPos : iter::range(tarStart_.size())){
			std::vector<uint32_t> hapCounts(numberOfHapsPerTarget_[tarPos], 0);
			uint32_t totalHaps = 0;
			for(const auto hapPos : iter::range(numberOfHapsPerTarget_[tarPos])){
				for(const auto sampPos : iter::range(hapsEncodeBySamp_.size())){
					hapCounts[hapPos] += hapsEncodeBySamp_[sampPos][tarStart_[tarPos] + hapPos];
					totalHaps += hapsEncodeBySamp_[sampPos][tarStart_[tarPos] + hapPos];
				}
			}
			for(const auto hapPos : iter::range(numberOfHapsPerTarget_[tarPos])){
				hapsProbs_[tarStart_[tarPos] + hapPos] = hapCounts[hapPos]/static_cast<double>(totalHaps);
			}
		}
	}


	void addMeta(const bfs::path & metaFnp){
		if(encodeKeysSet_){
			meta_ = std::make_shared<MultipleGroupMetaData>(metaFnp, std::set<std::string>(sampNames_.begin(), sampNames_.end()));
		}else{
			meta_ = std::make_shared<MultipleGroupMetaData>(metaFnp);
		}
	}

	void resetEncoding(){
		sampNamesVec_.clear();
		sampNamesKey_.clear();
		tarNamesVec_.clear();
		tarNameKey_.clear();
		hapNamesKey_.clear();
		numberOfHapsPerTarget_.clear();
		tarStart_.clear();
		hapsEncodeBySamp_.clear();
		targetsEncodeBySamp_.clear();
		totalHaps_ = 0;

	}
	void setEncodeKeys(){
		if(encodeKeysSet_){
			//reset
			resetEncoding();
		}
		encodeKeysSet_ = true;
		sampNamesVec_ = std::vector<std::string> (sampNames_.begin(), sampNames_.end());
		tarNamesVec_ = std::vector<std::string> (tarNames_.begin(), tarNames_.end());
		njh::sort(sampNamesVec_);
		njh::sort(tarNamesVec_);

		for(const auto pos : iter::range(tarNamesVec_.size())){
			tarNameKey_[tarNamesVec_[pos]] = pos;
		}

		for(const auto pos : iter::range(sampNamesVec_.size())){
			sampNamesKey_[sampNamesVec_[pos]] = pos;
		}

		for(const auto & hapName : hapNamesForTars_){
			std::vector<std::string> hapNamesVec{hapName.second.begin(), hapName.second.end()};
			for(const auto pos : iter::range(hapNamesVec.size())){
				hapNamesKey_[hapName.first][hapNamesVec[pos]] = pos;
				++totalHaps_;
			}
		}
		numberOfHapsPerTarget_ = std::vector<uint32_t> (tarNamesVec_.size());
		for(const auto & tarHaps : hapNamesForTars_){
			numberOfHapsPerTarget_[tarNameKey_[tarHaps.first]] = tarHaps.second.size();
		}
		tarStart_ = std::vector<uint32_t> (tarNamesVec_.size());

		{
			uint32_t pos = 0;
			for(const auto & tarPos : iter::range(tarNamesVec_.size())){
				tarStart_[tarPos] = pos;
				pos += numberOfHapsPerTarget_[tarPos];
			}
		}
		hapsEncodeBySamp_ = std::vector<std::vector<uint8_t>> (sampNames_.size());
		for(const auto & samp : sampNames_){
			hapsEncodeBySamp_[sampNamesKey_[samp]] = std::vector<uint8_t>(totalHaps_, 0);
		}

		targetsEncodeBySamp_ = std::vector<std::vector<uint8_t>> (sampNames_.size());
		for(const auto & samp : sampNames_){
			targetsEncodeBySamp_[sampNamesKey_[samp]] = std::vector<uint8_t>(tarNamesVec_.size(), 0);
		}
	}

	struct IndexResults{
		IndexResults(const uint64_t numOfSamps){
			byAllHaps = std::vector<std::vector<double>>(numOfSamps, std::vector<double>(numOfSamps));
			byHapsTarShared = std::vector<std::vector<double>> (numOfSamps, std::vector<double>(numOfSamps));

			avgJacard = std::vector<std::vector<double>> (numOfSamps, std::vector<double>(numOfSamps));
			byTarget = std::vector<std::vector<double>> (numOfSamps, std::vector<double>(numOfSamps));

			byHapsTarSharedWeighted = std::vector<std::vector<double>> (numOfSamps, std::vector<double>(numOfSamps));
			avgJacardWeighted = std::vector<std::vector<double>> (numOfSamps, std::vector<double>(numOfSamps));


			for(uint32_t pos = 0; pos < numOfSamps; ++pos){
				//set diagonal
				byHapsTarShared[pos][pos] = 1;
				byAllHaps[pos][pos] = 1;

				avgJacard[pos][pos] = 1;
				byTarget[pos][pos] = 1;

				byHapsTarSharedWeighted[pos][pos] = 1;
				avgJacardWeighted[pos][pos] = 1;
			}

		}
		std::vector<std::vector<double>> byAllHaps; //! jacard index by all input haplotpes
		std::vector<std::vector<double>> byHapsTarShared;//! jacard index for haplotypes for targets where both samples have data
		std::vector<std::vector<double>> byHapsTarSharedWeighted;//! jacard index for haplotypes for targets where both samples have data

		std::vector<std::vector<double>> avgJacard; //! averaged jacard distance for shared targets
		std::vector<std::vector<double>> avgJacardWeighted; //! averaged jacard distance for shared targets weighted by
		std::vector<std::vector<double>> byTarget; //! fraction of targets that have at least shared haplotype between samples

	};

	IndexResults genIndexMeasures(uint32_t numThreads, bool verbose = false) const{
		IndexResults ret(sampNames_.size());

		PairwisePairFactory pFactor(sampNames_.size());


		njh::ProgressBar progpar(pFactor.totalCompares_);
		if(verbose){
			std::cout << "totalComps: " << pFactor.totalCompares_ << std::endl;
		}

		std::function<void()> compSamps = [&pFactor,&progpar,
																			 this,
																			 &ret,
																			 &verbose](){
			PairwisePairFactory::PairwisePairVec pairVec;
			while(pFactor.setNextPairs(pairVec, 100)){
				for(const auto & pairPos : iter::range(pairVec.pairs_.size())){
					const auto & pair = pairVec.pairs_[pairPos];
					if(verbose){
						progpar.outputProgAdd(std::cout, 1, true);
					}

					//for shared targets
					{
						uint32_t totalSet = 0;
						uint32_t totalShared = 0;
						double totalSetWeighted = 0;
						double totalSharedWeighted = 0;
						uint32_t totalTarsWithDataForBoth = 0;
						uint32_t totalTarsWithAtLeastOneHapShared = 0;
						std::vector<double> jacardsByTarget;
						std::vector<double> jacardsByTargetWeighted;

						for(const auto tpos : iter::range(tarNamesVec_.size())){
							uint8_t tarRes = targetsEncodeBySamp_[pair.col_][tpos] + targetsEncodeBySamp_[pair.row_][tpos];
							//should be 2 if both samples have this target
							if(2 == tarRes){
								++totalTarsWithDataForBoth;
								uint32_t totalSetForTar = 0;
								uint32_t totalSharedForTar = 0;
								double totalWeightedSetForTar = 0;
								double totalWeightedSharedForTar = 0;
								for(const auto & hapPos : iter::range(numberOfHapsPerTarget_[tpos])){
									//position in the encoded vector should be the target start ranged over the possible haplotypes for that target
									uint8_t res = hapsEncodeBySamp_[pair.col_][tarStart_[tpos] + hapPos] + hapsEncodeBySamp_[pair.row_][tarStart_[tpos] + hapPos];
									//if res is 2 then haps are shared
									if(2 == res) {
										++totalShared;
										++totalSharedForTar;
										totalWeightedSharedForTar += (1.01 - hapsProbs_[tarStart_[tpos] + hapPos]);
										totalSharedWeighted += (1.01 - hapsProbs_[tarStart_[tpos] + hapPos]);
									}
									//if results is either 1 or 2 then at least one of them has this hap
									if(res > 0) {
										++totalSet;
										++totalSetForTar;
										totalSetWeighted += (1.01 - hapsProbs_[tarStart_[tpos] + hapPos]);
										totalWeightedSetForTar += (1.01 - hapsProbs_[tarStart_[tpos] + hapPos]);
									}
								}
								if(totalSharedForTar > 0){
									//at least one haplotype is shared between samps
									++totalTarsWithAtLeastOneHapShared;
								}
								jacardsByTarget.emplace_back(totalSharedForTar/static_cast<double>(totalSetForTar));
								jacardsByTargetWeighted.emplace_back(totalWeightedSharedForTar/(totalWeightedSetForTar));
							}
						}
						ret.byHapsTarShared[pair.row_][pair.col_] = totalShared/static_cast<double>(totalSet);
						ret.byHapsTarShared[pair.col_][pair.row_] = totalShared/static_cast<double>(totalSet);
						ret.byHapsTarSharedWeighted[pair.row_][pair.col_] = totalSharedWeighted/(totalSetWeighted);
						ret.byHapsTarSharedWeighted[pair.col_][pair.row_] = totalSharedWeighted/(totalSetWeighted);
						auto meanJacard = vectorMean(jacardsByTarget);;
						auto meanJacardWeighted = vectorMean(jacardsByTargetWeighted);;
						ret.avgJacard[pair.row_][pair.col_] = meanJacard;
						ret.avgJacard[pair.col_][pair.row_] = meanJacard;
						ret.avgJacardWeighted[pair.row_][pair.col_] = meanJacardWeighted;
						ret.avgJacardWeighted[pair.col_][pair.row_] = meanJacardWeighted;
						ret.byTarget[pair.row_][pair.col_] = totalTarsWithAtLeastOneHapShared/static_cast<double>(totalTarsWithDataForBoth);
						ret.byTarget[pair.col_][pair.row_] = totalTarsWithAtLeastOneHapShared/static_cast<double>(totalTarsWithDataForBoth);
					}
					//for all
					{
						uint32_t totalSet = 0;
						uint32_t totalShared = 0;
						for(const auto pos : iter::range(totalHaps_)){
							uint8_t res = hapsEncodeBySamp_[pair.col_][pos] + hapsEncodeBySamp_[pair.row_][pos];
							if(2 == res){
								++totalShared;
							}
							if(res > 0){
								++totalSet;
							}
						}
						ret.byAllHaps[pair.row_][pair.col_] = totalShared/static_cast<double>(totalSet);
						ret.byAllHaps[pair.col_][pair.row_] = totalShared/static_cast<double>(totalSet);
					}
				}
			}
		};


		njh::concurrent::runVoidFunctionThreaded(compSamps, numThreads);

		return ret;
	}

};



int genExpRunner::doPairwiseComparisonOnHapsSharing(const njh::progutils::CmdArgs & inputCommands){

	bfs::path metaFnp;
	bfs::path tableFnp = "";

	std::string sampleCol = "";
	std::string targetNameCol = "";
	std::string popIDCol = "";
	std::string relAbundCol = "";

	uint32_t numThreads = 1;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(metaFnp, "--metaFnp", "Table of meta data for samples, needs a column named sample and each additonal column will be the meta data associated with that sample");

	setUp.setOption(numThreads, "--numThreads", "number of cpus to use");
  setUp.setOption(tableFnp, "--tableFnp", "Table to read in (should be tab delimited)", true);

  setUp.setOption(sampleCol, "--sampleCol", "sampleCol", true);
  setUp.setOption(targetNameCol, "--targetNameCol", "targetNameCol", true);
  setUp.setOption(popIDCol, "--popIDCol", "popIDCol", true);
  setUp.setOption(relAbundCol, "--relAbundCol", "relAbundCol", true);

  setUp.processDirectoryOutputName(bfs::path(bfs::basename(tableFnp)).string() + "_doPairwiseComparisonOnHapsSharing_TODAY", true);
	setUp.finishSetUp(std::cout);

	setUp.startARunLog(setUp.pars_.directoryName_);

	if("" != metaFnp){
		MultipleGroupMetaData meta(metaFnp);
	}

	TableReader hapTab(TableIOOpts(InOptions(tableFnp), "\t", true));
	hapTab.header_.checkForColumnsThrow(VecStr{sampleCol, targetNameCol, popIDCol, relAbundCol}, __PRETTY_FUNCTION__);
	std::unordered_set<std::string> sampNames;
	std::unordered_map<std::string, uint32_t> sampNamesKey;

	std::unordered_set<std::string> tarNames;
	std::unordered_map<std::string, uint32_t> tarNameKey;

	std::unordered_map<std::string, std::unordered_set<std::string>> hapNamesForTars;
	std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>> hapNamesKey;



	uint64_t totalHaps{0};
	{
		//read in first to gather information on the table
		VecStr row;
		while(hapTab.getNextRow(row)){
			auto samp = row[hapTab.header_.getColPos(sampleCol)];
			auto tar = row[hapTab.header_.getColPos(targetNameCol)];
			auto hapName = row[hapTab.header_.getColPos(popIDCol)];
			sampNames.emplace(samp);
			tarNames.emplace(tar);
			hapNamesForTars[tar].emplace(hapName);
		}
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
	if(setUp.pars_.verbose_){
		std::cout << "totalHaps: " << totalHaps << std::endl;
	}
//	std::cout << "tarStart: " << njh::conToStr(tarStart, ",") << std::endl;

	std::vector<std::vector<uint8_t>> hapsEncodeBySamp(sampNames.size());
	for(const auto & samp : sampNames){
		hapsEncodeBySamp[sampNamesKey[samp]] = std::vector<uint8_t>(totalHaps, 0);
	}

	std::vector<std::vector<uint8_t>> targetsEncodeBySamp(sampNames.size());
	for(const auto & samp : sampNames){
		targetsEncodeBySamp[sampNamesKey[samp]] = std::vector<uint8_t>(tarNamesVec.size(), 0);
	}
	{
		//re-read to encode the data
		//read in for encoding
		TableReader reReadHapTab(TableIOOpts(InOptions(tableFnp), "\t", true));
		VecStr row;
		while(reReadHapTab.getNextRow(row)){
			auto samp = row[reReadHapTab.header_.getColPos(sampleCol)];
			auto tar = row[reReadHapTab.header_.getColPos(targetNameCol)];
			auto hapName = row[reReadHapTab.header_.getColPos(popIDCol)];
			auto rBund = row[reReadHapTab.header_.getColPos(relAbundCol)]; //doing nothing right now with this

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

int genExpRunner::doPairwiseComparisonOnHapsSharingDev(const njh::progutils::CmdArgs & inputCommands){

	bfs::path metaFnp;
	VecStr metaFieldsToCalcPopDiffs{};
	HapsEncodedMatrix::SetWithExternalPars pars;

	uint32_t numThreads = 1;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(metaFnp, "--metaFnp", "Table of meta data for samples, needs a column named sample and each additonal column will be the meta data associated with that sample");
	setUp.setOption(metaFieldsToCalcPopDiffs, "--metaFieldsToCalcPopDiffs", "Meta Fields To Calc Pop Diffs");

	setUp.setOption(numThreads, "--numThreads", "number of cpus to use");
  pars.setDefaults(setUp);

  setUp.processDirectoryOutputName(bfs::path(bfs::basename(pars.tableFnp)).string() + "_doPairwiseComparisonOnHapsSharing_TODAY", true);
	setUp.finishSetUp(std::cout);

	setUp.startARunLog(setUp.pars_.directoryName_);


	setUp.timer_.setLapName("initial");
	setUp.timer_.startNewLap("encode haplotypes");
	HapsEncodedMatrix haps(pars);
	if("" != metaFnp){
		haps.addMeta(metaFnp);
		if(!metaFieldsToCalcPopDiffs.empty()){
			haps.meta_->checkForFieldsThrow(metaFieldsToCalcPopDiffs);
		}
	}
	setUp.timer_.startNewLap("get hap probabilities");
	haps.calcHapProbs();
	setUp.timer_.startNewLap("get index measures");
	auto indexRes = haps.genIndexMeasures(numThreads, setUp.pars_.verbose_);

	OutputStream outSampNamesOut(njh::files::make_path(setUp.pars_.directoryName_, "sampleNames.tab.txt"));
	OutputStream byTargetOut(njh::files::make_path(setUp.pars_.directoryName_, "percOfTarSharingAtLeastOneHap.tab.txt.gz"));
	OutputStream byHapOut(njh::files::make_path(setUp.pars_.directoryName_, "jacardByAllHap.tab.txt.gz"));
	OutputStream byHapTarSharedOut(njh::files::make_path(setUp.pars_.directoryName_, "jacardByHapsTarShared.tab.txt.gz"));
	OutputStream avgHapOut(njh::files::make_path(setUp.pars_.directoryName_, "avgJacardPerTarget.tab.txt.gz"));

	OutputStream byHapTarSharedWeightedOut(njh::files::make_path(setUp.pars_.directoryName_, "jacardByHapsTarSharedWeighted.tab.txt.gz"));
	OutputStream avgHapWeightedOut(njh::files::make_path(setUp.pars_.directoryName_, "avgJacardPerTargetWeighted.tab.txt.gz"));



	outSampNamesOut << njh::conToStr(haps.sampNames_, "\n") << std::endl;
	for(const auto & it : indexRes.byTarget){
		byTargetOut << njh::conToStr(it, "\t") << std::endl;
	}
	for(const auto & ih : indexRes.byHapsTarShared){
		byHapTarSharedOut << njh::conToStr(ih, "\t") << std::endl;
	}
	for(const auto & ih : indexRes.byAllHaps){
		byHapOut << njh::conToStr(ih, "\t") << std::endl;
	}
	for(const auto & ih : indexRes.avgJacard){
		avgHapOut << njh::conToStr(ih, "\t") << std::endl;
	}
	for(const auto & ih : indexRes.byHapsTarSharedWeighted){
		byHapTarSharedWeightedOut << njh::conToStr(ih, "\t") << std::endl;
	}
	for(const auto & ih : indexRes.avgJacardWeighted){
		avgHapWeightedOut << njh::conToStr(ih, "\t") << std::endl;
	}

	{
		setUp.timer_.startNewLap("get population pairwise measures");
		std::vector<uint32_t> tarKeys(haps.tarNamesVec_.size());
		njh::iota<uint32_t>(tarKeys, 0);
		njh::concurrent::LockableQueue<uint32_t> tarQueue(tarKeys);
		OutputStream diversityMeasuresOut(njh::files::make_path(setUp.pars_.directoryName_, "diversityMeasuresPerTarget.tab.txt"));
		diversityMeasuresOut << "loci\tsampCount\ttotalHaps\tuniqueHaps\the\tsinglets\tdoublets\teffectiveNumOfAlleles\tShannonEntropyE" << '\n';
		std::mutex divOutMut;
		std::function<void()> getTargetInfo = [&tarQueue,&haps,&diversityMeasuresOut,&divOutMut](){

			uint32_t tarKey = std::numeric_limits<uint32_t>::max();
			while(tarQueue.getVal(tarKey)){
				std::vector<PopGenCalculator::PopHapInfo> hapsForTarget;
				for(const auto tarpos : iter::range(haps.numberOfHapsPerTarget_[tarKey])){
					hapsForTarget.emplace_back(PopGenCalculator::PopHapInfo(tarpos, 0));
				}
				uint32_t sampleCount = 0;
				for(const auto sampPos : iter::range(haps.hapsEncodeBySamp_.size())){
					if(haps.targetsEncodeBySamp_[sampPos][tarKey] > 0){
						++sampleCount;
					}
					for(const auto tarpos : iter::range(haps.numberOfHapsPerTarget_[tarKey])){
						if(haps.hapsEncodeBySamp_[sampPos][haps.tarStart_[tarKey] + tarpos] > 0){
							hapsForTarget[tarpos].count_ +=1;
						}
					}
				}
				auto diversityForTar = PopGenCalculator::getGeneralMeasuresOfDiversity(hapsForTarget);
				auto totalHaps = PopGenCalculator::PopHapInfo::getTotalPopCount(hapsForTarget);
				{
					std::lock_guard<std::mutex> lock(divOutMut);
					diversityMeasuresOut << haps.tarNamesVec_[tarKey]
															<< "\t" << sampleCount
															<< "\t" << totalHaps
															<< "\t" << diversityForTar.alleleNumber_
															<< "\t" << diversityForTar.heterozygostiy_
															<< "\t" << diversityForTar.singlets_
															<< "\t" << diversityForTar.doublets_
															<< "\t" << diversityForTar.effectiveNumOfAlleles_
															<< "\t" << diversityForTar.ShannonEntropyE_
															<< '\n';
				}
			}
		};

		njh::concurrent::runVoidFunctionThreaded(getTargetInfo, numThreads);
	}



	if("" != metaFnp && !metaFieldsToCalcPopDiffs.empty()){
		setUp.timer_.startNewLap("get population pairwise measures");

		auto popMeasuresDir = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar{"popDiffMeasures"});


		std::vector<uint32_t> tarKeys(haps.tarNamesVec_.size());
		njh::iota<uint32_t>(tarKeys, 0);
		for(const auto & field : metaFieldsToCalcPopDiffs){
			njh::concurrent::LockableQueue<uint32_t> tarQueue(tarKeys);
			OutputStream diversityMeasuresOut(njh::files::make_path(popMeasuresDir, njh::pasteAsStr(field, "_diversityMeasures.tab.txt.gz")));
			diversityMeasuresOut << field << "\tloci\tsampCount\ttotalHaps\tuniqueHaps\the\tsinglets\tdoublets\teffectiveNumOfAlleles\tShannonEntropyE\tIn" << '\n';
			OutputStream diffMeasuresOut(njh::files::make_path(popMeasuresDir, njh::pasteAsStr(field, "_diffMeasures.tab.txt.gz")));
			OutputStream pairwiseDiffMeasuresOut(njh::files::make_path(popMeasuresDir, njh::pasteAsStr(field, "_pairwiseDiffMeasures.tab.txt.gz")));
			diffMeasuresOut << "meta" << "\t"<< "loci"
					<<"\t"<<"totalHaps"
					<<"\t"<<"uniqueHaps"
					<<"\t"<<"nsamples"
					<<"\t"<<"HsSample"
					<<"\t"<<"HsEst"
					<<"\t"<<"HtSample"
					<<"\t"<<"HtEst"
					<<"\t"<<"Gst"
					<<"\t"<<"GstEst"
					<<"\t"<<"JostD"
					<<"\t"<<"JostDEst"
					<<"\t"<<"ChaoA"
					<<"\t"<<"ChaoB"
					<<"\t"<<"JostDChaoEst"
					<<"\t"<<"In"<< std::endl;

			pairwiseDiffMeasuresOut << "loci"
					<< "\t" << field << "1"
					<< "\t" << "popMeta" << "1_totalHaps"
					<< "\t" << "popMeta" << "1_uniqueHaps"
					<< "\t" << "popMeta" << "1_samples"
					<< "\t" << field << "2"
					<< "\t" << "popMeta" << "2_totalHaps"
					<< "\t" << "popMeta" << "2_uniqueHaps"
					<< "\t" << "popMeta" << "2_samples"
					<< "\t" << "HsSample"
									<<"\t"<<"HsEst"
									<<"\t"<<"HtSample"
									<<"\t"<<"HtEst"
									<<"\t"<<"Gst"
									<<"\t"<<"GstEst"
									<<"\t"<<"JostD"
									<<"\t"<<"JostDEst"
									<<"\t"<<"ChaoA"
									<<"\t"<<"ChaoB"
									<<"\t"<<"JostDChaoEst"
									<<"\t"<< "In"<< std::endl;

			std::mutex divOutMut;
			std::vector<std::string> sampleToMeta;
			std::unordered_set<std::string> subFields;
			for(const auto sampPos : iter::range(haps.sampNamesVec_.size())){

				sampleToMeta.emplace_back(haps.meta_->groupData_[field]->getGroupForSample(haps.sampNamesVec_[sampPos]));
				subFields.emplace(sampleToMeta.back());
			}
			std::function<void()> getPopDiffMeasures = [&tarQueue,&haps, &diversityMeasuresOut,&diffMeasuresOut,&pairwiseDiffMeasuresOut,&divOutMut,&sampleToMeta,&subFields,&field](){

				uint32_t tarKey = std::numeric_limits<uint32_t>::max();
				while(tarQueue.getVal(tarKey)){

					std::unordered_map<std::string, std::vector<PopGenCalculator::PopHapInfo>> hapsForTargetPerPopulationRaw;
					for(const auto & subField : subFields){
						for(const auto tarpos : iter::range(haps.numberOfHapsPerTarget_[tarKey])){
							hapsForTargetPerPopulationRaw[subField].emplace_back(PopGenCalculator::PopHapInfo(tarpos, 0));
						}
					}
					std::unordered_map<std::string, uint32_t> sampleCount;
					for(const auto sampPos : iter::range(haps.hapsEncodeBySamp_.size())){
						if(haps.targetsEncodeBySamp_[sampPos][tarKey] > 0){
							++sampleCount[sampleToMeta[sampPos]];
						}
						for(const auto tarpos : iter::range(haps.numberOfHapsPerTarget_[tarKey])){
							if(haps.hapsEncodeBySamp_[sampPos][haps.tarStart_[tarKey] + tarpos] > 0){
								hapsForTargetPerPopulationRaw[sampleToMeta[sampPos]][tarpos].count_ +=1;
							}
						}
					}
					std::unordered_map<std::string, std::vector<PopGenCalculator::PopHapInfo>> hapsForTargetPerPopulation;
					for(const auto & pop : hapsForTargetPerPopulationRaw){
						for(const auto & hap : pop.second){
							if(hap.count_ > 0){
								hapsForTargetPerPopulation[pop.first].emplace_back(hap);
							}
						}
					}
					auto generalDiff = PopGenCalculator::getOverallPopDiff(hapsForTargetPerPopulation);

					std::unordered_map<std::string, std::unordered_map<std::string, PopGenCalculator::PopDifferentiationMeasures>> pairwiseDiffs;

					if(hapsForTargetPerPopulation.size() > 1){
						pairwiseDiffs = PopGenCalculator::getPairwisePopDiff(hapsForTargetPerPopulation);
					}
					std::unordered_map<std::string, PopGenCalculator::DiversityMeasures> divMeausresPerPop;
					for(const auto & hapsForPop : hapsForTargetPerPopulation){
						divMeausresPerPop[hapsForPop.first] = PopGenCalculator::getGeneralMeasuresOfDiversity(hapsForPop.second);
					}
					{
						std::lock_guard<std::mutex> lock(divOutMut);
						uint32_t grandTotalHaps = 0;
						uint32_t grandTotalSamples = 0;
						std::unordered_map<std::string, uint32_t> totalHapsPerPop;
						for(const auto & popDiv : divMeausresPerPop){
							auto totalHaps = PopGenCalculator::PopHapInfo::getTotalPopCount(hapsForTargetPerPopulation[popDiv.first]);
							totalHapsPerPop[popDiv.first] = totalHaps;
							grandTotalHaps += totalHaps;
							grandTotalSamples += sampleCount[popDiv.first];
							diversityMeasuresOut
							<< popDiv.first
							<< "\t" << haps.tarNamesVec_[tarKey]
																	<< "\t" << sampleCount[popDiv.first]
																	<< "\t" << totalHaps
																	<< "\t" << popDiv.second.alleleNumber_
																	<< "\t" << popDiv.second.heterozygostiy_
																	<< "\t" << popDiv.second.singlets_
																	<< "\t" << popDiv.second.doublets_
																	<< "\t" << popDiv.second.effectiveNumOfAlleles_
																	<< "\t" << popDiv.second.ShannonEntropyE_
																	<< "\t" << generalDiff.informativenessForAssignPerPopulation_[popDiv.first]
																	<< '\n';
						}
						diffMeasuresOut << field << "\t"
								<< haps.tarNamesVec_[tarKey]
								<<"\t"<< grandTotalHaps
								<<"\t"<< haps.numberOfHapsPerTarget_[tarKey]
								<<"\t"<< grandTotalSamples
								<<"\t"<< generalDiff.hsSample_
								<<"\t"<< generalDiff.hsEst_
								<<"\t"<< generalDiff.htSample_
								<<"\t"<< generalDiff.htEst_
								<<"\t"<< generalDiff.gst_
								<<"\t"<< generalDiff.gstEst_
								<<"\t"<< generalDiff.jostD_
								<<"\t"<< generalDiff.jostDEst_
								<<"\t"<< generalDiff.chaoA_
								<<"\t"<< generalDiff.chaoB_
								<<"\t"<< generalDiff.jostDChaoEst_
								<<"\t"<< generalDiff.informativenessForAssign_<< std::endl;
						if(hapsForTargetPerPopulation.size() > 1){
							auto keys = getVectorOfMapKeys(pairwiseDiffs);
							njh::sort(keys);
							for(const auto & key : keys){
								auto subKeys = getVectorOfMapKeys(pairwiseDiffs.at(key));
								njh::sort(subKeys);
								for(const auto & subKey : subKeys){
									pairwiseDiffMeasuresOut << haps.tarNamesVec_[tarKey]
											<< "\t" << key
											<< "\t" << totalHapsPerPop[key]
											<< "\t" << divMeausresPerPop[key].alleleNumber_
											<< "\t" << sampleCount[key]
											<< "\t" << subKey
											<< "\t" << totalHapsPerPop[subKey]
											<< "\t" << divMeausresPerPop[subKey].alleleNumber_
											<< "\t" << sampleCount[subKey]
																<<"\t"<< pairwiseDiffs.at(key).at(subKey).hsSample_
																<<"\t"<< pairwiseDiffs.at(key).at(subKey).hsEst_
																<<"\t"<< pairwiseDiffs.at(key).at(subKey).htSample_
																<<"\t"<< pairwiseDiffs.at(key).at(subKey).htEst_
																<<"\t"<< pairwiseDiffs.at(key).at(subKey).gst_
																<<"\t"<< pairwiseDiffs.at(key).at(subKey).gstEst_
																<<"\t"<< pairwiseDiffs.at(key).at(subKey).jostD_
																<<"\t"<< pairwiseDiffs.at(key).at(subKey).jostDEst_
																<<"\t"<< pairwiseDiffs.at(key).at(subKey).chaoA_
																<<"\t"<< pairwiseDiffs.at(key).at(subKey).chaoB_
																<<"\t"<< pairwiseDiffs.at(key).at(subKey).jostDChaoEst_
																<<"\t"<< pairwiseDiffs.at(key).at(subKey).informativenessForAssign_<< std::endl;
								}
							}
						}
					}
				}
			};

			njh::concurrent::runVoidFunctionThreaded(getPopDiffMeasures, numThreads);
		}
	}



	setUp.timer_.logLapTimes(setUp.rLog_.runLogFile_, true, 6, true);
	return 0;
}

}  //namespace njhseq



