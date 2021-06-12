/*
 * HapsEncodedMatrix.cpp
 *
 *  Created on: Jun 11, 2021
 *      Author: nick
 */

#include "HapsEncodedMatrix.hpp"


namespace njhseq {




void HapsEncodedMatrix::SetWithExternalPars::setDefaults(seqSetUp & setUp){
  setUp.setOption(tableFnp, "--tableFnp", "Table to read in (should be tab delimited)", true);

  setUp.setOption(sampleCol, "--sampleCol", "sampleCol", true);
  setUp.setOption(targetNameCol, "--targetNameCol", "targetNameCol", true);
  setUp.setOption(popIDCol, "--popIDCol", "popIDCol", true);
  setUp.setOption(relAbundCol, "--relAbundCol", "relAbundCol", true);

  setUp.setOption(selectSamples, "--selectSamples", "Only analyze these select samples");
  setUp.setOption(selectTargets, "--selectTargets", "Only analyze these select targets");
	setUp.setOption(numThreads, "--numThreads", "number of cpus to use");
	setUp.setOption(majorOnly, "--calcMajorHapOnly", "calculate differences by major haplotype only");
}


HapsEncodedMatrix::HapsEncodedMatrix(const SetWithExternalPars & pars): pars_(pars){
	TableReader hapTab(TableIOOpts(InOptions(pars.tableFnp), "\t", true));
	hapTab.header_.checkForColumnsThrow(VecStr{pars.sampleCol, pars.targetNameCol, pars.popIDCol, pars.relAbundCol}, __PRETTY_FUNCTION__);

	//set all info
	{
		//read in first to gather information on the table
		VecStr row;
		while(hapTab.getNextRow(row)){

			auto samp = row[hapTab.header_.getColPos(pars.sampleCol)];
			auto tar = row[hapTab.header_.getColPos(pars.targetNameCol)];
			if(!pars.selectSamples.empty() && !njh::in(samp, pars.selectSamples)){
				continue;
			}
			if(!pars.selectTargets.empty() && !njh::in(tar, pars.selectTargets)){
				continue;
			}
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
			if(!pars.selectSamples.empty() && !njh::in(samp, pars.selectSamples)){
				continue;
			}
			if(!pars.selectTargets.empty() && !njh::in(tar, pars.selectTargets)){
				continue;
			}
			auto hapName = row[reReadHapTab.header_.getColPos(pars.popIDCol)];
			//auto rBund = njh::StrToNumConverter::stoToNum<double>(row[reReadHapTab.header_.getColPos(pars.relAbundCol)]); //doing nothing right now with this

			encodeSampTarHap(samp, tar, hapName);
		}
	}
	if(pars.majorOnly){
		std::vector<std::vector<double>> hapsEncodeBySampRelAbund = std::vector<std::vector<double>> (sampNames_.size());
		for(const auto & samp : sampNames_){
			hapsEncodeBySampRelAbund[sampNamesKey_[samp]] = std::vector<double>(totalHaps_, 0);
		}
		TableReader reReadHapTab(TableIOOpts(InOptions(pars.tableFnp), "\t", true));
		VecStr row;
		while(reReadHapTab.getNextRow(row)){
			auto samp = row[reReadHapTab.header_.getColPos(pars.sampleCol)];
			auto tar = row[reReadHapTab.header_.getColPos(pars.targetNameCol)];
			if(!pars.selectSamples.empty() && !njh::in(samp, pars.selectSamples)){
				continue;
			}
			if(!pars.selectTargets.empty() && !njh::in(tar, pars.selectTargets)){
				continue;
			}
			auto hapName = row[reReadHapTab.header_.getColPos(pars.popIDCol)];
			auto rBund = njh::StrToNumConverter::stoToNum<double>(row[reReadHapTab.header_.getColPos(pars.relAbundCol)]);
			auto tKey = tarNameKey_[tar];
			auto hKey = hapNamesKey_[tar][hapName];
//				if(rBund > 0 && rBund < 1){
//					std::cout << rBund << std::endl;
//				}
			hapsEncodeBySampRelAbund[sampNamesKey_[samp]][tarStart_[tKey] + hKey] = rBund;
		}
//			std::cout << __FILE__ << " " << __LINE__ << std::endl;
		//now determine the major hap per sample per target
		for(const auto sampIndex : iter::range(sampNamesVec_.size())){
			for(const auto targetIndex : iter::range(tarNamesVec_.size())){
				//check if sample has target
				if(targetsEncodeBySamp_[sampIndex][targetIndex] > 0){
					auto maxAbundance = std::numeric_limits<double>::min();
					uint32_t bestIndex = std::numeric_limits<uint32_t>::max();
					//determine the highest rel abundance
					for(const auto hapIndex : iter::range(tarStart_[targetIndex], tarStart_[targetIndex] + numberOfHapsPerTarget_[targetIndex]) ){
//							if(hapsEncodeBySampRelAbund[sampIndex][hapIndex] > 0 && hapsEncodeBySampRelAbund[sampIndex][hapIndex] < 1){
//								std::cout << hapsEncodeBySampRelAbund[sampIndex][hapIndex] << std::endl;
//							}
						if(hapsEncodeBySampRelAbund[sampIndex][hapIndex] > maxAbundance){
							maxAbundance = hapsEncodeBySampRelAbund[sampIndex][hapIndex];
							bestIndex = hapIndex;
						}
					}
					//set all other indexes for the sample target hap presence to be 0
					for(const auto hapIndex : iter::range(tarStart_[targetIndex], tarStart_[targetIndex] + numberOfHapsPerTarget_[targetIndex]) ){
						if(hapIndex != bestIndex){
							hapsEncodeBySamp_[sampIndex][hapIndex] = 0;
						}
					}
				}
			}
		}
	}
}

HapsEncodedMatrix::SampTarHap::SampTarHap(const std::string &samp, const std::string &tar,
		const std::string &hap) :
		samp_(samp), tar_(tar), hap_(hap) {

}


void HapsEncodedMatrix::addSampTarHapForEncoding(const std::string &samp, const std::string &tar,
		const std::string &hap){
	sampNames_.emplace(samp);
	tarNames_.emplace(tar);
	++hapNamesForTars_[tar][hap];
}

void HapsEncodedMatrix::addSampTarHapForEncoding(const SampTarHap & adding){
	addSampTarHapForEncoding(adding.samp_, adding.tar_, adding.hap_);
}

void HapsEncodedMatrix::encodeSampTarHap(const std::string &samp, const std::string &tar,
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

void HapsEncodedMatrix::encodeSampTarHap(const SampTarHap & adding){
	encodeSampTarHap(adding.samp_, adding.tar_, adding.hap_);
}


void HapsEncodedMatrix::calcHapProbs(){
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

table HapsEncodedMatrix::getNumberTargetsPerSample() const {
	table ret(VecStr { "sample", "targetCount" });
	for (const auto row : iter::range(targetsEncodeBySamp_.size())) {
		ret.addRow(sampNamesVec_[row], vectorSum(targetsEncodeBySamp_[row]));
	}
	return ret;
}

void HapsEncodedMatrix::addMeta(const bfs::path & metaFnp){
	if(encodeKeysSet_){
		meta_ = std::make_shared<MultipleGroupMetaData>(metaFnp, std::set<std::string>(sampNames_.begin(), sampNames_.end()));
	}else{
		meta_ = std::make_shared<MultipleGroupMetaData>(metaFnp);
	}
}

void HapsEncodedMatrix::resetEncoding(){
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
void HapsEncodedMatrix::setEncodeKeys(){
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
		std::vector<std::string> hapNamesVec = getVectorOfMapKeys(hapName.second);
		//sort by hap sample abundance
		njh::sort(hapNamesVec, [&hapName](const std::string & hap1, const std::string & hap2){
			return hapName.second.at(hap1) > hapName.second.at(hap2);
		});
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
		for(const auto tarPos : iter::range(tarNamesVec_.size())){
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

HapsEncodedMatrix::IndexResults::IndexResults(const uint64_t numOfSamps){

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



HapsEncodedMatrix::IndexResults HapsEncodedMatrix::genIndexMeasures(bool verbose) const{
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
			for(const auto pairPos : iter::range(pairVec.pairs_.size())){
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
							for(const auto hapPos : iter::range(numberOfHapsPerTarget_[tpos])){
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
//						if(("8034209115" == sampNamesVec_[pair.row_] && "8025874502" == sampNamesVec_[pair.col_]) ||
//							 ("8025874502" == sampNamesVec_[pair.row_] && "8034209115" == sampNamesVec_[pair.col_]) ){
//							std::cout << sampNamesVec_[pair.row_] << ":" << sampNamesVec_[pair.col_] << std::endl;
//							std::cout << "pair.row_: " << pair.row_ << " " << sampNamesVec_[pair.row_] << std::endl;
//							std::cout << "pair.col_: " << pair.col_ << " " << sampNamesVec_[pair.col_]<< std::endl;
//							std::cout << "totalShared: " << totalShared << std::endl;
//							std::cout << "totalSet: " << totalSet << std::endl;
//							std::cout << "totalShared/static_cast<double>(totalSet): " << totalShared/static_cast<double>(totalSet) << std::endl;
//						}
					ret.byAllHaps[pair.row_][pair.col_] = totalShared/static_cast<double>(totalSet);
					ret.byAllHaps[pair.col_][pair.row_] = totalShared/static_cast<double>(totalSet);
				}
			}
		}
	};


	njh::concurrent::runVoidFunctionThreaded(compSamps, pars_.numThreads);

	return ret;
}


}  // namespace njhseq
