/*
 * popGenExp_doPairwiseComparisonOnHapsSharing.cpp
 *
 *  Created on: Jun 11, 2021
 *      Author: nick
 */




#include "popGenExp.hpp"
#include "elucidator/objects/BioDataObject.h"
#include <njhseq/objects/counters/DNABaseCounter.hpp>
#include <njhseq/PopulationGenetics.h>


namespace njhseq {



int popGenExpRunner::doPairwiseComparisonOnHapsSharing(const njh::progutils::CmdArgs & inputCommands){
	bool writeOutDistMatrices = false;
	bfs::path metaFnp;
	VecStr metaFieldsToCalcPopDiffs{};
	HapsEncodedMatrix::SetWithExternalPars pars;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(writeOutDistMatrices, "--writeOutDistMatrices", "write Out Dist Matrices");
	setUp.setOption(metaFnp, "--metaFnp", "Table of meta data for samples, needs a column named sample and each additonal column will be the meta data associated with that sample");
	setUp.setOption(metaFieldsToCalcPopDiffs, "--metaFieldsToCalcPopDiffs", "Meta Fields To Calc Pop Diffs");
  pars.setDefaults(setUp);

  setUp.processDirectoryOutputName(bfs::path(bfs::basename(pars.tableFnp)).string() + "_doPairwiseComparisonOnHapsSharing_TODAY", true);
	setUp.finishSetUp(std::cout);

	setUp.startARunLog(setUp.pars_.directoryName_);


	setUp.timer_.setLapName("initial");
	setUp.timer_.startNewLap("encode haplotypes");
  HapsEncodedMatrix haps(pars);
  if ("" != metaFnp) {
    haps.addMeta(metaFnp);
    if (!metaFieldsToCalcPopDiffs.empty()) {
      haps.meta_->checkForFieldsThrow(metaFieldsToCalcPopDiffs);
    }
  } else if (!metaFieldsToCalcPopDiffs.empty()) {
    haps.addMetaWithInputTab(njh::vecToSet(metaFieldsToCalcPopDiffs));
  }
	setUp.timer_.startNewLap("get hap probabilities");
	haps.calcHapProbs();
	setUp.timer_.startNewLap("writing sample info");

	auto numTargetsPerSample = haps.getNumberTargetsPerSample();
	numTargetsPerSample.sortTable("sample", true);
	OutputStream outSamplesPerTarget(njh::files::make_path(setUp.pars_.directoryName_, "numTargetsPerSample.tab.txt"));
	numTargetsPerSample.outPutContents(outSamplesPerTarget, "\t");
	setUp.timer_.startNewLap("get index measures");
	auto indexRes = haps.genIndexMeasures(setUp.pars_.verbose_);
	if(writeOutDistMatrices){
		OutputStream outSampNamesOut(njh::files::make_path(setUp.pars_.directoryName_, "sampleNames.tab.txt"));
		OutputStream byTargetOut(njh::files::make_path(setUp.pars_.directoryName_, "percOfTarSharingAtLeastOneHap.tab.txt.gz"));
		OutputStream byHapOut(njh::files::make_path(setUp.pars_.directoryName_, "jacardByAllHap.tab.txt.gz"));
		OutputStream byHapTarSharedOut(njh::files::make_path(setUp.pars_.directoryName_, "jacardByHapsTarShared.tab.txt.gz"));
		OutputStream avgHapOut(njh::files::make_path(setUp.pars_.directoryName_, "avgJacardPerTarget.tab.txt.gz"));

		OutputStream byHapTarSharedWeightedOut(njh::files::make_path(setUp.pars_.directoryName_, "jacardByHapsTarSharedWeighted.tab.txt.gz"));
		OutputStream avgHapWeightedOut(njh::files::make_path(setUp.pars_.directoryName_, "avgJacardPerTargetWeighted.tab.txt.gz"));
		OutputStream targetsSharedBetweenSampsOut(njh::files::make_path(setUp.pars_.directoryName_, "targetsSharedBetweenSamps.tab.txt.gz"));



		outSampNamesOut << njh::conToStr(haps.sampNamesVec_, "\n") << std::endl;
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

		for(const auto & ih : indexRes.targetsShared){
			targetsSharedBetweenSampsOut << njh::conToStr(ih, "\t") << std::endl;
		}
	}

	{
		setUp.timer_.startNewLap("get population pairwise measures");
		std::vector<uint32_t> tarKeys(haps.tarNamesVec_.size());
		njh::iota<uint32_t>(tarKeys, 0);
		njh::concurrent::LockableQueue<uint32_t> tarQueue(tarKeys);
		OutputStream diversityMeasuresOut(njh::files::make_path(setUp.pars_.directoryName_, "diversityMeasuresPerTarget.tab.txt"));
		diversityMeasuresOut << "loci\tsampCount\ttotalHaps\tuniqueHaps\tSimpsonI\the\tExpP3\tExpP4\tExpP5\tsinglets\tdoublets\teffectiveNumOfAlleles\tShannonEntropyE" << '\n';
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
															<< "\t" << diversityForTar.simpsonIndex_
															<< "\t" << diversityForTar.heterozygostiy_
															<< "\t" << diversityForTar.ploidy3_.expectedCOIForPloidy_.at(3)
															<< "\t" << diversityForTar.ploidy4_.expectedCOIForPloidy_.at(4)
															<< "\t" << diversityForTar.ploidy5_.expectedCOIForPloidy_.at(5)
															<< "\t" << diversityForTar.singlets_
															<< "\t" << diversityForTar.doublets_
															<< "\t" << diversityForTar.effectiveNumOfAlleles_
															<< "\t" << diversityForTar.ShannonEntropyE_
															<< '\n';
				}
			}
		};

		njh::concurrent::runVoidFunctionThreaded(getTargetInfo, pars.numThreads);
	}



	if(!metaFieldsToCalcPopDiffs.empty()){
		setUp.timer_.startNewLap("get population pairwise measures");

		auto popMeasuresDir = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar{"popDiffMeasures"});


		std::vector<uint32_t> tarKeys(haps.tarNamesVec_.size());
		njh::iota<uint32_t>(tarKeys, 0);
		for(const auto & field : metaFieldsToCalcPopDiffs){
			njh::concurrent::LockableQueue<uint32_t> tarQueue(tarKeys);
			OutputStream diversityMeasuresOut(njh::files::make_path(popMeasuresDir, njh::pasteAsStr(field, "_diversityMeasures.tab.txt.gz")));
			diversityMeasuresOut << field << "\tloci\tsampCount\ttotalHaps\tuniqueHaps\tSimpsonI\the\tExpP3\tExpP4\tExpP5\tsinglets\tdoublets\teffectiveNumOfAlleles\tShannonEntropyE\tIn" << '\n';
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
					<< "\t" << "hapsOnlyIn_popMeta" << "1"
					<< "\t" << "hapsOnlyIn_popMeta" << "1CumFreq"
					<< "\t" << field << "2"
					<< "\t" << "popMeta" << "2_totalHaps"
					<< "\t" << "popMeta" << "2_uniqueHaps"
					<< "\t" << "popMeta" << "2_samples"
					<< "\t" << "hapsOnlyIn_popMeta" << "2"
					<< "\t" << "hapsOnlyIn_popMeta" << "2CumFreq"
					<< "\t" << "uniqHapsCombinedPops"
					<< "\t" << "uniqHapsSharedInPops"
					<< "\t" << "HsSample"
									<< "\t" << "HsEst"
									<< "\t" << "HtSample"
									<< "\t" << "HtEst"
									<< "\t" << "Gst"
									<< "\t" << "GstEst"
									<< "\t" << "JostD"
									<< "\t" << "JostDEst"
									<< "\t" << "ChaoA"
									<< "\t" << "ChaoB"
									<< "\t" << "JostDChaoEst"
									<< "\t" << "In"

									<< "\t" << "brayCurtisDissim"
									<< "\t" << "brayCurtisRelativeDissim"
									<< "\t" << "jaccardIndexDissim"
									<< "\t" << "sorensenDistance"
									<< "\t" << "RMSE"
									<< "\t" << "correlationDissim"
									<< "\t" << "matchingCoefficientDistance"
									<< "\t" << "plainAvalance"
									<< std::endl;

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

					PopGenCalculator::PopDifferentiationMeasures generalDiff;
					if(hapsForTargetPerPopulation.size() > 1){
						generalDiff = PopGenCalculator::getOverallPopDiff(hapsForTargetPerPopulation);
					}
					std::unordered_map<std::string, std::unordered_map<std::string, PopGenCalculator::PopDifferentiationMeasuresPairWise>> pairwiseDiffs;

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
																	<< "\t" << popDiv.second.simpsonIndex_
																	<< "\t" << popDiv.second.heterozygostiy_
																	<< "\t" << popDiv.second.ploidy3_.expectedCOIForPloidy_.at(3)
																	<< "\t" << popDiv.second.ploidy4_.expectedCOIForPloidy_.at(4)
																	<< "\t" << popDiv.second.ploidy5_.expectedCOIForPloidy_.at(5)
																	<< "\t" << popDiv.second.singlets_
																	<< "\t" << popDiv.second.doublets_
																	<< "\t" << popDiv.second.effectiveNumOfAlleles_
																	<< "\t" << popDiv.second.ShannonEntropyE_
																	<< "\t" << (hapsForTargetPerPopulation.size() > 1 ? generalDiff.informativenessForAssignPerPopulation_[popDiv.first] : 0)
																	<< '\n';
						}

						if(hapsForTargetPerPopulation.size() > 1){
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
											<< "\t" << pairwiseDiffs.at(key).at(subKey).uniqueHapsInPop1_
											<< "\t" << pairwiseDiffs.at(key).at(subKey).uniqueHapsInPop1CumFreq_
											<< "\t" << subKey
											<< "\t" << totalHapsPerPop[subKey]
											<< "\t" << divMeausresPerPop[subKey].alleleNumber_
											<< "\t" << sampleCount[subKey]
											<< "\t" << pairwiseDiffs.at(key).at(subKey).uniqueHapsInPop2_
											<< "\t" << pairwiseDiffs.at(key).at(subKey).uniqueHapsInPop2CumFreq_

											<< "\t" << pairwiseDiffs.at(key).at(subKey).uniqueHapsAll_
											<< "\t" << pairwiseDiffs.at(key).at(subKey).uniqueHapsShared_
																<<"\t"<< pairwiseDiffs.at(key).at(subKey).genDiffMeasures_.hsSample_
																<<"\t"<< pairwiseDiffs.at(key).at(subKey).genDiffMeasures_.hsEst_
																<<"\t"<< pairwiseDiffs.at(key).at(subKey).genDiffMeasures_.htSample_
																<<"\t"<< pairwiseDiffs.at(key).at(subKey).genDiffMeasures_.htEst_
																<<"\t"<< pairwiseDiffs.at(key).at(subKey).genDiffMeasures_.gst_
																<<"\t"<< pairwiseDiffs.at(key).at(subKey).genDiffMeasures_.gstEst_
																<<"\t"<< pairwiseDiffs.at(key).at(subKey).genDiffMeasures_.jostD_
																<<"\t"<< pairwiseDiffs.at(key).at(subKey).genDiffMeasures_.jostDEst_
																<<"\t"<< pairwiseDiffs.at(key).at(subKey).genDiffMeasures_.chaoA_
																<<"\t"<< pairwiseDiffs.at(key).at(subKey).genDiffMeasures_.chaoB_
																<<"\t"<< pairwiseDiffs.at(key).at(subKey).genDiffMeasures_.jostDChaoEst_
																<<"\t"<< pairwiseDiffs.at(key).at(subKey).genDiffMeasures_.informativenessForAssign_


																<< "\t" << pairwiseDiffs.at(key).at(subKey).brayCurtisDissim_
																<< "\t" << pairwiseDiffs.at(key).at(subKey).brayCurtisRelativeDissim_
																<< "\t" << pairwiseDiffs.at(key).at(subKey).jaccardIndexDissim_
																<< "\t" << pairwiseDiffs.at(key).at(subKey).sorensenDistance_
																<< "\t" << pairwiseDiffs.at(key).at(subKey).RMSE_
																<< "\t" << pairwiseDiffs.at(key).at(subKey).halfR_
																<< "\t" << pairwiseDiffs.at(key).at(subKey).matchingCoefficientDistance_
																<< "\t" << pairwiseDiffs.at(key).at(subKey).plainAvalance_

																<< std::endl;
								}
							}
						}
					}
				}
			};

			njh::concurrent::runVoidFunctionThreaded(getPopDiffMeasures, pars.numThreads);
		}
	}



	setUp.timer_.logLapTimes(setUp.rLog_.runLogFile_, true, 6, true);
	return 0;
}



} //namespace njhseq
