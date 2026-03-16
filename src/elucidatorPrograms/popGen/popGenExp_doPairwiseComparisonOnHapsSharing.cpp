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
#include <njhseq/objects/dataContainers/BasicPointMatrix.hpp>



namespace njhseq {




int popGenExpRunner::calc_pairwise_ccc_on_haps_sharing(const njh::progutils::CmdArgs & inputCommands){
	double minimumLociCoverageToKeepSamples = 0.90;
	HapsEncodedMatrix::SetWithExternalPars pars;
  uint32_t pairwise_factor_bin_size = 1000;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
  setUp.setOption(pairwise_factor_bin_size, "--pairwise_factor_bin_size", "pairwise_factor_bin_size");
	setUp.setOption(minimumLociCoverageToKeepSamples, "--minimumLociCoverageToKeepSamples", "minimum Loci Coverage To Keep Samples in post analysis steps, must have reads for at least this frction of the total loci");
  pars.setDefaults(setUp);

  setUp.processDirectoryOutputName(bfs::path(bfs::basename(pars.tableFnp)).string() + "_ccc_rmse_TODAY", true);
	setUp.finishSetUp(std::cout);

	setUp.startARunLog(setUp.pars_.directoryName_);


	setUp.timer_.setLapName("initial");
	setUp.timer_.startNewLap("encode haplotypes");
  HapsEncodedMatrix haps(pars);
	setUp.timer_.startNewLap("get hap probabilities");
	haps.calcHapProbs();
	setUp.timer_.startNewLap("add relative abundances");
  haps.add_relative_abundance();
  setUp.timer_.startNewLap("calc rmse and ccc");
	auto measures = haps.calc_ccc_rmse_measures(pairwise_factor_bin_size, setUp.pars_.verbose_);
  setUp.timer_.startNewLap("writing output matrices");
	OutputStream outSampNamesOut(njh::files::make_path(setUp.pars_.directoryName_, "sampleNames.tab.txt"));
	outSampNamesOut << njh::conToStr(haps.sampNamesVec_, "\n") << std::endl;
  {
	  auto ccc_out_fnp = njh::files::make_path(setUp.pars_.directoryName_, "ccc_on_targets_shared.tab.txt.gz");
	  OutputStream ccc_out(ccc_out_fnp);
	  for(const auto & ccc_row : measures.ccc){
	    ccc_out << njh::conToStr(ccc_row, "\t") << std::endl;
	  }
  }
  {
	  auto rmse_out_fnp = njh::files::make_path(setUp.pars_.directoryName_, "rmse_on_targets_shared.tab.txt.gz");
	  OutputStream rmse_out(rmse_out_fnp);
	  for(const auto & rmse_row : measures.rmse){
	    rmse_out << njh::conToStr(rmse_row, "\t") << std::endl;
	  }
  }

  {
	  auto targets_shared_out_fnp = njh::files::make_path(setUp.pars_.directoryName_, "targets_shared.tab.txt.gz");
	  OutputStream targets_shared_out(targets_shared_out_fnp);
	  for(const auto & targets_shared_row : measures.targets_shared){
	    targets_shared_out << njh::conToStr(targets_shared_row, "\t") << std::endl;
	  }
  }
  setUp.timer_.startNewLap("getting loci coverage info");

	std::unordered_map<std::string, double> lociCoveragePerSample = haps.getTargetCoveragePerSample();
	{
		table numTargetsPerSample = haps.getTableNumberTargetsPerSample(minimumLociCoverageToKeepSamples);
		OutputStream lociCoverageOut(njh::files::make_path(setUp.pars_.directoryName_, "loci_coverage_per_sample_info.tsv"));
		numTargetsPerSample.outPutContents(lociCoverageOut, "\t");
	}
	setUp.timer_.logLapTimes(setUp.rLog_.runLogFile_, true, 6, true);
	return 0;
}


int popGenExpRunner::cluster_samples_using_ccc_of_microhaps(const njh::progutils::CmdArgs & inputCommands){
	double minimumLociCoverageToKeepSamples = 0.90;
  double concordance_cut_off = 0.95;
  njhUndirWeightedGraph<double, std::vector<double>>::dbscanPars dbscanPars;
  // dbscanPars.eps_ = 0.50;
  dbscanPars.minEpNeighbors_ = 2;
	HapsEncodedMatrix::SetWithExternalPars pars;
  uint32_t pairwise_factor_bin_size = 1000;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
  setUp.setOption(pairwise_factor_bin_size, "--pairwise_factor_bin_size", "pairwise_factor_bin_size");
	setUp.setOption(minimumLociCoverageToKeepSamples, "--minimumLociCoverageToKeepSamples", "minimum Loci Coverage To Keep Samples in post analysis steps, must have reads for at least this frction of the total loci");
  setUp.setOption(concordance_cut_off, "--concordance_cut_off", "concordance cut off");
  dbscanPars.eps_ = 1 - concordance_cut_off;
  setUp.setOption(dbscanPars.minEpNeighbors_, "--min_group_size", "The minimum number of samples to group together");

  pars.setDefaults(setUp);

  setUp.processDirectoryOutputName(bfs::path(bfs::basename(pars.tableFnp)).string() + "_ccc_rmse_TODAY", true);
	setUp.finishSetUp(std::cout);

	setUp.startARunLog(setUp.pars_.directoryName_);


	setUp.timer_.setLapName("initial");
	setUp.timer_.startNewLap("encode haplotypes");
  HapsEncodedMatrix haps(pars);
	setUp.timer_.startNewLap("get hap probabilities");
	haps.calcHapProbs();
	setUp.timer_.startNewLap("add relative abundances");
  haps.add_relative_abundance();
  setUp.timer_.startNewLap("calc rmse and ccc");
	auto measures = haps.calc_ccc_rmse_measures(pairwise_factor_bin_size, setUp.pars_.verbose_);
  setUp.timer_.startNewLap("writing output matrices");
	OutputStream outSampNamesOut(njh::files::make_path(setUp.pars_.directoryName_, "sampleNames.tab.txt"));
	outSampNamesOut << njh::conToStr(haps.sampNamesVec_, "\n") << std::endl;
  {
	  auto ccc_out_fnp = njh::files::make_path(setUp.pars_.directoryName_, "ccc_on_targets_shared.tab.txt.gz");
	  OutputStream ccc_out(ccc_out_fnp);
	  for(const auto & ccc_row : measures.ccc){
	    ccc_out << njh::conToStr(ccc_row, "\t") << std::endl;
	  }
  }
  {
	  auto targets_shared_out_fnp = njh::files::make_path(setUp.pars_.directoryName_, "targets_shared.tab.txt.gz");
	  OutputStream targets_shared_out(targets_shared_out_fnp);
	  for(const auto & targets_shared_row : measures.targets_shared){
	    targets_shared_out << njh::conToStr(targets_shared_row, "\t") << std::endl;
	  }
  }
  setUp.timer_.startNewLap("getting loci coverage info");

	std::unordered_map<std::string, double> lociCoveragePerSample = haps.getTargetCoveragePerSample();
	{
		table numTargetsPerSample = haps.getTableNumberTargetsPerSample(minimumLociCoverageToKeepSamples);
		OutputStream lociCoverageOut(njh::files::make_path(setUp.pars_.directoryName_, "loci_coverage_per_sample_info.tsv"));
		numTargetsPerSample.outPutContents(lociCoverageOut, "\t");
	}
  setUp.timer_.setLapName("transforming matrix");
	{
	  //for the distance functions below to work, have to transform CCC so that the lower the better, CCC runs from -1 to 1, so below will transform it so it runs from 0 to 2 with 0 being CCC of 1, 1 being CCC 0, and 2 being CCC -2
	  PairwisePairFactory pairFactory(measures.ccc.size());
	  uint32_t pairBatchCount = 100000;
	  std::function<void()> transform_ccc =
    [&pairFactory,
      &pairBatchCount,
      &measures]() {
      PairwisePairFactory::PairwisePairVec pairs;
      while (pairFactory.setNextPairs(pairs, pairBatchCount)) {
        for (const auto & pair : pairs.pairs_) {
          measures.ccc[pair.row_][pair.col_] = -1 * (measures.ccc[pair.row_][pair.col_] - 1);
          measures.ccc[pair.col_][pair.row_] = measures.ccc[pair.row_][pair.col_];
        }
      }
    };
	  njh::concurrent::runVoidFunctionThreaded(transform_ccc, pars.numThreads);
    // fill the diagonal
	  for (uint32_t pos = 0; pos < measures.ccc.size(); ++pos) {
	    measures.ccc[pos][pos] = 0;
	  }
	}
  setUp.timer_.setLapName("building matrix");
  auto dist_graph = std::make_unique<njhUndirWeightedGraph<double, std::vector<double> > > ();
  for (const auto & pos : iter::range(measures.ccc.size())) {
    dist_graph->addNode(estd::to_string(pos), measures.ccc[pos]);
  }
	{
	  uint32_t belowEp = 0;
	  PairwisePairFactory pairFactory(measures.ccc.size());
	  uint32_t pairBatchCount = 100000;
	  std::mutex graphMut;
	  struct PairDist {
	    PairDist(const PairwisePairFactory::PairwisePair & pair, double dist) :
          pair_(pair), dist_(dist) {
	    }
	    PairwisePairFactory::PairwisePair pair_;
	    double dist_;
	  };

    std::function<void()> addToGraph =
        [&graphMut, &pairFactory,&pairBatchCount,&belowEp,
          &measures, &dbscanPars,
          &haps, &dist_graph,
          &lociCoveragePerSample, &minimumLociCoverageToKeepSamples]() {
      PairwisePairFactory::PairwisePairVec pairs;
      std::vector<PairDist> belowEps;
      while (pairFactory.setNextPairs(pairs, pairBatchCount)) {
        for (const auto &pair: pairs.pairs_) {
          if (lociCoveragePerSample[haps.sampNamesVec_[pair.row_]] < minimumLociCoverageToKeepSamples ||
              lociCoveragePerSample[haps.sampNamesVec_[pair.col_]] < minimumLociCoverageToKeepSamples) {
            continue;
          }
          auto dist = measures.ccc[pair.row_][pair.col_];
          if (dist < dbscanPars.eps_) {
            belowEps.emplace_back(PairDist{pair, dist});
          }
        }
      }
      if (!belowEps.empty()) {
        std::lock_guard<std::mutex> lock(graphMut);
        belowEp += belowEps.size();
        for (const auto &bEps: belowEps) {
          dist_graph->addEdge(estd::to_string(bEps.pair_.row_),
                              estd::to_string(bEps.pair_.col_),
                              bEps.dist_);
        }
      }
    };
	  njh::concurrent::runVoidFunctionThreaded(addToGraph, pars.numThreads);
	}

  setUp.timer_.startNewLap("dbscan");
  dist_graph->dbscan(dbscanPars);

  setUp.timer_.startNewLap("output");
  OutputStream outFile(OutOptions(njh::files::make_path(setUp.pars_.directoryName_, "clusters_by_ccc.tsv")));
  outFile << "sample\tgroup";
  outFile << std::endl;

  std::map<uint32_t, std::vector<uint32_t>> groupIndexes;
  std::vector<uint32_t> allGroupedIndices;
  for(const auto & n : iter::enumerate(dist_graph->nodes_)) {
    groupIndexes[n.element->group_].emplace_back(n.index);
    if (std::numeric_limits<uint32_t>::max() != n.element->group_) {
      allGroupedIndices.emplace_back(n.index);
    }
  }
  std::unordered_map<uint32_t, std::map<std::string, double> > groups_ccc_stats;
  for (const auto &group: groupIndexes) {
    std::vector<double> cccsWithinGroup; {
      PairwisePairFactory pfac(group.second.size());
      PairwisePairFactory::PairwisePair pair;
      while (pfac.setNextPair(pair)) {
        //have to transform back due to the previous transform
        cccsWithinGroup.emplace_back(measures.ccc[group.second[pair.col_]][group.second[pair.row_]] * -1 + 1);
      }
    }
    groups_ccc_stats[group.first] = getStatsOnVec(cccsWithinGroup);
    for (const auto &idx: group.second) {
      if (lociCoveragePerSample[haps.sampNamesVec_[idx]] < minimumLociCoverageToKeepSamples) {
        outFile << haps.sampNamesVec_[idx] << "\t" << "low_coverage_not_clustered";
      } else if (group.first == std::numeric_limits<uint32_t>::max()) {
        outFile << haps.sampNamesVec_[idx] << "\t" << "nogroup";
      } else {
        outFile << haps.sampNamesVec_[idx] << "\t" << group.first;
      }
      outFile << std::endl;
    }
  }
  OutputStream outGroupCountsFile(
    OutOptions(njh::files::make_path(setUp.pars_.directoryName_, "clusters_by_ccc_groupCounts.tsv")));
  outGroupCountsFile << "group\tsample_count";
  outGroupCountsFile << "\tmin_ccc\tmedian_ccc\tmean_ccc\tmax_ccc";
  outGroupCountsFile << std::endl;
  for (const auto &group: groupIndexes) {
    if (group.first == std::numeric_limits<uint32_t>::max()) {
      uint32_t low_coverage_not_clustered_cnt = 0;
      for (const auto &samp_idx: group.second) {
        if (lociCoveragePerSample[haps.sampNamesVec_[samp_idx]] < minimumLociCoverageToKeepSamples) {
          low_coverage_not_clustered_cnt++;
        }
      }
      outGroupCountsFile << "nogroup" << "\t" << group.second.size() - low_coverage_not_clustered_cnt;
      outGroupCountsFile << "\t" << "NA"
          << "\t" << "NA"
          << "\t" << "NA"
          << "\t" << "NA";
      outGroupCountsFile << std::endl;

      outGroupCountsFile << "low_coverage_not_clustered" << "\t" << low_coverage_not_clustered_cnt;
      outGroupCountsFile << "\t" << "NA"
          << "\t" << "NA"
          << "\t" << "NA"
          << "\t" << "NA";
      outGroupCountsFile << std::endl;
    } else {
      outGroupCountsFile << group.first << "\t" << group.second.size();
      outGroupCountsFile << "\t" << groups_ccc_stats[group.first]["min"]
          << "\t" << groups_ccc_stats[group.first]["median"]
          << "\t" << groups_ccc_stats[group.first]["mean"]
          << "\t" << groups_ccc_stats[group.first]["max"];
      outGroupCountsFile << std::endl;
    }
  }
	setUp.timer_.logLapTimes(setUp.rLog_.runLogFile_, true, 6, true);
	return 0;
}




int popGenExpRunner::doPairwiseComparisonOnHapsSharing(const njh::progutils::CmdArgs & inputCommands){
	bool writeOutDistMatrices = false;
	bool clusterOnJacardIndexShared = false;
	njhUndirWeightedGraph<double, std::shared_ptr<BasicPointMatrix<double>::BasicPoint>>::dbscanPars dbscanPars;
	// dbscanPars.eps_ = 0.50;
	dbscanPars.eps_ = 0.10;
	dbscanPars.minEpNeighbors_ = 2;
	double minimumLociCoverageToKeepSamples = 0.90;
	bfs::path metaFnp;
	VecStr metaFieldsToCalcPopDiffs{};
	HapsEncodedMatrix::SetWithExternalPars pars;
	bool writeOutTarsAbsoluteShared = false;
	bool doNotBreakWithRmse = false;
	double rmseCutOffToBreak = 0.10;
	bool doNotWriteOutGroupedRMSEs = false;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(minimumLociCoverageToKeepSamples, "--minimumLociCoverageToKeepSamples", "minimum Loci Coverage To Keep Samples in post analysis steps, must have reads for at least this frction of the total loci");

	setUp.setOption(clusterOnJacardIndexShared, "--clusterOnJacardIndexShared", "cluster On Jacard Index Shared");
	setUp.setOption(doNotBreakWithRmse, "--doNotBreakWithRmse", "do Not Break With Rmse");
	setUp.setOption(rmseCutOffToBreak, "--rmseCutOffToBreak", "rmse Cut Off To Break");
	setUp.setOption(doNotWriteOutGroupedRMSEs, "--doNotWriteOutGroupedRMSEs", "wriet Out Grouped RMSEs");
	bool writeOutGroupedRMSEs = !doNotWriteOutGroupedRMSEs;

	setUp.setOption(dbscanPars.eps_, "--eps", "Epsilon (distance sensitivity of algorithm)");
	setUp.setOption(dbscanPars.minEpNeighbors_, "--minpts", "The minimum number of epsilon neighbors");
	setUp.setOption(writeOutDistMatrices, "--writeOutDistMatrices", "write Out Dist Matrices");
	setUp.setOption(metaFnp, "--metaFnp", "Table of meta data for samples, needs a column named sample and each additional column will be the meta data associated with that sample");
	setUp.setOption(metaFieldsToCalcPopDiffs, "--metaFieldsToCalcPopDiffs", "Meta Fields To Calc Pop Diffs");
	setUp.setOption(writeOutTarsAbsoluteShared, "--writeOutTarsAbsoluteShared", "write Out Tars Absolute Shared");

  pars.setDefaults(setUp);

  setUp.processDirectoryOutputName(bfs::path(bfs::basename(pars.tableFnp)).string() + "_doPairwiseComparisonOnHapsSharing_TODAY", true);
	setUp.finishSetUp(std::cout);

	setUp.startARunLog(setUp.pars_.directoryName_);


	setUp.timer_.setLapName("initial");
	setUp.timer_.startNewLap("encode haplotypes");
  HapsEncodedMatrix haps(pars);

  if (!metaFnp.empty()) {
    haps.addMeta(metaFnp);
    if (!metaFieldsToCalcPopDiffs.empty()) {
      haps.meta_->checkForFieldsThrow(metaFieldsToCalcPopDiffs);
    }
  } else if (!metaFieldsToCalcPopDiffs.empty()) {
    haps.addMetaWithInputTab(njh::vecToSet(metaFieldsToCalcPopDiffs));
  }
	setUp.timer_.startNewLap("get hap probabilities");
	/**@todo look into whether or not this is being actually used */
	haps.calcHapProbs();

	haps.add_relative_abundance();
	setUp.timer_.startNewLap("writing sample info");


	setUp.timer_.startNewLap("get index measures");
	auto indexRes = haps.genIndexMeasures(setUp.pars_.verbose_);
	OutputStream outSampNamesOut(njh::files::make_path(setUp.pars_.directoryName_, "sampleNames.tab.txt"));
	outSampNamesOut << njh::conToStr(haps.sampNamesVec_, "\n") << std::endl;

	if(writeOutDistMatrices){
		OutputStream byTargetOut(njh::files::make_path(setUp.pars_.directoryName_, "percOfTarSharingAtLeastOneHap.tab.txt.gz"));
		OutputStream byHapOut(njh::files::make_path(setUp.pars_.directoryName_, "jacardByAllHap.tab.txt.gz"));
		OutputStream byHapTarSharedOut(njh::files::make_path(setUp.pars_.directoryName_, "jacardByHapsTarShared.tab.txt.gz"));
		OutputStream avgHapOut(njh::files::make_path(setUp.pars_.directoryName_, "avgJacardPerTarget.tab.txt.gz"));

		OutputStream byHapTarSharedWeightedOut(njh::files::make_path(setUp.pars_.directoryName_, "jacardByHapsTarSharedWeighted.tab.txt.gz"));
		OutputStream avgHapWeightedOut(njh::files::make_path(setUp.pars_.directoryName_, "avgJacardPerTargetWeighted.tab.txt.gz"));
		OutputStream targetsSharedBetweenSampsOut(njh::files::make_path(setUp.pars_.directoryName_, "targetsSharedBetweenSamps.tab.txt.gz"));



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

	std::unordered_map<std::string, double> lociCoveragePerSample = haps.getTargetCoveragePerSample();

	{
		table numTargetsPerSample = haps.getTableNumberTargetsPerSample(minimumLociCoverageToKeepSamples);
		OutputStream lociCoverageOut(njh::files::make_path(setUp.pars_.directoryName_, "loci_coverage_per_sample_info.tsv"));
		numTargetsPerSample.outPutContents(lociCoverageOut, "\t");
	}

	if(clusterOnJacardIndexShared) {
		auto distFnp = njh::files::make_path(setUp.pars_.directoryName_, "1MinusjacardByHapsTarShared.tab.txt.gz");
		{
			OutputStream byHapTarSharedOut(distFnp);
			for(const auto & ih : indexRes.byHapsTarShared){
				std::vector<double> outRow;
				outRow.reserve(ih.size());
				for(const auto & i : ih) {
					outRow.emplace_back(1 - i);
				}
				byHapTarSharedOut << njh::conToStr(outRow, "\t") << std::endl;
			}
		}



		OutputStream outFile(OutOptions(njh::files::make_path(setUp.pars_.directoryName_, "clusters_by_jacardTargetsShared.tsv")));
		njh::stopWatch watch;
		watch.setLapName("Reading in");
		auto mat = BasicPointMatrix<double>::readInBasicMatrix(distFnp, dbscanPars);
		watch.startNewLap("Adding nodes");
		std::vector<std::vector<double>> pairwiseRMSEs;
	  std::vector<std::vector<double>> pairwise_cccs;
		if(doNotBreakWithRmse) {
			// mat.setGraph(pars.numThreads, setUp.pars_.verbose_);
			mat.graph_ = std::make_unique<
					njhUndirWeightedGraph<double,
							std::shared_ptr<BasicPointMatrix<double>::BasicPoint>>>();

			for (const auto & pos : iter::range(mat.points_.size())) {
				mat.graph_->addNode(estd::to_string(pos), mat.points_[pos]);
			}

			uint32_t belowEp = 0;
			/**@todo this appears to be actually fairly slow, i think it's mostly because the eu calculations is so fast, perhaps a better way of multithreading this can be done
			 *
			 */
			PairwisePairFactory pairFactory(mat.points_.size());
			uint32_t pairBatchCount = 100000;
			std::mutex graphMut;
			struct PairDist {
				PairDist(const PairwisePairFactory::PairwisePair & pair, double dist) :
						pair_(pair), dist_(dist) {
				}
				PairwisePairFactory::PairwisePair pair_;
				double dist_;
			};

			std::function<void()> addToGraph =
					[&graphMut, &pairFactory,&pairBatchCount,&belowEp,
						&mat, &dbscanPars,
						&haps,
						&lociCoveragePerSample, &minimumLociCoverageToKeepSamples]() {
						PairwisePairFactory::PairwisePairVec pairs;
						std::vector<PairDist> belowEps;
						while(pairFactory.setNextPairs(pairs, pairBatchCount)) {
							for(const auto & pair : pairs.pairs_) {
								if (lociCoveragePerSample[haps.sampNamesVec_[pair.row_]] < minimumLociCoverageToKeepSamples ||
									lociCoveragePerSample[haps.sampNamesVec_[pair.col_]] < minimumLociCoverageToKeepSamples) {
									continue;
								}
								auto dist = mat.points_[pair.row_]->vals_[pair.col_];
								if (dist < dbscanPars.eps_) {
									belowEps.emplace_back(PairDist{pair, dist});
								}
							}
						}
						if(!belowEps.empty()) {
							std::lock_guard<std::mutex> lock(graphMut);
							belowEp += belowEps.size();
							for(const auto & bEps : belowEps) {
								mat.graph_->addEdge(estd::to_string(bEps.pair_.row_), estd::to_string(bEps.pair_.col_),
										bEps.dist_);
							}
						}
			};
			njh::concurrent::runVoidFunctionThreaded(addToGraph, pars.numThreads);
		} else {
			pairwiseRMSEs = std::vector<std::vector<double>>(haps.sampNames_.size(), std::vector<double>(haps.sampNames_.size(),1.0));
		  pairwise_cccs = std::vector<std::vector<double>>(haps.sampNames_.size(), std::vector<double>(haps.sampNames_.size(),0.0));
		  for(size_t pos = 0; pos < haps.sampNames_.size(); ++pos){
				//set diagonal
				pairwiseRMSEs[pos][pos] = 0;
		    pairwise_cccs[pos][pos] = 1.0;
			}

			if (0 == pars.numThreads) {
				pars.numThreads = 1;
			}

			mat.graph_ = std::make_unique<
					njhUndirWeightedGraph<double,
							std::shared_ptr<BasicPointMatrix<double>::BasicPoint>>>();

			for (const auto & pos : iter::range(mat.points_.size())) {
				mat.graph_->addNode(estd::to_string(pos), mat.points_[pos]);
			}

			uint32_t belowEp = 0;
			/**@todo this appears to be actually fairly slow, i think it's mostly because the eu calculations is so fast, perhaps a better way of multithreading this can be done
			 *
			 */
			PairwisePairFactory pairFactory(mat.points_.size());
			uint32_t pairBatchCount = 100000;
			std::mutex graphMut;
			struct PairDist {
				PairDist(const PairwisePairFactory::PairwisePair & pair, double dist) :
						pair_(pair), dist_(dist) {
				}
				PairwisePairFactory::PairwisePair pair_;
				double dist_;
			};

			std::function<void()> addToGraph =
					[&graphMut, &pairFactory,&pairBatchCount,&belowEp,&mat, &haps,&rmseCutOffToBreak,
						&pairwiseRMSEs, &pairwise_cccs, &lociCoveragePerSample,
						&minimumLociCoverageToKeepSamples]() {
						PairwisePairFactory::PairwisePairVec pairs;
						std::vector<PairDist> belowEps;
						while(pairFactory.setNextPairs(pairs, pairBatchCount)) {
							for (const auto &pair: pairs.pairs_) {
								if (lociCoveragePerSample[haps.sampNamesVec_[pair.row_]] < minimumLociCoverageToKeepSamples ||
								    lociCoveragePerSample[haps.sampNamesVec_[pair.col_]] < minimumLociCoverageToKeepSamples) {
									continue;
								}
								//auto dist = mat.points_[pair.row_]->euDist(*mat.points_[pair.col_]);
								auto dist = mat.points_[pair.row_]->vals_[pair.col_];
								if (dist < mat.dbscanPars_.eps_) {
								  std::vector<double> row_values;
								  std::vector<double> col_values;
									std::vector<double> rmses;
									double sum = 0;
									for(const auto tpos : iter::range(haps.tarNamesVec_.size())) {
										if(haps.targetsEncodeBySamp_[pair.col_][tpos]  + haps.targetsEncodeBySamp_[pair.row_][tpos] == 2) {
											double current_sum = 0;
											for(const auto hapPos : iter::range(haps.numberOfHapsPerTarget_[tpos])) {
												current_sum += std::pow(haps.hapsEncodeBySampRelAbund_[pair.col_][haps.tarStart_[tpos] + hapPos] - haps.hapsEncodeBySampRelAbund_[pair.row_][haps.tarStart_[tpos] + hapPos],2);
												sum +=         std::pow(haps.hapsEncodeBySampRelAbund_[pair.col_][haps.tarStart_[tpos] + hapPos] - haps.hapsEncodeBySampRelAbund_[pair.row_][haps.tarStart_[tpos] + hapPos],2);
											  row_values.emplace_back(haps.hapsEncodeBySampRelAbund_[pair.row_][haps.tarStart_[tpos] + hapPos]);
											  col_values.emplace_back(haps.hapsEncodeBySampRelAbund_[pair.col_][haps.tarStart_[tpos] + hapPos]);
											}
											rmses.emplace_back(std::sqrt(current_sum));
										}
									}
									//only add if mean RMSE is less than the cut off
									//if(vectorMean(rmses) < rmseCutOffToBreak) {
								  auto ccc_calc = ConcordanceCalculator::lins_ccc_with_ci(row_values, col_values);
									pairwiseRMSEs[pair.col_][pair.row_] = std::sqrt(sum/rmses.size());
									pairwiseRMSEs[pair.row_][pair.col_] = std::sqrt(sum/rmses.size());
								  pairwise_cccs[pair.col_][pair.row_] = ccc_calc.ccc;
								  pairwise_cccs[pair.row_][pair.col_] = ccc_calc.ccc;
									// std::cout << __FILE__ << " : " << __LINE__ << std::endl;
									// std::cout << "lociCoveragePerSample[haps.sampNamesVec_[pair.row_]]: " << lociCoveragePerSample[haps.sampNamesVec_[pair.row_]] << std::endl;
									// std::cout << "lociCoveragePerSample[haps.sampNamesVec_[pair.col_]]: " << lociCoveragePerSample[haps.sampNamesVec_[pair.col_]] << std::endl;
									// std::cout << "haps.sampNamesVec_[pair.row_]: " << haps.sampNamesVec_[pair.row_] << std::endl;
									// std::cout << "haps.sampNamesVec_[pair.col_]: " << haps.sampNamesVec_[pair.col_] << std::endl;
									// std::cout << "rmses.size(): " << rmses.size() << std::endl;
									// std::cout << "std::sqrt(sum/rmses.size()): " << std::sqrt(sum/rmses.size()) << std::endl;
									// std::cout << "rmseCutOffToBreak          : " << rmseCutOffToBreak << std::endl << std::endl;

									if(std::sqrt(sum/rmses.size()) < rmseCutOffToBreak){
										belowEps.emplace_back(PairDist{pair, dist});
									}
								}
							}
						}
						if(!belowEps.empty()) {
							std::lock_guard<std::mutex> lock(graphMut);
							belowEp += belowEps.size();
							for(const auto & bEps : belowEps) {
								mat.graph_->addEdge(estd::to_string(bEps.pair_.row_), estd::to_string(bEps.pair_.col_),
										bEps.dist_);
							}
						}
			};

			njh::concurrent::runVoidFunctionThreaded(addToGraph, pars.numThreads);

			if (setUp.pars_.verbose_) {
				std::cout << std::endl;
				std::cout << "below: " << belowEp << "/" << pairFactory.totalCompares_ << std::endl;
			}


		}


		watch.startNewLap("dbscan");
		mat.graph_->dbscan(dbscanPars);
		// if(!doNotBreakWithRmse){
		// 						//first fill the relative abundance vector with the input relative abundance
		// 	std::vector<std::vector<double>> hapsEncodeBySampRelAbund = std::vector<std::vector<double>> (haps.sampNames_.size());
		// 	for(const auto & samp : haps.sampNames_){
		// 		hapsEncodeBySampRelAbund[haps.sampNamesKey_[samp]] = std::vector<double>(haps.totalHaps_, 0);
		// 	}
		// 	TableReader reReadHapTab(TableIOOpts(InOptions(haps.pars_.tableFnp), "\t", true));
		// 	VecStr row;
		// 	while(reReadHapTab.getNextRow(row)){
		// 		const auto& samp = row[reReadHapTab.header_.getColPos(haps.pars_.sampleCol)];
		// 		const auto& tar = row[reReadHapTab.header_.getColPos(haps.pars_.targetNameCol)];
		// 		if(!haps.pars_.selectSamples.empty() && !njh::in(samp, haps.pars_.selectSamples)){
		// 			continue;
		// 		}
		// 		if(!haps.pars_.selectTargets.empty() && !njh::in(tar, haps.pars_.selectTargets)){
		// 			continue;
		// 		}
		// 		const auto& hapName = row[reReadHapTab.header_.getColPos(haps.pars_.popIDCol)];
		// 		auto rBund = njh::StrToNumConverter::stoToNum<double>(row[reReadHapTab.header_.getColPos(haps.pars_.relAbundCol)]);
		// 		auto tKey = haps.tarNameKey_[tar];
		// 		auto hKey = haps.hapNamesKey_[tar][hapName];
		// 		//				if(rBund > 0 && rBund < 1){
		// 		//					std::cout << rBund << std::endl;
		// 		//				}
		// 		hapsEncodeBySampRelAbund[haps.sampNamesKey_[samp]][haps.tarStart_[tKey] + hKey] = rBund;
		// 	}
		// 	//now recalculate the relative abundance to be 0-1
		// 	for(const auto pos : iter::range(haps.sampNames_.size())) {
		// 		for(const auto tpos : iter::range(haps.tarNamesVec_.size())) {
		// 			if(haps.targetsEncodeBySamp_[pos][tpos] == 1) {
		// 				double sum = 0;
		// 				for(const auto hapPos : iter::range(haps.numberOfHapsPerTarget_[tpos])) {
		// 					sum += hapsEncodeBySampRelAbund[pos][haps.tarStart_[tpos] + hapPos];
		// 				}
		// 				for(const auto hapPos : iter::range(haps.numberOfHapsPerTarget_[tpos])) {
		// 					hapsEncodeBySampRelAbund[pos][haps.tarStart_[tpos] + hapPos] = hapsEncodeBySampRelAbund[pos][haps.tarStart_[tpos] + hapPos]/sum;
		// 				}
		// 			}
		// 		}
		// 	}
		// 	// for(const auto pos : iter::range(haps.sampNamesKey_.size())) {
		// 	// 	std::cout << njh::conToStr( hapsEncodeBySampRelAbund[pos], "\t") << std::endl;
		// 	// }
		// 	std::map<uint32_t, std::vector<uint32_t>> groupIndexes;
		// 	for(const auto & n : iter::enumerate(mat.graph_->nodes_)) {
		// 		groupIndexes[n.element->group_].emplace_back(n.index);
		// 	}
		// 	for(const auto & group : groupIndexes) {
		// 		PairwisePairFactory pair_factory(group.second.size());
		// 		PairwisePairFactory::PairwisePair pair;
		// 		while(pair_factory.setNextPair(pair)) {
		// 			std::vector<double> rmses;
		// 			double sum = 0;
		// 			for(const auto tpos : iter::range(haps.tarNamesVec_.size())) {
		// 				if(haps.targetsEncodeBySamp_[group.second[pair.col_]][tpos]  + haps.targetsEncodeBySamp_[group.second[pair.row_]][tpos] == 2) {
		// 					double current_sum = 0;
		// 					for(const auto hapPos : iter::range(haps.numberOfHapsPerTarget_[tpos])) {
		// 						current_sum += std::pow(hapsEncodeBySampRelAbund[group.second[pair.col_]][haps.tarStart_[tpos] + hapPos] - hapsEncodeBySampRelAbund[group.second[pair.row_]][haps.tarStart_[tpos] + hapPos],2);
		// 						sum += std::pow(hapsEncodeBySampRelAbund[group.second[pair.col_]][haps.tarStart_[tpos] + hapPos] - hapsEncodeBySampRelAbund[group.second[pair.row_]][haps.tarStart_[tpos] + hapPos],2);
		// 					}
		// 					rmses.emplace_back(std::sqrt(current_sum));
		// 				}
		// 			}
		// 			std::cout << haps.sampNamesVec_[group.second[pair.col_]] << "\t" <<  haps.sampNamesVec_[group.second[pair.row_]] << "\t" << std::sqrt(sum) << "\t" << std::sqrt(sum)/rmses.size() << "\t" << std::sqrt(sum/rmses.size()) << "\t" << vectorMean(rmses) << std::endl;
		// 		}
		// 	}
		// }

		watch.startNewLap("output");
		//mat.writeGraph(outFile);
		outFile << "sample\tgroup";
		outFile << std::endl;
		std::map<uint32_t, std::vector<uint32_t>> groupIndexes;
		std::vector<uint32_t> allGroupedIndices;
		for(const auto & n : iter::enumerate(mat.graph_->nodes_)) {
			groupIndexes[n.element->group_].emplace_back(n.index);
			if (std::numeric_limits<uint32_t>::max() != n.element->group_) {
				allGroupedIndices.emplace_back(n.index);
			}
		}
		if (!doNotBreakWithRmse) {
			//fill in all RMSEs now for the grouped data
			watch.startNewLap("calculate all RMSEs");
			uint32_t pairBatchCount = 100;
			PairwisePairFactory allGroupedSamplesFactory(allGroupedIndices.size());
			std::function<void()> calcRMSEs =
					[&allGroupedSamplesFactory,
					  &pairwiseRMSEs,
					  &pairwise_cccs,
					  &pairBatchCount,
						&haps,&allGroupedIndices]() {
						PairwisePairFactory::PairwisePairVec pairs;
						while(allGroupedSamplesFactory.setNextPairs(pairs, pairBatchCount)) {
							for(const auto & groupedPair : pairs.pairs_) {
								auto colSamplePos = allGroupedIndices[groupedPair.col_];
								auto rowSamplePos = allGroupedIndices[groupedPair.row_];
								std::vector<double> rmses;
							  std::vector<double> row_values;
							  std::vector<double> col_values;
								double sum = 0;
								for(const auto tpos : iter::range(haps.tarNamesVec_.size())) {
									if(haps.targetsEncodeBySamp_[colSamplePos][tpos]  + haps.targetsEncodeBySamp_[rowSamplePos][tpos] == 2) {
										double current_sum = 0;
										for(const auto hapPos : iter::range(haps.numberOfHapsPerTarget_[tpos])) {
											current_sum += std::pow(haps.hapsEncodeBySampRelAbund_[colSamplePos][haps.tarStart_[tpos] + hapPos] - haps.hapsEncodeBySampRelAbund_[rowSamplePos][haps.tarStart_[tpos] + hapPos],2);
											sum +=         std::pow(haps.hapsEncodeBySampRelAbund_[colSamplePos][haps.tarStart_[tpos] + hapPos] - haps.hapsEncodeBySampRelAbund_[rowSamplePos][haps.tarStart_[tpos] + hapPos],2);
										  row_values.emplace_back(haps.hapsEncodeBySampRelAbund_[rowSamplePos][haps.tarStart_[tpos] + hapPos]);
										  col_values.emplace_back(haps.hapsEncodeBySampRelAbund_[colSamplePos][haps.tarStart_[tpos] + hapPos]);
										}
										rmses.emplace_back(std::sqrt(current_sum));
									}
								}
								//only add if mean RMSE is less than the cut off
								//if(vectorMean(rmses) < rmseCutOffToBreak) {
							  auto ccc_calc = ConcordanceCalculator::lins_ccc_with_ci(row_values, col_values);
								pairwiseRMSEs[colSamplePos][rowSamplePos] = std::sqrt(sum/rmses.size());
								pairwiseRMSEs[rowSamplePos][colSamplePos] = std::sqrt(sum/rmses.size());
							  pairwise_cccs[colSamplePos][rowSamplePos] = ccc_calc.ccc;
							  pairwise_cccs[rowSamplePos][colSamplePos] = ccc_calc.ccc;
							}
						}
			};

			njh::concurrent::runVoidFunctionThreaded(calcRMSEs, pars.numThreads);

		}
		std::unordered_map<uint32_t, std::map<std::string, double>> groups_jaccard_stats;
		std::unordered_map<uint32_t, std::map<std::string, double>> groups_rmse_stats;
	  std::unordered_map<uint32_t, std::map<std::string, double>> groups_ccc_stats;
		for(const auto & group : groupIndexes) {

			if (!doNotBreakWithRmse && std::numeric_limits<uint32_t>::max() != group.first) {
				std::vector<double> rmsesWithinGroup;
			  std::vector<double> cccsWithinGroup;
				PairwisePairFactory pfac(group.second.size());
				PairwisePairFactory::PairwisePair pair;
				// std::cout << "group: " << group.first << std::endl;
				while (pfac.setNextPair(pair)) {
					// std::cout << haps.sampNamesVec_[group.second[pair.col_]] << " vs " << haps.sampNamesVec_[group.second[pair.row_]] << " rmse: " << pairwiseRMSEs[group.second[pair.col_]][group.second[pair.row_]] << std::endl;
					rmsesWithinGroup.emplace_back(pairwiseRMSEs[group.second[pair.col_]][group.second[pair.row_]]);
				  cccsWithinGroup.emplace_back(pairwise_cccs[group.second[pair.col_]][group.second[pair.row_]]);
				}
				groups_rmse_stats[group.first] = getStatsOnVec(rmsesWithinGroup);
			  groups_ccc_stats[group.first]  = getStatsOnVec(cccsWithinGroup);
				// std::cout << njh::conToStr(rmsesWithinGroup, ",") << std::endl;
				// std::cout << "stats: " << njh::json::toJson(stats) << std::endl;
				// std::cout << std::endl;
			}


			if (std::numeric_limits<uint32_t>::max() != group.first) {
				std::vector<double> jacardWithinGroup;
				PairwisePairFactory pfac(group.second.size());
				PairwisePairFactory::PairwisePair pair;
				// std::cout << "group: " << group.first << std::endl;
				while (pfac.setNextPair(pair)) {
					// std::cout << haps.sampNamesVec_[group.second[pair.col_]] << " vs " << haps.sampNamesVec_[group.second[pair.row_]] << " rmse: " << pairwiseRMSEs[group.second[pair.col_]][group.second[pair.row_]] << std::endl;
					jacardWithinGroup.emplace_back(1 - mat.points_[group.second[pair.col_]]->vals_[group.second[pair.row_]]);
				}
				groups_jaccard_stats[group.first] = getStatsOnVec(jacardWithinGroup);
			}
			for (const auto &idx: group.second) {
				if (lociCoveragePerSample[haps.sampNamesVec_[idx]] < minimumLociCoverageToKeepSamples) {
					outFile << haps.sampNamesVec_[idx] << "\t" << "low_coverage_not_clustered";
				} else if (group.first == std::numeric_limits<uint32_t>::max()) {
					outFile << haps.sampNamesVec_[idx] << "\t" << "nogroup";
				} else {
					outFile << haps.sampNamesVec_[idx] << "\t" << group.first;
				}
				outFile << std::endl;
			}
		}
		if (writeOutGroupedRMSEs) {
			OutputStream outGroupCountsFile(OutOptions(njh::files::make_path(setUp.pars_.directoryName_, "rmses_of_grouped_samples.tsv.gz")));
			outGroupCountsFile << "sample";
			for (const auto & sample : allGroupedIndices) {
				outGroupCountsFile << "\t" << haps.sampNamesVec_[sample];
			}
			outGroupCountsFile << std::endl;
			for (const auto & sampleCol : allGroupedIndices) {
				outGroupCountsFile << haps.sampNamesVec_[sampleCol];
				for (const auto & sampleRow : allGroupedIndices) {
					outGroupCountsFile << "\t" << pairwiseRMSEs[sampleRow][sampleCol];
				}
				outGroupCountsFile << std::endl;
			}
		}

	  if (writeOutGroupedRMSEs) {
	    OutputStream outGroupCountsFile(OutOptions(njh::files::make_path(setUp.pars_.directoryName_, "ccc_of_grouped_samples.tsv.gz")));
	    outGroupCountsFile << "sample";
	    for (const auto & sample : allGroupedIndices) {
	      outGroupCountsFile << "\t" << haps.sampNamesVec_[sample];
	    }
	    outGroupCountsFile << std::endl;
	    for (const auto & sampleCol : allGroupedIndices) {
	      outGroupCountsFile << haps.sampNamesVec_[sampleCol];
	      for (const auto & sampleRow : allGroupedIndices) {
	        outGroupCountsFile << "\t" << pairwise_cccs[sampleRow][sampleCol];
	      }
	      outGroupCountsFile << std::endl;
	    }
	  }


		OutputStream outGroupCountsFile(OutOptions(njh::files::make_path(setUp.pars_.directoryName_, "clusters_by_jacardTargetsShared_groupCounts.tsv")));
		outGroupCountsFile << "group\tsampleCount";
		outGroupCountsFile << "\tmin_jaccard\tmedian_jaccard\tmean_jaccard\tmax_jaccard";
		if (!doNotBreakWithRmse) {
			outGroupCountsFile << "\tmin_rmse\tmedian_rmse\tmean_rmse\tmax_rmse";
		}
	  outGroupCountsFile << "\tmin_ccc\tmedian_ccc\tmean_ccc\tmax_ccc";
		outGroupCountsFile << std::endl;
		for(const auto & group : groupIndexes) {
			if(group.first == std::numeric_limits<uint32_t>::max()) {
				uint32_t low_coverage_not_clustered_cnt = 0;
				for(const auto & samp_idx : group.second) {
					if (lociCoveragePerSample[haps.sampNamesVec_[samp_idx]] < minimumLociCoverageToKeepSamples) {
						low_coverage_not_clustered_cnt++;
					}
				}
				outGroupCountsFile << "nogroup" << "\t" << group.second.size() - low_coverage_not_clustered_cnt;
				outGroupCountsFile << "\t" << "NA"
						<< "\t" << "NA"
						<< "\t" << "NA"
						<< "\t" << "NA";
        if (!doNotBreakWithRmse) {
          outGroupCountsFile << "\t" << "NA"
              << "\t" << "NA"
              << "\t" << "NA"
              << "\t" << "NA";
        }
        outGroupCountsFile << "\t" << "NA"
            << "\t" << "NA"
            << "\t" << "NA"
            << "\t" << "NA";
				outGroupCountsFile << std::endl;

				outGroupCountsFile << "low_coverage_not_clustered" << "\t" << low_coverage_not_clustered_cnt;
				outGroupCountsFile << "\t" << "NA"
						<< "\t" << "NA"
						<< "\t" << "NA"
						<< "\t" << "NA";
				if (!doNotBreakWithRmse) {
					outGroupCountsFile << "\t" << "NA"
							<< "\t" << "NA"
							<< "\t" << "NA"
							<< "\t" << "NA";
				}
        outGroupCountsFile << "\t" << "NA"
            << "\t" << "NA"
            << "\t" << "NA"
            << "\t" << "NA";
				outGroupCountsFile << std::endl;
			} else {
				outGroupCountsFile << group.first << "\t" << group.second.size();
				outGroupCountsFile << "\t" << groups_jaccard_stats[group.first]["min"]
						<< "\t" << groups_jaccard_stats[group.first]["median"]
						<< "\t" << groups_jaccard_stats[group.first]["mean"]
						<< "\t" << groups_jaccard_stats[group.first]["max"];
				if (!doNotBreakWithRmse) {
					outGroupCountsFile << "\t" << groups_rmse_stats[group.first]["min"]
							<< "\t" << groups_rmse_stats[group.first]["median"]
							<< "\t" << groups_rmse_stats[group.first]["mean"]
							<< "\t" << groups_rmse_stats[group.first]["max"];
				}
        outGroupCountsFile << "\t" << groups_ccc_stats[group.first]["min"]
            << "\t" << groups_ccc_stats[group.first]["median"]
            << "\t" << groups_ccc_stats[group.first]["mean"]
            << "\t" << groups_ccc_stats[group.first]["max"];
				outGroupCountsFile << std::endl;
			}
		}
		if(setUp.pars_.verbose_){
			watch.logLapTimes(std::cout, true, 6, true);
		}
	}

	if(writeOutTarsAbsoluteShared){
		OutOptions hapsAbsoluteSharedBetweenSampsBetweenTargetsOutopts(njh::files::make_path(setUp.pars_.directoryName_, "hapsAbsoluteSharedBetweenSampsBetweenTargetsOut.tab.txt.gz"));
		haps.writeAbsoluteHapSharedPerSamplePerTar(hapsAbsoluteSharedBetweenSampsBetweenTargetsOutopts, setUp.pars_.verbose_);
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
							hapsForTarget[tarpos].unweighted_count_ += 1;
							hapsForTarget[tarpos].weighted_count_ += haps.hapsEncodeBySampRelAbund_[sampPos][haps.tarStart_[tarKey] + tarpos];
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
															<< "\t" << (std::numeric_limits<long double>::max() == diversityForTar.expected_k_heterozygosities.at(3).k_heterozygosity_ ? "NA": estd::to_string(diversityForTar.expected_k_heterozygosities.at(3).k_heterozygosity_))
															<< "\t" << (std::numeric_limits<long double>::max() == diversityForTar.expected_k_heterozygosities.at(4).k_heterozygosity_ ? "NA": estd::to_string(diversityForTar.expected_k_heterozygosities.at(4).k_heterozygosity_))
															<< "\t" << (std::numeric_limits<long double>::max() == diversityForTar.expected_k_heterozygosities.at(5).k_heterozygosity_ ? "NA": estd::to_string(diversityForTar.expected_k_heterozygosities.at(5).k_heterozygosity_))
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
								hapsForTargetPerPopulationRaw[sampleToMeta[sampPos]][tarpos].unweighted_count_ +=1;
								hapsForTargetPerPopulationRaw[sampleToMeta[sampPos]][tarpos].weighted_count_ += haps.hapsEncodeBySampRelAbund_[sampPos][haps.tarStart_[tarKey] + tarpos];

							}
						}
					}
					std::unordered_map<std::string, std::vector<PopGenCalculator::PopHapInfo>> hapsForTargetPerPopulation;
					for(const auto & pop : hapsForTargetPerPopulationRaw){
						for(const auto & hap : pop.second){
							if(hap.unweighted_count_ > 0){
								hapsForTargetPerPopulation[pop.first].emplace_back(hap);
							}
						}
					}

					PopGenCalculator::PopDifferentiationMeasures generalDiff;
					if(hapsForTargetPerPopulation.size() > 1){
						generalDiff = PopGenCalculator::getOverallPopDiffWeighted(hapsForTargetPerPopulation);
					}
					std::unordered_map<std::string, std::unordered_map<std::string, PopGenCalculator::PopDifferentiationMeasuresPairWise>> pairwiseDiffs;

					if(hapsForTargetPerPopulation.size() > 1){
						pairwiseDiffs = PopGenCalculator::getPairwisePopDiffWeighted(hapsForTargetPerPopulation);
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
							<< "\t" << (std::numeric_limits<long double>::max() == popDiv.second.expected_k_heterozygosities.at(3).k_heterozygosity_ ? "NA": estd::to_string(popDiv.second.expected_k_heterozygosities.at(3).k_heterozygosity_))
							<< "\t" << (std::numeric_limits<long double>::max() == popDiv.second.expected_k_heterozygosities.at(4).k_heterozygosity_ ? "NA": estd::to_string(popDiv.second.expected_k_heterozygosities.at(4).k_heterozygosity_))
							<< "\t" << (std::numeric_limits<long double>::max() == popDiv.second.expected_k_heterozygosities.at(5).k_heterozygosity_ ? "NA": estd::to_string(popDiv.second.expected_k_heterozygosities.at(5).k_heterozygosity_))
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
