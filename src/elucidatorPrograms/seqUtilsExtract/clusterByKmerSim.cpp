//
// Created by Nicholas Hathaway on 4/9/22.
//

#include "seqUtilsExtractRunner.hpp"

#include "elucidator/utils/KmerUtils.hpp"
#include "elucidator/objects/seqObjects/seqKmers.h"
#include <njhseq/IO/SeqIO.h>
#include <njhseq/concurrency/PairwisePairFactory.hpp>
#include <njhseq/concurrency/AllByAllPairFactory.hpp>

#undef BOOST_HAS_THREADS
#include <boost/math/statistics/t_test.hpp>

namespace njhseq {


std::string distGraphRmdReport = R"(---
title: "Plot Dist Graph"
author: "Nicholas Hathaway"
output:
 html_document:
   highlight: textmate
   theme: flatly
   code_folding: hide
   toc: yes
   toc_float: yes
   fig_width: 12
   fig_height : 8
---


```{r setup, echo=FALSE, message=FALSE}
require(knitr)
require(tidyverse)
require(heatmaply)

myFormula= x~y
library(ggpmisc)
`%!in%` <- Negate(`%in%`)
opts_chunk$set(message=FALSE, warning=FALSE, comment = "", cache = F)
options(width = 200)
```
<style type="text/css">
div.main-container {
  max-width: 1800px;
  margin-left: auto;
  margin-right: auto;
}
</style>



```{r}
raw_dist = readr::read_tsv("distance.tab.txt", col_names = T)
dist = as.dist(raw_dist, diag = T)
dist_mat = as.matrix(dist)
dist_mat = dist_mat
heatmaply(dist_mat, plot_method = "plotly")

```

```{r}
parameters = readr::read_tsv("parameters.tab.txt")
parameters_eps = parameters %>%
  filter(parameter == "epsilon")
dist_mat[dist_mat > parameters_eps$value ] = 0
heatmaply(dist_mat, plot_method = "plotly")

```
)";

int seqUtilsExtractRunner::clusterByKmerSim(const njh::progutils::CmdArgs & inputCommands) {
	uint32_t kmerLength = 7;
	uint32_t numThreads = 1;
	bool HDBScountZeroNeighbors = true;
	bool HDBSCountSingletGroups = false;
	bool HDBSredetermineMaxEps = false;
	double HDBSmaxInitialEps = std::numeric_limits<double>::max();
	uint32_t proposedClusters = std::numeric_limits<uint32_t>::max();
	readDistGraph<double>::dbscanPars dbPars_;
	dbPars_.minEpNeighbors_ = 2;
	dbPars_.eps_ = .20;

	bool useHDBS = false;

	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processReadInNames(true);

	setUp.setOption(HDBSredetermineMaxEps, "--HDBSredetermineMaxEps", "HDBS redetermine Max Eps allowed in initial step based by setting it equal to mean of the non-same-group dist minus 2sd");
	bool doNotCountZeroNeighbors = false;
	setUp.setOption(doNotCountZeroNeighbors, "--HDBSdoNotCountZeroNeighbors", "HDBS when doing j-th nearest neighbor do Not Count Zero Neighbors");
	HDBScountZeroNeighbors = !doNotCountZeroNeighbors;
	setUp.setOption(HDBSmaxInitialEps, "--HDBSmaxInitialEps", "HDBS a hard cut off for max Initial Eps for initial DBSCAN step in H-DBSCAN");
	setUp.setOption(proposedClusters, "--HDBSproposedClusters", "HDBS proposed number of clusters");
	setUp.setOption(HDBSCountSingletGroups, "--HDBSCountSingletGroups", "For HD DBscan count Singlet Groups, by default these are not included in towards the proposed group counts");

	setUp.setOption(numThreads, "--numThreads", "numThreads");
	setUp.setOption(kmerLength, "--kmerLength", "kmer Length");
	setUp.setOption(useHDBS, "--useHDBS", "useHDBS");

	setUp.setOption(dbPars_.eps_, "--epsilon", "min distance for connecting neighbors");
	setUp.setOption(dbPars_.minEpNeighbors_, "--minEpNeighbors", "min epsilon neighbors");
	setUp.processDirectoryOutputName(true);
	setUp.finishSetUp(std::cout);

	setUp.startARunLog(setUp.pars_.directoryName_);
	njh::stopWatch watch;
	watch.setLapName("Index Kmers");

	auto reads = createKmerReadVec(setUp.pars_.ioOptions_, kmerLength, false);


	std::function<
					double(const std::unique_ptr<seqWithKmerInfo> &,
								 const std::unique_ptr<seqWithKmerInfo> &)> disFun =
					[](const std::unique_ptr<seqWithKmerInfo> & read1,
						 const std::unique_ptr<seqWithKmerInfo> & read2) {
						auto dist = read1->compareKmers(*read2);
						return 1 - dist.second;
					};
	watch.startNewLap("Distance");

	auto dist = getDistance(reads, numThreads, disFun);
	readDistGraph<double> distGraph(dist, reads);

	if(useHDBS){
		auto hdbsDir = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar{"HDBScan"});
		OutputStream hdbsRunInfoOut(njh::files::make_path(hdbsDir, "runInfo.tab.txt"));
		hdbsRunInfoOut << "Run\tnonSingletGroupCounts\ttotalGroupCount\tlowestCentroidDist" << std::endl;
		OutputStream nearestNeibhorDistsOut(njh::files::make_path(hdbsDir, "nearestNeighborsDists.tab.txt"));
		nearestNeibhorDistsOut << "name\tneighbor\tnonZeroNeighborPos\tdist" << std::endl;
		uint32_t minimum_minPts = 3;
		if(std::numeric_limits<uint32_t>::max() == proposedClusters){
			proposedClusters = reads.size()/minimum_minPts;
		}
		//HDBS

		uint32_t maximum_minPts = reads.size()/proposedClusters;

		double minEps = std::numeric_limits<double>::max();
		std::vector<double> neighbor2Dists;
		for(const auto & n : distGraph.nodes_){
			//sort so edges are sorted from lowest dist to highest dist
			n->sortEdges([](const std::shared_ptr<readDistGraph<double>::edge> & e1,
											const std::shared_ptr<readDistGraph<double>::edge> & e2){
				return *e1 < *e2;
			});

			uint32_t total = 0;
			uint32_t count = 0;
			for(const auto & e : n->edges_){
				++total;
				if(HDBScountZeroNeighbors || e->dist_ != 0){
					++count;
					nearestNeibhorDistsOut << n->value_->name_
																 << "\t" << total
																 << "\t" << count
																 << "\t" << e->dist_ << std::endl;
					if((minimum_minPts - 1) == count){
						neighbor2Dists.emplace_back(e->dist_);
						if(e->dist_ < minEps){
							minEps = e->dist_;
						}
						break;
					}
				}
			}
		}
		njh::sort(neighbor2Dists);
		std::vector<double> eps;
		uint32_t sq = std::round(std::sqrt(reads.size()));
		for(uint32_t i = sq; i < reads.size(); i += sq){
			eps.emplace_back(neighbor2Dists[i]);
		}
		double maxEps = std::numeric_limits<double>::lowest();
		for(const auto & n : distGraph.nodes_){
			//since all nodes are connected no need to check size edge
			if(n->edges_[maximum_minPts]->dist_ > maxEps){
				maxEps = n->edges_[maximum_minPts]->dist_;
			}
		}
		//std::cout << __FILE__ << " : " << __LINE__ << " : " << __PRETTY_FUNCTION__ << std::endl;
		if(setUp.pars_.debug_){
			std::cout << "Eps: " << std::endl;
			for(const auto & epEnum : iter::enumerate(eps)){
				std::cout << "\t"<< epEnum.index << ": " << epEnum.second << std::endl;
			}
			std::cout << "Max Eps: " << maxEps << std::endl;
			std::cout << "maximum_minPts: " << maximum_minPts << std::endl;
			std::cout << "medianNearestNeighbor: " << vectorMedianRef(neighbor2Dists) << std::endl;
			std::cout << "meanNearestNeighbor: " << vectorMean(neighbor2Dists) << std::endl;
			std::cout << "  sdNearestNeighbor: " << vectorStandardDeviationPop(neighbor2Dists) << std::endl;
		}
		//reset node's visited and group values
		distGraph.resetAllNodes();
		//turn off whole graph
		distGraph.turnOffAllCons();
		distGraph.turnOffAllNodes();
		//set number of groups to be 0
		distGraph.numberOfGroups_ = 0;
		for(const auto & epEnum : iter::enumerate(eps)){
			if(epEnum.second > HDBSmaxInitialEps){
				break;
			}
			readDistGraph<double>::dbscanPars currentPars;
			currentPars.eps_ = epEnum.second;
			currentPars.minEpNeighbors_ = minimum_minPts;
			for (auto & n : distGraph.nodes_) {
				//if the node has not be visited by an expand or spread try to expand it
				if (!n->visited_) {
					//std::cout << n->value_->name_ << std::endl;
					n->dbscanExpand(distGraph.numberOfGroups_, currentPars);
					//if it was assigned a group and expanded, increase group number
					if (std::numeric_limits<uint32_t>::max() != n->group_) {
						++distGraph.numberOfGroups_;
					}
				}
			}
			if(epEnum.index + 1 != eps.size() && eps[epEnum.first + 1] < HDBSmaxInitialEps ) {
				//reset unclustered nodes back on and unvisted for the next eps
				for(auto & n : distGraph.nodes_){
					if(!n->on_){
						n->visitedAmount_ = 0;
						n->visited_ = false;
						n->on_ = true;
					}
				}
			}
		}
		distGraph.assignNoiseNodesAGroup();
		if(HDBSredetermineMaxEps){
			std::vector<double> notSameGroupDists;
			for(const auto & n : distGraph.nodes_){
				for(const auto & e : n->edges_){
					if(e->nodeToNode_[n->name_].lock()->group_ != n->group_){
						notSameGroupDists.emplace_back(e->dist_);
					}
				}
			}
			auto newHDBSmaxInitialEps = vectorMean(notSameGroupDists) - 2 * vectorStandardDeviationPop(notSameGroupDists);

			if(setUp.pars_.debug_){
				std::cout << "medianNotSameGroupDists: " << vectorMedianRef(notSameGroupDists) << std::endl;
				std::cout << "meanNotSameGroupDists: " << vectorMean(notSameGroupDists) << std::endl;
				std::cout << "  sdNotSameGroupDists: " << vectorStandardDeviationPop(notSameGroupDists) << std::endl;
				std::cout << std::endl;
				std::cout << "max eps based on not same group =" << newHDBSmaxInitialEps << std::endl;
			}
			//reset node's visited and group values
			distGraph.resetAllNodes();
			//turn off whole graph
			distGraph.turnOffAllCons();
			distGraph.turnOffAllNodes();
			//set number of groups to be 0
			distGraph.numberOfGroups_ = 0;
			for(const auto & epEnum : iter::enumerate(eps)){
				if(epEnum.second > newHDBSmaxInitialEps){
					break;
				}
				readDistGraph<double>::dbscanPars currentPars;
				currentPars.eps_ = epEnum.second;
				currentPars.minEpNeighbors_ = minimum_minPts;
				for (auto & n : distGraph.nodes_) {
					//if the node has not be visited by an expand or spread try to expand it
					if (!n->visited_) {
						//std::cout << n->value_->name_ << std::endl;
						n->dbscanExpand(distGraph.numberOfGroups_, currentPars);
						//if it was assigned a group and expanded, increase group number
						if (std::numeric_limits<uint32_t>::max() != n->group_) {
							++distGraph.numberOfGroups_;
						}
					}
				}
				if(epEnum.index + 1 != eps.size() && eps[epEnum.first + 1] < newHDBSmaxInitialEps ) {
					//reset unclustered nodes back on and unvisted for the next eps
					for(auto & n : distGraph.nodes_){
						if(!n->on_){
							n->visitedAmount_ = 0;
							n->visited_ = false;
							n->on_ = true;
						}
					}
				}
			}
			distGraph.assignNoiseNodesAGroup();
		}
    std::vector<double> differentInitialGroupDists;

    double lowestCentroidDistInitial = std::numeric_limits<double>::max();
		std::vector<std::vector<double>> centroidDistances;
		for(const auto pos : iter::range(distGraph.numberOfGroups_)) {
			centroidDistances.emplace_back(std::vector<double>(pos + 1));
		}
		//set the zeros
		for(const auto pos : iter::range(distGraph.numberOfGroups_)) {
			centroidDistances[pos][pos] = 0;
		}
		{
			//initial centroid info
			//compute centroid distances
			std::map<uint32_t, std::vector<uint32_t>> groupNodes;
			for(const auto & node : distGraph.nodes_){
				groupNodes[node->group_].emplace_back(distGraph.nameToNodePos_[node->name_]);
			}
			auto groups = getVectorOfMapKeys(groupNodes);
			PairwisePairFactory pFac(groups.size());
			PairwisePairFactory::PairwisePair pair;
			while(pFac.setNextPair(pair)){
				//group1 == row_
				//group2 == col_
				uint32_t group1 = groups[pair.row_];
				uint32_t group2 = groups[pair.col_];
				uint32_t group1Size = groupNodes[group1].size();
				uint32_t group2Size = groupNodes[group2].size();
				double sumOfSquaresAll = 0;
				double sumOfSquaresGroup1 = 0;
				double sumOfSquaresGroup2 = 0;
				//group 1
				if(group1Size > 1){
					PairwisePairFactory group1_pFac(group1Size);
					PairwisePairFactory::PairwisePair group1_pair;
					while(group1_pFac.setNextPair(group1_pair)){

						uint32_t node1 = groupNodes[group1][group1_pair.col_];
						uint32_t node2 = groupNodes[group1][group1_pair.row_];
						uint32_t distRow = std::max(node1,  node2);
						uint32_t distCol = std::min(node1,  node2);
						double squareDist = std::pow(dist[distRow][distCol], 2.0);
						sumOfSquaresGroup1 += squareDist;
						sumOfSquaresAll += squareDist;
					}
				}
				//group 2
				if(group2Size > 1){
					PairwisePairFactory group2_pFac(group2Size);
					PairwisePairFactory::PairwisePair group2_pair;
					while(group2_pFac.setNextPair(group2_pair)){
						uint32_t node1 = groupNodes[group2][group2_pair.col_];
						uint32_t node2 = groupNodes[group2][group2_pair.row_];
						uint32_t distRow = std::max(node1,  node2);
						uint32_t distCol = std::min(node1,  node2);
						double squareDist = std::pow(dist[distRow][distCol], 2.0);
						sumOfSquaresGroup2 += squareDist;
						sumOfSquaresAll += squareDist;
					}
				}
				//between group
				for(const auto & group1Node : groupNodes[group1]){
					for(const auto & group2Node : groupNodes[group2]){
						uint32_t distRow = std::max(group1Node,  group2Node);
						uint32_t distCol = std::min(group1Node,  group2Node);
						double squareDist = std::pow(dist[distRow][distCol], 2.0);
            differentInitialGroupDists.emplace_back(dist[distRow][distCol]);
						sumOfSquaresAll += squareDist;
					}
				}
				double distBetweenCentroids = (sumOfSquaresAll - (group1Size + group2Size) * (sumOfSquaresGroup1 / group1Size + sumOfSquaresGroup2 / group2Size)) / (group1Size * group2Size);
				centroidDistances[std::max(group1,  group2)][std::min(group1,  group2)] = distBetweenCentroids;
				if(distBetweenCentroids < lowestCentroidDistInitial) {
					lowestCentroidDistInitial = distBetweenCentroids;
				}
			}
      {
        double differentInitialGroupDists_mean = vectorMean(differentInitialGroupDists);
        double differentInitialGroupDists_sd = vectorStandardDeviationPop(differentInitialGroupDists);
        boost::math::normal_distribution diffDistr(differentInitialGroupDists_mean, differentInitialGroupDists_sd);
        OutOptions tTestsOutOpts(njh::files::make_path(hdbsDir, "initialGroupsTTests.tab.txt"));
        OutputStream tTestsOut(tTestsOutOpts);
        tTestsOut << "group\tgroupSize\tt-statistic\tp-value\tdiffMean\tdiffSD\tgroupMean\tgroupSD" << std::endl;
        for(const auto & group : groups){
          uint32_t groupSize = groupNodes[group].size();
          if(groupSize > 1){
            std::vector<double> sameGroupDist;
            PairwisePairFactory group_pFac(groupSize);
            PairwisePairFactory::PairwisePair group_pair;
            while(group_pFac.setNextPair(group_pair)){

              uint32_t node1 = groupNodes[group][group_pair.col_];
              uint32_t node2 = groupNodes[group][group_pair.row_];
              uint32_t distRow = std::max(node1,  node2);
              uint32_t distCol = std::min(node1,  node2);
              sameGroupDist.emplace_back(dist[distRow][distCol]);

            }

            auto [t, p] = boost::math::statistics::two_sample_t_test(differentInitialGroupDists, sameGroupDist);
            tTestsOut << group << "\t" << groupSize << "\t" << t << "\t" << p
                << "\t" << differentInitialGroupDists_mean << "\t" << differentInitialGroupDists_sd
                << "\t" << vectorMean(sameGroupDist) << "\t" <<  vectorStandardDeviationSamp(sameGroupDist) << std::endl;
          } else {
            tTestsOut << group << "\t" << groupSize << "\t" << "NA" << "\t" << "NA"
                << "\t" << differentInitialGroupDists_mean << "\t" << differentInitialGroupDists_sd
                << "\t" << "NA" << "\t" << "NA"<< std::endl;
          }
        }
      }
			//
			{
				OutOptions initialGroupsOutOpts(njh::files::make_path(hdbsDir, "initialGroupNames.tab.txt"));
				OutputStream initialGroupsOut(initialGroupsOutOpts);
				initialGroupsOut << "name\tgroup" << std::endl;
				for(const auto & n : distGraph.nodes_){
					initialGroupsOut << n->name_ << "\t" << n->group_ << std::endl;
				}
				OutOptions distOutOpts(njh::files::make_path(hdbsDir, "initial_centroidDistances.tab.txt"));
				OutputStream distOut(distOutOpts);
				//write out distance matrix
				uint32_t rowCount = 0;
				printVector(njh::getVecOfMapKeys(groupNodes), "\t", distOut);
				for (const auto & row : centroidDistances) {
					distOut << njh::conToStr(row, "\t") << std::endl;
					++rowCount;
				}
			}
		}


		uint32_t numberOfNonSingletClusters = 0;
		{
			auto groupCounts = distGraph.getGroupCounts();
			for(const auto & count : groupCounts){
				if(HDBSCountSingletGroups || count.second > 1){
					++numberOfNonSingletClusters;
				}
			}
		}
		hdbsRunInfoOut << "0\t" << numberOfNonSingletClusters << "\t" << distGraph.numberOfGroups_ << "\t" << lowestCentroidDistInitial << std::endl;

		uint32_t runCount = 0;
		if(setUp.pars_.verbose_){
			std::cout << "Initial Number of Clusters: " << numberOfNonSingletClusters << std::endl;
			std::cout << "Proposed Clusters: " << proposedClusters << std::endl;
		}

		while(numberOfNonSingletClusters > proposedClusters){
			++runCount;
			if(setUp.pars_.verbose_){
				std::cout << "Run: " << runCount << std::endl;
				std::cout << "Current Number of Clusters: " << numberOfNonSingletClusters << std::endl;
				std::cout << "Proposed Clusters: " << proposedClusters << std::endl;
			}
			double lowestCentroidDist = std::numeric_limits<double>::max();
			std::set<uint32_t> groupsToCollapse;
			//get lowest centroid distances
			std::map<uint32_t, std::vector<uint32_t>> groupNodes;
			for(const auto & node : distGraph.nodes_){
				groupNodes[node->group_].emplace_back(distGraph.nameToNodePos_[node->name_]);
			}
			auto groups = getVectorOfMapKeys(groupNodes);
			PairwisePairFactory pFac(groups.size());
			PairwisePairFactory::PairwisePair pair;
			while(pFac.setNextPair(pair)){
				uint32_t group1 = groups[pair.row_];
				uint32_t group2 = groups[pair.col_];
				double distBetweenCentroids = centroidDistances[std::max(group1,  group2)][std::min(group1,  group2)];
				if(distBetweenCentroids < lowestCentroidDist){
					lowestCentroidDist = distBetweenCentroids;
					groupsToCollapse.clear();
					groupsToCollapse.emplace(groups[pair.row_]);
					groupsToCollapse.emplace(groups[pair.col_]);
				}else if(distBetweenCentroids == lowestCentroidDist){
					groupsToCollapse.emplace(groups[pair.row_]);
					groupsToCollapse.emplace(groups[pair.col_]);
				}
			}
			if(setUp.pars_.verbose_){
				std::cout << "lowestCentroidDist: " << lowestCentroidDist << std::endl;
			}
			//collapse group by setting group value
			for(const auto & group : groupsToCollapse){
				for(const auto & nodeIdx : groupNodes[group]){
					distGraph.nodes_[nodeIdx]->group_ = *groupsToCollapse.begin();
				}
				if(*groupsToCollapse.begin() != group){
					//add this group's nodes to the collapsed to group
					njh::addOtherVec(groupNodes[*groupsToCollapse.begin()], groupNodes[group]);
				}
			}
			numberOfNonSingletClusters = 0;
			{
				auto groupCounts = distGraph.getGroupCounts();
				distGraph.numberOfGroups_ = groupCounts.size();
				for(const auto & count : groupCounts){
					if(HDBSCountSingletGroups || count.second > 1){
						++numberOfNonSingletClusters;
					}
				}
			}
			//re-compute centroid distances now that clusters have been collapsed
			{
				uint32_t modifiedGroup = *groupsToCollapse.begin();
				njh::concurrent::LockableQueue<uint32_t> groupsQueue(groups);
				uint32_t otherGroup = std::numeric_limits<uint32_t>::max();
        //std::vector<double> allCurrentDists;
				while(groupsQueue.getVal(otherGroup)){
					if(!njh::in(otherGroup, groupsToCollapse)){
						uint32_t group1 = modifiedGroup;
						uint32_t group2 = otherGroup;
						uint32_t group1Size = groupNodes[group1].size();
						uint32_t group2Size = groupNodes[group2].size();
						double sumOfSquaresAll = 0;
						double sumOfSquaresGroup1 = 0;
						double sumOfSquaresGroup2 = 0;

						std::mutex sumsMut;
						//group 1
						if(group1Size > 1){
							PairwisePairFactory group1_pFac(group1Size);
							std::function<void()> computeSumsOfSqaures = [&group1_pFac,&sumsMut,&sumOfSquaresGroup1,&sumOfSquaresAll,&group1,&dist,&groupNodes](){
								PairwisePairFactory::PairwisePair group1_pair;
								while(group1_pFac.setNextPair(group1_pair)){
									uint32_t node1 = groupNodes[group1][group1_pair.col_];
									uint32_t node2 = groupNodes[group1][group1_pair.row_];
									uint32_t distRow = std::max(node1,  node2);
									uint32_t distCol = std::min(node1,  node2);
									double squareDist = std::pow(dist[distRow][distCol], 2.0);
									//allCurrentDists.emplace_back(dist[distRow][distCol]);
									{
										std::lock_guard<std::mutex> lock(sumsMut);
										sumOfSquaresGroup1 += squareDist;
										sumOfSquaresAll += squareDist;
									}
								}
							};
							njh::concurrent::runVoidFunctionThreaded(computeSumsOfSqaures, numThreads);
						}
						//group 2
						if(group2Size > 1){
							PairwisePairFactory group2_pFac(group2Size);
							std::function<void()> computeSumsOfSqaures = [&group2_pFac,&sumsMut,&sumOfSquaresGroup2,&sumOfSquaresAll,&group2,&dist,&groupNodes](){
								PairwisePairFactory::PairwisePair group2_pair;
								while(group2_pFac.setNextPair(group2_pair)){
									uint32_t node1 = groupNodes[group2][group2_pair.col_];
									uint32_t node2 = groupNodes[group2][group2_pair.row_];
									uint32_t distRow = std::max(node1,  node2);
									uint32_t distCol = std::min(node1,  node2);
									double squareDist = std::pow(dist[distRow][distCol], 2.0);
									//allCurrentDists.emplace_back(dist[distRow][distCol]);
									{
										std::lock_guard<std::mutex> lock(sumsMut);
										sumOfSquaresGroup2 += squareDist;
										sumOfSquaresAll += squareDist;
									}
								}
							};
							njh::concurrent::runVoidFunctionThreaded(computeSumsOfSqaures, numThreads);
						}
						//between group
						AllByAllPairFactory allFac(groupNodes[group1].size(), groupNodes[group2].size());
						std::function<void()> computeSumsOfSqaures = [&allFac,&sumsMut,&sumOfSquaresAll,&dist,&groupNodes,&group1,&group2](){
							AllByAllPairFactory::AllByAllPair allPair;
							while(allFac.setNextPair(allPair)){
								auto group1Node = groupNodes[group1][allPair.row_];
								auto group2Node = groupNodes[group2][allPair.col_];
								uint32_t distRow = std::max(group1Node,  group2Node);
								uint32_t distCol = std::min(group1Node,  group2Node);
								double squareDist = std::pow(dist[distRow][distCol], 2.0);
								//allCurrentDists.emplace_back(dist[distRow][distCol]);
								{
									std::lock_guard<std::mutex> lock(sumsMut);
									sumOfSquaresAll += squareDist;
								}
							}
						};
						njh::concurrent::runVoidFunctionThreaded(computeSumsOfSqaures, numThreads);
						double distBetweenCentroids = (sumOfSquaresAll - (group1Size + group2Size) * (sumOfSquaresGroup1 / group1Size + sumOfSquaresGroup2 / group2Size)) / (group1Size * group2Size);
						centroidDistances[std::max(group1,  group2)][std::min(group1,  group2)] = distBetweenCentroids;
					}
				}
//        {
//          auto [t, p] = boost::math::statistics::two_sample_t_test(differentInitialGroupDists, allCurrentDists);
//          if(setUp.pars_.verbose_){
//            auto diffDists_mean = vectorMean(differentInitialGroupDists);
//            auto diffDists_sd = vectorStandardDeviationSamp(differentInitialGroupDists);
//            boost::math::normal_distribution diffDistr(diffDists_mean, diffDists_sd);
//            auto allCurrentDists_mean = vectorMean(allCurrentDists);
//            auto allCurrentDists_sd = vectorStandardDeviationSamp(allCurrentDists);
//            boost::math::normal_distribution currentDistr(allCurrentDists_mean, allCurrentDists_sd);
//
//            std::cout << "differentInitialGroupDists.size(): " << differentInitialGroupDists.size() << std::endl;
//            std::cout << "allCurrentDists.size(): " << allCurrentDists.size() << std::endl;
//
//            std::cout << "\t" << vectorMean(differentInitialGroupDists) << "\t" <<  vectorStandardDeviationSamp(differentInitialGroupDists) << std::endl
//                    << "\t" << vectorMean(allCurrentDists) << "\t" <<  vectorStandardDeviationSamp(allCurrentDists) << std::endl;
//            std::cout << "t-statistic: " << t << ", p-value: " << p << ": " << njh::colorBool(p < 0.01) << std::endl;
//            std::cout << "in different overlap: " << cdf(diffDistr,allCurrentDists_mean + allCurrentDists_sd*2) << ": " << njh::colorBool(cdf(diffDistr,allCurrentDists_mean + allCurrentDists_sd*2) < 0.01) << std::endl;
//            std::cout << "in current overlap  : " << 1 - cdf(currentDistr, diffDists_mean - 2 * diffDists_sd) << std::endl;
//          }
//        }
			}

			hdbsRunInfoOut << runCount << "\t" << numberOfNonSingletClusters << "\t" << distGraph.numberOfGroups_ << "\t" << lowestCentroidDist << std::endl;
		}
		//final centroid
		{
			//compute centroid distances
			std::map<uint32_t, std::vector<uint32_t>> groupNodes;
			for(const auto & node : distGraph.nodes_){
				groupNodes[node->group_].emplace_back(distGraph.nameToNodePos_[node->name_]);
			}
			auto groups = getVectorOfMapKeys(groupNodes);
			//
			{
				OutOptions distOutOpts(njh::files::make_path(hdbsDir, "final_centroidDistances.tab.txt"));
				OutputStream distOut(distOutOpts);
				//write out distance matrix
				printVector(njh::getVecOfMapKeys(groupNodes), "\t", distOut);
				for(const auto & group1Idx : iter::range(groups.size())){
					std::vector<double> outVec;
					uint32_t group1 = groups[group1Idx];
					for(const auto & group2Idx : iter::range(0UL, group1Idx + 1)){
						uint32_t group2 = groups[group2Idx];
						outVec.emplace_back(centroidDistances[std::max(group1,  group2)][std::min(group1,  group2)] );
					}
					distOut << njh::conToStr(outVec, "\t") << std::endl;
				}
			}
		}
	} else {
		distGraph.dbscan(dbPars_);
		distGraph.assignNoiseNodesAGroup();
	}



	{
		OutOptions distOutOpts(njh::files::make_path(setUp.pars_.directoryName_, "distance.tab.txt"));
		OutputStream distOut(distOutOpts);
		//write out distance matrix
		uint32_t rowCount = 0;
		printVector(readVec::getNames(reads), "\t", distOut);
		for (const auto & row : dist) {
			distOut << njh::conToStr(row, "\t");
			if(rowCount > 0){
				distOut << "\t";
			}
			distOut << 0 << std::endl;
			++rowCount;
		}
	}

	//write out group counts;

	{
		//write out group counts;
		OutOptions parametersOutOpts(njh::files::make_path(setUp.pars_.directoryName_, "parameters.tab.txt"));
		OutputStream parametersOut(parametersOutOpts);
		parametersOut << "parameter\tvalue" << std::endl;
		parametersOut << "klen" << "\t" << kmerLength << std::endl;
		parametersOut << "epsilon" << "\t" << dbPars_.eps_ << std::endl;
		parametersOut << "minEpNeighbors" << "\t" << dbPars_.minEpNeighbors_ << std::endl;
	}
	{
		//write out group counts;
		OutOptions groupCountOutOpts(njh::files::make_path(setUp.pars_.directoryName_, "groupCounts.tab.txt"));
		OutputStream groupCountOut(groupCountOutOpts);
		std::map<uint32_t, uint32_t> groupCounts;
		for(const auto & n : distGraph.nodes_){
			++groupCounts[n->group_];
		}
		groupCountOut << "group\t" << "count\t" << "percent" << std::endl;
		for(const auto & count : groupCounts){
			groupCountOut << count.first
			<< "\t" << count.second
			<< "\t" << roundDecPlaces(static_cast<double>(count.second)/reads.size() * 100, 2) << std::endl;
		}
	}

	{
		OutOptions errorOutOpts(njh::files::make_path(setUp.pars_.directoryName_, "variance.tab.txt"));
		OutputStream errorOut(errorOutOpts);
		errorOut << "group\tnum\tmeanDiff\tvariance" << std::endl;

		std::map<uint32_t, std::vector<uint32_t>> groupNodes;
		for(const auto & node : distGraph.nodes_){
			groupNodes[node->group_].emplace_back(distGraph.nameToNodePos_[node->name_]);
		}
		for(const auto & group : groupNodes){
			if(group.second.size() < dbPars_.minEpNeighbors_){
				continue;
			}
			PairwisePairFactory pFactor(group.second.size());
			PairwisePairFactory::PairwisePair pair;
			std::vector<double> distancesWithinGroup;
			while(pFactor.setNextPair(pair)){
				uint32_t node1 = group.second[pair.col_];
				uint32_t node2 = group.second[pair.row_];
				uint32_t distRow = std::max(node1,  node2);
				uint32_t distCol = std::min(node1,  node2);
				distancesWithinGroup.emplace_back(dist[distRow][distCol]);
			}

			auto meanDist = vectorMean(distancesWithinGroup);
			auto var = vectorStandardDeviationPop(distancesWithinGroup);
			errorOut << group.first
			<< "\t" << group.second.size()
			<< "\t" << meanDist
			<< "\t" << var
			<< std::endl;
		}
	}

	{
		//write out groups
		MultiSeqIO seqWriter;
		std::map<uint32_t, uint32_t> groupCounts = distGraph.getGroupCounts();
		auto maxGroup = vectorMaximum(getVectorOfMapKeys(groupCounts));
		for (const auto group: njh::getVecOfMapKeys(groupCounts)) {
			seqWriter.addReader(njh::pasteAsStr(group), SeqIOOptions(
							njh::files::make_path(setUp.pars_.directoryName_, njh::leftPadNumStr<uint32_t>(group, maxGroup)),
							setUp.pars_.ioOptions_.outFormat_));
		}
		for(const auto & n : distGraph.nodes_){
			seqWriter.openWrite(estd::to_string(n->group_), n->value_);
		}
	}

	{
		OutOptions groupNamesOutOpts(njh::files::make_path(setUp.pars_.directoryName_, "groupNames.tab.txt"));
		OutputStream groupNamesOut(groupNamesOutOpts);
		groupNamesOut << "name\tgroup" << std::endl;
		for(const auto & n : distGraph.nodes_){
			groupNamesOut << n->name_ << "\t" << n->group_ << std::endl;
		}
	}

	{
		//write heatmap dist report
		OutOptions distGraphOutOpts(njh::files::make_path(setUp.pars_.directoryName_, "plotDistHeatmap.Rmd"));
		OutputStream distGraphOut(distGraphOutOpts);
		distGraphOut << distGraphRmdReport;
	}


	{
		auto subGroupVisDir = njh::files::make_path(setUp.pars_.directoryName_, "visulation");
		njh::files::makeDir(njh::files::MkdirPar{subGroupVisDir});
		VecStr groups = toVecStr(getVectorOfMapKeys(distGraph.getGroupCounts()));
		auto colors = getColorsForNames(groups);
		std::unordered_map<std::string, std::string> nameToColor;
		for(const auto & n : distGraph.nodes_){
			nameToColor[n->value_->name_] = colors[estd::to_string(n->group_)].getHexStr();
		}
		auto treeJson = distGraph.toJson(0, nameToColor);
		auto & nodes = treeJson["nodes"];
		for(auto & node : nodes){
			node["size"] = distGraph.nodes_[distGraph.nameToNodePos_[node["name"].asString()]]->value_->cnt_ * 50;
		}
		OutputStream treeJsonFile(njh::files::make_path(subGroupVisDir, "tree.json"));
		treeJsonFile << treeJson;

		OutputStream treeHtmlFile(njh::files::make_path(subGroupVisDir,  "tree.html"));
		genTreeHtml(treeHtmlFile, "tree.json", "tree.js");

		OutputStream treeJsFile(njh::files::make_path(subGroupVisDir, "tree.js"));
		genSimpleTreeJs(treeJsFile);
	}

	if(setUp.pars_.verbose_	){
		watch.logLapTimes(std::cout, true, 6, true);
	}

	return 0;
}


}  // namespace njhseq
