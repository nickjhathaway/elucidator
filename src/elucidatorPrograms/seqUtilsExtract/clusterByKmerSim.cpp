//
// Created by Nicholas Hathaway on 4/9/22.
//

#include "seqUtilsExtractRunner.hpp"

#include "elucidator/utils/KmerUtils.hpp"
#include "elucidator/objects/seqObjects/seqKmers.h"
#include <njhseq/IO/SeqIO.h>

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
	readDistGraph<double>::dbscanPars dbPars_;
	dbPars_.minEpNeighbors_ = 2;
	dbPars_.eps_ = .20;


	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processReadInNames(true);
	setUp.setOption(numThreads, "--numThreads", "numThreads");
	setUp.setOption(kmerLength, "--kmerLength", "kmer Length");

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
	distGraph.dbscan(dbPars_);
	distGraph.assignNoiseNodesAGroup();
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
		//write out groups
		MultiSeqIO seqWriter;
		for (const auto group: iter::range(0U, distGraph.numberOfGroups_)) {
			seqWriter.addReader(njh::pasteAsStr(group), SeqIOOptions(
							njh::files::make_path(setUp.pars_.directoryName_, njh::leftPadNumStr(group, distGraph.numberOfGroups_)),
							setUp.pars_.ioOptions_.outFormat_));
		}
		for(const auto & n : distGraph.nodes_){
			seqWriter.openWrite(estd::to_string(n->group_), n->value_);
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
		VecStr groups;
		for(const auto group : iter::range(distGraph.numberOfGroups_)){
			groups.emplace_back(estd::to_string(group));
		}
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
