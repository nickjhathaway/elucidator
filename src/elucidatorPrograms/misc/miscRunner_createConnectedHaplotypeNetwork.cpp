/*
 * miscRunner_createConnectedHaplotypeNetwork.cpp
 *
 *  Created on: Feb 6, 2018
 *      Author: nick
 */
#include "miscRunner.hpp"
#include "elucidator/simulation.h"
#include "elucidator/BamToolsUtils.h"
#include "elucidator/objects/seqContainers.h"
#include "elucidator/objects/BioDataObject.h"
#include "elucidator/utils/KmerUtils.hpp"
#include "elucidator/objects/dataContainers.h"

namespace njhseq {





int miscRunner::createConnectedHaplotypeNetwork(const njh::progutils::CmdArgs & inputCommands){
	bfs::path metafnp = "";
	bfs::path fieldColorFnp = "";
	ReadCompGraph::ConnectedHaplotypeNetworkPars netPars;
	seqSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processVerbose();
	setUp.setOption(netPars.noLabel, "--noLabel", "noLabel");
	setUp.setOption(netPars.matchPars.kmerLen_, "--compKmerLength", "Length of the k-mer to be used to skip comparisons");
	setUp.setOption(netPars.matchPars.kmerCutOff_, "--kmerCutOff", "k-mer similarity cut off to be used to skip comparisons");
	setUp.setOption(netPars.numThreads, "--numThreads", "Number of threads to use");
	setUp.setOption(netPars.minNumberOfEvents, "--minNumberOfDifferences", "Minimum number of differences to allow for connections ");
	setUp.setOption(netPars.minId,      "--minId", "minimum percent identity to for connections");
	setUp.setOption(netPars.setJustBest,"--setJustBest", "set just best connections");
	setUp.setOption(netPars.doTies,     "--doTies", "Allow ties while do doing setJustBest");
	setUp.setOption(metafnp,            "--metafnp", "meta file for adding additions to nodes");
	setUp.setOption(fieldColorFnp,      "--fieldColorFnp", "A table with at least two columns, the field name of the colorField and color, the hex color for the nodes");
	setUp.setOption(netPars.colorField, "--colorField", "The field to color by");
	setUp.setOption(netPars.labelField, "--labelField", "The field to add a label to the nodes for");
	setUp.processReadInNames();
	setUp.processAlignerDefualts();
	setUp.processDirectoryOutputName(true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);

	if ("" != metafnp) {
		netPars.loadMeta(metafnp);
		bfs::copy(bfs::absolute(metafnp), njh::files::make_path(setUp.pars_.directoryName_, "metaForSeqs.tab.txt"));
	}
	netPars.checkForColorField(__PRETTY_FUNCTION__);
	auto colorLookup = netPars.generateColorLookup(fieldColorFnp);
	if("" != netPars.colorField){
		OutputStream colorOut(OutOptions{njh::files::make_path(setUp.pars_.directoryName_, "fieldColorValues.tab.txt")});
		colorOut << netPars.colorField<< "\t" << "color" << std::endl;
		auto cKeys = njh::getVecOfMapKeys(colorLookup);
		njh::sort(cKeys);
		for(const auto & c : cKeys){
			colorOut << c
					<< "\t" << colorLookup.at(c).getHexStr() << std::endl;
		}
	}
	netPars.checkForLabelField(__PRETTY_FUNCTION__);


	OutOptions jsonOutOpts(njh::files::make_path(setUp.pars_.directoryName_, "tree"), ".json");
	OutputStream jsonOut(jsonOutOpts);
	OutOptions htmlOutOpts(njh::files::make_path(setUp.pars_.directoryName_, "tree"), ".html");
	OutputStream htmlOut(htmlOutOpts);

	htmlOut << ReadCompGraph::ConnectedHaplotypeNetworkPars::htmlPageForConnectHpaNet  << std::endl;

	uint64_t maxlen = 0;
	auto inReads = SeqInput::getSeqVec<seqInfo>(setUp.pars_.ioOptions_, maxlen);

	aligner alignerObj(maxlen, setUp.pars_.gapInfo_, setUp.pars_.scoring_,setUp.pars_.colOpts_.alignOpts_.countEndGaps_);
	concurrent::AlignerPool alnPool(alignerObj, netPars.numThreads);
	alnPool.inAlnDir_ = setUp.pars_.alnInfoDirName_;
	alnPool.outAlnDir_ = setUp.pars_.outAlnInfoDirName_;
	alnPool.initAligners();
	netPars.verbose = setUp.pars_.verbose_;

	ReadCompGraph graph(inReads);
	graph.addEdgesBasedOnIdOrMinDif(netPars, alnPool);
	Json::Value outputJson = graph.getSingleLineJsonOut(netPars, colorLookup);
	graph.writeAdjListPerId(
			OutOptions { njh::files::make_path(setUp.pars_.directoryName_, "links"),
					".txt" }, netPars);
	jsonOut << outputJson << std::endl;

	return 0;
}

}  // namespace njhseq
