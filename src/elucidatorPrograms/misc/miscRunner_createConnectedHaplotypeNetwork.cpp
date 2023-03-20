/*
 * miscRunner_createConnectedHaplotypeNetwork.cpp
 *
 *  Created on: Feb 6, 2018
 *      Author: nick
 */

// elucidator - A library for analyzing sequence data
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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
#include "miscRunner.hpp"
#include "elucidator/simulation.h"
#include "elucidator/BamToolsUtils.h"
#include "elucidator/objects/seqContainers.h"
#include "elucidator/objects/BioDataObject.h"
#include <njhseq/objects/kmer/KmerUtils.hpp>
#include "elucidator/objects/dataContainers.h"


#include <njhseq/objects/Meta.h>
#include <njhseq/objects/dataContainers/graphs/ReadCompGraph.hpp>

namespace njhseq {





int miscRunner::createConnectedHaplotypeNetwork(const njh::progutils::CmdArgs & inputCommands){
	bfs::path metafnp = "";
	bfs::path fieldColorFnp = "";
	std::string countField = "";
	ReadCompGraph::ConnectedHaplotypeNetworkPars netPars;
	seqSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processVerbose();
	setUp.setOption(netPars.noNodeLabel, "--noNodeLabel", "Put No Node Labels");
	setUp.setOption(netPars.noLinkLabel, "--noLinkLabel", "Put No Link Labels");

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
	setUp.setOption(countField, "--countField", "The field to use to determine the size of the node");

	setUp.processReadInNames();
	setUp.processAlignerDefualts();
	setUp.processDirectoryOutputName(true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);

	uint64_t maxlen = 0;
	auto inReads = SeqInput::getSeqVec<seqInfo>(setUp.pars_.ioOptions_, maxlen);

	if ("" != metafnp) {
		netPars.loadMeta(metafnp);
		bfs::copy_file(bfs::absolute(metafnp),
				njh::files::make_path(setUp.pars_.directoryName_,
						"metaForSeqs.tab.txt"));
	} else {
		if (!inReads.empty()
				&& MetaDataInName::nameHasMetaData(inReads.front().name_)) {
			auto metaInSeqs = seqsToMetaTable(inReads);
			metaInSeqs.columnNames_.front() = "sample";
			metaInSeqs.deleteColumn("seq");
			auto metaTabOutOpts = TableIOOpts::genTabFileOut(
					njh::files::make_path(setUp.pars_.directoryName_,
							"metaForSeqs.tab.txt"), true);
			metaInSeqs.outPutContents(metaTabOutOpts);
			netPars.loadMeta(metaTabOutOpts.out_.outName());
		}
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

	if("" != countField){
		if(nullptr == netPars.seqMeta){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", supplied a color field but no meta data" << "\n";
			throw std::runtime_error{ss.str()};
		}
		netPars.seqMeta->checkForFieldsThrow({countField});
		std::unordered_map<std::string, uint32_t > nameToPosition;
		uint32_t seqPos = 0;
		for(const auto & seq : inReads){
			nameToPosition[seq.name_] = seqPos;
			++seqPos;
		}
		double total = 0;
		for (const auto & samp : netPars.seqMeta->samples_) {
			inReads[nameToPosition[samp]].cnt_ =
					njh::StrToNumConverter::stoToNum<double>(
							njh::mapAt(netPars.seqMeta->groupData_, countField)->getGroupForSample(
									samp));
			total += inReads[nameToPosition[samp]].cnt_;
		}
		for (auto & seq : inReads) {
			seq.frac_ = seq.cnt_ / total;
		}
	}

	OutOptions jsonOutOpts(njh::files::make_path(setUp.pars_.directoryName_, "tree"), ".json");
	OutputStream jsonOut(jsonOutOpts);
	OutOptions htmlOutOpts(njh::files::make_path(setUp.pars_.directoryName_, "tree"), ".html");
	OutputStream htmlOut(htmlOutOpts);

	htmlOut << ReadCompGraph::ConnectedHaplotypeNetworkPars::htmlPageForConnectHpaNet  << std::endl;



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
