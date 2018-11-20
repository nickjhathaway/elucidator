/*
 * pacbioExp_sim.cpp
 *
 *  Created on: Oct 24, 2016
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

#include "pacbioExp.hpp"
#include "elucidator/simulation.h"

namespace njhseq {


int pacbioExpRunner::simPacbioPerRead(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path roundLikelihoodsFnp = "";
	bfs::path errorQualLikelihoodsFnp = "";
	bfs::path regQualLikelihoodsFnp = "";
	bfs::path actionLikelihoodsFnp = "";
	bfs::path deletionSizeLikelihoodsFnp = "";
	bfs::path insertionSizeLikelihoodsFnp = "";
	bfs::path baseMutLikelihoodsFnp = "";
	bfs::path profileDir = "";

	seqSetUp setUp(inputCommands);
	setUp.setOption(profileDir, "--profileDir", "Directory with error profile", true);
/*
	setUp.setOption(roundLikelihoodsFnp, "--roundLikelihoods", "Table with the likelihoods for number of rounds", true);
	setUp.setOption(errorQualLikelihoodsFnp, "--errorQualLikelihoodsFnp", "Table with the likelihoods for quality scores of the pacbio errors", true);
	setUp.setOption(regQualLikelihoodsFnp, "--regQualLikelihoodsFnp", "Table with the likelihoods for quality scores of the pacbio regular bases", true);
	setUp.setOption(actionLikelihoodsFnp, "--actionLikelihoodsFnp", "Table with the likelihoods for likelihood of match, snp, insertion, deletion", true);
	setUp.setOption(deletionSizeLikelihoodsFnp, "--deletionSizeLikelihoodsFnp", "Table with the likelihoods for size of deletions", true);
	setUp.setOption(insertionSizeLikelihoodsFnp, "--insertionSizeLikelihoodsFnp", "Table with the likelihoods for size of insertions", true);
	setUp.setOption(baseMutLikelihoodsFnp, "--baseMutLikelihoodsFnp", "Table with the likelihoods for base mutation rates", true);*/
	setUp.processDefaultReader(true);
	setUp.pars_.ioOptions_.outFormat_ = SeqIOOptions::outFormats::FASTQ;
	setUp.pars_.ioOptions_.out_.outExtention_ = ".fastq";
	if("out" == setUp.pars_.ioOptions_.out_.outFilename_){
		setUp.pars_.ioOptions_.out_.outFilename_ = bfs::path(setUp.pars_.ioOptions_.firstName_).replace_extension("").string() + "_simPacbio";
	}
	setUp.finishSetUp(std::cout);

	roundLikelihoodsFnp = njh::files::make_path(profileDir,
			"roundDist.tab.txt");
	errorQualLikelihoodsFnp = njh::files::make_path(profileDir,
			"errorQualDist.tab.txt");
	regQualLikelihoodsFnp = njh::files::make_path(profileDir,
			"regQualDist.tab.txt");
	actionLikelihoodsFnp = njh::files::make_path(profileDir,
			"baseActionDist.tab.txt");
	deletionSizeLikelihoodsFnp = njh::files::make_path(profileDir,
			"deletionSizeDist.tab.txt");
	insertionSizeLikelihoodsFnp = njh::files::make_path(profileDir,
			"insertionSizeDist.tab.txt");
	baseMutLikelihoodsFnp = njh::files::make_path(profileDir,
			"baseMutDist.tab.txt");


	//set up round likelihood
	table roundLikelihoodTab(roundLikelihoodsFnp.string(), "\t", true);
	std::vector<uint32_t> rounds  = njh::lexical_cast_con<VecStr, std::vector<uint32_t>>(roundLikelihoodTab.getColumn("rounds"));
	std::vector<double> roundLikelihoods = njh::lexical_cast_con<VecStr, std::vector<double>>(roundLikelihoodTab.getColumn("frac"));
	njh::randObjectGen <uint32_t, double> roundGenerator(rounds, roundLikelihoods);
	//set up qual error likelihoods
	table errorQualLikelihoodTab(errorQualLikelihoodsFnp.string(), "\t", true);
	auto errorQualSplits = errorQualLikelihoodTab.splitTableOnColumn("rounds");
	std::unordered_map<uint32_t, njh::randObjectGen <uint32_t, double> > errorQualGens;
	for(const auto & errorQual : errorQualSplits){
		std::vector<uint32_t> quals  = njh::lexical_cast_con<VecStr, std::vector<uint32_t>>(errorQual.second.getColumn("qual"));
		std::vector<double> qualsLikelihoods = njh::lexical_cast_con<VecStr, std::vector<double>>(errorQual.second.getColumn("rate"));
		errorQualGens.emplace(estd::stou(errorQual.first), njh::randObjectGen <uint32_t, double> (quals, qualsLikelihoods));
	}
	//set up qual non error likelihoods
	table regQualLikelihoodTab(regQualLikelihoodsFnp.string(), "\t", true);
	auto regQualSplits = regQualLikelihoodTab.splitTableOnColumn("rounds");
	std::unordered_map<uint32_t, njh::randObjectGen <uint32_t, double> > regQualGens;
	for(const auto & regQual : regQualSplits){
		std::vector<uint32_t> quals  = njh::lexical_cast_con<VecStr, std::vector<uint32_t>>(regQual.second.getColumn("qual"));
		std::vector<double> qualsLikelihoods = njh::lexical_cast_con<VecStr, std::vector<double>>(regQual.second.getColumn("rate"));
		regQualGens.emplace(estd::stou(regQual.first), njh::randObjectGen <uint32_t, double> (quals, qualsLikelihoods));
	}
	//set up base mut likelihoods
	table baseMutLikelihoodTab(baseMutLikelihoodsFnp.string(), "\t", true);
	auto baseMutSplits = baseMutLikelihoodTab.splitTableOnColumn("rounds");
	std::unordered_map<uint32_t, std::unordered_map<char, njh::randObjectGen <char, double>>> baseMutGen;
	for(const auto & baseMut : baseMutSplits){
		auto byBaseSplits = baseMut.second.splitTableOnColumn("refBase");
		std::unordered_map<char, njh::randObjectGen <char, double> > baseGen;
		for(const auto & base : byBaseSplits){
			VecStr basebaseMutsStr = base.second.getColumn("error");
			std::vector<char> basebaseMuts;
			for(const auto & b : basebaseMutsStr){
				basebaseMuts.emplace_back(b.front());
			}
			std::vector<double> basebaseMutsLikelihoods = njh::lexical_cast_con<VecStr, std::vector<double>>(base.second.getColumn("rate"));
			baseGen.emplace(base.first.front(), njh::randObjectGen <char, double> (basebaseMuts, basebaseMutsLikelihoods));
		}
		baseMutGen.emplace(estd::stou(baseMut.first), baseGen);
	}

	//set up insertion generator
	charCounter baseCounting{std::vector<char>{'A', 'C', 'G', 'T'}};
	{
		SeqIO reader(setUp.pars_.ioOptions_);
		reader.openIn();
		seqInfo seq;
		while(reader.readNextRead(seq)){
			baseCounting.increaseCountByString(seq.seq_);
		}
	}

	std::unordered_map<char, njh::randObjectGen<char, uint32_t>> charGen;
	charGen.emplace('T', njh::randObjectGen<char, uint32_t>(std::vector<char>{'A', 'C', 'G', 'T'},
			std::vector<uint32_t>{baseCounting.chars_['A'],baseCounting.chars_['C'],baseCounting.chars_['G'], baseCounting.chars_['T']*5}));
	charGen.emplace('C', njh::randObjectGen<char, uint32_t>(std::vector<char>{'A', 'C', 'G', 'T'},
			std::vector<uint32_t>{baseCounting.chars_['A'],baseCounting.chars_['C']*5,baseCounting.chars_['G'], baseCounting.chars_['T']}));
	charGen.emplace('G', njh::randObjectGen<char, uint32_t>(std::vector<char>{'A', 'C', 'G', 'T'},
			std::vector<uint32_t>{baseCounting.chars_['A'],baseCounting.chars_['C'],baseCounting.chars_['G']*5, baseCounting.chars_['T']}));
	charGen.emplace('A', njh::randObjectGen<char, uint32_t>(std::vector<char>{'A', 'C', 'G', 'T'},
			std::vector<uint32_t>{baseCounting.chars_['A']*5,baseCounting.chars_['C'],baseCounting.chars_['G'], baseCounting.chars_['T']}));
	//set up insertion size distribution
	table insertionSizeLikelihoodTab(insertionSizeLikelihoodsFnp.string(), "\t", true);
	auto insertionSizeSplits = insertionSizeLikelihoodTab.splitTableOnColumn("rounds");
	std::unordered_map<uint32_t, njh::randObjectGen <uint32_t, double> > insertionSizeGens;
	for(const auto & insertionSize : insertionSizeSplits){
		std::vector<uint32_t> insertionSizes  = njh::lexical_cast_con<VecStr, std::vector<uint32_t>>(insertionSize.second.getColumn("errorSize"));
		std::vector<double> insertionSizesLikelihoods = njh::lexical_cast_con<VecStr, std::vector<double>>(insertionSize.second.getColumn("rate"));
		insertionSizeGens.emplace(estd::stou(insertionSize.first), njh::randObjectGen <uint32_t, double> (insertionSizes, insertionSizesLikelihoods));
	}
	//set up deletion size distribution
	table deletionSizeLikelihoodTab(deletionSizeLikelihoodsFnp.string(), "\t", true);
	auto deletionSizeSplits = deletionSizeLikelihoodTab.splitTableOnColumn("rounds");
	std::unordered_map<uint32_t, njh::randObjectGen <uint32_t, double> > deletionSizeGens;
	for(const auto & deletionSize : deletionSizeSplits){
		std::vector<uint32_t> deletionSizes  = njh::lexical_cast_con<VecStr, std::vector<uint32_t>>(deletionSize.second.getColumn("errorSize"));
		std::vector<double> deletionSizesLikelihoods = njh::lexical_cast_con<VecStr, std::vector<double>>(deletionSize.second.getColumn("rate"));
		deletionSizeGens.emplace(estd::stou(deletionSize.first), njh::randObjectGen <uint32_t, double> (deletionSizes, deletionSizesLikelihoods));
	}
	//set up action distribution
	table actionLikelihoodTab(actionLikelihoodsFnp.string(), "\t", true);
	auto actionSplits = actionLikelihoodTab.splitTableOnColumn("rounds");
	std::unordered_map<uint32_t, std::unordered_map<char, njh::randObjectGen <std::string, double>>> actionGen;
	for(const auto & action : actionSplits){
		auto byBaseSplits = action.second.splitTableOnColumn("base");
		std::unordered_map<char, njh::randObjectGen <std::string, double> > baseGen;
		for(const auto & base : byBaseSplits){
			VecStr baseActions  = base.second.getColumn("action");
			std::vector<double> baseActionsLikelihoods = njh::lexical_cast_con<VecStr, std::vector<double>>(base.second.getColumn("rate"));
			baseGen.emplace(base.first.front(), njh::randObjectGen <std::string, double> (baseActions, baseActionsLikelihoods));
		}
		actionGen.emplace(estd::stou(action.first), baseGen);
	}

	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	reader.openOut();
	seqInfo seq;
	while(reader.readNextRead(seq)){
		auto round = roundGenerator.genObj();
		seqInfo outSeq(seq.name_);
		for(uint32_t basePos = 0; basePos < len(seq); ++basePos){
			char inputBase = seq.seq_[basePos];
			auto action = actionGen.at(round).at(inputBase).genObj();
			if("deletion" == action){
				//deletion, don't add base add size of deletion minus one to properly move the basePos
				//to delete the size of the deletion
				auto deletionSize = deletionSizeGens.at(round).genObj();
				basePos += deletionSize - 1;
			}else if("insertion" == action){
				//add the base and then insert
				outSeq.append(inputBase, regQualGens.at(round).genObj());
				auto insertSize = insertionSizeGens.at(round).genObj();
				for(uint32_t i = 0; i < insertSize; ++i){
					outSeq.append(charGen.at(inputBase).genObj(), errorQualGens.at(round).genObj());
				}
			}else if("match" == action){
				//just add regular base and add a reg qual
				outSeq.append(inputBase, regQualGens.at(round).genObj());
			}else if("snp" == action){
				//mutate the base and add errror qual
				outSeq.append(baseMutGen.at(round).at(inputBase).genObj(), errorQualGens.at(round).genObj());
			}else{
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ <<  ":  Error, unhandled action: " << action << "\n";
				throw std::runtime_error{ss.str()};
			}
		}
		reader.write(outSeq);
	}


	return 0;
}

} //namespace njhseq

