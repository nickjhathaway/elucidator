/*
 * miscRunner_calculateShannonEntropySimple.cpp
 *
 *  Created on: Sep 4, 2017
 *      Author: nicholashathaway
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


namespace njhseq {



int miscRunner::codeComparison(const njh::progutils::CmdArgs & inputCommands){
	OutOptions outOpts(bfs::path(""), ".tab.txt");
	seqSetUp setUp(inputCommands);
	setUp.processReadInNames(true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto allSeqs = reader.readAllReads<seqInfo>();
	if(0 == allSeqs.size()  || allSeqs.size() % 2 != 0){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error expected an even number of sequences" << "\n";
		throw std::runtime_error{ss.str()};
	}
	OutputStream codeOut(outOpts);

	std::unordered_map<std::string, std::vector<std::pair<std::string,std::vector<std::pair<uint32_t, uint32_t>>>>> codes;
	for(const auto pos : iter::range<uint32_t>(0, len(allSeqs) ,2)){
		if(len(allSeqs[pos]) != len(allSeqs[pos + 1])){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error sequences should be the same size, for seq: " << pos << "\n";
			ss << allSeqs[pos].name_ << " length: " << len(allSeqs[pos]) << "\n";
			ss << allSeqs[pos + 1].name_ << " length: " << len(allSeqs[pos + 1]) << "\n";
			throw std::runtime_error{ss.str()};
		}
		uint32_t firstNonGap = std::max(allSeqs[pos].seq_.find_first_not_of('-'), allSeqs[pos + 1].seq_.find_first_not_of('-'));
		uint32_t lastNonGap = std::min(allSeqs[pos].seq_.find_last_not_of('-'), allSeqs[pos + 1].seq_.find_last_not_of('-'));
		if(lastNonGap > firstNonGap){
			std::vector<std::pair<uint32_t, uint32_t>> codedRegions;
			for(const auto seqPos  : iter::range(firstNonGap, lastNonGap + 1)){
				if(allSeqs[pos].seq_[seqPos] != allSeqs[pos + 1].seq_[seqPos]){
					codedRegions.emplace_back(seqPos, 0);
				}else{
					codedRegions.emplace_back(seqPos, 1);
				}
			}
			codes[allSeqs[pos].name_].emplace_back(std::make_pair(allSeqs[pos + 1].name_, codedRegions));
		}
	}
	codeOut << "refName\tseqName\tpos\tcode\tidLevel" << std::endl;
	for(const auto & code : codes){
		std::unordered_map<std::string, uint32_t> levels;
		std::unordered_map<uint32_t, std::vector<uint32_t>> allTakenlevels;
		for(const auto & codedRegions : code.second){
			uint32_t level = 1;
			std::set<uint32_t> takenLevels;
			for(const auto & pos : codedRegions.second){
				if(njh::in(pos.first,allTakenlevels)){
					for(const auto & t : allTakenlevels[pos.first]){
						takenLevels.emplace(t);
					}
				}
			}
			while(njh::in(level, takenLevels)){
				++level;
			}
			for(const auto & pos : codedRegions.second){
				allTakenlevels[pos.first].emplace_back(level);
			}
			levels[codedRegions.first] = level;
		}
		for(const auto & codedRegions : code.second){
			for(const auto & pos : codedRegions.second){
				codeOut << code.first
						<< "\t" << codedRegions.first
						<< "\t" << pos.first
						<< "\t" << pos.second
						<< "\t" << levels[codedRegions.first] << std::endl;
			}
		}
	}
	return 0;
}




int miscRunner::getAlnPosToRealPosTable(const njh::progutils::CmdArgs & inputCommands){

	OutOptions outOpts(bfs::path(""));
	seqSetUp setUp(inputCommands);
	setUp.processWritingOptions(outOpts);
	setUp.processReadInNames({"--fasta"}, true);
	setUp.finishSetUp(std::cout);

	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	OutputStream out(outOpts);
	out << "name\talignPosition\trealPosition" << "\n";
	seqInfo seq;
	while(reader.readNextRead(seq)){

		auto start = seq.seq_.find_first_not_of("-");
		auto end = seq.seq_.find_last_not_of("-") + 1;
		for(const auto pos : iter::range(start, end)){
			out << seq.name_
					<< "\t" << pos
					<< "\t" << getRealPosForAlnPos(seq.seq_, pos);
			out << std::endl;
		}
	}

	return 0;
}

int miscRunner::calculateShannonEntropySimple(const njh::progutils::CmdArgs & inputCommands){

	OutOptions outOpts(bfs::path(""));
	seqSetUp setUp(inputCommands);
	setUp.processWritingOptions(outOpts);
	setUp.processReadInNames({"--fasta"}, true);
	setUp.finishSetUp(std::cout);

	SeqInput reader(setUp.pars_.ioOptions_);
	auto seqs = reader.readAllReads<seqInfo>();
	Muscler::checkAlignSeqsLensThrow(seqs, __PRETTY_FUNCTION__);
	OutputStream out(outOpts);

	std::unordered_map<uint32_t, charCounter> counters(len(seqs.front()));

	for(const auto & seq : seqs){
		for(const auto pos : iter::range(len(seq))){
			counters[pos].increaseCountOfBase(seq.seq_[pos], seq.cnt_);
		}
	}
	for(auto & counter : counters){
		counter.second.resetAlphabet(false);
		counter.second.setFractions();
	}

	out << "position\tshannonEntropy\the" << "\n";
	auto counterKeys = njh::getVecOfMapKeys(counters);
	njh::sort(counterKeys);
	for(const auto  counterKey : counterKeys){
		const auto & counter = counters.at(counterKey);
		double sum = 0;
		double sumOfSquares = 0;
		for(const auto & aa : counter.alphabet_){
			sumOfSquares += std::pow(counter.fractions_[aa], 2.0);
			sum += counter.fractions_[aa] * log2(counter.fractions_[aa]);
		}
		out << counterKey
				<< "\t" << (sum == 0 ? 0 : -sum)
				<< "\t" << 1 - sumOfSquares << std::endl;
	}

	return 0;
}

} // namespace njhseq
