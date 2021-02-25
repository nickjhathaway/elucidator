/*
 * kmerExp_getAvgKmerPositionPerSeq.cpp
 *
 *  Created on: Feb 24, 2021
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

#include "kmerExp.hpp"
#include "elucidator/utils/KmerUtils.hpp"
#include "elucidator/objects/seqObjects/seqKmers.h"
#include <njhseq/IO/SeqIO/SeqIO.hpp>



namespace njhseq {


int kmerExpRunner::getAvgKmerPosition(const njh::progutils::CmdArgs & inputCommands){
	uint32_t kmerLength = 21;
	uint32_t occurenceCutOff = 5;
	OutOptions outOpts(bfs::path(""), ".tab.txt");
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processReadInNames(true);
	setUp.processWritingOptions(outOpts);
	setUp.setOption(occurenceCutOff, "--occurenceCutOff", "occurence Cut Off");
	setUp.setOption(kmerLength, "--kmerLength", "kmer Length");


	setUp.finishSetUp(std::cout);

	OutputStream out(outOpts);
	out << "kmer\tcount\tavgPosition\tstd" << "\n";
	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	seqInfo seq;
	std::unordered_map<std::string, std::vector<uint32_t>> kmerPositons;
	while(reader.readNextReadLock(seq)){
		if(len(seq) > kmerLength){
			kmerInfo kInfo(seq.seq_, kmerLength, false);
			for(const auto & k : kInfo.kmers_){
				addOtherVec(kmerPositons[k.first], k.second.positions_);
			}
		}
	}

	for(const auto &kpos : kmerPositons){
		if(kpos.second.size() < occurenceCutOff){
			continue;
		}
		out << kpos.first
				<< "\t" << kpos.second.size()
				<< "\t" << vectorMean(kpos.second)
				<< "\t" << vectorStandardDeviationPop(kpos.second) << std::endl;
	}




	return 0;


}

int kmerExpRunner::getAvgKmerPositionPerSeq(const njh::progutils::CmdArgs & inputCommands){
	uint32_t kmerLength = 21;
	uint32_t occurenceCutOff = 5;
	OutOptions outOpts(bfs::path(""), ".tab.txt");
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processReadInNames(true);
	setUp.processWritingOptions(outOpts);
	setUp.setOption(occurenceCutOff, "--occurenceCutOff", "occurence Cut Off");
	setUp.setOption(kmerLength, "--kmerLength", "kmer Length");


	setUp.finishSetUp(std::cout);

	OutputStream out(outOpts);
	out << "seq\tseqPos\tkmer\tcount\tavgPosition\tstd" << "\n";
	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	seqInfo seq;
	std::unordered_map<std::string, std::vector<uint32_t>> kmerPositons;
	while(reader.readNextReadLock(seq)){
		if(len(seq) > kmerLength){
			kmerInfo kInfo(seq.seq_, kmerLength, false);
			for(const auto & k : kInfo.kmers_){
				addOtherVec(kmerPositons[k.first], k.second.positions_);
			}
		}
	}

	std::unordered_map<std::string, std::pair<double, double>> kmerPositionsAvgPosStds;

	reader.reOpenIn();



	for(const auto &kpos : kmerPositons){
		if(kpos.second.size() < occurenceCutOff){
			continue;
		}
		kmerPositionsAvgPosStds[kpos.first] = std::make_pair(vectorMean(kpos.second), vectorStandardDeviationPop(kpos.second));
	}

	while(reader.readNextReadLock(seq)){
		if(len(seq) > kmerLength){
			for(const auto pos : iter::range(len(seq) -kmerLength + 1)){
				auto k = seq.seq_.substr(pos, kmerLength);
				out << seq.name_
						<< "\t" << pos
						<< "\t" << k
						<< "\t" << kmerPositons[k].size()
						<< "\t" << kmerPositionsAvgPosStds[k].first
						<< "\t" << kmerPositionsAvgPosStds[k].second << std::endl;
			}
		}
	}



	return 0;


}



int kmerExpRunner::filterPerSeqAvgKmerPosition(const njh::progutils::CmdArgs & inputCommands){
	uint32_t kmerLength = 21;
	uint32_t occurenceCutOff = 5;
	double stdCutOff = 5;
	uint32_t within = 5;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processDefaultReader(true);

	setUp.setOption(occurenceCutOff, "--occurenceCutOff", "Occurrence Cut Off");
	setUp.setOption(kmerLength, "--kmerLength", "k-mer Length");
	setUp.setOption(stdCutOff, "--stdCutOff", "standard deviation Cut Off");
	setUp.setOption(within, "--within", "within this many bases from the edges");

	setUp.finishSetUp(std::cout);

	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	reader.openOut();

	SeqIOOptions filterOpts = setUp.pars_.ioOptions_;
	filterOpts.out_.outFilename_ = njh::files::prependFileBasename(filterOpts.out_.outFilename_, "filtered_");
	SeqOutput filterWriter(filterOpts);
	filterWriter.openOut();

	SeqIOOptions tempOpts = setUp.pars_.ioOptions_;
	tempOpts.out_.outFilename_ = njh::files::prependFileBasename(filterOpts.out_.outFilename_, "temp_");
	SeqOutput tempWriter(tempOpts);
	tempWriter.openOut();

	seqInfo seq;
	std::unordered_map<std::string, std::vector<uint32_t>> kmerPositons;
	std::unordered_map<std::string, std::vector<uint32_t>> kmerPositonsSecondPass;

	std::unordered_set<std::string> kmersBelowStdCutOff;

	while(reader.in_.readNextReadLock(seq)){
		if(len(seq) > kmerLength){
			kmerInfo kInfo(seq.seq_, kmerLength, false);
			for(const auto & k : kInfo.kmers_){
				addOtherVec(kmerPositons[k.first], k.second.positions_);
			}
		}
	}

	for(const auto &kpos : kmerPositons){
		if(kpos.second.size() < occurenceCutOff){
			continue;
		}
		if(vectorStandardDeviationPop(kpos.second) < stdCutOff){
			kmersBelowStdCutOff.emplace(kpos.first);
		}
	}

	reader.in_.reOpenIn();
	while(reader.in_.readNextReadLock(seq)){
		if(len(seq) > kmerLength){
			bool remove = false;
			for(const auto pos : iter::range(std::min<uint32_t>(within, len(seq) - kmerLength))){
				auto k = seq.seq_.substr(pos, kmerLength);
				if(njh::in(k, kmersBelowStdCutOff)){
					remove = true;
					break;
				}
			}
			if(!remove){
				for(const auto pos : iter::range(std::min<uint32_t>(within, len(seq) - kmerLength))){
					auto k = seq.seq_.substr(len(seq)-kmerLength - pos, kmerLength);
					if(njh::in(k, kmersBelowStdCutOff)){
						remove = true;
						break;
					}
				}
			}
			if(remove){
				filterWriter.write(seq);
			}else{
				kmerInfo kInfo(seq.seq_, kmerLength, false);
				for (const auto &k : kInfo.kmers_) {
					addOtherVec(kmerPositonsSecondPass[k.first], k.second.positions_);
				}
				tempWriter.write(seq);
			}
		}
	}
	for(const auto &kpos : kmerPositonsSecondPass){
		if(kpos.second.size() < occurenceCutOff){
			continue;
		}
		if(vectorStandardDeviationPop(kpos.second) < stdCutOff){
			kmersBelowStdCutOff.emplace(kpos.first);
		}
	}
	tempWriter.closeOut();
	SeqIOOptions secondPassRederOpts(tempOpts.out_.outName(), SeqIOOptions::getInFormat(tempOpts.outFormat_), false);
	SeqInput secondPassReder(secondPassRederOpts);
	secondPassReder.openIn();

	while(secondPassReder.readNextReadLock(seq)){
		if(len(seq) > kmerLength){
			bool remove = false;
			for(const auto pos : iter::range(std::min<uint32_t>(within, len(seq) - kmerLength))){
				auto k = seq.seq_.substr(pos, kmerLength);
				if(njh::in(k, kmersBelowStdCutOff)){
					remove = true;
					break;
				}
			}
			if(!remove){
				for(const auto pos : iter::range(std::min<uint32_t>(within, len(seq) - kmerLength))){
					auto k = seq.seq_.substr(len(seq)-kmerLength - pos, kmerLength);
					if(njh::in(k, kmersBelowStdCutOff)){
						remove = true;
						break;
					}
				}
			}
			if(remove){
				filterWriter.write(seq);
			}else{
				reader.out_.write(seq);
			}
		}
	}
	secondPassReder.closeIn();
	bfs::remove(tempOpts.out_.outName());

	return 0;
}

} //namespace njhseq

