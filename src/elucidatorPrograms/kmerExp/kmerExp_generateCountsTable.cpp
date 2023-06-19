/*
 * kmerExp_generateCountsTable.cpp
 *
 *  Created on: Nov 6, 2017
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
#include <njhseq/objects/kmer/KmerUtils.hpp>
#include "elucidator/objects/dataContainers.h"
#include "elucidator/simulation.h"
#include "elucidator/objects/seqObjects/seqKmers.h"

#include "elucidator/objects/MiscUtility/GenomeSeqSearch.hpp"


namespace njhseq {
int kmerExpRunner::generateCountsTable(const njh::progutils::CmdArgs & inputCommands){
	bfs::path inputFnp = "";
	uint32_t kmerLength = 7;
	uint32_t occurenceCutOff = 5;
	OutOptions outOpts(bfs::path("outmat.tab.txt"));
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processReadInNames(true);
	setUp.processWritingOptions(outOpts);
	setUp.setOption(occurenceCutOff, "--occurenceCutOff", "Occurence Cut Off");
	setUp.setOption(kmerLength, "--kmerLength", "kmer Length");
	setUp.finishSetUp(std::cout);

	OutputStream out(outOpts);


	SeqInput seqReader(setUp.pars_.ioOptions_);
	seqReader.openIn();
	seqInfo seq;
	std::map<std::string, uint32_t> kmers;
	while(seqReader.readNextRead(seq)){
		auto kInfo = seqWithKmerInfo(seq, kmerLength, false);
		for(const auto & k : kInfo.kInfo_.kmers_){
			++kmers[k.first];
		}
	}
	seqReader.reOpenIn();
	//collapseToUnique
	out << "sample" ;
	for(const auto & k : kmers){
		if(k.second >= occurenceCutOff){
			out << "\t" << k.first;
		}
	}
	out << "\n";

	while(seqReader.readNextRead(seq)){
		auto kInfo = seqWithKmerInfo(seq, kmerLength, false);
		out << seq.name_;
		for(const auto & k : kmers){
			if(k.second >= occurenceCutOff){
				if(njh::in(k.first, kInfo.kInfo_.kmers_)){
					out << "\t" << 1;
				}else{
					out << "\t" << 0;
				}
			}
		}
		out << "\n";
	}


	return 0;
}

}  //namespace njhseq

