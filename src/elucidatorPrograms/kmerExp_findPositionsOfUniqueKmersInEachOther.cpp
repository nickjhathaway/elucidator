/*
 * kmerExp_findPositionsOfUniqueKmersInEachOther.cpp
 *
 *  Created on: Jun 17, 2020
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
#include "elucidator/objects/dataContainers.h"
#include "elucidator/simulation.h"
#include "elucidator/objects/seqObjects/seqKmers.h"

#include "elucidator/objects/MiscUtility/GenomeSeqSearch.hpp"


namespace njhseq {


int kmerExpRunner::findPositionsOfUniqueKmersInEachOther(const njh::progutils::CmdArgs & inputCommands){
	seqInfo seq1("seq1");
	seqInfo seq2("seq2");
	uint32_t kmerLength = 7;
	OutOptions outOpts(bfs::path("out.tab.txt"));
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processWritingOptions(outOpts);
	setUp.setOption(kmerLength, "--kmerLength", "kmer Length");
	setUp.processSeq(seq1, "--seq1", "First sequence", true);
	setUp.processSeq(seq2, "--seq2", "Second sequence", true);
	setUp.finishSetUp(std::cout);

	OutputStream out(outOpts);


	kmerInfo kInfo1(seq1.seq_, kmerLength, false);
	kmerInfo kInfo2(seq2.seq_, kmerLength, false);

	std::unordered_set<std::string> sharedUniqueKmers;

	for(const auto & k :kInfo1.kmers_){
		if(1 == k.second.count_ && njh::in(k.first, kInfo2.kmers_) && 1 == kInfo2.kmers_[k.first].count_){
			sharedUniqueKmers.emplace(k.first);
		}
	}

	out << "kmer\tseq1Pos\tseq2Pos" << std::endl;
	VecStr sharedUniqueKmersVec(sharedUniqueKmers.begin(), sharedUniqueKmers.end());
	njh::sort(sharedUniqueKmersVec, [&kInfo1](const std::string & k1, const std::string & k2){
		return kInfo1.kmers_[k1].positions_.front() < kInfo1.kmers_[k2].positions_.front();
	});
	for(const auto & k : sharedUniqueKmersVec){
		out << k
				<< "\t" << kInfo1.kmers_[k].positions_.front()
				<< "\t" << kInfo2.kmers_[k].positions_.front() << std::endl;
	}


	return 0;
}

}  //namespace njhseq




