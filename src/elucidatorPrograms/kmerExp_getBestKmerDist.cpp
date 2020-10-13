/*
 * kmerExp_getBestKmerDist.cpp
 *
 *  Created on: Dec 3, 2017
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

int kmerExpRunner::getBestKmerDist(const njh::progutils::CmdArgs & inputCommands){
	uint32_t kmerLength = 7;
	OutOptions outOpts(bfs::path("out.tab.txt"));
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processReadInNames(true);
	setUp.processRefFilename(true);
	setUp.processWritingOptions(outOpts);
	setUp.setOption(kmerLength, "--kmerLength", "kmer Length");
	setUp.finishSetUp(std::cout);

	auto refSeqs = createKmerReadVec(setUp.pars_.refIoOptions_, kmerLength, false);

	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	seqInfo seq;
	OutputStream out(outOpts);
	out << "name\tBestRef\tdistance" << "\n";
	while(reader.readNextRead(seq)){
		kmerInfo kInfo(seq.seq_, kmerLength, false);
		double bestKmerDist = std::numeric_limits<double>::min();
		uint32_t index = std::numeric_limits<uint32_t>::max();
		for(const auto pos : iter::range(refSeqs.size())){
			const auto & refSeq = refSeqs[pos];
			auto dist = refSeq->kInfo_.compareKmers(kInfo);
			if(dist.second > bestKmerDist){
				bestKmerDist = dist.second;
				index = pos;
			}
		}
		if(std::numeric_limits<uint32_t>::max() != index){
			out << seq.name_
					<< "\t" << refSeqs[index]->seqBase_.name_
					<< "\t" << bestKmerDist
					<< "\n";
		}else{
			out << seq.name_
					<< "\t" << "*"
					<< "\t" << ""
					<< "\n";
		}
	}

	return 0;


}

} //namespace njhseq



