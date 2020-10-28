/*
 * kmerExp_kmerConnectionGraph.cpp
 *
 *  Created on: Oct 26, 2020
 *      Author: nick
 */




#include "kmerExp.hpp"
#include "elucidator/utils/KmerUtils.hpp"
#include "elucidator/objects/seqObjects/seqKmers.h"
#include <njhseq/IO/SeqIO/SeqIO.hpp>



namespace njhseq {




class KmerSimpleConGraph {




};


int kmerExpRunner::kmerConnectionGraph(const njh::progutils::CmdArgs & inputCommands){
	uint32_t kmerLength = 4;
	OutOptions outOpts(bfs::path("out.tab.txt"));
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processReadInNames(true);
	setUp.processWritingOptions(outOpts);
	setUp.setOption(kmerLength, "--kmerLength", "kmer Length");
	setUp.finishSetUp(std::cout);




	return 0;


}

int kmerExpRunner::countKmers(const njh::progutils::CmdArgs & inputCommands){



	//uint32_t numThreads = 1;
	uint32_t kmerLength = 4;
	OutOptions outOpts(bfs::path(""), ".tab.txt");
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processReadInNames(true);
	setUp.processWritingOptions(outOpts);
	setUp.setOption(kmerLength, "--kmerLength", "kmer Length");
	//setUp.setOption(numThreads, "--numThreads", "Number of threads to use");
	setUp.finishSetUp(std::cout);


	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();

	OutputStream out(outOpts);
	std::unordered_map<std::string, uint64_t> counts;

	auto addCounts = [&counts,&kmerLength](const seqInfo & seq){
		if(len(seq) >= kmerLength){
			for(const auto pos : iter::range(len(seq) + 1 - kmerLength)){
				++counts[seq.seq_.substr(pos, kmerLength)];
			}
		}
	};

	if(setUp.pars_.ioOptions_.isPairedIn()){
		PairedRead seq;
		while(reader.readNextRead(seq)){
			addCounts(seq.seqBase_);
			addCounts(seq.mateSeqBase_);
		}
	} else {
		seqInfo seq;
		while(reader.readNextRead(seq)){
			addCounts(seq);
		}
	}

	VecStr keys = njh::getVecOfMapKeys(counts);
	njh::sort(keys);
	for(const auto & key : keys){
		out << key << "\t" << counts[key] << "\n";
	}

	return 0;


}




} //namespace njhseq


