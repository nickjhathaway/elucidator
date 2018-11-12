/*
 * kmerExp_getBestKmerDist.cpp
 *
 *  Created on: Dec 3, 2017
 *      Author: nick
 */



#include "kmerExp.hpp"
#include "elucidator/utils/KmerUtils.hpp"
#include "elucidator/objects/seqObjects/seqKmers.h"

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



