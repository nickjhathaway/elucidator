/*
 * seqUtilsInfoRunner_getSlidingEntropy.cpp
 *
 *  Created on: Jun 27, 2020
 *      Author: nick
 */





#include "seqUtilsInfoRunner.hpp"
#include <njhseq/IO/SeqIO/SeqIO.hpp>
#include <njhseq/objects/seqObjects/seqKmers.h>


namespace njhseq {


int seqUtilsInfoRunner::getSlidingEntropy(const njh::progutils::CmdArgs & inputCommands) {

	OutOptions outOpts(bfs::path(""), ".tab.txt");
	uint32_t windowStep = 5;
	uint32_t windowSize = 50;
	uint32_t kLen = 1;
	seqSetUp setUp(inputCommands);
	setUp.description_ = "Trim front and back of sequences with a sliding window for low entropy, seqs will be removed if all is low entropy";

	setUp.processVerbose();
	setUp.processDebug();
	setUp.processWritingOptions(outOpts);
	setUp.setOption(windowStep, "--windowStep", "Window Step");
	setUp.setOption(windowSize, "--windowSize", "Window Size");
	setUp.setOption(kLen, "--kLen", "Kmer length for entropy calculation");

	setUp.processReadInNames(true);
	setUp.finishSetUp(std::cout);

	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	seqInfo seq;
	OutputStream out(outOpts);
	out << "name\tposition\tsubseq\tentropy" << std::endl;
	while(reader.readNextRead(seq)){
		if(len(seq) >= windowSize){
			for (auto pos : iter::range<uint32_t>(0, len(seq) - windowSize + 1, windowStep)) {
				kmerInfo kInfo(seq.seq_.substr(pos, windowSize), kLen, false);
				out << seq.name_
						<< "\t" << pos
						<< "\t" << seq.seq_.substr(pos, windowSize)
						<< "\t" << kInfo.computeKmerEntropy() << std::endl;
			}
		}
	}

	if(setUp.pars_.verbose_){
		setUp.logRunTime(std::cout);
	}
	return 0;
}



}  // namespace njhseq

