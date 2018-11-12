/*
 * miscRunner_calculateShannonEntropySimple.cpp
 *
 *  Created on: Sep 4, 2017
 *      Author: nicholashathaway
 */


#include "miscRunner.hpp"
#include "elucidator/simulation.h"
#include "elucidator/BamToolsUtils.h"

namespace njhseq {



int miscRunner::getSlidingQualityWindowMeans(const njh::progutils::CmdArgs & inputCommands){
	OutOptions outOpts(bfs::path(".tab.txt"));
	uint32_t windowSize = 10;
	uint32_t windowStep = 5;
	seqSetUp setUp(inputCommands);
	setUp.setOption(windowSize, "--windowSize", "Window Size");
	setUp.setOption(windowStep, "--windowStep", "Window Step");
	setUp.processReadInNames({"--fastq", "--fastqgz"}, true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	OutputStream out(outOpts);
	seqInfo seq;
	out << "name\tpos\tqualityMean" << "\n";
	while(reader.readNextRead(seq)){
		for(auto pos : iter::range<uint32_t>(0, len(seq) - windowSize + 1, windowStep)){
			out << seq.name_
					<< "\t" << pos
					<< "\t" << vectorMean(getSubVector(seq.qual_, pos, windowSize))
					<< "\n";
		}
	}

	return 0;
}



} // namespace njhseq
