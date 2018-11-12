/*
 * extractSeqsFromFile.cpp
 *
 *  Created on: Aug 6, 2017
 *      Author: nick
 */


#include "genExp.hpp"


namespace njhseq {

int genExpRunner::getSeqFromFile(const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	setUp.processDefaultReader(VecStr { "-fastq", "-fasta", "--fastagz", "-fastqgz", "--fastq1", "--fastq1gz"}, true);
	uint32_t pos = 0;
	uint32_t number = 1;
	setUp.setOption(pos, "--pos", "Read Position to Read", true);
	setUp.setOption(number, "--number", "Number of reads to read", true);
	setUp.processDebug();
	setUp.finishSetUp(std::cout);
	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	reader.openOut();
	reader.in_.loadIndex();
	reader.in_.seekToSeqIndex(pos);
	if (setUp.pars_.ioOptions_.isPairedIn()) {
		uint32_t count = 0;
		PairedRead seq;
		while (reader.in_.readNextRead(seq) && count < number) {
			reader.write(seq);
			++count;
		}
	} else {
		uint32_t count = 0;
		seqInfo seq;
		while (reader.in_.readNextRead(seq) && count < number) {
			reader.write(seq);
			++count;
		}
	}
	return 0;
}





} // namespace njhseq


