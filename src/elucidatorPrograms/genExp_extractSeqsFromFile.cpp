/*
 * extractSeqsFromFile.cpp
 *
 *  Created on: Aug 6, 2017
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


