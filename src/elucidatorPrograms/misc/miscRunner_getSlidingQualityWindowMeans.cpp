/*
 * miscRunner_calculateShannonEntropySimple.cpp
 *
 *  Created on: Sep 4, 2017
 *      Author: nicholashathaway
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
