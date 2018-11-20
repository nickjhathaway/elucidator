/*
 * parsingFileExp_parseSTOCKHOLM.cpp
 *
 *  Created on: Apr 20, 2018
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

#include "parsingFileExp.hpp"

namespace njhseq {


int parsingFileExpRunner::parseSTOCKHOLM(
		const njh::progutils::CmdArgs & inputCommands) {
	bfs::path fnp = "";
	OutOptions outOpts;
	seqSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processVerbose();
	setUp.processWritingOptions(outOpts);
	setUp.setOption(fnp, "--fnp", "STOCKHOLM alignment file", true);
	setUp.finishSetUp(std::cout);


	STOCKHOLMFileParser parser(fnp);
	parser.parseFile();

	parser.writeOutFile(outOpts);

	return 0;
}
int parsingFileExpRunner::parseSTOCKHOLMToFasta(
		const njh::progutils::CmdArgs & inputCommands) {
	bfs::path fnp = "";
	OutOptions outOpts;
	seqSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processVerbose();
	setUp.processWritingOptions(outOpts);
	setUp.setOption(fnp, "--fnp", "STOCKHOLM alignment file", true);
	setUp.finishSetUp(std::cout);


	STOCKHOLMFileParser parser(fnp);
	parser.parseFile();
	auto seqOut = SeqIOOptions::genFastaOut(outOpts.outName());
	seqOut.out_.transferOverwriteOpts(outOpts);
	parser.writeOutSeqFile(seqOut);

	return 0;
}
//


}  // namespace njhseq

