/*
 * parsingFileExp_parseSTOCKHOLM.cpp
 *
 *  Created on: Apr 20, 2018
 *      Author: nick
 */


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

