/*
 * parsingFileExp_parsehmmerDomainHitTab.cpp
 *
 *  Created on: Apr 20, 2018
 *      Author: nick
 */

#include "parsingFileExp.hpp"

namespace njhseq {


int parsingFileExpRunner::parsehmmerDomainHitTab(
		const njh::progutils::CmdArgs & inputCommands) {
	bfs::path fnp = "";
	OutOptions outOpts;
	seqSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processVerbose();
	setUp.processWritingOptions(outOpts);
	setUp.setOption(fnp, "--fnp", "hmmer domain hits file", true);
	setUp.finishSetUp(std::cout);

	OutputStream out(outOpts);
	BioDataFileIO<HmmerDomainHitTab> reader{IoOptions(InOptions(fnp))};
	HmmerDomainHitTab domain;
	reader.openIn();
	out << "[";
	while(reader.readNextRecord(domain)){
		out << domain.toJson() << std::endl;
	}
	out << "]" << std::endl;

	return 0;
}

}  // namespace njhseq


