/*
 * miscRunner_expandOutMaln.cpp
 *
 *  Created on: Sep 4, 2017
 *      Author: nicholashathaway
 */



#include "miscRunner.hpp"
#include "elucidator/simulation.h"
#include "elucidator/BamToolsUtils.h"

namespace njhseq {




int miscRunner::expandOutCollapsedToUnique(const njh::progutils::CmdArgs & inputCommands){
	bfs::path inputNames = "";
	seqSetUp setUp(inputCommands);
	setUp.pars_.ioOptions_.out_.outFilename_ = "";
	setUp.processDefaultReader({ "--fasta" }, true);
	setUp.setOption(inputNames, "--inputNames", "Input names", true);
	setUp.finishSetUp(std::cout);

	table namesTab(inputNames, "\t", true);
	namesTab.checkForColumnsThrow(VecStr{"name","reads"}, __PRETTY_FUNCTION__)	;

	std::unordered_map<std::string, VecStr> readNames;
	auto nameColPos = namesTab.getColPos("name");
	auto readsColPos = namesTab.getColPos("reads");
	for(const auto & row : namesTab.content_){
		readNames[row[nameColPos]] = tokenizeString(row[readsColPos], ",");
	}

	SeqIO reader(setUp.pars_.ioOptions_);

	reader.openIn();
	reader.openOut();

	seqInfo seq;
	while(reader.readNextRead(seq)){
		if(!njh::in(seq.name_, readNames)){
			reader.write(seq);
		}else{
			for(const auto & outName : readNames.at(seq.name_)){
				auto outSeq = seq;
				outSeq.name_ = outName;
				reader.write(outSeq);
			}
		}

	}
	return 0;
}

} // namespace njhseq
