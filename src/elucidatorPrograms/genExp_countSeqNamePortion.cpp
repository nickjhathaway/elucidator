/*
 * genExp_countSeqNamePortion.cpp
 *
 *  Created on: Jul 8, 2019
 *      Author: nicholashathaway
 */




#include "genExp.hpp"
#

namespace njhseq {

int genExpRunner::countSeqNamePortion(const njh::progutils::CmdArgs & inputCommands){
	std::string regexPatStr = "";
	uint32_t markCount = 1;
	OutOptions outOpts(bfs::path(""), ".tab.txt");
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(regexPatStr, "--regexPatStr", "Regex Pat Str", true);
	setUp.setOption(markCount, "--markCount", "The mark subexpression to count");
	setUp.processWritingOptions(outOpts);
	setUp.processReadInNames(true);
	setUp.finishSetUp(std::cout);

	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	OutputStream out(outOpts);
	seqInfo seq;
	std::unordered_map<std::string, uint32_t> counts;
	std::regex regexPat{regexPatStr};
	if(markCount > regexPat.mark_count()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "--markCount(" << markCount << ") can't be larger than the mark_count()(" << regexPat.mark_count() << ") of regular expression "<< "\n";
		throw std::runtime_error{ss.str()};
	}
	while(reader.readNextRead(seq)){
		std::smatch match;
		if(std::regex_match(seq.name_, match, regexPat)){
			++counts[match[markCount]];
		}
	}

	table outTab(counts, VecStr{"mark", "count"});
	outTab.sortTable("count", "mark", true);
	outTab.outPutContents(out, "\t");
	return 0;
}


}  //namespace njhseq


