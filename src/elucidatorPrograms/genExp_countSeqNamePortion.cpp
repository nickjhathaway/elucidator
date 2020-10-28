/*
 * genExp_countSeqNamePortion.cpp
 *
 *  Created on: Jul 8, 2019
 *      Author: nicholashathaway
 */




#include "genExp.hpp"
#include <SeekDeep/objects.h>

namespace njhseq {


int genExpRunner::countSeqNamePortion(const njh::progutils::CmdArgs & inputCommands){
	std::string regexPatStr = "([A-Za-z0-9_]+):([0-9]+):([A-Za-z0-9-]+):([0-9]+):([0-9]+):([0-9]+):([0-9]+) ([12]):([NY]):([0-9]):([AGCTN+]+)";
	uint32_t markCount = 11;
	std::string name = "";
	OutOptions outOpts(bfs::path(""), ".tab.txt");
	uint64_t numberOfReadsToCount = std::numeric_limits<uint64_t>::max();
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(numberOfReadsToCount, "--numberOfReadsToCount", "Number Of Reads To Count", false);
	setUp.setOption(regexPatStr, "--regexPatStr", "Regex Pat Str", false);
	setUp.setOption(markCount, "--markCount", "The mark subexpression to count");
	setUp.setOption(name, "--name", "name to add to output file");

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
	uint64_t totalCount = 0;


	while(reader.readNextRead(seq)){
		std::smatch match;
		if(std::regex_match(seq.name_, match, regexPat)){
			++counts[match[markCount]];
		}
		++totalCount;
		if(totalCount >= numberOfReadsToCount){
			break;
		}
	}

	table outTab(counts, VecStr{"mark", "count"});
	if("" != name){
		outTab.addColumn({name}, "name");
	}
	outTab.sortTable("count", "mark", true);
	outTab.outPutContents(out, "\t");
	return 0;
}

int genExpRunner::countIlluminaSampleNumber(const njh::progutils::CmdArgs & inputCommands){
	OutOptions outOpts(bfs::path(""), ".tab.txt");
	uint32_t testNumber = std::numeric_limits<uint32_t>::max();
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processWritingOptions(outOpts);
	setUp.processReadInNames(true);
	setUp.setOption(testNumber, "--testNumber", "Number of reads to read");
	setUp.finishSetUp(std::cout);

	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	OutputStream out(outOpts);
	seqInfo seq;
	std::unordered_map<std::string, uint32_t> counts;
	uint32_t totalCount = 0;
	while(reader.readNextRead(seq)){
		IlluminaNameFormatDecoder decoded(seq.name_);
		++counts[decoded.getSampleNumber()];
		++totalCount;
		if(totalCount > testNumber){
			break;
		}
	}

	table outTab(counts, VecStr{"sample", "count"});
	outTab.sortTable("count", "sample", true);
	outTab.outPutContents(out, "\t");
	return 0;
}




}  //namespace njhseq


