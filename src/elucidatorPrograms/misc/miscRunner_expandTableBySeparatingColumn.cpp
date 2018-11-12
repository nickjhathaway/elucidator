/*
 * miscRunner_expandTableBySeparatingColumn.cpp
 *
 *  Created on: Oct 29, 2017
 *      Author: nick
 */




#include "miscRunner.hpp"


namespace njhseq {



int miscRunner::expandTableBySeparatingColumn(const njh::progutils::CmdArgs & inputCommands){
	auto tabOpts = TableIOOpts::genTabFileOut("", true);
	tabOpts.inDelim_ = "\t";
	std::string columnName = "";
	seqSetUp setUp(inputCommands);
	setUp.setOption(tabOpts.in_.inFilename_, "--fnp", "File path for the table", true);
	setUp.setOption(tabOpts.hasHeader_, "--header", "Header is Present");
	setUp.setOption(tabOpts.inDelim_, "--delim", "File Delimiter");
	setUp.setOption(columnName, "--columnName", "Column Name to separate", true);
	if (tabOpts.inDelim_ == "tab") {
		tabOpts.inDelim_ = "\t";
	}
	tabOpts.outDelim_ = tabOpts.inDelim_;

	setUp.processWritingOptions(tabOpts.out_);
	setUp.finishSetUp(std::cout);

	table namesTab(tabOpts);
	namesTab.checkForColumnsThrow(VecStr{columnName}, __PRETTY_FUNCTION__)	;

	auto colPos = namesTab.getColPos(columnName);
	OutputStream out(tabOpts.out_);
	if(tabOpts.hasHeader_){
		out << njh::conToStr(namesTab.columnNames_, tabOpts.outDelim_) << "\n";
	}
	for(const auto & row : namesTab.content_){
		auto toks =  tokenizeString(row[colPos], ",");
		auto rowCopy = row;
		for(const auto & tok : toks){
			rowCopy[colPos] = tok;
			out << njh::conToStr(rowCopy, tabOpts.outDelim_) << "\n";
		}
	}
	return 0;
}

} // namespace njhseq

