/*
 * miscRunner_expandTableBySeparatingColumn.cpp
 *
 *  Created on: Oct 29, 2017
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
#include "miscRunner.hpp"

#include <njhseq/IO/SeqIO/SeqIO.hpp>


namespace njhseq {
int miscRunner::expandTableBySeparatingColumn(const njh::progutils::CmdArgs & inputCommands){
	auto tabOpts = TableIOOpts::genTabFileOut("", false);
	tabOpts.inDelim_ = "\t";
	std::string columnName = "";
	std::string sep =  ",";
	std::string append = "";
	std::string prepend = "";
	seqSetUp setUp(inputCommands);
	setUp.setOption(tabOpts.in_.inFilename_, "--fnp", "File path for the table", true);
	setUp.setOption(tabOpts.hasHeader_, "--header", "Header is Present");
	setUp.setOption(tabOpts.inDelim_, "--delim", "File Delimiter");
	setUp.setOption(sep, "--sep", "separator to expand column on");
	setUp.setOption(append, "--append", "append this onto the column elements");
	setUp.setOption(prepend, "--prepend", "prepend this onto the column elements");
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
		auto toks =  tokenizeString(row[colPos], sep);
		auto rowCopy = row;
		uint32_t count = 0;
		for(const auto & tok : toks){
			std::string element = tok;
			if(0 != count && "" != prepend){
				element = prepend + element;
			}
			if(toks.size() != (count + 1) && "" != append){
				element.append(append);
			}
			++count;
			rowCopy[colPos] = element;
			out << njh::conToStr(rowCopy, tabOpts.outDelim_) << "\n";
		}
	}
	return 0;
}
} // namespace njhseq
