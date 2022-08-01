//
// Created by Nicholas Hathaway on 7/31/22.
//



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
#include <njhseq/IO/OutputStream.hpp>

namespace njhseq {

int parsingFileExpRunner::parseMashTriangleResults(const njh::progutils::CmdArgs & inputCommands) {

	IoOptions ioOpts;
	ioOpts.out_.outExtention_ = ".tab.txt";
	seqSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processVerbose();
	setUp.processWritingOptions(ioOpts.out_);
	setUp.setOption(ioOpts.in_.inFilename_, "--triangleMashFile", "output from mash triangle", true);
	setUp.finishSetUp(std::cout);


	std::string line;
	InputStream  in(ioOpts.in_);
	OutputStream out(ioOpts.out_);
	njh::files::crossPlatGetline(in, line);
	uint32_t number = njh::StrToNumConverter::stoToNum<uint32_t>(line);
	std::cout << "number: " << number << std::endl;
	uint32_t count = 0;
	VecStr names;
	std::vector<std::vector<double>> data;
	while(count < number){
		++count;
		njh::files::crossPlatGetline(in, line);
		auto lineToks = njh::tokenizeString(line, "\t");
		if(lineToks.size() != count){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "expected like toks to be " << count << ", not " << lineToks.size() << "\n";
			throw std::runtime_error{ss.str()};
		}
		names.emplace_back(lineToks[0]);
		if(count >1){
			data.emplace_back(vecStrToVecNum<double>(getSubVector(lineToks,1)));
			data.back().emplace_back(0);
		}else{
			data.emplace_back(std::vector<double>{0});
		}
	}
	out << njh::conToStr(names, "\t") << std::endl;
	for(const auto row : iter::range(0UL, data.size())){
		out << njh::conToStr(data[row], "\t");
		for(const auto col : iter::range(row + 1, data.size())){
			out << "\t" << data[col][row];
		}
		out << std::endl;
	}
	return 0;
}


} //namespace njhseq

