/*
 * parsingFileExp_parsehmmerDomainHitTab.cpp
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
#include <njhseq/objects/BioDataObject/BioDataFileIO.hpp>

#include <njhseq/objects/BioDataObject/HmmerDomainHitTab.hpp>

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


