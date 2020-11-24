/*
 * benchMarkings.cpp
 *
 *  Created on: Feb 2, 2017
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

#include "benchMarking.hpp"
#include "elucidator/objects/BioDataObject.h"
#include "elucidator/BamToolsUtils.h"



namespace njhseq {


benchMarkingRunner::benchMarkingRunner()
    : njh::progutils::ProgramRunner(
          {
	 	 	 	 	 addFunc("benchGzWritingOneRead", benchGzWritingOneRead, false),
					 addFunc("benchGzWritingSetChunks", benchGzWritingSetChunks, false),
					 addFunc("benchGzWritingRdBuf", benchGzWritingRdBuf, false),
					 addFunc("mapLookUp", mapLookUp, false),
           },//
          "benchMarking") {}

int benchMarkingRunner::benchGzWritingOneRead(const njh::progutils::CmdArgs & inputCommands){
	bfs::path inputFnp = "";
	OutOptions outOpts;
	outOpts.outExtention_ = ".gz";
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processWritingOptions(outOpts);
	setUp.setOption(inputFnp, "--input", "Input file", true);
	setUp.finishSetUp(std::cout);

	njh::files::checkExistenceThrow(inputFnp, __PRETTY_FUNCTION__);
	if("" == outOpts.outFilename_ ){
		outOpts.outFilename_ = inputFnp.string() + ".gz";
		outOpts.outExtention_ = ".gz";
	}
	{
		njh::GZSTREAM::ogzstream outstream;
		outOpts.openGzFile(outstream);
		outstream << njh::files::get_file_contents(inputFnp, false);
	}
	return 0;
}

int benchMarkingRunner::benchGzWritingSetChunks(const njh::progutils::CmdArgs & inputCommands){
	bfs::path inputFnp = "";
	uint32_t chunckSize = 4096 * 10;
	OutOptions outOpts;
	outOpts.outExtention_ = ".gz";
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processWritingOptions(outOpts);
	setUp.setOption(inputFnp, "--input", "Input file", true);
	setUp.setOption(chunckSize, "--chunckSize", "chunck size");
	setUp.finishSetUp(std::cout);

	njh::files::checkExistenceThrow(inputFnp, __PRETTY_FUNCTION__);
	if("" == outOpts.outFilename_ ){
		outOpts.outFilename_ = inputFnp.string() + ".gz";
		outOpts.outExtention_ = ".gz";
	}
	{
		njh::GZSTREAM::ogzstream outstream;
		outOpts.openGzFile(outstream);
		std::ifstream infile(inputFnp.string(), std::ios::binary);
		std::vector<char> buffer(chunckSize);
		infile.read(buffer.data(), sizeof(char) * chunckSize);
		std::streamsize bytes = infile.gcount();

		while(bytes > 0){
			outstream.write(buffer.data(), bytes * sizeof(char));
			infile.read(buffer.data(), sizeof(char) * chunckSize);
			bytes = infile.gcount();
		}
	}
	return 0;
}

int benchMarkingRunner::benchGzWritingRdBuf(const njh::progutils::CmdArgs & inputCommands){
	bfs::path inputFnp = "";
	OutOptions outOpts;
	outOpts.outExtention_ = ".gz";
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processWritingOptions(outOpts);
	setUp.setOption(inputFnp, "--input", "Input file", true);
	setUp.finishSetUp(std::cout);

	njh::files::checkExistenceThrow(inputFnp, __PRETTY_FUNCTION__);
	if("" == outOpts.outFilename_ ){
		outOpts.outFilename_ = inputFnp.string() + ".gz";
		outOpts.outExtention_ = ".gz";
	}
	{
		njh::GZSTREAM::ogzstream outstream;
		outOpts.openGzFile(outstream);
		std::ifstream infile(inputFnp.string(), std::ios::binary);
		outstream << infile.rdbuf();
	}

	return 0;
}

} // namespace njhseq
