/*
 * genExp_concatLane.cpp
 *
 *  Created on: Jul 7, 2017
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

#include "genExp.hpp"
#include "elucidator/simulation.h"
#include "elucidator/utils/KmerUtils.hpp"
#include "elucidator/objects/BioDataObject.h"
#include <njhseq/GenomeUtils.h>

#include <TwoBit.h>



namespace njhseq {

int genExpRunner::concatenateDifferentLanes(const njh::progutils::CmdArgs & inputCommands){
	std::string regexPatStr = "(.*)(_L[0-9]{3}_)(.*)";
	bool removeOldFiles = false;
	bool overWrite = false;
	uint32_t numThreads = 1;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.setOption(numThreads, "--numThreads", "Number of cpus to use");
	setUp.setOption(overWrite, "--overWrite", "over write output files");
	setUp.setOption(removeOldFiles, "--removeOldFiles", "remove old files");
	setUp.finishSetUp(std::cout);
	std::regex regexPat{regexPatStr};

	auto fullFnps = njh::files::listAllFiles("./", false, {regexPat});

	std::vector<bfs::path> fnps;
	for(const auto & fnp : fullFnps){
		fnps.emplace_back(fnp.first.filename());
	}
	std::unordered_map<std::string, std::vector<bfs::path>> outputFiles;
	for(const auto & fnp : fnps){
		std::smatch match;
		if(std::regex_match(fnp.string(), match, regexPat)){
			if(4 == match.size()){
				std::stringstream ss;
				ss << match[1] << "_" << match[3];
				outputFiles[ss.str()].emplace_back(fnp);
			}
		}
	}

	//check if files exits already
	for (const auto & outs : outputFiles) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error the following files already exist, use --overWrite to overwrite it" << "\n";
		bool fail = false;
		OutOptions outOpts(bfs::path(outs.first));
		outOpts.overWriteFile_ = overWrite;
		if (outOpts.outExists() && !overWrite) {
			ss << outOpts.outName() << "\n";
			fail = true;
		}
		if(fail){
			throw std::runtime_error { ss.str() };
		}
	}
	for(const auto & outs : outputFiles){
		OutOptions outOpts(bfs::path(outs.first));
		outOpts.overWriteFile_ = overWrite;
		if(setUp.pars_.verbose_){
			std::cout << "Combining " << njh::conToStr(outs.second, ", ") << " into " << outOpts.outName();
		}
	}
	auto keys = getVectorOfMapKeys(outputFiles);

	njh::concurrent::LockableQueue<std::string> queue(keys);

	auto catFiles = [&queue,&outputFiles,&overWrite,&removeOldFiles](){
		std::string outName = "";
		while(queue.getVal(outName)){
			OutOptions outOpts{bfs::path(outName)};
			outOpts.overWriteFile_ = overWrite;
			concatenateFiles(outputFiles.at(outName), outOpts);
			if(removeOldFiles){
				for(const auto & fnp : outputFiles.at(outName)){
					bfs::remove(fnp);
				}
			}
		}
	};

	std::vector<std::thread> threads;
	for(uint32_t t = 0; t < numThreads; ++t){
		threads.emplace_back(catFiles);
	}
	for(auto & t : threads){
		t.join();
	}
	return 0;
}

} /* namespace njhseq */

