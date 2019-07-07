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
#include <SeekDeep/objects/ReadPairsOrganizer.hpp>


namespace njhseq {

std::vector<bfs::path> VecStrToVecPath(const VecStr & strs){
	std::vector<bfs::path> ret;
	std::transform(strs.begin(), strs.end(), std::back_inserter(ret), [](const std::string & str){return bfs::path(str);});
	return ret;
}

int genExpRunner::concatenateDifferentLanes(const njh::progutils::CmdArgs & inputCommands){
	std::string removeRegexPatStr = "";
	bool removeOldFiles = false;
	bool overWrite = false;
	uint32_t numThreads = 1;
	bfs::path dir = "./";
	bfs::path outDir = "./";
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(dir, "--dir", "directory to search");
	setUp.setOption(outDir, "--outDir", "output directory");
	setUp.setOption(removeRegexPatStr, "--removeRegexPatStr", "Remove Regex Pat Str");

	setUp.setOption(numThreads, "--numThreads", "Number of cpus to use");
	setUp.setOption(overWrite, "--overWrite", "over write output files");
	setUp.setOption(removeOldFiles, "--removeOldFiles", "remove old files");
	setUp.finishSetUp(std::cout);



	ReadPairsOrganizer pOrg{VecStr{}};

	std::regex regexPat{ReadPairsOrganizer::illuminaPat_};
	auto fullFnps = njh::files::listAllFiles("./", false, {regexPat});

	pOrg.processFiles(fullFnps);
	pOrg.readPairs_ = pOrg.readPairsUnrecognized_;
	auto pairs = pOrg.processReadPairs();

	if(setUp.pars_.debug_){
		std::cout << "fullFnps.size(): " << fullFnps.size() << std::endl;

		std::cout << "pairs.size(): " << pairs.size() << std::endl;
		for(const auto & pair : pairs){
			std::cout << pair.first << std::endl;
			if("" != removeRegexPatStr){
				std::string newName = std::regex_replace(pair.first, std::regex{removeRegexPatStr}, "");
				std::cout << newName << std::endl;
			}

			std::cout << "\tfirst" << std::endl;
			for(const auto & sub : pair.second.first){
				std::cout << "\t\t" << sub << std::endl;
				std::smatch mat;
				std::regex_match(sub, mat, ReadPairsOrganizer::illuminaPat_);
				std::cout << "\t\t" << "mat[mat.size() - 1] == \".gz\": " << njh::colorBool(mat[mat.size() - 1] == ".gz") << std::endl;
//				for(const auto & m : mat){
//					std::cout << "\t\t" << m << std::endl;
//				}
			}
			std::cout << "\tsecond" << std::endl;
			for(const auto & sub : pair.second.second){
				std::cout << "\t\t" << sub << std::endl;
				std::smatch mat;
				std::regex_match(sub, mat, ReadPairsOrganizer::illuminaPat_);
				std::cout << "\t\t" << "mat[mat.size() - 1] == \".gz\": " << njh::colorBool(mat[mat.size() - 1] == ".gz") << std::endl;
//				for(const auto & m : mat){
//					std::cout << "\t\t" << m << std::endl;
//				}
			}
		}
		return 1;
	}
	njh::files::makeDirP(njh::files::MkdirPar{outDir});

	//check if files exits already
	std::map<std::string, std::vector<bfs::path>> outputFiles;
	{
		bool fail = false;
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error the following files already exist, use --overWrite to overwrite it" << "\n";
		for (const auto & pair : pairs) {
			bool anyGz = false;
			for(const auto & sub : pair.second.first){
				std::smatch mat;
				std::regex_match(sub, mat, ReadPairsOrganizer::illuminaPat_);
				if(mat[mat.size() - 1] == ".gz"){
					anyGz = true;
					break;
				}
			}
			for(const auto & sub : pair.second.second){
				std::smatch mat;
				std::regex_match(sub, mat, ReadPairsOrganizer::illuminaPat_);
				if(mat[mat.size() - 1] == ".gz"){
					anyGz = true;
					break;
				}
			}
			std::string extension = ".fastq";
			if(anyGz){
				extension += ".gz";
			}
			std::string outStub = pair.first;
			if("" != removeRegexPatStr){
				std::string newName = std::regex_replace(pair.first, std::regex{removeRegexPatStr}, "");
				outStub = newName;
			}
			{
				OutOptions outOpts(bfs::path(njh::files::make_path(outDir, outStub + "_R1" + extension ) ) ) ;
				if (outOpts.outExists() && !overWrite) {
					ss << outOpts.outName() << "\n";
					fail = true;
				}
				outputFiles[outOpts.outName().string()] = VecStrToVecPath(pair.second.first);
			}
			{
				OutOptions outOpts(bfs::path(njh::files::make_path(outDir, outStub + "_R2" + extension ) ) ) ;
				if (outOpts.outExists() && !overWrite) {
					ss << outOpts.outName() << "\n";
					fail = true;
				}
				outputFiles[outOpts.outName().string()] = VecStrToVecPath(pair.second.second);

			}
		}
		if(fail){
			throw std::runtime_error { ss.str() };
		}
	}




	for(const auto & outs : outputFiles){
		OutOptions outOpts(bfs::path(outs.first));
		outOpts.overWriteFile_ = overWrite;
		if(setUp.pars_.verbose_){
			std::cout << "Combining " << njh::conToStr(outs.second, ", ") << " into " << outOpts.outName() << std::endl;
		}
	}

	auto keys = getVectorOfMapKeys(outputFiles);

	njh::concurrent::LockableQueue<std::string> queue(keys);

	std::function<void()> catFiles = [&queue,&outputFiles,&overWrite,&removeOldFiles](){
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
	njh::concurrent::runVoidFunctionThreaded(catFiles, numThreads);
	return 0;
}

} /* namespace njhseq */

