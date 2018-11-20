/*
 * SlimCounterMaster.cpp
 *
 *  Created on: Jun 16, 2016
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

#include "SlimCounterMaster.hpp"
#include <njhseq/BamToolsUtils/BamToolsUtils.hpp>
#include "elucidator/filesystem/GzSimpleBinFile.hpp"
#include "elucidator/BamToolsUtils/BamCountingUtils.hpp"


namespace njhseq {

SlimCounterMaster::SlimCounterMaster(bfs::path outputDir,
		const std::vector<BamTools::RefData> & refInfos) :
		mainOutputDir_(outputDir), countsOutputDir_(
				njh::files::make_path(outputDir, "counts")), refInfos_(refInfos) {
	if (bfs::exists(mainOutputDir_)) {
		std::stringstream ss;
		ss << "Error in " << __PRETTY_FUNCTION__ << ", directory " << mainOutputDir_
				<< " already exists" << std::endl;
		throw std::runtime_error { ss.str() };
	}
	njh::files::makeDirP(njh::files::MkdirPar(countsOutputDir_.string()));
	for (auto pos : iter::range<uint32_t>(refInfos_.size())) {
		counts_.emplace(pos,
				std::make_shared<SlimCounterRef>(pos,
						njh::files::make_path(countsOutputDir_, refInfos_[pos].RefName)));
	}
	auto refTab = refDataVecToTab(refInfos_);
	TableIOOpts outOpts = TableIOOpts::genTabFileOut(njh::files::make_path(mainOutputDir_, "refIdLookup.tab.txt").string(), true);
	refTab.outPutContents(outOpts);
}

std::shared_ptr<SlimCounterRef> SlimCounterMaster::getSlimCounterRef(uint32_t refId){
	std::lock_guard<std::mutex> lock(mut_);
	auto search = counts_.find(refId);
	if(search == counts_.end()){
		/**@todo should we just throw an exception instead? */
		return nullptr;
	}else{
		return search->second;
	}
}

std::shared_ptr<SlimCounterRef> SlimCounterMaster::getSlimCounterRef(const std::string & refName){
	uint32_t refId = std::numeric_limits<uint32_t>::max();
	for(const auto pos : iter::range(refInfos_.size())){
		if(refInfos_[pos].RefName == refName){
			refId = pos;
			break;
		}
	}
	return getSlimCounterRef(refId);
}

void SlimCounterMaster::increaseCount(uint32_t refId, uint32_t pos, char base,
		bool reverseStand, bool highQuality, uint32_t count) {
	counts_.at(refId)->increaseCountLockFree(pos, base, reverseStand, highQuality, count);
}


void SlimCounterMaster::checkDirectoryStructure(const bfs::path & dir){
	std::stringstream ss;
	bool failed = false;
	//check for main dir
	if(!bfs::exists(dir)){
		failed = true;
		ss << njh::bashCT::boldRed(dir.string()) << " doesn't exist" << std::endl;
	}
	//check for stats dir
	if(!bfs::exists(njh::files::make_path(dir, "stats"))){
		failed = true;
		ss << njh::bashCT::boldRed(njh::files::make_path(dir, "stats").string()) << " doesn't exist" << std::endl;
	}
	//check for data dir
	if(!bfs::exists(njh::files::make_path(dir, "data"))){
		failed = true;
		ss << njh::bashCT::boldRed(njh::files::make_path(dir, "data").string()) << " doesn't exist" << std::endl;
	}
	//check for counts dir
	if(!bfs::exists(njh::files::make_path(dir, "data", "counts"))){
		failed = true;
		ss << njh::bashCT::boldRed(njh::files::make_path(dir, "data", "counts").string()) << " doesn't exist" << std::endl;
	}

	//check for refIdLookupFile
	if(!bfs::exists(njh::files::make_path(dir, "data", "refIdLookup.tab.txt"))){
		failed = true;
		ss << njh::bashCT::boldRed(njh::files::make_path(dir, "data", "refIdLookup.tab.txt").string()) << " doesn't exist" << std::endl;
	}
	if(failed){
		throw std::runtime_error{ss.str()};
	}
}

std::unordered_map<uint64_t, uint64_t> SlimCounterMaster::readGzChromIndex(const bfs::path & filename){
	if(!bfs::exists(filename)){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": Error file " << filename << " doesn't exist" << std::endl;
		throw std::runtime_error{ss.str()};
	}
	std::unordered_map<uint64_t, uint64_t> ret;
	GzSimpleBinFile<uint64_t, SlimCounterRef::IndexNumOfElements> inFile(filename);
	auto d = inFile.genDataContainer();
	while(inFile.read(d)){
		ret[d.data_[0]] = d.data_[2];
	}
	return ret;
}



}  // namespace njhseq


