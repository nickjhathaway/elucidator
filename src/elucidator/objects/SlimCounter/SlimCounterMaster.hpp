#pragma once
/*
 * SlimCounterMaster.hpp
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


#include "elucidator/common.h"
#include "elucidator/objects/SlimCounter/SlimCounterRef.hpp"

namespace njhseq {
class SlimCounterMaster {
public:
	SlimCounterMaster(bfs::path outputDir,
			const std::vector<BamTools::RefData> & refInfos);

	const bfs::path mainOutputDir_;
	const bfs::path countsOutputDir_;

	std::vector<BamTools::RefData> refInfos_;

	std::unordered_map<uint32_t, std::shared_ptr<SlimCounterRef>> counts_;

	std::shared_ptr<SlimCounterRef> getSlimCounterRef(uint32_t refId);

	std::shared_ptr<SlimCounterRef> getSlimCounterRef(
			const std::string & refName);

	void increaseCount(uint32_t refId, uint32_t pos, char base, bool reverseStand,
			bool highQuality, uint32_t count);

	std::mutex mut_;


	static void checkDirectoryStructure(const bfs::path & dir);
	static std::unordered_map<uint64_t, uint64_t> readGzChromIndex(const bfs::path & filename);

};

}  // namespace njhseq





