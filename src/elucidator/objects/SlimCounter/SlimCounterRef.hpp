#pragma once
/*
 * SlimCounterRef.hpp
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
#include "elucidator/filesystem/wrappers/BGZFCPP.hpp"
#include "elucidator/objects/SlimCounter/SlimCounterPos.hpp"

namespace njhseq {
class SlimCounterRef {
public:
	SlimCounterRef(const uint32_t & refId, bfs::path outputFile);

	const uint32_t refId_;
	const bfs::path outputFile_;

	std::unique_ptr<BGZFCPP> bgzOutFile_{nullptr};

	const static uint32_t IndexNumOfElements = 3;

	std::unordered_map<uint32_t, SlimCounterPos> counts_;

	std::shared_timed_mutex mut_;
	//std::unique_lock<std::shared_timed_mutex> lock(mut_);
	//std::shared_lock<std::shared_timed_mutex> lock(mut_);

	void initiateResetPositions(uint32_t start, uint32_t end);

	void increaseCountLockFree(uint32_t pos, char base, bool reverseStand, bool highQuality,
				uint32_t count);

	void increaseInsertionCountLockFree(uint32_t pos, bool reverseStand, uint32_t count);

	void increaseDeletionCountLockFree(uint32_t pos, bool reverseStand, uint32_t count);

	uint32_t maxPos();

	uint32_t minPos();
	uint32_t minPosLockFree() const;

	uint32_t numOfCounts();
	uint32_t numOfCountsLockFree() const;

	bool flushCountsToFile(uint32_t fromPositionIncluding, uint32_t upToPosNotIncluding);

	void indexFile();

};



}  // namespace njhseq


