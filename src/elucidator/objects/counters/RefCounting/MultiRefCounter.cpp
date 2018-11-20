/*
 * MultiRefCounter.cpp
 *
 *  Created on: Jan 7, 2016
 *      Author: nick
 */
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

#include "MultiRefCounter.hpp"

namespace njhseq {

MultiRefCounter::MultiRefCounter(const std::string & twoBitFilename) :
		masterCounter_(twoBitFilename) {
}

void MultiRefCounter::addCounter(const std::string & refName,
		const RefCounter & counter, bool replace) {
	std::lock_guard<std::mutex> lock(mut_);
	auto search = allCounters_.find(refName);
	if (search == allCounters_.end()) {
		allCounters_.emplace(refName, counter);
	} else if (replace) {
		allCounters_.erase(search);
		allCounters_.emplace(refName, counter);
	} else {
		std::stringstream ss;
		ss << "Error in " << __PRETTY_FUNCTION__ << std::endl;
		ss << "already contains " << refName << " use replace=true to replace" << std::endl;
	}
}

bool MultiRefCounter::hasCounter(const std::string & refName) {
	std::lock_guard<std::mutex> lock(mut_);
	return allCounters_.find(refName) != allCounters_.end();
}

void MultiRefCounter::setMasterCounter() {
	std::lock_guard<std::mutex> lock(mut_);
	masterCounter_.resetCounts();
	for (const auto & counter : allCounters_) {
		masterCounter_.merge(counter.second);
	}
}

}  // namespace njhseq

