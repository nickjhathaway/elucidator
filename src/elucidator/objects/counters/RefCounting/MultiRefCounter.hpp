#pragma once
//
// njhseq - A library for analyzing sequence data
// Copyright (C) 2012-2016 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
// Jeffrey Bailey <Jeffrey.Bailey@umassmed.edu>
//
// This file is part of njhseq.
//
// njhseq is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// njhseq is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with njhseq.  If not, see <http://www.gnu.org/licenses/>.
//
/*
 * MultiRefCounter.hpp
 *
 *  Created on: Jan 7, 2016
 *      Author: nick
 */



#include "elucidator/objects/counters/RefCounting/RefCounter.hpp"

namespace njhseq {

class MultiRefCounter {
public:
	MultiRefCounter(const std::string & twoBitFilename);
	RefCounter masterCounter_;
	std::unordered_map<std::string, RefCounter> allCounters_;
	std::mutex mut_;

	void addCounter(const std::string & refName, const RefCounter & counter,
			bool replace = false);
	bool hasCounter(const std::string & refName);
	void setMasterCounter();
};

}  // namespace njhseq


