#pragma once
/*
 * randStrGen.hpp
 *
 *  Created on: Jul 27, 2014
 *      Author: nickhathaway
 */
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

#include <njhcpp/simulation/randObjGen.hpp>
#include "elucidator/simulation/simulationCommon.hpp"

namespace njhseq {

class randStrGen {
public:

	//constructor
	randStrGen(randomGenerator rGen,
			const std::vector<char> & letters):
				charGen_(njh::randObjectGen<char,uint32_t>(letters)),
				rGen_(rGen){}

	randStrGen(randomGenerator rGen,
			const std::vector<char> & letters,
			const std::vector<uint32_t> & counts):
					charGen_(njh::randObjectGen<char, uint32_t>(letters,counts)),
					rGen_(rGen){}

private:
	//members
	njh::randObjectGen<char, uint32_t> charGen_;
	randomGenerator rGen_;
public:
	//functions
	std::string rStr(uint64_t size);
	VecStr rStrs(uint64_t size, uint32_t num);
	std::string rStr(uint64_t minSize, uint64_t maxSize);
	VecStr rStrs(uint64_t minSize, uint64_t maxSize, uint32_t num);
};

} /* namespace njh */


