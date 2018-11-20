#pragma once
/*
 * BamUtilities.hpp
 *
 *  Created on: Jun 10, 2018
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

#include <njhseq/BamToolsUtils.h>

namespace njhseq {

struct CoverageFinderPars{
	OutOptions outOpts{bfs::path(""), ".tab.txt"};
	uint32_t numThreads = 1;
	std::string bams = "";
	std::string pat = ".*.bam$";
	uint32_t window = 500;
	uint32_t step = 100;
	uint32_t coverageCutOff = 1;
	uint32_t regionBatchSize = 10;
	bool byBases = false;



};
void RunCoverageFinder(const CoverageFinderPars & pars);

struct RegionRefinementPars{
	bfs::path bedFnp;
	bfs::path bamFnp;
	uint32_t numThreads = 1;
	OutOptions outOpts{"", ".bed"};

	bool reOrient = false;
};

void RunRegionRefinement(const RegionRefinementPars & pars);


}  // namespace njhseq




