#pragma once
/*
 * BamUtilities.hpp
 *
 *  Created on: Jun 10, 2018
 *      Author: nick
 */


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




