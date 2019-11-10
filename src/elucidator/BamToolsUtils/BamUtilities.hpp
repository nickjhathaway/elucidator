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
#include <SeekDeep/objects/PairedReadProcessor.hpp>

namespace njhseq {

struct RunCoverageFinderMultiPars{
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

void RunCoverageFinderMulti(const RunCoverageFinderMultiPars & pars);

struct RunCoverageFinderSinglePars{
	OutOptions outOpts{bfs::path(""), ".tab.txt"};
	uint32_t numThreads = 1;

	std::string bamFnp = "";
	uint32_t window = 500;
	uint32_t step = 100;
	uint32_t coverageCutOff = 1;
	uint32_t regionBatchSize = 100;
	bool byBases = false;



};

void RunCoverageFinderSingle(const RunCoverageFinderSinglePars & pars);

struct RegionRefinementPars{
	bfs::path bedFnp;
	bfs::path bamFnp;
	uint32_t numThreads = 1;
	OutOptions outOpts{"", ".bed"};

	bool reOrient = false;
};

std::vector<std::shared_ptr<Bed6RecordCore>> RunRegionRefinement(const RegionRefinementPars & pars);


struct BamCountSpecficRegionsPars{
	BamCountSpecficRegionsPars();
	uint32_t extendAmount = 30;
	uint32_t numThreads = 1;
	uint32_t mappingQuality = 20;
	uint32_t baseQuality = 25;
	double matchIDCutOff = 0.70;

	uint32_t totalCountCutOff = 5;
	uint32_t perBaseCountCutOff = 3;
	bool forcePlusStrand = false;
	bool countDuplicates = false;
	PairedReadProcessor::ProcessParams pairPars;

	void setDefaults(seqSetUp & setUp);
};

std::unordered_map<std::string, std::unordered_map<std::string, uint64_t>> BamCountSpecficRegions(
		const std::vector<GenomicRegion> & inputRegions,
		const std::unordered_map<std::string, seqInfo> & regionSeqs,
		const bfs::path bamFnp,
		const BamCountSpecficRegionsPars & pars);

}  // namespace njhseq




