#pragma once
/*
 * KmerVecUtils.hpp
 *
 *  Created on: May 24, 2016
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


namespace njhseq {

std::vector<std::vector<double>> getKmerAccerDistance(
		std::vector<std::unique_ptr<seqWithKmerInfo>>& reads, uint32_t kmerStart,
		uint32_t kmerStop, uint32_t numThreads, bool useKNumber, bool verbose);

readDistGraph<double> genKmerAccerDistGraphWithDbSmartBuild(
		std::vector<std::unique_ptr<seqWithKmerInfo>> & reads, uint32_t kmerStop,
		const readDistGraph<double>::dbscanPars & pars, bool breakLargeIndelCons,
		uint32_t largeIndel, aligner & alignerObj, bool doTies, uint32_t numThreads,
		bool verbose);

readDistGraph<double> genKmerAccerDistGraph(
		std::vector<std::unique_ptr<seqWithKmerInfo>> & reads, uint32_t kmerStart,
		uint32_t kmerStop, double avgKmerSimDecreasingRateCutOff,
		bool breakLargeIndelCons, uint32_t largeIndel, aligner & alignerObj,
		bool doTies, uint32_t numThreads, bool verbose);

readDistGraph<double> genKmerAccerDistGraphWithDb(
		std::vector<std::unique_ptr<seqWithKmerInfo>> & reads, uint32_t kmerStart,
		uint32_t kmerStop, const readDistGraph<double>::dbscanPars & pars,
		bool breakLargeIndelCons, uint32_t largeIndel, aligner & alignerObj,
		bool doTies, uint32_t numThreads, bool verbose);

}  // namespace njhseq
