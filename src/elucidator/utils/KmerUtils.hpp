#pragma once
/*
 * KmerUtils.hpp
 *
 *  Created on: Apr 18, 2016
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

#include <njhseq/objects/kmer.h>

#include <njhseq/objects/seqObjects/seqKmers/seqWithKmerInfo.hpp>

namespace njhseq {

std::vector<std::unique_ptr<seqWithKmerInfo>> createKmerReadVec(
		const SeqIOOptions & opts, uint32_t kLength, bool setReverse);

std::vector<std::vector<double>> readDistanceMatrix(std::istream & in);

void writeDistanceMatrix(std::ostream & out,
		const std::vector<std::vector<double>> & distances);

void writeDistanceMatrix(std::ostream & out,
		const std::vector<std::vector<double>> & distances, const VecStr & names);

table getKmerStatsOnFile(const SeqIOOptions & seqFile,
		const seqWithKmerInfo & compare);

}  // namespace njhseq



