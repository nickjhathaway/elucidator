/*
 * randomFileCreator.cpp
 *
 *  Created on: Mar 23, 2014
 *      Author: nickhathaway
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
#include "randomFileCreator.hpp"
#include "elucidator/simulation/randomStrGen.hpp"

namespace njhseq {

// functions


// functions
void randomFileCreator::randomFile(uint32_t lenStart, uint32_t lenStop,
		uint32_t numOfSeqs, bool processed, uint32_t bottomAmount,
		uint32_t topAmount, bool fastq, std::ostream& out) {
	if (fastq) {
		for (const auto pos : iter::range(numOfSeqs)) {
			auto len = rgen_.unifRand(lenStart, lenStop);
			seqInfo info("Seq." + njh::leftPadNumStr(pos, numOfSeqs),
					simulation::randStr(len, counter_, rgen_),
					rgen_.unifRandVector(qualStart_, qualStop_, len));
			if (processed) {
				info.cnt_ = rgen_.unifRand(bottomAmount, topAmount);
				info.updateName();
			}
			info.outPutFastq(out);
		}
	} else {
		for (const auto pos : iter::range(numOfSeqs)) {
			auto len = rgen_.unifRand(lenStart, lenStop);
			seqInfo info("Seq." + njh::leftPadNumStr(pos, numOfSeqs),
					simulation::randStr(len, counter_, rgen_));
			if (processed) {
				info.cnt_ = rgen_.unifRand(bottomAmount, topAmount);
				info.updateName();
			}
			info.outPutSeq(out);
		}
	}
}



} /* namespace njh */
