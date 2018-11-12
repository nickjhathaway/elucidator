#pragma once
/*
 * MinorVariantsMixture.hpp
 *
 *  Created on: Jan 9, 2016
 *      Author: nick
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
#include <njhcpp/simulation/randomGenerator.hpp>
#include "elucidator/common.h"


namespace njhseq {
class MinorVariantsMixture{
public:

	MinorVariantsMixture(const seqInfo & major, const std::string & mutationsString, size_t mutationOffSet);
	MinorVariantsMixture(const seqInfo & major);
	MinorVariantsMixture(size_t stringLength);
	void setNewMajor(size_t stringLength);
	std::shared_ptr<seqInfo> major_;
private:
	std::unordered_map<std::string, std::unordered_map<size_t, char>> variantsMutations_;
	std::unordered_map<size_t, char> alreadyMutatedPosBase_;
public:
	njh::randomGenerator gen_;
	std::vector<std::shared_ptr<seqInfo>> variants_;

	std::vector<char> dnaAlphabet_{'A', 'C', 'G', 'T'};
	std::vector<uint32_t> initialNumMutations_ {1,2,3,4,6,8,13};
	std::unordered_map<std::string, std::unordered_map<std::string, double>> mixturesAbundances_;

	uint32_t numOfVars() const;

	void addVariants(const std::string & mutationsString, size_t mutationOffSet);
	void addVariants(const std::vector<std::set<size_t>> & mutPositions);

	void createDisparatePairsMixture(uint32_t numberOfMinorStrains);
	void createGradientDisparatePairsMixture(uint32_t numberOfMinorStrains);
	void createVaryingMinorVariantsMixture(uint32_t numberOfMinorVariants,
			double oneVariantPerIdDiff = 0.04);

	void addAnEqualAbundanceMixture();
	void createAbundancesForMixtures(const std::string & str);
	void outputAbundanceFile(const std::string & abundanceFilename,
			bool dualReplicates = false) const;
	void outputVariantMutationInfofile(const std::string & filename) const;
	void writeFasta(const std::string & dirname) const;
	void writeFasta(const std::string & dirname, const std::string & filename) const;

	void prependPaddingSeq(const std::string & seq);
	void appendPaddingSeq(const std::string & seq);
};


}  // namespace njhseq




