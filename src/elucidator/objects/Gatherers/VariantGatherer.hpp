#pragma once
/*
 * SnpGatherer.hpp
 *
 *  Created on: Jun 20, 2016
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
#include "elucidator/objects/Gatherers/BaseGatherer.hpp"
#include "elucidator/objects/SlimCounter/SlimCounterPos.hpp"

namespace njhseq {

class VariantGatherer : public BaseGatherer{
public:

	struct GathVarPars {
		bool addIndels_ = false;
		bool addNonRefSnps_ = false;
		bool onlyBiallelicVariants_ = false;
	};

	struct VarPosRes {
		bool variant_;
		uint32_t numberOfAlleles_;
	};

	VariantGatherer(const std::vector<BamTools::RefData> & refInfos,
			double fracCutOff,
			uint32_t depthCutOff,
			double strandBiasCutOff);

	const double fracCutOff_;
	const uint32_t depthCutOff_;
	const double strandBiasCutOff_;



	table gatherVariantsFromChromFile(const std::string & filename,
			const std::string & sampleName,
			const GathVarPars pars,
			TwoBit::TwoBitFile & twoBitfile);

	VarPosRes determineIfPosVariant(const SlimCounterPos & pos,
			const char refBase,
			const GathVarPars pars);

};


}  // namespace njhseq




