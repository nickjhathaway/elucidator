#pragma once
/*
 * SnpGatherer.hpp
 *
 *  Created on: Jun 20, 2016
 *      Author: nick
 */

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




