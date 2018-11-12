#pragma once
/*
 * SnpGatherer.hpp
 *
 *  Created on: Jun 20, 2016
 *      Author: nick
 */

#include <TwoBit.h>
#include "elucidator/common.h"

namespace njhseq {

class BaseGatherer {
public:

	BaseGatherer(const std::vector<BamTools::RefData> & refInfos);

	const std::vector<BamTools::RefData> refInfos_;


};


}  // namespace njhseq




