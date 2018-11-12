/*
 * SnpGatherer.cpp
 *
 *  Created on: Jun 20, 2016
 *      Author: nick
 */




#include "BaseGatherer.hpp"



namespace njhseq {

BaseGatherer::BaseGatherer(const std::vector<BamTools::RefData> & refInfos) :
		refInfos_(refInfos) {

}


}  // namespace njhseq

