/*
 * ContigsCompareGraphDev.cpp
 *
 *  Created on: Nov 5, 2019
 *      Author: nicholashathaway
 */


#include "ContigsCompareGraph.hpp"

namespace njhseq {

ContigsCompareGraphDev::edge::edge(const std::shared_ptr<node> & head,
		const std::shared_ptr<node> & tail, uint32_t cnt, const std::vector<ConnectorInfo> & connectorInfos) :
		head_(head), tail_(tail), cnt_(cnt), connectorInfos_(connectorInfos) {
}



std::string ContigsCompareGraphDev::edge::createUid() const {
	return head_.lock()->uid_ + "_" + tail_.lock()->uid_;
}



}  // namespace njhseq
