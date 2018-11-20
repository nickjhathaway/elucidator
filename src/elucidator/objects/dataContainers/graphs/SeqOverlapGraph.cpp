/*
 * SeqOverlapGraph.cpp
 *
 *  Created on: Jul 18, 2016
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
#include "SeqOverlapGraph.hpp"

namespace njhseq {

SeqOverlapGraph::node::node(const std::shared_ptr<seqInfo> & val) :
		val_(val) {
}

void SeqOverlapGraph::node::resetVisitCount() {
	visitCount_ = 0;
}

void SeqOverlapGraph::node::addHead(const std::shared_ptr<edge> & e) {
	headEdges_.push_back(e);
}
void SeqOverlapGraph::node::addTail(const std::shared_ptr<edge> & e) {
	tailEdges_.push_back(e);
}

bool SeqOverlapGraph::node::headless() const {
	return headEdges_.empty();
}

bool SeqOverlapGraph::node::tailless() const {
	return tailEdges_.empty();
}

void SeqOverlapGraph::node::addToPath(std::ostream & out,
		std::string currentPath) {
	++visitCount_;
	currentPath += val_->name_;
	if (tailless()) {
		out << currentPath << std::endl;
		;
	}
	for (const auto & tail : tailEdges_) {
		tail->tail_.lock()->addToPath(out, currentPath + " -> ");
	}
}

SeqOverlapGraph::edge::edge(const std::shared_ptr<node> & head,
		const std::shared_ptr<node> & tail) :
		head_(head), tail_(tail) {

}

void SeqOverlapGraph::addNode(const std::shared_ptr<seqInfo> & n) {
	if (njh::has(nodes_, n->name_)) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": Error, alreay contains node with name "
				<< n->name_ << "\n";
		throw std::runtime_error { ss.str() };
	}
	nodes_.emplace(n->name_, std::make_shared<node>(n));
}

void SeqOverlapGraph::addEdge(const std::string & head,
		const std::string & tail) {
	/**@todo add way to check if edge already exists */
	auto headNode = nodes_.at(head);
	auto tailNode = nodes_.at(tail);
	std::shared_ptr<edge> e = std::make_shared<edge>(headNode, tailNode);
	headNode->addTail(e);
	tailNode->addHead(e);
}

void SeqOverlapGraph::writePaths(std::ostream & out) const {
	/**@todo add way to check if there are any cycles or no headless nodes */
	for (const auto & n : nodes_) {
		if (n.second->headless()) {
			n.second->addToPath(out, "");
		}
	}
	VecStr notVisitedNodes;
	for (const auto & n : nodes_) {
		if (n.second->visitCount_ == 0) {
			notVisitedNodes.emplace_back(n.first);
		}
	}
	if (!notVisitedNodes.empty()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__
				<< ": Error, the following nodes weren't visited:\n";
		ss << njh::conToStr(notVisitedNodes, ",") << "\n";
		throw std::runtime_error { ss.str() };
	}
}

Json::Value SeqOverlapGraph::createSankeyOutput() const {
	Json::Value ret;
	std::vector<std::shared_ptr<node>> nodesVec;
	std::unordered_map<std::string, uint32_t> nodePosition;
	for (const auto & n : nodes_) {
		nodePosition[n.first] = nodesVec.size();
		nodesVec.push_back(n.second);
	}
	auto &nodes = ret["nodes"];
	auto &links = ret["links"];

	for (const auto & n : nodesVec) {
		Json::Value nodeJson;
		nodeJson["name"] = n->val_->name_;
		nodeJson["cnt"] = n->val_->cnt_;
		nodeJson["frac"] = n->val_->frac_;
		nodes.append(nodeJson);
		double totalTail = 0;
		for (const auto & tl : n->tailEdges_) {
			totalTail += tl->tail_.lock()->val_->frac_;
		}
		for (const auto & tl : n->tailEdges_) {
			Json::Value linkJsons;
			linkJsons["source"] = nodePosition[tl->head_.lock()->val_->name_];
			linkJsons["target"] = nodePosition[tl->tail_.lock()->val_->name_];
			;
			linkJsons["value"] = tl->tail_.lock()->val_->frac_ / totalTail;
			links.append(linkJsons);
		}
	}
	return ret;
}

}  // namespace njhseq
