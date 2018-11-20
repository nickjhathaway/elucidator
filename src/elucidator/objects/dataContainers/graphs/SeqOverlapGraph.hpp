#pragma once
/*
 * SeqOverlapGraph.hpp
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

#include "elucidator/common.h"

namespace njhseq {
class SeqOverlapGraph {
public:

	class edge;
	class node {
	public:

		node(const std::shared_ptr<seqInfo> & val);
		std::shared_ptr<seqInfo> val_;

		std::vector<std::shared_ptr<edge>> headEdges_;
		std::vector<std::shared_ptr<edge>> tailEdges_;

		uint32_t visitCount_ = 0;

		void resetVisitCount();

		void addHead(const std::shared_ptr<edge> & e);
		void addTail(const std::shared_ptr<edge> & e);

		bool headless() const;

		bool tailless() const;

		void addToPath(std::ostream & out, std::string currentPath);
	};
	class edge {
	public:
		edge(const std::shared_ptr<node> & head,
				const std::shared_ptr<node> & tail);
		std::weak_ptr<node> head_;
		std::weak_ptr<node> tail_;

	};
	std::unordered_map<std::string, std::shared_ptr<node>> nodes_;

	void addNode(const std::shared_ptr<seqInfo> & n);

	void addEdge(const std::string & head, const std::string & tail);

	void writePaths(std::ostream & out) const;

	Json::Value createSankeyOutput() const;


};


}  // namespace njhseq




