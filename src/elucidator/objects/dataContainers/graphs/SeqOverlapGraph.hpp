#pragma once
/*
 * SeqOverlapGraph.hpp
 *
 *  Created on: Jul 18, 2016
 *      Author: nick
 */

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




