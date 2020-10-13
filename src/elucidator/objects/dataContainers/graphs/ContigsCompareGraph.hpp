#pragma once

/*
 * ContigsCompareGraphDev.hpp
 *
 *  Created on: Nov 5, 2019
 *      Author: nicholashathaway
 */

#include "elucidator/common.h"
#include <njhseq/alignment/alignerUtils/comparison.hpp>
#include <njhseq/objects/seqObjects/BaseObjects/seqInfo.hpp>


namespace njhseq {
class ContigsCompareGraphDev {
public:

	ContigsCompareGraphDev(uint32_t klen);
	ContigsCompareGraphDev(uint32_t klen, uint32_t occurenceCutOff);
	uint32_t klen_;
	uint32_t occurenceCutOff_{0};
	bool debug_ = false;
	bool verbose_ = false;

	bool hasNode(const std::string & nodeName) const;
	void addNode(const std::string & k, uint32_t cnt, uint32_t kLen);
	void populateNodesFromCounts();
	void setOccurenceCutOff(uint32_t cutOff);
	void resetNodeUids();
	void resetNodePositions();
	void removeOffNodes();


	void removeOffEdges();
	void removeOffNodesOffEdgesAndReset();
	void resetNodeVisitCounts();

	class edge;
	class node {
	public:

		node(const std::string & k, uint32_t cnt, uint32_t klen);
		std::string k_;
		std::string uid_;
		uint32_t cnt_;
		uint32_t kLen_;

		uint32_t group_ = std::numeric_limits<uint32_t>::max();

		std::vector<std::shared_ptr<edge>> headEdges_;
		std::vector<std::shared_ptr<edge>> tailEdges_;

		uint32_t visitCount_ = 0;
		bool on_ = true;

		std::set<std::string> inReadNamesIdx_;

		void resetVisitCount();

		void addHead(const std::shared_ptr<edge> & e);

		void addTail(const std::shared_ptr<edge> & e);

		bool headless() const;
		uint32_t headCount() const;
		std::shared_ptr<edge> getFirstOnHeadEdge() const;
		std::shared_ptr<edge> getLastOnHeadEdge() const;

		bool tailless() const;
		uint32_t tailCount() const;
		std::shared_ptr<edge> getFirstOnTailEdge() const;
		std::shared_ptr<edge> getLastOnTailEdge() const;

		std::shared_ptr<node> getHeadNode(const std::string & k) const;

		std::shared_ptr<node> getTailNode(const std::string & k) const;
	};

	class edge {
	public:


		struct ConnectorInfo {
			ConnectorInfo(const std::string & name, uint32_t head, uint32_t tail) :
					readName_(name), headPos_(head), tailPos_(tail) {

			}
			std::string readName_;
			uint32_t headPos_;
			uint32_t tailPos_;

			bool sameInfo(const ConnectorInfo & other)const{
				return readName_ == other.readName_&& headPos_ == other.headPos_ && tailPos_ == other.tailPos_;
			}
		};

		edge(const std::shared_ptr<node> & head,
				const std::shared_ptr<node> & tail,
				uint32_t cnt,
				const std::vector<ConnectorInfo> & connectorInfos);

		std::weak_ptr<node> head_;
		std::weak_ptr<node> tail_;
		uint32_t cnt_;
		std::vector<ConnectorInfo> connectorInfos_;
		bool on_ = true;

		std::string createUid() const;
	};



	std::vector<std::shared_ptr<node>> nodes_;
	std::unordered_map<std::string, uint32_t> nodePositions_;
	std::unordered_map<std::string, uint32_t> kCounts_;


	void increaseKCounts(const std::string & seq);
	void threadThroughSequence(
			const seqInfo & seq, const std::string & threadingSeqName);
	void writeRectangleDot(std::ostream & out, bool noLabels = false) const;
	void writeRectangleDotColorBySampleCount(std::ostream & out, bool noLabels = false) const;
	std::unordered_map<std::string, std::string> writeRectangleDotColorBySamples(std::ostream & out, bool noLabels = false) const;


	bool splitNodesWithRedundantKmers();

	void collapseSingleLinkedPathsSameReads();
	void collapseSingleLinkedPaths();


	bool collapseLowFreqNodes(const comparison & allowableError, uint32_t lowFreqCutOff);

	void removeNullNodes();


	std::vector<seqInfo> nodesToSeqs() const;



	void resetGroups() const ;


private:
};

} //namespace njhseq





