/*
 * ContigsCompareGraphDev.cpp
 *
 *  Created on: Nov 5, 2019
 *      Author: nicholashathaway
 */


#include "ContigsCompareGraph.hpp"
#include <njhseq/alignment.h>

#include <njhseq/readVectorManipulation/readVectorOperations/massGetters.hpp>
#include <njhseq/IO/SeqIO.h>
#include <njhseq/objects/BioDataObject/BioRecordsUtils/BedUtility.hpp>

namespace njhseq {




ContigsCompareGraphDev::ContigsCompareGraphDev(uint32_t klen) :
		klen_(klen), occurenceCutOff_(0) {
}
ContigsCompareGraphDev::ContigsCompareGraphDev(uint32_t klen,uint32_t occurenceCutOff) :
		klen_(klen), occurenceCutOff_(occurenceCutOff) {
}
void ContigsCompareGraphDev::setOccurenceCutOff(uint32_t cutOff){
	occurenceCutOff_ = cutOff;
}

void ContigsCompareGraphDev::increaseKCounts(const std::string & seq) {
	if(seq.size() > klen_){
		for (auto pos : iter::range(seq.size() - klen_ + 1)) {
			kCounts_[seq.substr(pos, klen_)] += 1;
		}
	}
}



bool ContigsCompareGraphDev::hasNode(const std::string & nodeName) const {
	for(const auto & node : nodes_){
		if(nodeName == node->k_){
			return true;
		}
	}
	return false;
}
void ContigsCompareGraphDev::populateNodesFromCounts() {
	nodes_.clear();
	nodePositions_.clear();
	for (const auto & kCount : kCounts_) {
		if (kCount.second > occurenceCutOff_) {
			//addNode(kCount.first, kCount.second, kCount.first.size());
			addNode(kCount.first, kCount.second, klen_);
		}
	}
	resetNodePositions();
}


void ContigsCompareGraphDev::addNode(const std::string & k, uint32_t cnt, uint32_t kLen){
	nodes_.emplace_back(std::make_shared<node>(k, cnt, kLen));
}




void ContigsCompareGraphDev::resetNodeUids(){
	std::unordered_map<std::string, uint32_t> previouslyFound;
	for( auto & node : nodes_){
		if(njh::in(node->k_, previouslyFound)){
			++previouslyFound[node->k_];
			node->uid_ = node->k_ + estd::to_string(previouslyFound[node->k_]);
		}else{
			previouslyFound[node->k_] = 0;
			node->uid_ = node->k_;
		}
	}
	nodePositions_.clear();
}



void ContigsCompareGraphDev::resetNodeVisitCounts() {
	for (auto & n : nodes_) {
		if (n->on_) {
			n->visitCount_ = 0;
		}
	}
}

void ContigsCompareGraphDev::removeOffNodesOffEdgesAndReset(){
	std::vector<uint32_t> toRemove;
	for(const auto nodePos : iter::range(nodes_.size())){
		const auto & n = nodes_[nodePos];
		if(!n->on_){
			toRemove.emplace_back(nodePos);
			//turn of any edges that might be pointing to this node
			for(auto & head : n->headEdges_){
				head->on_ = false;
			}
			for(auto & tail : n->tailEdges_){
				tail->on_ = false;
			}
		}
	}

	if(!toRemove.empty()){
		std::sort(toRemove.rbegin(), toRemove.rend());
		for(const auto & remove : toRemove){
			nodes_.erase(nodes_.begin() + remove);
		}
	}
	//always need to remove off nodes if going to name the function this way
	removeOffEdges();
	//always need to reset if going to name the function this way
	resetNodePositions();
}

void ContigsCompareGraphDev::resetNodePositions(){
	resetNodeUids();
	nodePositions_.clear();
	for(const auto pos : iter::range(nodes_.size())){
		nodePositions_[nodes_[pos]->uid_] = pos;
	}
}

void ContigsCompareGraphDev::removeOffNodes(){
	std::vector<uint32_t> toRemove;
	for(const auto nodePos : iter::range(nodes_.size())){
		const auto & n = nodes_[nodePos];
		if(!n->on_){
			toRemove.emplace_back(nodePos);
			//turn of any edges that might be pointing to this node
			for(auto & head : n->headEdges_){
				head->on_ = false;
			}
			for(auto & tail : n->tailEdges_){
				tail->on_ = false;
			}
		}
	}
	if(!toRemove.empty()){
		std::sort(toRemove.rbegin(), toRemove.rend());
		for(const auto & remove : toRemove){
			nodes_.erase(nodes_.begin() + remove);
		}
		removeOffEdges();
		resetNodePositions();
	}
}

void ContigsCompareGraphDev::removeOffEdges(){
	//remove off edges
	for(const auto &  n : nodes_){
		{
			//remove tails
			std::vector<uint32_t> toErase;
			for( const auto tailPos : iter::range(n->tailEdges_.size())){
				if(!n->tailEdges_[tailPos]->on_){
					toErase.emplace_back(tailPos);
				}
			}
			std::sort(toErase.rbegin(), toErase.rend());
			for(const auto pos : toErase){
				n->tailEdges_.erase(n->tailEdges_.begin() + pos);
			}
		}
		{
			//remove heads
			std::vector<uint32_t> toErase;
			for( const auto headPos : iter::range(n->headEdges_.size())){
				if(!n->headEdges_[headPos]->on_){
					toErase.emplace_back(headPos);
				}
			}
			std::sort(toErase.rbegin(), toErase.rend());
			for(const auto pos : toErase){
				n->headEdges_.erase(n->headEdges_.begin() + pos);
			}
		}
	}
}


bool ContigsCompareGraphDev::splitNodesWithRedundantKmers(){
	//std::cout << __FILE__ << " " << __LINE__ << std::endl;
	resetNodeVisitCounts();
	std::vector<uint32_t> nodesToProcess;
	for(const auto nodePos : iter::range(nodes_.size())){
		const auto & node = *nodes_[nodePos];
		if(node.headCount() > 1 || node.tailCount() > 1){
			std::unordered_map<std::string, std::unordered_map<uint32_t, uint32_t>> readNameCountsHeads;
			std::unordered_map<std::string, std::unordered_map<uint32_t, uint32_t>> readNameCountsTails;
			for(const auto headPos : iter::range(node.headEdges_.size())){
				for(const auto & con : node.headEdges_[headPos]->connectorInfos_){
					++readNameCountsHeads[con.readName_][headPos];
				}
			}
			for(const auto tailPos : iter::range(node.tailEdges_.size())){
				for(const auto & con : node.tailEdges_[tailPos]->connectorInfos_){
					++readNameCountsTails[con.readName_][tailPos];
				}
			}
			VecStr readWithMultiHeadsTails;
//			if("GATC" == node.uid_){
//				std::cout << __FILE__ << " " << __LINE__ << std::endl;
//				std::cout << "GATC" << std::endl;
//				std::cout << "readNameCountsHeads: " << readNameCountsHeads.size() << std::endl;
//			}
			for(const auto & readHead : readNameCountsHeads){
				if(readHead.second.size() > 1 || readNameCountsTails[readHead.first].size() > 1){
					readWithMultiHeadsTails.emplace_back(readHead.first);
				}
			}
			if(!readWithMultiHeadsTails.empty()){
				nodesToProcess.emplace_back(nodePos);
			}
		}
	}

	if(!nodesToProcess.empty()){
		for(const auto & nodePos : nodesToProcess){
			//std::cout << __FILE__ << " " << __LINE__ << std::endl;
			std::cout << nodePos << std::endl;
			std::cout << "\t" << nodes_[nodePos]->k_ << std::endl;
			//std::cout << __FILE__ << " " << __LINE__ << std::endl;
			std::cout << njh::bashCT::red << "heads(" << nodes_[nodePos]->headEdges_.size()<< "): " << std::endl;
			for(const auto & head : nodes_[nodePos]->headEdges_){
				std::cout << "\t" <<head->head_.lock()->k_ << std::endl;
				for(const auto & con : head->connectorInfos_){
					std::cout << "\t\t" << "con.readName_: " << con.readName_ << std::endl;
					std::cout << "\t\t" << "con.headPos_ : " << con.headPos_ << std::endl;
					std::cout << "\t\t" << "con.tailPos_ : " << con.tailPos_ << std::endl;
				}
			}
			//std::cout << __FILE__ << " " << __LINE__ << std::endl;
			std::cout << njh::bashCT::cyan  << "tails(" << nodes_[nodePos]->tailEdges_.size()<< "): " << std::endl;
			for(const auto & tail : nodes_[nodePos]->tailEdges_){
				//std::cout << __FILE__ << " " << __LINE__ << std::endl;
				std::cout << "tail->on_: " << njh::boolToStr(tail->on_) << std::endl;
				std::cout << "tail->connectorInfos_.size(): " << tail->connectorInfos_.size() << std::endl;
				//std::cout << __FILE__ << " " << __LINE__ << std::endl;
				std::cout << "\t" <<tail->tail_.lock()->k_ << std::endl;
				//std::cout << __FILE__ << " " << __LINE__ << std::endl;
				for(const auto & con : tail->connectorInfos_){
					std::cout << "\t\t" << "con.readName_: " << con.readName_ << std::endl;
					std::cout << "\t\t" << "con.headPos_ : " << con.headPos_ << std::endl;
					std::cout << "\t\t" << "con.tailPos_ : " << con.tailPos_ << std::endl;
				}
				//std::cout << __FILE__ << " " << __LINE__ << std::endl;
			}
			//std::cout << __FILE__ << " " << __LINE__ << std::endl;
			std::cout << njh::bashCT::reset;
		}
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
	}
	//std::cout << __FILE__ << " " << __LINE__ << std::endl;
	struct ConnectorInfoWithHeadTailInfo{
		ConnectorInfoWithHeadTailInfo(
				const ContigsCompareGraphDev::edge::ConnectorInfo & info,
				const std::string & head,
				const std::string & tail) :
				info_(info), head_(head), tail_(tail) {

		}
		ContigsCompareGraphDev::edge::ConnectorInfo info_;
		std::string head_;
		std::string tail_;
	};

	struct ConnectorInfoPath {
		ConnectorInfoPath(const ConnectorInfoWithHeadTailInfo & head,
				const ConnectorInfoWithHeadTailInfo & tail) :
				head_(head), tail_(tail) {

		}
		ConnectorInfoWithHeadTailInfo head_;
		ConnectorInfoWithHeadTailInfo tail_;

	};
	if(!nodesToProcess.empty()){
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		std::cout << njh::bashCT::red;
		for(const auto & n : nodes_){
			if(0 != n->visitCount_){
				std::cout << n->uid_ << ": visitCount: " << n->visitCount_ << std::endl;
			}
		}
		std::cout << njh::bashCT::reset;
		bool newNodessAdded = false;
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		for(const auto & nodePos : nodesToProcess){
			if(nodes_[nodePos]->visitCount_ > 0){
				continue;
			}
			std::unordered_map<std::string, std::vector<ConnectorInfoWithHeadTailInfo>> tailConnectorInfos;
			std::unordered_map<std::string, std::vector<ConnectorInfoWithHeadTailInfo>> headConnectorInfos;
			std::cout << __FILE__ << " " << __LINE__ << std::endl;
			std::cout << "nodes_[nodePos]->uid_: " << nodes_[nodePos]->uid_ << std::endl;
			std::cout << '\t' << "head cons: " << std::endl;
			for (const auto & head : nodes_[nodePos]->headEdges_) {
				for (const auto & con : head->connectorInfos_) {
					std::cout << "\t\tcon.headPos_ : " << con.headPos_ << std::endl;
					std::cout << "\t\tcon.tailPos_ : " << con.tailPos_ << std::endl;
					std::cout << "\t\tcon.readName_: " << con.readName_ << std::endl;
					headConnectorInfos[con.readName_].emplace_back(ConnectorInfoWithHeadTailInfo(con, head->head_.lock()->uid_, nodes_[nodePos]->uid_));
				}
			}
			std::cout << '\t' << "tail cons: " << std::endl;
			for (const auto & tail : nodes_[nodePos]->tailEdges_) {
				for (const auto & con : tail->connectorInfos_) {
					std::cout << "\t\tcon.headPos_ : " << con.headPos_ << std::endl;
					std::cout << "\t\tcon.tailPos_ : " << con.tailPos_ << std::endl;
					std::cout << "\t\tcon.readName_: " << con.readName_ << std::endl;

					tailConnectorInfos[con.readName_].emplace_back(ConnectorInfoWithHeadTailInfo(con, nodes_[nodePos]->uid_, tail->tail_.lock()->uid_));
				}
			}
			VecStr nonEvenConnectors;

			for( auto & head : headConnectorInfos){
				if(head.second.size() != tailConnectorInfos[head.first].size()){
					std::cout << __FILE__ << " " << __LINE__ << std::endl;
					std::cout << "\t" << "head.first: "<< head.first << std::endl;

					std::cout << "\t\t" << "head.second.size(): " << head.second.size() << std::endl;
					std::cout << "\t\t" << "tailConnectorInfos[head.first].size(): " << tailConnectorInfos[head.first].size() << std::endl;
					nonEvenConnectors.emplace_back(head.first);
				}
				njh::sort(head.second, [](const ConnectorInfoWithHeadTailInfo & p1, const ConnectorInfoWithHeadTailInfo & p2){
					return p1.info_.headPos_ < p2.info_.headPos_;
				});
			}
			for( auto & tail : tailConnectorInfos){

				if(tail.second.size() != headConnectorInfos[tail.first].size()){
					std::cout << __FILE__ << " " << __LINE__ << std::endl;
					std::cout << "\t" << "tail.first: " << tail.first << std::endl;
					std::cout << "\t\t" << "tail.second.size(): " << tail.second.size() << std::endl;
					std::cout << "\t\t" << "headConnectorInfos[head.first].size(): " << headConnectorInfos[tail.first].size() << std::endl;
					nonEvenConnectors.emplace_back(tail.first);
				}
				njh::sort(tail.second, [](const ConnectorInfoWithHeadTailInfo & p1, const ConnectorInfoWithHeadTailInfo & p2){
					return p1.info_.headPos_ < p2.info_.headPos_;
				});
			}
			//std::cout << __FILE__ << " " << __LINE__ << std::endl;
//			if(!nonEvenConnectors.empty()){
//				std::stringstream ss;
//				ss << __PRETTY_FUNCTION__ << ", error " << "uneven amount of connectors for node: " << nodes_[nodePos]->uid_ << " for reads:\n"
//						<< njh::conToStr(nonEvenConnectors, ",")<< "\n";
//				throw std::runtime_error{ss.str()};
//			}
			//std::cout << __FILE__ << " " << __LINE__ << std::endl;
			std::unordered_map<std::string, std::unordered_map<std::string, std::unordered_map<std::string, std::vector<ConnectorInfoPath>>>> paths;
//			if("AACAGAAGAAACAGAAGAAACAG" == nodes_[nodePos]->k_){
//				std::cout << njh::bashCT::cyan;
//			}
//			if("TAAA" == nodes_[nodePos]->uid_){
//				for(const auto & head : headConnectorInfos){
//					for(const auto & hCon : head.second){
//						for(const auto & tCon : tailConnectorInfos[head.first]){
//							std::cout << "positions difference              : " << tCon.info_.headPos_ - hCon.info_.headPos_ << std::endl;
//							std::cout << "nodes_[nodePos]->k_.size() - klen_: " << nodes_[nodePos]->k_.size() - klen_ + 1 << std::endl;
//							std::cout << "\t" << "hCon.info_.headPos_" << hCon.info_.headPos_ << std::endl;
//							std::cout << "\t" << "hCon.info_.headPos_" << hCon.info_.tailPos_ << std::endl;
//							std::cout << "\t" << "tCon.info_.headPos_" << tCon.info_.headPos_ << std::endl;
//							std::cout << "\t" << "tCon.info_.headPos_" << tCon.info_.tailPos_ << std::endl;
//							std::cout << "\t:pass: " << njh::colorBool((tCon.info_.headPos_ - hCon.info_.headPos_) == (nodes_[nodePos]->k_.size() - klen_ + 1)) << std::endl;
//							std::cout << std::endl;
//							if((tCon.info_.headPos_ - hCon.info_.headPos_) == (nodes_[nodePos]->k_.size() - klen_ + 1)) {
//								//paths[hCon.head_][tCon.tail_][hCon.info_.readName_].emplace_back(ConnectorInfoPath(hCon, tCon));
//							}
//						}
//					}
//				}
//				exit(1);
//			}
//			for(const auto & head : headConnectorInfos){
//
//
//				for(const auto pos : iter::range(head.second.size())){
//					std::cout << "positions difference              : " << tailConnectorInfos[head.first][pos].info_.headPos_ - head.second[pos].info_.headPos_ << std::endl;
//					std::cout << "nodes_[nodePos]->k_.size() - klen_: " << nodes_[nodePos]->k_.size() - klen_ + 1<< std::endl;
//					std::cout << std::endl;
//					if((tailConnectorInfos[head.first][pos].info_.headPos_ - head.second[pos].info_.headPos_) == (nodes_[nodePos]->k_.size() - klen_ + 1)) {
//						paths[head.second[pos].head_][tailConnectorInfos[head.first][pos].tail_][head.second[pos].info_.readName_].emplace_back(ConnectorInfoPath(head.second[pos], tailConnectorInfos[head.first][pos]));
//					}
//				}
//			}
			for(const auto & head : headConnectorInfos){
				for(const auto & hCon : head.second){
					for(const auto & tCon : tailConnectorInfos[head.first]){
						std::cout << "positions difference              : " << tCon.info_.headPos_ - hCon.info_.headPos_ << std::endl;
						std::cout << "nodes_[nodePos]->k_.size() - klen_: " << nodes_[nodePos]->k_.size() - klen_ + 1 << std::endl;
						std::cout << "\t" << "hCon.info_.headPos_" << hCon.info_.headPos_ << std::endl;
						std::cout << "\t" << "hCon.info_.headPos_" << hCon.info_.tailPos_ << std::endl;
						std::cout << "\t" << "tCon.info_.headPos_" << tCon.info_.headPos_ << std::endl;
						std::cout << "\t" << "tCon.info_.headPos_" << tCon.info_.tailPos_ << std::endl;
						std::cout << "\t:pass: " << njh::colorBool((tCon.info_.headPos_ - hCon.info_.headPos_) == (nodes_[nodePos]->k_.size() - klen_ + 1)) << std::endl;
						std::cout << std::endl;
						if((tCon.info_.headPos_ - hCon.info_.headPos_) == (nodes_[nodePos]->k_.size() - klen_ + 1)) {
							paths[hCon.head_][tCon.tail_][hCon.info_.readName_].emplace_back(ConnectorInfoPath(hCon, tCon));
						}
					}
				}
			}
			if("GATC" == nodes_[nodePos]->uid_){
				std::cout << __FILE__ << " " << __LINE__ << std::endl;
				//exit(1);
			}
			//std::cout << __FILE__ << " " << __LINE__ << std::endl;
			for(const auto & head : paths){
				std::cout << __FILE__ << " " << __LINE__ << std::endl;
				std::cout << "head: " << head.first << std::endl;
				std::cout << "      " << std::string(head.first.size() - klen_ + 1, ' ') << nodes_[nodePos]->k_ << std::endl;
				std::cout << __FILE__ << " " << __LINE__ << std::endl;
				for(const auto & tail : head.second){
					std::cout << "      " << std::string(head.first.size() - klen_ + 1 + (nodes_[nodePos]->k_.size() - klen_ + 1), ' ') << tail.first << std::endl;
					std::cout << "\ttail: " << std::endl;
					bool allOne = true;
					for(const auto & read : tail.second){
						if(1 != read.second.size()){
							allOne = false;
						}
						std::cout << "\t\tread: " << read.first << " " << read.second.size() << std::endl;
						std::cout << "\t\tread size: " << read.second.size() << " " << read.second.size() << std::endl;
						std::cout << "\t\tallOne: " << njh::colorBool(allOne) << std::endl;
						for(const auto & readheadPos : read.second){
							std::cout << "\t\t\t"<< readheadPos.head_.info_.headPos_ << ":" <<readheadPos.head_.info_.headPos_ + readheadPos.head_.head_.size() -klen_+ 1 << ":" << readheadPos.tail_.info_.tailPos_<< std::endl;
							std::cout << "\t\t\t" << "readheadPos.head_.info_.headPos_" << ":" <<readheadPos.head_.info_.headPos_ << std::endl;
							std::cout << "\t\t\t" << "readheadPos.head_.info_.tailPos_" << ":" <<readheadPos.head_.info_.tailPos_ << std::endl;
							std::cout << "\t\t\t" << "readheadPos.tail_.info_.headPos_" << ":" <<readheadPos.tail_.info_.headPos_ << std::endl;
							std::cout << "\t\t\t" << "readheadPos.tail_.info_.tailPos_" << ":" <<readheadPos.tail_.info_.tailPos_ << std::endl;
							std::cout << std::endl;;
						}
					}
					std::cout << "allOne: " << njh::colorBool(allOne) << std::endl;
//					if("AACAGAAGAAACAGAAGAAACAG" == nodes_[nodePos]->k_){
//						std::cout << njh::bashCT::reset;
//					}
//					if("GATC" == nodes_[nodePos]->uid_){
//						std::cout << __FILE__ << " " << __LINE__ << std::endl;
//						//exit(1);
//					}
					if(allOne){
						//std::cout << __FILE__ << " " << __LINE__ << std::endl;
						//grab head node and tail edges
						std::shared_ptr<edge> headEdge;
						for(const auto & headE  : nodes_[nodePos]->headEdges_){
							if(headE->head_.lock()->uid_ == head.first){
								headEdge = headE;
								break;
							}
						}
						if(nullptr == headEdge){
							std::stringstream ss;
							ss << __PRETTY_FUNCTION__ << ", error " << "couldn't find head node for " << head.first<< "\n";
							throw std::runtime_error{ss.str()};
						}
						std::shared_ptr<edge> tailEdge;
						for(const auto & tailE  : nodes_[nodePos]->tailEdges_){
							if(tailE->tail_.lock()->uid_ == tail.first){
								tailEdge = tailE;
								break;
							}
						}
						if(nullptr == tailEdge){
							std::stringstream ss;
							ss << __PRETTY_FUNCTION__ << ", error " << "couldn't find tail node for " << tail.first<< "\n";
							throw std::runtime_error{ss.str()};
						}
						auto headEdgeNode = headEdge->head_.lock();
						auto tailEdgeNode = tailEdge->tail_.lock();

						std::cout << "nodes_[nodePos]->visitCount_: " << nodes_[nodePos]->visitCount_ << std::endl;
						std::cout << "headEdgeNode->visitCount_   : " << headEdgeNode->visitCount_ << std::endl;
						std::cout << "tailEdgeNode->visitCount_   : " << tailEdgeNode->visitCount_ << std::endl;
//						if("GATC" == nodes_[nodePos]->uid_){
//							std::cout << __FILE__ << " " << __LINE__ << std::endl;
//							//exit(1);
//						}

						if(0 == nodes_[nodePos]->visitCount_ &&
							 0 == headEdgeNode->visitCount_ &&
							 0 == tailEdgeNode->visitCount_){
//						if(true){
							//can create a new node;
							newNodessAdded = true;
							auto newNode = std::make_shared<node>(nodes_[nodePos]->k_, nodes_[nodePos]->cnt_, klen_);
							newNode->uid_ = nodes_[nodePos]->uid_ + estd::to_string(nodes_.size());
							std::cout << njh::bashCT::purple;
							std::cout << "Adding new node " << newNode->uid_ << " from node position: " << nodePos << std::endl;
							std::cout << njh::bashCT::reset;
							nodePositions_[newNode->uid_] = nodes_.size();
							nodes_.emplace_back(newNode);

							//get the connector infos
							std::vector<ContigsCompareGraphDev::edge::ConnectorInfo> headConInfos;
							std::vector<ContigsCompareGraphDev::edge::ConnectorInfo> tailConInfos;
							for(const auto & read : tail.second){
								for(const auto & readheadPos : read.second){
									headConInfos.emplace_back(readheadPos.head_.info_);
									tailConInfos.emplace_back(readheadPos.tail_.info_);
									//add read names
									newNode->inReadNamesIdx_.emplace(readheadPos.head_.info_.readName_);
									newNode->inReadNamesIdx_.emplace(readheadPos.tail_.info_.readName_);
								}
							}

							//create new head node and tail edges
							//head
							auto newHeadEdge = std::make_shared<edge>(
									headEdge->head_.lock(),
									newNode,
									headConInfos.size(),
									headConInfos);
							//tail
							auto newTailEdge = std::make_shared<edge>(
									newNode,
									tailEdge->tail_.lock(),
									tailConInfos.size(),
									tailConInfos);
							//add edges
							headEdge->head_.lock()->tailEdges_.emplace_back(newHeadEdge);
							newNode->headEdges_.emplace_back(newHeadEdge);
							tailEdge->tail_.lock()->headEdges_.emplace_back(newTailEdge);
							newNode->tailEdges_.emplace_back(newTailEdge);
							//mark nodes as visited
							headEdge->head_.lock()->visitCount_ += 1;
							tailEdge->tail_.lock()->visitCount_ += 1;
							newNode->visitCount_ +=1;
							/**@todo do we need to remove the read names from the node? */
							//erase the connector nodes from the old edges, set edge to off if all connectors removed;
							//head
							std::vector<uint32_t> headConnectorPositionToRemove;
							for(const auto conPos : iter::range(headEdge->connectorInfos_.size())){
								for(const auto & currentCon : headConInfos){
									if(headEdge->connectorInfos_[conPos].sameInfo(currentCon)){
										headConnectorPositionToRemove.emplace_back(conPos);
										break;
									}
								}
							}
							std::sort(headConnectorPositionToRemove.rbegin(), headConnectorPositionToRemove.rend());
							for(const auto toRemovePos : headConnectorPositionToRemove){
								headEdge->connectorInfos_.erase(headEdge->connectorInfos_.begin() + toRemovePos);
							}
							if(headEdge->connectorInfos_.empty()){
								headEdge->on_ = false;
							}
							//tail
							std::vector<uint32_t> tailConnectorPositionToRemove;
							for(const auto conPos : iter::range(tailEdge->connectorInfos_.size())){
								for(const auto & currentCon : tailConInfos){
									if(tailEdge->connectorInfos_[conPos].sameInfo(currentCon)){
										tailConnectorPositionToRemove.emplace_back(conPos);
										break;
									}
								}
							}
							std::sort(tailConnectorPositionToRemove.rbegin(), tailConnectorPositionToRemove.rend());
							for(const auto toRemovePos : tailConnectorPositionToRemove){
								tailEdge->connectorInfos_.erase(tailEdge->connectorInfos_.begin() + toRemovePos);
							}
							if(tailEdge->connectorInfos_.empty()){
								tailEdge->on_ = false;
							}
						}
					}else{
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ", error " << "need to increase kmer size, can't disentagle non-redundant kmers at this size" << "\n";
						throw std::runtime_error{ss.str()};
						//exit(1);
					}
				}
			}

			//if all the edges have been turned off then this node can be turned off
			bool allOff = true;
			std::set<std::string> updateNames;
			for(const auto & tail : nodes_[nodePos]->tailEdges_){
				if(tail->on_){
					for(const auto & con : tail->connectorInfos_){
						updateNames.emplace(con.readName_);
					}
					allOff = false;
				}
			}
			for(const auto & head : nodes_[nodePos]->headEdges_){
				if(head->on_){
					allOff = false;
					for(const auto & con : head->connectorInfos_){
						updateNames.emplace(con.readName_);
					}
				}
			}
			if(allOff){
				nodes_[nodePos]->on_ = false;
			}else{
				//make sure the node only has names from the connections still on
				nodes_[nodePos]->inReadNamesIdx_ = updateNames;
			}

		}
		if(newNodessAdded){
			removeOffNodesOffEdgesAndReset();
		}
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		return newNodessAdded;
	}
	//std::cout << __FILE__ << " " << __LINE__ << std::endl;
	return false;
}



void ContigsCompareGraphDev::threadThroughSequence(
		const seqInfo & seq, const std::string & threadingSeqName) {
	//std::cout << __FILE__ << " " << __LINE__ << std::endl;
	if(seq.seq_.size() > klen_){
		std::string startKmer = seq.seq_.substr(0, klen_);
		auto lastPosIter = nodePositions_.find(startKmer);
		for (const auto pos : iter::range(seq.seq_.size() - klen_)) {
			std::string firstKmer = seq.seq_.substr(pos, klen_);
			std::string nextKmer = seq.seq_.substr(pos + 1, klen_);
			auto firstNodePosition = lastPosIter;
			auto nextNodePosition =  nodePositions_.find(nextKmer);
			//check if both nodes pass occurrence check
			bool passOccurenceCheck = nodePositions_.end() != firstNodePosition && nodePositions_.end() != nextNodePosition;
			//set the last position iterator as the next node iterator
			lastPosIter = nextNodePosition;
			if(passOccurenceCheck){
				auto headNode = nodes_[firstNodePosition->second];
				auto tailNode = nodes_[nextNodePosition->second];
				headNode->inReadNamesIdx_.emplace(threadingSeqName);
				tailNode->inReadNamesIdx_.emplace(threadingSeqName);
				bool foundEdge = false;
				edge::ConnectorInfo info(threadingSeqName, pos, pos+1);
				for(const auto & tail : headNode->tailEdges_){
					auto tailNode = tail->tail_.lock();
					if(nextKmer == tailNode->k_){
						foundEdge = true;
						tail->cnt_ += 1;
						tail->connectorInfos_.emplace_back(info);
						break;
					}
				}
				if(!foundEdge){
					std::shared_ptr<edge> e = std::make_shared<edge>(headNode, tailNode, 1, std::vector<edge::ConnectorInfo>{info});
					headNode->addTail(e);
					tailNode->addHead(e);
				}
			}
		}
	}
	//std::cout << __FILE__ << " " << __LINE__ << std::endl;
}


std::unordered_map<std::string, std::string> ContigsCompareGraphDev::writeRectangleDotColorBySamples(std::ostream & out, bool noLabels) const{
	VecStr moreColors = {"#006E82","#8214A0","#005AC8","#00A0FA","#FA78FA","#14D2DC","#AA0A3C","#FA7850","#0AB45A","#F0F032","#A0FA82","#FAE6BE",
			"#ff358f", "#07c652", "#a2009b", "#467f00", "#cb73fd",
			"#f4be45", "#0157d8", "#ff8d3b", "#0056a1", "#dd0d3a", "#01d5e4",
			"#b10049", "#7cda97", "#ff76e0", "#018a5a", "#ff87b8", "#4a5b00",
			"#664092", "#8f7400", "#02aee5", "#9e3500", "#8bd5b8", "#8a306b",
			"#e4b7f0", "#ff97a2" };
	out << "digraph graphname {" << std::endl;
	out << "\t" << "node [fixedsize=true, shape=rect]" << std::endl;
	double heightNormalizer = 250; // 250 bases will equate 1 inch
//	double widthNormalizer = 25; // per base coverage of 25 will equate 1 inch
	//double penWidthNormalizer = 100; // 50 will equate 1 inch
	double maxPenWidth = 12.5; // 50 will equate 1 inch, 25 == .5 inch
	std::vector<double> avgBaseCoverages;
	uint32_t maxSampleSize = 0;
	std::set<std::string> sampleCombos;
	for(const auto & node : nodes_){
		if(!node->on_){
			continue;
		}
		//double approxPerBaseCoverage = node->cnt_/(node->k_.length() - klen_ + 1);
		avgBaseCoverages.emplace_back(node->inReadNamesIdx_.size());
		sampleCombos.emplace(njh::conToStr(node->inReadNamesIdx_, ","));
	}

	maxSampleSize = sampleCombos.size();
	std::unordered_map<std::string, uint32_t> sampleComboToIndex;
	uint32_t sampleComboIndex = 0;
	for(const auto & sampCombo : sampleCombos){
		sampleComboToIndex[sampCombo] = sampleComboIndex++;
	}

	if(maxSampleSize > moreColors.size()){
		auto hColors = njh::heatmapColors(maxSampleSize + 1);
		njh::reverse(hColors);
		moreColors.clear();
		for(const auto & hColor : hColors){
			moreColors.emplace_back(hColor.getHexStr());
		}

		out << "\t"
				<< "graph [ bgcolor=black, resolution=128, fontname=Arial, fontcolor=white,  fontsize=12 ]; "
				<< std::endl;
		out << "\t" << "node [ fontname=Arial, fontcolor=white, fontsize=11];"
				<< std::endl;
		out << "\t" << "edge [ fontname=Helvetica, fontcolor=white, fontsize=10 ];"
				<< std::endl;
	}

	std::unordered_map<std::string, std::string> ret;
	for(const auto & sampCombo : sampleComboToIndex){
		ret[sampCombo.first] = moreColors.at(sampCombo.second);
	}

	auto maxCov = *std::max_element(avgBaseCoverages.begin(), avgBaseCoverages.end());

	uint32_t nodeIndex = 0;
	for(const auto & node : nodes_){
		if(!node->on_){
			continue;
		}

		uint32_t colorGroup = 0;
		std::unordered_set<std::string> sampleNames;
		for(const auto & name : node->inReadNamesIdx_){
			if(MetaDataInName::nameHasMetaData(name)){
				MetaDataInName meta(name);
				if(meta.containsMeta("sample")){
					sampleNames.emplace(meta.getMeta("sample"));
				}
			}
		}
		colorGroup = sampleNames.size();
		std::string nodeColor = "";
		nodeColor = moreColors.at(sampleComboToIndex.at(njh::conToStr(node->inReadNamesIdx_, ",")));

		double nheight = node->k_.length()/heightNormalizer;
		double approxPerBaseCoverage = node->inReadNamesIdx_.size();
		double nwidth = (approxPerBaseCoverage / maxCov) * 5;
		out << "\t" << node->uid_ << "[fixedsize=true,shape=rect,width=" << nwidth
				<< ",height=" << nheight << ",style=filled,fillcolor=\"" << nodeColor
//				<< "\", label=\""
//				<< (noLabels ?
//						"" :
//						njh::pasteAsStr("Len=", node->k_.size(), ";Cov=", roundDecPlaces(approxPerBaseCoverage, 2), ";", "sampleCount=", sampleNames.size(), ";"))
						<< "\", label=\""
						<< (noLabels ?
								"" :
								njh::pasteAsStr("Len=", node->k_.size(), ";\nCov=", roundDecPlaces(approxPerBaseCoverage, 2), ";\n", "sampleCount=", sampleNames.size(), ";\nnodeIndex:", nodeIndex))
				<< "\"]" << std::endl;
		;
		++nodeIndex;
	}

	uint32_t maxTailCnt = 0;
	for (const auto & node : nodes_){
		if(!node->on_){
			continue;
		}
		for(const auto & tail : node->tailEdges_){
			if(tail->on_){
				if(tail->cnt_ > maxTailCnt){
					maxTailCnt = tail->cnt_;
				}
			}
		}
	}

	scale<double> edgePenWidthScale(std::make_pair(0,maxTailCnt), std::make_pair(0,maxPenWidth));

	for (const auto & node : nodes_){
		if(!node->on_){
			continue;
		}
		for(const auto & tail : node->tailEdges_){
			if(tail->on_){
				std::set<std::string> sampleNamesTail;
				for(const auto & name : tail->connectorInfos_){
					sampleNamesTail.emplace(name.readName_);
				}

				out << "\t" << tail->head_.lock()->uid_
						<< " -> " << tail->tail_.lock()->uid_
						<< "[color=\""<< moreColors.at(sampleComboToIndex[njh::conToStr(sampleNamesTail, ",")])
						<<  "\", penwidth=" << edgePenWidthScale.get(tail->connectorInfos_.size())
						<< ", label=\"" << (noLabels ? "" : njh::pasteAsStr("cons:",tail->connectorInfos_.size(), ":samps", sampleNamesTail.size())) << "\""<< "]"   << std::endl;
			}
		}
	}
	out << "}" << std::endl;
	return ret;
}

void ContigsCompareGraphDev::writeRectangleDotColorBySampleCount(std::ostream & out, bool noLabels) const{
	VecStr moreColors = {"#006E82","#8214A0","#005AC8","#00A0FA","#FA78FA","#14D2DC","#AA0A3C","#FA7850","#0AB45A","#F0F032","#A0FA82","#FAE6BE",
			"#ff358f", "#07c652", "#a2009b", "#467f00", "#cb73fd",
			"#f4be45", "#0157d8", "#ff8d3b", "#0056a1", "#dd0d3a", "#01d5e4",
			"#b10049", "#7cda97", "#ff76e0", "#018a5a", "#ff87b8", "#4a5b00",
			"#664092", "#8f7400", "#02aee5", "#9e3500", "#8bd5b8", "#8a306b",
			"#e4b7f0", "#ff97a2" };
	out << "digraph graphname {" << std::endl;
	out << "\t" << "node [fixedsize=true, shape=rect]" << std::endl;
	double heightNormalizer = 250; // 250 bases will equate 1 inch
//	double widthNormalizer = 25; // per base coverage of 25 will equate 1 inch
	//double penWidthNormalizer = 100; // 50 will equate 1 inch
	double maxPenWidth = 12.5; // 50 will equate 1 inch, 25 == .5 inch
	std::vector<double> avgBaseCoverages;
	uint32_t maxSampleSize = 0;
	for(const auto & node : nodes_){
		if(!node->on_){
			continue;
		}
		//double approxPerBaseCoverage = node->cnt_/(node->k_.length() - klen_ + 1);
		avgBaseCoverages.emplace_back(node->inReadNamesIdx_.size());
		std::unordered_set<std::string> sampleNames;
		for(const auto & name : node->inReadNamesIdx_){
			if(MetaDataInName::nameHasMetaData(name)){
				MetaDataInName meta(name);
				if(meta.containsMeta("sample")){
					sampleNames.emplace(meta.getMeta("sample"));
				}
			}
			if(sampleNames.size() > maxSampleSize){
				maxSampleSize = sampleNames.size();
			}
		}

		for(const auto & tail : node->tailEdges_){
			std::unordered_set<std::string> sampleNamesTail;
			for(const auto & name : tail->connectorInfos_){
				if(MetaDataInName::nameHasMetaData(name.readName_)){
					MetaDataInName meta(name.readName_);
					if(meta.containsMeta("sample")){
						sampleNamesTail.emplace(meta.getMeta("sample"));
					}
				}
			}
			if(sampleNamesTail.size() > maxSampleSize){
				maxSampleSize = sampleNamesTail.size();
			}
		}
	}
	if(maxSampleSize > moreColors.size()){
		auto hColors = njh::heatmapColors(maxSampleSize + 1);
		njh::reverse(hColors);
		moreColors.clear();
		for(const auto & hColor : hColors){
			moreColors.emplace_back(hColor.getHexStr());
		}

		out << "\t"
				<< "graph [ bgcolor=black, resolution=128, fontname=Arial, fontcolor=white,  fontsize=12 ]; "
				<< std::endl;
		out << "\t" << "node [ fontname=Arial, fontcolor=white, fontsize=11];"
				<< std::endl;
		out << "\t" << "edge [ fontname=Helvetica, fontcolor=white, fontsize=10 ];"
				<< std::endl;
	}


	auto maxCov = *std::max_element(avgBaseCoverages.begin(), avgBaseCoverages.end());


	for(const auto & node : nodes_){
		if(!node->on_){
			continue;
		}
		uint32_t colorGroup = 0;
		std::unordered_set<std::string> sampleNames;
		for(const auto & name : node->inReadNamesIdx_){
			if(MetaDataInName::nameHasMetaData(name)){
				MetaDataInName meta(name);
				if(meta.containsMeta("sample")){
					sampleNames.emplace(meta.getMeta("sample"));
				}
			}
		}
		colorGroup = sampleNames.size();
		std::string nodeColor = "";
		nodeColor = moreColors.at(colorGroup);

		double nheight = node->k_.length()/heightNormalizer;
		double approxPerBaseCoverage = node->inReadNamesIdx_.size();
		double nwidth = (approxPerBaseCoverage / maxCov) * 5;
		out << "\t" << node->uid_ << "[fixedsize=true,shape=rect,width=" << nwidth
				<< ",height=" << nheight << ",style=filled,fillcolor=\"" << nodeColor
//				<< "\", label=\""
//				<< (noLabels ?
//						"" :
//						njh::pasteAsStr("Len=", node->k_.size(), ";Cov=", roundDecPlaces(approxPerBaseCoverage, 2), ";", "sampleCount=", sampleNames.size(), ";"))
						<< "\", label=\""
						<< (noLabels ?
								"" :
								njh::pasteAsStr("Len=", node->k_.size(), ";\nCov=", roundDecPlaces(approxPerBaseCoverage, 2), ";\n", "sampleCount=", sampleNames.size(), ";"))
				<< "\"]" << std::endl;
		;
	}

	uint32_t maxTailCnt = 0;
	for (const auto & node : nodes_){
		if(!node->on_){
			continue;
		}
		for(const auto & tail : node->tailEdges_){
			if(tail->on_){
				if(tail->cnt_ > maxTailCnt){
					maxTailCnt = tail->cnt_;
				}
			}
		}
	}

	scale<double> edgePenWidthScale(std::make_pair(0,maxTailCnt), std::make_pair(0,maxPenWidth));

	for (const auto & node : nodes_){
		if(!node->on_){
			continue;
		}
		for(const auto & tail : node->tailEdges_){
			if(tail->on_){
				uint32_t colorGroup = 0;
				std::unordered_set<std::string> sampleNamesTail;
				for(const auto & name : tail->connectorInfos_){
					if(MetaDataInName::nameHasMetaData(name.readName_)){
						MetaDataInName meta(name.readName_);
						if(meta.containsMeta("sample")){
							sampleNamesTail.emplace(meta.getMeta("sample"));
						}
					}
				}
				colorGroup = sampleNamesTail.size();
				out << "\t" << tail->head_.lock()->uid_
						<< " -> " << tail->tail_.lock()->uid_
						<< "[color=\""<< moreColors.at(colorGroup)
						<<  "\", penwidth=" << edgePenWidthScale.get(tail->connectorInfos_.size())
						<< ", label=\"" << (noLabels ? "" : njh::pasteAsStr("cons:",tail->connectorInfos_.size(), ":samps", sampleNamesTail.size())) << "\""<< "]"   << std::endl;
			}
		}
	}
	out << "}" << std::endl;
}


void ContigsCompareGraphDev::writeRectangleDot(std::ostream & out, bool noLabels) const{
	VecStr colors = {"#006E82","#8214A0","#005AC8","#00A0FA","#FA78FA","#14D2DC","#AA0A3C","#FA7850","#0AB45A","#F0F032","#A0FA82","#FAE6BE"};
	VecStr moreColors = { "#ff358f", "#07c652", "#a2009b", "#467f00", "#cb73fd",
			"#f4be45", "#0157d8", "#ff8d3b", "#0056a1", "#dd0d3a", "#01d5e4",
			"#b10049", "#7cda97", "#ff76e0", "#018a5a", "#ff87b8", "#4a5b00",
			"#664092", "#8f7400", "#02aee5", "#9e3500", "#8bd5b8", "#8a306b",
			"#e4b7f0", "#ff97a2" };
	out << "digraph graphname {" << std::endl;
	out << "\t" << "node [fixedsize=true, shape=rect]" << std::endl;
	double heightNormalizer = 250; // 250 bases will equate 1 inch
//	double widthNormalizer = 25; // per base coverage of 25 will equate 1 inch
	double penWidthNormalizer = 50; // 50 will equate 1 inch
	std::vector<double> avgBaseCoverages;
	for(const auto & node : nodes_){
		if(!node->on_){
			continue;
		}
		double approxPerBaseCoverage = node->cnt_/(node->k_.length() - klen_ + 1);
		avgBaseCoverages.emplace_back(approxPerBaseCoverage);
	}
	auto maxCov = *std::max_element(avgBaseCoverages.begin(), avgBaseCoverages.end());
	for(const auto & node : nodes_){
		if(!node->on_){
			continue;
		}
		uint32_t colorGroup = 0;
		if(node->headless() && node->tailless()){
			colorGroup = 1;
		}else if(node->headless()){
			colorGroup = 0;
		}else if(node->tailless()){
			colorGroup = 2;
		}else{
			colorGroup = 2 + node->tailCount();
		}
		std::string nodeColor = "";
		if(colorGroup >= colors.size()){
			if(colorGroup >= (colors.size() + moreColors.size())){
				std::stringstream ss;
				ss << "Error in " << __PRETTY_FUNCTION__ << "\n";
				ss << "Oh no, not enough colors!" << "\n";
				ss << "Need " << colorGroup << " but only have " << colors.size() + moreColors.size() << "\n";
				throw std::runtime_error{ss.str()};
			}else{
				nodeColor = moreColors.at(colorGroup - colors.size());
			}
		}else{
			nodeColor = colors.at(colorGroup);
		}
		std::unordered_set<std::string> sampleNames;
		for(const auto & name : node->inReadNamesIdx_){
			if(MetaDataInName::nameHasMetaData(name)){
				MetaDataInName meta(name);
				if(meta.containsMeta("sample")){
					sampleNames.emplace(meta.getMeta("sample"));
				}
			}
		}
		double nheight = node->k_.length()/heightNormalizer;
		double approxPerBaseCoverage = static_cast<double>(node->cnt_)/(node->k_.length() - klen_ + 1);
		double nwidth = (approxPerBaseCoverage / maxCov) * 5;
		out << "\t" << node->uid_ << "[fixedsize=true,shape=rect,width=" << nwidth
				<< ",height=" << nheight << ",style=filled,fillcolor=\"" << nodeColor
				<< "\", label=\""
				<< (noLabels ?
						"" :
						njh::pasteAsStr("Len=", node->k_.size(), ";\nCov=", roundDecPlaces(approxPerBaseCoverage, 2), ";\n",
								"sampleCount=", sampleNames.size(), ";"))
//				<< (noLabels ?
//						"" :
//						njh::pasteAsStr("Len=", node->k_.size(), ";\nCnt=", node->cnt_,
//								";\nCov=", roundDecPlaces(approxPerBaseCoverage, 2), "\n"))
				<< "\"]" << std::endl;
		;
	}
	for (const auto & node : nodes_){
		if(!node->on_){
			continue;
		}
		for(const auto & tail : node->tailEdges_){
			if(tail->on_){
				out << "\t" << tail->head_.lock()->uid_ << " -> " << tail->tail_.lock()->uid_ << "[penwidth=" << tail->cnt_/penWidthNormalizer << ", label=\"" << (noLabels ? "" : estd::to_string(tail->cnt_)) << "\""<< "]"   << std::endl;
			}
		}
	}
	out << "}" << std::endl;
}




void ContigsCompareGraphDev::collapseSingleLinkedPaths(){
	std::vector<uint32_t> nodePositionsToProcess;
	for (const auto nPos : iter::range(nodes_.size())) {
		if (1 == nodes_[nPos]->tailCount()) {
			nodePositionsToProcess.push_back(nPos);
		}
	}

	for(auto & nPos : nodePositionsToProcess){
		if(nullptr == nodes_[nPos]){
			continue;
		}
		nodes_[nPos]->visitCount_ += 1;
		std::vector<uint32_t> nodesToErase;
		auto next = nodes_[nPos]->getFirstOnTailEdge()->tail_.lock();
		//a tail count of one indicates that the next node should be added in
		//also head count should be less than 2
		while(next != nullptr && next->headCount() == 1 && nodes_[nPos]->uid_ != next->uid_){
			next->visitCount_ += 1;
			nodes_[nPos]->k_.append(next->k_.substr(klen_ - 1 ));
			nodes_[nPos]->cnt_+= next->cnt_;
			nodes_[nPos]->inReadNamesIdx_.insert(next->inReadNamesIdx_.begin(), next->inReadNamesIdx_.end());
			auto toErase = next->uid_;
			if (0 == next->tailCount()) {
				//head a tailless node, end of the line
				//need to erase tail edges
				nodes_[nPos]->tailEdges_.erase(nodes_[nPos]->tailEdges_.begin());
				next = nullptr;
			} else if (next->tailCount() > 1) {
				//head a multi tailed node, section collapse is done
				//need to add all the multiple tail edges
				nodes_[nPos]->tailEdges_.clear();
				for (const auto & tail : next->tailEdges_) {
					if(!tail->on_){
						continue;
					}
					nodes_[nPos]->tailEdges_.push_back(tail);
					nodes_[nPos]->tailEdges_.back()->head_ = nodes_[nPos];
				}
				next->tailEdges_.clear();
				next = nullptr;
			} else {
				//can still collapse more set next to the tail edge
				//will only hit here if tailEdges_.size() is 1
				nodes_[nPos]->tailEdges_.clear();
				for (const auto & tail : next->tailEdges_) {
					if (!tail->on_) {
						continue;
					}
					nodes_[nPos]->tailEdges_.push_back(tail);
					nodes_[nPos]->tailEdges_.back()->head_ = nodes_[nPos];
				}
				next = next->getFirstOnTailEdge()->tail_.lock();
			}
			nodesToErase.emplace_back(nodePositions_[toErase]);
		}
		if(!nodesToErase.empty()){
			std::sort(nodesToErase.rbegin(), nodesToErase.rend());
			for(const auto & remove : nodesToErase){
				nodes_[remove] = nullptr;
			}
		}
	}
	removeNullNodes();
}



bool ContigsCompareGraphDev::collapseLowFreqNodes(const comparison & allowableError, uint32_t lowFreqCutOff){
	std::map<std::string, std::vector<std::shared_ptr<node>>> groupedNodes;
	//gather together the nodes that have the same heads and tails and are also singlely linked forward and backwards
	for(const auto & n : nodes_){
		VecStr tailConnectorNames;
		VecStr headConnectorNames;
		if(( 1 == n->headCount() && 1== n->tailCount() ) ||
				(1 == n->headCount() && n->tailless() )      ||
				(n->headless()       && 1== n->tailCount() ) ){
			for(const auto & h : n->headEdges_){
				if(h->on_){
					headConnectorNames.emplace_back(h->head_.lock()->uid_);
				}
			}
			for(const auto & t : n->tailEdges_){
				if(t->on_){
					tailConnectorNames.emplace_back(t->tail_.lock()->uid_);
				}
			}
		}
		if(!headConnectorNames.empty() || !tailConnectorNames.empty() ){
			njh::sort(headConnectorNames);
			njh::sort(tailConnectorNames);
			groupedNodes[njh::pasteAsStr("head:",njh::conToStr(headConnectorNames, ","), ";tail:", njh::conToStr(tailConnectorNames, ","))].emplace_back(n);
		}
	}
	bool modifiedNodes = false;
	//iterate over the groups
	uint32_t groupCount = 0;
	for(auto & group : groupedNodes){
		// if there are only two nodes in the group, process them
		if(2 == group.second.size()){
			++groupCount;
			if(uAbsdiff(group.second.front()->k_.length(), group.second.back()->k_.length()) > 20){
				continue;
			}

			//sort the two nodes so that the front node is the higher read count than the 2nd one
			njh::sort(group.second, [](const std::shared_ptr<node> & node1, const std::shared_ptr<node> & node2){
				if(node1->inReadNamesIdx_.size() == node2->inReadNamesIdx_.size()){
					return node1->k_.length() < node2->k_.length();
				}else{
					return node1->inReadNamesIdx_.size() > node2->inReadNamesIdx_.size();
				}
			});
			auto node1 = group.second.front();
			auto node2 = group.second.back();


			//only process if node1 is above the cut off and if the second node is below the cut off
			if(node2->inReadNamesIdx_.size() > lowFreqCutOff || node1->inReadNamesIdx_.size() <= lowFreqCutOff){
				continue;
			}
			//compare the two nodes
			aligner alignerObj(std::max(node1->k_.length(), node2->k_.length()), gapScoringParameters(5,1,5,1,5,1));

			alignerObj.countEndGaps_ = true;
			alignerObj.weighHomopolymers_ = true;
			alignerObj.alignRegGlobal(seqInfo("1", node1->k_),
					seqInfo("2", node2->k_));
			alignerObj.profileAlignment(
					seqInfo("1", node1->k_),
					seqInfo("2", node2->k_), false, false, false);
//					//uncomment for debugging purposes
//					alignerObj.alignObjectA_.seqBase_.outPutSeqAnsi(std::cout);
//					alignerObj.alignObjectB_.seqBase_.outPutSeqAnsi(std::cout);
//					std::cout << "allowableErrorForHPIndexCollapse_.passErrorProfile(alignerObj.comp_): " << njh::colorBool(allowableErrorForHPIndexCollapse_.passErrorProfile(alignerObj.comp_)) << std::endl;
			if(allowableError.passErrorProfile(alignerObj.comp_)){
				modifiedNodes = true;
				//add in the other node's read names
				node1->inReadNamesIdx_.insert(
						node2->inReadNamesIdx_.begin(),
						node2->inReadNamesIdx_.end());
				if(!node1->headless()){
					auto node1HeadEdge = node1->getFirstOnHeadEdge();
					auto node2HeadEdge = node2->getFirstOnHeadEdge();
					node1HeadEdge->cnt_ += node2HeadEdge->cnt_;
					node1HeadEdge->connectorInfos_.insert(
							node1HeadEdge->connectorInfos_.end(),
							node2HeadEdge->connectorInfos_.begin(),
							node2HeadEdge->connectorInfos_.end());
				}
				if(!node1->tailless()){
					auto node1TailEdge = node1->getFirstOnTailEdge();
					auto node2TailEdge = node2->getFirstOnTailEdge();
					node1TailEdge->cnt_ += node2->getFirstOnTailEdge()->cnt_;
					node1TailEdge->connectorInfos_.insert(
							node1TailEdge->connectorInfos_.end(),
							node2TailEdge->connectorInfos_.begin(),
							node2TailEdge->connectorInfos_.end());
				}
				// now turn off the node and it's edges
				node2->on_ = false;
				if(!node2->headless()){
					node2->getFirstOnHeadEdge()->on_ = false;
				}
				if(!node2->tailless()){
					node2->getFirstOnTailEdge()->on_ = false;
				}
			}
		}
	}
//	//
	if (modifiedNodes) {
		//remove off nodes;
		removeOffNodes();
		removeOffEdges();
		resetNodePositions();
		return true;
	}
	return false;
}




void ContigsCompareGraphDev::collapseSingleLinkedPathsSameReads(){
	std::vector<uint32_t> nodePositionsToProcess;
	for (const auto nPos : iter::range(nodes_.size())) {
		if (1 == nodes_[nPos]->tailCount()) {
			nodePositionsToProcess.push_back(nPos);
		}
	}

	for(auto & nPos : nodePositionsToProcess){
		if(nullptr == nodes_[nPos]){
			continue;
		}
		nodes_[nPos]->visitCount_ += 1;
		std::vector<uint32_t> nodesToErase;
		auto next = nodes_[nPos]->getFirstOnTailEdge()->tail_.lock();
		//a tail count of one indicates that the next node should be added in
		//also head count should be less than 2
		bool containsNext = true;
		if (next != nullptr) {
			for (const auto & nextRead : next->inReadNamesIdx_) {
				if (!njh::in(nextRead, nodes_[nPos]->inReadNamesIdx_)) {
					containsNext = false;
					break;
				}
			}
			for (const auto & nextRead : nodes_[nPos]->inReadNamesIdx_ ) {
				if (!njh::in(nextRead, next->inReadNamesIdx_)) {
					containsNext = false;
					break;
				}
			}
		}
		while(next != nullptr && next->headCount() == 1 && containsNext && nodes_[nPos]->uid_ != next->uid_){
//			if("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT" == next->k_){
//				std::cout << __FILE__ << " " << __LINE__ << std::endl;
//				std::cout << " next->inReadNamesIdx_.size(): " << next->inReadNamesIdx_.size() << std::endl;
//				std::cout << "nodes_[nPos]->inReadNamesIdx_.size(): " << nodes_[nPos]->inReadNamesIdx_.size() << std::endl;
//				std::cout << "nodes_[nPos]->uid_: " << nodes_[nPos]->uid_ << std::endl;
//			}
//			if("ATATTTTAACATAAAAATAAGTACATATTT" == next->k_){
//				std::cout << __FILE__ << " " << __LINE__ << std::endl;
//				std::cout << " next->inReadNamesIdx_.size(): " << next->inReadNamesIdx_.size() << std::endl;
//				std::cout << "nodes_[nPos]->inReadNamesIdx_.size(): " << nodes_[nPos]->inReadNamesIdx_.size() << std::endl;
//				std::cout << "nodes_[nPos]->uid_: " << nodes_[nPos]->uid_ << std::endl;
//			}
//
//
//			if(next->k_.find("ATATTTTAACATAAAAATAAGTACATATTT") != std::string::npos ){
//				std::cout << __FILE__ << " " << __LINE__ << std::endl;
//				std::cout << " next->inReadNamesIdx_.size(): " << next->inReadNamesIdx_.size() << std::endl;
//				std::cout << "nodes_[nPos]->inReadNamesIdx_.size(): " << nodes_[nPos]->inReadNamesIdx_.size() << std::endl;
//				std::cout << "nodes_[nPos]->uid_: " << nodes_[nPos]->uid_ << std::endl;
//			}


			next->visitCount_ += 1;
			nodes_[nPos]->k_.append(next->k_.substr(klen_ - 1 ));
			nodes_[nPos]->cnt_+= next->cnt_;
			nodes_[nPos]->inReadNamesIdx_.insert(next->inReadNamesIdx_.begin(), next->inReadNamesIdx_.end());
			auto toErase = next->uid_;
			if (0 == next->tailCount()) {
				//head a tailless node, end of the line
				//need to erase tail edges
				nodes_[nPos]->tailEdges_.erase(nodes_[nPos]->tailEdges_.begin());
				next = nullptr;
			} else if (next->tailCount() > 1) {
				//head a multi tailed node, section collapse is done
				//need to add all the multiple tail edges
				nodes_[nPos]->tailEdges_.clear();
				for (const auto & tail : next->tailEdges_) {
					if(!tail->on_){
						continue;
					}
					nodes_[nPos]->tailEdges_.push_back(tail);
					nodes_[nPos]->tailEdges_.back()->head_ = nodes_[nPos];
				}
				next->tailEdges_.clear();
				next = nullptr;
			} else {
				//can still collapse more set next to the tail edge
				//will only hit here if tailEdges_.size() is 1
				nodes_[nPos]->tailEdges_.clear();
				for (const auto & tail : next->tailEdges_) {
					if (!tail->on_) {
						continue;
					}
					nodes_[nPos]->tailEdges_.push_back(tail);
					nodes_[nPos]->tailEdges_.back()->head_ = nodes_[nPos];
				}
				next = next->getFirstOnTailEdge()->tail_.lock();
				if (next != nullptr) {
					containsNext = true;
					for(const auto & nextRead : next->inReadNamesIdx_){
						if(!njh::in(nextRead, nodes_[nPos]->inReadNamesIdx_)){
							containsNext = false;
							break;
						}
					}
					for (const auto & nextRead : nodes_[nPos]->inReadNamesIdx_ ) {
						if (!njh::in(nextRead, next->inReadNamesIdx_)) {
							containsNext = false;
							break;
						}
					}
				}
			}
			nodesToErase.emplace_back(nodePositions_[toErase]);
		}
		if(!nodesToErase.empty()){
			std::sort(nodesToErase.rbegin(), nodesToErase.rend());
			for(const auto & remove : nodesToErase){
				nodes_[remove] = nullptr;
			}
		}
	}
	removeNullNodes();
}


void ContigsCompareGraphDev::removeNullNodes(){
	std::vector<uint32_t> toRemove;
	for(const auto nodePos : iter::range(nodes_.size())){
		if(nullptr == nodes_[nodePos]){
			toRemove.emplace_back(nodePos);
		}
	}
	if(!toRemove.empty()){
		std::sort(toRemove.rbegin(), toRemove.rend());
		for(const auto & remove : toRemove){
			nodes_.erase(nodes_.begin() + remove);
		}
		resetNodePositions();
	}
}

std::vector<seqInfo> ContigsCompareGraphDev::nodesToSeqs() const{
	std::vector<seqInfo> ret;
	for(const auto & node : nodes_){
		std::string seq = node->k_;
		ret.emplace_back(node->uid_, seq);
	}
	return ret;
}


void ContigsCompareGraphDev::resetGroups() const {
	//first reset all nodes back to max
	for(auto & n : nodes_){
		n->group_ = std::numeric_limits<uint32_t>::max();
	}
	std::vector<std::shared_ptr<node>> nodesToProcess;
	for(const auto & n : nodes_){
		if(n->headless()){
			nodesToProcess.emplace_back(n);
		}
	}
	std::function<void(const std::shared_ptr<node> &, uint32_t)> spreadGroup = [&spreadGroup](const std::shared_ptr<node> & spnode, uint32_t groupUID){
		//set self
		spnode->group_ = groupUID;
		//spread to heads
		for(const auto & head : spnode->headEdges_){
			if(head->on_){
				if(std::numeric_limits<uint32_t>::max() == head->head_.lock()->group_){
					spreadGroup(head->head_.lock(), groupUID);
				}
			}
		}
		//spread to tails
		for(const auto & tail : spnode->tailEdges_){
			if(tail->on_){
				if(std::numeric_limits<uint32_t>::max() == tail->tail_.lock()->group_){
					spreadGroup(tail->tail_.lock(), groupUID);
				}
			}
		}
	};
	uint32_t groupId = 0;
	for(const auto & n : nodesToProcess){
		//group hasn't been set yet
		if(std::numeric_limits<uint32_t>::max() == n->group_){
			spreadGroup(n, groupId);
			++groupId;
		}
	}
}




std::vector<seqInfo> ContigsCompareGraphDev::correctSeqsByGraph(const std::vector<seqInfo> & seqs, const correctSeqsByGraphPars & pars) {

	njh::stopWatch watch;
	//correction
	ContigsCompareGraphDev compGraph(pars.klen);
	compGraph.setOccurenceCutOff(pars.kmerOccurenceCutOff);
	for(const auto & seq : seqs){
		compGraph.increaseKCounts(seq.seq_);
	}
	compGraph.populateNodesFromCounts();
	for(const auto & seq : seqs){
		compGraph.threadThroughSequence(seq, seq.name_);
	}


	std::unordered_map<std::string, uint32_t> seqnameToPosition;
	{
		uint32_t seqCount = 0;
		for(const auto & seq : seqs){
			seqnameToPosition[seq.name_] = seqCount++;
		}
	}
	std::vector<seqInfo> correctedSeqs = seqs;
	for(const auto & seq : seqs){
		bool printInfo = pars.debug;
		std::vector<SubSeqment> subPositions;
		for(uint32_t pos = 0; pos < len(seq) + 1 - pars.klen; ++pos){
			std::string subK = seq.seq_.substr(pos, pars.klen);
			if(compGraph.kCounts_[subK] <= pars.correctionOccurenceCutOff){
				if(subPositions.empty() || (subPositions.back().pos_ + subPositions.back().size_ + 1 - pars.klen != pos)){
					subPositions.emplace_back(SubSeqment(pos, pars.klen));
				}else{
					subPositions.back().size_ += 1;
				}
			}
		}
		if(printInfo){
			std::cout << "subPositions.size(): " << subPositions.size() << std::endl;
		}
		if(!subPositions.empty()){
			if(printInfo){
				std::cout << seq.name_ << std::endl;
			}
			njh::sort(subPositions,[](const SubSeqment & s1, const SubSeqment & s2){
				return s1.pos_ > s2.pos_;
			});
			for(const auto & subPosition : subPositions){
				std::string subSegment = seq.seq_.substr(subPosition.pos_, subPosition.size_);
				if(printInfo){
					std::cout << "\t" <<"subsegment: " << subSegment << std::endl;
					std::cout << "\t" << "pos : " << subPosition.pos_ << std::endl;
					std::cout << "\t" << "size: " << subPosition.size_ << std::endl;

				}
				std::string head = "";
				std::string tail = "";
				if(subPosition.pos_ != 0){
					head = seq.seq_.substr(subPosition.pos_ - 1, pars.klen);
					if(printInfo){
						std::cout << "\t\t" << "head: " << std::endl;
						std::cout << "\t\t" << head << std::endl;
						std::cout << "\t\t" << compGraph.kCounts_[head] << std::endl;
					}
				}
				if(subPosition.pos_ + subPosition.size_ != len(seq)){
					tail = seq.seq_.substr(subPosition.pos_ + subPosition.size_ + 1 - pars.klen, pars.klen);
					if(printInfo){
						std::cout << "\t\t" << "tail: " << std::endl;
						std::cout << "\t\t" << tail << std::endl;
						std::cout << "\t\t" << compGraph.kCounts_[tail] << std::endl;
					}
				}
				std::unordered_map<std::string, uint32_t> headPositions;
				std::unordered_map<std::string, uint32_t> tailPositions;
				uint32_t replaceStart = subPosition.pos_;
				uint32_t replaceLen = subPosition.size_;
				if("" != head){
					--replaceStart;
					++replaceLen;
					for(const auto & tailEdge : compGraph.nodes_[compGraph.nodePositions_[head]]->tailEdges_){
						for(const auto & con : tailEdge->connectorInfos_){
							headPositions[con.readName_] = con.headPos_;
						}
					}
				}
				if("" != tail){
					++replaceLen;
					for(const auto & headEdge :compGraph.nodes_[compGraph.nodePositions_[tail]]->headEdges_){
						for(const auto & con : headEdge->connectorInfos_){
							tailPositions[con.readName_] = con.tailPos_;
						}
					}
				}
				if("" == head){
					for(const auto & tailPosition : tailPositions){
						headPositions[tailPosition.first] = 0;
					}
				}
				if("" == tail){
					for(const auto & headPosition : headPositions){
						tailPositions[headPosition.first] = len(seqs[seqnameToPosition[headPosition.first]]);
					}
				}
				std::unordered_map<std::string, uint32_t> subSeqCounts;
				for(const auto & headPosition : headPositions){
					if(headPosition.first == seq.name_){
						continue;
					}
					if(njh::in(headPosition.first, tailPositions) && tailPositions[headPosition.first] > headPosition.second){
						++subSeqCounts[seqs[seqnameToPosition[headPosition.first]].seq_.substr(headPosition.second, tailPositions[headPosition.first] - headPosition.second + pars.klen)];
					}
				}
				std::vector<std::string> subSeqsAbove;
				if(printInfo){
					std::cout << "\t" << "subSeqsCounts" << std::endl;
				}
				for(const auto & subSeq : subSeqCounts){
					if(printInfo){
						std::cout << "\t" << subSeq.first << ": " << subSeq.second << std::endl;
					}
					if(subSeq.second > std::max(pars.correctionOccurenceCutOff, pars.lowFreqCutOff)){
						subSeqsAbove.emplace_back(subSeq.first);
					}
				}
				if(1 == subSeqsAbove.size()){
					if(printInfo){
						std::cout << "replacing: " << replaceStart << "," << replaceLen << std::endl;
						std::cout << "old        : " << correctedSeqs[seqnameToPosition[seq.name_]].seq_.substr(replaceStart, replaceLen) << std::endl;
						std::cout << "replacement: " << subSeqsAbove.front() << std::endl;
					}

					correctedSeqs[seqnameToPosition[seq.name_]].seq_.replace(replaceStart, replaceLen, subSeqsAbove.front());
					correctedSeqs[seqnameToPosition[seq.name_]].qual_.erase(
							correctedSeqs[seqnameToPosition[seq.name_]].qual_.begin() + replaceStart,
							correctedSeqs[seqnameToPosition[seq.name_]].qual_.begin() + replaceStart + replaceLen);
					std::vector<uint32_t> insertQual(subSeqsAbove.front().size(), 40);
					correctedSeqs[seqnameToPosition[seq.name_]].qual_.insert(correctedSeqs[seqnameToPosition[seq.name_]].qual_.begin() + replaceStart,
							insertQual.begin(), insertQual.end());
				}
			}
		}
	}
	if(pars.debug){
		watch.logLapTimes(std::cout, true, 6, true);
	}
	return correctedSeqs;
}



uint32_t ContigsCompareGraphDev::findMinNonredundantKmer(
		uint32_t minimumKmerLen, const std::vector<seqInfo> &seqs,
		const std::string &funcName) {
	bool foundLength = false;
	uint32_t ret = minimumKmerLen;
	uint64_t maxLen = readVec::getMaxLength(seqs);

	while (ret < maxLen && !foundLength) {
		bool allPass = true;
		for (const auto &seq : seqs) {
			kmerInfo kinfo(seq.seq_, ret, false);
			for (const auto &k : kinfo.kmers_) {
				if (k.second.count_ > 1) {
					allPass = false;
					break;
				}
			}
		}
		if (allPass) {
			foundLength = true;
		} else {
			++ret;
		}
	}
	if (!foundLength) {
		std::stringstream ss;
		ss << funcName << ", error "
				<< " couldn't determine a minimum non-redundant kmer size" << "\n";
		throw std::runtime_error { ss.str() };
	}
	return ret;
}

std::vector<seqInfo> ContigsCompareGraphDev::filterSeqsOnLen(
		const std::vector<seqInfo> &seqs,
		const correctSeqsByGraphPars &graphCorrectingPars, const double lenFiltMultiplier,
		const SeqIOOptions &filteredSeqOpts) {
	auto seqLens = readVec::getLengths(seqs);
	uint32_t medianLen = std::round(vectorMedianRef(seqLens));
	uint32_t lenDiff = medianLen - std::round(medianLen * lenFiltMultiplier);
	std::vector<seqInfo> filterSeqs;
	std::vector<seqInfo> filteredOffSeqs;
	for (const auto &seq : seqs) {
		if (uAbsdiff(medianLen, len(seq)) > lenDiff) {
			filteredOffSeqs.emplace_back(seq);
		} else {
			filterSeqs.emplace_back(seq);
		}
	}
	if (!filteredOffSeqs.empty()
			&& filteredOffSeqs.size()
					<= graphCorrectingPars.correctionOccurenceCutOff) {
		SeqOutput::write(filteredOffSeqs, filteredSeqOpts);
	}
	return filterSeqs;
}

std::vector<seqInfo> ContigsCompareGraphDev::filterSeqsOnKmerSim(const std::vector<seqInfo> &seqs,
		const correctSeqsByGraphPars &graphCorrectingPars,
		const kmerInfo & refSeqKInfo,
		const double kSimCutOff,
		const SeqIOOptions & filteredSeqOpts){
	std::vector<seqInfo> filterSeqs;
	std::vector<seqInfo> filteredOffSeqs;
	for(const auto & seq : seqs){
		kmerInfo seqKInfo(seq.seq_, refSeqKInfo.kLen_, false);

		if(refSeqKInfo.compareKmers(seqKInfo).second < kSimCutOff){
			filteredOffSeqs.emplace_back(seq);
		}else{
			filterSeqs.emplace_back(seq);
		}
	}
	if(!filteredOffSeqs.empty() && filteredOffSeqs.size() <= graphCorrectingPars.correctionOccurenceCutOff){
		SeqOutput::write(filteredOffSeqs, filteredSeqOpts);
	}
	return filterSeqs;
}


std::vector<seqInfo> ContigsCompareGraphDev::readInSeqs(const SeqIOOptions & inOpts, seqInfo & refSeq, const std::string & refSeqName){
	std::vector<seqInfo> seqs;

	SeqInput reader(inOpts);
	reader.openIn();
	std::unordered_set<std::string> readNames;
	if("" != refSeq.name_){
		//supplying refname
		seqs.emplace_back(refSeq);
		readNames.emplace(refSeq.name_);
	}else if("" == refSeqName){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "if not supplying ref seq, nave to give a ref seq name to find within the input file: " << inOpts.firstName_ << "\n";
		throw std::runtime_error{ss.str()};
	}
	seqInfo seq;
	while(reader.readNextRead(seq)){
		bool matchingRefSeqName = false;
		if(seq.name_ == refSeqName){
			matchingRefSeqName = true;
		}else{
			if (MetaDataInName::nameHasMetaData(seq.name_)) {
				MetaDataInName meta(seq.name_);
				if (meta.containsMeta("sample")) {
					if(meta.getMeta("sample") == refSeqName){
						matchingRefSeqName = true;
					}
				}
			}
		}
		if(matchingRefSeqName){
			if("" != refSeq.name_){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "already found a seq maching reference name: " << refSeqName << "\n";
				ss << "previous: " << refSeq.name_ << "\n";
				ss << "current: " << seq.name_ << "\n";
				throw std::runtime_error{ss.str()};
			}else{
				refSeq = seq;
			}
		}
		if(njh::in(seq.name_, readNames)){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "already have seq name: " << seq.name_ << ", can't have repeat names"<< "\n";
			throw std::runtime_error{ss.str()};
		}
		seqs.emplace_back(seq);
		readNames.emplace(seq.name_);
	}
	if("" == refSeq.name_){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "didn't find a matching seq for reference: " << refSeqName<< "\n";
		throw std::runtime_error{ss.str()};
	}
	return seqs;
}


ContigsCompareGraphDev::processConservedNodesVecRes ContigsCompareGraphDev::processConservedNodesVec(const std::vector<std::shared_ptr<node>> & nodes){
	processConservedNodesVecRes ret;
	if(1 == nodes.size()){
		for(const auto & name : nodes.front()->inReadNamesIdx_){
			ret.nameToSubSegments[name].emplace_back(nodes.front()->k_, 0);
		}
	}else if(nodes.size() > 1){
		for(const auto & n : nodes){
			std::unordered_map<std::string, uint32_t> nodePositionPerSeq;
			for(const auto & h : n->headEdges_){
				for(const auto & con : h->connectorInfos_){
					nodePositionPerSeq[con.readName_] = con.tailPos_;
				}
			}
			for(const auto & t : n->tailEdges_){
				for(const auto & con : t->connectorInfos_){
					nodePositionPerSeq[con.readName_] = con.headPos_ - (n->k_.size() - n->kLen_);
				}
			}
			for(const auto & seq : nodePositionPerSeq){
				ret.nameToSubSegments[seq.first].emplace_back(n->k_, seq.second);
			}
		}
	}else{
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "nodes vector is empty, size:" << nodes.size()<< "\n";
		throw std::runtime_error{ss.str()};
	}

	std::unordered_set<std::string> paths;
	for(auto & subSeqs : ret.nameToSubSegments){
		njh::sort(subSeqs.second, [](const ContigsCompareGraphDev::SubSeqPos & p1, const ContigsCompareGraphDev::SubSeqPos & p2 ){
			return p1.pos_ < p2.pos_;
		});
		std::string path = "";
		for(const auto & sub : subSeqs.second){
			if("" != path){
				path += "-";
			}
			path += sub.subSeq_;
		}
		paths.emplace(path);
	}
	if(paths.size() > 1){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << " error more than one path , " << njh::conToStr(paths, ",") << "\n";
		throw std::runtime_error{ss.str()};
	}

	std::unordered_map<std::string, std::vector<uint32_t>> startingPositions;
	for(auto & subSeqs : ret.nameToSubSegments){
		for(const auto & sub : subSeqs.second){
			startingPositions[sub.subSeq_].emplace_back(sub.pos_);
		}
	}
	std::unordered_map<std::string, std::vector<uint32_t>> startingPositionsAvg;
	for(const auto & positions : startingPositions){
		startingPositionsAvg[positions.first].emplace_back(vectorMean(positions.second));
	}
	VecStr conservedSeqs;
	for(const auto & subNode : nodes){
		conservedSeqs.emplace_back(subNode->k_);
	}
	njh::sort(conservedSeqs, [&startingPositionsAvg](const std::string & seq1, const std::string & seq2){
		return startingPositionsAvg[seq1] < startingPositionsAvg[seq2];
	});
	for(auto & subSeqs : ret.nameToSubSegments){
		for(const auto & sub : subSeqs.second){
			auto outRegion = sub.genBed6Region(subSeqs.first);
			ret.nameToSubSegPositions[subSeqs.first].emplace_back(outRegion);
			ret.subSeqToNameToPos[sub.subSeq_].emplace(subSeqs.first, outRegion);
			//std::cout << "\t" << sub.pos_ << ": " << sub.subSeq_ << std::endl;
			//sharedLocsOut << outRegion.toDelimStr()<< std::endl;
		}
	}

	for(auto & subPositions : ret.nameToSubSegPositions){
		BedUtility::coordSort(subPositions.second, false);
	}


	uint32_t subseqCount = 0;

	for(const auto & subSeq : conservedSeqs){
		std::string ID = njh::pasteAsStr("subseq.", njh::leftPadNumStr<uint32_t>(subseqCount, nodes.size()));
		ret.subseqToIDKey[subSeq] = ID;
		ret.IDtoSubseqKey[ID] = subSeq;
		++subseqCount;
	}

	return ret;
}



} //namespace njhseq

