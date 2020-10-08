/*
 * ContigsCompareGraphDev.cpp
 *
 *  Created on: Nov 5, 2019
 *      Author: nicholashathaway
 */


#include "ContigsCompareGraph.hpp"


namespace njhseq {

ContigsCompareGraphDev::node::node(const std::string & k, uint32_t cnt, uint32_t kLen) :
		k_(k),uid_(k), cnt_(cnt), kLen_(kLen) {
}
void ContigsCompareGraphDev::node::resetVisitCount(){
	visitCount_ = 0;
}
void ContigsCompareGraphDev::node::addHead(const std::shared_ptr<edge> & e){
	headEdges_.push_back(e);
}
void ContigsCompareGraphDev::node::addTail(const std::shared_ptr<edge> & e){
	tailEdges_.push_back(e);
}
bool ContigsCompareGraphDev::node::headless() const {
	return 0 == headCount();
}
uint32_t ContigsCompareGraphDev::node::headCount() const {
	uint32_t headCount = 0;
	for (const auto & head : headEdges_) {
		if (head->on_) {
			++headCount;
		}
	}
	return headCount;
}
std::shared_ptr<ContigsCompareGraphDev::edge> ContigsCompareGraphDev::node::getFirstOnHeadEdge() const {
	for (const auto & head : headEdges_) {
		if (head->on_) {
			return head;
		}
	}
	return nullptr;
}
std::shared_ptr<ContigsCompareGraphDev::edge> ContigsCompareGraphDev::node::getLastOnHeadEdge() const {
	for (const auto & head : iter::reversed(headEdges_)) {
		if (head->on_) {
			return head;
		}
	}
	return nullptr;
}
bool ContigsCompareGraphDev::node::tailless() const{
	return 0 == tailCount();
}
uint32_t ContigsCompareGraphDev::node::tailCount() const{
	uint32_t tailCount = 0;
	for(const auto & tail : tailEdges_){
		if(tail->on_){
			++tailCount;
		}
	}
	return tailCount;
}
std::shared_ptr<ContigsCompareGraphDev::edge> ContigsCompareGraphDev::node::getFirstOnTailEdge() const {
	for (const auto & tail : tailEdges_) {
		if (tail->on_) {
			return tail;
		}
	}
	return nullptr;
}
std::shared_ptr<ContigsCompareGraphDev::edge> ContigsCompareGraphDev::node::getLastOnTailEdge() const {
	for (const auto & tail : iter::reversed(tailEdges_)) {
		if (tail->on_) {
			return tail;
		}
	}
	return nullptr;
}
std::shared_ptr<ContigsCompareGraphDev::node> ContigsCompareGraphDev::node::getHeadNode(const std::string & uid)const{
	std::shared_ptr<node> ret = nullptr;
	for(const auto & head : headEdges_){
		if(uid == head->head_.lock()->uid_){
			ret = head->head_.lock();
			break;
		}
	}
	return ret;
}
std::shared_ptr<ContigsCompareGraphDev::node> ContigsCompareGraphDev::node::getTailNode(const std::string & uid)const{
	std::shared_ptr<node> ret = nullptr;
	for(const auto & tail : tailEdges_){
		if(uid == tail->tail_.lock()->uid_){
			ret = tail->tail_.lock();
			break;
		}
	}
	return ret;
}


} // namespace njhseq

