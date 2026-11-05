#pragma once

// Created by Nicholas Hathaway on 7/12/23.
/* 
    
*/

#include <njhcpp/common.h>
#include <njhcpp/utils.h>
#include <njhseq/common/typedefs.hpp>
#include <njhcpp/simulation/randomGenerator.hpp>

namespace njhseq {
/**
 * @brief Undirected weighted graph, each node is a different region, each edge is a different set of primer pairs combo with a cost of having those primers together
 */
class PossibleAmpliconPanelGraph{

public:

	class edge;

	class node{
	public:
		explicit node(std::string region): region_(std::move(region)){

		}
		std::string region_;
		bool visited_{false};

		std::vector<std::shared_ptr<edge>> edges_;

		void turnOnAllEdges() {
			for (auto &e: edges_) {
				e->on_ = true;
			}
		}

		void reset() {
			visited_ = false;
			turnOnAllEdges();
		}

		void sortEdgesMaxTotalWeightTop(){
			njh::sort(edges_, [](const std::shared_ptr<edge> & e1, const std::shared_ptr<edge> & e2){
				return e1->getTotalWeights() > e2->getTotalWeights();
			});
		}

		void sortEdgesMaxRegionPrimerWeightTop(){
			njh::sort(edges_, [this](const std::shared_ptr<edge> & e1, const std::shared_ptr<edge> & e2){
				if(e1->primerPairIndvWeights_[e1->regionToPrimerName_[this->region_]] == e2->primerPairIndvWeights_[e2->regionToPrimerName_[this->region_]]){
					return e1->getTotalWeights() > e2->getTotalWeights();
				}else{
					return e1->primerPairIndvWeights_[e1->regionToPrimerName_[this->region_]] > e2->primerPairIndvWeights_[e2->regionToPrimerName_[this->region_]];
				}
			});
		}

		std::string getTopOnPrimerMaxTotalWeightOnTop(){
			sortEdgesMaxTotalWeightTop();
			if(edges_.empty()){
				return {""};
			}
			return edges_.front()->regionToPrimerName_[region_];
		}

		void toggleOnEdgesForPrimer(const std::string & name){
			for(const auto & e : edges_){
				e->on_ = e->regionToPrimerName_[region_] == name;
			}
		}

		void turnOffAllOtherPrimers(const std::string & name){
			for(const auto & e : edges_){
				if(e->regionToPrimerName_[region_] != name){
					e->on_ = false;
				}
			}
		}

		void turnOffPrimer(const std::string & name){
			for(const auto & e : edges_){
				if(e->regionToPrimerName_[region_] == name){
					e->on_ = false;
				}
			}
		}

		void turnOffEdge(const std::string & p1, const std::string & p2){
			for(const auto & e : edges_){
				auto primerPairs = njh::getVecOfMapValues(e->regionToPrimerName_);
				if(njh::in(p1, primerPairs) && njh::in(p2, primerPairs)){
					e->on_ = false;
				}
			}
		}

		/**
		 * @brief turn off all edges to the other region in the keepEdge expect for the keepEdge
		 * @param keepEdge the edge to keep, turn off all other edges going to this region
		 */
		void turnOfOtherEdgesForOtherRegion(const edge &keepEdge) {
			for (auto &e: edges_) {
				if (njh::in(keepEdge.regionToOtherRegion_.at(region_), e->regionToPrimerName_) &&
						e->regionToPrimerName_[keepEdge.regionToOtherRegion_.at(region_)] !=
						keepEdge.regionToPrimerName_.at(keepEdge.regionToOtherRegion_.at(region_))) {
					e->on_ = false;
				}
			}
		}

		[[nodiscard]] VecStr getUniquePrimersOn() const {
			auto onCounts = getTotalOnEdgesForPrimers();
			return njh::getVecOfMapKeys(onCounts);
		}

		[[nodiscard]] uint32_t getTotalUniquePrimersOn() const {
			auto onCounts = getTotalOnEdgesForPrimers();
			return onCounts.size();
		}

		[[nodiscard]] std::unordered_map<std::string, uint32_t> getTotalOnEdgesForPrimers() const {
			std::unordered_map<std::string, uint32_t> ret;
			for (const auto &e: edges_) {
				if (e->on_) {
					++ret[e->regionToPrimerName_[region_]];
				}
			}
			return ret;
		}

		[[nodiscard]] std::unordered_map<std::string, uint32_t> getTotalOnEdgesForOtherRegions() const {
			std::unordered_map<std::string, uint32_t> ret;
			for (const auto &e: edges_) {
				if (e->on_) {
					++ret[e->regionToOtherRegion_[region_]];
				}
			}
			return ret;
		}


	};

	class edge {
	public:
		std::unordered_map<std::string, std::string> regionToPrimerName_;/*< key is region name, key is primer pair name*/
		std::unordered_map<std::string, std::string> regionToOtherRegion_;/*< key is region name, key is the other region in this edge*/

		bool on_{true};
		double primerPairVsPrimerPairWeight_{0};

		double primerPairVsPrimerPairUnSpecAmpCnt_{0};
		/*< the same count from the sensitive hot-3'-end pass, where both ends of the product are hot */
		double primerPairVsPrimerPairHotUnSpecAmpCnt_{0};

		std::unordered_map<std::string, double> primerPairIndvWeights_;

		[[nodiscard]] std::string getUid() const{
			auto pairs = njh::getVecOfMapValues(regionToPrimerName_);
			njh::sort(pairs);
			return njh::pasteAsStr(pairs);
		}

		[[nodiscard]] double sumInvWeights() const {
			double ret = 0;
			for (const auto &pweight: primerPairIndvWeights_) {
				ret += pweight.second;
			}
			return ret;
		}

		[[nodiscard]] double getTotalWeights() const {
			return primerPairVsPrimerPairWeight_ + sumInvWeights();
		}
	};

	explicit PossibleAmpliconPanelGraph(const VecStr & regions){
		for(const auto & r : regions){
			nodes_.emplace_back(r);
		}
		setNodeIdx();
	}
	void setNodeIdx(){
		nodeIdx_.clear();
		for(const auto & n : iter::enumerate(nodes_)){
			nodeIdx_[n.element.region_] = n.index;
		}
	}

	void resetAll(){
		njh::for_each(nodes_,[](auto & n){
			n.reset();
		});
	}

	void turnOffPrimer(const std::string & name){
		for(auto & n : nodes_){
			n.turnOffPrimer(name);
		}
	}

	void resetAndGenerateRandomPool(){
		resetAll();
		njh::randomGenerator rgen;
		for(auto & n : nodes_){
			auto randomPrimer = rgen.unifRandSelection(n.getUniquePrimersOn());
			n.turnOffAllOtherPrimers(randomPrimer);
		}
	}

	std::vector<std::shared_ptr<PossibleAmpliconPanelGraph::edge>> greedyDetermineHeaviestPool(){
		auto keepGoingFunc = [this](){
			bool keepGoing = false;
			for(const auto & n : nodes_){
				if(n.getTotalUniquePrimersOn() > 1){
					keepGoing = true;
					break;
				}
			}
			return keepGoing;
		};
		std::vector<std::shared_ptr<PossibleAmpliconPanelGraph::edge>> topEdges;
		std::unordered_set<std::string> topEdgeUids;
		std::shared_ptr<PossibleAmpliconPanelGraph::edge> topEdge = nullptr;
		while(keepGoingFunc()){
			topEdge = nullptr;
			for(auto & n : nodes_){
				for(const auto & e : n.edges_){
					if(e->on_ && !njh::in(e->getUid(), topEdgeUids)){
						if(nullptr == topEdge || e->getTotalWeights() > topEdge->getTotalWeights()){
							topEdge = e;
						}
					}
				}
			}
			if(nullptr == topEdge){
				break;
			}
			//turn off all other primers
			nodes_[nodeIdx_[topEdge->regionToOtherRegion_.begin()->first]].turnOffAllOtherPrimers(topEdge->regionToPrimerName_[topEdge->regionToOtherRegion_.begin()->first]);
			nodes_[nodeIdx_[topEdge->regionToOtherRegion_[topEdge->regionToOtherRegion_.begin()->first]]].turnOffAllOtherPrimers(topEdge->regionToPrimerName_[topEdge->regionToOtherRegion_[topEdge->regionToOtherRegion_.begin()->first]]);
			//turn off all other connections between these two regions
			nodes_[nodeIdx_[topEdge->regionToOtherRegion_.begin()->first]].turnOfOtherEdgesForOtherRegion(*topEdge);
			nodes_[nodeIdx_[topEdge->regionToOtherRegion_[topEdge->regionToOtherRegion_.begin()->first]]].turnOfOtherEdgesForOtherRegion(*topEdge);
			topEdges.emplace_back(topEdge);
			topEdgeUids.emplace(topEdge->getUid());
		}
		return topEdges;
	}

	/**
	 * @brief Start from the heaviest edge, then keep adding the primer pair that scores best against every primer pair chosen so far
	 *
	 * A primer pair with an edge turned off to any chosen primer pair can't be added, and a region left with none ends up without a primer pair
	 */
	void greedyAccumulateHeaviestPool(){
		std::shared_ptr<edge> startEdge = nullptr;
		for(const auto & e : edges_){
			if(e->on_ && (nullptr == startEdge || e->getTotalWeights() > startEdge->getTotalWeights())){
				startEdge = e;
			}
		}
		std::unordered_map<std::string, std::string> chosenPrimerPairs;
		if(nullptr != startEdge){
			chosenPrimerPairs = startEdge->regionToPrimerName_;
		}
		while(nullptr != startEdge && chosenPrimerPairs.size() < nodes_.size()){
			bool found = false;
			double bestScore = 0;
			std::string bestRegion;
			std::string bestPrimerPair;
			for(const auto & n : nodes_){
				if(njh::in(n.region_, chosenPrimerPairs)){
					continue;
				}
				VecStr candidates;
				std::unordered_map<std::string, double> vsChosenWeights;
				std::unordered_map<std::string, uint32_t> vsChosenCounts;
				for(const auto & e : n.edges_){
					const auto & otherRegion = e->regionToOtherRegion_.at(n.region_);
					auto chosen = chosenPrimerPairs.find(otherRegion);
					if(e->on_ && chosenPrimerPairs.end() != chosen && e->regionToPrimerName_.at(otherRegion) == chosen->second){
						const auto & candidate = e->regionToPrimerName_.at(n.region_);
						if(!njh::in(candidate, vsChosenCounts)){
							candidates.emplace_back(candidate);
							vsChosenWeights[candidate] = e->primerPairIndvWeights_.at(candidate);
						}
						vsChosenWeights[candidate] += e->primerPairVsPrimerPairWeight_;
						++vsChosenCounts[candidate];
					}
				}
				for(const auto & candidate : candidates){
					if(vsChosenCounts[candidate] == chosenPrimerPairs.size() && (!found || vsChosenWeights[candidate] > bestScore)){
						found = true;
						bestScore = vsChosenWeights[candidate];
						bestRegion = n.region_;
						bestPrimerPair = candidate;
					}
				}
			}
			if(!found){
				break;
			}
			chosenPrimerPairs[bestRegion] = bestPrimerPair;
		}
		for(auto & n : nodes_){
			auto chosen = chosenPrimerPairs.find(n.region_);
			if(chosenPrimerPairs.end() != chosen){
				n.turnOffAllOtherPrimers(chosen->second);
			}else{
				for(auto & e : n.edges_){
					e->on_ = false;
				}
			}
		}
	}

	double getTotalOnWeight(){
		double ret = 0;
		for(const auto & n : nodes_){
			for(const auto & e : n.edges_){
				if(e->on_){
					ret += e->getTotalWeights();
				}
			}
		}
		return ret;
	}

	[[nodiscard]] std::set<std::string> getCurrentOnPrimerPairs() const {
		std::set < std::string > ret;
		for (const auto &n: nodes_) {
			for (const auto &e: n.edges_) {
				if (e->on_) {
					njh::addVecToSet(njh::getVecOfMapValues(e->regionToPrimerName_), ret);
				}
			}
		}
		return ret;
	}

	std::vector<std::shared_ptr<edge>> edges_;
	std::vector<node> nodes_;
	std::unordered_map<std::string, uint32_t> nodeIdx_;
};

}  // namespace njhseq

