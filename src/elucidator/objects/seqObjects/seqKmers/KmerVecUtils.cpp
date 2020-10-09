/*
 * KmerVecUtils.cpp
 *
 *  Created on: May 24, 2016
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
#include "KmerVecUtils.hpp"



namespace njhseq {



std::vector<std::vector<double>> getKmerAccerDistance(
		std::vector<std::unique_ptr<seqWithKmerInfo>>& reads, uint32_t kmerStart,
		uint32_t kmerStop, uint32_t numThreads, bool useKNumber, bool verbose) {
	//kmerStop is inclusive
	if(kmerStart > kmerStop){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": kmerStop is less than kmerStart, kmerStart: " << kmerStart <<
				", KmerStop: " << kmerStop;
		throw std::runtime_error{ss.str()};
	}
	std::vector<std::vector<double>> distances;
	if (useKNumber) {
		std::function<
				uint32_t(const std::unique_ptr<seqWithKmerInfo> &,
						const std::unique_ptr<seqWithKmerInfo> &)> disFun =
				[](const std::unique_ptr<seqWithKmerInfo> & read1,
						const std::unique_ptr<seqWithKmerInfo> & read2) {
					auto dist = read1->compareKmers(*read2);
					return dist.first;
				};
		std::unordered_map<uint32_t, std::vector<std::vector<uint32_t>>>distanceMaps;
		for (uint32_t k = kmerStart; k < kmerStop + 1; ++k) {
			if(verbose){
				std::cout << "K: " << k << std::endl;
				std::cout << "\tIndexing Kmers" << std::endl;
			}
			njh::stopWatch watch;
			allSetKmers(reads, k, false);

			if(verbose){
				std::cout << "\tIndexing Time: " << watch.totalTimeFormatted(0) << std::endl;
				std::cout << "\tCalculating Distances" << std::endl;
			}
			watch.reset();
			distanceMaps[k] = getDistance(reads, numThreads, disFun);
			if(verbose){
				std::cout << "\tCalculating Distances Time: " << watch.totalTimeFormatted(0) << std::endl;
			}
		}
		for (const auto rowPos : iter::range(distanceMaps[kmerStart].size())) {
			std::vector<double> temp;
			for (uint32_t i = 0; i < distanceMaps[kmerStart][rowPos].size(); ++i) {
				temp.emplace_back(0.00);
			}
			distances.emplace_back(temp);
		}
		for (const auto rowPos : iter::range(distances.size())) {
			for (const auto colPos : iter::range(distances[rowPos].size())) {
				std::vector<double> differences;
				for (uint32_t k = kmerStart; k < kmerStop; ++k) {
					differences.emplace_back(
							uAbsdiff(distanceMaps[k][rowPos][colPos],
									distanceMaps[k + 1][rowPos][colPos]));
				}
				distances[rowPos][colPos] = vectorMean(differences);
			}
		}
	} else {
		std::function<
				double(const std::unique_ptr<seqWithKmerInfo> &,
						const std::unique_ptr<seqWithKmerInfo> &)> disFun =
				[](const std::unique_ptr<seqWithKmerInfo> & read1,
						const std::unique_ptr<seqWithKmerInfo> & read2) {
					auto dist = read1->compareKmers(*read2);
					return dist.second;
				};
		std::unordered_map<uint32_t, std::vector<std::vector<double>>>distanceMaps;
		for (uint32_t k = kmerStart; k < kmerStop + 1; ++k) {
			if(verbose){
				std::cout << "K: " << k << std::endl;
				std::cout << "\tIndexing Kmers" << std::endl;
			}
			njh::stopWatch watch;
			allSetKmers(reads, k, false);

			if(verbose){
				std::cout << "\tIndexing Time: " << watch.totalTimeFormatted(0) << std::endl;
				std::cout << "\tCalculating Distances" << std::endl;
			}
			watch.reset();
			distanceMaps[k] = getDistance(reads, numThreads, disFun);
			if(verbose){
				std::cout << "\tCalculating Distances Time: " << watch.totalTimeFormatted(0) << std::endl;
			}
		}
		distances = distanceMaps[kmerStart];
		for (const auto rowPos : iter::range(distances.size())) {
			for (const auto colPos : iter::range(distances[rowPos].size())) {
				std::vector<double> differences;
				for (uint32_t k = kmerStart; k < kmerStop; ++k) {
					differences.emplace_back(
							std::abs(
									distanceMaps[k][rowPos][colPos]
											- distanceMaps[k + 1][rowPos][colPos]));
				}
				distances[rowPos][colPos] = vectorMean(differences);
			}
		}
	}
	return distances;
}


readDistGraph<double> genKmerAccerDistGraphWithDbSmartBuild(
		std::vector<std::unique_ptr<seqWithKmerInfo>> & reads,
		uint32_t kmerStop, const readDistGraph<double>::dbscanPars & pars,
		bool breakLargeIndelCons, uint32_t largeIndel, aligner & alignerObj,
		bool doTies, uint32_t numThreads, bool verbose) {
	if (verbose) {
		std::cout << njh::bashCT::bold << "Computing kmer distances and Building Graph"
				<< njh::bashCT::reset << std::endl;
	}
	if(kmerStop < 3 ){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": kmerStop is less than 3, should at least be 3 or greater " <<
				", KmerStop: " << kmerStop;
		throw std::runtime_error{ss.str()};
	}

	std::function<
			double(const std::unique_ptr<seqWithKmerInfo> &,
					const std::unique_ptr<seqWithKmerInfo> &)> disFun =
			[](const std::unique_ptr<seqWithKmerInfo> & read1,
					const std::unique_ptr<seqWithKmerInfo> & read2) {
				auto dist = read1->compareKmers(*read2);
				return dist.second;
			};
	std::unordered_map<uint32_t, std::vector<std::vector<double>>>distanceMaps;
	for (uint32_t k = 2; k < 4; ++k) {
		if(verbose){
			std::cout << "K: " << k << std::endl;
			std::cout << "\tIndexing Kmers" << std::endl;
		}
		njh::stopWatch watch;
		allSetKmers(reads, k, false);

		if(verbose){
			std::cout << "\tIndexing Time: " << watch.totalTimeFormatted(0) << std::endl;
			std::cout << "\tCalculating Distances" << std::endl;
		}
		watch.reset();
		distanceMaps[k] = getDistance(reads, numThreads, disFun);
		if(verbose){
			std::cout << "\tCalculating Distances Time: " << watch.totalTimeFormatted(0) << std::endl;
		}
	}
	std::vector<std::vector<double>> distances = distanceMaps[2];
	for (const auto rowPos : iter::range(distances.size())) {
		for (const auto colPos : iter::range(distances[rowPos].size())) {
			/*double res = std::abs(
					distanceMaps[2][rowPos][colPos]
							- distanceMaps[3][rowPos][colPos] );*/
			double res = distanceMaps[2][rowPos][colPos]
					- distanceMaps[3][rowPos][colPos];
			if(res <= pars.eps_){
				distances[rowPos][colPos] = res;
			}else{
				distances[rowPos][colPos] = std::numeric_limits<double>::max();
			}
		}
	}
	for(const auto k : iter::range<uint32_t>(4,kmerStop + 1)){
		if(verbose){
			std::cout << "K: " << k << std::endl;
			std::cout << "\tIndexing Kmers" << std::endl;
		}
		njh::stopWatch watch;
		allSetKmers(reads, k, false);

		if(verbose){
			std::cout << "\tIndexing Time: " << watch.totalTimeFormatted(0) << std::endl;
			std::cout << "\tCalculating Distances" << std::endl;
		}
		watch.reset();
		if(k > 5){
			for (const auto rowPos : iter::range(distances.size())) {
				for (const auto colPos : iter::range(distances[rowPos].size())) {
					std::vector<double> differences;
					for (uint32_t k = 2; k < kmerStop; ++k) {
						/*differences.emplace_back(
						 std::abs(
						 distanceMaps[k][rowPos][colPos]
						 - distanceMaps[k + 1][rowPos][colPos]));*/
						differences.emplace_back(
								distanceMaps[k][rowPos][colPos]
										- distanceMaps[k + 1][rowPos][colPos]);
					}
					if (std::numeric_limits<double>::max() != distances[rowPos][colPos]) {
						double res = vectorMean(differences);
						if(res <= pars.eps_){
							distances[rowPos][colPos] = res;
						}else{
							distances[rowPos][colPos] = std::numeric_limits<double>::max();
						}
					}
				}
			}
		}
		std::vector<std::vector<double>> currentKDists;
		std::vector<std::pair<uint32_t, uint32_t>> indices;
	  for(const auto pos : iter::range(reads.size())){
	  	currentKDists.emplace_back(std::vector<double>(pos));
	  	for(const auto secondPos : iter::range(pos)){
	  		if(std::numeric_limits<double>::max() != distances[pos][secondPos]){
	  			indices.emplace_back(pos, secondPos);
	  		}
	  	}
	  }
	  if(numThreads < 2 || numThreads >= reads.size()){
	  	paritialDis(reads, indices, currentKDists, disFun);
	  }else{
	  	std::vector<std::thread> threads;
	  	uint32_t step = std::round(indices.size()/static_cast<double>(numThreads));
	  	std::vector<std::vector<std::pair<uint32_t, uint32_t>>> indsSplit;
	  	for(const auto tNum : iter::range(numThreads - 1)){
	  		std::vector<std::pair<uint32_t, uint32_t>> temp {indices.begin() + tNum * step,
	  			indices.begin() + (tNum + 1)*step};
	  		indsSplit.emplace_back(temp);
	  	}
	  	std::vector<std::pair<uint32_t, uint32_t>> temp {indices.begin() + (numThreads - 1) * step,
	  	  			indices.end()};
	  	indsSplit.emplace_back(temp);
	  	for(const auto tNum : iter::range(numThreads)){
	  		threads.push_back(std::thread(paritialDis<std::unique_ptr<seqWithKmerInfo>,double>, std::cref(reads),
	    			indsSplit[tNum], std::ref(currentKDists), disFun));
	  	}
	  	for(auto & t : threads){
	  		t.join();
	  	}
	  }
		distanceMaps[k] = currentKDists;
		if(verbose){
			std::cout << "\tCalculating Distances Time: " << watch.totalTimeFormatted(0) << std::endl;
		}
	}



	readDistGraph<double> distanceGraph(reads);
	for (const auto rowPos : iter::range(distances.size())) {
		for (const auto colPos : iter::range(distances[rowPos].size())) {
			std::vector<double> differences;
			for (uint32_t k = 2; k < kmerStop; ++k) {
				/*differences.emplace_back(
				 std::abs(
				 distanceMaps[k][rowPos][colPos]
				 - distanceMaps[k + 1][rowPos][colPos]));*/
				differences.emplace_back(
						distanceMaps[k][rowPos][colPos]
								- distanceMaps[k + 1][rowPos][colPos]);
			}
			if (std::numeric_limits<double>::max() != distances[rowPos][colPos]) {
				double meanDecrease = vectorMean(differences);
				if (meanDecrease <= pars.eps_) {
					distanceGraph.addEdge(reads[rowPos]->seqBase_.name_,
							reads[colPos]->seqBase_.name_, meanDecrease);
				}
			}
		}
	}
	//distanceGraph.turnOffEdgesAbove(pars.eps_);
	//turn off connections if they have large indels
	if (breakLargeIndelCons) {
		if (verbose) {
			std::cout << njh::bashCT::bold << "Breaking connections with large indels"
					<< njh::bashCT::reset << std::endl;
		}
		for (auto & e : distanceGraph.edges_) {
			if (e->on_) {
				auto seq1 = e->nodeToNode_.begin()->second.lock()->value_;
				auto seq2 =
						distanceGraph.nodes_[distanceGraph.nameToNodePos_[e->nodeToNode_.begin()->first]]->value_;
				alignerObj.alignCacheGlobal(seq1, seq2);
				alignerObj.profilePrimerAlignment(seq1, seq2);
				bool foundLargeIndel = false;
				for (const auto & g : alignerObj.comp_.distances_.alignmentGaps_) {
					if (g.second.size_ >= largeIndel) {
						foundLargeIndel = true;
						break;
					}
				}
				if (foundLargeIndel) {
					e->on_ = false;
				}
			}
		}
	}

	distanceGraph.resetBestAndVistEdges();
	distanceGraph.resetVisitedNodes();
	distanceGraph.allDetermineLowestBest(doTies);
	distanceGraph.removeOffEdges();
	distanceGraph.dbscan(pars);
	distanceGraph.assignNoiseNodesAGroup();
	return distanceGraph;
}

readDistGraph<double> genKmerAccerDistGraph(
		std::vector<std::unique_ptr<seqWithKmerInfo>> & reads, uint32_t kmerStart,
		uint32_t kmerStop, double avgKmerSimDecreasingRateCutOff,
		bool breakLargeIndelCons, uint32_t largeIndel, aligner & alignerObj,
		bool doTies, uint32_t numThreads, bool verbose) {
	if (verbose) {
		std::cout << njh::bashCT::bold << "Computing kmer distances and Building Graph"
				<< njh::bashCT::reset << std::endl;
	}
	bool useKmerNum = false;
	readDistGraph<double> distanceGraph(getKmerAccerDistance(reads,
			kmerStart, kmerStop, numThreads, useKmerNum, verbose), reads);
	distanceGraph.turnOffEdgesAbove(avgKmerSimDecreasingRateCutOff);
	//turn off connections if they have large indels
	if (breakLargeIndelCons) {
		if (verbose) {
			std::cout << njh::bashCT::bold << "Breaking connections with large indels"
					<< njh::bashCT::reset << std::endl;
		}
		for (auto & e : distanceGraph.edges_) {
			if (e->on_) {
				auto seq1 = e->nodeToNode_.begin()->second.lock()->value_;
				auto seq2 =
						distanceGraph.nodes_[distanceGraph.nameToNodePos_[e->nodeToNode_.begin()->first]]->value_;
				alignerObj.alignCacheGlobal(seq1, seq2);
				alignerObj.profilePrimerAlignment(seq1, seq2);
				bool foundLargeIndel = false;
				for (const auto & g : alignerObj.comp_.distances_.alignmentGaps_) {
					if (g.second.size_ >= largeIndel) {
						foundLargeIndel = true;
						break;
					}
				}
				if (foundLargeIndel) {
					e->on_ = false;
				}
			}
		}
	}
	distanceGraph.resetBestAndVistEdges();
	distanceGraph.resetVisitedNodes();
	distanceGraph.allDetermineLowestBest(doTies);
	distanceGraph.determineGroups();
	return distanceGraph;
}

readDistGraph<double> genKmerAccerDistGraphWithDb(
		std::vector<std::unique_ptr<seqWithKmerInfo>> & reads, uint32_t kmerStart,
		uint32_t kmerStop, const readDistGraph<double>::dbscanPars & pars,
		bool breakLargeIndelCons, uint32_t largeIndel, aligner & alignerObj,
		bool doTies, uint32_t numThreads, bool verbose) {
	if (verbose) {
		std::cout << njh::bashCT::bold << "Computing kmer distances and Building Graph"
				<< njh::bashCT::reset << std::endl;
	}
	bool useKmerNum = false;
	auto kDists = getKmerAccerDistance(reads,
				kmerStart, kmerStop, numThreads, useKmerNum, verbose);
	readDistGraph<double> distanceGraph(reads);
  for(const auto pos : iter::range(kDists.size())){
  	for(const auto subPos : iter::range<uint64_t>(kDists[pos].size())){
			if (kDists[pos][subPos] <= pars.eps_) {
				distanceGraph.addEdge(reads[pos]->seqBase_.name_,
						reads[subPos]->seqBase_.name_, kDists[pos][subPos]);
			}
  	}
  }
	distanceGraph.turnOffEdgesAbove(pars.eps_);
	//turn off connections if they have large indels
	if (breakLargeIndelCons) {
		if (verbose) {
			std::cout << njh::bashCT::bold << "Breaking connections with large indels"
					<< njh::bashCT::reset << std::endl;
		}
		for (auto & e : distanceGraph.edges_) {
			if (e->on_) {
				auto seq1 = e->nodeToNode_.begin()->second.lock()->value_;
				auto seq2 =
						distanceGraph.nodes_[distanceGraph.nameToNodePos_[e->nodeToNode_.begin()->first]]->value_;
				alignerObj.alignCacheGlobal(seq1, seq2);
				alignerObj.profilePrimerAlignment(seq1, seq2);
				bool foundLargeIndel = false;
				for (const auto & g : alignerObj.comp_.distances_.alignmentGaps_) {
					if (g.second.size_ >= largeIndel) {
						foundLargeIndel = true;
						break;
					}
				}
				if (foundLargeIndel) {
					e->on_ = false;
				}
			}
		}
	}

	distanceGraph.resetBestAndVistEdges();
	distanceGraph.resetVisitedNodes();
	distanceGraph.allDetermineLowestBest(doTies);
	distanceGraph.removeOffEdges();
	distanceGraph.dbscan(pars);
	distanceGraph.assignNoiseNodesAGroup();
	return distanceGraph;
}


}  // namespace njhseq
