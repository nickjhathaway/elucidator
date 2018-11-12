#pragma once
//
//  debGraph.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 12/11/13.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
//
//
// njhseq - A library for analyzing sequence data
// Copyright (C) 2012-2016 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
// Jeffrey Bailey <Jeffrey.Bailey@umassmed.edu>
//
// This file is part of njhseq.
//
// njhseq is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// njhseq is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with njhseq.  If not, see <http://www.gnu.org/licenses/>.
//
#include "elucidator/common.h"
#include "elucidator/objects/dataContainers/graphs/nhGraph.hpp"
#include "elucidator/simulation/simulationCommon.hpp"

namespace njhseq {

class debGraph {

 public:
  // struc for edge
  struct edge {
    edge(const readObject& read, int minOverLap, int maxOverLap) : _read(read) {
      for (auto i : iter::range(minOverLap, maxOverLap + 1)) {
        _prefixes[i] = _read.seqBase_.seq_.substr(0, i);
        _suffixes[i] =
            _read.seqBase_.seq_.substr(read.seqBase_.seq_.length() - i);
      }
    }
    edge(const readObject& read,
         const std::unordered_map<int, std::string>& prefixes,
         const std::unordered_map<int, std::string>& suffixes)
        : _read(read), _prefixes(prefixes), _suffixes(suffixes) {}
    edge(const readObject& read) : _read(read) {}
    readObject _read;
    std::unordered_map<int, std::string> _prefixes;
    std::unordered_map<int, std::string> _suffixes;
  };
  // struct for node
  struct node {
    /*
    node(const std::string& value, const std::unordered_map<uint, uint> &
    headConnections,
         const std::unordered_map<uint, uint> & tailConnections)
    : _value(value), _nextConnections(headConnections),
    _tailConnections(tailConnections) {}*/
    node(const std::string& value)
        : _value(value),
          _hasUnExploredEdges(false),
          _numOfUnExploredEdges(0),
          _balanced(false) {}
    std::string _value;
    // forward, for this graph head overlap, first is position of edge and
    // second is position of connection
    // std::unordered_map<uint, uint> _previousConnections;
    std::vector<std::pair<uint, uint>> _previousConnections;
    // backward, for this graph tail overlap, first is position of edge and
    // second is position of connection
    // std::unordered_map<uint, uint> _nextConnections;
    std::vector<std::pair<uint, uint>> _nextConnections;
    //
    bool _hasUnExploredEdges;
    uint _numOfUnExploredEdges;
    bool _balanced;
    void updateNumOfExploredEdges(const std::vector<edge>& edges);
    bool travel(std::vector<edge>& edges, randomGenerator& gen,
                uint32_t& currentNode, bool& lastWasExplorable);
  };
  // debGraph(const std::vector<readObject> & reads, aligner& alignerObj,
  // int minOverLap, int maxOverLap);
  // constructors
  debGraph() {}
  debGraph(const nhGraph& regularGraph);
  debGraph(const std::map<int, std::vector<int>>& byNumbersGraph);
  debGraph(const std::map<std::string, VecStr>& byStringGraph);

  void defaultContructorOptions() { _gen = randomGenerator(); }
  // memebers
  std::vector<edge> _edges;
  std::vector<node> _nodes;
  std::unordered_map<std::string, uint> _nodesPositions;
  bool _unexploredEdgesAvailable;
  randomGenerator _gen;
  std::vector<uint32_t> _path;
  // functions
  void travelGraphEulerianCycleStyle();
  void travelGraphEulerianPathStyle();
  void addEdge(uint32_t previousNodePos, uint32_t nextNodePos,
               std::string edgeName = "", std::string edgeValue = "");
  // bool growPath(uint32_t & currentNode, std::vector<uint> & explorableNodes);
  void growPath(uint32_t& currentNode, uint64_t& explorableNodes);
  void printAdjacencyTable(std::ostream& out) const;
};
} // namespace njhseq




