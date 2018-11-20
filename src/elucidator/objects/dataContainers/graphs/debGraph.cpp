//
//  debGraph.cpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 12/11/13.
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
#include "debGraph.hpp"
namespace njhseq {

void debGraph::addEdge(uint32_t previousNodePos, uint32_t nextNodePos,
                       std::string edgeName, std::string edgeValue) {
  if ("" == edgeValue) {
    edgeValue = _nodes[previousNodePos]._value + _nodes[nextNodePos]._value;
  }
  if ("" == edgeName) {
    edgeName = _nodes[previousNodePos]._value + _nodes[nextNodePos]._value;
  }
  _edges.push_back(edge(readObject(seqInfo(edgeValue, edgeValue))));
  _nodes[previousNodePos]
      ._nextConnections.push_back({_edges.size() - 1, nextNodePos});
  _nodes[nextNodePos]
      ._previousConnections.push_back({_edges.size() - 1, previousNodePos});
  return;
}
debGraph::debGraph(const nhGraph& regularGraph) {
  defaultContructorOptions();
  for (const auto& currentNode : regularGraph._nodes) {
    _edges.push_back(edge(currentNode._readEnds._read,
                          currentNode._readEnds._prefixes,
                          currentNode._readEnds._suffixes));
    for (const auto& suf : currentNode._readEnds._suffixes) {
      if (_nodesPositions.find(suf.second) == _nodesPositions.end()) {
        _nodes.push_back(node(suf.second));
        _nodesPositions.insert({suf.second, _nodes.size() - 1});
      }
    }
    for (const auto& pre : currentNode._readEnds._prefixes) {
      if (_nodesPositions.find(pre.second) == _nodesPositions.end()) {
        _nodes.push_back(node(pre.second));
        _nodesPositions.insert({pre.second, _nodes.size() - 1});
      }
    }
  }
  uint32_t pos = 0;
  for (const auto& currentNode : regularGraph._nodes) {
    if (currentNode._tailConnections.empty()) {
      for (const auto& suf : currentNode._readEnds._suffixes) {
        for (const auto& pre : currentNode._readEnds._prefixes) {
          _nodes[_nodesPositions.at(suf.second)]._previousConnections.push_back(
              {pos, _nodesPositions.at(pre.second)});
        }
      }
    } else {
      for (const auto& tail : currentNode._tailConnections) {
        for (const auto& suf :
             regularGraph._nodes[tail.first]._readEnds._suffixes) {
          _nodes[_nodesPositions.at(tail.second)]._nextConnections.push_back(
              {tail.first, _nodesPositions.at(suf.second)});
          //_nodes[_nodesPositions.at(suf.second)]._previousConnections.insert({pos,
          //_nodesPositions.at(tail.second)});
        }
      }
    }
    if (currentNode._headConnections.empty()) {
      for (const auto& pre : currentNode._readEnds._prefixes) {
        for (const auto& suf : currentNode._readEnds._suffixes) {
          _nodes[_nodesPositions.at(pre.second)]._nextConnections.push_back(
              {pos, _nodesPositions.at(suf.second)});
        }
      }
    } else {
      for (const auto& head : currentNode._headConnections) {
        for (const auto& pre :
             regularGraph._nodes[head.first]._readEnds._prefixes) {
          _nodes[_nodesPositions.at(head.second)]
              ._previousConnections.push_back(
                   {head.first, _nodesPositions.at(pre.second)});
          //_nodes[_nodesPositions.at(pre.second)]._nextConnections.insert({pos,
          //_nodesPositions.at(head.second)});
        }
      }
    }
    ++pos;
  }
}
debGraph::debGraph(const std::map<int, std::vector<int>>& byNumbers) {
  defaultContructorOptions();
  for (const auto& num : byNumbers) {
    std::string currentNode = std::to_string(num.first);
    if (_nodesPositions.find(currentNode) == _nodesPositions.end()) {
      _nodes.push_back(node(currentNode));
      _nodesPositions.insert({currentNode, _nodes.size() - 1});
    }
    for (const auto& subNum : num.second) {
      std::string currentSubNode = std::to_string(subNum);
      if (_nodesPositions.find(currentSubNode) == _nodesPositions.end()) {
        _nodes.push_back(node(currentSubNode));
        _nodesPositions.insert({currentSubNode, _nodes.size() - 1});
      }
      addEdge(_nodesPositions[currentNode], _nodesPositions[currentSubNode]);
      /*
      _edges.push_back(edge(readObject(std::to_string(num.first)+std::to_string(subNum),
      "")));
      _nodes[_nodesPositions[currentNode]]._nextConnections.push_back({_edges.size()-1,
      _nodesPositions[currentSubNode]});
      _nodes[_nodesPositions[currentSubNode]]._previousConnections.push_back({_edges.size()-1,
      _nodesPositions[currentNode]});*/
    }
  }
}
debGraph::debGraph(const std::map<std::string, VecStr>& byStringGraph) {
  defaultContructorOptions();
  for (const auto& str : byStringGraph) {
    std::string currentNode = str.first;
    trimEndWhiteSpace(currentNode);
    if (_nodesPositions.find(currentNode) == _nodesPositions.end()) {
      _nodes.push_back(node(currentNode));
      _nodesPositions.insert({currentNode, _nodes.size() - 1});
    }
    for (const auto& subStr : str.second) {
      std::string currentSubNode = subStr;
      trimEndWhiteSpace(currentSubNode);
      if (_nodesPositions.find(currentSubNode) == _nodesPositions.end()) {
        _nodes.push_back(node(currentSubNode));
        _nodesPositions.insert({currentSubNode, _nodes.size() - 1});
      }
      _edges.push_back(
          edge(readObject(seqInfo(currentNode + currentSubNode, ""))));
      _nodes[_nodesPositions[currentNode]]._nextConnections.push_back(
          {_edges.size() - 1, _nodesPositions[currentSubNode]});
      _nodes[_nodesPositions[currentSubNode]]._previousConnections.push_back(
          {_edges.size() - 1, _nodesPositions[currentNode]});
    }
  }
}
void debGraph::printAdjacencyTable(std::ostream& out) const {
  std::multimap<std::string, std::string> ansOrganized;
  for (const auto& currentNode : _nodes) {
    if (currentNode._nextConnections.empty()) {
      continue;
    }
    VecStr cons;
    for (const auto& connection : currentNode._nextConnections) {
      cons.push_back(_nodes[connection.second]._value);
    }
    ansOrganized.insert({currentNode._value, currentNode._value + " -> " +
                                                 vectorToString(cons, ",")});
  }
  for (const auto& ans : ansOrganized) {
    out << ans.second << std::endl;
  }
}

void debGraph::travelGraphEulerianCycleStyle() {
  // assuming for the moment strongly connect and completely balanced
  // reset everything in the grpah
  uint64_t numOfUnExploredEdges = _edges.size();
  _path.clear();
  _path.reserve(numOfUnExploredEdges + 2);
  for (auto& ed : _edges) {
    ed._read.remove = false;
  }
  for (auto& no : _nodes) {
    no.updateNumOfExploredEdges(_edges);
  }
  // assuming node has an edge
  uint nodePos = 0;
  uint32_t currentNode = 0;
  for (const auto& n : _nodes) {
    std::cout << "node " << _nodes[nodePos]._value << " is balanaced "
              << convertBoolToString(_nodes[nodePos]._balanced) << std::endl;
    std::cout << "\tpConSize: " << n._previousConnections.size()
              << " nConSize: " << n._nextConnections.size() << std::endl;
    if (!_nodes[nodePos]._balanced) {
      if (_nodes[nodePos]._previousConnections.size() <
          _nodes[nodePos]._nextConnections.size()) {
        currentNode = nodePos;
      }
    }
    ++nodePos;
  }
  std::cout << "no previous connections nodePos: " << currentNode << std::endl;
  // std::cout<<"starting node: "<<_nodes[currentNode]._value<<std::endl;
  currentNode = _gen.unifRand(0, (int)_nodes.size());
  _path = {currentNode};

  // to store the explorable node positions in the growing path
  // std::vector<uint> explorableNodes;
  bool keepGrowing = true;
  while (keepGrowing) {
    std::cout << "numOfUnExploredEdges: " << numOfUnExploredEdges << std::endl;
    if (numOfUnExploredEdges == 0) {
      std::cout << "no more explorable edges left" << std::endl;
      keepGrowing = false;
    }
    if (!keepGrowing) {
      break;
    }
    // grow path
    growPath(currentNode, numOfUnExploredEdges);
    // collect exporable nodes
    std::vector<uint> realExplorableNodes;
    uint epos = 0;
    for (const auto& p : _path) {
      if (_nodes[p]._hasUnExploredEdges) {
        realExplorableNodes.push_back(epos);
      }
      ++epos;
    }
    if (realExplorableNodes.empty()) {
      std::cout << "no more nodes with explorable edges" << std::endl;
      keepGrowing = false;
    }
    if (!keepGrowing) {
      break;
    }
    _path.pop_back();
    uint nextNode = _gen.unifRand(0, (int)realExplorableNodes.size());
    // std::cout<<"current path "<<std::endl;
    /*for(const auto & p : _path){
     std::cout<<_nodes[p]._value<<"->";
     }
     std::cout<<std::endl;*/
    // printVector(_path,"->");
    // std::cout<<"NextNode: "<< nextNode<<" sizeOfexplorableNodes:
    // "<<realExplorableNodes.size()<<std::endl;
    auto pathPos = realExplorableNodes[nextNode];
    currentNode = _path[pathPos];
    // std::cout<<"pathPos: "<<pathPos<<std::endl;
    // std::cout<<"new current node with unexplored edge:
    // "<<currentNode<<std::endl;
    _path = concatVecs(getSubVector(_path, pathPos, _path.size() - pathPos),
                            getSubVector(_path, 0, pathPos + 1));
    // std::cout<<"new path orientation "<<std::endl;
    /*for(const auto & p : _path){
      std::cout<<_nodes[p]._value<<"->";
    }
    std::cout<<std::endl;*/
    // printVector(_path,"->");
    // std::cout<<std::endl;
  }
  VecStr ans;
  for (const auto& p : _path) {
    ans.push_back(_nodes[p]._value);
  }
  printVector(ans, "->");
  return;
}

void debGraph::travelGraphEulerianPathStyle() {
  // assuming for the moment strongly connect and completely balanced
  // reset everything in the grpah
  uint64_t numOfUnExploredEdges = _edges.size();
  _path.clear();
  _path.reserve(numOfUnExploredEdges + 2);
  for (auto& ed : _edges) {
    ed._read.remove = false;
  }
  for (auto& no : _nodes) {
    no.updateNumOfExploredEdges(_edges);
  }
  // assuming node has an edge
  uint32_t nodePos = 0;
  uint32_t currentNode = 0;
  uint32_t nodeNeedsPrevious = 0;
  uint32_t nodeNeedsNext = 0;
  for (const auto& n : _nodes) {
    std::cout << "node " << _nodes[nodePos]._value << " is balanaced "
              << convertBoolToString(_nodes[nodePos]._balanced) << std::endl;
    std::cout << "\tpConSize: " << n._previousConnections.size()
              << " nConSize: " << n._nextConnections.size() << std::endl;
    if (!_nodes[nodePos]._balanced) {
      if (_nodes[nodePos]._previousConnections.size() <
          _nodes[nodePos]._nextConnections.size()) {
        nodeNeedsPrevious = nodePos;
      } else if (_nodes[nodePos]._previousConnections.size() >
                 _nodes[nodePos]._nextConnections.size()) {
        nodeNeedsNext = nodePos;
      }
    }
    ++nodePos;
  }
  // now create a new edge to connect the ends and update
  addEdge(nodeNeedsNext, nodeNeedsPrevious);
  /*
  _edges.push_back(edge(readObject(_nodes[nodeNeedsNext]._value+_nodes[nodeNeedsPrevious]._value,
  "")));
  _nodes[nodeNeedsNext]._nextConnections.push_back({_edges.size()-1,
  nodeNeedsPrevious});
  _nodes[nodeNeedsPrevious]._previousConnections.push_back({_edges.size()-1,
  nodeNeedsNext});*/

  _nodes[nodeNeedsNext].updateNumOfExploredEdges(_edges);
  _nodes[nodeNeedsPrevious].updateNumOfExploredEdges(_edges);
  // std::cout<<"no previous connections nodePos: "<<currentNode<<std::endl;
  // std::cout<<"starting node: "<<_nodes[currentNode]._value<<std::endl;
  _gen.seed();
  currentNode = _gen.unifRand(0, (int)_nodes.size());
  _path = {currentNode};

  // to store the explorable node positions in the growing path
  // std::vector<uint> explorableNodes;
  bool keepGrowing = true;
  while (keepGrowing) {
    std::cout << "numOfUnExploredEdges: " << numOfUnExploredEdges << std::endl;
    if (numOfUnExploredEdges == 0) {
      std::cout << "no more explorable edges left" << std::endl;
      keepGrowing = false;
    }
    if (!keepGrowing) {
      break;
    }
    // grow path
    growPath(currentNode, numOfUnExploredEdges);
    // collect exporable nodes
    std::vector<uint> realExplorableNodes;
    uint32_t epos = 0;
    for (const auto& p : _path) {
      if (_nodes[p]._hasUnExploredEdges) {
        realExplorableNodes.push_back(epos);
      }
      ++epos;
    }
    if (realExplorableNodes.empty()) {
      std::cout << "no more nodes with explorable edges" << std::endl;
      keepGrowing = false;
    }
    if (!keepGrowing) {
      break;
    }
    _path.pop_back();
    uint32_t nextNode = _gen.unifRand<uint32_t>(0, realExplorableNodes.size());
    // std::cout<<"current path "<<std::endl;
    /*for(const auto & p : _path){
     std::cout<<_nodes[p]._value<<"->";
     }
     std::cout<<std::endl;*/
    // printVector(_path,"->");
    // std::cout<<"NextNode: "<< nextNode<<" sizeOfexplorableNodes:
    // "<<realExplorableNodes.size()<<std::endl;
    auto pathPos = realExplorableNodes[nextNode];
    currentNode = _path[pathPos];
    // std::cout<<"pathPos: "<<pathPos<<std::endl;
    // std::cout<<"new current node with unexplored edge:
    // "<<currentNode<<std::endl;
    _path = concatVecs(getSubVector(_path, pathPos, _path.size() - pathPos),
                            getSubVector(_path, 0, pathPos + 1));
    // std::cout<<"new path orientation "<<std::endl;
    /*for(const auto & p : _path){
     std::cout<<_nodes[p]._value<<"->";
     }
     std::cout<<std::endl;*/
    // printVector(_path,"->");
    // std::cout<<std::endl;
  }
  /*
  VecStr ans;
  for(const auto & p : _path){
    ans.push_back(_nodes[p]._value);
  }*/
  // printVector(ans, "->");
  // ans.clear();
  auto positions = getPositionsOfTarget(_path, nodeNeedsNext);
  uint32_t theEnd = 0;
  for (const auto& endPos : positions) {
    if (endPos + 1 < _path.size()) {
      if (_path[endPos + 1] == nodeNeedsPrevious) {
        theEnd = endPos;
      }
    }
  }
  // rearrange path
  _path =
      concatVecs(getSubVector(_path, theEnd + 1, _path.size() - 2 - theEnd),
                      getSubVector(_path, 0, theEnd + 1));
  /*for(const auto & p : _path){
    ans.push_back(_nodes[p]._value);
  }*/
  // printVector(ans, "->");
  // printVector(ans, "");
  return;
}

void debGraph::node::updateNumOfExploredEdges(const std::vector<edge>& edges) {
  if (_nextConnections.size() == _previousConnections.size()) {
    _balanced = true;
  } else {
    _balanced = false;
  }
  _numOfUnExploredEdges = 0;
  for (const auto& con : _nextConnections) {
    if (!edges.at(con.first)._read.remove) {
      _hasUnExploredEdges = true;
      ++_numOfUnExploredEdges;
    }
  }
}
void debGraph::growPath(uint32_t& currentNode, uint64_t& numOfUnExploredEdges) {
  bool lastWasExplorable = false;
  while (_nodes[currentNode]
             .travel(_edges, _gen, currentNode, lastWasExplorable)) {
    /*if(lastWasExplorable){
      std::cout<<"last one was explorable and now on
    "<<_nodes[currentNode]._value<<std::endl;
      explorableNodes.push_back(_path.size()-1);
      std::cout<<"Explorable nodes: ";printVector(explorableNodes);
    }*/
    --numOfUnExploredEdges;
    _path.push_back(currentNode);
    /*
    for(const auto & p : _path){
      std::cout<<_nodes[p]._value<<"->";
    }
    std::cout<<std::endl;
    printVector(_path,"->");
    std::cout<<std::endl;*/
  }
  return;
  // return !explorableNodes.empty();
}
bool debGraph::node::travel(std::vector<edge>& edges, randomGenerator& gen,
                            uint32_t& currentNode, bool& lastWasExplorable) {
  if (_hasUnExploredEdges) {
    bool keepSearching = true;
    while (keepSearching) {
      auto pos = gen.unifRand(0, (int)_nextConnections.size());
      if (!edges[_nextConnections[pos].first]._read.remove) {
        // stop searching, this is supper ineffecient to randomly find an
        // unexplored edge
        keepSearching = false;
        // turn off edge
        edges[_nextConnections[pos].first]._read.remove = true;
        // decrease edge count and check
        --_numOfUnExploredEdges;
        // std::cout<<"Node: "<<_value<<" now has "<<_numOfUnExploredEdges<<"
        // unexplored edges"<<std::endl;
        if (_numOfUnExploredEdges == 0) {
          // std::cout<<"Node: "<<_value<<" no longer has explorable
          // edges"<<std::endl;
          // std::cout<<"Number "
          _hasUnExploredEdges = false;
          lastWasExplorable = false;
        } else {
          lastWasExplorable = true;
        }
        // supply next node
        currentNode = _nextConnections[pos].second;
      }
    }
    return true;
  } else {
    lastWasExplorable = false;
    // can't travel
    return false;
  }
}
}
