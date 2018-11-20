
//
//  nhGraph.cpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 12/5/13.
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
#include "nhGraph.hpp"

namespace njhseq {
nhGraph::nhGraph(const std::vector<readObject>& reads, aligner& alignerObj,
                 int minOverLap, int maxOverLap) {
  // start nodes with no connections
  for (const auto& read : reads) {
    _nodes.emplace_back(node(read, minOverLap, maxOverLap));
    //_nodes.back()._read.setLetterCount();
  }
  // now connect the nodes
  std::map<int, int> countsOfOverlap;
  for (uint64_t pos : iter::range(_nodes.size())) {
    // bool foundAHeadConnection = false;
    // bool foundATailConnection = false;
    for (uint64_t search : iter::range(pos + 1, (uint64_t)_nodes.size())) {
      // need to find largest overLap
      std::string currentPrefix = "";
      std::string currentSuffix = "";
      bool foundBackConnection = false;
      bool foundForwardConnection = false;
      bool foundAnOverlap = false;
      for (auto i : iter::range(maxOverLap, minOverLap - 1, -1)) {
        if (_nodes[pos]._readEnds._suffixes[i] ==
            _nodes[search]._readEnds._prefixes[i]) {
          foundBackConnection = true;
          foundAnOverlap = true;
          currentPrefix = _nodes[search]._readEnds._prefixes[i];
          currentSuffix = _nodes[pos]._readEnds._suffixes[i];
          break;
        } else if (_nodes[pos]._readEnds._prefixes[i] ==
                   _nodes[search]._readEnds._suffixes[i]) {
          foundForwardConnection = true;
          foundAnOverlap = true;
          currentPrefix = _nodes[pos]._readEnds._prefixes[i];
          currentSuffix = _nodes[search]._readEnds._suffixes[i];
          break;
        }
      }
      if (!foundAnOverlap) {
        continue;
      }
      ++countsOfOverlap[alignerObj.comp_.highQualityMatches_];
      if (foundForwardConnection) {
        _nodes[pos]._headConnections.push_back({search, currentSuffix});
        _nodes[search]._tailConnections.push_back({pos, currentPrefix});
        // foundATailConnection = true;
      } else if (foundBackConnection) {
        _nodes[pos]._tailConnections.push_back({search, currentPrefix});
        _nodes[search]._headConnections.push_back({pos, currentSuffix});
        // foundAHeadConnection = true;
      }
    }
  }
  std::ofstream countsFile;
  openTextFile(countsFile, "counts", ".txt", true, false);
  countsFile << "lengthOfOverlap\tfreq" << std::endl;
  for (const auto& count : countsOfOverlap) {
    countsFile << count.first << "\t" << count.second << std::endl;
  }
}

void nhGraph::addNode(const readObject& read, aligner& alignerObj,
                      int minOverLap, int maxOverLap) {

  uint search = _nodes.size();
  _nodes.emplace_back(node(read, minOverLap, maxOverLap));
  for (uint64_t pos : iter::range(_nodes.size())) {
    // bool foundAHeadConnection = false;
    // bool foundATailConnection = false;
    std::string currentPrefix = "";
    std::string currentSuffix = "";
    bool foundBackConnection = false;
    bool foundForwardConnection = false;
    bool foundAnOverlap = false;
    for (auto i : iter::range(maxOverLap, minOverLap - 1, -1)) {
      if (_nodes[pos]._readEnds._suffixes[i] ==
          _nodes[search]._readEnds._prefixes[i]) {
        foundBackConnection = true;
        foundAnOverlap = true;
        currentPrefix = _nodes[search]._readEnds._prefixes[i];
        currentSuffix = _nodes[pos]._readEnds._suffixes[i];
        break;
      } else if (_nodes[pos]._readEnds._prefixes[i] ==
                 _nodes[search]._readEnds._suffixes[i]) {
        foundForwardConnection = true;
        foundAnOverlap = true;
        currentPrefix = _nodes[pos]._readEnds._prefixes[i];
        currentSuffix = _nodes[search]._readEnds._suffixes[i];
        break;
      }
    }
    if (!foundAnOverlap) {
      continue;
    }
    if (foundForwardConnection) {
      _nodes[pos]._headConnections.push_back({search, currentSuffix});
      _nodes[search]._tailConnections.push_back({pos, currentPrefix});
      // foundATailConnection = true;
    } else if (foundBackConnection) {
      _nodes[pos]._tailConnections.push_back({search, currentPrefix});
      _nodes[search]._headConnections.push_back({pos, currentSuffix});
      // foundAHeadConnection = true;
    }
  }
}
void nhGraph::printAdjacencyTable(std::ostream& out, bool printName) {
  std::multimap<std::string, std::string> ansOrganized;
  for (const auto& n : _nodes) {
    if (n._tailConnections.size() > 0) {
      for (auto i : n._tailConnections) {
        if (printName) {
          ansOrganized.insert({n._read.seqBase_.name_,
                               n._read.seqBase_.name_ + " " +
                                   _nodes[i.first]._read.seqBase_.name_});
        } else {
          ansOrganized.insert(
              {n._read.seqBase_.seq_, n._read.seqBase_.seq_ + " -> " +
                                          _nodes[i.first]._read.seqBase_.seq_});
        }
      }
    }
  }
  for (const auto& ans : ansOrganized) {
    out << ans.second << std::endl;
  }
}
std::multimap<uint, uint> nhGraph::getAdjacencyMap() {
  std::multimap<uint, uint> ansOrganized;
  uint pos = 0;
  for (const auto& n : _nodes) {
    if (n._tailConnections.size() > 0) {
      for (auto i : n._tailConnections) {
        ansOrganized.insert({pos, i.first});
      }
    }
    ++pos;
  }
  return ansOrganized;
}

}  // namespace njh
