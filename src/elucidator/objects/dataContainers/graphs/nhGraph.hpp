#pragma once
//
//  nhGraph.hpp
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


#include "elucidator/common.h"


#include <njhseq/alignment/aligner.h>
#include <njhseq/objects/seqObjects/readObject.hpp>


namespace njhseq {

class nhGraph {

 public:
  nhGraph(const std::vector<readObject>& reads, aligner& alignerObj,
          int minOverLap, int maxOverLap);
  struct readFixes {
    readFixes(const readObject& read, int minOverLap, int maxOverLap)
        : _read(read) {
      for (auto i : iter::range(minOverLap, maxOverLap + 1)) {
        _prefixes[i] = _read.seqBase_.seq_.substr(0, i);
        _suffixes[i] =
            _read.seqBase_.seq_.substr(read.seqBase_.seq_.length() - i);
      }
    }
    readObject _read;
    std::unordered_map<int, std::string> _prefixes;
    std::unordered_map<int, std::string> _suffixes;
  };
  struct node {
    node(const readObject& read, int minOverLap, int maxOverLap)
        : _read(read), _readEnds(readFixes(read, minOverLap, maxOverLap)) {}
    readObject _read;
    readFixes _readEnds;
    // forward, for this graph head overlap
    std::vector<std::pair<uint, std::string>> _headConnections;
    // backward, for this graph tail overlap
    std::vector<std::pair<uint, std::string>> _tailConnections;
  };
  struct connection {
    connection(const std::string& value, const std::pair<uint, uint>& firstCon)
        : _value(value), _nodePositions({firstCon}) {}
    std::string _value;
    // first being the tail connection possiton and the second being the head
    // connection position
    std::unordered_map<uint, uint> _nodePositions;
  };
  std::vector<node> _nodes;
  void addNode(const readObject& read, aligner& alignerObj, int minOverLap,
               int maxOverLap);
  void printAdjacencyTable(std::ostream& out, bool printName);
  std::multimap<uint, uint> getAdjacencyMap();
};

}  // namespace njhseq

