#pragma once
//
//  refMapContainer.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 3/7/14.
//  Copyright (c) 2014 Nicholas Hathaway. All rights reserved.
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
#include "elucidator/objects/seqContainers/baseContainer.hpp"

namespace njhseq {

template <typename T>
class refMapContainer : public baseContainer<T> {

 public:
  // Constructors
  refMapContainer() : baseContainer<T>() {}
  refMapContainer(const seqInfo& seqBase) : baseContainer<T>(seqBase) {
    initialize();
  }
  refMapContainer(const seqInfo& seqBase, const std::vector<T>& reads)
      : baseContainer<T>(seqBase, reads) {
    initialize();
  }

  void initialize() {
    // no reads mapped yet so initialize counts to zero
    this->seqBase_.cnt_ = 0;
    this->seqBase_.frac_ = 0;
  }
  // members

  // functions
  template<typename READS, typename REFS>
  static std::vector<refMapContainer> createContainers(const std::vector<REFS> & reads){
  	std::vector<refMapContainer> refContainers;
  	for(const auto & read : reads){
  		refContainers.emplace_back(read.seqBase_);
  	}
  	return refContainers;
  }
};

}  // namespace njhseq


