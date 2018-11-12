#pragma once
//
//  otuContainer.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 3/7/14.
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
#include "elucidator/objects/seqContainers/baseContainer.hpp"

namespace njhseq {

template <typename T>
class otuContainer : public baseContainer<T> {

 public:
  // contructors
  otuContainer() : baseContainer<T>() {}
  otuContainer(const seqInfo& seqBase) : baseContainer<T>(seqBase) {}
  otuContainer(const seqInfo& seqBase, const std::vector<T>& reads)
      : baseContainer<T>(seqBase, reads) {}
  otuContainer(const seqInfo& seqBase, double percentIdentity,
               double gapPercentCutOff)
      : percentIdentity_(percentIdentity),
        gapPercentCutOff_(gapPercentCutOff),
        baseContainer<T>(seqBase) {}
  otuContainer(const seqInfo& seqBase, const std::vector<T>& reads,
               double percentIdentity, double gapPercentCutOff)
      : percentIdentity_(percentIdentity),
        gapPercentCutOff_(gapPercentCutOff),
        baseContainer<T>(seqBase, reads) {}
  // members
  double percentIdentity_ = 0.97;
  double gapPercentCutOff_ = 0.03;

  // functions
};

}  // namespace njhseq


