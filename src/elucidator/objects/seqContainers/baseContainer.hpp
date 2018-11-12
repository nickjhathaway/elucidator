#pragma once
//
//  baseContainer.hpp
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

namespace njhseq {

template <typename T>
class baseContainer {

 public:
  // contructors
	baseContainer() :
			seqBase_(seqInfo()) {
	}

	baseContainer(const seqInfo& seqBase) :
			seqBase_(seqBase) {
	}

	baseContainer(const seqInfo& seqBase, const std::vector<T>& reads) :
			seqBase_(seqBase), reads_(reads) {
	}

	//copy constructor
	baseContainer(const baseContainer & other) :
			seqBase_(other.seqBase_), reads_(other.reads_) {
	}

	//move constructor
	baseContainer(baseContainer && other) :
			seqBase_(other.seqBase_), reads_(std::move(other.reads_)) {
	}

  // members
  seqInfo seqBase_;
  std::vector<T> reads_;
  std::mutex mut_;
  // functions
  // adding reads
	template<typename READ>
	void addRead(READ&& read) {
		std::lock_guard<std::mutex> lock(mut_);
		seqBase_.cnt_ += getSeqBase(read).cnt_;
		seqBase_.frac_ += getSeqBase(read).frac_;
		reads_.emplace_back(std::forward<READ>(read));
	}

	virtual void writeReads(SeqOutput & writer, bool writeBase) {
		std::lock_guard<std::mutex> lock(mut_);
		if (writeBase) {
			writer.openWrite(seqBase_);
		}
		for (const auto & read : reads_) {
			writer.openWrite(read);
		}
	}

	virtual ~baseContainer() {
	}
};

}  // namespace njhseq

