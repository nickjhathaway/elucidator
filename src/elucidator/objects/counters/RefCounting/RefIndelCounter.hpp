#pragma once
/*
 * RefIndelCounter.hpp
 *
 *  Created on: Dec 24, 2015
 *      Author: nick
 */
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

class RefIndelCounter {

public:
	class IndelCounter {
	public:
		IndelCounter();
		//strand (true = + strand), indel (true = insertion)
		std::unordered_map<bool, std::unordered_map<bool, strCounter>> counts_;

		void increaseCount(bool plusStrand, bool insertion,
				const std::string & gapStr, uint32_t count = 1);

		uint32_t getInsTotal() const;
		uint32_t getInsTotalPlusStrand() const;
		uint32_t getInsTotalNegStrand() const;
		double getInsFracPlusStrand() const;

		uint32_t getDelTotal() const;
		uint32_t getDelTotalPlusStrand() const;
		uint32_t getDelTotalNegStrand() const;
		double getDelFracPlusStrand() const;


		static std::string infoHeader(const std::string & delim);
		static std::string sizeInfoHeader(const std::string & delim);
		static std::string seqInfoHeader(const std::string & delim);

		void insertionInfo(std::ostream & out, const std::string delim) const;
		void insertionSizeInfo(std::ostream & out,const std::string refName,
				size_t refPos, char refBase,  const std::string delim) const;
		void insertionSeqInfo(std::ostream & out,const std::string refName,
				size_t refPos, char refBase,  const std::string delim) const;

		void deletionInfo(std::ostream & out, const std::string delim) const;
		void deletionSizeInfo(std::ostream & out,const std::string refName,
				size_t refPos, char refBase, const std::string delim) const;
		void deletionSeqInfo(std::ostream & out,const std::string refName,
				size_t refPos, char refBase,  const std::string delim) const;

		void merge(const IndelCounter & otherCounter);

	};

	//goes refId, refPos, indel count
	std::unordered_map<std::string, std::unordered_map<size_t, IndelCounter> > counts_;

	void increaseCount(const std::string & refName, size_t refPos,
			bool plusStrand, bool insertion, const std::string & gapStr,
			uint32_t count);
	VecStr getRefNames() const;
	std::vector<size_t> getPositionsForRef(const std::string & refName) const;
	void merge(const RefIndelCounter & otherCounter);
};




}  // namespace njhseq


