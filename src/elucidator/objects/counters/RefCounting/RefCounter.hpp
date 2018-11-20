#pragma once
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
/*
 * RefCounter.hpp
 *
 *  Created on: Dec 24, 2015
 *      Author: nick
 */

#include <TwoBit.h>
#include "elucidator/objects/counters/RefCounting/RefBaseCounter.hpp"
#include "elucidator/objects/counters/RefCounting/RefIndelCounter.hpp"

namespace njhseq {

class RefCounter {
public:
	RefCounter(const std::string & twoBitFilename);
	RefCounter(const RefCounter & that);
	TwoBit::TwoBitFile reader_;
	std::mutex mut_;
	uint32_t hardCutOff_ = 0; /**< Don't report anything on coverage of this or below */
private:
	RefBaseCounter baseCounter_;
	RefIndelCounter indelCounter_;
	std::unordered_map<std::string, std::string> referencesCache_;
public:
	void setHardCutOff(uint32_t hardCutOff);

	void increaseBaseCount(const std::string & refName, size_t refPos,
			bool plusStrand, bool highQuality, char base, uint32_t count);

	void increaseBaseCount(const std::string & refName, size_t refPos,
				bool plusStrand, char base, uint32_t count);

	void increaseIndelCount(const std::string & refName, size_t refPos,
			bool plusStrand, bool insertion, const std::string & gapStr,
			uint32_t count);

	void merge(const RefCounter & otherCounter);

	VecStr getRefNames() const;

	std::vector<size_t> getPositionsForRef(const std::string & refName) const;

	static std::string printInfoHeader() ;

	static std::string printInfoStratifiedByQualHeader();

	void printIndelSeqPosInfo(std::ostream & out, const std::string & refKey,
			size_t posKey);
	void printIndelSizePosInfo(std::ostream & out, const std::string & refKey,
			size_t posKey);

	void printPosInfo(std::ostream & out, const std::string & refKey,
			size_t posKey) ;

	void printPosInfoStratifiedByQuality(std::ostream & out, const std::string & refKey,
			size_t posKey);

	void printPosInfo(std::ostream & out, const std::string & refKey,
			size_t posKey, size_t depthCutOff, double freqCutOff);

	void printPosInfoStratifiedByQuality(std::ostream & out,
			const std::string & refKey, size_t posKey, size_t depthCutOff,
			double freqCutOff);

	void printInfo(std::ostream & out);

	void printInfoStratifiedByQual(std::ostream & out);

	void printInfo(std::ostream & out, const std::string & refKey, size_t start,
			size_t stop);

	void printInfoStratifiedByQual(std::ostream & out, const std::string & refKey,
			size_t start, size_t stop);

	void printInfo(std::ostream & out, size_t depthCutOff, double freqCutOff);

	void printInfoStratifiedByQual(std::ostream & out, size_t depthCutOff,
			double freqCutOff);

	void printInfo(std::ostream & out, const std::string & refKey, size_t start,
			size_t stop, size_t depthCutOff, double freqCutOff);

	void printInfoStratifiedByQual(std::ostream & out, const std::string & refKey,
			size_t start, size_t stop, size_t depthCutOff, double freqCutOff);

	static std::string printIndelSeqInfoHeader();
	void printIndelSeqInfo(std::ostream & out);
	void printIndelSeqInfo(std::ostream & out, const std::string & refKey,
			size_t start, size_t stop);

	static std::string printIndelSizeInfoHeader();
	void printIndelSizeInfo(std::ostream & out);
	void printIndelSizeInfo(std::ostream & out, const std::string & refKey,
			size_t start, size_t stop);

	substituteMatrix getRatesAgainstRef();

	substituteMatrix getRatesAgainstRef(const std::string & refName,
			size_t refStart, size_t refStop);

	substituteMatrix getRatesAgainstRefHq();

	substituteMatrix getRatesAgainstRefHq(const std::string & refName,
			size_t refStart, size_t refStop);

	void resetCounts();

};



}  // namespace njhseq


