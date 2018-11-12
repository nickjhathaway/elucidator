#pragma once
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
/*
 * RefBaseCounter.hpp
 *
 *  Created on: Dec 24, 2015
 *      Author: nick
 */

#include "elucidator/objects/counters/DNABaseCounter.hpp"

namespace njhseq {

class RefBaseCounter {
public:
	class BaseCount {
	public:
		BaseCount();
		/*
		enum class STRANDQUAL{
			POSHIGH,
			POSLOW,
			NEGHIGH,
			NEGLOW
		};*/

		//strand (true = plus strand), high vs low quality (true = high)
		std::unordered_map<bool, std::unordered_map<bool, DNABaseCounter>> counts_;



		void increaseCount(bool plusStrand, bool highQuality, char base,
				uint32_t count);
		void increaseCount(bool plusStrand, char base, uint32_t count);

		void resetAlphabets();
		void setFractions();
		std::vector<char> getAlphabet() const;

		uint64_t getTotalCount() const;

		uint64_t getHighQTotalCount() const;
		double getHighQFrac() const;

		uint64_t getLowQTotalCount() const;
		double getLowQFrac() const;

		uint64_t getPlusStrandTotalCount() const;
		double getPlusStrandFrac() const;

		uint64_t getNegStrandTotalCount() const;
		double getNegStrandFrac() const;

		uint64_t getHighQPlusStrandTotalCount() const;
		double getHighQFracPlusStrand() const;

		uint64_t getHighQNegStrandTotalCount() const;
		uint64_t getLowQPlusStrandTotalCount() const;
		double getHighQFracNegStrand() const;

		uint64_t getLowQNegStrandTotalCount() const;

		uint64_t getCount(char base) const;
		double getFrac(char base) const;

		uint64_t getHighQCount(char base) const;
		uint64_t getLowQCount(char base) const;
		double getHighQFrac(char base) const;

		double getBaseFracHighQ(char base) const;
		double getBaseFracLowQ(char base) const;

		uint64_t getPlusStrandCount(char base) const;
		uint64_t getNegStrandCount(char base) const;
		double getPlusStrandFrac(char base) const;

		uint64_t getHighQPlusStrandCount(char base) const;
		uint64_t getHighQNegStrandCount(char base) const;
		double getPlusStrandFracHighQ(char base) const;

		uint64_t getLowQPlusStrandCount(char base) const;
		uint64_t getLowQNegStrandCount(char base) const;
		double getPlusStrandFracLowQ(char base) const;

		static std::string infoStratifiedByQualityHeader(const std::string & delim);
		void infoStratifiedByQuality(std::ostream & out, char base,
				const std::string delim) const;
		static std::string infoHeader(const std::string & delim);
		void info(std::ostream & out, char base, const std::string delim) const;
		void merge(const BaseCount & otherCounter);

		void increaseRates(substituteMatrix & mat, char base) const;
		void increaseRatesHq(substituteMatrix & mat, char base) const;
		void increaseRatesLq(substituteMatrix & mat, char base) const;
		void increaseRatesForward(substituteMatrix & mat, char base) const;
		void increaseRatesReverse(substituteMatrix & mat, char base) const;
		void increaseRatesHqForward(substituteMatrix & mat, char base) const;
		void increaseRatesLqForward(substituteMatrix & mat, char base) const;
		void increaseRatesHqReverse(substituteMatrix & mat, char base) const;
		void increaseRatesLqReverse(substituteMatrix & mat, char base) const;

	};

	//goes refId, refPos, base count
	std::unordered_map<std::string, std::unordered_map<size_t, BaseCount>> counts_;
	bool needsToSet = true;
	void setAlphabetAndFractions();
	void increaseCount(const std::string & refName, size_t refPos,
			bool plusStrand, bool highQuality, char base, uint32_t count);
	void increaseCount(const std::string & refName, size_t refPos,
			bool plusStrand, char base, uint32_t count);
	VecStr getRefNames() const;
	std::vector<size_t> getPositionsForRef(const std::string & refName) const;

	void merge(const RefBaseCounter & otherCounter);
};

}  // namespace njhseq

