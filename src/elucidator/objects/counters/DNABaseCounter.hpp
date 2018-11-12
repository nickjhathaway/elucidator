#pragma once
/*
 * DNABaseCounter.hpp
 *
 *  Created on: Feb 26, 2016
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

/**@brief Class to hold counts for letters, by default upper case A to lowercase z
 *
 */
template<typename T, char START = 'A', char END = 'z'>
class LetArray {
	const char start_ = START;
	const char end_ = END;
	std::array<T, END - START> bases_;
public:
	// Element access.
	T& operator[](char base) noexcept {
		return bases_[base - start_];
	}

	constexpr const T& operator[](char base) const noexcept {
		return bases_[base - start_];
	}

	T& at(char base) {
		return bases_.at(base - start_);
	}

	constexpr const T& at(char base) const {
		return bases_.at(base - start_);
	}

  void fill(const T& val) {
		bases_.fill(val);
	}

	//delete some functions to prevent using things that can be cast to char
	T& operator[](int32_t base) noexcept = delete;
	constexpr const T& operator[](int32_t base) const noexcept = delete;
	T& at(int32_t base) = delete;
	constexpr const T& at(int32_t base) const = delete;
	T& operator[](uint32_t base) noexcept = delete;
	constexpr const T& operator[](uint32_t base) const noexcept = delete;
	T& at(uint32_t base) = delete;
	constexpr const T& at(uint32_t base) const = delete;
	T& operator[](size_t base) noexcept = delete;
	constexpr const T& operator[](size_t base) const noexcept = delete;
	T& at(size_t base) = delete;
	constexpr const T& at(size_t base) const = delete;


	static std::vector<char> bases(){
		std::vector<char> bases(END - START);
		std::iota(bases.begin(), bases.end(), START);
		return bases;
	}
};



class DNABaseCounter {
public:
	DNABaseCounter();
	DNABaseCounter(std::set<char> alphabet);

	LetArray<uint32_t> bases_;
	LetArray<double> baseFracs_;

	std::set<char> alphabet_ = { 'A', 'C', 'G', 'T', 'N' };

	const std::set<char> originalAlphabet_ = { 'A', 'C', 'G', 'T', 'N' };

	/**@brief convert to json representation
	 *
	 * @return Json::Value object
	 */
	Json::Value toJson() const;

	void resetAlphabet(bool keepOld);
	void reset();


	void increase(char base, uint32_t cnt = 1);
	void increase(const std::string &seq, uint32_t cnt = 1);
	//increase by seq portion, can be any range over a collection of chars
	template<class InputIt1>
	void increaseWithRange(InputIt1 first1, InputIt1 last1, uint32_t cnt = 1) {
		for (auto iter = first1; iter < last1; ++iter) {
			increase(*iter, cnt);
		}
	}
	void increaseWithSubStr(const std::string & str, size_t pos, size_t len,
			uint32_t cnt = 1);

	void setFractions();
	void setFractions(const std::set<char> & alphabet);
	void addOtherCounts(const DNABaseCounter & otherCounter, bool setFractions);
	uint32_t getTotalCount() const;
	uint32_t getTotalCount(const std::set<char> & alphabet) const;
	std::multimap<double, char, std::less<double>> createLikelihoodMaps(
			bool setFractionFirst);

	double calcGcContent() const;
	int getGcDifference() const;
	// compute entropy
	double computeEntrophy() const;

	char getDegenativeBase() const;
	// output data
	static std::string outPutInfoHeader();
	void outPutInfo(std::ostream &out, bool header) const;
	double getFracDifference(const DNABaseCounter & otherCounter,
			const std::vector<char> & alph) const;
	void increaseRates(substituteMatrix & mat, char refBase) const;

};

}  // namespace njhseq
