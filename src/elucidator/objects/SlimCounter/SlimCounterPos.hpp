#pragma once
/*
 * SlimCounterPos.hpp
 *
 *  Created on: Jun 16, 2016
 *      Author: nick
 */
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
#include "elucidator/objects/SlimCounter/SlimCounterCommon.hpp"

namespace njhseq {

class SlimCounterPos {
public:
	const static uint32_t NumOfElements = 26;
	SlimCounterPos();
	SlimCounterPos(const std::vector<uint32_t> & data);
	SlimCounterPos(const std::array<uint32_t, NumOfElements > & data);

	std::array<uint32_t, 20> baseCounts_;
	uint32_t fwInsertions_ = 0;
	uint32_t revInsertions_ = 0;
	uint32_t fwDeletions_ = 0;
	uint32_t revDeletions_ = 0;



	inline static uint32_t getBaseCountIndex(bool reverseStrand, bool highQuality) {
		return ((reverseStrand ? 1 : 0) + (highQuality ? 0 : 2)) * 5;
	}

	void increaseCount(char base, bool reverseStand, bool highQuality,
			uint32_t count);

	void increaseInsertionCount(bool reverseStrand,uint32_t count);

	void increaseDeletionCount(bool reverseStrand,uint32_t count);

	std::vector<uint32_t> genOutVec()const;

	static VecStr fullDetailHeader();
	static VecStr hqDetailHeader();

	void fullDetailOut(const std::string & refName, uint32_t pos,
			char refBase, std::ostream & out);

	std::vector<VecStr> fullDetailOutVec(const std::string & refName, uint32_t pos,
				char refBase);

	std::vector<VecStr> hqDetailOutVec(const std::string & refName, uint32_t pos,
					char refBase);

	uint32_t getEventsTotal() const;

	uint32_t getBaseTotal() const;

	//high quality

	uint32_t getHqBaseTotal() const;

	uint32_t getHqEventsTotal() const;

	uint32_t getHqFwBaseTotal() const;

	uint32_t getHqRevBaseTotal() const;

	uint32_t getHqFwBase(char base) const;

	uint32_t getHqRevBase(char base) const;

	uint32_t getHqBase(char base) const;

	double getHqFwBaseFrac(char base) const;

	double getHqBaseFrac(char base) const;

	//low quality

	uint32_t getLqBaseTotal() const;

	uint32_t getLqFwBaseTotal() const;

	uint32_t getLqRevBaseTotal() const;

	uint32_t getLqFwBase(char base) const ;

	uint32_t getLqRevBase(char base) const;

	uint32_t getLqBase(char base) const;

	double getLqFwBaseFrac(char base) const;

	double getLqBaseFrac(char base) const;

	double getHqFraction() const;

	//insertions

	uint32_t getFwInsertions() const;

	double getFwInsertionsFrac() const;

	uint32_t getRevInsertions() const;

	uint32_t getInsertions() const;

	double getInsertionsFrac() const;

	//deletions

	uint32_t getFwDeletions() const;

	double getFwDeletionsFrac() const;

	uint32_t getRevDeletions() const;

	uint32_t getDeletions() const;

	double getDeletionsFrac() const;


};

}  // namespace njhseq




