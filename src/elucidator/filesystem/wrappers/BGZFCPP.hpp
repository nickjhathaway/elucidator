#pragma once
/*
 * BGZFCPP.hpp
 *
 *  Created on: Jul 23, 2016
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


//due to a macro define clash with bamtools, htslib/bgzf.h must always come first in include orders
#include <htslib/kstring.h>
#include <htslib/bgzf.h>
#include "elucidator/common.h"


namespace njhseq {



class BGZFCPP {
	BGZF * file_;
	kstring_t kstr_ = { 0, 0, NULL };

public:

	const bfs::path fnp_;
	const std::string mode_;

	BGZFCPP(const bfs::path & fnp, const std::string & mode);

	template<typename T>
	bool readVec(std::vector<T>& data) {
		if (data.empty()) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__
					<< ": Error, trying to read into an empty vector, needs to have a size to determine how much data to read"
					<< "\n";
			throw std::runtime_error { ss.str() };
		}
		auto bytesRead = bgzf_read(file_, data.data(), sizeof(T) * data.size());
		if (bytesRead < 0) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ": Error in reading data from " << fnp_
					<< ", status: " << bytesRead << "\n";
			throw std::runtime_error { ss.str() };
		}
		if (bytesRead > 0) {
			if (static_cast<int64_t>(sizeof(T) * data.size()) != bytesRead) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ": Error in reading data from " << fnp_
						<< "\n";
				ss << "Read in " << bytesRead << " when expected "
						<< sizeof(T) * data.size() << "\n";
				throw std::runtime_error { ss.str() };
			}
			return true;
		}
		return false;
	}

	template<typename T>
	void writeVec(const std::vector<T>& data) {
		if(data.empty()){
			/**@todo maybe throw? */
			return;
		}
		auto bytesWritten = bgzf_write(file_, data.data(),
				sizeof(T) * data.size());
		if (bytesWritten < 0) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ": Error in writing data from " << fnp_
					<< ", status: " << bytesWritten << "\n";
			throw std::runtime_error { ss.str() };
		}
		if (static_cast<int64_t>(sizeof(T) * data.size()) != bytesWritten) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ": Error in writing data from " << fnp_
					<< "\n";
			ss << "Wrote " << bytesWritten << " when expected "
					<< sizeof(T) * data.size() << "\n";
			throw std::runtime_error { ss.str() };
		}
	}

	void writeStr(const std::string & str);

	long long int tell() const ;

	void seek(long long int pos) ;


	bool getline(std::string & line, int delim = '\n');

	~BGZFCPP();
};

}  // namespace njhseq



