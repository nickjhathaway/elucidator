/*
 * BGZFCPP.cpp
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


#include "BGZFCPP.hpp"

namespace njhseq {

BGZFCPP::BGZFCPP(const bfs::path & fnp, const std::string & mode) :
		fnp_(fnp), mode_(mode) {
	if (mode != "r" && mode != "w" && mode != "a") {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": Error, mode can only be r,w, or a not"
				<< mode << "\n";
		throw std::runtime_error { ss.str() };
	}
	file_ = bgzf_open(fnp_.c_str(), mode.c_str());
	if ("r" == mode) {
		auto status = bgzf_is_bgzf(fnp_.c_str());
		if (0 == status) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ": Error in reading from " << fnp_
					<< " it isn't in bgzf format or is malformed" << "\n";
			throw std::runtime_error { ss.str() };
		}
	}
}

void BGZFCPP::writeStr(const std::string & str) {
	if (str.empty()) {
		/**@todo maybe throw or at least warn?*/
		return;
	}
	auto bytesWritten = bgzf_write(file_, str.c_str(), strlen(str.c_str()));
	if (bytesWritten < 0) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": Error in writing string, " << str
				<< " from " << fnp_ << ", code: " << bytesWritten << "\n";
		throw std::runtime_error { ss.str() };
	} else if (strlen(str.c_str()) != static_cast<size_t>(bytesWritten)) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": Error in writing string, " << str
				<< " from " << fnp_ << ", bytes written doesn't equal the string length"
				<< "\n";
		ss << "strlen: " << strlen(str.c_str()) << ", bytesWritten: "
				<< bytesWritten << "\n";
		throw std::runtime_error { ss.str() };
	}
}

long long int BGZFCPP::tell() const {
	return bgzf_tell(file_);
}

void BGZFCPP::seek(long long int pos) {
	auto res = bgzf_seek(file_, pos, SEEK_SET);
	if (-1 == res) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": Error in seeking " << fnp_ << "\n";
		throw std::runtime_error { ss.str() };
	}
}


bool BGZFCPP::getline(std::string & line, int delim) {
	auto status = bgzf_getline(file_, '\n', &kstr_);
	if (status > 0) {
		line = std::string(kstr_.s);
	}
	if(status < -1){
		/** not sure why but it seems like status of -1
		 *  is return from bgzf_getline when at the end of the file*/
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": Error in reading " << fnp_  << ", code: " << status << "\n";
		throw std::runtime_error { ss.str() };
	}
	return status > 0;
}

BGZFCPP::~BGZFCPP() {
	free(kstr_.s);
	bgzf_close(file_);
}

}  // namespace njhseq
