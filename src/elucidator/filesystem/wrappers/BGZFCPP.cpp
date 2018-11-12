/*
 * BGZFCPP.cpp
 *
 *  Created on: Jul 23, 2016
 *      Author: nick
 */



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
