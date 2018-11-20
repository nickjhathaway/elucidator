/*
 * LockableJsonLog.cpp
 *
 *  Created on: Mar 18, 2017
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

#include "LockableJsonLog.hpp"

namespace njhseq {


LockableJsonLog::LockableJsonLog(const bfs::path & logFnp, bool overWrite) :
		logFileOpts_(logFnp) {
	logFileOpts_.overWriteFile_ = overWrite;
	log_["date"] = njh::json::toJson(njh::getCurrentDateFull());
}



void LockableJsonLog::addToLog(const std::string & uid, const Json::Value & adding){
	std::lock_guard<std::mutex> lock(logMut_);
	if(log_.isMember(uid)){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error, log already has info for " << uid << "\n";
		throw std::runtime_error{ss.str()};
	}
	log_[uid] = adding;
}


void LockableJsonLog::writeLog(){
	std::lock_guard<std::mutex> lock(logMut_);
	std::ofstream outFile;
	logFileOpts_.openFile(outFile);
	log_["totalTime"] = njh::json::toJson(watch_.totalTime());
	outFile << log_ << std::endl;;
}


} /* namespace njhseq */
