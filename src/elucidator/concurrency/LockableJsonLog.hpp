#pragma once
/*
 * LockableJsonLog.hpp
 *
 *  Created on: Mar 18, 2017
 *      Author: nick
 */

#include "elucidator/utils.h"

namespace njhseq {

class LockableJsonLog {

public:

	LockableJsonLog(const bfs::path & logFnp, bool overWrite = true);

	Json::Value log_;
	std::mutex logMut_;

	OutOptions logFileOpts_;

	njh::stopWatch watch_;

	void addToLog(const std::string & uid, const Json::Value & adding);

	void writeLog();
};

} /* namespace njhseq */

