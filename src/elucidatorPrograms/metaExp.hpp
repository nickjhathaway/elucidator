#pragma once
/*
 * metaExp.hpp
 *
 *  Created on: Mar 14, 2018
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


#include <njhcpp/progutils.h>
#include <njhseq.h>

namespace njhseq {

class metaExpRunner : public njh::progutils::ProgramRunner {
 public:
	metaExpRunner();

	static int addMetaBySampleName(const njh::progutils::CmdArgs & inputCommands);
	static int addMetaByMetaField(const njh::progutils::CmdArgs & inputCommands);
	static int excludeSeqsFileWithNumericMetaCutOff(const njh::progutils::CmdArgs & inputCommands);
	static int excludeSeqsFileWithMatchingMeta(const njh::progutils::CmdArgs & inputCommands);
	static int selectMetaFieldsToKeep(const njh::progutils::CmdArgs & inputCommands);
	static int splitSeqFileWithMeta(const njh::progutils::CmdArgs & inputCommands);
	static int splitSeqFileWithExternalMeta(const njh::progutils::CmdArgs & inputCommands);
  static int createTableFromSeqs(const njh::progutils::CmdArgs & inputCommands);
  static int printMetaFieldsFromSeqs(const njh::progutils::CmdArgs & inputCommands);

};

} //  namespace njhseq



