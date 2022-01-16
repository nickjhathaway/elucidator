#pragma once
/*
 * seqSearching.hpp
 *
 *  Created on: Nov 12, 2018
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
#include <njhseq/common.h>
#include <njhseq/programUtils/seqSetUp.hpp>

namespace njhseq {



class seqSearchingRunner : public njh::progutils::ProgramRunner {
 public:
	seqSearchingRunner();


	static int chopAndMap(const njh::progutils::CmdArgs & inputCommands);
	static int chopAndMapAndRefine(const njh::progutils::CmdArgs & inputCommands);
	static int chopAndMapAndRefineInvidual(const njh::progutils::CmdArgs & inputCommands);

  static int findHomopolymerLocations(const njh::progutils::CmdArgs & inputCommands);


  static int findMotifLocations(const njh::progutils::CmdArgs & inputCommands);
  static int findTandemMotifLocations(const njh::progutils::CmdArgs & inputCommands);


  static int findSimpleTandemRepeatLocations(const njh::progutils::CmdArgs & inputCommands);


};

}  // namespace njhseq


