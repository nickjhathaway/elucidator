#pragma once
/*
 * kmerExpRunner.h
 *
 *  Created on: Feb 4, 2015
 *      Author: nickhathaway
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

class kmerExpRunner : public njh::progutils::ProgramRunner {
 public:
	kmerExpRunner();


	static int convertKmerSearchToBinaryMatrix(const njh::progutils::CmdArgs & inputCommands);
	static int kmerSearch(const njh::progutils::CmdArgs & inputCommands);
	static int microsatsKmerSearch(const njh::progutils::CmdArgs & inputCommands);


	static int genomeKmerCompare(const njh::progutils::CmdArgs & inputCommands);

	static int getNewScanKmerDist(const njh::progutils::CmdArgs & inputCommands);

	static int readingDistanceCheck(const njh::progutils::CmdArgs & inputCommands);
	static int writingDistanceCheck(const njh::progutils::CmdArgs & inputCommands);

	static int writeKmerAccerlation(const njh::progutils::CmdArgs & inputCommands);
	static int writeKmerSimDistanceMatrix(const njh::progutils::CmdArgs & inputCommands);
	static int findUniqKmersBetweenSeqs(const njh::progutils::CmdArgs & inputCommands);


	static int kDistVsNucDist(const njh::progutils::CmdArgs & inputCommands);
	static int scoveViaKmers(const njh::progutils::CmdArgs & inputCommands);

	static int pidVsKmers(const njh::progutils::CmdArgs & inputCommands);

	static int randomSampleKmerCompare(const njh::progutils::CmdArgs & inputCommands);

	static int clostestKmerDist(const njh::progutils::CmdArgs & inputCommands);

	static int scaningKmerDist(const njh::progutils::CmdArgs & inputCommands);


	static int kmerRevVsForDist(const njh::progutils::CmdArgs & inputCommands);

	static int kDist(const njh::progutils::CmdArgs & inputCommands);
	static int getKmerDist(const njh::progutils::CmdArgs & inputCommands);

	static int getKmerDistStatsMultiple(const njh::progutils::CmdArgs & inputCommands);
	static int getKmerDistStats(const njh::progutils::CmdArgs & inputCommands);

	static int profileSharedKmerBlocks(const njh::progutils::CmdArgs & inputCommands);

	static int profileLargeKmerIndex(const njh::progutils::CmdArgs & inputCommands);
	static int profileScanningKmerDist(const njh::progutils::CmdArgs & inputCommands);
	static int profileKmerAccerlation(const njh::progutils::CmdArgs & inputCommands);

	static int findingMinimumKLenForNoRedundantKmers(const njh::progutils::CmdArgs & inputCommands);


	static int generateCountsTable(const njh::progutils::CmdArgs & inputCommands);

	static int getBestKmerDist(const njh::progutils::CmdArgs & inputCommands);




};
} /* namespace njhseq */


