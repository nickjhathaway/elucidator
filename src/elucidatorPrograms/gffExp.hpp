#pragma once
/*
 * gffExpRunner.h
 *
 *  Created on: October 5, 2015
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
#include <njhseq/common.h>
#include <njhseq/programUtils/seqSetUp.hpp>


namespace njhseq {

class gffExpRunner : public njh::progutils::ProgramRunner {
 public:
	gffExpRunner();
	static int gffFeatureCount(const njh::progutils::CmdArgs & inputCommands);
	static int gffDescriptionsCount(const njh::progutils::CmdArgs & inputCommands);
	static int gffAttributesCount(const njh::progutils::CmdArgs & inputCommands);
	static int gffCountAttribute(const njh::progutils::CmdArgs & inputCommands);
	static int gffToBed(const njh::progutils::CmdArgs & inputCommands);
	static int gffToBedByFeature(const njh::progutils::CmdArgs & inputCommands);
	static int gffToBedByAttribute(const njh::progutils::CmdArgs & inputCommands);
	static int gffToBedByDescription(const njh::progutils::CmdArgs & inputCommands);
	static int gffToBedByChrom(const njh::progutils::CmdArgs & inputCommands);
	static int gffToBedByBedLoc(const njh::progutils::CmdArgs & inputCommands);
	static int gffToBedByName(const njh::progutils::CmdArgs & inputCommands);
	static int gffToJsonByID(const njh::progutils::CmdArgs & inputCommands);

	static int extractGffFeature(const njh::progutils::CmdArgs & inputCommands);
	static int extractGffChrom(const njh::progutils::CmdArgs & inputCommands);



	static int gffSortInefficient(const njh::progutils::CmdArgs & inputCommands);
	static int roughGffConversionToOther(const njh::progutils::CmdArgs & inputCommands);


	static int extractGffRecordWithChildren(const njh::progutils::CmdArgs & inputCommands);
	static int gffToBedByAttributeIncludeExonInfo(const njh::progutils::CmdArgs & inputCommands);

	static int bedGetRegionsCompletelyInGenesInGff(const njh::progutils::CmdArgs & inputCommands);
	static int bedGetIntersectingGenesInGff(const njh::progutils::CmdArgs & inputCommands);
	static int bedGetIntersectingRecordsInGff(const njh::progutils::CmdArgs & inputCommands);

	static int reorientBedToIntersectingGeneInGff(const njh::progutils::CmdArgs & inputCommands);
	static int setBedPositionsToIntersectingGeneInGff(const njh::progutils::CmdArgs & inputCommands);

	static int gffTranscriptIDForGeneIDs(const njh::progutils::CmdArgs & inputCommands);
	static int gffGetNumOfTranscriptsForGenes(const njh::progutils::CmdArgs & inputCommands);
	static int gffPrintIds(const njh::progutils::CmdArgs & inputCommands);

	static int testingGffReading(const njh::progutils::CmdArgs & inputCommands);
	static int removeFastaFromGffFile(const njh::progutils::CmdArgs & inputCommands);


	static int aaPositionsToBed(const njh::progutils::CmdArgs & inputCommands);

	static int appendGff(const njh::progutils::CmdArgs & inputCommands);

	static int revCompGff(const njh::progutils::CmdArgs & inputCommands);




};
} /* namespace njhseq */


