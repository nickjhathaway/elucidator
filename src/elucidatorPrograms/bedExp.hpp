#pragma once
/*
 * bedExpRunner.h
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

class bedExpRunner : public njh::progutils::ProgramRunner {
 public:
	bedExpRunner();

	//extracting sequencing using bed file
	static int getFastaWithBed(const njh::progutils::CmdArgs & inputCommands);
	static int getSeqFromTwoBit(const njh::progutils::CmdArgs & inputCommands);

	//extracting sequence information per bed region
	static int getGCContentOfRegion(const njh::progutils::CmdArgs & inputCommands);
	static int getEntropyOfRegion(const njh::progutils::CmdArgs & inputCommands);
	static int getNumOfNsOfRegion(const njh::progutils::CmdArgs & inputCommands);
	static int getNumOfNsEntropyGCContentOfRegion(const njh::progutils::CmdArgs & inputCommands);

	static int getLongestHomopolymerLengthInRegion(const njh::progutils::CmdArgs & inputCommands);

	static int getAverageDistanceToOtherRegions(const njh::progutils::CmdArgs & inputCommands);
	static int getDistanceToClostestRegion(const njh::progutils::CmdArgs & inputCommands);


	//filtering of bed files
	static int getCloseBedRegions(const njh::progutils::CmdArgs & inputCommands);
	static int filterBedRecordsByLength(const njh::progutils::CmdArgs & inputCommands);
	static int getOverlappingBedRegions(const njh::progutils::CmdArgs & inputCommands);
	static int getNonOverlappingBedRegions(const njh::progutils::CmdArgs & inputCommands);
	static int getBestNonOverlapingRegions(const njh::progutils::CmdArgs & inputCommands);
	static int bedUnqiue(const njh::progutils::CmdArgs & inputCommands);
	static int getFirstRegionPerChrom(const njh::progutils::CmdArgs & inputCommands);
	static int getLastRegionPerChrom(const njh::progutils::CmdArgs & inputCommands);
	static int fillInRegionsByBestScore(const njh::progutils::CmdArgs & inputCommands);
	static int bedAddRegionsNonCompletelyInOther(const njh::progutils::CmdArgs & inputCommands);

	static int bedFilterRegionsCompletelyInOther(const njh::progutils::CmdArgs & inputCommands);
	static int bedKeepRegionsCompletelyInOther(const njh::progutils::CmdArgs & inputCommands);
	static int bedGetRegionsEncompassingOthers(const njh::progutils::CmdArgs & inputCommands);
	static int bedFilterRegionsStartingOrEndingInOther(const njh::progutils::CmdArgs & inputCommands);
	static int bedKeepRegionsStartingOrEndingInOther(const njh::progutils::CmdArgs & inputCommands);

	//getting info
	static int getDegreeOfOverlappingBedRegions(const njh::progutils::CmdArgs & inputCommands);
	static int bedGetOverlappingRegionsFromPanels(const njh::progutils::CmdArgs & inputCommands);
	static int bedGetOverlappinBetweenPanels(const njh::progutils::CmdArgs & inputCommands);


	//combing bed files
	static int bedCreateSpanningRegions(const njh::progutils::CmdArgs & inputCommands);
	static int bedRemoveOveringLappingRegions(const njh::progutils::CmdArgs & inputCommands);

	//adding a field to bed file
	static int bedAddSmartIDForPlotting(const njh::progutils::CmdArgs & inputCommands);
	static int bedBinCloseRegions(const njh::progutils::CmdArgs & inputCommands);

	//getting surrounding regions
	static int getUpstreamRegion(const njh::progutils::CmdArgs & inputCommands);
	static int getDownstreamRegion(const njh::progutils::CmdArgs & inputCommands);



	//modifying starts and ends of bed regions
	static int extendBedRegions(const njh::progutils::CmdArgs & inputCommands);
	static int trimBedRegions(const njh::progutils::CmdArgs & inputCommands);
	static int extendUpstreamRegion(const njh::progutils::CmdArgs & inputCommands);
	static int trimUpstreamRegion(const njh::progutils::CmdArgs & inputCommands);
	static int extendDownstreamRegion(const njh::progutils::CmdArgs & inputCommands);
	static int trimDownstreamRegion(const njh::progutils::CmdArgs & inputCommands);
	static int extendToEndOfChrom(const njh::progutils::CmdArgs & inputCommands);
	static int extendToStartOfChrom(const njh::progutils::CmdArgs & inputCommands);
	static int centerBedRegionWithFixSize(const njh::progutils::CmdArgs & inputCommands);
	static int reverseComplementRegion(const njh::progutils::CmdArgs & inputCommands);

	//modifying aspects of the bed regions
	static int bedToggleStrand(const njh::progutils::CmdArgs & inputCommands);
	static int bed3ToBed6(const njh::progutils::CmdArgs & inputCommands);
	static int bedRenameRepeatUids(const njh::progutils::CmdArgs & inputCommands);
	static int bedRenameWithCoords(const njh::progutils::CmdArgs & inputCommands);
	static int bedRenameChromosomes(const njh::progutils::CmdArgs & inputCommands);
	static int bedChangeScoreToLength(const njh::progutils::CmdArgs & inputCommands);
	static int bedCoordSort(const njh::progutils::CmdArgs & inputCommands);
	static int reorientBasedOnSingleReadsOrientationCounts(const njh::progutils::CmdArgs & inputCommands);

	//creating multiple bed files from one bed file
	static int splitBedFile(const njh::progutils::CmdArgs & inputCommands);
	static int separateOutRecordsInBedFile(const njh::progutils::CmdArgs & inputCommands);

	//extraction
	static int extractBedRecordsWithName(const njh::progutils::CmdArgs & inputCommands);

	//getting new regions
	static int getInterveningRegions(const njh::progutils::CmdArgs & inputCommands);
	static int createWindowsInbetweenRegions(const njh::progutils::CmdArgs & inputCommands);
	static int createWindowsInRegions(const njh::progutils::CmdArgs & inputCommands);

	static int getIntersectionBetweenTwoBedFiles(const njh::progutils::CmdArgs & inputCommands);
	static int removeSubRegionsFromBedFile(const njh::progutils::CmdArgs & inputCommands);

	//coverage
	static int roughSmoothingForBedCoverage(const njh::progutils::CmdArgs & inputCommands);

	//creating bed files
	static int fastaToBed(const njh::progutils::CmdArgs & inputCommands);
	static int differentSubRegionCombos(const njh::progutils::CmdArgs & inputCommands);

};
} /* namespace njhseq */


