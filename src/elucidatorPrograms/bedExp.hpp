#pragma once
/*
 * bedExpRunner.h
 *
 *  Created on: October 5, 2015
 *      Author: nickhathaway
 */

#include <njhcpp/progutils.h>
#include <njhseq.h>


namespace njhseq {

class bedExpRunner : public njh::progutils::ProgramRunner {
 public:
	bedExpRunner();

	//extracting sequencing using bed file
	static int getFastaWithBed(const njh::progutils::CmdArgs & inputCommands);

	//filtering of bed files
	static int getCloseBedRegions(const njh::progutils::CmdArgs & inputCommands);
	static int filterBedRecordsByLength(const njh::progutils::CmdArgs & inputCommands);
	static int getOverlappingBedRegions(const njh::progutils::CmdArgs & inputCommands);
	static int getNonOverlappingBedRegions(const njh::progutils::CmdArgs & inputCommands);
	static int getBestNonOverlapingRegions(const njh::progutils::CmdArgs & inputCommands);
	static int bedUnqiue(const njh::progutils::CmdArgs & inputCommands);

	//combing bed files
	static int bedCreateSpanningRegions(const njh::progutils::CmdArgs & inputCommands);

	//adding a field to bed file
	static int bedAddSmartIDForPlotting(const njh::progutils::CmdArgs & inputCommands);

	//getting surrounding regions
	static int getUpstreamRegion(const njh::progutils::CmdArgs & inputCommands);
	static int getDownstreamRegion(const njh::progutils::CmdArgs & inputCommands);

	//modifying starts and ends of bed regions
	static int extendBedRegions(const njh::progutils::CmdArgs & inputCommands);
	static int extendUpstreamRegion(const njh::progutils::CmdArgs & inputCommands);
	static int extendDownstreamRegion(const njh::progutils::CmdArgs & inputCommands);
	static int extendToEndOfChrom(const njh::progutils::CmdArgs & inputCommands);
	static int extendToStartOfChrom(const njh::progutils::CmdArgs & inputCommands);

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



};
} /* namespace njhseq */


