#pragma once
/*
 * bamExpRunner.h
 *
 *  Created on: Feb 4, 2015
 *      Author: nickhathaway
 */

#include <njhcpp/progutils.h>
#include <njhseq.h>

namespace njhseq {

class bamExpRunner : public njh::progutils::ProgramRunner {
 public:
	bamExpRunner();

	static int isBamSorted(const njh::progutils::CmdArgs & inputCommands);

	static int determineRegion(const njh::progutils::CmdArgs & inputCommands);
	static int determineRegionLastz(const njh::progutils::CmdArgs & inputCommands);

	static int bamReadProfiling(const njh::progutils::CmdArgs & inputCommands);

	static int bamToFastq(const njh::progutils::CmdArgs & inputCommands);
	static int getInsertSizesStats(const njh::progutils::CmdArgs & inputCommands);

	static int BamExtractReadsFromRegion(const njh::progutils::CmdArgs & inputCommands);

	static int BamGetFileIndexPositionOfName(const njh::progutils::CmdArgs & inputCommands);
	static int BamFindDifferenceInUnmappedFileIndexPosition(const njh::progutils::CmdArgs & inputCommands);





	static int getMateMapStatus(const njh::progutils::CmdArgs & inputCommands);
	static int getInsertSizeChanges(const njh::progutils::CmdArgs & inputCommands);
	static int getMapQualityCounts(const njh::progutils::CmdArgs & inputCommands);


	static int multiBamCoverageFinder(const njh::progutils::CmdArgs & inputCommands);
	static int bamMultiPairStats(const njh::progutils::CmdArgs & inputCommands);
	static int bamMulticov(const njh::progutils::CmdArgs & inputCommands);
	static int bamMulticovBases(const njh::progutils::CmdArgs & inputCommands);

	static int outputSoftClipCounts(const njh::progutils::CmdArgs & inputCommands);

	static int printBamRefIds(const njh::progutils::CmdArgs & inputCommands);
	static int bamToBed(const njh::progutils::CmdArgs & inputCommands);


	static int getBestBedRegionFromBam(const njh::progutils::CmdArgs & inputCommands);
	static int mergeBedRegionFromBam(const njh::progutils::CmdArgs & inputCommands);


	static int getUnmappedAlnsFromBam(const njh::progutils::CmdArgs & inputCommands);



	static int countUnMappedMateStatus(const njh::progutils::CmdArgs & inputCommands);

	static int genBedFromMappedGenome(const njh::progutils::CmdArgs & inputCommands);


	static int appendReadGroupToName(const njh::progutils::CmdArgs & inputCommands);
	static int BamRenameRefHeader(const njh::progutils::CmdArgs & inputCommands);


	static int GetBamLocationStats(const njh::progutils::CmdArgs & inputCommands);


	static int GetKmerCoverageForRegions(const njh::progutils::CmdArgs & inputCommands);
	static int testingBamToFastq(const njh::progutils::CmdArgs & inputCommands);


	static int printHeaderRefIndexes(const njh::progutils::CmdArgs & inputCommands);
	static int refineBedRegionFromBam(const njh::progutils::CmdArgs & inputCommands);


	static int BamFilterByChroms(const njh::progutils::CmdArgs & inputCommands);

};

} /* namespace njhseq */


