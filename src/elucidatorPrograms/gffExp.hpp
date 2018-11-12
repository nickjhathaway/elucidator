#pragma once
/*
 * gffExpRunner.h
 *
 *  Created on: October 5, 2015
 *      Author: nickhathaway
 */

#include <njhcpp/progutils.h>
#include <njhseq.h>


namespace njhseq {

class gffExpRunner : public njh::progutils::ProgramRunner {
 public:
	gffExpRunner();
	static int gffFeatureCount(const njh::progutils::CmdArgs & inputCommands);
	static int gffDescriptionsCount(const njh::progutils::CmdArgs & inputCommands);
	static int gffAttributesCount(const njh::progutils::CmdArgs & inputCommands);
	static int gffCountAttribute(const njh::progutils::CmdArgs & inputCommands);
	static int extractGffFeature(const njh::progutils::CmdArgs & inputCommands);
	static int gffToBed(const njh::progutils::CmdArgs & inputCommands);
	static int gffToBedByAttribute(const njh::progutils::CmdArgs & inputCommands);
	static int gffToBedByDescription(const njh::progutils::CmdArgs & inputCommands);
	static int gffToBedByChrom(const njh::progutils::CmdArgs & inputCommands);
	static int gffToBedByBedLoc(const njh::progutils::CmdArgs & inputCommands);
	static int gffToBedByName(const njh::progutils::CmdArgs & inputCommands);

	static int gffToJsonByID(const njh::progutils::CmdArgs & inputCommands);


	static int extractGffRecordWithChildren(const njh::progutils::CmdArgs & inputCommands);
	static int gffToBedByAttributeIncludeExonInfo(const njh::progutils::CmdArgs & inputCommands);

	static int bedGetIntersectingGenesInGff(const njh::progutils::CmdArgs & inputCommands);
	static int bedGetIntersectingRecordsInGff(const njh::progutils::CmdArgs & inputCommands);

	static int reorientBedToIntersectingGeneInGff(const njh::progutils::CmdArgs & inputCommands);
	static int setBedPositionsToIntersectingGeneInGff(const njh::progutils::CmdArgs & inputCommands);

	static int gffGetNumOfTranscriptsForGenes(const njh::progutils::CmdArgs & inputCommands);
	static int gffPrintIds(const njh::progutils::CmdArgs & inputCommands);

	static int testingGffReading(const njh::progutils::CmdArgs & inputCommands);
	static int removeFastaFromGffFile(const njh::progutils::CmdArgs & inputCommands);


	static int aaPositionsToBed(const njh::progutils::CmdArgs & inputCommands);

};
} /* namespace njhseq */


