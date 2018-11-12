#pragma once
/*
 * parsingFileExpRunner.h
 *
 *  Created on: July 4, 2017
 *      Author: nickhathaway
 */

#include <njhcpp/progutils.h>
#include <njhseq.h>


namespace njhseq {

class parsingFileExpRunner : public njh::progutils::ProgramRunner {
 public:
	parsingFileExpRunner();

	static int parsingBioSampleSetNCBIJson(const njh::progutils::CmdArgs & inputCommands);
	static int getStrainInfoTableBioSampleSetNCBIJson(const njh::progutils::CmdArgs & inputCommands);

	static int getAttributeLevelsBioSampleSetNCBIJson(const njh::progutils::CmdArgs & inputCommands);
	static int parseBlastpHitsTab(const njh::progutils::CmdArgs & inputCommands);

	static int parseSTOCKHOLM(const njh::progutils::CmdArgs & inputCommands);
	static int parseSTOCKHOLMToFasta(const njh::progutils::CmdArgs & inputCommands);

	static int parsehmmerDomainHitTab(const njh::progutils::CmdArgs & inputCommands);

};
} /* namespace njhseq */


