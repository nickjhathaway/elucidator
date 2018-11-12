#pragma once
/*
 * geneExpRunner.h
 *
 *  Created on: October 5, 2015
 *      Author: nickhathaway
 */

#include <njhcpp/progutils.h>
#include <njhseq.h>


namespace njhseq {

class geneExpRunner : public njh::progutils::ProgramRunner {
 public:
	geneExpRunner();

	static int cDNAPosTogDNAPos(const njh::progutils::CmdArgs & inputCommands);

	static int getBedOfAminoAcidPositions(const njh::progutils::CmdArgs & inputCommands);

	static int getBedOfCDnaPositions(const njh::progutils::CmdArgs & inputCommands);

	static int gffRecordIDToGeneInfo(const njh::progutils::CmdArgs & inputCommands);


};
} /* namespace njhseq */


