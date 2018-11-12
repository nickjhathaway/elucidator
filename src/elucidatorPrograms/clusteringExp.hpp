#pragma once
/*
 * clusteringExpRunner.h
 *
 *  Created on: March 5, 2015
 *      Author: nickhathaway
 */

#include <njhcpp/progutils.h>
#include <njhseq.h>

namespace njhseq {

class clusteringExpRunner : public njh::progutils::ProgramRunner {
 public:
	clusteringExpRunner();
	static int clusteringPairedEndReads(const njh::progutils::CmdArgs & inputCommands);

};
} /* namespace njhseq */


