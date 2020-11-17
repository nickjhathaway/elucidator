/*
 * genExp_clusterPointsIntoGroups.cpp
 *
 *  Created on: Nov 14, 2020
 *      Author: nick
 */


#include "genExp.hpp"

#include <njhseq/objects/dataContainers/graphs.h>

namespace njhseq {

int genExpRunner::clusterPointsIntoGroups(const njh::progutils::CmdArgs & inputCommands){
	//in progress
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();

	setUp.finishSetUp(std::cout);

	return 0;
}


}  // namespace njhseq
