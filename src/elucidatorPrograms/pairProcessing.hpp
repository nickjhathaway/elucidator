#pragma once
/*
 * pairProcessing.hpp
 *
 *  Created on: Nov 12, 2018
 *      Author: nick
 */
#include <njhcpp/progutils.h>
#include <njhseq.h>

namespace njhseq {

class pairProcessingRunner : public njh::progutils::ProgramRunner {
 public:
	pairProcessingRunner();



	static int detectPossiblePrimers(const njh::progutils::CmdArgs & inputCommands);
	static int StitchPairedReads(const njh::progutils::CmdArgs & inputCommands);


};

} //  namespace njhseq








