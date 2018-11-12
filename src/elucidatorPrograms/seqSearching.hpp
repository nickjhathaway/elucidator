#pragma once
/*
 * seqSearching.hpp
 *
 *  Created on: Nov 12, 2018
 *      Author: nick
 */


#include <njhcpp/progutils.h>
#include <njhseq.h>

namespace njhseq {



class seqSearchingRunner : public njh::progutils::ProgramRunner {
 public:
	seqSearchingRunner();


	static int chopAndMap(const njh::progutils::CmdArgs & inputCommands);
	static int chopAndMapAndRefine(const njh::progutils::CmdArgs & inputCommands);



};

}  // namespace njhseq


