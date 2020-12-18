#pragma once

/*
 * genomeExp.hpp
 *
 *  Created on: Dec 17, 2020
 *      Author: nick
 */







#include <njhcpp/progutils.h>
#include <njhseq/common.h>
#include <njhseq/programUtils/seqSetUp.hpp>

namespace njhseq {

class genomeExpRunner : public njh::progutils::ProgramRunner {
 public:
	genomeExpRunner();
	static int reorientToRefGenome(const njh::progutils::CmdArgs & inputCommands);

};
} /* namespace njhseq */



