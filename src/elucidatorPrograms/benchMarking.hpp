#pragma once
/*
 * benchMarkings.hpp
 *
 *  Created on: July 5, 2017
 *      Author: nick
 */



#include <njhcpp/progutils.h>
#include <njhseq.h>

namespace njhseq {

class benchMarkingRunner : public njh::progutils::ProgramRunner {
 public:
	benchMarkingRunner();


	static int benchGzWritingOneRead(const njh::progutils::CmdArgs & inputCommands);
	static int benchGzWritingSetChunks(const njh::progutils::CmdArgs & inputCommands);
	static int benchGzWritingRdBuf(const njh::progutils::CmdArgs & inputCommands);


};

} /* namespace njhseq */





