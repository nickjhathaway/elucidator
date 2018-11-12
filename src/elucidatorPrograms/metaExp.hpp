#pragma once
/*
 * metaExp.hpp
 *
 *  Created on: Mar 14, 2018
 *      Author: nick
 */



#include <njhcpp/progutils.h>
#include <njhseq.h>

namespace njhseq {

class metaExpRunner : public njh::progutils::ProgramRunner {
 public:
	metaExpRunner();

	static int addMetaBySampleName(const njh::progutils::CmdArgs & inputCommands);
	static int addMetaByMetaField(const njh::progutils::CmdArgs & inputCommands);
	static int excludeSeqsFileWithNumericMetaCutOff(const njh::progutils::CmdArgs & inputCommands);
	static int excludeSeqsFileWithMatchingMeta(const njh::progutils::CmdArgs & inputCommands);
	static int selectMetaFieldsToKeep(const njh::progutils::CmdArgs & inputCommands);
	static int splitSeqFileWithMeta(const njh::progutils::CmdArgs & inputCommands);
	static int splitSeqFileWithExternalMeta(const njh::progutils::CmdArgs & inputCommands);
  static int createTableFromSeqs(const njh::progutils::CmdArgs & inputCommands);
  static int printMetaFieldsFromSeqs(const njh::progutils::CmdArgs & inputCommands);

};

} //  namespace njhseq



