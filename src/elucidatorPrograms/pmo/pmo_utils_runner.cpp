//
// Created by Nicholas Hathaway on 8/20/25.
//

#include "pmo_utils_runner.hpp"


#include <njhseq/IO/SeqIO/SeqIO.hpp>
#include <njhseq/objects/BioDataObject/pmo.h>
#include <nlohmann/json.hpp>

namespace njhseq {
PMOUtilsRunner::PMOUtilsRunner()
        : njh::progutils::ProgramRunner(
        {
                addFunc("read_pmo", read_pmo, false),

        },//
        "PMOUtils") {}


int PMOUtilsRunner::read_pmo(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path pmo_fnp;
	OutOptions out_options;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(pmo_fnp, "--pmo", "pmo file to read", true);
	setUp.processWritingOptions(out_options);
	setUp.finishSetUp(std::cout);

	InputStream in(pmo_fnp);
	auto pmo_obj= pmo::PortableMicrohaplotypeObject::from_json(nlohmann::json::parse(in));
	pmo_obj.validate();
	OutputStream out(out_options);
	out << pmo_obj.to_json();

	return 0;
}


} //namespace njhseq


