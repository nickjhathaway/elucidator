
// Created on 2015/01/13
// main.cpp



#include "elucidatorPrograms.h"
#include <njhcpp/progutils/oneRing.hpp>

#include <seqServerPrograms.h>
#include <SeekDeepPrograms.h>
#include <TwoBitPrograms.h>
#include <njhseq/ProgramRunners.h>

namespace njhseq {
int testFunc(const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	setUp.finishSetUp(std::cout);
	std::cout << "This is a test!" << std::endl;
	return 0;
}
class elucidatorRunner: public njh::progutils::OneRing {
public:
	elucidatorRunner();
};
//,createReadConsensusRandomPicking  BamExtractPathawaysFromRegion
elucidatorRunner::elucidatorRunner() :
		njh::progutils::OneRing(
				{ addRing<pacbioExpRunner>(),      addRing<genExpRunner>(),
					addRing<kmerExpRunner>(),        addRing<clusteringExpRunner>(),
					addRing<fileFormatExpRunner>(),
					addRing<bamExpRunner>(),         addRing<bamCountingExpRunner>(),
					addRing<programWrapperRunner>(),
					addRing<seqUtilsTrimRunner>(),
					addRing<seqUtilsModRunner>(),
					addRing<gffExpRunner>(),   			 addRing<bedExpRunner>(),
					addRing<geneExpRunner>(),        addRing<SeqServerRunner>(),
					addRing<SeekDeepRunner>(),       addRing<TwoBit::TwoBitRunner>(),
					addRing<repelinRunner>(),        addRing<ampliconAnalysisRunner>(),
					addRing<readSimulatorRunner>(),  addRing<graphicsUtilsRunner>(),
					addRing<miscRunner>(),           addRing<printInfoRunner>(),
					addRing<seqUtilsRunner>(),       addRing<seqUtilsConverterRunner>(),
					addRing<seqUtilsExtractRunner>(),addRing<seqUtilsInfoRunner>(),
					addRing<simulatorRunner>(),      addRing<parsingFileExpRunner>(),
					addRing<benchMarkingRunner>(),   addRing<jsonExpRunner>(),
					addRing<ManipulateTableRunner>(),
					addRing<metaExpRunner>(),
					addRing<seqSearchingRunner>(),
					addRing<pairProcessingRunner>(),
				},//
				{ addFunc("testFunc", testFunc, false)
				}, "elucidator", "1", "0", "0-dev") {
}

} //namepsace njhseq



int main(int argc, char* argv[]) {
	try {
		njhseq::elucidatorRunner runner;
		return runner.run(argc, argv);
	} catch (std::exception & e) {
		std::cerr << e.what() << std::endl;
		return 1;
	}
}
