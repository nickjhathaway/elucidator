//
// Created by Nicholas H	athaway on 3/31/25.
//

#include "readSimulatorRunner.hpp"
#include "elucidator/simulation.h"
#include <SeekDeep/objects/TarAmpSetupUtils/PrimersAndMids.hpp>

namespace njhseq {



int readSimulatorRunner::simulateSpecificSamples(const njh::progutils::CmdArgs & inputCommands) {
	readSimulatorSetUp setUp(inputCommands);
	uint32_t numThreads = 2;
	uint32_t pcrNumThreads = 2;
	bool singleEnd = false;
	bool nonGz = false;
	bfs::path idFile = "";
	bfs::path referenceFile = "";
	bfs::path illuminaProfileDir = "";
	uint32_t defaultPcrRounds = 30;
	uint32_t initialPcrRounds = 10;
	std::map<uint32_t, uint32_t> initialPcrRoundsMap;
	bfs::path initialPcrRoundsMapFnp;
	long double errorRate = 3.5e-06;
	double pcrEfficiency = 0.85;
	bool keepPCRSeqs = false;
	uint32_t chimeraBasesIn = 5;
	uint32_t templateCap = 500000000;
	bool noChimeras = false;
	double finalReadAmountSDFrac = 0.1;

	uint32_t pairedEndLength = std::numeric_limits<uint32_t>::max();
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(pairedEndLength, "--pairedEndLength", "Paired End Length");
	setUp.setOption(finalReadAmountSDFrac, "--finalReadAmountSDFrac", "final Read Amount SD Frac", njh::progutils::ProgramSetUp::CheckCase::GREATERZERO);
	setUp.setOption(noChimeras, "--noChimeras", "Don't simulate chimeras");
	setUp.setOption(templateCap, "--templateCap", "Template Cap");
	setUp.setOption(chimeraBasesIn, "--chimeraBasesIn", "The number of bases needed for a template to lay down");
	setUp.setOption(keepPCRSeqs, "--keepPCRSeqs", "Keep PCR Seqs");
	setUp.setOption(errorRate, "--errorRate", "Polymerase Error Rate");
	setUp.setOption(pcrEfficiency, "--pcrEfficiency", "PCR Efficiency, between 0-1, chance a product gets amplified");
	setUp.setOption(defaultPcrRounds, "--pcrRounds", "Number of PCR rounds");
	setUp.setOption(initialPcrRounds, "--initialPcrRounds", "Number of Initial rounds of PCR before sampling");
	setUp.setOption(initialPcrRoundsMapFnp, "--initialPcrRoundsTable", "Number of Initial rounds of PCR before sampling per starting template amount, columns 1)template, 2) rounds");
	setUp.setOption(numThreads, "--numThreads", "Number of Threads to Use");
	setUp.setOption(pcrNumThreads, "--pcrNumThreads", "Number of Threads to Use for PCR sim");
	setUp.setOption(singleEnd, "--singleEnd", "Single End");
	setUp.setOption(nonGz, "--nonGz", "do not compress the output fastqs");
	setUp.setOption(illuminaProfileDir, "--illuminaProfileDir", "Illumina Profile Dir", true);
	setUp.processDirectoryOutputName("simulateSpecificSamples_TODAY", true);
	setUp.finishSetUp(std::cout);

	setUp.startARunLog(setUp.pars_.directoryName_);
	initialPcrRoundsMap[1] = initialPcrRounds;
	if(bfs::exists(initialPcrRoundsMapFnp)){
		initialPcrRoundsMap.clear();
		table initialPcrRoundsMapTab(initialPcrRoundsMapFnp, "\t", true);
		initialPcrRoundsMapTab.checkForColumnsThrow(VecStr{"template", "rounds"}, __PRETTY_FUNCTION__);
		for(const auto & row : initialPcrRoundsMapTab){
			initialPcrRoundsMap[njh::StrToNumConverter::stoToNum<uint32_t>(row[initialPcrRoundsMapTab.getColPos("template")])] =njh::StrToNumConverter::stoToNum<uint32_t>(row[initialPcrRoundsMapTab.getColPos("rounds")]);
		}
	}
	uint64_t intErrorRate = errorRate * std::numeric_limits<uint64_t>::max();

	setUp.writeParametersFile(setUp.pars_.directoryName_ + "parameters.tab.txt", false, true);
	std::unordered_map<std::string, std::shared_ptr<seqInfo>> refSeqs;
	std::unordered_map<std::string, std::vector<std::shared_ptr<seqInfo>>> refSeqsByPrimerName;
	std::unordered_map<std::string, std::set<std::string>> refSeqsNamesByPrimerName;

	return 0;
}

} // namespace njhseq

