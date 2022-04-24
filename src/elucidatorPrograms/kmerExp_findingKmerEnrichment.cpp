

#include "kmerExp.hpp"
#include "elucidator/utils/KmerUtils.hpp"
#include "elucidator/objects/dataContainers.h"
#include "elucidator/simulation.h"
#include "elucidator/objects/seqObjects/seqKmers.h"

#include "elucidator/objects/MiscUtility/GenomeSeqSearch.hpp"
#include "elucidator/objects/BioDataObject.h"

#include <njhseq/helpers.h>

namespace njhseq {

int kmerExpRunner::findingKmerEnrichment(const njh::progutils::CmdArgs & inputCommands) {

	std::vector<bfs::path> inputFiles;
	uint32_t kmerLength = 19;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(kmerLength, "--kmerLength", "kmer Length");
	setUp.setOption(inputFiles, "--inputFiles", "input files", true);

	setUp.processDirectoryOutputName("findingKmerEnrichment_TODAY", true);
	setUp.finishSetUp(std::cout);

	setUp.startARunLog(setUp.pars_.directoryName_);


	return 0;
}

} // namespace njhseq