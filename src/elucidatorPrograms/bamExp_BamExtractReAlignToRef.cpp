/*
 * bamExp_BamExtractReAlignToRef.cpp
 *
 *  Created on: Aug 16, 2019
 *      Author: nicholashathaway
 */



#include "bamExp.hpp"
#include <TwoBit.h>

#include "elucidator/BamToolsUtils.h"
#include <njhseq/objects/BioDataObject.h>
#include "elucidator/BioRecordsUtils/BedUtility.hpp"



namespace njhseq {





int bamExpRunner::BamExtractReAlignToRef(const njh::progutils::CmdArgs & inputCommands){
	ReAlignedSeq::genRealignmentPars ralnPars;
	std::string twoBitFnp = "";

	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(twoBitFnp, "--2bit", "2bit file of the genome", true);
	setUp.setOption(ralnPars.extendAmount, "--extendAmount", "Extend Amount", true);

	setUp.processReadInNames({"--bam"}, true);
	setUp.finishSetUp(std::cout);


	BamTools::BamReader bReader;
	bReader.Open(setUp.pars_.ioOptions_.firstName_.string());
	checkBamOpenThrow(bReader, setUp.pars_.ioOptions_.firstName_.string());

	auto refData = bReader.GetReferenceData();

	TwoBit::TwoBitFile tReader(twoBitFnp);

	BamTools::BamAlignment bAln;

	aligner alignerObj(500, gapScoringParameters(5,1,0,0,0,0));
	auto chromLengths = tReader.getSeqLens();

	while(bReader.GetNextAlignment(bAln)){
		if(bAln.IsMapped() && bAln.IsPrimaryAlignment()){
			auto realignment = ReAlignedSeq::genRealignment(bAln, refData, alignerObj, chromLengths, tReader,ralnPars);
		}
	}

	return 0;
}

}  // namespace njhseq











