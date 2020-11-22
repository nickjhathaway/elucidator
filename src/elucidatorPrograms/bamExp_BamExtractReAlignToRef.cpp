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
	OutOptions outOpts;
	gapScoringParameters gapPars(5,1,0,0,0,0);
	int32_t match = 2;
	int32_t mismatch = -2;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(twoBitFnp, "--2bit", "2bit file of the genome", true);
	setUp.setOption(ralnPars.extendAmount, "--extendAmount", "Extend Amount");
	setUp.setOption(match, "--match", "match score for alignment");
	setUp.setOption(mismatch, "--mismatch", "mismatch score for alignment");
	setUp.setOption(gapPars.gapOpen_, "--gapOpen", "gap open penalty for alignment");
	setUp.setOption(gapPars.gapExtend_, "--gapExtend", "gap extend penalty for alignment");

	setUp.processWritingOptions(outOpts);
	setUp.processReadInNames({"--bam"}, true);
	setUp.finishSetUp(std::cout);


	BamTools::BamReader bReader;
	bReader.Open(setUp.pars_.ioOptions_.firstName_.string());
	checkBamOpenThrow(bReader, setUp.pars_.ioOptions_.firstName_.string());

	auto refData = bReader.GetReferenceData();

	TwoBit::TwoBitFile tReader(twoBitFnp);

	BamTools::BamAlignment bAln;

//	aligner alignerObj(500, gapPars, substituteMatrix(match,mismatch));
//
//	{
//		aligner alignerObj(500, gapScoringParameters(5,1,0,0,0,0));
//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
//		std::cout << "alignerObj.parts_.scoring_.mat_['A']['A']" << alignerObj.parts_.scoring_.mat_['A']['A'] << std::endl;
//		std::cout << "alignerObj.parts_.scoring_.mat_['A']['G']" << alignerObj.parts_.scoring_.mat_['A']['G'] << std::endl;
//
//	}
	aligner alignerObj(500, gapPars, substituteMatrix(match,mismatch));
//	std::cout << __FILE__ << " " << __LINE__ << std::endl;
//	std::cout << "alignerObj.parts_.scoring_.mat_['A']['A']" << alignerObj.parts_.scoring_.mat_['A']['A'] << std::endl;
//	std::cout << "alignerObj.parts_.scoring_.mat_['A']['G']" << alignerObj.parts_.scoring_.mat_['A']['G'] << std::endl;
//


	auto chromLengths = tReader.getSeqLens();
	auto outSeqOpts = SeqIOOptions::genFastqOut(outOpts.outName());
	outSeqOpts.out_.transferOverwriteOpts(outOpts);
	SeqOutput writer(outSeqOpts);
	writer.openOut();
	while(bReader.GetNextAlignment(bAln)){
		if(bAln.IsMapped() && bAln.IsPrimaryAlignment()){
			auto realignment = ReAlignedSeq::genRealignment(bAln, refData, alignerObj, chromLengths, tReader, ralnPars);
			writer.write(realignment.alnRefSeq_);
			writer.write(realignment.querySeq_);
		}
	}
	return 0;
}

}  // namespace njhseq











