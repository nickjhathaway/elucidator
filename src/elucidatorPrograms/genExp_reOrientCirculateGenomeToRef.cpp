/*
 * genExp_reOrientCirculateGenomeToRef.cpp
 *
 *  Created on: Jun 27, 2020
 *      Author: nick
 */




#include "genExp.hpp"
#include <njhseq/readVectorManipulation/readVectorHelpers/readVecTrimmer.hpp>
#include <njhseq/IO/SeqIO/SeqIO.hpp>


namespace njhseq {





int genExpRunner::reOrientCirculateGenomeToRef(const njh::progutils::CmdArgs & inputCommands){

	readVecTrimmer::trimCircularGenomeToRefPars trimPars;

	auto outSeqOpts = SeqIOOptions::genFastaOut("");
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processWritingOptions(outSeqOpts.out_);
	setUp.setOption(trimPars.kmerLength_, "--kmerLength", "kmer Length");
	setUp.processReadInNames(true);
	setUp.processSeq(trimPars.refSeq_, "--ref",   "reference circulate genome", true);
	setUp.setOption(trimPars.doNotReOrientDirection_, "--doNotReOrientDirection", "do Not Re Orient Direction");
	setUp.setOption(trimPars.preferHeader_, "--preferHeader", "prefer Header");
	setUp.setOption(trimPars.padding_, "--padding", "padding");
	setUp.setOption(trimPars.mark_, "--mark", "mark sequence trim positions");

	setUp.pars_.gapLeft_ = "0,0";
	setUp.pars_.gapRight_ = "0,0";
	setUp.pars_.gapInfo_ = gapScoringParameters(5,1,0,0,0,0);
	setUp.pars_.colOpts_.iTOpts_.weighHomopolyer_ = false;
	setUp.processAlignerDefualts();
	setUp.processDirectoryOutputName(false);
  if(setUp.pars_.verbose_){
    std::cout << "Gap open: " << setUp.pars_.gapInfo_.gapOpen_ << std::endl;
    std::cout << "Gap extend: " << setUp.pars_.gapInfo_.gapExtend_ << std::endl;
    std::cout << "Gap Left Query open: " << setUp.pars_.gapInfo_.gapLeftQueryOpen_ << std::endl;
    std::cout << "Gap Left Query extend: " << setUp.pars_.gapInfo_.gapLeftQueryExtend_ << std::endl;
    std::cout << "Gap Right Query open: " << setUp.pars_.gapInfo_.gapRightQueryOpen_ << std::endl;
    std::cout << "Gap Right Query extend: " << setUp.pars_.gapInfo_.gapRightQueryExtend_ << std::endl;
    std::cout << "Gap Left Ref open: " << setUp.pars_.gapInfo_.gapLeftRefOpen_ << std::endl;
    std::cout << "Gap Left Ref extend: " << setUp.pars_.gapInfo_.gapLeftRefExtend_ << std::endl;
    std::cout << "Gap Right Ref open: " << setUp.pars_.gapInfo_.gapRightRefOpen_ << std::endl;
    std::cout << "Gap Right Ref extend: " << setUp.pars_.gapInfo_.gapRightRefExtend_ << std::endl;
  }
	setUp.finishSetUp(std::cout);

	SeqOutput writer(outSeqOpts);
	writer.openOut();

	uint64_t maxSize = 0;
	readVec::getMaxLength(trimPars.refSeq_, maxSize);
	{
		SeqInput reader(setUp.pars_.ioOptions_);
		seqInfo seq;
		reader.openIn();
		while(reader.readNextRead(seq)){
			readVec::getMaxLength(seq, maxSize);
		}
	}
	KmerMaps emptyMaps;
	aligner alignerObj(maxSize, setUp.pars_.gapInfo_,
			setUp.pars_.scoring_, emptyMaps, setUp.pars_.qScorePars_,
			setUp.pars_.colOpts_.alignOpts_.countEndGaps_,
			setUp.pars_.colOpts_.iTOpts_.weighHomopolyer_);
	alignerObj.processAlnInfoInput(setUp.pars_.alnInfoDirName_, setUp.pars_.verbose_);

	std::vector<seqInfo> outSeqs;
	SeqInput reader(setUp.pars_.ioOptions_);
	seqInfo seq;
	reader.openIn();
	while(reader.readNextRead(seq)){
		auto res = readVecTrimmer::trimCircularGenomeToRef(seq, trimPars, alignerObj);
		addOtherVec(outSeqs, res);
	}


	writer.write(outSeqs);
	alignerObj.processAlnInfoOutput(setUp.pars_.outAlnInfoDirName_, setUp.pars_.verbose_);

	return 0;
}



}  // namespace njhseq
