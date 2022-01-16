/*
 * seqUtilsTrimRunner_trimToRefWithGlobalAlignmentNonOverlappingRegions.cpp
 *
 *  Created on: Jan 15, 2022
 *      Author: nick
 */





#include "seqUtilsTrimRunner.hpp"
#include <njhseq/IO/SeqIO.h>
#include <njhseq/readVectorManipulation/readVectorHelpers/readVecTrimmer.hpp>
#include <njhseq/objects/Meta/MetaDataInName.hpp>

namespace njhseq {





int seqUtilsTrimRunner::trimToRefWithGlobalAlignmentNonOverlappingRegions(const njh::progutils::CmdArgs & inputCommands){
	//bool mark = false;
	bool getRevComp = false;
	seqUtilsTrimSetUp setUp(inputCommands);
	setUp.description_ = "Trim the left side of sequences to the mean alignment site of reference sequences";
	FullTrimReadsPars pars;
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processRefFilename(true);
	setUp.pars_.gapLeft_ = "0,0";
	setUp.pars_.gapRight_ = "0,0";
	setUp.pars_.gap_ = "5,1";
	setUp.setOption(getRevComp, "--getRevComp", "get Rev Comp");
	setUp.processAlignerDefualts();

	setUp.processIoOptions(setUp.readInFormatsAvailable_);
	setUp.finishSetUp(std::cout);

	SeqIO reader(setUp.pars_.ioOptions_	);
	reader.openIn();
	reader.openOut();
	uint64_t maxLen = 0;
	seqInfo seq;
	while(reader.readNextRead(seq)){
		readVec::getMaxLength(seq, maxLen);
	}

	auto refSeqs = SeqInput::getSeqVec<seqInfo>(setUp.pars_.refIoOptions_, maxLen);
	std::vector<seqInfo> revComp_refSeqs;

	std::vector<kmerInfo> refSeqsKInfos;
	std::vector<kmerInfo> revComp_refSeqsKInfos;

	if(refSeqs.size() != 1){
		for(const auto & rSeq : refSeqs){
			refSeqsKInfos.emplace_back(rSeq.seq_, 7, false);
		}
	}
	if(getRevComp){
		for(const auto & rSeq : refSeqs){
			revComp_refSeqs.emplace_back(rSeq);
			revComp_refSeqs.back().reverseComplementRead(false, true);

			revComp_refSeqsKInfos.emplace_back(revComp_refSeqs.back().seq_, 7, false);
		}
	}

	aligner alignerObj(maxLen, setUp.pars_.gapInfo_, setUp.pars_.scoring_, false);
	alignerObj.weighHomopolymers_ = false;
	alignerObj.processAlnInfoInputNoCheck(setUp.pars_.alnInfoDirName_, setUp.pars_.verbose_);

	reader.in_.reOpenIn();

	while(reader.readNextRead(seq)){
		std::vector<seqInfo> trimmed;
		if(getRevComp){
			trimmed = readVecTrimmer::trimSeqToRefByGlobalAlnBestNoOverlapIncludeRevComp(seq, refSeqs, revComp_refSeqs, refSeqsKInfos, revComp_refSeqsKInfos, alignerObj);;
		}else{
			trimmed = readVecTrimmer::trimSeqToRefByGlobalAlnBestNoOverlap(seq, refSeqs, refSeqsKInfos, alignerObj);;
		}
//		trimSeqToRefByGlobalAln(seq, refSeqs, refSeqsKInfos, alignerObj);
//		alignerObj.alignObjectA_.seqBase_.outPutSeq(std::cout);
//		alignerObj.alignObjectB_.seqBase_.outPutSeq(std::cout);
		for(const auto & trimmedSeq : trimmed){
			if(!pars.keepOnlyOn || trimmedSeq.on_){
				reader.write(trimmedSeq);
			}
		}
	}

	alignerObj.processAlnInfoOutputNoCheck(setUp.pars_.alnInfoDirName_, setUp.pars_.verbose_);
	if(setUp.pars_.verbose_){
		setUp.logRunTime(std::cout);
	}
	return 0;
}

} // namespace njhseq
