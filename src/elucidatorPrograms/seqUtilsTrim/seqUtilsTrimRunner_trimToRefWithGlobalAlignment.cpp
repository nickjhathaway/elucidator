/*
 * seqUtilsTrimRunner_trimToRefWithGlobalAlignment.cpp
 *
 *  Created on: Nov 26, 2017
 *      Author: nick
 */





#include "seqUtilsTrimRunner.hpp"


namespace njhseq {

int seqUtilsTrimRunner::trimToRefWithGlobalAlignment(const njh::progutils::CmdArgs & inputCommands) {
	seqUtilsTrimSetUp setUp(inputCommands);
	setUp.pars_.gapInfo_ = gapScoringParameters(5, 1, 0, 0, 0, 0);
	setUp.pars_.gapRight_ = "0,0";
	setUp.pars_.gapLeft_ = "0,0";
	FullTrimReadsPars pars;
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processGap();
	setUp.processScoringPars();
	setUp.processAlnInfoInput();
	setUp.processIoOptions();
	setUp.processRefFilename(true);
	setUp.setOption(pars.keepOnlyOn, "--keepOnlyOn", "Keep Only the Reads that are still on");
	setUp.finishSetUp(std::cout);

	SeqIO reader(setUp.pars_.ioOptions_	);
	reader.openIn();
	reader.openOut();
	uint64_t maxLen = 0;
	seqInfo seq;
	while(reader.readNextRead(seq)){
		readVec::getMaxLength(seq, maxLen);
	}

	auto refSeqs = SeqInput::getReferenceSeq(setUp.pars_.refIoOptions_, maxLen);
	std::vector<kmerInfo> refSeqsKInfos;
	if(refSeqs.size() != 1){
		for(const auto & rSeq : refSeqs){
			refSeqsKInfos.emplace_back(rSeq.seqBase_.seq_, 7, false);
		}
	}

	aligner alignerObj(maxLen, setUp.pars_.gapInfo_, setUp.pars_.scoring_, false);
	alignerObj.weighHomopolymers_ = false;
	alignerObj.processAlnInfoInputNoCheck(setUp.pars_.alnInfoDirName_, setUp.pars_.verbose_);

	reader.in_.reOpenIn();

	while(reader.readNextRead(seq)){
		readVecTrimmer::trimSeqToRefByGlobalAln(seq, refSeqs, refSeqsKInfos, alignerObj);
//		trimSeqToRefByGlobalAln(seq, refSeqs, refSeqsKInfos, alignerObj);
//		alignerObj.alignObjectA_.seqBase_.outPutSeq(std::cout);
//		alignerObj.alignObjectB_.seqBase_.outPutSeq(std::cout);
		if(!pars.keepOnlyOn || seq.on_){
			reader.write(seq);
		}
	}

	alignerObj.processAlnInfoOutputNoCheck(setUp.pars_.alnInfoDirName_, setUp.pars_.verbose_);

	return 0;

}




int seqUtilsTrimRunner::trimToRefWithGlobalAlignmentToRefPositions(const njh::progutils::CmdArgs & inputCommands) {
	FullTrimReadsPars pars;
	readVecTrimmer::GlobalAlnTrimPars gTrimPars;
	seqUtilsTrimSetUp setUp(inputCommands);
	setUp.pars_.gapInfo_ = gapScoringParameters(5, 1, 0, 0, 0, 0);
	setUp.pars_.gapRight_ = "0,0";
	setUp.pars_.gapLeft_ = "0,0";
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processGap();
	setUp.processScoringPars();
	setUp.processAlnInfoInput();
	setUp.processIoOptions();
	setUp.processRefFilename(true);
	setUp.setOption(gTrimPars.startInclusive_, "--startInclusive", "The start position of the given reference sequence to trim to, (inclusive, zero-based)", true);
	setUp.setOption(gTrimPars.endInclusive_,   "--endInclusive",   "The end position of the given reference sequence to trim to, (inclusive, zero-based)", true);
	setUp.setOption(pars.keepOnlyOn, "--keepOnlyOn",    "Keep Only the Reads that are still on");
	setUp.setOption(gTrimPars.needJustOneEnd_, "--needJustOneEnd",    "Need to just one end to be left on");
	setUp.finishSetUp(std::cout);

	SeqIO reader(setUp.pars_.ioOptions_	);
	reader.openIn();
	reader.openOut();
	uint64_t maxLen = 0;
	seqInfo seq;
	while(reader.readNextRead(seq)){
		readVec::getMaxLength(seq, maxLen);
	}

	auto refSeqs = SeqInput::getReferenceSeq(setUp.pars_.refIoOptions_, maxLen);

	aligner alignerObj(maxLen, setUp.pars_.gapInfo_, setUp.pars_.scoring_, false);
	alignerObj.weighHomopolymers_ = false;
	alignerObj.processAlnInfoInputNoCheck(setUp.pars_.alnInfoDirName_, setUp.pars_.verbose_);

	reader.in_.reOpenIn();
	uint32_t seqCount = 0;
	while (reader.readNextRead(seq)) {
		++seqCount;
		if(setUp.pars_.verbose_ && seqCount % 10  == 0){
			std::cout << '\r' << seqCount;
			std::cout.flush();
		}
		readVecTrimmer::trimSeqToRefByGlobalAln(seq, refSeqs.front(), gTrimPars,
				alignerObj);

		if (!pars.keepOnlyOn || seq.on_) {
			reader.write(seq);
		}
	}
	if(setUp.pars_.verbose_){
		std::cout << std::endl;
	}
	alignerObj.processAlnInfoOutputNoCheck(setUp.pars_.alnInfoDirName_, setUp.pars_.verbose_);
	return 0;

}

int seqUtilsTrimRunner::trimToRefWithGlobalAlignmentToRefMultiplePositions(const njh::progutils::CmdArgs & inputCommands) {
	FullTrimReadsPars pars;
	CutOutputRefPositions refPositionsSingular;
	bfs::path refPostionParsFnp = "";
	seqUtilsTrimSetUp setUp(inputCommands);
	setUp.processComparison(refPositionsSingular.comp_,   "--alignCutComp");
	setUp.setOption(refPositionsSingular.checkComp_, "--compAlignCut", "Compare on alignment cut");
	setUp.setOption(refPositionsSingular.refStopStr_, "--refStop", "ref stop");
	setUp.setOption(refPositionsSingular.refStartStr_, "--refStart", "ref start");

	setUp.pars_.gapInfo_ = gapScoringParameters(5, 1, 0, 0, 0, 0);
	setUp.pars_.gapRight_ = "0,0";
	setUp.pars_.gapLeft_ = "0,0";
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processGap();
	setUp.processScoringPars();
	setUp.processAlnInfoInput();
	setUp.processDefaultReader(VecStr{"--fasta", "--fastq"}, true);
	setUp.processRefFilename(true);
	setUp.setOption(refPostionParsFnp, "--refPostionParsFnp", "Ref Postion Pars Fnp of a json file with multiple positions");
	setUp.finishSetUp(std::cout);
	std::vector<CutOutputRefPositions> refPositionsVec;
	if ("" != refPostionParsFnp) {
		refPositionsVec = CutOutputRefPositions::readInPositions(refPostionParsFnp);
	} else {
		CutOutputRefPositions::setPositions(refPositionsSingular, __PRETTY_FUNCTION__);
		refPositionsVec = std::vector<CutOutputRefPositions> { refPositionsSingular };
	}
	if(setUp.pars_.debug_){
		std::cerr << "refPositionsVec.size(): " << refPositionsVec.size() << std::endl;
		std::cerr << njh::json::toJson(refPositionsVec) << std::endl;
	}
	SeqIO reader(setUp.pars_.ioOptions_	);
	reader.openIn();
	uint64_t maxLen = 0;
	seqInfo seq;
	while(reader.readNextRead(seq)){
		readVec::getMaxLength(seq, maxLen);
	}

	auto refSeqs = SeqInput::getSeqVec<seqInfo>(setUp.pars_.refIoOptions_, maxLen);
	if(refSeqs.empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << setUp.pars_.refIoOptions_.firstName_ << " was empty" << "\n";
		throw std::runtime_error{ss.str()};
	}
	auto refSeq = refSeqs.front();
	MultiSeqOutCache<seqInfo> seqOuts;
	/**@todo add a check to make sure ref positions aren't the same*/
	for(auto & refPositions : refPositionsVec){
		if(std::numeric_limits<uint32_t>::max() == refPositions.refStop_){
			refPositions.refStop_ = refSeq.seq_.size() - 1;
		}
		refPositions.checkPositionsThrow(__PRETTY_FUNCTION__, len(refSeq), refSeq.name_);
		SeqIOOptions refPosSeqOpts = setUp.pars_.ioOptions_;
		refPosSeqOpts.out_.outFilename_ = njh::files::prependFileBasename(
				refPosSeqOpts.out_.outFilename_,
				njh::pasteAsStr("trimmed_", refPositions.refStart_, "-",
						refPositions.refStop_, "_"));
		seqOuts.addReader(refPositions.getId(), refPosSeqOpts);
	}

	aligner alignerObj(maxLen, setUp.pars_.gapInfo_, setUp.pars_.scoring_, false);
	alignerObj.weighHomopolymers_ = false;
	alignerObj.processAlnInfoInputNoCheck(setUp.pars_.alnInfoDirName_, setUp.pars_.verbose_);

	reader.in_.reOpenIn();
	uint32_t seqCount = 0;


	while (reader.readNextRead(seq)) {
		++seqCount;
		if(setUp.pars_.verbose_ && seqCount % 10  == 0){
			std::cout << '\r' << seqCount;
			std::cout.flush();
		}
		//readVecTrimmer::trimSeqToRefByGlobalAln(seq, refSeqs[refPos], gTrimPars, alignerObj);
		alignerObj.alignCacheGlobal(refSeq, getSeqBase(seq));
		//aligned bases
		uint32_t firstAlignedBase = alignerObj.alignObjectB_.seqBase_.seq_.find_first_not_of('-');
		uint32_t lastAlignedBase = alignerObj.alignObjectB_.seqBase_.seq_.find_last_not_of('-');

		for(const auto & refPositions : refPositionsVec){
			//ref positions
			uint32_t refAlnStartPos = alignerObj.getAlignPosForSeqAPos(refPositions.refStart_);
			uint32_t refAlnStopPos = alignerObj.getAlignPosForSeqAPos(refPositions.refStop_);

			//front
			uint32_t trimFront = std::numeric_limits<uint32_t>::max();
			if (refAlnStartPos >= firstAlignedBase
					&& refAlnStartPos < lastAlignedBase) {
				trimFront = alignerObj.getSeqPosForAlnBPos(refAlnStartPos);
			}
			//back
			uint32_t trimBack = std::numeric_limits<uint32_t>::max();
			if (refAlnStopPos > firstAlignedBase && refAlnStopPos <= lastAlignedBase) {
				trimBack = alignerObj.getSeqPosForAlnBPos(refAlnStopPos);
			}

			//trim if found both
			if (std::numeric_limits<uint32_t>::max() != trimFront
					&& std::numeric_limits<uint32_t>::max() != trimBack) {
				seqInfo seqCopy = seq;

				seqCopy.setClip(trimFront, trimBack);
				auto malnRefStart = alignerObj.getAlignPosForSeqAPos(refPositions.refStart_);
				auto malnRefStop =  alignerObj.getAlignPosForSeqAPos(refPositions.refStop_);
				if(refPositions.refStartLength_ > 1){
					if(refPositions.checkComp_){
						alignerObj.profileAlignment(seq, refSeq, false, false, false, malnRefStart, malnRefStart + refPositions.refStartLength_);
						if(!refPositions.comp_.passErrorProfile(alignerObj.comp_)){
							seqCopy.on_ = false;
						}
					} else {
						for (const auto pos : iter::range(malnRefStart,
								malnRefStart + refPositions.refStartLength_)) {
							if ('-' == alignerObj.alignObjectB_.seqBase_.seq_[pos]) {
								seqCopy.on_ = false;
								break;
							}
						}
					}
				}
				if(seqCopy.on_){
					if(refPositions.refStopLength_ > 1){
						if(refPositions.checkComp_){
							alignerObj.profileAlignment(seq, refSeq, false, false, false, malnRefStop + 1 - refPositions.refStopLength_, malnRefStop + 1);
							if(!refPositions.comp_.passErrorProfile(alignerObj.comp_)){
								seqCopy.on_ = false;
							}
						}else{
							for(const auto pos : iter::range(malnRefStop + 1 - refPositions.refStopLength_, malnRefStop + 1)){
								if('-' == alignerObj.alignObjectB_.seqBase_.seq_[pos]){
									seqCopy.on_ = false;
									break;
								}
							}
						}
					}
				}
				if(seqCopy.on_){
					seqOuts.add(refPositions.getId(), seqCopy);
				}
			}
		}
	}

	if(setUp.pars_.verbose_){
		std::cout << std::endl;
	}
	alignerObj.processAlnInfoOutputNoCheck(setUp.pars_.alnInfoDirName_, setUp.pars_.verbose_);
	return 0;

}




} // namespace njhseq


