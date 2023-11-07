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
	//auto outSeqOpts = SeqIOOptions::genFastaOut("");
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();




	setUp.setOption(trimPars.extend_, "--extend", "extend this amount by adding this length to the front from the back and to the back from the front");

	setUp.setOption(trimPars.kmerLength_, "--kmerLength", "kmer Length");
//	setUp.processReadInNames(true);
//	outSeqOpts.out_.outFilename_ = njh::files::prependFileBasename(setUp.pars_.ioOptions_.firstName_, "trimmed_");;
//	setUp.processWritingOptions(outSeqOpts.out_);

	setUp.processDefaultReader(true);
	if(setUp.pars_.ioOptions_.out_.outFilename_ == "out"){
		setUp.pars_.ioOptions_.out_.outFilename_ = njh::files::prependFileBasename(setUp.pars_.ioOptions_.firstName_, "trimmed_");
	}

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
	//setUp.processDirectoryOutputName(false);
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

	SeqOutput writer(setUp.pars_.ioOptions_);
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
	maxSize += trimPars.extend_ * 2;
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
//		seqInfo originalSeq = seq;
//		if(trimPars.extend_ > 0){
//			if(len(seq) > trimPars.extend_){
//				auto front = seq.getSubRead(0, trimPars.extend_);
//				auto back = seq.getSubRead(len(seq) - trimPars.extend_);
//				seq.prepend(back);
//				seq.append(front);
//			}else{
//				seq.append(seq);
//				seq.prepend(seq);
//			}
//		}
		seqInfo originalSeq = seq;
		if(trimPars.extend_ > 0){
			seqInfo front;
			seqInfo back;
			if(len(seq) > trimPars.extend_){
				front = seq.getSubRead(0, trimPars.extend_);
				back = seq.getSubRead(len(seq) - trimPars.extend_);
			}else{
				front = seq;
				back = seq;
			}
			// since these are circular seqs check the fronts and back for same seq
			uint32_t checkLenTo = trimPars.extendSeqCheckLenTo_;
			//to adjust for seq len;
			checkLenTo = std::min<uint32_t>(checkLenTo, len(seq));
			uint32_t checkLenFrom = trimPars.extendSeqCheckLenFrom_;
			std::cout << __FILE__ << " " << __LINE__ << std::endl;
			std::cout << "checkLenFrom: " << checkLenFrom << std::endl;
			std::cout << "checkLenTo: " << checkLenTo << std::endl;

			if(checkLenFrom >= checkLenTo){
				checkLenFrom = checkLenTo - 1;
			}
			bool sameSeqFrontBack = false;
			uint32_t sameSeqSize = std::numeric_limits<uint32_t>::max();
			for(const auto pos : iter::range(checkLenFrom, checkLenTo)){
				if(std::equal(seq.seq_.begin(),                  seq.seq_.begin() + pos + 1,
						          seq.seq_.begin() + len(seq) - pos - 1 )){
					std::cout << "pos: " << pos << std::endl;
					std::cout << "pos + 1: " << pos + 1 << std::endl;
					std::cout << seq.seq_.substr(0, pos + 1) << std::endl;;
					std::cout << seq.seq_.substr(len(seq)-pos - 1, pos + 1) << std::endl;;
					sameSeqFrontBack = true;
					sameSeqSize = pos + 1;
					break;
				}
			}
			if(sameSeqFrontBack){
				if(len(front) > sameSeqSize){
					front = front.getSubRead(sameSeqSize);
					back = front.getSubRead(0, len(back) - sameSeqSize);
					seq.prepend(back);
					seq.append(front);
				}
			}else{
				seq.prepend(back);
				seq.append(front);
			}
		}

		auto res = readVecTrimmer::trimCircularGenomeToRef(seq, trimPars, alignerObj);

		if(trimPars.extend_ == 0){
			// if no extension just add results
			addOtherVec(outSeqs, res);
		}else{
			if(1 == res.size() && !res.front().on_){
				//no trim, adding original seq in case it was extended
				originalSeq.on_ = false;
				outSeqs.emplace_back(originalSeq);
			}else if(1 == res.size() && res.front().on_){
				//trimmed to one piece add results
				addOtherVec(outSeqs, res);
			} else {
				//there was trimming but not complete, since we extended just trim to the original seq
				auto res = readVecTrimmer::trimCircularGenomeToRef(originalSeq, trimPars, alignerObj);
				addOtherVec(outSeqs, res);
			}
		}
	}


	writer.write(outSeqs);
	alignerObj.processAlnInfoOutput(setUp.pars_.outAlnInfoDirName_, setUp.pars_.verbose_);

	return 0;
}



}  // namespace njhseq
