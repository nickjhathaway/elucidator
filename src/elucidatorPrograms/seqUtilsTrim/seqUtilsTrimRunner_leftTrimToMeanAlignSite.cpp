/*
 * seqUtilsTrimRunner_leftTrimToMeanAlignSite.cpp
 *
 *  Created on: Mar 4, 2020
 *      Author: nicholashathaway
 */


#include "seqUtilsTrimRunner.hpp"
#include <njhseq/IO/SeqIO/SeqIO.hpp>
#include <njhseq/alignment/aligner/aligner.hpp>

namespace njhseq {




int seqUtilsTrimRunner::leftTrimToMeanAlignSite(const njh::progutils::CmdArgs & inputCommands){
	bool mark = false;
	seqUtilsTrimSetUp setUp(inputCommands);
	setUp.description_ = "Trim the left side of sequences to the mean alignment site of reference sequences";
	FullTrimReadsPars pars;
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processRefFilename(true);
	setUp.pars_.gapLeft_ = "0,0";
	setUp.pars_.gapRight_ = "0,0";
	setUp.pars_.gap_ = "5,1";
	setUp.setOption(mark, "--mark", "Mark the trimmed sequences with the trim start determined");
	setUp.processAlignerDefualts();

	setUp.processIoOptions(setUp.readInFormatsAvailable_);
	setUp.finishSetUp(std::cout);

	if(setUp.pars_.debug_){
		std::cout << setUp.pars_.gapInfo_.toJson() << std::endl;
	}

	uint64_t maxLen = 0;
	{
		SeqIO reader(setUp.pars_.ioOptions_);
		reader.openIn();
		seqInfo seq;
		while(reader.readNextRead(seq)){
			readVec::getMaxLength(seq, maxLen);
		}
	}
	auto refSeqs = SeqInput::getSeqVec<seqInfo>(setUp.pars_.refIoOptions_, maxLen);

	aligner alignerObj(maxLen, setUp.pars_.gapInfo_, setUp.pars_.scoring_);
	alignerObj.processAlnInfoInputNoCheck(setUp.pars_.alnInfoDirName_, setUp.pars_.verbose_);
	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	reader.openOut();
	seqInfo seq;

	while(reader.readNextRead(seq)){
		std::vector<uint32_t> starts;
		if(setUp.pars_.debug_ || setUp.pars_.verbose_){
			std::cout << seq.name_ << std::endl;
		}
		for(const auto & ref : refSeqs){
			alignerObj.alignCacheGlobal(seq, ref);
			alignerObj.rearrangeObjsGlobal(seq, ref);
			uint32_t start = 0;
			if(alignerObj.alignObjectB_.seqBase_.seq_.front() == '-'){
				start = getRealPosForAlnPos(alignerObj.alignObjectA_.seqBase_.seq_, alignerObj.alignObjectB_.seqBase_.seq_.find_first_not_of("-"));
			}
			if(setUp.pars_.debug_){
				std::cout << '\t' << ref.name_ << std::endl;
				std::cout << "\t" << start << std::endl;
			}
			starts.emplace_back(start);
		}
		if(setUp.pars_.debug_){
			std::cout << seq.name_ << std::endl;
			std::cout << vectorMedianRef(starts) << " " << std::round(vectorMedianRef(starts)) << std::endl;
		}
		uint32_t start  = std::round(vectorMedianRef(starts));

		if(start > 0){
			seq = seq.getSubRead(start);
		}
		if (mark) {
			MetaDataInName seqMeta;
			if (MetaDataInName::nameHasMetaData(seq.name_)) {
				seqMeta = MetaDataInName(seq.name_);
			}
			seqMeta.addMeta("trimStart", start, true);
			seqMeta.resetMetaInName(seq.name_);
		}
		reader.write(seq);
	}

	alignerObj.processAlnInfoOutputNoCheck(setUp.pars_.outAlnInfoDirName_, setUp.pars_.verbose_);

	if(setUp.pars_.verbose_){
		setUp.logRunTime(std::cout);
	}
	return 0;
}


}  // namespace njhseq
