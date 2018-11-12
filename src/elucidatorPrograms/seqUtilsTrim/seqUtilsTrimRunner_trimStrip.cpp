/*
 * seqUtilsTrimRunner_trimStrip.cpp
 *
 *  Created on: Jun 8, 2017
 *      Author: nick
 */


#include "seqUtilsTrimRunner.hpp"


namespace njhseq {


int seqUtilsTrimRunner::breakupAtRegexPat(const njh::progutils::CmdArgs & inputCommands){
	std::string patStr = "N+";
	bool mark = false;
	seqUtilsTrimSetUp setUp(inputCommands);
	setUp.processDefaultReader();
	if("" == setUp.pars_.ioOptions_.out_.outFilename_){
		setUp.pars_.ioOptions_.out_.outFilename_ = njh::files::prependFileBasename(setUp.pars_.ioOptions_.firstName_, "breakup_");
	}else if("STDIN" == setUp.pars_.ioOptions_.firstName_){
		setUp.pars_.ioOptions_.out_.outFilename_ = "STDOUT";
	}
	setUp.setOption(patStr, "--pattern", "Pattern to break up seq on");
	setUp.setOption(mark, "--mark", "mark seq with postions");

	setUp.processVerbose();
	setUp.processDebug();
	setUp.finishSetUp(std::cout);

	std::regex pat{patStr};
	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	reader.openOut();

	seqInfo seq;

	struct PatPosSize{
		PatPosSize(const std::string & pat, size_t pos): pat_(pat), pos_(pos){

		}
		std::string pat_;
		size_t pos_;

		size_t end(){
			return pos_ + pat_.size();
		}
	};

	while(reader.readNextRead(seq)){

    std::sregex_iterator iter(seq.seq_.begin(), seq.seq_.end(), pat);
    std::sregex_iterator end;
  	std::vector<PatPosSize> pats;
		while (iter != end) {
			pats.emplace_back((*iter)[0], iter->position());
			++iter;
		}

		if(!pats.empty()){
			if(0 != pats.front().pos_ ){
				size_t start = 0;
				size_t end = pats.front().pos_;
				auto subSeq = seq.getSubRead(start, end - start);
				if(mark){
					subSeq.name_.append(njh::pasteAsStr("-s", start, "-e", end));
				}
				reader.write(subSeq);
			}
			if(pats.size() > 1){
				for(const auto & patPos : iter::range(pats.size() - 1)){
					const auto & p = pats[patPos];
					size_t start = p.pos_ + p.pat_.size();
					size_t end = pats[patPos + 1].pos_;
					auto subSeq = seq.getSubRead(start, end - start);
					if(mark){
						subSeq.name_.append(njh::pasteAsStr("-s", start, "-e", end));
					}
					reader.write(subSeq);
				}
			}
			if(seq.seq_.size() != pats.back().end()){
				size_t start = pats.back().end();
				size_t end = seq.seq_.size();
				auto subSeq = seq.getSubRead(start, end - start);
				if(mark){
					subSeq.name_.append(njh::pasteAsStr("-s", start, "-e", end));
				}
				reader.write(subSeq);
			}
		}else{
			reader.write(seq);
		}
	}


	return 0;
}

int seqUtilsTrimRunner::trimLstripBase(const njh::progutils::CmdArgs & inputCommands){

	seqUtilsTrimSetUp setUp(inputCommands);
	FullTrimReadsPars pars;

	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(pars.base, "--base", "Base to strip off the left side of the sequence", true);
	setUp.processIoOptions();
	setUp.setOption(pars.keepOnlyOn, "--keepOnlyOn",
			"Keep Only the Reads that are still on");
	setUp.finishSetUp(std::cout);

	SeqIO reader(setUp.pars_.ioOptions_);

	seqInfo seq;
	reader.openIn();
	reader.openOut();
	while(reader.readNextRead(seq)){
		readVecTrimmer::trimAtLstripBase(seq, pars.base);
		if(seq.seq_.empty()){
			seq.on_ = false;
		}
		if(!pars.keepOnlyOn || seq.on_){
			reader.write(seq);
		}
	}
	return 0;
}

int seqUtilsTrimRunner::trimRstripBase(const njh::progutils::CmdArgs & inputCommands){
	seqUtilsTrimSetUp setUp(inputCommands);
	FullTrimReadsPars pars;

	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(pars.base, "--base", "Base to strip off the right side of the sequence", true);
	setUp.processIoOptions();
	setUp.setOption(pars.keepOnlyOn, "--keepOnlyOn",
			"Keep Only the Reads that are still on");
	setUp.finishSetUp(std::cout);

	SeqIO reader(setUp.pars_.ioOptions_);

	seqInfo seq;
	reader.openIn();
	reader.openOut();
	while(reader.readNextRead(seq)){
		readVecTrimmer::trimAtRstripBase(seq, pars.base);
		if(seq.seq_.empty()){
			seq.on_ = false;
		}
		if(!pars.keepOnlyOn || seq.on_){
			reader.write(seq);
		}
	}
	return 0;
}

int seqUtilsTrimRunner::trimStripBase(const njh::progutils::CmdArgs & inputCommands){
	char rightBase = ' ';
	seqUtilsTrimSetUp setUp(inputCommands);
	FullTrimReadsPars pars;

	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(pars.base, "--leftBase", "Base to strip off the left side of the sequence", true);
	if(!setUp.setOption(rightBase, "--rightBase", "Base to strip off the right side of the sequence, defaults to leftBase")){
		rightBase = pars.base;
	}

	setUp.processIoOptions();
	setUp.setOption(pars.keepOnlyOn, "--keepOnlyOn",
			"Keep Only the Reads that are still on");
	setUp.finishSetUp(std::cout);

	SeqIO reader(setUp.pars_.ioOptions_);

	seqInfo seq;
	reader.openIn();
	reader.openOut();
	while(reader.readNextRead(seq)){
		readVecTrimmer::trimAtLstripBase(seq, pars.base);
		readVecTrimmer::trimAtRstripBase(seq, rightBase);
		if(seq.seq_.empty()){
			seq.on_ = false;
		}
		if(!pars.keepOnlyOn || seq.on_){
			reader.write(seq);
		}
	}
	return 0;
}

int seqUtilsTrimRunner::trimLstripQual(const njh::progutils::CmdArgs & inputCommands){
	seqUtilsTrimSetUp setUp(inputCommands);
	FullTrimReadsPars pars;

	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(pars.qual, "--qScore", "quality to strip off the left side of the sequence", true);
	setUp.processIoOptions();
	setUp.setOption(pars.keepOnlyOn, "--keepOnlyOn",
			"Keep Only the Reads that are still on");
	setUp.finishSetUp(std::cout);

	SeqIO reader(setUp.pars_.ioOptions_);

	seqInfo seq;
	reader.openIn();
	reader.openOut();
	while(reader.readNextRead(seq)){
		readVecTrimmer::trimAtLstripQualScore(seq, pars.qual);
		if(seq.seq_.empty()){
			seq.on_ = false;
		}
		if(!pars.keepOnlyOn || seq.on_){
			reader.write(seq);
		}
	}
	return 0;
}

int seqUtilsTrimRunner::trimRstripQual(const njh::progutils::CmdArgs & inputCommands){
	seqUtilsTrimSetUp setUp(inputCommands);
	FullTrimReadsPars pars;

	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(pars.qual, "--qScore", "quality to strip off the right side of the sequence", true);
	setUp.processIoOptions();
	setUp.setOption(pars.keepOnlyOn, "--keepOnlyOn",
			"Keep Only the Reads that are still on");
	setUp.finishSetUp(std::cout);

	SeqIO reader(setUp.pars_.ioOptions_);

	seqInfo seq;
	reader.openIn();
	reader.openOut();
	while(reader.readNextRead(seq)){
		readVecTrimmer::trimAtRstripQualScore(seq, pars.qual);
		if(seq.seq_.empty()){
			seq.on_ = false;
		}
		if(!pars.keepOnlyOn || seq.on_){
			reader.write(seq);
		}
	}
	return 0;
}

int seqUtilsTrimRunner::trimStripQual(const njh::progutils::CmdArgs & inputCommands){
	seqUtilsTrimSetUp setUp(inputCommands);
	FullTrimReadsPars pars;
	uint32_t rightQual = 2;
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(pars.qual, "--leftQScore", "quality to strip off the left side of the sequence", true);
	if(!setUp.setOption(rightQual, "--rightQual", "quality to strip off the right side of the sequence")){
		rightQual = pars.qual;
	}
	setUp.processIoOptions();
	setUp.setOption(pars.keepOnlyOn, "--keepOnlyOn",
			"Keep Only the Reads that are still on");
	setUp.finishSetUp(std::cout);

	SeqIO reader(setUp.pars_.ioOptions_);

	seqInfo seq;
	reader.openIn();
	reader.openOut();
	while(reader.readNextRead(seq)){
		readVecTrimmer::trimAtLstripQualScore(seq, pars.qual);
		readVecTrimmer::trimAtRstripQualScore(seq, rightQual);
		if(seq.seq_.empty()){
			seq.on_ = false;
		}
		if(!pars.keepOnlyOn || seq.on_){
			reader.write(seq);
		}
	}

	return 0;
}


} // namespace njhseq

