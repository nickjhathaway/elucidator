/*
 * seqUtilsTrimRunner_trimStrip.cpp
 *
 *  Created on: Jun 8, 2017
 *      Author: nick
 */
// elucidator - A library for analyzing sequence data
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
//
// This file is part of elucidator.
//
// elucidator is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// elucidator is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with elucidator.  If not, see <http://www.gnu.org/licenses/>.
//

#include "seqUtilsTrimRunner.hpp"


namespace njhseq {


int seqUtilsTrimRunner::breakupAtRegexPat(const njh::progutils::CmdArgs & inputCommands){
	std::string patStr = "N+";
	bool mark = false;
	seqUtilsTrimSetUp setUp(inputCommands);
	setUp.processDefaultReader();
	if ("" == setUp.pars_.ioOptions_.out_.outFilename_) {
		setUp.pars_.ioOptions_.out_.outFilename_ = njh::files::prependFileBasename(
				setUp.pars_.ioOptions_.firstName_, "breakup_");
	} else if ("STDIN" == setUp.pars_.ioOptions_.firstName_) {
		setUp.pars_.ioOptions_.out_.outFilename_ = "STDOUT";
	}
	setUp.setOption(patStr, "--pattern", "Pattern to break up seq on");
	setUp.setOption(mark, "--mark", "mark seq with postions");

	setUp.processVerbose();
	setUp.processDebug();
	setUp.finishSetUp(std::cout);

	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	reader.openOut();

	seqInfo seq;

	njh::PatPosFinder pFinder(patStr);
	while(reader.readNextRead(seq)){
		auto breaks = readVecTrimmer::breakUpSeqOnPat(seq, pFinder);
		for(auto & seqBreak : breaks){
			if(mark && "" != seqBreak.pat_ ){
				seqBreak.seqBase_.name_.append(njh::pasteAsStr("-s", seqBreak.start_, "-e", seqBreak.end_));
			}
			reader.write(seqBreak.seqBase_);
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

