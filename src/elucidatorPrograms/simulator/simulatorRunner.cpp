//
//  simulator.cpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 01/26/2013.
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

#include "simulatorRunner.hpp"
#include "elucidator/simulation.h"
#include "elucidator/seqToolsUtils/seqToolsUtils.hpp"
#include <njhseq/helpers/clusterCollapser.hpp>


namespace njhseq {



simulatorRunner::simulatorRunner()
    : njh::progutils::ProgramRunner({addFunc("randomStrings", randomStrings, false),
                     addFunc("randomStrsWithKmers", randomStrsWithKmers, false),
										 addFunc("randomSeqFile", randomSeqFile, false),
										 addFunc("randomSampleFile", randomSampleFile, false),
										 addFunc("randomSampleFast", randomSampleFast, false)},
    		//
                    "simulator") {}


template<typename T>
void randomSampleSeqFile(const SeqIOOptions & opts, const std::string& sample, bool verbose){

  uint64_t readNumber = countSeqs(opts, verbose);
  uint32_t sampleNum = processRunCutoff(sample, readNumber);
  SeqIO reader(opts);
  reader.openIn();
  reader.openOut();

  njh::randomGenerator gen;
  std::vector<uint64_t> allPositions(readNumber);
  njh::iota<uint64_t>(allPositions, 0);
  auto positions = gen.unifRandSelectionVec(allPositions, sampleNum, false);
  auto maxPos = *std::max_element(positions.begin(), positions.end());
  uint64_t pos = 0;
  uint32_t outCount = 0;
  if(verbose){
  	std::cout << "Sampling: " << sampleNum << " from " << opts.firstName_ << std::endl;
  	std::cout << "positions.size(): " << positions.size() << std::endl;
  	std::cout << "maxPos: " << maxPos << std::endl;
  }
	T seq;
  while(reader.readNextRead(seq)){
  	if(pos > maxPos){
  		break;
  	}
  	if(verbose && (outCount + 1) % 50 == 0){
  		std::cout << "\r" << outCount + 1 << ":" << sampleNum;
  	}
  	if(njh::in(pos, positions)){
  		++outCount;
  		reader.write(seq);
  	}
  	++pos;
  }
}



void randomSampleFastFile(const SeqIOOptions & opts,
		const std::string& sample){
	if (opts.inFormat_ != SeqIOOptions::inFormats::FASTQ
			&& opts.inFormat_ != SeqIOOptions::inFormats::FASTA) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": format should be fastq or fasta\n";
		throw std::runtime_error{ss.str()};
	}
  njh::randomGenerator gen;

  SeqIO reader(opts);
  reader.openIn();
  reader.in_.loadIndex();
  auto positions = reader.in_.randomlySampleIndex(gen,sample);
  reader.openOut();

	seqInfo seq;
	for(const auto pos : positions){
		reader.in_.seekgPri(pos);
		if(reader.readNextRead(seq)){
			reader.write(seq);
		}
	}
}


int simulatorRunner::randomSampleFast(const njh::progutils::CmdArgs & inputCommands){
	seqSetUp setUp(inputCommands);
	std::string sample = "";
	setUp.setOption(sample, "-n,--sample", "Sample, either absolute number or a percentage", true);
	setUp.processDefaultReader(VecStr{"-fasta", "-fastq"}, true);
	setUp.processVerbose();
	setUp.processDebug();
  setUp.finishSetUp(std::cout);
  randomSampleFastFile(setUp.pars_.ioOptions_, sample);
  return 0;
}

int simulatorRunner::randomSampleFile(const njh::progutils::CmdArgs & inputCommands){
	seqSetUp setUp(inputCommands);
	std::string sample = "";
	setUp.setOption(sample, "-n,--sample", "Sample, either absolute number or a percentage", true);
	setUp.processDefaultReader(true);
	setUp.processVerbose();
	setUp.processDebug();
  setUp.finishSetUp(std::cout);
  if(SeqIOOptions::inFormats::FASTQPAIRED == setUp.pars_.ioOptions_.inFormat_){
  	randomSampleSeqFile<PairedRead>(setUp.pars_.ioOptions_, sample, setUp.pars_.verbose_);
  }else{
  	randomSampleSeqFile<seqInfo>(setUp.pars_.ioOptions_, sample, setUp.pars_.verbose_);
  }
  if(setUp.pars_.verbose_){
  	std::cout << std::endl;
  }
  return 0;
}


int simulatorRunner::randomStrsWithKmers(const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	std::string kmersStrs = "AAAA,CCCC,GGGG,TTTT";
	bool header = false;
	bool ansi = false;
	uint32_t length = 10;
	uint32_t num = 10;
	setUp.setOption(ansi, "-ansi,-color", "ansiColorOutput");
	setUp.setOption(length, "-len,-length", "length");
	setUp.setOption(num, "-num,-number", "number");
	setUp.setOption(kmersStrs, "-kmers,-strs", "strings_kmers");
	setUp.setOption(setUp.pars_.ioOptions_.processed_, "-processed", "processed");
	setUp.setOption(header, "-header", "header");
	setUp.finishSetUp(std::cout);

	std::vector<readObject> inReads;

	VecStr kmers;
	if (bfs::exists(kmersStrs)) {
		std::string extention = njh::files::getExtension(kmersStrs);
		if (extention == "fasta" || extention == "fastq") {
			setUp.pars_.ioOptions_.firstName_ = kmersStrs;
			setUp.pars_.ioOptions_.inFormat_ = SeqIOOptions::getInFormat(extention);
			SeqInput reader(setUp.pars_.ioOptions_);
			reader.openIn();
			auto inReads = reader.readAllReads<readObject>();
		} else {
			table inTab(kmersStrs, "whitespace", header);
			if (inTab.columnNames_.size() == 2) {
				for (const auto & row : inTab.content_) {
					inReads.emplace_back(readObject(seqInfo(row[0], row[0])));
					inReads.back().seqBase_.cnt_ = std::stoi(row[1]);
				}
			} else {
				for (const auto & row : inTab.content_) {
					inReads.emplace_back(readObject(seqInfo(row[0], row[0])));
				}
			}
		}
	} else {
		kmers = tokenizeString(kmersStrs, ",");
		for (const auto & str : kmers) {
			inReads.emplace_back(readObject(seqInfo(str, str)));
		}
	}
	strCounter counter;
	for (const auto & read : inReads) {
		counter.increaseCountByString(read.seqBase_.seq_, read.seqBase_.cnt_);
	}
	counter.setFractions();
	njh::randomGenerator gen;
	VecStr rSeqs = simulation::randStrs(length, counter, gen, num);
	std::vector<readObject> rReads;
	uint32_t count = 0;
	for (const auto & rSeq : rSeqs) {
		rReads.emplace_back(
				readObject(seqInfo("seq." + estd::to_string(count), rSeq)));
		++count;
	}
	for (const auto & rRead : rReads) {
		if (ansi) {
			rRead.seqBase_.outPutSeqAnsi(std::cout);
		} else {
			rRead.seqBase_.outPutSeq(std::cout);
		}
	}
	return 0;
}

int simulatorRunner::randomStrings(const njh::progutils::CmdArgs & inputCommands) {
  simulatorSetUp setUp(inputCommands);
  uint32_t stringLength = 30;
  uint32_t stringNumber = 30;
  std::string logFilename = "";
  bool fastaOut = false;
  bool useHpRuns = false;

  setUp.setOption(fastaOut, "-fastaOut", "fastaOut");
  setUp.setOption(stringLength, "-stringLength,-len", "stringLength");
  setUp.setOption(stringNumber, "-stringNumber,-num", "stringNumber");
  setUp.setOption(logFilename, "-log,-logFileName", "logFileName");
  bool usePositionCounts = false;
  if(setUp.processDefaultReader(false)){
  	if(!setUp.setOption(usePositionCounts, "-usePositionCounts,-usePosition,-position", "usePositionCounts")){
  		setUp.setOption(useHpRuns, "-useHpRuns,-hRuns", "useHrRuns");
  	}
  }
  bool usingInputAlphabet = false;
  std::string alphabetString = "A,C,G,T";
  if(setUp.setOption(alphabetString, "-alph,-alphabet", "Alphabet")){
  	usingInputAlphabet = true;
  }
  setUp.finishSetUp(std::cout);
  njh::randomGenerator gen;
  std::vector<charCounter> counts;
  VecStr randoms;
  if(setUp.pars_.ioOptions_.firstName_ != ""){
  	SeqInput reader(setUp.pars_.ioOptions_);
  	reader.openIn();
  	auto inReads = reader.readAllReads<readObject>();
  	auto iden = clusterCollapser::collapseIdenticalReads(inReads, "median");
  	if(usePositionCounts){
  		uint64_t maxLen = 0;
  		readVec::getMaxLength(inReads, maxLen);
  		counts = std::vector<charCounter>(maxLen);
  		for(const auto & read : iden){
  			for(const auto pos : iter::range(read.seqBase_.seq_.size())){
  				counts[pos].increaseCountOfBase(read.seqBase_.seq_[pos], read.seqBase_.cnt_);
  			}
  		}
      for (uint32_t num = 0; num < stringNumber; ++num) {
        randoms.emplace_back(simulation::randStr<charCounter, char>(counts, gen));
      }
  	}else if(useHpRuns){
  		hrCounter counter;
  		readVec::allSetCondensedSeq(iden);
  		counter.inceaseCountByReads(iden);
  		counter.setFractions();
  		randoms = simulation::randStrs(stringLength, counter, gen, stringNumber);
  	}else {
    	charCounter counter;
    	for(const auto & read : iden){
    		counter.increaseCountByString(read.seqBase_.seq_, read.seqBase_.cnt_);
    	}
    	counter.resetAlphabet(false);
    	randoms = simulation::randStrs(stringLength, counter, gen, stringNumber);
  	}
  } else if (usingInputAlphabet) {
    VecStr alphbetStrs = processAlphStrVecStr(alphabetString, ",");
    charCounter counter;
    for(const auto & let : alphbetStrs){
    	if(let.size()> 1){
    		counter.increaseCountOfBase(let[0], std::stof(let.substr(1)));
    	}else{
    		counter.increaseCountOfBase(let[0]);
    	}
    }
    counter.resetAlphabet(false);
    counter.setFractions();
    //counter.outPutACGTFractionInfo(std::cout);
    randoms = simulation::randStrs(stringLength, counter, gen, stringNumber);
  } else {
    for (const auto i : iter::range<uint32_t>(0, stringLength)) {
    	charCounter testCounter =
    			charCounter(std::vector<char>{'A', 'G', 'T', 'C'});
      testCounter.chars_['A'] = 1 + i;
      testCounter.chars_['G'] = 1 + i;
      testCounter.chars_['C'] = stringLength + 1 - i;
      testCounter.chars_['T'] = stringLength + 1 - i;
      testCounter.setFractions();
      counts.emplace_back(testCounter);
    }
    {
      //simulation::TicToc t("testing", logFilename);
      for (uint32_t num = 0; num < stringNumber; ++num) {
        randoms.emplace_back(simulation::randStr<charCounter, char>(counts, gen));
      }
    }
  }
  if(fastaOut){
  	uint32_t count = 0;
  	for(const auto & str : randoms){
  		std::cout << ">Seq." << leftPadNumStr<uint32_t>(count, randoms.size()) << std::endl;
  		std::cout << str << std::endl;
  		++count;
  	}
  }else{
  	printVector(randoms, "\n");
  }

  return 0;
}


int simulatorRunner::randomSeqFile(const njh::progutils::CmdArgs & inputCommands) {
  simulatorSetUp setUp(inputCommands);
  uint32_t qualStart = 10;
  uint32_t qualStop = 40;
  std::string alphabetStr = "A,C,G,T";
  uint32_t len = 100;
  uint32_t lenStop = 100;
  uint32_t seqNum = 10;
  uint32_t topAmount = 100;
  uint32_t bottomAmount = 50;
  bool processed = false;
  bool fasta = false;
  uint32_t width = 50;
  setUp.setOption(width, "-width","Line length for when outputting fasta");
  setUp.setOption(fasta, "-fasta", "Output fasta instead of fastq");
  setUp.setOption(qualStart, "-qualStart", "qualStart");
  setUp.setOption(qualStop, "-qualStop", "qualStop");
  setUp.setOption(len, "-len", "len");
  setUp.setOption(lenStop, "-lenStop", "lenStop");
  setUp.setOption(seqNum, "-seqNum", "seqNum");
  setUp.setOption(processed, "-processed", "processed");
  setUp.setOption(topAmount, "-topAmount", "topAmount");
  setUp.setOption(bottomAmount, "-bottomAmount", "bottomAmount");
  setUp.setOption(alphabetStr, "-alphabet","alphabetStr");
  setUp.finishSetUp(std::cout);
  auto alphabetCounts = processAlphStrVecCharCounts(alphabetStr, ",");
  if(setUp.pars_.verbose_){
  	printVector(alphabetCounts.first, ",");
  	printVector(alphabetCounts.second, ",");
  }

  randomFileCreator fileGen(alphabetCounts.first, alphabetCounts.second, qualStart, qualStop);
  fileGen.randomFile(len, lenStop, seqNum, processed, bottomAmount, topAmount,!fasta, std::cout	);
  return 0;
}




}  // namespace njh
