
//  seqUtilsInfoRunner.cpp
//
//  Created by Nicholas Hathaway on 2015/05/28.
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
    
#include "seqUtilsInfoRunner.hpp"
#include "elucidator/simulation.h"
#include "elucidator/objects/counters/DNABaseCounter.hpp"
#include <njhseq/helpers.h>

    
namespace njhseq {

seqUtilsInfoRunner::seqUtilsInfoRunner()
    : njh::progutils::ProgramRunner({
	addFunc("getHpProfile", getHpProfile, false),
  addFunc("profileReadsToReference", profileReadsToReference, false),
  addFunc("findSeq", findSeq, false),
  addFunc("countSeqs", countSeqsFile, false),
  addFunc("fastaIdenticalInfo", fastaIdenticalInfo, true),
  addFunc("uniqueSeqInfo", fastaIdenticalInfo, false),
  addFunc("countLetters", countLetters, false),
  addFunc("countSeqPortion", countSeqPortion, false),
  addFunc("printTandems", printTandems, false),
  addFunc("countOtus", countOtus, false),
  addFunc("quickLenInfo", quickLenInfo, false),
  addFunc("quickMismatchDist", quickMismatchDist, false),
  addFunc("profileQualityScores", profileQualityScores, false),
  addFunc("countHPRuns", countHPRuns, false),
	addFunc("countKmersPlusStats", countKmersPlusStats, false),
	addFunc("profileErrors", profileErrors, false),
	addFunc("countAllSeqs", countAllSeqs, false),
  addFunc("printNames", printNames, false),
	addFunc("fracInfo", fracInfo, false),
	addFunc("genPsuedoAllMinTree", genPsuedoAllMinTree, false),
	addFunc("genPsuedoMismatchMinTree", genPsuedoMismatchMinTree, false),
	addFunc("qualCounts", qualCounts, false),
	addFunc("printSeqs", printSeqs, false),
	addFunc("quickHaplotypeInformation", quickHaplotypeInformation, false),
	addFunc("oldQuickHaplotypeInformationAndVariants", oldQuickHaplotypeInformationAndVariants, false),
	addFunc("getGCContent", getGCContent, false),
	addFunc("quickHaplotypeVariantsWithRegion", quickHaplotypeVariantsWithRegion, false),
	addFunc("readLengthDistribution", readLengthDistribution, false),
	addFunc("getReadLens", getReadLens, false),
	addFunc("multipleAlnProteinToPcaInput", multipleAlnProteinToPcaInput, false),
	addFunc("getHapPopDifAndVariantsInfo", getHapPopDifAndVariantsInfo, false),
	addFunc("getSlidingEntropy", getSlidingEntropy, false),
	},
                    "seqUtilsInfo") {}
//
//



int seqUtilsInfoRunner::multipleAlnProteinToPcaInput(const njh::progutils::CmdArgs & inputCommands) {
	bool setGapsAsZero = 0;
	seqSetUp setUp(inputCommands);
	setUp.setOption(setGapsAsZero, "--setGapsAsZero", "Set Gaps As Zero");
	setUp.processDefaultReader(true);
	setUp.pars_.ioOptions_.out_.outExtention_ = ".txt";
	setUp.finishSetUp(std::cout);

	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<readObject>();
  //auto reads = createKmerReadVec(inReads,setUp.pars_.colOpts_.kmerOpts_.kLength_, true);

	if(inReads.size() < 2){
		return 0;
	}
	for (const auto pos : iter::range(inReads.size() - 1)) {
		if (len(inReads[pos].seqBase_)
				!= len(inReads[pos + 1].seqBase_)) {
			std::stringstream ss;

			ss << inReads[pos].seqBase_.name_
					<< " is not the same size as " << inReads[pos + 1].seqBase_.name_
					<< std::endl;
			ss << inReads[pos].seqBase_.name_ << ": " << len(inReads[pos].seqBase_) << std::endl;
			ss << inReads[pos + 1].seqBase_.name_ << ": " << len(inReads[pos + 1].seqBase_) << std::endl;
			throw std::runtime_error{ss.str()};
		}
	}

	std::vector<charCounter> counters(len(inReads.front().seqBase_));
	for(const auto & read : inReads){
		for(const auto pos : iter::range(len(read))){
			counters[pos].increaseCountOfBase(read.seqBase_.seq_[pos]);
		}
	}

	for(auto & counter : counters){
		counter.resetAlphabet(false);
		counter.setFractions();
	}
	std::ofstream testFile;
	openTextFile(testFile, setUp.pars_.ioOptions_.out_);
	std::vector<std::map<char, uint32_t>> rankConverters;

	for(auto pos : iter::range(counters.size())){
		std::map<uint32_t, std::map<char, uint32_t>> ranks;
		if(setUp.pars_.verbose_){
			std::cout << njh::bashCT::boldBlack(estd::to_string(pos)) << std::endl;
		}
		for(const auto & aa : counters[pos].alphabet_){
			//gaps will be coded as zero
			if('-' == aa && setGapsAsZero){
				continue;
			}
			ranks[counters[pos].chars_[aa]][aa] = counters[pos].chars_[aa];
		}
		uint32_t count = 1;
		std::map<char, uint32_t> rankConver;
		for(const auto & occ : ranks){
			for(const auto & aa : occ.second){
				if(setUp.pars_.verbose_){
					std::cout << "\t" << aa.first << occ.first << std::endl;
				}
				rankConver[aa.first] = count;
				++count;
			}
		}
		for(const auto & rank : rankConver){
			if(setUp.pars_.verbose_){
				std::cout << "\t" << rank.first << ":" << rank.second << std::endl;
			}
		}
		rankConverters.emplace_back(rankConver);
	}
	std::vector<uint32_t> positions(counters.size());
	njh::iota<uint32_t>(positions, 0);
	testFile << "name\t" << vectorToString(positions, "\t") << std::endl;
	for(const auto & read: inReads){
		testFile << read.seqBase_.name_;
		for(const auto & pos : positions){
			//gaps are zero
			testFile << "\t" << ('-' == read.seqBase_.seq_[pos] && setGapsAsZero ? 0: rankConverters[pos][read.seqBase_.seq_[pos]]);
		}
		testFile << std::endl;
	}
	return 0;
}



int seqUtilsInfoRunner::readLengthDistribution(const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	bool stats = false;
	setUp.setOption(stats, "--getStats", "Print basic stats as well");
	setUp.processDefaultReader(true);
  setUp.pars_.ioOptions_.out_.outExtention_ = ".tab.txt";
  setUp.pars_.ioOptions_.out_.outFileFormat_ = "tab";
	setUp.finishSetUp(std::cout);

	seqInfo read;
	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	std::unordered_map<size_t,uint32_t> readLengths;
	while(reader.readNextRead(read)){
		readLengths[len(read)] += read.cnt_;
	}
	auto keys = getVectorOfMapKeys(readLengths);
	njh::sort(keys);
  table tab(VecStr{"length", "count"});
	for(auto key : keys){
		tab.content_.emplace_back(toVecStr(key, readLengths[key]));
	}
  TableIOOpts opts(setUp.pars_.ioOptions_.out_, "\t", true);
  tab.outPutContents(opts);
	return 0;
}

int seqUtilsInfoRunner::getReadLens(
		const njh::progutils::CmdArgs & inputCommands) {
	VecStr addHeader;
	OutOptions outOpts(bfs::path(""));
	outOpts.outExtention_ = ".tsv";
	seqSetUp setUp(inputCommands);
	setUp.processReadInNames(true);
	setUp.processWritingOptions(outOpts);
	setUp.setOption(addHeader, "--addHeader", "Add this header, must be two values separated by commas");
	if(!addHeader.empty() && 2 != addHeader.size()){
		setUp.failed_ = true;
		setUp.addWarning(njh::pasteAsStr("addHeader must be size 2, not ", addHeader.size()));
	}
	setUp.finishSetUp(std::cout);

	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	seqInfo read;
	OutputStream out(outOpts);
	if(!addHeader.empty()){
		out << njh::conToStr(addHeader, "\t") << std::endl;
	}
	while (reader.readNextRead(read)) {
		out << read.name_
				<< "\t" << len(read)
				<< std::endl;;
	}
	return 0;
}



int seqUtilsInfoRunner::getGCContent(const njh::progutils::CmdArgs & inputCommands) {
	auto outOpts = TableIOOpts::genTabFileOut(bfs::path(""), true);
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processWritingOptions(outOpts.out_);
	setUp.processReadInNames(true);
	setUp.finishSetUp(std::cout);
	seqInfo seq;

	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	OutputStream out(outOpts.out_);
	out << "name\tgcContent" << std::endl;
	while(reader.readNextRead(seq)){
		DNABaseCounter counter;
		counter.increase(seq.seq_);
		out << seq.name_
				<< "\t" <<counter.calcGcContent()
				<< std::endl;
	}
	return 0;
}
int seqUtilsInfoRunner::qualCounts(const njh::progutils::CmdArgs & inputCommands) {
	auto outOpts = TableIOOpts::genTabFileOut(bfs::path(""), true);
	seqSetUp setUp(inputCommands);
	setUp.processWritingOptions(outOpts.out_);
	setUp.processReadInNames(true);
	bool byBase = false;
	bool byPosition = false;
	setUp.setOption(byBase, "--byBase", "Break Down Qual Counts Per Base");
	setUp.setOption(byPosition, "--byPosition", "Break Down Qual Counts Per Positions");
	setUp.processVerbose();

	setUp.finishSetUp(std::cout);

	if ("STDOUT" != outOpts.out_.outFilename_
			&& "" != outOpts.out_.outFilename_
			&& outOpts.out_.outExists()
			&& !outOpts.out_.overWriteFile_) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << outOpts.out_.outName()
				<< " already exists, use --overWrite to over write" << "\n";
		throw std::runtime_error { ss.str() };
	}

	if(byBase && byPosition){
		std::unordered_map<uint32_t,std::unordered_map<char, std::unordered_map<uint32_t, uint32_t>>> qualCounts;
		SeqInput reader(setUp.pars_.ioOptions_);
		reader.openIn();
		seqInfo seq;
		uint32_t count = 0;
		while(reader.readNextRead(seq)){
			++count;
			if(setUp.pars_.verbose_ && count % 1000 == 0){
				std::cout << "\r" << count;
			}
			for(const auto pos : iter::range(len(seq))){
				qualCounts[pos][seq.seq_[pos]][seq.qual_[pos]] += seq.cnt_;
			}
		}
		if(setUp.pars_.verbose_ ){
			std::cout << std::endl;
		}
		table outTab(VecStr{"pos", "char", "qual", "count"});
		for(const auto & pos : qualCounts){
			for(const auto & base : pos.second){
				for(const auto & qual : base.second){
					outTab.addRow(pos.first, base.first, qual.first, qual.second);
				}
			}
		}
		outTab.sortTable("pos","char", "qual",   false);
		outTab.outPutContents(outOpts);
	}else if(byBase){
		std::unordered_map<char, std::unordered_map<uint32_t, uint32_t>> qualCounts;
		SeqInput reader(setUp.pars_.ioOptions_);
		reader.openIn();
		seqInfo seq;
		uint32_t count = 0;
		while(reader.readNextRead(seq)){
			++count;
			if(setUp.pars_.verbose_ && count % 1000 == 0){
				std::cout << "\r" << count;
			}
			for(const auto pos : iter::range(len(seq))){
				qualCounts[seq.seq_[pos]][seq.qual_[pos]] += seq.cnt_;
			}
		}
		if(setUp.pars_.verbose_ ){
			std::cout << std::endl;
		}
		table outTab(VecStr {"char", "qual", "count"});
		for(const auto & base : qualCounts) {
			for(const auto & qual : base.second) {
				outTab.addRow(base.first, qual.first, qual.second);
			}
		}
		outTab.sortTable("char", "qual",  false);
		outTab.outPutContents(outOpts);
	}else if(byPosition){
		std::unordered_map<uint32_t,std::unordered_map<uint32_t, uint32_t>> qualCounts;
		SeqInput reader(setUp.pars_.ioOptions_);
		reader.openIn();
		seqInfo seq;
		uint32_t count = 0;
		while(reader.readNextRead(seq)){
			++count;
			if(setUp.pars_.verbose_ && count % 1000 == 0){
				std::cout << "\r" << count;
			}
			for(const auto pos : iter::range(len(seq))){
				qualCounts[pos][seq.qual_[pos]] += seq.cnt_;
			}
		}
		if(setUp.pars_.verbose_ ){
			std::cout << std::endl;
		}
		table outTab(VecStr {"pos", "qual", "count"});
		for(const auto & pos : qualCounts) {
			for(const auto & qual : pos.second) {
				outTab.addRow(pos.first, qual.first, qual.second);
			}
		}
		outTab.sortTable( "pos", "qual", false);
		outTab.outPutContents(outOpts);
	} else {
		std::unordered_map<uint32_t, uint32_t> qualCounts;
		SeqInput reader(setUp.pars_.ioOptions_);
		reader.openIn();
		seqInfo seq;
		uint32_t count = 0;
		while(reader.readNextRead(seq)){
			++count;
			if(setUp.pars_.verbose_ && count % 1000 == 0){
				std::cout << "\r" << count;
			}
			for(const auto pos : iter::range(len(seq))){
				qualCounts[seq.qual_[pos]] += seq.cnt_;
			}
		}
		if(setUp.pars_.verbose_ ){
			std::cout << std::endl;
		}
		table outTab(VecStr {"qual", "count"});
		for(const auto & qual : qualCounts) {
			outTab.addRow(qual.first, qual.second);
		}
		outTab.sortTable("qual", false);
		outTab.outPutContents(outOpts);
	}
	return 0;
}
                    
int seqUtilsInfoRunner::countAllSeqs(const njh::progutils::CmdArgs & inputCommands){
	seqSetUp setUp(inputCommands);
	std::string ending = "";
	std::string format = "";
	bool processed = false;
	bool recursive = false;
	std::string dir = "./";
	OutOptions outOpts;
	setUp.setOption(ending, "--ending", "Ending of Files to Count Seqs", true);
	if(!setUp.setOption(format, "--format", "Format")){
		format = ending;
	}
	setUp.processWritingOptions(outOpts);
	setUp.setOption(processed, "--processed", "If the Read Name Contains Read Cnt Info");
	setUp.setOption(dir, "--dir", "dir");
	setUp.setOption(recursive, "--recursive", "recursive");
  setUp.finishSetUp(std::cout);
  uint64_t totalCount = 0;
  	OutputStream out(outOpts);
  auto files = njh::files::listAllFiles(dir, recursive, {std::regex{".*" + ending + "$"}});

  out << "filename\tcount" << "\n";
  for(const auto & f : files){
  		uint32_t currentCount = 0;
  		if(!f.second){
  			SeqInput reader(SeqIOOptions(f.first, SeqIOOptions::getInFormat(format), processed));
  			reader.openIn();
  			seqInfo seq;
  			while(reader.readNextRead(seq)){
  				++currentCount;
  			}
  		}
  		out << f.first << "\t" << currentCount << "\n";
  		totalCount += currentCount;
  }
  out << "total\t" << totalCount << "\n";
  return 0;
}



int seqUtilsInfoRunner::profileErrors(const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	setUp.processDefaultReader(true);
	setUp.processSeq(true);
	setUp.processDirectoryOutputName(true);
  setUp.processAlignerDefualts();
  setUp.finishSetUp(std::cout);
  setUp.startARunLog(setUp.pars_.directoryName_);
  //profileKmerAccerlation
  SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<readObject>();

	uint64_t maxReadLen = 0;
	readVec::getMaxLength(inReads, maxReadLen);
	readVec::getMaxLength(setUp.pars_.seqObj_, maxReadLen);
	KmerMaps kMaps = indexKmers(inReads, setUp.pars_.colOpts_.kmerOpts_.kLength_,
			setUp.pars_.colOpts_.kmerOpts_.runCutOff_,
			setUp.pars_.colOpts_.kmerOpts_.kmersByPosition_,
			setUp.pars_.expandKmerPos_, setUp.pars_.expandKmerSize_);
	aligner alignerObj(maxReadLen, setUp.pars_.gapInfo_, setUp.pars_.scoring_,
			kMaps, setUp.pars_.qScorePars_, setUp.pars_.colOpts_.alignOpts_.countEndGaps_,
			setUp.pars_.colOpts_.iTOpts_.weighHomopolyer_);

  alignerObj.processAlnInfoInput(setUp.pars_.alnInfoDirName_);
  std::ofstream misInfo;
  openTextFile(misInfo, setUp.pars_.directoryName_ + "misInfo.tab.txt",
  		".txt", setUp.pars_.ioOptions_.out_);
  std::ofstream disInfo;
  openTextFile(disInfo, setUp.pars_.directoryName_ + "disInfo.tab.txt",
  		".txt", setUp.pars_.ioOptions_.out_);
  std::map<uint32_t, uint32_t> qualCounts;
  readVec::updateQaulCountsMultiple(inReads, qualCounts);
  misInfo << "readName\trefBase\trefBasePos\tseqBase\tseqBasePos\tseqQual\tqualTotal\tleadingQual\ttrailingQual\tmisType\tkMerFreq\tkMerFreqByPos\n";
  disInfo << "readName\tpercentIdentity\tpercentGaps\tseqCoverage\trefCoverage\teventBasedIdentity\n";
  for(const auto & read : inReads){
  	alignerObj.alignCache(setUp.pars_.seqObj_.seqBase_, read.seqBase_,setUp.pars_.local_);
  	alignerObj.profilePrimerAlignment(setUp.pars_.seqObj_.seqBase_, read.seqBase_);
  	disInfo << read.seqBase_.name_
  			<< '\t' << alignerObj.comp_.distances_.percentMatch_
  			<< '\t' << alignerObj.comp_.distances_.percentGaps_
  			<< '\t' << alignerObj.comp_.distances_.query_.coverage_
  			<< '\t' << alignerObj.comp_.distances_.ref_.coverage_
  			<< '\t' << alignerObj.comp_.distances_.eventBasedIdentity_
  			<< '\n';
  	for(const auto & m : alignerObj.comp_.distances_.mismatches_){
  		std::string misType =  m.second.transition ? "transition" : "transversion";
  		misInfo << read.seqBase_.name_ << "\t" << m.second.refBase <<
  				"\t" << m.second.refBasePos << "\t" << m.second.seqBase <<
  				"\t" << m.second.seqBasePos << "\t" << m.second.seqQual <<
  				"\t" << qualCounts[m.second.seqQual]
  				<< "\t" << vectorToString(m.second.seqLeadingQual, ",")
  				<< "\t" << vectorToString(m.second.seqTrailingQual, ",")
  				<< "\t" << misType
  				<< "\t" << m.second.kMerFreq <<
  				"\t" << m.second.kMerFreqByPos << "\n";
  	}
  }
  if (setUp.pars_.writingOutAlnInfo_) {
    alignerObj.alnHolder_.write(setUp.pars_.outAlnInfoDirName_);
  }
  return 0;
}


int seqUtilsInfoRunner::countKmersPlusStats(const njh::progutils::CmdArgs & inputCommands){
	OutOptions outOpts;
	uint32_t kLen = 10;
	bool noHeader = false;
	seqSetUp setUp(inputCommands);
	setUp.setOption(noHeader, "--noHeader", "Don't print a header as well");
	setUp.setOption(kLen, "--kmerLength", "kmer Length");
	setUp.processReadInNames(true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	std::unordered_map<std::string, uint32_t> kmersCounts;
	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	charCounter counter;
	seqInfo seq;
	OutputStream out(outOpts);
	while(reader.readNextRead(seq)){
		counter.increaseCountByString(seq.seq_, seq.cnt_);
		for(auto pos : iter::range(seq.seq_.size() - kLen + 1)){
			kmersCounts[seq.seq_.substr(pos,kLen)]+= round(seq.cnt_);
		}
	}

	auto keys = getVectorOfMapKeys(kmersCounts);
	auto values = getVectorOfMapValues(kmersCounts);
	auto kSums = vectorSum(values);
	counter.resetAlphabet(false);
	counter.setFractions();

	njh::sort(keys);
	if(!noHeader){
		out << "kmer\tcount\tfraction\tchance\tobs/expected\n";
	}
	for(const auto & k : keys){
		double chance = 1;
		for(auto c : k){
			chance *=counter.fractions_[c];
		}
		out << k
				<< "\t" << kmersCounts[k]
				<< "\t" << kmersCounts[k]/static_cast<double>(kSums)
				<< "\t" << chance
				<< "\t" << kmersCounts[k]/(chance * kSums)
				<< "\n";
	}
	return 0;
}

int seqUtilsInfoRunner::countHPRuns(const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	bool plot = false;
	setUp.pars_.ioOptions_.out_.outFilename_ = "";
	setUp.setOption(plot, "-plot", "plot");

	setUp.processDefaultReader(true);
	setUp.finishSetUp(std::cout);
	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<readObject>();

	//std::unordered_map<char, std::unordered_map<uint32_t, uint32_t>> hCounts;
	readVec::allSetCondensedSeq(inReads);
	hrCounter counter;
	counter.inceaseCountByReads(inReads);
	table out(counter.hCounts_, VecStr{"char", "hpRunSize", "count"});
	counter.setFractions();
	table outFrac(counter.fractions_, VecStr{"char", "hpRunSize", "fractions"});
	out.sortTable("hpRunSize", true);
	out.sortTable("char", true);
	if(setUp.pars_.ioOptions_.out_.outFilename_ == ""){
		out.outPutContentOrganized(std::cout);
		outFrac.outPutContentOrganized(std::cout);
	}else{
		TableIOOpts outOptCount = TableIOOpts(OutOptions(bfs::path(bfs::basename(setUp.pars_.ioOptions_.out_.outFilename_) + "_count"),
				".tab.txt", "tab", false, setUp.pars_.ioOptions_.out_.overWriteFile_, false),"\t", out.hasHeader_);
		out.outPutContents(outOptCount);
		TableIOOpts outOptFrac = TableIOOpts(OutOptions(bfs::path(bfs::basename(setUp.pars_.ioOptions_.out_.outFilename_) + "_frac"),
				".tab.txt", "tab",false, setUp.pars_.ioOptions_.out_.overWriteFile_, false),"\t", outFrac.hasHeader_);
		outFrac.outPutContents(outOptFrac);
	}
//	if(plot){
//
//		njhRInside::OwnRInside rSes;
//		auto & R = rSes.get();
//		rSes.loadLib("ggplot2");
//		rSes.loadLib("colorspace");
//		rSes.openPdfDevice("hpCounts",11, 8.5);
//		{
//			Rcpp::DataFrame hpFrame = Rcpp::DataFrame::create(Rcpp::Named("char") = out.getColumn("char"),
//					Rcpp::Named("hpRunSize") = out.getColumn("hpRunSize"),
//					Rcpp::Named("count") = out.getColumn("count"));
//			R["hpCounts"] = hpFrame;
//			VecStr coms(3);
//			coms[0] = "hpCounts$char = as.factor(as.character(hpCounts$char))";
//			coms[1] = "hpCounts$hpRunSize = as.numeric(as.character(hpCounts$hpRunSize))";
//			coms[2] = "hpCounts$count = as.numeric(as.character(hpCounts$count))";
//			rSes.multipleParseEvalQ(coms);
//			R.parseEvalQ("print(ggplot(hpCounts, aes(x = hpRunSize, y = count, fill = as.factor(hpRunSize)) ) +geom_bar(stat = \"identity\")+ guides(fill=guide_legend(title=\"HP_Run_Size\")) + scale_fill_manual(values = heat_hcl(length(levels(as.factor(hpCounts$hpRunSize))), h = c(160, 420), c. = c(100, 100), l = c(40,70)))  + facet_grid(~char))");
//		}
//		{
//			Rcpp::DataFrame hpFrame = Rcpp::DataFrame::create(Rcpp::Named("char") = outFrac.getColumn("char"),
//					Rcpp::Named("hpRunSize") = outFrac.getColumn("hpRunSize"),
//					Rcpp::Named("count") = outFrac.getColumn("fractions"));
//			R["hpCounts"] = hpFrame;
//			VecStr coms(3);
//			coms[0] = "hpCounts$char = as.factor(as.character(hpCounts$char))";
//			coms[1] = "hpCounts$hpRunSize = as.numeric(as.character(hpCounts$hpRunSize))";
//			coms[2] = "hpCounts$count = as.numeric(as.character(hpCounts$count))";
//			rSes.multipleParseEvalQ(coms);
//			R.parseEvalQ("print(ggplot(hpCounts, aes(x = hpRunSize, y = count, fill = as.factor(hpRunSize)) ) +geom_bar(stat = \"identity\") + guides(fill=guide_legend(title=\"HP_Run_Size\")) + scale_fill_manual(values = heat_hcl(length(levels(as.factor(hpCounts$hpRunSize))), h = c(160, 420), c. = c(100, 100), l = c(40,70))) + facet_grid(~char))");
//		}
//		rSes.closeCurrentDevice();
//	}
	return 0;
}


int seqUtilsInfoRunner::profileQualityScores(const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	setUp.processDefaultReader(true);
	std::string qualWindow;
	uint32_t qualityWindowLength;
	uint32_t qualityWindowStep;
	uint32_t qualityWindowThres;
  if (setUp.setOption(qualWindow, "-qualWindow", "SlidingQualityWindow")) {
    seqUtil::processQualityWindowString(qualWindow, qualityWindowLength,
                                        qualityWindowStep, qualityWindowThres);
  } else {
    qualityWindowLength = 50;
    qualityWindowStep = 5;
    qualityWindowThres = 25;
  }
  uint32_t qualCheck = 30;
  setUp.setOption(qualCheck, "-qualCheck", "-qualCheck"	);
  setUp.processDirectoryOutputName(true);
	setUp.finishSetUp(std::cout);

	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<readObject>();


	std::unordered_map<uint32_t, uint32_t> qualWindowCounts;
	uint32_t count = 0;
	uint32_t tenPer = 0.1 * len(inReads);
	if(tenPer == 0){
		tenPer = 1;
	}
	readVec::allSetQualCheck(inReads, qualCheck);
	std::vector<double> qualChecks;
	for(const auto & read : inReads){
		if(count % tenPer == 0){
			std::cout << "Currently on " << count << " of " << len(inReads) << std::endl;
		}
		qualChecks.emplace_back(read.fractionAboveQualCheck_);
		++count;
	  uint32_t currentPos = 0;
	  while ((qualityWindowLength + currentPos) < read.seqBase_.qual_.size()) {
	    uint32_t sum = 0;
	    for (const auto qPos : iter::range(currentPos, currentPos + qualityWindowLength)) {
	      sum += read.seqBase_.qual_[qPos];
	    }
	    qualWindowCounts[sum]+=::round(read.seqBase_.cnt_);
	    currentPos += qualityWindowStep;
	  }
	}

	table out(VecStr{"qualWindowAverage", "count"});
	for(const auto & qCount : qualWindowCounts){
		out.content_.emplace_back(VecStr{estd::to_string(static_cast<double>(qCount.first)/qualityWindowLength),
				estd::to_string(qCount.second)});
	}
	out.sortTable("qualWindowAverage", true);
	TableIOOpts qualWindowOption(OutOptions(setUp.pars_.directoryName_ + "qualWindowAverage", ".tab.txt"), "\t", out.hasHeader_);
	out.outPutContents(qualWindowOption);
	std::map<uint32_t, uint32_t> qualCounts;
	readVec::updateQaulCountsMultiple(inReads, qualCounts);
	table qualCountsTab(qualCounts, VecStr{"qual", "counts"});
	out.outPutContentOrganized(std::cout);
	qualCountsTab.outPutContentOrganized(std::cout);
	TableIOOpts qualCountOption(OutOptions(setUp.pars_.directoryName_ + "qualCounts", ".tab.txt"), "\t", qualCountsTab.hasHeader_);
	qualCountsTab.outPutContents(qualCountOption);
	auto qualCheckStats = getStatsOnVec(qualChecks);
	table qualChecksTab(qualCheckStats, VecStr{"stat", "value"});
	qualChecksTab.outPutContentOrganized(std::cout);
	TableIOOpts qualChecksOpts(OutOptions(setUp.pars_.directoryName_ + "qualChecks", ".tab.txt"), "\t", qualChecksTab.hasHeader_);
	qualChecksTab.outPutContents(qualChecksOpts);
	return 0;
}
int seqUtilsInfoRunner::countSeqsFile(const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	setUp.processDefaultReader(true);
	setUp.processVerbose();
	setUp.finishSetUp(std::cout);
	uint32_t readAmount = countSeqs(setUp.pars_.ioOptions_, setUp.pars_.verbose_);
	std::cout << readAmount << std::endl;
	return 0;
}

int seqUtilsInfoRunner::printTandems(
		const njh::progutils::CmdArgs & inputCommands) {
	OutOptions outOpts;
	int32_t match = 2;
	int32_t mismatch = -2;
	int32_t gap = -7;
	int32_t minimumAlignScore = 50;
	seqSetUp setUp(inputCommands);
	setUp.setOption(match, "--matchScore", "Match score");
	setUp.setOption(mismatch, "--mismatchScore", "mismatch score");
	setUp.setOption(gap, "--gapScore", "gap score");
	setUp.setOption(minimumAlignScore, "--minimumAlignScore",
			"Minimum Align Score");
	setUp.processAlignerDefualts();
	if (!setUp.processSeq(false)) {
		setUp.processReadInNames(true);
	}
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);
	OutputStream out(outOpts);
	TandemRepeat::outPutInfoFormatedHeader(out, "\t");
	if ("" == setUp.pars_.seq_) {
		SeqInput reader(setUp.pars_.ioOptions_);
		reader.openIn();
		seqInfo seq;
		while (reader.readNextRead(seq)) {
			auto currentRepeats = aligner::findTandemRepeatsInSequence(seq.seq_,
					match, mismatch, gap, minimumAlignScore);
			for (const auto & rep : currentRepeats) {
				rep.outPutInfoFormated(out, seq.name_, "\t");
			}
		}
	} else {
		auto currentRepeats = aligner::findTandemRepeatsInSequence(setUp.pars_.seq_,
				match, mismatch, gap, minimumAlignScore);
		for (const auto & rep : currentRepeats) {
			rep.outPutInfoFormated(out, "seq", "\t");
		}
	}
	return 0;
}

int seqUtilsInfoRunner::countSeqPortion(const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	bool back = false;
	uint32_t size = 10;
	uint32_t occurenceCutOff = 1;
	uint32_t position = 0;
	setUp.processDefaultReader(true);
	setUp.setOption(size, "--size", "Size of the portion to count");
	setUp.setOption(occurenceCutOff, "--occurenceCutOff",
			"occurenceCutOff");
	setUp.setOption(back, "--back", "Count portion from the back of the sequence");
	setUp.setOption(position, "--position", "Position from which to count");
	setUp.finishSetUp(std::cout);
	std::function<std::string(const seqInfo &, size_t, uint32_t)> getSubStr;
	if (back) {
		if (position != 0 && size > position) {
			std::stringstream ss;
			ss << "Error, size must be less than or equal to position when counting from the back of the sequence"
					<< std::endl;
			ss << "position: " << position << ", size:" << size << std::endl;
			throw std::runtime_error { njh::bashCT::boldRed(ss.str()) };
		}
		getSubStr = [](const seqInfo & read, size_t position, uint32_t size){
			std::string ret = "";
			if(position < read.seq_.size()){
				if(position == 0){
					ret = read.seq_.substr(read.seq_.size() - size, size);
				}else{
					ret = read.seq_.substr(read.seq_.size() - position, size);
				}
			}
			return ret;
		};
	} else {
		getSubStr = [](const seqInfo & read, size_t position, uint32_t size){
			std::string ret = "";
			if(position + size <= read.seq_.size()){
				ret = read.seq_.substr(position, size);
			}
			return ret;
		};
	}

	if(setUp.pars_.ioOptions_.isPairedIn()){
		SeqIO reader(setUp.pars_.ioOptions_);
		reader.openIn();
		PairedRead seq;
		std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>> counts;
		double total = 0;
		while(reader.readNextRead(seq)){
			counts[getSubStr(seq.seqBase_, position, size)][getSubStr(seq.mateSeqBase_, position, size)] += seq.seqBase_.cnt_;
			total+= seq.seqBase_.cnt_;
		}

		table outTable(VecStr { "r1_str", "r2_str", "count", "fraction" });
		for (const auto & r1Str : counts) {
			for(const auto & r2Str : r1Str.second){
				if(r2Str.second > occurenceCutOff){
					outTable.addRow(r1Str.first, r2Str.first, r2Str.second, r2Str.second/total);
				}
			}
		}
		outTable.sortTable("count", true);
		outTable.outPutContentOrganized(std::cout);
	}else{
		SeqIO reader(setUp.pars_.ioOptions_);
		reader.openIn();
		seqInfo seq;
		strCounter counter;
		while(reader.readNextRead(seq)){
			counter.increaseCountByString(getSubStr(seq, position, size));
		}
		counter.setFractions();
		table outTable(VecStr { "str", "count", "fraction" });
		for (const auto & str : counter.counts_) {
			if (str.second > occurenceCutOff && str.first != "") {
				outTable.content_.emplace_back(
						toVecStr(str.first, str.second, counter.fractions_[str.first]));
			}
		}
		outTable.sortTable("count", true);
		outTable.outPutContentOrganized(std::cout);
	}

	return 0;
}
int seqUtilsInfoRunner::fastaIdenticalInfo(const njh::progutils::CmdArgs & inputCommands) {
  std::string qualRep = "median";
  int qualCheck = 30;

  seqUtilsInfoSetUp setUp (inputCommands);
  setUp.setUpFastaIdenticalInfo(qualRep, qualCheck);
  setUp.startARunLog(setUp.pars_.directoryName_);
  // read in sequences
  SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<readObject>();
  std::vector<readObject> inputReads = inReads;
  std::cout << "Read " << inputReads.size() << " reads" << std::endl;
  uint64_t maxReadLength = 0;
  readVec::getMaxLength(inputReads, maxReadLength);
  // info
  std::string info = "info.txt";
  auto inputVector = clusterCollapser::collapseIdenticalReads(
      inputReads, qualRep);
  readVecSorter::sort(inputVector);
  renameReadNames(inputVector, bfs::basename(setUp.pars_.ioOptions_.firstName_), true,
                  true);
	cluster::getInfo(inputVector, setUp.pars_.directoryName_, info, qualCheck);
	if (setUp.pars_.refIoOptions_.firstName_ != "") {
		std::vector<readObject> expectSeq = SeqInput::getReferenceSeq(
				setUp.pars_.refIoOptions_, maxReadLength);
		// aligner object
		auto scoringMatrixMap = substituteMatrix::createDegenScoreMatrix(1, -1);
		KmerMaps emptyMaps;
		gapScoringParameters gapPars(7, 1, 0, 0, 0, 0);
		aligner alignerObj(maxReadLength, gapPars, scoringMatrixMap, emptyMaps,
				setUp.pars_.qScorePars_, setUp.pars_.colOpts_.alignOpts_.countEndGaps_,
				setUp.pars_.colOpts_.iTOpts_.weighHomopolyer_);
		profiler::getFractionInfo(inputVector, setUp.pars_.directoryName_,
				"readsInfo", setUp.pars_.refIoOptions_.firstName_.string(), alignerObj, false);
	}
	setUp.logRunTime(std::cout);
	return 0;
}



int seqUtilsInfoRunner::countLetters(const njh::progutils::CmdArgs & inputCommands) {
	// parameters
	seqUtilsInfoSetUp setUp(inputCommands);
	if (setUp.needsHelp()) {
		std::cout << "countLetters" << std::endl;
		std::cout << "Required Options, order not necessary" << std::endl;
		//setUp.printInputUsage(std::cout);
		exit(1);
	}
	if (!setUp.processSeq(false)) {
		setUp.processDefaultReader(true);
	}
	bool resetAlph = false;
	bool keepDNA = false;
	std::string alphabet = "A,C,G,T";
	if (!setUp.setOption(alphabet, "-alphabet", "alphabet")) {
		setUp.setOption(resetAlph, "-resetAlphabet", "Reset Alphabet");
		setUp.setOption(keepDNA, "-keepDNA", "keepDNA");
	}
	setUp.finishSetUp(std::cout);
	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<readObject>();
	std::vector<char> alphVec = processAlphStrVecChar(alphabet, ",");
	charCounter counter(alphVec);
	if (setUp.pars_.seq_ == "") {
		SeqIO reader(setUp.pars_.ioOptions_);
		reader.openIn();
		seqInfo seq;
		while(reader.readNextRead(seq)){
			counter.increaseCountByString(seq.seq_);
		}
	} else {
		counter.increaseCountByString(setUp.pars_.seq_);
	}
	if (resetAlph) {
		counter.resetAlphabet(keepDNA);
	}
	counter.setFractions();
	counter.outPutInfo(std::cout, false);
	return 0;
}

int seqUtilsInfoRunner::printNames(const njh::progutils::CmdArgs & inputCommands){
	seqSetUp setUp(inputCommands);
	setUp.pars_.ioOptions_.out_.outFilename_ = "";
	setUp.processDefaultReader(true);
	setUp.pars_.ioOptions_.out_.outExtention_ = ".txt";
	setUp.finishSetUp(std::cout);
	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	seqInfo seq;
	OutputStream out(setUp.pars_.ioOptions_.out_);
	while (reader.readNextRead(seq)) {
		out << seq.name_ << "\n";
	}
	return 0;
}

int seqUtilsInfoRunner::printSeqs(const njh::progutils::CmdArgs & inputCommands){
	seqSetUp setUp(inputCommands);
	setUp.pars_.ioOptions_.out_.outFilename_ = "";
	setUp.processDefaultReader(true);
	setUp.pars_.ioOptions_.out_.outExtention_ = ".txt";
	setUp.finishSetUp(std::cout);
	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	seqInfo seq;
	OutputStream out(setUp.pars_.ioOptions_.out_);
	while (reader.readNextRead(seq)) {
		out << seq.seq_ << "\n";
	}
	return 0;
}


int seqUtilsInfoRunner::findSeq(const njh::progutils::CmdArgs & inputCommands) {
  // remove sequences very dissimilar to input compare sequence
  seqUtilsInfoSetUp setUp(inputCommands);
  std::string compareSeq = "";
  std::string singleSeq = "";
  double queryCutOff = 0.8;
  double percentIdentityCutoff = 0.8;
  double percentageGapsCutOff = 0.1;
  bool debug = false;
  bool extra = false;

  setUp.setOption(extra, "-extra", "extra");
  setUp.setOption(debug, "-debug", "debug");
  if(!setUp.processSeq(singleSeq, "-seq", "optionSequenceToSearch")){
  	setUp.processDefaultReader(true);
  }
  setUp.processSeq(compareSeq, "-target,-compare", "TargetSequenceToSearchFor", true);
  setUp.setOption(percentIdentityCutoff, "-idCutOff,-id", "PercentIdentityCutoff");
  setUp.setOption(queryCutOff, "-queryCutOff", "QueryCutOff");
  setUp.setOption(percentageGapsCutOff, "-gapCutoff,-gcutoff", "GapCutoff ");
  setUp.processAlignerDefualts();
  setUp.processVerbose();
  setUp.finishSetUp(std::cout);
  SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<readObject>();
  readObject compareObject(seqInfo("compare", compareSeq));
  readObject singleObject(seqInfo("singleObject", singleSeq));
  if(setUp.pars_.verbose_){
    std::cout << "QueryCutOff: " << queryCutOff << std::endl;
    std::cout << "IDCutOff: " << percentIdentityCutoff << std::endl;
    std::cout << "GapCutOff: " << percentageGapsCutOff << std::endl;
  }
  uint64_t maxReadLength = 0;
  if (setUp.pars_.ioOptions_.firstName_ != "") {

    readVec::getMaxLength(inReads, maxReadLength);
  }
  readVec::getMaxLength(compareObject, maxReadLength);
  readVec::getMaxLength(singleObject, maxReadLength);
  std::ofstream tempFile;
  if(extra){
  	openTextFile(tempFile, "tempFile.fastq", ".fastq",
  	  	setUp.pars_.ioOptions_.out_);
  }

	KmerMaps emptyMaps;
	auto gapPars = setUp.pars_.gapInfo_;
	aligner alignerObj(maxReadLength, gapPars, setUp.pars_.scoring_,
			setUp.pars_.colOpts_.alignOpts_.countEndGaps_);
	if (setUp.pars_.ioOptions_.firstName_ != "") {
		std::cout
				<< "readName\tscore\tposStart\tposStop\tcoverage\tidentity\tpercentGaps\teventBasedIdentity"
				<< std::endl;
		for (const auto &read : inReads) {
			auto positons = alignerObj.findReversePrimer(read, compareObject);
			alignerObj.rearrangeObjs(read, compareObject, true);
			// #! need to reorder seqs
			alignerObj.profilePrimerAlignment(read, compareObject);
			if (debug && extra) {
				alignerObj.alignObjectA_.seqBase_.outPutFastq(tempFile);
				alignerObj.alignObjectB_.seqBase_.outPutFastq(tempFile);
			}
      std::cout << read.seqBase_.name_ << "\t" << alignerObj.parts_.score_ << "\t"
      		<< positons.first << "\t" << positons.second
      		<< "\t" << alignerObj.comp_.distances_.query_.coverage_
      		<< "\t" << alignerObj.comp_.distances_.percentMatch_
      		<< "\t" << alignerObj.comp_.distances_.percentGaps_
      		<< "\t" << alignerObj.comp_.distances_.eventBasedIdentity_  << std::endl;
      if (alignerObj.comp_.distances_.query_.coverage_ >= queryCutOff &&
      		alignerObj.comp_.distances_.percentMatch_ >= percentIdentityCutoff &&
      		alignerObj.comp_.distances_.percentGaps_ <= percentageGapsCutOff &&
          extra) {
        read.seqBase_.outPutFastq(tempFile);
      }
    }
  } else {
    auto positons = alignerObj.findReversePrimer(singleObject, compareObject);
    alignerObj.rearrangeObjs(singleObject, compareObject, true);
    alignerObj.profilePrimerAlignment(singleObject, compareObject);
    std::cout << "readName\tscore\tposStart\tposStop\tcoverage\tidentity\tpercentGaps\teventBasedIdentity" << std::endl;
    std::cout << singleObject.seqBase_.name_ << "\t" << alignerObj.parts_.score_ << "\t"
    		<< positons.first << "\t" << positons.second
    		<< "\t" << alignerObj.comp_.distances_.query_.coverage_
    		<< "\t" << alignerObj.comp_.distances_.percentMatch_
    		<< "\t" << alignerObj.comp_.distances_.percentGaps_
    		<< "\t" << alignerObj.comp_.distances_.eventBasedIdentity_  << std::endl;
    //std::cout << positons.first << ":" << positons.second << std::endl;
    /*if(debug){

    	for(const auto & row : alignerObj.ScoreMatrix_){
    		for(const auto & col : row){
    			char temp = ' ';
    			double score = alignerObj.smithMaximum(col.upInherit, col.leftInherit, col.diagInherit, temp);
    			if(score > 0){
    				std::cout << changeBackground(46) << score << temp << "\t" << endAllAttributes();
    			}else{
    				std::cout << changeBackground(196) << score << temp << "\t" << endAllAttributes();
    			}
    		}
    		std::cout << std::endl;
    	}
    }*/

  }
  return 0;
}


int seqUtilsInfoRunner::profileReadsToReference(const njh::progutils::CmdArgs & inputCommands) {
  uint64_t maxReadLength = 0;

  seqUtilsInfoSetUp setUp(inputCommands);
  setUp.pars_.ioOptions_.lowerCaseBases_ = "upper";
  setUp.pars_.ioOptions_.out_.outFilename_ = "profile.tab.txt";
  setUp.processDefaultReader(true);
  setUp.processRefFilename(true);
  setUp.processKmerProfilingOptions();
  setUp.processGap();
  setUp.processQualThres();
  setUp.processScoringPars();
  setUp.finishSetUp(std::cout);

  SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<readObject>();

  readVec::getMaxLength(inReads, maxReadLength);
  std::vector<readObject> refSeqs = SeqInput::getReferenceSeq(
			setUp.pars_.refIoOptions_, maxReadLength);
	readVec::lowerCaseBasesToUpperCase(refSeqs);
	// create alinger class object
	KmerMaps kMaps = indexKmers(inReads, setUp.pars_.colOpts_.kmerOpts_.kLength_,
			setUp.pars_.colOpts_.kmerOpts_.runCutOff_,
			setUp.pars_.colOpts_.kmerOpts_.kmersByPosition_,
			setUp.pars_.expandKmerPos_, setUp.pars_.expandKmerSize_);
	auto gapPars = setUp.pars_.gapInfo_;
	// aligner alignerObj=aligner(maxReadLength, gapPars, scoringMatrixMap);
	aligner alignerObj(maxReadLength, gapPars, setUp.pars_.scoring_, kMaps,
			setUp.pars_.qScorePars_, setUp.pars_.colOpts_.alignOpts_.countEndGaps_,
			setUp.pars_.colOpts_.iTOpts_.weighHomopolyer_);
	readVec::allSetFractionByTotalCount(inReads);
	std::ofstream kMerInfo;
	openTextFile(kMerInfo, "kmerInfo.tab.txt", ".txt", true, true);
	//kmerMaps::outputKmerInfo(kMaps, kMerInfo);
	setUp.logRunTime(std::cout);
	return 0;
}




int seqUtilsInfoRunner::getHpProfile(const njh::progutils::CmdArgs & inputCommands) {
  seqUtilsInfoSetUp setUp(inputCommands);
  setUp.pars_.ioOptions_.out_.outFilename_ = "hpProfile.tab.txt";
  setUp.processDefaultReader(true);
  setUp.finishSetUp(std::cout);
  SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<readObject>();

  readVec::allSetFractionByTotalCount(inReads);
  readVec::lowerCaseBasesToUpperCase(inReads);
  std::map<std::string, std::map<uint32_t, hpRun>> runs =
      profiler::searchForHpRunsMultiple(inReads);
  std::ofstream profileFile;
  openTextFile(profileFile, setUp.pars_.ioOptions_.out_.outFilename_, ".tab.txt", false,
               false);

  profileFile << "ReadName\tclusterSize\tclusterFraction\tposition\thpPosition"
                 "\tbase\thpSize\tquals\tmeanQual" << std::endl;
  for (const auto &read : inReads) {
    for (auto &runIter : runs[read.seqBase_.name_]) {
      profileFile << read.seqBase_.name_ << "\t" << read.seqBase_.cnt_ << "\t"
                  << read.seqBase_.cnt_ << "\t";
      profileFile << runIter.second.getStringInfo() << std::endl;
    }
  }
  setUp.logRunTime(std::cout);
  return 0;
}



int seqUtilsInfoRunner::countOtus(const njh::progutils::CmdArgs & inputCommands) {
	std::cout << "Under construction" << std::endl;
	return 0;
  seqUtilsInfoSetUp setUp(inputCommands);
  /*
  int sizeCutOff = 1;
  bool condensed = false;
  bool grayScale = false;
  bool createDotFile = false;
  bool writeOutOtus = false;
  setUp.setOption(createDotFile, "-dot,-createDot", "CreateDotFile");
  setUp.setOption(grayScale, "-gray,-grayScale", "grayScale");
  setUp.setOption(condensed, "-condensed", "condensed");
  setUp.setOption(sizeCutOff, "-sizeCutOff,-cutOff", "SizeCutOff");
  setUp.pars_.ioOptions_.lowerCaseBases_ = "upper";
  setUp.processDefaultReader(true);
  setUp.setOption(writeOutOtus, "-write,-writeOut", "writeOutOtus");
  setUp.processQualThres();
  setUp.processGap();
  setUp.processKmerProfilingOptions();
  setUp.processScoringPars();
  // setUp.processSeq(true);
  double percentCutOff = 0.97;
  double gapCutOff = 0.01;
  double queryCutOff = 0.75;
  setUp.setOption(percentCutOff, "-idCutoff,-percentCutoff", "percentCutOff");
  setUp.setOption(gapCutOff, "-gapCutOff", "gapCutOff");
  setUp.setOption(queryCutOff, "-queryCutOff", "queryCutOff");
  setUp.processAlnInfoInput();
  setUp.finishSetUp(std::cout);

  SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<readObject>();

  uint64_t maxReadLength = 0;
  readVec::getMaxLength(inReads, maxReadLength);

  // create aligner class object
  kmerMaps kMaps;
  // kmerMaps kMaps = kmerCalculator::indexKmerMpas(inReads_,
  // setUp.pars_.kLength_,setUp.pars_.runCutoff_);
  auto gapPars = setUp.pars_.gapInfo_;
  // aligner alignerObj=aligner(maxReadLength, gapPars, scoringMatrixMap);
  aligner alignerObj(maxReadLength, gapPars,setUp.pars_.scoring_,
  		kMaps, setUp.pars_.primaryQual_,
                     setUp.pars_.secondaryQual_, setUp.pars_.qualThresWindow_,
                     setUp.pars_.countEndGaps_);
  alignerObj.processAlnInfoInput(setUp.pars_.alnInfoDirName_);
  std::vector<readObject> reads;
  if (condensed) {
    reads = createCondensedObjects(inReads);
  } else {
    reads = inReads;
  }

  std::pair<std::vector<readObject>, std::vector<readObject>>
      collapseReadsSplit;
  readVec::allSetFractionByTotalCount(inReads);
  if (setUp.pars_.ioOptions_.processed_) {
    collapseReadsSplit =
        readVecSplitter::splitVectorOnReadCount(inReads, sizeCutOff);
  } else {
    auto collapseReads = clusterCollapser::collapseIdenticalReads(
        reads, "median", setUp.pars_.ioOptions_.lowerCaseBases_);
    std::vector<readObject> convertedReads;
    for (const auto &cReads : collapseReads) {
      readObject tempObject(seqInfo(cReads.seqBase_.name_, cReads.seqBase_.seq_,
                                    cReads.seqBase_.qual_));
      tempObject.seqBase_.cnt_ = cReads.seqBase_.cnt_;
      tempObject.seqBase_.frac_ = cReads.seqBase_.frac_;
      tempObject.updateName();
      convertedReads.push_back(tempObject);
    }
    collapseReadsSplit =
        readVecSplitter::splitVectorOnReadCount(convertedReads, sizeCutOff);
  }

  // readVecSorter::sortReadVector(collapseReadsSplit.first, "size", true);
  // auto otus = roughCreateOTU<readObject, cluster>(
  //  collapseReadsSplit.first, alignerObj, percentCutOff, gapCutOff,
  // setUp.pars_.local_, setUp.pars_.weightHomopolymers_);
  std::vector<cluster> otus = roughCreateOTU<readObject, cluster>(
      collapseReadsSplit.first, alignerObj, percentCutOff, gapCutOff,queryCutOff,
      setUp.pars_.local_, setUp.pars_.weightHomopolymers_);
  if (createDotFile) {
    uint32_t count = 0;
    for (const auto &otu : otus) {
      ++count;
      std::cout << njh::bashCT::bold << "count: " + estd::to_string(count) << njh::bashCT::reset << std::endl;

      bestDistGraph testDist(otu.reads_, alignerObj, setUp.pars_.local_,
                             estd::to_string(count));
      std::ofstream outBestDistFile;
      openTextFile(outBestDistFile, "outBestDistFile", ".txt", true, false);
      outBestDistFile
          << "child\tparent\tchildFrac\tparentFrac\tdistance\tindelDi"
             "st\tmisDist\tidentity\tindels\tmismatches" << std::endl;
      testDist.printOutBestMatchInfos(outBestDistFile);

      std::unordered_map<std::string, njh::color> colorsForName;
      uint32_t count = 0;
      std::vector<njh::color> colors;
      if (grayScale) {
        std::vector<njh::color> tempColors =
            njh::evenHuesAll(0.75, 0.45, len(testDist.nodes_));
        for (const auto &c : tempColors) {
          colors.emplace_back(c.getGreyScaleColor());
        }
      } else {
        colors = njh::evenHuesAll(0.75, 0.45, len(testDist.nodes_));
      }
      for (const auto &n : testDist.nodes_) {
        colorsForName[n.read_.getReadId()] = colors[count];
        ++count;
      }
      testDist.createDotBestConnectedFile(setUp.pars_.directoryName_, colorsForName,
                                          false);
    }
  }
*/
  /*
  std::ofstream otusInfoFile;
  openTextFile(otusInfoFile, setUp.pars_.directoryName_ + "otusInfo.txt", ".txt",
  true, false);
  otusInfoFile << "otuNum\tread\tcnt\tbestRead\tbReadcCnt\tbestEditDistance" <<
  std::endl;
  uint32_t otuNum = 0;
  for(auto & otu : otus){
    readVec::allSetFractionByTotalCount(otu);
  }
  for(const auto & otu : otus){
    std::ofstream otusGraphFile;
    openTextFile(otusGraphFile, setUp.pars_.directoryName_ + std::to_string(otuNum) +
  "_gv.dot", ".dot", true, false);
    std::ofstream otusFirstLevelFile;
    openTextFile(otusFirstLevelFile, setUp.pars_.directoryName_ +
  std::to_string(otuNum) + "_firstLevel", ".txt", true, false);
    otuGraph currentGraph = otuGraph(otu, alignerObj, setUp.pars_.local_,
  std::to_string(otuNum));
    currentGraph.printInfo(otusInfoFile);
    currentGraph.printInfoGraphViz(otusGraphFile);
    currentGraph.printParentLevel(otusFirstLevelFile);
    ++otuNum;
  }*/

  // readVec::allUpdateName(otus);
  // readVecSorter::sort(otus);
  // reader.write(otus, "otus", setUp.pars_.ioOptions_.out_.outFormat,
  //           setUp.pars_.ioOptions_.overWriteFile,
  //         setUp.pars_.ioOptions_.exitOnFailureToWrite);
  /*
  if (writeOutOtus) {
    uint32_t otuCount = 0;
    std::string otuDirectory = njh::files::makeDir(
        "./", njh::files::getFileName(setUp.pars_.ioOptions_.firstName_) + "_otu_TODAY");
    for (const auto &otu : otus) {
      reader.write(otu.reads_, otuDirectory + estd::to_string(otuCount),
                   setUp.pars_.ioOptions_);
      ++otuCount;
    }
  }
  std::cout << "otus: " << otus.size() << std::endl;
  setUp.logRunTime(std::cout);
  */
  return 0;
}



int seqUtilsInfoRunner::quickLenInfo(const njh::progutils::CmdArgs & inputCommands) {
	VecStr additionalColumnsData;
	uint32_t testNumber = std::numeric_limits<uint32_t>::max();
  seqUtilsInfoSetUp setUp(inputCommands);
  setUp.processDebug();
  setUp.processVerbose();
  OutOptions outOpts;
  outOpts.outExtention_ = ".tab.txt";
  setUp.processWritingOptions(outOpts);
  setUp.processReadInNames();
  setUp.setOption(testNumber, "--testNumber", "testNumber");
  setUp.setOption(additionalColumnsData, "--additionalColumnsData", "Extra columns to add to the output table, comma separated values in the format of [COL_NAME]:[COL_VAL]");
  setUp.finishSetUp(std::cout);
	std::vector<uint64_t> readLengths;
	readObject read;
	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	OutputStream out(outOpts);
	uint32_t readCount = 0;
	while (reader.readNextRead(read)) {
		readLengths.emplace_back(len(read));
		++readCount;
		if(readCount > testNumber){
			break;
		}
	}
	auto stats = getStatsOnVec(readLengths);
	table statsTable;
	statsTable.columnNames_ = getVectorOfMapKeys(stats);
	statsTable.hasHeader_ = true;
	auto nums = getVectorOfMapValues(stats);
	statsTable.content_.emplace_back(numVecToVecStr(nums));
	statsTable.columnNames_.emplace_back("n");
	statsTable.content_.front().emplace_back(njh::pasteAsStr(readLengths.size()));
	for(const auto & d : additionalColumnsData){
		auto toks = tokenizeString(d, ":");
		if(2 != toks.size()){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "could not separate additional column on :, " << d << "\n";
			throw std::runtime_error{ss.str()};
		}
		statsTable.columnNames_.emplace_back(toks[0]);
		statsTable.content_.front().emplace_back(toks[1]);
	}
	statsTable.outPutContents(out, "\t");

//  if(plot){
//  	njhRInside::OwnRInside rSes;
//  	rSes.loadLib("ggplot2");
//  	auto & R = rSes.get();
//  	R["readLengths"] = readLengths;
//  	rSes.openPdfDevice(setUp.pars_.ioOptions_.out_.outFilename_.string(), 11.5, 8);
//  	if(std::numeric_limits<uint32_t>::max() != binSize  ){
//  		std::string plotCmd = "length_p = ggplot(data.frame(readLengths = readLengths), aes(x = readLengths) ) + stat_bin(binwidth=" + estd::to_string(binSize)+ " )";
//  		if(std::numeric_limits<uint32_t>::max() != start  && std::numeric_limits<uint32_t>::max() != stop){
//  			plotCmd.append(" + xlim(" + estd::to_string(start) + "," + estd::to_string(stop) + ")");
//  		}
//  		plotCmd.append("; print(length_p)");
//  		if(setUp.pars_.debug_){
//  			std::cout << plotCmd << std::endl;
//  		}
//  		R.parseEvalQ(plotCmd);
//  	}else{
//  		R.parseEvalQ("print(ggplot(data.frame(readLengths = readLengths), aes(x = readLengths) ) + geom_bar() )");
//  	}
//  	rSes.closeCurrentDevice();
//  }
  return 0;
}

int seqUtilsInfoRunner::quickMismatchDist(const njh::progutils::CmdArgs & inputCommands) {
  seqUtilsInfoSetUp setUp(inputCommands);
  setUp.pars_.ioOptions_.lowerCaseBases_ = "upper";
  setUp.processDefaultReader(true);

  setUp.processQualThres();
  setUp.processGap();
  setUp.processKmerProfilingOptions();
  setUp.processScoringPars();

  setUp.finishSetUp(std::cout);

  SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<readObject>();

  uint64_t maxReadLength = 0;
  readVec::getMaxLength(inReads, maxReadLength);

  // create aligner class object
  KmerMaps kMaps;
  // kmerMaps kMaps = kmerCalculator::indexKmerMpas(inReads_,
  // setUp.pars_.kLength_,setUp.pars_.runCutoff_);
  auto gapPars = setUp.pars_.gapInfo_;
  //auto gapPars = setUp.pars_.gapInfo_;
  // aligner alignerObj=aligner(maxReadLength, gapPars, scoringMatrixMap);
	aligner alignerObj(maxReadLength, gapPars, setUp.pars_.scoring_, kMaps,
			setUp.pars_.qScorePars_, setUp.pars_.colOpts_.alignOpts_.countEndGaps_,
			setUp.pars_.colOpts_.iTOpts_.weighHomopolyer_);
  std::ofstream distFile;
  openTextFile(distFile,
  		bfs::basename(setUp.pars_.ioOptions_.firstName_) + "_distFile.txt",
               ".txt", true, false);
  for (const auto pos : iter::range(inReads.size())) {
    distFile << " \t" << inReads[pos].seqBase_.name_;
  }
  distFile << std::endl;
  for (const auto pos : iter::range(inReads.size())) {
    std::cout << pos << ":" << inReads.size() << std::endl;
    distFile << inReads[pos].seqBase_.name_;
    for (const auto secondPos : iter::range(inReads.size())) {
      alignerObj.alignCache(inReads[pos], inReads[secondPos],
                          setUp.pars_.local_);
      alignerObj.profileAlignment(inReads[pos], inReads[secondPos],
                                  false, true, false);
      distFile << "\t" << alignerObj.comp_.hqMismatches_;
    }
    distFile << std::endl;
  }

  return 0;
}


int seqUtilsInfoRunner::fracInfo(const njh::progutils::CmdArgs & inputCommands) {
	//std::cout << "1" << std::endl;
  seqUtilsInfoSetUp setUp(inputCommands);
  setUp.processDefaultReader(true);
  setUp.finishSetUp(std::cout);
  SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<readObject>();

  std::ofstream outFile;
  if(setUp.pars_.ioOptions_.out_.outFilename_ == "out"){
  	setUp.pars_.ioOptions_.out_.outFilename_ = njh::files::removeExtension(bfs::basename(setUp.pars_.ioOptions_.firstName_)) + "Info";
  }
  openTextFile(outFile, setUp.pars_.ioOptions_.out_.outFilename_, ".tab.txt", true, false);
  profiler::getFractionInfo(inReads, outFile);
  return 0;
}




Json::Value toJsonMismatchGraph(readDistGraph<uint32_t> & graph,uint32_t groupCutOff,
			uint32_t mismatchesAllowed,
			njh::color backgroundColor, double hueStart, double hueStop,
	    double lumStart, double lumStop,
	    double satStart, double satStop, bool determineBest){

	Json::Value graphJson;
	graphJson["backgroundColor"] = "#" +  backgroundColor.hexStr_ ;
	auto & nodes = graphJson["nodes"];
	auto & links = graphJson["links"];
	uint32_t nCount = 0;
	std::unordered_map<std::string, uint64_t> nameToNewPos;
	graph.turnOffEdgesAbove(mismatchesAllowed);
	if(determineBest){
		graph.allDetermineLowestBest(true);
	}
	//graph.determineGroups();
	uint32_t groupCount = 0;
	for(auto & n : graph.nodes_){
		n->group_ = groupCount;
		++groupCount;
	}

	njh::randomGenerator gen;
	uint64_t pos = 0;
	double minReadCnt = std::numeric_limits<double>::max();
	double maxReadCnt = std::numeric_limits<double>::lowest();
	std::unordered_map<uint32_t, njh::color> groupColors;
  std::unordered_map<uint32_t, uint32_t> groupCounts;
  uint32_t numOfCutOffGroups = 0;
  for(const auto & n : graph.nodes_){
  	++groupCounts[n->group_];
  }
  std::vector<uint32_t> groups;
	for(const auto & g : groupCounts){
		if(g.second >= groupCutOff){
			++numOfCutOffGroups;
			groups.emplace_back(g.first);
		}
	}
	//printOutMapContents(groupCounts,"\t", std::cout);
	//printVector(groups);
	auto gColors = njh::getColsBetweenInc(hueStart, hueStop,
  		lumStart, lumStop,
  		satStart, satStop,
  		groups.size());

	for(const auto pos : iter::range(groups.size())){
		groupColors[groups[pos]] = gColors[pos];
	}

	for(const auto & n :graph.nodes_){
		if(n->value_->cnt_ < minReadCnt){
			minReadCnt = n->value_->cnt_;
		}
		if(n->value_->cnt_ > maxReadCnt){
			maxReadCnt = n->value_->cnt_;
		}
	}
	//scale<double> cntScale({minReadCnt, maxReadCnt},{50.0, 1000.0});
	scale<double> cntScale({0, maxReadCnt},{50.0, 1000.0});
	for(const auto & n : graph.nodes_){
  	if(groupCounts[n->group_] >= groupCutOff){
			nameToNewPos[n->name_] = pos;
			++pos;
			//std::cout << n->name_ << " : " << n->group_  << " : " << n->value_ << std::endl;
			nodes[nCount]["name"] = njh::json::toJson(n->name_);
			nodes[nCount]["group"] = njh::json::toJson(n->group_);
			nodes[nCount]["color"] = "#" + groupColors[n->group_].hexStr_;
			nodes[nCount]["size"] = cntScale.get(n->value_->cnt_);
			++nCount;
  	}
	}
	uint32_t lCount=0;
	for(const auto & e : graph.edges_){
		if(e->on_ && groupCounts[e->nodeToNode_.begin()->second.lock()->group_] >= groupCutOff){
			if(e->dist_ == 0){
				links[lCount]["source"] = njh::json::toJson(nameToNewPos[graph.nodes_[graph.nameToNodePos_[e->nodeToNode_.begin()->first]]->name_]);
				links[lCount]["target"] = njh::json::toJson(nameToNewPos[e->nodeToNode_.begin()->second.lock()->name_]);
				links[lCount]["value"] = njh::json::toJson(1);
				links[lCount]["on"] = njh::json::toJson(e->on_);
				auto lColor = getColsBetweenExcludeClosest(groupColors[graph.nodes_[graph.nameToNodePos_[e->nodeToNode_.begin()->first]]->group_],
						groupColors[e->nodeToNode_.begin()->second.lock()->group_], 1);
				links[lCount]["color"] = "#" + lColor.front().hexStr_;
				++lCount;
			}else{
				std::string lastName = graph.nodes_[graph.nameToNodePos_[e->nodeToNode_.begin()->first]]->name_;
				auto lColors = getColsBetweenExcludeClosest(groupColors[graph.nodes_[graph.nameToNodePos_[e->nodeToNode_.begin()->first]]->group_],
						groupColors[e->nodeToNode_.begin()->second.lock()->group_], e->dist_ + 1);
				for(const auto mis : iter::range(e->dist_)){
					std::string newName = graph.nodes_[graph.nameToNodePos_[e->nodeToNode_.begin()->first]]->name_
							+ estd::to_string(mis) + e->nodeToNode_.begin()->second.lock()->name_;
					nameToNewPos[newName] = pos;
					++pos;
					nodes[nCount]["name"] = njh::json::toJson(newName);
					nodes[nCount]["group"] = njh::json::toJson(e->nodeToNode_.begin()->second.lock()->group_);
					nodes[nCount]["color"] = "red";
					nodes[nCount]["size"] = 10;
					++nCount;
					links[lCount]["source"] = njh::json::toJson(nameToNewPos[lastName]);
					links[lCount]["target"] = njh::json::toJson(nameToNewPos[newName]);
					links[lCount]["value"] = njh::json::toJson(1);
					links[lCount]["on"] = njh::json::toJson(true);
					links[lCount]["color"] = "#" + lColors[mis].hexStr_;
					++lCount;
					lastName = newName;
				}

				links[lCount]["source"] = njh::json::toJson(nameToNewPos[lastName]);
				links[lCount]["target"] = njh::json::toJson(nameToNewPos[e->nodeToNode_.begin()->second.lock()->name_]);
				links[lCount]["value"] = njh::json::toJson(1);
				links[lCount]["on"] = njh::json::toJson(true);
				links[lCount]["color"] = "#" + lColors[e->dist_].hexStr_;
				++lCount;
			}
		}
	}
	return graphJson;
}



int seqUtilsInfoRunner::genPsuedoAllMinTree(const njh::progutils::CmdArgs & inputCommands){
  seqSetUp setUp(inputCommands);
  uint32_t numThreads = 2;
  setUp.setOption(numThreads, "--numThreads,-t", "Number of threads to use", false);
	comparison allowableErrors;
	bool settingEventsLimits = setUp.setOption(
			allowableErrors.distances_.overLappingEvents_, "--numberOfEvents",
			"Make a cut at this number of events and above");
	setUp.processDefaultReader(true);
  setUp.processAlignerDefualts();
  setUp.processDirectoryOutputName(true);
  setUp.processVerbose();
  setUp.processDebug();
  setUp.finishSetUp(std::cout);
  setUp.startARunLog(setUp.pars_.directoryName_);
  SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<readObject>();

  std::vector<std::shared_ptr<readObject>> reads;
  for(const auto & read : inReads){
  	reads.emplace_back(std::make_shared<readObject>(read));
  }
	uint64_t maxSize = 0;
	readVec::getMaxLength(reads, maxSize);
	aligner alignerObj = aligner(maxSize, setUp.pars_.gapInfo_,
			setUp.pars_.scoring_, setUp.pars_.colOpts_.alignOpts_.countEndGaps_);
  alignerObj.processAlnInfoInput(setUp.pars_.alnInfoDirName_);
  std::unordered_map<std::string, std::unique_ptr<aligner>> aligners;
  std::mutex alignerLock;
	auto graph = genReadComparisonGraph(reads, alignerObj, aligners, alignerLock,
			numThreads);
	std::vector<std::string> popNames;
	for (const auto & n : graph.nodes_) {
		if (n->on_) {
			popNames.emplace_back(n->name_);
		}
	}
	auto nameColors = getColorsForNames(popNames);
	if(settingEventsLimits){
		graph.turnOffEdgesWithComp(allowableErrors,
				[](const comparison & comp1, const comparison & cutOff){
			//std::cout << comp1.toJson() << std::endl;
			return comp1.distances_.getNumOfEvents(true) >= cutOff.distances_.overLappingEvents_;
		});
	}else{
		comparison maxEvents = graph.setMinimumEventConnections();
		if(setUp.pars_.debug_){
			std::cout << maxEvents.toJson() << std::endl;
		}
	}
	auto treeData = graph.toD3Json(njh::color("#000000"), nameColors);

	std::ofstream outJson(setUp.pars_.directoryName_ + "tree.json");
	std::ofstream outHtml(setUp.pars_.directoryName_ + "tree.html");
	outJson << treeData;
	auto outHtmlStr = genHtmlStrForPsuedoMintree("tree.json",
			"http://njh8.umassmed.edu/~hathawan/js/psuedoMinTreeWithIndels.js");
	outHtml << outHtmlStr;
  if(setUp.pars_.writingOutAlnInfo_){
  	for(const auto & alnObj : aligners){
  		if(setUp.pars_.debug_){
  			std::cout << alnObj.first << ": " << alnObj.second->numberOfAlingmentsDone_ << std::endl;
  		}
  		if(alnObj.second->numberOfAlingmentsDone_ > 0){
  			alignerObj.alnHolder_.mergeOtherHolder(alnObj.second->alnHolder_);
  		}
  	}
  	alignerObj.alnHolder_.write(setUp.pars_.outAlnInfoDirName_, setUp.pars_.verbose_);
  }
  return 0;
}

int seqUtilsInfoRunner::genPsuedoMismatchMinTree(const njh::progutils::CmdArgs & inputCommands){
	seqUtilsInfoSetUp setUp(inputCommands);
	uint32_t numThreads = 1;
	uint32_t strNum = 20;
	uint32_t strLen = 40;
	std::string backgroundColor = "#000000";
	double hueStart = 0;
	double hueStop = 360;
	double satStart = 1.0;
	double satStop = 1.0;
	double lumStart = 0.5;
	double lumStop = 0.5;
	uint32_t mismatches = 0;
	bool setMismatches = false;
	std::string colorFile = "";
	uint32_t groupCutOff = 2;
	bool mNoGroupColor = false;
	bool determineBest = false;
	setUp.setOption(determineBest, "--determineBest", "Determine Best Mismatches");
	setUp.setOption(mNoGroupColor, "--mNoGroupColor", "Don't do group coloring for graph");
	setMismatches = setUp.setOption(mismatches, "--mismatches", "Show Only up to this number of mismatches");
	setUp.setOption(groupCutOff, "--groupCutOff", "Group Size Cut Off To be included");
	setUp.setOption(numThreads, "-t,--numThreads", "Number of Threads to Use");
	setUp.setOption(strLen, "-l,--strLen", "Length of Strings to Generate");
	setUp.setOption(strNum, "-n,--strNum", "Number of Strings to Generate");
	setUp.setOption(backgroundColor, "-b,--backgroundColor", "Hex String for the color for the Background");
	setUp.setOption(hueStart, "--hueStart", "Hue Start for the color of Reads");
	setUp.setOption(hueStop, "--hueStop", "Hue Stop for the color of Reads");
	setUp.setOption(satStart, "--satStart", "Sat Start for the color of Reads");
	setUp.setOption(satStop, "--satStop", "Sat Stop for the color of Reads");
	setUp.setOption(lumStart, "--lumStart", "Lum Start for the color of Reads");
	setUp.setOption(lumStop, "--lumStop", "Lum Stop for the color of Reads");
	setUp.setOption(colorFile, "--colorFile", "Name of a color file, with first colum read name and second column is hex color string");
	setUp.processDefaultReader(false);
	setUp.processAlignerDefualts();
	setUp.processVerbose();
	setUp.processDirectoryOutputName(true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<readObject>();
	if(setUp.pars_.ioOptions_.firstName_ == ""){
		njh::randomGenerator gen;
		VecStr rSeqs = simulation::evenRandStrs(strLen, std::vector<char>{'A', 'C', 'G', 'T'}, gen, strNum);
		inReads = vecStrToReadObjs<readObject>(rSeqs, "Seq");
		for(auto & read : inReads){
			read.seqBase_.cnt_ = gen.unifRand(1.0,100.0);
			read.updateName();
		}
	}

  njh::stopWatch watch;
  uint64_t maxLength = 0;
  readVec::getMaxLength(inReads, maxLength);
  aligner alignerObj(maxLength, setUp.pars_.gapInfo_, setUp.pars_.scoring_);
  alignerObj.processAlnInfoInput(setUp.pars_.alnInfoDirName_);

  watch.startNewLap("regular distance");

  std::function<uint32_t(const readObject & ,
  		const readObject &, aligner)> misFun = getMismatches<readObject>;
	auto misDistances = getDistanceCopy(inReads, numThreads, misFun,
			alignerObj);
	watch.startNewLap("create graph");
  readDistGraph<uint32_t> graphMis(misDistances, inReads);
	std::vector<std::string> names;
  for(const auto & n : graphMis.nodes_){
  	names.emplace_back(n->name_);
  }
  if(hueStop == 360 && hueStart == 0){
  	hueStop = 360 - (360.0/names.size());
  }
	auto nameColors = njh::getColorsForNames(names, hueStart, hueStop,
			lumStart, lumStop, satStop, satStart);
	table colorToName;

	if(colorFile != ""){
		VecStr adding;
		colorToName = table(colorFile,"\t", true);
		if(colorToName.columnNames_.size() != 2){
			std::cerr << "Wrong number of columns, should only have two, "
					<< colorFile << " has " << colorToName.columnNames_.size()
					<< std::endl;
			exit(1);
		}
		for(const auto & row : colorToName.content_){
			if(stringToLowerReturn(row[0]) == "backgroundcolor"){
				backgroundColor = row[1];
				continue;
			}
			adding.emplace_back(row[0]);
			if(nameColors.find(row[0]) == nameColors.end()){
				std::cerr << "Warning, adding " << row[0] << " but it wasn't found in the reads" << std::endl;
			}
			nameColors[row[0]] = njh::color(row[1]);
		}
		for(const auto & n : names){
			if(!njh::in(n, adding)){
				std::cerr << "Warning, add a file to specify colors for reads but didn't have a color set for " << n << std::endl;
			}
		}
	}


  watch.startNewLap("determine connections");
  Json::Value graphJson;
  if(setMismatches){
  	if(mNoGroupColor){
  		graphJson = toJsonMismatchGraph(graphMis,groupCutOff,mismatches, njh::color(backgroundColor), hueStart, hueStop,
  		    			lumStart, lumStop, satStop, satStart, determineBest);
  	}else{
    	graphJson = graphMis.toJsonMismatchGraph(groupCutOff,mismatches, njh::color(backgroundColor), hueStart, hueStop,
    			lumStart, lumStop, satStop, satStart);
  	}

  }else{
  	graphJson = graphMis.toJsonMismatchGraphAll(njh::color(backgroundColor), nameColors);
  }
  if(setUp.pars_.verbose_){
  	watch.logLapTimes(std::cout, true, 6, true);
  }
	std::cout << graphMis.numberOfGroups_ << std::endl;
	std::ofstream outFile(setUp.pars_.directoryName_ + "psudoMinTree.json");
	outFile << graphJson << std::endl;
	std::string htmlOut = genHtmlStrForPsuedoMintree("psudoMinTree.json");
	std::ofstream outHtmlFile(setUp.pars_.directoryName_ + "psudoMinTree.html");
	outHtmlFile << htmlOut << std::endl;

	if(setUp.pars_.writingOutAlnInfo_){
		alignerObj.alnHolder_.write(setUp.pars_.outAlnInfoDirName_);
	}


	return 0;
}

} // namespace njhseq
