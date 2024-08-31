
//  seqUtilsModRunner.cpp
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
    
#include "seqUtilsModRunner.hpp"
#include <njhseq/objects/counters/DNABaseCounter.hpp>
#include <njhseq/IO/SeqIO/SeqIO.hpp>
#include <njhseq.h>


namespace njhseq {

seqUtilsModRunner::seqUtilsModRunner()
    : njh::progutils::ProgramRunner({
		addFunc("removeLowQualityBases", removeLowQualityBases, false),
		addFunc("translate", translate, false),
		addFunc("guessAProteinFromSeq", guessAProteinFromSeq, false),
		addFunc("revCompSeq", revCompSeq, false),
		addFunc("renameIDs", renameIDs, false),
		addFunc("sortReads", sortReads, false),
		addFunc("appendReads", appendReads, false),
		addFunc("prependReads", prependReads, false),
		addFunc("reOrientReads", reOrientReads, false),
		addFunc("collapseToUnique", collapseToUnique, false),
		addFunc("dereplicate", dereplicate, false),
		addFunc("prependNames", prependNames, false),
		addFunc("appendNames", appendNames, false),
		addFunc("changeLetterCase", changeLetterCase, false),
		addFunc("invertLetterCase", inverseLetterCase, false),
		addFunc("collapseToUniqueWithinMetaField", collapseToUniqueWithInMetaField, false),
		addFunc("expandOutCollapsedToUnique", expandOutCollapsedToUnique, false),
		addFunc("increaseQualityScores", increaseQualityScores, false),
		addFunc("sortReadsByEntropy", sortReadsByEntropy, false),
		addFunc("sortReadsByKmerEntropy", sortReadsByKmerEntropy, false),
		addFunc("sortReadsByKCompToTop", sortReadsByKCompToTop, false),
		addFunc("sortReadsByNameNaturalSort", sortReadsByNameNaturalSort, false),
		addFunc("renameSeqNameToUniqueNames", renameSeqNameToUniqueNames, false),
		addFunc("sortReadsPairedEnd", sortReadsPairedEnd, false),
		addFunc("compSeq", compSeq, false),
    	addFunc("changeLetterToOtherLetter", changeLetterToOtherLetter, false),


},//
                    "seqUtilsMod") {}

int seqUtilsModRunner::changeLetterToOtherLetter(const njh::progutils::CmdArgs & inputCommands){
	char base = ' ';
	char replacement = ' ';
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDefaultReader( true);
	setUp.setOption(base, "--base", "Base to replace", true);
	setUp.setOption(replacement, "--replacement", "replacement base", true);

	setUp.finishSetUp(std::cout);

	seqInfo seq;
	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	reader.openOut();
	std::regex baseRegex(std::string(1, base));
	auto replacementStr = std::string(1, replacement);
	while(reader.readNextRead(seq)){
		seq.seq_ = std::regex_replace(seq.seq_,baseRegex, replacementStr);
		reader.write(seq);
	}
	return 0;
}


int seqUtilsModRunner::increaseQualityScores(const njh::progutils::CmdArgs & inputCommands){
	uint32_t qualIncrease = 1;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.setOption(qualIncrease, "--qualIncrease", "Increase quality scores by this much");
	setUp.processDefaultReader({"--fastq", "--fastqgz"}, true);
	setUp.finishSetUp(std::cout);

	seqInfo seq;
	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	reader.openOut();
	while(reader.readNextRead(seq)){
		for(auto & q : seq.qual_){
			q += qualIncrease;
		}
		reader.write(seq);
	};
	return 0;
}


int seqUtilsModRunner::expandOutCollapsedToUnique(const njh::progutils::CmdArgs & inputCommands){
	bfs::path inputNames = "";
	std::string delim = ",";
	std::string readsColName = "reads";
	std::string nameColName = "name";
	seqSetUp setUp(inputCommands);
	setUp.pars_.ioOptions_.out_.outFilename_ = "";
	setUp.processDefaultReader({ "--fasta", "--fastagz", "--fastq", "--fastqgz"}, true);
	setUp.setOption(inputNames, "--inputNames", "Input names", true);
	setUp.setOption(delim, "--delim", "Delimiter", false);

	setUp.setOption(readsColName, "--readsColName", "Reads Col Name, containing the unique names", false);
	setUp.setOption(nameColName, "--nameColName", "Name Col Name, containing the name in the collapsed file", false);

	setUp.finishSetUp(std::cout);

	table namesTab(inputNames, "\t", true);
	namesTab.checkForColumnsThrow(VecStr{nameColName,readsColName}, __PRETTY_FUNCTION__)	;

	std::unordered_map<std::string, VecStr> readNames;
	auto nameColPos = namesTab.getColPos(nameColName);
	auto readsColPos = namesTab.getColPos(readsColName);
	for(const auto & row : namesTab.content_){
		addOtherVec(readNames[row[nameColPos]], tokenizeString(row[readsColPos], delim));
	}

	SeqIO reader(setUp.pars_.ioOptions_);

	reader.openIn();
	reader.openOut();

	seqInfo seq;
	while(reader.readNextRead(seq)){
		if(!njh::in(seq.name_, readNames)){
			reader.write(seq);
		}else{
			for(const auto & outName : readNames.at(seq.name_)){
				auto outSeq = seq;
				outSeq.name_ = outName;
				reader.write(outSeq);
			}
		}
	}
	return 0;
}

int seqUtilsModRunner::guessAProteinFromSeq(
		const njh::progutils::CmdArgs & inputCommands) {
	bool getLongest = false;
	bool removeTrailingStop = false;
	bool mark = false;
	bool overWriteMeta = false;
	bool forceMStart = false;
	seqSetUp setUp(inputCommands);
	setUp.processDefaultReader(true);
	if ("out" == setUp.pars_.ioOptions_.out_.outFilename_) {
		setUp.pars_.ioOptions_.out_.outFilename_ = njh::files::prependFileBasename(
				setUp.pars_.ioOptions_.firstName_, "translated_");
	}
	setUp.setOption(getLongest, "--getLongest", "Get the longest possible protein");
	setUp.setOption(mark, "--mark", "Mark the out seq with the frame used");
	setUp.setOption(overWriteMeta, "--overWriteMeta", "If marking whether or not to reset the meta already in the sequence name");

	setUp.setOption(removeTrailingStop, "--removeTrailingStop", "Remove Trailing Stop Codon");
	setUp.setOption(forceMStart, "--forceMStart", "force M Start");

	setUp.finishSetUp(std::cout);

	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	reader.openOut();

	seqInfo seq;
	njh::PatPosFinder pFinder(R"(\*+)");

	while (reader.readNextRead(seq)) {
		if (getLongest) {
			std::shared_ptr<readVecTrimmer::BreakUpRes> longest;
			uint32_t frameWithLongest = std::numeric_limits<uint32_t>::max();
			uint32_t lengthOfLongestFrame = std::numeric_limits<uint32_t>::max();
			bool isRevComp = false;
			for (auto comp : {false, true}){
				for (const auto start : iter::range(3)) {
					auto prot = seq.translateRet(comp, comp, start, forceMStart);
					auto possibleSubSeqs = readVecTrimmer::breakUpSeqOnPat(prot, pFinder);
					for (const auto & possible : possibleSubSeqs) {
						if (nullptr == longest) {
							longest = std::make_shared<readVecTrimmer::BreakUpRes>(possible);
							frameWithLongest = start;
							lengthOfLongestFrame = len(prot);
							isRevComp = comp;
						} else if (len(possible.seqBase_) > len(longest->seqBase_)) {
							longest = std::make_shared<readVecTrimmer::BreakUpRes>(possible);
							frameWithLongest = start;
							lengthOfLongestFrame = len(prot);
							isRevComp = comp;
						}
					}
				}
			}

			uint32_t seqStart = longest->start_ * 3 + frameWithLongest;
			uint32_t seqStop = longest->end_ * 3 + frameWithLongest;
			uint32_t aaStop = longest->end_;
			if(!removeTrailingStop && longest->end_ != lengthOfLongestFrame){

				//if the end doesn't equal the length of the protein that must mean that it broke on a stop codon
				longest->seqBase_.append("*");
				++aaStop;
				seqStop += 3;
//				std::cout << "longest->end_: " << longest->end_ << "  lengthOfLongestFrame: " << lengthOfLongestFrame<< std::endl;
//				auto protDebug = seq.translateRet(false, false, frameWithLongest);
//				std::cout << '\t' << protDebug.seq_.size() << " " << protDebug.seq_ << std::endl;
//				std::cout << "\t" << seq.seq_.substr(seqStart, seqStop - seqStart) << std::endl;
			}
			if(mark){
				MetaDataInName meta;
				bool nameHasMeta = false;
				if(MetaDataInName::nameHasMetaData(seq.name_)){
					meta = MetaDataInName(seq.name_);
					nameHasMeta = true;
				}
				meta.addMeta("frame", frameWithLongest, overWriteMeta);
				meta.addMeta("seqStart", seqStart, overWriteMeta);
				meta.addMeta("seqStop", seqStop, overWriteMeta);
				meta.addMeta("seqLen", len(seq), overWriteMeta);
				meta.addMeta("aaStart", longest->start_, overWriteMeta);
				meta.addMeta("aaStop", aaStop, overWriteMeta);
				meta.addMeta("aaLen", lengthOfLongestFrame, overWriteMeta);
				meta.addMeta("isRevComp", isRevComp);
				if(nameHasMeta){
					meta.resetMetaInName(longest->seqBase_.name_);
				}else{
					longest->seqBase_.name_.append(" " + meta.createMetaName());
				}
			}

			reader.write(longest);

		} else {
			std::unordered_map<bool, std::vector<seqInfo>>  translatedSeqs;
			translatedSeqs[false] =  std::vector<seqInfo>{
				seq.translateRet(false, false, 0, forceMStart),
				seq.translateRet(false, false, 1, forceMStart),
				seq.translateRet(false, false, 2, forceMStart)};
			translatedSeqs[true] =  std::vector<seqInfo>{
				seq.translateRet(true, true, 0, forceMStart),
				seq.translateRet(true, true, 1, forceMStart),
				seq.translateRet(true, true, 2, forceMStart)};


			uint32_t minPos = 0;
			uint32_t minStopCodon = std::numeric_limits<uint32_t>::max();
			bool isRevComp = false;
			for (auto comp : {false, true}){
				for(auto pos: iter::range(translatedSeqs.size())){
					auto stops =  countOccurences(translatedSeqs[comp][pos].seq_, "*");
					if(stops < minStopCodon){
						minStopCodon = stops;
						minPos = pos;
						isRevComp = comp;
					}
				}
			}

			seqInfo seqCopy = translatedSeqs[isRevComp][minPos];
			if (mark) {
				MetaDataInName meta;
				bool nameHasMeta = false;
				if (MetaDataInName::nameHasMetaData(seqCopy.name_)) {
					meta = MetaDataInName(seqCopy.name_);
					nameHasMeta = true;
				}
				meta.addMeta("frame", minPos);
				meta.addMeta("isRevComp", isRevComp);

				if (nameHasMeta) {
					meta.resetMetaInName(seqCopy.name_);
				} else {
					seqCopy.name_.append(" " + meta.createMetaName());
				}
			}

			if(removeTrailingStop ){
				readVecTrimmer::trimAtRstripBase(seqCopy, '*');
			}
			reader.write(seqCopy);
		}
	}
	return 0;
}

int seqUtilsModRunner::reOrientReads(
		const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	setUp.processDefaultReader(true);
	bool reOrientToBestWinner  = false;
	uint32_t numThreads = 1;
	if (setUp.pars_.ioOptions_.out_.outFilename_ == "out") {
		setUp.pars_.ioOptions_.out_.outFilename_ = njh::files::prependFileBasename(
				setUp.pars_.ioOptions_.firstName_, "reOriented_");
	}
	setUp.setOption(setUp.pars_.colOpts_.kmerOpts_.kLength_, "-k,--kLength",
			"kLength");
	setUp.setOption(reOrientToBestWinner, "---reOrientToBestWinner", "re-orient To Best Winner");
	setUp.setOption(numThreads, "---numThreads", "numThreads");

	setUp.processRefFilename(false);
	setUp.processSeq(false);
	setUp.finishSetUp(std::cout);

	if ("" != setUp.pars_.refIoOptions_.firstName_) {
		SeqInput refReader(setUp.pars_.refIoOptions_);
		auto refSeqs = refReader.readAllReads<seqInfo>();
		std::vector<std::unique_ptr<seqWithKmerInfo>> refKmerReads;
		for (const auto & seq : refSeqs) {
			refKmerReads.emplace_back(std::make_unique<seqWithKmerInfo>(seq));
		}
		allSetKmers(refKmerReads, setUp.pars_.colOpts_.kmerOpts_.kLength_, true);
		SeqIO reader(setUp.pars_.ioOptions_);
		reader.openIn();
		reader.openOut();

		std::mutex outMut;
		std::function<void()> reOrient = [&reader,&outMut,&setUp,&reOrientToBestWinner,&refKmerReads](){
			seqInfo seq;
			while (reader.in_.readNextReadLock(seq)) {
				auto seqKmer = std::make_unique<seqWithKmerInfo>(seq);
				seqKmer->setKmers(setUp.pars_.colOpts_.kmerOpts_.kLength_, true);
				uint32_t forwardWinners = 0;
				uint32_t revWinners = 0;
				if(reOrientToBestWinner){
					uint32_t best = 0;
					for (const auto & refSeq : refKmerReads) {
						auto forDist = refSeq->compareKmers(*seqKmer);
						auto revDist = refSeq->compareKmersRevComp(*seqKmer);
						if (forDist.first < revDist.first) {
							if(revDist.first > best){
								best = revDist.first;
								revWinners = 1;
								forwardWinners = 0;
							}
						} else {
							if(forDist.first > best){
								best = forDist.first;
								revWinners = 0;
								forwardWinners = 1;
							}
						}
					}
				}else{
					for (const auto & refSeq : refKmerReads) {
						auto forDist = refSeq->compareKmers(*seqKmer);
						auto revDist = refSeq->compareKmersRevComp(*seqKmer);
						if (forDist.first < revDist.first) {
							++revWinners;
						} else {
							++forwardWinners;
						}
					}
				}
				if (revWinners > forwardWinners) {
					seq.reverseComplementRead(true, true);
				}
				{
					std::lock_guard<std::mutex> lock(outMut);
					reader.write(seq);
				}
			}
		};

		njh::concurrent::runVoidFunctionThreaded(reOrient, numThreads);

	} else {
		seqInfo compareInfo;
		if ("" != setUp.pars_.seq_) {
			compareInfo = setUp.pars_.seqObj_.seqBase_;
		} else {
			SeqIO reader(setUp.pars_.ioOptions_);
			reader.openIn();
			seqInfo seq;
			bool success = reader.readNextRead(seq);
			if (success) {
				compareInfo = seq;
			} else {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error couldn't read any seqs from "
						<< setUp.pars_.ioOptions_.firstName_ << "\n";
				throw std::runtime_error { ss.str() };
			}
		}
		std::unique_ptr<seqWithKmerInfo> compareSeq = std::make_unique<
				seqWithKmerInfo>(compareInfo);
		compareSeq->setKmers(setUp.pars_.colOpts_.kmerOpts_.kLength_, true);
		SeqIO reader(setUp.pars_.ioOptions_);
		reader.openIn();
		reader.openOut();
		seqInfo seq;
		while (reader.readNextRead(seq)) {
			auto seqKmer = std::make_unique<seqWithKmerInfo>(seq);
			seqKmer->setKmers(setUp.pars_.colOpts_.kmerOpts_.kLength_, true);
			auto forDist = compareSeq->compareKmers(*seqKmer);
			auto revDist = compareSeq->compareKmersRevComp(*seqKmer);
			if (forDist.first < revDist.first) {
				seqKmer->seqBase_.reverseComplementRead(true);
			}
			reader.openWrite(seqKmer);
		}
	}

	return 0;
}



int seqUtilsModRunner::collapseToUniqueWithInMetaField(
		const njh::progutils::CmdArgs & inputCommands) {
	std::string fieldStr = "";
	bfs::path nameFilename = "";
	bool renameToField = false;

	seqUtilsModSetUp setUp(inputCommands);
	setUp.processDefaultReader(true);
	if (setUp.pars_.ioOptions_.out_.outFilename_ == "out") {
		setUp.pars_.ioOptions_.out_.outFilename_ = njh::files::prependFileBasename(setUp.pars_.ioOptions_.firstName_, "unique_");
	}
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(renameToField, "--renameToField",
			"Rename the collapse sequence to the field collapsed with");
	setUp.setOption(nameFilename, "--names",
			"A file to print the names that were collapsed to the main cluster",
			true);
	setUp.setOption(fieldStr, "--field", "The meta field to collapse within", true);
	setUp.finishSetUp(std::cout);
	VecStr fields = njh::tokenizeString(fieldStr, ",");
	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	reader.openOut();
	std::unordered_map<std::string, std::vector<identicalCluster>> ans;
	seqInfo seq;
	uint32_t seqCount = 0;
	while(reader.readNextRead(seq)) {
		readVec::handelLowerCaseBases(seq, setUp.pars_.ioOptions_.lowerCaseBases_);
		if(setUp.pars_.ioOptions_.removeGaps_){
			seq.removeGaps();
		}
		++seqCount;
		bool found = false;
		MetaDataInName seqMeta(seq.name_);
		VecStr metaFields;
		for(const auto & field : fields){
			metaFields.emplace_back(seqMeta.getMeta(field));
		}
		std::string field = njh::conToStr(metaFields, "-");
		for (auto &cIter : ans[field]) {
			if (cIter.seqBase_.seq_ == seq.seq_) {
				cIter.addRead(seq);
				found = true;
				break;
			}
		}
		if (!found) {
			ans[field].emplace_back(seq);
		}
	}
	for (auto &byField : ans) {
		if(renameToField){
			njh::sort(byField.second);
		}
		uint32_t withinFieldCount = 0;
		for(auto & cIter : byField.second){
			if(renameToField){
				cIter.seqBase_.name_ = njh::pasteAsStr(byField.first,".",withinFieldCount);
			}
			cIter.updateName();
			++withinFieldCount;
		}
		reader.write(byField.second);
	}

	if (setUp.pars_.verbose_) {
		std::cout << "Collapsed " << seqCount << " reads to " << ans.size()
				<< std::endl;
	}
	if (setUp.pars_.verbose_) {
		setUp.logRunTime(std::cout);
	}
	if (nameFilename != "") {
		std::ofstream out;
		openTextFile(out, nameFilename, ".tab.txt", setUp.pars_.ioOptions_.out_);
		out << "name\tnumber\treads\n";
		for (auto &byField : ans) {
			for (const auto & read : byField.second) {
				out << read.seqBase_.name_ << "\t" << read.reads_.size() << "\t"
						<< vectorToString(readVec::getNames(read.reads_), ",") << std::endl;
			}
		}
	}
	return 0;
}


int seqUtilsModRunner::collapseToUnique(
		const njh::progutils::CmdArgs & inputCommands) {
	bfs::path nameFilename = "";

	seqUtilsModSetUp setUp(inputCommands);
	setUp.processDefaultReader(true);
	if (setUp.pars_.ioOptions_.out_.outFilename_ == "out") {
		setUp.pars_.ioOptions_.out_.outFilename_ = njh::files::prependFileBasename(setUp.pars_.ioOptions_.firstName_, "unique_");
	}
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(nameFilename, "--names",
			"A file to print the names that were collapsed to the main cluster");
	if("" == nameFilename){
		nameFilename = njh::files::nameAppendBeforeExt(setUp.pars_.ioOptions_.out_.outFilename_, "_names");
	}
	setUp.finishSetUp(std::cout);

	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	reader.openOut();

	OutOptions nameOutOpts(bfs::path(nameFilename).replace_extension("tab.txt"), ".tab.txt");
	nameOutOpts.transferOverwriteOpts(setUp.pars_.ioOptions_.out_);
	OutputStream out(nameOutOpts);

	std::vector<identicalCluster> ans;
	seqInfo seq;
	uint32_t seqCount = 0;
	while(reader.readNextRead(seq)) {
		readVec::handelLowerCaseBases(seq, setUp.pars_.ioOptions_.lowerCaseBases_);
		if(setUp.pars_.ioOptions_.removeGaps_){
			seq.removeGaps();
		}
		++seqCount;
		bool found = false;
		for (auto &cIter : ans) {
			if (cIter.seqBase_.seq_ == seq.seq_) {
				cIter.addRead(seq);
				found = true;
				break;
			}
		}
		if (!found) {
			ans.emplace_back(seq);
		}
	}
	for (auto &cIter : ans) {
		//cIter.seqBase_.name_ = cIter.getStubName(false);
		cIter.updateName();
	}

	reader.write(ans);
	if (setUp.pars_.verbose_) {
		std::cout << "Collapsed " << seqCount << " reads to " << ans.size()
				<< std::endl;
	}
	if (setUp.pars_.verbose_) {
		setUp.logRunTime(std::cout);
	}
	out << "name\tnumber\treads\n";
	for (const auto & read : ans) {
		out << read.seqBase_.name_ << "\t" << read.reads_.size() << "\t"
				<< vectorToString(readVec::getNames(read.reads_), ",") << std::endl;
	}

	return 0;
}

int seqUtilsModRunner::translate(
		const njh::progutils::CmdArgs & inputCommands) {
	seqUtilsModSetUp setUp(inputCommands);
	bool complement = false;
	bool reverse = false;
	uint64_t start = 0;
	bool trimLastStopCodon = false;
	setUp.setOption(trimLastStopCodon, "--trimLastStopCodon", "Trim Last Stop Codon");

	setUp.setUpTranslate(start, complement, reverse);

	if (setUp.pars_.ioOptions_.firstName_ == "") {
		seqInfo tempRead("temp", setUp.pars_.seq_);
		tempRead.translate(complement, reverse, start);
		if(trimLastStopCodon && len(tempRead.seq_) > 0 && '*' == tempRead.seq_.back()){
			tempRead.trimBack(len(tempRead) - 1);
		}
		std::cout << tempRead.seq_ << std::endl;
	} else {
		SeqIO reader(setUp.pars_.ioOptions_);
		reader.openIn();
		reader.openOut();
		seqInfo seq;
		while(reader.readNextRead(seq)){
			if(setUp.pars_.ioOptions_.removeGaps_){
				seq.removeGaps();
			}
			seq.translate(complement, reverse, start);
			if(trimLastStopCodon && len(seq.seq_) > 0 && '*' == seq.seq_.back()){
				seq.trimBack(len(seq) - 1);
			}
			reader.write(seq);
		}
	}
	return 0;
}

int seqUtilsModRunner::compSeq(const njh::progutils::CmdArgs & inputCommands) {
	seqUtilsModSetUp setUp(inputCommands);
	std::string seqType = "DNA";
	if (!setUp.processDefaultReader(false)) {
		setUp.processSeq(true);
	} else {
		if (setUp.pars_.ioOptions_.out_.outFilename_ == "out") {
			setUp.pars_.ioOptions_.out_.outFilename_ = njh::files::prependFileBasename(njh::files::removeExtension(setUp.pars_.ioOptions_.firstName_), "comp_");
		}
	}
	setUp.setOption(seqType, "--seqType", "SequenceType");
	setUp.finishSetUp(std::cout);

	if (setUp.pars_.ioOptions_.firstName_.empty()) {
		std::cout << seqUtil::complement(setUp.pars_.seq_, seqType)
							<< std::endl;
	} else {
		SeqIO reader(setUp.pars_.ioOptions_);
		reader.openIn();
		reader.openOut();
		seqInfo read;
		while (reader.readNextRead(read)) {
			read.complementRead();
			reader.out_.writeNoCheck(read);
		}
	}
	return 0;
}
int seqUtilsModRunner::revCompSeq(const njh::progutils::CmdArgs & inputCommands) {
	seqUtilsModSetUp setUp(inputCommands);
  std::string seqType = "DNA";
  bool mark = false;
  setUp.setOption(mark, "--mark", "Mark the sequence name with a _Comp to indicate it's in the reverse complement now");
	setUp.setUpComplementSeq(seqType);
	if (setUp.pars_.ioOptions_.firstName_ == "") {
		std::cout << seqUtil::reverseComplement(setUp.pars_.seq_, seqType)
				<< std::endl;
	} else {
		SeqIO reader(setUp.pars_.ioOptions_);
		reader.openIn();
		reader.openOut();
		seqInfo read;
		while (reader.readNextRead(read)) {
			read.reverseComplementRead(mark, true);
			reader.out_.writeNoCheck(read);
		}
	}
  return 0;
}



int seqUtilsModRunner::appendReads(
		const njh::progutils::CmdArgs & inputCommands) {
	bool degen = false;
	seqSetUp setUp(inputCommands);
	setUp.setOption(degen, "--degen", "Create all possible degen strings if input seq has degenerative bases in it");
	setUp.pars_.ioOptions_.out_.outFilename_ = "";
	setUp.processDefaultReader(true);
	if("" == setUp.pars_.ioOptions_.out_.outFilename_){
		setUp.pars_.ioOptions_.out_.outFilename_ = njh::files::prependFileBasename(setUp.pars_.ioOptions_.firstName_, "appended_");
	}
	setUp.processSeq(true);
	setUp.processVerbose();
	setUp.finishSetUp(std::cout);
	VecStr ans;
	if(degen){
		ans = createDegenStrs(setUp.pars_.seq_);
	}else{
		ans = {setUp.pars_.seq_};
	}
	if (setUp.pars_.verbose_) {
		printVector(ans, "\n");
	}
	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	reader.openOut();
	seqInfo seq;
	while (reader.readNextRead(seq)) {
		if (len(ans) > 1) {
			uint32_t subCount = 0;
			for (const auto & str : ans) {
				auto copyObj = seq;
				copyObj.append(str);
				if(0 != subCount ){
					copyObj.name_.append("_" + estd::to_string(subCount));
				}
				++subCount;
				reader.write(copyObj);
			}
		} else {
			seq.append(ans.front());
			reader.write(seq);
		}
	}
	return 0;
}

int seqUtilsModRunner::prependReads(
		const njh::progutils::CmdArgs & inputCommands) {
	bool degen = false;
	seqSetUp setUp(inputCommands);
	setUp.setOption(degen, "--degen", "Create all possible degen strings if input seq has degenerative bases in it");
	setUp.pars_.ioOptions_.out_.outFilename_ = "";
	setUp.processDefaultReader(true);
	if("" == setUp.pars_.ioOptions_.out_.outFilename_){
		setUp.pars_.ioOptions_.out_.outFilename_ = njh::files::prependFileBasename(setUp.pars_.ioOptions_.firstName_, "prepended_");
	}
	setUp.processSeq(true);
	setUp.processVerbose();
	setUp.finishSetUp(std::cout);
	VecStr ans ;
	if(degen){
		ans = createDegenStrs(setUp.pars_.seq_);
	}else{
		ans = {setUp.pars_.seq_};
	}
	if (setUp.pars_.verbose_) {
		printVector(ans, "\n");
	}
	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	reader.openOut();
	seqInfo seq;
	while (reader.readNextRead(seq)) {
		if (len(ans) > 1) {
			uint32_t subCount = 0;
			for (const auto & str : ans) {
				auto copyObj = seq;
				copyObj.prepend(str);
				if(0 != subCount ){
					copyObj.name_.append("_" + estd::to_string(subCount));
				}
				++subCount;
				reader.write(copyObj);
			}
		} else {
			seq.prepend(ans.front());
			reader.write(seq);
		}
	}
	return 0;
}

int seqUtilsModRunner::appendNames(
		const njh::progutils::CmdArgs & inputCommands) {
	std::string suffix = "";
	seqSetUp setUp(inputCommands);
	setUp.pars_.ioOptions_.out_.outFilename_ = "";
	setUp.processDefaultReader(true);
	if("" == setUp.pars_.ioOptions_.out_.outFilename_){
		setUp.pars_.ioOptions_.out_.outFilename_ = njh::files::prependFileBasename(setUp.pars_.ioOptions_.firstName_, "appended_");
	}
	setUp.setOption(suffix, "--suffix", "Suffix to add to names", true);
	setUp.processVerbose();
	setUp.finishSetUp(std::cout);

	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	reader.openOut();
	seqInfo seq;
	while (reader.readNextRead(seq)) {
		readVec::handelLowerCaseBases(seq, setUp.pars_.ioOptions_.lowerCaseBases_);
		seq.name_ += suffix;
		reader.write(seq);
	}
	return 0;
}

int seqUtilsModRunner::prependNames(
		const njh::progutils::CmdArgs & inputCommands) {
	std::string prefix = "";
	seqSetUp setUp(inputCommands);
	setUp.pars_.ioOptions_.out_.outFilename_ = "";
	setUp.processDefaultReader(true);
	if("" == setUp.pars_.ioOptions_.out_.outFilename_){
		setUp.pars_.ioOptions_.out_.outFilename_ = njh::files::prependFileBasename(setUp.pars_.ioOptions_.firstName_, "prepended_");
	}
	setUp.setOption(prefix, "--prefix", "Prefix to add to names", true);
	setUp.processVerbose();
	setUp.finishSetUp(std::cout);

	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	reader.openOut();
	seqInfo seq;
	while (reader.readNextRead(seq)) {
		readVec::handelLowerCaseBases(seq, setUp.pars_.ioOptions_.lowerCaseBases_);
		seq.name_ = prefix + seq.name_;
		reader.write(seq);
	}
	return 0;
}





int seqUtilsModRunner::renameSeqNameToUniqueNames(const njh::progutils::CmdArgs & inputCommands) {
	// parameters
	std::string sortBy = "none";
	auto outOptsKey = TableIOOpts::genTabFileOut("", true);
	seqUtilsModSetUp setUp(inputCommands);
	setUp.setOption(outOptsKey.out_.outFilename_, "--keyOut", "A filename to write a key for the original name to");
	setUp.processVerbose();
	setUp.processDefaultReader(true);
	if (setUp.pars_.ioOptions_.out_.outFilename_ == "out") {
		setUp.pars_.ioOptions_.out_.outFilename_ = njh::files::prependFileBasename(njh::files::removeExtension(setUp.pars_.ioOptions_.firstName_), "renamed_");
	}
	setUp.setOption(sortBy, "--sortBy", "SortOption");
	outOptsKey.out_.overWriteFile_ = setUp.pars_.ioOptions_.out_.overWriteFile_;
	setUp.finishSetUp(std::cout);


	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	reader.openOut();

	auto inReads = reader.in_.readAllReads<seqInfo>();
	if (sortBy != "none") {
		readVecSorter::sortReadVectorSimple(inReads, sortBy);
	}

	VecStr originalNames = readVec::getNames(inReads);
	std::unordered_map<std::string, uint32_t> nameCounts;
	for(const auto & seqInfo : inReads){
		++nameCounts[seqInfo.name_];
	}

	if(nameCounts.size() != originalNames.size()){
		std::unordered_map<std::string, uint32_t> reNamedCounts;
		//there are non-unique names so let's rename;
		for(auto & seqInfo : inReads){
			if(nameCounts[seqInfo.name_] >1){
				std::string outName = seqInfo.name_;
				seqInfo.name_.append(njh::pasteAsStr(".", njh::leftPadNumStr(reNamedCounts[seqInfo.name_], nameCounts[seqInfo.name_])));
				++reNamedCounts[outName];
			}
		}
	}
	VecStr newNames = readVec::getNames(inReads);
	table nameKey(VecStr{"originalName", "newName"});
	for(const auto pos : iter::range(newNames.size())){
		nameKey.addRow(originalNames[pos], newNames[pos]);
	}
	if("" != outOptsKey.out_.outFilename_){
		nameKey.outPutContents(outOptsKey);
	}
	reader.openWrite(inReads);


	if(setUp.pars_.verbose_){
		setUp.logRunTime(std::cout);
	}
	return 0;
}


int seqUtilsModRunner::renameIDs(const njh::progutils::CmdArgs & inputCommands) {
  // parameters
  std::string sortBy = "none", stub = "seq";
  bool keepChimeraFlag = false;
  bool keepComplementFlag = false;
  auto outOptsKey = TableIOOpts::genTabFileOut("", true);
  bfs::path keyIn = "";
  seqUtilsModSetUp setUp(inputCommands);
  setUp.setOption(keyIn, "--keyIn", "A file with a key to rename seqs with, no header, first column is old name, second column is the new name");
  setUp.setOption(outOptsKey.out_.outFilename_, "--keyOut", "A filename to write a key for the original name to");
  setUp.setOption(keepComplementFlag, "--keepComplementFlag", "Keep any reads marked with _Comp");
  setUp.processVerbose();
  setUp.setUpRenameIDs(stub, sortBy, keepChimeraFlag);
  outOptsKey.out_.overWriteFile_ = setUp.pars_.ioOptions_.out_.overWriteFile_;
  SeqIO reader(setUp.pars_.ioOptions_);

	reader.openIn();
	reader.openOut();

	if("" != keyIn){
		table nameKey(keyIn, "\t", false);
		if(2 != nameKey.nCol()){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "name key has to be 2 columns no header, 1) old name, 2) new name, not size of : " << nameKey.nCol() << "\n";
			throw std::runtime_error{ss.str()};
		}
		std::unordered_map<std::string, std::string> nameKeyMap;
		for(const auto & row : nameKey){
			if(njh::in(row[0], nameKeyMap)){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "already have name: " << row[0] << "\n";
				throw std::runtime_error{ss.str()};
			}
			nameKeyMap[row[0]] = row[1];
		}
		seqInfo seq;
		while(reader.readNextRead(seq)){
			if(!njh::in(seq.name_, nameKeyMap)){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "couldn't find name: " << seq.name_ << " in key: " << keyIn << "\n";
				throw std::runtime_error{ss.str()};
			}
			seq.name_ = nameKeyMap[seq.name_];
			reader.write(seq);
		}
	} else {
		auto inReads = reader.in_.readAllReads<readObject>();
	  if (sortBy != "none") {
	    readVecSorter::sortReadVector(inReads, sortBy);
	  }

	  VecStr originalNames = readVec::getNames(inReads);
	  renameReadNames(inReads, stub, setUp.pars_.ioOptions_.processed_,
	                  keepChimeraFlag);
	  VecStr newNames = readVec::getNames(inReads);
	  table nameKey(VecStr{"originalName", "newName"});
	  for(const auto pos : iter::range(newNames.size())){
	  	nameKey.addRow(originalNames[pos], newNames[pos]);
	  }
	  if("" != outOptsKey.out_.outFilename_){
	  	nameKey.outPutContents(outOptsKey);
	  }
	  reader.openWrite(inReads);
	}

  if(setUp.pars_.verbose_){
  	setUp.logRunTime(std::cout);
  }
  return 0;
}


int seqUtilsModRunner::removeLowQualityBases(const njh::progutils::CmdArgs & inputCommands) {

  seqUtilsModSetUp setUp(inputCommands);
  int qualCutOff = 5;
  setUp.processDefaultReader(true);
  setUp.setOption(qualCutOff, "-qualCutoff", "QualityCutOff");
  setUp.finishSetUp(std::cout);
  SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.in_.readAllReads<readObject>();
  if (setUp.pars_.ioOptions_.out_.outFilename_ == "") {
    setUp.pars_.ioOptions_.out_.outFilename_ = "qualTrimmed";
  }
  readVec::allRemoveLowQualityBases(inReads, qualCutOff);
  reader.openWrite(inReads);
  return 0;
}

int seqUtilsModRunner::dereplicate(
		const njh::progutils::CmdArgs & inputCommands) {
	//std::cout << "1" << std::endl;
	seqUtilsModSetUp setUp(inputCommands);
	setUp.pars_.ioOptions_.out_.outFilename_ = "dereplicated";
	setUp.processDefaultReader(true);
	setUp.finishSetUp(std::cout);
	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.in_.readAllReads<readObject>();
	reader.openOut();
	for (auto & read : inReads) {
		for (uint32_t seqNum : iter::range(read.seqBase_.cnt_)) {
			read.seqBase_.name_ = read.getStubName(true) + "_"
					+ estd::to_string(seqNum);
			reader.out_.writeNoCheck(read);
		}
	}
	return 0;
}

int seqUtilsModRunner::changeLetterCase(
		const njh::progutils::CmdArgs & inputCommands) {
	uint32_t position = 0;
	uint32_t size = std::numeric_limits<uint32_t>::max();
	bool back = false;
	seqUtilsModSetUp setUp(inputCommands);
	setUp.setOption(position, "--position", "Position at which to start changing case");
	setUp.setOption(size,     "--size",     "The number of bases after position to change the case for");
	setUp.setOption(back,     "--back",     "Make size and position relative to the back of the sequence");

	setUp.pars_.ioOptions_.out_.outFilename_ = "";
	setUp.processDefaultReader(true);
	if("" == setUp.pars_.ioOptions_.out_.outFilename_){
		setUp.pars_.ioOptions_.out_.outFilename_ = njh::files::prependFileBasename(setUp.pars_.ioOptions_.firstName_, "changedCase_");
	}
	if("lower" != setUp.pars_.ioOptions_.lowerCaseBases_  &&
			"upper" != setUp.pars_.ioOptions_.lowerCaseBases_){
		setUp.failed_ = true;
		std::stringstream ss;
		ss << "--lower should be either lower or upper, not " << setUp.pars_.ioOptions_.lowerCaseBases_;
		setUp.addWarning(ss.str());
	}
	setUp.finishSetUp(std::cout);

	std::function<char(char)> modCaseFunc = [](char c ){
		return c;
	};
	if ("upper" == setUp.pars_.ioOptions_.lowerCaseBases_) {
		modCaseFunc = [](char c ) {
			return toupper(c);
		};
	}else if("lower" == setUp.pars_.ioOptions_.lowerCaseBases_){
		modCaseFunc = [](char c ) {
			return tolower(c);
		};
	}

	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	reader.openOut();
	seqInfo seq;
	while (reader.readNextRead(seq)) {
		if(back){
			auto beg = std::rbegin(seq.seq_);
			auto ending = std::rend(seq.seq_);
			if (position < len(seq)) {
				beg = std::rbegin(seq.seq_) + position;
				if (len(seq) - position > size) {
					ending = std::rbegin(seq.seq_) + position + size;
				}
			}
			std::transform(beg, ending, beg, modCaseFunc);
		}else{
			auto beg = std::begin(seq.seq_);
			auto ending = std::end(seq.seq_);
			if (position < len(seq)) {
				beg = std::begin(seq.seq_) + position;
				if (len(seq) - position > size) {
					ending = std::begin(seq.seq_) + position + size;
				}
			}
			std::transform(beg, ending, beg, modCaseFunc);
		}
		reader.out_.writeNoCheck(seq);
	}
	return 0;
}

int seqUtilsModRunner::inverseLetterCase(
		const njh::progutils::CmdArgs & inputCommands) {
	uint32_t position = 0;
	uint32_t size = std::numeric_limits<uint32_t>::max();
	bool back = false;
	seqUtilsModSetUp setUp(inputCommands);
	setUp.setOption(position, "--position", "Position at which to start changing case");
	setUp.setOption(size,     "--size",     "The number of bases after position to change the case for");
	setUp.setOption(back,     "--back",     "Make size and position relative to the back of the sequence");
	setUp.pars_.ioOptions_.out_.outFilename_ = "";
	setUp.processDefaultReader(true);
	if("" == setUp.pars_.ioOptions_.out_.outFilename_){
		setUp.pars_.ioOptions_.out_.outFilename_ = njh::files::prependFileBasename(setUp.pars_.ioOptions_.firstName_, "inversedCase_");
	}
	setUp.finishSetUp(std::cout);
	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	reader.openOut();
	seqInfo seq;
	while(reader.readNextRead(seq)){
		if(back){
			auto beg = std::rbegin(seq.seq_);
			auto ending = std::rend(seq.seq_);
			if (position < len(seq)) {
				beg = std::rbegin(seq.seq_) + position;
				if (len(seq) - position > size) {
					ending = std::rbegin(seq.seq_) + position + size;
				}
			}
			std::transform(beg, ending,
							beg, [](char c) {
								return std::islower(c) ? std::toupper(c) : std::tolower(c);
							});
		}else{
			auto beg = std::begin(seq.seq_);
			auto ending = std::end(seq.seq_);
			if (position < len(seq)) {
				beg = std::begin(seq.seq_) + position;
				if (len(seq) - position > size) {
					ending = std::begin(seq.seq_) + position + size;
				}
			}
			std::transform(beg, ending,
							beg, [](char c) {
								return std::islower(c) ? std::toupper(c) : std::tolower(c);
							});
		}

		reader.out_.writeNoCheck(seq);
	}
	return 0;
}


                    
} // namespace njhseq
