/*
 * pacbioExpRunner.cpp
 *
 *  Created on: Jan 25, 2015
 *      Author: nickhathaway
 */

#include "pacbioExp.hpp"

#include "elucidator/objects/seqContainers/refVariants.hpp"


namespace njhseq {
pacbioExpRunner::pacbioExpRunner()
    : njh::progutils::ProgramRunner(
          {
					 addFunc("profilePacbioReads", profilePacbioReads, false),
					 addFunc("createCssFromRawPacbio", createCssFromRawPacbio, false),
					 addFunc("testCssConsensusBuilding", testCssConsensusBuilding, false),
					 addFunc("checkCssReadsToRaw", checkCssReadsToRaw, false),
					 addFunc("indexErrorsRawToCcsProcessed", indexErrorsRawToCcsProcessed, false),
					 addFunc("indexErrorsRawToCcs", indexErrorsRawToCcs, false),
					 addFunc("processSnpTable", processSnpTable, false),
					 addFunc("testVariantCalling", testVariantCalling, false),
					 addFunc("printPacbioIDs", printPacbioIDs, false),
					 addFunc("getRoundsForPacbioIDs", getRoundsForPacbioIDs, false),
					 addFunc("simPacbioPerRead", simPacbioPerRead, false)
           },//
          "pacbioExp") {}

inline std::string getPacbioRunId(const std::string & str){
	auto toks = tokenizeString(str, "/");
	return toks.front();
}

inline uint64_t getPacbioGroupNumber(const std::string & str){
	auto toks = tokenizeString(str, "/");
	return estd::stou(toks[toks.size() - 2]);
}

inline std::string getPassRep(const std::string & str){
	auto toks = tokenizeString(str, "/");
	return toks.back();
}


template <typename T, typename FUNC>
std::string vectorToString(const std::vector<T>& vectorToConvert,
                           const std::string& delim, FUNC func) {
  if (vectorToConvert.empty()) {
    return "";
  }
  std::vector<std::string> tempVec(vectorToConvert.size());
  std::transform(vectorToConvert.begin(), vectorToConvert.end(), tempVec.begin(), func);
  return vectorToString(tempVec, delim);
}

int pacbioExpRunner::getRoundsForPacbioIDs(const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	setUp.processDefaultReader(true);
	if("out" == setUp.pars_.ioOptions_.out_.outFilename_){
		setUp.pars_.ioOptions_.out_.outFilename_ = bfs::path(setUp.pars_.ioOptions_.firstName_).replace_extension("").string() + "_RoundsForIds";
	}
	setUp.pars_.ioOptions_.out_.outExtention_ = ".txt";
	setUp.finishSetUp(std::cout);

	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	seqInfo seq;
	std::unordered_map<std::string, std::unordered_map<uint64_t, uint32_t>> roundCounts;

	while(reader.readNextRead(seq)){
		auto runId = getPacbioRunId(seq.name_);
		auto id = getPacbioGroupNumber(seq.name_);
		++roundCounts[runId][id];
	}
	std::ofstream outFile;
	openTextFile(outFile, setUp.pars_.ioOptions_.out_);
	auto keys = getVectorOfMapKeys(roundCounts);
	njh::sort(keys);
	outFile << "id\trounds\n";
	for(const auto key : keys){
		auto idKeys = getVectorOfMapKeys(roundCounts[key]);
		njh::sort(idKeys);
		for(const auto & idKey : idKeys){
			outFile << key << "/" << idKey << "\t" << roundCounts[key][idKey] << std::endl;
		}
	}
	return 0;
}

int pacbioExpRunner::printPacbioIDs(const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	setUp.processDefaultReader(true);
	if("out" == setUp.pars_.ioOptions_.out_.outFilename_){
		setUp.pars_.ioOptions_.out_.outFilename_ = bfs::path(setUp.pars_.ioOptions_.firstName_).replace_extension("").string() + "_ids";
	}
	setUp.pars_.ioOptions_.out_.outExtention_ = ".txt";
	setUp.finishSetUp(std::cout);

	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	seqInfo seq;
	std::unordered_map<std::string, std::vector<uint64_t>> ids;

	while(reader.readNextRead(seq)){
		auto runId = getPacbioRunId(seq.name_);
		auto id = getPacbioGroupNumber(seq.name_);
		if(njh::in(id, ids[runId])){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << " error, repeat id: " << id << "\n";
			throw std::runtime_error{ss.str()};
		}
		ids[runId].emplace_back(id);
	}
	std::ofstream outFile;
	openTextFile(outFile, setUp.pars_.ioOptions_.out_);
	auto keys = getVectorOfMapKeys(ids);
	njh::sort(keys);
	for(const auto & key : keys){
		for(const auto & id : ids[key]){
			outFile << key << "/" << id << "\n";
		}
	}
	return 0;
}

int pacbioExpRunner::testVariantCalling(const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	setUp.processDefaultReader(true);
	setUp.processAlignerDefualts();
	setUp.finishSetUp(std::cout);
	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<readObject>();


	std::vector<refVariants> rVars;
	for (const auto & read : inReads) {
		rVars.emplace_back(read.seqBase_);
	}
	uint64_t maxSize = 0;
	readVec::getMaxLength(inReads, maxSize);
	aligner alignerObj(maxSize, setUp.pars_.gapInfo_, setUp.pars_.scoring_,
			setUp.pars_.colOpts_.alignOpts_.countEndGaps_);
	alignerObj.processAlnInfoInput(setUp.pars_.alnInfoDirName_);

	for (const auto & readPos : iter::range(inReads.size())) {
		for (const auto & subReadPos : iter::range(inReads.size())) {
			if (readPos == subReadPos) {
				continue;
			}
			rVars[readPos].addVariant(inReads[subReadPos].seqBase_, alignerObj,
					false);
		}
	}
	std::ofstream outFile;
	openTextFile(outFile, setUp.pars_.ioOptions_.out_.outFilename_, ".txt", setUp.pars_.ioOptions_.out_);
	for(const auto & r : rVars){
		r.outPut(outFile);
	}
	alignerObj.alnHolder_.write(setUp.pars_.outAlnInfoDirName_);
	return 0;
}

int pacbioExpRunner::processSnpTable(const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	std::string filename = "";
	std::string outFilename = "";
	double flipDifference = 0.75;
	double freqOccurentCutOff = 0.20;
	setUp.setOption(filename, "--file", "Name of the snp file", true);
	setUp.setOption(outFilename, "--outFile", "Name of an out file, otherwise will print to stdout");
	setUp.setOption(flipDifference, "--fileDifference", "Difference in flipping occurrence to be counted as a flipping base");
	setUp.setOption(freqOccurentCutOff, "--freqOccurenceCutOff", "Cutoff for occurrence of the snp in the reads");
	setUp.processWritingOptions();
	setUp.finishSetUp(std::cout);
	table inTab(filename, "\t", true);
	table outTab(inTab.columnNames_);
	uint32_t forwFracPos = inTab.getColPos("forwFrac");
	uint32_t compFracPos = inTab.getColPos("compFrac");
	uint32_t repFractionPos = inTab.getColPos("repFraction");
	for (const VecStr & row : inTab) {
		if (std::abs(
				njh::lexical_cast<double>(row[forwFracPos])
						- njh::lexical_cast<double>(row[compFracPos])) > flipDifference
				&& njh::lexical_cast<double>(row[repFractionPos])
						> freqOccurentCutOff) {
			outTab.content_.emplace_back(row);
		}
	}

	outTab.outPutContents(
			TableIOOpts { OutOptions(outFilename, ".tab.txt", "tab", false,
					setUp.pars_.ioOptions_.out_.overWriteFile_, false), "\t",
					outTab.hasHeader_ });
	return 0;
}

int pacbioExpRunner::checkCssReadsToRaw(const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	std::string renameKeyFilename = "";
	setUp.setOption(renameKeyFilename, "--renameKeyFilename",
			"File with the first column is original pacbio name and the second is the renamed Mid name");
	setUp.processDefaultReader(true);
	setUp.processDirectoryOutputName(true);

	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<readObject>();

	std::unordered_map<uint64_t, uint32_t> nameCounts;
	std::unordered_map<uint64_t, std::vector<uint64_t>> readNamesPositions;
	for (const auto & readEnum : iter::enumerate(inReads)) {
		++nameCounts[getPacbioGroupNumber(readEnum.element.seqBase_.name_)];
		readNamesPositions[getPacbioGroupNumber(readEnum.element.seqBase_.name_)].emplace_back(
				readEnum.index);
	}
	table inTab { renameKeyFilename, "\t", true };
	std::unordered_map<std::string, std::vector<uint64_t>> outReads;
	std::unordered_map<std::string, std::vector<uint64_t>> passNumbers;
	table outInfoTab { VecStr { "midName", "rawReadNumber", "cssReadNumber",
			"medianPasses", "meanPasess", "minPasses", "maxPasses" } };
	for (const VecStr & row : inTab) {
		auto groupNum = getPacbioGroupNumber(row[0]);
		auto midName = njh::replaceString(row[1].substr(0, row[1].find(".")), "_Comp", "");
		addOtherVec(outReads[midName], readNamesPositions[groupNum]);
		passNumbers[midName].emplace_back(readNamesPositions[groupNum].size());
	}
	for (const auto & midReads : outReads) {
		std::ofstream midFile;
		openTextFile(midFile, setUp.pars_.directoryName_ + midReads.first, ".fastq",
				false, false);
		for (const auto & readPos : midReads.second) {
			inReads[readPos].seqBase_.outPutFastq(midFile);
		}
		auto stats = getStatsOnVec(passNumbers[midReads.first]);
		outInfoTab.content_.emplace_back(
				toVecStr(midReads.first, midReads.second.size(),
						passNumbers[midReads.first].size(), stats["median"], stats["mean"],
						stats["min"], stats["max"]));
	}
	outInfoTab.sortTable("midName", false);
	outInfoTab.outPutContents(
			TableIOOpts(OutOptions(setUp.pars_.directoryName_ + "info", ".tab.txt"), "\t", outInfoTab.hasHeader_));
	return 0;
}

int pacbioExpRunner::indexErrorsRawToCcsProcessed(const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	std::string renameKeyFilename = "";
	std::string ccsReads = "";
	setUp.setOption(ccsReads, "--ccsReadsFile", "File of the Renamed Circular Consensus Sequence(ccs) reads",true);
	setUp.setOption(renameKeyFilename, "--renameKeyFilename",
			"File with the first column is original pacbio name and the second is the renamed Mid name", true);
	setUp.processDefaultReader(true);
	setUp.processDirectoryOutputName(true);
	setUp.processSeq(true);
	setUp.pars_.gap_ = "5,1";
	setUp.pars_.gapLeft_ = "0,0";
	setUp.pars_.gapRight_ = "0,0";
	setUp.processAlignerDefualts();
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<readObject>();

	uint64_t maxSize = 0;
	readVec::getMaxLength(inReads, maxSize);
	std::unordered_map<uint64_t, uint32_t> nameCounts;
	std::unordered_map<uint64_t, std::vector<uint64_t>> readNamesPositions;
	for (const auto & readEnum : iter::enumerate(inReads)) {
		++nameCounts[getPacbioGroupNumber(readEnum.element.seqBase_.name_)];
		readNamesPositions[getPacbioGroupNumber(readEnum.element.seqBase_.name_)].emplace_back(
				readEnum.index);
	}
	SeqIOOptions ccsOpts = SeqIOOptions::genFastqIn(ccsReads);
	SeqInput ccsReader(ccsOpts);
	ccsReader.openIn();
	auto ccsinReads = ccsReader.readAllReads<readObject>();

	table inTab { renameKeyFilename, "\t", true };
	std::unordered_map<std::string, std::string> nameConverter;
	for(const auto & row : inTab){
		nameConverter[row[1]] = row[0];
	}
	readVec::getMaxLength(ccsinReads, maxSize);
	std::unordered_map<std::string, uint64_t> ccsReadPositions;
	for(const auto & readPos : iter::range(ccsinReads.size())){
		ccsReadPositions[ccsinReads[readPos].seqBase_.name_] = readPos;
	}

	aligner alignerObj(maxSize, setUp.pars_.gapInfo_, setUp.pars_.scoring_,false);
	alignerObj.processAlnInfoInput(setUp.pars_.alnInfoDirName_);
	bfs::path snpDir = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar("snpDir", false));
	bfs::path alnDir = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar("alnDir", false));
	bfs::path alnNoInsertsDir = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar("alnDirNoInserts", false));
	uint32_t count = 0;
	for (const VecStr & row : inTab) {
		auto groupNum = getPacbioGroupNumber(row[0]);
		if(ccsReadPositions.find(row[1]) == ccsReadPositions.end()
				|| readNamesPositions.find(groupNum) == readNamesPositions.end()){
			continue;
		}
		++count;
		if(count % 1 == 0){
			std::cout << "On " << count << "\r";
			std::cout.flush();
		}

		//auto midName = row[1].substr(0, row[1].find("."));
		std::vector<seqWithKmerInfo> currentGroupReads;
		currentGroupReads.reserve(readNamesPositions[groupNum].size());
		uint32_t index = 1;
		for (const auto & pos : readNamesPositions[groupNum]) {
			currentGroupReads.emplace_back(inReads[pos].seqBase_);
			if (index % 2 == 0) {
				currentGroupReads.back().seqBase_.reverseComplementRead(true, true);
			}
			currentGroupReads.back().setKmers(7, true);
			++index;
		}
		seqWithKmerInfo ccsRead = seqWithKmerInfo {
				ccsinReads[ccsReadPositions[row[1]]].seqBase_ };
		ccsRead.setKmers(7, true);

		for (auto & read : currentGroupReads) {
			auto forDist = ccsRead.compareKmers(read);
			auto revDist = ccsRead.compareKmersRevComp(read);
			if (forDist.first < revDist.first) {
				read.seqBase_.reverseComplementRead(true, true);
			}
		}
		std::unordered_map<uint32_t, std::unordered_map<char, VecStr>> mismatches;
		std::ofstream outAlnFile;
		openTextFile(outAlnFile, njh::files::make_path(alnDir, ccsRead.seqBase_.name_).string(), ".fastq", false, false);
		std::ofstream outAlnNoInsertsFile;
		openTextFile(outAlnNoInsertsFile, njh::files::make_path(alnNoInsertsDir,  ccsRead.seqBase_.name_).string(), ".fastq", false, false);
		ccsRead.seqBase_.outPutFastq(outAlnNoInsertsFile);
		for (const auto & subReadPos : iter::range(currentGroupReads.size())) {
			const auto & subRead = currentGroupReads[subReadPos];
			alignerObj.alignCache(ccsRead, subRead, false);
			std::vector<uint32_t> alnPositons(len(alignerObj.alignObjectA_));
			njh::iota<uint32_t>(alnPositons, 0);
			auto alnA = alignerObj.alignObjectA_.seqBase_;
			auto alnB = alignerObj.alignObjectB_.seqBase_;
			for(const auto & pos : iter::reversed(alnPositons)){
				if(alnA.seq_[pos] == '-'){
					alnA.removeBase(pos);
					alnB.removeBase(pos);
				}
			}
			//alnA.outPutFastq(outAlnNoInsertsFile);
			alnB.outPutFastq(outAlnNoInsertsFile);
			alignerObj.alignObjectA_.seqBase_.outPutFastq(outAlnFile);
			alignerObj.alignObjectB_.seqBase_.outPutFastq(outAlnFile);
			//count gaps and mismatches and get identity
			alignerObj.profilePrimerAlignment(ccsRead, subRead);
			for (const auto & m : alignerObj.comp_.distances_.mismatches_) {
				mismatches[m.second.refBasePos][m.second.seqBase].emplace_back(
						subRead.seqBase_.name_);
			}
		}
		table misTab { VecStr { "refPos", "refBase", "ccsPos", "ccsBase", "rawBase", "ccsQual", "misType", "freq",
				"fraction","ccsName", "seqs" } };
		double readTotal =
				std::accumulate(currentGroupReads.begin(), currentGroupReads.end(), 0.0,
						[](double res,const seqWithKmerInfo & seq) {return res + seq.seqBase_.cnt_;});
		alignerObj.alignCache(setUp.pars_.seqObj_, ccsRead, false);

		for (const auto & m : mismatches) {
			for (const auto & seqM : m.second) {
				auto ccsAlnPos = alignerObj.getAlignPosForSeqBPos(m.first);
				auto refPos = alignerObj.getSeqPosForAlnAPos(ccsAlnPos);
				misTab.content_.emplace_back(
						toVecStr(refPos, alignerObj.alignObjectA_.seqBase_.seq_[ccsAlnPos],
								m.first, ccsRead.seqBase_.seq_[m.first], seqM.first,
								ccsRead.seqBase_.qual_[m.first],
								(mismatch::isMismatchTransition(ccsRead.seqBase_.seq_[m.first],
										seqM.first) ? "transition" : "transversion"),
								seqM.second.size(), seqM.second.size() / readTotal,
								ccsRead.seqBase_.name_, vectorToString(seqM.second, ",")));
			}
		}
		misTab.sortTable("rawBase", false);
		misTab.sortTable("ccsPos", false);
		misTab.outPutContents(
				TableIOOpts(OutOptions(njh::files::make_path(snpDir,  ccsRead.seqBase_.name_).string(), ".tab.txt"), "\t", misTab.hasHeader_));
		//addOtherVec(outReads[midName], readNamesPositions[groupNum]);
		//passNumbers[midName].emplace_back(readNamesPositions[groupNum].size());
	}
	std::cout << std::endl;
	if (setUp.pars_.writingOutAlnInfo_) {
		alignerObj.alnHolder_.write(setUp.pars_.outAlnInfoDirName_);
	}
	return 0;
}

int pacbioExpRunner::indexErrorsRawToCcs(const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	std::string renameKeyFilename = "";
	std::string ccsReads = "";
	std::string namesFilename = "";
	setUp.setOption(ccsReads, "--ccsReadsFile", "File of the Renamed Circular Consensus Sequence(ccs) reads",true);
	setUp.setOption(renameKeyFilename, "--renameKeyFilename",
			"File with the first column is original pacbio name and the second is the renamed Mid name");
	setUp.setOption(namesFilename, "--names", "A file where the first column is the names of the reads to index");
	setUp.processDefaultReader(true);
	setUp.processDirectoryOutputName(true);
	setUp.processSeq(true);
	setUp.pars_.gap_ = "5,1";
	setUp.pars_.gapLeft_ = "0,0";
	setUp.pars_.gapRight_ = "0,0";
	setUp.processAlignerDefualts();
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<readObject>();

	uint64_t maxSize = 0;
	readVec::getMaxLength(inReads, maxSize);
	std::unordered_map<uint64_t, uint32_t> nameCounts;
	std::unordered_map<uint64_t, std::vector<uint64_t>> readNamesPositions;
	for (const auto & readEnum : iter::enumerate(inReads)) {
		++nameCounts[getPacbioGroupNumber(readEnum.element.seqBase_.name_)];
		readNamesPositions[getPacbioGroupNumber(readEnum.element.seqBase_.name_)].emplace_back(
				readEnum.index);
	}
	SeqIOOptions ccsOpts = SeqIOOptions::genFastqIn(ccsReads);
	SeqInput ccsReader(ccsOpts);
	auto ccsinReads = ccsReader.readAllReads<readObject>();
	std::unordered_map<std::string, std::string> cssNameToMidName;
	if(namesFilename != ""){
		table inTab { renameKeyFilename, "\t", true };
		std::unordered_map<std::string, std::string> nameConverter;
		for(const VecStr & row : inTab){
			nameConverter[njh::replaceString(row[0], "_Comp", "")] = row[1];
		}
		cssNameToMidName = nameConverter;
		table inNames(namesFilename);
		VecStr names;
		for(const VecStr & row : inNames){
			names.emplace_back(row[0]);
		}
		std::vector<readObject> specificReads;
		for(const auto & read : ccsinReads){
			if(njh::in(nameConverter[read.seqBase_.name_], names)){
				specificReads.emplace_back(read);
			}
		}
		ccsinReads = specificReads;
	}
	readVec::getMaxLength(ccsinReads, maxSize);
	std::unordered_map<std::string, uint64_t> ccsReadPositions;
	for(const auto & readPos : iter::range(ccsinReads.size())){
		ccsReadPositions[ccsinReads[readPos].seqBase_.name_] = readPos;
	}

	aligner alignerObj(maxSize, setUp.pars_.gapInfo_, setUp.pars_.scoring_,false);
	alignerObj.processAlnInfoInput(setUp.pars_.alnInfoDirName_);
	bfs::path snpDir = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar("snpDir", false));
	bfs::path alnDir = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar("alnDir", false));
	bfs::path alnNoInsertsDir = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar("alnDirNoInserts", false));
	bfs::path consensusBaseCountsdir = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar("consensusBaseCounts", false));
	uint32_t count = 0;
	seqWithKmerInfo refRead {setUp.pars_.seqObj_.seqBase_, 7, true};
	std::ofstream passInfoFile;
	openTextFile(passInfoFile, setUp.pars_.directoryName_ + "pasInfo", ".tab.txt", false, false);
	passInfoFile << "group\tpassNumber" << std::endl;
	for (const auto & read : ccsinReads) {

		seqWithKmerInfo ccsRead = seqWithKmerInfo {read.seqBase_};
		auto groupNum = getPacbioGroupNumber(ccsRead.seqBase_.name_);
		++count;
		std::cout << "On " << count << " of " << ccsinReads.size() << " ("
				<< roundDecPlaces(count / static_cast<double>(ccsinReads.size()),2) * 100 << "%)"
				<< "\r";
		std::cout.flush();
		std::vector<seqWithKmerInfo> currentGroupReads;
		currentGroupReads.reserve(readNamesPositions[groupNum].size());
		uint32_t index = 1;
		for (const auto & pos : readNamesPositions[groupNum]) {
			currentGroupReads.emplace_back(inReads[pos].seqBase_);
			if (index % 2 == 0) {
				currentGroupReads.back().seqBase_.reverseComplementRead(true, true);
			}
			currentGroupReads.back().setKmers(7, true);
			++index;
		}
		passInfoFile << groupNum << "\t" << currentGroupReads.size() << std::endl;
		ccsRead.setKmers(7, true);
		{
			auto forDist = ccsRead.compareKmers(refRead);
			auto revDist = ccsRead.compareKmersRevComp(refRead);
			if (forDist.first < revDist.first) {
				ccsRead.seqBase_.reverseComplementRead(true, true);
				ccsRead.setKmers(7, true);
			}
		}
		for (auto & groupRead : currentGroupReads) {
			auto forDist = ccsRead.compareKmers(groupRead);
			auto revDist = ccsRead.compareKmersRevComp(groupRead);
			if (forDist.first < revDist.first) {
				groupRead.seqBase_.reverseComplementRead(true, true);
			}
		}
		std::unordered_map<uint32_t, std::unordered_map<char, VecStr>> mismatches;
		std::ofstream outAlnFile;
		openTextFile(outAlnFile, njh::files::make_path(alnDir, estd::to_string(groupNum)).string(), ".fastq", false, false);
		std::ofstream outAlnNoInsertsFile;
		openTextFile(outAlnNoInsertsFile, njh::files::make_path(alnNoInsertsDir, estd::to_string(groupNum)).string(), ".fastq", false, false);
		ccsRead.seqBase_.outPutFastq(outAlnNoInsertsFile);
		std::vector<uint32_t> basePerPosCounts(len(ccsRead.seqBase_),0);
		std::ofstream consensusBaseCountFile;
		openTextFile(consensusBaseCountFile, njh::files::make_path(consensusBaseCountsdir, estd::to_string(groupNum)).string(), ".tab.txt", false, false);
		std::vector<charCounter> counters(len(ccsRead.seqBase_));
		for (const auto & subReadPos : iter::range(currentGroupReads.size())) {
			const auto & subRead = currentGroupReads[subReadPos];
			alignerObj.alignCache(ccsRead, subRead, false);
			std::vector<uint32_t> alnPositons(len(alignerObj.alignObjectA_));
			njh::iota<uint32_t>(alnPositons, 0);
			auto alnA = alignerObj.alignObjectA_.seqBase_;
			auto alnB = alignerObj.alignObjectB_.seqBase_;
			for(const auto & pos : iter::reversed(alnPositons)){
				if(alnA.seq_[pos] == '-'){
					alnA.removeBase(pos);
					alnB.removeBase(pos);
				}
			}
			for(const auto & pos : iter::range(len(alnB))){
				if(alnB.seq_[pos] != '-'){
					++basePerPosCounts[pos];
				}
			}
			//alnA.outPutFastq(outAlnNoInsertsFile);
			alnB.outPutFastq(outAlnNoInsertsFile);
			for(const auto & cPos : iter::range(len(alnB))){
				counters[cPos].increaseCountOfBase(alnB.seq_[cPos]);
			}
			alignerObj.alignObjectA_.seqBase_.outPutFastq(outAlnFile);
			alignerObj.alignObjectB_.seqBase_.outPutFastq(outAlnFile);
			//count gaps and mismatches and get identity
			alignerObj.profilePrimerAlignment(ccsRead, subRead);
			for (const auto & m : alignerObj.comp_.distances_.mismatches_) {
				mismatches[m.second.refBasePos][m.second.seqBase].emplace_back(
						subRead.seqBase_.name_);
			}
		}
		table misTab;
		if(namesFilename != ""){
			misTab = table{ VecStr { "refPos", "refBase", "ccsPos", "ccsBase", "rawBase", "ccsQual", "misType", "freq",
										"readfraction", "repFraction","ccsName", "midName", "MID","compCount", "compFrac", "forwCount", "forwFrac","passNumber", "seqs" } };
		}else{
			 misTab = table{ VecStr { "refPos", "refBase", "ccsPos", "ccsBase", "rawBase", "ccsQual", "misType", "freq",
							"readfraction", "repFraction","ccsName","compCount", "compFrac", "forwCount", "forwFrac","passNumber", "seqs" } };
		}
		double readTotal =
				std::accumulate(currentGroupReads.begin(), currentGroupReads.end(), 0.0,
						[](double res,const seqWithKmerInfo & seq) {return res + seq.seqBase_.cnt_;});
		alignerObj.alignCache(setUp.pars_.seqObj_, ccsRead, false);
		consensusBaseCountFile << "refPos\tccsPos\tgroupNum\ttotalReads\ttotalBases\tA\tC\tG\tT\n";
		for (const auto & cPos : iter::range(counters.size())) {
			auto ccsAlnPos = alignerObj.getAlignPosForSeqBPos(cPos);
			auto refPos = alignerObj.getSeqPosForAlnAPos(ccsAlnPos);
			counters[cPos].alphabet_ = {'A', 'C', 'G', 'T'};
			counters[cPos].setFractions();

			consensusBaseCountFile << refPos << "\t" << cPos << "\t" << groupNum
					<< "\t" << currentGroupReads.size() << "\t"
					<< counters[cPos].getTotalCount() << "\t"
					<< counters[cPos].fractions_['A'] << "\t"
					<< counters[cPos].fractions_['C'] << "\t"
					<< counters[cPos].fractions_['G'] << "\t"
					<< counters[cPos].fractions_['T'] << "\n";
		}
		for (const auto & m : mismatches) {
			for (const auto & seqM : m.second) {
				uint32_t compCount = 0;
				uint32_t forwCount = 0;
				for (const auto & name : seqM.second) {
					if (name.find("_Comp") != std::string::npos) {
						++compCount;
					} else {
						++forwCount;
					}
				}
				double totalNameCount = compCount + forwCount;
				auto ccsAlnPos = alignerObj.getAlignPosForSeqBPos(m.first);
				auto refPos = alignerObj.getSeqPosForAlnAPos(ccsAlnPos);
				auto getMid = [](std::string name){
					name = njh::replaceString(name, "_Comp", "");
					auto midPos = name.find("MID");
					return name.substr(midPos, name.find(".") - name.find("MID"));};
				if(namesFilename != ""){
					misTab.content_.emplace_back(
							toVecStr(refPos,
									alignerObj.alignObjectA_.seqBase_.seq_[ccsAlnPos], m.first,
									ccsRead.seqBase_.seq_[m.first], seqM.first,
									ccsRead.seqBase_.qual_[m.first],
									(mismatch::isMismatchTransition(
											ccsRead.seqBase_.seq_[m.first], seqM.first) ?
											"transition" : "transversion"), seqM.second.size(),
									seqM.second.size() / readTotal,
									seqM.second.size()
											/ static_cast<double>(basePerPosCounts[m.first]),
									getPacbioGroupNumber(ccsRead.seqBase_.name_),
									cssNameToMidName[njh::replaceString(ccsRead.seqBase_.name_, "_Comp", "")],
									getMid(cssNameToMidName[njh::replaceString(ccsRead.seqBase_.name_, "_Comp", "")]), compCount,
									compCount / totalNameCount, forwCount,
									forwCount / totalNameCount, currentGroupReads.size(),
									vectorToString(seqM.second, ",",
											[](const std::string & str) {return getPassRep(str);})));
				}else{
					misTab.content_.emplace_back(
							toVecStr(refPos, alignerObj.alignObjectA_.seqBase_.seq_[ccsAlnPos],
									m.first, ccsRead.seqBase_.seq_[m.first], seqM.first,
									ccsRead.seqBase_.qual_[m.first],
									(mismatch::isMismatchTransition(ccsRead.seqBase_.seq_[m.first],
											seqM.first) ? "transition" : "transversion"),
									seqM.second.size(), seqM.second.size() / readTotal,
									seqM.second.size()
											/ static_cast<double>(basePerPosCounts[m.first]),
									getPacbioGroupNumber(ccsRead.seqBase_.name_), compCount, compCount / totalNameCount,
									forwCount, forwCount / totalNameCount, currentGroupReads.size(),
									vectorToString(seqM.second, ",", [](const std::string & str){ return getPassRep(str);})));
				}

			}
		}
		misTab.sortTable("rawBase", false);
		misTab.sortTable("ccsPos", false);
		misTab.outPutContents(
				TableIOOpts(OutOptions(njh::files::make_path(snpDir, estd::to_string(groupNum)).string(), ".tab.txt"), "\t", misTab.hasHeader_));
	}
	std::cout << std::endl;
	if (setUp.pars_.writingOutAlnInfo_) {
		alignerObj.alnHolder_.write(setUp.pars_.outAlnInfoDirName_);
	}
	return 0;
}


int pacbioExpRunner::createCssFromRawPacbio(const njh::progutils::CmdArgs & inputCommands){
	seqSetUp setUp(inputCommands);
	setUp.pars_.ioOptions_.out_.outFilename_ = "out";
	setUp.processDefaultReader(true);
	setUp.processDirectoryOutputName(true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<readObject>();


	std::unordered_map<uint64_t, uint32_t> nameCounts;
	std::unordered_map<uint64_t, std::vector<uint64_t>> readNamesPositions;
	for(const auto & readEnum : iter::enumerate(inReads)){
		++nameCounts[getPacbioGroupNumber(readEnum.element.seqBase_.name_)];
		readNamesPositions[getPacbioGroupNumber(readEnum.element.seqBase_.name_)].emplace_back(readEnum.index);
	}
	std::unordered_map<uint32_t, uint32_t> countsDistribution;
	std::unordered_map<uint32_t, std::vector<uint64_t>> countsDistributionLens;
	for(const auto & count : nameCounts){
		++countsDistribution[count.second];
		for(const auto & pos : readNamesPositions[count.first]){
			countsDistributionLens[count.second].emplace_back(inReads[pos].seqBase_.seq_.length());
		}
	}
	table outTab{VecStr{"count","frequency", "readNumber", "meanLen", "medianLen", "minLen", "maxLen"}};
	for(const auto & cd : countsDistribution){
		auto stats = getStatsOnVec(countsDistributionLens[cd.first]);
		outTab.content_.emplace_back(toVecStr(cd.first, cd.second,
				getPercentageString(cd.first * cd.second, inReads.size()), stats["mean"],
				stats["median"], stats["min"], stats["max"]));
	}

	outTab.sortTable("count", false);
	outTab.outPutContentOrganized(std::cout);
	outTab.outPutContents(TableIOOpts(OutOptions(setUp.pars_.directoryName_ + "info", ".tab.txt"), "\t", outTab.hasHeader_));
	std::ofstream outFile;
	std::ofstream outFileTest;
	openTextFile(outFile,setUp.pars_.directoryName_ + setUp.pars_.ioOptions_.out_.outFilename_.string(), ".fastq", setUp.pars_.ioOptions_.out_);
	openTextFile(outFileTest,setUp.pars_.directoryName_ + setUp.pars_.ioOptions_.out_.outFilename_.string() + "_test", ".fastq", setUp.pars_.ioOptions_.out_);

	for(const auto positions : readNamesPositions){
		if(nameCounts[positions.first] > 2){
			auto paddedSizeName = leftPadNumStr(nameCounts[positions.first], estd::stou(outTab.content_.back()[0]));
			auto dirName = setUp.pars_.directoryName_ + paddedSizeName + "/";
			if(!njh::files::bfs::exists(dirName)){
				dirName = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar(paddedSizeName , false)).string();
			}
			inReads[positions.second.front()].seqBase_.outPutFastq(outFileTest);
			std::ofstream currentOutFile;
			openTextFile(currentOutFile,dirName + estd::to_string(positions.first), ".fastq", setUp.pars_.ioOptions_.out_);
			for(const auto & pos : positions.second){
				inReads[pos].seqBase_.outPutFastq(currentOutFile);
			}
		}
	}
	return 0;
}

int pacbioExpRunner::testCssConsensusBuilding(const njh::progutils::CmdArgs & inputCommands){
	seqSetUp setUp(inputCommands);
	setUp.pars_.ioOptions_.out_.outFilename_ = "out";
	setUp.processDefaultReader(true);
	uint32_t numThreads = 1;
	setUp.setOption(numThreads, "-t,--numThreads", "Number of Threads to Use");
	setUp.processAlignerDefualts();
	setUp.setOption(setUp.pars_.colOpts_.kmerOpts_.kLength_, "--kmerLength", "Length for kmers");
	setUp.processRefFilename(true);
	setUp.processVerbose();
	setUp.finishSetUp(std::cout);

	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<readObject>();


  std::vector<std::unique_ptr<seqWithKmerInfo>> reads;
  for(const auto & read : inReads){
  	reads.emplace_back(std::make_unique<seqWithKmerInfo>(read.seqBase_));
  }

  uint64_t maxLength = 0;
  readVec::getMaxLength(inReads, maxLength);
  auto refSeqs = SeqInput::getReferenceSeq(setUp.pars_.refIoOptions_, maxLength);
  std::vector<std::unique_ptr<seqWithKmerInfo>> refReads;
	for(const auto & read : refSeqs){
		refReads.emplace_back(std::make_unique<seqWithKmerInfo>(read.seqBase_));
	}

  allSetKmers(reads, setUp.pars_.colOpts_.kmerOpts_.kLength_,false);
  allSetKmers(refReads, setUp.pars_.colOpts_.kmerOpts_.kLength_,false);
  std::function<double(const std::unique_ptr<seqWithKmerInfo> & read1,
  		const std::unique_ptr<seqWithKmerInfo> & read2)> disFun =
  		[](const std::unique_ptr<seqWithKmerInfo> & read1,
  				const std::unique_ptr<seqWithKmerInfo> & read2){
  	auto dist = read1->compareKmers(*read2);
  	return dist.second;
  };

  auto distances = getDistance(reads, numThreads, disFun);

  njhUndirWeightedGraph<double, uint32_t> distGraph;
  for(const auto & pos : iter::range(reads.size())){
  	distGraph.addNode(reads[pos]->seqBase_.name_, pos);
  }
  for(const auto & pos : iter::range(distances.size())){
  	for(const auto & subPos : iter::range<uint64_t>(distances[pos].size())){
  		distGraph.addEdge(inReads[pos].seqBase_.name_,
  				inReads[subPos].seqBase_.name_, distances[pos][subPos]);
  	}
  }

  aligner alignerObj(maxLength, setUp.pars_.gapInfo_, setUp.pars_.scoring_);
  alignerObj.processAlnInfoInput(setUp.pars_.alnInfoDirName_);
  std::ofstream outFile;
  openTextFile(outFile, setUp.pars_.ioOptions_.out_.outFilename_, ".tab.txt", setUp.pars_.ioOptions_.out_);

  int counter = 0;
  bool eventBased = true;
  outFile
      << "ReadNumber\tReadId\tReadFraction\tBestRef\tavgDist\tscore\t1bIndel\t2bI"
         "ndel\t>2bIndel\tlqMismatch\thqMismatch" << std::endl;
  for (const auto& inputPos : iter::range(reads.size())) {
  	const auto & input = reads[inputPos];
    if ((counter + 1) % 5 == 0 && setUp.pars_.verbose_) {
      std::cout << "Currently on read " << counter + 1 << " out of "
                << reads.size() << std::endl;
    }
    double bestScore = std::numeric_limits<double>::lowest();
    std::vector<uint64_t> bestRead;
    for (const auto& refPos : iter::range(refReads.size())) {
    	const auto & ref = refReads[refPos];
      if (input->seqBase_.name_ == ref->seqBase_.name_) {
        continue;
      }
      alignerObj.alignCache(*ref, *input, false);
      double currentScore = 0;

      if(eventBased){
        alignerObj.profileAlignment(*ref, *input, false, true, false);
      	currentScore = alignerObj.comp_.distances_.eventBasedIdentity_;
      }else{
      	currentScore = alignerObj.parts_.score_;
      }

      if (currentScore == bestScore) {
        bestRead.push_back(refPos);
      }
      if (currentScore > bestScore) {
        bestRead.clear();
        bestRead.push_back(refPos);
        bestScore = currentScore;
      }
    }
    for (const auto& bestPos : bestRead) {
    	const auto & best = refReads[bestPos];
      alignerObj.alignCache(*best, *input, false);
      double score = 0;
      alignerObj.profileAlignment(*best, *input, false, true, false );
      if(eventBased){
      	score = alignerObj.comp_.distances_.eventBasedIdentity_;
      } else {
      	score = alignerObj.parts_.score_;
      }


      outFile << counter << "\t" << input->seqBase_.name_ << "\t"
                      << input->seqBase_.frac_ << "\t" << best->seqBase_.name_
											<< "\t" << distGraph.nodes_[inputPos]->getAverageDist()
                      << "\t" << score 	<< "\t"
                      << alignerObj.comp_.oneBaseIndel_ << "\t"
                      << alignerObj.comp_.twoBaseIndel_ << "\t"
                      << alignerObj.comp_.largeBaseIndel_ << "\t"
                      << alignerObj.comp_.lqMismatches_ << "\t"
                      << alignerObj.comp_.hqMismatches_ << std::endl;
    }
    counter++;
  }

	return 0;
}


int pacbioExpRunner::profilePacbioReads(const njh::progutils::CmdArgs & inputCommands){
  seqSetUp setUp(inputCommands);
  uint32_t windowSize = 50;
  uint32_t windowStep = 10;
	uint32_t kStart = 2;
	uint32_t kStop = 10;
	uint32_t readThres = 2;
	double snpFracThres = 0.001;
	uint32_t snpHardCut = 1;
	setUp.setOption(readThres, "--readThres", "Number of reads in order to be considered a new cluster");
	setUp.setOption(snpFracThres, "--snpFracThres", "Number of times a SNP has to appear in order to be considered high occurring, .001 corresponds to 0.1%");
	setUp.setOption(snpHardCut, "--snpHardCut", "Never consider SNP if it is this or less");
	setUp.setOption(kStart, "--kStart", "Starting kmer Analysis at this k len");
	setUp.setOption(kStop, "--kStop", "Stoping kmer Analysis at this k len");
  setUp.setOption(windowSize, "--windowSize", "Size of scanning kmer window");
  setUp.setOption(windowStep, "--windowStep", "size of the scanning kmer dist step");
  setUp.processDefaultReader(true);
  setUp.processSeq(true);
  setUp.pars_.gapLeft_= "0,0";
  setUp.pars_.gapRight_= "0,0";
  setUp.pars_.gapInfo_ = gapScoringParameters(setUp.pars_.gap_, setUp.pars_.gapLeft_, setUp.pars_.gapRight_);
  setUp.pars_.colOpts_.kmerOpts_.kLength_ = 10;
  setUp.processAlignerDefualts();
//profileSharedKmerBlocks
  setUp.processVerbose();
  setUp.processDirectoryOutputName(true);
  setUp.finishSetUp(std::cout);
  setUp.startARunLog(setUp.pars_.directoryName_);
  SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<readObject>();


  std::vector<std::unique_ptr<seqWithKmerInfo>> reads;
  for(const auto & read : inReads){
  	reads.emplace_back(std::make_unique<seqWithKmerInfo>(read.seqBase_));
  }
	uint64_t maxLength = 0;
	readVec::getMaxLength(inReads, maxLength);
	std::unique_ptr<seqWithKmerInfo> refRead = std::make_unique<seqWithKmerInfo>(
			setUp.pars_.seqObj_.seqBase_);
	readVec::getMaxLength(*refRead, maxLength);
	aligner alignerObj(maxLength, setUp.pars_.gapInfo_, setUp.pars_.scoring_,
			KmerMaps(), setUp.pars_.qScorePars_,
			setUp.pars_.colOpts_.alignOpts_.countEndGaps_,
			setUp.pars_.colOpts_.iTOpts_.weighHomopolyer_);
	alignerObj.processAlnInfoInput(setUp.pars_.alnInfoDirName_);
	std::ofstream outProfileFile;
  openTextFile(outProfileFile, setUp.pars_.directoryName_ + "info", ".tab.txt", false, false);
  allSetKmers(reads, setUp.pars_.colOpts_.kmerOpts_.kLength_, false);
  refRead->setKmers(setUp.pars_.colOpts_.kmerOpts_.kLength_, false);
	outProfileFile
			<< "ref\trefLen\tseq\tseqLen\tstart\tstop\tnumInserts\tbasesInInserts\tfracInsert\tnumDeletions"
					"\tbaseInDele\tfracDeletions\tnumMismatches\tfracSub\teventBasedIdentity\tklen\tkShared\tkDist\n";
	std::unordered_map<std::string, double> perIds;

	std::unordered_map<uint32_t, std::unordered_map<char, VecStr>> mismatches;

	for (const auto & readPos : iter::range(reads.size())) {
		if(readPos %100 == 0){
			std::cout << "On " << readPos << " of " << reads.size() << "\r";
			std::cout.flush();
		}
		const auto & read = reads[readPos];
		alignerObj.alignCache(*refRead, *read, false);
		//count gaps and mismatches and get identity
		alignerObj.profilePrimerAlignment(*refRead, *read);
		for(const auto & m : alignerObj.comp_.distances_.mismatches_){
			mismatches[m.second.refBasePos][m.second.seqBase].emplace_back(read->seqBase_.name_);
		}
		//calc indels
		uint32_t numOfInserts = 0;
		uint32_t basesInInserts = 0;
		uint32_t numOfDeletions = 0;
		uint32_t basesInDeletions = 0;
		for(const auto & g : alignerObj.comp_.distances_.alignmentGaps_){
			if(g.second.ref_){
				++numOfInserts;
				basesInInserts += g.second.size_;
			}else{
				++numOfDeletions;
				basesInDeletions += g.second.size_;
			}
		}
		//get kmer dist
		auto kDist = refRead->compareKmers(*read);
		//get start and end by counting ends gaps
		auto frontGapSize = countBeginChar(alignerObj.alignObjectB_.seqBase_.seq_);
		if(alignerObj.alignObjectB_.seqBase_.seq_.front() != '-'){
			frontGapSize = 0;
		}
		auto backGapSize = countEndChar(alignerObj.alignObjectB_.seqBase_.seq_);
		if(alignerObj.alignObjectB_.seqBase_.seq_.back() != '-'){
			backGapSize = 0;
		}
		perIds[read->seqBase_.name_] = alignerObj.comp_.distances_.eventBasedIdentity_;
		outProfileFile << refRead->seqBase_.name_
				<< "\t" << refRead->seqBase_.seq_.length()
				<< "\t" << read->seqBase_.name_
				<< "\t" << read->seqBase_.seq_.length()
				<< "\t" << frontGapSize
				<< "\t" << refRead->seqBase_.seq_.length() - backGapSize
				<< "\t" << numOfInserts
				<< "\t" << basesInInserts
				<< "\t" << basesInInserts/static_cast<double>(read->seqBase_.seq_.length())
				<< "\t" << numOfDeletions
				<< "\t" << basesInDeletions
				<< "\t" << basesInDeletions/static_cast<double>(refRead->seqBase_.seq_.length() - backGapSize - frontGapSize)
				<< "\t" << alignerObj.comp_.hqMismatches_
				<< "\t" << alignerObj.comp_.hqMismatches_/static_cast<double>(read->seqBase_.seq_.length() - basesInInserts)
				<< "\t" << alignerObj.comp_.distances_.eventBasedIdentity_
				<< "\t" << setUp.pars_.colOpts_.kmerOpts_.kLength_
				<< "\t" << kDist.first
				<< "\t" << kDist.second << "\n";
	}
	std::cout << std::endl;
	table misTab{VecStr{"refPos", "refBase", "seqBase","freq", "fraction", "seqs"}};
	double readTotal = readVec::getTotalReadCount(inReads);
	for (const auto & m : mismatches) {
		for (const auto & seqM : m.second) {
			misTab.content_.emplace_back(
					toVecStr(m.first, setUp.pars_.seqObj_.seqBase_.seq_[m.first], seqM.first,
							seqM.second.size(), seqM.second.size()/readTotal,
							vectorToString(seqM.second, ",")));
		}
	}
	misTab.sortTable("seqBase", false);
	misTab.sortTable("refPos", false);
	misTab.outPutContents(TableIOOpts(OutOptions(setUp.pars_.directoryName_ + "mismatches", ".tab.txt"), "\t", misTab.hasHeader_));
	std::unordered_map<std::string, readsWithSnps> erroneousReads;
	for (const auto & readPos : iter::range(reads.size())) {
		const auto & read = reads[readPos];
		alignerObj.alignCache(*refRead, *read, false);
		//count gaps and mismatches and get identity
		alignerObj.profilePrimerAlignment(*refRead, *read);
		std::map<uint32_t, mismatch> highMismatches;
		for(const auto & m : alignerObj.comp_.distances_.mismatches_){
			double mFrac = mismatches[m.second.refBasePos][m.second.seqBase].size()/readTotal;
			if(mFrac >= snpFracThres && mismatches[m.second.refBasePos][m.second.seqBase].size() > snpHardCut){
				highMismatches.emplace(m.second.refBasePos, m.second);
				highMismatches.at(m.second.refBasePos).freq = mismatches[m.second.refBasePos][m.second.seqBase].size();
			}
		}
		if(highMismatches.size() > 0){
			std::string mCode = "";
			for(const auto & hm : highMismatches){
				mCode.append(estd::to_string(hm.first) + estd::to_string(hm.second.seqBase));
			}
			auto search = erroneousReads.find(mCode);
			if(search == erroneousReads.end()){
				erroneousReads.emplace(mCode, readsWithSnps{highMismatches, readPos});
			}else{
				search->second.addReadPos(readPos);
			}
		}
	}
	if (erroneousReads.size() > 0) {
		bfs::path extracClustersDir = njh::files::makeDir(setUp.pars_.directoryName_,
				njh::files::MkdirPar("extraClusters", false));
		table outInfo { VecStr { "code", "freq", "frac" } };
		for (const auto & subReads : erroneousReads) {
			if(len(subReads.second) < readThres){
				continue;
			}
			outInfo.content_.emplace_back(
					toVecStr(subReads.first, len(subReads.second),
							len(subReads.second) / readTotal));
			std::ofstream outReadFile;
			openTextFile(outReadFile, njh::files::make_path(extracClustersDir, subReads.first).string(),
					".fastq", false, false);
			for (const auto & pos : subReads.second.readPositions_) {
				inReads[pos].seqBase_.outPutFastq(outReadFile);
			}
		}
		outInfo.sortTable("freq", false);
		outInfo.outPutContents(TableIOOpts{OutOptions(njh::files::make_path(extracClustersDir, "info").string(), ".tab.txt"), "\t", outInfo.hasHeader_});
	}

	std::ofstream kDistFile;
	openTextFile(kDistFile, setUp.pars_.directoryName_ + "kDistances.tab.txt", ".tab.txt", false, false);

	kDistFile << "ref\trefLen\tseq\tseqLen\teventBasedIdentity\tklen\tkShared\tkDist\n";
	for(const auto & k : iter::range(kStart,kStop + 1)){
		std::cout << "On k length of " << k << " of " << kStop << "\r";
		std::cout.flush();
		refRead->setKmers(k, false);
		allSetKmers(reads, k, false);
		for(const auto & read : reads){
			auto kDist = refRead->compareKmers(*read);
			kDistFile << refRead->seqBase_.name_
					<< "\t" << refRead->seqBase_.seq_.length()
					<< "\t" << read->seqBase_.name_
					<< "\t" << read->seqBase_.seq_.length()
					<< "\t" << perIds[read->seqBase_.name_]
					<< "\t" << k
					<< "\t" << kDist.first
					<< "\t" << kDist.second << "\n";
		}
	}
	std::cout << std::endl;
	if(setUp.pars_.writingOutAlnInfo_){
  	alignerObj.alnHolder_.write(setUp.pars_.alnInfoDirName_);
  }
  std::cout << alignerObj.numberOfAlingmentsDone_ << std::endl;
  setUp.logRunTime(std::cout);
	return 0;
}


} /* namespace njhseq */
