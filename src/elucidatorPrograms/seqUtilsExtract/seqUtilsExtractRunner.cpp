
//  seqUtilsExtractRunner.cpp
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
    
#include "seqUtilsExtractRunner.hpp"
    
    
namespace njhseq {

seqUtilsExtractRunner::seqUtilsExtractRunner()
    : njh::progutils::ProgramRunner({
	addFunc("extractBySeq", extractBySeq, false),
  addFunc("extractSameSeqs", extractSameSeqs, false),
  addFunc("extractByMIDs", extractByMIDs, false),
  addFunc("binOnNucComp", binOnNucComp, false),
  addFunc("binOnNucCompFaster", binOnNucCompFaster, false),
	addFunc("greedyKmerCluster", greedyKmerCluster, false),
  addFunc("extractSeqsBeginsWith", extractSeqsBeginsWith, false),
	addFunc("extractSeqsEndsWith", extractSeqsEndsWith, false),
	addFunc("extractSeqsBeginsWithEndsWith", extractSeqsBeginsWithEndsWith, false),
	addFunc("extractByIlluminaAaptors", extractByIlluminaAaptors, false),
	addFunc("extractByName", extractByName, false),
},
                    "seqUtilsExtract") {}
//
int seqUtilsExtractRunner::extractByName(const njh::progutils::CmdArgs & inputCommands){
	seqSetUp setUp(inputCommands);
	bfs::path filename = "";
	bool excluding = false;
	setUp.setOption(excluding, "--excluding", "Excluding These Names Instead, default is extract by these names");
	setUp.setOption(filename, "--file", "File Containing Seq Names in the first column column" , true);
	setUp.pars_.ioOptions_.out_.outFilename_ = "out";
	setUp.processDefaultReader(true);
	setUp.finishSetUp(std::cout);
	VecStr readNames;
	if(!bfs::exists(filename)){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << "Error, file " << filename << " doesn't exist" << "\n";
		throw std::runtime_error{ss.str()};
	}
	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	reader.openOut();
	std::string line = "";
	std::ifstream infile(filename.string());
	if(!infile){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << "Error, in opening " << filename << "" << "\n";
		throw std::runtime_error{ss.str()};
	}
	while(njh::files::crossPlatGetline(infile, line)){
		readNames.emplace_back(line);
	}
	readObject read;
	while(reader.readNextRead(read)){
		if(excluding){
			if(!njh::in(read.seqBase_.name_, readNames)){
				reader.write(read);
			}
		}else{
			if(njh::in(read.seqBase_.name_, readNames)){
				reader.write(read);
			}
		}
	}
	return 0;
}
int seqUtilsExtractRunner::extractSeqsBeginsWith(const njh::progutils::CmdArgs & inputCommands) {
	bool allowableErrors = 0;
	std::string beginsWith = "";
	bool writeOtherSeqs = false;
	seqUtilsExtractSetUp setUp(inputCommands);

	setUp.processVerbose();
	setUp.processDebug();
	setUp.processDefaultReader(true);
	setUp.setOption(beginsWith, "--beginsWith", "A motif for the sequences to begin with", true);
	setUp.setOption(writeOtherSeqs, "--writeOtherSeqs", "Write out other sequences as well");
	setUp.setOption(allowableErrors, "--allowableErrors", "Allowable Errors");

	setUp.finishSetUp(std::cout);

	motif beginsWithMotif(beginsWith);
	SeqIO readerWriter(setUp.pars_.ioOptions_);
	readerWriter.openIn();
	readerWriter.openOut();
	SeqIOOptions otherSeqsOpts = setUp.pars_.ioOptions_;
	otherSeqsOpts.out_.outFilename_ = njh::files::prependFileBasename(otherSeqsOpts.out_.outFilename_, "otherSeqs_");
	SeqOutput otherWriter(otherSeqsOpts);
	if(writeOtherSeqs){
		otherWriter.openOut();
	}
	seqInfo seq;

	while(readerWriter.readNextRead(seq)){
		if(len(seq) >= beginsWithMotif.size() &&
				beginsWithMotif.scoreMotif(seq.seq_.begin(), seq.seq_.begin() + beginsWithMotif.size()) + allowableErrors >= beginsWithMotif.size()){
			readerWriter.write(seq);
		}else if(writeOtherSeqs){
			otherWriter.write(seq);
		}
	}
	return 0;
}

int seqUtilsExtractRunner::extractSeqsEndsWith(const njh::progutils::CmdArgs & inputCommands) {
	bool allowableErrors = 0;
	std::string endsWith = "";
	bool writeOtherSeqs = false;
	seqUtilsExtractSetUp setUp(inputCommands);

	setUp.processVerbose();
	setUp.processDebug();
	setUp.processDefaultReader(true);
	setUp.setOption(endsWith, "--endsWith", "A motif for the sequences to ends with", true);
	setUp.setOption(writeOtherSeqs, "--writeOtherSeqs", "Write out other sequences as well");
	setUp.setOption(allowableErrors, "--allowableErrors", "Allowable Errors");

	setUp.finishSetUp(std::cout);

	motif endsWithMotif(endsWith);
	SeqIO readerWriter(setUp.pars_.ioOptions_);
	readerWriter.openIn();
	readerWriter.openOut();
	SeqIOOptions otherSeqsOpts = setUp.pars_.ioOptions_;
	otherSeqsOpts.out_.outFilename_ = njh::files::prependFileBasename(otherSeqsOpts.out_.outFilename_, "otherSeqs_");
	SeqOutput otherWriter(otherSeqsOpts);
	if(writeOtherSeqs){
		otherWriter.openOut();
	}
	seqInfo seq;

	while(readerWriter.readNextRead(seq)){
		if(len(seq) >= endsWithMotif.size() &&
				endsWithMotif.scoreMotif(seq.seq_.end() - endsWithMotif.size(), seq.seq_.end()) + allowableErrors >= endsWithMotif.size()){
			readerWriter.write(seq);
		}else if(writeOtherSeqs){
			otherWriter.write(seq);
		}
	}
	return 0;
}

int seqUtilsExtractRunner::extractSeqsBeginsWithEndsWith(const njh::progutils::CmdArgs & inputCommands) {
	bool allowableErrors = 0;
	std::string beginsWith = "";
	std::string endsWith = "";

	bool writeOtherSeqs = false;
	seqUtilsExtractSetUp setUp(inputCommands);

	setUp.processVerbose();
	setUp.processDebug();
	setUp.processDefaultReader(true);
	setUp.setOption(beginsWith, "--beginsWith", "A motif for the sequences to begin with", true);
	setUp.setOption(writeOtherSeqs, "--writeOtherSeqs", "Write out other sequences as well");
	setUp.setOption(endsWith, "--endsWith", "A motif for the sequences to ends with", true);
	setUp.setOption(allowableErrors, "--allowableErrors", "Allowable Errors");

	setUp.finishSetUp(std::cout);

	motif beginsWithMotif(beginsWith);
	motif endsWithMotif(endsWith);

	SeqIO readerWriter(setUp.pars_.ioOptions_);
	readerWriter.openIn();
	readerWriter.openOut();
	SeqIOOptions otherSeqsOpts = setUp.pars_.ioOptions_;
	otherSeqsOpts.out_.outFilename_ = njh::files::prependFileBasename(otherSeqsOpts.out_.outFilename_, "otherSeqs_");
	SeqOutput otherWriter(otherSeqsOpts);
	if(writeOtherSeqs){
		otherWriter.openOut();
	}

	if(setUp.pars_.ioOptions_.isPairedIn()){
		PairedRead seq;
		while(readerWriter.readNextRead(seq)){
			if(len(seq) >= beginsWithMotif.size() &&
					beginsWithMotif.scoreMotif(seq.seqBase_.seq_.begin(), seq.seqBase_.seq_.begin() + beginsWithMotif.size() ) + allowableErrors >= beginsWithMotif.size()&&
					endsWithMotif.scoreMotif(seq.mateSeqBase_.seq_.begin(), seq.mateSeqBase_.seq_.begin() + beginsWithMotif.size() ) + allowableErrors >= endsWithMotif.size()){
				readerWriter.write(seq);
			}else if(writeOtherSeqs){
				otherWriter.write(seq);
			}
		}
	}else{
		seqInfo seq;
		while(readerWriter.readNextRead(seq)){
			if(len(seq) >= beginsWithMotif.size() &&
					beginsWithMotif.scoreMotif(seq.seq_.begin(), seq.seq_.begin() + beginsWithMotif.size()) + allowableErrors >= beginsWithMotif.size()&&
					endsWithMotif.scoreMotif(seq.seq_.end() - endsWithMotif.size(), seq.seq_.end()) + allowableErrors >= endsWithMotif.size()){
				readerWriter.write(seq);
			}else if(writeOtherSeqs){
				otherWriter.write(seq);
			}
		}
	}

	return 0;
}


int seqUtilsExtractRunner::greedyKmerCluster(const njh::progutils::CmdArgs & inputCommands) {
	seqUtilsExtractSetUp setUp(inputCommands);
	bool checkComplement = false;
	double cutOff = 0.20;
	bool sort = false;
	setUp.setOption(sort, "--sort", "Sort the sequences first");

	setUp.setOption(cutOff, "--cutOff", "Kmer Similarity Score Cut off");
	setUp.setOption(checkComplement, "--checkComplement", "checkComplement");
	setUp.processKmerLenOptions();
	setUp.processDefaultReader(true);
	setUp.processDirectoryOutputName(true);
	setUp.processVerbose();
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);

	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<readObject>();
	if(sort){
		readVecSorter::sortReadVector(inReads, "seq", false);
	}
	std::vector<kmerCluster> kClusters = greedyKmerSimCluster(inReads,
			setUp.pars_.colOpts_.kmerOpts_.kLength_, cutOff, checkComplement,
			setUp.pars_.verbose_);
	std::ofstream infoFile;
	openTextFile(infoFile, setUp.pars_.directoryName_ + "info.tab.txt", ".txt",
			setUp.pars_.ioOptions_.out_);
	infoFile << "kmerCluster\treadNumber" << std::endl;
	std::sort(kClusters.begin(), kClusters.end(),
			[&](const kmerCluster & kClus1, const kmerCluster & kClus2) {
				return kClus1.reads_.size() > kClus2.reads_.size();});

	for (const auto & pos : iter::range(kClusters.size())) {
		auto outOpts = SeqIOOptions::genFastqOut(setUp.pars_.directoryName_
						+ estd::to_string(leftPadNumStr(pos, kClusters.size())));
		kClusters[pos].writeInfo(outOpts);
		infoFile << pos << "\t"
				<< getPercentageString(kClusters[pos].reads_.size(), inReads.size())
				<< std::endl;
	}
	setUp.logRunTime(std::cout);
	return 0;
}

int seqUtilsExtractRunner::binOnNucComp(const njh::progutils::CmdArgs & inputCommands) {
  seqSetUp setUp(inputCommands);
  double diffCutOff = 0.1;
  bool findBest = false;
  bool preSort = true;
  std::string alphabetStr = "A,C,G,T";
  uint32_t readBuffer = 10;
  bool write = false;
  uint32_t smallClusterSize = 1;
  setUp.setOption(smallClusterSize, "-smallClusterSize", "Small Cluster Size (inclusive)"	);
  setUp.setOption(write, "-write", "write");
  setUp.setOption(preSort, "-preSort", "preSort");
  setUp.setOption(readBuffer, "-readBuffer", "readBuffer");
  setUp.processAlignerDefualts();
  setUp.processDefaultReader(true);
  setUp.setOption(diffCutOff, "-diffCutOff","diffCutOff");
  setUp.setOption(findBest, "-findBest","findBest");
  setUp.processVerbose();
  if(write){
  	setUp.processDirectoryOutputName(true);
  }
  setUp.finishSetUp(std::cout);
  std::vector<char> alphabet = processAlphStrVecChar(alphabetStr, ",");
  std::vector<nucCompCluster> nucComps;
	SeqInput reader(setUp.pars_.ioOptions_);
		reader.openIn();
	auto inReads = reader.readAllReads<readObject>();
  std::vector<nucCompCluster> comps = clusterOnNucComp(inReads,
  		alphabet, diffCutOff, findBest, preSort, setUp.pars_.verbose_);
  table outInfo = getInfoNucCompVec(comps);
  outInfo.outPutContentOrganized(std::cout);

  if(write){
  	std::string smallDir = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar("smallClusters")).string();
  	TableIOOpts outInfoOpts(OutOptions(setUp.pars_.directoryName_ + "compInfo", ".tab.txt"), "\t", outInfo.hasHeader_);
  	outInfo.outPutContents(outInfoOpts)	;
  	for(const auto & compPos : iter::range(comps.size())){
  		std::vector<readObject> currentClusterReads = comps[compPos].getReads(inReads);
  		if(currentClusterReads.size() <= smallClusterSize){
    		SeqOutput::write(currentClusterReads, setUp.pars_.directoryName_ + leftPadNumStr(compPos, comps.size()),
    				setUp.pars_.ioOptions_);
  		}else{
    		SeqOutput::write(currentClusterReads, smallDir + leftPadNumStr(compPos, comps.size()),
    				setUp.pars_.ioOptions_);
  		}
  	}
  }
  setUp.logRunTime(std::cout);
	return 0;
}
int seqUtilsExtractRunner::binOnNucCompFaster(const njh::progutils::CmdArgs & inputCommands) {
  seqSetUp setUp(inputCommands);
  double diffCutOff = 0.1;
  bool findBest = false;
  std::string alphabetStr = "A,C,G,T";
  uint32_t readBuffer = 10;
  bool preSort = true;
  bool write = false;
  uint32_t smallClusterSize = 1;
  setUp.setOption(smallClusterSize, "-smallClusterSize", "Small Cluster Size (inclusive)"	);
  setUp.setOption(write, "-write", "write");
  setUp.setOption(readBuffer, "-readBuffer", "readBuffer");
  setUp.processAlignerDefualts();
  setUp.processDefaultReader(true);
  setUp.setOption(diffCutOff, "-diffCutOff","diffCutOff");
  setUp.setOption(findBest, "-findBest","findBest");
  setUp.processVerbose();
  setUp.processDirectoryOutputName(write);
  //setUp.processDirectoryOutputName(true);
  setUp.finishSetUp(std::cout);
  if(write){
  	 setUp.startARunLog(setUp.pars_.directoryName_);
  }
  nucCompCluster::readBufferMax_ = readBuffer;
  std::vector<char> alphabet = processAlphStrVecChar(alphabetStr, ",");
  std::vector<nucCompCluster> comps = clusterOnNucComp(setUp.pars_.ioOptions_,
  		alphabet, diffCutOff, findBest, preSort, setUp.pars_.verbose_);
  table outInfo = getInfoNucCompVec(comps);
  outInfo.outPutContentOrganized(std::cout);
  if(write){

  	std::string smallDir = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar("smallClusters")).string();
  	outInfo.outPutContents(TableIOOpts(OutOptions(setUp.pars_.directoryName_ + "compClusInfo.tab.txt",".txt"), "\t", outInfo.hasHeader_));
  	for(const auto & compPos : iter::range(comps.size())){
  		std::vector<readObject> currentClusterReads = comps[compPos].getReads<readObject>(setUp.pars_.ioOptions_);
  		if(currentClusterReads.size() <= smallClusterSize){
    		SeqOutput::write(currentClusterReads, smallDir + leftPadNumStr(compPos, comps.size()),
    				setUp.pars_.ioOptions_);
  		}else{
    		SeqOutput::write(currentClusterReads, setUp.pars_.directoryName_ + leftPadNumStr(compPos, comps.size()),
    				setUp.pars_.ioOptions_);
  		}
  	}
  }
  setUp.logRunTime(std::cout);
  return 0;
}
//
//int seqUtilsExtractRunner::extractByMID(const njh::progutils::CmdArgs & inputCommands) {
//  // parameters
//	ExtractByMIDPars pars;
//
//	seqUtilsExtractSetUp setUp(inputCommands);
//  setUp.setUpExtractByMID(pars);
//  setUp.startARunLog(setUp.pars_.directoryName_);
//  setUp.writeParametersFile(setUp.pars_.directoryName_ + "parametersUsed.tab.txt",
//                            true, false);
//  int32_t midSize;
//	table mids = seqUtil::readBarcodes(pars.idFilename.string(), pars.idFileDelim,
//			midSize);
//	std::unique_ptr<MidDeterminator> determinator = std::make_unique<MidDeterminator>(mids);
//	readObject read;
//  SeqInput reader(setUp.pars_.ioOptions_);
//  reader.openIn();
//
//	MultiSeqIO readerOuts;
//	for (const auto & mid : determinator->mids_) {
//		auto midOpts = setUp.pars_.ioOptions_;
//		midOpts.out_.outFilename_ = njh::files::make_path(setUp.pars_.directoryName_, mid.first).string();
//		if (setUp.pars_.debug_) {
//			std::cout << "Inserting: " << mid.first << std::endl;
//		}
//		readerOuts.addReader(mid.first, midOpts);
//	}
//	auto smallFragmentOpts = setUp.pars_.ioOptions_;
//	smallFragmentOpts.out_.outFilename_ = njh::files::make_path(setUp.pars_.directoryName_, "smallFragment").string();
//	readerOuts.addReader("smallFragment", smallFragmentOpts);
//	auto failureCases = MidDeterminator::midPos::getFailureCaseNames();
//	for(const auto & failureCase : failureCases){
//		std::string unRecName = "unrecognizedBarcode_" + failureCase;
//		auto midOpts = setUp.pars_.ioOptions_;
//		midOpts.out_.outFilename_ = setUp.pars_.directoryName_ + unRecName;
//		if (setUp.pars_.debug_){
//			std::cout << "Inserting: " << unRecName << std::endl;
//		}
//		readerOuts.addReader(unRecName, midOpts);
//	}
//
//	if(setUp.pars_.verbose_){
//		std::cout << njh::bashCT::boldGreen("Extracting on MIDs") << std::endl;
//	}
//	uint32_t count = 0;
//	uint32_t smallFragmentCount = 0;
//	uint32_t readsNotMatchedToBarcode = 0;
//	std::unordered_map<std::string, std::pair<uint32_t, uint32_t>> counts;
//	std::unordered_map<std::string, uint32_t>  failBarCodeCounts;
//	uint32_t found = 0;
//
//	while (reader.readNextRead(read)) {
//		++count;
//		if (setUp.pars_.verbose_ && count % 50 == 0) {
//			std::cout << "\r" << count ;
//			std::cout.flush();
//		}
//		readVec::handelLowerCaseBases(read, setUp.pars_.ioOptions_.lowerCaseBases_);
//
//		//possibly trim reads at low quality
//		if (len(read) < pars.smallFragmentCutoff) {
//			readerOuts.openWrite("smallFragment", read);
//			++smallFragmentCount;
//			continue;
//		}
//
//
//		std::pair<MidDeterminator::midPos, MidDeterminator::midPos> currentMid = determinator->fullDetermine(read, pars.mDetPars);
//		if ("unrecognized" != currentMid.first.midName_) {
//			++found;
//			if (read.seqBase_.name_.find("_Comp") != std::string::npos) {
//				++counts[currentMid.first.midName_].second;
//			} else {
//				++counts[currentMid.first.midName_].first;
//			}
//			/**@todo need to reorient the reads here before outputing if that's needed*/
//			readerOuts.openWrite(currentMid.first.midName_, read);
//		} else {
//			++readsNotMatchedToBarcode;
//			MidDeterminator::increaseFailedBarcodeCounts(currentMid.first, failBarCodeCounts);
//			std::string unRecName = "unrecognizedBarcode_" + MidDeterminator::midPos::getFailureCaseName(currentMid.first.fCase_);
//			readerOuts.openWrite(unRecName, read);
//		}
//
//	}
//	if(setUp.pars_.verbose_){
//		std::cout << std::endl;
//	}
//	//close mid outs;
//	readerOuts.closeOutAll();
//
//	table extractionPerMid(VecStr { "MID", "total", "forward", "reverse" });
//	table extractionOverall(VecStr { "totalInput", "extracted", "smallFragment",
//			"readsNotMatching" });
//	extractionOverall.addRow(count, getPercentageString(found, count),
//			getPercentageString(smallFragmentCount, count),
//			getPercentageString(readsNotMatchedToBarcode, count));
//	auto midCountKeys = getVectorOfMapKeys(counts);
//	njh::sort(midCountKeys);
//	for (const auto & midCountKey : midCountKeys) {
//		uint32_t midTotal = counts.at(midCountKey).first + counts.at(midCountKey).second;
//		extractionPerMid.addRow(midCountKey, midTotal,
//				getPercentageString(counts.at(midCountKey).first, midTotal),
//				getPercentageString(counts.at(midCountKey).second, midTotal));
//	}
//
//	extractionPerMid.outPutContents(TableIOOpts(OutOptions(bfs::path(setUp.pars_.directoryName_ + "extractionPerMid.tab.txt")), "\t", true));
//	extractionOverall.outPutContents(TableIOOpts(OutOptions(bfs::path(setUp.pars_.directoryName_ + "extractionOverall.tab.txt")), "\t", true));
//	std::ofstream failedBarcodeFile;
//	openTextFile(failedBarcodeFile, setUp.pars_.directoryName_ + "failedBarcode.tab.txt",
//			".txt", false, false);
//	failedBarcodeFile << "Reason\tcount"
//			<< std::endl;
//	for(const auto & count : failBarCodeCounts){
//		failedBarcodeFile << count.first << "\t" << getPercentageString(count.second, readsNotMatchedToBarcode)<< std::endl;
//	}
//  setUp.logRunTime(std::cout);
//  return 0;
//}

int seqUtilsExtractRunner::extractSameSeqs(const njh::progutils::CmdArgs & inputCommands) {
  std::string compareFileName = "";
  seqUtilsExtractSetUp setUp(inputCommands);
  setUp.setUpExtractSameSeqs(compareFileName);
  SeqInput reader(setUp.pars_.ioOptions_);
  reader.openIn();
	auto inReads = reader.readAllReads<readObject>();
  std::vector<readObject> inputReads = inReads;
  SeqIOOptions compareFileNameOpts;
  compareFileNameOpts.firstName_= compareFileName;
  compareFileNameOpts.inFormat_ = SeqIOOptions::inFormats::FASTA;
  SeqInput compareFileNameReader(setUp.pars_.ioOptions_);
  std::vector<readObject> compareReads = compareFileNameReader.readAllReads<readObject>();
  std::cout << "Read in " << inputReads.size() << " inputReads" << std::endl;
  std::cout << "Read in " << compareReads.size() << " compareReads"
            << std::endl;
  std::vector<readObject> ans =
      readVecExtractor::extractSameReadsByName(inputReads, compareReads);
  SeqOutput::write(ans, setUp.pars_.ioOptions_);
  setUp.logRunTime(std::cout);
  return 0;
}


int seqUtilsExtractRunner::extractBySeq(const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	uint32_t within = 100;
	bool checkComplement = false;
	double idCutoff = .80;
	double gapCutOff = .10;
	double coverage = .75;
	setUp.processAlignerDefualts();
	setUp.setOption(checkComplement, "-checkComplement", "checkComplement");
	setUp.setOption(within, "-within", "within");
	setUp.setOption(idCutoff, "-idCutoff", "idCutoff");
	setUp.setOption(gapCutOff, "-gapCutOff", "gapCutOff");
	setUp.setOption(coverage, "-coverage", "coverage");
	setUp.processDefaultReader(true);
	setUp.processSeq(true);
	setUp.processDirectoryOutputName(true);
	setUp.finishSetUp(std::cout);
	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<readObject>();

	substituteMatrix scoring = substituteMatrix::createDegenScoreMatrix(2, -2);
	gapScoringParameters gapPars = gapScoringParameters(5, 1);
	uint64_t maxReadLength = 0;
	readVec::getMaxLength(inReads, maxReadLength);
	aligner alignerObj(maxReadLength, gapPars, scoring);
	std::vector<readObject> readsMatched;
	std::vector<readObject> readsNotMatched;
	for (auto readPos : iter::range(inReads.size())) {
		auto & read = inReads[readPos];
		auto forPos = alignerObj.findReversePrimer(read.seqBase_.seq_, setUp.pars_.seq_);
		alignerObj.rearrangeSeq(read.seqBase_.seq_, setUp.pars_.seq_, true);
		alignerObj.profilePrimerAlignment(read, setUp.pars_.seqObj_);
		if (alignerObj.comp_.distances_.percentMatch_ > idCutoff
				&& alignerObj.comp_.distances_.query_.coverage_ > coverage
				&& alignerObj.comp_.distances_.percentGaps_ < gapCutOff
				&& forPos.first < within) {
			readsMatched.emplace_back(read);
		} else if (checkComplement) {
			read.seqBase_.reverseComplementRead(true);
			auto forPosComp = alignerObj.findReversePrimer(read.seqBase_.seq_,
					setUp.pars_.seq_);
			alignerObj.rearrangeSeq(read.seqBase_.seq_, setUp.pars_.seq_, true);
			alignerObj.profilePrimerAlignment(read, setUp.pars_.seqObj_);
			if (alignerObj.comp_.distances_.percentMatch_ > idCutoff
					&& alignerObj.comp_.distances_.query_.coverage_ > coverage
					&& alignerObj.comp_.distances_.percentGaps_ < gapCutOff
					&& forPosComp.first < within) {
				readsMatched.emplace_back(read);
			} else {
				readsNotMatched.emplace_back(read);
			}
		} else {
			readsNotMatched.emplace_back(read);
		}
	}
	SeqOutput::write(readsMatched, setUp.pars_.directoryName_ + "readsMatched",
			setUp.pars_.ioOptions_);
	SeqOutput::write(readsNotMatched, setUp.pars_.directoryName_ + "readsNotMatched",
			setUp.pars_.ioOptions_);
	return 0;
}

int seqUtilsExtractRunner::extractByMIDs(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path idFileFnp = "";
	MidDeterminator::MidDeterminePars mPars;

	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(idFileFnp, "--id", "ID file", true);
	setUp.setOption(mPars.checkComplement_, "--checkComplement", "Check Complement");
	setUp.setOption(mPars.checkForShorten_, "--checkForShorten", "Check For Shorten");
	setUp.setOption(mPars.searchStop_, "--searchStop", "Barcode must come at or before this position");
	setUp.setOption(mPars.searchStart_, "--searchStart", "Barcode search starts from this position on");
	setUp.setOption(mPars.allowableErrors_, "--allowableErrors", "Allowable Errors");
	setUp.processReadInNames(VecStr{"--fastq1", "--fastq1gz", "--fasta", "--fastq"});
	setUp.processDirectoryOutputName(true);

	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);


	MidDeterminator midDet(idFileFnp, mPars);

	std::unordered_map<std::string, uint32_t> forCounts;
	std::unordered_map<std::string, uint32_t> revCounts;

	std::unordered_map<std::string, uint32_t> failedCounts;

	if("" != setUp.pars_.ioOptions_.secondName_ ){
		//paired end
		PairedRead seq;
		MultiSeqOutCache<PairedRead> outs;
		for(const auto & mid : midDet.mids_){
			outs.addReader(mid.first, SeqIOOptions::genPairedOut(njh::files::make_path(setUp.pars_.directoryName_,mid.first)));
		}
		auto failureDir = njh::files::make_path(setUp.pars_.directoryName_, "failed");
		njh::files::makeDir(njh::files::MkdirPar{failureDir});

		for(const auto & caseName : MidDeterminator::ProcessedRes::getProcessedCaseNames()){
			outs.addReader(caseName, SeqIOOptions::genPairedOut(njh::files::make_path(failureDir,caseName)));
		}

		SeqInput reader(setUp.pars_.ioOptions_);
		reader.openIn();
		while (reader.readNextRead(seq)) {
			//std::cout << seq.seqBase_.name_ << std::endl;
			auto searchResult = midDet.searchPairedEndRead(seq);
			auto processedRes = midDet.processSearchPairedEndRead(seq, searchResult);
			if(MidDeterminator::ProcessedRes::PROCESSED_CASE::MATCH ==  processedRes.case_){
				if(processedRes.rcomplement_){
					++revCounts[processedRes.midName_];
				}else{
					++forCounts[processedRes.midName_];
				}
				outs.add(processedRes.midName_, seq);
			}else{
				outs.add(MidDeterminator::ProcessedRes::getProcessedCaseName(processedRes.case_), seq);
				++failedCounts[MidDeterminator::ProcessedRes::getProcessedCaseName(processedRes.case_)];
			}
		}
	}else{
		//paired end
		seqInfo seq;
		MultiSeqOutCache<seqInfo> outs;
		for(const auto & mid : midDet.mids_){
			outs.addReader(mid.first, SeqIOOptions(njh::files::make_path(setUp.pars_.directoryName_,mid.first), setUp.pars_.ioOptions_.outFormat_ ));
		}
		auto failureDir = njh::files::make_path(setUp.pars_.directoryName_, "failed");
		njh::files::makeDir(njh::files::MkdirPar{failureDir});

		for(const auto & caseName : MidDeterminator::ProcessedRes::getProcessedCaseNames()){
			outs.addReader(caseName, SeqIOOptions(njh::files::make_path(failureDir,caseName), setUp.pars_.ioOptions_.outFormat_ ));
		}

		SeqInput reader(setUp.pars_.ioOptions_);
		reader.openIn();
		while (reader.readNextRead(seq)) {
			//std::cout << seq.seqBase_.name_ << std::endl;
			auto searchResult = midDet.searchRead(seq);
			auto processedRes = midDet.processSearchRead(seq, searchResult);
			if(MidDeterminator::ProcessedRes::PROCESSED_CASE::MATCH ==  processedRes.case_){
				if(processedRes.rcomplement_){
					++revCounts[processedRes.midName_];
				}else{
					++forCounts[processedRes.midName_];
				}
				outs.add(processedRes.midName_, seq);
			}else{
				outs.add(MidDeterminator::ProcessedRes::getProcessedCaseName(processedRes.case_), seq);
				++failedCounts[MidDeterminator::ProcessedRes::getProcessedCaseName(processedRes.case_)];
			}
		}
	}

	std::set<std::string> midNames;
	njh::addVecToSet(getVectorOfMapKeys(midDet.mids_), midNames);
	table midCounts(VecStr{"MIDName", "ForCount", "RevCount"});
	for(const auto & mid : midNames){
		midCounts.addRow(mid, forCounts[mid], revCounts[mid]);
	}
	OutOptions midCountsOpts(njh::files::make_path(setUp.pars_.directoryName_, "midCounts.tab.txt"));
	OutputStream midCountsOut(midCountsOpts);
	midCounts.outPutContents(midCountsOut, "\t");

	table failedCountsTab(failedCounts, VecStr{"FailedCase", "Counts"});
	OutOptions failedCountsOpts(njh::files::make_path(setUp.pars_.directoryName_, "failedCounts.tab.txt"));
	OutputStream failedCountsOut(failedCountsOpts);
	failedCountsTab.outPutContents(failedCountsOut, "\t");

	return 0;
}

                    
} // namespace njhseq
