/*
 * genExpRunner.cpp
 *
 *  Created on: Jan 25, 2015
 *      Author: nickhathaway
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


#include "genExp.hpp"
#include "elucidator/simulation.h"
#include "elucidator/utils/KmerUtils.hpp"
#include "elucidator/objects/BioDataObject.h"
#include <njhseq/GenomeUtils.h>



#include <TwoBit.h>



/*
 *
 */


namespace njhseq {
genExpRunner::genExpRunner()
    : njh::progutils::ProgramRunner(
          {
					 addFunc("muscle", muscle, false),
					 addFunc("extractByName", extractByName, false),
					 addFunc("createReadConsensus", createReadConsensus, false),
					 addFunc("createReadConsensusRandomPicking", createReadConsensusRandomPicking, false),
					 addFunc("findVariableSites", findVariableSites, false),

					 addFunc("tableCountToDistGraph", tableCountToDistGraph, false),
					 addFunc("tableCountToFasta", tableCountToFasta, false),
					 addFunc("getSnpInfo", getSnpInfo, false),
					 addFunc("sortOnNucComp", sortOnNucComp, false),
					 addFunc("extractOnLen", extractOnLen, false),
					 addFunc("splitUpFile", splitUpFile, false),
					 addFunc("makeGraphFromAdjList", makeGraphFromAdjList, false),
					 addFunc("evenOutCompReads", evenOutCompReads, false),
					 addFunc("extractRefGenesFastaFiles", extractRefGenesFastaFiles, false),
					 addFunc("getMismatchDistances", getMismatchDistances, false),
					 addFunc("trimToMostProbableKmer", trimToMostProbableKmer, false),
					 addFunc("trimFromMostProbableKmer", trimFromMostProbableKmer, false),
					 addFunc("trimReadNameAtFirstWhiteSpace", trimReadNameAtFirstWhiteSpace, false),
					 addFunc("mapReads", mapReads, false),
					 addFunc("getIdDistances", getIdDistances, false),
					 addFunc("splitFile", splitFile, false),
					 addFunc("genomePrimerExtract", genomePrimerExtract, false),
					 addFunc("concatenateDifferentLanes", concatenateDifferentLanes, false),
					 addFunc("extractRefSeqsFromGenomes", extractRefSeqsFromGenomes, false),
					 addFunc("extractRefSeqsFromGenomesWithPrimers", extractRefSeqsFromGenomesWithPrimers, false),
					 addFunc("findOutliersWithMuscleToRefs", findOutliersWithMuscleToRefs, false),
					 addFunc("bioIndexGenomes", bioIndexGenomes, false),
					 addFunc("bioIndexGenome", bioIndexGenome, false),
					 addFunc("printMlnScores", printMlnScores, false),
					 addFunc("getSeqFromFile", getSeqFromFile, false),
					 addFunc("createAgreementMatrixWithMuscleRef", createAgreementMatrixWithMuscleRef, false),
					 addFunc("createAgreementSegmentsWithMuscleRef", createAgreementSegmentsWithMuscleRef, false),
					 addFunc("printMlnStreaks", printMlnStreaks, false),
					 addFunc("bowtie2ExtractAndCompare", bowtie2ExtractAndCompare, false),
					 addFunc("bowtie2ExtractAndCompareMultiple", bowtie2ExtractAndCompareMultiple, false),
					 addFunc("splitLibrary", splitLibrary, false),
					 addFunc("lastzExtractAndCompare", lastzExtractAndCompare, false),
					 addFunc("countSeqNamePortion", countSeqNamePortion, false),
           },//
          "genExp") {}





int genExpRunner::genomePrimerExtract(const njh::progutils::CmdArgs & inputCommands){
	bfs::path forwardBam = "";
	bfs::path reverseBam = "";
	bfs::path twoBitFnp = "";
	OutOptions outOpts(bfs::path(""));
	comparison allowableErrors;
	size_t insertSize = std::numeric_limits<size_t>::max();
	seqSetUp setUp(inputCommands);
	setUp.setOption(forwardBam, "--forwardBam", "Forward primer Bam", true);
	setUp.setOption(reverseBam, "--reverseBam", "Reverse primer Bam", true);
	setUp.processComparison(allowableErrors);
	setUp.setOption(twoBitFnp, "--2bit", "Twobit file", true);
	setUp.setOption(insertSize, "--insertSize", "Max Insert Size");
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	std::vector<std::shared_ptr<AlignmentResults>> forwardalnResults = gatherMapResults(
			forwardBam, twoBitFnp, allowableErrors);
	auto uniqueForwardResults = getUniqueLocationResults(forwardalnResults);

	std::vector<std::shared_ptr<AlignmentResults>> reversealnResults = gatherMapResults(
			reverseBam, twoBitFnp, allowableErrors);
	auto uniqueReverseResults = getUniqueLocationResults(reversealnResults);
	auto extracts = getPossibleGenomeExtracts(uniqueForwardResults, uniqueReverseResults, insertSize);
	std::ofstream outfile;
	std::ostream out(determineOutBuf(outfile, outOpts));

	for(const auto & extract : extracts){
		out << extract.gRegion_->genBedRecordCore().toDelimStr() << std::endl;
	}
	return 0;
}



int genExpRunner::mapReads(const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	comparison errors;
	bool byScore = false;
	bool writeOutNames = false;
	errors.distances_.eventBasedIdentity_ = 1.0;
	setUp.setOption(errors.distances_.eventBasedIdentity_, "--identity", "Percent identity to allow");
	setUp.setOption(errors.hqMismatches_, "--hqMismatches", "Number of high quality mismatches to allow");
	setUp.setOption(errors.lqMismatches_, "--lqMismatches", "Number of low quality mismatches to allow");
	setUp.setOption(errors.oneBaseIndel_, "--oneBaseIndel", "Number of one base indels to allow");
	setUp.setOption(errors.twoBaseIndel_, "--twoBaseIndel", "Number of two base indels to allow");
	setUp.setOption(errors.largeBaseIndel_, "--largeBaseIndel", "Number of large base indels to allow");
	setUp.setOption(writeOutNames, "--writeOutNames", "Write out the names for what maps to what");
	setUp.processDefaultReader(true);
	setUp.processRefFilename(true);
	setUp.processDirectoryOutputName(true);
	setUp.pars_.gapLeft_ = "0,0";
	setUp.pars_.gapRight_ = "0,0";
	setUp.pars_.gap_ = "5,1";
	setUp.processAlignerDefualts();
	setUp.setOption(byScore, "--byScore", "Determine best alignment by alignment score");
	setUp.processVerbose();
	setUp.processDebug();
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	uint64_t maxLen = 0;
	auto refSeqs = SeqInput::getReferenceSeq(setUp.pars_.refIoOptions_, maxLen);

	aligner alignerObj = aligner(maxLen, setUp.pars_.gapInfo_, setUp.pars_.scoring_,
			KmerMaps(), setUp.pars_.qScorePars_, setUp.pars_.colOpts_.alignOpts_.countEndGaps_,
			setUp.pars_.colOpts_.iTOpts_.weighHomopolyer_);
	alignerObj.processAlnInfoInput(setUp.pars_.alnInfoDirName_);
	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	readObject read;
	uint32_t pos = 0;
	std::unordered_map<std::string, std::vector<readObject>> mapCounts;

	while (reader.readNextRead(read)) {
		if(setUp.pars_.verbose_ && (pos + 1) % 100 == 0){
			std::cout << "\r" << pos + 1;
			std::cout.flush();
		}
		readVec::handelLowerCaseBases(read, setUp.pars_.ioOptions_.lowerCaseBases_);
		readVec::getMaxLength(read, maxLen);

		if (maxLen > alignerObj.parts_.maxSize_) {
			if(setUp.pars_.debug_){
				std::cout << std::endl << "Resizing alignment matrix"  << std::endl;;
			}
			alignerObj.parts_.setMaxSize(maxLen);
		}
		double bestScore = std::numeric_limits<double>::lowest();
		uint32_t bestRefPos = std::numeric_limits<uint32_t>::max();
		for (const auto & refSeqPos : iter::range(refSeqs.size())) {
			const auto & refSeq = refSeqs[refSeqPos];
			alignerObj.alignCache(refSeq, read, false);
			alignerObj.profileAlignment(refSeq, read, false, true, false);
			if (errors.passIdAndErrorThreshold(alignerObj.comp_)) {
				double currentScore = 0;
				if (byScore) {
					currentScore = alignerObj.parts_.score_;
				} else {
					currentScore = alignerObj.comp_.distances_.eventBasedIdentity_;
				}
				if (currentScore > bestScore) {
					bestScore = currentScore;
					bestRefPos = refSeqPos;
				}
			}
		}
		std::string bestRefName = "unmapped";
		if (bestRefPos != std::numeric_limits<uint32_t>::max()) {
			bestRefName = refSeqs[bestRefPos].seqBase_.name_;
		}
		mapCounts[bestRefName].emplace_back(read.seqBase_);
		++pos;
	}

	if(setUp.pars_.verbose_){
		std::cout << std::endl;
	}
	std::ofstream mapFreqInfoFile;
	openTextFile(mapFreqInfoFile, setUp.pars_.directoryName_ + "mappedFreqInfo", ".tab.txt", false, true);
	mapFreqInfoFile << "refSeq\treadsMapped\tclustersMapped\treadsFractionTotal\treadsFractionAdjusted" << std::endl;
	double totalCount = 0;
	double totalMappedCount = 0;
	for(const auto & map : mapCounts){
		totalCount += readVec::getTotalReadCount(map.second, true);
		if(map.first != "unmapped"){
			totalMappedCount += readVec::getTotalReadCount(map.second, true);
		}
	}
	for(const auto & refSeq : refSeqs){
		if(!njh::in(refSeq.seqBase_.name_, mapCounts)){
			mapCounts[refSeq.seqBase_.name_] = {};
		}
	}
	auto keys = getVectorOfMapKeys(mapCounts);
	njh::sort(keys);
	for(const auto & mapKey : keys){
		mapFreqInfoFile << mapKey
				<< "\t" << readVec::getTotalReadCount(mapCounts[mapKey],true)
				<< "\t" << mapCounts[mapKey].size()
				<< "\t" << readVec::getTotalReadCount(mapCounts[mapKey],true)/totalCount
				<< "\t" << readVec::getTotalReadCount(mapCounts[mapKey],true)/totalMappedCount << std::endl;
	}
	if(writeOutNames){
		bfs::path nameDir = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar("names", false));
		for(const auto & mapName : mapCounts){
			std::ofstream mapNamesFile;
			openTextFile(mapNamesFile, njh::files::make_path(nameDir, mapName.first).string(), ".txt", false, true);
			printVector(readVec::getNames(mapName.second), "\n", mapNamesFile);
			mapNamesFile << std::endl;
		}
	}
	if(setUp.pars_.writingOutAlnInfo_){
		alignerObj.alnHolder_.write(setUp.pars_.outAlnInfoDirName_, setUp.pars_.verbose_);
	}
	setUp.rLog_ << "Alignments Done: " << alignerObj.numberOfAlingmentsDone_ << "\n";
	return 0;
}

int genExpRunner::trimReadNameAtFirstWhiteSpace(const njh::progutils::CmdArgs & inputCommands) {
  seqSetUp setUp(inputCommands);
  setUp.processDefaultReader(true);
  setUp.finishSetUp(std::cout);
  SeqIO reader(setUp.pars_.ioOptions_);
  reader.openIn();
  reader.openOut();
  seqInfo read;
  while(reader.readNextRead(read)){
  	trimAtFirstWhitespace(read.name_);
  	reader.write(read);
  }
  return 0;
}






int genExpRunner::trimToMostProbableKmer(const njh::progutils::CmdArgs & inputCommands) {
  uint32_t kmerLength = 10;
  uint32_t trimLength = 0;
  uint32_t windowLength = 50;
  seqSetUp setUp(inputCommands);
  setUp.processDebug();
  setUp.processVerbose();
  setUp.processDefaultReader(true);
  /**@todo add a way to trim to similar condensed sequence */
  //bool useCondensed = false;
  //setUp.setOption(useCondensed, "--useCondensed", "Whether to use the condensed sequence instead");
  setUp.setOption(windowLength, "--windowLen,-w", "Window length to expand the kmer search at which to trim");
  setUp.setOption(trimLength,   "--trimLen,-t",   "Length to trim at approximately", true);
  setUp.setOption(kmerLength,   "--kmerLen,-k",   "Kmer Length");
  uint32_t testNum = 10;
  setUp.setOption(testNum,      "--testNum",   "Test number");

  setUp.finishSetUp(std::cout);
  if(windowLength > trimLength){
  	std::cerr << "Window length shouldn't be longer than trim length" << std::endl;
  	exit(1);
  }
  SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<readObject>();
  readVec::handelLowerCaseBases(inReads, setUp.pars_.ioOptions_.lowerCaseBases_);
  probabilityProfile profile(kmerLength);

  for(const auto & readPos : iter::range(len(inReads))){
  	auto & read = inReads[readPos];
  	//first turn off short sequences if they are too short and then continue on
  	if(len(read) < trimLength ){
  		read.seqBase_.on_ = false;
  		continue;
  	}
  	uint32_t startPos = trimLength - windowLength;
  	uint32_t stopPos = std::min<uint32_t>(len(read) - kmerLength + 1, trimLength + windowLength + 1);
  	for(const auto & pos : iter::range(startPos, stopPos)){
  		profile.add(read.seqBase_.seq_.substr(pos, kmerLength), false);
  	}
  }

  profile.updateProfile();
  uint32_t count = 0;
  SeqOutput writer(setUp.pars_.ioOptions_);
  writer.openOut();
  for(auto & read : inReads){
  	//skip if the length was too small
  	if(!read.seqBase_.on_){
  		continue;
  	}
  	++count;
  	uint32_t startPos = trimLength - windowLength;
  	uint32_t stopPos = std::min<uint32_t>(len(read) - kmerLength + 1, trimLength + windowLength + 1);
  	double bestProb = 0;
  	std::vector<uint32_t> bestPos;
  	for(const auto & pos : iter::range(startPos, stopPos)){
  		auto currentProb = roundDecPlaces(profile.getProbabilityOfKmer(read.seqBase_.seq_.substr(pos, kmerLength)), 10);
  		if(currentProb == bestProb){
  			bestProb = currentProb;
  			bestPos.emplace_back(pos);
  		}else if (currentProb > bestProb){
  			bestProb = currentProb;
  			bestPos.clear();
  			bestPos.emplace_back(pos);
  		}
  	}
  	auto maxLen = std::max_element(bestPos.begin(), bestPos.end());
  	if(setUp.pars_.debug_){
  		writer.writeNoCheck(read);
  	}
  	read.setClip(*maxLen);
  	writer.writeNoCheck(read);
  	if(setUp.pars_.debug_){
    	std::cout << read.seqBase_.name_ << std::endl;
    	std::cout << "bestProb: " << bestProb << std::endl;
    	std::cout << "bestPos: " << vectorToString(bestPos, ",") << std::endl;
  	}
  	if(setUp.pars_.debug_ && count > testNum){
  		break;
  	}
  }
  return 0;
}




int genExpRunner::trimFromMostProbableKmer(const njh::progutils::CmdArgs & inputCommands) {
  uint32_t kmerLength = 10;
  uint32_t windowLength = 50;
  bool keepSeq = false;
  seqSetUp setUp(inputCommands);
  setUp.processDebug();
  setUp.processVerbose();
  setUp.setOption(keepSeq, "--keepSeq", "Keep the sequence with the most probable kmer");
  setUp.processDefaultReader(true);
  /**@todo add a way to trim to similar condensed sequence */
  //bool useCondensed = false;
  //setUp.setOption(useCondensed, "--useCondensed", "Whether to use the condensed sequence instead");
  setUp.setOption(windowLength, "--windowLen,-w", "Window length to expand the kmer search at which to trim");
  setUp.setOption(kmerLength,   "--kmerLen,-k",   "Kmer Length");
  uint32_t testNum = 10;
  setUp.setOption(testNum,      "--testNum",   "Test number");

  setUp.finishSetUp(std::cout);

  SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<readObject>();
  readVec::handelLowerCaseBases(inReads, setUp.pars_.ioOptions_.lowerCaseBases_);
  probabilityProfile profile(kmerLength);

  for(const auto & readPos : iter::range(len(inReads))){
  	auto & read = inReads[readPos];
  	//first turn off short sequences if they are too short and then continue on
  	if(len(read) < windowLength ){
  		read.seqBase_.on_ = false;
  		continue;
  	}

  	uint32_t startPos = 0;
  	uint32_t stopPos = windowLength - kmerLength + 1;
  	for(const auto & pos : iter::range(startPos, stopPos)){
  		profile.add(read.seqBase_.seq_.substr(pos, kmerLength), false);
  	}
  }

  profile.updateProfile();
  uint32_t count = 0;
  SeqOutput writer(setUp.pars_.ioOptions_);
  writer.openOut();
  for(auto & read : inReads){
  	//skip if the length was too small
  	if(!read.seqBase_.on_){
  		continue;
  	}
  	++count;
  	uint32_t startPos = 0;
  	uint32_t stopPos = windowLength - kmerLength + 1;
  	double bestProb = 0;
  	std::vector<uint32_t> bestPos;
  	for(const auto & pos : iter::range(startPos, stopPos)){
  		auto currentProb = roundDecPlaces(profile.getProbabilityOfKmer(read.seqBase_.seq_.substr(pos, kmerLength)), 10);
  		if(currentProb == bestProb){
  			bestProb = currentProb;
  			if(keepSeq){
  				bestPos.emplace_back(pos);
  			}else{
  				bestPos.emplace_back(pos + kmerLength);
  			}

  		}else if (currentProb > bestProb){
  			bestProb = currentProb;
  			bestPos.clear();
  			if(keepSeq){
  				bestPos.emplace_back(pos);
  			}else{
  				bestPos.emplace_back(pos + kmerLength);
  			}
  		}
  	}
  	auto maxLen = std::max_element(bestPos.begin(), bestPos.end());
  	if(setUp.pars_.debug_){
  		writer.writeNoCheck(read);
  	}
  	readVecTrimmer::trimOffForwardBases(read,*maxLen);
  	writer.writeNoCheck(read);
  	if(setUp.pars_.debug_){
    	std::cout << read.seqBase_.name_ << std::endl;
    	std::cout << "bestProb: " << bestProb << std::endl;
    	std::cout << "bestPos: " << vectorToString(bestPos, ",") << std::endl;
  	}
  	if(setUp.pars_.debug_ && count > testNum){
  		break;
  	}
  }
  return 0;
}





int genExpRunner::getMismatchDistances(const njh::progutils::CmdArgs & inputCommands){
	seqSetUp setUp(inputCommands);
	setUp.processDefaultReader(true);
	setUp.processAlignerDefualts();
	setUp.finishSetUp(std::cout);
	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<readObject>();

	uint64_t maxSize = 0;
	readVec::getMaxLength(inReads, maxSize);
	aligner alignerObj(maxSize, setUp.pars_.gapInfo_, setUp.pars_.scoring_, setUp.pars_.colOpts_.alignOpts_.countEndGaps_);
	alignerObj.processAlnInfoInput(setUp.pars_.alnInfoDirName_);
	std::vector<uint32_t> mismatches;
	for(const auto & firstPos : iter::range(inReads.size())){
		for(const auto & secondPos : iter::range(firstPos)){
			std::cout << firstPos << ":" << secondPos << std::endl;
			alignerObj.alignCache(inReads[firstPos], inReads[secondPos], setUp.pars_.local_);
			alignerObj.profilePrimerAlignment(inReads[firstPos], inReads[secondPos]);
			mismatches.emplace_back(alignerObj.comp_.hqMismatches_ + alignerObj.comp_.lqMismatches_ + alignerObj.comp_.lowKmerMismatches_);
		}
	}
	auto stats = getStatsOnVec(mismatches);
	table mTab(stats,VecStr{"stat", "value"});
	mTab.outPutContentOrganized(std::cout);
	return 0;
}

int genExpRunner::getIdDistances(const njh::progutils::CmdArgs & inputCommands){
	seqSetUp setUp(inputCommands);
	uint32_t numThreads = 1;
	TableIOOpts outOptions = TableIOOpts::genTabFileOut("out.tab.txt", true);
	setUp.setOption(numThreads, "--numThreads", "Number of threads to use");
	setUp.processDefaultReader(true);
	setUp.processAlignerDefualts();
	setUp.pars_.ioOptions_.out_ = outOptions.out_;
	setUp.finishSetUp(std::cout);
	SeqInput reader(setUp.pars_.ioOptions_);
	auto inReads = reader.readAllReads<readObject>();

	uint64_t maxSize = 0;
	readVec::getMaxLength(inReads, maxSize);
	aligner alignerObj(maxSize, setUp.pars_.gapInfo_, setUp.pars_.scoring_,
			setUp.pars_.colOpts_.alignOpts_.countEndGaps_);
	alignerObj.processAlnInfoInput(setUp.pars_.alnInfoDirName_, setUp.pars_.verbose_);
	std::unordered_map<std::string, std::shared_ptr<aligner>> aligners;
	std::mutex alignerLock;
	std::function<
			double(const readObject &,
					const readObject &,
					std::unordered_map<std::string, std::shared_ptr<aligner>>&, aligner&)> alignFunc = [&alignerLock](const readObject & read1, const readObject & read2,
							std::unordered_map<std::string, std::shared_ptr<aligner>>& aligners,
							aligner &alignerObj) {
						alignerLock.lock();
						auto threadId = estd::to_string(std::this_thread::get_id());
						//std::cout << threadId<< std::endl;
						if(aligners.find(threadId) == aligners.end()) {
							aligners.emplace(threadId, std::make_unique<aligner>(alignerObj));
						}
						alignerLock.unlock();
						aligners.at(threadId)->alignCache(getSeqBase(read1),getSeqBase(read2), false);
						aligners.at(threadId)->profilePrimerAlignment(getSeqBase(read1), getSeqBase(read2));
						return 1 - aligners.at(threadId)->comp_.distances_.eventBasedIdentity_;
					};

	auto distances = getDistanceNonConst(inReads, numThreads, alignFunc,
						aligners, alignerObj);

	for(const auto & aligObj : aligners){
		alignerObj.alnHolder_.mergeOtherHolder(aligObj.second->alnHolder_);
	}
	alignerObj.processAlnInfoOutput(setUp.pars_.outAlnInfoDirName_, setUp.pars_.verbose_);

	std::ofstream outFile;
	openTextFile(outFile, setUp.pars_.ioOptions_.out_);

	writeDistanceMatrix(outFile, distances, readVec::getNames(inReads));

	return 0;
}


int genExpRunner::extractRefGenesFastaFiles(const njh::progutils::CmdArgs & inputCommands){
  seqSetUp setUp(inputCommands);
  std::string filename = "";
  std::string outDir = "";
  std::string twoBitFilename = "";
  setUp.setOption(twoBitFilename, "--twoBit", "Name of the two bit file", true);
  setUp.setOption(filename, "--file", "Filename of Mip Position File", true);
  setUp.setOption(outDir, "--outDir", "outDirectory", true);
  setUp.processWritingOptions();
  setUp.finishSetUp(std::cout);


  TwoBit::TwoBitFile twoBitFile(twoBitFilename);
  auto seqNames = twoBitFile.sequenceNames();
  auto refSeqInfo = getRefSeqRecs(filename, VecStr{});
  for(const auto & record : refSeqInfo){
  	auto chrName = record.second->chrom_;
  	if(!njh::in(chrName,seqNames)){
  		std::cerr << "chromosome name not found in seq names, skipping" << std::endl;
  		std::cerr << "chr: " << chrName << std::endl;
  		std::cerr << "possibleNames: " << vectorToString(seqNames,",") << std::endl;
  	}else{
  		std::string seq = "";
  		uint64_t start = record.second->txStart_;
  		uint64_t stop = record.second->txEnd_;
  		bool revrseStrand = (record.second->strand_ == '-');
  		twoBitFile[chrName]->getSequence(seq, start, stop);
  		if(revrseStrand){
  			seq = seqUtil::reverseComplement(seq, "DNA");
  		}
  		if(outDir.back() != '/'){
  			outDir.push_back('/');
  		}
  	  std::ofstream outFile;
  	  openTextFile(outFile, outDir + record.second->name_ + "_genomic", ".fasta", setUp.pars_.ioOptions_.out_);
  		outFile << ">" << record.second->name_ << std::endl;
  		outFile << seq << std::endl;
  	}
  }
  return 0;
}


int genExpRunner::evenOutCompReads(const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processDefaultReader(true);
	if(setUp.pars_.ioOptions_.out_.outFilename_ == "out"){
		setUp.pars_.ioOptions_.out_.outFilename_ = "evened_" + bfs::basename(setUp.pars_.ioOptions_.firstName_);
	}
	setUp.finishSetUp(std::cout);

	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<readObject>();

	std::vector<uint32_t> compPositions;
	std::vector<uint32_t> forPositons;
	for(const auto & readPos : iter::range(inReads.size())){
		if(njh::containsSubString(inReads[readPos].seqBase_.name_, "Comp")){
			compPositions.emplace_back(readPos);
		}else{
			forPositons.emplace_back(readPos);
		}
	}
	SeqOutput writer(setUp.pars_.ioOptions_);
	writer.openOut();
	std::vector<uint32_t> outPositons;
	njh::randomGenerator gen;
	if(setUp.pars_.debug_){
		gen.seedNum(0);
	}
	if(forPositons.size() == compPositions.size()){
		outPositons = concatVecs(compPositions, forPositons);
	}else if (forPositons.size() > compPositions.size()){
		outPositons = compPositions;
		addOtherVec(outPositons, gen.unifRandSelectionVec(forPositons, compPositions.size(), false));
	}else if (compPositions.size() > forPositons.size()){
		outPositons = forPositons;
		addOtherVec(outPositons, gen.unifRandSelectionVec(compPositions, forPositons.size(), false));
	}else{
		//shouldn't be here ...
		std::cout << "Well this shouldn't be happening" << std::endl;
	}
	for(const auto & pos : outPositons){
		writer.writeNoCheck(inReads[pos]);
	}
	return 0;
}






int genExpRunner::makeGraphFromAdjList(const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	std::string filename = "";
	std::string backgroundColor = "#000000";
	double hueStart = 120;
	double hueStop = 420;
	double lumStart = 0.40;
	double lumStop = 0.70;
	double satStart = 0.80;
	double satStop = 1.00;
	setUp.setOption(filename, "--filename", "Name of a file that contains an adjacency list");
	setUp.processDefaultReader(true);
	setUp.processDirectoryOutputName(true);
	setUp.setOption(backgroundColor, "-b,--backgroundColor", "Hex String for the color for the Background");
	setUp.setOption(hueStart, "--hueStart", "Hue Start for the color of Reads");
	setUp.setOption(hueStop, "--hueStop", "Hue Stop for the color of Reads");
	setUp.setOption(satStart, "--satStart", "Sat Start for the color of Reads");
	setUp.setOption(satStop, "--satStop", "Sat Stop for the color of Reads");
	setUp.setOption(lumStart, "--lumStart", "Lum Start for the color of Reads");
	setUp.setOption(lumStop, "--lumStop", "Lum Stop for the color of Reads");
	setUp.processDebug();
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<readObject>();

	std::ifstream inFile(filename);
	if(!inFile){
		std::cerr <<njh::bashCT::bold << "Error in opening " << njh::bashCT::red << filename << njh::bashCT::reset << std::endl;
		exit(1);
	}
	std::string line = "";
	std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>> readDists;
	while (std::getline(inFile, line)){
		auto lPos = line.find_last_not_of(" \f\n\r\t\v");
		auto sPos = line.find_last_of(" \f\n\r\t\v",lPos);
		auto distStr = line.substr(sPos);
		uint32_t dist = estd::stou(distStr);
		auto nodesStr = line.substr(0, sPos);
		auto toks = tokenizeString(nodesStr, "--");
		readDists[trimEndWhiteSpaceReturn(toks[0])][trimEndWhiteSpaceReturn(toks[1])] = dist;
		/**@todo do safe parse checking */
	}
	if(setUp.pars_.debug_){
		for(const auto & rn1 : readDists){
			for(const auto & rn2 : rn1.second){
				std::cout << rn1.first << " -- " << rn2.first << " " << rn2.second << std::endl;
			}
		}
		exit(1);
	}


  readDistGraph<uint32_t> graphMis(readDists, inReads);

	std::vector<std::string> names;
  for(const auto & n : graphMis.nodes_){
  	names.emplace_back(n->name_);
  }
  if(hueStop == 360 && hueStart == 0){
  	hueStop = 360 - (360.0/names.size());
  }
  for(auto pos : iter::range(graphMis.nodes_.size())){
  	graphMis.nodes_[pos]->group_ = pos;

  }
  auto graphJson = graphMis.toJsonMismatchGraph(njh::color(backgroundColor), hueStart, hueStop,
    			lumStart, lumStop, satStop, satStart);

	std::ofstream outFile(setUp.pars_.directoryName_ + "psudoMinTree.json");
	outFile << graphJson << std::endl;
	std::string htmlOut = genHtmlStrForPsuedoMintree("psudoMinTree.json");
	std::ofstream outHtmlFile(setUp.pars_.directoryName_ + "psudoMinTree.html");
	outHtmlFile << htmlOut << std::endl;
	return 0;
}

int genExpRunner::splitUpFile(const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	uint32_t sampleSize = 50000;
	setUp.setOption(sampleSize, "-sampleSize", "Sample Size");
	setUp.processDefaultReader(true);
	setUp.processDirectoryOutputName(true);
	setUp.finishSetUp(std::cout);
	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	readObject read;
	uint32_t count = 0;
	uint32_t sampleCount = 0;
	uint32_t numSamples = 0;
	MultiSeqIO multiWriter;
	multiWriter.addReader("sample_" + estd::to_string(numSamples),
			SeqIOOptions(
					setUp.pars_.directoryName_ + "sample_" + estd::to_string(numSamples),
					setUp.pars_.ioOptions_.outFormat_, setUp.pars_.ioOptions_.out_));
	std::string currentSample = "sample_" + estd::to_string(numSamples);
	while (reader.readNextRead(read)) {
		std::cout << count << "\r";
		std::cout.flush();
		multiWriter.openWrite(currentSample, read);
		++count;
		++sampleCount;
		if(sampleCount == sampleSize){
			sampleCount = 0;
			++numSamples;
			multiWriter.addReader("sample_" + estd::to_string(numSamples),
					SeqIOOptions(
							setUp.pars_.directoryName_ + "sample_" + estd::to_string(numSamples),
							setUp.pars_.ioOptions_.outFormat_, setUp.pars_.ioOptions_.out_));
			std::string currentSample = "sample_" + estd::to_string(numSamples);
		}
	}
	std::cout << std::endl;
	return 0;
}

int genExpRunner::extractOnLen(const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	uint32_t minLen = 100;
	uint32_t maxLen = 400;
	//uint32_t qualCheck = 30;
	//double qualCheckCutOff = .85;

	setUp.setOption(minLen, "--minLen", "Minimum Length");
	setUp.setOption(maxLen, "--maxLen", "Maximum Length");
	setUp.processDefaultReader(true);
	setUp.finishSetUp(std::cout);
	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	readObject read;
	MultiSeqIO multiWriter;
	multiWriter.addReader("within",
			SeqIOOptions(bfs::path(setUp.pars_.ioOptions_.out_.outFilename_.string() + "_within"),
					setUp.pars_.ioOptions_.outFormat_, setUp.pars_.ioOptions_.out_));
	multiWriter.addReader("out",
			SeqIOOptions(bfs::path(setUp.pars_.ioOptions_.out_.outFilename_.string() + "out"),
					setUp.pars_.ioOptions_.outFormat_, setUp.pars_.ioOptions_.out_));
	uint32_t count = 0;
	while(reader.readNextRead(read)){
		std::cout << count << "\r";
		std::cout.flush();
		if(read.seqBase_.seq_.length() >= minLen && read.seqBase_.seq_.length() <= maxLen){
			multiWriter.openWrite("within",read);
		}else{
			multiWriter.openWrite("out",read);
		}
		++count;
	}
	std::cout << std::endl;
	return 0;
}
int genExpRunner::sortOnNucComp(const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	setUp.processDefaultReader(true);
	setUp.finishSetUp(std::cout);
	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<readObject>();

	readVec::allSetLetterCount(inReads);
	uint64_t maxSize = 0;
	readVec::getMaxLength(inReads, maxSize);
	njh::randomGenerator gen;
	readObject rSeq = readObject(
			seqInfo("rSeq",
					simulation::evenRandStr(maxSize, { 'A', 'C', 'G', 'T' }, gen)));
	rSeq.counter_.setFractions();
	auto nucComp =
			[&rSeq](const readObject & rSeq1, const readObject & rSeq2) {return rSeq1.counter_.getFracDifference(rSeq.counter_, rSeq.counter_.alphabet_) < rSeq2.counter_.getFracDifference(rSeq.counter_, rSeq.counter_.alphabet_);};
	njh::sort(inReads, nucComp);
	SeqOutput::write(inReads, setUp.pars_.ioOptions_);
	return 0;
}


int genExpRunner::getSnpInfo(const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	uint32_t numThreads = 2;
	std::string testFastq = "back.fastq";
	setUp.setOption(numThreads, "-t,--numThreads", "Number of Threads to Use");
	setUp.setOption(testFastq, "--testFastq", "testFastq");
	setUp.processDefaultReader(true);
	setUp.processAlignerDefualts();
	setUp.finishSetUp(std::cout);
	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<readObject>();

	uint64_t maxLen = 0;
	readVec::getMaxLength(inReads, maxLen);
	SeqIOOptions testOpts = SeqIOOptions::genFastqIn(testFastq);
	SeqInput testReader(testOpts);
	testReader.openIn();
	auto testinReads = testReader.readAllReads<readObject>();
	readVec::getMaxLength(testinReads, maxLen);
	aligner alignerObj(maxLen, setUp.pars_.gapInfo_, setUp.pars_.scoring_,
			KmerMaps(), setUp.pars_.qScorePars_, setUp.pars_.colOpts_.alignOpts_.countEndGaps_,
			setUp.pars_.colOpts_.iTOpts_.weighHomopolyer_);
	alignerObj.processAlnInfoInput(setUp.pars_.alnInfoDirName_);
  std::function<uint32_t(const readObject & ,
  		const readObject &, aligner)> misFun = getMismatches<readObject>;
	auto misDistances = getDistanceCopy(inReads, numThreads, misFun,
			alignerObj);
  readDistGraph<uint32_t> graphMis(misDistances, inReads);
  graphMis.turnOffEdgesAbove(1);
  graphMis.determineGroups();
  graphMis.printAdjByGroup(std::cout);
  auto gValues = graphMis.getGroupValues();
  auto values = gValues[0];
  std::unordered_map<std::string,std::unordered_map<uint64_t,std::unordered_map<char, VecStr>>> snps;
  std::unordered_map<std::string, std::vector<uint64_t>> importantLocs;
  for(const auto & j : values){
  	for(const auto & i : values){
  		if(j != i){
  			alignerObj.alignCache(*j, *i, false);
  			alignerObj.profilePrimerAlignment(*j, *i);
  			for(const auto & m : alignerObj.comp_.distances_.mismatches_){
  				snps[j->name_][m.second.refBasePos][m.second.seqBase].emplace_back(i->name_);
  			}
  		}
  	}
  }
	for (const auto & fName : snps) {
		importantLocs[fName.first] = getVectorOfMapKeys(fName.second);
	}

	std::vector<readObject> refReads;
	for(const auto & v : values){
		refReads.emplace_back(*v);
	}

	table outExampleTab { VecStr { "ref", "read", "pos", "base", "otherRef" } };
	//std::unordered_map<std::string, std::vector<uint64_t>> readMatchingRefs;
	std::unordered_map<std::string, std::vector<std::string>> readMatchingRefs;
	for (const auto & readPos : iter::range(testinReads.size())) {
		const auto & read = testinReads[readPos];
		std::cout << read.seqBase_.name_ << " " << readPos << ":" << testinReads.size() << std::endl;;
		//std::cout.flush();
		//std::vector<uint32_t> matchingRefs;
		std::vector<std::string> matchingRefs;
		for (const auto & refPos : iter::range(refReads.size())) {
			const auto & ref = refReads[refPos];
			alignerObj.alignCache(ref, read, false);
			alignerObj.profilePrimerAlignment(ref, read);
			std::unordered_map<uint64_t, char> currentSnps;
			for (const auto & m : alignerObj.comp_.distances_.mismatches_) {
				currentSnps[m.second.refBasePos] = m.second.seqBase;
			}
			bool pass = true;
			for(const auto & loc : importantLocs[ref.seqBase_.name_]){
				auto search = currentSnps.find(loc);
				if((search == currentSnps.end() || snps[ref.seqBase_.name_][loc][search->second].empty()) &&
						alignerObj.alignObjectB_.seqBase_.seq_[alignerObj.getAlignPosForSeqAPos(loc)] != '-'){
				}else{
					pass = false;
				}
			}
			if(pass){
				matchingRefs.emplace_back(ref.seqBase_.name_);
			}
		}
		readMatchingRefs[read.seqBase_.name_] = matchingRefs;
	}
	std::cout << std::endl;
	table outTable { VecStr { "FirstName", "SecondName", "position", "base" } };
	for (const auto & fName : snps) {
		for (const auto & sName : fName.second) {
			for (const auto & pos : sName.second) {
				outTable.content_.emplace_back(
						toVecStr(fName.first, sName.first, pos.first,
								vectorToString(pos.second, ",")));
			}
		}
	}
	for (const auto & c : iter::reversed(outTable.columnNames_)) {
		outTable.sortTable(c, false);
	}
	outTable.outPutContentOrganized(std::cout);

	for (const auto & c : iter::reversed(outExampleTab.columnNames_)) {
		outExampleTab.sortTable(c, false);
	}
	outExampleTab.sortTable("read", false);
	outExampleTab.outPutContentOrganized(std::cout);
	table outMatchingRefsTab{VecStr{"read", "refs"}};
	for(const auto & r : readMatchingRefs){
		outMatchingRefsTab.content_.emplace_back(toVecStr(r.first, vectorToString(r.second, ",")));
	}
	outMatchingRefsTab.sortTable("read", false);
	outMatchingRefsTab.outPutContentOrganized(std::cout);
	return 0;
}

int genExpRunner::findVariableSites(const njh::progutils::CmdArgs & inputCommands){
  seqSetUp setUp(inputCommands);
  bool ignoreNs = false;
  bool ignoreGaps = false;
  bool writeFastaVars = false;
  std::string outFilename = "";
  setUp.processAlignerDefualts();
  setUp.pars_.ioOptions_.out_.outFilename_ = "out";
  setUp.processDefaultReader(true);
  setUp.setOption(ignoreNs, "-ignoreNs", "Ignore Ns");
  setUp.setOption(ignoreGaps, "-ignoreGaps", "Ignore Gaps");
  setUp.setOption(outFilename, "-varOutFilename", "Out Filename");
  setUp.setOption(writeFastaVars, "-writeFastaVars", "Write A Fasta File of only the Var Sites");
  setUp.finishSetUp(std::cout);
  SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<readObject>();
  uint32_t seqLength = len(inReads.front().seqBase_.seq_);
  std::vector<charCounter> counters(seqLength, charCounter());
  for(const auto & pos : iter::range(seqLength)){
  	for(const auto & read : inReads){
    	counters[pos].increaseCountOfBase(read.seqBase_.seq_[pos]);
  	}
  }
	charCounter allCounter;
	njh::for_each(inReads,
			[&](const readObject & read) {allCounter.increaseCountByString(read.seqBase_.seq_, read.seqBase_.cnt_);});
	if (ignoreNs) {
    for(const auto & pos : iter::range(seqLength)){
    	counters[pos].chars_['N'] = 0;
    	counters[pos].chars_['n'] = 0;
    }
    allCounter.chars_['N'] = 0;
    allCounter.chars_['n'] = 0;
  }
  if(ignoreGaps){
    for(const auto & pos : iter::range(seqLength)){
    	counters[pos].chars_['-'] = 0;
    	counters[pos].chars_['-'] = 0;
    }
    allCounter.chars_['-'] = 0;
    allCounter.chars_['-'] = 0;
  }
  allCounter.resetAlphabet(false);
  {
  	table outInfo(VecStr{toVecStr("pos", allCounter.alphabet_)});
    for(const auto & pos : iter::range(seqLength)){
    	counters[pos].resetAlphabet(false);
    	if(counters[pos].alphabet_.size() > 1){
      	VecStr currentRow;
      	currentRow.reserve(1 + counters[pos].alphabet_.size());
      	currentRow.emplace_back(estd::to_string(pos));
    		counters[pos].setFractions(allCounter.alphabet_);
  			for(const auto & c : allCounter.alphabet_){
  				currentRow.emplace_back(estd::to_string(counters[pos].fractions_[c]));
  			}
  			outInfo.content_.emplace_back(currentRow);
    	}
		}
		outInfo.outPutContents(
				TableIOOpts { OutOptions(outFilename, ".tab.txt", "tab", false,
						setUp.pars_.ioOptions_.out_.overWriteFile_, false), "\t",
						outInfo.hasHeader_ });

    if(writeFastaVars){
    	std::vector<uint32_t> positions = vecStrToVecNum<uint32_t>(outInfo.getColumn("pos"));
    	std::vector<readObject> outReads;
    	outReads.reserve(inReads.size());
    	for(auto & read : inReads){
    		seqInfo info{read.seqBase_.name_};
    		info.seq_.reserve(positions.size());
    		info.qual_.reserve(positions.size());

    		for(const auto & pos : positions){
    			info.append(read.seqBase_.seq_[pos], read.seqBase_.qual_[pos]);
    		}
    		outReads.emplace_back(info);
    	}
    	SeqOutput::write(outReads, setUp.pars_.ioOptions_);
    }
  }
  if(false){
    std::stringstream tempStream;
    tempStream << "pos\tbase\tbaseFrac" << std::endl;
    for(const auto & pos : iter::range(seqLength)){
    	counters[pos].resetAlphabet(false);
    	if(counters[pos].alphabet_.size() > 1){
    		counters[pos].setFractions(allCounter.alphabet_);
  			for(const auto & c : allCounter.alphabet_){
  				tempStream << pos << "\t" << c
  						<< "\t" << counters[pos].fractions_[c] << std::endl;
  			}
    	}
    }
    table outInfo(tempStream,"\t", true);
    outInfo.outPutContentOrganized(std::cout);
  }
  return 0;
}

int genExpRunner::extractByName(const njh::progutils::CmdArgs & inputCommands){
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


int genExpRunner::muscle(const njh::progutils::CmdArgs & inputCommands){
	seqSetUp setUp(inputCommands);
	setUp.pars_.ioOptions_.out_.outFilename_ = "out";
	setUp.processDefaultReader(true);
	setUp.finishSetUp(std::cout);
	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<readObject>();
	Muscler musclerOperator;
	musclerOperator.muscleSeqs(inReads);
	SeqOutput::write(inReads, setUp.pars_.ioOptions_);
	return 0;
}

seqInfo calculateConsensusToCurrent(const readObject & read,
		std::vector<readObject> & reads, aligner& alignerObj, bool setToConsensus) {
	// if the cluster is only one read, no need to create consensus
	// if (reads.size() <= 1) {
	//return;
	// }
	// create the map for letter counters for each position
	std::map<uint32_t, charCounter> counters;
	// create a map in case of insertions
	std::map<uint32_t, std::map<uint32_t, charCounter>> insertions;
	std::map<int32_t, charCounter> beginningGap;
	for (const auto & readPos : iter::range(reads.size())) {
		alignerObj.alignCache(read, reads[readPos], false);
		// the offset for the insertions
		uint32_t offSet = 0;
		uint32_t currentOffset = 1;
		uint32_t start = 0;
		//check to see if there is a gap at the beginning
		if (alignerObj.alignObjectA_.seqBase_.seq_.front() == '-') {
			start = alignerObj.alignObjectA_.seqBase_.seq_.find_first_not_of('-');
			for (uint32_t i = 0; i < start; ++i) {
				beginningGap[i - start].increaseCountOfBaseQual(
						alignerObj.alignObjectB_.seqBase_.seq_[i],
						alignerObj.alignObjectB_.seqBase_.qual_[i],
						reads[readPos].seqBase_.cnt_);
				++offSet;
			}
		}
		for (uint32_t i = start; i < len(alignerObj.alignObjectB_); ++i) {
			// if the longest reference has an insertion in it put it in the
			// insertions letter counter map
			if (alignerObj.alignObjectA_.seqBase_.seq_[i] == '-') {
				insertions[i - offSet][currentOffset].increaseCountOfBaseQual(
						alignerObj.alignObjectB_.seqBase_.seq_[i],
						alignerObj.alignObjectB_.seqBase_.qual_[i],
						reads[readPos].seqBase_.cnt_);

				++currentOffset;
				++offSet;
				continue;
			}
			currentOffset = 1;
			counters[i - offSet].increaseCountOfBaseQual(
					alignerObj.alignObjectB_.seqBase_.seq_[i],
					alignerObj.alignObjectB_.seqBase_.qual_[i],
					reads[readPos].seqBase_.cnt_);
		}
	}
	std::string calculatedConsensus = "";
	calculatedConsensus.reserve(len(read));
	std::vector<uint32_t> calculatedConsensusQuality;
	calculatedConsensusQuality.reserve(len(read));
	// first deal with any gaps in the beginning
	double fortyPercent = 0.40 * read.seqBase_.cnt_;
	for(const auto & bCount : beginningGap){
		uint32_t bestQuality = 0;
		char bestBase = ' ';
		bCount.second.getBest(bestBase, bestQuality);
		if (bestBase == '-' || bCount.second.getTotalCount() < fortyPercent) {
			continue;
		}
		calculatedConsensus.push_back(bestBase);
		calculatedConsensusQuality.emplace_back(bestQuality / bCount.second.getTotalCount());
	}
	//read.seqBase_.outPutFastq(std::cout);
	// the iterators to over the letter counter maps
	for (const auto & count : counters) {
		uint32_t bestQuality = 0;
		char bestBase = ' ';
		// if there is an insertion look at those if there is a majority of reads
		// with that insertion
		auto search = insertions.find(count.first);
		if (search != insertions.end()) {
			for (auto & counterInsert : search->second) {
				bestQuality = 0;
				bestBase = ' ';
				counterInsert.second.getBest(bestBase, bestQuality,
						std::round(read.seqBase_.cnt_));
				if (bestBase == ' ') {
					continue;
				} else {
					calculatedConsensus.push_back(bestBase);
					calculatedConsensusQuality.emplace_back(bestQuality);
				}
			}
		}
		count.second.getBest(bestBase, bestQuality);
		if (bestBase == '-' || count.second.getTotalCount() < fortyPercent) {
			continue;
		}
		calculatedConsensus.push_back(bestBase);
		calculatedConsensusQuality.emplace_back(
				bestQuality / count.second.getTotalCount());
	}
	//if (setToConsensus) {
	//  seqBase_.seq_ = calculatedConsensus;
	// seqBase_.qual_ = calculatedConsensusQuality;
	//}
	//previousErrorChecks.clear();
	//needToCalculateConsensus = false;
	return seqInfo { read.seqBase_.name_, calculatedConsensus,
			calculatedConsensusQuality, read.seqBase_.cnt_ };
}

int genExpRunner::createReadConsensus(const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	bool sort = false;
	std::string name = "Consensus";
	setUp.setOption(name, "-name", "Name of Consensus");
	setUp.setOption(sort, "-s,--sort", "Sort Reads First");
	setUp.processAlignerDefualts();
	setUp.processDefaultReader(true);
	setUp.finishSetUp(std::cout);


	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<readObject>();

	uint64_t maxSize = 0;
	readVec::getMaxLength(inReads, maxSize);
	std::vector<baseCluster> inClusters;
	if(sort){
		readVecSorter::sort(inReads, true);
	}
	readObject firstRead = readObject(inReads.begin()->seqBase_);
	baseCluster mainCluster(firstRead.seqBase_);
	for(const auto & read : inReads){
		inClusters.emplace_back(baseCluster(read.seqBase_));
	}
	for(const auto & clus : inClusters){
		mainCluster.addRead(clus);
	}
	aligner alignerObj = aligner(maxSize, setUp.pars_.gapInfo_,
			substituteMatrix(setUp.pars_.generalMatch_, setUp.pars_.generalMismatch_),
			KmerMaps(), setUp.pars_.qScorePars_,setUp.pars_.colOpts_.alignOpts_.countEndGaps_,
			setUp.pars_.colOpts_.iTOpts_.weighHomopolyer_);

	mainCluster.calculateConsensus(alignerObj, true);

	if (setUp.pars_.ioOptions_.out_.outFilename_ == "") {
		mainCluster.seqBase_.outPutFastq(std::cout);
	} else {
		SeqOutput writer(setUp.pars_.ioOptions_);
		writer.openWrite(mainCluster);
	}

  return 0;
}

int genExpRunner::createReadConsensusRandomPicking(const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	std::string name = "Consensus";
	uint32_t times = 10;
	setUp.setOption(times, "-times", "Number of times to randomly pick a consensus and create one");
	setUp.setOption(name, "-name", "Name of Consensus");
	setUp.processAlignerDefualts();

	setUp.processDefaultReader(true);
	setUp.finishSetUp(std::cout);


	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<readObject>();

	uint64_t maxSize = 0;
	readVec::getMaxLength(inReads, maxSize);
	aligner alignerObj = aligner(maxSize, setUp.pars_.gapInfo_,
			substituteMatrix(setUp.pars_.generalMatch_, setUp.pars_.generalMismatch_),
			KmerMaps(), setUp.pars_.qScorePars_,setUp.pars_.colOpts_.alignOpts_.countEndGaps_,
			setUp.pars_.colOpts_.iTOpts_.weighHomopolyer_);
	alignerObj.processAlnInfoInput(setUp.pars_.alnInfoDirName_);
	njh::randomGenerator gen;
	std::ofstream outFile;
	SeqOutput output(setUp.pars_.ioOptions_);
	output.openOut();
	for(uint32_t time = 0; time < times; ++time	){
		std::vector<baseCluster> inClusters;
		readObject firstRead;
		auto ranCon = gen.unifRand<uint64_t>(0, inReads.size());
		for(const auto & readPos : iter::range<uint32_t>(inReads.size())){
			if(readPos == ranCon){
				inClusters.emplace_back(baseCluster(inReads[readPos].seqBase_));
			}else{
				firstRead = readObject(inReads[readPos].seqBase_);
			}
		}
		baseCluster mainCluster(firstRead.seqBase_);
		for(const auto & clus : inClusters){
			mainCluster.addRead(clus);
		}
		mainCluster.calculateConsensusToCurrent(alignerObj, true);
		output.write(mainCluster.seqBase_);
	}
	if(setUp.pars_.writingOutAlnInfo_){
		alignerObj.alnHolder_.write(setUp.pars_.alnInfoDirName_);
	}

  return 0;
}//getKmerDistStats



} /* namespace njhseq */


namespace njhseq {






std::vector<readObject> seqMapCountsToSeqs(const std::map<std::string, uint32_t> & counts, const std::string & name){
	std::vector<readObject> ret;
	ret.reserve(counts.size());
	for(const auto & cEnum : iter::enumerate(counts)){
		ret.emplace_back(seqInfo(name + "." + estd::to_string(cEnum.index), cEnum.element.first));
		ret.back().seqBase_.cnt_ = cEnum.element.second;
		ret.back().updateName();
	}
	return ret;
}

std::vector<readObject> vecStrMultiplesToSeqs(const VecStr & strs, const std::string & name){
	auto counts = countVec(strs);
	return seqMapCountsToSeqs(counts, name);
}

int genExpRunner::tableCountToFasta(const njh::progutils::CmdArgs & inputCommands){
  seqSetUp setUp(inputCommands);
  std::string filename = "";
  std::string nameStub = "Bar";
  std::string column = "";
  setUp.setOption(filename, "-i,--inFile", "Filename of a table of first column seq and second column counts, or a table to count a column of", true);
  setUp.setOption(nameStub, "--name", "Name Stuf for seqs");
  setUp.setOption(column, "-c,--column", "name of the column to count and create seqs from");
  setUp.finishSetUp(std::cout);

  table inTab(filename, "\t", true);
  std::vector<readObject> seqs;
  if(column == ""){
    std::map<std::string, uint32_t> counts;
    for(const auto & row : inTab.content_){
    	counts[row[0]] = estd::stou(row[1]);
    }
    seqs = seqMapCountsToSeqs(counts, nameStub);
  }else{
  	seqs = vecStrMultiplesToSeqs(inTab.getColumn(column), nameStub);
  }
  readVec::allPrintSeqs(seqs);
	return 0;
}


uint32_t countBarcodeMismatches(const std::string & bar1,
		const std::string & bar2){
	uint32_t ret = 0;
	for(auto pos : iter::range(bar1.size())){
		if(bar1[pos] != bar2[pos]){
			++ret;
		}
	}
	return ret;
}
template<typename T>
uint32_t countBarcodeMismatches(const T & bar1,
		const T & bar2){
	uint32_t ret = 0;
	for(auto pos : iter::range(bar1.seqBase_.seq_.size())){
		if(bar1.seqBase_.seq_[pos] != bar2.seqBase_.seq_[pos]){
			++ret;
		}
	}
	return ret;
}


int genExpRunner::tableCountToDistGraph(const njh::progutils::CmdArgs & inputCommands){
  seqSetUp setUp(inputCommands);
  std::string filename = "";
	std::string backgroundColor = "#000000";
	uint32_t mismatches = 0;
  std::string nameStub = "Bar";
  std::string column = "";
	uint32_t numThreads = 1;
	double hueStart = 0;
	double hueStop = 360;
	double satStart = 1.0;
	double satStop = 1.0;
	double lumStart = 0.5;
	double lumStop = 0.5;
	uint32_t groupCutOff = 2;
	setUp.setOption(mismatches, "--mismatches", "Number of mismatches to allow");
	setUp.setOption(numThreads, "-t,--numThreads", "Number of Threads to Use");
  setUp.setOption(filename, "-i,--inFile", "Filename of a table of first column seq and second column counts, or a table to count a column of", true);
  setUp.setOption(nameStub, "--name", "Name Stuf for seqs");
  setUp.setOption(column, "-c,--column", "name of the column to count and create seqs from");
  setUp.processDirectoryOutputName(filename + "_" + inputCommands.getProgramName() + "_TODAY", true);
  setUp.processVerbose();
  setUp.finishSetUp(std::cout);

  table inTab(filename, "\t", true);
  std::vector<readObject> seqs;
  if(column == ""){
    std::map<std::string, uint32_t> counts;
    for(const auto & row : inTab.content_){
    	counts[row[0]] = estd::stou(row[1]);
    }
    seqs = seqMapCountsToSeqs(counts, nameStub);
  }else{
  	seqs = vecStrMultiplesToSeqs(inTab.getColumn(column), nameStub);
  }
  njh::stopWatch watch;

  watch.startNewLap("regular distance");

  std::function<uint32_t(const readObject & ,
  		const readObject &)> misFun = countBarcodeMismatches<readObject>;
	auto misDistances = getDistance(seqs, numThreads, misFun);
	watch.startNewLap("create graph");
  readDistGraph<uint32_t> graphMis(misDistances, seqs);
	std::vector<std::string> names;
  for(const auto & n : graphMis.nodes_){
  	names.emplace_back(n->name_);
  }

  if(hueStop == 360 && hueStart == 0){
  	hueStop = 360 - (360.0/names.size());
  }
	auto nameColors = njh::getColorsForNames(names, hueStart, hueStop,
			lumStart, lumStop, satStop, satStart);


  watch.startNewLap("determine connections");
  Json::Value graphJson = graphMis.toJsonMismatchGraph(groupCutOff, mismatches, njh::color(backgroundColor), hueStart, hueStop,
			lumStart, lumStop, satStop, satStart);
  if(setUp.pars_.verbose_){
    graphMis.printAdj(std::cout);
    auto gCounts = graphMis.getGroupCounts();
    uint32_t numOfGroups = 0;
    for(const auto & g : gCounts){
    	if(g.second >=groupCutOff){
    		++numOfGroups;
    	}
    }

    std::cout << "num of groups" <<  numOfGroups << std::endl;
  	watch.logLapTimes(std::cout, true, 6, true);
  }

	std::ofstream outFile(setUp.pars_.directoryName_ + "psudoMinTree.json");
	outFile << graphJson << std::endl;
	std::string htmlOut = genHtmlStrForPsuedoMintree("psudoMinTree.json");
	std::ofstream outHtmlFile(setUp.pars_.directoryName_ + "psudoMinTree.html");
	outHtmlFile << htmlOut << std::endl;
	return 0;
}




int genExpRunner::splitFile(const njh::progutils::CmdArgs & inputCommands){
	seqSetUp setUp(inputCommands);
	setUp.processDefaultReader(true);
	setUp.finishSetUp(std::cout);

	seqInfo seq;
	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();

	while(reader.readNextRead(seq)){
		SeqIOOptions outOpts = setUp.pars_.ioOptions_;
		outOpts.out_.outFilename_ = seq.name_;
		SeqOutput::write(std::vector<seqInfo>{seq}, outOpts);
	}

	return 0;
}

} /* namespace njhseq */
