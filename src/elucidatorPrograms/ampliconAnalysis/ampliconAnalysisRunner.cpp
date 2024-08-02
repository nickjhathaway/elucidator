#include "ampliconAnalysisRunner.hpp"
#include "elucidator/objects/dataContainers.h"
#include "elucidator/objects/seqObjects/seqKmers.h"
#include "elucidator/objects/seqContainers/refMapContainer.hpp"

#include <njhseq/IO/SeqIO/SeqIO.hpp>
#include <njhseq/objects/seqObjects/Clusters/clusterUtils.hpp>
#include <njhseq/objects/seqObjects/Clusters/cluster.hpp>
#include <njhseq/helpers/profiler.hpp>
#include <njhseq/helpers/clusterCollapser.hpp>
#include <njhseq/objects/collapseObjects/collapser.hpp>

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

namespace njhseq {

ampliconAnalysisRunner::ampliconAnalysisRunner()
    : njh::progutils::ProgramRunner({
                      addFunc("collapseTandems", collapseTandems, false),
                      addFunc("markChimeras", markChimeras, false),
                      addFunc("greedyCluster", greedyCluster, false),
                      addFunc("singleLinkageClusteringOnPerId", singleLinkageClusteringOnPerId, false),
                      addFunc("processRawExtractByKmerPathWeaverResults", processRawExtractByKmerPathWeaverResults, false),
                      addFunc("extractedTarAmpInfoFileToJson", extractedTarAmpInfoFileToJson, false),
                      addFunc("finalClustersFileToJson", finalClustersFileToJson, false),
                      addFunc("specimenInfoFileToJson", specimenInfoFileToJson, false),
                      addFunc("combingAllIntoPMOJson", combingAllIntoPMOJson, false),
                      addFunc("experimentInfoFileToJson", experimentInfoFileToJson, false),
                      addFunc("demultiplexedExperimentSampleFileToJson", demultiplexedExperimentSampleFileToJson, false),
                      addFunc("determinePossibleMaskFromSeqs", determinePossibleMaskFromSeqs, false),
                      addFunc("maskRegionBasedOnRefSubRegions", maskRegionBasedOnRefSubRegions, false),
										 },//
                    "ampliconAnalysis") {}

//


int ampliconAnalysisRunner::singleLinkageClusteringOnPerId(
		const njh::progutils::CmdArgs & inputCommands) {
	ampliconAnalysisSetUp setUp(inputCommands);
	double perIdCutOff = 0.97;
	uint32_t numThreads = 2;
	setUp.setOption(numThreads, "--threads,-t", "Number of Threads");
	setUp.setOption(perIdCutOff, "--perIdCutOff", "Per Id Cut Off");
	setUp.processDefaultReader(true);
	//setUp.processDirectoryOutputName(true);
	setUp.processAlignerDefualts();
	setUp.finishSetUp(std::cout);

	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto reads = reader.readAllReads<readObject>();
	uint64_t maxSize = 0;
	readVec::getMaxLength(reads, maxSize);
	aligner alignerObj(maxSize, setUp.pars_.gapInfo_, setUp.pars_.scoring_, KmerMaps(),
			setUp.pars_.qScorePars_,
			setUp.pars_.colOpts_.alignOpts_.countEndGaps_,
			setUp.pars_.colOpts_.iTOpts_.weighHomopolyer_);
	alignerObj.processAlnInfoInput(setUp.pars_.alnInfoDirName_);

	std::function<double(const readObject &, const readObject &, aligner)> getPerId =
			[](const readObject & read1,
					const readObject & read2, aligner alignerObj) {
		alignerObj.alignCache(read1, read2, false);
		alignerObj.profilePrimerAlignment(read1, read2);
		return alignerObj.comp_.distances_.eventBasedIdentity_;
	};
	auto misDistances = getDistanceCopy(reads, numThreads, getPerId,
			alignerObj);
	readDistGraph<double> graphMis(misDistances, reads);
	graphMis.turnOffEdgesUnder(perIdCutOff);
	graphMis.determineGroups();
	auto groupValues = graphMis.getGroupValues();
	for(const auto & g : groupValues){
		std::cout << g.first << std::endl;
		std::cout << "\t";
		for(const auto & r : g.second){
			std::cout << r->name_ << " ";
		}
		std::cout << std::endl;
	}
	return 0;
}

int ampliconAnalysisRunner::greedyCluster(const njh::progutils::CmdArgs & inputCommands) {
	uint32_t indelSize = 16;
	uint32_t runTimes = 2;
	bool findBest = false;
  bool additionalOut = false;
  std::string additionalOutLocationFile = "";
	seqSetUp setUp(inputCommands);
	setUp.pars_.ioOptions_.out_.outFilename_ = "output";
	setUp.processDefaultReader(true);
	setUp.setOption(indelSize, "--indelSize", "indelSize");
	setUp.setOption(runTimes, "--runTimes", "runTimes");
	additionalOut = setUp.setOption(additionalOutLocationFile, "--additionalOut",
      "AdditionalOutFilename", true);
	setUp.processRefFilename(false);
	setUp.processDirectoryOutputName(true);
	setUp.processAlignerDefualts();
	setUp.finishSetUp(std::cout);
	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto reads = reader.readAllReads<readObject>();;
	auto inReads = reader.readAllReads<readObject>();
	setUp.startARunLog(setUp.pars_.directoryName_);
	std::vector<cluster> clusters = readVec::convertVec<readObject, cluster>(inReads);
	njh::sort(clusters);
	uint64_t maxSize = 0;
	readVec::getMaxLength(clusters, maxSize);
	aligner alignerObj(maxSize, setUp.pars_.gapInfo_, setUp.pars_.scoring_, false);
	alignerObj.processAlnInfoInput(setUp.pars_.alnInfoDirName_);
	std::vector<uint64_t> clusPositions(clusters.size());
	njh::iota<uint64_t>(clusPositions, 0);
	for(uint32_t run = 0; run < runTimes; ++run){
		uint32_t count = 1;
		uint32_t clusCount = readVec::getReadVectorSize(clusters, false);
		for(const auto & readPos : iter::reversed(clusPositions)){
			if(clusters[readPos].remove){
				continue;
			}
			std::cout << "Currently on " << count << " of " << clusCount << std::endl;
			++count;
			uint64_t bestPosition = std::numeric_limits<uint64_t>::max();
			double bestScore = std::numeric_limits<double>::lowest();

			for(const auto searchPos : iter::range(clusters.size())){
				if(readPos == searchPos){
					break;
				}
				if(clusters[searchPos].remove){
					continue;
				}
				alignerObj.alignCache(clusters[searchPos].seqBase_,
						clusters[readPos].seqBase_, false);
				alignerObj.profileAlignment(clusters[searchPos].seqBase_,
						clusters[readPos].seqBase_, false, true, false);
				bool add = true;

				for (const auto & g : alignerObj.comp_.distances_.alignmentGaps_) {
					if (g.second.size_ >= indelSize){
						add = false;
						break;
					}
				}
				if(add){
					if(findBest){
						if(alignerObj.parts_.score_ > bestScore){
							bestPosition = searchPos;
							bestScore = alignerObj.parts_.score_;
						}
					}else{
						clusters[searchPos].addRead(clusters[readPos]);
						clusters[readPos].remove = true;
						break;
					}
				}
			}
			if(findBest){
				if(bestPosition != std::numeric_limits<uint64_t>::max()){
					clusters[bestPosition].addRead(clusters[readPos]);
					clusters[readPos].remove = true;
				}
			}
		}
		clusterVec::allCalculateConsensus(clusters, alignerObj, true);
		njh::sort(clusters);
    readVec::allUpdateName(clusters);
	}

	clusters = readVecSplitter::splitVectorOnRemove(clusters).first;
  std::string seqName = bfs::basename(setUp.pars_.ioOptions_.firstName_);
  renameReadNames(clusters, seqName, true, true);
  if (setUp.pars_.refIoOptions_.firstName_ == "") {
    profiler::getFractionInfoCluster(clusters, setUp.pars_.directoryName_, "outputInfo");
  } else {
    profiler::getFractionInfoCluster(clusters, setUp.pars_.directoryName_, "outputInfo",
    		setUp.pars_.refIoOptions_.firstName_.string(), alignerObj, setUp.pars_.local_);
  }
  bool containsCompReads = false;
  int compCount = 0;
  readVec::getCountOfReadNameContaining(inReads, "_Comp", compCount);
  if (compCount > 0) {
    containsCompReads = true;
  }
  std::ofstream compStats;
  if (containsCompReads) {
    openTextFile(compStats, setUp.pars_.directoryName_ + "compStats.tab.txt",
                 ".txt", false, false);
    compStats << "cluster\tcompAmount" << std::endl;
    for(const auto & clus : clusters){
      int currentCompAmount = 0;
      readVec::getCountOfReadNameContaining(clus.reads_, "_Comp",
                                            currentCompAmount);
      compStats << clus.seqBase_.name_ << "\t"
                << getPercentageString(currentCompAmount, clus.seqBase_.cnt_)
                << std::endl;
    }
  }
  if (setUp.pars_.writingOutAlnInfo_) {
  	//njh::scopedStopWatch writingTimer("Writing aln infos", true);
    alignerObj.alnHolder_.write(setUp.pars_.outAlnInfoDirName_);
  }
  SeqOutput::write(clusters, SeqIOOptions(setUp.pars_.directoryName_ + setUp.pars_.ioOptions_.out_.outFilename_.string(),
  		setUp.pars_.ioOptions_.outFormat_,setUp.pars_.ioOptions_.out_));
  std::string clusterDirectoryName =
      njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar("clusters")).string();
  clusterVec::allWriteClustersInDir(clusters, clusterDirectoryName, setUp.pars_.ioOptions_);
  if (additionalOut) {
    std::string additionalOutDir = findAdditonalOutLocation(
        additionalOutLocationFile, setUp.pars_.ioOptions_.firstName_.string());
    SeqOutput::write(clusters, SeqIOOptions(additionalOutDir + setUp.pars_.ioOptions_.out_.outFilename_.string(),
    		setUp.pars_.ioOptions_.outFormat_,setUp.pars_.ioOptions_.out_));
  }
  setUp.logRunTime(std::cout);
	return 0;
}

int ampliconAnalysisRunner::collapseTandems(const njh::progutils::CmdArgs & inputCommands) {
	clusterCollapser::collapseTandemsPars tpars;

  bool extra = false;
  bool checkingAgainstReference = false;
  bool additionalOut = false;
  std::string additionalOutName = "";
  uint64_t maxSize = 0;
  std::string additionalOutLocationFile = "";
  ampliconAnalysisSetUp setUp(inputCommands);
  setUp.setOption(additionalOutName, "--additionalOutName,-addName", "additionalOutName");
  setUp.setUpCollapseTandems(tpars, extra, additionalOut,
                             additionalOutLocationFile);
  // read in sequences
  SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto reads = reader.readAllReads<readObject>();;
  auto inReads = reader.readAllReads<readObject>();
  processRunCutoff(setUp.pars_.colOpts_.kmerOpts_.runCutOff_, setUp.pars_.colOpts_.kmerOpts_.runCutOffString_,
                   readVec::getTotalReadCount(inReads));
  std::vector<cluster> processedReads =
      baseCluster::convertVectorToClusterVector<cluster>(inReads);
  readVec::getMaxLength(processedReads, maxSize);
  // make runLog and optional log files
  setUp.startARunLog(setUp.pars_.directoryName_);
  setUp.rLog_ << "Read in " << processedReads.size() << " clusters\n";

  // create kmer map
  checkingAgainstReference = setUp.pars_.refIoOptions_.firstName_ != "";
  std::vector<readObject> refSeqs;
  if (checkingAgainstReference) {
    refSeqs = SeqInput::getReferenceSeq(
    		setUp.pars_.refIoOptions_, maxSize);
  }
  // construct aligner object

  // create aligner class object
	//auto scoringMatrixMap = createDegenScoreMatrix(1,-1);
	KmerMaps kMaps = indexKmers(processedReads, setUp.pars_.colOpts_.kmerOpts_.kLength_, setUp.pars_.colOpts_.kmerOpts_.runCutOff_,
			setUp.pars_.colOpts_.kmerOpts_.kmersByPosition_, setUp.pars_.expandKmerPos_, setUp.pars_.expandKmerSize_);
	gapScoringParameters gapPars(setUp.pars_.gapInfo_);
	aligner alignerObj(maxSize, gapPars, setUp.pars_.scoring_, kMaps,
			setUp.pars_.qScorePars_, setUp.pars_.colOpts_.alignOpts_.countEndGaps_,
			setUp.pars_.colOpts_.iTOpts_.weighHomopolyer_);
  clusterCollapser::collapseTandems(
      processedReads, alignerObj, tpars);
  processedReads = readVecSplitter::splitVectorOnRemove(processedReads).first;
  readVec::allUpdateName(processedReads);
  readVecSorter::sort(processedReads);
  SeqOutput::write(processedReads,
               setUp.pars_.directoryName_ + setUp.pars_.ioOptions_.out_.outFilename_.string(),
               setUp.pars_.ioOptions_);
  if (additionalOut) {
  	std::string additionalOutDir = "";
  	if(additionalOutName == ""){
  		additionalOutDir = findAdditonalOutLocation(
  		        additionalOutLocationFile, setUp.pars_.ioOptions_.firstName_.string());
  	}else{
  		additionalOutDir = findAdditonalOutLocation(
  		        additionalOutLocationFile, additionalOutName);
  	}
  	SeqOutput::write(processedReads,
                 SeqIOOptions(additionalOutDir + setUp.pars_.ioOptions_.out_.outFilename_.string(),
                		 setUp.pars_.ioOptions_.outFormat_,setUp.pars_.ioOptions_.out_));
  }
  if (checkingAgainstReference) {
    //alignerObj.CountEndGaps() = false;
    profiler::getFractionInfoCluster(processedReads, setUp.pars_.directoryName_,
                              "tandemCollapsedInfo", setUp.pars_.refIoOptions_.firstName_.string(),
                              alignerObj, setUp.pars_.local_);
  } else {
    profiler::getFractionInfoCluster(processedReads, setUp.pars_.directoryName_,
                              "tandemCollapsedInfo");
  }
  if (extra) {
    std::string clusDirName = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar("clusters")).string();
    clusterVec::allWriteClustersInDir(processedReads, clusDirName, setUp.pars_.ioOptions_);
    //clusterVec::allWriteOutAlignments(processedReads, setUp.pars_.directoryName_,
    //                                  alignerObj);
  }
  std::cout << "Read in " << inReads.size() << " clusters" << std::endl;
  std::cout << "Collapsed down to " << processedReads.size() << " clusters"
            << std::endl;
  return 0;
}

int ampliconAnalysisRunner::markChimeras(
		const njh::progutils::CmdArgs & inputCommands) {
	uint64_t maxSize = 0;
	bool checkingAgainstReference = false;

	ampliconAnalysisSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setUpMarkChimeras();
	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<readObject>();
	processRunCutoff(setUp.pars_.colOpts_.kmerOpts_.runCutOff_,
			setUp.pars_.colOpts_.kmerOpts_.runCutOffString_,
			readVec::getTotalReadCount(inReads));
	std::vector<cluster> processedReads =
			baseCluster::convertVectorToClusterVector<cluster>(inReads);
	readVec::getMaxLength(processedReads, maxSize);
	std::vector<readObject> refSeqs;
	if (setUp.pars_.refIoOptions_.firstName_ != "") {
		checkingAgainstReference = true;
		refSeqs = SeqInput::getReferenceSeq(setUp.pars_.refIoOptions_, maxSize);
	}
	KmerMaps kMaps = indexKmers(inReads, setUp.pars_.colOpts_.kmerOpts_.kLength_,
			setUp.pars_.colOpts_.kmerOpts_.runCutOff_,
			setUp.pars_.colOpts_.kmerOpts_.kmersByPosition_,
			setUp.pars_.expandKmerPos_, setUp.pars_.expandKmerSize_);
	gapScoringParameters gapPars(setUp.pars_.gapInfo_);
	aligner alignerObj = aligner(maxSize, gapPars, setUp.pars_.scoring_, kMaps,
			setUp.pars_.qScorePars_, setUp.pars_.colOpts_.alignOpts_.countEndGaps_,
			setUp.pars_.colOpts_.iTOpts_.weighHomopolyer_);


	readVecSorter::sort(processedReads);
	readVec::allSetFractionByTotalCount(processedReads);


	// mark the chimeras
	setUp.pars_.chiOpts_.chiOverlap_.largeBaseIndel_ = .99;
	collapser collapserObj(setUp.pars_.colOpts_);
	std::ofstream chimerasInfoFile;
	openTextFile(chimerasInfoFile,
			setUp.pars_.directoryName_ + "chimeraNumberInfo.txt", ".txt",
			setUp.pars_.ioOptions_.out_.overWriteFile_,
			setUp.pars_.ioOptions_.out_.exitOnFailureToWrite_);
	chimerasInfoFile << "#chimericClusters\t#chimericReads" << std::endl;
	setUp.pars_.chiOpts_.chiOverlap_.largeBaseIndel_ = .99;
	auto chiInfoTab = collapserObj.markChimeras(processedReads, alignerObj,
			setUp.pars_.chiOpts_);
	chiInfoTab.outPutContents(
			TableIOOpts(
					OutOptions(setUp.pars_.directoryName_ + "chiParentsInfo.txt",
							".txt", "tab", setUp.pars_.ioOptions_.out_.append_, setUp.pars_.ioOptions_.out_.overWriteFile_, setUp.pars_.ioOptions_.out_.exitOnFailureToWrite_), "\t", true));

	int clusterCount = 0;
	int readCount = 0;
	readVec::getCountOfReadNameContaining(processedReads, "CHI", clusterCount);
	readVec::getReadCountOfReadNameContaining(processedReads, "CHI", readCount);
	chimerasInfoFile << getPercentageString(clusterCount, processedReads.size())
			<< "\t"
			<< getPercentageString(readCount, readVec::getTotalReadCount(processedReads))
			<< std::endl;


	SeqOutput::write(processedReads,
			setUp.pars_.directoryName_ + setUp.pars_.ioOptions_.out_.outFilename_.string(),
			setUp.pars_.ioOptions_);
	if(bfs::exists(setUp.pars_.directoryName_ +
				setUp.pars_.ioOptions_.out_.outFilename_.string() + "Info.tab.txt") && setUp.pars_.ioOptions_.out_.overWriteFile_){
		bfs::remove(setUp.pars_.directoryName_ +
				setUp.pars_.ioOptions_.out_.outFilename_.string() + "Info.tab.txt");
	}
	if (checkingAgainstReference) {
		//alignerObj.CountEndGaps() = false;
		profiler::getFractionInfo(processedReads, setUp.pars_.directoryName_,
				setUp.pars_.ioOptions_.out_.outFilename_.string() + "Info",
				setUp.pars_.refIoOptions_.firstName_.string(), alignerObj, setUp.pars_.local_);
	} else {
		profiler::getFractionInfo(processedReads, setUp.pars_.directoryName_,
				setUp.pars_.ioOptions_.out_.outFilename_.string() + "Info");
	}
	if (setUp.pars_.verbose_) {
		std::cout << "Marked " << clusterCount << " as chimeric" << std::endl;
		setUp.logRunTime(std::cout);
	}
	return 0;
}

}  // namespace njh
