/*
 * clusteringExpRunner.cpp
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
#include "clusteringExp.hpp"


namespace njhseq {
clusteringExpRunner::clusteringExpRunner()
    : njh::progutils::ProgramRunner(
          {
					 addFunc("clusteringPairedEndReads",	clusteringPairedEndReads, false)
           },
          "clusteringExp") {}





namespace readVec {

template <>
void getMaxLength(const PairedRead& read, uint64_t & compare) {
	if (read.seqBase_.seq_.length() > compare) {
		compare = read.seqBase_.seq_.length();
	}
	if (read.mateSeqBase_.seq_.length() > compare) {
		compare = read.mateSeqBase_.seq_.length();
	}
}

template <>
void getMaxLength(const std::shared_ptr<PairedRead>& read, uint64_t & compare) {
	if (read->seqBase_.seq_.length() > compare) {
		compare = read->seqBase_.seq_.length();
	}
	if (read->mateSeqBase_.seq_.length() > compare) {
		compare = read->mateSeqBase_.seq_.length();
	}
}

}  // namespace readVec

int clusteringExpRunner::clusteringPairedEndReads(const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	setUp.processDefaultReader(VecStr{"--fastq1", "--fastq2", "--fastq1gz", "--fastq2gz"},true);
	setUp.processDirectoryOutputName(true);
	setUp.processAlignerDefualts();
	std::string qualRep = "median";
	setUp.setOption(setUp.pars_.ioOptions_.revComplMate_, "--reverseSecond", "Mate reads should be reverse complemented");
	setUp.setOption(qualRep, "--qualRep", "How to calculate the qual rep of identical clusters, options: median,max,mean");


  CollapseIterations iteratorMap;
	std::string parameters = "";

	if (setUp.setOption(parameters, "--par", "Parameters Filename")) {
		iteratorMap = setUp.processIteratorMap(parameters);
	}else{
		iteratorMap = CollapseIterations::genIlluminaDefaultParsCollapseHomopolymers(100);
	}

  std::string sortBy = "totalCount";
	setUp.processVerbose();
  setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	uint64_t maxSize = 0;
	if(setUp.pars_.verbose_){
		std::cout << "Read in from:" << std::endl;
		std::cout << setUp.pars_.ioOptions_.toJson() << std::endl;
	}
	SeqInput reader(setUp.pars_.ioOptions_);
	std::vector<PairedRead> reads = reader.readAllReads<PairedRead>();
	if(setUp.pars_.verbose_){
		std::cout << "Read in " << reads.size() << " reads" << std::endl;
	}
	readVec::getMaxLength(reads, maxSize);
	collapseIdenticalPairedClusters(reads, qualRep);
	if (setUp.pars_.verbose_) {
		std::cout << "Collapsed to  " << reads.size() << " clusters" << std::endl;
	}
	std::vector<PairedCluster> clusters;
	for(const auto & read : reads){
		clusters.emplace_back(PairedCluster(read));
	}
	// create aligner class object
	aligner alignerObj(maxSize, setUp.pars_.gapInfo_,
			setUp.pars_.scoring_, KmerMaps(), setUp.pars_.qScorePars_,
			setUp.pars_.colOpts_.alignOpts_.countEndGaps_,
			setUp.pars_.colOpts_.iTOpts_.weighHomopolyer_);
  alignerObj.processAlnInfoInput(setUp.pars_.alnInfoDirName_);
  // create collapser
  collapser collapserObj = collapser(setUp.pars_.colOpts_);

  clusters = collapserObj.runClustering(clusters, iteratorMap, alignerObj);

	if (setUp.pars_.writingOutAlnInfo_) {
		alignerObj.alnHolder_.write(setUp.pars_.outAlnInfoDirName_,
				setUp.pars_.verbose_);
	}


  renameReadNames(clusters,
			njh::replaceString(bfs::basename(setUp.pars_.ioOptions_.firstName_),
					"_R1", ""), true, true, true);
	std::ofstream outInfoFile;
	openTextFile(outInfoFile, setUp.pars_.directoryName_ + "outputInfo",
			".tab.txt", false, true);
	clusterVec::allSetFractionClusters(clusters);
	profiler::getFractionInfo(clusters, outInfoFile);

	std::ofstream outputFile1;
	openTextFile(outputFile1, setUp.pars_.directoryName_ + "output_R1",
			setUp.pars_.ioOptions_.out_.outExtention_, false, true);
	std::ofstream outputFile2;
	openTextFile(outputFile2, setUp.pars_.directoryName_ + "output_R2",
			setUp.pars_.ioOptions_.out_.outExtention_, false, true);
	auto clusDir = njh::files::makeDir(setUp.pars_.directoryName_,
			njh::files::MkdirPar("clusters", false));
	for (auto & clus : clusters) {
		clus.mateSeqBase_.name_ = clus.seqBase_.name_;
		clus.mateSeqBase_.frac_ = clus.seqBase_.frac_;
		clus.mateSeqBase_.cnt_ = clus.seqBase_.cnt_;
		if (setUp.pars_.ioOptions_.revComplMate_) {
			clus.mateSeqBase_.reverseComplementRead(false, true);
		}
		clus.outFastq(outputFile1, outputFile2);
		clus.writeOutClusters(clusDir.string(), setUp.pars_.ioOptions_);
	}
	setUp.rLog_ << alignerObj.numberOfAlingmentsDone_ << "\n";

	return 0;
}






} /* namespace njhseq */
