/*
 * seqUtilsRunner_compareToRef.cpp
 *
 *  Created on: Jun 21, 2017
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

#include "seqUtilsRunner.hpp"
#include "elucidator/simulation.h"
#include "elucidator/objects/seqContainers.h"
#include "elucidator/BamToolsUtils.h"
#include "elucidator/objects/seqObjects/seqKmers.h"
#include "elucidator/utils/KmerUtils.hpp"



namespace njhseq {

int seqUtilsRunner::compareToRef(const njh::progutils::CmdArgs & inputCommands) {
	uint32_t numThreads = 1;
	double kmerCutOff = 0.8;
	bool forceMatch = false;
	bool dontSkipSameName = false;
	OutOptions outOpts(bfs::path("refComparisonInfo.tab.txt"));
  seqUtilsSetUp setUp(inputCommands);
  setUp.pars_.colOpts_.kmerOpts_.kLength_ = 5;
  setUp.setOption(dontSkipSameName, "--dontSkipSameName", "By default reads with the same name are skipped, use this flag to compare reads even if they have the same name");
  setUp.setOption(kmerCutOff, "--kmerCutOff", "kmer Cut Off");
  setUp.setOption(forceMatch, "--forceMatch", "Force finding a match");
  setUp.setOption(numThreads, "--numThreads", "Number of threads to use when comparing");
  setUp.processWritingOptions(outOpts);
  setUp.setUpCompareToRef();

  // read in the clusters

	bool setReverse = false;
	auto inputSeqs = createKmerReadVec(setUp.pars_.ioOptions_, setUp.pars_.colOpts_.kmerOpts_.kLength_, setReverse);
	auto refSeqs = createKmerReadVec(setUp.pars_.refIoOptions_, setUp.pars_.colOpts_.kmerOpts_.kLength_, setReverse);
  uint64_t maxSize = 0;
  readVec::getMaxLength(inputSeqs, maxSize);
  readVec::getMaxLength(refSeqs, maxSize);
  if(setUp.pars_.verbose_){
		std::cout << "Read in " << inputSeqs.size()
				<< " reads to map to reference sequences" << std::endl;
    std::cout << "Read in " << refSeqs.size() << " reference sequences"
              << std::endl;
  }
	// aligner object
	KmerMaps emptyMaps;
	aligner alignerObjOrig(maxSize, setUp.pars_.gapInfo_,
			setUp.pars_.scoring_, emptyMaps, setUp.pars_.qScorePars_,
			setUp.pars_.colOpts_.alignOpts_.countEndGaps_,
			setUp.pars_.colOpts_.iTOpts_.weighHomopolyer_);
//	std::cout << setUp.pars_.gapInfo_.getIdentifer() << std::endl;
//	setUp.pars_.scoring_.printScores(std::cout);

	alignerObjOrig.processAlnInfoInput(setUp.pars_.alnInfoDirName_);

	uint32_t totalCount = 0;
	for(const auto & input : inputSeqs){
		totalCount += input->seqBase_.cnt_;
	}
	njh::for_each(inputSeqs, [&totalCount](std::unique_ptr<seqWithKmerInfo>& seq) { seq->seqBase_.setFractionByCount(totalCount); });

  OutputStream profileInfoFile(outOpts);

  std::ofstream tempFile;
  openTextFile(tempFile, setUp.pars_.directoryName_ + "tempFilealns.fastq", ".fastq", true,
               false);
  uint32_t counter = 0;
//  profileInfoFile
//      << "ReadNumber\tReadId\tReadFraction\tBestRef\tscore\t1bIndel\t2bI"
//         "ndel\t>2bIndel\tlqMismatch\thqMismatch" << std::endl;
	if (setUp.pars_.colOpts_.alignOpts_.eventBased_) {
		profileInfoFile
				<< "ReadNumber\tReadId\tReadFraction\tBestRef\tscore\talnScore\thqScore\tkDist-"
				<< setUp.pars_.colOpts_.kmerOpts_.kLength_ << "\t1bIndel\t2bI"
						"ndel\t>2bIndel\tlqMismatch\thqMismatch" << std::endl;
	} else {
		profileInfoFile
				<< "ReadNumber\tReadId\tReadFraction\tBestRef\tscore\tperId\thqScore\tkDist-"
				<< setUp.pars_.colOpts_.kmerOpts_.kLength_ << "\t1bIndel\t2bI"
						"ndel\t>2bIndel\tlqMismatch\thqMismatch" << std::endl;
	}
  std::mutex outMut;
  readVec::lowerCaseBasesToUpperCase(inputSeqs);
  readVec::lowerCaseBasesToUpperCase(refSeqs);



  //set up queue
  std::vector<uint32_t> positions(inputSeqs.size());
  njh::iota<uint32_t>(positions, 0);
  //njh::reverse(positions);
  njh::concurrent::LockableQueue<uint32_t> posQueue(positions);

  //set up progress bar
  njh::ProgressBar pBar(inputSeqs.size());


  //set up aligner pool
  concurrent::AlignerPool alnPool(alignerObjOrig, numThreads);
  //alnPool.inAlnDir_ = setUp.pars_.alnInfoDirName_;
  alnPool.outAlnDir_ = setUp.pars_.outAlnInfoDirName_;
  alnPool.initAligners();


  std::function<void()> compareInput = [&dontSkipSameName,&outMut,&profileInfoFile,&tempFile,&alnPool,&setUp,&inputSeqs,&refSeqs,&counter,&pBar,&posQueue,&kmerCutOff,&forceMatch](){
  	std::vector<uint32_t> subPositions;
  	auto curAligner = alnPool.popAligner();
  	while(posQueue.getVals(subPositions, 5	)){
  		std::unordered_map<uint32_t, std::vector<uint32_t>> bestRefsForPos;
		for(const auto pos : iter::reversed(subPositions)){
				const auto & input = inputSeqs[pos];
		    double bestScore = std::numeric_limits<double>::lowest();
		    std::vector<uint32_t> bestRefs;
		    double currentKmerCutOff = kmerCutOff;
		    bool run = true;
		    while(run){
			    for (const auto refPos : iter::range(refSeqs.size())) {
			      const auto & ref = refSeqs[refPos];
						if (!dontSkipSameName && ref->seqBase_.name_ == input->seqBase_.name_) {
							if(1 == refSeqs.size()) {
								run = false;
								currentKmerCutOff = -0.05;
							}
							continue;
						}
			      if(ref->compareKmers(*input).second < currentKmerCutOff){
			       	continue;
			      }
						curAligner->alignCache(ref, input, setUp.pars_.local_);
						double currentScore = 0;
						if(setUp.pars_.colOpts_.alignOpts_.eventBased_) {
							curAligner->profileAlignment(ref, input, false, true, false);
							currentScore = curAligner->comp_.distances_.eventBasedIdentity_;
						} else {
							currentScore = curAligner->parts_.score_;
						}
						if (currentScore == bestScore) {
							bestRefs.push_back(refPos);
						}
						if (currentScore > bestScore) {
							bestRefs.clear();
							bestRefs.push_back(refPos);
							bestScore = currentScore;
						}
					}
			    bestRefsForPos[pos] = bestRefs;
			    run = false;
			    if(bestRefs.empty() && forceMatch && currentKmerCutOff > 0){
			    		run = true;
			    }
			    currentKmerCutOff -= 0.05;
		    }
			}//quickHaplotypeInformationDeeper
			{
				std::lock_guard<std::mutex> lock(outMut);
				if(setUp.pars_.verbose_){
					pBar.outputProgAdd(std::cout, subPositions.size(), true);
				}
				for(const auto & bestRefs : bestRefsForPos) {
					const auto & input = inputSeqs[bestRefs.first];
					if(bestRefs.second.empty()) {
						profileInfoFile << counter << "\t" << input->seqBase_.name_ << "\t"
						<< input->seqBase_.frac_ << "\t" << "*"
						<< "\t" << "*";
						profileInfoFile << "\t" << "*";;
						profileInfoFile
						<< "\t" << "*"
						<< "\t"
						<< "*"<< "\t"
						<< "*" << "\t"
						<< "*" << "\t"
						<< "*" << "\t"
						<< "*" << "\t"
						<< "*" << std::endl;
					} else {
				    for (const auto& bestPos : bestRefs.second) {
				    	const auto & best = refSeqs[bestPos];
				      curAligner->alignCache(best, input, setUp.pars_.local_);
				      double score = 0;
				      curAligner->profileAlignment(best, input, false, true, false);
				      if(setUp.pars_.colOpts_.alignOpts_.eventBased_){
				      		score = curAligner->comp_.distances_.eventBasedIdentity_;
				      } else {
				      		score = curAligner->parts_.score_;
				      }
				      curAligner->alignObjectA_.seqBase_.outPutFastq(tempFile);
				      curAligner->alignObjectB_.seqBase_.outPutFastq(tempFile);
				      profileInfoFile << counter << "\t" << input->seqBase_.name_ << "\t"
									<< input->seqBase_.frac_ << "\t" << best->seqBase_.name_
									<< "\t" << score;
				      if(setUp.pars_.colOpts_.alignOpts_.eventBased_){
				      		profileInfoFile << "\t" << curAligner->parts_.score_;;
				      }else{
				      		profileInfoFile << "\t" << curAligner->comp_.distances_.eventBasedIdentity_;;
				      }
				      profileInfoFile
									<< "\t" << curAligner->comp_.distances_.eventBasedIdentityHq_
									<< "\t"
									<< input->compareKmers(*best).second << "\t"
									<< curAligner->comp_.oneBaseIndel_ << "\t"
									<< curAligner->comp_.twoBaseIndel_<< "\t"
									<< curAligner->comp_.largeBaseIndel_ << "\t"
									<< curAligner->comp_.lqMismatches_ << "\t"
									<< curAligner->comp_.hqMismatches_ << std::endl;
				    }
					}
			    ++counter;
				}
			}
  		}
	};
  njh::concurrent::runVoidFunctionThreaded(compareInput, numThreads);

	if (setUp.pars_.verbose_) {
		setUp.logRunTime(std::cout);
	}

	return 0;
}




}  // namespace njhseq

