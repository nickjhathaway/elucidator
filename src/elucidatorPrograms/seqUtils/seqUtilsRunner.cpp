
#include "seqUtilsRunner.hpp"
#include "elucidator/simulation.h"
#include "elucidator/objects/seqContainers.h"
#include "elucidator/BamToolsUtils.h"




namespace njhseq {

seqUtilsRunner::seqUtilsRunner()
    : njh::progutils::ProgramRunner(
          {
           addFunc("mapCount", mapCount, false),
           addFunc("compareToRef", compareToRef, false),
					 addFunc("compareAllByAll", compareAllByAll, false),
           addFunc("createConsensus", createConsensus, false),
           addFunc("createDegenerativeStr", createDegenerativeStr, false),
           addFunc("checkTwoReadFiles", checkTwoReadFiles, false),
           addFunc("alignToSequence", alignToSequence, false),

          },
          "seqUtils") {}














int seqUtilsRunner::createConsensus(const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	setUp.processAlignerDefualts();
	std::string filename = "";
	bool random = false;
	uint32_t runTimes = 10;
	setUp.setOption(runTimes, "-runTimes,-run", "runTimes");
	setUp.setOption(random, "-random", "randomSeqs");
	setUp.setOption(filename, "-file", "filename");
	setUp.finishSetUp(std::cout);
	aligner alignerObj = aligner(200, setUp.pars_.gapInfo_,
			substituteMatrix(setUp.pars_.generalMatch_, setUp.pars_.generalMismatch_), KmerMaps(),
			setUp.pars_.qScorePars_,setUp.pars_.colOpts_.alignOpts_.countEndGaps_,
			setUp.pars_.colOpts_.iTOpts_.weighHomopolyer_);
	std::vector<VecStr> randomStrs(runTimes);
	VecStr strs;
	if(random){
		njh::randomGenerator gen;
		strs = simulation::evenRandStrs(20, std::vector<char>{'A', 'C', 'G', 'T'}, gen, 10 );
		std::generate(randomStrs.begin(), randomStrs.end(),[&](){ return simulation::evenRandStrs(20, std::vector<char>{'A', 'C', 'G', 'T'}, gen, 10 );} );

	}else{
		table inSeqs(filename, "whitespace", false);
		for(const auto & row : inSeqs.content_){
			strs.emplace_back(row[0]);
		}
	}
	if(random){
		for(const auto & rStrs : randomStrs){
			//printVector(rStrs, "\n");
		  std::vector<baseCluster> inClusters;
		  readObject firstRead = readObject(seqInfo("str1", rStrs[0]));
		  //firstRead.seqBase_.cnt_ = 5;
		  baseCluster mainCluster(firstRead.seqBase_);
		  baseCluster mainClusterNew(firstRead.seqBase_);
		  uint32_t clusNum = 0;
		  for(const auto & st : rStrs){
		  	++clusNum;
		  	if(clusNum == 1){
		  		continue;
		  	}
		  	inClusters.emplace_back(baseCluster(seqInfo("str" + estd::to_string(clusNum), st)));
		  }
		  for(const auto & clus : inClusters){
		  	mainCluster.addRead(clus);
		  	mainClusterNew.addRead(clus);
		  }
		  mainCluster.calculateConsensus(alignerObj, true);
		  std::cout << convertBoolToString(mainCluster.seqBase_.seq_ == mainClusterNew.seqBase_.seq_) << std::endl;;
		  if(mainCluster.seqBase_.seq_ != mainClusterNew.seqBase_.seq_){
		  	printVector(rStrs, "\n");
		  	std::cout << mainCluster.seqBase_.seq_ << std::endl;
		  	std::cout << mainClusterNew.seqBase_.seq_ << std::endl;
		  	std::cout << std::endl;
		  	std::ofstream outDisFile;
		  	openTextFile(outDisFile, "outDis", ".fasta", true, false);
		  	//uint32_t seqNum = 0;
		  	for(const auto & pos : iter::range(len(mainCluster.reads_))){
		  		outDisFile << ">" << mainCluster.reads_[pos]->seqBase_.name_ << std::endl;
		  		outDisFile << ">" << mainCluster.reads_[pos]->seqBase_.name_ +"_ref" << std::endl;
		  	}
		  	exit(1);
		  }
		}
	}else{
		printVector(strs, "\n");
	  std::vector<baseCluster> inClusters;
	  readObject firstRead = readObject(seqInfo("str1", strs[0]));
	  //firstRead.seqBase_.cnt_ = 5;
	  baseCluster mainCluster(firstRead.seqBase_);
	  uint32_t clusNum = 0;
	  for(const auto & st : strs){
	  	++clusNum;
	  	if(clusNum == 1){
	  		continue;
	  	}
	  	inClusters.emplace_back(baseCluster(seqInfo("str" + estd::to_string(clusNum), st)));
	  }
	  for(const auto & clus : inClusters){
	  	mainCluster.addRead(clus);
	  }
	  mainCluster.calculateConsensus(alignerObj, true);
	  std::cout << mainCluster.seqBase_.cnt_ << std::endl;
	  std::cout << mainCluster.seqBase_.seq_ << std::endl;
	}
  return 0;
}




int seqUtilsRunner::createDegenerativeStr(const njh::progutils::CmdArgs & inputCommands) {
  seqSetUp setUp(inputCommands);
  setUp.processDefaultReader(true);
  setUp.finishSetUp(std::cout);

  SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<readObject>();
  uint32_t minLen = std::numeric_limits<uint32_t>::max();
  for(const auto & readPos : iter::range(len(inReads))){
  	if(inReads[readPos].seqBase_.seq_.length() < minLen){
  		minLen = len(inReads[readPos].seqBase_.seq_);
  	}
  }

  readVecTrimmer::trimToMaxLength(inReads, minLen);
  VecStr dnaStrings;
  for(const auto & readPos : iter::range(len(inReads))){
  	dnaStrings.emplace_back(inReads[readPos].seqBase_.seq_);
  }
  std::string ans = seqUtil::createDegenerativeString(dnaStrings);
  std::cout << ans << std::endl;
  return 0;
}






int seqUtilsRunner::checkTwoReadFiles(const njh::progutils::CmdArgs & inputCommands) {
  // parameters
	std::string filename1 = "";
	std::string filename2 = "";
  seqUtilsSetUp setUp(inputCommands);
  setUp.setOption(filename1, "-file1", "filename1", true);
  setUp.setOption(filename2, "-file2", "filename2", true);
  setUp.processDefaultReader(false);
  setUp.finishSetUp(std::cout);

  SeqIOOptions opts1;
  opts1.firstName_ = filename1;
  opts1.inFormat_ = SeqIOOptions::getInFormat(njh::files::getExtension(filename1));
  opts1.processed_ = setUp.pars_.ioOptions_.processed_;
  SeqInput reader(opts1);
	reader.openIn();
	auto inReads1 = reader.readAllReads<readObject>();
  SeqIOOptions opts2;
  opts2.firstName_ = filename2;
  opts2.inFormat_ = SeqIOOptions::getInFormat(njh::files::getExtension(filename2));
  opts2.processed_ = setUp.pars_.ioOptions_.processed_;
  SeqInput reader2(opts2);
  reader2.openIn();
	auto inReads2 = reader2.readAllReads<readObject>();

  bool equal = readVec::checkIfReadVecsAreSame(inReads1, inReads2);
  std::cout << njh::colorBool(equal) << std::endl;
  return 0;
}


//quickLenInfo


int seqUtilsRunner::compareAllByAll(const njh::progutils::CmdArgs & inputCommands) {
	uint32_t numThreads = 1;
	bool diagonal = false;
	seqSetUp setUp(inputCommands);
	setUp.pars_.ioOptions_.lowerCaseBases_ = "upper";
  setUp.processVerbose();
  setUp.processDebug();
  setUp.processDefaultReader(true);
  setUp.pars_.ioOptions_.out_.outExtention_ = ".tab.txt";
	setUp.setOption(numThreads, "--numThreads", "Number of Threads to Use");
	setUp.setOption(diagonal,   "--diagonal",   "Just solve a global diagonal");
  setUp.pars_.gapLeft_ = "0,0";
  setUp.processAlignerDefualts();
	setUp.pars_.colOpts_.iTOpts_.weighHomopolyer_ = false;
	setUp.setOption(setUp.pars_.colOpts_.iTOpts_.weighHomopolyer_,
			"--weighHomopolymerIndels", "Do Homopolymer Weighting");
  setUp.processDirectoryOutputName(false);
  if(setUp.pars_.verbose_){
    std::cout << "go: " << setUp.pars_.gapInfo_.gapOpen_ << std::endl;
    std::cout << "ge: " << setUp.pars_.gapInfo_.gapExtend_ << std::endl;
  }
  setUp.finishSetUp(std::cout);
  //read in
  SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<seqInfo>();

	// set up aligner object
  uint64_t maxSize = 0;
  readVec::getMaxLength(inReads, maxSize);

	aligner alignerObj(maxSize, setUp.pars_.gapInfo_, setUp.pars_.scoring_,
			KmerMaps(setUp.pars_.colOpts_.kmerOpts_.kLength_),
			setUp.pars_.qScorePars_,
			setUp.pars_.colOpts_.alignOpts_.countEndGaps_,
			setUp.pars_.colOpts_.iTOpts_.weighHomopolyer_);

	std::function<void(aligner &, const seqInfo&, const seqInfo&, bool)> alignFunc;
	if (diagonal) {
		alignFunc =
				[](aligner & alnObj, const seqInfo & refSeq, const seqInfo & querySeq, bool local) {
					alnObj.alignCacheGlobalDiag(refSeq, querySeq);
				};
	} else {
		alignFunc =
				[](aligner & alnObj, const seqInfo & refSeq, const seqInfo & querySeq, bool local) {
					alnObj.alignCache(refSeq, querySeq, local);
				};
	}

	double total = std::accumulate(inReads.begin(), inReads.end(), 0.0, [](double tots, const seqInfo & seq){
		return tots + seq.cnt_;
	});
	njh::for_each(inReads, [&total](seqInfo & seq){
		seq.frac_ = seq.cnt_/total;
	});

  std::ofstream profileInfoFile;
  setUp.pars_.ioOptions_.out_.openFile(profileInfoFile);

  std::ofstream tempFile;
  if(setUp.pars_.debug_){
    openTextFile(tempFile, "tempFilealns.fasta", ".fasta", setUp.pars_.ioOptions_.out_);
  }
  profileInfoFile
      << "ReadId\tReadFraction\tOtherReadId\talnScore\tperId\t1bIndel\t2bI"
         "ndel\t>2bIndel\tlqMismatch\thqMismatch\ttotalDiffs" << std::endl;
  readVec::handelLowerCaseBases(inReads, setUp.pars_.ioOptions_.lowerCaseBases_);
	PairwisePairFactory pairFac(len(inReads));
	if (setUp.pars_.verbose_) {
		std::cout << "Read in " << inReads.size() << " reads to compare" << std::endl;
		std::cout << "Total number of compares to do: " << pairFac.totalCompares_ << std::endl;
	}
	concurrent::AlignerPool alnPool(alignerObj, numThreads);
	alnPool.inAlnDir_ = setUp.pars_.alnInfoDirName_;
	alnPool.outAlnDir_ = setUp.pars_.outAlnInfoDirName_;
	alnPool.initAligners();
	std::mutex fileMut;
	njh::ProgressBar pbar(pairFac.totalCompares_);

	auto runCompare = [&pairFac,&fileMut,&profileInfoFile, &tempFile,&setUp,&inReads, &alnPool,&pbar,&alignFunc](){

		auto threadId = estd::to_string(std::this_thread::get_id());
		PairwisePairFactory::PairwisePair pair;
		std::stringstream ssProfile;
		std::stringstream ssTempFile;
		auto currentAligner = alnPool.popAligner();
		while(pairFac.setNextPair(pair)){
			if(setUp.pars_.verbose_){
				pbar.outputProgAdd(std::cout, 1, true);
			}
			const auto & ref = inReads[pair.row_];
			const auto & input = inReads[pair.col_];
			//currentAligner->alignCache(ref, input, setUp.pars_.local_);
			alignFunc(*currentAligner, ref, input, setUp.pars_.local_);
      if(setUp.pars_.debug_){
      	ssTempFile << ">" << ref.name_ << std::endl;
      	ssTempFile << currentAligner->alignObjectA_.seqBase_.seq_ << std::endl;
      	ssTempFile << ">" << input.name_ << std::endl;
      	ssTempFile << currentAligner->alignObjectB_.seqBase_.seq_ << std::endl;
      }
      currentAligner->profilePrimerAlignment(ref, input);
      ssProfile << input.name_
					<< "\t" << input.frac_
					<< "\t" << ref.name_
					<< "\t" << currentAligner->comp_.alnScore_
					<< "\t" << currentAligner->comp_.distances_.eventBasedIdentity_
					<< "\t" << currentAligner->comp_.oneBaseIndel_
					<< "\t" << currentAligner->comp_.twoBaseIndel_
					<< "\t" << currentAligner->comp_.largeBaseIndel_
					<< "\t" << currentAligner->comp_.lqMismatches_
					<< "\t" << currentAligner->comp_.hqMismatches_
					<< "\t" << currentAligner->comp_.distances_.getNumOfEvents(true) << std::endl;
		}
		{
			std::lock_guard<std::mutex> fileLock(fileMut);
			profileInfoFile << ssProfile.str();
			if(setUp.pars_.debug_){
				tempFile << ssTempFile.str();
			}
		}
	};
	std::vector<std::thread> threads;
	for(uint32_t t = 0; t < numThreads; ++t){
		threads.emplace_back(std::thread(runCompare));
	}
	njh::concurrent::joinAllThreads(threads);

	if(setUp.pars_.verbose_){
		 setUp.logRunTime(std::cout);
	}

  return 0;
}


std::vector<uint32_t> compareToRef(const seqInfo & input,
		const std::vector<seqInfo> & inputRefs,
		aligner& alignerObj,
		bool local,
		bool eventBased){
  double bestScore = std::numeric_limits<double>::lowest();
  std::vector<uint32_t> bestMatches;
  for (const auto& refPos : iter::range(inputRefs.size())) {
  	const auto & ref = inputRefs[refPos];
    if (input.name_ == ref.name_) {
      continue;
    }
    alignerObj.alignCache(ref, input, local);
    double currentScore = 0;
    if(eventBased){
      alignerObj.profileAlignment(ref, input, false, true, false);
    	currentScore = alignerObj.comp_.distances_.eventBasedIdentity_;
    }else{
    	currentScore = alignerObj.parts_.score_;
    }

    if (currentScore == bestScore) {
    	bestMatches.push_back(refPos);
    }
    if (currentScore > bestScore) {
    	bestMatches.clear();
    	bestMatches.push_back(refPos);
      bestScore = currentScore;
    }
  }
  return bestMatches;
}












/*
void mapRead(const std::vector<uint32_t> readPositions, const std::vector<identicalCluster> & reads, std::vector<identicalCluster> & ties){
  for (const auto & readPos : readPositions) {
  	auto & read = reads[readPos];
    std::vector<size_t> bestRefs;
    std::unordered_map<std::string, std::map<uint32_t, mismatch>>
        bestRefMismatches;
    if (count % 100 == 0) {
      std::cout << count << "/" << reads.size() << '\r';
      std::cout.flush();
    }
    ++count;
    double bestScore = 0.00;
    if(checkForChimeras){
    	bool foundChimera = checkPossibleChiByRefs(read, refSeqs, outInfo, alignerObj, chiOverlap, true, setUp.pars_.weightHomopolymers_);
  		if(foundChimera){
  			chiCount += read.seqBase_.cnt_;
  			read.seqBase_.outPutFastq(chiSeqs);
  			read.remove = true;
  			continue;
  		}
    }
    for (const auto &refReadPos : iter::range(refContainers.size())) {
    	auto & refRead = refContainers[refReadPos];
      alignerObj.alignCache(refRead.seqBase_, read.seqBase_, setUp.pars_.local_);
      double currentScore = std::numeric_limits<double>::lowest();
      if(eventBased){
      	alignerObj.profilePrimerAlignment(refRead.seqBase_, read.seqBase_, setUp.pars_.weightHomopolymers_);
      	currentScore = alignerObj.comp_.distances_.eventBasedIdentity_;
      }else{
      	currentScore = alignerObj.parts_.score_;
      }
      if (currentScore == bestScore && currentScore != 0.00) {
        bestRefs.emplace_back(refReadPos);
      } else if (currentScore > bestScore) {
        bestRefs.clear();
        bestScore = currentScore;
        bestRefs.emplace_back(refReadPos);
      }
    }
    if (bestRefs.size() == 1) {
      read.remove = true;
      bestReferenceCount[refContainers[bestRefs[0]].seqBase_.name_] += read.seqBase_.cnt_;
      refContainers[bestRefs[0]].addRead(read);
    }else{
    	ties.emplace_back(read);
    }
  }
}*/


std::array<double, 100> makeQualErrorArr() {
  std::array<double, 100> arr;
  for (auto i : iter::range(100)) {
    arr[i] = std::pow(10.0, (-i / 10.0));
  }
  return arr;
}

template<typename READ, typename REF>
bool checkPossibleChiByRefsForMapCount(const READ & read,
		const std::vector<REF> & refSeqs,
		table& outInfo,
		aligner & alignerObj,
		const comparison & chiOverlap,
		bool breakAtFirst,
		bool weightHomopolymers ){
	bool foundAnExactMatch = false;
	bool foundAChimera = false;
	std::string firstChiName = "";
	std::string secondChiName = "";
	uint32_t inflectionPoint = UINT32_MAX;
	uint32_t inflectionPointPar1 = UINT32_MAX;
	uint32_t inflectionPointPar2 = UINT32_MAX;
	for(const auto & refPos : iter::range(len(refSeqs))){
		auto & ref = refSeqs[refPos];
		alignerObj.alignCache(ref.seqBase_, read.seqBase_, false);
		alignerObj.profilePrimerAlignment(ref.seqBase_, read.seqBase_);
		if(alignerObj.comp_.distances_.mismatches_.empty()){
			foundAnExactMatch = true;
			foundAChimera = false;
			return foundAChimera;
		}
	}
	//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
	for(const auto & refPos : iter::range(len(refSeqs))){

		auto & ref = refSeqs[refPos];
		alignerObj.alignCache(ref.seqBase_, read.seqBase_, false);
		alignerObj.profilePrimerAlignment(ref.seqBase_, read.seqBase_);
		if(alignerObj.comp_.distances_.mismatches_.empty()){
			foundAnExactMatch = true;
			foundAChimera = false;
		} else if(alignerObj.comp_.distances_.mismatches_.size() > 0 && !foundAnExactMatch){
			auto savedMismatches = alignerObj.comp_.distances_.mismatches_;
			//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
			//check front
			bool passFront = false;
			if(0 != savedMismatches.begin()->second.seqBasePos){
				alignerObj.profileAlignment(ref.seqBase_, read.seqBase_, false, true,
						false, 0,
						alignerObj.getAlignPosForSeqBPos(
								savedMismatches.begin()->second.seqBasePos));
				passFront = chiOverlap.passErrorProfile(alignerObj.comp_);
			}
			//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
			//check back
			bool passBack = false;
			if(savedMismatches.rbegin()->second.seqBasePos + 1 < len(read.seqBase_)){
				alignerObj.profileAlignment(ref.seqBase_, read.seqBase_, false, true,
						false,
						alignerObj.getAlignPosForSeqBPos(
								savedMismatches.rbegin()->second.seqBasePos + 1));
				passBack = chiOverlap.passErrorProfile(alignerObj.comp_);
			}


			auto firstRefAlignA = alignerObj.alignObjectA_;
			auto firstRefAlignB = alignerObj.alignObjectB_;
//			std::cout << "pass front: " << njh::colorBool(passFront) << std::endl;
//			std::cout << "pass back:  " << njh::colorBool(passBack) << std::endl;

			if(passFront){
				for(const auto & secondRefPos : iter::range(refPos + 1, len(refSeqs))){
					auto & secondRef = refSeqs[secondRefPos];
					if(ref.seqBase_.name_ == secondRef.seqBase_.name_){
						continue;
					}

					alignerObj.alignCache(secondRef.seqBase_, read.seqBase_, false);
					//check to see if from the mismatch on from the other one ref matches
					//to another ref
					//when comparing to ref's need to make sure it's not the sure sequence actually ahah
					alignerObj.profileAlignment(secondRef.seqBase_, read.seqBase_, false, true, false);
					if(alignerObj.comp_.distances_.mismatches_.size() < 1){
						continue;
					}
					//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
					alignerObj.profileAlignment(secondRef.seqBase_, read.seqBase_, false, true, false,
												alignerObj.getAlignPosForSeqBPos(savedMismatches.begin()->second.seqBasePos));
					//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
					//alignerObj.errors_.printDescription(std::cout, true);
					if(chiOverlap.passErrorProfile(alignerObj.comp_)){
						firstChiName = ref.seqBase_.name_;
						secondChiName = secondRef.seqBase_.name_;
						inflectionPoint = savedMismatches.begin()->second.seqBasePos;
						inflectionPointPar1 = savedMismatches.begin()->second.refBasePos;
						std::string refTwoPortion = alignerObj.alignObjectA_.seqBase_.seq_.substr(0, alignerObj.getAlignPosForSeqAPos(savedMismatches.begin()->second.refBasePos) + 1);
						inflectionPointPar2 = removeCharReturn(refTwoPortion, '-').size() - 1;

						foundAChimera = true;
						outInfo.content_.emplace_back(VecStr{read.seqBase_.name_, estd::to_string(read.seqBase_.frac_),
							estd::to_string(inflectionPoint),
							firstChiName, "1", estd::to_string(inflectionPointPar1),
							secondChiName, "1", estd::to_string(inflectionPointPar2)});
						if(foundAChimera && breakAtFirst){
							//VecStr{"readName", "fraction", "par1", "par1Frac", "par2", "par2Frac", "inflectionPoint"}
							break;
						}
					}
				}
			}
			if(foundAChimera && breakAtFirst){
				break;
			}
			if(passBack){
				for(const auto & secondRefPos : iter::range(refPos + 1,len(refSeqs))){
					auto & secondRef = refSeqs[secondRefPos];
					if(ref.seqBase_.name_ == secondRef.seqBase_.name_){
						continue;
					}
					alignerObj.alignCache(secondRef.seqBase_, read.seqBase_, false);
					//when comparing to ref's need to make sure it's not the sure sequence actually ahah
					alignerObj.profileAlignment(secondRef.seqBase_, read.seqBase_, false, true, false);
					if(alignerObj.comp_.distances_.mismatches_.size() < 1){
						continue;
					}
					//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
					alignerObj.profileAlignment(secondRef.seqBase_, read.seqBase_, false, true, false,
												0, alignerObj.getAlignPosForSeqBPos(savedMismatches.rbegin()->second.seqBasePos + 1));
					//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
					if(chiOverlap.passErrorProfile(alignerObj.comp_)){
						foundAChimera = true;
						firstChiName = ref.seqBase_.name_;
						secondChiName = secondRef.seqBase_.name_;
						inflectionPoint = savedMismatches.rbegin()->second.seqBasePos;
						inflectionPointPar1 = savedMismatches.rbegin()->second.refBasePos;
						std::string refTwoPortion = alignerObj.alignObjectA_.seqBase_.seq_.substr(0, alignerObj.getAlignPosForSeqAPos(savedMismatches.rbegin()->second.refBasePos) + 1);
						inflectionPointPar2 = removeCharReturn(refTwoPortion, '-').size() - 1;
						outInfo.content_.emplace_back(VecStr{read.seqBase_.name_, estd::to_string(read.seqBase_.frac_),
							estd::to_string(inflectionPoint),
							firstChiName, "1", estd::to_string(inflectionPointPar1),
							secondChiName, "1", estd::to_string(inflectionPointPar2)});
						if(foundAChimera && breakAtFirst){
							break;
						}
					}
				}
			}
		}
		//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
	}
	//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
	return foundAChimera;
}



int seqUtilsRunner::mapCount(const njh::progutils::CmdArgs & inputCommands) {
	seqUtilsSetUp setUp(inputCommands);
	bool extra = false;
	bool weightTies = true;
	uint64_t maxSize = 0;
	std::string qualRep = "median";
	bool checkForChimeras = false;
	bool eventBased = false;
	std::string expectedBarcodesFile = "";
	bool throwOutTies = false;
	uint32_t numThreads = 1;
	//setUp.setOption(numThreads, "-numThreads", "Number of Threads to Use");
	setUp.setOption(throwOutTies, "-throwOutTies", "Throw Out Ties");
	setUp.setOption(checkForChimeras, "-checkForChimeras",
			"Check For Chimeras in the reads");
	setUp.setOption(expectedBarcodesFile, "-expect", "ExpectedBarcodesFile");
	bool notieweight = false;
	setUp.setOption(notieweight, "-notieweight", "weightTies");
	weightTies = !notieweight;
	setUp.setOption(qualRep, "-qualRep", "QualRep");
	setUp.setUpMapToReferenceCount(extra);


  setUp.writeParametersFile(setUp.pars_.directoryName_ + "parsUsed.tab.txt", false,
                            false);
  // read in sequences
  SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<readObject>();
  readVec::allSetFractionByTotalCount(inReads);
  readVec::getMaxLength(inReads, maxSize);
  // make runLog and optional log files
  setUp.startARunLog(setUp.pars_.directoryName_);
  setUp.rLog_.runLogFile_ << "Read in " << readVec::getReadVectorSize(inReads)
         << " clusters" << std::endl;
  auto refSeqs = SeqInput::getReferenceSeq(
  		setUp.pars_.refIoOptions_, maxSize);
  std::string seqName = bfs::basename(setUp.pars_.ioOptions_.firstName_);


  table inTab;
  std::string expected = "";
  std::string barocde = "";
  if (expectedBarcodesFile != "") {
  	inTab = table(expectedBarcodesFile, "\t", true);
    auto row = inTab.getRows("Id", seqName);
    expected = row.getColumn("Viruses").front();
    barocde = row.getColumn("Barcode").front();
    VecStr expToks = tokenizeString(expected, "_");
    refSeqs = readVecExtractor::extractReadsWithNames(refSeqs, expToks);
  }
  //make sure the ref seqs don't contain any lowercase
  readVec::lowerCaseBasesToUpperCase(refSeqs);
  std::vector<identicalCluster> reads;
  if (setUp.pars_.ioOptions_.processed_) {
    reads = baseCluster::convertVectorToClusterVector<identicalCluster>(
        inReads);
  } else {
    reads = clusterCollapser::collapseIdenticalReads(inReads, qualRep);
  }

	KmerMaps kMaps = indexKmers(reads, setUp.pars_.colOpts_.kmerOpts_.kLength_,
			setUp.pars_.colOpts_.kmerOpts_.runCutOff_,
			setUp.pars_.colOpts_.kmerOpts_.kmersByPosition_,
			setUp.pars_.expandKmerPos_, setUp.pars_.expandKmerSize_);
	aligner alignerObj(maxSize, setUp.pars_.gapInfo_, setUp.pars_.scoring_, kMaps,
			setUp.pars_.qScorePars_, setUp.pars_.colOpts_.alignOpts_.countEndGaps_,
			setUp.pars_.colOpts_.iTOpts_.weighHomopolyer_);
  alignerObj.processAlnInfoInput(setUp.pars_.alnInfoDirName_);
  njh::sort(reads);
  readVec::allSetFractionByTotalCount(reads);
  readVec::allUpdateName(reads);
  //std::vector<cluster> refClusters =
      //baseCluster::convertVectorToClusterVector<cluster>(refSeqs);
  std::vector<refMapContainer<identicalCluster>> refContainers = refMapContainer<identicalCluster>::createContainers<identicalCluster, readObject>(refSeqs);
  if(setUp.pars_.verbose_){
    std::cout << "inReads: " << readVec::getTotalReadCount(inReads)
              << std::endl;
    std::cout << "inReadsCluster: " << reads.size() << std::endl;
    std::cout << "refReads: " << refSeqs.size() << std::endl;
  }
  std::map<std::string, double> bestReferenceCount;
  std::vector<identicalCluster> ties;
  for (const auto & refCon : refContainers) {
    bestReferenceCount[refCon.seqBase_.name_] = 0;
  }
  table chiOutInfoTab(VecStr{"seqName", "fraction", "seqPoint", "par1", "par1Frac", "par1Point", "par2", "par2Frac", "par2Point"});
  comparison chiOverlap;
	chiOverlap.oneBaseIndel_ = 2;
	chiOverlap.twoBaseIndel_ = .99;
	chiOverlap.lowKmerMismatches_ = 0;


  std::ofstream chiSeqs;
  if(checkForChimeras){
  	openTextFile(chiSeqs, setUp.pars_.directoryName_ + "chiSeqs.fastq",
  	               ".fastq", false, false);
  }
  double chiCount = 0;

  uint32_t count = 1;
  if(numThreads > 1){

  }else{

    for (const auto & readPos : iter::range(reads.size())) {
    	auto & read = reads[readPos];
      std::vector<size_t> bestRefs;
      std::unordered_map<std::string, std::map<uint32_t, mismatch>>
          bestRefMismatches;
      if(count % 100 == 0 and setUp.pars_.verbose_){
        std::cout << '\r' << count << "/" << reads.size() ;
  			std::cout.flush();
      }
			++count;
			double bestScore = 0.00;
			if (checkForChimeras) {
				bool foundChimera = checkPossibleChiByRefsForMapCount(read, refSeqs, chiOutInfoTab,
						alignerObj, chiOverlap, true,
						setUp.pars_.colOpts_.iTOpts_.weighHomopolyer_);
				if (foundChimera) {
					chiCount += read.seqBase_.cnt_;
					read.seqBase_.outPutFastq(chiSeqs);
					read.remove = true;
					continue;
				}
			}
      for (const auto &refReadPos : iter::range(refContainers.size())) {
      	auto & refRead = refContainers[refReadPos];
        alignerObj.alignCache(refRead.seqBase_, read.seqBase_, setUp.pars_.local_);
        double currentScore = std::numeric_limits<double>::lowest();
        if(eventBased){
        	alignerObj.profilePrimerAlignment(refRead.seqBase_, read.seqBase_);
        	currentScore = alignerObj.comp_.distances_.eventBasedIdentity_;
        }else{
        	currentScore = alignerObj.parts_.score_;
        }
        if (currentScore == bestScore && currentScore != 0.00) {
          bestRefs.emplace_back(refReadPos);
        } else if (currentScore > bestScore) {
          bestRefs.clear();
          bestScore = currentScore;
          bestRefs.emplace_back(refReadPos);
        }
      }
      if (bestRefs.size() == 1) {
        read.remove = true;
        bestReferenceCount[refContainers[bestRefs[0]].seqBase_.name_] += read.seqBase_.cnt_;
        refContainers[bestRefs[0]].addRead(read);
      }else{
      	ties.emplace_back(read);
      }
    }
  }
  if(setUp.pars_.verbose_){
    std::cout << std::endl;;
  }
  count = 1;
	if (!throwOutTies) {
		for (const auto &read : ties) {
      if(count % 10 == 0 and setUp.pars_.verbose_){
        std::cout << '\r' << count << "/" << reads.size() ;
  			std::cout.flush();
      }
			++count;
			double bestScore = 0.00;
			std::vector<size_t> bestRefs;
			for (const auto &refReadPos : iter::range(refContainers.size())) {
				auto & refRead = refContainers[refReadPos];
				alignerObj.alignCache(refRead.seqBase_, read.seqBase_, setUp.pars_.local_);
				double currentScore = std::numeric_limits<double>::lowest();
				if (eventBased) {
					alignerObj.profilePrimerAlignment(refRead.seqBase_, read.seqBase_);
					currentScore = alignerObj.comp_.distances_.eventBasedIdentity_;
				} else {
					currentScore = alignerObj.parts_.score_;
				}

				if (currentScore == bestScore && currentScore != 0.00) {
					bestRefs.push_back(refReadPos);
				} else if (currentScore > bestScore) {
					bestRefs.clear();
					bestScore = currentScore;
					bestRefs.push_back(refReadPos);
				}
			}
			std::map<std::string, double> currentRefCounts;
			double total = 0.00;
			for (const auto & ref : bestRefs) {
				total += bestReferenceCount[refContainers[ref].seqBase_.name_];
			}
			for (const auto &refPos : bestRefs) {
				identicalCluster tempClus = read;
				// double totalCorrector = 0;
				if (total == 0 || !weightTies) {
					tempClus.seqBase_.cnt_ = tempClus.seqBase_.cnt_ / bestRefs.size();
					// totalCorrector = bestRefs.size();
				} else {
					tempClus.seqBase_.cnt_ = tempClus.seqBase_.cnt_
							* (bestReferenceCount[refContainers[refPos].seqBase_.name_]
									/ total);
					// totalCorrector = total;
				}
				tempClus.updateName();
				tempClus.reads_.front()->seqBase_.cnt_ = tempClus.seqBase_.cnt_;
				tempClus.reads_.front()->updateName();
				refContainers[refPos].addRead(tempClus);
			}
		}
	}
	if(setUp.pars_.verbose_){
		std::cout << std::endl;
	}
  std::string clusteredDirectory =
      njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar("clusters")).string();
  for (auto &ref : refContainers) {
    if (ref.reads_.size() > 1) {
			SeqIOOptions outOpts(clusteredDirectory + ref.seqBase_.name_,
					setUp.pars_.ioOptions_.outFormat_, setUp.pars_.ioOptions_.out_);
			SeqOutput writer(outOpts);
			ref.writeReads(writer, true);
    }
  }
	if (expectedBarcodesFile != "") {
		profiler::getMapInfo(refContainers, seqName, setUp.pars_.directoryName_,
				"mapFreqInfo", expected, barocde, alignerObj);
	} else {
		profiler::getMapInfo(refContainers, seqName, setUp.pars_.directoryName_,
				"mapFreqInfo", alignerObj);
	}
	if (!ties.empty()) {
		profiler::getFractionInfo(ties, setUp.pars_.directoryName_, "tiesInfo",
				setUp.pars_.refIoOptions_.firstName_.string(), alignerObj, setUp.pars_.local_);

		SeqOutput::write(ties, setUp.pars_.directoryName_ + "ties",
				setUp.pars_.ioOptions_);
	}
	std::unique_ptr<OutputStream>chiOut;
	if(checkForChimeras){
		chiOut = std::make_unique<OutputStream>(OutOptions(njh::files::make_path(setUp.pars_.directoryName_, "chiCountInfo.tab.txt")));
		(*chiOut) << "seqName\tchimeraCount\tchiPercentage" << std::endl;
	}
  std::ofstream errorsFile;
  openTextFile(errorsFile, setUp.pars_.directoryName_ + "errorsProfile.txt", ".txt", false, false);
  std::ofstream indelFile;
  openTextFile(indelFile, setUp.pars_.directoryName_ + "indelProfile.txt", ".txt", false, false);
  errorsFile << "readName\trefName\treadCnt\tinsertions\tdeletions\tmismatchs\tidentity" << std::endl;
  indelFile << "readName\trefName\treadCnt\tindel\tsize" << std::endl;
  for(const auto & ref : refContainers){
  	for(const auto & read : ref.reads_){
			alignerObj.alignCache(ref.seqBase_, read.seqBase_, false);
			alignerObj.profileAlignment(ref.seqBase_, read.seqBase_, false, true,
					false);
  		uint32_t insertions = 0;
  		uint32_t deletions = 0;
  		for(const auto & alnGap : alignerObj.comp_.distances_.alignmentGaps_){
  			if(alnGap.second.ref_){
  				++insertions;
  				indelFile << read.seqBase_.name_  << "\t" <<  ref.seqBase_.name_ <<
  						"\t" << read.seqBase_.cnt_ << "\t" << "insertion" <<
  						"\t" << alnGap.second.size_ << std::endl;
  			}else{
  				++deletions;
  				indelFile << ref.seqBase_.name_ << "\t" << read.seqBase_.name_ <<
  				  						"\t" << read.seqBase_.cnt_ << "\t" << "deletion" <<
  				  						"\t" << alnGap.second.size_ << std::endl;
  			}
  		}
  		errorsFile << read.seqBase_.name_
  				<< "\t" << ref.seqBase_.name_
					<< "\t" << read.seqBase_.cnt_
					<< "\t" << insertions
					<< "\t" << deletions
					<< "\t" << alignerObj.comp_.distances_.mismatches_.size()
					<< "\t" << alignerObj.comp_.distances_.eventBasedIdentity_
					<< std::endl;
  	}
  }
	if(checkForChimeras){
		(*chiOut) << "seqName"
				<< "\t" << chiCount
				<< "\t" << 100 *( chiCount/readVec::getTotalReadCount(inReads)) << std::endl;

	}

  if (setUp.pars_.writingOutAlnInfo_) {
    alignerObj.alnHolder_.write(setUp.pars_.outAlnInfoDirName_);
  }
  if(setUp.pars_.verbose_){
  	std::cout << "Number of alignments done: " <<  alignerObj.numberOfAlingmentsDone_ << std::endl;
  }
  setUp.rLog_.runLogFile_ << alignerObj.numberOfAlingmentsDone_ << std::endl;
  if(setUp.pars_.verbose_){
  	std::cout << "chiCount:" << getPercentageString(chiCount, readVec::getTotalReadCount(inReads)) << std::endl;
  }
  setUp.rLog_.runLogFile_ << "chiCount:" << getPercentageString(chiCount, readVec::getTotalReadCount(inReads)) << std::endl;
  TableIOOpts outWrite(OutOptions(setUp.pars_.directoryName_ + "chiParentInfo", ".tab.txt"), "\t", chiOutInfoTab.hasHeader_);
  chiOutInfoTab.outPutContents(outWrite);
  setUp.logRunTime(std::cout);
  return 0;
}



int seqUtilsRunner::alignToSequence(const njh::progutils::CmdArgs & inputCommands) {
  seqUtilsSetUp setUp(inputCommands);
  setUp.pars_.ioOptions_.out_.outFilename_ = "out";
  setUp.pars_.ioOptions_.lowerCaseBases_ = "upper";
  setUp.pars_.gapRight_ = "0,0";
  setUp.pars_.gapLeft_ = "0,0";
  setUp.pars_.gap_ = "5,1";
  setUp.processVerbose();
  setUp.processDefaultReader(true);
  setUp.processAlignerDefualts();
  setUp.processSeq(true);
  setUp.finishSetUp(std::cout);

  SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	reader.openOut();
	seqInfo seq;
  uint64_t maxReadLength = 0;
  while(reader.readNextRead(seq)){
    readVec::getMaxLength(seq, maxReadLength);
  }
  readVec::getMaxLength(setUp.pars_.seqObj_, maxReadLength);
  // create aligner class object
  aligner alignerObj(maxReadLength, setUp.pars_.gapInfo_,setUp.pars_.scoring_);
  alignerObj.processAlnInfoInput(setUp.pars_.alnInfoDirName_, setUp.pars_.verbose_);

	reader.in_.reOpenIn();

  while(reader.readNextRead(seq)){
  	alignerObj.alignCache(setUp.pars_.seqObj_.seqBase_, seq, setUp.pars_.local_);
  	reader.write(alignerObj.alignObjectA_.seqBase_);
  	reader.write(alignerObj.alignObjectB_.seqBase_);
  }
  alignerObj.processAlnInfoOutput(setUp.pars_.outAlnInfoDirName_, setUp.pars_.verbose_);
  if(setUp.pars_.verbose_){
  	setUp.logRunTime(std::cout);
  }
  return 0;
}







}  // namespace njh
