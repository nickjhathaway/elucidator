/*
 * kmerExpRunner.cpp
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

#include "kmerExp.hpp"
#include "elucidator/utils/KmerUtils.hpp"
#include "elucidator/objects/dataContainers.h"
#include "elucidator/simulation.h"
#include "elucidator/objects/seqObjects/seqKmers.h"

#include <njhseq/objects/seqObjects/seqKmers.h>

#include <njhseq/concurrency/AllByAllPairFactory.hpp>


namespace njhseq {
kmerExpRunner::kmerExpRunner()
    : njh::progutils::ProgramRunner(
          {

					 addFunc("profileScanningKmerDist", profileScanningKmerDist, false),
					 addFunc("profileKmerAccerlation", profileKmerAccerlation, false),
					 addFunc("getNewScanKmerDist", getNewScanKmerDist, false),
					 addFunc("readingDistanceCheck", readingDistanceCheck, false),
					 addFunc("getKmerDistAgainstRef", getKmerDistAgainstRef, false),
					 addFunc("writingDistanceCheck", writingDistanceCheck, false),
					 addFunc("writeKmerAccerlation",writeKmerAccerlation, false),
					 addFunc("kDistVsNucDist",kDistVsNucDist, false),
					 addFunc("clostestKmerDist",clostestKmerDist, false),
					 addFunc("getKmerDistStatsMultiple",getKmerDistStatsMultiple, false),
					 addFunc("getKmerDistStats",getKmerDistStats, false),
					 addFunc("scaningKmerDist",scaningKmerDist, false),
					 addFunc("getKmerDistTwoSeqs",getKmerDistTwoSeqs, false),
					 addFunc("scoveViaKmers",scoveViaKmers, false),
					 addFunc("kmerRevVsForDist",kmerRevVsForDist, false),
					 addFunc("profileLargeKmerIndex",profileLargeKmerIndex, false),
					 addFunc("profileSharedKmerBlocks",profileSharedKmerBlocks, false),
					 addFunc("pidVsKmers", pidVsKmers, false),
					 addFunc("randomSampleKmerCompare", randomSampleKmerCompare, false),
					 addFunc("findingMinimumKLenForNoRedundantKmers", findingMinimumKLenForNoRedundantKmers, false),
					 addFunc("microsatsKmerSearch", microsatsKmerSearch, false),
					 addFunc("genomeKmerCompare", genomeKmerCompare, false),
					 addFunc("kmerSearch", kmerSearch, false),
					 addFunc("convertKmerSearchToBinaryMatrix", convertKmerSearchToBinaryMatrix, false),
					 addFunc("generateCountsTable", generateCountsTable, false),
					 addFunc("getBestKmerDist", getBestKmerDist, false),
					 addFunc("writeKmerSimDistanceMatrix", writeKmerSimDistanceMatrix, false),
					 addFunc("findUniqKmersBetweenSeqs", findUniqKmersBetweenSeqs, false),
					 addFunc("kmerPositionQualCounts", kmerPositionQualCounts, false),
					 addFunc("kmerTestingGround", kmerTestingGround, false),
					 addFunc("pairwiseWithinComparisonOfUniqueKmers", pairwiseWithinComparisonOfUniqueKmers, false),
					 addFunc("getWithinGenomeUniqueKmers", getWithinGenomeUniqueKmers, false),
					 addFunc("getUniqKmerBlocksOnGenomeAgainstRef", getUniqKmerBlocksOnGenomeAgainstRef, false),
					 addFunc("findPositionsOfUniqueKmersInEachOther", findPositionsOfUniqueKmersInEachOther, false),
					 addFunc("kmerCompareTwoSetsOfContigs", kmerCompareTwoSetsOfContigs, false),
					 addFunc("kmerConnectionGraph", kmerConnectionGraph, false),
					 addFunc("countKmers", countKmers, false),
					 addFunc("getKmerSharedBlocksBetweenGenomes", getKmerSharedBlocksBetweenGenomes, false),
					 addFunc("allByAllComparisonOfUniqueKmers", allByAllComparisonOfUniqueKmers, false),
					 addFunc("getKmerSetDistBetween", getKmerSetDistBetween, false),
					 addFunc("getAvgKmerPosition", getAvgKmerPosition, false),
					 addFunc("getAvgKmerPositionPerSeq", getAvgKmerPositionPerSeq, false),
					 addFunc("filterPerSeqAvgKmerPosition", filterPerSeqAvgKmerPosition, false),
					 addFunc("findUniqKmersBetweenSeqSets", findUniqKmersBetweenSeqSets, false),
					 addFunc("findUniqKmersBetweenSeqSetsMulti", findUniqKmersBetweenSeqSetsMulti, false),
					 //addFunc("findUniqKmersBetweenSeqSetsMultiDev", findUniqKmersBetweenSeqSetsMultiDev, false),
					 addFunc("countingUniqKmersFromSets", countingUniqKmersFromSets, false),
					 addFunc("extractByCountingUniqKmersFromSets", extractByCountingUniqKmersFromSets, false),
					 addFunc("countingUniqKmersFromSetsPerRead", countingUniqKmersFromSetsPerRead, false),


					 addFunc("filterUniqueKmerSetForEntropy", filterUniqueKmerSetForEntropy, false),
					 addFunc("testingSimpleKmerHasher", testingSimpleKmerHasher, false),
					 addFunc("countingUniqKmersFromSetsInUnmappedAlns", countingUniqKmersFromSetsInUnmappedAlns, false),
					 addFunc("chromVsChromUniqueComp", chromVsChromUniqueComp, false),
					 addFunc("getKmerDetailedKmerDistAgainstRef", getKmerDetailedKmerDistAgainstRef, false),
					 addFunc("findingKmerEnrichment", findingKmerEnrichment, false),
					 addFunc("simpleHashKmer", simpleHashKmer, false),
					 addFunc("findKmersInSets", findKmersInSets, false),
					 //
           },
          "kmerExp") {}

int kmerExpRunner::findUniqKmersBetweenSeqs(const njh::progutils::CmdArgs & inputCommands) {
	seqInfo seq1("seq1");
	seqInfo seq2("seq2");
	OutOptions outOpts(bfs::path(""), ".txt");
	uint32_t klen = 40;
	bool sortByPosition = false;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processSeq(seq1, "--seq1", "seq1", true);
	setUp.processSeq(seq2, "--seq2", "seq2", true);
	setUp.processWritingOptions(outOpts);
	setUp.setOption(klen, "--klen", "Kmer Length");
	setUp.setOption(sortByPosition, "--sortByPosition", "Sort By Position");

	setUp.finishSetUp(std::cout);

	kmerInfo seq1Info(seq1.seq_, klen, false);
	kmerInfo seq2Info(seq2.seq_, klen, false);

	OutputStream out(outOpts);

	out << "seq\tuniqKmer\tposition" << std::endl;
	if(sortByPosition){
		{
			std::map<uint32_t, std::vector<std::string>> kmersByPositon;
			for (const auto & k : seq1Info.kmers_) {
				if (!njh::in(k.first, seq2Info.kmers_)) {
					for (const auto & pos : k.second.positions_) {
						kmersByPositon[pos].emplace_back(k.first);
					}
				}
			}
			for (const auto & pos : kmersByPositon) {
				for (const auto & k : pos.second) {
					out << seq1.name_ << "\t" << k << "\t" << pos.first << std::endl;
				}
			}
		}
		{
			std::map<uint32_t, std::vector<std::string>> kmersByPositon;
			for (const auto & k : seq2Info.kmers_) {
				if (!njh::in(k.first, seq1Info.kmers_)) {
					for (const auto & pos : k.second.positions_) {
						kmersByPositon[pos].emplace_back(k.first);
					}
				}
			}
			for (const auto & pos : kmersByPositon) {
				for (const auto & k : pos.second) {
					out << seq2.name_ << "\t" << k << "\t" << pos.first << std::endl;
				}
			}
		}
	}else{
		for(const auto & k : seq1Info.kmers_){
			if(!njh::in(k.first, seq2Info.kmers_)){
				for(const auto & pos : k.second.positions_){
					out << seq1.name_
							<< "\t" << k.first
							<< "\t"<< pos << std::endl;
				}
			}
		}
		for(const auto & k : seq2Info.kmers_){
			if(!njh::in(k.first, seq1Info.kmers_)){
				for(const auto & pos : k.second.positions_){
					out << seq2.name_
							<< "\t" << k.first
							<< "\t"<< pos << std::endl;
				}
			}
		}
	}


	return 0;

}



int kmerExpRunner::profileLargeKmerIndex(const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	uint32_t seqLen = 170000;
	uint32_t kmerLen = 18;
	setUp.setOption(seqLen, "--seqLen,-l", "Sequence Length");
	setUp.setOption(seqLen, "--kmerLen,-k", "Sequence Length");
	setUp.finishSetUp(std::cout);
	njh::randomGenerator gen;
	njh::stopWatch watch;
	watch.setLapName("gen");
	auto str = simulation::evenRandStr(seqLen, { 'A', 'C', 'G', 'T' }, gen);
	watch.startNewLap("kmer index");
	watch.logLapTimes(std::cout, true, 6, false);
	kmerInfo kInfo(str, kmerLen, false);
	watch.logLapTimes(std::cout, true, 6, true);
	return 0;
}




int kmerExpRunner::kmerRevVsForDist(const njh::progutils::CmdArgs & inputCommands){
	seqSetUp setUp(inputCommands);
	setUp.pars_.colOpts_.kmerOpts_.kLength_ = 10;
	uint32_t numThreads = 1;
	setUp.setOption(numThreads, "--threads,-t", "Number of Threads To use");
	setUp.processDefaultReader(true);
	setUp.processKmerLenOptions();
	setUp.processDebug();
	setUp.processVerbose();
  setUp.finishSetUp(std::cout);
  auto reads = createKmerReadVec(setUp.pars_.ioOptions_,setUp.pars_.colOpts_.kmerOpts_.kLength_, true);
	std::function<
			std::pair<double, double>(const std::unique_ptr<seqWithKmerInfo> &,
					const std::unique_ptr<seqWithKmerInfo> &)> disFun =
			[](const std::unique_ptr<seqWithKmerInfo> & read1,
					const std::unique_ptr<seqWithKmerInfo> & read2) {
				auto dist = read1->compareKmers(*read2);
				auto distRev = read1->compareKmersRevComp(*read2);
				return std::pair<double,double> {dist.second, distRev.second};
			};
	auto distances = getDistance(reads, numThreads, disFun);
	for (const auto pos : iter::range(reads.size())) {
		for (const auto subPos : iter::range(pos)) {
			std::cout << distances[pos][subPos].first << ","
					<< distances[pos][subPos].second << " ";
		}
		std::cout << '\n';
	}
  return 0;
}

//scaningKmerDist
int kmerExpRunner::randomSampleKmerCompare(const njh::progutils::CmdArgs & inputCommands) {
  uint32_t kmerStart = 2;
  uint32_t kmerStop = 5;
  uint32_t numThreads = 2;
  uint32_t sampleAmount = 600;
  uint64_t seed = std::numeric_limits<uint64_t>::max();
	seqSetUp setUp(inputCommands);
	setUp.pars_.ioOptions_.out_.outFilename_ = "pidVsKSim.tab.txt";
	setUp.pars_.ioOptions_.lowerCaseBases_ = "upper";
	setUp.processDefaultReader(true);
	setUp.pars_.ioOptions_.out_.outExtention_ = ".tab.txt";
	setUp.pars_.ioOptions_.out_.outFileFormat_ = "tab";
	setUp.processAlignerDefualts();
  setUp.setOption(numThreads, "--numThreads", "Number of Threads to Use");
  setUp.setOption(kmerStart, "--kmerLenStart", "Length for kmers to start at");
  setUp.setOption(kmerStop, "--kmerLenStop", "Length for kmers to stop at");
  setUp.setOption(sampleAmount, "--sampleAmount", "Sample Amount for the two groups");
  setUp.setOption(seed, "--seed", "seed for randomness");
  setUp.finishSetUp(std::cout);
  setUp.startARunLog(setUp.pars_.directoryName_);
	if (kmerStop < kmerStart) {
		std::stringstream ss;
		ss << "Kmer stop must be greater than kmerStart" << std::endl;
		ss << "KmerStart: " << kmerStart << std::endl;
		ss << "KmerStop: " << kmerStop << std::endl;
		throw std::runtime_error { ss.str() };
	}
	if (setUp.pars_.ioOptions_.out_.outExists()
			&& !setUp.pars_.ioOptions_.out_.overWriteFile_) {
		std::stringstream ss;
		ss << "File: "
				<< setUp.pars_.ioOptions_.out_.outName()
				<< " already exists, use -overWrite to over write it" << std::endl;
		throw std::runtime_error { ss.str() };
	}
  SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<readObject>();
	if(inReads.size() * 2 < sampleAmount){
		std::stringstream ss;
		ss << "Sample Amount is too great" << std::endl;
		ss << "sampleAmount: " << sampleAmount << std::endl;
		ss << "inReads.size() * 2: " << inReads.size() * 2 << std::endl;
		throw std::runtime_error { ss.str() };
	}
	std::vector<uint32_t> positions(inReads.size(), 0);
	njh::iota<uint32_t>(positions, 0);
	njh::randomGenerator gen;
	if(std::numeric_limits<uint64_t>::max() != seed){
		gen.seedNum(seed);
	}
	auto selPositions = gen.unifRandSelectionVec(positions,sampleAmount * 2, false);
	std::vector<readObject> group1;
	std::vector<readObject> group2;
	for(auto pos : selPositions){
		if(group1.size() < sampleAmount){
			group1.emplace_back(inReads[pos]);
		}else{
			group2.emplace_back(inReads[pos]);
		}
	}
  uint64_t maxSize = 0;
  readVec::getMaxLength(group1,maxSize);
  readVec::getMaxLength(group2,maxSize);
  std::vector<std::unique_ptr<seqWithKmerInfo>> reads1 = createKmerReadVec(group1);
  std::vector<std::unique_ptr<seqWithKmerInfo>> reads2 = createKmerReadVec(group2);
  VecStr names;
  aligner alignerObj(maxSize,setUp.pars_.gapInfo_, setUp.pars_.scoring_, setUp.pars_.colOpts_.alignOpts_.countEndGaps_);
  alignerObj.processAlnInfoInput(setUp.pars_.alnInfoDirName_, setUp.pars_.verbose_);
  table distTab{VecStr{"read1", "read2", "kLen", "kmerDist", "id"}};
	std::unordered_map<std::string, std::unique_ptr<aligner>> aligners;
	std::mutex  alignerLock;
	std::function<
			comparison(const std::unique_ptr<seqWithKmerInfo> &,
					const std::unique_ptr<seqWithKmerInfo> &,
					std::unordered_map<std::string, std::unique_ptr<aligner>>&, aligner&)> getMismatchesFunc =
			[&alignerLock](const std::unique_ptr<seqWithKmerInfo> & read1, const std::unique_ptr<seqWithKmerInfo> & read2,
					std::unordered_map<std::string, std::unique_ptr<aligner>>& aligners,
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
				return aligners.at(threadId)->comp_;
			};
	std::vector<uint32_t> group1Positions(group1.size(), 0);
	njh::iota<uint32_t>(group1Positions, 0);
	std::unordered_map<uint32_t, std::vector<comparison>> allComps;
	std::mutex allCompsMut;
	table timingTab(VecStr{"seqLength", "n", "type", "time"});
	{
		njh::stopWatch watch;
		njh::concurrent::LockableQueue<uint32_t> group1Queue(group1Positions);
		auto alignFunc = [&alignerLock, &allCompsMut, &allComps](njh::concurrent::LockableQueue<uint32_t> & posQueue,
				const std::vector<std::unique_ptr<seqWithKmerInfo>> & r1,
				const std::vector<std::unique_ptr<seqWithKmerInfo>> & r2,
				std::unordered_map<std::string, std::unique_ptr<aligner>>& aligners,
				aligner &alignerObj){
			alignerLock.lock();
			auto threadId = estd::to_string(std::this_thread::get_id());
			//std::cout << threadId<< std::endl;
			if(aligners.find(threadId) == aligners.end()) {
				aligners.emplace(threadId, std::make_unique<aligner>(alignerObj));
			}
			alignerLock.unlock();
			uint32_t pos = 0;
			std::unordered_map<uint32_t, std::vector<comparison>> comps;
			while(posQueue.getVal(pos)){
				std::vector<comparison> r2Comps;
				for(const auto r2Pos : iter::range(r2.size())){
					aligners.at(threadId)->alignCache(getSeqBase(r2[r2Pos]),getSeqBase(r1[pos]), false);
					aligners.at(threadId)->profilePrimerAlignment(r2[r2Pos], getSeqBase(r1[pos]));
					r2Comps.emplace_back(aligners.at(threadId)->comp_);
				}
				comps[pos] = r2Comps;
			}
			{
				std::lock_guard<std::mutex> lock(allCompsMut);
				for(const auto & comp : comps){
					allComps[comp.first] = comp.second;
				}
			}
		};
		std::vector<std::thread> threads;
		for (uint32_t t = 0; t < numThreads; ++t) {
			threads.emplace_back(
					std::thread(alignFunc, std::ref(group1Queue), std::cref(reads1),
							std::cref(reads2), std::ref(aligners), std::ref(alignerObj)));
		}
		for(auto & t : threads){
			t.join();
		}

		timingTab.addRow(maxSize, sampleAmount, "alignment", watch.totalTime());
	}
	std::mutex distTabMut;
	{
		auto kDistFunc = [&distTabMut,&distTab, &allComps](njh::concurrent::LockableQueue<uint32_t> & posQueue,
				const std::vector<std::unique_ptr<seqWithKmerInfo>> & r1,
				const std::vector<std::unique_ptr<seqWithKmerInfo>> & r2,
				uint32_t k){
			uint32_t pos = 0;
			std::unordered_map<uint32_t, std::vector<double>> comps;
			while(posQueue.getVal(pos)){
				std::vector<double> r2Comps;
				for(const auto r2Pos : iter::range(r2.size())){
					auto dist = r1[pos]->compareKmers(getRef(r2[r2Pos]));
					r2Comps.emplace_back(dist.second);
				}
				comps[pos] = r2Comps;
			}
			{
				std::lock_guard<std::mutex> lock(distTabMut);
				for(const auto & comp : comps){
					for(const auto compPos : iter::range(comp.second.size())){
						distTab.addRow(r1[comp.first]->seqBase_.name_,
																			r2[compPos]->seqBase_.name_, k,
																			comp.second[compPos],
																			allComps[comp.first][compPos].distances_.eventBasedIdentity_);
					}
				}
			}
		};
		for (uint32_t k = kmerStart; k < kmerStop + 1; ++k) {
			allSetKmers(reads1, k, false);
			njh::stopWatch watch;
			allSetKmers(reads2, k, false);
			njh::concurrent::LockableQueue<uint32_t> group1Queue(group1Positions);
			std::vector<std::thread> threads;
			for (uint32_t t = 0; t < numThreads; ++t) {
				threads.emplace_back(
						std::thread(kDistFunc, std::ref(group1Queue), std::cref(reads1),
								std::cref(reads2), k));
			}
			for(auto & t : threads){
				t.join();
			}
			timingTab.addRow(maxSize, sampleAmount, "k-" + estd::to_string(k), watch.totalTime());
		}
	}
	auto timmingOpts = setUp.pars_.ioOptions_.out_;
	timmingOpts.outFilename_ = timmingOpts.outFilename_.string() + "_timming";
  distTab.outPutContents(TableIOOpts(setUp.pars_.ioOptions_.out_, "\t", distTab.hasHeader_));
  timingTab.outPutContents(TableIOOpts(timmingOpts, "\t", distTab.hasHeader_));
  for(const auto & align : aligners){
  	alignerObj.alnHolder_.mergeOtherHolder(align.second->alnHolder_);
  }

  alignerObj.processAlnInfoOutput(setUp.pars_.alnInfoDirName_, setUp.pars_.verbose_);
  return 0;
}

int kmerExpRunner::pidVsKmers(const njh::progutils::CmdArgs & inputCommands) {
  uint32_t kmerStart = 2;
  uint32_t kmerStop = 4;
  uint32_t numThreads = 2;//profileErrors
	seqSetUp setUp(inputCommands);
	setUp.pars_.ioOptions_.lowerCaseBases_ = "upper";
	setUp.pars_.ioOptions_.out_.outFilename_ = "pidVsKSim.tab.txt";
	setUp.processDefaultReader(true);
	setUp.pars_.ioOptions_.out_.outExtention_ = ".tab.txt";
	setUp.pars_.ioOptions_.out_.outFileFormat_ = "tab";
	setUp.processAlignerDefualts();
  setUp.setOption(numThreads, "--numThreads", "Number of Threads to Use");
  setUp.setOption(kmerStart, "--kmerLenStart", "Length for kmers to start at");
  setUp.setOption(kmerStop, "--kmerLenStop", "Length for kmers to stop at");
  setUp.finishSetUp(std::cout);
  setUp.startARunLog(setUp.pars_.directoryName_);
	if (kmerStop < kmerStart) {
		std::stringstream ss;
		ss << "Kmer stop must be greater than kmerStart" << std::endl;
		ss << "KmerStart: " << kmerStart << std::endl;
		ss << "KmerStop: " << kmerStop << std::endl;
		throw std::runtime_error { ss.str() };
	}
	if (setUp.pars_.ioOptions_.out_.outExists()
			&& !setUp.pars_.ioOptions_.out_.overWriteFile_) {
		std::stringstream ss;
		ss << "File: "
				<< setUp.pars_.ioOptions_.out_.outName()
				<< " already exists, use -overWrite to over write it" << std::endl;
		throw std::runtime_error { ss.str() };
	}
  SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<readObject>();
  uint64_t maxSize = 0;
  readVec::getMaxLength(inReads,maxSize);
  std::vector<std::unique_ptr<seqWithKmerInfo>> reads;
  VecStr names;
  names.reserve(inReads.size() + 1);
  for(const auto & read : inReads){
  	reads.emplace_back(std::make_unique<seqWithKmerInfo>(read.seqBase_));
  }
  aligner alignerObj(maxSize,setUp.pars_.gapInfo_, setUp.pars_.scoring_, setUp.pars_.colOpts_.alignOpts_.countEndGaps_);

  table distTab{VecStr{"read1", "read2", "kLen", "kmerDist", "id"}};
	std::unordered_map<std::string, std::unique_ptr<aligner>> aligners;
	std::mutex alignerLock;
	std::function<
			comparison(const std::unique_ptr<seqWithKmerInfo> &,
					const std::unique_ptr<seqWithKmerInfo> &,
					std::unordered_map<std::string, std::unique_ptr<aligner>>&, aligner&)> getMismatchesFunc =
			[&alignerLock](const std::unique_ptr<seqWithKmerInfo> & read1, const std::unique_ptr<seqWithKmerInfo> & read2,
					std::unordered_map<std::string, std::unique_ptr<aligner>>& aligners,
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
				return aligners.at(threadId)->comp_;
			};

	auto distances = getDistanceNonConst(reads, numThreads, getMismatchesFunc,
			aligners, alignerObj);
	std::function<
			double(const std::unique_ptr<seqWithKmerInfo> &,
					const std::unique_ptr<seqWithKmerInfo> &)> disFun =
			[](const std::unique_ptr<seqWithKmerInfo> & read1,
					const std::unique_ptr<seqWithKmerInfo> & read2) {
				auto dist = read1->compareKmers(*read2);
				return dist.second;
			};
	for (uint32_t k = kmerStart; k < kmerStop + 1; ++k) {
		allSetKmers(reads, k, false);
		auto kDists = getDistance(reads, numThreads, disFun);
		for (const auto rowPos : iter::range(reads.size())) {
			for (const auto colPos : iter::range(rowPos)) {
				distTab.addRow(reads[rowPos]->seqBase_.name_,
						reads[colPos]->seqBase_.name_, k, kDists[rowPos][colPos],
						distances[rowPos][colPos].distances_.eventBasedIdentity_);
			}
		}
	}

  distTab.outPutContents(TableIOOpts(setUp.pars_.ioOptions_.out_, "\t", distTab.hasHeader_));
  return 0;
}

int kmerExpRunner::scoveViaKmers(const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	setUp.pars_.ioOptions_.out_.outFilename_ = "kDist.tab.txt";
	double cutOff = 0.75;
	uint32_t sizeCutOff = 1;
	setUp.setOption(sizeCutOff, "--sizeCutOff", "Cluster Size Cut Off");
	setUp.setOption(cutOff, "--cutOff", "cutOff");
	setUp.processDefaultReader(true);
	setUp.processDirectoryOutputName(true);
	setUp.processAlignerDefualts();
	setUp.processRefFilename(true);
  setUp.finishSetUp(std::cout);
  setUp.startARunLog(setUp.pars_.directoryName_);
  SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<readObject>();
  uint64_t maxSize = 0;
  readVec::getMaxLength(inReads,maxSize);
  auto refs = SeqInput::getReferenceSeq(setUp.pars_.refIoOptions_, maxSize);
  std::vector<std::unique_ptr<seqWithKmerInfo>> refReads;
  for(const auto & read : refs){
  	refReads.emplace_back(std::make_unique<seqWithKmerInfo>(read.seqBase_));
  }
  std::vector<std::unique_ptr<seqWithKmerInfo>> reads;
  VecStr names;
  names.reserve(inReads.size() + 1);
  for(const auto & read : inReads){
  	reads.emplace_back(std::make_unique<seqWithKmerInfo>(read.seqBase_));
  }
  allSetKmers(reads, setUp.pars_.colOpts_.kmerOpts_.kLength_,false);
  allSetKmers(refReads, setUp.pars_.colOpts_.kmerOpts_.kLength_, false);
  aligner alignerObj(maxSize,setUp.pars_.gapInfo_, setUp.pars_.scoring_, setUp.pars_.colOpts_.alignOpts_.countEndGaps_);

  table distTab{VecStr{"kmerDist", "alnScore", "color"}};
	std::vector<double> kDists;
	std::vector<int32_t> scores;
	kDists.reserve(refReads.size());
	scores.reserve(refReads.size());
  for(const auto readPos : iter::range(reads.size())){
  	if(readPos%10 == 0){
  		std::cout << "on " << readPos << " of " << reads.size() << "\r";
  		std::cout.flush();
  	}
  	auto & read = reads[readPos];
  	kDists.clear();
  	scores.clear();
  	for(const auto & ref : refReads){
  		auto dist = read->compareKmers(*ref);
  		alignerObj.alignCache(ref->seqBase_, read->seqBase_, setUp.pars_.local_);
  		kDists.emplace_back(dist.second);
  		scores.emplace_back(alignerObj.parts_.score_);
  	}
  	auto maxKDist = vectorMaximum(kDists);
  	auto maxScore = vectorMaximum(scores);
  	std::string outColor = "";
  	for(const auto pos : iter::range(refReads.size())){
  		if(kDists[pos] == maxKDist && scores[pos] == maxScore){
  			outColor = "#14B814";
  		}else if (kDists[pos] != maxKDist && scores[pos] == maxScore){
  			outColor = "#115BEE";
  		}else if (kDists[pos] == maxKDist && scores[pos] != maxScore){
  			outColor = "#F83AB9";
  		}else{
  			outColor = "#000000";
  		}
  		distTab.content_.emplace_back(toVecStr(kDists[pos], scores[pos], outColor));
  	}
  }
  std::cout << std::endl;
  distTab.outPutContents(TableIOOpts(OutOptions(setUp.pars_.ioOptions_.out_.outFilename_, ".tab.txt"), "\t", distTab.hasHeader_));
  return 0;
}//scoveViaKmers




int kmerExpRunner::getKmerSetDistBetween(const njh::progutils::CmdArgs & inputCommands){
	bfs::path file1 = "";
	bfs::path file2 = "";
	std::string name1 = "";
	std::string name2 = "";
	uint32_t kLength = 5;
	bool getRevComp = false;
	OutOptions outOpts(bfs::path(""), ".tab.txt");
	seqSetUp setUp(inputCommands);
	setUp.setOption(file1, "--file1", "file1", true);
	setUp.setOption(file2, "--file2", "file2", true);
	setUp.setOption(name1, "--name1", "name1");
	setUp.setOption(name2, "--name2", "name2");
	setUp.processWritingOptions(outOpts);
	setUp.setOption(kLength, "--kLength", "Kmer Length");
	setUp.setOption(getRevComp, "--getRevComp", "getRevComp");

  setUp.finishSetUp(std::cout);

  SeqIOOptions file1Opts(file1, SeqIOOptions::getInFormatFromFnp(file1));
  SeqIOOptions file2Opts(file2, SeqIOOptions::getInFormatFromFnp(file2));
  if("" == name1){
  	name1 = bfs::basename(file1);
  }
  if("" == name2){
  	name2 = bfs::basename(file2);
  }
  OutputStream out(outOpts);
  seqInfo seq;

  kmerInfo file1Info;
  file1Info.kLen_ = kLength;

  SeqInput file1Reader(file1Opts);
  file1Reader.openIn();
  while(file1Reader.readNextRead(seq)){
  	file1Info.updateKmers(seq.seq_, false);
  }

  kmerInfo file2Info;
  file2Info.kLen_ = kLength;

  SeqInput file2Reader(file2Opts);
  file2Reader.openIn();
  while(file2Reader.readNextRead(seq)){
  	file2Info.updateKmers(seq.seq_, getRevComp);
  }

  std::pair<uint32_t, double> forDist = file1Info.compareKmers(file2Info);
  std::pair<uint32_t, double> revDist;
  if(getRevComp){
  	revDist = file1Info.compareKmersRevComp(file2Info);
  }

  out << "set1\tset2\tkLength\tshared\tindex";
  if(getRevComp){
  	out << "\tsharedRev\tindexRev";
  }
  out << std::endl;
  out << name1
  		<< "\t" << name2
			<< "\t" << kLength
			<< "\t" << forDist.first
			<< "\t" << forDist.second;
  if(getRevComp){
    out << "\t" << revDist.first
  			<< "\t" << revDist.second;
  }
  out << std::endl;

  return 0;
}

int kmerExpRunner::getKmerDistTwoSeqs(const njh::progutils::CmdArgs & inputCommands){
	seqSetUp setUp(inputCommands);
	std::string seq1 = "";
	std::string seq2 = "";
	uint32_t kLength = 5;
	bool getReverseDistance = false;
	setUp.processSeq(seq1, "--seq1", "Sequence 1", true);
	setUp.processSeq(seq2, "--seq2", "Sequence 2", true);
	setUp.setOption(kLength, "--kLength", "Kmer Length");
	setUp.setOption(getReverseDistance, "--getReverseDistance", "Get Reverse Complement Distance");
  setUp.finishSetUp(std::cout);

  seqWithKmerInfo seq1Obj(seqInfo("seq1", seq1), kLength, getReverseDistance);
  seqWithKmerInfo seq2Obj(seqInfo("seq2", seq2), kLength, getReverseDistance);
  std::pair<uint32_t,double> distance;
  if(getReverseDistance){
  	distance = seq1Obj.compareKmers(seq2Obj);
  }else{
  	distance = seq1Obj.compareKmers(seq2Obj);
  }

  std::cout << "KmersShared:" << distance.first << "(" << std::min(seq1.size(), seq2.size()) - kLength + 1 << ")" << "\n";
  std::cout << "KmersDistance:" << distance.second << "\n";


  return 0;
}


class testRead : public seqWithKmerInfo {
public:
	testRead(const seqInfo & seqBase): seqWithKmerInfo(seqBase){

	}

	charCounter counter_;

	void setCount(const std::vector<char> & alph){
		counter_ = charCounter(alph);
		counter_.increaseCountByString(seqBase_.seq_, seqBase_.cnt_);
	}
};
/*
std::function<double(const std::unique_ptr<seqWithKmerInfo> & read1,
		const std::unique_ptr<seqWithKmerInfo> & read2,
		uint32_t windowSize, uint32_t windowStepSize)> disFun =
		[](const std::unique_ptr<seqWithKmerInfo> & read1,
				const std::unique_ptr<seqWithKmerInfo> & read2,
				uint32_t windowSize, uint32_t windowStepSize){
	auto scanDist = read1->slideCompareKmers(*read2, windowSize, windowStepSize);
	double averageDist = 0;
  for(const auto & distPos : iter::range(scanDist.size())){
  	averageDist += scanDist[distPos].second;
	}
	return averageDist/scanDist.size(); scoveViaKmers
};*/
int kmerExpRunner::scaningKmerDist(const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	uint32_t windowSize = 100;
	uint32_t windowStepSize = 10;
	uint32_t numThreads = 1;
	setUp.processDefaultReader(true);
	setUp.setOption(numThreads, "-numThreads", "numThreads");
	setUp.setOption(setUp.pars_.colOpts_.kmerOpts_.kLength_, "-k,-kLength", "kLength");
	setUp.setOption(windowSize, "-windowSize", "windowSize");
	setUp.setOption(windowStepSize, "-windowStepSize", "windowStepSize");
	//setUp.processDirectoryOutputName(true);
	setUp.processDebug();
	setUp.processSeq(true);
	setUp.finishSetUp(std::cout);
	//setUp.startARunLog(setUp.pars_.directoryName_);

	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<readObject>();
	auto reads = createKmerReadVec(inReads, setUp.pars_.colOpts_.kmerOpts_.kLength_, false);
	std::unique_ptr<seqWithKmerInfo> seq = std::make_unique<seqWithKmerInfo>(
			setUp.pars_.seqObj_.seqBase_, setUp.pars_.colOpts_.kmerOpts_.kLength_, false);
	seq->setKmers(setUp.pars_.colOpts_.kmerOpts_.kLength_, false);
	std::cout << "readName\taverageDist\tfullDist\n";
	for (const auto & read : reads) {
		auto scanDist = read->slideCompareKmers(*seq, windowSize, windowStepSize);
		double averageDist = 0;
		for (const auto distPos : iter::range(scanDist.size())) {
			averageDist += scanDist[distPos].second;
		}
		averageDist /= scanDist.size();
		auto fullDist = seq->compareKmers(*read);
		std::cout << read->seqBase_.name_ << "\t" << averageDist << "\t"
				<< fullDist.second << "\n";
	}
	return 0;
}


int kmerExpRunner::kDistVsNucDist(const njh::progutils::CmdArgs & inputCommands){
	seqSetUp setUp(inputCommands);
  uint32_t kLength = 10;
  uint32_t numThreads = 1;
	setUp.processDefaultReader(true);
	setUp.setOption(kLength, "-k,--kLength", "Kmer Length");
	setUp.setOption(numThreads, "-t,--numThreads", "Number of Threads to Use");
  setUp.finishSetUp(std::cout);

  SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<readObject>();
  std::vector<testRead> reads;
  for(const auto & read : inReads){
  	reads.emplace_back(read.seqBase_);
  }
  njh::for_each(reads, [&](testRead & read){
  	read.setCount({'A','C','G','T'});
  	read.counter_.setFractions();
  	read.setKmers(kLength, false);});

  std::function<std::pair<double, double>(const testRead & read1,
  		const testRead & read2)> disFun =
  		[](const testRead & read1,
  				const testRead & read2){
  	auto dist = read1.compareKmers(read2);
  	auto nucDist = read1.counter_.getFracDifference(read2.counter_, read1.counter_.alphabet_);
  	return std::pair<double,double>{nucDist, dist.second};
  	//return dist.second;
  };
  auto distances = getDistance(reads, numThreads, disFun);
  std::ofstream outInfoFile("nucVsNucDist_" + bfs::basename(setUp.pars_.ioOptions_.firstName_));
  outInfoFile << "read1\tread2\tnucDist\tDist\n";
  for(const auto pos : iter::range(distances.size())){
  	for(const auto secondPos : iter::range(pos)){
  		outInfoFile<< reads[pos].seqBase_.name_
  				<< "\t" << reads[secondPos].seqBase_.name_
  				<< "\t" << distances[pos][secondPos].first
  				<< "\t" << distances[pos][secondPos].second
  				<< "\n";
  	}
  }
  return 0;
}

struct kRefCluster{

	kRefCluster(const seqInfo & seqBase): refSeq_(std::make_unique<seqWithKmerInfo>(seqBase)){}
	std::unique_ptr<seqWithKmerInfo> refSeq_;
	std::vector<std::unique_ptr<seqWithKmerInfo>> reads_;

};

int kmerExpRunner::clostestKmerDist(const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	setUp.processDefaultReader(true);
	setUp.processRefFilename(true);
	setUp.processDirectoryOutputName(true);
	setUp.processVerbose();
	setUp.setOption(setUp.pars_.colOpts_.kmerOpts_.kLength_, "-k,-kLength", "kLength");
  setUp.finishSetUp(std::cout);
  uint64_t maxLength = 0;
  SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<readObject>();
  auto inRefSeqs = SeqInput::getReferenceSeq(setUp.pars_.refIoOptions_, maxLength);
  readVec::getMaxLength(inReads, maxLength);
  std::vector<kRefCluster> refSeqs;
  for(const auto & ref : inRefSeqs){
  	refSeqs.emplace_back(kRefCluster(ref.seqBase_));
  	refSeqs.back().refSeq_->setKmers(setUp.pars_.colOpts_.kmerOpts_.kLength_, false);
  }
  std::vector<std::unique_ptr<seqWithKmerInfo>> reads;
  for(const auto & read : inReads){
  	reads.emplace_back(std::make_unique<seqWithKmerInfo>(read.seqBase_));
  }
  allSetKmers(reads, setUp.pars_.colOpts_.kmerOpts_.kLength_, false);
  std::vector<std::unique_ptr<seqWithKmerInfo>> rejectReads;
  for(const auto readPos : iter::range(reads.size())){
  	if(readPos % 50 == 0){
  		std::cout << "Currently on " << readPos << " of " << reads.size() << "\n";
  	}
  	std::pair<uint32_t, double> best = {0,0};
  	uint32_t bestIndex = std::numeric_limits<uint32_t>::max();
  	for(const auto refPos : iter::range(refSeqs.size())){
  		auto current = refSeqs[refPos].refSeq_->compareKmers(*reads[readPos]);
  		if(current.second > best.second){
  			best = current;
  			bestIndex = refPos;
  		}
  	}
  	if(bestIndex != std::numeric_limits<uint32_t>::max()){
  		refSeqs[bestIndex].reads_.emplace_back(std::forward<std::unique_ptr<seqWithKmerInfo>>(reads[readPos]));
  	}else{
  		rejectReads.emplace_back(std::forward<std::unique_ptr<seqWithKmerInfo>>(reads[readPos]));
  	}
  }

  for(const auto & ref : refSeqs){
  	std::ofstream currentReadFile;
  	openTextFile(currentReadFile, setUp.pars_.directoryName_ + ref.refSeq_->seqBase_.name_,
  			"fastq", setUp.pars_.ioOptions_.out_);
  	for(const auto & read : ref.reads_){
  		read->seqBase_.outPutFastq(currentReadFile);
  	}
  }
  if(!rejectReads.empty()){
  	std::ofstream currentReadFile;
  	openTextFile(currentReadFile, setUp.pars_.directoryName_ + "rejected",
  			"fastq", setUp.pars_.ioOptions_.out_);
  	for(const auto & read : rejectReads){
  		read->seqBase_.outPutFastq(currentReadFile);
  	}
  }
  setUp.logRunTime(std::cout);
  return 0;
}



int kmerExpRunner::profileScanningKmerDist(const njh::progutils::CmdArgs & inputCommands){
  seqSetUp setUp(inputCommands);
  uint32_t windowSize = 50;
  uint32_t windowStep = 10;
  std::string outFilename = "";
  setUp.setOption(outFilename, "-o,--outFile","name of the output file", true );
  setUp.setOption(windowSize, "--windowSize", "Size of scanning kmer window");
  setUp.setOption(windowStep, "--windowStep", "size of the scanning kmer dist step");
  std::string seq2 = "";
  setUp.processSeq(seq2, "-seq2", "Seq2", true);
  setUp.processSeq(true);
  setUp.setOption(setUp.pars_.colOpts_.kmerOpts_.kLength_, "-k,--kmerLength", "Length for kmers");
  setUp.processWritingOptions();
  bool append = false;
  setUp.setOption(append, "--append", "Append File");
  setUp.processVerbose();
  setUp.finishSetUp(std::cout);

  seqWithKmerInfo seq1Obj(setUp.pars_.seqObj_.seqBase_);
  seqWithKmerInfo seq2Obj(seqInfo("Seq2", seq2));
  seq1Obj.setKmers(setUp.pars_.colOpts_.kmerOpts_.kLength_, false);
  seq2Obj.setKmers(setUp.pars_.colOpts_.kmerOpts_.kLength_, false);

  std::ofstream outFile;
  njh::appendAsNeeded(outFilename, ".tab.txt");
  bool fileExists = bfs::exists(outFilename);
	njh::files::openTextFile(outFile, outFilename,
			setUp.pars_.ioOptions_.out_.overWriteFile_, append,
			setUp.pars_.ioOptions_.out_.exitOnFailureToWrite_);
  if(!fileExists || setUp.pars_.ioOptions_.out_.overWriteFile_){
  	outFile << "fullSeq\tseqParts\tcase\tkLen\twindowNumber\tkShared\tkDist" << std::endl;
  }
  auto seq1Dist = seq1Obj.slideCompareSubKmersToFull(seq2Obj, windowSize, windowStep);
  auto seq2Dist = seq2Obj.slideCompareSubKmersToFull(seq1Obj, windowSize, windowStep);
  for(const auto & dist : seq2Dist){
  	outFile << "seq1\tseq2\tseq1vsseq2\t" << setUp.pars_.colOpts_.kmerOpts_.kLength_
  			<< "\t" << dist.first
  			<< "\t" << dist.second.first
				<< "\t" << dist.second.second
				<< "\n";
  }
  for(const auto & dist : seq2Dist){
  	outFile << "seq2\tseq1\tseq2vsseq1\t" << setUp.pars_.colOpts_.kmerOpts_.kLength_
  			<< "\t" << dist.first
  			<< "\t" << dist.second.first
				<< "\t" << dist.second.second
				<< "\n";
  }
	return 0;
}

int kmerExpRunner::profileKmerAccerlation(const njh::progutils::CmdArgs & inputCommands){
  seqSetUp setUp(inputCommands);
  uint32_t numThreads = 1;
  setUp.processDefaultReader(true);
  //setUp.processRefFilename(true);
  setUp.setOption(numThreads, "-t,--numThreads", "Number of Threads to Use");
  setUp.setOption(setUp.pars_.colOpts_.kmerOpts_.kLength_, "-k,--kmerLength", "Length for kmers");
  setUp.pars_.ioOptions_.out_.outExtention_ = ".tab.txt";

  setUp.finishSetUp(std::cout);
  SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<readObject>();

  njh::stopWatch watch;
  watch.setLapName("index kmers");
  std::vector<std::unique_ptr<seqWithKmerInfo>> reads;
  for(const auto & read : inReads){
  	reads.emplace_back(std::make_unique<seqWithKmerInfo>(read.seqBase_));
  }

  uint64_t maxLength = 0;
  readVec::getMaxLength(inReads, maxLength);


  std::ofstream outFile;
  openTextFile(outFile, setUp.pars_.ioOptions_.out_);
  //outFile << "seq1\tseq2\tkLen\tkDist\tkShared\n";
  outFile << njh::conToStr(VecStr{"seq1", "seq2","group", "kLen", "seq1Group", "seq2Group", "case", "kDist", "kShared"}, "\t") << std::endl;
  /*std::function<std::pair<uint32_t, double>(const std::unique_ptr<seqWithKmerInfo> & read1,
  		const std::unique_ptr<seqWithKmerInfo> & read2)> disFun =
  		[](const std::unique_ptr<seqWithKmerInfo> & read1,
  				const std::unique_ptr<seqWithKmerInfo> & read2){
  	auto dist = read1->compareKmers(*read2); profilePacbioReads
  	return dist;
  };
  auto distances = getDistance(reads, numThreads, disFun);*/
  std::function<std::pair<uint32_t,double>(const std::unique_ptr<seqWithKmerInfo> & read1,
			const std::unique_ptr<seqWithKmerInfo> & read2)> disFun =
			[](const std::unique_ptr<seqWithKmerInfo> & read1,
					const std::unique_ptr<seqWithKmerInfo> & read2){
		auto dist = read1->compareKmers(*read2);
		return dist;
	};
  auto processGroup = [](const std::string & name){
  	auto underPos = name.find("_");
  	if(underPos != std::string::npos){
  		return name.substr(0, underPos);
  	}else{
  		return name.substr(0, name.find("."));
  	}
  };
  for(const auto k : iter::range<uint32_t>(2, setUp.pars_.colOpts_.kmerOpts_.kLength_ + 1)){
  	std::cout << "Currently on Kmer Length " << k << '\r';
  	std::cout.flush();
  	allSetKmers(reads, k,false);
  	auto distances = getDistance(reads, numThreads, disFun);
  	for(const auto pos : iter::range(distances.size())){
  		for(const auto subPos : iter::range(distances[pos].size())){
  			auto seq1Name = reads[pos]->seqBase_.name_;
  			auto seq2Name = reads[subPos]->seqBase_.name_;
  			auto seq1Group = processGroup(seq1Name);
  			auto seq2Group = processGroup(seq2Name);
  			auto dist = distances[pos][subPos];
  			outFile << seq1Name
  					<< "\t" << seq2Name
						<< "\t" << seq1Name << "_v_" << seq2Name
						<< "\t" << k
						<< "\t" << seq1Group
						<< "\t" << seq2Group
						<< "\t" << (seq1Group == seq2Group ? "same" : "different")
						<< "\t" << dist.second
						<< "\t" << dist.first
						<< "\n";
  		}
  	}
  }
  std::cout << std::endl;
	return 0;
}



int kmerExpRunner::getNewScanKmerDist(const njh::progutils::CmdArgs & inputCommands){
  seqSetUp setUp(inputCommands);
  //uint32_t numThreads = 1; graphTest
  std::string filename = "";
  uint32_t windowSize = 50;
  uint32_t windowStep = 10;
  setUp.setOption(windowSize, "--windowSize", "Size of scanning kmer window");
  setUp.setOption(windowStep, "--windowStep", "size of the scanning kmer dist step");
  setUp.setOption(filename, "-o,--outFile", "Out Filename", true);
  setUp.processDefaultReader(true);
  setUp.processRefFilename(true);
  setUp.setOption(setUp.pars_.colOpts_.kmerOpts_.kLength_, "-k,--kmerLength", "Length for kmers");
  setUp.processVerbose();
  setUp.finishSetUp(std::cout);
  SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<readObject>();

  njh::stopWatch watch;
  watch.setLapName("index kmers");
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

  std::ofstream outFile;
  openTextFile(outFile, filename, ".tab.txt", setUp.pars_.ioOptions_.out_);
  outFile << "seq\tref\twindowNum\tkDist\tkShared\n";
	allSetKmers(refReads, setUp.pars_.colOpts_.kmerOpts_.kLength_,false);
	allSetKmers(reads, setUp.pars_.colOpts_.kmerOpts_.kLength_,false);
	for(const auto & read : reads){
		for(const auto & ref : refReads){
			auto dists = ref->slideCompareSubKmersToFull(*read, windowSize, windowStep);
			for(const auto & dist : dists){
				outFile << read->seqBase_.name_
						<< "\t" << ref->seqBase_.name_
						<< "\t" << dist.first
						<< "\t" << dist.second.first
						<< "\t" << dist.second.second
						<< "\n";
			}
		}
	}
	if(setUp.pars_.verbose_){
		watch.logLapTimes(std::cout, true, 6, true);
	}

	return 0;
}





int kmerExpRunner::getKmerDistAgainstRef(const njh::progutils::CmdArgs & inputCommands){
  bool dontSkipSameName = false;
  bool getRevComp = false;
  uint32_t numThreads = 1;
  uint32_t kLenStart = 2;
  uint32_t kLenStop = std::numeric_limits<uint32_t>::max();
  uint32_t kLenStep = 1;
	OutOptions outOpts(bfs::path(""), ".tab.txt");
  seqSetUp setUp(inputCommands);
	setUp.setOption(kLenStop, "--kLenStop", "kmer Length Stop", true);
	setUp.setOption(kLenStart, "--kLenStart", "kmer Length Start");
	setUp.setOption(kLenStep, "--kLenStep", "kmer Length Step");
	setUp.setOption(numThreads, "--numThreads", "number of threads to use");
	setUp.setOption(getRevComp, "--getRevComp", "Get Rev Comp");
	setUp.setOption(dontSkipSameName, "--dontSkipSameName", "Don't skip comparison if they have the same name");
	setUp.processVerbose();
  setUp.processReadInNames(true);
  setUp.processRefFilename(true);
  setUp.processKmerLenOptions();
  setUp.processWritingOptions(outOpts);
  setUp.finishSetUp(std::cout);
  SeqInput reader(setUp.pars_.ioOptions_);
	auto inReads = reader.readAllReads<readObject>();

  njh::stopWatch watch;
  watch.setLapName("index kmers");
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
  OutputStream outFile(outOpts);
  outFile << "seq\tref\tkLen\tkDist\tkShared\n";
  for(const auto k : iter::range<uint32_t>(kLenStart, kLenStop + 1, kLenStep)){
  	if(setUp.pars_.verbose_){
  		std::cerr << "Currently on Kmer Length " << k << '\r';
  		std::cerr.flush();
  	}
		allSetKmers(refReads, k, false);
		allSetKmers(reads, k, getRevComp);

		std::mutex outFileMut;

		AllByAllPairFactory allFactory(reads.size(), refReads.size());

		std::function<void()> getDist = [&allFactory,&getRevComp,&refReads,&reads,
										&outFileMut,&outFile,&k,&dontSkipSameName](){
			AllByAllPairFactory::AllByAllPair pair;
			while(allFactory.setNextPair(pair)){
				const auto & seq = reads[pair.row_];
				const auto & ref = refReads[pair.col_];
				if(seq->seqBase_.name_ == ref->seqBase_.name_ && !dontSkipSameName){
					continue;
				}
				auto dist = ref->compareKmers(*seq);
				std::pair<uint32_t,double> revDist;
				if(getRevComp){
					revDist = ref->compareKmersRevComp(*seq);
				}
				{
					std::lock_guard<std::mutex> lock(outFileMut);
					outFile << seq->seqBase_.name_
							<< "\t" << ref->seqBase_.name_
							<< "\t" << k
							<< "\t" << dist.second
							<< "\t" << dist.first << "\n";
					if(getRevComp){
		  			outFile << seq->seqBase_.name_ << "_revComp"
		  					<< "\t" << ref->seqBase_.name_
								<< "\t" << k
								<< "\t" << revDist.second
		  					<< "\t" << revDist.first << "\n";
					}
				}
			}
		};
		njh::concurrent::runVoidFunctionThreaded(getDist, numThreads);


  }
  if(setUp.pars_.verbose_){
  	std::cerr << std::endl;
  }
	return 0;
}


int kmerExpRunner::readingDistanceCheck(const njh::progutils::CmdArgs & inputCommands){
  seqSetUp setUp(inputCommands);
  uint32_t numThreads = 1;
  std::string filename = "";
  setUp.processDefaultReader(true);
  setUp.setOption(numThreads, "-t,--numThreads", "Number of Threads to Use");
  setUp.setOption(setUp.pars_.colOpts_.kmerOpts_.kLength_, "-k,--kmerLength", "Length for kmers");
  setUp.setOption(filename, "-i,--inFile", "In Filename", true);
  setUp.finishSetUp(std::cout);
  njh::stopWatch watch;
  watch.setLapName("index kmers");
  auto reads = createKmerReadVec(setUp.pars_.ioOptions_,setUp.pars_.colOpts_.kmerOpts_.kLength_,false );
  std::function<double(const std::unique_ptr<seqWithKmerInfo> & read1,
  		const std::unique_ptr<seqWithKmerInfo> & read2)> disFun =
  		[](const std::unique_ptr<seqWithKmerInfo> & read1,
  				const std::unique_ptr<seqWithKmerInfo> & read2){
  	auto dist = read1->compareKmers(*read2);
  	return dist.second;
  };
  watch.startNewLap("regular distance");
  auto distances = getDistance(reads, numThreads, disFun);
  std::ifstream inFile(filename);
  if(!inFile){
  	std::cerr << njh::bashCT::red <<
  			njh::bashCT::bold << "Error in opening " << filename
				<< njh::bashCT::reset << std::endl;
  	exit(1);
  }
  std::vector<std::vector<double>> inDist = readDistanceMatrix(inFile);
  bool match = true;
  for(const auto i : iter::range(inDist.size())){
  	for(const auto j : iter::range(inDist[i].size())){
  		if(roundDecPlaces(inDist[i][j], 4) != roundDecPlaces(distances[i][j], 4)){
  			std::cout << "inDist: " << inDist[i][j] << ":" << roundDecPlaces(inDist[i][j], 4)  << std::endl;
  			std::cout << "outDist: " << distances[i][j] << ":" << roundDecPlaces(distances[i][j], 4)  << std::endl;
  			match =false;
  			break;
  		}
  	}
  }
  if(match){
  	std::cout << njh::bashCT::green <<
  			njh::bashCT::bold << "match"<< filename
				<< njh::bashCT::reset << std::endl;
  }else{
  	std::cout << njh::bashCT::red <<
  			njh::bashCT::bold << "Don't match" << filename
				<< njh::bashCT::reset << std::endl;
  }

	return 0;
}


int kmerExpRunner::writingDistanceCheck(const njh::progutils::CmdArgs & inputCommands){
  seqSetUp setUp(inputCommands);
  uint32_t numThreads = 1;
  setUp.processDefaultReader(false);
  setUp.setOption(numThreads, "-t,--numThreads", "Number of Threads to Use");
  setUp.setOption(setUp.pars_.colOpts_.kmerOpts_.kLength_, "-k,--kmerLength", "Length for kmers");
  setUp.finishSetUp(std::cout);

	std::vector<readObject> inReads;
  if(setUp.pars_.ioOptions_.firstName_ == ""){
    njh::randomGenerator gen;
    VecStr rSeqs = simulation::evenRandStrs(40, std::vector<char>{'A', 'C', 'G', 'T'}, gen, 20);
    inReads = vecStrToReadObjs<readObject>(rSeqs, "Seq");
  }else{
  	SeqInput reader(setUp.pars_.ioOptions_);
  	reader.openIn();
  	inReads = reader.readAllReads<readObject>();
  }

  njh::stopWatch watch;
  watch.setLapName("index kmers");
  std::vector<std::unique_ptr<seqWithKmerInfo>> reads;
  for(const auto & read : inReads){
  	reads.emplace_back(std::make_unique<seqWithKmerInfo>(read.seqBase_));
  }
  allSetKmers(reads, setUp.pars_.colOpts_.kmerOpts_.kLength_,false);
  std::function<double(const std::unique_ptr<seqWithKmerInfo> & read1,
  		const std::unique_ptr<seqWithKmerInfo> & read2)> disFun =
  		[](const std::unique_ptr<seqWithKmerInfo> & read1,
  				const std::unique_ptr<seqWithKmerInfo> & read2){
  	auto dist = read1->compareKmers(*read2);
  	return dist.second;
  };
  watch.startNewLap("regular distance");
  auto distances = getDistance(reads, numThreads, disFun);
  std::ofstream outFile;
  openTextFile(outFile, setUp.pars_.ioOptions_.out_.outFilename_.string(), ".tab.txt", setUp.pars_.ioOptions_.out_);
  writeDistanceMatrix(outFile, distances);
  if(setUp.pars_.ioOptions_.firstName_ == ""){
    std::ofstream outSeqFile;
    openTextFile(outSeqFile, setUp.pars_.ioOptions_.out_.outFilename_.string() + ".fasta", ".fasta", setUp.pars_.ioOptions_.out_);
    for(const auto & read : reads){
    	read->seqBase_.outPutSeq(outSeqFile);
    }
  }
	return 0;
}
//scoveViaKmers

int kmerExpRunner::getKmerDistStats(
		const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	TableIOOpts tabOpts = TableIOOpts::genTabFileOut("", true);
	setUp.processWritingOptions(tabOpts.out_);
	setUp.processReadInNames(setUp.readInFormatsAvailable_);
	setUp.processSeq(true);
	setUp.processKmerLenOptions();
	setUp.finishSetUp(std::cout);
	seqWithKmerInfo compare(setUp.pars_.seqObj_.seqBase_, setUp.pars_.colOpts_.kmerOpts_.kLength_,
			false);
	auto outStatsTab = getKmerStatsOnFile(setUp.pars_.ioOptions_, compare);
	outStatsTab.outPutContents(tabOpts);
	return 0;
}

int kmerExpRunner::getKmerDistStatsMultiple(
		const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	std::string multipleFiles = "";
	TableIOOpts tabOpts = TableIOOpts::genTabFileOut("", true);
	setUp.processWritingOptions(tabOpts.out_);
	setUp.setOption(multipleFiles, "-files", "multipleFiles");
	setUp.processVerbose();
	setUp.processSeq(true);
	setUp.processKmerLenOptions();
	setUp.finishSetUp(std::cout);

	seqWithKmerInfo compare(setUp.pars_.seqObj_.seqBase_, setUp.pars_.colOpts_.kmerOpts_.kLength_,
			false);
	auto toks = tokenizeString(multipleFiles, ",");
	table outStats;
	for (const auto & file : toks) {
		SeqIOOptions seqOpts;
		seqOpts.firstName_ = file;
		seqOpts.inFormat_ = SeqIOOptions::getInFormat(njh::files::getExtension(file));
		if (outStats.content_.empty()) {
			outStats = getKmerStatsOnFile(seqOpts, compare);
		} else {
			outStats.rbind(getKmerStatsOnFile(seqOpts, compare), false);
		}
	}
	outStats.outPutContents(tabOpts);
	return 0;
}

} /* namespace njhseq */
