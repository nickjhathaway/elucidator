/*
 * pcrSimulationUtils.cpp
 *
 *  Created on: Jan 8, 2016
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
#include "pcrSimulationUtils.hpp"



namespace njhseq {
namespace sim {

ReadLenNormalDistribution::ReadLenNormalDistribution(double mean, double std,
		int32_t minLen) :
		minAllowed_(minLen), rgen_(std::random_device()()), ndis_(
				std::normal_distribution<>(mean, std)) {

}

ReadLenNormalDistribution::ReadLenNormalDistribution(double mean, double std,
		int32_t minLen, uint64_t seed) :
		minAllowed_(minLen), rgen_(seed), ndis_(
				std::normal_distribution<>(mean, std)) {

}

uint32_t ReadLenNormalDistribution::operator ()() {
	uint32_t res = std::max<double>(1, std::round(ndis_(rgen_)));
	while (res < minAllowed_) {
		res = std::round(ndis_(rgen_));
	}
	return res;
}

std::vector<std::pair<size_t, uint32_t>> genFragPosSizes(size_t seqLen,
		ReadLenNormalDistribution & uiNormDist, njh::randomGenerator & gen) {
	//std::cout << "__LINE__: " << __LINE__ << std::endl;
	std::vector<std::pair<size_t, uint32_t>> readPositionsLens;
	auto frontSeqSize = uiNormDist();
	auto backSeqSize = uiNormDist();
	//std::cout << "__LINE__: " << __LINE__ << std::endl;
	if(frontSeqSize <= seqLen){
		readPositionsLens.emplace_back(0, frontSeqSize);
	}
	//std::cout << "__LINE__: " << __LINE__ << std::endl;
	if(backSeqSize <= seqLen){
		readPositionsLens.emplace_back(seqLen - backSeqSize, backSeqSize);
	}
	//std::cout << "__LINE__: " << __LINE__ << std::endl;
	if(frontSeqSize + backSeqSize >= seqLen){
		return readPositionsLens;
	}
	//std::cout << "__LINE__: " << __LINE__ << std::endl;
	std::vector<size_t> positions(seqLen - frontSeqSize - backSeqSize);
	njh::iota<size_t>(positions, frontSeqSize);
	//std::cout << "__LINE__: " << __LINE__ << std::endl;
	auto readPos = gen.unifRandSelection(positions);
	//std::cout << "__LINE__: " << __LINE__ << std::endl;
	auto readLen = uiNormDist();
	//std::cout << "__LINE__: " << __LINE__ << std::endl;
	if (readLen + readPos < seqLen) {
		readPositionsLens.emplace_back(readPos, readLen);
	} else {
		readPositionsLens.emplace_back(readPos, seqLen - readPos);
	}
	//std::cout << "__LINE__: " << __LINE__ << std::endl;
	/**@todo need to fix below, currently shearing is using the fornt and end of the input seq twice */
	//forward
	{
		if (readPos + readLen < seqLen) {
			auto nextReadPos = readPos + readLen;
			auto nextReadLen = uiNormDist();
			while (nextReadPos + nextReadLen < seqLen) {
				readPositionsLens.emplace_back(nextReadPos, nextReadLen);
				nextReadPos = nextReadPos + nextReadLen;
				nextReadLen = uiNormDist();
			}
			readPositionsLens.emplace_back(nextReadPos, seqLen - nextReadPos);
		}
	}
	//std::cout << "__LINE__: " << __LINE__ << std::endl;
	//reverse
	{
		if (readPos > 0) {
			auto nextReadLen = uiNormDist();
			auto nextReadPos = readPos;
			while (nextReadPos > nextReadLen) {
				nextReadPos = nextReadPos - nextReadLen;
				readPositionsLens.emplace_back(nextReadPos, nextReadLen);
				nextReadLen = uiNormDist();
			}
			readPositionsLens.emplace_back(0, nextReadPos);
		}
	}
	//std::cout << "__LINE__: " << __LINE__ << std::endl;
	return readPositionsLens;
}

VecStr genFragments(const std::string & seq,
		const std::vector<std::pair<size_t, uint32_t>> & fragPosLens) {
	VecStr fragments;
	fragments.reserve(fragPosLens.size());
	for (const auto & posLen : fragPosLens) {
		fragments.emplace_back(seq.substr(posLen.first, posLen.second));
	}
	return fragments;
}

std::unordered_map<uint32_t, uint32_t> genMutCounts(std::mt19937_64 & gen,
		const uint64_t seqNumber,
		const uint64_t seqSize,
		const uint64_t intErrorRate) {
	std::unordered_map<uint32_t, uint32_t> currentMuts;
	for (uint32_t read = 0; read < seqNumber; ++read) {
		uint32_t count = 0;
		for (uint32_t base = 0; base < seqSize; ++base) {
			if (gen() <= intErrorRate) {
				++count;
			}
		}
		if (count != 0) {
			++currentMuts[count];
		}
	}
	return currentMuts;
}

std::vector<uint64_t> spreadReadNumAcrossThreads(const uint64_t readTotal,
		const uint32_t numThreads) {
	std::vector<uint64_t> tempAmounts;
	uint64_t tempAmount = readTotal / numThreads;
	uint64_t sum = 0;
	for (uint32_t t = 0; t < numThreads; ++t) {
		sum += tempAmount;
		tempAmounts.emplace_back(tempAmount);
	}
	tempAmounts.back() += (readTotal % sum);
	return tempAmounts;
}

std::pair<std::string, std::string> mutateSeq(std::string name, std::string seq,
		uint32_t mutNum, const std::vector<uint64_t>& readPositions,
		njh::randomGenerator & gen,
		std::unordered_map<char, njh::randObjectGen<char, uint32_t>>& charGen) {
	auto seqPositons = gen.unifRandSelectionVec(readPositions, mutNum, false);
	for (auto seqPos : seqPositons) {
		char base = charGen.at(seq[seqPos]).genObj();
		name.append("_" + estd::to_string(seqPos) + ":" + base);
		seq[seqPos] = base;
	}
	return {name, seq};
}



uint64_t runPcr(uint64_t intErrorRate, uint32_t numThreads, uint32_t pcrRounds,
		std::string seq, uint64_t startingTemplate, std::string name,
		std::unordered_map<std::string, uint64_t> & seqMap,
		std::mutex & seqMapLock, bool printTime) {
	std::random_device rd;
	njh::randomGenerator gen;
	std::vector<std::mt19937_64> gens;
	for(uint32_t t = 0; t < numThreads; ++t){
		gens.emplace_back(std::mt19937_64(rd()));
	}
	std::unordered_map<char, njh::randObjectGen<char, uint32_t>> charGen;
	charGen.emplace('T', njh::randObjectGen<char, uint32_t>(std::vector<char>{'A', 'C', 'G'}, std::vector<uint32_t>{1,80,1}));
	charGen.emplace('C', njh::randObjectGen<char, uint32_t>(std::vector<char>{'A', 'T', 'G'}, std::vector<uint32_t>{1,80,1}));
	charGen.emplace('G', njh::randObjectGen<char, uint32_t>(std::vector<char>{'C', 'A', 'T'}, std::vector<uint32_t>{1,80,1}));
	charGen.emplace('A', njh::randObjectGen<char, uint32_t>(std::vector<char>{'C', 'G', 'T'}, std::vector<uint32_t>{1,80,1}));

	std::vector<uint64_t> readPositions(seq.length());
	njh::iota<uint64_t>(readPositions, 0);
	for (uint32_t round = 1; round <= pcrRounds; ++round) {
		njh::stopWatch watch;
		uint32_t duplicatingAmount = startingTemplate;
		if(printTime){
			std::cout << "PCR Round " <<  round << std::endl;;
			std::cout << "\tDuplicatingAmount: " << duplicatingAmount << std::endl;
		}
		std::unordered_map<uint32_t, uint32_t> muts;
		std::mutex mutsLock;
		uint32_t numberMutated = 0;

		if(duplicatingAmount > numThreads){
			std::vector<uint64_t> tempAmounts = spreadReadNumAcrossThreads(duplicatingAmount, numThreads);
			auto mutate = [&mutsLock,&muts,&intErrorRate,&tempAmounts,&seq,&gens,&numberMutated](uint32_t threadNum) {
				auto currentMuts = genMutCounts(gens[threadNum], tempAmounts[threadNum], seq.size(), intErrorRate);
				{
					std::lock_guard<std::mutex> lock(mutsLock);
					for(const auto & mut : currentMuts){
						muts[mut.first] += mut.second;
						numberMutated += mut.second;
					}
				}
			};
			if(numThreads <=1){
				mutate(0);
			}else{
				std::vector<std::thread> threads;
				for (uint32_t thread = 0; thread < numThreads; ++thread) {
					threads.emplace_back(std::thread(mutate,thread));
				}
				for(auto & t : threads){
					t.join();
				}
			}
		} else {
			auto currentMuts = genMutCounts(gens.front(), duplicatingAmount, seq.size(),intErrorRate);
			for(const auto & mut : currentMuts){
				muts[mut.first] += mut.second;
				numberMutated += mut.second;
			}
		}
		if(printTime){
			std::cout << "\tmutation\tnumber" << std::endl;
			for(const auto & mutCount : muts){
				std::cout << "\t" << mutCount.first << "\t\t" <<mutCount.second << std::endl;
			}
		}

		//add on the number of reads that weren't mutated;
		startingTemplate += duplicatingAmount - numberMutated;
		std::vector<std::pair<std::string, std::string>> mutants;
		for(const auto & mut : muts){
			for(uint32_t num = 0; num < mut.second; ++num){
				auto mutSeq = mutateSeq(name,seq, mut.first, readPositions, gen, charGen);
				mutants.emplace_back(mutSeq);
			}
		}
		if(round == pcrRounds){
			std::lock_guard<std::mutex> outMutLock(seqMapLock);
			seqMap[seq] += startingTemplate;
			for(const auto & mutSeq : mutants){
				seqMap[mutSeq.second] += 1;
			}
		}else{
			if (mutants.size() > 1 && numThreads > 1) {
				auto runNextPcr =
						[&mutants,&seqMap,&seqMapLock,&intErrorRate](uint32_t pcrRound, std::vector<uint32_t> mutPositions) {
							for(const auto & mutPos : mutPositions) {
								runPcr(intErrorRate, 1, pcrRound, mutants[mutPos].second, 1, mutants[mutPos].first, seqMap, seqMapLock,false);
							}
						};
				std::unordered_map<uint32_t, std::vector<uint32_t>> mutantVecs;
				for (const auto & mutNamePos : iter::range(mutants.size())) {
					mutantVecs[mutNamePos % numThreads].emplace_back(mutNamePos);
				}
				std::vector<std::thread> threads;
				for (uint32_t threadNum = 0; threadNum < numThreads; ++threadNum) {
					if (njh::in(threadNum, mutantVecs)) {
						threads.emplace_back(
								std::thread(runNextPcr, pcrRounds - round,
										mutantVecs.at(threadNum)));
					}
				}
				for(auto & thread : threads){
					thread.join();
				}
			} else {
				for (const auto & mut : mutants) {
					runPcr(intErrorRate, 1, pcrRounds - round,mut.second, 1,
							mut.first, seqMap, seqMapLock, false);
				}
			}
		}
		if(printTime){
			std::cout << "\tRound " << round << " time: " << watch.totalTimeFormatted(6) << std::endl << std::endl;
		}
	}
	return startingTemplate;
}





std::pair<uint64_t,uint64_t> sampleReadsWithReplacement(const std::string & name, const std::string & seq,
		std::unordered_map<std::string, uint64_t> & seqCounts,
		uint64_t finalPerfectAmount, uint64_t finalReadAmount,
		std::ostream & sequenceOutFile, std::mutex & seqFileLock,
		uint32_t numThreads) {

	uint64_t mutatedOut = 0;
	uint64_t nonMutatedOut = 0;
	double finalPerfectAmountDbl = finalPerfectAmount;
	std::cout << "Sampling " << name << "..." << std::endl;
	auto sampleReads =
			[&seqFileLock,&sequenceOutFile,&seq,&seqCounts,&finalPerfectAmountDbl,&mutatedOut,&nonMutatedOut,&name](uint32_t sampleAmount) {
				std::random_device rd;
				std::mt19937_64 mtGen(rd());
				for(uint32_t read = 0; read < sampleAmount; ++read) {
					uint64_t randNum = mtGen();
					uint64_t sum = 0;
					for(const auto & seqCount : seqCounts) {
						sum += (seqCount.second/finalPerfectAmountDbl) * mtGen.max();
						if(sum >= randNum) {
							if(seqCount.first != seq) {
								std::lock_guard<std::mutex> fileLock(seqFileLock);
								++mutatedOut;
								sequenceOutFile << ">" << name << "_mut." << leftPadNumStr<uint64_t>(mutatedOut, finalPerfectAmountDbl) << std::endl;
								sequenceOutFile << seqCount.first << std::endl;
							} else {
								std::lock_guard<std::mutex> fileLock(seqFileLock);
								++nonMutatedOut;
								sequenceOutFile << ">" << name << "_seq." << leftPadNumStr<uint64_t>(nonMutatedOut, finalPerfectAmountDbl) << std::endl;
								sequenceOutFile << seqCount.first << std::endl;
							}
							break;
						}
					}
				}
			};
	std::vector<uint64_t> readAmounts;
	uint64_t threadReadAmount = finalReadAmount / numThreads;
	uint64_t threadReadSum = 0;
	for (uint32_t threadNum = 0; threadNum < numThreads; ++threadNum) {
		readAmounts.emplace_back(threadReadAmount);
		threadReadSum += threadReadAmount;
	}
	readAmounts.back() += finalReadAmount % threadReadSum;
	std::vector<std::thread> threads;
	for (const auto & amount : readAmounts) {
		threads.emplace_back(std::thread(sampleReads, amount));
	}
	for (auto & thread : threads) {
		thread.join();
	}
	return {nonMutatedOut, mutatedOut};
}

std::pair<uint64_t,uint64_t> sampleReadsWithoutReplacement(const std::string & name,
		const std::string & seq,
		std::unordered_map<std::string, uint64_t> & seqCounts,
		uint64_t finalReadAmount,
		std::ostream & sequenceOutFile, std::mutex & seqFileLock,
		bool verbose) {
	uint64_t mutatedOut = 0;
	uint64_t nonMutatedOut = 0;
	double totalTemplate = 0;
	for (const auto & seqCount : seqCounts) {
		totalTemplate += seqCount.second;
	}
	if(verbose){
		std::cout << "Sampling " << name << "..." << std::endl;
	}

	std::random_device rd;
	std::mt19937_64 mtGen(rd());
	for (uint32_t read = 0; read < finalReadAmount; ++read) {
		uint64_t randNumGen = mtGen();
		uint64_t randSel = (static_cast<double>(randNumGen)/mtGen.max()) * (totalTemplate - mutatedOut - nonMutatedOut);
		uint64_t sum = 0;
		for (auto & seqCount : seqCounts) {
			if (seqCount.second == 0) {
				continue;
			}
			sum += seqCount.second;
			if (sum >= randSel) {
				{
					std::lock_guard<std::mutex> fileLock(seqFileLock);
					--seqCount.second;
					if (seqCount.first != seq) {
						++mutatedOut;
						sequenceOutFile << ">" << name << "_mut."
								<< leftPadNumStr<uint64_t>(mutatedOut, finalReadAmount)
								<< std::endl;
						sequenceOutFile << seqCount.first << std::endl;
					} else {
						++nonMutatedOut;
						sequenceOutFile << ">" << name << "_seq."
								<< leftPadNumStr<uint64_t>(nonMutatedOut, finalReadAmount)
								<< std::endl;
						sequenceOutFile << seqCount.first << std::endl;
					}
				}
				break;
			}
		}
	}
	return {nonMutatedOut, mutatedOut};
}

std::unordered_map<std::string, std::pair<uint64_t, uint64_t>> sampleReadsWithoutReplacement(
		const std::unordered_map<std::string, std::string> & seqs,
		std::unordered_map<std::string, std::unordered_map < std::string, uint64_t>> & multipleSeqCounts,
		uint64_t finalReadAmount, std::ostream & sequenceOutFile,
		std::mutex & seqFileLock, bool verbose) {
	uint64_t mutatedOut = 0;
	uint64_t nonMutatedOut = 0;
	std::unordered_map<std::string, std::pair<uint64_t, uint64_t>> ret;

	auto seqNames = getVectorOfMapKeys(seqs);
	auto multipleSeqCountsNames = getVectorOfMapKeys(multipleSeqCounts);
	njh::sort(seqNames);
	njh::sort(multipleSeqCountsNames);
	if(seqNames != multipleSeqCountsNames){
		std::stringstream ss;
		ss << "Error in " << __PRETTY_FUNCTION__ << " names for seqNames and multipleSeqCountsNames counts don't match" << std::endl;
		ss << "seqNames: " << njh::conToStr(seqNames,",") << std::endl;
		ss << "multipleSeqCountsNames: " << njh::conToStr(multipleSeqCountsNames,",") << std::endl;
		throw std::runtime_error{ss.str()};
	}
	uint64_t totalTemplate = 0;
	for (const auto & seqCounts : multipleSeqCounts) {
		for (const auto & seqCount : seqCounts.second) {
			totalTemplate += seqCount.second;
		}
	}

	if(verbose){
		auto names = getVectorOfMapKeys(seqs);
		std::cout << "Sampling from " << njh::conToStr(names,",") << " ..." << std::endl;
	}

	std::random_device rd;
	std::mt19937_64 mtGen(rd());
	for (uint32_t read = 0; read < finalReadAmount; ++read) {
		uint64_t randNumGen = mtGen();
		uint64_t randSel = (static_cast<double>(randNumGen)/mtGen.max()) * (totalTemplate - mutatedOut - nonMutatedOut);
		uint64_t sum = 0;
		for(auto & seqCounts : multipleSeqCounts){
			bool foundSelection = false;
			for (auto & seqCount : seqCounts.second) {
				if (seqCount.second == 0) {
					continue;
				}
				sum += seqCount.second;
				if (sum >= randSel) {
					{
						std::lock_guard<std::mutex> fileLock(seqFileLock);
						--seqCount.second;
						if (seqCount.first != seqs.at(seqCounts.first)) {
							++mutatedOut;
							++ret[seqCounts.first].second;
							sequenceOutFile << ">" << seqCounts.first << "_mut."
									<< leftPadNumStr<uint64_t>(mutatedOut, finalReadAmount)
									<< std::endl;
							sequenceOutFile << seqCount.first << std::endl;
						} else {
							++nonMutatedOut;
							++ret[seqCounts.first].first;
							sequenceOutFile << ">" << seqCounts.first << "_seq."
									<< leftPadNumStr<uint64_t>(nonMutatedOut, finalReadAmount)
									<< std::endl;
							sequenceOutFile << seqCount.first << std::endl;
						}
					}
					foundSelection = true;
					break;
				}
			}
			if(foundSelection){
				break;
			}
		}
	}
	return ret;
}

std::string runPcrSingleTemplate(
		uint64_t intErrorRate,
		uint32_t roundsOfPcr,
		uint64_t randomNumberSelector,
		uint64_t randomNumberSelectorMax,
		std::string seq ){
	std::unordered_map<std::string, uint64_t> finalProducts;
	std::mutex finalProductsMutex;
	runPcr(intErrorRate, 1,
			roundsOfPcr, seq, 1, "singleTemplate",
			finalProducts, finalProductsMutex, false);
	uint64_t finalRandSel = (static_cast<double>(randomNumberSelector)/randomNumberSelectorMax) * (std::pow(2,roundsOfPcr));
	uint64_t finalSum = 0;
	for(const auto & final : finalProducts) {
		finalSum += final.second;
		if(finalSum >= finalRandSel) {
			seq = final.first;
			break;
		}
	}
	return seq;
}

std::unordered_map<std::string, std::pair<uint64_t, uint64_t>> sampleReadsWithoutReplacementFinishPCR(
		const std::unordered_map<std::string, std::string> & seqs,
		std::unordered_map<std::string, std::unordered_map<std::string, uint64_t>> & multipleSeqCounts,
		uint64_t finalReadAmount, std::ostream & sequenceOutFile,
		std::mutex & seqFileLock, uint32_t numberOfPCRRoundsLeft,
		uint64_t intErrorRate, uint32_t numThreads, bool verbose) {

	std::unordered_map<std::string, std::pair<uint64_t, uint64_t>> ret;
	auto seqNames = getVectorOfMapKeys(seqs);
	auto multipleSeqCountsNames = getVectorOfMapKeys(multipleSeqCounts);
	njh::sort(seqNames);
	njh::sort(multipleSeqCountsNames);
	if(seqNames != multipleSeqCountsNames){
		std::stringstream ss;
		ss << "Error in " << __PRETTY_FUNCTION__ << " names for seqNames and multipleSeqCountsNames counts don't match" << std::endl;
		ss << "seqNames: " << njh::conToStr(seqNames,",") << std::endl;
		ss << "multipleSeqCountsNames: " << njh::conToStr(multipleSeqCountsNames,",") << std::endl;
		throw std::runtime_error{ss.str()};
	}
	uint64_t totalTemplate = 0;
	for (const auto & seqCounts : multipleSeqCounts) {
		for (const auto & seqCount : seqCounts.second) {
			totalTemplate += seqCount.second;
		}
	}

	if(verbose){
		auto names = getVectorOfMapKeys(seqs);
		//std::cout << "Sampling from " << njh::conToStr(names,",") << " ..." << std::endl;
	}

	std::random_device rd;
	std::mt19937_64 mtGen(rd());
	std::unordered_map<uint32_t, std::vector<std::pair<std::string, std::string>>> allSampledSeqs;
	if(verbose){
		std::cout << "Initial Sampling" << std::endl;
	}
	for (uint32_t read = 0; read < finalReadAmount; ++read) {
		if(verbose){
			std::cout << '\r' << read + 1;
		}
		uint64_t randNumGen = mtGen();
		uint64_t randSel = (static_cast<double>(randNumGen)/mtGen.max()) * (totalTemplate - read);
		uint64_t sum = 0;
		for(auto & seqCounts : multipleSeqCounts){
			bool foundSelection = false;
			for (auto & seqCount : seqCounts.second) {
				if (seqCount.second == 0) {
					continue;
				}
				sum += seqCount.second;
				if (sum >= randSel) {
					--seqCount.second;
					allSampledSeqs[read % numThreads].emplace_back(std::pair<std::string,std::string>{seqCounts.first,seqCount.first});
					foundSelection = true;
					break;
				}
			}
			if(foundSelection){
				break;
			}
		}
	}
	if(verbose){
		std::cout << std::endl;;
	}
	auto finishPCR =
			[&allSampledSeqs,&intErrorRate,&numberOfPCRRoundsLeft,&seqs,&seqFileLock,&sequenceOutFile,&ret,&finalReadAmount,&verbose](uint32_t threadNumber ) {
				njh::stopWatch watch;
				std::vector<std::pair<std::string, std::string>> outputs;
				std::unordered_map<std::string, std::pair<uint64_t, uint64_t>> mutInfo;
				std::random_device rd;
				std::mt19937_64 mtGenFinal(rd());
				const uint32_t significantRounds = 10;
				for(const auto & namesSeqs : allSampledSeqs[threadNumber]) {
					std::string finalSeq = namesSeqs.second;
					uint32_t roundsStillLeft = numberOfPCRRoundsLeft;
					while(roundsStillLeft > significantRounds) {
						finalSeq = runPcrSingleTemplate(intErrorRate,
								significantRounds,mtGenFinal(), mtGenFinal.max(),
								finalSeq);
						roundsStillLeft -= significantRounds;
					}
					finalSeq = runPcrSingleTemplate(intErrorRate,
							roundsStillLeft,mtGenFinal(),mtGenFinal.max(),
							finalSeq);
					std::stringstream nameStream;
					MetaDataInName nameMeta;
					nameMeta.addMeta("hap", namesSeqs.first);

					if (finalSeq != seqs.at(namesSeqs.first)) {
						nameMeta.addMeta("mutated", true);
						nameMeta.addMeta("idNum", leftPadNumStr<uint64_t>(mutInfo[namesSeqs.first].second, finalReadAmount));
						++mutInfo[namesSeqs.first].second;
//						nameStream << ">" << namesSeqs.first << "_mut."
//						<< leftPadNumStr<uint64_t>(mutInfo[namesSeqs.first].second, finalReadAmount)
//						<< "_" << threadNumber;
					} else {
						++mutInfo[namesSeqs.first].first;
						nameMeta.addMeta("mutated", false);
						nameMeta.addMeta("idNum", leftPadNumStr<uint64_t>(mutInfo[namesSeqs.first].first, finalReadAmount));

//						nameStream << ">" << namesSeqs.first << "_seq."
//						<< leftPadNumStr<uint64_t>(mutInfo[namesSeqs.first].first, finalReadAmount)
//						<< "_" << threadNumber;
					}
//					outputs.emplace_back(std::make_pair(nameStream.str(), finalSeq));
					outputs.emplace_back(std::make_pair(">" + nameMeta.createMetaName(), finalSeq));

				}
				{
					std::lock_guard<std::mutex> fileLock(seqFileLock);
					if(verbose) {
						std::cout << "Thread " << threadNumber << " done" << std::endl;
						std::cout << "\tFinished PCR for "
								<< allSampledSeqs[threadNumber].size() << " reads" << std::endl;
						std::cout << "\tTime: " << watch.totalTimeFormatted(6) << std::endl;
					}
					for(const auto & out : outputs){
						sequenceOutFile << out.first << std::endl;
						sequenceOutFile << out.second << std::endl;
					}
					for(const auto & info : mutInfo){
						ret[info.first].first += info.second.first;
						ret[info.first].second += info.second.second;
					}
				}
			};

	if(numThreads <=1){
		finishPCR(0);
	}else{
		std::vector<std::thread> threads;
		for(uint32_t thread = 0; thread < numThreads; ++thread){
			threads.emplace_back(std::thread(finishPCR, thread));
		}

		for(auto & thread : threads){
			thread.join();
		}
	}
	return ret;
}

std::pair<uint64_t, uint64_t> sampleReadsWithoutReplacementFinishPCR(
		std::unordered_map<std::string, uint64_t> & seqCounts,
		uint64_t finalReadAmount, std::ostream & sequenceOutFile,
		std::mutex & seqFileLock, uint32_t numberOfPCRRoundsLeft,
		uint64_t intErrorRate, uint32_t numThreads, bool verbose) {
	std::pair<uint64_t, uint64_t> ret;
	uint64_t totalTemplate = 0;
	for (const auto & seqCount : seqCounts) {
		totalTemplate += seqCount.second;
	}
	std::random_device rd;
	std::mt19937_64 mtGen(rd());
	std::unordered_map<uint32_t, VecStr> allSampledSeqs;
	if(verbose){
		std::cout << "Initial Sampling from " << totalTemplate << " total template reads" << std::endl;
		std::cout << "Initial Sampling from " << seqCounts.size() << " total unique reads" << std::endl;
	}
	{
		std::vector<decltype(seqCounts.begin())> seqs;
		for(auto it = seqCounts.begin(); it != seqCounts.end(); ++it){
			if(it != seqCounts.end()){
				for(uint64_t i  = 0; i < it->second; ++i){
					seqs.push_back(it);
				}
			}
		}
		for (uint32_t read = 0; read < finalReadAmount; ++read) {
			uint64_t randSel = std::round(static_cast<double>(mtGen())/mtGen.max() * (totalTemplate - 1));
			if(verbose){
				std::cout << '\r' << read + 1 << " : " << randSel << " - " << seqs.size();
				std::cout.flush();
			}
			while(seqs[randSel]->second < 1){
				randSel = std::round(static_cast<double>(mtGen())/mtGen.max() * (totalTemplate - 1));
			}
			--(seqs[randSel]->second);
			allSampledSeqs[read % numThreads].emplace_back(seqs[randSel]->first);

			/*
			uint64_t randSel = (static_cast<double>(randNumGen)/mtGen.max()) * (totalTemplate - read);
			uint64_t sum = 0;
			for (auto & seqCount : seqCounts) {
				if (seqCount.second == 0) {
					continue;
				}
				sum += seqCount.second;
				if (sum >= randSel) {
					--seqCount.second;
					allSampledSeqs[read % numThreads].emplace_back(seqCount.first);
					break;
				}
			}*/
		}
		seqs.clear();
		seqs.resize(0);
	}

	if(verbose){
		std::cout << std::endl;;
	}
	auto finishPCR =
			[&allSampledSeqs,&intErrorRate,&numberOfPCRRoundsLeft,&seqFileLock,&sequenceOutFile,&ret,&finalReadAmount,&verbose](uint32_t threadNumber ) {
				njh::stopWatch watch;
				std::vector<std::pair<std::string, std::string>> outputs;
				std::pair<uint64_t, uint64_t> mutInfo;
				std::random_device rd;
				std::mt19937_64 mtGenFinal(rd());
				const uint32_t significantRounds = 10;
				for(const auto & seq : allSampledSeqs[threadNumber]) {
					std::string finalSeq = seq;
					uint32_t roundsStillLeft = numberOfPCRRoundsLeft;
					while(roundsStillLeft > significantRounds) {
						finalSeq = runPcrSingleTemplate(intErrorRate,
								significantRounds,mtGenFinal(), mtGenFinal.max(),
								finalSeq);
						roundsStillLeft -= significantRounds;
					}
					finalSeq = runPcrSingleTemplate(intErrorRate,
							roundsStillLeft,mtGenFinal(),mtGenFinal.max(),
							finalSeq);
					std::stringstream nameStream;
					if (finalSeq != seq) {
						++mutInfo.second;
						nameStream << ">" << seq.substr(0,8)<< "_mut."
						<< leftPadNumStr<uint64_t>(mutInfo.second, finalReadAmount)
						<< "_" << threadNumber;
					} else {
						++mutInfo.first;
						nameStream << ">" << seq.substr(0,8) << "_seq."
						<< leftPadNumStr<uint64_t>(mutInfo.first, finalReadAmount)
						<< "_" << threadNumber;
					}
					outputs.emplace_back(std::make_pair(nameStream.str(), finalSeq));
				}
				{
					std::lock_guard<std::mutex> fileLock(seqFileLock);
					if(verbose) {
						std::cout << "Thread " << threadNumber << " done" << std::endl;
						std::cout << "\tFinished PCR for "
								<< allSampledSeqs[threadNumber].size() << " reads" << std::endl;
						std::cout << "\tTime: " << watch.totalTimeFormatted(6) << std::endl;
					}
					for(const auto & out : outputs){
						sequenceOutFile << out.first << std::endl;
						sequenceOutFile << out.second << std::endl;
					}
					ret.first += mutInfo.first;
					ret.second += mutInfo.second;
				}
			};

	if(numThreads <=1){
		finishPCR(0);
	}else{
		std::vector<std::thread> threads;
		for(uint32_t thread = 0; thread < numThreads; ++thread){
			threads.emplace_back(std::thread(finishPCR, thread));
		}

		for(auto & thread : threads){
			thread.join();
		}
	}
	return ret;
}

void simLib(std::vector<std::shared_ptr<seqInfo> > & reads, const std::string & barcode,
		const std::string & workingDir, const std::string & libName,
		uint64_t intErrorRate, uint64_t startingTemplate, uint64_t finalReadAmount,
		uint32_t pcrRounds, uint32_t numThreads, bool verbose) {
	//auto finalPerfectAmount = static_cast<uint64_t>(startingTemplate * std::pow(2, pcrRounds));
	auto libDirName = njh::files::makeDir(workingDir, njh::files::MkdirPar(libName));
	std::ofstream libOutFile;
	openTextFile(libOutFile, OutOptions(bfs::path(libDirName.string() + "reads.fasta")));
	std::mutex seqFileLock;
	//check there is enough template and final read amount to make reads for the desired fractions
	std::unordered_map<std::string, uint64_t> templateAmountCounts;
	std::unordered_map<std::string, uint64_t> finalReadAmountCounts;
	uint64_t maxTemplateReads = 0;
	uint64_t readTemplateSum = 0;
	std::string maxTemplateRead = "";
	uint64_t maxFinalReads = 0;
	uint64_t readFinalSum = 0;
	std::string maxFinalRead = "";
	for(const auto & read : reads){
		uint64_t currentTemplateAmt = std::round(read->frac_ * startingTemplate);
		uint64_t currentFinalAmt = std::round(read->frac_ * finalReadAmount);
		if(currentTemplateAmt == 0){
			std::stringstream ss;
			ss << "Error in: " << __PRETTY_FUNCTION__
					<< ", not enough starting template reads, " << startingTemplate
					<< ", requested in order to simulate the desired read fraction:"
					<< read->frac_ << std::endl;
			throw std::runtime_error{njh::bashCT::boldRed(ss.str())};
		}
		if (currentFinalAmt == 0) {
			std::stringstream ss;
			ss << "Error in: " << __PRETTY_FUNCTION__
					<< ", not enough final reads, " << finalReadAmount
					<< ", requested in order to simulate the desired read fraction:"
					<< read->frac_ << std::endl;
			throw std::runtime_error{njh::bashCT::boldRed(ss.str())};
		}
		templateAmountCounts[read->name_] = currentTemplateAmt;
		finalReadAmountCounts[read->name_] = currentFinalAmt;
		readTemplateSum += currentTemplateAmt;
		readFinalSum += currentFinalAmt;
		if(currentTemplateAmt > maxTemplateReads){
			maxTemplateReads = currentTemplateAmt;
			maxTemplateRead = read->name_;
		}
		if(currentFinalAmt > maxFinalReads){
			maxFinalReads = currentFinalAmt;
			maxFinalRead = read->name_;
		}
	}
	if(readTemplateSum < startingTemplate){
		templateAmountCounts[maxTemplateRead] += startingTemplate - readTemplateSum;
	}

	if(readFinalSum < finalReadAmount){
		finalReadAmountCounts[maxFinalRead] += finalReadAmount - readFinalSum;
	}

	if(verbose){
		std::cout << "Simulating library " << libName << std::endl;
		std::cout << "Starting total template: " << startingTemplate << std::endl;
		for(const auto & tempAmount : templateAmountCounts){
			std::cout << tempAmount.first << ": " << tempAmount.second << std::endl;
		}
	}
	std::unordered_map<std::string,std::unordered_map<std::string, uint64_t>> allSeqCounts;
	std::unordered_map<std::string,std::string> barcodedSeqs;
	std::unordered_map<std::string,std::pair<uint64_t,uint64_t>> templateNonMutated;
	for (const auto & read : reads) {
		std::unordered_map<std::string, uint64_t> seqCounts;
		std::mutex seqMapLock;
		if(verbose){
			std::cout << "Simulating ref read " << read->name_ << " for library " << libName << std::endl;
			std::cout << "\tTemplate Amount:          " << templateAmountCounts[read->name_] << std::endl;
			std::cout << "\tFinal Read Sample amount: " << finalReadAmountCounts[read->name_] << std::endl;
			std::cout << "\tBarcode: " << barcode << std::endl;
		}
		std::string barcodedSeq = barcode + read->seq_;
		auto finalAmount = runPcr(intErrorRate, numThreads, pcrRounds, barcodedSeq,
				templateAmountCounts[read->name_], read->name_, seqCounts, seqMapLock,
				verbose);
		uint64_t finalPerfectAmount = templateAmountCounts[read->name_] * std::pow(2, pcrRounds);
		templateNonMutated[read->name_] = {finalAmount, finalPerfectAmount};
		barcodedSeqs[read->name_] = barcodedSeq;
		allSeqCounts[read->name_] = seqCounts;
	}
	auto sampleNumber = sampleReadsWithoutReplacement(barcodedSeqs,
			allSeqCounts, finalReadAmount,
			libOutFile, seqFileLock, verbose);
	if(verbose){
		{
			std::cout << "PCR amounts: " << std::endl;
			uint64_t nonMutated = 0;
			uint64_t mutated = 0;
			for (const auto & read : reads) {
				std::cout << read->name_ << std::endl;
				nonMutated += templateNonMutated[read->name_].first;
				mutated += templateNonMutated[read->name_].second
						- templateNonMutated[read->name_].first;
				std::cout << "\t"
						<< getPercentageString(templateNonMutated[read->name_].first,
								templateNonMutated[read->name_].second) << std::endl;
			}
			std::cout << "Total: " << nonMutated + mutated << std::endl;
			std::cout << "\t" << getPercentageString(nonMutated, nonMutated + mutated)
					<< std::endl;
			std::cout << "\t" << getPercentageString(mutated, nonMutated + mutated)
					<< std::endl;
		}
		std::cout << "Sampling Amounts:" << std::endl;
		uint64_t nonMutated = 0;
		uint64_t mutated = 0;
		for(const auto & read : reads){
			std::cout << read->name_ << std::endl;
			auto total = sampleNumber[read->name_].first + sampleNumber[read->name_].second;
			nonMutated+= sampleNumber[read->name_].first;
			mutated+= sampleNumber[read->name_].second;
			std::cout << "\tSampled     : " << getPercentageString(total, finalReadAmountCounts[read->name_])<< std::endl;
			std::cout << "\tNon-Mutated : " << getPercentageString(sampleNumber[read->name_].first, total) << std::endl;
			std::cout << "\tMutated     : " << getPercentageString(sampleNumber[read->name_].second, total) << std::endl;
		}
		std::cout << "Total\t       : " << nonMutated + mutated << std::endl;
		std::cout << "\tNon-Mutated : " << getPercentageString(nonMutated, nonMutated + mutated) << std::endl;
		std::cout << "\tMutated     : " << getPercentageString(mutated, nonMutated + mutated) << std::endl;

	}
	libOutFile.close();
}

void simLibFast(std::vector<std::shared_ptr<seqInfo> > & reads,
		const std::string & barcode, const std::string & workingDir,
		const std::string & libName, uint64_t intErrorRate,
		uint64_t startingTemplate, uint64_t finalReadAmount, uint32_t pcrRounds,
		uint32_t initialPcrRounds, uint32_t numThreads, bool verbose) {


	if(initialPcrRounds >= pcrRounds){
		std::stringstream ss;
		ss << "Error in " << __PRETTY_FUNCTION__ << std::endl;
		ss << "initialPcrRounds should be less than pcrRounds" << std::endl;
		ss << "pcrRounds:" << pcrRounds << std::endl;
		ss << "initialPcrRounds:" << initialPcrRounds << std::endl;
		throw std::runtime_error{ss.str()};
	}
	//auto finalPerfectAmount = static_cast<uint64_t>(startingTemplate * std::pow(2, pcrRounds));
	auto libDirName = njh::files::makeDir(workingDir, njh::files::MkdirPar(libName));
	std::ofstream libOutFile;
	openTextFile(libOutFile, OutOptions(bfs::path(libDirName.string() + "reads.fasta")));
	std::mutex seqFileLock;
	//check there is enough template and final read amount to make reads for the desired fractions
	std::unordered_map<std::string, uint64_t> templateAmountCounts;
	std::unordered_map<std::string, uint64_t> finalReadAmountCounts;
	uint64_t maxTemplateReads = 0;
	uint64_t readTemplateSum = 0;
	std::string maxTemplateRead = "";
	uint64_t maxFinalReads = 0;
	uint64_t readFinalSum = 0;
	std::string maxFinalRead = "";
	for(const auto & read : reads){
		uint64_t currentTemplateAmt = std::round(read->frac_ * startingTemplate);
		uint64_t currentFinalAmt = std::round(read->frac_ * finalReadAmount);
		if(currentTemplateAmt == 0){
			std::stringstream ss;
			ss << "Error in: " << __PRETTY_FUNCTION__
					<< ", not enough starting template reads, " << startingTemplate
					<< ", requested in order to simulate the desired read fraction:"
					<< read->frac_ << std::endl;
			throw std::runtime_error{njh::bashCT::boldRed(ss.str())};
		}
		if (currentFinalAmt == 0) {
			std::stringstream ss;
			ss << "Error in: " << __PRETTY_FUNCTION__
					<< ", not enough final reads, " << finalReadAmount
					<< ", requested in order to simulate the desired read fraction:"
					<< read->frac_ << std::endl;
			throw std::runtime_error{njh::bashCT::boldRed(ss.str())};
		}
		templateAmountCounts[read->name_] = currentTemplateAmt;
		finalReadAmountCounts[read->name_] = currentFinalAmt;
		readTemplateSum += currentTemplateAmt;
		readFinalSum += currentFinalAmt;
		if(currentTemplateAmt > maxTemplateReads){
			maxTemplateReads = currentTemplateAmt;
			maxTemplateRead = read->name_;
		}
		if(currentFinalAmt > maxFinalReads){
			maxFinalReads = currentFinalAmt;
			maxFinalRead = read->name_;
		}
	}
	if(readTemplateSum < startingTemplate){
		templateAmountCounts[maxTemplateRead] += startingTemplate - readTemplateSum;
	}

	if(readFinalSum < finalReadAmount){
		finalReadAmountCounts[maxFinalRead] += finalReadAmount - readFinalSum;
	}

	if(verbose){
		std::cout << "Simulating library " << libName << std::endl;
		std::cout << "Starting total template: " << startingTemplate << std::endl;
		for(const auto & tempAmount : templateAmountCounts){
			std::cout << tempAmount.first << ": " << tempAmount.second << std::endl;
		}
	}
	std::unordered_map<std::string,std::unordered_map<std::string, uint64_t>> allSeqCounts;
	std::unordered_map<std::string,std::string> barcodedSeqs;
	std::unordered_map<std::string,std::pair<uint64_t,uint64_t>> templateNonMutated;
	for (const auto & read : reads) {
		std::unordered_map<std::string, uint64_t> seqCounts;
		std::mutex seqMapLock;
		if(verbose){
			std::cout << "Simulating ref read " << read->name_ << " for library " << libName << std::endl;
			std::cout << "\tTemplate Amount:          " << templateAmountCounts[read->name_] << std::endl;
			std::cout << "\tFinal Read Sample amount: " << finalReadAmountCounts[read->name_] << std::endl;
			std::cout << "\tBarcode: " << barcode << std::endl;
		}
		std::string barcodedSeq = barcode + read->seq_;
		auto finalAmount = runPcr(intErrorRate, numThreads, initialPcrRounds, barcodedSeq,
				templateAmountCounts[read->name_], read->name_, seqCounts, seqMapLock,
				verbose);

		uint64_t finalPerfectAmount = templateAmountCounts[read->name_] * std::pow(2, initialPcrRounds);

		//minus off the starting template amount as this is just genomic DNA and won't be able to be sequenced

		if(seqCounts[barcodedSeq] <= templateAmountCounts[read->name_]){
			seqCounts.erase(barcodedSeq);
		}else{
			seqCounts[barcodedSeq] -= templateAmountCounts[read->name_];
		}
		if(!seqCounts.empty()){
			allSeqCounts[read->name_] = seqCounts;
			templateNonMutated[read->name_] = {finalAmount, finalPerfectAmount};
			barcodedSeqs[read->name_] = barcodedSeq;
		}
	}
	auto sampleNumber = sampleReadsWithoutReplacementFinishPCR(barcodedSeqs,
			allSeqCounts, finalReadAmount, libOutFile, seqFileLock,
			pcrRounds - initialPcrRounds, intErrorRate, numThreads, verbose);

	if(verbose){
		{
			std::cout << "PCR amounts: " << std::endl;
			uint64_t nonMutated = 0;
			uint64_t mutated = 0;
			for (const auto & read : reads) {
				std::cout << read->name_ << std::endl;
				nonMutated += templateNonMutated[read->name_].first;
				mutated += templateNonMutated[read->name_].second
						- templateNonMutated[read->name_].first;
				std::cout << "\t"
						<< getPercentageString(templateNonMutated[read->name_].first,
								templateNonMutated[read->name_].second) << std::endl;
			}
			std::cout << "Total: " << nonMutated + mutated << std::endl;
			std::cout << "\t" << getPercentageString(nonMutated, nonMutated + mutated)
					<< std::endl;
			std::cout << "\t" << getPercentageString(mutated, nonMutated + mutated)
					<< std::endl;
		}
		std::cout << "Sampling Amounts:" << std::endl;
		uint64_t nonMutated = 0;
		uint64_t mutated = 0;
		for(const auto & read : reads){
			std::cout << read->name_ << std::endl;
			auto total = sampleNumber[read->name_].first + sampleNumber[read->name_].second;
			nonMutated+= sampleNumber[read->name_].first;
			mutated+= sampleNumber[read->name_].second;
			std::cout << "\tSampled     : " << getPercentageString(total, finalReadAmountCounts[read->name_])<< std::endl;
			std::cout << "\tNon-Mutated : " << getPercentageString(sampleNumber[read->name_].first, total) << std::endl;
			std::cout << "\tMutated     : " << getPercentageString(sampleNumber[read->name_].second, total) << std::endl;
		}
		std::cout << "Total\t       : " << nonMutated + mutated << std::endl;
		std::cout << "\tNon-Mutated : " << getPercentageString(nonMutated, nonMutated + mutated) << std::endl;
		std::cout << "\tMutated     : " << getPercentageString(mutated, nonMutated + mutated) << std::endl;

	}
	libOutFile.close();
}

std::unordered_map<std::string, uint64_t> shearLibraryThreaded(
		std::unordered_map<std::string, std::unordered_map<std::string, uint64_t>> allSeqCounts,
		uint32_t mean, uint32_t std, uint32_t minLen, bool verbose) {
	std::unordered_map<std::string, uint64_t> ret;
	sim::ReadLenNormalDistribution uiNormDist(mean, std);
	njh::randomGenerator gen;
	if(verbose){
		std::cout << "Shearing library" << std::endl;
	}
	for (const auto & seqCount : allSeqCounts) {
		if(verbose){
			std::cout << "\tShearing " << seqCount.first << std::endl;
		}
		njh::stopWatch watch;
//		for (const auto & seq : seqCount.second) {
//			std::cout << seq.first << std::endl;
//		}
		for (const auto & seq : seqCount.second) {
//			std::cout << seq.first << std::endl;
			auto compSeq = seqUtil::reverseComplement(seq.first, "DNA");
			for (uint32_t i = 0; i < seq.second; ++i) {
//				std::cout << i << std::endl;
//				std::cout << seq.first.length() << std::endl;
				auto readPositionsLens = sim::genFragPosSizes(seq.first.length(),
						uiNormDist, gen);

				auto fragments = sim::genFragments(seq.first, readPositionsLens);
				auto compFragments = sim::genFragments(compSeq, readPositionsLens);
				for (const auto & frag : fragments) {
					if(frag.length() >=minLen){
						++ret[frag];
					}
				}
				for (const auto & frag : compFragments) {
					if(frag.length() >=minLen){
						++ret[frag];
					}
				}
			}
		}
		if(verbose){
			std::cout << "\t" << watch.totalTimeFormatted(6) << std::endl;
		}
	}
	return ret;
}



std::unordered_map<std::string, uint64_t> shearLibrary(
		std::unordered_map<std::string, std::unordered_map<std::string, uint64_t>> allSeqCounts,
		uint32_t mean, uint32_t std, uint32_t minLen, bool verbose) {
	std::unordered_map<std::string, uint64_t> ret;
	sim::ReadLenNormalDistribution uiNormDist(mean, std);
	njh::randomGenerator gen;
	if(verbose){
		std::cout << "Shearing library" << std::endl;
	}
	for (const auto & seqCount : allSeqCounts) {
		if(verbose){
			std::cout << "\tShearing " << seqCount.first << std::endl;
		}
//		for (const auto & seq : seqCount.second) {
//			std::cout << seq.first << std::endl;
//		}
		/*
		 *
		for (const auto & seq : seqCount.second) {
			std::cout << seq.first << std::endl;
			auto compSeq = seqUtil::reverseComplement(seq.first, "DNA");
			for (uint32_t i = 0; i < seq.second; ++i) {
				std::cout << i << std::endl;
				std::cout << seq.first.length() << std::endl;
				auto readPositionsLens = sim::genFragPosSizes(seq.first.length(),
						uiNormDist, gen);
		 */

		njh::stopWatch watch;
		for (const auto & seq : seqCount.second) {
			//std::cout << "__LINE__: " << __LINE__ << std::endl;
//			std::cout << "seq: " << seq.first << std::endl;
			auto compSeq = seqUtil::reverseComplement(seq.first, "DNA");
			//std::cout << "__LINE__: " << __LINE__ << std::endl;
			for (uint32_t i = 0; i < seq.second; ++i) {
//				std::cout << "i: " <<  i << std::endl;
//				std::cout << "seq.first.length(): " << seq.first.length() << std::endl;
				//std::cout << "__LINE__: " << __LINE__ << std::endl;
				auto readPositionsLens = sim::genFragPosSizes(seq.first.length(),
						uiNormDist, gen);
				//std::cout << "__LINE__: " << __LINE__ << std::endl;
				if(readPositionsLens.empty()){
					continue;
				}
				//std::cout << "__LINE__: " << __LINE__ << std::endl;
				auto fragments = sim::genFragments(seq.first, readPositionsLens);
				//std::cout << "__LINE__: " << __LINE__ << std::endl;
				auto compFragments = sim::genFragments(compSeq, readPositionsLens);
				//std::cout << "__LINE__: " << __LINE__ << std::endl;
				for (const auto & frag : fragments) {
					if(frag.length() >=minLen){
						++ret[frag];
					}
				}
				//std::cout << "__LINE__: " << __LINE__ << std::endl;
				for (const auto & frag : compFragments) {
					if(frag.length() >=minLen){
						++ret[frag];
					}
				}
				//std::cout << "__LINE__: " << __LINE__ << std::endl;
			}
		}
		if(verbose){
			std::cout << "\t" << watch.totalTimeFormatted(6) << std::endl;
		}
	}
	return ret;
}

void simShotgunLibFast(std::vector<std::shared_ptr<seqInfo> > & reads,
		const std::string & workingDir,
		const std::string & libName, uint64_t intErrorRate,
		uint64_t startingTemplate, uint64_t finalReadAmount, uint32_t pcrRounds,
		uint32_t initialPcrRounds, uint32_t numThreads, bool verbose, uint32_t mean,
		uint32_t std, uint32_t minLen) {
	if(initialPcrRounds >= pcrRounds){
		std::stringstream ss;
		ss << "Error in " << __PRETTY_FUNCTION__ << std::endl;
		ss << "initialPcrRounds should be less than pcrRounds" << std::endl;
		ss << "pcrRounds:" << pcrRounds << std::endl;
		ss << "initialPcrRounds:" << initialPcrRounds << std::endl;
		throw std::runtime_error{ss.str()};
	}
	//auto finalPerfectAmount = static_cast<uint64_t>(startingTemplate * std::pow(2, pcrRounds));
	auto libDirName = njh::files::makeDir(workingDir, njh::files::MkdirPar(libName));
	std::ofstream libOutFile;
	openTextFile(libOutFile, OutOptions(bfs::path(libDirName.string() + "reads.fasta")));
	std::mutex seqFileLock;
	//check there is enough template and final read amount to make reads for the desired fractions
	std::unordered_map<std::string, uint64_t> templateAmountCounts;
	std::unordered_map<std::string, uint64_t> finalReadAmountCounts;
	uint64_t maxTemplateReads = 0;
	uint64_t readTemplateSum = 0;
	std::string maxTemplateRead = "";
	uint64_t maxFinalReads = 0;
	uint64_t readFinalSum = 0;
	std::string maxFinalRead = "";
	for(const auto & read : reads){
		uint64_t currentTemplateAmt = std::round(read->frac_ * startingTemplate);
		uint64_t currentFinalAmt = std::round(read->frac_ * finalReadAmount);
		if(currentTemplateAmt == 0){
			std::stringstream ss;
			ss << "Error in: " << __PRETTY_FUNCTION__
					<< ", not enough starting template reads, " << startingTemplate
					<< ", requested in order to simulate the desired read fraction:"
					<< read->frac_ << std::endl;
			throw std::runtime_error{njh::bashCT::boldRed(ss.str())};
		}
		if (currentFinalAmt == 0) {
			std::stringstream ss;
			ss << "Error in: " << __PRETTY_FUNCTION__
					<< ", not enough final reads, " << finalReadAmount
					<< ", requested in order to simulate the desired read fraction:"
					<< read->frac_ << std::endl;
			throw std::runtime_error{njh::bashCT::boldRed(ss.str())};
		}
		templateAmountCounts[read->name_] = currentTemplateAmt;
		finalReadAmountCounts[read->name_] = currentFinalAmt;
		readTemplateSum += currentTemplateAmt;
		readFinalSum += currentFinalAmt;
		if(currentTemplateAmt > maxTemplateReads){
			maxTemplateReads = currentTemplateAmt;
			maxTemplateRead = read->name_;
		}
		if(currentFinalAmt > maxFinalReads){
			maxFinalReads = currentFinalAmt;
			maxFinalRead = read->name_;
		}
	}
	if(readTemplateSum < startingTemplate){
		templateAmountCounts[maxTemplateRead] += startingTemplate - readTemplateSum;
	}

	if(readFinalSum < finalReadAmount){
		finalReadAmountCounts[maxFinalRead] += finalReadAmount - readFinalSum;
	}

	if(verbose){
		std::cout << "Simulating library " << libName << std::endl;
		std::cout << "Starting total template: " << startingTemplate << std::endl;
		for(const auto & tempAmount : templateAmountCounts){
			std::cout << tempAmount.first << ": " << tempAmount.second << std::endl;
		}
	}
	std::unordered_map<std::string,std::unordered_map<std::string, uint64_t>> allSeqCounts;
	std::unordered_map<std::string,std::pair<uint64_t,uint64_t>> templateNonMutated;
	for (const auto & read : reads) {
		std::unordered_map<std::string, uint64_t> seqCounts;
		std::mutex seqMapLock;
		if(verbose){
			std::cout << "Simulating ref read " << read->name_ << " for library " << libName << std::endl;
			std::cout << "\tTemplate Amount:          " << templateAmountCounts[read->name_] << std::endl;
			std::cout << "\tFinal Read Sample amount: " << finalReadAmountCounts[read->name_] << std::endl;
		}
		auto finalAmount = runPcr(intErrorRate, numThreads, initialPcrRounds, read->seq_,
				templateAmountCounts[read->name_], read->name_, seqCounts, seqMapLock,
				verbose);
		uint64_t finalPerfectAmount = templateAmountCounts[read->name_] * std::pow(2, initialPcrRounds);
		templateNonMutated[read->name_] = {finalAmount, finalPerfectAmount};
		allSeqCounts[read->name_] = seqCounts;
	}
	//shear library
	auto shotgunReads = shearLibrary(allSeqCounts, mean, std, minLen, verbose);
	//sample and finish PCR
	auto sampleNumber = sampleReadsWithoutReplacementFinishPCR(shotgunReads,
			finalReadAmount, libOutFile, seqFileLock, pcrRounds - initialPcrRounds,
			intErrorRate, numThreads, verbose);
	if(verbose){
		{
			std::cout << "PCR amounts: " << std::endl;
			uint64_t nonMutated = 0;
			uint64_t mutated = 0;
			for (const auto & read : reads) {
				std::cout << read->name_ << std::endl;
				nonMutated += templateNonMutated[read->name_].first;
				mutated += templateNonMutated[read->name_].second
						- templateNonMutated[read->name_].first;
				std::cout << "\t"
						<< getPercentageString(templateNonMutated[read->name_].first,
								templateNonMutated[read->name_].second) << std::endl;
			}
			std::cout << "Total: " << nonMutated + mutated << std::endl;
			std::cout << "\t" << getPercentageString(nonMutated, nonMutated + mutated)
					<< std::endl;
			std::cout << "\t" << getPercentageString(mutated, nonMutated + mutated)
					<< std::endl;
		}
		std::cout << "Sampling Amounts:" << std::endl;
		std::cout << "Total\t       : " << sampleNumber.first + sampleNumber.second << std::endl;
		std::cout << "\tNon-Mutated : " << getPercentageString(sampleNumber.first, sampleNumber.first + sampleNumber.second) << std::endl;
		std::cout << "\tMutated     : " << getPercentageString(sampleNumber.second, sampleNumber.first + sampleNumber.second) << std::endl;
	}
	libOutFile.close();
}

std::unordered_map<uint32_t, std::unordered_map<std::string, double>> processLibraryAbundances(
		const std::string & abundanceFile, VecStr refNames,
		uint32_t maxLibraryAmount) {
	table abundTab(abundanceFile);
	std::unordered_map<uint32_t, std::unordered_map<std::string, double>> libraryAbundances;
	auto refNamesFromTable = abundTab.getColumn(0);
	for(const auto & refName : refNamesFromTable){
		if(!njh::in(refName, refNames)){
			std::stringstream ss;
			ss << "Error in abundance table, requesting read that isn't in supplied sequences" << std::endl;
			ss << "Requested: " << refName << std::endl;
			ss << "Options are: " << vectorToString(refNames) << std::endl;
			throw std::runtime_error{njh::bashCT::boldRed(ss.str())};
		}
	}
	if(maxLibraryAmount < abundTab.columnNames_.size() - 1){
		std::stringstream ss;
		ss << "Error in abundance table, requesting more barcode samples than there are barcodes supplied" << std::endl;
		ss << "Requested: " << abundTab.columnNames_.size() - 1 << std::endl;
		ss << "Barcode Number: " << maxLibraryAmount << std::endl;
		throw std::runtime_error{njh::bashCT::boldRed(ss.str())};
	}
	for(const auto pos : iter::range<size_t>(1, abundTab.columnNames_.size())){
		auto currentLibAbund = abundTab.getColumn(pos);
		if (!isVecOfDoubleStr(currentLibAbund)) {
			std::stringstream ss;
			ss << "Error in abundance table, column " << pos
					<< " has non numeric values" << std::endl;
			throw std::runtime_error { njh::bashCT::boldRed(ss.str()) };
		}
		double sum = 0;
		std::unordered_map<std::string, double> abundances;
		for(const auto abundPos : iter::range(currentLibAbund.size())){
			double abund =  njh::lexical_cast<double>(currentLibAbund[abundPos]);
			if(abund > 0){
				sum += abund;
				abundances[refNamesFromTable[abundPos]] = abund;
			}
		}
		for(auto & abund : abundances){
			abund.second /= sum;
		}
		libraryAbundances[pos] = abundances;
	}
	return libraryAbundances;
}

}  // namespace sim
}  // namespace njhseq

