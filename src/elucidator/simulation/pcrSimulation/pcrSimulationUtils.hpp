#pragma once
/*
 * pcrSimulationUtils.hpp
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
#include "elucidator/common.h"
namespace njhseq {
namespace sim {

class ReadLenNormalDistribution {
	uint32_t minAllowed_;
	std::mt19937 rgen_;
	std::normal_distribution<> ndis_;

public:

	ReadLenNormalDistribution(double mean, double std, int32_t minLen = 1);
	ReadLenNormalDistribution(double mean, double std, int32_t minLen, uint64_t seed);
	uint32_t operator ()();
};


VecStr genFragments(const std::string & seq,
		const std::vector<std::pair<size_t, uint32_t>> & fragPosLens);

std::vector<std::pair<size_t, uint32_t>> genFragPosSizes(size_t seqLen,
		ReadLenNormalDistribution & uiNormDist, njh::randomGenerator & gen);


std::unordered_map<uint32_t, uint32_t> genMutCounts(std::mt19937_64 & gen,
		const uint64_t seqNumber, const uint64_t seqSize,
		const uint64_t intErrorRate);

std::vector<uint64_t> spreadReadNumAcrossThreads(const uint64_t readTotal,
		const uint32_t numThreads);

std::pair<std::string, std::string> mutateSeq(std::string name, std::string seq,
		uint32_t mutNum, const std::vector<uint64_t>& readPositions,
		njh::randomGenerator & gen,
		std::unordered_map<char, njh::randObjectGen<char, uint32_t>>& charGen) ;



uint64_t runPcr(uint64_t intErrorRate, uint32_t numThreads, uint32_t pcrRounds,
		std::string seq, uint64_t startingTemplate, std::string name,
		std::unordered_map<std::string, uint64_t> & seqMap,
		std::mutex & seqMapLock, bool printTime);





std::pair<uint64_t,uint64_t> sampleReadsWithReplacement(const std::string & name, const std::string & seq,
		std::unordered_map<std::string, uint64_t> & seqCounts,
		uint64_t finalPerfectAmount, uint64_t finalReadAmount,
		std::ostream & sequenceOutFile, std::mutex & seqFileLock,
		uint32_t numThreads);

std::pair<uint64_t,uint64_t> sampleReadsWithoutReplacement(const std::string & name,
		const std::string & seq,
		std::unordered_map<std::string, uint64_t> & seqCounts,
		uint64_t finalReadAmount,
		std::ostream & sequenceOutFile, std::mutex & seqFileLock,
		bool verbose);

std::unordered_map<std::string, std::pair<uint64_t, uint64_t>> sampleReadsWithoutReplacement(
		const std::unordered_map<std::string, std::string> & seqs,
		std::unordered_map<std::string, std::unordered_map < std::string, uint64_t>> & multipleSeqCounts,
		uint64_t finalReadAmount, std::ostream & sequenceOutFile,
		std::mutex & seqFileLock, bool verbose);

std::string runPcrSingleTemplate(uint64_t intErrorRate,uint32_t roundsOfPcr,
		uint64_t randomNumberSelector,uint64_t randomNumberSelectorMax, std::string seq );

std::unordered_map<std::string, std::pair<uint64_t, uint64_t>> sampleReadsWithoutReplacementFinishPCR(
		const std::unordered_map<std::string, std::string> & seqs,
		std::unordered_map<std::string, std::unordered_map<std::string, uint64_t>> & multipleSeqCounts,
		uint64_t finalReadAmount, std::ostream & sequenceOutFile,
		std::mutex & seqFileLock, uint32_t numberOfPCRRoundsLeft,
		uint64_t intErrorRate, uint32_t numThreads, bool verbose);

std::pair<uint64_t, uint64_t> sampleReadsWithoutReplacementFinishPCR(
		std::unordered_map<std::string, uint64_t> & seqCounts,
		uint64_t finalReadAmount, std::ostream & sequenceOutFile,
		std::mutex & seqFileLock, uint32_t numberOfPCRRoundsLeft,
		uint64_t intErrorRate, uint32_t numThreads, bool verbose) ;

void simLib(std::vector<std::shared_ptr<seqInfo> > & reads, const std::string & barcode,
		const std::string & workingDir, const std::string & libName,
		uint64_t intErrorRate, uint64_t startingTemplate, uint64_t finalReadAmount,
		uint32_t pcrRounds, uint32_t numThreads, bool verbose);

void simLibFast(std::vector<std::shared_ptr<seqInfo> > & reads,
		const std::string & barcode, const std::string & workingDir,
		const std::string & libName, uint64_t intErrorRate,
		uint64_t startingTemplate, uint64_t finalReadAmount, uint32_t pcrRounds,
		uint32_t initialPcrRounds, uint32_t numThreads, bool verbose);

std::unordered_map<std::string, uint64_t> shearLibraryThreaded(
		std::unordered_map<std::string, std::unordered_map<std::string, uint64_t>> allSeqCounts,
		uint32_t mean, uint32_t std, uint32_t minLen, bool verbose);



std::unordered_map<std::string, uint64_t> shearLibrary(
		std::unordered_map<std::string, std::unordered_map<std::string, uint64_t>> allSeqCounts,
		uint32_t mean, uint32_t std, uint32_t minLen, bool verbose);

void simShotgunLibFast(std::vector<std::shared_ptr<seqInfo> > & reads,
		const std::string & workingDir,
		const std::string & libName, uint64_t intErrorRate,
		uint64_t startingTemplate, uint64_t finalReadAmount, uint32_t pcrRounds,
		uint32_t initialPcrRounds, uint32_t numThreads, bool verbose, uint32_t mean,
		uint32_t std, uint32_t minLen);

std::unordered_map<uint32_t, std::unordered_map<std::string, double>> processLibraryAbundances(
		const std::string & abundanceFile, VecStr refNames,
		uint32_t maxLibraryAmount) ;



}  // namespace sim
}  // namespace njhseq

