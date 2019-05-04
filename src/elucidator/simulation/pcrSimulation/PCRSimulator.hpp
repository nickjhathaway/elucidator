#pragma once

/*
 * PCRSimulator.hpp
 *
 *  Created on: May 2, 2019
 *      Author: nicholashathaway
 */


#include "elucidator/common.h"


namespace njhseq {

class PCRSimulator {
public:
	/**@brief Construct with the error rate (number between 0-1) converted into the uint64_t equalivalent
	 *
	 * @param intErrorRate error rate (rate times the uint64_t max, this will be compared to a random number generated between 0 and uin64_t max)
	 */
	PCRSimulator(uint64_t intErrorRate);

	/**@brief construct with PCR error rate, should be a number between 0 and 1
	 *
	 * @param errorRate the PCR error rate
	 */
	PCRSimulator(double errorRate);

	uint64_t intErrorRate_;/**< the error rate of PCR converted into uint64_t number, takes the rate and times std::numeric_limits<uint64_t>::max() */
	bool verbose_{false};
	double pcrEfficiency_{0.95}; /**< chance a product gets amplfied into the next round*/

	/**@brief Generate the number of mutated sequences with how many bases mutated
	 *
	 * @param gen the random number generator
	 * @param seqNumber the number of times this sequence should be attempted to be mutated
	 * @param seqSize the number of bases in this sequence
	 * @return an unordred_map, k: number of mutated bases, v: number of times this many mutates was observed
	 */
	std::unordered_map<uint32_t, uint32_t> genMutCounts(
			std::mt19937_64 & gen,
			const uint64_t seqNumber,
			const uint64_t seqSize) const;


	/**@brief Mutate a sequence a certain number of times
	 *
	 * @param seq the seq to mutate
	 * @param mutNum the number of bases to mutate
	 * @param readPositions positions that can be mutated
	 * @param gen a random number generator to pick positions to mutate
	 * @param charGen the mutation generator
	 * @return the mutated sequence
	 */
	static std::string mutateSeq(std::string seq,
			uint32_t mutNum,
			const std::vector<uint64_t>& readPositions,
			njh::randomGenerator & gen,
			std::unordered_map<char, njh::randObjectGen<char, uint32_t>>& charGen);

	/**@brief spread out the number of reads to simulate across the threads
	 *
	 * @param readTotal the total number of reads to simulate
	 * @param numThreads the number of threads to spread across
	 * @return a vector of size numThreads with readTotal evenly spread across it
	 */
	static std::vector<uint64_t> spreadReadNumAcrossThreads(const uint64_t readTotal,
			const uint32_t numThreads);

	struct PCRProduct {
		PCRProduct(const std::string & pcrSeq, const uint32_t pcrRoundsLeft, const uint32_t templateAmount);
		PCRProduct(const std::string & pcrSeq, const uint32_t pcrRoundsLeft);

		std::string pcrSeq_;
		uint32_t pcrRoundsLeft_;
		uint32_t templateAmount_{1};

	};


	void runPcr(uint32_t numThreads,
			uint32_t pcrRounds,
			std::string currentSeq,
			uint64_t currentStartingTemplate,
			std::unordered_map<std::string, uint64_t> & currentSeqMap,
			std::mutex & currentSeqMapLock) const ;

	void simLibFast(const std::vector<seqInfo> & reads,
			const OutOptions & outputFile,
			uint64_t startingTemplate,
			uint64_t finalReadAmount,
			uint32_t pcrRounds,
			uint32_t initialPcrRounds,
			uint32_t numThreads);

	std::unordered_map<std::string, std::pair<uint64_t, uint64_t>> sampleReadsWithoutReplacementFinishPCR(
			const std::unordered_map<std::string, std::string> & seqs,
			std::unordered_map<std::string, std::unordered_map<std::string, uint64_t>> & multipleSeqCounts,
			uint64_t finalReadAmount, std::ostream & sequenceOutFile,
			std::mutex & seqFileLock, uint32_t numberOfPCRRoundsLeft,
			uint32_t numThreads);

	std::string runPcrSampleSingleTemplate(
			uint32_t roundsOfPcr,
			uint64_t randomNumberSelector,
			uint64_t randomNumberSelectorMax,
			std::string seq );
};

}  // namespace njhseq

