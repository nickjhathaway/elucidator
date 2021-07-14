#pragma once

/*
 * KmerGatherer.hpp
 *
 *  Created on: Jul 13, 2021
 *      Author: nick
 */


#include <njhseq/common.h>
#include <njhseq/programUtils/seqSetUp.hpp>


namespace njhseq {

class KmerGatherer{
public:
	struct KmerGathererPars {
		KmerGathererPars(uint32_t kmerLength,
				bool noRevComp, uint32_t numThreads,
				std::vector<char> allowableCharacters);
		KmerGathererPars();
		uint32_t kmerLength_{19};
		bool noRevComp_{false};
		uint32_t numThreads_{1};
		std::vector<char> allowableCharacters_{'A', 'C', 'G', 'T'};

		void setOptions(seqSetUp & setUp);
	};

	KmerGatherer(const KmerGathererPars & pars);
	KmerGathererPars pars_;

	std::unordered_map<std::string, uint32_t> countGenomeKmers(const bfs::path & genomeFnp) const;

	std::unordered_set<std::string> getUniqueKmers(const bfs::path & genomeFnp) const;
	std::set<std::string> getUniqueKmersSet(const bfs::path & genomeFnp) const;

	struct TwobitFnpSeqNamePair{
		TwobitFnpSeqNamePair(const bfs::path twoBit, const std::string & seqName):twoBit_(twoBit), seqName_(seqName){

		}
		TwobitFnpSeqNamePair(){

		}
		bfs::path twoBit_;
		std::string seqName_;
	};
	std::unordered_map<std::string, std::set<std::string>> getUniqueKmersSet(const std::vector<bfs::path> & twobitFnps) const;

	std::unordered_map<std::string, std::set<uint64_t>> getUniqueKmersSetHash(const std::vector<bfs::path> & twobitFnps) const;
};


}  // namespace njhseq




