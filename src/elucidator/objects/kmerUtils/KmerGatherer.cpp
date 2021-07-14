/*
 * KmerGatherer.cpp
 *
 *  Created on: Jul 13, 2021
 *      Author: nick
 */

#include "KmerGatherer.hpp"
#include "elucidator/objects/kmerUtils/SimpleKmerHash.hpp"
#include "elucidator/utils/KmerUtils.hpp"
#include "elucidator/objects/dataContainers.h"
#include "elucidator/simulation.h"
#include "elucidator/objects/seqObjects/seqKmers.h"


namespace njhseq {

KmerGatherer::KmerGathererPars::KmerGathererPars(uint32_t kmerLength,
		bool noRevComp, uint32_t numThreads, std::set<char> allowableCharacters) :
		kmerLength_(kmerLength), noRevComp_(noRevComp), numThreads_(numThreads), allowableCharacters_(
				allowableCharacters) {

}
KmerGatherer::KmerGathererPars::KmerGathererPars() {

}

void KmerGatherer::KmerGathererPars::setOptions(seqSetUp &setUp) {
	setUp.setOption(kmerLength_, "--kmerLength", "kmer Length");
	setUp.setOption(numThreads_, "--numThreads", "number of threads");
	setUp.setOption(allowableCharacters_, "--allowableCharacters",
			"Only count kmers with these allowable Characters");
	setUp.setOption(noRevComp_, "--noRevComp", "noRevComp");
	setUp.setOption(entropyFilter_, "--entropyFilter", "entropy Filter cut off, exclusive, will only keep kmers abovet this entropy level");
}

KmerGatherer::KmerGatherer(const KmerGathererPars &pars) :
		pars_(pars) {

}

std::unordered_map<std::string, uint32_t> KmerGatherer::countGenomeKmers(
		const bfs::path &genomeFnp) const {
	std::unordered_map < std::string, uint32_t > ret;
	SeqInput reader(SeqIOOptions::genFastaIn(genomeFnp));
	reader.openIn();
	std::mutex genomeCountsMut;
	std::function < void() > countKmers =
			[&reader, this, &genomeCountsMut, &ret]() {
				seqInfo seq;

				while (reader.readNextReadLock(seq)) {
					std::unordered_map < std::string, uint32_t > genomeCountsCurrent;
					if (len(seq) > pars_.kmerLength_) {
						for (uint32_t pos = 0; pos < len(seq) - pars_.kmerLength_ + 1;
								++pos) {
							++genomeCountsCurrent[seq.seq_.substr(pos, pars_.kmerLength_)];
						}
						if (!pars_.noRevComp_) {
							seq.reverseComplementRead(false, true);
							for (uint32_t pos = 0; pos < len(seq) - pars_.kmerLength_ + 1;
									++pos) {
								++genomeCountsCurrent[seq.seq_.substr(pos, pars_.kmerLength_)];
							}
						}
					}
					{
						std::lock_guard < std::mutex > lock(genomeCountsMut);
						for (const auto &count : genomeCountsCurrent) {
							ret[count.first] += count.second;
						}
					}
				}
			};
	njh::concurrent::runVoidFunctionThreaded(countKmers, pars_.numThreads_);

	return ret;
}

std::unordered_set<std::string> KmerGatherer::getUniqueKmers(
		const bfs::path &genomeFnp) const {
	std::unordered_set < std::string > ret;
	SeqInput reader(SeqIOOptions::genFastaIn(genomeFnp));
	reader.openIn();
	std::mutex genomeKmersMut;
	std::function < void() > gatherKmers =
			[&reader, this, &genomeKmersMut, &ret]() {
				seqInfo seq;

				while (reader.readNextReadLock(seq)) {
					std::unordered_set < std::string > genomeKmersCurrent;
					if (len(seq) > pars_.kmerLength_) {
						for (uint32_t pos = 0; pos < len(seq) - pars_.kmerLength_ + 1;
								++pos) {
							genomeKmersCurrent.emplace(
									seq.seq_.substr(pos, pars_.kmerLength_));
						}
						if (!pars_.noRevComp_) {
							seq.reverseComplementRead(false, true);
							for (uint32_t pos = 0; pos < len(seq) - pars_.kmerLength_ + 1;
									++pos) {
								genomeKmersCurrent.emplace(
										seq.seq_.substr(pos, pars_.kmerLength_));
							}
						}
					}
					{
						std::lock_guard < std::mutex > lock(genomeKmersMut);
						ret.insert(genomeKmersCurrent.begin(), genomeKmersCurrent.end());
					}
				}
			};
	njh::concurrent::runVoidFunctionThreaded(gatherKmers, pars_.numThreads_);
	return ret;
}
std::set<std::string> KmerGatherer::getUniqueKmersSet(
		const bfs::path &genomeFnp) const {
	std::set < std::string > ret;
	SeqInput reader(SeqIOOptions::genFastaIn(genomeFnp));
	reader.openIn();
	std::mutex genomeKmersMut;
	std::function < void() > gatherKmers =
			[&reader, this, &genomeKmersMut, &ret]() {
				seqInfo seq;
				while (reader.readNextReadLock(seq)) {
					std::set < std::string > genomeKmersCurrent;
					if (len(seq) > pars_.kmerLength_) {
						for (uint32_t pos = 0; pos < len(seq) - pars_.kmerLength_ + 1;
								++pos) {
							genomeKmersCurrent.emplace(
									seq.seq_.substr(pos, pars_.kmerLength_));
						}
						if (!pars_.noRevComp_) {
							seq.reverseComplementRead(false, true);
							for (uint32_t pos = 0; pos < len(seq) - pars_.kmerLength_ + 1;
									++pos) {
								genomeKmersCurrent.emplace(
										seq.seq_.substr(pos, pars_.kmerLength_));
							}
						}
					}
					{
						std::lock_guard < std::mutex > lock(genomeKmersMut);
						ret.insert(genomeKmersCurrent.begin(), genomeKmersCurrent.end());
					}
				}
			};
	njh::concurrent::runVoidFunctionThreaded(gatherKmers, pars_.numThreads_);
	return ret;
}

std::unordered_map<std::string, std::set<std::string>> KmerGatherer::getUniqueKmersSet(
		const std::vector<bfs::path> &twobitFnps) const {
	std::vector < TwobitFnpSeqNamePair > pairs;

	for (const auto &twoBit : twobitFnps) {
		TwoBit::TwoBitFile tReader(twoBit);
		auto seqNames = tReader.sequenceNames();
		for (const auto &seqName : seqNames) {
			pairs.emplace_back(TwobitFnpSeqNamePair(twoBit, seqName));
		}
	}
	njh::concurrent::LockableQueue < TwobitFnpSeqNamePair > pairsQueue(pairs);

	std::unordered_map<std::string, std::set<std::string>> ret;
	std::mutex mut;

	std::function < void() > gatherKmers = [&pairsQueue, this, &ret, &mut]() {

		TwobitFnpSeqNamePair pair;
		std::unordered_map<std::string, std::set<std::string>> current;
		while (pairsQueue.getVal(pair)) {
			std::set < std::string > genomeKmersCurrent;
			TwoBit::TwoBitFile tReader(pair.twoBit_);
			std::string buffer;
			tReader[pair.seqName_]->getSequence(buffer);
			for (uint32_t pos = 0; pos < len(buffer) - pars_.kmerLength_ + 1; ++pos) {
				genomeKmersCurrent.emplace(buffer.substr(pos, pars_.kmerLength_));
			}
			if (!pars_.noRevComp_) {
				buffer = seqUtil::reverseComplement(buffer, "DNA");
				for (uint32_t pos = 0; pos < len(buffer) - pars_.kmerLength_ + 1;
						++pos) {
					genomeKmersCurrent.emplace(buffer.substr(pos, pars_.kmerLength_));
				}
			}
			current[pair.twoBit_.string()].insert(genomeKmersCurrent.begin(),
					genomeKmersCurrent.end());
		}
		{
			std::lock_guard < std::mutex > lock(mut);
			for (const auto &kmerSet : current) {
				ret[kmerSet.first].insert(kmerSet.second.begin(), kmerSet.second.end());
			}
		}
	};
	njh::concurrent::runVoidFunctionThreaded(gatherKmers, pars_.numThreads_);
	return ret;
}



std::unordered_map<std::string, std::set<uint64_t>> KmerGatherer::getUniqueKmersSetHashWithFilters(
		const std::vector<bfs::path> &twobitFnps) const {
	std::vector < TwobitFnpSeqNamePair > pairs;
	if (pars_.kmerLength_ > 19) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error "
				<< "cannot do kmer lengths greater than 19" << "\n";
		throw std::runtime_error { ss.str() };
	}
	for (const auto &twoBit : twobitFnps) {
		TwoBit::TwoBitFile tReader(twoBit);
		auto seqNames = tReader.sequenceNames();
		for (const auto &seqName : seqNames) {
			pairs.emplace_back(TwobitFnpSeqNamePair(twoBit, seqName));
		}
	}
	njh::concurrent::LockableQueue < TwobitFnpSeqNamePair > pairsQueue(pairs);

	std::unordered_map<std::string, std::set<uint64_t>> ret;
	std::mutex mut;
	std::function<bool(const std::string&)> seqCheck = [this](const std::string & k){
		return std::all_of(k.begin(), k.end(), [this](char base){return njh::in(base, pars_.allowableCharacters_);});
	};

	std::function < void() > gatherKmers = [&pairsQueue, this, &ret, &mut,&seqCheck]() {
		SimpleKmerHash hasher;
		TwobitFnpSeqNamePair pair;
		std::unordered_map<std::string, std::set<uint64_t>> current;
		while (pairsQueue.getVal(pair)) {
			std::set < uint64_t > genomeKmersCurrent;
			TwoBit::TwoBitFile tReader(pair.twoBit_);
			std::string buffer;
			tReader[pair.seqName_]->getSequence(buffer);
			for (uint32_t pos = 0; pos < len(buffer) - pars_.kmerLength_ + 1; ++pos) {
				auto k = buffer.substr(pos, pars_.kmerLength_);
				if(seqCheck(k)){
					DNABaseCounter counter(pars_.allowableCharacters_);
					counter.increase(k);
					if(counter.computeEntrophy() > pars_.entropyFilter_){
						genomeKmersCurrent.emplace(hasher.hash(k));
					}
				}
			}
			if (!pars_.noRevComp_) {
				//buffer = seqUtil::reverseComplement(buffer, "DNA");
				for (uint32_t pos = 0; pos < len(buffer) - pars_.kmerLength_ + 1;
						++pos) {
					auto k = buffer.substr(pos, pars_.kmerLength_);
					if(seqCheck(k)){
						DNABaseCounter counter(pars_.allowableCharacters_);
						counter.increase(k);
						if(counter.computeEntrophy() > pars_.entropyFilter_){
							genomeKmersCurrent.emplace(hasher.revCompHash(k));
						}
					}
				}
			}
			current[pair.twoBit_.string()].insert(genomeKmersCurrent.begin(),genomeKmersCurrent.end());
		}
		{
			std::lock_guard < std::mutex > lock(mut);
			for (const auto &kmerSet : current) {
				ret[kmerSet.first].insert(kmerSet.second.begin(), kmerSet.second.end());
			}
		}
	};
	njh::concurrent::runVoidFunctionThreaded(gatherKmers, pars_.numThreads_);
	return ret;
}

std::unordered_map<std::string, std::set<uint64_t>> KmerGatherer::getUniqueKmersSetHash(
		const std::vector<bfs::path> &twobitFnps) const {
	std::vector < TwobitFnpSeqNamePair > pairs;
	if (pars_.kmerLength_ > 19) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error "
				<< "cannot do kmer lengths greater than 19" << "\n";
		throw std::runtime_error { ss.str() };
	}
	for (const auto &twoBit : twobitFnps) {
		TwoBit::TwoBitFile tReader(twoBit);
		auto seqNames = tReader.sequenceNames();
		for (const auto &seqName : seqNames) {
			pairs.emplace_back(TwobitFnpSeqNamePair(twoBit, seqName));
		}
	}
	njh::concurrent::LockableQueue < TwobitFnpSeqNamePair > pairsQueue(pairs);

	std::unordered_map<std::string, std::set<uint64_t>> ret;
	std::mutex mut;

	std::function < void() > gatherKmers = [&pairsQueue, this, &ret, &mut]() {
		SimpleKmerHash hasher;
		TwobitFnpSeqNamePair pair;
		std::unordered_map<std::string, std::set<uint64_t>> current;
		while (pairsQueue.getVal(pair)) {
			std::set < uint64_t > genomeKmersCurrent;
			TwoBit::TwoBitFile tReader(pair.twoBit_);
			std::string buffer;
			tReader[pair.seqName_]->getSequence(buffer);
			for (uint32_t pos = 0; pos < len(buffer) - pars_.kmerLength_ + 1; ++pos) {
				genomeKmersCurrent.emplace(
						hasher.hash(buffer.substr(pos, pars_.kmerLength_)));
			}
			if (!pars_.noRevComp_) {
				buffer = seqUtil::reverseComplement(buffer, "DNA");
				for (uint32_t pos = 0; pos < len(buffer) - pars_.kmerLength_ + 1;
						++pos) {
					genomeKmersCurrent.emplace(
							hasher.hash(buffer.substr(pos, pars_.kmerLength_)));
				}
			}
			current[pair.twoBit_.string()].insert(genomeKmersCurrent.begin(),
					genomeKmersCurrent.end());
		}
		{
			std::lock_guard < std::mutex > lock(mut);
			for (const auto &kmerSet : current) {
				ret[kmerSet.first].insert(kmerSet.second.begin(), kmerSet.second.end());
			}
		}
	};
	njh::concurrent::runVoidFunctionThreaded(gatherKmers, pars_.numThreads_);
	return ret;
}

}  // namespace njhseq

