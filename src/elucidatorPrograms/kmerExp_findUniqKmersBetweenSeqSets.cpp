/*
 * kmerExp_findUniqKmersBetweenSeqSets.cpp
 *
 *  Created on: Jul 10, 2021
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

#include "kmerExp.hpp"
#include "elucidator/utils/KmerUtils.hpp"
#include "elucidator/objects/dataContainers.h"
#include "elucidator/simulation.h"
#include "elucidator/objects/seqObjects/seqKmers.h"

#include "elucidator/objects/MiscUtility/GenomeSeqSearch.hpp"


namespace njhseq {

class SimpleKmerHash{

public:
	SimpleKmerHash(){
		hasher_ = std::vector<char>(255, '5');
		hasher_['A'] = '1';
		hasher_['C'] = '2';
		hasher_['G'] = '3';
		hasher_['T'] = '4';
		hasher_['N'] = '5';
		reverseHasher_ = std::vector(255, 'N');
		reverseHasher_['1'] = 'A';
		reverseHasher_['2'] = 'C';
		reverseHasher_['3'] = 'G';
		reverseHasher_['4'] = 'T';
		reverseHasher_['5'] = 'N';
	}
	std::vector<char> hasher_;

	std::vector<char> reverseHasher_;


	uint64_t hash(const std::string & str) const{
		std::string convert;
		for(size_t pos = 0; pos < std::min<size_t>(20, str.size()); ++pos){
			convert.push_back(hasher_[str[pos]]);
		}
		return njh::StrToNumConverter::stoToNum<uint64_t>(convert);
	}

	std::string reverseHash(uint64_t hash) const {

		std::string hashStr = estd::to_string(hash);
		std::string back;
		back.reserve(hashStr.size());
		for(const auto pos : iter::range(hashStr.size())){
			back.push_back(reverseHasher_[hashStr[pos]]);
		}
		return back;
	}
};


class KmerGatherer{
public:
	struct KmerGathererPars {
		KmerGathererPars(uint32_t kmerLength,
				bool noRevComp, uint32_t numThreads,
				std::vector<char> allowableCharacters) :
				kmerLength_(kmerLength), noRevComp_(noRevComp), numThreads_(
						numThreads),
						allowableCharacters_(allowableCharacters){

		}
		KmerGathererPars(){

		}
		uint32_t kmerLength_{21};
		bool noRevComp_{false};
		uint32_t numThreads_{1};
		std::vector<char> allowableCharacters_{'A', 'C', 'G', 'T'};

		void setOptions(seqSetUp & setUp){
			setUp.setOption(kmerLength_, "--kmerLength", "kmer Length");
			setUp.setOption(numThreads_, "--numThreads", "number of threads");
			setUp.setOption(allowableCharacters_, "--allowableCharacters", "Only count kmers with these allowable Characters");
			setUp.setOption(noRevComp_, "--noRevComp", "noRevComp");
		}
	};

	KmerGatherer(const KmerGathererPars & pars):pars_(pars){

	}
	KmerGathererPars pars_;

	std::unordered_map<std::string, uint32_t> countGenomeKmers(const bfs::path & genomeFnp) const {
		std::unordered_map<std::string, uint32_t> ret;
		SeqInput reader(SeqIOOptions::genFastaIn(genomeFnp));
		reader.openIn();
		std::mutex genomeCountsMut;
		std::function<void()> countKmers =[&reader,this,&genomeCountsMut,&ret](){
			seqInfo seq;

			while(reader.readNextReadLock(seq)){
				std::unordered_map<std::string, uint32_t> genomeCountsCurrent;
				if(len(seq) > pars_.kmerLength_){
					for(uint32_t pos = 0; pos < len(seq) - pars_.kmerLength_ + 1; ++pos){
						++genomeCountsCurrent[seq.seq_.substr(pos, pars_.kmerLength_)];
					}
					if(!pars_.noRevComp_){
						seq.reverseComplementRead(false, true);
						for(uint32_t pos = 0; pos < len(seq) - pars_.kmerLength_ + 1; ++pos){
							++genomeCountsCurrent[seq.seq_.substr(pos, pars_.kmerLength_)];
						}
					}
				}
				{
					std::lock_guard<std::mutex> lock(genomeCountsMut);
					for(const auto & count : genomeCountsCurrent){
						ret[count.first] += count.second;
					}
				}
			}
		};
		njh::concurrent::runVoidFunctionThreaded(countKmers, pars_.numThreads_);

		return ret;
	}

	std::unordered_set<std::string> getUniqueKmers(const bfs::path & genomeFnp) const {
		std::unordered_set<std::string> ret;
		SeqInput reader(SeqIOOptions::genFastaIn(genomeFnp));
		reader.openIn();
		std::mutex genomeKmersMut;
		std::function<void()> gatherKmers =[&reader,this,&genomeKmersMut,&ret](){
			seqInfo seq;

			while(reader.readNextReadLock(seq)){
				std::unordered_set<std::string>  genomeKmersCurrent;
				if(len(seq) > pars_.kmerLength_){
					for(uint32_t pos = 0; pos < len(seq) - pars_.kmerLength_ + 1; ++pos){
						genomeKmersCurrent.emplace(seq.seq_.substr(pos, pars_.kmerLength_));
					}
					if(!pars_.noRevComp_){
						seq.reverseComplementRead(false, true);
						for(uint32_t pos = 0; pos < len(seq) - pars_.kmerLength_ + 1; ++pos){
							genomeKmersCurrent.emplace(seq.seq_.substr(pos, pars_.kmerLength_));
						}
					}
				}
				{
					std::lock_guard<std::mutex> lock(genomeKmersMut);
					ret.insert(genomeKmersCurrent.begin(), genomeKmersCurrent.end());
				}
			}
		};
		njh::concurrent::runVoidFunctionThreaded(gatherKmers, pars_.numThreads_);
		return ret;
	}
	std::set<std::string> getUniqueKmersSet(const bfs::path & genomeFnp) const {
		std::set<std::string> ret;
		SeqInput reader(SeqIOOptions::genFastaIn(genomeFnp));
		reader.openIn();
		std::mutex genomeKmersMut;
		std::function<void()> gatherKmers =[&reader,this,&genomeKmersMut,&ret](){
			seqInfo seq;
			while(reader.readNextReadLock(seq)){
				std::set<std::string>  genomeKmersCurrent;
				if(len(seq) > pars_.kmerLength_){
					for(uint32_t pos = 0; pos < len(seq) - pars_.kmerLength_ + 1; ++pos){
						genomeKmersCurrent.emplace(seq.seq_.substr(pos, pars_.kmerLength_));
					}
					if(!pars_.noRevComp_){
						seq.reverseComplementRead(false, true);
						for(uint32_t pos = 0; pos < len(seq) - pars_.kmerLength_ + 1; ++pos){
							genomeKmersCurrent.emplace(seq.seq_.substr(pos, pars_.kmerLength_));
						}
					}
				}
				{
					std::lock_guard<std::mutex> lock(genomeKmersMut);
					ret.insert(genomeKmersCurrent.begin(), genomeKmersCurrent.end());
				}
			}
		};
		njh::concurrent::runVoidFunctionThreaded(gatherKmers, pars_.numThreads_);
		return ret;
	}

	struct TwobitFnpSeqNamePair{
		TwobitFnpSeqNamePair(const bfs::path twoBit, const std::string & seqName):twoBit_(twoBit), seqName_(seqName){

		}
		TwobitFnpSeqNamePair(){

		}
		bfs::path twoBit_;
		std::string seqName_;
	};
	std::unordered_map<std::string, std::set<std::string>> getUniqueKmersSet(const std::vector<bfs::path> & twobitFnps) const {
		std::vector<TwobitFnpSeqNamePair> pairs;

		for(const auto & twoBit : twobitFnps){
			TwoBit::TwoBitFile tReader(twoBit);
			auto seqNames = tReader.sequenceNames();
			for(const auto & seqName : seqNames){
				pairs.emplace_back(TwobitFnpSeqNamePair(twoBit, seqName));
			}
		}
		njh::concurrent::LockableQueue<TwobitFnpSeqNamePair> pairsQueue(pairs);

		std::unordered_map<std::string, std::set<std::string>> ret;
		std::mutex mut;

		std::function<void()> gatherKmers=[&pairsQueue,this,&ret,&mut](){

			TwobitFnpSeqNamePair pair;
			std::unordered_map<std::string, std::set<std::string>> current;
			while(pairsQueue.getVal(pair)){
				std::set<std::string>  genomeKmersCurrent;
				TwoBit::TwoBitFile tReader(pair.twoBit_);
				std::string buffer;
				tReader[pair.seqName_]->getSequence(buffer);
				for(uint32_t pos = 0; pos < len(buffer) - pars_.kmerLength_ + 1; ++pos){
					genomeKmersCurrent.emplace(buffer.substr(pos, pars_.kmerLength_));
				}
				if(!pars_.noRevComp_){
					buffer = seqUtil::reverseComplement(buffer, "DNA");
					for(uint32_t pos = 0; pos < len(buffer) - pars_.kmerLength_ + 1; ++pos){
						genomeKmersCurrent.emplace(buffer.substr(pos, pars_.kmerLength_));
					}
				}
				current[pair.twoBit_.string()].insert(genomeKmersCurrent.begin(), genomeKmersCurrent.end());
			}
			{
				std::lock_guard<std::mutex> lock(mut);
				for(const auto & kmerSet : current){
					ret[kmerSet.first].insert(kmerSet.second.begin(), kmerSet.second.end());
				}
			}
		};
		njh::concurrent::runVoidFunctionThreaded(gatherKmers, pars_.numThreads_);
		return ret;
	}

	std::unordered_map<std::string, std::set<uint64_t>> getUniqueKmersSetHash(const std::vector<bfs::path> & twobitFnps) const {
		std::vector<TwobitFnpSeqNamePair> pairs;
		if(pars_.kmerLength_ > 19){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "cannot do kmer lengths greater than 19" << "\n";
			throw std::runtime_error{ss.str()};
		}
		for(const auto & twoBit : twobitFnps){
			TwoBit::TwoBitFile tReader(twoBit);
			auto seqNames = tReader.sequenceNames();
			for(const auto & seqName : seqNames){
				pairs.emplace_back(TwobitFnpSeqNamePair(twoBit, seqName));
			}
		}
		njh::concurrent::LockableQueue<TwobitFnpSeqNamePair> pairsQueue(pairs);

		std::unordered_map<std::string, std::set<uint64_t>> ret;
		std::mutex mut;

		std::function<void()> gatherKmers=[&pairsQueue,this,&ret,&mut](){
			SimpleKmerHash hasher;
			TwobitFnpSeqNamePair pair;
			std::unordered_map<std::string, std::set<uint64_t>> current;
			while(pairsQueue.getVal(pair)){
				std::set<uint64_t>  genomeKmersCurrent;
				TwoBit::TwoBitFile tReader(pair.twoBit_);
				std::string buffer;
				tReader[pair.seqName_]->getSequence(buffer);
				for(uint32_t pos = 0; pos < len(buffer) - pars_.kmerLength_ + 1; ++pos){
					genomeKmersCurrent.emplace(hasher.hash(buffer.substr(pos, pars_.kmerLength_)));
				}
				if(!pars_.noRevComp_){
					buffer = seqUtil::reverseComplement(buffer, "DNA");
					for(uint32_t pos = 0; pos < len(buffer) - pars_.kmerLength_ + 1; ++pos){
						genomeKmersCurrent.emplace(hasher.hash(buffer.substr(pos, pars_.kmerLength_)));
					}
				}
				current[pair.twoBit_.string()].insert(genomeKmersCurrent.begin(), genomeKmersCurrent.end());
			}
			{
				std::lock_guard<std::mutex> lock(mut);
				for(const auto & kmerSet : current){
					ret[kmerSet.first].insert(kmerSet.second.begin(), kmerSet.second.end());
				}
			}
		};
		njh::concurrent::runVoidFunctionThreaded(gatherKmers, pars_.numThreads_);
		return ret;
	}
};



namespace StrToNumConverter {

/**@brief Function for converting a string to a number, which is just njh::lexical_cast by default and then several specific int conversions are defined for faster converting
 *
 * @param str the string to convert
 * @return the string convert to a number
 */
	template<typename T>
	T stoToNum(const std::string & str){
		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		return njh::lexical_cast<T>(str);
	}

	template<>
	unsigned short stoToNum(const std::string & str){
		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		return estd::stous(str);
	}

	template<>
	unsigned stoToNum(const std::string & str){
		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		return estd::stou(str);
	}

	template<>
	unsigned long stoToNum(const std::string & str){
		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		return std::stoul(str);
	}

	template<>
	unsigned long long stoToNum(const std::string & str){
		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		return std::stoull(str);
	}

	template<>
	short stoToNum(const std::string & str){
		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		return estd::stos(str);
	}

	template<>
	int stoToNum(const std::string & str){
		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		return std::stoi(str);
	}

	template<>
	long int stoToNum(const std::string & str){
		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		return std::stol(str);
	}

	template<>
	long long int stoToNum(const std::string & str){
		return std::stoll(str);
	}

}  // namespace StrToNumConverter




int testingOfWeirdStod(){

	{
		std::cout << "13313441321414123" << std::endl;
		double testDoub = 13313441321414123;
		std::cout << std::setprecision(20) << testDoub << std::endl;
		std::cout << "13313441321414124" << std::endl;
		double testDoub2 = 13313441321414124;
		std::cout << std::setprecision(20) << testDoub2 << std::endl;
		std::cout << njh::colorBool(testDoub == testDoub2) << std::endl;
		std::cout << __FILE__ << " " << __LINE__ << std::endl;
	}
	{
		std::cout << "13313441321414125" << std::endl;
		double testDoub = 13313441321414125;
		std::cout << std::setprecision(20) << testDoub << std::endl;
		std::cout << "13313441321414126" << std::endl;
		double testDoub2 = 13313441321414126;
		std::cout << std::setprecision(20) << testDoub2 << std::endl;
		std::cout << njh::colorBool(testDoub == testDoub2) << std::endl;
		std::cout << __FILE__ << " " << __LINE__ << std::endl;
	}
	{
		std::cout << "13313441321414127" << std::endl;
		double testDoub = 13313441321414127;
		std::cout << std::setprecision(20) << testDoub << std::endl;
		std::cout << "13313441321414128" << std::endl;
		double testDoub2 = 13313441321414128;
		std::cout << std::setprecision(20) << testDoub2 << std::endl;
		std::cout << njh::colorBool(testDoub == testDoub2) << std::endl;
		std::cout << __FILE__ << " " << __LINE__ << std::endl;
	}


	std::cout << std::endl;

	{
		std::string testStr = "3313441321414121";
		std::cout << "places: " << testStr.size() << std::endl;
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "3313441321414122";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "3313441321414123";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "3313441321414124";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "3313441321414125";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "3313441321414126";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "3313441321414127";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "3313441321414129";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "3313441321414129";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}

	std::cout << __FILE__ << " " << __LINE__ << std::endl;

	{
		std::string testStr = "13313441321414121";
		std::cout << "places: " << testStr.size() << std::endl;
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "13313441321414122";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "13313441321414123";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "13313441321414124";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "13313441321414125";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "13313441321414126";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "13313441321414127";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "13313441321414128";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "13313441321414129";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	std::cout << __FILE__ << " " << __LINE__ << std::endl;

	std::cout << std::endl;

	{
		std::string testStr = "113313441321414121";
		std::cout << "places: " << testStr.size() << std::endl;
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "113313441321414122";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "113313441321414123";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "113313441321414124";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "113313441321414125";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "113313441321414126";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "113313441321414127";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "113313441321414128";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "113313441321414129";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}



	std::cout << __FILE__ << " " << __LINE__ << std::endl;

	{
		std::string testStr = "1213313441321414121";
		std::cout << "places: " << testStr.size() << std::endl;
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "1213313441321414122";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "1213313441321414123";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "1213313441321414124";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "1213313441321414125";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "1213313441321414126";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "1213313441321414127";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "1213313441321414128";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "1213313441321414129";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << __FILE__ << " " << __LINE__ << std::endl;


	SimpleKmerHash hashifier;

	{
		std::string kmerTest = "AGGAGTTAGCATATACN";
		uint64_t hash = hashifier.hash(kmerTest);
		std::string hashback = hashifier.reverseHash(hash);
		std::cout << "hash    : " << hash << std::endl;
		std::cout << "original: " << kmerTest << std::endl;
		std::cout << "convert : " << hashback << std::endl;

	}

	return 0;

}



int kmerExpRunner::findUniqKmersBetweenSeqSetsMultiDev(const njh::progutils::CmdArgs & inputCommands){
	KmerGatherer::KmerGathererPars countPars;
	bfs::path seqSetTableFnp = "";
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.description_ = "Get kmers that appear within all the sequences supplied and are unique to that set";
	setUp.setOption(seqSetTableFnp, "--seqSetTableFnp", "Seq Set Table, 2 columns, 1)set,2)fasta", true);
	countPars.setOptions(setUp);

	setUp.processDirectoryOutputName(njh::pasteAsStr(bfs::basename(seqSetTableFnp), "_TODAY"), true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);

	KmerGatherer kGather(countPars);

	std::unordered_map<std::string, std::set<std::string>> fastasForSet;
	table input(seqSetTableFnp, "\t", true);
	input.checkForColumnsThrow(VecStr{"set", "2bit"}, __PRETTY_FUNCTION__);
	for(const auto & row : input){
		fastasForSet[row[input.getColPos("set")]].emplace(row[input.getColPos("2bit")]);
	}
	setUp.rLog_.setCurrentLapName("initial");

	std::vector<bfs::path> twoBitFiles;
	for(const auto & seqSet : fastasForSet){
		for(const auto & fnp : seqSet.second){
			twoBitFiles.emplace_back(fnp);
		}
	}
	std::map<std::string, std::set<uint64_t>> kmersPerSet;

	std::function<bool(const std::string&)> seqCheck = [&countPars](const std::string & k){
		return std::all_of(k.begin(), k.end(), [&countPars](char base){return njh::in(base, countPars.allowableCharacters_);});
	};


	{
		setUp.rLog_.logCurrentTime("count_all");
		setUp.rLog_.runLogFile_.flush();
		auto allKmers = kGather.getUniqueKmersSetHash(twoBitFiles);
		setUp.rLog_.logCurrentTime("condense");
		setUp.rLog_.runLogFile_.flush();
		njh::concurrent::LockableQueue<std::string> seqSetNamesQueue(getVectorOfMapKeys(fastasForSet));
		for(const auto & name : fastasForSet){
			kmersPerSet[name.first] = std::set<uint64_t>{};
		}
		std::function<void()> condenseKmers = [&seqSetNamesQueue,&allKmers,&fastasForSet,&kmersPerSet,&seqCheck](){
			std::string name;
			while(seqSetNamesQueue.getVal(name)){
				SimpleKmerHash hasher;
				for(const auto & twobit : fastasForSet.at(name)){
					for(const auto & k : allKmers.at(twobit)){
						if(seqCheck(hasher.reverseHash(k))){
							kmersPerSet[name].emplace(k);
						}
					}
				}
			}
		};
		njh::concurrent::runVoidFunctionThreaded(condenseKmers, countPars.numThreads_);

	}
	std::map<std::string, std::set<uint64_t>> uniqueKmersFinal;
	setUp.rLog_.logCurrentTime("compare");
	setUp.rLog_.runLogFile_.flush();
	for(const auto & kmersForSet : kmersPerSet){
		uniqueKmersFinal[kmersForSet.first] = std::set<uint64_t>{};
	}
	{
		njh::concurrent::LockableQueue<std::string> seqSetNamesQueue(getVectorOfMapKeys(kmersPerSet));
		std::function<void()> compareKmers = [&seqSetNamesQueue,&kmersPerSet,&uniqueKmersFinal](){
			std::string name;
			while(seqSetNamesQueue.getVal(name)){
				std::set<uint64_t> uniqueKmers;
				uint32_t count = 0;
				for(const auto & otherSet : kmersPerSet){
					if(otherSet.first == name){
						continue;
					}
					if(0 == count){
						std::vector<uint64_t> notShared;
						std::set_difference(
								kmersPerSet.at(name).begin(), kmersPerSet.at(name).end(),
								otherSet.second.begin(), otherSet.second.end(),
								std::back_inserter(notShared));
						uniqueKmers = njh::vecToSet(notShared);
					}else{
						std::vector<uint64_t> notShared;
						std::set_difference(
								uniqueKmers.begin(), uniqueKmers.end(),
								otherSet.second.begin(), otherSet.second.end(),
								std::back_inserter(notShared));
						uniqueKmers = njh::vecToSet(notShared);
					}
					++count;
				}
				uniqueKmersFinal[name] = uniqueKmers;
			}
		};
		njh::concurrent::runVoidFunctionThreaded(compareKmers, countPars.numThreads_);
	}

	SimpleKmerHash hasher;
	OutputStream out(njh::files::make_path(setUp.pars_.directoryName_, "uniqueKmers.tab.txt"));
	for(const auto & kmersForSet : uniqueKmersFinal){
		for(const auto & kmer : kmersForSet.second){
			out << kmersForSet.first
					<< "\t" << hasher.reverseHash(kmer) << "\n";
		}
	}
	return 0;
}

int kmerExpRunner::findUniqKmersBetweenSeqSetsMulti(const njh::progutils::CmdArgs & inputCommands){
	KmerGatherer::KmerGathererPars countPars;
	bfs::path seqSetTableFnp = "";
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.description_ = "Get kmers that appear within all the sequences supplied and are unique to that set";
	setUp.setOption(seqSetTableFnp, "--seqSetTableFnp", "Seq Set Table, 2 columns, 1)set,2)fasta", true);
	countPars.setOptions(setUp);

	setUp.processDirectoryOutputName(njh::pasteAsStr(bfs::basename(seqSetTableFnp), "_TODAY"), true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);

	KmerGatherer kGather(countPars);

	std::unordered_map<std::string, std::set<std::string>> fastasForSet;
	table input(seqSetTableFnp, "\t", true);
	input.checkForColumnsThrow(VecStr{"set", "2bit"}, __PRETTY_FUNCTION__);
	for(const auto & row : input){
		fastasForSet[row[input.getColPos("set")]].emplace(row[input.getColPos("2bit")]);
	}
	setUp.rLog_.setCurrentLapName("initial");

	std::vector<bfs::path> twoBitFiles;
	for(const auto & seqSet : fastasForSet){
		for(const auto & fnp : seqSet.second){
			twoBitFiles.emplace_back(fnp);
		}
	}
	std::map<std::string, std::set<std::string>> kmersPerSet;

	std::function<bool(const std::string&)> seqCheck = [&countPars](const std::string & k){
		return std::all_of(k.begin(), k.end(), [&countPars](char base){return njh::in(base, countPars.allowableCharacters_);});
	};


	{
		setUp.rLog_.logCurrentTime("count_all");
		auto allKmers = kGather.getUniqueKmersSet(twoBitFiles);
		setUp.rLog_.logCurrentTime("condense");
		std::mutex setMut;
		njh::concurrent::LockableQueue<std::string> seqSetNamesQueue(getVectorOfMapKeys(fastasForSet));
		for(const auto & name : fastasForSet){
			kmersPerSet[name.first] = std::set<std::string>{};
		}
		std::function<void()> condenseKmers = [&seqSetNamesQueue,&allKmers,&fastasForSet,&kmersPerSet,&seqCheck](){
			std::string name;
			while(seqSetNamesQueue.getVal(name)){

				for(const auto & twobit : fastasForSet.at(name)){
					for(const auto & k : allKmers.at(twobit)){
						if(seqCheck(k)){
							kmersPerSet[name].emplace(k);
						}
					}
				}
			}
		};
		njh::concurrent::runVoidFunctionThreaded(condenseKmers, countPars.numThreads_);

	}
		std::map<std::string, std::set<std::string>> uniqueKmersFinal;

		for(const auto & kmersForSet : kmersPerSet){
			std::set<std::string> uniqueKmers;
			uint32_t count = 0;
			for(const auto & otherSet : kmersPerSet){
				if(otherSet.first == kmersForSet.first){
					continue;
				}
				if(0 == count){
					std::vector<std::string> notShared;
					std::set_difference(
							kmersForSet.second.begin(), kmersForSet.second.end(),
							otherSet.second.begin(), otherSet.second.end(),
							std::back_inserter(notShared));
					uniqueKmers = njh::vecToSet(notShared);
				}else{
					std::vector<std::string> notShared;
					std::set_difference(
							uniqueKmers.begin(), uniqueKmers.end(),
							otherSet.second.begin(), otherSet.second.end(),
							std::back_inserter(notShared));
					uniqueKmers = njh::vecToSet(notShared);
				}
				++count;
			}
			uniqueKmersFinal[kmersForSet.first] = uniqueKmers;
		}
		OutputStream out(njh::files::make_path(setUp.pars_.directoryName_, "uniqueKmers.tab.txt"));
		for(const auto & kmersForSet : uniqueKmersFinal){
			for(const auto & kmer : kmersForSet.second){
				out << kmersForSet.first
						<< "\t" << kmer << "\n";
			}
		}
	return 0;
}

int kmerExpRunner::findUniqKmersBetweenSeqSets(const njh::progutils::CmdArgs & inputCommands){
	KmerGatherer::KmerGathererPars countPars;
	bfs::path genomeFnp1 = "";
	bfs::path genomeFnp2 = "";
	bool writeOutAllCounts = false;

	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(genomeFnp1, "--genomeFnp1", "genomeFnp1", true);
	setUp.setOption(genomeFnp2, "--genomeFnp2", "genomeFnp2", true);
	countPars.setOptions(setUp);
	setUp.setOption(writeOutAllCounts, "--writeOutAllCounts", "writeOutAllCounts");
	std::string genomeFnp1Base = bfs::basename(njh::files::removeExtension(genomeFnp1));
	std::string genomeFnp2Base = bfs::basename(njh::files::removeExtension(genomeFnp2));

	setUp.processDirectoryOutputName(njh::pasteAsStr(genomeFnp1Base, "_vs_", genomeFnp2Base, "_TODAY"), true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);

	KmerGatherer kGather(countPars);

	setUp.rLog_.setCurrentLapName("initial");

	setUp.rLog_.logCurrentTime("genome1_counting-" + genomeFnp1Base);
	std::unordered_map<std::string, uint32_t> genome1Counts = kGather.countGenomeKmers(genomeFnp1);
	setUp.rLog_.logCurrentTime("genome2_counting-" + genomeFnp2Base);
	std::unordered_map<std::string, uint32_t> genome2Counts = kGather.countGenomeKmers(genomeFnp2);


	if(writeOutAllCounts){
		setUp.rLog_.logCurrentTime("genome1_writing_all-" + genomeFnp1Base);

		OutputStream genome1CountsOut(njh::files::make_path(setUp.pars_.directoryName_, genomeFnp1Base + "_counts.tab.txt.gz"));
		for(const auto & count : genome1Counts){
			genome1CountsOut << count.first << "\t" << count.second << "\n";
		}
		setUp.rLog_.logCurrentTime("genome2_writing_all-" + genomeFnp2Base);

		OutputStream genome2CountsOut(njh::files::make_path(setUp.pars_.directoryName_, genomeFnp2Base + "_counts.tab.txt.gz"));
		for(const auto & count : genome2Counts){
			genome2CountsOut << count.first << "\t" << count.second << "\n";
		}
	}


	std::function<bool(char)> charCheck = [&countPars](char base){
		return njh::in(base, countPars.allowableCharacters_);
	};
	setUp.rLog_.logCurrentTime("genome1_unique_deterination-" + genomeFnp1Base);

	OutputStream genome1UniqueCountsOut(njh::files::make_path(setUp.pars_.directoryName_, genomeFnp1Base + "_uniqueKmers_counts.tab.txt.gz"));
	for(const auto & count : genome1Counts){
		if(std::all_of(count.first.begin(), count.first.end(), charCheck) && !njh::in(count.first, genome2Counts) ){
			genome1UniqueCountsOut << count.first << "\t" << count.second << "\n";
		}
	}
	setUp.rLog_.logCurrentTime("genome2_unique_deterination-" + genomeFnp2Base);
	OutputStream genome2UniqueCountsOut(njh::files::make_path(setUp.pars_.directoryName_, genomeFnp2Base + "_uniqueKmers_counts.tab.txt.gz"));
	for(const auto & count : genome2Counts){
		if(std::all_of(count.first.begin(), count.first.end(), charCheck) && !njh::in(count.first, genome1Counts) ){
			genome2UniqueCountsOut << count.first << "\t" << count.second << "\n";
		}
	}
	setUp.rLog_.logCurrentTime("end");

	return 0;
}

}  //namespace njhseq

