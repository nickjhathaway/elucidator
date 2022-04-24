/*
 * SimpleKmerHash.cpp
 *
 *  Created on: Jul 13, 2021
 *      Author: nick
 */


#include "SimpleKmerHash.hpp"


namespace njhseq {


SimpleKmerHash::SimpleKmerHash(){
	hasher_ = std::vector<char>(255, '5'); //by defaulting to 5 everything but A, C, G, T will be turned into an N
	hasher_['A'] = '1';
	hasher_['C'] = '2';
	hasher_['G'] = '3';
	hasher_['T'] = '4';
	hasher_['N'] = '5';

	revCompHasher_ = std::vector<char>(255, '5'); //by defaulting to 5 everything but A, C, G, T will be turned into an N
	revCompHasher_['A'] = '4';//T
	revCompHasher_['C'] = '3';//G
	revCompHasher_['G'] = '2';//C
	revCompHasher_['T'] = '1';//A
	revCompHasher_['N'] = '5';//stil N

	reverseHasher_ = std::vector(255, 'N');
	reverseHasher_['1'] = 'A';
	reverseHasher_['2'] = 'C';
	reverseHasher_['3'] = 'G';
	reverseHasher_['4'] = 'T';
	reverseHasher_['5'] = 'N';

	revCompReverseHasher_ = std::vector(255, 'N');
	revCompReverseHasher_['1'] = 'T';//1 = A, now T
	revCompReverseHasher_['2'] = 'G';//2 = C, now G
	revCompReverseHasher_['3'] = 'C';//3 = G, now C
	revCompReverseHasher_['4'] = 'A';//4 = T, now A
	revCompReverseHasher_['5'] = 'N';//still N
}



uint64_t SimpleKmerHash::hash(const std::string & str) const{
	std::string convert;
	for(size_t pos = 0; pos < std::min<size_t>(20, str.size()); ++pos){
		convert.push_back(hasher_[str[pos]]);
	}
	return njh::StrToNumConverter::stoToNum<uint64_t>(convert);
}

uint64_t SimpleKmerHash::hash(const std::string & str, uint32_t start, uint32_t size) const{
	std::string convert;

//	std::cout << "start + size: " << start + size << std::endl;
//	std::cout << "str.size(): " << str.size()  << std::endl;
	for(size_t pos = start; pos < std::min<size_t>(start + size, str.size()); ++pos){
		convert.push_back(hasher_[str[pos]]);
	}
//	std::cout << "convert: " << convert << std::endl;
	return njh::StrToNumConverter::stoToNum<uint64_t>(convert);
}




std::string SimpleKmerHash::reverseHash(uint64_t hash) const {
	std::string hashStr = estd::to_string(hash);
	std::string back;
	back.reserve(hashStr.size());
	for(const auto pos : iter::range(hashStr.size())){
		back.push_back(reverseHasher_[hashStr[pos]]);
	}
	return back;
}



uint64_t SimpleKmerHash::revCompHash(const std::string & str) const{
	std::string convert;
	//go over backwards to reverse complement
	for(size_t pos = std::min<size_t>(20, str.size()); pos >0 ; --pos){
		convert.push_back(revCompHasher_[str[pos - 1]]);
	}
	return njh::StrToNumConverter::stoToNum<uint64_t>(convert);
}

uint64_t SimpleKmerHash::revCompHash(const std::string & str, uint32_t start, uint32_t size) const{
	std::string convert;
	//go over backwards to reverse complement
	for(size_t pos = std::min<size_t>(start + size, str.size()); pos >start ; --pos){
		convert.push_back(revCompHasher_[str[pos - 1]]);
	}
	return njh::StrToNumConverter::stoToNum<uint64_t>(convert);
}

std::string SimpleKmerHash::revCompReverseHash(uint64_t hash) const {
	std::string hashStr = estd::to_string(hash);
	std::string back;
	back.reserve(hashStr.size());
	for(const auto pos : iter::range(hashStr.size())){
		//go from the back to reverse complement
		back.push_back(revCompReverseHasher_[hashStr[hashStr.size() - 1 - pos]]);
	}
	return back;
}


}  // namespace njhseq

