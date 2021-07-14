#pragma once

/*
 * SimpleKmerHash.hpp
 *
 *  Created on: Jul 13, 2021
 *      Author: nick
 */

#include <njhseq/common.h>


namespace njhseq {

class SimpleKmerHash{

public:
	SimpleKmerHash();
	std::vector<char> hasher_;
	std::vector<char> revCompHasher_;
	std::vector<char> reverseHasher_;
	std::vector<char> revCompReverseHasher_;


	uint64_t hash(const std::string & str) const;
	std::string reverseHash(uint64_t hash) const;

	uint64_t revCompHash(const std::string & str) const;
	std::string revCompReverseHash(uint64_t hash) const;
};

}  // namespace njhseq




