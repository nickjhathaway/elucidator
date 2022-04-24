#pragma once

/*
 * SimpleKmerHash.hpp
 *
 *  Created on: Jul 13, 2021
 *      Author: nick
 */

#include <njhseq/common.h>


namespace njhseq {

/**
 * @brief A simple DNA k-mer hasher, max up to 19mers, has added benefit of being a reversible hash
 */
class SimpleKmerHash{
	std::vector<char> hasher_;
	std::vector<char> revCompHasher_;
	std::vector<char> reverseHasher_;
	std::vector<char> revCompReverseHasher_;
public:
	SimpleKmerHash();

	/**
	 * @brief hash up to the first 19 characters of string, hash can only do 19 characters, no check on len is done but only up to the first 19 characters are hashed no matter the input size
	 * @param str the string to hash
	 * @return an uint64_t with the hash for the input str
	 */
	[[nodiscard]] uint64_t hash(const std::string & str) const;
	/**
	 * @brief reverse a the hash to the original string
	 * @param hash the hash to reverse
	 * @return the original kmer
	 */
	[[nodiscard]] std::string reverseHash(uint64_t hash) const;

	/**
	 * @brief create a hash at the same time as reverse complementing, this way the original string doesn't need to be reverse complement
	 * @param str a kmer that will be reverse complemented and hashed
	 * @return a uint64_t hash representing the reverse complement of the input str
	 */
	[[nodiscard]] uint64_t revCompHash(const std::string & str) const;

	/**
	 * @brief reverse a hash and reverse complement
	 * @param hash the hash to reverse and reverse complement
	 * @return the reverse complement of the reversed hash
	 */
	[[nodiscard]] std::string revCompReverseHash(uint64_t hash) const;
};

}  // namespace njhseq




