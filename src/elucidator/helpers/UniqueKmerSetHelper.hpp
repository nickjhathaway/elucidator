#pragma once

// Created by Nicholas Hathaway on 5/25/23.
/* 
    
*/


#include <njhcpp/common.h>
#include <njhseq/objects/kmer.h>
#include <njhseq/objects/dataContainers/tables/TableReader.hpp>


namespace njhseq {


class UniqueKmerSetHelper {
public:

	/**
	 * @brief Get kmer length used to create unique kmer table
	 * @param uniqueKmerTableFnp the file path to the table of unique kmers, headerless column 1) set name, 2) kmer
	 * @return the kmer length of the table
	 *
	 * Expects the file to be a two column table with kmers in the second column, kmers can't be longer than 19 in length
	 *
	 */
	static uint32_t getKmerLenFromUniqueKmerTable(const bfs::path &uniqueKmerTableFnp);

	/**
	 * @brief Read in the kmers in the file and store thier hashed values in map, key=set name, value=std::set of hashed kmer
	 * @param uniqueKmerTableFnp the file path to the table of unique kmers, headerless column 1) set name, 2) kmer
	 * @return unordered map key=set name, value=std::set of hashed kmer
	 *
	 * Expects the file to be a two column table with kmers in the second column, kmers can't be longer than 19 in length
	 */
	static std::unordered_map<std::string, std::unordered_set<uint64_t>> readInUniqueKmerTablePerSet(const bfs::path &uniqueKmerTableFnp);

	/**
	 * @brief Read in the second column of colum of kmers and don't store the sets they are coming from
	 * @param uniqueKmerTableFnp the file path to the table of unique kmers, headerless column 1) set name, 2) kmer
	 * @return an std::unordered_set of kmers
	 *
	 * Expects the file to be a two column table with kmers in the second column, kmers can't be longer than 19 in length
	 */
	static std::unordered_set<uint64_t> readInUniqueKmerTableSetsCollapsed(const bfs::path &uniqueKmerTableFnp);

};



}  // namespace njhseq

