#pragma once

// Created by Nicholas Hathaway on 5/25/23.
/* 
    
*/


#include <njhcpp/common.h>
#include <njhseq/objects/kmer.h>
#include <njhseq/objects/dataContainers/tables/TableReader.hpp>
#include <njhseq/IO/SeqIO/MultiSeqIO.hpp>


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

	struct CompareReadToSetPars {
		bool includeRevComp = false;
		uint32_t klen = 19;
		std::string sampleName;
		bool pairsSeparate = false;
		uint32_t hardCountOff = 0;
		double fracCutOff = 0;

		std::set<std::string> excludeSetNames{"genomeRest"};

		uint32_t initialExcludeHardCountOff = 0;
		double initialExcludeFracCutOff = 0;

		uint32_t kmerLengthForEntropyCalc_{1};
		double entropyFilter_{1.95};
	};

	struct CompareReadToSetRes {
		std::unordered_map<uint64_t, uint64_t> hashedInputKmers;
		std::unordered_map<uint64_t, uint64_t> hashedInputKmersRevComp;
		std::unordered_map<std::string, uint32_t> foundPerSet;
		std::unordered_map<std::string, uint32_t> foundPerSetRevComp;

		std::string winnerSet = "undetermined";
		double bestFrac = 0;
		bool winnerRevComp = false;

		uint32_t getTotalDetermined();

		void zeroFillFoundSets(const VecStr &setNames);
		static VecStr genOutputHeader(const CompareReadToSetPars &compPars);
		static void writeOutputHeader(std::ostream & out, const CompareReadToSetPars &compPars, const std::string & delim = "\t");

		std::vector<std::vector<std::string>> genOutput(const seqInfo &seq,
										 const std::unordered_map<std::string, std::unordered_set<uint64_t>> &uniqueKmersPerSet,
										 const CompareReadToSetPars &pars);

		void writeOutput(std::ostream & out, const seqInfo &seq,
										 const std::unordered_map<std::string, std::unordered_set<uint64_t>> &uniqueKmersPerSet,
										 const CompareReadToSetPars &pars, const std::string & delim = "\t");

	};


	static CompareReadToSetRes compareReadToSetRes(PairedRead &pseq,
																								 const std::unordered_map<std::string, std::unordered_set<uint64_t>> &uniqueKmersPerSet,
																								 const CompareReadToSetPars &compPars, const SimpleKmerHash &hasher);

	static CompareReadToSetRes compareReadToSetRes(seqInfo &seq,
																								 const std::unordered_map<std::string, std::unordered_set<uint64_t>> &uniqueKmersPerSet,
																								 const CompareReadToSetPars &compPars, const SimpleKmerHash &hasher);



	struct ProcessReadForExtractingPars {
		CompareReadToSetPars compPars;
		uint32_t smallLenCutOff = 0;
		bool writeOutExclude = false;
		bool doNotWriteUndetermined = false;
		bool doReCheckExcludeSets = false;
		uint32_t addingInKmersCountCutOff = 3;

		bool markReadsPerIteration = false;
		bool writeOutFinalKmerSets = false;
	};

	struct ProcessReadForExtractingCounts{
//		std::unordered_map<std::string, uint32_t> readsPerSet;
//		std::unordered_map<std::string, uint32_t> readsPerSetRevComp;
		ProcessReadForExtractingCounts(){
			readCountsPerSet[false] = std::unordered_map<std::string, uint32_t>{};
			readCountsPerSet[true] = std::unordered_map<std::string, uint32_t>{};
		}
		std::unordered_map<bool, std::unordered_map<std::string, uint32_t>> readCountsPerSet;
		uint32_t smallLenCutOffCount = 0;

		uint32_t filteredDissimilarCount = 0;//not to be included final counts

		void addOtherCounts(const ProcessReadForExtractingCounts & otherCounts);

		[[nodiscard]] uint64_t getTotalCounts() const;
		[[nodiscard]] uint64_t genTotalUndeterminedCount() const;
		[[nodiscard]] uint64_t genTotalDeterminedCount() const;
		static VecStr genOutCountsHeader(const ProcessReadForExtractingPars & extractingPars);
		static void writeOutCountsHeader(std::ostream & out, const ProcessReadForExtractingPars & extractingPars, const std::string & delim = "\t");


		[[nodiscard]] std::vector<std::vector<std::string>> genOutCounts(const ProcessReadForExtractingPars & extractingPars,
												const std::unordered_map<std::string, std::unordered_set<uint64_t>> &uniqueKmersPerSet,
												const std::string & iterName) const;

		void writeOutCounts(std::ostream & out, const ProcessReadForExtractingPars & extractingPars,
											const std::unordered_map<std::string, std::unordered_set<uint64_t>> &uniqueKmersPerSet,
												const std::string & iterName, const std::string & delim = "\t") const;



	};

	static void processReadForExtracting(PairedRead &pseq,
																			 const std::unordered_map<std::string, std::unordered_set<uint64_t>> &uniqueKmersPerSet,
																			 const ProcessReadForExtractingPars &extractingPars, const SimpleKmerHash &hasher,
																			 MultiSeqIO &seqOut, ProcessReadForExtractingCounts &counts,
																			 const std::string & iterationName);

	static void processReadForFilteringPairsSeparate(PairedRead &pseq,
																										const std::unordered_map<std::string, std::unordered_set<uint64_t>> &uniqueKmersPerSet,
																										const ProcessReadForExtractingPars &extractingPars,
																										const SimpleKmerHash &hasher,
																										MultiSeqIO &seqOut,
																										ProcessReadForExtractingCounts &counts,
																										const std::string & iterationName,
																										const std::string & filterName,
																										const std::string & keepName);

	static void processReadForFilteringPairsTogether(PairedRead &pseq,
																									 const std::unordered_map<std::string, std::unordered_set<uint64_t>> &uniqueKmersPerSet,
																									 const ProcessReadForExtractingPars &extractingPars,
																									 const SimpleKmerHash &hasher,
																									 MultiSeqIO &seqOut,
																									 ProcessReadForExtractingCounts &counts,
																									 const std::string & iterationName,
																									 const std::string & filterName,
																									 const std::string & keepName);

	static void processReadForFiltering(seqInfo &seq,
																									 const std::unordered_map<std::string, std::unordered_set<uint64_t>> &uniqueKmersPerSet,
																									 const ProcessReadForExtractingPars &extractingPars,
																									 const SimpleKmerHash &hasher,
																									 MultiSeqIO &seqOut,
																									 ProcessReadForExtractingCounts &counts,
																									 const std::string & iterationName,
																									 const std::string & filterName,
																									 const std::string & keepName);

	static void processReadForExtractingPairsSeparate(PairedRead &pseq,
																			 const std::unordered_map<std::string, std::unordered_set<uint64_t>> &uniqueKmersPerSet,
																			 const ProcessReadForExtractingPars &extractingPars, const SimpleKmerHash &hasher,
																			 MultiSeqIO &seqOut, ProcessReadForExtractingCounts &counts,
																			 const std::string & iterationName);

	static void processReadForExtractingPairsTogether(PairedRead &pseq,
																			 const std::unordered_map<std::string, std::unordered_set<uint64_t>> &uniqueKmersPerSet,
																			 const ProcessReadForExtractingPars &extractingPars, const SimpleKmerHash &hasher,
																			 MultiSeqIO &seqOut, ProcessReadForExtractingCounts &counts,
																			 const std::string & iterationName);

	static void processReadForExtracting(seqInfo &seq,
																			 const std::unordered_map<std::string, std::unordered_set<uint64_t>> &uniqueKmersPerSet,
																			 const ProcessReadForExtractingPars &extractingPars, const SimpleKmerHash &hasher,
																			 MultiSeqIO &seqOut, ProcessReadForExtractingCounts &counts,
																			 const std::string & iterationName);


	static std::unordered_map<std::string, std::unordered_map<uint64_t, uint32_t>>
	readInNewKmersFromExtractedReads(const bfs::path &directoryIn,
																	 const VecStr &kmerSets,
																	 const ProcessReadForExtractingPars &extractingPars);
	struct FilePositons {
		FilePositons() = default;
		std::ifstream::pos_type r1FnpEnd = 0;
		std::ifstream::pos_type r2FnpEnd = 0;
		std::ifstream::pos_type singleFnpEnd = 0;
	};
	static std::unordered_map<std::string, std::unordered_map<uint64_t, uint32_t>>
	readInNewKmersFromExtractedReads(const bfs::path &directoryIn,
																	 const VecStr &kmerSets,
																	 const ProcessReadForExtractingPars &extractingPars,
																	 const std::unordered_map<std::string, UniqueKmerSetHelper::FilePositons> & positionAfterLastIteration);

	static std::unordered_map<std::string, std::unordered_set<uint64_t>> filterReExtractedKmersForNonUnique(
					const std::unordered_map<std::string, std::unordered_map<uint64_t, uint32_t>> &rawKmersPerInput,
					const ProcessReadForExtractingPars &extractingPars,
					const std::unordered_map<std::string, std::unordered_set<uint64_t>> &uniqueKmersPerSet,
					std::unordered_set<uint64_t> &nonUniqueKmersPerSet) ;

	static std::unordered_map<std::string, std::unordered_set<uint64_t>> filterReExtractedKmersForNonUniqueIncludeExcludedSets(
					const std::unordered_map<std::string, std::unordered_map<uint64_t, uint32_t>> &rawKmersPerInput,
					const ProcessReadForExtractingPars &extractingPars,
					const std::unordered_map<std::string, std::unordered_set<uint64_t>> &uniqueKmersPerSet,
					std::unordered_set<uint64_t> &nonUniqueKmersPerSet) ;


};



}  // namespace njhseq

