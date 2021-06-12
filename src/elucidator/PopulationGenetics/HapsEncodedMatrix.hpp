#pragma once

/*
 * HapsEncodedMatrix.hpp
 *
 *  Created on: Jun 11, 2021
 *      Author: nick
 */

#include "elucidator/common.h"

#include <njhseq/programUtils/seqSetUp.hpp>
#include <njhseq/objects/dataContainers/tables/TableReader.hpp>
#include <njhseq/concurrency/PairwisePairFactory.hpp>
#include <njhseq/objects/Meta/MultipleGroupMetaData.hpp>


namespace njhseq {



class HapsEncodedMatrix{

public:


	struct SetWithExternalPars{
		bfs::path tableFnp = "";

		std::string sampleCol = "";
		std::string targetNameCol = "";
		std::string popIDCol = "";
		std::string relAbundCol = "";
		std::unordered_set<std::string> selectSamples{};
		std::unordered_set<std::string> selectTargets{};
		void setDefaults(seqSetUp & setUp);
		uint32_t numThreads = 1;
		bool majorOnly = false;
	};

	HapsEncodedMatrix(const SetWithExternalPars & pars);

	SetWithExternalPars pars_;

	std::unordered_set<std::string> sampNames_; //! used during the encoding phase
	std::vector<std::string> sampNamesVec_; //! sample names, index is the key to the sample
	std::unordered_map<std::string, uint32_t> sampNamesKey_;//! key to sample name to sample index

	std::unordered_set<std::string> tarNames_; //! used during the encoding phase
	std::vector<std::string> tarNamesVec_; //! target names, index is the key to the targets
	std::unordered_map<std::string, uint32_t> tarNameKey_; //! key to target name to target index

	std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>> hapNamesForTars_; //! all haplotype names for each target
	std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>> hapNamesKey_; //! key for each haplotype within each target, will range from 0 to haplotype count for target


	std::vector<uint32_t> numberOfHapsPerTarget_; //! number of haplotypes per target, index is the target as denoted from tarNameKey_
	std::vector<uint32_t> tarStart_; //! the starting position of target in targetsEncodeBySamp_ columns

	std::vector<std::vector<uint8_t>> hapsEncodeBySamp_; //! each row is a sample, each column is hap, 0 for not present, 1 for present
	std::vector<std::vector<uint8_t>> targetsEncodeBySamp_; //! each row is a sample, each column is a target, 0 if sample has no data for target, 1 for has data

	std::vector<double> hapsProbs_;

	uint64_t totalHaps_{0};

	bool encodeKeysSet_{false};

	std::shared_ptr<MultipleGroupMetaData> meta_;


	struct SampTarHap {
		SampTarHap(const std::string &samp, const std::string &tar,
				const std::string &hap);
		std::string samp_;
		std::string tar_;
		std::string hap_;

	};
	void addSampTarHapForEncoding(const std::string &samp, const std::string &tar,
			const std::string &hap);

	void addSampTarHapForEncoding(const SampTarHap & adding);

	void encodeSampTarHap(const std::string &samp, const std::string &tar,
			const std::string &hap);

	void encodeSampTarHap(const SampTarHap & adding);


	void calcHapProbs();

	table getNumberTargetsPerSample() const;

	void addMeta(const bfs::path & metaFnp);

	void resetEncoding();
	void setEncodeKeys();

	struct IndexResults{
		IndexResults(const uint64_t numOfSamps);
		std::vector<std::vector<double>> byAllHaps; //! jacard index by all input haplotpes
		std::vector<std::vector<double>> byHapsTarShared;//! jacard index for haplotypes for targets where both samples have data
		std::vector<std::vector<double>> byHapsTarSharedWeighted;//! jacard index for haplotypes for targets where both samples have data

		std::vector<std::vector<double>> avgJacard; //! averaged jacard distance for shared targets
		std::vector<std::vector<double>> avgJacardWeighted; //! averaged jacard distance for shared targets weighted by
		std::vector<std::vector<double>> byTarget; //! fraction of targets that have at least shared haplotype between samples

	};

	IndexResults genIndexMeasures(bool verbose = false) const;

};




}  // namespace njhseq



