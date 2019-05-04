#pragma once

/*
 * RoughIlluminaSimulator.hpp
 *
 *  Created on: May 4, 2019
 *      Author: nicholashathaway
 */

#include <njhcpp.h>
#include <njhseq/IO.h>


namespace njhseq {




class RoughIlluminaSimulator {
public:

	RoughIlluminaSimulator(const bfs::path & profileDir, double errorRateCorrection = 0);

	bfs::path profileDir_;

	struct ReadProfile{

		struct ReadProfileFnps {
			bfs::path base_substitution_rates_fnp;
			bfs::path positional_error_rate_fnp;
			bfs::path quality_distribution_for_correct_calls_fnp;
			bfs::path quality_distribution_for_error_calls_fnp;
			bfs::path overhang_profile_fnp;
		};
		ReadProfile();

		ReadProfile(const ReadProfileFnps & fnps, double errorRateCorrection = 0);

		std::unordered_map<char, njh::randObjectGen<char, uint32_t>> errorBaseGen_;

		std::vector<double> errorRates_; //chance of error, position in vector is position of read;

		std::vector<njh::randObjectGen<uint32_t, uint32_t>> errorQualityDist_;
		std::vector<njh::randObjectGen<uint32_t, uint32_t>> regularQualityDist_;

		std::vector<njh::randObjectGen<char, uint32_t>> overHangBaseGen_;
		std::vector<std::unordered_map<char, njh::randObjectGen<uint32_t, uint32_t>>> overHangQualForBaseGen_;
		njh::randomGenerator rGen_;

		seqInfo simRead(seqInfo input, uint32_t length);
	};

	ReadProfile r1Profile_;
	ReadProfile r2Profile_;

	njh::randomGenerator rGen_;


	seqInfo simR1(const seqInfo & input, uint32_t r1Length);

	seqInfo simR2(const seqInfo & input, uint32_t r2Length);

};




}  // namespace njhseq




