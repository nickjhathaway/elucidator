/*
 * RoughIlluminaSimulator.cpp
 *
 *  Created on: May 4, 2019
 *      Author: nicholashathaway
 */



#include "RoughIlluminaSimulator.hpp"

namespace njhseq {

RoughIlluminaSimulator::RoughIlluminaSimulator(const bfs::path & profileDir, double errorRateCorrection) :
		profileDir_(profileDir) {
	//check for files

	std::vector<bfs::path> inputProfileFiles{
		"r1_base_substitution_rates.tab.txt",
		"r1_positional_error_rate.tab.txt",
		"r1_quality_distribution_for_correct_calls.tab.txt",
		"r1_quality_distribution_for_error_calls.tab.txt",
		"r1_overhang_profile.tab.txt",
		"r2_base_substitution_rates.tab.txt",
		"r2_positional_error_rate.tab.txt",
		"r2_quality_distribution_for_correct_calls.tab.txt",
		"r2_quality_distribution_for_error_calls.tab.txt",
		"r2_overhang_profile.tab.txt"
	};

	std::vector<bfs::path> inputProfileFilesFull;
	for(const auto & fnp : inputProfileFiles){
		inputProfileFilesFull.emplace_back(njh::files::make_path(profileDir_, fnp));
	}
	njh::files::checkExistenceThrow(inputProfileFilesFull, __PRETTY_FUNCTION__);


	ReadProfile::ReadProfileFnps r1Fnps;
	r1Fnps.base_substitution_rates_fnp = njh::files::make_path(profileDir_, "r1_base_substitution_rates.tab.txt");
	r1Fnps.overhang_profile_fnp = njh::files::make_path(profileDir_, "r1_overhang_profile.tab.txt");
	r1Fnps.positional_error_rate_fnp = njh::files::make_path(profileDir_, "r1_positional_error_rate.tab.txt");
	r1Fnps.quality_distribution_for_correct_calls_fnp = njh::files::make_path(profileDir_, "r1_quality_distribution_for_correct_calls.tab.txt");
	r1Fnps.quality_distribution_for_error_calls_fnp = njh::files::make_path(profileDir_, "r1_quality_distribution_for_error_calls.tab.txt");
	r1Profile_ = ReadProfile(r1Fnps, errorRateCorrection);

	ReadProfile::ReadProfileFnps r2Fnps;
	r2Fnps.base_substitution_rates_fnp = njh::files::make_path(profileDir_, "r2_base_substitution_rates.tab.txt");
	r2Fnps.overhang_profile_fnp = njh::files::make_path(profileDir_, "r2_overhang_profile.tab.txt");
	r2Fnps.positional_error_rate_fnp = njh::files::make_path(profileDir_, "r2_positional_error_rate.tab.txt");
	r2Fnps.quality_distribution_for_correct_calls_fnp = njh::files::make_path(profileDir_, "r2_quality_distribution_for_correct_calls.tab.txt");
	r2Fnps.quality_distribution_for_error_calls_fnp = njh::files::make_path(profileDir_, "r2_quality_distribution_for_error_calls.tab.txt");
	r2Profile_ = ReadProfile(r2Fnps, errorRateCorrection);

}



RoughIlluminaSimulator::ReadProfile::ReadProfile(){

}

RoughIlluminaSimulator::ReadProfile::ReadProfile(const ReadProfileFnps & fnps, double errorRateCorrection){
	//base sub r
	{
		table rBases(fnps.base_substitution_rates_fnp, "\t", true);
		rBases.checkForColumnsThrow(VecStr{"ref", "seq", "count"}, __PRETTY_FUNCTION__);
		std::unordered_map<char, std::pair<std::vector<char>, std::vector<uint32_t>>> baseCounts;
		for(const auto & row : rBases){
			baseCounts[row[rBases.getColPos("ref")].front()].first.emplace_back(row[rBases.getColPos("seq")].front());
			baseCounts[row[rBases.getColPos("ref")].front()].second.emplace_back(njh::StrToNumConverter::stoToNum<uint32_t>(row[rBases.getColPos("count")] ));
		}
		for(const auto & baseCount : baseCounts){
			errorBaseGen_.emplace(baseCount.first, njh::randObjectGen<char, uint32_t>(baseCount.second.first,baseCount.second.second ));
		}
	}
	//positional error rate
	{
		table rPositionError(fnps.positional_error_rate_fnp, "\t", true);
		rPositionError.checkForColumnsThrow(VecStr{"position", "errors", "total", "rate"}, __PRETTY_FUNCTION__);
		uint32_t lastPosition = 0;
		for(const auto & row : rPositionError){
			uint32_t position = njh::StrToNumConverter::stoToNum<uint32_t>(row[rPositionError.getColPos("position")]);
			if (errorRates_.empty()) {
				if (0 != position) {
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__
							<< "error, first position should be zero, not " << position
							<< "\n";
					throw std::runtime_error { ss.str() };
				}
			} else {
				if (lastPosition + 1 != position) {
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__
							<< "error, position should be one greater than the last position, position: "
							<< position << ", lastPosition: " << lastPosition << "\n";
					throw std::runtime_error { ss.str() };
				}
			}
			//should make sure this is between 0 and 1...
			errorRates_.emplace_back(njh::StrToNumConverter::stoToNum<double>(row[rPositionError.getColPos("rate")]));
			for(auto & errorRate : errorRates_){
				errorRate = std::max(errorRate - errorRateCorrection, 0.0);
			}
			lastPosition = position;
		}
	}
	//positional qual calls
	{

		table rCorrectQuals(fnps.quality_distribution_for_correct_calls_fnp, "\t", true);
		rCorrectQuals.checkForColumnsThrow(VecStr{"position", "quality", "count"}, __PRETTY_FUNCTION__);
		std::unordered_map<uint32_t, std::pair<std::vector<uint32_t>, std::vector<uint32_t>>> qualPerPositionCalls;
		for(const auto & row : rCorrectQuals){
			uint32_t position = njh::StrToNumConverter::stoToNum<uint32_t>(row[rCorrectQuals.getColPos("position")]);
			qualPerPositionCalls[position].first.emplace_back(njh::StrToNumConverter::stoToNum<uint32_t>(row[rCorrectQuals.getColPos("quality")]));
			qualPerPositionCalls[position].second.emplace_back(njh::StrToNumConverter::stoToNum<uint32_t>(row[rCorrectQuals.getColPos("count")]));
		}
		for (const auto pos : iter::range(errorRates_.size())) {
			if (!njh::in(pos, qualPerPositionCalls)) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << "error, pos " << pos
						<< " missing from quality dist in "
						<< fnps.quality_distribution_for_correct_calls_fnp << "\n";
				throw std::runtime_error { ss.str() };
			}
			regularQualityDist_.emplace_back(
					njh::randObjectGen<uint32_t, uint32_t>(
							qualPerPositionCalls[pos].first,
							qualPerPositionCalls[pos].second));
		}
	}

	//positional error qual calls
	{
		table rErrorQuals(fnps.quality_distribution_for_error_calls_fnp, "\t", true);
		rErrorQuals.checkForColumnsThrow(VecStr{"position", "quality", "count"}, __PRETTY_FUNCTION__);
		std::unordered_map<uint32_t, std::pair<std::vector<uint32_t>, std::vector<uint32_t>>> qualPerPositionCalls;
		for(const auto & row : rErrorQuals){
			uint32_t position = njh::StrToNumConverter::stoToNum<uint32_t>(row[rErrorQuals.getColPos("position")]);
			qualPerPositionCalls[position].first.emplace_back(njh::StrToNumConverter::stoToNum<uint32_t>(row[rErrorQuals.getColPos("quality")]));
			qualPerPositionCalls[position].second.emplace_back(njh::StrToNumConverter::stoToNum<uint32_t>(row[rErrorQuals.getColPos("count")]));
		}
		for (const auto pos : iter::range(errorRates_.size())) {
			if (!njh::in(pos, qualPerPositionCalls)) {
				if(errorQualityDist_.empty()){
					if(qualPerPositionCalls.empty()){
						errorQualityDist_.emplace_back(
								njh::randObjectGen<uint32_t, uint32_t>(
										{15,16,17,18,19},
										{1,1,1,1,1}));
					}else{
						auto maxPosition = vectorMaximum(getVectorOfMapKeys(qualPerPositionCalls));
						if(pos > maxPosition){
							errorQualityDist_.emplace_back(
														njh::randObjectGen<uint32_t, uint32_t>(
																qualPerPositionCalls[maxPosition].first,
																qualPerPositionCalls[maxPosition].second));
						}else{
							uint32_t nextPos = pos + 1;
							while(!njh::in(nextPos, qualPerPositionCalls)){
								++nextPos;
							}
							errorQualityDist_.emplace_back(
														njh::randObjectGen<uint32_t, uint32_t>(
																qualPerPositionCalls[nextPos].first,
																qualPerPositionCalls[nextPos].second));
						}
					}
				}else{
					errorQualityDist_.emplace_back(errorQualityDist_.back());
				}
			}else{
				errorQualityDist_.emplace_back(
						njh::randObjectGen<uint32_t, uint32_t>(
								qualPerPositionCalls[pos].first,
								qualPerPositionCalls[pos].second) );
			}
		}
	}

	//overhang generator
	{
		table rOverhangProfile(fnps.overhang_profile_fnp, "\t", true);
		rOverhangProfile.checkForColumnsThrow(VecStr{"pos", "char", "qual", "count"}, __PRETTY_FUNCTION__);
		std::unordered_map<uint32_t, std::unordered_map<char, std::unordered_map<uint32_t, uint32_t>>> counts;

		for(const auto & row : rOverhangProfile){

			uint32_t position = njh::StrToNumConverter::stoToNum<uint32_t>(row[rOverhangProfile.getColPos("pos")]);
			char base = row[rOverhangProfile.getColPos("char")].front();
			uint32_t qual = njh::StrToNumConverter::stoToNum<uint32_t>(row[rOverhangProfile.getColPos("qual")]);
			uint32_t count = njh::StrToNumConverter::stoToNum<uint32_t>(row[rOverhangProfile.getColPos("count")]);
			counts[position][base][qual] = count;
		}
		auto maxPosition = vectorMaximum(getVectorOfMapKeys(counts));
		for(uint32_t pos = 0; pos < maxPosition; ++pos){
			if(!njh::in(pos, counts)){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << "error, pos " << pos
						<< " missing from overhang profile "
						<< fnps.overhang_profile_fnp << "\n";
				throw std::runtime_error { ss.str() };
			}
		}
		for(uint32_t pos = 0; pos < maxPosition; ++pos){
			std::unordered_map<char, uint32_t> baseCounts;
			std::unordered_map<char, njh::randObjectGen<uint32_t, uint32_t>> qualForBase;
			for(const auto & base : counts[pos]){
				qualForBase.emplace(base.first, njh::randObjectGen<uint32_t, uint32_t>(getVectorOfMapKeys(base.second), getVectorOfMapValues(base.second)));
				for(const auto & qual : base.second){
					baseCounts[base.first] += qual.second;
				}
			}
			overHangBaseGen_.emplace_back(getVectorOfMapKeys(baseCounts), getVectorOfMapValues(baseCounts));
			overHangQualForBaseGen_.emplace_back(qualForBase);
		}
	}
}

seqInfo RoughIlluminaSimulator::ReadProfile::simRead(seqInfo input, uint32_t length){ //const
	if(length > errorRates_.size()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error length is longer than can be simulated" << "\n";
		throw std::runtime_error{ss.str()};
	}
	for(const auto pos : iter::range(len(input))){
		//decide if error
		if(rGen_() <= errorRates_[pos]){
			//error
			input.seq_[pos] = errorBaseGen_.at(input.seq_[pos]).genObj();
			input.qual_[pos] = errorQualityDist_[pos].genObj();
		}else{
			//no error
			input.qual_[pos] = regularQualityDist_[pos].genObj();
		}
	}
	//add overhang
	if(length > len(input)){
		for(uint32_t pos = len(input); pos < length; ++pos){
			uint32_t overHangPos = pos - len(input);
			if(overHangPos > overHangBaseGen_.size()){
				char overhangBase = overHangBaseGen_.back().genObj();
				uint32_t qual = overHangQualForBaseGen_.back().at(overhangBase).genObj();
				input.append(overhangBase, qual);
			}else{
				char overhangBase = overHangBaseGen_[overHangPos].genObj();
				uint32_t qual = overHangQualForBaseGen_[overHangPos].at(overhangBase).genObj();
				input.append(overhangBase, qual);
			}
		}
	}
	return input;
}

seqInfo RoughIlluminaSimulator::simR1(const seqInfo & input, uint32_t r1Length){ //const
	return r1Profile_.simRead(input, r1Length);
}

seqInfo RoughIlluminaSimulator::simR2(const seqInfo & input, uint32_t r2Length){ //const
	return r2Profile_.simRead(input, r2Length);
}

}  // namespace njhseq
