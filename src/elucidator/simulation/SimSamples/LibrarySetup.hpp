#pragma once

/*
 * LibrarySetup.hpp
 *
 *  Created on: May 2, 2019
 *      Author: nicholashathaway
 */

#include <SeekDeep/objects/PrimersAndMids.hpp>

#include "elucidator/common.h"
#include "elucidator/simulation/pcrSimulation/PCRSimulator.hpp"

namespace njhseq {




class ControlPopulation {
public:
	class Strain {
	public:
		Strain(const std::string & name,
				const std::map<std::string, std::string> & genomicRegionToHapNames) :
				name_(name), genomicRegionToHapNames_(genomicRegionToHapNames) {

		}
		Strain(const std::string & name,
				const std::map<std::string, std::string> & genomicRegionToHapNames,
				double relativeAbundance) :
				name_(name), genomicRegionToHapNames_(genomicRegionToHapNames), relativeAbundance_(
						relativeAbundance) {

		}
		std::string name_;
		std::map<std::string, std::string> genomicRegionToHapNames_;
		double relativeAbundance_ { 1.0 };
		MetaDataInName meta_;

		Json::Value toJson() const {
			Json::Value ret;
			ret["class"] = njh::json::toJson(njh::getTypeName(*this));
			ret["name_"] = njh::json::toJson(name_);
			ret["genomicRegionToHapNames_"] = njh::json::toJson(genomicRegionToHapNames_);
			ret["relativeAbundance_"] = njh::json::toJson(relativeAbundance_);
			ret["meta_"] = njh::json::toJson(meta_);
			return ret;
		}

		static VecStr genJsonMemNames() {
			return VecStr { "name_", "genomicRegionToHapNames_", "relativeAbundance_", "meta_" };
		}

		static Strain genFromJson(const Json::Value & val){
			njh::json::MemberChecker memCheck(val);
			memCheck.failMemberCheckThrow(genJsonMemNames(), __PRETTY_FUNCTION__);
			Strain ret(val["sampName_"].asString(),
					njh::json::JsonToMap<std::string, std::string>(val["genomicRegionToHapNames_"]),
							val["relativeAbundance_"].asDouble());
			ret.meta_ = MetaDataInName::genMetaFromJson(val["meta_"]);
			return ret;
		}
	};

	class GenomicRegionAmp {
	public:
		GenomicRegionAmp(const std::string & tarName,
				const std::map<std::string, double>& hapAbundSampled) :
				name_(tarName), hapAbundSampled_(hapAbundSampled) {

		}

		std::string name_;
		std::map<std::string, double> hapAbundSampled_;
		std::map<std::string, double> hapAbundSequenced_;
		std::map<std::string, double> hapAbundObserved_;

		MetaDataInName meta_;

		Json::Value toJson() const {
			Json::Value ret;
			ret["class"] = njh::json::toJson(njh::getTypeName(*this));
			ret["name_"] = njh::json::toJson(name_);
			ret["hapAbundSampled_"] = njh::json::toJson(hapAbundSampled_);
			ret["hapAbundSequenced_"] = njh::json::toJson(hapAbundSequenced_);
			ret["hapAbundObserved_"] = njh::json::toJson(hapAbundObserved_);
			ret["meta_"] = njh::json::toJson(meta_);
			return ret;
		}

		static VecStr genJsonMem() {
			return VecStr { "name_", "hapAbundSampled_", "hapAbundSequenced_", "hapAbundObserved_", "meta_" };
		}

		static GenomicRegionAmp genFromJson(const Json::Value & val){
			njh::json::MemberChecker memCheck(val);
			memCheck.failMemberCheckThrow(genJsonMem(), __PRETTY_FUNCTION__);
			GenomicRegionAmp ret(val["name_"].asString(), njh::json::JsonToMap<std::string, double>(val["hapAbundSampled_"]));
			ret.hapAbundSequenced_ = njh::json::JsonToMap<std::string, double>(val["hapAbundSequenced_"]);
			ret.hapAbundObserved_ = njh::json::JsonToMap<std::string, double>(val["hapAbundObserved_"]);
			ret.meta_ = MetaDataInName::genMetaFromJson(val["meta_"]);
			return ret;
		}
	};

	struct SeqSampledAmounts{
		uint32_t sequencedReadAmount_{0};
		uint32_t totalGenomesSampled_{0};

		Json::Value toJson() const {
			Json::Value ret;
			ret["class"] = njh::json::toJson(njh::getTypeName(*this));
			ret["sequencedReadAmount_"] = njh::json::toJson(sequencedReadAmount_);
			ret["totalGenomesSampled_"] = njh::json::toJson(totalGenomesSampled_);
			return ret;
		}
		static VecStr genJsonMem(){
			return VecStr{"sequencedReadAmount_", "totalGenomesSampled_"};
		}

		static SeqSampledAmounts genFromJson(const Json::Value & val){
			njh::json::MemberChecker memCheck(val);
			memCheck.failMemberCheckThrow(genJsonMem(), __PRETTY_FUNCTION__);
			SeqSampledAmounts ret;
			ret.sequencedReadAmount_ = val["sequencedReadAmount_"].asUInt();
			ret.totalGenomesSampled_ = val["totalGenomesSampled_"].asUInt();
			return ret;
		}
	};

	class Sample {
	public:
		Sample(const std::string & name,
				const std::map<std::string, Strain> & strainsExpected) :
				sampName_(name), strainsExpected_(strainsExpected) {

		}



		class ExperimentRun {
		public:

			ExperimentRun(const std::string & name, const SeqSampledAmounts & expAmounts) :
					runName_(name), expAmounts_(expAmounts) {

			}
			std::string runName_;
			std::map<std::string, Strain> hapAbundGenomesSampled_;
			std::map<std::string, GenomicRegionAmp> hapRegionAmplified_;
			SeqSampledAmounts expAmounts_;
			MetaDataInName meta_;

			void setGenomesSampled(const std::map<std::string, Strain> & strainsExpected){
				hapAbundGenomesSampled_.clear();
				hapRegionAmplified_.clear();
				std::vector<seqInfo> seqs;
				for(const auto & strain : strainsExpected){
					seqInfo seq(strain.second.name_);
					seq.cnt_ = strain.second.relativeAbundance_;
					seq.frac_ = strain.second.relativeAbundance_;
					seqs.emplace_back(seq);
				}
				auto genomesSampled = PCRSimulator::randomlySampleGenomes(seqs,
						expAmounts_.totalGenomesSampled_);
				std::map<std::string, std::map<std::string, double>> subRegionsSampled;
				for (const auto & strainSampled : genomesSampled) {
					Strain strain(strainSampled.seqBase_.name_,
							strainsExpected.at(strainSampled.seqBase_.name_).genomicRegionToHapNames_,
							strainSampled.genomeCnt_);
					strain.meta_ = strainsExpected.at(strainSampled.seqBase_.name_).meta_;
					hapAbundGenomesSampled_.emplace(strainSampled.seqBase_.name_, strain);
					for(const auto & hap : strain.genomicRegionToHapNames_){
						subRegionsSampled[hap.first][hap.second] += strain.relativeAbundance_;
					}
				}
				for(const auto & subRegionSampled : subRegionsSampled){
					hapRegionAmplified_.emplace(subRegionSampled.first, GenomicRegionAmp(subRegionSampled.first,subRegionSampled.second));
				}
			}

			Json::Value toJson() const {
				Json::Value ret;
				ret["class"] = njh::json::toJson(njh::getTypeName(*this));
				ret["runName_"] = njh::json::toJson(runName_);
				ret["hapAbundGenomesSampled_"] = njh::json::toJson(hapAbundGenomesSampled_);
				ret["hapRegionAmplified_"] = njh::json::toJson(hapRegionAmplified_);
				ret["expAmounts_"] = njh::json::toJson(expAmounts_);
				ret["meta_"] = njh::json::toJson(meta_);
				return ret;
			}

			static VecStr genJsonMem(){
				return VecStr{"runName_", "hapAbundGenomesSampled_", "expAmounts_", "hapRegionAmplified_", "meta_"};
			}

			static ExperimentRun genFromJson(const Json::Value & val){
				njh::json::MemberChecker memCheck(val);
				memCheck.failMemberCheckThrow(genJsonMem(), __PRETTY_FUNCTION__);
				ExperimentRun ret(val["runName_"].asString(), SeqSampledAmounts::genFromJson(val["expAmounts_"]));
				for(const auto & sName : val["hapAbundGenomesSampled_"].getMemberNames()){
					ret.hapAbundGenomesSampled_.emplace(sName, Strain::genFromJson(val["hapAbundGenomesSampled_"][sName]));
				}
				for(const auto & region : val["hapRegionAmplified_"].getMemberNames()){
					ret.hapRegionAmplified_.emplace(region, GenomicRegionAmp::genFromJson(val["hapRegionAmplified_"][region]));
				}
				ret.meta_ = MetaDataInName::genMetaFromJson(val["meta_"]);
				return ret;
			}

			static std::map<std::string, Strain> averageStrainsSampled(const std::vector<ExperimentRun> & runs, const uint32_t runsRequired){
				std::map<std::string, Strain> ret;

				std::map<std::string, std::vector<Strain>> acrossRuns;
				for(const auto & run : runs){
					for(const auto & strain : run.hapAbundGenomesSampled_){
						acrossRuns[strain.first].emplace_back(strain.second);
					}
				}
				for (const auto & strain : acrossRuns) {
					if (runsRequired == strain.second.size()) {
						double sum = 0;
						for (const auto & s : strain.second) {
							sum += s.relativeAbundance_;
						}
						auto nonZeroStrain = strain.second.front();
						nonZeroStrain.relativeAbundance_ = sum/runs.size();
						ret.emplace(nonZeroStrain.name_, nonZeroStrain);
					}else{
						auto zeroStrain = strain.second.front();
						zeroStrain.relativeAbundance_ = 0;
						ret.emplace(zeroStrain.name_, zeroStrain);
					}
				}
				//normalize fractions
				double total = 0;
				for(const auto & s : ret){
					total += s.second.relativeAbundance_;
				}
				for(auto & s : ret){
					s.second.relativeAbundance_ = s.second.relativeAbundance_/total;
				}
				return ret;
			}

			static std::map<std::string,std::map<std::string, double>> averageHapAbundSampled(const std::vector<ExperimentRun> & runs, const uint32_t runsRequired){
				std::map<std::string,std::map<std::string, double>> ret;
				std::map<std::string,std::map<std::string, std::vector<double>>> acrossRuns;
				for(const auto & run : runs){
					for(const auto & region : run.hapRegionAmplified_){
						for(const auto & hap : region.second.hapAbundSampled_){
							if(hap.second > 0){
								acrossRuns[region.first][hap.first].emplace_back(hap.second);
							}
						}
					}
				}
				for(const auto & region : acrossRuns){
					for(const auto & hap : region.second){
						if(hap.second.size() == runsRequired){
							auto sum = vectorSum(hap.second);
							ret[region.first][hap.first] = sum/runs.size();
						}else{
							ret[region.first][hap.first] = 0;
						}
					}
				}
				//normalize fractions
				for (auto & region : ret) {
					double total = 0;
					for (auto & hap : region.second) {
						total += hap.second;
					}
					for (auto & hap : region.second) {
						hap.second = hap.second / total;
					}
				}
				return ret;
			}

			static std::map<std::string,std::map<std::string, double>> averageHapAbundSequenced(const std::vector<ExperimentRun> & runs, const uint32_t runsRequired){
				std::map<std::string,std::map<std::string, double>> ret;
				std::map<std::string,std::map<std::string, std::vector<double>>> acrossRuns;
				for(const auto & run : runs){
					for(const auto & region : run.hapRegionAmplified_){
						for(const auto & hap : region.second.hapAbundSequenced_){
							if(hap.second > 0){
								acrossRuns[region.first][hap.first].emplace_back(hap.second);
							}
						}
					}
				}
				for(const auto & region : acrossRuns){
					for(const auto & hap : region.second){
						if(hap.second.size() == runsRequired){
							auto sum = vectorSum(hap.second);
							ret[region.first][hap.first] = sum/runs.size();
						}else{
							ret[region.first][hap.first] = 0;
						}
					}
				}
				//normalize fractions
				for (auto & region : ret) {
					double total = 0;
					for (auto & hap : region.second) {
						total += hap.second;
					}
					for (auto & hap : region.second) {
						hap.second = hap.second / total;
					}
				}
				return ret;
			}

			static std::map<std::string,std::map<std::string, double>> averageHapAbundObserved(const std::vector<ExperimentRun> & runs, const uint32_t runsRequired){
				std::map<std::string,std::map<std::string, double>> ret;
				std::map<std::string,std::map<std::string, std::vector<double>>> acrossRuns;
				for(const auto & run : runs){
					for(const auto & region : run.hapRegionAmplified_){
						for(const auto & hap : region.second.hapAbundObserved_){
							if(hap.second > 0){
								acrossRuns[region.first][hap.first].emplace_back(hap.second);
							}
						}
					}
				}
				for(const auto & region : acrossRuns){
					for(const auto & hap : region.second){
						if(hap.second.size() == runsRequired){
							auto sum = vectorSum(hap.second);
							ret[region.first][hap.first] = sum/runs.size();
						}else{
							ret[region.first][hap.first] = 0;
						}
					}
				}
				//normalize fractions
				for (auto & region : ret) {
					double total = 0;
					for (auto & hap : region.second) {
						total += hap.second;
					}
					for (auto & hap : region.second) {
						hap.second = hap.second / total;
					}
				}
				return ret;
			}

			static std::map<std::string, Strain> averageStrainsSampled(const std::vector<ExperimentRun> & runs){
				return averageStrainsSampled(runs, runs.size());
			}

			static std::map<std::string,std::map<std::string, double>> averageHapAbundSampled(const std::vector<ExperimentRun> & runs){
				return averageHapAbundSampled(runs, runs.size());
			}

			static std::map<std::string,std::map<std::string, double>> averageHapAbundSequenced(const std::vector<ExperimentRun> & runs){
				return averageHapAbundSequenced(runs, runs.size());
			}

			static std::map<std::string,std::map<std::string, double>> averageHapAbundObserved(const std::vector<ExperimentRun> & runs){
				return averageHapAbundObserved(runs, runs.size());
			}
		};

		std::string sampName_;
		std::map<std::string, std::vector<ExperimentRun>> expRuns_; /**< replicates will be placed in the same sub-vector*/
		std::map<std::string, Strain> strainsExpected_;
		MetaDataInName meta_;

		Json::Value toJson() const {
			Json::Value ret;
			ret["class"] = njh::json::toJson(njh::getTypeName(*this));
			ret["sampName_"] = njh::json::toJson(sampName_);
			ret["expRuns_"] = njh::json::toJson(expRuns_);
			ret["strainsExpected_"] = njh::json::toJson(strainsExpected_);
			ret["meta_"] = njh::json::toJson(meta_);
			return ret;
		}

		static VecStr genJsonMem() {
			return VecStr { "sampName_", "strainsExpected_", "expRuns_", "meta_" };
		}

		static Sample genFromJson(const Json::Value & val){
			njh::json::MemberChecker memCheck(val);
			memCheck.failMemberCheckThrow(genJsonMem(), __PRETTY_FUNCTION__);
			std::map<std::string, Strain> strainsExpected;
			for(const auto & strain : val["strainsExpected_"].getMemberNames()){
				strainsExpected.emplace(strain, Strain::genFromJson((val["strainsExpected_"][strain])));
			}
			Sample ret(val["sampName_"].asString(), strainsExpected);
			for(const auto & run : val["expRuns_"].getMemberNames()){
				std::vector<ExperimentRun> runs;
				for(const auto & subRun : val["expRuns_"][run] ){
					runs.emplace_back(ExperimentRun::genFromJson(subRun)) ;
				}
				ret.expRuns_.emplace(run, runs);
			}
			ret.meta_ = MetaDataInName::genMetaFromJson(val["meta_"]);
			return ret;
		}

		bool hasExpRun(const std::string & run){
			return njh::in(run, expRuns_);
		}

		void addExpRun(const std::string & run, const SeqSampledAmounts & amounts, const uint32_t numOfRuns = 1) {
			if (hasExpRun(run)) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "already have experiment "
						<< run << "\n";
				throw std::runtime_error { ss.str() };
			}

			std::vector<ExperimentRun> runs;
			for(uint32_t runNum = 1; runNum <= numOfRuns; ++runNum){
				std::string eRunName = run;
				if(numOfRuns > 1){
					eRunName += njh::pasteAsStr("-rep", runNum);
				}
				ExperimentRun eRun(eRunName, amounts);
				eRun.setGenomesSampled(strainsExpected_);
				runs.emplace_back(eRun);
			}
			expRuns_.emplace(run, runs);
		}
	};

	ControlPopulation(const std::string & name): populationName_(name){

	}
	std::string populationName_;
	std::map<std::string, Sample> samples_;
	MetaDataInName meta_;
	Json::Value toJson() const {
		Json::Value ret;
		ret["class"] = njh::json::toJson(njh::getTypeName(*this));
		ret["populationName_"] = njh::json::toJson(populationName_);
		ret["samples_"] = njh::json::toJson(samples_);
		ret["meta_"] = njh::json::toJson(meta_);
		return ret;
	}

	static VecStr genJsonMem() {
		return VecStr { "populationName_", "samples_", "meta_" };
	}

	static ControlPopulation genFromJson(const Json::Value & val){
		njh::json::MemberChecker memCheck(val);
		memCheck.failMemberCheckThrow(genJsonMem(), __PRETTY_FUNCTION__);
		ControlPopulation ret(val["populationName_"].asString());
		for(const auto & reg : val["samples_"].getMemberNames()){
			ret.samples_.emplace(reg, Sample::genFromJson((val["samples_"][reg])));
		}
		ret.meta_ = MetaDataInName::genMetaFromJson(val["meta_"]);
		return ret;
	}

	bool hasSample(const std::string & samp) {
		return njh::in(samp, samples_);
	}
	void addSample(const Sample & samp) {
		if (hasSample(samp.sampName_)) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "already have sample: " << samp.sampName_
					<< "\n";
			throw std::runtime_error { ss.str() };
		}
		samples_.emplace(samp.sampName_, samp);
	}
	void addSample(const std::string & samp,
			const std::map<std::string, Strain> & strainsExpected) {
		addSample(Sample(samp, strainsExpected));
	}
};



class MixtureSetUp {
public:

	struct PrimerSetup{
		PrimerSetup(const std::string & name,
				const std::string & forward,
				const std::string & reverse): name_(name),
						forward_(forward),
						forward_3_5_(seqUtil::reverseComplement(forward, "DNA")),
						reverse_(reverse),
						reverse_3_5_(seqUtil::reverseComplement(reverse, "DNA")){

		}
		std::string name_;
		std::string forward_;    //5`-3` direction
		std::string forward_3_5_;
		std::string reverse_; //5`-3` direction
		std::string reverse_3_5_;

		uint32_t forward_randomPrecedingBases_{0}; //this will insert up to this number of random bases before
		uint32_t reverse_randomPrecedingBases_{0}; //this will insert up to this number of random bases before
		bool randomOneLength_{false}; //random number of bases preceding is just one length



		Json::Value toJson() const {
			Json::Value ret;
			ret["name"] = njh::json::toJson(name_);
			ret["forward"] = njh::json::toJson(forward_);
			ret["reverse"] = njh::json::toJson(reverse_);
			ret["forward_randomPrecedingBases"] = njh::json::toJson(forward_randomPrecedingBases_);
			ret["reverse_randomPrecedingBases"] = njh::json::toJson(reverse_randomPrecedingBases_);
			ret["randomOneLength"] = njh::json::toJson(randomOneLength_);
			return ret;
		}
	};

	struct MidSetup{
		MidSetup(const std::string & name,
				const std::string & barcode): name_(name),
						barcode_(barcode),
						barcode_3_5_(seqUtil::reverseComplement(barcode, "DNA")){

		}
		std::string name_;
		std::string barcode_; //5`-3`
		std::string barcode_3_5_;
		uint32_t randomPrecedingBases_{0}; //this will insert up to this number of random bases before MID
		bool randomOneLength_{false}; //random number of bases preceding is just one length


		Json::Value toJson() const {
			Json::Value ret;
			ret["name"] = njh::json::toJson(name_);
			ret["barcode"] = njh::json::toJson(barcode_);
			ret["randomPrecedingBases"] = njh::json::toJson(randomPrecedingBases_);
			ret["randomOneLength"] = njh::json::toJson(randomOneLength_);
			return ret;
		}
	};

	MixtureSetUp(const std::string & name) :
			name_(name) {
	}

	std::string name_;

	std::shared_ptr<PrimerSetup> primers_;

	std::shared_ptr<MidSetup> forwardBarcode_;
	std::shared_ptr<MidSetup> reverseBarcode_;

	std::map<std::string, double> expectedAbundances_;

	std::map<std::string, uint64_t> genomeCounts_;

	uint32_t startingTemplateAmount_{3000};
	uint32_t finalReadAmount_{1000};

	std::unique_ptr<MetaDataInName> meta_;

	Json::Value toJson() const {
		Json::Value ret;
		ret["class"] = njh::json::toJson(njh::getTypeName(*this));
		if (nullptr != primers_) {
			ret["primers"] = njh::json::toJson(primers_);
		}
		if (nullptr != forwardBarcode_) {
			ret["forwardBarcode"] = njh::json::toJson(*forwardBarcode_);
		}
		if (nullptr != reverseBarcode_) {
			ret["reverseBarcode"] = njh::json::toJson(*reverseBarcode_);
		}
		ret["name"] = njh::json::toJson(name_);
		ret["expectedAbundances_"] = njh::json::toJson(expectedAbundances_);
		ret["genomeCounts_"] = njh::json::toJson(genomeCounts_);

		ret["startingTemplateAmount_"] = njh::json::toJson(startingTemplateAmount_);
		ret["finalReadAmount_"] = njh::json::toJson(finalReadAmount_);
		if(nullptr != meta_){
			ret["meta_"] = njh::json::toJson(*meta_);
		}
		return ret;
	}

	void addAbundance(const std::string & seqName, double abundnace){
		if(njh::in(seqName, expectedAbundances_)){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", already have abundance for " << seqName << " in mixutre: " << name_ << "\n";
			throw std::runtime_error{ss.str()};
		}
		if(abundnace <=0){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error in adding abundance for  " << seqName << " in mixutre: " << name_ << " abundance must be greater than 0, not " << abundnace << "\n";
			throw std::runtime_error{ss.str()};
		}
		expectedAbundances_[seqName] = abundnace;
	}

	void addGenomeCount(const std::string & seqName, uint32_t genomeCount){
		if(njh::in(seqName, genomeCounts_)){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", already have abundance for " << seqName << " in mixutre: " << name_ << "\n";
			throw std::runtime_error{ss.str()};
		}
		if(!njh::in(seqName, expectedAbundances_)){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", should have an expected abundance for " << seqName << " in mixutre: " << name_ << " when adding a genome count" << "\n";
			throw std::runtime_error{ss.str()};
		}
		genomeCounts_[seqName] = genomeCount;
	}


	void setPrimers(const std::string & primerName, const std::string & forward, const std::string & reverse){
		if(nullptr != primers_){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error, already set primers for mixture " <<  name_ << "\n";
			throw std::runtime_error{ss.str()};
		}
		primers_ = std::make_shared<PrimerSetup>(primerName, forward, reverse);
	}

	void setForwardBarcode(const std::string & barcodeName, const std::string & barcode){
		if(nullptr != forwardBarcode_){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error, already set forward barcode for mixture " << name_  << "\n";
			throw std::runtime_error{ss.str()};
		}
		if(nullptr != reverseBarcode_ && barcodeName != reverseBarcode_->name_){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error in setting barcode " << barcodeName << " for forward barcode, reverse barcode set but name doesn't match " << ", reverse name: " << reverseBarcode_->name_<< "\n";
			throw std::runtime_error{ss.str()};
		}
		if (barcode.empty() || barcodeName.empty()) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__
					<< ", error, barcode name or barcode cannot be empty for " << barcodeName << " in mixture: " << name_
					<< "\n";
			throw std::runtime_error { ss.str() };
		}
		forwardBarcode_ = std::make_shared<MidSetup>(barcodeName, barcode);
	}

	void setReverseBarcode(const std::string & barcodeName, const std::string & barcode){
		if(nullptr != reverseBarcode_){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error, already set reverse barcode for mixture " << name_  << "\n";
			throw std::runtime_error{ss.str()};
		}
		if(nullptr != forwardBarcode_ && barcodeName != forwardBarcode_->name_){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error in setting barcode " << barcodeName << " for reverse barcode, forward barcode set but name doesn't match " << ", forward name: " << forwardBarcode_->name_<< "\n";
			throw std::runtime_error{ss.str()};
		}
		if (barcode.empty() || barcodeName.empty()) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__
					<< ", error, barcode name or barcode cannot be empty for " << barcodeName << " in mixture: " << name_
					<< "\n";
			throw std::runtime_error { ss.str() };
		}
		reverseBarcode_ = std::make_shared<MidSetup>(barcodeName, barcode);
	}
};

class SampleSetup{

public:
	SampleSetup(const std::string & name): name_(name){

	}

	std::string name_;
	std::map<std::string, std::shared_ptr<MixtureSetUp>> mixtures_;

	void addMixture(const std::shared_ptr<MixtureSetUp> & mixture){
		if(hasMixture(mixture->name_)){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error, already have mixture " << mixture->name_;
			throw std::runtime_error{ss.str()};
		}
		mixtures_[mixture->name_] = mixture;
	}
	bool hasMixture(const std::string & name) const{
		return njh::in(name, mixtures_);
	}
	Json::Value toJson() const{
		Json::Value ret;
		ret["class"] = njh::json::toJson(njh::getTypeName(*this));
		ret["name"] = njh::json::toJson(name_);
		ret["mixtures"] = njh::json::toJson(mixtures_);
		return ret;
	}

};

class LibrarySetup {

public:

	struct SimLibrarySetupPars{

		SimLibrarySetupPars(){}
		explicit SimLibrarySetupPars(const Json::Value & val){
			njh::json::MemberChecker checker(val);
			checker.failMemberCheckThrow(VecStr{
				"pairedEndLength",
				"addReverseComplement",
				"addBluntEndingArtifact",
				"bluntEndingArtifactChance",
				"barcodeRandomPrecedingBases",
				"primerRandomPrecedingBases",
			  "noAddPrimers"}, __PRETTY_FUNCTION__);
			pairedEndLength_ = val["pairedEndLength"].asUInt();
			addReverseComplement_ = val["addReverseComplement"].asBool();
			addBluntEndingArtifact_ = val["addBluntEndingArtifact"].asBool();
			bluntEndingArtifactChance_ = val["bluntEndingArtifactChance"].asDouble();
			barcodeRandomPrecedingBases_ = val["barcodeRandomPrecedingBases"].asUInt();
			primerRandomPrecedingBases_ = val["primerRandomPrecedingBases"].asUInt();

			noAddPrimers_ = val["noAddPrimers"].asBool();
		}

		uint32_t pairedEndLength_ = 250;
		bool addReverseComplement_ = false;
		bool addBluntEndingArtifact_ = false;
		double bluntEndingArtifactChance_ = .10;
		uint32_t barcodeRandomPrecedingBases_ = 0;
		uint32_t primerRandomPrecedingBases_ = 0;
		bool noAddPrimers_ = false;

		Json::Value toJson() const{
			Json::Value ret;
			ret["class"] = njh::json::toJson(njh::getTypeName(*this));
			ret["pairedEndLength"] = njh::json::toJson(pairedEndLength_);
			ret["addReverseComplement"] = njh::json::toJson(addReverseComplement_);
			ret["addBluntEndingArtifact"] = njh::json::toJson(addBluntEndingArtifact_);
			ret["bluntEndingArtifactChance"] = njh::json::toJson(bluntEndingArtifactChance_);
			ret["barcodeRandomPrecedingBases"] = njh::json::toJson(barcodeRandomPrecedingBases_);
			ret["primerRandomPrecedingBases"] = njh::json::toJson(primerRandomPrecedingBases_);


			ret["noAddPrimers"] = njh::json::toJson(noAddPrimers_);

			return ret;
		}
	};

	LibrarySetup(const std::string & name, const SimLibrarySetupPars & pars) :
			name_(name), pars_(pars), pop_(name) {

	}

	LibrarySetup(const std::string & name, const SimLibrarySetupPars & pars,
			const Json::Value & librarySetUp,
			const VecStr & refSeqsNames) :
			name_(name), pars_(pars), pop_(name) {

		njh::json::MemberChecker memcheck(librarySetUp);
		memcheck.failMemberCheckThrow(jsonMembers(), __PRETTY_FUNCTION__);

		ids_ = std::make_unique<PrimersAndMids>(std::unordered_map<std::string, PrimersAndMids::Target>{});
		for(const auto & sample : librarySetUp["samples"]){
			njh::json::MemberChecker sampleChecker(sample);
			sampleChecker.failMemberCheckThrow({"name", "mixtures"}, __PRETTY_FUNCTION__);
			auto sampleSet = std::make_shared<SampleSetup>(sample["name"].asString());
			for(const auto & mixture : sample["mixtures"]){
				njh::json::MemberChecker mixtureChecker(mixture);
				mixtureChecker.failMemberCheckThrow({"name", "expectedAbundances_"}, __PRETTY_FUNCTION__);
				mixtureChecker.failMemberCheckThrow({"name", "genomeCounts_"}, __PRETTY_FUNCTION__);

				auto mixtureSet = std::make_shared<MixtureSetUp>(mixture["name"].asString());
				{
					auto memebers = mixture["expectedAbundances_"].getMemberNames();
					for(const auto & member : memebers){
						if(!njh::in(member, refSeqsNames)){
							std::stringstream ss;
							ss << __PRETTY_FUNCTION__ << ", error, " << " no reference seq named: " << member << "\n";
							ss << "options are: " << njh::conToStr(refSeqsNames) << "\n";
							throw std::runtime_error{ss.str()};
						}
						mixtureSet->addAbundance(member, mixture["expectedAbundances_"][member].asDouble());
					}
				}
				{
					auto memebers = mixture["genomeCounts_"].getMemberNames();
					for(const auto & member : memebers){
						if(!njh::in(member, refSeqsNames)){
							std::stringstream ss;
							ss << __PRETTY_FUNCTION__ << ", error, " << " no reference seq named: " << member << "\n";
							ss << "options are: " << njh::conToStr(refSeqsNames) << "\n";
							throw std::runtime_error{ss.str()};
						}
						mixtureSet->addGenomeCount(member, mixture["genomeCounts_"][member].asUInt64());
					}
				}
				if(mixture.isMember("meta_")){
					std::unordered_map<std::string, std::string> mixMeta;
					auto keys = mixture["meta_"]["meta_"].getMemberNames();
					for(const auto & k : keys){
						mixMeta[k] = mixture["meta_"]["meta_"][k].asString();
					}
					mixtureSet->meta_ = std::make_unique<MetaDataInName>();
					mixtureSet->meta_->meta_ = mixMeta;
				}
				if(mixture.isMember("primers")){
					if(ids_->hasTarget(mixture["primers"]["name"].asString())){
						if (ids_->targets_.at(mixture["primers"]["name"].asString()).info_.forwardPrimer_
								!= mixture["primers"]["forward"].asString()) {
							std::stringstream ss;
							ss << __PRETTY_FUNCTION__ << ", error, target "
									<< mixture["primers"]["name"].asString()
									<< " already added with a different forward primer " << "\n";
							ss << "New Forward: " << mixture["primers"]["forward"].asString()
									<< " , Current Forward: "
									<< ids_->targets_.at(mixture["primers"]["name"].asString()).info_.forwardPrimer_
									<< "\n";
							throw std::runtime_error { ss.str() };
						}
						if (ids_->targets_.at(mixture["primers"]["name"].asString()).info_.reversePrimer_
								!= mixture["primers"]["reverse"].asString()) {
							std::stringstream ss;
							ss << __PRETTY_FUNCTION__ << ", error, target "
									<< mixture["primers"]["name"].asString()
									<< " already added with a different forward primer " << "\n";
							ss << "New Reverse: " << mixture["primers"]["reverse"].asString()
									<< " , Current Reverse: "
									<< ids_->targets_.at(mixture["primers"]["name"].asString()).info_.reversePrimer_
									<< "\n";
							throw std::runtime_error { ss.str() };
						}
					}else{
						ids_->addTarget(
							mixture["primers"]["name"].asString(),
							mixture["primers"]["forward"].asString(),
							mixture["primers"]["reverse"].asString());
					}
					mixtureSet->setPrimers(
							mixture["primers"]["name"].asString(),
							mixture["primers"]["forward"].asString(),
							mixture["primers"]["reverse"].asString()
							);

					if(mixture["primers"].isMember("forward_randomPrecedingBases")){
						mixtureSet->primers_->forward_randomPrecedingBases_ = mixture["primers"]["forward_randomPrecedingBases"].asUInt();
					}
					if(mixture["primers"].isMember("reverse_randomPrecedingBases")){
						mixtureSet->primers_->reverse_randomPrecedingBases_ = mixture["primers"]["reverse_randomPrecedingBases"].asUInt();
					}
				}
				mixtureSet->startingTemplateAmount_ = mixture["startingTemplateAmount_"].asUInt();
				mixtureSet->finalReadAmount_ = mixture["finalReadAmount_"].asUInt();

				std::string barcodeName = "";
				std::string forBar = "";
				std::string revBar = "";
				if(mixture.isMember("forwardBarcode")){
					mixtureSet->setForwardBarcode(
							mixture["forwardBarcode"]["name"].asString(),
							mixture["forwardBarcode"]["barcode"].asString()
							);
					barcodeName = mixture["forwardBarcode"]["name"].asString();
					forBar = mixture["forwardBarcode"]["barcode"].asString();
					if(mixture["forwardBarcode"].isMember("randomPrecedingBases")){
						mixtureSet->forwardBarcode_->randomPrecedingBases_ = mixture["forwardBarcode"]["randomPrecedingBases"].asUInt();
					}
				}
				if(mixture.isMember("reverseBarcode")){
					mixtureSet->setReverseBarcode(
							mixture["reverseBarcode"]["name"].asString(),
							mixture["reverseBarcode"]["barcode"].asString()
							);
					revBar = mixture["reverseBarcode"]["barcode"].asString();
					if(mixture["reverseBarcode"].isMember("randomPrecedingBases")){
						mixtureSet->reverseBarcode_->randomPrecedingBases_ = mixture["reverseBarcode"]["randomPrecedingBases"].asUInt();
					}
				}
				if("" != barcodeName){
					if(njh::in(barcodeName, ids_->mids_)){
						if(ids_->mids_.at(barcodeName).forwardBar_->bar_->motifOriginal_ != forBar){
							std::stringstream ss;
							ss << __PRETTY_FUNCTION__ << ", error, MID "
									<< barcodeName
									<< " already added with a different forward barcode " << "\n";
							ss << "New Barcode: " << forBar
									<< " , Current Reverse: "
									<< ids_->mids_.at(barcodeName).forwardBar_->bar_->motifOriginal_
									<< "\n";
							throw std::runtime_error { ss.str() };
						}
						if("" != revBar && nullptr == ids_->mids_.at(barcodeName).reverseBar_ ){
							std::stringstream ss;
							ss << __PRETTY_FUNCTION__ << ", error, MID " << barcodeName
									<< " already added without a reverwe barcode but new barcode has one, "
									<< revBar << "\n";
							throw std::runtime_error { ss.str() };
						}else if ("" == revBar && nullptr != ids_->mids_.at(barcodeName).reverseBar_ ){
							std::stringstream ss;
							ss << __PRETTY_FUNCTION__ << ", error, MID " << barcodeName
									<< " already added has a reverwe barcode but new barcode has doesnt have one, have "
									<< ids_->mids_.at(barcodeName).reverseBar_->bar_->motifOriginal_ << "\n";
							throw std::runtime_error { ss.str() };
						}else if("" != revBar && nullptr != ids_->mids_.at(barcodeName).reverseBar_ ){
							if(ids_->mids_.at(barcodeName).reverseBar_->bar_->motifOriginal_ != revBar){
								std::stringstream ss;
								ss << __PRETTY_FUNCTION__ << ", error, MID "
										<< barcodeName
										<< " already added with a different forward barcode " << "\n";
								ss << "New Barcode: " << revBar
										<< " , Current Reverse: "
										<< ids_->mids_.at(barcodeName).reverseBar_->bar_->motifOriginal_
										<< "\n";
								throw std::runtime_error { ss.str() };
							}
						}
					}else{
						if("" != revBar){
							ids_->addMid(barcodeName, forBar, revBar);
						}else{
							ids_->addMid(barcodeName, forBar);
						}
					}
				}
				sampleSet->addMixture(mixtureSet);
			}
			addSample(sampleSet);
		}
		pop_ = ControlPopulation::genFromJson(librarySetUp["pop_"]);
	}
	std::string name_;
	std::unordered_map<std::string, std::shared_ptr<SampleSetup>> samples_;

	SimLibrarySetupPars pars_;

	std::unique_ptr<PrimersAndMids> ids_;
	ControlPopulation pop_;

	void addSample(const std::shared_ptr<SampleSetup> & sample){
		if(hasSample(sample->name_)){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error, already have sample " << sample->name_;
			throw std::runtime_error{ss.str()};
		}
		samples_[sample->name_] = sample;
	}
	bool hasSample(const std::string & name) const{
		return njh::in(name, samples_);
	}

	static VecStr jsonMembers(){
		return {
		"name",
		"samples",
		"pars",
		"pop_"};
	}
	Json::Value toJson() const{
		Json::Value ret;
		ret["class"] = njh::json::toJson(njh::getTypeName(*this));
		ret["name"] = njh::json::toJson(name_);
		ret["samples"] = njh::json::toJson(samples_);
		ret["pars"] = njh::json::toJson(pars_);
		ret["pop_"] = njh::json::toJson(pop_);
		return ret;
	}


	table createAbundanceTable() const{
		table ret(VecStr{"library", "sample", "mix", "hap", "abundanceRaw", "totalRawAbundance", "Percentage"});
		for(const auto & sample : samples_){
			for(const auto & mix : sample.second->mixtures_){
				double totalRawAbundance = 0.0;
				for(const auto & hap : mix.second->expectedAbundances_){
					totalRawAbundance += hap.second;
				}
				for(const auto & hap : mix.second->expectedAbundances_){
					ret.addRow(name_, sample.second->name_, mix.second->name_, hap.first, hap.second, totalRawAbundance, 100 * (hap.second/totalRawAbundance));
				}
			}
		}
		return ret;
	}

};




}  // namespace njhseq




