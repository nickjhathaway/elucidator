#pragma once

/*
 * LibrarySetup.hpp
 *
 *  Created on: May 2, 2019
 *      Author: nicholashathaway
 */

#include <SeekDeep/objects/PrimersAndMids.hpp>

#include "elucidator/common.h"

namespace njhseq {

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

	std::unordered_map<std::string, double> abundances_;

	uint32_t startingTemplateAmount_{3000};
	uint32_t finalReadAmount_{1000};


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
		ret["abundances"] = njh::json::toJson(abundances_);
		ret["startingTemplateAmount_"] = njh::json::toJson(startingTemplateAmount_);
		ret["finalReadAmount_"] = njh::json::toJson(finalReadAmount_);

		return ret;
	}

	void addAbundance(const std::string & seqName, double abundnace){
		if(njh::in(seqName, abundances_)){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", already have abundance for " << seqName << " in mixutre: " << name_ << "\n";
			throw std::runtime_error{ss.str()};
		}
		if(abundnace <=0){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error in adding abundance for  " << seqName << " in mixutre: " << name_ << " abundance must be greater than 0, not " << abundnace << "\n";
			throw std::runtime_error{ss.str()};
		}
		abundances_[seqName] = abundnace;
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
	std::unordered_map<std::string, std::shared_ptr<MixtureSetUp>> mixtures_;

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
				"sampleReadAmount",
				"pairedEndLength",
				"addReverseComplement",
				"addBluntEndingArtifact",
				"bluntEndingArtifactChance",
				"barcodeRandomPrecedingBases",
				"primerRandomPrecedingBases",
			  "noAddPrimers"}, __PRETTY_FUNCTION__);
			sampleReadAmount_ = val["sampleReadAmount"].asUInt();
			pairedEndLength_ = val["pairedEndLength"].asUInt();
			addReverseComplement_ = val["addReverseComplement"].asBool();
			addBluntEndingArtifact_ = val["addBluntEndingArtifact"].asBool();
			bluntEndingArtifactChance_ = val["bluntEndingArtifactChance"].asDouble();
			barcodeRandomPrecedingBases_ = val["barcodeRandomPrecedingBases"].asUInt();
			primerRandomPrecedingBases_ = val["primerRandomPrecedingBases"].asUInt();

			noAddPrimers_ = val["noAddPrimers"].asBool();
		}
		uint32_t sampleReadAmount_ = 5000;
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
			ret["sampleReadAmount"] = njh::json::toJson(sampleReadAmount_);
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
			name_(name), pars_(pars) {

	}

	LibrarySetup(const std::string & name,
			const SimLibrarySetupPars & pars,
			const Json::Value & librarySetUp,
			const std::unordered_map<std::string, seqInfo> & refSeqs
			) :
			name_(name), pars_(pars) {

		ids_ = std::make_unique<PrimersAndMids>(std::unordered_map<std::string, PrimersAndMids::Target>{});
		for(const auto & sample : librarySetUp["samples"]){
			njh::json::MemberChecker sampleChecker(sample);
			sampleChecker.failMemberCheckThrow({"name", "mixtures"}, __PRETTY_FUNCTION__);
			auto sampleSet = std::make_shared<SampleSetup>(sample["name"].asString());
			for(const auto & mixture : sample["mixtures"]){
				njh::json::MemberChecker mixtureChecker(mixture);
				mixtureChecker.failMemberCheckThrow({"name", "abundances"}, __PRETTY_FUNCTION__);
				auto mixtureSet = std::make_shared<MixtureSetUp>(mixture["name"].asString());
				auto memebers = mixture["abundances"].getMemberNames();
				for(const auto & member : memebers){
					if(!njh::in(member, refSeqs)){
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ", error, " << " no reference seq named: " << member << "\n";
						ss << "options are: " << njh::conToStr(getVectorOfMapKeys(refSeqs)) << "\n";
						throw std::runtime_error{ss.str()};
					}
					mixtureSet->addAbundance(member, mixture["abundances"][member].asDouble());
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
				if(mixture.isMember("startingTemplateAmount_")){
					mixtureSet->startingTemplateAmount_ = mixture["startingTemplateAmount_"].asUInt();
				}
				if(mixture.isMember("finalReadAmount_")){
					mixtureSet->finalReadAmount_ = mixture["finalReadAmount_"].asUInt();
				}else{
					mixtureSet->finalReadAmount_ = pars_.sampleReadAmount_;
				}
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

	}
	std::string name_;
	std::unordered_map<std::string, std::shared_ptr<SampleSetup>> samples_;

	SimLibrarySetupPars pars_;

	std::unique_ptr<PrimersAndMids> ids_;

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
		"pars"};
	}
	Json::Value toJson() const{
		Json::Value ret;
		ret["class"] = njh::json::toJson(njh::getTypeName(*this));
		ret["name"] = njh::json::toJson(name_);
		ret["samples"] = njh::json::toJson(samples_);
		ret["pars"] = njh::json::toJson(pars_);
		return ret;
	}

	table createAbundanceTable() const{
		table ret(VecStr{"library", "sample", "mix", "hap", "abundanceRaw", "totalRawAbundance", "Percentage"});
		for(const auto & sample : samples_){
			for(const auto & mix : sample.second->mixtures_){
				double totalRawAbundance = 0.0;
				for(const auto & hap : mix.second->abundances_){
					totalRawAbundance += hap.second;
				}
				for(const auto & hap : mix.second->abundances_){
					ret.addRow(name_, sample.second->name_, mix.second->name_, hap.first, hap.second, totalRawAbundance, 100 * (hap.second/totalRawAbundance));
				}
			}
		}
		return ret;
	}

};




}  // namespace njhseq




