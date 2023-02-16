/*
 * createLibrarySimMultipleMixtureDrugResistant.cpp
 *
 *  Created on: Oct 7, 2019
 *      Author: nicholashathaway
 */




#include "readSimulatorRunner.hpp"
#include "elucidator/simulation.h"
#include <SeekDeep/objects/TarAmpSetupUtils/PrimersAndMids.hpp>

namespace njhseq {



int readSimulatorRunner::createLibrarySimMultipleMixtureDrugResistant(
		const njh::progutils::CmdArgs & inputCommands) {
	LibrarySetup::SimLibrarySetupPars simPars;
	bfs::path coiTableFnp = "";
	bfs::path haplotypeInfo = "";
	bfs::path primerMidFnp = "";
	bfs::path patientSetupFile = "";
	std::string libraryName = "";
	bool simReplicates = false;
	bool rawPatientSetupFile = false;
	PCRAmountPars pcrNumbers;
	std::vector<uint32_t> timePoints;
	double adjustAroundGivenFrac = 0.5;
	double twoSdFrac = std::numeric_limits<double>::max();

	readSimulatorSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processVerbose();
	setUp.setOption(pcrNumbers.startingTemplateAmounts_, "--startingTemplateAmount", "Starting PCR Template Amount", njh::progutils::ProgramSetUp::ConCheckCase::NONZERO);
	setUp.setOption(pcrNumbers.finalReadAmount_, "--perMixtureReadAmount", "Final Read Amount to create per mixture", njh::progutils::ProgramSetUp::ConCheckCase::NONZERO);
	setUp.setOption(simReplicates, "--replicates", "Replicates");
	setUp.setOption(adjustAroundGivenFrac, "--adjustAroundGivenFrac", "Adjust Around Given Frac", njh::progutils::ProgramSetUp::CheckCase::GREATERZERO);
	adjustAroundGivenFrac = std::min(adjustAroundGivenFrac, 0.99);
	setUp.setOption(twoSdFrac, "--twoSdFrac", "Two Sd Fracs around genome starting template amount", njh::progutils::ProgramSetUp::CheckCase::GREATERZERO);
	setUp.setOption(timePoints, "--timePoints", "Additional time points, will automatically have time point 00", njh::progutils::ProgramSetUp::ConCheckCase::NONZERO);

	setUp.setOption(rawPatientSetupFile, "--rawPatientSetupFile", "Patient Setup File designates individual samples");

	setUp.setOption(patientSetupFile, "--patientSetupFile", "Patient Setup File", true);
	setUp.setOption(coiTableFnp, "--coiTable", "COI Table", true);
	setUp.setOption(haplotypeInfo, "--haplotypeInfo", "Haplotype Info", true);
	setUp.setOption(primerMidFnp, "--primerMidFnp", "Primer MID Fnp", true);
	setUp.setOption(libraryName, "--libraryName", "Library Name", true);

	setUp.setOption(simPars.pairedEndLength_, "--pairedEndLength", "Paired End Length");
	setUp.setOption(simPars.noAddPrimers_, "--noAddPrimers", "Primers are already present");
	setUp.setOption(simPars.barcodeRandomPrecedingBases_, "--barcodeRandomPrecedingBases", "Barcode Random Preceding Bases");
	setUp.setOption(simPars.primerRandomPrecedingBases_, "--primerRandomPrecedingBases", "Primer Random Preceding Bases");
	setUp.setOption(simPars.addReverseComplement_, "--addReverseComplement", "Add Reverse Complement");
	setUp.setOption(simPars.addBluntEndingArtifact_, "--addBluntEndingArtifact", "Add Blunt Ending Artifact");
	setUp.setOption(simPars.bluntEndingArtifactChance_, "--bluntEndingArtifactChance", "Blunt Ending Artifact Chance");

	setUp.processDirectoryOutputName(libraryName, true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_ );
	njh::sort(timePoints);
	uint32_t maxTime = 0;
	if(!timePoints.empty()){
		maxTime = *std::max_element(timePoints.begin(), timePoints.end());
	}
	LibrarySetup lSetup(libraryName, simPars);

	table coiTable;;
	table haplotypeInfoTable(haplotypeInfo, "\t", true);
	table patientSetupTable(patientSetupFile, "\t", true);
	lSetup.ids_ = std::make_unique<PrimersAndMids>(primerMidFnp);
	lSetup.ids_->initPrimerDeterminator();
	if(lSetup.ids_->containsMids()){
		lSetup.ids_->initMidDeterminator(MidDeterminator::MidDeterminePars{});
	}
	bool patientSetUpHasCOICol = njh::in<std::string>("COI", patientSetupTable.columnNames_);

	if (rawPatientSetupFile) {

	} else {
		//checks
		patientSetupTable.checkForColumnsThrow(VecStr{"number"}, __PRETTY_FUNCTION__);
		if(patientSetupTable.nCol() < 2){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ <<  " error, there should be more than 1 column" << "\n";
			throw std::runtime_error{ss.str()};
		}
	}
	if(!patientSetUpHasCOICol){
		coiTable = table(coiTableFnp, "\t", true);
	}

	//check for header in other tables
	VecStr metaLevels = patientSetupTable.columnNames_;
	if(rawPatientSetupFile){
		if(njh::in(std::string("name"), metaLevels ) ){
			removeElement<std::string>(metaLevels, "name");
		}
		if(njh::in(std::string("replicate"), metaLevels ) ){
			removeElement<std::string>(metaLevels, "replicate");
		}
		if(njh::in(std::string("COI"), metaLevels ) ){
			removeElement<std::string>(metaLevels, "COI");
		}
    if(njh::in(std::string("startingTemplate"), metaLevels ) ){
      removeElement<std::string>(metaLevels, "startingTemplate");
    }
	}else{
		removeElement<std::string>(metaLevels, "number");
	}
	haplotypeInfoTable.checkForColumnsThrow(njh::catVecs(metaLevels, VecStr{"frac"}), __PRETTY_FUNCTION__);
	if(!patientSetUpHasCOICol){
		coiTable.checkForColumnsThrow(njh::catVecs(metaLevels, VecStr{"COI", "n"}), __PRETTY_FUNCTION__);
	}
	//checks for uniqueness
	auto checkForUniqLevels = [](const table & tab, const VecStr & header, const bfs::path & fnp){
		std::unordered_set<std::string> levels;
		for (const auto & row : tab) {
			VecStr currentLevelVec;
			for (const auto & col : header) {
				currentLevelVec.emplace_back(row[tab.getColPos(col)]);
			}
			std::string currentLevel = njh::pasteAsStr(currentLevelVec);
			if (njh::in(currentLevel, levels)) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << " error, there should be unique levels in "
						<< fnp << "\n";
				ss << "Found " << njh::conToStr(currentLevelVec, ",") << " already "
						<< "\n";
				throw std::runtime_error { ss.str() };
			} else {
				levels.emplace(currentLevel);
			}
		}
	};
	auto getLevelsForHeader = [](const table & tab, VecStr header){
		std::unordered_set<std::string> levels;
		for (const auto & row : tab) {
			VecStr currentLevelVec;
			for (const auto & col : header) {
				currentLevelVec.emplace_back(row[tab.getColPos(col)]);
			}
			std::string currentLevel = njh::pasteAsStr(currentLevelVec);
			levels.emplace(currentLevel);
		}
		return levels;
	};
	auto checkForLevelsForHeader = [](const table & tab, const std::unordered_set<std::string> & levels, VecStr header, const bfs::path & fnp){
		std::unordered_set<std::string> currentLevels;
		for (const auto & row : tab) {
			VecStr currentLevelVec;
			for (const auto & col : header) {
				currentLevelVec.emplace_back(row[tab.getColPos(col)]);
			}
			std::string currentLevel = njh::pasteAsStr(currentLevelVec);
			currentLevels.emplace(currentLevel);
		}
		std::set<std::string> missing;
		for(const auto & lev :levels){
			if(!njh::in(lev, currentLevels)){
				missing.emplace(lev);
			}
		}
		if(!missing.empty()){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << " error, missing levels "
					<<  njh::conToStr(missing, ",") << " from " << fnp << "\n";
			throw std::runtime_error { ss.str() };
		}
	};
	if(!rawPatientSetupFile){
		checkForUniqLevels(patientSetupTable, metaLevels, patientSetupFile);
	}
	if(!patientSetUpHasCOICol){
		checkForUniqLevels(coiTable, njh::catVecs(metaLevels, VecStr{"COI"}), coiTableFnp);
	}
	auto levels = getLevelsForHeader(patientSetupTable, metaLevels);
	checkForLevelsForHeader(haplotypeInfoTable, levels, metaLevels, haplotypeInfo);
	if(!patientSetUpHasCOICol){
		checkForLevelsForHeader(coiTable, levels, metaLevels, coiTableFnp);
	}
	//primer check
	VecStr ignoreColumns = metaLevels;
	addOtherVec(ignoreColumns, VecStr{"n","total","frac","resistant","resistantFactor"});
	VecStr primerPairNames;
	for(const auto & col : haplotypeInfoTable.columnNames_){
		if(!njh::in(col, ignoreColumns)){
			primerPairNames.emplace_back(col);
		}
	}
	if(primerPairNames.empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << " error, couldn't determine primer pair names in "
				<<  haplotypeInfo << "\n";
		throw std::runtime_error { ss.str() };
	}
	VecStr missingPrimerPairs;
	for(const auto & primerPair : primerPairNames){
		if(!njh::in(primerPair, lSetup.ids_->pDeterminator_->primers_)){
			missingPrimerPairs.emplace_back(primerPair);
		}
	}
	if (!missingPrimerPairs.empty()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << " error, missing primer pair names in "
				<< primerMidFnp << "\n";
		ss << "Missing: " << njh::conToStr(missingPrimerPairs) << "\n";
		throw std::runtime_error { ss.str() };
	}

	checkForUniqLevels(haplotypeInfoTable, njh::catVecs(metaLevels, primerPairNames), haplotypeInfo);

	//create COI generator
	std::unordered_map<std::string, std::pair<std::vector<uint32_t>, std::vector<uint32_t>>> coiCounts;
	std::unordered_map<std::string, njh::randObjectGen<uint32_t, uint32_t>> coiRGenerator;
	if(!patientSetUpHasCOICol){
		for(const auto & row : coiTable){
			VecStr currentLevelVec;
			for (const auto & col : metaLevels) {
				currentLevelVec.emplace_back(row[coiTable.getColPos(col)]);
			}
			std::string currentLevel = njh::pasteAsStr(currentLevelVec);
			uint32_t coi = njh::StrToNumConverter::stoToNum<uint32_t>(row[coiTable.getColPos("COI")]);
			uint32_t n = njh::StrToNumConverter::stoToNum<uint32_t>(row[coiTable.getColPos("n")]);
			coiCounts[currentLevel].first.emplace_back(coi);
			coiCounts[currentLevel].second.emplace_back(n);
		}
		for(const auto & coi : coiCounts){
			coiRGenerator.emplace(coi.first, njh::randObjectGen<uint32_t, uint32_t>(coi.second.first, coi.second.second));
		}
	}
	//create haplotype chooser
	//   create a key for identifier
	//   uid, primer pair name, hap for primer pair
	std::unordered_map<std::string, std::map<std::string, std::string>> haplotypeKey;
	auto justHapNames_from_haplotypeInfoTable = haplotypeInfoTable.getColumns(primerPairNames).getUniqueRows();
	justHapNames_from_haplotypeInfoTable.setColNamePositions();

	for(const auto & row : justHapNames_from_haplotypeInfoTable){
		std::map<std::string, std::string> haps;
		VecStr hapNames;
		for(const auto & primerpair : primerPairNames){
			haps[primerpair] = row[justHapNames_from_haplotypeInfoTable.getColPos(primerpair)];
			hapNames.emplace_back(row[justHapNames_from_haplotypeInfoTable.getColPos(primerpair)]);
		}
		njh::sort(hapNames);
		haplotypeKey.emplace(njh::conToStr(hapNames, "__"), haps);
	}

	// resistance factor
	std::unordered_map<std::string, double> haplotypeKeyTorResistanceFactor;
	auto haplotypeInfoTableRes = haplotypeInfoTable.getColumns(njh::catVecs(primerPairNames, VecStr{"resistantFactor"})).getUniqueRows();
	checkForUniqLevels(haplotypeInfoTableRes, njh::catVecs(primerPairNames, VecStr{"resistantFactor"}), "resistance_" + haplotypeInfo.string());
	for(const auto & row : haplotypeInfoTableRes){
		VecStr hapNames;
		for(const auto & primerpair : primerPairNames){
			hapNames.emplace_back(row[haplotypeInfoTableRes.getColPos(primerpair)]);
		}
		njh::sort(hapNames);
		double res = njh::StrToNumConverter::stoToNum<double>(row[haplotypeInfoTableRes.getColPos("resistantFactor")]);
		haplotypeKeyTorResistanceFactor[njh::conToStr(hapNames, "__")] = res;
	}
	//  hap chooser gen
	std::unordered_map<std::string, std::unordered_map<std::string, double>> hapFracForLevels;
	std::unordered_map<std::string, std::pair<std::vector<std::string>, std::vector<double>>> hapCounts;
	std::unordered_map<std::string, njh::randObjectGen<std::string, double>> hapRGenerator;


	for(const auto & row : haplotypeInfoTable){
		VecStr currentLevelVec;
		for (const auto & col : metaLevels) {
			currentLevelVec.emplace_back(row[haplotypeInfoTable.getColPos(col)]);
		}
		std::string currentLevel = njh::pasteAsStr(currentLevelVec);
		VecStr hapNames;
		for(const auto & primerpair : primerPairNames){
			hapNames.emplace_back(row[haplotypeInfoTable.getColPos(primerpair)] );
		}
		njh::sort(hapNames);
		std::string hapUID = njh::conToStr(hapNames, "__");
		double frac = njh::StrToNumConverter::stoToNum<double>(row[haplotypeInfoTable.getColPos("frac")]);
		hapFracForLevels[currentLevel][hapUID] = frac;
		hapCounts[currentLevel].first.emplace_back(hapUID);
		hapCounts[currentLevel].second.emplace_back(frac);
	}
	for(const auto & hap : hapCounts){
		hapRGenerator.emplace(hap.first, njh::randObjectGen<std::string, double>(hap.second.first, hap.second.second));
	}
	struct TimePoint {
		uint32_t time_{std::numeric_limits<uint32_t>::max()};
    uint32_t startingTemplate_{std::numeric_limits<uint32_t>::max()};
    std::unordered_map<std::string, double> hapFracs_;

		Json::Value toJson() const{
			Json::Value ret;
			ret["class"] = njh::json::toJson(njh::getTypeName(*this));
			ret["time_"] = njh::json::toJson(time_);
			ret["hapFracs_"] = njh::json::toJson(hapFracs_);
			return ret;
		}
	};
	struct DrugResPatient {
		DrugResPatient(const std::string & name, const MetaDataInName & meta) :
				name_(name), meta_(meta) {
		}
		std::string name_;
		MetaDataInName meta_;
		std::vector<TimePoint> timePoints_;
		bool replicate_{false};
		uint32_t initalCOI_{1};
    uint32_t initialStartingTemplate_{std::numeric_limits<uint32_t>::max()};
		Json::Value toJson() const{
			Json::Value ret;
			ret["class"] = njh::json::toJson(njh::getTypeName(*this));
			ret["name_"] = njh::json::toJson(name_);
			ret["meta_"] = njh::json::toJson(meta_.meta_);
			ret["timePoints_"] = njh::json::toJson(timePoints_);
			ret["replicate_"] = njh::json::toJson(replicate_);
			ret["initalCOI_"] = njh::json::toJson(initalCOI_);
      ret["initialStartingTemplate_"] = njh::json::toJson(initialStartingTemplate_);
			return ret;
		}
	};
	//std::vector<DrugResPatient> patients;
	uint32_t totalPatients = 0;
	for(const auto & row : patientSetupTable){
		if(rawPatientSetupFile){
			++totalPatients;
		}else{
			totalPatients += njh::StrToNumConverter::stoToNum<uint32_t>(row[patientSetupTable.getColPos("number")]);
		}
	}
	njh::randomGenerator rGen;
	OutOptions coiOutOpts(njh::files::make_path(setUp.pars_.directoryName_, "coiCounts.tab.txt"));
	OutputStream coiCout(coiOutOpts);
	coiCout << "Patient\tinitialCOI" << std::endl;
	std::vector<DrugResPatient> patientSetUpPars;
	//
	if(rawPatientSetupFile){
		uint32_t patientCount = 1;
		bool hasNameCol = njh::in<std::string>("name", patientSetupTable.columnNames_);
		bool hasReplicateCol = njh::in<std::string>("replicate", patientSetupTable.columnNames_);
    bool hasStartingTemplateCol = njh::in<std::string>("startingTemplate", patientSetupTable.columnNames_);

    if(hasStartingTemplateCol){
      pcrNumbers.startingTemplateAmounts_ = {2048};
    }

		std::set<std::string> alreadyAddedNames;
		for(const auto & row : patientSetupTable){
			MetaDataInName meta;
			for (const auto & col : metaLevels) {
				meta.addMeta(col, row[patientSetupTable.getColPos(col)]);
			}
			DrugResPatient currentPatient(njh::pasteAsStr("P", leftPadNumStr(patientCount, totalPatients)), meta);
			if(hasNameCol){
				currentPatient.name_ = row[patientSetupTable.getColPos("name")];
				if(njh::in(currentPatient.name_, alreadyAddedNames)){
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", already added " << currentPatient.name_ << "\n";
					throw std::runtime_error{ss.str()};
				}
				alreadyAddedNames.emplace(currentPatient.name_);
			}
			auto currentLevel = currentPatient.meta_.pasteLevels(metaLevels);
			uint32_t initalCOI;
			if(patientSetUpHasCOICol){
				initalCOI = njh::StrToNumConverter::stoToNum<uint32_t>(row[patientSetupTable.getColPos("COI")]);
			}else{
				initalCOI = njh::mapAt(coiRGenerator, currentLevel).genObj();
			}
			coiCout << currentPatient.name_ << "\t" << initalCOI << std::endl;
			currentPatient.initalCOI_ = initalCOI;
			if (simReplicates && !hasReplicateCol) {
				currentPatient.replicate_ = true;
			} else if (hasReplicateCol) {
				currentPatient.replicate_ = "true"
						== njh::strToLowerRet(row[patientSetupTable.getColPos("replicate")])
						|| "t"
								== njh::strToLowerRet(
										row[patientSetupTable.getColPos("replicate")] ) ;
			}
      if(hasStartingTemplateCol){
        currentPatient.initialStartingTemplate_ = njh::StrToNumConverter::stoToNum<uint32_t>(row[patientSetupTable.getColPos("startingTemplate")]);
      }
			patientSetUpPars.emplace_back(currentPatient);
			++patientCount;
		}
	} else {
		uint32_t patientCount = 1;
		for(const auto & row : patientSetupTable){
			uint32_t number = njh::StrToNumConverter::stoToNum<uint32_t>(row[patientSetupTable.getColPos("number")]);
			for(uint32_t num = 0; num < number; ++num) {
				MetaDataInName meta;
				for (const auto & col : metaLevels) {
					meta.addMeta(col, row[patientSetupTable.getColPos(col)]);
				}
				DrugResPatient currentPatient(njh::pasteAsStr("P", leftPadNumStr(patientCount, totalPatients)), meta);
				auto currentLevel = currentPatient.meta_.pasteLevels(metaLevels);
				uint32_t initalCOI = njh::mapAt(coiRGenerator, currentLevel).genObj();
				initalCOI = std::min<uint32_t>(initalCOI, njh::mapAt(hapRGenerator,currentLevel).objs().size() );
				coiCout << currentPatient.name_ << "\t" << initalCOI << std::endl;
				currentPatient.initalCOI_ = initalCOI;
				if(simReplicates){
					currentPatient.replicate_ = true;
				}
				patientSetUpPars.emplace_back(currentPatient);
				++patientCount;
			}
		}
	}
	//add time points
	for(auto & currentPatient : patientSetUpPars){
		if(setUp.pars_.verbose_){
			std::cout << currentPatient.name_ << std::endl;
		}
		auto currentLevel = currentPatient.meta_.pasteLevels(metaLevels);
		std::vector<std::string> haplotypes{njh::mapAt(hapRGenerator, currentLevel).genObj()};
		if(currentPatient.initalCOI_ > 1){
			if(currentPatient.initalCOI_ >= njh::mapAt(hapRGenerator, currentLevel).objs().size()){
				if (setUp.pars_.verbose_) {
					std::cout << "For patient " << currentPatient.name_ << " COI of "
							<< currentPatient.initalCOI_
							<< " is greater than the number of haplotypes supplied, setting COI to "
							<< njh::mapAt(hapRGenerator, currentLevel).objs().size()
							<< std::endl;
					;
				}
				currentPatient.initalCOI_ = njh::mapAt(hapRGenerator, currentLevel).objs().size();
				haplotypes = njh::mapAt(hapRGenerator, currentLevel).objs();
			}else{
				for(uint32_t hapNum = 1; hapNum < currentPatient.initalCOI_; ++hapNum){
					if(setUp.pars_.verbose_){
						std::cout << "Simulating hapNum: " << hapNum << std::endl;
					}
					std::string currentHap = njh::mapAt(hapRGenerator, currentLevel).genObj();
					while(njh::in(currentHap, haplotypes)){
						currentHap = njh::mapAt(hapRGenerator,currentLevel).genObj();
					}
					if(setUp.pars_.verbose_){
						std::cout << "\tPicked: " << currentHap << std::endl;
					}
					haplotypes.emplace_back(currentHap);
				}
			}

			TimePoint initialTimePoint;
			initialTimePoint.time_ = 0;

			for(const auto & hap : haplotypes){
				double hapFrac = 100 * rGen.unifRand(1-adjustAroundGivenFrac, 1+adjustAroundGivenFrac) * hapFracForLevels[currentLevel][hap];
				initialTimePoint.hapFracs_[hap] = hapFrac;
			}
			currentPatient.timePoints_.emplace_back(initialTimePoint);

			for(const auto & time : timePoints){
				TimePoint tPoint;
				tPoint.time_ = time;
				double diffInTime = time - currentPatient.timePoints_.back().time_;
				double maxFrac = 0;
				std::string prevHapMax = "";
				for(const auto & prevHap : currentPatient.timePoints_.back().hapFracs_){
					if(prevHap.second > maxFrac){
						maxFrac = prevHap.second;
						prevHapMax = prevHap.first;
					}
					double nextFrac = std::max(0.0, prevHap.second + (prevHap.second * ((diffInTime/24) * (haplotypeKeyTorResistanceFactor[prevHap.first] * rGen.unifRand(0.75,1.25)))  ));;
					if(nextFrac  > 0){
						tPoint.hapFracs_[prevHap.first] = nextFrac;
					}
				}
				if(tPoint.hapFracs_.empty()){
					tPoint.hapFracs_[prevHapMax] = maxFrac;
				}
				currentPatient.timePoints_.emplace_back(tPoint);
			}
		}else{
			TimePoint initialTimePoint;
			initialTimePoint.time_ = 0;
			initialTimePoint.hapFracs_[haplotypes.front()] = 1;
			currentPatient.timePoints_.emplace_back(initialTimePoint);
			for(const auto & time : timePoints){
				TimePoint tPoint;
				tPoint.time_ = time;
				tPoint.hapFracs_[haplotypes.front()] = 1;
				currentPatient.timePoints_.emplace_back(tPoint);
			}
		}
	}
	OutOptions sampleNamesOutOpts(
			njh::files::make_path(setUp.pars_.directoryName_, "sampleNames.tab.txt"));
	VecStr sampleNamesTabHeader{};
	if (!lSetup.ids_->containsMids()) {
		sampleNamesTabHeader = {"#target", "sample", "run1"};
	}else{
		sampleNamesTabHeader = {"#file", "sample", "run1"};
	}
	bool anyReps = false;
	for (const auto & patient : patientSetUpPars) {
		if(patient.replicate_){
			anyReps = true;
			break;
		}
	}
	if(anyReps){
		sampleNamesTabHeader.emplace_back("run2");
	}


	auto maxStartingTemplateAmounts = vectorMaximum(pcrNumbers.startingTemplateAmounts_);
	auto maxFinalReadAmount = vectorMaximum(pcrNumbers.finalReadAmount_);
	auto genTimePointNameSection = [&maxTime,&timePoints](uint32_t time){
		std::string name;
		if(!timePoints.empty()){
			name = njh::pasteAsStr("-TP", njh::leftPadNumStr(time, maxTime));
		}
		return name;
	};

	auto genStartingTempNameSection = [&pcrNumbers,&maxStartingTemplateAmounts](uint32_t startingTemplateAmount){
		std::string name;
		if(pcrNumbers.startingTemplateAmounts_.size() > 1){
			name = njh::pasteAsStr("-ST", njh::leftPadNumStr(startingTemplateAmount, maxStartingTemplateAmounts));
		}
		return name;
	};
	auto genFinalReadNameSection = [&pcrNumbers,&maxFinalReadAmount](uint32_t finalReadAmount){
		std::string name;
		if(pcrNumbers.finalReadAmount_.size() > 1){
			name = njh::pasteAsStr("-RD", njh::leftPadNumStr(finalReadAmount, maxFinalReadAmount));
		}
		return name;
	};
	for (const auto & patient : patientSetUpPars) {
		for (const auto & tp : patient.timePoints_) {
			std::map<std::string, ControlPopulation::Strain> strains;
			double totalFracs = 0;
			for (const auto & hap : tp.hapFracs_) {
				totalFracs += hap.second;
			}
			for (const auto & hap : tp.hapFracs_) {
				strains.emplace(hap.first,
						ControlPopulation::Strain(hap.first, haplotypeKey[hap.first],
								hap.second/totalFracs));
			}
			ControlPopulation::Sample currentSamp(patient.name_ + genTimePointNameSection(tp.time_), strains);
			currentSamp.meta_.addMeta(patient.meta_, false);
			currentSamp.meta_.addMeta("Patient", patient.name_);
			if(!timePoints.empty()){
				currentSamp.meta_.addMeta("TimePoint", tp.time_);
			}

			for (const auto & startingTemplateAmount : pcrNumbers.startingTemplateAmounts_){
				for (const auto & finalReadAmount : pcrNumbers.finalReadAmount_){
					uint32_t numRuns = patient.replicate_ ? 2 : 1;
          uint32_t expectedStartingTemplateAmount = startingTemplateAmount;
          if(std::numeric_limits<uint32_t>::max() != patient.initialStartingTemplate_){
            expectedStartingTemplateAmount = patient.initialStartingTemplate_;
          }
					uint32_t startingTemplateAmountForSamp = expectedStartingTemplateAmount;
					if(std::numeric_limits<double>::max() != twoSdFrac){
						double startingTemplateAmountForSampNew = 0;
						std::normal_distribution<double> ndist(startingTemplateAmountForSamp, (startingTemplateAmountForSamp * twoSdFrac)/2);
				    // random device class instance, source of 'true' randomness for initializing random seed
				    std::random_device rd;
				    // Mersenne twister PRNG, initialized with seed from previous random device instance
				    std::mt19937 rgen(rd());
				    while(startingTemplateAmountForSampNew <= 1){
				    	startingTemplateAmountForSampNew = std::round(ndist(rgen));
				    }
				    startingTemplateAmountForSamp = static_cast<uint32_t>(startingTemplateAmountForSampNew);
					}
					currentSamp.addExpRun(
							genStartingTempNameSection(expectedStartingTemplateAmount)
									+ genFinalReadNameSection(finalReadAmount),
							ControlPopulation::SeqSampledAmounts { finalReadAmount,
						startingTemplateAmountForSamp }, numRuns);

				}
			}
			lSetup.pop_.addSample(currentSamp);
		}
	}

	table sampleNamesTables(sampleNamesTabHeader);
	if (!lSetup.ids_->containsMids()) {
		auto sampleKeys = njh::getVecOfMapKeys(lSetup.pop_.samples_);
		njh::sort(sampleKeys);
		for(const auto & sampleKey : sampleKeys){
			const auto & sample = lSetup.pop_.samples_.at(sampleKey);
			for(const auto & experiment : sample.expRuns_){
				std::string outSampName = sample.sampName_ + experiment.first;
				//add to sample name table
				for(const auto & primerPair : primerPairNames){
					VecStr addingRow{primerPair, outSampName};
					for(const auto & run : experiment.second){
						addingRow.emplace_back(sample.sampName_ + run.runName_);
					}
					sampleNamesTables.addRow(addingRow);
				}
				//add experiments to sequence setup;
				uint32_t runNumber = 0;
				for(const auto & run : experiment.second){
					auto sampleSet = std::make_shared<SampleSetup>(sample.sampName_ + run.runName_);
					std::unordered_map<std::string, std::shared_ptr<MixtureSetUp> > mixtures;
					for (const auto & primerPair : primerPairNames) {
						mixtures[primerPair] = std::make_shared<MixtureSetUp>(primerPair);
						mixtures[primerPair]->meta_ = std::make_unique<MetaDataInName>();
						mixtures[primerPair]->meta_->addMeta("PrimerPair", primerPair);
						mixtures[primerPair]->meta_->addMeta("PatientSample", sample.sampName_ + experiment.first);
						//add in identifying information
						mixtures[primerPair]->meta_->addMeta("sampName_", sample.sampName_);
						mixtures[primerPair]->meta_->addMeta("experiment", experiment.first );
						mixtures[primerPair]->meta_->addMeta("runNumber", runNumber);

						mixtures[primerPair]->startingTemplateAmount_ = run.expAmounts_.totalGenomesSampled_;
						mixtures[primerPair]->finalReadAmount_ = run.expAmounts_.sequencedReadAmount_;
						mixtures[primerPair]->setPrimers(primerPair,
								njh::mapAt(lSetup.ids_->targets_, primerPair).info_.forwardPrimerRaw_,
								njh::mapAt(lSetup.ids_->targets_, primerPair).info_.reversePrimerRaw_);
						mixtures[primerPair]->primers_->reverse_randomPrecedingBases_ = simPars.primerRandomPrecedingBases_;
						mixtures[primerPair]->primers_->forward_randomPrecedingBases_ = simPars.primerRandomPrecedingBases_;
					}
//					std::cout << sampleKey << std::endl;
//					std::cout << "sample.strainsExpected_: " << sample.strainsExpected_.size() << std::endl;
//					std::cout << sample.sampName_ << run.runName_ << std::endl;
//					std::cout << "\tsample.strainsExpected_.size():" << sample.strainsExpected_.size() << std::endl;
					for (const auto & strain : sample.strainsExpected_) {
//						std::cout << strain.first << std::endl;
//						std::cout << "\tstrain.second.name_: " << strain.second.name_ << std::endl;
//						std::cout << "\tstrain.second.relativeAbundance_: " << strain.second.relativeAbundance_ << std::endl;
//						std::cout << "\tstrain.second.genomicRegionToHapNames_.size(): " << strain.second.genomicRegionToHapNames_.size() << std::endl;
						uint32_t genomeCountsForStrain = run.hapAbundGenomesSampled_.at(strain.second.name_).relativeAbundance_;
//						std::cout << strain.first << "\t" << genomeCountsForStrain << std::endl;
						for (const auto & subHapName : strain.second.genomicRegionToHapNames_) {
//							std::cout << strain.first << std::endl;
//							std::cout << '\t' << subHapName.first << std::endl;
//							std::cout << "\t\t" << subHapName.second << std::endl;
//							std::cout << "\t\t" << "strain.second.relativeAbundance_: " << strain.second.relativeAbundance_ << std::endl;
							if (njh::in(subHapName.second, mixtures[subHapName.first]->expectedAbundances_)) {
								mixtures[subHapName.first]->expectedAbundances_[subHapName.second] += strain.second.relativeAbundance_;
								mixtures[subHapName.first]->genomeCounts_[subHapName.second] += genomeCountsForStrain;
							} else {
								mixtures[subHapName.first]->addAbundance(subHapName.second, strain.second.relativeAbundance_);
								mixtures[subHapName.first]->addGenomeCount(subHapName.second, genomeCountsForStrain);
							}
						}
					}
					for (const auto & mix : mixtures) {
//						std::cout << sample.sampName_ << run.runName_ << "-" << mix.first << std::endl;
//						for(const auto & hap : mix.second->genomeCounts_){
//							std::cout << "\t" << hap.first << "\t" << hap.second << std::endl;
//						}
						sampleSet->addMixture(mix.second);
					}
					lSetup.addSample(sampleSet);
					++runNumber;
				}
			}
		}

	}else{
		uint32_t finalSubSetAmount = 0;
		{
			uint32_t midCount = 0;
			for (const auto & patient : lSetup.pop_.samples_) {
				for(const auto & experiment : patient.second.expRuns_){
					for(uint32_t run = 0; run < experiment.second.size(); ++run){
						++midCount;
						if(midCount >= lSetup.ids_->mids_.size()){
							midCount = 0;
							++finalSubSetAmount;
						}
					}
				}
			}
		}
		uint32_t indexCount = 0;
		uint32_t midCount = 0;
		auto sampleSet = std::make_shared<SampleSetup>(njh::pasteAsStr("Subset-", njh::leftPadNumStr(indexCount, finalSubSetAmount) ) );
		auto midNames = getVectorOfMapKeys(lSetup.ids_->mids_);
		njh::sort(midNames);
		auto sampleKeys = njh::getVecOfMapKeys(lSetup.pop_.samples_);
		njh::sort(sampleKeys);

		for(const auto & sampleKey : sampleKeys){
			const auto & sample = lSetup.pop_.samples_.at(sampleKey);
			for(const auto & experiment : sample.expRuns_){
				std::string outSampName = sample.sampName_ + experiment.first;
				std::unordered_map<std::string, VecStr> replicateNames;
				uint32_t runNumber = 0;
				for(const auto & run : experiment.second){
					std::string midName =  midNames[midCount];
					std::string indexName = sampleSet->name_;
					replicateNames[indexName].emplace_back(midName);
					std::unordered_map<std::string, std::shared_ptr<MixtureSetUp> > mixtures;
					for (const auto & primerPair : primerPairNames) {
						mixtures[primerPair] = std::make_shared<MixtureSetUp>(njh::pasteAsStr(primerPair, "-", midName));
						mixtures[primerPair]->meta_ = std::make_unique<MetaDataInName>();
						mixtures[primerPair]->meta_->addMeta("PrimerPair", primerPair);
						mixtures[primerPair]->meta_->addMeta("PatientSample", outSampName);
						//add in identifying information
						mixtures[primerPair]->meta_->addMeta("sampName_", sample.sampName_);
						mixtures[primerPair]->meta_->addMeta("experiment", experiment.first );
						mixtures[primerPair]->meta_->addMeta("runNumber", runNumber);

						mixtures[primerPair]->startingTemplateAmount_ = run.expAmounts_.totalGenomesSampled_;
						mixtures[primerPair]->finalReadAmount_ = run.expAmounts_.sequencedReadAmount_;
						mixtures[primerPair]->setPrimers(primerPair,
								njh::mapAt(lSetup.ids_->targets_, primerPair).info_.forwardPrimerRaw_,
								njh::mapAt(lSetup.ids_->targets_, primerPair).info_.reversePrimerRaw_);
						mixtures[primerPair]->primers_->reverse_randomPrecedingBases_ = simPars.primerRandomPrecedingBases_;
						mixtures[primerPair]->primers_->forward_randomPrecedingBases_ = simPars.primerRandomPrecedingBases_;
					}
					for (const auto & strain : sample.strainsExpected_) {
						for (const auto & subHapName : strain.second.genomicRegionToHapNames_) {
							uint32_t genomeCountsForSubRegion = run.hapRegionAmplified_.at(subHapName.first).hapAbundSampled_.at(subHapName.second);
							if (njh::in(subHapName.second, mixtures[subHapName.first]->expectedAbundances_)) {
								mixtures[subHapName.first]->expectedAbundances_[subHapName.second] += strain.second.relativeAbundance_;
								mixtures[subHapName.first]->genomeCounts_[subHapName.second] += genomeCountsForSubRegion;
							} else {
								mixtures[subHapName.first]->addAbundance(subHapName.second, strain.second.relativeAbundance_);
								mixtures[subHapName.first]->addGenomeCount(subHapName.second, genomeCountsForSubRegion);
							}
						}
					}
					for (auto & mix : mixtures) {
						if(nullptr != lSetup.ids_->mids_.at(midName).forwardBar_){
							mix.second->setForwardBarcode(midName, lSetup.ids_->mids_.at(midName).forwardBar_->bar_->motifOriginal_);
							mix.second->forwardBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
						}
						if(nullptr != lSetup.ids_->mids_.at(midName).reverseBar_){
							mix.second->setReverseBarcode(midName, lSetup.ids_->mids_.at(midName).reverseBar_->bar_->motifOriginal_);
							mix.second->reverseBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
						}
						sampleSet->addMixture(mix.second);
					}
					++midCount;
					if(midCount >= lSetup.ids_->mids_.size()){
						lSetup.addSample(sampleSet);
						midCount = 0;
						++indexCount;
						sampleSet = std::make_shared<SampleSetup>(njh::pasteAsStr("Subset-", njh::leftPadNumStr(indexCount, finalSubSetAmount) ) );
					}
					++runNumber;
				}
				for(const auto & index : replicateNames){
					VecStr addingRow{index.first, outSampName};
					for(const auto & run : index.second){
						addingRow.emplace_back(run);
					}
					sampleNamesTables.addRow(addingRow);
				}
			}
		}
		if(midCount > 0){
			lSetup.addSample(sampleSet);
		}
	}
	njh::sort(sampleNamesTables.content_, [](const VecStr & row1, const VecStr & row2){
		if(row1[0] == row2[0]){
			if(row1[1] == row2[1]){
				return row1[2] < row2[2];
			}else{
				return row1[1] < row2[1];
			}
		}else{
			return row1[0] < row2[0];
		}
	});

	OutputStream sampleNamesOut(sampleNamesOutOpts);
	sampleNamesTables.outPutContents(sampleNamesOut, "\t");

	VecStr metaHeader{"sample"};
	if(!timePoints.empty() || pcrNumbers.finalReadAmount_.size() > 1 || pcrNumbers.startingTemplateAmounts_.size() > 1){
		metaHeader.emplace_back("Patient");
	}
	if(!timePoints.empty()){
		metaHeader.emplace_back("TimePoint");
	}
	metaHeader.emplace_back("PCRStartingTemplate");
	metaHeader.emplace_back("FinalSamplingReadAmount");
	addOtherVec(metaHeader, metaLevels);
	table metaTable(metaHeader);
	for(const auto & patient : patientSetUpPars){
		for(const auto & tp : patient.timePoints_){
			for (const auto & iterStartingTemplateAmount : pcrNumbers.startingTemplateAmounts_){
				for (const auto & finalReadAmount : pcrNumbers.finalReadAmount_){
          uint32_t startingTemplateAmount = iterStartingTemplateAmount;
          if(std::numeric_limits<uint32_t>::max() != patient.initialStartingTemplate_){
            startingTemplateAmount = patient.initialStartingTemplate_;
          }

					MetaDataInName sampMeta = patient.meta_;
					std::string outputSampName = patient.name_
							+ genTimePointNameSection(tp.time_)
							+ genStartingTempNameSection(startingTemplateAmount)
							+ genFinalReadNameSection(finalReadAmount);
					VecStr row;
					row.emplace_back(outputSampName);
					if(!timePoints.empty() || pcrNumbers.finalReadAmount_.size() > 1 || pcrNumbers.startingTemplateAmounts_.size() > 1){
						sampMeta.addMeta("Patient", patient.name_);
						row.emplace_back(estd::to_string(patient.name_));
					}
					if(!timePoints.empty()){
						sampMeta.addMeta("TimePoint", tp.time_);
						row.emplace_back(estd::to_string(tp.time_));
					}
					sampMeta.addMeta("PCRStartingTemplate", startingTemplateAmount);
					sampMeta.addMeta("FinalSamplingReadAmount", finalReadAmount);
					row.emplace_back(estd::to_string(startingTemplateAmount));
					row.emplace_back(estd::to_string(finalReadAmount));
					for(const auto & head : metaLevels){
						row.emplace_back(sampMeta.getMeta(head));
					}
					metaTable.addRow(row);
				}
			}
		}
	}
	OutOptions libSetUpOpts(njh::files::make_path(setUp.pars_.directoryName_, "librarySetup.json"));
	OutputStream libSetUpOut(libSetUpOpts);
	libSetUpOut << njh::json::toJson(lSetup) << std::endl;

	OutOptions metaOutOpts(njh::files::make_path(setUp.pars_.directoryName_, "metaData.tab.txt"));
	OutputStream metaOut(metaOutOpts);
	metaTable.sortTable("sample", false);
	metaTable.outPutContents(metaOut, "\t");

	auto abundTab = lSetup.createAbundanceTable();
	OutOptions abundTabOutOpts(njh::files::make_path(setUp.pars_.directoryName_, "abundanceTable.tab.txt"));
	OutputStream abundTabOut(abundTabOutOpts);
	njh::sort(abundTab.content_, [](const VecStr & row1, const VecStr & row2){
		if(row1[1] == row2[1]){
			if(row1[2] == row2[2]){
				return row1[3] < row2[3];
			}else{
				return row1[2] < row2[2];
			}
		}else{
			return row1[1] < row2[1];
		}
	});
	abundTab.outPutContents(abundTabOut, "\t");
	return 0;
}






} // namespace njhseq


