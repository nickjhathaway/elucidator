/*
 * readSimulatorRunner_simMixture.cpp
 *
 *  Created on: Sep 16, 2018
 *      Author: nick
 */
// elucidator - A library for analyzing sequence data
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
//
// This file is part of elucidator.
//
// elucidator is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// elucidator is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with elucidator.  If not, see <http://www.gnu.org/licenses/>.
//
#include "readSimulatorRunner.hpp"
#include "elucidator/simulation.h"
#include <SeekDeep/objects/PrimersAndMids.hpp>
namespace njhseq {


int readSimulatorRunner::createIlluminaErrorProfile(
		const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bedFnp = "";
	bfs::path twoBitFnp = "";
	uint32_t numThreads = 1;
	uint32_t lengthLimit = 0;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processReadInNames(VecStr{"--bam"});
	setUp.processDirectoryOutputName(true);
	setUp.setOption(lengthLimit, "--lengthLimit", "Length Limit");
	setUp.setOption(numThreads, "--numThreads", "Number Threads");
	setUp.setOption(bedFnp, "--bed", "Bed file");
	setUp.setOption(twoBitFnp, "--twoBitFnp", "2bit file for genome aligned to", true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	BamTools::BamReader bReader;
	bReader.Open(setUp.pars_.ioOptions_.firstName_.string());
	checkBamOpenThrow(bReader, setUp.pars_.ioOptions_.firstName_);
	loadBamIndexThrow(bReader);
	auto refs = bReader.GetReferenceData();
	BamAlnsCache alnCache;
	auto refData = bReader.GetReferenceData();
	std::unordered_map<std::string, uint32_t> refNameToId;
	for (auto pos : iter::range(refData.size())) {
		refNameToId[refData[pos].RefName] = pos;
	}

	TwoBit::TwoBitFile tReader(twoBitFnp);
	RoughIlluminaProfiler profiler;
	if("" != bedFnp){
		auto regions = bedPtrsToGenomicRegs(getBeds(bedFnp));
		njh::concurrent::LockableQueue<GenomicRegion> regionsQueue(regions);
		concurrent::BamReaderPool bamPools(setUp.pars_.ioOptions_.firstName_, numThreads);
		bamPools.openBamFile();
		std::mutex mut;
		auto increaseCount = [&regionsQueue,&bamPools,&mut,&profiler,&twoBitFnp,&refData,&setUp,&lengthLimit](){
			GenomicRegion region;
			TwoBit::TwoBitFile tReader(twoBitFnp);
			auto bReader = bamPools.popReader();
			BamTools::BamAlignment bAln;
			RoughIlluminaProfiler currentProfiler;
			while(regionsQueue.getVal(region)){
				if(setUp.pars_.verbose_){
					std::lock_guard<std::mutex> lock(mut);
					std::cout << region.uid_ << std::endl;
				}
				setBamFileRegionThrow(*bReader, region);
				while (bReader->GetNextAlignment(bAln)) {
					if(bAln.IsMapped() && bAln.IsPrimaryAlignment()){
						if(bAln.GetEndPosition() - bAln.Position > lengthLimit){
							currentProfiler.increaseCounts(bAln, refData, tReader);
						}
					}
				}
			}
			{
				std::lock_guard<std::mutex> lock(mut);
				profiler.addOther(currentProfiler);
			}
		};
		std::vector<std::thread> threads;
		for(uint32_t t = 0; t < numThreads; ++t){
			threads.emplace_back(std::thread(increaseCount));
		}
		njh::concurrent::joinAllJoinableThreads(threads);
	}else{
		BamTools::BamAlignment bAln;
		while (bReader.GetNextAlignment(bAln)) {
			if(bAln.IsMapped() && bAln.IsPrimaryAlignment()){
				if(bAln.GetEndPosition() - bAln.Position > lengthLimit){
					profiler.increaseCounts(bAln, refData, tReader);
				}
			}
		}
	}
	profiler.r1_counts.writeProfiles(njh::files::make_path(setUp.pars_.directoryName_, "r1").string(), false);
	profiler.r2_counts.writeProfiles(njh::files::make_path(setUp.pars_.directoryName_, "r2").string(), false);
	return 0;
}

struct PCRAmountPars {
	PCRAmountPars() {
	}
	PCRAmountPars(const std::vector<uint32_t> & startingTemplateAmounts,
			const std::vector<uint32_t> & finalReadAmount) :
			startingTemplateAmounts_(startingTemplateAmounts),
			finalReadAmount_(finalReadAmount) {
	}
	PCRAmountPars(const uint32_t & startingTemplateAmount,
			const uint32_t & finalReadAmount) :
			startingTemplateAmounts_({startingTemplateAmount}),
			finalReadAmount_({finalReadAmount}) {
	}
	std::vector<uint32_t> startingTemplateAmounts_ { std::vector<uint32_t>{ 2048 } };
	std::vector<uint32_t> finalReadAmount_ { std::vector<uint32_t>{ 5000 } };
};
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
	readSimulatorSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processVerbose();
	setUp.setOption(pcrNumbers.startingTemplateAmounts_, "--startingTemplateAmount", "Starting PCR Template Amount", njh::progutils::ProgramSetUp::ConCheckCase::NONZERO);
	setUp.setOption(pcrNumbers.finalReadAmount_, "--perMixtureReadAmount", "Final Read Amount to create per mixture", njh::progutils::ProgramSetUp::ConCheckCase::NONZERO);
	setUp.setOption(simReplicates, "--replicates", "Replicates");
	setUp.setOption(adjustAroundGivenFrac, "--adjustAroundGivenFrac", "Adjust Around Given Frac", njh::progutils::ProgramSetUp::CheckCase::GREATERZERO);
	adjustAroundGivenFrac = std::min(adjustAroundGivenFrac, 0.99);
	setUp.setOption(timePoints, "--timePoints", "Additional time points, will automatically have time point 00", njh::progutils::ProgramSetUp::ConCheckCase::NONZERO);
	setUp.setOption(patientSetupFile, "--patientSetupFile", "Patient Setup File", true);
	setUp.setOption(rawPatientSetupFile, "--rawPatientSetupFile", "Patient Setup File designates individual samples");
	setUp.setOption(coiTableFnp, "--coiTable", "COI Table", true);
	setUp.setOption(haplotypeInfo, "--haplotypeInfo", "Haplotype Info", true);
	setUp.setOption(primerMidFnp, "--primerMidFnp", "Primer MID Fnp", true);
	setUp.setOption(simPars.pairedEndLength_, "--pairedEndLength", "Paired End Length");
	setUp.setOption(simPars.noAddPrimers_, "--noAddPrimers", "Primers are already present");
	setUp.setOption(simPars.barcodeRandomPrecedingBases_, "--barcodeRandomPrecedingBases", "Barcode Random Preceding Bases");
	setUp.setOption(simPars.primerRandomPrecedingBases_, "--primerRandomPrecedingBases", "Primer Random Preceding Bases");
	setUp.setOption(simPars.addReverseComplement_, "--addReverseComplement", "Add Reverse Complement");
	setUp.setOption(simPars.addBluntEndingArtifact_, "--addBluntEndingArtifact", "Add Blunt Ending Artifact");
	setUp.setOption(simPars.bluntEndingArtifactChance_, "--bluntEndingArtifactChance", "Blunt Ending Artifact Chance");
	setUp.setOption(libraryName, "--libraryName", "Library Name", true);
	setUp.processDirectoryOutputName(libraryName, true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_ );
	njh::sort(timePoints);
	uint32_t maxTime = 0;
	if(!timePoints.empty()){
		maxTime = *std::max_element(timePoints.begin(), timePoints.end());
	}
	LibrarySetup lSetup(libraryName, simPars);

	table coiTable(coiTableFnp, "\t", true);
	table haplotypeInfoTable(haplotypeInfo, "\t", true);
	table patientSetupTable(patientSetupFile, "\t", true);
	lSetup.ids_ = std::make_unique<PrimersAndMids>(primerMidFnp);
	lSetup.ids_->initPrimerDeterminator();
	if(lSetup.ids_->containsMids()){
		lSetup.ids_->initMidDeterminator(MidDeterminator::MidDeterminePars{});
	}
	if(rawPatientSetupFile){
	}else{
		//checks
		patientSetupTable.checkForColumnsThrow(VecStr{"number"}, __PRETTY_FUNCTION__);
		if(patientSetupTable.nCol() < 2){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ <<  " error, there should be more than 1 column" << "\n";
			throw std::runtime_error{ss.str()};
		}
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
	}else{
		removeElement<std::string>(metaLevels, "number");
	}
	haplotypeInfoTable.checkForColumnsThrow(njh::catVecs(metaLevels, VecStr{"frac"}), __PRETTY_FUNCTION__);
	coiTable.checkForColumnsThrow(njh::catVecs(metaLevels, VecStr{"COI", "n"}), __PRETTY_FUNCTION__);
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
	checkForUniqLevels(coiTable, njh::catVecs(metaLevels, VecStr{"COI"}), coiTableFnp);
	auto levels = getLevelsForHeader(patientSetupTable, metaLevels);
	checkForLevelsForHeader(haplotypeInfoTable, levels, metaLevels, haplotypeInfo);
	checkForLevelsForHeader(coiTable, levels, metaLevels, coiTableFnp);
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

//	for(const auto pos : iter::range(haplotypeInfoTable.columnNames_.size())){
//		std::cout << "pos: " << pos << " vs " << haplotypeInfoTable.getColPos(haplotypeInfoTable.columnNames_[pos]) << " for " << haplotypeInfoTable.columnNames_[pos] << std::endl;
//	}

	for(const auto & row : haplotypeInfoTable){
		VecStr currentLevelVec;
		for (const auto & col : metaLevels) {
			currentLevelVec.emplace_back(row[haplotypeInfoTable.getColPos(col)]);
		}
		std::string currentLevel = njh::pasteAsStr(currentLevelVec);
		VecStr hapNames;
		for(const auto & primerpair : primerPairNames){
//			std::cout << "primerpair: " << primerpair << std::endl;
//			std::cout << "\t" << "haplotypeInfoTableNames.getColPos(primerpair): " << haplotypeInfoTable.getColPos(primerpair) << std::endl;
//			std::cout << "\t" << "row[haplotypeInfoTableNames.getColPos(primerpair)] : " << row[haplotypeInfoTable.getColPos(primerpair)] << std::endl;
			hapNames.emplace_back(row[haplotypeInfoTable.getColPos(primerpair)] );
		}
		njh::sort(hapNames);
		std::string hapUID = njh::conToStr(hapNames, "__");
//		std::cout << "hapUID: " << hapUID << std::endl;

		//exit(1);
		double frac = njh::StrToNumConverter::stoToNum<double>(row[haplotypeInfoTable.getColPos("frac")]);
		hapFracForLevels[currentLevel][hapUID] = frac;
		hapCounts[currentLevel].first.emplace_back(hapUID);
		hapCounts[currentLevel].second.emplace_back(frac);
	}
	for(const auto & hap : hapCounts){
		hapRGenerator.emplace(hap.first, njh::randObjectGen<std::string, double>(hap.second.first, hap.second.second));
	}
	struct TimePoint {
		uint32_t time_;
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
		Json::Value toJson() const{
			Json::Value ret;
			ret["class"] = njh::json::toJson(njh::getTypeName(*this));
			ret["name_"] = njh::json::toJson(name_);
			ret["meta_"] = njh::json::toJson(meta_.meta_);
			ret["timePoints_"] = njh::json::toJson(timePoints_);
			ret["replicate_"] = njh::json::toJson(replicate_);
			ret["initalCOI_"] = njh::json::toJson(initalCOI_);
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
		bool hasNameCol = false;
		if(njh::in<std::string>("name", patientSetupTable.columnNames_)){
			hasNameCol = true;
		}
		bool hasReplicateCol = false;
		if(njh::in<std::string>("replicate", patientSetupTable.columnNames_)){
			hasReplicateCol = true;
		}
		bool hasCOICol = false;
		if(njh::in<std::string>("COI", patientSetupTable.columnNames_)){
			hasCOICol = true;
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
			uint32_t initalCOI = njh::mapAt(coiRGenerator, currentLevel).genObj();
			if(hasCOICol){
				initalCOI = njh::StrToNumConverter::stoToNum<uint32_t>(row[patientSetupTable.getColPos("COI")]);
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
//				std::cout << "hap.first: " << hap.first << std::endl;
//				std::cout << "hap.second/totalFracs: " << hap.second/totalFracs << std::endl;
//				std::cout << "njh::in(hap.first, haplotypeKey): " << njh::colorBool(njh::in(hap.first, haplotypeKey)) << std::endl;
//				std::cout << "haplotypeKey.size(): " << haplotypeKey.size() << std::endl;
//				std::cout << "hap.first.size(): " << hap.first.size() << std::endl;
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
					currentSamp.addExpRun(
							genStartingTempNameSection(startingTemplateAmount)
									+ genFinalReadNameSection(finalReadAmount),
							ControlPopulation::SeqSampledAmounts { finalReadAmount,
									startingTemplateAmount }, numRuns);
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
								njh::mapAt(lSetup.ids_->targets_, primerPair).info_.forwardPrimer_,
								njh::mapAt(lSetup.ids_->targets_, primerPair).info_.reversePrimer_);
						mixtures[primerPair]->primers_->reverse_randomPrecedingBases_ = simPars.primerRandomPrecedingBases_;
						mixtures[primerPair]->primers_->forward_randomPrecedingBases_ = simPars.primerRandomPrecedingBases_;
					}
//					std::cout << sampleKey << std::endl;
//					std::cout << "sample.strainsExpected_: " << sample.strainsExpected_.size() << std::endl;

					for (const auto & strain : sample.strainsExpected_) {
//						std::cout << strain.first << std::endl;
//						std::cout << "\tstrain.second.name_: " << strain.second.name_ << std::endl;
//						std::cout << "\tstrain.second.relativeAbundance_: " << strain.second.relativeAbundance_ << std::endl;
//						std::cout << "\tstrain.second.genomicRegionToHapNames_.size(): " << strain.second.genomicRegionToHapNames_.size() << std::endl;
						for (const auto & subHapName : strain.second.genomicRegionToHapNames_) {
							uint32_t genomeCountsForSubRegion = run.hapRegionAmplified_.at(subHapName.first).hapAbundSampled_.at(subHapName.second);
//							std::cout << strain.first << std::endl;
//							std::cout << '\t' << subHapName.first << std::endl;
//							std::cout << "\t" << genomeCountsForSubRegion << std::endl;

							if (njh::in(subHapName.second, mixtures[subHapName.first]->expectedAbundances_)) {
								mixtures[subHapName.first]->expectedAbundances_[subHapName.second] += strain.second.relativeAbundance_;
								mixtures[subHapName.first]->genomeCounts_[subHapName.second] += genomeCountsForSubRegion;
							} else {
								mixtures[subHapName.first]->addAbundance(subHapName.second, strain.second.relativeAbundance_);
								mixtures[subHapName.first]->addGenomeCount(subHapName.second, genomeCountsForSubRegion);
							}
						}
					}
					for (const auto & mix : mixtures) {
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
								njh::mapAt(lSetup.ids_->targets_, primerPair).info_.forwardPrimer_,
								njh::mapAt(lSetup.ids_->targets_, primerPair).info_.reversePrimer_);
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
			for (const auto & startingTemplateAmount : pcrNumbers.startingTemplateAmounts_){
				for (const auto & finalReadAmount : pcrNumbers.finalReadAmount_){
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


int createLibrarySimMultipleMixtureExampleTesting(
		const njh::progutils::CmdArgs & inputCommands) {
	LibrarySetup::SimLibrarySetupPars simPars;
	OutOptions outOpts(bfs::path(""), ".json");
	readSimulatorSetUp setUp(inputCommands);
	setUp.processWritingOptions(outOpts);
	setUp.setOption(simPars.barcodeRandomPrecedingBases_, "--barcodeRandomPrecedingBases", "Barcode Random Preceding Bases");
	setUp.setOption(simPars.addReverseComplement_, "--addReverseComplement", "Add Reverse Complement");
	setUp.setOption(simPars.addBluntEndingArtifact_, "--addBluntEndingArtifact", "Add Blunt Ending Artifact");
	setUp.setOption(simPars.bluntEndingArtifactChance_, "--bluntEndingArtifactChance", "Blunt Ending Artifact Chance");
	setUp.finishSetUp(std::cout);
	OutputStream out(outOpts);
	LibrarySetup lSetup("testing_barcode_scheme", simPars);
	{
		auto sampleSet = std::make_shared<SampleSetup>("Sample1-FrontBarcode");
		{
			auto mixture = std::make_shared<MixtureSetUp>("Mixture1");
			mixture->setPrimers("ama1_S0_Sub0_mip2", "CAATATTTAAAAGATGGAGGTTTTGCT",
					"ATCTTCACAATTTCCATCGACCCAT");
			mixture->addAbundance("PfML01", 10);
			mixture->addAbundance("PfDd2", 10);
			mixture->addAbundance("PfGB4-PfTG01", 10);
			mixture->addAbundance("PfKH02", 10);
			mixture->setForwardBarcode("MID07", "CGTGTCTCTA");
			mixture->forwardBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
			sampleSet->addMixture(mixture);
		}
		{
			auto mixture = std::make_shared<MixtureSetUp>("Mixture2");
			mixture->setPrimers("ama1_S0_Sub0_mip2", "CAATATTTAAAAGATGGAGGTTTTGCT",
					"ATCTTCACAATTTCCATCGACCCAT");
			mixture->addAbundance("PfML01", 10);
			mixture->addAbundance("PfDd2", 20);
			mixture->addAbundance("PfGB4-PfTG01", 30);
			mixture->addAbundance("PfKH02", 40);
			mixture->setForwardBarcode("MID08", "CTCGCGTGTC");
			mixture->forwardBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
			sampleSet->addMixture(mixture);
		}

		{
			auto mixture = std::make_shared<MixtureSetUp>("Mixture3");
			mixture->setPrimers("ama1_S0_Sub0_mip2", "CAATATTTAAAAGATGGAGGTTTTGCT",
					"ATCTTCACAATTTCCATCGACCCAT");
			mixture->addAbundance("PfML01", 10);
			mixture->addAbundance("PfDd2", 20);
			mixture->addAbundance("PfGB4-PfTG01", 30);
			mixture->addAbundance("PfKH02", 40);
			mixture->setForwardBarcode("MID10", "TCTCTATGCG");
			mixture->forwardBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
			sampleSet->addMixture(mixture);
		}
		{
			auto mixture = std::make_shared<MixtureSetUp>("Mixture4");
			mixture->setPrimers("ama1_S0_Sub0_mip2", "CAATATTTAAAAGATGGAGGTTTTGCT",
					"ATCTTCACAATTTCCATCGACCCAT");
			mixture->addAbundance("PfML01", 10);
			mixture->addAbundance("PfDd2", 20);
			mixture->addAbundance("PfGB4-PfTG01", 30);
			mixture->addAbundance("PfKH02", 40);
			mixture->setForwardBarcode("MID11", "TGATACGTCT");
			mixture->forwardBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
			sampleSet->addMixture(mixture);
		}
		{
			auto mixture = std::make_shared<MixtureSetUp>("Mixture5");
			mixture->setPrimers("ama1_S0_Sub0_mip2", "CAATATTTAAAAGATGGAGGTTTTGCT",
					"ATCTTCACAATTTCCATCGACCCAT");
			mixture->addAbundance("PfML01", 2);
			mixture->addAbundance("PfDd2", 3);
			mixture->addAbundance("PfGB4-PfTG01", 5);
			mixture->addAbundance("PfKH02", 90);
			mixture->setForwardBarcode("MID01", "ACGAGTGCGT");
			mixture->forwardBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
			sampleSet->addMixture(mixture);
		}
		{
			auto mixture = std::make_shared<MixtureSetUp>("Mixture6");
			mixture->setPrimers("ama1_S0_Sub0_mip2", "CAATATTTAAAAGATGGAGGTTTTGCT",
					"ATCTTCACAATTTCCATCGACCCAT");
			mixture->addAbundance("PfML01", 1);
			mixture->addAbundance("PfDd2", 1);
			mixture->addAbundance("PfGB4-PfTG01", 1);
			mixture->addAbundance("PfKH02", 99);
			mixture->setForwardBarcode("MID29", "ACTGTACAGT");
			mixture->forwardBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
			sampleSet->addMixture(mixture);
		}
		lSetup.addSample(sampleSet);
	}
//	{
//		auto sampleSet = std::make_shared<SampleSetup>("Sample2-ReverseBarcode");
//		{
//			auto mixture = std::make_shared<MixtureSetUp>("Mixture1");
//			mixture->setPrimers("ama1_S0_Sub0_mip2", "CAATATTTAAAAGATGGAGGTTTTGCT",
//					"ATCTTCACAATTTCCATCGACCCAT");
//			mixture->addAbundance("PfML01", 10);
//			mixture->addAbundance("PfDd2", 10);
//			mixture->addAbundance("PfGB4-PfTG01", 10);
//			mixture->addAbundance("PfKH02", 10);
//
//			mixture->setReverseBarcode("MID07", "CGTGTCTCTA");
//			sampleSet->addMixture(mixture);
//		}
//
//		{
//			auto mixture = std::make_shared<MixtureSetUp>("Mixture2");
//			mixture->setPrimers("ama1_S0_Sub0_mip2", "CAATATTTAAAAGATGGAGGTTTTGCT",
//					"ATCTTCACAATTTCCATCGACCCAT");
//			mixture->addAbundance("PfML01", 10);
//			mixture->addAbundance("PfDd2", 20);
//			mixture->addAbundance("PfGB4-PfTG01", 30);
//			mixture->addAbundance("PfKH02", 40);
//			mixture->setReverseBarcode("MID08", "CTCGCGTGTC");
//			sampleSet->addMixture(mixture);
//		}
//		{
//			auto mixture = std::make_shared<MixtureSetUp>("Mixture3");
//			mixture->setPrimers("ama1_S0_Sub0_mip2", "CAATATTTAAAAGATGGAGGTTTTGCT",
//					"ATCTTCACAATTTCCATCGACCCAT");
//			mixture->addAbundance("PfML01", 10);
//			mixture->addAbundance("PfDd2", 20);
//			mixture->addAbundance("PfGB4-PfTG01", 30);
//			mixture->addAbundance("PfKH02", 40);
//			mixture->setReverseBarcode("MID10", "TCTCTATGCG");
//
//			sampleSet->addMixture(mixture);
//		}
//
//		{
//			auto mixture = std::make_shared<MixtureSetUp>("Mixture4");
//			mixture->setPrimers("ama1_S0_Sub0_mip2", "CAATATTTAAAAGATGGAGGTTTTGCT",
//					"ATCTTCACAATTTCCATCGACCCAT");
//			mixture->addAbundance("PfML01", 10);
//			mixture->addAbundance("PfDd2", 20);
//			mixture->addAbundance("PfGB4-PfTG01", 30);
//			mixture->addAbundance("PfKH02", 40);
//			mixture->setReverseBarcode("MID11", "TGATACGTCT");
//			sampleSet->addMixture(mixture);
//		}
//		lSetup.addSample(sampleSet);
//	}
	{
		auto sampleSet = std::make_shared<SampleSetup>("Sample3-BarcodeBothEndsDifferent");
		{
			auto mixture = std::make_shared<MixtureSetUp>("Mixture1");
			mixture->setPrimers("ama1_S0_Sub0_mip2", "CAATATTTAAAAGATGGAGGTTTTGCT",
					"ATCTTCACAATTTCCATCGACCCAT");
			mixture->addAbundance("PfML01", 10);
			mixture->addAbundance("PfDd2", 10);
			mixture->addAbundance("PfGB4-PfTG01", 10);
			mixture->addAbundance("PfKH02", 10);
			mixture->setForwardBarcode("MID07", "CGAGAGATAC");
			mixture->setReverseBarcode("MID07", "CGTGTCTCTA");
			mixture->forwardBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
			mixture->reverseBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
			sampleSet->addMixture(mixture);
		}
		{
			auto mixture = std::make_shared<MixtureSetUp>("Mixture2");
			mixture->setPrimers("ama1_S0_Sub0_mip2", "CAATATTTAAAAGATGGAGGTTTTGCT",
					"ATCTTCACAATTTCCATCGACCCAT");
			mixture->addAbundance("PfML01", 10);
			mixture->addAbundance("PfDd2", 20);
			mixture->addAbundance("PfGB4-PfTG01", 30);
			mixture->addAbundance("PfKH02", 40);
			mixture->setForwardBarcode("MID08", "TCACGTACTA");
			mixture->setReverseBarcode("MID08", "CTCGCGTGTC");
			mixture->forwardBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
			mixture->reverseBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
			sampleSet->addMixture(mixture);
		}
		{
			auto mixture = std::make_shared<MixtureSetUp>("Mixture3");
			mixture->setPrimers("ama1_S0_Sub0_mip2", "CAATATTTAAAAGATGGAGGTTTTGCT",
					"ATCTTCACAATTTCCATCGACCCAT");
			mixture->addAbundance("PfML01", 10);
			mixture->addAbundance("PfDd2", 20);
			mixture->addAbundance("PfGB4-PfTG01", 30);
			mixture->addAbundance("PfKH02", 40);
			mixture->setForwardBarcode("MID10", "CGTCTAGTAC");
			mixture->setReverseBarcode("MID10", "TCTCTATGCG");
			mixture->forwardBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
						mixture->reverseBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
			sampleSet->addMixture(mixture);
		}
		{
			auto mixture = std::make_shared<MixtureSetUp>("Mixture4");
			mixture->setPrimers("ama1_S0_Sub0_mip2", "CAATATTTAAAAGATGGAGGTTTTGCT",
					"ATCTTCACAATTTCCATCGACCCAT");
			mixture->addAbundance("PfML01", 10);
			mixture->addAbundance("PfDd2", 20);
			mixture->addAbundance("PfGB4-PfTG01", 30);
			mixture->addAbundance("PfKH02", 40);
			mixture->setForwardBarcode("MID11", "TCTACGTAGC");
			mixture->setReverseBarcode("MID11", "TGATACGTCT");
			mixture->forwardBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
						mixture->reverseBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
			sampleSet->addMixture(mixture);
		}
		{
			auto mixture = std::make_shared<MixtureSetUp>("Mixture5");
			mixture->setPrimers("ama1_S0_Sub0_mip2", "CAATATTTAAAAGATGGAGGTTTTGCT",
					"ATCTTCACAATTTCCATCGACCCAT");
			mixture->addAbundance("PfML01", 2);
			mixture->addAbundance("PfDd2", 3);
			mixture->addAbundance("PfGB4-PfTG01", 5);
			mixture->addAbundance("PfKH02", 90);
			mixture->setForwardBarcode("MID01", "ACGAGTGCGT");
			mixture->setReverseBarcode("MID01", "GAGTGCGTCT");
			mixture->forwardBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
			sampleSet->addMixture(mixture);
		}

		{
			auto mixture = std::make_shared<MixtureSetUp>("Mixture6");
			mixture->setPrimers("ama1_S0_Sub0_mip2", "CAATATTTAAAAGATGGAGGTTTTGCT",
					"ATCTTCACAATTTCCATCGACCCAT");
			mixture->addAbundance("PfML01", 1);
			mixture->addAbundance("PfDd2", 1);
			mixture->addAbundance("PfGB4-PfTG01", 1);
			mixture->addAbundance("PfKH02", 99);
			mixture->setForwardBarcode("MID29", "ACTGTACAGT");
			mixture->setReverseBarcode("MID29", "CACGCTACGT");
			mixture->forwardBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
			sampleSet->addMixture(mixture);
		}
		lSetup.addSample(sampleSet);
	}
	{
		auto sampleSet = std::make_shared<SampleSetup>("Sample4-BarcodeBothEndsSame");
		{
			auto mixture = std::make_shared<MixtureSetUp>("Mixture1");
			mixture->setPrimers("ama1_S0_Sub0_mip2", "CAATATTTAAAAGATGGAGGTTTTGCT",
					"ATCTTCACAATTTCCATCGACCCAT");
			mixture->addAbundance("PfML01", 10);
			mixture->addAbundance("PfDd2", 10);
			mixture->addAbundance("PfGB4-PfTG01", 10);
			mixture->addAbundance("PfKH02", 10);
			mixture->setForwardBarcode("MID07", "CGTGTCTCTA");
			mixture->setReverseBarcode("MID07", "CGTGTCTCTA");
			mixture->forwardBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
						mixture->reverseBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
			sampleSet->addMixture(mixture);
		}
		{
			auto mixture = std::make_shared<MixtureSetUp>("Mixture2");
			mixture->setPrimers("ama1_S0_Sub0_mip2", "CAATATTTAAAAGATGGAGGTTTTGCT",
					"ATCTTCACAATTTCCATCGACCCAT");
			mixture->addAbundance("PfML01", 10);
			mixture->addAbundance("PfDd2", 20);
			mixture->addAbundance("PfGB4-PfTG01", 30);
			mixture->addAbundance("PfKH02", 40);
			mixture->setForwardBarcode("MID08", "CTCGCGTGTC");
			mixture->setReverseBarcode("MID08", "CTCGCGTGTC");
			mixture->forwardBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
						mixture->reverseBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
			sampleSet->addMixture(mixture);
		}
		{
			auto mixture = std::make_shared<MixtureSetUp>("Mixture3");
			mixture->setPrimers("ama1_S0_Sub0_mip2", "CAATATTTAAAAGATGGAGGTTTTGCT",
					"ATCTTCACAATTTCCATCGACCCAT");
			mixture->addAbundance("PfML01", 10);
			mixture->addAbundance("PfDd2", 20);
			mixture->addAbundance("PfGB4-PfTG01", 30);
			mixture->addAbundance("PfKH02", 40);
			mixture->setForwardBarcode("MID10", "TCTCTATGCG");
			mixture->setReverseBarcode("MID10", "TCTCTATGCG");
			mixture->forwardBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
						mixture->reverseBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
			sampleSet->addMixture(mixture);
		}
		{
			auto mixture = std::make_shared<MixtureSetUp>("Mixture4");
			mixture->setPrimers("ama1_S0_Sub0_mip2", "CAATATTTAAAAGATGGAGGTTTTGCT",
					"ATCTTCACAATTTCCATCGACCCAT");
			mixture->addAbundance("PfML01", 10);
			mixture->addAbundance("PfDd2", 20);
			mixture->addAbundance("PfGB4-PfTG01", 30);
			mixture->addAbundance("PfKH02", 40);
			mixture->setForwardBarcode("MID11", "TGATACGTCT");
			mixture->setReverseBarcode("MID11", "TGATACGTCT");
			mixture->forwardBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
						mixture->reverseBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
			sampleSet->addMixture(mixture);
		}
		{
			auto mixture = std::make_shared<MixtureSetUp>("Mixture5");
			mixture->setPrimers("ama1_S0_Sub0_mip2", "CAATATTTAAAAGATGGAGGTTTTGCT",
					"ATCTTCACAATTTCCATCGACCCAT");
			mixture->addAbundance("PfML01", 2);
			mixture->addAbundance("PfDd2", 3);
			mixture->addAbundance("PfGB4-PfTG01", 5);
			mixture->addAbundance("PfKH02", 90);
			mixture->setForwardBarcode("MID01", "ACGAGTGCGT");
			mixture->setReverseBarcode("MID01", "ACGAGTGCGT");
			mixture->forwardBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
			sampleSet->addMixture(mixture);
		}
		{
			auto mixture = std::make_shared<MixtureSetUp>("Mixture6");
			mixture->setPrimers("ama1_S0_Sub0_mip2", "CAATATTTAAAAGATGGAGGTTTTGCT",
					"ATCTTCACAATTTCCATCGACCCAT");
			mixture->addAbundance("PfML01", 1);
			mixture->addAbundance("PfDd2", 1);
			mixture->addAbundance("PfGB4-PfTG01", 1);
			mixture->addAbundance("PfKH02", 99);
			mixture->setForwardBarcode("MID29", "ACTGTACAGT");
			mixture->setReverseBarcode("MID29", "ACTGTACAGT");
			mixture->forwardBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
			sampleSet->addMixture(mixture);
		}
		lSetup.addSample(sampleSet);
	}

	{
		auto sampleSet = std::make_shared<SampleSetup>("Sample5-BarcodeBothEndsRComp");
		{
			auto mixture = std::make_shared<MixtureSetUp>("Mixture1");
			mixture->setPrimers("ama1_S0_Sub0_mip2", "CAATATTTAAAAGATGGAGGTTTTGCT",
					"ATCTTCACAATTTCCATCGACCCAT");
			mixture->addAbundance("PfML01", 10);
			mixture->addAbundance("PfDd2", 10);
			mixture->addAbundance("PfGB4-PfTG01", 10);
			mixture->addAbundance("PfKH02", 10);
			mixture->setForwardBarcode("MID07", "CGTGTCTCTA");
			mixture->setReverseBarcode("MID07", "TAGAGACACG");
			mixture->forwardBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
						mixture->reverseBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
			sampleSet->addMixture(mixture);
		}
		{
			auto mixture = std::make_shared<MixtureSetUp>("Mixture2");
			mixture->setPrimers("ama1_S0_Sub0_mip2", "CAATATTTAAAAGATGGAGGTTTTGCT",
					"ATCTTCACAATTTCCATCGACCCAT");
			mixture->addAbundance("PfML01", 10);
			mixture->addAbundance("PfDd2", 20);
			mixture->addAbundance("PfGB4-PfTG01", 30);
			mixture->addAbundance("PfKH02", 40);
			mixture->setForwardBarcode("MID08", "CTCGCGTGTC");
			mixture->setReverseBarcode("MID08", "GACACGCGAG");
			mixture->forwardBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
						mixture->reverseBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
			sampleSet->addMixture(mixture);
		}
		{
			auto mixture = std::make_shared<MixtureSetUp>("Mixture3");
			mixture->setPrimers("ama1_S0_Sub0_mip2", "CAATATTTAAAAGATGGAGGTTTTGCT",
					"ATCTTCACAATTTCCATCGACCCAT");
			mixture->addAbundance("PfML01", 10);
			mixture->addAbundance("PfDd2", 20);
			mixture->addAbundance("PfGB4-PfTG01", 30);
			mixture->addAbundance("PfKH02", 40);
			mixture->setForwardBarcode("MID10", "TCTCTATGCG");
			mixture->setReverseBarcode("MID10", "CGCATAGAGA");
			mixture->forwardBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
						mixture->reverseBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
			sampleSet->addMixture(mixture);
		}
		{
			auto mixture = std::make_shared<MixtureSetUp>("Mixture4");
			mixture->setPrimers("ama1_S0_Sub0_mip2", "CAATATTTAAAAGATGGAGGTTTTGCT",
					"ATCTTCACAATTTCCATCGACCCAT");
			mixture->addAbundance("PfML01", 10);
			mixture->addAbundance("PfDd2", 20);
			mixture->addAbundance("PfGB4-PfTG01", 30);
			mixture->addAbundance("PfKH02", 40);
			mixture->setForwardBarcode("MID11", "TGATACGTCT");
			mixture->setReverseBarcode("MID11", "AGACGTATCA");
			mixture->forwardBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
						mixture->reverseBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
			sampleSet->addMixture(mixture);
		}
		{
			auto mixture = std::make_shared<MixtureSetUp>("Mixture5");
			mixture->setPrimers("ama1_S0_Sub0_mip2", "CAATATTTAAAAGATGGAGGTTTTGCT",
					"ATCTTCACAATTTCCATCGACCCAT");
			mixture->addAbundance("PfML01", 2);
			mixture->addAbundance("PfDd2", 3);
			mixture->addAbundance("PfGB4-PfTG01", 5);
			mixture->addAbundance("PfKH02", 90);
			mixture->setForwardBarcode("MID01", "ACGAGTGCGT");
			mixture->setReverseBarcode("MID01", "ACGCACTCGT");
			mixture->forwardBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
			sampleSet->addMixture(mixture);
		}
		{
			auto mixture = std::make_shared<MixtureSetUp>("Mixture6");
			mixture->setPrimers("ama1_S0_Sub0_mip2", "CAATATTTAAAAGATGGAGGTTTTGCT",
					"ATCTTCACAATTTCCATCGACCCAT");
			mixture->addAbundance("PfML01", 1);
			mixture->addAbundance("PfDd2", 1);
			mixture->addAbundance("PfGB4-PfTG01", 1);
			mixture->addAbundance("PfKH02", 99);
			mixture->setForwardBarcode("MID29", "ACTGTACAGT");
			mixture->setReverseBarcode("MID29", "ACTGTACAGT");
			mixture->forwardBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
			sampleSet->addMixture(mixture);
		}

		lSetup.addSample(sampleSet);
	}
	out << lSetup.toJson() << std::endl;

//	PfML01 vs PfDd2 1
//	PfGB4-PfTG01 vs PfKH02 3
//	PfGA01 vs Pf7G8 5
	return 0;
}

//
//
//int readSimulatorRunner::simMultipleMixture(const njh::progutils::CmdArgs & inputCommands) {
//
//  readSimulatorSetUp setUp(inputCommands);
//
//	uint32_t numThreads = 2;
//	bool singleEnd = false;
//	bool nonGz = false;
//	bfs::path idFile = "";
//	bfs::path referenceFile = "";
//	bfs::path librarySetUpFile = "";
//
//	bfs::path illuminaProfileDir = "";
//	setUp.processVerbose();
//	setUp.processDebug();
//	setUp.setOption(librarySetUpFile, "--librarySetUpFile", "Library Set Up File", true);
//	setUp.processReadInNames(VecStr{"--fasta", "--fastagz"},true);
//	setUp.setOption(numThreads, "--numThreads", "Number of Threads to Use");
//	setUp.setOption(singleEnd, "--singleEnd", "Single End");
//	setUp.setOption(nonGz, "--nonGz", "do not compress the output fastqs");
//	setUp.setOption(illuminaProfileDir, "--illuminaProfileDir", "Illumina Profile Dir", true);
//	setUp.processDirectoryOutputName("simMultipleMixture_TODAY", true);
//	setUp.finishSetUp(std::cout);
//	setUp.startARunLog(setUp.pars_.directoryName_);
//	setUp.writeParametersFile(setUp.pars_.directoryName_ + "parameters.tab.txt", false, true);
//
//	std::unordered_map<std::string, std::shared_ptr<seqInfo>> refSeqs;
//
//
//	{
//		SeqInput reader(setUp.pars_.ioOptions_);
//		reader.openIn();
//		seqInfo seq;
//		while(reader.readNextRead(seq)){
//			refSeqs[seq.name_] = std::make_shared<seqInfo>(seq);
//		}
//	}
//	Json::Value librarySetUp = njh::json::parseFile(librarySetUpFile.string());
//
//	njh::json::MemberChecker libraryChecker(librarySetUp);
//	libraryChecker.failMemberCheckThrow(LibrarySetup::jsonMembers(), __PRETTY_FUNCTION__);
//	LibrarySetup::SimLibrarySetupPars simPars(librarySetUp["pars"]);
//	LibrarySetup lSetup(librarySetUp["name"].asString(), simPars, librarySetUp, getVectorOfMapKeys(refSeqs));
//
//
//	OutOptions idOutOpts(njh::files::make_path(setUp.pars_.directoryName_, "ids.tab.txt"));
//
//	lSetup.ids_->writeIdFile(OutOptions(njh::files::make_path(setUp.pars_.directoryName_, "ids.tab.txt")));
//
//	OutOptions outOptsSetup(njh::files::make_path(setUp.pars_.directoryName_, "setup.json"));
//	OutputStream outOptsSetupOut(outOptsSetup);
//	outOptsSetupOut << lSetup.toJson() << std::endl;
//
//	bfs::path fastqDirectory =njh::files::make_path(setUp.pars_.directoryName_, "fastq");
//	njh::files::makeDir(njh::files::MkdirPar{fastqDirectory});
//	auto sampleKeys = njh::getVecOfMapKeys(lSetup.samples_);
//	njh::sort(sampleKeys);
//
//	std::atomic_uint sampleCountAtom{0};
//
//	njh::concurrent::LockableQueue<std::string> samplesQueue(sampleKeys);
//	std::string extension = nonGz ? ".fastq" : ".fastq.gz";
//
//	auto simSample = [&samplesQueue,&extension,&sampleCountAtom,&fastqDirectory,&lSetup,&illuminaProfileDir,&singleEnd,&refSeqs](){
//		std::string sampleKey = "";
//		njh::randObjectGen<char,uint32_t> baseRGen({'A', 'C', 'G', 'T'}, {1,1,1,1});
//		RoughIlluminaSimulator simulator(illuminaProfileDir);
//		njh::randomGenerator rGen;
//
//		while(samplesQueue.getVal(sampleKey)){
//			uint32_t sampleCount = sampleCountAtom++;
//			OutOptions targetOutOpts(njh::files::make_path(fastqDirectory, njh::pasteAsStr(sampleKey, "_S", sampleCount + 1, extension)));
//			OutOptions r1OutOpts(njh::files::make_path(fastqDirectory, njh::pasteAsStr(sampleKey, "_S", sampleCount + 1, "_R1_001" + extension)));
//			OutOptions r2OutOpts(njh::files::make_path(fastqDirectory, njh::pasteAsStr(sampleKey, "_S", sampleCount + 1, "_R2_001" + extension)));
//			std::shared_ptr<OutputStream> targetOut;
//
//			std::shared_ptr<OutputStream> r1Out;
//			std::shared_ptr<OutputStream> r2Out;
//			if(singleEnd) {
//				targetOut = std::make_shared<OutputStream>(targetOutOpts);
//			} else {
//				r1Out = std::make_shared<OutputStream>(r1OutOpts);
//				r2Out = std::make_shared<OutputStream>(r2OutOpts);
//			}
//
//
//			/**@todo add back in the ability to have different Illumina barcodes in overhangs */
//	//		std::string sampleAdapter1 = defaultAdapter1;
//	//		std::string sampleAdapter2 = defaultAdapter2;
//	//		for(const auto pos : iter::range(sampleAdapter1.size())){
//	//			if('N' == sampleAdapter1[pos]){
//	//				sampleAdapter1[pos] = baseRGen.genObj();
//	//			}
//	//		}
//	//		for(const auto pos : iter::range(sampleAdapter2.size())){
//	//			if('N' == sampleAdapter2[pos]){
//	//				sampleAdapter2[pos] = baseRGen.genObj();
//	//			}
//	//		}
//			for(const auto & mixture : lSetup.samples_.at(sampleKey)->mixtures_ ){
//				VecStr refNames;
//				std::vector<double> amounts;
//				for(const auto & abun : mixture.second->expectedAbundances_){
//					refNames.emplace_back(abun.first);
//					amounts.emplace_back(abun.second);
//				}
//				njh::randObjectGen<std::string,double> refNameGen(refNames, amounts);
//				for(uint32_t readNumber = 0; readNumber < mixture.second->finalReadAmount_; ++readNumber){
//				//for(uint32_t readNumber = 0; readNumber < lSetup.pars_.sampleReadAmount_; ++readNumber){
//					MetaDataInName nameMeta;
//					auto refName = refNameGen.genObj();
//
//					nameMeta.addMeta("Mixture", mixture.second->name_);
//					nameMeta.addMeta("Sample", sampleKey);
//					nameMeta.addMeta("readNumber", readNumber);
//					nameMeta.addMeta("refName", refName);
//					seqInfo targetSeq("", refSeqs.at(refName)->seq_);
//
//					uint32_t forwardPrimerPosition = 0;
//					uint32_t reversePrimerPosition = 0;
//
//					uint32_t forwardBarcodePosition = 0;
//					uint32_t reverseBarcodePosition = 0;
//
//					if(nullptr != mixture.second->primers_){
////						if(!njh::in(mixture.second->primers_->name_, primerNames)){
////							primerOut << mixture.second->primers_->name_
////									<< "\t" << mixture.second->primers_->forward_
////									<< "\t" << mixture.second->primers_->reverse_
////									<< std::endl;
////							primerNames.emplace_back(mixture.second->primers_->name_);
////						}
//						if(!lSetup.pars_.noAddPrimers_){
//							/**@todo add allowing for ambigious bases in primers*/
//							targetSeq.prepend(mixture.second->primers_->forward_);
//							reversePrimerPosition = targetSeq.seq_.size();
//							targetSeq.append(mixture.second->primers_->reverse_3_5_);
//						}
//
//
//
//						if(mixture.second->primers_->forward_randomPrecedingBases_ > 0){
//							auto minBaseAmount = mixture.second->primers_->randomOneLength_ ? mixture.second->primers_->forward_randomPrecedingBases_: 0;
//							uint32_t numFBase = rGen.unifRand<uint32_t>(minBaseAmount, mixture.second->primers_->forward_randomPrecedingBases_ + 1);
//							if(numFBase > 0){
//								targetSeq.prepend(njh::pasteAsStr(baseRGen.genObjs(numFBase)));
//								reversePrimerPosition += numFBase;
//								forwardPrimerPosition += numFBase;
//							}
//						}
//						if(mixture.second->primers_->reverse_randomPrecedingBases_ > 0){
//							auto minBaseAmount = mixture.second->primers_->randomOneLength_ ? mixture.second->primers_->reverse_randomPrecedingBases_: 0;
//							uint32_t numRBase = rGen.unifRand<uint32_t>(minBaseAmount, mixture.second->primers_->reverse_randomPrecedingBases_ + 1);
//							if(numRBase > 0){
//								targetSeq.append(njh::pasteAsStr(baseRGen.genObjs(numRBase)));
//							}
//						}
//					}
//					if(nullptr != mixture.second->forwardBarcode_ || nullptr != mixture.second->reverseBarcode_){
//						std::string midName = nullptr != mixture.second->forwardBarcode_ ? mixture.second->forwardBarcode_->name_ : mixture.second->reverseBarcode_->name_;
////						if(!njh::in(midName, midNames)){
////							midNames.emplace_back(midName);
////							midOut << midName;
////							if(nullptr != mixture.second->forwardBarcode_){
////								midOut << "\t" << mixture.second->forwardBarcode_->barcode_;
////							}else{
////								midOut << "\t";
////							}
////							if(nullptr != mixture.second->reverseBarcode_){
////								midOut << "\t" << mixture.second->reverseBarcode_->barcode_;
////							}
////							midOut << std::endl;
////						}
//					}
//					if(nullptr != mixture.second->forwardBarcode_){
//						targetSeq.prepend(mixture.second->forwardBarcode_->barcode_);
//						reversePrimerPosition += mixture.second->forwardBarcode_->barcode_.size();
//						forwardPrimerPosition += mixture.second->forwardBarcode_->barcode_.size();
//						if(mixture.second->forwardBarcode_->randomPrecedingBases_ > 0){
//							auto minBaseAmount = mixture.second->forwardBarcode_->randomOneLength_ ? mixture.second->forwardBarcode_->randomPrecedingBases_: 0;
//							uint32_t numFBase = rGen.unifRand<uint32_t>(minBaseAmount, mixture.second->forwardBarcode_->randomPrecedingBases_ + 1);
//							if(numFBase > 0){
//								targetSeq.prepend(njh::pasteAsStr(baseRGen.genObjs(numFBase)));
//								reversePrimerPosition += numFBase;
//								forwardPrimerPosition += numFBase;
//								forwardBarcodePosition+= numFBase;
//							}
//						}
//					}
//					if(nullptr != mixture.second->reverseBarcode_){
//						reverseBarcodePosition = len(targetSeq);
//						targetSeq.append(mixture.second->reverseBarcode_->barcode_3_5_);
//						if(mixture.second->reverseBarcode_->randomPrecedingBases_ > 0){
//							auto minBaseAmount = mixture.second->reverseBarcode_->randomOneLength_ ? mixture.second->reverseBarcode_->randomPrecedingBases_: 0;
//							uint32_t numFBase = rGen.unifRand<uint32_t>(minBaseAmount, mixture.second->reverseBarcode_->randomPrecedingBases_ + 1);
//							if(numFBase > 0){
//								targetSeq.append(njh::pasteAsStr(baseRGen.genObjs(numFBase)));
//							}
//						}
//					}
//					//blunt ending
//					if(lSetup.pars_.addBluntEndingArtifact_){
//						if('A' == targetSeq.seq_[0] && rGen() < lSetup.pars_.bluntEndingArtifactChance_){
//							readVecTrimmer::trimOffForwardBases(targetSeq, 1);
//							forwardPrimerPosition = forwardPrimerPosition == 0 ? forwardPrimerPosition: forwardPrimerPosition - 1;
//							forwardBarcodePosition = forwardBarcodePosition == 0 ? forwardBarcodePosition: forwardBarcodePosition - 1;
//							reverseBarcodePosition -= 1;
//							reversePrimerPosition -= 1;
//							nameMeta.addMeta("frontEndBluntEndArtifact", true);
//						}else{
//							nameMeta.addMeta("frontEndBluntEndArtifact", false);
//						}
//						if ('T' == targetSeq.seq_.back()
//								&& rGen() < lSetup.pars_.bluntEndingArtifactChance_) {
//							readVecTrimmer::trimOffEndBases(targetSeq, 1);
//							nameMeta.addMeta("backEndBluntEndArtifact", true);
//						} else {
//							nameMeta.addMeta("backEndBluntEndArtifact", false);
//						}
//					}
//					//positions
//					if(nullptr != mixture.second->primers_){
//						nameMeta.addMeta("forwardPrimerPosition", forwardPrimerPosition);
//						nameMeta.addMeta("reversePrimerPosition", reversePrimerPosition);
//					}
//					if(nullptr != mixture.second->forwardBarcode_){
//						nameMeta.addMeta("forwardBarcodePosition", forwardBarcodePosition);
//					}
//					if(nullptr != mixture.second->reverseBarcode_){
//						nameMeta.addMeta("reverseBarcodePosition", reverseBarcodePosition);
//					}
//					//complement
//					if(lSetup.pars_.addReverseComplement_){
//						bool complement = rGen.unifRand(0,2) == 0;
//						nameMeta.addMeta("complement", complement);
//						if(complement){
//							targetSeq.reverseComplementRead(false, true);
//						}
//					}
//					//length
//					nameMeta.addMeta("targetSeqLength", len(targetSeq));
//
//					targetSeq.name_ = nameMeta.createMetaName();
//
//					if(singleEnd){
//						//target
//						targetSeq.outPutSeq(*targetOut);
//					}else{
//						{
//							//r1
//							auto subSeq = len(targetSeq) > lSetup.pars_.pairedEndLength_? targetSeq.getSubRead(0, lSetup.pars_.pairedEndLength_) : targetSeq;
//							simulator.simR1(subSeq, lSetup.pars_.pairedEndLength_).outPutFastq(*r1Out);
//						}
//						{
//							//r2
//							targetSeq.reverseComplementRead(false, true);
//							auto subSeq = len(targetSeq) > lSetup.pars_.pairedEndLength_? targetSeq.getSubRead(0, lSetup.pars_.pairedEndLength_) : targetSeq;
//							simulator.simR2(subSeq, lSetup.pars_.pairedEndLength_).outPutFastq(*r2Out);
//						}
//					}
//				}
//			}
//		}
//	};
//
//	std::vector<std::thread> threads;
//	for(uint32_t t = 0; t < numThreads; ++t){
//		threads.emplace_back(std::thread(simSample));
//	}
//	njh::concurrent::joinAllJoinableThreads(threads);
//
//	return 0;
//}
//





int readSimulatorRunner::simMultipleMixtureSimPCR(const njh::progutils::CmdArgs & inputCommands) {
  readSimulatorSetUp setUp(inputCommands);
	uint32_t numThreads = 2;
	uint32_t pcrNumThreads = 2;
	bool singleEnd = false;
	bool nonGz = false;
	bfs::path idFile = "";
	bfs::path referenceFile = "";
	bfs::path librarySetUpFile = "";
	bfs::path illuminaProfileDir = "";
	uint32_t pcrRounds = 30;
	uint32_t initialPcrRounds = 10;
	std::map<uint32_t, uint32_t> initialPcrRoundsMap;
	bfs::path initialPcrRoundsMapFnp;
	long double errorRate = 3.5e-06;
	double pcrEfficiency = 0.95;
	bool keepPCRSeqs = false;
	uint32_t chimeraBasesIn = 5;
	uint32_t templateCap = 500000000;
	bool noChimeras = false;
	uint32_t pairedEndLength = std::numeric_limits<uint32_t>::max();
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(pairedEndLength, "--pairedEndLength", "Paired End Length");

	setUp.setOption(noChimeras, "--noChimeras", "Don't simulate chimeras");
	setUp.setOption(templateCap, "--templateCap", "Template Cap");
	setUp.setOption(chimeraBasesIn, "--chimeraBasesIn", "The number of bases needed for a template to lay down");
	setUp.setOption(keepPCRSeqs, "--keepPCRSeqs", "Keep PCR Seqs");
	setUp.setOption(errorRate, "--errorRate", "Polymerase Error Rate");
	setUp.setOption(pcrEfficiency, "--pcrEfficiency", "PCR Efficiency, between 0-1, chance a product gets amplified");
	setUp.setOption(pcrRounds, "--pcrRounds", "Number of PCR rounds");
	setUp.setOption(initialPcrRounds, "--initialPcrRounds", "Number of Initial rounds of PCR before sampling");
	setUp.setOption(initialPcrRoundsMapFnp, "--initialPcrRoundsTable", "Number of Initial rounds of PCR before sampling per starting template amount, columns 1)template, 2) rounds");

	setUp.setOption(librarySetUpFile, "--librarySetUpFile", "Library Set Up File", true);
	setUp.processReadInNames(VecStr{"--fasta", "--fastagz"},true);
	setUp.setOption(numThreads, "--numThreads", "Number of Threads to Use");
	setUp.setOption(pcrNumThreads, "--pcrNumThreads", "Number of Threads to Use for PCR sim");
	setUp.setOption(singleEnd, "--singleEnd", "Single End");
	setUp.setOption(nonGz, "--nonGz", "do not compress the output fastqs");
	setUp.setOption(illuminaProfileDir, "--illuminaProfileDir", "Illumina Profile Dir", true);
	setUp.processDirectoryOutputName("simMultipleMixture_TODAY", true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	initialPcrRoundsMap[1] = initialPcrRounds;
	if(bfs::exists(initialPcrRoundsMapFnp)){
		initialPcrRoundsMap.clear();
		table initialPcrRoundsMapTab(initialPcrRoundsMapFnp, "\t", true);
		initialPcrRoundsMapTab.checkForColumnsThrow(VecStr{"template", "rounds"}, __PRETTY_FUNCTION__);
		for(const auto & row : initialPcrRoundsMapTab){
			initialPcrRoundsMap[njh::StrToNumConverter::stoToNum<uint32_t>(row[initialPcrRoundsMapTab.getColPos("template")])] =njh::StrToNumConverter::stoToNum<uint32_t>(row[initialPcrRoundsMapTab.getColPos("rounds")]);
		}
	}
	uint64_t intErrorRate = errorRate * std::numeric_limits<uint64_t>::max();

	setUp.writeParametersFile(setUp.pars_.directoryName_ + "parameters.tab.txt", false, true);
	std::unordered_map<std::string, std::shared_ptr<seqInfo>> refSeqs;
	std::unordered_map<std::string, std::vector<std::shared_ptr<seqInfo>>> refSeqsByPrimerName;
	std::unordered_map<std::string, std::set<std::string>> refSeqsNamesByPrimerName;
	uint64_t maxLen = 0;
	{
		SeqInput reader(setUp.pars_.ioOptions_);
		reader.openIn();
		seqInfo seq;
		while(reader.readNextRead(seq)){
			readVec::getMaxLength(seq, maxLen);
			refSeqs[seq.name_] = std::make_shared<seqInfo>(seq);
		}
	}
	Json::Value librarySetUp = njh::json::parseFile(librarySetUpFile.string());
	njh::json::MemberChecker libraryChecker(librarySetUp);
	libraryChecker.failMemberCheckThrow(LibrarySetup::jsonMembers(), __PRETTY_FUNCTION__);
	LibrarySetup::SimLibrarySetupPars simPars(librarySetUp["pars"]);
	if(std::numeric_limits<uint32_t>::max() != pairedEndLength){
		simPars.pairedEndLength_ = pairedEndLength;
	}
	LibrarySetup lSetup(librarySetUp["name"].asString(), simPars, librarySetUp, getVectorOfMapKeys(refSeqs));

	OutOptions idOutOpts(njh::files::make_path(setUp.pars_.directoryName_, "ids.tab.txt"));
	lSetup.ids_->writeIdFile(OutOptions(njh::files::make_path(setUp.pars_.directoryName_, "ids.tab.txt")));
	for(const auto & samp : lSetup.samples_){
		for(const auto & mix : samp.second->mixtures_){
			std::string primerName = "NoPrimers";
			if(nullptr != mix.second->primers_){
				primerName = mix.second->primers_->name_;
			}
			for(const auto & hap : mix.second->expectedAbundances_){
				refSeqsNamesByPrimerName[primerName].emplace(hap.first);
			}
		}
	}
	for(const auto & refSeqNames : refSeqsNamesByPrimerName){
		for(const auto & refSeqName : refSeqNames.second){
			refSeqsByPrimerName[refSeqNames.first].push_back(njh::mapAt(refSeqs, refSeqName));
		}
	}
	comparison allowableErrors;
	allowableErrors.hqMismatches_ = 3;
	FullTrimReadsPars::trimSeqPars trimPars;
	trimPars.alwaysTrim = true;
	trimPars.includeSequence_ = false;

	bfs::path refSeqsDir = njh::files::make_path(setUp.pars_.directoryName_, "refSeqs");
	lSetup.ids_->initPrimerDeterminator();
	njh::files::makeDir(njh::files::MkdirPar{refSeqsDir});
	aligner alignerObj(maxLen, gapScoringParameters(5,1,0,0,0,0));
	for(const auto & refSeqs : refSeqsByPrimerName){
		SeqOutput seqWriter(SeqIOOptions::genFastaOut(njh::files::make_path(refSeqsDir, refSeqs.first)));
		if("NoPrimers" == refSeqs.first){
			for(const auto & refSeq : refSeqs.second){
				seqWriter.openWrite(refSeq);
			}
		}else{
			for(const auto & refSeq : refSeqs.second){
				seqInfo copySeqToTrim = *refSeq;
				readVecTrimmer::trimBetweenSequences(copySeqToTrim,
						lSetup.ids_->pDeterminator_->primers_.at(refSeqs.first).forwardPrimerInfo_,
						lSetup.ids_->pDeterminator_->primers_.at(refSeqs.first).reversePrimerInfo_,
						alignerObj, allowableErrors, trimPars);
				seqWriter.openWrite(copySeqToTrim);
			}
		}
	}

	OutOptions outOptsSetup(njh::files::make_path(setUp.pars_.directoryName_, "setup.json"));
	OutputStream outOptsSetupOut(outOptsSetup);
	outOptsSetupOut << lSetup.toJson() << std::endl;
	bfs::path fastqDirectory =njh::files::make_path(setUp.pars_.directoryName_, "fastq");
	bfs::path chimeraDirectory =njh::files::make_path(setUp.pars_.directoryName_, "chimeras");
	njh::files::makeDir(njh::files::MkdirPar{fastqDirectory});
	if(!noChimeras){
		njh::files::makeDir(njh::files::MkdirPar{chimeraDirectory});
	}
	auto sampleKeys = njh::getVecOfMapKeys(lSetup.samples_);
	njh::sort(sampleKeys);
	std::atomic_uint sampleCountAtom{0};
	njh::concurrent::LockableQueue<std::string> samplesQueue(sampleKeys);
	std::string extension = nonGz ? ".fastq" : ".fastq.gz";
	bool verbose = setUp.pars_.verbose_;
	auto simAmountsDir = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar{"simulation_sampling_info"});

	auto simSample = [&samplesQueue,&extension,&sampleCountAtom,&fastqDirectory,
										&lSetup,&illuminaProfileDir,&singleEnd,&refSeqs,
										&intErrorRate, &pcrRounds, &initialPcrRounds,&verbose,
										&keepPCRSeqs, &pcrNumThreads,&pcrEfficiency,&simAmountsDir,
										&chimeraDirectory,&templateCap,&noChimeras,&chimeraBasesIn,&initialPcrRoundsMap](){
		std::string sampleKey = "";
		njh::randObjectGen<char,uint32_t> baseRGen({'A', 'C', 'G', 'T'}, {1,1,1,1});
		RoughIlluminaSimulator simulator(illuminaProfileDir);
		njh::randomGenerator rGen;
		while(samplesQueue.getVal(sampleKey)){
			uint32_t sampleCount = sampleCountAtom++;
			OutOptions targetOutOpts(njh::files::make_path(fastqDirectory, njh::pasteAsStr(sampleKey, "_S", sampleCount + 1, extension)));
			OutOptions r1OutOpts(njh::files::make_path(fastqDirectory, njh::pasteAsStr(sampleKey, "_S", sampleCount + 1, "_R1_001" + extension)));
			OutOptions r2OutOpts(njh::files::make_path(fastqDirectory, njh::pasteAsStr(sampleKey, "_S", sampleCount + 1, "_R2_001" + extension)));
			std::shared_ptr<OutputStream> targetOut;
			std::shared_ptr<OutputStream> r1Out;
			std::shared_ptr<OutputStream> r2Out;
			if(singleEnd) {
				targetOut = std::make_shared<OutputStream>(targetOutOpts);
			} else {
				r1Out = std::make_shared<OutputStream>(r1OutOpts);
				r2Out = std::make_shared<OutputStream>(r2OutOpts);
			}

			/**@todo add back in the ability to have different Illumina barcodes in overhangs */
	//		std::string sampleAdapter1 = defaultAdapter1;
	//		std::string sampleAdapter2 = defaultAdapter2;
	//		for(const auto pos : iter::range(sampleAdapter1.size())){
	//			if('N' == sampleAdapter1[pos]){
	//				sampleAdapter1[pos] = baseRGen.genObj();
	//			}
	//		}
	//		for(const auto pos : iter::range(sampleAdapter2.size())){
	//			if('N' == sampleAdapter2[pos]){
	//				sampleAdapter2[pos] = baseRGen.genObj();
	//			}
	//		}
			for(const auto & mixture : lSetup.samples_.at(sampleKey)->mixtures_ ){
				VecStr refNames;
				std::vector<double> amounts;
				double totalAbundance = 0;
				std::vector<seqInfo> currentSeqs;
				for(const auto & abun : mixture.second->expectedAbundances_){
					totalAbundance+=abun.second;
				}
				for(const auto & abun : mixture.second->expectedAbundances_){
					std::string seqName = abun.first;
					std::string seq = refSeqs.at(abun.first)->seq_;
					/**@todo add back in randomizing this per seq*/
					if(nullptr != mixture.second->primers_){
						if(!lSetup.pars_.noAddPrimers_){
							/**@todo add allowing for ambigious bases in primers*/
							seq = mixture.second->primers_->forward_ + seq;
							seq.append(mixture.second->primers_->reverse_3_5_);
						}
						if(mixture.second->primers_->forward_randomPrecedingBases_ > 0){
							auto minBaseAmount = mixture.second->primers_->randomOneLength_ ? mixture.second->primers_->forward_randomPrecedingBases_: 0;
							uint32_t numFBase = rGen.unifRand<uint32_t>(minBaseAmount, mixture.second->primers_->forward_randomPrecedingBases_ + 1);
							if(numFBase > 0){
								seq = njh::pasteAsStr(baseRGen.genObjs(numFBase)) + seq;
							}
						}
						if(mixture.second->primers_->reverse_randomPrecedingBases_ > 0){
							auto minBaseAmount = mixture.second->primers_->randomOneLength_ ? mixture.second->primers_->reverse_randomPrecedingBases_: 0;
							uint32_t numRBase = rGen.unifRand<uint32_t>(minBaseAmount, mixture.second->primers_->reverse_randomPrecedingBases_ + 1);
							if(numRBase > 0){
								seq.append(njh::pasteAsStr(baseRGen.genObjs(numRBase)));
							}
						}
					}
					if(nullptr != mixture.second->forwardBarcode_ || nullptr != mixture.second->reverseBarcode_){
						std::string midName = nullptr != mixture.second->forwardBarcode_ ? mixture.second->forwardBarcode_->name_ : mixture.second->reverseBarcode_->name_;
					}
					if(nullptr != mixture.second->forwardBarcode_){
						seq =mixture.second->forwardBarcode_->barcode_ + seq;
						if(mixture.second->forwardBarcode_->randomPrecedingBases_ > 0){
							auto minBaseAmount = mixture.second->forwardBarcode_->randomOneLength_ ? mixture.second->forwardBarcode_->randomPrecedingBases_: 0;
							uint32_t numFBase = rGen.unifRand<uint32_t>(minBaseAmount, mixture.second->forwardBarcode_->randomPrecedingBases_ + 1);
							if(numFBase > 0){
								seq = njh::pasteAsStr(baseRGen.genObjs(numFBase)) + seq;
							}
						}
					}
					if(nullptr != mixture.second->reverseBarcode_){
						seq.append(mixture.second->reverseBarcode_->barcode_3_5_);
						if(mixture.second->reverseBarcode_->randomPrecedingBases_ > 0){
							auto minBaseAmount = mixture.second->reverseBarcode_->randomOneLength_ ? mixture.second->reverseBarcode_->randomPrecedingBases_: 0;
							uint32_t numFBase = rGen.unifRand<uint32_t>(minBaseAmount, mixture.second->reverseBarcode_->randomPrecedingBases_ + 1);
							if(numFBase > 0){
								seq.append(njh::pasteAsStr(baseRGen.genObjs(numFBase)));
							}
						}
					}
					currentSeqs.emplace_back(seqInfo{seqName, seq});
					currentSeqs.back().frac_ = abun.second/totalAbundance;
					currentSeqs.back().cnt_ = abun.second/totalAbundance;
				}
//				std::cout << sampleKey + mixture.second->name_ << std::endl;
//				std::cout << "\tmixture.second->startingTemplateAmount_: " << mixture.second->startingTemplateAmount_ << std::endl;
				PCRSimulator pcrSim(intErrorRate);
				pcrSim.verbose_ = verbose;
				pcrSim.pcrEfficiency_ = pcrEfficiency;
				pcrSim.templateCap_ = templateCap;
				pcrSim.noChimeras_ = noChimeras;
				pcrSim.chimeraPad_ = chimeraBasesIn;
				auto libDirName = njh::files::makeDir(njh::files::make_path(fastqDirectory).string(), njh::files::MkdirPar(sampleKey + mixture.second->name_));
				auto pcrReadsFnp = njh::files::make_path(fastqDirectory, sampleKey + mixture.second->name_, "reads.fasta");
				OutOptions pcrReadsOpts(pcrReadsFnp);
				std::vector<PCRSimulator::SeqGenomeCnt> seqGCounts;
				uint32_t totalGenomes = 0;
				for(const auto & currentSeq : currentSeqs){
					if(!njh::in(currentSeq.name_, mixture.second->genomeCounts_)){
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ", error " << " need to supply genome counts " << " none fond for " << currentSeq.name_ << "\n";
						ss << "Found: " << njh::conToStrEndSpecial(njh::getVecOfMapKeys(mixture.second->genomeCounts_), ", ", " and ") << "\n";
						throw std::runtime_error{ss.str()};
					}
					if(mixture.second->genomeCounts_.at(currentSeq.name_) > 0){
						totalGenomes += mixture.second->genomeCounts_.at(currentSeq.name_);
						seqGCounts.emplace_back(currentSeq, mixture.second->genomeCounts_.at(currentSeq.name_));
					}
				}
				uint32_t currentInitial = initialPcrRounds;
				for(const auto & initialRound : initialPcrRoundsMap){
					if(initialRound.first >= totalGenomes){
						currentInitial = initialRound.second;
						break;
					}
				}
				auto pcrSimAmounts = pcrSim.simLibFast(seqGCounts, pcrReadsOpts,mixture.second->finalReadAmount_, pcrRounds, currentInitial, pcrNumThreads);
				//PrimerPair	experiment	runNumber	sampName_
				std::string PrimerPair = mixture.second->meta_->getMeta("PrimerPair");
				std::string experiment = mixture.second->meta_->getMeta("experiment");
				std::string sampName   = mixture.second->meta_->getMeta("sampName_");
				uint32_t runNumber     = mixture.second->meta_->getMeta<uint32_t>("runNumber");
				for(const auto & sequenced : pcrSimAmounts.sampledForSequencing_){
					njh::mapAt(njh::mapAt(njh::mapAt(lSetup.pop_.samples_, sampName).expRuns_,experiment)[runNumber].hapRegionAmplified_, PrimerPair).hapAbundSequenced_[sequenced.first] = sequenced.second.mutated_ + sequenced.second.nonMutated_;
					//lSetup.pop_.samples_.at(sampName).expRuns_.at(experiment)[runNumber].hapRegionAmplified_.at(PrimerPair).hapAbundSequenced_.at(sequenced.first) = sequenced.second.mutated_ + sequenced.second.nonMutated_;
				}
				OutputStream simAmountsOut(njh::files::make_path(simAmountsDir, sampleKey + mixture.second->name_ + "_simAmounts.tab.txt"));
				simAmountsOut << "Sample\tMixture\tHapName\tExpectedAbund\tGenomesSampled\tFinalSampledMutated\tFinalSampledTotal";
				VecStr metalevels;
				if(nullptr != mixture.second->meta_){
					metalevels = getVectorOfMapKeys(mixture.second->meta_->meta_);
					njh::sort(metalevels);
					for(const auto & m : metalevels){
						simAmountsOut << "\t" << m;
					}
				}
				simAmountsOut << std::endl;
				uint32_t totalNonChimera = 0;
				for(const auto & seq : currentSeqs){
					simAmountsOut << sampleKey
							<< "\t" << mixture.second->name_
							<< "\t" << seq.name_
							<< "\t" << seq.frac_
							<< "\t" << pcrSimAmounts.genomesSampled_[seq.name_]
							<< "\t" << pcrSimAmounts.sampledForSequencing_[seq.name_].mutated_
							<< "\t" << pcrSimAmounts.sampledForSequencing_[seq.name_].nonMutated_ + pcrSimAmounts.sampledForSequencing_[seq.name_].mutated_;
					totalNonChimera+=  pcrSimAmounts.sampledForSequencing_[seq.name_].nonMutated_ + pcrSimAmounts.sampledForSequencing_[seq.name_].mutated_;
					if(nullptr != mixture.second->meta_){
						for(const auto & m : metalevels){
							simAmountsOut << "\t" << mixture.second->meta_->getMeta(m);
						}
					}
					simAmountsOut << std::endl;
				}
				if(!pcrSimAmounts.chimeraSeqs_.empty()){
					SeqOutput::write(pcrSimAmounts.chimeraSeqs_, SeqIOOptions::genFastaOut(njh::files::make_path(chimeraDirectory, sampleKey + mixture.second->name_ + "_chimerasGenerated.fasta")));
					OutputStream simAmountsChimeraOut(njh::files::make_path(simAmountsDir, sampleKey + mixture.second->name_ + "_chimeraSimAmounts.tab.txt"));
					simAmountsChimeraOut << "Sample\tMixture\tHapName\tFinalSampledMutated\tFinalSampledTotal";
					if(nullptr != mixture.second->meta_){
						for(const auto & m : metalevels){
							simAmountsChimeraOut << "\t" << m;
						}
					}
					simAmountsChimeraOut << std::endl;
					uint64_t totalChimera = 0;
					for(const auto & seq : pcrSimAmounts.chimeraSeqs_){
						simAmountsChimeraOut << sampleKey
								<< "\t" << mixture.second->name_
								<< "\t" << seq.name_
								<< "\t" << pcrSimAmounts.chimerasSampledForSequencing_[seq.name_].mutated_
								<< "\t" << pcrSimAmounts.chimerasSampledForSequencing_[seq.name_].nonMutated_ + pcrSimAmounts.chimerasSampledForSequencing_[seq.name_].mutated_;
						if(nullptr != mixture.second->meta_){
							for(const auto & m : metalevels){
								simAmountsChimeraOut << "\t" << mixture.second->meta_->getMeta(m);
							}
						}
						simAmountsChimeraOut << std::endl;
						totalChimera += pcrSimAmounts.chimerasSampledForSequencing_[seq.name_].nonMutated_ + pcrSimAmounts.chimerasSampledForSequencing_[seq.name_].mutated_;
					}
					OutputStream simAmountsChimeraTotCountOut(njh::files::make_path(simAmountsDir, sampleKey + mixture.second->name_ + "_chimeraSimFinals.tab.txt"));//
					simAmountsChimeraTotCountOut << "Sample\tMixture\tTotal\tchimera\tfullTemplate\tpercentChimera";
					if(nullptr != mixture.second->meta_){
						for(const auto & m : metalevels){
							simAmountsChimeraTotCountOut << "\t" << m;
						}
					}
					simAmountsChimeraTotCountOut<< std::endl;
					simAmountsChimeraTotCountOut << sampleKey
							<< "\t" << mixture.second->name_
							<< "\t" << totalNonChimera + totalChimera
							<< "\t" << totalChimera
							<< "\t" << totalNonChimera
							<< "\t" << getPercentageString(totalChimera, totalChimera + totalNonChimera);
					if(nullptr != mixture.second->meta_){
						for(const auto & m : metalevels){
							simAmountsChimeraTotCountOut << "\t" << mixture.second->meta_->getMeta(m);
						}
					}
					simAmountsChimeraTotCountOut << std::endl;
				}
				auto pcrSeqsInOpts = SeqIOOptions::genFastaIn(pcrReadsFnp);
				SeqInput pcrReader(pcrSeqsInOpts);
				pcrReader.openIn();
				seqInfo targetSeq;
				while(pcrReader.readNextRead(targetSeq)){
//					MetaDataInName nameMeta;
//
//					nameMeta.addMeta("Mixture", mixture.second->name_);
//					nameMeta.addMeta("Sample", sampleKey);
//					nameMeta.addMeta("readNumber", readNumber);
//					nameMeta.addMeta("refName", refName);
//					seqInfo targetSeq("", refSeqs.at(refName).seq_);
					//blunt ending
					if(lSetup.pars_.addBluntEndingArtifact_){
						if('A' == targetSeq.seq_[0] && rGen() < lSetup.pars_.bluntEndingArtifactChance_){
							readVecTrimmer::trimOffForwardBases(targetSeq, 1);
//							forwardPrimerPosition = forwardPrimerPosition == 0 ? forwardPrimerPosition: forwardPrimerPosition - 1;
//							forwardBarcodePosition = forwardBarcodePosition == 0 ? forwardBarcodePosition: forwardBarcodePosition - 1;
//							reverseBarcodePosition -= 1;
//							reversePrimerPosition -= 1;
//							nameMeta.addMeta("frontEndBluntEndArtifact", true);
						}else{
//							nameMeta.addMeta("frontEndBluntEndArtifact", false);
						}
						if ('T' == targetSeq.seq_.back()
								&& rGen() < lSetup.pars_.bluntEndingArtifactChance_) {
							readVecTrimmer::trimOffEndBases(targetSeq, 1);
//							nameMeta.addMeta("backEndBluntEndArtifact", true);
						} else {
//							nameMeta.addMeta("backEndBluntEndArtifact", false);
						}
					}
					//positions
//					if(nullptr != mixture.second->primers_){
//						nameMeta.addMeta("forwardPrimerPosition", forwardPrimerPosition);
//						nameMeta.addMeta("reversePrimerPosition", reversePrimerPosition);
//					}
//					if(nullptr != mixture.second->forwardBarcode_){
//						nameMeta.addMeta("forwardBarcodePosition", forwardBarcodePosition);
//					}
//					if(nullptr != mixture.second->reverseBarcode_){
//						nameMeta.addMeta("reverseBarcodePosition", reverseBarcodePosition);
//					}
					//complement
					if(lSetup.pars_.addReverseComplement_){
						bool complement = rGen.unifRand(0,2) == 0;
//						nameMeta.addMeta("complement", complement);
						if(complement){
							targetSeq.reverseComplementRead(false, true);
						}
					}
					//length
//					nameMeta.addMeta("targetSeqLength", len(targetSeq));
//					targetSeq.name_ = nameMeta.createMetaName();
					if(singleEnd){
						//target
						targetSeq.outPutSeq(*targetOut);
					}else{
						{
							//r1
							auto subSeq = len(targetSeq) > lSetup.pars_.pairedEndLength_? targetSeq.getSubRead(0, lSetup.pars_.pairedEndLength_) : targetSeq;
							simulator.simR1(subSeq, lSetup.pars_.pairedEndLength_).outPutFastq(*r1Out);
						}
						{
							//r2
							targetSeq.reverseComplementRead(false, true);
							auto subSeq = len(targetSeq) > lSetup.pars_.pairedEndLength_? targetSeq.getSubRead(0, lSetup.pars_.pairedEndLength_) : targetSeq;
							simulator.simR2(subSeq, lSetup.pars_.pairedEndLength_).outPutFastq(*r2Out);
						}
					}
				}
				pcrReader.closeIn();
				if(!keepPCRSeqs){
					njh::files::rmDirForce(njh::files::make_path(fastqDirectory, sampleKey + mixture.second->name_));
				}
			}
		}
	};
	std::vector<std::thread> threads;
	for(uint32_t t = 0; t < numThreads; ++t){
		threads.emplace_back(std::thread(simSample));
	}
	njh::concurrent::joinAllJoinableThreads(threads);
	//bind together the sampling information
	auto bindTables = [&verbose](const std::vector<bfs::path> & files, const bfs::path & outFnp){
		OutputStream tabOut(outFnp);
		table mainTable;
		uint32_t count = 0;
		for (const auto &file : files) {
			if (verbose) {
				std::cout << file.string() << std::endl;
			}
			if (njh::files::bfs::is_directory(file)) {
				if (verbose) {
					std::cout << "Skipping directory: " << file.string() << std::endl;
				}
				continue;
			}
			table inTab(file.string(), "\t", true);
			if (count == 0) {
				mainTable = inTab;
			} else {
				try {
					mainTable.rbind(inTab, false);
				}catch (std::exception & e) {
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", failed to add table from " << file << "\n";
					ss << e.what();
					throw std::runtime_error{ss.str()};
				}
			}
			++count;
		}
		mainTable.outPutContents(tabOut, "\t");
	};
	{
		auto allFiles = njh::files::gatherFiles(simAmountsDir, "_simAmounts.tab.txt", false);
		bfs::path outFnp(njh::files::make_path(setUp.pars_.directoryName_, "simulation_sampling_info.tab.txt"));
		bindTables(allFiles, outFnp);
	}
	{
		auto allFiles = njh::files::gatherFiles(simAmountsDir, "_chimeraSimAmounts.tab.txt", false);
		bfs::path outFnp(njh::files::make_path(chimeraDirectory, "simulation_chimeraSimAmounts_info.tab.txt"));
		bindTables(allFiles, outFnp);
	}
	{
		auto allFiles = njh::files::gatherFiles(simAmountsDir, "_chimeraSimFinals.tab.txt", false);
		bfs::path outFnp(njh::files::make_path(chimeraDirectory, "simulation_chimeraSimFinals_info.tab.txt"));
		bindTables(allFiles, outFnp);
	}
	{
		bfs::path benchmarkingFilesDir = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar{"benchmarkingFilesDir"});

		//sample to mixture
		//
		table sampleToMixtureTab(VecStr{"sample", "MixName"});

		std::unordered_map<std::string, table> sampleToMixtureRepsSeparateTables;
		std::unordered_map<std::string, table> sampleToMixtureTables;

		for(const auto & sample : lSetup.pop_.samples_){
			for(const auto & experiment : sample.second.expRuns_){
				std::string sampleName = sample.first + experiment.first;
				for(const auto & run : experiment.second){
					for(const auto & primerPair : run.hapRegionAmplified_){
						if(!njh::in(primerPair.first, sampleToMixtureTables)){
							sampleToMixtureTables[primerPair.first] = sampleToMixtureTab;
						}
						if(!njh::in(primerPair.first, sampleToMixtureRepsSeparateTables)){
							sampleToMixtureRepsSeparateTables[primerPair.first] = sampleToMixtureTab;
						}
						sampleToMixtureTables[primerPair.first].addRow(sampleName, sampleName + "--" + primerPair.first);
						sampleToMixtureRepsSeparateTables[primerPair.first].addRow(sample.first + run.runName_, sample.first + run.runName_ + "--" + primerPair.first);
					}
				}
			}
		}
		for(const auto & sampleToMixtureRepsSeparateTable : sampleToMixtureRepsSeparateTables){
			auto outTab = sampleToMixtureRepsSeparateTable.second.getUniqueRows();
			OutputStream sampToMixRepSepOut(njh::files::make_path(benchmarkingFilesDir, "samplesToMix_repsSeparate_" + sampleToMixtureRepsSeparateTable.first + ".tab.txt"));
			outTab.outPutContents(sampToMixRepSepOut, "\t");
		}

		for(const auto & sampleToMixtureTable : sampleToMixtureTables){
			auto outTab = sampleToMixtureTable.second.getUniqueRows();
			OutputStream sampToMixOut(njh::files::make_path(benchmarkingFilesDir, "samplesToMix_" + sampleToMixtureTable.first + ".tab.txt"));
			outTab.outPutContents(sampToMixOut, "\t");
		}
		//mixture setup
		//
		table mixtureSetUpTab(VecStr{"MixName", "strain", "relative_abundance"});
		std::unordered_map<std::string, table> mixtureSetUpRepsSeparateExpectedTables;
		std::unordered_map<std::string, table> mixtureSetUpExpectedTables;
		std::unordered_map<std::string, table> mixtureSetUpRepsSeparateSequencedTables;
		std::unordered_map<std::string, table> mixtureSetUpSequencedTables;
		for(const auto & sample : lSetup.pop_.samples_){

			std::map<std::string,std::map<std::string, double>> hapExpected;
			for(const auto & region : sample.second.strainsExpected_){
				for(const auto & hap : region.second.genomicRegionToHapNames_){
					hapExpected[hap.first][hap.second] += region.second.relativeAbundance_;
				}
			}
			for(const auto & experiment : sample.second.expRuns_){
				std::string sampleName = sample.first + experiment.first;
				for(const auto & region : hapExpected){
					if(!njh::in(region.first, mixtureSetUpExpectedTables)){
						mixtureSetUpExpectedTables[region.first] = mixtureSetUpTab;
					}
					for(const auto & hap : region.second){
						mixtureSetUpExpectedTables[region.first].addRow(sampleName + "--" + region.first, hap.first, hap.second);
					}
				}

				auto combined = ControlPopulation::Sample::ExperimentRun::averageHapAbundSequenced(experiment.second);
				for(const auto & region : combined){
					double total = 0;
					for(const auto & hap : region.second){
						total += hap.second;
					}
					for(const auto & hap : region.second){
						if(!njh::in(region.first, mixtureSetUpSequencedTables)){
							mixtureSetUpSequencedTables[region.first] = mixtureSetUpTab;
						}
						mixtureSetUpSequencedTables[region.first].addRow(sampleName + "--" + region.first, hap.first, hap.second/total);
					}
				}

				for(const auto & run : experiment.second){
					for(const auto & region : hapExpected){
						if(!njh::in(region.first, mixtureSetUpRepsSeparateExpectedTables)){
							mixtureSetUpRepsSeparateExpectedTables[region.first] = mixtureSetUpTab;
						}
						for(const auto & hap : region.second){
							mixtureSetUpRepsSeparateExpectedTables[region.first].addRow(sample.first + run.runName_ + "--" + region.first, hap.first, hap.second);
						}
					}
					for(const auto & primerPair : run.hapRegionAmplified_){
						if(!njh::in(primerPair.first, mixtureSetUpRepsSeparateSequencedTables)){
							mixtureSetUpRepsSeparateSequencedTables[primerPair.first] = mixtureSetUpTab;
						}
						std::string mixname = sample.first + run.runName_ + "--" + primerPair.first;
						double total = 0;
						for(const auto & hap : primerPair.second.hapAbundSequenced_){
							total += hap.second;
						}
						for(const auto & hap : primerPair.second.hapAbundSequenced_){
							mixtureSetUpRepsSeparateSequencedTables[primerPair.first].addRow(mixname, hap.first, hap.second/total);
						}
					}
				}
			}
		}

		for(const auto & mixtureSetUpRepsSeparateExpectedTable : mixtureSetUpRepsSeparateExpectedTables){
			auto outTab = mixtureSetUpRepsSeparateExpectedTable.second.getUniqueRows();
			OutputStream sampToMixRepSepOut(njh::files::make_path(benchmarkingFilesDir, "mixtureSetUpRepsSeparateExpected_" + mixtureSetUpRepsSeparateExpectedTable.first + ".tab.txt"));
			outTab.outPutContents(sampToMixRepSepOut, "\t");
		}

		for(const auto & mixtureSetUpExpectedTable : mixtureSetUpExpectedTables){
			auto outTab = mixtureSetUpExpectedTable.second.getUniqueRows();
			OutputStream sampToMixOut(njh::files::make_path(benchmarkingFilesDir, "mixtureSetUpExpected_" + mixtureSetUpExpectedTable.first + ".tab.txt"));
			outTab.outPutContents(sampToMixOut, "\t");
		}

		for(const auto & mixtureSetUpRepsSeparateSequencedTable : mixtureSetUpRepsSeparateSequencedTables){
			auto outTab = mixtureSetUpRepsSeparateSequencedTable.second.getUniqueRows();
			OutputStream sampToMixRepSepOut(njh::files::make_path(benchmarkingFilesDir, "mixtureSetUpRepsSeparateSequenced_" + mixtureSetUpRepsSeparateSequencedTable.first + ".tab.txt"));
			outTab.outPutContents(sampToMixRepSepOut, "\t");
		}

		for(const auto & mixtureSetUpSequencedTable : mixtureSetUpSequencedTables){
			auto outTab = mixtureSetUpSequencedTable.second.getUniqueRows();
			OutputStream sampToMixOut(njh::files::make_path(benchmarkingFilesDir, "mixtureSetUpSequenced_" + mixtureSetUpSequencedTable.first + ".tab.txt"));
			outTab.outPutContents(sampToMixOut, "\t");
		}
	}
	//sampleToMixture,


	njh::files::rmDirForce(simAmountsDir);
	return 0;
}

} // namespace njhseq
