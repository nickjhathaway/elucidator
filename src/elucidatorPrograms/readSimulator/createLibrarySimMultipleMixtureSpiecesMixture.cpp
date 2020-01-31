/*
 * createLibrarySimMultipleMixtureSpiecesMixture.cpp
 *
 *  Created on: Jan 30, 2020
 *      Author: nicholashathaway
 */




#include "readSimulatorRunner.hpp"
#include "elucidator/simulation.h"
#include <SeekDeep/objects/TarAmpSetupUtils/PrimersAndMids.hpp>

namespace njhseq {



int readSimulatorRunner::createLibrarySimMultipleMixtureSpiecesMixture(
		const njh::progutils::CmdArgs & inputCommands) {

	/**@todo this approach needs a re-vamp to handle what is being attempt which is having multiple regions (like 18S rRNA, etc) amplfied using one primer pair */
//	LibrarySetup::SimLibrarySetupPars simPars;
//	bfs::path primerMidFnp = "";
//
//	bfs::path mixtureSetupFile = "";
//	bfs::path strainSeqsFnp = "";
//
//	std::string libraryName = "";
//	bool simReplicates = false;
//	PCRAmountPars pcrNumbers;
//	double adjustAroundGivenFrac = 0.5;
//	double twoSdFrac = std::numeric_limits<double>::max();
//
//	readSimulatorSetUp setUp(inputCommands);
//	setUp.processDebug();
//	setUp.processVerbose();
//	setUp.setOption(pcrNumbers.startingTemplateAmounts_, "--startingTemplateAmount", "Starting PCR Template Amount", njh::progutils::ProgramSetUp::ConCheckCase::NONZERO);
//	setUp.setOption(pcrNumbers.finalReadAmount_, "--perMixtureReadAmount", "Final Read Amount to create per mixture", njh::progutils::ProgramSetUp::ConCheckCase::NONZERO);
//	setUp.setOption(simReplicates, "--replicates", "Replicates");
//	setUp.setOption(adjustAroundGivenFrac, "--adjustAroundGivenFrac", "Adjust Around Given Frac", njh::progutils::ProgramSetUp::CheckCase::GREATERZERO);
//	adjustAroundGivenFrac = std::min(adjustAroundGivenFrac, 0.99);
//	setUp.setOption(twoSdFrac, "--twoSdFrac", "Two Sd Fracs around genome starting template amount", njh::progutils::ProgramSetUp::CheckCase::GREATERZERO);
//
//	setUp.setOption(mixtureSetupFile, "--mixtureSetupFile", "The mixtrues of strains, 3 columns. MixName, strain, RelativeAbundance", true);
//	setUp.setOption(strainSeqsFnp, "--strainSeqs", "Table of sequences for given strain, 4 columns needed, strain, target, seq, name");
//
//
//	setUp.setOption(primerMidFnp, "--primerMidFnp", "Primer MID Fnp", true);
//
//	setUp.setOption(simPars.pairedEndLength_, "--pairedEndLength", "Paired End Length");
//	setUp.setOption(simPars.noAddPrimers_, "--noAddPrimers", "Primers are already present");
//	setUp.setOption(simPars.barcodeRandomPrecedingBases_, "--barcodeRandomPrecedingBases", "Barcode Random Preceding Bases");
//	setUp.setOption(simPars.primerRandomPrecedingBases_, "--primerRandomPrecedingBases", "Primer Random Preceding Bases");
//	setUp.setOption(simPars.addReverseComplement_, "--addReverseComplement", "Add Reverse Complement");
//	setUp.setOption(simPars.addBluntEndingArtifact_, "--addBluntEndingArtifact", "Add Blunt Ending Artifact");
//	setUp.setOption(simPars.bluntEndingArtifactChance_, "--bluntEndingArtifactChance", "Blunt Ending Artifact Chance");
//	setUp.setOption(libraryName, "--libraryName", "Library Name", true);
//	setUp.processDirectoryOutputName(libraryName, true);
//	setUp.finishSetUp(std::cout);
//	setUp.startARunLog(setUp.pars_.directoryName_ );
//
//	LibrarySetup lSetup(libraryName, simPars);
//
//	table coiTable;;
//	table mixtureSetupTab(mixtureSetupFile, "\t", true);
//	table strainSeqsTab(strainSeqsFnp, "\t", true);
//
//	//key1 = target, key2 = seq, value = input names
//	std::unordered_map<std::string, std::unordered_map<std::string, VecStr>> hapsByTarBySeq_;
//	//key1 = strain, key2 = target, value = vector of seqs for target for strain
//	std::unordered_map<std::string, std::unordered_map<std::string, VecStr>> byStrainByTarSeqs_;
//	for(const auto & row : strainSeqsTab){
//		auto target = row[strainSeqsTab.getColPos("target")];
//		auto seq = row[strainSeqsTab.getColPos("seq")];
//		auto strain = row[strainSeqsTab.getColPos("strain")];
//		auto name = row[strainSeqsTab.getColPos("name")];
//		byStrainByTarSeqs_[strain][target].emplace_back(seq);
//		hapsByTarBySeq_[target][seq].emplace_back(njh::pasteAsStr(strain, "__", name));
//	}
//
//	//key1 = target, key2 = seq, value = output name
//	std::unordered_map<std::string, std::unordered_map<std::string, std::string>> byTarSeqToName_;
//
//	{
//		//writing out seqs;
//
//		SeqOutput allSeqsWriter(SeqIOOptions::genFastaOut(njh::files::make_path(setUp.pars_.directoryName_, "allSeqs.fasta")));
//		allSeqsWriter.openOut();
//		for(const auto & tar : hapsByTarBySeq_){
//			VecStr seqs = getVectorOfMapValues(tar.second);
//			njh::sort(seqs,[&tar](const std::string & seq1, const std::string & seq2){
//				if(tar.second.at(seq1).size() == tar.second.at(seq2).size()){
//					return seq1 < seq2;
//				}else{
//					return tar.second.at(seq1).size() > tar.second.at(seq2).size();
//				}
//			});
//			uint32_t seqCount = 0;
//			for(const auto & seq : seqs){
//				std::string newName = njh::pasteAsStr(tar.first, ".", njh::leftPadNumStr<uint32_t>(seqCount, seqs.size()));
//				++seqCount;
//				allSeqsWriter.write(seqInfo(newName, seq));
//				byTarSeqToName_[tar.first][seq] = newName;
//			}
//		}
//	}
//
//
//
//
//	lSetup.ids_ = std::make_unique<PrimersAndMids>(primerMidFnp);
//	lSetup.ids_->initPrimerDeterminator();
//	if(lSetup.ids_->containsMids()){
//		lSetup.ids_->initMidDeterminator(MidDeterminator::MidDeterminePars{});
//	}
//
//	njh::for_each(mixtureSetupTab.columnNames_, [](std::string & col){
//		njh::strToLower(col);
//	});
//	mixtureSetupTab.setColNamePositions();
//	njh::for_each(strainSeqsTab.columnNames_, [](std::string & col){
//		njh::strToLower(col);
//	});
//	strainSeqsTab.setColNamePositions();
//	mixtureSetupTab.checkForColumnsThrow(VecStr{"MixName", "strain", "RelativeAbundance"}, __PRETTY_FUNCTION__);
//	strainSeqsTab.checkForColumnsThrow(VecStr{"strain", "target", "seq", "name"}, __PRETTY_FUNCTION__);
//
//	auto targetsFromSeqTab = strainSeqsTab.getColumnLevels("target");
//
//	//primer check
//	VecStr missingPrimerPairs;
//	for(const auto & primerPair : targetsFromSeqTab){
//		if(!njh::in(primerPair, lSetup.ids_->pDeterminator_->primers_)){
//			missingPrimerPairs.emplace_back(primerPair);
//		}
//	}
//	if (!missingPrimerPairs.empty()) {
//		std::stringstream ss;
//		ss << __PRETTY_FUNCTION__ << " error, missing primer pair names in "
//				<< primerMidFnp << "\n";
//		ss << "Missing: " << njh::conToStr(missingPrimerPairs) << "\n";
//		throw std::runtime_error { ss.str() };
//	}
//
//
//
//	struct StrainMixture {
//		StrainMixture(const std::string & name,
//				const std::map<std::string, double> & hapFracs,
//				const MetaDataInName & meta) :
//				name_(name), hapFracs_(hapFracs), meta_(meta) {
//		}
//		std::string name_;
//		std::unordered_map<std::string, double> hapFracs_;
//		MetaDataInName meta_;
//		bool replicate_{false};
//		Json::Value toJson() const{
//			Json::Value ret;
//			ret["class"] = njh::json::toJson(njh::getTypeName(*this));
//			ret["name_"] = njh::json::toJson(name_);
//			ret["meta_"] = njh::json::toJson(meta_.meta_);
//			ret["hapFracs_"] = njh::json::toJson(hapFracs_);
//			ret["replicate_"] = njh::json::toJson(replicate_);
//			return ret;
//		}
//	};
//	std::vector<StrainMixture> mixtureSetUpPars;
//	std::unordered_map<std::string, std::map<std::string, double>> hapMixInfos;
//	for(const auto & row : mixtureSetupTab){
//		auto name = row[mixtureSetupTab.getColPos("mixname")];
//		auto strain = row[mixtureSetupTab.getColPos("strain")];
//		auto amount = row[mixtureSetupTab.getColPos("relativeabundance")];
//		hapMixInfos[name][strain] = njh::StrToNumConverter::stoToNum<double>(amount);
//	}
//	uint32_t totalMixture = 0;
//
//	for(const auto & hapMix : hapMixInfos){
//		mixtureSetUpPars.emplace_back(StrainMixture(hapMix.first, hapMix.second, MetaDataInName()));
//		++totalMixture;
//	}
//
//	OutOptions sampleNamesOutOpts(
//			njh::files::make_path(setUp.pars_.directoryName_, "sampleNames.tab.txt"));
//	VecStr sampleNamesTabHeader{};
//	if (!lSetup.ids_->containsMids()) {
//		sampleNamesTabHeader = {"#target", "sample", "run1"};
//	}else{
//		sampleNamesTabHeader = {"#file", "sample", "run1"};
//	}
//	bool anyReps = false;
//	for (const auto & patient : mixtureSetUpPars) {
//		if(patient.replicate_){
//			anyReps = true;
//			break;
//		}
//	}
//	if(anyReps){
//		sampleNamesTabHeader.emplace_back("run2");
//	}
//
//
//	auto maxStartingTemplateAmounts = vectorMaximum(pcrNumbers.startingTemplateAmounts_);
//	auto maxFinalReadAmount = vectorMaximum(pcrNumbers.finalReadAmount_);
//
//
//	auto genStartingTempNameSection = [&pcrNumbers,&maxStartingTemplateAmounts](uint32_t startingTemplateAmount){
//		std::string name;
//		if(pcrNumbers.startingTemplateAmounts_.size() > 1){
//			name = njh::pasteAsStr("-ST", njh::leftPadNumStr(startingTemplateAmount, maxStartingTemplateAmounts));
//		}
//		return name;
//	};
//	auto genFinalReadNameSection = [&pcrNumbers,&maxFinalReadAmount](uint32_t finalReadAmount){
//		std::string name;
//		if(pcrNumbers.finalReadAmount_.size() > 1){
//			name = njh::pasteAsStr("-RD", njh::leftPadNumStr(finalReadAmount, maxFinalReadAmount));
//		}
//		return name;
//	};
//	for (const auto & mixture : mixtureSetUpPars) {
//
//		std::map<std::string, ControlPopulation::Strain> strains;
//		double totalFracs = 0;
//		for (const auto & hap : mixture.hapFracs_) {
//			totalFracs += hap.second;
//		}
//		for (const auto & hap : mixture.hapFracs_) {
//			strains.emplace(hap.first, ControlPopulation::Strain(hap.first, haplotypeKey[hap.first], hap.second/totalFracs));
//		}
//		ControlPopulation::Sample currentSamp(mixture.name_, strains);
//		currentSamp.meta_.addMeta(mixture.meta_, false);
//		currentSamp.meta_.addMeta("MixName", mixture.name_);
//
//		for (const auto & startingTemplateAmount : pcrNumbers.startingTemplateAmounts_){
//			for (const auto & finalReadAmount : pcrNumbers.finalReadAmount_){
//				uint32_t numRuns = mixture.replicate_ ? 2 : 1;
//				uint32_t startingTemplateAmountForSamp = startingTemplateAmount;
//				if(std::numeric_limits<double>::max() != twoSdFrac){
//					double startingTemplateAmountForSampNew = 0;
//					std::normal_distribution<double> ndist(startingTemplateAmountForSamp, (startingTemplateAmountForSamp * twoSdFrac)/2);
//			    // random device class instance, source of 'true' randomness for initializing random seed
//			    std::random_device rd;
//			    // Mersenne twister PRNG, initialized with seed from previous random device instance
//			    std::mt19937 rgen(rd());
//			    while(startingTemplateAmountForSampNew <= 1){
//			    	startingTemplateAmountForSampNew = std::round(ndist(rgen));
//			    }
//			    startingTemplateAmountForSamp = startingTemplateAmountForSampNew;
//				}
//				currentSamp.addExpRun(
//						genStartingTempNameSection(startingTemplateAmount) + genFinalReadNameSection(finalReadAmount),
//						ControlPopulation::SeqSampledAmounts { finalReadAmount,startingTemplateAmountForSamp }, numRuns);
//			}
//		}
//		lSetup.pop_.addSample(currentSamp);
//	}
//
//	table sampleNamesTables(sampleNamesTabHeader);
//	if (!lSetup.ids_->containsMids()) {
//		auto sampleKeys = njh::getVecOfMapKeys(lSetup.pop_.samples_);
//		njh::sort(sampleKeys);
//		for(const auto & sampleKey : sampleKeys){
//			const auto & sample = lSetup.pop_.samples_.at(sampleKey);
//			for(const auto & experiment : sample.expRuns_){
//				std::string outSampName = sample.sampName_ + experiment.first;
//				//add to sample name table
//				for(const auto & primerPair : targetsFromSeqTab){
//					VecStr addingRow{primerPair, outSampName};
//					for(const auto & run : experiment.second){
//						addingRow.emplace_back(sample.sampName_ + run.runName_);
//					}
//					sampleNamesTables.addRow(addingRow);
//				}
//				//add experiments to sequence setup;
//				uint32_t runNumber = 0;
//				for(const auto & run : experiment.second){
//					auto sampleSet = std::make_shared<SampleSetup>(sample.sampName_ + run.runName_);
//					std::unordered_map<std::string, std::shared_ptr<MixtureSetUp> > mixtures;
//					for (const auto & primerPair : targetsFromSeqTab) {
//						mixtures[primerPair] = std::make_shared<MixtureSetUp>(primerPair);
//						mixtures[primerPair]->meta_ = std::make_unique<MetaDataInName>();
//						mixtures[primerPair]->meta_->addMeta("PrimerPair", primerPair);
//						mixtures[primerPair]->meta_->addMeta("PatientSample", sample.sampName_ + experiment.first);
//						//add in identifying information
//						mixtures[primerPair]->meta_->addMeta("sampName_", sample.sampName_);
//						mixtures[primerPair]->meta_->addMeta("experiment", experiment.first );
//						mixtures[primerPair]->meta_->addMeta("runNumber", runNumber);
//
//						mixtures[primerPair]->startingTemplateAmount_ = run.expAmounts_.totalGenomesSampled_;
//						mixtures[primerPair]->finalReadAmount_ = run.expAmounts_.sequencedReadAmount_;
//						mixtures[primerPair]->setPrimers(primerPair,
//								njh::mapAt(lSetup.ids_->targets_, primerPair).info_.forwardPrimerRaw_,
//								njh::mapAt(lSetup.ids_->targets_, primerPair).info_.reversePrimerRaw_);
//						mixtures[primerPair]->primers_->reverse_randomPrecedingBases_ = simPars.primerRandomPrecedingBases_;
//						mixtures[primerPair]->primers_->forward_randomPrecedingBases_ = simPars.primerRandomPrecedingBases_;
//					}
////					std::cout << sampleKey << std::endl;
////					std::cout << "sample.strainsExpected_: " << sample.strainsExpected_.size() << std::endl;
////					std::cout << sample.sampName_ << run.runName_ << std::endl;
////					std::cout << "\tsample.strainsExpected_.size():" << sample.strainsExpected_.size() << std::endl;
//					for (const auto & strain : sample.strainsExpected_) {
////						std::cout << strain.first << std::endl;
////						std::cout << "\tstrain.second.name_: " << strain.second.name_ << std::endl;
////						std::cout << "\tstrain.second.relativeAbundance_: " << strain.second.relativeAbundance_ << std::endl;
////						std::cout << "\tstrain.second.genomicRegionToHapNames_.size(): " << strain.second.genomicRegionToHapNames_.size() << std::endl;
//						uint32_t genomeCountsForStrain = run.hapAbundGenomesSampled_.at(strain.second.name_).relativeAbundance_;
////						std::cout << strain.first << "\t" << genomeCountsForStrain << std::endl;
//						for (const auto & subHapName : strain.second.genomicRegionToHapNames_) {
////							std::cout << strain.first << std::endl;
////							std::cout << '\t' << subHapName.first << std::endl;
////							std::cout << "\t\t" << subHapName.second << std::endl;
////							std::cout << "\t\t" << "strain.second.relativeAbundance_: " << strain.second.relativeAbundance_ << std::endl;
//							if (njh::in(subHapName.second, mixtures[subHapName.first]->expectedAbundances_)) {
//								mixtures[subHapName.first]->expectedAbundances_[subHapName.second] += strain.second.relativeAbundance_;
//								mixtures[subHapName.first]->genomeCounts_[subHapName.second] += genomeCountsForStrain;
//							} else {
//								mixtures[subHapName.first]->addAbundance(subHapName.second, strain.second.relativeAbundance_);
//								mixtures[subHapName.first]->addGenomeCount(subHapName.second, genomeCountsForStrain);
//							}
//						}
//					}
//					for (const auto & mix : mixtures) {
////						std::cout << sample.sampName_ << run.runName_ << "-" << mix.first << std::endl;
////						for(const auto & hap : mix.second->genomeCounts_){
////							std::cout << "\t" << hap.first << "\t" << hap.second << std::endl;
////						}
//						sampleSet->addMixture(mix.second);
//					}
//					lSetup.addSample(sampleSet);
//					++runNumber;
//				}
//			}
//		}
//
//	}else{
//		uint32_t finalSubSetAmount = 0;
//		{
//			uint32_t midCount = 0;
//			for (const auto & patient : lSetup.pop_.samples_) {
//				for(const auto & experiment : patient.second.expRuns_){
//					for(uint32_t run = 0; run < experiment.second.size(); ++run){
//						++midCount;
//						if(midCount >= lSetup.ids_->mids_.size()){
//							midCount = 0;
//							++finalSubSetAmount;
//						}
//					}
//				}
//			}
//		}
//		uint32_t indexCount = 0;
//		uint32_t midCount = 0;
//		auto sampleSet = std::make_shared<SampleSetup>(njh::pasteAsStr("Subset-", njh::leftPadNumStr(indexCount, finalSubSetAmount) ) );
//		auto midNames = getVectorOfMapKeys(lSetup.ids_->mids_);
//		njh::sort(midNames);
//		auto sampleKeys = njh::getVecOfMapKeys(lSetup.pop_.samples_);
//		njh::sort(sampleKeys);
//
//		for(const auto & sampleKey : sampleKeys){
//			const auto & sample = lSetup.pop_.samples_.at(sampleKey);
//			for(const auto & experiment : sample.expRuns_){
//				std::string outSampName = sample.sampName_ + experiment.first;
//				std::unordered_map<std::string, VecStr> replicateNames;
//				uint32_t runNumber = 0;
//				for(const auto & run : experiment.second){
//					std::string midName =  midNames[midCount];
//					std::string indexName = sampleSet->name_;
//					replicateNames[indexName].emplace_back(midName);
//					std::unordered_map<std::string, std::shared_ptr<MixtureSetUp> > mixtures;
//					for (const auto & primerPair : targetsFromSeqTab) {
//						mixtures[primerPair] = std::make_shared<MixtureSetUp>(njh::pasteAsStr(primerPair, "-", midName));
//						mixtures[primerPair]->meta_ = std::make_unique<MetaDataInName>();
//						mixtures[primerPair]->meta_->addMeta("PrimerPair", primerPair);
//						mixtures[primerPair]->meta_->addMeta("PatientSample", outSampName);
//						//add in identifying information
//						mixtures[primerPair]->meta_->addMeta("sampName_", sample.sampName_);
//						mixtures[primerPair]->meta_->addMeta("experiment", experiment.first );
//						mixtures[primerPair]->meta_->addMeta("runNumber", runNumber);
//
//						mixtures[primerPair]->startingTemplateAmount_ = run.expAmounts_.totalGenomesSampled_;
//						mixtures[primerPair]->finalReadAmount_ = run.expAmounts_.sequencedReadAmount_;
//						mixtures[primerPair]->setPrimers(primerPair,
//								njh::mapAt(lSetup.ids_->targets_, primerPair).info_.forwardPrimerRaw_,
//								njh::mapAt(lSetup.ids_->targets_, primerPair).info_.reversePrimerRaw_);
//						mixtures[primerPair]->primers_->reverse_randomPrecedingBases_ = simPars.primerRandomPrecedingBases_;
//						mixtures[primerPair]->primers_->forward_randomPrecedingBases_ = simPars.primerRandomPrecedingBases_;
//					}
//					for (const auto & strain : sample.strainsExpected_) {
//						for (const auto & subHapName : strain.second.genomicRegionToHapNames_) {
//							uint32_t genomeCountsForSubRegion = run.hapRegionAmplified_.at(subHapName.first).hapAbundSampled_.at(subHapName.second);
//							if (njh::in(subHapName.second, mixtures[subHapName.first]->expectedAbundances_)) {
//								mixtures[subHapName.first]->expectedAbundances_[subHapName.second] += strain.second.relativeAbundance_;
//								mixtures[subHapName.first]->genomeCounts_[subHapName.second] += genomeCountsForSubRegion;
//							} else {
//								mixtures[subHapName.first]->addAbundance(subHapName.second, strain.second.relativeAbundance_);
//								mixtures[subHapName.first]->addGenomeCount(subHapName.second, genomeCountsForSubRegion);
//							}
//						}
//					}
//					for (auto & mix : mixtures) {
//						if(nullptr != lSetup.ids_->mids_.at(midName).forwardBar_){
//							mix.second->setForwardBarcode(midName, lSetup.ids_->mids_.at(midName).forwardBar_->bar_->motifOriginal_);
//							mix.second->forwardBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
//						}
//						if(nullptr != lSetup.ids_->mids_.at(midName).reverseBar_){
//							mix.second->setReverseBarcode(midName, lSetup.ids_->mids_.at(midName).reverseBar_->bar_->motifOriginal_);
//							mix.second->reverseBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
//						}
//						sampleSet->addMixture(mix.second);
//					}
//					++midCount;
//					if(midCount >= lSetup.ids_->mids_.size()){
//						lSetup.addSample(sampleSet);
//						midCount = 0;
//						++indexCount;
//						sampleSet = std::make_shared<SampleSetup>(njh::pasteAsStr("Subset-", njh::leftPadNumStr(indexCount, finalSubSetAmount) ) );
//					}
//					++runNumber;
//				}
//				for(const auto & index : replicateNames){
//					VecStr addingRow{index.first, outSampName};
//					for(const auto & run : index.second){
//						addingRow.emplace_back(run);
//					}
//					sampleNamesTables.addRow(addingRow);
//				}
//			}
//		}
//		if(midCount > 0){
//			lSetup.addSample(sampleSet);
//		}
//	}
//	njh::sort(sampleNamesTables.content_, [](const VecStr & row1, const VecStr & row2){
//		if(row1[0] == row2[0]){
//			if(row1[1] == row2[1]){
//				return row1[2] < row2[2];
//			}else{
//				return row1[1] < row2[1];
//			}
//		}else{
//			return row1[0] < row2[0];
//		}
//	});
//
//	OutputStream sampleNamesOut(sampleNamesOutOpts);
//	sampleNamesTables.outPutContents(sampleNamesOut, "\t");
//
//	VecStr metaHeader{"sample"};
//	if(pcrNumbers.finalReadAmount_.size() > 1 || pcrNumbers.startingTemplateAmounts_.size() > 1){
//		metaHeader.emplace_back("MixName");
//	}
//	metaHeader.emplace_back("PCRStartingTemplate");
//	metaHeader.emplace_back("FinalSamplingReadAmount");
//
//	//addOtherVec(metaHeader, metaLevels);
//	table metaTable(metaHeader);
//	for(const auto & mixture : mixtureSetUpPars){
//
//		for (const auto & startingTemplateAmount : pcrNumbers.startingTemplateAmounts_){
//			for (const auto & finalReadAmount : pcrNumbers.finalReadAmount_){
//				MetaDataInName sampMeta = mixture.meta_;
//				std::string outputSampName = mixture.name_
//						+ genStartingTempNameSection(startingTemplateAmount)
//						+ genFinalReadNameSection(finalReadAmount);
//				VecStr row;
//				row.emplace_back(outputSampName);
//				if(pcrNumbers.finalReadAmount_.size() > 1 || pcrNumbers.startingTemplateAmounts_.size() > 1){
//					sampMeta.addMeta("MixName", mixture.name_);
//					row.emplace_back(estd::to_string(mixture.name_));
//				}
//				sampMeta.addMeta("PCRStartingTemplate", startingTemplateAmount);
//				sampMeta.addMeta("FinalSamplingReadAmount", finalReadAmount);
//				row.emplace_back(estd::to_string(startingTemplateAmount));
//				row.emplace_back(estd::to_string(finalReadAmount));
////				for(const auto & head : metaLevels){
////					row.emplace_back(sampMeta.getMeta(head));
////				}
//				metaTable.addRow(row);
//			}
//		}
//
//	}
//	OutOptions libSetUpOpts(njh::files::make_path(setUp.pars_.directoryName_, "librarySetup.json"));
//	OutputStream libSetUpOut(libSetUpOpts);
//	libSetUpOut << njh::json::toJson(lSetup) << std::endl;
//
//	OutOptions metaOutOpts(njh::files::make_path(setUp.pars_.directoryName_, "metaData.tab.txt"));
//	OutputStream metaOut(metaOutOpts);
//	metaTable.sortTable("sample", false);
//	metaTable.outPutContents(metaOut, "\t");
//
//	auto abundTab = lSetup.createAbundanceTable();
//	OutOptions abundTabOutOpts(njh::files::make_path(setUp.pars_.directoryName_, "abundanceTable.tab.txt"));
//	OutputStream abundTabOut(abundTabOutOpts);
//	njh::sort(abundTab.content_, [](const VecStr & row1, const VecStr & row2){
//		if(row1[1] == row2[1]){
//			if(row1[2] == row2[2]){
//				return row1[3] < row2[3];
//			}else{
//				return row1[2] < row2[2];
//			}
//		}else{
//			return row1[1] < row2[1];
//		}
//	});
//	abundTab.outPutContents(abundTabOut, "\t");
	return 0;
}






} // namespace njhseq


