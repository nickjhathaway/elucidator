/*
 * simMultipleMixtureSimPCR.cpp
 *
 *  Created on: Oct 7, 2019
 *      Author: nicholashathaway
 */




#include "readSimulatorRunner.hpp"
#include "elucidator/simulation.h"
#include <SeekDeep/objects/TarAmpSetupUtils/PrimersAndMids.hpp>
#include <njhseq/IO/SeqIO/SeqIO.hpp>
#include <njhseq/readVectorManipulation/readVectorHelpers/readVecTrimmer.hpp>

namespace njhseq {




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
	uint32_t defaultPcrRounds = 30;
	uint32_t initialPcrRounds = 10;
	std::map<uint32_t, uint32_t> initialPcrRoundsMap;
	bfs::path initialPcrRoundsMapFnp;
	long double errorRate = 3.5e-06;
	double pcrEfficiency = 0.85;
	bool keepPCRSeqs = false;
	uint32_t chimeraBasesIn = 5;
	uint32_t templateCap = 500000000;
	bool noChimeras = false;
	double finalReadAmountSDFrac = 0.1;

	uint32_t pairedEndLength = std::numeric_limits<uint32_t>::max();
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(pairedEndLength, "--pairedEndLength", "Paired End Length");
	setUp.setOption(finalReadAmountSDFrac, "--finalReadAmountSDFrac", "final Read Amount SD Frac", njh::progutils::ProgramSetUp::CheckCase::GREATERZERO);


	setUp.setOption(noChimeras, "--noChimeras", "Don't simulate chimeras");
	setUp.setOption(templateCap, "--templateCap", "Template Cap");
	setUp.setOption(chimeraBasesIn, "--chimeraBasesIn", "The number of bases needed for a template to lay down");
	setUp.setOption(keepPCRSeqs, "--keepPCRSeqs", "Keep PCR Seqs");
	setUp.setOption(errorRate, "--errorRate", "Polymerase Error Rate");
	setUp.setOption(pcrEfficiency, "--pcrEfficiency", "PCR Efficiency, between 0-1, chance a product gets amplified");
	setUp.setOption(defaultPcrRounds, "--pcrRounds", "Number of PCR rounds");
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
		if("NoPrimers" == refSeqs.first || !lSetup.pars_.noAddPrimers_){
			for(const auto & refSeq : refSeqs.second){
				seqWriter.openWrite(refSeq);
			}
		}else{
			for(const auto & refSeq : refSeqs.second){
				seqInfo copySeqToTrim = *refSeq;
				readVecTrimmer::trimBetweenSequences(copySeqToTrim,
						lSetup.ids_->pDeterminator_->primers_.at(refSeqs.first).fwds_.front().info_ ,
						lSetup.ids_->pDeterminator_->primers_.at(refSeqs.first).revs_.front().infoRC_,
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

	std::function<void()> simSample = [&samplesQueue,&extension,&sampleCountAtom,&fastqDirectory,
										&lSetup,&illuminaProfileDir,&singleEnd,&refSeqs,
										&intErrorRate, &defaultPcrRounds, &initialPcrRounds,&verbose,
										&keepPCRSeqs, &pcrNumThreads,&pcrEfficiency,&simAmountsDir,
										&chimeraDirectory,&templateCap,&noChimeras,&chimeraBasesIn,&initialPcrRoundsMap,
										&finalReadAmountSDFrac](){
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
				uint64_t readSimAmount = mixture.second->finalReadAmount_;
				if(std::numeric_limits<double>::max() != finalReadAmountSDFrac){
					std::normal_distribution<double> ndist(readSimAmount, finalReadAmountSDFrac * readSimAmount);
					double newReadSimAmount = 0;
					while(newReadSimAmount <= 0){
						newReadSimAmount = std::round(ndist(rGen.mtGen_));
					}
					readSimAmount = newReadSimAmount;
				}
				uint32_t pcrRoundsForMixture = defaultPcrRounds;
				if(std::numeric_limits<uint32_t>::max() != mixture.second->pcrRounds_) {
					pcrRoundsForMixture = mixture.second->pcrRounds_;
				}
				auto pcrSimAmounts = pcrSim.simLibFast(seqGCounts, pcrReadsOpts, readSimAmount, pcrRoundsForMixture, currentInitial, pcrNumThreads);
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
	njh::concurrent::runVoidFunctionThreaded(simSample, numThreads);
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
