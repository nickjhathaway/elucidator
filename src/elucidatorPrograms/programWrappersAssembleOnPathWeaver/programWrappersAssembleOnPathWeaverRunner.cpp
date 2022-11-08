/*
 * programWrappersAssembleOnPathWeaverRunner.cpp
 *
 *  Created on: Sep 8, 2021
 *      Author: nick
 */


#include "elucidator/objects/seqContainers.h"
#include "elucidator/objects/BioDataObject.h"
#include "elucidator/seqToolsUtils/seqToolsUtils.hpp"


#include <njhseq/objects/seqObjects/seqKmers.h>
#include <njhseq/system.h>
#include <PathWeaver/objects/bam/RegionInvestigatorInBam.hpp>
#include <utility>

#include "programWrappersAssembleOnPathWeaverRunner.hpp"
#include "otherAssemblersUtils.hpp"

namespace njhseq {

programWrappersAssembleOnPathWeaverRunner::programWrappersAssembleOnPathWeaverRunner()
    : njh::progutils::ProgramRunner(
          {

					 addFunc("runSpadesOnPathWeaverRegions", runSpadesOnPathWeaverRegions, false),
					 addFunc("runUnicyclerOnPathWeaverRegions", runUnicyclerOnPathWeaverRegions, false),
					 addFunc("runMegahitOnPathWeaverRegions", runMegahitOnPathWeaverRegions, false),
					 addFunc("runPRICEOnPathWeaverRegions", runPRICEOnPathWeaverRegions, false),
					 addFunc("runSavageOnPathWeaverRegions", runSavageOnPathWeaverRegions, false),
					 addFunc("runPolyteOnPathWeaverRegions", runPolyteOnPathWeaverRegions, false),

					 addFunc("runVelvetOptimizerAndMetaVelvetOnPathWeaverRegions", runVelvetOptimizerAndMetaVelvetOnPathWeaverRegions, false),
					 addFunc("runTrinityOnPathWeaverRegions", runTrinityOnPathWeaverRegions, false),
					 addFunc("runIDBAUDOnPathWeaverRegions", runIDBAUDOnPathWeaverRegions, false),
					 addFunc("runRayOnPathWeaverRegions", runRayOnPathWeaverRegions, false),
					 addFunc("runMIRAOnPathWeaverRegions", runMIRAOnPathWeaverRegions, false),
					 addFunc("runFermiLiteOnPathWeaverRegions", runFermiLiteOnPathWeaverRegions, false),

					 addFunc("runSpadesOnPathWeaverRegionsAndUnmapped", runSpadesOnPathWeaverRegionsAndUnmapped, false),
					 addFunc("runUnicyclerOnPathWeaverRegionsAndUnmapped", runUnicyclerOnPathWeaverRegionsAndUnmapped, false),
					 addFunc("runMegahitOnPathWeaverRegionsAndUnmapped", runMegahitOnPathWeaverRegionsAndUnmapped, false),
					 addFunc("runPRICEOnPathWeaverRegionsAndUnmapped", runPRICEOnPathWeaverRegionsAndUnmapped, false),
					 addFunc("runSavageOnPathWeaverRegionsAndUnmapped", runSavageOnPathWeaverRegionsAndUnmapped, false),
					 addFunc("runPolyteOnPathWeaverRegionsAndUnmapped", runPolyteOnPathWeaverRegionsAndUnmapped, false),

					 addFunc("runVelvetOptimizerAndMetaVelvetOnPathWeaverRegionsAndUnmapped", runVelvetOptimizerAndMetaVelvetOnPathWeaverRegionsAndUnmapped, false),
					 addFunc("runTrinityOnPathWeaverRegionsAndUnmapped", runTrinityOnPathWeaverRegionsAndUnmapped, false),
					 addFunc("runIDBAUDOnPathWeaverRegionsAndUnmapped", runIDBAUDOnPathWeaverRegionsAndUnmapped, false),
					 addFunc("runRayOnPathWeaverRegionsAndUnmapped", runRayOnPathWeaverRegionsAndUnmapped, false),
					 addFunc("runMIRAOnPathWeaverRegionsAndUnmapped", runMIRAOnPathWeaverRegionsAndUnmapped, false),
					 addFunc("runFermiLiteOnPathWeaverRegionsAndUnmapped", runFermiLiteOnPathWeaverRegionsAndUnmapped, false),
           },//,
          "programWrappersAssembleOnPathWeaverRunner") {

}


int programWrappersAssembleOnPathWeaverRunner::runMIRAOnPathWeaverRegions(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bedFile = "";
	bfs::path pwOutputDir = "";
	std::string sample;

	uint32_t MIRANumThreads = 1;
	std::string extraMIRAOptions;
	uint32_t reOrientingKmerLength = 9;
	uint32_t minFinalLength = 40;
	bfs::path MIRAOutDir = "MIRAOut";
	uint32_t numThreads = 1;
	uint32_t miraAttempts = 3 ;
	seqSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processVerbose();
	setUp.setOption(bedFile, "--bed", "The Regions to analyze", true);
	setUp.setOption(pwOutputDir, "--pwOutputDir", "The PathWeaver directory", true);
	setUp.setOption(sample, "--sample", "sample name", true);

	setUp.setOption(numThreads, "--numThreads", "num Threads");
	setUp.setOption(miraAttempts, "--miraAttempts", "mira Attempts");


	setUp.setOption(MIRANumThreads, "--MIRANumThreads", "MIRA Num Threads");
	setUp.setOption(extraMIRAOptions, "--extraMIRAOptions", "extra MIRA Options");

	setUp.setOption(minFinalLength, "--minFinalLength", "min Final Length");
	setUp.setOption(reOrientingKmerLength, "--reOrientingKmerLength", "re-orienting K-mer Length");

	setUp.setOption(MIRAOutDir,     "--MIRAOutDir",     "MIRA Out Directory name, will be relative to final pass directory");


	setUp.processDirectoryOutputName(njh::pasteAsStr(bfs::basename(pwOutputDir), "_MIRA_TODAY"), true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	njh::sys::requireExternalProgramThrow("mira");
	njh::sys::requireExternalProgramThrow("miraconvert");

	auto inputRegions = gatherRegions(bedFile.string(), "", setUp.pars_.verbose_);
	sortGRegionsByStart(inputRegions);

	std::set<std::string> regionNames;
	for(const auto & reg : inputRegions){
		regionNames.emplace(reg.uid_);
	}
	//njh::sort(regionNames);
	njh::concurrent::LockableQueue<std::string> regionsQueue(regionNames);

	bfs::path finalDirectory = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar("final"));
	bfs::path partialDirectory = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar("partial"));
	auto allFinalSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(finalDirectory, "allFinal.fasta"));
	auto allPartialSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(partialDirectory, "allPartial.fasta"));
	SeqOutput allFinalWriter(allFinalSeqOpts);
	SeqOutput allPartialWriter(allPartialSeqOpts);
	allFinalWriter.openOut();
	allPartialWriter.openOut();
	std::mutex allFinalWriterMut;
	std::mutex allPartialWriterMut;


	std::unordered_map<std::string,
			std::vector<std::shared_ptr<BamRegionInvestigator::RegionInfo>> >regInfosByUID;
	for (const auto & reg : inputRegions) {
		regInfosByUID[reg.uid_].emplace_back(std::make_shared<BamRegionInvestigator::RegionInfo>(reg));
	}

	std::unordered_map<std::string, std::string> exceptions;
	std::mutex exceptionsMut;

	std::function<void()> runMIRAOnRegion = [&](){
		std::string regionUid;
		while(regionsQueue.getVal(regionUid)){
			const auto & regInfo = njh::mapAt(regInfosByUID, regionUid);
			auto regionOutputDir = njh::files::make_path(setUp.pars_.directoryName_, regionUid, sample);
			njh::files::makeDirP(njh::files::MkdirPar{regionOutputDir});
			bfs::path refFnp = njh::files::make_path(pwOutputDir, regionUid, "allRefs.fasta");
			//first extract the reads
			bfs::path extractBam = njh::files::make_path(pwOutputDir, regionUid, sample + "_extraction", "extracted.bam");
			OutOptions outOpts(njh::files::make_path(regionOutputDir, "extracted"));
			auto readCounts = rawWriteExtractReadsFromBamOnlyMapped(extractBam, outOpts);
			bfs::path pairedR1 = njh::files::make_path(regionOutputDir, "extracted_R1.fastq");
			bfs::path pairedR2 = njh::files::make_path(regionOutputDir, "extracted_R2.fastq");
			bfs::path singles =  njh::files::make_path(regionOutputDir, "extracted.fastq");
			for(auto & reg : regInfo){
				reg->totalPairedReads_ = readCounts.pairedReads_;
				reg->totalReads_ = readCounts.pairedReads_ + readCounts.unpaiedReads_ + readCounts.orphans_;
				reg->totalFinalReads_ = readCounts.pairedReads_ + readCounts.unpaiedReads_ + readCounts.orphans_;
			}
			try {
				if(!exists(pairedR1) && !exists(singles)){
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", couldn't find " << pairedR1 << " or " << singles << ", need to have at least one of them" << "\n";
					throw std::runtime_error{ss.str()};
				}
				std::stringstream MIRACmdStream;
				MIRACmdStream << "cd " << regionOutputDir;
				MIRACmdStream << " && mira ";

				{
					OutputStream miramanifestOutput(njh::files::make_path(regionOutputDir, "mira_manifest.txt"));
					miramanifestOutput << "project = " << MIRAOutDir.filename().string() << std::endl;
					miramanifestOutput << "job = genome,denovo,accurate"<< std::endl;


					miramanifestOutput << "parameters = -CO:force_nonIUPACconsensus_perseqtype=yes -GENERAL:number_of_threads=" << MIRANumThreads << " COMMON_SETTINGS -NW:cmrnl=no -NW:cac=warn -NW:csrn=no -NW:cdrn=no"<< std::endl;
					//-EDIT:edit_homopolymer_overcalls=yes
					if(exists(pairedR1)){
						if(!exists(pairedR2)){
							std::stringstream ss;
							ss << __PRETTY_FUNCTION__ << ", found: " << pairedR1 << " but cound't find it's mate file: " << pairedR2 << "\n";
							throw std::runtime_error{ss.str()};
						}
						bfs::path pairedR1_app = njh::files::make_path(regionOutputDir, "appended_extracted_R1.fastq");
						bfs::path pairedR2_app = njh::files::make_path(regionOutputDir, "appended_extracted_R2.fastq");
						SeqOutput pairedR1_app_writer(SeqIOOptions::genFastqOut(pairedR1_app));
						pairedR1_app_writer.openOut();
						SeqOutput pairedR2_app_writer(SeqIOOptions::genFastqOut(pairedR2_app));
						pairedR2_app_writer.openOut();
						seqInfo seq;
						SeqInput pairedR1_reading(SeqIOOptions::genFastqIn(pairedR1));
						pairedR1_reading.openIn();
						while(pairedR1_reading.readNextRead(seq)){
							seq.name_.append("/1");
							pairedR1_app_writer.write(seq);
						}
						SeqInput pairedR2_reading(SeqIOOptions::genFastqIn(pairedR2));
						pairedR2_reading.openIn();
						while(pairedR2_reading.readNextRead(seq)){
							seq.name_.append("/2");
							pairedR2_app_writer.write(seq);
						}

						miramanifestOutput << "readgroup = " << sample << "--" << regionUid << std::endl;
						miramanifestOutput << "autopairing"<< std::endl;

						miramanifestOutput << "data = appended_extracted_R1.fastq appended_extracted_R2.fastq"<< std::endl;
						miramanifestOutput << "technology = solexa"<< std::endl;
						miramanifestOutput << "template_size = 50 1000 autorefine"<< std::endl;
						miramanifestOutput << "segment_placement = ---> <---"<< std::endl;
					}
					if(exists(singles)){
						miramanifestOutput << "readgroup = " << sample << "--" << regionUid << "-single" << std::endl;
						miramanifestOutput << "data = extracted.fastq"<< std::endl;
						miramanifestOutput << "technology = solexa"<< std::endl;
					}
				}

				MIRACmdStream  << " mira_manifest.txt "
												<< " >> MIRARunLog_" << njh::getCurrentDate() << ".txt 2>&1";
				auto MIRAFullOutputDir = njh::files::make_path(regionOutputDir, MIRAOutDir.filename().string() + "_assembly");

				auto MIRARunOutput = njh::sys::run({MIRACmdStream.str()});
				uint32_t currentMiraAttempt = 1;
				while(!MIRARunOutput.success_ && currentMiraAttempt <= miraAttempts){
					MIRARunOutput = njh::sys::run({MIRACmdStream.str()});

					++currentMiraAttempt;
				}
				BioCmdsUtils::checkRunOutThrow(MIRARunOutput, __PRETTY_FUNCTION__);

				OutOptions MIRARunOutputLogOpts(njh::files::make_path(MIRAFullOutputDir, "MIRARunOutput.json"));
				OutputStream MIRARunOutputLogOut(MIRARunOutputLogOpts);
				MIRARunOutputLogOut << njh::json::toJson(MIRARunOutput) << std::endl;

				auto contigsFnp = njh::files::make_path(regionOutputDir, MIRAOutDir.filename().string() + "_assembly/", MIRAOutDir.filename().string() + "_d_results", MIRAOutDir.filename().string() + "_out.unpadded.fasta");


				auto contigsSeqIoOpts = SeqIOOptions::genFastaIn(contigsFnp);
//				contigsSeqIoOpts.includeWhiteSpaceInName_ = false;
				contigsSeqIoOpts.lowerCaseBases_ = "upper";
				SeqInput contigsReader(contigsSeqIoOpts);
				auto contigsSeqs = contigsReader.readAllReads<seqInfo>();
				std::vector<std::shared_ptr<seqWithKmerInfo>> contigsKmerReads;
				contigsKmerReads.reserve(contigsSeqs.size());
				for (const auto & seq : contigsSeqs) {
					contigsKmerReads.emplace_back(std::make_shared<seqWithKmerInfo>(seq));
				}
				allSetKmers(contigsKmerReads, reOrientingKmerLength, true);

				SeqInput refReader(SeqIOOptions::genFastaIn(refFnp));
				auto refSeqs = refReader.readAllReads<seqInfo>();
				std::vector<std::shared_ptr<seqWithKmerInfo>> refKmerReads;
				std::vector<seqInfo> revComp_refSeqs;
				std::vector<kmerInfo> revComp_refSeqsKInfos;
				refKmerReads.reserve(refSeqs.size());
				for (const auto & seq : refSeqs) {
					refKmerReads.emplace_back(std::make_shared<seqWithKmerInfo>(seq));
				}
				allSetKmers(refKmerReads, reOrientingKmerLength, true);
				for(const auto & rSeq : refSeqs){
					revComp_refSeqs.emplace_back(rSeq);
					revComp_refSeqs.back().reverseComplementRead(false, true);
					revComp_refSeqsKInfos.emplace_back(revComp_refSeqs.back().seq_, 7, false);
				}
				for(const auto & seqKmer : contigsKmerReads) {
					uint32_t forwardWinners = 0;
					uint32_t revWinners = 0;
					for (const auto & refSeq : refKmerReads) {
						auto forDist = refSeq->compareKmers(*seqKmer);
						auto revDist = refSeq->compareKmersRevComp(*seqKmer);
						if (forDist.first < revDist.first) {
							++revWinners;
						} else {
							++forwardWinners;
						}
					}
					if (revWinners > forwardWinners) {
						seqKmer->seqBase_.reverseComplementRead(true, true);
					}
				}

				//sort by sequence length;
				njh::sort(contigsKmerReads, [](const std::shared_ptr<seqWithKmerInfo> & seq1, const std::shared_ptr<seqWithKmerInfo> & seq2){
					return len(seq1->seqBase_) > len(seq2->seqBase_);
				});

				OutOptions contigInfoOpts(njh::files::make_path(MIRAFullOutputDir, "contigs_outputInfo.tab.txt"));
				OutputStream contigInfoOut(contigInfoOpts);
				contigInfoOut << "name\tlength\tcoverage" << std::endl;

				for(const auto & contigsKmerRead : contigsKmerReads){
					auto assembleInfo = MIRAAssembleNameInfo(contigsKmerRead->seqBase_.name_);

					contigInfoOut << contigsKmerRead->seqBase_.name_
							<< "\t" << len(contigsKmerRead->seqBase_)
							<< "\t" << assembleInfo.coverage_ << std::endl;
				}
				auto reOrientedContigsFnp = njh::files::make_path(MIRAFullOutputDir, "reOriented_contigs.fasta");

				SeqOutput::write(contigsKmerReads, SeqIOOptions::genFastaOut(reOrientedContigsFnp));

				uint64_t maxLen = 0;
				readVec::getMaxLength(refSeqs, maxLen);
				readVec::getMaxLength(contigsKmerReads, maxLen);
				aligner alignerObj(maxLen, gapScoringParameters(5,1,0,0,0,0), substituteMatrix(2,-2), false);
				//alignerObj.processAlnInfoInputNoCheck(njh::files::make_path(resultsDirectory, "trimAlnCache").string(), setUp.pars_.verbose_);
				std::vector<kmerInfo> refSeqsKmerInfos;
				refSeqsKmerInfos.reserve(refSeqs.size());
				for(const auto & input : refSeqs){
					refSeqsKmerInfos.emplace_back(input.seq_, 7, false);
				}
				//readVecTrimmer::trimSeqToRefByGlobalAln(contigsKmerReads, refSeqs, refSeqsKmerInfos, alignerObj);
				std::vector<std::shared_ptr<seqWithKmerInfo>> finalSeqs;
				std::unordered_map<std::string, uint32_t> finalSeqCounts;
				std::unordered_map<std::string, uint32_t> finalSeqCountsWritten;
				for(const auto & seq : contigsKmerReads){
					auto trimmed = readVecTrimmer::trimSeqToRefByGlobalAlnBestNoOverlapIncludeRevComp(seq, refSeqs, revComp_refSeqs, refSeqsKmerInfos, revComp_refSeqsKInfos, alignerObj, false);
					for(auto & trimmedSeq : trimmed){
						bool found = false;
						for(const auto & finalSeq : finalSeqs){
							if(finalSeq->seqBase_.seq_ == trimmedSeq.seq_){
								found = true;
								break;
							}
						}
						if(!found){
							finalSeqs.emplace_back(std::make_shared<seqWithKmerInfo>(trimmedSeq, 7, false));
							++finalSeqCounts[trimmedSeq.name_];
						}
					}
				}
				double totalCoverage = 0;
				for(auto & seq : finalSeqs){
					auto assembleInfo = MIRAAssembleNameInfo(seq->seqBase_.name_);
					totalCoverage += assembleInfo.coverage_;
				}

				for(auto & seq : finalSeqs){
					auto assembleInfo = MIRAAssembleNameInfo(seq->seqBase_.name_);
					MetaDataInName seqMeta;
					seqMeta.addMeta("trimmedLength", len(seq->seqBase_));
					seqMeta.addMeta("estimatedPerBaseCoverage", assembleInfo.coverage_);
					seqMeta.addMeta("trimStatus", seq->seqBase_.on_);
					seqMeta.addMeta("regionUID", regionUid);
					seqMeta.addMeta("sample", sample);
					if(finalSeqCounts[seq->seqBase_.name_] > 1){
						seqMeta.addMeta("seqTrimmedCount", finalSeqCountsWritten[seq->seqBase_.name_]);
						++finalSeqCountsWritten[seq->seqBase_.name_];
					}
					seqMeta.resetMetaInName(seq->seqBase_.name_);
					seq->seqBase_.cnt_ = assembleInfo.seqNumber_;
					seq->seqBase_.name_ += njh::pasteAsStr("_t", assembleInfo.seqNumber_);
				}

				OutOptions trimmedContigInfoOpts(njh::files::make_path(MIRAFullOutputDir, "trimmed_reOriented_contigs_outputInfo.tab.txt"));
				OutputStream trimmedContigInfoOut(trimmedContigInfoOpts);
				trimmedContigInfoOut << "name\tlength\tcoverage" << std::endl;
				auto trimmedReOrientedContigsFnp = njh::files::make_path(MIRAFullOutputDir, "trimmed_reOriented_contigs.fasta");
				SeqOutput outputWriter(SeqIOOptions::genFastaOut(trimmedReOrientedContigsFnp));
				auto trimmedReOrientedContigsFnp_belowCutOff = njh::files::make_path(MIRAFullOutputDir, "trimmed_reOriented_contigs_belowCutOff.fasta");
				SeqOutput belowCutOffOutputWriter(SeqIOOptions::genFastaOut(trimmedReOrientedContigsFnp_belowCutOff));

				uint32_t belowCutOff = 0;
				uint32_t aboveCutOff = 0;
				bool allPassTrim = true;
				for (const auto & contigsKmerRead : finalSeqs) {
					if (len(contigsKmerRead->seqBase_) < minFinalLength) {
						++belowCutOff;
						belowCutOffOutputWriter.openWrite(contigsKmerRead);
						contigsKmerRead->seqBase_.on_ = false;
					} else {
						MetaDataInName seqMeta(contigsKmerRead->seqBase_.name_);
						trimmedContigInfoOut << contigsKmerRead->seqBase_.name_
								<< "\t" << len(contigsKmerRead->seqBase_)
								<< "\t" << seqMeta.getMeta("estimatedPerBaseCoverage")
								<< std::endl;
						if(!contigsKmerRead->seqBase_.on_){
							allPassTrim = false;
						}else{
							++aboveCutOff;
						}
						outputWriter.openWrite(contigsKmerRead);
					}
				}
				if(allPassTrim){
					std::lock_guard<std::mutex> lock(allFinalWriterMut);
					for (const auto & contigsKmerRead : finalSeqs) {
						if (len(contigsKmerRead->seqBase_) >= minFinalLength) {
							allFinalWriter.write(contigsKmerRead);
						}
					}
					for(auto & reg : regInfo){
						reg->infoCalled_ = true;
						reg->uniqHaps_ = aboveCutOff;
					}
				}else{
					std::lock_guard<std::mutex> lock(allPartialWriterMut);
					for (const auto & contigsKmerRead : finalSeqs) {
						if (len(contigsKmerRead->seqBase_) >= minFinalLength) {
							allPartialWriter.write(contigsKmerRead);
						}
					}
					for(auto & reg : regInfo){
						reg->infoCalled_ = false;
						reg->uniqHaps_ = 0;
					}
				}
			} catch (std::exception & e) {
				std::lock_guard<std::mutex> lock(exceptionsMut);
				exceptions[regionUid] = e.what();
				for(auto & reg : regInfo){
					reg->infoCalled_ = false;
					reg->uniqHaps_ = 0;
				}
			}
		}
	};


	njh::concurrent::runVoidFunctionThreaded(runMIRAOnRegion, numThreads);
	allFinalWriter.closeOut();
	allPartialWriter.closeOut();
	//sample,readTotal,readTotalUsed, success, name
	//
	OutputStream basicInfo(njh::files::make_path(finalDirectory, "basicInfoPerRegion.tab.txt"));

	basicInfo << "#chrom\tstart\tend\tname\tlength\tstrand\tsuccess\tuniqHaps\treadTotal\treadTotalUsed\ttotalPairedReads";
	basicInfo << "\tsample";
	uint32_t maxExtraFields = 0;
	for(const auto & p : inputRegions){

		auto bedOut = p.genBedRecordCore();
		if(bedOut.extraFields_.size() > maxExtraFields){
			maxExtraFields = bedOut.extraFields_.size();
		}
	}
	for(uint32_t t = 0; t < maxExtraFields; ++t){
		basicInfo << "\textraField"<<t;
	}
	basicInfo << "\n";

	std::map<uint32_t, uint32_t> coiCounts;

	for (const auto & reg : inputRegions) {
		const auto & regInfos = njh::mapAt(regInfosByUID, reg.uid_);
		++coiCounts[regInfos.front()->uniqHaps_];
		for(auto & regInfo : regInfos){
			auto bedOut = regInfo->region_.genBedRecordCore();
			basicInfo << bedOut.toDelimStr();
			basicInfo << "\t" << njh::boolToStr(regInfo->infoCalled_)
								<< "\t" << regInfo->uniqHaps_
								<< "\t" << regInfo->totalReads_
								<< "\t" << regInfo->totalFinalReads_
								<< "\t" << regInfo->totalPairedReads_
								<< "\t" << sample;
			for(const auto & extra : bedOut.extraFields_){
				basicInfo << "\t" << extra;
			}
			basicInfo << std::endl;
		}
	}

	OutputStream coiOut(njh::files::make_path(finalDirectory, "coiCounts.tab.txt"));
	coiOut << "coi\tcount" << std::endl;
	for(const auto & count : coiCounts){
		coiOut << count.first << "\t" << count.second << std::endl;
	}

	OutputStream exceptionsOut(njh::files::make_path(finalDirectory, "exceptionsMessages.tab.txt"));
	exceptionsOut << "regionUID\tmessage" << std::endl;
	for(const auto & exp : exceptions){
		exceptionsOut << exp.first << "\t" << exp.second << std::endl;
	}
	return 0;
}


int programWrappersAssembleOnPathWeaverRunner::runFermiLiteOnPathWeaverRegions(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bedFile = "";
	bfs::path pwOutputDir = "";
	std::string sample;

	uint32_t fermiLiteNumThreads = 1;
	std::string extraFermiLiteOptions;
	uint32_t reOrientingKmerLength = 9;
	uint32_t minFinalLength = 40;
	uint32_t numThreads = 1;
	seqSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processVerbose();
	setUp.setOption(bedFile, "--bed", "The Regions to analyze", true);
	setUp.setOption(pwOutputDir, "--pwOutputDir", "The PathWeaver directory", true);
	setUp.setOption(sample, "--sample", "sample name", true);

	setUp.setOption(numThreads, "--numThreads", "num Threads");


	setUp.setOption(fermiLiteNumThreads, "--fermiLiteNumThreads", "fermiLite Num Threads");
	setUp.setOption(extraFermiLiteOptions, "--extraFermiLiteOptions", "extra FermiLite Options");

	setUp.setOption(minFinalLength, "--minFinalLength", "min Final Length");
	setUp.setOption(reOrientingKmerLength, "--reOrientingKmerLength", "re-orienting K-mer Length");

	setUp.processDirectoryOutputName(njh::pasteAsStr(bfs::basename(pwOutputDir), "_fermiLite_TODAY"), true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	njh::sys::requireExternalProgramThrow("fml-asm");

	auto inputRegions = gatherRegions(bedFile.string(), "", setUp.pars_.verbose_);
	sortGRegionsByStart(inputRegions);

	std::set<std::string> regionNames;
	for(const auto & reg : inputRegions){
		regionNames.emplace(reg.uid_);
	}
	//njh::sort(regionNames);
	njh::concurrent::LockableQueue<std::string> regionsQueue(regionNames);

	bfs::path finalDirectory = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar("final"));
	bfs::path partialDirectory = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar("partial"));
	auto allFinalSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(finalDirectory, "allFinal.fasta"));
	auto allPartialSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(partialDirectory, "allPartial.fasta"));
	SeqOutput allFinalWriter(allFinalSeqOpts);
	SeqOutput allPartialWriter(allPartialSeqOpts);
	allFinalWriter.openOut();
	allPartialWriter.openOut();
	std::mutex allFinalWriterMut;
	std::mutex allPartialWriterMut;


	std::unordered_map<std::string,
					std::vector<std::shared_ptr<BamRegionInvestigator::RegionInfo>> >regInfosByUID;
	for (const auto & reg : inputRegions) {
		regInfosByUID[reg.uid_].emplace_back(std::make_shared<BamRegionInvestigator::RegionInfo>(reg));
	}

	std::unordered_map<std::string, std::string> exceptions;
	std::mutex exceptionsMut;

	std::function<void()> runFermiLiteOnRegion = [&](){
		std::string regionUid;
		while(regionsQueue.getVal(regionUid)){
			const auto & regInfo = njh::mapAt(regInfosByUID, regionUid);
			auto regionOutputDir = njh::files::make_path(setUp.pars_.directoryName_, regionUid, sample);
			njh::files::makeDirP(njh::files::MkdirPar{regionOutputDir});

			//

			bfs::path refFnp = njh::files::make_path(pwOutputDir, regionUid, "allRefs.fasta");
			//first extract the reads
			bfs::path pwRunDir = njh::files::make_path(pwOutputDir, regionUid, sample);
			bfs::path filtStiched_pairedR1Fnp_ = njh::files::make_path(pwRunDir, "filteredExtractedPairs_R1.fastq");
			bfs::path filtStiched_pairedR2Fnp_ = njh::files::make_path(pwRunDir, "filteredExtractedPairs_R2.fastq");
			bfs::path filtStiched_singlesFnp_ =  njh::files::make_path(pwRunDir, "filteredSingles.fastq");

			try {
				//concatenate into 1 file
				std::vector<bfs::path> filesToCollapse;
				if(bfs::exists(filtStiched_pairedR1Fnp_)){
					filesToCollapse.emplace_back(filtStiched_pairedR1Fnp_);
					filesToCollapse.emplace_back(filtStiched_pairedR2Fnp_);
				}
				if(bfs::exists(filtStiched_singlesFnp_)){
					filesToCollapse.emplace_back(filtStiched_singlesFnp_);
				}
				auto inputFnp = njh::files::make_path(regionOutputDir, "input.fastq.gz");
				auto outputFnp = njh::files::make_path(regionOutputDir, "raw_output.fastq");
				concatenateFiles(filesToCollapse, OutOptions(inputFnp));
				uint64_t totalReads = 0;
				double medianReadLength = 0;
				{
					std::vector<uint32_t> readLens;
					seqInfo inputSeq;
					SeqInput reader(SeqIOOptions::genFastqIn(inputFnp));
					reader.openIn();
					while(reader.readNextRead(inputSeq)){
						readLens.emplace_back(len(inputSeq));
						++totalReads;
					}
					medianReadLength = vectorMedianRef(readLens);
				}
				uint64_t pairedReads = countSeqs(SeqIOOptions::genFastqIn(filtStiched_pairedR1Fnp_), false);

				for(auto & reg : regInfo){
					reg->totalPairedReads_ = pairedReads;
					reg->totalReads_ = totalReads;
					reg->totalFinalReads_ = totalReads;
				}
				if(!exists(filtStiched_pairedR1Fnp_) && !exists(filtStiched_singlesFnp_)){
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", couldn't find " << filtStiched_pairedR1Fnp_ << " or " << filtStiched_singlesFnp_ << ", need to have at least one of them" << "\n";
					throw std::runtime_error{ss.str()};
				}

				std::stringstream raw_fermiLiteCmdStream;
				raw_fermiLiteCmdStream << "cd " << regionOutputDir;
				raw_fermiLiteCmdStream << " && fml-asm ";

				raw_fermiLiteCmdStream  << " -t " << fermiLiteNumThreads
																<< " " << extraFermiLiteOptions
																<< " input.fastq.gz "
																<< " > " << "raw_output.fastq";
				std::string raw_fermiLiteCmd = raw_fermiLiteCmdStream.str();
				std::stringstream fermiLiteCmdStream;
				fermiLiteCmdStream << raw_fermiLiteCmd << " 2> fermiLiteRunLog_" << njh::getCurrentDate() << ".txt";
				const auto & fermiLiteFullOutputDir = regionOutputDir;

				auto fermiLiteRunOutput = njh::sys::run({fermiLiteCmdStream.str()});

				OutOptions fermiLiteRunOutputLogOpts(njh::files::make_path(fermiLiteFullOutputDir, "fermiLiteRunOutput.json"));
				OutputStream fermiLiteRunOutputLogOut(fermiLiteRunOutputLogOpts);
				fermiLiteRunOutputLogOut << njh::json::toJson(fermiLiteRunOutput) << std::endl;

				auto contigsFnp = njh::files::make_path(fermiLiteFullOutputDir, "raw_output.fastq");

				auto contigsSeqIoOpts = SeqIOOptions::genFastqIn(contigsFnp);
//				contigsSeqIoOpts.includeWhiteSpaceInName_ = false;
				contigsSeqIoOpts.lowerCaseBases_ = "upper";
				SeqInput contigsReader(contigsSeqIoOpts);
				auto contigsSeqs = contigsReader.readAllReads<seqInfo>();
				std::vector<std::shared_ptr<seqWithKmerInfo>> contigsKmerReads;
				contigsKmerReads.reserve(contigsSeqs.size());
				for (const auto & seq : contigsSeqs) {
					contigsKmerReads.emplace_back(std::make_shared<seqWithKmerInfo>(seq));
				}
				allSetKmers(contigsKmerReads, reOrientingKmerLength, true);

				RefSeqsWithKmers refSeqs(refFnp, reOrientingKmerLength);

				readVec::reorientSeqs(contigsKmerReads, refSeqs.refKmerReads_);
				//sort by sequence length;
				njh::sort(contigsKmerReads, [](const std::shared_ptr<seqWithKmerInfo> & seq1, const std::shared_ptr<seqWithKmerInfo> & seq2){
					return len(seq1->seqBase_) > len(seq2->seqBase_);
				});

				OutOptions contigInfoOpts(njh::files::make_path(fermiLiteFullOutputDir, "contigs_outputInfo.tab.txt"));
				OutputStream contigInfoOut(contigInfoOpts);
				contigInfoOut << "name\tlength\tcoverage" << std::endl;

				for(const auto & contigsKmerRead : contigsKmerReads){
					auto assembleInfo = FermiLiteNameParse(contigsKmerRead->seqBase_.name_);
					contigsKmerRead->seqBase_.name_ = assembleInfo.modFullname_; //get rid of the \t characters in the name
					assembleInfo.coverage_ = (assembleInfo.coverage_ * medianReadLength)/len(contigsKmerRead->seqBase_);
					contigInfoOut << contigsKmerRead->seqBase_.name_
												<< "\t" << len(contigsKmerRead->seqBase_)
												<< "\t" << assembleInfo.coverage_ << std::endl;
				}
				auto reOrientedContigsFnp = njh::files::make_path(fermiLiteFullOutputDir, "reOriented_contigs.fasta");

				SeqOutput::write(contigsKmerReads, SeqIOOptions::genFastaOut(reOrientedContigsFnp));

				std::vector<std::shared_ptr<seqWithKmerInfo>> finalSeqs = trimToFinalSeqs(contigsKmerReads, refSeqs);
				std::unordered_map<std::string, uint32_t> finalSeqCounts;
				for(const auto & seq : finalSeqs){
					++finalSeqCounts[seq->seqBase_.name_];
				}

				double totalCoverage = 0;
				for(auto & seq : finalSeqs){
					auto assembleInfo = FermiLiteNameParse(seq->seqBase_.name_);
					assembleInfo.coverage_ = (assembleInfo.coverage_ * medianReadLength)/len(seq->seqBase_);
					totalCoverage += assembleInfo.coverage_;
				}
				std::unordered_map<std::string, uint32_t> finalSeqCountsWritten;
				for(auto & seq : finalSeqs){
					auto assembleInfo = FermiLiteNameParse(seq->seqBase_.name_);
					assembleInfo.coverage_ = (assembleInfo.coverage_ * medianReadLength)/len(seq->seqBase_);
					MetaDataInName seqMeta;
					seqMeta.addMeta("trimmedLength", len(seq->seqBase_));
					seqMeta.addMeta("estimatedPerBaseCoverage", assembleInfo.coverage_);
					seqMeta.addMeta("trimStatus", seq->seqBase_.on_);
					seqMeta.addMeta("regionUID", regionUid);
					seqMeta.addMeta("sample", sample);
					if(finalSeqCounts[seq->seqBase_.name_] > 1){
						seqMeta.addMeta("seqTrimmedCount", finalSeqCountsWritten[seq->seqBase_.name_]);
						++finalSeqCountsWritten[seq->seqBase_.name_];
					}
					seqMeta.resetMetaInName(seq->seqBase_.name_);
					seq->seqBase_.cnt_ = (assembleInfo.coverage_/totalCoverage) * (totalReads);
					seq->seqBase_.name_ += njh::pasteAsStr("_t", seq->seqBase_.cnt_);
				}

				OutOptions trimmedContigInfoOpts(njh::files::make_path(fermiLiteFullOutputDir, "trimmed_reOriented_contigs_outputInfo.tab.txt"));
				OutputStream trimmedContigInfoOut(trimmedContigInfoOpts);
				trimmedContigInfoOut << "name\tlength\tcoverage" << std::endl;
				auto trimmedReOrientedContigsFnp = njh::files::make_path(fermiLiteFullOutputDir, "trimmed_reOriented_contigs.fasta");
				SeqOutput outputWriter(SeqIOOptions::genFastaOut(trimmedReOrientedContigsFnp));
				auto trimmedReOrientedContigsFnp_belowCutOff = njh::files::make_path(fermiLiteFullOutputDir, "trimmed_reOriented_contigs_belowCutOff.fasta");
				SeqOutput belowCutOffOutputWriter(SeqIOOptions::genFastaOut(trimmedReOrientedContigsFnp_belowCutOff));

				uint32_t belowCutOff = 0;
				uint32_t aboveCutOff = 0;
				bool allPassTrim = true;
				for (const auto & contigsKmerRead : finalSeqs) {
					if (len(contigsKmerRead->seqBase_) < minFinalLength) {
						++belowCutOff;
						belowCutOffOutputWriter.openWrite(contigsKmerRead);
						contigsKmerRead->seqBase_.on_ = false;
					} else {
						MetaDataInName seqMeta(contigsKmerRead->seqBase_.name_);
						trimmedContigInfoOut << contigsKmerRead->seqBase_.name_
																 << "\t" << len(contigsKmerRead->seqBase_)
																 << "\t" << seqMeta.getMeta("estimatedPerBaseCoverage")
																 << std::endl;
						if(!contigsKmerRead->seqBase_.on_){
							allPassTrim = false;
						}else{
							++aboveCutOff;
						}
						outputWriter.openWrite(contigsKmerRead);
					}
				}
				if(allPassTrim && !finalSeqs.empty()){
					std::lock_guard<std::mutex> lock(allFinalWriterMut);
					for (const auto & contigsKmerRead : finalSeqs) {
						if (len(contigsKmerRead->seqBase_) >= minFinalLength) {
							allFinalWriter.write(contigsKmerRead);
						}
					}
					for(auto & reg : regInfo){
						reg->infoCalled_ = true;
						reg->uniqHaps_ = aboveCutOff;
					}
				}else{
					for(auto & reg : regInfo){
						reg->infoCalled_ = false;
						reg->uniqHaps_ = 0;
					}
					std::lock_guard<std::mutex> lock(allPartialWriterMut);
					for (const auto & contigsKmerRead : finalSeqs) {
						if (len(contigsKmerRead->seqBase_) >= minFinalLength) {
							allPartialWriter.write(contigsKmerRead);
						}
					}
					for(auto & reg : regInfo){
						reg->infoCalled_ = false;
						reg->uniqHaps_ = 0;
					}
				}
			} catch (std::exception & e) {
				std::lock_guard<std::mutex> lock(exceptionsMut);
				exceptions[regionUid] = e.what();
				for(auto & reg : regInfo){
					reg->infoCalled_ = false;
					reg->uniqHaps_ = 0;
				}
			}
		}
	};


	njh::concurrent::runVoidFunctionThreaded(runFermiLiteOnRegion, numThreads);
	allFinalWriter.closeOut();
	allPartialWriter.closeOut();
	//sample,readTotal,readTotalUsed, success, name
	//
	OutputStream basicInfo(njh::files::make_path(finalDirectory, "basicInfoPerRegion.tab.txt"));

	basicInfo << "#chrom\tstart\tend\tname\tlength\tstrand\tsuccess\tuniqHaps\treadTotal\treadTotalUsed\ttotalPairedReads";
	basicInfo << "\tsample";
	uint32_t maxExtraFields = 0;
	for(const auto & p : inputRegions){

		auto bedOut = p.genBedRecordCore();
		if(bedOut.extraFields_.size() > maxExtraFields){
			maxExtraFields = bedOut.extraFields_.size();
		}
	}
	for(uint32_t t = 0; t < maxExtraFields; ++t){
		basicInfo << "\textraField"<<t;
	}
	basicInfo << "\n";

	std::map<uint32_t, uint32_t> coiCounts;

	for (const auto & reg : inputRegions) {
		const auto & regInfos = njh::mapAt(regInfosByUID, reg.uid_);
		++coiCounts[regInfos.front()->uniqHaps_];
		for(auto & regInfo : regInfos){
			auto bedOut = regInfo->region_.genBedRecordCore();
			basicInfo << bedOut.toDelimStr();
			basicInfo << "\t" << njh::boolToStr(regInfo->infoCalled_)
								<< "\t" << regInfo->uniqHaps_
								<< "\t" << regInfo->totalReads_
								<< "\t" << regInfo->totalFinalReads_
								<< "\t" << regInfo->totalPairedReads_
								<< "\t" << sample;
			for(const auto & extra : bedOut.extraFields_){
				basicInfo << "\t" << extra;
			}
			basicInfo << std::endl;
		}
	}

	OutputStream coiOut(njh::files::make_path(finalDirectory, "coiCounts.tab.txt"));
	coiOut << "coi\tcount" << std::endl;
	for(const auto & count : coiCounts){
		coiOut << count.first << "\t" << count.second << std::endl;
	}

	OutputStream exceptionsOut(njh::files::make_path(finalDirectory, "exceptionsMessages.tab.txt"));
	exceptionsOut << "regionUID\tmessage" << std::endl;
	for(const auto & exp : exceptions){
		exceptionsOut << exp.first << "\t" << exp.second << std::endl;
	}



	return 0;
}


//
int programWrappersAssembleOnPathWeaverRunner::runUnicyclerOnPathWeaverRegions(const njh::progutils::CmdArgs & inputCommands) {
		bfs::path bedFile = "";
		bfs::path pwOutputDir = "";
		std::string sample;

		uint32_t unicyclerNumThreads = 1;
		std::string extraUnicyclerOptions = "--no_pilon";
		uint32_t reOrientingKmerLength = 9;
		uint32_t minFinalLength = 40;
		uint32_t numThreads = 1;
		bfs::path unicyclerOutDir = "unicyclerOut";
		seqSetUp setUp(inputCommands);
		setUp.processDebug();
		setUp.processVerbose();
		setUp.setOption(bedFile, "--bed", "The Regions to analyze", true);
		setUp.setOption(pwOutputDir, "--pwOutputDir", "The PathWeaver directory", true);
		setUp.setOption(sample, "--sample", "sample name", true);

		setUp.setOption(numThreads, "--numThreads", "num Threads");


		setUp.setOption(unicyclerNumThreads, "--unicyclerNumThreads", "unicycler Num Threads");
		setUp.setOption(extraUnicyclerOptions, "--extraUnicyclerOptions", "extra Unicycler Options");

		setUp.setOption(minFinalLength, "--minFinalLength", "min Final Length");
		setUp.setOption(reOrientingKmerLength, "--reOrientingKmerLength", "re-orienting K-mer Length");
		setUp.setOption(unicyclerOutDir,     "--unicyclerOutDir",     "unicycler Out Directory name, will be relative to final pass directory");


		setUp.processDirectoryOutputName(njh::pasteAsStr(bfs::basename(pwOutputDir), "_unicycler_TODAY"), true);
		setUp.finishSetUp(std::cout);
		setUp.startARunLog(setUp.pars_.directoryName_);
		njh::sys::requireExternalProgramThrow("unicycler");

		auto inputRegions = gatherRegions(bedFile.string(), "", setUp.pars_.verbose_);
		sortGRegionsByStart(inputRegions);

		std::set<std::string> regionNames;
		for(const auto & reg : inputRegions){
			regionNames.emplace(reg.uid_);
		}
		//njh::sort(regionNames);
		njh::concurrent::LockableQueue<std::string> regionsQueue(regionNames);

		bfs::path finalDirectory = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar("final"));
		bfs::path partialDirectory = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar("partial"));
		auto allFinalSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(finalDirectory, "allFinal.fasta"));
		auto allPartialSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(partialDirectory, "allPartial.fasta"));
		SeqOutput allFinalWriter(allFinalSeqOpts);
		SeqOutput allPartialWriter(allPartialSeqOpts);
		allFinalWriter.openOut();
		allPartialWriter.openOut();
		std::mutex allFinalWriterMut;
		std::mutex allPartialWriterMut;


		std::unordered_map<std::string,
						std::vector<std::shared_ptr<BamRegionInvestigator::RegionInfo>> >regInfosByUID;
		for (const auto & reg : inputRegions) {
			regInfosByUID[reg.uid_].emplace_back(std::make_shared<BamRegionInvestigator::RegionInfo>(reg));
		}

		std::unordered_map<std::string, std::string> exceptions;
		std::mutex exceptionsMut;

		std::function<void()> runUnicyclerOnRegion = [&](){
			std::string regionUid;
			while(regionsQueue.getVal(regionUid)){
				const auto & regInfo = njh::mapAt(regInfosByUID, regionUid);
				auto regionOutputDir = njh::files::make_path(setUp.pars_.directoryName_, regionUid, sample);
				njh::files::makeDirP(njh::files::MkdirPar{regionOutputDir});

				//

				bfs::path refFnp = njh::files::make_path(pwOutputDir, regionUid, "allRefs.fasta");


				//first extract the reads
				bfs::path extractBam = njh::files::make_path(pwOutputDir, regionUid, sample + "_extraction", "extracted.bam");
				OutOptions outOpts(njh::files::make_path(regionOutputDir, "extracted"));
				auto readCounts = rawWriteExtractReadsFromBamOnlyMapped(extractBam, outOpts);
				bfs::path pairedR1 = njh::files::make_path(regionOutputDir, "extracted_R1.fastq");
				bfs::path pairedR2 = njh::files::make_path(regionOutputDir, "extracted_R2.fastq");
				bfs::path singles =  njh::files::make_path(regionOutputDir, "extracted.fastq");
				uint64_t totalReads = readCounts.pairedReads_ + readCounts.unpaiedReads_ + readCounts.orphans_;
				uint32_t minRegionSize = std::numeric_limits<uint32_t>::max();

				for(auto & reg : regInfo){
					reg->totalPairedReads_ = readCounts.pairedReads_;
					reg->totalReads_ = totalReads;
					reg->totalFinalReads_ = totalReads;
					if(reg->region_.getLen() < minRegionSize){
						minRegionSize = reg->region_.getLen();
					}
				}


				try {
					if(!exists(pairedR1) && !exists(singles)){
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ", couldn't find " << pairedR1 << " or " << singles << ", need to have at least one of them" << "\n";
						throw std::runtime_error{ss.str()};
					}
					std::stringstream raw_unicyclerCmdStream;
					raw_unicyclerCmdStream << "cd " << regionOutputDir;
					raw_unicyclerCmdStream << " && unicycler ";

					if(exists(pairedR1)){
						if(!exists(pairedR2)){
							std::stringstream ss;
							ss << __PRETTY_FUNCTION__ << ", found: " << pairedR1 << " but cound't find it's mate file: " << pairedR2 << "\n";
							throw std::runtime_error{ss.str()};
						}else{
							raw_unicyclerCmdStream << " -1 " << pairedR1.filename() << " -2 " << pairedR2.filename() << " ";
						}
					}
					if(exists(singles)){
						raw_unicyclerCmdStream << " -s  " << singles.filename();
					}
					raw_unicyclerCmdStream  << " -t " << unicyclerNumThreads
													 << " " << extraUnicyclerOptions
													 << " -o " << unicyclerOutDir << " ";
					if(minRegionSize < 1000){
						uint32_t minComponentSize = static_cast<uint32_t>(std::max(minRegionSize * .10, 1.0));
						raw_unicyclerCmdStream << " --min_fasta_length " << minComponentSize << " --min_component_size " << minComponentSize << " --min_dead_end_size " << minComponentSize << " ";
					}
					std::string raw_unicyclerCmd = raw_unicyclerCmdStream.str();
					std::stringstream unicyclerCmdStream;
					unicyclerCmdStream << raw_unicyclerCmd << " > unicyclerRunLog_" << njh::getCurrentDate() << ".txt 2>&1";
					auto unicyclerFullOutputDir = njh::files::make_path(regionOutputDir, unicyclerOutDir);

					auto unicyclerRunOutput = njh::sys::run({unicyclerCmdStream.str()});
					if(!unicyclerRunOutput.success_){
						std::stringstream unicyclerCmdStream_reAttempt;
						//sometimes unicycler fails at the spades correction step, so try re-runing without correcting
						unicyclerCmdStream_reAttempt << raw_unicyclerCmd << " --no_correct  > unicyclerReRunLog_" << njh::getCurrentDate() << ".txt 2>&1";
						if(bfs::exists(unicyclerFullOutputDir)){
							bfs::rename(unicyclerFullOutputDir, njh::files::make_path(regionOutputDir, "failed_" + unicyclerOutDir.string()));
						}
						auto unicyclerReRunOutput = njh::sys::run({unicyclerCmdStream_reAttempt.str()});
						if(!unicyclerReRunOutput.success_ && !unicyclerRunOutput.success_){
							//throw if they both failed
							std::stringstream ss;
							ss << __PRETTY_FUNCTION__ << ", error:\n" << unicyclerRunOutput.stdErr_ << "\n" << unicyclerReRunOutput.stdErr_ << "\n";
							throw std::runtime_error{ss.str()};
						}
					}


					OutOptions unicyclerRunOutputLogOpts(njh::files::make_path(unicyclerFullOutputDir, "unicyclerRunOutput.json"));
					OutputStream unicyclerRunOutputLogOut(unicyclerRunOutputLogOpts);
					unicyclerRunOutputLogOut << njh::json::toJson(unicyclerRunOutput) << std::endl;

					auto contigsFnp = njh::files::make_path(unicyclerFullOutputDir, "assembly.fasta");
					std::regex unicyclerNamePat{R"(^(\d+) length=(\d+) depth=([0-9.]+)x.*)"};

					auto contigsSeqIoOpts = SeqIOOptions::genFastaIn(contigsFnp);
//				contigsSeqIoOpts.includeWhiteSpaceInName_ = false;
					contigsSeqIoOpts.lowerCaseBases_ = "upper";
					SeqInput contigsReader(contigsSeqIoOpts);
					auto contigsSeqs = contigsReader.readAllReads<seqInfo>();
					std::vector<std::shared_ptr<seqWithKmerInfo>> contigsKmerReads;
					contigsKmerReads.reserve(contigsSeqs.size());
					for (const auto & seq : contigsSeqs) {
						contigsKmerReads.emplace_back(std::make_shared<seqWithKmerInfo>(seq));
					}
					allSetKmers(contigsKmerReads, reOrientingKmerLength, true);

					RefSeqsWithKmers refSeqs(refFnp, reOrientingKmerLength);

					readVec::reorientSeqs(contigsKmerReads, refSeqs.refKmerReads_);
					//sort by sequence length;
					njh::sort(contigsKmerReads, [](const std::shared_ptr<seqWithKmerInfo> & seq1, const std::shared_ptr<seqWithKmerInfo> & seq2){
						return len(seq1->seqBase_) > len(seq2->seqBase_);
					});

					OutOptions contigInfoOpts(njh::files::make_path(unicyclerFullOutputDir, "contigs_outputInfo.tab.txt"));
					OutputStream contigInfoOut(contigInfoOpts);
					contigInfoOut << "name\tlength\tcoverage" << std::endl;

					for(const auto & contigsKmerRead : contigsKmerReads){
						auto assembleInfo = DefaultAssembleNameInfo(contigsKmerRead->seqBase_.name_, unicyclerNamePat);
						assembleInfo.coverage_ = (assembleInfo.coverage_ * totalReads * readCounts.medianReadLength_)/assembleInfo.len_;
						contigInfoOut << contigsKmerRead->seqBase_.name_
													<< "\t" << len(contigsKmerRead->seqBase_)
													<< "\t" << assembleInfo.coverage_ << std::endl;
					}
					auto reOrientedContigsFnp = njh::files::make_path(unicyclerFullOutputDir, "reOriented_contigs.fasta");

					SeqOutput::write(contigsKmerReads, SeqIOOptions::genFastaOut(reOrientedContigsFnp));

					std::vector<std::shared_ptr<seqWithKmerInfo>> finalSeqs = trimToFinalSeqs(contigsKmerReads, refSeqs);
					std::unordered_map<std::string, uint32_t> finalSeqCounts;
					for(const auto & seq : finalSeqs){
						++finalSeqCounts[seq->seqBase_.name_];
					}

					double totalCoverage = 0;
					for(auto & seq : finalSeqs){
						auto assembleInfo = DefaultAssembleNameInfo(seq->seqBase_.name_, unicyclerNamePat);
						assembleInfo.coverage_ = (assembleInfo.coverage_ * totalReads * readCounts.medianReadLength_)/assembleInfo.len_;
						totalCoverage += assembleInfo.coverage_;
					}
					std::unordered_map<std::string, uint32_t> finalSeqCountsWritten;
					for(auto & seq : finalSeqs){
						auto assembleInfo = DefaultAssembleNameInfo(seq->seqBase_.name_, unicyclerNamePat);
						assembleInfo.coverage_ = (assembleInfo.coverage_ * totalReads * readCounts.medianReadLength_)/assembleInfo.len_;
						MetaDataInName seqMeta;
						seqMeta.addMeta("trimmedLength", len(seq->seqBase_));
						seqMeta.addMeta("estimatedPerBaseCoverage", assembleInfo.coverage_);
						seqMeta.addMeta("trimStatus", seq->seqBase_.on_);
						seqMeta.addMeta("regionUID", regionUid);
						seqMeta.addMeta("sample", sample);
						if(finalSeqCounts[seq->seqBase_.name_] > 1){
							seqMeta.addMeta("seqTrimmedCount", finalSeqCountsWritten[seq->seqBase_.name_]);
							++finalSeqCountsWritten[seq->seqBase_.name_];
						}
						seqMeta.resetMetaInName(seq->seqBase_.name_);
						seq->seqBase_.cnt_ = (assembleInfo.coverage_/totalCoverage) * (readCounts.pairedReads_ + readCounts.unpaiedReads_ + readCounts.orphans_);
						seq->seqBase_.name_ += njh::pasteAsStr("_t", seq->seqBase_.cnt_);
					}

					OutOptions trimmedContigInfoOpts(njh::files::make_path(unicyclerFullOutputDir, "trimmed_reOriented_contigs_outputInfo.tab.txt"));
					OutputStream trimmedContigInfoOut(trimmedContigInfoOpts);
					trimmedContigInfoOut << "name\tlength\tcoverage" << std::endl;
					auto trimmedReOrientedContigsFnp = njh::files::make_path(unicyclerFullOutputDir, "trimmed_reOriented_contigs.fasta");
					SeqOutput outputWriter(SeqIOOptions::genFastaOut(trimmedReOrientedContigsFnp));
					auto trimmedReOrientedContigsFnp_belowCutOff = njh::files::make_path(unicyclerFullOutputDir, "trimmed_reOriented_contigs_belowCutOff.fasta");
					SeqOutput belowCutOffOutputWriter(SeqIOOptions::genFastaOut(trimmedReOrientedContigsFnp_belowCutOff));

					uint32_t belowCutOff = 0;
					uint32_t aboveCutOff = 0;
					bool allPassTrim = true;
					for (const auto & contigsKmerRead : finalSeqs) {
						if (len(contigsKmerRead->seqBase_) < minFinalLength) {
							++belowCutOff;
							belowCutOffOutputWriter.openWrite(contigsKmerRead);
							contigsKmerRead->seqBase_.on_ = false;
						} else {
							MetaDataInName seqMeta(contigsKmerRead->seqBase_.name_);
							trimmedContigInfoOut << contigsKmerRead->seqBase_.name_
																	 << "\t" << len(contigsKmerRead->seqBase_)
																	 << "\t" << seqMeta.getMeta("estimatedPerBaseCoverage")
																	 << std::endl;
							if(!contigsKmerRead->seqBase_.on_){
								allPassTrim = false;
							}else{
								++aboveCutOff;
							}
							outputWriter.openWrite(contigsKmerRead);
						}
					}
					if(allPassTrim && !finalSeqs.empty()){
						std::lock_guard<std::mutex> lock(allFinalWriterMut);
						for (const auto & contigsKmerRead : finalSeqs) {
							if (len(contigsKmerRead->seqBase_) >= minFinalLength) {
								allFinalWriter.write(contigsKmerRead);
							}
						}
						for(auto & reg : regInfo){
							reg->infoCalled_ = true;
							reg->uniqHaps_ = aboveCutOff;
						}
					}else{
						std::lock_guard<std::mutex> lock(allPartialWriterMut);
						for (const auto & contigsKmerRead : finalSeqs) {
							if (len(contigsKmerRead->seqBase_) >= minFinalLength) {
								allPartialWriter.write(contigsKmerRead);
							}
						}
						for(auto & reg : regInfo){
							reg->infoCalled_ = false;
							reg->uniqHaps_ = 0;
						}
					}
				} catch (std::exception & e) {
					std::lock_guard<std::mutex> lock(exceptionsMut);
					exceptions[regionUid] = e.what();
					for(auto & reg : regInfo){
						reg->infoCalled_ = false;
						reg->uniqHaps_ = 0;
					}
				}
			}
		};


		njh::concurrent::runVoidFunctionThreaded(runUnicyclerOnRegion, numThreads);
		allFinalWriter.closeOut();
		allPartialWriter.closeOut();
		//sample,readTotal,readTotalUsed, success, name
		//
		OutputStream basicInfo(njh::files::make_path(finalDirectory, "basicInfoPerRegion.tab.txt"));

		basicInfo << "#chrom\tstart\tend\tname\tlength\tstrand\tsuccess\tuniqHaps\treadTotal\treadTotalUsed\ttotalPairedReads";
		basicInfo << "\tsample";
		uint32_t maxExtraFields = 0;
		for(const auto & p : inputRegions){

			auto bedOut = p.genBedRecordCore();
			if(bedOut.extraFields_.size() > maxExtraFields){
				maxExtraFields = bedOut.extraFields_.size();
			}
		}
		for(uint32_t t = 0; t < maxExtraFields; ++t){
			basicInfo << "\textraField"<<t;
		}
		basicInfo << "\n";

		std::map<uint32_t, uint32_t> coiCounts;

		for (const auto & reg : inputRegions) {
			const auto & regInfos = njh::mapAt(regInfosByUID, reg.uid_);
			++coiCounts[regInfos.front()->uniqHaps_];
			for(auto & regInfo : regInfos){
				auto bedOut = regInfo->region_.genBedRecordCore();
				basicInfo << bedOut.toDelimStr();
				basicInfo << "\t" << njh::boolToStr(regInfo->infoCalled_)
									<< "\t" << regInfo->uniqHaps_
									<< "\t" << regInfo->totalReads_
									<< "\t" << regInfo->totalFinalReads_
									<< "\t" << regInfo->totalPairedReads_
									<< "\t" << sample;
				for(const auto & extra : bedOut.extraFields_){
					basicInfo << "\t" << extra;
				}
				basicInfo << std::endl;
			}
		}

		OutputStream coiOut(njh::files::make_path(finalDirectory, "coiCounts.tab.txt"));
		coiOut << "coi\tcount" << std::endl;
		for(const auto & count : coiCounts){
			coiOut << count.first << "\t" << count.second << std::endl;
		}

		OutputStream exceptionsOut(njh::files::make_path(finalDirectory, "exceptionsMessages.tab.txt"));
		exceptionsOut << "regionUID\tmessage" << std::endl;
		for(const auto & exp : exceptions){
			exceptionsOut << exp.first << "\t" << exp.second << std::endl;
		}



		return 0;
	}

int programWrappersAssembleOnPathWeaverRunner::runSpadesOnPathWeaverRegions(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bedFile = "";
	bfs::path pwOutputDir = "";
	std::string sample;

	uint32_t spadesNumThreads = 1;
	std::string extraSpadesOptions;
	uint32_t reOrientingKmerLength = 9;
	uint32_t minFinalLength = 40;
	bool runMeta = false;
	bfs::path spadesOutDir = "spadesOut";
	uint32_t numThreads = 1;
	seqSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processVerbose();
	setUp.setOption(bedFile, "--bed", "The Regions to analyze", true);
	setUp.setOption(pwOutputDir, "--pwOutputDir", "The PathWeaver directory", true);
	setUp.setOption(sample, "--sample", "sample name", true);

	setUp.setOption(numThreads, "--numThreads", "num Threads");


	setUp.setOption(spadesNumThreads, "--spadesNumThreads", "spades Num Threads");
	setUp.setOption(extraSpadesOptions, "--extraSpadesOptions", "extra Spades Options");

	setUp.setOption(minFinalLength, "--minFinalLength", "min Final Length");
	setUp.setOption(runMeta, "--runMeta", "Run Meta");
	setUp.setOption(reOrientingKmerLength, "--reOrientingKmerLength", "re-orienting K-mer Length");
	if(runMeta){
		spadesOutDir = "metaspadesOut";
	}
	setUp.setOption(spadesOutDir,     "--spadesOutDir",     "spades Out Directory name, will be relative to final pass directory");


	setUp.processDirectoryOutputName(njh::pasteAsStr(bfs::basename(pwOutputDir), "_spades_TODAY"), true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	njh::sys::requireExternalProgramThrow("spades.py");

	auto inputRegions = gatherRegions(bedFile.string(), "", setUp.pars_.verbose_);
	sortGRegionsByStart(inputRegions);

	std::set<std::string> regionNames;
	for(const auto & reg : inputRegions){
		regionNames.emplace(reg.uid_);
	}
	//njh::sort(regionNames);
	njh::concurrent::LockableQueue<std::string> regionsQueue(regionNames);

	bfs::path finalDirectory = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar("final"));
	bfs::path partialDirectory = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar("partial"));
	auto allFinalSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(finalDirectory, "allFinal.fasta"));
	auto allPartialSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(partialDirectory, "allPartial.fasta"));
	SeqOutput allFinalWriter(allFinalSeqOpts);
	SeqOutput allPartialWriter(allPartialSeqOpts);
	allFinalWriter.openOut();
	allPartialWriter.openOut();
	std::mutex allFinalWriterMut;
	std::mutex allPartialWriterMut;


	std::unordered_map<std::string,
			std::vector<std::shared_ptr<BamRegionInvestigator::RegionInfo>> >regInfosByUID;
	for (const auto & reg : inputRegions) {
		regInfosByUID[reg.uid_].emplace_back(std::make_shared<BamRegionInvestigator::RegionInfo>(reg));
	}

	std::unordered_map<std::string, std::string> exceptions;
	std::mutex exceptionsMut;

	std::function<void()> runSpadesOnRegion = [&](){
		std::string regionUid;
		while(regionsQueue.getVal(regionUid)){
			const auto & regInfo = njh::mapAt(regInfosByUID, regionUid);
			auto regionOutputDir = njh::files::make_path(setUp.pars_.directoryName_, regionUid, sample);
			njh::files::makeDirP(njh::files::MkdirPar{regionOutputDir});

			//

			bfs::path refFnp = njh::files::make_path(pwOutputDir, regionUid, "allRefs.fasta");


			//first extract the reads
			bfs::path extractBam = njh::files::make_path(pwOutputDir, regionUid, sample + "_extraction", "extracted.bam");
			OutOptions outOpts(njh::files::make_path(regionOutputDir, "extracted"));
			auto readCounts = rawWriteExtractReadsFromBamOnlyMapped(extractBam, outOpts);
			bfs::path pairedR1 = njh::files::make_path(regionOutputDir, "extracted_R1.fastq");
			bfs::path pairedR2 = njh::files::make_path(regionOutputDir, "extracted_R2.fastq");
			bfs::path singles =  njh::files::make_path(regionOutputDir, "extracted.fastq");
			for(auto & reg : regInfo){
				reg->totalPairedReads_ = readCounts.pairedReads_;
				reg->totalReads_ = readCounts.pairedReads_ + readCounts.unpaiedReads_ + readCounts.orphans_;
				reg->totalFinalReads_ = readCounts.pairedReads_ + readCounts.unpaiedReads_ + readCounts.orphans_;
			}
			try {
				if(!exists(pairedR1) && !exists(singles)){
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", couldn't find " << pairedR1 << " or " << singles << ", need to have at least one of them" << "\n";
					throw std::runtime_error{ss.str()};
				}
				std::stringstream spadesCmdStream;
				spadesCmdStream << "cd " << regionOutputDir;
				if(runMeta){
					spadesCmdStream << " && metaspades.py ";
				}else{
					spadesCmdStream << " && spades.py ";
				}

				if(exists(pairedR1)){
					if(!exists(pairedR2)){
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ", found: " << pairedR1 << " but cound't find it's mate file: " << pairedR2 << "\n";
						throw std::runtime_error{ss.str()};
					}else{
						spadesCmdStream << " -1 " << pairedR1.filename() << " -2 " << pairedR2.filename() << " ";
					}
				}
				if(exists(singles)){
					spadesCmdStream << " -s  " << singles.filename();
				}
				spadesCmdStream  << " -t " << spadesNumThreads
												<< " " << extraSpadesOptions
												<< " -o " << spadesOutDir
												<< " > spadesRunLog_" << njh::getCurrentDate() << ".txt 2>&1";
				auto spadesFullOutputDir = njh::files::make_path(regionOutputDir, spadesOutDir);

				auto spadesRunOutput = njh::sys::run({spadesCmdStream.str()});

				BioCmdsUtils::checkRunOutThrow(spadesRunOutput, __PRETTY_FUNCTION__);

				OutOptions spadesRunOutputLogOpts(njh::files::make_path(spadesFullOutputDir, "spadesRunOutput.json"));
				OutputStream spadesRunOutputLogOut(spadesRunOutputLogOpts);
				spadesRunOutputLogOut << njh::json::toJson(spadesRunOutput) << std::endl;

				auto contigsFnp = njh::files::make_path(spadesFullOutputDir, "contigs.fasta");


				auto contigsSeqIoOpts = SeqIOOptions::genFastaIn(contigsFnp);
//				contigsSeqIoOpts.includeWhiteSpaceInName_ = false;
				contigsSeqIoOpts.lowerCaseBases_ = "upper";
				SeqInput contigsReader(contigsSeqIoOpts);
				auto contigsSeqs = contigsReader.readAllReads<seqInfo>();
				std::vector<std::shared_ptr<seqWithKmerInfo>> contigsKmerReads;
				contigsKmerReads.reserve(contigsSeqs.size());
				for (const auto & seq : contigsSeqs) {
					contigsKmerReads.emplace_back(std::make_shared<seqWithKmerInfo>(seq));
				}
				allSetKmers(contigsKmerReads, reOrientingKmerLength, true);

				SeqInput refReader(SeqIOOptions::genFastaIn(refFnp));
				auto refSeqs = refReader.readAllReads<seqInfo>();
				std::vector<std::shared_ptr<seqWithKmerInfo>> refKmerReads;
				std::vector<seqInfo> revComp_refSeqs;
				std::vector<kmerInfo> revComp_refSeqsKInfos;
				refKmerReads.reserve(refSeqs.size());
				for (const auto & seq : refSeqs) {
					refKmerReads.emplace_back(std::make_shared<seqWithKmerInfo>(seq));
				}
				allSetKmers(refKmerReads, reOrientingKmerLength, true);
				for(const auto & rSeq : refSeqs){
					revComp_refSeqs.emplace_back(rSeq);
					revComp_refSeqs.back().reverseComplementRead(false, true);
					revComp_refSeqsKInfos.emplace_back(revComp_refSeqs.back().seq_, 7, false);
				}
				for(const auto & seqKmer : contigsKmerReads) {
					uint32_t forwardWinners = 0;
					uint32_t revWinners = 0;
					for (const auto & refSeq : refKmerReads) {
						auto forDist = refSeq->compareKmers(*seqKmer);
						auto revDist = refSeq->compareKmersRevComp(*seqKmer);
						if (forDist.first < revDist.first) {
							++revWinners;
						} else {
							++forwardWinners;
						}
					}
					if (revWinners > forwardWinners) {
						seqKmer->seqBase_.reverseComplementRead(true, true);
					}
				}

				//sort by sequence length;
				njh::sort(contigsKmerReads, [](const std::shared_ptr<seqWithKmerInfo> & seq1, const std::shared_ptr<seqWithKmerInfo> & seq2){
					return len(seq1->seqBase_) > len(seq2->seqBase_);
				});

				OutOptions contigInfoOpts(njh::files::make_path(spadesFullOutputDir, "contigs_outputInfo.tab.txt"));
				OutputStream contigInfoOut(contigInfoOpts);
				contigInfoOut << "name\tlength\tcoverage" << std::endl;

				for(const auto & contigsKmerRead : contigsKmerReads){
					auto assembleInfo = DefaultAssembleNameInfo(contigsKmerRead->seqBase_.name_);
					contigInfoOut << contigsKmerRead->seqBase_.name_
							<< "\t" << len(contigsKmerRead->seqBase_)
							<< "\t" << assembleInfo.coverage_ << std::endl;
				}
				auto reOrientedContigsFnp = njh::files::make_path(spadesFullOutputDir, "reOriented_contigs.fasta");

				SeqOutput::write(contigsKmerReads, SeqIOOptions::genFastaOut(reOrientedContigsFnp));

				uint64_t maxLen = 0;
				readVec::getMaxLength(refSeqs, maxLen);
				readVec::getMaxLength(contigsKmerReads, maxLen);
				aligner alignerObj(maxLen, gapScoringParameters(5,1,0,0,0,0), substituteMatrix(2,-2), false);
				//alignerObj.processAlnInfoInputNoCheck(njh::files::make_path(resultsDirectory, "trimAlnCache").string(), setUp.pars_.verbose_);
				std::vector<kmerInfo> refSeqsKmerInfos;
				refSeqsKmerInfos.reserve(refSeqs.size());
				for(const auto & input : refSeqs){
					refSeqsKmerInfos.emplace_back(input.seq_, 7, false);
				}

				//readVecTrimmer::trimSeqToRefByGlobalAln(contigsKmerReads, refSeqs, refSeqsKmerInfos, alignerObj);
				std::vector<std::shared_ptr<seqWithKmerInfo>> finalSeqs;
				std::unordered_map<std::string, uint32_t> finalSeqCounts;
				std::unordered_map<std::string, uint32_t> finalSeqCountsWritten;
				for(const auto & seq : contigsKmerReads){
					auto trimmed = readVecTrimmer::trimSeqToRefByGlobalAlnBestNoOverlapIncludeRevComp(seq, refSeqs, revComp_refSeqs, refSeqsKmerInfos, revComp_refSeqsKInfos, alignerObj, false);
					for(auto & trimmedSeq : trimmed){
						bool found = false;
						for(const auto & finalSeq : finalSeqs){
							if(finalSeq->seqBase_.seq_ == trimmedSeq.seq_){
								found = true;
								break;
							}
						}
						if(!found){
							finalSeqs.emplace_back(std::make_shared<seqWithKmerInfo>(trimmedSeq, 7, false));
							++finalSeqCounts[trimmedSeq.name_];
						}
					}
				}
				double totalCoverage = 0;
				for(auto & seq : finalSeqs){
					auto assembleInfo = DefaultAssembleNameInfo(seq->seqBase_.name_);
					totalCoverage += assembleInfo.coverage_;
				}

				for(auto & seq : finalSeqs){
					auto assembleInfo = DefaultAssembleNameInfo(seq->seqBase_.name_);
					MetaDataInName seqMeta;
					seqMeta.addMeta("trimmedLength", len(seq->seqBase_));
					seqMeta.addMeta("estimatedPerBaseCoverage", assembleInfo.coverage_);
					seqMeta.addMeta("trimStatus", seq->seqBase_.on_);
					seqMeta.addMeta("regionUID", regionUid);
					seqMeta.addMeta("sample", sample);
					if(finalSeqCounts[seq->seqBase_.name_] > 1){
						seqMeta.addMeta("seqTrimmedCount", finalSeqCountsWritten[seq->seqBase_.name_]);
						++finalSeqCountsWritten[seq->seqBase_.name_];
					}
					seqMeta.resetMetaInName(seq->seqBase_.name_);
					seq->seqBase_.cnt_ = (assembleInfo.coverage_/totalCoverage) * (readCounts.pairedReads_ + readCounts.unpaiedReads_ + readCounts.orphans_);
					seq->seqBase_.name_ += njh::pasteAsStr("_t", seq->seqBase_.cnt_);
				}

				OutOptions trimmedContigInfoOpts(njh::files::make_path(spadesFullOutputDir, "trimmed_reOriented_contigs_outputInfo.tab.txt"));
				OutputStream trimmedContigInfoOut(trimmedContigInfoOpts);
				trimmedContigInfoOut << "name\tlength\tcoverage" << std::endl;
				auto trimmedReOrientedContigsFnp = njh::files::make_path(spadesFullOutputDir, "trimmed_reOriented_contigs.fasta");
				SeqOutput outputWriter(SeqIOOptions::genFastaOut(trimmedReOrientedContigsFnp));
				auto trimmedReOrientedContigsFnp_belowCutOff = njh::files::make_path(spadesFullOutputDir, "trimmed_reOriented_contigs_belowCutOff.fasta");
				SeqOutput belowCutOffOutputWriter(SeqIOOptions::genFastaOut(trimmedReOrientedContigsFnp_belowCutOff));

				uint32_t belowCutOff = 0;
				uint32_t aboveCutOff = 0;
				bool allPassTrim = true;
				for (const auto & contigsKmerRead : finalSeqs) {
					if (len(contigsKmerRead->seqBase_) < minFinalLength) {
						++belowCutOff;
						belowCutOffOutputWriter.openWrite(contigsKmerRead);
						contigsKmerRead->seqBase_.on_ = false;
					} else {
						MetaDataInName seqMeta(contigsKmerRead->seqBase_.name_);
						trimmedContigInfoOut << contigsKmerRead->seqBase_.name_
								<< "\t" << len(contigsKmerRead->seqBase_)
								<< "\t" << seqMeta.getMeta("estimatedPerBaseCoverage")
								<< std::endl;
						if(!contigsKmerRead->seqBase_.on_){
							allPassTrim = false;
						}else{
							++aboveCutOff;
						}
						outputWriter.openWrite(contigsKmerRead);
					}
				}
				if(allPassTrim && !finalSeqs.empty()){
					std::lock_guard<std::mutex> lock(allFinalWriterMut);
					for (const auto & contigsKmerRead : finalSeqs) {
						if (len(contigsKmerRead->seqBase_) >= minFinalLength) {
							allFinalWriter.write(contigsKmerRead);
						}
					}
					for(auto & reg : regInfo){
						reg->infoCalled_ = true;
						reg->uniqHaps_ = aboveCutOff;
					}
				}else{
					std::lock_guard<std::mutex> lock(allPartialWriterMut);
					for (const auto & contigsKmerRead : finalSeqs) {
						if (len(contigsKmerRead->seqBase_) >= minFinalLength) {
							allPartialWriter.write(contigsKmerRead);
						}
					}
					for(auto & reg : regInfo){
						reg->infoCalled_ = false;
						reg->uniqHaps_ = 0;
					}
				}
			} catch (std::exception & e) {
				std::lock_guard<std::mutex> lock(exceptionsMut);
				exceptions[regionUid] = e.what();
				for(auto & reg : regInfo){
					reg->infoCalled_ = false;
					reg->uniqHaps_ = 0;
				}
			}
		}
	};


	njh::concurrent::runVoidFunctionThreaded(runSpadesOnRegion, numThreads);
	allFinalWriter.closeOut();
	allPartialWriter.closeOut();
	//sample,readTotal,readTotalUsed, success, name
	//
	OutputStream basicInfo(njh::files::make_path(finalDirectory, "basicInfoPerRegion.tab.txt"));

	basicInfo << "#chrom\tstart\tend\tname\tlength\tstrand\tsuccess\tuniqHaps\treadTotal\treadTotalUsed\ttotalPairedReads";
	basicInfo << "\tsample";
	uint32_t maxExtraFields = 0;
	for(const auto & p : inputRegions){

		auto bedOut = p.genBedRecordCore();
		if(bedOut.extraFields_.size() > maxExtraFields){
			maxExtraFields = bedOut.extraFields_.size();
		}
	}
	for(uint32_t t = 0; t < maxExtraFields; ++t){
		basicInfo << "\textraField"<<t;
	}
	basicInfo << "\n";

	std::map<uint32_t, uint32_t> coiCounts;

	for (const auto & reg : inputRegions) {
		const auto & regInfos = njh::mapAt(regInfosByUID, reg.uid_);
		++coiCounts[regInfos.front()->uniqHaps_];
		for(auto & regInfo : regInfos){
			auto bedOut = regInfo->region_.genBedRecordCore();
			basicInfo << bedOut.toDelimStr();
			basicInfo << "\t" << njh::boolToStr(regInfo->infoCalled_)
								<< "\t" << regInfo->uniqHaps_
								<< "\t" << regInfo->totalReads_
								<< "\t" << regInfo->totalFinalReads_
								<< "\t" << regInfo->totalPairedReads_
								<< "\t" << sample;
			for(const auto & extra : bedOut.extraFields_){
				basicInfo << "\t" << extra;
			}
			basicInfo << std::endl;
		}
	}

	OutputStream coiOut(njh::files::make_path(finalDirectory, "coiCounts.tab.txt"));
	coiOut << "coi\tcount" << std::endl;
	for(const auto & count : coiCounts){
		coiOut << count.first << "\t" << count.second << std::endl;
	}

	OutputStream exceptionsOut(njh::files::make_path(finalDirectory, "exceptionsMessages.tab.txt"));
	exceptionsOut << "regionUID\tmessage" << std::endl;
	for(const auto & exp : exceptions){
		exceptionsOut << exp.first << "\t" << exp.second << std::endl;
	}



	return 0;
}

//


//


int programWrappersAssembleOnPathWeaverRunner::runRayOnPathWeaverRegions(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bedFile = "";
	bfs::path pwOutputDir = "";
	std::string sample;

	//uint32_t RayNumThreads = 1;
	std::string extraRayOptions;
	uint32_t reOrientingKmerLength = 9;
	uint32_t minFinalLength = 40;
	uint32_t RayKmerLength = std::numeric_limits<uint32_t>::max();
	bfs::path RayOutDir = "RayOut";
	uint32_t numThreads = 1;
	seqSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processVerbose();
	setUp.setOption(bedFile, "--bed", "The Regions to analyze", true);
	setUp.setOption(pwOutputDir, "--pwOutputDir", "The PathWeaver directory", true);
	setUp.setOption(sample, "--sample", "sample name", true);

	setUp.setOption(RayKmerLength, "--RayKmerLength", "Ray Kmer Length");

	setUp.setOption(numThreads, "--numThreads", "num Threads");


	//setUp.setOption(RayNumThreads, "--RayNumThreads", "Ray Num Threads");
	setUp.setOption(extraRayOptions, "--extraRayOptions", "extra Ray Options");

	setUp.setOption(minFinalLength, "--minFinalLength", "min Final Length");
	setUp.setOption(reOrientingKmerLength, "--reOrientingKmerLength", "re-orienting K-mer Length");

	setUp.setOption(RayOutDir,     "--RayOutDir",     "Ray Out Directory name, will be relative to final pass directory");


	setUp.processDirectoryOutputName(njh::pasteAsStr(bfs::basename(pwOutputDir), "_Ray_TODAY"), true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	njh::sys::requireExternalProgramThrow("Ray");

	auto inputRegions = gatherRegions(bedFile.string(), "", setUp.pars_.verbose_);
	sortGRegionsByStart(inputRegions);

	std::set<std::string> regionNames;
	for(const auto & reg : inputRegions){
		regionNames.emplace(reg.uid_);
	}
	//njh::sort(regionNames);
	njh::concurrent::LockableQueue<std::string> regionsQueue(regionNames);

	bfs::path finalDirectory = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar("final"));
	bfs::path partialDirectory = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar("partial"));
	auto allFinalSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(finalDirectory, "allFinal.fasta"));
	auto allPartialSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(partialDirectory, "allPartial.fasta"));
	SeqOutput allFinalWriter(allFinalSeqOpts);
	SeqOutput allPartialWriter(allPartialSeqOpts);
	allFinalWriter.openOut();
	allPartialWriter.openOut();
	std::mutex allFinalWriterMut;
	std::mutex allPartialWriterMut;


	std::unordered_map<std::string,
			std::vector<std::shared_ptr<BamRegionInvestigator::RegionInfo>> >regInfosByUID;
	for (const auto & reg : inputRegions) {
		regInfosByUID[reg.uid_].emplace_back(std::make_shared<BamRegionInvestigator::RegionInfo>(reg));
	}

	std::unordered_map<std::string, std::string> exceptions;
	std::mutex exceptionsMut;

	std::function<void()> runRayOnRegion = [&](){
		std::string regionUid;
		while(regionsQueue.getVal(regionUid)){
			const auto & regInfo = njh::mapAt(regInfosByUID, regionUid);
			auto regionOutputDir = njh::files::make_path(setUp.pars_.directoryName_, regionUid, sample);
			njh::files::makeDirP(njh::files::MkdirPar{regionOutputDir});

			//

			bfs::path refFnp = njh::files::make_path(pwOutputDir, regionUid, "allRefs.fasta");


			//first extract the reads
			bfs::path extractBam = njh::files::make_path(pwOutputDir, regionUid, sample + "_extraction", "extracted.bam");
			OutOptions outOpts(njh::files::make_path(regionOutputDir, "extracted"));
			auto readCounts = rawWriteExtractReadsFromBamOnlyMapped(extractBam, outOpts);
			bfs::path pairedR1 = njh::files::make_path(regionOutputDir, "extracted_R1.fastq");
			bfs::path pairedR2 = njh::files::make_path(regionOutputDir, "extracted_R2.fastq");
			bfs::path singles =  njh::files::make_path(regionOutputDir, "extracted.fastq");
			for(auto & reg : regInfo){
				reg->totalPairedReads_ = readCounts.pairedReads_;
				reg->totalReads_ = readCounts.pairedReads_ + readCounts.unpaiedReads_ + readCounts.orphans_;
				reg->totalFinalReads_ = readCounts.pairedReads_ + readCounts.unpaiedReads_ + readCounts.orphans_;
			}
			try {
				bfs::path optimJsonFnp = njh::files::make_path(pwOutputDir, regionUid, sample, "optimizationInfoBest.json");
				Json::Value optimJson = njh::json::parseFile(optimJsonFnp.string());
				if(std::numeric_limits<uint32_t>::max() == RayKmerLength){
					RayKmerLength = optimJson["runParams_"]["klen_"].asUInt64() ;
				}
//				std::cout << njh::conToStr(optimJson.getMemberNames(), "\n") << std::endl;
//
//				std::cout << "optimJson[\"runParams_\"][\"klen_\"].asUInt64(): " << optimJson["runParams_"]["klen_"].asUInt64() << std::endl;
//
//				exit(1);
				if(!exists(pairedR1) && !exists(singles)){
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", couldn't find " << pairedR1 << " or " << singles << ", need to have at least one of them" << "\n";
					throw std::runtime_error{ss.str()};
				}
				std::stringstream RayCmdStream;
				RayCmdStream << "cd " << regionOutputDir << " && Ray ";
				if(!exists(pairedR1)){
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error " << "Ray requires paired reads"<< "\n";
					throw std::runtime_error{ss.str()};
				}
//
//
				RayCmdStream    << " -k " << RayKmerLength
				                << " -p " << pairedR1.filename() << " " <<  pairedR2.filename()
												<< " " << extraRayOptions
												<< " -o " << RayOutDir
												<< " > RayRunLog_" << njh::getCurrentDate() << ".txt 2>&1";
				auto RayFullOutputDir = njh::files::make_path(regionOutputDir, RayOutDir);

				auto RayRunOutput = njh::sys::run({RayCmdStream.str()});

				BioCmdsUtils::checkRunOutThrow(RayRunOutput, __PRETTY_FUNCTION__);

				OutOptions RayRunOutputLogOpts(njh::files::make_path(RayFullOutputDir, "RayRunOutput.json"));
				OutputStream RayRunOutputLogOut(RayRunOutputLogOpts);
				RayRunOutputLogOut << njh::json::toJson(RayRunOutput) << std::endl;

				auto contigsFnp = njh::files::make_path(RayFullOutputDir, "Contigs.fasta");


				std::unordered_map<std::string, double> kmerCoverage;
				auto coverageInfo = njh::files::make_path(RayFullOutputDir, "BiologicalAbundances/_DeNovoAssembly/Contigs.tsv");
				table covTab(coverageInfo, "\t", true);
				for(const auto & row : covTab){
					kmerCoverage[row[covTab.getColPos("#Contig name")]] = njh::StrToNumConverter::stoToNum<double>(row[covTab.getColPos("Mode k-mer coverage depth")]);
				}
				auto contigsSeqIoOpts = SeqIOOptions::genFastaIn(contigsFnp);
				contigsSeqIoOpts.includeWhiteSpaceInName_ = false;
				contigsSeqIoOpts.lowerCaseBases_ = "upper";
				SeqInput contigsReader(contigsSeqIoOpts);
				auto contigsSeqs = contigsReader.readAllReads<seqInfo>();
				std::vector<std::shared_ptr<seqWithKmerInfo>> contigsKmerReads;
				contigsKmerReads.reserve(contigsSeqs.size());
				for (const auto & seq : contigsSeqs) {
					contigsKmerReads.emplace_back(std::make_shared<seqWithKmerInfo>(seq));
				}
				allSetKmers(contigsKmerReads, reOrientingKmerLength, true);

				SeqInput refReader(SeqIOOptions::genFastaIn(refFnp));
				auto refSeqs = refReader.readAllReads<seqInfo>();
				std::vector<std::shared_ptr<seqWithKmerInfo>> refKmerReads;
				std::vector<seqInfo> revComp_refSeqs;
				std::vector<kmerInfo> revComp_refSeqsKInfos;
				refKmerReads.reserve(refSeqs.size());
				for (const auto & seq : refSeqs) {
					refKmerReads.emplace_back(std::make_shared<seqWithKmerInfo>(seq));
				}
				allSetKmers(refKmerReads, reOrientingKmerLength, true);
				for(const auto & rSeq : refSeqs){
					revComp_refSeqs.emplace_back(rSeq);
					revComp_refSeqs.back().reverseComplementRead(false, true);
					revComp_refSeqsKInfos.emplace_back(revComp_refSeqs.back().seq_, 7, false);
				}
				for(const auto & seqKmer : contigsKmerReads) {
					uint32_t forwardWinners = 0;
					uint32_t revWinners = 0;
					for (const auto & refSeq : refKmerReads) {
						auto forDist = refSeq->compareKmers(*seqKmer);
						auto revDist = refSeq->compareKmersRevComp(*seqKmer);
						if (forDist.first < revDist.first) {
							++revWinners;
						} else {
							++forwardWinners;
						}
					}
					if (revWinners > forwardWinners) {
						auto oldName = seqKmer->seqBase_.name_;
						seqKmer->seqBase_.reverseComplementRead(true, true);
						kmerCoverage[seqKmer->seqBase_.name_] = kmerCoverage[oldName];
					}
				}

				//sort by sequence length;
				njh::sort(contigsKmerReads, [](const std::shared_ptr<seqWithKmerInfo> & seq1, const std::shared_ptr<seqWithKmerInfo> & seq2){
					return len(seq1->seqBase_) > len(seq2->seqBase_);
				});

				OutOptions contigInfoOpts(njh::files::make_path(RayFullOutputDir, "contigs_outputInfo.tab.txt"));
				OutputStream contigInfoOut(contigInfoOpts);
				contigInfoOut << "name\tlength\tcoverage" << std::endl;

				for(const auto & contigsKmerRead : contigsKmerReads){
					//auto assembleInfo = DefaultAssembleNameInfo(contigsKmerRead->seqBase_.name_, true);
					contigInfoOut << contigsKmerRead->seqBase_.name_
							<< "\t" << len(contigsKmerRead->seqBase_)
							<< "\t" << kmerCoverage[contigsKmerRead->seqBase_.name_] << std::endl;
				}
				auto reOrientedContigsFnp = njh::files::make_path(RayFullOutputDir, "reOriented_contigs.fasta");

				SeqOutput::write(contigsKmerReads, SeqIOOptions::genFastaOut(reOrientedContigsFnp));

				uint64_t maxLen = 0;
				readVec::getMaxLength(refSeqs, maxLen);
				readVec::getMaxLength(contigsKmerReads, maxLen);
				aligner alignerObj(maxLen, gapScoringParameters(5,1,0,0,0,0), substituteMatrix(2,-2), false);
				//alignerObj.processAlnInfoInputNoCheck(njh::files::make_path(resultsDirectory, "trimAlnCache").string(), setUp.pars_.verbose_);
				std::vector<kmerInfo> refSeqsKmerInfos;
				refSeqsKmerInfos.reserve(refSeqs.size());
				for(const auto & input : refSeqs){
					refSeqsKmerInfos.emplace_back(input.seq_, 7, false);
				}

				//readVecTrimmer::trimSeqToRefByGlobalAln(contigsKmerReads, refSeqs, refSeqsKmerInfos, alignerObj);
				std::vector<std::shared_ptr<seqWithKmerInfo>> finalSeqs;
				std::unordered_map<std::string, uint32_t> finalSeqCounts;
				std::unordered_map<std::string, uint32_t> finalSeqCountsWritten;
				for(const auto & seq : contigsKmerReads){
					auto trimmed = readVecTrimmer::trimSeqToRefByGlobalAlnBestNoOverlapIncludeRevComp(seq, refSeqs, revComp_refSeqs, refSeqsKmerInfos, revComp_refSeqsKInfos, alignerObj, false);
					for(auto & trimmedSeq : trimmed){
						bool found = false;
						for(const auto & finalSeq : finalSeqs){
							if(finalSeq->seqBase_.seq_ == trimmedSeq.seq_){
								found = true;
								break;
							}
						}
						if(!found){
							finalSeqs.emplace_back(std::make_shared<seqWithKmerInfo>(trimmedSeq, 7, false));
							++finalSeqCounts[trimmedSeq.name_];
						}
					}
				}

				double totalCoverage = 0;
				for(auto & seq : finalSeqs){
					//auto assembleInfo = DefaultAssembleNameInfo(seq->seqBase_.name_, true);
					totalCoverage += kmerCoverage[seq->seqBase_.name_];
				}

				for(auto & seq : finalSeqs){
					//auto assembleInfo = DefaultAssembleNameInfo(seq->seqBase_.name_, true);
					MetaDataInName seqMeta;
					seqMeta.addMeta("trimmedLength", len(seq->seqBase_));
					seqMeta.addMeta("estimatedPerBaseCoverage", kmerCoverage[seq->seqBase_.name_]);
					seqMeta.addMeta("trimStatus", seq->seqBase_.on_);
					seqMeta.addMeta("regionUID", regionUid);
					seqMeta.addMeta("sample", sample);
					if(finalSeqCounts[seq->seqBase_.name_] > 1){
						seqMeta.addMeta("seqTrimmedCount", finalSeqCountsWritten[seq->seqBase_.name_]);
						++finalSeqCountsWritten[seq->seqBase_.name_];
					}
					//seq->seqBase_.cnt_ = (defaultCoverage/totalCoverage) * (readCounts.pairedReads_ + readCounts.unpaiedReads_ + readCounts.orphans_);
//					std::cout << "kmerCoverage[seq->seqBase_.name_]: " << kmerCoverage[seq->seqBase_.name_] << std::endl;
//					std::cout << "totalCoverage: " << totalCoverage << std::endl;
//					std::cout << "reads: " << readCounts.pairedReads_ + readCounts.unpaiedReads_ + readCounts.orphans_ << std::endl;
					seq->seqBase_.cnt_ = (kmerCoverage[seq->seqBase_.name_]/totalCoverage) * (readCounts.pairedReads_ + readCounts.unpaiedReads_ + readCounts.orphans_);
					seqMeta.resetMetaInName(seq->seqBase_.name_);
					seq->seqBase_.name_ += njh::pasteAsStr("_t", seq->seqBase_.cnt_);
				}

				OutOptions trimmedContigInfoOpts(njh::files::make_path(RayFullOutputDir, "trimmed_reOriented_contigs_outputInfo.tab.txt"));
				OutputStream trimmedContigInfoOut(trimmedContigInfoOpts);
				trimmedContigInfoOut << "name\tlength\tcoverage" << std::endl;
				auto trimmedReOrientedContigsFnp = njh::files::make_path(RayFullOutputDir, "trimmed_reOriented_contigs.fasta");
				SeqOutput outputWriter(SeqIOOptions::genFastaOut(trimmedReOrientedContigsFnp));
				auto trimmedReOrientedContigsFnp_belowCutOff = njh::files::make_path(RayFullOutputDir, "trimmed_reOriented_contigs_belowCutOff.fasta");
				SeqOutput belowCutOffOutputWriter(SeqIOOptions::genFastaOut(trimmedReOrientedContigsFnp_belowCutOff));

				uint32_t belowCutOff = 0;
				uint32_t aboveCutOff = 0;
				bool allPassTrim = true;
				for (const auto & contigsKmerRead : finalSeqs) {
					if (len(contigsKmerRead->seqBase_) < minFinalLength) {
						++belowCutOff;
						belowCutOffOutputWriter.openWrite(contigsKmerRead);
						contigsKmerRead->seqBase_.on_ = false;
					} else {
						MetaDataInName seqMeta(contigsKmerRead->seqBase_.name_);
						trimmedContigInfoOut << contigsKmerRead->seqBase_.name_
								<< "\t" << len(contigsKmerRead->seqBase_)
								<< "\t" << seqMeta.getMeta("estimatedPerBaseCoverage")
								<< std::endl;
						if(!contigsKmerRead->seqBase_.on_){
							allPassTrim = false;
						}else{
							++aboveCutOff;
						}
						outputWriter.openWrite(contigsKmerRead);
					}
				}
				if(allPassTrim && !finalSeqs.empty()){
					std::lock_guard<std::mutex> lock(allFinalWriterMut);
					for (const auto & contigsKmerRead : finalSeqs) {
						if (len(contigsKmerRead->seqBase_) >= minFinalLength) {
							allFinalWriter.write(contigsKmerRead);
						}
					}
					for(auto & reg : regInfo){
						reg->infoCalled_ = true;
						reg->uniqHaps_ = aboveCutOff;
					}
				}else{
					std::lock_guard<std::mutex> lock(allPartialWriterMut);
					for (const auto & contigsKmerRead : finalSeqs) {
						if (len(contigsKmerRead->seqBase_) >= minFinalLength) {
							allPartialWriter.write(contigsKmerRead);
						}
					}
					for(auto & reg : regInfo){
						reg->infoCalled_ = false;
						reg->uniqHaps_ = 0;
					}
				}
			} catch (std::exception & e) {
				std::lock_guard<std::mutex> lock(exceptionsMut);
				exceptions[regionUid] = e.what();
				for(auto & reg : regInfo){
					reg->infoCalled_ = false;
					reg->uniqHaps_ = 0;
				}
			}
		}
	};


	njh::concurrent::runVoidFunctionThreaded(runRayOnRegion, numThreads);
	allFinalWriter.closeOut();
	allPartialWriter.closeOut();
	//sample,readTotal,readTotalUsed, success, name
	//
	OutputStream basicInfo(njh::files::make_path(finalDirectory, "basicInfoPerRegion.tab.txt"));

	basicInfo << "#chrom\tstart\tend\tname\tlength\tstrand\tsuccess\tuniqHaps\treadTotal\treadTotalUsed\ttotalPairedReads";
	basicInfo << "\tsample";
	uint32_t maxExtraFields = 0;
	for(const auto & p : inputRegions){

		auto bedOut = p.genBedRecordCore();
		if(bedOut.extraFields_.size() > maxExtraFields){
			maxExtraFields = bedOut.extraFields_.size();
		}
	}
	for(uint32_t t = 0; t < maxExtraFields; ++t){
		basicInfo << "\textraField"<<t;
	}
	basicInfo << "\n";

	std::map<uint32_t, uint32_t> coiCounts;

	for (const auto & reg : inputRegions) {
		const auto & regInfos = njh::mapAt(regInfosByUID, reg.uid_);
		++coiCounts[regInfos.front()->uniqHaps_];
		for(auto & regInfo : regInfos){
			auto bedOut = regInfo->region_.genBedRecordCore();
			basicInfo << bedOut.toDelimStr();
			basicInfo << "\t" << njh::boolToStr(regInfo->infoCalled_)
								<< "\t" << regInfo->uniqHaps_
								<< "\t" << regInfo->totalReads_
								<< "\t" << regInfo->totalFinalReads_
								<< "\t" << regInfo->totalPairedReads_
								<< "\t" << sample;
			for(const auto & extra : bedOut.extraFields_){
				basicInfo << "\t" << extra;
			}
			basicInfo << std::endl;
		}
	}

	OutputStream coiOut(njh::files::make_path(finalDirectory, "coiCounts.tab.txt"));
	coiOut << "coi\tcount" << std::endl;
	for(const auto & count : coiCounts){
		coiOut << count.first << "\t" << count.second << std::endl;
	}

	OutputStream exceptionsOut(njh::files::make_path(finalDirectory, "exceptionsMessages.tab.txt"));
	exceptionsOut << "regionUID\tmessage" << std::endl;
	for(const auto & exp : exceptions){
		exceptionsOut << exp.first << "\t" << exp.second << std::endl;
	}



	return 0;
}


int programWrappersAssembleOnPathWeaverRunner::runIDBAUDOnPathWeaverRegions(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bedFile = "";
	bfs::path pwOutputDir = "";
	std::string sample;

	uint32_t IDBAUDNumThreads = 1;
	std::string extraIDBAUDOptions;
	uint32_t reOrientingKmerLength = 9;
	uint32_t minFinalLength = 40;
	uint32_t IDBAUDKmerLength = std::numeric_limits<uint32_t>::max();
	bfs::path IDBAUDOutDir = "IDBAUDOut";
	uint32_t numThreads = 1;
	seqSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processVerbose();
	setUp.setOption(bedFile, "--bed", "The Regions to analyze", true);
	setUp.setOption(pwOutputDir, "--pwOutputDir", "The PathWeaver directory", true);
	setUp.setOption(sample, "--sample", "sample name", true);

	setUp.setOption(IDBAUDKmerLength, "--IDBAUDKmerLength", "IDBAUD Kmer Length");

	setUp.setOption(numThreads, "--numThreads", "num Threads");


	setUp.setOption(IDBAUDNumThreads, "--IDBAUDNumThreads", "IDBAUD Num Threads");
	setUp.setOption(extraIDBAUDOptions, "--extraIDBAUDOptions", "extra IDBAUD Options");

	setUp.setOption(minFinalLength, "--minFinalLength", "min Final Length");
	setUp.setOption(reOrientingKmerLength, "--reOrientingKmerLength", "re-orienting K-mer Length");

	setUp.setOption(IDBAUDOutDir,     "--IDBAUDOutDir",     "IDBAUD Out Directory name, will be relative to final pass directory");


	setUp.processDirectoryOutputName(njh::pasteAsStr(bfs::basename(pwOutputDir), "_IDBAUD_TODAY"), true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	njh::sys::requireExternalProgramThrow("fq2fa");
	njh::sys::requireExternalProgramThrow("idba_ud");

	auto inputRegions = gatherRegions(bedFile.string(), "", setUp.pars_.verbose_);
	sortGRegionsByStart(inputRegions);

	std::set<std::string> regionNames;
	for(const auto & reg : inputRegions){
		regionNames.emplace(reg.uid_);
	}
	//njh::sort(regionNames);
	njh::concurrent::LockableQueue<std::string> regionsQueue(regionNames);

	bfs::path finalDirectory = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar("final"));
	bfs::path partialDirectory = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar("partial"));
	auto allFinalSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(finalDirectory, "allFinal.fasta"));
	auto allPartialSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(partialDirectory, "allPartial.fasta"));
	SeqOutput allFinalWriter(allFinalSeqOpts);
	SeqOutput allPartialWriter(allPartialSeqOpts);
	allFinalWriter.openOut();
	allPartialWriter.openOut();
	std::mutex allFinalWriterMut;
	std::mutex allPartialWriterMut;


	std::unordered_map<std::string,
			std::vector<std::shared_ptr<BamRegionInvestigator::RegionInfo>> >regInfosByUID;
	for (const auto & reg : inputRegions) {
		regInfosByUID[reg.uid_].emplace_back(std::make_shared<BamRegionInvestigator::RegionInfo>(reg));
	}

	std::unordered_map<std::string, std::string> exceptions;
	std::mutex exceptionsMut;

	std::function<void()> runIDBAUDOnRegion = [&](){
		std::string regionUid;
		while(regionsQueue.getVal(regionUid)){
			const auto & regInfo = njh::mapAt(regInfosByUID, regionUid);
			auto regionOutputDir = njh::files::make_path(setUp.pars_.directoryName_, regionUid, sample);
			njh::files::makeDirP(njh::files::MkdirPar{regionOutputDir});

			//

			bfs::path refFnp = njh::files::make_path(pwOutputDir, regionUid, "allRefs.fasta");


			//first extract the reads
			bfs::path extractBam = njh::files::make_path(pwOutputDir, regionUid, sample + "_extraction", "extracted.bam");
			OutOptions outOpts(njh::files::make_path(regionOutputDir, "extracted"));
			auto readCounts = rawWriteExtractReadsFromBamOnlyMapped(extractBam, outOpts);
			bfs::path pairedR1 = njh::files::make_path(regionOutputDir, "extracted_R1.fastq");
			bfs::path pairedR2 = njh::files::make_path(regionOutputDir, "extracted_R2.fastq");
			bfs::path singles =  njh::files::make_path(regionOutputDir, "extracted.fastq");
			for(auto & reg : regInfo){
				reg->totalPairedReads_ = readCounts.pairedReads_;
				reg->totalReads_ = readCounts.pairedReads_ + readCounts.unpaiedReads_ + readCounts.orphans_;
				reg->totalFinalReads_ = readCounts.pairedReads_ + readCounts.unpaiedReads_ + readCounts.orphans_;
			}
			try {
				bfs::path optimJsonFnp = njh::files::make_path(pwOutputDir, regionUid, sample, "optimizationInfoBest.json");
				Json::Value optimJson = njh::json::parseFile(optimJsonFnp.string());
				if(std::numeric_limits<uint32_t>::max() == IDBAUDKmerLength){
					IDBAUDKmerLength = optimJson["runParams_"]["klen_"].asUInt64() ;
				}
//				std::cout << njh::conToStr(optimJson.getMemberNames(), "\n") << std::endl;
//
//				std::cout << "optimJson[\"runParams_\"][\"klen_\"].asUInt64(): " << optimJson["runParams_"]["klen_"].asUInt64() << std::endl;
//
//				exit(1);
				if(!exists(pairedR1) && !exists(singles)){
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", couldn't find " << pairedR1 << " or " << singles << ", need to have at least one of them" << "\n";
					throw std::runtime_error{ss.str()};
				}
				std::stringstream IDBAUDCmdStream;

				IDBAUDCmdStream << "cd " << regionOutputDir ;

				IDBAUDCmdStream << " && fq2fa --merge extracted_R1.fastq extracted_R2.fastq extracted.fasta";

				IDBAUDCmdStream << " && idba_ud ";
				if(!exists(pairedR1)){
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error " << "IDBAUD requires paired reads"<< "\n";
					throw std::runtime_error{ss.str()};
				}
//
//
				//
				//idba_ud  -o ibda_testing

				IDBAUDCmdStream
				                << " -r extracted.fasta "
												<< " --num_threads " << IDBAUDNumThreads
												<< " " << extraIDBAUDOptions
												<< " -o " << IDBAUDOutDir
												<< " > IDBAUDRunLog_" << njh::getCurrentDate() << ".txt 2>&1";
				auto IDBAUDFullOutputDir = njh::files::make_path(regionOutputDir, IDBAUDOutDir);

				auto IDBAUDRunOutput = njh::sys::run({IDBAUDCmdStream.str()});

				BioCmdsUtils::checkRunOutThrow(IDBAUDRunOutput, __PRETTY_FUNCTION__);

				OutOptions IDBAUDRunOutputLogOpts(njh::files::make_path(IDBAUDFullOutputDir, "IDBAUDRunOutput.json"));
				OutputStream IDBAUDRunOutputLogOut(IDBAUDRunOutputLogOpts);
				IDBAUDRunOutputLogOut << njh::json::toJson(IDBAUDRunOutput) << std::endl;

				auto contigsFnp = njh::files::make_path(IDBAUDFullOutputDir, "contig.fa");


				auto contigsSeqIoOpts = SeqIOOptions::genFastaIn(contigsFnp);
				contigsSeqIoOpts.includeWhiteSpaceInName_ = false;
				contigsSeqIoOpts.lowerCaseBases_ = "upper";
				SeqInput contigsReader(contigsSeqIoOpts);
				auto contigsSeqs = contigsReader.readAllReads<seqInfo>();
				std::vector<std::shared_ptr<seqWithKmerInfo>> contigsKmerReads;
				contigsKmerReads.reserve(contigsSeqs.size());
				for (const auto & seq : contigsSeqs) {
					contigsKmerReads.emplace_back(std::make_shared<seqWithKmerInfo>(seq));
				}
				allSetKmers(contigsKmerReads, reOrientingKmerLength, true);

				SeqInput refReader(SeqIOOptions::genFastaIn(refFnp));
				auto refSeqs = refReader.readAllReads<seqInfo>();
				std::vector<std::shared_ptr<seqWithKmerInfo>> refKmerReads;
				std::vector<seqInfo> revComp_refSeqs;
				std::vector<kmerInfo> revComp_refSeqsKInfos;
				refKmerReads.reserve(refSeqs.size());
				for (const auto & seq : refSeqs) {
					refKmerReads.emplace_back(std::make_shared<seqWithKmerInfo>(seq));
				}
				allSetKmers(refKmerReads, reOrientingKmerLength, true);
				for(const auto & rSeq : refSeqs){
					revComp_refSeqs.emplace_back(rSeq);
					revComp_refSeqs.back().reverseComplementRead(false, true);
					revComp_refSeqsKInfos.emplace_back(revComp_refSeqs.back().seq_, 7, false);
				}
				for(const auto & seqKmer : contigsKmerReads) {
					uint32_t forwardWinners = 0;
					uint32_t revWinners = 0;
					for (const auto & refSeq : refKmerReads) {
						auto forDist = refSeq->compareKmers(*seqKmer);
						auto revDist = refSeq->compareKmersRevComp(*seqKmer);
						if (forDist.first < revDist.first) {
							++revWinners;
						} else {
							++forwardWinners;
						}
					}
					if (revWinners > forwardWinners) {
						auto oldName = seqKmer->seqBase_.name_;
						seqKmer->seqBase_.reverseComplementRead(true, true);
//						kmerCoverage[seqKmer->seqBase_.name_] = kmerCoverage[oldName];
					}
				}

				uint32_t defaultCoverage = 10;

				//sort by sequence length;
				njh::sort(contigsKmerReads, [](const std::shared_ptr<seqWithKmerInfo> & seq1, const std::shared_ptr<seqWithKmerInfo> & seq2){
					return len(seq1->seqBase_) > len(seq2->seqBase_);
				});

				OutOptions contigInfoOpts(njh::files::make_path(IDBAUDFullOutputDir, "contigs_outputInfo.tab.txt"));
				OutputStream contigInfoOut(contigInfoOpts);
				contigInfoOut << "name\tlength\tcoverage" << std::endl;

				for(const auto & contigsKmerRead : contigsKmerReads){
					//auto assembleInfo = DefaultAssembleNameInfo(contigsKmerRead->seqBase_.name_, true);
					contigInfoOut << contigsKmerRead->seqBase_.name_
							<< "\t" << len(contigsKmerRead->seqBase_)
							<< "\t" << defaultCoverage << std::endl;
				}
				auto reOrientedContigsFnp = njh::files::make_path(IDBAUDFullOutputDir, "reOriented_contigs.fasta");

				SeqOutput::write(contigsKmerReads, SeqIOOptions::genFastaOut(reOrientedContigsFnp));

				uint64_t maxLen = 0;
				readVec::getMaxLength(refSeqs, maxLen);
				readVec::getMaxLength(contigsKmerReads, maxLen);
				aligner alignerObj(maxLen, gapScoringParameters(5,1,0,0,0,0), substituteMatrix(2,-2), false);
				//alignerObj.processAlnInfoInputNoCheck(njh::files::make_path(resultsDirectory, "trimAlnCache").string(), setUp.pars_.verbose_);
				std::vector<kmerInfo> refSeqsKmerInfos;
				refSeqsKmerInfos.reserve(refSeqs.size());
				for(const auto & input : refSeqs){
					refSeqsKmerInfos.emplace_back(input.seq_, 7, false);
				}
				//readVecTrimmer::trimSeqToRefByGlobalAln(contigsKmerReads, refSeqs, refSeqsKmerInfos, alignerObj);
				std::vector<std::shared_ptr<seqWithKmerInfo>> finalSeqs;
				std::unordered_map<std::string, uint32_t> finalSeqCounts;
				std::unordered_map<std::string, uint32_t> finalSeqCountsWritten;
				for(const auto & seq : contigsKmerReads){
					auto trimmed = readVecTrimmer::trimSeqToRefByGlobalAlnBestNoOverlapIncludeRevComp(seq, refSeqs, revComp_refSeqs, refSeqsKmerInfos, revComp_refSeqsKInfos, alignerObj, false);
					for(auto & trimmedSeq : trimmed){
						bool found = false;
						for(const auto & finalSeq : finalSeqs){
							if(finalSeq->seqBase_.seq_ == trimmedSeq.seq_){
								found = true;
								break;
							}
						}
						if(!found){
							finalSeqs.emplace_back(std::make_shared<seqWithKmerInfo>(trimmedSeq, 7, false));
							++finalSeqCounts[trimmedSeq.name_];
						}
					}
				}
				double totalCoverage = defaultCoverage * finalSeqs.size();
//				for(auto & seq : finalSeqs){
//					//auto assembleInfo = DefaultAssembleNameInfo(seq->seqBase_.name_, true);
//					totalCoverage += kmerCoverage[seq->seqBase_.name_];
//				}

				for(auto & seq : finalSeqs){
					//auto assembleInfo = DefaultAssembleNameInfo(seq->seqBase_.name_, true);
					MetaDataInName seqMeta;
					seqMeta.addMeta("trimmedLength", len(seq->seqBase_));
//					seqMeta.addMeta("estimatedPerBaseCoverage", kmerCoverage[seq->seqBase_.name_]);
					seqMeta.addMeta("estimatedPerBaseCoverage", defaultCoverage);

					seqMeta.addMeta("trimStatus", seq->seqBase_.on_);
					seqMeta.addMeta("regionUID", regionUid);
					seqMeta.addMeta("sample", sample);
					if(finalSeqCounts[seq->seqBase_.name_] > 1){
						seqMeta.addMeta("seqTrimmedCount", finalSeqCountsWritten[seq->seqBase_.name_]);
						++finalSeqCountsWritten[seq->seqBase_.name_];
					}
					//seq->seqBase_.cnt_ = (defaultCoverage/totalCoverage) * (readCounts.pairedReads_ + readCounts.unpaiedReads_ + readCounts.orphans_);
//					std::cout << "kmerCoverage[seq->seqBase_.name_]: " << kmerCoverage[seq->seqBase_.name_] << std::endl;
//					std::cout << "totalCoverage: " << totalCoverage << std::endl;
//					std::cout << "reads: " << readCounts.pairedReads_ + readCounts.unpaiedReads_ + readCounts.orphans_ << std::endl;
					seq->seqBase_.cnt_ = (defaultCoverage/totalCoverage) * (readCounts.pairedReads_ + readCounts.unpaiedReads_ + readCounts.orphans_);
					seqMeta.resetMetaInName(seq->seqBase_.name_);
					seq->seqBase_.name_ += njh::pasteAsStr("_t", seq->seqBase_.cnt_);
				}

				OutOptions trimmedContigInfoOpts(njh::files::make_path(IDBAUDFullOutputDir, "trimmed_reOriented_contigs_outputInfo.tab.txt"));
				OutputStream trimmedContigInfoOut(trimmedContigInfoOpts);
				trimmedContigInfoOut << "name\tlength\tcoverage" << std::endl;
				auto trimmedReOrientedContigsFnp = njh::files::make_path(IDBAUDFullOutputDir, "trimmed_reOriented_contigs.fasta");
				SeqOutput outputWriter(SeqIOOptions::genFastaOut(trimmedReOrientedContigsFnp));
				auto trimmedReOrientedContigsFnp_belowCutOff = njh::files::make_path(IDBAUDFullOutputDir, "trimmed_reOriented_contigs_belowCutOff.fasta");
				SeqOutput belowCutOffOutputWriter(SeqIOOptions::genFastaOut(trimmedReOrientedContigsFnp_belowCutOff));

				uint32_t belowCutOff = 0;
				uint32_t aboveCutOff = 0;
				bool allPassTrim = true;
				for (const auto & contigsKmerRead : finalSeqs) {
					if (len(contigsKmerRead->seqBase_) < minFinalLength) {
						++belowCutOff;
						belowCutOffOutputWriter.openWrite(contigsKmerRead);
						contigsKmerRead->seqBase_.on_ = false;
					} else {
						MetaDataInName seqMeta(contigsKmerRead->seqBase_.name_);
						trimmedContigInfoOut << contigsKmerRead->seqBase_.name_
								<< "\t" << len(contigsKmerRead->seqBase_)
								<< "\t" << seqMeta.getMeta("estimatedPerBaseCoverage")
								<< std::endl;
						if(!contigsKmerRead->seqBase_.on_){
							allPassTrim = false;
						}else{
							++aboveCutOff;
						}
						outputWriter.openWrite(contigsKmerRead);
					}
				}
				if(allPassTrim && !finalSeqs.empty()){
					std::lock_guard<std::mutex> lock(allFinalWriterMut);
					for (const auto & contigsKmerRead : finalSeqs) {
						if (len(contigsKmerRead->seqBase_) >= minFinalLength) {
							allFinalWriter.write(contigsKmerRead);
						}
					}
					for(auto & reg : regInfo){
						reg->infoCalled_ = true;
						reg->uniqHaps_ = aboveCutOff;
					}
				}else{
					std::lock_guard<std::mutex> lock(allPartialWriterMut);
					for (const auto & contigsKmerRead : finalSeqs) {
						if (len(contigsKmerRead->seqBase_) >= minFinalLength) {
							allPartialWriter.write(contigsKmerRead);
						}
					}
					for(auto & reg : regInfo){
						reg->infoCalled_ = false;
						reg->uniqHaps_ = 0;
					}
				}
			} catch (std::exception & e) {
				std::lock_guard<std::mutex> lock(exceptionsMut);
				exceptions[regionUid] = e.what();
				for(auto & reg : regInfo){
					reg->infoCalled_ = false;
					reg->uniqHaps_ = 0;
				}
			}
		}
	};


	njh::concurrent::runVoidFunctionThreaded(runIDBAUDOnRegion, numThreads);
	allFinalWriter.closeOut();
	allPartialWriter.closeOut();
	//sample,readTotal,readTotalUsed, success, name
	//
	OutputStream basicInfo(njh::files::make_path(finalDirectory, "basicInfoPerRegion.tab.txt"));

	basicInfo << "#chrom\tstart\tend\tname\tlength\tstrand\tsuccess\tuniqHaps\treadTotal\treadTotalUsed\ttotalPairedReads";
	basicInfo << "\tsample";
	uint32_t maxExtraFields = 0;
	for(const auto & p : inputRegions){

		auto bedOut = p.genBedRecordCore();
		if(bedOut.extraFields_.size() > maxExtraFields){
			maxExtraFields = bedOut.extraFields_.size();
		}
	}
	for(uint32_t t = 0; t < maxExtraFields; ++t){
		basicInfo << "\textraField"<<t;
	}
	basicInfo << "\n";

	std::map<uint32_t, uint32_t> coiCounts;

	for (const auto & reg : inputRegions) {
		const auto & regInfos = njh::mapAt(regInfosByUID, reg.uid_);
		++coiCounts[regInfos.front()->uniqHaps_];
		for(auto & regInfo : regInfos){
			auto bedOut = regInfo->region_.genBedRecordCore();
			basicInfo << bedOut.toDelimStr();
			basicInfo << "\t" << njh::boolToStr(regInfo->infoCalled_)
								<< "\t" << regInfo->uniqHaps_
								<< "\t" << regInfo->totalReads_
								<< "\t" << regInfo->totalFinalReads_
								<< "\t" << regInfo->totalPairedReads_
								<< "\t" << sample;
			for(const auto & extra : bedOut.extraFields_){
				basicInfo << "\t" << extra;
			}
			basicInfo << std::endl;
		}
	}

	OutputStream coiOut(njh::files::make_path(finalDirectory, "coiCounts.tab.txt"));
	coiOut << "coi\tcount" << std::endl;
	for(const auto & count : coiCounts){
		coiOut << count.first << "\t" << count.second << std::endl;
	}

	OutputStream exceptionsOut(njh::files::make_path(finalDirectory, "exceptionsMessages.tab.txt"));
	exceptionsOut << "regionUID\tmessage" << std::endl;
	for(const auto & exp : exceptions){
		exceptionsOut << exp.first << "\t" << exp.second << std::endl;
	}



	return 0;
}



int programWrappersAssembleOnPathWeaverRunner::runTrinityOnPathWeaverRegions(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bedFile = "";
	bfs::path pwOutputDir = "";
	std::string sample;

	uint32_t TrinityNumThreads = 1;
	uint32_t TrinityMaxMemory = 10;
	std::string extraTrinityOptions;
	uint32_t reOrientingKmerLength = 9;
	uint32_t minFinalLength = 40;
	bfs::path TrinityOutDir = "TrinityOut";
	uint32_t numThreads = 1;
	seqSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processVerbose();
	setUp.setOption(bedFile, "--bed", "The Regions to analyze", true);
	setUp.setOption(pwOutputDir, "--pwOutputDir", "The PathWeaver directory", true);
	setUp.setOption(sample, "--sample", "sample name", true);

	setUp.setOption(numThreads, "--numThreads", "num Threads");
	setUp.setOption(TrinityMaxMemory, "--TrinityMaxMemory", "Trinity Max Memory (in gigabytes)");


	setUp.setOption(TrinityNumThreads, "--TrinityNumThreads", "Trinity Num Threads");
	setUp.setOption(extraTrinityOptions, "--extraTrinityOptions", "extra Trinity Options");

	setUp.setOption(minFinalLength, "--minFinalLength", "min Final Length");
	setUp.setOption(reOrientingKmerLength, "--reOrientingKmerLength", "re-orienting K-mer Length");

	setUp.setOption(TrinityOutDir,     "--TrinityOutDir",     "Trinity Out Directory name, will be relative to final pass directory");


	setUp.processDirectoryOutputName(njh::pasteAsStr(bfs::basename(pwOutputDir), "_Trinity_TODAY"), true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	njh::sys::requireExternalProgramThrow("Trinity");

	auto inputRegions = gatherRegions(bedFile.string(), "", setUp.pars_.verbose_);
	sortGRegionsByStart(inputRegions);

	std::set<std::string> regionNames;
	for(const auto & reg : inputRegions){
		regionNames.emplace(reg.uid_);
	}
	//njh::sort(regionNames);
	njh::concurrent::LockableQueue<std::string> regionsQueue(regionNames);

	bfs::path finalDirectory = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar("final"));
	bfs::path partialDirectory = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar("partial"));
	auto allFinalSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(finalDirectory, "allFinal.fasta"));
	auto allPartialSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(partialDirectory, "allPartial.fasta"));
	SeqOutput allFinalWriter(allFinalSeqOpts);
	SeqOutput allPartialWriter(allPartialSeqOpts);
	allFinalWriter.openOut();
	allPartialWriter.openOut();
	std::mutex allFinalWriterMut;
	std::mutex allPartialWriterMut;


	std::unordered_map<std::string,
			std::vector<std::shared_ptr<BamRegionInvestigator::RegionInfo>> >regInfosByUID;
	for (const auto & reg : inputRegions) {
		regInfosByUID[reg.uid_].emplace_back(std::make_shared<BamRegionInvestigator::RegionInfo>(reg));
	}

	std::unordered_map<std::string, std::string> exceptions;
	std::mutex exceptionsMut;

	std::function<void()> runTrinityOnRegion = [&](){
		std::string regionUid;
		while(regionsQueue.getVal(regionUid)){
			const auto & regInfo = njh::mapAt(regInfosByUID, regionUid);
			auto regionOutputDir = njh::files::make_path(setUp.pars_.directoryName_, regionUid, sample);
			njh::files::makeDirP(njh::files::MkdirPar{regionOutputDir});

			//

			bfs::path refFnp = njh::files::make_path(pwOutputDir, regionUid, "allRefs.fasta");


			//first extract the reads
			bfs::path extractBam = njh::files::make_path(pwOutputDir, regionUid, sample + "_extraction", "extracted.bam");
			OutOptions outOpts(njh::files::make_path(regionOutputDir, "extracted"));
			auto readCounts = rawWriteExtractReadsFromBamOnlyMapped(extractBam, outOpts);
			bfs::path pairedR1 = njh::files::make_path(regionOutputDir, "extracted_R1.fastq");
			bfs::path pairedR2 = njh::files::make_path(regionOutputDir, "extracted_R2.fastq");
			bfs::path singles =  njh::files::make_path(regionOutputDir, "extracted.fastq");
			for(auto & reg : regInfo){
				reg->totalPairedReads_ = readCounts.pairedReads_;
				reg->totalReads_ = readCounts.pairedReads_ + readCounts.unpaiedReads_ + readCounts.orphans_;
				reg->totalFinalReads_ = readCounts.pairedReads_ + readCounts.unpaiedReads_ + readCounts.orphans_;
			}
			try {
				if(!exists(pairedR1) && !exists(singles)){
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", couldn't find " << pairedR1 << " or " << singles << ", need to have at least one of them" << "\n";
					throw std::runtime_error{ss.str()};
				}
				std::stringstream TrinityCmdStream;
				TrinityCmdStream << "cd " << regionOutputDir << " && Trinity ";

				if(exists(pairedR1)){
					if(!exists(pairedR2)){
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ", found: " << pairedR1 << " but cound't find it's mate file: " << pairedR2 << "\n";
						throw std::runtime_error{ss.str()};
					}else{
						TrinityCmdStream << " --left " << pairedR1.filename() << " --right " << pairedR2.filename() << " ";
					}
				}else if(exists(singles)){
					TrinityCmdStream << " --single  " << singles.filename();
				}

				TrinityCmdStream  << " --seqType fq  --CPU " << TrinityNumThreads
											  << " --max_memory " << TrinityMaxMemory << "G"
												<< " --min_contig_length " << minFinalLength
												<< " " << extraTrinityOptions
												<< " --output " << TrinityOutDir
												<< " > TrinityRunLog_" << njh::getCurrentDate() << ".txt 2>&1";
				auto TrinityFullOutputDir = njh::files::make_path(regionOutputDir, TrinityOutDir);

				auto TrinityRunOutput = njh::sys::run({TrinityCmdStream.str()});

				BioCmdsUtils::checkRunOutThrow(TrinityRunOutput, __PRETTY_FUNCTION__);

				OutOptions TrinityRunOutputLogOpts(njh::files::make_path(TrinityFullOutputDir, "TrinityRunOutput.json"));
				OutputStream TrinityRunOutputLogOut(TrinityRunOutputLogOpts);
				TrinityRunOutputLogOut << njh::json::toJson(TrinityRunOutput) << std::endl;

				auto contigsFnp = njh::files::make_path(TrinityFullOutputDir, "Trinity.fasta");
        auto tmpContigsFnp = njh::files::make_path(regionOutputDir, TrinityOutDir.string() + "." + "Trinity.fasta");
        if(bfs::exists(tmpContigsFnp) && !bfs::exists(contigsFnp)){
          bfs::copy_file(tmpContigsFnp, contigsFnp);
        }
				auto contigsSeqIoOpts = SeqIOOptions::genFastaIn(contigsFnp);
				contigsSeqIoOpts.includeWhiteSpaceInName_ = false;
				contigsSeqIoOpts.lowerCaseBases_ = "upper";
				SeqInput contigsReader(contigsSeqIoOpts);
				auto contigsSeqs = contigsReader.readAllReads<seqInfo>();
				std::vector<std::shared_ptr<seqWithKmerInfo>> contigsKmerReads;
				contigsKmerReads.reserve(contigsSeqs.size());
				for (const auto & seq : contigsSeqs) {
					contigsKmerReads.emplace_back(std::make_shared<seqWithKmerInfo>(seq));
				}
				allSetKmers(contigsKmerReads, reOrientingKmerLength, true);

				SeqInput refReader(SeqIOOptions::genFastaIn(refFnp));
				auto refSeqs = refReader.readAllReads<seqInfo>();
				std::vector<std::shared_ptr<seqWithKmerInfo>> refKmerReads;
				std::vector<seqInfo> revComp_refSeqs;
				std::vector<kmerInfo> revComp_refSeqsKInfos;
				refKmerReads.reserve(refSeqs.size());
				for (const auto & seq : refSeqs) {
					refKmerReads.emplace_back(std::make_shared<seqWithKmerInfo>(seq));
				}
				allSetKmers(refKmerReads, reOrientingKmerLength, true);
				for(const auto & rSeq : refSeqs){
					revComp_refSeqs.emplace_back(rSeq);
					revComp_refSeqs.back().reverseComplementRead(false, true);
					revComp_refSeqsKInfos.emplace_back(revComp_refSeqs.back().seq_, 7, false);
				}
				for(const auto & seqKmer : contigsKmerReads) {
					uint32_t forwardWinners = 0;
					uint32_t revWinners = 0;
					for (const auto & refSeq : refKmerReads) {
						auto forDist = refSeq->compareKmers(*seqKmer);
						auto revDist = refSeq->compareKmersRevComp(*seqKmer);
						if (forDist.first < revDist.first) {
							++revWinners;
						} else {
							++forwardWinners;
						}
					}
					if (revWinners > forwardWinners) {
						seqKmer->seqBase_.reverseComplementRead(true, true);
					}
				}

				//sort by sequence length;
				njh::sort(contigsKmerReads, [](const std::shared_ptr<seqWithKmerInfo> & seq1, const std::shared_ptr<seqWithKmerInfo> & seq2){
					return len(seq1->seqBase_) > len(seq2->seqBase_);
				});

				OutOptions contigInfoOpts(njh::files::make_path(TrinityFullOutputDir, "contigs_outputInfo.tab.txt"));
				OutputStream contigInfoOut(contigInfoOpts);
				contigInfoOut << "name\tlength" << std::endl;
				uint32_t defaultCoverage = 10;
				for(const auto & contigsKmerRead : contigsKmerReads){
					//auto assembleInfo = DefaultAssembleNameInfo(contigsKmerRead->seqBase_.name_, true);
					contigInfoOut << contigsKmerRead->seqBase_.name_
							<< "\t" << len(contigsKmerRead->seqBase_)
							<< std::endl;
							//<< "\t" << assembleInfo.coverage_ << std::endl;
				}
				auto reOrientedContigsFnp = njh::files::make_path(TrinityFullOutputDir, "reOriented_contigs.fasta");

				SeqOutput::write(contigsKmerReads, SeqIOOptions::genFastaOut(reOrientedContigsFnp));

				uint64_t maxLen = 0;
				readVec::getMaxLength(refSeqs, maxLen);
				readVec::getMaxLength(contigsKmerReads, maxLen);
				aligner alignerObj(maxLen, gapScoringParameters(5,1,0,0,0,0), substituteMatrix(2,-2), false);
				//alignerObj.processAlnInfoInputNoCheck(njh::files::make_path(resultsDirectory, "trimAlnCache").string(), setUp.pars_.verbose_);
				std::vector<kmerInfo> refSeqsKmerInfos;
				refSeqsKmerInfos.reserve(refSeqs.size());
				for(const auto & input : refSeqs){
					refSeqsKmerInfos.emplace_back(input.seq_, 7, false);
				}
				//readVecTrimmer::trimSeqToRefByGlobalAln(contigsKmerReads, refSeqs, refSeqsKmerInfos, alignerObj);
				std::vector<std::shared_ptr<seqWithKmerInfo>> finalSeqs;
				std::unordered_map<std::string, uint32_t> finalSeqCounts;
				std::unordered_map<std::string, uint32_t> finalSeqCountsWritten;
				for(const auto & seq : contigsKmerReads){
					auto trimmed = readVecTrimmer::trimSeqToRefByGlobalAlnBestNoOverlapIncludeRevComp(seq, refSeqs, revComp_refSeqs, refSeqsKmerInfos, revComp_refSeqsKInfos, alignerObj, false);
					for(auto & trimmedSeq : trimmed){
						bool found = false;
						for(const auto & finalSeq : finalSeqs){
							if(finalSeq->seqBase_.seq_ == trimmedSeq.seq_){
								found = true;
								break;
							}
						}
						if(!found){
							finalSeqs.emplace_back(std::make_shared<seqWithKmerInfo>(trimmedSeq, 7, false));
							++finalSeqCounts[trimmedSeq.name_];
						}
					}
				}
				double totalCoverage = finalSeqs.size() * defaultCoverage;
//				for(auto & seq : finalSeqs){
//					//auto assembleInfo = DefaultAssembleNameInfo(seq->seqBase_.name_, true);
//					totalCoverage += assembleInfo.coverage_;
//				}

				for(auto & seq : finalSeqs){
					//auto assembleInfo = DefaultAssembleNameInfo(seq->seqBase_.name_, true);
					MetaDataInName seqMeta;
					seqMeta.addMeta("trimmedLength", len(seq->seqBase_));
					//seqMeta.addMeta("estimatedPerBaseCoverage", assembleInfo.coverage_);
					seqMeta.addMeta("estimatedPerBaseCoverage", defaultCoverage);
					seqMeta.addMeta("trimStatus", seq->seqBase_.on_);
					seqMeta.addMeta("regionUID", regionUid);
					seqMeta.addMeta("sample", sample);
					if(finalSeqCounts[seq->seqBase_.name_] > 1){
						seqMeta.addMeta("seqTrimmedCount", finalSeqCountsWritten[seq->seqBase_.name_]);
						++finalSeqCountsWritten[seq->seqBase_.name_];
					}
					seqMeta.resetMetaInName(seq->seqBase_.name_);
					seq->seqBase_.cnt_ = (defaultCoverage/totalCoverage) * (readCounts.pairedReads_ + readCounts.unpaiedReads_ + readCounts.orphans_);
					//seq->seqBase_.cnt_ = (assembleInfo.coverage_/totalCoverage) * (readCounts.pairedReads_ + readCounts.unpaiedReads_ + readCounts.orphans_);
					seq->seqBase_.name_ += njh::pasteAsStr("_t", seq->seqBase_.cnt_);
				}

				OutOptions trimmedContigInfoOpts(njh::files::make_path(TrinityFullOutputDir, "trimmed_reOriented_contigs_outputInfo.tab.txt"));
				OutputStream trimmedContigInfoOut(trimmedContigInfoOpts);
				trimmedContigInfoOut << "name\tlength\tcoverage" << std::endl;
				auto trimmedReOrientedContigsFnp = njh::files::make_path(TrinityFullOutputDir, "trimmed_reOriented_contigs.fasta");
				SeqOutput outputWriter(SeqIOOptions::genFastaOut(trimmedReOrientedContigsFnp));
				auto trimmedReOrientedContigsFnp_belowCutOff = njh::files::make_path(TrinityFullOutputDir, "trimmed_reOriented_contigs_belowCutOff.fasta");
				SeqOutput belowCutOffOutputWriter(SeqIOOptions::genFastaOut(trimmedReOrientedContigsFnp_belowCutOff));

				uint32_t belowCutOff = 0;
				uint32_t aboveCutOff = 0;
				bool allPassTrim = true;
				for (const auto & contigsKmerRead : finalSeqs) {
					if (len(contigsKmerRead->seqBase_) < minFinalLength) {
						++belowCutOff;
						belowCutOffOutputWriter.openWrite(contigsKmerRead);
						contigsKmerRead->seqBase_.on_ = false;
					} else {
						MetaDataInName seqMeta(contigsKmerRead->seqBase_.name_);
						trimmedContigInfoOut << contigsKmerRead->seqBase_.name_
								<< "\t" << len(contigsKmerRead->seqBase_)
								<< "\t" << seqMeta.getMeta("estimatedPerBaseCoverage")
								<< std::endl;
						if(!contigsKmerRead->seqBase_.on_){
							allPassTrim = false;
						}else{
							++aboveCutOff;
						}
						outputWriter.openWrite(contigsKmerRead);
					}
				}
				if(allPassTrim && !finalSeqs.empty()){
					std::lock_guard<std::mutex> lock(allFinalWriterMut);
					for (const auto & contigsKmerRead : finalSeqs) {
						if (len(contigsKmerRead->seqBase_) >= minFinalLength) {
							allFinalWriter.write(contigsKmerRead);
						}
					}
					for(auto & reg : regInfo){
						reg->infoCalled_ = true;
						reg->uniqHaps_ = aboveCutOff;
					}
				}else{
					std::lock_guard<std::mutex> lock(allPartialWriterMut);
					for (const auto & contigsKmerRead : finalSeqs) {
						if (len(contigsKmerRead->seqBase_) >= minFinalLength) {
							allPartialWriter.write(contigsKmerRead);
						}
					}
					for(auto & reg : regInfo){
						reg->infoCalled_ = false;
						reg->uniqHaps_ = 0;
					}
				}
			} catch (std::exception & e) {
				std::lock_guard<std::mutex> lock(exceptionsMut);
				exceptions[regionUid] = e.what();
				for(auto & reg : regInfo){
					reg->infoCalled_ = false;
					reg->uniqHaps_ = 0;
				}
			}
		}
	};


	njh::concurrent::runVoidFunctionThreaded(runTrinityOnRegion, numThreads);
	allFinalWriter.closeOut();
	allPartialWriter.closeOut();
	//sample,readTotal,readTotalUsed, success, name
	//
	OutputStream basicInfo(njh::files::make_path(finalDirectory, "basicInfoPerRegion.tab.txt"));

	basicInfo << "#chrom\tstart\tend\tname\tlength\tstrand\tsuccess\tuniqHaps\treadTotal\treadTotalUsed\ttotalPairedReads";
	basicInfo << "\tsample";
	uint32_t maxExtraFields = 0;
	for(const auto & p : inputRegions){

		auto bedOut = p.genBedRecordCore();
		if(bedOut.extraFields_.size() > maxExtraFields){
			maxExtraFields = bedOut.extraFields_.size();
		}
	}
	for(uint32_t t = 0; t < maxExtraFields; ++t){
		basicInfo << "\textraField"<<t;
	}
	basicInfo << "\n";

	std::map<uint32_t, uint32_t> coiCounts;

	for (const auto & reg : inputRegions) {
		const auto & regInfos = njh::mapAt(regInfosByUID, reg.uid_);
		++coiCounts[regInfos.front()->uniqHaps_];
		for(auto & regInfo : regInfos){
			auto bedOut = regInfo->region_.genBedRecordCore();
			basicInfo << bedOut.toDelimStr();
			basicInfo << "\t" << njh::boolToStr(regInfo->infoCalled_)
								<< "\t" << regInfo->uniqHaps_
								<< "\t" << regInfo->totalReads_
								<< "\t" << regInfo->totalFinalReads_
								<< "\t" << regInfo->totalPairedReads_
								<< "\t" << sample;
			for(const auto & extra : bedOut.extraFields_){
				basicInfo << "\t" << extra;
			}
			basicInfo << std::endl;
		}
	}

	OutputStream coiOut(njh::files::make_path(finalDirectory, "coiCounts.tab.txt"));
	coiOut << "coi\tcount" << std::endl;
	for(const auto & count : coiCounts){
		coiOut << count.first << "\t" << count.second << std::endl;
	}

	OutputStream exceptionsOut(njh::files::make_path(finalDirectory, "exceptionsMessages.tab.txt"));
	exceptionsOut << "regionUID\tmessage" << std::endl;
	for(const auto & exp : exceptions){
		exceptionsOut << exp.first << "\t" << exp.second << std::endl;
	}



	return 0;
}




int programWrappersAssembleOnPathWeaverRunner::runMegahitOnPathWeaverRegions(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bedFile = "";
	bfs::path pwOutputDir = "";
	std::string sample;

	uint32_t megahitNumThreads = 1;
	std::string extraMegahitOptions;
	uint32_t reOrientingKmerLength = 9;
	uint32_t minFinalLength = 40;
	bfs::path megahitOutDir = "megahitOut";
	uint32_t numThreads = 1;
	seqSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processVerbose();
	setUp.setOption(bedFile, "--bed", "The Regions to analyze", true);
	setUp.setOption(pwOutputDir, "--pwOutputDir", "The PathWeaver directory", true);
	setUp.setOption(sample, "--sample", "sample name", true);

	setUp.setOption(numThreads, "--numThreads", "num Threads");


	setUp.setOption(megahitNumThreads, "--megahitNumThreads", "megahit Num Threads");
	setUp.setOption(extraMegahitOptions, "--extraMegahitOptions", "extra Megahit Options");

	setUp.setOption(minFinalLength, "--minFinalLength", "min Final Length");
	setUp.setOption(reOrientingKmerLength, "--reOrientingKmerLength", "re-orienting K-mer Length");

	setUp.setOption(megahitOutDir,     "--megahitOutDir",     "megahit Out Directory name, will be relative to final pass directory");


	setUp.processDirectoryOutputName(njh::pasteAsStr(bfs::basename(pwOutputDir), "_megahit_TODAY"), true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	njh::sys::requireExternalProgramThrow("megahit");

	auto inputRegions = gatherRegions(bedFile.string(), "", setUp.pars_.verbose_);
	sortGRegionsByStart(inputRegions);

	std::set<std::string> regionNames;
	for(const auto & reg : inputRegions){
		regionNames.emplace(reg.uid_);
	}
	//njh::sort(regionNames);
	njh::concurrent::LockableQueue<std::string> regionsQueue(regionNames);

	bfs::path finalDirectory = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar("final"));
	bfs::path partialDirectory = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar("partial"));
	auto allFinalSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(finalDirectory, "allFinal.fasta"));
	auto allPartialSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(partialDirectory, "allPartial.fasta"));
	SeqOutput allFinalWriter(allFinalSeqOpts);
	SeqOutput allPartialWriter(allPartialSeqOpts);
	allFinalWriter.openOut();
	allPartialWriter.openOut();
	std::mutex allFinalWriterMut;
	std::mutex allPartialWriterMut;


	std::unordered_map<std::string,
			std::vector<std::shared_ptr<BamRegionInvestigator::RegionInfo>> >regInfosByUID;
	for (const auto & reg : inputRegions) {
		regInfosByUID[reg.uid_].emplace_back(std::make_shared<BamRegionInvestigator::RegionInfo>(reg));
	}

	std::unordered_map<std::string, std::string> exceptions;
	std::mutex exceptionsMut;

	std::function<void()> runMegahitOnRegion = [&](){
		std::string regionUid;
		while(regionsQueue.getVal(regionUid)){
			const auto & regInfo = njh::mapAt(regInfosByUID, regionUid);
			auto regionOutputDir = njh::files::make_path(setUp.pars_.directoryName_, regionUid, sample);
			njh::files::makeDirP(njh::files::MkdirPar{regionOutputDir});

			//

			bfs::path refFnp = njh::files::make_path(pwOutputDir, regionUid, "allRefs.fasta");


			//first extract the reads
			bfs::path extractBam = njh::files::make_path(pwOutputDir, regionUid, sample + "_extraction", "extracted.bam");
			OutOptions outOpts(njh::files::make_path(regionOutputDir, "extracted"));
			auto readCounts = rawWriteExtractReadsFromBamOnlyMapped(extractBam, outOpts);
			bfs::path pairedR1 = njh::files::make_path(regionOutputDir, "extracted_R1.fastq");
			bfs::path pairedR2 = njh::files::make_path(regionOutputDir, "extracted_R2.fastq");
			bfs::path singles =  njh::files::make_path(regionOutputDir, "extracted.fastq");
			for(auto & reg : regInfo){
				reg->totalPairedReads_ = readCounts.pairedReads_;
				reg->totalReads_ = readCounts.pairedReads_ + readCounts.unpaiedReads_ + readCounts.orphans_;
				reg->totalFinalReads_ = readCounts.pairedReads_ + readCounts.unpaiedReads_ + readCounts.orphans_;
			}
			try {
				if(!exists(pairedR1) && !exists(singles)){
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", couldn't find " << pairedR1 << " or " << singles << ", need to have at least one of them" << "\n";
					throw std::runtime_error{ss.str()};
				}
				std::stringstream megahitCmdStream;
				megahitCmdStream << "cd " << regionOutputDir << " && megahit ";

				if(exists(pairedR1)){
					if(!exists(pairedR2)){
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ", found: " << pairedR1 << " but cound't find it's mate file: " << pairedR2 << "\n";
						throw std::runtime_error{ss.str()};
					}else{
						megahitCmdStream << " -1 " << pairedR1.filename() << " -2 " << pairedR2.filename() << " ";
					}
				}else if(exists(singles)){
					megahitCmdStream << " -r  " << singles.filename();
				}
				megahitCmdStream  << " -t " << megahitNumThreads
												<< " " << extraMegahitOptions
												<< " -o " << megahitOutDir
												<< " > megahitRunLog_" << njh::getCurrentDate() << ".txt 2>&1";
				auto megahitFullOutputDir = njh::files::make_path(regionOutputDir, megahitOutDir);

				auto megahitRunOutput = njh::sys::run({megahitCmdStream.str()});

				BioCmdsUtils::checkRunOutThrow(megahitRunOutput, __PRETTY_FUNCTION__);

				OutOptions megahitRunOutputLogOpts(njh::files::make_path(megahitFullOutputDir, "megahitRunOutput.json"));
				OutputStream megahitRunOutputLogOut(megahitRunOutputLogOpts);
				megahitRunOutputLogOut << njh::json::toJson(megahitRunOutput) << std::endl;

				auto contigsFnp = njh::files::make_path(megahitFullOutputDir, "final.contigs.fa");

				auto contigsSeqIoOpts = SeqIOOptions::genFastaIn(contigsFnp);
//				contigsSeqIoOpts.includeWhiteSpaceInName_ = false;
				contigsSeqIoOpts.lowerCaseBases_ = "upper";
				SeqInput contigsReader(contigsSeqIoOpts);
				auto contigsSeqs = contigsReader.readAllReads<seqInfo>();
				std::vector<std::shared_ptr<seqWithKmerInfo>> contigsKmerReads;
				contigsKmerReads.reserve(contigsSeqs.size());
				for (const auto & seq : contigsSeqs) {
					contigsKmerReads.emplace_back(std::make_shared<seqWithKmerInfo>(seq));
				}
				allSetKmers(contigsKmerReads, reOrientingKmerLength, true);

				SeqInput refReader(SeqIOOptions::genFastaIn(refFnp));
				auto refSeqs = refReader.readAllReads<seqInfo>();
				std::vector<std::shared_ptr<seqWithKmerInfo>> refKmerReads;
				std::vector<seqInfo> revComp_refSeqs;
				std::vector<kmerInfo> revComp_refSeqsKInfos;
				refKmerReads.reserve(refSeqs.size());
				for (const auto & seq : refSeqs) {
					refKmerReads.emplace_back(std::make_shared<seqWithKmerInfo>(seq));
				}
				allSetKmers(refKmerReads, reOrientingKmerLength, true);
				for(const auto & rSeq : refSeqs){
					revComp_refSeqs.emplace_back(rSeq);
					revComp_refSeqs.back().reverseComplementRead(false, true);
					revComp_refSeqsKInfos.emplace_back(revComp_refSeqs.back().seq_, 7, false);
				}
				for(const auto & seqKmer : contigsKmerReads) {
					uint32_t forwardWinners = 0;
					uint32_t revWinners = 0;
					for (const auto & refSeq : refKmerReads) {
						auto forDist = refSeq->compareKmers(*seqKmer);
						auto revDist = refSeq->compareKmersRevComp(*seqKmer);
						if (forDist.first < revDist.first) {
							++revWinners;
						} else {
							++forwardWinners;
						}
					}
					if (revWinners > forwardWinners) {
						seqKmer->seqBase_.reverseComplementRead(true, true);
					}
				}

				//sort by sequence length;
				njh::sort(contigsKmerReads, [](const std::shared_ptr<seqWithKmerInfo> & seq1, const std::shared_ptr<seqWithKmerInfo> & seq2){
					return len(seq1->seqBase_) > len(seq2->seqBase_);
				});

				OutOptions contigInfoOpts(njh::files::make_path(megahitFullOutputDir, "contigs_outputInfo.tab.txt"));
				OutputStream contigInfoOut(contigInfoOpts);
				contigInfoOut << "name\tlength\tcoverage" << std::endl;

				for(const auto & contigsKmerRead : contigsKmerReads){
					auto assembleInfo = DefaultAssembleNameInfo(contigsKmerRead->seqBase_.name_, true);
					contigInfoOut << contigsKmerRead->seqBase_.name_
							<< "\t" << len(contigsKmerRead->seqBase_)
							<< "\t" << assembleInfo.coverage_ << std::endl;
				}
				auto reOrientedContigsFnp = njh::files::make_path(megahitFullOutputDir, "reOriented_contigs.fasta");

				SeqOutput::write(contigsKmerReads, SeqIOOptions::genFastaOut(reOrientedContigsFnp));

				uint64_t maxLen = 0;
				readVec::getMaxLength(refSeqs, maxLen);
				readVec::getMaxLength(contigsKmerReads, maxLen);
				aligner alignerObj(maxLen, gapScoringParameters(5,1,0,0,0,0), substituteMatrix(2,-2), false);
				//alignerObj.processAlnInfoInputNoCheck(njh::files::make_path(resultsDirectory, "trimAlnCache").string(), setUp.pars_.verbose_);
				std::vector<kmerInfo> refSeqsKmerInfos;
				refSeqsKmerInfos.reserve(refSeqs.size());
				for(const auto & input : refSeqs){
					refSeqsKmerInfos.emplace_back(input.seq_, 7, false);
				}
				//readVecTrimmer::trimSeqToRefByGlobalAln(contigsKmerReads, refSeqs, refSeqsKmerInfos, alignerObj);
				std::vector<std::shared_ptr<seqWithKmerInfo>> finalSeqs;
				std::unordered_map<std::string, uint32_t> finalSeqCounts;
				std::unordered_map<std::string, uint32_t> finalSeqCountsWritten;
				for(const auto & seq : contigsKmerReads){
					auto trimmed = readVecTrimmer::trimSeqToRefByGlobalAlnBestNoOverlapIncludeRevComp(seq, refSeqs, revComp_refSeqs, refSeqsKmerInfos, revComp_refSeqsKInfos, alignerObj, false);
					for(auto & trimmedSeq : trimmed){
						bool found = false;
						for(const auto & finalSeq : finalSeqs){
							if(finalSeq->seqBase_.seq_ == trimmedSeq.seq_){
								found = true;
								break;
							}
						}
						if(!found){
							finalSeqs.emplace_back(std::make_shared<seqWithKmerInfo>(trimmedSeq, 7, false));
							++finalSeqCounts[trimmedSeq.name_];
						}
					}
				}
				double totalCoverage = 0;
				for(auto & seq : finalSeqs){
					auto assembleInfo = DefaultAssembleNameInfo(seq->seqBase_.name_, true);
					totalCoverage += assembleInfo.coverage_;
				}

				for(auto & seq : finalSeqs){
					auto assembleInfo = DefaultAssembleNameInfo(seq->seqBase_.name_, true);
					MetaDataInName seqMeta;
					seqMeta.addMeta("trimmedLength", len(seq->seqBase_));
					seqMeta.addMeta("estimatedPerBaseCoverage", assembleInfo.coverage_);
					seqMeta.addMeta("trimStatus", seq->seqBase_.on_);
					seqMeta.addMeta("regionUID", regionUid);
					seqMeta.addMeta("sample", sample);
					if(finalSeqCounts[seq->seqBase_.name_] > 1){
						seqMeta.addMeta("seqTrimmedCount", finalSeqCountsWritten[seq->seqBase_.name_]);
						++finalSeqCountsWritten[seq->seqBase_.name_];
					}
					seqMeta.resetMetaInName(seq->seqBase_.name_);
					seq->seqBase_.cnt_ = (assembleInfo.coverage_/totalCoverage) * (readCounts.pairedReads_ + readCounts.unpaiedReads_ + readCounts.orphans_);
					seq->seqBase_.name_ += njh::pasteAsStr("_t", seq->seqBase_.cnt_);
				}

				OutOptions trimmedContigInfoOpts(njh::files::make_path(megahitFullOutputDir, "trimmed_reOriented_contigs_outputInfo.tab.txt"));
				OutputStream trimmedContigInfoOut(trimmedContigInfoOpts);
				trimmedContigInfoOut << "name\tlength\tcoverage" << std::endl;
				auto trimmedReOrientedContigsFnp = njh::files::make_path(megahitFullOutputDir, "trimmed_reOriented_contigs.fasta");
				SeqOutput outputWriter(SeqIOOptions::genFastaOut(trimmedReOrientedContigsFnp));
				auto trimmedReOrientedContigsFnp_belowCutOff = njh::files::make_path(megahitFullOutputDir, "trimmed_reOriented_contigs_belowCutOff.fasta");
				SeqOutput belowCutOffOutputWriter(SeqIOOptions::genFastaOut(trimmedReOrientedContigsFnp_belowCutOff));

				uint32_t belowCutOff = 0;
				uint32_t aboveCutOff = 0;
				bool allPassTrim = true;
				for (const auto & contigsKmerRead : finalSeqs) {
					if (len(contigsKmerRead->seqBase_) < minFinalLength) {
						++belowCutOff;
						belowCutOffOutputWriter.openWrite(contigsKmerRead);
						contigsKmerRead->seqBase_.on_ = false;
					} else {
						MetaDataInName seqMeta(contigsKmerRead->seqBase_.name_);
						trimmedContigInfoOut << contigsKmerRead->seqBase_.name_
								<< "\t" << len(contigsKmerRead->seqBase_)
								<< "\t" << seqMeta.getMeta("estimatedPerBaseCoverage")
								<< std::endl;
						if(!contigsKmerRead->seqBase_.on_){
							allPassTrim = false;
						}else{
							++aboveCutOff;
						}
						outputWriter.openWrite(contigsKmerRead);
					}
				}
				if(allPassTrim && !finalSeqs.empty()){
					std::lock_guard<std::mutex> lock(allFinalWriterMut);
					for (const auto & contigsKmerRead : finalSeqs) {
						if (len(contigsKmerRead->seqBase_) >= minFinalLength) {
							allFinalWriter.write(contigsKmerRead);
						}
					}
					for(auto & reg : regInfo){
						reg->infoCalled_ = true;
						reg->uniqHaps_ = aboveCutOff;
					}
				}else{
					std::lock_guard<std::mutex> lock(allPartialWriterMut);
					for (const auto & contigsKmerRead : finalSeqs) {
						if (len(contigsKmerRead->seqBase_) >= minFinalLength) {
							allPartialWriter.write(contigsKmerRead);
						}
					}
					for(auto & reg : regInfo){
						reg->infoCalled_ = false;
						reg->uniqHaps_ = 0;
					}
				}
			} catch (std::exception & e) {
				std::lock_guard<std::mutex> lock(exceptionsMut);
				exceptions[regionUid] = e.what();
				for(auto & reg : regInfo){
					reg->infoCalled_ = false;
					reg->uniqHaps_ = 0;
				}
			}
		}
	};


	njh::concurrent::runVoidFunctionThreaded(runMegahitOnRegion, numThreads);
	allFinalWriter.closeOut();
	allPartialWriter.closeOut();
	//sample,readTotal,readTotalUsed, success, name
	//
	OutputStream basicInfo(njh::files::make_path(finalDirectory, "basicInfoPerRegion.tab.txt"));

	basicInfo << "#chrom\tstart\tend\tname\tlength\tstrand\tsuccess\tuniqHaps\treadTotal\treadTotalUsed\ttotalPairedReads";
	basicInfo << "\tsample";
	uint32_t maxExtraFields = 0;
	for(const auto & p : inputRegions){

		auto bedOut = p.genBedRecordCore();
		if(bedOut.extraFields_.size() > maxExtraFields){
			maxExtraFields = bedOut.extraFields_.size();
		}
	}
	for(uint32_t t = 0; t < maxExtraFields; ++t){
		basicInfo << "\textraField"<<t;
	}
	basicInfo << "\n";

	std::map<uint32_t, uint32_t> coiCounts;

	for (const auto & reg : inputRegions) {
		const auto & regInfos = njh::mapAt(regInfosByUID, reg.uid_);
		++coiCounts[regInfos.front()->uniqHaps_];
		for(auto & regInfo : regInfos){
			auto bedOut = regInfo->region_.genBedRecordCore();
			basicInfo << bedOut.toDelimStr();
			basicInfo << "\t" << njh::boolToStr(regInfo->infoCalled_)
								<< "\t" << regInfo->uniqHaps_
								<< "\t" << regInfo->totalReads_
								<< "\t" << regInfo->totalFinalReads_
								<< "\t" << regInfo->totalPairedReads_
								<< "\t" << sample;
			for(const auto & extra : bedOut.extraFields_){
				basicInfo << "\t" << extra;
			}
			basicInfo << std::endl;
		}
	}

	OutputStream coiOut(njh::files::make_path(finalDirectory, "coiCounts.tab.txt"));
	coiOut << "coi\tcount" << std::endl;
	for(const auto & count : coiCounts){
		coiOut << count.first << "\t" << count.second << std::endl;
	}

	OutputStream exceptionsOut(njh::files::make_path(finalDirectory, "exceptionsMessages.tab.txt"));
	exceptionsOut << "regionUID\tmessage" << std::endl;
	for(const auto & exp : exceptions){
		exceptionsOut << exp.first << "\t" << exp.second << std::endl;
	}



	return 0;
}




int programWrappersAssembleOnPathWeaverRunner::runSavageOnPathWeaverRegions(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bedFile = "";
	bfs::path pwOutputDir = "";
	std::string sample;

	bfs::path haploconductPath = "/home/hathawan/sourceCodes/savage/HaploConduct/haploconduct";


	uint32_t savageNumThreads = 1;
	std::string extraSavageOptions;
	uint32_t reOrientingKmerLength = 9;
	uint32_t minFinalLength = 40;
	bfs::path savageOutDir = "savageOut";
	uint32_t numThreads = 1;
	seqSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processVerbose();
	setUp.setOption(bedFile, "--bed", "The Regions to analyze", true);
	setUp.setOption(pwOutputDir, "--pwOutputDir", "The PathWeaver directory", true);
	setUp.setOption(sample, "--sample", "sample name", true);

	setUp.setOption(numThreads, "--numThreads", "num Threads");
	setUp.setOption(haploconductPath, "--haploconductPath", "haploconductPath", true);


	setUp.setOption(savageNumThreads, "--savageNumThreads", "savage Num Threads");
	setUp.setOption(extraSavageOptions, "--extraSavageOptions", "extra Savage Options");

	setUp.setOption(minFinalLength, "--minFinalLength", "min Final Length");
	setUp.setOption(reOrientingKmerLength, "--reOrientingKmerLength", "re-orienting K-mer Length");

	setUp.setOption(savageOutDir,     "--savageOutDir",     "savage Out Directory name, will be relative to final pass directory");


	setUp.processDirectoryOutputName(njh::pasteAsStr(bfs::basename(pwOutputDir), "_savage_TODAY"), true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	//njh::sys::requireExternalProgramThrow("savage");

	auto inputRegions = gatherRegions(bedFile.string(), "", setUp.pars_.verbose_);
	sortGRegionsByStart(inputRegions);

	std::set<std::string> regionNames;
	for(const auto & reg : inputRegions){
		regionNames.emplace(reg.uid_);
	}
	//njh::sort(regionNames);
	njh::concurrent::LockableQueue<std::string> regionsQueue(regionNames);

	bfs::path finalDirectory = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar("final"));
	bfs::path partialDirectory = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar("partial"));
	auto allFinalSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(finalDirectory, "allFinal.fasta"));
	auto allPartialSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(partialDirectory, "allPartial.fasta"));
	SeqOutput allFinalWriter(allFinalSeqOpts);
	SeqOutput allPartialWriter(allPartialSeqOpts);
	allFinalWriter.openOut();
	allPartialWriter.openOut();
	std::mutex allFinalWriterMut;
	std::mutex allPartialWriterMut;


	std::unordered_map<std::string,
			std::vector<std::shared_ptr<BamRegionInvestigator::RegionInfo>> >regInfosByUID;
	for (const auto & reg : inputRegions) {
		regInfosByUID[reg.uid_].emplace_back(std::make_shared<BamRegionInvestigator::RegionInfo>(reg));
	}

	std::unordered_map<std::string, std::string> exceptions;
	std::mutex exceptionsMut;

	std::function<void()> runSavageOnRegion = [&](){
		std::string regionUid;
		while(regionsQueue.getVal(regionUid)){
			const auto & regInfo = njh::mapAt(regInfosByUID, regionUid);
			auto regionOutputDir = njh::files::make_path(setUp.pars_.directoryName_, regionUid, sample);
			njh::files::makeDirP(njh::files::MkdirPar{regionOutputDir});

			//

			bfs::path refFnp = njh::files::make_path(pwOutputDir, regionUid, "allRefs.fasta");


			//first extract the reads
			bfs::path extractBam = njh::files::make_path(pwOutputDir, regionUid, sample + "_extraction", "extracted.bam");
			OutOptions outOpts(njh::files::make_path(regionOutputDir, "extracted"));
			auto readCounts = rawWriteExtractReadsFromBamOnlyMapped(extractBam, outOpts, true);
			bfs::path pairedR1 = njh::files::make_path(regionOutputDir, "extracted_R1.fastq");
			bfs::path pairedR2 = njh::files::make_path(regionOutputDir, "extracted_R2.fastq");
			bfs::path singles =  njh::files::make_path(regionOutputDir, "extracted.fastq");
			for(auto & reg : regInfo){
				reg->totalPairedReads_ = readCounts.pairedReads_;
				reg->totalReads_ = readCounts.pairedReads_ + readCounts.unpaiedReads_ + readCounts.orphans_;
				reg->totalFinalReads_ = readCounts.pairedReads_ + readCounts.unpaiedReads_ + readCounts.orphans_;
			}
			try {
				if(!exists(pairedR1) && !exists(singles)){
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", couldn't find " << pairedR1 << " or " << singles << ", need to have at least one of them" << "\n";
					throw std::runtime_error{ss.str()};
				}
				std::stringstream savageCmdStream;
				savageCmdStream << "cd " << regionOutputDir << " && " << haploconductPath << " savage ";

				if(exists(pairedR1)){
					if(!exists(pairedR2)){
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ", found: " << pairedR1 << " but cound't find it's mate file: " << pairedR2 << "\n";
						throw std::runtime_error{ss.str()};
					}else{
						savageCmdStream << " -p1 " << pairedR1.filename() << " -p2 " << pairedR2.filename() << " ";
					}
				}
				if(exists(singles)){
					savageCmdStream << " -s  " << singles.filename();
				}
				savageCmdStream  << " -t " << savageNumThreads
											  << " --split 1 "
												<< " " << extraSavageOptions
												<< " -o " << savageOutDir
												<< " > savageRunLog_" << njh::getCurrentDate() << ".txt 2>&1";

				// ~/sourceCodes/savage/HaploConduct/haploconduct  -p1 extracted_R1.fastq -p2 extracted_R2.fastq -s extracted.fastq -t 10 --split 1


				auto savageFullOutputDir = njh::files::make_path(regionOutputDir, savageOutDir);

				auto savageRunOutput = njh::sys::run({savageCmdStream.str()});

				BioCmdsUtils::checkRunOutThrow(savageRunOutput, __PRETTY_FUNCTION__);

				OutOptions savageRunOutputLogOpts(njh::files::make_path(savageFullOutputDir, "savageRunOutput.json"));
				OutputStream savageRunOutputLogOut(savageRunOutputLogOpts);
				savageRunOutputLogOut << njh::json::toJson(savageRunOutput) << std::endl;

				auto contigsFnp = njh::files::make_path(savageFullOutputDir, "contigs_stage_c.fasta");

				auto contigsSeqIoOpts = SeqIOOptions::genFastaIn(contigsFnp);
//				contigsSeqIoOpts.includeWhiteSpaceInName_ = false;
				contigsSeqIoOpts.lowerCaseBases_ = "upper";
				SeqInput contigsReader(contigsSeqIoOpts);
				auto contigsSeqs = contigsReader.readAllReads<seqInfo>();
				std::vector<std::shared_ptr<seqWithKmerInfo>> contigsKmerReads;
				contigsKmerReads.reserve(contigsSeqs.size());
				for (const auto & seq : contigsSeqs) {
					contigsKmerReads.emplace_back(std::make_shared<seqWithKmerInfo>(seq));
				}
				allSetKmers(contigsKmerReads, reOrientingKmerLength, true);

				SeqInput refReader(SeqIOOptions::genFastaIn(refFnp));
				auto refSeqs = refReader.readAllReads<seqInfo>();
				std::vector<std::shared_ptr<seqWithKmerInfo>> refKmerReads;
				std::vector<seqInfo> revComp_refSeqs;
				std::vector<kmerInfo> revComp_refSeqsKInfos;
				refKmerReads.reserve(refSeqs.size());
				for (const auto & seq : refSeqs) {
					refKmerReads.emplace_back(std::make_shared<seqWithKmerInfo>(seq));
				}
				allSetKmers(refKmerReads, reOrientingKmerLength, true);
				for(const auto & rSeq : refSeqs){
					revComp_refSeqs.emplace_back(rSeq);
					revComp_refSeqs.back().reverseComplementRead(false, true);
					revComp_refSeqsKInfos.emplace_back(revComp_refSeqs.back().seq_, 7, false);
				}
				for(const auto & seqKmer : contigsKmerReads) {
					uint32_t forwardWinners = 0;
					uint32_t revWinners = 0;
					for (const auto & refSeq : refKmerReads) {
						auto forDist = refSeq->compareKmers(*seqKmer);
						auto revDist = refSeq->compareKmersRevComp(*seqKmer);
						if (forDist.first < revDist.first) {
							++revWinners;
						} else {
							++forwardWinners;
						}
					}
					if (revWinners > forwardWinners) {
						seqKmer->seqBase_.reverseComplementRead(true, true);
					}
				}

				//sort by sequence length;
				njh::sort(contigsKmerReads, [](const std::shared_ptr<seqWithKmerInfo> & seq1, const std::shared_ptr<seqWithKmerInfo> & seq2){
					return len(seq1->seqBase_) > len(seq2->seqBase_);
				});

				OutOptions contigInfoOpts(njh::files::make_path(savageFullOutputDir, "contigs_outputInfo.tab.txt"));
				OutputStream contigInfoOut(contigInfoOpts);
				contigInfoOut << "name\tlength" << std::endl;

				for(const auto & contigsKmerRead : contigsKmerReads){
					//auto assembleInfo = DefaultAssembleNameInfo(contigsKmerRead->seqBase_.name_, true);
					contigInfoOut << contigsKmerRead->seqBase_.name_
							<< "\t" << len(contigsKmerRead->seqBase_) << std::endl;
				}
				auto reOrientedContigsFnp = njh::files::make_path(savageFullOutputDir, "reOriented_contigs.fasta");

				SeqOutput::write(contigsKmerReads, SeqIOOptions::genFastaOut(reOrientedContigsFnp));

				uint64_t maxLen = 0;
				readVec::getMaxLength(refSeqs, maxLen);
				readVec::getMaxLength(contigsKmerReads, maxLen);
				aligner alignerObj(maxLen, gapScoringParameters(5,1,0,0,0,0), substituteMatrix(2,-2), false);
				//alignerObj.processAlnInfoInputNoCheck(njh::files::make_path(resultsDirectory, "trimAlnCache").string(), setUp.pars_.verbose_);
				std::vector<kmerInfo> refSeqsKmerInfos;
				refSeqsKmerInfos.reserve(refSeqs.size());
				for(const auto & input : refSeqs){
					refSeqsKmerInfos.emplace_back(input.seq_, 7, false);
				}
				//readVecTrimmer::trimSeqToRefByGlobalAln(contigsKmerReads, refSeqs, refSeqsKmerInfos, alignerObj);
				std::vector<std::shared_ptr<seqWithKmerInfo>> finalSeqs;
				std::unordered_map<std::string, uint32_t> finalSeqCounts;
				std::unordered_map<std::string, uint32_t> finalSeqCountsWritten;
				for(const auto & seq : contigsKmerReads){
					auto trimmed = readVecTrimmer::trimSeqToRefByGlobalAlnBestNoOverlapIncludeRevComp(seq, refSeqs, revComp_refSeqs, refSeqsKmerInfos, revComp_refSeqsKInfos, alignerObj, false);
					for(auto & trimmedSeq : trimmed){
						bool found = false;
						for(const auto & finalSeq : finalSeqs){
							if(finalSeq->seqBase_.seq_ == trimmedSeq.seq_){
								found = true;
								break;
							}
						}
						if(!found){
							finalSeqs.emplace_back(std::make_shared<seqWithKmerInfo>(trimmedSeq, 7, false));
							++finalSeqCounts[trimmedSeq.name_];
						}
					}
				}
				double totalCoverage = finalSeqs.size();

				for(auto & seq : finalSeqs){
					//auto assembleInfo = DefaultAssembleNameInfo(seq->seqBase_.name_, true);
					MetaDataInName seqMeta;
					seqMeta.addMeta("trimmedLength", len(seq->seqBase_));
					seqMeta.addMeta("estimatedPerBaseCoverage", 10);
					seqMeta.addMeta("trimStatus", seq->seqBase_.on_);
					seqMeta.addMeta("regionUID", regionUid);
					seqMeta.addMeta("sample", sample);
					if(finalSeqCounts[seq->seqBase_.name_] > 1){
						seqMeta.addMeta("seqTrimmedCount", finalSeqCountsWritten[seq->seqBase_.name_]);
						++finalSeqCountsWritten[seq->seqBase_.name_];
					}
					seqMeta.resetMetaInName(seq->seqBase_.name_);
					seq->seqBase_.cnt_ = (1/totalCoverage) * (readCounts.pairedReads_ + readCounts.unpaiedReads_ + readCounts.orphans_);
					seq->seqBase_.name_ += njh::pasteAsStr("_t", seq->seqBase_.cnt_);
				}

				OutOptions trimmedContigInfoOpts(njh::files::make_path(savageFullOutputDir, "trimmed_reOriented_contigs_outputInfo.tab.txt"));
				OutputStream trimmedContigInfoOut(trimmedContigInfoOpts);
				trimmedContigInfoOut << "name\tlength\tcoverage" << std::endl;
				auto trimmedReOrientedContigsFnp = njh::files::make_path(savageFullOutputDir, "trimmed_reOriented_contigs.fasta");
				SeqOutput outputWriter(SeqIOOptions::genFastaOut(trimmedReOrientedContigsFnp));
				auto trimmedReOrientedContigsFnp_belowCutOff = njh::files::make_path(savageFullOutputDir, "trimmed_reOriented_contigs_belowCutOff.fasta");
				SeqOutput belowCutOffOutputWriter(SeqIOOptions::genFastaOut(trimmedReOrientedContigsFnp_belowCutOff));

				uint32_t belowCutOff = 0;
				uint32_t aboveCutOff = 0;
				bool allPassTrim = true;
				for (const auto & contigsKmerRead : finalSeqs) {
					if (len(contigsKmerRead->seqBase_) < minFinalLength) {
						++belowCutOff;
						belowCutOffOutputWriter.openWrite(contigsKmerRead);
						contigsKmerRead->seqBase_.on_ = false;
					} else {
						MetaDataInName seqMeta(contigsKmerRead->seqBase_.name_);
						trimmedContigInfoOut << contigsKmerRead->seqBase_.name_
								<< "\t" << len(contigsKmerRead->seqBase_)
								<< "\t" << seqMeta.getMeta("estimatedPerBaseCoverage")
								<< std::endl;
						if(!contigsKmerRead->seqBase_.on_){
							allPassTrim = false;
						}else{
							++aboveCutOff;
						}
						outputWriter.openWrite(contigsKmerRead);
					}
				}
				if(allPassTrim && !finalSeqs.empty()){
					std::lock_guard<std::mutex> lock(allFinalWriterMut);
					for (const auto & contigsKmerRead : finalSeqs) {
						if (len(contigsKmerRead->seqBase_) >= minFinalLength) {
							allFinalWriter.write(contigsKmerRead);
						}
					}
					for(auto & reg : regInfo){
						reg->infoCalled_ = true;
						reg->uniqHaps_ = aboveCutOff;
					}
				}else{
					std::lock_guard<std::mutex> lock(allPartialWriterMut);
					for (const auto & contigsKmerRead : finalSeqs) {
						if (len(contigsKmerRead->seqBase_) >= minFinalLength) {
							allPartialWriter.write(contigsKmerRead);
						}
					}
					for(auto & reg : regInfo){
						reg->infoCalled_ = false;
						reg->uniqHaps_ = 0;
					}
				}
			} catch (std::exception & e) {
				std::lock_guard<std::mutex> lock(exceptionsMut);
				exceptions[regionUid] = e.what();
				for(auto & reg : regInfo){
					reg->infoCalled_ = false;
					reg->uniqHaps_ = 0;
				}
			}
		}
	};


	njh::concurrent::runVoidFunctionThreaded(runSavageOnRegion, numThreads);
	allFinalWriter.closeOut();
	allPartialWriter.closeOut();
	//sample,readTotal,readTotalUsed, success, name
	//
	OutputStream basicInfo(njh::files::make_path(finalDirectory, "basicInfoPerRegion.tab.txt"));

	basicInfo << "#chrom\tstart\tend\tname\tlength\tstrand\tsuccess\tuniqHaps\treadTotal\treadTotalUsed\ttotalPairedReads";
	basicInfo << "\tsample";
	uint32_t maxExtraFields = 0;
	for(const auto & p : inputRegions){

		auto bedOut = p.genBedRecordCore();
		if(bedOut.extraFields_.size() > maxExtraFields){
			maxExtraFields = bedOut.extraFields_.size();
		}
	}
	for(uint32_t t = 0; t < maxExtraFields; ++t){
		basicInfo << "\textraField"<<t;
	}
	basicInfo << "\n";

	std::map<uint32_t, uint32_t> coiCounts;

	for (const auto & reg : inputRegions) {
		const auto & regInfos = njh::mapAt(regInfosByUID, reg.uid_);
		++coiCounts[regInfos.front()->uniqHaps_];
		for(auto & regInfo : regInfos){
			auto bedOut = regInfo->region_.genBedRecordCore();
			basicInfo << bedOut.toDelimStr();
			basicInfo << "\t" << njh::boolToStr(regInfo->infoCalled_)
								<< "\t" << regInfo->uniqHaps_
								<< "\t" << regInfo->totalReads_
								<< "\t" << regInfo->totalFinalReads_
								<< "\t" << regInfo->totalPairedReads_
								<< "\t" << sample;
			for(const auto & extra : bedOut.extraFields_){
				basicInfo << "\t" << extra;
			}
			basicInfo << std::endl;
		}
	}

	OutputStream coiOut(njh::files::make_path(finalDirectory, "coiCounts.tab.txt"));
	coiOut << "coi\tcount" << std::endl;
	for(const auto & count : coiCounts){
		coiOut << count.first << "\t" << count.second << std::endl;
	}

	OutputStream exceptionsOut(njh::files::make_path(finalDirectory, "exceptionsMessages.tab.txt"));
	exceptionsOut << "regionUID\tmessage" << std::endl;
	for(const auto & exp : exceptions){
		exceptionsOut << exp.first << "\t" << exp.second << std::endl;
	}



	return 0;
}


int programWrappersAssembleOnPathWeaverRunner::runPolyteOnPathWeaverRegions(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bedFile = "";
	bfs::path pwOutputDir = "";
	std::string sample;
	uint32_t hardInsertSizeCutOff = 10000;
	uint32_t mapQualityCutOff = 20;
	bfs::path haploconductPath = "/home/hathawan/sourceCodes/savage/HaploConduct/haploconduct";


	uint32_t polyteNumThreads = 1;
	std::string extraSavageOptions;
	uint32_t reOrientingKmerLength = 9;
	uint32_t minFinalLength = 40;
	uint32_t numThreads = 1;
	seqSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processVerbose();
	setUp.setOption(bedFile, "--bed", "The Regions to analyze", true);
	setUp.setOption(pwOutputDir, "--pwOutputDir", "The PathWeaver directory", true);
	setUp.setOption(sample, "--sample", "sample name", true);

	setUp.setOption(numThreads, "--numThreads", "num Threads");
	setUp.setOption(haploconductPath, "--haploconductPath", "haploconductPath", true);


	setUp.setOption(polyteNumThreads, "--polyteNumThreads", "polyte Num Threads");
	setUp.setOption(extraSavageOptions, "--extraSavageOptions", "extra Savage Options");

	setUp.setOption(minFinalLength, "--minFinalLength", "min Final Length");
	setUp.setOption(reOrientingKmerLength, "--reOrientingKmerLength", "re-orienting K-mer Length");

	setUp.processDirectoryOutputName(njh::pasteAsStr(bfs::basename(pwOutputDir), "_polyte_TODAY"), true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	//njh::sys::requireExternalProgramThrow("polyte");

	auto inputRegions = gatherRegions(bedFile.string(), "", setUp.pars_.verbose_);
	sortGRegionsByStart(inputRegions);

	std::set<std::string> regionNames;
	for(const auto & reg : inputRegions){
		regionNames.emplace(reg.uid_);
	}
	//njh::sort(regionNames);
	njh::concurrent::LockableQueue<std::string> regionsQueue(regionNames);

	bfs::path finalDirectory = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar("final"));
	bfs::path partialDirectory = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar("partial"));
	auto allFinalSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(finalDirectory, "allFinal.fasta"));
	auto allPartialSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(partialDirectory, "allPartial.fasta"));
	SeqOutput allFinalWriter(allFinalSeqOpts);
	SeqOutput allPartialWriter(allPartialSeqOpts);
	allFinalWriter.openOut();
	allPartialWriter.openOut();
	std::mutex allFinalWriterMut;
	std::mutex allPartialWriterMut;

	std::unordered_map<std::string, double> coverageForRegion;
	{
		auto basicInfoFnp = njh::files::make_path(pwOutputDir, "final", "basicInfoPerRegion.tab.txt");
		table basicInfo(basicInfoFnp, "\t", true);
		for(const auto & row : basicInfo){
			double perBaseCoverage = njh::StrToNumConverter::stoToNum<double>(row[basicInfo.getColPos("perBaseCoverage")]);
			if("TRUE" == row[basicInfo.getColPos("success")]){
				perBaseCoverage /= njh::StrToNumConverter::stoToNum<double>(row[basicInfo.getColPos("uniqHaps")]);
			}
			coverageForRegion[row[basicInfo.getColPos("name")]] = perBaseCoverage;
		}
	}

	std::unordered_map<std::string,
					std::vector<std::shared_ptr<BamRegionInvestigator::RegionInfo>> >regInfosByUID;
	for (const auto & reg : inputRegions) {
		regInfosByUID[reg.uid_].emplace_back(std::make_shared<BamRegionInvestigator::RegionInfo>(reg));
	}

	std::unordered_map<std::string, std::string> exceptions;
	std::mutex exceptionsMut;

	std::function<void()> runSavageOnRegion = [&](){
		std::string regionUid;
		while(regionsQueue.getVal(regionUid)){
			const auto & regInfo = njh::mapAt(regInfosByUID, regionUid);
			auto regionOutputDir = njh::files::make_path(setUp.pars_.directoryName_, regionUid, sample);
			njh::files::makeDirP(njh::files::MkdirPar{regionOutputDir});

			//

			bfs::path refFnp = njh::files::make_path(pwOutputDir, regionUid, "allRefs.fasta");


			//first extract the reads
			bfs::path extractBam = njh::files::make_path(pwOutputDir, regionUid, sample + "_extraction", "extracted.bam");
			OutOptions outOpts(njh::files::make_path(regionOutputDir, "extracted"));

			// get insert size
			uint32_t insertSize = 0;
			double insertSizeSD = 0;

			{
				BamTools::BamReader bReader;
				bReader.Open(extractBam.string());
				//checkBamOpenThrow(bReader, extractionBam);
				std::vector<uint32_t> currentInsertSizes;
				BamTools::BamAlignment aln;

				while (bReader.GetNextAlignmentCore(aln)) {
					//skip alignments that don't start in this region
					//this way if regions are close to each other it will avoid counting twice
					if (aln.IsMapped() &&
							aln.IsMateMapped() &&
							aln.IsPrimaryAlignment() &&
							aln.IsReverseStrand() != aln.IsMateReverseStrand()) {
						if (aln.RefID == aln.MateRefID) {
							if(std::abs(aln.InsertSize) > hardInsertSizeCutOff || aln.MapQuality < mapQualityCutOff){
								continue;
							}
							currentInsertSizes.emplace_back(std::abs(aln.InsertSize));
						}
					}
				}
				insertSize = vectorMedianRef(currentInsertSizes);
				insertSizeSD = vectorStandardDeviationPop(currentInsertSizes);
			}

			auto readCounts = rawWriteExtractReadsFromBamOnlyMapped(extractBam, outOpts, true);
			bfs::path pairedR1 = njh::files::make_path(regionOutputDir, "extracted_R1.fastq");
			bfs::path pairedR2 = njh::files::make_path(regionOutputDir, "extracted_R2.fastq");
			bfs::path singles =  njh::files::make_path(regionOutputDir, "extracted.fastq");
			for(auto & reg : regInfo){
				reg->totalPairedReads_ = readCounts.pairedReads_;
				reg->totalReads_ = readCounts.pairedReads_ + readCounts.unpaiedReads_ + readCounts.orphans_;
				reg->totalFinalReads_ = readCounts.pairedReads_ + readCounts.unpaiedReads_ + readCounts.orphans_;
			}
			try {
				if(!exists(pairedR1) && !exists(singles)){
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", couldn't find " << pairedR1 << " or " << singles << ", need to have at least one of them" << "\n";
					throw std::runtime_error{ss.str()};
				}
				std::stringstream polyteCmdStream;
				polyteCmdStream << "cd " << regionOutputDir << " && " << haploconductPath << " polyte ";

				if(exists(pairedR1)){
					if(!exists(pairedR2)){
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ", found: " << pairedR1 << " but cound't find it's mate file: " << pairedR2 << "\n";
						throw std::runtime_error{ss.str()};
					}else{
						polyteCmdStream << " -p1 " << pairedR1.filename() << " -p2 " << pairedR2.filename() << " ";
					}
				}
				if(exists(singles)){
					polyteCmdStream << " -s  " << singles.filename();
				}
				polyteCmdStream  << " -t " << polyteNumThreads
				 								 << " --hap_cov " << njh::mapAt(coverageForRegion, regionUid) << " --insert_size " << insertSize << "  --stddev " << insertSizeSD
												 << " " << extraSavageOptions
												 << " > polyteRunLog_" << njh::getCurrentDate() << ".txt 2>&1";

				// ~/sourceCodes/polyte/HaploConduct/haploconduct  -p1 extracted_R1.fastq -p2 extracted_R2.fastq -s extracted.fastq -t 10 --split 1


				auto polyteFullOutputDir = njh::files::make_path(regionOutputDir);

				auto polyteRunOutput = njh::sys::run({polyteCmdStream.str()});

				BioCmdsUtils::checkRunOutThrow(polyteRunOutput, __PRETTY_FUNCTION__);

				OutOptions polyteRunOutputLogOpts(njh::files::make_path(polyteFullOutputDir, "polyteRunOutput.json"));
				OutputStream polyteRunOutputLogOut(polyteRunOutputLogOpts);
				polyteRunOutputLogOut << njh::json::toJson(polyteRunOutput) << std::endl;

				auto contigsFnp = njh::files::make_path(polyteFullOutputDir, "contigs.fasta");

				auto contigsSeqIoOpts = SeqIOOptions::genFastaIn(contigsFnp);
//				contigsSeqIoOpts.includeWhiteSpaceInName_ = false;
				contigsSeqIoOpts.lowerCaseBases_ = "upper";
				SeqInput contigsReader(contigsSeqIoOpts);
				auto contigsSeqs = contigsReader.readAllReads<seqInfo>();
				std::vector<std::shared_ptr<seqWithKmerInfo>> contigsKmerReads;
				contigsKmerReads.reserve(contigsSeqs.size());
				for (const auto & seq : contigsSeqs) {
					contigsKmerReads.emplace_back(std::make_shared<seqWithKmerInfo>(seq));
				}
				allSetKmers(contigsKmerReads, reOrientingKmerLength, true);

				SeqInput refReader(SeqIOOptions::genFastaIn(refFnp));
				auto refSeqs = refReader.readAllReads<seqInfo>();
				std::vector<std::shared_ptr<seqWithKmerInfo>> refKmerReads;
				std::vector<seqInfo> revComp_refSeqs;
				std::vector<kmerInfo> revComp_refSeqsKInfos;
				refKmerReads.reserve(refSeqs.size());
				for (const auto & seq : refSeqs) {
					refKmerReads.emplace_back(std::make_shared<seqWithKmerInfo>(seq));
				}
				allSetKmers(refKmerReads, reOrientingKmerLength, true);
				for(const auto & rSeq : refSeqs){
					revComp_refSeqs.emplace_back(rSeq);
					revComp_refSeqs.back().reverseComplementRead(false, true);
					revComp_refSeqsKInfos.emplace_back(revComp_refSeqs.back().seq_, 7, false);
				}
				for(const auto & seqKmer : contigsKmerReads) {
					uint32_t forwardWinners = 0;
					uint32_t revWinners = 0;
					for (const auto & refSeq : refKmerReads) {
						auto forDist = refSeq->compareKmers(*seqKmer);
						auto revDist = refSeq->compareKmersRevComp(*seqKmer);
						if (forDist.first < revDist.first) {
							++revWinners;
						} else {
							++forwardWinners;
						}
					}
					if (revWinners > forwardWinners) {
						seqKmer->seqBase_.reverseComplementRead(true, true);
					}
				}

				//sort by sequence length;
				njh::sort(contigsKmerReads, [](const std::shared_ptr<seqWithKmerInfo> & seq1, const std::shared_ptr<seqWithKmerInfo> & seq2){
					return len(seq1->seqBase_) > len(seq2->seqBase_);
				});

				OutOptions contigInfoOpts(njh::files::make_path(polyteFullOutputDir, "contigs_outputInfo.tab.txt"));
				OutputStream contigInfoOut(contigInfoOpts);
				contigInfoOut << "name\tlength" << std::endl;

				for(const auto & contigsKmerRead : contigsKmerReads){
					//auto assembleInfo = DefaultAssembleNameInfo(contigsKmerRead->seqBase_.name_, true);
					contigInfoOut << contigsKmerRead->seqBase_.name_
												<< "\t" << len(contigsKmerRead->seqBase_) << std::endl;
				}
				auto reOrientedContigsFnp = njh::files::make_path(polyteFullOutputDir, "reOriented_contigs.fasta");

				SeqOutput::write(contigsKmerReads, SeqIOOptions::genFastaOut(reOrientedContigsFnp));

				uint64_t maxLen = 0;
				readVec::getMaxLength(refSeqs, maxLen);
				readVec::getMaxLength(contigsKmerReads, maxLen);
				aligner alignerObj(maxLen, gapScoringParameters(5,1,0,0,0,0), substituteMatrix(2,-2), false);
				//alignerObj.processAlnInfoInputNoCheck(njh::files::make_path(resultsDirectory, "trimAlnCache").string(), setUp.pars_.verbose_);
				std::vector<kmerInfo> refSeqsKmerInfos;
				refSeqsKmerInfos.reserve(refSeqs.size());
				for(const auto & input : refSeqs){
					refSeqsKmerInfos.emplace_back(input.seq_, 7, false);
				}
				//readVecTrimmer::trimSeqToRefByGlobalAln(contigsKmerReads, refSeqs, refSeqsKmerInfos, alignerObj);
				std::vector<std::shared_ptr<seqWithKmerInfo>> finalSeqs;
				std::unordered_map<std::string, uint32_t> finalSeqCounts;
				std::unordered_map<std::string, uint32_t> finalSeqCountsWritten;
				for(const auto & seq : contigsKmerReads){
					auto trimmed = readVecTrimmer::trimSeqToRefByGlobalAlnBestNoOverlapIncludeRevComp(seq, refSeqs, revComp_refSeqs, refSeqsKmerInfos, revComp_refSeqsKInfos, alignerObj, false);
					for(auto & trimmedSeq : trimmed){
						bool found = false;
						for(const auto & finalSeq : finalSeqs){
							if(finalSeq->seqBase_.seq_ == trimmedSeq.seq_){
								found = true;
								break;
							}
						}
						if(!found){
							finalSeqs.emplace_back(std::make_shared<seqWithKmerInfo>(trimmedSeq, 7, false));
							++finalSeqCounts[trimmedSeq.name_];
						}
					}
				}
				double totalCoverage = finalSeqs.size();

				for(auto & seq : finalSeqs){
					//auto assembleInfo = DefaultAssembleNameInfo(seq->seqBase_.name_, true);
					MetaDataInName seqMeta;
					seqMeta.addMeta("trimmedLength", len(seq->seqBase_));
					seqMeta.addMeta("estimatedPerBaseCoverage", 10);
					seqMeta.addMeta("trimStatus", seq->seqBase_.on_);
					seqMeta.addMeta("regionUID", regionUid);
					seqMeta.addMeta("sample", sample);
					if(finalSeqCounts[seq->seqBase_.name_] > 1){
						seqMeta.addMeta("seqTrimmedCount", finalSeqCountsWritten[seq->seqBase_.name_]);
						++finalSeqCountsWritten[seq->seqBase_.name_];
					}
					seqMeta.resetMetaInName(seq->seqBase_.name_);
					seq->seqBase_.cnt_ = (1/totalCoverage) * (readCounts.pairedReads_ + readCounts.unpaiedReads_ + readCounts.orphans_);
					seq->seqBase_.name_ += njh::pasteAsStr("_t", seq->seqBase_.cnt_);
				}

				OutOptions trimmedContigInfoOpts(njh::files::make_path(polyteFullOutputDir, "trimmed_reOriented_contigs_outputInfo.tab.txt"));
				OutputStream trimmedContigInfoOut(trimmedContigInfoOpts);
				trimmedContigInfoOut << "name\tlength\tcoverage" << std::endl;
				auto trimmedReOrientedContigsFnp = njh::files::make_path(polyteFullOutputDir, "trimmed_reOriented_contigs.fasta");
				SeqOutput outputWriter(SeqIOOptions::genFastaOut(trimmedReOrientedContigsFnp));
				auto trimmedReOrientedContigsFnp_belowCutOff = njh::files::make_path(polyteFullOutputDir, "trimmed_reOriented_contigs_belowCutOff.fasta");
				SeqOutput belowCutOffOutputWriter(SeqIOOptions::genFastaOut(trimmedReOrientedContigsFnp_belowCutOff));

				uint32_t belowCutOff = 0;
				uint32_t aboveCutOff = 0;
				bool allPassTrim = true;
				for (const auto & contigsKmerRead : finalSeqs) {
					if (len(contigsKmerRead->seqBase_) < minFinalLength) {
						++belowCutOff;
						belowCutOffOutputWriter.openWrite(contigsKmerRead);
						contigsKmerRead->seqBase_.on_ = false;
					} else {
						MetaDataInName seqMeta(contigsKmerRead->seqBase_.name_);
						trimmedContigInfoOut << contigsKmerRead->seqBase_.name_
																 << "\t" << len(contigsKmerRead->seqBase_)
																 << "\t" << seqMeta.getMeta("estimatedPerBaseCoverage")
																 << std::endl;
						if(!contigsKmerRead->seqBase_.on_){
							allPassTrim = false;
						}else{
							++aboveCutOff;
						}
						outputWriter.openWrite(contigsKmerRead);
					}
				}
				if(allPassTrim && !finalSeqs.empty()){
					std::lock_guard<std::mutex> lock(allFinalWriterMut);
					for (const auto & contigsKmerRead : finalSeqs) {
						if (len(contigsKmerRead->seqBase_) >= minFinalLength) {
							allFinalWriter.write(contigsKmerRead);
						}
					}
					for(auto & reg : regInfo){
						reg->infoCalled_ = true;
						reg->uniqHaps_ = aboveCutOff;
					}
				}else{
					std::lock_guard<std::mutex> lock(allPartialWriterMut);
					for (const auto & contigsKmerRead : finalSeqs) {
						if (len(contigsKmerRead->seqBase_) >= minFinalLength) {
							allPartialWriter.write(contigsKmerRead);
						}
					}
					for(auto & reg : regInfo){
						reg->infoCalled_ = false;
						reg->uniqHaps_ = 0;
					}
				}
			} catch (std::exception & e) {
				std::lock_guard<std::mutex> lock(exceptionsMut);
				exceptions[regionUid] = e.what();
				for(auto & reg : regInfo){
					reg->infoCalled_ = false;
					reg->uniqHaps_ = 0;
				}
			}
		}
	};


	njh::concurrent::runVoidFunctionThreaded(runSavageOnRegion, numThreads);
	allFinalWriter.closeOut();
	allPartialWriter.closeOut();
	//sample,readTotal,readTotalUsed, success, name
	//
	OutputStream basicInfo(njh::files::make_path(finalDirectory, "basicInfoPerRegion.tab.txt"));

	basicInfo << "#chrom\tstart\tend\tname\tlength\tstrand\tsuccess\tuniqHaps\treadTotal\treadTotalUsed\ttotalPairedReads";
	basicInfo << "\tsample";
	uint32_t maxExtraFields = 0;
	for(const auto & p : inputRegions){

		auto bedOut = p.genBedRecordCore();
		if(bedOut.extraFields_.size() > maxExtraFields){
			maxExtraFields = bedOut.extraFields_.size();
		}
	}
	for(uint32_t t = 0; t < maxExtraFields; ++t){
		basicInfo << "\textraField"<<t;
	}
	basicInfo << "\n";

	std::map<uint32_t, uint32_t> coiCounts;

	for (const auto & reg : inputRegions) {
		const auto & regInfos = njh::mapAt(regInfosByUID, reg.uid_);
		++coiCounts[regInfos.front()->uniqHaps_];
		for(auto & regInfo : regInfos){
			auto bedOut = regInfo->region_.genBedRecordCore();
			basicInfo << bedOut.toDelimStr();
			basicInfo << "\t" << njh::boolToStr(regInfo->infoCalled_)
								<< "\t" << regInfo->uniqHaps_
								<< "\t" << regInfo->totalReads_
								<< "\t" << regInfo->totalFinalReads_
								<< "\t" << regInfo->totalPairedReads_
								<< "\t" << sample;
			for(const auto & extra : bedOut.extraFields_){
				basicInfo << "\t" << extra;
			}
			basicInfo << std::endl;
		}
	}

	OutputStream coiOut(njh::files::make_path(finalDirectory, "coiCounts.tab.txt"));
	coiOut << "coi\tcount" << std::endl;
	for(const auto & count : coiCounts){
		coiOut << count.first << "\t" << count.second << std::endl;
	}

	OutputStream exceptionsOut(njh::files::make_path(finalDirectory, "exceptionsMessages.tab.txt"));
	exceptionsOut << "regionUID\tmessage" << std::endl;
	for(const auto & exp : exceptions){
		exceptionsOut << exp.first << "\t" << exp.second << std::endl;
	}



	return 0;
}




int programWrappersAssembleOnPathWeaverRunner::runVelvetOptimizerAndMetaVelvetOnPathWeaverRegions(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bedFile = "";
	bfs::path pwOutputDir = "";
	std::string sample;

	uint32_t velvetStartKmer = 31;
	uint32_t velvetEndKmer = 71;
	uint32_t velvetKmerStep = 10;
//	std::string optFuncKmer = "n50*Lcon/tbp+log(Lbp)";
//	std::string optFuncCov = "n50*Lcon/tbp+log(Lbp)";

	std::string optFuncKmer = "n50";
	std::string optFuncCov = "n50";

	uint32_t velvetNumOfThreads = 1;
	VecStr optimizerFuncsAvail{"LNbp","Lbp","Lcon","max","n50","ncon","tbp"};

	bool overWriteDir = false;
	bool breakUpAmbigousContigs = true;
	double coverageCutOff = 2.0;

	std::string extraVelvetOptimiserOptions;
	bfs::path VelvetOptimiserOutDir = "VelvetOptimiserOutDir";


	uint32_t reOrientingKmerLength = 9;
	uint32_t minFinalLength = 40;
	uint32_t numThreads = 1;
	seqSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processVerbose();
	setUp.setOption(bedFile, "--bed", "The Regions to analyze", true);
	setUp.setOption(pwOutputDir, "--pwOutputDir", "The PathWeaver directory", true);
	setUp.setOption(sample, "--sample", "sample name", true);

	setUp.setOption(velvetStartKmer, "--velvetStartKmer", "velvet Start Kmer Size");
	setUp.setOption(velvetEndKmer, "--velvetEndKmer", "velvet End Kmer");
	setUp.setOption(velvetKmerStep, "--velvetKmerStep", "velvet Kmer Step");
	setUp.setOption(optFuncKmer, "--optFuncKmer", "optimizer Func Kmer");
	setUp.setOption(optFuncCov, "--optFuncCov", "optimizer Func Coverage");


	bool doNotBreakUpAmbigousContigs = false;
	setUp.setOption(doNotBreakUpAmbigousContigs, "--doBreakUpAmbigousContigs", "Do not Break Up Ambigous Contigs at Ns");
	breakUpAmbigousContigs = !doNotBreakUpAmbigousContigs;

	setUp.setOption(coverageCutOff, "--coverageCutOff", "Don't include these sequences in final output");
//	if(!njh::in(optFuncKmer, optimizerFuncsAvail)){
//		setUp.failed_ = true;
//		setUp.addWarning("Error for --optFuncKmer, value was set as " + optFuncKmer + " but doesn't match available options " + njh::conToStr(optimizerFuncsAvail, ", "));
//	}
//	if(!njh::in(optFuncCov, optimizerFuncsAvail)){
//		setUp.failed_ = true;
//		setUp.addWarning("Error for --optFuncCov, value was set as " + optFuncCov + " but doesn't match available options " + njh::conToStr(optimizerFuncsAvail, ", "));
//	}


	setUp.setOption(minFinalLength, "--minFinalLength", "min Final Length");
	setUp.setOption(reOrientingKmerLength, "--reOrientingKmerLength", "re-orienting K-mer Length");
	setUp.setOption(numThreads, "--numThreads", "number of threads to use");
	setUp.setOption(velvetNumOfThreads, "--velvetNumOfThreads", "velvet Num Of Threads");


	setUp.setOption(extraVelvetOptimiserOptions, "--extraVelvetOptimiserOptions", "Extra options to give to spades");
	setUp.setOption(VelvetOptimiserOutDir,     "--VelvetOptimiserOutDir",     "VelvetOptimiser.pl Out Directory name, will be relative to final pass directory");
	if("VelvetOptimiserOutDir" == VelvetOptimiserOutDir){
		if(!njh::in(optFuncKmer, optimizerFuncsAvail) || !njh::in(optFuncCov, optimizerFuncsAvail)){
			VelvetOptimiserOutDir = VelvetOptimiserOutDir.string() + "_complex";
		}else{
			VelvetOptimiserOutDir = VelvetOptimiserOutDir.string() + "_kfunc" + optFuncKmer + "_covfunc" + optFuncCov;
		}
	}

	setUp.processDirectoryOutputName(njh::pasteAsStr(bfs::basename(pwOutputDir), "_Velvet_TODAY"), true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);

	bfs::path metaVelvetFinalDir = njh::replaceString(setUp.pars_.directoryName_, "Velvet", "MetaVelvet");
	njh::files::MkdirPar metaVelvetFinalDirPar(metaVelvetFinalDir, setUp.pars_.overWriteDir_);
	if(bfs::exists(metaVelvetFinalDir) && setUp.pars_.overWriteDir_){
		njh::files::rmDirForce(metaVelvetFinalDir);
	}
	njh::files::makeDirP(metaVelvetFinalDirPar);

	njh::sys::requireExternalProgramThrow("meta-velvetg");
	njh::sys::requireExternalProgramThrow("velveth");
	njh::sys::requireExternalProgramThrow("velvetg");
	njh::sys::requireExternalProgramThrow("VelvetOptimiser.pl");


	auto inputRegions = gatherRegions(bedFile.string(), "", setUp.pars_.verbose_);
	sortGRegionsByStart(inputRegions);

	std::set<std::string> regionNames;
	for(const auto & reg : inputRegions){
		regionNames.emplace(reg.uid_);
	}
	//njh::sort(regionNames);
	njh::concurrent::LockableQueue<std::string> regionsQueue(regionNames);

	bfs::path finalDirectory = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar("final"));
	bfs::path partialDirectory = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar("partial"));
	auto allFinalSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(finalDirectory, "allFinal.fasta"));
	auto allPartialSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(partialDirectory, "allPartial.fasta"));
	SeqOutput allFinalWriter(allFinalSeqOpts);
	SeqOutput allPartialWriter(allPartialSeqOpts);
	allFinalWriter.openOut();
	allPartialWriter.openOut();
	std::mutex allFinalWriterMut;
	std::mutex allPartialWriterMut;




	bfs::path metav_finalDirectory = njh::files::makeDir(metaVelvetFinalDir, njh::files::MkdirPar("final"));
	bfs::path metav_partialDirectory = njh::files::makeDir(metaVelvetFinalDir, njh::files::MkdirPar("partial"));
	auto metav_allFinalSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(metav_finalDirectory, "allFinal.fasta"));
	auto metav_allPartialSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(metav_partialDirectory, "allPartial.fasta"));
	SeqOutput metav_allFinalWriter(metav_allFinalSeqOpts);
	SeqOutput metav_allPartialWriter(metav_allPartialSeqOpts);
	metav_allFinalWriter.openOut();
	metav_allPartialWriter.openOut();
	std::mutex metav_allFinalWriterMut;
	std::mutex metav_allPartialWriterMut;




	std::unordered_map<std::string,
			std::vector<std::shared_ptr<BamRegionInvestigator::RegionInfo>> >regInfosByUID;

	std::unordered_map<std::string,
			std::vector<std::shared_ptr<BamRegionInvestigator::RegionInfo>> > metav_regInfosByUID;

	for (const auto & reg : inputRegions) {
		regInfosByUID[reg.uid_].emplace_back(std::make_shared<BamRegionInvestigator::RegionInfo>(reg));
		metav_regInfosByUID[reg.uid_].emplace_back(std::make_shared<BamRegionInvestigator::RegionInfo>(reg));
	}

	std::unordered_map<std::string, std::string> exceptions;
	std::unordered_map<std::string, std::string> metav_exceptions;

	std::mutex exceptionsMut;

	std::function<void()> runVelvetsOnRegion = [&](){
		std::string regionUid;
		while(regionsQueue.getVal(regionUid)){
			const auto & metav_regInfo = njh::mapAt(metav_regInfosByUID, regionUid);

			const auto & regInfo = njh::mapAt(regInfosByUID, regionUid);
			auto regionOutputDir = njh::files::make_path(setUp.pars_.directoryName_, regionUid, sample);
			njh::files::makeDirP(njh::files::MkdirPar{regionOutputDir});

			//

			bfs::path refFnp = njh::files::make_path(pwOutputDir, regionUid, "allRefs.fasta");


			//first extract the reads
			bfs::path extractBam = njh::files::make_path(pwOutputDir, regionUid, sample + "_extraction", "extracted.bam");
			OutOptions outOpts(njh::files::make_path(regionOutputDir, "extracted"));
			auto readCounts = rawWriteExtractReadsFromBamOnlyMapped(extractBam, outOpts);
			bfs::path pairedR1 = njh::files::make_path(regionOutputDir, "extracted_R1.fastq");
			bfs::path pairedR2 = njh::files::make_path(regionOutputDir, "extracted_R2.fastq");
			bfs::path singles =  njh::files::make_path(regionOutputDir, "extracted.fastq");

			for(auto & reg : regInfo){
				reg->totalPairedReads_ = readCounts.pairedReads_;
				reg->totalReads_ = readCounts.pairedReads_ + readCounts.unpaiedReads_ + readCounts.orphans_;
				reg->totalFinalReads_ = readCounts.pairedReads_ + readCounts.unpaiedReads_ + readCounts.orphans_;
			}
			for(auto & reg : metav_regInfo){
				reg->totalPairedReads_ = readCounts.pairedReads_;
				reg->totalReads_ = readCounts.pairedReads_ + readCounts.unpaiedReads_ + readCounts.orphans_;
				reg->totalFinalReads_ = readCounts.pairedReads_ + readCounts.unpaiedReads_ + readCounts.orphans_;
			}


			auto vOptFullOutputDir = njh::files::make_path(regionOutputDir, VelvetOptimiserOutDir);

			try {

				if(!exists(pairedR1) && !exists(singles)){
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", couldn't find " << pairedR1 << " or " << singles << ", need to have at least one of them" << "\n";
					throw std::runtime_error{ss.str()};
				}

				std::stringstream vOptCmdStream;
				vOptCmdStream << "cd " << regionOutputDir << " && VelvetOptimiser.pl -f '";
				if(exists(pairedR1)){
					if(!exists(pairedR2)){
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ", found: " << pairedR1 << " but cound't find it's mate file: " << pairedR2 << "\n";
						throw std::runtime_error{ss.str()};
					}else{
						vOptCmdStream << " -fastq -shortPaired -separate " << pairedR1.filename() << " " <<  pairedR2.filename() <<  " ";
					}
				}
				if(exists(singles)){
					vOptCmdStream << " -fastq -short " << singles.filename() << " " ;
				}
				vOptCmdStream << "'";
				//VelvetOptimiser.pl  -x 2 -f ' '  --d withRevComp_optimize_n50

				vOptCmdStream  << " -t " << velvetNumOfThreads
											<< " -s " << velvetStartKmer
											<< " -e " << velvetEndKmer
											<< " -x " << velvetKmerStep
											<< " -optFuncKmer '" << optFuncKmer << "'"
											<< " -optFuncCov '" << optFuncCov << "'"
											<< " " << extraVelvetOptimiserOptions
											<< " --d " << VelvetOptimiserOutDir
											<< " -o '-exp_cov auto -cov_cutoff 300' ";
				std::string vOptCmdPreCovCutOff = vOptCmdStream.str();
				vOptCmdStream
											<< " -m " << coverageCutOff
											<< " > VelvetOptimzerRunLog_" << njh::getCurrentDate() << ".txt 2>&1 ";
				if(exists(vOptFullOutputDir)){
					if(overWriteDir){
						njh::files::rmDirForce(vOptFullOutputDir);
					}else{
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ", error " << vOptFullOutputDir << " already exists, use --overWriteDir to overwrite" << "\n";
						throw std::runtime_error{ss.str()};
					}
				}

				auto vOptRunOutput = njh::sys::run({vOptCmdStream.str()});
				auto firstRunOutputJson = njh::json::toJson(vOptRunOutput);
				bool failedFirstRun = false;
				if(!vOptRunOutput.success_){
					failedFirstRun = true;
					if(exists(vOptFullOutputDir)){
						bfs::rename(vOptFullOutputDir, njh::files::prependFileBasename(vOptFullOutputDir, "failed_" +  njh::getCurrentDate() + "_"));
					}
					std::stringstream vOptCmdStreamAgain;
					vOptCmdStreamAgain << 	vOptCmdPreCovCutOff
							<< " > VelvetOptimzerRunLog_" << njh::getCurrentDate() << ".txt 2>&1";
					vOptRunOutput = njh::sys::run({vOptCmdStreamAgain.str()});
				}

				//failed a second time
				auto secondRunOutputJson = njh::json::toJson(vOptRunOutput);
				bool failedSecondRun = false;
				if (!vOptRunOutput.success_) {
					failedSecondRun = true;
					if (exists(vOptFullOutputDir)) {
						bfs::rename(vOptFullOutputDir,
								njh::files::prependFileBasename(vOptFullOutputDir,
										"failed2_" + njh::getCurrentDate() + "_"));
					}
					std::stringstream vOptCmdStreamAgain;
					vOptCmdStreamAgain << 	vOptCmdPreCovCutOff
							<< " -m " << "-1" << " > VelvetOptimzerRunLog_"
							<< njh::getCurrentDate() << ".txt 2>&1";
					vOptRunOutput = njh::sys::run( { vOptCmdStreamAgain.str() });
				}


				BioCmdsUtils::checkRunOutThrow(vOptRunOutput, __PRETTY_FUNCTION__);
				//setUp.startARunLog(vOptFullOutputDir.string());

				OutOptions spadesRunOutputLogOpts(njh::files::make_path(vOptFullOutputDir, "VelvetOptimiserRunOutput.json"));
				OutputStream spadesRunOutputLogOut(spadesRunOutputLogOpts);
				spadesRunOutputLogOut << njh::json::toJson(vOptRunOutput) << std::endl;

				if(failedFirstRun){
					OutOptions failedVelvetRunOutputLogOpts(njh::files::make_path(vOptFullOutputDir, "failed_VelvetOptimiserRunOutput.json"));
					OutputStream failedVelvetRunOutputLogOut(failedVelvetRunOutputLogOpts);
					failedVelvetRunOutputLogOut << firstRunOutputJson << std::endl;
				}

				if(failedSecondRun){
					OutOptions failedVelvetRunOutputLogOpts(njh::files::make_path(vOptFullOutputDir, "failed_VelvetOptimiserSecondTimeRunOutput.json"));
					OutputStream failedVelvetRunOutputLogOut(failedVelvetRunOutputLogOpts);
					failedVelvetRunOutputLogOut << secondRunOutputJson << std::endl;
				}

				{ //velvet optimizer run
					auto contigsFnp = njh::files::make_path(vOptFullOutputDir, "contigs.fasta");
					{
						//make contigs upper case
						auto originalContigsFnp = njh::files::make_path(vOptFullOutputDir, "contigs.fa");
						auto originalContigsOpts = SeqIOOptions::genFastaInOut(originalContigsFnp, contigsFnp);
						originalContigsOpts.lowerCaseBases_ = "upper";
						SeqInput originalContigsReader(originalContigsOpts);
						auto originalContigsSeqs = originalContigsReader.readAllReads<seqInfo>();
						SeqOutput::write(originalContigsSeqs, originalContigsOpts);
					}

					if (breakUpAmbigousContigs) {
						std::regex pat{"N+"};
						bool mark = true;
						//break up
						auto beforeBreakupContigsFnp = njh::files::make_path(vOptFullOutputDir, "beforeBreakUp_contigs.fasta");
						bfs::copy_file(contigsFnp, beforeBreakupContigsFnp);
						auto contigsOpts = SeqIOOptions::genFastaInOut(beforeBreakupContigsFnp, contigsFnp);
						contigsOpts.out_.overWriteFile_ = true;
						SeqIO contigsReader(contigsOpts);
						auto contigsSeqs = contigsReader.in_.readAllReads<seqInfo>();
						contigsReader.out_.openOut();
						for (const auto & seq : contigsSeqs) {
							auto trimmedSeqs = readVecTrimmer::breakUpSeqOnPat(seq, pat);
							for(auto & trimmedSeq : trimmedSeqs){
								if(mark && 0 != trimmedSeq.start_ && len(seq) != trimmedSeq.end_){
									trimmedSeq.seqBase_.name_.append(njh::pasteAsStr("-s", trimmedSeq.start_, "-e", trimmedSeq.end_));
								}
								contigsReader.write(trimmedSeq.seqBase_);
							}
						}
					}


//					auto contigsOpts = SeqIOOptions::genFastaIn(contigsFnp);
					auto contigsSeqIoOpts = SeqIOOptions::genFastaIn(contigsFnp);
//					contigsSeqIoOpts.includeWhiteSpaceInName_ = false;
					contigsSeqIoOpts.lowerCaseBases_ = "upper";
					SeqInput contigsReader(contigsSeqIoOpts);
					auto contigsSeqs = contigsReader.readAllReads<seqInfo>();

					std::vector<std::shared_ptr<seqWithKmerInfo>> contigsKmerReads;
					contigsKmerReads.reserve(contigsSeqs.size());
					for (const auto & seq : contigsSeqs) {
						contigsKmerReads.emplace_back(std::make_shared<seqWithKmerInfo>(seq));
					}

					allSetKmers(contigsKmerReads, reOrientingKmerLength, true);

					SeqInput refReader(SeqIOOptions::genFastaIn( refFnp ) );
					auto refSeqs = refReader.readAllReads<seqInfo>();
					std::vector<std::shared_ptr<seqWithKmerInfo>> refKmerReads;
					std::vector<seqInfo> revComp_refSeqs;
					std::vector<kmerInfo> revComp_refSeqsKInfos;
					refKmerReads.reserve(refSeqs.size());
					for (const auto & seq : refSeqs) {
						refKmerReads.emplace_back(std::make_shared<seqWithKmerInfo>(seq));
					}
					allSetKmers(refKmerReads, reOrientingKmerLength, true);
					for(const auto & rSeq : refSeqs){
						revComp_refSeqs.emplace_back(rSeq);
						revComp_refSeqs.back().reverseComplementRead(false, true);
						revComp_refSeqsKInfos.emplace_back(revComp_refSeqs.back().seq_, 7, false);
					}
					for(const auto & seqKmer : contigsKmerReads) {
						uint32_t forwardWinners = 0;
						uint32_t revWinners = 0;
						for (const auto & refSeq : refKmerReads) {
							auto forDist = refSeq->compareKmers(*seqKmer);
							auto revDist = refSeq->compareKmersRevComp(*seqKmer);
							if (forDist.first < revDist.first) {
								++revWinners;
							} else {
								++forwardWinners;
							}
						}
						if (revWinners > forwardWinners) {
							seqKmer->seqBase_.reverseComplementRead(true, true);
						}
					}

					//sort by sequence length;
					njh::sort(contigsKmerReads, [](const std::shared_ptr<seqWithKmerInfo> & seq1, const std::shared_ptr<seqWithKmerInfo> & seq2){
						return len(seq1->seqBase_) > len(seq2->seqBase_);
					});




					OutOptions contigInfoOpts(njh::files::make_path(vOptFullOutputDir, "contigs_outputInfo.tab.txt"));
					OutputStream contigInfoOut(contigInfoOpts);
					contigInfoOut << "name\tlength\tcoverage" << std::endl;

					for(const auto & contigsKmerRead : contigsKmerReads){
						auto assembleInfo = DefaultAssembleNameInfo(contigsKmerRead->seqBase_.name_);
						contigInfoOut << contigsKmerRead->seqBase_.name_
								<< "\t" << len(contigsKmerRead->seqBase_)
								<< "\t" << assembleInfo.coverage_ << std::endl;
					}

					auto reOrientedContigsFnp = njh::files::make_path(vOptFullOutputDir, "reOriented_contigs.fasta");

					SeqOutput::write(contigsKmerReads, SeqIOOptions::genFastaOut(reOrientedContigsFnp));

					uint64_t maxLen = 0;
					readVec::getMaxLength(refSeqs, maxLen);
					readVec::getMaxLength(contigsKmerReads, maxLen);
					aligner alignerObj(maxLen, gapScoringParameters(5,1,0,0,0,0), substituteMatrix(2,-2), false);
					//alignerObj.processAlnInfoInputNoCheck(njh::files::make_path(resultsDirectory, "trimAlnCache").string(), setUp.pars_.verbose_);
					std::vector<kmerInfo> refSeqsKmerInfos;
					refSeqsKmerInfos.reserve(refSeqs.size());
					for(const auto & input : refSeqs){
						refSeqsKmerInfos.emplace_back(input.seq_, 7, false);
					}
					//readVecTrimmer::trimSeqToRefByGlobalAln(contigsKmerReads, refSeqs, refSeqsKmerInfos, alignerObj);
					std::vector<std::shared_ptr<seqWithKmerInfo>> finalSeqs;
					std::unordered_map<std::string, uint32_t> finalSeqCounts;
					std::unordered_map<std::string, uint32_t> finalSeqCountsWritten;
					for(const auto & seq : contigsKmerReads){
						auto trimmed = readVecTrimmer::trimSeqToRefByGlobalAlnBestNoOverlapIncludeRevComp(seq, refSeqs, revComp_refSeqs, refSeqsKmerInfos, revComp_refSeqsKInfos, alignerObj, false);
						for(auto & trimmedSeq : trimmed){
							bool found = false;
							for(const auto & finalSeq : finalSeqs){
								if(finalSeq->seqBase_.seq_ == trimmedSeq.seq_){
									found = true;
									break;
								}
							}
							if(!found){
								finalSeqs.emplace_back(std::make_shared<seqWithKmerInfo>(trimmedSeq, 7, false));
								++finalSeqCounts[trimmedSeq.name_];
							}
						}
					}

					double totalCoverage = 0;
					for(auto & seq : finalSeqs){
						auto assembleInfo = DefaultAssembleNameInfo(seq->seqBase_.name_);
						totalCoverage += assembleInfo.coverage_;
					}

					for(auto & seq : finalSeqs){
						auto assembleInfo = DefaultAssembleNameInfo(seq->seqBase_.name_);
						MetaDataInName seqMeta;
						seqMeta.addMeta("trimmedLength", len(seq->seqBase_));
						seqMeta.addMeta("estimatedPerBaseCoverage", assembleInfo.coverage_);
						seqMeta.addMeta("trimStatus", seq->seqBase_.on_);
						seqMeta.addMeta("regionUID", regionUid);
						seqMeta.addMeta("sample", sample);
						if(finalSeqCounts[seq->seqBase_.name_] > 1){
							seqMeta.addMeta("seqTrimmedCount", finalSeqCountsWritten[seq->seqBase_.name_]);
							++finalSeqCountsWritten[seq->seqBase_.name_];
						}
						seqMeta.resetMetaInName(seq->seqBase_.name_);
						seq->seqBase_.cnt_ = (assembleInfo.coverage_/totalCoverage) * (readCounts.pairedReads_ + readCounts.unpaiedReads_ + readCounts.orphans_);
						seq->seqBase_.name_ += njh::pasteAsStr("_t", seq->seqBase_.cnt_);
					}

					OutOptions trimmedContigInfoOpts(njh::files::make_path(vOptFullOutputDir, "trimmed_reOriented_contigs_outputInfo.tab.txt"));
					OutputStream trimmedContigInfoOut(trimmedContigInfoOpts);
					trimmedContigInfoOut << "name\tlength\tcoverage" << std::endl;
					auto trimmedReOrientedContigsFnp = njh::files::make_path(vOptFullOutputDir, "trimmed_reOriented_contigs.fasta");
					SeqOutput outputWriter(SeqIOOptions::genFastaOut(trimmedReOrientedContigsFnp));
					auto trimmedReOrientedContigsFnp_belowCutOff = njh::files::make_path(vOptFullOutputDir, "trimmed_reOriented_contigs_belowCutOff.fasta");
					SeqOutput belowCutOffOutputWriter(SeqIOOptions::genFastaOut(trimmedReOrientedContigsFnp_belowCutOff));

					uint32_t belowCutOff = 0;
					uint32_t aboveCutOff = 0;
					bool allPassTrim = true;
					for (const auto & contigsKmerRead : finalSeqs) {
						MetaDataInName seqMeta(contigsKmerRead->seqBase_.name_);

						if (len(contigsKmerRead->seqBase_) < minFinalLength || seqMeta.getMeta<double>("estimatedPerBaseCoverage") < coverageCutOff) {
							++belowCutOff;
							belowCutOffOutputWriter.openWrite(contigsKmerRead);
							contigsKmerRead->seqBase_.on_ = false;
						} else {

							trimmedContigInfoOut << contigsKmerRead->seqBase_.name_
									<< "\t" << len(contigsKmerRead->seqBase_)
									<< "\t" << seqMeta.getMeta("estimatedPerBaseCoverage")
									<< std::endl;
							if(!contigsKmerRead->seqBase_.on_){
								allPassTrim = false;
							}else{
								++aboveCutOff;
							}
							outputWriter.openWrite(contigsKmerRead);
						}
					}
					if(allPassTrim && !finalSeqs.empty()){
						std::lock_guard<std::mutex> lock(allFinalWriterMut);
						for (const auto & contigsKmerRead : finalSeqs) {
							if (len(contigsKmerRead->seqBase_) >= minFinalLength) {
								allFinalWriter.write(contigsKmerRead);
							}
						}
						for(auto & reg : regInfo){
							reg->infoCalled_ = true;
							reg->uniqHaps_ = aboveCutOff;
						}
					}else{
						std::lock_guard<std::mutex> lock(allPartialWriterMut);
						for (const auto & contigsKmerRead : finalSeqs) {
							if (len(contigsKmerRead->seqBase_) >= minFinalLength) {
								allPartialWriter.write(contigsKmerRead);
							}
						}
						for(auto & reg : regInfo){
							reg->infoCalled_ = false;
							reg->uniqHaps_ = 0;
						}
					}
				}
			} catch (std::exception &e) {
			std::lock_guard<std::mutex> lock(exceptionsMut);
			exceptions[regionUid] = e.what();
			for (auto &reg : regInfo) {
				reg->infoCalled_ = false;
				reg->uniqHaps_ = 0;
			}
			for (auto &reg : metav_regInfo) {
				reg->infoCalled_ = false;
				reg->uniqHaps_ = 0;
			}
		}

			try {



				// now run meta-velvetg
				bfs::path MetaVelvetOutDir = njh::files::make_path(vOptFullOutputDir, "MetaVelvetOutDir");
				njh::files::makeDir(njh::files::MkdirPar{MetaVelvetOutDir});
				auto logfiles = njh::files::gatherFiles(vOptFullOutputDir, "_Logfile.txt");
				if(logfiles.empty() || 1 != logfiles.size()){
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error " << "couldn't determine log file in " << vOptFullOutputDir << "\n";
					if(logfiles.size() > 1){
						ss << "Multiple log files found: " << njh::conToStr(logfiles, ", ") << "\n";
					}
					throw std::runtime_error{ss.str()};
				}
				InputStream logFileIn(InOptions(logfiles.front()));
				std::string line;
				bool encounteredFinalInfo = false;
				std::string velvethArgs;
				std::string velvetgArgs;
				while(njh::files::crossPlatGetline(logFileIn, line)){
					if(njh::beginsWith(line, "Final optimised assembly details")){
						encounteredFinalInfo = true;
					}
					if(encounteredFinalInfo && std::string::npos != line.find(':')){
						auto toks = tokenizeString(line, ":");
						if(toks.size() == 2 && "Velveth parameter string" == toks[0]){
							njh::trim(toks[1]);
							velvethArgs = toks[1];
						}
						if(toks.size() == 2 && "Velvetg parameter string" == toks[0]){
							njh::trim(toks[1]);
							velvetgArgs = njh::replaceString(toks[1], "-clean yes", "");
						}
					}
				}

				if(velvethArgs.empty()){
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error " << "couldn't determine velveth args from " << logfiles.front()<< "\n";
					throw std::runtime_error{ss.str()};
				}

				if(velvetgArgs.empty()){
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error " << "couldn't determine velvetg args from " << logfiles.front()<< "\n";
					throw std::runtime_error{ss.str()};
				}
			//	std::string velvetOutputDirStr = velvetgArgs;
			//	trimAtFirstWhitespace(velvetOutputDirStr);
				std::string velvetOutputDirStr = "optimized_velvet_run";
				bfs::path velvetOutputDir = njh::files::make_path(VelvetOptimiserOutDir, "MetaVelvetOutDir", velvetOutputDirStr);

				std::stringstream velvethCmd;
				velvethCmd << "cd " << vOptFullOutputDir << "/../ && velveth " << velvetOutputDir << " " << velvethArgs.substr(velvethArgs.find(' '))
								<< " > VelvetOptimzer_velveth_RunLog_" << njh::getCurrentDate() << ".txt 2>&1 ";

				std::stringstream velvetgCmd;
				velvetgCmd << "cd " << vOptFullOutputDir << "/../ && velvetg " << velvetOutputDir << " " << velvetgArgs.substr(velvetgArgs.find(' '))
								<< " > VelvetOptimzer_velvetg_RunLog_" << njh::getCurrentDate() << ".txt 2>&1 ";

				std::stringstream metavelvetgCmd;
				metavelvetgCmd << "cd " << vOptFullOutputDir << "/../ && meta-velvetg " << velvetOutputDir
								<< " > VelvetOptimzer_meta-velvetg_RunLog_" << njh::getCurrentDate() << ".txt 2>&1 ";

				auto velvethCmdOutut = njh::sys::run(VecStr{velvethCmd.str()});
				BioCmdsUtils::checkRunOutThrow(velvethCmdOutut, __PRETTY_FUNCTION__);

				auto velvetgCmdOutut = njh::sys::run(VecStr{velvetgCmd.str()});
				BioCmdsUtils::checkRunOutThrow(velvetgCmdOutut, __PRETTY_FUNCTION__);

				auto metavelvetgCmdOutput = njh::sys::run(VecStr{metavelvetgCmd.str()});
				BioCmdsUtils::checkRunOutThrow(metavelvetgCmdOutput, __PRETTY_FUNCTION__);

				Json::Value runOutputs;
				runOutputs["velvethCmd"] = velvethCmdOutut.toJson();
				runOutputs["velvetgCmd"] = velvetgCmdOutut.toJson();
				runOutputs["metavelvetgCmdOutput"] = metavelvetgCmdOutput.toJson();
				{
					OutputStream runOutputsLogOut(njh::files::make_path(vOptFullOutputDir, "MetaVelvetOutDir", velvetOutputDirStr, "runOutputsLog.json"));
					runOutputsLogOut << runOutputs << std::endl;
				}



				{
					auto metaVelvetDir = njh::files::make_path(vOptFullOutputDir, "MetaVelvetOutDir", velvetOutputDirStr);
					bfs::path metaVelvetContigsfile = njh::files::make_path(metaVelvetDir, "meta-velvetg.contigs.fasta");
					{
						//make contigs upper case
						auto originalContigsFnp = njh::files::make_path(metaVelvetDir, "meta-velvetg.contigs.fa");
						auto originalContigsOpts = SeqIOOptions::genFastaInOut(originalContigsFnp, metaVelvetContigsfile);
						originalContigsOpts.lowerCaseBases_ = "upper";
						SeqInput originalContigsReader(originalContigsOpts);
						auto originalContigsSeqs = originalContigsReader.readAllReads<seqInfo>();
						SeqOutput::write(originalContigsSeqs, originalContigsOpts);
					}

					if (breakUpAmbigousContigs) {
						std::regex pat{"N+"};
						bool mark = true;
						//break up
						auto beforeBreakupContigsFnp = njh::files::make_path(metaVelvetDir, "beforeBreakUp_meta-velvetg.contigs.fasta");
						bfs::copy_file(metaVelvetContigsfile, beforeBreakupContigsFnp);
						auto contigsOpts = SeqIOOptions::genFastaInOut(beforeBreakupContigsFnp, metaVelvetContigsfile);
						contigsOpts.out_.overWriteFile_ = true;
						SeqIO contigsReader(contigsOpts);
						auto contigsSeqs = contigsReader.in_.readAllReads<seqInfo>();
						contigsReader.out_.openOut();
						for (const auto & seq : contigsSeqs) {
							auto trimmedSeqs = readVecTrimmer::breakUpSeqOnPat(seq, pat);
							for(auto & trimmedSeq : trimmedSeqs){
								if(mark && 0 != trimmedSeq.start_ && len(seq) != trimmedSeq.end_){
									trimmedSeq.seqBase_.name_.append(njh::pasteAsStr("-s", trimmedSeq.start_, "-e", trimmedSeq.end_));
								}
								contigsReader.write(trimmedSeq.seqBase_);
							}
						}
					}


//					auto contigsOpts = SeqIOOptions::genFastaIn(metaVelvetContigsfile);
					auto contigsSeqIoOpts = SeqIOOptions::genFastaIn(metaVelvetContigsfile);
//					contigsSeqIoOpts.includeWhiteSpaceInName_ = false;
					contigsSeqIoOpts.lowerCaseBases_ = "upper";
					SeqInput contigsReader(contigsSeqIoOpts);
					auto contigsSeqs = contigsReader.readAllReads<seqInfo>();

					std::vector<std::shared_ptr<seqWithKmerInfo>> contigsKmerReads;
					contigsKmerReads.reserve(contigsSeqs.size());
					for (const auto & seq : contigsSeqs) {
						contigsKmerReads.emplace_back(std::make_shared<seqWithKmerInfo>(seq));
					}

					allSetKmers(contigsKmerReads, reOrientingKmerLength, true);

					SeqInput refReader(SeqIOOptions::genFastaIn( refFnp ) );
					auto refSeqs = refReader.readAllReads<seqInfo>();
					std::vector<std::shared_ptr<seqWithKmerInfo>> refKmerReads;
					std::vector<seqInfo> revComp_refSeqs;
					std::vector<kmerInfo> revComp_refSeqsKInfos;
					refKmerReads.reserve(refSeqs.size());
					for (const auto & seq : refSeqs) {
						refKmerReads.emplace_back(std::make_shared<seqWithKmerInfo>(seq));
					}
					allSetKmers(refKmerReads, reOrientingKmerLength, true);
					for(const auto & rSeq : refSeqs){
						revComp_refSeqs.emplace_back(rSeq);
						revComp_refSeqs.back().reverseComplementRead(false, true);
						revComp_refSeqsKInfos.emplace_back(revComp_refSeqs.back().seq_, 7, false);
					}
					for(const auto & seqKmer : contigsKmerReads) {
						uint32_t forwardWinners = 0;
						uint32_t revWinners = 0;
						for (const auto & refSeq : refKmerReads) {
							auto forDist = refSeq->compareKmers(*seqKmer);
							auto revDist = refSeq->compareKmersRevComp(*seqKmer);
							if (forDist.first < revDist.first) {
								++revWinners;
							} else {
								++forwardWinners;
							}
						}
						if (revWinners > forwardWinners) {
							seqKmer->seqBase_.reverseComplementRead(true, true);
						}
					}

					//sort by sequence length;
					njh::sort(contigsKmerReads, [](const std::shared_ptr<seqWithKmerInfo> & seq1, const std::shared_ptr<seqWithKmerInfo> & seq2){
						return len(seq1->seqBase_) > len(seq2->seqBase_);
					});



					OutOptions contigInfoOpts(njh::files::make_path(metaVelvetDir, "meta-velvetg.contigs_outputInfo.tab.txt"));
					OutputStream contigInfoOut(contigInfoOpts);
					contigInfoOut << "name\tlength\tcoverage" << std::endl;

					for(const auto & contigsKmerRead : contigsKmerReads){
						auto assembleInfo = DefaultAssembleNameInfo(contigsKmerRead->seqBase_.name_);
						contigInfoOut << contigsKmerRead->seqBase_.name_
								<< "\t" << len(contigsKmerRead->seqBase_)
								<< "\t" << assembleInfo.coverage_ << std::endl;
					}

					auto reOrientedContigsFnp = njh::files::make_path(metaVelvetDir, "reOriented_meta-velvetg.contigs.fasta");

					SeqOutput::write(contigsKmerReads, SeqIOOptions::genFastaOut(reOrientedContigsFnp));

					uint64_t maxLen = 0;
					readVec::getMaxLength(refSeqs, maxLen);
					readVec::getMaxLength(contigsKmerReads, maxLen);
					aligner alignerObj(maxLen, gapScoringParameters(5,1,0,0,0,0), substituteMatrix(2,-2), false);
					//alignerObj.processAlnInfoInputNoCheck(njh::files::make_path(resultsDirectory, "trimAlnCache").string(), setUp.pars_.verbose_);
					std::vector<kmerInfo> refSeqsKmerInfos;
					refSeqsKmerInfos.reserve(refSeqs.size());
					for(const auto & input : refSeqs){
						refSeqsKmerInfos.emplace_back(input.seq_, 7, false);
					}
					//readVecTrimmer::trimSeqToRefByGlobalAln(contigsKmerReads, refSeqs, refSeqsKmerInfos, alignerObj);
					std::vector<std::shared_ptr<seqWithKmerInfo>> finalSeqs;
					std::unordered_map<std::string, uint32_t> finalSeqCounts;
					std::unordered_map<std::string, uint32_t> finalSeqCountsWritten;
					for(const auto & seq : contigsKmerReads){
						auto trimmed = readVecTrimmer::trimSeqToRefByGlobalAlnBestNoOverlapIncludeRevComp(seq, refSeqs, revComp_refSeqs, refSeqsKmerInfos, revComp_refSeqsKInfos, alignerObj, false);
						for(auto & trimmedSeq : trimmed){
							bool found = false;
							for(const auto & finalSeq : finalSeqs){
								if(finalSeq->seqBase_.seq_ == trimmedSeq.seq_){
									found = true;
									break;
								}
							}
							if(!found){
								finalSeqs.emplace_back(std::make_shared<seqWithKmerInfo>(trimmedSeq, 7, false));
								++finalSeqCounts[trimmedSeq.name_];
							}
						}
					}

					double totalCoverage = 0;
					for(auto & seq : finalSeqs){
						auto assembleInfo = DefaultAssembleNameInfo(seq->seqBase_.name_);
						totalCoverage += assembleInfo.coverage_;
					}

					for(auto & seq : finalSeqs){
						auto assembleInfo = DefaultAssembleNameInfo(seq->seqBase_.name_);
						MetaDataInName seqMeta;
						seqMeta.addMeta("trimmedLength", len(seq->seqBase_));
						seqMeta.addMeta("estimatedPerBaseCoverage", assembleInfo.coverage_);
						seqMeta.addMeta("trimStatus", seq->seqBase_.on_);
						seqMeta.addMeta("regionUID", regionUid);
						seqMeta.addMeta("sample", sample);
						if(finalSeqCounts[seq->seqBase_.name_] > 1){
							seqMeta.addMeta("seqTrimmedCount", finalSeqCountsWritten[seq->seqBase_.name_]);
							++finalSeqCountsWritten[seq->seqBase_.name_];
						}
						seqMeta.resetMetaInName(seq->seqBase_.name_);
						seq->seqBase_.cnt_ = (assembleInfo.coverage_/totalCoverage) * (readCounts.pairedReads_ + readCounts.unpaiedReads_ + readCounts.orphans_);
						seq->seqBase_.name_ += njh::pasteAsStr("_t", seq->seqBase_.cnt_);
					}

					OutOptions trimmedContigInfoOpts(njh::files::make_path(metaVelvetDir, "trimmed_reOriented_meta-velvetg.contigs_outputInfo.tab.txt"));
					OutputStream trimmedContigInfoOut(trimmedContigInfoOpts);
					trimmedContigInfoOut << "name\tlength\tcoverage" << std::endl;
					auto trimmedReOrientedContigsFnp = njh::files::make_path(metaVelvetDir, "trimmed_reOriented_meta-velvetg.contigs.fasta");
					SeqOutput outputWriter(SeqIOOptions::genFastaOut(trimmedReOrientedContigsFnp));
					auto trimmedReOrientedContigsFnp_belowCutOff = njh::files::make_path(metaVelvetDir, "trimmed_reOriented_meta-velvetg.contigs_belowCutOff.fasta");
					SeqOutput belowCutOffOutputWriter(SeqIOOptions::genFastaOut(trimmedReOrientedContigsFnp_belowCutOff));

					uint32_t belowCutOff = 0;
					uint32_t aboveCutOff = 0;
					bool allPassTrim = true;
					for (const auto & contigsKmerRead : finalSeqs) {
						MetaDataInName seqMeta(contigsKmerRead->seqBase_.name_);
						if (len(contigsKmerRead->seqBase_) < minFinalLength || seqMeta.getMeta<double>("estimatedPerBaseCoverage") < coverageCutOff) {
							++belowCutOff;
							belowCutOffOutputWriter.openWrite(contigsKmerRead);
							contigsKmerRead->seqBase_.on_ = false;
						} else {
							trimmedContigInfoOut << contigsKmerRead->seqBase_.name_
									<< "\t" << len(contigsKmerRead->seqBase_)
									<< "\t" << seqMeta.getMeta("estimatedPerBaseCoverage")
									<< std::endl;
							if(!contigsKmerRead->seqBase_.on_){
								allPassTrim = false;
							}else{
								++aboveCutOff;
							}
							outputWriter.openWrite(contigsKmerRead);
						}
					}
					if(allPassTrim && !finalSeqs.empty()){
						std::lock_guard<std::mutex> lock(metav_allFinalWriterMut);
						for (const auto & contigsKmerRead : finalSeqs) {
							if (len(contigsKmerRead->seqBase_) >= minFinalLength) {
								metav_allFinalWriter.write(contigsKmerRead);
							}
						}
						for(auto & reg : metav_regInfo){
							reg->infoCalled_ = true;
							reg->uniqHaps_ = aboveCutOff;
						}
					}else{
						std::lock_guard<std::mutex> lock(metav_allPartialWriterMut);
						for (const auto & contigsKmerRead : finalSeqs) {
							if (len(contigsKmerRead->seqBase_) >= minFinalLength) {
								metav_allPartialWriter.write(contigsKmerRead);
							}
						}
						for(auto & reg : regInfo){
							reg->infoCalled_ = false;
							reg->uniqHaps_ = 0;
						}
					}
				}
			} catch (std::exception & e) {
				std::lock_guard<std::mutex> lock(exceptionsMut);
				metav_exceptions[regionUid] = e.what();
				for (auto &reg : metav_regInfo) {
					reg->infoCalled_ = false;
					reg->uniqHaps_ = 0;
				}
			}
		}
	};


	njh::concurrent::runVoidFunctionThreaded(runVelvetsOnRegion, numThreads);
	allFinalWriter.closeOut();
	allPartialWriter.closeOut();

	metav_allFinalWriter.closeOut();
	metav_allPartialWriter.closeOut();
	//sample,readTotal,readTotalUsed, success, name
	//

	{
		OutputStream basicInfo(njh::files::make_path(finalDirectory, "basicInfoPerRegion.tab.txt"));
		basicInfo << "#chrom\tstart\tend\tname\tlength\tstrand\tsuccess\tuniqHaps\treadTotal\treadTotalUsed\ttotalPairedReads";
		basicInfo << "\tsample";
		uint32_t maxExtraFields = 0;
		for(const auto & p : inputRegions){

			auto bedOut = p.genBedRecordCore();
			if(bedOut.extraFields_.size() > maxExtraFields){
				maxExtraFields = bedOut.extraFields_.size();
			}
		}
		for(uint32_t t = 0; t < maxExtraFields; ++t){
			basicInfo << "\textraField"<<t;
		}
		basicInfo << "\n";

		std::map<uint32_t, uint32_t> coiCounts;

		for (const auto & reg : inputRegions) {
			const auto & regInfos = njh::mapAt(regInfosByUID, reg.uid_);
			++coiCounts[regInfos.front()->uniqHaps_];
			for(auto & regInfo : regInfos){
				auto bedOut = regInfo->region_.genBedRecordCore();
				basicInfo << bedOut.toDelimStr();
				basicInfo << "\t" << njh::boolToStr(regInfo->infoCalled_)
									<< "\t" << regInfo->uniqHaps_
									<< "\t" << regInfo->totalReads_
									<< "\t" << regInfo->totalFinalReads_
									<< "\t" << regInfo->totalPairedReads_
									<< "\t" << sample;
				for(const auto & extra : bedOut.extraFields_){
					basicInfo << "\t" << extra;
				}
				basicInfo << std::endl;
			}
		}

		OutputStream coiOut(njh::files::make_path(finalDirectory, "coiCounts.tab.txt"));
		coiOut << "coi\tcount" << std::endl;
		for(const auto & count : coiCounts){
			coiOut << count.first << "\t" << count.second << std::endl;
		}

		OutputStream exceptionsOut(njh::files::make_path(finalDirectory, "exceptionsMessages.tab.txt"));
		exceptionsOut << "regionUID\tmessage" << std::endl;
		for(const auto & exp : exceptions){
			exceptionsOut << exp.first << "\t" << exp.second << std::endl;
		}
	}

	{
		OutputStream basicInfo(njh::files::make_path(metav_finalDirectory, "basicInfoPerRegion.tab.txt"));
		basicInfo << "#chrom\tstart\tend\tname\tlength\tstrand\tsuccess\tuniqHaps\treadTotal\treadTotalUsed\ttotalPairedReads";
		basicInfo << "\tsample";
		uint32_t maxExtraFields = 0;
		for(const auto & p : inputRegions){

			auto bedOut = p.genBedRecordCore();
			if(bedOut.extraFields_.size() > maxExtraFields){
				maxExtraFields = bedOut.extraFields_.size();
			}
		}
		for(uint32_t t = 0; t < maxExtraFields; ++t){
			basicInfo << "\textraField"<<t;
		}
		basicInfo << "\n";

		std::map<uint32_t, uint32_t> coiCounts;

		for (const auto & reg : inputRegions) {
			const auto & regInfos = njh::mapAt(metav_regInfosByUID, reg.uid_);
			++coiCounts[regInfos.front()->uniqHaps_];
			for(auto & regInfo : regInfos){
				auto bedOut = regInfo->region_.genBedRecordCore();
				basicInfo << bedOut.toDelimStr();
				basicInfo << "\t" << njh::boolToStr(regInfo->infoCalled_)
									<< "\t" << regInfo->uniqHaps_
									<< "\t" << regInfo->totalReads_
									<< "\t" << regInfo->totalFinalReads_
									<< "\t" << regInfo->totalPairedReads_
									<< "\t" << sample;
				for(const auto & extra : bedOut.extraFields_){
					basicInfo << "\t" << extra;
				}
				basicInfo << std::endl;
			}
		}

		OutputStream coiOut(njh::files::make_path(metav_finalDirectory, "coiCounts.tab.txt"));
		coiOut << "coi\tcount" << std::endl;
		for(const auto & count : coiCounts){
			coiOut << count.first << "\t" << count.second << std::endl;
		}

		OutputStream exceptionsOut(njh::files::make_path(metav_finalDirectory, "exceptionsMessages.tab.txt"));
		exceptionsOut << "regionUID\tmessage" << std::endl;
		for(const auto & exp : metav_exceptions){
			exceptionsOut << exp.first << "\t" << exp.second << std::endl;
		}
	}

	return 0;
}






int programWrappersAssembleOnPathWeaverRunner::runPRICEOnPathWeaverRegions(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bedFile = "";
	bfs::path pwOutputDir = "";
	std::string sample;

	uint32_t PRICENumThreads = 1;
	std::string extraPRICEOptions;
	uint32_t reOrientingKmerLength = 9;
	uint32_t minFinalLength = 40;

	/// add options
	uint32_t hardInsertSizeCutOff = 10000;
	uint32_t mapQualityCutOff = 20;
	uint32_t numberOfCycles = 10;
	bfs::path priceOutputDir = "price_output";
	///

	bfs::path PRICEOutDir = "PRICEOut";
	uint32_t numThreads = 1;
	seqSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processVerbose();
	setUp.setOption(bedFile, "--bed", "The Regions to analyze", true);
	setUp.setOption(pwOutputDir, "--pwOutputDir", "The PathWeaver directory", true);
	setUp.setOption(sample, "--sample", "sample name", true);

	setUp.setOption(numThreads, "--numThreads", "num Threads");


	setUp.setOption(PRICENumThreads, "--PRICENumThreads", "PRICE Num Threads");
	setUp.setOption(extraPRICEOptions, "--extraPRICEOptions", "extra PRICE Options");

	setUp.setOption(minFinalLength, "--minFinalLength", "min Final Length");
	setUp.setOption(reOrientingKmerLength, "--reOrientingKmerLength", "re-orienting K-mer Length");

	setUp.setOption(PRICEOutDir,     "--PRICEOutDir",     "PRICE Out Directory name, will be relative to final pass directory");


	setUp.processDirectoryOutputName(njh::pasteAsStr(bfs::basename(pwOutputDir), "_PRICE_TODAY"), true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	njh::sys::requireExternalProgramThrow("PriceTI");

	auto inputRegions = gatherRegions(bedFile.string(), "", setUp.pars_.verbose_);
	sortGRegionsByStart(inputRegions);

	std::set<std::string> regionNames;
	for(const auto & reg : inputRegions){
		regionNames.emplace(reg.uid_);
	}
	//njh::sort(regionNames);
	njh::concurrent::LockableQueue<std::string> regionsQueue(regionNames);

	bfs::path finalDirectory = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar("final"));
	bfs::path partialDirectory = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar("partial"));
	auto allFinalSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(finalDirectory, "allFinal.fasta"));
	auto allPartialSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(partialDirectory, "allPartial.fasta"));
	SeqOutput allFinalWriter(allFinalSeqOpts);
	SeqOutput allPartialWriter(allPartialSeqOpts);
	allFinalWriter.openOut();
	allPartialWriter.openOut();
	std::mutex allFinalWriterMut;
	std::mutex allPartialWriterMut;


	std::unordered_map<std::string,
			std::vector<std::shared_ptr<BamRegionInvestigator::RegionInfo>> >regInfosByUID;
	for (const auto & reg : inputRegions) {
		regInfosByUID[reg.uid_].emplace_back(std::make_shared<BamRegionInvestigator::RegionInfo>(reg));
	}

	std::unordered_map<std::string, std::string> exceptions;
	std::mutex exceptionsMut;

	std::function<void()> runPRICEOnRegion = [&](){
		std::string regionUid;
		while(regionsQueue.getVal(regionUid)){
			const auto & regInfo = njh::mapAt(regInfosByUID, regionUid);
			auto regionOutputDir = njh::files::make_path(setUp.pars_.directoryName_, regionUid, sample);
			njh::files::makeDirP(njh::files::MkdirPar{regionOutputDir});

			//

			bfs::path refFnp = njh::files::make_path(pwOutputDir, regionUid, "allRefs.fasta");


			//first extract the reads
			bfs::path extractBam = njh::files::make_path(pwOutputDir, regionUid, sample + "_extraction", "extracted.bam");
			OutOptions outOpts(njh::files::make_path(regionOutputDir, "extracted"));
			auto readCounts = rawWriteExtractReadsFromBamOnlyMapped(extractBam, outOpts);

			// get insert size
			uint32_t insertSize = 0;
			{

				BamTools::BamReader bReader;
				bReader.Open(extractBam.string());
				//checkBamOpenThrow(bReader, extractionBam);
				std::vector<uint32_t> currentInsertSizes;
				BamTools::BamAlignment aln;

				while (bReader.GetNextAlignmentCore(aln)) {
					//skip alignments that don't start in this region
					//this way if regions are close to each other it will avoid counting twice
					if (aln.IsMapped() &&
							aln.IsMateMapped() &&
							aln.IsPrimaryAlignment() &&
							aln.IsReverseStrand() != aln.IsMateReverseStrand()) {
						if (aln.RefID == aln.MateRefID) {
							if(std::abs(aln.InsertSize) > hardInsertSizeCutOff || aln.MapQuality < mapQualityCutOff){
								continue;
							}
							currentInsertSizes.emplace_back(std::abs(aln.InsertSize));
						}
					}
				}
				insertSize = vectorMedianRef(currentInsertSizes);
			}

			bfs::path pairedR1 = njh::files::make_path(regionOutputDir, "extracted_R1.fastq");
			bfs::path pairedR2 = njh::files::make_path(regionOutputDir, "extracted_R2.fastq");
			bfs::path singles =  njh::files::make_path(regionOutputDir, "extracted.fastq");
			for(auto & reg : regInfo){
				reg->totalPairedReads_ = readCounts.pairedReads_;
				reg->totalReads_ = readCounts.pairedReads_ + readCounts.unpaiedReads_ + readCounts.orphans_;
				reg->totalFinalReads_ = readCounts.pairedReads_ + readCounts.unpaiedReads_ + readCounts.orphans_;
			}
			try {



				if(!exists(pairedR1) && !exists(singles)){
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", couldn't find " << pairedR1 << " or " << singles << ", need to have at least one of them" << "\n";
					throw std::runtime_error{ss.str()};
				}
				std::stringstream PriceTICmdStream;
				PriceTICmdStream << "cd " << regionOutputDir << " && PriceTI ";
				if (!exists(pairedR1)){
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error " << "PriceTI only works on paired end reads"<< "\n";
					throw std::runtime_error{ss.str()};
				}
				if(numThreads > 1){
					PriceTICmdStream << " -a " << PRICENumThreads << " ";
				}
				if (!exists(pairedR2)) {
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", found: " << pairedR1
							<< " but cound't find it's mate file: " << pairedR2 << "\n";
					throw std::runtime_error { ss.str() };
				} else {
					PriceTICmdStream
							<< " -fpp " << pairedR1.filename() << " " << pairedR2.filename() << " " << insertSize << " 95 ";
					PriceTICmdStream << " -nc " << numberOfCycles;
				}
				if(exists(singles)){
					PriceTICmdStream << " -icf " << singles.filename() << " 2 1 1 ";
				}else{
					PriceTICmdStream << " -picf 400 " << pairedR1.filename() << " 1 1 1 -picf 400 " << pairedR2.filename() << " 1 1 1 ";
				}
				PriceTICmdStream << " -o " << priceOutputDir << "/price_out.fasta";



				auto PRICEFullOutputDir = njh::files::make_path(regionOutputDir, priceOutputDir);

				PriceTICmdStream
												<< " > " << priceOutputDir << "/priceOutputDirRunLog_" << njh::getCurrentDate() << ".txt 2>&1";
				njh::files::makeDir(njh::files::MkdirPar{PRICEFullOutputDir});

				auto PRICERunOutput = njh::sys::run({PriceTICmdStream.str()});

				BioCmdsUtils::checkRunOutThrow(PRICERunOutput, __PRETTY_FUNCTION__);

				OutOptions PRICERunOutputLogOpts(njh::files::make_path(PRICEFullOutputDir, "PRICERunOutput.json"));
				OutputStream PRICERunOutputLogOut(PRICERunOutputLogOpts);
				PRICERunOutputLogOut << njh::json::toJson(PRICERunOutput) << std::endl;

				auto contigsFnp = njh::files::make_path(PRICEFullOutputDir, njh::pasteAsStr("price_out.cycle", numberOfCycles, ".fasta"));

				auto contigsSeqIoOpts = SeqIOOptions::genFastaIn(contigsFnp);
//				contigsSeqIoOpts.includeWhiteSpaceInName_ = false;
				contigsSeqIoOpts.lowerCaseBases_ = "upper";
				SeqInput contigsReader(contigsSeqIoOpts);
				auto contigsSeqs = contigsReader.readAllReads<seqInfo>();
				std::vector<std::shared_ptr<seqWithKmerInfo>> contigsKmerReads;
				contigsKmerReads.reserve(contigsSeqs.size());
				for (const auto & seq : contigsSeqs) {
					contigsKmerReads.emplace_back(std::make_shared<seqWithKmerInfo>(seq));
				}
				allSetKmers(contigsKmerReads, reOrientingKmerLength, true);

				SeqInput refReader(SeqIOOptions::genFastaIn(refFnp));
				auto refSeqs = refReader.readAllReads<seqInfo>();
				std::vector<std::shared_ptr<seqWithKmerInfo>> refKmerReads;
				std::vector<seqInfo> revComp_refSeqs;
				std::vector<kmerInfo> revComp_refSeqsKInfos;
				refKmerReads.reserve(refSeqs.size());
				for (const auto & seq : refSeqs) {
					refKmerReads.emplace_back(std::make_shared<seqWithKmerInfo>(seq));
				}
				allSetKmers(refKmerReads, reOrientingKmerLength, true);
				for(const auto & rSeq : refSeqs){
					revComp_refSeqs.emplace_back(rSeq);
					revComp_refSeqs.back().reverseComplementRead(false, true);
					revComp_refSeqsKInfos.emplace_back(revComp_refSeqs.back().seq_, 7, false);
				}
				for(const auto & seqKmer : contigsKmerReads) {
					uint32_t forwardWinners = 0;
					uint32_t revWinners = 0;
					for (const auto & refSeq : refKmerReads) {
						auto forDist = refSeq->compareKmers(*seqKmer);
						auto revDist = refSeq->compareKmersRevComp(*seqKmer);
						if (forDist.first < revDist.first) {
							++revWinners;
						} else {
							++forwardWinners;
						}
					}
					if (revWinners > forwardWinners) {
						seqKmer->seqBase_.reverseComplementRead(true, true);
					}
				}

				//sort by sequence length;
				njh::sort(contigsKmerReads, [](const std::shared_ptr<seqWithKmerInfo> & seq1, const std::shared_ptr<seqWithKmerInfo> & seq2){
					return len(seq1->seqBase_) > len(seq2->seqBase_);
				});

				OutOptions contigInfoOpts(njh::files::make_path(PRICEFullOutputDir, "contigs_outputInfo.tab.txt"));
				OutputStream contigInfoOut(contigInfoOpts);
				contigInfoOut << "name\tlength" << std::endl;
				uint32_t defaultCoverage = 10;

				for(const auto & contigsKmerRead : contigsKmerReads){
					//auto assembleInfo = DefaultAssembleNameInfo(contigsKmerRead->seqBase_.name_, true);
					contigInfoOut << contigsKmerRead->seqBase_.name_
							<< "\t" << len(contigsKmerRead->seqBase_)
							<< std::endl;
				}
				auto reOrientedContigsFnp = njh::files::make_path(PRICEFullOutputDir, "reOriented_contigs.fasta");

				SeqOutput::write(contigsKmerReads, SeqIOOptions::genFastaOut(reOrientedContigsFnp));

				uint64_t maxLen = 0;
				readVec::getMaxLength(refSeqs, maxLen);
				readVec::getMaxLength(contigsKmerReads, maxLen);
				aligner alignerObj(maxLen, gapScoringParameters(5,1,0,0,0,0), substituteMatrix(2,-2), false);
				//alignerObj.processAlnInfoInputNoCheck(njh::files::make_path(resultsDirectory, "trimAlnCache").string(), setUp.pars_.verbose_);
				std::vector<kmerInfo> refSeqsKmerInfos;
				refSeqsKmerInfos.reserve(refSeqs.size());
				for(const auto & input : refSeqs){
					refSeqsKmerInfos.emplace_back(input.seq_, 7, false);
				}
				//readVecTrimmer::trimSeqToRefByGlobalAln(contigsKmerReads, refSeqs, refSeqsKmerInfos, alignerObj);
				std::vector<std::shared_ptr<seqWithKmerInfo>> finalSeqs;
				std::unordered_map<std::string, uint32_t> finalSeqCounts;
				std::unordered_map<std::string, uint32_t> finalSeqCountsWritten;
				for(const auto & seq : contigsKmerReads){
					auto trimmed = readVecTrimmer::trimSeqToRefByGlobalAlnBestNoOverlapIncludeRevComp(seq, refSeqs, revComp_refSeqs, refSeqsKmerInfos, revComp_refSeqsKInfos, alignerObj, false);
					for(auto & trimmedSeq : trimmed){
						bool found = false;
						for(const auto & finalSeq : finalSeqs){
							if(finalSeq->seqBase_.seq_ == trimmedSeq.seq_){
								found = true;
								break;
							}
						}
						if(!found){
							finalSeqs.emplace_back(std::make_shared<seqWithKmerInfo>(trimmedSeq, 7, false));
							++finalSeqCounts[trimmedSeq.name_];
						}
					}
				}
				double totalCoverage = finalSeqs.size() * defaultCoverage;
//				for(auto & seq : finalSeqs){
//					//auto assembleInfo = DefaultAssembleNameInfo(seq->seqBase_.name_, true);
//					//totalCoverage += 10;
//				}

				for(auto & seq : finalSeqs){
					//auto assembleInfo = DefaultAssembleNameInfo(seq->seqBase_.name_, true);
					MetaDataInName seqMeta;
					seqMeta.addMeta("trimmedLength", len(seq->seqBase_));
					seqMeta.addMeta("estimatedPerBaseCoverage", defaultCoverage);
					seqMeta.addMeta("trimStatus", seq->seqBase_.on_);
					seqMeta.addMeta("regionUID", regionUid);
					seqMeta.addMeta("sample", sample);
					if(finalSeqCounts[seq->seqBase_.name_] > 1){
						seqMeta.addMeta("seqTrimmedCount", finalSeqCountsWritten[seq->seqBase_.name_]);
						++finalSeqCountsWritten[seq->seqBase_.name_];
					}
					seqMeta.resetMetaInName(seq->seqBase_.name_);
					seq->seqBase_.cnt_ = (defaultCoverage/totalCoverage) * (readCounts.pairedReads_ + readCounts.unpaiedReads_ + readCounts.orphans_);
					seq->seqBase_.name_ += njh::pasteAsStr("_t", seq->seqBase_.cnt_);
				}

				OutOptions trimmedContigInfoOpts(njh::files::make_path(PRICEFullOutputDir, "trimmed_reOriented_contigs_outputInfo.tab.txt"));
				OutputStream trimmedContigInfoOut(trimmedContigInfoOpts);
				trimmedContigInfoOut << "name\tlength\tcoverage" << std::endl;
				auto trimmedReOrientedContigsFnp = njh::files::make_path(PRICEFullOutputDir, "trimmed_reOriented_contigs.fasta");
				SeqOutput outputWriter(SeqIOOptions::genFastaOut(trimmedReOrientedContigsFnp));
				auto trimmedReOrientedContigsFnp_belowCutOff = njh::files::make_path(PRICEFullOutputDir, "trimmed_reOriented_contigs_belowCutOff.fasta");
				SeqOutput belowCutOffOutputWriter(SeqIOOptions::genFastaOut(trimmedReOrientedContigsFnp_belowCutOff));

				uint32_t belowCutOff = 0;
				uint32_t aboveCutOff = 0;
				bool allPassTrim = true;
				for (const auto & contigsKmerRead : finalSeqs) {
					if (len(contigsKmerRead->seqBase_) < minFinalLength) {
						++belowCutOff;
						belowCutOffOutputWriter.openWrite(contigsKmerRead);
						contigsKmerRead->seqBase_.on_ = false;
					} else {
						MetaDataInName seqMeta(contigsKmerRead->seqBase_.name_);
						trimmedContigInfoOut << contigsKmerRead->seqBase_.name_
								<< "\t" << len(contigsKmerRead->seqBase_)
								<< "\t" << seqMeta.getMeta("estimatedPerBaseCoverage")
								<< std::endl;
						if(!contigsKmerRead->seqBase_.on_){
							allPassTrim = false;
						}else{
							++aboveCutOff;
						}
						outputWriter.openWrite(contigsKmerRead);
					}
				}
				if(allPassTrim && !finalSeqs.empty()){
					std::lock_guard<std::mutex> lock(allFinalWriterMut);
					for (const auto & contigsKmerRead : finalSeqs) {
						if (len(contigsKmerRead->seqBase_) >= minFinalLength) {
							allFinalWriter.write(contigsKmerRead);
						}
					}
					for(auto & reg : regInfo){
						reg->infoCalled_ = true;
						reg->uniqHaps_ = aboveCutOff;
					}
				}else{
					std::lock_guard<std::mutex> lock(allPartialWriterMut);
					for (const auto & contigsKmerRead : finalSeqs) {
						if (len(contigsKmerRead->seqBase_) >= minFinalLength) {
							allPartialWriter.write(contigsKmerRead);
						}
					}
					for(auto & reg : regInfo){
						reg->infoCalled_ = false;
						reg->uniqHaps_ = 0;
					}
				}
			} catch (std::exception & e) {
				std::lock_guard<std::mutex> lock(exceptionsMut);
				exceptions[regionUid] = e.what();
				for(auto & reg : regInfo){
					reg->infoCalled_ = false;
					reg->uniqHaps_ = 0;
				}
			}
		}
	};


	njh::concurrent::runVoidFunctionThreaded(runPRICEOnRegion, numThreads);
	allFinalWriter.closeOut();
	allPartialWriter.closeOut();
	//sample,readTotal,readTotalUsed, success, name
	//
	OutputStream basicInfo(njh::files::make_path(finalDirectory, "basicInfoPerRegion.tab.txt"));

	basicInfo << "#chrom\tstart\tend\tname\tlength\tstrand\tsuccess\tuniqHaps\treadTotal\treadTotalUsed\ttotalPairedReads";
	basicInfo << "\tsample";
	uint32_t maxExtraFields = 0;
	for(const auto & p : inputRegions){

		auto bedOut = p.genBedRecordCore();
		if(bedOut.extraFields_.size() > maxExtraFields){
			maxExtraFields = bedOut.extraFields_.size();
		}
	}
	for(uint32_t t = 0; t < maxExtraFields; ++t){
		basicInfo << "\textraField"<<t;
	}
	basicInfo << "\n";

	std::map<uint32_t, uint32_t> coiCounts;

	for (const auto & reg : inputRegions) {
		const auto & regInfos = njh::mapAt(regInfosByUID, reg.uid_);
		++coiCounts[regInfos.front()->uniqHaps_];
		for(auto & regInfo : regInfos){
			auto bedOut = regInfo->region_.genBedRecordCore();
			basicInfo << bedOut.toDelimStr();
			basicInfo << "\t" << njh::boolToStr(regInfo->infoCalled_)
								<< "\t" << regInfo->uniqHaps_
								<< "\t" << regInfo->totalReads_
								<< "\t" << regInfo->totalFinalReads_
								<< "\t" << regInfo->totalPairedReads_
								<< "\t" << sample;
			for(const auto & extra : bedOut.extraFields_){
				basicInfo << "\t" << extra;
			}
			basicInfo << std::endl;
		}
	}

	OutputStream coiOut(njh::files::make_path(finalDirectory, "coiCounts.tab.txt"));
	coiOut << "coi\tcount" << std::endl;
	for(const auto & count : coiCounts){
		coiOut << count.first << "\t" << count.second << std::endl;
	}

	OutputStream exceptionsOut(njh::files::make_path(finalDirectory, "exceptionsMessages.tab.txt"));
	exceptionsOut << "regionUID\tmessage" << std::endl;
	for(const auto & exp : exceptions){
		exceptionsOut << exp.first << "\t" << exp.second << std::endl;
	}



	return 0;
}





}  // namespace njhseq

