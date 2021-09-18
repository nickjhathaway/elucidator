/*
 * programWrappersAssembleOnPathWeaverRunner.cpp
 *
 *  Created on: Sep 8, 2021
 *      Author: nick
 */

#include "elucidator/simulation.h"
#include "elucidator/BamToolsUtils.h"
#include "elucidator/objects/seqContainers.h"
#include "elucidator/objects/BioDataObject.h"
#include "elucidator/seqToolsUtils/seqToolsUtils.hpp"


#include <njhseq/objects/seqObjects/seqKmers.h>
#include <njhseq/system.h>
#include <njhseq/BamToolsUtils.h>
#include <PathWeaver/objects/bam/RegionInvestigatorInBam.hpp>

#include "programWrappersAssembleOnPathWeaverRunner.hpp"

namespace njhseq {






programWrappersAssembleOnPathWeaverRunner::programWrappersAssembleOnPathWeaverRunner()
    : njh::progutils::ProgramRunner(
          {
					 addFunc("runSpadesOnPathWeaverRegions", runSpadesOnPathWeaverRegions, false),
					 addFunc("runMegahitOnPathWeaverRegions", runMegahitOnPathWeaverRegions, false),
					 addFunc("runPRICEOnPathWeaverRegions", runPRICEOnPathWeaverRegions, false),



					 addFunc("runSavageOnPathWeaverRegions", runSavageOnPathWeaverRegions, false),
					 addFunc("runVelvetOptimizerAndMetaVelvetOnPathWeaverRegions", runVelvetOptimizerAndMetaVelvetOnPathWeaverRegions, false),

           },//
          "programWrappersAssembleOnPathWeaverRunner") {}





BamExtractor::ExtractCounts rawWriteExtractReadsFromBamOnlyMapped(const bfs::path & bamFnp,
		const OutOptions & outOpts,
		bool keepRefStrainDir = false) {
	BamExtractor::ExtractCounts ret;
	BamTools::BamReader bReader;
	BamTools::BamAlignment bAln;
	bReader.Open(bamFnp.string());
	checkBamOpenThrow(bReader, bamFnp);
	//loadBamIndexThrow(bReader);
	auto refs = bReader.GetReferenceData();
	BamAlnsCache alnCache;
	auto refData = bReader.GetReferenceData();
	std::unordered_map<std::string, uint32_t> refNameToId;
	for (auto pos : iter::range(refData.size())) {
		refNameToId[refData[pos].RefName] = pos;
	}
	//pair writer
	SeqIOOptions outOptsPaired(outOpts.outFilename_,
			SeqIOOptions::outFormats::FASTQPAIRED, outOpts);
	SeqOutput pairWriter(outOptsPaired);

	//non paired writer
	SeqIOOptions outOptsSingle(outOpts.outFilename_,
			SeqIOOptions::outFormats::FASTQ, outOpts);
	SeqOutput writer(outOptsSingle);

	while (bReader.GetNextAlignment(bAln)) {
		//skip secondary alignments
		if (!bAln.IsPrimaryAlignment()) {
			continue;
		}
		//won't handled sequences with pairs that have a unmmaped mate, will make them orphan reads
		//should improve upon
		if(!bAln.IsMapped()){
			continue;
		}
		/**@todo consider skipping duplicates */
		if (bAln.IsPaired()) {
			if (!alnCache.has(bAln.Name)) {
				//pair hasn't been added to cache yet so add to cache
				//this only works if mate and first mate have the same name
				alnCache.add(bAln);
				continue;
			} else {
				auto search = alnCache.get(bAln.Name);
				if (bAln.IsFirstMate()) {
					pairWriter.openWrite(
							PairedRead(bamAlnToSeqInfo(bAln, keepRefStrainDir), bamAlnToSeqInfo(*search, keepRefStrainDir),
									false));
				} else {
					pairWriter.openWrite(
							PairedRead(bamAlnToSeqInfo(*search, keepRefStrainDir), bamAlnToSeqInfo(bAln, keepRefStrainDir),
									false));
				}
				++ret.pairedReads_;
				// now that operations have been computed, remove first mate found from cache
				alnCache.remove(search->Name);
				continue;
			}
		} else {
			//unpaired read
			++ret.unpaiedReads_;
			if (!bAln.IsMapped()) {
				// do non-mapping operation
				writer.openWrite(bamAlnToSeqInfo(bAln, keepRefStrainDir));
			} else {
				//do unpaired read operation
				writer.openWrite(bamAlnToSeqInfo(bAln, keepRefStrainDir));
			}
		}
	}
	//save the orphans;
	if (len(alnCache) > 0) {
		auto names = alnCache.getNames();
		for (const auto & name : names) {
			++ret.orphans_;
			auto search = alnCache.get(name);
			writer.openWrite(bamAlnToSeqInfo(*search, keepRefStrainDir));
			alnCache.remove(name);
		}
	}
	return ret;
}




struct DefaultAssembleNameInfo{

//	DefaultAssembleNameInfo(const std::string & fullname):fullname_(fullname){
//		setInfoFromName();
//	}

	DefaultAssembleNameInfo(const std::string & fullname, bool megahit = false):fullname_(fullname){
		if(megahit){
			setInfoFromNameMegahit();
		}else{
			setInfoFromName();
		}
	}

	std::string fullname_;

	void setInfoFromName(){
		std::smatch match;
		std::regex pat{R"((NODE_\d+)_length_(\d+)_cov_([0-9.]+).*)"};
		if(!std::regex_match(fullname_, match, pat)){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error in processing " << fullname_ << " for basic assembly info" << "\n";
			throw std::runtime_error{ss.str()};
		}
		nodeName_ = match[1];
		len_ =  njh::StrToNumConverter::stoToNum<uint32_t>(match[2]);
		coverage_ =  njh::StrToNumConverter::stoToNum<double>(match[3]);
	}
	//k99_0 flag=1 multi=94.5593 len=1027
	void setInfoFromNameMegahit(){
		std::smatch match;
		std::regex pat{R"((k[0-9]+_\d+) flag=\d+ multi=([0-9.]+) len=(\d+).*)"};
		if(!std::regex_match(fullname_, match, pat)){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error in processing " << fullname_ << " for basic assembly info" << "\n";
			throw std::runtime_error{ss.str()};
		}
		nodeName_ = match[1];

		coverage_ =  njh::StrToNumConverter::stoToNum<double>(match[2]);

		len_ =  njh::StrToNumConverter::stoToNum<uint32_t>(match[3]);
	}
	std::string nodeName_;
	uint32_t len_;
	double coverage_;

};


int programWrappersAssembleOnPathWeaverRunner::runSpadesOnPathWeaverRegions(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bedFile = "";
	bfs::path pwOutputDir = "";
	std::string sample = "";

	uint32_t spadesNumThreads = 1;
	std::string extraSpadesOptions = "";
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

	VecStr regionNames;
	for(const auto & reg : inputRegions){
		regionNames.emplace_back(reg.uid_);
	}
	njh::sort(regionNames);
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
		std::string regionUid = "";
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


				SeqInput contigsReader(SeqIOOptions::genFastaIn(contigsFnp));
				auto contigsSeqs = contigsReader.readAllReads<seqInfo>();
				std::vector<std::shared_ptr<seqWithKmerInfo>> contigsKmerReads;
				for (const auto & seq : contigsSeqs) {
					contigsKmerReads.emplace_back(std::make_shared<seqWithKmerInfo>(seq));
				}
				allSetKmers(contigsKmerReads, reOrientingKmerLength, true);

				SeqInput refReader(SeqIOOptions::genFastaIn(refFnp));
				auto refSeqs = refReader.readAllReads<seqInfo>();
				std::vector<std::shared_ptr<seqWithKmerInfo>> refKmerReads;
				for (const auto & seq : refSeqs) {
					refKmerReads.emplace_back(std::make_shared<seqWithKmerInfo>(seq));
				}
				allSetKmers(refKmerReads, reOrientingKmerLength, true);

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
				std::vector<kmerInfo> inputSeqsKmerInfos;
				for(const auto & input : refSeqs){
					inputSeqsKmerInfos.emplace_back(input.seq_, 7, false);
				}
				readVecTrimmer::trimSeqToRefByGlobalAln(contigsKmerReads, refSeqs, inputSeqsKmerInfos, alignerObj	);
				//alignerObj.processAlnInfoOutputNoCheck(njh::files::make_path(resultsDirectory, "trimAlnCache").string(), setUp.pars_.verbose_);

				std::vector<std::shared_ptr<seqWithKmerInfo>> finalSeqs;
				for(auto & seq : contigsKmerReads){
					bool found = false;
					for(const auto & finalSeq : finalSeqs){
						if(finalSeq->seqBase_.seq_ == seq->seqBase_.seq_){
							found = true;
							break;
						}
					}
					if(!found){
						finalSeqs.emplace_back(seq);
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
					seqMeta.resetMetaInName(seq->seqBase_.name_);
					seq->seqBase_.cnt_ = (assembleInfo.coverage_/totalCoverage) * (readCounts.pairedReads_ + readCounts.unpaiedReads_ + readCounts.orphans_);
					seq->seqBase_.name_ += njh::pasteAsStr("_t", seq->seqBase_.cnt_);
				}

				OutOptions trimmedContigInfoOpts(njh::files::make_path(spadesFullOutputDir, "trimmed_reOriented_contigs_outputInfo.tab.txt"));
				OutputStream trimmedContigInfoOut(trimmedContigInfoOpts);
				trimmedContigInfoOut << "name\tlength\tcoverage" << std::endl;
				auto trimmkedReOrientedContigsFnp = njh::files::make_path(spadesFullOutputDir, "trimmed_reOriented_contigs.fasta");
				SeqOutput outputWriter(SeqIOOptions::genFastaOut(trimmkedReOrientedContigsFnp));
				auto trimmkedReOrientedContigsFnp_belowCutOff = njh::files::make_path(spadesFullOutputDir, "trimmed_reOriented_contigs_belowCutOff.fasta");
				SeqOutput belowCutOffOutputWriter(SeqIOOptions::genFastaOut(trimmkedReOrientedContigsFnp_belowCutOff));

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
		const auto & regInfo = njh::mapAt(regInfosByUID, reg.uid_);
		++coiCounts[regInfo.front()->uniqHaps_];
		for(auto & reg : regInfo){
			auto bedOut = reg->region_.genBedRecordCore();
			basicInfo << bedOut.toDelimStr();
			basicInfo << "\t" << njh::boolToStr(reg->infoCalled_)
								<< "\t" << reg->uniqHaps_
								<< "\t" << reg->totalReads_
								<< "\t" << reg->totalFinalReads_
								<< "\t" << reg->totalPairedReads_
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
	std::string sample = "";

	uint32_t megahitNumThreads = 1;
	std::string extraMegahitOptions = "";
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

	VecStr regionNames;
	for(const auto & reg : inputRegions){
		regionNames.emplace_back(reg.uid_);
	}
	njh::sort(regionNames);
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
		std::string regionUid = "";
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

				SeqInput contigsReader(SeqIOOptions::genFastaIn(contigsFnp));
				auto contigsSeqs = contigsReader.readAllReads<seqInfo>();
				std::vector<std::shared_ptr<seqWithKmerInfo>> contigsKmerReads;
				for (const auto & seq : contigsSeqs) {
					contigsKmerReads.emplace_back(std::make_shared<seqWithKmerInfo>(seq));
				}
				allSetKmers(contigsKmerReads, reOrientingKmerLength, true);

				SeqInput refReader(SeqIOOptions::genFastaIn(refFnp));
				auto refSeqs = refReader.readAllReads<seqInfo>();
				std::vector<std::shared_ptr<seqWithKmerInfo>> refKmerReads;
				for (const auto & seq : refSeqs) {
					refKmerReads.emplace_back(std::make_shared<seqWithKmerInfo>(seq));
				}
				allSetKmers(refKmerReads, reOrientingKmerLength, true);

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
				std::vector<kmerInfo> inputSeqsKmerInfos;
				for(const auto & input : refSeqs){
					inputSeqsKmerInfos.emplace_back(input.seq_, 7, false);
				}
				readVecTrimmer::trimSeqToRefByGlobalAln(contigsKmerReads, refSeqs, inputSeqsKmerInfos, alignerObj	);
				//alignerObj.processAlnInfoOutputNoCheck(njh::files::make_path(resultsDirectory, "trimAlnCache").string(), setUp.pars_.verbose_);

				std::vector<std::shared_ptr<seqWithKmerInfo>> finalSeqs;
				for(auto & seq : contigsKmerReads){
					bool found = false;
					for(const auto & finalSeq : finalSeqs){
						if(finalSeq->seqBase_.seq_ == seq->seqBase_.seq_){
							found = true;
							break;
						}
					}
					if(!found){
						finalSeqs.emplace_back(seq);
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
					seqMeta.resetMetaInName(seq->seqBase_.name_);
					seq->seqBase_.cnt_ = (assembleInfo.coverage_/totalCoverage) * (readCounts.pairedReads_ + readCounts.unpaiedReads_ + readCounts.orphans_);
					seq->seqBase_.name_ += njh::pasteAsStr("_t", seq->seqBase_.cnt_);
				}

				OutOptions trimmedContigInfoOpts(njh::files::make_path(megahitFullOutputDir, "trimmed_reOriented_contigs_outputInfo.tab.txt"));
				OutputStream trimmedContigInfoOut(trimmedContigInfoOpts);
				trimmedContigInfoOut << "name\tlength\tcoverage" << std::endl;
				auto trimmkedReOrientedContigsFnp = njh::files::make_path(megahitFullOutputDir, "trimmed_reOriented_contigs.fasta");
				SeqOutput outputWriter(SeqIOOptions::genFastaOut(trimmkedReOrientedContigsFnp));
				auto trimmkedReOrientedContigsFnp_belowCutOff = njh::files::make_path(megahitFullOutputDir, "trimmed_reOriented_contigs_belowCutOff.fasta");
				SeqOutput belowCutOffOutputWriter(SeqIOOptions::genFastaOut(trimmkedReOrientedContigsFnp_belowCutOff));

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
		const auto & regInfo = njh::mapAt(regInfosByUID, reg.uid_);
		++coiCounts[regInfo.front()->uniqHaps_];
		for(auto & reg : regInfo){
			auto bedOut = reg->region_.genBedRecordCore();
			basicInfo << bedOut.toDelimStr();
			basicInfo << "\t" << njh::boolToStr(reg->infoCalled_)
								<< "\t" << reg->uniqHaps_
								<< "\t" << reg->totalReads_
								<< "\t" << reg->totalFinalReads_
								<< "\t" << reg->totalPairedReads_
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
	std::string sample = "";

	bfs::path haploconductPath = "/home/hathawan/sourceCodes/savage/HaploConduct/haploconduct";


	uint32_t savageNumThreads = 1;
	std::string extraSavageOptions = "";
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

	VecStr regionNames;
	for(const auto & reg : inputRegions){
		regionNames.emplace_back(reg.uid_);
	}
	njh::sort(regionNames);
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
		std::string regionUid = "";
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

				SeqInput contigsReader(SeqIOOptions::genFastaIn(contigsFnp));
				auto contigsSeqs = contigsReader.readAllReads<seqInfo>();
				std::vector<std::shared_ptr<seqWithKmerInfo>> contigsKmerReads;
				for (const auto & seq : contigsSeqs) {
					contigsKmerReads.emplace_back(std::make_shared<seqWithKmerInfo>(seq));
				}
				allSetKmers(contigsKmerReads, reOrientingKmerLength, true);

				SeqInput refReader(SeqIOOptions::genFastaIn(refFnp));
				auto refSeqs = refReader.readAllReads<seqInfo>();
				std::vector<std::shared_ptr<seqWithKmerInfo>> refKmerReads;
				for (const auto & seq : refSeqs) {
					refKmerReads.emplace_back(std::make_shared<seqWithKmerInfo>(seq));
				}
				allSetKmers(refKmerReads, reOrientingKmerLength, true);

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
				std::vector<kmerInfo> inputSeqsKmerInfos;
				for(const auto & input : refSeqs){
					inputSeqsKmerInfos.emplace_back(input.seq_, 7, false);
				}
				readVecTrimmer::trimSeqToRefByGlobalAln(contigsKmerReads, refSeqs, inputSeqsKmerInfos, alignerObj	);
				//alignerObj.processAlnInfoOutputNoCheck(njh::files::make_path(resultsDirectory, "trimAlnCache").string(), setUp.pars_.verbose_);

				std::vector<std::shared_ptr<seqWithKmerInfo>> finalSeqs;
				for(auto & seq : contigsKmerReads){
					bool found = false;
					for(const auto & finalSeq : finalSeqs){
						if(finalSeq->seqBase_.seq_ == seq->seqBase_.seq_){
							found = true;
							break;
						}
					}
					if(!found){
						finalSeqs.emplace_back(seq);
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
					seqMeta.resetMetaInName(seq->seqBase_.name_);
					seq->seqBase_.cnt_ = (1/totalCoverage) * (readCounts.pairedReads_ + readCounts.unpaiedReads_ + readCounts.orphans_);
					seq->seqBase_.name_ += njh::pasteAsStr("_t", seq->seqBase_.cnt_);
				}

				OutOptions trimmedContigInfoOpts(njh::files::make_path(savageFullOutputDir, "trimmed_reOriented_contigs_outputInfo.tab.txt"));
				OutputStream trimmedContigInfoOut(trimmedContigInfoOpts);
				trimmedContigInfoOut << "name\tlength\tcoverage" << std::endl;
				auto trimmkedReOrientedContigsFnp = njh::files::make_path(savageFullOutputDir, "trimmed_reOriented_contigs.fasta");
				SeqOutput outputWriter(SeqIOOptions::genFastaOut(trimmkedReOrientedContigsFnp));
				auto trimmkedReOrientedContigsFnp_belowCutOff = njh::files::make_path(savageFullOutputDir, "trimmed_reOriented_contigs_belowCutOff.fasta");
				SeqOutput belowCutOffOutputWriter(SeqIOOptions::genFastaOut(trimmkedReOrientedContigsFnp_belowCutOff));

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
		const auto & regInfo = njh::mapAt(regInfosByUID, reg.uid_);
		++coiCounts[regInfo.front()->uniqHaps_];
		for(auto & reg : regInfo){
			auto bedOut = reg->region_.genBedRecordCore();
			basicInfo << bedOut.toDelimStr();
			basicInfo << "\t" << njh::boolToStr(reg->infoCalled_)
								<< "\t" << reg->uniqHaps_
								<< "\t" << reg->totalReads_
								<< "\t" << reg->totalFinalReads_
								<< "\t" << reg->totalPairedReads_
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
	std::string sample = "";

	uint32_t velvetStartKmer = 31;
	uint32_t velvetEndKmer = 71;
	uint32_t velvetKmerStep = 10;
	std::string optFuncKmer = "n50*Lcon/tbp+log(Lbp)";
	std::string optFuncCov = "n50*Lcon/tbp+log(Lbp)";

	uint32_t velvetNumOfThreads = 1;
	VecStr optimizerFuncsAvail{"LNbp","Lbp","Lcon","max","n50","ncon","tbp"};

	bool overWriteDir = false;
	bool breakUpAmbigousContigs = true;
	double coverageCutOff = 2.0;

	std::string extraVelvetOptimiserOptions = "";
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

	VecStr regionNames;
	for(const auto & reg : inputRegions){
		regionNames.emplace_back(reg.uid_);
	}
	njh::sort(regionNames);
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
		std::string regionUid = "";
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
											<< " -o '-exp_cov auto -cov_cutoff 300 -read_trkg yes' ";
				std::string vOptCmdPreCovCutOff = vOptCmdStream.str();
				vOptCmdStream
											<< " -m " << coverageCutOff
											<< " > VelvetOptimzerRunLog_" << njh::getCurrentDate() << ".txt 2>&1";
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
						struct PatPosSize{
							PatPosSize(const std::string & pat, size_t pos): pat_(pat), pos_(pos){

							}
							std::string pat_;
							size_t pos_;

							size_t end(){
								return pos_ + pat_.size();
							}
						};

						for (const auto & seq : contigsSeqs) {
							std::sregex_iterator iter(seq.seq_.begin(), seq.seq_.end(), pat);
							std::sregex_iterator end;
							std::vector<PatPosSize> pats;
							while (iter != end) {
								pats.emplace_back((*iter)[0], iter->position());
								++iter;
							}

							if(!pats.empty()){
								if(0 != pats.front().pos_ ){
									size_t start = 0;
									size_t end = pats.front().pos_;
									auto subSeq = seq.getSubRead(start, end - start);
									if(mark){
										subSeq.name_.append(njh::pasteAsStr("-s", start, "-e", end));
									}
									contigsReader.write(subSeq);
								}
								if(pats.size() > 1){
									for(const auto patPos : iter::range(pats.size() - 1)){
										const auto & p = pats[patPos];
										size_t start = p.pos_ + p.pat_.size();
										size_t end = pats[patPos + 1].pos_;
										auto subSeq = seq.getSubRead(start, end - start);
										if(mark){
											subSeq.name_.append(njh::pasteAsStr("-s", start, "-e", end));
										}
										contigsReader.write(subSeq);
									}
								}
								if(seq.seq_.size() != pats.back().end()){
									size_t start = pats.back().end();
									size_t end = seq.seq_.size();
									auto subSeq = seq.getSubRead(start, end - start);
									if(mark){
										subSeq.name_.append(njh::pasteAsStr("-s", start, "-e", end));
									}
									contigsReader.write(subSeq);
								}
							}else{
								contigsReader.write(seq);
							}
						}
					}


					auto contigsOpts = SeqIOOptions::genFastaIn(contigsFnp);
					SeqInput contigsReader(contigsOpts);
					auto contigsSeqs = contigsReader.readAllReads<seqInfo>();

					std::vector<std::shared_ptr<seqWithKmerInfo>> contigsKmerReads;
					for (const auto & seq : contigsSeqs) {
						contigsKmerReads.emplace_back(std::make_shared<seqWithKmerInfo>(seq));
					}

					allSetKmers(contigsKmerReads, reOrientingKmerLength, true);

					SeqInput refReader(SeqIOOptions::genFastaIn( refFnp ) );
					auto refSeqs = refReader.readAllReads<seqInfo>();
					std::vector<std::shared_ptr<seqWithKmerInfo>> refKmerReads;
					for (const auto & seq : refSeqs) {
						refKmerReads.emplace_back(std::make_shared<seqWithKmerInfo>(seq));
					}
					allSetKmers(refKmerReads, reOrientingKmerLength, true);

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
					std::vector<kmerInfo> inputSeqsKmerInfos;
					for(const auto & input : refSeqs){
						inputSeqsKmerInfos.emplace_back(input.seq_, 7, false);
					}
					readVecTrimmer::trimSeqToRefByGlobalAln(contigsKmerReads, refSeqs, inputSeqsKmerInfos, alignerObj	);
					//alignerObj.processAlnInfoOutputNoCheck(njh::files::make_path(resultsDirectory, "trimAlnCache").string(), setUp.pars_.verbose_);
					std::vector<std::shared_ptr<seqWithKmerInfo>> finalSeqs;
					for(auto & seq : contigsKmerReads){
						bool found = false;
						for(const auto & finalSeq : finalSeqs){
							if(finalSeq->seqBase_.seq_ == seq->seqBase_.seq_){
								found = true;
								break;
							}
						}
						if(!found){
							finalSeqs.emplace_back(seq);
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
						seqMeta.resetMetaInName(seq->seqBase_.name_);
						seq->seqBase_.cnt_ = (assembleInfo.coverage_/totalCoverage) * (readCounts.pairedReads_ + readCounts.unpaiedReads_ + readCounts.orphans_);
						seq->seqBase_.name_ += njh::pasteAsStr("_t", seq->seqBase_.cnt_);
					}

					OutOptions trimmedContigInfoOpts(njh::files::make_path(vOptFullOutputDir, "trimmed_reOriented_contigs_outputInfo.tab.txt"));
					OutputStream trimmedContigInfoOut(trimmedContigInfoOpts);
					trimmedContigInfoOut << "name\tlength\tcoverage" << std::endl;
					auto trimmkedReOrientedContigsFnp = njh::files::make_path(vOptFullOutputDir, "trimmed_reOriented_contigs.fasta");
					SeqOutput outputWriter(SeqIOOptions::genFastaOut(trimmkedReOrientedContigsFnp));
					auto trimmkedReOrientedContigsFnp_belowCutOff = njh::files::make_path(vOptFullOutputDir, "trimmed_reOriented_contigs_belowCutOff.fasta");
					SeqOutput belowCutOffOutputWriter(SeqIOOptions::genFastaOut(trimmkedReOrientedContigsFnp_belowCutOff));

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
				std::string velvethArgs = "";
				std::string velvetgArgs = "";
				while(njh::files::crossPlatGetline(logFileIn, line)){
					if(njh::beginsWith(line, "Final optimised assembly details")){
						encounteredFinalInfo = true;
					}
					if(encounteredFinalInfo && std::string::npos != line.find(":")){
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

				if("" == velvethArgs){
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error " << "couldn't determine velveth args from " << logfiles.front()<< "\n";
					throw std::runtime_error{ss.str()};
				}

				if("" == velvetgArgs){
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error " << "couldn't determine velvetg args from " << logfiles.front()<< "\n";
					throw std::runtime_error{ss.str()};
				}
			//	std::string velvetOutputDirStr = velvetgArgs;
			//	trimAtFirstWhitespace(velvetOutputDirStr);
				std::string velvetOutputDirStr = "optimized_velvet_run";
				bfs::path velvetOutputDir = njh::files::make_path(VelvetOptimiserOutDir, "MetaVelvetOutDir", velvetOutputDirStr);

				std::stringstream velvethCmd;
				velvethCmd << "cd " << vOptFullOutputDir << "/../ && velveth " << velvetOutputDir << " " << velvethArgs.substr(velvethArgs.find(" "));

				std::stringstream velvetgCmd;
				velvetgCmd << "cd " << vOptFullOutputDir << "/../ && velvetg " << velvetOutputDir << " " << velvetgArgs.substr(velvetgArgs.find(" "));

				std::stringstream metavelvetgCmd;
				metavelvetgCmd << "cd " << vOptFullOutputDir << "/../ && meta-velvetg " << velvetOutputDir;

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
						struct PatPosSize{
							PatPosSize(const std::string & pat, size_t pos): pat_(pat), pos_(pos){

							}
							std::string pat_;
							size_t pos_;

							size_t end(){
								return pos_ + pat_.size();
							}
						};

						for (const auto & seq : contigsSeqs) {
							std::sregex_iterator iter(seq.seq_.begin(), seq.seq_.end(), pat);
							std::sregex_iterator end;
							std::vector<PatPosSize> pats;
							while (iter != end) {
								pats.emplace_back((*iter)[0], iter->position());
								++iter;
							}

							if(!pats.empty()){
								if(0 != pats.front().pos_ ){
									size_t start = 0;
									size_t end = pats.front().pos_;
									auto subSeq = seq.getSubRead(start, end - start);
									if(mark){
										subSeq.name_.append(njh::pasteAsStr("-s", start, "-e", end));
									}
									contigsReader.write(subSeq);
								}
								if(pats.size() > 1){
									for(const auto patPos : iter::range(pats.size() - 1)){
										const auto & p = pats[patPos];
										size_t start = p.pos_ + p.pat_.size();
										size_t end = pats[patPos + 1].pos_;
										auto subSeq = seq.getSubRead(start, end - start);
										if(mark){
											subSeq.name_.append(njh::pasteAsStr("-s", start, "-e", end));
										}
										contigsReader.write(subSeq);
									}
								}
								if(seq.seq_.size() != pats.back().end()){
									size_t start = pats.back().end();
									size_t end = seq.seq_.size();
									auto subSeq = seq.getSubRead(start, end - start);
									if(mark){
										subSeq.name_.append(njh::pasteAsStr("-s", start, "-e", end));
									}
									contigsReader.write(subSeq);
								}
							}else{
								contigsReader.write(seq);
							}
						}
					}


					auto contigsOpts = SeqIOOptions::genFastaIn(metaVelvetContigsfile);
					SeqInput contigsReader(contigsOpts);
					auto contigsSeqs = contigsReader.readAllReads<seqInfo>();

					std::vector<std::shared_ptr<seqWithKmerInfo>> contigsKmerReads;
					for (const auto & seq : contigsSeqs) {
						contigsKmerReads.emplace_back(std::make_shared<seqWithKmerInfo>(seq));
					}

					allSetKmers(contigsKmerReads, reOrientingKmerLength, true);

					SeqInput refReader(SeqIOOptions::genFastaIn( refFnp ) );
					auto refSeqs = refReader.readAllReads<seqInfo>();
					std::vector<std::shared_ptr<seqWithKmerInfo>> refKmerReads;
					for (const auto & seq : refSeqs) {
						refKmerReads.emplace_back(std::make_shared<seqWithKmerInfo>(seq));
					}
					allSetKmers(refKmerReads, reOrientingKmerLength, true);

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
					std::vector<kmerInfo> inputSeqsKmerInfos;
					for(const auto & input : refSeqs){
						inputSeqsKmerInfos.emplace_back(input.seq_, 7, false);
					}
					readVecTrimmer::trimSeqToRefByGlobalAln(contigsKmerReads, refSeqs, inputSeqsKmerInfos, alignerObj	);
					//alignerObj.processAlnInfoOutputNoCheck(njh::files::make_path(resultsDirectory, "trimAlnCache").string(), setUp.pars_.verbose_);
					std::vector<std::shared_ptr<seqWithKmerInfo>> finalSeqs;
					for(auto & seq : contigsKmerReads){
						bool found = false;
						for(const auto & finalSeq : finalSeqs){
							if(finalSeq->seqBase_.seq_ == seq->seqBase_.seq_){
								found = true;
								break;
							}
						}
						if(!found){
							finalSeqs.emplace_back(seq);
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
						seqMeta.resetMetaInName(seq->seqBase_.name_);
						seq->seqBase_.cnt_ = (assembleInfo.coverage_/totalCoverage) * (readCounts.pairedReads_ + readCounts.unpaiedReads_ + readCounts.orphans_);
						seq->seqBase_.name_ += njh::pasteAsStr("_t", seq->seqBase_.cnt_);
					}

					OutOptions trimmedContigInfoOpts(njh::files::make_path(metaVelvetDir, "trimmed_reOriented_meta-velvetg.contigs_outputInfo.tab.txt"));
					OutputStream trimmedContigInfoOut(trimmedContigInfoOpts);
					trimmedContigInfoOut << "name\tlength\tcoverage" << std::endl;
					auto trimmkedReOrientedContigsFnp = njh::files::make_path(metaVelvetDir, "trimmed_reOriented_meta-velvetg.contigs.fasta");
					SeqOutput outputWriter(SeqIOOptions::genFastaOut(trimmkedReOrientedContigsFnp));
					auto trimmkedReOrientedContigsFnp_belowCutOff = njh::files::make_path(metaVelvetDir, "trimmed_reOriented_meta-velvetg.contigs_belowCutOff.fasta");
					SeqOutput belowCutOffOutputWriter(SeqIOOptions::genFastaOut(trimmkedReOrientedContigsFnp_belowCutOff));

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
			const auto & regInfo = njh::mapAt(regInfosByUID, reg.uid_);
			++coiCounts[regInfo.front()->uniqHaps_];
			for(auto & reg : regInfo){
				auto bedOut = reg->region_.genBedRecordCore();
				basicInfo << bedOut.toDelimStr();
				basicInfo << "\t" << njh::boolToStr(reg->infoCalled_)
									<< "\t" << reg->uniqHaps_
									<< "\t" << reg->totalReads_
									<< "\t" << reg->totalFinalReads_
									<< "\t" << reg->totalPairedReads_
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
			const auto & regInfo = njh::mapAt(metav_regInfosByUID, reg.uid_);
			++coiCounts[regInfo.front()->uniqHaps_];
			for(auto & reg : regInfo){
				auto bedOut = reg->region_.genBedRecordCore();
				basicInfo << bedOut.toDelimStr();
				basicInfo << "\t" << njh::boolToStr(reg->infoCalled_)
									<< "\t" << reg->uniqHaps_
									<< "\t" << reg->totalReads_
									<< "\t" << reg->totalFinalReads_
									<< "\t" << reg->totalPairedReads_
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
	std::string sample = "";

	uint32_t PRICENumThreads = 1;
	std::string extraPRICEOptions = "";
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

	VecStr regionNames;
	for(const auto & reg : inputRegions){
		regionNames.emplace_back(reg.uid_);
	}
	njh::sort(regionNames);
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
		std::string regionUid = "";
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

				SeqInput contigsReader(SeqIOOptions::genFastaIn(contigsFnp));
				auto contigsSeqs = contigsReader.readAllReads<seqInfo>();
				std::vector<std::shared_ptr<seqWithKmerInfo>> contigsKmerReads;
				for (const auto & seq : contigsSeqs) {
					contigsKmerReads.emplace_back(std::make_shared<seqWithKmerInfo>(seq));
				}
				allSetKmers(contigsKmerReads, reOrientingKmerLength, true);

				SeqInput refReader(SeqIOOptions::genFastaIn(refFnp));
				auto refSeqs = refReader.readAllReads<seqInfo>();
				std::vector<std::shared_ptr<seqWithKmerInfo>> refKmerReads;
				for (const auto & seq : refSeqs) {
					refKmerReads.emplace_back(std::make_shared<seqWithKmerInfo>(seq));
				}
				allSetKmers(refKmerReads, reOrientingKmerLength, true);

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
				contigInfoOut << "name\tlength\tcoverage" << std::endl;

				for(const auto & contigsKmerRead : contigsKmerReads){
					auto assembleInfo = DefaultAssembleNameInfo(contigsKmerRead->seqBase_.name_, true);
					contigInfoOut << contigsKmerRead->seqBase_.name_
							<< "\t" << len(contigsKmerRead->seqBase_)
							<< "\t" << assembleInfo.coverage_ << std::endl;
				}
				auto reOrientedContigsFnp = njh::files::make_path(PRICEFullOutputDir, "reOriented_contigs.fasta");

				SeqOutput::write(contigsKmerReads, SeqIOOptions::genFastaOut(reOrientedContigsFnp));

				uint64_t maxLen = 0;
				readVec::getMaxLength(refSeqs, maxLen);
				readVec::getMaxLength(contigsKmerReads, maxLen);
				aligner alignerObj(maxLen, gapScoringParameters(5,1,0,0,0,0), substituteMatrix(2,-2), false);
				//alignerObj.processAlnInfoInputNoCheck(njh::files::make_path(resultsDirectory, "trimAlnCache").string(), setUp.pars_.verbose_);
				std::vector<kmerInfo> inputSeqsKmerInfos;
				for(const auto & input : refSeqs){
					inputSeqsKmerInfos.emplace_back(input.seq_, 7, false);
				}
				readVecTrimmer::trimSeqToRefByGlobalAln(contigsKmerReads, refSeqs, inputSeqsKmerInfos, alignerObj	);
				//alignerObj.processAlnInfoOutputNoCheck(njh::files::make_path(resultsDirectory, "trimAlnCache").string(), setUp.pars_.verbose_);

				std::vector<std::shared_ptr<seqWithKmerInfo>> finalSeqs;
				for(auto & seq : contigsKmerReads){
					bool found = false;
					for(const auto & finalSeq : finalSeqs){
						if(finalSeq->seqBase_.seq_ == seq->seqBase_.seq_){
							found = true;
							break;
						}
					}
					if(!found){
						finalSeqs.emplace_back(seq);
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
					seqMeta.resetMetaInName(seq->seqBase_.name_);
					seq->seqBase_.cnt_ = (assembleInfo.coverage_/totalCoverage) * (readCounts.pairedReads_ + readCounts.unpaiedReads_ + readCounts.orphans_);
					seq->seqBase_.name_ += njh::pasteAsStr("_t", seq->seqBase_.cnt_);
				}

				OutOptions trimmedContigInfoOpts(njh::files::make_path(PRICEFullOutputDir, "trimmed_reOriented_contigs_outputInfo.tab.txt"));
				OutputStream trimmedContigInfoOut(trimmedContigInfoOpts);
				trimmedContigInfoOut << "name\tlength\tcoverage" << std::endl;
				auto trimmkedReOrientedContigsFnp = njh::files::make_path(PRICEFullOutputDir, "trimmed_reOriented_contigs.fasta");
				SeqOutput outputWriter(SeqIOOptions::genFastaOut(trimmkedReOrientedContigsFnp));
				auto trimmkedReOrientedContigsFnp_belowCutOff = njh::files::make_path(PRICEFullOutputDir, "trimmed_reOriented_contigs_belowCutOff.fasta");
				SeqOutput belowCutOffOutputWriter(SeqIOOptions::genFastaOut(trimmkedReOrientedContigsFnp_belowCutOff));

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
		const auto & regInfo = njh::mapAt(regInfosByUID, reg.uid_);
		++coiCounts[regInfo.front()->uniqHaps_];
		for(auto & reg : regInfo){
			auto bedOut = reg->region_.genBedRecordCore();
			basicInfo << bedOut.toDelimStr();
			basicInfo << "\t" << njh::boolToStr(reg->infoCalled_)
								<< "\t" << reg->uniqHaps_
								<< "\t" << reg->totalReads_
								<< "\t" << reg->totalFinalReads_
								<< "\t" << reg->totalPairedReads_
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

