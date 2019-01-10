/*
 * programWrappers_runAdapterRemoval.cpp
 *
 *  Created on: Aug 9, 2017
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

#include "programWrappers.hpp"
#include "elucidator/objects/BioDataObject.h"
#include "elucidator/BamToolsUtils.h"



namespace njhseq {

struct RunAdapterRemovalPars {
	SeqIOOptions inOpts;
	bool gzip = false;
	std::string adapter1 = "";
	std::string adapter2 = "";
	bfs::path adapterList = "";
	uint32_t numThreads = 1;
	bool combineSingles = false;
	bool noTrimmingNsAndQuals = false;
	OutOptions outOpts{bfs::path("")};
	bool force = false;
	bool doCollapse = false;
	uint32_t maxns = 2;
	std::string extraArgs = "";
	bool keepSinglesAfterCombining = false;
	uint32_t maxQual = 52;

	void processArgs(seqSetUp & setUp){
		setUp.setOption(adapterList, "--adapterList", "Adapter List");
		setUp.setOption(force, "--force", "force");
		setUp.setOption(keepSinglesAfterCombining, "--keepSinglesAfterCombining", "Keep Singles After Combining, otherwise they are deleted");
		setUp.setOption(maxns, "--maxns", "max ns");
		setUp.setOption(doCollapse, "--doCollapse", "Merge (collapse) overlapping reads");
		setUp.setOption(noTrimmingNsAndQuals, "--noTrimmingNsAndQuals", "Don't Trim Ns And Quals");
		setUp.setOption(combineSingles, "--combineSingles", "Combine single files all into 1");
		setUp.setOption(gzip, "--gzip", "Make output gzip");
		setUp.setOption(numThreads, "--numThreads", "Number of threads to utilize");
		setUp.setOption(adapter1, "--adapter1", "The forward adapter");
		setUp.setOption(adapter2, "--adapter2", "The forward adapter");
		setUp.setOption(extraArgs, "--extraArgs", "Extra Arguments to AdapterRemoval");
		setUp.setOption(maxQual, "--qualityMax", "quality max for cut adapter, otherwise will throw if a quality is found above this score");

		setUp.processWritingOptions(outOpts);
		setUp.processReadInNames(VecStr{"--fastq1", "--fastq2"});
		inOpts = setUp.pars_.ioOptions_;
	}
	void processArgsSE(seqSetUp & setUp){
		setUp.setOption(adapterList, "--adapterList", "Adapter List");
		setUp.setOption(force, "--force", "force");
		setUp.setOption(maxns, "--maxns", "max ns");
		setUp.setOption(noTrimmingNsAndQuals, "--noTrimmingNsAndQuals", "Don't Trim Ns And Quals");
		setUp.setOption(gzip, "--gzip", "Make output gzip");
		setUp.setOption(numThreads, "--numThreads", "Number of threads to utilize");
		setUp.setOption(adapter1, "--adapter1", "The forward adapter");
		setUp.setOption(extraArgs, "--extraArgs", "Extra Arguments to AdapterRemoval");
		setUp.setOption(maxQual, "--qualityMax", "quality max for cut adapter, otherwise will throw if a quality is found above this score");
		setUp.processWritingOptions(outOpts);
		setUp.processReadInNames(VecStr{"--fastq", "--fastqgz"});
		inOpts = setUp.pars_.ioOptions_;
		if(njh::endsWith(inOpts.firstName_.string(), ".gz")){
			gzip = true;
		}
	}
};


njh::sys::RunOutput RunAdapterRemoval(RunAdapterRemovalPars adapRemovPars){
	njh::sys::requireExternalProgramThrow("AdapterRemoval");
	std::string extention = adapRemovPars.gzip ? ".fastq.gz" : ".fastq";

	njh::files::checkExistenceThrow(adapRemovPars.inOpts.firstName_, __PRETTY_FUNCTION__);
	njh::files::checkExistenceThrow(adapRemovPars.inOpts.secondName_, __PRETTY_FUNCTION__);
	if("" == adapRemovPars.outOpts.outFilename_){
		auto outName = bfs::path(bfs::basename(adapRemovPars.inOpts.firstName_)).filename().string();
		adapRemovPars.outOpts.outFilename_ = "trimmed_" + outName.substr(0, outName.find("_"));
	}
	//check if output exists already

	OutOptions output1(bfs::path(adapRemovPars.outOpts.outFilename_.string() + "_1" + extention));
	OutOptions outputCollapse(bfs::path(adapRemovPars.outOpts.outFilename_.string() + "" + extention));
	output1.transferOverwriteOpts(adapRemovPars.outOpts);
	outputCollapse.transferOverwriteOpts(adapRemovPars.outOpts);
	output1.throwIfOutExistsNoOverWrite(__PRETTY_FUNCTION__);
	outputCollapse.throwIfOutExistsNoOverWrite(__PRETTY_FUNCTION__);


	std::stringstream adapterRemvalCmdStream;
	adapterRemvalCmdStream << "AdapterRemoval --threads " << adapRemovPars.numThreads << " ";
	if("" != adapRemovPars.adapterList){
		njh::files::checkExistenceThrow(adapRemovPars.adapterList, __PRETTY_FUNCTION__);
		adapterRemvalCmdStream << "--adapter-list " << adapRemovPars.adapterList << " ";
	}else{
		if("" != adapRemovPars.adapter1){
			adapterRemvalCmdStream<< "--adapter1 " << adapRemovPars.adapter1 << " ";
		}
		if("" != adapRemovPars.adapter2){
			adapterRemvalCmdStream<< "--adapter2 " << adapRemovPars.adapter2 << " ";
		}
	}
	if(!adapRemovPars.noTrimmingNsAndQuals){
		adapterRemvalCmdStream << "--trimns --trimqualities ";
	}
	adapterRemvalCmdStream << "--maxns " << adapRemovPars.maxns << " ";
	adapterRemvalCmdStream
			<< "--file1 " << adapRemovPars.inOpts.firstName_ << " "
			<< "--file2 " << adapRemovPars.inOpts.secondName_ << " ";
	adapterRemvalCmdStream << "--qualitymax " << adapRemovPars.maxQual << " ";
	if(adapRemovPars.gzip){
		adapterRemvalCmdStream << "--gzip ";
	}
	if(adapRemovPars.doCollapse){
		adapterRemvalCmdStream << "--collapse ";
	}
	adapterRemvalCmdStream
			<< "--output1 " << adapRemovPars.outOpts.outFilename_.string() << "_1" << extention << " "
			<< "--output2 " << adapRemovPars.outOpts.outFilename_.string() << "_2" << extention << " "
			<< "--outputcollapsed " << adapRemovPars.outOpts.outFilename_.string() << "" << extention << " "
			<< "--outputcollapsedtruncated " << adapRemovPars.outOpts.outFilename_.string() << "_trunc" << extention << " "
			<< "--singleton " << adapRemovPars.outOpts.outFilename_.string() << "_mateLost" << extention << " "
			<< "--discarded " << adapRemovPars.outOpts.outFilename_.string() << "_discarded" << extention << " "
			<< "--settings " << adapRemovPars.outOpts.outFilename_.string() << "_AdapterRemoval.log ";

	//add any extra commands
	if("" != adapRemovPars.extraArgs){
		adapterRemvalCmdStream << adapRemovPars.extraArgs << " ";
	}
	adapterRemvalCmdStream << " >> " << adapRemovPars.outOpts.outFilename_.string() << "_AdapterRemoval.runlog 2>&1 ";

	bool needToRun = adapRemovPars.force;
	if(!adapRemovPars.force){
		bfs::path checkFile = adapRemovPars.outOpts.outFilename_.string() + "_AdapterRemoval.log";
		if( !bfs::exists(checkFile) ||
				njh::files::firstFileIsOlder(checkFile, adapRemovPars.inOpts.firstName_) ||
				njh::files::firstFileIsOlder(checkFile, adapRemovPars.inOpts.secondName_)){
			needToRun = true;
		}
	}
	if(needToRun){
		bfs::path runLogFnp = adapRemovPars.outOpts.outFilename_.string() + "_AdapterRemoval.runlog";
		std::string adapterRemoveCmd = adapterRemvalCmdStream.str();
		{
			OutOptions runLogOpts(runLogFnp);
			runLogOpts.overWriteFile_ = true;
			OutputStream runLogOut(runLogOpts);
			runLogOut << "Ran on: " << njh::getCurrentDate() << std::endl;
			runLogOut << "Ran from: " << bfs::current_path() << std::endl;
			runLogOut << adapterRemoveCmd << std::endl;
		}
		auto runOutput = njh::sys::run({adapterRemoveCmd});
		BioCmdsUtils::checkRunOutThrow(runOutput, __PRETTY_FUNCTION__);

		{
			//remove any empty files as AdapterRemoval just creates files and might not put anything in certain ones if certain conditions don't occur
			std::vector<bfs::path> checkFiles;
			bfs::path check_output1 = adapRemovPars.outOpts.outFilename_.string() + "_1" + extention ;
			bfs::path check_output2 = adapRemovPars.outOpts.outFilename_.string() + "_2" + extention ;
			bfs::path check_outputcollapsed = adapRemovPars.outOpts.outFilename_.string() + "" + extention ;
			bfs::path check_outputcollapsedtruncated = adapRemovPars.outOpts.outFilename_.string() + "_trunc" + extention ;
			bfs::path check_singleton = adapRemovPars.outOpts.outFilename_.string() + "_mateLost" + extention ;
			bfs::path check_discarded = adapRemovPars.outOpts.outFilename_.string() + "_discarded" + extention ;
			checkFiles.emplace_back(check_output1);
			checkFiles.emplace_back(check_output2);
			checkFiles.emplace_back(check_outputcollapsed);
			checkFiles.emplace_back(check_outputcollapsedtruncated);
			checkFiles.emplace_back(check_singleton);
			checkFiles.emplace_back(check_discarded);
			for(const auto & fnp : checkFiles){
				if(bfs::exists(fnp) && 0 == bfs::file_size(fnp)){
					bfs::remove(fnp);
				}
			}
		}
		//combine the single outputs into one to make it easier to map latter;
		if(adapRemovPars.combineSingles){
			bfs::path singletonFnp = adapRemovPars.outOpts.outFilename_.string() + "_mateLost" + extention;
			bfs::path collapdedFnp = adapRemovPars.outOpts.outFilename_.string() + "" + extention;
			bfs::path truncCollapdedFnp = adapRemovPars.outOpts.outFilename_.string() + "_trunc" + extention;
			std::vector<bfs::path> collapseSingleFnps;
			auto addIfExists = [&collapseSingleFnps](const bfs::path & fnp){
				if(bfs::exists(fnp)){
					collapseSingleFnps.emplace_back(fnp);
				}
			};
			addIfExists(singletonFnp);
			addIfExists(collapdedFnp);
			addIfExists(truncCollapdedFnp);

			if(!collapseSingleFnps.empty()){
				OutOptions collapsedSinglesOpts(bfs::path(adapRemovPars.outOpts.outFilename_.string() + "_singles" + extention));
				collapsedSinglesOpts.transferOverwriteOpts(adapRemovPars.outOpts);
				concatenateFiles(collapseSingleFnps, collapsedSinglesOpts);
				for(const auto & fnp : collapseSingleFnps){
					bfs::remove(fnp);
				}
			}
		}
		return runOutput;
	}else{
		return njh::sys::RunOutput{};
	}
}

njh::sys::RunOutput RunAdapterRemovalSE(RunAdapterRemovalPars adapRemovPars){
	njh::sys::requireExternalProgramThrow("AdapterRemoval");
	std::string extention = adapRemovPars.gzip ? ".fastq.gz" : ".fastq";

	njh::files::checkExistenceThrow(adapRemovPars.inOpts.firstName_, __PRETTY_FUNCTION__);
	if("" == adapRemovPars.outOpts.outFilename_){
		auto outName = bfs::path(bfs::basename(adapRemovPars.inOpts.firstName_)).filename().string();
		adapRemovPars.outOpts.outFilename_ = "trimmed_" + outName.substr(0, outName.find("_"));
	}
	//check if output exists already

	OutOptions output(bfs::path(adapRemovPars.outOpts.outFilename_.string() + "" + extention));
	output.transferOverwriteOpts(adapRemovPars.outOpts);
	output.throwIfOutExistsNoOverWrite(__PRETTY_FUNCTION__);
//AdapterRemoval --file1 ../SRR066764_1.fastq.gz   --output1 trimmed_SRR066764.fastq --settings trimmed_SRR066764.log --discarded trimmed_SRR066764_discarded.fastq  --trimns --trimqualities --maxns 2  --threads 10

	std::stringstream adapterRemvalCmdStream;
	adapterRemvalCmdStream << "AdapterRemoval --threads " << adapRemovPars.numThreads << " ";
	if("" != adapRemovPars.adapterList){
		njh::files::checkExistenceThrow(adapRemovPars.adapterList, __PRETTY_FUNCTION__);
		adapterRemvalCmdStream << "--adapter-list " << adapRemovPars.adapterList << " ";
	}else{
		if("" != adapRemovPars.adapter1){
			adapterRemvalCmdStream<< "--adapter1 " << adapRemovPars.adapter1 << " ";
		}
	}
	if(!adapRemovPars.noTrimmingNsAndQuals){
		adapterRemvalCmdStream << "--trimns --trimqualities ";
	}
	adapterRemvalCmdStream << "--maxns " << adapRemovPars.maxns << " ";
	adapterRemvalCmdStream
			<< "--file1 " << adapRemovPars.inOpts.firstName_ << " ";
	adapterRemvalCmdStream << "--qualitymax " << adapRemovPars.maxQual << " ";
	if(adapRemovPars.gzip){
		adapterRemvalCmdStream << "--gzip ";
	}
	adapterRemvalCmdStream
			<< "--output1 " << adapRemovPars.outOpts.outFilename_.string() << "" << extention << " "
			<< "--discarded " << adapRemovPars.outOpts.outFilename_.string() << "_discarded" << extention << " "
			<< "--settings " << adapRemovPars.outOpts.outFilename_.string() << "_AdapterRemoval.log ";
	// add extra commands
	if("" != adapRemovPars.extraArgs){
		adapterRemvalCmdStream << adapRemovPars.extraArgs << " ";
	}
	adapterRemvalCmdStream  << ">> " << adapRemovPars.outOpts.outFilename_.string() << "_AdapterRemoval.runlog 2>&1 ";

	bool needToRun = adapRemovPars.force;
	if(!adapRemovPars.force){
		bfs::path checkFile = adapRemovPars.outOpts.outFilename_.string() + "_AdapterRemoval.log";
		if( !bfs::exists(checkFile) ||
				njh::files::firstFileIsOlder(checkFile, adapRemovPars.inOpts.firstName_)){
			needToRun = true;
		}
	}
	if(needToRun){
		bfs::path runLogFnp = adapRemovPars.outOpts.outFilename_.string() + "_AdapterRemoval.runlog";
		std::string adapterRemoveCmd = adapterRemvalCmdStream.str();
		{
			OutOptions runLogOpts(runLogFnp);
			runLogOpts.overWriteFile_ = true;
			OutputStream runLogOut(runLogOpts);
			runLogOut << "Ran on: " << njh::getCurrentDate() << std::endl;
			runLogOut << "Ran from: " << bfs::current_path() << std::endl;
			runLogOut << adapterRemoveCmd << std::endl;
		}
		auto runOutput = njh::sys::run({adapterRemoveCmd});
		BioCmdsUtils::checkRunOutThrow(runOutput, __PRETTY_FUNCTION__);

		{
			//remove any empty files as AdapterRemoval just creates files and might not put anything in certain ones if certain conditions don't occur
			std::vector<bfs::path> checkFiles;
			bfs::path check_output1 = adapRemovPars.outOpts.outFilename_.string() + "_1" + extention ;
			bfs::path check_discarded = adapRemovPars.outOpts.outFilename_.string() + "_discarded" + extention ;
			checkFiles.emplace_back(check_output1);
			checkFiles.emplace_back(check_discarded);
			for(const auto & fnp : checkFiles){
				if(bfs::exists(fnp) && 0 == bfs::file_size(fnp)){
					bfs::remove(fnp);
				}
			}
		}
		return runOutput;
	}else{
		return njh::sys::RunOutput{};
	}
}



int programWrapperRunner::setUpRunAdapterRemoval(const njh::progutils::CmdArgs & inputCommands){
	RunAdapterRemovalPars adapRemovPars;
	bfs::path detectPrimersDir = "";
	seqSetUp setUp(inputCommands);
	adapRemovPars.processArgs(setUp);
	setUp.setOption(detectPrimersDir, "--detectPrimersDir", "Detect Primers Dir");
	setUp.processVerbose();
	setUp.processDebug();
	setUp.finishSetUp(std::cout);
	njh::files::checkExistenceThrow(detectPrimersDir, __PRETTY_FUNCTION__);
	//check if consensus sequences were detected

	std::string defaultAdapter1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG";
	std::string defaultAdapter2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT";
	auto r1Consensus = njh::files::make_path(detectPrimersDir, "r1_consensus.fasta");
	auto r2Consensus = njh::files::make_path(detectPrimersDir, "r2_consensus.fasta");
	VecStr acceptableSeqs{"TGATACGG","ATCGGAAG","ATGTGTAT","AGCAGAAG","GTCTCTTA","CTCGTGGG","GTCGGCAG"};

	std::vector<seqInfo> acceptedR1Cons;

	if(bfs::exists(r1Consensus)){
		auto r1ConOpts = SeqIOOptions::genFastaIn(r1Consensus);
		seqInfo seq;
		SeqInput r1Reader(r1ConOpts);
		r1Reader.openIn();
		while(r1Reader.readNextRead(seq)){
			if(std::string::npos != seq.seq_.find("N")){
				seq.trimBack(seq.seq_.find("N"));
			}
			for(const auto & accepSeq : acceptableSeqs){
				if(std::string::npos != seq.seq_.find(accepSeq)){
					acceptedR1Cons.emplace_back(seq);
				}
			}
		}
	}else{
		acceptedR1Cons.emplace_back("defaultAdapter1", defaultAdapter1);
	}
	if(acceptedR1Cons.empty()){
		acceptedR1Cons.emplace_back("defaultAdapter1", defaultAdapter1);
	}

	std::vector<seqInfo> acceptedR2Cons;
	if(bfs::exists(r2Consensus)){
		auto r2ConOpts = SeqIOOptions::genFastaIn(r2Consensus);
		seqInfo seq;
		SeqInput r2Reader(r2ConOpts);
		r2Reader.openIn();
		while(r2Reader.readNextRead(seq)){
			if(std::string::npos != seq.seq_.find("N")){
				seq.trimBack(seq.seq_.find("N"));
			}
			for(const auto & accepSeq : acceptableSeqs){
				if(std::string::npos != seq.seq_.find(accepSeq)){
					acceptedR2Cons.emplace_back(seq);
				}
			}
		}
	}else{
		acceptedR2Cons.emplace_back("defaultAdapter2", defaultAdapter2);
	}
	if(acceptedR2Cons.empty()){
		acceptedR2Cons.emplace_back("defaultAdapter2", defaultAdapter2);
	}


	if("" == adapRemovPars.outOpts.outFilename_){
		auto outName = bfs::path(bfs::basename(adapRemovPars.inOpts.firstName_)).filename().string();
		adapRemovPars.outOpts.outFilename_ = "trimmed_" + outName.substr(0, outName.find("_"));
	}

	if(1 == acceptedR1Cons.size() && 1 == acceptedR2Cons.size()){
		adapRemovPars.adapter1 = acceptedR1Cons.front().seq_;
		adapRemovPars.adapter2 = acceptedR2Cons.front().seq_;
	}else{
		OutOptions adapterListOpts(bfs::path(adapRemovPars.outOpts.outFilename_.string() + "_adapterList.txt"));
		adapterListOpts.transferOverwriteOpts(adapRemovPars.outOpts);
		OutputStream adapterListOut(adapterListOpts);
		for(const auto & r1Con : acceptedR1Cons){
			for( auto r2Con : acceptedR2Cons){
				r2Con.reverseComplementRead(false, true);
				adapterListOut << r1Con.seq_ << "\t" << r2Con.seq_ << std::endl;;
			}
		}
		adapRemovPars.adapterList = adapterListOpts.outName();
	}
	RunAdapterRemoval(adapRemovPars);
	return 0;
}




int programWrapperRunner::runAdapterRemoval(const njh::progutils::CmdArgs & inputCommands){
	RunAdapterRemovalPars adapRemovPars;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	adapRemovPars.processArgs(setUp);
	setUp.finishSetUp(std::cout);

	RunAdapterRemoval(adapRemovPars);

	return 0;
}

int programWrapperRunner::runAdapterRemovalSE(const njh::progutils::CmdArgs & inputCommands){
	RunAdapterRemovalPars adapRemovPars;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	adapRemovPars.processArgsSE(setUp);
	setUp.finishSetUp(std::cout);

	RunAdapterRemovalSE(adapRemovPars);

	return 0;
}




int programWrapperRunner::runBwaOnAdapterReomvalOutputSE(const njh::progutils::CmdArgs & inputCommands){

	bfs::path outputDir = "";
	bfs::path trimStub = "";
	bool force = false;
	uint32_t numThreads = 1;
	std::string sampName = "";
	bool removeIntermediateFiles = false;
	bfs::path genomeFnp = "";
	std::string extraBwaArgs = "";
	bfs::path outputFnp = "";
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(genomeFnp, "--genomeFnp", "genome file to align to", true);
	setUp.setOption(trimStub, "--trimStub", "AdapterRemoval output stub", true);
	setUp.setOption(sampName, "--sampName", "Sample Name to give to final bam", true);
	setUp.setOption(force, "--force", "force run even if file already exists");
	setUp.setOption(outputFnp, "--outputFnp", "output name, will default to sampName.sorted.bam");
	setUp.setOption(extraBwaArgs, "--extraBwaArgs", "extra bwa arguments");
	setUp.setOption(numThreads, "--numThreads", "Number of threads");
	setUp.setOption(outputDir, "--outputDir", "output directory");
	setUp.setOption(removeIntermediateFiles, "--removeIntermediateFiles", "remove the intermediate bam files and just keep the ");
	setUp.finishSetUp(std::cout);
	BioCmdsUtils bioRunner(setUp.pars_.verbose_);
	njh::sys::requireExternalProgramsThrow(VecStr{"bwa", "bamtools", "samtools"});
	bfs::path inputSingles = njh::files::make_path(trimStub.string() + ".fastq");
	njh::files::checkExistenceThrow(genomeFnp,__PRETTY_FUNCTION__);
	bioRunner.RunBwaIndex(genomeFnp);
	if (!bfs::exists(inputSingles)) {
		inputSingles = njh::files::make_path(trimStub.string() + ".fastq.gz");
		if(!bfs::exists(inputSingles)){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error, one of the following pairs  " << "\n"
					<< " " << njh::files::make_path(trimStub.string() + ".fastq")  << "\n"
					<< " or " << "\n"
					<< " " << njh::files::make_path(trimStub.string() + ".fastq.gz")
					<< "\n" << "have to exist" << "\n";
			throw std::runtime_error { ss.str() };
		}
	}
	if("" == outputFnp){
		outputFnp = njh::files::make_path(outputDir, sampName + ".sorted.bam");
	}
	bfs::path outputFnpBai = outputFnp.string() + ".bai";
	bool needToRun = true;
	if(bfs::exists(outputFnp) &&
			bfs::exists(outputFnpBai) &&
			njh::files::firstFileIsOlder(inputSingles, outputFnp) ){
		needToRun = false;
	}
	if(force){
		needToRun = true;
	}
	bfs::path singlesSortedBam = njh::files::make_path(outputDir, trimStub.filename().string() + ".sorted.bam");


	std::string bNameStub = trimStub.filename().string();
	std::stringstream singlesCmd;
	singlesCmd << "bwa mem  -M -t " << numThreads
			<< " -R " << R"("@RG\tID:)" << bNameStub << "" << R"(\tSM:)"
			<< sampName << R"(")"
			<< " "   << extraBwaArgs
			<< " "   << genomeFnp
			<< " "   << inputSingles
			<< " 2> " << bfs::path(singlesSortedBam.string() + ".bwa.log")
			<< " | samtools sort -@ " << numThreads << " -o " << outputFnp
			<< " && samtools index " << outputFnp;

	if(needToRun){
		bfs::path logFnp = njh::files::make_path(outputDir, "alignTrimoOutputs_" + sampName + "_" + njh::getCurrentDate() + "_log.json");
		logFnp = njh::files::findNonexitantFile(logFnp);
		OutOptions logOpts(logFnp);
		std::ofstream logFile;
		logOpts.openFile(logFile);
		std::unordered_map<std::string, njh::sys::RunOutput> runOutputs;
		if(bfs::exists(inputSingles)){
			auto singlesRunOutput = njh::sys::run({singlesCmd.str()});
			BioCmdsUtils::checkRunOutThrow(singlesRunOutput, __PRETTY_FUNCTION__);
			runOutputs["bwa-singles"] = singlesRunOutput;
		}
		logFile << njh::json::toJson(runOutputs) << std::endl;
	}
	return 0;
}


int programWrapperRunner::runBwa(const njh::progutils::CmdArgs & inputCommands){
	bfs::path outputDir = "";
	bfs::path pairR1 = "";
	bfs::path pairR2 = "";
	bfs::path singles = "";

	bool force = false;
	uint32_t numThreads = 1;
	std::string sampName = "";
	bool removeIntermediateFiles = false;
	bfs::path genomeFnp = "";
	std::string extraBwaArgs = "";
	bfs::path outputFnp = "";
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(genomeFnp, "--genomeFnp", "genome file to align to", true);
	bool r1Set = setUp.setOption(pairR1, "--pairR1", "AdapterRemoval output stub");
	setUp.setOption(pairR2, "--pairR2", "AdapterRemoval output stub", r1Set);
	setUp.setOption(singles, "--single", "AdapterRemoval output stub", !r1Set);

	setUp.setOption(sampName, "--sampName", "Sample Name to give to final bam", true);
	setUp.setOption(force, "--force", "force run even if file already exists");
	setUp.setOption(outputFnp, "--outputFnp", "output name, will default to sampName.sorted.bam");
	setUp.setOption(extraBwaArgs, "--extraBwaArgs", "extra bwa arguments");
	setUp.setOption(numThreads, "--numThreads", "Number of threads");
	setUp.setOption(outputDir, "--outputDir", "output directory");
	setUp.setOption(removeIntermediateFiles, "--removeIntermediateFiles", "remove the intermediate bam files and just keep the ");
	setUp.finishSetUp(std::cout);
	BioCmdsUtils bioRunner(setUp.pars_.verbose_);
	njh::sys::requireExternalProgramsThrow(VecStr{"bwa", "bamtools", "samtools"});
	bfs::path inputSingles = singles;
	bfs::path inputPairedFirstMates = pairR1;
	bfs::path inputPairedSecondMates = pairR2;
	njh::files::checkExistenceThrow(genomeFnp,__PRETTY_FUNCTION__);
	bioRunner.RunBwaIndex(genomeFnp);

	if("" == outputFnp){
		outputFnp = njh::files::make_path(outputDir, sampName + ".sorted.bam");
	}
	bfs::path outputFnpBai = outputFnp.string() + ".bai";
	bool needToRun = true;
	if(bfs::exists(outputFnp) &&
			bfs::exists(outputFnpBai) &&
			njh::files::firstFileIsOlder(inputPairedFirstMates, outputFnp) &&
			njh::files::firstFileIsOlder(inputPairedSecondMates, outputFnp)){
		needToRun = false;
	}
	if(force){
		needToRun = true;
	}
	bfs::path singlesSortedBam = njh::files::make_path(outputDir, singles.string() + ".sorted.bam");
	bfs::path pairedSortedBam = njh::files::make_path(outputDir,  pairR1.string() + ".sorted.bam");


	std::stringstream singlesCmd;
	singlesCmd << "bwa mem  -M -t " << numThreads
			<< " -R " << R"("@RG\tID:)" << bfs::basename(singles) << "" << R"(\tSM:)"
			<< sampName << R"(")"
			<< " "   << extraBwaArgs
			<< " "   << genomeFnp
			<< " "   << inputSingles
			<< " 2> " << bfs::path(singlesSortedBam.string() + ".bwa.log")
			<< " | samtools sort -@ " << numThreads << " -o " << singlesSortedBam;

	std::stringstream pairedCmd;
	pairedCmd << "bwa mem  -M -t " << numThreads
			<< " -R " << R"("@RG\tID:)" << bfs::basename(pairR1) << "" << R"(\tSM:)"
			<< sampName << R"(")"
			<< " " << extraBwaArgs
			<< " " << genomeFnp
			<< " " << inputPairedFirstMates
			<< " " << inputPairedSecondMates
			<< " 2> " << bfs::path(pairedSortedBam.string() + ".bwa.log")
			<< " | samtools sort -@ " << numThreads << " -o " << pairedSortedBam;


	std::stringstream bamtoolsMergeAndIndexCmd;
	bamtoolsMergeAndIndexCmd << "bamtools merge " << " -in " << pairedSortedBam;
	if(bfs::exists(inputSingles)){
		bamtoolsMergeAndIndexCmd << " -in " <<  singlesSortedBam;
	}
	bamtoolsMergeAndIndexCmd << " -out " << outputFnp
			<< " && samtools index " << outputFnp;

	if(needToRun){
		bfs::path logFnp = njh::files::make_path(outputDir, "alignTrimoOutputs_" + sampName + "_" + njh::getCurrentDate() + "_log.json");
		logFnp = njh::files::findNonexitantFile(logFnp);
		OutOptions logOpts(logFnp);
		std::ofstream logFile;
		logOpts.openFile(logFile);
		std::unordered_map<std::string, njh::sys::RunOutput> runOutputs;
		if(bfs::exists(inputSingles)){
			auto singlesRunOutput = njh::sys::run({singlesCmd.str()});
			BioCmdsUtils::checkRunOutThrow(singlesRunOutput, __PRETTY_FUNCTION__);
			runOutputs["bwa-singles"] = singlesRunOutput;
		}

		auto pairedRunOutput = njh::sys::run({pairedCmd.str()});
		BioCmdsUtils::checkRunOutThrow(pairedRunOutput, __PRETTY_FUNCTION__);
		runOutputs["bwa-paired"] = pairedRunOutput;
		if(bfs::exists(singlesSortedBam) ){
			auto bamtoolsMergeAndIndexRunOutput = njh::sys::run({bamtoolsMergeAndIndexCmd.str()});
			BioCmdsUtils::checkRunOutThrow(bamtoolsMergeAndIndexRunOutput, __PRETTY_FUNCTION__);
			runOutputs["bamtools-merge-index"] = bamtoolsMergeAndIndexRunOutput;
		}else{
			std::stringstream ss;
			ss << "samtools index " << outputFnp;
			bfs::rename(pairedSortedBam, outputFnp);
			auto indexRunOutput = njh::sys::run({ss.str()});
			BioCmdsUtils::checkRunOutThrow(indexRunOutput, __PRETTY_FUNCTION__);
			runOutputs["index"] = indexRunOutput;
		}
		logFile << njh::json::toJson(runOutputs) << std::endl;
		if(removeIntermediateFiles){
			if(bfs::exists(singlesSortedBam)){
				bfs::remove(pairedSortedBam);
				bfs::remove(singlesSortedBam);
			}
		}
	}
	return 0;
}

int programWrapperRunner::runBwaOnAdapterReomvalOutputSinglesCombined(const njh::progutils::CmdArgs & inputCommands){
	bfs::path outputDir = "";
	bfs::path trimStub = "";
	bool force = false;
	uint32_t numThreads = 1;
	std::string sampName = "";
	bool removeIntermediateFiles = false;
	bfs::path genomeFnp = "";
	std::string extraBwaArgs = "";
	bfs::path outputFnp = "";
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(genomeFnp, "--genomeFnp", "genome file to align to", true);
	setUp.setOption(trimStub, "--trimStub", "AdapterRemoval output stub", true);
	setUp.setOption(sampName, "--sampName", "Sample Name to give to final bam", true);
	setUp.setOption(force, "--force", "force run even if file already exists");
	setUp.setOption(outputFnp, "--outputFnp", "output name, will default to sampName.sorted.bam");
	setUp.setOption(extraBwaArgs, "--extraBwaArgs", "extra bwa arguments");
	setUp.setOption(numThreads, "--numThreads", "Number of threads");
	setUp.setOption(outputDir, "--outputDir", "output directory");
	setUp.setOption(removeIntermediateFiles, "--removeIntermediateFiles", "remove the intermediate bam files and just keep the ");
	setUp.finishSetUp(std::cout);
	BioCmdsUtils bioRunner(setUp.pars_.verbose_);
	njh::sys::requireExternalProgramsThrow(VecStr{"bwa", "bamtools", "samtools"});
	bfs::path inputSingles = njh::files::make_path(trimStub.string() + "_singles.fastq");
	bfs::path inputPairedFirstMates = njh::files::make_path(trimStub.string() + "_1.fastq");
	bfs::path inputPairedSecondMates = njh::files::make_path(trimStub.string() + "_2.fastq");
	njh::files::checkExistenceThrow(genomeFnp,__PRETTY_FUNCTION__);
	bioRunner.RunBwaIndex(genomeFnp);

	if (!bfs::exists(inputSingles)
			&& !bfs::exists(inputPairedFirstMates)
			&& !bfs::exists(inputPairedSecondMates)) {
		inputPairedFirstMates = njh::files::make_path(trimStub.string() + "_1.fastq.gz");
		inputPairedSecondMates = njh::files::make_path(trimStub.string() + "_2.fastq.gz");
		inputSingles = njh::files::make_path(trimStub.string() + "_singles.fastq.gz");
		if(!bfs::exists(inputPairedFirstMates) &&
				!bfs::exists(inputPairedSecondMates)){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error, one of the following pairs  " << "\n"
					<< " " << njh::files::make_path(trimStub.string() + "_1.fastq") << " "
					       << njh::files::make_path(trimStub.string() + "_2.fastq") << "\n"
					<< " or " << "\n"
					<< " " << njh::files::make_path(trimStub.string() + "_1.fastq.gz")
					<< " " << njh::files::make_path(trimStub.string() + "_2.fastq.gz")
					<< "\n" << "have to exist" << "\n";
			throw std::runtime_error { ss.str() };
		}
	}
	if("" == outputFnp){
		outputFnp = njh::files::make_path(outputDir, sampName + ".sorted.bam");
	}
	bfs::path outputFnpBai = outputFnp.string() + ".bai";
	bool needToRun = true;
	if(bfs::exists(outputFnp) &&
			bfs::exists(outputFnpBai) &&
			njh::files::firstFileIsOlder(inputPairedFirstMates, outputFnp) &&
			njh::files::firstFileIsOlder(inputPairedSecondMates, outputFnp)){
		needToRun = false;
	}
	if(force){
		needToRun = true;
	}
	bfs::path singlesSortedBam = njh::files::make_path(outputDir, trimStub.filename().string() + "_singles.sorted.bam");
	bfs::path pairedSortedBam = njh::files::make_path(outputDir,trimStub.filename().string() + ".sorted.bam");


	std::string bNameStub = trimStub.filename().string();
	std::stringstream singlesCmd;
	singlesCmd << "bwa mem  -M -t " << numThreads
			<< " -R " << R"("@RG\tID:)" << bNameStub << "_singles" << R"(\tSM:)"
			<< sampName << R"(")"
			<< " "   << extraBwaArgs
			<< " "   << genomeFnp
			<< " "   << inputSingles
			<< " 2> " << bfs::path(singlesSortedBam.string() + ".bwa.log")
			<< " | samtools sort -@ " << numThreads << " -o " << singlesSortedBam;

	std::stringstream pairedCmd;
	pairedCmd << "bwa mem  -M -t " << numThreads
			<< " -R " << R"("@RG\tID:)" << bNameStub << "" << R"(\tSM:)"
			<< sampName << R"(")"
			<< " " << extraBwaArgs
			<< " " << genomeFnp
			<< " " << inputPairedFirstMates
			<< " " << inputPairedSecondMates
			<< " 2> " << bfs::path(pairedSortedBam.string() + ".bwa.log")
			<< " | samtools sort -@ " << numThreads << " -o " << pairedSortedBam;


	std::stringstream bamtoolsMergeAndIndexCmd;
	bamtoolsMergeAndIndexCmd << "bamtools merge " << " -in " << pairedSortedBam;
	if(bfs::exists(inputSingles)){
		bamtoolsMergeAndIndexCmd << " -in " <<  singlesSortedBam;
	}
	bamtoolsMergeAndIndexCmd << " -out " << outputFnp
			<< " && samtools index " << outputFnp;

	if(needToRun){
		bfs::path logFnp = njh::files::make_path(outputDir, "alignTrimoOutputs_" + sampName + "_" + njh::getCurrentDate() + "_log.json");
		logFnp = njh::files::findNonexitantFile(logFnp);
		OutOptions logOpts(logFnp);
		std::ofstream logFile;
		logOpts.openFile(logFile);
		std::unordered_map<std::string, njh::sys::RunOutput> runOutputs;
		if(bfs::exists(inputSingles)){
			auto singlesRunOutput = njh::sys::run({singlesCmd.str()});
			BioCmdsUtils::checkRunOutThrow(singlesRunOutput, __PRETTY_FUNCTION__);
			runOutputs["bwa-singles"] = singlesRunOutput;
		}

		auto pairedRunOutput = njh::sys::run({pairedCmd.str()});
		BioCmdsUtils::checkRunOutThrow(pairedRunOutput, __PRETTY_FUNCTION__);
		runOutputs["bwa-paired"] = pairedRunOutput;
		if(bfs::exists(singlesSortedBam) ){
			auto bamtoolsMergeAndIndexRunOutput = njh::sys::run({bamtoolsMergeAndIndexCmd.str()});
			BioCmdsUtils::checkRunOutThrow(bamtoolsMergeAndIndexRunOutput, __PRETTY_FUNCTION__);
			runOutputs["bamtools-merge-index"] = bamtoolsMergeAndIndexRunOutput;
		}else{
			std::stringstream ss;
			ss << "samtools index " << outputFnp;
			bfs::rename(pairedSortedBam, outputFnp);
			auto indexRunOutput = njh::sys::run({ss.str()});
			BioCmdsUtils::checkRunOutThrow(indexRunOutput, __PRETTY_FUNCTION__);
			runOutputs["index"] = indexRunOutput;
		}
		logFile << njh::json::toJson(runOutputs) << std::endl;
		if(removeIntermediateFiles){
			if(bfs::exists(singlesSortedBam)){
				bfs::remove(pairedSortedBam);
				bfs::remove(singlesSortedBam);
			}
		}
	}
	return 0;
}


int programWrapperRunner::runBowtieOnAdapterReomvalOutputSinglesCombined(const njh::progutils::CmdArgs & inputCommands){
	bfs::path outputDir = "";
	bfs::path trimStub = "";
	bool force = false;
	uint32_t numThreads = 1;
	std::string sampName = "";
	bool removeIntermediateFiles = false;
	bfs::path genomePrefix = "";
	std::string extraBowtie2Args = "";
	bfs::path outputFnp = "";
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(genomePrefix, "--genomePrefix", "genome prefix for bowtie2", true);
	setUp.setOption(trimStub, "--trimStub", "AdapterRemoval output stub", true);
	setUp.setOption(sampName, "--sampName", "Sample Name to give to final bam", true);
	setUp.setOption(force, "--force", "force run even if file already exists");
	setUp.setOption(outputFnp, "--outputFnp", "output name, will default to sampName.sorted.bam");
	setUp.setOption(extraBowtie2Args, "--extraBowtie2Args", "extra bowtie2 arguments");
	setUp.setOption(numThreads, "--numThreads", "Number of threads");
	setUp.setOption(outputDir, "--outputDir", "output directory");
	setUp.setOption(removeIntermediateFiles, "--removeIntermediateFiles", "remove the intermediate bam files and just keep the ");
	setUp.finishSetUp(std::cout);
	BioCmdsUtils bioRunner(setUp.pars_.verbose_);
	njh::sys::requireExternalProgramsThrow(VecStr{"bowtie2", "bamtools", "samtools"});
	bfs::path inputSingles = njh::files::make_path(trimStub.string() + "_singles.fastq");
	bfs::path inputPairedFirstMates = njh::files::make_path(trimStub.string() + "_1.fastq");
	bfs::path inputPairedSecondMates = njh::files::make_path(trimStub.string() + "_2.fastq");
	bfs::path genomeFnp = genomePrefix.string() + ".fasta";
	njh::files::checkExistenceThrow(genomeFnp,__PRETTY_FUNCTION__);
	if(setUp.pars_.debug_){
		std::cout << "genomeFnp: " << genomeFnp << std::endl;
	}
	bioRunner.RunBowtie2Index(genomeFnp);

	if (!bfs::exists(inputSingles)
			&& !bfs::exists(inputPairedFirstMates)
			&& !bfs::exists(inputPairedSecondMates)) {
		inputPairedFirstMates = njh::files::make_path(trimStub.string() + "_1.fastq.gz");
		inputPairedSecondMates = njh::files::make_path(trimStub.string() + "_2.fastq.gz");
		inputSingles = njh::files::make_path(trimStub.string() + "_singles.fastq.gz");
		if(!bfs::exists(inputPairedFirstMates) &&
				!bfs::exists(inputPairedSecondMates)){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error, one of the following pairs  " << "\n"
					<< " " << njh::files::make_path(trimStub.string() + "_1.fastq") << " "
					       << njh::files::make_path(trimStub.string() + "_2.fastq") << "\n"
					<< " or " << "\n"
					<< " " << njh::files::make_path(trimStub.string() + "_1.fastq.gz")
					<< " " << njh::files::make_path(trimStub.string() + "_2.fastq.gz")
					<< "\n" << "have to exist" << "\n";
			throw std::runtime_error { ss.str() };
		}
	}
	if("" == outputFnp){
		outputFnp = njh::files::make_path(outputDir, sampName + ".sorted.bam");
	}
	bfs::path outputFnpBai = outputFnp.string() + ".bai";
	bool needToRun = true;
	if(bfs::exists(outputFnp) &&
			bfs::exists(outputFnpBai) &&
			njh::files::firstFileIsOlder(inputPairedFirstMates, outputFnp) &&
			njh::files::firstFileIsOlder(inputPairedSecondMates, outputFnp)){
		needToRun = false;
	}
	if(force){
		needToRun = true;
	}
	bfs::path singlesSortedBam = njh::files::make_path(outputDir, trimStub.filename().string() + "_singles.sorted.bam");
	bfs::path pairedSortedBam = njh::files::make_path(outputDir,trimStub.filename().string() + ".sorted.bam");


	std::string bNameStub = trimStub.filename().string();
	std::stringstream singlesCmd;
	singlesCmd << "bowtie2 "
			<< " -threads " << numThreads
			<< " --rg-id " << bNameStub << "_singles "
			<< " --rg \"SM:" << sampName << "\""
			<< " "   << extraBowtie2Args
			<< " -x "   << genomePrefix
			<< " -U "   << inputSingles
			<< " 2> " << bfs::path(singlesSortedBam.string() + ".bowtie2.log")
			<< " | samtools sort -@ " << numThreads << " -o " << singlesSortedBam;

	std::stringstream pairedCmd;
	pairedCmd << "bowtie2 "
			<< " -threads " << numThreads
			<< " --rg-id " << bNameStub
			<< " --rg \"SM:" << sampName << "\""
			<< " "   << extraBowtie2Args
			<< " -x "   << genomePrefix
			<< " -1 "   << inputPairedFirstMates
			<< " -2 "   << inputPairedSecondMates
			<< " 2> " << bfs::path(pairedSortedBam.string() + ".bowtie2.log")
			<< " | samtools sort -@ " << numThreads << " -o " << pairedSortedBam;


	std::stringstream bamtoolsMergeAndIndexCmd;
	bamtoolsMergeAndIndexCmd << "bamtools merge " << " -in " << pairedSortedBam;
	bool inputSingleEmpty = false;
	if(bfs::exists(inputSingles)){
		uint32_t count  = 0;
		InputStream singlesInput(inputSingles);
		std::string line = "";
		while(njh::files::crossPlatGetline(singlesInput, line)){
			++count;
			break;
		}
		if(0 == count){
			inputSingleEmpty = true;
		}
	}
	if(bfs::exists(inputSingles) && !inputSingleEmpty){
		bamtoolsMergeAndIndexCmd << " -in " <<  singlesSortedBam;
	}
	bamtoolsMergeAndIndexCmd << " -out " << outputFnp
			<< " && samtools index " << outputFnp;

	if(needToRun){
		bfs::path logFnp = njh::files::make_path(outputDir, "alignTrimoOutputs_" + sampName + "_" + njh::getCurrentDate() + "_log.json");
		logFnp = njh::files::findNonexitantFile(logFnp);
		OutOptions logOpts(logFnp);
		std::ofstream logFile;
		logOpts.openFile(logFile);
		std::unordered_map<std::string, njh::sys::RunOutput> runOutputs;
		if(bfs::exists(inputSingles) && !inputSingleEmpty){
			auto singlesRunOutput = njh::sys::run({singlesCmd.str()});
			BioCmdsUtils::checkRunOutThrow(singlesRunOutput, __PRETTY_FUNCTION__);
			runOutputs["bowtie2-singles"] = singlesRunOutput;
		}

		auto pairedRunOutput = njh::sys::run({pairedCmd.str()});
		BioCmdsUtils::checkRunOutThrow(pairedRunOutput, __PRETTY_FUNCTION__);
		runOutputs["bowtie2-paired"] = pairedRunOutput;
		if(bfs::exists(singlesSortedBam) ){
			auto bamtoolsMergeAndIndexRunOutput = njh::sys::run({bamtoolsMergeAndIndexCmd.str()});
			BioCmdsUtils::checkRunOutThrow(bamtoolsMergeAndIndexRunOutput, __PRETTY_FUNCTION__);
			runOutputs["bamtools-merge-index"] = bamtoolsMergeAndIndexRunOutput;
		}else{
			std::stringstream ss;
			ss << "samtools index " << outputFnp;
			bfs::rename(pairedSortedBam, outputFnp);
			auto indexRunOutput = njh::sys::run({ss.str()});
			BioCmdsUtils::checkRunOutThrow(indexRunOutput, __PRETTY_FUNCTION__);
			runOutputs["index"] = indexRunOutput;
		}
		logFile << njh::json::toJson(runOutputs) << std::endl;
		if(removeIntermediateFiles){
			if(bfs::exists(singlesSortedBam)){
				bfs::remove(pairedSortedBam);
				bfs::remove(singlesSortedBam);
			}
		}
	}
	return 0;
}

} // namespace njhseq
