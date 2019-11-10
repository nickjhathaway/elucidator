/*
 * readSimulatorRunner_simIllumina.cpp
 *
 *  Created on: Aug 24, 2019
 *      Author: nicholashathaway
 */


#include "readSimulatorRunner.hpp"
#include "elucidator/simulation.h"
#include <SeekDeep/objects/PrimersAndMids.hpp>


namespace njhseq {




int readSimulatorRunner::simIlluminaOnSeqs(const njh::progutils::CmdArgs & inputCommands) {
  readSimulatorSetUp setUp(inputCommands);
	bool singleEnd = false;
	bool nonGz = false;
	bfs::path illuminaProfileDir = "";
	OutOptions outOpts(bfs::path("out"));
	uint32_t outLength = 150;
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processReadInNames(VecStr{"--fasta", "--fastagz", "--fastq", "--fastqgz"},true);
	setUp.setOption(singleEnd, "--singleEnd", "Single End");
	setUp.setOption(nonGz, "--nonGz", "do not compress the output fastqs");
	setUp.setOption(illuminaProfileDir, "--illuminaProfileDir", "Illumina Profile Dir", true);
	setUp.setOption(outLength, "--outLength", "Illumina Length to simulate", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	RoughIlluminaSimulator simulator(illuminaProfileDir);

	if(outLength > simulator.r1Profile_.errorRates_.size()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "out length: " << outLength << " longer than able to simulate, " << simulator.r1Profile_.errorRates_.size()<< "\n";
		throw std::runtime_error{ss.str()};
	}
	seqInfo seq;
	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	SeqIOOptions outSeqOpts;
	if(singleEnd){
		if(nonGz){
			outSeqOpts = SeqIOOptions::genFastqOut(outOpts.outName());
			outSeqOpts.out_.transferOverwriteOpts(outOpts);
		}else{
			outSeqOpts = SeqIOOptions::genFastqOutGz(outOpts.outName());
			outSeqOpts.out_.transferOverwriteOpts(outOpts);
		}
	}else{
		if(nonGz){
			outSeqOpts = SeqIOOptions::genPairedOut(outOpts.outName());
			outSeqOpts.out_.transferOverwriteOpts(outOpts);
		}else{
			outSeqOpts = SeqIOOptions::genPairedOutGz(outOpts.outName());
			outSeqOpts.out_.transferOverwriteOpts(outOpts);
		}
	}
	SeqOutput writer(outSeqOpts);
	writer.openOut();
	while(reader.readNextRead(seq)){
		auto r1Seq = simulator.simR1(seq, outLength);
		if(singleEnd){
			writer.write(r1Seq);
		}else{
			seq.reverseComplementRead(false, true);
			auto r2Seq = simulator.simR2(seq, outLength);
			writer.write(PairedRead(r1Seq, r2Seq, false));

		}
	}

	return 0;
}


int readSimulatorRunner::shearSequences(const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	uint32_t mean = 240;
	uint32_t std = 81;//245.431 140.977
	uint32_t minLen = 8;
	uint32_t readNumber = 100;
	setUp.processDefaultReader();
	setUp.setOption(mean, "--mean", "Mean");
	setUp.setOption(minLen, "--minLen", "Minimum Length");
	setUp.setOption(std, "--std", "Standard Deviation");
	setUp.setOption(readNumber, "--readNumber", "Number of fragments to create");

	setUp.finishSetUp(std::cout);
	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	reader.openOut();
	sim::ReadLenNormalDistribution uiNormDist(mean, std);
	njh::randomGenerator gen;
	seqInfo seq;
	while(reader.readNextRead(seq)){
		for (uint32_t i = 0; i < readNumber; ++i) {
			auto readPositionsLens = sim::genFragPosSizes(len(seq), uiNormDist, gen);
			auto fragments = sim::genFragments(seq.seq_, readPositionsLens);
			for(const auto & fragPos : iter::range(fragments.size())){
				if(readPositionsLens[fragPos].second >= minLen){
					seqInfo fragInfo(
							njh::pasteAsStr(seq.name_, "_id", i ,"_pos" ,readPositionsLens[fragPos].first ,"_size" ,readPositionsLens[fragPos].second),
							fragments[fragPos]);
					reader.write(fragInfo);
				}
			}
		}
	}
	return 0;
}

int readSimulatorRunner::shearSimIlluminaAlign(const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	uint32_t mean = 240;
	uint32_t std = 81;//245.431 140.977
	uint32_t minLen = 8;
	uint32_t readNumber = 100;
	bfs::path illuminaProfileDir = njh::files::make_path(elucidator_INSTALLDIR, "etc/illumina_profiles/miseq_250");
	uint32_t outLength = 150;
	uint32_t numThreads = 1;
	bfs::path genomeFnp;
	std::string sampleName = "";
	setUp.setOption(illuminaProfileDir, "--illuminaProfileDir", njh::pasteAsStr("Illumina Profile Dir, by default will simulate miseq 250 paired end, see ", illuminaProfileDir, " for examples on the input profile files"));
	setUp.setOption(outLength, "--outLength", "Illumina Length to simulate", true);
	setUp.setOption(sampleName, "--sampleName", "Sample Name to Give Bam", true);
	setUp.setOption(genomeFnp, "--genomeFnp", "Genome File name to align to", true);

	setUp.processReadInNames();
	setUp.setOption(numThreads, "--numThreads", "num Threads to use");

	setUp.setOption(mean, "--mean", "Mean Fragment Size");
	setUp.setOption(minLen, "--minLen", "Minimum Length for the fragments");
	setUp.setOption(std, "--std", "Standard Deviation of the fragment size");
	setUp.setOption(readNumber, "--readNumber", "Number of reads to create fragments off of per sequence, proxy for read depth");
	setUp.processDirectoryOutputName(true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);

	njh::sys::requireExternalProgramThrow("AdapterRemoval");
	njh::sys::requireExternalProgramThrow("elucidator");
	njh::sys::requireExternalProgramThrow("bwa");
	njh::sys::requireExternalProgramThrow("picard");


	RoughIlluminaSimulator simulator(illuminaProfileDir);


	if(outLength > simulator.r1Profile_.errorRates_.size()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "out length: " << outLength << " longer than able to simulate, " << simulator.r1Profile_.errorRates_.size()<< "\n";
		throw std::runtime_error{ss.str()};
	}

	//shearing
	auto outSheared = SeqIOOptions::genFastaOut(njh::files::make_path(setUp.pars_.directoryName_, "shearedSeqs"));
	{
		SeqInput reader(setUp.pars_.ioOptions_);
		reader.openIn();
		SeqOutput writer(outSheared);
		writer.openOut();
		sim::ReadLenNormalDistribution uiNormDist(mean, std);
		njh::randomGenerator gen;
		seqInfo seq;
		while(reader.readNextRead(seq)){
			for (uint32_t i = 0; i < readNumber; ++i) {
				auto readPositionsLens = sim::genFragPosSizes(len(seq), uiNormDist, gen);
				auto fragments = sim::genFragments(seq.seq_, readPositionsLens);
				for(const auto & fragPos : iter::range(fragments.size())){
					if(readPositionsLens[fragPos].second >= minLen){
						seqInfo fragInfo(
								njh::pasteAsStr(seq.name_, "_id", i ,"_pos" ,readPositionsLens[fragPos].first ,"_size" ,readPositionsLens[fragPos].second),
								fragments[fragPos]);
						if(gen.unifRand() > 0.5){
							fragInfo.reverseComplementRead(false, true);
						}
						writer.write(fragInfo);
					}
				}
			}
		}
	}

	//illumina
	auto illuminaPairedOutStub = njh::files::make_path(setUp.pars_.directoryName_, "fastq/seqeuncedShearedSeqs");
	njh::files::makeDir(njh::files::MkdirPar{njh::files::make_path(setUp.pars_.directoryName_, "fastq")});
	{
		seqInfo seq;
		SeqInput reader(SeqIOOptions::genFastaIn(outSheared.out_.outName()));
		reader.openIn();
		SeqIOOptions outSeqOpts = SeqIOOptions::genPairedOutGz(illuminaPairedOutStub);
		SeqOutput writer(outSeqOpts);
		writer.openOut();
		while(reader.readNextRead(seq)){
			auto r1Seq = simulator.simR1(seq, outLength);
			seq.reverseComplementRead(false, true);
			auto r2Seq = simulator.simR2(seq, outLength);
			writer.write(PairedRead(r1Seq, r2Seq, false));
		}
	}
	auto tempDir = njh::files::make_path(setUp.pars_.directoryName_, "temp");
	njh::files::makeDir(njh::files::MkdirPar{tempDir});

	auto bamsDir = njh::files::make_path(setUp.pars_.directoryName_, "bams");
	njh::files::makeDir(njh::files::MkdirPar{bamsDir});

	std::stringstream detectPrimersCmd;
	detectPrimersCmd << "elucidator detectPossiblePrimers --overWriteDir --fastq1gz " << illuminaPairedOutStub << "_R1.fastq.gz --fastq2gz " << illuminaPairedOutStub << "_R2.fastq.gz --dout " << tempDir << "/seqeuncedShearedSeqs_detectPrimers";
	auto detectPrimersCmd_runOut = njh::sys::run({detectPrimersCmd.str()});
	BioCmdsUtils::checkRunOutThrow(detectPrimersCmd_runOut, __PRETTY_FUNCTION__);

	std::stringstream setUpTrimmingCmd;
	setUpTrimmingCmd << "cd " << tempDir << " && elucidator setUpRunAdapterRemoval --extraArgs=\"--trim3p 1 --minlength 30\" --numThreads " << numThreads << " --fastq1 ../fastq/seqeuncedShearedSeqs_R1.fastq.gz --fastq2 ../fastq/seqeuncedShearedSeqs_R2.fastq.gz --overWrite --combineSingles --gzip --detectPrimersDir seqeuncedShearedSeqs_detectPrimers";
	auto setUpTrimmingCmd_runOut = njh::sys::run({setUpTrimmingCmd.str()});
	BioCmdsUtils::checkRunOutThrow(setUpTrimmingCmd_runOut, __PRETTY_FUNCTION__);

	std::stringstream runAlignmentCmd;
	runAlignmentCmd << "cd " << tempDir << " && elucidator runBwaOnAdapterReomvalOutputSinglesCombined --removeIntermediateFiles --trimStub trimmed_seqeuncedShearedSeqs --genomeFnp "<< genomeFnp << " --sampName " << sampleName <<" --numThreads " << numThreads << "";
	auto runAlignmentCmd_runOut = njh::sys::run({runAlignmentCmd.str()});
	BioCmdsUtils::checkRunOutThrow(runAlignmentCmd_runOut, __PRETTY_FUNCTION__);

	std::stringstream appendNameCmd;
	appendNameCmd << "cd " << tempDir << " && elucidator appendReadGroupToName --bam " << sampleName << ".sorted.bam";
	auto appendNameCmd_runOut = njh::sys::run({appendNameCmd.str()});
	BioCmdsUtils::checkRunOutThrow(appendNameCmd_runOut, __PRETTY_FUNCTION__);

	std::stringstream markDupsCmds;
	markDupsCmds << "cd " << tempDir << " && elucidator runPicardMarkDups --bamFnp renamed_" << sampleName << ".sorted.bam --outDir ../bams/ --outNameStub " << sampleName << " --logDir ./";
	auto markDupsCmds_runOut = njh::sys::run({markDupsCmds.str()});
	BioCmdsUtils::checkRunOutThrow(markDupsCmds_runOut, __PRETTY_FUNCTION__);

	return 0;
}




}  // namespace njhseq
