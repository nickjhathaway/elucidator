/*
 * programWrappers.cpp
 *
 *  Created on: Feb 2, 2017
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


programWrapperRunner::programWrapperRunner()
    : njh::progutils::ProgramRunner(
          {
					 addFunc("runTrinityOnRegion", runTrinityOnRegion, false),
					 addFunc("runSSAKEOnRegion", runSSAKEOnRegion, false),
					 addFunc("runLastz", runLastz, false),
					 addFunc("sraIsPairedEnd", sraIsPairedEnd, false),
					 addFunc("sraFastqDump", sraFastqDump, false),
					 addFunc("runTrimmomatic", runTrimmomatic, false),
					 addFunc("runBwaOnTrimmomaticOutputPE", runBwaOnTrimmomaticOutputPE, false),
					 addFunc("parsePrimer3OutputToJson", parsePrimer3OutputToJson, false),
					 addFunc("parsePrimer3OutputToBed", parsePrimer3OutputToBed, false),
					 addFunc("parsePrimer3OutputToPossibleMipArms", parsePrimer3OutputToPossibleMipArms, false),
					 addFunc("findNonUniquePrimerArms", findNonUniquePrimerArms, false),
					 addFunc("convertShorahSupportSeqs", convertShorahSupportSeqs, false),
					 addFunc("runShorahAmplian", runShorahAmplian, false),
					 addFunc("runSamtoolsFlagStat", runSamtoolsFlagStat, false),
					 addFunc("runAdapterRemoval", runAdapterRemoval, false),
					 addFunc("setUpRunAdapterRemoval", setUpRunAdapterRemoval, false),
					 addFunc("runAdapterRemovalSE", runAdapterRemovalSE, false),
					 addFunc("runBwaOnAdapterReomvalOutputSinglesCombined", runBwaOnAdapterReomvalOutputSinglesCombined, false),
					 addFunc("runBwaOnAdapterReomvalOutputSE", runBwaOnAdapterReomvalOutputSE, false),
					 addFunc("runBwa", runBwa, false),
					 addFunc("generatingPrime3TemplatesBasedOnMALN", generatingPrime3TemplatesBasedOnMALN, false),
					 addFunc("testHasProgram", testHasProgram, false),
					 addFunc("runBowtieOnAdapterReomvalOutputSinglesCombined", runBowtieOnAdapterReomvalOutputSinglesCombined, false),
           },//
          "programWrapper") {}
int programWrapperRunner::generatingPrime3TemplatesBasedOnMALN(const njh::progutils::CmdArgs & inputCommands) {
	bool muscle = false;
	Muscler::TrimWithMusclePars mPars;
	mPars.baseCutOff = 1;       // all bases in this position
	mPars.hardGapCutOff = 0;    // no gaps
	mPars.spanningCutOff = 1;   // covered by all seqs
	mPars.streakLenCutOff = 18; // at least 18 adjacent positions
	bool noExclusion = false;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.setOption(noExclusion, "--noExclusion", "No Exclusion");
	setUp.setOption(muscle, "--muscle", "Muscle the sequences if they haven't already been aligned");
	setUp.processReadInNames(true);
	setUp.processDirectoryOutputName(true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);

	SeqInput reader(setUp.pars_.ioOptions_);

	auto seqs = reader.readAllReads<seqInfo>();

	if(muscle){
		Muscler mRunner;
		readVec::removeGapsFromReads(seqs);
		mRunner.muscleSeqs(seqs);
	}

	auto totalInputSeqs = len(seqs);
	std::function<bool(const std::shared_ptr<Muscler::AlnPosScore> &)> scorePred = [&mPars,&totalInputSeqs](const std::shared_ptr<Muscler::AlnPosScore> & score){
		if(score->getBaseSpannedPerc()  >= mPars.baseCutOff &&
				score->gapCount_ <= mPars.hardGapCutOff &&
				score->getPercentOfSequencesSpanningPosition(totalInputSeqs) >= mPars.spanningCutOff){
			return true;
		}else{
			return false;
		}//" | "
	};

	auto streaks = Muscler::getAlignmentStreaksPositions(seqs, scorePred, mPars.streakLenCutOff);
	OutputStream outStreaks(OutOptions(njh::files::make_path(setUp.pars_.directoryName_, "streaks")));
	outStreaks << "start\tend\tlength" << std::endl;
	for(const auto & streak : streaks){
		outStreaks << streak.start_
				<< "\t" << streak.end_
				<< "\t" << streak.end_ - streak.start_
				<< std::endl;
	}
	auto alnSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(setUp.pars_.directoryName_, "aln_seqs.fasta"));
	SeqOutput::write(seqs, alnSeqOpts);
	auto unalignedSeqs = seqs;
	readVec::removeGapsFromReads(unalignedSeqs);
	if (streaks.size() > 1) {
		for (const auto & seqPos : iter::range(seqs.size())) {
			OutputStream outprimer3Template(OutOptions(njh::files::make_path(setUp.pars_.directoryName_, seqs[seqPos].name_ + "_primer3_template.txt")));
			outprimer3Template << "SEQUENCE_ID=" << seqs[seqPos].name_ << std::endl;// ideel-barcode-chr8-1100058_S0_Sub0_ext,CHR8:1099906-1100208
			outprimer3Template << "SEQUENCE_TEMPLATE=" << unalignedSeqs[seqPos].seq_ << std::endl;
			outprimer3Template << "SEQUENCE_TARGET=" << std::endl;
			std::string exclusion = "";
//SEQUENCE_TEMPLATE=TGTGATTTATGTAATAATCCTATAAGTCCGTTATGTTATGTGTATGAATGTAACATATGTGATAATTTTGCTTTGTGTAAGAAATGTTATAAAAAAAATAAGCATGAACACAATTTAAAAAAAATATTAGTACCCAGACATTGCATACCGCCACAGGATTATCAAAATGAAGAATTGATAGCAAAAGATGGACAAATTAATAATAATAATAATAATAATAATAATAATAATAATAATAATAATATATATGATAATAATTTAATGGAATACCTAGAAAATGATTCAAGTGCATATGAAGATCTG
//SEQUENCE_TARGET=
//SEQUENCE_EXCLUDED_REGION=38,1 59,1 116,1 149,18 197,31"
			for (const auto & streakPos : iter::range<uint32_t>(0, streaks.size() - 1)) {
				auto start = getRealPosForAlnPos(seqs[seqPos].seq_,
						streaks[streakPos].end_);
				auto end = getRealPosForAlnPos(seqs[seqPos].seq_,
						streaks[streakPos + 1].start_);
				if("" != exclusion){
					exclusion.append(" ");
				}
				exclusion.append(njh::pasteAsStr(start, ",",end - start));
				changeSubStrToLower(unalignedSeqs[seqPos].seq_, start, end - start);
			}
			if(!noExclusion){
				outprimer3Template << "SEQUENCE_EXCLUDED_REGION=" << exclusion<< std::endl;
			}
			outprimer3Template << "=" << std::endl;
		}
	}
	auto unalnSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(setUp.pars_.directoryName_, "unaln_seqs.fasta"));
	SeqOutput::write(unalignedSeqs, unalnSeqOpts);



	return 0;
}



int programWrapperRunner::testHasProgram(const njh::progutils::CmdArgs & inputCommands) {
	std::string program = "";
	seqSetUp setUp(inputCommands);
	setUp.setOption(program, "--program", "program to test for");
	setUp.finishSetUp(std::cout);
	auto output = njh::sys::run({"which", program});

	std::cout << njh::colorBool(njh::sys::hasSysCommand(program)) << std::endl;
	std::cout << output.toJson() << std::endl;

	return 0;
}


int programWrapperRunner::runSamtoolsFlagStat(
		const njh::progutils::CmdArgs & inputCommands) {
	bool keepIndividualFiles = false;
	OutOptions jsonOutOpts(bfs::path("samtoolsFlagStats"), ".json");
	std::vector<bfs::path> bams;
	uint32_t numThreads = 1;
	seqSetUp setUp(inputCommands);
	setUp.pars_.directoryName_ = "";
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processWritingOptions(jsonOutOpts);
	setUp.setOption(bams, "--bams", "Bam files to run on", true);
	setUp.setOption(numThreads, "--numThreads", "Number of threads for samtools to use");
	setUp.setOption(keepIndividualFiles, "--keepIndividualFiles", "Keep Individual Files");
	setUp.finishSetUp(std::cout);

	njh::sys::requireExternalProgramThrow("samtools");
	jsonOutOpts.throwIfOutExistsNoOverWrite(__PRETTY_FUNCTION__);
	OutputStream jsonOut(jsonOutOpts);
	Json::Value allValues;
	std::vector<bfs::path> samtoolsFlagFiles;
	for(const auto & bam : bams){
		njh::files::checkExistenceThrow(bam, __PRETTY_FUNCTION__);
		OutOptions samtoolsOutputOpts(bfs::path("samtools_flagstat_out_for_" + bam.filename().string() ));
		samtoolsOutputOpts.transferOverwriteOpts(jsonOutOpts);
		std::stringstream samtoolsCmd;
		samtoolsCmd << "samtools flagstat -@ " << numThreads << " " << bam << " > " << samtoolsOutputOpts.outName();
		samtoolsFlagFiles.emplace_back(samtoolsOutputOpts.outName());

		auto runOutput = njh::sys::run({samtoolsCmd.str()});
		BioCmdsUtils::checkRunOutThrow(runOutput, __PRETTY_FUNCTION__	);

		InputStream in(InOptions(samtoolsOutputOpts.outName()));
		std::string line = "";
		Json::Value outVal;
		outVal["bamFnp"] = njh::json::toJson(njh::files::normalize(bam));
		std::regex linePat(R"(([0-9]+) \+ ([0-9]+) (.*))");
		std::regex percPat(R"(\((.+) : (.+)\))");
		while(njh::files::crossPlatGetline(in, line)){
			std::smatch match;
			if(std::regex_match(line, match, linePat)){
				if(4 != match.size()){
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error in processing line : " << "\n";
					ss << line << "\n";
					ss << "Pattern should have matched 4 elements, not " << match.size() << "\n";
					throw std::runtime_error{ss.str()};
				}
				uint32_t qcPass = njh::StrToNumConverter::stoToNum<uint32_t>(match[1]);
				uint32_t qcFail = njh::StrToNumConverter::stoToNum<uint32_t>(match[2]);
				Json::Value val;
				std::string valTitle = match[3].str();
				val["QC-passed reads"] = njh::json::toJson(qcPass);
				val["QC-failed reads"] = njh::json::toJson(qcFail);;
				if(std::string::npos != valTitle.find(":")){
					std::string rest = valTitle.substr(valTitle.find(" (") + 1);
					std::smatch percMatch;
					if(std::regex_match(rest, percMatch, percPat)){
						if(3 != percMatch.size()){
							std::stringstream ss;
							ss << __PRETTY_FUNCTION__ << ", error in processing substr : " << "\n";
							ss << rest << "\n";
							ss << "Pattern should have matched 3 elements, not " << percMatch.size() << "\n";
							throw std::runtime_error{ss.str()};
						}
						if(std::string::npos != percMatch[1].str().find("%")){
							val["QC-passed reads percentage"] = njh::json::toJson(njh::StrToNumConverter::stoToNum<double>(percMatch[1].str().substr(0, percMatch[1].str().find("%"))));
						}else{
							val["QC-passed reads percentage"] = njh::json::toJson(percMatch[1].str());
						}
						if(std::string::npos != percMatch[2].str().find("%")){
							val["QC-failed reads percentage"] = njh::json::toJson(njh::StrToNumConverter::stoToNum<double>(percMatch[2].str().substr(0, percMatch[2].str().find("%"))));
						}else{
							val["QC-failed reads percentage"] = njh::json::toJson(percMatch[2].str());
						}
					}else{
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ", error in processing sustr : " << "\n";
						ss << rest << "\n";
						throw std::runtime_error{ss.str()};
					}
					valTitle = valTitle.substr(0, valTitle.find(" ("));

				}else if(std::string::npos != valTitle.find("QC-passed reads")){
					valTitle = valTitle.substr(0, valTitle.find(" ("));
				}
				outVal[valTitle] = val;
			} else {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error in processing line : " << "\n";
				ss << line << "\n";
				throw std::runtime_error{ss.str()};
			}
		}
		allValues[bam.filename().string()] = outVal;
	}

	jsonOut << allValues << std::endl;
	if(!keepIndividualFiles){
		for(const auto & fnp : samtoolsFlagFiles){
			bfs::remove(fnp);
		}
	}
	return 0;
}


int programWrapperRunner::runTrinityOnRegion(
		const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bedFile = "";
	std::string maxMemory = "10G";
	uint32_t ncpus = 4;
	double percInRegion = 0.80;
	seqSetUp setUp(inputCommands);
	setUp.pars_.directoryName_ = "";
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(percInRegion, "--percInRegion", "Percent of bases n Region");
	setUp.setOption(bedFile, "--bed", "Bed file, only first entry is used", true);
	setUp.setOption(maxMemory, "--maxMemory", "The maximum amount of memory to allow Trinity to use");
	setUp.setOption(ncpus, "--ncpus", "Number of cpus to use");
	setUp.processDefaultReader( { "--bam" }, true);
	setUp.setOption(setUp.pars_.overWriteDir_, "--overWriteDir", "If the directory already exists over write it");
	setUp.setOption(setUp.pars_.directoryName_, "--dout", "Output Directory Name");
	setUp.finishSetUp(std::cout);
	BamExtractor bExtractor(setUp.pars_.verbose_);
	njh::sys::requireExternalProgramThrow("Trinity");

	if (!setUp.pars_.ioOptions_.inExists()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << setUp.pars_.ioOptions_.firstName_
				<< " doesn't exist";
		throw std::runtime_error { ss.str() };
	}
	auto regions = gatherRegions(bedFile.string(), "", setUp.pars_.verbose_);
	std::string defaultName = bfs::basename(setUp.pars_.ioOptions_.firstName_)
			+ "_" + njh::replaceString(setUp.commands_.getProgramName(), " ", "-")
			+ "_" + regions.front().uid_ + "_" + getCurrentDate();
	if ("" != setUp.pars_.directoryName_) {
		defaultName = setUp.pars_.directoryName_;
	}
	std::string outDirName = njh::replaceString(defaultName, "TODAY",
			getCurrentDate()) + "/";
	bool needsToRun = true;
	if(bfs::exists(outDirName)){
		if(setUp.pars_.overWriteDir_){
			//calling to remove directory no matter what
			njh::files::rmDirForce(outDirName);
		}else{
			if(njh::files::firstFileIsOlder(outDirName, setUp.pars_.ioOptions_.firstName_) ||
			njh::files::firstFileIsOlder(outDirName,bedFile)){
				//out direcotry exists and the input files have changed since the run time so need to re-run again
				njh::files::rmDirForce(outDirName);
			}else{
				//out directory exists and the input files are older than the directory meaning nothing has changed no need to run
				needsToRun = false;
			}
		}
	}
	if(needsToRun){
		setUp.pars_.directoryName_ = njh::files::makeDir("./",
				njh::files::MkdirPar(defaultName, setUp.pars_.overWriteDir_)).string();

		setUp.startARunLog(setUp.pars_.directoryName_);

		OutOptions seqBaseOpts(njh::files::make_path(setUp.pars_.directoryName_,
						"trinityInput_" + bfs::basename(setUp.pars_.ioOptions_.firstName_) + "_"
								+ regions.front().uid_));

		bExtractor.writeExtractReadsFromBamRegion(setUp.pars_.ioOptions_.firstName_,
				regions.front(), percInRegion, seqBaseOpts);

		auto seqOpts = SeqIOOptions::genPairedOut(
				"trinityInput_" + bfs::basename(setUp.pars_.ioOptions_.firstName_) + "_"
						+ regions.front().uid_);

		std::stringstream trinityCmdTemplate;
		bfs::path trinityOutputDir = njh::files::make_path(
				"trinity_" + bfs::basename(setUp.pars_.ioOptions_.firstName_) + "_"
						+ regions.front().uid_);
		trinityCmdTemplate << "Trinity --seqType fq --max_memory " << maxMemory << " --left "
				<< seqOpts.getPriamryOutName() << " --right "
				<< seqOpts.getSecondaryOutName() << " --CPU " << ncpus << " --output "
				<< trinityOutputDir << " > trinityLog 2>&1";

		{
			std::ofstream outTrinityCmdFile;
			auto trinityOutOpts = OutOptions(njh::files::make_path(setUp.pars_.directoryName_, "runTrinity.sh"));
			trinityOutOpts.openExecutableFile(outTrinityCmdFile);

			outTrinityCmdFile << "#!/usr/bin/env bash" << std::endl;
			outTrinityCmdFile << trinityCmdTemplate.str() << std::endl;
		}

		auto trinityOut = njh::sys::run(VecStr{"cd " + setUp.pars_.directoryName_, "&&", "./runTrinity.sh"});

		auto trinityLogOutOpts = OutOptions(njh::files::make_path(setUp.pars_.directoryName_, "runTrinityLog.json"));
		auto outTrinityRunLog = trinityLogOutOpts.openFile();

		(*outTrinityRunLog) << trinityOut.toJson() << std::endl;

		if(trinityOut.success_){
			std::string seqName = "trinityOutput_" + bfs::basename(setUp.pars_.ioOptions_.firstName_) + "_"
										+ regions.front().uid_;
			auto trinityOpts = SeqIOOptions::genFastaInOut(njh::files::make_path(setUp.pars_.directoryName_, trinityOutputDir, "Trinity.fasta"),
					njh::files::make_path(setUp.pars_.directoryName_, seqName));
			trinityOpts.out_.overWriteFile_ = true;
			SeqIO reader(trinityOpts);
			auto trinitySeqs = reader.in_.readAllReads<seqInfo>();
			readVecSorter::sortBySeqSize(trinitySeqs, false);
			seqInfo::size_type count = 0;
			for(auto & seq : trinitySeqs){
				seq.name_ = seqName + "." + leftPadNumStr(count, len(trinitySeqs));
				++count;
			}
			std::ofstream outLenInfo;
			OutOptions outLenOpts(njh::files::make_path(setUp.pars_.directoryName_, "trinityOutput_" + bfs::basename(setUp.pars_.ioOptions_.firstName_) + "_"
					+ regions.front().uid_ + "_lenInfo.tab.txt"));
			outLenOpts.overWriteFile_ = true;
			outLenOpts.openFile(outLenInfo);
			outLenInfo << "SeqName\tlen" << '\n';
			for(const auto & seq : trinitySeqs){
				outLenInfo << seq.name_ << "\t" << len(seq) << "\n";
			}

			reader.out_.openWrite(trinitySeqs);
		}
	}
	return 0;
}

int programWrapperRunner::runSSAKEOnRegion(
		const njh::progutils::CmdArgs & inputCommands) {
	bool stitch = false;
	bfs::path bedFile = "";
	double percInRegion = 0.80;
	std::string extraSSAKEopts = "";
	seqSetUp setUp(inputCommands);
	setUp.pars_.directoryName_ = "";
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(percInRegion, "--percInRegion", "Percent of bases n Region");
	setUp.setOption(bedFile, "--bed", "Bed file, only first entry is used", true);
	setUp.processDefaultReader( { "--bam" }, true);
	setUp.setOption(setUp.pars_.overWriteDir_, "--overWriteDir", "If the directory already exists over write it");
	setUp.setOption(setUp.pars_.directoryName_, "--dout", "Output Directory Name");
	setUp.setOption(stitch, "--stitch", "Try to stitch the output pairs if there are any");
	setUp.setOption(extraSSAKEopts, "--extraSSAKEopts", "Extra options to give to SSAKE");

	setUp.finishSetUp(std::cout);
	BamExtractor bExtractor(setUp.pars_.verbose_);
	njh::sys::requireExternalProgramThrow("SSAKE");
	if (!setUp.pars_.ioOptions_.inExists()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << setUp.pars_.ioOptions_.firstName_
				<< " doesn't exist";
		throw std::runtime_error { ss.str() };
	}
	auto regions = gatherRegions(bedFile.string(), "", setUp.pars_.verbose_);
	std::string defaultName = bfs::basename(setUp.pars_.ioOptions_.firstName_)
			+ "_" + njh::replaceString(setUp.commands_.getProgramName(), " ", "-")
			+ "_" + regions.front().uid_ + "_" + getCurrentDate();
	if ("" != setUp.pars_.directoryName_) {
		defaultName = setUp.pars_.directoryName_;
	}
	std::string outDirName = njh::replaceString(defaultName, "TODAY",
			getCurrentDate()) + "/";
	if(bfs::exists(outDirName)){
		if(setUp.pars_.overWriteDir_){
			//calling to remove directory no matter what
			njh::files::rmDirForce(outDirName);
		}else{
			if(njh::files::firstFileIsOlder(outDirName, setUp.pars_.ioOptions_.firstName_) ||
			njh::files::firstFileIsOlder(outDirName,bedFile)){
				//out direcotry exists and the input files have changed since the run time so need to re-run again
				njh::files::rmDirForce(outDirName);
			}else{
				//out directory exists and the input files are older than the directory meaning nothing has changed no need to run
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << outDirName
						<< " already exists, use --overWriteDir to over write it"  << "\n";
				throw std::runtime_error { ss.str() };
			}
		}
	}

	setUp.pars_.directoryName_ = njh::files::makeDir("./",
			njh::files::MkdirPar(defaultName, setUp.pars_.overWriteDir_)).string();

	setUp.startARunLog(setUp.pars_.directoryName_);

	auto ssakeOutSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(setUp.pars_.directoryName_,
			"ssakeInput_" + bfs::basename(setUp.pars_.ioOptions_.firstName_) + "_"
					+ regions.front().uid_));
	SeqOutput writer(ssakeOutSeqOpts);
	writer.openOut();

	auto extractedBase = njh::files::make_path(setUp.pars_.directoryName_,
			"extracted_" + bfs::basename(setUp.pars_.ioOptions_.firstName_) + "_"
					+ regions.front().uid_);
	OutOptions seqBaseOpts(extractedBase);

	if(stitch){
		auto opts = bExtractor.writeExtractReadsFromBamRegionStitch(setUp.pars_.ioOptions_.firstName_,
				regions.front(),percInRegion, seqBaseOpts);
		//singles
		if (opts.inUnpaired_.inExists()) {
			seqInfo seq;
			SeqInput reader(opts.inUnpaired_);
			reader.openIn();
			while (reader.readNextRead(seq)) {
				writer.write(seq);
			}
		}
		//stitched
		if (opts.stitched_.inExists()) {
			seqInfo seq;
			SeqInput reader(opts.stitched_);
			reader.openIn();
			while (reader.readNextRead(seq)) {
				writer.write(seq);
			}
		}
		//un-stitched
		if (opts.notCombinedPairs_.inExists()) {
			PairedRead seq;
			SeqInput reader(opts.notCombinedPairs_);
			reader.openIn();
			while (reader.readNextRead(seq)) {
				writer.write(seq.seqBase_);
				writer.write(seq.mateSeqBase_);
			}
		}
	}else{
		bExtractor.writeExtractReadsFromBamRegion(setUp.pars_.ioOptions_.firstName_,
				regions.front(),percInRegion, seqBaseOpts);
		auto pairedSeqOpts = SeqIOOptions::genPairedOut(extractedBase);
		auto unpaired = SeqIOOptions::genFastqOut(extractedBase);
		if(unpaired.outExists()){
			auto inOpts = SeqIOOptions::genFastqIn(unpaired.getPriamryOutName());
			seqInfo seq;
			SeqInput reader(inOpts);
			reader.openIn();
			while(reader.readNextRead(seq)){
				writer.write(seq);
			}
		}
		if(pairedSeqOpts.outExists()){
			auto inOpts = SeqIOOptions::genPairedIn(pairedSeqOpts.getPriamryOutName(), pairedSeqOpts.getSecondaryOutName());
			PairedRead seq;
			SeqInput reader(inOpts);
			reader.openIn();
			while(reader.readNextRead(seq)){
				writer.write(seq.seqBase_);
				writer.write(seq.mateSeqBase_);
			}
		}
	}


	writer.closeOut();

	std::stringstream ssakeCmdTemplate;
	ssakeCmdTemplate << "SSAKE -w 1 -f " << ssakeOutSeqOpts.getPriamryOutName().filename()  << " -b ssakeOutput " << extraSSAKEopts << " >  SSAKELog 2>&1";

	{
		std::ofstream outSSAKECmdFile;
		auto ssakeOutOpts = OutOptions(njh::files::make_path(setUp.pars_.directoryName_, "runSSAKE.sh"));
		ssakeOutOpts.openExecutableFile(outSSAKECmdFile);

		outSSAKECmdFile << "#!/usr/bin/env bash" << std::endl;
		outSSAKECmdFile << ssakeCmdTemplate.str() << std::endl;
	}

	auto SSAKEOut = njh::sys::run(VecStr{"cd " + setUp.pars_.directoryName_, "&&", "./runSSAKE.sh"});

	auto ssakeLogOutOpts = OutOptions(njh::files::make_path(setUp.pars_.directoryName_, "runSSAKELog.json"));
	auto outSSAKERunLog = ssakeLogOutOpts.openFile();

	(*outSSAKERunLog) << ssakeLogOutOpts.toJson() << std::endl;

	if(SSAKEOut.success_){
		auto ssakeContigsOpts = SeqIOOptions::genFastaInOut(njh::files::make_path(setUp.pars_.directoryName_, "ssakeOutput.contigs"),
				njh::files::make_path(setUp.pars_.directoryName_, "output.fasta"));
		SeqInput reader(ssakeContigsOpts);
		auto seqs = reader.readAllReads<seqInfo>();
		std::regex ssakeNamePat{"([a-z]+)([-+]?[0-9]*\\.?[0-9]+)"};
		double total = 0;
		for(auto & seq : seqs){
			auto toks = njh::tokenizeString(seq.name_, "|");
			MetaDataInName meta;
			for(const auto & tok : toks){
				std::smatch match;
				if(std::regex_match(tok, match, ssakeNamePat)){
					meta.addMeta(match[1], match[2]);
				}else{
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error name should contain a pattern of letters followed by numbers not " << tok << "\n";
					throw std::runtime_error{ss.str()};
				}
			}
			if(!meta.containsMeta("contig") || !meta.containsMeta("read")){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << "Error" << "\n";
				if(!meta.containsMeta("contig")){
					ss << "couldn't find field contig in " << seq.name_ << "\n";
				}
				if(!meta.containsMeta("read")){
					ss << "couldn't find field read in " << seq.name_ << "\n";
				}
				throw std::runtime_error{ss.str()};
			}
			seq.name_ = "contig." + leftPadNumStr(meta.getMeta<size_t>("contig"), seqs.size());
			seq.cnt_ = meta.getMeta<double>("read");
			total+= seq.cnt_;
			meta.resetMetaInName(seq.name_);
			seq.updateName();

		}
		for(auto & seq : seqs){
			MetaDataInName meta(seq.name_);
			seq.setFractionByCount(total);
			meta.addMeta("frac", seq.frac_);
			meta.resetMetaInName(seq.name_);
		}
		SeqOutput::write(seqs, ssakeContigsOpts);
		auto infoTabOpts = TableIOOpts::genTabFileOut(njh::files::make_path(setUp.pars_.directoryName_, "outputInfo.tab.txt"));
		table infoTab(VecStr{"seqName", "reads", "fraction", "len"});
		for(const auto & seq : seqs){
			infoTab.addRow(seq.name_, seq.cnt_, seq.frac_, len(seq));
		}
		infoTab.outPutContents(infoTabOpts);
		/*
		std::string seqName = "trinityOutput_" + bfs::basename(setUp.pars_.ioOptions_.firstName_) + "_"
									+ regions.front().uid_;
		auto trinityOpts = SeqIOOptions::genFastaInOut(njh::files::make_path(setUp.pars_.directoryName_, trinityOutputDir, "Trinity.fasta"),
				njh::files::make_path(setUp.pars_.directoryName_, seqName));
		trinityOpts.out_.overWriteFile_ = true;
		SeqIO reader(trinityOpts);
		auto trinitySeqs = reader.in_.readAllReads<seqInfo>();
		readVecSorter::sortBySeqSize(trinitySeqs, false);
		seqInfo::size_type count = 0;
		for(auto & seq : trinitySeqs){
			seq.name_ = seqName + "." + leftPadNumStr(count, len(trinitySeqs));
			++count;
		}
		std::ofstream outLenInfo;
		OutOptions outLenOpts(njh::files::make_path(setUp.pars_.directoryName_, "trinityOutput_" + bfs::basename(setUp.pars_.ioOptions_.firstName_) + "_"
				+ regions.front().uid_ + "_lenInfo.tab.txt"));
		outLenOpts.overWriteFile_ = true;
		outLenOpts.openFile(outLenInfo);
		outLenInfo << "SeqName\tlen" << '\n';
		for(const auto & seq : trinitySeqs){
			outLenInfo << seq.name_ << "\t" << len(seq) << "\n";
		}

		reader.out_.openWrite(trinitySeqs);
		*/
	}
	return 0;
}


int programWrapperRunner::runLastz(const njh::progutils::CmdArgs & inputCommands){
	double coverage = 90;
	double identity = 50;
	bfs::path genomeFnp = "";
	std::string outFormat = "SAM";
	std::string extraLastzArgs = "";
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(outFormat, "--lastzOutFormat", "lastz Out Format");
	setUp.setOption(extraLastzArgs, "--extraLastzArgs", "extra lastz Args");
	setUp.processDefaultReader(VecStr{"--fasta"});
	setUp.pars_.ioOptions_.out_.outExtention_ = ".sorted.bam";
	setUp.setOption(genomeFnp, "--genomeFnp", "Filename path to genome", true);
	setUp.setOption(coverage, "--coverage", "Coverage");
	setUp.setOption(identity, "--identity", "Identity");

	setUp.finishSetUp(std::cout);

	if(setUp.pars_.ioOptions_.out_.outExists() && !setUp.pars_.ioOptions_.out_.overWriteFile_){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error, file "
				<< setUp.pars_.ioOptions_.out_.outName()
				<< " already exists use --overWrite to over it" << "\n";
		throw std::runtime_error{ss.str() };
	}
	njh::sys::requireExternalProgramThrow("lastz");
	njh::sys::requireExternalProgramThrow("samtools");
	//determineRegionLastz
	std::stringstream lastzCmd;
	lastzCmd << "lastz " << genomeFnp << "[multiple] "
			<< setUp.pars_.ioOptions_.firstName_ << " --format=" << outFormat
			<< " --coverage=" << coverage << " --identity=" << identity << " --ambiguous=iupac "
			<< extraLastzArgs << " | samtools view - -b | samtools sort - -o "
			<< setUp.pars_.ioOptions_.out_.outName() << " " << "&& samtools index "
			<< setUp.pars_.ioOptions_.out_.outName();

	auto runOutput = njh::sys::run(VecStr{lastzCmd.str()});

	if(!runOutput.success_ || setUp.pars_.debug_){
		std::cout << runOutput.toJson() << std::endl;
	}

	return 0;
}


int programWrapperRunner::sraIsPairedEnd(const njh::progutils::CmdArgs & inputCommands){
	std::string fastqDumpCmd = "fastq-dump";
	bfs::path sraFnp = "";
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(fastqDumpCmd, "--fastqDumpCmd", "fastq Dump Cmd");
	setUp.setOption(sraFnp, "--sraFnp", "Filename path to SRA file", true);
	setUp.finishSetUp(std::cout);

	if(!bfs::exists(sraFnp)){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << sraFnp << " doesn't exist" << "\n";
		throw std::runtime_error{ss.str()};
	}

	std::stringstream cmd;
	cmd <<  fastqDumpCmd << " --log-level 0 -X 1 -Z --split-spot " << sraFnp;
	auto cmdOutput = njh::sys::run({cmd.str()});
	BioCmdsUtils::checkRunOutThrow(cmdOutput, __PRETTY_FUNCTION__);
	//njh::sys::run trim end white space so have to add one;
	uint32_t newLines = countOccurences(cmdOutput.stdOut_, "\n") + 1;
	if(4 == newLines){
		std::cout << "false" << std::endl;
	}else if(8 == newLines || 12 == newLines){
		//check for 8 or 12, 12 means triple file _1 forward mate, _2 barcode, _3 reverse mate
		std::cout << "true" << std::endl;
	}else{
		std::cout << "indeterminate" << std::endl;
		if(setUp.pars_.verbose_){
			std::cout << "indeterminate, had " << newLines << ", was expecting 4, 8 or 12" << std::endl;
		}
	}

	return 0;
}

int programWrapperRunner::sraFastqDump(const njh::progutils::CmdArgs & inputCommands){
	std::string fastqDumpCmd = "fastq-dump";
	bfs::path sraFnp = "";
	BioCmdsUtils::FastqDumpPars sraPars;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(fastqDumpCmd, "--fastqDumpCmd", "fastq Dump Cmd");
	setUp.setOption(sraPars.outputDir_, "--outputDir", "Output directory");
	setUp.setOption(sraPars.exportBarCode_, "--exportBarCode", "export Bar Code file if present");
	setUp.setOption(sraPars.gzip_, "--gzip", "Make output gzip compressed");
	setUp.setOption(sraPars.force_, "--force", "Even if the output files exist and appear to be up to date, extract again");
	setUp.setOption(sraPars.sraFnp_, "--sraFnp", "Filename path to SRA file", true);
	setUp.finishSetUp(std::cout);

	BioCmdsUtils bioRunner(setUp.pars_.verbose_);
	auto results = bioRunner.runFastqDump(sraPars);
	if (setUp.pars_.verbose_) {
		std::cout << njh::json::toJson(results) << std::endl;
	}
	return 0;
}

int programWrapperRunner::runTrimmomatic(const njh::progutils::CmdArgs & inputCommands){
	std::string trimmomaticCmd = "trimmomatic";
	bfs::path adaptersDir = "";
	std::string adaptor = "";
	bfs::path outputDir = "";
	bfs::path fastq1 = "";
	bfs::path fastq2 = "";
	bfs::path fastqStub = "";
	bool force = false;
	uint32_t numThreads = 1;
	bool createTrimmomaticLog = false;
	bool keepBoth = false;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	bool stubGiven = setUp.setOption(fastqStub, "--fastqStub", "fastq stub");
	setUp.setOption(force, "--force", "force run");
	setUp.setOption(createTrimmomaticLog, "--createTrimmomaticLog", "create trimmomatic log file");
	setUp.setOption(numThreads, "--numThreads", "Number of threads");
	auto adaptorGiven = setUp.setOption(adaptor, "--adaptor", "adaptor");
	setUp.setOption(adaptersDir, "--adaptorsDir", "adapters directory", adaptorGiven);
	setUp.setOption(outputDir, "--outputDir", "output directory");
	setUp.setOption(trimmomaticCmd, "--trimmomaticCmd", "trimmomatic Cmd");
	setUp.setOption(fastq1, "--fastq1", "Filename path to fastq1", !stubGiven);
	setUp.setOption(fastq2, "--fastq2", "Filename path to fastq2");
	setUp.setOption(keepBoth, "--keepBoth", "keep both reads when there is read through");
	setUp.finishSetUp(std::cout);
	if(stubGiven){
		fastq1 = fastqStub.string() + "_R1.fastq";
		fastq2 = fastqStub.string() + "_R2.fastq";
		if(!bfs::exists(fastq1) && !bfs::exists(fastq2)){
			fastq1 = fastqStub.string() + "_R1.fastq.gz";
			fastq2 = fastqStub.string() + "_R2.fastq.gz";
			if(!bfs::exists(fastq1) && !bfs::exists(fastq2)){
				fastq1 = fastqStub.string() + "_1.fastq";
				fastq2 = fastqStub.string() + "_2.fastq";
			}
			if(!bfs::exists(fastq1) && !bfs::exists(fastq2)){
				fastq1 = fastqStub.string() + "_1.fastq.gz";
				fastq2 = fastqStub.string() + "_2.fastq.gz";
			}
		}
	}
	bool pairedEnd = false;
	if ("" != fastq2) {
		pairedEnd = true;
	}
	njh::files::checkExistenceThrow(fastq1, __PRETTY_FUNCTION__);
	if (pairedEnd) {
		njh::files::checkExistenceThrow(fastq2, __PRETTY_FUNCTION__);
	}
	if ("" != adaptersDir) {
		njh::files::checkExistenceThrow(adaptersDir, __PRETTY_FUNCTION__);
		njh::files::checkExistenceThrow(njh::files::make_path(adaptersDir, adaptor),
				__PRETTY_FUNCTION__);
	}
	auto getFqBaseName = [](const bfs::path & fnp){
		auto bName = bfs::basename(fnp);
		if(njh::endsWith(bName, ".fq")){
			bName = bName.substr(0, bName.rfind(".fq"));
		}else if(njh::endsWith(bName, ".fastq")){
			bName = bName.substr(0, bName.rfind(".fastq"));
		}
		return bName;
	};
	std::stringstream cmd;
	bool needToRun = true;
	if(pairedEnd){


		auto outFirstPrimary = njh::files::make_path(outputDir, "trimmed_" + getFqBaseName(fastq1) + ".fastq");
		if(njh::endsWith(fastq1.string(), ".gz")){
			outFirstPrimary = outFirstPrimary.string() + ".gz";
		}
		auto outFirstUnpaired = njh::files::make_path(outputDir, "trimmed_" + getFqBaseName(fastq1) + "_unpaired.fastq");
		if(njh::endsWith(fastq1.string(), ".gz")){
			outFirstUnpaired = outFirstUnpaired.string() + ".gz";
		}
		auto outSecondPrimary = njh::files::make_path(outputDir, "trimmed_" + getFqBaseName(fastq2) + ".fastq");
		if(njh::endsWith(fastq1.string(), ".gz")){
			outSecondPrimary = outSecondPrimary.string() + ".gz";
		}
		auto outSecondUnpaired = njh::files::make_path(outputDir, "trimmed_" + getFqBaseName(fastq2) + "_unpaired.fastq");
		if(njh::endsWith(fastq1.string(), ".gz")){
			outSecondUnpaired = outSecondUnpaired.string() + ".gz";
		}

		cmd << trimmomaticCmd << " PE "
				<< " -threads " << numThreads;
		if(createTrimmomaticLog){
			cmd << " -trimlog " << njh::files::make_path(outputDir, "trimmomatic_" + getFqBaseName(fastq1) + "_log.txt");
		}
		cmd << " " << fastq1
				<< " " << fastq2
				<< " " << outFirstPrimary
				<< " " << outFirstUnpaired
				<< " " << outSecondPrimary
				<< " " << outSecondUnpaired;
		if ("" != adaptor) {
			cmd << " ILLUMINACLIP:"<< njh::files::make_path(adaptersDir, adaptor)<<":2:30:10";
			if(keepBoth){
				cmd << ":TRUE";
			}
		}
		cmd	<< " " << "LEADING:3"
				<< " " << "TRAILING:3"
				<< " " << "SLIDINGWINDOW:4:15"
				<< " " << "MINLEN:36";
		if(bfs::exists(outFirstPrimary) &&
				njh::files::firstFileIsOlder(fastq1, outFirstPrimary) &&
				njh::files::firstFileIsOlder(fastq1, outFirstUnpaired) &&
				bfs::exists(outSecondPrimary) &&
								njh::files::firstFileIsOlder(fastq2, outSecondPrimary) &&
								njh::files::firstFileIsOlder(fastq2, outSecondUnpaired)){
			needToRun = false;
		}
	} else {
		auto outFirstPrimary = njh::files::make_path(outputDir, "trimmed_" + getFqBaseName(fastq1) + ".fastq");
		if(njh::endsWith(fastq1.string(), ".gz")){
			outFirstPrimary = outFirstPrimary.string() + ".gz";
		}
		cmd << trimmomaticCmd << " SE "
				<< " -threads " << numThreads;
		if(createTrimmomaticLog){
			cmd << " -trimlog " << njh::files::make_path(outputDir,"trimmomatic_" + getFqBaseName(fastq1) + "_log.txt");
		}
		cmd << " " << fastq1
				<< " " << outFirstPrimary;
		if ("" != adaptor) {
			cmd << " ILLUMINACLIP:"<< njh::files::make_path(adaptersDir, adaptor)<<":2:30:10";
		}
		cmd		  << " " << "LEADING:3"
						<< " " << "TRAILING:3"
						<< " " << "SLIDINGWINDOW:4:15"
						<< " " << "MINLEN:36";
		if(bfs::exists(outFirstPrimary) &&
				njh::files::firstFileIsOlder(fastq1, outFirstPrimary)){
			needToRun = false;
		}
	}
	if(force){
		needToRun = true;
	}
	if(needToRun){
		auto cmdOutput = njh::sys::run({cmd.str()});
		BioCmdsUtils::checkRunOutThrow(cmdOutput, __PRETTY_FUNCTION__);
		if(setUp.pars_.verbose_){
			std::cout << cmdOutput.toJson() << std::endl;
		}
	}
	return 0;
}


int programWrapperRunner::runBwaOnTrimmomaticOutputPE(const njh::progutils::CmdArgs & inputCommands){
	bfs::path outputDir = "";
	bfs::path trimStub = "";
	bool force = false;
	uint32_t numThreads = 1;
	std::string sampName = "";
	bool removeIntermediateFiles = false;
	bfs::path genomeFnp = "";
	std::string extraBwaArgs = "";
	bfs::path outputFnp = "";
	bool leaveOutUnpaired = false;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(genomeFnp, "--genomeFnp", "genome file to align to", true);
	setUp.setOption(trimStub, "--trimStub", "Trimmomatic output stub", true);
	setUp.setOption(sampName, "--sampName", "Sample Name to give to final bam", true);
	setUp.setOption(force, "--force", "force run");
	setUp.setOption(leaveOutUnpaired, "--leaveOutUnpaired", "leave out unpaired sequences");
	setUp.setOption(outputFnp, "--outputFnp", "output name, will default to sampName.sorted.bam");
	setUp.setOption(extraBwaArgs, "--extraBwaArgs", "extra bwa arguments");
	setUp.setOption(numThreads, "--numThreads", "Number of threads");
	setUp.setOption(outputDir, "--outputDir", "output directory");
	setUp.setOption(removeIntermediateFiles, "--removeIntermediateFiles", "remove the intermediate bam files and just keep the ");
	setUp.finishSetUp(std::cout);
	BioCmdsUtils bioRunner(setUp.pars_.verbose_);
	njh::sys::requireExternalProgramsThrow(VecStr{"bwa", "bamtools", "samtools"});
	bfs::path inputUnpairedFirstMates = njh::files::make_path(trimStub.string() + "_1_unpaired.fastq");
	bfs::path inputUnpairedSecondMates = njh::files::make_path(trimStub.string() + "_2_unpaired.fastq");
	bfs::path inputPairedFirstMates = njh::files::make_path(trimStub.string() + "_1.fastq");
	bfs::path inputPairedSecondMates = njh::files::make_path(trimStub.string() + "_2.fastq");
	njh::files::checkExistenceThrow(genomeFnp,__PRETTY_FUNCTION__);
	bioRunner.RunBwaIndex(genomeFnp);



	if (!bfs::exists(inputPairedFirstMates)
			&& !bfs::exists(inputPairedSecondMates)) {
		inputPairedFirstMates = njh::files::make_path(trimStub.string() + "_R1.fastq");
		inputPairedSecondMates = njh::files::make_path(trimStub.string() + "_R2.fastq");
		inputUnpairedFirstMates = njh::files::make_path(trimStub.string() + "_R1_unpaired.fastq");
		inputUnpairedSecondMates = njh::files::make_path(trimStub.string() + "_R2_unpaired.fastq");
		if(!bfs::exists(inputPairedFirstMates) &&
				!bfs::exists(inputPairedSecondMates)){
			inputPairedFirstMates = njh::files::make_path(trimStub.string() + "_1.fastq.gz");
			inputPairedSecondMates = njh::files::make_path(trimStub.string() + "_2.fastq.gz");
			inputUnpairedFirstMates = njh::files::make_path(trimStub.string() + "_1_unpaired.fastq.gz");
			inputUnpairedSecondMates = njh::files::make_path(trimStub.string() + "_2_unpaired.fastq.gz");

			if (!bfs::exists(inputPairedFirstMates)
					&& !bfs::exists(inputPairedSecondMates)) {
				inputPairedFirstMates = njh::files::make_path(trimStub.string() + "_R1.fastq.gz");
				inputPairedSecondMates = njh::files::make_path(trimStub.string() + "_R2.fastq.gz");
				inputUnpairedFirstMates = njh::files::make_path(trimStub.string() + "_R1_unpaired.fastq.gz");
				inputUnpairedSecondMates = njh::files::make_path(trimStub.string() + "_R2_unpaired.fastq.gz");
				if (!bfs::exists(inputPairedFirstMates)
						&& !bfs::exists(inputPairedSecondMates)) {
					inputPairedFirstMates = njh::files::make_path(trimStub.string() + "_R1_001.fastq.gz");
					inputPairedSecondMates = njh::files::make_path(trimStub.string() + "_R2_001.fastq.gz");
					inputUnpairedFirstMates = njh::files::make_path(trimStub.string() + "_R1_001_unpaired.fastq.gz");
					inputUnpairedSecondMates = njh::files::make_path(trimStub.string() + "_R2_001_unpaired.fastq.gz");
					if (!bfs::exists(inputPairedFirstMates)
							&& !bfs::exists(inputPairedSecondMates)) {
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ", error, one of the following pairs  " << "\n"
								<< " " << njh::files::make_path(trimStub.string() + "_1.fastq") << " "
								<< njh::files::make_path(trimStub.string() + "_2.fastq") << "\n"
								<< " or " << "\n"
								<< " " << njh::files::make_path(trimStub.string() + "_1.fastq.gz")
								<< " " << njh::files::make_path(trimStub.string() + "_2.fastq.gz")
								<< " or " << "\n"
								<< " " << njh::files::make_path(trimStub.string() + "_R1.fastq")
								<< " " << njh::files::make_path(trimStub.string() + "_R2.fastq")
								<< " or " << "\n"
								<< " " << njh::files::make_path(trimStub.string() + "_R1.fastq.gz")
								<< " " << njh::files::make_path(trimStub.string() + "_R2.fastq.gz")
								<< "\n" << "have to exist" << "\n";
						throw std::runtime_error { ss.str() };
					}
				}
			}
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
	bfs::path unpairedFirstMatesSortedBam = njh::files::make_path(outputDir, trimStub.filename().string() + "_1_unpaired.sorted.bam");
	bfs::path unpairedSecondMatesSortedBam = njh::files::make_path(outputDir,trimStub.filename().string() + "_2_unpaired.sorted.bam");
	bfs::path pairedSortedBam = njh::files::make_path(outputDir,trimStub.filename().string() + ".sorted.bam");


	std::string bNameStub = trimStub.filename().string();
	std::stringstream firstMateUnpairedCmd;
	firstMateUnpairedCmd << "bwa mem  -M -t " << numThreads
			<< " -R " << R"("@RG\tID:)" << bNameStub << "_1_unpaired" << R"(\tSM:)"
			<< sampName << R"(")"
			<< " "   << extraBwaArgs
			<< " "   << genomeFnp
			<< " "   << inputUnpairedFirstMates
			<< " 2> " << bfs::path(unpairedFirstMatesSortedBam.string() + ".bwa.log")
			<< " | samtools sort -@ " << numThreads << " -o " << unpairedFirstMatesSortedBam;

	std::stringstream secondMateUnpairedCmd;
	secondMateUnpairedCmd << "bwa mem  -M -t " << numThreads
			<< " -R " << R"("@RG\tID:)" << bNameStub << "_2_unpaired" << R"(\tSM:)"
			<< sampName << R"(")"
			<< " " << extraBwaArgs
			<< " " << genomeFnp
			<< " " << inputUnpairedSecondMates
			<< " 2> " << bfs::path(unpairedSecondMatesSortedBam.string() + ".bwa.log")
			<< " | samtools sort -@ " << numThreads << " -o " << unpairedSecondMatesSortedBam;

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
	if(bfs::exists(inputUnpairedFirstMates)){
		bamtoolsMergeAndIndexCmd << " -in " <<  unpairedFirstMatesSortedBam;
	}
	if(bfs::exists(inputUnpairedSecondMates)){
		bamtoolsMergeAndIndexCmd << " -in " <<  unpairedSecondMatesSortedBam;
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
		if(bfs::exists(inputUnpairedFirstMates) && !leaveOutUnpaired){
			auto unpairedFirstMatesRunOutput = njh::sys::run({firstMateUnpairedCmd.str()});
			BioCmdsUtils::checkRunOutThrow(unpairedFirstMatesRunOutput, __PRETTY_FUNCTION__);
			runOutputs["bwa-unpairedFirstMates"] = unpairedFirstMatesRunOutput;
		}

		if(bfs::exists(inputUnpairedSecondMates) && !leaveOutUnpaired){
			auto unpairedSecondMatesRunOutput = njh::sys::run({secondMateUnpairedCmd.str()});
			BioCmdsUtils::checkRunOutThrow(unpairedSecondMatesRunOutput, __PRETTY_FUNCTION__);
			runOutputs["bwa-unpairedSecondMates"] = unpairedSecondMatesRunOutput;
		}

		auto pairedRunOutput = njh::sys::run({pairedCmd.str()});
		BioCmdsUtils::checkRunOutThrow(pairedRunOutput, __PRETTY_FUNCTION__);
		runOutputs["bwa-paired"] = pairedRunOutput;
		if(bfs::exists(unpairedFirstMatesSortedBam) || bfs::exists(unpairedSecondMatesSortedBam)){
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
			if(bfs::exists(unpairedFirstMatesSortedBam) || bfs::exists(unpairedSecondMatesSortedBam)){
				bfs::remove(pairedSortedBam);
				if(bfs::exists(unpairedFirstMatesSortedBam)){
					bfs::remove(unpairedFirstMatesSortedBam);
				}
				if(bfs::exists(unpairedSecondMatesSortedBam)){
					bfs::remove(unpairedSecondMatesSortedBam);
				}
			}
		}
	}
	return 0;
}


} // namespace njhseq
