/*
 * programWrappers_runPicardMarkDups.cpp
 *
 *  Created on: May 27, 2019
 *      Author: nicholashathaway
 */




#include "programWrappers.hpp"
#include "elucidator/objects/BioDataObject.h"
#include "elucidator/BamToolsUtils.h"



namespace njhseq {


int programWrapperRunner::runPicardMarkDups(const njh::progutils::CmdArgs & inputCommands){
	bfs::path bamFnp = "";
	bool overWrite = false;
	std::string extraArgs = "";
	bfs::path outDir = "./";
	std::string outNameStub = "";
	seqSetUp setUp(inputCommands);
	setUp.description_ = "Wrapper for running picard MarkDuplicates on a coordinate sorted bam file";
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(bamFnp, "--bamFnp", "Bam file",true);
	setUp.setOption(extraArgs, "--extraArgs", "Extra arguments to give to picard");
	setUp.setOption(overWrite, "--overWrite", "over write existing files");
	setUp.setOption(outDir, "--outDir", "output Directory");
	bfs::path logDir = outDir;
	setUp.setOption(logDir, "--logDir", "logs Directory");
	setUp.setOption(outNameStub, "--outNameStub", "Out Name Stub");

	setUp.finishSetUp(std::cout);
	setUp.startARunLog(logDir.string());
	njh::sys::requireExternalProgramThrow("picard");
	if("" == outNameStub){
		outNameStub = bamFnp.filename().string();
		if(njh::endsWith(outNameStub, ".sorted.bam")){
			outNameStub = outNameStub.substr(0, outNameStub.rfind(".sorted.bam"));
		}else if(njh::endsWith(outNameStub, ".bam")){
			outNameStub = outNameStub.substr(0, outNameStub.rfind(".bam"));
		}
	}
	bfs::path outBamFnp = njh::files::make_path(outDir,
			outNameStub + ".mdups.sorted.bam");
	bfs::path outBamBaiFnp = njh::files::make_path(outDir,
			outNameStub + ".mdups.sorted.bam.bai");
	bfs::path outBamBaiPicardFnp = njh::files::make_path(outDir,
			outNameStub + ".mdups.sorted.bai");
	bfs::path outMetricsFnp = njh::files::make_path(logDir,
			outNameStub + ".mdups_MarkDuplicates.metrics.txt");
	;
	bfs::path outRunLog = njh::files::make_path(logDir,
			outNameStub + "_runPicardMarkDups_" + njh::getCurrentDate() + ".txt");
	if (!overWrite) {
		VecStr fnpExists;
		auto checkExists = [&fnpExists](const bfs::path & fnp) {
			if(bfs::exists(fnp)) {
				fnpExists.emplace_back(fnp.string());
			}
		};
		checkExists(outBamFnp);
		checkExists(outBamBaiFnp);
		checkExists(outMetricsFnp);
		checkExists(outRunLog);
		if (!fnpExists.empty()) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error "
					<< "the following files exist, use --overWrite to over write them: "
					<< njh::conToStr(fnpExists, ", ") << "\n";
			throw std::runtime_error { ss.str() };
		}
	}

	std::stringstream cmd;
	cmd << "picard MarkDuplicates "
			<< " I=" << bamFnp
			<< " O=" << outBamFnp
			<< " M=" << outMetricsFnp
			<< " CREATE_INDEX=true ";
	if("" != extraArgs){
		cmd << " " << extraArgs;
	}
	cmd << " > " << outRunLog << " 2>&1 ";

	auto runCmdLog = njh::sys::run(VecStr{cmd.str()});
	setUp.rLog_ << "\n" << runCmdLog.toJson() << "\n";

	if(runCmdLog.success_){
		if(bfs::exists(outBamBaiPicardFnp)){
			bfs::rename(outBamBaiPicardFnp, outBamBaiFnp);
		}
	}
	return 0;

}

}// namespace njhseq

