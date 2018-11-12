/*
 * genExp_splitLibrary.cpp
 *
 *  Created on: May 21, 2018
 *      Author: nick
 */

#include "genExp.hpp"


namespace njhseq {

int genExpRunner::splitLibrary(const njh::progutils::CmdArgs & inputCommands){
	bool nocompress = false;
	uint32_t numThreads = 1;
	seqSetUp setUp(inputCommands);
	setUp.setOption(nocompress,   "--noCompression", "Don't compress the output files afterwards using gzip");
	setUp.setOption(numThreads, "--numThreads", "Number of CPUs allowed to be utilized");
	setUp.processDefaultReader(VecStr{"--fastq1", "--fastq2", "--fastq1gz", "--fastq2gz"}, true);
	setUp.processDirectoryOutputName("split_library");
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);

	MultiSeqIO outs;

	PairedRead seq;
	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();

	while(reader.readNextRead(seq)){
		auto toks = tokenizeString(seq.seqBase_.name_, "whitespace");
		if(toks.size() < 2){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error in processing, " << seq.seqBase_.name_ << " was expecting at least 2 sections separated by whitespace" << "\n";
			throw std::runtime_error{ss.str()};
		}
		auto subToks = tokenizeString(toks[1], ":");
		if(4 > subToks.size()){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error in processing, " << toks[1] << " was expecting at least 4 sections separated by :" << "\n";
			throw std::runtime_error{ss.str()};
		}
		if(!outs.containsReader(subToks[3])){
			auto outOpts = setUp.pars_.ioOptions_;
			outOpts.out_.outFilename_ = njh::files::make_path(setUp.pars_.directoryName_, subToks[3]).string();
			outs.addReader(subToks[3], outOpts);
		}
		outs.openWrite(subToks[3], seq);
	}
	outs.closeOutAll();

	if(!nocompress){
		VecStr cmds;
		auto outFiles = njh::files::listAllFiles(setUp.pars_.directoryName_, false, {std::regex(R"(.*\.fastq$)")});
		for(const auto & f : outFiles){
			cmds.emplace_back("gzip " + f.first.string());
		}
		auto outs = njh::sys::runCmdsThreaded(cmds, numThreads, setUp.pars_.verbose_, false);
		bool fail = false;
		std::stringstream ss;
		for(const auto & out : outs){
			if(!out.success_){
				fail = true;
				ss << __PRETTY_FUNCTION__ << " Error in running " << out.cmd_ <<"\n";
			}
		}
		if(fail){
			throw std::runtime_error{ss.str()};
		}
	}

	return 0;
}



}  // namespace njhseq
