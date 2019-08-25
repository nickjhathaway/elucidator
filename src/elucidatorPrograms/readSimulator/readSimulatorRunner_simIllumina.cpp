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


}  // namespace njhseq
