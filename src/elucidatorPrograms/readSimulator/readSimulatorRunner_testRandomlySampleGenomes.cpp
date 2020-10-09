/*
 * readSimulatorRunner_testRandomlySampleGenomes.cpp
 *
 *  Created on: Oct 7, 2019
 *      Author: nicholashathaway
 */





#include "readSimulatorRunner.hpp"
#include "elucidator/simulation.h"
#include <SeekDeep/objects/TarAmpSetupUtils/PrimersAndMids.hpp>

namespace njhseq {



int readSimulatorRunner::testSDFinalReadAmount(
		const njh::progutils::CmdArgs & inputCommands) {

	uint64_t finalReadAmount = 5000;
	double sdFrac = 0.389;
	uint32_t simAmount = 10;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(simAmount, "--simAmount", "Number of sim to conduct");
	setUp.setOption(sdFrac, "--sdFrac", "sdFrac", true, njh::progutils::ProgramSetUp::CheckCase::GREATERZERO);
	setUp.setOption(finalReadAmount, "--finalReadAmount", "finalReadAmount", true, njh::progutils::ProgramSetUp::CheckCase::GREATERZERO);

	setUp.finishSetUp(std::cout);
  std::random_device rd;
  // Mersenne twister PRNG, initialized with seed from previous random device instance
  std::mt19937 gen(rd());

	std::normal_distribution<double> ndist(finalReadAmount, sdFrac * finalReadAmount);
	for(uint32_t sim = 0; sim < simAmount; ++sim){
		std::cout << std::round(ndist(gen)) << std::endl;
	}
	return 0;
}



int readSimulatorRunner::testRandomlySampleGenomes(
		const njh::progutils::CmdArgs & inputCommands) {

	uint64_t totalGenomes = 6000;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processReadInNames(true);
	setUp.setOption(totalGenomes, "--totalGenomes", "Total Starting Genomes", true);
	setUp.finishSetUp(std::cout);

	auto seqs = SeqInput::getSeqVec<seqInfo>(setUp.pars_.ioOptions_);
	auto results = PCRSimulator::randomlySampleGenomes(seqs, totalGenomes);
	for(const auto & res : results){
		std::cout << res.seqBase_.name_ << "\t" << res.genomeCnt_ << std::endl;
	}

	return 0;
}

int readSimulatorRunner::simRandomlySampleGenomes(
		const njh::progutils::CmdArgs & inputCommands) {
	uint32_t sims = 10;
	std::vector<uint64_t> totalGenomesCounts{6,60,600,6000};
	std::vector<double> minorAlleleFreqsVec(49);
	njh::iota(minorAlleleFreqsVec, 1.0);
	std::set<double>minorAlleleFreqs(minorAlleleFreqsVec.begin(), minorAlleleFreqsVec.end());

	OutOptions outOpts(bfs::path(""), ".tab.txt");

	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(totalGenomesCounts, "--totalGenomes", "Total Starting Genomes", njh::progutils::ProgramSetUp::ConCheckCase::GREATERZERO);
	setUp.setOption(minorAlleleFreqs, "--minorAlleleFreqs", "Minor Allele Freqs", njh::progutils::ProgramSetUp::ConCheckCase::GREATERZERO);
	setUp.setOption(sims, "--numberOfSims", "Number Of Sims");
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	for(const auto & minorAF : minorAlleleFreqs){
		if(minorAF >=50){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "minor allele frequency can't be greater than 50"<< "\n";
			throw std::runtime_error{ss.str()};
		}
	}
	OutputStream out(outOpts);
	out << "TotalGenomes\tsim\texp\tobs" << std::endl;
	for(const auto sim : iter::range(sims)){
		for(const auto & gCount : totalGenomesCounts){
			for(const auto & maf : minorAlleleFreqs){
				seqInfo major{"major"};
				seqInfo minor{"minor"};
				major.cnt_ = 100 - maf;
				major.frac_ = major.cnt_/100;
				minor.cnt_ = maf;
				minor.frac_ = minor.cnt_/100;
				std::vector<seqInfo> seqs{major, minor};
				auto results = PCRSimulator::randomlySampleGenomes(seqs, gCount);
				for(const auto & res : results){
					out << gCount
							<< "\t" << sim
							<< "\t" << res.seqBase_.cnt_
							<< "\t" << (static_cast<double>(res.genomeCnt_)/gCount) * 100
							<< std::endl;
				}
			}
		}
	}

	return 0;
}




} // namespace njhseq
