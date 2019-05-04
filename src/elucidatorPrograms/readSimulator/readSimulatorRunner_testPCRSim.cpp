/*
 * readSimulatorRunner_simMixture.cpp
 *
 *  Created on: Sep 16, 2018
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

#include "readSimulatorRunner.hpp"
#include "elucidator/simulation.h"
#include <SeekDeep/objects/PrimersAndMids.hpp>

namespace njhseq {




int readSimulatorRunner::testPCRSim(
		const njh::progutils::CmdArgs & inputCommands) {
	uint32_t pcrRounds = 30;
	uint32_t initialPcrRounds = 10;
	long double errorRate = 3.5e-06;
	uint32_t numThreads = 1;
	uint32_t startingTemplate = 1000;
	uint32_t finalReadAmount = 5000;
	double pcrEfficiency = 0.95;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processSeq(true);
	setUp.processDirectoryOutputName(setUp.pars_.seqObj_.seqBase_.name_ + "_testPCRSim_TODAY", true );
	setUp.setOption(errorRate, "--errorRate", "Polymerase Error Rate");
	setUp.setOption(numThreads, "--numThreads", "Number Threads");
	setUp.setOption(pcrRounds, "--pcrRounds", "Total Number of PCR Rounds");
	setUp.setOption(initialPcrRounds, "--initialPcrRounds", "Number of Initial PCR Rounds");

	setUp.setOption(startingTemplate, "--startingTemplate", "Starting Template");
	setUp.setOption(finalReadAmount, "--finalReadAmount", "Final Read Amount");
	setUp.setOption(pcrEfficiency, "--pcrEfficiency", "PCR Efficiency, between 0-1, chance a product gets amplified");

	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);

	uint64_t intErrorRate = errorRate * std::numeric_limits<uint64_t>::max();
	PCRSimulator pcrSim(intErrorRate);
	pcrSim.verbose_ = setUp.pars_.verbose_;
	pcrSim.pcrEfficiency_ = pcrEfficiency;

	OutOptions pcrSeqFileOpts(
			njh::files::make_path(setUp.pars_.directoryName_, "pcrSeqs.fasta"));
	setUp.pars_.seqObj_.seqBase_.cnt_ = 1;
	setUp.pars_.seqObj_.seqBase_.frac_ = 1;

	pcrSim.simLibFast(std::vector<seqInfo> { setUp.pars_.seqObj_.seqBase_ },
			pcrSeqFileOpts, startingTemplate, finalReadAmount, pcrRounds,
			initialPcrRounds, numThreads);


	return 0;
}




} // namespace njhseq
