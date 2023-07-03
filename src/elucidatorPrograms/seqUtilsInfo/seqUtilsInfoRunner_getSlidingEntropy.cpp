/*
 * seqUtilsInfoRunner_getSlidingEntropy.cpp
 *
 *  Created on: Jun 27, 2020
 *      Author: nick
 */





#include "seqUtilsInfoRunner.hpp"
#include <njhseq/IO/SeqIO/SeqIO.hpp>
#include <njhseq/objects/seqObjects/seqKmers.h>
#include <njhseq/objects/BioDataObject/reading.hpp>


namespace njhseq {


int seqUtilsInfoRunner::getSlidingEntropy(const njh::progutils::CmdArgs & inputCommands) {

	OutOptions outOpts(bfs::path(""), ".tab.txt");
	uint32_t windowStep = 5;
	uint32_t windowSize = 50;
	uint32_t kLen = 1;
	seqSetUp setUp(inputCommands);
	setUp.description_ = "Trim front and back of sequences with a sliding window for low entropy, seqs will be removed if all is low entropy";

	setUp.processVerbose();
	setUp.processDebug();
	setUp.processWritingOptions(outOpts);
	setUp.setOption(windowStep, "--windowStep", "Window Step");
	setUp.setOption(windowSize, "--windowSize", "Window Size");
	setUp.setOption(kLen, "--kLen", "Kmer length for entropy calculation");

	setUp.processReadInNames(true);
	setUp.finishSetUp(std::cout);

	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	seqInfo seq;
	OutputStream out(outOpts);
	out << "name\tposition\tsubseq\tentropy" << std::endl;
	while(reader.readNextRead(seq)){
		if(len(seq) >= windowSize){
			for (auto pos : iter::range<uint32_t>(0, len(seq) - windowSize + 1, windowStep)) {
				kmerInfo kInfo(seq.seq_.substr(pos, windowSize), kLen, false);
				out << seq.name_
						<< "\t" << pos
						<< "\t" << seq.seq_.substr(pos, windowSize)
						<< "\t" << kInfo.computeKmerEntropy() << std::endl;
			}
		}
	}

	if(setUp.pars_.verbose_){
		setUp.logRunTime(std::cout);
	}
	return 0;
}


int seqUtilsInfoRunner::getSlidingEntropyGenomicRegion(const njh::progutils::CmdArgs & inputCommands) {

	OutOptions outOpts(bfs::path(""), ".tab.txt");
	uint32_t windowStep = 5;
	uint32_t windowSize = 50;
	uint32_t kLen = 1;
	bfs::path bedFnp;
	bfs::path twoBitFnp;
	seqSetUp setUp(inputCommands);
	setUp.description_ = "Trim front and back of sequences with a sliding window for low entropy, seqs will be removed if all is low entropy";

	setUp.processVerbose();
	setUp.processDebug();
	setUp.processWritingOptions(outOpts);
	setUp.setOption(windowStep, "--windowStep", "Window Step");
	setUp.setOption(windowSize, "--windowSize", "Window Size");
	setUp.setOption(kLen, "--kLen", "Kmer length for entropy calculation");

	setUp.setOption(bedFnp, "--bed", "bedFnp", true);
	setUp.setOption(twoBitFnp, "--twoBitFnp", "twoBitFnp", true);

	setUp.finishSetUp(std::cout);
	TwoBit::TwoBitFile treader(twoBitFnp);
	auto inputRegions = getBed3s(bedFnp);

	OutputStream out(outOpts);
	out << "#chrom\tstart\tend\tname\tscore\tstrand" << std::endl;
	for(const auto & region : inputRegions){
		auto gRegion = GenomicRegion(*region);
		gRegion.reverseSrand_ = false; //make plus strand since entropy won't change based on whether it's reverse comp or not
		auto seq = gRegion.extractSeq(treader);
		if(len(seq) >= windowSize){
			for (auto pos : iter::range<uint32_t>(0, len(seq) - windowSize + 1, windowStep)) {
				kmerInfo kInfo(seq.seq_.substr(pos, windowSize), kLen, false);
				out << region->chrom_
						<< "\t" << region->chromStart_ + pos
						<< "\t" << region->chromStart_ + pos + windowSize
						<< "\t" << seq.seq_.substr(pos, windowSize)
						<< "\t" << kInfo.computeKmerEntropy()
						<< "\t" << '+' << std::endl;
			}
		}
	}

	if(setUp.pars_.verbose_){
		setUp.logRunTime(std::cout);
	}
	return 0;
}





}  // namespace njhseq

