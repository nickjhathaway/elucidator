//
// Created by Nicholas Hathaway on 8/20/25.
//



#include "seqUtilsInfoRunner.hpp"
#include <njhseq/IO/SeqIO/SeqIO.hpp>
#include <njhseq/seqToolsUtils/tandemRepeatUtils.hpp>
#include <njhseq/concurrency/PairwisePairFactory.hpp>
#include <njhseq/concurrency/pools/AlignerPool.hpp>


namespace njhseq {

int seqUtilsInfoRunner::countEditDistancesSameLenSeqs(const njh::progutils::CmdArgs &inputCommands) {
	OutOptions outOpts(bfs::path(""), ".tsv");

	uint32_t numThreads = 1;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processReadInNames(true);
	setUp.processWritingOptions(outOpts);
	setUp.setOption(numThreads, "--numThreads", "number of threads");
	setUp.finishSetUp(std::cout);


	auto input = SeqInput::getSeqVec<seqInfo>(setUp.pars_.ioOptions_);
	for (const auto & seq : input) {
		if (len(seq) != len(input.front())) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << ", error "
			<< ", all input has to be the same length, seq: " << seq.name_
			<< " has len: " << len(seq) << " which is different from "
			<< input.front().name_ << " which has len " << len(input.front()) << "\n";
			throw std::runtime_error{ss.str()};
		}
	}
	OutputStream out(outOpts);
	PairwisePairFactory pairFactory(input.size());
	//set up progress bar
	njh::ProgressBar pBar(pairFactory.totalCompares_);

	njhseq::concurrent::AlignerPool aligner_pool(
		len(input.front()), gapScoringParameters(10, 1), substituteMatrix::createScoreMatrix(1, 0, false, false, true), numThreads
	);
	aligner_pool.initAligners();
	std::unordered_map<uint32_t, uint32_t> counts;
	std::mutex counts_mut;
	std::function<void()> getEditDistances = [&counts,&counts_mut,&pairFactory,&input,&aligner_pool,&setUp, &pBar]() {
		PairwisePairFactory::PairwisePair pair;
		auto current_aligner = aligner_pool.popAligner();
		std::unordered_map<uint32_t, uint32_t> current_counts;
		while (pairFactory.setNextPair(pair)) {
			if (setUp.pars_.verbose_) {
				pBar.outputProgAdd(std::cout, 1, true);
			}
			current_aligner->noAlignSetAndScore(input[pair.col_], input[pair.row_]);
			//with match score being 1 and mismatch being 0, the edit distance (number of snps) will be length minus score
			++current_counts[input[pair.col_].seq_.size() - current_aligner->parts_.score_];
			// current_aligner->profilePrimerAlignment(input[pair.col_], input[pair.row_]);
			// ++current_counts[current_aligner->comp_.hqMismatches_];
		}
		{
			std::lock_guard lock(counts_mut);
			for (const auto & count : current_counts) {
				counts[count.first] += count.second;
			}
		}
	};

	njh::concurrent::runVoidFunctionThreaded(getEditDistances, numThreads);
	out << "edit_distance\tcount" << std::endl;
	auto counts_key = njh::getSetOfMapKeys(counts);
	for (const auto & dist : counts_key) {
		out << dist << "\t" << counts[dist] << std::endl;
	}
	return 0;
}



} //namespace njhseq


