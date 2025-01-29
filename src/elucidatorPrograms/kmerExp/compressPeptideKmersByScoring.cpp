//
// Created by Nicholas Hathaway on 1/28/25.
//
#include <njhseq/alignment/aligner/aligner.hpp>
#include <njhseq/concurrency/PairwisePairFactory.hpp>
#include <njhseq/concurrency/pools/AlignerPool.hpp>

#include "kmerExp.hpp"
#include "elucidator/objects/BioDataObject.h"
#include <njhseq/objects/dataContainers/graphs/UndirWeightedGraph.hpp>
#include <njhseq/objects/seqObjects/seqKmers/KmerVecUtils.hpp>
#include <njhseq/PopulationGenetics/PopGenCalcs.hpp>
#include <njhseq/objects/helperObjects/PeptideLibraryReducer.hpp>

namespace njhseq {



int kmerExpRunner::compressPeptidesByReducedLibrary(const njh::progutils::CmdArgs &inputCommands) {
	std::string reduction = "UNIPROT18";
	OutOptions keyOutOpts("", ".tsv");
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processDefaultReader(true);
	setUp.setOption(reduction, "--reduction", "reduction library to use, options are: "  + njh::conToStr(PeptideLibraryReducer::availableReductions_, ","));
	setUp.setOption(keyOutOpts.outFilename_, "--outKey", "output a file with the reduction key");
	keyOutOpts.transferOverwriteOpts(setUp.pars_.ioOptions_.out_);
	setUp.finishSetUp(std::cout);

	PeptideLibraryReducer reducer;
	reducer.setReduction(reduction);
	if (setUp.pars_.debug_) {
		std::cout << "residueToReduction: " << std::endl;
		for (const auto & red : reducer.residueToReduction) {
			std::cout << red.first << " " << red.second << std::endl;
		}
		std::cout << "reductionToResidues: " << std::endl;
		for (const auto & red : reducer.reductionToResidues) {
			std::cout << red.first << " " << njh::conToStr(red.second, ",") << std::endl;
		}
		std::cout << std::endl;
	}
	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	reader.openOut();

	std::unordered_map<std::string, VecStr> reducedPeptideKey;
	std::unique_ptr<OutputStream> out;
	if (!keyOutOpts.outFilename_.empty()) {
		out = std::make_unique<OutputStream>(keyOutOpts);
	}

	seqInfo seq;
	while (reader.readNextRead(seq)) {
		auto reduced = reducer.reverseSimpleFirst(reducer.reducePeptide(seq.seq_));
		if (!keyOutOpts.outFilename_.empty()) {
			reducedPeptideKey[reduced].emplace_back(seq.seq_);
		}
		seq.seq_ = reduced;
		reader.write(seq);
	}

	if (!keyOutOpts.outFilename_.empty()) {
		*out << "reduced_peptide\tinput_peptides" << std::endl;
		for (const auto & key : reducedPeptideKey) {
			*out << key.first << '\t' << njh::conToStr(key.second, ",") << std::endl;
		}
	}
	return 0;
}

/*
int kmerExpRunner::compressPeptideKmersByReducedLibrary(const njh::progutils::CmdArgs &inputCommands) {
  uint32_t kmerLength = 16;
	uint32_t numThreads = 1;
	uint32_t batchNumber = 100;
	njhUndirWeightedGraph<int32_t, std::string>::dbscanPars dbscan_pars;
	dbscan_pars.eps_ = 0;
	dbscan_pars.minEpNeighbors_ = 5;

	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(batchNumber, "--batchNumber", "batchNumber");

	setUp.setOption(kmerLength, "--kmerLength", "kmer Length");
	setUp.processReadInNames(true);
	setUp.setOption(numThreads, "--numThreads", "number of threads");
	setUp.setOption(dbscan_pars.eps_, "--epsilon", "epislon distance for making connections, will make connections if less than < -BLOSUM62 score since distance algorithm makes connections with a distance < a score");
	setUp.setOption(dbscan_pars.minEpNeighbors_, "--minNeighbors", "minimum number of of neighbors");


	setUp.finishSetUp(std::cout);


	auto inputSeqs = createKmerReadVec(SeqInput::getSeqVec<readObject>(setUp.pars_.ioOptions_), kmerLength, false);

	std::unordered_set<std::string> allKmers;
	for (const auto & seq : inputSeqs) {
		for (const auto & kmer : seq->kInfo_.kmers_) {
			allKmers.emplace(kmer.first);
		}
	}

	njhUndirWeightedGraph<int32_t, std::string> kmerGraph;

	for (const auto & kmer : allKmers) {
		kmerGraph.nodes_.emplace_back(std::make_shared<njhUndirWeightedGraph<int32_t, std::string>::node>(estd::to_string(kmerGraph.nodes_.size()), kmer));
	}

	PairwisePairFactory pairFactory(kmerGraph.nodes_.size());

	njhseq::concurrent::AlignerPool aligner_pool(
		kmerLength * 2, gapScoringParameters(10, 1), substituteMatrix::createBlosum62(), numThreads
	);
	aligner_pool.initAligners();
	std::mutex kmerGraph_mtx;
	njh::ProgressBar progress_bar(pairFactory.totalCompares_);

	std::function<void()> fillInDistMatrix = [&kmerGraph,&kmerGraph_mtx,&pairFactory,&aligner_pool,&batchNumber, &dbscan_pars,&progress_bar,&setUp]() {
		PairwisePairFactory::PairwisePairVec pairs;
		auto current_aligner = aligner_pool.popAligner();
		while (pairFactory.setNextPairs(pairs, batchNumber)) {
			PairwisePairFactory::PairwisePairVec passing_pairs;
			for (const auto & pair : pairs.pairs_) {
				if (setUp.pars_.verbose_) {
					progress_bar.outputProgAdd(std::cout, 1, true);
				}
				current_aligner->alignObjectA_.seqBase_.seq_ = kmerGraph.nodes_[pair.col_]->value_;
				current_aligner->alignObjectA_.seqBase_.seq_ = kmerGraph.nodes_[pair.row_]->value_;
				current_aligner->scoreAlignment(true);
				//negate the score so a better score will be smaller which is what is expected for the distance algo
				current_aligner->parts_.score_ = -current_aligner->parts_.score_;
				if (current_aligner->parts_.score_ <= dbscan_pars.eps_) {
					passing_pairs.pairs_.emplace_back(pair);
				}
			}
			if (!passing_pairs.pairs_.empty()) {
				std::lock_guard lock(kmerGraph_mtx);
				for (const auto & p : passing_pairs.pairs_) {
					kmerGraph.addEdge(kmerGraph.nodes_[p.col_]->name_, kmerGraph.nodes_[p.row_]->name_, current_aligner->parts_.score_);
				}
			}
		}
	};

	njh::concurrent::runVoidFunctionThreaded(fillInDistMatrix, numThreads);

	return 0;
}

int kmerExpRunner::compressPeptideKmersByScoring(const njh::progutils::CmdArgs &inputCommands) {
  uint32_t kmerLength = 16;
	uint32_t numThreads = 1;
	uint32_t batchNumber = 100;
	njhUndirWeightedGraph<int32_t, std::string>::dbscanPars dbscan_pars;
	dbscan_pars.eps_ = 0;
	dbscan_pars.minEpNeighbors_ = 5;

	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(batchNumber, "--batchNumber", "batchNumber");

	setUp.setOption(kmerLength, "--kmerLength", "kmer Length");
	setUp.processReadInNames(true);
	setUp.setOption(numThreads, "--numThreads", "number of threads");
	setUp.setOption(dbscan_pars.eps_, "--epsilon", "epislon distance for making connections, will make connections if less than < -BLOSUM62 score since distance algorithm makes connections with a distance < a score");
	setUp.setOption(dbscan_pars.minEpNeighbors_, "--minNeighbors", "minimum number of of neighbors");


	setUp.finishSetUp(std::cout);


	auto inputSeqs = createKmerReadVec(SeqInput::getSeqVec<readObject>(setUp.pars_.ioOptions_), kmerLength, false);

	std::unordered_set<std::string> allKmers;
	for (const auto & seq : inputSeqs) {
		for (const auto & kmer : seq->kInfo_.kmers_) {
			allKmers.emplace(kmer.first);
		}
	}

	njhUndirWeightedGraph<int32_t, std::string> kmerGraph;

	for (const auto & kmer : allKmers) {
		kmerGraph.nodes_.emplace_back(std::make_shared<njhUndirWeightedGraph<int32_t, std::string>::node>(estd::to_string(kmerGraph.nodes_.size()), kmer));
	}

	PairwisePairFactory pairFactory(kmerGraph.nodes_.size());

	njhseq::concurrent::AlignerPool aligner_pool(
		kmerLength * 2, gapScoringParameters(10, 1), substituteMatrix::createBlosum62(), numThreads
	);
	aligner_pool.initAligners();
	std::mutex kmerGraph_mtx;
	njh::ProgressBar progress_bar(pairFactory.totalCompares_);

	std::function<void()> fillInDistMatrix = [&kmerGraph,&kmerGraph_mtx,&pairFactory,&aligner_pool,&batchNumber, &dbscan_pars,&progress_bar,&setUp]() {
		PairwisePairFactory::PairwisePairVec pairs;
		auto current_aligner = aligner_pool.popAligner();
		while (pairFactory.setNextPairs(pairs, batchNumber)) {
			PairwisePairFactory::PairwisePairVec passing_pairs;
			for (const auto & pair : pairs.pairs_) {
				if (setUp.pars_.verbose_) {
					progress_bar.outputProgAdd(std::cout, 1, true);
				}
				current_aligner->alignObjectA_.seqBase_.seq_ = kmerGraph.nodes_[pair.col_]->value_;
				current_aligner->alignObjectA_.seqBase_.seq_ = kmerGraph.nodes_[pair.row_]->value_;
				current_aligner->scoreAlignment(true);
				//negate the score so a better score will be smaller which is what is expected for the distance algo
				current_aligner->parts_.score_ = -current_aligner->parts_.score_;
				if (current_aligner->parts_.score_ <= dbscan_pars.eps_) {
					passing_pairs.pairs_.emplace_back(pair);
				}
			}
			if (!passing_pairs.pairs_.empty()) {
				std::lock_guard lock(kmerGraph_mtx);
				for (const auto & p : passing_pairs.pairs_) {
					kmerGraph.addEdge(kmerGraph.nodes_[p.col_]->name_, kmerGraph.nodes_[p.row_]->name_, current_aligner->parts_.score_);
				}
			}
		}
	};

	njh::concurrent::runVoidFunctionThreaded(fillInDistMatrix, numThreads);

	return 0;
}
*/

int kmerExpRunner::generatePeptideKmerDistMatrixByScore(const njh::progutils::CmdArgs &inputCommands) {
	OutOptions outOpts(bfs::path(""), ".tsv");

	bool noHeader = false;
	uint32_t numThreads = 1;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processReadInNames(true);

	setUp.setOption(noHeader, "--noHeader", "no header");
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
	std::vector<std::vector<int64_t>> distMatrix(input.size(), std::vector<int64_t>(input.size(), 0));
	PairwisePairFactory pairFactory(input.size());

	njhseq::concurrent::AlignerPool aligner_pool(
		len(input.front()), gapScoringParameters(10, 1), substituteMatrix::createBlosum62(), numThreads
	);
	aligner_pool.initAligners();
	std::function<void()> fillInDistMatrix = [&distMatrix,&pairFactory,&input,&aligner_pool]() {
		PairwisePairFactory::PairwisePair pair;
		auto current_aligner = aligner_pool.popAligner();
		while (pairFactory.setNextPair(pair)) {
			current_aligner->noAlignSetAndScore(input[pair.col_], input[pair.row_]);
			distMatrix[pair.col_][pair.row_] = current_aligner->parts_.score_;
			distMatrix[pair.row_][pair.col_] = current_aligner->parts_.score_;
		}
	};

	njh::concurrent::runVoidFunctionThreaded(fillInDistMatrix, numThreads);

	{
		auto current_aligner = aligner_pool.popAligner();

		//fill in own score
		for (uint32_t i = 0; i < input.size(); ++i) {
			current_aligner->alignObjectA_ = input[i];
			current_aligner->alignObjectB_ = input[i];
			current_aligner->scoreAlignment(true);
			distMatrix[i][i] = current_aligner->parts_.score_;
			distMatrix[i][i] = current_aligner->parts_.score_;
		}
	}
	if (!noHeader) {
		out << "name" << '\t' << njh::conToStr(readVec::getNames(input), "\t") << std::endl;
	}
	for (const auto & rowEnum : iter::enumerate(distMatrix)) {
		if (!noHeader) {
			out << input[rowEnum.index].name_ << '\t';
		}
		out << njh::conToStr(rowEnum.element, "\t") << std::endl;
	}
	return 0;
}



}// namespace njhseq
