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

int kmerExpRunner::comparePeptidesDistMatrixByBlosumScore(const njh::progutils::CmdArgs &inputCommands) {
	OutOptions outOpts(bfs::path(""), ".tsv");

	bool noHeader = false;
	uint32_t numThreads = 1;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processReadInNames(true);
	setUp.processWritingOptions(outOpts);
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
