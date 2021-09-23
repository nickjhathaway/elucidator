/*
 * popGenExp_quickHaplotypeInformation.cpp
 *
 *  Created on: Jun 23, 2021
 *      Author: nick
 */
#include "popGenExp.hpp"
#include "elucidator/objects/BioDataObject.h"
#include "elucidator/objects/counters/DNABaseCounter.hpp"
#include "elucidator/PopulationGenetics.h"

#include "elucidator/objects/seqContainers/CollapsedHaps.hpp"
#include <njhseq/concurrency/PairwisePairFactory.hpp>

namespace njhseq {

int popGenExpRunner::quickHaplotypeInformation(const njh::progutils::CmdArgs & inputCommands) {
	std::string identifier = "region";
	uint32_t numThreads = 1;
	bool getPairwiseComps = false;
	seqSetUp setUp(inputCommands);
	setUp.pars_.ioOptions_.out_.outFilename_ = "";
	setUp.pars_.ioOptions_.out_.outFileFormat_ = "tab";
	setUp.pars_.ioOptions_.out_.outExtention_ = ".tab.txt";
	setUp.processReadInNames(true);
	setUp.setOption(getPairwiseComps, "--getPairwiseComps", "Get Pairwise Comps");
	setUp.setOption(identifier, "--identifier", "Give a identifier name for info");
	setUp.processWritingOptions(setUp.pars_.ioOptions_.out_);
	setUp.processVerbose();
	setUp.processAlignerDefualts();
	setUp.finishSetUp(std::cout);
	bool writeHeader = true;
	if(bfs::exists(setUp.pars_.ioOptions_.out_.outName()) && setUp.pars_.ioOptions_.out_.append_){
		writeHeader = false;
	}
	auto inputSeqs = CollapsedHaps::readInReads(setUp.pars_.ioOptions_);
	inputSeqs.verbose_ = setUp.pars_.verbose_;
	OutputStream out(setUp.pars_.ioOptions_.out_);
	CollapsedHaps::AvgPairwiseMeasures avgPMeasures;
	if(getPairwiseComps){
		uint64_t maxLen = readVec::getMaxLength(inputSeqs.seqs_);
		aligner alignerObj(maxLen, setUp.pars_.gapInfo_, setUp.pars_.scoring_);
		alignerObj.processAlnInfoInput(setUp.pars_.alnInfoDirName_, setUp.pars_.verbose_);
		auto allComps = inputSeqs.getPairwiseComps(alignerObj, numThreads);
		alignerObj.processAlnInfoOutput(setUp.pars_.outAlnInfoDirName_, setUp.pars_.verbose_);
		avgPMeasures = inputSeqs.getAvgPairwiseMeasures(allComps);
	}

	if(writeHeader){
		out << "id\ttotalHaplotypes\tuniqueHaplotypes\tsinglets\tdoublets\texpShannonEntropy\tShannonEntropyE\teffectiveNumOfAlleles\the\tlengthPolymorphism" ;
		if(getPairwiseComps){
			out << "\tavgPercentID\tavgNumOfDiffs";
		}
		out << std::endl;
	}
	std::unordered_map<uint32_t, uint32_t> readLens = inputSeqs.getReadLenMap();
	inputSeqs.setFrequencies();
	auto divMeasures = PopGenCalculator::getGeneralMeasuresOfDiversity(inputSeqs.seqs_);
	out << identifier
			<< "\t" << inputSeqs.getTotalHapCount()
			<< "\t" << inputSeqs.seqs_.size()
			<< "\t" << divMeasures.singlets_
			<< "\t" << divMeasures.doublets_
			<< "\t" << divMeasures.expShannonEntropy_
			<< "\t" << divMeasures.ShannonEntropyE_
			<< "\t" << divMeasures.effectiveNumOfAlleles_
			<< "\t" << divMeasures.heterozygostiy_
			<< "\t" << (readLens.size() > 1 ? "true" : "false");
	if(getPairwiseComps){
		out << "\t" << avgPMeasures.avgPercentId << "\t" << avgPMeasures.avgNumOfDiffs;
	}
	out << std::endl;

	return 0;
}


}  // namespace njhseq

