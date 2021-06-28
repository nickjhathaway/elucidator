/*
 * popGenExp_callVariantsAgainstRefSeqIndividual.cpp
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
#include <njhseq/GenomeUtils.h>


namespace njhseq {

int popGenExpRunner::callVariantsAgainstRefSeqIndividual(const njh::progutils::CmdArgs & inputCommands) {
	OutOptions outOpts(bfs::path(""), ".tab.txt");
	bfs::path bedFnp = "";
	uint32_t outwardsExpand = 5;
	bfs::path genomeFnp = "";
	uint32_t numThreads = 1;
	MultiGenomeMapper::inputParameters gPars;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processReadInNames(true);
	setUp.processWritingOptions(outOpts);
	setUp.processAlnInfoInput();
	setUp.setOption(outwardsExpand, "--outwardsExpand", "The amount to expand outwards from given region when determining variants positions with extracted ref seq");
	setUp.setOption(numThreads, "--numThreads", "number of Threads");

	setUp.setOption(bedFnp,    "--bed",    "A bed file of the location for the extraction", true);
	setUp.setOption(genomeFnp, "--genome", "A reference genome to compare against", true);
	if(!bfs::is_regular_file(genomeFnp)){
		setUp.failed_ = true;
		setUp.addWarning(njh::pasteAsStr(genomeFnp, " should be a file, not a directory"));
	}
	gPars.genomeDir_ = njh::files::normalize(genomeFnp.parent_path());
	gPars.primaryGenome_ = bfs::basename(genomeFnp);
	gPars.primaryGenome_ = gPars.primaryGenome_.substr(0, gPars.primaryGenome_.rfind("."));
	gPars.selectedGenomes_ = {gPars.primaryGenome_};

	setUp.finishSetUp(std::cout);


	njh::files::checkExistenceThrow(bedFnp,    __PRETTY_FUNCTION__);
	njh::files::checkExistenceThrow(genomeFnp, __PRETTY_FUNCTION__);
	OutputStream out(outOpts);
	auto beds = getBeds(bedFnp);
	if(beds.empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error no records found in " << bedFnp << "\n";
		throw std::runtime_error{ss.str()};
	}
	MultiGenomeMapper gMapper(gPars);
	gMapper.loadInGenomes();
	gMapper.setUpGenomes();
	auto gPos = beds.front();
	auto inputSeqs = CollapsedHaps::readInReads(setUp.pars_.ioOptions_);
	uint64_t maxLen = readVec::getMaxLength(inputSeqs.seqs_);
	inputSeqs.setFrequencies();
	if (gPos->reverseStrand()) {
		inputSeqs.revCompSeqs();
	}

	uint32_t oldLen = gPos->length();
	BedUtility::extendLeftRight(*gPos, outwardsExpand, outwardsExpand,
			gMapper.genomes_.at(gPars.primaryGenome_)->chromosomeLengths_.at(
					gPos->chrom_));
	if(oldLen == gPos->score_){
		gPos->score_ = gPos->length();
	}
	GenomicRegion refRegion(*gPos);
	//put genomic location into reference orientation
	refRegion.reverseSrand_ = false;
	TwoBit::TwoBitFile tReader(gMapper.genomes_.at(gPars.primaryGenome_)->fnpTwoBit_);
	auto refSeq = refRegion.extractSeq(tReader);
	readVec::getMaxLength(refSeq, maxLen);
	aligner alignerObj(maxLen, gapScoringParameters(5,1,0,0,0,0), substituteMatrix(2,-2), false);
	alignerObj.weighHomopolymers_ = false;
	alignerObj.processAlnInfoInput(setUp.pars_.alnInfoDirName_, setUp.pars_.verbose_);

	out << "#" << njh::conToStr(DistanceComp::BasicInfoHeader(), "\t") << "\n";
	auto cmps = inputSeqs.getCompsAgainstRef(refSeq, alignerObj, numThreads);
	for(const auto pos : iter::range(inputSeqs.seqs_.size())){
		for(const auto & name : inputSeqs.names_[pos]){
			cmps[pos].distances_.writeBasicInfo(out, *gPos, name);
		}
	}
	alignerObj.processAlnInfoOutput(setUp.pars_.outAlnInfoDirName_, setUp.pars_.verbose_);
	return 0;
}




} //namespace njhseq

