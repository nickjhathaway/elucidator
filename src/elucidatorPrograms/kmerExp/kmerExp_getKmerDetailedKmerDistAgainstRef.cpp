/*
 * kmerExp_getKmerDetailedKmerDistAgainstRef.cpp
 *
 *  Created on: Oct 30, 2021
 *      Author: nick
 */




#include "kmerExp.hpp"
#include <njhseq/objects/kmer/KmerUtils.hpp>
#include "elucidator/objects/dataContainers.h"
#include "elucidator/simulation.h"
#include "elucidator/objects/seqObjects/seqKmers.h"

#include "elucidator/objects/MiscUtility/GenomeSeqSearch.hpp"
#include "elucidator/objects/BioDataObject.h"

#include <njhseq/helpers.h>

namespace njhseq {


int kmerExpRunner::getBestKmerDetailedKmerDistAgainstRef(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path inputFnp = "";
	OutOptions outOpts(bfs::path(""), ".tab.txt");
	uint32_t numThreads = 1;
	uint32_t kmerLength = 7;
	bool doNotSkipSameName = false;
	bool getRevComp = false;
	VecStr measureOpts{"similarity", "similarityLenAdjust", "uniqSimilarity", "uniqSimilarityLenAdjust"};
	std::string measureOpt = measureOpts.front();
	bool alignBest = false;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processReadInNames(true);
	setUp.processRefFilename(true);
	setUp.processWritingOptions(outOpts);
	setUp.setOption(doNotSkipSameName, "--doNotSkipSameName", "do Not Skip Same Name");
	setUp.setOption(alignBest, "--alignBest", "align the best hits");
	if(alignBest) {
		setUp.processAlignerDefualts();
	}
	setUp.setOption(measureOpt, "--measureOpt", "Measure option to compare on, options are: " + njh::conToStrEndSpecial(measureOpts, ", ", " or "));
	if(njh::notIn(measureOpt, measureOpts)) {
		setUp.addWarning(njh::pasteAsStr("Measure opts must be one of the following ", njh::conToStrEndSpecial(measureOpts, ", ", " or "), ", not: ", measureOpt));
		setUp.failed_ = true;
	}
	setUp.setOption(numThreads, "--numThreads", "number of threads to use");
	setUp.setOption(kmerLength, "--kmerLength", "kmer Length");
	setUp.setOption(getRevComp, "--getRevComp", "Compare Reverse Complement as well");

	setUp.finishSetUp(std::cout);

	setUp.pars_.refIoOptions_.includeWhiteSpaceInName_ = setUp.pars_.ioOptions_.includeWhiteSpaceInName_;
	auto refSeqs = createKmerReadVec(setUp.pars_.refIoOptions_, kmerLength, getRevComp);
	uint64_t maxLen = readVec::getMaxLength(refSeqs	);
	{
		SeqInput reader(setUp.pars_.ioOptions_);
		reader.openIn();
		seqInfo seq;
		while(reader.readNextRead(seq)) {
			readVec::getMaxLength(seq, maxLen	);
		}
	}
	std::unique_ptr<aligner> alignerPtr;
	if(alignBest) {
		alignerPtr = std::make_unique<aligner>(maxLen,setUp.pars_.gapInfo_, setUp.pars_.scoring_);
	}
	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();

	std::unique_ptr<OutputStream> alignOut;
	if(alignBest) {
		auto bestAlignsFnp = njh::files::prependFileBasename(outOpts.outName(), "bestAligns_");
		alignOut = std::make_unique<OutputStream>(bfs::path(bestAlignsFnp	).replace_extension(".fastq.gz"));
	}

	OutputStream out(outOpts);
	out << "name\tref\trevComp\ttotalShared\tsimilarity\tsimilarityLenAdjust\ttotalKmersIn1\ttotalKmersIn2\ttotalUniqueShared\tuniqSimilarity\tuniqSimilarityLenAdjust\ttotalUniqKmersIn1\ttotalUniqKmersIn2\ttotalUniq";
	if(alignBest) {
		out << "\tpercentIdentity";
	}
	out << "\n";

	std::function<bool(kmerInfo::DetailedKmerDist &, const kmerInfo::DetailedKmerDist &)> compToBest;


	if("similarity" == measureOpt) {
		compToBest = [](kmerInfo::DetailedKmerDist & bestSimScore, const kmerInfo::DetailedKmerDist & comp) {
			if(comp.getDistTotalShared() > bestSimScore.getDistTotalShared()) {
				bestSimScore = comp;
				return true;
			}
			return false;
		};
	} else if ("similarityLenAdjust" == measureOpt) {
		compToBest = [](kmerInfo::DetailedKmerDist & bestSimScore, const kmerInfo::DetailedKmerDist & comp) {
			if(comp.getDistTotalSharedLenAdjusted() > bestSimScore.getDistTotalSharedLenAdjusted()) {
				bestSimScore = comp;
				return true;
			}
			return false;
		};
	} else if ("uniqSimilarity" == measureOpt) {
		compToBest = [](kmerInfo::DetailedKmerDist & bestSimScore, const kmerInfo::DetailedKmerDist & comp) {
			if(comp.getDistUniqueShared() > bestSimScore.getDistUniqueShared()) {
				bestSimScore = comp;
				return true;
			}
			return false;
		};
	} else if ("uniqSimilarityLenAdjust" == measureOpt) {
		compToBest = [](kmerInfo::DetailedKmerDist & bestSimScore, const kmerInfo::DetailedKmerDist & comp) {
			if(comp.getDistUniqueSharedLenAdjusted() > bestSimScore.getDistUniqueSharedLenAdjusted()) {
				bestSimScore = comp;
				return true;
			}
			return false;
		};
	}


	std::mutex outMut;
	std::function<void()> getBestSimilarity = [&reader,&out,&outMut,
																			 &getRevComp,&kmerLength,
																			 &refSeqs, &compToBest,
																			 &doNotSkipSameName,
																			 &alignBest,
																			 &alignerPtr,
																			 &alignOut](){
		seqInfo seq;

		while(reader.readNextReadLock(seq)){
			kmerInfo::DetailedKmerDist bestSimilarityScore;
			bool bestSimSocreRevComp = false;
			uint32_t bestRefPos = std::numeric_limits<uint32_t>::max();
			kmerInfo kInfo(seq.seq_, kmerLength, false);
			for(const auto pos : iter::range(refSeqs.size())){
				const auto & refSeq = refSeqs[pos];
				if(!doNotSkipSameName && refSeq->seqBase_.name_ == seq.name_) {
					continue;
				}
				{
					auto similarity = kInfo.compareKmersDetailed(refSeq->kInfo_);
					if(compToBest(bestSimilarityScore, similarity)) {
						bestSimSocreRevComp = false;
						bestRefPos = pos;
					}
				}
				if(getRevComp){
					auto revDist = kInfo.compareKmersRevCompDetailed(refSeq->kInfo_);
					if(compToBest(bestSimilarityScore, revDist)) {
						bestSimSocreRevComp = true;
						bestRefPos = pos;
					}
				}
			}
			{
				std::lock_guard<std::mutex> lock(outMut);
				if(alignBest) {
					alignerPtr->alignCacheGlobal(refSeqs[bestRefPos]->seqBase_, seq);
					alignerPtr->profileAlignment(refSeqs[bestRefPos]->seqBase_, seq,false, false,false);
					alignerPtr->alignObjectA_.seqBase_.outPutFastq(*alignOut);
					alignerPtr->alignObjectB_.seqBase_.outPutFastq(*alignOut);
				}

				out << seq.name_
						<< "\t" << (bestRefPos != std::numeric_limits<uint32_t>::max() ? refSeqs[bestRefPos]->seqBase_.name_ : std::string("NA") )
						<< "\t" << njh::boolToStr(bestSimSocreRevComp)
						<< "\t" << bestSimilarityScore.totalShared_
						<< "\t" << bestSimilarityScore.getDistTotalShared()
						<< "\t" << bestSimilarityScore.getDistTotalSharedLenAdjusted()
						<< "\t" << bestSimilarityScore.totalKmersIn1_
						<< "\t" << bestSimilarityScore.totalKmersIn2_
						<< "\t" << bestSimilarityScore.totalUniqShared_
						<< "\t" << bestSimilarityScore.getDistUniqueShared()
						<< "\t" << bestSimilarityScore.getDistUniqueSharedLenAdjusted()
						<< "\t" << bestSimilarityScore.totalUniqKmersIn1_
						<< "\t" << bestSimilarityScore.totalUniqKmersIn2_
						<< "\t" << bestSimilarityScore.totalUniqBetween_;
				if(alignBest) {
					out << "\t" << alignerPtr->comp_.distances_.eventBasedIdentityHq_;
				}
				out
						<< "\n";
			}
		}
	};

	njh::concurrent::runVoidFunctionThreaded(getBestSimilarity, numThreads);



	return 0;
}

int kmerExpRunner::getKmerDetailedKmerDistAgainstRef(const njh::progutils::CmdArgs & inputCommands){
	bfs::path inputFnp = "";
	OutOptions outOpts(bfs::path(""), ".tab.txt");
	uint32_t numThreads = 1;
	uint32_t kmerLength = 7;

	bool getRevComp = false;


	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processReadInNames(true);
	setUp.processRefFilename(true);
	setUp.processWritingOptions(outOpts);

	setUp.setOption(numThreads, "--numThreads", "number of threads to use");
	setUp.setOption(kmerLength, "--kmerLength", "kmer Length");
	setUp.setOption(getRevComp, "--getRevComp", "Compare Reverse Complement as well");

	setUp.finishSetUp(std::cout);

	setUp.pars_.refIoOptions_.includeWhiteSpaceInName_ = setUp.pars_.ioOptions_.includeWhiteSpaceInName_;
	auto refSeqs = createKmerReadVec(setUp.pars_.refIoOptions_, kmerLength, getRevComp);

	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();

	OutputStream out(outOpts);
	out << "name\tref\trevComp\ttotalShared\tdist\tdistLenAdjust\ttotalKmersIn1\ttotalKmersIn2\ttotalUniqueShared\tuniqDist\tuniqDistLenAdjust\ttotalUniqKmersIn1\ttotalUniqKmersIn2\ttotalUniq" << "\n";

	std::mutex outMut;
	std::function<void()> getBestDist = [&reader,&out,&outMut,
																			 &getRevComp,&kmerLength,
																			 &refSeqs](){
		seqInfo seq;
		while(reader.readNextReadLock(seq)){

			kmerInfo kInfo(seq.seq_, kmerLength, false);
			for(const auto pos : iter::range(refSeqs.size())){
				const auto & refSeq = refSeqs[pos];
				auto dist = kInfo.compareKmersDetailed(refSeq->kInfo_);
				{
					std::lock_guard<std::mutex> lock(outMut);
					out << seq.name_
							<< "\t" << refSeq->seqBase_.name_
							<< "\t" << njh::boolToStr(false)
							<< "\t" << dist.totalShared_
							<< "\t" << dist.getDistTotalShared()
							<< "\t" << dist.getDistTotalSharedLenAdjusted()
							<< "\t" << dist.totalKmersIn1_
							<< "\t" << dist.totalKmersIn2_
							<< "\t" << dist.totalUniqShared_
							<< "\t" << dist.getDistUniqueShared()
							<< "\t" << dist.getDistUniqueSharedLenAdjusted()
							<< "\t" << dist.totalUniqKmersIn1_
							<< "\t" << dist.totalUniqKmersIn2_
							<< "\t" << dist.totalUniqBetween_
							<< std::endl;;
				}
				if(getRevComp){
					auto revDist = kInfo.compareKmersRevCompDetailed(refSeq->kInfo_);
					{
						std::lock_guard<std::mutex> lock(outMut);
						out << seq.name_
								<< "\t" << refSeq->seqBase_.name_
								<< "\t" << njh::boolToStr(true)
								<< "\t" << revDist.totalShared_
								<< "\t" << revDist.getDistTotalShared()
								<< "\t" << revDist.getDistTotalSharedLenAdjusted()
								<< "\t" << revDist.totalKmersIn1_
								<< "\t" << revDist.totalKmersIn2_
								<< "\t" << revDist.totalUniqShared_
								<< "\t" << revDist.getDistUniqueShared()
								<< "\t" << revDist.getDistUniqueSharedLenAdjusted()
								<< "\t" << revDist.totalUniqKmersIn1_
								<< "\t" << revDist.totalUniqKmersIn2_
								<< "\t" << revDist.totalUniqBetween_
								<< "\n";
					}
				}
			}
		}
	};

	njh::concurrent::runVoidFunctionThreaded(getBestDist, numThreads);



	return 0;
}


} // namespace njhseq



