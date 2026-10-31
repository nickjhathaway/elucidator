/*
 * kmerSetExp_findUniqKmersBetweenSeqSets.cpp
 *
 *  Created on: Jul 10, 2021
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

#include "kmerSetExp.hpp"
#include <njhseq/objects/kmer/KmerUtils.hpp>
#include <njhseq/objects/kmer.h>

#include "elucidator/simulation.h"

#include "elucidator/objects/MiscUtility/GenomeSeqSearch.hpp"

#include <njhseq/objects/dataContainers/tables/TableReader.hpp>
#include <njhseq/IO/SeqIO/MultiSeqOutCache.hpp>
#include <njhseq/objects/BioDataObject/reading.hpp>
#include <njhseq/BamToolsUtils/LockedBamReader.hpp>
#include "elucidator/helpers/UniqueKmerSetHelper.hpp"

namespace njhseq {




int kmerSetExpRunner::filterUniqueKmerSetForEntropy(const njh::progutils::CmdArgs & inputCommands){
	bfs::path countTable = "";

	KmerGatherer::KmerGathererPars countPars;
	OutOptions outOpts(bfs::path(""), ".tab.txt.gz");
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(countTable, "--countTable", "countTable, 1)set,2)kmer", true);
	setUp.setOption(countPars.entropyFilter_, "--entropyCutOff", "entropy Cut Off to keep kmers (exclusive)");
	setUp.setOption(countPars.kmerLengthForEntropyCalc_, "--kmerLengthForEntropyCalc", "kmer Length For Entropy Calc");
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);
	njh::stopWatch watch;
	watch.setLapName("initial");
	watch.startNewLap("reading in unique kmer table");
	OutputStream out(outOpts);

	{
		TableReader uniqKmers(TableIOOpts::genTabFileIn(countTable, false));
		if (uniqKmers.header_.nCol() != 2) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "need to have 2 columns" << "\n";
			throw std::runtime_error{ss.str()};
		}
		VecStr row;
		while (uniqKmers.getNextRow(row)) {
			kmerInfo kinfo(row[1], countPars.kmerLengthForEntropyCalc_, false);
			if (kinfo.computeKmerEntropy() > countPars.entropyFilter_) {
				out << row[0]
						<< "\t" << row[1];
				if (setUp.pars_.debug_) {
					out << "\t" << kinfo.computeKmerEntropy();
				}

				out << std::endl;
			}
		}
	}
	return 0;
}

int kmerSetExpRunner::countingUniqKmersFromSetsInUnmappedAlns(const njh::progutils::CmdArgs & inputCommands){
	uint32_t numThreads = 1;
	bfs::path countTable = "";
	std::string sampleName;
	bool includeRevComp = false;
	OutOptions outOpts(bfs::path(""), ".tab.txt.gz");
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.pars_.ioOptions_.revComplMate_ = true;
	setUp.processReadInNames(VecStr{"--bam"}, true);
	setUp.setOption(countTable, "--countTable", "countTable, 1)set,2)kmer", true);
	setUp.setOption(sampleName, "--sampleName", "Name to add to output file", true);

	setUp.setOption(numThreads, "--numThreads", "numThreads");
	bool noInvertedKmerIndex = false;
	setUp.setOption(noInvertedKmerIndex, "--noInvertedKmerIndex",
			"Probe every k-mer set for every read k-mer instead of using one k-mer -> set(s) index. Identical results either way; this exists to fall back if that is ever in doubt");
	const bool useInvertedKmerIndex = !noInvertedKmerIndex;
	setUp.setOption(includeRevComp, "--includeRevComp", "include Rev Comp of the input seqs");




	setUp.processWritingOptions(outOpts);
	//setUp.processDirectoryOutputName("true");
	setUp.finishSetUp(std::cout);
	//setUp.startARunLog(setUp.pars_.directoryName_);
	njh::stopWatch watch;
	watch.setLapName("initial");

	auto klen =	KmerCollectionIO::getKmerLength(countTable);

	std::mutex mut;

	watch.startNewLap("reading in unique kmer table");
	std::unordered_map<std::string, std::unordered_set<uint64_t>> uniqueKmersPerSet = KmerCollectionIO::readHashedCollections(countTable);
	if(setUp.pars_.verbose_){
		std::cout << watch.getLapName() << "\t" << watch.timeLapFormatted() <<std::endl;
	}
	UniqueKmerSetHelper::KmerSetTally tally;
	OutputStream out(outOpts);
	VecStr names = getVectorOfMapKeys(uniqueKmersPerSet);

	LockedBamReader bamReader(setUp.pars_.ioOptions_.firstName_);

	const auto invIndex = UniqueKmerSetHelper::getUsableInvertedIndex(uniqueKmersPerSet, useInvertedKmerIndex);
	std::function<void()> readInComp = [&bamReader, &uniqueKmersPerSet, &tally, &mut, &klen, &includeRevComp, &invIndex]() {
		SimpleKmerHash hasher;
		UniqueKmerSetHelper::KmerSetTally currentTally;
		std::vector<BamTools::BamAlignment> alns;
		while(bamReader.readNextAlnBatch(alns, 10000)){
			for(const auto & bAln : alns){
				if(bAln.IsPrimaryAlignment() && !bAln.IsMapped()){
					std::unordered_map<uint64_t, uint64_t> hashedInputKmers;
					//unmapped reads are stored as sequenced, so the query bases are used as is
					UniqueKmerSetHelper::countReadKmers(bAln.QueryBases, klen, includeRevComp, hasher, hashedInputKmers);
					currentTally.addRead(hashedInputKmers, uniqueKmersPerSet, invIndex);
				}
			}
		}
		std::lock_guard<std::mutex> lock(mut);
		tally.add(currentTally);
	};
	njh::concurrent::runVoidFunctionThreaded(readInComp, numThreads);
	UniqueKmerSetHelper::KmerSetTally::writeHeader(out);
	njh::sort(names);
	for(const auto & name : names){
		tally.writeSetRow(out, sampleName, name, uniqueKmersPerSet.at(name));
	}

	return 0;
}

int kmerSetExpRunner::countingUniqKmersFromSets(const njh::progutils::CmdArgs & inputCommands){
	uint32_t numThreads = 1;
	bfs::path countTable = "";
	std::string sampleName;
	bool includeRevComp = false;
	OutOptions outOpts(bfs::path(""), ".tab.txt.gz");
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.pars_.ioOptions_.revComplMate_ = true;
	setUp.processReadInNames(true);
	setUp.setOption(countTable, "--countTable", "countTable, 1)set,2)kmer", true);
	setUp.setOption(sampleName, "--sampleName", "Name to add to output file", true);

	setUp.setOption(numThreads, "--numThreads", "numThreads");
	bool noInvertedKmerIndex = false;
	setUp.setOption(noInvertedKmerIndex, "--noInvertedKmerIndex",
			"Probe every k-mer set for every read k-mer instead of using one k-mer -> set(s) index. Identical results either way; this exists to fall back if that is ever in doubt");
	const bool useInvertedKmerIndex = !noInvertedKmerIndex;
	setUp.setOption(includeRevComp, "--includeRevComp", "include Rev Comp of the input seqs");




	setUp.processWritingOptions(outOpts);
	//setUp.processDirectoryOutputName("true");
	setUp.finishSetUp(std::cout);
	//setUp.startARunLog(setUp.pars_.directoryName_);
	njh::stopWatch watch;
	watch.setLapName("initial");
	auto klen =	KmerCollectionIO::getKmerLength(countTable);

	std::mutex mut;

	watch.startNewLap("reading in unique kmer table");
	std::unordered_map<std::string, std::unordered_set<uint64_t>> uniqueKmersPerSet = KmerCollectionIO::readHashedCollections(countTable);
	if(setUp.pars_.verbose_){
		std::cout << watch.getLapName() << "\t" << watch.timeLapFormatted() <<std::endl;
	}
	UniqueKmerSetHelper::KmerSetTally tally;
	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();

	OutputStream out(outOpts);
	const auto invIndex = UniqueKmerSetHelper::getUsableInvertedIndex(uniqueKmersPerSet, useInvertedKmerIndex);
	VecStr names = getVectorOfMapKeys(uniqueKmersPerSet);
	const bool pairedIn = setUp.pars_.ioOptions_.isPairedIn();

	std::function<void()> readInComp = [&reader, &uniqueKmersPerSet, &tally, &mut, &klen, &includeRevComp, &invIndex, &pairedIn]() {
		SimpleKmerHash hasher;
		UniqueKmerSetHelper::KmerSetTally currentTally;
		std::unordered_map<uint64_t, uint64_t> hashedInputKmers;
		if (pairedIn) {
			PairedRead pseq;
			while(reader.readNextReadLock(pseq)){
				hashedInputKmers.clear();
				UniqueKmerSetHelper::countPairedReadKmers(pseq, klen, includeRevComp, hasher, hashedInputKmers);
				currentTally.addRead(hashedInputKmers, uniqueKmersPerSet, invIndex);
			}
		} else {
			seqInfo seq;
			while(reader.readNextReadLock(seq)){
				hashedInputKmers.clear();
				UniqueKmerSetHelper::countReadKmers(seq.seq_, klen, includeRevComp, hasher, hashedInputKmers);
				currentTally.addRead(hashedInputKmers, uniqueKmersPerSet, invIndex);
			}
		}
		std::lock_guard<std::mutex> lock(mut);
		tally.add(currentTally);
	};
	njh::concurrent::runVoidFunctionThreaded(readInComp, numThreads);
	UniqueKmerSetHelper::KmerSetTally::writeHeader(out);
	njh::sort(names);
	for(const auto & name : names){
		tally.writeSetRow(out, sampleName, name, uniqueKmersPerSet.at(name));
	}
	return 0;
}



int kmerSetExpRunner::countingUniqKmersFromSetsInRegionsAlnsBestSet(const njh::progutils::CmdArgs & inputCommands){
	uint32_t numThreads = 1;
//	uint32_t numThreadsPerBam = 1;
	bfs::path countTable = "";

	UniqueKmerSetHelper::ProcessReadForExtractingPars extractingPars;
	extractingPars.compPars.hardCountOff = 10;
	extractingPars.compPars.fracCutOff = 0;
	extractingPars.compPars.kmerLengthForEntropyCalc_ = 2;
	extractingPars.compPars.entropyFilter_ = 1.20;
	extractingPars.compPars.initialExcludeHardCountOff = 60;
	extractingPars.compPars.initialExcludeFracCutOff = 0.25;
	OutOptions outOpts(bfs::path(""), ".tab.txt.gz");

	std::set <bfs::path> bams;
	bfs::path bedFnp;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.pars_.ioOptions_.revComplMate_ = true;
	setUp.setOption(bedFnp, "--bedFnp", "regions to count from", true);
	setUp.setOption(bams, "--bams", "input bams to work on", true);
	setUp.setOption(countTable, "--countTable", "countTable, 1)set,2)kmer", true);
	//setUp.setOption(extractingPars.compPars.sampleName, "--sampleName", "Name to add to output file", true);

	setUp.setOption(extractingPars.compPars.includeRevComp, "--includeRevComp", "include Rev Comp of the input seqs");
	setUp.setOption(extractingPars.compPars.entropyFilter_, "--entropyFilter", "entropy Filter_");
	setUp.setOption(extractingPars.compPars.kmerLengthForEntropyCalc_, "--kmerLengthForEntropyCalc", "kmerLength For Entropy Calculation");
	setUp.setOption(extractingPars.compPars.fracCutOff, "--readKmerFracCutOff", "per read Kmer Frac Cut Off");
	setUp.setOption(extractingPars.compPars.hardCountOff, "--hardCountOff", "hard Count Off, do not count sets unless greater thant his number");
	setUp.setOption(extractingPars.compPars.excludeSetNames, "--excludeSetNames", "names of sets of unique kmers to ignore from the --kmerTable");

	setUp.setOption(numThreads, "--numThreads", "numThreads");
//	setUp.setOption(numThreadsPerBam, "--numThreadsPerBam", "numThreadsPerBam");

	setUp.processWritingOptions(outOpts);
	//setUp.processDirectoryOutputName("true");
	setUp.finishSetUp(std::cout);
	//setUp.startARunLog(setUp.pars_.directoryName_);
	njh::stopWatch watch;
	watch.setLapName("initial");

	auto regions = getBeds(bedFnp);
	checkBamFilesForIndexesAndAbilityToOpen(std::vector<bfs::path> {bams.begin(), bams.end()}, numThreads  );
	OutputStream out(outOpts);

	extractingPars.compPars.klen =	KmerCollectionIO::getKmerLength(countTable);

	watch.startNewLap("reading in unique kmer table");
	std::unordered_map<std::string, std::unordered_set<uint64_t>> uniqueKmersPerSet = KmerCollectionIO::readHashedCollections(countTable);
	if(setUp.pars_.verbose_){
		std::cout << watch.getLapName() << "\t" << watch.timeLapFormatted() <<std::endl;
	}

	VecStr names = getVectorOfMapKeys(uniqueKmersPerSet);
	njh::naturalSortNameSet(names);



	out << "sample\ttotal\tset\treads\treadsFracTotal\treadsInForward\tforwardFrac" << std::endl;
	std::mutex outMut;
	njh::concurrent::LockableQueue<bfs::path> bamsQueue(bams);

	std::function<void()> countBam = [&bamsQueue, &out, &outMut,
					&uniqueKmersPerSet,
					&regions,
					&extractingPars, &names, &setUp]() {
		bfs::path bam;
		while(bamsQueue.getVal(bam)){
			MultiSeqIO seqOut;
			if(setUp.pars_.debug_){
				for (const auto &name : names) {
//					auto paired_seqOutOpts = SeqIOOptions::genPairedOutGz(njh::pasteAsStr(basename(bam), "-", name, "-paired"));
//					paired_seqOutOpts.out_.overWriteFile_ = true;
//					seqOut.addReader(name + "-paired", paired_seqOutOpts);

					auto single_seqOutOpts = SeqIOOptions::genFastqOutGz(njh::pasteAsStr(basename(bam), "-", name));
					single_seqOutOpts.out_.overWriteFile_ = true;
					seqOut.addReader(name, single_seqOutOpts);
				}
			}
			BamTools::BamReader bamReader;
			bamReader.Open(bam.string());
			bamReader.LocateIndex();
			checkBamOpenThrow(bamReader, bam);
			std::unordered_map<bool, std::unordered_map<std::string, uint64_t>> matchingCounts;
			uint64_t totalInput = 0;
			SimpleKmerHash hasher;
			BamTools::BamAlignment bAln;
			for(const auto & reg : regions){
				setBamFileRegionThrow(bamReader, GenomicRegion(*reg));
				while(bamReader.GetNextAlignment(bAln)){
					if (bAln.IsPrimaryAlignment() && bAln.IsMapped()) {
						seqInfo seq = bamAlnToSeqInfo(bAln);
						auto compRes = UniqueKmerSetHelper::compareReadToSets(seq, uniqueKmersPerSet, extractingPars.compPars, hasher);
						compRes.clearWinnerIfOtherSetsHit();
						if(setUp.pars_.debug_ && "undetermined" != compRes.winnerSet){
							seqOut.openWrite(compRes.winnerSet, seq);
						}
						++matchingCounts[compRes.winnerRevComp][compRes.winnerSet];
						++totalInput;
					}
				}
			}
			{
				std::lock_guard<std::mutex> lock(outMut);
				for (const auto &name: names) {
					uint64_t total = matchingCounts[true][name] + matchingCounts[false][name];
					uint64_t totalForward = matchingCounts[false][name];
					out << bfs::basename(bam)
							<< "\t" << totalInput
							<< "\t" << name
							<< "\t" << total
							<< "\t" << (totalInput > 0 ? total / static_cast<long double>(totalInput) : 0)
							<< "\t" << totalForward
							<< "\t" << (total > 0 ? totalForward / static_cast<long double>(total) : 0)
							<< std::endl;
				}
			}
		}
	};
	njh::concurrent::runVoidFunctionThreaded(countBam, numThreads);

	return 0;
}

int kmerSetExpRunner::countingUniqKmersFromSetsInUnmappedAlnsBestSet(const njh::progutils::CmdArgs & inputCommands){
	uint32_t numThreads = 1;
	uint32_t numThreadsPerBam = 1;
	bfs::path countTable = "";

	UniqueKmerSetHelper::ProcessReadForExtractingPars extractingPars;
	extractingPars.compPars.hardCountOff = 10;
	extractingPars.compPars.fracCutOff = 0;
	extractingPars.compPars.kmerLengthForEntropyCalc_ = 2;
	extractingPars.compPars.entropyFilter_ = 1.20;
	extractingPars.compPars.initialExcludeHardCountOff = 60;
	extractingPars.compPars.initialExcludeFracCutOff = 0.25;
	OutOptions outOpts(bfs::path(""), ".tab.txt.gz");

	std::set <bfs::path> bams;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.pars_.ioOptions_.revComplMate_ = true;
	setUp.setOption(bams, "--bams", "input bams to work on", true);
	setUp.setOption(countTable, "--countTable", "countTable, 1)set,2)kmer", true);
	//setUp.setOption(extractingPars.compPars.sampleName, "--sampleName", "Name to add to output file", true);

	setUp.setOption(extractingPars.compPars.includeRevComp, "--includeRevComp", "include Rev Comp of the input seqs");
	setUp.setOption(extractingPars.compPars.entropyFilter_, "--entropyFilter", "entropy Filter_");
	setUp.setOption(extractingPars.compPars.kmerLengthForEntropyCalc_, "--kmerLengthForEntropyCalc", "kmerLength For Entropy Calculation");
	setUp.setOption(extractingPars.compPars.fracCutOff, "--readKmerFracCutOff", "per read Kmer Frac Cut Off");
	setUp.setOption(extractingPars.compPars.hardCountOff, "--hardCountOff", "hard Count Off, do not count sets unless greater thant his number");
	setUp.setOption(extractingPars.compPars.excludeSetNames, "--excludeSetNames", "names of sets of unique kmers to ignore from the --kmerTable");

	setUp.setOption(numThreads, "--numThreads", "numThreads");
	setUp.setOption(numThreadsPerBam, "--numThreadsPerBam", "numThreadsPerBam");

	setUp.processWritingOptions(outOpts);
	//setUp.processDirectoryOutputName("true");
	setUp.finishSetUp(std::cout);
	//setUp.startARunLog(setUp.pars_.directoryName_);
	njh::stopWatch watch;
	watch.setLapName("initial");

	checkBamFilesForIndexesAndAbilityToOpen(std::vector<bfs::path> {bams.begin(), bams.end()}, numThreads * numThreadsPerBam );
	OutputStream out(outOpts);

	extractingPars.compPars.klen =	KmerCollectionIO::getKmerLength(countTable);

	watch.startNewLap("reading in unique kmer table");
	std::unordered_map<std::string, std::unordered_set<uint64_t>> uniqueKmersPerSet = KmerCollectionIO::readHashedCollections(countTable);
	if(setUp.pars_.verbose_){
		std::cout << watch.getLapName() << "\t" << watch.timeLapFormatted() <<std::endl;
	}

	VecStr names = getVectorOfMapKeys(uniqueKmersPerSet);
	njh::naturalSortNameSet(names);

	out << "sample\ttotal\ttotalUnmapped\tset\treads\treadsFracTotal\treadsFracUnmapped\treadsInForward\tforwardFrac" << std::endl;
	std::mutex outMut;
	njh::concurrent::LockableQueue<bfs::path> bamsQueue(bams);

	std::function<void()> countBam = [&bamsQueue, &out, &outMut,
																		&numThreadsPerBam, &uniqueKmersPerSet,
																		&extractingPars, &names]() {
		bfs::path bam;
		while(bamsQueue.getVal(bam)){
			LockedBamReader bamReader(bam);
			std::unordered_map<bool, std::unordered_map<std::string, uint64_t>> matchingCounts;
			uint64_t totalInput = 0;
			uint64_t totalUnmapped = 0;
			std::mutex matchingCountsMut;
			std::function<void()> readInComp = [&bamReader, &uniqueKmersPerSet, &matchingCounts,&matchingCountsMut,&extractingPars,
							&totalInput, &totalUnmapped]() {
				SimpleKmerHash hasher;
				std::unordered_map<bool, std::unordered_map<std::string, uint64_t>> currentMatchingCounts;
				uint64_t currentTotalInput = 0;
				uint64_t currentTotalUnmapped = 0;
				std::vector<BamTools::BamAlignment> alns;
				while (bamReader.readNextAlnBatch(alns, 10000)) {
					for (const auto &bAln: alns) {
						if (bAln.IsPrimaryAlignment() && !bAln.IsMapped()) {
							seqInfo seq = bamAlnToSeqInfo(bAln);
							auto compRes = UniqueKmerSetHelper::compareReadToSets(seq, uniqueKmersPerSet, extractingPars.compPars,
																																		hasher);
							compRes.clearWinnerIfOtherSetsHit();
							++currentMatchingCounts[compRes.winnerRevComp][compRes.winnerSet];
							++currentTotalInput;
							++currentTotalUnmapped;
						} else if (bAln.IsPrimaryAlignment() && bAln.IsMapped()) {
							++currentTotalInput;
						}
					}
				}
				{
					std::lock_guard<std::mutex> lock(matchingCountsMut);
					totalInput += currentTotalInput;
					totalUnmapped += currentTotalUnmapped;
					for(const auto & foundPerSet : currentMatchingCounts){
						for(const auto & count : foundPerSet.second){
							matchingCounts[foundPerSet.first][count.first] += count.second;
						}
					}
				}
			};
			njh::concurrent::runVoidFunctionThreaded(readInComp, numThreadsPerBam);
			{
				std::lock_guard<std::mutex> lock(outMut);
				for (const auto &name: names) {
					uint64_t total = matchingCounts[true][name] + matchingCounts[false][name];
					uint64_t totalForward = matchingCounts[false][name];
					out << bfs::basename(bam)
							<< "\t" << totalInput
							<< "\t" << totalUnmapped
							<< "\t" << name
							<< "\t" << total
							<< "\t" << (totalInput > 0 ? total / static_cast<long double>(totalInput) : 0)
							<< "\t" << (totalUnmapped > 0 ? total / static_cast<long double>(totalUnmapped) : 0)
							<< "\t" << totalForward
							<< "\t" << (total > 0 ? totalForward / static_cast<long double>(total) : 0)
							<< std::endl;
				}
			}
		}
	};
	njh::concurrent::runVoidFunctionThreaded(countBam, numThreads);

	return 0;
}

int kmerSetExpRunner::countingUniqKmersFromSetsBestSet(const njh::progutils::CmdArgs & inputCommands){
	uint32_t numThreads = 1;
	bfs::path countTable = "";

	UniqueKmerSetHelper::ProcessReadForExtractingPars extractingPars;
	extractingPars.compPars.hardCountOff = 10;
	extractingPars.compPars.fracCutOff = 0;
	extractingPars.compPars.kmerLengthForEntropyCalc_ = 2;
	extractingPars.compPars.entropyFilter_ = 1.20;
	extractingPars.compPars.initialExcludeHardCountOff = 60;
	extractingPars.compPars.initialExcludeFracCutOff = 0.25;
	bool looseCounting = false;
	OutOptions outOpts(bfs::path(""), ".tab.txt.gz");
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.pars_.ioOptions_.revComplMate_ = true;
	bool pairedEndSet = setUp.processReadInNames(njhseq::seqSetUp::pairedReadInFormatsAvailable_, false);

	auto singlesOption = SeqIOOptions();
	setUp.processJustReadInNames(singlesOption, njhseq::seqSetUp::singleInFormatsAvailable_, !pairedEndSet);
	singlesOption.lowerCaseBases_ = setUp.pars_.ioOptions_.lowerCaseBases_;
	singlesOption.removeGaps_ = setUp.pars_.ioOptions_.removeGaps_;
	singlesOption.includeWhiteSpaceInName_ = setUp.pars_.ioOptions_.includeWhiteSpaceInName_;
	setUp.setOption(countTable, "--kmerTable", "kmer set table, 1)set,2)kmer", true);
	setUp.setOption(extractingPars.compPars.sampleName, "--sampleName", "Name to add to output file", true);
	setUp.setOption(looseCounting, "--looseCounting", "by default will only count reads if there is only one set with reads, by adding this flag will count just best winner");

	setUp.setOption(extractingPars.compPars.includeRevComp, "--includeRevComp", "include Rev Comp of the input seqs");
//	setUp.setOption(extractingPars.compPars.entropyFilter_, "--entropyFilter", "entropy Filter_");
//	setUp.setOption(extractingPars.compPars.kmerLengthForEntropyCalc_, "--kmerLengthForEntropyCalc", "kmerLength For Entropy Calculation");
	setUp.setOption(extractingPars.compPars.fracCutOff, "--readKmerFracCutOff", "per read Kmer Frac Cut Off");
	setUp.setOption(extractingPars.compPars.hardCountOff, "--hardCountOff", "hard Count Off, do not count sets unless greater thant his number");
	setUp.setOption(extractingPars.compPars.excludeSetNames, "--excludeSetNames", "names of sets of unique kmers to ignore from the --kmerTable");

	setUp.setOption(numThreads, "--numThreads", "numThreads");


	setUp.processWritingOptions(outOpts);
	//setUp.processDirectoryOutputName("true");
	setUp.finishSetUp(std::cout);
	//setUp.startARunLog(setUp.pars_.directoryName_);
	njh::stopWatch watch;
	watch.setLapName("initial");
	auto klen = KmerCollectionIO::getKmerLength(countTable);
	extractingPars.compPars.klen = klen;


	watch.startNewLap("reading in unique kmer table");
	std::unordered_map<std::string, std::unordered_set<uint64_t>> uniqueKmersPerSet = KmerCollectionIO::readHashedCollections(countTable);
	if(setUp.pars_.verbose_){
		std::cout << watch.getLapName() << "\t" << watch.timeLapFormatted() <<std::endl;
	}

	VecStr names = getVectorOfMapKeys(uniqueKmersPerSet);


	OutputStream out(outOpts);
	std::function<void()> readInComp;
//	MultiSeqIO seqOut;
	uint64_t totalInputReads = 0;
	std::unordered_map<bool, std::unordered_map<std::string, uint64_t>> matchingCounts;
	std::mutex matchingCountsMut;
	if (!setUp.pars_.ioOptions_.firstName_.empty()) {
		SeqInput reader(setUp.pars_.ioOptions_);
		reader.openIn();
//		for(const auto & name : names){
//			auto seqOutOpts = SeqIOOptions::genPairedOutGz(name);
//			seqOutOpts.out_.overWriteFile_ = true;
//			seqOut.addReader(name, seqOutOpts);
//		}
//		readInComp = [&reader, &uniqueKmersPerSet, &uniqueKmersFoundPerSet,&kmersFoundPerSeq,&mut,&klen,&includeRevComp,&seqOut]() {
		readInComp = [&reader, &uniqueKmersPerSet,&extractingPars,
									&matchingCounts,&matchingCountsMut,
//									&seqOut,
									&totalInputReads , &looseCounting]() {

			SimpleKmerHash hasher;
			PairedRead pseq;

			std::unordered_map<bool, std::unordered_map<std::string, uint64_t>> currentMatchingCounts;
			uint64_t currentTotalInputReads = 0;
			while(reader.readNextReadLock(pseq)){
				auto compRes = UniqueKmerSetHelper::compareReadToSets(pseq, uniqueKmersPerSet, extractingPars.compPars, hasher);
				if (!looseCounting) {
					compRes.clearWinnerIfOtherSetsHit();
				}
				++currentMatchingCounts[compRes.winnerRevComp][compRes.winnerSet];
				++currentTotalInputReads;
//				if("undetermined" != compRes.winnerSet){
//					seqOut.openWrite(compRes.winnerSet, pseq);
//				}
			}
			{
				std::lock_guard<std::mutex> lock(matchingCountsMut);
				totalInputReads+= currentTotalInputReads;
				for(const auto & foundPerSet : currentMatchingCounts){
					for(const auto & count : foundPerSet.second){
						matchingCounts[foundPerSet.first][count.first] += count.second;
					}
				}
			}
		};
		njh::concurrent::runVoidFunctionThreaded(readInComp, numThreads);
	}
	if(!singlesOption.firstName_.empty()){
		SeqInput reader(singlesOption);
		reader.openIn();
//		for(const auto & name : names){
//			auto seqOutOpts = SeqIOOptions::genFastqOutGz(name);
//			seqOutOpts.out_.overWriteFile_ = true;
//			seqOut.addReader(name, seqOutOpts);
//		}
//		readInComp = [&reader, &uniqueKmersPerSet, &uniqueKmersFoundPerSet,&kmersFoundPerSeq,&mut,&klen,&includeRevComp,&seqOut]() {
		readInComp = [&reader, &uniqueKmersPerSet,&extractingPars,
						&matchingCounts,&matchingCountsMut,
//						&seqOut,
						&totalInputReads, &looseCounting]() {

			SimpleKmerHash hasher;
			seqInfo seq;

			std::unordered_map<bool, std::unordered_map<std::string, uint64_t>> currentMatchingCounts;
			uint64_t currentTotalInputReads = 0;
			while(reader.readNextReadLock(seq)){
				auto compRes = UniqueKmerSetHelper::compareReadToSets(seq, uniqueKmersPerSet, extractingPars.compPars,
																															hasher);
				if (!looseCounting) {
					compRes.clearWinnerIfOtherSetsHit();
				}
				++currentTotalInputReads;
				++currentMatchingCounts[compRes.winnerRevComp][compRes.winnerSet];
//				if("undetermined" != compRes.winnerSet){
//					seqOut.openWrite(compRes.winnerSet, seq);
//				}
			}
			{
				std::lock_guard<std::mutex> lock(matchingCountsMut);
				totalInputReads+= currentTotalInputReads;
				for(const auto & foundPerSet : currentMatchingCounts){
					for(const auto & count : foundPerSet.second){
						matchingCounts[foundPerSet.first][count.first] += count.second;
					}
				}
			}
		};
		njh::concurrent::runVoidFunctionThreaded(readInComp, numThreads);
	}

	out << "sample\ttotalInputReads\tset\treads\treadsFrac\treadsInForward\tforwardFrac" << std::endl;
	njh::sort(names);
	for(const auto & name : names){
		uint64_t total = matchingCounts[true][name] + matchingCounts[false][name];
		uint64_t totalForward = matchingCounts[false][name];
		out << extractingPars.compPars.sampleName
		<< "\t" << totalInputReads
				<< "\t" << name
				<< "\t" << total
				<< "\t" << (totalInputReads > 0 ? total/static_cast<long double>(totalInputReads) : 0)
				<< "\t" << totalForward
				<< "\t" << (total > 0 ? totalForward/static_cast<long double>(total) : 0)
				<< std::endl;
	}
	return 0;
}



int kmerSetExpRunner::findKmersUniqueAndConservedToSet(const njh::progutils::CmdArgs & inputCommands){
	KmerGatherer::KmerGathererPars countPars;
	bfs::path seqSetTableFnp = "";
	bool fasta = false;
	bool noFilters = false;
	bool exportNonUniqueKmers = false;
	uint32_t topThreadCount = 1;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.description_ = "Get kmers that appear within the sequences supplied and are unique to that set compared to the other sets";
	setUp.setOption(fasta, "--fasta", "file contains fasta files instead of 2bit files");
	setUp.setOption(noFilters, "--noFilters", "Don't do filtering on kmers sets");
	setUp.setOption(exportNonUniqueKmers, "--exportNonUniqueKmers", "export Non Unique Kmers");
	setUp.setOption(topThreadCount, "--topThreadCount", "number of sets to gather at once, --numThreads is split between them");

	if(noFilters){
		countPars.entropyFilter_ = 0;
		countPars.allowableCharacters_ = njh::genSetOfAnsiPrintable();
	}
	std::string columnsHelp = njh::pasteAsStr("1)set,2)", fasta? "fasta": "2bit");
	setUp.setOption(seqSetTableFnp, "--seqSetTableFnp", "Seq Set Table, 2 columns, " + columnsHelp, true);
	countPars.setOptions(setUp);
	if(0 == topThreadCount){
		setUp.failed_ = true;
		setUp.addWarning("--topThreadCount must be at least 1");
	}
	setUp.processDirectoryOutputName(njh::pasteAsStr(bfs::basename(seqSetTableFnp), "_TODAY"), true);
	setUp.finishSetUp(std::cout);

	setUp.startARunLog(setUp.pars_.directoryName_);

	VecStr columnRequired{"set", "2bit"};
	if(fasta){
		columnRequired = {"set", "fasta"};
	}
	KmerGatherer::KmerGathererPars perSetPars = countPars;
	perSetPars.numThreads_ = std::max<uint32_t>(1, countPars.numThreads_ / topThreadCount);

	std::unordered_map<std::string, std::set<bfs::path>> seqFilesForSet;
	{
		table input(seqSetTableFnp, "\t", true);
		input.checkForColumnsThrow(columnRequired, __PRETTY_FUNCTION__);
		for(const auto & row : input){
			seqFilesForSet[row[input.getColPos("set")]].emplace(row[input.getColPos(columnRequired[1])]);
		}
	}



	setUp.rLog_.setCurrentLapName("initial");

	std::vector<bfs::path> seqFiles;
	for(const auto & seqSet : seqFilesForSet){
		for(const auto & fnp : seqSet.second){
			seqFiles.emplace_back(fnp);
		}
	}
	std::map<std::string, std::set<uint64_t>> kmersPerSet;

	{
		setUp.rLog_.logCurrentTime("count_and_condense_per_set");
		setUp.rLog_.runLogFile_.flush();
		njh::concurrent::LockableQueue<std::string> seqSetNamesQueue(getVectorOfMapKeys(seqFilesForSet));
		for(const auto & name : seqFilesForSet){
			kmersPerSet[name.first] = std::set<uint64_t>{};
		}
		OutputStream countOut_per_file(njh::files::make_path(setUp.pars_.directoryName_, "counts_per_input.tsv.gz"));
		countOut_per_file << "set\tcount" << std::endl;
		std::mutex out_all_kmer_count_mutex;
		std::function<void()> condenseKmers = [&seqSetNamesQueue,&seqFilesForSet,&kmersPerSet,
			&perSetPars,&countOut_per_file, &fasta,
			&out_all_kmer_count_mutex](){
			std::string name;
			while(seqSetNamesQueue.getVal(name)){
				std::unordered_map<std::string, std::set<uint64_t>> allKmers_per_set;
				KmerGatherer kGather(perSetPars);
				if (!fasta) {
					allKmers_per_set = kGather.getFilteredHashedKmerCollectionsFromTwobits(std::vector(seqFilesForSet.at(name).begin(), seqFilesForSet.at(name).end()));
				} else {
					allKmers_per_set = kGather.getFilteredHashedKmerCollectionsFromFastas(std::vector(seqFilesForSet.at(name).begin(), seqFilesForSet.at(name).end()) );
				}
				{
					std::lock_guard<std::mutex> lock(out_all_kmer_count_mutex);
					for(const auto & allKmers_for_file : allKmers_per_set){
						countOut_per_file << allKmers_for_file.first << "\t" << allKmers_for_file.second.size() << std::endl;
					}
				}
				std::unordered_map<uint64_t, uint32_t> counts;
				for(const auto & fnp : seqFilesForSet.at(name)){
					for (const auto k : allKmers_per_set.at(fnp.string())) {
						++counts[k];
					}
				}
				//right now looking for perfectly conserved but could expand to doing a percentage
				for (const auto & count : counts) {
					if (count.second == seqFilesForSet.at(name).size()) {
						kmersPerSet.at(name).emplace(count.first);
					}
				}
			}
		};
		njh::concurrent::runVoidFunctionThreaded(condenseKmers, topThreadCount);

		OutputStream countOut_per_set(njh::files::make_path(setUp.pars_.directoryName_, "counts_per_set.tsv.gz"));
		countOut_per_set << "set\tcount" << std::endl;
		for(const auto & allKmers_per_set : kmersPerSet){
			countOut_per_set << allKmers_per_set.first << "\t" << allKmers_per_set.second.size() << std::endl;
		}
	}
	setUp.rLog_.logCurrentTime("compare");
	setUp.rLog_.runLogFile_.flush();
	auto uniqueKmersFinal = KmerGatherer::filterToKmersUniqueToEachCollection(kmersPerSet);
	std::set<uint64_t> nonUniqueKmers;
	if(exportNonUniqueKmers){
		nonUniqueKmers = KmerGatherer::getKmersSharedBetweenCollections(kmersPerSet);
	}

	std::map<std::string, uint64_t> kmerPerSetCounts;
	OutputStream out(njh::files::make_path(setUp.pars_.directoryName_, "kmerSets.tsv.gz"));
	KmerCollectionIO::writeHashedCollections(out, uniqueKmersFinal);
	for(const auto & kmersForSet : uniqueKmersFinal){
		kmerPerSetCounts[kmersForSet.first] = kmersForSet.second.size();
	}
	if(exportNonUniqueKmers){
		std::string nonUniqueRegionName = "NON_UNIQUE";
		OutputStream nonUniquOut(njh::files::make_path(setUp.pars_.directoryName_, "nonUniqueKmers_kmerSets.tsv.gz"));
		kmerPerSetCounts[nonUniqueRegionName] = nonUniqueKmers.size();
		KmerCollectionIO::writeHashedCollection(nonUniquOut, nonUniqueRegionName, nonUniqueKmers);
	}
	OutputStream countOut(njh::files::make_path(setUp.pars_.directoryName_, "counts.tsv.gz"));
	countOut << "set\tcount" << std::endl;
	for(const auto & kmerPerSetCount : kmerPerSetCounts){
		countOut << kmerPerSetCount.first << "\t" << kmerPerSetCount.second << std::endl;
	}
	return 0;
}


int kmerSetExpRunner::findUniqKmersBetweenSeqSetsMulti(const njh::progutils::CmdArgs & inputCommands){
	KmerGatherer::KmerGathererPars countPars;
	bfs::path seqSetTableFnp = "";
	bfs::path seqSetSuppFastaTableFnp = "";
	bool fasta = false;
	bool noFilters = false;
	bool exportNonUniqueKmers = false;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.description_ = "Get kmers that appear within the sequences supplied and are unique to that set compared to the other sets";
	setUp.setOption(fasta, "--fasta", "file contains fasta files instead of 2bit files");
	setUp.setOption(noFilters, "--noFilters", "Don't do filtering on kmers sets");
	setUp.setOption(exportNonUniqueKmers, "--exportNonUniqueKmers", "export Non Unique Kmers");

	if(noFilters){
		countPars.entropyFilter_ = 0;
		countPars.allowableCharacters_ = njh::genSetOfAnsiPrintable();
	}
	std::string columnsHelp = njh::pasteAsStr("1)set,2)", fasta? "fasta": "2bit");
	setUp.setOption(seqSetTableFnp, "--seqSetTableFnp", "Seq Set Table, 2 columns, " + columnsHelp, true);
	setUp.setOption(seqSetSuppFastaTableFnp, "--seqSetSuppFastaTableFnp", "Seq Set Table supplement small fasta files, 2 columns, 1)set,2)fasta");
	countPars.setOptions(setUp);
	setUp.processDirectoryOutputName(njh::pasteAsStr(bfs::basename(seqSetTableFnp), "_TODAY"), true);
	setUp.finishSetUp(std::cout);

	setUp.startARunLog(setUp.pars_.directoryName_);

	VecStr columnRequired{"set", "2bit"};
	if(fasta){
		columnRequired = {"set", "fasta"};
	}
	KmerGatherer kGather(countPars);

	std::unordered_map<std::string, std::set<std::string>> seqFilesForSet;
	{
		table input(seqSetTableFnp, "\t", true);
		input.checkForColumnsThrow(columnRequired, __PRETTY_FUNCTION__);
		for(const auto & row : input){
			seqFilesForSet[row[input.getColPos("set")]].emplace(row[input.getColPos(columnRequired[1])]);
		}
	}

	table inputFastaSupp;
	if(!seqSetSuppFastaTableFnp.empty()){
		inputFastaSupp= table(seqSetSuppFastaTableFnp, "\t", true);
		inputFastaSupp.checkForColumnsThrow(VecStr{"set", "fasta"}, __PRETTY_FUNCTION__);
	}

	setUp.rLog_.setCurrentLapName("initial");

	std::vector<bfs::path> seqFiles;
	for(const auto & seqSet : seqFilesForSet){
		for(const auto & fnp : seqSet.second){
			seqFiles.emplace_back(fnp);
		}
	}
	std::map<std::string, std::set<uint64_t>> kmersPerSet;

	{
		setUp.rLog_.logCurrentTime("count_all");
		setUp.rLog_.runLogFile_.flush();
		std::unordered_map<std::string, std::set<uint64_t>> allKmers;
		if (!fasta) {
			allKmers = kGather.getFilteredHashedKmerCollectionsFromTwobits(seqFiles);
		} else {
			allKmers = kGather.getFilteredHashedKmerCollectionsFromFastas(seqFiles);
		}
		if(!seqSetSuppFastaTableFnp.empty()){
			std::vector<bfs::path> fastaFiles;
			for(const auto & row : inputFastaSupp){
				seqFilesForSet[row[inputFastaSupp.getColPos("set")]].emplace(row[inputFastaSupp.getColPos("fasta")]);
				fastaFiles.emplace_back(row[inputFastaSupp.getColPos("fasta")]);
			}
			auto suppKmers = kGather.getFilteredHashedKmerCollectionsFromFastas(fastaFiles);
			for(const auto & supp : suppKmers){
				allKmers[supp.first].insert(supp.second.begin(), supp.second.end());
			}
		}
		setUp.rLog_.logCurrentTime("condense");
		setUp.rLog_.runLogFile_.flush();
		kmersPerSet = KmerGatherer::condenseCollections(allKmers, seqFilesForSet);
	}
	setUp.rLog_.logCurrentTime("compare");
	setUp.rLog_.runLogFile_.flush();
	auto uniqueKmersFinal = KmerGatherer::filterToKmersUniqueToEachCollection(kmersPerSet);
	std::set<uint64_t> nonUniqueKmers;
	if(exportNonUniqueKmers){
		nonUniqueKmers = KmerGatherer::getKmersSharedBetweenCollections(kmersPerSet);
	}

	std::map<std::string, uint64_t> kmerPerSetCounts;
	OutputStream out(njh::files::make_path(setUp.pars_.directoryName_, "kmerSets.tsv.gz"));
	KmerCollectionIO::writeHashedCollections(out, uniqueKmersFinal);
	for(const auto & kmersForSet : uniqueKmersFinal){
		kmerPerSetCounts[kmersForSet.first] = kmersForSet.second.size();
	}
	if(exportNonUniqueKmers){
		std::string nonUniqueRegionName = "NON_UNIQUE";
		OutputStream nonUniquOut(njh::files::make_path(setUp.pars_.directoryName_, "nonUniqueKmers_kmerSets.tsv.gz"));
		kmerPerSetCounts[nonUniqueRegionName] = nonUniqueKmers.size();
		KmerCollectionIO::writeHashedCollection(nonUniquOut, nonUniqueRegionName, nonUniqueKmers);
	}
	OutputStream countOut(njh::files::make_path(setUp.pars_.directoryName_, "counts.tsv.gz"));
	countOut << "set\tcount" << std::endl;
	for(const auto & kmerPerSetCount : kmerPerSetCounts){
		countOut << kmerPerSetCount.first << "\t" << kmerPerSetCount.second << std::endl;
	}
	return 0;
}

int kmerSetExpRunner::findKmersInSets(const njh::progutils::CmdArgs & inputCommands){
	KmerGatherer::KmerGathererPars countPars;
	std::vector<bfs::path> seqSetFnps;
	uint32_t minOccurrences = 1;
	bfs::path countTable;
	OutOptions outOpts(bfs::path(""), ".tsv.gz");
	seqSetUp setUp(inputCommands);
	setUp.description_ = "Take a sequence file and report for each kmer stepped along the seq if that kmer can be found within a set of other fasta files";
	setUp.processVerbose();
	setUp.processDebug();
	bool seqSetFnpsSet = setUp.setOption(seqSetFnps, "--seqSetFnps", "seqSetFnps", false);
	setUp.setOption(minOccurrences, "--minOccurrences", "minOccurrences");

	setUp.setOption(countTable, "--kmerTable", "countTable, 1)set,2)kmer", !seqSetFnpsSet);

	countPars.setOptions(setUp);
	setUp.processReadInNames(true);
	setUp.processWritingOptions(outOpts);
	// setUp.processDirectoryOutputName(true);
	setUp.finishSetUp(std::cout);
	// setUp.startARunLog(setUp.pars_.directoryName_);

	std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>> allCounts;
	// setUp.rLog_.setCurrentLapName("initial");

	if (countTable.empty()) {
		KmerGatherer kGather(countPars);
		for(const auto & seqSetFnp : seqSetFnps){
			std::string genomeFnp1Base = bfs::basename(njh::files::removeExtension(seqSetFnp));
			// setUp.rLog_.logCurrentTime("genome1_counting-" + genomeFnp1Base);
			allCounts[genomeFnp1Base] = kGather.countKmers(seqSetFnp);
		}
	}	else {
		TableReader counts_table_reader(TableIOOpts::genTabFileIn(countTable, false));
		VecStr row;
		if (counts_table_reader.header_.nCol() != 2) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "expected 2 columns for " << countTable << " but found " << counts_table_reader.header_.nCol() << "\n";
			throw std::runtime_error{ss.str()};
		}
		while (counts_table_reader.getNextRow(row)) {
			++allCounts[row[0]][row[1]];
			countPars.kmerLength_ = row[1].size();
		}
		if (allCounts.empty()) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "no kmers read from " << countTable << "\n";
			throw std::runtime_error{ss.str()};
		}
	}


	seqInfo seq;
	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn()	;
	// OutputStream kmerClassifierOut(njh::files::make_path(setUp.pars_.directoryName_, "seqsKmerClassified.tab.txt.gz"));
	OutputStream kmerClassifierOut(outOpts);
	kmerClassifierOut << "name\tpos\tkmer\tfoundIn" << std::endl;
	while(reader.readNextRead(seq)){
		if(len(seq) >= countPars.kmerLength_){
			for(auto pos = 0; pos < (len(seq) + 1 - countPars.kmerLength_); ++pos){
				auto currentK = seq.seq_.substr(pos, countPars.kmerLength_);
				std::set<std::string> setsFounds;
				for(auto & seqSet : allCounts){
					if(njh::in(currentK, seqSet.second) && seqSet.second[currentK] >= minOccurrences){
						setsFounds.emplace(seqSet.first);
					}
				}
				if (!countPars.noRevComp_) {
					auto revComp_currentK = seqUtil::reverseComplement(currentK, "DNA");
					for(auto & seqSet : allCounts){
						if(njh::in(revComp_currentK, seqSet.second) && seqSet.second[revComp_currentK] >= minOccurrences){
							setsFounds.emplace(njh::pasteAsStr(seqSet.first, "-revComp"));
						}
					}
				}
				kmerClassifierOut
								<< seq.name_
								<< "\t" << pos
								<< "\t" << currentK
								<< "\t";
				if (setsFounds.empty()) {
					kmerClassifierOut << "none";
				} else {
					kmerClassifierOut << njh::conToStr(setsFounds, ",");
				}
				kmerClassifierOut << std::endl;
			}
		}
	}

	// setUp.rLog_.logCurrentTime("end");

	return 0;
}





}  //namespace njhseq

