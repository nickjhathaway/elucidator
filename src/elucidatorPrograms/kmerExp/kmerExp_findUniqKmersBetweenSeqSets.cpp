/*
 * kmerExp_findUniqKmersBetweenSeqSets.cpp
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

#include "kmerExp.hpp"
#include <njhseq/objects/kmer/KmerUtils.hpp>
#include <njhseq/objects/kmer.h>

#include "elucidator/simulation.h"

#include "elucidator/objects/MiscUtility/GenomeSeqSearch.hpp"

#include <njhseq/objects/dataContainers/tables/TableReader.hpp>
#include <njhseq/IO/SeqIO/MultiSeqOutCache.hpp>
#include "elucidator/helpers/UniqueKmerSetHelper.hpp"

namespace njhseq {




int kmerExpRunner::filterUniqueKmerSetForEntropy(const njh::progutils::CmdArgs & inputCommands){
	bfs::path countTable = "";
	double entropyCutOff = 1;
	OutOptions outOpts(bfs::path(""), ".tab.txt.gz");
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(countTable, "--countTable", "countTable, 1)set,2)kmer", true);
	setUp.setOption(entropyCutOff, "--entropyCutOff", "entropy Cut Off to keep kmers (exclusive)");

	setUp.processWritingOptions(outOpts);
	//setUp.processDirectoryOutputName("true");
	setUp.finishSetUp(std::cout);
	//setUp.startARunLog(setUp.pars_.directoryName_);
	njh::stopWatch watch;
	watch.setLapName("initial");
	watch.startNewLap("reading in unique kmer table");
	OutputStream out(outOpts);

	{
		TableReader uniqKmers(TableIOOpts::genTabFileIn(countTable, false));
		if(uniqKmers.header_.nCol() != 2){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "need to have 2 columns" << "\n";
			throw std::runtime_error{ss.str()};
		}
		VecStr row;
		while(uniqKmers.getNextRow(row)){
			DNABaseCounter counter;
			counter.increase(row[1]);
			if(counter.computeEntrophy() >entropyCutOff){
				out << row[0]
					 << "\t" << row[1]
					 << "\t" << counter.computeEntrophy()
					 << std::endl;
			}
		}
	}

	return 0;
}

int kmerExpRunner::countingUniqKmersFromSetsInUnmappedAlns(const njh::progutils::CmdArgs & inputCommands){
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
	setUp.setOption(includeRevComp, "--includeRevComp", "include Rev Comp of the input seqs");




	setUp.processWritingOptions(outOpts);
	//setUp.processDirectoryOutputName("true");
	setUp.finishSetUp(std::cout);
	//setUp.startARunLog(setUp.pars_.directoryName_);
	njh::stopWatch watch;
	watch.setLapName("initial");

	auto klen =	UniqueKmerSetHelper::getKmerLenFromUniqueKmerTable(countTable);

	std::mutex mut;

	watch.startNewLap("reading in unique kmer table");
	std::unordered_map<std::string, std::unordered_set<uint64_t>> uniqueKmersPerSet = UniqueKmerSetHelper::readInUniqueKmerTablePerSet(countTable);
	if(setUp.pars_.verbose_){
		std::cout << watch.getLapName() << "\t" << watch.timeLapFormatted() <<std::endl;
	}
	std::unordered_map<std::string, std::unordered_map<uint64_t, uint64_t>> uniqueKmersFoundPerSet;
	std::unordered_map<std::string, std::vector<uint64_t>> kmersFoundPerSeq;
	OutputStream out(outOpts);
	VecStr names = getVectorOfMapKeys(uniqueKmersPerSet);
//	MultiSeqIO seqOut;
//	if(setUp.pars_.debug_){
//		for (const auto &name : names) {
//			auto seqOutOpts = SeqIOOptions::genPairedOutGz(njh::pasteAsStr(basename(setUp.pars_.ioOptions_.firstName_), "-", name, "-paired"));
//			seqOutOpts.out_.overWriteFile_ = true;
//			seqOut.addReader(name + "-paired", seqOutOpts);
//
//			auto seqOutOpts = SeqIOOptions::genFastqOutGz(njh::pasteAsStr(basename(setUp.pars_.ioOptions_.firstName_), "-", name, "-single"));
//			seqOutOpts.out_.overWriteFile_ = true;
//			seqOut.addReader(name+ "-single", seqOutOpts);
//		}
//	}
	struct LockedBamReader {

		LockedBamReader(const bfs::path & bamFnp){
			bReader_.Open(bamFnp.string());
			checkBamOpenThrow(bReader_, bamFnp);
		}
		BamTools::BamReader bReader_;

		std::mutex mut_;

		bool readNextAlnLockFree(BamTools::BamAlignment & bAln){
			return bReader_.GetNextAlignment(bAln);
		}

		bool readNextAlnLock(BamTools::BamAlignment & bAln){
			std::lock_guard<std::mutex> lock(mut_);
			return bReader_.GetNextAlignment(bAln);
		}

		bool readNextAlnLockBatch(std::vector<BamTools::BamAlignment> & alns, uint32_t batchSize){
			std::lock_guard<std::mutex> lock(mut_);
			alns.clear();
			uint32_t count = 0;
			BamTools::BamAlignment baln;
			while(count < batchSize && readNextAlnLockFree(baln)){
				alns.emplace_back(baln);
				++count;
			}
			return !alns.empty();
		}
	};


	LockedBamReader bamReader(setUp.pars_.ioOptions_.firstName_);

	std::function<void()> readInComp = [&bamReader, &uniqueKmersPerSet, &uniqueKmersFoundPerSet,&kmersFoundPerSeq,&mut,&klen,&includeRevComp]() {
		SimpleKmerHash hasher;
		std::unordered_map<std::string, std::unordered_map<uint64_t, uint64_t>> uniqueKmersFoundPerSetCurrent;
		std::unordered_map<std::string, std::vector<uint64_t>> kmersFoundPerSeqCurrent;
		std::vector<BamTools::BamAlignment> alns;
		while(bamReader.readNextAlnLockBatch(alns, 10000)){
			for(const auto & bAln : alns){
				if(bAln.IsPrimaryAlignment() && !bAln.IsMapped()){

					std::unordered_map<uint64_t, uint64_t> hashedInputKmers;

					//since unmapped reads are modifed can just take the query bases and don't have to worry about revComping
					if(len(bAln.QueryBases) > klen){
						for(uint32_t pos = 0; pos < len(bAln.QueryBases) - klen + 1; ++pos){
							auto hash = hasher.hash(bAln.QueryBases.substr(pos, klen));
							++hashedInputKmers[hash];
						}
					}
					if(includeRevComp){
						//bAln.QueryBases = seqUtil::reverseComplement(bAln.QueryBases, "DNA");
						if(len(bAln.QueryBases) > klen){
							for(uint32_t pos = 0; pos < len(bAln.QueryBases) - klen + 1; ++pos){
								auto hash = hasher.revCompHash(bAln.QueryBases.substr(pos, klen));
								++hashedInputKmers[hash];
							}
						}
					}
					std::unordered_map<std::string, uint32_t> foundPerSet;
					for(const auto & hashedKmer : hashedInputKmers){
						for(const auto & uniqueKmers : uniqueKmersPerSet){
							if(njh::in(hashedKmer.first, uniqueKmers.second)){
								uniqueKmersFoundPerSetCurrent[uniqueKmers.first][hashedKmer.first]+= hashedKmer.second;
								++foundPerSet[uniqueKmers.first];
								break;
							}
						}
					}
					for(const auto & found : foundPerSet){
						kmersFoundPerSeqCurrent[found.first].emplace_back(found.second);
		//					seqOut.openWrite(found.first, pseq);
					}
				}
			}
		}
		{
			std::lock_guard<std::mutex> lock(mut);
			for(const auto & foundPerSet : uniqueKmersFoundPerSetCurrent){
				for(const auto & count : foundPerSet.second){
					uniqueKmersFoundPerSet[foundPerSet.first][count.first] += count.second;
				}
			}
			for(const auto & foundPerSeq : kmersFoundPerSeqCurrent){
				addOtherVec(kmersFoundPerSeq[foundPerSeq.first], foundPerSeq.second);
			}
		}
	};
	njh::concurrent::runVoidFunctionThreaded(readInComp, numThreads);
	out << "sample\tset\treads\tmeanPerRead\ttotal\ttotal1\ttotal2\tunique\tunique1\tunique2\tuniqueInSet\tmeanOccurrence\tmeanOccurrence1\tmeanOccurrence2\tfracUniqFound\tfracUniqFound1\tfracUniqFound2" << std::endl;
	njh::sort(names);
	for(const auto & name : names){
		uint64_t total = 0;
		uint64_t totalMoreThanOnce = 0;
		uint64_t totalMoreThanTwice = 0;
		std::unordered_set<uint64_t> moreThanOnce;
		std::unordered_set<uint64_t> moreThanTwice;
		for(const auto & countPerSet : uniqueKmersFoundPerSet[name]){
			total += countPerSet.second;
			if(countPerSet.second>1){
				totalMoreThanOnce += countPerSet.second;
				moreThanOnce.emplace(countPerSet.first);
			}
			if(countPerSet.second>2){
				totalMoreThanTwice += countPerSet.second;
				moreThanTwice.emplace(countPerSet.first);
			}
		}
		uint64_t readCount = kmersFoundPerSeq[name].size();
		long double meanPerSeq = vectorMean(kmersFoundPerSeq[name]);
		auto occMean =  static_cast<long double>(total)/uniqueKmersFoundPerSet[name].size();
		auto occMean1 = static_cast<long double>(totalMoreThanOnce)/moreThanOnce.size();
		auto occMean2 = static_cast<long double>(totalMoreThanTwice)/moreThanTwice.size();

		out << sampleName
				<< "\t" << name
				<< "\t" << readCount
				<< "\t" << meanPerSeq
				<< "\t" << total
				<< "\t" << totalMoreThanOnce
				<< "\t" << totalMoreThanTwice
				<< "\t" << uniqueKmersFoundPerSet[name].size()
				<< "\t" << moreThanOnce.size()
				<< "\t" << moreThanTwice.size()
				<< "\t" << uniqueKmersPerSet[name].size()
				<< "\t" << occMean
				<< "\t" << occMean1
				<< "\t" << occMean2
				<< "\t" << static_cast<long double>(uniqueKmersFoundPerSet[name].size())/uniqueKmersPerSet[name].size()
				<< "\t" << static_cast<long double>(moreThanOnce.size())/uniqueKmersPerSet[name].size()
				<< "\t" << static_cast<long double>(moreThanTwice.size())/uniqueKmersPerSet[name].size()
				<< std::endl;
	}

	return 0;
}


int kmerExpRunner::countingUniqKmersFromSets(const njh::progutils::CmdArgs & inputCommands){
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
	setUp.setOption(includeRevComp, "--includeRevComp", "include Rev Comp of the input seqs");




	setUp.processWritingOptions(outOpts);
	//setUp.processDirectoryOutputName("true");
	setUp.finishSetUp(std::cout);
	//setUp.startARunLog(setUp.pars_.directoryName_);
	njh::stopWatch watch;
	watch.setLapName("initial");
	auto klen =	UniqueKmerSetHelper::getKmerLenFromUniqueKmerTable(countTable);

	std::mutex mut;

	watch.startNewLap("reading in unique kmer table");
	std::unordered_map<std::string, std::unordered_set<uint64_t>> uniqueKmersPerSet = UniqueKmerSetHelper::readInUniqueKmerTablePerSet(countTable);
	if(setUp.pars_.verbose_){
		std::cout << watch.getLapName() << "\t" << watch.timeLapFormatted() <<std::endl;
	}
	std::unordered_map<std::string, std::unordered_map<uint64_t, uint64_t>> uniqueKmersFoundPerSet;
	std::unordered_map<std::string, std::vector<uint64_t>> kmersFoundPerSeq;
	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();

	OutputStream out(outOpts);
	std::function<void()> readInComp;
	VecStr names = getVectorOfMapKeys(uniqueKmersPerSet);
//	MultiSeqIO seqOut;


	if (setUp.pars_.ioOptions_.isPairedIn()) {
//		for(const auto & name : names){
//			auto seqOutOpts = SeqIOOptions::genPairedOutGz(name);
//			seqOutOpts.out_.overWriteFile_ = true;
//			seqOut.addReader(name, seqOutOpts);
//		}
//		readInComp = [&reader, &uniqueKmersPerSet, &uniqueKmersFoundPerSet,&kmersFoundPerSeq,&mut,&klen,&includeRevComp,&seqOut]() {
		readInComp = [&reader, &uniqueKmersPerSet, &uniqueKmersFoundPerSet,&kmersFoundPerSeq,&mut,&klen,&includeRevComp]() {

			SimpleKmerHash hasher;
			PairedRead pseq;
			std::unordered_map<std::string, std::unordered_map<uint64_t, uint64_t>> uniqueKmersFoundPerSetCurrent;
			std::unordered_map<std::string, std::vector<uint64_t>> kmersFoundPerSeqCurrent;

			while(reader.readNextReadLock(pseq)){
				std::unordered_map<uint64_t, uint64_t> hashedInputKmers;

				std::unordered_set<uint64_t> hashedInputKmersInFirstMate;
				if(len(pseq.seqBase_.seq_) > klen){
					for(uint32_t pos = 0; pos < len(pseq.seqBase_.seq_) - klen + 1; ++pos){
						auto hash = hasher.hash(pseq.seqBase_.seq_.substr(pos, klen));
						hashedInputKmersInFirstMate.emplace(hash);
						++hashedInputKmers[hash];
					}
				}
				if(len(pseq.mateSeqBase_.seq_) > klen){
					for(uint32_t pos = 0; pos < len(pseq.mateSeqBase_.seq_) - klen + 1; ++pos){
						auto hash = hasher.hash(pseq.mateSeqBase_.seq_.substr(pos, klen));
						if(!njh::in(hash, hashedInputKmersInFirstMate)){
							++hashedInputKmers[hash];
						}
					}
				}
				if(includeRevComp){
					//pseq.seqBase_.seq_ = seqUtil::reverseComplement(pseq.seqBase_.seq_, "DNA");
					if(len(pseq.seqBase_.seq_) > klen){
						for(uint32_t pos = 0; pos < len(pseq.seqBase_.seq_) - klen + 1; ++pos){
							auto hash = hasher.revCompHash(pseq.seqBase_.seq_.substr(pos, klen));
							hashedInputKmersInFirstMate.emplace(hash);
							++hashedInputKmers[hash];
						}
					}
					//pseq.mateSeqBase_.seq_ = seqUtil::reverseComplement(pseq.mateSeqBase_.seq_, "DNA");
					if(len(pseq.mateSeqBase_.seq_) > klen){
						for(uint32_t pos = 0; pos < len(pseq.mateSeqBase_.seq_) - klen + 1; ++pos){
							auto hash = hasher.revCompHash(pseq.mateSeqBase_.seq_.substr(pos, klen));
							if(!njh::in(hash, hashedInputKmersInFirstMate)){
								++hashedInputKmers[hash];
							}
						}
					}
				}
				std::unordered_map<std::string, uint32_t> foundPerSet;
				for(const auto & hashedKmer : hashedInputKmers){
					for(const auto & uniqueKmers : uniqueKmersPerSet){
						if(njh::in(hashedKmer.first, uniqueKmers.second)){
							uniqueKmersFoundPerSetCurrent[uniqueKmers.first][hashedKmer.first]+= hashedKmer.second;
							++foundPerSet[uniqueKmers.first];
							break;
						}
					}
				}
				for(const auto & found : foundPerSet){
					kmersFoundPerSeqCurrent[found.first].emplace_back(found.second);
//					seqOut.openWrite(found.first, pseq);
				}
			}
			{
				std::lock_guard<std::mutex> lock(mut);
				for(const auto & foundPerSet : uniqueKmersFoundPerSetCurrent){
					for(const auto & count : foundPerSet.second){
						uniqueKmersFoundPerSet[foundPerSet.first][count.first] += count.second;
					}
				}
				for(const auto & foundPerSeq : kmersFoundPerSeqCurrent){
					addOtherVec(kmersFoundPerSeq[foundPerSeq.first], foundPerSeq.second);
				}
			}
		};
	} else {
		readInComp = [&reader, &uniqueKmersPerSet, &uniqueKmersFoundPerSet,&kmersFoundPerSeq,&mut,&klen,&includeRevComp]() {
			SimpleKmerHash hasher;
			seqInfo seq;
			std::unordered_map<std::string, std::unordered_map<uint64_t, uint64_t>> uniqueKmersFoundPerSetCurrent;
			std::unordered_map<std::string, std::vector<uint64_t>> kmersFoundPerSeqCurrent;

			while(reader.readNextReadLock(seq)){
				std::unordered_map<uint64_t, uint64_t> hashedInputKmers;
				if(len(seq.seq_) > klen){
					for(uint32_t pos = 0; pos < len(seq.seq_) - klen + 1; ++pos){
						auto hash = hasher.hash(seq.seq_.substr(pos, klen));
						++hashedInputKmers[hash];
					}
				}
				if(includeRevComp){
					//seq.seq_ = seqUtil::reverseComplement(seq.seq_, "DNA");
					if(len(seq.seq_) > klen){
						for(uint32_t pos = 0; pos < len(seq.seq_) - klen + 1; ++pos){
							auto hash = hasher.revCompHash(seq.seq_.substr(pos, klen));
							++hashedInputKmers[hash];
						}
					}
				}
				std::unordered_map<std::string, uint32_t> foundPerSet;
				for(const auto & hashedKmer : hashedInputKmers){
					for(const auto & uniqueKmers : uniqueKmersPerSet){
						if(njh::in(hashedKmer.first, uniqueKmers.second)){
							uniqueKmersFoundPerSetCurrent[uniqueKmers.first][hashedKmer.first]+= hashedKmer.second;
							++foundPerSet[uniqueKmers.first];
							break;
						}
					}
				}
				for(const auto & found : foundPerSet){
					kmersFoundPerSeqCurrent[found.first].emplace_back(found.second);
//					seqOut.openWrite(found.first, pseq);
				}
			}
			{
				std::lock_guard<std::mutex> lock(mut);
				for(const auto & foundPerSet : uniqueKmersFoundPerSetCurrent){
					for(const auto & count : foundPerSet.second){
						uniqueKmersFoundPerSet[foundPerSet.first][count.first] += count.second;
					}
				}
				for(const auto & foundPerSeq : kmersFoundPerSeqCurrent){
					addOtherVec(kmersFoundPerSeq[foundPerSeq.first], foundPerSeq.second);
				}
			}
		};
	}
	njh::concurrent::runVoidFunctionThreaded(readInComp, numThreads);
	out << "sample\tset\treads\tmeanPerRead\ttotal\ttotal1\ttotal2\tunique\tunique1\tunique2\tuniqueInSet\tmeanOccurrence\tmeanOccurrence1\tmeanOccurence2\tfracUniqFound\tfracUniqFound1\tfracUniqFound2" << std::endl;
	njh::sort(names);
	for(const auto & name : names){
		uint64_t total = 0;
		uint64_t totalMoreThanOnce = 0;
		uint64_t totalMoreThanTwice = 0;
		std::unordered_set<uint64_t> moreThanOnce;
		std::unordered_set<uint64_t> moreThanTwice;
		for(const auto & countPerSet : uniqueKmersFoundPerSet[name]){
			total += countPerSet.second;
			if(countPerSet.second>1){
				totalMoreThanOnce += countPerSet.second;
				moreThanOnce.emplace(countPerSet.first);
			}
			if(countPerSet.second>2){
				totalMoreThanTwice += countPerSet.second;
				moreThanTwice.emplace(countPerSet.first);
			}
		}
		uint64_t readCount = kmersFoundPerSeq[name].size();
		long double meanPerSeq = vectorMean(kmersFoundPerSeq[name]);
		auto occMean =  static_cast<long double>(total)/uniqueKmersFoundPerSet[name].size();
		auto occMean1 = static_cast<long double>(totalMoreThanOnce)/moreThanOnce.size();
		auto occMean2 = static_cast<long double>(totalMoreThanTwice)/moreThanTwice.size();

		out << sampleName
				<< "\t" << name
				<< "\t" << readCount
				<< "\t" << meanPerSeq
				<< "\t" << total
				<< "\t" << totalMoreThanOnce
				<< "\t" << totalMoreThanTwice
				<< "\t" << uniqueKmersFoundPerSet[name].size()
				<< "\t" << moreThanOnce.size()
				<< "\t" << moreThanTwice.size()
				<< "\t" << uniqueKmersPerSet[name].size()
				<< "\t" << occMean
				<< "\t" << occMean1
				<< "\t" << occMean2
				<< "\t" << static_cast<long double>(uniqueKmersFoundPerSet[name].size())/uniqueKmersPerSet[name].size()
				<< "\t" << static_cast<long double>(moreThanOnce.size())/uniqueKmersPerSet[name].size()
				<< "\t" << static_cast<long double>(moreThanTwice.size())/uniqueKmersPerSet[name].size()
				<< std::endl;
	}

	return 0;

}




int kmerExpRunner::findUniqKmersBetweenSeqSetsMulti(const njh::progutils::CmdArgs & inputCommands){
	KmerGatherer::KmerGathererPars countPars;

	bfs::path seqSetTableFnp = "";
	bfs::path seqSetSuppFastaTableFnp = "";
	bool fasta = false;
	bool noFilters = false;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.description_ = "Get kmers that appear within the sequences supplied and are unique to that set compared to the other sets";
	setUp.setOption(fasta, "--fasta", "file contains fasta files instead of 2bit files");
	setUp.setOption(noFilters, "--noFilters", "Don't do filtering on kmers sets");
	if(noFilters){
		countPars.entropyFilter_ = 0;
		countPars.allowableCharacters_ = njh::genSetOfAnsiPrintable();
	}

	std::string columnsHelp = njh::pasteAsStr("1)set,2)", (fasta? "fasta": "2bit"));
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

	std::function<bool(const std::string&)> seqCheck = [&countPars](const std::string & k){
		return std::all_of(k.begin(), k.end(), [&countPars](char base){return njh::in(base, countPars.allowableCharacters_);});
	};


	{
		setUp.rLog_.logCurrentTime("count_all");
		setUp.rLog_.runLogFile_.flush();
		std::unordered_map<std::string, std::set<uint64_t>> allKmers;
		if (!fasta) {
			allKmers = kGather.getUniqueKmersSetHashWithFilters(seqFiles);
		} else {
			allKmers = kGather.getUniqueKmersSetHashWithFiltersFromFastas(seqFiles);
		}
		if(!seqSetSuppFastaTableFnp.empty()){
			std::vector<bfs::path> fastaFiles;
			for(const auto & row : inputFastaSupp){
				seqFilesForSet[row[inputFastaSupp.getColPos("set")]].emplace(row[inputFastaSupp.getColPos("fasta")]);
				fastaFiles.emplace_back(row[inputFastaSupp.getColPos("fasta")]);
			}
			auto suppKmers = kGather.getUniqueKmersSetHashWithFiltersFromFastas(fastaFiles);
			for(const auto & supp : suppKmers){
				allKmers.emplace(supp);
			}
		}
		setUp.rLog_.logCurrentTime("condense");
		setUp.rLog_.runLogFile_.flush();
		njh::concurrent::LockableQueue<std::string> seqSetNamesQueue(getVectorOfMapKeys(seqFilesForSet));
		for(const auto & name : seqFilesForSet){
			kmersPerSet[name.first] = std::set<uint64_t>{};
		}
		std::function<void()> condenseKmers = [&seqSetNamesQueue,&allKmers,&seqFilesForSet,&kmersPerSet](){
			std::string name;
			while(seqSetNamesQueue.getVal(name)){
				for(const auto & twobit : seqFilesForSet.at(name)){
//					std::cout << __FILE__ << " " << __LINE__ << std::endl;
//					std::cout << "twobit: " << twobit << std::endl;
//					std::cout << allKmers.at(twobit).size() << std::endl;
					kmersPerSet[name].insert(allKmers.at(twobit).begin(), allKmers.at(twobit).end());
//					for(const auto & k : allKmers.at(twobit)){
//						kmersPerSet[name].emplace(k);
//						if(seqCheck(hasher.reverseHash(k))){
//							kmersPerSet[name].emplace(k);
//						}
//					}
				}
			}
		};
		njh::concurrent::runVoidFunctionThreaded(condenseKmers, countPars.numThreads_);
	}
	std::map<std::string, std::set<uint64_t>> uniqueKmersFinal;
	setUp.rLog_.logCurrentTime("compare");
	setUp.rLog_.runLogFile_.flush();
	for(const auto & kmersForSet : kmersPerSet){
		uniqueKmersFinal[kmersForSet.first] = std::set<uint64_t>{};
	}
	{

		auto namesFound = getVectorOfMapKeys(kmersPerSet);
		if(namesFound.size() > 1){
			njh::concurrent::LockableQueue<std::string> seqSetNamesQueue(namesFound);
			std::function<void()> compareKmers = [&seqSetNamesQueue,&kmersPerSet,&uniqueKmersFinal](){
				std::string name;
				while(seqSetNamesQueue.getVal(name)){
					std::set<uint64_t> uniqueKmers;
					uint32_t count = 0;
					for(const auto & otherSet : kmersPerSet){
						if(otherSet.first == name){
							continue;
						}
						if(0 == count){
							std::vector<uint64_t> notShared;
							std::set_difference(
									kmersPerSet.at(name).begin(), kmersPerSet.at(name).end(),
									otherSet.second.begin(), otherSet.second.end(),
									std::back_inserter(notShared));
							uniqueKmers = njh::vecToSet(notShared);
						}else{
							std::vector<uint64_t> notShared;
							std::set_difference(
									uniqueKmers.begin(), uniqueKmers.end(),
									otherSet.second.begin(), otherSet.second.end(),
									std::back_inserter(notShared));
							uniqueKmers = njh::vecToSet(notShared);
						}
						++count;
					}
					uniqueKmersFinal[name] = uniqueKmers;
				}
			};
			njh::concurrent::runVoidFunctionThreaded(compareKmers, countPars.numThreads_);
		}else if( namesFound.size() == 1){
			uniqueKmersFinal[namesFound.front()] = kmersPerSet[namesFound.front()];
		}
	}

	SimpleKmerHash hasher;
	OutputStream out(njh::files::make_path(setUp.pars_.directoryName_, "uniqueKmers.tab.txt.gz"));
	for(const auto & kmersForSet : uniqueKmersFinal){
		for(const auto & kmer : kmersForSet.second){
			out << kmersForSet.first
					<< "\t" << hasher.reverseHash(kmer) << "\n";
		}
	}
	return 0;
}

int kmerExpRunner::findKmersInSets(const njh::progutils::CmdArgs & inputCommands){
	KmerGatherer::KmerGathererPars countPars;
	std::vector<bfs::path> seqSetFnps;
	uint32_t minOccurrences = 1;
	seqSetUp setUp(inputCommands);
	setUp.description_ = "Take a sequence file and report for each kmer stepped along the seq if that kmer can be found within a set of other fasta files";
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(seqSetFnps, "--seqSetFnps", "seqSetFnps", true);
	setUp.setOption(minOccurrences, "--minOccurrences", "minOccurrences");

	countPars.setOptions(setUp);
	setUp.processReadInNames(true);

	setUp.processDirectoryOutputName(true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);

	KmerGatherer kGather(countPars);
	std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>> allCounts;
	setUp.rLog_.setCurrentLapName("initial");
	for(const auto & seqSetFnp : seqSetFnps){
		std::string genomeFnp1Base = bfs::basename(njh::files::removeExtension(seqSetFnp));
		setUp.rLog_.logCurrentTime("genome1_counting-" + genomeFnp1Base);
		allCounts[genomeFnp1Base] = kGather.countGenomeKmers(seqSetFnp);
	}

	seqInfo seq;
	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn()	;
	OutputStream kmerClassifierOut(njh::files::make_path(setUp.pars_.directoryName_, "seqsKmerClassified.tab.txt.gz"));
	kmerClassifierOut << "name\tpos\tkmer\tfoundIn" << std::endl;
	while(reader.readNextRead(seq)){
		if(len(seq) >= countPars.kmerLength_){
			for(auto pos = 0; pos < (len(seq) + 1 - countPars.kmerLength_); ++pos){
				auto currentK = seq.seq_.substr(pos, countPars.kmerLength_);
				std::set<std::string> setsFounds;
				for( auto & seqSet : allCounts){
					if(njh::in(currentK, seqSet.second) && seqSet.second[currentK] >= minOccurrences){
						setsFounds.emplace(seqSet.first);
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

	setUp.rLog_.logCurrentTime("end");

	return 0;
}





}  //namespace njhseq
