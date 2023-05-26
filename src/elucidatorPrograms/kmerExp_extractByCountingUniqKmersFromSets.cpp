//
// Created by Nicholas Hathaway on 5/26/23.
//

#include "kmerExp.hpp"
#include <njhseq/IO/SeqIO/MultiSeqIO.hpp>
#include "elucidator/helpers/UniqueKmerSetHelper.hpp"

namespace njhseq {

int kmerExpRunner::countingUniqKmersFromSetsPerRead(const njh::progutils::CmdArgs & inputCommands){
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
	setUp.setOption(includeRevComp, "--includeRevComp", "include Rev Comp of the input seqs");




	setUp.processWritingOptions(outOpts);
	//setUp.processDirectoryOutputName("true");
	setUp.finishSetUp(std::cout);
	//setUp.startARunLog(setUp.pars_.directoryName_);
	njh::stopWatch watch;
	watch.setLapName("initial");
	std::unordered_map<std::string, std::unordered_set<uint64_t>> uniqueKmersPerSet;

	uint32_t klen = 0;
	watch.startNewLap("reading in unique kmer table");
	{
		SimpleKmerHash hasher;
		TableReader uniqKmers(TableIOOpts::genTabFileIn(countTable, false));
		if(uniqKmers.header_.nCol() < 2){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "need to have 2 columns" << "\n";
			throw std::runtime_error{ss.str()};
		}
		VecStr row;
		while(uniqKmers.getNextRow(row)){
			klen = row[1].size();
			uniqueKmersPerSet[row[0]].emplace(hasher.hash(row[1]));
		}
	}
	if(setUp.pars_.verbose_){
		std::cout << watch.getLapName() << "\t" << watch.timeLapFormatted() <<std::endl;
	}
	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();

	OutputStream out(outOpts);
	std::mutex outMut;
	std::function<void()> readInComp;
	VecStr names = getVectorOfMapKeys(uniqueKmersPerSet);
	out << "sample\tname\ttotalUniqueKmers";
	out << "\tset";
	out << "\ttotalUniqueInSet";
	out << "\tkmersFoundFromSet";
	out << "\tfracInSet";
	out << "\tfracOfSetFound";
	if(includeRevComp){
		out << "\tkmersFoundFromSetRevComp";
		out << "\tfracInSetRevComp";
		out << "\tfracOfSetFoundRevComp";
	}
	out << "\twinnerSet\twinnerSetFrac\twinnerSetRecComp";
	out << std::endl;

	if (setUp.pars_.ioOptions_.isPairedIn()) {
		readInComp = [&reader, &uniqueKmersPerSet,&klen,&includeRevComp,&out,&outMut,&sampleName]() {
			//readInComp = [&reader, &uniqueKmersPerSet, &uniqueKmersFoundPerSet,&kmersFoundPerSeq,&mut,&klen,&includeRevComp]() {

			SimpleKmerHash hasher;
			PairedRead pseq;

			while(reader.readNextReadLock(pseq)){


				std::unordered_map<uint64_t, uint64_t> hashedInputKmers;
				std::unordered_map<uint64_t, uint64_t> hashedInputKmersRevComp;
				if(len(pseq.seqBase_.seq_) > klen){
					for(uint32_t pos = 0; pos < len(pseq.seqBase_.seq_) - klen + 1; ++pos){
						auto hash = hasher.hash(pseq.seqBase_.seq_.substr(pos, klen));
						++hashedInputKmers[hash];
					}
				}
				if(len(pseq.mateSeqBase_.seq_) > klen){
					for(uint32_t pos = 0; pos < len(pseq.mateSeqBase_.seq_) - klen + 1; ++pos){
						auto hash = hasher.hash(pseq.mateSeqBase_.seq_.substr(pos, klen));
						++hashedInputKmers[hash];
					}
				}
				if(includeRevComp){
					if(len(pseq.seqBase_.seq_) > klen){
						for(uint32_t pos = 0; pos < len(pseq.seqBase_.seq_) - klen + 1; ++pos){
							auto hash = hasher.revCompHash(pseq.seqBase_.seq_.substr(pos, klen));
							++hashedInputKmersRevComp[hash];
						}
					}
					if(len(pseq.mateSeqBase_.seq_) > klen){
						for(uint32_t pos = 0; pos < len(pseq.mateSeqBase_.seq_) - klen + 1; ++pos){
							auto hash = hasher.revCompHash(pseq.mateSeqBase_.seq_.substr(pos, klen));
							++hashedInputKmersRevComp[hash];
						}
					}
				}

				std::unordered_map<std::string, uint32_t> foundPerSet;
				std::unordered_map<std::string, uint32_t> foundPerSetRevComp;
				for(const auto & setName  : njh::getVecOfMapKeys(uniqueKmersPerSet)){
					foundPerSet[setName] = 0;
					foundPerSetRevComp[setName] = 0;
				}
				for(const auto & hashedKmer : hashedInputKmers){
					for(const auto & uniqueKmers : uniqueKmersPerSet){
						if(njh::in(hashedKmer.first, uniqueKmers.second)){
							++foundPerSet[uniqueKmers.first];
						}
					}
				}

				if(includeRevComp){
					for(const auto & hashedKmer : hashedInputKmersRevComp){
						for(const auto & uniqueKmers : uniqueKmersPerSet){
							if(njh::in(hashedKmer.first, uniqueKmers.second)){
								++foundPerSetRevComp[uniqueKmers.first];
							}
						}
					}
				}
				{
					std::lock_guard<std::mutex> lockGuard(outMut);

					std::string winnerSet = "undetermined";
					double bestFrac = 0;
					bool winnerRevComp = false;

					for(const auto & setName  : njh::getVecOfMapKeys(uniqueKmersPerSet)){
						if(static_cast<double>(foundPerSet[setName])/static_cast<double>(hashedInputKmers.size()) > bestFrac){
							bestFrac = static_cast<double>(foundPerSet[setName])/static_cast<double>(hashedInputKmers.size());
							winnerSet = setName;
							winnerRevComp = false;
						}
						if (includeRevComp) {
							if(static_cast<double>(foundPerSetRevComp[setName])/static_cast<double>(hashedInputKmers.size()) > bestFrac){
								bestFrac = static_cast<double>(foundPerSetRevComp[setName])/static_cast<double>(hashedInputKmers.size());
								winnerSet = setName;
								winnerRevComp = true;
							}
						}
					}

					for(const auto & setName  : njh::getVecOfMapKeys(uniqueKmersPerSet)){
						out << sampleName
								<< "\t" << pseq.seqBase_.name_
								<< "\t" << hashedInputKmers.size()
								<< "\t" << setName
								<< "\t" << uniqueKmersPerSet[setName].size()
								<< "\t" << foundPerSet[setName]
								<< "\t" << static_cast<double>(foundPerSet[setName])/static_cast<double>(hashedInputKmers.size())
								<< "\t" << static_cast<double>(foundPerSet[setName])/static_cast<double>(uniqueKmersPerSet[setName].size());
						if (includeRevComp) {
							out << "\t" << foundPerSetRevComp[setName]
									<< "\t" << static_cast<double>(foundPerSetRevComp[setName])/static_cast<double>(hashedInputKmers.size())
									<< "\t" << static_cast<double>(foundPerSetRevComp[setName])/static_cast<double>(uniqueKmersPerSet[setName].size());
						}
						out << "\t" << winnerSet
								<< "\t" << bestFrac
								<< "\t" << njh::boolToStr(winnerRevComp);
						out << std::endl;
					}
				}
			}
		};
	} else {
		readInComp = [&reader, &uniqueKmersPerSet,&klen,&includeRevComp,&out,&outMut,&sampleName]() {
			SimpleKmerHash hasher;
			seqInfo seq;

			while(reader.readNextReadLock(seq)){

				std::unordered_map<uint64_t, uint64_t> hashedInputKmers;
				std::unordered_map<uint64_t, uint64_t> hashedInputKmersRevComp;
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
							++hashedInputKmersRevComp[hash];
						}
					}
				}

				std::unordered_map<std::string, uint32_t> foundPerSet;
				std::unordered_map<std::string, uint32_t> foundPerSetRevComp;
				for(const auto & setName  : njh::getVecOfMapKeys(uniqueKmersPerSet)){
					foundPerSet[setName] = 0;
					foundPerSetRevComp[setName] = 0;
				}
				for(const auto & hashedKmer : hashedInputKmers){
					for(const auto & uniqueKmers : uniqueKmersPerSet){
						if(njh::in(hashedKmer.first, uniqueKmers.second)){
							++foundPerSet[uniqueKmers.first];
						}
					}
				}
				if(includeRevComp){
					for(const auto & hashedKmer : hashedInputKmersRevComp){
						for(const auto & uniqueKmers : uniqueKmersPerSet){
							if(njh::in(hashedKmer.first, uniqueKmers.second)){
								++foundPerSetRevComp[uniqueKmers.first];
							}
						}
					}
				}
				{
					std::lock_guard<std::mutex> lockGuard(outMut);

					std::string winnerSet = "undetermined";
					double bestFrac = 0;
					bool winnerRevComp = false;

					for(const auto & setName  : njh::getVecOfMapKeys(uniqueKmersPerSet)){
						if(static_cast<double>(foundPerSet[setName])/static_cast<double>(hashedInputKmers.size()) > bestFrac){
							bestFrac = static_cast<double>(foundPerSet[setName])/static_cast<double>(hashedInputKmers.size());
							winnerSet = setName;
							winnerRevComp = false;
						}
						if (includeRevComp) {
							if(static_cast<double>(foundPerSetRevComp[setName])/static_cast<double>(hashedInputKmers.size()) > bestFrac){
								bestFrac = static_cast<double>(foundPerSetRevComp[setName])/static_cast<double>(hashedInputKmers.size());
								winnerSet = setName;
								winnerRevComp = true;
							}
						}
					}
					for(const auto & setName  : njh::getVecOfMapKeys(uniqueKmersPerSet)){
						out << sampleName
								<< "\t" << seq.name_
								<< "\t" << hashedInputKmers.size()
								<< "\t" << setName
								<< "\t" << uniqueKmersPerSet[setName].size()
								<< "\t" << foundPerSet[setName]
								<< "\t" << static_cast<double>(foundPerSet[setName])/static_cast<double>(hashedInputKmers.size())
								<< "\t" << static_cast<double>(foundPerSet[setName])/static_cast<double>(uniqueKmersPerSet[setName].size());
						if (includeRevComp) {
							out << "\t" << foundPerSetRevComp[setName]
									<< "\t" << static_cast<double>(foundPerSetRevComp[setName])/static_cast<double>(hashedInputKmers.size())
									<< "\t" << static_cast<double>(foundPerSetRevComp[setName])/static_cast<double>(uniqueKmersPerSet[setName].size());
						}
						out << "\t" << winnerSet
								<< "\t" << bestFrac
								<< "\t" << njh::boolToStr(winnerRevComp);
						out << std::endl;
					}
				}
			}
		};
	}

	njh::concurrent::runVoidFunctionThreaded(readInComp, numThreads);
	return 0;

}


int kmerExpRunner::extractByCountingUniqKmersFromSets(const njh::progutils::CmdArgs & inputCommands){

//	VecStr testReads{};

	uint32_t numThreads = 1;
	std::set<std::string> excludeSetNames;
	bfs::path countTable;
	bfs::path excludesCountTable;
	bfs::path nonUniqueKmerTable;
	uint32_t maxIterations = std::numeric_limits<uint32_t>::max();
	std::string sampleName;
	bool includeRevComp = false;
	bool doNotWriteUndetermined = false;
	bool writeOutExclude = false;
	//bool extractAfterIterating = false;
	bool keepTemporaryFiles = false;

	bool doNotDoFinalExtract = false;
	uint32_t hardCountOff = 30;
	uint32_t addingInKmersCountCutOff = 0;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.pars_.ioOptions_.revComplMate_ = true;
	setUp.processReadInNames(true);
	setUp.setOption(keepTemporaryFiles, "--keepTemporaryFiles", "Keep Temporary Files");
	setUp.setOption(addingInKmersCountCutOff, "--addingInKmersCountCutOff", "adding In Kmers Count Cut Off");

	setUp.setOption(excludesCountTable, "--excludesCountTable", "excludes Count Table, 1)set,2)kmer");
	setUp.setOption(excludeSetNames, "--excludeSetNames", "names of sets of unqiue kmers to ignore from the --kmerTable");
	bool doReCheckExcludeSets = false;
	setUp.setOption(doReCheckExcludeSets, "--doReCheckExcludeSets", "do Re Check Exclude Sets on iterations, could lead to some reads being recruit away but will drastically speed up run time");

	bool doNotFinalRecheckOfExcludeSeqs = false;
	setUp.setOption(doNotFinalRecheckOfExcludeSeqs, "--doNotFinalRecheckOfExcludeSeqs", "redo a unique kmer filter on the exclusion sets during the final re-extraction process");

	bool doFinalRecheckOfExcludeSeqs = !doNotFinalRecheckOfExcludeSeqs;
//	setUp.setOption(doFinalRecheckOfExcludeSeqs, "--doFinalRecheckOfExcludeSeqs", "Normal the exclusion kmer set is left alone, if this flag is set, it will be filter during the final re-extraction process");

	setUp.setOption(doNotDoFinalExtract, "--doNotDoFinalExtract", "do Not Do Final Extract");



	setUp.setOption(nonUniqueKmerTable, "--nonUniqueKmerTable", "non-unique Kmer Table, 1)set,2)kmer");
	//setUp.setOption(extractAfterIterating, "--extractAfterIterating", "Extract After Iterating");
	setUp.setOption(maxIterations, "--maxIterations", "max Iterations to perform");


	setUp.setOption(countTable, "--kmerTable,--countTable", "countTable, 1)set,2)kmer", true);
	setUp.setOption(sampleName, "--sampleName", "Name to add to output file", true);
	setUp.setOption(doNotWriteUndetermined, "--doNotWriteUndetermined", "do Not Write Undetermined");
	setUp.setOption(writeOutExclude, "--writeOutExclude", "write Out Excluded reads that match the excluded kmer sets");
	uint32_t klen = 19;
	if(bfs::exists(countTable)){
		auto row = tokenizeString(njh::files::getFirstLine(countTable), "\t");
		klen = row[1].size();
	}
	setUp.setOption(hardCountOff, "--hardCountOff", "hard Count Off, do not count sets unless greater thant his number");
	uint32_t smallLenCutOff = hardCountOff + klen + 1;
	setUp.setOption(smallLenCutOff, "--smallLenCutOff", "small Len Cut Off of input sequences, if less than this size will skip over");

	setUp.setOption(numThreads, "--numThreads", "numThreads");
	setUp.setOption(includeRevComp, "--includeRevComp", "include Rev Comp of the input seqs");
	setUp.processDirectoryOutputName(true);

	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	njh::stopWatch watch;
	watch.setLapName("initial");
	std::map<std::string, std::unordered_set<uint64_t>> uniqueKmersPerSet;

	std::map<std::string, uint32_t> readsPerSet;
	std::map<std::string, uint32_t> readsPerSetRevComp;
	std::mutex mut;

	watch.startNewLap("reading in unique kmer table");


	{
		SimpleKmerHash hasher;

		TableReader uniqKmers(TableIOOpts::genTabFileIn(countTable, false));
		if(uniqKmers.header_.nCol() < 2){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "need to have 2 columns" << "\n";
			throw std::runtime_error{ss.str()};
		}
		VecStr row;
		while(uniqKmers.getNextRow(row)){
			uniqueKmersPerSet[row[0]].emplace(hasher.hash(row[1]));
		}
	}
	if(!excludesCountTable.empty()){
		SimpleKmerHash hasher;

		TableReader uniqKmers(TableIOOpts::genTabFileIn(excludesCountTable, false));
		if(uniqKmers.header_.nCol() < 2){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "need to have 2 columns" << "\n";
			throw std::runtime_error{ss.str()};
		}
		VecStr row;
		while(uniqKmers.getNextRow(row)){
			excludeSetNames.emplace(row[0]);
			klen = row[1].size();
			uniqueKmersPerSet[row[0]].emplace(hasher.hash(row[1]));
		}
	}

	std::unordered_set<uint64_t> nonUniqueKmersPerSet;
	std::string nonUniqueRegionName = "NON_UNIQUE";
	if(!nonUniqueKmerTable.empty()){

		SimpleKmerHash hasher;
		TableReader uniqKmers(TableIOOpts::genTabFileIn(nonUniqueKmerTable, false));
		if(uniqKmers.header_.nCol() < 2){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "need to have 2 columns" << "\n";
			throw std::runtime_error{ss.str()};
		}
		VecStr row;
		while(uniqKmers.getNextRow(row)){
			klen = row[1].size();
			//currently this will overwrite whatever the original non-unique kmer set name was
			nonUniqueKmersPerSet.emplace(hasher.hash(row[1]));
		}
	}

	if(setUp.pars_.verbose_){
		std::cout << watch.getLapName() << "\t" << watch.timeLapFormatted() <<std::endl;
	}
	MultiSeqIO seqOut;
	watch.startNewLap("initial scan");
	uint32_t initialSmallLenCutOffCount = 0;
	OutOptions outCountsOpts(njh::files::make_path(setUp.pars_.directoryName_, "extractionCounts.tab.txt"));
	OutputStream outCounts(outCountsOpts);
	{//initial

		if(setUp.pars_.ioOptions_.isPairedIn()) {
			setUp.pars_.ioOptions_.revComplMate_ = true;
			//to make below count correctly, make sure mate is automatically reversed complemented
		}

		SeqInput reader(setUp.pars_.ioOptions_);
		reader.openIn();


		std::function<void()> readInComp;
		VecStr names = getVectorOfMapKeys(uniqueKmersPerSet);
//		std::cout << "sample\tname\ttotalUniqueKmers";
//		std::cout << "\tset";
//		std::cout << "\ttotalUniqueInSet";
//		std::cout << "\tkmersFoundFromSet";
//		std::cout << "\tfracInSet";
//		std::cout << "\tfracOfSetFound";
//		if(includeRevComp){
//			std::cout << "\tkmersFoundFromSetRevComp";
//			std::cout << "\tfracInSetRevComp";
//			std::cout << "\tfracOfSetFoundRevComp";
//		}
//		std::cout << "\twinnerSet\twinnerSetFrac\twinnerSetRecComp";
//		std::cout << std::endl;

		if (setUp.pars_.ioOptions_.isPairedIn()) {
			for(const auto & name : names){
				auto seqOutOpts = SeqIOOptions::genPairedOutGz(njh::files::make_path(setUp.pars_.directoryName_, name));
				seqOut.addReader(name, seqOutOpts);
			}
			seqOut.addReader("undetermined", SeqIOOptions::genPairedOutGz(njh::files::make_path(setUp.pars_.directoryName_, "undetermined")));
			readInComp = [&reader, &uniqueKmersPerSet, &readsPerSet, &readsPerSetRevComp,
							&mut, &klen, &includeRevComp, &seqOut, &doNotWriteUndetermined,
							&hardCountOff, &excludeSetNames, &writeOutExclude,
							&initialSmallLenCutOffCount,
							&smallLenCutOff
//							, &testReads
			]() {
				//readInComp = [&reader, &uniqueKmersPerSet, &uniqueKmersFoundPerSet,&kmersFoundPerSeq,&mut,&klen,&includeRevComp]() {
				uint32_t currentInitialSmallLenCutOffCount = 0;

				SimpleKmerHash hasher;
				PairedRead pseq;
				std::unordered_map<std::string, uint32_t> readsPerSetCurrent;
				std::unordered_map<std::string, uint32_t> readsPerSetRevCompCurrent;
				while(reader.readNextReadLock(pseq)){
					std::unordered_map<uint64_t, uint64_t> hashedInputKmers;
					std::unordered_map<uint64_t, uint64_t> hashedInputKmersRevComp;

					if(len(pseq.seqBase_.seq_) < smallLenCutOff && len(pseq.mateSeqBase_.seq_) < smallLenCutOff){
						++currentInitialSmallLenCutOffCount;
						continue;
					}

					if(len(pseq.seqBase_.seq_) > klen){
						for(uint32_t pos = 0; pos < len(pseq.seqBase_.seq_) - klen + 1; ++pos){
							auto hash = hasher.hash(pseq.seqBase_.seq_.substr(pos, klen));
							++hashedInputKmers[hash];
						}
					}
					if(len(pseq.mateSeqBase_.seq_) > klen){
						for(uint32_t pos = 0; pos < len(pseq.mateSeqBase_.seq_) - klen + 1; ++pos){
							auto hash = hasher.hash(pseq.mateSeqBase_.seq_.substr(pos, klen));
							++hashedInputKmers[hash];
						}
					}
					if(includeRevComp){
						if(len(pseq.seqBase_.seq_) > klen){
							for(uint32_t pos = 0; pos < len(pseq.seqBase_.seq_) - klen + 1; ++pos){
								auto hash = hasher.revCompHash(pseq.seqBase_.seq_.substr(pos, klen));
								++hashedInputKmersRevComp[hash];
							}
						}
						if(len(pseq.mateSeqBase_.seq_) > klen){
							for(uint32_t pos = 0; pos < len(pseq.mateSeqBase_.seq_) - klen + 1; ++pos){
								auto hash = hasher.revCompHash(pseq.mateSeqBase_.seq_.substr(pos, klen));
								++hashedInputKmersRevComp[hash];
							}
						}
					}

					std::unordered_map<std::string, uint32_t> foundPerSet;
					std::unordered_map<std::string, uint32_t> foundPerSetRevComp;
					for(const auto & setName  : njh::getVecOfMapKeys(uniqueKmersPerSet)){
						foundPerSet[setName] = 0;
						foundPerSetRevComp[setName] = 0;
					}
					for(const auto & hashedKmer : hashedInputKmers){
						for(const auto & uniqueKmers : uniqueKmersPerSet){
							if(njh::in(hashedKmer.first, uniqueKmers.second)){
								++foundPerSet[uniqueKmers.first];
							}
						}
					}

					if(includeRevComp){
						for(const auto & hashedKmer : hashedInputKmersRevComp){
							for(const auto & uniqueKmers : uniqueKmersPerSet){
								if(njh::in(hashedKmer.first, uniqueKmers.second)){
									++foundPerSetRevComp[uniqueKmers.first];
								}
							}
						}
					}
					//set a hard cut off, reset counts to zero
					for(auto & perSet : foundPerSet){
						if(perSet.second < hardCountOff){
							perSet.second = 0;
						}
					}
					for(auto & perSet : foundPerSetRevComp){
						if(perSet.second < hardCountOff){
							perSet.second = 0;
						}
					}
					std::string winnerSet = "undetermined";
					double bestFrac = 0;
					bool winnerRevComp = false;

					for(const auto & setName  : njh::getVecOfMapKeys(uniqueKmersPerSet)){
						if(static_cast<double>(foundPerSet[setName])/static_cast<double>(hashedInputKmers.size()) > bestFrac){
							bestFrac = static_cast<double>(foundPerSet[setName])/static_cast<double>(hashedInputKmers.size());
							winnerSet = setName;
							winnerRevComp = false;
						}
						if (includeRevComp) {
							if(static_cast<double>(foundPerSetRevComp[setName])/static_cast<double>(hashedInputKmers.size()) > bestFrac){
								bestFrac = static_cast<double>(foundPerSetRevComp[setName])/static_cast<double>(hashedInputKmers.size());
								winnerSet = setName;
								winnerRevComp = true;
							}
						}
					}
					if (winnerRevComp) {
						++readsPerSetRevCompCurrent[winnerSet];
						pseq.seqBase_.reverseComplementRead(false, true);
						pseq.mateSeqBase_.reverseComplementRead(false, true);
					} else {
						++readsPerSetCurrent[winnerSet];
					}
//
//					if(njh::in(pseq.seqBase_.name_, testReads)){
//						for(const auto & setName  : njh::getVecOfMapKeys(uniqueKmersPerSet)) {
//							std::cout << "test"
//												<< "\t" << pseq.seqBase_.name_
//												<< "\t" << hashedInputKmers.size()
//												<< "\t" << setName
//												<< "\t" << uniqueKmersPerSet[setName].size()
//												<< "\t" << foundPerSet[setName]
//												<< "\t"
//												<< static_cast<double>(foundPerSet[setName]) / static_cast<double>(hashedInputKmers.size())
//												<< "\t" << static_cast<double>(foundPerSet[setName]) /
//																	 static_cast<double>(uniqueKmersPerSet[setName].size());
//							if (includeRevComp) {
//								std::cout << "\t" << foundPerSetRevComp[setName]
//													<< "\t" << static_cast<double>(foundPerSetRevComp[setName]) /
//																		 static_cast<double>(hashedInputKmers.size())
//													<< "\t" << static_cast<double>(foundPerSetRevComp[setName]) /
//																		 static_cast<double>(uniqueKmersPerSet[setName].size());
//							}
//							std::cout << "\t" << winnerSet
//												<< "\t" << bestFrac
//												<< "\t" << njh::boolToStr(winnerRevComp);
//							std::cout << std::endl;
//						}
//					}

					if(!doNotWriteUndetermined || winnerSet != "undetermined"){
						if(writeOutExclude || !njh::in(winnerSet, excludeSetNames)){
							seqOut.openWrite(winnerSet, pseq);
						}
					}
				}
				{
					std::lock_guard<std::mutex> lockGuard(mut);
					initialSmallLenCutOffCount += currentInitialSmallLenCutOffCount;
					for(const auto & readsPerSetCount : readsPerSetCurrent){
						readsPerSet[readsPerSetCount.first] += readsPerSetCount.second;
					}
					for(const auto & readsPerSetRevCompCount : readsPerSetRevCompCurrent){
						readsPerSetRevComp[readsPerSetRevCompCount.first] += readsPerSetRevCompCount.second;
					}
				}
			};
		} else {
			for(const auto & name : names){
				auto seqOutOpts = SeqIOOptions::genFastqOutGz(njh::files::make_path(setUp.pars_.directoryName_, name) );
				seqOut.addReader(name, seqOutOpts);
			}
			seqOut.addReader("undetermined", SeqIOOptions::genFastqOutGz(njh::files::make_path(setUp.pars_.directoryName_, "undetermined")));
			readInComp = [&reader, &uniqueKmersPerSet, &readsPerSet,&readsPerSetRevComp,
							&mut,&klen,&includeRevComp,&seqOut,&doNotWriteUndetermined,
							&hardCountOff, &excludeSetNames, &writeOutExclude,
							&smallLenCutOff,
							&initialSmallLenCutOffCount]() {
//		readInComp = [&reader, &uniqueKmersPerSet, &uniqueKmersFoundPerSet,&kmersFoundPerSeq,&mut,&klen,&includeRevComp]() {
				SimpleKmerHash hasher;
				seqInfo seq;
				std::unordered_map<std::string, uint32_t> readsPerSetCurrent;
				std::unordered_map<std::string, uint32_t> readsPerSetRevCompCurrent;
				uint32_t currentInitialSmallLenCutOffCount = 0;
				while(reader.readNextReadLock(seq)){
					if(len(seq) < smallLenCutOff){
						++currentInitialSmallLenCutOffCount;
					}
					std::unordered_map<uint64_t, uint64_t> hashedInputKmers;
					std::unordered_map<uint64_t, uint64_t> hashedInputKmersRevComp;
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
								++hashedInputKmersRevComp[hash];
							}
						}
					}

					std::unordered_map<std::string, uint32_t> foundPerSet;
					std::unordered_map<std::string, uint32_t> foundPerSetRevComp;
					for(const auto & setName  : njh::getVecOfMapKeys(uniqueKmersPerSet)){
						foundPerSet[setName] = 0;
						foundPerSetRevComp[setName] = 0;
					}
					for(const auto & hashedKmer : hashedInputKmers){
						for(const auto & uniqueKmers : uniqueKmersPerSet){
							if(njh::in(hashedKmer.first, uniqueKmers.second)){
								++foundPerSet[uniqueKmers.first];
							}
						}
					}
					if(includeRevComp){
						for(const auto & hashedKmer : hashedInputKmersRevComp){
							for(const auto & uniqueKmers : uniqueKmersPerSet){
								if(njh::in(hashedKmer.first, uniqueKmers.second)){
									++foundPerSetRevComp[uniqueKmers.first];
								}
							}
						}
					}
					//set a hard cut off, reset counts to zero
					for(auto & perSet : foundPerSet){
						if(perSet.second < hardCountOff){
							perSet.second = 0;
						}
					}
					for(auto & perSet : foundPerSetRevComp){
						if(perSet.second < hardCountOff){
							perSet.second = 0;
						}
					}
					std::string winnerSet = "undetermined";
					double bestFrac = 0;
					bool winnerRevComp = false;

					for(const auto & setName  : njh::getVecOfMapKeys(uniqueKmersPerSet)){
						if(static_cast<double>(foundPerSet[setName])/static_cast<double>(hashedInputKmers.size()) > bestFrac){
							bestFrac = static_cast<double>(foundPerSet[setName])/static_cast<double>(hashedInputKmers.size());
							winnerSet = setName;
							winnerRevComp = false;
						}
						if (includeRevComp) {
							if(static_cast<double>(foundPerSetRevComp[setName])/static_cast<double>(hashedInputKmers.size()) > bestFrac){
								bestFrac = static_cast<double>(foundPerSetRevComp[setName])/static_cast<double>(hashedInputKmers.size());
								winnerSet = setName;
								winnerRevComp = true;
							}
						}
					}

					if(winnerRevComp){
						++readsPerSetRevCompCurrent[winnerSet];
						seq.reverseComplementRead(true, true);
					}else{
						++readsPerSetCurrent[winnerSet];
					}
					if(!doNotWriteUndetermined || winnerSet != "undetermined"){
						if(writeOutExclude || !njh::in(winnerSet, excludeSetNames)){
							seqOut.openWrite(winnerSet, seq);
						}
					}
				}
				{
					std::lock_guard<std::mutex> lockGuard(mut);

					initialSmallLenCutOffCount += currentInitialSmallLenCutOffCount;
					for(const auto & readsPerSetCount : readsPerSetCurrent){
						readsPerSet[readsPerSetCount.first] += readsPerSetCount.second;
					}
					for(const auto & readsPerSetRevCompCount : readsPerSetRevCompCurrent){
						readsPerSetRevComp[readsPerSetRevCompCount.first] += readsPerSetRevCompCount.second;
					}
				}
			};
		}
		njh::concurrent::runVoidFunctionThreaded(readInComp, numThreads);


		uint64_t totalReads = initialSmallLenCutOffCount;
		for(const auto & readsPerSetCount : readsPerSet){
			totalReads += readsPerSetCount.second;
		}
		for(const auto & readsPerSetRevCompCount : readsPerSetRevComp){
			totalReads += readsPerSetRevCompCount.second;
		}

		if (includeRevComp) {
			outCounts << "iteration\tsample\ttotalReads\ttarget\tcount\tfrac\tforwardCount\tfracForward";
			outCounts << std::endl;
			for(const auto & setName : uniqueKmersPerSet){
				uint64_t totalExtracted = readsPerSet[setName.first] + readsPerSetRevComp[setName.first];
				outCounts << "initial"
									<< "\t" << sampleName
									<< "\t" << totalReads
									<< "\t" << setName.first
									<< "\t" << totalExtracted
									<< "\t" << static_cast<double>(totalExtracted) / static_cast<double>(totalReads)
									<< "\t" << readsPerSet[setName.first]
									<< "\t" << (totalExtracted > 0 ? static_cast<double>(readsPerSet[setName.first]) / static_cast<double>(totalExtracted) : 0)
									<< std::endl;
			}
			{
				uint64_t totalExtracted = readsPerSet["undetermined"] + readsPerSetRevComp["undetermined"];
				outCounts << "initial"
									<< "\t" << sampleName
									<< "\t" << totalReads
									<< "\t" << "undetermined"
									<< "\t" << totalExtracted
									<< "\t" << static_cast<double>(totalExtracted) / static_cast<double>(totalReads)
									<< "\t" << readsPerSet["undetermined"]
									<< "\t" << (totalExtracted > 0 ?static_cast<double>(readsPerSet["undetermined"]) / static_cast<double>(totalExtracted) : 0)
									<< std::endl;
			}
			{
				outCounts << "initial"
									<< "\t" << sampleName
									<< "\t" << totalReads
									<< "\t" << "smallLenReads"
									<< "\t" << initialSmallLenCutOffCount
									<< "\t" << static_cast<double>(initialSmallLenCutOffCount) / static_cast<double>(totalReads)
									<< "\t" << initialSmallLenCutOffCount
									<< "\t" << 0
									<< std::endl;
			}
		} else {
			outCounts << "iteration\tsample\ttotalReads\ttarget\tcount\tfrac";
			outCounts << std::endl;
			for(const auto & setName : uniqueKmersPerSet){
				outCounts << "initial"
									<< "\t" << sampleName
									<< "\t" << totalReads
									<< "\t" << setName.first
									<< "\t" << readsPerSet[setName.first]
									<< "\t" << static_cast<double>(readsPerSet[setName.first]) / static_cast<double>(totalReads)
									<< std::endl;
			}
			{
				outCounts << "initial"
									<< "\t" << sampleName
									<< "\t" << totalReads
									<< "\t" << "undetermined"
									<< "\t" << readsPerSet["undetermined"]
									<< "\t" << static_cast<double>(readsPerSet["undetermined"]) / static_cast<double>(totalReads)
									<< std::endl;
			}
			{
				outCounts << "initial"
									<< "\t" << sampleName
									<< "\t" << totalReads
									<< "\t" << "smallLenReads"
									<< "\t" << initialSmallLenCutOffCount
									<< "\t" << static_cast<double>(initialSmallLenCutOffCount) / static_cast<double>(totalReads)
									<< std::endl;
			}
		}
	}


	//read in the extract reads and add to the unique kmer sets
	bool iterate = true;
	uint64_t totalUndeterminedExtracted = readsPerSet["undetermined"] + readsPerSetRevComp["undetermined"];
	uint32_t iterNumber = 0;
	if(0 == totalUndeterminedExtracted){
		iterate = false;
	}
	auto initial_readsPerSet = readsPerSet;
	auto initial_readsPerSetRevComp = readsPerSetRevComp;
	readsPerSet["undetermined"] = 0;
	readsPerSetRevComp["undetermined"] = 0;
	uint64_t totalDeterminedReads = 0;
	for(const auto & readsPerSetCount : readsPerSet){
		totalDeterminedReads += readsPerSetCount.second;
	}
	for(const auto & readsPerSetRevCompCount : readsPerSetRevComp){
		totalDeterminedReads += readsPerSetRevCompCount.second;
	}


	while(iterate && iterNumber < maxIterations){

		watch.startNewLap(njh::pasteAsStr(iterNumber, "- read in new kmers"));

		seqOut.closeOutForReopeningAll();
		std::map<std::string, std::unordered_map<uint64_t, uint32_t>> rawKmersPerInput;
		for (const auto &kmerSet: uniqueKmersPerSet) {
			auto inOpts = SeqIOOptions::genPairedInGz(
							njh::files::make_path(setUp.pars_.directoryName_, kmerSet.first + "_R1.fastq.gz"),
							njh::files::make_path(setUp.pars_.directoryName_, kmerSet.first + "_R2.fastq.gz"));
			inOpts.revComplMate_ = true;
			SimpleKmerHash hasher;
			if (inOpts.inExists()) {
				PairedRead pseq;
				SeqInput extractedReader(inOpts);
				extractedReader.openIn();
				while (extractedReader.readNextRead(pseq)) {
					for (uint32_t pos = 0; pos < len(pseq.seqBase_.seq_) - klen + 1; ++pos) {
						++rawKmersPerInput[kmerSet.first][hasher.hash(pseq.seqBase_.seq_.substr(pos, klen))];
					}
					for (uint32_t pos = 0; pos < len(pseq.mateSeqBase_.seq_) - klen + 1; ++pos) {
						++rawKmersPerInput[kmerSet.first][hasher.hash(pseq.mateSeqBase_.seq_.substr(pos, klen))];
					}
				}
			}
		}

		std::map<std::string, std::unordered_set<uint64_t>> outputUniqueKmersPerSet;
		watch.startNewLap(njh::pasteAsStr(iterNumber, "- check new kmers against current set"));
		for (const auto & finalNewKmerSet : rawKmersPerInput){
			for (const auto finalNewKmer: finalNewKmerSet.second) {
				if(finalNewKmer.second <addingInKmersCountCutOff){
					continue; //skip if low count, helps to avoid just PCR and sequencing noisy kmers
				}
				bool pass = true;
				if (njh::in(finalNewKmer.first, nonUniqueKmersPerSet)) {
					pass = false;
				} else {
					for (const auto &set: uniqueKmersPerSet) {
						if (set.first == finalNewKmerSet.first) {
							continue;
						}
						if (njh::in(finalNewKmer.first, set.second)) {
							pass = false;
							nonUniqueKmersPerSet.emplace(finalNewKmer.first);
							break;
						}
					}
				}
				if (pass) {
					outputUniqueKmersPerSet[finalNewKmerSet.first].emplace(finalNewKmer.first);
				}
			}
		}
		watch.startNewLap(njh::pasteAsStr(iterNumber, "- check new kmers against current set - other sets"));
		//filter other sets
		for (const auto &set: uniqueKmersPerSet) {
			if(!doReCheckExcludeSets && njh::in(set.first, excludeSetNames)){
				continue;
			}
			for(const auto & finalKmer : set.second){
				bool pass = true;
				for(const auto & rawInput : rawKmersPerInput){
					if (set.first == rawInput.first) {
						pass = false;
						continue;
					}
					if(njh::in(finalKmer, rawInput.second)){
						pass = false;
						break;
					}
				}
				if(pass){
					outputUniqueKmersPerSet[set.first].emplace(finalKmer);
				}
			}
		}
		if(!doReCheckExcludeSets){
			for(const auto & excludeName : excludeSetNames){
				outputUniqueKmersPerSet[excludeName] = std::move(uniqueKmersPerSet[excludeName]);
			}
		}

		watch.startNewLap(njh::pasteAsStr(iterNumber, "- copy new set"));
		uniqueKmersPerSet = std::move(outputUniqueKmersPerSet);
		watch.startNewLap(njh::pasteAsStr(iterNumber, "- scan again"));
		bfs::path orig_undeterminedInput = njh::files::make_path(setUp.pars_.directoryName_, "undetermined.fastq.gz");
		bfs::path orig_undeterminedInputR1 = njh::files::make_path(setUp.pars_.directoryName_, "undetermined_R1.fastq.gz");
		bfs::path orig_undeterminedInputR2 = njh::files::make_path(setUp.pars_.directoryName_, "undetermined_R2.fastq.gz");

		bfs::path undeterminedInput = njh::files::make_path(setUp.pars_.directoryName_, njh::pasteAsStr(iterNumber, "_undetermined.fastq.gz"));
		bfs::path undeterminedInputR1 = njh::files::make_path(setUp.pars_.directoryName_, njh::pasteAsStr(iterNumber, "_undetermined_R1.fastq.gz"));
		bfs::path undeterminedInputR2 = njh::files::make_path(setUp.pars_.directoryName_, njh::pasteAsStr(iterNumber, "_undetermined_R2.fastq.gz"));

		SeqIOOptions undeterminedInOpts;
		if (setUp.pars_.ioOptions_.isPairedIn()) {
			bfs::rename(orig_undeterminedInputR1, undeterminedInputR1);
			bfs::rename(orig_undeterminedInputR2, undeterminedInputR2);
			undeterminedInOpts.firstName_ = undeterminedInputR1;
			undeterminedInOpts.secondName_ = undeterminedInputR2;
			undeterminedInOpts.inFormat_ = SeqIOOptions::inFormats::FASTQPAIREDGZ;
			undeterminedInOpts.revComplMate_ = true;
			//to make below count correctly, make sure mate is automatically reversed complemented
		} else {
			bfs::rename(orig_undeterminedInput, undeterminedInput);
			undeterminedInOpts.firstName_ = undeterminedInput;
			undeterminedInOpts.inFormat_ = SeqIOOptions::inFormats::FASTQGZ;
		}

		SeqInput reader(undeterminedInOpts);
		reader.openIn();

		std::function<void()> readInComp;
		VecStr names = getVectorOfMapKeys(uniqueKmersPerSet);


		if (setUp.pars_.ioOptions_.isPairedIn()) {
			readInComp = [&reader, &uniqueKmersPerSet, &readsPerSet, &readsPerSetRevComp,
							&mut, &klen, &includeRevComp, &seqOut, &doNotWriteUndetermined,
							&hardCountOff, &excludeSetNames, &writeOutExclude
							//,&testReads
							//,&iterNumber
			]() {
				//readInComp = [&reader, &uniqueKmersPerSet, &uniqueKmersFoundPerSet,&kmersFoundPerSeq,&mut,&klen,&includeRevComp]() {

				SimpleKmerHash hasher;
				PairedRead pseq;
				std::unordered_map<std::string, uint32_t> readsPerSetCurrent;
				std::unordered_map<std::string, uint32_t> readsPerSetRevCompCurrent;
				while(reader.readNextReadLock(pseq)){
					std::unordered_map<uint64_t, uint64_t> hashedInputKmers;
					std::unordered_map<uint64_t, uint64_t> hashedInputKmersRevComp;
					if(len(pseq.seqBase_.seq_) > klen){
						for(uint32_t pos = 0; pos < len(pseq.seqBase_.seq_) - klen + 1; ++pos){
							auto hash = hasher.hash(pseq.seqBase_.seq_.substr(pos, klen));
							++hashedInputKmers[hash];
						}
					}
					if(len(pseq.mateSeqBase_.seq_) > klen){
						for(uint32_t pos = 0; pos < len(pseq.mateSeqBase_.seq_) - klen + 1; ++pos){
							auto hash = hasher.hash(pseq.mateSeqBase_.seq_.substr(pos, klen));
							++hashedInputKmers[hash];
						}
					}
					if(includeRevComp){
						if(len(pseq.seqBase_.seq_) > klen){
							for(uint32_t pos = 0; pos < len(pseq.seqBase_.seq_) - klen + 1; ++pos){
								auto hash = hasher.revCompHash(pseq.seqBase_.seq_.substr(pos, klen));
								++hashedInputKmersRevComp[hash];
							}
						}
						if(len(pseq.mateSeqBase_.seq_) > klen){
							for(uint32_t pos = 0; pos < len(pseq.mateSeqBase_.seq_) - klen + 1; ++pos){
								auto hash = hasher.revCompHash(pseq.mateSeqBase_.seq_.substr(pos, klen));
								++hashedInputKmersRevComp[hash];
							}
						}
					}

					std::unordered_map<std::string, uint32_t> foundPerSet;
					std::unordered_map<std::string, uint32_t> foundPerSetRevComp;
					for(const auto & setName  : njh::getVecOfMapKeys(uniqueKmersPerSet)){
						foundPerSet[setName] = 0;
						foundPerSetRevComp[setName] = 0;
					}
					for(const auto & hashedKmer : hashedInputKmers){
						for(const auto & uniqueKmers : uniqueKmersPerSet){
							if(njh::in(hashedKmer.first, uniqueKmers.second)){
								++foundPerSet[uniqueKmers.first];
							}
						}
					}

					if(includeRevComp){
						for(const auto & hashedKmer : hashedInputKmersRevComp){
							for(const auto & uniqueKmers : uniqueKmersPerSet){
								if(njh::in(hashedKmer.first, uniqueKmers.second)){
									++foundPerSetRevComp[uniqueKmers.first];
								}
							}
						}
					}
					//set a hard cut off, reset counts to zero
					for(auto & perSet : foundPerSet){
						if(perSet.second < hardCountOff){
							perSet.second = 0;
						}
					}
					for(auto & perSet : foundPerSetRevComp){
						if(perSet.second < hardCountOff){
							perSet.second = 0;
						}
					}
					std::string winnerSet = "undetermined";
					double bestFrac = 0;
					bool winnerRevComp = false;

					for(const auto & setName  : njh::getVecOfMapKeys(uniqueKmersPerSet)){
						if(static_cast<double>(foundPerSet[setName])/static_cast<double>(hashedInputKmers.size()) > bestFrac){
							bestFrac = static_cast<double>(foundPerSet[setName])/static_cast<double>(hashedInputKmers.size());
							winnerSet = setName;
							winnerRevComp = false;
						}
						if (includeRevComp) {
							if(static_cast<double>(foundPerSetRevComp[setName])/static_cast<double>(hashedInputKmers.size()) > bestFrac){
								bestFrac = static_cast<double>(foundPerSetRevComp[setName])/static_cast<double>(hashedInputKmers.size());
								winnerSet = setName;
								winnerRevComp = true;
							}
						}
					}

//					if(njh::in(pseq.seqBase_.name_, testReads)){
//						for(const auto & setName  : njh::getVecOfMapKeys(uniqueKmersPerSet)) {
//							std::cout << "iter:" << iterNumber
//												<< "\t" << pseq.seqBase_.name_
//												<< "\t" << hashedInputKmers.size()
//												<< "\t" << setName
//												<< "\t" << uniqueKmersPerSet[setName].size()
//												<< "\t" << foundPerSet[setName]
//												<< "\t"
//												<< static_cast<double>(foundPerSet[setName]) / static_cast<double>(hashedInputKmers.size())
//												<< "\t" << static_cast<double>(foundPerSet[setName]) /
//																	 static_cast<double>(uniqueKmersPerSet[setName].size());
//							if (includeRevComp) {
//								std::cout << "\t" << foundPerSetRevComp[setName]
//													<< "\t" << static_cast<double>(foundPerSetRevComp[setName]) /
//																		 static_cast<double>(hashedInputKmers.size())
//													<< "\t" << static_cast<double>(foundPerSetRevComp[setName]) /
//																		 static_cast<double>(uniqueKmersPerSet[setName].size());
//							}
//							std::cout << "\t" << winnerSet
//												<< "\t" << bestFrac
//												<< "\t" << njh::boolToStr(winnerRevComp);
//							std::cout << std::endl;
//						}
//					}

					if (winnerRevComp) {
						++readsPerSetRevCompCurrent[winnerSet];
						pseq.seqBase_.reverseComplementRead(false, true);
						pseq.mateSeqBase_.reverseComplementRead(false, true);
					} else {
						++readsPerSetCurrent[winnerSet];
					}
//					std::cout << std::this_thread::get_id() << " " << winnerSet << std::endl;
					if(!doNotWriteUndetermined || winnerSet != "undetermined"){
						if(writeOutExclude || !njh::in(winnerSet, excludeSetNames)){
							seqOut.openWrite(winnerSet, pseq);
						}
					}
				}
				{
					std::lock_guard<std::mutex> lockGuard(mut);
					for(const auto & readsPerSetCount : readsPerSetCurrent){
//						std::cout << "readsPerSetCount.first: " << readsPerSetCount.first << ", " << readsPerSetCount.second << std::endl;
						readsPerSet[readsPerSetCount.first] += readsPerSetCount.second;
					}
					for(const auto & readsPerSetRevCompCount : readsPerSetRevCompCurrent){
//						std::cout << "readsPerSetRevCompCount.first: " << readsPerSetRevCompCount.first << ", " << readsPerSetRevCompCount.second << std::endl;
						readsPerSetRevComp[readsPerSetRevCompCount.first] += readsPerSetRevCompCount.second;
					}
				}
			};
		} else {
			readInComp = [&reader, &uniqueKmersPerSet, &readsPerSet,&readsPerSetRevComp,
							&mut,&klen,&includeRevComp,&seqOut,&doNotWriteUndetermined,
							&hardCountOff, &excludeSetNames, &writeOutExclude]() {
//		readInComp = [&reader, &uniqueKmersPerSet, &uniqueKmersFoundPerSet,&kmersFoundPerSeq,&mut,&klen,&includeRevComp]() {
				SimpleKmerHash hasher;
				seqInfo seq;
				std::unordered_map<std::string, uint32_t> readsPerSetCurrent;
				std::unordered_map<std::string, uint32_t> readsPerSetRevCompCurrent;

				while(reader.readNextReadLock(seq)){

					std::unordered_map<uint64_t, uint64_t> hashedInputKmers;
					std::unordered_map<uint64_t, uint64_t> hashedInputKmersRevComp;
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
								++hashedInputKmersRevComp[hash];
							}
						}
					}

					std::unordered_map<std::string, uint32_t> foundPerSet;
					std::unordered_map<std::string, uint32_t> foundPerSetRevComp;
					for(const auto & setName  : njh::getVecOfMapKeys(uniqueKmersPerSet)){
						foundPerSet[setName] = 0;
						foundPerSetRevComp[setName] = 0;
					}
					for(const auto & hashedKmer : hashedInputKmers){
						for(const auto & uniqueKmers : uniqueKmersPerSet){
							if(njh::in(hashedKmer.first, uniqueKmers.second)){
								++foundPerSet[uniqueKmers.first];
							}
						}
					}
					if(includeRevComp){
						for(const auto & hashedKmer : hashedInputKmersRevComp){
							for(const auto & uniqueKmers : uniqueKmersPerSet){
								if(njh::in(hashedKmer.first, uniqueKmers.second)){
									++foundPerSetRevComp[uniqueKmers.first];
								}
							}
						}
					}
					//set a hard cut off, reset counts to zero
					for(auto & perSet : foundPerSet){
						if(perSet.second < hardCountOff){
							perSet.second = 0;
						}
					}
					for(auto & perSet : foundPerSetRevComp){
						if(perSet.second < hardCountOff){
							perSet.second = 0;
						}
					}
					std::string winnerSet = "undetermined";
					double bestFrac = 0;
					bool winnerRevComp = false;

					for(const auto & setName  : njh::getVecOfMapKeys(uniqueKmersPerSet)){
						if(static_cast<double>(foundPerSet[setName])/static_cast<double>(hashedInputKmers.size()) > bestFrac){
							bestFrac = static_cast<double>(foundPerSet[setName])/static_cast<double>(hashedInputKmers.size());
							winnerSet = setName;
							winnerRevComp = false;
						}
						if (includeRevComp) {
							if(static_cast<double>(foundPerSetRevComp[setName])/static_cast<double>(hashedInputKmers.size()) > bestFrac){
								bestFrac = static_cast<double>(foundPerSetRevComp[setName])/static_cast<double>(hashedInputKmers.size());
								winnerSet = setName;
								winnerRevComp = true;
							}
						}
					}
					if(winnerRevComp){
						++readsPerSetRevCompCurrent[winnerSet];
						seq.reverseComplementRead(true, true);
					}else{
						++readsPerSetCurrent[winnerSet];
					}
					if(!doNotWriteUndetermined || winnerSet != "undetermined"){
						if(writeOutExclude || !njh::in(winnerSet, excludeSetNames)){
							seqOut.openWrite(winnerSet, seq);
						}
					}
				}
				{
					std::lock_guard<std::mutex> lockGuard(mut);
					for(const auto & readsPerSetCount : readsPerSetCurrent){
						readsPerSet[readsPerSetCount.first] += readsPerSetCount.second;
					}
					for(const auto & readsPerSetRevCompCount : readsPerSetRevCompCurrent){
						readsPerSetRevComp[readsPerSetRevCompCount.first] += readsPerSetRevCompCount.second;
					}
				}
			};
		}
		njh::concurrent::runVoidFunctionThreaded(readInComp, numThreads);


		uint64_t totalReads = 0;
		for(const auto & readsPerSetCount : readsPerSet){
			totalReads += readsPerSetCount.second;
		}
		for(const auto & readsPerSetRevCompCount : readsPerSetRevComp){
			totalReads += readsPerSetRevCompCount.second;
		}

		if (includeRevComp) {

			for(const auto & setName : uniqueKmersPerSet){
				uint64_t totalExtracted = readsPerSet[setName.first] + readsPerSetRevComp[setName.first];
				outCounts << iterNumber
									<< "\t" << sampleName
									<< "\t" << totalReads
									<< "\t" << setName.first
									<< "\t" << totalExtracted
									<< "\t" << static_cast<double>(totalExtracted) / static_cast<double>(totalReads)
									<< "\t" << readsPerSet[setName.first]
									<< "\t" << (totalExtracted > 0 ? static_cast<double>(readsPerSet[setName.first]) / static_cast<double>(totalExtracted) : 0)
									<< std::endl;
			}
			{
				uint64_t totalExtracted = readsPerSet["undetermined"] + readsPerSetRevComp["undetermined"];
				outCounts << iterNumber
									<< "\t" << sampleName
									<< "\t" << totalReads
									<< "\t" << "undetermined"
									<< "\t" << totalExtracted
									<< "\t" << static_cast<double>(totalExtracted) / static_cast<double>(totalReads)
									<< "\t" << readsPerSet["undetermined"]
									<< "\t" << (totalExtracted > 0 ?static_cast<double>(readsPerSet["undetermined"]) / static_cast<double>(totalExtracted) : 0)
									<< std::endl;
			}
		} else {
			for(const auto & setName : uniqueKmersPerSet){
				outCounts << iterNumber
									<< "\t" << sampleName
									<< "\t" << totalReads
									<< "\t" << setName.first
									<< "\t" << readsPerSet[setName.first]
									<< "\t" << static_cast<double>(readsPerSet[setName.first]) / static_cast<double>(totalReads)
									<< std::endl;
			}
			{
				outCounts << iterNumber
									<< "\t" << sampleName
									<< "\t" << totalReads
									<< "\t" << "undetermined"
									<< "\t" << readsPerSet["undetermined"]
									<< "\t" << static_cast<double>(readsPerSet["undetermined"]) / static_cast<double>(totalReads)
									<< std::endl;
			}
		}
		uint64_t currentTotalUndeterminedReads = readsPerSet["undetermined"] + readsPerSetRevComp["undetermined"];
		readsPerSet["undetermined"] = 0;
		readsPerSetRevComp["undetermined"] = 0;
		uint64_t currentTotalDeterminedReads = 0;
		for(const auto & readsPerSetCount : readsPerSet){
			currentTotalDeterminedReads += readsPerSetCount.second;
		}
		for(const auto & readsPerSetRevCompCount : readsPerSetRevComp){
			currentTotalDeterminedReads += readsPerSetRevCompCount.second;
		}
		if(setUp.pars_.verbose_){
			std::cout << "iterNumber: " << iterNumber << std::endl;
			std::cout << "currentTotalUndeterminedReads: " << currentTotalUndeterminedReads << std::endl;
			std::cout << "totalDeterminedReads: " << totalDeterminedReads << std::endl;
			std::cout << "currentTotalDeterminedReads: " << currentTotalDeterminedReads << std::endl;
		}
		if(currentTotalUndeterminedReads == 0 || totalDeterminedReads == currentTotalDeterminedReads){
			//if there's no more undetermined reads or if the number determined this round was zero then stop iterating
			iterate = false;
		}
		totalDeterminedReads = currentTotalDeterminedReads;
		reader.closeIn();
		if(!keepTemporaryFiles){
			bfs::remove(undeterminedInOpts.firstName_);
			if (setUp.pars_.ioOptions_.isPairedIn()) {
				bfs::remove(undeterminedInOpts.secondName_);
			}
		}
		++iterNumber;
	}

//	if(extractAfterIterating){
	if(!doNotDoFinalExtract){
		seqOut.closeOutForReopeningAll();

		if(!doReCheckExcludeSets && doFinalRecheckOfExcludeSeqs){
			watch.startNewLap(njh::pasteAsStr(iterNumber, "- final read in new kmers"));

			std::map<std::string, std::unordered_map<uint64_t, uint32_t>> rawKmersPerInput;
			for (const auto &kmerSet: uniqueKmersPerSet) {
				auto inOpts = SeqIOOptions::genPairedInGz(
								njh::files::make_path(setUp.pars_.directoryName_, kmerSet.first + "_R1.fastq.gz"),
								njh::files::make_path(setUp.pars_.directoryName_, kmerSet.first + "_R2.fastq.gz"));
				inOpts.revComplMate_ = true;
				SimpleKmerHash hasher;
				if (inOpts.inExists()) {
					PairedRead pseq;
					SeqInput extractedReader(inOpts);
					extractedReader.openIn();
					while (extractedReader.readNextRead(pseq)) {
						for (uint32_t pos = 0; pos < len(pseq.seqBase_.seq_) - klen + 1; ++pos) {
							++rawKmersPerInput[kmerSet.first][hasher.hash(pseq.seqBase_.seq_.substr(pos, klen))];
						}
						for (uint32_t pos = 0; pos < len(pseq.mateSeqBase_.seq_) - klen + 1; ++pos) {
							++rawKmersPerInput[kmerSet.first][hasher.hash(pseq.mateSeqBase_.seq_.substr(pos, klen))];
						}
					}
				}
			}

			std::map<std::string, std::unordered_set<uint64_t>> outputUniqueKmersPerSet;
			watch.startNewLap(njh::pasteAsStr(iterNumber, "- final check new kmers against current set"));
			for (const auto & finalNewKmerSet : rawKmersPerInput){
				for (const auto finalNewKmer: finalNewKmerSet.second) {
					if(finalNewKmer.second <addingInKmersCountCutOff){
						continue; //skip if low count, helps to avoid just PCR and sequencing noisy kmers
					}
					bool pass = true;
					if (njh::in(finalNewKmer.first, nonUniqueKmersPerSet)) {
						pass = false;
					} else {
						for (const auto &set: uniqueKmersPerSet) {
							if (set.first == finalNewKmerSet.first) {
								continue;
							}
							if (njh::in(finalNewKmer.first, set.second)) {
								pass = false;
								nonUniqueKmersPerSet.emplace(finalNewKmer.first);
								break;
							}
						}
					}
					if (pass) {
						outputUniqueKmersPerSet[finalNewKmerSet.first].emplace(finalNewKmer.first);
					}
				}
			}
			watch.startNewLap(njh::pasteAsStr(iterNumber, "- final check new kmers against current set - other sets"));
			//filter other sets
			for (const auto &set: uniqueKmersPerSet) {

				for(const auto & finalKmer : set.second){
					bool pass = true;
					for(const auto & rawInput : rawKmersPerInput){
						if (set.first == rawInput.first) {
							pass = false;
							continue;
						}
						if(njh::in(finalKmer, rawInput.second)){
							pass = false;
							break;
						}
					}
					if(pass){
						outputUniqueKmersPerSet[set.first].emplace(finalKmer);
					}
				}
			}
			watch.startNewLap(njh::pasteAsStr(iterNumber, "- final copy new set"));
			uniqueKmersPerSet = std::move(outputUniqueKmersPerSet);
		}
		uint32_t finalSmallLenCutOffCount = 0;

		watch.startNewLap(njh::pasteAsStr(iterNumber, "- scan whole set again"));
		MultiSeqIO finalSeqOut;
		std::map<std::string, uint32_t> final_readsPerSet;
		std::map<std::string, uint32_t> final_readsPerSetRevComp;

		//initial
		auto finalExtractionDir = njh::files::make_path(setUp.pars_.directoryName_, "finalExtraction");
		njh::files::makeDir(njh::files::MkdirPar{finalExtractionDir});
		if(setUp.pars_.ioOptions_.isPairedIn()) {
			setUp.pars_.ioOptions_.revComplMate_ = true;
			//to make below count correctly, make sure mate is automatically reversed complemented
		}

		SeqInput reader(setUp.pars_.ioOptions_);
		reader.openIn();

		std::function<void()> readInComp;
		VecStr names = getVectorOfMapKeys(uniqueKmersPerSet);


		if (setUp.pars_.ioOptions_.isPairedIn()) {
			for(const auto & name : names){
				auto seqOutOpts = SeqIOOptions::genPairedOutGz(njh::files::make_path(finalExtractionDir, name));
				finalSeqOut.addReader(name, seqOutOpts);
			}
			finalSeqOut.addReader("undetermined", SeqIOOptions::genPairedOutGz(njh::files::make_path(finalExtractionDir, "undetermined")));
			readInComp = [&reader, &uniqueKmersPerSet, &final_readsPerSet, &final_readsPerSetRevComp,
							&mut, &klen, &includeRevComp, &finalSeqOut, &doNotWriteUndetermined,
							&hardCountOff, &excludeSetNames, &writeOutExclude,
							&smallLenCutOff,
							&finalSmallLenCutOffCount
							//, &testReads
			]() {
				//readInComp = [&reader, &uniqueKmersPerSet, &uniqueKmersFoundPerSet,&kmersFoundPerSeq,&mut,&klen,&includeRevComp]() {

				SimpleKmerHash hasher;
				PairedRead pseq;
				std::unordered_map<std::string, uint32_t> readsPerSetCurrent;
				std::unordered_map<std::string, uint32_t> readsPerSetRevCompCurrent;
				uint32_t currentFinalSmallLenCutOffCount = 0;

				while(reader.readNextReadLock(pseq)){
					if(len(pseq.seqBase_.seq_) < smallLenCutOff && len(pseq.mateSeqBase_.seq_) < smallLenCutOff){
						++currentFinalSmallLenCutOffCount;
						continue;
					}
					std::unordered_map<uint64_t, uint64_t> hashedInputKmers;
					std::unordered_map<uint64_t, uint64_t> hashedInputKmersRevComp;
					if(len(pseq.seqBase_.seq_) > klen){
						for(uint32_t pos = 0; pos < len(pseq.seqBase_.seq_) - klen + 1; ++pos){
							auto hash = hasher.hash(pseq.seqBase_.seq_.substr(pos, klen));
							++hashedInputKmers[hash];
						}
					}
					if(len(pseq.mateSeqBase_.seq_) > klen){
						for(uint32_t pos = 0; pos < len(pseq.mateSeqBase_.seq_) - klen + 1; ++pos){
							auto hash = hasher.hash(pseq.mateSeqBase_.seq_.substr(pos, klen));
							++hashedInputKmers[hash];
						}
					}
					if(includeRevComp){
						if(len(pseq.seqBase_.seq_) > klen){
							for(uint32_t pos = 0; pos < len(pseq.seqBase_.seq_) - klen + 1; ++pos){
								auto hash = hasher.revCompHash(pseq.seqBase_.seq_.substr(pos, klen));
								++hashedInputKmersRevComp[hash];
							}
						}
						if(len(pseq.mateSeqBase_.seq_) > klen){
							for(uint32_t pos = 0; pos < len(pseq.mateSeqBase_.seq_) - klen + 1; ++pos){
								auto hash = hasher.revCompHash(pseq.mateSeqBase_.seq_.substr(pos, klen));
								++hashedInputKmersRevComp[hash];
							}
						}
					}

					std::unordered_map<std::string, uint32_t> foundPerSet;
					std::unordered_map<std::string, uint32_t> foundPerSetRevComp;
					for(const auto & setName  : njh::getVecOfMapKeys(uniqueKmersPerSet)){
						foundPerSet[setName] = 0;
						foundPerSetRevComp[setName] = 0;
					}
					for(const auto & hashedKmer : hashedInputKmers){
						for(const auto & uniqueKmers : uniqueKmersPerSet){
							if(njh::in(hashedKmer.first, uniqueKmers.second)){
								++foundPerSet[uniqueKmers.first];
							}
						}
					}

					if(includeRevComp){
						for(const auto & hashedKmer : hashedInputKmersRevComp){
							for(const auto & uniqueKmers : uniqueKmersPerSet){
								if(njh::in(hashedKmer.first, uniqueKmers.second)){
									++foundPerSetRevComp[uniqueKmers.first];
								}
							}
						}
					}
					//set a hard cut off, reset counts to zero
					for(auto & perSet : foundPerSet){
						if(perSet.second < hardCountOff){
							perSet.second = 0;
						}
					}
					for(auto & perSet : foundPerSetRevComp){
						if(perSet.second < hardCountOff){
							perSet.second = 0;
						}
					}
					std::string winnerSet = "undetermined";
					double bestFrac = 0;
					bool winnerRevComp = false;

					for(const auto & setName  : njh::getVecOfMapKeys(uniqueKmersPerSet)){
						if(static_cast<double>(foundPerSet[setName])/static_cast<double>(hashedInputKmers.size()) > bestFrac){
							bestFrac = static_cast<double>(foundPerSet[setName])/static_cast<double>(hashedInputKmers.size());
							winnerSet = setName;
							winnerRevComp = false;
						}
						if (includeRevComp) {
							if(static_cast<double>(foundPerSetRevComp[setName])/static_cast<double>(hashedInputKmers.size()) > bestFrac){
								bestFrac = static_cast<double>(foundPerSetRevComp[setName])/static_cast<double>(hashedInputKmers.size());
								winnerSet = setName;
								winnerRevComp = true;
							}
						}
					}
//					if(njh::in(pseq.seqBase_.name_, testReads)){
//						for(const auto & setName  : njh::getVecOfMapKeys(uniqueKmersPerSet)) {
//							std::cout << "final"
//												<< "\t" << pseq.seqBase_.name_
//												<< "\t" << hashedInputKmers.size()
//												<< "\t" << setName
//												<< "\t" << uniqueKmersPerSet[setName].size()
//												<< "\t" << foundPerSet[setName]
//												<< "\t"
//												<< static_cast<double>(foundPerSet[setName]) / static_cast<double>(hashedInputKmers.size())
//												<< "\t" << static_cast<double>(foundPerSet[setName]) /
//																	 static_cast<double>(uniqueKmersPerSet[setName].size());
//							if (includeRevComp) {
//								std::cout << "\t" << foundPerSetRevComp[setName]
//													<< "\t" << static_cast<double>(foundPerSetRevComp[setName]) /
//																		 static_cast<double>(hashedInputKmers.size())
//													<< "\t" << static_cast<double>(foundPerSetRevComp[setName]) /
//																		 static_cast<double>(uniqueKmersPerSet[setName].size());
//							}
//							std::cout << "\t" << winnerSet
//												<< "\t" << bestFrac
//												<< "\t" << njh::boolToStr(winnerRevComp);
//							std::cout << std::endl;
//						}
//					}
					if (winnerRevComp) {
						++readsPerSetRevCompCurrent[winnerSet];
						pseq.seqBase_.reverseComplementRead(false, true);
						pseq.mateSeqBase_.reverseComplementRead(false, true);
					} else {
						++readsPerSetCurrent[winnerSet];
					}
					if(!doNotWriteUndetermined || winnerSet != "undetermined"){
						if(writeOutExclude || !njh::in(winnerSet, excludeSetNames)){
							finalSeqOut.openWrite(winnerSet, pseq);
						}
					}
				}
				{
					std::lock_guard<std::mutex> lockGuard(mut);
					finalSmallLenCutOffCount += currentFinalSmallLenCutOffCount;
					for(const auto & readsPerSetCount : readsPerSetCurrent){
						final_readsPerSet[readsPerSetCount.first] += readsPerSetCount.second;
					}
					for(const auto & readsPerSetRevCompCount : readsPerSetRevCompCurrent){
						final_readsPerSetRevComp[readsPerSetRevCompCount.first] += readsPerSetRevCompCount.second;
					}
				}
			};
		} else {
			for(const auto & name : names){
				auto seqOutOpts = SeqIOOptions::genFastqOutGz(njh::files::make_path(finalExtractionDir, name) );
				finalSeqOut.addReader(name, seqOutOpts);
			}
			finalSeqOut.addReader("undetermined", SeqIOOptions::genFastqOutGz(njh::files::make_path(finalExtractionDir, "undetermined")));
			readInComp = [&reader, &uniqueKmersPerSet, &final_readsPerSet,&final_readsPerSetRevComp,
							&mut,&klen,&includeRevComp,&finalSeqOut,&doNotWriteUndetermined,
							&hardCountOff, &excludeSetNames, &writeOutExclude,
							&smallLenCutOff,
							&finalSmallLenCutOffCount]() {
//		readInComp = [&reader, &uniqueKmersPerSet, &uniqueKmersFoundPerSet,&kmersFoundPerSeq,&mut,&klen,&includeRevComp]() {
				SimpleKmerHash hasher;
				seqInfo seq;
				std::unordered_map<std::string, uint32_t> readsPerSetCurrent;
				std::unordered_map<std::string, uint32_t> readsPerSetRevCompCurrent;
				uint32_t currentFinalSmallLenCutOffCount = 0;
				while(reader.readNextReadLock(seq)){
					if(len(seq) < smallLenCutOff){
						++currentFinalSmallLenCutOffCount;
					}
					std::unordered_map<uint64_t, uint64_t> hashedInputKmers;
					std::unordered_map<uint64_t, uint64_t> hashedInputKmersRevComp;
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
								++hashedInputKmersRevComp[hash];
							}
						}
					}

					std::unordered_map<std::string, uint32_t> foundPerSet;
					std::unordered_map<std::string, uint32_t> foundPerSetRevComp;
					for(const auto & setName  : njh::getVecOfMapKeys(uniqueKmersPerSet)){
						foundPerSet[setName] = 0;
						foundPerSetRevComp[setName] = 0;
					}
					for(const auto & hashedKmer : hashedInputKmers){
						for(const auto & uniqueKmers : uniqueKmersPerSet){
							if(njh::in(hashedKmer.first, uniqueKmers.second)){
								++foundPerSet[uniqueKmers.first];
							}
						}
					}
					if(includeRevComp){
						for(const auto & hashedKmer : hashedInputKmersRevComp){
							for(const auto & uniqueKmers : uniqueKmersPerSet){
								if(njh::in(hashedKmer.first, uniqueKmers.second)){
									++foundPerSetRevComp[uniqueKmers.first];
								}
							}
						}
					}
					//set a hard cut off, reset counts to zero
					for(auto & perSet : foundPerSet){
						if(perSet.second < hardCountOff){
							perSet.second = 0;
						}
					}
					for(auto & perSet : foundPerSetRevComp){
						if(perSet.second < hardCountOff){
							perSet.second = 0;
						}
					}
					std::string winnerSet = "undetermined";
					double bestFrac = 0;
					bool winnerRevComp = false;

					for(const auto & setName  : njh::getVecOfMapKeys(uniqueKmersPerSet)){
						if(static_cast<double>(foundPerSet[setName])/static_cast<double>(hashedInputKmers.size()) > bestFrac){
							bestFrac = static_cast<double>(foundPerSet[setName])/static_cast<double>(hashedInputKmers.size());
							winnerSet = setName;
							winnerRevComp = false;
						}
						if (includeRevComp) {
							if(static_cast<double>(foundPerSetRevComp[setName])/static_cast<double>(hashedInputKmers.size()) > bestFrac){
								bestFrac = static_cast<double>(foundPerSetRevComp[setName])/static_cast<double>(hashedInputKmers.size());
								winnerSet = setName;
								winnerRevComp = true;
							}
						}
					}
					if(winnerRevComp){
						++readsPerSetRevCompCurrent[winnerSet];
						seq.reverseComplementRead(true, true);
					}else{
						++readsPerSetCurrent[winnerSet];
					}
					if(!doNotWriteUndetermined || winnerSet != "undetermined"){
						if(writeOutExclude || !njh::in(winnerSet, excludeSetNames)){
							finalSeqOut.openWrite(winnerSet, seq);
						}
					}
				}
				{
					std::lock_guard<std::mutex> lockGuard(mut);
					finalSmallLenCutOffCount += currentFinalSmallLenCutOffCount;
					for(const auto & readsPerSetCount : readsPerSetCurrent){
						final_readsPerSet[readsPerSetCount.first] += readsPerSetCount.second;
					}
					for(const auto & readsPerSetRevCompCount : readsPerSetRevCompCurrent){
						final_readsPerSetRevComp[readsPerSetRevCompCount.first] += readsPerSetRevCompCount.second;
					}
				}
			};
		}
		njh::concurrent::runVoidFunctionThreaded(readInComp, numThreads);


		uint64_t totalReads = 0;
		for(const auto & readsPerSetCount : final_readsPerSet){
			totalReads += readsPerSetCount.second;
		}
		for(const auto & readsPerSetRevCompCount : final_readsPerSetRevComp){
			totalReads += readsPerSetRevCompCount.second;
		}

		if (includeRevComp) {
			for(const auto & setName : uniqueKmersPerSet){
				uint64_t totalExtracted = final_readsPerSet[setName.first] + final_readsPerSetRevComp[setName.first];
				outCounts << "final"
									<< "\t" << sampleName
									<< "\t" << totalReads
									<< "\t" << setName.first
									<< "\t" << totalExtracted
									<< "\t" << static_cast<double>(totalExtracted) / static_cast<double>(totalReads)
									<< "\t" << final_readsPerSet[setName.first]
									<< "\t" << (totalExtracted > 0 ? static_cast<double>(final_readsPerSet[setName.first]) / static_cast<double>(totalExtracted) : 0)
									<< std::endl;
			}
			{
				uint64_t totalExtracted = final_readsPerSet["undetermined"] + final_readsPerSetRevComp["undetermined"];
				outCounts << "final"
									<< "\t" << sampleName
									<< "\t" << totalReads
									<< "\t" << "undetermined"
									<< "\t" << totalExtracted
									<< "\t" << static_cast<double>(totalExtracted) / static_cast<double>(totalReads)
									<< "\t" << final_readsPerSet["undetermined"]
									<< "\t" << (totalExtracted > 0 ?static_cast<double>(final_readsPerSet["undetermined"]) / static_cast<double>(totalExtracted) : 0)
									<< std::endl;
			}
			{
				outCounts << "final"
									<< "\t" << sampleName
									<< "\t" << totalReads
									<< "\t" << "smallLenReads"
									<< "\t" << finalSmallLenCutOffCount
									<< "\t" << static_cast<double>(finalSmallLenCutOffCount) / static_cast<double>(totalReads)
									<< "\t" << finalSmallLenCutOffCount
									<< "\t" << 0
									<< std::endl;
			}
		} else {

			for(const auto & setName : uniqueKmersPerSet){
				outCounts << "final"
									<< "\t" << sampleName
									<< "\t" << totalReads
									<< "\t" << setName.first
									<< "\t" << final_readsPerSet[setName.first]
									<< "\t" << static_cast<double>(final_readsPerSet[setName.first]) / static_cast<double>(totalReads)
									<< std::endl;
			}
			{
				outCounts << "final"
									<< "\t" << sampleName
									<< "\t" << totalReads
									<< "\t" << "undetermined"
									<< "\t" << final_readsPerSet["undetermined"]
									<< "\t" << static_cast<double>(final_readsPerSet["undetermined"]) / static_cast<double>(totalReads)
									<< std::endl;
			}
			{
				outCounts << "final"
									<< "\t" << sampleName
									<< "\t" << totalReads
									<< "\t" << "smallLenReads"
									<< "\t" << finalSmallLenCutOffCount
									<< "\t" << static_cast<double>(finalSmallLenCutOffCount) / static_cast<double>(totalReads)
									<< std::endl;
			}
		}
		if(setUp.pars_.verbose_){
			std::cout << "iterNumber: " << iterNumber << std::endl;
			std::cout << "totalReads: " << totalReads << std::endl;
			std::cout << "totalDeterminedReads: " << totalReads - final_readsPerSet["undetermined"] + final_readsPerSetRevComp["undetermined"] << std::endl;
			std::cout << "totalUndeterminedReads: " << final_readsPerSet["undetermined"] + final_readsPerSetRevComp["undetermined"] << std::endl;
		}

		if(!keepTemporaryFiles){
			//remove all fast[aq] files and fast[aq].gz files within directory
			auto readFiles = njh::files::listAllFiles(setUp.pars_.directoryName_,false,{std::regex{R"((.*)(\.(fast[aq])(\.gz)?)$)"}});
			for(const auto & f : readFiles){
				bfs::remove(f.first);
			}
		}
	}
	watch.logLapTimes(setUp.rLog_.runLogFile_, true, 6, true);
	return 0;
}


}  // namespace njhseq

