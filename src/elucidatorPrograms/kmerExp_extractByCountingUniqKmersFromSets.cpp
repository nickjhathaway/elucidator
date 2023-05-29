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
	UniqueKmerSetHelper::CompareReadToSetPars compPars;
	OutOptions outOpts(bfs::path(""), ".tab.txt.gz");
	seqSetUp setUp(inputCommands);
	//default options
	setUp.processVerbose();
	setUp.processDebug();
	//read in options
	setUp.pars_.ioOptions_.revComplMate_ = true;
	bool pairedEndSet = setUp.processReadInNames(njhseq::seqSetUp::pairedReadInFormatsAvailable_, false);
	auto singlesOption = SeqIOOptions();
	setUp.processJustReadInNames(singlesOption, njhseq::seqSetUp::singleInFormatsAvailable_, !pairedEndSet);
	singlesOption.lowerCaseBases_ = setUp.pars_.ioOptions_.lowerCaseBases_;
	singlesOption.removeGaps_ = setUp.pars_.ioOptions_.removeGaps_;
	singlesOption.includeWhiteSpaceInName_ = setUp.pars_.ioOptions_.includeWhiteSpaceInName_;

	//comp results
	setUp.setOption(countTable, "--kmerTable", "kmerTable, no header, columns 1)set,2)kmer", true);
	setUp.setOption(compPars.sampleName, "--sampleName", "Name to add to output file", true);
	setUp.setOption(compPars.includeRevComp, "--includeRevComp", "include Rev Comp of the input seqs");
	setUp.setOption(compPars.pairsSeparate, "--pairsSeparate", "count the pairs separately");
	//run options
	setUp.setOption(numThreads, "--numThreads", "numThreads");
	//write options
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);
	njh::stopWatch watch;
	watch.setLapName("initial");

	watch.startNewLap("reading in unique kmer table");
	compPars.klen = UniqueKmerSetHelper::getKmerLenFromUniqueKmerTable(countTable);
	std::unordered_map<std::string, std::unordered_set<uint64_t>> uniqueKmersPerSet = UniqueKmerSetHelper::readInUniqueKmerTablePerSet(countTable);
	if(setUp.pars_.verbose_){
		std::cout << watch.getLapName() << "\t" << watch.timeLapFormatted() <<std::endl;
	}

	OutputStream out(outOpts);
	UniqueKmerSetHelper::CompareReadToSetRes::writeOutputHeader(out, compPars);

	std::mutex outMut;
	std::function<void()> readInComp;
	if(!setUp.pars_.ioOptions_.firstName_.empty()){
		SeqInput reader(setUp.pars_.ioOptions_);
		reader.openIn();
		readInComp = [&reader, &uniqueKmersPerSet,&compPars,&out,&outMut]() {
			SimpleKmerHash hasher;
			PairedRead pseq;
			while(reader.readNextReadLock(pseq)){
				if(compPars.pairsSeparate){
					auto compResFirstMate = UniqueKmerSetHelper::compareReadToSetRes(pseq.seqBase_, uniqueKmersPerSet, compPars, hasher);
					pseq.seqBase_.name_.append("_firstMate");
					auto compResSecondMate = UniqueKmerSetHelper::compareReadToSetRes(pseq.mateSeqBase_, uniqueKmersPerSet, compPars, hasher);
					pseq.mateSeqBase_.name_.append("_secondMate");
					{
						std::lock_guard<std::mutex> lockGuard(outMut);
						compResFirstMate.writeOutput(out,pseq.seqBase_,uniqueKmersPerSet, compPars);
						compResSecondMate.writeOutput(out,pseq.mateSeqBase_,uniqueKmersPerSet, compPars);
					}
				}else{
					auto compRes = UniqueKmerSetHelper::compareReadToSetRes(pseq, uniqueKmersPerSet, compPars, hasher);
					{
						std::lock_guard<std::mutex> lockGuard(outMut);
						compRes.writeOutput(out,pseq.seqBase_,uniqueKmersPerSet, compPars);
					}
				}
			}
		};
		njh::concurrent::runVoidFunctionThreaded(readInComp, numThreads);
	}
	if(!singlesOption.firstName_.empty()){
		SeqInput reader(singlesOption);
		reader.openIn();
		readInComp = [&reader, &uniqueKmersPerSet,&compPars,&out,&outMut]() {
			SimpleKmerHash hasher;
			seqInfo seq;
			while(reader.readNextReadLock(seq)){
				auto compRes = UniqueKmerSetHelper::compareReadToSetRes(seq, uniqueKmersPerSet, compPars, hasher);
				{
					std::lock_guard<std::mutex> lockGuard(outMut);
					compRes.writeOutput(out,seq,uniqueKmersPerSet, compPars);
				}
			}
		};
		njh::concurrent::runVoidFunctionThreaded(readInComp, numThreads);
	}
	return 0;
}


int kmerExpRunner::extractByCountingUniqKmersFromSets(const njh::progutils::CmdArgs & inputCommands){

	uint32_t numThreads = 1;
	bfs::path countTable;
	bfs::path excludesCountTable;
	bfs::path nonUniqueKmerTable;
	uint32_t maxIterations = std::numeric_limits<uint32_t>::max();

	bool keepTemporaryFiles = false;
	bool doNotDoFinalExtract = false;

	UniqueKmerSetHelper::ProcessReadForExtractingPars extractingPars;
	extractingPars.compPars.hardCountOff = 30;

	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	//read in options
	setUp.pars_.ioOptions_.revComplMate_ = true;
	bool pairedEndSet = setUp.processReadInNames(njhseq::seqSetUp::pairedReadInFormatsAvailable_, false);
	auto singlesOption = SeqIOOptions();
	setUp.processJustReadInNames(singlesOption, njhseq::seqSetUp::singleInFormatsAvailable_, !pairedEndSet);
	singlesOption.lowerCaseBases_ = setUp.pars_.ioOptions_.lowerCaseBases_;
	singlesOption.removeGaps_ = setUp.pars_.ioOptions_.removeGaps_;
	singlesOption.includeWhiteSpaceInName_ = setUp.pars_.ioOptions_.includeWhiteSpaceInName_;

	setUp.setOption(countTable, "--kmerTable", "countTable, 1)set,2)kmer", true);
	setUp.setOption(extractingPars.compPars.sampleName, "--sampleName", "Name to add to output file", true);
	setUp.setOption(extractingPars.compPars.includeRevComp, "--includeRevComp", "include Rev Comp of the input seqs");
	setUp.setOption(extractingPars.compPars.pairsSeparate, "--pairsSeparate", "count the pairs separately");
	if(bfs::exists(countTable)){
		auto row = tokenizeString(njh::files::getFirstLine(countTable), "\t");
		extractingPars.compPars.klen = row[1].size();
	}
	setUp.setOption(extractingPars.compPars.hardCountOff, "--hardCountOff", "hard Count Off, do not count sets unless greater thant his number");
	extractingPars.smallLenCutOff = extractingPars.compPars.hardCountOff + extractingPars.compPars.klen + 11;
	setUp.setOption(extractingPars.smallLenCutOff, "--smallLenCutOff", "small Len Cut Off of input sequences, if less than this size will skip over");

	setUp.setOption(excludesCountTable, "--excludesCountTable", "excludes Count Table, 1)set,2)kmer");
	setUp.setOption(extractingPars.excludeSetNames, "--excludeSetNames", "names of sets of unqiue kmers to ignore from the --kmerTable");
	//
	setUp.setOption(keepTemporaryFiles, "--keepTemporaryFiles", "Keep Temporary Files");
	setUp.setOption(extractingPars.addingInKmersCountCutOff, "--addingInKmersCountCutOff", "adding In Kmers Count Cut Off");


	setUp.setOption(extractingPars.doReCheckExcludeSets, "--doReCheckExcludeSets", "do Re Check Exclude Sets on iterations, could lead to some reads being recruit away but will drastically speed up run time");

	bool doNotFinalRecheckOfExcludeSeqs = false;
	setUp.setOption(doNotFinalRecheckOfExcludeSeqs, "--doNotFinalRecheckOfExcludeSeqs", "redo a unique kmer filter on the exclusion sets during the final re-extraction process");

	bool doFinalRecheckOfExcludeSeqs = !doNotFinalRecheckOfExcludeSeqs;
//	setUp.setOption(doFinalRecheckOfExcludeSeqs, "--doFinalRecheckOfExcludeSeqs", "Normal the exclusion kmer set is left alone, if this flag is set, it will be filter during the final re-extraction process");

	setUp.setOption(doNotDoFinalExtract, "--doNotDoFinalExtract", "do Not Do Final Extract");
	setUp.setOption(nonUniqueKmerTable, "--nonUniqueKmerTable", "non-unique Kmer Table, 1)set,2)kmer");
	//setUp.setOption(extractAfterIterating, "--extractAfterIterating", "Extract After Iterating");
	setUp.setOption(maxIterations, "--maxIterations", "max Iterations to perform");


	setUp.setOption(extractingPars.doNotWriteUndetermined, "--doNotWriteUndetermined", "do Not Write Undetermined");
	setUp.setOption(extractingPars.writeOutExclude, "--writeOutExclude", "write Out Excluded reads that match the excluded kmer sets");



	setUp.setOption(numThreads, "--numThreads", "numThreads");
	setUp.processDirectoryOutputName(true);

	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	njh::stopWatch watch;
	watch.setLapName("initial");


	UniqueKmerSetHelper::ProcessReadForExtractingCounts initialCounts;

	std::mutex mut;

	watch.startNewLap("reading in unique kmer table");


	std::unordered_map<std::string, std::unordered_set<uint64_t>> uniqueKmersPerSet = UniqueKmerSetHelper::readInUniqueKmerTablePerSet(
					countTable);
	std::string nonUniqueRegionName = "NON_UNIQUE";
	std::unordered_set<uint64_t> nonUniqueKmersPerSet;
	if (!nonUniqueKmerTable.empty()) {
		nonUniqueKmersPerSet = UniqueKmerSetHelper::readInUniqueKmerTableSetsCollapsed(nonUniqueKmerTable);
	}

	if(setUp.pars_.verbose_){
		std::cout << watch.getLapName() << "\t" << watch.timeLapFormatted() <<std::endl;
		std::cout << "klen: " << extractingPars.compPars.klen << std::endl;
		for(const auto & uni : uniqueKmersPerSet){
			std::cout << uni.first << "\t" << uni.second.size() << std::endl;
		}
		std::cout << "NON_UNIQUE" << "\t" << nonUniqueKmersPerSet.size() << std::endl;
	}
	MultiSeqIO initialSeqOut;
	watch.startNewLap("initial scan");

	OutOptions outCountsOpts(njh::files::make_path(setUp.pars_.directoryName_, "extractionCounts.tab.txt"));
	OutputStream outCounts(outCountsOpts);
	UniqueKmerSetHelper::ProcessReadForExtractingCounts::writeOutCountsHeader(outCounts, extractingPars);

	{//initial

		//add outputs
		VecStr names = getVectorOfMapKeys(uniqueKmersPerSet);
		for(const auto & name : names){
			auto seqOutOpts = SeqIOOptions::genPairedOutGz(njh::files::make_path(setUp.pars_.directoryName_, name));
			initialSeqOut.addReader(njh::pasteAsStr(name, "-paired"), seqOutOpts);
		}
		initialSeqOut.addReader("undetermined-paired", SeqIOOptions::genPairedOutGz(njh::files::make_path(setUp.pars_.directoryName_, "undetermined")));
		for(const auto & name : names){
			auto seqOutOpts = SeqIOOptions::genFastqOutGz(njh::files::make_path(setUp.pars_.directoryName_, name));
			initialSeqOut.addReader(njh::pasteAsStr(name, "-single"), seqOutOpts);
		}
		initialSeqOut.addReader("undetermined-single", SeqIOOptions::genFastqOutGz(njh::files::make_path(setUp.pars_.directoryName_, "undetermined")));

		//to make below count correctly, make sure mate is automatically reversed complemented
		setUp.pars_.ioOptions_.revComplMate_ = true;

		if (!setUp.pars_.ioOptions_.firstName_.empty()) {
			SeqInput reader(setUp.pars_.ioOptions_);
			reader.openIn();
			std::function<void()> readInComp = [&reader, &uniqueKmersPerSet,
							&mut, &extractingPars, &initialSeqOut, &initialCounts
//							, &testReads
			]() {
				SimpleKmerHash hasher;
				PairedRead pseq;
				UniqueKmerSetHelper::ProcessReadForExtractingCounts currentCounts;
				while (reader.readNextReadLock(pseq)) {
					UniqueKmerSetHelper::processReadForExtracting(pseq, uniqueKmersPerSet,
																												extractingPars, hasher, initialSeqOut,
																												currentCounts);
				}
				{
					std::lock_guard<std::mutex> lockGuard(mut);
					initialCounts.addOtherCounts(currentCounts);
				}
			};
			njh::concurrent::runVoidFunctionThreaded(readInComp, numThreads);
		}
		if (!singlesOption.firstName_.empty()) {
			SeqInput reader(singlesOption);
			reader.openIn();
			std::function<void()> readInComp = [&reader, &uniqueKmersPerSet,
							&mut, &extractingPars, &initialSeqOut, &initialCounts
//							, &testReads
			]() {
				SimpleKmerHash hasher;
				seqInfo seq;
				UniqueKmerSetHelper::ProcessReadForExtractingCounts currentCounts;
				while (reader.readNextReadLock(seq)) {
					UniqueKmerSetHelper::processReadForExtracting(seq, uniqueKmersPerSet,
																												extractingPars, hasher, initialSeqOut,
																												currentCounts);
				}
				{
					std::lock_guard<std::mutex> lockGuard(mut);
					initialCounts.addOtherCounts(currentCounts);
				}
			};
			njh::concurrent::runVoidFunctionThreaded(readInComp, numThreads);

		}
		initialCounts.writeOutCounts(outCounts,extractingPars,uniqueKmersPerSet,"initial");
	}

	//
	if(setUp.pars_.verbose_){
		std::cout << "iterNumber: " << "initial" << std::endl;
		std::cout << "Total Initial Reads: " << initialCounts.getTotalCounts() << std::endl;
		std::cout << "Total Initial Undetermined Reads: " << initialCounts.genTotalUndeterminedCount() << " (" << initialCounts.genTotalUndeterminedCount() * 100 /static_cast<double>(initialCounts.getTotalCounts())<< "%)" << std::endl;
		std::cout << "Total Initial Determined Reads: " << initialCounts.genTotalDeterminedCount() << " (" << initialCounts.genTotalDeterminedCount() * 100 /static_cast<double>(initialCounts.getTotalCounts())<< "%)"<< std::endl;
		std::cout << "Total Initial Less Than Small Len Cut Off Count: " << initialCounts.smallLenCutOffCount << " (" << initialCounts.smallLenCutOffCount * 100 /static_cast<double>(initialCounts.getTotalCounts())<< "%)"<< std::endl;

	}
	//read in the extract reads and add to the unique kmer sets
	bool iterate = true;
	uint64_t totalUndeterminedExtracted = initialCounts.genTotalUndeterminedCount();
	uint32_t iterNumber = 0;
	if(0 == totalUndeterminedExtracted){
		iterate = false;
	}
	auto original_initialCounts = initialCounts;
	initialCounts.readsPerSet["undetermined"] = 0;
	initialCounts.readsPerSetRevComp["undetermined"] = 0;
	uint64_t totalDeterminedReads = initialCounts.genTotalDeterminedCount();
//	std::cout << __FILE__ << " " << __LINE__ << std::endl;
	while(iterate && iterNumber < maxIterations){
//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		watch.startNewLap(njh::pasteAsStr(iterNumber, "- read in new kmers"));
//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		initialSeqOut.closeOutForReopeningAll();
		auto rawUniqKmerSets = njh::getVecOfMapKeys(uniqueKmersPerSet);
		removeElement(rawUniqKmerSets, std::string("undetermined"));
		std::unordered_map<std::string, std::unordered_map<uint64_t, uint32_t>> rawKmersPerInput = UniqueKmerSetHelper::readInNewKmersFromExtractedReads(
						setUp.pars_.directoryName_, rawUniqKmerSets, extractingPars);

//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		watch.startNewLap(njh::pasteAsStr(iterNumber, "- check new kmers against current set"));
		std::unordered_map<std::string, std::unordered_set<uint64_t>> outputUniqueKmersPerSet = UniqueKmerSetHelper::filterReExtractedKmersForNonUnique(rawKmersPerInput,
																																																																										extractingPars,
																																																																										uniqueKmersPerSet,
																																																																										nonUniqueKmersPerSet);
//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		if(!extractingPars.doReCheckExcludeSets){
			for(const auto & excludeName : extractingPars.excludeSetNames){
				outputUniqueKmersPerSet[excludeName] = std::move(uniqueKmersPerSet[excludeName]);
			}
		}
//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		watch.startNewLap(njh::pasteAsStr(iterNumber, "- copy new set"));
		uniqueKmersPerSet = std::move(outputUniqueKmersPerSet);
		watch.startNewLap(njh::pasteAsStr(iterNumber, "- scan again"));
//		std::cout << __FILE__ << " " << __LINE__ << std::endl;

		bfs::path orig_undeterminedInput = njh::files::make_path(setUp.pars_.directoryName_, "undetermined.fastq.gz");
		bfs::path orig_undeterminedInputR1 = njh::files::make_path(setUp.pars_.directoryName_, "undetermined_R1.fastq.gz");
		bfs::path orig_undeterminedInputR2 = njh::files::make_path(setUp.pars_.directoryName_, "undetermined_R2.fastq.gz");

		bfs::path undeterminedInput = njh::files::make_path(setUp.pars_.directoryName_, njh::pasteAsStr(iterNumber, "_undetermined.fastq.gz"));
		bfs::path undeterminedInputR1 = njh::files::make_path(setUp.pars_.directoryName_, njh::pasteAsStr(iterNumber, "_undetermined_R1.fastq.gz"));
		bfs::path undeterminedInputR2 = njh::files::make_path(setUp.pars_.directoryName_, njh::pasteAsStr(iterNumber, "_undetermined_R2.fastq.gz"));

		SeqIOOptions undeterminedPairedInOpts;
		SeqIOOptions undeterminedSingleInOpts;
		bool pairedInUndeterminedExists = bfs::exists(orig_undeterminedInputR1);
		bool singleInUndeterminedExists = bfs::exists(orig_undeterminedInput);
		if (pairedInUndeterminedExists) {
			bfs::rename(orig_undeterminedInputR1, undeterminedInputR1);
			bfs::rename(orig_undeterminedInputR2, undeterminedInputR2);
			undeterminedPairedInOpts.firstName_ = undeterminedInputR1;
			undeterminedPairedInOpts.secondName_ = undeterminedInputR2;
			undeterminedPairedInOpts.inFormat_ = SeqIOOptions::inFormats::FASTQPAIREDGZ;
			undeterminedPairedInOpts.revComplMate_ = true;
			//to make below count correctly, make sure mate is automatically reversed complemented
		}
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		if(singleInUndeterminedExists) {
			bfs::rename(orig_undeterminedInput, undeterminedInput);
			undeterminedSingleInOpts.firstName_ = undeterminedInput;
			undeterminedSingleInOpts.inFormat_ = SeqIOOptions::inFormats::FASTQGZ;
		}
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		if (pairedInUndeterminedExists) {
			//std::cout << __FILE__ << " " << __LINE__ << std::endl;
			//re-extract paired in
			SeqInput reader(undeterminedPairedInOpts);
			reader.openIn();
			//std::cout << __FILE__ << " " << __LINE__ << std::endl;
			std::function<void()> readInComp = [&reader, &uniqueKmersPerSet,
							&mut, &extractingPars, &initialSeqOut, &initialCounts
			]() {
				//std::cout << __FILE__ << " " << __LINE__ << std::endl;
				SimpleKmerHash hasher;
				PairedRead pseq;
				UniqueKmerSetHelper::ProcessReadForExtractingCounts currentCounts;
				//std::cout << __FILE__ << " " << __LINE__ << std::endl;
				while (reader.readNextReadLock(pseq)) {
					//std::cout << __FILE__ << " " << __LINE__ << std::endl;
					UniqueKmerSetHelper::processReadForExtracting(pseq, uniqueKmersPerSet,
																												extractingPars, hasher, initialSeqOut,
																												currentCounts);
					//std::cout << __FILE__ << " " << __LINE__ << std::endl;
				}
				{
					std::lock_guard<std::mutex> lockGuard(mut);
					//std::cout << __FILE__ << " " << __LINE__ << std::endl;
					initialCounts.addOtherCounts(currentCounts);
					//std::cout << __FILE__ << " " << __LINE__ << std::endl;

				}
			};
			//std::cout << __FILE__ << " " << __LINE__ << std::endl;
			njh::concurrent::runVoidFunctionThreaded(readInComp, numThreads);
			//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		}
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		if(singleInUndeterminedExists) {
			//std::cout << __FILE__ << " " << __LINE__ << std::endl;
			//re-extract single in
			SeqInput reader(undeterminedSingleInOpts);
			reader.openIn();
			//std::cout << __FILE__ << " " << __LINE__ << std::endl;
			std::function<void()> readInComp = [&reader, &uniqueKmersPerSet,
							&mut, &extractingPars, &initialSeqOut, &initialCounts
//							, &testReads
			]() {
				SimpleKmerHash hasher;
				seqInfo seq;
				UniqueKmerSetHelper::ProcessReadForExtractingCounts currentCounts;
				while (reader.readNextReadLock(seq)) {
					UniqueKmerSetHelper::processReadForExtracting(seq, uniqueKmersPerSet,
																												extractingPars, hasher, initialSeqOut,
																												currentCounts);
				}
				{
					std::lock_guard<std::mutex> lockGuard(mut);
					initialCounts.addOtherCounts(currentCounts);
				}
			};
			//std::cout << __FILE__ << " " << __LINE__ << std::endl;
			njh::concurrent::runVoidFunctionThreaded(readInComp, numThreads);
			//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		}
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		initialCounts.writeOutCounts(outCounts,extractingPars,uniqueKmersPerSet,njh::pasteAsStr(iterNumber));
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		uint64_t currentTotalUndeterminedReads = initialCounts.genTotalUndeterminedCount();
		uint64_t currentTotalReads = initialCounts.getTotalCounts();
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		initialCounts.readsPerSet["undetermined"] = 0;
		initialCounts.readsPerSetRevComp["undetermined"] = 0;
		uint64_t currentTotalDeterminedReads = initialCounts.genTotalDeterminedCount();
		if(setUp.pars_.verbose_){

			std::cout << "iterNumber: " << iterNumber << std::endl;
			std::cout << "Determined: " << totalDeterminedReads - currentTotalDeterminedReads << " more reads" <<std::endl;
			std::cout << "Current counts: " << std::endl;
			std::cout << "\tTotal Reads: " << currentTotalReads << std::endl;
			std::cout << "\tTotal Undetermined Reads: " << currentTotalUndeterminedReads << " (" << currentTotalUndeterminedReads * 100 /static_cast<double>(currentTotalReads)<< "%)" << std::endl;
			std::cout << "\tTotal Determined Reads: " << currentTotalDeterminedReads << " (" << currentTotalDeterminedReads * 100 /static_cast<double>(currentTotalReads)<< "%)"<< std::endl;
			std::cout << "\tTotal Less Than Small Len Cut Off Count: " << initialCounts.smallLenCutOffCount << " (" << initialCounts.smallLenCutOffCount * 100 /static_cast<double>(currentTotalReads)<< "%)"<< std::endl;
			std::cout << std::endl;
		}
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		if(currentTotalUndeterminedReads == 0 || totalDeterminedReads == currentTotalDeterminedReads){
			//if there's no more undetermined reads or if the number determined this round was zero then stop iterating
			iterate = false;
		}
		totalDeterminedReads = currentTotalDeterminedReads;
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		if (!keepTemporaryFiles) {
			//std::cout << __FILE__ << " " << __LINE__ << std::endl;
			if (bfs::exists(undeterminedPairedInOpts.firstName_)) {
				bfs::remove(undeterminedPairedInOpts.firstName_);
				bfs::remove(undeterminedPairedInOpts.secondName_);
			}
			//std::cout << __FILE__ << " " << __LINE__ << std::endl;
			if (bfs::exists(undeterminedSingleInOpts.firstName_)) {
				bfs::remove(undeterminedSingleInOpts.firstName_);
			}
			//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		}
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		++iterNumber;
	}
	//std::cout << __FILE__ << " " << __LINE__ << std::endl;
//	if(extractAfterIterating){
	if(!doNotDoFinalExtract){
		initialSeqOut.closeOutForReopeningAll();

		if(!extractingPars.doReCheckExcludeSets && doFinalRecheckOfExcludeSeqs){
			watch.startNewLap(njh::pasteAsStr(iterNumber, "- final read in new kmers"));
			auto rawUniqKmerSets = njh::getVecOfMapKeys(uniqueKmersPerSet);
			removeElement(rawUniqKmerSets, std::string("undetermined"));
			std::unordered_map<std::string, std::unordered_map<uint64_t, uint32_t>> rawKmersPerInput = UniqueKmerSetHelper::readInNewKmersFromExtractedReads(
							setUp.pars_.directoryName_, rawUniqKmerSets, extractingPars);


			watch.startNewLap(njh::pasteAsStr(iterNumber, "- final check new kmers against current set"));
			std::unordered_map<std::string, std::unordered_set<uint64_t>> outputUniqueKmersPerSet = UniqueKmerSetHelper::filterReExtractedKmersForNonUniqueIncludeExcludedSets(
							rawKmersPerInput,
							extractingPars,
							uniqueKmersPerSet,
							nonUniqueKmersPerSet);

			watch.startNewLap(njh::pasteAsStr(iterNumber, "- final copy new set"));
			uniqueKmersPerSet = std::move(outputUniqueKmersPerSet);
		}

		UniqueKmerSetHelper::ProcessReadForExtractingCounts finalCounts;


		watch.startNewLap(njh::pasteAsStr(iterNumber, "- scan whole set again"));


		auto finalExtractionDir = njh::files::make_path(setUp.pars_.directoryName_, "finalExtraction");
		njh::files::makeDir(njh::files::MkdirPar{finalExtractionDir});

		//add outputs
		MultiSeqIO finalSeqOut;

		VecStr names = getVectorOfMapKeys(uniqueKmersPerSet);
		for(const auto & name : names){
			auto seqOutOpts = SeqIOOptions::genPairedOutGz(njh::files::make_path(finalExtractionDir, name));
			finalSeqOut.addReader(njh::pasteAsStr(name, "-paired"), seqOutOpts);
		}
		finalSeqOut.addReader("undetermined-paired", SeqIOOptions::genPairedOutGz(njh::files::make_path(finalExtractionDir, "undetermined")));
		for(const auto & name : names){
			auto seqOutOpts = SeqIOOptions::genFastqOutGz(njh::files::make_path(finalExtractionDir, name));
			finalSeqOut.addReader(njh::pasteAsStr(name, "-single"), seqOutOpts);
		}
		finalSeqOut.addReader("undetermined-single", SeqIOOptions::genFastqOutGz(njh::files::make_path(finalExtractionDir, "undetermined")));


		if (!setUp.pars_.ioOptions_.firstName_.empty()) {
			SeqInput reader(setUp.pars_.ioOptions_);
			reader.openIn();
			std::function<void()> readInComp = [&reader, &uniqueKmersPerSet,
							&mut, &extractingPars, &finalSeqOut, &finalCounts
//							, &testReads
			]() {
				SimpleKmerHash hasher;
				PairedRead pseq;
				UniqueKmerSetHelper::ProcessReadForExtractingCounts currentCounts;
				while (reader.readNextReadLock(pseq)) {
					UniqueKmerSetHelper::processReadForExtracting(pseq, uniqueKmersPerSet,
																												extractingPars, hasher, finalSeqOut,
																												currentCounts);
				}
				{
					std::lock_guard<std::mutex> lockGuard(mut);
					finalCounts.addOtherCounts(currentCounts);
				}
			};
			njh::concurrent::runVoidFunctionThreaded(readInComp, numThreads);
		}
		if (!singlesOption.firstName_.empty()) {
			SeqInput reader(singlesOption);
			reader.openIn();
			std::function<void()> readInComp = [&reader, &uniqueKmersPerSet,
							&mut, &extractingPars, &finalSeqOut, &finalCounts
//							, &testReads
			]() {
				SimpleKmerHash hasher;
				seqInfo seq;
				UniqueKmerSetHelper::ProcessReadForExtractingCounts currentCounts;
				while (reader.readNextReadLock(seq)) {
					UniqueKmerSetHelper::processReadForExtracting(seq, uniqueKmersPerSet,
																												extractingPars, hasher, finalSeqOut,
																												currentCounts);
				}
				{
					std::lock_guard<std::mutex> lockGuard(mut);
					finalCounts.addOtherCounts(currentCounts);
				}
			};
			njh::concurrent::runVoidFunctionThreaded(readInComp, numThreads);
		}
		finalCounts.writeOutCounts(outCounts,extractingPars,uniqueKmersPerSet,"final");
		if(setUp.pars_.verbose_){
			std::cout << "iterNumber: " << "final" << std::endl;
			std::cout << "Total Reads: " << finalCounts.getTotalCounts() << std::endl;
			std::cout << "Total Undetermined Reads: " << finalCounts.genTotalUndeterminedCount() << " (" << finalCounts.genTotalUndeterminedCount() * 100 /static_cast<double>(finalCounts.getTotalCounts())<< "%)" << std::endl;
			std::cout << "Total Determined Reads: " << finalCounts.genTotalDeterminedCount() << " (" << finalCounts.genTotalDeterminedCount() * 100 /static_cast<double>(finalCounts.getTotalCounts())<< "%)"<< std::endl;
			std::cout << "Total Less Than Small Len Cut Off Count: " << finalCounts.smallLenCutOffCount << " (" << finalCounts.smallLenCutOffCount * 100 /static_cast<double>(finalCounts.getTotalCounts())<< "%)"<< std::endl;
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

