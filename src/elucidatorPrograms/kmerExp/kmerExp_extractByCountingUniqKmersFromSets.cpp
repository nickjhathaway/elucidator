//
// Created by Nicholas Hathaway on 5/26/23.
//

#include "kmerExp.hpp"
#include <njhseq/IO/SeqIO/MultiSeqIO.hpp>
#include <njhseq/BamToolsUtils/bamExtractUtils.hpp>
#include "elucidator/helpers/UniqueKmerSetHelper.hpp"

namespace njhseq {

int kmerExpRunner::countingUniqKmersFromSetsPerRead(const njh::progutils::CmdArgs & inputCommands){
	uint32_t numThreads = 1;
	bfs::path countTable = "";
	UniqueKmerSetHelper::CompareReadToSetPars compPars;
	compPars.kmerLengthForEntropyCalc_ = 2;
	compPars.entropyFilter_ = 1.20;
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
	setUp.setOption(compPars.kmerLengthForEntropyCalc_, "--kmerLengthForEntropyCalc", "kmer Length For Entropy Calc");
	setUp.setOption(compPars.entropyFilter_, "--entropyFilter", "entropy Filter");
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
	if (!singlesOption.firstName_.empty()) {
		SeqInput reader(singlesOption);
		reader.openIn();
		readInComp = [&reader, &uniqueKmersPerSet, &compPars, &out, &outMut]() {
			SimpleKmerHash hasher;
			seqInfo seq;
			while (reader.readNextReadLock(seq)) {
				auto compRes = UniqueKmerSetHelper::compareReadToSetRes(seq, uniqueKmersPerSet, compPars, hasher);
				{
					std::lock_guard<std::mutex> lockGuard(outMut);
					compRes.writeOutput(out, seq, uniqueKmersPerSet, compPars);
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

	uint32_t iterationsToFilterDissimilar = std::numeric_limits<uint32_t>::max();
	UniqueKmerSetHelper::CompareReadToSetPars dissimilarFilterPars;
	dissimilarFilterPars.hardCountOff = 0;
	dissimilarFilterPars.fracCutOff = 0;

	UniqueKmerSetHelper::ProcessReadForExtractingPars extractingPars;
	extractingPars.compPars.hardCountOff = 20;
	extractingPars.compPars.fracCutOff = 0.12;
	extractingPars.compPars.kmerLengthForEntropyCalc_ = 2;
	extractingPars.compPars.entropyFilter_ = 1.20;

	extractingPars.compPars.initialExcludeHardCountOff = 60;
	extractingPars.compPars.initialExcludeFracCutOff = 0.25;
	uint32_t iterationToLowerExcludeCutOff = 6;

	bool finalExtractionPairsSeparate = false;



	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();

	//read in options
	bfs::path bamFnp;
	bfs::path bamExtractBedFnp;
	double bamExtractBedPercInRegion = 0.50;
	setUp.setOption(bamExtractBedPercInRegion, "--bamExtractBedPercInRegion", "bam Extract Bed Perc In Region", false);
	bool bamFnpSet = setUp.setOption(bamFnp, "--bamFnp", "bam Fnp", false);
	setUp.setOption(bamExtractBedFnp, "--bamExtractBedFnp", "bam Extract Bed Fnp regions", false);

	setUp.pars_.ioOptions_.revComplMate_ = true;
	bool pairedEndSet = setUp.processReadInNames(njhseq::seqSetUp::pairedReadInFormatsAvailable_, false);
	auto singlesOption = SeqIOOptions();
	bool singlesEndSet = setUp.processJustReadInNames(singlesOption, njhseq::seqSetUp::singleInFormatsAvailable_, !pairedEndSet && !bamFnpSet);
	singlesOption.lowerCaseBases_ = setUp.pars_.ioOptions_.lowerCaseBases_;
	singlesOption.removeGaps_ = setUp.pars_.ioOptions_.removeGaps_;
	singlesOption.includeWhiteSpaceInName_ = setUp.pars_.ioOptions_.includeWhiteSpaceInName_;

	if((singlesEndSet || pairedEndSet) && bamFnpSet){
		setUp.failed_ = true;
		setUp.addWarning(njh::pasteAsStr("can't set both bam file and fastq inputs"));
	}

	setUp.setOption(countTable, "--kmerTable", "countTable, 1)set,2)kmer", true);
	{
		auto possibleNonUniqueTable = njh::files::prependFileBasename(countTable, "nonUniqueKmers_");
		if(exists(possibleNonUniqueTable)){
			nonUniqueKmerTable = possibleNonUniqueTable;
		}
	}
	setUp.setOption(nonUniqueKmerTable, "--nonUniqueKmerTable", "non-unique Kmer Table, 1)set,2)kmer");
	setUp.setOption(extractingPars.compPars.sampleName, "--sampleName", "Name to add to output file", true);
	bool doNotIncludeRevComp = false;
	setUp.setOption(doNotIncludeRevComp, "--doNotIncludeRevComp", "do not include Rev Comp of the input seqs");
	extractingPars.compPars.includeRevComp = !doNotIncludeRevComp;
	//setUp.setOption(extractingPars.compPars.includeRevComp, "--includeRevComp", "include Rev Comp of the input seqs");

//	bool pairsTogether = false;
//	setUp.setOption(pairsTogether, "--pairsTogether", "extract the pairs together");
//	extractingPars.compPars.pairsSeparate = !pairsTogether;
	setUp.setOption(extractingPars.compPars.pairsSeparate, "--pairsSeparate", "extract the pairs separately");

	bool finalExtractionPairsTogether = false;
	setUp.setOption(finalExtractionPairsTogether, "--finalExtractionPairsTogether", "final extraction the pairs together");
	finalExtractionPairsSeparate = !finalExtractionPairsTogether;
	//setUp.setOption(finalExtractionPairsSeparate, "--finalExtractionPairsSeparate", "extract the pairs separately for final extraction step");



	//setUp.setOption(extractingPars.compPars.pairsSeparate, "--pairsSeparate", "extract the pairs separately");
	setUp.setOption(extractingPars.compPars.entropyFilter_, "--entropyFilter", "entropy Filter_");
	setUp.setOption(extractingPars.compPars.kmerLengthForEntropyCalc_, "--kmerLengthForEntropyCalc", "kmerLength For Entropy Calculation");


	setUp.setOption(dissimilarFilterPars.hardCountOff, "--hardCountOffForDissimilarFiltering", "hard Count Off For Dissimilar Filtering");
	setUp.setOption(dissimilarFilterPars.fracCutOff, "--fracCutOffForDissimilarFiltering", "frac Cut Off For Dissimilar Filtering");
	setUp.setOption(iterationsToFilterDissimilar, "--iterationsToFilterDissimilar", "iterations To Filter Dissimilar");

	setUp.setOption(extractingPars.compPars.initialExcludeFracCutOff, "--initialExcludeFracCutOff", "initial for exclusion sets per read Kmer Frac Cut Off");
	setUp.setOption(extractingPars.compPars.initialExcludeHardCountOff, "--initialExcludeHardCountOff", "initial for exclusion sets hard Count Off, do not count sets unless greater thant his number");
	setUp.setOption(iterationToLowerExcludeCutOff, "--iterationToLowerExcludeCutOff", "iterations To Lower Exclusion set Cut Off to the nomral cut offs");

	if(bfs::exists(countTable)){
		auto row = tokenizeString(njh::files::getFirstLine(countTable), "\t");
		extractingPars.compPars.klen = row[1].size();
	}
	setUp.setOption(extractingPars.compPars.fracCutOff, "--readKmerFracCutOff", "per read Kmer Frac Cut Off");
	setUp.setOption(extractingPars.compPars.hardCountOff, "--hardCountOff", "hard Count Off, do not count sets unless greater thant his number");
	extractingPars.smallLenCutOff = extractingPars.compPars.hardCountOff + extractingPars.compPars.klen + 11;
	setUp.setOption(extractingPars.smallLenCutOff, "--smallLenCutOff", "small Len Cut Off of input sequences, if less than this size will skip over");

	setUp.setOption(extractingPars.markReadsPerIteration, "--markReadsPerIteration", "mark Reads Per Iteration");
	setUp.setOption(extractingPars.writeOutFinalKmerSets, "--writeOutFinalKmerSets", "write Out Final Kmer Sets");

	setUp.setOption(excludesCountTable, "--excludesCountTable", "excludes Count Table, 1)set,2)kmer");
	setUp.setOption(extractingPars.compPars.excludeSetNames, "--excludeSetNames", "names of sets of unique kmers to ignore from the --kmerTable");
	//
	setUp.setOption(keepTemporaryFiles, "--keepTemporaryFiles", "Keep Temporary Files");
	setUp.setOption(extractingPars.addingInKmersCountCutOff, "--addingInKmersCountCutOff", "adding In Kmers Count Cut Off");


	setUp.setOption(extractingPars.doReCheckExcludeSets, "--doReCheckExcludeSets", "do Re Check Exclude Sets on iterations, could lead to some reads being recruit away but will drastically speed up run time");

	bool doNotFinalRecheckOfExcludeSeqs = false;
	setUp.setOption(doNotFinalRecheckOfExcludeSeqs, "--doNotFinalRecheckOfExcludeSeqs", "redo a unique kmer filter on the exclusion sets during the final re-extraction process");

	bool doFinalRecheckOfExcludeSeqs = !doNotFinalRecheckOfExcludeSeqs;
//	setUp.setOption(doFinalRecheckOfExcludeSeqs, "--doFinalRecheckOfExcludeSeqs", "Normal the exclusion kmer set is left alone, if this flag is set, it will be filter during the final re-extraction process");

	setUp.setOption(doNotDoFinalExtract, "--doNotDoFinalExtract", "do Not Do Final Extract");
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
	OutOptions bamExtractOut(njh::files::make_path(setUp.pars_.directoryName_, "bamExtract"));
	bamExtractOut.append_ = true;
	if (bamFnpSet) {
		watch.startNewLap(njh::pasteAsStr("bam extract files"));
		njh::files::checkExistenceThrow(bamFnp, __PRETTY_FUNCTION__);
		BamExtractor bExtractor(setUp.pars_.verbose_);
		if (bamExtractBedFnp.empty()) {
			bExtractor.writeExtractReadsFromBam(bamFnp, bamExtractOut);
		} else {
			auto regions = gatherRegions(bamExtractBedFnp.string(), "", setUp.pars_.verbose_);
			for (const auto &region: regions) {
				bExtractor.writeExtractReadsFromBamRegion(bamFnp,
																									region, bamExtractBedPercInRegion, bamExtractOut);
			}
		}

		bfs::path bamExtractR1Fnp = bamExtractOut.outFilename_.string() + "_R1.fastq.gz";
		bfs::path bamExtractR2Fnp = bamExtractOut.outFilename_.string() + "_R2.fastq.gz";
		bfs::path bamExtractFnp = bamExtractOut.outFilename_.string() + ".fastq.gz";
		if (exists(bamExtractR1Fnp)) {
			setUp.pars_.ioOptions_ = SeqIOOptions::genPairedInGz(bamExtractR1Fnp, bamExtractR2Fnp);
//			setUp.pars_.ioOptions_.inFormat_ = SeqIOOptions::inFormats::FASTQPAIREDGZ;
//			setUp.pars_.ioOptions_.outFormat_ = SeqIOOptions::outFormats::FASTQPAIREDGZ;
//			setUp.pars_.ioOptions_.firstName_ = bamExtractR1Fnp;
//			setUp.pars_.ioOptions_.firstName_ = bamExtractR2Fnp;
		}
		if (exists(bamExtractFnp)) {
			singlesOption = SeqIOOptions::genFastqInGz(bamExtractFnp);
//			singlesOption.inFormat_ = SeqIOOptions::inFormats::FASTQGZ;
//			singlesOption.outFormat_ = SeqIOOptions::outFormats::FASTQGZ;
//			singlesOption.firstName_ = bamExtractFnp;
		}
		if(setUp.pars_.verbose_){
			std::cout << watch.getLapName() << "\t" << watch.timeLapFormatted() <<std::endl;
		}
	}

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
		auto setNames = njh::getSetOfMapKeys(uniqueKmersPerSet);
		for(const auto & setName : setNames){
			std::cout << setName << "\t" << uniqueKmersPerSet[setName].size() << std::endl;
		}
		std::cout << "NON_UNIQUE" << "\t" << nonUniqueKmersPerSet.size() << std::endl;
	}
	MultiSeqIO initialSeqOut;
	watch.startNewLap("initial scan");

	OutOptions outCountsOpts(njh::files::make_path(setUp.pars_.directoryName_, "extractionCounts.tab.txt"));
	OutputStream outCounts(outCountsOpts);
	UniqueKmerSetHelper::ProcessReadForExtractingCounts::writeOutCountsHeader(outCounts, extractingPars);




	std::unordered_map<std::string, UniqueKmerSetHelper::FilePositons> positionsAfterLastIteration;
	{//initial

		//add outputs
		VecStr names = getVectorOfMapKeys(uniqueKmersPerSet);
		for(const auto & name : names){
			auto seqOutOpts = SeqIOOptions::genPairedOut(njh::files::make_path(setUp.pars_.directoryName_, name));
			initialSeqOut.addReader(njh::pasteAsStr(name, "-paired"), seqOutOpts);
		}
		initialSeqOut.addReader("undetermined-paired", SeqIOOptions::genPairedOut(njh::files::make_path(setUp.pars_.directoryName_, "undetermined")));
		for(const auto & name : names){
			auto seqOutOpts = SeqIOOptions::genFastqOut(njh::files::make_path(setUp.pars_.directoryName_, name));
			initialSeqOut.addReader(njh::pasteAsStr(name, "-single"), seqOutOpts);
		}
		initialSeqOut.addReader("undetermined-single", SeqIOOptions::genFastqOut(njh::files::make_path(setUp.pars_.directoryName_, "undetermined")));

		//to make below count correctly, make sure mate is automatically reversed complemented
		setUp.pars_.ioOptions_.revComplMate_ = true;

		if (!setUp.pars_.ioOptions_.firstName_.empty()) {
			SeqInput reader(setUp.pars_.ioOptions_);
			reader.openIn();
			std::function<void()> readInComp;

			if(extractingPars.compPars.pairsSeparate){
				readInComp = [&reader, &uniqueKmersPerSet,
								&mut, &extractingPars, &initialSeqOut, &initialCounts
				]() {
					SimpleKmerHash hasher;
					PairedRead pseq;
					UniqueKmerSetHelper::ProcessReadForExtractingCounts currentCounts;
					while (reader.readNextReadLock(pseq)) {
						UniqueKmerSetHelper::processReadForExtractingPairsSeparate(pseq, uniqueKmersPerSet,
																																			 extractingPars, hasher, initialSeqOut,
																																			 currentCounts,
																																			 "initial");
					}
					{
						std::lock_guard<std::mutex> lockGuard(mut);
						initialCounts.addOtherCounts(currentCounts);
					}
				};
			}else{
				readInComp = [&reader, &uniqueKmersPerSet,
								&mut, &extractingPars, &initialSeqOut, &initialCounts
				]() {
					SimpleKmerHash hasher;
					PairedRead pseq;
					UniqueKmerSetHelper::ProcessReadForExtractingCounts currentCounts;
					while (reader.readNextReadLock(pseq)) {
						UniqueKmerSetHelper::processReadForExtractingPairsTogether(pseq, uniqueKmersPerSet,
																													extractingPars, hasher, initialSeqOut,
																													currentCounts,
																													"initial");
					}
					{
						std::lock_guard<std::mutex> lockGuard(mut);
						initialCounts.addOtherCounts(currentCounts);
					}
				};
			}
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
																												currentCounts, "initial");
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
	initialSeqOut.closeOutForReopeningAll();

	{
		auto kmerSetNames = njh::getSetOfMapKeys(uniqueKmersPerSet);
		for(const auto & kmerSetName : kmerSetNames){
			if(njh::in(kmerSetName, extractingPars.compPars.excludeSetNames)){
				continue;
			}
			UniqueKmerSetHelper::FilePositons positions;
			auto r1Fnp = njh::files::make_path(setUp.pars_.directoryName_, kmerSetName + "_R1.fastq");
			auto r2Fnp =  njh::files::make_path(setUp.pars_.directoryName_, kmerSetName + "_R2.fastq");
			auto singleFnp = njh::files::make_path(setUp.pars_.directoryName_, kmerSetName + ".fastq");
			positions.r1FnpEnd = njh::files::getFilePositionsEndOfFile(r1Fnp);
			positions.r2FnpEnd = njh::files::getFilePositionsEndOfFile(r2Fnp);
			positions.singleFnpEnd = njh::files::getFilePositionsEndOfFile(singleFnp);
			positionsAfterLastIteration[kmerSetName] = positions;
		}
	}
	//
	if(setUp.pars_.verbose_) {
		std::cout << "iterNumber: " << "initial" << std::endl;
		std::cout << "Total Initial Reads: " << initialCounts.getTotalCounts() << std::endl;
		std::cout << "Total Initial Undetermined Reads: " << initialCounts.genTotalUndeterminedCount() << " ("
							<< initialCounts.genTotalUndeterminedCount() * 100 / static_cast<double>(initialCounts.getTotalCounts())
							<< "%)" << std::endl;
		std::cout << "Total Initial Determined Reads: " << initialCounts.genTotalDeterminedCount() << " ("
							<< initialCounts.genTotalDeterminedCount() * 100 / static_cast<double>(initialCounts.getTotalCounts())
							<< "%)" << std::endl;
		std::cout << "Total Initial Less Than Small Len Cut Off Count: " << initialCounts.smallLenCutOffCount << " ("
							<< initialCounts.smallLenCutOffCount * 100 / static_cast<double>(initialCounts.getTotalCounts()) << "%)"
							<< std::endl;
		std::cout << watch.getLapName() << "\t" << watch.timeLapFormatted() << std::endl;
	}
	//read in the extract reads and add to the unique kmer sets
	bool iterate = true;
	uint64_t totalUndeterminedExtracted = initialCounts.genTotalUndeterminedCount();
	uint32_t iterNumber = 0;
	if(0 == totalUndeterminedExtracted){
		iterate = false;
	}




	auto original_initialCounts = initialCounts;
	initialCounts.readCountsPerSet[false]["undetermined"] = 0;
	initialCounts.readCountsPerSet[true]["undetermined"] = 0;
	uint64_t totalDeterminedReads = initialCounts.genTotalDeterminedCount();
//	std::cout << __FILE__ << " " << __LINE__ << std::endl;
	while(iterate && iterNumber < maxIterations){
//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		watch.startNewLap(njh::pasteAsStr(iterNumber, "- read in new kmers"));
//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		initialSeqOut.closeOutForReopeningAll();
		auto rawUniqKmerSets = njh::getVecOfMapKeys(uniqueKmersPerSet);
		removeElement(rawUniqKmerSets, std::string("undetermined"));
		removeElements(rawUniqKmerSets, VecStr(extractingPars.compPars.excludeSetNames.begin(), extractingPars.compPars.excludeSetNames.end()));
		std::unordered_map<std::string, std::unordered_map<uint64_t, uint32_t>> rawKmersPerInput;
		if (iterNumber > 0) {
			rawKmersPerInput = UniqueKmerSetHelper::readInNewKmersFromExtractedReads(
							setUp.pars_.directoryName_, rawUniqKmerSets, extractingPars, positionsAfterLastIteration);
		} else {
			rawKmersPerInput = UniqueKmerSetHelper::readInNewKmersFromExtractedReads(
							setUp.pars_.directoryName_, rawUniqKmerSets, extractingPars);
		}
		if (setUp.pars_.verbose_ && setUp.pars_.debug_) {
			std::cout << "new kmer counts for: " << rawKmersPerInput.size() << " sets" << std::endl;
			auto setNames = njh::getSetOfMapKeys(rawKmersPerInput);
			for (const auto &setName: setNames) {
				std::cout << setName << "\t" << rawKmersPerInput[setName].size() << std::endl;
			}
			std::cout << "NON_UNIQUE" << "\t" << nonUniqueKmersPerSet.size() << std::endl << std::endl;
		}

//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		watch.startNewLap(njh::pasteAsStr(iterNumber, "- check new kmers against current set"));
		std::unordered_map<std::string, std::unordered_set<uint64_t>> outputUniqueKmersPerSet = UniqueKmerSetHelper::filterReExtractedKmersForNonUnique(rawKmersPerInput,
																																																																										extractingPars,
																																																																										uniqueKmersPerSet,
																																																																										nonUniqueKmersPerSet);
//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		if(!extractingPars.doReCheckExcludeSets){
			for(const auto & excludeName : extractingPars.compPars.excludeSetNames){
				outputUniqueKmersPerSet[excludeName] = std::move(uniqueKmersPerSet[excludeName]);
			}
		}

//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		watch.startNewLap(njh::pasteAsStr(iterNumber, "- copy new set"));
		uniqueKmersPerSet = std::move(outputUniqueKmersPerSet);
		watch.startNewLap(njh::pasteAsStr(iterNumber, "- scan again"));
//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		if (setUp.pars_.verbose_ && setUp.pars_.debug_) {
			auto setNames = njh::getSetOfMapKeys(uniqueKmersPerSet);
			for (const auto &setName: setNames) {
				std::cout << setName << "\t" << uniqueKmersPerSet[setName].size() << std::endl;
			}
			std::cout << "NON_UNIQUE" << "\t" << nonUniqueKmersPerSet.size() << std::endl;
		}

		if(iterNumber >= iterationToLowerExcludeCutOff){
			extractingPars.compPars.initialExcludeHardCountOff = extractingPars.compPars.hardCountOff;
			extractingPars.compPars.initialExcludeFracCutOff = extractingPars.compPars.fracCutOff;
		}

		bfs::path orig_undeterminedInput = njh::files::make_path(setUp.pars_.directoryName_, "undetermined.fastq");
		bfs::path orig_undeterminedInputR1 = njh::files::make_path(setUp.pars_.directoryName_, "undetermined_R1.fastq");
		bfs::path orig_undeterminedInputR2 = njh::files::make_path(setUp.pars_.directoryName_, "undetermined_R2.fastq");

		bfs::path undeterminedInput = njh::files::make_path(setUp.pars_.directoryName_, njh::pasteAsStr(iterNumber, "_undetermined.fastq"));
		bfs::path undeterminedInputR1 = njh::files::make_path(setUp.pars_.directoryName_, njh::pasteAsStr(iterNumber, "_undetermined_R1.fastq"));
		bfs::path undeterminedInputR2 = njh::files::make_path(setUp.pars_.directoryName_, njh::pasteAsStr(iterNumber, "_undetermined_R2.fastq"));

		SeqIOOptions undeterminedPairedInOpts;
		SeqIOOptions undeterminedSingleInOpts;
		bool pairedInUndeterminedExists = bfs::exists(orig_undeterminedInputR1);
		bool singleInUndeterminedExists = bfs::exists(orig_undeterminedInput);

		if (iterNumber == iterationsToFilterDissimilar) {
			//add outputs
			auto extractingParsForFiltering = extractingPars;
			extractingParsForFiltering.compPars.hardCountOff = dissimilarFilterPars.hardCountOff;
			extractingParsForFiltering.compPars.fracCutOff = dissimilarFilterPars.fracCutOff;
			extractingParsForFiltering.compPars.initialExcludeHardCountOff = dissimilarFilterPars.hardCountOff;
			extractingParsForFiltering.compPars.initialExcludeFracCutOff = dissimilarFilterPars.fracCutOff;

			extractingParsForFiltering.compPars.excludeSetNames = extractingPars.compPars.excludeSetNames;
			MultiSeqIO filteringSeqOut;
			std::string keep = "undetermined";
			std::string filter = "dissimilar";
			filteringSeqOut.addReader("undetermined-paired", SeqIOOptions::genPairedOut(njh::files::make_path(setUp.pars_.directoryName_, njh::pasteAsStr(iterNumber, "_undetermined"))));
			filteringSeqOut.addReader("undetermined-single", SeqIOOptions::genFastqOut(njh::files::make_path(setUp.pars_.directoryName_, njh::pasteAsStr(iterNumber, "_undetermined"))));
			filteringSeqOut.addReader("dissimilar-paired", SeqIOOptions::genPairedOut(njh::files::make_path(setUp.pars_.directoryName_, njh::pasteAsStr(iterNumber, "_dissimilar"))));
			filteringSeqOut.addReader("dissimilar-single", SeqIOOptions::genFastqOut(njh::files::make_path(setUp.pars_.directoryName_, njh::pasteAsStr(iterNumber, "_dissimilar"))));
			if (pairedInUndeterminedExists) {
				SeqIOOptions undeterminedPairedInOptsForFiltering;
				undeterminedPairedInOptsForFiltering.firstName_ = orig_undeterminedInputR1;
				undeterminedPairedInOptsForFiltering.secondName_ = orig_undeterminedInputR2;
				undeterminedPairedInOptsForFiltering.inFormat_ = SeqIOOptions::inFormats::FASTQPAIRED;
				undeterminedPairedInOptsForFiltering.revComplMate_ = true;


				//re-extract paired in
				SeqInput reader(undeterminedPairedInOptsForFiltering);
				reader.openIn();
				std::function<void()> readInComp;

				if(extractingPars.compPars.pairsSeparate){
					readInComp = [&reader, &uniqueKmersPerSet,
									&mut, &extractingParsForFiltering, &filteringSeqOut, &initialCounts,
									&iterNumber,&keep,&filter
					]() {
						SimpleKmerHash hasher;
						PairedRead pseq;
						UniqueKmerSetHelper::ProcessReadForExtractingCounts currentCounts;
						while (reader.readNextReadLock(pseq)) {
							UniqueKmerSetHelper::processReadForFilteringPairsSeparate(pseq, uniqueKmersPerSet,
																																				extractingParsForFiltering, hasher, filteringSeqOut,
																																				 currentCounts,
																																				 njh::pasteAsStr(iterNumber),
																																				 filter,
																																				 keep);
						}
						{
							std::lock_guard<std::mutex> lockGuard(mut);
							initialCounts.addOtherCounts(currentCounts);
						}
					};
				} else {
					readInComp = [&reader, &uniqueKmersPerSet,
									&mut, &extractingParsForFiltering, &filteringSeqOut, &initialCounts,
									&iterNumber,&keep,&filter
					]() {
						SimpleKmerHash hasher;
						PairedRead pseq;
						UniqueKmerSetHelper::ProcessReadForExtractingCounts currentCounts;
						while (reader.readNextReadLock(pseq)) {
							UniqueKmerSetHelper::processReadForFilteringPairsTogether(pseq, uniqueKmersPerSet,
																																				extractingParsForFiltering, hasher, filteringSeqOut,
																																				 currentCounts,
																																				 njh::pasteAsStr(iterNumber),
																																				filter,
																																				keep);
						}
						{
							std::lock_guard<std::mutex> lockGuard(mut);
							initialCounts.addOtherCounts(currentCounts);
						}
					};
				}
				njh::concurrent::runVoidFunctionThreaded(readInComp, numThreads);
			}
			if (singleInUndeterminedExists) {
				//re-extract single in
				SeqIOOptions undeterminedSingleInOptsForFiltering;
				undeterminedSingleInOptsForFiltering.firstName_ = orig_undeterminedInput;
				undeterminedSingleInOptsForFiltering.inFormat_ = SeqIOOptions::inFormats::FASTQ;
				SeqInput reader(undeterminedSingleInOptsForFiltering);
				reader.openIn();
				std::function<void()> readInComp = [&reader, &uniqueKmersPerSet,
								&mut, &extractingParsForFiltering, &filteringSeqOut, &initialCounts
								, &iterNumber,&keep,&filter]() {
					SimpleKmerHash hasher;
					seqInfo seq;
					UniqueKmerSetHelper::ProcessReadForExtractingCounts currentCounts;
					while (reader.readNextReadLock(seq)) {
						UniqueKmerSetHelper::processReadForFiltering(seq, uniqueKmersPerSet,
																												 extractingParsForFiltering, hasher, filteringSeqOut,
																												 currentCounts,
																												 njh::pasteAsStr(iterNumber),
																												 filter,
																												 keep);
					}
					{
						std::lock_guard<std::mutex> lockGuard(mut);
						initialCounts.addOtherCounts(currentCounts);
					}
				};
				njh::concurrent::runVoidFunctionThreaded(readInComp, numThreads);
			}
			if (pairedInUndeterminedExists) {
				bfs::remove(orig_undeterminedInputR1);
				bfs::remove(orig_undeterminedInputR2);
				undeterminedPairedInOpts.firstName_ = undeterminedInputR1;
				undeterminedPairedInOpts.secondName_ = undeterminedInputR2;
				undeterminedPairedInOpts.inFormat_ = SeqIOOptions::inFormats::FASTQPAIRED;
				undeterminedPairedInOpts.revComplMate_ = true;
			}
			if (singleInUndeterminedExists) {
				bfs::remove(orig_undeterminedInput);
				undeterminedSingleInOpts.firstName_ = undeterminedInput;
				undeterminedSingleInOpts.inFormat_ = SeqIOOptions::inFormats::FASTQ;
			}
		} else {
			if (pairedInUndeterminedExists) {
				bfs::rename(orig_undeterminedInputR1, undeterminedInputR1);
				bfs::rename(orig_undeterminedInputR2, undeterminedInputR2);
				undeterminedPairedInOpts.firstName_ = undeterminedInputR1;
				undeterminedPairedInOpts.secondName_ = undeterminedInputR2;
				undeterminedPairedInOpts.inFormat_ = SeqIOOptions::inFormats::FASTQPAIRED;
				undeterminedPairedInOpts.revComplMate_ = true;
				//to make below count correctly, make sure mate is automatically reversed complemented
			}
			//std::cout << __FILE__ << " " << __LINE__ << std::endl;
			if (singleInUndeterminedExists) {
				bfs::rename(orig_undeterminedInput, undeterminedInput);
				undeterminedSingleInOpts.firstName_ = undeterminedInput;
				undeterminedSingleInOpts.inFormat_ = SeqIOOptions::inFormats::FASTQ;
			}
		}

		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		{
			//determine positions before extracting more
			auto kmerSetNames = njh::getSetOfMapKeys(uniqueKmersPerSet);
			for(const auto & kmerSetName : kmerSetNames){
				if(njh::in(kmerSetName, extractingPars.compPars.excludeSetNames)){
					continue;
				}
				UniqueKmerSetHelper::FilePositons positions;
				auto r1Fnp = njh::files::make_path(setUp.pars_.directoryName_, kmerSetName + "_R1.fastq");
				auto r2Fnp =  njh::files::make_path(setUp.pars_.directoryName_, kmerSetName + "_R2.fastq");
				auto singleFnp = njh::files::make_path(setUp.pars_.directoryName_, kmerSetName + ".fastq");
				positions.r1FnpEnd = njh::files::getFilePositionsEndOfFile(r1Fnp);
				positions.r2FnpEnd = njh::files::getFilePositionsEndOfFile(r2Fnp);
				positions.singleFnpEnd = njh::files::getFilePositionsEndOfFile(singleFnp);
				positionsAfterLastIteration[kmerSetName] = positions;
			}
		}
		if (pairedInUndeterminedExists) {
			//std::cout << __FILE__ << " " << __LINE__ << std::endl;
			//re-extract paired in
			SeqInput reader(undeterminedPairedInOpts);
			reader.openIn();
			std::function<void()> readInComp;

			if(extractingPars.compPars.pairsSeparate){
				readInComp = [&reader, &uniqueKmersPerSet,
								&mut, &extractingPars, &initialSeqOut, &initialCounts,
								&iterNumber
				]() {
					SimpleKmerHash hasher;
					PairedRead pseq;
					UniqueKmerSetHelper::ProcessReadForExtractingCounts currentCounts;
					while (reader.readNextReadLock(pseq)) {
						UniqueKmerSetHelper::processReadForExtractingPairsSeparate(pseq, uniqueKmersPerSet,
																																			 extractingPars, hasher, initialSeqOut,
																																			 currentCounts,
																																			 njh::pasteAsStr(iterNumber));
					}
					{
						std::lock_guard<std::mutex> lockGuard(mut);
						initialCounts.addOtherCounts(currentCounts);
					}
				};
			}else{
				readInComp = [&reader, &uniqueKmersPerSet,
								&mut, &extractingPars, &initialSeqOut, &initialCounts,
								&iterNumber
				]() {
					SimpleKmerHash hasher;
					PairedRead pseq;
					UniqueKmerSetHelper::ProcessReadForExtractingCounts currentCounts;
					while (reader.readNextReadLock(pseq)) {
						UniqueKmerSetHelper::processReadForExtractingPairsTogether(pseq, uniqueKmersPerSet,
																																			 extractingPars, hasher, initialSeqOut,
																																			 currentCounts,
																																			 njh::pasteAsStr(iterNumber));
					}
					{
						std::lock_guard<std::mutex> lockGuard(mut);
						initialCounts.addOtherCounts(currentCounts);
					}
				};
			}
			njh::concurrent::runVoidFunctionThreaded(readInComp, numThreads);
		}
		if(singleInUndeterminedExists) {
			//re-extract single in
			SeqInput reader(undeterminedSingleInOpts);
			reader.openIn();
			std::function<void()> readInComp = [&reader, &uniqueKmersPerSet,
							&mut, &extractingPars, &initialSeqOut, &initialCounts
//							, &testReads
			, &iterNumber]() {
				SimpleKmerHash hasher;
				seqInfo seq;
				UniqueKmerSetHelper::ProcessReadForExtractingCounts currentCounts;
				while (reader.readNextReadLock(seq)) {
					UniqueKmerSetHelper::processReadForExtracting(seq, uniqueKmersPerSet,
																												extractingPars, hasher, initialSeqOut,
																												currentCounts,
																												njh::pasteAsStr(iterNumber));
				}

				{
					std::lock_guard<std::mutex> lockGuard(mut);
					initialCounts.addOtherCounts(currentCounts);
				}
			};

			njh::concurrent::runVoidFunctionThreaded(readInComp, numThreads);
		}
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		initialCounts.writeOutCounts(outCounts,extractingPars,uniqueKmersPerSet,njh::pasteAsStr(iterNumber));
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		uint64_t currentTotalUndeterminedReads = initialCounts.genTotalUndeterminedCount();
		uint64_t currentTotalReads = initialCounts.getTotalCounts();
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		initialCounts.readCountsPerSet[false]["undetermined"] = 0;
		initialCounts.readCountsPerSet[true]["undetermined"] = 0;
		uint64_t currentTotalDeterminedReads = initialCounts.genTotalDeterminedCount();
		if(setUp.pars_.verbose_){
			std::cout << "iterNumber: " << iterNumber << std::endl;
			std::cout << "Determined: " << currentTotalDeterminedReads - totalDeterminedReads << " more reads" <<std::endl;
			std::cout << "Current counts: " << std::endl;
			std::cout << "\tTotal Reads: " << currentTotalReads << std::endl;
			std::cout << "\tTotal Undetermined Reads: " << currentTotalUndeterminedReads << " (" << currentTotalUndeterminedReads * 100 /static_cast<double>(currentTotalReads)<< "%)" << std::endl;
			std::cout << "\tTotal Determined Reads: " << currentTotalDeterminedReads << " (" << currentTotalDeterminedReads * 100 /static_cast<double>(currentTotalReads)<< "%)"<< std::endl;
			std::cout << "\tTotal Less Than Small Len Cut Off Count: " << initialCounts.smallLenCutOffCount << " (" << initialCounts.smallLenCutOffCount * 100 /static_cast<double>(currentTotalReads)<< "%)"<< std::endl;
			std::cout << "\tDissimilar Filter: " << initialCounts.filteredDissimilarCount << " (" << initialCounts.filteredDissimilarCount * 100 /static_cast<double>(currentTotalReads)<< "%)"<< std::endl;
			std::cout << std::endl;
		}
		if(currentTotalUndeterminedReads == 0 || totalDeterminedReads == currentTotalDeterminedReads){
			//if there's no more undetermined reads or if the number determined this round was zero then stop iterating
			iterate = false;
		}
		totalDeterminedReads = currentTotalDeterminedReads;
		if (!keepTemporaryFiles) {
			if (bfs::exists(undeterminedPairedInOpts.firstName_)) {
				bfs::remove(undeterminedPairedInOpts.firstName_);
				bfs::remove(undeterminedPairedInOpts.secondName_);
			}
			if (bfs::exists(undeterminedSingleInOpts.firstName_)) {
				bfs::remove(undeterminedSingleInOpts.firstName_);
			}
		}
		++iterNumber;
	}
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

		extractingPars.compPars.pairsSeparate = finalExtractionPairsSeparate;
		if (!setUp.pars_.ioOptions_.firstName_.empty()) {
			SeqInput reader(setUp.pars_.ioOptions_);
			reader.openIn();

			std::function<void()> readInComp;

			if(extractingPars.compPars.pairsSeparate){
				readInComp = [&reader, &uniqueKmersPerSet,
								&mut, &extractingPars, &finalSeqOut, &finalCounts
				]() {
					SimpleKmerHash hasher;
					PairedRead pseq;
					UniqueKmerSetHelper::ProcessReadForExtractingCounts currentCounts;
					while (reader.readNextReadLock(pseq)) {
						UniqueKmerSetHelper::processReadForExtractingPairsSeparate(pseq, uniqueKmersPerSet,
																																			 extractingPars, hasher, finalSeqOut,
																																			 currentCounts,
																																			 "final");
					}
					{
						std::lock_guard<std::mutex> lockGuard(mut);
						finalCounts.addOtherCounts(currentCounts);
					}
				};
			}else{
				readInComp = [&reader, &uniqueKmersPerSet,
								&mut, &extractingPars, &finalSeqOut, &finalCounts
				]() {
					SimpleKmerHash hasher;
					PairedRead pseq;
					UniqueKmerSetHelper::ProcessReadForExtractingCounts currentCounts;
					while (reader.readNextReadLock(pseq)) {
						UniqueKmerSetHelper::processReadForExtractingPairsTogether(pseq, uniqueKmersPerSet,
																																			 extractingPars, hasher, finalSeqOut,
																																			 currentCounts,
																																			 "final");
					}
					{
						std::lock_guard<std::mutex> lockGuard(mut);
						finalCounts.addOtherCounts(currentCounts);
					}
				};
			}
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
																												currentCounts,
																												"final");
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
			if (setUp.pars_.debug_) {
				std::cout << watch.getLapName() << "\t" << watch.timeLapFormatted() << std::endl;
				std::cout << "klen: " << extractingPars.compPars.klen << std::endl;
				auto setNames = njh::getSetOfMapKeys(uniqueKmersPerSet);
				for (const auto &setName: setNames) {
					std::cout << setName << "\t" << uniqueKmersPerSet[setName].size() << std::endl;
				}
				std::cout << "NON_UNIQUE" << "\t" << nonUniqueKmersPerSet.size() << std::endl;
			}
		}
		if(!keepTemporaryFiles){
			//remove all fast[aq] files and fast[aq].gz files within directory
			auto readFiles = njh::files::listAllFiles(setUp.pars_.directoryName_,false,{std::regex{R"((.*)(\.(fast[aq])(\.gz)?)$)"}});
			for(const auto & f : readFiles){
				bfs::remove(f.first);
			}
		} else {
			//zip up any temp files if keeping
			auto readFiles = njh::files::listAllFiles(setUp.pars_.directoryName_,false,{std::regex{R"((.*)(\.(fast[aq]))$)"}});
			for(const auto & f : readFiles){
				njh::gzZipFile(njh::IoOptions{njh::InOptions(f.first), njh::OutOptions(bfs::path(f.first.string() + ".gz"))});
				bfs::remove(f.first);
			}
		}
	} else {
		auto readFiles = njh::files::listAllFiles(setUp.pars_.directoryName_,false,{std::regex{R"((.*)(\.(fast[aq]))$)"}});
		for(const auto & f : readFiles){
			njh::gzZipFile(njh::IoOptions{njh::InOptions(f.first), njh::OutOptions(bfs::path(f.first.string() + ".gz"))});
			bfs::remove(f.first);
		}
	}

	if(extractingPars.writeOutFinalKmerSets){
		SimpleKmerHash hasher;
		watch.startNewLap(njh::pasteAsStr("write Out Final Kmer Sets"));
		OutputStream finalUniqueKmerSetsOut(njh::files::make_path(setUp.pars_.directoryName_, "finalUniqueKmerSets.tsv.gz"));
		for(const auto & uniqKSet : uniqueKmersPerSet){
			for(const auto & uniqK : uniqKSet.second){
				finalUniqueKmerSetsOut << uniqKSet.first << "\t" << hasher.reverseHash(uniqK)  << "\n";
			}
		}
		OutputStream nonUniqueKmerSets(njh::files::make_path(setUp.pars_.directoryName_, "nonUniqueKmers_finalUniqueKmerSets.tsv.gz"));
		for(const auto & nonUniqK : nonUniqueKmersPerSet){
			nonUniqueKmerSets << nonUniqueRegionName << "\t" << hasher.reverseHash(nonUniqK) << "\n";
		}
	}
	if(!keepTemporaryFiles){
		if (bamFnpSet) {
			watch.startNewLap(njh::pasteAsStr("remove bam extracted files"));

			bfs::path bamExtractR1Fnp = bamExtractOut.outFilename_.string() + "_R1.fastq.gz";
			bfs::path bamExtractR2Fnp = bamExtractOut.outFilename_.string() + "_R2.fastq.gz";
			bfs::path bamExtractFnp = bamExtractOut.outFilename_.string() + ".fastq.gz";
			if (exists(bamExtractR1Fnp)) {
				bfs::remove(bamExtractR1Fnp);
				bfs::remove(bamExtractR2Fnp);
			}
			if (exists(bamExtractFnp)) {
				bfs::remove(bamExtractFnp);
			}
		}
	}
	watch.logLapTimes(setUp.rLog_.runLogFile_, true, 6, true);
	return 0;
}


}  // namespace njhseq

