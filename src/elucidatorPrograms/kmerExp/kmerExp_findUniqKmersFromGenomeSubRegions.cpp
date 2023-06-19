//
// Created by Nicholas Hathaway on 5/6/23.
//

#include "kmerExp.hpp"
#include <njhseq/objects/kmer/KmerUtils.hpp>
#include <njhseq/objects/kmer.h>


#include "elucidator/objects/MiscUtility/GenomeSeqSearch.hpp"
#include "elucidator/helpers/UniqueKmerSetHelper.hpp"

#include <njhseq/objects/dataContainers/tables/TableReader.hpp>
#include <njhseq/objects/BioDataObject/reading.hpp>
#include <njhseq/objects/BioDataObject/BioRecordsUtils/BedUtility.hpp>
#include <njhseq/objects/counters/DNABaseCounter.hpp>

namespace njhseq {



int kmerExpRunner::uniqKmersSetToFasta(const njh::progutils::CmdArgs & inputCommands){

	bfs::path countTable;
	bfs::path nonUniqueKmerTable;

	OutOptions outOpts;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.description_ = "Add unique kmers to set of already determined unique kmers";
	setUp.setOption(countTable, "--kmerTable", "unique kmer sets, 1)set,2)kmer", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	njh::files::checkExistenceThrow(countTable, __PRETTY_FUNCTION__ );
	if(!nonUniqueKmerTable.empty()){
		njh::files::checkExistenceThrow(nonUniqueKmerTable, __PRETTY_FUNCTION__ );
	}

	std::unordered_map<std::string, std::unordered_set<uint64_t>> uniqueKmersPerSet = UniqueKmerSetHelper::readInUniqueKmerTablePerSet(
					countTable);
	SimpleKmerHash hasher;
	for (const auto &set: uniqueKmersPerSet) {
		bfs::path outFnp = set.first;
		if (!outOpts.outFilename_.empty()) {
			outFnp = njh::files::nameAppendBeforeExt(outOpts.outFilename_, set.first);
		}
		auto seqOutForSet = SeqIOOptions::genFastaOutGz(outFnp);
		seqOutForSet.out_.transferOverwriteOpts(outOpts);

		SeqOutput writer(seqOutForSet);
		writer.openOut();
		uint32_t kmerCount = 0;
		for(const auto & k : set.second){
			writer.write(seqInfo(estd::to_string(kmerCount), hasher.reverseHash(k)));
			++kmerCount;
		}
	}

	return 0;
}

int kmerExpRunner::reportOnUniqKmersSet(const njh::progutils::CmdArgs & inputCommands){

	bfs::path countTable;
	bfs::path nonUniqueKmerTable;

	OutOptions outOpts;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.description_ = "Add unique kmers to set of already determined unique kmers";
	setUp.setOption(nonUniqueKmerTable, "--nonUniqueKmerTable", "non-unique Kmer Table, 1)set,2)kmer");
	setUp.setOption(countTable, "--kmerTable,--countTable", "unique kmer sets, 1)set,2)kmer", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	njh::files::checkExistenceThrow(countTable, __PRETTY_FUNCTION__ );
	if(!nonUniqueKmerTable.empty()){
		njh::files::checkExistenceThrow(nonUniqueKmerTable, __PRETTY_FUNCTION__ );
	}
	OutputStream out(outOpts);
	uint32_t klen = UniqueKmerSetHelper::getKmerLenFromUniqueKmerTable(countTable);
	std::unordered_map<std::string, std::unordered_set<uint64_t>> uniqueKmersPerSet = UniqueKmerSetHelper::readInUniqueKmerTablePerSet(
					countTable);
	std::string nonUniqueRegionName = "NON_UNIQUE";
	std::unordered_set<uint64_t> nonUniqueKmersPerSet;
	if (!nonUniqueKmerTable.empty()) {
		nonUniqueKmersPerSet = UniqueKmerSetHelper::readInUniqueKmerTableSetsCollapsed(nonUniqueKmerTable);
	}
	out << "set\tuniquerKmers" << std::endl;
	out << "klen\t" << klen << std::endl;
	auto setNames = njh::getSetOfMapKeys(uniqueKmersPerSet);
	for (const auto &setName: setNames) {
		out << setName << "\t" << uniqueKmersPerSet[setName].size() << std::endl;
	}
	out << nonUniqueRegionName << "\t" << nonUniqueKmersPerSet.size() << std::endl;
	return 0;
}

int kmerExpRunner::addToUniqKmersSet(const njh::progutils::CmdArgs & inputCommands){

	bfs::path countTable;
	bfs::path nonUniqueKmerTable;
	bool getRevComp = false;
	std::string regionName;
	OutOptions outOpts;
	KmerGatherer::KmerGathererPars countPars;
	countPars.allUpper_ = false;
	countPars.entropyFilter_ = 0;

	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.description_ = "Add unique kmers to set of already determined unique kmers";
	setUp.setOption(nonUniqueKmerTable, "--nonUniqueKmerTable", "non-unique Kmer Table, 1)set,2)kmer");
	setUp.setOption(getRevComp, "--getRevComp", "Add Rev Comp kmers");

	setUp.setOption(countPars.allowableCharacters_, "--allowableCharacters",
									"Only count kmers with these allowable Characters");
	setUp.setOption(countPars.entropyFilter_, "--entropyFilter", "entropy Filter cut off, exclusive, will only keep kmers above this entropy level");
	setUp.setOption(countPars.allUpper_, "--changeSeqsToUpper", "change all sequences to upper case, otherwise they will filter off");

	setUp.setOption(countTable, "--kmerTable,--countTable", "unique kmer sets, 1)set,2)kmer", true);
	setUp.setOption(regionName, "--regionName", "region name for the input, can be a name already in --kmerTable", true);
	setUp.processReadInNames(VecStr{"--fasta", "--fastagz", "--fastq", "--fastqgz"}, true);
	setUp.processWritingOptions(outOpts);
	bfs::path nonUniqueOutputFnp = njh::files::prependFileBasename(outOpts.outName(), "nonUniqueKmers_");
	setUp.setOption(nonUniqueOutputFnp, "--nonUniqueOutputFnp", "non Unique Output Fnp");
	OutOptions nonUniqKmersOutOpts(nonUniqueOutputFnp);
	nonUniqKmersOutOpts.transferOverwriteOpts(outOpts);
	setUp.finishSetUp(std::cout);

	njh::files::checkExistenceThrow(countTable, __PRETTY_FUNCTION__ );
	if(!nonUniqueKmerTable.empty()){
		njh::files::checkExistenceThrow(nonUniqueKmerTable, __PRETTY_FUNCTION__ );
	}
	OutputStream out(outOpts);
	OutputStream nonUniqKmersOut(nonUniqKmersOutOpts);

	uint32_t klen = UniqueKmerSetHelper::getKmerLenFromUniqueKmerTable(countTable);
	std::unordered_map<std::string, std::unordered_set<uint64_t>> uniqueKmersPerSet = UniqueKmerSetHelper::readInUniqueKmerTablePerSet(
					countTable);
	std::string nonUniqueRegionName = "NON_UNIQUE";
	std::unordered_set<uint64_t> nonUniqueKmersPerSet;
	if (!nonUniqueKmerTable.empty()) {
		nonUniqueKmersPerSet = UniqueKmerSetHelper::readInUniqueKmerTableSetsCollapsed(nonUniqueKmerTable);
	}


	seqInfo seq;
	if(countPars.allUpper_){
		setUp.pars_.ioOptions_.lowerCaseBases_ = "upper";
	}
	std::function<bool(const std::string&)> seqCheck = [&countPars](const std::string & k){
		return std::all_of(k.begin(), k.end(), [&countPars](char base){return njh::in(base, countPars.allowableCharacters_);});
	};
	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	SimpleKmerHash hasher;
	std::unordered_set<uint64_t> rawKmersPerInput;
	while (reader.readNextRead(seq)) {
		for (uint32_t pos = 0; pos < len(seq.seq_) - klen + 1; ++pos) {

			auto k = seq.seq_.substr(pos, klen);
			if (seqCheck(k)) {
				DNABaseCounter counter(countPars.allowableCharacters_);
				counter.increase(k);
				if (counter.computeEntrophy() > countPars.entropyFilter_) {
					rawKmersPerInput.emplace(hasher.hash(k));
				}
			}
		}
		if (getRevComp) {
			for (uint32_t pos = 0; pos < len(seq.seq_) - klen + 1;
					 ++pos) {
				auto k = seq.seq_.substr(pos, klen);
				if (seqCheck(k)) {
					DNABaseCounter counter(countPars.allowableCharacters_);
					counter.increase(k);
					if (counter.computeEntrophy() > countPars.entropyFilter_) {
						rawKmersPerInput.emplace(
										hasher.revCompHash(k));
					}
				}
			}
		}
	}
	std::unordered_map<std::string, std::unordered_set<uint64_t>> outputUniqueKmersPerSet;

	if(njh::in(regionName, uniqueKmersPerSet)){
		outputUniqueKmersPerSet[regionName] = uniqueKmersPerSet[regionName];
	}
	for (const auto finalNewKmer: rawKmersPerInput) {
		bool pass = true;
		if (njh::in(finalNewKmer, nonUniqueKmersPerSet)) {
			pass = false;
		} else {
			for (const auto &set: uniqueKmersPerSet) {
				if (set.first == regionName) {
					continue;
				}
				if (njh::in(finalNewKmer, set.second)) {
					pass = false;
					nonUniqueKmersPerSet.emplace(finalNewKmer);
					break;
				}
			}
		}
		if (pass) {
			outputUniqueKmersPerSet[regionName].emplace(finalNewKmer);
		}
	}

	//filter other sets
	for (const auto &set: uniqueKmersPerSet) {
		if (set.first == regionName) {
			continue;
		}
		for(const auto & finalKmer : set.second){
			if(!njh::in(finalKmer, rawKmersPerInput)){
				outputUniqueKmersPerSet[set.first].emplace(finalKmer);
			}
		}
	}

	for(const auto & finalKmerSets : outputUniqueKmersPerSet){
		for(const auto & finalKmer : finalKmerSets.second){
			out << finalKmerSets.first << "\t" << hasher.reverseHash(finalKmer) << std::endl;
		}
	}

	for(const auto & finalKmer : nonUniqueKmersPerSet){
		nonUniqKmersOut << nonUniqueRegionName << "\t" << hasher.reverseHash(finalKmer) << std::endl;
	}


	return 0;
}

int kmerExpRunner::findUniqKmersFromGenomeSubRegionsMultiple(const njh::progutils::CmdArgs & inputCommands){
	std::string nonUniqueRegionName = "NON_UNIQUE";
	KmerGatherer::KmerGathererPars countPars;
	countPars.entropyFilter_ = 0;
	countPars.kmerLengthForEntropyCalc_ = 2;
	bfs::path genomeFnp;
	bfs::path regionTableFnp;
	bool getReverseCompOfInputRegions = false;
	bool getReverseCompOfGenomeRegions = false;
	bool separateGenomeOutput = false;
	OutOptions outOpts;
	bool filterGenomeRestForEntropyToo = false;
	uint32_t expand = 0;

	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.description_ = "Get the unique kmers that appear within the bed region compared to the rest of the genome";
	setUp.setOption(expand, "--expand", "expand this amount out around the regions when creating the rest of genome unique kmers");
	setUp.setOption(filterGenomeRestForEntropyToo, "--filterGenomeRestForEntropyToo", "filter Genome Rest For low Entropy kmers Too");

	setUp.setOption(genomeFnp, "--genomeFnp", "genome file to extract from, a .2bit file needs to exist for the supplied genome", true);
	setUp.setOption(regionTableFnp, "--bedFnp", "sub regions to extract to compare to the rest of the genome", true);


	setUp.setOption(getReverseCompOfInputRegions, "--getReverseCompOfInputRegions", "get Reverse Comp Of Input Regions");
	setUp.setOption(getReverseCompOfGenomeRegions, "--getReverseCompOfGenomeRegions", "get Reverse Comp Of Genome Regions");

	setUp.processWritingOptions(outOpts);
	bfs::path nonUniqueOutputFnp = njh::files::prependFileBasename(outOpts.outName(), "nonUniqueKmers_");
	setUp.setOption(nonUniqueOutputFnp, "--nonUniqueOutputFnp", "non Unique Output Fnp");
	OutOptions nonUniqKmersOutOpts(nonUniqueOutputFnp);
	nonUniqKmersOutOpts.transferOverwriteOpts(outOpts);

	bfs::path restOfGenomeKmersOutputFnp = njh::files::prependFileBasename(outOpts.outName(), "restOfGenomeKmers_");

	if(setUp.setOption(separateGenomeOutput, "--separateGenomeOutput", "create a separate file for the genome output")) {
		setUp.setOption(restOfGenomeKmersOutputFnp, "--restGenomeOutputFnp", "rest of genome Output Fnp");
	}
	OutOptions restOfGenomeKmersOutOpts(restOfGenomeKmersOutputFnp);
	restOfGenomeKmersOutOpts.transferOverwriteOpts(outOpts);
	std::string restOfGenomeRegionName = "genomeRest";
	setUp.setOption(restOfGenomeRegionName, "--restOfGenomeRegionName", "rest of genome kmer set name");
	if(njh::containsSpecialChars(restOfGenomeRegionName)){
		setUp.failed_ = true;
		setUp.addWarning("restOfGenomeRegionName can't have special characters");
	}
	if(njh::strHasWhitesapce(restOfGenomeRegionName)){
		setUp.failed_ = true;
		setUp.addWarning("restOfGenomeRegionName can't have whitespace characters");
	}

	setUp.setOption(nonUniqueRegionName, "--nonUniqueRegionName", "non unique kmer set name");
	if(njh::containsSpecialChars(nonUniqueRegionName)){
		setUp.failed_ = true;
		setUp.addWarning("nonUniqueRegionName can't have special characters");
	}
	if(njh::strHasWhitesapce(nonUniqueRegionName)){
		setUp.failed_ = true;
		setUp.addWarning("nonUniqueRegionName can't have whitespace characters");
	}
	countPars.setOptions(setUp);
	setUp.finishSetUp(std::cout);

	KmerGatherer kGather(countPars);

	//set up 2bit file
	bfs::path genome2bitFnp = bfs::path(genomeFnp).replace_extension("").string() + ".2bit";
	njh::files::checkExistenceThrow(genome2bitFnp, __PRETTY_FUNCTION__ );
	TwoBit::TwoBitFile tReader(genome2bitFnp);
	auto chromLengths = tReader.getSeqLens();

	std::unordered_map<std::string, std::vector<std::shared_ptr<Bed6RecordCore>>> allBedsByRegion;
	//get input regions
	std::vector<std::shared_ptr<Bed6RecordCore>> allBeds;
	table regionTab(regionTableFnp, "\t", true);
	regionTab.checkForColumnsThrow(VecStr{"set", "bed"},__PRETTY_FUNCTION__ );
	{
		auto setColPos = regionTab.getColPos("set");
		auto bedColPos = regionTab.getColPos("bed");
		VecStr warnings;
		for(const auto & row : regionTab){
			auto regionName = row[setColPos];
			if(njh::containsSpecialChars(regionName)){
				warnings.emplace_back("region name can't have special characters");
			}
			if(njh::strHasWhitesapce(regionName)){
				warnings.emplace_back("region name can't have whitespace characters");
			}
			if(njh::in(regionName, allBedsByRegion)){
				warnings.emplace_back(njh::pasteAsStr("set name : " , regionName, " already read in"));
			}
			auto currentBeds = getBeds(row[bedColPos]);
			allBedsByRegion[regionName] = currentBeds;
			njh::addOtherVec(allBeds, currentBeds);
		}
		if(!warnings.empty()){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error: "  << "\n";
			ss << njh::conToStr(warnings, "\n");
			throw std::runtime_error{ss.str()};
		}
	}

	BedUtility::coordSort(allBeds);

	////check in put
	////// check for overlaping
	VecStr overlappingRegionWarnings;
	if(allBeds.size() > 1){
		for(const auto pos : iter::range(1UL, allBeds.size())){
			if(allBeds[pos]->overlaps(*allBeds[pos - 1], 1)){
				overlappingRegionWarnings.emplace_back(njh::pasteAsStr(allBeds[pos]->name_, " overlaps ", allBeds[pos - 1]->name_));
			}
		}
	}
	if(!overlappingRegionWarnings.empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "found the following overlaps" << "\n";
		ss << njh::conToStr(overlappingRegionWarnings, "\n") << "\n";
		throw std::runtime_error{ss.str()};
	}

	////// chroms
	VecStr chromsDoNotExist;
	for(const auto & b : allBeds){
		if(!njh::in(b->chrom_, chromLengths)){
			chromsDoNotExist.emplace_back(b->chrom_);
		}
	}
	njh::sort(chromsDoNotExist);
	if(!chromsDoNotExist.empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "the following chroms could not be found within " << genome2bitFnp << "\n";
		ss << njh::conToStr(chromsDoNotExist, ", ") << "\n";
		throw std::runtime_error{ss.str()};
	}
	std::unordered_map<std::string, std::vector<std::shared_ptr<Bed6RecordCore>>> bedsByChromForSubstracting;

	std::vector<std::shared_ptr<Bed6RecordCore>> allBedsForSubtracting;
	for (const auto &region: allBeds) {
		Bed6RecordCore regionToAdd = *region;
		if (expand > 0) {
			BedUtility::extendLeftRight(regionToAdd, expand, expand, chromLengths.at(regionToAdd.chrom_));
		}
		if (!allBedsForSubtracting.empty() && allBedsForSubtracting.back()->overlaps(regionToAdd, 1)) {
			allBedsForSubtracting.back()->chromEnd_ = std::max(allBedsForSubtracting.back()->chromEnd_,
																												 regionToAdd.chromEnd_);
		} else {
			allBedsForSubtracting.emplace_back(std::make_shared<Bed6RecordCore>(regionToAdd));
		}
	}

	for(const auto & bedRegion : allBedsForSubtracting){
		bedsByChromForSubstracting[bedRegion->chrom_].emplace_back(bedRegion);
	}

	//// subtract region to create other regions
	std::vector<Bed6RecordCore> inbetweenRegions;
	for (const auto &chromKey: njh::getSetOfMapKeys(chromLengths)) {

		if (njh::in(chromKey, bedsByChromForSubstracting)) {
			//add front
			if (0 != bedsByChromForSubstracting[chromKey].front()->chromStart_) {
				inbetweenRegions.emplace_back(chromKey, 0, bedsByChromForSubstracting[chromKey].front()->chromStart_,
																			njh::pasteAsStr(chromKey, "-", 0, "-",
																											bedsByChromForSubstracting[chromKey].front()->chromStart_),
																			bedsByChromForSubstracting[chromKey].front()->chromStart_, '+');
			}
			if(bedsByChromForSubstracting[chromKey].size() > 1){
				for(const auto pos : iter::range(bedsByChromForSubstracting[chromKey].size() - 1)){
					auto start = bedsByChromForSubstracting[chromKey][pos]->chromEnd_;
					auto end = bedsByChromForSubstracting[chromKey][pos + 1]->chromStart_;
					if(start != end){
						inbetweenRegions.emplace_back(chromKey, start, end,
																					njh::pasteAsStr(chromKey, "-", start, "-",
																													end),
																					end - start, '+');
					}
				}
			}
			//add back
			if (chromLengths[chromKey] != bedsByChromForSubstracting[chromKey].back()->chromEnd_) {
				inbetweenRegions.emplace_back(chromKey, bedsByChromForSubstracting[chromKey].back()->chromEnd_, chromLengths[chromKey],
																			njh::pasteAsStr(chromKey, "-", bedsByChromForSubstracting[chromKey].back()->chromEnd_, "-",
																											chromLengths[chromKey]),
																			chromLengths[chromKey] - bedsByChromForSubstracting[chromKey].back()->chromEnd_, '+');
			}
		} else {
			//no bed regions for this chrom
			inbetweenRegions.emplace_back(chromKey, 0, chromLengths[chromKey], chromKey, chromLengths[chromKey], '+');
		}
	}

	//// output
	OutputStream out(outOpts);

	OutputStream nonUniqKmersOut(nonUniqKmersOutOpts);
	//OutputStream restOfGenomeKmersOut(restOfGenomeKmersOutOpts);
	OutputStreamWrap restOfGenomeKmersOut(restOfGenomeKmersOutOpts);
	if(separateGenomeOutput){
		restOfGenomeKmersOut.openOut();
	}
	std::unordered_set<uint64_t> nonUniqueKmers;
	std::unordered_set<uint64_t> restOfGenomeUniqueKmers;

	//// hasher
	SimpleKmerHash hasher;
	//// get kmers for input bed regions
	std::map<std::string, std::unordered_set<uint64_t>> rawKmersPerInput;
	for (const auto &regionBeds: allBedsByRegion) {
		for (const auto &region: regionBeds.second) {
			auto currentSeq = GenomicRegion(*region).extractSeq(tReader);
			for (uint32_t pos = 0; pos < len(currentSeq.seq_) - countPars.kmerLength_ + 1; ++pos) {
				rawKmersPerInput[regionBeds.first].emplace(hasher.hash(currentSeq.seq_.substr(pos, countPars.kmerLength_)));
//				auto k = currentSeq.seq_.substr(pos, countPars.kmerLength_);
//				kmerInfo kInfo(k, countPars.kmerLengthForEntropyCalc_, false);
//				if(kInfo.computeKmerEntropy() > countPars.entropyFilter_){
//					rawKmersPerInput[regionBeds.first].emplace(hasher.hash(k));
//				}
			}
			if (getReverseCompOfInputRegions) {
				for (uint32_t pos = 0; pos < len(currentSeq.seq_) - countPars.kmerLength_ + 1;
						 ++pos) {
					rawKmersPerInput[regionBeds.first].emplace(hasher.revCompHash(currentSeq.seq_.substr(pos, countPars.kmerLength_)));
//					auto k = currentSeq.seq_.substr(pos, countPars.kmerLength_);
//					kmerInfo kInfo(k, countPars.kmerLengthForEntropyCalc_, false);
//					if(kInfo.computeKmerEntropy() > countPars.entropyFilter_){
//						rawKmersPerInput[regionBeds.first].emplace(hasher.revCompHash(k));
//					}
				}
			}
		}
	}
	std::map<std::string, std::unordered_set<uint64_t>> finalKmersPerInput;
	//initial filter to other groups
	if (rawKmersPerInput.size() > 1) {
		finalKmersPerInput = rawKmersPerInput;
	} else {
		for(const auto & kmerSet : rawKmersPerInput){
			for(const auto & kmer : kmerSet.second){
				bool pass = true;
				for(const auto & otherKmerSet : rawKmersPerInput){
					if(otherKmerSet.first != kmerSet.first){
						if(njh::in(kmer, otherKmerSet.second)){
							pass = false;
							nonUniqueKmers.emplace(kmer);
							break;
						}
					}
				}
				if(pass){
					finalKmersPerInput[kmerSet.first].emplace(kmer);
				}
			}
		}
	}

	//// get kmers for the rest of the genome

	for(const auto & region : inbetweenRegions){
		if(setUp.pars_.verbose_){
			std::cout << region.name_ << std::endl;
			for(const auto & finalKmerSet : finalKmersPerInput){
				std::cout << "\t" << "prior to filter for set: " << finalKmerSet.first  << " " << finalKmerSet.second.size() << std::endl;
			}
			std::cout << "\t" << "nonUniqueKmers: " << nonUniqueKmers.size()  << std::endl;
			std::cout << "\t" << "restOfGenomeUniqueKmers: " << restOfGenomeUniqueKmers.size() << std::endl << std::endl;
		}
		std::set<uint64_t> kmersPerInbetween;
		auto currentSeq = GenomicRegion(region).extractSeq(tReader);

		for (uint32_t pos = 0; pos < len(currentSeq.seq_) - countPars.kmerLength_ + 1; ++pos) {
			kmersPerInbetween.emplace(hasher.hash(currentSeq.seq_.substr(pos, countPars.kmerLength_)));
//			auto k = currentSeq.seq_.substr(pos, countPars.kmerLength_);
//			kmerInfo kInfo(k, countPars.kmerLengthForEntropyCalc_, false);
//			if(kInfo.computeKmerEntropy() > countPars.entropyFilter_){
//				kmersPerInbetween.emplace(hasher.hash(k));
//			}
		}
		if (getReverseCompOfGenomeRegions) {
			for (uint32_t pos = 0; pos < len(currentSeq.seq_) - countPars.kmerLength_ + 1;
					 ++pos) {
				kmersPerInbetween.emplace(hasher.revCompHash(currentSeq.seq_.substr(pos, countPars.kmerLength_)));
//				auto k = currentSeq.seq_.substr(pos, countPars.kmerLength_);
//				kmerInfo kInfo(k, countPars.kmerLengthForEntropyCalc_, false);
//				if(kInfo.computeKmerEntropy() > countPars.entropyFilter_){
//					kmersPerInbetween.emplace(hasher.revCompHash(k));
//				}
			}
		}

		//filter kmers
		std::map<std::string, std::unordered_set<uint64_t> > filterKmers;
		for (const auto &finalKmerSet: finalKmersPerInput) {
			for (const auto &finalKmer: finalKmerSet.second) {
				if (!njh::in(finalKmer, kmersPerInbetween)) {
					filterKmers[finalKmerSet.first].emplace(finalKmer);
				} else {
					nonUniqueKmers.emplace(finalKmer);
				}
			}
		}
		finalKmersPerInput = filterKmers;
		for(const auto & finalKmer: kmersPerInbetween){
			if(!njh::in(finalKmer, nonUniqueKmers)){
				restOfGenomeUniqueKmers.emplace(finalKmer);
			}
		}
		if(setUp.pars_.verbose_){
			for(const auto & finalKmerSet : finalKmersPerInput){
				std::cout << "\t" << "set: " << finalKmerSet.first  << " " << finalKmerSet.second.size() << std::endl;
			}
			std::cout << "\t" << "nonUniqueKmers: " << nonUniqueKmers.size()  << std::endl;
			std::cout << "\t" << "restOfGenomeUniqueKmers: " << restOfGenomeUniqueKmers.size() << std::endl << std::endl;

		}
	}
	if(separateGenomeOutput){
		for(const auto & finalKmer : restOfGenomeUniqueKmers){
			if(filterGenomeRestForEntropyToo){
				auto k = hasher.reverseHash(finalKmer);
				kmerInfo kInfo(k, countPars.kmerLengthForEntropyCalc_, false);
				if (kInfo.computeKmerEntropy() > countPars.entropyFilter_) {
					*restOfGenomeKmersOut.out_ << restOfGenomeRegionName << "\t" << k << std::endl;
				}
			} else {
				*restOfGenomeKmersOut.out_ << restOfGenomeRegionName << "\t" << hasher.reverseHash(finalKmer) << std::endl;
			}
		}
	} else {
		for(const auto & finalKmer : restOfGenomeUniqueKmers){
			if(filterGenomeRestForEntropyToo){
				auto k = hasher.reverseHash(finalKmer);
				kmerInfo kInfo(k, countPars.kmerLengthForEntropyCalc_, false);
				if (kInfo.computeKmerEntropy() > countPars.entropyFilter_) {
					out << restOfGenomeRegionName << "\t" << k << std::endl;
				}
			} else {
				out << restOfGenomeRegionName << "\t" << hasher.reverseHash(finalKmer) << std::endl;
			}
		}
	}
	for(const auto & finalKmerSet : finalKmersPerInput){
		for(const auto & finalKmer : finalKmerSet.second){
			auto k = hasher.reverseHash(finalKmer);
			kmerInfo kInfo(k, countPars.kmerLengthForEntropyCalc_, false);
			if (kInfo.computeKmerEntropy() > countPars.entropyFilter_) {
				out << finalKmerSet.first << "\t" << k << std::endl;
			}
		}
	}
	for(const auto & finalKmer : nonUniqueKmers){
		nonUniqKmersOut << nonUniqueRegionName << "\t" << hasher.reverseHash(finalKmer) << std::endl;
	}

	return 0;
}


//
//int kmerExpRunner::findUniqKmersFromGenomeSubRegions(const njh::progutils::CmdArgs & inputCommands){
//	std::string nonUniqueRegionName = "NON_UNIQUE";
//	KmerGatherer::KmerGathererPars countPars;
//	countPars.entropyFilter_ = 1.20;
//	countPars.kmerLengthForEntropyCalc_ = 2;
//	bfs::path genomeFnp;
//	bfs::path bedFnp;
//	std::string regionName;
//	bool getReverseCompOfInputRegions = false;
//	bool getReverseCompOfGenomeRegions = false;
//	OutOptions outOpts;
//	seqSetUp setUp(inputCommands);
//	setUp.processVerbose();
//	setUp.processDebug();
//	setUp.description_ = "Get the unique kmers that appear within the bed region compared to the rest of the genome";
//
//	setUp.setOption(genomeFnp, "--genomeFnp", "genome file to extract from, a .2bit file needs to exist for the supplied genome", true);
//	setUp.setOption(bedFnp, "--bedFnp", "sub regions to extract to compare to the rest of the genome", true);
//	setUp.setOption(regionName, "--regionName", "region name to give to the unique kmer set", true);
//	if(njh::containsSpecialChars(regionName)){
//		setUp.failed_ = true;
//		setUp.addWarning("region name can't have special characters");
//	}
//	if(njh::strHasWhitesapce(regionName)){
//		setUp.failed_ = true;
//		setUp.addWarning("region name can't have whitespace characters");
//	}
//	setUp.setOption(getReverseCompOfInputRegions, "--getReverseCompOfInputRegions", "get Reverse Comp Of Input Regions");
//	setUp.setOption(getReverseCompOfGenomeRegions, "--getReverseCompOfGenomeRegions", "get Reverse Comp Of Genome Regions");
//
//	setUp.processWritingOptions(outOpts);
//	bfs::path nonUniqueOutputFnp = njh::files::prependFileBasename(outOpts.outName(), "nonUniqueKmers_");
//	setUp.setOption(nonUniqueOutputFnp, "--nonUniqueOutputFnp", "non Unique Output Fnp");
//	OutOptions nonUniqKmersOutOpts(nonUniqueOutputFnp);
//	nonUniqKmersOutOpts.transferOverwriteOpts(outOpts);
//
//	bfs::path restOfGenomeKmersOutputFnp = njh::files::prependFileBasename(outOpts.outName(), "restOfGenomeKmers_");
//	setUp.setOption(restOfGenomeKmersOutputFnp, "--restGenomeOutputFnp", "rest of genome Output Fnp");
//	OutOptions restOfGenomeKmersOutOpts(restOfGenomeKmersOutputFnp);
//	restOfGenomeKmersOutOpts.transferOverwriteOpts(outOpts);
//	std::string restOfGenomeRegionName = "genomeMinus_" + regionName;
//	setUp.setOption(restOfGenomeRegionName, "--restOfGenomeRegionName", "rest of genome kmer set name");
//	if(njh::containsSpecialChars(restOfGenomeRegionName)){
//		setUp.failed_ = true;
//		setUp.addWarning("restOfGenomeRegionName can't have special characters");
//	}
//	if(njh::strHasWhitesapce(restOfGenomeRegionName)){
//		setUp.failed_ = true;
//		setUp.addWarning("restOfGenomeRegionName can't have whitespace characters");
//	}
//
//	setUp.setOption(nonUniqueRegionName, "--nonUniqueRegionName", "non unique kmer set name");
//	if(njh::containsSpecialChars(nonUniqueRegionName)){
//		setUp.failed_ = true;
//		setUp.addWarning("nonUniqueRegionName can't have special characters");
//	}
//	if(njh::strHasWhitesapce(nonUniqueRegionName)){
//		setUp.failed_ = true;
//		setUp.addWarning("nonUniqueRegionName can't have whitespace characters");
//	}
//
//	countPars.setOptions(setUp);
//	setUp.finishSetUp(std::cout);
//
//	KmerGatherer kGather(countPars);
//
//	//set up 2bit file
//	bfs::path genome2bitFnp = bfs::path(genomeFnp).replace_extension("").string() + ".2bit";
//	njh::files::checkExistenceThrow(genome2bitFnp, __PRETTY_FUNCTION__ );
//	TwoBit::TwoBitFile tReader(genome2bitFnp);
//	auto chromLengths = tReader.getSeqLens();
//
//	//get input regions
//	auto bedRegions = getBeds(bedFnp);
//	BedUtility::coordSort(bedRegions);
//
//	////check in put
//	////// check for overlaping
//	VecStr overlappingRegionWarnings;
//	if(bedRegions.size() > 1){
//		for(const auto pos : iter::range(1UL, bedRegions.size())){
//			if(bedRegions[pos]->overlaps(*bedRegions[pos - 1], 1)){
//				overlappingRegionWarnings.emplace_back(njh::pasteAsStr(bedRegions[pos]->name_, " overlaps ", bedRegions[pos - 1]->name_));
//			}
//		}
//	}
//	if(!overlappingRegionWarnings.empty()){
//		std::stringstream ss;
//		ss << __PRETTY_FUNCTION__ << ", error " << "found the following overlaps" << "\n";
//		ss << njh::conToStr(overlappingRegionWarnings, "\n") << "\n";
//		throw std::runtime_error{ss.str()};
//	}
//
//	////// chroms
//	VecStr chromsDoNotExist;
//	for(const auto & b : bedRegions){
//		if(!njh::in(b->chrom_, chromLengths)){
//			chromsDoNotExist.emplace_back(b->chrom_);
//		}
//	}
//	if(!chromsDoNotExist.empty()){
//		std::stringstream ss;
//		ss << __PRETTY_FUNCTION__ << ", error " << "the following chroms could not be found within " << genome2bitFnp << "\n";
//		ss << njh::conToStr(chromsDoNotExist, ", ") << "\n";
//		throw std::runtime_error{ss.str()};
//	}
//	std::unordered_map<std::string, std::vector<std::shared_ptr<Bed6RecordCore>>> bedsByChrom;
//	for(const auto & bedRegion : bedRegions){
//		bedsByChrom[bedRegion->chrom_].emplace_back(bedRegion);
//	}
//
//	//// subtract region to create other regions
//	std::vector<Bed6RecordCore> inbetweenRegions;
//	for (const auto &chrom: chromLengths) {
//		if (njh::in(chrom.first, bedsByChrom)) {
//			//add front
//			if (0 != bedsByChrom[chrom.first].front()->chromStart_) {
//				inbetweenRegions.emplace_back(chrom.first, 0, bedsByChrom[chrom.first].front()->chromStart_,
//																			njh::pasteAsStr(chrom.first, "-", 0, "-",
//																											bedsByChrom[chrom.first].front()->chromStart_),
//																			bedsByChrom[chrom.first].front()->chromStart_, '+');
//			}
//			if(bedsByChrom[chrom.first].size() > 1){
//				for(const auto pos : iter::range(bedsByChrom[chrom.first].size() - 1)){
//					auto start = bedsByChrom[chrom.first][pos]->chromEnd_;
//					auto end = bedsByChrom[chrom.first][pos + 1]->chromStart_;
//					if(start != end){
//						inbetweenRegions.emplace_back(chrom.first, start, end,
//																					njh::pasteAsStr(chrom.first, "-", start, "-",
//																													end),
//																					end - start, '+');
//					}
//				}
//			}
//			//add back
//			if (chrom.second != bedsByChrom[chrom.first].back()->chromEnd_) {
//				inbetweenRegions.emplace_back(chrom.first, bedsByChrom[chrom.first].back()->chromEnd_, chrom.second,
//																			njh::pasteAsStr(chrom.first, "-", bedsByChrom[chrom.first].back()->chromEnd_, "-",
//																											chrom.second),
//																			chrom.second - bedsByChrom[chrom.first].back()->chromEnd_, '+');
//			}
//		} else {
//			//no bed regions for this chrom
//			inbetweenRegions.emplace_back(chrom.first, 0, chrom.second, chrom.first, chrom.second, '+');
//		}
//	}
//
//	//// output
//	OutputStream out(outOpts);
//
//	OutputStream nonUniqKmersOut(nonUniqKmersOutOpts);
//	OutputStream restOfGenomeKmersOut(restOfGenomeKmersOutOpts);
//
//	std::unordered_set<uint64_t> nonUniqueKmers;
//	std::unordered_set<uint64_t> restOfGenomeUniqueKmers;
//
//	//// hasher
//	SimpleKmerHash hasher;
//
//	//// get kmers for input bed regions
//	std::set<uint64_t> rawKmersPerInput;
//	for (const auto &region: bedRegions) {
//		auto currentSeq = GenomicRegion(*region).extractSeq(tReader);
//
//		for (uint32_t pos = 0; pos < len(currentSeq.seq_) - countPars.kmerLength_ + 1; ++pos) {
//			auto k = currentSeq.seq_.substr(pos, countPars.kmerLength_);
//			kmerInfo kInfo(k, countPars.kmerLengthForEntropyCalc_, false);
//			if(kInfo.computeKmerEntropy() > countPars.entropyFilter_){
//				rawKmersPerInput.emplace(hasher.hash(k));
//			}
//		}
//		if (getReverseCompOfInputRegions) {
//			for (uint32_t pos = 0; pos < len(currentSeq.seq_) - countPars.kmerLength_ + 1;
//					 ++pos) {
//				auto k = currentSeq.seq_.substr(pos, countPars.kmerLength_);
//				kmerInfo kInfo(k, countPars.kmerLengthForEntropyCalc_, false);
//				if(kInfo.computeKmerEntropy() > countPars.entropyFilter_){
//					rawKmersPerInput.emplace(hasher.revCompHash(k));
//				}
//			}
//		}
//	}
//
//	std::set<uint64_t> finalKmersPerInput = rawKmersPerInput;
//
//	//// get kmers for the rest of the genome
//
//	for(const auto & region : inbetweenRegions){
//		if(setUp.pars_.verbose_){
//			std::cout << region.name_ << std::endl;
//			std::cout << "\t" << finalKmersPerInput.size() << std::endl;
//		}
//		std::set<uint64_t> kmersPerInbetween;
//		auto currentSeq = GenomicRegion(region).extractSeq(tReader);
//
//		for (uint32_t pos = 0; pos < len(currentSeq.seq_) - countPars.kmerLength_ + 1; ++pos) {
//			auto k = currentSeq.seq_.substr(pos, countPars.kmerLength_);
//			kmerInfo kInfo(k, countPars.kmerLengthForEntropyCalc_, false);
//			if(kInfo.computeKmerEntropy() > countPars.entropyFilter_){
//				kmersPerInbetween.emplace(hasher.hash(k));
//			}
//		}
//		if (getReverseCompOfGenomeRegions) {
//			for (uint32_t pos = 0; pos < len(currentSeq.seq_) - countPars.kmerLength_ + 1;
//					 ++pos) {
//				auto k = currentSeq.seq_.substr(pos, countPars.kmerLength_);
//				kmerInfo kInfo(k, countPars.kmerLengthForEntropyCalc_, false);
//				if(kInfo.computeKmerEntropy() > countPars.entropyFilter_){
//					kmersPerInbetween.emplace(hasher.revCompHash(k));
//				}
//			}
//		}
//
//		//filter kmers
//		std::set < uint64_t > filterKmers;
//		for (const auto &finalKmer: finalKmersPerInput) {
//			if (!njh::in(finalKmer, kmersPerInbetween)) {
//				filterKmers.emplace(finalKmer);
//			} else {
//				nonUniqueKmers.emplace(finalKmer);
//			}
//		}
//		finalKmersPerInput = filterKmers;
//		for(const auto & finalKmer: kmersPerInbetween){
//			if(!njh::in(finalKmer, nonUniqueKmers)){
//				restOfGenomeUniqueKmers.emplace(finalKmer);
//			}
//		}
//		if(setUp.pars_.verbose_){
//			std::cout << "\t" << "finalKmersPerInput: " << finalKmersPerInput.size() << std::endl << std::endl;
//			std::cout << "\t" << "restOfGenomeUniqueKmers: " << restOfGenomeUniqueKmers.size() << std::endl << std::endl;
//		}
//	}
//
//	for(const auto & finalKmer : restOfGenomeUniqueKmers){
//		restOfGenomeKmersOut << restOfGenomeRegionName << "\t" << hasher.reverseHash(finalKmer) << std::endl;
//	}
//	for(const auto & finalKmer : finalKmersPerInput){
//		out << regionName << "\t" << hasher.reverseHash(finalKmer) << std::endl;
//	}
//	for(const auto & finalKmer : nonUniqueKmers){
//		nonUniqKmersOut << nonUniqueRegionName << "\t" << hasher.reverseHash(finalKmer) << std::endl;
//	}
//
//	return 0;
//}

}  //namespace njhseq

