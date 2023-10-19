//
// Created by Nicholas Hathaway on 5/25/23.
//

#include <njhseq/readVectorManipulation/readVectorHelpers/readVecChecker.hpp>
#include "UniqueKmerSetHelper.hpp"


namespace njhseq {


UniqueKmerSetHelper::CompareReadToSetRes UniqueKmerSetHelper::compareReadToSets(PairedRead &pseq,
																																								const std::unordered_map<std::string, std::unordered_set<uint64_t>> &uniqueKmersPerSet,
																																								const CompareReadToSetPars &compPars,
																																								const SimpleKmerHash &hasher) {
	CompareReadToSetRes ret;
	if (len(pseq.seqBase_.seq_) > compPars.klen) {
		for (uint32_t pos = 0; pos < len(pseq.seqBase_.seq_) - compPars.klen + 1; ++pos) {
			++ret.hashedInputKmers[hasher.hash(pseq.seqBase_.seq_.substr(pos, compPars.klen))];
		}
	}
	if (len(pseq.mateSeqBase_.seq_) > compPars.klen) {
		for (uint32_t pos = 0; pos < len(pseq.mateSeqBase_.seq_) - compPars.klen + 1; ++pos) {
			++ret.hashedInputKmers[hasher.hash(pseq.mateSeqBase_.seq_.substr(pos, compPars.klen))];
		}
	}
	if (compPars.includeRevComp) {
		if (len(pseq.seqBase_.seq_) > compPars.klen) {
			for (uint32_t pos = 0; pos < len(pseq.seqBase_.seq_) - compPars.klen + 1; ++pos) {
				++ret.hashedInputKmersRevComp[hasher.revCompHash(pseq.seqBase_.seq_.substr(pos, compPars.klen))];
			}
		}
		if (len(pseq.mateSeqBase_.seq_) > compPars.klen) {
			for (uint32_t pos = 0; pos < len(pseq.mateSeqBase_.seq_) - compPars.klen + 1; ++pos) {
				++ret.hashedInputKmersRevComp[hasher.revCompHash(pseq.mateSeqBase_.seq_.substr(pos, compPars.klen))];
			}
		}
	}


	for (const auto &hashedKmer: ret.hashedInputKmers) {
		for (const auto &uniqueKmers: uniqueKmersPerSet) {
			if (njh::in(hashedKmer.first, uniqueKmers.second)) {
				++ret.foundPerSet[uniqueKmers.first];
			}
		}
	}

	if (compPars.includeRevComp) {
		for (const auto &hashedKmer: ret.hashedInputKmersRevComp) {
			for (const auto &uniqueKmers: uniqueKmersPerSet) {
				if (njh::in(hashedKmer.first, uniqueKmers.second)) {
					++ret.foundPerSetRevComp[uniqueKmers.first];
				}
			}
		}
	}

	//set a hard cut off, reset counts to zero
	for (auto &perSet: ret.foundPerSet) {
		if (perSet.second < compPars.hardCountOff) {
			perSet.second = 0;
		}
	}
	for (auto &perSet: ret.foundPerSetRevComp) {
		if (perSet.second < compPars.hardCountOff) {
			perSet.second = 0;
		}
	}

	uint32_t kmersPossible = (pseq.seqBase_.seq_.size() - compPars.klen + 1) + (pseq.mateSeqBase_.seq_.size() - compPars.klen + 1);
	//set a hard cut off or the fraction of possible kmers is less than frac cut off, reset counts to zero
	for (auto &perSet: ret.foundPerSet) {
		if (perSet.second < compPars.hardCountOff || perSet.second / static_cast<double>(kmersPossible) < compPars.fracCutOff) {
			perSet.second = 0;
		}
	}
	for (auto &perSet: ret.foundPerSetRevComp) {
		if (perSet.second < compPars.hardCountOff || perSet.second / static_cast<double>(kmersPossible) < compPars.fracCutOff) {
			perSet.second = 0;
		}
	}

	for (const auto &setName: njh::getVecOfMapKeys(uniqueKmersPerSet)) {
		if (static_cast<double>(ret.foundPerSet[setName]) / static_cast<double>(ret.hashedInputKmers.size()) >
				ret.bestFrac) {
			ret.bestFrac = static_cast<double>(ret.foundPerSet[setName]) / static_cast<double>(ret.hashedInputKmers.size());
			ret.winnerSet = setName;
			ret.winnerRevComp = false;
		}
		if (compPars.includeRevComp) {
			if (static_cast<double>(ret.foundPerSetRevComp[setName]) / static_cast<double>(ret.hashedInputKmers.size()) >
					ret.bestFrac) {
				ret.bestFrac = static_cast<double>(ret.foundPerSetRevComp[setName]) /
											 static_cast<double>(ret.hashedInputKmers.size());
				ret.winnerSet = setName;
				ret.winnerRevComp = true;
			}
		}
	}

	if (ret.winnerRevComp) {
		pseq.seqBase_.reverseComplementRead(false, true);
		pseq.mateSeqBase_.reverseComplementRead(false, true);
	}
	return ret;
}

UniqueKmerSetHelper::CompareReadToSetRes UniqueKmerSetHelper::compareReadToSets(seqInfo &seq,
																																								const std::unordered_map<std::string, std::unordered_set<uint64_t>> &uniqueKmersPerSet,
																																								const CompareReadToSetPars &compPars,
																																								const SimpleKmerHash &hasher) {
	CompareReadToSetRes ret;
	if (len(seq.seq_) > compPars.klen) {
		for (uint32_t pos = 0; pos < len(seq.seq_) - compPars.klen + 1; ++pos) {
			++ret.hashedInputKmers[hasher.hash(seq.seq_.substr(pos, compPars.klen))];
		}
		if (compPars.includeRevComp) {
			for (uint32_t pos = 0; pos < len(seq.seq_) - compPars.klen + 1; ++pos) {
				++ret.hashedInputKmersRevComp[hasher.revCompHash(seq.seq_.substr(pos, compPars.klen))];
			}
		}
	}


	for (const auto &hashedKmer: ret.hashedInputKmers) {
		for (const auto &uniqueKmers: uniqueKmersPerSet) {
			if (njh::in(hashedKmer.first, uniqueKmers.second)) {
				++ret.foundPerSet[uniqueKmers.first];
			}
		}
	}

	if (compPars.includeRevComp) {
		for (const auto &hashedKmer: ret.hashedInputKmersRevComp) {
			for (const auto &uniqueKmers: uniqueKmersPerSet) {
				if (njh::in(hashedKmer.first, uniqueKmers.second)) {
					++ret.foundPerSetRevComp[uniqueKmers.first];
				}
			}
		}
	}

//	bool print = false;
//	VecStr testReadNames = {"ERR980507.sra.1612064_secondMate","ERR980507.sra.506751_secondMate","ERR980507.sra.1358147_firstMate","ERR980507.sra.101811_secondMate","ERR980507.sra.101803_secondMate"};
//	if(njh::in(seq.name_, testReadNames)){
//		print = true;
//	}

	uint32_t kmersPossible = seq.seq_.size() - compPars.klen + 1;
	//set a hard cut-off or the fraction of possible kmers is less than frac cut off, reset counts to zero

//	if(print){
//		auto foundNames = njh::getSetOfMapKeys(ret.foundPerSet);
//		std::cout << seq.name_ << std::endl;
//		std::cout << "forward" << std::endl;
//		for(const auto & name : foundNames){
//			std::cout << name << " " << ret.foundPerSet.at(name) << " " << ret.foundPerSet.at(name)/static_cast<double>(kmersPossible) << std::endl;
//			std::cout << "\t" << ret.foundPerSet.at(name)/static_cast<double>(ret.hashedInputKmers.size()) << std::endl;
//		}
//		std::cout << "reverse" << std::endl;
//		foundNames = njh::getSetOfMapKeys(ret.foundPerSetRevComp);
//		for(const auto & name : foundNames){
//			std::cout << name << " " << ret.foundPerSetRevComp.at(name) << " " << ret.foundPerSetRevComp.at(name)/static_cast<double>(kmersPossible) << std::endl;
//			std::cout << "\t" << ret.foundPerSetRevComp.at(name)/static_cast<double>(ret.hashedInputKmers.size()) << std::endl;
//		}
//	}

	for (auto &perSet: ret.foundPerSet) {
		if (njh::in(perSet.first, compPars.excludeSetNames)) {
			if (perSet.second < compPars.initialExcludeHardCountOff || perSet.second / static_cast<double>(kmersPossible) < compPars.initialExcludeFracCutOff) {
				perSet.second = 0;
			}
		} else {
			if (perSet.second < compPars.hardCountOff || perSet.second / static_cast<double>(kmersPossible) < compPars.fracCutOff) {
				perSet.second = 0;
			}
		}
	}
	for (auto &perSet: ret.foundPerSetRevComp) {
		if (njh::in(perSet.first, compPars.excludeSetNames)) {
			if (perSet.second < compPars.initialExcludeHardCountOff || perSet.second / static_cast<double>(kmersPossible) < compPars.initialExcludeFracCutOff) {
				perSet.second = 0;
			}
		} else {
			if (perSet.second < compPars.hardCountOff || perSet.second / static_cast<double>(kmersPossible) < compPars.fracCutOff) {
				perSet.second = 0;
			}
		}
	}
//	if(print){
//		auto foundNames = njh::getSetOfMapKeys(ret.foundPerSet);
//		std::cout << "forward" << std::endl;
//		foundNames = njh::getSetOfMapKeys(ret.foundPerSet);
//		for(const auto & name : foundNames){
//			std::cout << name << " " << ret.foundPerSet.at(name) << " " << ret.foundPerSet.at(name)/static_cast<double>(kmersPossible) << std::endl;
//			std::cout << "\t" << ret.foundPerSet.at(name)/static_cast<double>(ret.hashedInputKmers.size()) << std::endl;
//		}
//		std::cout << "reverse" << std::endl;
//		foundNames = njh::getSetOfMapKeys(ret.foundPerSetRevComp);
//		for(const auto & name : foundNames){
//			std::cout << name << " " << ret.foundPerSetRevComp.at(name) << " " << ret.foundPerSetRevComp.at(name)/static_cast<double>(kmersPossible) << std::endl;
//			std::cout << "\t" << ret.foundPerSetRevComp.at(name)/static_cast<double>(ret.hashedInputKmers.size()) << std::endl;
//		}
//	}

	for (const auto &setName: njh::getVecOfMapKeys(uniqueKmersPerSet)) {
		if (static_cast<double>(ret.foundPerSet[setName]) / static_cast<double>(ret.hashedInputKmers.size()) >
				ret.bestFrac) {
			ret.bestFrac = static_cast<double>(ret.foundPerSet[setName]) / static_cast<double>(ret.hashedInputKmers.size());
			ret.winnerSet = setName;
			ret.winnerRevComp = false;
		}
		if (compPars.includeRevComp) {
			if (static_cast<double>(ret.foundPerSetRevComp[setName]) / static_cast<double>(ret.hashedInputKmers.size()) >
					ret.bestFrac) {
				ret.bestFrac = static_cast<double>(ret.foundPerSetRevComp[setName]) /
											 static_cast<double>(ret.hashedInputKmers.size());
				ret.winnerSet = setName;
				ret.winnerRevComp = true;
			}
		}
	}
//	if(print){
//		std::cout << "ret.winnerSet: " << ret.winnerSet << std::endl;
//		std::cout << std::endl;
//	}

	if (ret.winnerRevComp) {
		seq.reverseComplementRead(false, true);
	}
	return ret;
}




void UniqueKmerSetHelper::CompareReadToSetRes::zeroFillFoundSets(const VecStr &setNames) {
	for (const auto &setName: setNames) {
		foundPerSet[setName] = 0;
		foundPerSetRevComp[setName] = 0;
	}
}


uint32_t UniqueKmerSetHelper::CompareReadToSetRes::getTotalDetermined() {
	uint32_t count = 0;
	for (auto &perSet: foundPerSet) {
		count += perSet.second;
	}
	for (auto &perSet: foundPerSetRevComp) {
		count += perSet.second;
	}
	return count;
}



VecStr UniqueKmerSetHelper::CompareReadToSetRes::genOutputHeader(const CompareReadToSetPars &compPars) {
	VecStr header{"sample", "name", "totalUniqueKmers", "set", "totalUniqueInSet",
								"kmersFoundFromSet", "fracInSet", "fracOfSetFound"};
	if (compPars.includeRevComp) {
		addOtherVec(header, VecStr{"kmersFoundFromSetRevComp", "fracInSetRevComp",
															 "fracOfSetFoundRevComp"});
	}
	addOtherVec(header, VecStr{"winnerSet", "winnerSetFrac", "winnerSetRecComp"});
	return header;
}

void
UniqueKmerSetHelper::CompareReadToSetRes::writeOutputHeader(std::ostream &out, const CompareReadToSetPars &compPars,
																														const std::string &delim) {
	out << njh::conToStr(genOutputHeader(compPars), delim) << "\n";
}


std::vector<std::vector<std::string>> UniqueKmerSetHelper::CompareReadToSetRes::genOutput(const seqInfo &seq,
																																													const std::unordered_map<std::string, std::unordered_set<uint64_t>> &uniqueKmersPerSet,
																																													const CompareReadToSetPars &compPars) {
	std::vector<std::vector<std::string>> content;
	for (const auto &setName: uniqueKmersPerSet) {
		VecStr row = toVecStr(
						compPars.sampleName,
						seq.name_,
						hashedInputKmers.size(),
						setName.first,
						setName.second.size(),
						(njh::in(setName.first, foundPerSet) ? foundPerSet.at(setName.first) : 0),
						static_cast<double>(njh::in(setName.first, foundPerSet) ? foundPerSet.at(setName.first) : 0) /
						static_cast<double>(hashedInputKmers.size()),
						static_cast<double>(njh::in(setName.first, foundPerSet) ? foundPerSet.at(setName.first) : 0) /
						static_cast<double>(setName.second.size()));
		if (compPars.includeRevComp) {
			addOtherVec(row, toVecStr(
							(njh::in(setName.first, foundPerSetRevComp) ? foundPerSetRevComp.at(setName.first) : 0),
							static_cast<double>(njh::in(setName.first, foundPerSetRevComp) ? foundPerSetRevComp.at(setName.first)
																																						 : 0) /
							static_cast<double>(hashedInputKmersRevComp.size()),
							static_cast<double>(njh::in(setName.first, foundPerSetRevComp) ? foundPerSetRevComp.at(setName.first)
																																						 : 0) /
							static_cast<double>(setName.second.size()))

			);
		}
		addOtherVec(row, toVecStr(
						winnerSet,
						bestFrac,
						njh::boolToStr(winnerRevComp)));
		content.emplace_back(row);
	}
	return content;
}

void UniqueKmerSetHelper::CompareReadToSetRes::writeOutput(std::ostream &out, const seqInfo &seq,
																													 const std::unordered_map<std::string, std::unordered_set<uint64_t>> &uniqueKmersPerSet,
																													 const CompareReadToSetPars &compPars,
																													 const std::string &delim) {
	auto content = genOutput(seq, uniqueKmersPerSet, compPars);
	for (const auto &row: content) {
		out << njh::conToStr(row, delim) << "\n";
	}
}


uint32_t UniqueKmerSetHelper::getKmerLenFromUniqueKmerTable(const bfs::path &uniqueKmerTable) {
	auto toks = tokenizeString(njh::files::getFirstLine(uniqueKmerTable), "\t");
	return toks.at(1).size();
}


std::unordered_map<std::string, std::unordered_set<uint64_t>>
UniqueKmerSetHelper::readInUniqueKmerTablePerSet(const bfs::path &uniqueKmerTableFnp) {
	std::unordered_map<std::string, std::unordered_set<uint64_t>> ret;
	SimpleKmerHash hasher;
	TableReader uniqKmers(TableIOOpts::genTabFileIn(uniqueKmerTableFnp, false));
	if (uniqKmers.header_.nCol() < 2) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "need to have 2 columns" << "\n";
		throw std::runtime_error{ss.str()};
	}
	VecStr row;
	while (uniqKmers.getNextRow(row)) {
		ret[row[0]].emplace(hasher.hash(row[1]));
	}
	return ret;
}

std::unordered_set<uint64_t>
UniqueKmerSetHelper::readInUniqueKmerTableSetsCollapsed(const bfs::path &uniqueKmerTableFnp) {
	std::unordered_set<uint64_t> ret;
	SimpleKmerHash hasher;
	TableReader uniqKmers(TableIOOpts::genTabFileIn(uniqueKmerTableFnp, false));
	if (uniqKmers.header_.nCol() < 2) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "need to have 2 columns" << "\n";
		throw std::runtime_error{ss.str()};
	}
	VecStr row;
	while (uniqKmers.getNextRow(row)) {
		ret.emplace(hasher.hash(row[1]));
	}
	return ret;
}


void
UniqueKmerSetHelper::ProcessReadForExtractingCounts::addOtherCounts(const ProcessReadForExtractingCounts &otherCounts) {
	smallLenCutOffCount += otherCounts.smallLenCutOffCount;
	poorQualityCount += otherCounts.poorQualityCount;
	containsNs += otherCounts.containsNs;
	filteredDissimilarCount += otherCounts.filteredDissimilarCount;
	for (const auto &readsPerSetCount: otherCounts.readCountsPerSet.at(false)) {
		readCountsPerSet[false][readsPerSetCount.first] += readsPerSetCount.second;
	}
	for (const auto &readsPerSetRevCompCount: otherCounts.readCountsPerSet.at(true)) {
		readCountsPerSet[true][readsPerSetRevCompCount.first] += readsPerSetRevCompCount.second;
	}
}

VecStr UniqueKmerSetHelper::ProcessReadForExtractingCounts::genOutCountsHeader(
				const ProcessReadForExtractingPars &extractingPars) {
	VecStr header;
	if (extractingPars.compPars.includeRevComp) {
		header = toVecStr("iteration", "sample", "totalReads", "target", "count", "frac", "forwardCount", "fracForward");
	} else {
		header = toVecStr("iteration", "sample", "totalReads", "target", "count", "frac");
	}
	return header;
}


[[nodiscard]] std::vector<std::vector<std::string>> UniqueKmerSetHelper::ProcessReadForExtractingCounts::genOutCounts(const ProcessReadForExtractingPars & extractingPars,
																																 const std::unordered_map<std::string, std::unordered_set<uint64_t>> &uniqueKmersPerSet,
																																																											const std::string & iterName) const {
	std::vector<std::vector<std::string>> content;
	uint64_t totalReads = getTotalCounts();
	if (extractingPars.compPars.includeRevComp) {
		auto setNames = njh::getVecOfMapKeys(uniqueKmersPerSet);
		njh::sort(setNames);
		for(const auto & setName : setNames){
			uint64_t totalExtracted = (njh::in(setName, readCountsPerSet.at(false)) ? readCountsPerSet.at(false).at(setName) : 0) + (njh::in(setName, readCountsPerSet.at(true)) ? readCountsPerSet.at(true).at(setName) : 0);
			content.emplace_back(toVecStr(iterName
							, extractingPars.compPars.sampleName
							, totalReads
							, setName
							, totalExtracted
							, static_cast<double>(totalExtracted) / static_cast<double>(totalReads)
							, (njh::in(setName, readCountsPerSet.at(false)) ? readCountsPerSet.at(false).at(setName) : 0)
							, (totalExtracted > 0 ? static_cast<double>((njh::in(setName, readCountsPerSet.at(false)) ? readCountsPerSet.at(false).at(setName) : 0)) / static_cast<double>(totalExtracted) : 0)
			));
		}
		{
			uint64_t totalUndeterminedExtracted = (njh::in("undetermined", readCountsPerSet.at(false)) ? readCountsPerSet.at(false).at("undetermined") : 0) + (njh::in("undetermined", readCountsPerSet.at(true)) ? readCountsPerSet.at(true).at("undetermined") : 0);
			content.emplace_back(toVecStr( iterName
							, extractingPars.compPars.sampleName
							, totalReads
							, "undetermined"
							, totalUndeterminedExtracted
							, static_cast<double>(totalUndeterminedExtracted) / static_cast<double>(totalReads)
							, (njh::in("undetermined", readCountsPerSet.at(false)) ? readCountsPerSet.at(false).at("undetermined") : 0)
							, (totalUndeterminedExtracted > 0 ?static_cast<double>(readCountsPerSet.at(false).at("undetermined")) / static_cast<double>(totalUndeterminedExtracted) : 0)
			));
		}
		{
			content.emplace_back(toVecStr( iterName
							, extractingPars.compPars.sampleName
							, totalReads
							, "smallLenReads"
							, smallLenCutOffCount
							, static_cast<double>(smallLenCutOffCount) / static_cast<double>(totalReads)
							, smallLenCutOffCount
							, 0
			));

			content.emplace_back(toVecStr( iterName
							, extractingPars.compPars.sampleName
							, totalReads
							, "poorQualityCount"
							, poorQualityCount
							, static_cast<double>(poorQualityCount) / static_cast<double>(totalReads)
							, poorQualityCount
							, 0
			));

			content.emplace_back(toVecStr( iterName
							, extractingPars.compPars.sampleName
							, totalReads
							, "containsNs"
							, containsNs
							, static_cast<double>(containsNs) / static_cast<double>(totalReads)
							, containsNs
							, 0
			));




		}
	} else {
		auto setNames = njh::getVecOfMapKeys(uniqueKmersPerSet);
		njh::sort(setNames);
		for(const auto & setName : setNames){
			content.emplace_back(toVecStr( iterName
							, extractingPars.compPars.sampleName
							, totalReads
							, setName
							, (njh::in(setName, readCountsPerSet.at(false)) ? readCountsPerSet.at(false).at(setName) : 0)
							, static_cast<double>((njh::in(setName, readCountsPerSet.at(false)) ? readCountsPerSet.at(false).at(setName) : 0)) / static_cast<double>(totalReads)
			));
		}
		{
			content.emplace_back(toVecStr( iterName
							, extractingPars.compPars.sampleName
							, totalReads
							, "undetermined"
							, (njh::in("undetermined", readCountsPerSet.at(false)) ? readCountsPerSet.at(false).at("undetermined") : 0)
							, static_cast<double>((njh::in("undetermined", readCountsPerSet.at(false)) ? readCountsPerSet.at(false).at("undetermined") : 0)) / static_cast<double>(totalReads)
			));
		}
		{
			content.emplace_back(toVecStr( iterName
							, extractingPars.compPars.sampleName
							, totalReads
							, "smallLenReads"
							, smallLenCutOffCount
							, static_cast<double>(smallLenCutOffCount) / static_cast<double>(totalReads)
			));

			content.emplace_back(toVecStr( iterName
							, extractingPars.compPars.sampleName
							, totalReads
							, "poorQualityCount"
							, poorQualityCount
							, static_cast<double>(poorQualityCount) / static_cast<double>(totalReads)
			));

			content.emplace_back(toVecStr( iterName
							, extractingPars.compPars.sampleName
							, totalReads
							, "containsNs"
							, containsNs
							, static_cast<double>(containsNs) / static_cast<double>(totalReads)
			));
		}
	}
	return content;
}


void UniqueKmerSetHelper::ProcessReadForExtractingCounts::writeOutCounts(std::ostream &out,
																																				 const ProcessReadForExtractingPars &extractingPars,
																																				 const std::unordered_map<std::string, std::unordered_set<uint64_t>> &uniqueKmersPerSet,
																																				 const std::string &iterName,
																																				 const std::string &delim) const {
	auto content = genOutCounts(extractingPars, uniqueKmersPerSet, iterName);
	for (const auto &row: content) {
		out << njh::conToStr(row, delim) << std::endl;
	}
}




void UniqueKmerSetHelper::ProcessReadForExtractingCounts::writeOutCountsHeader(std::ostream &out,
																																							 const ProcessReadForExtractingPars &extractingPars,
																																							 const std::string &delim) {
	auto header = genOutCountsHeader(extractingPars);
	out << njh::conToStr(header, delim) << std::endl;
}


uint64_t UniqueKmerSetHelper::ProcessReadForExtractingCounts::getTotalCounts() const {
	uint64_t totalReads = poorQualityCount + containsNs + smallLenCutOffCount + genTotalUndeterminedCount() + genTotalDeterminedCount();
	return totalReads;
}

uint64_t UniqueKmerSetHelper::ProcessReadForExtractingCounts::genTotalUndeterminedCount() const {
	return (njh::in("undetermined", readCountsPerSet.at(false)) ? readCountsPerSet.at(false).at("undetermined") : 0) +
				 (njh::in("undetermined", readCountsPerSet.at(true)) ? readCountsPerSet.at(true).at("undetermined") : 0);
}

uint64_t UniqueKmerSetHelper::ProcessReadForExtractingCounts::genTotalDeterminedCount() const {
	uint64_t totalDeterminedCount = 0;
	for (const auto &readsPerSetCount: readCountsPerSet.at(false)) {
		if(readsPerSetCount.first != "undetermined"){
			totalDeterminedCount += readsPerSetCount.second;
		}
	}
	for (const auto &readsPerSetRevCompCount: readCountsPerSet.at(true)) {
		if(readsPerSetRevCompCount.first != "undetermined"){
			totalDeterminedCount += readsPerSetRevCompCount.second;
		}
	}
	return totalDeterminedCount;
}


void UniqueKmerSetHelper::processReadForExtracting(PairedRead &pseq,
																									 const std::unordered_map<std::string, std::unordered_set<uint64_t>> &uniqueKmersPerSet,
																									 const ProcessReadForExtractingPars &extractingPars,
																									 const SimpleKmerHash &hasher,
																									 MultiSeqIO &seqOut, ProcessReadForExtractingCounts &counts,
																									 const std::string & iterationName) {

	if(extractingPars.compPars.pairsSeparate){
		processReadForExtractingPairsSeparate(pseq, uniqueKmersPerSet, extractingPars, hasher, seqOut, counts, iterationName);
	}else{
		processReadForExtractingPairsTogether(pseq, uniqueKmersPerSet, extractingPars, hasher, seqOut, counts, iterationName);
	}
}

void UniqueKmerSetHelper::processReadForExtractingPairsSeparate(PairedRead &pseq,
																																const std::unordered_map<std::string, std::unordered_set<uint64_t>> &uniqueKmersPerSet,
																																const ProcessReadForExtractingPars &extractingPars,
																																const SimpleKmerHash &hasher,
																																MultiSeqIO &seqOut, ProcessReadForExtractingCounts &counts,
																																const std::string & iterationName) {

	if(len(pseq.seqBase_.seq_) < extractingPars.smallLenCutOff && len(pseq.mateSeqBase_.seq_) < extractingPars.smallLenCutOff){
		++counts.smallLenCutOffCount;
		return;
	}
	ReadCheckerQualCheck checker(extractingPars.qPars.qualCheck_, extractingPars.qPars.qualCheckCutOff_, false);

	if(extractingPars.qPars.checkingQFrac_ && !checker.checkRead(pseq)){
		++counts.poorQualityCount;
		return;
	}

	if(extractingPars.filterOnNs && (std::string::npos != pseq.seqBase_.seq_.find('N') || std::string::npos != pseq.mateSeqBase_.seq_.find('N')) ){
		++counts.containsNs;
		return;
	}

	MetaDataInName meta;
	if(extractingPars.markReadsPerIteration){
		if(MetaDataInName::nameHasMetaData(pseq.seqBase_.name_)){
			meta = MetaDataInName(pseq.seqBase_.name_);
		}
		meta.addMeta("iteration", iterationName, true);
		meta.resetMetaInName(pseq.seqBase_.name_);
		meta.resetMetaInName(pseq.mateSeqBase_.name_);
	}
	auto compResFirstMate = UniqueKmerSetHelper::compareReadToSets(pseq.seqBase_, uniqueKmersPerSet,
																																 extractingPars.compPars, hasher);

	auto compResSecondMate = UniqueKmerSetHelper::compareReadToSets(pseq.mateSeqBase_, uniqueKmersPerSet,
																																	extractingPars.compPars, hasher);

	if (compResFirstMate.winnerSet == compResSecondMate.winnerSet &&
			compResFirstMate.winnerRevComp == compResSecondMate.winnerRevComp &&
			(len(pseq.seqBase_.seq_) >= extractingPars.smallLenCutOff && len(pseq.mateSeqBase_.seq_) >= extractingPars.smallLenCutOff)) {
		++counts.readCountsPerSet[compResFirstMate.winnerRevComp][compResFirstMate.winnerSet];
		if (!extractingPars.doNotWriteUndetermined || compResFirstMate.winnerSet != "undetermined") {
			if (extractingPars.writeOutExclude || !njh::in(compResFirstMate.winnerSet, extractingPars.compPars.excludeSetNames)) {
				seqOut.openWrite(njh::pasteAsStr(compResFirstMate.winnerSet, "-paired"), pseq);
			}
		}
	} else {

		if (len(pseq.seqBase_.seq_) < extractingPars.smallLenCutOff) {
			++counts.smallLenCutOffCount;
		} else if (extractingPars.qPars.checkingQFrac_ && !checker.checkRead(pseq.seqBase_)) {
			++counts.poorQualityCount;
		} else if (extractingPars.filterOnNs && (std::string::npos != pseq.seqBase_.seq_.find('N'))) {
			++counts.containsNs;
		} else {
			pseq.seqBase_.name_.append("_firstMate");
			++counts.readCountsPerSet[compResFirstMate.winnerRevComp][compResFirstMate.winnerSet];

			if (!extractingPars.doNotWriteUndetermined || compResFirstMate.winnerSet != "undetermined") {
				if (extractingPars.writeOutExclude || !njh::in(compResFirstMate.winnerSet, extractingPars.compPars.excludeSetNames)) {
					seqOut.openWrite(njh::pasteAsStr(compResFirstMate.winnerSet, "-single"), pseq.seqBase_);
				}
			}
		}

		if (len(pseq.mateSeqBase_.seq_) < extractingPars.smallLenCutOff) {
			++counts.smallLenCutOffCount;
		} else if (extractingPars.qPars.checkingQFrac_ && !checker.checkRead(pseq.mateSeqBase_)) {
			++counts.poorQualityCount;
		} else if (extractingPars.filterOnNs && (std::string::npos != pseq.mateSeqBase_.seq_.find('N'))) {
			++counts.containsNs;
		} else {
			pseq.mateSeqBase_.name_.append("_secondMate");
			++counts.readCountsPerSet[compResSecondMate.winnerRevComp][compResSecondMate.winnerSet];
			if (!extractingPars.doNotWriteUndetermined || compResSecondMate.winnerSet != "undetermined") {
				if (extractingPars.writeOutExclude || !njh::in(compResSecondMate.winnerSet, extractingPars.compPars.excludeSetNames)) {
					seqOut.openWrite(njh::pasteAsStr(compResSecondMate.winnerSet, "-single"), pseq.mateSeqBase_);
				}
			}
		}
	}
}



void UniqueKmerSetHelper::processReadForFilteringPairsSeparate(PairedRead &pseq,
																									 const std::unordered_map<std::string, std::unordered_set<uint64_t>> &uniqueKmersPerSet,
																									 const ProcessReadForExtractingPars &extractingPars,
																									 const SimpleKmerHash &hasher,
																									 MultiSeqIO &seqOut,
																									 ProcessReadForExtractingCounts & counts,
																									 const std::string & iterationName,
																															 const std::string & filterName,
																															 const std::string & keepName) {

	MetaDataInName meta;
	if(extractingPars.markReadsPerIteration){
		if(MetaDataInName::nameHasMetaData(pseq.seqBase_.name_)){
			meta = MetaDataInName(pseq.seqBase_.name_);
		}
		meta.addMeta("iteration", iterationName, true);
		meta.resetMetaInName(pseq.seqBase_.name_);
		meta.resetMetaInName(pseq.mateSeqBase_.name_);
	}
	auto compResFirstMate = UniqueKmerSetHelper::compareReadToSets(pseq.seqBase_, uniqueKmersPerSet,
																																 extractingPars.compPars, hasher);

	auto compResSecondMate = UniqueKmerSetHelper::compareReadToSets(pseq.mateSeqBase_, uniqueKmersPerSet,
																																	extractingPars.compPars, hasher);

	auto firstMateTotalDetermined = compResFirstMate.getTotalDetermined();
	auto secondMateTotalDetermined = compResSecondMate.getTotalDetermined();
	if (firstMateTotalDetermined > 0 && secondMateTotalDetermined > 0) {
		seqOut.openWrite(njh::pasteAsStr(keepName, "-paired"), pseq);
	} else if (firstMateTotalDetermined == 0 && secondMateTotalDetermined == 0) {
		if (extractingPars.writeOutExclude) {
			seqOut.openWrite(njh::pasteAsStr(filterName, "-paired"), pseq);
		}
		counts.filteredDissimilarCount += 1;
	} else if (firstMateTotalDetermined == 0) {
		seqOut.openWrite(njh::pasteAsStr(keepName, "-single"), pseq.mateSeqBase_);
		counts.filteredDissimilarCount += 1;
		if (extractingPars.writeOutExclude) {
			seqOut.openWrite(njh::pasteAsStr(filterName, "-single"), pseq.seqBase_);
		}
	} else if (secondMateTotalDetermined == 0) {
		seqOut.openWrite(njh::pasteAsStr(keepName, "-single"), pseq.seqBase_);
		counts.filteredDissimilarCount += 1;
		if (extractingPars.writeOutExclude) {
			seqOut.openWrite(njh::pasteAsStr(filterName, "-single"), pseq.mateSeqBase_);
		}
	}
}



void UniqueKmerSetHelper::processReadForFilteringPairsTogether(PairedRead &pseq,
																															 const std::unordered_map<std::string, std::unordered_set<uint64_t>> &uniqueKmersPerSet,
																															 const ProcessReadForExtractingPars &extractingPars,
																															 const SimpleKmerHash &hasher,
																															 MultiSeqIO &seqOut,
																															 ProcessReadForExtractingCounts & counts,
																															 const std::string & iterationName,
																															 const std::string & filterName,
																															 const std::string & keepName) {

	MetaDataInName meta;
	if(extractingPars.markReadsPerIteration){
		if(MetaDataInName::nameHasMetaData(pseq.seqBase_.name_)){
			meta = MetaDataInName(pseq.seqBase_.name_);
		}
		meta.addMeta("iteration", iterationName, true);
		meta.resetMetaInName(pseq.seqBase_.name_);
		meta.resetMetaInName(pseq.mateSeqBase_.name_);
	}
	auto compResFirstMate = UniqueKmerSetHelper::compareReadToSets(pseq.seqBase_, uniqueKmersPerSet,
																																 extractingPars.compPars, hasher);

	auto compResSecondMate = UniqueKmerSetHelper::compareReadToSets(pseq.mateSeqBase_, uniqueKmersPerSet,
																																	extractingPars.compPars, hasher);

	auto firstMateTotalDetermined = compResFirstMate.getTotalDetermined();
	auto secondMateTotalDetermined = compResSecondMate.getTotalDetermined();
	if (firstMateTotalDetermined == 0 && secondMateTotalDetermined == 0) {
		if (extractingPars.writeOutExclude) {
			seqOut.openWrite(njh::pasteAsStr(filterName, "-paired"), pseq);
		}
		counts.filteredDissimilarCount += 1;
	} else {
		seqOut.openWrite(njh::pasteAsStr(keepName, "-paired"), pseq);
	}
}


void UniqueKmerSetHelper::processReadForFiltering
				(seqInfo &seq,
				 const std::unordered_map<std::string, std::unordered_set<uint64_t>> &uniqueKmersPerSet,
				 const ProcessReadForExtractingPars &extractingPars,
				 const SimpleKmerHash &hasher,
				 MultiSeqIO &seqOut,
				 ProcessReadForExtractingCounts &counts,
				 const std::string &iterationName,
				 const std::string &filterName,
				 const std::string &keepName) {

	MetaDataInName meta;
	if(extractingPars.markReadsPerIteration){
		if(MetaDataInName::nameHasMetaData(seq.name_)){
			meta = MetaDataInName(seq.name_);
		}
		meta.addMeta("iteration", iterationName, true);
		meta.resetMetaInName(seq.name_);
	}
	auto compRes = UniqueKmerSetHelper::compareReadToSets(seq, uniqueKmersPerSet, extractingPars.compPars, hasher);

	auto totalDetermined = compRes.getTotalDetermined();
	if (totalDetermined == 0) {
		if (extractingPars.writeOutExclude) {
			seqOut.openWrite(njh::pasteAsStr(filterName, "-single"), seq);
		}
		counts.filteredDissimilarCount += 1;
	} else {
		seqOut.openWrite(njh::pasteAsStr(keepName, "-single"), seq);
	}
}




void UniqueKmerSetHelper::processReadForExtractingPairsTogether(PairedRead &pseq,
																									 const std::unordered_map<std::string, std::unordered_set<uint64_t>> &uniqueKmersPerSet,
																									 const ProcessReadForExtractingPars &extractingPars,
																									 const SimpleKmerHash &hasher,
																									 MultiSeqIO &seqOut, ProcessReadForExtractingCounts &counts,
																									 const std::string & iterationName) {
	if(len(pseq.seqBase_.seq_) < extractingPars.smallLenCutOff && len(pseq.mateSeqBase_.seq_) < extractingPars.smallLenCutOff){
		++counts.smallLenCutOffCount;
		return;
	}

	ReadCheckerQualCheck checker(extractingPars.qPars.qualCheck_, extractingPars.qPars.qualCheckCutOff_, false);

	if(extractingPars.qPars.checkingQFrac_ && !checker.checkRead(pseq)){
		++counts.poorQualityCount;
		return;
	}

	if(extractingPars.filterOnNs && (std::string::npos != pseq.seqBase_.seq_.find('N') || std::string::npos != pseq.mateSeqBase_.seq_.find('N')) ){
		++counts.containsNs;
		return;
	}


	MetaDataInName meta;
	if(extractingPars.markReadsPerIteration){
		if(MetaDataInName::nameHasMetaData(pseq.seqBase_.name_)){
			meta = MetaDataInName(pseq.seqBase_.name_);
		}
		meta.addMeta("iteration", iterationName, true);
		meta.resetMetaInName(pseq.seqBase_.name_);
		meta.resetMetaInName(pseq.mateSeqBase_.name_);
	}

	auto compRes = UniqueKmerSetHelper::compareReadToSets(pseq, uniqueKmersPerSet, extractingPars.compPars, hasher);
	++counts.readCountsPerSet[compRes.winnerRevComp][compRes.winnerSet];
	if(!extractingPars.doNotWriteUndetermined || compRes.winnerSet != "undetermined"){
		if(extractingPars.writeOutExclude || !njh::in(compRes.winnerSet, extractingPars.compPars.excludeSetNames)){
			seqOut.openWrite(njh::pasteAsStr(compRes.winnerSet, "-paired"), pseq);
		}
	}
}

void UniqueKmerSetHelper::processReadForExtracting(seqInfo &seq,
																									 const std::unordered_map<std::string, std::unordered_set<uint64_t>> &uniqueKmersPerSet,
																									 const ProcessReadForExtractingPars & extractingPars, const SimpleKmerHash &hasher,
																									 MultiSeqIO & seqOut, ProcessReadForExtractingCounts & counts,
																									 const std::string & iterationName){
	if(len(seq.seq_) < extractingPars.smallLenCutOff){
		++counts.smallLenCutOffCount;
		return;
	}
	ReadCheckerQualCheck checker(extractingPars.qPars.qualCheck_, extractingPars.qPars.qualCheckCutOff_, false);

	if(extractingPars.qPars.checkingQFrac_ && !checker.checkRead(seq)){
		++counts.poorQualityCount;
		return;
	}

	if(extractingPars.filterOnNs && (std::string::npos != seq.seq_.find('N')) ){
		++counts.containsNs;
		return;
	}

	MetaDataInName meta;
	if(extractingPars.markReadsPerIteration){
		if(MetaDataInName::nameHasMetaData(seq.name_)){
			meta = MetaDataInName(seq.name_);
		}
		meta.addMeta("iteration", iterationName, true);
		meta.resetMetaInName(seq.name_);
	}
	auto compRes = UniqueKmerSetHelper::compareReadToSets(seq, uniqueKmersPerSet, extractingPars.compPars, hasher);
	++counts.readCountsPerSet[compRes.winnerRevComp][compRes.winnerSet];
	if(!extractingPars.doNotWriteUndetermined || compRes.winnerSet != "undetermined"){
		if(extractingPars.writeOutExclude || !njh::in(compRes.winnerSet, extractingPars.compPars.excludeSetNames)){
			seqOut.openWrite(njh::pasteAsStr(compRes.winnerSet, "-single"), seq);
		}
	}
}


std::unordered_map<std::string, std::unordered_map<uint64_t, uint32_t>>
UniqueKmerSetHelper::readInNewKmersFromExtractedReads(
				const bfs::path &directoryIn,
				const VecStr &kmerSets,
				const ProcessReadForExtractingPars &extractingPars) {
	std::unordered_map<std::string, std::unordered_map<uint64_t, uint32_t>> rawKmersPerInput;
	SimpleKmerHash hasher;

	for (const auto &kmerSetName: kmerSets) {
		{//paired in
			auto inOpts = SeqIOOptions::genPairedIn(
							njh::files::make_path(directoryIn, kmerSetName + "_R1.fastq"),
							njh::files::make_path(directoryIn, kmerSetName + "_R2.fastq"));
			inOpts.revComplMate_ = true;
			if (inOpts.inExists()) {
				PairedRead pseq;
				SeqInput extractedReader(inOpts);
				extractedReader.openIn();
				while (extractedReader.readNextRead(pseq)) {
					if(len(pseq.seqBase_.seq_) >= extractingPars.compPars.klen){
						for (uint32_t pos = 0; pos < len(pseq.seqBase_.seq_) - extractingPars.compPars.klen + 1; ++pos) {
							++rawKmersPerInput[kmerSetName][hasher.hash(pseq.seqBase_.seq_.substr(pos, extractingPars.compPars.klen))];
						}
					}
					if(len(pseq.mateSeqBase_.seq_) >= extractingPars.compPars.klen){
						for (uint32_t pos = 0; pos < len(pseq.mateSeqBase_.seq_) - extractingPars.compPars.klen + 1; ++pos) {
							++rawKmersPerInput[kmerSetName][hasher.hash(pseq.mateSeqBase_.seq_.substr(pos, extractingPars.compPars.klen))];
						}
					}
				}
			}
		}

		{//singles in
			auto inOpts = SeqIOOptions::genFastqIn(
							njh::files::make_path(directoryIn, kmerSetName + ".fastq"));
			if (inOpts.inExists()) {
				seqInfo seq;
				SeqInput extractedReader(inOpts);
				extractedReader.openIn();

				while (extractedReader.readNextRead(seq)) {
					if(len(seq.seq_) >= extractingPars.compPars.klen){
						for (uint32_t pos = 0; pos < len(seq.seq_) - extractingPars.compPars.klen + 1; ++pos) {
							++rawKmersPerInput[kmerSetName][hasher.hash(seq.seq_.substr(pos, extractingPars.compPars.klen))];
						}
					}
				}
			}
		}
	}
	return rawKmersPerInput;
}

std::unordered_map<std::string, std::unordered_map<uint64_t, uint32_t>>
UniqueKmerSetHelper::readInNewKmersFromExtractedReads(
				const bfs::path &directoryIn,
				const VecStr &kmerSets,
				const ProcessReadForExtractingPars &extractingPars,
				const std::unordered_map<std::string, UniqueKmerSetHelper::FilePositons> & positionAfterLastIteration) {
	std::unordered_map<std::string, std::unordered_map<uint64_t, uint32_t>> rawKmersPerInput;
	SimpleKmerHash hasher;

	for (const auto &kmerSetName: kmerSets) {
		{//paired in
			auto inOpts = SeqIOOptions::genPairedIn(
							njh::files::make_path(directoryIn, kmerSetName + "_R1.fastq"),
							njh::files::make_path(directoryIn, kmerSetName + "_R2.fastq"));
			inOpts.revComplMate_ = true;
			if (inOpts.inExists()) {
				bool moreToRead = true;
				if(njh::in(kmerSetName, positionAfterLastIteration)){
					if(njh::files::getFilePositionsEndOfFile(inOpts.firstName_) == positionAfterLastIteration.at(kmerSetName).r1FnpEnd){
						moreToRead = false;
					}
				}
				if(moreToRead){
					PairedRead pseq;
					SeqInput extractedReader(inOpts);
					extractedReader.openIn();
					if(njh::in(kmerSetName, positionAfterLastIteration)){
						extractedReader.seekgPri(positionAfterLastIteration.at(kmerSetName).r1FnpEnd);
						extractedReader.seekgSec(positionAfterLastIteration.at(kmerSetName).r2FnpEnd);
					}
					while (extractedReader.readNextRead(pseq)) {
						if(len(pseq.seqBase_.seq_) >= extractingPars.compPars.klen){
							for (uint32_t pos = 0; pos < len(pseq.seqBase_.seq_) - extractingPars.compPars.klen + 1; ++pos) {
								++rawKmersPerInput[kmerSetName][hasher.hash(pseq.seqBase_.seq_.substr(pos, extractingPars.compPars.klen))];
							}
						}
						if(len(pseq.mateSeqBase_.seq_) >= extractingPars.compPars.klen){
							for (uint32_t pos = 0; pos < len(pseq.mateSeqBase_.seq_) - extractingPars.compPars.klen + 1; ++pos) {
								++rawKmersPerInput[kmerSetName][hasher.hash(pseq.mateSeqBase_.seq_.substr(pos, extractingPars.compPars.klen))];
							}
						}
					}
				}
			}
		}

		{//singles in
			auto inOpts = SeqIOOptions::genFastqIn(
							njh::files::make_path(directoryIn, kmerSetName + ".fastq"));
			if (inOpts.inExists()) {
				bool moreToRead = true;
				if(njh::in(kmerSetName, positionAfterLastIteration)){
					if(njh::files::getFilePositionsEndOfFile(inOpts.firstName_) == positionAfterLastIteration.at(kmerSetName).singleFnpEnd){
						moreToRead = false;
					}
				}
				if(moreToRead){
					seqInfo seq;
					SeqInput extractedReader(inOpts);
					extractedReader.openIn();
					if(njh::in(kmerSetName, positionAfterLastIteration)){
						extractedReader.seekgPri(positionAfterLastIteration.at(kmerSetName).singleFnpEnd);
					}
					while (extractedReader.readNextRead(seq)) {
						if(len(seq.seq_) >= extractingPars.compPars.klen){
							for (uint32_t pos = 0; pos < len(seq.seq_) - extractingPars.compPars.klen + 1; ++pos) {
								++rawKmersPerInput[kmerSetName][hasher.hash(seq.seq_.substr(pos, extractingPars.compPars.klen))];
							}
						}
					}
				}
			}
		}
	}
	return rawKmersPerInput;
}




std::unordered_map<std::string, std::unordered_set<uint64_t>> UniqueKmerSetHelper::filterReExtractedKmersForNonUnique(
				const std::unordered_map<std::string, std::unordered_map<uint64_t, uint32_t>> &rawKmersPerInput,
				const ProcessReadForExtractingPars &extractingPars,
				const std::unordered_map<std::string, std::unordered_set<uint64_t>> &uniqueKmersPerSet,
				std::unordered_set<uint64_t> &nonUniqueKmersPerSet) {

	std::unordered_map<std::string, std::unordered_set<uint64_t>> outputUniqueKmersPerSet;
	SimpleKmerHash hasher;
	//have to adjust for the fact we are going to just use the unhashed num string;
	std::vector<char> newAllowableCharacters;
	newAllowableCharacters.reserve(extractingPars.compPars.allowableCharacters_.size());
	for (const auto c: extractingPars.compPars.allowableCharacters_) {
		newAllowableCharacters.push_back(hasher.hashBase(c));
	}
	std::function<bool(const std::string&)> hashedSeqCheck = [&newAllowableCharacters](const std::string & k){
		return std::all_of(k.begin(), k.end(), [&newAllowableCharacters](char base){return njh::in(base, newAllowableCharacters);});
	};

	for (const auto &finalNewKmerSet: rawKmersPerInput) {
		for (const auto finalNewKmer: finalNewKmerSet.second) {
			if (finalNewKmer.second < extractingPars.addingInKmersCountCutOff) {
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
				auto k = estd::to_string(finalNewKmer.first);
				if (hashedSeqCheck(k) && kmerInfo(k, extractingPars.compPars.kmerLengthForEntropyCalc_, false).computeKmerEntropy() > extractingPars.compPars.entropyFilter_) {
					outputUniqueKmersPerSet[finalNewKmerSet.first].emplace(finalNewKmer.first);
				}
			}
		}
	}
	for (const auto &set: uniqueKmersPerSet) {
		if(!extractingPars.doReCheckExcludeSets && njh::in(set.first, extractingPars.compPars.excludeSetNames)){
			continue;
		}
		for(const auto & finalKmer : set.second){
			bool pass = true;
			for(const auto & rawInput : rawKmersPerInput){
				if (set.first == rawInput.first) {
					//don't filter on self
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
	return outputUniqueKmersPerSet;
}


std::unordered_map<std::string, std::unordered_set<uint64_t>> UniqueKmerSetHelper::filterReExtractedKmersForNonUniqueIncludeExcludedSets(
				const std::unordered_map<std::string, std::unordered_map<uint64_t, uint32_t>> &rawKmersPerInput,
				const ProcessReadForExtractingPars &extractingPars,
				const std::unordered_map<std::string, std::unordered_set<uint64_t>> &uniqueKmersPerSet,
				std::unordered_set<uint64_t> &nonUniqueKmersPerSet) {
	std::unordered_map<std::string, std::unordered_set<uint64_t>> outputUniqueKmersPerSet;
	for (const auto &finalNewKmerSet: rawKmersPerInput) {
		for (const auto finalNewKmer: finalNewKmerSet.second) {
			if (finalNewKmer.second < extractingPars.addingInKmersCountCutOff) {
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

	for (const auto &set: uniqueKmersPerSet) {
		for(const auto & finalKmer : set.second){
			bool pass = true;
			for(const auto & rawInput : rawKmersPerInput){
				if (set.first == rawInput.first) {
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
	return outputUniqueKmersPerSet;
}



}  // namespace njhseq

