//
// Created by Nicholas Hathaway on 5/25/23.
//

#include "UniqueKmerSetHelper.hpp"


namespace njhseq {


UniqueKmerSetHelper::CompareReadToSetRes UniqueKmerSetHelper::compareReadToSetRes(PairedRead &pseq,
																																									const std::unordered_map<std::string, std::unordered_set<uint64_t>> &uniqueKmersPerSet,
																																									const CompareReadToSetPars &pars,
																																									const SimpleKmerHash &hasher) {
	CompareReadToSetRes ret;
	if (len(pseq.seqBase_.seq_) > pars.klen) {
		for (uint32_t pos = 0; pos < len(pseq.seqBase_.seq_) - pars.klen + 1; ++pos) {
			auto hash = hasher.hash(pseq.seqBase_.seq_.substr(pos, pars.klen));
			++ret.hashedInputKmers[hash];
		}
	}
	if (len(pseq.mateSeqBase_.seq_) > pars.klen) {
		for (uint32_t pos = 0; pos < len(pseq.mateSeqBase_.seq_) - pars.klen + 1; ++pos) {
			auto hash = hasher.hash(pseq.mateSeqBase_.seq_.substr(pos, pars.klen));
			++ret.hashedInputKmers[hash];
		}
	}
	if (pars.includeRevComp) {
		if (len(pseq.seqBase_.seq_) > pars.klen) {
			for (uint32_t pos = 0; pos < len(pseq.seqBase_.seq_) - pars.klen + 1; ++pos) {
				auto hash = hasher.revCompHash(pseq.seqBase_.seq_.substr(pos, pars.klen));
				++ret.hashedInputKmersRevComp[hash];
			}
		}
		if (len(pseq.mateSeqBase_.seq_) > pars.klen) {
			for (uint32_t pos = 0; pos < len(pseq.mateSeqBase_.seq_) - pars.klen + 1; ++pos) {
				auto hash = hasher.revCompHash(pseq.mateSeqBase_.seq_.substr(pos, pars.klen));
				++ret.hashedInputKmersRevComp[hash];
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

	if (pars.includeRevComp) {
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
		if (perSet.second < pars.hardCountOff) {
			perSet.second = 0;
		}
	}
	for (auto &perSet: ret.foundPerSetRevComp) {
		if (perSet.second < pars.hardCountOff) {
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
		if (pars.includeRevComp) {
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

UniqueKmerSetHelper::CompareReadToSetRes UniqueKmerSetHelper::compareReadToSetRes(seqInfo &seq,
																																									const std::unordered_map<std::string, std::unordered_set<uint64_t>> &uniqueKmersPerSet,
																																									const CompareReadToSetPars &pars,
																																									const SimpleKmerHash &hasher) {
	CompareReadToSetRes ret;
	if (len(seq.seq_) > pars.klen) {
		for (uint32_t pos = 0; pos < len(seq.seq_) - pars.klen + 1; ++pos) {
			auto hash = hasher.hash(seq.seq_.substr(pos, pars.klen));
			++ret.hashedInputKmers[hash];
		}
	}
	if (pars.includeRevComp) {
		if (len(seq.seq_) > pars.klen) {
			for (uint32_t pos = 0; pos < len(seq.seq_) - pars.klen + 1; ++pos) {
				auto hash = hasher.revCompHash(seq.seq_.substr(pos, pars.klen));
				++ret.hashedInputKmersRevComp[hash];
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

	if (pars.includeRevComp) {
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
		if (perSet.second < pars.hardCountOff) {
			perSet.second = 0;
		}
	}
	for (auto &perSet: ret.foundPerSetRevComp) {
		if (perSet.second < pars.hardCountOff) {
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
		if (pars.includeRevComp) {
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
	for (const auto &readsPerSetCount: otherCounts.readsPerSet) {
		readsPerSet[readsPerSetCount.first] += readsPerSetCount.second;
	}
	for (const auto &readsPerSetRevCompCount: otherCounts.readsPerSetRevComp) {
		readsPerSetRevComp[readsPerSetRevCompCount.first] += readsPerSetRevCompCount.second;
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
			uint64_t totalExtracted = (njh::in(setName, readsPerSet) ? readsPerSet.at(setName) : 0) + (njh::in(setName, readsPerSetRevComp) ? readsPerSetRevComp.at(setName) : 0);
			content.emplace_back(toVecStr(iterName
							, extractingPars.compPars.sampleName
							, totalReads
							, setName
							, totalExtracted
							, static_cast<double>(totalExtracted) / static_cast<double>(totalReads)
							, (njh::in(setName, readsPerSet) ? readsPerSet.at(setName) : 0)
							, (totalExtracted > 0 ? static_cast<double>((njh::in(setName, readsPerSet) ? readsPerSet.at(setName) : 0)) / static_cast<double>(totalExtracted) : 0)
			));
		}
		{
			uint64_t totalUndeterminedExtracted = (njh::in("undetermined", readsPerSet) ? readsPerSet.at("undetermined") : 0) + (njh::in("undetermined", readsPerSetRevComp) ? readsPerSetRevComp.at("undetermined") : 0);
			content.emplace_back(toVecStr( iterName
							, extractingPars.compPars.sampleName
							, totalReads
							, "undetermined"
							, totalUndeterminedExtracted
							, static_cast<double>(totalUndeterminedExtracted) / static_cast<double>(totalReads)
							, (njh::in("undetermined", readsPerSet) ? readsPerSet.at("undetermined") : 0)
							, (totalUndeterminedExtracted > 0 ?static_cast<double>(readsPerSet.at("undetermined")) / static_cast<double>(totalUndeterminedExtracted) : 0)
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
		}
	} else {
		auto setNames = njh::getVecOfMapKeys(uniqueKmersPerSet);
		njh::sort(setNames);
		for(const auto & setName : setNames){
			content.emplace_back(toVecStr( iterName
							, extractingPars.compPars.sampleName
							, totalReads
							, setName
							, (njh::in(setName, readsPerSet) ? readsPerSet.at(setName) : 0)
							, static_cast<double>((njh::in(setName, readsPerSet) ? readsPerSet.at(setName) : 0)) / static_cast<double>(totalReads)
			));
		}
		{
			content.emplace_back(toVecStr( iterName
							, extractingPars.compPars.sampleName
							, totalReads
							, "undetermined"
							, (njh::in("undetermined", readsPerSet) ? readsPerSet.at("undetermined") : 0)
							, static_cast<double>((njh::in("undetermined", readsPerSet) ? readsPerSet.at("undetermined") : 0)) / static_cast<double>(totalReads)
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
		out << njh::conToStr(row) << std::endl;
	}
}




void UniqueKmerSetHelper::ProcessReadForExtractingCounts::writeOutCountsHeader(std::ostream &out,
																																							 const ProcessReadForExtractingPars &extractingPars,
																																							 const std::string &delim) {
	auto header = genOutCountsHeader(extractingPars);
	out << njh::conToStr(header, delim) << std::endl;
}


uint64_t UniqueKmerSetHelper::ProcessReadForExtractingCounts::getTotalCounts() const {
	uint64_t totalReads = smallLenCutOffCount + genTotalUndeterminedCount() + genTotalDeterminedCount();
	return totalReads;
}

uint64_t UniqueKmerSetHelper::ProcessReadForExtractingCounts::genTotalUndeterminedCount() const {
	return (njh::in("undetermined", readsPerSet) ? readsPerSet.at("undetermined") : 0) +
				 (njh::in("undetermined", readsPerSetRevComp) ? readsPerSetRevComp.at("undetermined") : 0);
}

uint64_t UniqueKmerSetHelper::ProcessReadForExtractingCounts::genTotalDeterminedCount() const {
	uint64_t totalDeterminedCount = 0;
	for (const auto &readsPerSetCount: readsPerSet) {
		if(readsPerSetCount.first != "undetermined"){
			totalDeterminedCount += readsPerSetCount.second;
		}
	}
	for (const auto &readsPerSetRevCompCount: readsPerSetRevComp) {
		if(readsPerSetRevCompCount.first != "undetermined"){
			totalDeterminedCount += readsPerSetRevCompCount.second;
		}
	}
	return totalDeterminedCount;
}



void UniqueKmerSetHelper::processReadForExtracting(PairedRead &pseq,
															const std::unordered_map<std::string, std::unordered_set<uint64_t>> &uniqueKmersPerSet,
															const ProcessReadForExtractingPars & extractingPars, const SimpleKmerHash &hasher,
															MultiSeqIO & seqOut, ProcessReadForExtractingCounts & counts){
//	std::cout << __FILE__ << " " << __LINE__ << std::endl;
	if(len(pseq.seqBase_.seq_) < extractingPars.smallLenCutOff && len(pseq.mateSeqBase_.seq_) < extractingPars.smallLenCutOff){
		++counts.smallLenCutOffCount;
		return;
	}
//	std::cout << __FILE__ << " " << __LINE__ << std::endl;

	if(extractingPars.compPars.pairsSeparate){
//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		auto compResFirstMate = UniqueKmerSetHelper::compareReadToSetRes(pseq.seqBase_, uniqueKmersPerSet, extractingPars.compPars,
																																		 hasher);
//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		auto compResSecondMate = UniqueKmerSetHelper::compareReadToSetRes(pseq.mateSeqBase_, uniqueKmersPerSet,
																																			extractingPars.compPars, hasher);
//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		if (compResFirstMate.winnerSet == compResSecondMate.winnerSet &&
				compResFirstMate.winnerRevComp == compResSecondMate.winnerRevComp &&
				(len(pseq.seqBase_.seq_) >= extractingPars.smallLenCutOff && len(pseq.mateSeqBase_.seq_) >= extractingPars.smallLenCutOff)) {
//			std::cout << __FILE__ << " " << __LINE__ << std::endl;
			if (compResFirstMate.winnerRevComp) {
				++counts.readsPerSetRevComp[compResFirstMate.winnerSet];
			} else {
				++counts.readsPerSet[compResFirstMate.winnerSet];
			}
//			std::cout << __FILE__ << " " << __LINE__ << std::endl;
			if (!extractingPars.doNotWriteUndetermined || compResFirstMate.winnerSet != "undetermined") {
				if (extractingPars.writeOutExclude || !njh::in(compResFirstMate.winnerSet, extractingPars.excludeSetNames)) {
//					std::cout << __FILE__ << " " << __LINE__ << std::endl;
//					std::cout << "compResFirstMate.winnerSet: " << compResFirstMate.winnerSet << std::endl;
//					std::cout << "njh::pasteAsStr(compResFirstMate.winnerSet, \"-paired\"): " << njh::pasteAsStr(compResFirstMate.winnerSet, "-paired") << std::endl;
//					std::cout << "pseq.seqBase_.name_ " << pseq.seqBase_.name_ << std::endl;
					seqOut.openWrite(njh::pasteAsStr(compResFirstMate.winnerSet, "-paired"), pseq);
//					std::cout << __FILE__ << " " << __LINE__ << std::endl;
				}
			}
//			std::cout << __FILE__ << " " << __LINE__ << std::endl;
		} else {
//			std::cout << __FILE__ << " " << __LINE__ << std::endl;
			if (len(pseq.seqBase_.seq_) < extractingPars.smallLenCutOff) {
				++counts.smallLenCutOffCount;
			} else {
//				std::cout << __FILE__ << " " << __LINE__ << std::endl;
				pseq.seqBase_.name_.append("_firstMate");
				if (compResFirstMate.winnerRevComp) {
					++counts.readsPerSetRevComp[compResFirstMate.winnerSet];
				} else {
					++counts.readsPerSet[compResFirstMate.winnerSet];
				}
//				std::cout << __FILE__ << " " << __LINE__ << std::endl;
				if (!extractingPars.doNotWriteUndetermined || compResFirstMate.winnerSet != "undetermined") {
					if (extractingPars.writeOutExclude || !njh::in(compResFirstMate.winnerSet, extractingPars.excludeSetNames)) {
//						std::cout << __FILE__ << " " << __LINE__ << std::endl;
						seqOut.openWrite(njh::pasteAsStr(compResFirstMate.winnerSet, "-single"), pseq.seqBase_);
//						std::cout << __FILE__ << " " << __LINE__ << std::endl;
					}
				}
			}
//			std::cout << __FILE__ << " " << __LINE__ << std::endl;
			if (len(pseq.mateSeqBase_.seq_) < extractingPars.smallLenCutOff) {
				++counts.smallLenCutOffCount;
			} else {
//				std::cout << __FILE__ << " " << __LINE__ << std::endl;
				pseq.mateSeqBase_.name_.append("_secondMate");
				if (compResSecondMate.winnerRevComp) {
					++counts.readsPerSetRevComp[compResSecondMate.winnerSet];
				} else {
					++counts.readsPerSet[compResSecondMate.winnerSet];
				}
//				std::cout << __FILE__ << " " << __LINE__ << std::endl;
				if (!extractingPars.doNotWriteUndetermined || compResSecondMate.winnerSet != "undetermined") {
					if (extractingPars.writeOutExclude || !njh::in(compResSecondMate.winnerSet, extractingPars.excludeSetNames)) {
//						std::cout << __FILE__ << " " << __LINE__ << std::endl;
						pseq.mateSeqBase_.reverseComplementRead(false, true);
//						std::cout << __FILE__ << " " << __LINE__ << std::endl;
						seqOut.openWrite(njh::pasteAsStr(compResSecondMate.winnerSet, "-single"), pseq.mateSeqBase_);
//						std::cout << __FILE__ << " " << __LINE__ << std::endl;
					}
				}
			}
//			std::cout << __FILE__ << " " << __LINE__ << std::endl;
		}
//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
	} else {
		auto compRes = UniqueKmerSetHelper::compareReadToSetRes(pseq, uniqueKmersPerSet, extractingPars.compPars, hasher);
		if (compRes.winnerRevComp) {
			++counts.readsPerSetRevComp[compRes.winnerSet];
		} else {
			++counts.readsPerSet[compRes.winnerSet];
		}
		if(!extractingPars.doNotWriteUndetermined || compRes.winnerSet != "undetermined"){
			if(extractingPars.writeOutExclude || !njh::in(compRes.winnerSet, extractingPars.excludeSetNames)){
				seqOut.openWrite(njh::pasteAsStr(compRes.winnerSet, "-paired"), pseq);
			}
		}
	}
}

void UniqueKmerSetHelper::processReadForExtracting(seqInfo &seq,
																									 const std::unordered_map<std::string, std::unordered_set<uint64_t>> &uniqueKmersPerSet,
																									 const ProcessReadForExtractingPars & extractingPars, const SimpleKmerHash &hasher,
																									 MultiSeqIO & seqOut, ProcessReadForExtractingCounts & counts){
	if(len(seq.seq_) < extractingPars.smallLenCutOff){
		++counts.smallLenCutOffCount;
		return;
	}

	auto compRes = UniqueKmerSetHelper::compareReadToSetRes(seq, uniqueKmersPerSet, extractingPars.compPars, hasher);
	if (compRes.winnerRevComp) {
		++counts.readsPerSetRevComp[compRes.winnerSet];
	} else {
		++counts.readsPerSet[compRes.winnerSet];
	}
	if(!extractingPars.doNotWriteUndetermined || compRes.winnerSet != "undetermined"){
		if(extractingPars.writeOutExclude || !njh::in(compRes.winnerSet, extractingPars.excludeSetNames)){
			seqOut.openWrite(njh::pasteAsStr(compRes.winnerSet, "-single"), seq);
		}
	}
}


std::unordered_map<std::string, std::unordered_map<uint64_t, uint32_t>> UniqueKmerSetHelper::readInNewKmersFromExtractedReads(const bfs::path & directoryIn,
																																																				 const VecStr & kmerSets,
																																																				 const ProcessReadForExtractingPars & extractingPars){
	std::unordered_map<std::string, std::unordered_map<uint64_t, uint32_t>> rawKmersPerInput;
	SimpleKmerHash hasher;

	for (const auto &kmerSetName: kmerSets) {
		{//paired in
			auto inOpts = SeqIOOptions::genPairedInGz(
							njh::files::make_path(directoryIn, kmerSetName + "_R1.fastq.gz"),
							njh::files::make_path(directoryIn, kmerSetName + "_R2.fastq.gz"));
			inOpts.revComplMate_ = true;
			if (inOpts.inExists()) {
				PairedRead pseq;
				SeqInput extractedReader(inOpts);
				extractedReader.openIn();
				while (extractedReader.readNextRead(pseq)) {
					for (uint32_t pos = 0; pos < len(pseq.seqBase_.seq_) - extractingPars.compPars.klen + 1; ++pos) {
						++rawKmersPerInput[kmerSetName][hasher.hash(pseq.seqBase_.seq_.substr(pos, extractingPars.compPars.klen))];
					}
					for (uint32_t pos = 0; pos < len(pseq.mateSeqBase_.seq_) - extractingPars.compPars.klen + 1; ++pos) {
						++rawKmersPerInput[kmerSetName][hasher.hash(
										pseq.mateSeqBase_.seq_.substr(pos, extractingPars.compPars.klen))];
					}
				}
			}
		}

		{//singles in
			auto inOpts = SeqIOOptions::genFastqInGz(
							njh::files::make_path(directoryIn, kmerSetName + ".fastq.gz"));
			if (inOpts.inExists()) {
				seqInfo seq;
				SeqInput extractedReader(inOpts);
				extractedReader.openIn();
				while (extractedReader.readNextRead(seq)) {
					for (uint32_t pos = 0; pos < len(seq.seq_) - extractingPars.compPars.klen + 1; ++pos) {
						++rawKmersPerInput[kmerSetName][hasher.hash(seq.seq_.substr(pos, extractingPars.compPars.klen))];
					}
				}
			}
		}
	}
	return rawKmersPerInput;
};


std::unordered_map<std::string, std::unordered_set<uint64_t>> UniqueKmerSetHelper::filterReExtractedKmersForNonUnique(
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
		if(!extractingPars.doReCheckExcludeSets && njh::in(set.first, extractingPars.excludeSetNames)){
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
	return outputUniqueKmersPerSet;
}



}  // namespace njhseq

