//
// Created by Nicholas Hathaway on 10/15/24.
//
#include <njhseq/readVectorManipulation/readVectorHelpers/readVecChecker.hpp>
#include "UniqueKmerSetHelper.hpp"


namespace njhseq {


namespace {

/**
 * @brief This thread's matcher for these sets, rebuilt when the sets or settings change
 */
const KmerCollectionReadMatcher & getMatcher(const std::unordered_map<std::string, std::unordered_set<uint64_t>> &uniqueKmersPerSet,
																						 const UniqueKmerSetHelper::CompareReadToSetPars &compPars) {
	static thread_local std::unordered_map<const void *, std::unique_ptr<KmerCollectionReadMatcher>> matchers;
	auto & matcher = matchers[&uniqueKmersPerSet];
	if (nullptr == matcher || !matcher->usable(uniqueKmersPerSet, compPars.klen, compPars.useInvertedKmerIndex_, false)) {
		matcher = std::make_unique<KmerCollectionReadMatcher>(uniqueKmersPerSet, compPars.klen, compPars.useInvertedKmerIndex_, false);
	}
	return *matcher;
}

/**
 * @brief Count the read's kmers in each set, zero sets under the cut offs (excluded sets have their own initial cut offs) and pick the winner
 */
UniqueKmerSetHelper::CompareReadToSetRes compareKmersToSets(const KmerCollectionReadMatcher &matcher,
																														KmerCollectionReadMatcher::Match &match,
																														const UniqueKmerSetHelper::CompareReadToSetPars &compPars) {
	matcher.countFound(match);
	const auto & names = matcher.names();
	const auto kmersPossible = static_cast<double>(match.kmers.size());
	for (uint32_t idx = 0; idx < names.size(); ++idx) {
		const bool excluded = njh::in(names[idx], compPars.excludeSetNames);
		const uint32_t hardCountOff = excluded ? compPars.initialExcludeHardCountOff : compPars.hardCountOff;
		const double fracCutOff = excluded ? compPars.initialExcludeFracCutOff : compPars.fracCutOff;
		for (auto * count : {&match.found[idx], &match.foundRevComp[idx]}) {
			if (*count < hardCountOff || *count / kmersPossible < fracCutOff) {
				*count = 0;
			}
		}
	}
	matcher.determineWinner(match);

	UniqueKmerSetHelper::CompareReadToSetRes ret;
	if (match.hasWinner()) {
		ret.winnerSet = names[match.winnerIdx];
		ret.winnerRevComp = match.winnerRevComp;
		ret.bestFrac = match.bestFrac;
	}
	for (uint32_t idx = 0; idx < names.size(); ++idx) {
		ret.foundPerSet[names[idx]] = match.found[idx];
		ret.foundPerSetRevComp[names[idx]] = match.foundRevComp[idx];
		if (match.found[idx] > 0 || match.foundRevComp[idx] > 0) {
			ret.allHits.emplace(names[idx]);
		}
	}
	ret.hashedInputKmers = std::move(match.kmers);
	ret.hashedInputKmersRevComp = std::move(match.kmersRevComp);
	return ret;
}

}  // namespace

UniqueKmerSetHelper::CompareReadToSetRes UniqueKmerSetHelper::compareReadToSets(PairedRead &pseq,
																																								const std::unordered_map<std::string, std::unordered_set<uint64_t>> &uniqueKmersPerSet,
																																								const CompareReadToSetPars &compPars,
																																								const SimpleKmerHash &hasher) {
	const auto & matcher = getMatcher(uniqueKmersPerSet, compPars);
	auto match = matcher.genMatch();
	//both mates' kmers are pooled
	matcher.addKmers(pseq.seqBase_.seq_, compPars.includeRevComp, match);
	matcher.addKmers(pseq.mateSeqBase_.seq_, compPars.includeRevComp, match);
	auto ret = compareKmersToSets(matcher, match, compPars);
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
	const auto & matcher = getMatcher(uniqueKmersPerSet, compPars);
	auto match = matcher.genMatch();
	matcher.addKmers(seq.seq_, compPars.includeRevComp, match);
	auto ret = compareKmersToSets(matcher, match, compPars);
	if (ret.winnerRevComp) {
		seq.reverseComplementRead(false, true);
	}
	return ret;
}




bool UniqueKmerSetHelper::readPassesQCChecks(seqInfo &seq,
																 const ProcessReadForExtractingPars &extractingPars,
																 ProcessReadForExtractingCounts &counts
	) {
	if(len(seq.seq_) < extractingPars.smallLenCutOff){
		++counts.smallLenCutOffCount;
		return false;
	}
	if(extractingPars.qPars.checkingQFrac_ && !extractingPars.qual_checker->checkRead(seq)){
		++counts.poorQualityCount;
		return false;
	}
	if(extractingPars.filterOnNs && std::string::npos != seq.seq_.find('N') ){
		++counts.containsNs;
		return false;
	}
	return true;
}

void UniqueKmerSetHelper::processReadForExtracting(seqInfo &seq,
																									 const std::unordered_map<std::string, std::unordered_set<uint64_t>> &uniqueKmersPerSet,
																									 const ProcessReadForExtractingPars & extractingPars, const SimpleKmerHash &hasher,
																									 MultiSeqIO & seqOut, ProcessReadForExtractingCounts & counts
																									 ){


	auto compRes = UniqueKmerSetHelper::compareReadToSets(seq, uniqueKmersPerSet, extractingPars.compPars, hasher);
	extractingPars.correctForSet(compRes.winnerSet, seq);
	// if (seq.name_.find("PfSN01_13_2727936_2728809_PfSN01_130074700_id48_pos433_size440_firstMate") != std::string::npos) {
	// 	std::cout << __FILE__ << " : " << __LINE__ << std::endl;
	// 	std::cout << "extractingPars.compPars.hardCountOff: " << extractingPars.compPars.hardCountOff << std::endl;
	// 	std::cout << "extractingPars.compPars.finalHardCountOff: " << extractingPars.compPars.finalHardCountOff << std::endl;
	// 	std::cout << "extractingPars.compPars.fracCutOff: " << extractingPars.compPars.fracCutOff << std::endl;
	// 	std::cout << "extractingPars.compPars.finalFracCutOff: " << extractingPars.compPars.finalFracCutOff << std::endl;
	// 	compRes.writeOutputHeader(std::cout, extractingPars.compPars);
	// 	compRes.writeOutput(std::cout, seq, uniqueKmersPerSet, extractingPars.compPars);
	// 	// auto compJson = njh::json::toJson(compRes);
	// 	// std::cout <<  compJson<< std::endl;
	// }



	if (extractingPars.filterMultiHitReads && compRes.allHits.size() > 1) {
		++counts.multiHitReadCounts[njh::conToStr(compRes.allHits, "::")];
		if(!extractingPars.doNotWriteUndetermined) {
			if (extractingPars.writeMultiHitSeparate) {
				seq.name_.append(njh::pasteAsStr("[hits=", njh::conToStr(compRes.allHits, "::"), ";]"));
				seqOut.openWrite(njh::pasteAsStr("multihit", "-single"), seq);
			} else {
				seqOut.openWrite(njh::pasteAsStr("undetermined", "-single"), seq);
			}
		}
	} else {
		++counts.readCountsPerSet[compRes.winnerRevComp][compRes.winnerSet];
		if(!extractingPars.doNotWriteUndetermined || compRes.winnerSet != "undetermined"){
			if(extractingPars.writeOutExclude || !njh::in(compRes.winnerSet, extractingPars.compPars.excludeSetNames)){
				seqOut.openWrite(njh::pasteAsStr(compRes.winnerSet, "-single"), seq);
			}
		}
	}
}

void UniqueKmerSetHelper::processReadForExtracting(PairedRead &pseq,
																									 const std::unordered_map<std::string, std::unordered_set<uint64_t>> &uniqueKmersPerSet,
																									 const ProcessReadForExtractingPars &extractingPars,
																									 const SimpleKmerHash &hasher,
																									 MultiSeqIO &seqOut, ProcessReadForExtractingCounts &counts
																									 ) {
	if(extractingPars.compPars.pairsSeparate){
		processReadForExtractingPairsSeparate(pseq, uniqueKmersPerSet, extractingPars, hasher, seqOut, counts);
	}else{
		processReadForExtractingPairsTogether(pseq, uniqueKmersPerSet, extractingPars, hasher, seqOut, counts);
	}
}

void UniqueKmerSetHelper::processReadForExtractingPairsTogether(PairedRead &pseq,
																									 const std::unordered_map<std::string, std::unordered_set<uint64_t>> &uniqueKmersPerSet,
																									 const ProcessReadForExtractingPars &extractingPars,
																									 const SimpleKmerHash &hasher,
																									 MultiSeqIO &seqOut, ProcessReadForExtractingCounts &counts
																									 ) {

	// auto compRes = UniqueKmerSetHelper::compareReadToSets(pseq, uniqueKmersPerSet, extractingPars.compPars, hasher);
	//
	// ++counts.readCountsPerSet[compRes.winnerRevComp][compRes.winnerSet];
	// ++counts.readCountsPerSet[compRes.winnerRevComp][compRes.winnerSet];
	// if(!extractingPars.doNotWriteUndetermined || compRes.winnerSet != "undetermined"){
	// 	if(extractingPars.writeOutExclude || !njh::in(compRes.winnerSet, extractingPars.compPars.excludeSetNames)){
	// 		seqOut.openWrite(njh::pasteAsStr(compRes.winnerSet, "-paired"), pseq);
	// 	}
	// }

	// auto compResFirstMate = UniqueKmerSetHelper::compareReadToSets(pseq.seqBase_, uniqueKmersPerSet,
	// 																															 extractingPars.compPars, hasher);
	//
	// auto compResSecondMate = UniqueKmerSetHelper::compareReadToSets(pseq.mateSeqBase_, uniqueKmersPerSet,
	// 																																extractingPars.compPars, hasher);
	//
	//
	// if (compResFirstMate.winnerSet == compResSecondMate.winnerSet &&
	// 		compResFirstMate.winnerRevComp == compResSecondMate.winnerRevComp ) {
	// 	++counts.readCountsPerSet[compResFirstMate.winnerRevComp][compResFirstMate.winnerSet];
	// 	++counts.readCountsPerSet[compResFirstMate.winnerRevComp][compResFirstMate.winnerSet];
	// 	if (!extractingPars.doNotWriteUndetermined || compResFirstMate.winnerSet != "undetermined") {
	// 		if (extractingPars.writeOutExclude || !njh::in(compResFirstMate.winnerSet, extractingPars.compPars.excludeSetNames)) {
	// 			seqOut.openWrite(njh::pasteAsStr(compResFirstMate.winnerSet, "-paired"), pseq);
	// 		}
	// 	}
	// } else {
	// 	if (compResFirstMate.winnerSet == compResSecondMate.winnerSet &&
	// 					 compResFirstMate.winnerRevComp != compResSecondMate.winnerRevComp &&
	// 					 !njh::in(compResFirstMate.winnerSet, extractingPars.compPars.excludeSetNames)) {
	// 		++counts.inversePairReadCountsPerSet[compResFirstMate.winnerSet];
	// 		++counts.inversePairReadCountsPerSet[compResFirstMate.winnerSet];
	// 	}
	// 	//mates didn't match into same category, leave in undetermined
	// 	++counts.readCountsPerSet[false]["undetermined"];
	// 	++counts.readCountsPerSet[false]["undetermined"];
	// 	if (!extractingPars.doNotWriteUndetermined) {
	// 		if(compResFirstMate.winnerRevComp != compResSecondMate.winnerRevComp) {
	// 			if(compResSecondMate.winnerRevComp) {
	// 				pseq.mateSeqBase_.reverseComplementRead(false, true);
	// 			} else {
	// 				pseq.seqBase_.reverseComplementRead(false, true);
	// 			}
	// 		}
	// 		seqOut.openWrite(njh::pasteAsStr("undetermined", "-paired"), pseq);
	// 	}
	// }

	auto compResFirstMate = UniqueKmerSetHelper::compareReadToSets(pseq.seqBase_, uniqueKmersPerSet,
																																 extractingPars.compPars, hasher);

	auto compResSecondMate = UniqueKmerSetHelper::compareReadToSets(pseq.mateSeqBase_, uniqueKmersPerSet,
																																	extractingPars.compPars, hasher);

	//both mates are now in their winning set's orientation, so strand-bias corrections (built by the
	//collection pass) can be applied before any write path sees them
	extractingPars.correctForSet(compResFirstMate.winnerSet, pseq.seqBase_);
	extractingPars.correctForSet(compResSecondMate.winnerSet, pseq.mateSeqBase_);

	if (compResFirstMate.winnerSet == compResSecondMate.winnerSet &&
						 compResFirstMate.winnerRevComp != compResSecondMate.winnerRevComp &&
						 !njh::in(compResFirstMate.winnerSet, extractingPars.compPars.excludeSetNames)) {
		++counts.inversePairReadCountsPerSet[compResFirstMate.winnerSet];
		++counts.inversePairReadCountsPerSet[compResFirstMate.winnerSet];
		//mates didn't match into same category, leave in undetermined
		++counts.readCountsPerSet[false]["undetermined"];
		++counts.readCountsPerSet[false]["undetermined"];
		if (!extractingPars.doNotWriteUndetermined) {
			// if(compResFirstMate.winnerRevComp != compResSecondMate.winnerRevComp) {
			// 	if(compResSecondMate.winnerRevComp) {
			// 		pseq.mateSeqBase_.reverseComplementRead(false, true);
			// 	} else {
			// 		pseq.seqBase_.reverseComplementRead(false, true);
			// 	}
			// }
			if(compResFirstMate.winnerRevComp) {
				pseq.seqBase_.reverseComplementRead(false, true);
			}
			if(compResSecondMate.winnerRevComp) {
				pseq.mateSeqBase_.reverseComplementRead(false, true);
			}
			seqOut.openWrite(njh::pasteAsStr("undetermined", "-paired"), pseq);
		}
	} else {
		if(compResFirstMate.winnerRevComp) {
			pseq.seqBase_.reverseComplementRead(false, true);
		}
		if(compResSecondMate.winnerRevComp) {
			pseq.mateSeqBase_.reverseComplementRead(false, true);
		}
		auto compRes = UniqueKmerSetHelper::compareReadToSets(pseq, uniqueKmersPerSet, extractingPars.compPars, hasher);
		if (extractingPars.filterMultiHitReads && compRes.allHits.size() > 1) {
			++counts.multiHitReadCounts[njh::conToStr(compRes.allHits, "::")];
			++counts.multiHitReadCounts[njh::conToStr(compRes.allHits, "::")];
			if(!extractingPars.doNotWriteUndetermined) {
				if (extractingPars.writeMultiHitSeparate) {
					pseq.seqBase_.name_.append(njh::pasteAsStr("[firstMateHits=", njh::conToStr(compResFirstMate.allHits, "::"), ";", "secondMateHits=", njh::conToStr(compResSecondMate.allHits, "::"), ";]"));
					pseq.mateSeqBase_.name_.append(njh::pasteAsStr("[firstMateHits=", njh::conToStr(compResFirstMate.allHits, "::"), ";", "secondMateHits=", njh::conToStr(compResSecondMate.allHits, "::"), ";]"));
					seqOut.openWrite(njh::pasteAsStr("multihit", "-paired"), pseq);
				} else {
					seqOut.openWrite(njh::pasteAsStr("undetermined", "-paired"), pseq);
				}
			}
		} else {
			++counts.readCountsPerSet[compRes.winnerRevComp][compRes.winnerSet];
			++counts.readCountsPerSet[compRes.winnerRevComp][compRes.winnerSet];
			if(!extractingPars.doNotWriteUndetermined || compRes.winnerSet != "undetermined"){
				if(extractingPars.writeOutExclude || !njh::in(compRes.winnerSet, extractingPars.compPars.excludeSetNames)){
					seqOut.openWrite(njh::pasteAsStr(compRes.winnerSet, "-paired"), pseq);
				}
			}
		}
	}
}

void UniqueKmerSetHelper::processReadForExtractingPairsSeparate(PairedRead &pseq,
																																const std::unordered_map<std::string, std::unordered_set<uint64_t>> &uniqueKmersPerSet,
																																const ProcessReadForExtractingPars &extractingPars,
																																const SimpleKmerHash &hasher,
																																MultiSeqIO &seqOut, ProcessReadForExtractingCounts &counts
																																) {
	// if (pseq.seqBase_.name_.find("PfSN01_080036600.1-exon2_id49_pos520_size170") != std::string::npos) {
	// 	std::cout << __FILE__ << " " << __LINE__ << std::endl;
	// 	std::cout << extractingPars.compPars.hardCountOff << std::endl;
	// 	std::cout << extractingPars.compPars.fracCutOff << std::endl;
	// }

	auto compResFirstMate = UniqueKmerSetHelper::compareReadToSets(pseq.seqBase_, uniqueKmersPerSet,
																																 extractingPars.compPars, hasher);

	auto compResSecondMate = UniqueKmerSetHelper::compareReadToSets(pseq.mateSeqBase_, uniqueKmersPerSet,
																																	extractingPars.compPars, hasher);

	//both mates are now in their winning set's orientation, so strand-bias corrections (built by the
	//collection pass) can be applied before any write path sees them
	extractingPars.correctForSet(compResFirstMate.winnerSet, pseq.seqBase_);
	extractingPars.correctForSet(compResSecondMate.winnerSet, pseq.mateSeqBase_);
	// if (pseq.seqBase_.name_.find("PfSN01_080036600.1-exon2_id49_pos520_size170") != std::string::npos) {
	// 	std::cout << "compResFirstMate: " << std::endl;
	// 	std::cout << njh::json::JsonConversion::toJson(compResFirstMate) << std::endl;
	// 	compResFirstMate.writeOutputHeader(std::cout, extractingPars.compPars);
	// 	compResFirstMate.writeOutput(std::cout, pseq.seqBase_, uniqueKmersPerSet, extractingPars.compPars);
	// 	std::cout << "compResSecondMate: " << std::endl;
	// 	std::cout << njh::json::JsonConversion::toJson(compResSecondMate) << std::endl;
	// 	compResSecondMate.writeOutputHeader(std::cout, extractingPars.compPars);
	// 	compResSecondMate.writeOutput(std::cout, pseq.mateSeqBase_, uniqueKmersPerSet, extractingPars.compPars);
	// }

	if (compResFirstMate.winnerSet == compResSecondMate.winnerSet &&
			compResFirstMate.winnerRevComp == compResSecondMate.winnerRevComp &&
			(!extractingPars.filterMultiHitReads || (compResFirstMate.allHits.size() == 1 && compResSecondMate.allHits.size() == 1) )) {
		++counts.readCountsPerSet[compResFirstMate.winnerRevComp][compResFirstMate.winnerSet];
		++counts.readCountsPerSet[compResFirstMate.winnerRevComp][compResFirstMate.winnerSet];
		if (!extractingPars.doNotWriteUndetermined || compResFirstMate.winnerSet != "undetermined") {
			if (extractingPars.writeOutExclude || !njh::in(compResFirstMate.winnerSet, extractingPars.compPars.excludeSetNames)) {
				seqOut.openWrite(njh::pasteAsStr(compResFirstMate.winnerSet, "-paired"), pseq);
			}
		}
	} else if (compResFirstMate.winnerSet == compResSecondMate.winnerSet &&
	           compResFirstMate.winnerRevComp != compResSecondMate.winnerRevComp &&
	           !njh::in(compResFirstMate.winnerSet, extractingPars.compPars.excludeSetNames)) {
		//one mate matched one direction of the input while the other mate matched the other direction, for now will put into "undetermined"
		++counts.readCountsPerSet[false]["undetermined"];
		++counts.readCountsPerSet[false]["undetermined"];
		++counts.inversePairReadCountsPerSet[compResFirstMate.winnerSet];
		++counts.inversePairReadCountsPerSet[compResFirstMate.winnerSet];
		if (!extractingPars.doNotWriteUndetermined) {
			if(compResSecondMate.winnerRevComp) {
				pseq.mateSeqBase_.reverseComplementRead(false, true);
			} else {
				pseq.seqBase_.reverseComplementRead(false, true);
			}
			seqOut.openWrite(njh::pasteAsStr("undetermined", "-paired"), pseq);
		}
	} else {
		pseq.seqBase_.name_.append("_firstMate");
		if (extractingPars.filterMultiHitReads && compResFirstMate.allHits.size() > 1) {
			++counts.multiHitReadCounts[njh::conToStr(compResFirstMate.allHits, "::")];
			if(!extractingPars.doNotWriteUndetermined) {
				if (extractingPars.writeMultiHitSeparate) {
					pseq.seqBase_.name_.append(njh::pasteAsStr("[hits=", njh::conToStr(compResFirstMate.allHits, "::"), ";]"));
					seqOut.openWrite(njh::pasteAsStr("multihit", "-single"), pseq.seqBase_);
				} else {
					seqOut.openWrite(njh::pasteAsStr("undetermined", "-single"), pseq.seqBase_);
				}
			}
		} else {
			++counts.readCountsPerSet[compResFirstMate.winnerRevComp][compResFirstMate.winnerSet];
			if (!extractingPars.doNotWriteUndetermined || compResFirstMate.winnerSet != "undetermined") {
				if (extractingPars.writeOutExclude || !njh::in(compResFirstMate.winnerSet, extractingPars.compPars.excludeSetNames)) {
					seqOut.openWrite(njh::pasteAsStr(compResFirstMate.winnerSet, "-single"), pseq.seqBase_);
				}
			}
		}

		pseq.mateSeqBase_.name_.append("_secondMate");
		if (extractingPars.filterMultiHitReads && compResSecondMate.allHits.size() > 1) {
			++counts.multiHitReadCounts[njh::conToStr(compResSecondMate.allHits, "::")];
			if(!extractingPars.doNotWriteUndetermined) {
				if (extractingPars.writeMultiHitSeparate) {
					pseq.mateSeqBase_.name_.append(njh::pasteAsStr("[hits=", njh::conToStr(compResSecondMate.allHits, "::"), ";]"));
					seqOut.openWrite(njh::pasteAsStr("multihit", "-single"), pseq.mateSeqBase_);
				} else {
					seqOut.openWrite(njh::pasteAsStr("undetermined", "-single"), pseq.mateSeqBase_);
				}
			}
		} else {
			++counts.readCountsPerSet[compResSecondMate.winnerRevComp][compResSecondMate.winnerSet];
			if (!extractingPars.doNotWriteUndetermined || compResSecondMate.winnerSet != "undetermined") {
				if (extractingPars.writeOutExclude || !njh::in(compResSecondMate.winnerSet, extractingPars.compPars.excludeSetNames)) {
					seqOut.openWrite(njh::pasteAsStr(compResSecondMate.winnerSet, "-single"), pseq.mateSeqBase_);
				}
			}
		}
	}
}





//filtering

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
	} else {
		//will always be true to reach here if secondMateTotalDetermined == 0
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




}  //namespace njhseq


