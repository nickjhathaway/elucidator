//
// Created by Nicholas Hathaway on 10/15/24.
//
#include <njhseq/readVectorManipulation/readVectorHelpers/readVecChecker.hpp>
#include "UniqueKmerSetHelper.hpp"


namespace njhseq {


UniqueKmerSetHelper::CompareReadToSetRes UniqueKmerSetHelper::compareReadToSets(PairedRead &pseq,
																																								const std::unordered_map<std::string, std::unordered_set<uint64_t>> &uniqueKmersPerSet,
																																								const CompareReadToSetPars &compPars,
																																								const SimpleKmerHash &hasher) {
	CompareReadToSetRes ret;
	// hash kmers
	if(compPars.includeRevComp) {
		ret.increaseHashedInputKmerCountsBoth(pseq.seqBase_, compPars, hasher);
		ret.increaseHashedInputKmerCountsBoth(pseq.mateSeqBase_, compPars, hasher);
	} else {
		ret.increaseHashedInputKmerCounts(pseq.seqBase_, compPars, hasher);
		ret.increaseHashedInputKmerCounts(pseq.mateSeqBase_, compPars, hasher);
	}

	// check for presence of kmers in sets
	ret.checkForPresenceOfHashedInputKmers(uniqueKmersPerSet);

	uint32_t kmersPossible = (pseq.seqBase_.seq_.size() - compPars.klen + 1) + (pseq.mateSeqBase_.seq_.size() - compPars.klen + 1);

	ret.preprocessedFoundKmers(compPars, kmersPossible);



	//set a hard cut off or the fraction of possible kmers is less than frac cut off, reset counts to zero
	// for (auto &perSet: ret.foundPerSet) {
	// 	if (perSet.second < compPars.hardCountOff || perSet.second / static_cast<double>(kmersPossible) < compPars.fracCutOff) {
	// 		perSet.second = 0;
	// 	}
	// }
	// for (auto &perSet: ret.foundPerSetRevComp) {
	// 	if (perSet.second < compPars.hardCountOff || perSet.second / static_cast<double>(kmersPossible) < compPars.fracCutOff) {
	// 		perSet.second = 0;
	// 	}
	// }

	ret.determineWinner(uniqueKmersPerSet);

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

	// hash kmers
	if(compPars.includeRevComp) {
		ret.increaseHashedInputKmerCountsBoth(seq, compPars, hasher);
	} else {
		ret.increaseHashedInputKmerCounts(seq, compPars, hasher);
	}

  // check for presence of kmers in sets
	ret.checkForPresenceOfHashedInputKmers(uniqueKmersPerSet);

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
	ret.preprocessedFoundKmers(compPars, kmersPossible);
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

	ret.determineWinner(uniqueKmersPerSet);

//	if(print){
//		std::cout << "ret.winnerSet: " << ret.winnerSet << std::endl;
//		std::cout << std::endl;
//	}

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

	++counts.readCountsPerSet[compRes.winnerRevComp][compRes.winnerSet];

	if(!extractingPars.doNotWriteUndetermined || compRes.winnerSet != "undetermined"){
		if(extractingPars.writeOutExclude || !njh::in(compRes.winnerSet, extractingPars.compPars.excludeSetNames)){
			seqOut.openWrite(njh::pasteAsStr(compRes.winnerSet, "-single"), seq);
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

		++counts.readCountsPerSet[compRes.winnerRevComp][compRes.winnerSet];
		++counts.readCountsPerSet[compRes.winnerRevComp][compRes.winnerSet];
		if(!extractingPars.doNotWriteUndetermined || compRes.winnerSet != "undetermined"){
			if(extractingPars.writeOutExclude || !njh::in(compRes.winnerSet, extractingPars.compPars.excludeSetNames)){
				seqOut.openWrite(njh::pasteAsStr(compRes.winnerSet, "-paired"), pseq);
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
	auto compResFirstMate = UniqueKmerSetHelper::compareReadToSets(pseq.seqBase_, uniqueKmersPerSet,
																																 extractingPars.compPars, hasher);

	auto compResSecondMate = UniqueKmerSetHelper::compareReadToSets(pseq.mateSeqBase_, uniqueKmersPerSet,
																																	extractingPars.compPars, hasher);
	if (compResFirstMate.winnerSet == compResSecondMate.winnerSet &&
			compResFirstMate.winnerRevComp == compResSecondMate.winnerRevComp ) {
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
		++counts.readCountsPerSet[compResFirstMate.winnerRevComp][compResFirstMate.winnerSet];
		if (!extractingPars.doNotWriteUndetermined || compResFirstMate.winnerSet != "undetermined") {
			if (extractingPars.writeOutExclude || !njh::in(compResFirstMate.winnerSet, extractingPars.compPars.excludeSetNames)) {
				seqOut.openWrite(njh::pasteAsStr(compResFirstMate.winnerSet, "-single"), pseq.seqBase_);
			}
		}
		pseq.mateSeqBase_.name_.append("_secondMate");
		++counts.readCountsPerSet[compResSecondMate.winnerRevComp][compResSecondMate.winnerSet];
		if (!extractingPars.doNotWriteUndetermined || compResSecondMate.winnerSet != "undetermined") {
			if (extractingPars.writeOutExclude || !njh::in(compResSecondMate.winnerSet, extractingPars.compPars.excludeSetNames)) {
				seqOut.openWrite(njh::pasteAsStr(compResSecondMate.winnerSet, "-single"), pseq.mateSeqBase_);
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


