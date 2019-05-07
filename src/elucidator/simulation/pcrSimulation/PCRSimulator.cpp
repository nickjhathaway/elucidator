/*
 * PCRSimulator.cpp
 *
 *  Created on: May 2, 2019
 *      Author: nicholashathaway
 */


#include "PCRSimulator.hpp"

namespace njhseq {


PCRSimulator::PCRSimulator(uint64_t intErrorRate) :
		intErrorRate_(intErrorRate) {

}


std::unordered_map<uint32_t, uint32_t> PCRSimulator::genMutCounts(
		std::mt19937_64 & gen,
		const uint64_t seqNumber,
		const uint64_t seqSize) const {
	std::unordered_map<uint32_t, uint32_t> currentMuts;
	for (uint32_t read = 0; read < seqNumber; ++read) {
		uint32_t count = 0;
		for (uint32_t base = 0; base < seqSize; ++base) {
			if (gen() <= intErrorRate_) {
				++count;
			}
		}
		if (count != 0) {
			++currentMuts[count];
		}
	}
	return currentMuts;
}

std::string PCRSimulator::mutateSeq(std::string seq,
		uint32_t mutNum,
		const std::vector<uint64_t>& readPositions,
		njh::randomGenerator & gen,
		std::unordered_map<char, njh::randObjectGen<char, uint32_t>>& charGen) {
	auto seqPositons = gen.unifRandSelectionVec(readPositions, mutNum, false);
	for (auto seqPos : seqPositons) {
		char base = charGen.at(seq[seqPos]).genObj();
		seq[seqPos] = base;
	}
	return seq;
}


std::vector<uint64_t> PCRSimulator::spreadReadNumAcrossThreads(const uint64_t readTotal,
		const uint32_t numThreads) {
	std::vector<uint64_t> tempAmounts;
	uint64_t tempAmount = readTotal / numThreads;
	uint64_t sum = 0;
	for (uint32_t t = 0; t < numThreads; ++t) {
		sum += tempAmount;
		tempAmounts.emplace_back(tempAmount);
	}
	tempAmounts.back() += (readTotal % sum);
	return tempAmounts;
}


PCRSimulator::PCRProduct::PCRProduct(const std::string & pcrSeq,
		const uint32_t pcrRoundsLeft, const uint32_t templateAmount) :
		pcrSeq_(pcrSeq), pcrRoundsLeft_(pcrRoundsLeft), templateAmount_(
				templateAmount) {

}

PCRSimulator::PCRProduct::PCRProduct(const std::string & pcrSeq,
		const uint32_t pcrRoundsLeft) :
		pcrSeq_(pcrSeq), pcrRoundsLeft_(pcrRoundsLeft) {

}

void PCRSimulator::runPcr(uint32_t numThreads,
		uint32_t pcrRounds,
		std::string currentSeq,
		uint64_t currentStartingTemplate,
		std::unordered_map<std::string, uint64_t> & currentSeqMap,
		std::mutex & currentSeqMapLock) const {
	std::vector<njh::randomGenerator> gens;
	std::vector<std::unordered_map<char, njh::randObjectGen<char, uint32_t>>> charGens;
	for(uint32_t t = 0; t < numThreads; ++t){
		gens.emplace_back(njh::randomGenerator());
		std::unordered_map<char, njh::randObjectGen<char, uint32_t>> charGen;
		charGen.emplace('T', njh::randObjectGen<char, uint32_t>(std::vector<char>{'A', 'C', 'G'}, std::vector<uint32_t>{1,80,1}));
		charGen.emplace('C', njh::randObjectGen<char, uint32_t>(std::vector<char>{'A', 'T', 'G'}, std::vector<uint32_t>{1,80,1}));
		charGen.emplace('G', njh::randObjectGen<char, uint32_t>(std::vector<char>{'C', 'A', 'T'}, std::vector<uint32_t>{1,80,1}));
		charGen.emplace('A', njh::randObjectGen<char, uint32_t>(std::vector<char>{'C', 'G', 'T'}, std::vector<uint32_t>{1,80,1}));
		charGens.emplace_back(charGen);
	}

	std::vector<uint64_t> readPositions(currentSeq.length());
	njh::iota<uint64_t>(readPositions, 0);
	std::stack<PCRProduct> pcrProducts;
	pcrProducts.push(PCRProduct(currentSeq, pcrRounds, currentStartingTemplate));
	while(!pcrProducts.empty()){
		auto currentPCRProduct = pcrProducts.top();
		pcrProducts.pop();
		for(uint32_t round = 1; round <= currentPCRProduct.pcrRoundsLeft_; ++round){
			uint32_t duplicatingAmount = 0;
			for(uint32_t temp = 0; temp < currentPCRProduct.templateAmount_; ++temp){
				if(gens.front()() < pcrEfficiency_){
					++duplicatingAmount;
				}
			}

			std::unordered_map<uint32_t, uint32_t> muts;
			std::vector<std::string> mutants;
			std::mutex mutsLock;
			uint32_t numberMutated = 0;
			if(duplicatingAmount > 0){
				//generate the number of mutants and the number of bases mutated
				std::vector<uint64_t> tempAmounts;
				if(duplicatingAmount > numThreads){
					tempAmounts = spreadReadNumAcrossThreads(duplicatingAmount, numThreads);
				}else{
					tempAmounts.emplace_back(duplicatingAmount);
					for(uint32_t t = 1; t < numThreads; ++t){
						tempAmounts.emplace_back(0);
					}
				}
				auto mutate = [this, &readPositions,&mutsLock,&muts,&tempAmounts,&currentPCRProduct,&gens,&charGens,&numberMutated,&mutants](uint32_t threadNum) {

					if(tempAmounts[threadNum] > 0){
						std::unordered_map<uint32_t, uint32_t> currentMuts;
						std::vector<std::string> currentMutSeqs;
						for (uint32_t read = 0; read < tempAmounts[threadNum]; ++read) {
							uint32_t count = 0;
							for (uint32_t base = 0; base < currentPCRProduct.pcrSeq_.size(); ++base) {
								if (gens[threadNum].mtGen_() <= intErrorRate_) {
									++count;
								}
							}
							if (count != 0) {
								++currentMuts[count];
								currentMutSeqs.emplace_back(mutateSeq(currentPCRProduct.pcrSeq_, count, readPositions, gens[threadNum], charGens[threadNum]));
							}
						}

						//auto currentMuts = genMutCounts(gens[threadNum].mtGen_, tempAmounts[threadNum], currentPCRProduct.pcrSeq_.size());
						{
							std::lock_guard<std::mutex> lock(mutsLock);
							for(const auto & mut : currentMuts){
								muts[mut.first] += mut.second;
								numberMutated += mut.second;
							}
							addOtherVec(mutants, currentMutSeqs);
						}
					}
				};
				std::vector<std::thread> threads;
				for (uint32_t thread = 0; thread < numThreads; ++thread) {
					threads.emplace_back(std::thread(mutate,thread));
				}
				for(auto & t : threads){
					t.join();
				}
				if(verbose_){
					std::cout << "Round: " << round << std::endl;
					std::cout << '\t' << "seq: " << currentPCRProduct.pcrSeq_ << std::endl;
					std::cout << '\t' << "duplicatingAmount: " << duplicatingAmount << std::endl;
					std::cout << "\t" << "pcrRoundsLeft: " << currentPCRProduct.pcrRoundsLeft_ << std::endl;
					std::cout << "\t" << "templateAmount: " << currentPCRProduct.templateAmount_ << std::endl;
					std::cout << "\t" << "numberMutated: " << numberMutated << std::endl;
					std::cout << "\t" << "mutants.size(): " << mutants.size() << std::endl;
					std::cout << "\t" << "pcrProducts Left: " << pcrProducts.size() << std::endl;
					std::cout << "\t" << "muts " << std::endl;
					for(const auto & mut : muts){
						std::cout << "\t\t" << mut.first << "\t" << mut.second << std::endl;
					}
				}
			}
			//add on the number of seqs that weren't mutated;
			currentPCRProduct.templateAmount_ += duplicatingAmount - numberMutated;
			//check if this is the last round or not
			if (round == currentPCRProduct.pcrRoundsLeft_) {
				std::lock_guard<std::mutex> outMutLock(currentSeqMapLock);
				currentSeqMap[currentPCRProduct.pcrSeq_] +=
						currentPCRProduct.templateAmount_;
				for (const auto & mutSeq : mutants) {
					currentSeqMap[mutSeq] += 1;
				}
				//exit(1);
			} else {
				for (const auto & mutSeq : mutants) {
					//add on the new mutants to go through the rest of the PCR rounds
					pcrProducts.push(PCRProduct(mutSeq, currentPCRProduct.pcrRoundsLeft_ - round));
				}
			}
		}
	}
}

PCRSimulator::SeqGenomeCnt::SeqGenomeCnt(const seqInfo & seqBase,
		const uint64_t genomeCnt) :
		seqBase_(seqBase), genomeCnt_(genomeCnt) {
}


std::vector<PCRSimulator::SeqGenomeCnt> PCRSimulator::randomlySampleGenomes(
		const std::vector<seqInfo> & seqs, uint64_t totalGenomes) {
	std::vector<PCRSimulator::SeqGenomeCnt> ret;

	std::vector<uint32_t> seqPositions;
	std::vector<double> seqFracs;
	std::unordered_map<uint32_t, uint64_t> genomeCounts;
	for(const auto & seqPos : iter::range(seqs.size())){
		seqPositions.emplace_back(seqPos);
		seqFracs.emplace_back(seqs[seqPos].frac_);
		//to make sure if none is sampled the haplotype name still appears
		genomeCounts[seqPos] = 0;
	}
	njh::randObjectGen<uint32_t, double> templateSampler(seqPositions, seqFracs);


	for(uint32_t tempSampleCnt = 0; tempSampleCnt < totalGenomes; ++tempSampleCnt){
		auto hapSampled = templateSampler.genObj();
		++genomeCounts[hapSampled];
	}

	for(const auto & genomeCount : genomeCounts){
		ret.emplace_back(PCRSimulator::SeqGenomeCnt(seqs[genomeCount.first], genomeCount.second));
	}

	return ret;
}


PCRSimulator::SimHapCounts PCRSimulator::simLibFast(const std::vector<SeqGenomeCnt> & seqs,
		const OutOptions & outputFileOpts,
		uint64_t finalReadAmount,
		uint32_t pcrRounds,
		uint32_t initialPcrRounds,
		uint32_t numThreads){
	SimHapCounts ret;
	if(initialPcrRounds >= pcrRounds){
		std::stringstream ss;
		ss << "Error in " << __PRETTY_FUNCTION__ << std::endl;
		ss << "initialPcrRounds should be less than pcrRounds" << std::endl;
		ss << "pcrRounds:" << pcrRounds << std::endl;
		ss << "initialPcrRounds:" << initialPcrRounds << std::endl;
		throw std::runtime_error{ss.str()};
	}
	//auto finalPerfectAmount = static_cast<uint64_t>(startingTemplate * std::pow(2, pcrRounds));
	OutputStream libOutFile(outputFileOpts);
	std::mutex seqFileLock;

	//account for double stranded
	std::unordered_map<std::string, uint64_t> totalTemplateStrandCounts;
	std::unordered_map<std::string, uint64_t> genomeCounts;

	uint64_t totalTemplateStrandCount = 0;
	for(const auto & seq : seqs){
		genomeCounts[seq.seqBase_.name_] = seq.genomeCnt_;
		//account for double stranded
		totalTemplateStrandCounts[seq.seqBase_.name_] = seq.genomeCnt_ * 2;
		totalTemplateStrandCount += seq.genomeCnt_ * 2;
	}

	std::unordered_map<std::string,std::unordered_map<std::string, uint64_t>> allSeqCounts;
	std::unordered_map<std::string,std::string> barcodedSeqs;
//	std::unordered_map<std::string,std::pair<uint64_t,uint64_t>> templateNonMutated;
	for (const auto & seq : seqs) {
		std::unordered_map<std::string, uint64_t> seqCounts;
		std::mutex seqMapLock;

		runPcr(numThreads, initialPcrRounds, seq.seqBase_.seq_,
				totalTemplateStrandCounts[seq.seqBase_.name_], seqCounts, seqMapLock);

//		uint64_t finalPerfectAmount = templateAmountCounts[seq.name_] * std::pow(2, initialPcrRounds);

		//minus off the starting template amount as this is just genomic DNA and won't be able to be sequenced

		if(seqCounts[seq.seqBase_.seq_] <= totalTemplateStrandCounts[seq.seqBase_.name_]){
			seqCounts.erase(seq.seqBase_.seq_);
		}else{
			seqCounts[seq.seqBase_.seq_] -= totalTemplateStrandCounts[seq.seqBase_.name_];
		}
		if(!seqCounts.empty()){
			allSeqCounts[seq.seqBase_.name_] = seqCounts;
//			templateNonMutated[seq.name_] = {seqCounts[seq.seq_], finalPerfectAmount};
			barcodedSeqs[seq.seqBase_.name_] = seq.seqBase_.seq_;
		}
	}
	if(verbose_){
		uint64_t templateAmountAfterInitialPCR = 0;
		for(const auto  & allSeqCount : allSeqCounts){
			for(const auto & seqCount : allSeqCount.second){
				templateAmountAfterInitialPCR += seqCount.second;
			}
		}
		std::cout << "templateAmountAfterInitialPCR: " << templateAmountAfterInitialPCR << std::endl;
		std::cout << "ideal templateAmountAfterInitialPCR: " << totalTemplateStrandCount * std::pow(2, initialPcrRounds) - totalTemplateStrandCount << std::endl;
	}

	//now sample the rest of PCR
	auto sampleNumber = sampleReadsWithoutReplacementFinishPCR(barcodedSeqs,
			allSeqCounts, finalReadAmount, libOutFile, seqFileLock,
			pcrRounds - initialPcrRounds, numThreads);
	ret.genomesSampled_ = genomeCounts;
	ret.sampledForSequencing_ = sampleNumber;
	return ret;

//	if(verbose){
//		{
//			std::cout << "PCR amounts: " << std::endl;
//			uint64_t nonMutated = 0;
//			uint64_t mutated = 0;
//			for (const auto & read : seqs) {
//				std::cout << seq.name_ << std::endl;
//				nonMutated += templateNonMutated[seq.name_].first;
//				mutated += templateNonMutated[seq.name_].second
//						- templateNonMutated[seq.name_].first;
//				std::cout << "\t"
//						<< getPercentageString(templateNonMutated[seq.name_].first,
//								templateNonMutated[seq.name_].second) << std::endl;
//			}
//			std::cout << "Total: " << nonMutated + mutated << std::endl;
//			std::cout << "\t" << getPercentageString(nonMutated, nonMutated + mutated)
//					<< std::endl;
//			std::cout << "\t" << getPercentageString(mutated, nonMutated + mutated)
//					<< std::endl;
//		}
//		std::cout << "Sampling Amounts:" << std::endl;
//		uint64_t nonMutated = 0;
//		uint64_t mutated = 0;
//		for(const auto & seq : seqs){
//			std::cout << seq.name_ << std::endl;
//			auto total = sampleNumber[seq.name_].first + sampleNumber[seq.name_].second;
//			nonMutated+= sampleNumber[seq.name_].first;
//			mutated+= sampleNumber[seq.name_].second;
//			std::cout << "\tSampled     : " << getPercentageString(total, finalReadAmountCounts[seq.name_])<< std::endl;
//			std::cout << "\tNon-Mutated : " << getPercentageString(sampleNumber[seq.name_].first, total) << std::endl;
//			std::cout << "\tMutated     : " << getPercentageString(sampleNumber[seq.name_].second, total) << std::endl;
//		}
//		std::cout << "Total\t       : " << nonMutated + mutated << std::endl;
//		std::cout << "\tNon-Mutated : " << getPercentageString(nonMutated, nonMutated + mutated) << std::endl;
//		std::cout << "\tMutated     : " << getPercentageString(mutated, nonMutated + mutated) << std::endl;
//
//	}
}


PCRSimulator::SimHapCounts PCRSimulator::simLibFast(
		const std::vector<seqInfo> & seqs,
		const OutOptions & outputFileOpts,
		uint64_t startingTemplate,
		uint64_t finalReadAmount,
		uint32_t pcrRounds,
		uint32_t initialPcrRounds,
		uint32_t numThreads){
	auto seqsWithGenomes = randomlySampleGenomes(seqs, startingTemplate);
	return simLibFast(seqsWithGenomes, outputFileOpts, finalReadAmount, pcrRounds, initialPcrRounds, numThreads);
}



std::unordered_map<std::string, PCRSimulator::SimHapCounts::MutInfo> PCRSimulator::sampleReadsWithoutReplacementFinishPCR(
		const std::unordered_map<std::string, std::string> & seqs,
		std::unordered_map<std::string, std::unordered_map<std::string, uint64_t>> & multipleSeqCounts,
		uint64_t finalReadAmount, std::ostream & sequenceOutFile,
		std::mutex & seqFileLock, uint32_t numberOfPCRRoundsLeft,
		uint32_t numThreads) {

	std::unordered_map<std::string, PCRSimulator::SimHapCounts::MutInfo> ret;
	auto seqNames = getVectorOfMapKeys(seqs);
	auto multipleSeqCountsNames = getVectorOfMapKeys(multipleSeqCounts);
	njh::sort(seqNames);
	njh::sort(multipleSeqCountsNames);
	if(seqNames != multipleSeqCountsNames){
		std::stringstream ss;
		ss << "Error in " << __PRETTY_FUNCTION__ << " names for seqNames and multipleSeqCountsNames counts don't match" << std::endl;
		ss << "seqNames: " << njh::conToStr(seqNames,",") << std::endl;
		ss << "multipleSeqCountsNames: " << njh::conToStr(multipleSeqCountsNames,",") << std::endl;
		throw std::runtime_error{ss.str()};
	}
	uint64_t totalTemplate = 0;
	for (const auto & seqCounts : multipleSeqCounts) {
		for (const auto & seqCount : seqCounts.second) {
			totalTemplate += seqCount.second;
		}
	}

	if(verbose_){
		auto names = getVectorOfMapKeys(seqs);
		//std::cout << "Sampling from " << njh::conToStr(names,",") << " ..." << std::endl;
	}

	std::random_device rd;
	std::mt19937_64 mtGen(rd());
	std::unordered_map<uint32_t, std::vector<std::pair<std::string, std::string>>> allSampledSeqs;
	if(verbose_){
		std::cout << "Initial Sampling" << std::endl;
	}
	for (uint32_t read = 0; read < finalReadAmount; ++read) {
		if(verbose_){
			std::cout << '\r' << read + 1;
		}
		uint64_t randNumGen = mtGen();
		uint64_t randSel = (static_cast<double>(randNumGen)/mtGen.max()) * (totalTemplate - read);
		uint64_t sum = 0;
		for(auto & seqCounts : multipleSeqCounts){
			bool foundSelection = false;
			for (auto & seqCount : seqCounts.second) {
				if (seqCount.second == 0) {
					continue;
				}
				sum += seqCount.second;
				if (sum >= randSel) {
					--seqCount.second; //decrease to represent sampling without replacement
					allSampledSeqs[read % numThreads].emplace_back(std::pair<std::string,std::string>{seqCounts.first,seqCount.first});
					foundSelection = true;
					break;
				}
			}
			if(foundSelection){
				break;
			}
		}
	}
	if(verbose_){
		std::cout << std::endl;;
	}
	auto finishPCR =
			[this,&allSampledSeqs,&numberOfPCRRoundsLeft,&seqs,&seqFileLock,&sequenceOutFile,&ret,&finalReadAmount](uint32_t threadNumber ) {
				njh::stopWatch watch;
				std::vector<std::pair<std::string, std::string>> outputs;
				std::unordered_map<std::string, PCRSimulator::SimHapCounts::MutInfo> mutInfo;
				std::random_device rd;
				std::mt19937_64 mtGenFinal(rd());
				const uint32_t significantRounds = 10;
				for(const auto & namesSeqs : allSampledSeqs[threadNumber]) {
					std::string finalSeq = namesSeqs.second;
					uint32_t roundsStillLeft = numberOfPCRRoundsLeft;
					while(roundsStillLeft > significantRounds) {
						finalSeq = runPcrSampleSingleTemplate(
								significantRounds,mtGenFinal(), mtGenFinal.max(),
								finalSeq);
						roundsStillLeft -= significantRounds;
					}
					finalSeq = runPcrSampleSingleTemplate(
							roundsStillLeft,mtGenFinal(),mtGenFinal.max(),
							finalSeq);
					MetaDataInName nameMeta;
					nameMeta.addMeta("hap", namesSeqs.first);

					if (finalSeq != seqs.at(namesSeqs.first)) {
						nameMeta.addMeta("mutated", true);
						nameMeta.addMeta("idNum", leftPadNumStr<uint64_t>(mutInfo[namesSeqs.first].mutated_, finalReadAmount));
						++mutInfo[namesSeqs.first].mutated_;
					} else {
						nameMeta.addMeta("mutated", false);
						nameMeta.addMeta("idNum", leftPadNumStr<uint64_t>(mutInfo[namesSeqs.first].nonMutated_, finalReadAmount));
						++mutInfo[namesSeqs.first].nonMutated_;
					}
					outputs.emplace_back(std::make_pair(">" + nameMeta.createMetaName(), finalSeq));
				}
				{
					std::lock_guard<std::mutex> fileLock(seqFileLock);
					if(verbose_) {
						std::cout << "Thread " << threadNumber << " done" << std::endl;
						std::cout << "\tFinished PCR for "
								<< allSampledSeqs[threadNumber].size() << " reads" << std::endl;
						std::cout << "\tTime: " << watch.totalTimeFormatted(6) << std::endl;
					}
					for(const auto & out : outputs){
						sequenceOutFile << out.first << std::endl;
						sequenceOutFile << out.second << std::endl;
					}
					for(const auto & info : mutInfo){
						ret[info.first].mutated_ += info.second.mutated_;
						ret[info.first].nonMutated_ += info.second.nonMutated_;
					}
				}
			};

	std::vector<std::thread> threads;
	for(uint32_t thread = 0; thread < numThreads; ++thread){
		threads.emplace_back(std::thread(finishPCR, thread));
	}

	for(auto & thread : threads){
		thread.join();
	}
	return ret;
}


std::string PCRSimulator::runPcrSampleSingleTemplate(
		uint32_t roundsOfPcr,
		uint64_t randomNumberSelector,
		uint64_t randomNumberSelectorMax,
		std::string seq ){
	std::unordered_map<std::string, uint64_t> finalProducts;
	std::mutex finalProductsMutex;
	runPcr(1, roundsOfPcr, seq, 1, finalProducts, finalProductsMutex);
	uint64_t finalProductSum = 0;
	for(const auto & final : finalProducts) {
		finalProductSum += final.second;
	}
	uint64_t finalRandSel = (static_cast<double>(randomNumberSelector)/randomNumberSelectorMax) * (finalProductSum );
	uint64_t finalSum = 0;
	for(const auto & final : finalProducts) {
		finalSum += final.second;
		if(finalSum >= finalRandSel) {
			seq = final.first;
			break;
		}
	}
	return seq;
}


}  // namespace njhseq


