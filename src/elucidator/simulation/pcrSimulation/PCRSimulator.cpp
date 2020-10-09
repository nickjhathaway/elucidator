/*
 * PCRSimulator.cpp
 *
 *  Created on: May 2, 2019
 *      Author: nicholashathaway
 */


#include "PCRSimulator.hpp"

namespace njhseq {


class ChimeraProfiler{
public:
	struct SharedSegments {

		struct Segment{
			Segment(const uint32_t starta, const uint32_t startb, uint32_t size) :
					starta_(starta), startb_(startb), size_(size) {

			}
			uint32_t starta_;
			uint32_t startb_;
			uint32_t size_;
			uint32_t enda() const{
				return starta_ + size_;
			}
			uint32_t endb() const{
				return startb_ + size_;
			}

			bool doesSegsOverLapInA(const Segment & otherSeg){
				if((otherSeg.starta_ > starta_ && otherSeg.starta_ < enda()) || (otherSeg.enda() > starta_ && otherSeg.enda() < enda())){
					return true;
				}
				return false;
			}
			bool doesSegsOverLapInB(const Segment & otherSeg){
				if((otherSeg.startb_ > startb_ && otherSeg.startb_ < endb()) || (otherSeg.endb() > startb_ && otherSeg.endb() < endb())){
					return true;
				}
				return false;
			}
		};

		SharedSegments(const std::string & namea, const std::string & nameb) :
				seqNameA_(namea), seqNameB_(nameb) {

		}

		void addSegment(Segment newSeg){
			//check to make sure there's no overlap
			for(const auto & seg : segs_){
				if(newSeg.doesSegsOverLapInA(seg)){
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error " << "adding segment,"
							<< "start: " << newSeg.starta_ << " size: " << newSeg.size_
							<< " that overlaps with "
							<<"start: " <<  seg.starta_<< " size: " << seg.size_<< "\n";
					throw std::runtime_error{ss.str()};
				}
				if(newSeg.doesSegsOverLapInA(seg)){
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error " << "adding segment,"
							<< "start: " << newSeg.starta_ << " size: " << newSeg.size_
							<< " that overlaps with "
							<<"start: " <<  seg.starta_<< " size: " << seg.size_<< "\n";
					throw std::runtime_error{ss.str()};
				}
			}
			segs_.emplace_back(newSeg);
		}
		void addSegment(uint32_t starta, uint32_t startb, uint32_t size){
			addSegment(Segment{starta, startb, size});

		}

		std::string seqNameA_;
		std::string seqNameB_;

		std::vector<Segment> segs_;

	};

	struct SegmentsForSeqs {
		SegmentsForSeqs(const seqInfo & seqBase) :
			seqBase_(seqBase) {

		}
		seqInfo seqBase_;
		std::unordered_map<std::string, SharedSegments> sharedSegs_;

		uint32_t minSegmentStart(){
			uint32_t ret = std::numeric_limits<uint32_t>::max();
			for(const auto & seg : sharedSegs_){
				if(seg.second.segs_.front().starta_ < ret){
					ret = seg.second.segs_.front().starta_;
				}
			}
			return ret;
		}
		uint32_t maxSegmentEnd(){
			uint32_t ret = 0;
			for(const auto & seg : sharedSegs_){
				if(seg.second.segs_.back().enda() > ret){
					ret = seg.second.segs_.back().enda();
				}
			}
			return ret;
		}

		void addSegments(const seqInfo & alnA, const seqInfo & alnB, uint32_t padding){
			if(njh::in(alnB.name_, sharedSegs_)){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << " error, already have segments for " << alnB.name_<< "\n";
				throw std::runtime_error{ss.str()};
			}

			uint32_t aLen = len(alnA) - countOccurences(alnA.seq_, "-");
			uint32_t bLen = len(alnB) - countOccurences(alnB.seq_, "-");

			uint32_t firstNonGap = std::max(alnA.seq_.find_first_not_of("-"), alnB.seq_.find_first_not_of("-"));
			uint32_t lastNonGap = std::max(alnA.seq_.find_last_not_of("-"), alnB.seq_.find_last_not_of("-"));
			uint32_t currentSegLength = 0;
			uint32_t currentSegStart = std::numeric_limits<uint32_t>::max();
			SharedSegments segs(alnA.name_, alnB.name_);

			for(const auto pos : iter::range(firstNonGap, lastNonGap + 1)){
				if(alnA.seq_[pos] != '-' && alnA.seq_[pos] == alnB.seq_[pos]){
					if(0 == currentSegLength){
						currentSegLength = 1;
						currentSegStart = pos;
					}else{
						++currentSegLength;
					}
				}else{
					if(currentSegLength > 0 && currentSegLength > padding ){
						uint32_t startA = getRealPosForAlnPos(alnA.seq_, currentSegStart);
						uint32_t startB = getRealPosForAlnPos(alnB.seq_, currentSegStart);
//						std::cout << "pos: " << pos << std::endl;
//						std::cout <<"startA + currentSegLength: " << startA + currentSegLength << std::endl;
//						std::cout <<"aLen                     : " << aLen << std::endl;
//						std::cout <<"startB + currentSegLength: " << startB + currentSegLength << std::endl;
//						std::cout <<"bLen                     : " << bLen << std::endl;
//						std::cout << std::endl;
						if(!(startA + currentSegLength == aLen && startB + currentSegLength == bLen) && !(0 == startA && 0 == startB)){
							segs.addSegment(startA + padding, startB + padding, currentSegLength - padding);
						}
					}
					currentSegLength = 0;
					currentSegStart = std::numeric_limits<uint32_t>::max();
				}
			}
			if(currentSegLength > 0 && currentSegLength > padding ){
				uint32_t startA = getRealPosForAlnPos(alnA.seq_, currentSegStart);
				uint32_t startB = getRealPosForAlnPos(alnB.seq_, currentSegStart);
//				std::cout << "after" << std::endl;
//				std::cout <<"startA + currentSegLength: " << startA + currentSegLength << std::endl;
//				std::cout <<"aLen                     : " << aLen << std::endl;
//				std::cout <<"startB + currentSegLength: " << startB + currentSegLength << std::endl;
//				std::cout <<"bLen                     : " << bLen << std::endl;

				if(!(startA + currentSegLength == aLen && startB + currentSegLength == bLen) && !(0 == startA && 0 == startB)){
					segs.addSegment(startA + padding, startB + padding, currentSegLength - padding);
				}
			}
			if(!segs.segs_.empty()){
				sharedSegs_.emplace(alnB.name_, segs);
			}
		}

		struct SeqNameSeg {
			SeqNameSeg(const std::string & name, uint32_t pos):name_(name), pos_(pos){

			}
			std::string name_;
			uint32_t pos_;

		};

		std::vector<SeqNameSeg> genPosChimeras(uint32_t pos){
			std::vector<SeqNameSeg> ret;
			for(const auto & seq : sharedSegs_){
				for(const auto & seg : seq.second.segs_){
					if(pos >= seg.starta_ && pos < seg.enda()){
						ret.emplace_back(seq.first, seg.startb_ + (pos - seg.starta_));
					}
				}
			}
			return ret;
		}

	};

	std::unordered_map<std::string, std::shared_ptr<SegmentsForSeqs>> segsForSeqs_;


	template<typename T>
	void addSegsForSeqs(const std::vector<T> & seqs, aligner & alignerObj, uint32_t padding){
		segsForSeqs_.clear();
		for(const auto & seq : seqs){
			segsForSeqs_[getSeqBase(seq).name_] = std::make_shared<SegmentsForSeqs>(getSeqBase(seq));
		}
		for (auto pos1 : iter::range(seqs.size())) {
			for (auto pos2 : iter::range(pos1 + 1, seqs.size())) {
				alignerObj.alignCacheGlobal(getSeqBase(seqs[pos1]), getSeqBase(seqs[pos2]));
				alignerObj.profileAlignment(getSeqBase(seqs[pos1]), getSeqBase(seqs[pos2]), false, false, false);
				segsForSeqs_[getSeqBase(seqs[pos1]).name_]->addSegments(alignerObj.alignObjectA_.seqBase_,
						alignerObj.alignObjectB_.seqBase_, padding);
				segsForSeqs_[getSeqBase(seqs[pos2]).name_]->addSegments(alignerObj.alignObjectB_.seqBase_,
						alignerObj.alignObjectA_.seqBase_, padding);
			}
		}
	}

	template<typename T>
	void addSegsForSeq(const T & seq, aligner & alignerObj, uint32_t padding){
		if(njh::in(getSeqBase(seq).name_, segsForSeqs_)){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << " already have seq with name " << getSeqBase(seq).name_ << "\n";
			throw std::runtime_error{ss.str()};
		}
		segsForSeqs_[getSeqBase(seq).name_] = std::make_shared<SegmentsForSeqs>(getSeqBase(seq));
		for (const auto & otherSeq : segsForSeqs_) {
			if(otherSeq.first == getSeqBase(seq).name_){
				continue;
			}
			alignerObj.alignCacheGlobal(otherSeq.second->seqBase_, getSeqBase(seq));
			alignerObj.profileAlignment(otherSeq.second->seqBase_, getSeqBase(seq), false, false, false);
			segsForSeqs_[otherSeq.second->seqBase_.name_]->addSegments(alignerObj.alignObjectA_.seqBase_,
					alignerObj.alignObjectB_.seqBase_, padding);
			segsForSeqs_[getSeqBase(seq).name_]->addSegments(alignerObj.alignObjectB_.seqBase_,
					alignerObj.alignObjectA_.seqBase_, padding);
		}
	}

};



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
				if(numThreads <=1){
					mutate(0);
				}else{
					std::vector<std::thread> threads;
					for (uint32_t thread = 0; thread < numThreads; ++thread) {
						threads.emplace_back(std::thread(mutate,thread));
					}
					for(auto & t : threads){
						t.join();
					}
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
	bool allFracZero = std::all_of(seqs.begin(), seqs.end(), [](const seqInfo & seq){
		return 0 == seq.frac_;
	});

//	bool allCntZero = std::all_of(seqs.begin(), seqs.end(), [](const seqInfo & seq){
//		return 0 == seq.cnt_;
//	});

	if(allFracZero){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << " all the fractions can't be set to zero" << "\n";
		ss << "Fracs: ";
		njh::for_each(seqs, [&ss](const seqInfo & seq){
			ss << "\n" << seq.name_ << ":" << seq.frac_;
		});
		ss << "\n";
		throw std::runtime_error{ss.str()};
	}
	for(const auto seqPos : iter::range(seqs.size())){

//		std::cout << "seqs[seqPos].frac_: " << seqs[seqPos].frac_ << std::endl;
//		std::cout << "seqs[seqPos].cnt_: " << seqs[seqPos].cnt_ << std::endl;

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



	uint32_t numberOfSeqsSampled = 0;
	std::vector<SeqGenomeCnt> seqsSampled;
	for (const auto & seq : seqs) {
		//if the gnomes sampled is zero don't simulate any seqs
		if(genomeCounts[seq.seqBase_.name_] > 0){
			++numberOfSeqsSampled;
			seqsSampled.emplace_back(seq);
		}
	}

	std::unordered_map<std::string,std::unordered_map<std::string, uint64_t>> allSeqCounts;
	std::unordered_map<std::string,std::string> barcodedSeqs;
//	std::unordered_map<std::string,std::pair<uint64_t,uint64_t>> templateNonMutated;
	for (const auto & seq : seqsSampled) {
		//if the gnomes sampled is zero don't simulate any seqs
		if(0 == genomeCounts[seq.seqBase_.name_]){
			continue;
		}
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

	uint64_t totalFinalChimeraCount = 0;
	std::vector<seqInfo> chimeraSeqs;

	if(!noChimeras_){
		if(numberOfSeqsSampled > 1){
			auto maxLen = readVec::getMaxLength(seqsSampled);
			ChimeraProfiler chiProf;
			aligner alignerObj(maxLen, gapScoringParameters(5,1,0,0,0,0));
			chiProf.addSegsForSeqs(seqsSampled, alignerObj, chimeraPad_);
			uint32_t seqsWithSegs = 0;
			for(const auto & seq : chiProf.segsForSeqs_){
				if(!seq.second->sharedSegs_.empty()){
					++seqsWithSegs;
				}
			}

//			std::cout << "seqsWithSegs: " << seqsWithSegs << std::endl;
			if(seqsWithSegs > 0 ){
				std::unordered_map<std::string, uint32_t> seqPositions;
				for(const auto pos : iter::range(seqsSampled.size())){
					seqPositions[seqsSampled[pos].seqBase_.name_] = pos;
				}


//				for(const auto & seq : chiProf.segsForSeqs_){
//					std::cout << seq.first << std::endl;
//					for(const auto & otherSeq : seq.second->sharedSegs_){
//						std::cout << "\t" << otherSeq.first << std::endl;
//						for(const auto & seg : otherSeq.second.segs_){
//							std::cout << "\t\tseg.starta_:" << seg.starta_ << std::endl;
//							std::cout << "\t\tseg.startb_:" << seg.startb_ << std::endl;
//							std::cout << "\t\tseg.size_:" << seg.size_ << std::endl;
//							seqsSampled[seqPositions[seq.first]].seqBase_.getSubRead(seg.starta_, seg.size_).outPutSeqAnsi(std::cout);
//							seqsSampled[seqPositions[otherSeq.first]].seqBase_.getSubRead(seg.startb_, seg.size_).outPutSeqAnsi(std::cout);
//						}
//					}
//				}
				njh::randomGenerator rGen;
				std::unordered_map<std::string, uint64_t> finishedAmount;
				std::unordered_map<std::string, uint64_t> partialTotalAmount;
				uint64_t intialRoundChiTotal = 0;
				std::unordered_map<std::string, uint64_t> chimeras;
				std::mutex chimerasMut;

				{

					table roundCountsTab(VecStr{"round", "seqName", "total", "fullTemplateFinished", "fullTemplateThisRound", "partialThisRound", "fracPartialThisRound", "partialTotalAmount", "partialTotalAmountFrac"});
					for(const auto & seqSampled : seqsSampled){
						finishedAmount[seqSampled.seqBase_.name_] = seqSampled.genomeCnt_ * 2;
						partialTotalAmount[seqSampled.seqBase_.name_] = 0;
					}

					for(auto round : iter::range<uint32_t>(1, initialPcrRounds + 1)){
//						std::cout << "round: " << round << std::endl;
						for(const auto & seqSampled : seqsSampled){

//							std::cout << "seq: " << seqSampled.seqBase_.name_ << std::endl;
							auto chiStart = chiProf.segsForSeqs_[seqSampled.seqBase_.name_]->minSegmentStart();
							auto chiStop = chiProf.segsForSeqs_[seqSampled.seqBase_.name_]->maxSegmentEnd();
							uint64_t finished = 0;
							uint64_t partial = 0;
							if (round > 10) {
								finished += std::min<uint64_t>(std::round(pcrEfficiency_ * finishedAmount[seqSampled.seqBase_.name_]), templateCap_);
		 						partial +=  std::min<uint64_t>(std::round((1 - pcrEfficiency_) * finishedAmount[seqSampled.seqBase_.name_]), templateCap_);
							} else {
								for(uint64_t temp = 0; temp < finishedAmount[seqSampled.seqBase_.name_]; ++temp){
									if(rGen.unifRand() < pcrEfficiency_){
										++finished;
									}else{
										++partial;
									}
								}
							}
							if(partial > 0){
								std::vector<uint64_t> partials;
								if(partial > numThreads){
									partials = spreadReadNumAcrossThreads(partial, numThreads);
								}else{
									partials.emplace_back(partial);
									for(uint32_t t = 1; t < numThreads; ++t){
										partials.emplace_back(0);
									}
								}
								auto genPartials = [&partials,&seqSampled,&chiProf,&chiStart,&chiStop,&chimerasMut,&chimeras,&finishedAmount](const uint32_t threadNum){
									njh::randomGenerator ranGen;
									std::unordered_map<std::string, uint32_t> threadChimeras;
									uint32_t chisGenerated = 0;
									for(uint64_t temp = 0; temp < partials[threadNum]; ++temp){
										//ranGen.unifRand();
										auto partLength = ranGen.unifRand<uint32_t>(0, seqSampled.seqBase_.seq_.size());
//											std::cout << "chiStart: " << chiStart  << " chiStop:" << chiStop << std::endl;
//											std::cout << "partLength: " << partLength << " chi pos? " << njh::colorBool(partLength >= chiStart && partLength <chiStop) << std::endl;
										if(partLength >= chiStart && partLength <chiStop){
											auto posChis = chiProf.segsForSeqs_[seqSampled.seqBase_.name_]->genPosChimeras(partLength);
											VecStr names;
											names.emplace_back(seqSampled.seqBase_.name_);
//												std::cout << "Pos chis" << std::endl;
											std::unordered_map<std::string, uint32_t> otherTempPos;
											for(const auto & posChi : posChis){
												names.emplace_back(posChi.name_);
//													std::cout << posChi.name_ << " " << posChi.pos_ << std::endl;
												otherTempPos[posChi.name_] = posChi.pos_;
											}
											std::vector<double> fracs;
											for(const auto & name : names){
												fracs.emplace_back(finishedAmount[name]);
											}
											njh::randObjectGen gen(names,fracs);
											std::string templatedLandedOn = gen.genObj();
//												std::cout << "templatedLandedOn different? " << njh::colorBool(templatedLandedOn != seqSampled.seqBase_.name_) << std::endl;
											if(templatedLandedOn != seqSampled.seqBase_.name_){
												++chisGenerated;
												++threadChimeras[seqSampled.seqBase_.seq_.substr(0, partLength) + chiProf.segsForSeqs_[templatedLandedOn]->seqBase_.seq_.substr(otherTempPos[templatedLandedOn])];
											}
										}
									}
									{
										std::lock_guard<std::mutex> lock(chimerasMut);
										for(const auto & chi : threadChimeras){
											chimeras[chi.first] += chi.second;
										}
									}
//										std::cout << "chisGenerated: " << chisGenerated << std::endl;
//										std::cout << "partials[threadNum]: " << partials[threadNum] << std::endl;
								};

								{
									if(numThreads <=1){
										genPartials(0);
									}else{
										std::vector<std::thread> threads;
										for(uint32_t t = 0; t < numThreads; ++t){
											threads.emplace_back(genPartials,t);
										}
										njh::concurrent::joinAllThreads(threads);
									}
								}
							}
							partialTotalAmount[seqSampled.seqBase_.name_] += partial;
							finishedAmount[seqSampled.seqBase_.name_] += finished;
							roundCountsTab.addRow(round, seqSampled.seqBase_.name_, partialTotalAmount[seqSampled.seqBase_.name_] + finishedAmount[seqSampled.seqBase_.name_], finishedAmount[seqSampled.seqBase_.name_],
									finished, partial, static_cast<double>(partial)/(finished + partial),
									partialTotalAmount[seqSampled.seqBase_.name_], static_cast<double>(partialTotalAmount[seqSampled.seqBase_.name_])/(partialTotalAmount[seqSampled.seqBase_.name_] + finishedAmount[seqSampled.seqBase_.name_]));
						}

						for(const auto & chiSeq : chimeraSeqs){

//							std::cout << "seq: " << chiSeq.name_ << std::endl;
							auto chiStart = chiProf.segsForSeqs_[chiSeq.name_]->minSegmentStart();
							auto chiStop = chiProf.segsForSeqs_[chiSeq.name_]->maxSegmentEnd();
							uint64_t finished = 0;
							uint64_t partial = 0;
							if (round > 10) {
								finished += std::min<uint64_t>(std::round(pcrEfficiency_ * finishedAmount[chiSeq.name_]), templateCap_);
		 						partial +=  std::min<uint64_t>(std::round((1 - pcrEfficiency_) * finishedAmount[chiSeq.name_]), templateCap_);
							} else {
								for(uint64_t temp = 0; temp < finishedAmount[chiSeq.name_]; ++temp){
									if(rGen.unifRand() < pcrEfficiency_){
										++finished;
									}else{
										++partial;
									}
								}
							}
							if(partial > 0){
								std::vector<uint64_t> partials;
								if(partial > numThreads){
									partials = spreadReadNumAcrossThreads(partial, numThreads);
								}else{
									partials.emplace_back(partial);
									for(uint32_t t = 1; t < numThreads; ++t){
										partials.emplace_back(0);
									}
								}
								auto genPartials = [&partials,&chiSeq,&chiProf,&chiStart,&chiStop,&chimerasMut,&chimeras,&finishedAmount](const uint32_t threadNum){
									//std::unordered_map<std::string, uint32_t> chimeras;
									njh::randomGenerator ranGen;
									std::unordered_map<std::string, uint32_t> threadChimeras;
									uint32_t chisGenerated = 0;
									for(uint64_t temp = 0; temp < partials[threadNum]; ++temp){
										//ranGen.unifRand();
										auto partLength = ranGen.unifRand<uint32_t>(0, chiSeq.seq_.size());
//											std::cout << "chiStart: " << chiStart  << " chiStop:" << chiStop << std::endl;
//											std::cout << "partLength: " << partLength << " chi pos? " << njh::colorBool(partLength >= chiStart && partLength <chiStop) << std::endl;
										if(partLength >= chiStart && partLength <chiStop){
											auto posChis = chiProf.segsForSeqs_[chiSeq.name_]->genPosChimeras(partLength);
											VecStr names;
											names.emplace_back(chiSeq.name_);
//												std::cout << "Pos chis" << std::endl;
											std::unordered_map<std::string, uint32_t> otherTempPos;
											for(const auto & posChi : posChis){
												names.emplace_back(posChi.name_);
//													std::cout << posChi.name_ << " " << posChi.pos_ << std::endl;
												otherTempPos[posChi.name_] = posChi.pos_;
											}
											std::vector<double> fracs;
											for(const auto & name : names){
												fracs.emplace_back(finishedAmount[name]);
											}
											njh::randObjectGen gen(names,fracs);
											std::string templatedLandedOn = gen.genObj();
//												std::cout << "templatedLandedOn different? " << njh::colorBool(templatedLandedOn != seqSampled.seqBase_.name_) << std::endl;
											if(templatedLandedOn != chiSeq.name_){
												++chisGenerated;
												++threadChimeras[chiSeq.seq_.substr(0, partLength) + chiProf.segsForSeqs_[templatedLandedOn]->seqBase_.seq_.substr(otherTempPos[templatedLandedOn])];
											}
										}
									}
									{
										std::lock_guard<std::mutex> lock(chimerasMut);
										for(const auto & chi : threadChimeras){
											chimeras[chi.first] += chi.second;
										}
									}
//										std::cout << "chisGenerated: " << chisGenerated << std::endl;
//										std::cout << "partials[threadNum]: " << partials[threadNum] << std::endl;
								};

								{
									if(numThreads <=1){
										genPartials(0);
									}else{
										std::vector<std::thread> threads;
										for(uint32_t t = 0; t < numThreads; ++t){
											threads.emplace_back(genPartials,t);
										}
										njh::concurrent::joinAllThreads(threads);
									}
								}

							}

							partialTotalAmount[chiSeq.name_] += partial;
							finishedAmount[chiSeq.name_] += finished;

							roundCountsTab.addRow(round, chiSeq.name_, partialTotalAmount[chiSeq.name_] + finishedAmount[chiSeq.name_], finishedAmount[chiSeq.name_],
									finished, partial, static_cast<double>(partial)/(finished + partial),
									partialTotalAmount[chiSeq.name_], static_cast<double>(partialTotalAmount[chiSeq.name_])/(partialTotalAmount[chiSeq.name_] + finishedAmount[chiSeq.name_]));
						}
						//add in chimeras
						for(const auto & chi : chimeras){
							bool matchesInputSeq = false;
							for(const auto & seqSampled : seqsSampled){
								if(seqSampled.seqBase_.seq_ == chi.first){
									matchesInputSeq = true;
									finishedAmount[seqSampled.seqBase_.name_] += chi.second;
									break;
								}
							}
							bool matchesPreviousChimera = false;
							for(const auto & chimeraSeq : chimeraSeqs){
								if(chimeraSeq.seq_ == chi.first){
									matchesPreviousChimera = true;
									finishedAmount[chimeraSeq.name_] += chi.second;
									break;
								}
							}
							if(!matchesInputSeq && !matchesPreviousChimera){
								std::string name =njh::pasteAsStr("Chi.", chimeraSeqs.size(), "[PCRRound=", round, "]");
								auto backSeq = chimeraSeqs.emplace_back(seqInfo(name, chi.first));
								backSeq.cnt_ = chi.second;
								chiProf.addSegsForSeq(backSeq, alignerObj, chimeraPad_);
								finishedAmount[backSeq.name_] = chi.second;
								partialTotalAmount[backSeq.name_] = 0;
							}
						}
					}
//					roundCountsTab.outPutContentOrganized(std::cout);
					for(const auto & chi : chimeras){
						bool matchesInputSeq = false;
						for(const auto & seqSampled : seqsSampled){
							if(seqSampled.seqBase_.seq_ == chi.first){
								matchesInputSeq = true;
								finishedAmount[seqSampled.seqBase_.name_] += chi.second;
								allSeqCounts[seqSampled.seqBase_.name_][seqSampled.seqBase_.seq_]+= chi.second;
								break;
							}
						}
						bool matchesPreviousChimera = false;
						if(!matchesInputSeq){
							intialRoundChiTotal+= chi.second;
							for(const auto & chimeraSeq : chimeraSeqs){
								if(chimeraSeq.seq_ == chi.first){
									matchesPreviousChimera = true;
									finishedAmount[chimeraSeq.name_] += chi.second;
									allSeqCounts[chimeraSeq.name_][chimeraSeq.seq_] = chi.second;
						//			templateNonMutated[seq.name_] = {seqCounts[seq.seq_], finalPerfectAmount};
									barcodedSeqs[chimeraSeq.name_] = chimeraSeq.seq_;
									break;
								}
							}
						}
						if(!matchesInputSeq && !matchesPreviousChimera){
							std::string name =njh::pasteAsStr("Chi.", chimeraSeqs.size(), "[PCRRound=", initialPcrRounds, "]");
							auto backSeq = chimeraSeqs.emplace_back(seqInfo(name, chi.first));
							backSeq.cnt_ = chi.second;
							chiProf.addSegsForSeq(backSeq, alignerObj, chimeraPad_);
							finishedAmount[backSeq.name_] = chi.second;
							partialTotalAmount[backSeq.name_] = 0;
							allSeqCounts[backSeq.name_][backSeq.seq_] = chi.second;
				//			templateNonMutated[seq.name_] = {seqCounts[seq.seq_], finalPerfectAmount};
							barcodedSeqs[backSeq.name_] = backSeq.seq_;
						}

//						if(matchesInputSeq){
//							std::cout << njh::bashCT::green << std::endl;
//						}else if(matchesPreviousChimera){
//							std::cout << njh::bashCT::red << std::endl;
//						}
//
//						std::cout << chi.first << std::endl;
//						std::cout << chi.second << std::endl << std::endl;
//						if(matchesInputSeq || matchesPreviousChimera){
//							std::cout << njh::bashCT::reset << std::endl;
//						}
					}
					uint64_t all_finishedAmount = 0;
					for(const auto & finised : finishedAmount){
						all_finishedAmount += finised.second;
					}
//					std::cout << "chiTotal: " << intialRoundChiTotal << " " << static_cast<double>(intialRoundChiTotal)/(all_finishedAmount + intialRoundChiTotal) << std::endl;
					//exit(1);
				}
				{
					uint64_t all_intialFinishedAmount = 0;
					for(const auto & finised : finishedAmount){
						all_intialFinishedAmount += finised.second;
					}

					table roundCountsTab(VecStr{"round", "seqName", "total", "fullTemplateFinished", "fullTemplateThisRound", "partialThisRound", "fracPartialThisRound", "partialTotalAmount", "partialTotalAmountFrac"});
					std::unordered_map<uint32_t, uint64_t> roundPartialCounts;

					for(auto round : iter::range<uint32_t>(initialPcrRounds + 1, pcrRounds + 1)){
						uint64_t roundPartials = 0;
//						std::cout << "round: " << round << std::endl;
						for(const auto & seqSampled : seqsSampled){
//							std::cout << "seq: " << seqSampled.seqBase_.name_ << std::endl;
							uint64_t finished = 0;
							uint64_t partial = 0;
							if (round > 10) {
								finished += std::min<uint64_t>(std::round(pcrEfficiency_ * finishedAmount[seqSampled.seqBase_.name_]), templateCap_);
		 						partial +=  std::min<uint64_t>(std::round((1 - pcrEfficiency_) * finishedAmount[seqSampled.seqBase_.name_]), templateCap_);
							} else {
								for(uint64_t temp = 0; temp < finishedAmount[seqSampled.seqBase_.name_]; ++temp){
									if(rGen.unifRand() < pcrEfficiency_){
										++finished;
									}else{
										++partial;
									}
								}
							}
							roundPartials+= partial;
							partialTotalAmount[seqSampled.seqBase_.name_] += partial;
							finishedAmount[seqSampled.seqBase_.name_] += finished;
							roundCountsTab.addRow(round, seqSampled.seqBase_.name_, partialTotalAmount[seqSampled.seqBase_.name_] + finishedAmount[seqSampled.seqBase_.name_], finishedAmount[seqSampled.seqBase_.name_],
									finished, partial, static_cast<double>(partial)/(finished + partial),
									partialTotalAmount[seqSampled.seqBase_.name_], static_cast<double>(partialTotalAmount[seqSampled.seqBase_.name_])/(partialTotalAmount[seqSampled.seqBase_.name_] + finishedAmount[seqSampled.seqBase_.name_]));
						}
						roundPartialCounts[round] = roundPartials;
					}
//					roundCountsTab.outPutContentOrganized(std::cout);
					uint64_t allFinishedPlusPartial = 0;
					std::unordered_map<std::string, double> percPartial;
					for(const auto & seqSampled : seqsSampled){
						allFinishedPlusPartial += finishedAmount[seqSampled.seqBase_.name_] + partialTotalAmount[seqSampled.seqBase_.name_];
						percPartial[seqSampled.seqBase_.name_] = static_cast<double>(partialTotalAmount[seqSampled.seqBase_.name_])/(partialTotalAmount[seqSampled.seqBase_.name_] + finishedAmount[seqSampled.seqBase_.name_]);
					}
					VecStr seqNames;
					std::vector<double> seqPercs;
					for(const auto & seqSampled: seqsSampled){
						double perc = static_cast<double>(finishedAmount[seqSampled.seqBase_.name_] + partialTotalAmount[seqSampled.seqBase_.name_])/allFinishedPlusPartial;
						seqPercs.emplace_back(perc);
						seqNames.emplace_back(seqSampled.seqBase_.name_);
					}
					njh::randObjectGen rObj(seqNames, seqPercs);



					std::vector<uint32_t> rounds;
					std::vector<uint64_t> roundsCounts;
					for(const auto & roundPar : roundPartialCounts){
						rounds.emplace_back(roundPar.first);
						roundsCounts.emplace_back(roundPar.second);
					}
					njh::randObjectGen roundGen(rounds, roundsCounts);



					std::unordered_map<std::string, uint64_t> finalPartials;
					for(uint32_t r = 0; r < finalReadAmount - (finalReadAmount *(static_cast<double>(intialRoundChiTotal)/(all_intialFinishedAmount + intialRoundChiTotal))); ++r){
						auto seqSampled = rObj.genObj();
						if(rGen.unifRand() < percPartial[seqSampled]){
							++finalPartials[seqSampled];
						}
					}

					std::unordered_map<uint32_t, std::unordered_map<std::string, uint64_t>> finalChimerasPerRound;
					std::unordered_map<std::string, uint64_t> finalChimeras;
					for(const auto & finalPartial : finalPartials){
//						std::cout << "seq: " << finalPartial.first << " final partial: " << finalPartial.second << std::endl;
						auto chiStart = chiProf.segsForSeqs_[finalPartial.first]->minSegmentStart();
						auto chiStop = chiProf.segsForSeqs_[finalPartial.first]->maxSegmentEnd();

						std::vector<uint64_t> partials;
						if(finalPartial.second > numThreads){
							partials = spreadReadNumAcrossThreads(finalPartial.second, numThreads);
						}else{
							partials.emplace_back(finalPartial.second);
							for(uint32_t t = 1; t < numThreads; ++t){
								partials.emplace_back(0);
							}
						}

						auto genPartials = [&finalPartial,&partials,&chiProf,&chiStart,&chiStop,&chimerasMut,&finishedAmount,&finalChimeras,&roundGen,&finalChimerasPerRound](const uint32_t threadNum){
							njh::randomGenerator ranGen;
							std::unordered_map<std::string, uint32_t> threadChimeras;
							uint32_t chisGenerated = 0;
							for(uint64_t temp = 0; temp < partials[threadNum]; ++temp){
								//ranGen.unifRand();
								auto partLength = ranGen.unifRand<uint32_t>(0, chiProf.segsForSeqs_[finalPartial.first]->seqBase_.seq_.size());
//											std::cout << "chiStart: " << chiStart  << " chiStop:" << chiStop << std::endl;
//											std::cout << "partLength: " << partLength << " chi pos? " << njh::colorBool(partLength >= chiStart && partLength <chiStop) << std::endl;
								if(partLength >= chiStart && partLength <chiStop){
									auto posChis = chiProf.segsForSeqs_[chiProf.segsForSeqs_[finalPartial.first]->seqBase_.name_]->genPosChimeras(partLength);
									VecStr names;
									names.emplace_back(chiProf.segsForSeqs_[finalPartial.first]->seqBase_.name_);
//												std::cout << "Pos chis" << std::endl;
									std::unordered_map<std::string, uint32_t> otherTempPos;
									for(const auto & posChi : posChis){
										names.emplace_back(posChi.name_);
//													std::cout << posChi.name_ << " " << posChi.pos_ << std::endl;
										otherTempPos[posChi.name_] = posChi.pos_;
									}
									std::vector<double> fracs;
									for(const auto & name : names){
										fracs.emplace_back(finishedAmount[name]);
									}
									njh::randObjectGen gen(names,fracs);
									std::string templatedLandedOn = gen.genObj();
//												std::cout << "templatedLandedOn different? " << njh::colorBool(templatedLandedOn != chiProf.segsForSeqs_[finalPartial.first]->seqBase_.name_) << std::endl;
									if(templatedLandedOn != chiProf.segsForSeqs_[finalPartial.first]->seqBase_.name_){
										++chisGenerated;
										++threadChimeras[chiProf.segsForSeqs_[finalPartial.first]->seqBase_.seq_.substr(0, partLength) + chiProf.segsForSeqs_[templatedLandedOn]->seqBase_.seq_.substr(otherTempPos[templatedLandedOn])];
									}
								}
							}
							{
								std::lock_guard<std::mutex> lock(chimerasMut);
								for(const auto & chi : threadChimeras){
									finalChimeras[chi.first] += chi.second;
									for(uint32_t count = 0; count < chi.second; ++count){
										auto round = roundGen.genObj();
										++finalChimerasPerRound[round][chi.first];
									}
								}
							}
//										std::cout << "chisGenerated: " << chisGenerated << std::endl;
//										std::cout << "partials[threadNum]: " << partials[threadNum] << std::endl;
						};

						{
							if(numThreads <=1){
								genPartials(0);
							}else{
								std::vector<std::thread> threads;
								for(uint32_t t = 0; t < numThreads; ++t){
									threads.emplace_back(genPartials,t);
								}
								njh::concurrent::joinAllThreads(threads);
							}
						}
					}
					std::unordered_map<uint32_t, uint32_t> totalsForRounds;
					for(const auto & finalChimeraForRound : finalChimerasPerRound){
//						std::cout << "finalChimeraForRound.first: " << finalChimeraForRound.first << std::endl;
						uint32_t totalForthisRound = 0;
						for(const auto & finalChimera : finalChimeraForRound.second){
							bool matchesInputSeq = false;
							for(const auto & seqSampled : seqsSampled){
								if(seqSampled.seqBase_.seq_ == finalChimera.first){
									matchesInputSeq = true;
									break;
								}
							}
							if(!matchesInputSeq){
								totalForthisRound += finalChimera.second;
								bool matchesPreviousChimera = false;
								if(!matchesInputSeq){
									for(const auto & chimeraSeq : chimeraSeqs){
										if(chimeraSeq.seq_ == finalChimera.first){
											matchesPreviousChimera = true;
											break;
										}
									}
								}
								if(!matchesInputSeq && !matchesPreviousChimera){
									std::string name =njh::pasteAsStr("Chi.", chimeraSeqs.size(), "[PCRRound=", finalChimeraForRound.first, "]");
									auto backSeq = chimeraSeqs.emplace_back(seqInfo(name, finalChimera.first));
									backSeq.cnt_ = finalChimera.second;
								}
							}
						}
						totalsForRounds[finalChimeraForRound.first] = totalForthisRound;
//						std::cout << "\t" << "totalForthisRound: " << totalForthisRound << std::endl;
					}
					std::unordered_map<std::string, uint32_t> chiSeqPositions;
					for(const auto pos : iter::range(chimeraSeqs.size())){
						chiSeqPositions[chimeraSeqs[pos].seq_] = pos;
					}

					for(const auto & finalChimeraForRound : finalChimerasPerRound){
						std::unordered_map<std::string,std::unordered_map<std::string, uint64_t>> currentChiAllSeqCounts;
						std::unordered_map<std::string, std::string> barcodedSeqsCurrent;
						for(const auto & finalChi : finalChimeraForRound.second){
							bool matchesInputSeq = false;
							for(const auto & seqSampled : seqsSampled){
								if(seqSampled.seqBase_.seq_ == finalChi.first){
									matchesInputSeq = true;
									break;
								}
							}
							if(!matchesInputSeq){
								const auto & finalChimeraSeq = chimeraSeqs[chiSeqPositions[finalChi.first]];
								currentChiAllSeqCounts[finalChimeraSeq.name_][finalChimeraSeq.seq_] = finalChi.second;
								barcodedSeqsCurrent[finalChimeraSeq.name_] = finalChimeraSeq.seq_;
							}
						}
						if(!barcodedSeqsCurrent.empty()){
							auto sampleNumberAll = sampleReadsWithoutReplacementFinishPCR(barcodedSeqsCurrent,
									currentChiAllSeqCounts, totalsForRounds[finalChimeraForRound.first], libOutFile, seqFileLock,
									pcrRounds - finalChimeraForRound.first, numThreads);
							for(const auto & info : sampleNumberAll){
								ret.chimerasSampledForSequencing_[info.first].mutated_ += info.second.mutated_;
								ret.chimerasSampledForSequencing_[info.first].nonMutated_ += info.second.nonMutated_;
							}
						}
					}

					for(const auto & finalChimera : finalChimeras){
						bool matchesInputSeq = false;
						for(const auto & seqSampled : seqsSampled){
							if(seqSampled.seqBase_.seq_ == finalChimera.first){
								matchesInputSeq = true;
								break;
							}
						}
//						if(matchesInputSeq){
//							std::cout << njh::bashCT::green << std::endl;
//						}
//						std::cout << finalChimera.first << std::endl;
//						std::cout << finalChimera.second << std::endl;
//
//						std::cout << std::endl;
//						if(matchesInputSeq ){
//							std::cout << njh::bashCT::reset << std::endl;
//						}
						if(!matchesInputSeq){
							totalFinalChimeraCount+= finalChimera.second;
						}
					}
//					std::cout << "totalFinalChimera: " << totalFinalChimeraCount << " " << getPercentageString(totalFinalChimeraCount, finalReadAmount) << std::endl;
				}
				//can only sim chimeras if there are more than 1 sequence
			}
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
	auto sampleNumberAll = sampleReadsWithoutReplacementFinishPCR(barcodedSeqs,
			allSeqCounts, finalReadAmount - totalFinalChimeraCount, libOutFile, seqFileLock,
			pcrRounds - initialPcrRounds, numThreads);
	std::unordered_map<std::string, PCRSimulator::SimHapCounts::MutInfo> sampleNumber;
	for(const auto & samp : sampleNumberAll){
		if(!njh::beginsWith(samp.first, "Chi.")){
			sampleNumber.emplace(samp);
		}else{
			ret.chimerasSampledForSequencing_[samp.first].mutated_ += samp.second.mutated_;
			ret.chimerasSampledForSequencing_[samp.first].nonMutated_ += samp.second.nonMutated_;
		}
	}
	ret.genomesSampled_ = genomeCounts;
	ret.sampledForSequencing_ = sampleNumber;
	ret.chimeraSeqs_ = chimeraSeqs;

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

	numThreads = std::min<uint32_t>(numThreads, finalReadAmount);
	if(verbose_){
		auto names = getVectorOfMapKeys(seqs);
		std::cout << "Sampling from " << njh::conToStr(names,",") << " ..." << std::endl;
		std::cout << "numThreads: " << numThreads << std::endl;
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

//	std::cout << "sampled" << std::endl;
//	uint32_t totalSampled = 0;
//	for(const auto & sampled : allSampledSeqs){
//		std::cout << sampled.first << std::endl;
//		std::cout << "\t" << sampled.second.size() << std::endl;
//		totalSampled += sampled.second.size();
//	}
//	std::cout << "totalSampled:   " << totalSampled << std::endl;
//	std::cout << "finalReadAmount:" << finalReadAmount << std::endl;


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
					nameMeta.addMeta("threadNumber", threadNumber);
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

	if(numThreads <=1){
		finishPCR(0);
	}else{
		std::vector<std::thread> threads;
		for(uint32_t thread = 0; thread < numThreads; ++thread){
			threads.emplace_back(std::thread(finishPCR, thread));
		}

		for(auto & thread : threads){
			thread.join();
		}
	}
	return ret;
}


std::string PCRSimulator::runPcrSampleSingleTemplate(
		uint32_t roundsOfPcr,
		uint64_t randomNumberSelector,
		uint64_t randomNumberSelectorMax,
		std::string seq ){
	if(0 == roundsOfPcr){
		return seq;
	}
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


