/*
 * CollapsedHaps.cpp
 *
 *  Created on: Jun 12, 2021
 *      Author: nick
 */


#include "CollapsedHaps.hpp"

#include <njhseq/concurrency/PairwisePairFactory.hpp>

namespace njhseq {


void CollapsedHaps::setSubNamesToMainSeqPos(){
	subNamesToMainSeqPos_.clear();
	for(const auto & names : iter::enumerate(names_)){
		for(const auto & name : names.second){
			subNamesToMainSeqPos_[name] = names.first;
		}
	}
}

void CollapsedHaps::revCompSeqs(){
	njh::for_each(seqs_, [](std::shared_ptr<seqInfo> & seq){ seq->reverseComplementRead(false, true);});
}

void CollapsedHaps::setFrequencies(uint32_t total){
	for(auto & seq : seqs_){
		seq->frac_ = seq->cnt_/total;
	}
}

void CollapsedHaps::setFrequencies(){
	auto total = getTotalHapCount();
	setFrequencies(total);
}

uint32_t CollapsedHaps::getTotalHapCount() const{
	uint32_t ret = 0;
	for(const auto & seq : seqs_){
		ret += seq->cnt_;
	}
	return ret;
}

uint32_t CollapsedHaps::getTotalUniqueHapCount() const{
	return seqs_.size();
}

std::vector<uint32_t> CollapsedHaps::getReadLenVec() const{
	std::vector<uint32_t> ret;
	for(const auto & seq : seqs_){
		for(uint32_t count = 0; count < seq->cnt_; ++count){
			ret.emplace_back(len(*seq));
		}
	}
	return ret;
}

std::unordered_map<uint32_t, uint32_t> CollapsedHaps::getReadLenMap() const{
	std::unordered_map<uint32_t, uint32_t> ret;
	for(const auto & seq : seqs_){
		ret[len(*seq)] += seq->cnt_;
	}
	return ret;
}

CollapsedHaps CollapsedHaps::readInReads(const SeqIOOptions & inOpts, std::unique_ptr<MultipleGroupMetaData> meta,
		std::unordered_map<std::string, std::string> metaValuesToAvoid){
	CollapsedHaps ret;
	SeqInput reader(inOpts);
	reader.openIn();
	seqInfo seq;
	uint32_t seqCount = 0;
	std::unordered_set<std::string> allNames;

	while(reader.readNextRead(seq)) {
		if(nullptr != meta) {
			meta->attemptToAddSeqMeta(seq);
		}
		//get meta keys if available
		if(MetaDataInName::nameHasMetaData(getSeqBase(seq).name_)){
			MetaDataInName metaData(getSeqBase(seq).name_);
			bool skip = false;
			for(const auto & ignoreField : metaValuesToAvoid){
				if(metaData.containsMeta(ignoreField.first) && metaData.getMeta(ignoreField.first) == ignoreField.second){
					skip = true;
					break;
					//skip this seq
				}
			}
			if(skip){
				continue;
			}
		}
		if(njh::in(seq.name_, allNames)){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "can't have seqs with the same name in input"<< "\n";
			ss << seq.name_ << " found more than once" << "\n";
			throw std::runtime_error{ss.str()};
		}
		allNames.emplace(seq.name_);
		readVec::handelLowerCaseBases(seq, inOpts.lowerCaseBases_);
		if(inOpts.removeGaps_){
			seq.removeGaps();
		}
		seqCount+= std::round(seq.cnt_);
		bool found = false;
		for (const auto & pos : iter::range(ret.seqs_.size())) {
			const auto & otherSeq = ret.seqs_[pos];
			if (otherSeq->seq_ == seq.seq_) {
				otherSeq->cnt_ += 1;
				ret.names_[pos].emplace(seq.name_);
				found = true;
				break;
			}
		}
		if (!found) {
			//since we are considering counts here as haplotypes, will make sure cnt_ is set to 1
			seq.cnt_ = 1;
			ret.seqs_.emplace_back(std::make_shared<seqInfo>(seq));
			ret.names_.emplace_back(std::unordered_set<std::string>{seq.name_});
		}
	}
	ret.setFrequencies();
	ret.setSubNamesToMainSeqPos();
	return ret;
}

std::vector<comparison> CollapsedHaps::getCompsAgainstRef(const seqInfo & refSeq, aligner & alignerObj, uint32_t numThreads) const {
	std::vector<comparison> ret(seqs_.size());
	concurrent::AlignerPool alnPool(alignerObj, numThreads);
	alnPool.initAligners();
	std::vector<uint32_t> positions;
	njh::concurrent::LockableQueue<uint32_t> posQueue(positions);

	std::mutex mut;
	std::function<void()> alignComp = [this,&posQueue,&refSeq,&alnPool,&ret,&mut,&alignerObj](){
		auto currentAligner = alnPool.popAligner();
		uint32_t pos = std::numeric_limits<uint32_t>::max();
		while(posQueue.getVal(pos)){
			currentAligner->alignCacheGlobal(refSeq, seqs_[pos]);
			currentAligner->profileAlignment(refSeq, seqs_[pos], false, false, false);
			ret[pos] = currentAligner->comp_;
		}
		{
			std::lock_guard<std::mutex> lock(mut);
			alignerObj.alnHolder_.mergeOtherHolder(currentAligner->alnHolder_);
		}
	};
	njh::concurrent::runVoidFunctionThreaded(alignComp, numThreads);
	return ret;
}

std::vector<std::vector<comparison>> CollapsedHaps::getPairwiseComps(aligner & alignerObj, uint32_t numThreads) const{
	PairwisePairFactory pFac(seqs_.size());
	std::vector<std::vector<comparison>> ret;
	for(const auto pos : iter::range(seqs_.size())){
		ret.emplace_back(std::vector<comparison>(pos));
	}
	concurrent::AlignerPool alnPool(alignerObj, numThreads);
	alnPool.initAligners();
	njh::ProgressBar pBar(pFac.totalCompares_);
	std::mutex mut;
	std::function<void()> alignComp = [this,&pFac,&alnPool,&ret,&mut,&alignerObj,&pBar](){
		auto currentAligner = alnPool.popAligner();
		PairwisePairFactory::PairwisePair pair;
		while(pFac.setNextPair(pair)){
			currentAligner->alignCacheGlobal(seqs_[pair.row_], seqs_[pair.col_]);
			currentAligner->profileAlignment(seqs_[pair.row_], seqs_[pair.col_], false, false, false);
			ret[pair.row_][pair.col_] = currentAligner->comp_;
			if(verbose_){
				pBar.outputProgAdd(std::cout, 1, true);
			}
		}
		{
			std::lock_guard<std::mutex> lock(mut);
			alignerObj.alnHolder_.mergeOtherHolder(currentAligner->alnHolder_);
		}
	};
	njh::concurrent::runVoidFunctionThreaded(alignComp, numThreads);
	return ret;
}

CollapsedHaps::AvgPairwiseMeasures CollapsedHaps::getAvgPairwiseMeasures(const std::vector<std::vector<comparison>> & allComps) const{
	AvgPairwiseMeasures ret;
	PairwisePairFactory pFac(seqs_.size());
	PairwisePairFactory::PairwisePair pair;
	while(pFac.setNextPair(pair)){
		uint32_t toalSeqInComp = seqs_[pair.row_]->cnt_ + seqs_[pair.col_]->cnt_;

		uint32_t totalBetweenComps =
				  PairwisePairFactory::getTotalPairwiseComps(toalSeqInComp)
				- PairwisePairFactory::getTotalPairwiseComps(seqs_[pair.row_]->cnt_)
				- PairwisePairFactory::getTotalPairwiseComps(seqs_[pair.col_]->cnt_);

		ret.avgPercentId  += allComps[pair.row_][pair.col_].distances_.eventBasedIdentityHq_ * totalBetweenComps;
		ret.avgNumOfDiffs += allComps[pair.row_][pair.col_].distances_.getNumOfEvents(true) * totalBetweenComps;
	}
	uint64_t total = getTotalHapCount();
	for(const auto pos : iter::range(getTotalUniqueHapCount() )){
		//add in the 100% percent matches between the same exact seqs for the pairwise comps
		ret.avgPercentId += PairwisePairFactory::getTotalPairwiseComps(seqs_[pos]->cnt_);
	}
	ret.avgPercentId  /= PairwisePairFactory::getTotalPairwiseComps(total);
	ret.avgNumOfDiffs /= PairwisePairFactory::getTotalPairwiseComps(total);
	return ret;
}


}  // namespace njhseq

