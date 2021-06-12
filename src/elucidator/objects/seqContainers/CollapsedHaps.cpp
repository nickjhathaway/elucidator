/*
 * CollapsedHaps.cpp
 *
 *  Created on: Jun 12, 2021
 *      Author: nick
 */


#include "CollapsedHaps.hpp"

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

CollapsedHaps CollapsedHaps::readInReads(const SeqIOOptions & inOpts){
	CollapsedHaps ret;
	SeqInput reader(inOpts);
	reader.openIn();
	seqInfo seq;
	uint32_t seqCount = 0;
	std::unordered_set<std::string> allNames;

	while(reader.readNextRead(seq)) {
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


}  // namespace njhseq

