/*
 * GenomeSeqSearch.cpp
 *
 *  Created on: Aug 26, 2017
 *      Author: nick
 */


#include "GenomeSeqSearch.hpp"

namespace njhseq {

GenomeSeqSearch::SeqStreak::SeqStreak(size_t seqPos, const GenomicRegion & region) :
		seqPos_(seqPos), region_(region) {
	region_.uid_ = estd::to_string(seqPos_);
}



bool GenomeSeqSearch::SeqStreak::isNextPos(const std::string & chrom, const size_t gPos,
		const bool isRevCompl, const size_t seqPos, const uint32_t kLen) {
	return region_.chrom_ == chrom
			&& region_.end_ + 1 == (gPos + kLen)
			&& region_.reverseSrand_ == isRevCompl
			&& seqPos_ + gPos - region_.start_  == seqPos;
}

void GenomeSeqSearch::SeqStreak::increaseEndPos(uint32_t size) {
	growing_ = true;
	region_.end_ += size;
}




GenomeSeqSearch::GenomeSeqSearch(const OutOptions & outOpts) {
	out_ = outOpts.openFile();
}


void GenomeSeqSearch::growStreaks(const std::string & kmer, size_t seqPos, bool isRevComp, const std::unordered_map<std::string, std::unique_ptr<kmerInfo>> & infos ){
	auto isHomopolymer =
			[](const std::string & k) {
				return std::all_of(k.begin(), k.end(),[&k](const char c) {return k.front() == c;});
			};

	std::vector<SeqStreak> adding;
	//set all growing streaks to false until a streak is found
	for(auto & streak : growingStreaks_){
		streak.growing_ = false;
	}
	bool onlyGrow = false;
	if(isHomopolymer(kmer)){
		//don't count strings of all 1 letter as this will have a massive amount of occurrences;
		onlyGrow = true;
	}
	uint32_t totalCount = 0;
	for(const auto & chrom : infos){
		if(chrom.second->kmers_.end() != chrom.second->kmers_.find(kmer)){
			totalCount += chrom.second->kmers_.at(kmer).count_;
		}
	}
//		std::cout << std::endl;
//		std::cout << totalCount << std::endl;
	if(totalCount > countCapForGrowing_){
		onlyGrow = true;
	}
	for(const auto & chrom : infos){
		if(chrom.second->kmers_.end() != chrom.second->kmers_.find(kmer)){
			for(const auto & gPos : chrom.second->kmers_.at(kmer).positions_){
				bool foundStreak = false;
				for(auto & streak : growingStreaks_){
					if(streak.isNextPos(chrom.first, gPos, isRevComp, seqPos, kmer.size())){
						foundStreak = true;
						streak.increaseEndPos();
					}
				}
				if(!foundStreak && !onlyGrow){
					adding.emplace_back(SeqStreak(seqPos, GenomicRegion("", chrom.first, gPos, gPos + kmer.size(), isRevComp)));
				}
			}
		}
	}
	addOtherVec(growingStreaks_, adding);
}

void GenomeSeqSearch::pruneStreaks(uint32_t minSize){
	std::lock_guard<std::mutex> lock(mut_);
	std::vector<size_t> positions;
	for(const auto pos : iter::range(growingStreaks_.size())){
		if(!growingStreaks_[pos].growing_){
			positions.emplace_back(pos);
		}
	}
	//sort backwards so not to invalidate positions while erasing
	std::sort(positions.rbegin(), positions.rend());
	for(const auto pos : positions){
		//write out streaks if they are longer than the minSize
		if(growingStreaks_[pos].region_.getLen() >= minSize){
			(*out_) << growingStreaks_[pos].region_.genBedRecordCore().toDelimStr() << std::endl;
		}
		growingStreaks_.erase(growingStreaks_.begin() + pos);
	}
}

void GenomeSeqSearch::purgeAllStreaks(uint32_t minSize){
	for(const auto& streak : growingStreaks_){
		//write out streaks if they are longer than the minSize
		if(streak.region_.getLen() >= minSize){
			(*out_) << streak.region_.genBedRecordCore().toDelimStr() << std::endl;
		}
	}
	growingStreaks_.clear();
}


}  // namespace njhseq

