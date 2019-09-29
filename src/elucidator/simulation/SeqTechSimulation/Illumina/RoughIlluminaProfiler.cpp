/*
 * RoughIlluminaProfiler.cpp
 *
 *  Created on: May 2, 2019
 *      Author: nicholashathaway
 */



#include "RoughIlluminaProfiler.hpp"

namespace njhseq {

void RoughIlluminaProfiler::Counts::addOtherCounts(const Counts & other){
	for(const auto & pos : other.positionErrorCounts){
		positionErrorCounts[pos.first] += pos.second;
	}
	for(const auto & pos : other.positionTotalCounts){
		positionTotalCounts[pos.first] += pos.second;
	}

	for(const auto & refBase : other.baseChangeCounts){
		for(const auto & seqBase : refBase.second){
			baseChangeCounts[refBase.first][seqBase.first] += seqBase.second;
		}
	}

	for(const auto & pos : other.baseChangeCountsPerPosition){
		for(const auto & refBase : pos.second){
			for(const auto & seqBase : refBase.second){
				baseChangeCountsPerPosition[pos.first][refBase.first][seqBase.first] += seqBase.second;
			}
		}
	}

	for(const auto & pos : other.qualCounts){
		for(const auto & qual : pos.second){
			qualCounts[pos.first][qual.first] += qual.second;
		}
	}

	for(const auto & pos : other.qualErrorsCounts){
		for(const auto & qual : pos.second){
			qualErrorsCounts[pos.first][qual.first] += qual.second;
		}
	}
	for(const auto & pos : other.deletions_){
		addOtherVec(deletions_[pos.first], pos.second);
	}

	for(const auto & pos : other.insertions_){
		addOtherVec(insertions_[pos.first], pos.second);
	}
	perfectHits_ += other.perfectHits_;
	addOtherVec(percentIds_, other.percentIds_);
}

void RoughIlluminaProfiler::Counts::increaseCounts(const seqInfo & refAln,
		const seqInfo & queryAln,
		const comparison & comp){

	if(refAln.seq_.size() != queryAln.seq_.size()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error ref alignment seq size, " << refAln.seq_.size() << " does not match query aln size, " <<  queryAln.seq_.size()<< "\n";
		throw std::runtime_error{ss.str()};
	}
	percentIds_.emplace_back(comp.distances_.eventBasedIdentity_);
	if(comp.distances_.eventBasedIdentity_ >= 1){
		++perfectHits_;
	}
	uint32_t alnSeqPosOffSet = 0;
	for(const auto pos : iter::range(queryAln.seq_.size())){
		if('-' == queryAln.seq_[pos]){
			++alnSeqPosOffSet;
		}
		uint32_t realPosition = pos - alnSeqPosOffSet;
		if('-' != queryAln.seq_[pos] && '-' != refAln.seq_[pos]){
			++positionTotalCounts[realPosition];
			if(queryAln.seq_[pos] != refAln.seq_[pos]){
				//error
				++positionErrorCounts[realPosition];
				++qualErrorsCounts[realPosition][queryAln.qual_[pos]];
				++baseChangeCounts[refAln.seq_[pos]][queryAln.seq_[pos]];
				++baseChangeCountsPerPosition[realPosition][refAln.seq_[pos]][queryAln.seq_[pos]];
			}else{
				//correct
				++qualCounts[realPosition][queryAln.qual_[pos]];
			}
		}
	}
	auto isHomopolymer =
			[](const std::string & k) {
				return std::all_of(k.begin(), k.end(),[&k](const char c) {return k.front() == c;});
			};
	for(const auto & gap : comp.distances_.alignmentGaps_){
		std::string refHomopolymer;
		if(isHomopolymer(gap.second.gapedSequence_)){
			{//search back
				auto seqAlnPosition = gap.first;
				while(seqAlnPosition > 0){
					--seqAlnPosition;
					if(gap.second.gapedSequence_.front() == queryAln.seq_[seqAlnPosition] &&
						 gap.second.gapedSequence_.front() == refAln.seq_[seqAlnPosition]){
						refHomopolymer.push_back(gap.second.gapedSequence_.front());
					}else{
						break;
					}
				}
			}
			{//search forward
				auto seqAlnPosition = gap.first + gap.second.gapedSequence_.size();
				while(seqAlnPosition < (queryAln.seq_.size() - 1) ){
					++seqAlnPosition;
					if(gap.second.gapedSequence_.front() == queryAln.seq_[seqAlnPosition] &&
						 gap.second.gapedSequence_.front() == refAln.seq_[seqAlnPosition]){
						refHomopolymer.push_back(gap.second.gapedSequence_.front());
					}else{
						break;
					}
				}
			}
		}
		if(gap.second.ref_){
			insertions_[gap.second.seqPos_].emplace_back(Indel(gap.second, refHomopolymer));
		}else{
			deletions_[gap.second.seqPos_].emplace_back(Indel(gap.second, refHomopolymer));
		}
	}
}



void RoughIlluminaProfiler::Counts::increaseCounts(const ReAlignedSeq & res){
	if(getSoftClipAmount(res.bAln_) < softClipCutOff_){
		increaseCounts(res.alnRefSeq_, res.alnQuerySeq_, res.comp_);
	}
}

void RoughIlluminaProfiler::Counts::increaseCounts(const AlignmentResults & res){
	if(nullptr == res.alnSeqAligned_ || nullptr == res.refSeqAligned_){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error alignment sequences don't appear to be set" << "\n";
		throw std::runtime_error{ss.str()};
	}
	if(getSoftClipAmount(res.bAln_) < softClipCutOff_){
		increaseCounts(*res.alnSeqAligned_, *res.refSeqAligned_, res.comp_);
	}
}




void RoughIlluminaProfiler::Counts::writeIndels(const std::string & prefix, bool overWrite){
	//positional error rate
	OutOptions positional_error_rateOpts(bfs::path(prefix + "_indels.tab.txt"));
	positional_error_rateOpts.overWriteFile_ = overWrite;
	OutputStream positional_error_rateOut(positional_error_rateOpts);
	positional_error_rateOut << "position\tindel\tgapSeq\trefHomopolymer" << std::endl;
	{
		auto posKeys = getVectorOfMapKeys(deletions_);
		njh::sort(posKeys);
		for(const auto & pos : posKeys){
			for(const auto & del : deletions_[pos]){
				positional_error_rateOut << pos
						<< "\t" << "deletion"
						<< "\t" << del.gapinfo_.gapedSequence_
						<< "\t" << del.refHomopolymer_ << std::endl;
			}
		}
	}

	{
		auto posKeys = getVectorOfMapKeys(insertions_);
		njh::sort(posKeys);
		for(const auto & pos : posKeys){
			for(const auto & ins : insertions_[pos]){
				positional_error_rateOut << pos
						<< "\t" << "insertion"
						<< "\t" << ins.gapinfo_.gapedSequence_
						<< "\t" << ins.refHomopolymer_ << std::endl;
			}
		}
	}
}

void RoughIlluminaProfiler::Counts::writeProfiles(const std::string & prefix, bool overWrite){

	//scores
	OutOptions scoresOpts(bfs::path(prefix + "_scores.tab.txt"));
	scoresOpts.overWriteFile_ = overWrite;
	OutputStream scoresOut(scoresOpts);
	scoresOut << "score" << std::endl;
	{
		scoresOut << njh::conToStr(percentIds_, "\n") << std::endl;
	}
	//percent hits
	OutOptions hitsOpts(bfs::path(prefix + "_hits.tab.txt"));
	hitsOpts.overWriteFile_ = overWrite;
	OutputStream hitsOut(hitsOpts);
	hitsOut << "totalReads\tperfectHits\tperfectHitsRatio" << std::endl;
	{
		hitsOut << percentIds_.size()
				<< "\t" << perfectHits_
				<< "\t" << perfectHits_/static_cast<double>(percentIds_.size()) << std::endl;
	}


	//positional error rate
	OutOptions positional_error_rateOpts(bfs::path(prefix + "_positional_error_rate.tab.txt"));
	positional_error_rateOpts.overWriteFile_ = overWrite;
	OutputStream positional_error_rateOut(positional_error_rateOpts);
	positional_error_rateOut << "position\terrors\ttotal\trate" << std::endl;
	{
		auto posKeys = getVectorOfMapKeys(positionTotalCounts);
		njh::sort(posKeys);
		for(const auto & pos : posKeys){
			positional_error_rateOut << pos
					<< "\t" << positionErrorCounts[pos]
					<< "\t" << positionTotalCounts[pos]
					<< "\t" << positionErrorCounts[pos]/static_cast<double>(positionTotalCounts[pos]) << std::endl;
		}
	}


	//base substitution rates
	OutOptions base_substitution_ratesOpts(bfs::path(prefix + "_base_substitution_rates.tab.txt"));
	base_substitution_ratesOpts.overWriteFile_ = overWrite;
	OutputStream base_substitution_ratesOut(base_substitution_ratesOpts);
	base_substitution_ratesOut << "ref\tseq\tcount" << std::endl;
	std::vector<char> bases{'A', 'C', 'G', 'T'};
	for(const auto & refBase : bases){
		double total = 0;
		for(const auto & seqBase : bases){
			if(seqBase == refBase){
				continue;
			}
			total += baseChangeCounts[refBase][seqBase];
		}
		for(const auto & seqBase : bases){
			if(seqBase == refBase){
				continue;
			}
			base_substitution_ratesOut << refBase
				 << "\t" << seqBase
				 << "\t" << baseChangeCounts[refBase][seqBase] << std::endl;
		}
	}

	//base substitution rates per position
	OutOptions base_substitution_rates_perPositionOpts(bfs::path(prefix + "_base_substitution_rates_perPosition.tab.txt"));
	base_substitution_ratesOpts.overWriteFile_ = overWrite;
	OutputStream base_substitution_rates_perPositionOut(base_substitution_rates_perPositionOpts);
	base_substitution_rates_perPositionOut << "pos\tref\tseq\tcount" << std::endl;
	auto allPositions = getVectorOfMapKeys(baseChangeCountsPerPosition);
	njh::sort(allPositions);
	for(const auto & pos :  allPositions){
		for(const auto & refBase : bases){
			double total = 0;
			for(const auto & seqBase : bases){
				if(seqBase == refBase){
					continue;
				}
				total += baseChangeCountsPerPosition[pos][refBase][seqBase];
			}
			for(const auto & seqBase : bases){
				if(seqBase == refBase){
					continue;
				}
				base_substitution_rates_perPositionOut << pos
				   << "\t" << refBase
					 << "\t" << seqBase
					 << "\t" << baseChangeCountsPerPosition[pos][refBase][seqBase] << std::endl;
			}
		}
	}


	//quality distribution for errors
	OutOptions quality_distribution_for_error_callsOpts(bfs::path(prefix + "_quality_distribution_for_error_calls.tab.txt"));
	quality_distribution_for_error_callsOpts.overWriteFile_ = overWrite;
	OutputStream quality_distribution_for_error_callsOut(quality_distribution_for_error_callsOpts);
	quality_distribution_for_error_callsOut << "position\tquality\tcount" << std::endl;
	{
		auto posKeys = getVectorOfMapKeys(qualErrorsCounts);
		njh::sort(posKeys);
		for(const auto pos : posKeys){
			auto qualKeys = getVectorOfMapKeys(qualErrorsCounts[pos]);
			njh::sort(qualKeys);
			for(const auto qual : qualKeys){
				quality_distribution_for_error_callsOut << pos
						<< "\t" << qual
						<< "\t" << qualErrorsCounts[pos][qual] << std::endl;
			}
		}
	}


	//quality distribution for correct calls
	OutOptions quality_distribution_for_correct_callsOpts(bfs::path(prefix + "_quality_distribution_for_correct_calls.tab.txt"));
	quality_distribution_for_correct_callsOpts.overWriteFile_ = overWrite;
	OutputStream quality_distribution_for_correct_callsOut(quality_distribution_for_correct_callsOpts);
	quality_distribution_for_correct_callsOut << "position\tquality\tcount" << std::endl;
	{
		auto posKeys = getVectorOfMapKeys(qualCounts);
		njh::sort(posKeys);
		for(const auto pos : posKeys){
			auto qualKeys = getVectorOfMapKeys(qualCounts[pos]);
			njh::sort(qualKeys);
			for(const auto qual : qualKeys){
				quality_distribution_for_correct_callsOut << pos
						<< "\t" << qual
						<< "\t" << qualCounts[pos][qual] << std::endl;
			}
		}
	}
}



void RoughIlluminaProfiler::addOther(const RoughIlluminaProfiler & other){
	r1_counts.addOtherCounts(other.r1_counts);
	r2_counts.addOtherCounts(other.r2_counts);

}

void RoughIlluminaProfiler::increaseCounts(
		const BamTools::BamAlignment & bAln,
		const BamTools::RefVector & refData,
		TwoBit::TwoBitFile & tReader) {
	AlignmentResults res(bAln, refData);
	res.setRefSeq(tReader);
	res.setAlignedObjects();
	if(bAln.IsFirstMate()){
		r1_counts.increaseCounts(res);
	}else{
		r2_counts.increaseCounts(res);
	}
}

void RoughIlluminaProfiler::increaseCounts(const ReAlignedSeq & res) {
	if (res.bAln_.IsFirstMate()) {
		r1_counts.increaseCounts(res);
	} else {
		r2_counts.increaseCounts(res);
	}
}




}  // namespace njhseq


