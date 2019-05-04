/*
 * IlluminaRoughProfiler.cpp
 *
 *  Created on: May 2, 2019
 *      Author: nicholashathaway
 */



#include "IlluminaRoughProfiler.hpp"

namespace njhseq {

void IlluminaRoughProfiler::Counts::addOtherCounts(const Counts & other){
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
}

void IlluminaRoughProfiler::Counts::increaseCounts(const AlignmentResults & res){
	if(nullptr == res.alnSeqAligned_ || nullptr == res.refSeqAligned_){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error alignment sequences don't appear to be set" << "\n";
		throw std::runtime_error{ss.str()};
	}
	uint32_t alnSeqPosOffSet = 0;
	for(const auto pos : iter::range(len(*res.alnSeqAligned_))){
		if('-' == res.alnSeqAligned_->seq_[pos]){
			++alnSeqPosOffSet;
		}
		uint32_t realPosition = pos - alnSeqPosOffSet;
		if('-' != res.alnSeqAligned_->seq_[pos] && '-' != res.refSeqAligned_->seq_[pos]){
			++positionTotalCounts[realPosition];
			if(res.alnSeqAligned_->seq_[pos] != res.refSeqAligned_->seq_[pos]){
				//error
				++positionErrorCounts[realPosition];
				++qualErrorsCounts[realPosition][res.alnSeqAligned_->qual_[pos]];
				++baseChangeCounts[res.refSeqAligned_->seq_[pos]][res.alnSeqAligned_->seq_[pos]];
				++baseChangeCountsPerPosition[realPosition][res.refSeqAligned_->seq_[pos]][res.alnSeqAligned_->seq_[pos]];
			}else{
				//correct
				++qualCounts[realPosition][res.alnSeqAligned_->qual_[pos]];
			}
		}
	}
}

void IlluminaRoughProfiler::Counts::writeProfiles(const std::string & prefix, bool overWrite){
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



void IlluminaRoughProfiler::addOther(const IlluminaRoughProfiler & other){
	r1_counts.addOtherCounts(other.r1_counts);
	r2_counts.addOtherCounts(other.r2_counts);
}

void IlluminaRoughProfiler::increaseCounts(
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

}  // namespace njhseq


