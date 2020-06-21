/*
 * kmerExp_findPositionsOfUniqueKmersInEachOther.cpp
 *
 *  Created on: Jun 17, 2020
 *      Author: nick
 */


// elucidator - A library for analyzing sequence data
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
//
// This file is part of elucidator.
//
// elucidator is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// elucidator is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with elucidator.  If not, see <http://www.gnu.org/licenses/>.
//

#include "kmerExp.hpp"
#include "elucidator/utils/KmerUtils.hpp"
#include "elucidator/objects/dataContainers.h"
#include "elucidator/simulation.h"
#include "elucidator/objects/seqObjects/seqKmers.h"

#include "elucidator/objects/MiscUtility/GenomeSeqSearch.hpp"


namespace njhseq {


int kmerExpRunner::kmerCompareTwoSetsOfContigs(const njh::progutils::CmdArgs & inputCommands){

	bfs::path refContigs = "";
	uint32_t kmerLength = 31;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(kmerLength, "--kmerLength", "kmer Length");
	uint32_t minSegmentLen = kmerLength + 1;
	setUp.setOption(minSegmentLen, "--minSegmentLen", "min Segment Len");

	setUp.processReadInNames(true);
	setUp.processRefFilename(true);
	setUp.processDirectoryOutputName(true);
	setUp.finishSetUp(std::cout);

	setUp.startARunLog(setUp.pars_.directoryName_);
	auto inputSeqs = SeqInput::getSeqVec<seqInfo>(setUp.pars_.ioOptions_);
	auto refSeqs = SeqInput::getSeqVec<seqInfo>(setUp.pars_.refIoOptions_);


	std::unordered_map<std::string, kmerInfo> inputSeqKmerInfos;
	std::unordered_map<std::string, kmerInfo> refSeqKmerInfos;

	for(const auto & seq : inputSeqs){
		if(len(seq) >= kmerLength){
			inputSeqKmerInfos.emplace(seq.name_, kmerInfo(seq.seq_, kmerLength, false));
		}
	}

	for(const auto & seq : refSeqs){
		if(len(seq) >= kmerLength){
			refSeqKmerInfos.emplace(seq.name_, kmerInfo(seq.seq_, kmerLength, false));
		}
	}


	struct GrowingSegment{
		GrowingSegment(const Bed3RecordCore & queryRegion, const Bed3RecordCore & refRegion):queryRegion_(queryRegion), refRegion_(refRegion){}
		Bed3RecordCore queryRegion_;
		Bed3RecordCore refRegion_;
		bool grewThisRound_{true};
	};


	{
		OutputStream inputSeqsLensOut(njh::files::make_path(setUp.pars_.directoryName_, "inputSeqLengths.tab.txt"));
		inputSeqsLensOut << "name\tlength" << std::endl;


		OutputStream inputSharedRegionsOut(njh::files::make_path(setUp.pars_.directoryName_, "sharedRegions.tab.txt"));
		inputSharedRegionsOut << "#name\tstart\tend\trefName\trefStart\trefEnd\tlength" << std::endl;


		for(const auto & seq : inputSeqs){
			inputSeqsLensOut << seq.name_ << '\t' << len(seq) << std::endl;
			if(len(seq) > kmerLength){
				std::vector<GrowingSegment> growingSegments;
				for(const auto pos : iter::range(len(seq) + 1 - kmerLength)){
					for(auto & grow : growingSegments){
						grow.grewThisRound_  = false;
					}
					std::string k = seq.seq_.substr(pos, kmerLength);
					std::vector<Bed3RecordCore> kPositions;
					for(const auto & refSeqKmerInfo : refSeqKmerInfos){
						if(njh::in(k, refSeqKmerInfo.second.kmers_)){
							for(const auto & kPos : refSeqKmerInfo.second.kmers_.at(k).positions_){
								kPositions.emplace_back(Bed3RecordCore(refSeqKmerInfo.first, kPos, kPos + kmerLength));
							}
						}
					}
					for(const auto & kPos : kPositions){
						bool found = false;
						for( auto & grow : growingSegments){
							if(grow.refRegion_.chrom_ == kPos.chrom_ && grow.refRegion_.chromEnd_ + 1 == kPos.chromEnd_){
								found = true;
								++grow.refRegion_.chromEnd_;
								++grow.queryRegion_.chromEnd_;
								grow.grewThisRound_ = true;
							}
						}
						if(!found){
							growingSegments.emplace_back(Bed3RecordCore(seq.name_, pos, pos + kmerLength), kPos);
						}
					}
					//prune
					std::vector<uint32_t> toErase;
					for(const auto gPos : iter::range(growingSegments.size())){
						if(!growingSegments[gPos].grewThisRound_){
							if(growingSegments[gPos].queryRegion_.length() >= minSegmentLen){
								inputSharedRegionsOut << growingSegments[gPos].queryRegion_.chrom_
										<< "\t" << growingSegments[gPos].queryRegion_.chromStart_
										<< "\t" << growingSegments[gPos].queryRegion_.chromEnd_
										<< "\t" << growingSegments[gPos].refRegion_.chrom_
										<< "\t" << growingSegments[gPos].refRegion_.chromStart_
										<< "\t" << growingSegments[gPos].refRegion_.chromEnd_
										<< "\t" << growingSegments[gPos].refRegion_.length() << std::endl;
							}
							toErase.emplace_back(gPos);
						}
					}
					for(const auto & toErasePos : iter::reversed(toErase)){
						growingSegments.erase(growingSegments.begin() + toErasePos);
					}
				}
				//prune
				for(auto & grow : growingSegments){
					grow.grewThisRound_  = false;
				}
				std::vector<uint32_t> toErase;
				for(const auto gPos : iter::range(growingSegments.size())){
					if(!growingSegments[gPos].grewThisRound_){
						if(growingSegments[gPos].queryRegion_.length() >= minSegmentLen){
							inputSharedRegionsOut << growingSegments[gPos].queryRegion_.chrom_
									<< "\t" << growingSegments[gPos].queryRegion_.chromStart_
									<< "\t" << growingSegments[gPos].queryRegion_.chromEnd_
									<< "\t" << growingSegments[gPos].refRegion_.chrom_
									<< "\t" << growingSegments[gPos].refRegion_.chromStart_
									<< "\t" << growingSegments[gPos].refRegion_.chromEnd_
									<< "\t" << growingSegments[gPos].refRegion_.length() << std::endl;
						}
						toErase.emplace_back(gPos);
					}
				}
				for(const auto & toErasePos : iter::reversed(toErase)){
					growingSegments.erase(growingSegments.begin() + toErasePos);
				}
			}
		}
	}
	{
		OutputStream refSeqsLensOut(njh::files::make_path(setUp.pars_.directoryName_, "refSeqLengths.tab.txt"));
		refSeqsLensOut << "name\tlength" << std::endl;
		for(const auto & seq : refSeqs){
			refSeqsLensOut << seq.name_ << '\t' << len(seq) << std::endl;
		}
	}
	return 0;
}





int kmerExpRunner::findPositionsOfUniqueKmersInEachOther(const njh::progutils::CmdArgs & inputCommands){
	seqInfo seq1("seq1");
	seqInfo seq2("seq2");
	uint32_t kmerLength = 7;
	OutOptions outOpts(bfs::path("out.tab.txt"));
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processWritingOptions(outOpts);
	setUp.setOption(kmerLength, "--kmerLength", "kmer Length");
	setUp.processSeq(seq1, "--seq1", "First sequence", true);
	setUp.processSeq(seq2, "--seq2", "Second sequence", true);
	setUp.finishSetUp(std::cout);

	OutputStream out(outOpts);


	kmerInfo kInfo1(seq1.seq_, kmerLength, false);
	kmerInfo kInfo2(seq2.seq_, kmerLength, false);

	std::unordered_set<std::string> sharedUniqueKmers;

	for(const auto & k :kInfo1.kmers_){
		if(1 == k.second.count_ && njh::in(k.first, kInfo2.kmers_) && 1 == kInfo2.kmers_[k.first].count_){
			sharedUniqueKmers.emplace(k.first);
		}
	}

	out << "kmer\tseq1Pos\tseq2Pos" << std::endl;
	VecStr sharedUniqueKmersVec(sharedUniqueKmers.begin(), sharedUniqueKmers.end());
	njh::sort(sharedUniqueKmersVec, [&kInfo1](const std::string & k1, const std::string & k2){
		return kInfo1.kmers_[k1].positions_.front() < kInfo1.kmers_[k2].positions_.front();
	});
	for(const auto & k : sharedUniqueKmersVec){
		out << k
				<< "\t" << kInfo1.kmers_[k].positions_.front()
				<< "\t" << kInfo2.kmers_[k].positions_.front() << std::endl;
	}


	return 0;
}

}  //namespace njhseq




