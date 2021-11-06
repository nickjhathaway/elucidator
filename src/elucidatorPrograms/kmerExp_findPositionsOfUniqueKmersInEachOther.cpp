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


#include <njhseq/concurrency/PairwisePairFactory.hpp>

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

class ChromVsChromGraph {
public:

	struct node{
		node(const Bed6RecordCore & chrom1Reg, const Bed6RecordCore & chrom2Reg):chrom1Reg_(chrom1Reg), chrom2Reg_(chrom2Reg){

		}
		Bed6RecordCore chrom1Reg_;
		Bed6RecordCore chrom2Reg_;

		std::vector<uint32_t> edges_;

		uint32_t group_{std::numeric_limits<uint32_t>::max()};
	};


	std::vector<node> nodes_;

	void connectNodes(uint32_t minDist){
		PairwisePairFactory pFac(nodes_.size());
		PairwisePairFactory::PairwisePair pair;

		std::function<bool(const node & n1, const node & n2)> nodesPass = [&minDist](const node & n1, const node & n2){
			if((n1.chrom1Reg_.overlaps(n2.chrom1Reg_, 1) && n1.chrom2Reg_.overlaps(n2.chrom2Reg_, 1)) ||
				 (n1.chrom1Reg_.chrom_ == n2.chrom1Reg_.chrom_ &&
				  n1.chrom2Reg_.chrom_ == n2.chrom2Reg_.chrom_ &&
					(uAbsdiff(n1.chrom1Reg_.chromStart_, n2.chrom1Reg_.chromEnd_) <= minDist ||
					 uAbsdiff(n1.chrom1Reg_.chromEnd_,   n2.chrom1Reg_.chromStart_) <= minDist) &&
					(uAbsdiff(n1.chrom2Reg_.chromStart_, n2.chrom2Reg_.chromEnd_) <= minDist ||
					 uAbsdiff(n1.chrom2Reg_.chromEnd_,   n2.chrom2Reg_.chromStart_) <= minDist)) ){
				return true;
			}
			return false;
		};

		while(pFac.setNextPair(pair)){
			if(nodesPass(nodes_[pair.col_], nodes_[pair.row_]) ){
				nodes_[pair.col_].edges_.emplace_back(pair.row_);
				nodes_[pair.row_].edges_.emplace_back(pair.col_);
			}
		}
	}

	void determineGroups(){
		uint32_t currentGroup = 0;
		for(const auto pos : iter::range(nodes_.size())){
			//skip if a group already assigned
			if(std::numeric_limits<uint32_t>::max() == nodes_[pos].group_){
				std::deque<uint32_t> toSpread;
				toSpread.emplace_back(pos);
				while(!toSpread.empty()){
					uint32_t nextPos = toSpread.front();
					toSpread.pop_front();
					nodes_[nextPos].group_ = currentGroup;
					for(const auto e : nodes_[nextPos].edges_){
						//spread if no group assigned
						if(std::numeric_limits<uint32_t>::max() == nodes_[e].group_){
							toSpread.emplace_back(e);
						}
					}
				}
				++currentGroup;
			}
		}


	}

	std::unordered_set<uint32_t> getGroupNumbers() const{
		std::unordered_set<uint32_t> groups;
		for(const auto & node : nodes_){
			groups.emplace(node.group_);
		}
		return groups;
	}
	struct GenomicRanges{
		GenomicRanges() {

		}
		std::string chrom1 = "";
		uint32_t minChrom1 = std::numeric_limits<uint32_t>::max();
		uint32_t maxChrom1 = 0;

		std::string chrom2 = "";
		uint32_t minChrom2 = std::numeric_limits<uint32_t>::max();
		uint32_t maxChrom2 = 0;

		uint32_t numberOfRevComp = 0;
		uint32_t totalSubRegions = 0;

		double getFracRevComp() const {
			if (totalSubRegions > 0) {
				return static_cast<double>(numberOfRevComp) / totalSubRegions;
			}
			return 0;
		}

		uint32_t chrom1Len() const {
			return maxChrom1 - minChrom1;
		}
		uint32_t chrom2Len() const {
			return maxChrom2 - minChrom2;
		}

		std::string genGroupName() const{
			return njh::pasteAsStr(
								chrom1, "-", minChrom1, "-", maxChrom1, getFracRevComp() > 0.5 ? "-rev": "-for",
								"--",
								chrom2, "-", minChrom2, "-", maxChrom2, "-for");
		}
	};

	std::unordered_map<uint32_t, GenomicRanges> genGroupRanges() const {
		std::unordered_map<uint32_t, GenomicRanges> ranges;
		for(const auto & node : nodes_){
			ranges[node.group_].chrom1 = node.chrom1Reg_.chrom_;
			ranges[node.group_].minChrom1 = std::min(node.chrom1Reg_.chromStart_, ranges[node.group_].minChrom1);
			ranges[node.group_].maxChrom1 = std::max(node.chrom1Reg_.chromEnd_,   ranges[node.group_].maxChrom1);

			ranges[node.group_].chrom2 = node.chrom2Reg_.chrom_;
			ranges[node.group_].minChrom2 = std::min(node.chrom2Reg_.chromStart_, ranges[node.group_].minChrom2);
			ranges[node.group_].maxChrom2 = std::max(node.chrom2Reg_.chromEnd_,   ranges[node.group_].maxChrom2);
			++ranges[node.group_].totalSubRegions;
			if(node.chrom1Reg_.reverseStrand()){
				++ranges[node.group_].numberOfRevComp;
			}
		}
		return ranges;
	}

	std::unordered_map<uint32_t, std::string> determineGroupNames() const {
		std::unordered_map<uint32_t, std::string> ret;


		auto ranges = genGroupRanges();
		for(const auto & range : ranges){
			ret[range.first] = range.second.genGroupName();
		}

		return ret;
	}

	std::vector<Bed6RecordCore> generateGroupRegions() const{
		std::vector<Bed6RecordCore> regions;
		auto ranges = genGroupRanges();
		for(const auto & range : ranges){
			std::string groupName = range.second.genGroupName();
			Bed6RecordCore region1(range.second.chrom1, range.second.minChrom1, range.second.maxChrom1, groupName, range.second.chrom1Len(),range.second.getFracRevComp() > 0.50 ? '-' : '+');
			Bed6RecordCore region2(range.second.chrom2, range.second.minChrom2, range.second.maxChrom2, groupName, range.second.chrom2Len(), '+');
			regions.emplace_back(region1);
			regions.emplace_back(region2);
		}
		return regions;
	}

	void sortByGroup() {
		njh::sort(nodes_, [](const node & n1, const node & n2){
			return n1.group_ < n2.group_;
		});
	}

private:
};

int kmerExpRunner::chromVsChromUniqueComp(const njh::progutils::CmdArgs & inputCommands){
	//Currently inefficient
	bfs::path genomeFnp = "";
	uint32_t numThreads = 1;
	bool checkComplement = false;
	uint32_t kmerLength = 31;
	uint32_t minSize = 100;
	uint32_t minDisForGrouping = 1000;
	OutOptions outOpts(bfs::path("out"), ".tab.txt");
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processWritingOptions(outOpts);
	setUp.setOption(minSize, "--minSize", "min Size");
	setUp.setOption(kmerLength, "--kmerLength", "kmer Length");
	setUp.setOption(checkComplement, "--checkComplement", "Check Complement");
	setUp.setOption(minDisForGrouping, "--minDisForGrouping", "min Dis For Grouping");

	setUp.setOption(numThreads, "--numThreads", "num Threads");
	setUp.setOption(genomeFnp, "--genome2bit", "genome 2bit", true);
	setUp.finishSetUp(std::cout);

	OutOptions outBedOpts(outOpts.outFilename_.string(), ".bed");
	OutOptions outGroupOpts(outOpts.outFilename_.string(), "_groupRegions.bed");


	OutputStream out(outOpts);
	outBedOpts.transferOverwriteOpts(outOpts);
	OutputStream outBed(outBedOpts);

	outGroupOpts.transferOverwriteOpts(outOpts);
	OutputStream outGroup(outGroupOpts);


	out << "chrom1\tchrom1Pos";
	out << "\tchrom2\tchrom2Pos\tsize\tchrom1RevComp\tgroup" << std::endl;
	std::mutex outMut;

	TwoBit::TwoBitFile topTReader(genomeFnp);
	VecStr chromNames = topTReader.sequenceNames();

	PairwisePairFactory pfactor(chromNames.size());

	std::function<void()> compSeqs = [&pfactor,&chromNames,&out,
																		&outMut,&kmerLength,
																		&checkComplement,&minSize,
																		&genomeFnp, &outBed,
																		&minDisForGrouping,
																		&outGroup](){

		PairwisePairFactory::PairwisePair pair;
		TwoBit::TwoBitFile tReader(genomeFnp);

		while(pfactor.setNextPair(pair)){

			auto seq1Chrom = chromNames[pair.row_];
			auto seq2Chrom = chromNames[pair.col_];

			ChromVsChromGraph graph;
			std::string seq1Seq;
			tReader[seq1Chrom]->getSequence(seq1Seq);

			std::string seq2Seq;
			tReader[seq2Chrom]->getSequence(seq2Seq);
			{
				kmerInfo kInfo1(seq1Seq, kmerLength, false);
				auto seq2Block = KmersSharedBlocks(seq2Chrom, seq2Seq, kmerLength);
				for (const auto pos : iter::range(seq1Seq.size() - kmerLength + 1)) {
					auto sub = seq1Seq.substr(pos, kmerLength);

					auto search = seq2Block.kInfo_->kmers_.find(sub);
					//if both sequences contain the kmer and the kmer is unique in both sequences, add them
					if (   search != seq2Block.kInfo_->kmers_.end()
							&& search->second.positions_.size() == 1
							&& kInfo1.kmers_[search->first].positions_.size() == 1) {
						seq2Block.addComp(kInfo1.kmers_[search->first].positions_.front(), search->second.positions_.front());
					}
				}
				//
				seq2Block.finish();
				auto positionsKeys = getVectorOfMapKeys(seq2Block.kComps_);
				njh::sort(positionsKeys);

				for(const auto & pos : positionsKeys){
					uint32_t actualSize = seq2Block.kComps_[pos].size_ + kmerLength - 1;
					if(actualSize > minSize){
						uint32_t chrom1End = seq2Block.kComps_[pos].refStart_ + actualSize;
						uint32_t chrom2End = seq2Block.kComps_[pos].start_ + actualSize;

						std::string chrom1OutName = njh::pasteAsStr(seq1Chrom,"-", seq2Block.kComps_[pos].refStart_, "-", chrom1End, "-for");
						std::string chrom2OutName = njh::pasteAsStr(seq2Chrom,"-", seq2Block.kComps_[pos].start_,    "-", chrom2End, "-for");

						std::string combinedName = njh::pasteAsStr(chrom1OutName, "--", chrom2OutName);

						Bed6RecordCore chrom1Reg(seq1Chrom, seq2Block.kComps_[pos].refStart_, chrom1End, chrom1OutName, actualSize, '+');
						Bed6RecordCore chrom2Reg(seq2Chrom, seq2Block.kComps_[pos].start_, chrom2End, chrom2OutName, actualSize, '+');
						graph.nodes_.emplace_back(ChromVsChromGraph::node(chrom1Reg, chrom2Reg));
					}
				}
			}
			if(checkComplement){
				seq1Seq = seqUtil::reverseComplement(seq1Seq, "DNA");
				kmerInfo kInfo1(seq1Seq, kmerLength, false);
				auto seq2Block = KmersSharedBlocks(seq2Chrom, seq2Seq, kmerLength);
				for (const auto pos : iter::range(seq1Seq.size() - kmerLength + 1)) {
					auto sub = seq1Seq.substr(pos, kmerLength);

					auto search = seq2Block.kInfo_->kmers_.find(sub);
					//if both sequences contain the kmer and the kmer is unique in both sequences, add them
					if (   search != seq2Block.kInfo_->kmers_.end()
							&& search->second.positions_.size() == 1
							&& kInfo1.kmers_[search->first].positions_.size() == 1) {
						seq2Block.addComp(kInfo1.kmers_[search->first].positions_.front(), search->second.positions_.front());
					}
				}
				//
				seq2Block.finish();
				auto positionsKeys = getVectorOfMapKeys(seq2Block.kComps_);
				njh::sort(positionsKeys);

				for(const auto & pos : positionsKeys){
					uint32_t actualSize = seq2Block.kComps_[pos].size_ + kmerLength - 1;
					if(actualSize > minSize){
						uint32_t chrom1End = seq1Seq.size() - seq2Block.kComps_[pos].refStart_;
						uint32_t chrom1Start = seq1Seq.size() - (seq2Block.kComps_[pos].refStart_ + actualSize);


						uint32_t chrom2End = seq2Block.kComps_[pos].start_ + actualSize;

						std::string chrom1OutName = njh::pasteAsStr(seq1Chrom,"-", chrom1Start,                      "-", chrom1End, "-rev");
						std::string chrom2OutName = njh::pasteAsStr(seq2Chrom,"-", seq2Block.kComps_[pos].start_,    "-", chrom2End, "-for");

						std::string combinedName = njh::pasteAsStr(chrom1OutName, "--", chrom2OutName);

						Bed6RecordCore chrom1Reg(seq1Chrom, chrom1Start, chrom1End, chrom1OutName, actualSize, '-');
						Bed6RecordCore chrom2Reg(seq2Chrom, seq2Block.kComps_[pos].start_, chrom2End, chrom2OutName, actualSize, '+');
						graph.nodes_.emplace_back(ChromVsChromGraph::node(chrom1Reg, chrom2Reg));
					}
				}
			}
			graph.connectNodes(minDisForGrouping);
			graph.determineGroups();
			graph.sortByGroup();
			auto groupNames = graph.determineGroupNames();
			auto groupRegions = graph.generateGroupRegions();
			{
				std::lock_guard<std::mutex> lock(outMut);
				for(const auto & node : graph.nodes_){

					outBed << node.chrom1Reg_.toDelimStr() << "\t" << groupNames[node.group_]<< std::endl;
					outBed << node.chrom2Reg_.toDelimStr() << "\t" << groupNames[node.group_]<< std::endl;
					out << node.chrom1Reg_.chrom_
							<< "\t" << node.chrom1Reg_.chromStart_
							<< "\t" << node.chrom2Reg_.chrom_
							<< "\t" << node.chrom2Reg_.chromStart_
							<< "\t" << node.chrom1Reg_.length()
							<< "\t" << njh::boolToStr('-' == node.chrom1Reg_.strand_)
							<< "\t" << groupNames[node.group_]
							<< std::endl;


				}
				for(const auto & groupRegion : groupRegions){
					outGroup << groupRegion.toDelimStrWithExtra() << std::endl;
				}
			}
		}
	};

	njh::concurrent::runVoidFunctionThreaded(compSeqs, numThreads);

	return 0;
}


int kmerExpRunner::findPositionsOfUniqueKmersInEachOther(const njh::progutils::CmdArgs & inputCommands){
	seqInfo seq1("seq1");
	seqInfo seq2("seq2");
	bool collapse = false;
	bool checkComplement = false;
	uint32_t kmerLength = 31;
	uint32_t minSize = 100;
	OutOptions outOpts(bfs::path("out"), ".tab.txt");
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processWritingOptions(outOpts);
	setUp.setOption(minSize, "--minSize", "min Size");
	setUp.setOption(kmerLength, "--kmerLength", "kmer Length");
	setUp.setOption(collapse, "--collapse", "collapse");
	setUp.setOption(checkComplement, "--checkComplement", "Check Complement");

	setUp.processSeq(seq1, "--seq1", "First sequence", true);
	setUp.processSeq(seq2, "--seq2", "Second sequence", true);
	setUp.finishSetUp(std::cout);

	OutputStream out(outOpts);




	if(collapse){
		out << "seq1\tseq1Pos";
		if(checkComplement){
			out << "\trevComp";
		}
		out << "\tseq2\tseq2Pos\tsize" << std::endl;

		{
			kmerInfo kInfo1(seq1.seq_, kmerLength, false);
			auto seq2Block = KmersSharedBlocks(seq2.name_, seq2.seq_, kmerLength);
			for (const auto pos : iter::range(seq1.seq_.size() - kmerLength + 1)) {
				auto sub = seq1.seq_.substr(pos, kmerLength);

				auto search = seq2Block.kInfo_->kmers_.find(sub);
				//if both sequences contain the kmer and the kmer is unique in both sequences, add them
				if (   search != seq2Block.kInfo_->kmers_.end()
						&& search->second.positions_.size() == 1
						&& kInfo1.kmers_[search->first].positions_.size() == 1) {
					seq2Block.addComp(kInfo1.kmers_[search->first].positions_.front(), search->second.positions_.front());
				}
			}
			//
			seq2Block.finish();
			auto positionsKeys = getVectorOfMapKeys(seq2Block.kComps_);
			njh::sort(positionsKeys);

			for(const auto & pos : positionsKeys){
				uint32_t actualSize = seq2Block.kComps_[pos].size_ + kmerLength - 1;
				if(actualSize > minSize){
					out << seq1.name_
					    << "\t" << seq2Block.kComps_[pos].refStart_;
					if(checkComplement){
						out << "\t" << njh::boolToStr(false);
					}
					out << "\t" << seq2.name_
							<< "\t" << seq2Block.kComps_[pos].start_
							<< "\t" << actualSize
							<< std::endl;
				}
			}
		}
		if(checkComplement){
			seq1.reverseComplementRead(false, true);
			kmerInfo kInfo1(seq1.seq_, kmerLength, false);
			auto seq2Block = KmersSharedBlocks(seq2.name_, seq2.seq_, kmerLength);
			for (const auto pos : iter::range(seq1.seq_.size() - kmerLength + 1)) {
				auto sub = seq1.seq_.substr(pos, kmerLength);

				auto search = seq2Block.kInfo_->kmers_.find(sub);
				//if both sequences contain the kmer and the kmer is unique in both sequences, add them
				if (   search != seq2Block.kInfo_->kmers_.end()
						&& search->second.positions_.size() == 1
						&& kInfo1.kmers_[search->first].positions_.size() == 1) {
					seq2Block.addComp(kInfo1.kmers_[search->first].positions_.front(), search->second.positions_.front());
				}
			}
			//
			seq2Block.finish();
			auto positionsKeys = getVectorOfMapKeys(seq2Block.kComps_);
			njh::sort(positionsKeys);

			for(const auto & pos : positionsKeys){
				uint32_t actualSize = seq2Block.kComps_[pos].size_ + kmerLength - 1;
				if(actualSize > minSize){
					out << seq1.name_
					    << "\t" << seq2Block.kComps_[pos].refStart_
							<< "\t" << njh::boolToStr(true)
							<< "\t" << seq2.name_
							<< "\t" << seq2Block.kComps_[pos].start_
							<< "\t" << actualSize
							<< std::endl;
				}
			}
		}

	}else{
		kmerInfo kInfo2(seq2.seq_, kmerLength, false);
		kmerInfo kInfo1(seq1.seq_, kmerLength, checkComplement);
		{
			std::unordered_set<std::string> sharedUniqueKmers;
			//


			out << "kmer\tseq1Pos\tseq2Pos";
			if(checkComplement){
				out << "\trevComp";
			}
			out << std::endl;
			VecStr sharedUniqueKmersVec(sharedUniqueKmers.begin(), sharedUniqueKmers.end());
			njh::sort(sharedUniqueKmersVec, [&kInfo1](const std::string & k1, const std::string & k2){
				return kInfo1.kmers_[k1].positions_.front() < kInfo1.kmers_[k2].positions_.front();
			});
			for(const auto & k : sharedUniqueKmersVec){
				out << k
						<< "\t" << kInfo1.kmers_[k].positions_.front()
						<< "\t" << kInfo2.kmers_[k].positions_.front();
				if(checkComplement){
					out << "\t" << njh::boolToStr(false);
				}
				out << std::endl;
			}
		}
		if(checkComplement){
			std::unordered_set<std::string> sharedUniqueKmers;
			//
			for(const auto & k :kInfo1.kmersRevComp_){
				if(1 == k.second.count_ && njh::in(k.first, kInfo2.kmers_) && 1 == kInfo2.kmers_[k.first].count_){
					sharedUniqueKmers.emplace(k.first);
				}
			}
			VecStr sharedUniqueKmersVec(sharedUniqueKmers.begin(), sharedUniqueKmers.end());
			njh::sort(sharedUniqueKmersVec, [&kInfo1](const std::string & k1, const std::string & k2){
				return kInfo1.kmersRevComp_[k1].positions_.front() < kInfo1.kmersRevComp_[k2].positions_.front();
			});
			for(const auto & k : sharedUniqueKmersVec){
				out << k
						<< "\t" << kInfo1.kmersRevComp_[k].positions_.front()
						<< "\t" << kInfo2.kmers_[k].positions_.front();
				out << "\t" << njh::boolToStr(true);
				out << std::endl;
			}
		}
	}


	return 0;
}

}  //namespace njhseq




