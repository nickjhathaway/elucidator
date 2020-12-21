/*
 * createSharedPathwaysFromContigs.cpp
 *
 *  Created on: Nov 5, 2019
 *      Author: nicholashathaway
 */



#include <njhseq/GenomeUtils.h>

#include "miscRunner.hpp"
#include "elucidator/simulation.h"
#include "elucidator/BamToolsUtils.h"
#include "elucidator/objects/dataContainers/graphs/ContigsCompareGraph.hpp"
#include "elucidator/BioRecordsUtils/BedUtility.hpp"


#include "elucidator/PopulationGenetics.h"

namespace njhseq {




int miscRunner::createSharedPathwaysFromReads(const njh::progutils::CmdArgs & inputCommands){

	std::vector<std::string> contigFnps;
	std::string contigFnpsStr;
	uint32_t klen = 25;
	uint32_t kmerOccurenceCutOff = 0;
	bool writeBasesShared = false;
	seqSetUp setUp(inputCommands);
	setUp.setOption(contigFnpsStr, "--contigFnps", "contigFnps", true);
	setUp.setOption(klen, "--klen", "k-mer length", true);
	setUp.setOption(writeBasesShared, "--writeBasesShared", "Write Bases Shared");
	setUp.setOption(kmerOccurenceCutOff, "--kmerOccurenceCutOff", "K-mer Occurence Cut Off");

	setUp.processDirectoryOutputName("createSharedPathwaysFromContigs_TODAY", true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	OutOptions outOpts(njh::files::make_path(setUp.pars_.directoryName_, "graph"), ".dot");

	contigFnps = tokenizeString(contigFnpsStr, ",");
	std::vector<std::vector<seqInfo>> inputs;
	for(const auto & contigFnp : contigFnps){
		inputs.emplace_back(SeqInput::getSeqVec<seqInfo>(SeqIOOptions::genFastaIn(contigFnp)));
	}



	ContigsCompareGraphDev compGraph(klen);
	compGraph.setOccurenceCutOff(kmerOccurenceCutOff);
	for(const auto & contigs : inputs){
		for(const auto & contig : contigs){
			compGraph.increaseKCounts(contig.seq_);
		}
	}
	compGraph.populateNodesFromCounts();
	for(const auto & contigs : inputs){
		for(const auto & contig : contigs){
			compGraph.threadThroughSequence(contig, contig.name_);
		}
	}
	compGraph.collapseSingleLinkedPaths();
	uint32_t graphCount = 0;
	{
		auto outOptsCurrent = outOpts;
		outOptsCurrent.outFilename_ = njh::files::prependFileBasename(outOptsCurrent.outFilename_,
				njh::leftPadNumStr<uint32_t>(graphCount++, 10000) + "_");
		OutputStream out(outOptsCurrent);
		compGraph.writeRectangleDotColorBySampleCount(out);
		auto outSeqs = compGraph.nodesToSeqs();
		auto seqsOutOpts = SeqIOOptions::genFastaOut(outOptsCurrent.outFilename_);
		seqsOutOpts.out_.transferOverwriteOpts(outOpts);
		SeqOutput::write(outSeqs, seqsOutOpts);
	}

	uint32_t splitCount = 0;
	while(compGraph.splitNodesWithRedundantKmers()){
		{
			++splitCount;
			auto outOptsCurrent = outOpts;
			outOptsCurrent.outFilename_ = njh::files::prependFileBasename(outOptsCurrent.outFilename_,
					njh::pasteAsStr(njh::leftPadNumStr<uint32_t>(graphCount++, 10000), "_splitCount_", splitCount, "_"));
			OutputStream out(outOptsCurrent);
			compGraph.writeRectangleDotColorBySampleCount(out);
			auto outSeqs = compGraph.nodesToSeqs();
			auto seqsOutOpts = SeqIOOptions::genFastaOut(outOptsCurrent.outFilename_);
			seqsOutOpts.out_.transferOverwriteOpts(outOpts);
			SeqOutput::write(outSeqs, seqsOutOpts);
		}
		compGraph.collapseSingleLinkedPaths();

		{
			auto outOptsCurrent = outOpts;
			outOptsCurrent.outFilename_ = njh::files::prependFileBasename(outOptsCurrent.outFilename_,
					njh::pasteAsStr(njh::leftPadNumStr<uint32_t>(graphCount, 10000), "_splitCount_", splitCount, "_collpased_"));
			OutputStream out(outOptsCurrent);
			compGraph.writeRectangleDotColorBySampleCount(out);
			auto outSeqs = compGraph.nodesToSeqs();
			auto seqsOutOpts = SeqIOOptions::genFastaOut(outOptsCurrent.outFilename_);
			seqsOutOpts.out_.transferOverwriteOpts(outOpts);
			SeqOutput::write(outSeqs, seqsOutOpts);
		}
	}



	if(writeBasesShared){
		std::map<std::string, std::map<std::string, uint32_t>> basesShared;

		for(const auto & n : compGraph.nodes_){
			if(n->on_){
				VecStr sampNames;
				for(const auto & readName : n->inReadNamesIdx_){
					MetaDataInName meta(readName);
					sampNames.emplace_back(meta.getMeta("sample"));
				}
				for(const auto pos1 : iter::range(sampNames.size())){
					for(const auto pos2 : iter::range(sampNames.size())){
						if(pos1 != pos2 && sampNames[pos1] == sampNames[pos2]){
							continue;
						}
						basesShared[sampNames[pos1]][ sampNames[pos2]] += n->k_.size();
						if(n->tailCount() > 0){
							//need to take into account the amount of overlap
							basesShared[sampNames[pos1]][ sampNames[pos2]] -=(n->kLen_ - 1);
						}
					}
				}
			}
		}
		OutputStream sharedBase(njh::files::make_path(setUp.pars_.directoryName_, "basesShared.tab.txt"));
		sharedBase << "sample1\tsample2\tbaseShared\tfrac1\tfrac2" << std::endl;
		for(const auto & samp1 : basesShared){
			for(const auto & samp2 : samp1.second){
				sharedBase << samp1.first
						<< "\t" << samp2.first
						<< "\t" << samp2.second
						<< "\t" << samp2.second/static_cast<double>(basesShared[samp1.first][samp1.first])
						<< "\t" << samp2.second/static_cast<double>(basesShared[samp2.first][samp2.first])<< std::endl;
			}
		}
	}
	{
		OutputStream contigsWtihSamplesOut (njh::files::make_path(setUp.pars_.directoryName_, "contigsWtihSamples.tab.txt"));
		contigsWtihSamplesOut << "contig\tsampleName" << std::endl;
		uint32_t nodePos = 0;
		for(const auto & n : compGraph.nodes_){
			for(const auto & readName : n->inReadNamesIdx_){
				MetaDataInName meta(readName);
				contigsWtihSamplesOut << nodePos
						<< "\t" << meta.getMeta("sample") << std::endl;
			}
			++nodePos;
		}

	}

	return 0;
}



int miscRunner::createSharedPathwaysFromContigs(const njh::progutils::CmdArgs & inputCommands){
	std::vector<std::string> contigFnps;
	std::string contigFnpsStr;
	uint32_t klen = 25;
	uint32_t kmerOccurenceCutOff = 0;
	bool writeBasesShared = false;
	seqSetUp setUp(inputCommands);
	setUp.setOption(contigFnpsStr, "--contigFnps", "contigFnps", true);
	setUp.setOption(klen, "--klen", "k-mer length", true);
	setUp.setOption(writeBasesShared, "--writeBasesShared", "Write Bases Shared");
	setUp.setOption(kmerOccurenceCutOff, "--kmerOccurenceCutOff", "K-mer Occurence Cut Off");

	setUp.processDirectoryOutputName("createSharedPathwaysFromContigs_TODAY", true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	OutOptions outOpts(njh::files::make_path(setUp.pars_.directoryName_, "graph"), ".dot");

	contigFnps = tokenizeString(contigFnpsStr, ",");
	std::vector<std::vector<seqInfo>> inputs;
	for(const auto & contigFnp : contigFnps){
		inputs.emplace_back(SeqInput::getSeqVec<seqInfo>(SeqIOOptions::genFastaIn(contigFnp)));
	}



	ContigsCompareGraphDev compGraph(klen);
	compGraph.setOccurenceCutOff(kmerOccurenceCutOff);
	for(const auto & contigs : inputs){
		for(const auto & contig : contigs){
			compGraph.increaseKCounts(contig.seq_);
		}
	}
	compGraph.populateNodesFromCounts();
	for(const auto & contigs : inputs){
		for(const auto & contig : contigs){
			compGraph.threadThroughSequence(contig, contig.name_);
		}
	}
	compGraph.collapseSingleLinkedPathsSameReads();
	uint32_t graphCount = 0;
	{
		auto outOptsCurrent = outOpts;
		outOptsCurrent.outFilename_ = njh::files::prependFileBasename(outOptsCurrent.outFilename_,
				njh::leftPadNumStr<uint32_t>(graphCount++, 10000) + "_");
		OutputStream out(outOptsCurrent);
		compGraph.writeRectangleDotColorBySampleCount(out);
		auto outSeqs = compGraph.nodesToSeqs();
		auto seqsOutOpts = SeqIOOptions::genFastaOut(outOptsCurrent.outFilename_);
		seqsOutOpts.out_.transferOverwriteOpts(outOpts);
		SeqOutput::write(outSeqs, seqsOutOpts);
	}

	uint32_t splitCount = 0;
	while(compGraph.splitNodesWithRedundantKmers()){
		{
			++splitCount;
			auto outOptsCurrent = outOpts;
			outOptsCurrent.outFilename_ = njh::files::prependFileBasename(outOptsCurrent.outFilename_,
					njh::pasteAsStr(njh::leftPadNumStr<uint32_t>(graphCount++, 10000), "_splitCount_", splitCount, "_"));
			OutputStream out(outOptsCurrent);
			compGraph.writeRectangleDotColorBySampleCount(out);
			auto outSeqs = compGraph.nodesToSeqs();
			auto seqsOutOpts = SeqIOOptions::genFastaOut(outOptsCurrent.outFilename_);
			seqsOutOpts.out_.transferOverwriteOpts(outOpts);
			SeqOutput::write(outSeqs, seqsOutOpts);
		}
		compGraph.collapseSingleLinkedPathsSameReads();

		{
			auto outOptsCurrent = outOpts;
			outOptsCurrent.outFilename_ = njh::files::prependFileBasename(outOptsCurrent.outFilename_,
					njh::pasteAsStr(njh::leftPadNumStr<uint32_t>(graphCount, 10000), "_splitCount_", splitCount, "_collpased_"));
			OutputStream out(outOptsCurrent);
			compGraph.writeRectangleDotColorBySampleCount(out);
			auto outSeqs = compGraph.nodesToSeqs();
			auto seqsOutOpts = SeqIOOptions::genFastaOut(outOptsCurrent.outFilename_);
			seqsOutOpts.out_.transferOverwriteOpts(outOpts);
			SeqOutput::write(outSeqs, seqsOutOpts);
		}
	}



	if(writeBasesShared){
		std::map<std::string, std::map<std::string, uint32_t>> basesShared;

		for(const auto & n : compGraph.nodes_){
			if(n->on_){
				VecStr sampNames;
				for(const auto & readName : n->inReadNamesIdx_){
					MetaDataInName meta(readName);
					sampNames.emplace_back(meta.getMeta("sample"));
				}
				for(const auto pos1 : iter::range(sampNames.size())){
					for(const auto pos2 : iter::range(sampNames.size())){
						if(pos1 != pos2 && sampNames[pos1] == sampNames[pos2]){
							continue;
						}
						basesShared[sampNames[pos1]][ sampNames[pos2]] += n->k_.size();
						if(n->tailCount() > 0){
							//need to take into account the amount of overlap
							basesShared[sampNames[pos1]][ sampNames[pos2]] -=(n->kLen_ - 1);
						}
					}
				}
			}
		}
		OutputStream sharedBase(njh::files::make_path(setUp.pars_.directoryName_, "basesShared.tab.txt"));
		sharedBase << "sample1\tsample2\tbaseShared\tfrac1\tfrac2" << std::endl;
		for(const auto & samp1 : basesShared){
			for(const auto & samp2 : samp1.second){
				sharedBase << samp1.first
						<< "\t" << samp2.first
						<< "\t" << samp2.second
						<< "\t" << samp2.second/static_cast<double>(basesShared[samp1.first][samp1.first])
						<< "\t" << samp2.second/static_cast<double>(basesShared[samp2.first][samp2.first])<< std::endl;
			}
		}
	}
	{
		OutputStream contigsWtihSamplesOut (njh::files::make_path(setUp.pars_.directoryName_, "contigsWtihSamples.tab.txt"));
		contigsWtihSamplesOut << "contig\tsampleName" << std::endl;
		uint32_t nodePos = 0;
		for(const auto & n : compGraph.nodes_){
			for(const auto & readName : n->inReadNamesIdx_){
				MetaDataInName meta(readName);
				contigsWtihSamplesOut << nodePos
						<< "\t" << meta.getMeta("sample") << std::endl;
			}
			++nodePos;
		}
	}

	return 0;
}



int miscRunner::createSharedPathwaysFromRefSeqs(const njh::progutils::CmdArgs & inputCommands){
	uint32_t klen = 25;
	uint32_t kmerOccurenceCutOff = 0;
	uint32_t correctionOccurenceCutOff = 2;
	comparison allowableError;
	allowableError.hqMismatches_ = 1;
	allowableError.oneBaseIndel_ = 1;
	uint32_t lowFreqCutOff = 3;
	bool collapseLowFreqNodes = false;
	bool writeBasesShared = false;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(klen, "--klen", "k-mer length", true);
	setUp.setOption(writeBasesShared, "--writeBasesShared", "Write Bases Shared");
	setUp.setOption(kmerOccurenceCutOff, "--kmerOccurenceCutOff", "K-mer Occurence Cut Off");

	setUp.setOption(collapseLowFreqNodes, "--collapseLowFreqNodes", "collapse Low Freq Nodes");
	setUp.setOption(lowFreqCutOff, "--lowFreqCutOff", "Low Freq Cut Off");
	setUp.setOption(correctionOccurenceCutOff, "--correctionOccurenceCutOff", "correction Occurrence Cut Off");

	setUp.processComparison(allowableError);


	setUp.processReadInNames(true);
	setUp.processDirectoryOutputName(true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	OutOptions outOpts(njh::files::make_path(setUp.pars_.directoryName_, "graph"), ".dot");


	std::vector<seqInfo> seqs;
	{
		seqInfo seq;
		SeqInput reader(setUp.pars_.ioOptions_);
		reader.openIn();
		std::unordered_set<std::string> readNames;
		while(reader.readNextRead(seq)){
//			if(MetaDataInName::nameHasMetaData(seq.name_)){
//				MetaDataInName meta(seq.name_);
//				if(!meta.containsMeta("sample")){
//					auto sampleMeta = seq.name_;
//					MetaDataInName::removeMetaDataInName(sampleMeta);
//					meta.addMeta("sample", sampleMeta);
//					meta.resetMetaInName(seq.name_);
//				}
//			}else{
//				MetaDataInName meta;
//				auto sampleMeta = seq.name_;
//				meta.addMeta("sample", sampleMeta);
//				meta.resetMetaInName(seq.name_);
//			}
			if(njh::in(seq.name_, readNames)){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "already have seq name: " << seq.name_ << ", can't have repeat names"<< "\n";
				throw std::runtime_error{ss.str()};
			}
			seqs.emplace_back(seq);
			readNames.emplace(seq.name_);
		}
	}

	if(setUp.pars_.debug_){
		std::cout << "Read in" << std::endl;
		for(const auto & seq : seqs){
			std::cout << seq.name_ << std::endl;
		}
	}



	{
		njh::stopWatch watch;
		//correction
		ContigsCompareGraphDev compGraph(klen);
		compGraph.setOccurenceCutOff(kmerOccurenceCutOff);
		for(const auto & seq : seqs){
			compGraph.increaseKCounts(seq.seq_);
		}
		compGraph.populateNodesFromCounts();
		for(const auto & seq : seqs){
			compGraph.threadThroughSequence(seq, seq.name_);
		}

		struct SubSeqment{
			SubSeqment(uint32_t pos, uint32_t size):pos_(pos), size_(size){

			}
			std::string head_;
			std::string tail_;

			uint32_t pos_;
			uint32_t size_;

		};
		std::unordered_map<std::string, uint32_t> seqnameToPosition;
		{
			uint32_t seqCount = 0;
			for(const auto & seq : seqs){
				seqnameToPosition[seq.name_] = seqCount++;
			}
		}
		std::vector<seqInfo> correctedSeqs = seqs;
		for(const auto & seq : seqs){
			bool printInfo = seq.name_ == "PE0362-C-0[sample=PE0362-C-0]";
			std::vector<SubSeqment> subPositions;
			for(uint32_t pos = 0; pos < len(seq) + 1 - klen; ++pos){
				std::string subK = seq.seq_.substr(pos, klen);
				if(compGraph.kCounts_[subK] <= correctionOccurenceCutOff){
					if(subPositions.empty() || (subPositions.back().pos_ + subPositions.back().size_ + 1 - klen != pos)){
						subPositions.emplace_back(SubSeqment(pos, klen));
					}else{
						subPositions.back().size_ += 1;
					}
				}
			}
			if(printInfo){
				std::cout << "subPositions.size(): " << subPositions.size() << std::endl;
			}
			if(!subPositions.empty()){
				if(printInfo){
					std::cout << seq.name_ << std::endl;
				}
				njh::sort(subPositions,[](const SubSeqment & s1, const SubSeqment & s2){
					return s1.pos_ > s2.pos_;
				});
				for(const auto & subPosition : subPositions){
					std::string subSegment = seq.seq_.substr(subPosition.pos_, subPosition.size_);
					if(printInfo){
						std::cout << "\t" <<"subsegment: " << subSegment << std::endl;
						std::cout << "\t" << "pos : " << subPosition.pos_ << std::endl;
						std::cout << "\t" << "size: " << subPosition.size_ << std::endl;

					}
					std::string head = "";
					std::string tail = "";
					if(subPosition.pos_ != 0){
						head = seq.seq_.substr(subPosition.pos_ - 1, klen);
						if(printInfo){
							std::cout << "\t\t" << "head: " << std::endl;
							std::cout << "\t\t" << head << std::endl;
							std::cout << "\t\t" << compGraph.kCounts_[head] << std::endl;
						}
					}
					if(subPosition.pos_ + subPosition.size_ != len(seq)){
						tail = seq.seq_.substr(subPosition.pos_ + subPosition.size_ + 1 - klen, klen);
						if(printInfo){
							std::cout << "\t\t" << "tail: " << std::endl;
							std::cout << "\t\t" << tail << std::endl;
							std::cout << "\t\t" << compGraph.kCounts_[tail] << std::endl;
						}
					}
					std::unordered_map<std::string, uint32_t> headPositions;
					std::unordered_map<std::string, uint32_t> tailPositions;
					uint32_t replaceStart = subPosition.pos_;
					uint32_t replaceLen = subPosition.size_;
					if("" != head){
						--replaceStart;
						++replaceLen;
						for(const auto & tailEdge : compGraph.nodes_[compGraph.nodePositions_[head]]->tailEdges_){
							for(const auto & con : tailEdge->connectorInfos_){
								headPositions[con.readName_] = con.headPos_;
							}
						}
					}
					if("" != tail){
						++replaceLen;
						for(const auto & headEdge :compGraph.nodes_[compGraph.nodePositions_[tail]]->headEdges_){
							for(const auto & con : headEdge->connectorInfos_){
								tailPositions[con.readName_] = con.tailPos_;
							}
						}
					}
					if("" == head){
						for(const auto & tailPosition : tailPositions){
							headPositions[tailPosition.first] = 0;
						}
					}
					if("" == tail){
						for(const auto & headPosition : headPositions){
							tailPositions[headPosition.first] = len(seqs[seqnameToPosition[headPosition.first]]);
						}
					}
					std::unordered_map<std::string, uint32_t> subSeqCounts;
					for(const auto & headPosition : headPositions){
						if(headPosition.first == seq.name_){
							continue;
						}
						if(njh::in(headPosition.first, tailPositions) && tailPositions[headPosition.first] > headPosition.second){
							++subSeqCounts[seqs[seqnameToPosition[headPosition.first]].seq_.substr(headPosition.second, tailPositions[headPosition.first] - headPosition.second + klen)];
						}
					}
					std::vector<std::string> subSeqsAbove;
					if(printInfo){
						std::cout << "\t" << "subSeqsCounts" << std::endl;
					}
					for(const auto & subSeq : subSeqCounts){
						if(printInfo){
							std::cout << "\t" << subSeq.first << ": " << subSeq.second << std::endl;
						}
						if(subSeq.second > std::max(correctionOccurenceCutOff, lowFreqCutOff)){
							subSeqsAbove.emplace_back(subSeq.first);
						}
					}
					if(1 == subSeqsAbove.size()){
						if(printInfo){
							std::cout << "replacing: " << replaceStart << "," << replaceLen << std::endl;
							std::cout << "old        : " << correctedSeqs[seqnameToPosition[seq.name_]].seq_.substr(replaceStart, replaceLen) << std::endl;
							std::cout << "replacement: " << subSeqsAbove.front() << std::endl;
						}
						correctedSeqs[seqnameToPosition[seq.name_]].seq_.replace(replaceStart, replaceLen, subSeqsAbove.front());
					}
				}
			}
		}
		if(setUp.pars_.verbose_){
			watch.logLapTimes(std::cout, true, 6, true);
		}
		SeqOutput::write(correctedSeqs, SeqIOOptions::genFastaOutGz(njh::files::make_path(setUp.pars_.directoryName_, "correctedSeqs")));
		seqs = correctedSeqs;
	}


	ContigsCompareGraphDev compGraph(klen);
	compGraph.setOccurenceCutOff(kmerOccurenceCutOff);
	for(const auto & seq : seqs){
		compGraph.increaseKCounts(seq.seq_);
	}
	compGraph.populateNodesFromCounts();
	for(const auto & seq : seqs){
		compGraph.threadThroughSequence(seq, seq.name_);
	}




	compGraph.collapseSingleLinkedPathsSameReads();

	uint32_t graphCount = 0;
	{
		auto outOptsCurrent = outOpts;
		outOptsCurrent.outFilename_ = njh::files::prependFileBasename(outOptsCurrent.outFilename_,
				njh::leftPadNumStr<uint32_t>(graphCount++, 10000) + "_");
		OutputStream out(outOptsCurrent);
		compGraph.writeRectangleDotColorBySampleCount(out);
		auto outSeqs = compGraph.nodesToSeqs();
		auto seqsOutOpts = SeqIOOptions::genFastaOut(outOptsCurrent.outFilename_);
		seqsOutOpts.out_.transferOverwriteOpts(outOpts);
		//SeqOutput::write(outSeqs, seqsOutOpts);
	}
	if(collapseLowFreqNodes){
		uint32_t collapseCount = 0;
		while(compGraph.collapseLowFreqNodes(allowableError, lowFreqCutOff)){
			{
				auto outOptsCurrent = outOpts;
				outOptsCurrent.outFilename_ = njh::files::prependFileBasename(outOptsCurrent.outFilename_,
						njh::leftPadNumStr<uint32_t>(graphCount++, 10000) + "-collapseLowFreq-" + estd::to_string(collapseCount) + "_") ;
				OutputStream out(outOptsCurrent);
				compGraph.writeRectangleDotColorBySampleCount(out);
				auto outSeqs = compGraph.nodesToSeqs();
				auto seqsOutOpts = SeqIOOptions::genFastaOut(outOptsCurrent.outFilename_);
				seqsOutOpts.out_.transferOverwriteOpts(outOpts);
				//SeqOutput::write(outSeqs, seqsOutOpts);
			}
			compGraph.collapseSingleLinkedPathsSameReads();
			{
				auto outOptsCurrent = outOpts;
				outOptsCurrent.outFilename_ = njh::files::prependFileBasename(outOptsCurrent.outFilename_,
						njh::leftPadNumStr<uint32_t>(graphCount++, 10000)+ "-collapseLowFreq-" + estd::to_string(collapseCount) + "-collpased" + "_");
				OutputStream out(outOptsCurrent);
				compGraph.writeRectangleDotColorBySampleCount(out);
				auto outSeqs = compGraph.nodesToSeqs();
				auto seqsOutOpts = SeqIOOptions::genFastaOut(outOptsCurrent.outFilename_);
				seqsOutOpts.out_.transferOverwriteOpts(outOpts);
				//SeqOutput::write(outSeqs, seqsOutOpts);
			}
			++collapseCount;
		}
	}
	uint32_t splitCount = 0;
	while(compGraph.splitNodesWithRedundantKmers()){
		{
			++splitCount;
			auto outOptsCurrent = outOpts;
			outOptsCurrent.outFilename_ = njh::files::prependFileBasename(outOptsCurrent.outFilename_,
					njh::pasteAsStr(njh::leftPadNumStr<uint32_t>(graphCount++, 10000), "_splitCount_", splitCount, "_"));
			OutputStream out(outOptsCurrent);
			compGraph.writeRectangleDotColorBySampleCount(out);
			auto outSeqs = compGraph.nodesToSeqs();
			auto seqsOutOpts = SeqIOOptions::genFastaOut(outOptsCurrent.outFilename_);
			seqsOutOpts.out_.transferOverwriteOpts(outOpts);
			//SeqOutput::write(outSeqs, seqsOutOpts);
		}
		compGraph.collapseSingleLinkedPathsSameReads();

		{
			auto outOptsCurrent = outOpts;
			outOptsCurrent.outFilename_ = njh::files::prependFileBasename(outOptsCurrent.outFilename_,
					njh::pasteAsStr(njh::leftPadNumStr<uint32_t>(graphCount, 10000), "_splitCount_", splitCount, "_collpased_"));
			OutputStream out(outOptsCurrent);
			compGraph.writeRectangleDotColorBySampleCount(out);
			auto outSeqs = compGraph.nodesToSeqs();
			auto seqsOutOpts = SeqIOOptions::genFastaOut(outOptsCurrent.outFilename_);
			seqsOutOpts.out_.transferOverwriteOpts(outOpts);
			//SeqOutput::write(outSeqs, seqsOutOpts);
		}
	}

	{
		{
			auto outOptsCurrent = outOpts;
			outOptsCurrent.outFilename_ = njh::files::prependFileBasename(outOptsCurrent.outFilename_,"final_");
			OutputStream out(outOptsCurrent);
			compGraph.writeRectangleDotColorBySampleCount(out);
		}
		if(seqs.size() < 5){
			auto outOptsCurrent = outOpts;
			outOptsCurrent.outFilename_ = njh::files::prependFileBasename(outOptsCurrent.outFilename_,"final_colorBySampleCombos_");
			OutputStream out(outOptsCurrent);
			auto colors = compGraph.writeRectangleDotColorBySamples(out);
			OutputStream outColorKey(njh::files::make_path(setUp.pars_.directoryName_, "final_colorBySampleCombos_colorKey.tab.txt"));
			outColorKey << "sampleCombo\tcolor" << std::endl;
			for(const auto & color : colors){
				outColorKey << color.first << "\t" << color.second << std::endl;
			}
		}
	}

	{
		std::unordered_map<uint32_t, std::vector<std::shared_ptr<ContigsCompareGraphDev::node>>> nodesWithFull;
		compGraph.resetGroups();

		//get the max count of contigs per group to get the max
		std::unordered_map<uint32_t, std::unordered_set<std::string>> contigCountsPerGroup;
		for(const auto & n : compGraph.nodes_){
			contigCountsPerGroup[n->group_].insert(n->inReadNamesIdx_.begin(), n->inReadNamesIdx_.end());
		}

		for(const auto & n : compGraph.nodes_){
			if(n->inReadNamesIdx_.size() >= contigCountsPerGroup[n->group_].size() ){
				nodesWithFull[n->group_].emplace_back(n);
			}
		}
		OutputStream groupContigCountsOut(njh::files::make_path(setUp.pars_.directoryName_, "groupContigCountsOut.tab.txt"));
		groupContigCountsOut << "group\tcount" << std::endl;
		for(const auto & group : contigCountsPerGroup){
			groupContigCountsOut << group.first << '\t' << group.second.size() << std::endl;
		}
		OutputStream conservedNodeOut(njh::files::make_path(setUp.pars_.directoryName_, "conservedNodesPerGroup.tab.txt"));

		conservedNodeOut << "group\tkmer\tseqCount\ttailCount\theadCount" << std::endl;
		for(const auto & nodes : nodesWithFull){
			for(const auto & n : nodes.second){
				conservedNodeOut << nodes.first
						<< "\t" << n->k_
						<< "\t" << n->inReadNamesIdx_.size()
						<< "\t" << n->tailCount()
						<< "\t" << n->headCount()
						<< std::endl;

			}
		}
		OutputStream conservedNodeDetailedOut(njh::files::make_path(setUp.pars_.directoryName_, "conservedNodeDetailedOut.tab.txt"));
		conservedNodeDetailedOut << "group\tnode\theadOrTail\tedgeCount\theadPos\treadName\ttailPos" << std::endl;
		for(const auto & nodes : nodesWithFull){
			for(const auto & n : nodes.second){
				uint32_t headCount = 0;
				for(const auto & h : n->headEdges_){
					for(const auto & con : h->connectorInfos_){
						conservedNodeDetailedOut << nodes.first
								<< "\t" << n->k_
								<< "\t" << "head"
								<< "\t" << headCount
								<< "\t" << con.headPos_
								<< "\t" << con.readName_
								<< "\t" << con.tailPos_
								<< std::endl;
					}
					++headCount;
				}
				uint32_t tailCount = 0;

				for(const auto & t : n->tailEdges_){
					for(const auto & con : t->connectorInfos_){
						conservedNodeDetailedOut << nodes.first
								<< "\t" << n->k_
								<< "\t" << "tail"
								<< "\t" << tailCount
								<< "\t" << con.headPos_
								<< "\t" << con.readName_
								<< "\t" << con.tailPos_
								<< std::endl;
					}
					++tailCount;
				}
			}
		}
		struct SubSeqPos{
			SubSeqPos(const std::string & subSeq, const uint32_t pos): subSeq_(subSeq), pos_(pos){

			}
			std::string subSeq_;
			uint32_t pos_;
		};
		auto bedDir = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar{"sharedLocs"});
		for(const auto & nodes : nodesWithFull){
			std::unordered_map<std::string, std::vector<SubSeqPos>> seqToNodePositions;

			for(const auto & n : nodes.second){
				std::unordered_map<std::string, uint32_t> nodePositionPerSeq;
				for(const auto & h : n->headEdges_){
					for(const auto & con : h->connectorInfos_){
						nodePositionPerSeq[con.readName_] = con.tailPos_;
					}
				}
				for(const auto & t : n->tailEdges_){
					for(const auto & con : t->connectorInfos_){
						nodePositionPerSeq[con.readName_] = con.headPos_ - (n->k_.size() - n->kLen_);
					}
				}
				for(const auto & seq : nodePositionPerSeq){
					seqToNodePositions[seq.first].emplace_back(n->k_, seq.second);
				}
			}
			std::unordered_set<std::string> paths;

			for(auto & subSeqs : seqToNodePositions){
				njh::sort(subSeqs.second, [](const SubSeqPos & p1, const SubSeqPos & p2 ){
					return p1.pos_ < p2.pos_;
				});
				std::string path = "";

				for(const auto & sub : subSeqs.second){
					if("" != path){
						path += "-";
					}
					path += sub.subSeq_;
				}
				paths.emplace(path);
			}
			if(paths.size() > 1){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << " error more than one path , " << njh::conToStr(paths, ",") << "\n";
				throw std::runtime_error{ss.str()};
			}

			OutputStream sharedLocsOut(njh::files::make_path(bedDir, njh::pasteAsStr(nodes.first, "_sharedLocs.bed")));

			for(auto & subSeqs : seqToNodePositions){
				//std::cout << subSeqs.first <<std::endl;
				for(const auto & sub : subSeqs.second){
					//std::cout << "\t" << sub.pos_ << ": " << sub.subSeq_ << std::endl;
					sharedLocsOut << subSeqs.first
							<< "\t" << sub.pos_
							<< "\t" << sub.pos_ + sub.subSeq_.size()
							<< "\t" << sub.subSeq_
							<< "\t" << sub.subSeq_.size()
							<< "\t" << "+" << std::endl;
				}
			}
		}
	}


	if(writeBasesShared){
		std::map<std::string, std::map<std::string, uint32_t>> basesShared;

		for(const auto & n : compGraph.nodes_){
			if(n->on_){
				VecStr sampNames;
				for(const auto & readName : n->inReadNamesIdx_){
					MetaDataInName meta(readName);
					sampNames.emplace_back(meta.getMeta("sample"));
				}
				for(const auto pos1 : iter::range(sampNames.size())){
					for(const auto pos2 : iter::range(sampNames.size())){
						if(pos1 != pos2 && sampNames[pos1] == sampNames[pos2]){
							continue;
						}
						basesShared[sampNames[pos1]][ sampNames[pos2]] += n->k_.size();
						if(n->tailCount() > 0){
							//need to take into account the amount of overlap
							basesShared[sampNames[pos1]][ sampNames[pos2]] -=(n->kLen_ - 1);
						}
					}
				}
			}
		}
		OutputStream sharedBase(njh::files::make_path(setUp.pars_.directoryName_, "basesShared.tab.txt"));
		sharedBase << "sample1\tsample2\tbaseShared\tfrac1\tfrac2" << std::endl;
		for(const auto & samp1 : basesShared){
			for(const auto & samp2 : samp1.second){
				sharedBase << samp1.first
						<< "\t" << samp2.first
						<< "\t" << samp2.second
						<< "\t" << samp2.second/static_cast<double>(basesShared[samp1.first][samp1.first])
						<< "\t" << samp2.second/static_cast<double>(basesShared[samp2.first][samp2.first])<< std::endl;
			}
		}
	}
//	{
//		OutputStream contigsWtihSamplesOut (njh::files::make_path(setUp.pars_.directoryName_, "contigsWtihSamples.tab.txt"));
//		contigsWtihSamplesOut << "contig\tsampleName" << std::endl;
//		uint32_t nodePos = 0;
//		for(const auto & n : compGraph.nodes_){
//			for(const auto & readName : n->inReadNamesIdx_){
//				MetaDataInName meta(readName);
//				contigsWtihSamplesOut << nodePos << "\t" << meta.getMeta("sample") << std::endl;
//			}
//			++nodePos;
//		}
//	}
	std::unordered_map<std::string, uint32_t> seqLens;
	for(const auto  & seq : seqs){
		seqLens[seq.name_] = len(seq);
	}
	{
		OutputStream outDetailedNodes(njh::files::make_path(setUp.pars_.directoryName_, "nodeDetails.tab.txt"));
		outDetailedNodes << "node\tread\tstart\tend\treadCombo" << std::endl;
		uint32_t nodePos = 0;
		for(const auto & n : compGraph.nodes_){
			if (n->on_) {
				std::unordered_map<std::string, uint32_t> starts;
				std::unordered_map<std::string, uint32_t> ends;
				for(const auto & name : n->inReadNamesIdx_){
					starts[name] = 0;
					ends[name] = seqLens[name];
				}
				for(const auto & edge : n->headEdges_){
					for(const auto & con : edge->connectorInfos_){
						starts[con.readName_] = con.tailPos_ + klen - 1;
					}
				}
				for(const auto & edge : n->tailEdges_){
					for(const auto & con : edge->connectorInfos_){
						ends[con.readName_] = con.headPos_ + klen;
					}
				}
				for(const auto & name : n->inReadNamesIdx_){
					outDetailedNodes << nodePos
							<< "\t" << name
							<< "\t" << starts[name]
							<< "\t" << ends[name]

							<< "\t" << njh::conToStr(n->inReadNamesIdx_, ",") << std::endl;
				}
				++nodePos;
			}
		}
	}

//	if(writeOutForRD3SankeyGraph){
//		auto sankeyOutputDirName = njh::files::makeDir(setUp.pars_.directoryName_,njh::files::MkdirPar("sankeyGraphOut"));
//
//		OutputStream linksOuts(njh::files::make_path(sankeyOutputDirName, "links.tab.txt"));
//		OutputStream nodesOuts(njh::files::make_path(sankeyOutputDirName, "nodes.tab.txt"));
//		OutputStream nodeInfoDetailed(njh::files::make_path(sankeyOutputDirName, "nodes.tab.txt"));
//		nodesOuts << "contig\tsamples" << "\n";
//		linksOuts << "source\ttarget" << "\n";
//		uint32_t nodePos = 0;
//		for(const auto & n : compGraph.nodes_){
//			if(n->on_){
//				++nodePos;
//			}
//		}
//	}
	return 0;
}






int miscRunner::createSharedSubSegmentsFromRefSeqs(const njh::progutils::CmdArgs & inputCommands){
	uint32_t minimumKlen = 12;
	uint32_t kmerOccurenceCutOff = 0;
	uint32_t correctionOccurenceCutOff = 2;
	comparison allowableError;
	allowableError.hqMismatches_ = 1;
	allowableError.oneBaseIndel_ = 1;
	uint32_t lowFreqCutOff = 3;
	bool doNotCollapseLowFreqNodes = false;
	std::string refSeqName = "";
	bfs::path genomeFnp = "";
	bool filterRegionsWithinRegions = false;
	BedUtility::genSubRegionCombosPars subRegionComboPars;
	bool lenFilter = false;
	double lenFiltMultiplier = 0.91;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(lenFilter, "--lenFilter", "filter input sequences on length for artifacts");
	setUp.setOption(lenFiltMultiplier, "--lenFiltMultiplier", "length Filter Multiplier");

	setUp.setOption(subRegionComboPars.includeFromBack, "--includeFromBack", "include from back if the last region isn't shared");
	setUp.setOption(subRegionComboPars.includeFromFront, "--includeFromFront", "include from front if the last region isn't shared");
	setUp.setOption(subRegionComboPars.includeFromSubRegionSize, "--includeFromSubRegionSize", "the maximum to include from the conserved regions");
	setUp.setOption(subRegionComboPars.justToNextRegion, "--justToNextRegion", "Just do the adjacent regions rather than all possible regions");
	setUp.setOption(subRegionComboPars.maxLen, "--maxSubRegionLen", "Maximum subregion size");
	setUp.setOption(subRegionComboPars.minLen, "--minSubRegionLen", "Minimum subregion size");
	setUp.setOption(subRegionComboPars.minBlockRegionLen, "--minBlockRegionLen", "Minimum subregion to start a block from");

	setUp.setOption(filterRegionsWithinRegions, "--filterRegionsWithinRegions", "Filter Regions Within Regions");


	setUp.setOption(genomeFnp, "--genome",
			"Path to a genome to determine the location in", true);
	setUp.setOption(minimumKlen, "--minimumKlen", "minimum k-mer length");
	setUp.setOption(kmerOccurenceCutOff, "--kmerOccurenceCutOff", "K-mer Occurence Cut Off");
	setUp.setOption(doNotCollapseLowFreqNodes, "--doNotCollapseLowFreqNodes", "don't collapse Low Freq Nodes");
	setUp.setOption(lowFreqCutOff, "--lowFreqCutOff", "Low Freq Cut Off");
	setUp.setOption(correctionOccurenceCutOff, "--correctionOccurenceCutOff", "correction Occurrence Cut Off");
	setUp.setOption(refSeqName, "--refSeqName", "The sample name of the sequence to output reference to", true);

	setUp.processComparison(allowableError);


	setUp.processReadInNames(true);
	setUp.processDirectoryOutputName(true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	OutOptions outOpts(njh::files::make_path(setUp.pars_.directoryName_, "graph"), ".dot");

	std::vector<seqInfo> seqs;
	seqInfo refSeq;
	std::vector<uint32_t> seqLens ;
	{
		seqInfo seq;
		SeqInput reader(setUp.pars_.ioOptions_);
		reader.openIn();
		std::unordered_set<std::string> readNames;
		while(reader.readNextRead(seq)){
			seqLens.emplace_back(len(seq));
			bool matchingRefSeqName = false;
			if(seq.name_ == refSeqName){
				matchingRefSeqName = true;
				refSeq = seq;
			}else{
				if (MetaDataInName::nameHasMetaData(seq.name_)) {
					MetaDataInName meta(seq.name_);
					if (meta.containsMeta("sample")) {
						if(meta.getMeta("sample") == refSeqName){
							matchingRefSeqName = true;
						}
					}
				}
			}
			if(matchingRefSeqName){
				if("" != refSeq.name_){
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error " << "already found a seq maching reference name: " << refSeqName << "\n";
					ss << "previous: " << refSeq.name_ << "\n";
					ss << "current: " << seq.name_ << "\n";
					throw std::runtime_error{ss.str()};
				}else{
					refSeq = seq;
				}
			}
			if(njh::in(seq.name_, readNames)){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "already have seq name: " << seq.name_ << ", can't have repeat names"<< "\n";
				throw std::runtime_error{ss.str()};
			}
			seqs.emplace_back(seq);
			readNames.emplace(seq.name_);
		}
	}
	if("" == refSeq.name_){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "didn't find a matching seq for reference: " << refSeqName<< "\n";
		throw std::runtime_error{ss.str()};
	}

	if(lenFilter){
		uint32_t medianLen = std::round(vectorMedianRef(seqLens));
		uint32_t lenDiff = medianLen - std::round(medianLen * lenFiltMultiplier);
		std::vector<seqInfo> filterSeqs;
		std::vector<seqInfo> filteredOffSeqs;
		for(const auto & seq : seqs){
			if(uAbsdiff(medianLen, len(seq)) > lenDiff){
				filteredOffSeqs.emplace_back(seq);
			}else{
				filterSeqs.emplace_back(seq);
			}
		}
//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
//		std::cout << "filteredOffSeqs.size(): " << filteredOffSeqs.size() << std::endl;
		if(!filteredOffSeqs.empty() && filteredOffSeqs.size() <= correctionOccurenceCutOff){
			seqs = filterSeqs;
			auto filtOffOpts = SeqIOOptions::genFastaOut(njh::files::make_path(setUp.pars_.directoryName_, "filteredSeqs.fasta.gz"));
			SeqOutput::write(filteredOffSeqs, filtOffOpts);
		}
	}

	//refSeq.name_ = refSeqName;
	{

		auto refSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(setUp.pars_.directoryName_, "refSeq.fasta"));
		SeqOutput::write(std::vector<seqInfo>{refSeq}, refSeqOpts);
	}
	GenomicRegion refSeqLoc;
	{
		//determine refseq location

		auto genomeName = bfs::basename(genomeFnp);
		if (std::string::npos != genomeName.rfind(".")) {
			genomeName = genomeName.substr(0, genomeName.rfind("."));
		}

		MultiGenomeMapper genomeMapper(genomeFnp.parent_path(),
				genomeName);
		genomeMapper.setSelectedGenomes(VecStr { genomeName });
		genomeMapper.loadInGenomes();
		TwoBit::TwoBitFile tReader(genomeMapper.genomes_[genomeName]->fnpTwoBit_);
		//open out file
		OutputStream out(njh::files::make_path(setUp.pars_.directoryName_, "refSeq.bed"));
		//map reads
		auto outputs = genomeMapper.alignToGenomes(SeqIOOptions::genFastaIn(njh::files::make_path(setUp.pars_.directoryName_, "refSeq.fasta")), setUp.pars_.directoryName_);
		std::unordered_map<std::string, bfs::path> bamFnps;
		for (const auto & output : outputs) {
			bamFnps[output.first] = output.second.alignedFnp_;
		}
		auto bamFnp = bamFnps.begin()->second;
		BamTools::BamReader bReader;
		BamTools::BamAlignment bAln;
		bReader.Open(bamFnp.string());
		checkBamOpenThrow(bReader, bamFnp.string());
		auto refDAta = bReader.GetReferenceData();
		std::vector<AlignmentResults> bestAlnResults;
		double bestAlnScore =std::numeric_limits<double>::min();
		while(bReader.GetNextAlignment(bAln)){
			AlignmentResults alnRes(bAln, refDAta);
			alnRes.setRefSeq(tReader);
			alnRes.setComparison(true);
			if(alnRes.comp_.distances_.eventBasedIdentity_ > bestAlnScore){
				bestAlnResults.clear();
				bestAlnResults.emplace_back(alnRes);
				bestAlnScore = alnRes.comp_.distances_.eventBasedIdentity_;
			}else if(alnRes.comp_.distances_.eventBasedIdentity_ == bestAlnScore){
				bestAlnResults.emplace_back(alnRes);
			}
		}
		if(bestAlnResults.empty()){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "couldn't determine " << refSeq.name_ << " in " << genomeFnp << "\n";
			throw std::runtime_error{ss.str()};
		}else if(bestAlnResults.size() > 1){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << " found multiple possible matches "<< "\n";
			throw std::runtime_error{ss.str()};
		}
		refSeqLoc = bestAlnResults.front().gRegion_;
		out << bestAlnResults.front().gRegion_.genBedRecordCore().toDelimStr() << std::endl;
	}
	uint64_t maxLen = readVec::getMaxLength(seqs);
	uint32_t klen = minimumKlen;

	bool foundLength = false;
	while (klen < maxLen && !foundLength) {
		if (setUp.pars_.verbose_) {
			std::cout << klen << std::endl;
		}
		bool allPass = true;
		for(const auto & seq : seqs) {
			kmerInfo kinfo(seq.seq_, klen, false);
			for (const auto &k : kinfo.kmers_) {
				if (k.second.count_ > 1) {
					allPass = false;
					if (setUp.pars_.verbose_) {
						std::cout << "\t" << seq.name_ << std::endl;
					}
					break;
				}
			}
		}
		if (allPass) {
			foundLength = true;
		} else {
			++klen;
		}
	}
	if(!foundLength){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << " couldn't determine a minimum non-redundant kmer size"<< "\n";
		throw std::runtime_error{ss.str()};
	}

	if(setUp.pars_.debug_){
		std::cout << "Read in" << std::endl;
		for(const auto & seq : seqs){
			std::cout << seq.name_ << std::endl;
		}
	}
	if(std::numeric_limits<uint32_t>::max() == subRegionComboPars.includeFromSubRegionSize){
		subRegionComboPars.includeFromSubRegionSize = klen;
	}


	std::unordered_map<std::string, uint32_t> seqKey;
	for(const auto seq: iter::enumerate(seqs)){
		seqKey[seq.element.name_] = seq.index;
	}

	{
		njh::stopWatch watch;
		//correction
		ContigsCompareGraphDev compGraph(klen);
		compGraph.setOccurenceCutOff(kmerOccurenceCutOff);
		for(const auto & seq : seqs){
			compGraph.increaseKCounts(seq.seq_);
		}
		compGraph.populateNodesFromCounts();
		for(const auto & seq : seqs){
			compGraph.threadThroughSequence(seq, seq.name_);
		}

		struct SubSeqment{
			SubSeqment(uint32_t pos, uint32_t size):pos_(pos), size_(size){

			}
			std::string head_;
			std::string tail_;

			uint32_t pos_;
			uint32_t size_;

		};
		std::unordered_map<std::string, uint32_t> seqnameToPosition;
		{
			uint32_t seqCount = 0;
			for(const auto & seq : seqs){
				seqnameToPosition[seq.name_] = seqCount++;
			}
		}
		std::vector<seqInfo> correctedSeqs = seqs;
		for(const auto & seq : seqs){
			bool printInfo = setUp.pars_.debug_;
			std::vector<SubSeqment> subPositions;
			for(uint32_t pos = 0; pos < len(seq) + 1 - klen; ++pos){
				std::string subK = seq.seq_.substr(pos, klen);
				if(compGraph.kCounts_[subK] <= correctionOccurenceCutOff){
					if(subPositions.empty() || (subPositions.back().pos_ + subPositions.back().size_ + 1 - klen != pos)){
						subPositions.emplace_back(SubSeqment(pos, klen));
					}else{
						subPositions.back().size_ += 1;
					}
				}
			}
			if(printInfo){
				std::cout << "subPositions.size(): " << subPositions.size() << std::endl;
			}
			if(!subPositions.empty()){
				if(printInfo){
					std::cout << seq.name_ << std::endl;
				}
				njh::sort(subPositions,[](const SubSeqment & s1, const SubSeqment & s2){
					return s1.pos_ > s2.pos_;
				});
				for(const auto & subPosition : subPositions){
					std::string subSegment = seq.seq_.substr(subPosition.pos_, subPosition.size_);
					if(printInfo){
						std::cout << "\t" <<"subsegment: " << subSegment << std::endl;
						std::cout << "\t" << "pos : " << subPosition.pos_ << std::endl;
						std::cout << "\t" << "size: " << subPosition.size_ << std::endl;

					}
					std::string head = "";
					std::string tail = "";
					if(subPosition.pos_ != 0){
						head = seq.seq_.substr(subPosition.pos_ - 1, klen);
						if(printInfo){
							std::cout << "\t\t" << "head: " << std::endl;
							std::cout << "\t\t" << head << std::endl;
							std::cout << "\t\t" << compGraph.kCounts_[head] << std::endl;
						}
					}
					if(subPosition.pos_ + subPosition.size_ != len(seq)){
						tail = seq.seq_.substr(subPosition.pos_ + subPosition.size_ + 1 - klen, klen);
						if(printInfo){
							std::cout << "\t\t" << "tail: " << std::endl;
							std::cout << "\t\t" << tail << std::endl;
							std::cout << "\t\t" << compGraph.kCounts_[tail] << std::endl;
						}
					}
					std::unordered_map<std::string, uint32_t> headPositions;
					std::unordered_map<std::string, uint32_t> tailPositions;
					uint32_t replaceStart = subPosition.pos_;
					uint32_t replaceLen = subPosition.size_;
					if("" != head){
						--replaceStart;
						++replaceLen;
						for(const auto & tailEdge : compGraph.nodes_[compGraph.nodePositions_[head]]->tailEdges_){
							for(const auto & con : tailEdge->connectorInfos_){
								headPositions[con.readName_] = con.headPos_;
							}
						}
					}
					if("" != tail){
						++replaceLen;
						for(const auto & headEdge :compGraph.nodes_[compGraph.nodePositions_[tail]]->headEdges_){
							for(const auto & con : headEdge->connectorInfos_){
								tailPositions[con.readName_] = con.tailPos_;
							}
						}
					}
					if("" == head){
						for(const auto & tailPosition : tailPositions){
							headPositions[tailPosition.first] = 0;
						}
					}
					if("" == tail){
						for(const auto & headPosition : headPositions){
							tailPositions[headPosition.first] = len(seqs[seqnameToPosition[headPosition.first]]);
						}
					}
					std::unordered_map<std::string, uint32_t> subSeqCounts;
					for(const auto & headPosition : headPositions){
						if(headPosition.first == seq.name_){
							continue;
						}
						if(njh::in(headPosition.first, tailPositions) && tailPositions[headPosition.first] > headPosition.second){
							++subSeqCounts[seqs[seqnameToPosition[headPosition.first]].seq_.substr(headPosition.second, tailPositions[headPosition.first] - headPosition.second + klen)];
						}
					}
					std::vector<std::string> subSeqsAbove;
					if(printInfo){
						std::cout << "\t" << "subSeqsCounts" << std::endl;
					}
					for(const auto & subSeq : subSeqCounts){
						if(printInfo){
							std::cout << "\t" << subSeq.first << ": " << subSeq.second << std::endl;
						}
						if(subSeq.second > std::max(correctionOccurenceCutOff, lowFreqCutOff)){
							subSeqsAbove.emplace_back(subSeq.first);
						}
					}
					if(1 == subSeqsAbove.size()){
						if(printInfo){
							std::cout << "replacing: " << replaceStart << "," << replaceLen << std::endl;
							std::cout << "old        : " << correctedSeqs[seqnameToPosition[seq.name_]].seq_.substr(replaceStart, replaceLen) << std::endl;
							std::cout << "replacement: " << subSeqsAbove.front() << std::endl;
						}
						correctedSeqs[seqnameToPosition[seq.name_]].seq_.replace(replaceStart, replaceLen, subSeqsAbove.front());
					}
				}
			}
		}
		if(setUp.pars_.verbose_){
			watch.logLapTimes(std::cout, true, 6, true);
		}
		SeqOutput::write(correctedSeqs, SeqIOOptions::genFastaOutGz(njh::files::make_path(setUp.pars_.directoryName_, "correctedSeqs")));
		seqs = correctedSeqs;
	}
	//replace refseq with any error correction done
	refSeq = seqs[seqKey[refSeq.name_]];

	ContigsCompareGraphDev compGraph(klen);
	compGraph.setOccurenceCutOff(kmerOccurenceCutOff);
	for(const auto & seq : seqs){
		compGraph.increaseKCounts(seq.seq_);
	}
	compGraph.populateNodesFromCounts();
	for(const auto & seq : seqs){
		compGraph.threadThroughSequence(seq, seq.name_);
	}




	compGraph.collapseSingleLinkedPathsSameReads();

	uint32_t graphCount = 0;
	{
		auto outOptsCurrent = outOpts;
		outOptsCurrent.outFilename_ = njh::files::prependFileBasename(outOptsCurrent.outFilename_,
				njh::leftPadNumStr<uint32_t>(graphCount++, 10000) + "_");
		OutputStream out(outOptsCurrent);
		compGraph.writeRectangleDotColorBySampleCount(out);
		auto outSeqs = compGraph.nodesToSeqs();
		auto seqsOutOpts = SeqIOOptions::genFastaOut(outOptsCurrent.outFilename_);
		seqsOutOpts.out_.transferOverwriteOpts(outOpts);
		//SeqOutput::write(outSeqs, seqsOutOpts);
	}
	if(!doNotCollapseLowFreqNodes){
		uint32_t collapseCount = 0;
		while(compGraph.collapseLowFreqNodes(allowableError, lowFreqCutOff)){
			{
				auto outOptsCurrent = outOpts;
				outOptsCurrent.outFilename_ = njh::files::prependFileBasename(outOptsCurrent.outFilename_,
						njh::leftPadNumStr<uint32_t>(graphCount++, 10000) + "-collapseLowFreq-" + estd::to_string(collapseCount) + "_") ;
				OutputStream out(outOptsCurrent);
				compGraph.writeRectangleDotColorBySampleCount(out);
				auto outSeqs = compGraph.nodesToSeqs();
				auto seqsOutOpts = SeqIOOptions::genFastaOut(outOptsCurrent.outFilename_);
				seqsOutOpts.out_.transferOverwriteOpts(outOpts);
				//SeqOutput::write(outSeqs, seqsOutOpts);
			}
			compGraph.collapseSingleLinkedPathsSameReads();
			{
				auto outOptsCurrent = outOpts;
				outOptsCurrent.outFilename_ = njh::files::prependFileBasename(outOptsCurrent.outFilename_,
						njh::leftPadNumStr<uint32_t>(graphCount++, 10000)+ "-collapseLowFreq-" + estd::to_string(collapseCount) + "-collpased" + "_");
				OutputStream out(outOptsCurrent);
				compGraph.writeRectangleDotColorBySampleCount(out);
				auto outSeqs = compGraph.nodesToSeqs();
				auto seqsOutOpts = SeqIOOptions::genFastaOut(outOptsCurrent.outFilename_);
				seqsOutOpts.out_.transferOverwriteOpts(outOpts);
				//SeqOutput::write(outSeqs, seqsOutOpts);
			}
			++collapseCount;
		}
	}
	uint32_t splitCount = 0;
	while(compGraph.splitNodesWithRedundantKmers()){
		{
			++splitCount;
			auto outOptsCurrent = outOpts;
			outOptsCurrent.outFilename_ = njh::files::prependFileBasename(outOptsCurrent.outFilename_,
					njh::pasteAsStr(njh::leftPadNumStr<uint32_t>(graphCount++, 10000), "_splitCount_", splitCount, "_"));
			OutputStream out(outOptsCurrent);
			compGraph.writeRectangleDotColorBySampleCount(out);
			auto outSeqs = compGraph.nodesToSeqs();
			auto seqsOutOpts = SeqIOOptions::genFastaOut(outOptsCurrent.outFilename_);
			seqsOutOpts.out_.transferOverwriteOpts(outOpts);
			//SeqOutput::write(outSeqs, seqsOutOpts);
		}
		compGraph.collapseSingleLinkedPathsSameReads();

		{
			auto outOptsCurrent = outOpts;
			outOptsCurrent.outFilename_ = njh::files::prependFileBasename(outOptsCurrent.outFilename_,
					njh::pasteAsStr(njh::leftPadNumStr<uint32_t>(graphCount, 10000), "_splitCount_", splitCount, "_collpased_"));
			OutputStream out(outOptsCurrent);
			compGraph.writeRectangleDotColorBySampleCount(out);
			auto outSeqs = compGraph.nodesToSeqs();
			auto seqsOutOpts = SeqIOOptions::genFastaOut(outOptsCurrent.outFilename_);
			seqsOutOpts.out_.transferOverwriteOpts(outOpts);
			//SeqOutput::write(outSeqs, seqsOutOpts);
		}
	}

	{
		{
			auto outOptsCurrent = outOpts;
			outOptsCurrent.outFilename_ = njh::files::prependFileBasename(outOptsCurrent.outFilename_,"final_");
			OutputStream out(outOptsCurrent);
			compGraph.writeRectangleDotColorBySampleCount(out);
		}
	}

	{
		std::unordered_map<uint32_t, std::vector<std::shared_ptr<ContigsCompareGraphDev::node>>> nodesWithFull;
		compGraph.resetGroups();

		//get the max count of contigs per group to get the max
		std::unordered_map<uint32_t, std::unordered_set<std::string>> contigCountsPerGroup;
		for(const auto & n : compGraph.nodes_){
			contigCountsPerGroup[n->group_].insert(n->inReadNamesIdx_.begin(), n->inReadNamesIdx_.end());
		}

		for(const auto & n : compGraph.nodes_){
			if(n->inReadNamesIdx_.size() >= contigCountsPerGroup[n->group_].size() ){
				nodesWithFull[n->group_].emplace_back(n);
			}
		}
		OutputStream groupContigCountsOut(njh::files::make_path(setUp.pars_.directoryName_, "groupContigCountsOut.tab.txt"));
		groupContigCountsOut << "group\tcount" << std::endl;
		for(const auto & group : contigCountsPerGroup){
			groupContigCountsOut << group.first << '\t' << group.second.size() << std::endl;
		}
		OutputStream conservedNodeOut(njh::files::make_path(setUp.pars_.directoryName_, "conservedNodesPerGroup.tab.txt"));
		conservedNodeOut << "group\tsubSegment\tseqCount\ttailCount\theadCount" << std::endl;

		for(const auto & nodes : nodesWithFull){
			std::vector<seqInfo> conservedSeqs;

			for(const auto & n : nodes.second){
				conservedSeqs.emplace_back(seqInfo(n->k_, n->k_));
				conservedNodeOut << nodes.first
						<< "\t" << n->k_
						<< "\t" << n->inReadNamesIdx_.size()
						<< "\t" << n->tailCount()
						<< "\t" << n->headCount()
						<< std::endl;

			}
			auto conservedSeqsOpts = SeqIOOptions::genFastaOut(njh::files::make_path(setUp.pars_.directoryName_, njh::pasteAsStr(nodes.first, "_conservedNodes.fasta")));
			SeqOutput::write(conservedSeqs, conservedSeqsOpts);
		}

		struct SubSeqPos{
			SubSeqPos(const std::string & subSeq, const uint32_t pos): subSeq_(subSeq), pos_(pos){

			}
			std::string subSeq_;
			uint32_t pos_;

			Bed6RecordCore genBedRegion(const std::string & chrom)const{
				return Bed6RecordCore(chrom, pos_, pos_ + subSeq_.size(),subSeq_,subSeq_.size(), '+' );
			}
		};
		auto bedDir = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar{"sharedLocs"});
		auto seqsDir = njh::files::makeDir(bedDir, njh::files::MkdirPar{"subSeqs"});
		OutputStream outDivMeasures(njh::files::make_path(bedDir, "divMeasuresPerSubRegion.tab.txt"));
		outDivMeasures << "target\ttotalHaps\tuniqueHaps\the\tsinglets\tdoublets\teffectiveNumOfAlleles\tShannonEntropyE" << std::endl;

		//first full div
		{
			std::unordered_map<std::string, uint32_t> popCounts;
			auto divMeasures = PopGenCalculator::getGeneralMeasuresOfDiversityRawInput(seqs);
			outDivMeasures << refSeqLoc.createUidFromCoordsStrand()
													<< "\t" << seqs.size()
													<< "\t" << divMeasures.alleleNumber_
													<< "\t" << divMeasures.heterozygostiy_
													<< "\t" << divMeasures.singlets_
													<< "\t" << divMeasures.doublets_
													<< "\t" << divMeasures.effectiveNumOfAlleles_
													<< "\t" << divMeasures.ShannonEntropyE_
													<< std::endl;
		}
		OutputStream sharedLocsOutRefSeq(njh::files::make_path(bedDir, njh::pasteAsStr("refSeqGenomicLocs.bed")));
		table outSubSeqs(VecStr{"target", "frontSeq", "endSeq"});

		for(const auto & nodes : nodesWithFull){
			std::vector<GenomicRegion> refSeqPositionsRelativeToCurentSeq;
			std::unordered_map<std::string, std::vector<SubSeqPos>> seqToNodePositions;
			if(1 == nodes.second.size()){
				for(const auto & name : nodes.second.front()->inReadNamesIdx_){
					seqToNodePositions[name].emplace_back(nodes.second.front()->k_, 0);
				}
			}else if(nodes.second.size() > 2){
				for(const auto & n : nodes.second){
					std::unordered_map<std::string, uint32_t> nodePositionPerSeq;
					for(const auto & h : n->headEdges_){
						for(const auto & con : h->connectorInfos_){
							nodePositionPerSeq[con.readName_] = con.tailPos_;
						}
					}
					for(const auto & t : n->tailEdges_){
						for(const auto & con : t->connectorInfos_){
							nodePositionPerSeq[con.readName_] = con.headPos_ - (n->k_.size() - n->kLen_);
						}
					}
					for(const auto & seq : nodePositionPerSeq){
						seqToNodePositions[seq.first].emplace_back(n->k_, seq.second);
					}
				}
			}else{
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "nodes for " << nodes.first << " is empty, size:" << nodes.second.size()<< "\n";
				throw std::runtime_error{ss.str()};
			}

			std::unordered_set<std::string> paths;
			for(auto & subSeqs : seqToNodePositions){
				njh::sort(subSeqs.second, [](const SubSeqPos & p1, const SubSeqPos & p2 ){
					return p1.pos_ < p2.pos_;
				});
				std::string path = "";

				for(const auto & sub : subSeqs.second){
					if("" != path){
						path += "-";
					}
					path += sub.subSeq_;
				}
				paths.emplace(path);
			}
			if(paths.size() > 1){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << " error more than one path , " << njh::conToStr(paths, ",") << "\n";
				throw std::runtime_error{ss.str()};
			}

			//OutputStream sharedLocsOut(njh::files::make_path(bedDir, njh::pasteAsStr(nodes.first, "_sharedLocs.bed")));
			OutputStream refSharedLocsOut(njh::files::make_path(bedDir, njh::pasteAsStr(nodes.first, "_ref_sharedLocs.bed")));


			std::unordered_map<std::string, std::unordered_map<std::string, Bed6RecordCore>> subSeqToReadNameToPos;

			for(auto & subSeqs : seqToNodePositions){
				//std::cout << subSeqs.first <<std::endl;
				for(const auto & sub : subSeqs.second){
					auto outRegion = sub.genBedRegion(subSeqs.first);
					if(subSeqs.first == refSeq.name_){
						refSeqPositionsRelativeToCurentSeq.emplace_back(GenomicRegion(outRegion));
					}
					subSeqToReadNameToPos[sub.subSeq_].emplace(subSeqs.first, outRegion);
					//std::cout << "\t" << sub.pos_ << ": " << sub.subSeq_ << std::endl;
					//sharedLocsOut << outRegion.toDelimStr()<< std::endl;
				}
			}
			std::unordered_map<std::string, uint32_t> lengths;
			lengths[refSeq.name_] = len(refSeq);
			std::vector<Bed6RecordCore> refSeqPositionsRelativeToCurentSeqBeds;
			std::vector<Bed6RecordCore> refSeqPositionsBeds;

			for(const auto & g : refSeqPositionsRelativeToCurentSeq){
				auto bedReg = g.genBedRecordCore();
				refSeqPositionsRelativeToCurentSeqBeds.emplace_back(bedReg);
				refSharedLocsOut << bedReg.toDelimStr() << std::endl;
			}


			std::vector<BedUtility::SubRegionCombo> subRegions = BedUtility::genSubRegionCombos(refSeqPositionsRelativeToCurentSeqBeds, lengths, subRegionComboPars);
			OutputStream refSubRegionsLocsOut(njh::files::make_path(bedDir, njh::pasteAsStr(nodes.first, "_ref_subRegions.bed")));


			//key = sub seq UID, val = 1st node seq start, 2nd node seq stop
			std::unordered_map<std::string, std::pair<std::string, std::string>> subPortionsByConservedSeq;

			njh::sort(subRegions,[](const BedUtility::SubRegionCombo & com1, const BedUtility::SubRegionCombo & com2){
				if(com1.startRegion_.chromStart_ == com2.startRegion_.chromStart_){
					return com1.endRegion_.chromEnd_ > com2.endRegion_.chromEnd_;
				}else{
					return com1.startRegion_.chromStart_ < com2.startRegion_.chromStart_;
				}
			});
			if(!refSeqPositionsRelativeToCurentSeqBeds.empty()){
				//getting the region with the last and first conserved regions
				BedUtility::coordSort(refSeqPositionsRelativeToCurentSeqBeds, false);
				Bed6RecordCore genomicLoc = refSeqPositionsRelativeToCurentSeqBeds.front();
				genomicLoc.chromEnd_ = refSeqPositionsRelativeToCurentSeqBeds.back().chromEnd_;
				genomicLoc.score_ = genomicLoc.length();
				Bed6RecordCore relativeBedReg = genomicLoc;
				genomicLoc.chrom_ = refSeqLoc.chrom_;
				if(refSeqLoc.reverseSrand_){
					genomicLoc.strand_ = '-';
					genomicLoc.chromStart_ = refSeqLoc.start_ +(len(refSeq) - 1 - relativeBedReg.chromEnd_ );
					genomicLoc.chromEnd_ = refSeqLoc.start_ + (len(refSeq) - 1 - relativeBedReg.chromStart_ );
				} else {
					genomicLoc.chromStart_ = refSeqLoc.start_ + relativeBedReg.chromStart_;
					genomicLoc.chromEnd_ = refSeqLoc.start_ + relativeBedReg.chromEnd_;
					genomicLoc.strand_ = '+';
				}
				GenomicRegion gRegion(genomicLoc);
				auto uid = gRegion.createUidFromCoordsStrand();
				genomicLoc.name_ = uid;
				OutputStream refSharedLocsOut(njh::files::make_path(bedDir, njh::pasteAsStr(nodes.first, "_refFirstToLastSharedRegion.bed")));
				refSharedLocsOut << genomicLoc.toDelimStr() << std::endl;
			}
			if(refSeqPositionsRelativeToCurentSeqBeds.size() > 1 ){
				//getting the region with the last and first conserved regions
				BedUtility::coordSort(refSeqPositionsRelativeToCurentSeqBeds, false);
				BedUtility::SubRegionCombo subReg(refSeqPositionsRelativeToCurentSeqBeds.front(), refSeqPositionsRelativeToCurentSeqBeds.back());
				auto bedReg = subReg.genOut(subRegionComboPars.includeFromSubRegionSize);;
//				std::cout << bedReg.toDelimStr() << std::endl;
				Bed6RecordCore genomicLoc = bedReg;
				genomicLoc.chrom_ = refSeqLoc.chrom_;
				auto startSubReg = subReg.genSubStart(subRegionComboPars.includeFromSubRegionSize);;
				auto endSubReg   = subReg.genSubEnd(subRegionComboPars.includeFromSubRegionSize);;
				seqInfo startSubSeq = refSeq.getSubRead(startSubReg.chromStart_, startSubReg.length());
				seqInfo endSubSeq = refSeq.getSubRead(endSubReg.chromStart_, endSubReg.length());

				if(refSeqLoc.reverseSrand_){
					genomicLoc.strand_ = '-';
					genomicLoc.chromStart_ = refSeqLoc.start_ +(len(refSeq) - 1 - bedReg.chromEnd_ );
					genomicLoc.chromEnd_ = refSeqLoc.start_ + (len(refSeq) - 1 - bedReg.chromStart_ );
				}else{
					genomicLoc.chromStart_ = refSeqLoc.start_ + bedReg.chromStart_;
					genomicLoc.chromEnd_ = refSeqLoc.start_ + bedReg.chromEnd_;
					genomicLoc.strand_ = '+';
				}

				GenomicRegion gRegion(genomicLoc);
				auto uid = gRegion.createUidFromCoordsStrand();
				seqInfo startSubSeqFull = refSeq.getSubRead(subReg.startRegion_.chromStart_, subReg.startRegion_.length());
				seqInfo endSubSeqFull = refSeq.getSubRead(subReg.endRegion_.chromStart_, subReg.endRegion_.length());
				genomicLoc.name_ = uid;
				OutputStream refSharedLocsOut(njh::files::make_path(bedDir, njh::pasteAsStr(nodes.first, "_refFirstToLastSharedRegionSlimmer.bed")));
				refSharedLocsOut << genomicLoc.toDelimStr() << std::endl;
			}
			if(filterRegionsWithinRegions){
				std::vector<BedUtility::SubRegionCombo> filteredRegions;

				for(const auto & subRegion : subRegions){
					bool foundWithInAnother = false;
					for(const auto & other : filteredRegions){
						if(subRegion.startRegion_.chromStart_ >= other.startRegion_.chromStart_
						&& subRegion.endRegion_.chromEnd_     <= other.endRegion_.chromEnd_){
							foundWithInAnother = true;
							break;
						}
					}
					if(!foundWithInAnother){
						filteredRegions.emplace_back(subRegion);
					}
				}
				subRegions = filteredRegions;
			}

			for(const auto & subReg : subRegions){
				auto bedReg = subReg.genOut(subRegionComboPars.includeFromSubRegionSize);;
				Bed6RecordCore genomicLoc = bedReg;
				genomicLoc.chrom_ = refSeqLoc.chrom_;
				auto startSubReg = subReg.genSubStart(subRegionComboPars.includeFromSubRegionSize);;
				auto endSubReg   = subReg.genSubEnd(subRegionComboPars.includeFromSubRegionSize);;
				seqInfo startSubSeq = refSeq.getSubRead(startSubReg.chromStart_, startSubReg.length());
				seqInfo endSubSeq = refSeq.getSubRead(endSubReg.chromStart_, endSubReg.length());

				if(refSeqLoc.reverseSrand_){
					genomicLoc.strand_ = '-';
					genomicLoc.chromStart_ = refSeqLoc.start_ +(len(refSeq) - 1 - bedReg.chromEnd_ );
					genomicLoc.chromEnd_ = refSeqLoc.start_ + (len(refSeq) - 1 - bedReg.chromStart_ );
				}else{
					genomicLoc.chromStart_ = refSeqLoc.start_ + bedReg.chromStart_;
					genomicLoc.chromEnd_ = refSeqLoc.start_ + bedReg.chromEnd_;
					genomicLoc.strand_ = '+';
				}

				GenomicRegion gRegion(genomicLoc);
				auto uid = gRegion.createUidFromCoordsStrand();
				seqInfo startSubSeqFull = refSeq.getSubRead(subReg.startRegion_.chromStart_, subReg.startRegion_.length());
				seqInfo endSubSeqFull = refSeq.getSubRead(subReg.endRegion_.chromStart_, subReg.endRegion_.length());

				genomicLoc.name_ = uid;

				subPortionsByConservedSeq[uid] = std::make_pair(startSubSeqFull.seq_, endSubSeqFull.seq_);
				outSubSeqs.addRow(uid, startSubSeq.seq_, endSubSeq.seq_);
				refSeqPositionsBeds.emplace_back(genomicLoc);
				refSubRegionsLocsOut << bedReg.toDelimStr() << std::endl;
			}

			for(const auto & subReg : refSeqPositionsBeds){
				sharedLocsOutRefSeq << subReg.toDelimStr()<< std::endl;
			}

			//cut input to subsegments
			std::unordered_map<std::string, std::vector<seqInfo>> subRegionsSeqs;
//			std::cout << "subPortionsByConservedSeq.size(): " << subPortionsByConservedSeq.size() << std::endl;
//			std::cout << "subSeqToReadNameToPos.size(): " << subSeqToReadNameToPos.size() << std::endl;
//			std::cout << njh::conToStr(getVectorOfMapKeys(subSeqToReadNameToPos), "\n") << std::endl;

			for(const auto & subRegion : subPortionsByConservedSeq){
//				std::cout << subRegion.second.first << std::endl;

				for(const auto & frontRegion : subSeqToReadNameToPos[subRegion.second.first]){
//					std::cout << "\t" << subRegion.second.first << std::endl;
//					std::cout << "\t" << subRegion.second.second << std::endl;

					auto startRegion = frontRegion.second;
					auto endRegion = njh::mapAt(njh::mapAt(subSeqToReadNameToPos,subRegion.second.second), frontRegion.first);
					BedUtility::SubRegionCombo subRegionForSeq(startRegion, endRegion);
					auto outRegion = subRegionForSeq.genOut(subRegionComboPars.includeFromSubRegionSize);
					seqInfo subSeq(seqs[seqKey[frontRegion.first]].getSubRead(outRegion.chromStart_, outRegion.length()));
					MetaDataInName seqMeta(subSeq.name_);
					seqMeta.addMeta("SubRegUID", subRegion.first, true);
					seqMeta.addMeta("SeqStart", outRegion.chromStart_, true);
					seqMeta.addMeta("SeqStop", outRegion.chromEnd_, true);
					seqMeta.resetMetaInName(subSeq.name_);
					subRegionsSeqs[subRegion.first].emplace_back(subSeq);
				}
			}
//			std::cout << "subRegionsSeqs.size(): " << subRegionsSeqs.size() << std::endl;

			for(const auto & subSeqs : subRegionsSeqs){
				auto subSeqOpts = SeqIOOptions::genFastaOutGz(njh::files::make_path(seqsDir, subSeqs.first + ".fasta.gz"));
				SeqOutput::write(subSeqs.second, subSeqOpts);
				auto divMeasures = PopGenCalculator::getGeneralMeasuresOfDiversityRawInput(subSeqs.second);
				outDivMeasures << subSeqs.first
														<< "\t" << subSeqs.second.size()
														<< "\t" << divMeasures.alleleNumber_
														<< "\t" << divMeasures.heterozygostiy_
														<< "\t" << divMeasures.singlets_
														<< "\t" << divMeasures.doublets_
														<< "\t" << divMeasures.effectiveNumOfAlleles_
														<< "\t" << divMeasures.ShannonEntropyE_
														<< '\n';
			}
		}
		OutputStream outSubSeqsOut(njh::files::make_path(bedDir, "startEndSeqs.tab.txt"));
		outSubSeqs.outPutContents(outSubSeqsOut, "\t");


	}
	return 0;
}




} // namespace njhseq

