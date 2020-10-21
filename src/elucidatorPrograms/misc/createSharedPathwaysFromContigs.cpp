/*
 * createSharedPathwaysFromContigs.cpp
 *
 *  Created on: Nov 5, 2019
 *      Author: nicholashathaway
 */




#include "miscRunner.hpp"
#include "elucidator/simulation.h"
#include "elucidator/BamToolsUtils.h"
#include "elucidator/objects/dataContainers/graphs/ContigsCompareGraph.hpp"

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
	uint32_t lowFreqCutOff = 0;
//	bool writeOutForRD3SankeyGraph = false;
	bool collapseLowFreqNodes = false;
	bool writeBasesShared = false;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(klen, "--klen", "k-mer length", true);
	setUp.setOption(writeBasesShared, "--writeBasesShared", "Write Bases Shared");
	setUp.setOption(kmerOccurenceCutOff, "--kmerOccurenceCutOff", "K-mer Occurence Cut Off");
//	setUp.setOption(writeOutForRD3SankeyGraph, "--writeOutForRD3SankeyGraph", "writeOutForRD3SankeyGraph");

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
		SeqOutput::write(correctedSeqs, SeqIOOptions::genFastaOut(njh::files::make_path(setUp.pars_.directoryName_, "correctedSeqs")));
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
		SeqOutput::write(outSeqs, seqsOutOpts);
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
				SeqOutput::write(outSeqs, seqsOutOpts);
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
				SeqOutput::write(outSeqs, seqsOutOpts);
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

} // namespace njhseq

