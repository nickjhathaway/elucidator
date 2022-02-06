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


#include <njhseq/PopulationGenetics.h>
#include <njhseq/objects/seqContainers/CollapsedHaps.hpp>
#include <njhseq/objects/Gene/TranslatorByAlignment.hpp>
#include <boost/filesystem.hpp>


namespace njhseq {
namespace bfs = boost::filesystem;




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
			if(len(seq) <= klen){
				continue;
			}
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




GenomicRegion determineRefSeqLocation(const seqInfo & refSeq,
		const bfs::path & refSeqFnp,
		const bfs::path & genomeFnp){

	SeqOutput::write(std::vector<seqInfo>{refSeq}, SeqIOOptions::genFastaOut(refSeqFnp));
	GenomicRegion refSeqLoc;
	{
		//determine refseq location

		auto genomeName = bfs::basename(genomeFnp);
		if (std::string::npos != genomeName.rfind(".")) {
			genomeName = genomeName.substr(0, genomeName.rfind("."));
		}

		MultiGenomeMapper genomeMapper(genomeFnp.parent_path(), genomeName);
		genomeMapper.setSelectedGenomes(VecStr { genomeName });
		genomeMapper.loadInGenomes();
		TwoBit::TwoBitFile tReader(genomeMapper.genomes_[genomeName]->fnpTwoBit_);
		//map reads
		auto outputs = genomeMapper.alignToGenomes(SeqIOOptions::genFastaIn(refSeqFnp), refSeqFnp.parent_path().string() + "/");
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
		//out << bestAlnResults.front().gRegion_.genBedRecordCore().toDelimStr() << std::endl;
	}
	return refSeqLoc;
}


int miscRunner::createSharedSubSegmentsFromRefSeqs(const njh::progutils::CmdArgs & inputCommands){
	uint32_t minimumKlen = 12;


//	uint32_t lowFreqCutOff = 3;
//	bool doNotCollapseLowFreqNodes = false;
	//	uint32_t kmerOccurenceCutOff = 0;
	//	uint32_t correctionOccurenceCutOff = 2;
	ContigsCompareGraphDev::correctSeqsByGraphPars graphCorrectingPars;

	std::string refSeqName = "";
	bfs::path genomeFnp = "";
	seqInfo refSeq;
	//bool filterRegionsWithinRegions = false;
	BedUtility::genSubRegionCombosPars subRegionComboPars;

	bool lenFilter = false;
	bool kSimFilter = false;
	uint32_t kLenForFilter= 5;
	double kSimCutOff =  0.40;
	double lenFiltMultiplier = 0.91;

	uint32_t numThreads = 1;

	TranslatorByAlignment::RunPars variantCallerRunPars;
	CollapsedHaps::GenPopMeasuresPar calcPopMeasuresPars;
	calcPopMeasuresPars.getPairwiseComps = true;

	variantCallerRunPars.lowVariantCutOff = 0.005;
	variantCallerRunPars.occurrenceCutOff = 1;
	calcPopMeasuresPars.lowVarFreq = variantCallerRunPars.lowVariantCutOff;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();

	bool refNameRequired = setUp.processSeq(refSeq, "--refSeq", "Reference seq", false);
	setUp.setOption(numThreads, "--numThreads", "num Threads");

	setUp.setOption(kSimFilter, "--kSimFilter", "filter input sequences on kSimFilter for artifacts");
	setUp.setOption(kLenForFilter, "--kLenForFilter", "filter kLenForFilter for artifacts");
	setUp.setOption(kSimCutOff, "--kSimCutOff", "filter kSimCutOff for artifacts");


	setUp.setOption(lenFilter, "--lenFilter", "filter input sequences on length for artifacts");
	setUp.setOption(lenFiltMultiplier, "--lenFiltMultiplier", "length Filter Multiplier");

	setUp.setOption(subRegionComboPars.includeFromBack, "--includeFromBack", "include from back if the last region isn't shared");
	setUp.setOption(subRegionComboPars.includeFromFront, "--includeFromFront", "include from front if the last region isn't shared");
	setUp.setOption(subRegionComboPars.includeFromSubRegionSize, "--includeFromSubRegionSize", "the maximum to include from the conserved regions");
	setUp.setOption(subRegionComboPars.justToNextRegion, "--justToNextRegion", "Just do the adjacent regions rather than all possible regions");
	setUp.setOption(subRegionComboPars.maxLen, "--maxSubRegionLen", "Maximum subregion size");
	setUp.setOption(subRegionComboPars.minLen, "--minSubRegionLen", "Minimum subregion size");
	setUp.setOption(subRegionComboPars.minBlockRegionLen, "--minBlockRegionLen", "Minimum subregion to start a block from");
	setUp.setOption(subRegionComboPars.filterFinalRegionsWithinRegions, "--filterRegionsWithinRegions", "Filter Regions Within Regions");
	subRegionComboPars.filterFinalRegionsToNonoverlapingOrAdjacentRegions = true;
	//setUp.setOption(subRegionComboPars.filterFinalRegionsToNonoverlapingOrAdjacentRegions, "--filterRegionsWithinRegions", "Filter Regions Within Regions");


	setUp.setOption(genomeFnp, "--genome",
			"Path to a genome to determine the location in", true);
	setUp.setOption(minimumKlen, "--minimumKlen", "minimum k-mer length");
	setUp.setOption(graphCorrectingPars.kmerOccurenceCutOff, "--kmerOccurenceCutOff", "K-mer Occurence Cut Off");
	setUp.setOption(graphCorrectingPars.doNotCollapseLowFreqNodes, "--doNotCollapseLowFreqNodes", "don't collapse Low Freq Nodes");
	setUp.setOption(graphCorrectingPars.lowFreqCutOff, "--lowFreqCutOff", "Low Freq Cut Off");
	setUp.setOption(graphCorrectingPars.correctionOccurenceCutOff, "--correctionOccurenceCutOff", "correction Occurrence Cut Off");
	if(setUp.setOption(refSeqName, "--refSeqName", "The sample name of the sequence to output reference to", !refNameRequired)){
		refSeq = seqInfo();//create empty refSeq if setting a ref seqname
	}

	setUp.processComparison(graphCorrectingPars.allowableError);


	setUp.processReadInNames(true);
	setUp.processDirectoryOutputName(true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);

	auto conservedRegionInfoDir = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar{"subRegionInfo"});
	auto subsegmentInfoDir = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar{"subsegmentInfo"});

	if ("" != refSeq.name_) {
		MetaDataInName refSeqMeta;
		if (MetaDataInName::nameHasMetaData(refSeq.name_)) {
			refSeqMeta = MetaDataInName(refSeq.name_);
		}
		refSeqMeta.addMeta("sample", refSeq.name_, true);
		refSeqMeta.resetMetaInName(refSeq.name_);
	}
	std::vector<seqInfo> seqs = ContigsCompareGraphDev::readInSeqs(setUp.pars_.ioOptions_, refSeq, refSeqName);
	if (lenFilter) {
		seqs = ContigsCompareGraphDev::filterSeqsOnLen(seqs, graphCorrectingPars,
				lenFiltMultiplier,
				SeqIOOptions::genFastaOut(
						njh::files::make_path(setUp.pars_.directoryName_,
								"filteredSeqs.fasta.gz")));
	}

	if (kSimFilter) {
		kmerInfo refSeqKInfo(refSeq.seq_, kLenForFilter, false);
		seqs = ContigsCompareGraphDev::filterSeqsOnKmerSim(seqs,
				graphCorrectingPars, refSeqKInfo, kSimCutOff,
				SeqIOOptions::genFastaOut(
						njh::files::make_path(setUp.pars_.directoryName_,
								"kSimCutOffFilteredSeqs.fasta.gz")));
	}

	auto refSeqOutFnp = njh::files::make_path(setUp.pars_.directoryName_, "refSeq.fasta");
	GenomicRegion refSeqLoc = determineRefSeqLocation(refSeq, refSeqOutFnp, genomeFnp);
	{
		//open out file
		OutputStream out(njh::files::make_path(setUp.pars_.directoryName_, "refSeq.bed"));
		out << refSeqLoc.genBedRecordCore().toDelimStrWithExtra() << std::endl;
	}

	//uint32_t klen = minimumKlen;
	graphCorrectingPars.klen = ContigsCompareGraphDev::findMinNonredundantKmer(minimumKlen, seqs, __PRETTY_FUNCTION__);

	if(setUp.pars_.debug_){
		std::cout << "Read in" << std::endl;
		for(const auto & seq : seqs){
			std::cout << seq.name_ << std::endl;
		}
	}
	{
		OutputStream klenOut(njh::files::make_path(setUp.pars_.directoryName_, "nonRedundantKlen.txt"));
		klenOut << graphCorrectingPars.klen << std::endl;
	}

	if(std::numeric_limits<uint32_t>::max() == subRegionComboPars.includeFromSubRegionSize){
		subRegionComboPars.includeFromSubRegionSize = graphCorrectingPars.klen;
	}
	seqInfo refCorrectedInfo = refSeq;
	if(graphCorrectingPars.correctionOccurenceCutOff > 0){
		auto correctedSeqs = ContigsCompareGraphDev::correctSeqsByGraph(seqs, graphCorrectingPars);
		for(const auto & cseq : correctedSeqs){
			if(cseq.name_ == refSeq.name_){
				refCorrectedInfo = cseq;
				break;
			}
		}
		SeqOutput::write(correctedSeqs, SeqIOOptions::genFastaOutGz(njh::files::make_path(setUp.pars_.directoryName_, "correctedSeqs")));
		seqs = correctedSeqs;

		if(lenFilter){
			seqs = ContigsCompareGraphDev::correctSeqsByGraph(seqs, graphCorrectingPars);
		}
		graphCorrectingPars.klen = ContigsCompareGraphDev::findMinNonredundantKmer(minimumKlen, seqs, __PRETTY_FUNCTION__);
		{
			OutputStream klenOut(njh::files::make_path(setUp.pars_.directoryName_, "nonRedundantKlenAfterCorrection.txt"));
			klenOut << graphCorrectingPars.klen << std::endl;
		}
		if(std::numeric_limits<uint32_t>::max() == subRegionComboPars.includeFromSubRegionSize){
			subRegionComboPars.includeFromSubRegionSize = graphCorrectingPars.klen;
		}
	}


	//diversity
	{
		OutputStream outDivMeasures(njh::files::make_path(setUp.pars_.directoryName_, "divMeasuresFullRegion.tab.txt"));
		outDivMeasures << "target\tlength\ttotalHaps\tuniqueHaps\the\texp3\texp4\texp5\tsinglets\tdoublets\teffectiveNumOfAlleles\tShannonEntropyE" << std::endl;

		//first full div
		{
			std::unordered_map<std::string, uint32_t> popCounts;
			auto divMeasures = PopGenCalculator::getGeneralMeasuresOfDiversityRawInput(seqs);
			outDivMeasures << refSeqLoc.createUidFromCoordsStrand()
													<< "\t" << refSeqLoc.getLen()
													<< "\t" << seqs.size()
													<< "\t" << divMeasures.alleleNumber_
													<< "\t" << divMeasures.heterozygostiy_
													<< "\t" << divMeasures.ploidy3_.expectedPolyClonal_
													<< "\t" << divMeasures.ploidy4_.expectedPolyClonal_
													<< "\t" << divMeasures.ploidy5_.expectedPolyClonal_
													<< "\t" << divMeasures.singlets_
													<< "\t" << divMeasures.doublets_
													<< "\t" << divMeasures.effectiveNumOfAlleles_
													<< "\t" << divMeasures.ShannonEntropyE_
													<< std::endl;
		}
	}


	std::unordered_map<std::string, uint32_t> seqKey;
	for(const auto seq: iter::enumerate(seqs)){
		seqKey[seq.element.name_] = seq.index;
	}
	//replace refseq with any error correction done
	refSeq = seqs[seqKey[refSeq.name_]];

	ContigsCompareGraphDev compGraph(graphCorrectingPars.klen);
	compGraph.setOccurenceCutOff(graphCorrectingPars.kmerOccurenceCutOff);
	for(const auto & seq : seqs){
		compGraph.increaseKCounts(seq.seq_);
	}
	compGraph.populateNodesFromCounts();
	for(const auto & seq : seqs){
		compGraph.threadThroughSequence(seq, seq.name_);
	}

	if(setUp.pars_.debug_){
		OutputStream out(njh::files::make_path(setUp.pars_.directoryName_, "initialGraph.dot"));
		compGraph.writeRectangleDotColorBySampleCount(out);

	}


	compGraph.collapseSingleLinkedPathsSameReads();
	OutOptions outOpts(njh::files::make_path(setUp.pars_.directoryName_, "graph"), ".dot");

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
	if(!graphCorrectingPars.doNotCollapseLowFreqNodes){
		uint32_t collapseCount = 0;
		while(compGraph.collapseLowFreqNodes(graphCorrectingPars.allowableError, graphCorrectingPars.lowFreqCutOff)){
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
		{
			auto outOptsCurrent = outOpts;
			outOptsCurrent.outFilename_ = njh::files::prependFileBasename(outOptsCurrent.outFilename_,"noLabels_final_");
			OutputStream out(outOptsCurrent);
			compGraph.writeRectangleDotColorBySampleCount(out, true);
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
		OutputStream groupContigCountsOut(njh::files::make_path(conservedRegionInfoDir, "groupContigCountsOut.tab.txt"));
		groupContigCountsOut << "group\tcount" << std::endl;
		for(const auto & group : contigCountsPerGroup){
			groupContigCountsOut << group.first << '\t' << group.second.size() << std::endl;
		}




		auto seqsDir = njh::files::makeDir(subsegmentInfoDir, njh::files::MkdirPar{"subSeqs"});
		OutputStream outDivMeasures(njh::files::make_path(subsegmentInfoDir, "divMeasuresPerSubRegion.tab.txt"));
		outDivMeasures << "target\tlength\ttotalHaps\tuniqueHaps\the\texp3\texp4\texp5\tsinglets\tdoublets\teffectiveNumOfAlleles\tShannonEntropyE" << std::endl;

		//first full div
		OutputStream sharedLocsOutRefSeq(njh::files::make_path(subsegmentInfoDir, njh::pasteAsStr("refSeqGenomicLocs.bed")));
		table outSubSeqs(VecStr{"target", "frontSeq", "endSeq"});

		aligner alignerObj(std::max(refSeq.seq_.size(), refCorrectedInfo.seq_.size()), gapScoringParameters(5,1,0,0,0,0), substituteMatrix::createDegenScoreMatrixLessN(2, -2));
		alignerObj.alignCacheGlobal(refSeq, refCorrectedInfo);
		alignerObj.rearrangeObjsGlobal(refSeq, refCorrectedInfo);

		auto transfromCorrectedRefPosition = [&alignerObj](
				uint32_t correctSeqPos) {
			return getRealPosForAlnPos(alignerObj.alignObjectA_.seqBase_.seq_,
					getAlnPosForRealPos(alignerObj.alignObjectB_.seqBase_.seq_,
							correctSeqPos));
		};

		auto modSubSegmentToCorrectedRef = [&transfromCorrectedRefPosition](
				Bed3RecordCore & region) {
			region.chromStart_ = transfromCorrectedRefPosition(region.chromStart_);
			region.chromEnd_ = transfromCorrectedRefPosition(region.chromEnd_ -1 ) + 1;
		};
		OutputStream conservedNodeOut(njh::files::make_path(conservedRegionInfoDir, "conservedNodesPerGroup.tab.txt"));
		conservedNodeOut << "group\tsubSegment\tsubSegmentID\tseqCount\ttailCount\theadCount" << std::endl;
		for(const auto & nodes : nodesWithFull){


			auto processedNodes = ContigsCompareGraphDev::processConservedNodesVec(nodes.second);


			//output seqs
			{
				std::vector<seqInfo> conservedSeqs;

				for(const auto & n : nodes.second){
					auto ID = processedNodes.subseqToIDKey[n->k_];

					conservedSeqs.emplace_back(seqInfo(ID, n->k_));
					conservedNodeOut << nodes.first
							<< "\t" << n->k_
							<< "\t" << ID
							<< "\t" << n->inReadNamesIdx_.size()
							<< "\t" << n->tailCount()
							<< "\t" << n->headCount()
							<< std::endl;
				}
				auto conservedSeqsOpts = SeqIOOptions::genFastaOut(njh::files::make_path(conservedRegionInfoDir, njh::pasteAsStr(nodes.first, "_conservedNodes.fasta")));
				SeqOutput::write(conservedSeqs, conservedSeqsOpts);
			}

			//output conserved seq key
			{
				OutputStream keyOut(njh::files::make_path(conservedRegionInfoDir, njh::pasteAsStr(nodes.first, "_conservedNodesKey.tab.txt")));
				keyOut << "seq\tID" << std::endl;
				for(const auto & seqKey : processedNodes.subseqToIDKey){
					keyOut << seqKey.first << "\t" << seqKey.second << std::endl;
				}
			}

			//output reference locations relative to current seq
			{
				OutputStream refLocsOut(njh::files::make_path(conservedRegionInfoDir, njh::pasteAsStr(nodes.first, "_ref_sharedLocs.bed")));
				for(auto loc : processedNodes.nameToSubSegPositions_raw.at(refSeq.name_)){
					loc.name_ = processedNodes.subseqToIDKey[loc.name_];
					refLocsOut << loc.toDelimStrWithExtra() << std::endl;
				}
			}
			//output reference locations relative to genomic location
			{
				OutputStream refLocsOut(njh::files::make_path(conservedRegionInfoDir, njh::pasteAsStr(nodes.first, "_ref_sharedLocs_genomic.bed")));
				for(auto loc : processedNodes.nameToSubSegPositions_raw.at(refSeq.name_)){
					modSubSegmentToCorrectedRef(loc);
					std::string newName = processedNodes.subseqToIDKey[loc.name_];
					loc = refSeqLoc.genBedRecordCore().adjustSubRegionToRelativePosition(loc);
					loc.name_ = newName;
					refLocsOut << loc.toDelimStrWithExtra() << std::endl;
				}
			}

			uint32_t klen_minus1 = graphCorrectingPars.klen -1;
			auto createSubRegion = [&refSeq,&refSeqLoc,&modSubSegmentToCorrectedRef](uint32_t start, uint32_t end, const std::string & varName,
					std::vector<Bed6RecordCore> & relative, std::vector<Bed6RecordCore> & genomic){
				Bed3RecordCore varSubRegion(refSeq.name_, start, end);
				relative.emplace_back(GenomicRegion(varSubRegion).genBedRecordCore());
				relative.back().name_ = varName;
				modSubSegmentToCorrectedRef(varSubRegion);
				Bed6RecordCore varRegion = refSeqLoc.genBedRecordCore().adjustSubRegionToRelativePosition(varSubRegion);
				varRegion.name_ = varName;
				genomic.emplace_back(varRegion);
			};

			auto renameSeq = [](seqInfo & subSeq, uint32_t start, uint32_t end, const std::string & varName){
				MetaDataInName seqMeta;
				if(MetaDataInName::nameHasMetaData(subSeq.name_)){
					seqMeta = MetaDataInName(subSeq.name_);
				}
				seqMeta.addMeta("SubRegUID", varName, true);
				seqMeta.addMeta("SeqStart", start, true);
				seqMeta.addMeta("SeqStop", end, true);
				seqMeta.resetMetaInName(subSeq.name_);
			};
			//output variable region info
			{
				//key1 = input name; key2 = variable region name; val = sequence;

				std::unordered_map<std::string, std::unordered_map<std::string, std::string>> variableRegions;

				auto seqsDir = njh::files::makeDir(conservedRegionInfoDir, njh::files::MkdirPar{"subSeqsVariableRegions"});
				OutputStream outDivMeasures(njh::files::make_path(conservedRegionInfoDir, "divMeasuresPerVarRegion.tab.txt"));
//				outDivMeasures << "target\tlength\ttotalHaps\tuniqueHaps\the\texp3\texp4\texp5\tsinglets\tdoublets\teffectiveNumOfAlleles\tShannonEntropyE" << std::endl;
				outDivMeasures << njh::conToStr(calcPopMeasuresPars.genHeader(), "\t") << std::endl;
				OutputStream variableRegionsRelOut(njh::files::make_path(conservedRegionInfoDir, njh::pasteAsStr(nodes.first, "_ref_variable.bed")));
				OutputStream variableRegionsGenomicOut(njh::files::make_path(conservedRegionInfoDir, njh::pasteAsStr(nodes.first, "_ref_variable_genomic.bed")));
				std::vector<Bed6RecordCore> variableRegionsRelative;
				std::vector<Bed6RecordCore> variableRegionsGenomic;

				OutputStream variableRegionsRelExpanedOut(njh::files::make_path(conservedRegionInfoDir, njh::pasteAsStr(nodes.first, "_ref_variable_expaned.bed")));
				OutputStream variableRegionsGenomicExpanedOut(njh::files::make_path(conservedRegionInfoDir, njh::pasteAsStr(nodes.first, "_ref_variable_expanded_genomic.bed")));
				std::vector<Bed6RecordCore> variableRegionsRelativeExpaned;
				std::vector<Bed6RecordCore> variableRegionsGenomicExpaned;


				std::map<std::string, std::vector<seqInfo>> subRegionsSeqs;
				std::map<std::string, seqInfo> subRegionsRefSeqs;
				std::map<std::string, GenomicRegion> subRegionsRefSeqsLoc;
				std::map<std::string, CollapsedHaps> subRegionsSeqsUniq;
				//key1 = var region name; key2 = sequence; val = unique seq key;
				std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>> varRegionSeqKey;

				uint32_t varCount = 0;
				uint32_t totalVar = processedNodes.nameToSubSegPositions_filt.at(refSeq.name_).size() - 1;
				if(0 != processedNodes.nameToSubSegPositions_filt.at(refSeq.name_).front().chromStart_){
					++totalVar;
				}
				if(len(refCorrectedInfo) !=  processedNodes.nameToSubSegPositions_filt.at(refSeq.name_).back().chromEnd_){
					++totalVar;
				}
//				std::cout << __FILE__ << " " << __LINE__ << std::endl;



				//check front
				if(0 != processedNodes.nameToSubSegPositions_filt.at(refSeq.name_).front().chromStart_){
					auto varName = njh::pasteAsStr("var.", njh::leftPadNumStr<uint32_t>(varCount, totalVar));;
					varName = njh::pasteAsStr(refSeqLoc.createUidFromCoordsStrand(), "__", varName);
					uint32_t start = 0;
					uint32_t end = processedNodes.nameToSubSegPositions_filt.at(refSeq.name_).front().chromStart_;
					createSubRegion(start, end, varName, variableRegionsRelative, variableRegionsGenomic);
					//adjust to expanded
					//start -= klen_minus1;
					end += klen_minus1;
					createSubRegion(start, end, varName, variableRegionsRelativeExpaned, variableRegionsGenomicExpaned);
					subRegionsRefSeqsLoc[varName] = variableRegionsGenomicExpaned.back();
					++varCount;
					//getting diversity
					std::vector<seqInfo> subSeqs;
					for(const auto & endRegion : processedNodes.subSeqToNameToPos[processedNodes.nameToSubSegPositions_filt.at(refSeq.name_).front().name_]){
						uint32_t start = 0;
						uint32_t end = endRegion.second.chromStart_ + klen_minus1;
						auto subSeq = seqs[seqKey[endRegion.first]].getSubRead(start, end);
						if(subSeq.name_ == refSeqName){
							subRegionsRefSeqs[varName] = subSeq;
						}
						variableRegions[subSeq.name_][varName] = subSeq.seq_;
						renameSeq(subSeq, start, end, varName);
						subSeqs.emplace_back(subSeq);
					}
					subRegionsSeqs[varName] = subSeqs;
				}
//				std::cout << __FILE__ << " " << __LINE__ << std::endl;
				//middle
				if(processedNodes.nameToSubSegPositions_filt.at(refSeq.name_).size() > 1){
//					std::cout << "processedNodes.nameToSubSegPositions_filt.at(refSeq.name_).size(): " << processedNodes.nameToSubSegPositions_filt.at(refSeq.name_).size() << std::endl;
					for(auto pos : iter::range(processedNodes.nameToSubSegPositions_filt.at(refSeq.name_).size() - 1)){
//						std::cout << "pos: " << pos << std::endl;
						auto varName = njh::pasteAsStr("var.", njh::leftPadNumStr<uint32_t>(varCount, totalVar));
						varName = njh::pasteAsStr(refSeqLoc.createUidFromCoordsStrand(), "__", varName);
//						std::cout << varName << std::endl;
						uint32_t start = processedNodes.nameToSubSegPositions_filt.at(refSeq.name_)[pos].chromEnd_;
						uint32_t end = processedNodes.nameToSubSegPositions_filt.at(refSeq.name_)[pos + 1].chromStart_;
						if(end <= start){
							std::stringstream ss;
							ss << __PRETTY_FUNCTION__ << ", error " << "end: " << end << " is less than start: " << start << " for " << refSeqName << "\n";
							throw std::runtime_error{ss.str()};
						}
						if(end > refCorrectedInfo.seq_.size()){
							std::stringstream ss;
							ss << __PRETTY_FUNCTION__ << ", error " << "end: " << end << " is past the end of the seq len: " << refCorrectedInfo.seq_.size() << "\n";
							throw std::runtime_error{ss.str()};
						}
//						std::cout << "start: " << start << std::endl;
//						std::cout << "end: " << end << std::endl;
						createSubRegion(start, end, varName, variableRegionsRelative, variableRegionsGenomic);
						//adjust to expanded
						start -= klen_minus1;
						end += klen_minus1;
						createSubRegion(start, end, varName, variableRegionsRelativeExpaned, variableRegionsGenomicExpaned);
						subRegionsRefSeqsLoc[varName] = variableRegionsGenomicExpaned.back();
//						Bed3RecordCore varSubRegion(refSeq.name_, start, end);
//						variableRegionsRelative.emplace_back(GenomicRegion(varSubRegion).genBedRecordCore());
//						variableRegionsRelative.back().name_ = varName;
//						modSubSegmentToCorrectedRef(varSubRegion);
//						Bed6RecordCore varRegion = refSeqLoc.genBedRecordCore().adjustSubRegionToRelativePosition(varSubRegion);
//						varRegion.name_ = varName;
//						variableRegionsGenomic.emplace_back(varRegion);
						++varCount;
//						std::cout << __FILE__ << " " << __LINE__ << std::endl;
						//getting diversity
						std::vector<seqInfo> subSeqs;
						for(const auto & frontRegion : processedNodes.subSeqToNameToPos[processedNodes.nameToSubSegPositions_filt.at(refSeq.name_)[pos].name_]){

							auto endRegion = processedNodes.subSeqToNameToPos[processedNodes.nameToSubSegPositions_filt.at(refSeq.name_)[pos + 1].name_][frontRegion.first];
//							std::cout << "seqs[seqKey[frontRegion.first]].name_: " << seqs[seqKey[frontRegion.first]].name_ << std::endl;
							uint32_t start = frontRegion.second.chromEnd_ - klen_minus1;
							uint32_t end = endRegion.chromStart_ + klen_minus1;
							if(end <= start){
								std::stringstream ss;
								ss << __PRETTY_FUNCTION__ << ", error " << "end: " << end << " is less than start: " << start << " for " << seqs[seqKey[frontRegion.first]].name_ << "\n";
								throw std::runtime_error{ss.str()};
							}
							if(end > seqs[seqKey[frontRegion.first]].seq_.size()){
								std::stringstream ss;
								ss << __PRETTY_FUNCTION__ << ", error " << "end: " << end << " is past the end of the seq len: " << seqs[seqKey[frontRegion.first]].seq_.size() << "\n";
								throw std::runtime_error{ss.str()};
							}
//							std::cout << "start: " << start << std::endl;
//							std::cout << "end: " << end << std::endl;

							auto subSeq = seqs[seqKey[frontRegion.first]].getSubRead(start, end - start);
							if(subSeq.name_ == refSeqName){
								subRegionsRefSeqs[varName] = subSeq;
							}
							variableRegions[subSeq.name_][varName] = subSeq.seq_;
							renameSeq(subSeq, start, end, varName);
							subSeqs.emplace_back(subSeq);
						}
						subRegionsSeqs[varName] = subSeqs;
					}
				}
//				std::cout << __FILE__ << " " << __LINE__ << std::endl;
				//check back
				if(len(refCorrectedInfo) !=  processedNodes.nameToSubSegPositions_filt.at(refSeq.name_).back().chromEnd_){
					auto varName = njh::pasteAsStr("var.", njh::leftPadNumStr<uint32_t>(varCount, totalVar));
					varName = njh::pasteAsStr(refSeqLoc.createUidFromCoordsStrand(), "__", varName);
					uint32_t start = processedNodes.nameToSubSegPositions_filt.at(refSeq.name_).back().chromEnd_;
					uint32_t end = len(refCorrectedInfo);
					createSubRegion(start, end, varName, variableRegionsRelative, variableRegionsGenomic);
					//adjust to expanded
					start -= klen_minus1;
					//end += klen_minus1;
					createSubRegion(start, end, varName, variableRegionsRelativeExpaned, variableRegionsGenomicExpaned);
					subRegionsRefSeqsLoc[varName] = variableRegionsGenomicExpaned.back();
					++varCount;

					//getting diversity
					std::vector<seqInfo> subSeqs;
					for(const auto & frontRegion : processedNodes.subSeqToNameToPos[processedNodes.nameToSubSegPositions_filt.at(refSeq.name_).back().name_]){
						uint32_t start = frontRegion.second.chromEnd_ - klen_minus1;
						uint32_t end = len(seqs[seqKey[frontRegion.first]]);
						auto subSeq = seqs[seqKey[frontRegion.first]].getSubRead(start);
						if(subSeq.name_ == refSeqName){
							subRegionsRefSeqs[varName] = subSeq;
						}
						variableRegions[subSeq.name_][varName] = subSeq.seq_;
						renameSeq(subSeq, start, end, varName);
						subSeqs.emplace_back(subSeq);
					}
					subRegionsSeqs[varName] = subSeqs;
				}
//				std::cout << __FILE__ << " " << __LINE__ << std::endl;
				//genomic out
				for(auto & varRegion : variableRegionsGenomic){
					variableRegionsGenomicOut << varRegion.toDelimStrWithExtra() << std::endl;
				}
				for(auto & varRegion : variableRegionsGenomicExpaned){
					variableRegionsGenomicExpanedOut << varRegion.toDelimStrWithExtra() << std::endl;
				}
				//relative out
				for(auto & varRegion : variableRegionsRelative){
					variableRegionsRelOut << varRegion.toDelimStrWithExtra() << std::endl;
				}
				for(auto & varRegion : variableRegionsRelativeExpaned){
					variableRegionsRelExpanedOut << varRegion.toDelimStrWithExtra() << std::endl;
				}
				OutputStream outSubregionCodingKeyOut(njh::files::make_path(conservedRegionInfoDir, "subRegionCodedNameKey.tab.txt"));
				outSubregionCodingKeyOut << "variableRegName\tname\tnameID" << std::endl;
				//process info on variable regions
				{
					for(const auto & subSeqs : subRegionsSeqs){
						auto subSeqOpts = SeqIOOptions::genFastaOutGz(njh::files::make_path(seqsDir, subSeqs.first + ".fasta.gz"));
						SeqOutput::write(subSeqs.second, subSeqOpts);
						auto uniqSeqs = CollapsedHaps::collapseReads(subSeqs.second);
						//auto identifier = njh::pasteAsStr(refSeqLoc.createUidFromCoordsStrand(), "__", subSeqs.first);
						auto identifier = subSeqs.first;
						uniqSeqs.renameBaseOnFreq(identifier);

						for(const auto pos : uniqSeqs.getOrderByTopCnt()){
							varRegionSeqKey[subSeqs.first][uniqSeqs.seqs_[pos]->seq_] = pos;
							outSubregionCodingKeyOut << subSeqs.first
									<< "\t" << uniqSeqs.seqs_[pos]->name_
									<< "\t" << pos << std::endl;
						}
						subRegionsSeqsUniq[subSeqs.first] = uniqSeqs;
						SeqOutput::write(uniqSeqs.seqs_, SeqIOOptions::genFastaOutGz(njh::files::make_path(seqsDir, "uniq_" + subSeqs.first + ".fasta.gz")));
						//samples names
						auto sampNamesPerSeq = uniqSeqs.getSampleNamesPerSeqs();
						auto allSamples = uniqSeqs.getAllSampleNames();
						uint64_t maxLen = readVec::getMaxLength(uniqSeqs.seqs_);
						readVec::getMaxLength(refSeq, maxLen);

						std::shared_ptr<aligner> alignerObj = std::make_shared<aligner>(maxLen, gapScoringParameters(5,1,0,0,0,0), substituteMatrix(2,-2), false);
						alignerObj->weighHomopolymers_ = false;
						alignerObj->processAlnInfoInput(setUp.pars_.alnInfoDirName_, setUp.pars_.verbose_);
						//set up variant info
						auto idSeq = refSeq;
						idSeq.name_ = identifier;
						TranslatorByAlignment::VariantsInfo varInfo(subRegionsRefSeqsLoc[subSeqs.first].genBed3RecordCore(), idSeq);
						//get variant info
						auto refComps = uniqSeqs.getCompsAgainstRef(refSeq, *alignerObj, numThreads);
						for(const auto pos :  iter::range(uniqSeqs.size())){
							varInfo.addVariantInfo(
									refComps[pos].refAlnSeq_,
									refComps[pos].queryAlnSeq_,
									uniqSeqs.seqs_[pos]->cnt_,
									sampNamesPerSeq[pos],
									refComps[pos].comp_,
									subRegionsRefSeqsLoc[subSeqs.first].start_);
						}
						varInfo.setFinals(variantCallerRunPars);
						{
							calcPopMeasuresPars.numSegSites_ = varInfo.getFinalNumberOfSegratingSites();
							auto divMeasures = uniqSeqs.getGeneralMeasuresOfDiversity(
									calcPopMeasuresPars, alignerObj);
//							OutputStream divMeasuresOut(outOpts);
							outDivMeasures << njh::conToStr(divMeasures.getOut(uniqSeqs, identifier, calcPopMeasuresPars), "\t") << std::endl;
//							divMeasures.writeDivMeasures(
//									njh::files::make_path(conservedRegionInfoDir,
//											subSeqs.first + "_divMeasures.tab.txt"), uniqSeqs, identifier, calcPopMeasuresPars);
						}


//						auto divMeasures = PopGenCalculator::getGeneralMeasuresOfDiversityRawInput(subSeqs.second);
//						outDivMeasures << subSeqs.first
//								<< "\t" << vectorMean(readVec::getLengths(subSeqs.second))
//																<< "\t" << subSeqs.second.size()
//																<< "\t" << divMeasures.alleleNumber_
//																<< "\t" << divMeasures.heterozygostiy_
//																<< "\t" << divMeasures.ploidy3_.expectedCOIForPloidy_[3]
//																<< "\t" << divMeasures.ploidy4_.expectedCOIForPloidy_[4]
//																<< "\t" << divMeasures.ploidy5_.expectedCOIForPloidy_[5]
//																<< "\t" << divMeasures.singlets_
//																<< "\t" << divMeasures.doublets_
//																<< "\t" << divMeasures.effectiveNumOfAlleles_
//																<< "\t" << divMeasures.ShannonEntropyE_
//																<< '\n';
					}
				}
				OutputStream outSubregionCoding(njh::files::make_path(conservedRegionInfoDir, "subRegionCoded.tab.txt"));
				outSubregionCoding << "name\t" << njh::conToStr(njh::getVecOfMapKeys(subRegionsSeqs), "\t") << std::endl;
				for(const auto & name : variableRegions){
					outSubregionCoding << name.first;
					for(const auto & subRegion : name.second){
						outSubregionCoding << "\t" << varRegionSeqKey[subRegion.first][subRegion.second];
					}
					outSubregionCoding << std::endl;
				}
			}


			if(!processedNodes.nameToSubSegPositions_filt.at(refSeq.name_).empty()){
				//getting the region with the last and first conserved regions
				Bed3RecordCore firstToLastShared(processedNodes.nameToSubSegPositions_filt.at(refSeq.name_).front());
				firstToLastShared.chromEnd_ = processedNodes.nameToSubSegPositions_filt.at(refSeq.name_).back().chromEnd_;
				modSubSegmentToCorrectedRef(firstToLastShared);
				auto firstToLastShared_genomic = refSeqLoc.genBedRecordCore().adjustSubRegionToRelativePosition(firstToLastShared);
				OutputStream refSharedLocsOut(njh::files::make_path(subsegmentInfoDir, njh::pasteAsStr(nodes.first, "_refFirstToLastSharedRegion.bed")));
				refSharedLocsOut << firstToLastShared_genomic.toDelimStr() << std::endl;
			}
			if(processedNodes.nameToSubSegPositions_filt.at(refSeq.name_).size() > 1 ){
				//getting the region with the last and first conserved regions
				BedUtility::SubRegionCombo subReg(processedNodes.nameToSubSegPositions_filt.at(refSeq.name_).front(), processedNodes.nameToSubSegPositions_filt.at(refSeq.name_).back());
				auto bedReg = subReg.genOut(subRegionComboPars.includeFromSubRegionSize);
				modSubSegmentToCorrectedRef(bedReg);
				auto bedReg_genomic = refSeqLoc.genBedRecordCore().adjustSubRegionToRelativePosition(bedReg);
				OutputStream refSharedLocsOut(njh::files::make_path(subsegmentInfoDir, njh::pasteAsStr(nodes.first, "_refFirstToLastSharedRegionSlimmer.bed")));
				refSharedLocsOut << bedReg_genomic.toDelimStrWithExtra() << std::endl;
			}

			std::unordered_map<std::string, uint32_t>lengths;
			lengths[refSeq.name_] = len(refCorrectedInfo);
			std::vector<BedUtility::SubRegionCombo> subRegions = BedUtility::genSubRegionCombos(processedNodes.nameToSubSegPositions_filt.at(refSeq.name_), lengths, subRegionComboPars);


			OutputStream refSubRegionsLocsRelOut(njh::files::make_path(subsegmentInfoDir, njh::pasteAsStr(nodes.first, "_ref_subRegions.bed")));


			//key = sub seq UID, val = 1st node seq start, 2nd node seq stop
			std::unordered_map<std::string, std::pair<std::string, std::string>> subPortionsByConservedSeq;
			//
			if(processedNodes.nameToSubSegPositions_filt.at(refSeq.name_).size() > 1 ){
				BedUtility::SubRegionCombo subReg(processedNodes.nameToSubSegPositions_filt.at(refSeq.name_).front(), processedNodes.nameToSubSegPositions_filt.at(refSeq.name_).back());
				auto bedReg = subReg.genOut(subRegionComboPars.includeFromSubRegionSize);
				modSubSegmentToCorrectedRef(bedReg);
				auto bedReg_genomic = refSeqLoc.genBedRecordCore().adjustSubRegionToRelativePosition(bedReg);

				auto startSubReg = subReg.genSubStart(subRegionComboPars.includeFromSubRegionSize);;
				auto endSubReg   = subReg.genSubEnd(subRegionComboPars.includeFromSubRegionSize);;
				seqInfo startSubSeq = refCorrectedInfo.getSubRead(startSubReg.chromStart_, startSubReg.length());
				seqInfo endSubSeq = refCorrectedInfo.getSubRead(endSubReg.chromStart_, endSubReg.length());

//				seqInfo startSubSeqFull = refCorrectedInfo.getSubRead(subReg.startRegion_.chromStart_, subReg.startRegion_.length());
//				seqInfo endSubSeqFull = refCorrectedInfo.getSubRead(subReg.endRegion_.chromStart_, subReg.endRegion_.length());


				endSubSeq.reverseComplementRead(false, true);
				outSubSeqs.addRow(bedReg_genomic.name_, startSubSeq.seq_, endSubSeq.seq_);
				std::vector<seqInfo> seqsFromLargestSubRegion;
				for(const auto & frontRegion : processedNodes.subSeqToNameToPos[subReg.startRegion_.name_]){
					auto startRegion = frontRegion.second;
					auto endRegion = njh::mapAt(njh::mapAt(processedNodes.subSeqToNameToPos,subReg.endRegion_.name_), frontRegion.first);
					BedUtility::SubRegionCombo subRegionForSeq(startRegion, endRegion);
					auto outRegion = subRegionForSeq.genOut(subRegionComboPars.includeFromSubRegionSize);
	//						std::cout << __FILE__ << " " << __LINE__ << std::endl;
	//						std::cout << "frontRegion.first: " << frontRegion.first << std::endl;
	//						std::cout << "seqKey[frontRegion.first]: " << seqKey[frontRegion.first] << std::endl;
	//						std::cout << "outRegion.chromStart_: " << outRegion.chromStart_ << std::endl;
	//						std::cout << "outRegion.length: " << outRegion.length() << std::endl;
	//						std::cout << "seqs[seqKey[frontRegion.first]].size(): " << seqs[seqKey[frontRegion.first]].seq_.size()<< std::endl;

					seqInfo subSeq(seqs[seqKey[frontRegion.first]].getSubRead(outRegion.chromStart_, outRegion.length()));
	//						std::cout << __FILE__ << " " << __LINE__ << std::endl;
					renameSeq(subSeq, outRegion.chromStart_, outRegion.chromEnd_, bedReg_genomic.name_);
//					MetaDataInName seqMeta(subSeq.name_);
//					seqMeta.addMeta("SubRegUID", bedReg_genomic.name_, true);
//					seqMeta.addMeta("SeqStart", outRegion.chromStart_, true);
//					seqMeta.addMeta("SeqStop", outRegion.chromEnd_, true);
//	//						std::cout << __FILE__ << " " << __LINE__ << std::endl;
//
//					seqMeta.resetMetaInName(subSeq.name_);
					seqsFromLargestSubRegion.emplace_back(subSeq);
				}
				// diversity
				{
					auto divMeasures = PopGenCalculator::getGeneralMeasuresOfDiversityRawInput(seqsFromLargestSubRegion);
					outDivMeasures << bedReg_genomic.name_
							<< "\t" << bedReg_genomic.length()
															<< "\t" << seqsFromLargestSubRegion.size()
															<< "\t" << divMeasures.alleleNumber_
															<< "\t" << divMeasures.heterozygostiy_
															<< "\t" << divMeasures.ploidy3_.expectedCOIForPloidy_[3]
															<< "\t" << divMeasures.ploidy4_.expectedCOIForPloidy_[4]
															<< "\t" << divMeasures.ploidy5_.expectedCOIForPloidy_[5]
															<< "\t" << divMeasures.singlets_
															<< "\t" << divMeasures.doublets_
															<< "\t" << divMeasures.effectiveNumOfAlleles_
															<< "\t" << divMeasures.ShannonEntropyE_
															<< std::endl;
				}
			}

//			std::cout << __FILE__ << " " << __LINE__ << std::endl;
			std::vector<Bed6RecordCore> subRegionsGenomicLocs;
			if(processedNodes.nameToSubSegPositions_filt.at(refSeq.name_).size() > 1 ){
				for(const auto & subReg : subRegions){

					auto bedReg = subReg.genOut(subRegionComboPars.includeFromSubRegionSize);
					modSubSegmentToCorrectedRef(bedReg);
					auto genomicLoc = refSeqLoc.genBedRecordCore().adjustSubRegionToRelativePosition(bedReg);

					GenomicRegion gRegion(genomicLoc);
					auto uid = gRegion.createUidFromCoordsStrand();

//					seqInfo startSubSeqFull = refCorrectedInfo.getSubRead(subReg.startRegion_.chromStart_, subReg.startRegion_.length());
//					seqInfo endSubSeqFull = refCorrectedInfo.getSubRead(subReg.endRegion_.chromStart_, subReg.endRegion_.length());
					auto startSubReg = subReg.genSubStart(subRegionComboPars.includeFromSubRegionSize);;
					auto endSubReg   = subReg.genSubEnd(subRegionComboPars.includeFromSubRegionSize);;
					seqInfo startSubSeq = refCorrectedInfo.getSubRead(startSubReg.chromStart_, startSubReg.length());
					seqInfo endSubSeq = refCorrectedInfo.getSubRead(endSubReg.chromStart_, endSubReg.length());
					genomicLoc.name_ = uid;

					subPortionsByConservedSeq[uid] = std::make_pair(subReg.startRegion_.name_, subReg.endRegion_.name_);
					endSubSeq.reverseComplementRead(false, true);
					outSubSeqs.addRow(uid, startSubSeq.seq_, endSubSeq.seq_);
					subRegionsGenomicLocs.emplace_back(genomicLoc);
					bedReg.name_ = njh::pasteAsStr(
							processedNodes.subseqToIDKey[subReg.startRegion_.name_], "-",
							processedNodes.subseqToIDKey[subReg.endRegion_.name_]);
					refSubRegionsLocsRelOut << bedReg.toDelimStr() << std::endl;
				}
//				std::cout << __FILE__ << " " << __LINE__ << std::endl;

				for(const auto & subReg : subRegionsGenomicLocs){
					sharedLocsOutRefSeq << subReg.toDelimStr()<< std::endl;
				}
//				std::cout << __FILE__ << " " << __LINE__ << std::endl;
				//cut input to subsegments
				std::unordered_map<std::string, std::vector<seqInfo>> subRegionsSeqs;
	//			std::cout << "subPortionsByConservedSeq.size(): " << subPortionsByConservedSeq.size() << std::endl;
	//			std::cout << "processedNodes.subSeqToNameToPos.size(): " << processedNodes.subSeqToNameToPos.size() << std::endl;
	//			std::cout << njh::conToStr(getVectorOfMapKeys(processedNodes.subSeqToNameToPos), "\n") << std::endl;

				for(const auto & subRegion : subPortionsByConservedSeq){
	//				std::cout << subRegion.second.first << std::endl;

					for(const auto & frontRegion : processedNodes.subSeqToNameToPos[subRegion.second.first]){
	//					std::cout << "\t" << subRegion.second.first << std::endl;
	//					std::cout << "\t" << subRegion.second.second << std::endl;
//						std::cout << __FILE__ << " " << __LINE__ << std::endl;

						auto startRegion = frontRegion.second;
						auto endRegion = njh::mapAt(njh::mapAt(processedNodes.subSeqToNameToPos,subRegion.second.second), frontRegion.first);
						BedUtility::SubRegionCombo subRegionForSeq(startRegion, endRegion);
						auto outRegion = subRegionForSeq.genOut(subRegionComboPars.includeFromSubRegionSize);
//						std::cout << __FILE__ << " " << __LINE__ << std::endl;
//						std::cout << "frontRegion.first: " << frontRegion.first << std::endl;
//						std::cout << "seqKey[frontRegion.first]: " << seqKey[frontRegion.first] << std::endl;
//						std::cout << "outRegion.chromStart_: " << outRegion.chromStart_ << std::endl;
//						std::cout << "outRegion.length: " << outRegion.length() << std::endl;
//						std::cout << "seqs[seqKey[frontRegion.first]].size(): " << seqs[seqKey[frontRegion.first]].seq_.size()<< std::endl;

						seqInfo subSeq(seqs[seqKey[frontRegion.first]].getSubRead(outRegion.chromStart_, outRegion.length()));
//						std::cout << __FILE__ << " " << __LINE__ << std::endl;
						renameSeq(subSeq, outRegion.chromStart_, outRegion.chromEnd_, subRegion.first);
//						MetaDataInName seqMeta(subSeq.name_);
//						seqMeta.addMeta("SubRegUID", subRegion.first, true);
//						seqMeta.addMeta("SeqStart", outRegion.chromStart_, true);
//						seqMeta.addMeta("SeqStop", outRegion.chromEnd_, true);
//						std::cout << __FILE__ << " " << __LINE__ << std::endl;
//						seqMeta.resetMetaInName(subSeq.name_);
						subRegionsSeqs[subRegion.first].emplace_back(subSeq);
					}
				}
	//			std::cout << "subRegionsSeqs.size(): " << subRegionsSeqs.size() << std::endl;
//				std::cout << __FILE__ << " " << __LINE__ << std::endl;

				for(const auto & subSeqs : subRegionsSeqs){
					auto subSeqOpts = SeqIOOptions::genFastaOutGz(njh::files::make_path(seqsDir, subSeqs.first + ".fasta.gz"));
					SeqOutput::write(subSeqs.second, subSeqOpts);
					auto divMeasures = PopGenCalculator::getGeneralMeasuresOfDiversityRawInput(subSeqs.second);
					outDivMeasures << subSeqs.first
															<< "\t" << vectorMean(readVec::getLengths(subSeqs.second))
															<< "\t" << subSeqs.second.size()
															<< "\t" << divMeasures.alleleNumber_
															<< "\t" << divMeasures.heterozygostiy_
															<< "\t" << divMeasures.ploidy3_.expectedCOIForPloidy_[3]
															<< "\t" << divMeasures.ploidy4_.expectedCOIForPloidy_[4]
															<< "\t" << divMeasures.ploidy5_.expectedCOIForPloidy_[5]
															<< "\t" << divMeasures.singlets_
															<< "\t" << divMeasures.doublets_
															<< "\t" << divMeasures.effectiveNumOfAlleles_
															<< "\t" << divMeasures.ShannonEntropyE_
															<< '\n';
				}
			}
		}
		OutputStream outSubSeqsOut(njh::files::make_path(subsegmentInfoDir, "startEndSeqs.tab.txt"));
		outSubSeqs.outPutContents(outSubSeqsOut, "\t");


	}
	return 0;
}




} // namespace njhseq

