/*
 * kmerExp_add.cpp
 *
 *  Created on: Dec 18, 2016
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
#include "elucidator/objects/BioDataObject.h"



namespace njhseq {





int kmerExpRunner::convertKmerSearchToBinaryMatrix(const njh::progutils::CmdArgs & inputCommands){
	bfs::path inputFnp = "";
	OutOptions outOpts(bfs::path("outmat.tab.txt"));
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processSeq(true);
	setUp.setOption(inputFnp, "--inputFnp", "Input Fnp", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	std::set<std::string> inputSeqNames;
	{
		Bed6RecordCore bRecord;
		BioDataFileIO<Bed6RecordCore> bedReader(IoOptions{InOptions{inputFnp}});
		bedReader.openIn();
		while(bedReader.readNextRecord(bRecord)){
			inputSeqNames.insert(bRecord.chrom_);
		}
	}

	std::unordered_map<std::string, uint32_t> seqIndex;
	{
		uint32_t pos = 0;
		for(const auto name : inputSeqNames){
			seqIndex[name] = pos;
			++pos;
		}
	}


	Bed6RecordCore bRecord;
	BioDataFileIO<Bed6RecordCore> bedReader(IoOptions{InOptions{inputFnp}});
	bedReader.openIn();
	std::vector<std::vector<uint16_t>> occurenceMatrix(inputSeqNames.size(),
			std::vector<uint16_t>(setUp.pars_.seqObj_.seqBase_.seq_.size()));

	while(bedReader.readNextRecord(bRecord)){
		uint32_t pos = njh::StrToNumConverter::stoToNum<uint32_t>(bRecord.name_);
		for(const auto refSeqPos : iter::range<uint32_t>(pos, pos + bRecord.score_)){
			occurenceMatrix[seqIndex[bRecord.chrom_]][refSeqPos] = 1;
		}
	}

	OutputStream out(outOpts);
	for (const auto & inputName : inputSeqNames) {
		out << inputName << "\t"
				<< njh::conToStr(occurenceMatrix[seqIndex[inputName]], "\t")
				<< std::endl;
	}

	return 0;
}

int kmerExpRunner::kmerSearch(const njh::progutils::CmdArgs & inputCommands){
	uint32_t kLen = 11;
	bool noReverse = false;
	bfs::path genome = "";
	OutOptions outOpts(bfs::path("out.tab.txt"));
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processSeq(true);
	setUp.setOption(kLen, "--kLen", "Kmer Length");
	setUp.setOption(noReverse, "--noReverse", "No Reverse");
	uint32_t minStreakSize = kLen + 1;
	setUp.setOption(minStreakSize, "--minStreakSize", "Min Streak Size");
	setUp.setOption(genome, "--genome", "Genome Fasta file", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	if(setUp.pars_.debug_){
		setUp.pars_.seqObj_.seqBase_.outPutSeq(std::cout);
	}

	GenomeSeqSearch searcher(outOpts);

	auto genomeOpts = SeqIOOptions::genFastaIn(genome.string());

	std::unordered_map<std::string, std::unique_ptr<kmerInfo>> genomeInfos;
	seqInfo seq;
	njh::stopWatch watch;
	watch.setLapName("Genome Index");
	SeqInput reader(genomeOpts);
	reader.openIn();
	while(reader.readNextRead(seq)){
		auto chromName = seq.name_;
		trimAtFirstWhitespace(chromName);
		genomeInfos[chromName] = std::make_unique<kmerInfo>(seq.seq_, kLen, !noReverse);
	}
	watch.startNewLap("Forward search");
	for (const auto & seqPos : iter::range<size_t>(0,
			len(setUp.pars_.seqObj_) - kLen + 1)) {
		auto currentK = setUp.pars_.seqObj_.seqBase_.seq_.substr(seqPos, kLen);
		if(setUp.pars_.verbose_){
			std::cout << "\r" << seqPos << " : " << currentK;
			std::cout.flush();
		}
		searcher.growStreaks(currentK, seqPos, false, genomeInfos);
		searcher.pruneStreaks(minStreakSize);
	}
	searcher.purgeAllStreaks(minStreakSize);
	if(setUp.pars_.verbose_){
		std::cout << std::endl;
	}
	if(!noReverse){
		setUp.pars_.seqObj_.seqBase_.reverseComplementRead(false, true);
		for (const auto & seqPos : iter::range<size_t>(0,
				len(setUp.pars_.seqObj_) - kLen + 1)) {
			auto currentK = setUp.pars_.seqObj_.seqBase_.seq_.substr(seqPos, kLen);
			if(setUp.pars_.verbose_){
				std::cout << "\r" << seqPos << " : " << currentK;
				std::cout.flush();
			}

			uint32_t forwardSeqPos = len(setUp.pars_.seqObj_) - kLen - seqPos;
			searcher.growStreaks(currentK, forwardSeqPos, true, genomeInfos);
			searcher.pruneStreaks(minStreakSize);
		}
		searcher.purgeAllStreaks(minStreakSize);
	}
	watch.startNewLap("Reverse search");



	if(setUp.pars_.verbose_){
		std::cout << std::endl;
		watch.logLapTimes(std::cout, true, 6, true);
	}

	return 0;
}

int kmerExpRunner::microsatsKmerSearch(const njh::progutils::CmdArgs & inputCommands){
	bfs::path genome = "";
	OutOptions outOpts(bfs::path("out.tab.txt"));
	size_t minRepeatCount = 3;
	size_t minLength = 6;
	std::string kLens = "2,3,4,5";
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(genome, "--genome", "Genome Fasta file", true);
	setUp.setOption(minRepeatCount, "--minRepeatCount", "Min Repeat Count");
	setUp.setOption(minLength, "--minLength", "Minimum Length in base pairs");
	setUp.setOption(kLens, "--kLens", "Kmer Lengths to search for");
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	std::vector<uint32_t> kLensNum = vecStrToVecNum<uint32_t>(tokenizeString(kLens, ","));

	auto genomeOpts = SeqIOOptions::genFastaIn(genome.string());

	auto isHomopolymer =
			[](const std::string & k) {
				return std::all_of(k.begin(), k.end(),[&k](const char c) {return k.front() == c;});
			};
	seqInfo seq;
	njh::stopWatch watch;
	std::ofstream out;
	outOpts.openFile(out);
	//out << "chrom\tstart\tend\tname\tscore\tstrand" << std::endl;
	for(auto kLen : kLensNum){
		if(setUp.pars_.verbose_){
			std::cout << "On klen: " << kLen << std::endl;
		}
		watch.startNewLap("Genome Index : " + estd::to_string(kLen));
		std::unordered_map<std::string, std::unique_ptr<kmerInfo>> genomeInfos;
		SeqInput reader(genomeOpts);
		reader.openIn();
		while(reader.readNextRead(seq)){
			auto chromName = seq.name_;
			trimAtFirstWhitespace(chromName);
			genomeInfos[chromName] = std::make_unique<kmerInfo>(seq.seq_, kLen, true);
		}
		watch.startNewLap("Search counts : " + estd::to_string(kLen));
		for(const auto & chrom : genomeInfos){
			for(const auto & k : chrom.second->kmers_){
				if(isHomopolymer(k.first)){
					continue;
				}
				if(k.second.positions_.size() <= 1){
					continue;
				}
				if(kLen > 3 && kLen % 2 == 0){
					if(k.first.substr(0,kLen/2) == k.first.substr(kLen/2)){
						continue;
					}
				}
				bool growing = false;
				size_t startPos = std::numeric_limits<size_t>::max();
				uint32_t repeatNumber = 1;
				for(const auto & pos : iter::range(k.second.positions_.size() - 1)){
					if( (k.second.positions_[pos + 1] - k.second.positions_[pos]) == kLen){
						++repeatNumber;
						if (!growing) {
							growing = true;
							startPos = k.second.positions_[pos];
						}
					}else{
						if(growing){
							if(repeatNumber >= minRepeatCount && k.first.size() * repeatNumber >=minLength){
								out << njh::conToStr(toVecStr(chrom.first, startPos, startPos + repeatNumber * kLen, k.first, repeatNumber, "+"), "\t") << std::endl;;//(, k.first, startPos, repeatNumber, repeatNumber * kLen);
							}
						}
						startPos = std::numeric_limits<size_t>::max();
						repeatNumber = 1;
						growing = false;
					}
				}
				if(growing){
					if(repeatNumber >= minRepeatCount && k.first.size() * repeatNumber >=minLength){
						out << njh::conToStr(toVecStr(chrom.first, startPos, startPos + repeatNumber * kLen, k.first, repeatNumber, "+"), "\t") << std::endl;;//(, k.first, startPos, repeatNumber, repeatNumber * kLen);
					}
				}
			}
		}
	}


	if(setUp.pars_.verbose_){
		std::cout << std::endl;
		watch.logLapTimes(std::cout, true, 6, true);
	}
	return 0;
}



int kmerExpRunner::genomeKmerCompare(const njh::progutils::CmdArgs & inputCommands){
	uint32_t kLen = 10;
	uint32_t testNumber = 10;
	bfs::path genome1;
	bfs::path genome2;
	seqSetUp setUp(inputCommands);
	setUp.setOption(kLen, "--kLen", "Kmer Length");
	setUp.setOption(genome1, "--genome1", "Genome 1", true);
	setUp.setOption(genome2, "--genome2", "Genome 2", true);
	setUp.finishSetUp(std::cout);

	auto genome1Opts = SeqIOOptions::genFastaIn(genome1.string());
	auto genome2Opts = SeqIOOptions::genFastaIn(genome2.string());

	std::unordered_map<std::string, std::unique_ptr<kmerInfo>> genome1Infos;
	std::unordered_map<std::string, std::unique_ptr<kmerInfo>> genome2Infos;
	seqInfo seq;
	njh::stopWatch watch;
	watch.setLapName("Genome1 Index");
	{

		SeqInput reader(genome1Opts);
		reader.openIn();
		while(reader.readNextRead(seq)){
			genome1Infos[seq.name_] = std::make_unique<kmerInfo>(seq.seq_, kLen, true);
		}
	}

	{
		watch.startNewLap("Genome2 Index");
		SeqInput reader(genome2Opts);
		reader.openIn();
		while(reader.readNextRead(seq)){
			genome2Infos[seq.name_] = std::make_unique<kmerInfo>(seq.seq_, kLen, true);
		}
	}
	watch.startNewLap("Search");
	SeqInput reader(genome1Opts);
	reader.openIn();
	while(reader.readNextRead(seq)){
		uint32_t num = 0;
		for(const auto pos : iter::range(len(seq) - kLen + 1)){
			auto kmer = seq.seq_.substr(pos, kLen);
			bool found = false;
			for(const auto & kInfo : genome2Infos){
				if(kInfo.second->kmers_.end() != kInfo.second->kmers_.find(kmer) ||
						kInfo.second->kmersRevComp_.end() != kInfo.second->kmersRevComp_.find(kmer)){
					found = true;
					break;
				}
			}
			if(!found){
				std::cout << seq.name_ << "\t" << pos << "\t" << kmer << std::endl;
				++num;
				if(num >= testNumber){
					break;
				}
			}
		}
	}
	watch.logLapTimes(std::cout, true, 6,true);
	return 0;
}


int kmerExpRunner::findingMinimumKLenForNoRedundantKmers(const njh::progutils::CmdArgs & inputCommands){

	seqSetUp setUp(inputCommands);
  setUp.processVerbose();
  setUp.processReadInNames(true);
  setUp.finishSetUp(std::cout);


	uint64_t maxLen = 0;
	seqInfo seq;
	SeqInput reader(setUp.pars_.ioOptions_);

	reader.openIn();
	while (reader.readNextRead(seq)) {
		if(setUp.pars_.verbose_){
			std::cout << seq.name_ << std::endl;
		}
		readVec::getMaxLength(seq, maxLen);
	}


  uint32_t klen = 2;
  bool foundLength = false;
  while(klen < maxLen && !foundLength){
  		if(setUp.pars_.verbose_){
  			std::cout  << klen  << std::endl;
  		}
  		reader.reOpenIn();
  	  bool allPass = true;
  		while(reader.readNextRead(seq)){
  			kmerInfo kinfo(seq.seq_, klen, false);
  			for(const auto & k : kinfo.kmers_){
  				if(k.second.count_ > 1){
  					allPass = false;
  					break;
  				}
  			}
  	  }
  		if(allPass){
  			foundLength = true;
  		}else{
  			++klen;
  		}
  }

  if(foundLength){
  		std::cout << klen << std::endl;

  }else{
  		std::cout << "Found no kmer length that would recreate non-redundant kmers for all input seqs" << std::endl;
  }

  return 0;
}

int kmerExpRunner::writeKmerAccerlation(const njh::progutils::CmdArgs & inputCommands){
  seqSetUp setUp(inputCommands);
  uint32_t numThreads = 1;
  bool doTies = false;
  uint32_t kmerStart = 2;
  uint32_t kmerStop = 5;
	readDistGraph<double>::dbscanPars dbpars;
	dbpars.eps_ =  0.01;
	dbpars.minEpNeighbors_ = 5;
  bool useKNumber = false;
  uint32_t nodeSize = 50;
  uint32_t clustersizeCutOff = 0;
  bool useDbscan = false;
  setUp.processVerbose();
  setUp.setOption(useDbscan, "--useDbscan", "useDbscan");
  setUp.setOption(useKNumber, "--useKNumber", "Use number of decreasing kmers rather than distance metric");
  setUp.setOption(doTies, "--doTies", "Do ties");
  setUp.setOption(dbpars.eps_, "-cutOff", "cutOff");
  setUp.setOption(nodeSize, "--nodeSize", "nodeSize");
  setUp.setOption(clustersizeCutOff, "-clustersizeCutOff", "clustersizeCutOff");
  setUp.setOption(dbpars.minEpNeighbors_, "--minNeighbors", "minNeighbors");
  setUp.processDefaultReader(true);
  setUp.setOption(numThreads, "--numThreads,-t", "Number of Threads to Use");
  setUp.setOption(kmerStart, "--kmerLenStart", "Length for kmers to start at");
  setUp.setOption(kmerStop, "--kmerLenStop", "Length for kmers to stop at");
  setUp.processDirectoryOutputName(true);
  setUp.processRefFilename(false);
  setUp.processAlignerDefualts();

  setUp.finishSetUp(std::cout);
  setUp.startARunLog(setUp.pars_.directoryName_);
  njh::stopWatch watch;
  watch.setLapName("Index Kmers");
  auto reads = createKmerReadVec(setUp.pars_.ioOptions_, kmerStart, false);
	setUp.rLog_.logCurrentTime("All by All Comparison");
	std::vector<std::vector<double>> distances = getKmerAccerDistance(reads, kmerStart, kmerStop, numThreads, useKNumber, setUp.pars_.verbose_);
	setUp.rLog_.logCurrentTime("Graph Creation");
	readDistGraph<double> distanceGraph(distances, reads);

	if(useDbscan){
		distanceGraph.resetBestAndVistEdges();
		distanceGraph.resetVisitedNodes();
		distanceGraph.allDetermineLowestBest(doTies);
		distanceGraph.removeOffEdges();
		distanceGraph.dbscan(dbpars);
		distanceGraph.assignNoiseNodesAGroup();
	}else{
		distanceGraph.turnOffEdgesAbove(dbpars.eps_);
		distanceGraph.resetBestAndVistEdges();
		distanceGraph.resetVisitedNodes();
		distanceGraph.allDetermineLowestBest(doTies);
		distanceGraph.determineGroups();
	}


	if(setUp.pars_.debug_){
		distanceGraph.printAdjByGroup(std::cout);
	}
	std::vector<readObject> refSeqs;
	uint64_t maxReadLength = 0;
	if("" != setUp.pars_.refIoOptions_.firstName_){
		refSeqs = SeqInput::getReferenceSeq(setUp.pars_.refIoOptions_, maxReadLength);
	}
	//output visualization where clusters are written and colored either by a comparison to a reference sequence
	//or to the best matching consensus sequence
	//@todo change the color to the read it gets mapped to if remapping
	setUp.rLog_.logCurrentTime("Visualizing");
	bfs::path visDir = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar("vis"));
	Json::Value treeJson;


	if(!refSeqs.empty()){
		std::unordered_map<std::string, std::string> nameToColor;
		std::mutex nameToColorMut;
		std::vector<uint32_t> positions(reads.size());
		njh::iota(positions,0U);
		njh::concurrent::LockableQueue<uint32_t> positionsQueue(positions);

		//create aligner with gap info from ref so a different gap scoring can be used for mapping determination
		aligner alignerObjMapper(maxReadLength, setUp.pars_.gapInfoRef_,
				setUp.pars_.scoring_, KmerMaps(), setUp.pars_.qScorePars_,
				setUp.pars_.colOpts_.alignOpts_.countEndGaps_,
				setUp.pars_.colOpts_.iTOpts_.weighHomopolyer_);
		alignerObjMapper.processAlnInfoInput(setUp.pars_.alnInfoDirName_);
		concurrent::AlignerPool alignPool(alignerObjMapper, numThreads);
		alignPool.initAligners();
		alignPool.outAlnDir_ = setUp.pars_.outAlnInfoDirName_;

		auto compareReads = [&alignPool, &refSeqs,&nameToColorMut,&nameToColor,&reads,&positionsQueue](){
			uint32_t pos = std::numeric_limits<uint32_t>::max();
			std::unordered_map<std::string, std::string> nameToColorCurrent;
			auto currentAligner = alignPool.popAligner();
			while(positionsQueue.getVal(pos)){
				const auto & seq = reads[pos];
				auto popColors = getColorsForNames(readVec::getNames(refSeqs));
				auto bestRef = profiler::compareToRefSingle(refSeqs, seq, *currentAligner,
						false, false);
				auto bestRefToks = tokenizeString(bestRef.front(), ",");
				nameToColorCurrent[seq->seqBase_.name_] = popColors[bestRefToks[0]].getHexStr();
			}

			{
				std::lock_guard<std::mutex> lock(nameToColorMut);
				for(const auto & name : nameToColorCurrent){
					nameToColor[name.first] = name.second;
				}
			}
		};

		std::vector<std::thread> threads;
		for(uint32_t threadNum = 0; threadNum < numThreads; ++threadNum){
			threads.emplace_back(std::thread(compareReads));
		}
		for(auto & t : threads){
			t.join();
		}
		//auto gJson = distanceGraph.toJson(2);
		treeJson = distanceGraph.toJson(clustersizeCutOff, nameToColor);
	}else{
		treeJson = distanceGraph.toJson(clustersizeCutOff);
	}

	auto & nodes = treeJson["nodes"];
	for(auto & node : nodes){
		node["size"] = distanceGraph.nodes_[distanceGraph.nameToNodePos_[node["name"].asString()]]->value_->cnt_ * nodeSize;
	}

	std::ofstream treeJsonFile(visDir.string() + "tree.json");
	treeJsonFile << treeJson;

	std::ofstream treeHtmlFile(visDir.string() + "tree.html");
	genTreeHtml(treeHtmlFile, "tree.json", "tree.js");

	std::ofstream treeJsFile(visDir.string() + "tree.js");
	genSimpleTreeJs(treeJsFile);

	OutOptions groupFileOpts(njh::files::make_path(setUp.pars_.directoryName_, "groupInfo.tab.txt"));
	std::ofstream groupFile;
	groupFileOpts.openFile(groupFile);
	groupFile << "readName\tgroup" << "\n";
	for(const auto & seq : distanceGraph.nodes_){
		groupFile << seq->name_ << "\t" << seq->group_ << "\n";
	}

  std::ofstream outFile;
  openTextFile(outFile, setUp.pars_.directoryName_ + "distanceMatrix", ".tab.txt", setUp.pars_.ioOptions_.out_);
  writeDistanceMatrix(outFile, distances, readVec::getNames(reads));

	return 0;
}



}  // namespace njhseq

