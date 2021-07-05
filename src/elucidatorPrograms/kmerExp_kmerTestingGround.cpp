/*
 * kmerExp_kmerTestingGround.cpp
 *
 *  Created on: May 24, 2020
 *      Author: nick
 */




#include "kmerExp.hpp"
#include "elucidator/utils/KmerUtils.hpp"
#include "elucidator/objects/seqObjects/seqKmers.h"

#include <njhseq/concurrency/AllByAllPairFactory.hpp>
#include <njhseq/concurrency/PairwisePairFactory.hpp>

#include <njhseq/GenomeUtils/GenomeMapping/MultiGenomeMapper.hpp>
#include <njhseq/IO/SeqIO/SeqIO.hpp>



namespace njhseq {

std::string perfectKmerHtmlIndexFilePerRecStr = "<!DOCTYPE html>\n"
		"<html>\n"
		"<head>\n"
		"<meta charset=\"utf-8\">\n"
		"<title>Genome Kmer Comparison</title>\n"
		"<style>\n"
		"\n"
		"body {\n"
		"  font: 10px sans-serif;\n"
		"}\n"
		"\n"
		".axis path,\n"
		".axis line {\n"
		"  fill: none;\n"
		"  stroke: #000;\n"
		"  shape-rendering: crispEdges;\n"
		"}\n"
		"\n"
		"</style>\n"
		"\n"
		"<script src=\"http://cdnjs.cloudflare.com/ajax/libs/d3/3.4.13/d3.min.js\"></script>\n"
		"\n"
		"\n"
		"<body>\n"
		"\n"
		"<div id=\"figure\" style=\"margin-bottom: 50px;\"></div>\n"
		"\n"
		"\n"
		"\n"
		"<script>\n"
		"\n"
		"var margin = {top: 50, right: 10, bottom: 10, left: 200},\n"
		"    width = window.innerWidth *20 - margin.left - margin.right,\n"
		"    height = window.innerHeight - margin.top - margin.bottom;\n"
		"\n"
		"var y = d3.scale.ordinal()\n"
		"    .rangeRoundBands([0, height], .3);\n"
		"\n"
		"var x = d3.scale.linear()\n"
		"    .rangeRound([0, width]);\n"
		"\n"
		"var color = d3.scale.ordinal()\n"
		"    .range([\"#c7001e\", \"#f6a580\", \"#cccccc\", \"#92c6db\", \"#086fad\"]);\n"
		"\n"
		"var xAxis = d3.svg.axis()\n"
		"    .scale(x)\n"
		"    .orient(\"top\")\n"
		"    .ticks(100);\n"
		"\n"
		"var yAxis = d3.svg.axis()\n"
		"    .scale(y)\n"
		"    .orient(\"left\")\n"
		"\n"
		"var svg = d3.select(\"#figure\").append(\"svg\")\n"
		"    .attr(\"width\", width + margin.left + margin.right)\n"
		"    .attr(\"height\", height + margin.top + margin.bottom)\n"
		"    .attr(\"id\", \"d3-plot\")\n"
		"  .append(\"g\")\n"
		"    .attr(\"transform\", \"translate(\" + margin.left + \",\" + margin.top + \")\");\n"
		"\n"
		"  color.domain([\"Strongly disagree\", \"Disagree\", \"Neither agree nor disagree\", \"Agree\", \"Strongly agree\"]);\n"
		"\n"
		"  d3.json(\"INPUTJSONFNP\", function(error, data) {\n"
		"\n"
		"\n"
		"  x.domain([0, MAXPOSITION]).nice();\n"
		"\n"
		"  y.domain(data.map(function(d){ return d.name_}));\n"
		"\n"
		"  svg.append(\"g\")\n"
		"      .attr(\"class\", \"x axis\")\n"
		"      .call(xAxis);\n"
		"\n"
		"  svg.append(\"g\")\n"
		"      .attr(\"class\", \"y axis\")\n"
		"      .call(yAxis)\n"
		"  var genomes = svg.selectAll(\".genomes\")\n"
		"      .data(data)\n"
		"    .enter().append(\"g\")\n"
		"      .attr(\"class\", \"bar\")\n"
		"      .attr(\"id\", function(d) {return d.name_;})\n"
		"      .attr(\"transform\", function(d) {return \"translate(0,\" + y(d.name_) + \")\"; });\n"
		"\n"
		"  var bars = genomes.selectAll(\"rect\")\n"
		"      .data(function(d) {return d.kComps_; })\n"
		"    .enter().append(\"g\").attr(\"class\", \"subbar\");\n"
		"    \n"
		"  var tooltip = d3.select(\"body\")\n"
		"				.append(\"div\")\n"
		"				.style(\"position\", \"absolute\")\n"
		"				.style(\"visibility\", \"hidden\")\n"
		"				.style(\"background-color\", \"#88aaaa\")\n"
		"				.style(\"width\", \"300\")\n"
		"				.attr(\"id\", \"popHover\");\n"
		"  bars.append(\"rect\")\n"
		"      .attr(\"height\", y.rangeBand())\n"
		"      .attr(\"x\", function(d) { return x(d.start_); })\n"
		"      .attr(\"width\", function(d) { return x(d.size_); })\n"
		"      .attr(\"id\", function(d) { return \"ref\" + d.refStart_.toString(); })\n"
		"      .on(\"mouseover\", function(d,i,j){\n"
		"      		  var idName = \"#\"  + \"ref\" + d.refStart_.toString();\n"
		"      		  d3.selectAll(idName).style(\"fill\", \"#ffad33\");\n"
		"      		  var refStop = d.refStart_ + d.size_;\n"
		"      		  var seqStop = d.realStart_+ d.size_;\n"
		"	    	  tooltip.node().innerHTML = \"name: \"  + bars[j].parentNode.__data__.name_ + \"<br>refRange: \" + d.refStart_.toString() + \":\"  + refStop.toString() + \"<br>seqRange: \"  + d.realStart_.toString() + \":\" + seqStop.toString() + \"(size=\" +d.size_.toString()  +\")\" ;\n"
		"	    	  tooltip.style(\"visibility\", \"visible\");\n"
		"	      	 return;})\n"
		"	  .on(\"mousemove\", function(){\n"
		"	  	tooltip.style(\"top\", (d3.event.layerY-10)+\"px\").style(\"left\",(d3.event.layerX+10)+\"px\")\n"
		"	  	return ;})\n"
		"	  .on(\"mouseout\", function(d){\n"
		"	  	var idName = \"#\"  + \"ref\" + d.refStart_.toString();\n"
		"	  	d3.selectAll(idName).style(\"fill\", \"#000000\");\n"
		"	  	tooltip.style(\"visibility\", \"hidden\");\n"
		"	  return ;});\n"
		"\n"
		"\n"
		"  svg.append(\"g\")\n"
		"      .attr(\"class\", \"y axis\")\n"
		"  .append(\"line\")\n"
		"      .attr(\"x1\", x(0))\n"
		"      .attr(\"x2\", x(0))\n"
		"      .attr(\"y2\", height);\n"
		"\n"
		"\n"
		"  d3.selectAll(\".axis path\")\n"
		"      .style(\"fill\", \"none\")\n"
		"      .style(\"stroke\", \"#000\")\n"
		"      .style(\"shape-rendering\", \"crispEdges\");\n"
		"\n"
		"  d3.selectAll(\".axis line\")\n"
		"      .style(\"fill\", \"none\")\n"
		"      .style(\"stroke\", \"#000\")\n"
		"      .style(\"shape-rendering\", \"crispEdges\");\n"
		"});\n"
		"</script>\n"
		"\n"
		"</body>\n"
		"</html>\n";

//on 2020-11-24 kmerTestingGround changed to getKmerSharedBlocksBetweenGenomes
int kmerExpRunner::kmerTestingGround(const njh::progutils::CmdArgs & inputCommands){
	std::cout << "Nothing in testing ground currently, a previous version of getKmerSharedBlocksBetweenGenomes was until moved on 2020-11-24" << std::endl;
	return 0;
//	uint32_t klen = 35;
//	MultiGenomeMapper::inputParameters pars;
//	uint32_t allowableMissing = 1;
//
//
//
//	seqSetUp setUp(inputCommands);
//	setUp.processVerbose();
//	setUp.processDebug();
//
//	setUp.setOption(allowableMissing, "--allowableMissing", "allowable genomes missing");
//
//
//	setUp.setOption(pars.numThreads_, "--numThreads", "Number of CPUs to utilize");
//	setUp.setOption(pars.genomeDir_, "--genomeDir", "Name of the genome file fnp", true);
//	setUp.setOption(pars.primaryGenome_, "--primaryGenome", "The primary reference genome");
//  setUp.setOption(pars.selectedGenomes_, "--selectedGenomes", "Name of the other genomes in --genomeDir to be read in, leave blank to just do all fastas");
//
//
//	setUp.setOption(pars.gffDir_, "--gffDir", "A directory with a gff for the genomes in --genomeDir, should be named GENOME.gff (for GENOME.fasta)");
//  setUp.setOption(pars.gffIntersectPars_.extraAttributes_, "--gffExtraAttributes", "Extra attributes to add to genome that has an accompany gff");
//
//
//	setUp.setOption(klen, "--klen", "kmer Length", true);
//	setUp.processDirectoryOutputName(pars.primaryGenome_ + "_kmerTestingGround" + "_TODAY", true);
//	setUp.finishSetUp(std::cout);
//
//	setUp.startARunLog(setUp.pars_.directoryName_);
//  njh::stopWatch watch;
//
//  watch.setLapName("Set up genomes");
//
//	//set up genome mapper;
//	auto gMapper = std::make_unique<MultiGenomeMapper>(pars);
//	//set up selected genomes
//	if(!pars.selectedGenomes_.empty()){
//		gMapper->setSelectedGenomes(pars.selectedGenomes_);
//	}
//	//init is threaded
//	gMapper->init();
//
//	auto genomeNames = getVectorOfMapKeys(gMapper->genomes_);
//
//	std::unordered_map<std::string, std::unordered_map<std::string, std::shared_ptr<KmersSharedBlocks>>> uniqueKmersPerGenomePerRecord;
//	std::mutex uniqKsMut;
//
//	njh::concurrent::LockableQueue<std::string> genomeNamesQueue(genomeNames);
//	std::function<void()> indexGenomeForUniqueKmers = [&uniqueKmersPerGenomePerRecord,&uniqKsMut,&genomeNamesQueue,&gMapper,&klen](){
//		std::string genome = "";
//		while(genomeNamesQueue.getVal(genome)){
//		  std::unordered_map<std::string, std::shared_ptr<KmersSharedBlocks>> uniqueKmersPerRecord;
//		  {
//		  	std::unordered_map<std::string, uint32_t> genomeWideKmers;
//		  	auto genomeFnpInOpts = SeqIOOptions::genFastaIn(gMapper->genomes_.at(genome)->fnp_);
//		  	genomeFnpInOpts.includeWhiteSpaceInName_ = false;
//		  	{
//		  		//get genome wide kmers
//		    	SeqInput reader(genomeFnpInOpts);
//		    	reader.openIn();
//		    	seqInfo seq;
//		    	while (reader.readNextRead(seq)) {
//		    		if (len(seq) > klen + 1) {
//		    			for (uint32_t pos = 0; pos < len(seq) - klen + 1; ++pos) {
//		    				++genomeWideKmers[seq.seq_.substr(pos, klen)];
//		    			}
//		    		}
//		    	}
//		  	}
//		  	{
//		  		//get genome wide kmers
//		  		SeqInput reader(genomeFnpInOpts);
//		    	reader.openIn();
//		    	seqInfo seq;
//		    	while (reader.readNextRead(seq)) {
//		    		if (len(seq) > klen + 1) {
//		    			kmerInfo uniqueKmers;
//		    			uniqueKmers.kLen_ = klen;
//		    			for (uint32_t pos = 0; pos < len(seq) - klen + 1; ++pos) {
//		    				auto k = seq.seq_.substr(pos, klen);
//		    				if(1 == genomeWideKmers[k]){
//		    					uniqueKmers.kmers_[k] = kmer(k,pos);
//		    				}
//		    			}
//		    			uniqueKmersPerRecord[seq.name_] = std::make_shared<KmersSharedBlocks>(seq, uniqueKmers);
//		    		}
//		    	}
//		  	}
//		  }
//		  {
//		  	std::lock_guard<std::mutex> lock(uniqKsMut);
//		  	uniqueKmersPerGenomePerRecord[genome] = uniqueKmersPerRecord;
//		  }
//		}
//	};
//	watch.startNewLap("Get unique kmers for all genomes");
//	njh::concurrent::runVoidFunctionThreaded(indexGenomeForUniqueKmers, gMapper->pars_.numThreads_);
//
//
//	auto perGenomeUniCountsDir = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar{"uniqueKmersCounts"});
//	watch.startNewLap("Outputing info per genome");
//
//	for(const auto & genome : uniqueKmersPerGenomePerRecord){
//		OutputStream out(njh::files::make_path(perGenomeUniCountsDir, genome.first + "_uniqueKmerNumbers.tab.txt"));
//		out << "record1\tnumberOfUniqueKmers" << std::endl;
//		auto recordNames = getVectorOfMapKeys(genome.second);
//		njh::sort(recordNames);
//		for(const auto & name : recordNames){
//			out << name
//					<< "\t" << genome.second.at(name)->kInfo_->kmers_.size() << std::endl;
//		}
//
//		OutputStream outBed(njh::files::make_path(perGenomeUniCountsDir, genome.first + "_uniqueKmerNumbers.bed"));
//		for(const auto & name : recordNames){
//			std::vector<uint32_t> positions;
//			for(const auto & k : genome.second.at(name)->kInfo_->kmers_){
//				positions.emplace_back(k.second.positions_.front());
//			}
//			njh::sort(positions);
//			if(positions.size() > 0){
//				std::vector<Bed6RecordCore> regions;
//				regions.emplace_back(Bed6RecordCore(name, positions[0], positions[0] +1, "", 1, '+'));
//				if(positions.size() > 1){
//					for(uint32_t pos = 1; pos < positions.size(); ++pos){
//						auto gPos = positions[pos];
//						if(regions.back().chromEnd_ == gPos){
//							regions.back().chromEnd_ += 1;
//						}else{
//							regions.emplace_back(Bed6RecordCore(name, positions[pos], positions[pos] +1, "", 1, '+'));
//						}
//					}
//				}
//				njh::for_each(regions,[&klen](Bed6RecordCore & record){
//					record.chromEnd_ += klen - 1;
//					record.score_ = record.length();
//					record.name_ = record.genUIDFromCoords();
//				});
//				for(const auto & rec : regions){
//					outBed << rec.toDelimStr() << std::endl;
//				}
//			}
//		}
//	}
//	watch.startNewLap("Compare Against Ref");
//
//	auto refRecNames = getVectorOfMapKeys(uniqueKmersPerGenomePerRecord[pars.primaryGenome_]);
//	struct UniKmersCompRes {
//		UniKmersCompRes(const std::string &gRec, const std::string &rRec,
//				uint32_t genomeUniqKs, uint32_t refUniKs) :
//				genomeRec_(gRec), refRec_(rRec), genomeUniqKs_(genomeUniqKs), refUniKs_(
//						refUniKs) {
//
//		}
//		std::string genomeRec_;
//		std::string refRec_;
//		uint32_t genomeUniqKs_ { 0 };
//		uint32_t refUniKs_ { 0 };
//		uint32_t kmersShared_ { 0 };
//
//		double proportionOfGenomeUniqShared() const {
//			return static_cast<double>(kmersShared_) / genomeUniqKs_;
//		}
//		double indexScore() const {
//			return static_cast<double>(kmersShared_)
//					/ (genomeUniqKs_ + refUniKs_ - kmersShared_);
//		}
//
//	};
//	std::unordered_map<std::string, std::vector<UniKmersCompRes>> compsAgainstRef;
//	std::unordered_map<std::string, std::vector<std::shared_ptr<KmersSharedBlocks>>> blocksByRefRecord;
//
//	for(const auto & genome : uniqueKmersPerGenomePerRecord){
//		//if(genome.first != pars.primaryGenome_){
//		if(true){
//			OutputStream genomeCompFile(njh::files::make_path(setUp.pars_.directoryName_, genome.first + "_againstRefGenome.tab.txt"));
//			genomeCompFile << "record\trefRecord\ttotalUniqueKs\ttotalUniqueKsRef\tkShared\tindexScore\tproportionOfKsShared\tbestRecord" << std::endl;
//			std::mutex genomeCompMut;
//			auto genomeRecNames = getVectorOfMapKeys(genome.second);
//			AllByAllPairFactory pairFactory(genome.second.size(), uniqueKmersPerGenomePerRecord.at(pars.primaryGenome_).size());
//			std::function<void()> compToRef = [&refRecNames,&genomeRecNames,&pairFactory,&genome,&genomeCompMut,&compsAgainstRef,&pars,&uniqueKmersPerGenomePerRecord](){
//				AllByAllPairFactory::AllByAllPair pair;
//				while(pairFactory.setNextPair(pair)){
//					std::string refRec = refRecNames[pair.col_];
//					std::string genomeRec = genomeRecNames[pair.row_];
//					UniKmersCompRes compRes(genomeRec, refRec, genome.second.at(genomeRec)->kInfo_->kmers_.size(), uniqueKmersPerGenomePerRecord[pars.primaryGenome_][refRec]->kInfo_->kmers_.size());
//					for(const auto & k : genome.second.at(genomeRec)->kInfo_->kmers_){
//						if(njh::in(k.first, uniqueKmersPerGenomePerRecord[pars.primaryGenome_][refRec]->kInfo_->kmers_)){
//							++compRes.kmersShared_;
//						}
//					}
//					{
//						std::lock_guard<std::mutex> lock(genomeCompMut);
//						compsAgainstRef[genome.first].emplace_back(compRes);
//					}
//				}
//			};
//			njh::concurrent::runVoidFunctionThreaded(compToRef, pars.numThreads_);
//			njh::sort(compsAgainstRef[genome.first], [](const UniKmersCompRes & comp1, const UniKmersCompRes & comp2){
//				if(comp1.genomeRec_ == comp2.genomeRec_){
//					return comp1.refRec_ < comp2.refRec_;
//				}
//				return comp1.genomeRec_ < comp2.genomeRec_;
//			});
//			std::unordered_map<std::string, std::pair<std::string, double>> highestMatchingRefRec;
//			for(const auto & comp : compsAgainstRef[genome.first]){
//				if(comp.indexScore() > highestMatchingRefRec[comp.genomeRec_].second){
//					highestMatchingRefRec[comp.genomeRec_] = std::make_pair(comp.refRec_, comp.indexScore());
//				}
//			}
//			if(genome.first != pars.primaryGenome_){
//				for(const auto & comp : compsAgainstRef[genome.first]){
//					genomeCompFile << comp.genomeRec_
//							<< "\t" << comp.refRec_
//							<< "\t" << comp.genomeUniqKs_
//							<< "\t" << comp.refUniKs_
//							<< "\t" << comp.kmersShared_
//							<< "\t" << comp.indexScore()
//							<< "\t" << comp.proportionOfGenomeUniqShared()
//							<< "\t" << highestMatchingRefRec[comp.genomeRec_].first<< std::endl;
//				}
//			}
//			for(const auto & high : highestMatchingRefRec){
//				blocksByRefRecord[high.second.first].emplace_back(uniqueKmersPerGenomePerRecord[genome.first][high.first]);
//			}
//		}
//	}
//	{
//		watch.startNewLap("Search All");
//
//		njh::concurrent::LockableQueue<std::string> refNameQueue(refRecNames);
//		std::function<void()> profileRecord = [&refNameQueue,&blocksByRefRecord,&gMapper,&uniqueKmersPerGenomePerRecord,&klen](){
//			std::string recName = "";
//			TwoBit::TwoBitFile refReader(gMapper->genomes_.at(gMapper->pars_.primaryGenome_)->fnpTwoBit_);
//			while(refNameQueue.getVal(recName)){
//				std::string refSeq;
//				refReader[recName]->getSequence(refSeq);
//				for (const auto pos : iter::range(refSeq.size() - klen + 1)) {
//					auto sub = refSeq.substr(pos, klen);
//					for (auto & seq : blocksByRefRecord.at(recName)) {
//						auto search = seq->kInfo_->kmers_.find(sub);
//						//if both sequences contain the kmer and the kmer is unique in both sequences, add them
//						if (search != seq->kInfo_->kmers_.end()
//								&& search->second.positions_.size() == 1
//								&& njh::in(search->first, uniqueKmersPerGenomePerRecord.at(gMapper->pars_.primaryGenome_).at(recName)->kInfo_->kmers_)) {
//							seq->addComp(uniqueKmersPerGenomePerRecord.at(gMapper->pars_.primaryGenome_).at(recName)->kInfo_->kmers_[search->first].positions_.front(), search->second.positions_.front());
//						}
//					}
//				}
//				for (auto & seq : blocksByRefRecord.at(recName)) {
//					seq->finish();
//				}
//			}
//		};
//		njh::concurrent::runVoidFunctionThreaded(profileRecord, gMapper->pars_.numThreads_);
//	}
//
//
//	{
//		watch.startNewLap("Count Shared Kmers");
//
//		njh::concurrent::LockableQueue<std::string> refNameQueue(refRecNames);
//
//
//		auto sharedInfoDir = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar{"sharedInfo"});
//
//		std::function<void()> profileRecord = [&refNameQueue,&blocksByRefRecord,&gMapper,&klen,
//																					 &allowableMissing,&sharedInfoDir](){
//			std::string recName = "";
//			TwoBit::TwoBitFile refReader(gMapper->genomes_.at(gMapper->pars_.primaryGenome_)->fnpTwoBit_);
//			while(refNameQueue.getVal(recName)){
//
//
//				OutputStream perfectCompResultsJsonFile(njh::files::make_path(sharedInfoDir, recName + "_" + "perfect.json"));
//				OutputStream perfectCompResultsHtmlFile(njh::files::make_path(sharedInfoDir, recName + "_" + "perfect.html"));
//				OutputStream perfectCompResultsFile(njh::files::make_path(sharedInfoDir, recName + "_" + "perfect.fasta"));
//				OutputStream nearPerfectCompResultsFile(njh::files::make_path(sharedInfoDir, recName + "_" + "nearPerfect.fasta"));
//				OutputStream compResultsJsonFile(njh::files::make_path(sharedInfoDir, recName + "_" + "fractions.json"));
//				OutputStream compResultsFile(njh::files::make_path(sharedInfoDir, recName + "_" + "fractions.tab.txt"));
//
//
//
//
//				std::string refSeq;
//				refReader[recName]->getSequence(refSeq);
//				;
//				perfectCompResultsHtmlFile << njh::replaceString(njh::replaceString(perfectKmerHtmlIndexFilePerRecStr, "MAXPOSITION", estd::to_string(refSeq.size()))
//				, "INPUTJSONFNP", recName + "_" + "perfect.json");
//
//
//				//count up the shared kmers at the reference positions
//				std::vector<uint32_t> counts(refSeq.size(), 0);
//				for(const auto & g : blocksByRefRecord.at(recName)){
//					for(const auto & k : g->kComps_){
//						for(const auto pos : iter::range(k.second.refStart_, k.second.size_ + k.second.refStart_)){
//							++counts[pos];
//						}
//					}
//				}
//
//
//				std::unordered_map<std::string, std::shared_ptr<KmersSharedBlocks>> perfectSeqs;
//				std::unordered_map<std::string, std::shared_ptr<KmersSharedBlocks>> nearPerfectSeqs;
//
//				for(const auto & seq : blocksByRefRecord.at(recName)){
//					perfectSeqs[seq->seqBase_->name_] = std::make_shared<KmersSharedBlocks>(seqInfo(seq->seqBase_->name_), kmerInfo());
//					nearPerfectSeqs[seq->seqBase_->name_] = std::make_shared<KmersSharedBlocks>(seqInfo(seq->seqBase_->name_), kmerInfo());
//				}
//				for(const auto & seq : blocksByRefRecord.at(recName)){
//					for(const auto & k : seq->kComps_){
//						for(const auto pos : iter::range(k.second.size_)){
//							if(counts[k.second.refStart_ + pos] == blocksByRefRecord.at(recName).size()){
//								perfectSeqs.at(seq->seqBase_->name_)->addComp(k.second.refStart_ + pos,k.second.start_ + pos);
//							}
//							if(counts[k.second.refStart_ + pos] >= blocksByRefRecord.at(recName).size() - allowableMissing){
//								nearPerfectSeqs.at(seq->seqBase_->name_)->addComp(k.second.refStart_ + pos,k.second.start_ + pos);
//							}
//						}
//					}
//				}
//				for (auto & seq : perfectSeqs) {
//					seq.second->finish();
//				}
//				for (auto & seq : nearPerfectSeqs) {
//					seq.second->finish();
//				}
//
//
//				{
//					//per seq
//					for(const auto & seq : perfectSeqs){
//						OutputStream perfectCompBedPerGenomeFile(
//								njh::files::make_path(sharedInfoDir,
//										njh::pasteAsStr(recName, "_", seq.first, "_","perfect.bed")));
//						for (const auto &g : seq.second->kComps_){
//							auto seqStop = g.second.size_ + g.second.start_ + klen - 1;
//							perfectCompBedPerGenomeFile << seq.first
//									<< "\t" << g.second.start_
//									<< "\t" << seqStop
//									<< "\t" << njh::pasteAsStr(recName, "-", g.second.start_, "-", seqStop )
//									<< "\t" << g.second.size_ + klen - 1
//									<< "\t" << "+" << std::endl;
//						}
//					}
//				}
//
//				Json::Value perfectResults;
//				uint64_t maxMinStartAll = 0;
//				for (const auto &g : perfectSeqs) {
//					uint64_t minStart = std::numeric_limits<uint64_t>::max();
//					for (const auto &k : g.second->kComps_) {
//						if (k.second.start_ < minStart) {
//							minStart = k.second.start_;
//						}
//					}
//					if(minStart > maxMinStartAll){
//						maxMinStartAll = minStart;
//					}
//				}
//
//
//				for(const auto & g : perfectSeqs){
//					Json::Value currentGen;
//					currentGen["name_"] = g.first;
//					currentGen["kLen_"] = klen;
//					Json::Value & kComps = currentGen["kComps_"];
//					uint64_t minStart = std::numeric_limits<uint64_t>::max();
//					for(const auto & k : g.second->kComps_ ){
//						if(k.second.start_ < minStart){
//							minStart = k.second.start_;
//						}
//					}
//					for(const auto & k : g.second->kComps_ ){
//						auto outJsonForBlock = njh::json::toJson(k.second);
//						outJsonForBlock["realStart_"] = outJsonForBlock["start_"].asUInt64();
//						outJsonForBlock["start_"] = outJsonForBlock["start_"].asUInt64() - minStart + maxMinStartAll;
//						outJsonForBlock["size_"] = outJsonForBlock["size_"].asUInt64() + klen - 1;
//						kComps.append(outJsonForBlock);
//					}
//					perfectResults.append(currentGen);
//				}
//				perfectCompResultsJsonFile << njh::json::writeAsOneLine(perfectResults) << std::endl;
//				{
//					auto posKeys = getVectorOfMapKeys(perfectSeqs.begin()->second->kComps_) ;
//					njh::sort(posKeys);
//					for (const auto & posKey : posKeys) {
//						const auto & k = perfectSeqs.begin()->second->kComps_[posKey];
//						perfectCompResultsFile << ">" << k.refStart_ << std::endl;
//						perfectCompResultsFile
//								<< refSeq.substr(k.refStart_, k.size_ + klen - 1)
//								<< std::endl;
//					}
//				}
//
//				{
//					auto posKeys = getVectorOfMapKeys(nearPerfectSeqs.begin()->second->kComps_) ;
//					njh::sort(posKeys);
//					for (const auto & posKey : posKeys) {
//						const auto & k = nearPerfectSeqs.begin()->second->kComps_[posKey];
//						nearPerfectCompResultsFile << ">" << k.refStart_ << std::endl;
//						nearPerfectCompResultsFile
//								<< refSeq.substr(k.refStart_, k.size_ + klen - 1)
//								<< std::endl;
//					}
//				}
//				//fractional results, json
//				Json::Value fracResults;
//				fracResults["name_"] = recName;
//				fracResults["kLen_"] = klen;
//				Json::Value fracs;
//				for(const auto & c : counts){
//					fracs.append(static_cast<double>(c)/perfectSeqs.size());
//				}
//				fracResults["fracs"] = fracs;
//				compResultsJsonFile << njh::json::writeAsOneLine(fracResults) << std::endl;
//
//				compResultsFile << "refPos\trefSeq\tfractionShared\tsharing\tnotSharing\ttotalInput" << std::endl;
//				for(auto pos : iter::range(counts.size() - klen + 1)){
//					compResultsFile << pos
//							<< "\t" << refSeq.substr(pos,klen)
//							<< "\t" << static_cast<double>(counts[pos])/perfectSeqs.size()
//							<< "\t" << counts[pos]
//							<< "\t" << perfectSeqs.size() - counts[pos]
//							<< "\t" << perfectSeqs.size() << std::endl;
//				}
//			}
//		};
//		njh::concurrent::runVoidFunctionThreaded(profileRecord, gMapper->pars_.numThreads_);
//	}
//
//	/*
//	 * 	watch.startNewLap("Count All Shared Kmers");
//	//count up the shared kmers at the reference positions
//	std::vector<uint32_t> counts(len(setUp.pars_.seqObj_), 0);
//	for(const auto & g : seqs){
//		for(const auto & k : g.second.kComps_){
//			for(const auto & pos : iter::range(k.second.refStart_, k.second.size_ + k.second.refStart_)){
//				++counts[pos];
//			}
//		}
//	}
//	watch.startNewLap("Determine Shared Kmers Between All");
//	for(const auto & seq : seqs){
//		for(const auto & k : seq.second.kComps_){
//			for(const auto & pos : iter::range(k.second.size_)){
//				if(counts[k.second.refStart_ + pos] == seqs.size()){
//					perfectSeqs.at(seq.first).addComp(k.second.refStart_ + pos,k.second.start_ + pos);
//				}
//				if(counts[k.second.refStart_ + pos] >= seqs.size() - allowableMissing){
//					nearPerfectSeqs.at(seq.first).addComp(k.second.refStart_ + pos,k.second.start_ + pos);
//				}
//			}
//		}
//	}
//
//	for (auto & seq : perfectSeqs) {
//		seq.second.finish();
//	}
//
//	for (auto & seq : nearPerfectSeqs) {
//		seq.second.finish();
//	}
//
//	watch.startNewLap("Outputing");
//	Json::Value perfectResults;
//	uint64_t maxMinStartAll = 0;
//	for (const auto &g : perfectSeqs.begin()->second.kComps_){
//		auto refStop = g.second.size_ + g.second.refStart_ + kLen - 1;
//		perfectCompBedFile << setUp.pars_.seqObj_.seqBase_.name_
//				<< "\t" << g.second.refStart_
//				<< "\t" << refStop
//				<< "\t" << njh::pasteAsStr(setUp.pars_.seqObj_.seqBase_.name_, "-", g.second.refStart_, "-", refStop )
//				<< "\t" << g.second.size_
//				<< "\t" << "+" << std::endl;
//	}
//	for (const auto &g : perfectSeqs) {
//		uint64_t minStart = std::numeric_limits<uint64_t>::max();
//		for (const auto &k : g.second.kComps_) {
//			if (k.second.start_ < minStart) {
//				minStart = k.second.start_;
//			}
//		}
//		if(minStart > maxMinStartAll){
//			maxMinStartAll = minStart;
//		}
//	}
//	for(const auto & g : perfectSeqs){
//		Json::Value currentGen;
//		currentGen["name_"] = g.first;
//		currentGen["kLen_"] = kLen;
//		Json::Value & kComps = currentGen["kComps_"];
//		uint64_t minStart = std::numeric_limits<uint64_t>::max();
//		for(const auto & k : g.second.kComps_ ){
//			if(k.second.start_ < minStart){
//				minStart = k.second.start_;
//			}
//		}
//		for(const auto & k : g.second.kComps_ ){
//			auto outJsonForBlock = njh::json::toJson(k.second);
//			outJsonForBlock["realStart_"] = outJsonForBlock["start_"].asUInt64();
//			outJsonForBlock["start_"] = outJsonForBlock["start_"].asUInt64() - minStart + maxMinStartAll;
//			outJsonForBlock["size_"] = outJsonForBlock["size_"].asUInt64() + kLen - 1;
//			kComps.append(outJsonForBlock);
//		}
//		perfectResults.append(currentGen);
//	}
//	perfectCompResultsJsonFile << njh::json::writeAsOneLine(perfectResults) << std::endl;
//	{
//		auto posKeys = getVectorOfMapKeys(perfectSeqs.begin()->second.kComps_) ;
//		njh::sort(posKeys);
//		for (const auto & posKey : posKeys) {
//			const auto & k = perfectSeqs.begin()->second.kComps_[posKey];
//			perfectCompResultsFile << ">" << k.refStart_ << std::endl;
//			perfectCompResultsFile
//					<< setUp.pars_.seqObj_.seqBase_.seq_.substr(k.refStart_, k.size_ + kLen - 1)
//					<< std::endl;
//		}
//	}
//
//	{
//		auto posKeys = getVectorOfMapKeys(nearPerfectSeqs.begin()->second.kComps_) ;
//		njh::sort(posKeys);
//		for (const auto & posKey : posKeys) {
//			const auto & k = nearPerfectSeqs.begin()->second.kComps_[posKey];
//			nearPerfectCompResultsFile << ">" << k.refStart_ << std::endl;
//			nearPerfectCompResultsFile
//					<< setUp.pars_.seqObj_.seqBase_.seq_.substr(k.refStart_, k.size_ + kLen - 1)
//					<< std::endl;
//		}
//	}
//	//fractional results, json
//	Json::Value fracResults;
//	fracResults["name_"] = setUp.pars_.seqObj_.seqBase_.name_;
//	fracResults["kLen_"] = kLen;
//	Json::Value fracs;
//	for(const auto & c : counts){
//		fracs.append(static_cast<double>(c)/seqs.size());
//	}
//	fracResults["fracs"] = fracs;
//	compResultsJsonFile << njh::json::writeAsOneLine(fracResults) << std::endl;
//
//	compResultsFile << "refPos\trefSeq\tfractionShared\tsharing\tnotSharing\ttotalInput" << std::endl;
//	for(auto pos : iter::range(counts.size() - kLen + 1)){
//		compResultsFile << pos
//				<< "\t" << setUp.pars_.seqObj_.seqBase_.seq_.substr(pos,kLen)
//				<< "\t" << static_cast<double>(counts[pos])/seqs.size()
//				<< "\t" << counts[pos]
//				<< "\t" << seqs.size() - counts[pos]
//				<< "\t" << seqs.size() << std::endl;
//	}
//	 */
//
//  if(setUp.pars_.verbose_){
//  	watch.logLapTimes(std::cout, true, 6, true);
//  	std::cout << "totalTime: " << watch.totalTimeFormatted(6) << std::endl;
//  }else{
//  	OutputStream watchLogTime(njh::files::make_path(setUp.pars_.directoryName_, "time.txt"));
//  	watch.logLapTimes(watchLogTime, true, 6, true);
//  	watchLogTime << "totalTime: " << watch.totalTimeFormatted(6) << std::endl;
//  }
//
//	return 0;
}

int kmerExpRunner::getKmerSharedBlocksBetweenGenomes(const njh::progutils::CmdArgs & inputCommands){
	uint32_t klen = 35;
	MultiGenomeMapper::inputParameters pars;
	uint32_t allowableMissing = 1;



	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();

	setUp.setOption(allowableMissing, "--allowableMissing", "allowable genomes missing");


	setUp.setOption(pars.numThreads_, "--numThreads", "Number of CPUs to utilize");
	setUp.setOption(pars.genomeDir_, "--genomeDir", "Name of the genome file fnp", true);
	setUp.setOption(pars.primaryGenome_, "--primaryGenome", "The primary reference genome", true);
  setUp.setOption(pars.selectedGenomes_, "--selectedGenomes", "Name of the other genomes in --genomeDir to be read in, leave blank to just do all fastas");


	setUp.setOption(pars.gffDir_, "--gffDir", "A directory with a gff for the genomes in --genomeDir, should be named GENOME.gff (for GENOME.fasta)");
  setUp.setOption(pars.gffIntersectPars_.extraAttributes_, "--gffExtraAttributes", "Extra attributes to add to genome that has an accompany gff");


	setUp.setOption(klen, "--klen", "kmer Length", true);
	setUp.processDirectoryOutputName(pars.primaryGenome_ + "_kmerTestingGround" + "_TODAY", true);
	setUp.finishSetUp(std::cout);

	setUp.startARunLog(setUp.pars_.directoryName_);
  njh::stopWatch watch;

  watch.setLapName("Set up genomes");

	//set up genome mapper;
	auto gMapper = std::make_unique<MultiGenomeMapper>(pars);
	//set up selected genomes
	if(!pars.selectedGenomes_.empty()){
		gMapper->setSelectedGenomes(pars.selectedGenomes_);
	}
	//init is threaded
	gMapper->init();

	auto genomeNames = getVectorOfMapKeys(gMapper->genomes_);
	std::unordered_map<std::string, std::unordered_map<std::string, std::shared_ptr<KmersSharedBlocks>>> uniqueKmersPerGenomePerRecord;
	std::mutex uniqKsMut;

	njh::concurrent::LockableQueue<std::string> genomeNamesQueue(genomeNames);
	std::function<void()> indexGenomeForUniqueKmers = [&uniqueKmersPerGenomePerRecord,&uniqKsMut,&genomeNamesQueue,&gMapper,&klen](){
		std::string genome = "";
		while(genomeNamesQueue.getVal(genome)){
		  std::unordered_map<std::string, std::shared_ptr<KmersSharedBlocks>> uniqueKmersPerRecord;
		  {
		  	std::unordered_map<std::string, uint32_t> genomeWideKmers;
		  	auto genomeFnpInOpts = SeqIOOptions::genFastaIn(gMapper->genomes_.at(genome)->fnp_);
		  	genomeFnpInOpts.includeWhiteSpaceInName_ = false;
		  	{
		  		//get genome wide kmers
		    	SeqInput reader(genomeFnpInOpts);
		    	reader.openIn();
		    	seqInfo seq;
		    	while (reader.readNextRead(seq)) {
		    		if (len(seq) > klen + 1) {
		    			for (uint32_t pos = 0; pos < len(seq) - klen + 1; ++pos) {
		    				++genomeWideKmers[seq.seq_.substr(pos, klen)];
		    			}
		    		}
		    	}
		  	}
		  	{
		  		//get genome wide kmers
		  		SeqInput reader(genomeFnpInOpts);
		    	reader.openIn();
		    	seqInfo seq;
		    	while (reader.readNextRead(seq)) {
		    		if (len(seq) > klen + 1) {
		    			kmerInfo uniqueKmers;
		    			uniqueKmers.kLen_ = klen;
		    			for (uint32_t pos = 0; pos < len(seq) - klen + 1; ++pos) {
		    				auto k = seq.seq_.substr(pos, klen);
		    				if(1 == genomeWideKmers[k]){
		    					uniqueKmers.kmers_[k] = kmer(k,pos);
		    				}
		    			}
		    			uniqueKmersPerRecord[seq.name_] = std::make_shared<KmersSharedBlocks>(seq, uniqueKmers);
		    		}
		    	}
		  	}
		  }
		  {
		  	std::lock_guard<std::mutex> lock(uniqKsMut);
		  	uniqueKmersPerGenomePerRecord[genome] = uniqueKmersPerRecord;
		  }
		}
	};
	watch.startNewLap("Get unique kmers for all genomes");
	njh::concurrent::runVoidFunctionThreaded(indexGenomeForUniqueKmers, gMapper->pars_.numThreads_);

	auto perGenomeUniCountsDir = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar{"uniqueKmersCounts"});
	watch.startNewLap("Outputing info per genome");
	for(const auto & genome : uniqueKmersPerGenomePerRecord){
		OutputStream out(njh::files::make_path(perGenomeUniCountsDir, genome.first + "_uniqueKmerNumbers.tab.txt"));
		out << "record1\tnumberOfUniqueKmers" << std::endl;
		auto recordNames = getVectorOfMapKeys(genome.second);
		njh::sort(recordNames);
		for(const auto & name : recordNames){
			out << name
					<< "\t" << genome.second.at(name)->kInfo_->kmers_.size() << std::endl;
		}

		OutputStream outBed(njh::files::make_path(perGenomeUniCountsDir, genome.first + "_uniqueKmerNumbers.bed"));
		for(const auto & name : recordNames){
			std::vector<uint32_t> positions;
			for(const auto & k : genome.second.at(name)->kInfo_->kmers_){
				positions.emplace_back(k.second.positions_.front());
			}
			njh::sort(positions);
			if(positions.size() > 0){
				std::vector<Bed6RecordCore> regions;
				regions.emplace_back(Bed6RecordCore(name, positions[0], positions[0] +1, "", 1, '+'));
				if(positions.size() > 1){
					for(uint32_t pos = 1; pos < positions.size(); ++pos){
						auto gPos = positions[pos];
						if(regions.back().chromEnd_ == gPos){
							regions.back().chromEnd_ += 1;
						}else{
							regions.emplace_back(Bed6RecordCore(name, positions[pos], positions[pos] +1, "", 1, '+'));
						}
					}
				}
				njh::for_each(regions,[&klen](Bed6RecordCore & record){
					record.chromEnd_ += klen - 1;
					record.score_ = record.length();
					record.name_ = record.genUIDFromCoords();
				});
				for(const auto & rec : regions){
					outBed << rec.toDelimStr() << std::endl;
				}
			}
		}
	}
	watch.startNewLap("Compare Against Ref");
	auto refRecNames = getVectorOfMapKeys(uniqueKmersPerGenomePerRecord[pars.primaryGenome_]);
	struct UniKmersCompRes {
		UniKmersCompRes(const std::string &gRec, const std::string &rRec,
				uint32_t genomeUniqKs, uint32_t refUniKs) :
				genomeRec_(gRec), refRec_(rRec), genomeUniqKs_(genomeUniqKs), refUniKs_(
						refUniKs) {

		}
		std::string genomeRec_;
		std::string refRec_;
		uint32_t genomeUniqKs_ { 0 };
		uint32_t refUniKs_ { 0 };
		uint32_t kmersShared_ { 0 };

		double proportionOfGenomeUniqShared() const {
			return static_cast<double>(kmersShared_) / genomeUniqKs_;
		}
		double indexScore() const {
			return static_cast<double>(kmersShared_)
					/ (genomeUniqKs_ + refUniKs_ - kmersShared_);
		}

	};
	std::unordered_map<std::string, std::vector<UniKmersCompRes>> compsAgainstRef;
	std::unordered_map<std::string, std::vector<std::shared_ptr<KmersSharedBlocks>>> blocksByRefRecord;
	for(const auto & genome : uniqueKmersPerGenomePerRecord){
		//if(genome.first != pars.primaryGenome_){
		if(true){
			OutputStream genomeCompFile(njh::files::make_path(setUp.pars_.directoryName_, genome.first + "_againstRefGenome.tab.txt"));
			genomeCompFile << "record\trefRecord\ttotalUniqueKs\ttotalUniqueKsRef\tkShared\tindexScore\tproportionOfKsShared\tbestRecord" << std::endl;
			std::mutex genomeCompMut;
			auto genomeRecNames = getVectorOfMapKeys(genome.second);
			AllByAllPairFactory pairFactory(genome.second.size(), uniqueKmersPerGenomePerRecord.at(pars.primaryGenome_).size());
			std::function<void()> compToRef = [&refRecNames,&genomeRecNames,&pairFactory,&genome,&genomeCompMut,&compsAgainstRef,&pars,&uniqueKmersPerGenomePerRecord](){
				AllByAllPairFactory::AllByAllPair pair;
				while(pairFactory.setNextPair(pair)){
					std::string refRec = refRecNames[pair.col_];
					std::string genomeRec = genomeRecNames[pair.row_];
					UniKmersCompRes compRes(genomeRec, refRec, genome.second.at(genomeRec)->kInfo_->kmers_.size(), uniqueKmersPerGenomePerRecord[pars.primaryGenome_][refRec]->kInfo_->kmers_.size());
					for(const auto & k : genome.second.at(genomeRec)->kInfo_->kmers_){
						if(njh::in(k.first, uniqueKmersPerGenomePerRecord[pars.primaryGenome_][refRec]->kInfo_->kmers_)){
							++compRes.kmersShared_;
						}
					}
					{
						std::lock_guard<std::mutex> lock(genomeCompMut);
						compsAgainstRef[genome.first].emplace_back(compRes);
					}
				}
			};
			njh::concurrent::runVoidFunctionThreaded(compToRef, pars.numThreads_);
			njh::sort(compsAgainstRef[genome.first], [](const UniKmersCompRes & comp1, const UniKmersCompRes & comp2){
				if(comp1.genomeRec_ == comp2.genomeRec_){
					return comp1.refRec_ < comp2.refRec_;
				}
				return comp1.genomeRec_ < comp2.genomeRec_;
			});
			std::unordered_map<std::string, std::pair<std::string, double>> highestMatchingRefRec;
			for(const auto & comp : compsAgainstRef[genome.first]){
				if(comp.indexScore() > highestMatchingRefRec[comp.genomeRec_].second){
					highestMatchingRefRec[comp.genomeRec_] = std::make_pair(comp.refRec_, comp.indexScore());
				}
			}
			if(genome.first != pars.primaryGenome_){
				for(const auto & comp : compsAgainstRef[genome.first]){
					genomeCompFile << comp.genomeRec_
							<< "\t" << comp.refRec_
							<< "\t" << comp.genomeUniqKs_
							<< "\t" << comp.refUniKs_
							<< "\t" << comp.kmersShared_
							<< "\t" << comp.indexScore()
							<< "\t" << comp.proportionOfGenomeUniqShared()
							<< "\t" << highestMatchingRefRec[comp.genomeRec_].first<< std::endl;
				}
			}
			for(const auto & high : highestMatchingRefRec){
				blocksByRefRecord[high.second.first].emplace_back(uniqueKmersPerGenomePerRecord[genome.first][high.first]);
			}
		}
	}
	{
		watch.startNewLap("Search All");

		njh::concurrent::LockableQueue<std::string> refNameQueue(refRecNames);
		std::function<void()> profileRecord = [&refNameQueue,&blocksByRefRecord,&gMapper,&uniqueKmersPerGenomePerRecord,&klen](){
			std::string recName = "";
			TwoBit::TwoBitFile refReader(gMapper->genomes_.at(gMapper->pars_.primaryGenome_)->fnpTwoBit_);
			while(refNameQueue.getVal(recName)){
				std::string refSeq;
				refReader[recName]->getSequence(refSeq);
				for (const auto pos : iter::range(refSeq.size() - klen + 1)) {
					auto sub = refSeq.substr(pos, klen);
					for (auto & seq : blocksByRefRecord.at(recName)) {
						auto search = seq->kInfo_->kmers_.find(sub);
						//if both sequences contain the kmer and the kmer is unique in both sequences, add them
						if (search != seq->kInfo_->kmers_.end()
								&& search->second.positions_.size() == 1
								&& njh::in(search->first, uniqueKmersPerGenomePerRecord.at(gMapper->pars_.primaryGenome_).at(recName)->kInfo_->kmers_)) {
							seq->addComp(uniqueKmersPerGenomePerRecord.at(gMapper->pars_.primaryGenome_).at(recName)->kInfo_->kmers_[search->first].positions_.front(), search->second.positions_.front());
						}
					}
				}
				for (auto & seq : blocksByRefRecord.at(recName)) {
					seq->finish();
				}
			}
		};
		njh::concurrent::runVoidFunctionThreaded(profileRecord, gMapper->pars_.numThreads_);
	}

	{
		watch.startNewLap("Count Shared Kmers");

		njh::concurrent::LockableQueue<std::string> refNameQueue(refRecNames);


		auto sharedInfoDir = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar{"sharedInfo"});

		std::function<void()> profileRecord = [&refNameQueue,&blocksByRefRecord,&gMapper,&klen,
																					 &allowableMissing,&sharedInfoDir](){
			std::string recName = "";
			TwoBit::TwoBitFile refReader(gMapper->genomes_.at(gMapper->pars_.primaryGenome_)->fnpTwoBit_);
			while(refNameQueue.getVal(recName)){


				OutputStream perfectCompResultsJsonFile(njh::files::make_path(sharedInfoDir, recName + "_" + "perfect.json"));
				OutputStream perfectCompResultsHtmlFile(njh::files::make_path(sharedInfoDir, recName + "_" + "perfect.html"));
				OutputStream perfectCompResultsFile(njh::files::make_path(sharedInfoDir, recName + "_" + "perfect.fasta"));
				OutputStream nearPerfectCompResultsFile(njh::files::make_path(sharedInfoDir, recName + "_" + "nearPerfect.fasta"));
				OutputStream compResultsJsonFile(njh::files::make_path(sharedInfoDir, recName + "_" + "fractions.json"));
				OutputStream compResultsFile(njh::files::make_path(sharedInfoDir, recName + "_" + "fractions.tab.txt"));




				std::string refSeq;
				refReader[recName]->getSequence(refSeq);
				;
				perfectCompResultsHtmlFile << njh::replaceString(njh::replaceString(perfectKmerHtmlIndexFilePerRecStr, "MAXPOSITION", estd::to_string(refSeq.size()))
				, "INPUTJSONFNP", recName + "_" + "perfect.json");


				//count up the shared kmers at the reference positions
				std::vector<uint32_t> counts(refSeq.size(), 0);
				for(const auto & g : blocksByRefRecord.at(recName)){
					for(const auto & k : g->kComps_){
						for(const auto pos : iter::range(k.second.refStart_, k.second.size_ + k.second.refStart_)){
							++counts[pos];
						}
					}
				}


				std::unordered_map<std::string, std::shared_ptr<KmersSharedBlocks>> perfectSeqs;
				std::unordered_map<std::string, std::shared_ptr<KmersSharedBlocks>> nearPerfectSeqs;

				for(const auto & seq : blocksByRefRecord.at(recName)){
					perfectSeqs[seq->seqBase_->name_] = std::make_shared<KmersSharedBlocks>(seqInfo(seq->seqBase_->name_), kmerInfo());
					nearPerfectSeqs[seq->seqBase_->name_] = std::make_shared<KmersSharedBlocks>(seqInfo(seq->seqBase_->name_), kmerInfo());
				}
				for(const auto & seq : blocksByRefRecord.at(recName)){
					for(const auto & k : seq->kComps_){
						for(const auto pos : iter::range(k.second.size_)){
							if(counts[k.second.refStart_ + pos] == blocksByRefRecord.at(recName).size()){
								perfectSeqs.at(seq->seqBase_->name_)->addComp(k.second.refStart_ + pos,k.second.start_ + pos);
							}
							if(counts[k.second.refStart_ + pos] >= blocksByRefRecord.at(recName).size() - allowableMissing){
								nearPerfectSeqs.at(seq->seqBase_->name_)->addComp(k.second.refStart_ + pos,k.second.start_ + pos);
							}
						}
					}
				}
				for (auto & seq : perfectSeqs) {
					seq.second->finish();
				}
				for (auto & seq : nearPerfectSeqs) {
					seq.second->finish();
				}


				{
					//per seq
					for(const auto & seq : perfectSeqs){
						OutputStream perfectCompBedPerGenomeFile(
								njh::files::make_path(sharedInfoDir,
										njh::pasteAsStr(recName, "_", seq.first, "_","perfect.bed")));
						for (const auto &g : seq.second->kComps_){
							auto seqStop = g.second.size_ + g.second.start_ + klen - 1;
							perfectCompBedPerGenomeFile << seq.first
									<< "\t" << g.second.start_
									<< "\t" << seqStop
									<< "\t" << njh::pasteAsStr(recName, "-", g.second.refStart_, "-", g.second.refStart_ + (g.second.size_ + klen - 1) )
									<< "\t" << g.second.size_ + klen - 1
									<< "\t" << "+" << std::endl;
						}
					}
				}

				Json::Value perfectResults;
				uint64_t maxMinStartAll = 0;
				for (const auto &g : perfectSeqs) {
					uint64_t minStart = std::numeric_limits<uint64_t>::max();
					for (const auto &k : g.second->kComps_) {
						if (k.second.start_ < minStart) {
							minStart = k.second.start_;
						}
					}
					if(minStart > maxMinStartAll){
						maxMinStartAll = minStart;
					}
				}


				for(const auto & g : perfectSeqs){
					Json::Value currentGen;
					currentGen["name_"] = g.first;
					currentGen["kLen_"] = klen;
					Json::Value & kComps = currentGen["kComps_"];
					uint64_t minStart = std::numeric_limits<uint64_t>::max();
					for(const auto & k : g.second->kComps_ ){
						if(k.second.start_ < minStart){
							minStart = k.second.start_;
						}
					}
					for(const auto & k : g.second->kComps_ ){
						auto outJsonForBlock = njh::json::toJson(k.second);
						outJsonForBlock["realStart_"] = outJsonForBlock["start_"].asUInt64();
						outJsonForBlock["start_"] = outJsonForBlock["start_"].asUInt64() - minStart + maxMinStartAll;
						outJsonForBlock["size_"] = outJsonForBlock["size_"].asUInt64() + klen - 1;
						kComps.append(outJsonForBlock);
					}
					perfectResults.append(currentGen);
				}
				perfectCompResultsJsonFile << njh::json::writeAsOneLine(perfectResults) << std::endl;
				{
					auto posKeys = getVectorOfMapKeys(perfectSeqs.begin()->second->kComps_) ;
					njh::sort(posKeys);
					for (const auto & posKey : posKeys) {
						const auto & k = perfectSeqs.begin()->second->kComps_[posKey];
						perfectCompResultsFile << ">" << k.refStart_ << std::endl;
						perfectCompResultsFile
								<< refSeq.substr(k.refStart_, k.size_ + klen - 1)
								<< std::endl;
					}
				}

				{
					auto posKeys = getVectorOfMapKeys(nearPerfectSeqs.begin()->second->kComps_) ;
					njh::sort(posKeys);
					for (const auto & posKey : posKeys) {
						const auto & k = nearPerfectSeqs.begin()->second->kComps_[posKey];
						nearPerfectCompResultsFile << ">" << k.refStart_ << std::endl;
						nearPerfectCompResultsFile
								<< refSeq.substr(k.refStart_, k.size_ + klen - 1)
								<< std::endl;
					}
				}
				//fractional results, json
				Json::Value fracResults;
				fracResults["name_"] = recName;
				fracResults["kLen_"] = klen;
				Json::Value fracs;
				for(const auto & c : counts){
					fracs.append(static_cast<double>(c)/perfectSeqs.size());
				}
				fracResults["fracs"] = fracs;
				compResultsJsonFile << njh::json::writeAsOneLine(fracResults) << std::endl;

				compResultsFile << "refPos\trefSeq\tfractionShared\tsharing\tnotSharing\ttotalInput" << std::endl;
				for(auto pos : iter::range(counts.size() - klen + 1)){
					compResultsFile << pos
							<< "\t" << refSeq.substr(pos,klen)
							<< "\t" << static_cast<double>(counts[pos])/perfectSeqs.size()
							<< "\t" << counts[pos]
							<< "\t" << perfectSeqs.size() - counts[pos]
							<< "\t" << perfectSeqs.size() << std::endl;
				}
			}
		};
		njh::concurrent::runVoidFunctionThreaded(profileRecord, gMapper->pars_.numThreads_);
	}

	/*
	 * 	watch.startNewLap("Count All Shared Kmers");
	//count up the shared kmers at the reference positions
	std::vector<uint32_t> counts(len(setUp.pars_.seqObj_), 0);
	for(const auto & g : seqs){
		for(const auto & k : g.second.kComps_){
			for(const auto & pos : iter::range(k.second.refStart_, k.second.size_ + k.second.refStart_)){
				++counts[pos];
			}
		}
	}
	watch.startNewLap("Determine Shared Kmers Between All");
	for(const auto & seq : seqs){
		for(const auto & k : seq.second.kComps_){
			for(const auto & pos : iter::range(k.second.size_)){
				if(counts[k.second.refStart_ + pos] == seqs.size()){
					perfectSeqs.at(seq.first).addComp(k.second.refStart_ + pos,k.second.start_ + pos);
				}
				if(counts[k.second.refStart_ + pos] >= seqs.size() - allowableMissing){
					nearPerfectSeqs.at(seq.first).addComp(k.second.refStart_ + pos,k.second.start_ + pos);
				}
			}
		}
	}

	for (auto & seq : perfectSeqs) {
		seq.second.finish();
	}

	for (auto & seq : nearPerfectSeqs) {
		seq.second.finish();
	}

	watch.startNewLap("Outputing");
	Json::Value perfectResults;
	uint64_t maxMinStartAll = 0;
	for (const auto &g : perfectSeqs.begin()->second.kComps_){
		auto refStop = g.second.size_ + g.second.refStart_ + kLen - 1;
		perfectCompBedFile << setUp.pars_.seqObj_.seqBase_.name_
				<< "\t" << g.second.refStart_
				<< "\t" << refStop
				<< "\t" << njh::pasteAsStr(setUp.pars_.seqObj_.seqBase_.name_, "-", g.second.refStart_, "-", refStop )
				<< "\t" << g.second.size_
				<< "\t" << "+" << std::endl;
	}
	for (const auto &g : perfectSeqs) {
		uint64_t minStart = std::numeric_limits<uint64_t>::max();
		for (const auto &k : g.second.kComps_) {
			if (k.second.start_ < minStart) {
				minStart = k.second.start_;
			}
		}
		if(minStart > maxMinStartAll){
			maxMinStartAll = minStart;
		}
	}
	for(const auto & g : perfectSeqs){
		Json::Value currentGen;
		currentGen["name_"] = g.first;
		currentGen["kLen_"] = kLen;
		Json::Value & kComps = currentGen["kComps_"];
		uint64_t minStart = std::numeric_limits<uint64_t>::max();
		for(const auto & k : g.second.kComps_ ){
			if(k.second.start_ < minStart){
				minStart = k.second.start_;
			}
		}
		for(const auto & k : g.second.kComps_ ){
			auto outJsonForBlock = njh::json::toJson(k.second);
			outJsonForBlock["realStart_"] = outJsonForBlock["start_"].asUInt64();
			outJsonForBlock["start_"] = outJsonForBlock["start_"].asUInt64() - minStart + maxMinStartAll;
			outJsonForBlock["size_"] = outJsonForBlock["size_"].asUInt64() + kLen - 1;
			kComps.append(outJsonForBlock);
		}
		perfectResults.append(currentGen);
	}
	perfectCompResultsJsonFile << njh::json::writeAsOneLine(perfectResults) << std::endl;
	{
		auto posKeys = getVectorOfMapKeys(perfectSeqs.begin()->second.kComps_) ;
		njh::sort(posKeys);
		for (const auto & posKey : posKeys) {
			const auto & k = perfectSeqs.begin()->second.kComps_[posKey];
			perfectCompResultsFile << ">" << k.refStart_ << std::endl;
			perfectCompResultsFile
					<< setUp.pars_.seqObj_.seqBase_.seq_.substr(k.refStart_, k.size_ + kLen - 1)
					<< std::endl;
		}
	}

	{
		auto posKeys = getVectorOfMapKeys(nearPerfectSeqs.begin()->second.kComps_) ;
		njh::sort(posKeys);
		for (const auto & posKey : posKeys) {
			const auto & k = nearPerfectSeqs.begin()->second.kComps_[posKey];
			nearPerfectCompResultsFile << ">" << k.refStart_ << std::endl;
			nearPerfectCompResultsFile
					<< setUp.pars_.seqObj_.seqBase_.seq_.substr(k.refStart_, k.size_ + kLen - 1)
					<< std::endl;
		}
	}
	//fractional results, json
	Json::Value fracResults;
	fracResults["name_"] = setUp.pars_.seqObj_.seqBase_.name_;
	fracResults["kLen_"] = kLen;
	Json::Value fracs;
	for(const auto & c : counts){
		fracs.append(static_cast<double>(c)/seqs.size());
	}
	fracResults["fracs"] = fracs;
	compResultsJsonFile << njh::json::writeAsOneLine(fracResults) << std::endl;

	compResultsFile << "refPos\trefSeq\tfractionShared\tsharing\tnotSharing\ttotalInput" << std::endl;
	for(auto pos : iter::range(counts.size() - kLen + 1)){
		compResultsFile << pos
				<< "\t" << setUp.pars_.seqObj_.seqBase_.seq_.substr(pos,kLen)
				<< "\t" << static_cast<double>(counts[pos])/seqs.size()
				<< "\t" << counts[pos]
				<< "\t" << seqs.size() - counts[pos]
				<< "\t" << seqs.size() << std::endl;
	}
	 */

  if(setUp.pars_.verbose_){
  	watch.logLapTimes(std::cout, true, 6, true);
  	std::cout << "totalTime: " << watch.totalTimeFormatted(6) << std::endl;
  }else{
  	OutputStream watchLogTime(njh::files::make_path(setUp.pars_.directoryName_, "time.txt"));
  	watch.logLapTimes(watchLogTime, true, 6, true);
  	watchLogTime << "totalTime: " << watch.totalTimeFormatted(6) << std::endl;
  }

	return 0;
}




int kmerExpRunner::getWithinGenomeUniqueKmers(const njh::progutils::CmdArgs & inputCommands){
	uint32_t klen = 35;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processReadInNames(true);
	setUp.setOption(klen, "--klen", "kmer Length", true);
	setUp.processDirectoryOutputName(true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
  njh::stopWatch watch;
  OutputStream out(njh::files::make_path(setUp.pars_.directoryName_, "uniqueKmerNumbers.tab.txt"));
  std::unordered_map<std::string, std::shared_ptr<KmersSharedBlocks>> uniqueKmersPerRecord;
  {
  	std::unordered_map<std::string, uint32_t> genomeWideKmers;
  	watch.setLapName("Get genome wide kmers");
  	{
  		//get genome wide kmers
    	SeqInput reader(setUp.pars_.ioOptions_);
    	reader.openIn();
    	seqInfo seq;
    	while (reader.readNextRead(seq)) {
    		if (len(seq) > klen + 1) {
    			for (uint32_t pos = 0; pos < len(seq) - klen + 1; ++pos) {
    				++genomeWideKmers[seq.seq_.substr(pos, klen)];
    			}
    		}
    	}
  	}
  	watch.startNewLap("Get unique kmers");
  	{
  		//get genome wide kmers
    	SeqInput reader(setUp.pars_.ioOptions_);
    	reader.openIn();
    	seqInfo seq;
    	while (reader.readNextRead(seq)) {
    		if (len(seq) > klen + 1) {
    			kmerInfo uniqueKmers;
    			uniqueKmers.kLen_ = klen;
    			for (uint32_t pos = 0; pos < len(seq) - klen + 1; ++pos) {
    				auto k = seq.seq_.substr(pos, klen);
    				if(1 == genomeWideKmers[k]){
    					uniqueKmers.kmers_[k] = kmer(k,pos);
    				}
    			}
    			uniqueKmersPerRecord[seq.name_] = std::make_shared<KmersSharedBlocks>(seq, uniqueKmers);
    		}
    	}
  	}
  }


	out << "record1\tnumberOfUniqueKmers" << std::endl;
	auto recordNames = getVectorOfMapKeys(uniqueKmersPerRecord);
	njh::sort(recordNames);
	for(const auto & name : recordNames){
		out << name
				<< "\t" << uniqueKmersPerRecord[name]->kInfo_->kmers_.size() << std::endl;
	}

  OutputStream outBed(njh::files::make_path(setUp.pars_.directoryName_, "uniqueKmerNumbers.bed"));


	for(const auto & name : recordNames){
		std::vector<uint32_t> positions;
		for(const auto & k : uniqueKmersPerRecord[name]->kInfo_->kmers_){
			positions.emplace_back(k.second.positions_.front());
		}
		njh::sort(positions);
		if(positions.size() > 0){
			std::vector<Bed6RecordCore> regions;
			regions.emplace_back(Bed6RecordCore(name, positions[0], positions[0] +1, "", 1, '+'));
			if(positions.size() > 1){
				for(uint32_t pos = 1; pos < positions.size(); ++pos){
					auto gPos = positions[pos];
					if(regions.back().chromEnd_ == gPos){
						regions.back().chromEnd_ += 1;
					}else{
						regions.emplace_back(Bed6RecordCore(name, positions[pos], positions[pos] +1, "", 1, '+'));
					}
				}
			}
			njh::for_each(regions,[&klen](Bed6RecordCore & record){
				record.chromEnd_ += klen - 1;
				record.score_ = record.length();
				record.name_ = record.genUIDFromCoords();
			});
			for(const auto & rec : regions){
				outBed << rec.toDelimStr() << std::endl;
			}
		}
	}


  if(setUp.pars_.verbose_){
  	watch.logLapTimes(std::cout, true, 6, true);
  	std::cout << "totalTime: " << watch.totalTimeFormatted(6) << std::endl;
  }
	return 0;
}



int kmerExpRunner::allByAllComparisonOfUniqueKmers(const njh::progutils::CmdArgs & inputCommands){
	uint32_t klen = 35;
	uint32_t numThreads = 1;
	bool getReverseComp = false;
	OutOptions outOpts(bfs::path(""), ".tab.txt");
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processReadInNames(true);
	setUp.processWritingOptions(outOpts);
	setUp.processRefFilename(true);
	setUp.setOption(numThreads, "--numThreads", "Number of threads to use", njh::progutils::ProgramSetUp::CheckCase::GTEQ1);
	setUp.setOption(klen, "--klen", "kmer Length", true);
	setUp.setOption(getReverseComp, "--getReverseComp", "Compare the reverse complement as well");

	setUp.finishSetUp(std::cout);

  njh::stopWatch watch;

  OutputStream out(outOpts);

	std::unordered_map<std::string, std::unordered_map<std::string, kmer>> uniqueKmersPerRecord;
	std::unordered_map<std::string, std::unordered_map<std::string, kmer>> uniqueKmersPerRecordRevComplement;


	{
	  watch.startNewLap("Get Unique K-mers");
		SeqInput reader(setUp.pars_.ioOptions_);
		reader.openIn();
		std::mutex uniqueKmersPerRecord_mut;
		std::function<void()> getUniqueKmers = [&reader,&uniqueKmersPerRecord, &uniqueKmersPerRecordRevComplement,
																						&uniqueKmersPerRecord_mut, &klen, &getReverseComp](){
			seqInfo seq;
			while (reader.readNextReadLock(seq)) {
				std::unordered_map<std::string, kmer> allKmers;
				std::unordered_map<std::string, kmer> allKmersRevComp;
				std::unordered_map<std::string, kmer> uniqueKmersForRecord;
				std::unordered_map<std::string, kmer> uniqueKmersForRecordRevComplement;
				if (len(seq) > klen + 1) {
					for (uint32_t pos = 0; pos < len(seq) - klen + 1; ++pos) {
						auto k = seq.seq_.substr(pos, klen);
						if(njh::in(k, allKmers)){
							allKmers[k].addPosition(pos);
						}else{
							allKmers[k] = kmer(k, pos);
						}
					}
				}
				for (const auto &k : allKmers) {
					if (1 == k.second.count_) {
						uniqueKmersForRecord.emplace(k);
					}
				}
				if(getReverseComp){
					seq.reverseComplementRead(false, true);
					if (len(seq) > klen + 1) {
						for (uint32_t pos = 0; pos < len(seq) - klen + 1; ++pos) {
							auto k = seq.seq_.substr(pos, klen);
							if(njh::in(k, allKmersRevComp)){
								allKmersRevComp[k].addPosition(pos);
							}else{
								allKmersRevComp[k] = kmer(k, pos);
							}
						}
					}
					for (const auto &k : allKmersRevComp) {
						if (1 == k.second.count_) {
							uniqueKmersForRecordRevComplement.emplace(k);
						}
					}
				}
				{
					std::lock_guard<std::mutex> lock(uniqueKmersPerRecord_mut);
					uniqueKmersPerRecord[seq.name_] = uniqueKmersForRecord;
					if(getReverseComp){
						uniqueKmersPerRecordRevComplement[seq.name_] = uniqueKmersForRecordRevComplement;
					}
				}
			}
		};

		njh::concurrent::runVoidFunctionThreaded(getUniqueKmers, numThreads);

	}
	std::unordered_map<std::string, std::unordered_map<std::string, kmer>> ref_uniqueKmersPerRecord;
	std::unordered_map<std::string, std::unordered_map<std::string, kmer>> ref_uniqueKmersPerRecordRevComplement;
	{
		//ref
	  watch.startNewLap("Get Unique K-mers");
		SeqInput reader(setUp.pars_.refIoOptions_);
		reader.openIn();
		std::mutex ref_uniqueKmersPerRecord_mut;
		std::function<void()> getUniqueKmers = [&reader,&ref_uniqueKmersPerRecord, &ref_uniqueKmersPerRecordRevComplement,
																						&ref_uniqueKmersPerRecord_mut, &klen, &getReverseComp](){
			seqInfo seq;
			while (reader.readNextReadLock(seq)) {
				std::unordered_map<std::string, kmer> allKmers;
				std::unordered_map<std::string, kmer> allKmersRevComp;
				std::unordered_map<std::string, kmer> uniqueKmersForRecord;
				std::unordered_map<std::string, kmer> uniqueKmersForRecordRevComplement;
				if (len(seq) > klen + 1) {
					for (uint32_t pos = 0; pos < len(seq) - klen + 1; ++pos) {
						auto k = seq.seq_.substr(pos, klen);
						if(njh::in(k, allKmers)){
							allKmers[k].addPosition(pos);
						}else{
							allKmers[k] = kmer(k, pos);
						}
					}
				}
				for (const auto &k : allKmers) {
					if (1 == k.second.count_) {
						uniqueKmersForRecord.emplace(k);
					}
				}
				if(getReverseComp){
					seq.reverseComplementRead(false, true);
					if (len(seq) > klen + 1) {
						for (uint32_t pos = 0; pos < len(seq) - klen + 1; ++pos) {
							auto k = seq.seq_.substr(pos, klen);
							if(njh::in(k, allKmersRevComp)){
								allKmersRevComp[k].addPosition(pos);
							}else{
								allKmersRevComp[k] = kmer(k, pos);
							}
						}
					}
					for (const auto &k : allKmersRevComp) {
						if (1 == k.second.count_) {
							uniqueKmersForRecordRevComplement.emplace(k);
						}
					}
				}
				{
					std::lock_guard<std::mutex> lock(ref_uniqueKmersPerRecord_mut);
					ref_uniqueKmersPerRecord[seq.name_] = uniqueKmersForRecord;
					if(getReverseComp){
						ref_uniqueKmersPerRecordRevComplement[seq.name_] = uniqueKmersForRecordRevComplement;
					}
				}
			}
		};

		njh::concurrent::runVoidFunctionThreaded(getUniqueKmers, numThreads);

	}
	watch.startNewLap("comparing");

	AllByAllPairFactory allFac(ref_uniqueKmersPerRecord.size(), uniqueKmersPerRecord.size());
	out << "query_record\tref_record\tr1TotalUniqueKmers\tr2TotalUniqueKmers\ttotalUniqueKmersShared\tshareScore\tshareScore2" << std::endl;
	auto query_recordNames = getVectorOfMapKeys(uniqueKmersPerRecord);
	njh::sort(query_recordNames);

	auto ref_recordNames = getVectorOfMapKeys(ref_uniqueKmersPerRecord);
	njh::sort(ref_recordNames);

	std::mutex outMut;

	std::function<void()> runComp = [&out,&outMut,&allFac,
																	 &uniqueKmersPerRecord,
																	 &ref_uniqueKmersPerRecord, &ref_uniqueKmersPerRecordRevComplement,
																	 &getReverseComp,&query_recordNames,
																	 &ref_recordNames](){
		AllByAllPairFactory::AllByAllPair between_pair;
		while(allFac.setNextPair(between_pair)){
			{
				uint32_t shared = 0;
				for(const auto & k :uniqueKmersPerRecord[query_recordNames[between_pair.col_]]){
					if(njh::in(k.first, ref_uniqueKmersPerRecord[ref_recordNames[between_pair.row_]])){
						++shared;
					}
				}
				{
					std::lock_guard<std::mutex> lock(outMut);
					out << query_recordNames[between_pair.col_]
							<< "\t" << ref_recordNames[between_pair.row_]
							<< "\t" << uniqueKmersPerRecord[query_recordNames[between_pair.col_]].size()
							<< "\t" << ref_uniqueKmersPerRecord[ref_recordNames[between_pair.row_]].size()
							<< "\t" << shared
							<< "\t" <<  static_cast<double>(shared)/std::min(uniqueKmersPerRecord[query_recordNames[between_pair.col_]].size(), ref_uniqueKmersPerRecord[ref_recordNames[between_pair.row_]].size())
					    << "\t" << (static_cast<double>(shared))/(uniqueKmersPerRecord[query_recordNames[between_pair.col_]].size() +  ref_uniqueKmersPerRecord[ref_recordNames[between_pair.row_]].size() - shared)
					<< std::endl;
				}
			}
			if(getReverseComp){
				uint32_t shared = 0;
				for(const auto & k :uniqueKmersPerRecord[query_recordNames[between_pair.col_]]){
					if(njh::in(k.first, ref_uniqueKmersPerRecordRevComplement[ref_recordNames[between_pair.row_]])){
						++shared;
					}
				}
				{
					std::lock_guard<std::mutex> lock(outMut);
					out << query_recordNames[between_pair.col_]
							<< "\t" << ref_recordNames[between_pair.row_] << "_revComp"
							<< "\t" << uniqueKmersPerRecord[query_recordNames[between_pair.col_]].size()
							<< "\t" << ref_uniqueKmersPerRecordRevComplement[ref_recordNames[between_pair.row_]].size()
							<< "\t" << shared
							<< "\t" <<  static_cast<double>(shared)/std::min(uniqueKmersPerRecord[query_recordNames[between_pair.col_]].size(), ref_uniqueKmersPerRecordRevComplement[ref_recordNames[between_pair.row_]].size())
					    << "\t" << (static_cast<double>(shared))/(uniqueKmersPerRecord[query_recordNames[between_pair.col_]].size() +  ref_uniqueKmersPerRecordRevComplement[ref_recordNames[between_pair.row_]].size() - shared)
					<< std::endl;
				}
			}
		}
	};

	njh::concurrent::runVoidFunctionThreaded(runComp, numThreads);



  if(setUp.pars_.verbose_){
  	watch.logLapTimes(std::cout, true, 6, true);
  	std::cout << "totalTime: " << watch.totalTimeFormatted(6) << std::endl;
  }
	return 0;
}


int kmerExpRunner::pairwiseWithinComparisonOfUniqueKmers(const njh::progutils::CmdArgs & inputCommands){
	uint32_t klen = 35;
	uint32_t numThreads = 1;
	bool getReverseComp = false;
	OutOptions outOpts(bfs::path(""), ".tab.txt");
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processReadInNames(true);
	setUp.processWritingOptions(outOpts);
	setUp.setOption(numThreads, "--numThreads", "Number of threads to use", njh::progutils::ProgramSetUp::CheckCase::GTEQ1);
	setUp.setOption(klen, "--klen", "kmer Length", true);
	setUp.setOption(getReverseComp, "--getReverseComp", "Compare the reverse complement as well");

	setUp.finishSetUp(std::cout);

  njh::stopWatch watch;
  watch.setLapName("Get Unique K-mers");
  OutputStream out(outOpts);

	SeqInput reader(setUp.pars_.ioOptions_);
	std::unordered_map<std::string, std::unordered_map<std::string, kmer>> uniqueKmersPerRecord;
	std::unordered_map<std::string, std::unordered_map<std::string, kmer>> uniqueKmersPerRecordRevComplement;

	reader.openIn();


	std::mutex uniqueKmersPerRecord_mut;

	std::function<void()> getUniqueKmers = [&reader,&uniqueKmersPerRecord, &uniqueKmersPerRecordRevComplement,
																					&uniqueKmersPerRecord_mut, &klen, &getReverseComp](){
		seqInfo seq;


		while (reader.readNextReadLock(seq)) {
			std::unordered_map<std::string, kmer> allKmers;
			std::unordered_map<std::string, kmer> allKmersRevComp;
			std::unordered_map<std::string, kmer> uniqueKmersForRecord;
			std::unordered_map<std::string, kmer> uniqueKmersForRecordRevComplement;
			if (len(seq) > klen + 1) {
				for (uint32_t pos = 0; pos < len(seq) - klen + 1; ++pos) {
					auto k = seq.seq_.substr(pos, klen);
					if(njh::in(k, allKmers)){
						allKmers[k].addPosition(pos);
					}else{
						allKmers[k] = kmer(k, pos);
					}
				}
			}
			for (const auto &k : allKmers) {
				if (1 == k.second.count_) {
					uniqueKmersForRecord.emplace(k);
				}
			}
			if(getReverseComp){
				seq.reverseComplementRead(false, true);
				if (len(seq) > klen + 1) {
					for (uint32_t pos = 0; pos < len(seq) - klen + 1; ++pos) {
						auto k = seq.seq_.substr(pos, klen);
						if(njh::in(k, allKmersRevComp)){
							allKmersRevComp[k].addPosition(pos);
						}else{
							allKmersRevComp[k] = kmer(k, pos);
						}
					}
				}
				for (const auto &k : allKmersRevComp) {
					if (1 == k.second.count_) {
						uniqueKmersForRecordRevComplement.emplace(k);
					}
				}
			}
			{
				std::lock_guard<std::mutex> lock(uniqueKmersPerRecord_mut);
				uniqueKmersPerRecord[seq.name_] = uniqueKmersForRecord;
				if(getReverseComp){
					uniqueKmersPerRecordRevComplement[seq.name_] = uniqueKmersForRecordRevComplement;
				}
			}
		}
	};

	njh::concurrent::runVoidFunctionThreaded(getUniqueKmers, numThreads);

	watch.startNewLap("comparing");
	PairwisePairFactory pFac(uniqueKmersPerRecord.size());


	out << "record1\trecord2\tr1TotalUniqueKmers\tr2TotalUniqueKmers\ttotalUniqueKmersShared\tshareScore\tshareScore2" << std::endl;
	auto recordNames = getVectorOfMapKeys(uniqueKmersPerRecord);
	njh::sort(recordNames);

	std::mutex outMut;

	std::function<void()> runComp = [&out,&outMut,&pFac,&uniqueKmersPerRecord, &uniqueKmersPerRecordRevComplement, &getReverseComp,&recordNames](){
		PairwisePairFactory::PairwisePair pair;
		while(pFac.setNextPair(pair)){

			{
				uint32_t shared = 0;
				for(const auto & k :uniqueKmersPerRecord[recordNames[pair.col_]]){
					if(njh::in(k.first, uniqueKmersPerRecord[recordNames[pair.row_]])){
						++shared;
					}
				}
				{
					std::lock_guard<std::mutex> lock(outMut);
					out << recordNames[pair.col_]
							<< "\t" << recordNames[pair.row_]
							<< "\t" << uniqueKmersPerRecord[recordNames[pair.col_]].size()
							<< "\t" << uniqueKmersPerRecord[recordNames[pair.row_]].size()
							<< "\t" << shared
							<< "\t" <<  static_cast<double>(shared)/std::min(uniqueKmersPerRecord[recordNames[pair.col_]].size(), uniqueKmersPerRecord[recordNames[pair.row_]].size())
					    << "\t" << (static_cast<double>(shared))/(uniqueKmersPerRecord[recordNames[pair.col_]].size() +  uniqueKmersPerRecord[recordNames[pair.row_]].size() - shared)
					<< std::endl;
				}
			}
			if(getReverseComp){
				uint32_t shared = 0;
				for(const auto & k :uniqueKmersPerRecord[recordNames[pair.col_]]){
					if(njh::in(k.first, uniqueKmersPerRecordRevComplement[recordNames[pair.row_]])){
						++shared;
					}
				}
				{
					std::lock_guard<std::mutex> lock(outMut);
					out << recordNames[pair.col_]
							<< "\t" << recordNames[pair.row_] << "_revComp"
							<< "\t" << uniqueKmersPerRecord[recordNames[pair.col_]].size()
							<< "\t" << uniqueKmersPerRecordRevComplement[recordNames[pair.row_]].size()
							<< "\t" << shared
							<< "\t" <<  static_cast<double>(shared)/std::min(uniqueKmersPerRecord[recordNames[pair.col_]].size(), uniqueKmersPerRecordRevComplement[recordNames[pair.row_]].size())
					    << "\t" << (static_cast<double>(shared))/(uniqueKmersPerRecord[recordNames[pair.col_]].size() +  uniqueKmersPerRecordRevComplement[recordNames[pair.row_]].size() - shared)
					<< std::endl;
				}
			}
		}
	};

	njh::concurrent::runVoidFunctionThreaded(runComp, numThreads);



  if(setUp.pars_.verbose_){
  	watch.logLapTimes(std::cout, true, 6, true);
  	std::cout << "totalTime: " << watch.totalTimeFormatted(6) << std::endl;
  }
	return 0;
}



}  // namespace njhseq
