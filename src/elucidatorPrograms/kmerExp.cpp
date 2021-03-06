/*
 * kmerExpRunner.cpp
 *
 *  Created on: Jan 25, 2015
 *      Author: nickhathaway
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


namespace njhseq {
kmerExpRunner::kmerExpRunner()
    : njh::progutils::ProgramRunner(
          {
					 addFunc("profileScanningKmerDist", profileScanningKmerDist, false),
					 addFunc("profileKmerAccerlation", profileKmerAccerlation, false),
					 addFunc("getNewScanKmerDist", getNewScanKmerDist, false),
					 addFunc("readingDistanceCheck", readingDistanceCheck, false),
					 addFunc("getKmerDist", getKmerDist, false),
					 addFunc("writingDistanceCheck", writingDistanceCheck, false),
					 addFunc("writeKmerAccerlation",writeKmerAccerlation, false),
					 addFunc("kDistVsNucDist",kDistVsNucDist, false),
					 addFunc("clostestKmerDist",clostestKmerDist, false),
					 addFunc("getKmerDistStatsMultiple",getKmerDistStatsMultiple, false),
					 addFunc("getKmerDistStats",getKmerDistStats, false),
					 addFunc("scaningKmerDist",scaningKmerDist, false),
					 addFunc("kDist",kDist, false),
					 addFunc("scoveViaKmers",scoveViaKmers, false),
					 addFunc("kmerRevVsForDist",kmerRevVsForDist, false),
					 addFunc("profileLargeKmerIndex",profileLargeKmerIndex, false),
					 addFunc("profileSharedKmerBlocks",profileSharedKmerBlocks, false),
					 addFunc("pidVsKmers", pidVsKmers, false),
					 addFunc("randomSampleKmerCompare", randomSampleKmerCompare, false),
					 addFunc("findingMinimumKLenForNoRedundantKmers", findingMinimumKLenForNoRedundantKmers, false),
					 addFunc("microsatsKmerSearch", microsatsKmerSearch, false),
					 addFunc("genomeKmerCompare", genomeKmerCompare, false),
					 addFunc("kmerSearch", kmerSearch, false),
					 addFunc("convertKmerSearchToBinaryMatrix", convertKmerSearchToBinaryMatrix, false),
					 addFunc("generateCountsTable", generateCountsTable, false),
					 addFunc("getBestKmerDist", getBestKmerDist, false),
					 addFunc("writeKmerSimDistanceMatrix", writeKmerSimDistanceMatrix, false),
					 addFunc("findUniqKmersBetweenSeqs", findUniqKmersBetweenSeqs, false),
					 //
           },
          "kmerExp") {}

int kmerExpRunner::findUniqKmersBetweenSeqs(const njh::progutils::CmdArgs & inputCommands) {
	seqInfo seq1("seq1");
	seqInfo seq2("seq2");
	OutOptions outOpts(bfs::path(""), ".txt");
	uint32_t klen = 40;
	bool sortByPosition = false;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processSeq(seq1, "--seq1", "seq1", true);
	setUp.processSeq(seq2, "--seq2", "seq2", true);
	setUp.processWritingOptions(outOpts);
	setUp.setOption(klen, "--klen", "Kmer Length");
	setUp.setOption(sortByPosition, "--sortByPosition", "Sort By Position");

	setUp.finishSetUp(std::cout);

	kmerInfo seq1Info(seq1.seq_, klen, false);
	kmerInfo seq2Info(seq2.seq_, klen, false);

	OutputStream out(outOpts);

	out << "seq\tuniqKmer\tposition" << std::endl;
	if(sortByPosition){
		{
			std::map<uint32_t, std::vector<std::string>> kmersByPositon;
			for (const auto & k : seq1Info.kmers_) {
				if (!njh::in(k.first, seq2Info.kmers_)) {
					for (const auto & pos : k.second.positions_) {
						kmersByPositon[pos].emplace_back(k.first);
					}
				}
			}
			for (const auto & pos : kmersByPositon) {
				for (const auto & k : pos.second) {
					out << seq1.name_ << "\t" << k << "\t" << pos.first << std::endl;
				}
			}
		}
		{
			std::map<uint32_t, std::vector<std::string>> kmersByPositon;
			for (const auto & k : seq2Info.kmers_) {
				if (!njh::in(k.first, seq1Info.kmers_)) {
					for (const auto & pos : k.second.positions_) {
						kmersByPositon[pos].emplace_back(k.first);
					}
				}
			}
			for (const auto & pos : kmersByPositon) {
				for (const auto & k : pos.second) {
					out << seq2.name_ << "\t" << k << "\t" << pos.first << std::endl;
				}
			}
		}
	}else{
		for(const auto & k : seq1Info.kmers_){
			if(!njh::in(k.first, seq2Info.kmers_)){
				for(const auto & pos : k.second.positions_){
					out << seq1.name_
							<< "\t" << k.first
							<< "\t"<< pos << std::endl;
				}
			}
		}
		for(const auto & k : seq2Info.kmers_){
			if(!njh::in(k.first, seq1Info.kmers_)){
				for(const auto & pos : k.second.positions_){
					out << seq2.name_
							<< "\t" << k.first
							<< "\t"<< pos << std::endl;
				}
			}
		}
	}


	return 0;

}


std::string perfectKmerHtmlIndexFileStr = "<!DOCTYPE html>\n"
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
		"  d3.json(\"perfect.json\", function(error, data) {\n"
		"\n"
		"\n"
		"  x.domain([0, 180000]).nice();\n"
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
		"      		  var seqStop = d.start_+ d.size_;\n"
		"	    	  tooltip.node().innerHTML = \"name: \"  + bars[j].parentNode.__data__.name_ + \"<br>refRange: \" + d.refStart_.toString() + \":\"  + refStop.toString() + \"<br>seqRange: \"  + d.start_.toString() + \":\" + seqStop.toString() ;\n"
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

int kmerExpRunner::profileSharedKmerBlocks(const njh::progutils::CmdArgs & inputCommands){
	seqSetUp setUp(inputCommands);
	uint32_t kLen = 18;
	std::string seqsDir = "";
	uint32_t allowableMissing = 1;
	setUp.processVerbose();
	setUp.setOption(kLen, "--kLen", "kmer Len");
	setUp.setOption(allowableMissing, "--allowableMissing", "allowable Missing for near perfect outputs");
	setUp.processSeq(true);
	bool seqsDirPresent = setUp.setOption(seqsDir, "--seqsDir", "Name of the directory that holds several fasta files to profile");
	setUp.processDefaultReader(seqsDirPresent);
	setUp.processDirectoryOutputName("_profileSharedKmerBlocks_TODAY", true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	std::ofstream perfectCompResultsJsonFile;
	openTextFile(perfectCompResultsJsonFile, setUp.pars_.directoryName_ + "perfect", ".json",
			setUp.pars_.ioOptions_.out_);
	std::ofstream perfectCompResultsHtmlFile;
	openTextFile(perfectCompResultsHtmlFile, setUp.pars_.directoryName_ + "perfect", ".html",
			setUp.pars_.ioOptions_.out_);
	perfectCompResultsHtmlFile << perfectKmerHtmlIndexFileStr;
	std::ofstream perfectCompResultsFile;
	openTextFile(perfectCompResultsFile, setUp.pars_.directoryName_ + "perfect", ".fasta",
			setUp.pars_.ioOptions_.out_);

	std::ofstream nearPerfectCompResultsFile;
	openTextFile(nearPerfectCompResultsFile, setUp.pars_.directoryName_ + "nearPerfect", ".fasta",
			setUp.pars_.ioOptions_.out_);

	std::ofstream compResultsJsonFile;
	openTextFile(compResultsJsonFile, setUp.pars_.directoryName_ + "fractions", ".json",
			setUp.pars_.ioOptions_.out_);
	std::ofstream compResultsFile;
	openTextFile(compResultsFile, setUp.pars_.directoryName_ + "fractions", ".tab.txt",
			setUp.pars_.ioOptions_.out_);

	njh::stopWatch watch;
	watch.setLapName("kmer index");
	kmerInfo refGenInfo = kmerInfo(setUp.pars_.seqObj_.seqBase_.seq_, kLen, false);
	std::unordered_map<std::string, KmersSharedBlocks> seqs;
	std::unordered_map<std::string, KmersSharedBlocks> perfectSeqs;
	std::unordered_map<std::string, KmersSharedBlocks> nearPerfectSeqs;
	if(seqsDirPresent){
		auto files = njh::files::listAllFiles(seqsDir, false, {std::regex{".*.fasta"}});
		for (const auto & f : files) {
			SeqIOOptions opts = SeqIOOptions::genFastaIn(f.first.string());
			SeqInput reader(opts);
			reader.openIn();
			auto inReads = reader.readAllReads<readObject>();
			seqs.emplace(inReads.front().seqBase_.name_,
					KmersSharedBlocks(inReads.front().seqBase_.name_,
							inReads.front().seqBase_.seq_, kLen));
		}
	}else{
		SeqInput reader(setUp.pars_.ioOptions_);
		reader.openIn();
		auto inReads = reader.readAllReads<readObject>();
		for(const auto & seq : inReads){
			seqs.emplace(seq.seqBase_.name_,
					KmersSharedBlocks(seq.seqBase_.name_,
							seq.seqBase_.seq_, kLen));
		}
	}


	perfectSeqs = seqs;
	nearPerfectSeqs = seqs;
	watch.startNewLap("Search All");
	for (const auto & pos : iter::range(
			setUp.pars_.seqObj_.seqBase_.seq_.size() - kLen + 1)) {
		auto sub = setUp.pars_.seqObj_.seqBase_.seq_.substr(pos, kLen);
		for (auto & seq : seqs) {
			auto search = seq.second.kInfo_->kmers_.find(sub);
			//if both sequences contain the kmer and the kmer is unique in both sequences, add them
			if (search != seq.second.kInfo_->kmers_.end()
					&& search->second.positions_.size() == 1
					&& refGenInfo.kmers_[search->first].positions_.size() == 1) {
				seq.second.addComp(refGenInfo.kmers_[search->first].positions_.front(),
										search->second.positions_.front());
			}
		}
	}
	//
	for (auto & seq : seqs) {
		seq.second.finish();
	}


	watch.startNewLap("Print Test Shared Kmers");
	if(setUp.pars_.verbose_){
		for(const auto & g : seqs){
			auto keys = getVectorOfMapKeys(g.second.kComps_);
			njh::sort(keys);
			uint32_t matchCount = 0;
			for (const auto & compKey : keys) {
				std::cout << g.first << std::endl;
				std::cout << "\t" << g.second.kComps_.at(compKey).refStart_ << std::endl;
				std::cout << "\t" << g.second.kComps_.at(compKey).start_ << std::endl;
				std::cout << "\t" << g.second.kComps_.at(compKey).size_ << std::endl;
				++matchCount;
			}
		}
	}
	watch.startNewLap("Count All Shared Kmers");
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
	for(const auto & g : perfectSeqs){
		Json::Value currentGen;
		currentGen["name_"] = g.first;
		currentGen["kLen_"] = kLen;
		Json::Value & kComps = currentGen["kComps_"];
		for(const auto & k : g.second.kComps_ ){
			kComps.append(njh::json::toJson(k.second));
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
	watch.logLapTimes(std::cout, true, 6, true);
	return 0;
}

int kmerExpRunner::profileLargeKmerIndex(const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	uint32_t seqLen = 170000;
	uint32_t kmerLen = 18;
	setUp.setOption(seqLen, "--seqLen,-l", "Sequence Length");
	setUp.setOption(seqLen, "--kmerLen,-k", "Sequence Length");
	setUp.finishSetUp(std::cout);
	njh::randomGenerator gen;
	njh::stopWatch watch;
	watch.setLapName("gen");
	auto str = simulation::evenRandStr(seqLen, { 'A', 'C', 'G', 'T' }, gen);
	watch.startNewLap("kmer index");
	watch.logLapTimes(std::cout, true, 6, false);
	kmerInfo kInfo(str, kmerLen, false);
	watch.logLapTimes(std::cout, true, 6, true);
	return 0;
}




int kmerExpRunner::kmerRevVsForDist(const njh::progutils::CmdArgs & inputCommands){
	seqSetUp setUp(inputCommands);
	setUp.pars_.colOpts_.kmerOpts_.kLength_ = 10;
	uint32_t numThreads = 1;
	setUp.setOption(numThreads, "--threads,-t", "Number of Threads To use");
	setUp.processDefaultReader(true);
	setUp.processKmerLenOptions();
	setUp.processDebug();
	setUp.processVerbose();
  setUp.finishSetUp(std::cout);
  auto reads = createKmerReadVec(setUp.pars_.ioOptions_,setUp.pars_.colOpts_.kmerOpts_.kLength_, true);
	std::function<
			std::pair<double, double>(const std::unique_ptr<seqWithKmerInfo> &,
					const std::unique_ptr<seqWithKmerInfo> &)> disFun =
			[](const std::unique_ptr<seqWithKmerInfo> & read1,
					const std::unique_ptr<seqWithKmerInfo> & read2) {
				auto dist = read1->compareKmers(*read2);
				auto distRev = read1->compareKmersRevComp(*read2);
				return std::pair<double,double> {dist.second, distRev.second};
			};
	auto distances = getDistance(reads, numThreads, disFun);
	for (const auto & pos : iter::range(reads.size())) {
		for (const auto & subPos : iter::range(pos)) {
			std::cout << distances[pos][subPos].first << ","
					<< distances[pos][subPos].second << " ";
		}
		std::cout << '\n';
	}
  return 0;
}

//scaningKmerDist
int kmerExpRunner::randomSampleKmerCompare(const njh::progutils::CmdArgs & inputCommands) {
  uint32_t kmerStart = 2;
  uint32_t kmerStop = 5;
  uint32_t numThreads = 2;
  uint32_t sampleAmount = 600;
  uint64_t seed = std::numeric_limits<uint64_t>::max();
	seqSetUp setUp(inputCommands);
	setUp.pars_.ioOptions_.out_.outFilename_ = "pidVsKSim.tab.txt";
	setUp.pars_.ioOptions_.lowerCaseBases_ = "upper";
	setUp.processDefaultReader(true);
	setUp.pars_.ioOptions_.out_.outExtention_ = ".tab.txt";
	setUp.pars_.ioOptions_.out_.outFileFormat_ = "tab";
	setUp.processAlignerDefualts();
  setUp.setOption(numThreads, "--numThreads", "Number of Threads to Use");
  setUp.setOption(kmerStart, "--kmerLenStart", "Length for kmers to start at");
  setUp.setOption(kmerStop, "--kmerLenStop", "Length for kmers to stop at");
  setUp.setOption(sampleAmount, "--sampleAmount", "Sample Amount for the two groups");
  setUp.setOption(seed, "--seed", "seed for randomness");
  setUp.finishSetUp(std::cout);
  setUp.startARunLog(setUp.pars_.directoryName_);
	if (kmerStop < kmerStart) {
		std::stringstream ss;
		ss << "Kmer stop must be greater than kmerStart" << std::endl;
		ss << "KmerStart: " << kmerStart << std::endl;
		ss << "KmerStop: " << kmerStop << std::endl;
		throw std::runtime_error { ss.str() };
	}
	if (setUp.pars_.ioOptions_.out_.outExists()
			&& !setUp.pars_.ioOptions_.out_.overWriteFile_) {
		std::stringstream ss;
		ss << "File: "
				<< setUp.pars_.ioOptions_.out_.outName()
				<< " already exists, use -overWrite to over write it" << std::endl;
		throw std::runtime_error { ss.str() };
	}
  SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<readObject>();
	if(inReads.size() * 2 < sampleAmount){
		std::stringstream ss;
		ss << "Sample Amount is too great" << std::endl;
		ss << "sampleAmount: " << sampleAmount << std::endl;
		ss << "inReads.size() * 2: " << inReads.size() * 2 << std::endl;
		throw std::runtime_error { ss.str() };
	}
	std::vector<uint32_t> positions(inReads.size(), 0);
	njh::iota<uint32_t>(positions, 0);
	njh::randomGenerator gen;
	if(std::numeric_limits<uint64_t>::max() != seed){
		gen.seedNum(seed);
	}
	auto selPositions = gen.unifRandSelectionVec(positions,sampleAmount * 2, false);
	std::vector<readObject> group1;
	std::vector<readObject> group2;
	for(auto pos : selPositions){
		if(group1.size() < sampleAmount){
			group1.emplace_back(inReads[pos]);
		}else{
			group2.emplace_back(inReads[pos]);
		}
	}
  uint64_t maxSize = 0;
  readVec::getMaxLength(group1,maxSize);
  readVec::getMaxLength(group2,maxSize);
  std::vector<std::unique_ptr<seqWithKmerInfo>> reads1 = createKmerReadVec(group1);
  std::vector<std::unique_ptr<seqWithKmerInfo>> reads2 = createKmerReadVec(group2);
  VecStr names;
  aligner alignerObj(maxSize,setUp.pars_.gapInfo_, setUp.pars_.scoring_, setUp.pars_.colOpts_.alignOpts_.countEndGaps_);
  alignerObj.processAlnInfoInput(setUp.pars_.alnInfoDirName_, setUp.pars_.verbose_);
  table distTab{VecStr{"read1", "read2", "kLen", "kmerDist", "id"}};
	std::unordered_map<std::string, std::unique_ptr<aligner>> aligners;
	std::mutex  alignerLock;
	std::function<
			comparison(const std::unique_ptr<seqWithKmerInfo> &,
					const std::unique_ptr<seqWithKmerInfo> &,
					std::unordered_map<std::string, std::unique_ptr<aligner>>&, aligner&)> getMismatchesFunc =
			[&alignerLock](const std::unique_ptr<seqWithKmerInfo> & read1, const std::unique_ptr<seqWithKmerInfo> & read2,
					std::unordered_map<std::string, std::unique_ptr<aligner>>& aligners,
					aligner &alignerObj) {
				alignerLock.lock();
				auto threadId = estd::to_string(std::this_thread::get_id());
				//std::cout << threadId<< std::endl;
				if(aligners.find(threadId) == aligners.end()) {
					aligners.emplace(threadId, std::make_unique<aligner>(alignerObj));
				}
				alignerLock.unlock();
				aligners.at(threadId)->alignCache(getSeqBase(read1),getSeqBase(read2), false);
				aligners.at(threadId)->profilePrimerAlignment(getSeqBase(read1), getSeqBase(read2));
				return aligners.at(threadId)->comp_;
			};
	std::vector<uint32_t> group1Positions(group1.size(), 0);
	njh::iota<uint32_t>(group1Positions, 0);
	std::unordered_map<uint32_t, std::vector<comparison>> allComps;
	std::mutex allCompsMut;
	table timingTab(VecStr{"seqLength", "n", "type", "time"});
	{
		njh::stopWatch watch;
		njh::concurrent::LockableQueue<uint32_t> group1Queue(group1Positions);
		auto alignFunc = [&alignerLock, &allCompsMut, &allComps](njh::concurrent::LockableQueue<uint32_t> & posQueue,
				const std::vector<std::unique_ptr<seqWithKmerInfo>> & r1,
				const std::vector<std::unique_ptr<seqWithKmerInfo>> & r2,
				std::unordered_map<std::string, std::unique_ptr<aligner>>& aligners,
				aligner &alignerObj){
			alignerLock.lock();
			auto threadId = estd::to_string(std::this_thread::get_id());
			//std::cout << threadId<< std::endl;
			if(aligners.find(threadId) == aligners.end()) {
				aligners.emplace(threadId, std::make_unique<aligner>(alignerObj));
			}
			alignerLock.unlock();
			uint32_t pos = 0;
			std::unordered_map<uint32_t, std::vector<comparison>> comps;
			while(posQueue.getVal(pos)){
				std::vector<comparison> r2Comps;
				for(const auto & r2Pos : iter::range(r2.size())){
					aligners.at(threadId)->alignCache(getSeqBase(r2[r2Pos]),getSeqBase(r1[pos]), false);
					aligners.at(threadId)->profilePrimerAlignment(r2[r2Pos], getSeqBase(r1[pos]));
					r2Comps.emplace_back(aligners.at(threadId)->comp_);
				}
				comps[pos] = r2Comps;
			}
			{
				std::lock_guard<std::mutex> lock(allCompsMut);
				for(const auto & comp : comps){
					allComps[comp.first] = comp.second;
				}
			}
		};
		std::vector<std::thread> threads;
		for (uint32_t t = 0; t < numThreads; ++t) {
			threads.emplace_back(
					std::thread(alignFunc, std::ref(group1Queue), std::cref(reads1),
							std::cref(reads2), std::ref(aligners), std::ref(alignerObj)));
		}
		for(auto & t : threads){
			t.join();
		}

		timingTab.addRow(maxSize, sampleAmount, "alignment", watch.totalTime());
	}
	std::mutex distTabMut;
	{
		auto kDistFunc = [&distTabMut,&distTab, &allComps](njh::concurrent::LockableQueue<uint32_t> & posQueue,
				const std::vector<std::unique_ptr<seqWithKmerInfo>> & r1,
				const std::vector<std::unique_ptr<seqWithKmerInfo>> & r2,
				uint32_t k){
			uint32_t pos = 0;
			std::unordered_map<uint32_t, std::vector<double>> comps;
			while(posQueue.getVal(pos)){
				std::vector<double> r2Comps;
				for(const auto & r2Pos : iter::range(r2.size())){
					auto dist = r1[pos]->compareKmers(getRef(r2[r2Pos]));
					r2Comps.emplace_back(dist.second);
				}
				comps[pos] = r2Comps;
			}
			{
				std::lock_guard<std::mutex> lock(distTabMut);
				for(const auto & comp : comps){
					for(const auto & compPos : iter::range(comp.second.size())){
						distTab.addRow(r1[comp.first]->seqBase_.name_,
																			r2[compPos]->seqBase_.name_, k,
																			comp.second[compPos],
																			allComps[comp.first][compPos].distances_.eventBasedIdentity_);
					}
				}
			}
		};
		for (uint32_t k = kmerStart; k < kmerStop + 1; ++k) {
			allSetKmers(reads1, k, false);
			njh::stopWatch watch;
			allSetKmers(reads2, k, false);
			njh::concurrent::LockableQueue<uint32_t> group1Queue(group1Positions);
			std::vector<std::thread> threads;
			for (uint32_t t = 0; t < numThreads; ++t) {
				threads.emplace_back(
						std::thread(kDistFunc, std::ref(group1Queue), std::cref(reads1),
								std::cref(reads2), k));
			}
			for(auto & t : threads){
				t.join();
			}
			timingTab.addRow(maxSize, sampleAmount, "k-" + estd::to_string(k), watch.totalTime());
		}
	}
	auto timmingOpts = setUp.pars_.ioOptions_.out_;
	timmingOpts.outFilename_ = timmingOpts.outFilename_.string() + "_timming";
  distTab.outPutContents(TableIOOpts(setUp.pars_.ioOptions_.out_, "\t", distTab.hasHeader_));
  timingTab.outPutContents(TableIOOpts(timmingOpts, "\t", distTab.hasHeader_));
  for(const auto & align : aligners){
  	alignerObj.alnHolder_.mergeOtherHolder(align.second->alnHolder_);
  }

  alignerObj.processAlnInfoOutput(setUp.pars_.alnInfoDirName_, setUp.pars_.verbose_);
  return 0;
}

int kmerExpRunner::pidVsKmers(const njh::progutils::CmdArgs & inputCommands) {
  uint32_t kmerStart = 2;
  uint32_t kmerStop = 4;
  uint32_t numThreads = 2;//profileErrors
	seqSetUp setUp(inputCommands);
	setUp.pars_.ioOptions_.lowerCaseBases_ = "upper";
	setUp.pars_.ioOptions_.out_.outFilename_ = "pidVsKSim.tab.txt";
	setUp.processDefaultReader(true);
	setUp.pars_.ioOptions_.out_.outExtention_ = ".tab.txt";
	setUp.pars_.ioOptions_.out_.outFileFormat_ = "tab";
	setUp.processAlignerDefualts();
  setUp.setOption(numThreads, "--numThreads", "Number of Threads to Use");
  setUp.setOption(kmerStart, "--kmerLenStart", "Length for kmers to start at");
  setUp.setOption(kmerStop, "--kmerLenStop", "Length for kmers to stop at");
  setUp.finishSetUp(std::cout);
  setUp.startARunLog(setUp.pars_.directoryName_);
	if (kmerStop < kmerStart) {
		std::stringstream ss;
		ss << "Kmer stop must be greater than kmerStart" << std::endl;
		ss << "KmerStart: " << kmerStart << std::endl;
		ss << "KmerStop: " << kmerStop << std::endl;
		throw std::runtime_error { ss.str() };
	}
	if (setUp.pars_.ioOptions_.out_.outExists()
			&& !setUp.pars_.ioOptions_.out_.overWriteFile_) {
		std::stringstream ss;
		ss << "File: "
				<< setUp.pars_.ioOptions_.out_.outName()
				<< " already exists, use -overWrite to over write it" << std::endl;
		throw std::runtime_error { ss.str() };
	}
  SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<readObject>();
  uint64_t maxSize = 0;
  readVec::getMaxLength(inReads,maxSize);
  std::vector<std::unique_ptr<seqWithKmerInfo>> reads;
  VecStr names;
  names.reserve(inReads.size() + 1);
  for(const auto & read : inReads){
  	reads.emplace_back(std::make_unique<seqWithKmerInfo>(read.seqBase_));
  }
  aligner alignerObj(maxSize,setUp.pars_.gapInfo_, setUp.pars_.scoring_, setUp.pars_.colOpts_.alignOpts_.countEndGaps_);

  table distTab{VecStr{"read1", "read2", "kLen", "kmerDist", "id"}};
	std::unordered_map<std::string, std::unique_ptr<aligner>> aligners;
	std::mutex alignerLock;
	std::function<
			comparison(const std::unique_ptr<seqWithKmerInfo> &,
					const std::unique_ptr<seqWithKmerInfo> &,
					std::unordered_map<std::string, std::unique_ptr<aligner>>&, aligner&)> getMismatchesFunc =
			[&alignerLock](const std::unique_ptr<seqWithKmerInfo> & read1, const std::unique_ptr<seqWithKmerInfo> & read2,
					std::unordered_map<std::string, std::unique_ptr<aligner>>& aligners,
					aligner &alignerObj) {
				alignerLock.lock();
				auto threadId = estd::to_string(std::this_thread::get_id());
				//std::cout << threadId<< std::endl;
				if(aligners.find(threadId) == aligners.end()) {
					aligners.emplace(threadId, std::make_unique<aligner>(alignerObj));
				}
				alignerLock.unlock();
				aligners.at(threadId)->alignCache(getSeqBase(read1),getSeqBase(read2), false);
				aligners.at(threadId)->profilePrimerAlignment(getSeqBase(read1), getSeqBase(read2));
				return aligners.at(threadId)->comp_;
			};

	auto distances = getDistanceNonConst(reads, numThreads, getMismatchesFunc,
			aligners, alignerObj);
	std::function<
			double(const std::unique_ptr<seqWithKmerInfo> &,
					const std::unique_ptr<seqWithKmerInfo> &)> disFun =
			[](const std::unique_ptr<seqWithKmerInfo> & read1,
					const std::unique_ptr<seqWithKmerInfo> & read2) {
				auto dist = read1->compareKmers(*read2);
				return dist.second;
			};
	for (uint32_t k = kmerStart; k < kmerStop + 1; ++k) {
		allSetKmers(reads, k, false);
		auto kDists = getDistance(reads, numThreads, disFun);
		for (const auto & rowPos : iter::range(reads.size())) {
			for (const auto & colPos : iter::range(rowPos)) {
				distTab.addRow(reads[rowPos]->seqBase_.name_,
						reads[colPos]->seqBase_.name_, k, kDists[rowPos][colPos],
						distances[rowPos][colPos].distances_.eventBasedIdentity_);
			}
		}
	}

  distTab.outPutContents(TableIOOpts(setUp.pars_.ioOptions_.out_, "\t", distTab.hasHeader_));
  return 0;
}
//tempGetGV
int kmerExpRunner::scoveViaKmers(const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	setUp.pars_.ioOptions_.out_.outFilename_ = "kDist.tab.txt";
	double cutOff = 0.75;
	uint32_t sizeCutOff = 1;
	setUp.setOption(sizeCutOff, "--sizeCutOff", "Cluster Size Cut Off");
	setUp.setOption(cutOff, "--cutOff", "cutOff");
	setUp.processDefaultReader(true);
	setUp.processDirectoryOutputName(true);
	setUp.processAlignerDefualts();
	setUp.processRefFilename(true);
  setUp.finishSetUp(std::cout);
  setUp.startARunLog(setUp.pars_.directoryName_);
  SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<readObject>();
  uint64_t maxSize = 0;
  readVec::getMaxLength(inReads,maxSize);
  auto refs = SeqInput::getReferenceSeq(setUp.pars_.refIoOptions_, maxSize);
  std::vector<std::unique_ptr<seqWithKmerInfo>> refReads;
  for(const auto & read : refs){
  	refReads.emplace_back(std::make_unique<seqWithKmerInfo>(read.seqBase_));
  }
  std::vector<std::unique_ptr<seqWithKmerInfo>> reads;
  VecStr names;
  names.reserve(inReads.size() + 1);
  for(const auto & read : inReads){
  	reads.emplace_back(std::make_unique<seqWithKmerInfo>(read.seqBase_));
  }
  allSetKmers(reads, setUp.pars_.colOpts_.kmerOpts_.kLength_,false);
  allSetKmers(refReads, setUp.pars_.colOpts_.kmerOpts_.kLength_, false);
  aligner alignerObj(maxSize,setUp.pars_.gapInfo_, setUp.pars_.scoring_, setUp.pars_.colOpts_.alignOpts_.countEndGaps_);

  table distTab{VecStr{"kmerDist", "alnScore", "color"}};
	std::vector<double> kDists;
	std::vector<int32_t> scores;
	kDists.reserve(refReads.size());
	scores.reserve(refReads.size());
  for(const auto & readPos : iter::range(reads.size())){
  	if(readPos%10 == 0){
  		std::cout << "on " << readPos << " of " << reads.size() << "\r";
  		std::cout.flush();
  	}
  	auto & read = reads[readPos];
  	kDists.clear();
  	scores.clear();
  	for(const auto & ref : refReads){
  		auto dist = read->compareKmers(*ref);
  		alignerObj.alignCache(ref->seqBase_, read->seqBase_, setUp.pars_.local_);
  		kDists.emplace_back(dist.second);
  		scores.emplace_back(alignerObj.parts_.score_);
  	}
  	auto maxKDist = vectorMaximum(kDists);
  	auto maxScore = vectorMaximum(scores);
  	std::string outColor = "";
  	for(const auto & pos : iter::range(refReads.size())){
  		if(kDists[pos] == maxKDist && scores[pos] == maxScore){
  			outColor = "#14B814";
  		}else if (kDists[pos] != maxKDist && scores[pos] == maxScore){
  			outColor = "#115BEE";
  		}else if (kDists[pos] == maxKDist && scores[pos] != maxScore){
  			outColor = "#F83AB9";
  		}else{
  			outColor = "#000000";
  		}
  		distTab.content_.emplace_back(toVecStr(kDists[pos], scores[pos], outColor));
  	}
  }
  std::cout << std::endl;
  distTab.outPutContents(TableIOOpts(OutOptions(setUp.pars_.ioOptions_.out_.outFilename_, ".tab.txt"), "\t", distTab.hasHeader_));
  return 0;
}//scoveViaKmers

int kmerExpRunner::kDist(const njh::progutils::CmdArgs & inputCommands){
	seqSetUp setUp(inputCommands);
	std::string seq1 = "";
	std::string seq2 = "";
	uint32_t kLength = 5;
	bool getReverseDistance = false;
	setUp.processSeq(seq1, "-seq1", "Sequence 1", true);
	setUp.processSeq(seq2, "-seq2", "Sequence 2", true);
	setUp.setOption(kLength, "-kLength", "Kmer Length");
	setUp.setOption(getReverseDistance, "-getReverseDistance", "Get Reverse Complement Distance");
  setUp.finishSetUp(std::cout);

  seqWithKmerInfo seq1Obj(seqInfo("seq1", seq1), kLength, getReverseDistance);

  seqWithKmerInfo seq2Obj(seqInfo("seq2", seq2), kLength, getReverseDistance);
  std::pair<uint32_t,double> distance;
  if(getReverseDistance){
  	distance = seq1Obj.compareKmers(seq2Obj);
  }else{
  	distance = seq1Obj.compareKmers(seq2Obj);
  }

  std::cout << "KmersShared:" << distance.first << "("
  		<< std::min(seq1.size(), seq2.size()) - kLength + 1 << ")" << "\n";
  std::cout << "KmersDistance:" << distance.second << "\n";


  return 0;
}
class testRead : public seqWithKmerInfo {
public:
	testRead(const seqInfo & seqBase): seqWithKmerInfo(seqBase){

	}

	charCounter counter_;

	void setCount(const std::vector<char> & alph){
		counter_ = charCounter(alph);
		counter_.increaseCountByString(seqBase_.seq_, seqBase_.cnt_);
	}
};
/*
std::function<double(const std::unique_ptr<seqWithKmerInfo> & read1,
		const std::unique_ptr<seqWithKmerInfo> & read2,
		uint32_t windowSize, uint32_t windowStepSize)> disFun =
		[](const std::unique_ptr<seqWithKmerInfo> & read1,
				const std::unique_ptr<seqWithKmerInfo> & read2,
				uint32_t windowSize, uint32_t windowStepSize){
	auto scanDist = read1->slideCompareKmers(*read2, windowSize, windowStepSize);
	double averageDist = 0;
  for(const auto & distPos : iter::range(scanDist.size())){
  	averageDist += scanDist[distPos].second;
	}
	return averageDist/scanDist.size(); scoveViaKmers
};*/
int kmerExpRunner::scaningKmerDist(const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	uint32_t windowSize = 100;
	uint32_t windowStepSize = 10;
	uint32_t numThreads = 1;
	setUp.processDefaultReader(true);
	setUp.setOption(numThreads, "-numThreads", "numThreads");
	setUp.setOption(setUp.pars_.colOpts_.kmerOpts_.kLength_, "-k,-kLength", "kLength");
	setUp.setOption(windowSize, "-windowSize", "windowSize");
	setUp.setOption(windowStepSize, "-windowStepSize", "windowStepSize");
	//setUp.processDirectoryOutputName(true);
	setUp.processDebug();
	setUp.processSeq(true);
	setUp.finishSetUp(std::cout);
	//setUp.startARunLog(setUp.pars_.directoryName_);

	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<readObject>();
	auto reads = createKmerReadVec(inReads, setUp.pars_.colOpts_.kmerOpts_.kLength_, false);
	std::unique_ptr<seqWithKmerInfo> seq = std::make_unique<seqWithKmerInfo>(
			setUp.pars_.seqObj_.seqBase_, setUp.pars_.colOpts_.kmerOpts_.kLength_, false);
	seq->setKmers(setUp.pars_.colOpts_.kmerOpts_.kLength_, false);
	std::cout << "readName\taverageDist\tfullDist\n";
	for (const auto & read : reads) {
		auto scanDist = read->slideCompareKmers(*seq, windowSize, windowStepSize);
		double averageDist = 0;
		for (const auto & distPos : iter::range(scanDist.size())) {
			averageDist += scanDist[distPos].second;
		}
		averageDist /= scanDist.size();
		auto fullDist = seq->compareKmers(*read);
		std::cout << read->seqBase_.name_ << "\t" << averageDist << "\t"
				<< fullDist.second << "\n";
	}
	return 0;
}


int kmerExpRunner::kDistVsNucDist(const njh::progutils::CmdArgs & inputCommands){
	seqSetUp setUp(inputCommands);
  uint32_t kLength = 10;
  uint32_t numThreads = 1;
	setUp.processDefaultReader(true);
	setUp.setOption(kLength, "-k,--kLength", "Kmer Length");
	setUp.setOption(numThreads, "-t,--numThreads", "Number of Threads to Use");
  setUp.finishSetUp(std::cout);

  SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<readObject>();
  std::vector<testRead> reads;
  for(const auto & read : inReads){
  	reads.emplace_back(read.seqBase_);
  }
  njh::for_each(reads, [&](testRead & read){
  	read.setCount({'A','C','G','T'});
  	read.counter_.setFractions();
  	read.setKmers(kLength, false);});

  std::function<std::pair<double, double>(const testRead & read1,
  		const testRead & read2)> disFun =
  		[](const testRead & read1,
  				const testRead & read2){
  	auto dist = read1.compareKmers(read2);
  	auto nucDist = read1.counter_.getFracDifference(read2.counter_, read1.counter_.alphabet_);
  	return std::pair<double,double>{nucDist, dist.second};
  	//return dist.second;
  };
  auto distances = getDistance(reads, numThreads, disFun);
  std::ofstream outInfoFile("nucVsNucDist_" + bfs::basename(setUp.pars_.ioOptions_.firstName_));
  outInfoFile << "read1\tread2\tnucDist\tDist\n";
  for(const auto & pos : iter::range(distances.size())){
  	for(const auto & secondPos : iter::range(pos)){
  		outInfoFile<< reads[pos].seqBase_.name_
  				<< "\t" << reads[secondPos].seqBase_.name_
  				<< "\t" << distances[pos][secondPos].first
  				<< "\t" << distances[pos][secondPos].second
  				<< "\n";
  	}
  }
  return 0;
}

struct kRefCluster{

	kRefCluster(const seqInfo & seqBase): refSeq_(std::make_unique<seqWithKmerInfo>(seqBase)){}
	std::unique_ptr<seqWithKmerInfo> refSeq_;
	std::vector<std::unique_ptr<seqWithKmerInfo>> reads_;

};

int kmerExpRunner::clostestKmerDist(const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	setUp.processDefaultReader(true);
	setUp.processRefFilename(true);
	setUp.processDirectoryOutputName(true);
	setUp.processVerbose();
	setUp.setOption(setUp.pars_.colOpts_.kmerOpts_.kLength_, "-k,-kLength", "kLength");
  setUp.finishSetUp(std::cout);
  uint64_t maxLength = 0;
  SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<readObject>();
  auto inRefSeqs = SeqInput::getReferenceSeq(setUp.pars_.refIoOptions_, maxLength);
  readVec::getMaxLength(inReads, maxLength);
  std::vector<kRefCluster> refSeqs;
  for(const auto & ref : inRefSeqs){
  	refSeqs.emplace_back(kRefCluster(ref.seqBase_));
  	refSeqs.back().refSeq_->setKmers(setUp.pars_.colOpts_.kmerOpts_.kLength_, false);
  }
  std::vector<std::unique_ptr<seqWithKmerInfo>> reads;
  for(const auto & read : inReads){
  	reads.emplace_back(std::make_unique<seqWithKmerInfo>(read.seqBase_));
  }
  allSetKmers(reads, setUp.pars_.colOpts_.kmerOpts_.kLength_, false);
  std::vector<std::unique_ptr<seqWithKmerInfo>> rejectReads;
  for(const auto & readPos : iter::range(reads.size())){
  	if(readPos % 50 == 0){
  		std::cout << "Currently on " << readPos << " of " << reads.size() << "\n";
  	}
  	std::pair<uint32_t, double> best = {0,0};
  	uint32_t bestIndex = std::numeric_limits<uint32_t>::max();
  	for(const auto & refPos : iter::range(refSeqs.size())){
  		auto current = refSeqs[refPos].refSeq_->compareKmers(*reads[readPos]);
  		if(current.second > best.second){
  			best = current;
  			bestIndex = refPos;
  		}
  	}
  	if(bestIndex != std::numeric_limits<uint32_t>::max()){
  		refSeqs[bestIndex].reads_.emplace_back(std::forward<std::unique_ptr<seqWithKmerInfo>>(reads[readPos]));
  	}else{
  		rejectReads.emplace_back(std::forward<std::unique_ptr<seqWithKmerInfo>>(reads[readPos]));
  	}
  }

  for(const auto & ref : refSeqs){
  	std::ofstream currentReadFile;
  	openTextFile(currentReadFile, setUp.pars_.directoryName_ + ref.refSeq_->seqBase_.name_,
  			"fastq", setUp.pars_.ioOptions_.out_);
  	for(const auto & read : ref.reads_){
  		read->seqBase_.outPutFastq(currentReadFile);
  	}
  }
  if(!rejectReads.empty()){
  	std::ofstream currentReadFile;
  	openTextFile(currentReadFile, setUp.pars_.directoryName_ + "rejected",
  			"fastq", setUp.pars_.ioOptions_.out_);
  	for(const auto & read : rejectReads){
  		read->seqBase_.outPutFastq(currentReadFile);
  	}
  }
  setUp.logRunTime(std::cout);
  return 0;
}



int kmerExpRunner::profileScanningKmerDist(const njh::progutils::CmdArgs & inputCommands){
  seqSetUp setUp(inputCommands);
  uint32_t windowSize = 50;
  uint32_t windowStep = 10;
  std::string outFilename = "";
  setUp.setOption(outFilename, "-o,--outFile","name of the output file", true );
  setUp.setOption(windowSize, "--windowSize", "Size of scanning kmer window");
  setUp.setOption(windowStep, "--windowStep", "size of the scanning kmer dist step");
  std::string seq2 = "";
  setUp.processSeq(seq2, "-seq2", "Seq2", true);
  setUp.processSeq(true);
  setUp.setOption(setUp.pars_.colOpts_.kmerOpts_.kLength_, "-k,--kmerLength", "Length for kmers");
  setUp.processWritingOptions();
  bool append = false;
  setUp.setOption(append, "--append", "Append File");
  setUp.processVerbose();
  setUp.finishSetUp(std::cout);

  seqWithKmerInfo seq1Obj(setUp.pars_.seqObj_.seqBase_);
  seqWithKmerInfo seq2Obj(seqInfo("Seq2", seq2));
  seq1Obj.setKmers(setUp.pars_.colOpts_.kmerOpts_.kLength_, false);
  seq2Obj.setKmers(setUp.pars_.colOpts_.kmerOpts_.kLength_, false);

  std::ofstream outFile;
  njh::appendAsNeeded(outFilename, ".tab.txt");
  bool fileExists = bfs::exists(outFilename);
	njh::files::openTextFile(outFile, outFilename,
			setUp.pars_.ioOptions_.out_.overWriteFile_, append,
			setUp.pars_.ioOptions_.out_.exitOnFailureToWrite_);
  if(!fileExists || setUp.pars_.ioOptions_.out_.overWriteFile_){
  	outFile << "fullSeq\tseqParts\tcase\tkLen\twindowNumber\tkShared\tkDist" << std::endl;
  }
  auto seq1Dist = seq1Obj.slideCompareSubKmersToFull(seq2Obj, windowSize, windowStep);
  auto seq2Dist = seq2Obj.slideCompareSubKmersToFull(seq1Obj, windowSize, windowStep);
  for(const auto & dist : seq2Dist){
  	outFile << "seq1\tseq2\tseq1vsseq2\t" << setUp.pars_.colOpts_.kmerOpts_.kLength_
  			<< "\t" << dist.first
  			<< "\t" << dist.second.first
				<< "\t" << dist.second.second
				<< "\n";
  }
  for(const auto & dist : seq2Dist){
  	outFile << "seq2\tseq1\tseq2vsseq1\t" << setUp.pars_.colOpts_.kmerOpts_.kLength_
  			<< "\t" << dist.first
  			<< "\t" << dist.second.first
				<< "\t" << dist.second.second
				<< "\n";
  }
	return 0;
}

int kmerExpRunner::profileKmerAccerlation(const njh::progutils::CmdArgs & inputCommands){
  seqSetUp setUp(inputCommands);
  uint32_t numThreads = 1;
  setUp.processDefaultReader(true);
  //setUp.processRefFilename(true);
  setUp.setOption(numThreads, "-t,--numThreads", "Number of Threads to Use");
  setUp.setOption(setUp.pars_.colOpts_.kmerOpts_.kLength_, "-k,--kmerLength", "Length for kmers");
  setUp.pars_.ioOptions_.out_.outExtention_ = ".tab.txt";

  setUp.finishSetUp(std::cout);
  SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<readObject>();

  njh::stopWatch watch;
  watch.setLapName("index kmers");
  std::vector<std::unique_ptr<seqWithKmerInfo>> reads;
  for(const auto & read : inReads){
  	reads.emplace_back(std::make_unique<seqWithKmerInfo>(read.seqBase_));
  }

  uint64_t maxLength = 0;
  readVec::getMaxLength(inReads, maxLength);


  std::ofstream outFile;
  openTextFile(outFile, setUp.pars_.ioOptions_.out_);
  //outFile << "seq1\tseq2\tkLen\tkDist\tkShared\n";
  outFile << njh::conToStr(VecStr{"seq1", "seq2","group", "kLen", "seq1Group", "seq2Group", "case", "kDist", "kShared"}, "\t") << std::endl;
  /*std::function<std::pair<uint32_t, double>(const std::unique_ptr<seqWithKmerInfo> & read1,
  		const std::unique_ptr<seqWithKmerInfo> & read2)> disFun =
  		[](const std::unique_ptr<seqWithKmerInfo> & read1,
  				const std::unique_ptr<seqWithKmerInfo> & read2){
  	auto dist = read1->compareKmers(*read2); profilePacbioReads
  	return dist;
  };
  auto distances = getDistance(reads, numThreads, disFun);*/
  std::function<std::pair<uint32_t,double>(const std::unique_ptr<seqWithKmerInfo> & read1,
			const std::unique_ptr<seqWithKmerInfo> & read2)> disFun =
			[](const std::unique_ptr<seqWithKmerInfo> & read1,
					const std::unique_ptr<seqWithKmerInfo> & read2){
		auto dist = read1->compareKmers(*read2);
		return dist;
	};
  auto processGroup = [](const std::string & name){
  	auto underPos = name.find("_");
  	if(underPos != std::string::npos){
  		return name.substr(0, underPos);
  	}else{
  		return name.substr(0, name.find("."));
  	}
  };
  for(const auto & k : iter::range<uint32_t>(2, setUp.pars_.colOpts_.kmerOpts_.kLength_ + 1)){
  	std::cout << "Currently on Kmer Length " << k << '\r';
  	std::cout.flush();
  	allSetKmers(reads, k,false);
  	auto distances = getDistance(reads, numThreads, disFun);
  	for(const auto & pos : iter::range(distances.size())){
  		for(const auto & subPos : iter::range(distances[pos].size())){
  			auto seq1Name = reads[pos]->seqBase_.name_;
  			auto seq2Name = reads[subPos]->seqBase_.name_;
  			auto seq1Group = processGroup(seq1Name);
  			auto seq2Group = processGroup(seq2Name);
  			auto dist = distances[pos][subPos];
  			outFile << seq1Name
  					<< "\t" << seq2Name
						<< "\t" << seq1Name << "_v_" << seq2Name
						<< "\t" << k
						<< "\t" << seq1Group
						<< "\t" << seq2Group
						<< "\t" << (seq1Group == seq2Group ? "same" : "different")
						<< "\t" << dist.second
						<< "\t" << dist.first
						<< "\n";
  		}
  	}
  }
  std::cout << std::endl;
	return 0;
}



int kmerExpRunner::getNewScanKmerDist(const njh::progutils::CmdArgs & inputCommands){
  seqSetUp setUp(inputCommands);
  //uint32_t numThreads = 1; graphTest
  std::string filename = "";
  uint32_t windowSize = 50;
  uint32_t windowStep = 10;
  setUp.setOption(windowSize, "--windowSize", "Size of scanning kmer window");
  setUp.setOption(windowStep, "--windowStep", "size of the scanning kmer dist step");
  setUp.setOption(filename, "-o,--outFile", "Out Filename", true);
  setUp.processDefaultReader(true);
  setUp.processRefFilename(true);
  setUp.setOption(setUp.pars_.colOpts_.kmerOpts_.kLength_, "-k,--kmerLength", "Length for kmers");
  setUp.processVerbose();
  setUp.finishSetUp(std::cout);
  SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<readObject>();

  njh::stopWatch watch;
  watch.setLapName("index kmers");
  std::vector<std::unique_ptr<seqWithKmerInfo>> reads;
  for(const auto & read : inReads){
  	reads.emplace_back(std::make_unique<seqWithKmerInfo>(read.seqBase_));
  }

  uint64_t maxLength = 0;
  readVec::getMaxLength(inReads, maxLength);
  auto refSeqs = SeqInput::getReferenceSeq(setUp.pars_.refIoOptions_, maxLength);
  std::vector<std::unique_ptr<seqWithKmerInfo>> refReads;
  for(const auto & read : refSeqs){
  	refReads.emplace_back(std::make_unique<seqWithKmerInfo>(read.seqBase_));
  }

  std::ofstream outFile;
  openTextFile(outFile, filename, ".tab.txt", setUp.pars_.ioOptions_.out_);
  outFile << "seq\tref\twindowNum\tkDist\tkShared\n";
	allSetKmers(refReads, setUp.pars_.colOpts_.kmerOpts_.kLength_,false);
	allSetKmers(reads, setUp.pars_.colOpts_.kmerOpts_.kLength_,false);
	for(const auto & read : reads){
		for(const auto & ref : refReads){
			auto dists = ref->slideCompareSubKmersToFull(*read, windowSize, windowStep);
			for(const auto & dist : dists){
				outFile << read->seqBase_.name_
						<< "\t" << ref->seqBase_.name_
						<< "\t" << dist.first
						<< "\t" << dist.second.first
						<< "\t" << dist.second.second
						<< "\n";
			}
		}
	}
	if(setUp.pars_.verbose_){
		watch.logLapTimes(std::cout, true, 6, true);
	}

	return 0;
}





int kmerExpRunner::getKmerDist(const njh::progutils::CmdArgs & inputCommands){
  seqSetUp setUp(inputCommands);
  setUp.processDefaultReader(true);
  setUp.processRefFilename(true);
  setUp.processKmerLenOptions();
  setUp.finishSetUp(std::cout);
  SeqInput reader(setUp.pars_.ioOptions_);
	auto inReads = reader.readAllReads<readObject>();

  njh::stopWatch watch;
  watch.setLapName("index kmers");
  std::vector<std::unique_ptr<seqWithKmerInfo>> reads;
  for(const auto & read : inReads){
  	reads.emplace_back(std::make_unique<seqWithKmerInfo>(read.seqBase_));
  }

  uint64_t maxLength = 0;
  readVec::getMaxLength(inReads, maxLength);
  auto refSeqs = SeqInput::getReferenceSeq(setUp.pars_.refIoOptions_, maxLength);
  std::vector<std::unique_ptr<seqWithKmerInfo>> refReads;
  for(const auto & read : refSeqs){
  	refReads.emplace_back(std::make_unique<seqWithKmerInfo>(read.seqBase_));
  }

  std::ofstream outFile;
  openTextFile(outFile, setUp.pars_.ioOptions_.out_.outFilename_, ".tab.txt", setUp.pars_.ioOptions_.out_);
  outFile << "seq\tref\tkLen\tkDist\tkShared\n";
  for(const auto & k : iter::range<uint32_t>(2, setUp.pars_.colOpts_.kmerOpts_.kLength_ + 1)){
  	std::cout << "Currently on Kmer Length " << k << '\r';
  	std::cout.flush();
  	allSetKmers(refReads, k,false);
  	allSetKmers(reads, k,false);
  	for(const auto & read : reads){
  		for(const auto & ref : refReads){
  			auto dist = ref->compareKmers(*read);
  			outFile << read->seqBase_.name_
  					<< "\t" << ref->seqBase_.name_
						<< "\t" << k
						<< "\t" << dist.second
  					<< "\t" << dist.first << "\n";
  		}
  	}
  }
  std::cout << std::endl;
	return 0;
}


int kmerExpRunner::readingDistanceCheck(const njh::progutils::CmdArgs & inputCommands){
  seqSetUp setUp(inputCommands);
  uint32_t numThreads = 1;
  std::string filename = "";
  setUp.processDefaultReader(true);
  setUp.setOption(numThreads, "-t,--numThreads", "Number of Threads to Use");
  setUp.setOption(setUp.pars_.colOpts_.kmerOpts_.kLength_, "-k,--kmerLength", "Length for kmers");
  setUp.setOption(filename, "-i,--inFile", "In Filename", true);
  setUp.finishSetUp(std::cout);
  njh::stopWatch watch;
  watch.setLapName("index kmers");
  auto reads = createKmerReadVec(setUp.pars_.ioOptions_,setUp.pars_.colOpts_.kmerOpts_.kLength_,false );
  std::function<double(const std::unique_ptr<seqWithKmerInfo> & read1,
  		const std::unique_ptr<seqWithKmerInfo> & read2)> disFun =
  		[](const std::unique_ptr<seqWithKmerInfo> & read1,
  				const std::unique_ptr<seqWithKmerInfo> & read2){
  	auto dist = read1->compareKmers(*read2);
  	return dist.second;
  };
  watch.startNewLap("regular distance");
  auto distances = getDistance(reads, numThreads, disFun);
  std::ifstream inFile(filename);
  if(!inFile){
  	std::cerr << njh::bashCT::red <<
  			njh::bashCT::bold << "Error in opening " << filename
				<< njh::bashCT::reset << std::endl;
  	exit(1);
  }
  std::vector<std::vector<double>> inDist = readDistanceMatrix(inFile);
  bool match = true;
  for(const auto & i : iter::range(inDist.size())){
  	for(const auto & j : iter::range(inDist[i].size())){
  		if(roundDecPlaces(inDist[i][j], 4) != roundDecPlaces(distances[i][j], 4)){
  			std::cout << "inDist: " << inDist[i][j] << ":" << roundDecPlaces(inDist[i][j], 4)  << std::endl;
  			std::cout << "outDist: " << distances[i][j] << ":" << roundDecPlaces(distances[i][j], 4)  << std::endl;
  			match =false;
  			break;
  		}
  	}
  }
  if(match){
  	std::cout << njh::bashCT::green <<
  			njh::bashCT::bold << "match"<< filename
				<< njh::bashCT::reset << std::endl;
  }else{
  	std::cout << njh::bashCT::red <<
  			njh::bashCT::bold << "Don't match" << filename
				<< njh::bashCT::reset << std::endl;
  }

	return 0;
}


int kmerExpRunner::writingDistanceCheck(const njh::progutils::CmdArgs & inputCommands){
  seqSetUp setUp(inputCommands);
  uint32_t numThreads = 1;
  setUp.processDefaultReader(false);
  setUp.setOption(numThreads, "-t,--numThreads", "Number of Threads to Use");
  setUp.setOption(setUp.pars_.colOpts_.kmerOpts_.kLength_, "-k,--kmerLength", "Length for kmers");
  setUp.finishSetUp(std::cout);

	std::vector<readObject> inReads;
  if(setUp.pars_.ioOptions_.firstName_ == ""){
    njh::randomGenerator gen;
    VecStr rSeqs = simulation::evenRandStrs(40, std::vector<char>{'A', 'C', 'G', 'T'}, gen, 20);
    inReads = vecStrToReadObjs<readObject>(rSeqs, "Seq");
  }else{
  	SeqInput reader(setUp.pars_.ioOptions_);
  	reader.openIn();
  	inReads = reader.readAllReads<readObject>();
  }

  njh::stopWatch watch;
  watch.setLapName("index kmers");
  std::vector<std::unique_ptr<seqWithKmerInfo>> reads;
  for(const auto & read : inReads){
  	reads.emplace_back(std::make_unique<seqWithKmerInfo>(read.seqBase_));
  }
  allSetKmers(reads, setUp.pars_.colOpts_.kmerOpts_.kLength_,false);
  std::function<double(const std::unique_ptr<seqWithKmerInfo> & read1,
  		const std::unique_ptr<seqWithKmerInfo> & read2)> disFun =
  		[](const std::unique_ptr<seqWithKmerInfo> & read1,
  				const std::unique_ptr<seqWithKmerInfo> & read2){
  	auto dist = read1->compareKmers(*read2);
  	return dist.second;
  };
  watch.startNewLap("regular distance");
  auto distances = getDistance(reads, numThreads, disFun);
  std::ofstream outFile;
  openTextFile(outFile, setUp.pars_.ioOptions_.out_.outFilename_.string(), ".tab.txt", setUp.pars_.ioOptions_.out_);
  writeDistanceMatrix(outFile, distances);
  if(setUp.pars_.ioOptions_.firstName_ == ""){
    std::ofstream outSeqFile;
    openTextFile(outSeqFile, setUp.pars_.ioOptions_.out_.outFilename_.string() + ".fasta", ".fasta", setUp.pars_.ioOptions_.out_);
    for(const auto & read : reads){
    	read->seqBase_.outPutSeq(outSeqFile);
    }
  }
	return 0;
}
//scoveViaKmers

int kmerExpRunner::getKmerDistStats(
		const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	TableIOOpts tabOpts = TableIOOpts::genTabFileOut("", true);
	setUp.processWritingOptions(tabOpts.out_);
	setUp.processReadInNames(setUp.readInFormatsAvailable_);
	setUp.processSeq(true);
	setUp.processKmerLenOptions();
	setUp.finishSetUp(std::cout);
	seqWithKmerInfo compare(setUp.pars_.seqObj_.seqBase_, setUp.pars_.colOpts_.kmerOpts_.kLength_,
			false);
	auto outStatsTab = getKmerStatsOnFile(setUp.pars_.ioOptions_, compare);
	outStatsTab.outPutContents(tabOpts);
	return 0;
}

int kmerExpRunner::getKmerDistStatsMultiple(
		const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	std::string multipleFiles = "";
	TableIOOpts tabOpts = TableIOOpts::genTabFileOut("", true);
	setUp.processWritingOptions(tabOpts.out_);
	setUp.setOption(multipleFiles, "-files", "multipleFiles");
	setUp.processVerbose();
	setUp.processSeq(true);
	setUp.processKmerLenOptions();
	setUp.finishSetUp(std::cout);

	seqWithKmerInfo compare(setUp.pars_.seqObj_.seqBase_, setUp.pars_.colOpts_.kmerOpts_.kLength_,
			false);
	auto toks = tokenizeString(multipleFiles, ",");
	table outStats;
	for (const auto & file : toks) {
		SeqIOOptions seqOpts;
		seqOpts.firstName_ = file;
		seqOpts.inFormat_ = SeqIOOptions::getInFormat(njh::files::getExtension(file));
		if (outStats.content_.empty()) {
			outStats = getKmerStatsOnFile(seqOpts, compare);
		} else {
			outStats.rbind(getKmerStatsOnFile(seqOpts, compare), false);
		}
	}
	outStats.outPutContents(tabOpts);
	return 0;
}

} /* namespace njhseq */
