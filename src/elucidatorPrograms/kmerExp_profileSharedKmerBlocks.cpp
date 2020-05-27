/*
 * kmerExp_profileSharedKmerBlocks.cpp
 *
 *  Created on: May 25, 2020
 *      Author: nick
 */




#include "kmerExp.hpp"
#include "elucidator/utils/KmerUtils.hpp"
#include "elucidator/objects/seqObjects/seqKmers.h"

namespace njhseq {
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

	std::ofstream perfectCompResultsFile;
	openTextFile(perfectCompResultsFile, setUp.pars_.directoryName_ + "perfect", ".fasta",
			setUp.pars_.ioOptions_.out_);
	std::ofstream perfectCompBedFile;
	openTextFile(perfectCompBedFile, setUp.pars_.directoryName_ + "perfect", ".bed",
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
	kmerInfo refGenInfo(setUp.pars_.seqObj_.seqBase_.seq_, kLen, false);

	std::unordered_map<std::string, KmersSharedBlocks> seqs;
	std::unordered_map<std::string, KmersSharedBlocks> perfectSeqs;
	std::unordered_map<std::string, KmersSharedBlocks> nearPerfectSeqs;
	uint64_t maxPosition = len(setUp.pars_.seqObj_);

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
			readVec::getMaxLength(inReads.front(), maxPosition);
		}
	}else{
		SeqInput reader(setUp.pars_.ioOptions_);
		reader.openIn();
		auto inReads = reader.readAllReads<readObject>();
		readVec::getMaxLength(inReads.front(), maxPosition);
		for(const auto & seq : inReads){
			seqs.emplace(seq.seqBase_.name_,
					KmersSharedBlocks(seq.seqBase_.name_,
							seq.seqBase_.seq_, kLen));
		}
	}


	perfectCompResultsHtmlFile << njh::replaceString(perfectKmerHtmlIndexFileStr, "MAXPOSITION", estd::to_string(maxPosition));
	perfectSeqs = seqs;
	nearPerfectSeqs = seqs;
	watch.startNewLap("Search All");
	for (const auto & pos : iter::range(setUp.pars_.seqObj_.seqBase_.seq_.size() - kLen + 1)) {
		auto sub = setUp.pars_.seqObj_.seqBase_.seq_.substr(pos, kLen);
		for (auto & seq : seqs) {
			auto search = seq.second.kInfo_->kmers_.find(sub);
			//if both sequences contain the kmer and the kmer is unique in both sequences, add them
			if (search != seq.second.kInfo_->kmers_.end()
					&& search->second.positions_.size() == 1
					&& refGenInfo.kmers_[search->first].positions_.size() == 1) {
				seq.second.addComp(refGenInfo.kmers_[search->first].positions_.front(), search->second.positions_.front());
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
	uint64_t maxMinStartAll = 0;
	for (const auto &g : perfectSeqs.begin()->second.kComps_){
		auto refStop = g.second.size_ + g.second.refStart_ + kLen - 1;
		perfectCompBedFile << setUp.pars_.seqObj_.seqBase_.name_
				<< "\t" << g.second.refStart_
				<< "\t" << refStop
				<< "\t" << njh::pasteAsStr(setUp.pars_.seqObj_.seqBase_.name_, "-", g.second.refStart_, "-", refStop )
				<< "\t" << g.second.size_ + kLen - 1
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

	watch.logLapTimes(std::cout, true, 6, true);
	return 0;
}

} // namespace njhseq

