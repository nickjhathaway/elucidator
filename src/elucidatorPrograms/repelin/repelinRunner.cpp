/*
 * repelinRunner.cpp
 *
 *  Created on: Jan 11, 2015
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
#include "repelinRunner.hpp"
#include <TwoBit.h>
#include <seqServer.h>
#include "elucidator/objects/BioDataObject.h"

#include <njhseq.h>



namespace njhseq {
namespace bfs = boost::filesystem;

repelinRunner::repelinRunner()
    : njh::progutils::ProgramRunner(

    		{	 addFunc("parseRepeatMaskerOutputUnitTest", parseRepeatMaskerOutputUnitTest, false),
					 addFunc("parseRepeatMaskerOutputForElement", parseRepeatMaskerOutputForElement, false),
					 addFunc("extractFullElementSequences", extractFullElementSequences, false),
					 addFunc("extractElementSequences", extractElementSequences, false),
					 addFunc("parseRMForElementWithoutIntervening", parseRMForElementWithoutIntervening, false),
					 addFunc("parseTandemRepeatFinderOutputUnitTest", parseTandemRepeatFinderOutputUnitTest, false),
					 addFunc("TandemRepeatFinderOutputToBed", TandemRepeatFinderOutputToBed, false),
					 addFunc("runTRF", runTRF, false)
           },//
          "repelin") {}

std::unordered_map<uint32_t, std::vector<RepeatMaskerRecord>> readRepeatMaskerFile(const std::string & filename, bool verbose){
	std::ifstream inFile(filename);
	if (!inFile) {
		std::stringstream ss;
		ss << "Error in opening " << njh::bashCT::boldRed(filename) << std::endl;
		throw std::runtime_error{ss.str()};
	}
	std::unordered_map<uint32_t, std::vector<RepeatMaskerRecord>> ret;
	std::string line;
	uint32_t lineNumber = 0;
	if(verbose){
		std::cout << "Parsing file: " << filename << std::endl;
	}
	while (njh::files::crossPlatGetline(inFile, line)) {
		if (njh::beginsWith(line, "   SW")) {
			continue;
		}
		if (njh::beginsWith(line, "score")) {
			continue;
		}
		if(line.empty() || allWhiteSpaceStr(line)){
			continue;
		}
		if (verbose && lineNumber % 1000 == 0) {
			std::cout << "\r" << "On " << lineNumber;
			std::cout.flush();
		}
		try {
			auto element = RepeatMaskerRecord(line);
			ret[element.regionSegment_].emplace_back(element);
		} catch (std::exception & e) {
			std::stringstream ss;
			ss << std::endl;
			ss << e.what() << std::endl;
			ss << "With parsing line: " << lineNumber << std::endl;
			ss << line << std::endl;
			throw std::runtime_error{ss.str()};
		}
		++lineNumber;
	}
	if(verbose){
		std::cout << std::endl;
	}
	return ret;
}

std::unordered_map<uint32_t, std::vector<std::shared_ptr<RepeatMaskerRecord>>> readRepeatMaskerFilePtrs(const std::string & filename, bool verbose){
	std::ifstream inFile(filename);
	if (!inFile) {
		std::stringstream ss;
		ss << "Error in opening " << njh::bashCT::boldRed(filename) << std::endl;
		throw std::runtime_error{ss.str()};
	}
	std::unordered_map<uint32_t, std::vector<std::shared_ptr<RepeatMaskerRecord>>> ret;
	std::string line;
	uint32_t lineNumber = 0;
	if(verbose){
		std::cout << "Parsing file: " << filename << std::endl;
	}
	while (njh::files::crossPlatGetline(inFile, line)) {
		if (njh::beginsWith(line, "   SW")) {
			continue;
		}
		if (njh::beginsWith(line, "score")) {
			continue;
		}
		if(line.empty() || allWhiteSpaceStr(line)){
			continue;
		}
		if (verbose && lineNumber % 1000 == 0) {
			std::cout << "\r" << "On " << lineNumber;
			std::cout.flush();
		}
		try {
			auto element = std::make_shared<RepeatMaskerRecord>(line);
			ret[element->regionSegment_].emplace_back(element);
		} catch (std::exception & e) {
			std::stringstream ss;
			ss << std::endl;
			ss << e.what() << std::endl;
			ss << __PRETTY_FUNCTION__ << std::endl;
			ss << "With parsing line: " << lineNumber << std::endl;
			ss << line << std::endl;
			throw std::runtime_error{ss.str()};
		}
		++lineNumber;
	}
	if(verbose){
		std::cout << std::endl;
	}
	return ret;
}

int repelinRunner::parseRepeatMaskerOutputUnitTest(const njh::progutils::CmdArgs & inputCommands){
  seqSetUp setUp(inputCommands);
  std::string filename = "";
  setUp.setOption(filename, "--file,-f", "Name of repeat masker output file", true);
  setUp.processVerbose();
  setUp.finishSetUp(std::cout);

  auto elements = readRepeatMaskerFile(filename, setUp.pars_.verbose_);

  return 0;
}

int repelinRunner::parseRepeatMaskerOutputForElement(const njh::progutils::CmdArgs & inputCommands){
  seqSetUp setUp(inputCommands);
  std::string filename = "";
  TableIOOpts outOpts = TableIOOpts::genTabFileOut("out", false);
  std::string elementName = "";
  setUp.setOption(filename, "-i,--inFile", "Name of repeat masker output file", true);
  setUp.setOption(elementName, "-e,--elementName", "Name of the element to extract for in regex format", true);
  setUp.processWritingOptions(outOpts.out_);
  setUp.processVerbose();
  setUp.finishSetUp(std::cout);
  std::regex pat{elementName};
  auto elements = readRepeatMaskerFile(filename, setUp.pars_.verbose_);
  if (setUp.pars_.verbose_){
  	std::cout << "Parsing elements for elements containing: " << elementName << std::endl;
  }
  std::unordered_map<uint32_t, std::vector<RepeatMaskerRecord>> ret;
	for (const auto regionEn : iter::enumerate(elements)) {
		if (setUp.pars_.verbose_ && regionEn.index % 1000 == 0) {
			std::cout << "\r" << "On " << regionEn.index << " of " << elements.size() ;
			std::cout.flush();
		}
		for (const auto & element : regionEn.element.second) {
			if (std::regex_match(element.nameOfMatchedSeq_, pat)) {
				ret.emplace(regionEn.element);
				break;
			}
		}
	}
	if (setUp.pars_.verbose_){
		std::cout << std::endl;
	}
  std::ofstream outFile;
  openTextFile(outFile,  outOpts.out_);
  std::vector<uint32_t> keys = getVectorOfMapKeys(ret);
  njh::sort(keys);
  for(const auto & k: keys){
  	const auto & region = ret[k];
  	for(const auto & element : region){
  		outFile << element.getDelimitedInfoStr("\t") << std::endl;
  	}
  }
  return 0;
}

int repelinRunner::parseRMForElementWithoutIntervening(const njh::progutils::CmdArgs & inputCommands){
  seqSetUp setUp(inputCommands);
  std::string filename = "";
  std::string outFilename = "";
  std::string elementName = "";
  std::string twoBitFilename = "";
  setUp.processVerbose();
  setUp.setOption(filename, "-i,--inFile", "Name of repeat masker output file", true);
  setUp.setOption(outFilename, "-o,--outFile", "Name of the outfile to print to", true);
  setUp.setOption(elementName, "-e,--elementName", "Name of the element to extract for in regex format", true);
  setUp.setOption(twoBitFilename, "--twoBitFile", "Filename of two bit file", true);
  setUp.processWritingOptions();
  setUp.finishSetUp(std::cout);

  if (!bfs::exists(twoBitFilename)) {
  	std::cout << "Filename: " << twoBitFilename << " doesn't exist"
  			<< std::endl;
  	exit(1);
  }
  TwoBit::TwoBitFile f(twoBitFilename);
  std::string buffer;

  std::ifstream inFile(filename);
  if(!inFile){
  	std::cerr << "Error in opening " << filename << std::endl;
  	exit(1);
  }
  std::regex pat{elementName};
  std::string line;
  auto elements = readRepeatMaskerFilePtrs(filename, setUp.pars_.verbose_);
  std::unordered_map<std::string, std::unordered_map<uint32_t,std::shared_ptr<RepeatMaskerRecord>>> elementsByPositions;
  for(const auto & rmUnit : elements){
  	for(const auto & subUnits : rmUnit.second){
  		elementsByPositions[subUnits->nameOfQuery_][subUnits->start_] = subUnits;
  	}
  }
  std::cout << std::endl;
  std::cout << "Parsing elements for elements containing: " << elementName << std::endl;
  std::unordered_map<uint32_t, std::vector<std::shared_ptr<RepeatMaskerRecord>>> ret;
  for(const auto regionEn : iter::enumerate(elements)){
  	if(regionEn.index % 100  == 0){
  		std::cout << "On " << regionEn.index << " of " << elements.size()<< "\r";
  		std::cout.flush();
  	}
  	for(const auto & element : regionEn.element.second){
  		if(std::regex_match(element->nameOfMatchedSeq_,pat)){
  			ret.emplace(regionEn.element);
  			break;
  		}
  	}
  }
  std::cout << std::endl;
  std::ofstream outFile;
  std::ofstream outPlusFile;
  openTextFile(outFile, outFilename, ".tab.txt", setUp.pars_.ioOptions_.out_);
  openTextFile(outPlusFile, outFilename, "Plus.tab.txt", setUp.pars_.ioOptions_.out_);
  std::vector<uint32_t> keys = getVectorOfMapKeys(ret);
  njh::sort(keys);
  uint32_t kCount = 0;
  for(const auto & k: keys){
  	++kCount;
  	if(kCount % 100  == 0){
  		std::cout << "On " << kCount << " of " << keys.size()<< "\r";
  		std::cout.flush();
  	}
  	const auto & region = ret[k];
  	uint32_t minPos = std::numeric_limits<uint32_t>::max();
  	uint32_t maxPos = 0;
  	for(const auto & element : region){
  		outFile << element->getDelimitedInfoStr("\t") << std::endl;
  		if(element->start_ < minPos){
  			minPos = element->start_;
  		}
  		if(element->end_ > maxPos){
  			maxPos = element->end_;
  		}
  	}
		std::unordered_map<uint32_t,
				std::vector<std::shared_ptr<RepeatMaskerRecord>>>inBetween;
		for (const auto pos : iter::range(minPos, maxPos)) {
			auto search = elementsByPositions[region.front()->nameOfQuery_].find(pos);
			if (search != elementsByPositions[region.front()->nameOfQuery_].end()
					&& search->second->regionSegment_ != region.front()->regionSegment_) {
				inBetween.emplace(search->second->regionSegment_,
						elements[search->second->regionSegment_]);
			}
		}
  	if(!inBetween.empty()){
    	for(const auto & element : region){
    		outPlusFile << element->getDelimitedInfoStr("\t") << std::endl;
    	}
    	for(const auto & ib : inBetween){
    		for(const auto & ele : ib.second){
    			outPlusFile << ele->getDelimitedInfoStr("\t") << std::endl;
    		}
    	}
  	}
  }
  std::cout << std::endl;
  return 0;
}

void checkFileExistsThrow(bfs::path file) {
	if (!bfs::exists(file)) {
		std::stringstream ss;
		ss << "Filename: " << njh::bashCT::boldRed(file.string())
				<< " doesn't exist\n";
		throw std::runtime_error { ss.str() };
	}
}

void genElementsHtml(std::ostream & out){
  out << "<!DOCTYPE html>" << std::endl;
  out << "<html>" << std::endl;
  out << "<head>" << std::endl;
  out << "<meta charset=\"utf-8\">" << std::endl;
  out << "<title>Repeat elements</title>" << std::endl;
  out << "<style>" << std::endl;
  out << "" << std::endl;
  out << "body {" << std::endl;
  out << "  font: 10px sans-serif;" << std::endl;
  out << "}" << std::endl;
  out << "" << std::endl;
  out << ".axis path," << std::endl;
  out << ".axis line {" << std::endl;
  out << "  fill: none;" << std::endl;
  out << "  stroke: #000;" << std::endl;
  out << "  shape-rendering: crispEdges;" << std::endl;
  out << "}" << std::endl;
  out << "" << std::endl;
  out << "</style>" << std::endl;
  out << "" << std::endl;
  out << "<script src=\"http://cdnjs.cloudflare.com/ajax/libs/d3/3.4.13/d3.min.js\"></script>" << std::endl;
  out << "<script>" << std::endl;
  out << "var colorbrewer={YlGn:{3:[\"#f7fcb9\",\"#addd8e\",\"#31a354\"],4:[\"#ffffcc\",\"#c2e699\",\"#78c679\",\"#238443\"],5:[\"#ffffcc\",\"#c2e699\",\"#78c679\",\"#31a354\",\"#006837\"],6:[\"#ffffcc\",\"#d9f0a3\",\"#addd8e\",\"#78c679\",\"#31a354\",\"#006837\"],7:[\"#ffffcc\",\"#d9f0a3\",\"#addd8e\",\"#78c679\",\"#41ab5d\",\"#238443\",\"#005a32\"],8:[\"#ffffe5\",\"#f7fcb9\",\"#d9f0a3\",\"#addd8e\",\"#78c679\",\"#41ab5d\",\"#238443\",\"#005a32\"],9:[\"#ffffe5\",\"#f7fcb9\",\"#d9f0a3\",\"#addd8e\",\"#78c679\",\"#41ab5d\",\"#238443\",\"#006837\",\"#004529\"]},YlGnBu:{3:[\"#edf8b1\",\"#7fcdbb\",\"#2c7fb8\"],4:[\"#ffffcc\",\"#a1dab4\",\"#41b6c4\",\"#225ea8\"],5:[\"#ffffcc\",\"#a1dab4\",\"#41b6c4\",\"#2c7fb8\",\"#253494\"],6:[\"#ffffcc\",\"#c7e9b4\",\"#7fcdbb\",\"#41b6c4\",\"#2c7fb8\",\"#253494\"],7:[\"#ffffcc\",\"#c7e9b4\",\"#7fcdbb\",\"#41b6c4\",\"#1d91c0\",\"#225ea8\",\"#0c2c84\"],8:[\"#ffffd9\",\"#edf8b1\",\"#c7e9b4\",\"#7fcdbb\",\"#41b6c4\",\"#1d91c0\",\"#225ea8\",\"#0c2c84\"],9:[\"#ffffd9\",\"#edf8b1\",\"#c7e9b4\",\"#7fcdbb\",\"#41b6c4\",\"#1d91c0\",\"#225ea8\",\"#253494\",\"#081d58\"]},GnBu:{3:[\"#e0f3db\",\"#a8ddb5\",\"#43a2ca\"],4:[\"#f0f9e8\",\"#bae4bc\",\"#7bccc4\",\"#2b8cbe\"],5:[\"#f0f9e8\",\"#bae4bc\",\"#7bccc4\",\"#43a2ca\",\"#0868ac\"],6:[\"#f0f9e8\",\"#ccebc5\",\"#a8ddb5\",\"#7bccc4\",\"#43a2ca\",\"#0868ac\"],7:[\"#f0f9e8\",\"#ccebc5\",\"#a8ddb5\",\"#7bccc4\",\"#4eb3d3\",\"#2b8cbe\",\"#08589e\"],8:[\"#f7fcf0\",\"#e0f3db\",\"#ccebc5\",\"#a8ddb5\",\"#7bccc4\",\"#4eb3d3\",\"#2b8cbe\",\"#08589e\"],9:[\"#f7fcf0\",\"#e0f3db\",\"#ccebc5\",\"#a8ddb5\",\"#7bccc4\",\"#4eb3d3\",\"#2b8cbe\",\"#0868ac\",\"#084081\"]},BuGn:{3:[\"#e5f5f9\",\"#99d8c9\",\"#2ca25f\"],4:[\"#edf8fb\",\"#b2e2e2\",\"#66c2a4\",\"#238b45\"],5:[\"#edf8fb\",\"#b2e2e2\",\"#66c2a4\",\"#2ca25f\",\"#006d2c\"],6:[\"#edf8fb\",\"#ccece6\",\"#99d8c9\",\"#66c2a4\",\"#2ca25f\",\"#006d2c\"],7:[\"#edf8fb\",\"#ccece6\",\"#99d8c9\",\"#66c2a4\",\"#41ae76\",\"#238b45\",\"#005824\"],8:[\"#f7fcfd\",\"#e5f5f9\",\"#ccece6\",\"#99d8c9\",\"#66c2a4\",\"#41ae76\",\"#238b45\",\"#005824\"],9:[\"#f7fcfd\",\"#e5f5f9\",\"#ccece6\",\"#99d8c9\",\"#66c2a4\",\"#41ae76\",\"#238b45\",\"#006d2c\",\"#00441b\"]},PuBuGn:{3:[\"#ece2f0\",\"#a6bddb\",\"#1c9099\"],4:[\"#f6eff7\",\"#bdc9e1\",\"#67a9cf\",\"#02818a\"],5:[\"#f6eff7\",\"#bdc9e1\",\"#67a9cf\",\"#1c9099\",\"#016c59\"],6:[\"#f6eff7\",\"#d0d1e6\",\"#a6bddb\",\"#67a9cf\",\"#1c9099\",\"#016c59\"],7:[\"#f6eff7\",\"#d0d1e6\",\"#a6bddb\",\"#67a9cf\",\"#3690c0\",\"#02818a\",\"#016450\"],8:[\"#fff7fb\",\"#ece2f0\",\"#d0d1e6\",\"#a6bddb\",\"#67a9cf\",\"#3690c0\",\"#02818a\",\"#016450\"],9:[\"#fff7fb\",\"#ece2f0\",\"#d0d1e6\",\"#a6bddb\",\"#67a9cf\",\"#3690c0\",\"#02818a\",\"#016c59\",\"#014636\"]},PuBu:{3:[\"#ece7f2\",\"#a6bddb\",\"#2b8cbe\"],4:[\"#f1eef6\",\"#bdc9e1\",\"#74a9cf\",\"#0570b0\"],5:[\"#f1eef6\",\"#bdc9e1\",\"#74a9cf\",\"#2b8cbe\",\"#045a8d\"],6:[\"#f1eef6\",\"#d0d1e6\",\"#a6bddb\",\"#74a9cf\",\"#2b8cbe\",\"#045a8d\"],7:[\"#f1eef6\",\"#d0d1e6\",\"#a6bddb\",\"#74a9cf\",\"#3690c0\",\"#0570b0\",\"#034e7b\"],8:[\"#fff7fb\",\"#ece7f2\",\"#d0d1e6\",\"#a6bddb\",\"#74a9cf\",\"#3690c0\",\"#0570b0\",\"#034e7b\"],9:[\"#fff7fb\",\"#ece7f2\",\"#d0d1e6\",\"#a6bddb\",\"#74a9cf\",\"#3690c0\",\"#0570b0\",\"#045a8d\",\"#023858\"]},BuPu:{3:[\"#e0ecf4\",\"#9ebcda\",\"#8856a7\"],4:[\"#edf8fb\",\"#b3cde3\",\"#8c96c6\",\"#88419d\"],5:[\"#edf8fb\",\"#b3cde3\",\"#8c96c6\",\"#8856a7\",\"#810f7c\"],6:[\"#edf8fb\",\"#bfd3e6\",\"#9ebcda\",\"#8c96c6\",\"#8856a7\",\"#810f7c\"],7:[\"#edf8fb\",\"#bfd3e6\",\"#9ebcda\",\"#8c96c6\",\"#8c6bb1\",\"#88419d\",\"#6e016b\"],8:[\"#f7fcfd\",\"#e0ecf4\",\"#bfd3e6\",\"#9ebcda\",\"#8c96c6\",\"#8c6bb1\",\"#88419d\",\"#6e016b\"],9:[\"#f7fcfd\",\"#e0ecf4\",\"#bfd3e6\",\"#9ebcda\",\"#8c96c6\",\"#8c6bb1\",\"#88419d\",\"#810f7c\",\"#4d004b\"]},RdPu:{3:[\"#fde0dd\",\"#fa9fb5\",\"#c51b8a\"],4:[\"#feebe2\",\"#fbb4b9\",\"#f768a1\",\"#ae017e\"],5:[\"#feebe2\",\"#fbb4b9\",\"#f768a1\",\"#c51b8a\",\"#7a0177\"],6:[\"#feebe2\",\"#fcc5c0\",\"#fa9fb5\",\"#f768a1\",\"#c51b8a\",\"#7a0177\"],7:[\"#feebe2\",\"#fcc5c0\",\"#fa9fb5\",\"#f768a1\",\"#dd3497\",\"#ae017e\",\"#7a0177\"],8:[\"#fff7f3\",\"#fde0dd\",\"#fcc5c0\",\"#fa9fb5\",\"#f768a1\",\"#dd3497\",\"#ae017e\",\"#7a0177\"],9:[\"#fff7f3\",\"#fde0dd\",\"#fcc5c0\",\"#fa9fb5\",\"#f768a1\",\"#dd3497\",\"#ae017e\",\"#7a0177\",\"#49006a\"]},PuRd:{3:[\"#e7e1ef\",\"#c994c7\",\"#dd1c77\"],4:[\"#f1eef6\",\"#d7b5d8\",\"#df65b0\",\"#ce1256\"],5:[\"#f1eef6\",\"#d7b5d8\",\"#df65b0\",\"#dd1c77\",\"#980043\"],6:[\"#f1eef6\",\"#d4b9da\",\"#c994c7\",\"#df65b0\",\"#dd1c77\",\"#980043\"],7:[\"#f1eef6\",\"#d4b9da\",\"#c994c7\",\"#df65b0\",\"#e7298a\",\"#ce1256\",\"#91003f\"],8:[\"#f7f4f9\",\"#e7e1ef\",\"#d4b9da\",\"#c994c7\",\"#df65b0\",\"#e7298a\",\"#ce1256\",\"#91003f\"],9:[\"#f7f4f9\",\"#e7e1ef\",\"#d4b9da\",\"#c994c7\",\"#df65b0\",\"#e7298a\",\"#ce1256\",\"#980043\",\"#67001f\"]},OrRd:{3:[\"#fee8c8\",\"#fdbb84\",\"#e34a33\"],4:[\"#fef0d9\",\"#fdcc8a\",\"#fc8d59\",\"#d7301f\"],5:[\"#fef0d9\",\"#fdcc8a\",\"#fc8d59\",\"#e34a33\",\"#b30000\"],6:[\"#fef0d9\",\"#fdd49e\",\"#fdbb84\",\"#fc8d59\",\"#e34a33\",\"#b30000\"],7:[\"#fef0d9\",\"#fdd49e\",\"#fdbb84\",\"#fc8d59\",\"#ef6548\",\"#d7301f\",\"#990000\"],8:[\"#fff7ec\",\"#fee8c8\",\"#fdd49e\",\"#fdbb84\",\"#fc8d59\",\"#ef6548\",\"#d7301f\",\"#990000\"],9:[\"#fff7ec\",\"#fee8c8\",\"#fdd49e\",\"#fdbb84\",\"#fc8d59\",\"#ef6548\",\"#d7301f\",\"#b30000\",\"#7f0000\"]},YlOrRd:{3:[\"#ffeda0\",\"#feb24c\",\"#f03b20\"],4:[\"#ffffb2\",\"#fecc5c\",\"#fd8d3c\",\"#e31a1c\"],5:[\"#ffffb2\",\"#fecc5c\",\"#fd8d3c\",\"#f03b20\",\"#bd0026\"],6:[\"#ffffb2\",\"#fed976\",\"#feb24c\",\"#fd8d3c\",\"#f03b20\",\"#bd0026\"],7:[\"#ffffb2\",\"#fed976\",\"#feb24c\",\"#fd8d3c\",\"#fc4e2a\",\"#e31a1c\",\"#b10026\"],8:[\"#ffffcc\",\"#ffeda0\",\"#fed976\",\"#feb24c\",\"#fd8d3c\",\"#fc4e2a\",\"#e31a1c\",\"#b10026\"],9:[\"#ffffcc\",\"#ffeda0\",\"#fed976\",\"#feb24c\",\"#fd8d3c\",\"#fc4e2a\",\"#e31a1c\",\"#bd0026\",\"#800026\"]},YlOrBr:{3:[\"#fff7bc\",\"#fec44f\",\"#d95f0e\"],4:[\"#ffffd4\",\"#fed98e\",\"#fe9929\",\"#cc4c02\"],5:[\"#ffffd4\",\"#fed98e\",\"#fe9929\",\"#d95f0e\",\"#993404\"],6:[\"#ffffd4\",\"#fee391\",\"#fec44f\",\"#fe9929\",\"#d95f0e\",\"#993404\"],7:[\"#ffffd4\",\"#fee391\",\"#fec44f\",\"#fe9929\",\"#ec7014\",\"#cc4c02\",\"#8c2d04\"],8:[\"#ffffe5\",\"#fff7bc\",\"#fee391\",\"#fec44f\",\"#fe9929\",\"#ec7014\",\"#cc4c02\",\"#8c2d04\"],9:[\"#ffffe5\",\"#fff7bc\",\"#fee391\",\"#fec44f\",\"#fe9929\",\"#ec7014\",\"#cc4c02\",\"#993404\",\"#662506\"]},Purples:{3:[\"#efedf5\",\"#bcbddc\",\"#756bb1\"],4:[\"#f2f0f7\",\"#cbc9e2\",\"#9e9ac8\",\"#6a51a3\"],5:[\"#f2f0f7\",\"#cbc9e2\",\"#9e9ac8\",\"#756bb1\",\"#54278f\"],6:[\"#f2f0f7\",\"#dadaeb\",\"#bcbddc\",\"#9e9ac8\",\"#756bb1\",\"#54278f\"],7:[\"#f2f0f7\",\"#dadaeb\",\"#bcbddc\",\"#9e9ac8\",\"#807dba\",\"#6a51a3\",\"#4a1486\"],8:[\"#fcfbfd\",\"#efedf5\",\"#dadaeb\",\"#bcbddc\",\"#9e9ac8\",\"#807dba\",\"#6a51a3\",\"#4a1486\"],9:[\"#fcfbfd\",\"#efedf5\",\"#dadaeb\",\"#bcbddc\",\"#9e9ac8\",\"#807dba\",\"#6a51a3\",\"#54278f\",\"#3f007d\"]},Blues:{3:[\"#deebf7\",\"#9ecae1\",\"#3182bd\"],4:[\"#eff3ff\",\"#bdd7e7\",\"#6baed6\",\"#2171b5\"],5:[\"#eff3ff\",\"#bdd7e7\",\"#6baed6\",\"#3182bd\",\"#08519c\"],6:[\"#eff3ff\",\"#c6dbef\",\"#9ecae1\",\"#6baed6\",\"#3182bd\",\"#08519c\"],7:[\"#eff3ff\",\"#c6dbef\",\"#9ecae1\",\"#6baed6\",\"#4292c6\",\"#2171b5\",\"#084594\"],8:[\"#f7fbff\",\"#deebf7\",\"#c6dbef\",\"#9ecae1\",\"#6baed6\",\"#4292c6\",\"#2171b5\",\"#084594\"],9:[\"#f7fbff\",\"#deebf7\",\"#c6dbef\",\"#9ecae1\",\"#6baed6\",\"#4292c6\",\"#2171b5\",\"#08519c\",\"#08306b\"]},Greens:{3:[\"#e5f5e0\",\"#a1d99b\",\"#31a354\"],4:[\"#edf8e9\",\"#bae4b3\",\"#74c476\",\"#238b45\"],5:[\"#edf8e9\",\"#bae4b3\",\"#74c476\",\"#31a354\",\"#006d2c\"],6:[\"#edf8e9\",\"#c7e9c0\",\"#a1d99b\",\"#74c476\",\"#31a354\",\"#006d2c\"],7:[\"#edf8e9\",\"#c7e9c0\",\"#a1d99b\",\"#74c476\",\"#41ab5d\",\"#238b45\",\"#005a32\"],8:[\"#f7fcf5\",\"#e5f5e0\",\"#c7e9c0\",\"#a1d99b\",\"#74c476\",\"#41ab5d\",\"#238b45\",\"#005a32\"],9:[\"#f7fcf5\",\"#e5f5e0\",\"#c7e9c0\",\"#a1d99b\",\"#74c476\",\"#41ab5d\",\"#238b45\",\"#006d2c\",\"#00441b\"]},Oranges:{3:[\"#fee6ce\",\"#fdae6b\",\"#e6550d\"],4:[\"#feedde\",\"#fdbe85\",\"#fd8d3c\",\"#d94701\"],5:[\"#feedde\",\"#fdbe85\",\"#fd8d3c\",\"#e6550d\",\"#a63603\"],6:[\"#feedde\",\"#fdd0a2\",\"#fdae6b\",\"#fd8d3c\",\"#e6550d\",\"#a63603\"],7:[\"#feedde\",\"#fdd0a2\",\"#fdae6b\",\"#fd8d3c\",\"#f16913\",\"#d94801\",\"#8c2d04\"],8:[\"#fff5eb\",\"#fee6ce\",\"#fdd0a2\",\"#fdae6b\",\"#fd8d3c\",\"#f16913\",\"#d94801\",\"#8c2d04\"],9:[\"#fff5eb\",\"#fee6ce\",\"#fdd0a2\",\"#fdae6b\",\"#fd8d3c\",\"#f16913\",\"#d94801\",\"#a63603\",\"#7f2704\"]},Reds:{3:[\"#fee0d2\",\"#fc9272\",\"#de2d26\"],4:[\"#fee5d9\",\"#fcae91\",\"#fb6a4a\",\"#cb181d\"],5:[\"#fee5d9\",\"#fcae91\",\"#fb6a4a\",\"#de2d26\",\"#a50f15\"],6:[\"#fee5d9\",\"#fcbba1\",\"#fc9272\",\"#fb6a4a\",\"#de2d26\",\"#a50f15\"],7:[\"#fee5d9\",\"#fcbba1\",\"#fc9272\",\"#fb6a4a\",\"#ef3b2c\",\"#cb181d\",\"#99000d\"],8:[\"#fff5f0\",\"#fee0d2\",\"#fcbba1\",\"#fc9272\",\"#fb6a4a\",\"#ef3b2c\",\"#cb181d\",\"#99000d\"],9:[\"#fff5f0\",\"#fee0d2\",\"#fcbba1\",\"#fc9272\",\"#fb6a4a\",\"#ef3b2c\",\"#cb181d\",\"#a50f15\",\"#67000d\"]},Greys:{3:[\"#f0f0f0\",\"#bdbdbd\",\"#636363\"],4:[\"#f7f7f7\",\"#cccccc\",\"#969696\",\"#525252\"],5:[\"#f7f7f7\",\"#cccccc\",\"#969696\",\"#636363\",\"#252525\"],6:[\"#f7f7f7\",\"#d9d9d9\",\"#bdbdbd\",\"#969696\",\"#636363\",\"#252525\"],7:[\"#f7f7f7\",\"#d9d9d9\",\"#bdbdbd\",\"#969696\",\"#737373\",\"#525252\",\"#252525\"],8:[\"#ffffff\",\"#f0f0f0\",\"#d9d9d9\",\"#bdbdbd\",\"#969696\",\"#737373\",\"#525252\",\"#252525\"],9:[\"#ffffff\",\"#f0f0f0\",\"#d9d9d9\",\"#bdbdbd\",\"#969696\",\"#737373\",\"#525252\",\"#252525\",\"#000000\"]},PuOr:{3:[\"#f1a340\",\"#f7f7f7\",\"#998ec3\"],4:[\"#e66101\",\"#fdb863\",\"#b2abd2\",\"#5e3c99\"],5:[\"#e66101\",\"#fdb863\",\"#f7f7f7\",\"#b2abd2\",\"#5e3c99\"],6:[\"#b35806\",\"#f1a340\",\"#fee0b6\",\"#d8daeb\",\"#998ec3\",\"#542788\"],7:[\"#b35806\",\"#f1a340\",\"#fee0b6\",\"#f7f7f7\",\"#d8daeb\",\"#998ec3\",\"#542788\"],8:[\"#b35806\",\"#e08214\",\"#fdb863\",\"#fee0b6\",\"#d8daeb\",\"#b2abd2\",\"#8073ac\",\"#542788\"],9:[\"#b35806\",\"#e08214\",\"#fdb863\",\"#fee0b6\",\"#f7f7f7\",\"#d8daeb\",\"#b2abd2\",\"#8073ac\",\"#542788\"],10:[\"#7f3b08\",\"#b35806\",\"#e08214\",\"#fdb863\",\"#fee0b6\",\"#d8daeb\",\"#b2abd2\",\"#8073ac\",\"#542788\",\"#2d004b\"],11:[\"#7f3b08\",\"#b35806\",\"#e08214\",\"#fdb863\",\"#fee0b6\",\"#f7f7f7\",\"#d8daeb\",\"#b2abd2\",\"#8073ac\",\"#542788\",\"#2d004b\"]},BrBG:{3:[\"#d8b365\",\"#f5f5f5\",\"#5ab4ac\"],4:[\"#a6611a\",\"#dfc27d\",\"#80cdc1\",\"#018571\"],5:[\"#a6611a\",\"#dfc27d\",\"#f5f5f5\",\"#80cdc1\",\"#018571\"],6:[\"#8c510a\",\"#d8b365\",\"#f6e8c3\",\"#c7eae5\",\"#5ab4ac\",\"#01665e\"],7:[\"#8c510a\",\"#d8b365\",\"#f6e8c3\",\"#f5f5f5\",\"#c7eae5\",\"#5ab4ac\",\"#01665e\"],8:[\"#8c510a\",\"#bf812d\",\"#dfc27d\",\"#f6e8c3\",\"#c7eae5\",\"#80cdc1\",\"#35978f\",\"#01665e\"],9:[\"#8c510a\",\"#bf812d\",\"#dfc27d\",\"#f6e8c3\",\"#f5f5f5\",\"#c7eae5\",\"#80cdc1\",\"#35978f\",\"#01665e\"],10:[\"#543005\",\"#8c510a\",\"#bf812d\",\"#dfc27d\",\"#f6e8c3\",\"#c7eae5\",\"#80cdc1\",\"#35978f\",\"#01665e\",\"#003c30\"],11:[\"#543005\",\"#8c510a\",\"#bf812d\",\"#dfc27d\",\"#f6e8c3\",\"#f5f5f5\",\"#c7eae5\",\"#80cdc1\",\"#35978f\",\"#01665e\",\"#003c30\"]},PRGn:{3:[\"#af8dc3\",\"#f7f7f7\",\"#7fbf7b\"],4:[\"#7b3294\",\"#c2a5cf\",\"#a6dba0\",\"#008837\"],5:[\"#7b3294\",\"#c2a5cf\",\"#f7f7f7\",\"#a6dba0\",\"#008837\"],6:[\"#762a83\",\"#af8dc3\",\"#e7d4e8\",\"#d9f0d3\",\"#7fbf7b\",\"#1b7837\"],7:[\"#762a83\",\"#af8dc3\",\"#e7d4e8\",\"#f7f7f7\",\"#d9f0d3\",\"#7fbf7b\",\"#1b7837\"],8:[\"#762a83\",\"#9970ab\",\"#c2a5cf\",\"#e7d4e8\",\"#d9f0d3\",\"#a6dba0\",\"#5aae61\",\"#1b7837\"],9:[\"#762a83\",\"#9970ab\",\"#c2a5cf\",\"#e7d4e8\",\"#f7f7f7\",\"#d9f0d3\",\"#a6dba0\",\"#5aae61\",\"#1b7837\"],10:[\"#40004b\",\"#762a83\",\"#9970ab\",\"#c2a5cf\",\"#e7d4e8\",\"#d9f0d3\",\"#a6dba0\",\"#5aae61\",\"#1b7837\",\"#00441b\"],11:[\"#40004b\",\"#762a83\",\"#9970ab\",\"#c2a5cf\",\"#e7d4e8\",\"#f7f7f7\",\"#d9f0d3\",\"#a6dba0\",\"#5aae61\",\"#1b7837\",\"#00441b\"]},PiYG:{3:[\"#e9a3c9\",\"#f7f7f7\",\"#a1d76a\"],4:[\"#d01c8b\",\"#f1b6da\",\"#b8e186\",\"#4dac26\"],5:[\"#d01c8b\",\"#f1b6da\",\"#f7f7f7\",\"#b8e186\",\"#4dac26\"],6:[\"#c51b7d\",\"#e9a3c9\",\"#fde0ef\",\"#e6f5d0\",\"#a1d76a\",\"#4d9221\"],7:[\"#c51b7d\",\"#e9a3c9\",\"#fde0ef\",\"#f7f7f7\",\"#e6f5d0\",\"#a1d76a\",\"#4d9221\"],8:[\"#c51b7d\",\"#de77ae\",\"#f1b6da\",\"#fde0ef\",\"#e6f5d0\",\"#b8e186\",\"#7fbc41\",\"#4d9221\"],9:[\"#c51b7d\",\"#de77ae\",\"#f1b6da\",\"#fde0ef\",\"#f7f7f7\",\"#e6f5d0\",\"#b8e186\",\"#7fbc41\",\"#4d9221\"],10:[\"#8e0152\",\"#c51b7d\",\"#de77ae\",\"#f1b6da\",\"#fde0ef\",\"#e6f5d0\",\"#b8e186\",\"#7fbc41\",\"#4d9221\",\"#276419\"],11:[\"#8e0152\",\"#c51b7d\",\"#de77ae\",\"#f1b6da\",\"#fde0ef\",\"#f7f7f7\",\"#e6f5d0\",\"#b8e186\",\"#7fbc41\",\"#4d9221\",\"#276419\"]},RdBu:{3:[\"#ef8a62\",\"#f7f7f7\",\"#67a9cf\"],4:[\"#ca0020\",\"#f4a582\",\"#92c5de\",\"#0571b0\"],5:[\"#ca0020\",\"#f4a582\",\"#f7f7f7\",\"#92c5de\",\"#0571b0\"],6:[\"#b2182b\",\"#ef8a62\",\"#fddbc7\",\"#d1e5f0\",\"#67a9cf\",\"#2166ac\"],7:[\"#b2182b\",\"#ef8a62\",\"#fddbc7\",\"#f7f7f7\",\"#d1e5f0\",\"#67a9cf\",\"#2166ac\"],8:[\"#b2182b\",\"#d6604d\",\"#f4a582\",\"#fddbc7\",\"#d1e5f0\",\"#92c5de\",\"#4393c3\",\"#2166ac\"],9:[\"#b2182b\",\"#d6604d\",\"#f4a582\",\"#fddbc7\",\"#f7f7f7\",\"#d1e5f0\",\"#92c5de\",\"#4393c3\",\"#2166ac\"],10:[\"#67001f\",\"#b2182b\",\"#d6604d\",\"#f4a582\",\"#fddbc7\",\"#d1e5f0\",\"#92c5de\",\"#4393c3\",\"#2166ac\",\"#053061\"],11:[\"#67001f\",\"#b2182b\",\"#d6604d\",\"#f4a582\",\"#fddbc7\",\"#f7f7f7\",\"#d1e5f0\",\"#92c5de\",\"#4393c3\",\"#2166ac\",\"#053061\"]},RdGy:{3:[\"#ef8a62\",\"#ffffff\",\"#999999\"],4:[\"#ca0020\",\"#f4a582\",\"#bababa\",\"#404040\"],5:[\"#ca0020\",\"#f4a582\",\"#ffffff\",\"#bababa\",\"#404040\"],6:[\"#b2182b\",\"#ef8a62\",\"#fddbc7\",\"#e0e0e0\",\"#999999\",\"#4d4d4d\"],7:[\"#b2182b\",\"#ef8a62\",\"#fddbc7\",\"#ffffff\",\"#e0e0e0\",\"#999999\",\"#4d4d4d\"],8:[\"#b2182b\",\"#d6604d\",\"#f4a582\",\"#fddbc7\",\"#e0e0e0\",\"#bababa\",\"#878787\",\"#4d4d4d\"],9:[\"#b2182b\",\"#d6604d\",\"#f4a582\",\"#fddbc7\",\"#ffffff\",\"#e0e0e0\",\"#bababa\",\"#878787\",\"#4d4d4d\"],10:[\"#67001f\",\"#b2182b\",\"#d6604d\",\"#f4a582\",\"#fddbc7\",\"#e0e0e0\",\"#bababa\",\"#878787\",\"#4d4d4d\",\"#1a1a1a\"],11:[\"#67001f\",\"#b2182b\",\"#d6604d\",\"#f4a582\",\"#fddbc7\",\"#ffffff\",\"#e0e0e0\",\"#bababa\",\"#878787\",\"#4d4d4d\",\"#1a1a1a\"]},RdYlBu:{3:[\"#fc8d59\",\"#ffffbf\",\"#91bfdb\"],4:[\"#d7191c\",\"#fdae61\",\"#abd9e9\",\"#2c7bb6\"],5:[\"#d7191c\",\"#fdae61\",\"#ffffbf\",\"#abd9e9\",\"#2c7bb6\"],6:[\"#d73027\",\"#fc8d59\",\"#fee090\",\"#e0f3f8\",\"#91bfdb\",\"#4575b4\"],7:[\"#d73027\",\"#fc8d59\",\"#fee090\",\"#ffffbf\",\"#e0f3f8\",\"#91bfdb\",\"#4575b4\"],8:[\"#d73027\",\"#f46d43\",\"#fdae61\",\"#fee090\",\"#e0f3f8\",\"#abd9e9\",\"#74add1\",\"#4575b4\"],9:[\"#d73027\",\"#f46d43\",\"#fdae61\",\"#fee090\",\"#ffffbf\",\"#e0f3f8\",\"#abd9e9\",\"#74add1\",\"#4575b4\"],10:[\"#a50026\",\"#d73027\",\"#f46d43\",\"#fdae61\",\"#fee090\",\"#e0f3f8\",\"#abd9e9\",\"#74add1\",\"#4575b4\",\"#313695\"],11:[\"#a50026\",\"#d73027\",\"#f46d43\",\"#fdae61\",\"#fee090\",\"#ffffbf\",\"#e0f3f8\",\"#abd9e9\",\"#74add1\",\"#4575b4\",\"#313695\"]},Spectral:{3:[\"#fc8d59\",\"#ffffbf\",\"#99d594\"],4:[\"#d7191c\",\"#fdae61\",\"#abdda4\",\"#2b83ba\"],5:[\"#d7191c\",\"#fdae61\",\"#ffffbf\",\"#abdda4\",\"#2b83ba\"],6:[\"#d53e4f\",\"#fc8d59\",\"#fee08b\",\"#e6f598\",\"#99d594\",\"#3288bd\"],7:[\"#d53e4f\",\"#fc8d59\",\"#fee08b\",\"#ffffbf\",\"#e6f598\",\"#99d594\",\"#3288bd\"],8:[\"#d53e4f\",\"#f46d43\",\"#fdae61\",\"#fee08b\",\"#e6f598\",\"#abdda4\",\"#66c2a5\",\"#3288bd\"],9:[\"#d53e4f\",\"#f46d43\",\"#fdae61\",\"#fee08b\",\"#ffffbf\",\"#e6f598\",\"#abdda4\",\"#66c2a5\",\"#3288bd\"],10:[\"#9e0142\",\"#d53e4f\",\"#f46d43\",\"#fdae61\",\"#fee08b\",\"#e6f598\",\"#abdda4\",\"#66c2a5\",\"#3288bd\",\"#5e4fa2\"],11:[\"#9e0142\",\"#d53e4f\",\"#f46d43\",\"#fdae61\",\"#fee08b\",\"#ffffbf\",\"#e6f598\",\"#abdda4\",\"#66c2a5\",\"#3288bd\",\"#5e4fa2\"]},RdYlGn:{3:[\"#fc8d59\",\"#ffffbf\",\"#91cf60\"],4:[\"#d7191c\",\"#fdae61\",\"#a6d96a\",\"#1a9641\"],5:[\"#d7191c\",\"#fdae61\",\"#ffffbf\",\"#a6d96a\",\"#1a9641\"],6:[\"#d73027\",\"#fc8d59\",\"#fee08b\",\"#d9ef8b\",\"#91cf60\",\"#1a9850\"],7:[\"#d73027\",\"#fc8d59\",\"#fee08b\",\"#ffffbf\",\"#d9ef8b\",\"#91cf60\",\"#1a9850\"],8:[\"#d73027\",\"#f46d43\",\"#fdae61\",\"#fee08b\",\"#d9ef8b\",\"#a6d96a\",\"#66bd63\",\"#1a9850\"],9:[\"#d73027\",\"#f46d43\",\"#fdae61\",\"#fee08b\",\"#ffffbf\",\"#d9ef8b\",\"#a6d96a\",\"#66bd63\",\"#1a9850\"],10:[\"#a50026\",\"#d73027\",\"#f46d43\",\"#fdae61\",\"#fee08b\",\"#d9ef8b\",\"#a6d96a\",\"#66bd63\",\"#1a9850\",\"#006837\"],11:[\"#a50026\",\"#d73027\",\"#f46d43\",\"#fdae61\",\"#fee08b\",\"#ffffbf\",\"#d9ef8b\",\"#a6d96a\",\"#66bd63\",\"#1a9850\",\"#006837\"]},Accent:{3:[\"#7fc97f\",\"#beaed4\",\"#fdc086\"],4:[\"#7fc97f\",\"#beaed4\",\"#fdc086\",\"#ffff99\"],5:[\"#7fc97f\",\"#beaed4\",\"#fdc086\",\"#ffff99\",\"#386cb0\"],6:[\"#7fc97f\",\"#beaed4\",\"#fdc086\",\"#ffff99\",\"#386cb0\",\"#f0027f\"],7:[\"#7fc97f\",\"#beaed4\",\"#fdc086\",\"#ffff99\",\"#386cb0\",\"#f0027f\",\"#bf5b17\"],8:[\"#7fc97f\",\"#beaed4\",\"#fdc086\",\"#ffff99\",\"#386cb0\",\"#f0027f\",\"#bf5b17\",\"#666666\"]},Dark2:{3:[\"#1b9e77\",\"#d95f02\",\"#7570b3\"],4:[\"#1b9e77\",\"#d95f02\",\"#7570b3\",\"#e7298a\"],5:[\"#1b9e77\",\"#d95f02\",\"#7570b3\",\"#e7298a\",\"#66a61e\"],6:[\"#1b9e77\",\"#d95f02\",\"#7570b3\",\"#e7298a\",\"#66a61e\",\"#e6ab02\"],7:[\"#1b9e77\",\"#d95f02\",\"#7570b3\",\"#e7298a\",\"#66a61e\",\"#e6ab02\",\"#a6761d\"],8:[\"#1b9e77\",\"#d95f02\",\"#7570b3\",\"#e7298a\",\"#66a61e\",\"#e6ab02\",\"#a6761d\",\"#666666\"]},Paired:{3:[\"#a6cee3\",\"#1f78b4\",\"#b2df8a\"],4:[\"#a6cee3\",\"#1f78b4\",\"#b2df8a\",\"#33a02c\"],5:[\"#a6cee3\",\"#1f78b4\",\"#b2df8a\",\"#33a02c\",\"#fb9a99\"],6:[\"#a6cee3\",\"#1f78b4\",\"#b2df8a\",\"#33a02c\",\"#fb9a99\",\"#e31a1c\"],7:[\"#a6cee3\",\"#1f78b4\",\"#b2df8a\",\"#33a02c\",\"#fb9a99\",\"#e31a1c\",\"#fdbf6f\"],8:[\"#a6cee3\",\"#1f78b4\",\"#b2df8a\",\"#33a02c\",\"#fb9a99\",\"#e31a1c\",\"#fdbf6f\",\"#ff7f00\"],9:[\"#a6cee3\",\"#1f78b4\",\"#b2df8a\",\"#33a02c\",\"#fb9a99\",\"#e31a1c\",\"#fdbf6f\",\"#ff7f00\",\"#cab2d6\"],10:[\"#a6cee3\",\"#1f78b4\",\"#b2df8a\",\"#33a02c\",\"#fb9a99\",\"#e31a1c\",\"#fdbf6f\",\"#ff7f00\",\"#cab2d6\",\"#6a3d9a\"],11:[\"#a6cee3\",\"#1f78b4\",\"#b2df8a\",\"#33a02c\",\"#fb9a99\",\"#e31a1c\",\"#fdbf6f\",\"#ff7f00\",\"#cab2d6\",\"#6a3d9a\",\"#ffff99\"],12:[\"#a6cee3\",\"#1f78b4\",\"#b2df8a\",\"#33a02c\",\"#fb9a99\",\"#e31a1c\",\"#fdbf6f\",\"#ff7f00\",\"#cab2d6\",\"#6a3d9a\",\"#ffff99\",\"#b15928\"]},Pastel1:{3:[\"#fbb4ae\",\"#b3cde3\",\"#ccebc5\"],4:[\"#fbb4ae\",\"#b3cde3\",\"#ccebc5\",\"#decbe4\"],5:[\"#fbb4ae\",\"#b3cde3\",\"#ccebc5\",\"#decbe4\",\"#fed9a6\"],6:[\"#fbb4ae\",\"#b3cde3\",\"#ccebc5\",\"#decbe4\",\"#fed9a6\",\"#ffffcc\"],7:[\"#fbb4ae\",\"#b3cde3\",\"#ccebc5\",\"#decbe4\",\"#fed9a6\",\"#ffffcc\",\"#e5d8bd\"],8:[\"#fbb4ae\",\"#b3cde3\",\"#ccebc5\",\"#decbe4\",\"#fed9a6\",\"#ffffcc\",\"#e5d8bd\",\"#fddaec\"],9:[\"#fbb4ae\",\"#b3cde3\",\"#ccebc5\",\"#decbe4\",\"#fed9a6\",\"#ffffcc\",\"#e5d8bd\",\"#fddaec\",\"#f2f2f2\"]},Pastel2:{3:[\"#b3e2cd\",\"#fdcdac\",\"#cbd5e8\"],4:[\"#b3e2cd\",\"#fdcdac\",\"#cbd5e8\",\"#f4cae4\"],5:[\"#b3e2cd\",\"#fdcdac\",\"#cbd5e8\",\"#f4cae4\",\"#e6f5c9\"],6:[\"#b3e2cd\",\"#fdcdac\",\"#cbd5e8\",\"#f4cae4\",\"#e6f5c9\",\"#fff2ae\"],7:[\"#b3e2cd\",\"#fdcdac\",\"#cbd5e8\",\"#f4cae4\",\"#e6f5c9\",\"#fff2ae\",\"#f1e2cc\"],8:[\"#b3e2cd\",\"#fdcdac\",\"#cbd5e8\",\"#f4cae4\",\"#e6f5c9\",\"#fff2ae\",\"#f1e2cc\",\"#cccccc\"]},Set1:{3:[\"#e41a1c\",\"#377eb8\",\"#4daf4a\"],4:[\"#e41a1c\",\"#377eb8\",\"#4daf4a\",\"#984ea3\"],5:[\"#e41a1c\",\"#377eb8\",\"#4daf4a\",\"#984ea3\",\"#ff7f00\"],6:[\"#e41a1c\",\"#377eb8\",\"#4daf4a\",\"#984ea3\",\"#ff7f00\",\"#ffff33\"],7:[\"#e41a1c\",\"#377eb8\",\"#4daf4a\",\"#984ea3\",\"#ff7f00\",\"#ffff33\",\"#a65628\"],8:[\"#e41a1c\",\"#377eb8\",\"#4daf4a\",\"#984ea3\",\"#ff7f00\",\"#ffff33\",\"#a65628\",\"#f781bf\"],9:[\"#e41a1c\",\"#377eb8\",\"#4daf4a\",\"#984ea3\",\"#ff7f00\",\"#ffff33\",\"#a65628\",\"#f781bf\",\"#999999\"]},Set2:{3:[\"#66c2a5\",\"#fc8d62\",\"#8da0cb\"],4:[\"#66c2a5\",\"#fc8d62\",\"#8da0cb\",\"#e78ac3\"],5:[\"#66c2a5\",\"#fc8d62\",\"#8da0cb\",\"#e78ac3\",\"#a6d854\"],6:[\"#66c2a5\",\"#fc8d62\",\"#8da0cb\",\"#e78ac3\",\"#a6d854\",\"#ffd92f\"],7:[\"#66c2a5\",\"#fc8d62\",\"#8da0cb\",\"#e78ac3\",\"#a6d854\",\"#ffd92f\",\"#e5c494\"],8:[\"#66c2a5\",\"#fc8d62\",\"#8da0cb\",\"#e78ac3\",\"#a6d854\",\"#ffd92f\",\"#e5c494\",\"#b3b3b3\"]},Set3:{3:[\"#8dd3c7\",\"#ffffb3\",\"#bebada\"],4:[\"#8dd3c7\",\"#ffffb3\",\"#bebada\",\"#fb8072\"],5:[\"#8dd3c7\",\"#ffffb3\",\"#bebada\",\"#fb8072\",\"#80b1d3\"],6:[\"#8dd3c7\",\"#ffffb3\",\"#bebada\",\"#fb8072\",\"#80b1d3\",\"#fdb462\"],7:[\"#8dd3c7\",\"#ffffb3\",\"#bebada\",\"#fb8072\",\"#80b1d3\",\"#fdb462\",\"#b3de69\"],8:[\"#8dd3c7\",\"#ffffb3\",\"#bebada\",\"#fb8072\",\"#80b1d3\",\"#fdb462\",\"#b3de69\",\"#fccde5\"],9:[\"#8dd3c7\",\"#ffffb3\",\"#bebada\",\"#fb8072\",\"#80b1d3\",\"#fdb462\",\"#b3de69\",\"#fccde5\",\"#d9d9d9\"],10:[\"#8dd3c7\",\"#ffffb3\",\"#bebada\",\"#fb8072\",\"#80b1d3\",\"#fdb462\",\"#b3de69\",\"#fccde5\",\"#d9d9d9\",\"#bc80bd\"],11:[\"#8dd3c7\",\"#ffffb3\",\"#bebada\",\"#fb8072\",\"#80b1d3\",\"#fdb462\",\"#b3de69\",\"#fccde5\",\"#d9d9d9\",\"#bc80bd\",\"#ccebc5\"],12:[\"#8dd3c7\",\"#ffffb3\",\"#bebada\",\"#fb8072\",\"#80b1d3\",\"#fdb462\",\"#b3de69\",\"#fccde5\",\"#d9d9d9\",\"#bc80bd\",\"#ccebc5\",\"#ffed6f\"]}};" << std::endl;
  out << "" << std::endl;
  out << "</script>" << std::endl;
  out << "" << std::endl;
  out << "<body>" << std::endl;
  out << "" << std::endl;
  out << "<div id=\"figure\" style=\"margin-bottom: 50px;\"></div>" << std::endl;
  out << "" << std::endl;
  out << "" << std::endl;
  out << "" << std::endl;
  out << "<script>" << std::endl;
  out << "" << std::endl;
  out << "var margin = {top: 50, right: 10, bottom: 10, left: 200}," << std::endl;
  out << "    width = window.innerWidth * 1 - margin.left - margin.right," << std::endl;
  out << "    height = 10000 - margin.top - margin.bottom;" << std::endl;
  out << "" << std::endl;
  out << "var y = d3.scale.ordinal()" << std::endl;
  out << "    .rangeRoundBands([0, height], .3);" << std::endl;
  out << "" << std::endl;
  out << "var x = d3.scale.linear()" << std::endl;
  out << "    .rangeRound([0, width]);" << std::endl;
  out << "" << std::endl;
  out << "" << std::endl;
  out << "var xAxis = d3.svg.axis()" << std::endl;
  out << "    .scale(x)" << std::endl;
  out << "    .orient(\"top\")" << std::endl;
  out << "    .ticks(100);" << std::endl;
  out << "" << std::endl;
  out << "var yAxis = d3.svg.axis()" << std::endl;
  out << "    .scale(y)" << std::endl;
  out << "    .orient(\"left\")" << std::endl;
  out << "" << std::endl;
  out << "var svg = d3.select(\"#figure\").append(\"svg\")" << std::endl;
  out << "    .attr(\"width\", width + margin.left + margin.right)" << std::endl;
  out << "    .attr(\"height\", height + margin.top + margin.bottom)" << std::endl;
  out << "    .attr(\"id\", \"d3-plot\")" << std::endl;
  out << "  .append(\"g\")" << std::endl;
  out << "    .attr(\"transform\", \"translate(\" + margin.left + \",\" + margin.top + \")\");" << std::endl;
  out << "" << std::endl;
  out << "" << std::endl;
  out << "  d3.json(\"elements.json\", function(error, data) {" << std::endl;
  out << "" << std::endl;
  out << "" << std::endl;
  out << "  x.domain([0, 10000]).nice();" << std::endl;
  out << "" << std::endl;
  out << "  y.domain(data.map(function(d){return d.name}));" << std::endl;
  out << "var colors = d3.scale.ordinal()" << std::endl;
  out << "    .domain(data.map(function(d){return d.elements.map(function(subd){return subd.nameOfQuery_;});}))" << std::endl;
  out << "    .range(d3.scale.category20b().range().concat(d3.scale.category20c().range()));" << std::endl;
  out << "    //.range(colorbrewer.Set1[9]);" << std::endl;
  out << "" << std::endl;
  out << "  svg.append(\"g\")" << std::endl;
  out << "      .attr(\"class\", \"x axis\")" << std::endl;
  out << "      .call(xAxis);" << std::endl;
  out << "" << std::endl;
  out << "  svg.append(\"g\")" << std::endl;
  out << "      .attr(\"class\", \"y axis\")" << std::endl;
  out << "      .call(yAxis)" << std::endl;
  out << "  var genomes = svg.selectAll(\".genomes\")" << std::endl;
  out << "      .data(data)" << std::endl;
  out << "    .enter().append(\"g\")" << std::endl;
  out << "      .attr(\"class\", \"bar\")" << std::endl;
  out << "      .attr(\"id\", function(d) {return d.name;})" << std::endl;
  out << "      .attr(\"transform\", function(d) {return \"translate(0,\" + y(d.name) + \")\"; });" << std::endl;
  out << "	//console.log(d.name.toString()); console.log(y(d.name)); " << std::endl;
  out << "  var bars = genomes.selectAll(\"rect\")" << std::endl;
  out << "      .data(function(d) {return d.elements; })" << std::endl;
  out << "    .enter()" << std::endl;
  out << "    	.append(\"g\")" << std::endl;
  out << "    	.attr(\"class\", function(d) {return \"subbar \" + d.nameOfMatchedSeq_; } );" << std::endl;
  out << "  d3.selectAll(\".LTR7\").remove()" << std::endl;
  out << "  d3.selectAll(\".LTR7B\").remove()" << std::endl;
  out << "  d3.selectAll(\".LTR7C\").remove()" << std::endl;
  out << "  d3.selectAll(\".LTR7Y\").remove()" << std::endl;
  out << "" << std::endl;
  out << "  var tooltip = d3.select(\"body\")" << std::endl;
  out << "				.append(\"div\")" << std::endl;
  out << "				.style(\"position\", \"absolute\")" << std::endl;
  out << "				.style(\"visibility\", \"hidden\")" << std::endl;
  out << "				.style(\"background-color\", \"#88aaaa\")" << std::endl;
  out << "				.style(\"width\", \"300\")" << std::endl;
  out << "				.attr(\"id\", \"popHover\");" << std::endl;
  out << "  bars.append(\"rect\")" << std::endl;
  out << "      .attr(\"height\", y.rangeBand())" << std::endl;
  out << "      .attr(\"x\", function(d) { return x(d.startInMatch_); })" << std::endl;
  out << "      .attr(\"width\", function(d) { return x(d.endInMatch_ - d.startInMatch_ + 1); })" << std::endl;
  out << "      .attr(\"id\", function(d) { return \"ref\" + d.startInMatch_.toString(); })" << std::endl;
  out << "      //.attr(\"fill\", function(d){ return colors(d.nameOfMatchedSeq_);})" << std::endl;
  out << "      .attr(\"fill\", function(d){ return colors(d.nameOfQuery_);})" << std::endl;
  out << "      .on(\"mouseover\", function(d,i,j){" << std::endl;
  out << "	    	  tooltip.node().innerHTML = \"name: \"  + d.nameOfMatchedSeq_ + \"<br>refRange: \" + d.startInMatch_.toString() + \":\"  + d.endInMatch_.toString() + \"<br>genome: \" + d.nameOfQuery_ + \":\" + d.start_.toString() + \":\" + d.end_.toString() ;" << std::endl;
  out << "	    	  tooltip.style(\"visibility\", \"visible\");" << std::endl;
  out << "	      	 return;})" << std::endl;
  out << "	  .on(\"mousemove\", function(){" << std::endl;
  out << "	  	tooltip.style(\"top\", (d3.event.layerY-10)+\"px\").style(\"left\",(d3.event.layerX+10)+\"px\")" << std::endl;
  out << "	  	return ;})" << std::endl;
  out << "	  .on(\"mouseout\", function(d){" << std::endl;
  out << "	  	tooltip.style(\"visibility\", \"hidden\");" << std::endl;
  out << "	  return ;});" << std::endl;
  out << "" << std::endl;
  out << "  svg.append(\"g\")" << std::endl;
  out << "      .attr(\"class\", \"y axis\")" << std::endl;
  out << "  .append(\"line\")" << std::endl;
  out << "      .attr(\"x1\", x(0))" << std::endl;
  out << "      .attr(\"x2\", x(0))" << std::endl;
  out << "      .attr(\"y2\", height);" << std::endl;
  out << "" << std::endl;
  out << "" << std::endl;
  out << "  d3.selectAll(\".axis path\")" << std::endl;
  out << "      .style(\"fill\", \"none\")" << std::endl;
  out << "      .style(\"stroke\", \"#000\")" << std::endl;
  out << "      .style(\"shape-rendering\", \"crispEdges\")" << std::endl;
  out << "" << std::endl;
  out << "  d3.selectAll(\".axis line\")" << std::endl;
  out << "      .style(\"fill\", \"none\")" << std::endl;
  out << "      .style(\"stroke\", \"#000\")" << std::endl;
  out << "      .style(\"shape-rendering\", \"crispEdges\")" << std::endl;
  out << "" << std::endl;
  out << "});" << std::endl;
  out << "</script>" << std::endl;
  out << "" << std::endl;
  out << "</body>" << std::endl;
  out << "</html>" << std::endl;
}

int repelinRunner::extractElementSequences(const njh::progutils::CmdArgs & inputCommands){
  seqSetUp setUp(inputCommands);
  std::string repeatFilename = "";
  std::string twoBitFilename = "";
  setUp.setOption(repeatFilename, "--repeatFilename", "Filename of repeat masker output", true);
  setUp.setOption(twoBitFilename, "--twoBitFile", "Filename of two bit file", true);
  setUp.pars_.ioOptions_.firstName_ = repeatFilename;
  setUp.processDirectoryOutputName(true);
  setUp.processVerbose();
  setUp.finishSetUp(std::cout);
  setUp.startARunLog(setUp.pars_.directoryName_);
  std::string eachElementDir = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar("repeatElements")).string();
	std::ofstream outFile;
	MultiSeqIO writer;
	checkFileExistsThrow(twoBitFilename);
	TwoBit::TwoBitFile f(twoBitFilename);
	std::string buffer = "";
	std::unordered_map<uint32_t, std::vector<RepeatMaskerRecord>> elements =
			readRepeatMaskerFile(repeatFilename, setUp.pars_.verbose_);
	auto eleKeys = getVectorOfMapKeys(elements);
	njh::sort(eleKeys);

	Json::Value outVisJson;
	for(const auto & eKey : eleKeys){
		const auto & e = elements[eKey];
		Json::Value eJson;
		eJson["name"] = estd::to_string(eKey);
		eJson["chr"] = estd::to_string(e.front().nameOfQuery_);
		auto & subEJsons = eJson["elements"];
		for(const auto & subE : e){
			subEJsons.append(subE.toJson());
			njh::files::makeDirP(eachElementDir, njh::files::MkdirPar(subE.nameOfMatchedSeq_));
			if (!writer.containsReader(subE.nameOfMatchedSeq_)) {
				auto outOpts = SeqIOOptions::genFastaOut(
						njh::files::make_path(eachElementDir, subE.nameOfMatchedSeq_,
								subE.nameOfMatchedSeq_).string());
				writer.addReader(subE.nameOfMatchedSeq_, outOpts);
			}
			f[subE.nameOfQuery_]->getSequence(buffer,subE.start_ - 1, subE.end_);
			if(subE.reverseStrand_){
				buffer = seqUtil::reverseComplement(buffer, "DNA");
			}
			std::stringstream outName;
			outName << subE.regionSegment_ << "." << subE.nameOfMatchedSeq_ << "["
					<< "repeatMaskerId=" << subE.regionSegment_ << ";"
					<< "repeatType=" << subE.nameOfMatchedSeq_ << ";"
					<< "rStart=" << subE.startInMatch_.first - 1 << ";"
					<< "rStop=" << subE.endInMatch_.first << ";"
					<< "chr=" << subE.nameOfQuery_ << ";"
					<< "gStart=" << subE.start_ - 1 << ";"
					<< "gStop=" << subE.end_ << ";"
					<< "strand=" << (subE.reverseStrand_ ? '-' : '+') << ";"
					<< "]";
			seqInfo outSeq = seqInfo(outName.str(), stringToUpperReturn(buffer));
			writer.openWrite(subE.nameOfMatchedSeq_, outSeq);
		}
		outVisJson.append(eJson);
	}
	std::ofstream outJson;
	openTextFile(outJson, setUp.pars_.directoryName_ + "elements.json", ".json", false, false);
	outJson << outVisJson;

	std::ofstream outHtml;
	openTextFile(outHtml, setUp.pars_.directoryName_ + "elements.html", ".html", false, false);
	genElementsHtml(outHtml);

	return 0;
}



int repelinRunner::extractFullElementSequences(const njh::progutils::CmdArgs & inputCommands){
  seqSetUp setUp(inputCommands);
  std::string repeatFilename = "";
  std::string twoBitFilename = "";
  bool toUpper = false;
  setUp.setOption(repeatFilename, "--repeatFilename", "Filename of repeat masker output", true);
  setUp.setOption(twoBitFilename, "--twoBitFile", "Filename of two bit file", true);
  setUp.setOption(toUpper, "--toUpper", "Make all Bases Upper case");
  setUp.pars_.ioOptions_.firstName_ = repeatFilename;
  setUp.processDirectoryOutputName(true);
  setUp.processVerbose();
  setUp.finishSetUp(std::cout);
  setUp.startARunLog(setUp.pars_.directoryName_);
  std::string eachElementDir = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar("eachElement")).string();
  std::ofstream allFasta;
  openTextFile(allFasta, setUp.pars_.directoryName_ + "all", ".fasta", false, false);
	checkFileExistsThrow(twoBitFilename);
	TwoBit::TwoBitFile f(twoBitFilename);
	std::string buffer;
	std::unordered_map<uint32_t, std::vector<RepeatMaskerRecord>> elements =
			readRepeatMaskerFile(repeatFilename, setUp.pars_.verbose_);
	auto eleKeys = getVectorOfMapKeys(elements);
	njh::sort(eleKeys);
	for(const auto & eKey : eleKeys){
		const auto & e = elements[eKey];
		std::vector<seqInfo> seqs;
		std::string buffer = "";
		uint32_t minPos = std::numeric_limits<uint32_t>::max();
		uint32_t maxPos = std::numeric_limits<uint32_t>::min();
		for(const auto & subE : e){
			f[subE.nameOfQuery_]->getSequence(buffer,subE.start_,subE.end_);
			if(minPos > subE.start_){
				minPos = subE.start_;
			}
			if(maxPos < subE.end_){
				maxPos = subE.end_;
			}
			if(subE.reverseStrand_){
				buffer = seqUtil::reverseComplement(buffer, "DNA");
			}
			if(toUpper){
				stringToUpper(buffer);
			}
			seqs.emplace_back(njh::err::F() << subE.nameOfMatchedSeq_ << ";" << subE.startInMatch_.first<< ";" << subE.endInMatch_.first<< ";" , buffer);
		}

		f[e.front().nameOfQuery_]->getSequence(buffer,minPos,maxPos);
		if(e.front().reverseStrand_){
			buffer = seqUtil::reverseComplement(buffer, "DNA");
		}
		if(toUpper){
			stringToUpper(buffer);
		}

		std::ofstream outFile;
		openTextFile(outFile,
				eachElementDir + estd::to_string(e.front().regionSegment_), ".fasta",
				false, false);
		outFile << ">" << e.front().regionSegment_ << "_full" << std::endl;
		outFile << buffer << std::endl;
		for(const auto & seq : seqs){
			seq.outPutSeq(outFile);
		}
		allFasta << ">" << e.front().regionSegment_ << std::endl;
		allFasta << buffer << std::endl;
	}
	return 0;
}

int repelinRunner::parseTandemRepeatFinderOutputUnitTest(const njh::progutils::CmdArgs & inputCommands){
  seqSetUp setUp(inputCommands);
  bfs::path repeatFilename = "";
  setUp.setOption(repeatFilename, "--repeatFilename", "Filename of trf output", true);
  setUp.processVerbose();
  setUp.finishSetUp(std::cout);

  auto trfRecrods = getTandemRepeatFinderRecords(repeatFilename);
  if(setUp.pars_.verbose_){
  	std::cout << "Read : " << trfRecrods.size() << " records" << std::endl;
  	std::map<std::string, uint32_t> recrodCountsBySeqName;
  	for(const auto & record : trfRecrods){
  		++recrodCountsBySeqName[record->seqName_];
  	}
  	printOutMapContents(recrodCountsBySeqName, "\t", std::cout);
  }
  return 0;
}

int repelinRunner::TandemRepeatFinderOutputToBed(const njh::progutils::CmdArgs & inputCommands){
  seqSetUp setUp(inputCommands);
  bfs::path repeatFilename = "";
  OutOptions outOpts{bfs::path("")};
  outOpts.outExtention_ = ".bed";
  setUp.processVerbose();
  setUp.setOption(repeatFilename, "--repeatFilename", "Filename of trf output", true);
  setUp.processWritingOptions(outOpts);
  setUp.finishSetUp(std::cout);

  std::ofstream outFile;
  std::ostream out(determineOutBuf(outFile, outOpts));
	BioDataFileIO<TandemRepeatFinderRecord> reader{IoOptions(InOptions(repeatFilename))};
	reader.openIn();
	std::string currentSeqName = "";
	std::string line = "";
	while (njh::files::crossPlatGetline(*reader.inFile_, line)) {
		if (njh::beginsWith(line, "Sequence")) {
			auto toks = njh::tokenizeString(line, "whitespace");
			if (toks.size() < 2) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error in processing line " << line
						<< ", was expecting at least two values separated by white sapce but found "
						<< toks.size() << " instead " << "\n";
				throw std::runtime_error { ss.str() };
			}
			currentSeqName = toks[1];
		}
		//skip lines that don't start with a digit or is a blank line, it seems to be the only indicator that the file is on a data line, this format is stupid
		if("" == line || allWhiteSpaceStr(line) || !::isdigit(line.front())){
			continue;
		}
		std::shared_ptr<TandemRepeatFinderRecord> record = std::make_shared<TandemRepeatFinderRecord>(line);
		record->setSeqName(currentSeqName);
		out << GenomicRegion(*record).genBedRecordCore().toDelimStr() << std::endl;
	}
  return 0;
}


int repelinRunner::runTRF(const njh::progutils::CmdArgs & inputCommands){
	SimpleTRFinderLocsPars pars;
	bfs::path genomicLocation;
	pars.maxRepeatUnitSize = 6;
	uint32_t match{2};
	uint32_t mismatch{7};
	uint32_t delta{7};
	uint32_t PM{80};
	uint32_t PI{10};
	uint32_t Minscore{50}; //with a match of 2 and a min score of 50, the smallest repeat size would be 25bps
	uint32_t MaxPeriod{1000};

	bool supplement = false;

  seqSetUp setUp(inputCommands);

  setUp.processVerbose();
  setUp.setOption(supplement, "--supplement", "supplement TRF output with simple repeat determination which can sometimes be missed by TRF");
  setUp.processReadInNames({"--fasta", "--fastq", "--fastagz", "--fastqgz"});
  //TRF options
  setUp.setOption(match,     "--match", "match");
  setUp.setOption(mismatch,  "--mismatch", "mismatch");
  setUp.setOption(delta,     "--delta", "delta");
  setUp.setOption(PM,        "--PM", "PM");
  setUp.setOption(PI,        "--PI", "PI");
  setUp.setOption(Minscore,  "--Minscore", "Minscore");
  setUp.setOption(MaxPeriod, "--MaxPeriod", "MaxPeriod");
  setUp.setOption(genomicLocation, "--genomicLocation", "genomic Location of the template sequence to adjust to create genomic bed file location of ouput");



  setUp.processDirectoryOutputName(true);
  //supplemental options
  pars.verbose = setUp.pars_.verbose_;
	pars.setDefaultOpts(setUp);
  setUp.finishSetUp();

  setUp.startARunLog(setUp.pars_.directoryName_);
  pars.outOpts = OutOptions(njh::files::make_path(setUp.pars_.directoryName_, njh::pasteAsStr("repeats", pars.minRepeatUnitSize, "-", pars.maxRepeatUnitSize, ".bed")));

  njh::sys::requireExternalProgramThrow("trf");

  GenomicRegion gRegion;
  if("" != genomicLocation){
  	auto regions = bedPtrsToGenomicRegs(getBeds(genomicLocation));
  	gRegion = regions.front();//just use the very first region
  }

	{
		SeqInput reader(setUp.pars_.ioOptions_);
		auto fastaInputOut = SeqIOOptions::genFastaOut(
				njh::files::make_path(setUp.pars_.directoryName_, "input.fasta"));
		SeqOutput writer(fastaInputOut);
		seqInfo seq;
		reader.openIn();
		writer.openOut();

		while (reader.readNextRead(seq)) {
			writer.write(seq);
		}
		writer.closeOut();

		std::stringstream cmdStream;
		cmdStream << "cd " << setUp.pars_.directoryName_ << " && " << "trf " << "input.fasta "
				<< " " << match
				<< " " << mismatch
				<< " " << delta
				<< " " << PM
				<< " " << PI
				<< " " << Minscore
				<< " " << MaxPeriod
				<< " -f -d -m -h > trf.log 2>&1";
		auto cmdOut = njh::sys::run(VecStr{cmdStream.str()});
	  OutOptions trfOutRunLogOpts{njh::files::make_path(setUp.pars_.directoryName_, "trfOutRunLog.json")};
	  OutputStream trfOutRunLogOut(trfOutRunLogOpts);

	  trfOutRunLogOut << cmdOut.toJson() << std::endl;
	}

	std::stringstream repeatFnpStream;
	repeatFnpStream << "input.fasta"
			<< "." << match
			<< "." << mismatch
			<< "." << delta
			<< "." << PM
			<< "." << PI
			<< "." << Minscore
			<< "." << MaxPeriod << ".dat";


  OutOptions outOpts{njh::files::make_path(setUp.pars_.directoryName_, "repeats.bed")};
  OutOptions outGenomicOpts{njh::files::make_path(setUp.pars_.directoryName_, "repeats_genomic.bed")};

	{
	  bfs::path repeatFilename = njh::files::make_path(setUp.pars_.directoryName_, repeatFnpStream.str());
	  OutputStream outFile(outOpts);
	  std::unique_ptr<OutputStream> outGenomicLocationOut;
	  if("" != genomicLocation){
	  	outGenomicLocationOut = std::make_unique<OutputStream>(outGenomicOpts);
	  }
		BioDataFileIO<TandemRepeatFinderRecord> reader{IoOptions(InOptions(repeatFilename))};
		reader.openIn();
		std::string currentSeqName = "";
		std::string line = "";
		while (njh::files::crossPlatGetline(*reader.inFile_, line)) {
			if (njh::beginsWith(line, "Sequence")) {
				auto toks = njh::tokenizeString(line, "whitespace");
				if (toks.size() < 2) {
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error in processing line " << line
							<< ", was expecting at least two values separated by white sapce but found "
							<< toks.size() << " instead " << "\n";
					throw std::runtime_error { ss.str() };
				}
				currentSeqName = toks[1];
			}
			//skip lines that don't start with a digit or is a blank line, it seems to be the only indicator that the file is on a data line, this format is stupid
			if("" == line || allWhiteSpaceStr(line) || !::isdigit(line.front())){
				continue;
			}
			std::shared_ptr<TandemRepeatFinderRecord> record = std::make_shared<TandemRepeatFinderRecord>(line);
			if(record->numberOfAlignedRepeats_ < pars.minNumRepeats ||
					record->fullSeq_.size() < pars.lengthCutOff){
				continue;
			}
			record->setSeqName(currentSeqName);
			outFile << GenomicRegion(*record).genBedRecordCore().toDelimStr() << std::endl;
			if(outGenomicLocationOut){
				auto relativeToTemplateLoc = GenomicRegion(*record).genBedRecordCore();
				auto outRegion = relativeToTemplateLoc;
				outRegion.chrom_ = gRegion.chrom_;
				outRegion.strand_ = (gRegion.reverseSrand_ ? '-' : '+');
				outRegion.chromStart_ = gRegion.getRelativePositionFromStartStrandAware(relativeToTemplateLoc.chromStart_);
				outRegion.chromEnd_ = gRegion.getRelativePositionFromStartStrandAware(relativeToTemplateLoc.chromEnd_);
				*outGenomicLocationOut << outRegion.toDelimStrWithExtra() << std::endl;
			}
		}
	}

	if(supplement){
		runSimpleTRFinderLocsPars(setUp.pars_.ioOptions_, pars);
		//combine

		std::vector<Bed6RecordCore> combinedRepeats;
		{
			BioDataFileIO<Bed6RecordCore> trfRepeatsReader{IoOptions(InOptions(outOpts.outFilename_))};
			trfRepeatsReader.openIn();
			Bed6RecordCore rec;
			while(trfRepeatsReader.readNextRecord(rec)){
				bool found = false;
				for(const auto & other : combinedRepeats){
					if(other.sameRegion(rec) && other.name_ == rec.name_){
						found = true;
						break;
					}
				}
				if(!found){
					combinedRepeats.emplace_back(rec);
				}
			}
		}
		{
			BioDataFileIO<Bed6RecordCore> strRepeatsReader{IoOptions(InOptions(pars.outOpts.outFilename_))};
			strRepeatsReader.openIn();
			Bed6RecordCore rec;
			while(strRepeatsReader.readNextRecord(rec)){
				bool found = false;
				for(const auto & other : combinedRepeats){
					if(other.sameRegion(rec) && other.name_ == rec.name_){
						found = true;
						break;
					}
				}
				if(!found && rec.score_ * match >= Minscore){
					combinedRepeats.emplace_back(rec);
				}
			}
		}

		BedUtility::coordSort(combinedRepeats);

		OutputStream combinedOut(njh::files::make_path(setUp.pars_.directoryName_, "combined.bed"));

		std::unique_ptr<OutputStream> combinedGenomicOut;


		if("" != genomicLocation){
			combinedGenomicOut = std::make_unique<OutputStream>((njh::files::make_path(setUp.pars_.directoryName_, "genomic_combined.bed")));
		}
		for(const auto & reg : combinedRepeats){
			combinedOut << reg.toDelimStr() << std::endl;
			if("" != genomicLocation){
				auto outRegion = reg;
				outRegion.chrom_ = gRegion.chrom_;
				outRegion.strand_ = (gRegion.reverseSrand_ ? '-' : '+');
				if(gRegion.reverseSrand_){
					outRegion.chromStart_ = gRegion.getRelativePositionFromStartStrandAware(reg.chromEnd_) + 1;
					outRegion.chromEnd_ = gRegion.getRelativePositionFromStartStrandAware(reg.chromStart_) + 1;
				}else{
					outRegion.chromStart_ = gRegion.getRelativePositionFromStartStrandAware(reg.chromStart_);
					outRegion.chromEnd_ = gRegion.getRelativePositionFromStartStrandAware(reg.chromEnd_);
				}
				*combinedGenomicOut << outRegion.toDelimStrWithExtra() << std::endl;
			}
		}
	}

  return 0;
}




} /* namespace njhseq */
