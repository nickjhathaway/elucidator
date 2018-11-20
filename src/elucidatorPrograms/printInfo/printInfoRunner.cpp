
//  printInfoRunner.cpp
//
//  Created by Nicholas Hathaway on 2015/05/29.
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

    
#include "printInfoRunner.hpp"
    
    
namespace njhseq {

printInfoRunner::printInfoRunner()
    : njh::progutils::ProgramRunner({
	addFunc("printFastqAscII", printFastqAscII, false),
  addFunc("printDegen", printDegen, false),
  addFunc("printAminoAcidInfo", printAminoAcidInfo, false)},
                    "printInfo") {}
                    
int printInfoRunner::printFastqAscII(const njh::progutils::CmdArgs & inputCommands){
	seqSetUp setUp(inputCommands);
	bool illuminaOffset = false;
	setUp.setOption(illuminaOffset, "--illumina", "Use illumina offset");
	setUp.finishSetUp(std::cout);
	table outTable{VecStr{"Char", "Qual"}};
	for(const auto & pos : iter::range(93)){
		if(illuminaOffset){
			outTable.content_.emplace_back(toVecStr(static_cast<char>(pos + IlluminaQualOffset), pos));
		}else{
			outTable.content_.emplace_back(toVecStr(static_cast<char>(pos + SangerQualOffset), pos));
		}
	}
	outTable.outPutContentOrganized(std::cout);
	return 0;
}

int printInfoRunner::printAminoAcidInfo(const njh::progutils::CmdArgs & inputCommands){
	seqSetUp setUp(inputCommands);
	bool colorOutput = false;
	setUp.setOption(colorOutput, "-color", "colorOutput");
	setUp.finishSetUp(std::cout);
	table outInfo(VecStr{"dnaCodons", "rnaCodons", "numCodes", "letCode", "triCode",
		"fullName", "classification", "weight", "acidHydrophobicity"});

	for(const auto & aa : aminoAcidInfo::infos::allInfo){
		outInfo.content_.emplace_back(aa.second.getInfo());
	}
	outInfo.sortTable("letCode", false);
	std::stringstream infoStream;
	outInfo.outPutContentOrganized(infoStream);
	VecStr lines = streamToVecStr(infoStream);
	VecStr classificationVec = outInfo.getColumn("classification");
	for(const auto & pos : iter::range(len(outInfo.content_))){
		if(pos > 0 && colorOutput){
			njh::color currentColor = njh::color(aminoAcidInfo::infos::aaClassColorCode.at(classificationVec[pos -1]));
					std::cout <<
					njh::bashCT::addBGColor(getClosetAnsiColor(currentColor).first)
					<< lines[pos] << njh::bashCT::reset << std::endl;
		}else{
			std::cout << lines[pos] << std::endl;
		}
	}
	return 0;
}

int printInfoRunner::printDegen(const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	setUp.finishSetUp(std::cout);
	substituteMatrix mat = substituteMatrix::createDegenScoreMatrix(1, 0);
	std::unordered_map<char, std::vector<char>> acceptable;
	for (const auto & pos : iter::range(len(mat.mat_))) {
		for (const auto & subPos : iter::range(len(mat.mat_[pos]))) {
			if (pos == subPos) {
				continue;
			}
			if (pos == 'N' || subPos == 'N') {
				continue;
			}
			if (njh::in<char>(pos, std::vector<char> { 'A', 'C', 'G', 'T' })) {
				continue;
			}
			if (!njh::in<char>(subPos, std::vector<char> { 'A', 'C', 'G', 'T' })) {
				continue;
			}
			if (mat.mat_[pos][subPos] > 0) {
				acceptable[pos].emplace_back(subPos);
			}
		}
	}
	std::vector<char> keys = getVectorOfMapKeys(acceptable);
	njh::sort(keys);
	struct ownComp {
		bool operator()(const std::string & str1, const std::string & str2) const {
			if (str1.length() > str2.length()) {
				return false;
			} else if (str1.length() == str2.length()) {
				if (str1 < str2) {
					return true;
				} else {
					return false;
				}
			} else {
				return true;
			}
		}
	};
	std::multimap<std::string, char, ownComp> byBases;
	for (const auto & base : keys) {
		std::stringstream tempStream;
		tempStream << vectorToString(acceptable[base], ",");
		byBases.insert( { tempStream.str(), base });
	}
	std::unordered_map<char, std::string> descriptions = { { 'W',
			njh::bashCT::boldBlack("W") + "eak" }, { 'S', njh::bashCT::boldBlack("S")
			+ "trong" }, { 'M', "a" + njh::bashCT::boldBlack("M") + "ino" }, { 'K',
			njh::bashCT::boldBlack("K") + "eto" }, { 'R', "pu"
			+ njh::bashCT::boldBlack("R") + "ine" }, { 'Y', "p"
			+ njh::bashCT::boldBlack("Y") + "rimidine" }, { 'B', "not A ("
			+ njh::bashCT::boldBlack("B") + " comes after A)" }, { 'D', "not C ("
			+ njh::bashCT::boldBlack("D") + " comes after C)" }, { 'H', "not G ("
			+ njh::bashCT::boldBlack("H") + " comes after G)" }, { 'V', "not T ("
			+ njh::bashCT::boldBlack("V") + " comes after T and U)" } };

	table outTab(VecStr { "Symbol", "Bases", "Description" });
	for (const auto & final : byBases) {
		outTab.content_.emplace_back(
				VecStr { estd::to_string(final.second), final.first,
						descriptions[final.second] });
	}
	outTab.outPutContentOrganized(std::cout);
	return 0;
}
                    
} // namespace njhseq
