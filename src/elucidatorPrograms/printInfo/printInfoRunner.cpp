
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
    	addFunc("printAminoAcidInfo", printAminoAcidInfo, false),
			addFunc("printCPPNumericalLimits", printCPPNumericalLimits, false),
    	addFunc("testingBoostFilesystem", testingBoostFilesystem, false)

    },
                    "printInfo") {}
                    //
int printInfoRunner::printFastqAscII(const njh::progutils::CmdArgs & inputCommands){
	seqSetUp setUp(inputCommands);
	bool illuminaOffset = false;
	setUp.setOption(illuminaOffset, "--illumina", "Use illumina offset");
	setUp.finishSetUp(std::cout);
	table outTable{VecStr{"Char", "Qual"}};
	for(const auto pos : iter::range(93)){
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
	infoStream.flush();
	VecStr lines = streamToVecStr(infoStream);

	VecStr classificationVec = outInfo.getColumn("classification");
	for(const auto pos : iter::range(len(outInfo.content_) +1)){
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
	for (const auto pos : iter::range(len(mat.mat_))) {
		for (const auto subPos : iter::range(len(mat.mat_[pos]))) {
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



int printInfoRunner::printCPPNumericalLimits(const njh::progutils::CmdArgs & inputCommands){
	njh::progutils::ProgramSetUp setUp(inputCommands);
	setUp.finishSetUp(std::cout);
	std::cout << "Unsigned integers (only numbers >= 0)" << std::endl;
	std::cout << "               uint8_t low: " << static_cast<uint16_t>(std::numeric_limits<uint8_t>::lowest()) << std::endl;
	std::cout << "               uint8_t min: " << static_cast<uint16_t>(std::numeric_limits<uint8_t>::min()) << std::endl;
	std::cout << "               uint8_t max: " << static_cast<uint16_t>(std::numeric_limits<uint8_t>::max()) << std::endl;
	std::cout << "            sizeof uint8_t: " << sizeof (uint8_t) << std::endl;
	std::cout << "    actual type of uint8_t: " << njh::typeStr<uint8_t>() << std::endl << std::endl;
	std::cout << "              uint16_t low: " << std::numeric_limits<uint16_t>::lowest() << std::endl;
	std::cout << "              uint16_t min: " << std::numeric_limits<uint16_t>::min() << std::endl;
	std::cout << "              uint16_t max: " << std::numeric_limits<uint16_t>::max() << std::endl;
	std::cout << "           sizeof uint16_t: " << sizeof (uint16_t) << std::endl;
	std::cout << "   actual type of uint16_t: " << njh::typeStr<uint16_t>() << std::endl << std::endl;
	std::cout << "              uint32_t low: " << std::numeric_limits<uint32_t>::lowest() << std::endl;
	std::cout << "              uint32_t min: " << std::numeric_limits<uint32_t>::min() << std::endl;
	std::cout << "              uint32_t max: " << std::numeric_limits<uint32_t>::max() << std::endl;
	std::cout << "           sizeof uint32_t: " << sizeof (uint32_t) << std::endl;
	std::cout << "   actual type of uint32_t: " << njh::typeStr<uint32_t>() << std::endl << std::endl;
	std::cout << "                  uint min: " << std::numeric_limits<uint>::lowest() << std::endl;
	std::cout << "                  uint min: " << std::numeric_limits<uint>::min() << std::endl;
	std::cout << "                  uint max: " << std::numeric_limits<uint>::max() << std::endl;
	std::cout << "               sizeof uint: " << sizeof (uint) << std::endl;
	std::cout << "       actual type of uint: " << njh::typeStr<uint>() << std::endl << std::endl;
	std::cout << "              uint64_t low: " << std::numeric_limits<uint64_t>::lowest() << std::endl;
	std::cout << "              uint64_t min: " << std::numeric_limits<uint64_t>::min() << std::endl;
	std::cout << "              uint64_t max: " << std::numeric_limits<uint64_t>::max() << std::endl;
	std::cout << "           sizeof uint64_t: " << sizeof (uint64_t) << std::endl;
	std::cout << "   actual type of uint64_t: " << njh::typeStr<uint64_t>() << std::endl << std::endl;
	std::cout << "                size_t low: " << std::numeric_limits<size_t>::lowest() << std::endl;
	std::cout << "                size_t min: " << std::numeric_limits<size_t>::min() << std::endl;
	std::cout << "                size_t max: " << std::numeric_limits<size_t>::max() << std::endl;
	std::cout << "             sizeof size_t: " << sizeof (size_t) << std::endl;
	std::cout << "     actual type of size_t: " << njh::typeStr<size_t>() << std::endl << std::endl;
	std::cout << "Signed Integers" << std::endl;
	std::cout << "                int8_t low: " << static_cast<int16_t>(std::numeric_limits<int8_t>::lowest()) << std::endl;
	std::cout << "                int8_t min: " << static_cast<int16_t>(std::numeric_limits<int8_t>::min()) << std::endl;
	std::cout << "                int8_t max: " << static_cast<int16_t>(std::numeric_limits<int8_t>::max()) << std::endl;
	std::cout << "             sizeof int8_t: " << sizeof (int8_t) << std::endl;
	std::cout << "     actual type of int8_t: " << njh::typeStr<int8_t>() << std::endl << std::endl;
	std::cout << "               int16_t low: " << std::numeric_limits<int16_t>::lowest() << std::endl;
	std::cout << "               int16_t min: " << std::numeric_limits<int16_t>::min() << std::endl;
	std::cout << "               int16_t max: " << std::numeric_limits<int16_t>::max() << std::endl;
	std::cout << "            sizeof int16_t: " << sizeof (int16_t) << std::endl;
	std::cout << "    actual type of int16_t: " << njh::typeStr<int16_t>() << std::endl << std::endl;
	std::cout << "               int32_t min: " << std::numeric_limits<int32_t>::lowest() << std::endl;
	std::cout << "               int32_t min: " << std::numeric_limits<int32_t>::min() << std::endl;
	std::cout << "               int32_t max: " << std::numeric_limits<int32_t>::max() << std::endl;
	std::cout << "            sizeof int32_t: " << sizeof (int32_t) << std::endl;
	std::cout << "    actual type of int32_t: " << njh::typeStr<int32_t>() << std::endl << std::endl;
	std::cout << "                   int min: " << std::numeric_limits<int>::lowest() << std::endl;
	std::cout << "                   int min: " << std::numeric_limits<int>::min() << std::endl;
	std::cout << "                   int max: " << std::numeric_limits<int>::max() << std::endl;
	std::cout << "                sizeof int: " << sizeof (int) << std::endl;
	std::cout << "        actual type of int: " << njh::typeStr<int>() << std::endl << std::endl;
	std::cout << "               int64_t low: " << std::numeric_limits<int64_t>::lowest() << std::endl;
	std::cout << "               int64_t min: " << std::numeric_limits<int64_t>::min() << std::endl;
	std::cout << "               int64_t max: " << std::numeric_limits<int64_t>::max() << std::endl;
	std::cout << "            sizeof int64_t: " << sizeof (int64_t) << std::endl;
	std::cout << "    actual type of int64_t: " << njh::typeStr<int64_t>() << std::endl << std::endl;
	std::cout << "Floating Point Numbers" << std::endl;
	std::cout << "                 float low: " << std::numeric_limits<float>::lowest() << std::endl;
	std::cout << "                 float min: " << std::numeric_limits<float>::min() << std::endl;
	std::cout << "                 float max: " << std::numeric_limits<float>::max() << std::endl;
	std::cout << "              sizeof float: " << sizeof (float) << std::endl;
	std::cout << "      actual type of float: " << njh::typeStr<float>() << std::endl << std::endl;
	std::cout << "                double low: " << std::numeric_limits<double>::lowest() << std::endl;
	std::cout << "                double min: " << std::numeric_limits<double>::min() << std::endl;
	std::cout << "                double max: " << std::numeric_limits<double>::max() << std::endl;
	std::cout << "             sizeof double: " << sizeof (double) << std::endl;
	std::cout << "     actual type of double: " << njh::typeStr<double>() << std::endl << std::endl;
	std::cout << "           long double low: " << std::numeric_limits<long double>::lowest() << std::endl;
	std::cout << "           long double min: " << std::numeric_limits<long double>::min() << std::endl;
	std::cout << "           long double max: " << std::numeric_limits<long double>::max() << std::endl;
	std::cout << "        sizeof long double: " << sizeof (long double) << std::endl;
	std::cout << "actual type of long double: " << njh::typeStr<long double>() << std::endl << std::endl;
	return 0;
}

int printInfoRunner::testingBoostFilesystem(const njh::progutils::CmdArgs & inputCommands){
  std::string filename;
  std::string filename2;
  njh::progutils::ProgramSetUp setUp(inputCommands);
  setUp.setOption(filename, "--file", "Filename", true);
  setUp.setOption(filename2, "--file2", "Second Filename");
  setUp.finishSetUp(std::cout);
  boost::filesystem::path testPath(filename);

  std::stringstream ss;
  ss << "returnType command returnValue" << std::endl;;
  ss << "path path " << testPath << std::endl;
  ss << "string path.filename().string() " << testPath.filename().string() << std::endl;
  ss << "string bfs::basename(path) " << bfs::basename(filename) << std::endl;
  ss << "path path.branch_path() " << testPath.branch_path() << std::endl;
  ss << "path path.extension() " << testPath.extension() << std::endl;
  ss << "path njh::files::prependFileBasename(path,string) " << njh::files::prependFileBasename(testPath, "pre_") << std::endl;
  auto testPath2 = testPath;
  ss << "path& path.replace_extension(\"\") " << testPath2.replace_extension("") << std::endl;
  ss << "string njh::files::removeExtension(string) " << njh::files::removeExtension(filename) << std::endl;
  ss << "path path.parent_path() " << testPath.parent_path()<< std::endl;
  ss << "path path.relative_path() " << testPath.relative_path() << std::endl;
  ss << "path path.root_directory() " << testPath.root_directory() << std::endl;
  ss << "path path.root_name() " << testPath.root_name() << std::endl;
  ss << "path path.root_path() " << testPath.root_path() << std::endl;
  ss << "string njh::files::getExtension(string) " << njh::files::getExtension(filename) << std::endl;
  ss << "path boost::filesystem::absolute(path) " << boost::filesystem::absolute(testPath) << std::endl;
  ss << "path njh::files::normalize(path) " << njh::files::normalize(testPath) << std::endl;
  ss << "path boost::filesystem::canonical(path) " << boost::filesystem::canonical(testPath) << std::endl;
  ss << "path testPath.lexically_relative(bfs::current_path()) " << testPath.lexically_relative(bfs::current_path()) << std::endl;
  ;
  auto testPath3 = testPath;
  ss << "path& path.replace_extension(\"other\") " << testPath3.replace_extension("other") << std::endl;
  if("" != filename2){
    bfs::path testPath4(filename2);
    //ss << "std::chrono::time_point<std::chrono::system_clock> njh::files::last_write_time(file1) " << njh::files::last_write_time(testPath).time_since_epoch() << std::endl;
    //ss << "std::chrono::time_point<std::chrono::system_clock> njh::files::last_write_time(file2) " << njh::files::last_write_time(testPath4) << std::endl;
    ss << "bool njh::files::last_write_time(file1)>njh::files::last_write_time(file2) " << njh::colorBool(njh::files::last_write_time(testPath) > njh::files::last_write_time(testPath4)) << std::endl;
    ss << "bool njh::files::firstFileIsOlder(testPath,testPath4) " << njh::colorBool(njh::files::firstFileIsOlder(testPath, testPath4)) << std::endl;
  }
  table tab = table(ss, " ", true);
  tab.outPutContentOrganized(std::cout);
  return 0;
}

} // namespace njhseq
