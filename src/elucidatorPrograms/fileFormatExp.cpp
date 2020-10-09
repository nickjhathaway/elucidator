/*
 * fileFormatExpRunner.cpp
 *
 *  Created on: May 18, 2015
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


#include "fileFormatExp.hpp"
#include "elucidator/objects/BioDataObject.h"


namespace njhseq {
fileFormatExpRunner::fileFormatExpRunner()
    : njh::progutils::ProgramRunner(
          {
					 addFunc("extractRefSeqRecords", extractRefSeqRecords, false),
					 addFunc("parsePf3kEmblFilesToGff3", parsePf3kEmblFilesToGff3, false)
           },//
          "fileFormatExp") {}



class EmblEntry {
public:
	EmblEntry(const std::string & idLine){
		if(idLine.size() < 2 || "ID" != idLine.substr(0,2)){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error in processing line " << idLine << "\n";
			ss << "Line should start with ID" << "\n";
			throw std::runtime_error{ss.str()};
		}

		auto toks = tokenizeString(idLine.substr(2), ";");
		if(7 != toks.size()){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error in processing line " << idLine << "\n";
			ss << "ID line should be 7 values seperated by semicolons" << "\n";
			throw std::runtime_error{ss.str()};
		}
		acc_ = trimEndWhiteSpaceReturn(toks[0]);
		sequence_version_number_ = trimEndWhiteSpaceReturn(toks[1]);
		topology_ = trimEndWhiteSpaceReturn(toks[2]);
		molecular_type_ = trimEndWhiteSpaceReturn(toks[3]);
		data_class_ = trimEndWhiteSpaceReturn(toks[4]);
		taxonomic_division_ = trimEndWhiteSpaceReturn(toks[5]);
		toks[6] = trimEndWhiteSpaceReturn(toks[6]);
		trimAtFirstWhitespace(toks[6]);
		sequence_length_ = njh::StrToNumConverter::stoToNum<uint32_t>(toks[6]);

	}

	std::string acc_;
	std::string sequence_version_number_;
	std::string topology_;
	std::string molecular_type_;
	std::string data_class_;
	std::string taxonomic_division_;
	uint32_t sequence_length_{std::numeric_limits<uint32_t>::max()};

	void addInfo(const std::string & line){
		if(line.size() < 2){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error all lines should at least contain two characters" << "\n";
			throw std::runtime_error{ss.str()};
		}
		std::string begin = line.substr(0,2);
		std::string rest = line.substr(2);
		if ("AC" == begin) {

		} else if ("PR" == begin) {

		} else if ("DT" == begin) {

		} else if ("DE" == begin) {

		} else if ("KW" == begin) {

		} else if ("OS" == begin) {

		} else if ("OC" == begin) {

		} else if ("OG" == begin) {

		} else if ("RN" == begin) {

		} else if ("RC" == begin) {

		} else if ("RP" == begin) {

		} else if ("RX" == begin) {

		} else if ("RG" == begin) {

		} else if ("RA" == begin) {

		} else if ("RT" == begin) {

		} else if ("RL" == begin) {

		} else if ("DR" == begin) {

		} else if ("CC" == begin) {

		} else if ("AS" == begin) {

		} else if ("FH" == begin) {

		} else if ("FT" == begin) {

		} else if ("SQ" == begin) {

		} else if ("CO" == begin) {

		} else if ("bb" == begin) {

		} else {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error, unhandled line: " << begin << "\n";
			throw std::runtime_error{ss.str()};
		}
	}

};

class EMBLFeature{
public:
	EMBLFeature(VecStr  lines){
		if(lines.empty()){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", input vector shouldn't be empty" << "\n";
			throw std::runtime_error{ss.str()};
		}
		auto firstLineToks = tokenizeString(lines.front(), "whitespace");
		if(3 != firstLineToks.size()){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", first line should be three items separated by whitespace, not: " << lines.front() << "\n";
			throw std::runtime_error{ss.str()};
		}
		type_ = firstLineToks[1];
		location_str_ = firstLineToks[2];
		if(lines.size() > 1){
			//process lines
			for(const auto pos : iter::range<uint32_t>(1, lines.size())){
				if(lines[pos].size() < 3){
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", qualifier lines should be greater than 3 characters: " << lines[pos] << "\n";
					throw std::runtime_error{ss.str()};
				}
				auto firstNonSpace = lines[pos].find_first_not_of(" ", 2);
				if(std::string::npos == firstNonSpace){
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", qualifier lines should have information after spaces: " << lines[pos] << "\n";
					throw std::runtime_error{ss.str()};
				}
				lines[pos] = lines[pos].substr(firstNonSpace);
			}
			VecStr trueLines;
			for(const auto pos : iter::range<uint32_t>(1, lines.size())){
				if(lines[pos].front() == '/'){
					trueLines.emplace_back(lines[pos]);
				}else{
					if(trueLines.empty()){
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ", first qualifier line should start with /: " << lines[pos] << "\n";
						throw std::runtime_error{ss.str()};
					}
					trueLines.back().append(lines[pos]);
				}
			}
			for(const auto & tline : trueLines){
				if(tline.size() < 1 || '/' != tline.front() ){
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", qualifier line should start with /: " << tline << "\n";
					throw std::runtime_error{ss.str()};
				}
				auto equalSignPos = tline.find("=");
				if(std::string::npos == equalSignPos){
					std::string front = tline.substr(1, equalSignPos-1);
					qualifiers_[front] = "";
				}else{
					std::string front = tline.substr(1, equalSignPos-1);
					std::string rest = tline.substr(equalSignPos + 1);
					if(rest.size() < 1){
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ", qualifier line should have data have = : " << tline << ", rest: " << rest << "\n";
						throw std::runtime_error{ss.str()};
					}
					if(rest.size() > 2 && '"' == rest.front()){
						if('"' != rest.back()){
							std::stringstream ss;
							ss << __PRETTY_FUNCTION__ << ", qualifier lines that start with a quote need to end with a quote = : " << tline << ", rest: " << rest << "\n";
							throw std::runtime_error{ss.str()};
						}
						rest = rest.substr(1, rest.size() - 2);
					}
					qualifiers_[front] = rest;
				}
			}
		}
	}

	std::string type_;
	std::string location_str_;
	std::unordered_map<std::string, std::string> qualifiers_;

	Json::Value toJson() const{
		Json::Value ret;
		ret["class"] = njh::getTypeName(*this);
		ret["type_"] = njh::json::toJson(type_);
		ret["location_str_"] = njh::json::toJson(location_str_);
		ret["qualifiers_"] = njh::json::toJson(qualifiers_);
		return ret;
	}
};

int fileFormatExpRunner::parsePf3kEmblFilesToGff3(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path emblInput = "";
	OutOptions gffOutOpts;
	gffOutOpts.outExtention_ = ".gff3";
	seqSetUp setUp(inputCommands);
	setUp.processWritingOptions(gffOutOpts);
	setUp.setOption(emblInput, "--embl", "The embl file supplied by Pf3k", true);
	setUp.finishSetUp(std::cout);
	std::regex locPat{R"(([0-9]+)(\.\.)([0-9]+))"};

	InputStream in(emblInput);
	OutputStream gffOut(gffOutOpts);
	gffOut << "##gff-version   3" << std::endl;
	std::string line = "";
	std::vector<std::shared_ptr<EmblEntry>> emblEnteries;
	std::shared_ptr<EmblEntry> entery;
	while(njh::files::crossPlatGetline(in, line)){
		if(njh::beginsWith(line, "XX")){
			continue;
		}
		if(line.size() < 2){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error all lines should at least contain two characters" << "\n";
			throw std::runtime_error{ss.str()};
		}
		if(njh::beginsWith(line, "ID")){
			if(nullptr != entery){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", new ID line found by the last entery was not terminated" << "\n";
				throw std::runtime_error{ss.str()};
			}
			entery = std::make_shared<EmblEntry>(line);
		}else if(njh::beginsWith(line, "//")){
			emblEnteries.emplace_back(entery);
			entery = nullptr;
		}else{
			if(nullptr == entery){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", new entery not started yet but next line is: " << line << "\n";
				throw std::runtime_error{ss.str()};
			}else{
				if(njh::beginsWith(line, "FH   Key             Location/Qualifiers")){
					if(!njh::files::crossPlatGetline(in, line)){
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ", there should be a blank line following the beginning of the feature table" << "\n";
						throw std::runtime_error{ss.str()};
					}
					if("FH" != line){
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ", there should be a blank line following the beginning of the feature table" << "\n";
						throw std::runtime_error{ss.str()};
					}
					VecStr currentFeature;
					std::vector<std::shared_ptr<EMBLFeature>> features;
					while(njh::files::crossPlatGetline(in, line)){

						if(!njh::beginsWith(line, "FT")){
							if(!currentFeature.empty()){
								features.emplace_back(std::make_shared<EMBLFeature>(currentFeature));
								currentFeature.clear();
							}
							break;
						}
						if(line.size() > 6){
							if(' ' != line[6]){
								if(!currentFeature.empty()){
									features.emplace_back(std::make_shared<EMBLFeature>(currentFeature));
									currentFeature.clear();
								}
							}
							currentFeature.emplace_back(line);
						}else{
							std::stringstream ss;
							ss << __PRETTY_FUNCTION__ << "FT lines should be more than 6 characters" << "\n";
							throw std::runtime_error{ss.str()};
						}
					}
					std::unordered_map<std::string, uint32_t> counts;
					std::unordered_map<std::string, uint32_t> locusTagNamesCounts;
					struct ChromeRecord{
						ChromeRecord(
								const std::string & name,
								uint32_t start,
								uint32_t end): name_(name), start_(start), end_(end){

						}
						std::string name_;
						uint32_t start_;
						uint32_t end_;
					};
					std::vector<ChromeRecord> chromRecs;
					for(const auto & feature : features){
						if(("source" == feature->type_ || "fasta_record" == feature->type_) && njh::in("origid", feature->qualifiers_)){
							std::smatch match;
							bool matched = std::regex_match(feature->location_str_, match, locPat);
							if(!matched){
								std::stringstream ss;
								ss << __PRETTY_FUNCTION__ << "chrome records need to have positions created START..END only" << "\n";
								ss << njh::json::toJson(feature) << "\n";
								throw std::runtime_error{ss.str()};

							}
							chromRecs.emplace_back(feature->qualifiers_["origid"],
									njh::StrToNumConverter::stoToNum<uint32_t>(match[1]),
									njh::StrToNumConverter::stoToNum<uint32_t>(match[3]) );
						}
						njh::sort(chromRecs, [](const ChromeRecord & rec1, const ChromeRecord & rec2){
							return rec1.start_ < rec2.start_;
						});
						if(chromRecs.size() > 1){
							for(const auto pos : iter::range<uint32_t>(1, chromRecs.size())){
								if(chromRecs[pos-1].end_ >= chromRecs[pos].start_){
									std::stringstream ss;
									ss << __PRETTY_FUNCTION__ << " positions shouldn't be overlapping : " << "chromRecs[pos-1].end_: " << chromRecs[pos-1].end_ <<
											" positions[pos].start_: " << chromRecs[pos].start_
											<< std::endl
											<< njh::json::toJson(feature) << "\n";

									throw std::runtime_error{ss.str()};
								}else if(chromRecs[pos-1].end_ + 1 != chromRecs[pos].start_){
									std::stringstream ss;
									ss << __PRETTY_FUNCTION__ << " chrome positions should be one plus the end of the previous chromosome : " << "chromRecs[pos-1].end_ + 1: " << chromRecs[pos-1].end_ + 1 <<
											" positions[pos].start_: " << chromRecs[pos].start_
											<< std::endl
											<< njh::json::toJson(feature) << "\n";
									throw std::runtime_error{ss.str()};
								}
							}
						}
					}
					std::unordered_map<std::string, ChromeRecord> chromosomes;
					for(const auto & chromRec : chromRecs){
						chromosomes.emplace(chromRec.name_, chromRec);
					}

					for(const auto & feature : features){

						if(("source" == feature->type_ || "fasta_record" == feature->type_)){
							continue;
						}
						if("misc_feature" == feature->type_){
							continue;
						}

						if(!njh::in("locus_tag", feature->qualifiers_)){
							continue;
						}

						std::string loc = feature->location_str_;
						VecStr actions;
						while(std::string::npos != loc.find("(")){
							if(')' != feature->location_str_.back()){
								std::stringstream ss;
								ss << __PRETTY_FUNCTION__ << "location string should end with ): " << feature->location_str_ << "\n";
								throw std::runtime_error{ss.str()};
							}
							std::string action = loc.substr(0, loc.find("("));
							if("join" != action && "complement" != action){
								std::stringstream ss;
								ss << __PRETTY_FUNCTION__ << "in location string, don't know how to deal with the action: " << action << "\n";
								throw std::runtime_error{ss.str()};
							}
							actions.emplace_back(action);
							loc = loc.substr(loc.find("(") + 1, loc.find(")") - 1 - loc.find("("));
						}
						//std::cout << loc << std::endl;
						auto locToks = tokenizeString(loc, ",");
						bool fail = false;
						struct Pos{
							Pos(uint32_t start, uint32_t end): start_(start), end_(end){}
							uint32_t start_;
							uint32_t end_;
						};


						std::vector<Pos> positions;
						for(const auto & tok : locToks){
							std::smatch match;
							bool matched = std::regex_match(tok, match, locPat);
//							if(feature->qualifiers_["locus_tag"] == "Pf7G8_000015300"){
//								std::cout << "Pf7G8_000015300" << std::endl;
//								std::cout << "tok: " << tok << std::endl;
//								std::cout << njh::colorBool(matched) << std::endl;
//							}
							if(!matched){
								fail = true;
								break;
							} else {
								positions.emplace_back(
										njh::StrToNumConverter::stoToNum<uint32_t>(match[1]),
										njh::StrToNumConverter::stoToNum<uint32_t>(match[3]) );
							}
						}
						if(positions.size() > 1){
							for(const auto pos : iter::range<uint32_t>(1, positions.size())){
								if(positions[pos-1].end_ >= positions[pos].start_){
//									std::stringstream ss;
//									ss << __PRETTY_FUNCTION__ << " positions shouldn't be overlapping : " << "positions[pos-1].end_: " << positions[pos-1].end_ <<
//											" positions[pos].start_: " << positions[pos].start_
//											<< std::endl
//											<< njh::json::toJson(feature) << "\n";
//
//									throw std::runtime_error{ss.str()};
									fail = true;
								}
							}
						}

						if(!fail){
							if(njh::in("locus_tag", feature->qualifiers_)){
								++locusTagNamesCounts[feature->qualifiers_["locus_tag"]];
							}
							//bool reverseStrand =
							njh::sort(positions, [](const Pos & p1, const Pos & p2){
								return p1.start_ < p2.start_;
							});

							if(!njh::in("locus_tag", feature->qualifiers_)){
								//std::stringstream ss;
								std::cout << __PRETTY_FUNCTION__ << "no locus tag : " << std::endl << njh::json::toJson(feature) << "\n";
								//throw std::runtime_error{ss.str()};
							}

							GFFCore outGff;

							outGff.start_ = positions.front().start_;
							outGff.end_ = positions.back().end_;
							for(const auto & chromRec : chromRecs){
								if(outGff.start_ >= chromRec.start_ && outGff.start_ < chromRec.end_){
									if(outGff.end_ > chromRec.end_){
										std::stringstream ss;
										ss << __PRETTY_FUNCTION__ << "feature not bounded to chromosome: " << chromRec.name_ << ": " << chromRec.start_ << "-" << chromRec.end_ << "\n";
										ss << outGff.start_ << "-" << outGff.end_;
										throw std::runtime_error{ss.str()};
									}
									outGff.seqid_ = chromRec.name_;
									outGff.start_ = outGff.start_ + 1 - chromRec.start_;
									outGff.end_ = outGff.end_ + 1 - chromRec.start_;
								}
							}
							if("." == outGff.seqid_){
								std::stringstream ss;
								ss << __PRETTY_FUNCTION__ << "feature could be determined on a chromosme: "  << "\n";
								ss << outGff.start_ << "-" << outGff.end_;
								throw std::runtime_error{ss.str()};
							}
							outGff.source_ = "pf3k";

							if( njh::in(std::string("complement"), actions)){
								outGff.strand_ = '-';
							}else{
								outGff.strand_ = '+';
							}

							if(njh::in("pseudo", feature->qualifiers_)){
								outGff.type_ = "pseudogene";
							}else{
								outGff.type_ = "gene";
							}

							for(const auto & qual : feature->qualifiers_){
								if("product" == qual.first){
									outGff.attributes_["description"] = qual.second;
								}else if("locus_tag" == qual.first){
									if(locusTagNamesCounts[qual.second]> 1){
										outGff.attributes_["ID"] = njh::pasteAsStr(qual.second, "-", locusTagNamesCounts[qual.second]);
									}else{
										outGff.attributes_["ID"] = qual.second;
									}
								}else if("" == qual.second) {
									outGff.attributes_[qual.first] = "true";
								}	else{
									outGff.attributes_[qual.first] = qual.second;
								}
							}

							outGff.writeGffRecord(gffOut);
							auto subGffRecord = outGff;
							subGffRecord.type_ = feature->type_;
							if("CDS" == feature->type_ ){
								subGffRecord.type_ = "mRNA";
							} else if("pseudogene" == outGff.type_){
								subGffRecord.type_ = "pseudogenic_transcript";
							}

							subGffRecord.attributes_.clear();
							subGffRecord.attributes_["ID"] = outGff.attributes_["ID"] + ".1";
							subGffRecord.attributes_["Parent"] = outGff.attributes_["ID"];
							subGffRecord.writeGffRecord(gffOut);

							uint32_t subCount = 1;
							for(const auto & pos : positions){
								auto subGff = outGff;
								subGff.start_ = pos.start_;
								subGff.end_ = pos.end_;
								subGff.start_ = subGff.start_ + 1 - chromosomes.at(outGff.seqid_).start_;
								subGff.end_ = subGff.end_ + 1 - chromosomes.at(outGff.seqid_).start_;
								subGff.attributes_.clear();
								subGff.type_ = "CDS";
								if("pseudogene" == outGff.type_) {
									subGff.type_ = "pseudogenic_exon";
								}
								subGff.attributes_["ID"] = njh::pasteAsStr(subGffRecord.attributes_["ID"], "-","exon", "-",('-' == outGff.strand_? positions.size() + 1 - subCount :subCount) );
								subGff.attributes_["Parent"] = subGffRecord.attributes_["ID"];
								subGff.writeGffRecord(gffOut);
								++subCount;
							}
						}
					}
					if(njh::beginsWith(line, "SQ")){
						break;
					}
				}else{
					entery->addInfo(line);
				}
			}

		}

	}

	return 0;
}


int fileFormatExpRunner::extractRefSeqRecords(const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	std::string input = "";
	std::string aliasDictFile = "";
	std::string extractNames = "";

	setUp.setOption(input, "--file", "Input File Name", true);
	setUp.setOption(aliasDictFile, "--nameAlias", "Alias name file in json");
	setUp.setOption(extractNames, "--names", "Names of genes to extract");
	setUp.finishSetUp(std::cout);
	VecStr nameVec;
	if(extractNames != ""){
		nameVec = tokenizeString(extractNames, ",");
	}
	auto records = getRefSeqRecs(input, nameVec, aliasDictFile);

	for(const auto & record : records){
		std::cout << record.second->toStrLine() << std::endl;
	}
	//std::cout << records.front().toStrLine() << std::endl;

	return 0;
}






} /* namespace njhseq */
