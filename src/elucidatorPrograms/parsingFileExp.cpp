/*
 * parsingFileExpRunner.cpp
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


#include "parsingFileExp.hpp"


namespace njhseq {
parsingFileExpRunner::parsingFileExpRunner()
    : njh::progutils::ProgramRunner(
          {
					 addFunc("parsingBioSampleSetNCBIJson", parsingBioSampleSetNCBIJson, false),
					 addFunc("getAttributeLevelsBioSampleSetNCBIJson", getAttributeLevelsBioSampleSetNCBIJson, false),
					 addFunc("getStrainInfoTableBioSampleSetNCBIJson", getStrainInfoTableBioSampleSetNCBIJson, false),
					 addFunc("parseBlastpHitsTab", parseBlastpHitsTab, false),
					 addFunc("parseSTOCKHOLM", parseSTOCKHOLM, false),
					 addFunc("parseSTOCKHOLMToFasta", parseSTOCKHOLMToFasta, false),
					 addFunc("parsehmmerDomainHitTab", parsehmmerDomainHitTab, false),
           },
          "parsingFileExp") {}
//

class BioSampleSetNCBIJsonParse {
public:
	BioSampleSetNCBIJsonParse(const bfs::path & fnp) :
			fnp_(fnp), info_(njh::json::parseFile(fnp_.string())) {
	}

	const bfs::path fnp_;

	const Json::Value info_;

	std::set<std::string> attributeLevels_;

	void setAttributeLevelsAll(){
		attributeLevels_.clear();
		const auto & bioSampleList = info_["BioSampleSet"]["BioSample"];
		for(const auto & sample : bioSampleList){
			auto levels = getAttrLevels(sample);
			attributeLevels_.insert(levels.begin(), levels.end());
		}
	}

	static VecStr getAttrLevels(const Json::Value & sample){
		VecStr ret;
		if(sample.isMember("Attributes") && sample["Attributes"].isMember("Attribute") && sample["Attributes"]["Attribute"].isArray()){
			for(const auto & attr : sample["Attributes"]["Attribute"]){
				if(attr.isMember("attribute_name")){
					ret.emplace_back(attr["attribute_name"].asString());
				}
			}
		}
		return ret;
	}

	static std::string getLevelForAttrName(const Json::Value & sample, const std::string & attrName){
		std::string ret = "NA";
		if(sample.isMember("Attributes") && sample["Attributes"].isMember("Attribute") && sample["Attributes"]["Attribute"].isArray()){
			for(const auto & attr : sample["Attributes"]["Attribute"]){
				if(attr.isMember("attribute_name")){
					if(attrName == attr["attribute_name"].asString()){
						ret = attr["attribute_name"].asString();
					}
				}
			}
		}
		return ret;
	};

	static std::string getMemeberValue(const Json::Value & sample,
			const std::string & member) {
		std::string ret = "NA";
		if(sample.isMember(member)){
			ret = sample[member].asString();
		}
		return ret;
	}

	static std::string getMemeberValueNested1(const Json::Value & sample,
			const std::string & member, const std::string & nestedMember) {
		std::string ret = "NA";
		if(sample.isMember(member) && sample[member].isMember(nestedMember)){
			ret = sample[member][nestedMember].asString();
		}
		return ret;
	}

	static std::string getAttributeValues(const Json::Value & sample,
			std::function<std::string(const Json::Value &)> extractorFunc) {

		std::string ret = "NA";

		if (sample.isMember("Attributes")
				&& sample["Attributes"].isMember("Attribute")
				&& sample["Attributes"]["Attribute"].isArray()) {
			ret = extractorFunc(sample["Attributes"]["Attribute"]);
		}
		return ret;
	}

};

int parsingFileExpRunner::getStrainInfoTableBioSampleSetNCBIJson(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path jsonFile = "";
	auto outOpts = TableIOOpts::genTabFileOut("out.tab.txt", true);
	seqSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processVerbose();
	setUp.processWritingOptions(outOpts.out_);
	setUp.setOption(jsonFile, "--jsonFile",
				"A JSON file convert by xml2json from xml downloaded from NCBI search", true);
	setUp.finishSetUp(std::cout);

	BioSampleSetNCBIJsonParse parser(jsonFile);
	const auto & bioSampleList = parser.info_["BioSampleSet"]["BioSample"];
	table outTable(VecStr{"id", "accession", "Title", "strain"});

	for(const auto & sample : bioSampleList){
		outTable.addRow(
									 BioSampleSetNCBIJsonParse::getMemeberValue(sample, "id"),
									 BioSampleSetNCBIJsonParse::getMemeberValue(sample, "accession"),
									 BioSampleSetNCBIJsonParse::getMemeberValueNested1(sample, "Description", "Title"),
									 BioSampleSetNCBIJsonParse::getAttributeValues(sample, [](const Json::Value & val){
			std::string ret = "NA";
			for(const auto & attr : val){
				if(attr.isMember("attribute_name") && "strain" == njh::strToLowerRet(attr["attribute_name"].asString())){
					if(attr.isMember("$t")){
						if("NA" == ret){
							ret = attr["$t"].asString();
						}else{
							ret.append(",");
							ret.append(attr["$t"].asString());
						}
					}
				}
			}
			return ret;
		}));
	}
	outTable.outPutContents(outOpts);
	return 0;
}

int parsingFileExpRunner::parsingBioSampleSetNCBIJson(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path jsonFile = "";
	auto outOpts = TableIOOpts::genTabFileOut("out.tab.txt", true);
	seqSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processVerbose();
	setUp.processWritingOptions(outOpts.out_);
	setUp.setOption(jsonFile, "--jsonFile",
				"A JSON file convert by xml2json from xml downloaded from NCBI search", true);
	setUp.finishSetUp(std::cout);

	BioSampleSetNCBIJsonParse parser(jsonFile);
	const auto & bioSampleList = parser.info_["BioSampleSet"]["BioSample"];
	uint32_t noAccession = 0;
	uint32_t noId = 0;
	uint32_t noDescription = 0;
	for(const auto & sample : bioSampleList){
		if(sample.isMember("accession")){
			std::cout << "accession" << std::endl;
			std::cout << sample["accession"] << std::endl;
		}else{
			++	noAccession;
		}
		if(sample.isMember("id")){
			std::cout << "id" << std::endl;
			std::cout << sample["id"] << std::endl;
		}else{
			++	noId;
		}
		if(sample.isMember("Description")){
			std::cout << "Description" << std::endl;
			std::cout << sample["Description"] << std::endl;
		}else{
			++	noDescription;
		}
	}
	std::cout << "noId: " << noId << std::endl;
	std::cout << "noAccession: " << noAccession << std::endl;
	std::cout << "noDescription: " << noDescription << std::endl;
	return 0;
}



int parsingFileExpRunner::getAttributeLevelsBioSampleSetNCBIJson(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path jsonFile = "";
	auto outOpts = TableIOOpts::genTabFileOut("out.tab.txt", true);
	seqSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processVerbose();
	setUp.processWritingOptions(outOpts.out_);
	setUp.setOption(jsonFile, "--jsonFile",
				"A JSON file convert by xml2json from xml downloaded from NCBI search", true);
	setUp.finishSetUp(std::cout);

	BioSampleSetNCBIJsonParse parser(jsonFile);
	parser.setAttributeLevelsAll();
	std::cout << njh::conToStr(parser.attributeLevels_, "\n") << std::endl;
	return 0;
}

int parsingFileExpRunner::parseBlastpHitsTab(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path jsonFile = "";
	bfs::path hitFnp = "";
	IoOptions ioOpts;
	ioOpts.out_.outExtention_ = ".json";
	seqSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processVerbose();
	setUp.processWritingOptions(ioOpts.out_);
	setUp.setOption(ioOpts.in_.inFilename_, "--hitFnp", "output from blastp filename", true);
	setUp.finishSetUp(std::cout);

	BioDataFileIO<BLASTHitTab> reader(ioOpts);
	reader.openIn();
	reader.openOut();
	BLASTHitTab hit;
	while(reader.readNextRecord(hit)){
		(*reader.out_) << hit.toJson() << std::endl;
	}

	return 0;
}


} /* namespace njhseq */
