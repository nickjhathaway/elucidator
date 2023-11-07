/*
 * jsonExpRunner.cpp
 *
 *  Created on: Dec 19, 2017
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

#include "jsonExp.hpp"
#include "elucidator/objects/BioDataObject.h"


namespace njhseq {
jsonExpRunner::jsonExpRunner()
    : njh::progutils::ProgramRunner(
          {
					 addFunc("jsonPrintFieldInArray", jsonPrintFieldInArray, false),
					 addFunc("jsonExtractInArrayMatchingField", jsonExtractInArrayMatchingField, false),
           addFunc("jsonPrintFirstLevelNames", jsonPrintFirstLevelNames, false),
           },
          "jsonExp") {}


int jsonExpRunner::jsonPrintFirstLevelNames(const njh::progutils::CmdArgs & inputCommands) {
  bfs::path jsonFnp = "";
  seqSetUp setUp(inputCommands);
  setUp.setOption(jsonFnp, "--jsonFnp", "JSON file to read from", true);
  setUp.finishSetUp(std::cout);

  auto inputValues = njh::json::parseFile(jsonFnp.string());
  std::cout << njh::conToStr(inputValues.getMemberNames(), ",") << std::endl;

  return 0;
}

int jsonExpRunner::jsonPrintFieldInArray(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path jsonFnp = "";
	std::string field = "";
	OutOptions outOpts;
	seqSetUp setUp(inputCommands);
	setUp.setOption(jsonFnp, "--jsonFnp", "JSON file to read from", true);
	setUp.setOption(field, "--field", "field to print", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	auto inputValues = njh::json::parseFile(jsonFnp.string());
	if(!inputValues.isArray()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error was expecting " << jsonFnp << " to be an array" << "\n";
		throw std::runtime_error{ss.str()};
	}
	OutputStream out(outOpts);
	for(const auto & val : inputValues){
		if(!val.isMember(field)){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error val didn't have " << field << "\n";
			ss << "fields: " << njh::conToStr(val.getMemberNames(), ", ") << "\n";
			throw std::runtime_error{ss.str()};
		}
		out << val[field].asString() << std::endl;
	}

	return 0;
}

int jsonExpRunner::jsonExtractInArrayMatchingField(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path jsonFnp = "";
	std::string field = "";
	std::string valuesStr = "";
	OutOptions outOpts;
	outOpts.outExtention_ = ".json";
	seqSetUp setUp(inputCommands);
	setUp.setOption(jsonFnp, "--jsonFnp", "JSON file to read from", true);
	setUp.setOption(field, "--field", "field to print", true);
	setUp.setOption(valuesStr, "--values", "Values that field has to match to be extracted", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	auto inputValues = njh::json::parseFile(jsonFnp.string());
	if(!inputValues.isArray()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error was expecting " << jsonFnp << " to be an array" << "\n";
		throw std::runtime_error{ss.str()};
	}
	auto matchValues = getInputValues(valuesStr, ",");
	OutputStream out(outOpts);
	Json::Value outVals;
	for(const auto & val : inputValues){
		if(!val.isMember(field)){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error val didn't have " << field << "\n";
			ss << "fields: " << njh::conToStr(val.getMemberNames(), ", ") << "\n";
			throw std::runtime_error{ss.str()};
		}
		if(njh::in(val[field].asString(), matchValues)){
			outVals.append(val);
		}
	}
	out << outVals << std::endl;

	return 0;
}




} /* namespace njhseq */
