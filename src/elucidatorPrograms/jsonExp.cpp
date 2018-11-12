/*
 * jsonExpRunner.cpp
 *
 *  Created on: Dec 19, 2017
 *      Author: nickhathaway
 */

#include "jsonExp.hpp"
#include "elucidator/objects/BioDataObject.h"


namespace njhseq {
jsonExpRunner::jsonExpRunner()
    : njh::progutils::ProgramRunner(
          {
					 addFunc("jsonPrintFieldInArray", jsonPrintFieldInArray, false),
					 addFunc("jsonExtractInArrayMatchingField", jsonExtractInArrayMatchingField, false),
           },
          "jsonExp") {}



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
