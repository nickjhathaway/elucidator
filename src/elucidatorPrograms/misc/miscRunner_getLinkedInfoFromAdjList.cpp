/*
 * miscRunner_getLinkedInfoFromAdjList.cpp
 *
 *  Created on: Feb 15, 2018
 *      Author: nick
 */




#include "miscRunner.hpp"
#include "elucidator/simulation.h"
#include "elucidator/BamToolsUtils.h"
#include "elucidator/objects/seqContainers.h"
#include "elucidator/objects/BioDataObject.h"



namespace njhseq {



int miscRunner::getLinkedInfoFromAdjList(const njh::progutils::CmdArgs & inputCommands){
	bfs::path inputLinks = "";
	std::string metaGrouping = "";

	OutOptions outOpts(bfs::path(""));
	//bool withId = false; /**@todo implement with id done links*/
	seqSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.setOption(inputLinks, "--inputLinks", "Name of the links file");
	setUp.setOption(metaGrouping, "--metaGrouping", "Meta Grouping to count for grouping");
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	njh::files::checkExistenceThrow(inputLinks, __PRETTY_FUNCTION__);

	InputStream in(InOptions{inputLinks});
	std::string line = "";
	OutputStream out(outOpts);
	out << "node1\tnode2\tnode1Grouping\tnode2Grouping\tnode1count\tnode2Count" << std::endl;

	while(njh::files::crossPlatGetline(in, line)){
		auto lineToks = njh::tokenizeString(line, " -- ");
		if(lineToks.size() != 2){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error in processing line " << line << "\n";
			ss << "Was expecting to to find one  -- " << "\n";
			throw std::runtime_error{ss.str()};
		}
		std::string firstNode = lineToks[0];
		std::string secondNode = lineToks[1].substr(0, lineToks[1].rfind(" "));

		std::string firstNodeGrouping;
		std::string secondNodeGrouping;
		if("" == metaGrouping){
			firstNodeGrouping = firstNode.substr(0, firstNode.rfind("."));
			secondNodeGrouping = secondNode.substr(0, secondNode.rfind("."));
		}else{
			MetaDataInName firstMeta(firstNode);
			MetaDataInName secondMeta(secondNode);
			firstNodeGrouping = firstMeta.getMeta(metaGrouping);
			secondNodeGrouping = secondMeta.getMeta(metaGrouping);
		}


		uint32_t firstNodeCount = njh::StrToNumConverter::stoToNum<uint32_t>(firstNode.substr(firstNode.rfind("_t") + 2));
		uint32_t secondNodeCount = njh::StrToNumConverter::stoToNum<uint32_t>(secondNode.substr(secondNode.rfind("_t") + 2));

		out << firstNode
				<< "\t" << secondNode
				<< "\t" << firstNodeGrouping
				<< "\t" << secondNodeGrouping
				<< "\t" << firstNodeCount
				<< "\t" << secondNodeCount
				<< std::endl;


	}

	return 0;
}

}
