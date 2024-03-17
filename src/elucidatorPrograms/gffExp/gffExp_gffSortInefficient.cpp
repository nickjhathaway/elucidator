/*
 * gffExp_gffSortInefficient.cpp
 *
 *  Created on: Dec 8, 2019
 *      Author: nicholashathaway
 */


#include "elucidatorPrograms/gffExp/gffExp.hpp"
#include "elucidator/objects/BioDataObject.h"
#include "elucidator/seqToolsUtils/seqToolsUtils.hpp"

#include <njhseq/objects/Gene.h>
#include <TwoBit.h>




namespace njhseq {

//sort function
auto sortGffs = [](const std::shared_ptr<GFFCore> & p1Ptr, const std::shared_ptr<GFFCore> & p2Ptr){
	const auto & p1 = *p1Ptr;
	const auto & p2 = *p2Ptr;
	if(p1.seqid_ == p2.seqid_){
		if(p1.start_ == p2.start_){
			if(p1.end_ == p2.end_){
				if(p1.type_ == p2.type_){
					return p1.getIDAttr() < p2.getIDAttr();
				}else{
					if("mRNA" == p1.type_ ){
						return true;
					}else{
						return false;
					}
				}
			}else{
				return p1.end_ > p2.end_;
			}
		}else{
			return p1.start_ < p2.start_;
		}
	}else{
		return p1.seqid_ < p2.seqid_;
	}
};

int gffExpRunner::gffSortInefficient(const njh::progutils::CmdArgs & inputCommands){
	bfs::path inputFile;
	OutOptions outOpts(bfs::path("out.gff"));
	seqSetUp setUp(inputCommands);
	setUp.setOption(inputFile, "--gff", "Input gff file", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	BioDataFileIO<GFFCore> reader{(IoOptions(InOptions(inputFile)))};
	reader.openIn();
	// uint32_t count = 0;
	std::string line = "";

	std::ofstream outFile;
	outOpts.openFile(outFile);
	{
		//write header
		std::ifstream infile(inputFile.string());
		while('#' == infile.peek()){
			njh::files::crossPlatGetline(infile, line);
			outFile << line << std::endl;
		}
	}
	std::vector<std::shared_ptr<GFFCore>> allRecords;
	std::shared_ptr<GFFCore> gRecord = reader.readNextRecord();
	std::stringstream endOfFile;
	while(nullptr != gRecord) {
		allRecords.emplace_back(gRecord);
		bool end = false;
		while ('#' == reader.inFile_->peek()) {
			if (njh::files::nextLineBeginsWith(*reader.inFile_, "##FASTA")) {
				//write out the fasta if there
				while('#' == reader.inFile_->peek()){
					njh::files::crossPlatGetline(*reader.inFile_, line);
					endOfFile << line << "\n";
				}
				end = true;
				break;
			}
			njh::files::crossPlatGetline(*reader.inFile_, line);
		}
		if (end) {
			break;
		}
		gRecord = reader.readNextRecord();
		// ++count;
	}

	//now sort records by utilizing parents

	struct node {
		std::vector<uint32_t> children_;
		uint32_t parent_{std::numeric_limits<uint32_t>::max()};
		uint32_t group_{std::numeric_limits<uint32_t>::max()};
	};
	std::vector<node> gffNodes;


	//first create an index of ID to vector position
	std::unordered_map<std::string, uint32_t> IDIndex;
	for(const auto pos : iter::range(allRecords.size())){
		IDIndex[allRecords[pos]->getIDAttr()] = pos;
		gffNodes.emplace_back(node{});
	}
	//set children and parents

	for(const auto pos : iter::range(allRecords.size())){
		if(allRecords[pos]->hasAttr("Parent")){
			auto parentIDs = tokenizeString(allRecords[pos]->getAttr("Parent"), ",");
			for(const auto & parentID : parentIDs){
				if(!njh::in(parentID, IDIndex)){
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error " << "couldn't find index for parent: " << parentID<< "\n";
					throw std::runtime_error{ss.str()};
				}
				auto parentIDIndex = IDIndex[parentID];
				gffNodes[pos].parent_ = parentIDIndex;
				gffNodes[parentIDIndex].children_.emplace_back(pos);
			}
		}
	}




	struct GffNodeGroup{
		GffNodeGroup(const std::vector<std::shared_ptr<GFFCore>> & allRecords,
				const std::vector<uint32_t> & groupNodePositions,
				uint32_t groupID):groupID_(groupID) {
			for(const auto & pos : groupNodePositions){
				if(allRecords[pos]->hasAttr("Parent")){
					childrenNodes_.emplace_back(allRecords[pos]);
				}else{
					parentLessNodes_.emplace_back(allRecords[pos]);
				}
			}
			njh::sort(parentLessNodes_, sortGffs);
			njh::sort(childrenNodes_,   sortGffs);
			if(parentLessNodes_.empty()){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << " no parentless nodes for group: " << groupID_<< "\n";
				throw std::runtime_error{ss.str()};
			}
		}

		uint32_t groupID_;
		std::vector<std::shared_ptr<GFFCore>> parentLessNodes_;
		std::vector<std::shared_ptr<GFFCore>> childrenNodes_;
	};




	//set group

	auto addNodePosToStack = [&gffNodes](std::stack<uint32_t> & nodesToSpreadTo, uint32_t nodePos){
		if(std::numeric_limits<uint32_t>::max() == gffNodes[nodePos].group_){
			nodesToSpreadTo.push(nodePos);
		}
	};
	auto setGroupAndSpread = [&gffNodes,&addNodePosToStack](std::stack<uint32_t> & nodesToSpreadTo,
			uint32_t nodePos, uint32_t nodeGroup){
		gffNodes[nodePos].group_ = nodeGroup;
		if(std::numeric_limits<uint32_t>::max() != gffNodes[nodePos].parent_){
			addNodePosToStack(nodesToSpreadTo, gffNodes[nodePos].parent_);
		}
		for(const auto & child : gffNodes[nodePos].children_){
			addNodePosToStack(nodesToSpreadTo, child);
		}
	};
	uint32_t nodeGroup = 0;
	for(const auto nodePos : iter::range(gffNodes.size())){
		if(std::numeric_limits<uint32_t>::max() == gffNodes[nodePos].group_){
			std::stack<uint32_t> nodesToSpreadTo;
			setGroupAndSpread(nodesToSpreadTo, nodePos, nodeGroup);
			while(!nodesToSpreadTo.empty()){
				auto currentNode = nodesToSpreadTo.top();
				nodesToSpreadTo.pop();
				setGroupAndSpread(nodesToSpreadTo, currentNode, nodeGroup);
			}
			++nodeGroup;
		}
	}


	std::vector<GffNodeGroup> groups;
	std::unordered_map<uint32_t, std::vector<uint32_t>> groupPositions;
	for(const auto pos : iter::range(gffNodes.size())){
		groupPositions[gffNodes[pos].group_].emplace_back(pos);
	}

	for(const auto & groupPos : groupPositions){
		groups.emplace_back(GffNodeGroup(allRecords, groupPos.second, groupPos.first));
	}

	njh::sort(groups, [](const GffNodeGroup & g1, const GffNodeGroup & g2){
		return sortGffs(g1.parentLessNodes_.front(), g2.parentLessNodes_.front());
	});

	for(const auto & group : groups){
		for(const auto & parentLess : group.parentLessNodes_){
			parentLess->writeGffRecord(outFile);
		}
		for(const auto & child : group.childrenNodes_){
			child->writeGffRecord(outFile);
		}
	}

	outFile << endOfFile.str();

	return 0;
}





int gffExpRunner::gffRenameChroms(const njh::progutils::CmdArgs & inputCommands){
	bfs::path inputFile;
	bfs::path keyIn = "";

	OutOptions outOpts(bfs::path("out.gff"));
	seqSetUp setUp(inputCommands);
	setUp.setOption(inputFile, "--gff", "Input gff file", true);
	setUp.setOption(keyIn, "--keyIn", "A file with a key to rename seqs with");

	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);
	table nameKey(keyIn, "\t", false);
	if(2 != nameKey.nCol()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "name key has to be 2 columns no header, 1) old name, 2) new name, not size of : " << nameKey.nCol() << "\n";
		throw std::runtime_error{ss.str()};
	}
	std::unordered_map<std::string, std::string> nameKeyMap;
	for(const auto & row : nameKey){
		if(njh::in(row[0], nameKeyMap)){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "already have name: " << row[0] << "\n";
			throw std::runtime_error{ss.str()};
		}
		nameKeyMap[row[0]] = row[1];
	}
	BioDataFileIO<GFFCore> reader{(IoOptions(InOptions(inputFile)))};
	reader.openIn();
	// uint32_t count = 0;
	std::string line = "";

	std::ofstream outFile;
	outOpts.openFile(outFile);
	{
		//write header
		std::ifstream infile(inputFile.string());
		while('#' == infile.peek()){
			njh::files::crossPlatGetline(infile, line);
			if (njh::beginsWith(line, "##sequence-region ")) {
				auto lineToks = tokenizeString(line, "whitespace");
				if (lineToks.size() == 4 && njh::in(lineToks[1], nameKeyMap)) {
					lineToks[1] = nameKeyMap[lineToks[1]];
				}
				line = njh::conToStr(lineToks, " ");
			}
			outFile << line << std::endl;
		}
	}
	std::shared_ptr<GFFCore> gRecord = reader.readNextRecord();
	std::stringstream endOfFile;
	while(nullptr != gRecord) {

		if(!njh::in(gRecord->seqid_, nameKeyMap)){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "name: " << gRecord->seqid_ << " not in " << keyIn << "\n";
			throw std::runtime_error{ss.str()};
		}
		gRecord->seqid_ = nameKeyMap[gRecord->seqid_];
		gRecord->writeGffRecord(outFile);
		bool end = false;

		while ('#' == reader.inFile_->peek()) {
			if (njh::files::nextLineBeginsWith(*reader.inFile_, "##FASTA")) {
				//write out the fasta if there
				while('#' == reader.inFile_->peek()){
					njh::files::crossPlatGetline(*reader.inFile_, line);
					endOfFile << line << "\n";
				}
				end = true;
				break;
			}


			njh::files::crossPlatGetline(*reader.inFile_, line);
		}
		if (end) {
			break;
		}
		gRecord = reader.readNextRecord();
		// ++count;
	}


	return 0;
}



}  // namespace njhseq
