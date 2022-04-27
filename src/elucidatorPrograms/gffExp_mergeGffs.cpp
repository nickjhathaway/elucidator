/*
 * gffExp_mergeGffs.cpp
 *
 *  Created on: Jun 29, 2020
 *      Author: nick
 */




#include "gffExp.hpp"
#include "elucidator/objects/BioDataObject.h"
#include "elucidator/seqToolsUtils/seqToolsUtils.hpp"

#include <njhseq/objects/Gene.h>
#include <TwoBit.h>




namespace njhseq {


int gffExpRunner::combineGffs(const njh::progutils::CmdArgs & inputCommands){
	std::vector<bfs::path> inputFiles;
	VecStr recordsToRemove{};
	OutOptions outOpts(bfs::path("out.gff"));
	seqSetUp setUp(inputCommands);
	setUp.setOption(inputFiles, "--gffs", "Input gff files", true);
	setUp.setOption(recordsToRemove, "--recordsToRemove", "records To Remove from original file if replacing chromosomes");

	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);


	std::string line = "";

	OutputStream outFile(outOpts);
	uint32_t headerCount = 0;
	for(const auto & inputFile : inputFiles){
		//write header
		std::ifstream infile(inputFile.string());
		while('#' == infile.peek()){
			njh::files::crossPlatGetline(infile, line);
			bool skip = false;
			if(line.find("##sequence-region") != std::string::npos){
				for(const auto & rec : recordsToRemove){
					if(line.find(rec) != std::string::npos){
						skip =true;
						break;
					}
				}
			}
			if(headerCount >0 && line.find("##gff-version") != std::string::npos){
				skip = true;
			}
			if(!skip){
				outFile << line << std::endl;
			}
		}
		++headerCount;
	}

	std::stringstream endOfFile;
	for(const auto & inputFile : inputFiles){

		{
			BioDataFileIO<GFFCore> reader{(IoOptions(InOptions(inputFile)))};
			reader.openIn();
			uint32_t count = 0;
			std::shared_ptr<GFFCore> gRecord = reader.readNextRecord();

			while(nullptr != gRecord) {
				if(!njh::in(gRecord->seqid_, recordsToRemove)){
					gRecord->writeGffRecord(outFile);
				}
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
				++count;
			}
		}
	}



	outFile << endOfFile.str();

	return 0;
}


int gffExpRunner::appendGff(const njh::progutils::CmdArgs & inputCommands){
	bfs::path inputFile;
	bfs::path appendingFile;
	VecStr recordsToRemove{};
	OutOptions outOpts(bfs::path("out.gff"));
	seqSetUp setUp(inputCommands);
	setUp.setOption(inputFile, "--gff", "Input gff file", true);
	setUp.setOption(appendingFile, "--appendingFile", "appending gff file", true);
	setUp.setOption(recordsToRemove, "--recordsToRemove", "records To Remove from original file if replacing chromosomes");

	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);


	std::string line = "";

	OutputStream outFile(outOpts);
	{
		//write header
		std::ifstream infile(inputFile.string());
		while('#' == infile.peek()){
			njh::files::crossPlatGetline(infile, line);
			bool skip = false;
			if(line.find("##sequence-region") != std::string::npos){
				for(const auto & rec : recordsToRemove){
					if(line.find(rec) != std::string::npos){
						skip =true;
						break;
					}
				}
			}
			if(!skip){
				outFile << line << std::endl;
			}
		}
	}
	{
		//write header
		std::ifstream infile(appendingFile.string());
		while('#' == infile.peek()){
			njh::files::crossPlatGetline(infile, line);
			bool skip = false;
			if(line.find("##gff-version") != std::string::npos){
				skip =true;
			}
			if(!skip){
				outFile << line << std::endl;
			}
		}
	}

	std::stringstream endOfFile;

	{
		BioDataFileIO<GFFCore> reader{(IoOptions(InOptions(inputFile)))};
		reader.openIn();
		uint32_t count = 0;
		std::shared_ptr<GFFCore> gRecord = reader.readNextRecord();

		while(nullptr != gRecord) {
			if(!njh::in(gRecord->seqid_, recordsToRemove)){
				gRecord->writeGffRecord(outFile);
			}
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
			++count;
		}
	}
	{
		BioDataFileIO<GFFCore> reader{(IoOptions(InOptions(appendingFile)))};
		reader.openIn();
		uint32_t count = 0;
		std::shared_ptr<GFFCore> gRecord = reader.readNextRecord();

		while(nullptr != gRecord) {
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
			++count;
		}
	}


	outFile << endOfFile.str();

	return 0;
}


}  // namespace njhseq
