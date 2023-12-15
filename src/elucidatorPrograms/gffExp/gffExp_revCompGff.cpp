/*
 * gffExp_revCompGff.cpp
 *
 *  Created on: Dec 17, 2020
 *      Author: nick
 */





#include "elucidatorPrograms/gffExp/gffExp.hpp"
#include "elucidator/objects/BioDataObject.h"
#include "elucidator/seqToolsUtils/seqToolsUtils.hpp"

#include <njhseq/objects/Gene.h>
#include <TwoBit.h>




namespace njhseq {



int gffExpRunner::revCompGff(const njh::progutils::CmdArgs & inputCommands){
	//currently just turns phase into NA since it's hard to determine what rev complementing will do to the phase without knowing all genes
	bfs::path inputFile;
	bfs::path genomeTwobit;
	OutOptions outOpts(bfs::path("out.gff"));
	seqSetUp setUp(inputCommands);
	setUp.setOption(inputFile, "--gff", "Input gff file", true);
	setUp.setOption(genomeTwobit, "--genome2bit", "2bit genome file", true);


	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);


	OutputStream outFile(outOpts);

	BioDataFileIO<GFFCore> reader{(IoOptions(InOptions(inputFile)))};
	reader.openIn();
	// uint32_t count = 0;
	std::string line = "";
	TwoBit::TwoBitFile tReader(genomeTwobit);
	auto chromLens = tReader.getSeqLens();

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
	while(nullptr != gRecord) {
		gRecord->revCompRecord(chromLens[gRecord->seqid_]);
		gRecord->writeGffRecord(outFile);
		bool end = false;
		while ('#' == reader.inFile_->peek()) {
			if (njh::files::nextLineBeginsWith(*reader.inFile_, "##FASTA")) {
				//write out the fasta if there
				while('#' == reader.inFile_->peek()){
					njh::files::crossPlatGetline(*reader.inFile_, line);
					outFile << line << "\n";
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


