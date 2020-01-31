/*
 * renameDownloadedGenome.cpp
 *
 *  Created on: Jan 29, 2020
 *      Author: nicholashathaway
 */

#include "miscRunner.hpp"


namespace njhseq {


int miscRunner::renameDownloadedGenome(const njh::progutils::CmdArgs & inputCommands){
	std::string genomeName = "";
	bfs::path genomeFnp = "";
	bfs::path gffFnp = "";
	std::string lower = "upper";
	bool overWrite = false;
	seqSetUp setUp(inputCommands);
	setUp.setOption(genomeName, "--name", "Name to give genome", true);
	setUp.setOption(gffFnp, "--gff", "GFF file", true);
	setUp.setOption(genomeFnp, "--genome", "Genome fasta file", true);
	setUp.setOption(overWrite, "--overWrite", "Over Write");
	setUp.setOption(lower, "--lower", "How to handle lower case in the genome");

	setUp.finishSetUp(std::cout);

	OutOptions genomeOutOps(bfs::path(njh::pasteAsStr(genomeName, ".fasta")));
	genomeOutOps.overWriteFile_ = overWrite;
	OutputStream genomeOut(genomeOutOps);

	InOptions gffInOpts(gffFnp);
	OutOptions gffOutOps(bfs::path(njh::pasteAsStr(genomeName, ".gff")));
	gffOutOps.overWriteFile_ = overWrite;
	BioDataFileIO<GFFCore> gffReader(IoOptions(gffInOpts, gffOutOps));
	gffReader.openIn();
	gffReader.openOut();


	OutOptions nameKeyOutOpts(bfs::path(njh::pasteAsStr(genomeName, "_renameKey.txt")));
	nameKeyOutOpts.overWriteFile_ = overWrite;


	OutputStream nameKeyOut(nameKeyOutOpts);


	std::map<std::string, std::string> renameKey;

	SeqIOOptions genomeInOpts(genomeFnp, SeqIOOptions::getInFormatFromFnp(genomeFnp), false);
	SeqInput genomeReader(genomeInOpts);
	genomeReader.openIn();
	seqInfo seq;

	std::regex chromosomePat{".*chromosome: ([0-9]+).*"};
	std::regex contigPat{".*contig: ([0-9_A-z]+).*"};
	std::regex organellePlastidPat{".*organelle: plastid:([0-9_A-z]+).*"};
	std::regex organellePat{".*organelle: ([0-9_A-z]+).*"};

	while(genomeReader.readNextRead(seq)){
		std::string name = seq.name_;
		trimAtFirstWhitespace(name);
		if(njh::in(name, renameKey)){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "can't have multiple of the same name in genome, repeat: " << name<< "\n";
			throw std::runtime_error{ss.str()};
		}
		std::string newName = name;
		if(seq.name_.find("chromosome:") != std::string::npos){
			std::smatch match;
			if(std::regex_match(seq.name_, match, chromosomePat)){
				auto chromNum = njh::StrToNumConverter::stoToNum<uint32_t>(match[1]);
				newName = njh::pasteAsStr(genomeName, "_", njh::leftPadNumStr(chromNum, 10U));
			}else{
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "couldn't match chromosome pattern to " << seq.name_<< "\n";
				throw std::runtime_error{ss.str()};
			}
		}else if(seq.name_.find("contig:") != std::string::npos){
			std::smatch match;
			if(std::regex_match(seq.name_, match, contigPat)){
				newName = match[1];
			}else{
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "couldn't match contig pattern to " << seq.name_<< "\n";
				throw std::runtime_error{ss.str()};
			}
		}else if(seq.name_.find("organelle: plastid:") != std::string::npos){
			std::smatch match;
			if(std::regex_match(seq.name_, match, organellePlastidPat)){
				newName = njh::pasteAsStr(genomeName, "_", match[1]);
			}else{
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "couldn't match organelle plastid pattern to " << seq.name_<< "\n";
				throw std::runtime_error{ss.str()};
			}
		}else if(seq.name_.find("organelle:") != std::string::npos){
			std::smatch match;
			if(std::regex_match(seq.name_, match, organellePat)){
				newName = njh::pasteAsStr(genomeName, "_", match[1]);
			}else{
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "couldn't match organelle pattern to " << seq.name_<< "\n";
				throw std::runtime_error{ss.str()};
			}
		}else{
			//unhandled case, no renaming;

		}
		//update key
		renameKey[name] = newName;

		//lower case handling
		readVec::handelLowerCaseBases(seq, lower);
		seq.name_ = njh::pasteAsStr(newName, " ", seq.name_);
		seq.outPutSeq(genomeOut);
	}


	//write rename key
	nameKeyOut << "oldName\tnewName" << std::endl;
	for(const auto & nameKey : renameKey){
		nameKeyOut << nameKey.first << "\t" << nameKey.second << std::endl;
	}

	GFFCore record;
	while(gffReader.readNextRecord(record)){
		if(!njh::in(record.seqid_, renameKey)){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "couldn't find " << record.seqid_ << " in rename key"<< "\n";
			throw std::runtime_error{ss.str()};
		}
		record.seqid_ = renameKey[record.seqid_];
		record.writeGffRecord(*gffReader.out_);
	}
	return 0;
}






}  // namespace njhseq
