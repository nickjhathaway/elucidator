
//  seqUtilsConverterRunner.cpp
//
//  Created by Nicholas Hathaway on 2015/05/28.
//  Copyright (c) 2015 Nicholas Hathaway. All rights reserved.
//

    
#include "seqUtilsConverterRunner.hpp"
    
    
namespace njhseq {

seqUtilsConverterRunner::seqUtilsConverterRunner()
    : njh::progutils::ProgramRunner({
	addFunc("convertFiles", convertFiles, false),
  addFunc("convertTxtToFasta", convertTxtToFasta, false),
  addFunc("convertFastaToTxt", convertFastaToTxt, false),
  addFunc("convertTxtOneLineToFasta", convertTxtOneLineToFasta, false),
	addFunc("sffInfo", sffInfo, false),
  addFunc("plasmoDBTxtToFasta", plasmoDBTxtToFasta, false)},
                    "seqUtilsConverter") {}

int seqUtilsConverterRunner::plasmoDBTxtToFasta(const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	IoOptions opts;
	setUp.setOption(opts.in_.inFilename_, "--in", "Input file name", true);
	setUp.setOption(opts.out_.overWriteFile_,"--overWrite",  "Overwrite output file");
  setUp.finishSetUp(std::cout);
  table inTab(opts.in_.inFilename_, "whitespace", false);

  std::unordered_map<std::string, std::string> seqs;
  for(const auto & row : inTab.content_){
  	for(const auto & colPos : iter::range(len(row))){
  		if(colPos <=1){
  			continue;
  		}
  		seqs[row[0]].append(row[colPos]);
  	}
  }
  std::vector<readObject> outSeqs;
  for(const auto & s : seqs){
  	if(s.first == ""){
  		continue;
  	}
  	outSeqs.emplace_back(readObject(seqInfo(s.first, njh::replaceString(s.second, ".", ""))));
  }
  SeqIOOptions options;
  options.out_.outFilename_ = bfs::basename(opts.in_.inFilename_);
  options.outFormat_ = SeqIOOptions::outFormats::FASTA;
  options.out_.overWriteFile_ = opts.out_.overWriteFile_;
  SeqOutput::write(outSeqs, options);

  return 0;
}
                    
int seqUtilsConverterRunner::convertTxtOneLineToFasta(const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	std::string filename = "";
	setUp.pars_.ioOptions_.out_.outFilename_ = "out";
	setUp.setOption(filename, "--file,-f", "filename", true);
	setUp.processWritingOptions();
	setUp.finishSetUp(std::cout);
	table inTab(filename, "whitespace", false);
	std::vector<readObject> reads;
	for(const auto & row : inTab.content_){
		reads.emplace_back(readObject(seqInfo(row[0], row[1])));
	}
	auto options = setUp.pars_.ioOptions_;
  options.outFormat_ = SeqIOOptions::outFormats::FASTA;
	SeqOutput::write(reads, options);
	return 0;
}



int seqUtilsConverterRunner::convertFiles(const njh::progutils::CmdArgs & inputCommands) {
  seqUtilsConverterSetUp setUp(inputCommands);
  setUp.setUpConvertFiles();
  // read in the sequences
  SeqIO reader(setUp.pars_.ioOptions_);
  reader.openIn();
  reader.openOut();
  seqInfo read;
  while(reader.readNextRead(read)){
  	reader.write(read);
  }
  setUp.logRunTime(std::cout);
  return 0;
}


int seqUtilsConverterRunner::convertTxtToFasta(const njh::progutils::CmdArgs & inputCommands) {
  seqUtilsConverterSetUp setUp(inputCommands);
  std::string fileName = "";
  std::string outFileName = "";
  std::string stub = "Seq";
  bool overWrite = false;
  setUp.setOption(fileName, "-file", "Filename", true);
  if (!setUp.setOption(outFileName, "-out", "Oufilename")) {
    outFileName = bfs::basename(fileName) + ".fasta";
  }
  setUp.setOption(overWrite, "-overWrite", "overWriteCurrentFile");
  setUp.finishSetUp(std::cout);
  VecStr dnaStrings;
  table inTab(fileName);
  std::ofstream outFile;
  openTextFile(outFile, outFileName, ".fasta", overWrite, false);
  int count = 0;
  for (const auto &line : inTab.content_) {
    outFile << ">" << stub << "." << count << std::endl;
    outFile << line[0] << std::endl;
    ++count;
  }
  return 0;
}

int seqUtilsConverterRunner::convertFastaToTxt(const njh::progutils::CmdArgs & inputCommands) {
  bool addingName = false;
  bool addNameAfter = false;
  OutOptions outOpts;
  outOpts.outExtention_ = ".txt";
  seqUtilsConverterSetUp setUp(inputCommands);
  setUp.setOption(addingName, "--addName", "Whether adding the name to the output or not");
  setUp.setOption(addNameAfter, "--addNameAfter", "When adding name, add it after the sequence");

  setUp.pars_.ioOptions_.out_.outFilename_ = "";
  setUp.processReadInNames(true);
  setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);
	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	OutputStream out(outOpts);
	seqInfo seq;
	while(reader.readNextRead(seq)){
		if(addingName && !addNameAfter){
			out << seq.name_ << "\t";
		}
		out << seq.seq_;
		if(addingName && addNameAfter){
			out << "\t"<< seq.name_;
		}
		out << std::endl;
	}
  return 0;//multiBamCoverageFinder
}

int seqUtilsConverterRunner::sffInfo(const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	std::string filename = "";
	uint32_t testLimit = std::numeric_limits<uint32_t>::max();

	setUp.setOption(testLimit, "-testLimit,-limit", "testLimit");
	setUp.setOption(setUp.pars_.ioOptions_.firstName_, "-file,-filename", "filename", true);
	setUp.pars_.ioOptions_.inFormat_ = SeqIOOptions::inFormats::SFFBIN;
	bfs::path outFilename = filename + ".txt";
	setUp.setOption(outFilename, "-out", "Out filename");
	setUp.finishSetUp(std::cout);

	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();



	outFilename = njh::files::findNonexitantFile(outFilename);

	std::ofstream outSfftxt(outFilename.string());
	reader.sffBinHeader_->printSffTxtStyle(outSfftxt);
	seqInfo read;
	uint32_t count = 0;
	while (reader.readNextRead(read)) {
		if(count > testLimit){
			break;
		}

		bool okay = reader.lastSffRead_->sanityCheck();
		if (!okay) {
			break;
		}
		reader.lastSffRead_->printHeaderSffTxtStyle(outSfftxt);
		reader.lastSffRead_->printSffTxtSeqData(outSfftxt);
		++count;
		//report progress
		if((count + 1) % 10000 == 0){
			std::cout << count + 1 << std::endl;
		}
		if (count >= reader.sffBinHeader_->numReads) {
			break;
		}
	}

	return 0;
}
                    
} // namespace njhseq
