/*
 * programWrappers_runDada2.cpp
 *
 *  Created on: Nov 26, 2019
 *      Author: nicholashathaway
 */


#include "programWrappers.hpp"


namespace njhseq {




void genDada2RScriptSingleSamplePaired(std::ostream & out,
		const std::string & inputDir,
		const std::string & outputDir,
		const std::string & filePat,
		const std::string & outName
		){
  out << "#!/usr/bin/env Rscript" << std::endl;
  out << "suppressMessages(library(dada2)); #packageVersion(\"dada2\")" << std::endl;
  out << "" << std::endl;
  out << "path <- \"" << inputDir << "\"" << std::endl;
  out << "filtpath <- \"" << outputDir << "\"" << std::endl;
  out << "outputdir <- \"" << outputDir << "\"" << std::endl;
  out << "filePattern <- \"" << filePat << "\"" << std::endl;
  out << "" << std::endl;
  out << "fns <- sort(list.files(path)) # Sort should keep them paired in order" << std::endl;
  out << "fastqs <- fns[grepl(filePattern, fns)]" << std::endl;
  out << "fnFs <- fastqs[grepl(\"_R1\", fastqs)]" << std::endl;
  out << "fnRs <- fastqs[grepl(\"_R2\", fastqs)]" << std::endl;
  out << "sam_names <- sapply(strsplit(fnFs, \"_\"), `[`, 1)" << std::endl;
  out << "" << std::endl;
  out << "" << std::endl;
  out << "filtFs <- paste0(filtpath, sapply(strsplit(fnFs, \"" << R"(\\.)"<< "\"), `[`, 1), \"_filt.fastq.gz\")" << std::endl;
  out << "filtRs <- paste0(filtpath, sapply(strsplit(fnRs, \"" << R"(\\.)"<< "\"), `[`, 1), \"_filt.fastq.gz\")" << std::endl;
  out << "for(i in seq_along(fnFs)) {" << std::endl;
  out << "  fastqPairedFilter(paste0(path, c(fnFs[i], fnRs[i])), c(filtFs[i], filtRs[i]), maxN=0, maxEE=2, truncQ=2, compress=TRUE, verbose=TRUE)" << std::endl;
  out << "}" << std::endl;
  out << "" << std::endl;
  out << "derepFs <- derepFastq(filtFs, verbose=TRUE)" << std::endl;
  out << "derepRs <- derepFastq(filtRs, verbose=TRUE)" << std::endl;
  out << "" << std::endl;
  out << "#names(derepFs) <- sam_names" << std::endl;
  out << "#names(derepRs) <- sam_names" << std::endl;
  out << "" << std::endl;
  out << "dadaFs.part <- dada(derepFs, err=inflateErr(tperr1,3), selfConsist = TRUE)" << std::endl;
  out << "dadaRs.part <- dada(derepRs, err=inflateErr(tperr1,3), selfConsist = TRUE)" << std::endl;
  out << "" << std::endl;
  out << "errF <- dadaFs.part$err_out" << std::endl;
  out << "errR <- dadaRs.part$err_out" << std::endl;
  out << "" << std::endl;
  out << "dadaFs <- dada(derepFs, err=errF, pool = TRUE, selfConsist=FALSE) # 9m" << std::endl;
  out << "dadaRs <- dada(derepRs, err=errR, pool = TRUE, selfConsist=FALSE) # 6m" << std::endl;
  out << "" << std::endl;
  out << "mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)" << std::endl;
  out << "" << std::endl;
  out << "" << std::endl;
  out << "seqtab.all <- makeSequenceTable(mergers)" << std::endl;
  out << "bim <- isBimeraDenovo(seqtab.all,minFoldParentOverAbundance = 3,  verbose=TRUE)" << std::endl;
  out << "seqtab <- seqtab.all[,!bim]" << std::endl;
  out << "" << std::endl;
  out << "chr <- function(n) { rawToChar(as.raw(n)) }" << std::endl;
  out << "writeFastqFromDada<-function(dadaSingleSample, filename, bim){" << std::endl;
  out << "  fastqLines = c()" << std::endl;
  out << "  sampleName = basename(filename)" << std::endl;
  out << "  for (seqNum in 1:length(dadaSingleSample$sequence)){" << std::endl;
  out << "    if(!missing(bim)){" << std::endl;
  out << "      if(bim[names(bim) == dadaSingleSample$sequence[seqNum]]){" << std::endl;
  out << "        fastqLines = c(fastqLines, as.character(paste(\"@CHI_\", sampleName, \".\", as.character(seqNum), \"_t\", as.character(dadaSingleSample$clustering$abundance[seqNum]), sep = \"\")))" << std::endl;
  out << "      }else{" << std::endl;
  out << "        fastqLines = c(fastqLines, as.character(paste(\"@\", sampleName, \".\", as.character(seqNum), \"_t\", as.character(dadaSingleSample$clustering$abundance[seqNum]), sep = \"\")))" << std::endl;
  out << "      }" << std::endl;
  out << "      fastqLines = c(fastqLines,  as.character(dadaSingleSample$sequence[seqNum]))" << std::endl;
  out << "      fastqLines = c(fastqLines,\"+\")" << std::endl;
  out << "      print (dadaSingleSample$quality)" << std::endl;
  out << "      fastqLines = c(fastqLines,as.character(chr(round(dadaSingleSample$quality + 33)[seqNum,1:nchar(dadaSingleSample$sequence[seqNum])])))" << std::endl;
  out << "    }else{" << std::endl;
  out << "      fastqLines = c(fastqLines, as.character(paste(\"@\", sampleName, \".\", as.character(seqNum), \"_t\", as.character(dadaSingleSample$clustering$abundance[seqNum]), sep = \"\")))" << std::endl;
  out << "      fastqLines = c(fastqLines,  as.character(dadaSingleSample$sequence[seqNum]))" << std::endl;
  out << "      fastqLines = c(fastqLines,\"+\")" << std::endl;
  out << "      fastqLines = c(fastqLines,as.character(chr(round(dadaSingleSample$quality + 33)[seqNum,1:nchar(dadaSingleSample$sequence[seqNum])])))" << std::endl;
  out << "    }" << std::endl;
  out << "  }" << std::endl;
  out << "  fileConn<-file(filename)" << std::endl;
  out << "  writeLines(fastqLines, fileConn)" << std::endl;
  out << "  close(fileConn)" << std::endl;
  out << "}" << std::endl;
  out << "writeFastaFromDadaMerged<-function(dadaSingleSample, filename, bim){" << std::endl;
  out << "  fastqLines = c()" << std::endl;
  out << "  sampleName = basename(filename)" << std::endl;
  out << "  for (seqNum in 1:length(dadaSingleSample$sequence)){" << std::endl;
  out << "    if(!missing(bim)){" << std::endl;
  out << "      if(bim[names(bim) == dadaSingleSample$sequence[seqNum]]){" << std::endl;
  out << "        fastqLines = c(fastqLines, as.character(paste(\">CHI_\", sampleName, \".\", as.character(seqNum), \"_t\", as.character(dadaSingleSample$abundance[seqNum]), sep = \"\")))" << std::endl;
  out << "      }else{" << std::endl;
  out << "        fastqLines = c(fastqLines, as.character(paste(\">\", sampleName, \".\", as.character(seqNum), \"_t\", as.character(dadaSingleSample$abundance[seqNum]), sep = \"\")))" << std::endl;
  out << "      }" << std::endl;
  out << "      fastqLines = c(fastqLines,  as.character(dadaSingleSample$sequence[seqNum]))" << std::endl;
  out << "    }else{" << std::endl;
  out << "      fastqLines = c(fastqLines, as.character(paste(\">\", sampleName, \".\", as.character(seqNum), \"_t\", as.character(dadaSingleSample$abundance[seqNum]), sep = \"\")))" << std::endl;
  out << "      fastqLines = c(fastqLines,  as.character(dadaSingleSample$sequence[seqNum]))" << std::endl;
  out << "    }" << std::endl;
  out << "  }" << std::endl;
  out << "  fileConn<-file(filename)" << std::endl;
  out << "  writeLines(fastqLines, fileConn)" << std::endl;
  out << "  close(fileConn)" << std::endl;
  out << "}" << std::endl;
  out << "" << std::endl;
  out << "currentBimHaps = bim" << std::endl;
  out << "writeFastaFromDadaMerged(mergers, paste0(outputdir, \"dada2_\", \"" << outName << "\", \".fasta\"), currentBimHaps)" << std::endl;

}


void genDada2RScript(std::ostream & out,
		const std::string & inputDir,
		const std::string & outputDir,
		const std::string & filePat,
		uint32_t trimLeft,
		uint32_t trimRight,
		bool on454,
		bool checkBimeras){
  out << "#!/usr/bin/env Rscript" << std::endl;
  out << "suppressMessages(library(dada2)); #packageVersion(\"dada2\")" << std::endl;
  out << "chr <- function(n) { rawToChar(as.raw(n)) }" << std::endl;
  out << "writeFastqFromDada<-function(dadaSingleSample, filename, bim){" << std::endl;
  out << "  fastqLines = c()" << std::endl;
  out << "  sampleName = basename(filename)" << std::endl;
  out << "  for (seqNum in 1:length(dadaSingleSample$sequence)){" << std::endl;
  out << "    if(!missing(bim)){" << std::endl;
  out << "      if(bim[names(bim) == dadaSingleSample$sequence[seqNum]]){" << std::endl;
  out << "        fastqLines = c(fastqLines, as.character(paste(\"@CHI_\", sampleName, \".\", as.character(seqNum), \"_t\", as.character(dadaSingleSample$clustering$abundance[seqNum]), sep = \"\")))" << std::endl;
  out << "      }else{" << std::endl;
  out << "        fastqLines = c(fastqLines, as.character(paste(\"@\", sampleName, \".\", as.character(seqNum), \"_t\", as.character(dadaSingleSample$clustering$abundance[seqNum]), sep = \"\")))" << std::endl;
  out << "      }" << std::endl;
  out << "      fastqLines = c(fastqLines,  as.character(dadaSingleSample$sequence[seqNum]))" << std::endl;
  out << "      fastqLines = c(fastqLines,\"+\")" << std::endl;
  out << "      fastqLines = c(fastqLines,as.character(chr(round(dadaSingleSample$quality + 33)[seqNum,1:nchar(dadaSingleSample$sequence[seqNum])])))" << std::endl;
  out << "    }else{" << std::endl;
  out << "      fastqLines = c(fastqLines, as.character(paste(\"@\", sampleName, \".\", as.character(seqNum), \"_t\", as.character(dadaSingleSample$clustering$abundance[seqNum]), sep = \"\")))" << std::endl;
  out << "      fastqLines = c(fastqLines,  as.character(dadaSingleSample$sequence[seqNum]))" << std::endl;
  out << "      fastqLines = c(fastqLines,\"+\")" << std::endl;
  out << "      fastqLines = c(fastqLines,as.character(chr(round(dadaSingleSample$quality + 33)[seqNum,1:nchar(dadaSingleSample$sequence[seqNum])])))" << std::endl;
  out << "    }" << std::endl;
  out << "  }" << std::endl;
  out << "  fileConn<-file(filename)" << std::endl;
  out << "  writeLines(fastqLines, fileConn)" << std::endl;
  out << "  close(fileConn)" << std::endl;
  out << "}" << std::endl;
  out << "runDada<-function(inputPath,  outputDir, filePattern, trimLeft, trimRight, on454){" << std::endl;
  out << "  fns <- list.files(inputPath)" << std::endl;
  out << "  fastqs <- fns[grepl(filePattern, fns)]" << std::endl;
  out << "  sample_names <- sapply(strsplit(fastqs, \".\",fixed = T), `[`, 1)" << std::endl;
  out << "  fastqs <- paste0(inputPath, fastqs)" << std::endl;
  out << "  # Make filenames for the filtered fastq files" << std::endl;
  out << "  filtFs <- paste0(outputDir, sample_names, \"_filt.fastq.gz\")" << std::endl;
  out << "  # Filter" << std::endl;
  out << "  for(i in seq_along(fastqs)) {" << std::endl;
  out << "    fastqFilter(fastqs[i], filtFs[i]," << std::endl;
  out << "                trimLeft=trimLeft, truncLen=trimRight, " << std::endl;
  out << "                maxN=0, maxEE=Inf, truncQ=0, " << std::endl;
  out << "                compress=TRUE, verbose=TRUE)" << std::endl;
  out << "  }" << std::endl;
  out << "  dereps <- derepFastq(filtFs, verbose=TRUE)" << std::endl;
  out << "  # Name the derep-class objects by the sample names" << std::endl;
  out << "  if(length(filtFs) > 1){" << std::endl;
  out << "    names(dereps) <- sample_names" << std::endl;
  out << "  }" << std::endl;
  out << "  if(on454){" << std::endl;
  out << "    dadaFsBoth <- dada(dereps, err=inflateErr(tperr1,3), selfConsist = TRUE, HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32)" << std::endl;
  out << "  }else{" << std::endl;
  out << "    dadaFsBoth <- dada(dereps, err=inflateErr(tperr1,3), selfConsist = TRUE)" << std::endl;
  out << "  }" << std::endl;
  out << "  return(dadaFsBoth)" << std::endl;
  out << "}" << std::endl;
  out << "path <- \"" << inputDir << "\" " << std::endl;
  out << "filePat <- \"" << filePat << "\"" << std::endl;
  out << "outputdir <- \"" << outputDir << "\"" << std::endl;
  out << "trimLeft <- " << trimLeft << " " << std::endl;
  out << "trimRight <- " << trimRight << "" << std::endl;
  out << "on454 <- " << (on454 ? 'T' : 'F') << "" << std::endl;
  out << "checkBimeras <- " << (checkBimeras ? 'T' : 'F') << "" << std::endl;
  out << "output = runDada(path,  outputdir, filePat, trimLeft, trimRight, on454)" << std::endl;
  out << "if (checkBimeras) {" << std::endl;
  out << "  seqtab.all <- makeSequenceTable(output)" << std::endl;
  out << "  bim <- isBimeraDenovo(seqtab.all, verbose=TRUE)" << std::endl;
  out << "  if(class(output)[1] == \"list\"){" << std::endl;
  out << "    for (samp in names(output)){" << std::endl;
  out << "      currentSampHaps = seqtab.all[rownames(seqtab.all) == samp, ]" << std::endl;
  out << "      currentBimHaps = bim[(1:length(currentSampHaps))[currentSampHaps > 0]]" << std::endl;
  out << "      writeFastqFromDada(output[[samp]], paste0(outputdir, \"dada2_\", samp, \".fastq\"), currentBimHaps)" << std::endl;
  out << "    }" << std::endl;
  out << "  }else{" << std::endl;
  out << "    fns <- list.files(path)" << std::endl;
  out << "    fastqs <- fns[grepl(filePat, fns)]" << std::endl;
  out << "    sample_names <- sapply(strsplit(fastqs, \".\",fixed = T), `[`, 1)" << std::endl;
  out << "    currentSampHaps = seqtab.all[1, ]" << std::endl;
  out << "    currentBimHaps = bim[(1:length(currentSampHaps))[currentSampHaps > 0]]" << std::endl;
  out << "    writeFastqFromDada(output, paste0(outputdir, \"dada2_\", sample_names[1], \".fastq\"), currentBimHaps)" << std::endl;
  out << "  }" << std::endl;
  out << "" << std::endl;
  out << "} else {" << std::endl;
  out << "  if(class(output)[1] == \"list\"){" << std::endl;
  out << "    for (samp in names(output)){" << std::endl;
  out << "      writeFastqFromDada(output[[samp]], paste0(outputdir, \"dada2_\", samp, \".fastq\"))" << std::endl;
  out << "    }" << std::endl;
  out << "  }else{" << std::endl;
  out << "    fns <- list.files(path)" << std::endl;
  out << "    fastqs <- fns[grepl(filePat, fns)]" << std::endl;
  out << "    sample_names <- sapply(strsplit(fastqs, \".\",fixed = T), `[`, 1)" << std::endl;
  out << "    writeFastqFromDada(output, paste0(outputdir, \"dada2_\", sample_names[1], \".fastq\"))" << std::endl;
  out << "  }" << std::endl;
  out << "}" << std::endl;
}




int programWrapperRunner::runDada2SingleSamplePaired(const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	std::string inputDir = "";
	std::string outputDir = "";
	std::string filePat = "[0-9]+.fastq.gz$";
	std::string additionalOut = "";
	bool overWriteDir = false;
	setUp.setOption(inputDir,     "--inputDir",     "Name of the input directory", true);
	setUp.setOption(outputDir,    "--outputDir",    "Name of the output directory", true);
	setUp.setOption(overWriteDir, "--overWriteDir", "Whether to Overwrite the --outputDir");
	setUp.setOption(filePat,      "--filePat",      "The file pattern to run dada2 on in the --inputDir", true);
	setUp.setOption(additionalOut,"--additionalOut","Additional out location file");
	setUp.processVerbose();
	setUp.finishSetUp(std::cout);
	if (!bfs::exists(inputDir)) {
		std::stringstream ss;
		ss << "Input directory " << njh::bashCT::bold << inputDir
				<< njh::bashCT::boldRed(" doesn't exist") << std::endl;
		throw std::runtime_error { ss.str() };
	}
	if (bfs::exists(outputDir) && !overWriteDir) {
		std::stringstream ss;
		ss << "Output directory " << njh::bashCT::bold << outputDir
				<< njh::bashCT::boldRed(" already exists ") << std::endl;
		ss << "Use " << njh::bashCT::bold << "--overWriteDir" << njh::bashCT::reset
				<< " to overwrite" << std::endl;
		throw std::runtime_error { ss.str() };
	}
	njh::files::makeDir(njh::files::MkdirPar(outputDir, overWriteDir));
	outputDir = bfs::canonical(outputDir).string();
	inputDir = bfs::canonical(inputDir).string();
	njh::appendAsNeeded(outputDir, "/");
	njh::appendAsNeeded(inputDir, "/");
	setUp.startARunLog(outputDir);


	{
		auto dada2Script = njh::files::join(outputDir, "runDada2.R");
		std::ofstream dataScript(dada2Script.string());
		::chmod(dada2Script.c_str(),
		S_IWUSR | S_IRUSR | S_IXUSR | S_IRGRP | S_IWGRP | S_IXGRP | S_IROTH);
		genDada2RScriptSingleSamplePaired(dataScript, inputDir, outputDir, filePat, filePat);
	}
	std::string runDada2 = "cd " + outputDir + " && ./runDada2.R ";

	auto runOut = njh::sys::run({runDada2});

	std::ofstream dataLogFile(njh::files::join(outputDir, "dada2Log.json").string());
	dataLogFile << runOut.toJson() << std::endl;
	if ("" != additionalOut) {
		auto files = njh::files::listAllFiles(outputDir, false, { std::regex {R"(dada2_.*\.fasta)"} });
		//std::cout << "files.size(): " << files.size() << std::endl;
		for (const auto & file : files) {
			auto opts = SeqIOOptions::genFastaIn(file.first.string(), true);
			std::string additionalOutDir = findAdditonalOutLocation(additionalOut,
					opts.firstName_.string());
			opts.out_.outFilename_ = njh::files::join(additionalOutDir, "output.fasta");
			if(bfs::exists(opts.out_.outFilename_)){
				bfs::remove(opts.out_.outFilename_);
			}
			bfs::copy_file(opts.firstName_, opts.out_.outFilename_);
		}
	}

	return 0;
}


int programWrapperRunner::runDada2(const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	std::string inputDir = "";
	std::string outputDir = "";
	std::string filePat = "[0-9]+.fastq.gz$";
	std::string additionalOut = "";
	uint32_t trimLeft = 0;
	uint32_t trimRight = 0;
	bool trimToMin = false;
	bool on454 = false;
	bool checkBimeras = false;
	bool overWriteDir = false;
	setUp.setOption(trimToMin,    "--trimToMin",    "Trim to the minimum length of the input");
	setUp.setOption(inputDir,     "--inputDir",     "Name of the input directory", true);
	setUp.setOption(outputDir,    "--outputDir",    "Name of the output directory", true);
	setUp.setOption(overWriteDir, "--overWriteDir", "Whether to Overwrite the --outputDir");
	setUp.setOption(filePat,      "--filePat",      "The file pattern to run dada2 on in the --inputDir");
	setUp.setOption(trimLeft,     "--trimLeft",     "Trim left position",true);
	setUp.setOption(trimRight,    "--trimRight",    "Trim Right position", !trimToMin);
	setUp.setOption(on454,        "--on454",        "Whether the data is IonTorrent or 454");
	setUp.setOption(checkBimeras, "--checkBimeras", "Whether to check for bimeras or not");
	setUp.setOption(additionalOut,"--additionalOut","Additional out location file");
	setUp.processVerbose();
	setUp.finishSetUp(std::cout);
	if (!bfs::exists(inputDir)) {
		std::stringstream ss;
		ss << "Input directory " << njh::bashCT::bold << inputDir
				<< njh::bashCT::boldRed(" doesn't exist") << std::endl;
		throw std::runtime_error { ss.str() };
	}
	if (bfs::exists(outputDir)) {
		if(!overWriteDir){
			std::stringstream ss;
			ss << "Output directory " << njh::bashCT::bold << outputDir
					<< njh::bashCT::boldRed(" already exists ") << std::endl;
			ss << "Use " << njh::bashCT::bold << "--overWriteDir" << njh::bashCT::reset
					<< " to overwrite" << std::endl;
			throw std::runtime_error { ss.str() };
		}else{
			njh::files::rmDirForce(outputDir);
		}
	}
	njh::files::makeDir(njh::files::MkdirPar(outputDir, overWriteDir));
	outputDir = bfs::canonical(outputDir).string();
	inputDir = bfs::canonical(inputDir).string();
	njh::appendAsNeeded(outputDir, "/");
	njh::appendAsNeeded(inputDir, "/");
	setUp.startARunLog(outputDir);
	if(trimToMin){
		auto files = njh::files::listAllFiles(inputDir, false, {std::regex(filePat)});
		uint64_t minLen = std::numeric_limits<uint64_t>::max();
		seqInfo seq;
		for(const auto & f : files){
			auto inputOpts = SeqIOOptions(f.first.string(), SeqIOOptions::getInFormat(njh::files::getExtension(f.first.string())), false);
			SeqInput reader(inputOpts);
			reader.openIn();
			while(reader.readNextRead(seq)){
				readVec::getMinLength(seq, minLen);
			}
		}
		trimRight = minLen;
	}

	{
		auto dada2Script = njh::files::join(outputDir, "runDada2.R");
		std::ofstream dataScript(dada2Script.string());
		::chmod(dada2Script.c_str(),
		S_IWUSR | S_IRUSR | S_IXUSR | S_IRGRP | S_IWGRP | S_IXGRP | S_IROTH);
		genDada2RScript(dataScript, inputDir, outputDir, filePat, trimLeft, trimRight, on454, checkBimeras);
	}
	//std::string runDada2 = "cd " + outputDir + " && /usr/bin/time -o dada2Time.txt --portability ./runDada2.R ";
	std::string runDada2 = "cd " + outputDir + " && ./runDada2.R ";

	auto runOut = njh::sys::run({runDada2});

	std::ofstream dataLogFile(njh::files::join(outputDir, "dada2Log.json").string());
	dataLogFile << runOut.toJson() << std::endl;
	if ("" != additionalOut) {
		auto files = njh::files::listAllFiles(outputDir, false, { std::regex {R"(dada2_.*\.fastq)"} });
		//std::cout << "files.size(): " << files.size() << std::endl;
		for (const auto & file : files) {
			auto opts = SeqIOOptions::genFastqIn(file.first.string(), true);
			std::string additionalOutDir = findAdditonalOutLocation(additionalOut,
					opts.firstName_.filename().string());
			if ("" == additionalOutDir) {
				std::cerr << njh::bashCT::red << njh::bashCT::bold;
				std::cerr << "No additional out directory found for: "
						<< opts.firstName_.filename() << std::endl;
				std::cerr << njh::bashCT::reset;
				table inTab(additionalOut, "\t");
			  MapStrStr additionalOutNames;
			  for (const auto& fIter : inTab.content_) {
			    additionalOutNames[makeIDNameComparable(fIter[0])] = fIter[1];
			  }
			  std::cerr << "Options are: " << njh::conToStr(getVectorOfMapKeys(additionalOutNames), ", ") << std::endl;

			} else {
				opts.out_.outFilename_ = njh::files::join(additionalOutDir, "output.fastq");
				if(bfs::exists(opts.out_.outFilename_)){
					bfs::remove(opts.out_.outFilename_);
				}
				bfs::copy_file(opts.firstName_, opts.out_.outFilename_);
			}
		}
	}

	return 0;
}




}  // namespace njhseq




