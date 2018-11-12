//
//  readSimulator.cpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 12/23/15.
//  Copyright (c) 2015 Nicholas Hathaway. All rights reserved.
//

#include "readSimulatorRunner.hpp"
#include "elucidator/simulation.h"
#include <SeekDeep/objects/PrimersAndMids.hpp>

namespace njhseq {


readSimulatorRunner::readSimulatorRunner()
    : njh::progutils::ProgramRunner({
	addFunc("simSingleSequencePCR", simSingleSequencePCR, false),
	addFunc("createRandomSequenceMixtures", createRandomSequenceMixtures, false),
	addFunc("simPcrShotgunSequences", simPcrShotgunSequences, false),
	addFunc("chimeraSim", chimeraSim, false),
	addFunc("createFractionAbundanceFile", createFractionAbundanceFile, false),
	addFunc("simMultipleMixturePCR", simMultipleMixturePCR, false),
	addFunc("createMinorVariantsFromKnown", createMinorVariantsFromKnown, false),
	addFunc("effectsOfPcrErrorRatePerRounds", effectsOfPcrErrorRatePerRounds, false),
	addFunc("simMultipleMixture", simMultipleMixture, false),
	addFunc("createLibrarySimMultipleMixture", createLibrarySimMultipleMixture, false),
	addFunc("createIlluminaErrorProfile", createIlluminaErrorProfile, false),
	addFunc("createLibrarySimMultipleMixtureDrugResistant", createLibrarySimMultipleMixtureDrugResistant, false),
},
                    "readSimulator") {}
//


int readSimulatorRunner::createFractionAbundanceFile(const njh::progutils::CmdArgs & inputCommands) {
  seqSetUp setUp(inputCommands);
  std::string abundanceFile = "";
  std::string barcodeFile = "";
  std::string outFile = "";
  setUp.setOption(abundanceFile, "-a,--abundanceFile", "The name of the abundance File created by createSimReads", true);
  setUp.setOption(barcodeFile, "-b,--barcodes", "Barcodes file used to create the abundanceFile", true);
  setUp.setOption(outFile, "-o,--outFile", "Name of an output file");
  setUp.finishSetUp(std::cout);
  //read in abundance table
  table abundTab(abundanceFile, "whitespace", false);
  //read in the barcodes for their names
  SeqIOOptions barInfo = SeqIOOptions::genFastaIn(barcodeFile);
  SeqInput reader(barInfo);
  reader.openIn();
  auto inReads = reader.readAllReads<readObject>();

  //get barcdoe names
  VecStr barcodeNames;
  for(const auto & read : inReads){
  	barcodeNames.emplace_back(read.seqBase_.name_);
  }

  //sum up the abundances
  std::vector<double> refSums{0};
  for(const auto & colPos : iter::range<uint64_t>(1,abundTab.columnNames_.size())){
  	uint32_t sum = 0;
  	for(const auto & rowPos : iter::range(abundTab.content_.size())){
  		sum+= estd::stou(abundTab.content_[rowPos][colPos]);
  	}
  	refSums.emplace_back(sum);
  }

  //calculate the fractional abundances
  table outTab(concatVecs(VecStr{"ref"}, barcodeNames));
  for(const auto & row : abundTab.content_){
  	VecStr rowOut{row.front()};
  	for(const auto & colPos : iter::range<uint64_t>(1,abundTab.columnNames_.size())){
  		rowOut.emplace_back(estd::to_string(estd::stou(row[colPos])/refSums[colPos]));
  	}
  	outTab.content_.emplace_back(rowOut);
  }
  //output table
  outTab.outPutContents(TableIOOpts(OutOptions(outFile,".tab.txt"), "\t", outTab.hasHeader_));
  return 0;
}




int readSimulatorRunner::createRandomSequenceMixtures(const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
  uint32_t stringLength = 300;
  uint32_t forwardPrimerSize = 20;
  uint32_t reversePrimerSize = 20;
  uint32_t libNum = 10;
  bool disparatePairs = false;
  bool gradientDisparatePairs = false;
  double startingMinorStrainPercent = 10;
  uint32_t numberOfMinorStrains = 7;
  std::string barcodesFile = "";
  std::string simName = "RSEQ";
  std::string abundancesStr = "";
  uint64_t seed = 0;
  bool useSeed = setUp.setOption(seed, "--seed", "Use this seed to generate the same exact library");
  setUp.setOption(libNum, "-libNum", "Number of libraries to sim, for sim reps, out libs will be twice this");
  setUp.setOption(abundancesStr, "-abundances", "A comma separated string indicated the minor strain freqs");
	std::vector<double> minorAbundances;
	if (abundancesStr != "") {
		minorAbundances = njh::lexical_cast_con<VecStr, std::vector<double>>(
				tokenizeString(abundancesStr, ","));
		libNum = minorAbundances.size() + 1;
	}
	setUp.setOption(forwardPrimerSize, "-forwardPrimerSize",
			"Forward Primer Size");
	setUp.setOption(reversePrimerSize, "-reversePrimerSize",
			"Reverse Primer Size");
	setUp.setOption(barcodesFile, "-barcodes", "Barcodes File", true);
	setUp.setOption(stringLength, "-len", "Length of Random Reference Sequences");
	setUp.setOption(disparatePairs, "-disparatePairs",
			"Whether to simulate one base pair difference which are distant from the major strain");
	setUp.setOption(gradientDisparatePairs, "-gradientDisparatePairs",
			"Whether to simulate minor strains with gradient base pair differences which are distant from the major strain");
	setUp.setOption(startingMinorStrainPercent, "-startingMinorStrainPercent",
			"Starting Minor Strain Percent");
	setUp.setOption(numberOfMinorStrains, "-numberOfMinorStrains",
			"Number Of Minor Strains");
	if (libNum > 36) {
		std::cerr << "Only library sizes 36 or less are supported" << std::endl;
		setUp.failed_ = true;
	}
  setUp.processDirectoryOutputName("simReads_TODAY", true);
  setUp.finishSetUp(std::cout);
  setUp.startARunLog(setUp.pars_.directoryName_);

  //sim the minor strains off the major
  MinorVariantsMixture mixture(stringLength);
  if(useSeed){
  	mixture.gen_.seedNum(seed);
  	mixture.setNewMajor(stringLength);
  }
  std::ofstream seedFile;
  openTextFile(seedFile, setUp.pars_.directoryName_ + "seed","",false, true);
  seedFile << mixture.gen_.currentSeed_;

  if(disparatePairs){
  	mixture.createDisparatePairsMixture(numberOfMinorStrains);
  }else if(gradientDisparatePairs){
  	mixture.createGradientDisparatePairsMixture(numberOfMinorStrains);
  }else{
  	mixture.createVaryingMinorVariantsMixture(numberOfMinorStrains);
  }
  //sim the primers
  std::string forward = simulation::evenRandStr(forwardPrimerSize, mixture.dnaAlphabet_ , mixture.gen_);
  std::string reverse = simulation::evenRandStr(reversePrimerSize, mixture.dnaAlphabet_ , mixture.gen_);
  std::string reverseComp = seqUtil::reverseComplement(reverse, "DNA");

  //primer file
  std::ofstream outPrimerFile;
  openTextFile(outPrimerFile, setUp.pars_.directoryName_ + "primers.fasta", ".fasta", false, false);
  outPrimerFile << ">for\n";
  outPrimerFile << forward << "\n";
  outPrimerFile << ">rev\n";
  outPrimerFile << reverseComp << "\n";
  mixture.prependPaddingSeq(forward);
  mixture.appendPaddingSeq(reverse);
  //ref file
  mixture.writeFasta(setUp.pars_.directoryName_, "refFile.fasta");

  //create abundances
  mixture.addAnEqualAbundanceMixture();
  mixture.createAbundancesForMixtures(abundancesStr);

	//dual replicates
  //id file
  SeqIOOptions barInfo = SeqIOOptions::genFastaIn(barcodesFile);
  SeqInput reader(barInfo);
  reader.openIn();
  auto inReads = reader.readAllReads<readObject>();
	//dirnames
	std::ofstream directoryNamesFile;
	openTextFile(directoryNamesFile, setUp.pars_.directoryName_ + "dualRepDirNames.txt", ".txt", true, false);
	//barcodes
	std::ofstream barcodesOutFile;
	openTextFile(barcodesOutFile, setUp.pars_.directoryName_ + "dualRepBarcodesFile.fasta", ".fasta", true, false);
	std::ofstream idFile;
	openTextFile(idFile, setUp.pars_.directoryName_ + "dualRepIdFile.txt", ".txt", true, false);
	idFile << "gene\tforwardPrimer\treversePrimer" << std::endl;
	idFile << simName << "\t" << forward << "\t" <<  reverseComp << std::endl;
	idFile << "id\tbarcode" << std::endl;
	for (const auto & barPos : iter::range(libNum * 2)) {
		if (barPos % 2 == 0) {
			directoryNamesFile << "dualReps\t" << inReads[barPos].seqBase_.name_
					<< njh::replaceString(inReads[barPos + 1].seqBase_.name_, "MID", "")
					<< "\t" << inReads[barPos].seqBase_.name_;
		} else {
			directoryNamesFile << "\t" << inReads[barPos].seqBase_.name_ << "\n";
		}
		inReads[barPos].seqBase_.outPutSeq(barcodesOutFile);
		idFile << inReads[barPos].seqBase_.name_ << "\t"
				<< inReads[barPos].seqBase_.seq_ << std::endl;
	}
	//dual abundance file
	mixture.outputAbundanceFile(setUp.pars_.directoryName_ + "dualRepAbundanceFile.tab.txt", true);

	//single replicates
	//id file
	//dirnames
	std::ofstream singleRepDirectoryNamesFile;
	openTextFile(singleRepDirectoryNamesFile, setUp.pars_.directoryName_ + "singleRepDirNames.txt", ".txt", true, false);
	//barcodes
	std::ofstream singleRepBarcodesOutFile;
	openTextFile(singleRepBarcodesOutFile, setUp.pars_.directoryName_ + "singleRepBarcodesFile.fasta", ".fasta", true, false);
	std::ofstream singleRepIdFile;
	openTextFile(singleRepIdFile, setUp.pars_.directoryName_ + "singleRepIdFile.txt", ".txt", true, false);
	singleRepIdFile << "gene\tforwardPrimer\treversePrimer" << std::endl;
	singleRepIdFile << simName << "\t" << forward << "\t" <<  reverseComp << std::endl;
	singleRepIdFile << "id\tbarcode" << std::endl;
	for(const auto & barPos : iter::range(libNum)){
		singleRepDirectoryNamesFile << "singleReps\t" << inReads[barPos].seqBase_.name_
				<< "\t" << inReads[barPos].seqBase_.name_ << "\n";
		inReads[barPos].seqBase_.outPutSeq(singleRepBarcodesOutFile);
		singleRepIdFile << inReads[barPos].seqBase_.name_ << "\t" << inReads[barPos].seqBase_.seq_ << std::endl;
	}
	//single replicate abundance file
	mixture.outputAbundanceFile(setUp.pars_.directoryName_ + "singleRepAbundanceFile.tab.txt", false);

	//output variant info
	mixture.outputVariantMutationInfofile(setUp.pars_.directoryName_ + "variantMutationsInfo.tab.txt");

	return 0;
}

int readSimulatorRunner::chimeraSim(const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	uint32_t kLength = 10;
	setUp.setOption(kLength, "-kLenth,-k", "kLength");
	setUp.processDefaultReader(true);

	setUp.finishSetUp(std::cout);

	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<readObject>();

	std::string seq1 = inReads[0].seqBase_.seq_;
	std::string seq2 = inReads[1].seqBase_.seq_;

	std::vector<uint32_t> positions;
	uint64_t maxLength = 0;
	if(seq1.size() > maxLength){
		maxLength = seq1.size();
	}
	if(seq2.size() > maxLength){
		maxLength = seq2.size();
	}
	for(const auto & pos : iter::range(maxLength - kLength)){
		if(std::equal(seq1.begin() + pos, seq1.begin() + pos + kLength, seq2.begin() + pos)){
			positions.emplace_back(pos);
		}
	}
	printVector(positions);

	return 0;
}




int readSimulatorRunner::simMultipleMixturePCR(const njh::progutils::CmdArgs & inputCommands) {
  readSimulatorSetUp setUp(inputCommands);
  double mismatchDensity = 0.25;
  uint32_t minOverLap = 10;
  uint32_t maxOverlap = 200;
	uint64_t startingTemplate = 3000;
	uint64_t finalReadAmount = 5000;
	uint32_t pcrRounds = 20;
	uint32_t initialPcrRounds = 10;
	long double errorRate = 3.5e-06;
	uint32_t numThreads = 2;
	std::string barcodesFile = "";
	std::string referenceFile = "";
	std::string abundanceFile = "";
	bool simIllumina = false;
	bool sim454 = false;
	bool keepIndividualLibReads = false;
	uint32_t pairedEndReadLength = 200;
	setUp.setOption(pairedEndReadLength, "--pairedEndReadLength", "Paired End Read Length For illumina simulation");
	setUp.setOption(keepIndividualLibReads, "--keepIndividualLibReads", "Keep the individual library fasta files which have been combined into all.fasta");
	setUp.setOption(simIllumina, "--simIllumina", "Simulate Illumina reads");
	setUp.setOption(sim454, "--sim454", "Simulate both 454");
	setUp.setOption(barcodesFile, "--barcodesFile", "Barcodes fasta File", true);
	setUp.setOption(referenceFile, "--referenceFile", "Reference sequences fasta File", true);
	setUp.setOption(abundanceFile, "--abundanceFile", "Abundance File", true);
	setUp.setOption(startingTemplate, "--startingTemplate", "Staring Template Amount");
	setUp.setOption(pcrRounds, "--pcrRounds", "Number of PCR rounds");
	setUp.setOption(initialPcrRounds, "--initialPcrRounds", "Number of Initial rounds of PCR before sampling");
	setUp.setOption(errorRate, "--errorRate", "Polymerase Error Rate");
	setUp.setOption(numThreads, "--numThreads", "Number of Threads to Use");
	setUp.setOption(finalReadAmount, "--finalReadAmount", "Final Read Amount to sample");
	setUp.processDirectoryOutputName("simPCR_TODAY", true);
	setUp.processVerbose();
	setUp.processWritingOptions();
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	setUp.writeParametersFile(setUp.pars_.directoryName_ + "parameters.tab.txt", false, true);
	njh::randomGenerator gen;
	uint64_t intErrorRate = errorRate * std::numeric_limits<uint64_t>::max();
  SeqIOOptions refInfo = SeqIOOptions::genFastaIn(referenceFile);
  SeqInput reader(refInfo);
  auto refReads = reader.readAllReads<readObject>();
  SeqIOOptions barInfo = SeqIOOptions::genFastaIn(barcodesFile);
  SeqInput barReader(barInfo);
  auto barcodes = barReader.readAllReads<readObject>();

	auto refNamesFromSeqs = readVec::getNames(refReads);
	auto libraryAbundances = sim::processLibraryAbundances(abundanceFile, refNamesFromSeqs, barcodes.size());
	VecStr libNames;
	VecStr libDirNames;
	for (const auto & lib : libraryAbundances) {
		std::vector<std::shared_ptr<seqInfo>> libReads;
		for (const auto & names : lib.second) {
			libReads.emplace_back(
					std::make_shared<seqInfo>(
							readVec::getReadByName(refReads, names.first).seqBase_));
			libReads.back()->frac_ = names.second;
		}
		auto libName = "lib"
				+ leftPadNumStr<size_t>(lib.first, libraryAbundances.size());
		sim::simLibFast(libReads, barcodes[lib.first - 1].seqBase_.seq_,
				setUp.pars_.directoryName_, libName, intErrorRate, startingTemplate,
				finalReadAmount, pcrRounds,initialPcrRounds, numThreads, setUp.pars_.verbose_);
		libNames.emplace_back(libName);
		libDirNames.emplace_back(setUp.pars_.directoryName_ + libName);
	}
	SeqIOOptions allOpts = SeqIOOptions::genFastaOut(njh::files::join(setUp.pars_.directoryName_,"all.fasta"));
	SeqIO allReaderOut(allOpts);
	allReaderOut.openOut();
	SeqIOOptions libGenOpts;
	libGenOpts.inFormat_ = SeqIOOptions::inFormats::FASTA;
	seqInfo read;
	for(const auto & lib : libNames ){
		auto libOpts = libGenOpts;
		libOpts.firstName_ = setUp.pars_.directoryName_ + lib + "/reads.fasta";
		SeqIO libReader(libOpts);
		libReader.openIn();
		while(libReader.readNextRead(read)){
			read.name_ = lib + "_" + read.name_;
			allReaderOut.write(read);
		}
	}
	if(!keepIndividualLibReads){
		for(const auto & libDir : libDirNames){
			njh::files::rmDirForce(libDir);
		}
	}
	//sim reads


	std::ofstream runLogs;
	if(simIllumina || sim454){
		openTextFile(runLogs, OutOptions(bfs::path(setUp.pars_.directoryName_ + "simProgramLogs.json")));
	}
	//sim 454
	if (sim454) {
		std::string simCmd454 = "454sim -d $(echo $(dirname $(which 454sim))/gen) "
				+ setUp.pars_.directoryName_ + "all.fasta -o " + setUp.pars_.directoryName_
				+ "all.sff";
		auto simOutPut454 = njh::sys::run(VecStr{simCmd454});
		runLogs << simOutPut454.toJson();
	}
	//sim illumina
	if (simIllumina) {
		uint32_t illuminaAttempts = 10; //art fails for no reason sometimes
		std::string simCmdIllumina = "art_illumina -amp -p -na -i "
				+ setUp.pars_.directoryName_ + "all.fasta -l " + estd::to_string(pairedEndReadLength) + " -f 1 -o "
				+ setUp.pars_.directoryName_ + "amplicon_pair_dat";
		auto simOutPutIllumina = njh::sys::run(VecStr{simCmdIllumina});
		runLogs << simOutPutIllumina.toJson();
		uint32_t numberOfAttempts = 1;
		while(!simOutPutIllumina.success_ && numberOfAttempts <= illuminaAttempts){
			simOutPutIllumina = njh::sys::run(VecStr{simCmdIllumina});
			runLogs << simOutPutIllumina.toJson();
			++numberOfAttempts;
		}
		std::string flashCmdTemplate = "flash "
				+ setUp.pars_.directoryName_ + "amplicon_pair_dat1.fq "
				+ setUp.pars_.directoryName_ + "amplicon_pair_dat2.fq "
				+ "-o " + setUp.pars_.directoryName_ + "amplicon_pair_dat "
				" --min-overlap " + estd::to_string(minOverLap) +
				" --max-overlap " + estd::to_string(maxOverlap) +
				" -x " + estd::to_string(mismatchDensity) +
				" -t " + estd::to_string(numThreads) +
				" 2>&1 | tee "
				+ setUp.pars_.directoryName_ + "flash.log";

		auto pearOutPut = njh::sys::run(VecStr { flashCmdTemplate });
		runLogs << pearOutPut.toJson();
	}
	std::cout << "Done" << std::endl;
	setUp.logRunTime(std::cout);
	return 0;
}

int readSimulatorRunner::effectsOfPcrErrorRatePerRounds(const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	uint32_t seqLen = 350;
	uint32_t seqNum = 1000;
	long double errorRate = 3.5e-06;
	uint32_t rounds = 20;
	setUp.setOption(seqNum, "--startingTemplate", "Staring Template Amount");
	setUp.setOption(seqLen, "--len", "Template Length");
	setUp.setOption(rounds, "--pcrRounds", "Number of PCR rounds");
	setUp.setOption(errorRate, "--errorRate,-e", "Error Rate");
	setUp.finishSetUp(std::cout);
	long double finalAmount = seqNum;
	long double correctRate = std::pow(1 - errorRate, seqLen);
	table roundAmounts(VecStr { "round", "Amount Correct", "percentCorrect" });
	for (uint32_t round = 1; round <= rounds; ++round) {
		finalAmount += finalAmount * correctRate;
		roundAmounts.content_.emplace_back(
				toVecStr(round, finalAmount,
						finalAmount / (seqNum * std::pow(2, round))));
	}
	roundAmounts.outPutContentOrganized(std::cout);
	std::cout << "Approximate percent of a jackpot event in first round" << std::endl;
	std::cout << std::fixed << std::pow(2, rounds - 1)/((seqNum + seqNum - 1) * std::pow(2, rounds - 1 )) << std::endl;
	return 0;
}

int readSimulatorRunner::simSingleSequencePCR(const njh::progutils::CmdArgs & inputCommands) {
  readSimulatorSetUp setUp(inputCommands);
	uint32_t len = 300;
	uint64_t startingTemplate = 3000;
	uint32_t pcrRounds = 20;
	uint32_t initialPcrRounds = 10;
	long double errorRate = 3.5e-06;
	uint32_t numThreads = 2;
	uint32_t finalReadAmount = 5000;

	setUp.setOption(startingTemplate, "--startingTemplate", "Staring Template Amount");
	setUp.setOption(len, "--len", "Template Length");
	setUp.setOption(pcrRounds, "--pcrRounds", "Number of PCR rounds");
	setUp.setOption(initialPcrRounds, "--initialPcrRounds", "Initial number of PCR rounds");
	setUp.processSeq(false);
	if(initialPcrRounds > pcrRounds){
		std::stringstream ss;
		ss << "Error in " << __PRETTY_FUNCTION__ << std::endl;
		ss << "initialPcrRounds should be less than pcrRounds" << std::endl;
		ss << "pcrRounds:" << pcrRounds << std::endl;
		ss << "initialPcrRounds:" << initialPcrRounds << std::endl;
		throw std::runtime_error{ss.str()};
	}
	setUp.setOption(errorRate, "--errorRate", "Polymerase Error Rate");
	setUp.setOption(numThreads, "--numThreads", "Number of Threads to Use");
	setUp.setOption(finalReadAmount, "--finalReadAmount", "Final Read Amount to sample");
	setUp.processDirectoryOutputName("simPCR_TODAY", true);
	setUp.processVerbose();
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	njh::randomGenerator gen;
	uint64_t intErrorRate = errorRate * std::numeric_limits<uint64_t>::max();
	std::string seq = "";
	if("" == setUp.pars_.seq_){
		seq = simulation::evenRandStr(len, std::vector<char>{'A', 'C', 'G', 'T'}, gen);
	}else{
		seq = setUp.pars_.seq_;
	}
	std::unordered_map<std::string, uint64_t> seqCounts;
	std::mutex seqMapLock;
	auto finalPerfectAmount = static_cast<uint64_t>(startingTemplate
			* std::pow(2, initialPcrRounds));
	auto finalAmount = sim::runPcr(intErrorRate, numThreads, initialPcrRounds, seq,
			startingTemplate, "Seq", seqCounts, seqMapLock, true);
	std::cout << "Initial Pcr amplified Round " << initialPcrRounds  << std::endl;
	std::cout << "Final Amount If No Errors Round " << initialPcrRounds << ": " << finalPerfectAmount << std::endl;
	std::cout << "Final Amount Observed Round " << initialPcrRounds << ": "
			<< getPercentageString(finalAmount, finalPerfectAmount) << std::endl;
	std::ofstream logFile;
	openTextFile(logFile, OutOptions(bfs::path(setUp.pars_.directoryName_ + "info.tab.txt")));
	logFile << "Initial Pcr amplified Round " << initialPcrRounds  << std::endl;
	logFile << "Final Amount If No Errors Round " << initialPcrRounds << ": " << finalPerfectAmount << std::endl;
	logFile << "Final Amount Observed Round " << initialPcrRounds << ": "
			<< getPercentageString(finalAmount, finalPerfectAmount) << std::endl;
	std::ofstream sequenceRefOutFile;
	openTextFile(sequenceRefOutFile, OutOptions(bfs::path(setUp.pars_.directoryName_ + "ref.fasta")));
	sequenceRefOutFile << ">ref" << std::endl;
	sequenceRefOutFile << seq << std::endl;
	std::ofstream sequenceOutFile;
	openTextFile(sequenceOutFile, OutOptions(bfs::path(setUp.pars_.directoryName_ + "output.fasta")));
	std::cout << std::endl;
	std::cout << "Sampling..." << std::endl;
	std::mutex seqFileLock;
	std::unordered_map<std::string, std::string> seqs;
	seqs["Seq"] = seq;
	std::unordered_map<std::string, std::unordered_map < std::string, uint64_t>>  multipleSeqCounts;
	multipleSeqCounts["Seq"] = seqCounts;
	std::unordered_map<std::string, std::pair<uint64_t, uint64_t>> sampleNumber;
	if (initialPcrRounds == pcrRounds) {
		sampleNumber["Seq"] = sim::sampleReadsWithoutReplacement("Seq", seq, seqCounts,
				finalReadAmount, sequenceOutFile, seqFileLock, setUp.pars_.verbose_);

	} else {
		sampleNumber = sim::sampleReadsWithoutReplacementFinishPCR(seqs,
				multipleSeqCounts, finalReadAmount, sequenceOutFile, seqFileLock,
				pcrRounds - initialPcrRounds, intErrorRate, numThreads, setUp.pars_.verbose_);
	}
	std::cout << std::endl;
	std::cout << "Sampled" << std::endl;
	std::cout << finalReadAmount << std::endl;
	std::cout << getPercentageString(sampleNumber["Seq"].first, finalReadAmount) << std::endl;
	logFile << "Final Sampled Amount With No Errors: "
			<< getPercentageString(sampleNumber["Seq"].first, finalReadAmount) << std::endl;
	logFile << "Final Sampled Amount With Errors: "
			<< getPercentageString(sampleNumber["Seq"].second, finalReadAmount) << std::endl;
	return 0;
}





int readSimulatorRunner::createMinorVariantsFromKnown(const njh::progutils::CmdArgs & inputCommands) {
  seqSetUp setUp(inputCommands);
  //bool allowRepeatMutations = false;
  bool gradientDisparatePairs = false;
  uint32_t pairNumber = 7;
  uint32_t mutationOffSet = 100;
  std::string mutations = "1";   // format if just number,
  														   // do a random mutation for that number, or 1024:A which means to mutation base 1024 to A,
  														   // separate these by ; to do multiple mutations, eg. 1000:A, 1409:C
  std::string abundances = "10"; // will determine the abundance of the minor alleles, can be separated by : to set for different alleles but must equal number of minor alleles
  															 // eg. 10:10,20:30,
  setUp.setOption(mutationOffSet, "--mutationOffSet", "The number of bases from the end to simulate");
  setUp.setOption(mutations, "--mutations", "The mutations to simulate");
  setUp.setOption(gradientDisparatePairs, "--gradientDisparatePairs", "Create pairs of variants that are distantly related");
  setUp.setOption(pairNumber, "--pairNumber", "Number of pairs to create if --gradientDisparatePairs is used");
  setUp.setOption(abundances, "--abundances", "The abundances of the minor alleles");
  uint64_t seed = 0;
  bool useSeed = setUp.setOption(seed, "--seed", "Use this seed to generate the same exact library");
  //setUp.setOption(allowRepeatMutations, "--allowRepeatMutations", "Allow different variants to have same mutations");
  setUp.processSeq(true); //read in one seq record
  setUp.processDirectoryOutputName("simReads_createMinorVariantsFromKnown_TODAY", true);
  setUp.processWritingOptions();
  setUp.processVerbose();
  setUp.processDebug();
  setUp.finishSetUp(std::cout);
  setUp.startARunLog(setUp.pars_.directoryName_);
  MinorVariantsMixture mixture(setUp.pars_.seqObj_.seqBase_);
  if(useSeed){
  	mixture.gen_.seedNum(seed);
  }
  std::ofstream seedFile;
  openTextFile(seedFile, setUp.pars_.directoryName_ + "seed","",false, true);
  seedFile << mixture.gen_.currentSeed_;
  if(gradientDisparatePairs){
  	mixture.createGradientDisparatePairsMixture(pairNumber);
  }else{
  	 mixture.addVariants(mutations, mutationOffSet);
  }
	//process abundance
  mixture.createAbundancesForMixtures(abundances);
	//write out variants
  mixture.writeFasta(setUp.pars_.directoryName_);
  //write out variant mutation info
  mixture.outputVariantMutationInfofile(setUp.pars_.directoryName_ + "variantsInfo.tab.txt");
  //write out abundance info
  mixture.outputAbundanceFile(setUp.pars_.directoryName_ + "abundance.tab.txt");

  return 0;
}

int readSimulatorRunner::simPcrShotgunSequences(const njh::progutils::CmdArgs & inputCommands) {
  readSimulatorSetUp setUp(inputCommands);
	uint32_t mean = 240;
	uint32_t std = 81;
	uint32_t minLen = 8;
	uint64_t startingTemplate = 3000;
	uint64_t finalReadAmount = 5000;
	uint32_t pcrRounds = 20;
	uint32_t initialPcrRounds = 10;
	long double errorRate = 3.5e-06;
	uint32_t numThreads = 2;
	std::string abundanceFile = "";
	bool simIllumina = false;
	bool sim454 = false;
	bool both = false;
  double mismatchDensity = 0.25;
  uint32_t minOverLap = 10;
  size_t pairedEndLength = 200;
  setUp.setOption(pairedEndLength, "--pairedEndLength", "Paired End Length For Illumina Simulation");
	setUp.setOption(minLen, "--minLen", "Minimum Length to allow in the simulation");
	setUp.setOption(std, "--std", "Standard Deviation of the mean read length");
  setUp.setOption(mean, "--readLength,-l", "The mean read length to simulate, the actual lengths will be distributed around this");
	setUp.setOption(both, "--both", "Simulate both Illumina and 454");
	if(both){
		sim454 = true;
		simIllumina = true;
	}
  setUp.setOption(simIllumina, "--simIllumina", "simulate illumina reads");
	setUp.setOption(sim454, "--sim454", "simulate 454 reads");

	//setUp.setOption(referenceFile, "--referenceFile", "Reference sequences fasta File", true);
	setUp.processDefaultReader(VecStr{"-fasta", "-fastq"},true);
	setUp.setOption(abundanceFile, "--abundanceFile", "Abundance File", true);
	setUp.setOption(startingTemplate, "--startingTemplate", "Staring Template Amount");
	setUp.setOption(pcrRounds, "--pcrRounds", "Number of PCR rounds");
	setUp.setOption(initialPcrRounds, "--initialPcrRounds", "Number of Initial rounds of PCR before sampling");
	setUp.setOption(errorRate, "--errorRate", "Polymerase Error Rate");
	setUp.setOption(numThreads, "--numThreads", "Number of Threads to Use");
	setUp.setOption(finalReadAmount, "--finalReadAmount", "Final Read Amount to sample");
	setUp.processDirectoryOutputName("simPCR_TODAY", true);
	setUp.processVerbose();
	//setUp.processWritingOptions();
	setUp.finishSetUp(std::cout);
	uint32_t maxOverlap = pairedEndLength;
	setUp.startARunLog(setUp.pars_.directoryName_);
	setUp.writeParametersFile(setUp.pars_.directoryName_ + "parameters.tab.txt", false, true);
	njh::randomGenerator gen;
	uint64_t intErrorRate = errorRate * std::numeric_limits<uint64_t>::max();
  SeqInput reader(setUp.pars_.ioOptions_);
  reader.openIn();
  auto inReads = reader.readAllReads<readObject>();
	std::vector<readObject> refReads = inReads;
	auto refNamesFromSeqs = readVec::getNames(refReads);
	auto libraryAbundances = sim::processLibraryAbundances(abundanceFile, refNamesFromSeqs, std::numeric_limits<uint32_t>::max());
	VecStr libNames;
	auto allDir = njh::files::makeDir(setUp.pars_.directoryName_,njh::files::MkdirPar( "all"));
	for (const auto & lib : libraryAbundances) {
		std::vector<std::shared_ptr<seqInfo>> libReads;
		for (const auto & names : lib.second) {
			libReads.emplace_back(
					std::make_shared<seqInfo>(
							readVec::getReadByName(refReads, names.first).seqBase_));
			libReads.back()->frac_ = names.second;
		}
		auto libName = "lib"
				+ leftPadNumStr<size_t>(lib.first, libraryAbundances.size());
		sim::simShotgunLibFast(libReads, setUp.pars_.directoryName_, libName, intErrorRate,
				startingTemplate, finalReadAmount, pcrRounds, initialPcrRounds,
				numThreads, setUp.pars_.verbose_, mean, std, minLen);
		libNames.emplace_back(libName);
		std::ofstream runLogs;
		auto libDir = njh::appendAsNeededRet(njh::files::join(setUp.pars_.directoryName_, libName).string(), "/");
		if(sim454 || simIllumina){
			openTextFile(runLogs, OutOptions(bfs::path(libDir + "simProgramLogs.json")));
		}
		if (simIllumina) {
			uint32_t illuminaAttempts = 10; //art fails for no reason sometimes
			//pstreams freezes if art_illumina reports warnings which include length warnings
			SeqIOOptions opts = SeqIOOptions::genFastaIn(njh::files::join(libDir, "reads.fasta"));
			opts.out_.outFilename_ = libDir + "longReads.fasta";
			SeqIO reader(opts);
			reader.openIn();
			reader.openOut();
			seqInfo seq;

			while (reader.readNextRead(seq)) {
				if (len(seq) >= pairedEndLength) {
					reader.write(seq);
				}
			}
			std::string simCmdIllumina = "art_illumina -amp -p -na -i "
					+ libDir + "longReads.fasta -l " + estd::to_string(pairedEndLength) + " -f 1 -o "
					+ libDir + "amplicon_pair_dat";
			auto simOutPutIllumina = njh::sys::run(VecStr{simCmdIllumina});
			if(setUp.pars_.debug_){
				std::cout << "That finished" << std::endl;
			}
			runLogs << simOutPutIllumina.toJson();
			uint32_t numberOfAttempts = 1;
			while(!simOutPutIllumina.success_ && numberOfAttempts <= illuminaAttempts){
				simOutPutIllumina = njh::sys::run(VecStr{simCmdIllumina});
				runLogs << simOutPutIllumina.toJson();
				++numberOfAttempts;
			}
			std::string flashCmdTemplate = "flash "
							+ libDir + "amplicon_pair_dat1.fq "
							+ libDir + "amplicon_pair_dat2.fq "
							+ "-o " + libDir + "amplicon_pair_dat "
							" --min-overlap " + estd::to_string(minOverLap) +
							" --max-overlap " + estd::to_string(maxOverlap) +
							" -x " + estd::to_string(mismatchDensity) +
							" -t " + estd::to_string(numThreads) +
							" 2>&1 | tee "
							+ libDir + "flash.log";
			auto simOutPutFlash = njh::sys::run(VecStr { flashCmdTemplate });
			runLogs << simOutPutFlash.toJson();
			njh::files::bfs::copy(
					njh::files::make_path(libDir,
							"amplicon_pair_dat.extendedFrags.fastq"),
					njh::files::make_path(allDir, libName + "_stitchedIllumina.fastq"));
		}
		if(sim454) {
			std::string simCmd454 = "454sim -d $(echo $(dirname $(which 454sim))/gen) "
					+ libDir + "reads.fasta -o " + libDir
					+ "reads.sff";
			auto simOutPut454 = njh::sys::run(VecStr{simCmd454});
			runLogs << simOutPut454.toJson();
			auto sffOpts = SeqIOOptions::genFastqOut(libDir + "454_reads.fastq");
			sffOpts.inFormat_ = SeqIOOptions::inFormats::SFFBIN;
			sffOpts.firstName_ = libDir + "reads.sff";
		  SeqIO reader(sffOpts);
		  reader.openIn();
		  reader.openOut();
		  seqInfo sffRead;
		  while(reader.readNextRead(sffRead)){
		  	reader.write(sffRead);
		  }
			njh::files::bfs::copy(
								njh::files::make_path(libDir,
										"454_reads.fastq"),
								njh::files::make_path(allDir, libName + "_454.fastq"));
		}
	}
	std::cout << "Done" << std::endl;
	return 0;
}









}  // namespace njh
