//
// Created by Nicholas Hathaway on 7/31/23.
//
#include <njhseq/GenomeUtils.h>

#include "miscRunner.hpp"
#include "elucidator/simulation.h"
#include "elucidator/BamToolsUtils.h"
#include "elucidator/objects/dataContainers/graphs/ContigsCompareGraph.hpp"
#include "elucidator/BioRecordsUtils/BedUtility.hpp"


#include <njhseq/PopulationGenetics.h>
#include <njhseq/objects/seqContainers/CollapsedHaps.hpp>
#include <njhseq/objects/Gene/TranslatorByAlignment.hpp>
#include <boost/filesystem.hpp>


namespace njhseq {

int miscRunner::countPWExtractedReadsWithPattern(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path inputDirectory = "./";
	bfs::path bedFnp = "";
	std::set<std::string> samples;
	std::string pat;
	std::set<std::string> seqPats;
	uint32_t numThreads = 1;
	uint32_t minPatPerRead = 1;
	uint32_t minReadCounts = 2;
	OutOptions outOpts;

	VecStr filesToInvestigation = {"extracted.fastq",
																 "extracted_R1.fastq",
																 "extracted_R2.fastq",
																 "filteredPairs_extracted_R1.fastq",
																 "filteredPairs_extracted_R2.fastq",
																 "filteredSingles_extracted.fastq",
																 "thrownAwayMate_extracted.fastq"};

	seqSetUp setUp(inputCommands);

	setUp.processDebug();
	setUp.processVerbose();

	setUp.setOption(bedFnp, "--bedFnp", "regions", true);
	setUp.setOption(pat, "--pat", "file pattern", true);
	setUp.setOption(seqPats, "--seqPat", "sequence patterns to count", true);
	setUp.setOption(inputDirectory, "--inputDirectory", "Input Directory to search");
	setUp.setOption(samples, "--samples", "Process input from only these samples");
	setUp.setOption(minReadCounts, "--minReadCounts", "min Read Counts to count a file");
	setUp.setOption(minPatPerRead, "--minPatPerRead", "min count of patterns per Read Counts to count a read");
	setUp.setOption(filesToInvestigation, "--filesToInvestigation", "files To Investigation");
	setUp.setOption(numThreads, "--numThreads", "number of threads");

	setUp.processWritingOptions(outOpts);

	setUp.finishSetUp(std::cout);
//	setUp.startARunLog(setUp.pars_.directoryName_);


	std::vector<bfs::path> directories;
	if (samples.empty()) {
		auto allFiles = njh::files::listAllFiles(inputDirectory, false, {std::regex { ".*" + pat + "$" } });
		for (const auto &f : allFiles) {
			if (f.second) {
				directories.emplace_back(f.first);
			}
		}
	} else {
		for (const auto &samp : samples) {
			directories.emplace_back(njh::files::make_path(inputDirectory, samp + pat));
		}
	}

	if(directories.empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error no directories found in " << inputDirectory << " ending with " << pat << "\n";
		throw std::runtime_error{ss.str()};
	}

	auto beds = getBeds(bedFnp);

	OutputStream out(outOpts);
	std::mutex outMut;
	out << "sample\tregion\tfile\tseqPat\tcount" << std::endl;
	
	njh::concurrent::LockableVec<bfs::path> dirQueue(directories);
	std::function<void()> countSample = [&dirQueue,&beds,&filesToInvestigation,&out,&outMut, &minReadCounts,
																			 &seqPats, &setUp, &pat, &minPatPerRead](){
		bfs::path dir;
		while(dirQueue.getVal(dir)){
			std::string sample = std::regex_replace(dir.filename().string(), std::regex(pat), "");
			for(const auto & bed : beds){
				for(const auto & f : filesToInvestigation){
					bfs::path inputFnp = njh::files::make_path(dir, bed->name_, sample + "_extraction", f);
					if(setUp.pars_.debug_){
						std::cout << "dir: " << dir << ", " << "region: " << bed->name_ << ", " << "f: " << f << std::endl;
						std::cout << "\t" << inputFnp << std::endl;
					}
					if(bfs::exists(inputFnp)){
						for(const auto & seqPat : seqPats){
							std::regex seqPatReg{seqPat};
							seqInfo seq;
							auto inputFnpOpts = SeqIOOptions::genFastqIn(inputFnp);
							SeqInput reader(inputFnpOpts);
							reader.openIn();
							std::stringstream ss;
							uint32_t readCounts = 0;
							while(reader.readNextRead(seq)){
								std::ptrdiff_t const match_count(std::distance(
												std::sregex_iterator(seq.seq_.begin(), seq.seq_.end(), seqPatReg),
												std::sregex_iterator()));
								if(match_count >= minPatPerRead){
									++readCounts;
								}
							}
							if(readCounts >= minReadCounts){
								std::lock_guard<std::mutex> lock(outMut);
								out << sample
										<< "\t" << bed->name_
										<< "\t" << f
										<< "\t" << seqPat
										<< "\t" << readCounts << std::endl;
							}
						}
					}
				}
			}
		}
	};

	njh::concurrent::runVoidFunctionThreaded(countSample, numThreads);

	return 0;

}

}  //namespace njhseq

