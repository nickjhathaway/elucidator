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
	std::string seqPat;

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
	setUp.setOption(seqPat, "--seqPat", "seqPat", true);
	setUp.setOption(inputDirectory, "--inputDirectory", "Input Directory to search");
	setUp.setOption(samples, "--samples", "Process input from only these samples");
	setUp.setOption(minReadCounts, "--minReadCounts", "min Read Counts to count a file");
	setUp.setOption(filesToInvestigation, "--filesToInvestigation", "files To Investigation");

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
	out << "sample\tregion\tfile\tcount" << std::endl;
	std::regex seqPatReg{seqPat};
	seqInfo seq;

	for(const auto & dir : directories){
		std::string sample = std::regex_replace(dir.string(), std::regex(pat), "");
		for(const auto & bed : beds){
			for(const auto & f : filesToInvestigation){
				bfs::path inputFnp = njh::files::make_path(dir, bed->name_, sample + "_extraction", f);
				if(setUp.pars_.debug_){
					std::cout << "dir: " << dir << ", " << "region: " << bed->name_ << ", " << "f: " << f << std::endl;
					std::cout << "\t" << inputFnp << std::endl;
				}
				if(bfs::exists(inputFnp)){
					auto inputFnpOpts = SeqIOOptions::genFastqIn(inputFnp);
					SeqInput reader(inputFnpOpts);
					reader.openIn();

					while(reader.readNextRead(seq)){
						std::ptrdiff_t const match_count(std::distance(
										std::sregex_iterator(seq.seq_.begin(), seq.seq_.end(), seqPatReg),
										std::sregex_iterator()));

						if(match_count >= minReadCounts){
							out << sample
							<< "\t" << bed->name_
							<< "\t" << f
							<< "\t" << match_count << std::endl;
						}
					}
				}
			}
		}
	}

	return 0;

}

}  //namespace njhseq

