//
// Created by Nicholas Hathaway on 6/2/23.
//

#include <njhseq/IO/OutputStream.hpp>
#include "ampliconAnalysisRunner.hpp"

namespace njhseq {


int ampliconAnalysisRunner::processRawExtractByKmerPathWeaverResults(
				const njh::progutils::CmdArgs & inputCommands) {
	std::string PathWeaverFinalFastaFnp = "output_aboveCutOff.fasta";
	std::string PathWeaverResDirEnding = "_PathWeaver";
	VecStr extractionFilter{"initial", "final"};
	std::vector<bfs::path> inputDirectories;
	ampliconAnalysisSetUp setUp(inputCommands);
	uint32_t numThreads = 2;
	setUp.setOption(inputDirectories, "--inputDirectories", "Input Directories to analyze");
	setUp.setOption(extractionFilter, "--extractionFilter", "when combining extractions only take extractions with matching these values");
	setUp.setOption(numThreads, "--threads,-t", "Number of Threads");
	setUp.setOption(PathWeaverResDirEnding, "--PathWeaverResDirEnding", "PathWeaver Res Dir Ending");
	setUp.setOption(PathWeaverFinalFastaFnp, "--PathWeaverFinalFastaFnp", "PathWeaver Final Fasta Fnp");


	setUp.processDirectoryOutputName(true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	njh::files::checkExistenceThrow(inputDirectories, __PRETTY_FUNCTION__ );

	auto reportsDir = njh::files::make_path(setUp.pars_.directoryName_, "reports");
	njh::files::makeDir(njh::files::MkdirPar{reportsDir});

	{
		//concatenating extraction counts
		VecStr expectedHeaders{"iteration","sample","totalReads","target","count","frac","forwardCount","fracForward"};
		OutputStream bindExtractionsOut(njh::files::make_path(reportsDir, "allExtractions.tsv.gz"));
		bindExtractionsOut << "inputName\t" << njh::conToStr(expectedHeaders, "\t") << std::endl;
		for(const auto & fnp : inputDirectories){
			//table extractionCounts(TableIOOpts::genTabFileIn(njh::files::make_path(fnp, "extractionCounts.tab.txt")));
			table extractionCounts(njh::files::make_path(fnp, "extractionCounts.tab.txt"), "whitespace", true);
			extractionCounts.checkForColumnsThrow(expectedHeaders,__PRETTY_FUNCTION__ );
			for(const auto & row : extractionCounts){
				if(njh::in(row[extractionCounts.getColPos("iteration")], extractionFilter)){
					bindExtractionsOut<< bfs::basename(fnp);
					for(const auto & col : expectedHeaders){
						bindExtractionsOut << "\t" << row[extractionCounts.getColPos(col)];
					}
					bindExtractionsOut << std::endl;
				}
			}
		}
	}

	{
		for (const auto &inputDir: inputDirectories) {
			auto pathWeaverResDirs = njh::files::listAllFiles(njh::files::make_path(inputDir, "finalExtraction"), false,
																												{std::regex{njh::pasteAsStr(".*",PathWeaverResDirEnding, "$")}});
			std::cout << inputDir << std::endl;
			for(const auto & d : pathWeaverResDirs){
				std::vector<bfs::path> fnps;
				for (const auto& e : njh::files::dir(d.first)) {
					auto outputResults = njh::files::make_path(e.path(), PathWeaverFinalFastaFnp);
					if(bfs::is_directory(e) && bfs::exists(outputResults)){
						fnps.emplace_back(outputResults);
					}
				}
				std::cout << "\t" << njh::conToStr(fnps, ",") << std::endl;
			}
		}
	}



	return 0;
}



}  // namespace njhseq
