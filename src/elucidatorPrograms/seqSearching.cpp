/*
 * seqSearching.cpp
 *
 *  Created on: Nov 12, 2018
 *      Author: nick
 */


#include "seqSearching.hpp"

#include "elucidator/objects/BioDataObject.h"
#include "elucidator/BamToolsUtils/BamUtilities.hpp"
#include "elucidator/BioRecordsUtils/BedUtility.hpp"

namespace njhseq {


seqSearchingRunner::seqSearchingRunner()
    : njh::progutils::ProgramRunner(
          {

					 addFunc("chopAndMap", chopAndMap, false),
					 addFunc("chopAndMapAndRefine", chopAndMapAndRefine, false),
           },
          "seqSearching") {}




struct ChopAndMapPars{
	uint32_t windowLength = 100;
	uint32_t windowStep = 10;
	uint32_t numThreads = 1;
	uint32_t perFragmentCount = 10;
	bfs::path genomeFnp = "";
	SeqIOOptions inOpts;
	bfs::path outputDirectory;
	bool debug = false;
};

void runChopAndMap(const ChopAndMapPars & pars){
	njh::sys::requireExternalProgramThrow("bwa");
	njh::files::checkExistenceThrow(pars.genomeFnp, __PRETTY_FUNCTION__);
	auto fragOutOpts = SeqIOOptions::genFastaOut(njh::files::make_path(pars.outputDirectory, "fragments.fasta"));

	{
		SeqOutput fragWriter(fragOutOpts);
		fragWriter.openOut();
		seqInfo seq;
		SeqInput reader(pars.inOpts);
		reader.openIn();
		while (reader.readNextRead(seq)) {
			if (len(seq) < pars.windowLength ) {
				std::cerr << "Seq: " << seq.name_ << " is too short" << std::endl;
			} else {
				for(auto const pos : iter::range<uint32_t>(0, len(seq) + 1 - pars.windowLength, pars.windowStep)){
					auto fragment = seq.getSubRead(pos, pars.windowLength);
					for(uint32_t fragCount = 0; fragCount <= pars.perFragmentCount; ++ fragCount){
						fragWriter.write(fragment);
					}
				}
			}
		}
	}

	bfs::path fragmentBamFnp =  njh::files::make_path(pars.outputDirectory, "fragments.sorted.bam");
	std::stringstream bwaMappCmd;
	bwaMappCmd << "bwa mem -M -t " << pars.numThreads << " "
			<< pars.genomeFnp  << " "
			<< fragOutOpts.out_.outName() << " "
			<< " 2> " << njh::files::make_path(pars.outputDirectory,"bwa.log.txt")
			<< " | samtools sort - -o " << fragmentBamFnp
			<< " && samtools index " << fragmentBamFnp ;
	auto bwaRunOutput = njh::sys::run({bwaMappCmd.str()});
	BioCmdsUtils::checkRunOutThrow(bwaRunOutput, __PRETTY_FUNCTION__);

}

int seqSearchingRunner::chopAndMap(const njh::progutils::CmdArgs & inputCommands) {
	ChopAndMapPars chopPars;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	chopPars.debug = setUp.pars_.debug_;
	setUp.setOption(chopPars.numThreads, "--numThreads", "Number of threads to use");
	setUp.setOption(chopPars.genomeFnp, "--genomeFnp", "Genome to map to", true);
	setUp.setOption(chopPars.perFragmentCount, "--perFragmentCount", "perFragmentCount");
	setUp.setOption(chopPars.windowLength, "--windowLength", "windowLength");
	setUp.setOption(chopPars.windowStep, "--windowStep", "windowStep");
	setUp.processReadInNames(VecStr{"--fasta", "--fastq"});
	setUp.processDirectoryOutputName(true);
	chopPars.inOpts = setUp.pars_.ioOptions_;
	chopPars.outputDirectory = setUp.pars_.directoryName_;
	setUp.finishSetUp(std::cout);

	runChopAndMap(chopPars);

	return 0;
}

int seqSearchingRunner::chopAndMapAndRefine(const njh::progutils::CmdArgs & inputCommands) {
	ChopAndMapPars chopPars;
	CoverageFinderPars covPars;
	RegionRefinementPars refinePars;
	uint32_t expandLeft = 0;
	uint32_t expandRight = 0;
	uint32_t minLength = 0;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	chopPars.debug = setUp.pars_.debug_;
	setUp.setOption(chopPars.numThreads, "--numThreads", "Number of threads to use");
	setUp.setOption(chopPars.genomeFnp, "--genomeFnp", "Genome to map to", true);
	setUp.setOption(chopPars.perFragmentCount, "--perFragmentCount", "perFragmentCount");
	setUp.setOption(chopPars.windowLength, "--windowLength", "windowLength");
	setUp.setOption(chopPars.windowStep, "--windowStep", "windowStep");
	setUp.setOption(refinePars.reOrient, "--reOrient", "reOrient");

	setUp.setOption(expandLeft, "--expandLeft", "expandLeft");
	setUp.setOption(expandRight, "--expandRight", "expandRight");

	setUp.setOption(minLength, "--minLength", "minLength");

	setUp.processReadInNames(VecStr{"--fasta", "--fastq"});
	setUp.processDirectoryOutputName(true);
	chopPars.inOpts = setUp.pars_.ioOptions_;
	chopPars.outputDirectory = setUp.pars_.directoryName_;
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	//chop
	runChopAndMap(chopPars);
	//determine coverage
	covPars.numThreads = chopPars.numThreads;
	covPars.coverageCutOff = 5;
	covPars.window = chopPars.windowLength;
	covPars.step = chopPars.windowStep;
	covPars.bams = njh::files::make_path(chopPars.outputDirectory, "fragments.sorted.bam").string();
	covPars.outOpts = OutOptions(njh::files::make_path(chopPars.outputDirectory, "coverage.bed"));
	RunCoverageFinder(covPars);
	//merge regions
	auto beds = getBed3s(covPars.outOpts.outName());
	njh::for_each(beds,
			[&expandLeft,&expandRight](const std::shared_ptr<Bed3RecordCore> & region) {
				if(0 != expandLeft) {
					BedUtility::extendLeft(*region, expandLeft);
				}
				if(0 != expandRight) {
					BedUtility::extendRight(*region, expandRight);
				}
			});
	njh::sort(beds,
			[](const std::shared_ptr<Bed3RecordCore> & region1, const std::shared_ptr<Bed3RecordCore> & region2) {
				return region1->chrom_ == region2->chrom_ ? (region1->chromStart_ == region2->chromStart_ ? region1->chromEnd_ < region2->chromEnd_ : region1->chromStart_ < region2->chromStart_): region1->chrom_ < region2->chrom_;
			});


	OutOptions mergedCoverageOpts(njh::files::make_path(chopPars.outputDirectory, "merged_coverage.bed"));
	{
		OutputStream mergedCoverageOut(mergedCoverageOpts);
		Bed3RecordCore currentRegion = *beds.front();
		for(const auto & regPos : iter::range<uint32_t>(1, beds.size())){
			if(currentRegion.chrom_ == beds[regPos]->chrom_ && currentRegion.overlaps(*beds[regPos], 1)){
				currentRegion.chromEnd_ = std::max(currentRegion.chromEnd_, beds[regPos]->chromEnd_);
			}else{
				Bed6RecordCore regOut(currentRegion.chrom_, currentRegion.chromStart_, currentRegion.chromEnd_, "", currentRegion.length(), '+');
				regOut.name_ = GenomicRegion(regOut).createUidFromCoords();

				mergedCoverageOut << regOut.toDelimStr() << std::endl;
				currentRegion = *beds[regPos];
			}
		}
		Bed6RecordCore regOut(currentRegion.chrom_, currentRegion.chromStart_, currentRegion.chromEnd_, "", currentRegion.length(), '+');
		regOut.name_ = GenomicRegion(regOut).createUidFromCoords();
		mergedCoverageOut << regOut.toDelimStr() << std::endl;
	}

	//refine regions
	refinePars.bedFnp = mergedCoverageOpts.outName();
	refinePars.bamFnp = njh::files::make_path(chopPars.outputDirectory, "fragments.sorted.bam");
	refinePars.outOpts = OutOptions(njh::files::make_path(chopPars.outputDirectory, "refined_merged.bed"));
	RunRegionRefinement(refinePars);

	if(0 != minLength){
		auto bedAgain = getBed3s(refinePars.outOpts.outName());
		OutOptions filteredOpts(njh::files::make_path(chopPars.outputDirectory, "filtered_refined_merged.bed"));
		OutputStream filteredOut(filteredOpts);
		for(const auto & b : bedAgain){
			if(b->length() >= minLength){
				filteredOut << b->toDelimStrWithExtra() << std::endl;
			}
		}
	}

	return 0;
}


}  // namespace njhseq

