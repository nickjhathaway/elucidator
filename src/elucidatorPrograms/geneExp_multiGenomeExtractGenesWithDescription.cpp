//
// Created by Nicholas Hathaway on 11/7/24.
//
#include "geneExp.hpp"
#include <njhseq/objects/Gene/GeneFromGffs.hpp>
#include "elucidator/objects/BioDataObject.h"
#include "elucidator/seqToolsUtils/seqToolsUtils.hpp"

#include <TwoBit.h>
#include <njhseq/GenomeUtils/GenomeMapping/MultiGenomeMapper.hpp>
#include <njhseq/objects/helperObjects/AminoAcidPositionInfo.hpp>


namespace njhseq {




int geneExpRunner::multiGenomeExtractGenesWithDescription(const njh::progutils::CmdArgs & inputCommands){
	MultiGenomeMapper::inputParameters genomePars;
	bool doNotSkipEmpty = false;
  std::set<std::string> descriptions;
	std::string geneName;

	seqSetUp setUp(inputCommands);
	setUp.setOption(descriptions, "--descriptions", "descriptions", true);
	setUp.setOption(geneName, "--geneName", "geneName", true);
	setUp.setOption(genomePars.genomeDir_, "--genomeDir", "a directory containing genomes", true);
	setUp.setOption(genomePars.gffDir_, "--gffDir", "a directory containing annotation files", true);
	// setUp.setOption(genomePars.gffIntersectPars_.filterSubRegionFeatures_, "--excludeFeatures", "gene features to exclude");
	setUp.setOption(genomePars.gffIntersectPars_.selectFeatures_, "--features", "gene features to extract");
	setUp.setOption(genomePars.acceptableGenomeExtensions_, "--acceptableGenomeExtensions", "acceptable Genome Extensions");
	setUp.setOption(genomePars.acceptableGffExtensions_, "--acceptableGffExtensions", "acceptable Gff Extensions");
	setUp.setOption(genomePars.selectedGenomes_, "--selectGenomes", "selectGenomes");
	setUp.setOption(doNotSkipEmpty, "--doNotSkipEmpty", "do Not Skip Empty");


	setUp.processVerbose();
	setUp.processDirectoryOutputName(geneName, true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);

	//check genome and gff annotation directories

	auto genomeFnps = njh::files::gatherFiles(genomePars.genomeDir_, genomePars.acceptableGenomeExtensions_);
	auto gffFnps = njh::files::gatherFiles(genomePars.gffDir_, genomePars.acceptableGffExtensions_);
	std::set<std::string> genomeNames;
	for (const auto & g : genomeFnps) {
		genomeNames.emplace(bfs::basename(g));
	}

	std::set<std::string> gffNames;
	for (const auto & g : gffFnps) {
		gffNames.emplace(bfs::basename(g));
	}
	std::vector<std::string> uniqueGenomeNames;
	std::vector<std::string> uniqueGffNames;
	std::vector<std::string> inBoth;
	njh::decompose_sets(genomeNames.begin(), genomeNames.end(),
		gffNames.begin(), gffNames.end(),
		std::back_inserter(uniqueGenomeNames),
		std::back_inserter(uniqueGffNames),
		std::back_inserter(inBoth));

	if (!uniqueGenomeNames.empty() || !uniqueGffNames.empty()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << ", error " << "\n"
				<< "uniqueGenomeNames: " << njh::conToStr(uniqueGenomeNames, ",") << "\n"
				<< "uniqueGffNames: " << njh::conToStr(uniqueGffNames, ",") << "\n";
		throw std::runtime_error{ss.str()};
	}

	MultiGenomeMapper genomes(genomePars);
	genomes.loadInGenomes();
	genomes.loadGffFnps();
	std::set<std::string> allGeneIDs;

	for (const auto & g : genomes.genomes_) {
		std::set<std::string> geneIDs;
		auto regions = MultiGenomeMapper::gatherGffRegionsWithDescriptions(g.second->gffFnp_, descriptions, genomePars);
		for (const auto & region : regions) {
			geneIDs.emplace(region.meta_.getMeta("ID"));
			// regionsOut << region.genBedRecordCore().toDelimStrWithExtra() << std::endl;
		}
		if (doNotSkipEmpty && geneIDs.empty()) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << ", error " << "no genes extracted for " << g.first << "\n";
			throw std::runtime_error{ss.str()};
		} else if (!geneIDs.empty()) {
			auto outputDir = njh::files::make_path(setUp.pars_.directoryName_, njh::pasteAsStr(g.first, "_", geneName, "GeneInfos"));
			njh::files::makeDir(outputDir);
			GeneFromGffs::gffRecordIDsToGeneInfoPars pars;
			pars.inputFile = g.second->gffFnp_;
			pars.twoBitFnp = g.second->fnpTwoBit_;
			pars.outOpts.outFilename_ = njh::files::make_path(outputDir, g.first);
			pars.ids = geneIDs;

			GeneFromGffs::gffRecordIDsToGeneInfo(pars);
			bfs::copy_file(njh::files::make_path(outputDir, g.first + "_allTranscripts.bed"),njh::files::make_path(setUp.pars_.directoryName_, njh::pasteAsStr(g.first, "_", geneName, "Genes.bed")));

			concatenateFiles(njh::files::gatherFiles(outputDir, "_cDNA.fasta", true), njh::files::make_path(outputDir, "allCDNA.fasta"));
			concatenateFiles(njh::files::gatherFiles(outputDir, "_gDNA.fasta", true), njh::files::make_path(outputDir, "allGDNA.fasta"));
			concatenateFiles(njh::files::gatherFiles(outputDir, "_protein.fasta", true), njh::files::make_path(outputDir, "allProtein.fasta"));
			concatenateFiles(njh::files::gatherFiles(outputDir, "_withUTR.bed", true), njh::files::make_path(outputDir, "allWithUTR.bed"));
			concatenateFiles(njh::files::gatherFiles(outputDir, "_exonIntronPositions.bed", true), njh::files::make_path(outputDir, "allExonIntronPositions.bed"));
			allGeneIDs.insert(geneIDs.begin(), geneIDs.end());
		}

	}

	concatenateFiles(njh::files::gatherFiles(setUp.pars_.directoryName_, "_cDNA.fasta", true), njh::files::make_path(setUp.pars_.directoryName_, "allCDNA.fasta"));
	concatenateFiles(njh::files::gatherFiles(setUp.pars_.directoryName_, "_gDNA.fasta", true), njh::files::make_path(setUp.pars_.directoryName_, "allGDNA.fasta"));
	concatenateFiles(njh::files::gatherFiles(setUp.pars_.directoryName_, "_protein.fasta", true), njh::files::make_path(setUp.pars_.directoryName_, "allProtein.fasta"));
	concatenateFiles(njh::files::gatherFiles(setUp.pars_.directoryName_, "_withUTR.bed", true), njh::files::make_path(setUp.pars_.directoryName_, "allWithUTR.bed"));
	concatenateFiles(njh::files::gatherFiles(setUp.pars_.directoryName_, "_exonIntronPositions.bed", true), njh::files::make_path(setUp.pars_.directoryName_, "allExonIntronPositions.bed"));
	concatenateFiles(njh::files::gatherFiles(setUp.pars_.directoryName_, "_allTranscripts.bed", true), njh::files::make_path(setUp.pars_.directoryName_, "all.bed"));
	return 0;
}




} /* namespace njhseq */


