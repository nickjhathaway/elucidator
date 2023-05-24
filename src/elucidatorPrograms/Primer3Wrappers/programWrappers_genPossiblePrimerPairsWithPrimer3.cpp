//
// Created by Nicholas Hathaway on 5/21/23.
//

#include "elucidatorPrograms/programWrappers.hpp"
#include "elucidator/objects/BioDataObject.h"
#include "elucidator/BamToolsUtils.h"

#include <SeekDeep/objects/TarAmpSetupUtils/PrimersAndMids.hpp>

#include <njhseq/concurrency/PairwisePairFactory.hpp>
#include <njhseq/objects/seqContainers/CollapsedHaps.hpp>
#include <njhseq/GenomeUtils/GenomeExtraction/ParsingAlignmentInfo/GenomeExtractResult.hpp>

namespace njhseq {

std::string primer3RoughQmdFile = R"(---
title: "Primer3 Primer Pair Locations"
author: "Nicholas Hathaway"
format:
  html:
    theme: cosmo
    fig-height: 10
    fig-width: 15
    toc: true
    toc-depth: 4 # default is 3
    toc-title: Contents
    number-sections: true
    anchor-sections: true
    smooth-scroll: true
    highlight-style: textmate
    page-layout: full
    code-fold: true
    code-tools: true
    code-link: true
---


```{r setup, echo=FALSE, message=FALSE}
require(knitr)
require(DT)
require(tidyverse)
require(stringr)
require(dbscan)
require(ape)
require(rwantshue)
require(tsne)
require(ggforce)
require(GGally)
require(plotly)
require(heatmaply)
require(ComplexHeatmap)
require(fastmatch)
require(HaplotypeRainbows)
require(ggtree)
require(ggthemes)
#turn off messages and warnings and make it so output isn't prefixed by anything,
#default is to put "##" in front of all output for some reason
#also set tidy to true so code is wrapped properly
opts_chunk$set(message=FALSE, warning=FALSE, comment = "", cache = F)
options(width = 200)
`%!in%` <- Negate(`%in%`)
scheme <- iwanthue(seed = 42, force_init = TRUE)

myFormula= x~y
library(ggpmisc)
`%!in%` <- Negate(`%in%`)

sofonias_theme = theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank() )+
  theme(axis.line.x = element_line(color="black", linewidth = 0.3),axis.line.y =
          element_line(color="black", linewidth = 0.3))+
  theme(text=element_text(size=12, family="Helvetica"))+
  theme(axis.text.y = element_text(size=12))+
  theme(axis.text.x = element_text(size=12)) +
   theme(legend.position = "bottom") +
   theme(plot.title = element_text(hjust = 0.5))

sofonias_theme_xRotate = sofonias_theme +
  theme(axis.text.x = element_text(size=12, angle = -90, vjust = 0.5, hjust = 0))
create_dt <- function(x){
  DT::datatable(x,
                extensions = 'Buttons',
                options = list(dom = 'Blfrtip',
                               buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                               lengthMenu = list(c(10,25,50,-1),
                                                 c(10,25,50,"All"))),
                filter = "top")
}
colorPalette_08 = c("#2271B2","#F748A5","#359B73","#F0E442","#D55E00","#3DB7E9","#E69F00","#000000")
colorPalette_12 = c("#E20134","#FF6E3A","#008DF9","#8400CD","#FFC33B","#9F0162","#009F81","#FF5AAF","#00FCCF","#00C2F9","#FFB2FD","#A40122")
colorPalette_15 = c("#F60239","#003C86","#EF0096","#9400E6","#009FFA","#008169","#68023F","#00DCB5","#FFCFE2","#FF71FD","#7CFFFA","#6A0213","#008607","#00E307","#FFDC3D")

createColorListFromDf <- function(df, colorPalette = colorPalette_12, iwanthudSeed = rnorm(1) * 100){
  colorList = list()
  for (dfColname in colnames(df)) {
    levels = sort(unique(df[[dfColname]]))
    scheme <- iwanthue(seed = iwanthudSeed, force_init = TRUE)
    if (length(levels) <= length(colorPalette_12)) {
      levelsCols = colorPalette_12[1:length(levels)]
    } else{
      levelsCols = scheme$hex(length(levels))
    }
    names(levelsCols) = levels
    colorList[[dfColname]] = levelsCols
  }
  return(colorList)
}
```

<style type="text/css">
div.content {
  max-width: 1800px;
  margin-left: auto;
  margin-right: auto;
}
</style>

```{r}
diversity = readr::read_tsv("diversityPerPair/diversityPerPrimerPair.tab.txt")

primer3_results_ampLocs = readr::read_tsv("primer3_results_ampLocs.bed", col_names = F) %>%
  rename(SeqIDPrimerPairName = X4) %>%
  left_join(diversity)%>%
  arrange(desc(he)) %>%
  mutate(id = row_number())
refRegion = readr::read_tsv("../refSeq.bed", col_names = F)
variableRegions_divs = readr::read_tsv("../subRegionInfo/divMeasuresPerVarRegion.tab.txt", col_names = T)

variableRegions = readr::read_tsv("../subRegionInfo/0_ref_variable_genomic.bed", col_names = F)%>%
  left_join(variableRegions_divs %>%
              rename(X4 = id))
conservedRegions = readr::read_tsv("../subRegionInfo/0_ref_sharedLocs_genomic.bed", col_names = F)

primer3Results = readr::read_tsv("primer3_results.tab.txt")
create_dt(primer3Results)
```

Visualization of the diversity of the primer pairs created, bars are sorted by diversity, on the bottom is the input region and the variable and conserved regions around which the primers are designed.

```{r}
#| column: screen-inset-shaded
#| fig-column: screen-inset-shaded
ggplotly(ggplot() +
  geom_rect(aes(
    xmin = X2,
    xmax = X3,
    ymin = 0,
    ymax = 1,
    fill = "conservedRegions"
  ),
  data = conservedRegions) +
  geom_rect(aes(
    xmin = X2,
    xmax = X3,
    ymin = 1,
    ymax = 2,
    he = he,
    ExpP5 = ExpP5,
    fill = "variableRegions"
  ),
  data = variableRegions) +
  geom_rect(aes(
    xmin = X2,
    xmax = X3,
    ymin = -1,
    ymax = 0,
    fill = "refRegion"
  ),
  data = refRegion)  +
  # scale_fill_tableau() +
  # ggnewscale::new_scale_fill() +
  geom_rect(aes(
    xmin = X2,
    xmax = X3,
    ymin = 3 + id - 0.9,
    ymax = 3 + id - 0.2,
    fill = factor(ExpP5),
    SeqIDPrimerPairName = SeqIDPrimerPairName,
    he = he, len = len, start = start, end = end
  ),
  data = primer3_results_ampLocs %>% mutate(len = X5, start = X2, end = X3)) +
  scale_fill_manual(values = createColorListFromDf(tibble(factorExpP5 = factor(c("conservedRegions", "variableRegions", "refRegion",
                                                                                 primer3_results_ampLocs$ExpP5))))$factorExpP5) +
  sofonias_theme )
```
)";


int programWrapperRunner::genPossiblePrimerPairsWithPrimer3(const njh::progutils::CmdArgs &inputCommands) {


	Primer3Runner::Primer3Options p3Opts;

	uint32_t errorAllowed = 0;
	double majorAFVCFFreq = 0.99;
	bfs::path vcfFnp;
	uint32_t blastExpandSize = 10;
	uint32_t numThreads = 1;
	bfs::path excludeRegionsFnp;
	bfs::path regionsOfInterestFnp;
	bfs::path bedFnp;
	bfs::path twoBitFnp;

	bfs::path createSharedSubSegmentsFromRefSeqsDirFnp;

	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();

	setUp.setOption(majorAFVCFFreq, "--majorAFVCFFreq",
									"If supplying a VCF file than change SNPs with this frequency or more before designing primers");
	setUp.setOption(vcfFnp, "--vcfFnp", "VCF Fnp");
	setUp.setOption(createSharedSubSegmentsFromRefSeqsDirFnp, "--createSharedSubSegmentsFromRefSeqsDirFnp", "createSharedSubSegmentsFromRefSeqs directory to compute diversity results of possible primers");

	setUp.setOption(excludeRegionsFnp, "--excludeRegions", "exclude Regions");
	setUp.setOption(regionsOfInterestFnp, "--regionsOfInterest", "regions Of Interest");
	setUp.setOption(numThreads, "--numThreads", "numThreads");

	p3Opts.setInsertSizeOptions(setUp);
	p3Opts.setPrimaryOptions(setUp);
	p3Opts.setPrimerSizeOpts(setUp);
	p3Opts.setReturnOptions(setUp);

	//setUp.setOption(task, "--task", "primer picking task, examples include:generic(default), pick_sequencing_primers, pick_primer_list");


	uint32_t extendRegion = std::max<uint32_t>(p3Opts.PRIMER_MAX_SIZE + 5, std::round(p3Opts.maxSize * 1.5));

	setUp.setOption(extendRegion, "--extendRegion", "extend region by this amount");

	setUp.setOption(bedFnp, "--bedFnp", "genomic locations", true);
	setUp.setOption(twoBitFnp, "--twoBit", "two Bit file", true);

	setUp.processDirectoryOutputName("genPrimers_TODAY", true);

	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);

	std::vector<std::shared_ptr<Bed6RecordCore>> excludeRegions;
	std::vector<std::shared_ptr<Bed6RecordCore>> regionsOfInterest;

	if (bfs::exists(regionsOfInterestFnp)) {
		regionsOfInterest = getBeds(regionsOfInterestFnp);
	}

	if (bfs::exists(excludeRegionsFnp)) {
		excludeRegions = getBeds(excludeRegionsFnp);
	}
	std::vector<VCFVariant> allVariants;
	if (bfs::exists(vcfFnp)) {
		InputStream vcfIn(vcfFnp);
		std::string line;
		while (njh::files::crossPlatGetline(vcfIn, line)) {
			if (njh::beginsWith(line, "#")) {
				continue;
			}
			auto variants = VCFVariant::readVCFLine(line);
			njh::addConToVec(allVariants, variants);
		}
	}


	auto regions = getBeds(bedFnp);

	std::unordered_map<std::string, std::shared_ptr<Bed6RecordCore>> nameToRegion;
	std::unordered_map<std::string, std::shared_ptr<Bed6RecordCore>> templateToRegion;

	TwoBit::TwoBitFile tReader(twoBitFnp);
	auto chromLens = tReader.getSeqLens();
	std::set<std::string> regionNames;
	for (const auto &reg: regions) {
		regionNames.emplace(reg->name_);
	}

	BedUtility::coordSort(excludeRegions, false);
	BedUtility::coordSort(regions, false);
	BedUtility::coordSort(regionsOfInterest, false);

	njh::sort(allVariants, [](const VCFVariant &reg1, const VCFVariant &reg2) {
		if (reg1.region_.chrom_ == reg2.region_.chrom_) {
			if (reg1.region_.start_ == reg2.region_.start_) {
				return reg1.region_.end_ < reg2.region_.end_;
			} else {
				return reg1.region_.start_ < reg2.region_.start_;
			}
		} else {
			return reg1.region_.chrom_ < reg2.region_.chrom_;
		}
	});

	OutputStream inputRegionsExpanded(njh::files::make_path(setUp.pars_.directoryName_, "inputRegionsExpanded.bed"));


	//SEQUENCE_TARGET, SEQUENCE_EXCLUDED_REGION
	bfs::path regionsOfInterestDir = njh::files::make_path(setUp.pars_.directoryName_, "regionsOfInterest");
	if (!regionsOfInterest.empty()) {

	}
	njh::files::makeDir(njh::files::MkdirPar{regionsOfInterestDir});
	{
		OutputStream primer3Input(njh::files::make_path(setUp.pars_.directoryName_, "primer3File.txt"));
		primer3Input << "PRIMER_THERMODYNAMIC_PARAMETERS_PATH=" << p3Opts.primer3ConfigPath.string() << std::endl;

		for (const auto &originalReg: regions) {
			std::shared_ptr<Bed6RecordCore> reg = std::make_shared<Bed6RecordCore>(*originalReg);
			reg->strand_ = '+';
			std::vector<Bed6RecordCore> regionsOfInterest_withinRegion;
			std::vector<Bed6RecordCore> regionsOfInterest_relToRegion;

			std::vector<Bed6RecordCore> excludeRegions_withinRegion;
			std::vector<Bed6RecordCore> excludeRegions_relToRegion;

			std::vector<VCFVariant> variants_withinRegion;
			std::vector<VCFVariant> variants_relToRegion;
			{
				for (const auto &inter: regionsOfInterest) {
					if (inter->chrom_ < reg->chrom_) {
						//bother regions are sorted if we haven't reached this region's chromosome yet
						continue;
					}
					if (inter->chrom_ > reg->chrom_) {
						//both regions are sorted so if we run into this we can break
						break;
					}
					if (inter->chrom_ == reg->chrom_ && inter->chromStart_ >= reg->chromEnd_) {
						//both regions are sorted so if we run into this we can break
						break;
					}
					if (inter->chrom_ == reg->chrom_ && inter->chromEnd_ < reg->chromStart_) {
						continue;
					}
					if (reg->overlaps(*inter, 1) && inter->chromStart_ >= reg->chromStart_ &&
							inter->chromEnd_ <= reg->chromEnd_) {
						//take only completely within this region
						regionsOfInterest_withinRegion.emplace_back(*inter);
					}
				}
			}
			OutputStream regionsOfInterest_withinRegion_bedOut(
							njh::files::make_path(regionsOfInterestDir, reg->name_ + ".bed"));
			for (const auto &bed: regionsOfInterest_withinRegion) {
				regionsOfInterest_withinRegion_bedOut << bed.toDelimStrWithExtra() << std::endl;
			}
			bool lenIsScore = reg->score_ == reg->length();
			BedUtility::extendLeftRight(*reg, extendRegion, extendRegion, njh::mapAt(chromLens, reg->chrom_));
			if (lenIsScore) {
				reg->score_ = reg->length();
			}
			inputRegionsExpanded << reg->toDelimStrWithExtra() << std::endl;
			nameToRegion[reg->name_] = reg;
			{
				for (const auto &inter: excludeRegions) {
					if (inter->chrom_ < reg->chrom_) {
						//bother regions are sorted if we haven't reached this region's chromosome yet
						continue;
					}
					if (inter->chrom_ > reg->chrom_) {
						//both regions are sorted so if we run into this we can break
						break;
					}
					if (inter->chrom_ == reg->chrom_ && inter->chromStart_ >= reg->chromEnd_) {
						//both regions are sorted so if we run into this we can break
						break;
					}
					if (inter->chrom_ == reg->chrom_ && inter->chromEnd_ < reg->chromStart_) {
						continue;
					}
					if (reg->overlaps(*inter, 1)) {
						//check first if we already have this exact region for exclusion, this can happen if input generated from variant calls and there is more than biallelic SNPs etc
						bool alreadyHave = false;
						for (const auto &already: excludeRegions_withinRegion) {
							if (already.genUIDFromCoords() == inter->genUIDFromCoords()) {
								alreadyHave = true;
								break;
							}
						}

						if (!alreadyHave) {
							excludeRegions_withinRegion.emplace_back(*inter);
						}
					}
				}
				//book end the exclusion regions
				for (auto &exclude: excludeRegions_withinRegion) {
					exclude.chromStart_ = std::max(exclude.chromStart_, reg->chromStart_);
					exclude.chromEnd_ = std::min(exclude.chromEnd_, reg->chromEnd_);
				}
			}
			//variant regions;
			{
				for (const auto &var: allVariants) {
					if (var.region_.chrom_ < reg->chrom_) {
						//bother regions are sorted if we haven't reached this region's chromosome yet
						continue;
					}
					if (var.region_.chrom_ > reg->chrom_) {
						//both regions are sorted so if we run into this we can break
						break;
					}
					if (var.region_.chrom_ == reg->chrom_ && var.region_.start_ >= reg->chromEnd_) {
						//both regions are sorted so if we run into this we can break
						break;
					}
					if (var.region_.chrom_ == reg->chrom_ && var.region_.end_ < reg->chromStart_) {
						continue;
					}
					if (reg->overlaps(var.region_.genBed3RecordCore(), 1)) {
						variants_withinRegion.emplace_back(var);
					}
				}
				//book end the exclusion regions
				for (auto &var: variants_withinRegion) {
					var.region_.start_ = std::max<size_t>(var.region_.start_, reg->chromStart_);
					var.region_.end_ = std::min<size_t>(var.region_.end_, reg->chromEnd_);
				}
			}
			//convert to regions relative to the template region
			for (const auto &roi: regionsOfInterest_withinRegion) {
				auto modRegion = roi;
				modRegion.chromStart_ = modRegion.chromStart_ - reg->chromStart_;
				modRegion.chromEnd_ = modRegion.chromEnd_ - reg->chromStart_;
				regionsOfInterest_relToRegion.emplace_back(modRegion);
			}
			for (const auto &var: variants_withinRegion) {
				auto modRegion = var;
				modRegion.region_.start_ = modRegion.region_.start_ - reg->chromStart_;
				modRegion.region_.end_ = modRegion.region_.end_ - reg->chromStart_;
				variants_relToRegion.emplace_back(modRegion);
				if (var.freq_ < majorAFVCFFreq) {
					//add variants
					excludeRegions_relToRegion.emplace_back(modRegion.region_.genBedRecordCore());
				}
			}

			for (const auto &exclude: excludeRegions_withinRegion) {
				bool withinVCF = false;
				for (const auto &var: variants_withinRegion) {
					if (var.region_.createUidFromCoords() == exclude.genUIDFromCoords()) {
						withinVCF = true;
						break;
					}
				}
				if (!withinVCF) {
					auto modRegion = exclude;
					modRegion.chromStart_ = modRegion.chromStart_ - reg->chromStart_;
					modRegion.chromEnd_ = modRegion.chromEnd_ - reg->chromStart_;
					excludeRegions_relToRegion.emplace_back(modRegion);
				}
			}
			{ //collapse down regions that fall completely within other regions
				//internal filtration
				auto bedCoordSorterFunc =
								[](const Bed6RecordCore & reg1, const Bed6RecordCore & reg2) {
									if(reg1.chrom_ == reg2.chrom_) {
										if(reg1.chromStart_ == reg2.chromStart_) {
											if(reg1.chromEnd_ == reg2.chromEnd_){
												return reg1.name_ < reg2.name_;
											}else{
												return reg1.chromEnd_ > reg2.chromEnd_;
											}
										} else {
											return reg1.chromStart_ < reg2.chromStart_;
										}
									} else {
										return reg1.chrom_ < reg2.chrom_;
									}
								};
				njh::sort(excludeRegions_withinRegion, bedCoordSorterFunc);
				std::vector<Bed6RecordCore> excludeRegions_withinRegion_filter;
				for(const auto & excludeRegion : excludeRegions_withinRegion){
					if(excludeRegions_withinRegion_filter.empty()){
						excludeRegions_withinRegion_filter.emplace_back(excludeRegion);
					}else{
						auto completelyInBackRegion = excludeRegion.chromStart_ >= excludeRegions_withinRegion_filter.back().chromStart_ && excludeRegion.chromEnd_ <= excludeRegions_withinRegion_filter.back().chromEnd_;
						if(!completelyInBackRegion){
							excludeRegions_withinRegion_filter.emplace_back(excludeRegion);
						}
					}
				}
				excludeRegions_withinRegion = excludeRegions_withinRegion_filter;
				njh::sort(excludeRegions_relToRegion, bedCoordSorterFunc);
				std::vector<Bed6RecordCore> excludeRegions_relToRegion_filter;
				for (const auto &excludeRegion: excludeRegions_relToRegion) {
					if (excludeRegions_relToRegion_filter.empty()) {
						excludeRegions_relToRegion_filter.emplace_back(excludeRegion);
					} else {
						auto completelyInBackRegion =
										excludeRegion.chromStart_ >= excludeRegions_relToRegion_filter.back().chromStart_ &&
										excludeRegion.chromEnd_ <= excludeRegions_relToRegion_filter.back().chromEnd_;
						if (!completelyInBackRegion) {
							excludeRegions_relToRegion_filter.emplace_back(excludeRegion);
						}
					}
				}
				excludeRegions_relToRegion = excludeRegions_relToRegion_filter;
			}

			auto seqTemplate = GenomicRegion(*reg).extractSeq(tReader);
			//make any changes to the templates as needed
			for (const auto &var: variants_relToRegion) {
				//only doing SNPs right now cause otherwise handling INDELs and the changes to the genomic locations would be kind of a nightmare
				if (var.vtype_ == VCFVariant::VarType::SNP && var.freq_ >= majorAFVCFFreq) {
					seqTemplate.seq_[var.region_.start_] = var.variant_[0]; //change base
				}
			}
			std::string seqID = reg->name_;
			std::stringstream defaultPars;
			auto p3OptsJson = p3Opts.toJson();
			defaultPars << "SEQUENCE_TEMPLATE=" << seqTemplate.seq_ << std::endl;
			defaultPars << "SEQUENCE_ID=" << seqID << std::endl;
			defaultPars << "PRIMER_FIRST_BASE_INDEX=0" << std::endl;
			defaultPars << "PRIMER_TASK=" << p3Opts.task << std::endl;
			defaultPars << "PRIMER_SEQUENCING_LEAD=5" << std::endl;
			defaultPars << "PRIMER_EXPLAIN_FLAG=1" << std::endl;
			defaultPars << "PRIMER_PRODUCT_SIZE_RANGE=" << p3Opts.minSize << "-" << p3Opts.maxSize << std::endl;
			defaultPars << "PRIMER_NUM_RETURN=" << p3Opts.PRIMER_NUM_RETURN << std::endl;
			defaultPars << "PRIMER_MAX_SIZE=" << p3Opts.PRIMER_MAX_SIZE << std::endl;
			defaultPars << "PRIMER_MIN_SIZE=" << p3Opts.PRIMER_MIN_SIZE << std::endl;
			defaultPars << "PRIMER_OPT_SIZE=" << p3Opts.PRIMER_OPT_SIZE << std::endl;
			defaultPars << "PRIMER_OPT_GC_PERCENT=" << p3Opts.PRIMER_OPT_GC_PERCENT << std::endl;
			defaultPars << "PRIMER_INTERNAL_MIN_GC=" << p3Opts.PRIMER_INTERNAL_MIN_GC << std::endl;
			for(const auto & add : p3Opts.additionalOpts){
				if(!njh::in(add.first, p3OptsJson.getMemberNames())){
					defaultPars << add.first << "=" << add.second << std::endl;
				}
			}
			if (!excludeRegions_relToRegion.empty()) {
				std::string SEQUENCE_EXCLUDED_REGION;
				for (const auto &exclude: excludeRegions_relToRegion) {
					if (!SEQUENCE_EXCLUDED_REGION.empty()) {
						SEQUENCE_EXCLUDED_REGION += " ";
					}
					SEQUENCE_EXCLUDED_REGION += njh::pasteAsStr(exclude.chromStart_, ",",
																											exclude.chromEnd_ - exclude.chromStart_);
				}
				defaultPars << "SEQUENCE_EXCLUDED_REGION=" << SEQUENCE_EXCLUDED_REGION << std::endl;
			}
			if (!regionsOfInterest_relToRegion.empty()) {
				for (const auto &roi: regionsOfInterest_relToRegion) {
					std::string SEQUENCE_TARGET = njh::pasteAsStr(roi.chromStart_, ",", roi.chromEnd_ - roi.chromStart_);
					primer3Input << defaultPars.str();
					primer3Input << "SEQUENCE_TARGET=" << SEQUENCE_TARGET << std::endl;
					primer3Input << "=" << std::endl;
				}
			} else {
				primer3Input << defaultPars.str();
				primer3Input << "=" << std::endl;
			}
		}
	}

	std::unordered_map<std::string, std::string> primerNameToFor;
	std::unordered_map<std::string, std::string> primerNameToRev;


	std::string primer3Cmd = njh::pasteAsStr("cd ", setUp.pars_.directoryName_, " && ", "primer3_core ",
																					 "primer3File.txt > primer3_output.txt 2> primer3_error.txt");
	auto runOutput = njh::sys::run({primer3Cmd});
	OutputStream primer3RunLog(njh::files::make_path(setUp.pars_.directoryName_, "primer3RunLog.json"));
	primer3RunLog << runOutput.toJson() << std::endl;

	auto primer3ResultsFnp = njh::files::make_path(setUp.pars_.directoryName_, "primer3_output.txt");
	auto results = Primer3Runner::Primer3ResultsGeneric::parsePrimer3OutputResults(primer3ResultsFnp, true);

	OutputStream primer3ResultsOut(njh::files::make_path(setUp.pars_.directoryName_, "primer3_results.tab.txt"));
	OutputStream primer3ResultsPrimerTableOut(njh::files::make_path(setUp.pars_.directoryName_, "primer3_results_primerTable.tab.txt"));

	OutputStream primer3ResultsPrimerLocs(
					njh::files::make_path(setUp.pars_.directoryName_, "primer3_results_primerLocs.bed"));
	OutputStream primer3ResultsInsertLocs(
					njh::files::make_path(setUp.pars_.directoryName_, "primer3_results_insertLocs.bed"));
	OutputStream primer3ResultsAmpLocs(njh::files::make_path(setUp.pars_.directoryName_, "primer3_results_ampLocs.bed"));
	auto bestDir = njh::files::make_path(setUp.pars_.directoryName_, "bestPrimerPerTarget");
	njh::files::makeDir(njh::files::MkdirPar{bestDir});
	OutputStream best_primer3ResultsPrimerLocs(njh::files::make_path(bestDir, "primer3_results_primerLocs.bed"));
	OutputStream best_primer3ResultsInsertLocs(njh::files::make_path(bestDir, "primer3_results_insertLocs.bed"));
	OutputStream best_primer3ResultsAmpLocs(njh::files::make_path(bestDir, "primer3_results_ampLocs.bed"));
	OutputStream best_primersTable(njh::files::make_path(bestDir, "primer3_results_primers.tab.txt"));
	best_primersTable << "target\tforward\treverse" << std::endl;
	primer3ResultsPrimerTableOut << "target\tforward\treverse" << std::endl;
	primer3ResultsOut << njh::conToStr(
					VecStr{"seqID", "SeqIDPrimerPairName", "chrom", "ampStart", "ampEnd", "primerPairName", "compl_any_th", "compl_end_th",
								 "pair_penalty", "pair_penalty_noSize", "product_size", "product_gc_percent",
								 "left_name", "left_seq", "left_start", "left_end", "left_size", "left_end_stability",
								 "left_gc_percent", "left_hairpin_th", "left_penalty", "left_penalty_noSize", "left_self_any_th",
								 "left_self_end_th", "left_tm", "left_tm_hairpin_diff", "left_problems",
								 "right_name", "right_seq", "right_start", "right_end", "right_size", "right_end_stability",
								 "right_gc_percent", "right_hairpin_th", "right_penalty", "right_penalty_noSize", "right_self_any_th",
								 "right_self_end_th", "right_tm", "right_tm_hairpin_diff", "right_problems"
					}, "\t") << std::endl;
	VecStr noTargetsFor;
	std::stringstream noTargetsBed;

	uint32_t maxChromLen = 0;
	for (const auto &chrom: chromLens) {
		if (chrom.second > maxChromLen) {
			maxChromLen = chrom.second;
		}
	}

	std::unordered_set<std::string> allPrimers;


	std::vector<std::shared_ptr<njhseq::Primer3Runner::PrimerPairGeneric>> allBestPrimers;
	std::unordered_set<std::string> alreadyHaveUID;
	for (const auto &res: results) {
		const auto &region = nameToRegion[res->sequence_id_];
		bool targeted = !res->sequence_target_.empty();
		if (0 == res->primer_pair_num_returned_) {
			//should log which attempts had 0 returned
			std::string targetedName;
			if (targeted) {
				targetedName = njh::pasteAsStr(region->chrom_, "-",
																			 region->chromStart_ + res->sequence_target_.front().start_, "-",
																			 region->chromStart_ + res->sequence_target_.front().start_ +
																			 res->sequence_target_.front().size_);
				noTargetsFor.emplace_back(targetedName);
				noTargetsBed << region->chrom_
										 << "\t" << region->chromStart_ + res->sequence_target_.front().start_
										 << "\t"
										 << region->chromStart_ + res->sequence_target_.front().start_ + res->sequence_target_.front().size_
										 << "\t" << targetedName
										 << "\t" << res->sequence_target_.front().size_
										 << "\t" << "+" << std::endl;
			}
			continue;
		}

		if (targeted) {
			auto bestPrimer = res->getLowestPenaltyPair();
			bool passBestPrimer = true;
			for (const auto &bestPrimerAlready: allBestPrimers) {
				if (bestPrimer->left_.seq_ == bestPrimerAlready->left_.seq_ &&
						bestPrimer->right_.seq_ == bestPrimerAlready->right_.seq_) {
					passBestPrimer = false;
					break;
				}
			}
			if (passBestPrimer) {
				allBestPrimers.emplace_back(std::make_shared<njhseq::Primer3Runner::PrimerPairGeneric>(*bestPrimer));
			} else {
				continue;
			}

			std::string targetedName;
			if (targeted) {
				targetedName = njh::pasteAsStr(region->chrom_, "-",
																			 region->chromStart_ + res->sequence_target_.front().start_, "-",
																			 region->chromStart_ + res->sequence_target_.front().start_
																			 + res->sequence_target_.front().size_);
			}
			auto primerPairName = njh::pasteAsStr(res->sequence_id_ + (targeted ? "--" + targetedName : ""), "--PP",
																						njh::replaceString(bestPrimer->name_, "PRIMER_PAIR_", ""));

			GenomicRegion insertRegion(
							res->sequence_id_,
							region->chrom_,
							region->chromStart_ + bestPrimer->left_.forwardOrientationPos_.start_ +
							bestPrimer->left_.forwardOrientationPos_.size_,
							region->chromStart_ + bestPrimer->right_.forwardOrientationPos_.start_,
							false);
//			insertRegion.uid_ = njh::pasteAsStr(insertRegion.chrom_, "-",
//																					njh::leftPadNumStr<size_t>(insertRegion.start_, maxChromLen), "-",
//																					njh::leftPadNumStr<size_t>(insertRegion.end_, maxChromLen));
			insertRegion.uid_ = primerPairName;
			uint32_t uidCount = 0;
			while (njh::in(insertRegion.uid_, alreadyHaveUID)) {
				++uidCount;
				insertRegion.uid_ = njh::pasteAsStr(insertRegion.chrom_, "-",
																						njh::leftPadNumStr<size_t>(insertRegion.start_, maxChromLen), "-",
																						njh::leftPadNumStr<size_t>(insertRegion.end_, maxChromLen), "--", uidCount);
			}
			best_primer3ResultsInsertLocs << insertRegion.genBedRecordCore().toDelimStrWithExtra() << std::endl;
			GenomicRegion ampRegion(insertRegion.uid_,
															region->chrom_,
															region->chromStart_ + bestPrimer->getStart(),
															region->chromStart_ + bestPrimer->getEnd(),
															false);

			best_primer3ResultsAmpLocs << ampRegion.genBedRecordCore().toDelimStrWithExtra() << std::endl;

			GenomicRegion leftPrimerRegion(
							insertRegion.uid_ + "-forwardPrimer",
							region->chrom_,
							region->chromStart_ + bestPrimer->left_.forwardOrientationPos_.start_,
							region->chromStart_ + bestPrimer->left_.forwardOrientationPos_.start_ +
							bestPrimer->left_.forwardOrientationPos_.size_,
							false);
			GenomicRegion rightPrimerRegion(
							insertRegion.uid_ + "-reversePrimer",
							region->chrom_,
							region->chromStart_ + bestPrimer->right_.forwardOrientationPos_.start_,
							region->chromStart_ + bestPrimer->right_.forwardOrientationPos_.start_ +
							bestPrimer->right_.forwardOrientationPos_.size_,
							true);
			best_primer3ResultsPrimerLocs << leftPrimerRegion.genBedRecordCore().toDelimStrWithExtra() << std::endl;
			best_primer3ResultsPrimerLocs << rightPrimerRegion.genBedRecordCore().toDelimStrWithExtra() << std::endl;

			best_primersTable << insertRegion.uid_
												<< "\t" << bestPrimer->left_.seq_
												<< "\t" << bestPrimer->right_.seq_ << std::endl;
		} else {
			auto bestPrimers = res->getLowestPenaltyNonOverlappingPairs();
			for (const auto &bestPrimer: bestPrimers) {
				bool passBestPrimer = true;
				for (const auto &bestPrimerAlready: allBestPrimers) {
					if (bestPrimer->left_.seq_ == bestPrimerAlready->left_.seq_ &&
							bestPrimer->right_.seq_ == bestPrimerAlready->right_.seq_) {
						passBestPrimer = false;
						break;
					}
				}
				if (passBestPrimer) {
					allBestPrimers.emplace_back(std::make_shared<njhseq::Primer3Runner::PrimerPairGeneric>(*bestPrimer));
				} else {
					continue;
				}
				GenomicRegion insertRegion(
								res->sequence_id_,
								region->chrom_,
								region->chromStart_ + bestPrimer->left_.forwardOrientationPos_.start_ +
								bestPrimer->left_.forwardOrientationPos_.size_,
								region->chromStart_ + bestPrimer->right_.forwardOrientationPos_.start_,
								false);
				insertRegion.uid_ = njh::pasteAsStr(insertRegion.chrom_, "-",
																						njh::leftPadNumStr<size_t>(insertRegion.start_, maxChromLen), "-",
																						njh::leftPadNumStr<size_t>(insertRegion.end_, maxChromLen));
				uint32_t uidCount = 0;
				while (njh::in(insertRegion.uid_, alreadyHaveUID)) {
					++uidCount;
					insertRegion.uid_ = njh::pasteAsStr(insertRegion.chrom_, "-",
																							njh::leftPadNumStr<size_t>(insertRegion.start_, maxChromLen), "-",
																							njh::leftPadNumStr<size_t>(insertRegion.end_, maxChromLen), "--",
																							uidCount);
				}
				best_primer3ResultsInsertLocs << insertRegion.genBedRecordCore().toDelimStrWithExtra() << std::endl;
				GenomicRegion ampRegion(insertRegion.uid_,
																region->chrom_,
																region->chromStart_ + bestPrimer->getStart(),
																region->chromStart_ + bestPrimer->getEnd(),
																false);
				best_primer3ResultsAmpLocs << ampRegion.genBedRecordCore().toDelimStrWithExtra() << std::endl;

				GenomicRegion leftPrimerRegion(
								insertRegion.uid_ + "-forwardPrimer",
								region->chrom_,
								region->chromStart_ + bestPrimer->left_.forwardOrientationPos_.start_,
								region->chromStart_ + bestPrimer->left_.forwardOrientationPos_.start_ +
								bestPrimer->left_.forwardOrientationPos_.size_,
								false);
				GenomicRegion rightPrimerRegion(
								insertRegion.uid_ + "-reversePrimer",
								region->chrom_,
								region->chromStart_ + bestPrimer->right_.forwardOrientationPos_.start_,
								region->chromStart_ + bestPrimer->right_.forwardOrientationPos_.start_ +
								bestPrimer->right_.forwardOrientationPos_.size_,
								true);
				best_primer3ResultsPrimerLocs << leftPrimerRegion.genBedRecordCore().toDelimStrWithExtra() << std::endl;
				best_primer3ResultsPrimerLocs << rightPrimerRegion.genBedRecordCore().toDelimStrWithExtra() << std::endl;

				best_primersTable << insertRegion.uid_
													<< "\t" << bestPrimer->left_.seq_
													<< "\t" << bestPrimer->right_.seq_ << std::endl;
			}
		}

		for (const auto &primerPair: res->primerPairs_) {
			uint32_t targetStart = primerPair->left_.forwardOrientationPos_.start_;
			uint32_t targetEnd =
							primerPair->right_.forwardOrientationPos_.start_ + primerPair->right_.forwardOrientationPos_.size_;
			std::string target = res->sequence_template_.substr(targetStart, targetEnd - targetStart);
			uint32_t chromStart = region->chromStart_ + targetStart;
			uint32_t chromEnd = region->chromStart_ + targetEnd;

			DNABaseCounter counter;
			counter.increase(target);
			counter.setFractions();

			//since above we split out targets to be just a single at a time will just take the front
			if (res->sequence_target_.size() > 1) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error "
					 << ", target size shouldn't be more than 1, res->sequence_target_.size(): " << res->sequence_target_.size()
					 << "\n";
				throw std::runtime_error{ss.str()};
			}
			std::string targetedName;
			if (targeted) {
				targetedName = njh::pasteAsStr(region->chrom_, "-",
																			 region->chromStart_ + res->sequence_target_.front().start_, "-",
																			 region->chromStart_ + res->sequence_target_.front().start_
																			 + res->sequence_target_.front().size_);
			}
			auto primerPairName = njh::pasteAsStr(res->sequence_id_ + (targeted ? "--" + targetedName : ""), "--PP",
																						njh::replaceString(primerPair->name_, "PRIMER_PAIR_", ""));
			{
				MetaDataInName ampMeta;
				ampMeta.addMeta("SeqID", res->sequence_id_);
				ampMeta.addMeta("PrimerID", primerPair->name_);
				if (targeted) {
					ampMeta.addMeta("target", targetedName);
				}
				primer3ResultsAmpLocs << region->chrom_
															<< "\t" << region->chromStart_ + targetStart
															<< "\t" << region->chromStart_ + targetEnd
															<< "\t" << primerPairName
															<< "\t" << targetEnd - targetStart
															<< "\t" << "+"
															<< "\t" << ampMeta.createMetaName() << std::endl;
			}

			{
				MetaDataInName insertMeta;
				insertMeta.addMeta("SeqID", res->sequence_id_);
				insertMeta.addMeta("PrimerID", primerPair->name_);
				if (targeted) {
					insertMeta.addMeta("target", targetedName);
				}
				primer3ResultsInsertLocs << region->chrom_
																 << "\t" << region->chromStart_ + primerPair->left_.forwardOrientationPos_.start_ +
																						primerPair->left_.forwardOrientationPos_.size_
																 << "\t" << region->chromStart_ + primerPair->right_.forwardOrientationPos_.start_
																 << "\t" << primerPairName
																 << "\t" << primerPair->right_.forwardOrientationPos_.start_ -
																						(primerPair->left_.forwardOrientationPos_.start_ +
																						 primerPair->left_.forwardOrientationPos_.size_)
																 << "\t" << "+"
																 << "\t" << insertMeta.createMetaName()
																 << std::endl;
			}

			{
				MetaDataInName leftMeta;
				leftMeta.addMeta("SeqID", res->sequence_id_);
				leftMeta.addMeta("PrimerID", primerPair->name_);
				leftMeta.addMeta("PrimerName", primerPair->left_.name_);
				if (targeted) {
					leftMeta.addMeta("target", targetedName);
				}
				primer3ResultsPrimerLocs << region->chrom_
																 << "\t" << region->chromStart_ + primerPair->left_.forwardOrientationPos_.start_
																 << "\t" << region->chromStart_ + primerPair->left_.forwardOrientationPos_.start_ +
																						primerPair->left_.forwardOrientationPos_.size_
																 << "\t" << primerPairName
																 << "\t" << primerPair->left_.forwardOrientationPos_.size_
																 << "\t" << "+"
																 << "\t" << leftMeta.createMetaName() << std::endl;

				MetaDataInName rightMeta;
				rightMeta.addMeta("SeqID", res->sequence_id_);
				rightMeta.addMeta("PrimerID", primerPair->name_);
				rightMeta.addMeta("PrimerName", primerPair->right_.name_);
				if (targeted) {
					rightMeta.addMeta("target", targetedName);
				}
				primer3ResultsPrimerLocs << region->chrom_
																 << "\t" << region->chromStart_ + primerPair->right_.forwardOrientationPos_.start_
																 << "\t" << region->chromStart_ + primerPair->right_.forwardOrientationPos_.start_ +
																						primerPair->right_.forwardOrientationPos_.size_
																 << "\t" << primerPairName
																 << "\t" << primerPair->right_.forwardOrientationPos_.size_
																 << "\t" << "-"
																 << "\t" << rightMeta.createMetaName() << std::endl;
			}

			double leftPenaltyWithOutSize = primerPair->left_.penalty_ -
																			uAbsdiff(primerPair->left_.forwardOrientationPos_.size_, p3Opts.PRIMER_OPT_SIZE);
			double rightPenaltyWithOutSize = primerPair->right_.penalty_ -
																			 uAbsdiff(primerPair->right_.forwardOrientationPos_.size_,
																								p3Opts.PRIMER_OPT_SIZE);
			double pairPenaltyWithOutSize = leftPenaltyWithOutSize + rightPenaltyWithOutSize;

			primer3ResultsPrimerTableOut
							<< primerPairName
							<< "\t" << primerPair->left_.seq_
							<< "\t" << primerPair->right_.seq_ << std::endl;
			allPrimers.emplace(primerPair->left_.seq_);
			allPrimers.emplace(primerPair->right_.seq_);
			primerNameToFor[primerPairName] = primerPair->left_.seq_;
			primerNameToRev[primerPairName] = primerPair->right_.seq_;

			primer3ResultsOut << njh::conToStr(toVecStr(
							res->sequence_id_ + (targeted ? "--" + targetedName : ""),
							primerPairName,
							region->chrom_,


							chromStart,
							chromEnd,
							primerPair->name_,
							primerPair->compl_any_th_,
							primerPair->compl_end_th_,
							primerPair->penalty_,
							pairPenaltyWithOutSize,
							primerPair->product_size_,
							roundDecPlaces(counter.calcGcContent() * 100, 2),
							primerPair->left_.name_,
							primerPair->left_.seq_,
							region->chromStart_ + primerPair->left_.forwardOrientationPos_.start_,
							region->chromStart_ + primerPair->left_.forwardOrientationPos_.start_ +
							primerPair->left_.forwardOrientationPos_.size_,
							primerPair->left_.forwardOrientationPos_.size_,
							primerPair->left_.end_stability_,
							primerPair->left_.gc_percent_,
							primerPair->left_.hairpin_th_,
							primerPair->left_.penalty_,
							leftPenaltyWithOutSize,
							primerPair->left_.self_any_th_,
							primerPair->left_.self_end_th_,
							primerPair->left_.tm_,
							primerPair->left_.tm_ - primerPair->left_.hairpin_th_,
							njh::conToStr(primerPair->left_.problems_, ";"),

							primerPair->right_.name_,
							primerPair->right_.seq_,
							region->chromStart_ + primerPair->right_.forwardOrientationPos_.start_,
							region->chromStart_ + primerPair->right_.forwardOrientationPos_.start_ +
							primerPair->right_.forwardOrientationPos_.size_,
							primerPair->right_.forwardOrientationPos_.size_,
							primerPair->right_.end_stability_,
							primerPair->right_.gc_percent_,
							primerPair->right_.hairpin_th_,
							primerPair->right_.penalty_,
							rightPenaltyWithOutSize,
							primerPair->right_.self_any_th_,
							primerPair->right_.self_end_th_,
							primerPair->right_.tm_,
							primerPair->right_.tm_ - primerPair->right_.hairpin_th_,
							njh::conToStr(primerPair->right_.problems_, ";")
			), "\t") << std::endl;
		}

	}
	std::vector<std::string> allPrimersVec{allPrimers.begin(), allPrimers.end()};
	{
		OutputStream allPrimersOut(njh::files::make_path(setUp.pars_.directoryName_, "allPrimers.fasta"));
		for(const auto & item : iter::enumerate(allPrimersVec)){
			allPrimersOut << ">" << item.index << std::endl;
			allPrimersOut << item.element << std::endl;
		}
	}

	if (!noTargetsFor.empty()) {
		OutputStream noResultsForRegions(
						njh::files::make_path(setUp.pars_.directoryName_, "noTargetsForRegionsOfInterest.bed"));
		noResultsForRegions << noTargetsBed.str();
	}

	OutputStream primer3ResultsSummaryOut(
					njh::files::make_path(setUp.pars_.directoryName_,
																"primer3_results_summary.tab.txt"));
	primer3ResultsSummaryOut << njh::conToStr(VecStr{"seqID", "chrom", "start",
																									 "end", "len", "targeted", "targetedID", "targetedChrom",
																									 "targetedStart",
																									 "targetedEnd", "targetedSize", "excludedRegions",
																									 "excludedRegionsIDs", "primerTask",
																									 "pairNumReturned", "leftNumReturned", "rightNumReturned",
																									 "pairExplained",
																									 "leftExplained", "rightExplained", "warnings"}, "\t") << std::endl;

	for (const auto &res: results) {
		const auto &region = nameToRegion[res->sequence_id_];
		bool targeted = !res->sequence_target_.empty();
		std::string targetedName;
		if (targeted) {
			targetedName = njh::pasteAsStr(region->chrom_, "-",
																		 region->chromStart_ + res->sequence_target_.front().start_, "-",
																		 region->chromStart_ + res->sequence_target_.front().start_
																		 + res->sequence_target_.front().size_);
		}
		std::string seqExcludedRegions;
		std::string seqExcludedRegionsIDs;
		for (const auto &exclude: res->sequence_excluded_region_) {
			std::string excludeRel = njh::pasteAsStr(exclude.start_, exclude.size_, ";");
			std::string excludeID = njh::pasteAsStr(region->chrom_, "-",
																							region->chromStart_ + exclude.start_, "-",
																							region->chromStart_ + exclude.start_ + exclude.size_, ";");
			seqExcludedRegions += excludeRel;
			seqExcludedRegionsIDs += excludeID;
		}

		primer3ResultsSummaryOut << njh::conToStr(toVecStr(
						res->sequence_id_,
						region->chrom_,
						region->chromStart_,
						region->chromEnd_,
						region->length(),
						njh::boolToStr(targeted),
						(targeted ? targetedName : std::string("NA")),
						(targeted ? region->chrom_ : std::string("NA")),
						(targeted ? estd::to_string(region->chromStart_ + res->sequence_target_.front().start_) : std::string(
										"NA")),
						(targeted ? estd::to_string(
										region->chromStart_ + res->sequence_target_.front().start_ + res->sequence_target_.front().size_)
											: std::string("NA")),
						(targeted ? estd::to_string(res->sequence_target_.front().size_) : std::string("NA")),
						seqExcludedRegions,
						seqExcludedRegionsIDs,
						res->primer_task_,
						res->primer_pair_num_returned_,
						res->primer_left_num_returned_,
						res->primer_right_num_returned_,
						res->primer_pair_explain_,
						res->primer_left_explain_,
						res->primer_right_explain_,
						njh::conToStr(res->warnings_, ";")
		), "\t") << std::endl;
	}

	if(bfs::exists(createSharedSubSegmentsFromRefSeqsDirFnp)){
		ReAlignedSeq::genRealignmentPars reAlnPars;
		reAlnPars.extendAmount = blastExpandSize;
		reAlnPars.adjustForSoftClipping = true;

		auto diversityPerPairDir = njh::files::make_path(setUp.pars_.directoryName_, "diversityPerPair");
		njh::files::makeDir(njh::files::MkdirPar{diversityPerPairDir});

		//make a file of forward and reverse primers
		//read in unqiue sequences
		bfs::path uniqueSeqFnp = njh::files::make_path(createSharedSubSegmentsFromRefSeqsDirFnp, "uniqueSeqs.fasta.gz");
		std::vector<seqInfo> seqs = SeqInput::getSeqVec<seqInfo>(SeqIOOptions::genFastaInGz(uniqueSeqFnp));
		std::unordered_map<std::string, std::unordered_set<std::string>> seqNames;
		{
			TableReader nameReader(TableIOOpts::genTabFileIn(njh::files::make_path(createSharedSubSegmentsFromRefSeqsDirFnp, "uniqueSeqs_names.tab.txt.gz")));
			VecStr row;
			while(nameReader.getNextRow(row)){
				//name	number	inputNames
				auto names = njh::tokenizeString(row[nameReader.header_.getColPos("inputNames")], ";;");
				auto clusterSize = njh::StrToNumConverter::stoToNum<uint32_t>(row[nameReader.header_.getColPos("number")]);
				if(names.size() != clusterSize){
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error " << "for " << row[nameReader.header_.getColPos("name")] << " the number of input doesn't match cluster size, likily naems contains commas and can't be delimited well" << "\n";
					throw std::runtime_error{ss.str()};
				}
				auto namesSet = njh::vecToUOSet(names);
				seqNames[row[nameReader.header_.getColPos("name")]] = namesSet;
			}
		}
		std::vector<std::unordered_set<std::string>> setNames;
		setNames.reserve(seqs.size());
		for (auto &seq: seqs) {
			setNames.emplace_back(njh::mapAt(seqNames, seq.name_));
			seq.cnt_ = static_cast<double>(njh::mapAt(seqNames, seq.name_).size());
//			std::cout << seq.name_ << std::endl;
//			std::cout << "\tnjh::mapAt(seqNames, seq.name_).size(): " << njh::mapAt(seqNames, seq.name_).size() << std::endl;
		}
		auto uniqCorrectedSeqs = CollapsedHaps(seqs, setNames);
//		for(const auto & nameEnum : iter::enumerate(uniqCorrectedSeqs.names_)){
//			std::cout << uniqCorrectedSeqs.seqs_[nameEnum.index]->name_ << std::endl;
//			std::cout << "\t " << nameEnum.element.size() << std::endl;
//		}
		auto outputUniqSeqsFasta = njh::files::make_path(diversityPerPairDir, "uniqueSeqs.fasta");
		SeqOutput::write(seqs, SeqIOOptions::genFastaOut(outputUniqSeqsFasta));
//		std::cout << "uniqCorrectedSeqs.size(): " << uniqCorrectedSeqs.size() << std::endl;
//		std::cout << "uniqCorrectedSeqs.getTotalHapCount(): " << uniqCorrectedSeqs.getTotalHapCount() << std::endl;
//		exit(1);
		//make a blast database of the unique seqs
		std::string makedbCmd = "makeblastdb -in " + outputUniqSeqsFasta.string() + " -dbtype nucl -out " + njh::files::make_path(diversityPerPairDir, "uniqueSeqs").string();
		auto makedbCmdRunOut = njh::sys::run(VecStr{makedbCmd});
		if(!makedbCmdRunOut.success_){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error running " << makedbCmd << "\n";
			ss  << makedbCmdRunOut.stdErr_ << "\n";
			throw std::runtime_error{ss.str()};
		}
		std::string runBlastCmd =
						"blastn -query " + njh::files::make_path(setUp.pars_.directoryName_, "allPrimers.fasta").string() +
						" -perc_identity 90 -qcov_hsp_perc 90 -task blastn-short -word_size 6 -db " +
						njh::files::make_path(diversityPerPairDir, "uniqueSeqs").string() + " -outfmt 6 -max_target_seqs 100000000 > " + njh::files::make_path(diversityPerPairDir, "primerHits.tsv").string();
		auto runBlastCmdRunOut = njh::sys::run(VecStr{runBlastCmd});
		if(!runBlastCmdRunOut.success_){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error running " << runBlastCmd << "\n";
			ss  << makedbCmdRunOut.stdErr_ << "\n";
			throw std::runtime_error{ss.str()};
		}
		auto blastHitsFnp = njh::files::make_path(diversityPerPairDir, "primerHits.tsv");
		std::unordered_map<std::string, std::vector<std::shared_ptr<BLASTHitTab>>> blastHitsPerPrimers;

		std::vector<std::shared_ptr<BLASTHitTab>> allBlastHits;
		{
			BioDataFileIO<BLASTHitTab> reader{IoOptions(InOptions(blastHitsFnp))};
			reader.openIn();
			BLASTHitTab hit;
			while(reader.readNextRecord(hit)){
				allBlastHits.emplace_back(std::make_shared<BLASTHitTab>(hit));
				blastHitsPerPrimers[allPrimersVec[njh::StrToNumConverter::stoToNum<uint32_t>(hit.queryName_)]].emplace_back(allBlastHits.back());
			}
		}

		//subset unique to in between the primers and get diversity
		struct GenExtracRes{
			uint32_t forwardHits_{0};
			uint32_t reverseHits_{0};
			uint32_t extractCounts_{0};
		};

		std::unordered_map<std::string, GenExtracRes> genomeExtractionsResults;
		std::unordered_map<std::string, uint32_t> nameToPosUniqSeqs;

		for(const auto & seqEnum : iter::enumerate(uniqCorrectedSeqs.seqs_)){
			genomeExtractionsResults[seqEnum.element->name_] = GenExtracRes{};
			nameToPosUniqSeqs[seqEnum.element->name_] = seqEnum.index;
		}


		comparison allowableErrors;
		allowableErrors.hqMismatches_ = errorAllowed;
		njh::concurrent::LockableQueue<std::string> targetsQueue(njh::getVecOfMapKeys(primerNameToFor));
		std::mutex divMut;

		table divTable(
						VecStr{"SeqIDPrimerPairName", "totalHaplotypes", "uniqueHaplotypes", "singlets", "doublets", "expShannonEntropy",
									 "ShannonEntropyE", "effectiveNumOfAlleles", "SimpsonIndex", "he", "ExpP3", "ExpP4", "ExpP5",
									 "lengthPolymorphism"});

		OutputStream visOut(njh::files::make_path(setUp.pars_.directoryName_, "vis.qmd"));
		visOut << primer3RoughQmdFile << std::endl;
//		table performanceTab(
//						VecStr{"genome", "forwardPrimerHits", "reversePrimerHits", "extractionCounts", "target"});
		std::function<void()> extractPathway =
						[&targetsQueue, &blastHitsPerPrimers, &allowableErrors, &primerNameToFor, &primerNameToRev,
						 &allPrimers, &uniqCorrectedSeqs, &reAlnPars, &genomeExtractionsResults, &nameToPosUniqSeqs,
										&divTable, &divMut
//										,&performanceTab
										]() {

							std::string target;
							while (targetsQueue.getVal(target)) {

								uint64_t maxLen = 0;
								for (const auto &primer: allPrimers) {
									if (len(primer) > maxLen) {
										maxLen = len(primer);
									}
								}

								aligner alignerObj(maxLen, gapScoringParameters(5, 1, 0, 0, 5, 1, 0, 0, 5, 1), substituteMatrix(2, -2));
								std::unordered_map<std::string, std::vector<GenomeExtractResult>> genomeExtracts;
								for (const auto &seq: uniqCorrectedSeqs.seqs_) {

									std::vector<std::shared_ptr<BLASTHitTab>> forBlastHits;
									std::vector<std::shared_ptr<ReAlignedSeq>> forResults;
									{
										for (const auto &hit: blastHitsPerPrimers.at(primerNameToFor.at(target))) {
											if (hit->subjectName_ == seq->name_) {
												forBlastHits.emplace_back(hit);
												auto realignedSeq = ReAlignedSeq::genRealignment(*hit, primerNameToFor.at(target), alignerObj,
																																				 *seq, reAlnPars);
												if (1 == realignedSeq.comp_.distances_.query_.coverage_ &&
														allowableErrors.passErrorProfile(realignedSeq.comp_)) {
													forResults.emplace_back(std::make_shared<ReAlignedSeq>(realignedSeq));
												}
											}
										}
									}

									std::vector<std::shared_ptr<BLASTHitTab>> revBlastHits;
									std::vector<std::shared_ptr<ReAlignedSeq>> revResults;
									{
										for (const auto &hit: blastHitsPerPrimers.at(primerNameToRev.at(target))) {
											if (hit->subjectName_ == seq->name_) {
												revBlastHits.emplace_back(hit);
												auto realignedSeq = ReAlignedSeq::genRealignment(*hit, primerNameToRev.at(target), alignerObj,
																																				 *seq, reAlnPars);
												if (1 == realignedSeq.comp_.distances_.query_.coverage_ &&
														allowableErrors.passErrorProfile(realignedSeq.comp_)) {
													revResults.emplace_back(std::make_shared<ReAlignedSeq>(realignedSeq));
												}
											}
										}
									}

									forResults = ReAlignedSeq::getUniqueLocationResults(forResults);
									revResults = ReAlignedSeq::getUniqueLocationResults(revResults);


									genomeExtractionsResults[seq->name_].forwardHits_ = forResults.size();
									genomeExtractionsResults[seq->name_].reverseHits_ = revResults.size();
									if (!forResults.empty() && !revResults.empty()) {
										auto uniForRes = ReAlignedSeq::getUniqueLocationResults(forResults);
										auto uniRevRes = ReAlignedSeq::getUniqueLocationResults(revResults);
										genomeExtracts[seq->name_] = getPossibleGenomeExtracts(uniForRes, uniRevRes,
																																					 std::numeric_limits<uint32_t>::max());
									}
								}

								std::vector<seqInfo> subSeqs;
								std::vector<std::unordered_set<std::string>> subSeqsNames;
								for (auto &genome: genomeExtracts) {
									if (genome.second.empty()) {
										continue;
									}
									genomeExtractionsResults[genome.first].extractCounts_ = genome.second.size();\
                  for (auto &extract: genome.second) {
										extract.setRegion();
										subSeqs.emplace_back(uniqCorrectedSeqs.seqs_[nameToPosUniqSeqs.at(genome.first)]->getSubRead(
														extract.gRegionInner_->start_, extract.gRegionInner_->getLen()));
										subSeqsNames.emplace_back(uniqCorrectedSeqs.names_[nameToPosUniqSeqs.at(genome.first)]);
									}
								}
								auto subSeqsCollapsed = CollapsedHaps::collapseReads(subSeqs, subSeqsNames);
								auto divMeasuresForSubRegion = subSeqsCollapsed.getGeneralMeasuresOfDiversity(
												CollapsedHaps::GenPopMeasuresPar{});

								{
									std::lock_guard<std::mutex> lock(divMut);
									divTable.addRow(
													divMeasuresForSubRegion.getOut(subSeqsCollapsed, target, CollapsedHaps::GenPopMeasuresPar{})
									);
//									auto genomeKeys = getVectorOfMapKeys(genomeExtractionsResults);
//									njh::sort(genomeKeys);
//									for (const auto &genomeKey: genomeKeys) {
//										performanceTab.addRow(genomeKey,
//																					genomeExtractionsResults[genomeKey].forwardHits_,
//																					genomeExtractionsResults[genomeKey].reverseHits_,
//																					genomeExtractionsResults[genomeKey].extractCounts_,
//																					target);
//									}
								}


							}
						};

		njh::concurrent::runVoidFunctionThreaded(extractPathway, numThreads);
//		auto perTabOpts = TableIOOpts::genTabFileOut(njh::files::make_path(diversityPerPairDir, "extractionCounts"), true);
//		performanceTab.outPutContents(perTabOpts);

		auto divTabOpts = TableIOOpts::genTabFileOut(njh::files::make_path(diversityPerPairDir, "diversityPerPrimerPair"), true);
		divTable.naturlSortTable("SeqIDPrimerPairName", false);
		divTable.outPutContents(divTabOpts);
	}


	return 0;
}


} // namespace njhseq
