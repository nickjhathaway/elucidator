/*
 * programWrappers_runnhmmscan.cpp
 *
 *  Created on: Jan 15, 2022
 *      Author: nick
 */






#include "programWrappers.hpp"
#include <njhseq/IO/SeqIO/SeqIO.hpp>
#include <njhseq/objects/BioDataObject/reading.hpp>
#include <njhseq/objects/BioDataObject/BedRecordCore.hpp>
#include <njhseq/objects/BioDataObject/BioRecordsUtils/HmmerUtility.hpp>
#include <njhseq/objects/BioDataObject/BioRecordsUtils/BedUtility.hpp>

namespace njhseq {

std::string plotDomainsQmd = R"(---
title: "Plotting Domain Hits"
format:
  html:
    theme: cosmo
    fig-height: 10
    fig-width: 15
    fig-align: center
    toc: true
    toc-depth: 4 # default is 3
    toc-title: Contents
    number-sections: true
    anchor-sections: true
    smooth-scroll: true
    highlight-style: dracula
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
nhmmscan_hits_filtered_merged_noOverlap_table_hits = readr::read_tsv("nhmmscan_hits_filtered_merged_noOverlap_table_hits.tab.txt")


nhmmscan_hits_filtered_merged_noOverlap_table_hits_qLens = nhmmscan_hits_filtered_merged_noOverlap_table_hits %>%
  select(`#chrom`, queryLen) %>%
  unique() %>%
  mutate(`#chrom` = factor(`#chrom`),
         shortName = factor(gsub("\\[.*\\]", "", `#chrom`)))

nhmmscan_hits_filtered_merged_noOverlap_table_hits = nhmmscan_hits_filtered_merged_noOverlap_table_hits %>%
  mutate(`#chrom` = factor(`#chrom`, levels = levels(nhmmscan_hits_filtered_merged_noOverlap_table_hits_qLens$`#chrom`)))%>%
  mutate(shortName = factor(gsub("\\[.*\\]", "", `#chrom`), levels = levels(nhmmscan_hits_filtered_merged_noOverlap_table_hits_qLens$shortName) ))

domainPlot = ggplot() +
  geom_rect(aes(
    xmin = 0,
    xmax = queryLen,
    ymin = as.numeric(shortName) - 0.0,
    ymax = as.numeric(shortName) + 0.4
  ),
  fill = "#2271B2", data = nhmmscan_hits_filtered_merged_noOverlap_table_hits_qLens) +

    geom_rect(aes(
    xmin = env_from,
    xmax = env_to,
    ymin = as.numeric(shortName) - 0.4,
    ymax = as.numeric(shortName) + 0.4,
    fill = name,
    env_from = env_from,
    env_to = env_to,
    env_len = env_len,
    acc = acc,
    evalue = evalue,
    score = score,
    scoreOverLen = scoreOverLen,
    hmm_covered = hmm_covered,
    hmm_From = hmm_From, 
    hmm_to = hmm_to,
    hmm_edges = hmm_edges

  ),
  color = "black",
  alpha = 0.75,
   data = nhmmscan_hits_filtered_merged_noOverlap_table_hits) +

  scale_fill_manual(values = scheme$hex(length(nhmmscan_hits_filtered_merged_noOverlap_table_hits$name))) +

  scale_y_continuous(
    labels = nhmmscan_hits_filtered_merged_noOverlap_table_hits_qLens$shortName,
    breaks = 1:nrow(nhmmscan_hits_filtered_merged_noOverlap_table_hits_qLens)
  ) +
  sofonias_theme
```


```{r}
#| column: screen-inset-shaded
create_dt(nhmmscan_hits_filtered_merged_noOverlap_table_hits)
```


```{r}
#| column: screen-inset-shaded
#| fig-height: !expr 'nrow(nhmmscan_hits_filtered_merged_noOverlap_table_hits) * 0.1'
ggplotly(domainPlot)
```
)";


int programWrapperRunner::runnhmmscan(const njh::progutils::CmdArgs & inputCommands){

	bfs::path subRegions;
	uint32_t extendSubRegions = 0;
	std::string defaultParameters = "--nonull2 --incT 50 --incdomT 50 -T 50 --notextw";
	bfs::path hmmModel;

  double fullDomainTrimHmmCoverageCutOff = 0.90;

	nhmmscanOutput::PostProcessHitsPars postProcessPars;
  postProcessPars.accCutOff = 0.80;
	postProcessPars.scoreCutOff = 200;
	postProcessPars.evalueCutOff = 1e-100;
	postProcessPars.scoreNormCutOff = 0.50;
	postProcessPars.softModelCovergeCutOff = 0.80;
	postProcessPars.hardScoreNormCutOff = 0.20;
	nhmmscanOutput::QueryResults::mergeOverlapingHitsPars mergePars;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processDefaultReader(VecStr{"fasta", "fastagz", "fastq", "fastqgz"}, true);
  setUp.setOption(fullDomainTrimHmmCoverageCutOff, "--fullDomainTrimHmmCoverageCutOff", "full Domain Trim Hmm Coverage Cut Off for when trimming from first to last determined region");

	setUp.setOption(defaultParameters, "--defaultParameters", "The default parameters given to hmmsearch");
	setUp.setOption(extendSubRegions, "--extendSubRegions", "Extend the sub regions by this much in both directions");

	setUp.setOption(hmmModel, "--hmmModel", "hmm model database, created by hmmbuild", true);
	postProcessPars.minOverlapFilt_ = 60;
	setUp.setOption(postProcessPars.minOverlapFilt_, "--minOverlapLenForFilt", "minimum Overlap length for the non-overlapping regions");

	setUp.setOption(postProcessPars.hmmStartFilter, "--hmmStartFilter", "Filter partial hmms domain hits if they start or end this far into the model");
	setUp.setOption(postProcessPars.minLength, "--minLength", "Minimum output domain hit length");
	setUp.setOption(subRegions, "--subRegions", "Run on the subregions as defined by this bed file, needs to have a 2bit file named with same prefix as the input file");

	setUp.setOption(postProcessPars.accCutOff, "--accCutOff", "soft accuracy cut off");
	setUp.setOption(postProcessPars.scoreCutOff, "--scoreCutOff", "soft score cut off");
	setUp.setOption(postProcessPars.evalueCutOff, "--evalueCutOff", "soft evalue cut off");
	setUp.setOption(postProcessPars.scoreNormCutOff, "--scoreNormCutOff", "soft scoreNorm cut off");

	setUp.setOption(postProcessPars.hardAccCutOff, "--hardAccCutOff", "hard accuracy cut off");
	if(postProcessPars.hardAccCutOff > postProcessPars.accCutOff){
		postProcessPars.accCutOff = postProcessPars.hardAccCutOff;
	}

	setUp.setOption(postProcessPars.hardScoreCutOff, "--hardScoreCutOff", "hard score cut off");
	if(postProcessPars.hardScoreCutOff > postProcessPars.scoreCutOff){
		postProcessPars.scoreCutOff = postProcessPars.hardScoreCutOff;
	}

	setUp.setOption(postProcessPars.hardEvalueCutOff, "--hardEvalueCutOff", "hard evalue cut off");
	if(postProcessPars.hardEvalueCutOff < postProcessPars.evalueCutOff){
		postProcessPars.evalueCutOff = postProcessPars.hardEvalueCutOff;
	}

	setUp.setOption(postProcessPars.hardScoreNormCutOff, "--hardScoreNormCutOff", "hard scoreNorm cut off");
	if(postProcessPars.hardScoreNormCutOff > postProcessPars.scoreNormCutOff){
		postProcessPars.scoreNormCutOff = postProcessPars.hardScoreNormCutOff;
	}

	setUp.setOption(postProcessPars.hardModelCovergeCutOff, "--hardModelCovergeCutOff", "hard mdoel coverage cut off");
	if(postProcessPars.hardModelCovergeCutOff > postProcessPars.softModelCovergeCutOff){
		postProcessPars.softModelCovergeCutOff = postProcessPars.hardModelCovergeCutOff;
	}
	setUp.processDirectoryOutputName(true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);

	//convert to fasta and non-gz hmm model

	auto inputSeqFnp = njh::files::make_path(setUp.pars_.directoryName_, "inputSeqs.fasta");

  auto noOverlapFiltAllModelsTrimFnp = njh::files::make_path(setUp.pars_.directoryName_, "noOverlapFiltAllModelsTrim.fasta");
  auto noOverlapFiltFnp = njh::files::make_path(setUp.pars_.directoryName_, "noOverlapFiltHits.fasta");

  auto noOverlapMergedFiltAllModelsTrimFnp = njh::files::make_path(setUp.pars_.directoryName_, "noOverlapMergedFiltAllModelsTrim.fasta");
	auto noOverlapMergedFiltFnp = njh::files::make_path(setUp.pars_.directoryName_, "noOverlapMergedFiltHits.fasta");
	auto hmmModelFnp = njh::files::make_path(setUp.pars_.directoryName_, "hmmModel.txt");

	std::unordered_map<std::string, GenomicRegion> regionsByName;
	std::unordered_map<std::string, uint32_t> realQueryLens;
	std::vector<GenomicRegion> regions;
	auto inputOpts = setUp.pars_.ioOptions_;
	if(!subRegions.empty()){
		auto twoBitFnp = bfs::path(setUp.pars_.ioOptions_.firstName_).replace_extension(".2bit");
		TwoBit::TwoBitFile tReader(twoBitFnp);
		realQueryLens = tReader.getSeqLens();
		auto regionsOutOpts = SeqIOOptions::genFastaOut(njh::files::make_path(setUp.pars_.directoryName_, "subRegions.fasta"));
		regions = bed3PtrsToGenomicRegs(getBed3s(subRegions));
		if(0 != extendSubRegions){
			njh::for_each(regions, [&realQueryLens,&extendSubRegions](GenomicRegion & region){
				BedUtility::extendLeftRight(region, extendSubRegions, extendSubRegions, njh::mapAt(realQueryLens, region.chrom_));
			});
			sortGRegionsByStart(regions);
			std::vector<GenomicRegion> mergedRegions;
			for(const auto & reg : regions){
				bool overlaps = false;
				if(!mergedRegions.empty()){
					if(mergedRegions.back().overlaps(reg, 1)) {
						overlaps = true;
						mergedRegions.back().end_ = std::max(mergedRegions.back().end_, reg.end_);
					}
				}
				if(!overlaps){
					mergedRegions.emplace_back(reg);
				}
			}
			OutputStream mergedBedOut(njh::files::make_path(setUp.pars_.directoryName_, "extendedSubRegions.bed"));
			njh::for_each(mergedRegions,[&mergedBedOut](GenomicRegion & region){
				region.uid_ = region.createUidFromCoords();
				mergedBedOut << region.genBedRecordCore().toDelimStr() << std::endl;
			});
			regions = mergedRegions;
		}


		SeqOutput writer(regionsOutOpts);
		writer.openOut();
		// uint32_t regionCount = 0;
		for(const auto & region : regions){
			regionsByName[region.uid_] = region;
			// ++regionCount;
			auto subRegion = region.extractSeq(tReader);
			subRegion.name_ = region.uid_;
			writer.write(subRegion);
		}
		inputOpts = SeqIOOptions::genFastaIn(regionsOutOpts.out_.outName());
	}

	auto seqKey = SeqIO::rewriteSeqsWithIndexAsName(
					inputOpts,
					SeqIOOptions::genFastaOut(inputSeqFnp),
					njh::files::make_path(setUp.pars_.directoryName_, "inputSeqNameKey.tab.txt"));

	//write out the model, this way you can supply it in gzip format but this will unzip it
	njh::files::reWriteFile(hmmModel, hmmModelFnp);
	nhmmscanOutput::run_hmmpress_ifNeed(hmmModelFnp);

	std::stringstream cmdSs;
	cmdSs << "nhmmscan " << defaultParameters
			<< " " << "--tblout " << "raw_all_hits_table.txt"
			<< " " << "hmmModel.txt"
			<< " " << "inputSeqs.fasta"
			<< " " << " > " << "nhmmscan_raw_output.txt";
	std::string cdCmd = "cd " + setUp.pars_.directoryName_ + " && ";
	auto cmdOutput = njh::sys::run( { cdCmd, cmdSs.str() });
	if (!cmdOutput.success_) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "failed to run hmmsearch "
				<< "\n";
		ss << cmdOutput.toJson() << "\n";
		throw std::runtime_error { ss.str() };
	} else {
		OutputStream nhmmscanOutput(
				njh::files::make_path(setUp.pars_.directoryName_,
						"nhmmscanCmdRunDetails.json"));
		nhmmscanOutput << cmdOutput.toJson() << std::endl;
	}



	auto nhmmscan_raw_outputFnp = njh::files::make_path(setUp.pars_.directoryName_, "nhmmscan_raw_output.txt");
	auto seqsWithNoDomainHitsFnp = njh::files::make_path(setUp.pars_.directoryName_, "seqsWithNoDomainHits.tab.txt");
	nhmmscanOutput outputParsed;
	if(regions.empty()){
		outputParsed = nhmmscanOutput::parseRawOutput(nhmmscan_raw_outputFnp, seqKey);
	}else{
		outputParsed = nhmmscanOutput::parseRawOutput(nhmmscan_raw_outputFnp, seqKey, regionsByName, realQueryLens);
	}
	auto postProcessResults = outputParsed.postProcessHits(postProcessPars);
	//convert hits table into a real table and
	outputParsed.writeInfoFiles(postProcessResults, setUp.pars_.directoryName_);

	if(setUp.pars_.debug_){
		//merging debugging
		for(const auto & filteredHits : postProcessResults.filteredHitsByQuery_){
			auto mergedResults = nhmmscanOutput::QueryResults::mergeOverlapingHits(filteredHits.second, mergePars);
			for(const auto & merged : mergedResults){
				if(merged.hits_.size() > 1){
					std::cout << merged.region_.genBedRecordCore().toDelimStrWithExtra() << std::endl;
					for(const auto & hit : merged.hits_){
						std::cout << "\t" << hit.genBed6_env().toDelimStrWithExtra() << std::endl;
					}
					std::cout << std::endl;
				}
			}
		}
	}

	//trim seqs to best overlapping positions sub seqs
	{
		VecStr seqsWithNoDomains;
		seqInfo seq;
		SeqInput reader(setUp.pars_.ioOptions_);

		SeqOutput noOverlapFiltWriter(SeqIOOptions::genFastaOut(noOverlapFiltFnp));
    SeqOutput noOverlapFiltAllModelsTrimWriter(SeqIOOptions::genFastaOut(noOverlapFiltAllModelsTrimFnp));

    SeqOutput noOverlapMergedFiltWriter(SeqIOOptions::genFastaOut(noOverlapMergedFiltFnp));
    SeqOutput noOverlapMergedFiltAllModelsTrimWriter(SeqIOOptions::genFastaOut(noOverlapMergedFiltAllModelsTrimFnp));


		reader.openIn();
		noOverlapFiltWriter.openOut();
		noOverlapMergedFiltWriter.openOut();
    noOverlapFiltAllModelsTrimWriter.openOut();
    noOverlapMergedFiltAllModelsTrimWriter.openOut();

		while(reader.readNextRead(seq)){

			if(njh::in(seq.name_, postProcessResults.filteredHitsMergedNonOverlapByQuery_)){
        std::vector<Bed6RecordCore> allRegions;

				for(const auto & hitGroup : postProcessResults.filteredHitsMergedNonOverlapByQuery_[seq.name_]){


          Bed6RecordCore region = hitGroup.region_.genBedRecordCore();
          for(const auto & hit : hitGroup.hits_){
            if(hit.modelCoverage() >= fullDomainTrimHmmCoverageCutOff){
              allRegions.emplace_back(region);
              allRegions.back().name_ = hit.targetName_;
              break;
            }
          }
					auto subSeq = seq.getSubRead(region.chromStart_, region.length());
					if(region.reverseStrand()){
						subSeq.reverseComplementRead(false, true);
					}
					MetaDataInName meta = hitGroup.genOutRegion().meta_;
					meta.addMeta("trimStart", region.chromStart_, true);
					meta.addMeta("trimEnd", region.chromEnd_, true);
					meta.addMeta("trimLen", region.length(), true);
					meta.addMeta("trimCov", region.length()/static_cast<double>(len(seq)), true);
					meta.addMeta("revStrand", region.reverseStrand(), true);
					meta.resetMetaInName(subSeq.name_);
					noOverlapMergedFiltWriter.write(subSeq);
				}

        if(!allRegions.empty()){
          BedUtility::coordSort(allRegions);
          auto start = allRegions.front().chromStart_;
          auto end = allRegions.back().chromEnd_;
          auto outSeq = seq.getSubRead(start, end - start);
          std::string models;
          for(const auto & r : allRegions){
            if(!models.empty()){
              models += "--";
            }
            models += r.name_;
          }
          MetaDataInName meta;
          if(MetaDataInName::nameHasMetaData(outSeq.name_)){
            meta = MetaDataInName(outSeq.name_);
          }
          meta.addMeta("models", models);
          if(meta.containsMeta("length")){
            meta.addMeta("length", len(outSeq), true);
          }
          meta.resetMetaInName(outSeq.name_);
          noOverlapMergedFiltAllModelsTrimWriter.write(outSeq);
        }
			}

			if(njh::in(seq.name_, postProcessResults.filteredNonOverlapHitsByQuery_)){
        std::vector<Bed6RecordCore> allRegions;
//        std::cout << seq.name_ << std::endl;
//        std::cout << "\tpostProcessResults.filteredHitsByQuery_[seq.name_].size():                 " << postProcessResults.filteredHitsByQuery_[seq.name_].size() << std::endl;
//        std::cout << "\tpostProcessResults.filteredNonOverlapHitsByQuery_[seq.name_].size():       " << postProcessResults.filteredNonOverlapHitsByQuery_[seq.name_].size() << std::endl;
//        std::cout << "\tpostProcessResults.filteredHitsMergedByQuery_[seq.name_].size():           " << postProcessResults.filteredHitsMergedByQuery_[seq.name_].size() << std::endl;
//        std::cout << "\tpostProcessResults.filteredHitsMergedNonOverlapByQuery_[seq.name_].size(): " << postProcessResults.filteredHitsMergedNonOverlapByQuery_[seq.name_].size() << std::endl;

				for(const auto & hit : postProcessResults.filteredNonOverlapHitsByQuery_[seq.name_]){
					Bed6RecordCore region = hit.genBed6_env();
          if(hit.modelCoverage() >=fullDomainTrimHmmCoverageCutOff){
            allRegions.emplace_back(region);
            allRegions.back().name_ = hit.targetName_;
          }

					auto subSeq = seq.getSubRead(region.chromStart_, region.length());
					if(region.reverseStrand()){
						subSeq.reverseComplementRead(false, true);
					}
					MetaDataInName meta;
					meta.addMeta("acc", hit.acc_, true);
					meta.addMeta("hmmFrom", hit.hmmFrom_, true);
					meta.addMeta("hmmTo", hit.hmmTo_, true);
					meta.addMeta("hmmCovered", hit.modelCoverage());
					meta.addMeta("trimStart", region.chromStart_, true);
					meta.addMeta("trimEnd", region.chromEnd_, true);
					meta.addMeta("trimLen", region.length(), true);
					meta.addMeta("trimCov", region.length()/static_cast<double>(len(seq)), true);
					meta.addMeta("revStrand", region.reverseStrand(), true);
					meta.addMeta("score", hit.modelScore_, true);
					meta.addMeta("scoreNorm", hit.modelScore_/region.length(), true);
					meta.addMeta("evalue", hit.modelEvalue_, true);
					meta.addMeta("model", hit.targetName_);
					meta.addMeta("ID", hit.targetDesc_);
					meta.resetMetaInName(subSeq.name_);
					noOverlapFiltWriter.write(subSeq);
				}
        if(!allRegions.empty()){
          BedUtility::coordSort(allRegions);
          auto start = allRegions.front().chromStart_;
          auto end = allRegions.back().chromEnd_;
          auto outSeq = seq.getSubRead(start, end - start);
          std::string models;
          for(const auto & r : allRegions){
            if(!models.empty()){
              models += "--";
            }
            models += r.name_;
          }
          MetaDataInName meta;
          if(MetaDataInName::nameHasMetaData(outSeq.name_)){
            meta = MetaDataInName(outSeq.name_);
          }
          meta.addMeta("models", models);
          if(meta.containsMeta("length")){
            meta.addMeta("length", len(outSeq), true);
          }
          meta.resetMetaInName(outSeq.name_);
          noOverlapFiltAllModelsTrimWriter.write(outSeq);
        }
			} else {
				seqsWithNoDomains.emplace_back(seq.name_);
			}
		}
		OutputStream noDomainHits(seqsWithNoDomainHitsFnp);
		for(const auto & name : seqsWithNoDomains){
			noDomainHits << name << std::endl;
		}
	}

	{
		OutputStream qmdPlotDomainsOut(njh::files::make_path(setUp.pars_.directoryName_, "plotDomainsHits.qmd"));
		qmdPlotDomainsOut << plotDomainsQmd << std::endl;
	}

	return 0;
}


} // namespace njhseq

