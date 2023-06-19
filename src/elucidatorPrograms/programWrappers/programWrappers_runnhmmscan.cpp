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



int programWrapperRunner::runnhmmscan(const njh::progutils::CmdArgs & inputCommands){

	bfs::path subRegions;
	uint32_t extendSubRegions = 0;
	std::string defaultParameters = "--nonull2 --incT 50 --incdomT 50 -T 50 --notextw";
	bfs::path hmmModel;
	nhmmscanOutput::PostProcessHitsPars postProcessPars;
	postProcessPars.accCutOff = 0.80;
	postProcessPars.scoreCutOff = 200;
	postProcessPars.evalueCutOff = 1e-100;
	postProcessPars.scoreNormCutOff = 0.50;
	nhmmscanOutput::QueryResults::mergeOverlapingHitsPars mergePars;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processDefaultReader(VecStr{"fasta", "fastagz", "fastq", "fastqgz"}, true);
	setUp.setOption(defaultParameters, "--defaultParameters", "The default parameters given to hmmsearch");
	setUp.setOption(extendSubRegions, "--extendSubRegions", "Extend the sub regions by this much in both directions");

	setUp.setOption(hmmModel, "--hmmModel", "hmm model database, created by hmmbuild", true);
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
	if(postProcessPars.hardEvalueCutOff > postProcessPars.evalueCutOff){
		postProcessPars.evalueCutOff = postProcessPars.hardEvalueCutOff;
	}
	setUp.setOption(postProcessPars.hardScoreNormCutOff, "--hardScoreNormCutOff", "hard scoreNorm cut off");
	if(postProcessPars.hardScoreNormCutOff > postProcessPars.scoreNormCutOff){
		postProcessPars.scoreNormCutOff = postProcessPars.hardScoreNormCutOff;
	}
	setUp.processDirectoryOutputName(true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);

	//convert to fasta and non-gz hmm model

	auto inputSeqFnp = njh::files::make_path(setUp.pars_.directoryName_, "inputSeqs.fasta");
	auto noOverlapFiltFnp = njh::files::make_path(setUp.pars_.directoryName_, "noOverlapFiltHits.fasta");
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
		uint32_t regionCount = 0;
		for(const auto & region : regions){
			regionsByName[region.uid_] = region;
			++regionCount;
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
		SeqOutput noOverlapMergedFiltWriter(SeqIOOptions::genFastaOut(noOverlapMergedFiltFnp));

		reader.openIn();
		noOverlapFiltWriter.openOut();
		noOverlapMergedFiltWriter.openOut();

		while(reader.readNextRead(seq)){
			if(njh::in(seq.name_, postProcessResults.filteredHitsMergedNonOverlapByQuery_)){
				for(const auto & hitGroup : postProcessResults.filteredHitsMergedNonOverlapByQuery_[seq.name_]){
					Bed6RecordCore region = hitGroup.region_.genBedRecordCore();
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
			}
			if(njh::in(seq.name_, postProcessResults.filteredNonOverlapHitsByQuery_)){
				for(const auto & hit : postProcessResults.filteredNonOverlapHitsByQuery_[seq.name_]){
					Bed6RecordCore region = hit.genBed6_env();
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
			} else {
				seqsWithNoDomains.emplace_back(seq.name_);
			}
		}
		OutputStream noDomainHits(seqsWithNoDomainHitsFnp);
		for(const auto & name : seqsWithNoDomains){
			noDomainHits << name << std::endl;
		}
	}

	return 0;
}


} // namespace njhseq

