/*
 * popGenExp_callVariantsAgainstRefSeq.cpp
 *
 *  Created on: Jun 23, 2021
 *      Author: nick
 */


#include "popGenExp.hpp"
#include "elucidator/objects/BioDataObject.h"
#include "elucidator/objects/counters/DNABaseCounter.hpp"
#include "elucidator/PopulationGenetics.h"

#include "elucidator/objects/seqContainers/CollapsedHaps.hpp"
#include <njhseq/concurrency/PairwisePairFactory.hpp>
#include <njhseq/GenomeUtils.h>
#include <njhseq/objects/Gene/GenomicAminoAcidPositionTyper.hpp>
#include <njhseq/objects/Gene/TranslatorByAlignment.hpp>


#include <njhseq/objects/seqObjects/Clusters/identicalCluster.hpp>


namespace njhseq {




int popGenExpRunner::callVariantsAgainstRefSeq(const njh::progutils::CmdArgs & inputCommands) {

	/**@todo should add the following 1) average pairwise difference, 2) make into a function, 3) generate connected hap map, 4) unique haps to region count, 5) doing multiple pop fields at once */

  TranslatorByAlignment::TranslatorByAlignmentPars transPars;
  TranslatorByAlignment::RunPars variantCallerRunPars;
  bfs::path knownAminoAcidChangesFnp;

	std::string identifier = "";
	variantCallerRunPars.lowVariantCutOff = 0.005;
	variantCallerRunPars.occurrenceCutOff = 1;

	bfs::path bedFnp = "";
	uint32_t outwardsExpand = 5;

	uint32_t numThreads = 1;
	bool getPairwiseComps = false;
	bfs::path metaFnp;
	bool keepNonFieldSamples = false;
	std::set<std::string> ignoreSubFields;
	MultiGenomeMapper::inputParameters gPars;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processReadInNames(true);
	setUp.processDirectoryOutputName(true);
	setUp.processAlnInfoInput();
	setUp.setOption(numThreads, "--numThreads", "number of threads");
	setUp.setOption(getPairwiseComps, "--getPairwiseComps", "get Pairwise comparison metrics");

	setUp.setOption(variantCallerRunPars.occurrenceCutOff, "--occurrenceCutOff", "Occurrence Cut Off, don't report variants below this count");
	setUp.setOption(variantCallerRunPars.lowVariantCutOff, "--lowVariantCutOff", "Low Variant Cut Off, don't report variants below this fraction");
	setUp.setOption(outwardsExpand, "--outwardsExpand", "The amount to expand outwards from given region when determining variants positions with extracted ref seq");
	setUp.setOption(metaFnp,    "--metaFnp",    "Meta data to add to sequences");

	setUp.setOption(bedFnp,    "--bed",    "A bed file of the location for the extraction", true);
	transPars.setOptions(setUp);

	gPars.genomeDir_ = njh::files::normalize(transPars.lzPars_.genomeFnp.parent_path());
	gPars.primaryGenome_ = bfs::basename(transPars.lzPars_.genomeFnp);
	gPars.primaryGenome_ = gPars.primaryGenome_.substr(0, gPars.primaryGenome_.rfind("."));
	gPars.selectedGenomes_ = {gPars.primaryGenome_};

	setUp.setOption(ignoreSubFields, "--ignoreSubFields", "Meta Sub Field values to ignore when calculating variants, e.g. --ignoreSubFields \"isFieldSample:TRUE,PreferredSample:FALSE\"");

	setUp.setOption(identifier, "--identifier", "Give a identifier name for info");
	setUp.setOption(knownAminoAcidChangesFnp, "--proteinMutantTypingFnp", "Protein Mutant Typing Fnp, columns should be ID=gene id in gff, AAPosition=amino acid position", false);
	setUp.setOption(keepNonFieldSamples, "--keepNonFieldSamples", "Keep Non Field Samples for population stats");
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);

	njh::files::checkExistenceThrow(bedFnp, __PRETTY_FUNCTION__);
	njh::files::checkExistenceThrow(transPars.lzPars_.genomeFnp, __PRETTY_FUNCTION__);
	std::unique_ptr<TranslatorByAlignment> translator;
	if("" != knownAminoAcidChangesFnp){
		if("" == transPars.gffFnp_){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "if supplying known amino acid positions than must also supply --gffFnp file"<< "\n";
			throw std::runtime_error{ss.str()};
		}
	}
	if("" != transPars.gffFnp_){
		translator = std::make_unique<TranslatorByAlignment>(transPars);
		if("" == transPars.lzPars_.genomeFnp){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << " if supplying gff file, must also supply --genomeFnp"<< "\n";
			throw std::runtime_error{ss.str()};
		}
	}
	std::unordered_map<std::string, std::set<uint32_t>> knownMutationsLocationsMap;
	if("" != knownAminoAcidChangesFnp){
		knownMutationsLocationsMap = TranslatorByAlignment::readInAAPositions(knownAminoAcidChangesFnp);
	}
	auto beds = getBeds(bedFnp);
	if(beds.empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error no records found in " << bedFnp << "\n";
		throw std::runtime_error{ss.str()};
	}
	MultiGenomeMapper gMapper(gPars);
	gMapper.loadInGenomes();
	gMapper.setUpGenomes();
	auto gPos = beds.front();

	//get identifer for sequences
	if("" == identifier){
		identifier = gPos->name_;
	}
	//get reference seq
	//get region and ref seq for mapping of variants
	OutputStream inputRegionBeforeExpandOut(OutOptions(njh::files::make_path(setUp.pars_.directoryName_, "inputRegionBeforeExpand.bed")));
	inputRegionBeforeExpandOut << gPos->toDelimStrWithExtra() << std::endl;

	uint32_t oldLen = gPos->length();
	BedUtility::extendLeftRight(*gPos, outwardsExpand, outwardsExpand,
			gMapper.genomes_.at(gPars.primaryGenome_)->chromosomeLengths_.at(
					gPos->chrom_));
	if(oldLen == gPos->score_){
		gPos->score_ = gPos->length();
	}

	OutputStream bedOut(OutOptions(njh::files::make_path(setUp.pars_.directoryName_, "inputRegion.bed")));
	bedOut << gPos->toDelimStrWithExtra() << std::endl;

	GenomicRegion refRegion(*gPos);
	refRegion.reverseSrand_ = false;
	TwoBit::TwoBitFile tReader(gMapper.genomes_.at(gPars.primaryGenome_)->fnpTwoBit_);
	auto refSeq = refRegion.extractSeq(tReader);

	SeqOutput::write(std::vector<seqInfo> { refSeq },
			SeqIOOptions::genFastaOut(
					njh::files::make_path(setUp.pars_.directoryName_,
							"inputRegion.fasta")));

	//read in meta if available
	std::unique_ptr<MultipleGroupMetaData> meta;
	if("" != metaFnp){
		meta = std::make_unique<MultipleGroupMetaData>(metaFnp);
	}

	std::unordered_map<std::string, std::set<std::string>> metaValuesToAvoid = njh::progutils::CmdArgs::sepSubArgsMulti<std::string, std::string>(ignoreSubFields);


	auto inputSeqs = CollapsedHaps::readInReads(setUp.pars_.ioOptions_, meta, metaValuesToAvoid);
	uint64_t maxLen = readVec::getMaxLength(inputSeqs.seqs_);
	inputSeqs.setFrequencies();
	if (gPos->reverseStrand()) {
		inputSeqs.revCompSeqs();
	}
	//samples names
	auto sampNamesPerSeq = inputSeqs.getSampleNamesPerSeqs();
	auto allSamples = inputSeqs.getAllSampleNames();
	//rename based on freq
	std::vector<uint32_t> orderByCnt = inputSeqs.getOrderByTopCnt();
	uint32_t seqId = 0;
	for(const auto pos : orderByCnt){
		MetaDataInName popMeta;
		popMeta.addMeta("HapPopUIDCount", static_cast<uint32_t>(std::round(inputSeqs.seqs_[pos]->cnt_)));
		inputSeqs.seqs_[pos]->name_ = njh::pasteAsStr(identifier, ".", leftPadNumStr<uint32_t>(seqId, inputSeqs.size()),popMeta.createMetaName());
		++seqId;
	}
	readVec::getMaxLength(refSeq, maxLen);
	aligner alignerObj(maxLen, gapScoringParameters(5,1,0,0,0,0), substituteMatrix(2,-2), false);
	alignerObj.weighHomopolymers_ = false;
	alignerObj.processAlnInfoInput(setUp.pars_.alnInfoDirName_, setUp.pars_.verbose_);
	//set up variant info
	auto idSeq = refSeq;
	idSeq.name_ = identifier;
	TranslatorByAlignment::VariantsInfo varInfo(refRegion.genBed3RecordCore(), idSeq);
	//get variant info
	auto refComps = inputSeqs.getCompsAgainstRef(refSeq, alignerObj, numThreads);
	for(const auto pos :  iter::range(inputSeqs.size())){
		varInfo.addVariantInfo(
				refComps[pos].refAlnSeq_,
				refComps[pos].queryAlnSeq_,
				inputSeqs.seqs_[pos]->cnt_,
				sampNamesPerSeq[pos],
				refComps[pos].comp_,
				refRegion.start_);
	}
	varInfo.setFinals(variantCallerRunPars);
	//variable location
	std::string variableRegionID = "";
	if(!varInfo.variablePositons_.empty()){
		GenomicRegion variableRegion(varInfo.getVariableRegion());
		OutputStream bedVariableRegionOut(OutOptions(njh::files::make_path(setUp.pars_.directoryName_, "variableRegion.bed")));
		bedVariableRegionOut << variableRegion.genBedRecordCore().toDelimStrWithExtra() << std::endl;
		variableRegionID = njh::pasteAsStr(variableRegion.chrom_, "::", variableRegion.start_, "::", variableRegion.end_);
	}
	auto variantCallsDir = njh::files::make_path(setUp.pars_.directoryName_, "variantCalls");
	njh::files::makeDir(njh::files::MkdirPar{variantCallsDir});
	if(!varInfo.snpsFinal.empty() || ! varInfo.deletionsFinal.empty() || !varInfo.insertionsFinal.empty()){
		varInfo.writeVCF(njh::files::make_path(variantCallsDir, "allVariants.vcf"));
	}

	{
		OutputStream divMeasuresOut(njh::files::make_path(setUp.pars_.directoryName_, "divMeasures.tab.txt"));
		CollapsedHaps::AvgPairwiseMeasures avgPMeasures;
		if(getPairwiseComps && inputSeqs.size() > 1 ){
			auto allComps = inputSeqs.getPairwiseComps(alignerObj, numThreads);
			avgPMeasures = inputSeqs.getAvgPairwiseMeasures(allComps);
		}
		divMeasuresOut << "id\tname\ttotalHaplotypes\tuniqueHaplotypes\tsinglets\tdoublets\texpShannonEntropy\tShannonEntropyE\teffectiveNumOfAlleles\the\tlengthPolymorphism" ;
		if(getPairwiseComps){
			divMeasuresOut << "\tavgPercentID\tavgNumOfDiffs";
			divMeasuresOut << "\tnSegratingSites";

			divMeasuresOut << "\tTajimaD\tTajimaDPVal";
		}
		divMeasuresOut << std::endl;
		std::unordered_map<uint32_t, uint32_t> readLens = inputSeqs.getReadLenMap();
		inputSeqs.setFrequencies();
		auto divMeasures = PopGenCalculator::getGeneralMeasuresOfDiversity(inputSeqs.seqs_);
		divMeasuresOut << identifier
				<< "\t" << bfs::basename(setUp.pars_.ioOptions_.firstName_)
				<< "\t" << inputSeqs.getTotalHapCount()
				<< "\t" << inputSeqs.seqs_.size()
				<< "\t" << divMeasures.singlets_
				<< "\t" << divMeasures.doublets_
				<< "\t" << divMeasures.expShannonEntropy_
				<< "\t" << divMeasures.ShannonEntropyE_
				<< "\t" << divMeasures.effectiveNumOfAlleles_
				<< "\t" << divMeasures.heterozygostiy_
				<< "\t" << (readLens.size() > 1 ? "true" : "false");
		if(getPairwiseComps){
			if(inputSeqs.size() > 1){
				divMeasuresOut << "\t" << avgPMeasures.avgPercentId
										<< "\t" << avgPMeasures.avgNumOfDiffs
										<< "\t" << varInfo.getFinalNumberOfSegratingSites();
				if(0 == varInfo.getFinalNumberOfSegratingSites()){
					divMeasuresOut
							<< "\t" << 0
							<< "\t" << 1;
				}else{
					try {
						auto tajimad = PopGenCalculator::calcTajimaTest(inputSeqs.getTotalHapCount(), varInfo.getFinalNumberOfSegratingSites(), avgPMeasures.avgNumOfDiffs);
						divMeasuresOut
								<< "\t" << tajimad.d_
								<< "\t" << tajimad.pval_beta_;
					} catch (std::exception & e) {
						divMeasuresOut
								<< "\t" << "NA"
								<< "\t" << "NA";
					}

				}
			}else{
				divMeasuresOut << "\t" << 1
										<< "\t" << avgPMeasures.avgNumOfDiffs
										<< "\t" << 0;
				divMeasuresOut
						<< "\t" << 0
						<< "\t" << 1;
			}
		}
		divMeasuresOut << std::endl;
	}
	//std::cout << __FILE__ << " " << __LINE__ << std::endl;

	{
		//write out seqs
		auto uniqueSeqsOpts = SeqIOOptions::genFastaOutGz(njh::files::make_path(setUp.pars_.directoryName_, "uniqueSeqs.fasta.gz"));
		{
			SeqOutput uniqueWriter(uniqueSeqsOpts);
			uniqueWriter.openOut();
			for(const auto pos : orderByCnt){
				uniqueWriter.write(inputSeqs.seqs_[pos]);
			}
			uniqueWriter.closeOut();
		}

		OutputStream nameOut(OutOptions(njh::files::make_path(setUp.pars_.directoryName_, "uniqueSeqs_names.tab.txt.gz")));
		nameOut << "name\tnumber\tinputNames"	<< std::endl;
		OutputStream metaLabNamesOut(OutOptions(njh::files::make_path(setUp.pars_.directoryName_, "uniqueSeqs_nonFieldSampleNames.tab.txt.gz")));
		metaLabNamesOut << "name\tsamples" << std::endl;

		for(const auto pos : orderByCnt){
			nameOut << inputSeqs.seqs_[pos]->name_
					<< "\t" << inputSeqs.seqs_[pos]->cnt_
					<< "\t" << njh::conToStr(inputSeqs.names_[pos], ",") << std::endl;
			VecStr nonFieldSampleNames{};
			for(const auto & name : inputSeqs.names_[pos]){
				if(MetaDataInName::nameHasMetaData(name)){
					MetaDataInName meta(name);
					std::string sampleField = "sample";
					if(meta.containsMeta("BiologicalSample")){
						sampleField = "BiologicalSample";
					}
					if(meta.containsMeta(sampleField) && meta.containsMeta("IsFieldSample") && !meta.getMeta<bool>("IsFieldSample")){
						if(meta.containsMeta("site") && "LabIsolate" == meta.getMeta("site")){
							nonFieldSampleNames.emplace_back(meta.getMeta(sampleField));
						}
					}
				}
			}
			metaLabNamesOut << inputSeqs.seqs_[pos]->name_ << "\t" << njh::conToStr(nonFieldSampleNames, ",") << std::endl;
		}

		std::map<std::string, std::map<std::string, MetaDataInName>> knownAAMeta;
		std::map<std::string, MetaDataInName> fullAATyped;
		std::map<std::string, std::vector<TranslatorByAlignment::AAInfo>> fullAATypedWithCodonInfo;
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		if("" != transPars.gffFnp_){
			//std::cout << __FILE__ << " " << __LINE__ << std::endl;
			auto uniqueSeqInOpts = SeqIOOptions::genFastaIn(uniqueSeqsOpts.out_.outName());
			auto variantInfoDir =  njh::files::make_path(setUp.pars_.directoryName_, "proteinVariantInfo");
			//std::cout << __FILE__ << " " << __LINE__ << std::endl;
			njh::files::makeDir(njh::files::MkdirPar{variantInfoDir});
			translator->pars_.keepTemporaryFiles_ = true;
			translator->pars_.workingDirtory_ = variantInfoDir;
			//std::cout << __FILE__ << " " << __LINE__ << std::endl;
			std::unordered_map<std::string, std::unordered_set<std::string>> sampNamesForPopHaps;
			for(const auto pos : iter::range(inputSeqs.size())){
				sampNamesForPopHaps[inputSeqs.seqs_[pos]->name_] = sampNamesPerSeq[pos];
			}
			//std::cout << __FILE__ << " " << __LINE__ << std::endl;
			auto translatedRes = translator->run(uniqueSeqInOpts, sampNamesForPopHaps, variantCallerRunPars);
			//std::cout << __FILE__ << " " << __LINE__ << std::endl;
			OutputStream popBedLocs(njh::files::make_path(variantInfoDir, "inputSeqs.bed"));
			translatedRes.writeSeqLocations(popBedLocs);
			if(!translatedRes.translations_.empty()){
				//std::cout << __FILE__ << " " << __LINE__ << std::endl;
				SeqOutput transwriter(SeqIOOptions::genFastaOut(njh::files::make_path(variantInfoDir, "translatedInput.fasta")));
				for(const auto & seqName : translatedRes.translations_){
					for(const auto & transcript : seqName.second){
						transwriter.openWrite(transcript.second.translation_);
					}
				}

				OutputStream transBedLocs(njh::files::make_path(variantInfoDir, "translatedInput.bed"));
				translatedRes.writeSeqLocationsTranslation(transBedLocs);
				std::unordered_map<std::string, std::set<uint32_t>> knownAAMutsChromPositions;
				//std::cout << __FILE__ << " " << __LINE__ << std::endl;
				{
					//protein
					for(auto & varPerTrans : translatedRes.proteinVariants_){
						//std::cout << __FILE__ << " " << __LINE__ << std::endl;

						varPerTrans.second.writeVCF(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerTrans.first +  "-protein.vcf")));
						//std::cout << __FILE__ << " " << __LINE__ << std::endl;
						varPerTrans.second.writeOutSNPsFinalInfo(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerTrans.first +  "-protein_aminoAcidVariable.tab.txt")), varPerTrans.first, true);
						//std::cout << __FILE__ << " " << __LINE__ << std::endl;
						std::set<uint32_t> knownMutationsLocations;
						for(const auto & snpPos : varPerTrans.second.allBases){
							if(njh::in(snpPos.first + 1, knownMutationsLocationsMap[varPerTrans.first])){
								knownMutationsLocations.emplace(snpPos.first);
								auto genomicLocationForAAPos = translatedRes.translationInfoForTranscirpt_.at(varPerTrans.first)->genBedFromAAPositions(snpPos.first, snpPos.first + 1);
								for(const auto gPos : iter::range(genomicLocationForAAPos.chromStart_, genomicLocationForAAPos.chromEnd_)){
									knownAAMutsChromPositions[genomicLocationForAAPos.chrom_].emplace(gPos);
								}
							}
						}
						//std::cout << __FILE__ << " " << __LINE__ << std::endl;
						varPerTrans.second.writeOutSNPsAllInfo(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerTrans.first +  "-protein_aminoAcidsAll.tab.txt")), varPerTrans.first, true);
						//std::cout << __FILE__ << " " << __LINE__ << std::endl;
						if(!varPerTrans.second.variablePositons_.empty()){
							GenomicRegion variableRegion = varPerTrans.second.getVariableRegion();
							variableRegion.start_ += 1; //do one based positioning
							OutputStream bedVariableRegionOut(OutOptions(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerTrans.first +  "-protein_variableRegion.bed"))));
							bedVariableRegionOut << variableRegion.genBedRecordCore().toDelimStrWithExtra() << std::endl;
						}
						std::set<uint32_t> allLocations(knownMutationsLocations.begin(), knownMutationsLocations.end());
						for(const auto & variablePos : varPerTrans.second.snpsFinal){
							allLocations.emplace(variablePos.first);
						}
						//std::cout << __FILE__ << " " << __LINE__ << std::endl;
						for (auto & seqName : translatedRes.translations_) {
							if (njh::in(varPerTrans.first, seqName.second)) {
								std::string popName = seqName.first.substr(0, seqName.first.rfind("_f"));
								std::string transcript = varPerTrans.first;
								for (const auto & loc : allLocations) {
									//location is not within the aligned translation
									if(loc < std::get<0>(seqName.second[varPerTrans.first].firstAminoInfo_).aaPos_ || loc > std::get<0>(seqName.second[varPerTrans.first].lastAminoInfo_).aaPos_){
										continue;
									}

									auto codon = seqName.second[varPerTrans.first].getCodonForAARefPos(loc);

									fullAATyped[popName].addMeta(njh::pasteAsStr(varPerTrans.first, "-", loc + 1), codon.aa_);
									fullAATypedWithCodonInfo[popName].emplace_back(
											TranslatorByAlignment::AAInfo(varPerTrans.first, loc, codon,
													njh::in(loc, knownMutationsLocations)));
								}
							}
						}
						if(!knownMutationsLocations.empty()){
							varPerTrans.second.writeOutSNPsInfo(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerTrans.first +  "-protein_aminoAcidKnownMutations.tab.txt")), varPerTrans.first, knownMutationsLocations, true);
						}
					}
				}
				//std::cout << __FILE__ << " " << __LINE__ << std::endl;
				{
					//snps
					for( auto & varPerChrom : translatedRes.seqVariants_){
						//std::cout << __FILE__ << " " << __LINE__ << std::endl;
						varPerChrom.second.writeVCF(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerChrom.first +  "-genomic.vcf")));
						//std::cout << __FILE__ << " " << __LINE__ << std::endl;
						varPerChrom.second.writeOutSNPsFinalInfo(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerChrom.first +  "-SNPs.tab.txt")), varPerChrom.first);
						//std::cout << __FILE__ << " " << __LINE__ << std::endl;
						if(!knownAAMutsChromPositions[varPerChrom.first].empty()){
							varPerChrom.second.writeOutSNPsInfo(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerChrom.first +  "-knownAA_SNPs.tab.txt")), varPerChrom.first, knownAAMutsChromPositions[varPerChrom.first], true);
						}
						//std::cout << __FILE__ << " " << __LINE__ << std::endl;
						varPerChrom.second.writeOutSNPsAllInfo(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerChrom.first +  "-allBases.tab.txt")), varPerChrom.first);
						//std::cout << __FILE__ << " " << __LINE__ << std::endl;
						if(!varPerChrom.second.variablePositons_.empty()){
							//std::cout << __FILE__ << " " << __LINE__ << std::endl;
							GenomicRegion variableRegion = varPerChrom.second.getVariableRegion();
							//std::cout << __FILE__ << " " << __LINE__ << std::endl;
							OutputStream bedVariableRegionOut(OutOptions(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerChrom.first +  "-chromosome_variableRegion.bed"))));
							//std::cout << __FILE__ << " " << __LINE__ << std::endl;
							bedVariableRegionOut << variableRegion.genBedRecordCore().toDelimStrWithExtra() << std::endl;
						}
						//std::cout << __FILE__ << " " << __LINE__ << std::endl;
					}
				}
			}
			//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		}
	}


	//std::cout << __FILE__ << " " << __LINE__ << std::endl;
	alignerObj.processAlnInfoOutput(setUp.pars_.outAlnInfoDirName_, setUp.pars_.verbose_);
//	std::cout << "varInfo.snpsFinal.size(): " << varInfo.snpsFinal.size() << std::endl;
//	std::cout << "varInfo.deletionsFinal.size(): " << varInfo.deletionsFinal.size() << std::endl;
//	std::cout << "varInfo.insertionsFinal.size(): " << varInfo.insertionsFinal.size() << std::endl;
//


	return 0;
}



} //namespace njhseq

