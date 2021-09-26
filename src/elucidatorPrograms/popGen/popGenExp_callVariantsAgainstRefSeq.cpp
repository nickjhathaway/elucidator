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


int popGenExpRunner::callVariantsAgainstRefGenome(const njh::progutils::CmdArgs & inputCommands) {

	/**@todo should add the following  2) make into a function, 3) generate connected hap map, 4) unique haps to region count, 5) doing multiple pop fields at once */

  TranslatorByAlignment::TranslatorByAlignmentPars transPars;
  TranslatorByAlignment::RunPars variantCallerRunPars;
  CollapsedHaps::GenPopMeasuresPar calcPopMeasuresPars;


	std::string identifier = "";
	variantCallerRunPars.lowVariantCutOff = 0.005;
	variantCallerRunPars.occurrenceCutOff = 1;

	uint32_t numThreads = 1;
	bool noDiagAlnPairwiseComps = false;


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
	calcPopMeasuresPars.numThreads = numThreads;
	setUp.setOption(calcPopMeasuresPars.getPairwiseComps, "--getPairwiseComps", "get Pairwise comparison metrics");
	setUp.setOption(noDiagAlnPairwiseComps, "--noDiagAlnPairwiseComps", "Use diagonal Alignment for Pairwise Comparisons");
	calcPopMeasuresPars.diagAlnPairwiseComps = !noDiagAlnPairwiseComps;

	setUp.setOption(variantCallerRunPars.occurrenceCutOff, "--occurrenceCutOff", "Occurrence Cut Off, don't report variants below this count");
	setUp.setOption(variantCallerRunPars.lowVariantCutOff, "--lowVariantCutOff", "Low Variant Cut Off, don't report variants below this fraction");
	calcPopMeasuresPars.lowVarFreq = variantCallerRunPars.lowVariantCutOff;

	setUp.setOption(metaFnp,    "--metaFnp",    "Meta data to add to sequences");
	transPars.setOptions(setUp, true);

	gPars.genomeDir_ = njh::files::normalize(transPars.lzPars_.genomeFnp.parent_path());
	gPars.primaryGenome_ = bfs::basename(transPars.lzPars_.genomeFnp);
	gPars.primaryGenome_ = gPars.primaryGenome_.substr(0, gPars.primaryGenome_.rfind("."));
	gPars.selectedGenomes_ = {gPars.primaryGenome_};

	setUp.setOption(ignoreSubFields, "--ignoreSubFields", "Meta Sub Field values to ignore when calculating variants, e.g. --ignoreSubFields \"isFieldSample:TRUE,PreferredSample:FALSE\"");

	setUp.setOption(identifier, "--identifier", "Give a identifier name for info", true);
	setUp.setOption(transPars.knownAminoAcidMutationsFnp_, "--proteinMutantTypingFnp", "Protein Mutant Typing Fnp, columns should be ID=gene id in gff, AAPosition=amino acid position", false);
	setUp.setOption(keepNonFieldSamples, "--keepNonFieldSamples", "Keep Non Field Samples for population stats");
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);



	njh::files::checkExistenceThrow(transPars.lzPars_.genomeFnp, __PRETTY_FUNCTION__);
	std::unique_ptr<TranslatorByAlignment> translator = std::make_unique<TranslatorByAlignment>(transPars);

	MultiGenomeMapper gMapper(gPars);
	gMapper.loadInGenomes();
	gMapper.setUpGenomes();

	//read in meta if available
	std::unique_ptr<MultipleGroupMetaData> meta;
	if("" != metaFnp){
		meta = std::make_unique<MultipleGroupMetaData>(metaFnp);
	}

	std::unordered_map<std::string, std::set<std::string>> metaValuesToAvoid = njh::progutils::CmdArgs::sepSubArgsMulti<std::string, std::string>(ignoreSubFields);
	auto inputSeqs = CollapsedHaps::readInReads(setUp.pars_.ioOptions_, meta, metaValuesToAvoid);
	inputSeqs.setFrequencies();
	//samples names
	auto sampNamesPerSeq = inputSeqs.getSampleNamesPerSeqs();
	auto allSamples = inputSeqs.getAllSampleNames();
	//rename based on freq
	inputSeqs.renameBaseOnFreq(identifier);

	uint64_t maxLen = readVec::getMaxLength(inputSeqs.seqs_);

	std::shared_ptr<aligner> alignerObj = std::make_shared<aligner>(maxLen, gapScoringParameters(5,1,0,0,0,0), substituteMatrix(2,-2), false);
	alignerObj->weighHomopolymers_ = false;
	alignerObj->processAlnInfoInput(setUp.pars_.alnInfoDirName_, setUp.pars_.verbose_);


	//std::cout << __FILE__ << " " << __LINE__ << std::endl;
	std::unordered_map<std::string, uint32_t> seqNameKey = inputSeqs.genSeqNameKey();
	//write out seqs
	auto uniqueSeqsOpts = SeqIOOptions::genFastaOutGz(njh::files::make_path(setUp.pars_.directoryName_, "uniqueSeqs.fasta.gz"));
	inputSeqs.writeOutSeqsOrdCnt(uniqueSeqsOpts);
	inputSeqs.writeNames(OutOptions(njh::files::make_path(setUp.pars_.directoryName_, "uniqueSeqs_names.tab.txt.gz")));
	inputSeqs.writeOutMetaFields(OutOptions(njh::files::make_path(setUp.pars_.directoryName_, "uniqueSeqs_meta.tab.txt.gz")));
	inputSeqs.writeLabIsolateNames(OutOptions(njh::files::make_path(setUp.pars_.directoryName_, "uniqueSeqs_labIsolateNames.tab.txt.gz")));


	auto uniqueSeqInOpts = SeqIOOptions::genFastaIn(uniqueSeqsOpts.out_.outName());
	auto variantInfoDir =  njh::files::make_path(setUp.pars_.directoryName_, "variantCalls");
	njh::files::makeDir(njh::files::MkdirPar{variantInfoDir});
	translator->pars_.keepTemporaryFiles_ = true;
	translator->pars_.workingDirtory_ = variantInfoDir;
	std::unordered_map<std::string, std::unordered_set<std::string>> sampNamesForPopHaps;
	for(const auto pos : iter::range(inputSeqs.size())){
		sampNamesForPopHaps[inputSeqs.seqs_[pos]->name_] = sampNamesPerSeq[pos];
	}
	auto translatedRes = translator->run(uniqueSeqInOpts, sampNamesForPopHaps, variantCallerRunPars);
	OutputStream popBedLocs(njh::files::make_path(variantInfoDir, "inputSeqs.bed"));
	translatedRes.writeOutSeqAlnIndvVars(njh::files::make_path(variantInfoDir, "variantsPerSeqAln.tab.txt"));
	translatedRes.writeSeqLocations(popBedLocs);
	std::unordered_map<std::string, std::set<uint32_t>> knownAAMutsChromPositions;

	if(!translatedRes.translations_.empty()){
		SeqOutput transwriter(SeqIOOptions::genFastaOut(njh::files::make_path(variantInfoDir, "translatedInput.fasta")));
		std::vector<seqInfo> translatedSeqs;
		std::vector<std::unordered_set<std::string>> translatedSeqInputNames;
		auto seqNames = njh::getVecOfMapKeys(translatedRes.translations_);
		njh::sort(seqNames);
		for(const auto & seqName : seqNames){
			for(const auto & transcript : translatedRes.translations_.at(seqName)){
				transwriter.openWrite(transcript.second.translation_);
				translatedSeqs.emplace_back(transcript.second.translation_);
				translatedSeqs.back().cnt_ = inputSeqs.names_[seqNameKey[seqName]].size();
				translatedSeqInputNames.emplace_back(inputSeqs.names_[seqNameKey[seqName]]);
			}
		}
		auto inputTranslatedSeq = CollapsedHaps::collapseReads(translatedSeqs, translatedSeqInputNames);
		inputTranslatedSeq.renameBaseOnFreq(identifier);
		//write out seqs
		auto uniqueSeqsOpts = SeqIOOptions::genFastaOutGz(njh::files::make_path(variantInfoDir, "uniqueTranslatedSeqs.fasta.gz"));
		inputTranslatedSeq.writeOutSeqsOrdCnt(uniqueSeqsOpts);
		inputTranslatedSeq.writeNames(OutOptions(njh::files::make_path(variantInfoDir, "uniqueTranslatedSeqs_names.tab.txt.gz")));
		inputTranslatedSeq.writeOutMetaFields(OutOptions(njh::files::make_path(variantInfoDir, "uniqueTranslatedSeqs_meta.tab.txt.gz")));
		inputTranslatedSeq.writeLabIsolateNames(OutOptions(njh::files::make_path(variantInfoDir, "uniqueTranslatedSeqs_labIsolateNames.tab.txt.gz")));
		calcPopMeasuresPars.numSegSites_ =
				translatedRes.proteinVariants_.begin()->second.getFinalNumberOfSegratingSites();
		auto divMeasures = inputTranslatedSeq.getGeneralMeasuresOfDiversity(
				calcPopMeasuresPars, alignerObj);
		divMeasures.writeDivMeasures(
				njh::files::make_path(variantInfoDir,
						"translatedDivMeasures.tab.txt"), inputTranslatedSeq,
				identifier, calcPopMeasuresPars);

		OutputStream transBedLocs(njh::files::make_path(variantInfoDir, "translatedInput.bed"));
		translatedRes.writeSeqLocationsTranslation(transBedLocs);
		translatedRes.writeOutTranslatedIndvVars(njh::files::make_path(variantInfoDir, "variantsPerTranslatedSeq.tab.txt"));
		{
			//protein
			for(auto & varPerTrans : translatedRes.proteinVariants_){
				varPerTrans.second.writeVCF(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerTrans.first +  "-protein.vcf")));
				varPerTrans.second.writeOutSNPsFinalInfo(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerTrans.first +  "-protein_aminoAcidVariable.tab.txt")), varPerTrans.first, true);
				std::set<uint32_t> knownMutationsLocationsZeroBased;
				for(const auto & snpPos : varPerTrans.second.allBases){
					if(njh::in(snpPos.first + 1, translator->knownAminoAcidPositions_[varPerTrans.first])){
						knownMutationsLocationsZeroBased.emplace(snpPos.first);
						auto genomicLocationForAAPos = translatedRes.translationInfoForTranscirpt_.at(varPerTrans.first)->genBedFromAAPositions(snpPos.first, snpPos.first + 1);
						for(const auto gPos : iter::range(genomicLocationForAAPos.chromStart_, genomicLocationForAAPos.chromEnd_)){
							knownAAMutsChromPositions[genomicLocationForAAPos.chrom_].emplace(gPos);
						}
					}
				}
				varPerTrans.second.writeOutSNPsAllInfo(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerTrans.first +  "-protein_aminoAcidsAll.tab.txt")), varPerTrans.first, true);
				if(!varPerTrans.second.variablePositons_.empty()){
					GenomicRegion variableRegion = varPerTrans.second.getVariableRegion();
					variableRegion.start_ += 1; //do one based positioning
					OutputStream bedVariableRegionOut(OutOptions(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerTrans.first +  "-protein_variableRegion.bed"))));
					bedVariableRegionOut << variableRegion.genBedRecordCore().toDelimStrWithExtra() << std::endl;
				}
				if(!knownMutationsLocationsZeroBased.empty()){
					varPerTrans.second.writeOutSNPsInfo(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerTrans.first +  "-protein_aminoAcidKnownMutations.tab.txt")), varPerTrans.first, knownMutationsLocationsZeroBased, true);
				}
			}
		}
	}


	//snps
	uint32_t maxSeqCount = 0;

	for(auto & varPerChrom : translatedRes.seqVariants_){
		for(const auto & count : varPerChrom.second.depthPerPosition){
			if(count.second > maxSeqCount){
				//cheap way of doing this for now
				maxSeqCount = count.second;
				calcPopMeasuresPars.numSegSites_ = varPerChrom.second.getFinalNumberOfSegratingSites();
			}
		}
		varPerChrom.second.writeVCF(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerChrom.first +  "-genomic.vcf")));
		varPerChrom.second.writeOutSNPsFinalInfo(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerChrom.first +  "-SNPs.tab.txt")), varPerChrom.first);
		if(!knownAAMutsChromPositions[varPerChrom.first].empty()){
			varPerChrom.second.writeOutSNPsInfo(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerChrom.first +  "-knownAA_SNPs.tab.txt")), varPerChrom.first, knownAAMutsChromPositions[varPerChrom.first]);
		}
		varPerChrom.second.writeOutSNPsAllInfo(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerChrom.first +  "-allBases.tab.txt")), varPerChrom.first);
		if(!varPerChrom.second.variablePositons_.empty()){
			GenomicRegion variableRegion = varPerChrom.second.getVariableRegion();
			OutputStream bedVariableRegionOut(OutOptions(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerChrom.first +  "-chromosome_variableRegion.bed"))));
			bedVariableRegionOut << variableRegion.genBedRecordCore().toDelimStrWithExtra() << std::endl;
		}
	}

	auto divMeasures = inputSeqs.getGeneralMeasuresOfDiversity(
			calcPopMeasuresPars, alignerObj);
	divMeasures.writeDivMeasures(
			njh::files::make_path(setUp.pars_.directoryName_,
					"divMeasures.tab.txt"), inputSeqs, identifier, calcPopMeasuresPars);

	alignerObj->processAlnInfoOutput(setUp.pars_.outAlnInfoDirName_, setUp.pars_.verbose_);


	return 0;
}




int popGenExpRunner::callVariantsAgainstRefSeq(const njh::progutils::CmdArgs & inputCommands) {

	/**@todo should add the following 1) average pairwise difference, 2) make into a function, 3) generate connected hap map, 4) unique haps to region count, 5) doing multiple pop fields at once */
	bfs::path genomeFnp = "";
  TranslatorByAlignment::RunPars variantCallerRunPars;
  CollapsedHaps::GenPopMeasuresPar calcPopMeasuresPars;


	std::string identifier = "";
	variantCallerRunPars.lowVariantCutOff = 0.005;
	variantCallerRunPars.occurrenceCutOff = 1;

	bfs::path bedFnp = "";
	uint32_t outwardsExpand = 5;

	uint32_t numThreads = 1;
	bool noDiagAlnPairwiseComps = false;


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
	calcPopMeasuresPars.numThreads = numThreads;
	setUp.setOption(calcPopMeasuresPars.getPairwiseComps, "--getPairwiseComps", "get Pairwise comparison metrics");
	setUp.setOption(noDiagAlnPairwiseComps, "--noDiagAlnPairwiseComps", "Use diagonal Alignment for Pairwise Comparisons");
	calcPopMeasuresPars.diagAlnPairwiseComps = !noDiagAlnPairwiseComps;

	setUp.setOption(variantCallerRunPars.occurrenceCutOff, "--occurrenceCutOff", "Occurrence Cut Off, don't report variants below this count");
	setUp.setOption(variantCallerRunPars.lowVariantCutOff, "--lowVariantCutOff", "Low Variant Cut Off, don't report variants below this fraction");
	calcPopMeasuresPars.lowVarFreq = variantCallerRunPars.lowVariantCutOff;

	setUp.setOption(outwardsExpand, "--outwardsExpand", "The amount to expand outwards from given region when determining variants positions with extracted ref seq");
	setUp.setOption(metaFnp,    "--metaFnp",    "Meta data to add to sequences");

	bool directRefSeqInput = setUp.processSeq(false);
	if(!directRefSeqInput){
		setUp.setOption(genomeFnp, "--genome", "Genome file to extract ref seq from", true);
		setUp.setOption(bedFnp,    "--bed",    "A bed file of the location for the extraction", true);
		gPars.genomeDir_ = njh::files::normalize(genomeFnp.parent_path());
		gPars.primaryGenome_ = bfs::basename(genomeFnp);
		gPars.primaryGenome_ = gPars.primaryGenome_.substr(0, gPars.primaryGenome_.rfind("."));
		gPars.selectedGenomes_ = {gPars.primaryGenome_};
	}

	setUp.setOption(ignoreSubFields, "--ignoreSubFields", "Meta Sub Field values to ignore when calculating variants, e.g. --ignoreSubFields \"isFieldSample:TRUE,PreferredSample:FALSE\"");

	setUp.setOption(identifier, "--identifier", "Give a identifier name for info", false);
	setUp.setOption(keepNonFieldSamples, "--keepNonFieldSamples", "Keep Non Field Samples for population stats");
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);


	//gPos, refRegion (always forward strand), refSeq

	seqInfo refSeq;
	std::shared_ptr<Bed6RecordCore> gPos;
	GenomicRegion refRegion;

	if(directRefSeqInput){
		refSeq = setUp.pars_.seqObj_.seqBase_;
		gPos = std::make_shared<Bed6RecordCore>(refSeq.name_, 0, len(refSeq), refSeq.name_, len(refSeq), '+');
		refRegion = GenomicRegion(*gPos);
	}else{
		njh::files::checkExistenceThrow(bedFnp, __PRETTY_FUNCTION__);
		njh::files::checkExistenceThrow(genomeFnp, __PRETTY_FUNCTION__);

		MultiGenomeMapper gMapper(gPars);
		gMapper.loadInGenomes();
		gMapper.setUpGenomes();

		auto beds = getBeds(bedFnp);
		if(beds.empty()){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error no records found in " << bedFnp << "\n";
			throw std::runtime_error{ss.str()};
		}

		gPos = beds.front();
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

		refRegion = GenomicRegion(*gPos);
		refRegion.reverseSrand_ = false;
		TwoBit::TwoBitFile tReader(gMapper.genomes_.at(gPars.primaryGenome_)->fnpTwoBit_);



		refSeq = refRegion.extractSeq(tReader);

		SeqOutput::write(std::vector<seqInfo> { refSeq },
				SeqIOOptions::genFastaOut(
						njh::files::make_path(setUp.pars_.directoryName_,
								"inputRegion.fasta")));
	}
	//get identifer for sequences
	if("" == identifier){
		identifier = gPos->name_;
	}
	//read in meta if available
	std::unique_ptr<MultipleGroupMetaData> meta;
	if("" != metaFnp){
		meta = std::make_unique<MultipleGroupMetaData>(metaFnp);
	}

	std::unordered_map<std::string, std::set<std::string>> metaValuesToAvoid = njh::progutils::CmdArgs::sepSubArgsMulti<std::string, std::string>(ignoreSubFields);


	auto inputSeqs = CollapsedHaps::readInReads(setUp.pars_.ioOptions_, meta, metaValuesToAvoid);

	inputSeqs.setFrequencies();
	if (gPos->reverseStrand()) {
		inputSeqs.revCompSeqs();
	}
	//samples names
	auto sampNamesPerSeq = inputSeqs.getSampleNamesPerSeqs();
	auto allSamples = inputSeqs.getAllSampleNames();
	//rename based on freq
	inputSeqs.renameBaseOnFreq(identifier);

	uint64_t maxLen = readVec::getMaxLength(inputSeqs.seqs_);
	readVec::getMaxLength(refSeq, maxLen);

	std::shared_ptr<aligner> alignerObj = std::make_shared<aligner>(maxLen, gapScoringParameters(5,1,0,0,0,0), substituteMatrix(2,-2), false);
	alignerObj->weighHomopolymers_ = false;
	alignerObj->processAlnInfoInput(setUp.pars_.alnInfoDirName_, setUp.pars_.verbose_);
	//set up variant info
	auto idSeq = refSeq;
	idSeq.name_ = identifier;
	TranslatorByAlignment::VariantsInfo varInfo(refRegion.genBed3RecordCore(), idSeq);
	//get variant info
	auto refComps = inputSeqs.getCompsAgainstRef(refSeq, *alignerObj, numThreads);
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
	{
		calcPopMeasuresPars.numSegSites_ = varInfo.getFinalNumberOfSegratingSites();
		auto divMeasures = inputSeqs.getGeneralMeasuresOfDiversity(
				calcPopMeasuresPars, alignerObj);
		divMeasures.writeDivMeasures(
				njh::files::make_path(setUp.pars_.directoryName_,
						"divMeasures.tab.txt"), inputSeqs, identifier, calcPopMeasuresPars);
	}

	//std::cout << __FILE__ << " " << __LINE__ << std::endl;
	std::unordered_map<std::string, uint32_t> seqNameKey = inputSeqs.genSeqNameKey();
	//write out seqs
	auto uniqueSeqsOpts = SeqIOOptions::genFastaOutGz(njh::files::make_path(setUp.pars_.directoryName_, "uniqueSeqs.fasta.gz"));
	inputSeqs.writeOutSeqsOrdCnt(uniqueSeqsOpts);
	inputSeqs.writeNames(OutOptions(njh::files::make_path(setUp.pars_.directoryName_, "uniqueSeqs_names.tab.txt.gz")));
	inputSeqs.writeOutMetaFields(OutOptions(njh::files::make_path(setUp.pars_.directoryName_, "uniqueSeqs_meta.tab.txt.gz")));
	inputSeqs.writeLabIsolateNames(OutOptions(njh::files::make_path(setUp.pars_.directoryName_, "uniqueSeqs_labIsolateNames.tab.txt.gz")));
	auto variantCallsDir = njh::files::make_path(setUp.pars_.directoryName_, "variantCalls");
	njh::files::makeDir(njh::files::MkdirPar{variantCallsDir});
	if(!varInfo.snpsFinal.empty() || ! varInfo.deletionsFinal.empty() || !varInfo.insertionsFinal.empty()){
		varInfo.writeVCF(njh::files::make_path(variantCallsDir, "variants.vcf"));
		GenomicRegion variableRegion(varInfo.getVariableRegion());
		OutputStream bedVariableRegionOut(OutOptions(njh::files::make_path(variantCallsDir, "variableRegion.bed")));
		bedVariableRegionOut << variableRegion.genBedRecordCore().toDelimStrWithExtra() << std::endl;
		std::string variableRegionID = njh::pasteAsStr(variableRegion.chrom_, "::", variableRegion.start_, "::", variableRegion.end_);
		OutputStream variantsPerSeqOut(OutOptions(njh::files::make_path(variantCallsDir, "variantsPerSeq.tab.txt")));
		auto order = inputSeqs.getOrderByTopCnt();
		for(const auto & pos : order){
			refComps[pos].comp_.distances_.writeBasicInfo(variantsPerSeqOut, varInfo.region_, inputSeqs.seqs_[pos]->name_);
		}
	}

	alignerObj->processAlnInfoOutput(setUp.pars_.outAlnInfoDirName_, setUp.pars_.verbose_);


	return 0;
}



} //namespace njhseq

