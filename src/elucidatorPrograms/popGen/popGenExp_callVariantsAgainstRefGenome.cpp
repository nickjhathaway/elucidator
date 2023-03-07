/*
 * popGenExp_callVariantsAgainstRefGenome.cpp
 *
 *  Created on: Sep 26, 2021
 *      Author: nicholashathaway
 */





#include "popGenExp.hpp"
#include <njhseq/objects/counters/DNABaseCounter.hpp>
#include <njhseq/PopulationGenetics.h>
#include <njhseq/PopulationGeneticsUtils.h>

#include <njhseq/objects/seqContainers/CollapsedHaps.hpp>
#include <njhseq/concurrency/PairwisePairFactory.hpp>
#include <njhseq/GenomeUtils.h>
#include <njhseq/objects/Gene/GenomicAminoAcidPositionTyper.hpp>
#include <njhseq/objects/Gene/TranslatorByAlignment.hpp>


#include <njhseq/objects/seqObjects/Clusters/identicalCluster.hpp>


namespace njhseq {


int popGenExpRunner::callVariantsAgainstRefGenome(const njh::progutils::CmdArgs & inputCommands) {

	/**@todo should add the following  2) make into a function, 3) generate connected hap map, 4) unique haps to region count, 5) doing multiple pop fields at once */


	bool noDiagAlnPairwiseComps = false;
	CollapseAndCallVariantsPars pars;

	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processReadInNames(true);
	pars.inOpts = setUp.pars_.ioOptions_;



	setUp.processAlnInfoInput();
	pars.alnCacheDir = setUp.pars_.alnInfoDirName_;

	setUp.setOption(pars.numThreads, "--numThreads", "number of threads");
	pars.calcPopMeasuresPars.numThreads = pars.numThreads;
	setUp.setOption(pars.calcPopMeasuresPars.getPairwiseComps, "--getPairwiseComps", "get Pairwise comparison metrics");
	setUp.setOption(pars.noDiagAlnPairwiseComps, "--noDiagAlnPairwiseComps", "Use diagonal Alignment for Pairwise Comparisons");
	pars.calcPopMeasuresPars.diagAlnPairwiseComps = !noDiagAlnPairwiseComps;

	setUp.setOption(pars.variantCallerRunPars.occurrenceCutOff, "--occurrenceCutOff", "Occurrence Cut Off, don't report variants below this count");
	setUp.setOption(pars.variantCallerRunPars.lowVariantCutOff, "--lowVariantCutOff", "Low Variant Cut Off, don't report variants below this fraction");
	pars.calcPopMeasuresPars.lowVarFreq = pars.variantCallerRunPars.lowVariantCutOff;
	pars.transPars.setOptions(setUp, true);

	setUp.setOption(pars.metaFnp,    "--metaFnp",    "Meta data to add to sequences");
	setUp.setOption(pars.ignoreSubFields, "--ignoreSubFields", "Meta Sub Field values to ignore when calculating variants, e.g. --ignoreSubFields \"isFieldSample:TRUE,PreferredSample:FALSE\"");

	setUp.setOption(pars.identifier, "--identifier", "Give a identifier name for info", true);
	setUp.setOption(pars.metaFieldsToCalcPopDiffs, "--metaFieldsToCalcPopDiffs", "Meta Fields To Calculate Pop Diffs");
	setUp.setOption(pars.outputDirectory, "--dout", "Output directory", true);
	setUp.setOption(pars.overWriteDirectory, "--overWriteDir", "over write output Directory");

	setUp.finishSetUp(std::cout);

	collapseAndCallVariants(pars);
	setUp.startARunLog(pars.outputDirectory.string());

//
//
//	//read in meta if available
//	std::unique_ptr<MultipleGroupMetaData> meta;
//	if("" != metaFnp){
//		meta = std::make_unique<MultipleGroupMetaData>(metaFnp);
//	}
//
//	std::unordered_map<std::string, std::set<std::string>> metaValuesToAvoid = njh::progutils::CmdArgs::sepSubArgsMulti<std::string, std::string>(ignoreSubFields);
//
//
//	auto inputSeqs = CollapsedHaps::readInReads(setUp.pars_.ioOptions_, meta, metaValuesToAvoid);
//	//samples names
//	auto sampNamesPerSeq = inputSeqs.getSampleNamesPerSeqs();
//	auto allSamples = inputSeqs.getAllSampleNames();
//	//rename based on freq
//	inputSeqs.renameBaseOnFreq(identifier);
//	//write out seqs
//	auto uniqueSeqsOpts = SeqIOOptions::genFastaOutGz(njh::files::make_path(setUp.pars_.directoryName_, "uniqueSeqs.fasta.gz"));
//	inputSeqs.writeOutAll(setUp.pars_.directoryName_, "uniqueSeqs");
//
//	uint64_t maxLen = readVec::getMaxLength(inputSeqs.seqs_);
//
//	std::shared_ptr<aligner> alignerObj = std::make_shared<aligner>(maxLen, gapScoringParameters(5,1,0,0,0,0), substituteMatrix(2,-2), false);
//	alignerObj->weighHomopolymers_ = false;
//	alignerObj->processAlnInfoInput(setUp.pars_.alnInfoDirName_, setUp.pars_.verbose_);
//
//	std::unordered_map<std::string, uint32_t> seqNameKey = inputSeqs.genSeqNameKey();
//
//
//	auto variantInfoDir =  njh::files::make_path(setUp.pars_.directoryName_, "variantCalls");
//	njh::files::makeDir(njh::files::MkdirPar{variantInfoDir});
//	std::unique_ptr<TranslatorByAlignment> translator = std::make_unique<TranslatorByAlignment>(transPars);
//	translator->pars_.keepTemporaryFiles_ = true;
//	translator->pars_.workingDirtory_ = variantInfoDir;
//	std::unordered_map<std::string, std::unordered_set<std::string>> sampNamesForPopHaps;
//	for(const auto pos : iter::range(inputSeqs.size())){
//		sampNamesForPopHaps[inputSeqs.seqs_[pos]->name_] = sampNamesPerSeq[pos];
//	}
//	auto translatedRes = translator->run(SeqIOOptions::genFastaIn(uniqueSeqsOpts.out_.outName()), sampNamesForPopHaps, variantCallerRunPars);
//	OutputStream popBedLocs(njh::files::make_path(variantInfoDir, "inputSeqs.bed"));
//	translatedRes.writeOutSeqAlnIndvVars(njh::files::make_path(variantInfoDir, "variantsPerSeqAln.tab.txt.gz"));
//	translatedRes.writeSeqLocations(popBedLocs);
//
//	std::unordered_map<std::string, std::set<uint32_t>> knownAAMutsChromPositions;
//
//	OutputStream seqsUnableToBeMappedOut(njh::files::make_path(variantInfoDir, "seqsUnableToBeMapped.txt"));
//	seqsUnableToBeMappedOut << njh::conToStr(translatedRes.seqsUnableToBeMapped_, "\n") << std::endl;
//	OutputStream seqsTranslationFilteredOut(njh::files::make_path(variantInfoDir, "seqsTranslationFiltered.txt"));
//	seqsTranslationFilteredOut << njh::conToStr(translatedRes.seqsTranslationFiltered_, "\n") << std::endl;
//
//
//	if(!translatedRes.translations_.empty()){
//		SeqOutput transwriter(SeqIOOptions::genFastaOut(njh::files::make_path(variantInfoDir, "translatedInput.fasta.gz")));
//		std::vector<seqInfo> translatedSeqs;
//		std::vector<std::unordered_set<std::string>> translatedSeqInputNames;
//		auto seqNames = njh::getVecOfMapKeys(translatedRes.translations_);
//		njh::sort(seqNames);
//		for(const auto & seqName : seqNames){
//			for(const auto & transcript : translatedRes.translations_.at(seqName)){
//				transwriter.openWrite(transcript.second.translation_);
//				translatedSeqs.emplace_back(transcript.second.translation_);
//				translatedSeqs.back().cnt_ = inputSeqs.names_[seqNameKey[seqName]].size();
//				translatedSeqInputNames.emplace_back(inputSeqs.names_[seqNameKey[seqName]]);
//			}
//		}
//		auto inputTranslatedSeq = CollapsedHaps::collapseReads(translatedSeqs, translatedSeqInputNames);
//		std::string identifierTranslated = njh::pasteAsStr(identifier, "-translated");
//		inputTranslatedSeq.renameBaseOnFreq(identifierTranslated);
//		//write out seqs
//		inputTranslatedSeq.writeOutAll(variantInfoDir, "uniqueTranslatedSeqs");
//		//get div measures
//		calcPopMeasuresPars.numSegSites_ = translatedRes.proteinVariants_.begin()->second.getFinalNumberOfSegratingSites();
//		auto divMeasures = inputTranslatedSeq.getGeneralMeasuresOfDiversity(calcPopMeasuresPars, alignerObj);
//		divMeasures.writeDivMeasures(
//				njh::files::make_path(variantInfoDir, "translatedDivMeasures.tab.txt"),
//				inputTranslatedSeq,
//				identifierTranslated, calcPopMeasuresPars);
//
//		OutputStream transBedLocs(njh::files::make_path(variantInfoDir, "translatedInput.bed"));
//		translatedRes.writeSeqLocationsTranslation(transBedLocs);
//		translatedRes.writeOutTranslatedIndvVars(njh::files::make_path(variantInfoDir, "variantsPerTranslatedSeq.tab.txt.gz"), translator->knownAminoAcidPositions_);
//		{
//			//protein
//			for(auto & varPerTrans : translatedRes.proteinVariants_){
//				varPerTrans.second.writeVCF(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerTrans.first +  "-protein.vcf")));
//				varPerTrans.second.writeOutSNPsFinalInfo(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerTrans.first +  "-protein_aminoAcidVariable.tab.txt.gz")), varPerTrans.first, true);
//				std::set<uint32_t> knownMutationsLocationsZeroBased;
//				for(const auto & snpPos : varPerTrans.second.allBases){
//					if(njh::in(snpPos.first + 1, translator->knownAminoAcidPositions_[varPerTrans.first])){
//						knownMutationsLocationsZeroBased.emplace(snpPos.first);
//						auto genomicLocationForAAPos = translatedRes.translationInfoForTranscirpt_.at(varPerTrans.first)->genBedFromAAPositions(snpPos.first, snpPos.first + 1);
//						for(const auto gPos : iter::range(genomicLocationForAAPos.chromStart_, genomicLocationForAAPos.chromEnd_)){
//							knownAAMutsChromPositions[genomicLocationForAAPos.chrom_].emplace(gPos);
//						}
//					}
//				}
//				varPerTrans.second.writeOutSNPsAllInfo(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerTrans.first +  "-protein_aminoAcidsAll.tab.txt.gz")), varPerTrans.first, true);
//				if(!varPerTrans.second.variablePositons_.empty()){
//					GenomicRegion variableRegion = varPerTrans.second.getVariableRegion();
//					variableRegion.start_ += 1; //do one based positioning
//					OutputStream bedVariableRegionOut(OutOptions(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerTrans.first +  "-protein_variableRegion.bed"))));
//					bedVariableRegionOut << variableRegion.genBedRecordCore().toDelimStrWithExtra() << std::endl;
//				}
//				if(!knownMutationsLocationsZeroBased.empty()){
//					varPerTrans.second.writeOutSNPsInfo(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerTrans.first +  "-protein_aminoAcidKnownMutations.tab.txt.gz")), varPerTrans.first, knownMutationsLocationsZeroBased, true);
//				}
//			}
//
//			translatedRes.writeOutAATypedInfo(njh::files::make_path(variantInfoDir, "seqKnownAATyped.tab.txt.gz"));
//		}
//	}
//
//
//	//snps
//	uint32_t maxSeqCount = 0;
//
//	for(auto & varPerChrom : translatedRes.seqVariants_){
//		for(const auto & count : varPerChrom.second.depthPerPosition){
//			if(count.second > maxSeqCount){
//				//cheap way of doing this for now
//				maxSeqCount = count.second;
//				calcPopMeasuresPars.numSegSites_ = varPerChrom.second.getFinalNumberOfSegratingSites();
//			}
//		}
//		varPerChrom.second.writeVCF(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerChrom.first +  "-genomic.vcf")));
//		varPerChrom.second.writeOutSNPsFinalInfo(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerChrom.first +  "-SNPs.tab.txt.gz")), varPerChrom.first);
//		if(!knownAAMutsChromPositions[varPerChrom.first].empty()){
//			varPerChrom.second.writeOutSNPsInfo(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerChrom.first +  "-knownAA_SNPs.tab.txt.gz")), varPerChrom.first, knownAAMutsChromPositions[varPerChrom.first]);
//		}
//		varPerChrom.second.writeOutSNPsAllInfo(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerChrom.first +  "-allBases.tab.txt.gz")), varPerChrom.first);
//		if(!varPerChrom.second.variablePositons_.empty()){
//			GenomicRegion variableRegion = varPerChrom.second.getVariableRegion();
//			OutputStream bedVariableRegionOut(OutOptions(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerChrom.first +  "-chromosome_variableRegion.bed"))));
//			bedVariableRegionOut << variableRegion.genBedRecordCore().toDelimStrWithExtra() << std::endl;
//		}
//	}
//
//	{
//		auto snpTyped = translatedRes.genSeqSNPTypedStr();
//		OutputStream snpTypedOut(njh::files::make_path(variantInfoDir, "seqSnpTyped.tab.txt.gz"));
//		snpTypedOut << "name\tsnpTyped" << std::endl;
//		for(const auto & seqPos : inputSeqs.getOrderByTopCnt()){
//			snpTypedOut << inputSeqs.seqs_[seqPos]->name_
//					<< "\t" << snpTyped[inputSeqs.seqs_[seqPos]->name_];
//			snpTypedOut << std::endl;
//		}
//	}
//
//	auto divMeasures = inputSeqs.getGeneralMeasuresOfDiversity(calcPopMeasuresPars, alignerObj);
//	divMeasures.writeDivMeasures(
//			njh::files::make_path(setUp.pars_.directoryName_,
//					"divMeasures.tab.txt"), inputSeqs, identifier, calcPopMeasuresPars);
//
//	alignerObj->processAlnInfoOutput(setUp.pars_.outAlnInfoDirName_, setUp.pars_.verbose_);


	return 0;
}



} // namespace njhseq

