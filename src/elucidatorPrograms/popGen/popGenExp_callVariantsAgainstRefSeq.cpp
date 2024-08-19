/*
 * popGenExp_callVariantsAgainstRefSeq.cpp
 *
 *  Created on: Jun 23, 2021
 *      Author: nick
 */


#include "popGenExp.hpp"
#include "elucidator/objects/BioDataObject.h"
#include <njhseq/objects/counters/DNABaseCounter.hpp>
#include <njhseq/objects/seqContainers/CollapsedHaps.hpp>

#include <njhseq/concurrency/PairwisePairFactory.hpp>
#include <njhseq/GenomeUtils.h>
#include <njhseq/objects/Gene/GenomicAminoAcidPositionTyper.hpp>
#include <njhseq/objects/Gene/TranslatorByAlignment.hpp>


#include <njhseq/objects/seqObjects/Clusters/identicalCluster.hpp>


namespace njhseq {


int popGenExpRunner::callVariantsAgainstRefSeq(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path genomeFnp = "";
	bfs::path gffFnp = "";
	TranslatorByAlignment::RunPars variantCallerRunPars;
	TranslatorByAlignment::TranslatorByAlignmentPars transPars;
	CollapsedHaps::GenPopMeasuresPar calcPopMeasuresPars;


	std::string identifier;
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
	setUp.setOption(noDiagAlnPairwiseComps, "--noDiagAlnPairwiseComps", "Use diagonal Alignment for Pairwise Comparisons");
	calcPopMeasuresPars.diagAlnPairwiseComps = !noDiagAlnPairwiseComps;

	setUp.setOption(variantCallerRunPars.occurrenceCutOff, "--occurrenceCutOff", "Occurrence Cut Off, don't report variants below this count");
	setUp.setOption(variantCallerRunPars.lowVariantCutOff, "--lowVariantCutOff", "Low Variant Cut Off, don't report variants below this fraction");
	calcPopMeasuresPars.lowVarFreq = variantCallerRunPars.lowVariantCutOff;

	setUp.setOption(outwardsExpand, "--outwardsExpand", "The amount to expand outwards from given region when determining variants positions with extracted ref seq");
	setUp.setOption(metaFnp,    "--metaFnp",    "Meta data to add to sequences");

	bool directRefSeqInput = setUp.processSeq(false);
	if(!directRefSeqInput){
		setUp.setOption(gffFnp, "--gffFnp", "gff Fnp for translating", false);
		setUp.setOption(genomeFnp, "--genome", "Genome file to extract ref seq from", true);
		setUp.setOption(bedFnp,    "--bed",    "A bed file of the location for the extraction", true);
		gPars.genomeDir_ = njh::files::normalize(genomeFnp.parent_path());
		gPars.primaryGenome_ = bfs::basename(genomeFnp);
		gPars.primaryGenome_ = gPars.primaryGenome_.substr(0, gPars.primaryGenome_.rfind('.'));
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
		if(identifier.empty()){
			gPos = std::make_shared<Bed6RecordCore>(*beds.front());
		} else {
			bool found = false;
			for(const auto & region : beds){
				if(region->name_ == identifier){
					gPos = std::make_shared<Bed6RecordCore>(*region);
					found = true;
					break;
				}
			}
			if(!found){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "didn't find " << identifier << " in " << bedFnp << "\n";
				throw std::runtime_error{ss.str()};
			}
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

		refRegion = GenomicRegion(*gPos);
		refRegion.reverseSrand_ = false;
		TwoBit::TwoBitFile tReader(gMapper.genomes_.at(gPars.primaryGenome_)->fnpTwoBit_);



		refSeq = refRegion.extractSeq(tReader);

		SeqOutput::write(std::vector<seqInfo> { refSeq },
				SeqIOOptions::genFastaOut(
						njh::files::make_path(setUp.pars_.directoryName_,
								"inputRegion.fasta")));
	}
	//get identifier for sequences
	if(identifier.empty()){
		identifier = gPos->name_;
	}
	//read in meta if available
	std::unique_ptr<MultipleGroupMetaData> meta;
	if(!metaFnp.empty()){
		meta = std::make_unique<MultipleGroupMetaData>(metaFnp);
	}

	std::unordered_map<std::string, std::set<std::string>> metaValuesToAvoid = njh::progutils::CmdArgs::sepSubArgsMulti<std::string, std::string>(ignoreSubFields);


	auto inputSeqs = CollapsedHaps::readInReads(setUp.pars_.ioOptions_, meta, metaValuesToAvoid);

	inputSeqs.setFrequencies();
	if (gPos->reverseStrand()) {
		inputSeqs.revCompSeqs();
	}
	//samples names
	auto sampReadCountsPerSeq = inputSeqs.getSampleReadCntsPerSeqs();

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
		// seqInfo(refComps[pos].comp_.refName_, refComps[pos].refAlnSeq_).outPutSeqAnsi(std::cout);
		// seqInfo(refComps[pos].comp_.queryName_, refComps[pos].queryAlnSeq_).outPutSeqAnsi(std::cout);
		varInfo.addVariantInfo(
				refComps[pos].refAlnSeq_,
				refComps[pos].queryAlnSeq_,
				inputSeqs.seqs_[pos]->cnt_,
				sampReadCountsPerSeq[pos],
				refComps[pos].comp_,
				refRegion.start_);
	}
	varInfo.setFinals(variantCallerRunPars);


	////std::cout << __FILE__ << " " << __LINE__ << std::endl;
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
	//std::cout << __FILE__ << " " << __LINE__ << std::endl;
	{
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		calcPopMeasuresPars.numSegSites_ = varInfo.getFinalNumberOfSegratingSites();
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		auto divMeasures = inputSeqs.getGeneralMeasuresOfDiversity(
				calcPopMeasuresPars, alignerObj);
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		divMeasures.writeDivMeasures(
				njh::files::make_path(setUp.pars_.directoryName_,
						"divMeasures.tab.txt"), inputSeqs, identifier, calcPopMeasuresPars);
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
	}
	if(!gffFnp.empty()) {
		if (gPos->reverseStrand()) {
			inputSeqs.revCompSeqs();
		}
		MultiGenomeMapper gMapper(gPars);
		gMapper.loadInGenomes();
		gMapper.setUpGenomes();
		TwoBit::TwoBitFile tReader(gMapper.genomes_.at(gPars.primaryGenome_)->fnpTwoBit_);

		//get gene for position
		auto geneInfoDir = njh::files::make_path(setUp.pars_.directoryName_, "geneInfos");
		njh::files::makeDir(njh::files::MkdirPar{geneInfoDir});
		OutOptions outOpts(njh::files::make_path(geneInfoDir, "gene"));
		auto inputRegions = std::vector<GenomicRegion>{GenomicRegion(*gPos)};
		auto geneIds = getFeatureIdsFromOverlappingRegions(inputRegions,
		                                                   gffFnp);

		// auto geneIds = getFeatureIdsFromOverlappingRegions(inputRegions,
		//                                                    std::function<const GenomicRegion&(const GenomicRegion&)>(
		// 	                                                   [](const GenomicRegion& reg) { return getRef(reg); }),
		//                                                    gffFnp);


		if(!geneIds.empty()) {
			SeqOutput transwriter(SeqIOOptions::genFastaOutGz(njh::files::make_path(variantCallsDir, "translatedInput.fasta.gz")));
			std::unordered_map<std::string, VecStr> idToTranscriptName;
			std::unordered_map<std::string, std::shared_ptr<GeneFromGffs>> rawGenes = GeneFromGffs::getGenesFromGffForIds(gffFnp, geneIds);
	//		std::cout << "rawGenes.size(): " << rawGenes.size() << std::endl;

			std::unordered_map<std::string, std::shared_ptr<GeneFromGffs>> genes;
			std::unordered_map<std::string, std::unordered_map<std::string, std::shared_ptr<GeneSeqInfo>>>  transcriptInfosForGene;
			std::unordered_map<std::string, std::shared_ptr<GeneSeqInfo>>  translationInfoForTranscirpt;
			std::unordered_map<std::string, TranslatorByAlignment::VariantsInfo> proteinVariants;
			std::unordered_map<std::string, std::unordered_map<std::string, TranslatorByAlignment::TranslateSeqRes>> translations;
			std::unordered_map<std::string, std::unordered_map<std::string, TranslatorByAlignment::TranslateSeqRes>> filteredOffTranslations;

			for(const auto & gene : rawGenes){
				bool failFilter = false;
				for(const auto & transcript : gene.second->mRNAs_){
					if(njh::in(njh::strToLowerRet(transcript->type_), VecStr{"rrna", "trna", "snorna","snrna","ncrna"}) ){
						////std::cout << __FILE__ << " " << __LINE__ << std::endl;
						transcript->writeGffRecord(std::cout);
						failFilter = true;
						break;
					}
				}

				if(!failFilter){
					 ////std::cout << __FILE__ << " " << __LINE__ << std::endl;
					genes[gene.first] = gene.second;
				}
			}
	//		std::cout << "genes.size(): " << genes.size() << std::endl;

			//std::unordered_map<std::string, std::unordered_map<std::string, std::shared_ptr<GeneSeqInfo>>> geneTranscriptInfos;
			uint64_t proteinMaxLen = readVec::getMaxLength(inputSeqs.seqs_);
			std::unordered_map<std::string, std::unordered_map<std::string, std::shared_ptr<GeneFromGffs>>> genesByChrom;
			// ////std::cout << __FILE__ << " " << __LINE__ << std::endl;
			for(const auto & gene : genes){
				genesByChrom[gene.second->gene_->seqid_].emplace(gene.first, gene.second);
				for(const auto & transcript : gene.second->mRNAs_){
					idToTranscriptName[gene.second->gene_->getIDAttr()].emplace_back(transcript->getIDAttr());
				}
				transcriptInfosForGene[gene.first] = gene.second->generateGeneSeqInfo(tReader, false);
				for(const auto & transcriptInfo : transcriptInfosForGene[gene.first]){
					translationInfoForTranscirpt[transcriptInfo.first] = transcriptInfo.second;
					readVec::getMaxLength(transcriptInfo.second->protein_, proteinMaxLen);
		//			ret.proteinForTranscript_[transcriptInfo.first] = transcriptInfo.second->protein_.seq_;
					proteinVariants.emplace(transcriptInfo.first,
											TranslatorByAlignment::VariantsInfo{Bed3RecordCore(transcriptInfo.first, 0, transcriptInfo.second->protein_.seq_.size()), transcriptInfo.second->protein_});
				}
				gene.second->writeOutGeneInfo(tReader, outOpts);
			}
			// ////std::cout << __FILE__ << " " << __LINE__ << std::endl;
			aligner alignObjForProtein(proteinMaxLen, gapScoringParameters(6,1,0,0,0,0), substituteMatrix(2,-2));
			//aligner alignObjSeq(seqMaxLen + rPars.realnPars.extendAmount * 2, gapScoringParameters(5,1,0,0,0,0), substituteMatrix(2,-2));

			std::unordered_map<std::string, std::unordered_map<std::string, std::set<std::string>>> regionsToGeneIds;
			//targetName, GeneID, AA Position
			std::unordered_map<std::string, std::vector<std::string>> alnRegionToGeneIds;
			for(const auto & region : inputRegions){
				for (const auto & g : genesByChrom[region.chrom_]) {
					for(const auto & t : g.second->mRNAs_){
						regionsToGeneIds[region.createUidFromCoords()][g.first].emplace(t->getIDAttr());
						alnRegionToGeneIds[region.createUidFromCoords()].emplace_back(g.first);
					}
				}
			}

			std::set<std::string> unmappable;
			std::set<std::string> untranslatable;
			ReAlignedSeq::genRealignmentPars realingPars;
			auto chromLengths = tReader.getSeqLens();

			GenomicRegion refSeqRegion(*gPos);
			realingPars.extendAmount = 0;//set to zero since above we already expanded
			for(const auto pos :  iter::range(inputSeqs.size())) {
				// seqInfo(refComps[pos].comp_.refName_, refComps[pos].refAlnSeq_).outPutSeqAnsi(std::cout);
				// seqInfo(refComps[pos].comp_.queryName_, refComps[pos].queryAlnSeq_).outPutSeqAnsi(std::cout);
				// std::cout << std::endl;

				if(refComps[pos].refAlnSeq_.front() == '-' || refComps[pos].refAlnSeq_.back() == '-') {

					unmappable.emplace(inputSeqs.seqs_[pos]->name_);

				} else {
					for(const auto & g : geneIds) {
						const auto & currentGene = njh::mapAt(genes, g);
						////std::cout << __FILE__ << " " << __LINE__ << std::endl;
						const auto & currentGeneInfo = njh::mapAt(transcriptInfosForGene, g);
						std::unordered_map<std::string, TranslatorByAlignment::TranslateSeqRes> currentTranslations;
						//            std::cout << __PRETTY_FUNCTION__  << " " << __LINE__ << std::endl;
						auto results = ReAlignedSeq::genRealignment(*inputSeqs.seqs_[pos], refSeqRegion, *alignerObj, chromLengths, tReader, realingPars);
						currentTranslations = TranslatorByAlignment::translateBasedOnAlignment(results, *currentGene, currentGeneInfo, alignObjForProtein, transPars);
						for(const auto & trans : currentTranslations){
							auto queryTransStart = trans.second.queryAlnTranslation_.seq_.find_first_not_of('-');
							if('-' == trans.second.refAlnTranslation_.seq_[queryTransStart] || countOccurences(trans.second.queryAlnTranslation_.seq_, "*") > transPars.allowableStopCodons_){
								//probably should do a more intensive check here fo
								filteredOffTranslations[inputSeqs.seqs_[pos]->name_].emplace(trans);
								//                std::cout << __PRETTY_FUNCTION__  << " " << __LINE__ << std::endl;
								untranslatable.emplace(inputSeqs.seqs_[pos]->name_);
							} else {
								translations[inputSeqs.seqs_[pos]->name_].emplace(trans);
								//                std::cout << __PRETTY_FUNCTION__  << " " << __LINE__ << std::endl;
							}
						}
					}
				}
			}
			OutputStream unmmapableSeqsOut(njh::files::make_path(variantCallsDir, "unmappable.txt"));
			unmmapableSeqsOut << njh::conToStr(unmappable, "\n") << std::endl;
			OutputStream untranslatableSeqsOut(njh::files::make_path(variantCallsDir, "untranslatable.txt"));
			untranslatableSeqsOut << njh::conToStr(untranslatable, "\n") << std::endl;

			//index amino acid changes per transcript
			for(const auto & seqName : translations){
				for(const auto & transcript : seqName.second){
					if(countOccurences(transcript.second.queryAlnTranslation_.seq_, "*") > transPars.allowableStopCodons_){
						//should log which ones have messed up translations
					}else{
						auto seqKey = inputSeqs.genSeqNameKey();
						auto sampNamesPerKey = inputSeqs.getSampleReadCntsPerSeqs();
						auto popCount = sampNamesPerKey.at(njh::mapAt(seqKey, seqName.first)).size();
						transwriter.openWrite(transcript.second.translation_);
						proteinVariants.at(transcript.first).addVariantInfo(
								transcript.second.refAlnTranslation_.seq_,
								transcript.second.queryAlnTranslation_.seq_,
								popCount,
								sampNamesPerKey.at(njh::mapAt(seqKey, seqName.first)),
								transcript.second.comp_,
								0);
					}
				}
			}
			////std::cout << __FILE__ << " " << __LINE__ << std::endl;
			for(auto & varPerTrans : proteinVariants){
				varPerTrans.second.setFinals(variantCallerRunPars);
				varPerTrans.second.writeVCF(njh::files::make_path(variantCallsDir, varPerTrans.first + ".vcf"));
			}
		}
	}



	//std::cout << __FILE__ << " " << __LINE__ << std::endl;
	alignerObj->processAlnInfoOutput(setUp.pars_.outAlnInfoDirName_, setUp.pars_.verbose_);


	return 0;
}



} //namespace njhseq

