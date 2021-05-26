/*
 * callVariantsAgainstRefSeq.cpp
 *
 *  Created on: Feb 24, 2021
 *      Author: nick
 *
 *
 */





#include "seqUtilsInfoRunner.hpp"
#include "elucidator/seqToolsUtils/seqToolsUtils.hpp"
#include "elucidator/objects/BioDataObject.h"
#include <njhseq/GenomeUtils.h>
#include <njhseq/objects/Gene/GeneFromGffs.hpp>

#include "elucidator/BioRecordsUtils/BedUtility.hpp"
#include "elucidator/PopulationGenetics.h"
#include <njhseq/objects/Gene/GenomicAminoAcidPositionTyper.hpp>

#include <njhseq/objects/Gene/TranslatorByAlignment.hpp>


#include <njhseq.h>



namespace njhseq {




int seqUtilsInfoRunner::callVariantsAgainstRefSeqIndividual(const njh::progutils::CmdArgs & inputCommands) {




	OutOptions outOpts(bfs::path(""), ".tab.txt");

	bfs::path bedFnp = "";
	uint32_t outwardsExpand = 5;

	bfs::path genomeFnp = "";



	MultiGenomeMapper::inputParameters gPars;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processReadInNames(true);
	setUp.processWritingOptions(outOpts);
	setUp.processAlnInfoInput();
	setUp.setOption(outwardsExpand, "--outwardsExpand", "The amount to expand outwards from given region when determining variants positions with extracted ref seq");

	setUp.setOption(bedFnp,    "--bed",    "A bed file of the location for the extraction", true);
	setUp.setOption(genomeFnp, "--genome", "A reference genome to compare against", true);
	if(!bfs::is_regular_file(genomeFnp)){
		setUp.failed_ = true;
		setUp.addWarning(njh::pasteAsStr(genomeFnp, " should be a file, not a directory"));
	}
	gPars.genomeDir_ = njh::files::normalize(genomeFnp.parent_path());
	gPars.primaryGenome_ = bfs::basename(genomeFnp);
	gPars.primaryGenome_ = gPars.primaryGenome_.substr(0, gPars.primaryGenome_.rfind("."));
	gPars.selectedGenomes_ = {gPars.primaryGenome_};

	setUp.finishSetUp(std::cout);


	njh::files::checkExistenceThrow(bedFnp, __PRETTY_FUNCTION__);
	njh::files::checkExistenceThrow(genomeFnp, __PRETTY_FUNCTION__);
	OutputStream out(outOpts);
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

	uint64_t maxLen = 0;
	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();

	std::vector<identicalCluster> originalOrientationClusters;
	seqInfo seq;



	// read in reads and collapse to unique
	while(reader.readNextRead(seq)) {
		//get meta keys if available
		readVec::handelLowerCaseBases(seq, setUp.pars_.ioOptions_.lowerCaseBases_);
		if(setUp.pars_.ioOptions_.removeGaps_){
			seq.removeGaps();
		}
		bool found = false;
		for (auto &cIter : originalOrientationClusters) {
			if (cIter.seqBase_.seq_ == seq.seq_) {
				cIter.addRead(seq);
				found = true;
				break;
			}
		}
		if (!found) {
			originalOrientationClusters.emplace_back(seq);
		}
	}
	std::vector<identicalCluster> forwardStrandClusters;
	if (gPos->reverseStrand()) {
		for (const auto &clus : originalOrientationClusters) {
			forwardStrandClusters.emplace_back(clus);
			forwardStrandClusters.back().seqBase_.reverseComplementRead(false, true);
		}
	} else {
		forwardStrandClusters = originalOrientationClusters;
	}

	njh::sort(forwardStrandClusters);


	uint32_t oldLen = gPos->length();
	BedUtility::extendLeftRight(*gPos, outwardsExpand, outwardsExpand,
			gMapper.genomes_.at(gPars.primaryGenome_)->chromosomeLengths_.at(
					gPos->chrom_));
	if(oldLen == gPos->score_){
		gPos->score_ = gPos->length();
	}

	GenomicRegion refRegion(*gPos);

	TwoBit::TwoBitFile tReader(
			gMapper.genomes_.at(gPars.primaryGenome_)->fnpTwoBit_);
	auto regionOrientationReferenceSeq = refRegion.extractSeq(tReader);
	auto forwardStrandRefSeq = regionOrientationReferenceSeq;
	if(refRegion.reverseSrand_){
		forwardStrandRefSeq = regionOrientationReferenceSeq	;
		forwardStrandRefSeq.reverseComplementRead(false, true);
	}

	readVec::getMaxLength(forwardStrandRefSeq, maxLen);
	aligner alignerObj(maxLen, gapScoringParameters(5,1,0,0,0,0), substituteMatrix(2,-2), false);
	alignerObj.weighHomopolymers_ = false;
	alignerObj.processAlnInfoInput(setUp.pars_.alnInfoDirName_, setUp.pars_.verbose_);

	out << "#chrom\tstart\tend\tname\ttype\tref\tquery" << "\n";
	for(const auto & seq : forwardStrandClusters){

		alignerObj.alignCacheGlobal(forwardStrandRefSeq, seq);
		alignerObj.profileAlignment(forwardStrandRefSeq, seq, false, false, false);
		if(setUp.pars_.debug_){
			alignerObj.alignObjectA_.seqBase_.outPutSeq(std::cout);
			alignerObj.alignObjectB_.seqBase_.outPutSeq(std::cout);

		}
		//iterate over the sub clusters
		for(const auto & subName : seq.reads_){
			//mismatches
			for(const auto & m : alignerObj.comp_.distances_.mismatches_){
				out << njh::conToStr(toVecStr(
						gPos->chrom_,
						gPos->chromStart_ + m.second.refBasePos,
						gPos->chromStart_ + m.second.refBasePos + 1,
						subName->seqBase_.name_,
						"SNP",
						m.second.refBase,
						m.second.seqBase), "\t") << "\n";
			}
			//INDELs
			for(const auto & indel: alignerObj.comp_.distances_.alignmentGaps_){
				if(indel.second.ref_){
					//gap is in reference so a insertion
					out << njh::conToStr(toVecStr(
							gPos->chrom_,
							gPos->chromStart_ + indel.second.refPos_,
							gPos->chromStart_ + indel.second.refPos_ + 1,
							subName->seqBase_.name_,
							"insertion",
							std::string(indel.second.gapedSequence_.size(), '-'),
							indel.second.gapedSequence_), "\t") << "\n";
				}else{
					//gap is in query so a deletion
					out << njh::conToStr(toVecStr(
							gPos->chrom_,
							gPos->chromStart_ + indel.second.refPos_,
							gPos->chromStart_ + indel.second.refPos_ + indel.second.gapedSequence_.size(),
							subName->seqBase_.name_,
							"deletion",
							indel.second.gapedSequence_,
							std::string(indel.second.gapedSequence_.size(), '-')), "\t") << "\n";
				}
			}
		}
	}
	return 0;
}


int seqUtilsInfoRunner::callVariantsAgainstRefSeq(const njh::progutils::CmdArgs & inputCommands) {

	/**@todo should add the following 1) average pairwise difference, 2) make into a function, 3) generate connected hap map, 4) unique haps to region count, 5) doing multiple pop fields at once */

  TranslatorByAlignment::TranslatorByAlignmentPars transPars;
  TranslatorByAlignment::RunPars variantCallerRunPars;
  bfs::path knownAminoAcidChangesFnp;

	std::string identifier = "";
	variantCallerRunPars.lowVariantCutOff = 0.005;
	variantCallerRunPars.occurrenceCutOff = 1;

	bfs::path bedFnp = "";
	uint32_t outwardsExpand = 5;

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
	setUp.setOption(variantCallerRunPars.occurrenceCutOff, "--occurrenceCutOff", "Occurrence Cut Off, don't report variants below this count");
	setUp.setOption(variantCallerRunPars.lowVariantCutOff, "--lowVariantCutOff", "Low Variant Cut Off, don't report variants below this fraction");
	setUp.setOption(outwardsExpand, "--outwardsExpand", "The amount to expand outwards from given region when determining variants positions with extracted ref seq");
	setUp.setOption(metaFnp,    "--metaFnp",    "Meta data to add to sequences");

	setUp.setOption(bedFnp,    "--bed",    "A bed file of the location for the extraction", true);
	setUp.setOption(transPars.lzPars_.genomeFnp, "--genome", "A reference genome to compare against", true);
	setUp.setOption(transPars.gffFnp_,    "--gff", "A gff3 file for genome file");
	if(!bfs::is_regular_file(transPars.lzPars_.genomeFnp)){
		setUp.failed_ = true;
		setUp.addWarning(njh::pasteAsStr(transPars.lzPars_.genomeFnp, " should be a file, not a directory"));
	}
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
	if("" != transPars.gffFnp_){
		translator = std::make_unique<TranslatorByAlignment>(transPars);
		if("" == transPars.lzPars_.genomeFnp){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << " if supplying gff file, must also supply --genomeFnp"<< "\n";
			throw std::runtime_error{ss.str()};
		}
	}
	table knownAminoAcidChanges;
	if("" != knownAminoAcidChangesFnp){
		if("" == transPars.gffFnp_){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "if supplying known amino acid positions than must also supply --gffFnp file"<< "\n";
			throw std::runtime_error{ss.str()};
		}
		knownAminoAcidChanges = table(knownAminoAcidChangesFnp, "\t", true);
		njh::for_each(knownAminoAcidChanges.columnNames_, [](std::string & col){
			njh::strToLower(col);
		});
		knownAminoAcidChanges.setColNamePositions();
		knownAminoAcidChanges.checkForColumnsThrow(VecStr{"transcriptid", "aaposition"}, __PRETTY_FUNCTION__);
	}
	std::unordered_map<std::string, std::set<uint32_t>> knownMutationsLocationsMap;
	if(knownAminoAcidChanges.nRow() > 0){
		for(const auto & row : knownAminoAcidChanges){
			if(std::all_of(row.begin(), row.end(), [](const std::string & element){
				return "" ==element;
			})){
				continue;
			}
			knownMutationsLocationsMap[row[knownAminoAcidChanges.getColPos("transcriptid")]].emplace(njh::StrToNumConverter::stoToNum<uint32_t>(row[knownAminoAcidChanges.getColPos("aaposition")]));
		}
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
	uint64_t maxLen = 0;
	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();

	std::unique_ptr<MultipleGroupMetaData> meta;
	if("" != metaFnp){
		meta = std::make_unique<MultipleGroupMetaData>(metaFnp);
	}

	std::vector<identicalCluster> originalOrientationClusters;
	seqInfo seq;
	uint32_t totalInputSeqs = 0;
	std::set<std::string> samples;

	std::unordered_map<std::string, std::string> metaValuesToAvoid;
	for(const auto & fieldValue : ignoreSubFields){
		auto toks = tokenizeString(fieldValue, ":");
		if(2 != toks.size()){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "meta field values to ignore need to be in format of METAFIELD:VALUE not " << fieldValue<< "\n";
			throw std::runtime_error{ss.str()};
		}
		metaValuesToAvoid[toks[0]] = toks[1];
	}


	std::unordered_map<std::string, std::unordered_set<std::string>> samplesPerSeq;
	// read in reads and collapse to unique
	while(reader.readNextRead(seq)) {
		if ("" != metaFnp) {
			meta->attemptToAddSeqMeta(seq);
		}
		//get meta keys if available
		if(MetaDataInName::nameHasMetaData(getSeqBase(seq).name_)){
			MetaDataInName metaData(getSeqBase(seq).name_);
			bool skip = false;
			for(const auto & ignoreField : metaValuesToAvoid){
				if(metaData.containsMeta(ignoreField.first) && metaData.getMeta(ignoreField.first) == ignoreField.second){
					skip = true;
					break;
					//skip this seq
				}
			}
			if(skip){
				continue;
			}
		}

		readVec::handelLowerCaseBases(seq, setUp.pars_.ioOptions_.lowerCaseBases_);
		if(setUp.pars_.ioOptions_.removeGaps_){
			seq.removeGaps();
		}
		readVec::getMaxLength(seq, maxLen);
		std::string sample = seq.name_;

		if(seq.nameHasMetaData()){
			MetaDataInName seqMeta(seq.name_);
			if(seqMeta.containsMeta("sample")){
				sample = seqMeta.getMeta("sample");
			}else if(seqMeta.containsMeta("BiologicalSample")){
				sample = seqMeta.getMeta("BiologicalSample");
			}
		}
		samples.emplace(sample);
		samplesPerSeq[seq.seq_].emplace(sample);

		++totalInputSeqs;
		bool found = false;
		for (auto &cIter : originalOrientationClusters) {
			if (cIter.seqBase_.seq_ == seq.seq_) {
				cIter.addRead(seq);
				found = true;
				break;
			}
		}
		if (!found) {
			originalOrientationClusters.emplace_back(seq);
		}
	}
	std::vector<identicalCluster> forwardStrandClusters;


	uint32_t samplesCalled = totalInputSeqs;
	if(!samples.empty()){
		samplesCalled = samples.size();
	}

	if(gPos->reverseStrand()){
		for(const auto & clus : originalOrientationClusters){
			forwardStrandClusters.emplace_back(clus);
			forwardStrandClusters.back().seqBase_.reverseComplementRead(false, true);
		}
	}else{
		forwardStrandClusters = originalOrientationClusters;
	}

	njh::sort(forwardStrandClusters);
	uint32_t seqId = 0;
	std::unordered_map<std::string, std::string> nameLookUp;
	std::unordered_map<std::string, uint32_t> sampCountsForPopHaps;
	std::unordered_map<std::string, std::unordered_set<std::string>> sampNamesForPopHaps;

	for (auto &cIter : forwardStrandClusters) {
		MetaDataInName popMeta;
		popMeta.addMeta("HapPopUIDCount", static_cast<uint32_t>(std::round(cIter.seqBase_.cnt_)));
		cIter.seqBase_.name_ = njh::pasteAsStr(identifier, ".", leftPadNumStr<uint32_t>(seqId, forwardStrandClusters.size()),popMeta.createMetaName());
		sampCountsForPopHaps[cIter.seqBase_.name_] = cIter.seqBase_.cnt_;
		sampNamesForPopHaps[cIter.seqBase_.name_] = samplesPerSeq[cIter.seqBase_.seq_];
		nameLookUp[cIter.seqBase_.seq_] = cIter.seqBase_.name_;
		++seqId;
	}


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

	TwoBit::TwoBitFile tReader(
			gMapper.genomes_.at(gPars.primaryGenome_)->fnpTwoBit_);
	auto regionOrientationReferenceSeq = refRegion.extractSeq(tReader);
	auto forwardStrandRefSeq = regionOrientationReferenceSeq;
	std::unordered_map<uint32_t, char> baseForPosition;
	if(refRegion.reverseSrand_){
		forwardStrandRefSeq = regionOrientationReferenceSeq	;
		forwardStrandRefSeq.reverseComplementRead(false, true);
	}

	for(uint32_t seqPos = 0; seqPos < forwardStrandRefSeq.seq_.size(); ++ seqPos){
		baseForPosition[refRegion.start_ + seqPos] = forwardStrandRefSeq.seq_[seqPos];
	}

	SeqOutput::write(std::vector<seqInfo> { regionOrientationReferenceSeq },
			SeqIOOptions::genFastaOut(
					njh::files::make_path(setUp.pars_.directoryName_,
							"inputRegion.fasta")));

	readVec::getMaxLength(forwardStrandRefSeq, maxLen);
	aligner alignerObj(maxLen, gapScoringParameters(5,1,0,0,0,0), substituteMatrix(2,-2), false);
	alignerObj.weighHomopolymers_ = false;
	alignerObj.processAlnInfoInput(setUp.pars_.alnInfoDirName_, setUp.pars_.verbose_);

	auto idSeq = forwardStrandRefSeq;
	idSeq.name_ = identifier;
	TranslatorByAlignment::VariantsInfo varInfo(refRegion.genBed3RecordCore(), idSeq);

	for(const auto & seq : forwardStrandClusters){

		alignerObj.alignCacheGlobal(forwardStrandRefSeq, seq);
		alignerObj.profileAlignment(forwardStrandRefSeq, seq, false, false, false);
		varInfo.addVariantInfo(
				alignerObj.alignObjectA_.seqBase_.seq_,
				alignerObj.alignObjectB_.seqBase_.seq_,
				seq.seqBase_.cnt_,
				samplesPerSeq[seq.seqBase_.seq_],
				alignerObj.comp_,
				refRegion.start_);
	}
	varInfo.setFinals(variantCallerRunPars);


	std::string variableRegionID = "";
	if(!varInfo.variablePositons_.empty()){
		uint32_t variableStart = vectorMinimum(std::vector<uint32_t>(varInfo.variablePositons_.begin(), varInfo.variablePositons_.end()));
		uint32_t variableStop = vectorMaximum(std::vector<uint32_t>(varInfo.variablePositons_.begin(), varInfo.variablePositons_.end()));

		size_t genomeVarStart = variableStart;
		size_t genomeVarStop = variableStop;

		GenomicRegion variableRegion = refRegion;
		variableRegion.start_ = genomeVarStart;
		variableRegion.end_ = genomeVarStop + 1;
		OutputStream bedVariableRegionOut(OutOptions(njh::files::make_path(setUp.pars_.directoryName_, "variableRegion.bed")));
		bedVariableRegionOut << variableRegion.genBedRecordCore().toDelimStrWithExtra() << std::endl;
		variableRegionID = njh::pasteAsStr(variableRegion.chrom_, "::", variableRegion.start_, "::", variableRegion.end_);
	}



	{
//		uint32_t numVariants = varInfo.snpsFinal.size() + varInfo.insertionsFinal.size() + varInfo.deletionsFinal.size();
//		OutputStream out(njh::files::make_path(setUp.pars_.directoryName_, "basicInfo.tab.txt"));
//		out << "id\ttotalHaplotypes\tuniqueHaplotypes\tnsamples\tsingletons\tdoublets\texpShannonEntropy\tShannonEntropyE\teffectiveNumOfAlleles\the\tlengthPolymorphism\tnvariantsAbove"
//				<< variantCallerRunPars.lowVariantCutOff * 100 << "%" << "\tvariableRegionID\t" << "avgGCContent";
		//out << "\tAvgPairwiseDistWeighted";
//		out << std::endl;
		std::map<uint32_t, uint32_t> readLens;
		//writing out unique sequences
		auto uniqueSeqsOpts = SeqIOOptions::genFastaOut(njh::files::make_path(setUp.pars_.directoryName_, "uniqueSeqs.fasta"));
		OutputStream nameOut(OutOptions(njh::files::make_path(setUp.pars_.directoryName_, "uniqueSeqs_names.tab.txt")));
		nameOut << "name\tnumber\treads"	<< std::endl;
		{
			SeqOutput uniqueWriter(uniqueSeqsOpts);
			uniqueWriter.openOut();
			uniqueWriter.write(originalOrientationClusters);
			uniqueWriter.closeOut();
		}


		OutputStream metaLabNamesOut(OutOptions(njh::files::make_path(setUp.pars_.directoryName_, "uniqueSeqs_nonFieldSampleNames.tab.txt")));
		metaLabNamesOut << "name\tsamples" << std::endl;
		for(const auto & cIter : forwardStrandClusters){
			auto inputNames = readVec::getNames(cIter.reads_);
			nameOut << cIter.seqBase_.name_
					<< "\t" << cIter.seqBase_.cnt_
					<< "\t" << njh::conToStr(inputNames, ",") << std::endl;
			readLens[len(cIter.seqBase_)]+= cIter.seqBase_.cnt_;
			VecStr nonFieldSampleNames{};
			for(const auto & name : inputNames){
				//std::cout << name << std::endl;

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
			metaLabNamesOut << cIter.seqBase_.name_ << "\t" << njh::conToStr(nonFieldSampleNames, ",") << std::endl;
		}

		std::map<std::string, std::map<std::string, MetaDataInName>> knownAAMeta;

		std::map<std::string, MetaDataInName> fullAATyped;
		std::map<std::string, std::vector<TranslatorByAlignment::AAInfo>> fullAATypedWithCodonInfo;
		if("" != transPars.gffFnp_){
			auto uniqueSeqInOpts = SeqIOOptions::genFastaIn(uniqueSeqsOpts.out_.outName());
			auto variantInfoDir =  njh::files::make_path(setUp.pars_.directoryName_, "proteinVariantInfo");

			njh::files::makeDir(njh::files::MkdirPar{variantInfoDir});
			translator->pars_.keepTemporaryFiles_ = true;
			translator->pars_.workingDirtory_ = variantInfoDir;
			auto translatedRes = translator->run(uniqueSeqInOpts, sampNamesForPopHaps, variantCallerRunPars);
			SeqOutput transwriter(SeqIOOptions::genFastaOut(njh::files::make_path(variantInfoDir, "translatedInput.fasta")));
			for(const auto & seqName : translatedRes.translations_){
				for(const auto & transcript : seqName.second){
					transwriter.openWrite(transcript.second.translation_);
				}
			}
			SeqInput popReader(uniqueSeqInOpts);
			auto popSeqs = popReader.readAllReads<seqInfo>();
			std::unordered_map<std::string, uint32_t> popSeqsPosition;
			for(const auto popPos : iter::range(popSeqs.size())){
				popSeqsPosition[popSeqs[popPos].name_] = popPos;
			}
			OutputStream popBedLocs(njh::files::make_path(variantInfoDir, "uniqueSeqs.bed"));
			for(const auto & seqLocs : translatedRes.seqAlns_){
				for(const auto & loc : seqLocs.second){
					popBedLocs << loc.gRegion_.genBedRecordCore().toDelimStrWithExtra() << std::endl;
				}
			}
			for(const auto & pop : popSeqs){
				if(!njh::in(pop.name_, translatedRes.seqAlns_)){
					popBedLocs << "*"
							<< "\t" << "*"
							<< "\t" << "*"
							<< "\t" << pop.name_
							<< "\t" << "*"
							<< "\t" << "*" << std::endl;
				}
			}
			std::unordered_map<std::string, std::set<uint32_t>> knownAAMutsChromPositions;
			{
				//protein
				for(auto & varPerTrans : translatedRes.proteinVariants_){
					auto snpsPositions = getVectorOfMapKeys(varPerTrans.second.snpsFinal);
					njh::sort(snpsPositions);
					OutputStream snpTabOut(OutOptions(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerTrans.first +  "-protein_aminoAcidVariable.tab.txt"))));
					snpTabOut << "transcript\tposition(1-based)\trefAA\tAA\tcount\tfraction\talleleDepth\tsamples" << std::endl;
					for(const auto & snpPos : snpsPositions){
						for(const auto & aa : varPerTrans.second.allBases[snpPos]){
							snpTabOut << varPerTrans.first
									<< "\t" << snpPos + 1
									<< "\t" << translatedRes.proteinForTranscript_[varPerTrans.first][snpPos]
									<< "\t" << aa.first
									<< "\t" << aa.second
									<< "\t" << aa.second/static_cast<double>(varPerTrans.second.depthPerPosition[snpPos])
									<< "\t" << varPerTrans.second.depthPerPosition[snpPos]
									<< "\t" << samplesCalled << std::endl;
						}
					}
					std::set<uint32_t> knownMutationsLocations;
					OutputStream allAATabOut(OutOptions(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerTrans.first +  "-protein_aminoAcidsAll.tab.txt"))));
					allAATabOut << "transcript\tposition(1-based)\trefAA\tAA\tcount\tfraction\talleleDepth\tsamples" << std::endl;
					for(const auto & snpPos : varPerTrans.second.allBases){
						if(njh::in(snpPos.first + 1, knownMutationsLocationsMap[varPerTrans.first])){
							knownMutationsLocations.emplace(snpPos.first);
							auto genomicLocationForAAPos = translatedRes.translationInfoForTranscirpt_.at(varPerTrans.first)->genBedFromAAPositions(snpPos.first, snpPos.first + 1);
							for(const auto gPos : iter::range(genomicLocationForAAPos.chromStart_, genomicLocationForAAPos.chromEnd_)){
								knownAAMutsChromPositions[genomicLocationForAAPos.chrom_].emplace(gPos);
							}
						}
						for(const auto & aa : snpPos.second){
							allAATabOut << varPerTrans.first
									<< "\t" << snpPos.first + 1
									<< "\t" << translatedRes.proteinForTranscript_[varPerTrans.first][snpPos.first]
									<< "\t" << aa.first
									<< "\t" << aa.second
									<< "\t" << aa.second/static_cast<double>(varPerTrans.second.depthPerPosition[snpPos.first])
									<< "\t" << varPerTrans.second.depthPerPosition[snpPos.first]
									<< "\t" << samplesCalled << std::endl;
						}
					}
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

					for (auto & seqName : translatedRes.translations_) {
						if (njh::in(varPerTrans.first, seqName.second)) {
							VecStr allAAPosCoded;
							std::string popName = seqName.first.substr(0, seqName.first.rfind("_f"));
							std::string transcript = varPerTrans.first;
							for (const auto & loc : allLocations) {
								//location is not within the aligned translation
								if(loc < std::get<0>(seqName.second[varPerTrans.first].firstAminoInfo_).aaPos_ || loc > std::get<0>(seqName.second[varPerTrans.first].lastAminoInfo_).aaPos_){
									continue;
								}
								auto aa = seqName.second[varPerTrans.first].queryAlnTranslation_.seq_[getAlnPosForRealPos(seqName.second[varPerTrans.first].refAlnTranslation_.seq_,loc)];
								fullAATyped[popName].addMeta(njh::pasteAsStr(varPerTrans.first, "-", loc + 1), aa);
//								std::cout << "loc: " << loc << std::endl;
								uint32_t refCDnaPos = loc * 3;
//								std::cout << "cDnaPos: " << cDnaPos << std::endl;
//								std::cout << "std::get<0>(seqName.second[varPerTrans.first].firstAminoInfo_).cDNAPos_: " << std::get<0>(seqName.second[varPerTrans.first].firstAminoInfo_).cDNAPos_ << std::endl;
								//subtract off the start location
								refCDnaPos -= std::get<0>(seqName.second[varPerTrans.first].firstAminoInfo_).cDNAPos_;
//								seqName.second[varPerTrans.first].cDna_.outPutSeqAnsi(std::cout);
//								seqName.second[varPerTrans.first].refAlnTranslation_.outPutSeq(std::cout);
//								seqName.second[varPerTrans.first].queryAlnTranslation_.outPutSeq(std::cout);
//								std::cout << "std::get<0>(seqName.second[varPerTrans.first].firstAminoInfo_).cDNAPos_;: " << std::get<0>(seqName.second[varPerTrans.first].firstAminoInfo_).cDNAPos_ << std::endl;
//								std::cout << "std::get<0>(seqName.second[varPerTrans.first].lastAminoInfo_).cDNAPos_;: " << std::get<0>(seqName.second[varPerTrans.first].lastAminoInfo_).cDNAPos_ << std::endl;
//								std::cout << "std::get<0>(seqName.second[varPerTrans.first].firstAminoInfo_).aaPos_;: " << std::get<0>(seqName.second[varPerTrans.first].firstAminoInfo_).aaPos_ << std::endl;
//								std::cout << "std::get<0>(seqName.second[varPerTrans.first].lastAminoInfo_).aaPos_;: " << std::get<0>(seqName.second[varPerTrans.first].lastAminoInfo_).aaPos_ << std::endl;
								uint32_t alnPosForLoc = getAlnPosForRealPos(seqName.second[varPerTrans.first].refAlnTranslation_.seq_,loc);
								uint32_t queryCDNAPos = getRealPosForAlnPos(seqName.second[varPerTrans.first].queryAlnTranslation_.seq_,alnPosForLoc ) * 3;

								if(queryCDNAPos +2 >=seqName.second[varPerTrans.first].cDna_.seq_.size()){
									std::stringstream ss;
									ss << __PRETTY_FUNCTION__ << ", error " << "" << queryCDNAPos +2 << " is greater than seqName.second[varPerTrans.first].cDna_.seq_.size(): " << seqName.second[varPerTrans.first].cDna_.seq_.size()<< "\n";
									throw std::runtime_error{ss.str()};
								}
								fullAATypedWithCodonInfo[popName].emplace_back(
										TranslatorByAlignment::AAInfo(varPerTrans.first, aa, loc,
												std::tie(
														seqName.second[varPerTrans.first].cDna_.seq_[queryCDNAPos+ 0],
														seqName.second[varPerTrans.first].cDna_.seq_[queryCDNAPos+ 1],
														seqName.second[varPerTrans.first].cDna_.seq_[queryCDNAPos+ 2]),
														njh::in(loc, knownMutationsLocations)));
							}
						}
					}
					if(!knownMutationsLocations.empty()){
						OutputStream knownTabOut(OutOptions(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerTrans.first +  "-protein_aminoAcidKnownMutations.tab.txt"))));
						knownTabOut << "transcript\tposition(1-based)\trefAA\tAA\tcount\tfraction\talleleDepth\tsamples" << std::endl;
						for(const auto & snpPos : knownMutationsLocations){
							for(const auto & aa : varPerTrans.second.allBases[snpPos]){
								knownTabOut << varPerTrans.first
										<< "\t" << snpPos + 1
										<< "\t" << translatedRes.proteinForTranscript_[varPerTrans.first][snpPos]
										<< "\t" << aa.first
										<< "\t" << aa.second
										<< "\t" << aa.second/static_cast<double>(varPerTrans.second.depthPerPosition[snpPos])
										<< "\t" << varPerTrans.second.depthPerPosition[snpPos]
										<< "\t" << samplesCalled << std::endl;
							}
						}
					}
				}
			}

			{
				//snps
				for( auto & varPerChrom : translatedRes.seqVariants_){
					auto snpsPositions = getVectorOfMapKeys(varPerChrom.second.snpsFinal);
					njh::sort(snpsPositions);
					OutputStream snpTabOut(OutOptions(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerChrom.first +  "-SNPs.tab.txt"))));
					snpTabOut << "chromosome\tposition(0-based)\trefBase\tbase\tcount\tfraction\talleleDepth\tsamples" << std::endl;
					for(const auto & snpPos : snpsPositions){
						for(const auto & base : varPerChrom.second.allBases[snpPos]){
							snpTabOut << varPerChrom.first
									<< "\t" << snpPos
									<< "\t" << translatedRes.baseForPosition_[varPerChrom.first][snpPos]
									<< "\t" << base.first
									<< "\t" << base.second
									<< "\t" << base.second/static_cast<double>(varPerChrom.second.depthPerPosition[snpPos])
									<< "\t" << varPerChrom.second.depthPerPosition[snpPos]
									<< "\t" << samplesCalled << std::endl;
						}
					}
					if(!knownAAMutsChromPositions[varPerChrom.first].empty()){
						OutputStream knownAASNPTabOut(OutOptions(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerChrom.first +  "-knownAA_SNPs.tab.txt"))));
						knownAASNPTabOut << "chromosome\tposition(0-based)\trefBase\tbase\tcount\tfraction\talleleDepth\tsamples" << std::endl;
						for(const auto & snpPos : knownAAMutsChromPositions[varPerChrom.first]){
							for(const auto & base : varPerChrom.second.allBases[snpPos]){
								knownAASNPTabOut << varPerChrom.first
										<< "\t" << snpPos
//<< "\t" <<
										<< "\t" << translatedRes.baseForPosition_[varPerChrom.first][snpPos]
										<< "\t" << base.first
										<< "\t" << base.second
										<< "\t" << base.second/static_cast<double>(varPerChrom.second.depthPerPosition[snpPos])
										<< "\t" << varPerChrom.second.depthPerPosition[snpPos]
										<< "\t" << samplesCalled << std::endl;
							}
						}
					}
					OutputStream allBasesTabOut(OutOptions(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerChrom.first +  "-allBases.tab.txt"))));
					allBasesTabOut << "chromosome\tposition(0-based)\trefBase\tbase\tcount\tfraction\talleleDepth\tsamples" << std::endl;
					for(const auto & snpPos : varPerChrom.second.allBases){
						for(const auto & base : snpPos.second){
							allBasesTabOut << varPerChrom.first
									<< "\t" << snpPos.first
									<< "\t" << translatedRes.baseForPosition_[varPerChrom.first][snpPos.first]
									<< "\t" << base.first
									<< "\t" << base.second
									<< "\t" << base.second/static_cast<double>(varPerChrom.second.depthPerPosition[snpPos.first])
									<< "\t" << varPerChrom.second.depthPerPosition[snpPos.first]
									<< "\t" << samplesCalled << std::endl;
						}
					}
					if(!varPerChrom.second.variablePositons_.empty()){
						GenomicRegion variableRegion = varPerChrom.second.getVariableRegion();
						OutputStream bedVariableRegionOut(OutOptions(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerChrom.first +  "-chromosome_variableRegion.bed"))));
						bedVariableRegionOut << variableRegion.genBedRecordCore().toDelimStrWithExtra() << std::endl;
					}
				}
			}
		}
	}

	auto variantCallsDir = njh::files::make_path(setUp.pars_.directoryName_, "variantCalls");
	njh::files::makeDir(njh::files::MkdirPar{variantCallsDir});


	alignerObj.processAlnInfoOutput(setUp.pars_.outAlnInfoDirName_, setUp.pars_.verbose_);
//	std::cout << "varInfo.snpsFinal.size(): " << varInfo.snpsFinal.size() << std::endl;
//	std::cout << "varInfo.deletionsFinal.size(): " << varInfo.deletionsFinal.size() << std::endl;
//	std::cout << "varInfo.insertionsFinal.size(): " << varInfo.insertionsFinal.size() << std::endl;
//

	if(!varInfo.snpsFinal.empty() || ! varInfo.deletionsFinal.empty() || !varInfo.insertionsFinal.empty()){
		varInfo.writeVCF(njh::files::make_path(variantCallsDir, "allVariants.vcf"));
		varInfo.writeSNPTable(njh::files::make_path(variantCallsDir, "allSNPs.tab.txt"));

	}
	return 0;
}





}  //namespace njhseq



