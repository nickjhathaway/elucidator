/*
 * popGenExp_variantCalling.cpp
 *
 *  Created on: Jun 11, 2021
 *      Author: nick
 */



#include "popGenExp.hpp"
#include "elucidator/objects/BioDataObject.h"
#include <njhseq/objects/counters/DNABaseCounter.hpp>


#include "elucidator/seqToolsUtils/seqToolsUtils.hpp"
#include <njhseq/GenomeUtils.h>
#include <njhseq/objects/Gene/GeneFromGffs.hpp>

#include "elucidator/BioRecordsUtils/BedUtility.hpp"
#include <njhseq/objects/Gene/GenomicAminoAcidPositionTyper.hpp>
#include <njhseq/objects/Gene/TranslatorByAlignment.hpp>

#include <njhseq/objects/seqObjects/Clusters/identicalCluster.hpp>
#include <njhseq/concurrency/pools/BamReaderPool.hpp>
#include <njhseq/concurrency/pools/AlignerPool.hpp>
#include <njhseq/PopulationGenetics.h>
#include <njhseq/objects/seqContainers/CollapsedHaps.hpp>



namespace njhseq {



int popGenExpRunner::quickHaplotypeVariantsWithRegion(const njh::progutils::CmdArgs & inputCommands) {
	std::string identifier = "";
	std::string fstMeta = "";
	double lowVariantCutOff = 0.005;
	uint32_t occurrenceCutOff = 5;
	uint32_t lengthDiffForLengthPoly = 6;
	bfs::path bedFnp = "";
	bfs::path intersectBedFnp = "";

	bfs::path genomefnp = "";
	uint32_t outwardsExpand = 5;
	OutOptions outOpts(bfs::path(""));
	outOpts.outExtention_ = ".tab.txt";
	MultiGenomeMapper::inputParameters gPars;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processReadInNames(true);
	setUp.processWritingOptions(outOpts);
	setUp.processAlnInfoInput();
	setUp.setOption(lengthDiffForLengthPoly, "--lengthDiffForLengthPoly", "length Diff For Length Poly");
	setUp.setOption(occurrenceCutOff, "--occurrenceCutOff", "Occurrence Cut Off, don't report variants above this percentage");
	setUp.setOption(outwardsExpand, "--outwardsExpand", "The amount to expand outwards from given region when determining variants positions with extracted ref seq");
	setUp.setOption(intersectBedFnp,    "--intersectBed",    "A bed to get intersecting variants for", true);
	setUp.setOption(bedFnp,    "--bed",    "A bed file of the location for the extraction", true);
	setUp.setOption(genomefnp, "--genome", "A reference genome to compare against", true);
	gPars.genomeDir_ = genomefnp.parent_path();
	gPars.primaryGenome_ = bfs::basename(genomefnp);
	gPars.primaryGenome_ = gPars.primaryGenome_.substr(0, gPars.primaryGenome_.rfind("."));

	setUp.setOption(identifier, "--identifier", "Give a identifier name for info");
	setUp.setOption(fstMeta,    "--fstMeta",    "Meta field to calculate fsts by", true);
	setUp.finishSetUp(std::cout);


	njh::files::checkExistenceThrow(intersectBedFnp, __PRETTY_FUNCTION__);
	njh::files::checkExistenceThrow(bedFnp,          __PRETTY_FUNCTION__);
	njh::files::checkExistenceThrow(genomefnp,       __PRETTY_FUNCTION__);

	uint64_t maxLen = 0;
	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();

	std::unordered_map<std::string, std::shared_ptr<std::vector<identicalCluster>>> uniqueSeqsByMeta;
	std::vector<identicalCluster> clusters;
	seqInfo seq;
	// uint32_t totalInput = 0;
	bool calculatingFst = "" != fstMeta;
	std::set<std::string> samples;
	std::set<std::string> allMetaKeys;
	std::unordered_map<std::string, std::set<std::string>> samplesByMeta;
	// read in reads and collapse to unique
	while(reader.readNextRead(seq)) {
		//get meta keys if available
		if(MetaDataInName::nameHasMetaData(getSeqBase(seq).name_)){
			MetaDataInName metaData(getSeqBase(seq).name_);
			for(const auto & meta : metaData.meta_){
				allMetaKeys.emplace(meta.first);
			}
		}

		readVec::handelLowerCaseBases(seq, setUp.pars_.ioOptions_.lowerCaseBases_);
		if(setUp.pars_.ioOptions_.removeGaps_){
			seq.removeGaps();
		}
		readVec::getMaxLength(seq, maxLen);
		if(seq.nameHasMetaData()){
			MetaDataInName seqMeta(seq.name_);
			if(seqMeta.containsMeta("sample")){
				samples.emplace(seqMeta.getMeta("sample"));
			}
		}
		// ++totalInput;
		bool found = false;
		for (auto &cIter : clusters) {
			if (cIter.seqBase_.seq_ == seq.seq_) {
				cIter.addRead(seq);
				found = true;
				break;
			}
		}
		if (!found) {
			clusters.emplace_back(seq);
		}
		//if calculating sub population differences collect unique sequences for the given field
		if(calculatingFst){
			MetaDataInName seqMeta(seq.name_);
			bool found = false;
			if(!njh::in(seqMeta.getMeta(fstMeta),  uniqueSeqsByMeta)){
				uniqueSeqsByMeta[seqMeta.getMeta(fstMeta)] = std::make_shared<std::vector<identicalCluster>>();
			}
			if(seqMeta.containsMeta("sample")){
				samplesByMeta[seqMeta.getMeta(fstMeta)].emplace(seqMeta.getMeta("sample"));
			}

			for (auto &cIter : *uniqueSeqsByMeta[seqMeta.getMeta(fstMeta)]) {
				if (cIter.seqBase_.seq_ == seq.seq_) {
					cIter.addRead(seq);
					found = true;
					break;
				}
			}
			if (!found) {
				uniqueSeqsByMeta[seqMeta.getMeta(fstMeta)]->emplace_back(seq);
			}
		}
	}

	VecStr allMetaKeysVec(allMetaKeys.begin(), allMetaKeys.end());

	// uint32_t samplesCalled = totalInput;
	// if(!samples.empty()){
	// 	samplesCalled = samples.size();
	// }
	auto beds = getBeds(bedFnp);
	if(beds.empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error no records found in " << bedFnp << "\n";
		throw std::runtime_error{ss.str()};
	}
	auto intersectingBeds = getBeds(intersectBedFnp);
	if(intersectingBeds.empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error no records found in " << intersectBedFnp << "\n";
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
	njh::sort(clusters);
	uint32_t seqId = 0;
	std::unordered_map<std::string, std::string> nameLookUp;
	for (auto &cIter : clusters) {
		MetaDataInName popMeta;
		popMeta.addMeta("HapPopUIDCount", static_cast<uint32_t>(std::round(cIter.seqBase_.cnt_)));
		cIter.seqBase_.name_ = njh::pasteAsStr(identifier, ".", leftPadNumStr<uint32_t>(seqId, clusters.size()),popMeta.createMetaName());
		nameLookUp[cIter.seqBase_.seq_] = cIter.seqBase_.name_;
		++seqId;
	}
	if(calculatingFst){
		std::unordered_map<std::string, uint32_t> seqIdForField;
		for(auto & field : uniqueSeqsByMeta){
			readVec::allSetFractionByTotalCount(*field.second);
			for (auto &cIter : *field.second) {
				MetaDataInName popMeta;
				std::string popName = nameLookUp[cIter.seqBase_.seq_];
				MetaDataInName popNameMeta(popName);
				MetaDataInName::removeMetaDataInName(popName);
				popMeta.addMeta("HapPopUID", popName);
				popMeta.addMeta(popNameMeta, true);
				popMeta.addMeta("HapUniqueCount", static_cast<uint32_t>(std::round(cIter.seqBase_.cnt_)));
				cIter.seqBase_.name_ = njh::pasteAsStr(field.first, ".", leftPadNumStr<uint32_t>(seqIdForField[field.first], field.second->size()), popMeta.createMetaName());
				++seqIdForField[field.first];
			}
		}
	}


	//get region and ref seq for mapping of variants
	BedUtility::extendLeftRight(*gPos, outwardsExpand, outwardsExpand,
			gMapper.genomes_.at(gPars.primaryGenome_)->chromosomeLengths_.at(
					gPos->chrom_));
	GenomicRegion refRegion(*gPos);

	TwoBit::TwoBitFile tReader(
			gMapper.genomes_.at(gPars.primaryGenome_)->fnpTwoBit_);
	auto refSeq = refRegion.extractSeq(tReader);

	readVec::getMaxLength(refSeq, maxLen);
	aligner alignerObj(maxLen, gapScoringParameters(5,1,0,0,0,0), substituteMatrix(2,-2), false);
	alignerObj.weighHomopolymers_ =false;
	alignerObj.processAlnInfoInput(setUp.pars_.alnInfoDirName_, setUp.pars_.verbose_);

	table outputTab(VecStr{"CHROM", "POS", "ID", "REF", "ALT", fstMeta, "count", "total"});
	if(calculatingFst){
		for(auto & field : uniqueSeqsByMeta){
			njh::sort(*field.second);
			double seqCount = readVec::getTotalReadCount(*field.second);

			std::unordered_map<uint32_t, std::unordered_map<char,        uint32_t>> snps;
			std::unordered_map<uint32_t, std::unordered_map<std::string, uint32_t>> insertions;
			std::unordered_map<uint32_t, std::unordered_map<std::string, uint32_t>> deletions;

			std::unordered_map<uint32_t, std::unordered_map<char,        uint32_t>> snpsFinal;
			std::unordered_map<uint32_t, std::unordered_map<std::string, uint32_t>> insertionsFinal;
			std::unordered_map<uint32_t, std::unordered_map<std::string, uint32_t>> deletionsFinal;

			for(const auto & seq : *field.second){
				alignerObj.alignCacheGlobal(refSeq, seq);
				alignerObj.profileAlignment(refSeq, seq, false, false, false);
				for(const auto & m : alignerObj.comp_.distances_.mismatches_){
					snps[m.second.refBasePos][m.second.seqBase]+= seq.seqBase_.cnt_;
				}
				for(const auto & g : alignerObj.comp_.distances_.alignmentGaps_){
					if(g.second.ref_){
						//insertion
						insertions[g.second.refPos_ - 1][g.second.gapedSequence_]++;
					}else{
						//deletion
						deletions[g.second.refPos_ - 1][g.second.gapedSequence_]++;
					}
				}
			}
			std::vector<uint32_t> variablePositions;
			//filter snps and indels by occurrence cut off
			for(const auto & snp : snps){
				for(const auto & b : snp.second){
					if(b.second < occurrenceCutOff){
						continue;
					}
					snpsFinal[snp.first][b.first] = b.second;
					if(b.second/static_cast<double>(seqCount) > lowVariantCutOff){
						variablePositions.emplace_back(snp.first);
					}
				}
			}
			for(const auto & del : deletions){
				for(const auto & d : del.second){
					if(d.second < occurrenceCutOff){
						continue;
					}
					deletionsFinal[del.first][d.first] = d.second;
					if(d.second/static_cast<double>(seqCount) > lowVariantCutOff){
						variablePositions.emplace_back(del.first);
					}
				}
			}
			for(const auto & ins : insertions){
				for(const auto & i : ins.second){
					if(i.second < occurrenceCutOff){
						continue;
					}
					insertionsFinal[ins.first][i.first] = i.second;
					if(i.second/static_cast<double>(seqCount) > lowVariantCutOff){
						variablePositions.emplace_back(ins.first);
					}
				}
			}
			if(!snpsFinal.empty() || ! deletionsFinal.empty() || !insertionsFinal.empty()){

				std::unordered_set<uint32_t> positionsSet;
				for(const auto & snps : snpsFinal){
					positionsSet.emplace(snps.first);
				}
				for(const auto & ins : insertionsFinal){
					positionsSet.emplace(ins.first);
				}
				for(const auto & del : deletionsFinal){
					positionsSet.emplace(del.first);
				}
				std::vector<uint32_t> positions(positionsSet.begin(), positionsSet.end());
				// uint32_t samplesCalledForField = seqCount;
				// if(!samplesByMeta[field.first].empty()){
				// 	samplesCalledForField = samplesByMeta[field.first].size();
				// }
				njh::sort(positions);
				if(refRegion.reverseSrand_){
					njh::reverse(positions);
				}
				for(const auto & pos : positions){
					auto relativePos = refRegion.getRelativePositionFromStartStrandAware(pos);
					if (intersectingBeds.front()->chrom_ == refRegion.chrom_
							&& relativePos < intersectingBeds.front()->chromEnd_
							&& relativePos >= intersectingBeds.front()->chromStart_) {
						//	table outputTab(VecStr{"CHROM", "POS", "ID", "REF", "ALT", fstMeta, "count", "total"});
						//intersectingBeds.front()->chrom_,relativePos,identifier,
						std::string ref = njh::pasteAsStr(refSeq.seq_[pos]);
						if(refRegion.reverseSrand_){
							ref = seqUtil::reverseComplement(std::string(1, refSeq.seq_[pos]), "DNA");
						}
						if(njh::in(pos, snpsFinal)){
							for(const auto & b : snpsFinal[pos]){
								outputTab.addRow(intersectingBeds.front()->chrom_,relativePos,identifier, ref, b.first, field.first, b.second, seqCount);
							}
						}
					}
				}
			}
		}
	}
	if(outputTab.nRow() > 0){
		OutputStream out(outOpts);
		outputTab.outPutContents(out, "\t");
	}
	return 0;
}


int popGenExpRunner::getHapPopDifAndVariantsInfo(const njh::progutils::CmdArgs & inputCommands) {
	return 0;
//
//	/**@todo should add the following 1) average pairwise difference, 2) make into a function, 3) generate connected hap map, 4) unique haps to region count, 5) doing multiple pop fields at once */
//
//  TranslatorByAlignment::TranslatorByAlignmentPars transPars;
//  TranslatorByAlignment::RunPars variantCallerRunPars;
//  bfs::path knownAminoAcidChangesFnp;
//
//	std::string identifier = "";
//	std::string popMeta = "";
//	variantCallerRunPars.lowVariantCutOff = 0.005;
//	variantCallerRunPars.occurrenceCutOff = 1;
//
//	uint32_t lengthDiffForLengthPoly = 6;
//	bfs::path bedFnp = "";
//	uint32_t outwardsExpand = 5;
//	bool zeroBased = false;
//	bool keepNonFieldSamples = false;
//	std::set<std::string> ignoreSubFields;
//	MultiGenomeMapper::inputParameters gPars;
//	seqSetUp setUp(inputCommands);
//	setUp.processVerbose();
//	setUp.processDebug();
//	setUp.processReadInNames(true);
//	setUp.processDirectoryOutputName(true);
//	setUp.processAlnInfoInput();
//	setUp.setOption(lengthDiffForLengthPoly, "--lengthDiffForLengthPoly", "If there is a length difference of this or above, label as there being length polymorphism");
//	setUp.setOption(variantCallerRunPars.occurrenceCutOff, "--occurrenceCutOff", "Occurrence Cut Off, don't report variants below this count");
//	setUp.setOption(variantCallerRunPars.lowVariantCutOff, "--lowVariantCutOff", "Low Variant Cut Off, don't report variants below this fraction");
//	setUp.setOption(outwardsExpand, "--outwardsExpand", "The amount to expand outwards from given region when determining variants positions with extracted ref seq");
//
//	setUp.setOption(bedFnp,    "--bed",    "A bed file of the location for the extraction", true);
//	setUp.setOption(transPars.lzPars_.genomeFnp, "--genome", "A reference genome to compare against", true);
//	setUp.setOption(transPars.gffFnp_,    "--gff", "A gff3 file for genome file");
//	if(!bfs::is_regular_file(transPars.lzPars_.genomeFnp)){
//		setUp.failed_ = true;
//		setUp.addWarning(njh::pasteAsStr(transPars.lzPars_.genomeFnp, " should be a file, not a directory"));
//	}
//	gPars.genomeDir_ = njh::files::normalize(transPars.lzPars_.genomeFnp.parent_path());
//	gPars.primaryGenome_ = bfs::basename(transPars.lzPars_.genomeFnp);
//	gPars.primaryGenome_ = gPars.primaryGenome_.substr(0, gPars.primaryGenome_.rfind("."));
//	gPars.selectedGenomes_ = {gPars.primaryGenome_};
//
//	setUp.setOption(ignoreSubFields, "--ignoreSubFields", "Meta Sub Fields to ignore when doing population analysis");
//
//	setUp.setOption(identifier, "--identifier", "Give a identifier name for info");
//	setUp.setOption(popMeta,    "--popMeta",    "Meta field to calculate population stats by");
//	setUp.setOption(knownAminoAcidChangesFnp, "--proteinMutantTypingFnp", "Protein Mutant Typing Fnp, columns should be ID=gene id in gff, AAPosition=amino acid position", false);
//	setUp.setOption(zeroBased, "--zeroBased", "If the positions in the proteinMutantTypingFnp are zero based");
//	setUp.setOption(keepNonFieldSamples, "--keepNonFieldSamples", "Keep Non Field Samples for population stats");
//
//	setUp.finishSetUp(std::cout);
//	setUp.startARunLog(setUp.pars_.directoryName_);
//
//	njh::files::checkExistenceThrow(bedFnp, __PRETTY_FUNCTION__);
//	njh::files::checkExistenceThrow(transPars.lzPars_.genomeFnp, __PRETTY_FUNCTION__);
//
//	std::unique_ptr<TranslatorByAlignment> translator;
//	if("" != transPars.gffFnp_){
//		translator = std::make_unique<TranslatorByAlignment>(transPars);
//		if("" == transPars.lzPars_.genomeFnp){
//			std::stringstream ss;
//			ss << __PRETTY_FUNCTION__ << ", error " << " if supplying gff file, must also supply --genomeFnp"<< "\n";
//			throw std::runtime_error{ss.str()};
//		}
//	}
//	table knownAminoAcidChanges;
//	if("" != knownAminoAcidChangesFnp){
//		if("" == transPars.gffFnp_){
//			std::stringstream ss;
//			ss << __PRETTY_FUNCTION__ << ", error " << "if supplying known amino acid positions than must also supply --gffFnp file"<< "\n";
//			throw std::runtime_error{ss.str()};
//		}
//		knownAminoAcidChanges = table(knownAminoAcidChangesFnp, "\t", true);
//		njh::for_each(knownAminoAcidChanges.columnNames_, [](std::string & col){
//			njh::strToLower(col);
//		});
//		knownAminoAcidChanges.setColNamePositions();
//		knownAminoAcidChanges.checkForColumnsThrow(VecStr{"transcriptid", "aaposition"}, __PRETTY_FUNCTION__);
//	}
//	std::unordered_map<std::string, std::set<uint32_t>> knownMutationsLocationsMap;
//	if(knownAminoAcidChanges.nRow() > 0){
//		for(const auto & row : knownAminoAcidChanges){
//			if(std::all_of(row.begin(), row.end(), [](const std::string & element){
//				return "" ==element;
//			})){
//				continue;
//			}
//			knownMutationsLocationsMap[row[knownAminoAcidChanges.getColPos("transcriptid")]].emplace(njh::StrToNumConverter::stoToNum<uint32_t>(row[knownAminoAcidChanges.getColPos("aaposition")]));
//		}
//	}
//	auto beds = getBeds(bedFnp);
//	if(beds.empty()){
//		std::stringstream ss;
//		ss << __PRETTY_FUNCTION__ << ", error no records found in " << bedFnp << "\n";
//		throw std::runtime_error{ss.str()};
//	}
//	MultiGenomeMapper gMapper(gPars);
//	gMapper.loadInGenomes();
//	gMapper.setUpGenomes();
//	auto gPos = beds.front();
//
//	//get identifer for sequences
//	if("" == identifier){
//		identifier = gPos->name_;
//	}
//	uint64_t maxLen = 0;
//	SeqInput reader(setUp.pars_.ioOptions_);
//	reader.openIn();
//
//	std::unordered_map<std::string, std::shared_ptr<std::vector<identicalCluster>>> uniqueSeqsByMetaForwardStrand;
//	std::unordered_map<std::string,std::unordered_map<std::string, std::unordered_set<std::string>>> samplesPerSeqPerMeta;
//
//	std::vector<identicalCluster> originalOrientationClusters;
//	seqInfo seq;
//	uint32_t totalInputSeqs = 0;
//	bool calculatingFst = "" != popMeta;
//	std::set<std::string> samples;
//	std::set<std::string> allMetaKeys;
//	std::unordered_map<std::string, std::set<std::string>> samplesByMeta;
//	std::unordered_map<std::string, std::unordered_set<std::string>> samplesPerSeq;
//
//	// read in reads and collapse to unique
//	while(reader.readNextRead(seq)) {
//		//get meta keys if available
//		if(MetaDataInName::nameHasMetaData(getSeqBase(seq).name_)){
//			MetaDataInName metaData(getSeqBase(seq).name_);
//			for(const auto & meta : metaData.meta_){
//				allMetaKeys.emplace(meta.first);
//			}
//		}
//
//		readVec::handelLowerCaseBases(seq, setUp.pars_.ioOptions_.lowerCaseBases_);
//		if(setUp.pars_.ioOptions_.removeGaps_){
//			seq.removeGaps();
//		}
//		readVec::getMaxLength(seq, maxLen);
//		std::string sample = seq.name_;
//
//		if(seq.nameHasMetaData()){
//			MetaDataInName seqMeta(seq.name_);
//			if(seqMeta.containsMeta("sample")){
//				sample = seqMeta.getMeta("sample");
//			}else if(seqMeta.containsMeta("BiologicalSample")){
//				sample = seqMeta.getMeta("BiologicalSample");
//			}
//		}
//		samples.emplace(sample);
//		samplesPerSeq[seq.seq_].emplace(sample);
//
//		++totalInputSeqs;
//		bool found = false;
//		for (auto &cIter : originalOrientationClusters) {
//			if (cIter.seqBase_.seq_ == seq.seq_) {
//				cIter.addRead(seq);
//				found = true;
//				break;
//			}
//		}
//		if (!found) {
//			originalOrientationClusters.emplace_back(seq);
//		}
//
//		//if calculating sub population differences collect unique sequences for the given field
//		if(calculatingFst){
//			MetaDataInName seqMeta(seq.name_);
//			bool found = false;
//			//skip ignore fields
//			if(njh::in(seqMeta.getMeta(popMeta), ignoreSubFields)){
//				continue;
//			}
//			//skip non field samples
//			if(!keepNonFieldSamples && seqMeta.containsMeta("IsFieldSample") && !seqMeta.getMeta<bool>("IsFieldSample")){
//				continue;
//			}
//			if(!njh::in(seqMeta.getMeta(popMeta),  uniqueSeqsByMetaForwardStrand)){
//				uniqueSeqsByMetaForwardStrand[seqMeta.getMeta(popMeta)] = std::make_shared<std::vector<identicalCluster>>();
//
//			}
//			if(seqMeta.containsMeta("sample")){
//				samplesByMeta[seqMeta.getMeta(popMeta)].emplace(seqMeta.getMeta("sample"));
//			}
//
//			for (auto &cIter : *uniqueSeqsByMetaForwardStrand[seqMeta.getMeta(popMeta)]) {
//				if (cIter.seqBase_.seq_ == seq.seq_) {
//					cIter.addRead(seq);
//					found = true;
//					break;
//				}
//			}
//			samplesPerSeqPerMeta[seqMeta.getMeta(popMeta)][seq.name_].emplace(sample);
//			if (!found) {
//				uniqueSeqsByMetaForwardStrand[seqMeta.getMeta(popMeta)]->emplace_back(seq);
//			}
//		}
//	}
//	std::vector<identicalCluster> forwardStrandClusters;
//
//
//
//	double sumGCContent = 0;
//	for(auto & clus : originalOrientationClusters){
//		clus.setLetterCount();
//		clus.counter_.calcGcContent();
//		sumGCContent += clus.counter_.gcContent_ * clus.seqBase_.cnt_;
//	}
//	double avgGCContent = sumGCContent/totalInputSeqs;
//
//
//	VecStr allMetaKeysVec(allMetaKeys.begin(), allMetaKeys.end());
//
//	uint32_t samplesCalled = totalInputSeqs;
//	if(!samples.empty()){
//		samplesCalled = samples.size();
//	}
//
//
//
//
//
//	if(gPos->reverseStrand()){
//		for( auto & field : uniqueSeqsByMetaForwardStrand){
//			for( auto & seq : *field.second){
//				seq.seqBase_.reverseComplementRead(false, true);
//			}
//		}
//		for(const auto & clus : originalOrientationClusters){
//			forwardStrandClusters.emplace_back(clus);
//			forwardStrandClusters.back().seqBase_.reverseComplementRead(false, true);
//		}
//	}else{
//		forwardStrandClusters = originalOrientationClusters;
//	}
//
//	njh::sort(forwardStrandClusters);
//	uint32_t seqId = 0;
//	std::unordered_map<std::string, std::string> nameLookUp;
//	std::unordered_map<std::string, uint32_t> sampCountsForPopHaps;
//	std::unordered_map<std::string, std::unordered_set<std::string>> sampNamesForPopHaps;
//
//	for (auto &cIter : forwardStrandClusters) {
//		MetaDataInName popMeta;
//		popMeta.addMeta("HapPopUIDCount", static_cast<uint32_t>(std::round(cIter.seqBase_.cnt_)));
//		cIter.seqBase_.name_ = njh::pasteAsStr(identifier, ".", leftPadNumStr<uint32_t>(seqId, forwardStrandClusters.size()),popMeta.createMetaName());
//		sampCountsForPopHaps[cIter.seqBase_.name_] = cIter.seqBase_.cnt_;
//		sampNamesForPopHaps[cIter.seqBase_.name_] = samplesPerSeq[cIter.seqBase_.seq_];
//
//		nameLookUp[cIter.seqBase_.seq_] = cIter.seqBase_.name_;
//		++seqId;
//	}
//
//	if(calculatingFst){
//		std::unordered_map<std::string, uint32_t> seqIdForField;
//		for(auto & field : uniqueSeqsByMetaForwardStrand){
//			readVec::allSetFractionByTotalCount(*field.second);
//			for (auto &cIter : *field.second) {
//				MetaDataInName popMeta;
//				std::string popName = nameLookUp[cIter.seqBase_.seq_];
//				MetaDataInName popNameMeta(popName);
//				MetaDataInName::removeMetaDataInName(popName);
//				popMeta.addMeta("HapPopUID", popName);
//				popMeta.addMeta(popNameMeta, true);
//				popMeta.addMeta("HapUniqueCount", static_cast<uint32_t>(std::round(cIter.seqBase_.cnt_)));
//				cIter.seqBase_.name_ = njh::pasteAsStr(field.first, ".", leftPadNumStr<uint32_t>(seqIdForField[field.first], field.second->size()), popMeta.createMetaName());
//				++seqIdForField[field.first];
//			}
//		}
//	}
//
//	//get region and ref seq for mapping of variants
//	OutputStream inputRegionBeforeExpandOut(OutOptions(njh::files::make_path(setUp.pars_.directoryName_, "inputRegionBeforeExpand.bed")));
//	inputRegionBeforeExpandOut << gPos->toDelimStrWithExtra() << std::endl;
//
//	uint32_t oldLen = gPos->length();
//	BedUtility::extendLeftRight(*gPos, outwardsExpand, outwardsExpand,
//			gMapper.genomes_.at(gPars.primaryGenome_)->chromosomeLengths_.at(
//					gPos->chrom_));
//	if(oldLen == gPos->score_){
//		gPos->score_ = gPos->length();
//	}
//
//	OutputStream bedOut(OutOptions(njh::files::make_path(setUp.pars_.directoryName_, "inputRegion.bed")));
//	bedOut << gPos->toDelimStrWithExtra() << std::endl;
//
//	GenomicRegion refRegion(*gPos);
//
//	TwoBit::TwoBitFile tReader(
//			gMapper.genomes_.at(gPars.primaryGenome_)->fnpTwoBit_);
//	auto regionOrientationReferenceSeq = refRegion.extractSeq(tReader);
//	auto forwardStrandRefSeq = regionOrientationReferenceSeq;
//	std::unordered_map<uint32_t, char> baseForPosition;
//	if(refRegion.reverseSrand_){
//		forwardStrandRefSeq = regionOrientationReferenceSeq	;
//		forwardStrandRefSeq.reverseComplementRead(false, true);
//	}
//
//	for(uint32_t seqPos = 0; seqPos < forwardStrandRefSeq.seq_.size(); ++ seqPos){
//		baseForPosition[refRegion.start_ + seqPos] = forwardStrandRefSeq.seq_[seqPos];
//	}
//
//	SeqOutput::write(std::vector<seqInfo> { regionOrientationReferenceSeq },
//			SeqIOOptions::genFastaOut(
//					njh::files::make_path(setUp.pars_.directoryName_,
//							"inputRegion.fasta")));
//
//	readVec::getMaxLength(forwardStrandRefSeq, maxLen);
//	aligner alignerObj(maxLen, gapScoringParameters(5,1,0,0,0,0), substituteMatrix(2,-2), false);
//	alignerObj.weighHomopolymers_ = false;
//	alignerObj.processAlnInfoInput(setUp.pars_.alnInfoDirName_, setUp.pars_.verbose_);
//	auto seqWithID = forwardStrandRefSeq;
//	seqWithID.name_ = identifier;
//	TranslatorByAlignment::VariantsInfo varInfo(refRegion.genBed3RecordCore(), seqWithID);
//
////	std::unordered_map<uint32_t, std::unordered_map<char, uint32_t>> snps;
////	std::unordered_map<uint32_t, std::unordered_map<std::string,uint32_t>> insertions;
////	std::unordered_map<uint32_t, std::unordered_map<std::string,uint32_t>> deletions;
////
////	std::unordered_map<uint32_t, std::unordered_map<char, uint32_t>> snpsFinal;
////	std::unordered_map<uint32_t, std::unordered_map<std::string,uint32_t>> insertionsFinal;
////	std::unordered_map<uint32_t, std::unordered_map<std::string,uint32_t>> deletionsFinal;
//	//
//
//	for(const auto & seq : forwardStrandClusters){
//
//		alignerObj.alignCacheGlobal(forwardStrandRefSeq, seq);
//		alignerObj.profileAlignment(forwardStrandRefSeq, seq, false, false, false);
//		varInfo.addVariantInfo(
//				alignerObj.alignObjectA_.seqBase_.seq_,
//				alignerObj.alignObjectB_.seqBase_.seq_,
//				seq.seqBase_.cnt_,
//				samplesPerSeq[seq.seqBase_.seq_],
//				alignerObj.comp_,
//				refRegion.start_);
//	}
//	varInfo.setFinals(variantCallerRunPars);
//
//
//	std::string variableRegionID = "";
//	if(!varInfo.variablePositons_.empty()){
//		uint32_t variableStart = vectorMinimum(std::vector<uint32_t>(varInfo.variablePositons_.begin(), varInfo.variablePositons_.end()));
//		uint32_t variableStop = vectorMaximum(std::vector<uint32_t>(varInfo.variablePositons_.begin(), varInfo.variablePositons_.end()));
//
//		size_t genomeVarStart = variableStart;
//		size_t genomeVarStop = variableStop;
//
//		GenomicRegion variableRegion = refRegion;
//		variableRegion.start_ = genomeVarStart;
//		variableRegion.end_ = genomeVarStop + 1;
//		OutputStream bedVariableRegionOut(OutOptions(njh::files::make_path(setUp.pars_.directoryName_, "variableRegion.bed")));
//		bedVariableRegionOut << variableRegion.genBedRecordCore().toDelimStrWithExtra() << std::endl;
//		variableRegionID = njh::pasteAsStr(variableRegion.chrom_, "::", variableRegion.start_, "::", variableRegion.end_);
//	}
//
//
//
//	{
//		uint32_t numVariants = varInfo.snpsFinal.size() + varInfo.insertionsFinal.size() + varInfo.deletionsFinal.size();
//		OutputStream out(njh::files::make_path(setUp.pars_.directoryName_, "basicInfo.tab.txt"));
//		out << "id\ttotalHaplotypes\tuniqueHaplotypes\tnsamples\tsingletons\tdoublets\texpShannonEntropy\tShannonEntropyE\teffectiveNumOfAlleles\the\tlengthPolymorphism\tnvariantsAbove"
//				<< variantCallerRunPars.lowVariantCutOff * 100 << "%" << "\tvariableRegionID\t" << "avgGCContent";
//		//out << "\tAvgPairwiseDistWeighted";
//		out << std::endl;
//		std::map<uint32_t, uint32_t> readLens;
//		//writing out unique sequences
//		auto uniqueSeqsOpts = SeqIOOptions::genFastaOut(njh::files::make_path(setUp.pars_.directoryName_, "uniqueSeqs.fasta"));
//		OutputStream nameOut(OutOptions(njh::files::make_path(setUp.pars_.directoryName_, "uniqueSeqs_names.tab.txt")));
//		nameOut << "name\tnumber\treads"	<< std::endl;
//		{
//			SeqOutput uniqueWriter(uniqueSeqsOpts);
//			uniqueWriter.openOut();
//			uniqueWriter.write(originalOrientationClusters);
//			uniqueWriter.closeOut();
//		}
//
//
//		OutputStream metaLabNamesOut(OutOptions(njh::files::make_path(setUp.pars_.directoryName_, "uniqueSeqs_nonFieldSampleNames.tab.txt")));
//		metaLabNamesOut << "name\tsamples" << std::endl;
//		for(const auto & cIter : forwardStrandClusters){
//			auto inputNames = readVec::getNames(cIter.reads_);
//			nameOut << cIter.seqBase_.name_
//					<< "\t" << cIter.seqBase_.cnt_
//					<< "\t" << njh::conToStr(inputNames, ",") << std::endl;
//			readLens[len(cIter.seqBase_)]+= cIter.seqBase_.cnt_;
//			VecStr nonFieldSampleNames{};
//			for(const auto & name : inputNames){
//				if(MetaDataInName::nameHasMetaData(name)){
//					MetaDataInName meta(name);
//					if(meta.containsMeta("BiologicalSample") && meta.containsMeta("IsFieldSample") && !meta.getMeta<bool>("IsFieldSample")){
//						if(!meta.containsMeta("country") || "NA" != meta.getMeta("country")){
//							nonFieldSampleNames.emplace_back(meta.getMeta("BiologicalSample"));
//						}
//					}
//				}
//			}
//			metaLabNamesOut << cIter.seqBase_.name_ << "\t" << njh::conToStr(nonFieldSampleNames, ",") << std::endl;
//		}
//
//
//		//       seqName               transcript   amino acid positions and amino acid
//		//std::map<std::string, std::map<std::string, MetaDataInName>> fullAATyped;
//		//seqname meta
//		std::map<std::string, std::map<std::string, MetaDataInName>> knownAAMeta;
//		std::map<std::string, MetaDataInName> fullAATyped;
//		std::map<std::string, std::vector<TranslatorByAlignment::AAInfo>> fullAATypedWithCodonInfo;
//		if("" != transPars.gffFnp_){
//			auto uniqueSeqInOpts = SeqIOOptions::genFastaIn(uniqueSeqsOpts.out_.outName());
//			auto variantInfoDir =  njh::files::make_path(setUp.pars_.directoryName_, "proteinVariantInfo");
//
//			njh::files::makeDir(njh::files::MkdirPar{variantInfoDir});
//			translator->pars_.keepTemporaryFiles_ = true;
//			translator->pars_.workingDirtory_ = variantInfoDir;
//			auto translatedRes = translator->run(uniqueSeqInOpts, sampNamesForPopHaps, variantCallerRunPars);
//			//write out the translation
//			SeqOutput transwriter(SeqIOOptions::genFastaOut(njh::files::make_path(variantInfoDir, "translatedInput.fasta")));
//			for(const auto & seqName : translatedRes.translations_){
//				for(const auto & transcript : seqName.second){
//					transwriter.openWrite(transcript.second.translation_);
//				}
//			}
//			SeqInput popReader(uniqueSeqInOpts);
//			auto popSeqs = popReader.readAllReads<seqInfo>();
//			OutputStream popBedLocs(njh::files::make_path(variantInfoDir, "uniqueSeqs.bed"));
//			for(const auto & seqLocs : translatedRes.seqAlns_){
//				for(const auto & loc : seqLocs.second){
//					popBedLocs << loc.gRegion_.genBedRecordCore().toDelimStrWithExtra() << std::endl;
//				}
//			}
//			for(const auto & pop : popSeqs){
//				if(!njh::in(pop.name_, translatedRes.seqAlns_)){
//					popBedLocs << "*"
//							<< "\t" << "*"
//							<< "\t" << "*"
//							<< "\t" << pop.name_
//							<< "\t" << "*"
//							<< "\t" << "*" << std::endl;
//				}
//			}
//			std::unordered_map<std::string, std::set<uint32_t>> knownAAMutsChromPositions;
//			{
//				//protein
//				for(auto & varPerTrans : translatedRes.proteinVariants_){
//					auto snpsPositions = getVectorOfMapKeys(varPerTrans.second.snpsFinal);
//					njh::sort(snpsPositions);
//					OutputStream snpTabOut(OutOptions(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerTrans.first +  "-protein_aminoAcidVariable.tab.txt"))));
//					snpTabOut << "transcript\tposition(1-based)\trefAA\tAA\tcount\tfraction\talleleDepth\tsamples" << std::endl;
//					for(const auto & snpPos : snpsPositions){
//						for(const auto & aa : varPerTrans.second.allBases[snpPos]){
//							snpTabOut << varPerTrans.first
//									<< "\t" << snpPos + 1
//									<< "\t" << translatedRes.proteinForTranscript_[varPerTrans.first][snpPos]
//									<< "\t" << aa.first
//									<< "\t" << aa.second
//									<< "\t" << aa.second/static_cast<double>(varPerTrans.second.depthPerPosition[snpPos])
//									<< "\t" << varPerTrans.second.depthPerPosition[snpPos]
//									<< "\t" << varPerTrans.second.samplesPerPosition[snpPos].size() << std::endl;
//						}
//					}
//					std::set<uint32_t> knownMutationsLocations;
//					OutputStream allAATabOut(OutOptions(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerTrans.first +  "-protein_aminoAcidsAll.tab.txt"))));
//					allAATabOut << "transcript\tposition(1-based)\trefAA\tAA\tcount\tfraction\talleleDepth\tsamples" << std::endl;
//					for(const auto & snpPos : varPerTrans.second.allBases){
//						if(njh::in(snpPos.first + 1, knownMutationsLocationsMap[varPerTrans.first])){
//							knownMutationsLocations.emplace(snpPos.first);
//							auto genomicLocationForAAPos = translatedRes.translationInfoForTranscirpt_.at(varPerTrans.first)->genBedFromAAPositions(snpPos.first, snpPos.first + 1);
//							for(const auto gPos : iter::range(genomicLocationForAAPos.chromStart_, genomicLocationForAAPos.chromEnd_)){
//								knownAAMutsChromPositions[genomicLocationForAAPos.chrom_].emplace(gPos);
//							}
//						}
//						for(const auto & aa : snpPos.second){
//							allAATabOut << varPerTrans.first
//									<< "\t" << snpPos.first + 1
//									<< "\t" << translatedRes.proteinForTranscript_[varPerTrans.first][snpPos.first]
//									<< "\t" << aa.first
//									<< "\t" << aa.second
//									<< "\t" << aa.second/static_cast<double>(varPerTrans.second.depthPerPosition[snpPos.first])
//									<< "\t" << varPerTrans.second.depthPerPosition[snpPos.first]
//									<< "\t" << samplesCalled << std::endl;
//						}
//					}
//					if(!varPerTrans.second.variablePositons_.empty()){
//						GenomicRegion variableRegion = varPerTrans.second.getVariableRegion();
//						variableRegion.start_ += 1; //do one based positioning
//						OutputStream bedVariableRegionOut(OutOptions(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerTrans.first +  "-protein_variableRegion.bed"))));
//						bedVariableRegionOut << variableRegion.genBedRecordCore().toDelimStrWithExtra() << std::endl;
//					}
//					std::set<uint32_t> allLocations(knownMutationsLocations.begin(), knownMutationsLocations.end());
//					for(const auto & variablePos : varPerTrans.second.snpsFinal){
//						allLocations.emplace(variablePos.first);
//					}
//
//					for (auto & seqName : translatedRes.translations_) {
//						if (njh::in(varPerTrans.first, seqName.second)) {
//							VecStr allAAPosCoded;
//							std::string popName = seqName.first.substr(0, seqName.first.rfind("_f"));
//							std::string transcript = varPerTrans.first;
//							for (const auto & loc : allLocations) {
//								//location is not within the aligned translation
//								if(loc < std::get<0>(seqName.second[varPerTrans.first].firstAminoInfo_).aaPos_ || loc > std::get<0>(seqName.second[varPerTrans.first].lastAminoInfo_).aaPos_){
//									continue;
//								}
//								auto aa = seqName.second[varPerTrans.first].queryAlnTranslation_.seq_[getAlnPosForRealPos(seqName.second[varPerTrans.first].refAlnTranslation_.seq_,loc)];
//								fullAATyped[popName].addMeta(njh::pasteAsStr(varPerTrans.first, "-", loc + 1), aa);
////								std::cout << "loc: " << loc << std::endl;
//								uint32_t refCDnaPos = loc * 3;
////								std::cout << "cDnaPos: " << cDnaPos << std::endl;
////								std::cout << "std::get<0>(seqName.second[varPerTrans.first].firstAminoInfo_).cDNAPos_: " << std::get<0>(seqName.second[varPerTrans.first].firstAminoInfo_).cDNAPos_ << std::endl;
//								//subtract off the start location
//								refCDnaPos -= std::get<0>(seqName.second[varPerTrans.first].firstAminoInfo_).cDNAPos_;
////								seqName.second[varPerTrans.first].cDna_.outPutSeqAnsi(std::cout);
////								seqName.second[varPerTrans.first].refAlnTranslation_.outPutSeq(std::cout);
////								seqName.second[varPerTrans.first].queryAlnTranslation_.outPutSeq(std::cout);
////								std::cout << "std::get<0>(seqName.second[varPerTrans.first].firstAminoInfo_).cDNAPos_;: " << std::get<0>(seqName.second[varPerTrans.first].firstAminoInfo_).cDNAPos_ << std::endl;
////								std::cout << "std::get<0>(seqName.second[varPerTrans.first].lastAminoInfo_).cDNAPos_;: " << std::get<0>(seqName.second[varPerTrans.first].lastAminoInfo_).cDNAPos_ << std::endl;
////								std::cout << "std::get<0>(seqName.second[varPerTrans.first].firstAminoInfo_).aaPos_;: " << std::get<0>(seqName.second[varPerTrans.first].firstAminoInfo_).aaPos_ << std::endl;
////								std::cout << "std::get<0>(seqName.second[varPerTrans.first].lastAminoInfo_).aaPos_;: " << std::get<0>(seqName.second[varPerTrans.first].lastAminoInfo_).aaPos_ << std::endl;
//								uint32_t alnPosForLoc = getAlnPosForRealPos(seqName.second[varPerTrans.first].refAlnTranslation_.seq_,loc);
//								uint32_t queryCDNAPos = getRealPosForAlnPos(seqName.second[varPerTrans.first].queryAlnTranslation_.seq_,alnPosForLoc ) * 3;
//
//								if(queryCDNAPos +2 >=seqName.second[varPerTrans.first].cDna_.seq_.size()){
//									std::stringstream ss;
//									ss << __PRETTY_FUNCTION__ << ", error " << "" << queryCDNAPos +2 << " is greater than seqName.second[varPerTrans.first].cDna_.seq_.size(): " << seqName.second[varPerTrans.first].cDna_.seq_.size()<< "\n";
//									throw std::runtime_error{ss.str()};
//								}
//								fullAATypedWithCodonInfo[popName].emplace_back(
//										TranslatorByAlignment::AAInfo(varPerTrans.first, aa, loc,
//												std::tie(
//														seqName.second[varPerTrans.first].cDna_.seq_[queryCDNAPos+ 0],
//														seqName.second[varPerTrans.first].cDna_.seq_[queryCDNAPos+ 1],
//														seqName.second[varPerTrans.first].cDna_.seq_[queryCDNAPos+ 2]),
//														njh::in(loc, knownMutationsLocations)));
//							}
//						}
//					}
//					if(!knownMutationsLocations.empty()){
//						OutputStream knownTabOut(OutOptions(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerTrans.first +  "-protein_aminoAcidKnownMutations.tab.txt"))));
//						knownTabOut << "transcript\tposition(1-based)\trefAA\tAA\tcount\tfraction\talleleDepth\tsamples" << std::endl;
//						for(const auto & snpPos : knownMutationsLocations){
//							for(const auto & aa : varPerTrans.second.allBases[snpPos]){
//								knownTabOut << varPerTrans.first
//										<< "\t" << snpPos + 1
//										<< "\t" << translatedRes.proteinForTranscript_[varPerTrans.first][snpPos]
//										<< "\t" << aa.first
//										<< "\t" << aa.second
//										<< "\t" << aa.second/static_cast<double>(varPerTrans.second.depthPerPosition[snpPos])
//										<< "\t" << varPerTrans.second.depthPerPosition[snpPos]
//										<< "\t" << samplesCalled << std::endl;
//							}
//						}
//					}
//				}
//			}
//
//			{
//				//snps
//				for( auto & varPerChrom : translatedRes.seqVariants_){
//					auto snpsPositions = getVectorOfMapKeys(varPerChrom.second.snpsFinal);
//					njh::sort(snpsPositions);
//					OutputStream snpTabOut(OutOptions(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerChrom.first +  "-SNPs.tab.txt"))));
//					snpTabOut << "chromosome\tposition(0-based)\trefBase\tbase\tcount\tfraction\talleleDepth\tsamples" << std::endl;
//					for(const auto & snpPos : snpsPositions){
//						for(const auto & base : varPerChrom.second.allBases[snpPos]){
//							snpTabOut << varPerChrom.first
//									<< "\t" << snpPos
//									<< "\t" << translatedRes.baseForPosition_[varPerChrom.first][snpPos]
//									<< "\t" << base.first
//									<< "\t" << base.second
//									<< "\t" << base.second/static_cast<double>(varPerChrom.second.depthPerPosition[snpPos])
//									<< "\t" << varPerChrom.second.depthPerPosition[snpPos]
//									<< "\t" << samplesCalled << std::endl;
//						}
//					}
//					if(!knownAAMutsChromPositions[varPerChrom.first].empty()){
//						OutputStream knownAASNPTabOut(OutOptions(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerChrom.first +  "-knownAA_SNPs.tab.txt"))));
//						knownAASNPTabOut << "chromosome\tposition(0-based)\trefBase\tbase\tcount\tfraction\talleleDepth\tsamples" << std::endl;
//						for(const auto & snpPos : knownAAMutsChromPositions[varPerChrom.first]){
//							for(const auto & base : varPerChrom.second.allBases[snpPos]){
//								knownAASNPTabOut << varPerChrom.first
//										<< "\t" << snpPos
//										<< "\t" << translatedRes.baseForPosition_[varPerChrom.first][snpPos]
//										<< "\t" << base.first
//										<< "\t" << base.second
//										<< "\t" << base.second/static_cast<double>(varPerChrom.second.depthPerPosition[snpPos])
//										<< "\t" << varPerChrom.second.depthPerPosition[snpPos]
//										<< "\t" << samplesCalled << std::endl;
//							}
//						}
//					}
//					OutputStream allBasesTabOut(OutOptions(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerChrom.first +  "-allBases.tab.txt"))));
//					allBasesTabOut << "chromosome\tposition(0-based)\trefBase\tbase\tcount\tfraction\talleleDepth\tsamples" << std::endl;
//					for(const auto & snpPos : varPerChrom.second.allBases){
//						for(const auto & base : snpPos.second){
//							allBasesTabOut << varPerChrom.first
//									<< "\t" << snpPos.first
//									<< "\t" << translatedRes.baseForPosition_[varPerChrom.first][snpPos.first]
//									<< "\t" << base.first
//									<< "\t" << base.second
//									<< "\t" << base.second/static_cast<double>(varPerChrom.second.depthPerPosition[snpPos.first])
//									<< "\t" << varPerChrom.second.depthPerPosition[snpPos.first]
//									<< "\t" << samplesCalled << std::endl;
//						}
//					}
//					if(!varPerChrom.second.variablePositons_.empty()){
//						GenomicRegion variableRegion = varPerChrom.second.getVariableRegion();
//						OutputStream bedVariableRegionOut(OutOptions(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerChrom.first +  "-chromosome_variableRegion.bed"))));
//						bedVariableRegionOut << variableRegion.genBedRecordCore().toDelimStrWithExtra() << std::endl;
//					}
//				}
//			}
//		}
//
//		if(!allMetaKeys.empty()){
//			OutputStream metaOutPut(OutOptions(njh::files::make_path(setUp.pars_.directoryName_, "popHapMeta.tab.txt")));
//			VecStr columnNames = allMetaKeysVec;
//			metaOutPut << "Identifier\tPopName\tHapPopUIDCount" << "\t" << njh::conToStr(columnNames, "\t");
//			metaOutPut << std::endl;
//			for(const auto & popSeq : forwardStrandClusters){
//				for(const auto & inputSeq : popSeq.reads_){
//					MetaDataInName outMeta;
//					VecStr missingMetaFields;
//					if (MetaDataInName::nameHasMetaData(inputSeq->seqBase_.name_)) {
//						MetaDataInName seqMeta(inputSeq->seqBase_.name_);
//						outMeta = seqMeta;
//						for (const auto & mf : allMetaKeysVec) {
//							if (!njh::in(mf, seqMeta.meta_)) {
//								missingMetaFields.emplace_back(mf);
//							}
//						}
//					} else {
//						missingMetaFields = allMetaKeysVec;
//					}
//					for(const auto & mf : missingMetaFields){
//						outMeta.addMeta(mf, "NA");
//					}
//					metaOutPut << identifier
//							<< "\t" << popSeq.seqBase_.name_;
//					MetaDataInName popSeqMeta(popSeq.seqBase_.name_);
//					metaOutPut << "\t" << popSeqMeta.getMeta("HapPopUIDCount");
//					for(const auto & mf : allMetaKeysVec){
//						metaOutPut << '\t' << outMeta.meta_[mf];
//					}
//
//					metaOutPut<< std::endl;
//				}
//			}
//		}
//
//		//amino typing
//		{
//			OutputStream aminoMetaOutPut(OutOptions(njh::files::make_path(setUp.pars_.directoryName_, "popHapMetaWithAminoChanges.tab.txt")));
//
//			std::set<std::string> aminoTypingFields;
//			for (const auto & popHapTyped : fullAATyped) {
//				for (const auto & aminoMetaLevel : popHapTyped.second.meta_) {
//					aminoTypingFields.emplace(aminoMetaLevel.first);
//				}
//			}
//			if (!aminoTypingFields.empty()) {
//
//				aminoMetaOutPut << "PopName\tHapPopUIDCount";
//				if(!allMetaKeys.empty()){
//					aminoMetaOutPut << "\t" << njh::conToStr(allMetaKeys, "\t");
//				}
//				//aminoMetaOutPut << "\t" << njh::conToStr(aminoTypingFields, "\t");
//				aminoMetaOutPut << "\t" << "transcript"
//						<< "\t" << "AA"
//						<< "\t" << "AAPos(1-based)"
//						<< "\t" << "codon"
//						<< "\t" << "knownMutation(based-on-supplied)";
//
//				aminoMetaOutPut << std::endl;
//
//				for(const auto & popSeq : forwardStrandClusters){
//					for(const auto & inputSeq : popSeq.reads_){
//						MetaDataInName outMeta;
//						if(!allMetaKeys.empty()){
//							VecStr missingMetaFields;
//							if (MetaDataInName::nameHasMetaData(inputSeq->seqBase_.name_)) {
//								MetaDataInName seqMeta(inputSeq->seqBase_.name_);
//								outMeta = seqMeta;
//								for (const auto & mf : allMetaKeysVec) {
//									if (!njh::in(mf, seqMeta.meta_)) {
//										missingMetaFields.emplace_back(mf);
//									}
//								}
//							} else {
//								missingMetaFields = allMetaKeysVec;
//							}
//							for(const auto & mf : missingMetaFields){
//								outMeta.addMeta(mf, "NA");
//							}
//						}
//						for(const auto & aminoInfo : fullAATypedWithCodonInfo[popSeq.seqBase_.name_]){
//							aminoMetaOutPut << popSeq.seqBase_.name_;
//							MetaDataInName popSeqMeta(popSeq.seqBase_.name_);
//							aminoMetaOutPut << "\t" << popSeqMeta.getMeta("HapPopUIDCount");
//							if(!allMetaKeys.empty()){
//								for(const auto & mf : allMetaKeysVec){
//									aminoMetaOutPut << '\t' << outMeta.meta_[mf];
//								}
//							}
//							aminoMetaOutPut << "\t" << aminoInfo.transcriptName_
//									<< "\t" << aminoInfo.aa_
//									<< "\t" << aminoInfo.zeroBasedPos_ + 1
//									<< "\t" << std::get<0>(aminoInfo.codon_) << std::get<1>(aminoInfo.codon_) << std::get<2>(aminoInfo.codon_)
//									<< "\t" << njh::boolToStr(aminoInfo.knownMut_);
//							aminoMetaOutPut<< std::endl;
//						}
//					}
//				}
//			}
//		}
//
//
//
//
//		table readLensTab(readLens, VecStr{"length", "count"});
//		OutputStream readLensOut(njh::files::make_path(setUp.pars_.directoryName_, "readLengthCounts.tab.txt"));
//		readLensTab.outPutContents(readLensOut, "\t");
//		uint32_t mostCommonReadLen = 0;
//		uint32_t mostCommonReadLenCount = 0;
//		for(const auto & rlen : readLens){
//			if(rlen.second > mostCommonReadLenCount){
//				mostCommonReadLen = rlen.first;
//				mostCommonReadLenCount = rlen.second;
//			}
//		}
//		uint32_t countOfNotMostCommonReadLen = 0;
//		for(const auto & rlen : readLens){
//			if(rlen.first != mostCommonReadLen && uAbsdiff(mostCommonReadLen , rlen.first) >= lengthDiffForLengthPoly){
//				countOfNotMostCommonReadLen += rlen.second;
//			}
//		}
//		bool lengthPoly = false;
//		if((countOfNotMostCommonReadLen/static_cast<double>(totalInputSeqs)) > variantCallerRunPars.lowVariantCutOff){
//			lengthPoly = true;
//		}
//		readVec::allSetFractionByTotalCount(forwardStrandClusters);
//		PopGenCalculator::DiversityMeasures divMeasures  = PopGenCalculator::getGeneralMeasuresOfDiversity(forwardStrandClusters);
//		//PairwisePairFactory pFacIfTotal(totalInput);
//		out << identifier
//				<< "\t" << totalInputSeqs
//				<< "\t" << forwardStrandClusters.size()
//				<< "\t" << samplesCalled
//				<< "\t" << divMeasures.singlets_
//				<< "\t" << divMeasures.doublets_
//				<< "\t" << divMeasures.expShannonEntropy_
//				<< "\t" << divMeasures.ShannonEntropyE_
//				<< "\t" << divMeasures.effectiveNumOfAlleles_
//				<< "\t" << divMeasures.heterozygostiy_
//				<< "\t" << (lengthPoly ? "true" : "false")
//				<< "\t" << numVariants
//				<< "\t" << variableRegionID
//				<< "\t" << avgGCContent
//				//<< "\t" << sumOfPairwiseDist/pFac.totalCompares_
//				<< std::endl;
//				//<< "\t" << sumOfPairwiseDistWeighted/pFacIfTotal.totalCompares_<< std::endl;
//	}
//
//	auto variantCallsDir = njh::files::make_path(setUp.pars_.directoryName_, "variantCalls");
//	njh::files::makeDir(njh::files::MkdirPar{variantCallsDir});
//
//	if(calculatingFst){
//		OutputStream out(njh::files::make_path(setUp.pars_.directoryName_, popMeta + "_basicInfo.tab.txt"));
//		out << "id\t" << popMeta << "\ttotalHaplotypes\tuniqueHaplotypes\tnsamples\tsingletons\tdoublets\texpShannonEntropy\tShannonEntropyE\teffectiveNumOfAlleles\the\tlengthPolymorphism\tnvariantsAbove" << variantCallerRunPars.lowVariantCutOff * 100 << "%" << std::endl;
//		std::unordered_map<std::string, double> heForFields;
//
//		for(auto & field : uniqueSeqsByMetaForwardStrand){
//			double sumOfSquares = 0;
//			std::map<uint32_t, uint32_t> readLens;
//			njh::sort(*field.second);
//			uint32_t singletons = 0;
//			uint32_t doublets = 0;
//			double seqCountPerField = readVec::getTotalReadCount(*field.second);
//			for(const auto & cIter : *field.second){
//				if(1 == cIter.seqBase_.cnt_){
//					++singletons;
//				}
//				if(2 == cIter.seqBase_.cnt_){
//					++doublets;
//				}
//				sumOfSquares += std::pow(cIter.seqBase_.cnt_/seqCountPerField, 2);
//				readLens[len(cIter.seqBase_)]+= cIter.seqBase_.cnt_;
//			}
//			TranslatorByAlignment::VariantsInfo varInfo(refRegion.genBed3RecordCore(), seqWithID);
//
//		//	std::unordered_map<uint32_t, std::unordered_map<char, uint32_t>> snps;
//		//	std::unordered_map<uint32_t, std::unordered_map<std::string,uint32_t>> insertions;
//		//	std::unordered_map<uint32_t, std::unordered_map<std::string,uint32_t>> deletions;
//		//
//		//	std::unordered_map<uint32_t, std::unordered_map<char, uint32_t>> snpsFinal;
//		//	std::unordered_map<uint32_t, std::unordered_map<std::string,uint32_t>> insertionsFinal;
//		//	std::unordered_map<uint32_t, std::unordered_map<std::string,uint32_t>> deletionsFinal;
//			//
//
//			for(const auto & seq : *field.second){
//
//				alignerObj.alignCacheGlobal(forwardStrandRefSeq, seq);
//				alignerObj.profileAlignment(forwardStrandRefSeq, seq, false, false, false);
//				varInfo.addVariantInfo(
//						alignerObj.alignObjectA_.seqBase_.seq_,
//						alignerObj.alignObjectB_.seqBase_.seq_,
//						seq.seqBase_.cnt_,
//						njh::mapAt(samplesPerSeqPerMeta, field.first)[seq.seqBase_.seq_],
//						alignerObj.comp_,
//						refRegion.start_);
//			}
//			varInfo.setFinals(variantCallerRunPars);
//
//
//			uint32_t mostCommonReadLen = 0;
//			uint32_t mostCommonReadLenCount = 0;
//			for(const auto & rlen : readLens){
//				if(rlen.second > mostCommonReadLenCount){
//					mostCommonReadLen = rlen.first;
//					mostCommonReadLenCount = rlen.second;
//				}
//			}
//			uint32_t countOfNotMostCommonReadLen = 0;
//			for(const auto & rlen : readLens){
//				if(rlen.first != mostCommonReadLen && uAbsdiff(mostCommonReadLen , rlen.first) >= lengthDiffForLengthPoly){
//					countOfNotMostCommonReadLen += rlen.second;
//				}
//			}
//
//			bool lengthPoly = false;
//			if((countOfNotMostCommonReadLen/static_cast<double>(seqCountPerField)) > variantCallerRunPars.lowVariantCutOff){
//				lengthPoly = true;
//			}
//
//			double he = 1 - sumOfSquares;
//			heForFields[field.first] = he;
//			auto divMeasures = PopGenCalculator::getGeneralMeasuresOfDiversity(*field.second);
//			out << identifier
//					<< "\t" << field.first
//					<< "\t" << seqCountPerField
//					<< "\t" << field.second->size()
//					<< "\t" << samplesByMeta.at(field.first).size()
//					<< "\t" << divMeasures.singlets_
//					<< "\t" << divMeasures.doublets_
//					<< "\t" << divMeasures.expShannonEntropy_
//					<< "\t" << divMeasures.ShannonEntropyE_
//					<< "\t" << divMeasures.effectiveNumOfAlleles_
//					<< "\t" << divMeasures.heterozygostiy_
//					<< "\t" << (lengthPoly ? "true" : "false")
//					<< "\t" << varInfo.snpsFinal.size() + varInfo.deletionsFinal.size() + varInfo.insertionsFinal.size() << std::endl;
//
//			if(!varInfo.snpsFinal.empty() || ! varInfo.deletionsFinal.empty() || !varInfo.insertionsFinal.empty()){
//				OutputStream vcfOut(njh::files::make_path(variantCallsDir, "variants_" + field.first + ".vcf"));
//				std::unordered_set<uint32_t> positionsSet;
//				for(const auto & snps : varInfo.snpsFinal){
//					positionsSet.emplace(snps.first);
//				}
//				for(const auto & ins : varInfo.insertionsFinal){
//					positionsSet.emplace(ins.first);
//				}
//				for(const auto & del : varInfo.deletionsFinal){
//					positionsSet.emplace(del.first);
//				}
//				vcfOut << "##fileformat=VCFv4.0" << std::endl;
//				vcfOut << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Allele Depth\">" << std::endl;
//				vcfOut << "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">" << std::endl;
//				vcfOut << "##INFO=<ID=AF,Number=.,Type=Float,Description=\"Allele Frequency\">" << std::endl;
//				vcfOut << "##INFO=<ID=AC,Number=.,Type=Integer,Description=\"Allele Count\">" << std::endl;
//				vcfOut << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" << std::endl;
//				std::vector<uint32_t> positions(positionsSet.begin(), positionsSet.end());
//				uint32_t samplesCalledForField = seqCountPerField;
//				if(!samplesByMeta[field.first].empty()){
//					samplesCalledForField = samplesByMeta[field.first].size();
//				}
//				njh::sort(positions);
////				if(refRegion.reverseSrand_){
////					njh::reverse(positions);
////				}
//				for(const auto & pos : positions){
//					vcfOut
//							<<  refRegion.chrom_
//							<< "\t" << pos + 1
//							<< "\t" << "."
//							<< "\t";
//					std::vector<std::string> alts;
//					std::vector<uint32_t> altsCounts;
//					std::vector<double> altsFreqs;
//					vcfOut << baseForPosition[pos] << "\t";
//					if(njh::in(pos, varInfo.snpsFinal)){
//						for(const auto & b : varInfo.snpsFinal[pos]){
//							alts.emplace_back(std::string(1, b.first));
//							altsCounts.emplace_back(b.second);
//							altsFreqs.emplace_back(b.second/static_cast<double>(varInfo.depthPerPosition[pos]));
//						}
//					}
//					if (njh::in(pos, varInfo.insertionsFinal)) {
//						for (const auto & ins : varInfo.insertionsFinal[pos]) {
//							alts.emplace_back(ins.first);
//							altsCounts.emplace_back(ins.second);
//							altsFreqs.emplace_back(ins.second/static_cast<double>(varInfo.depthPerPosition[pos]));
//						}
//					}
//					vcfOut << njh::conToStr(alts, ",")
//					<< "\t40\tPASS\t";
//					vcfOut
//							<< "DP=" << varInfo.depthPerPosition[pos] << ";"
//							<< "NS=" << samplesCalledForField << ";"
//							<< "AC=" << njh::conToStr(altsCounts, ",") << ";"
//							<< "AF=" << njh::conToStr(altsFreqs, ",")
//					<< std::endl;
//					if (njh::in(pos, varInfo.deletionsFinal)) {
//						for (const auto & d : varInfo.deletionsFinal[pos]) {
//							vcfOut <<  refRegion.chrom_
//									<< "\t" << pos + 1
//									<< "\t" << "."
//									<< "\t";
//							vcfOut << baseForPosition[pos] << d.first
//							<< "\t" << baseForPosition[pos] << "\t";
//							vcfOut << "40\tPASS\t";
//							vcfOut
//									<< "DP=" << varInfo.depthPerPosition[pos] << ";"
//									<< "NS=" << samplesCalledForField << ";"
//									<< "AC=" << d.second << ";"
//									<< "AF=" << d.second/static_cast<double>(varInfo.depthPerPosition[pos])
//							<< std::endl;
//						}
//					}
//				}
//			}
//		}
//
//
//		OutputStream pairwiseDiffsOut(njh::files::make_path(setUp.pars_.directoryName_, "pairwisePopDiffMeasures.tab.txt"));
//		pairwiseDiffsOut << "identifier"
//				<< "\t" << popMeta << "1"
//				<< "\t" << "popMeta" << "1_totalHaps"
//				<< "\t" << "popMeta" << "1_uniqueHaps"
//				<< "\t" << "popMeta" << "1_samples"
//				<< "\t" << popMeta << "2"
//				<< "\t" << "popMeta" << "2_totalHaps"
//				<< "\t" << "popMeta" << "2_uniqueHaps"
//				<< "\t" << "popMeta" << "2_samples"
//				<< "\t" << "HsSample"
//								<<"\t"<<"HsEst"
//								<<"\t"<<"HtSample"
//								<<"\t"<<"HtEst"
//								<<"\t"<<"Gst"
//								<<"\t"<<"GstEst"
//								<<"\t"<<"JostD"
//								<<"\t"<<"JostDEst"
//								<<"\t"<<"ChaoA"
//								<<"\t"<<"ChaoB"
//								<<"\t"<<"JostDChaoEst" << std::endl;
//
//		auto pairwiseMeasurements = PopGenCalculator::getPairwisePopDiff(uniqueSeqsByMetaForwardStrand);
//		for(const auto & field1 : pairwiseMeasurements){
//			for(const auto & field2 : field1.second){
//				uint32_t total1 = 0;
//				uint32_t total2 = 0;
//				for(const auto & haps1 : *uniqueSeqsByMetaForwardStrand.at(field1.first)){
//					total1 += haps1.seqBase_.cnt_;
//				}
//				for(const auto & haps2 : *uniqueSeqsByMetaForwardStrand.at(field2.first)){
//					total2 += haps2.seqBase_.cnt_;
//				}
//				pairwiseDiffsOut << identifier
//						<< "\t" << field1.first
//						<< "\t" << total1
//						<< "\t" << uniqueSeqsByMetaForwardStrand.at(field1.first)->size()
//						<< "\t" << samplesByMeta.at(field1.first).size()
//						<< "\t" << field2.first
//						<< "\t" << total2
//						<< "\t" << uniqueSeqsByMetaForwardStrand.at(field2.first)->size()
//						<< "\t" << samplesByMeta.at(field2.first).size()
//				    << "\t" << field2.second.hsSample_
//						<< "\t" << field2.second.hsEst_
//						<< "\t" << field2.second.htSample_
//						<< "\t" << field2.second.htEst_
//						<< "\t" << field2.second.gst_
//						<< "\t" << field2.second.gstEst_
//						<< "\t" << field2.second.jostD_
//						<< "\t" << field2.second.jostDEst_
//						<< "\t" << field2.second.chaoA_
//						<< "\t" << field2.second.chaoB_
//						<< "\t" << field2.second.jostDChaoEst_
//						<< std::endl;
//			}
//		}
//		OutputStream popDiffOut(njh::files::make_path(setUp.pars_.directoryName_, "popDiffMeasures.tab.txt"));
//		popDiffOut << "identifier"
//				<<"\t"<<"totalHaps"
//				<<"\t"<<"uniqueHaps"
//				<<"\t"<<"nsamples"
//				<<"\t"<<"HsSample"
//				<<"\t"<<"HsEst"
//				<<"\t"<<"HtSample"
//				<<"\t"<<"HtEst"
//				<<"\t"<<"Gst"
//				<<"\t"<<"GstEst"
//				<<"\t"<<"JostD"
//				<<"\t"<<"JostDEst"
//				<<"\t"<<"ChaoA"
//				<<"\t"<<"ChaoB"
//				<<"\t"<<"JostDChaoEst" << std::endl;
//
//		//PopGenCalculator::PopDifferentiationMeasures temp;
//		auto diffMeasures = PopGenCalculator::getOverallPopDiffForSeqs(uniqueSeqsByMetaForwardStrand);
//
//		popDiffOut << identifier
//
//				<< "\t" << totalInputSeqs
//				<< "\t" << forwardStrandClusters.size()
//				<< "\t" << samplesCalled
//		    << "\t" << diffMeasures.hsSample_
//				<< "\t" << diffMeasures.hsEst_
//				<< "\t" << diffMeasures.htSample_
//				<< "\t" << diffMeasures.htEst_
//				<< "\t" << diffMeasures.gst_
//				<< "\t" << diffMeasures.gstEst_
//				<< "\t" << diffMeasures.jostD_
//				<< "\t" << diffMeasures.jostDEst_
//				<< "\t" << diffMeasures.chaoA_
//				<< "\t" << diffMeasures.chaoB_
//				<< "\t" << diffMeasures.jostDChaoEst_
//				<< std::endl;
//	}
//
//	alignerObj.processAlnInfoOutput(setUp.pars_.outAlnInfoDirName_, setUp.pars_.verbose_);
//	if(!varInfo.snpsFinal.empty() || ! varInfo.deletionsFinal.empty() || !varInfo.insertionsFinal.empty()){
//		OutputStream snpTabOut(njh::files::make_path(variantCallsDir, "allSNPs.tab.txt"));
//		OutputStream vcfOut(njh::files::make_path(variantCallsDir, "allVariants.vcf"));
//		std::unordered_set<uint32_t> positionsSet;
//		for(const auto & snps : varInfo.snpsFinal){
//			positionsSet.emplace(snps.first);
//		}
//		for(const auto & ins : varInfo.insertionsFinal){
//			positionsSet.emplace(ins.first);
//		}
//		for(const auto & del : varInfo.deletionsFinal){
//			positionsSet.emplace(del.first);
//		}
//		vcfOut << "##fileformat=VCFv4.0" << std::endl;
//		vcfOut << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Allele Depth\">" << std::endl;
//		vcfOut << "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">" << std::endl;
//		vcfOut << "##INFO=<ID=AF,Number=.,Type=Float,Description=\"Allele Frequency\">" << std::endl;
//		vcfOut << "##INFO=<ID=AC,Number=.,Type=Integer,Description=\"Allele Count\">" << std::endl;
//		vcfOut << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" << std::endl;
//		snpTabOut << "chrom\tposition\tref\tvariant\tcount\tfrequency\talleleDepth\tsamples" << std::endl;
//
//		std::vector<uint32_t> positions(positionsSet.begin(), positionsSet.end());
//
//		njh::sort(positions);
////		if(refRegion.reverseSrand_){
////			njh::reverse(positions);
////		}
//		for(const auto & pos : positions){
//			if (njh::in(pos, varInfo.insertionsFinal) || njh::in(pos, varInfo.snpsFinal)) {
//				vcfOut <<  refRegion.chrom_
//						<< "\t" << pos+ 1
//						<< "\t" << "."
//						<< "\t";
//				std::vector<std::string> alts;
//				std::vector<uint32_t> altsCounts;
//				std::vector<double> altsFreqs;
//
//
//				vcfOut << baseForPosition[pos] << "\t";
//				if(njh::in(pos, varInfo.snpsFinal)){
//					uint32_t snpCount = 0;
//					for(const auto & b : varInfo.snpsFinal[pos]){
//						snpTabOut << refRegion.chrom_
//								<< "\t" << pos
//								<< "\t" << baseForPosition[pos]
//								<< "\t" << std::string(1, b.first)
//								<< "\t" << b.second
//								<< "\t" << b.second/static_cast<double>(varInfo.depthPerPosition[pos])
//								<< "\t" << varInfo.depthPerPosition[pos]
//								<< "\t" << samplesCalled << std::endl;
//						snpCount+= b.second;
//						alts.emplace_back(std::string(1, b.first));
//						altsCounts.emplace_back(b.second);
//						altsFreqs.emplace_back(b.second/static_cast<double>(varInfo.depthPerPosition[pos]));
//					}
//					snpTabOut << refRegion.chrom_
//							<< "\t" << pos
//							<< "\t" << baseForPosition[pos]
//							<< "\t" << baseForPosition[pos]
//							<< "\t" << varInfo.depthPerPosition[pos] - snpCount
//							<< "\t" << (varInfo.depthPerPosition[pos] - snpCount)/static_cast<double>(varInfo.depthPerPosition[pos])
//							<< "\t" << varInfo.depthPerPosition[pos]
//							<< "\t" << samplesCalled << std::endl;
//				}
//				if (njh::in(pos, varInfo.insertionsFinal)) {
//					for (const auto & ins : varInfo.insertionsFinal[pos]) {
//						alts.emplace_back(ins.first);
//						altsCounts.emplace_back(ins.second);
//						altsFreqs.emplace_back(ins.second/static_cast<double>(varInfo.depthPerPosition[pos]));
//					}
//				}
//
//				vcfOut << njh::conToStr(alts, ",")
//				<< "\t40\tPASS\t";
//				vcfOut
//						<< "DP=" << varInfo.depthPerPosition[pos] << ";"
//						<< "NS=" << samplesCalled << ";"
//						<< "AC=" << njh::conToStr(altsCounts, ",") << ";"
//						<< "AF=" << njh::conToStr(altsFreqs, ",")
//				<< std::endl;
//			}
//			if (njh::in(pos, varInfo.deletionsFinal)) {
//				for (const auto & d : varInfo.deletionsFinal[pos]) {
//					vcfOut <<  refRegion.chrom_
//							<< "\t" << pos+ 1
//							<< "\t" << "."
//							<< "\t";
//					vcfOut << baseForPosition[pos] << d.first
//					<< "\t" << baseForPosition[pos] << "\t";
//					vcfOut << "40\tPASS\t";
//					vcfOut
//							<< "DP=" << varInfo.depthPerPosition[pos] << ";"
//							<< "NS=" << samplesCalled << ";"
//							<< "AC=" << d.second << ";"
//							<< "AF=" << d.second/static_cast<double>(varInfo.depthPerPosition[pos])
//					<< std::endl;
//				}
//			}
//		}
//	}
//	return 0;
}


int popGenExpRunner::oldQuickHaplotypeInformationAndVariants(const njh::progutils::CmdArgs & inputCommands) {
//
//	/**@todo should add the following 1) average pairwise difference, 2) make into a function, 3) generate connected hap map, 4) unique haps to region count, 5) doing multiple pop fields at once */
//
//	std::string identifier = "";
//	std::string popMeta = "";
//	double lowVariantCutOff = 0.005;
//	uint32_t occurrenceCutOff = 1;
//	uint32_t lengthDiffForLengthPoly = 6;
//	bfs::path bedFnp = "";
//	bfs::path genomefnp = "";
//	bfs::path gffFnp = "";
//	uint32_t outwardsExpand = 5;
//	bfs::path proteinMutantTypingFnp = "";
//	bool zeroBased = false;
//	bool keepNonFieldSamples = false;
//	std::set<std::string> ignoreSubFields;
//	MultiGenomeMapper::inputParameters gPars;
//	seqSetUp setUp(inputCommands);
//	setUp.processVerbose();
//	setUp.processDebug();
//	setUp.processReadInNames(true);
//	setUp.processDirectoryOutputName(true);
//	setUp.processAlnInfoInput();
//	setUp.setOption(lengthDiffForLengthPoly, "--lengthDiffForLengthPoly", "If there is a length difference of this or above, label as there being length polymorphism");
//	setUp.setOption(occurrenceCutOff, "--occurrenceCutOff", "Occurrence Cut Off, don't report variants below this count");
//	setUp.setOption(lowVariantCutOff, "--lowVariantCutOff", "Low Variant Cut Off, don't report variants below this fraction");
//
//	setUp.setOption(outwardsExpand, "--outwardsExpand", "The amount to expand outwards from given region when determining variants positions with extracted ref seq");
//
//	setUp.setOption(bedFnp,    "--bed",    "A bed file of the location for the extraction", true);
//	setUp.setOption(genomefnp, "--genome", "A reference genome to compare against", true);
//	setUp.setOption(gffFnp, "--gff", "A gff3 file for genome file");
//	if(!bfs::is_regular_file(genomefnp)){
//		setUp.failed_ = true;
//		setUp.addWarning(njh::pasteAsStr(genomefnp, " should be a file, not a directory"));
//	}
//	gPars.genomeDir_ = njh::files::normalize(genomefnp.parent_path());
//	gPars.primaryGenome_ = bfs::basename(genomefnp);
//	gPars.primaryGenome_ = gPars.primaryGenome_.substr(0, gPars.primaryGenome_.rfind("."));
//	gPars.selectedGenomes_ = {gPars.primaryGenome_};
//
//	setUp.setOption(ignoreSubFields, "--ignoreSubFields", "Meta Sub Fields to ignore when doing population analysis");
//
//	setUp.setOption(identifier, "--identifier", "Give a identifier name for info");
//	setUp.setOption(popMeta,    "--popMeta",    "Meta field to calculate population stats by");
//	setUp.setOption(proteinMutantTypingFnp, "--proteinMutantTypingFnp", "Protein Mutant Typing Fnp, columns should be ID=gene id in gff, AAPosition=amino acid position", "" != gffFnp);
//	setUp.setOption(zeroBased, "--zeroBased", "If the positions in the proteinMutantTypingFnp are zero based");
//	setUp.setOption(keepNonFieldSamples, "--keepNonFieldSamples", "Keep Non Field Samples for population stats");
//
//	setUp.finishSetUp(std::cout);
//	setUp.startARunLog(setUp.pars_.directoryName_);
//
//	njh::files::checkExistenceThrow(bedFnp, __PRETTY_FUNCTION__);
//	njh::files::checkExistenceThrow(genomefnp, __PRETTY_FUNCTION__);
//
//	uint64_t maxLen = 0;
//	SeqInput reader(setUp.pars_.ioOptions_);
//	reader.openIn();
//
//	std::unordered_map<std::string, std::shared_ptr<std::vector<identicalCluster>>> uniqueSeqsByMeta;
//	std::vector<identicalCluster> clusters;
//	seqInfo seq;
//	uint32_t totalInput = 0;
//	bool calculatingFst = "" != popMeta;
//	std::set<std::string> samples;
//	std::set<std::string> allMetaKeys;
//	std::unordered_map<std::string, std::set<std::string>> samplesByMeta;
//
//	// read in reads and collapse to unique
//	while(reader.readNextRead(seq)) {
//		//get meta keys if available
//		if(MetaDataInName::nameHasMetaData(getSeqBase(seq).name_)){
//			MetaDataInName metaData(getSeqBase(seq).name_);
//			for(const auto & meta : metaData.meta_){
//				allMetaKeys.emplace(meta.first);
//			}
//		}
//
//		readVec::handelLowerCaseBases(seq, setUp.pars_.ioOptions_.lowerCaseBases_);
//		if(setUp.pars_.ioOptions_.removeGaps_){
//			seq.removeGaps();
//		}
//		readVec::getMaxLength(seq, maxLen);
//		if(seq.nameHasMetaData()){
//			MetaDataInName seqMeta(seq.name_);
//			if(seqMeta.containsMeta("sample")){
//				samples.emplace(seqMeta.getMeta("sample"));
//			}
//		}
//		++totalInput;
//		bool found = false;
//		for (auto &cIter : clusters) {
//			if (cIter.seqBase_.seq_ == seq.seq_) {
//				cIter.addRead(seq);
//				found = true;
//				break;
//			}
//		}
//		if (!found) {
//			clusters.emplace_back(seq);
//		}
//
//		//if calculating sub population differences collect unique sequences for the given field
//		if(calculatingFst){
//			MetaDataInName seqMeta(seq.name_);
//			bool found = false;
//			//skip ignore fields
//			if(njh::in(seqMeta.getMeta(popMeta), ignoreSubFields)){
//				continue;
//			}
//			//skip non field samples
//			if(!keepNonFieldSamples && seqMeta.containsMeta("IsFieldSample") && !seqMeta.getMeta<bool>("IsFieldSample")){
//				continue;
//			}
//			if(!njh::in(seqMeta.getMeta(popMeta),  uniqueSeqsByMeta)){
//				uniqueSeqsByMeta[seqMeta.getMeta(popMeta)] = std::make_shared<std::vector<identicalCluster>>();
//			}
//			if(seqMeta.containsMeta("sample")){
//				samplesByMeta[seqMeta.getMeta(popMeta)].emplace(seqMeta.getMeta("sample"));
//			}
//
//			for (auto &cIter : *uniqueSeqsByMeta[seqMeta.getMeta(popMeta)]) {
//				if (cIter.seqBase_.seq_ == seq.seq_) {
//					cIter.addRead(seq);
//					found = true;
//					break;
//				}
//			}
//			if (!found) {
//				uniqueSeqsByMeta[seqMeta.getMeta(popMeta)]->emplace_back(seq);
//			}
//		}
//	}
//
//	double sumGCContent = 0;
//	for(auto & clus : clusters){
//		clus.setLetterCount();
//		clus.counter_.calcGcContent();
//		sumGCContent += clus.counter_.gcContent_ * clus.seqBase_.cnt_;
//	}
//	double avgGCContent = sumGCContent/totalInput;
//
//
//	VecStr allMetaKeysVec(allMetaKeys.begin(), allMetaKeys.end());
//
//	uint32_t samplesCalled = totalInput;
//	if(!samples.empty()){
//		samplesCalled = samples.size();
//	}
//
//
//
//	auto beds = getBeds(bedFnp);
//	if(beds.empty()){
//		std::stringstream ss;
//		ss << __PRETTY_FUNCTION__ << ", error no records found in " << bedFnp << "\n";
//		throw std::runtime_error{ss.str()};
//	}
//	MultiGenomeMapper gMapper(gPars);
//	gMapper.loadInGenomes();
//	gMapper.setUpGenomes();
//	auto gPos = beds.front();
//
//	//get identifer for sequences
//	if("" == identifier){
//		identifier = gPos->name_;
//	}
//	njh::sort(clusters);
//	uint32_t seqId = 0;
//	std::unordered_map<std::string, std::string> nameLookUp;
//	for (auto &cIter : clusters) {
//		MetaDataInName popMeta;
//		popMeta.addMeta("HapPopUIDCount", static_cast<uint32_t>(std::round(cIter.seqBase_.cnt_)));
//		cIter.seqBase_.name_ = njh::pasteAsStr(identifier, ".", leftPadNumStr<uint32_t>(seqId, clusters.size()),popMeta.createMetaName());
//		nameLookUp[cIter.seqBase_.seq_] = cIter.seqBase_.name_;
//		++seqId;
//	}
//
//	if(calculatingFst){
//		std::unordered_map<std::string, uint32_t> seqIdForField;
//		for(auto & field : uniqueSeqsByMeta){
//			readVec::allSetFractionByTotalCount(*field.second);
//			for (auto &cIter : *field.second) {
//				MetaDataInName popMeta;
//				std::string popName = nameLookUp[cIter.seqBase_.seq_];
//				MetaDataInName popNameMeta(popName);
//				MetaDataInName::removeMetaDataInName(popName);
//				popMeta.addMeta("HapPopUID", popName);
//				popMeta.addMeta(popNameMeta, true);
//				popMeta.addMeta("HapUniqueCount", static_cast<uint32_t>(std::round(cIter.seqBase_.cnt_)));
//				cIter.seqBase_.name_ = njh::pasteAsStr(field.first, ".", leftPadNumStr<uint32_t>(seqIdForField[field.first], field.second->size()), popMeta.createMetaName());
//				++seqIdForField[field.first];
//			}
//		}
//	}
//
//	//get region and ref seq for mapping of variants
//	OutputStream inputRegionBeforeExpandOut(OutOptions(njh::files::make_path(setUp.pars_.directoryName_, "inputRegionBeforeExpand.bed")));
//	inputRegionBeforeExpandOut << gPos->toDelimStrWithExtra() << std::endl;
//
//	uint32_t oldLen = gPos->length();
//	BedUtility::extendLeftRight(*gPos, outwardsExpand, outwardsExpand,
//			gMapper.genomes_.at(gPars.primaryGenome_)->chromosomeLengths_.at(
//					gPos->chrom_));
//	if(oldLen == gPos->score_){
//		gPos->score_ = gPos->length();
//	}
//
//	OutputStream bedOut(OutOptions(njh::files::make_path(setUp.pars_.directoryName_, "inputRegion.bed")));
//	bedOut << gPos->toDelimStrWithExtra() << std::endl;
//
//	GenomicRegion refRegion(*gPos);
//
//	TwoBit::TwoBitFile tReader(
//			gMapper.genomes_.at(gPars.primaryGenome_)->fnpTwoBit_);
//	auto refSeq = refRegion.extractSeq(tReader);
//	SeqOutput::write(std::vector<seqInfo> { refSeq },
//			SeqIOOptions::genFastaOut(
//					njh::files::make_path(setUp.pars_.directoryName_,
//							"inputRegion.fasta")));
//
//	readVec::getMaxLength(refSeq, maxLen);
//	aligner alignerObj(maxLen, gapScoringParameters(5,1,0,0,0,0), substituteMatrix(2,-2), false);
//	alignerObj.weighHomopolymers_ = false;
//	alignerObj.processAlnInfoInput(setUp.pars_.alnInfoDirName_, setUp.pars_.verbose_);
//
//	std::unordered_map<uint32_t, std::unordered_map<char, uint32_t>> snps;
//	std::unordered_map<uint32_t, std::unordered_map<std::string,uint32_t>> insertions;
//	std::unordered_map<uint32_t, std::unordered_map<std::string,uint32_t>> deletions;
//
//	std::unordered_map<uint32_t, std::unordered_map<char, uint32_t>> snpsFinal;
//	std::unordered_map<uint32_t, std::unordered_map<std::string,uint32_t>> insertionsFinal;
//	std::unordered_map<uint32_t, std::unordered_map<std::string,uint32_t>> deletionsFinal;
//	//
//
//	for(const auto & seq : clusters){
//		alignerObj.alignCacheGlobal(refSeq, seq);
//		alignerObj.profileAlignment(refSeq, seq, false, false, false);
//		for(const auto & m : alignerObj.comp_.distances_.mismatches_){
//			snps[m.second.refBasePos][m.second.seqBase]+= seq.seqBase_.cnt_;
//		}
//		for(const auto & g : alignerObj.comp_.distances_.alignmentGaps_){
//			if(g.second.ref_){
//				//insertion
//				insertions[g.second.refPos_ - 1][g.second.gapedSequence_]++;
//			}else{
//				//deletion
//				deletions[g.second.refPos_ - 1][g.second.gapedSequence_]++;
//			}
//		}
//	}
//
//	std::vector<uint32_t> variablePositions;
//	//filter snps and indels by occurrence cut off
//	for(const auto & snp : snps){
//		for(const auto & b : snp.second){
//			if(b.second < occurrenceCutOff){
//				continue;
//			}
//			snpsFinal[snp.first][b.first] = b.second;
//			if(b.second/static_cast<double>(totalInput) > lowVariantCutOff){
//				variablePositions.emplace_back(snp.first);
//			}
//		}
//	}
//	for(const auto & del : deletions){
//		for(const auto & d : del.second){
//			if(d.second < occurrenceCutOff){
//				continue;
//			}
//			deletionsFinal[del.first][d.first] = d.second;
//			if(d.second/static_cast<double>(totalInput) > lowVariantCutOff){
//				variablePositions.emplace_back(del.first);
//			}
//		}
//	}
//	for(const auto & ins : insertions){
//		for(const auto & i : ins.second){
//			if(i.second < occurrenceCutOff){
//				continue;
//			}
//			insertionsFinal[ins.first][i.first] = i.second;
//			if(i.second/static_cast<double>(totalInput) > lowVariantCutOff){
//				variablePositions.emplace_back(ins.first);
//			}
//		}
//	}
//	std::string variableRegionID = "";
//	if(!variablePositions.empty()){
//		uint32_t variableStart = vectorMinimum(variablePositions);
//		uint32_t variableStop = vectorMaximum(variablePositions);
//
//		size_t genomeVarStart = refRegion.getRelativePositionFromStartStrandAware(variableStart);
//		size_t genomeVarStop = refRegion.getRelativePositionFromStartStrandAware(variableStop);
//
//		GenomicRegion variableRegion = refRegion;
//		if(variableRegion.reverseSrand_){
//			variableRegion.start_ = genomeVarStop;
//			variableRegion.end_ = genomeVarStart + 1;
//		}else{
//			variableRegion.start_ = genomeVarStart;
//			variableRegion.end_ = genomeVarStop + 1;
//		}
//		OutputStream bedVariableRegionOut(OutOptions(njh::files::make_path(setUp.pars_.directoryName_, "variableRegion.bed")));
//		bedVariableRegionOut << variableRegion.genBedRecordCore().toDelimStrWithExtra() << std::endl;
//		variableRegionID = njh::pasteAsStr(variableRegion.chrom_, "::", variableRegion.start_, "::", variableRegion.end_);
//	}
//
//	{
//		uint32_t numSNPs = 0;
//		for(const auto & snp : snpsFinal){
//			uint32_t above = 0;
//			for(const auto & b : snp.second){
//				if(b.second/static_cast<double>(totalInput) > lowVariantCutOff){
//					++above;
//				}
//			}
//			if(above >= 1){
//				++numSNPs;
//			}
//		}
//		OutputStream out(njh::files::make_path(setUp.pars_.directoryName_, "basicInfo.tab.txt"));
//		out << "id\ttotalHaplotypes\tuniqueHaplotypes\tnsamples\tsingletons\tdoublets\texpShannonEntropy\tShannonEntropyE\teffectiveNumOfAlleles\the\tlengthPolymorphism\tnSNPsAbove"
//				<< lowVariantCutOff * 100 << "%" << "\tvariableRegionID\t" << "avgGCContent";
//		//out << "\tAvgPairwiseDistWeighted";
//		out << std::endl;
//		std::map<uint32_t, uint32_t> readLens;
//		//writing out unique sequences
//		auto uniqueSeqsOpts = SeqIOOptions::genFastaOut(njh::files::make_path(setUp.pars_.directoryName_, "uniqueSeqs.fasta"));
//		OutputStream nameOut(OutOptions(njh::files::make_path(setUp.pars_.directoryName_, "uniqueSeqs_names.tab.txt")));
//		nameOut << "name\tnumber\treads"	<< std::endl;
//		{
//			SeqOutput uniqueWriter(uniqueSeqsOpts);
//			uniqueWriter.openOut();
//			uniqueWriter.write(clusters);
//			uniqueWriter.closeOut();
//		}
//
////		PairwisePairFactory pFac(clusters.size());
////		PairwisePairFactory::PairwisePair pPair;
////		double sumOfPairwiseDist = 0;
////		double sumOfPairwiseDistWeighted = 0;
////
////		while (pFac.setNextPair(pPair)) {
////			alignerObj.alignCacheGlobal(clusters[pPair.col_], clusters[pPair.row_]);
////			alignerObj.profileAlignment(clusters[pPair.col_], clusters[pPair.row_],
////					false, false, false);
////			sumOfPairwiseDist += alignerObj.comp_.distances_.eventBasedIdentity_;
////			sumOfPairwiseDistWeighted +=
////					alignerObj.comp_.distances_.eventBasedIdentity_
////							* (clusters[pPair.col_].seqBase_.cnt_
////									* clusters[pPair.row_].seqBase_.cnt_);
////		}
//
//
//		std::unordered_map<std::string, std::unordered_map<std::string, GenomicAminoAcidPositionTyper::GeneAminoTyperInfo>> popHapsTyped;
//		std::unordered_map<std::string, std::unordered_map<std::string, GenomicAminoAcidPositionTyper::GeneAminoTyperInfo>> popHapsTypeddenovo;
//
//		if("" != gffFnp){
//			BioCmdsUtils bRunner(setUp.pars_.verbose_);
//			bRunner.RunFaToTwoBit(genomefnp);
//			bRunner.RunBowtie2Index(genomefnp);
//
//			auto gprefix = bfs::path(genomefnp).replace_extension("");
//			auto twoBitFnp = gprefix.string() + ".2bit";
//			auto uniqueSeqInOpts = SeqIOOptions::genFastaIn(uniqueSeqsOpts.out_.outName());
//			uniqueSeqInOpts.out_.outFilename_ = njh::files::make_path(setUp.pars_.directoryName_, "aligned_uniqueSeqs.sorted.bam");
//			uniqueSeqInOpts.out_.outExtention_ = ".sorted.bam";
//			auto bowtieRunOut = bRunner.bowtie2Align(uniqueSeqInOpts, genomefnp);
//
//			auto regionsCounter = GenomicRegionCounter::countRegionsInBam(uniqueSeqInOpts.out_.outName());
//			auto ids = regionsCounter.getIntersectingGffIds(gffFnp);
//			OutOptions idOpts(njh::files::make_path(setUp.pars_.directoryName_, "geneIds.txt"));
//			OutputStream idOut(idOpts);
//			idOut << njh::conToStr(ids, "\n") << std::endl;
//
//			//auto genes = GeneFromGffs::getGenesFromGffForIds(gffFnp, ids);
//
//
//			//initiate typer
//			GenomicAminoAcidPositionTyper aaTyper(proteinMutantTypingFnp, zeroBased);
//
//			// get gene information
//			auto geneInfoDir = njh::files::make_path(setUp.pars_.directoryName_, "geneInfos");
//			njh::files::makeDir(njh::files::MkdirPar{geneInfoDir});
//
//			OutOptions outOpts(njh::files::make_path(geneInfoDir, "gene"));
//
//			std::unordered_map<std::string, VecStr> idToTranscriptName;
//			std::unordered_map<std::string, std::shared_ptr<GeneFromGffs>> genes = GeneFromGffs::getGenesFromGffForIds(gffFnp, ids);
//			std::unordered_map<std::string, std::shared_ptr<GeneSeqInfo>> geneInfos;
//
//			uint64_t proteinMaxLen = 0;
//			std::unordered_map<std::string, std::unordered_map<std::string, std::shared_ptr<GeneFromGffs>>> genesByChrom;
//
//			for(const auto & gene : genes){
//				genesByChrom[gene.second->gene_->seqid_].emplace(gene.first, gene.second);
//				for(const auto & transcript : gene.second->mRNAs_){
//					idToTranscriptName[gene.second->gene_->getIDAttr()].emplace_back(transcript->getIDAttr());
//				}
//				geneInfos[gene.first] = gene.second->generateGeneSeqInfo(tReader, false).begin()->second;
//				readVec::getMaxLength(geneInfos[gene.first]->protein_, proteinMaxLen);
//				gene.second->writeOutGeneInfo(tReader, outOpts);
//			}
//			//check for multiple transcripts
//			/**@todo add inputing what transcript to avoid this check, cannot be handled now because the input positions can only refer to one transcript*/
//			bool failed = false;
//			VecStr idsWithMoreThanOneTranscript;
//			for(const auto & idTrans : idToTranscriptName){
//				if(idTrans.second.size() > 1){
//					failed = true;
//					idsWithMoreThanOneTranscript.emplace_back(idTrans.first);
//				}
//			}
//			if (failed) {
//				std::stringstream ss;
//				ss << __PRETTY_FUNCTION__
//						<< ", error the following ids were found to have more than 1 transcript which can't be handled in the current version"
//						<< "\n";
//				ss << "ids: " << njh::conToStr(idsWithMoreThanOneTranscript, ", ") << "\n";
//				throw std::runtime_error { ss.str() };
//			}
////			//check for missing ids from the input query
////			VecStr idsMissing;
////			for(const auto & id : aaTyper.getGeneIds()){
////				if(!njh::in(id, idToTranscriptName)){
////					idsMissing.emplace_back(id);
////					failed = true;
////				}
////			}
////
////			if (failed) {
////				std::stringstream ss;
////				ss << __PRETTY_FUNCTION__
////						<< ", the following ids were not found in the gff file, " << gffFnp
////						<< "\n";
////				ss << "ids: " << njh::conToStr(idsMissing, ", ") << "\n";
////				throw std::runtime_error { ss.str() };
////			}
//			//check to make sure the asked for positons are not out of range
//			//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
//			bool failedPosition = false;
//			std::stringstream ssAminoAcidPosCheck;
//			ssAminoAcidPosCheck << __PRETTY_FUNCTION__ <<  "\n" ;
//			for (const auto & geneId : aaTyper.aminoPositionsForTyping_) {
//				if(!njh::in(geneId.first, ids)){
//					continue;
//				}
//				std::unordered_map<std::string, std::shared_ptr<GeneSeqInfo>> gsInfos = genes[geneId.first]->generateGeneSeqInfo(tReader, false);
//				auto gsInfo = gsInfos.begin()->second;
//				std::map<uint32_t, char> refAminoAcids;
//				for (const auto & aaPos : geneId.second) {
//					if (aaPos >= gsInfo->protein_.seq_.size()) {
//						ssAminoAcidPosCheck << "amino acid position "
//								<< (zeroBased ? aaPos : aaPos + 1)
//								<< " is out of range of the gene id " << geneId.first
//								<< ", max position: "
//								<< (zeroBased ?
//										gsInfo->protein_.seq_.size() - 1 : gsInfo->protein_.seq_.size())
//								<< '\n';
//					} else {
//						refAminoAcids[aaPos] = gsInfo->protein_.seq_[aaPos];
//					}
//				}
//				aaTyper.aminoPositionsForTypingWithInfo_.emplace(geneId.first, GenomicAminoAcidPositionTyper::GeneAminoTyperInfo(geneId.first, refAminoAcids));
//				if(aaTyper.altNamesForIds_[geneId.first].size() == 1){
//					aaTyper.aminoPositionsForTypingWithInfo_.at(geneId.first).altId_ = *aaTyper.altNamesForIds_[geneId.first].begin();
//				}
//			}
//			if(failedPosition){
//				throw std::runtime_error{ssAminoAcidPosCheck.str()};
//			}
//			//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
//			aligner alignObj(proteinMaxLen, gapScoringParameters(5,1,0,0,0,0), substituteMatrix(2,-2));
//
//			//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
//			std::unordered_map<std::string, std::set<std::string>> regionsToGeneIds;
//			//targetName, GeneID, AA Position
//			std::unordered_map<std::string, std::vector<std::string>> alnRegionToGeneIds;
//			//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
//			for(const auto & gCount : regionsCounter.counts_){
//				for (const auto & g : genesByChrom[gCount.second.region_.chrom_]) {
//					if (gCount.second.region_.overlaps(*g.second->gene_)) {
//						alnRegionToGeneIds[gCount.second.region_.createUidFromCoords()].emplace_back(g.first);
//					}
//				}
//			}
//			{
//				//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
//					//typing by protein changes
//					BamTools::BamReader bReader;
//					//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
//					bReader.Open(uniqueSeqInOpts.out_.outName().string());
//					checkBamOpenThrow(bReader, uniqueSeqInOpts.out_.outName());
//					//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
//					auto refData = bReader.GetReferenceData();
//					//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
//					uint32_t numOfThreadsToUse = std::min<uint32_t>(refData.size(), 1);
//					//create bam readers
//					//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
//					concurrent::BamReaderPool bamPool(uniqueSeqInOpts.out_.outName(), numOfThreadsToUse);
//					bamPool.openBamFile();
//					//create alingers
//					concurrent::AlignerPool alnPool(alignObj, numOfThreadsToUse);
//					alnPool.initAligners();
//					//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
//					MultiSeqOutCache<seqInfo> proteinSeqOuts;
//					for(const auto & g :genes){
//						proteinSeqOuts.addReader(g.first, SeqIOOptions::genFastaOut(njh::files::make_path(setUp.pars_.directoryName_, g.first)));
//					}
//					//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
//					std::mutex proteinSeqOutsMut;
//					std::mutex transferInfoMut;
//					//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
//					auto chromRegions = genGenRegionsFromRefData(refData);
//					njh::concurrent::LockableVec<GenomicRegion> chromRegionsVec(chromRegions);
//
//					struct VariantsInfo {
//						std::unordered_map<uint32_t, std::unordered_map<char, uint32_t>> snps;
//						std::unordered_map<uint32_t, std::unordered_map<std::string,uint32_t>> insertions;
//						std::unordered_map<uint32_t, std::unordered_map<std::string,uint32_t>> deletions;
//
//						std::unordered_map<uint32_t, std::unordered_map<char, uint32_t>> snpsFinal;
//						std::unordered_map<uint32_t, std::unordered_map<std::string,uint32_t>> insertionsFinal;
//						std::unordered_map<uint32_t, std::unordered_map<std::string,uint32_t>> deletionsFinal;
//					};
//					std::unordered_map<std::string, VariantsInfo> variantInfoPerGene;
//
//					//
//					//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
//					std::function<void()> typeOnChrom = [&bamPool,&alnPool,&proteinSeqOuts,&proteinSeqOutsMut,
//															&chromRegionsVec,&refData,&genes,
//															&geneInfos,&alnRegionToGeneIds,&transferInfoMut,
//															&setUp,&twoBitFnp,&aaTyper,
//															&popHapsTyped,&variantInfoPerGene
//															](){
//
//						GenomicRegion currentChrom;
//						auto curAligner = alnPool.popAligner();
//						auto curBReader = bamPool.popReader();
//						TwoBit::TwoBitFile tReader(twoBitFnp);
//						BamTools::BamAlignment bAln;
//						//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
//						std::unordered_map<std::string, std::unordered_map<std::string, GenomicAminoAcidPositionTyper::GeneAminoTyperInfo>> currentPopHapsTyped;
//
//						while(chromRegionsVec.getVal(currentChrom)){
//							//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
//							setBamFileRegionThrow(*curBReader, currentChrom);
//							//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
//							while (curBReader->GetNextAlignment(bAln)) {
//								if (bAln.IsMapped() && bAln.IsPrimaryAlignment()) {
//									//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
//									auto balnGenomicRegion = GenomicRegion(bAln, refData);
//
//									if (setUp.pars_.verbose_) {
//										std::cout << balnGenomicRegion.genBedRecordCore().toDelimStr() << std::endl;
//									}
//									if(!njh::in(balnGenomicRegion.createUidFromCoords(), alnRegionToGeneIds)){
//										continue;
//									}
//									//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
//									for (const auto & g : alnRegionToGeneIds.at(balnGenomicRegion.createUidFromCoords())) {
//										//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
//										const auto & currentGene = genes.at(g);
//										const auto & currentGeneInfo = geneInfos.at(g);
//										//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
//										auto aminoTyping = aaTyper.typeAlignment(bAln, *currentGene, *currentGeneInfo, tReader, *curAligner, refData);
//										//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
//										MetaDataInName seqMeta(bAln.Name);
//										auto popCount = seqMeta.getMeta<uint32_t>("HapPopUIDCount");
//										for(const auto & m : curAligner->comp_.distances_.mismatches_){
//											variantInfoPerGene[g].snps[m.second.refBasePos][m.second.seqBase]+= popCount;
//										}
//										for(const auto & gap : curAligner->comp_.distances_.alignmentGaps_){
//											if(gap.second.ref_){
//												//insertion
//												variantInfoPerGene[g].insertions[gap.second.refPos_ - 1][gap.second.gapedSequence_]+=popCount;
//											}else{
//												//deletion
//												variantInfoPerGene[g].deletions[gap.second.refPos_ - 1][gap.second.gapedSequence_]+=popCount;
//											}
//										}
//										if (!aminoTyping.empty() && njh::in(g, aaTyper.aminoPositionsForTyping_)) {
//											//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
//											currentPopHapsTyped[bAln.Name].emplace(g, GenomicAminoAcidPositionTyper::GeneAminoTyperInfo(g, aminoTyping));
//											//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
//											auto typeMeta = MetaDataInName::mapToMeta(aminoTyping);
//											//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
//											curAligner->alignObjectB_.seqBase_.name_.append(typeMeta.createMetaName());
//										}
//										//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
//										if(setUp.pars_.debug_){
//											std::lock_guard<std::mutex> lock(proteinSeqOutsMut);
//											proteinSeqOuts.add(g, curAligner->alignObjectA_.seqBase_);
//											proteinSeqOuts.add(g, curAligner->alignObjectB_.seqBase_);
//										}
//									}
//									//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
//								}
//							}
//							{
//								std::lock_guard<std::mutex> lock(transferInfoMut);
//								for(const auto & popHapTyped : currentPopHapsTyped){
//									popHapsTyped.emplace(popHapTyped.first, popHapTyped.second);
//								}
//							}
//						}
//					};
//
//					njh::concurrent::runVoidFunctionThreaded(typeOnChrom, numOfThreadsToUse);
//
//					for(auto & geneVarInfo : variantInfoPerGene){
//						//std::cout << geneVarInfo.first << " : "<< "geneVarInfo.second.snps: " << geneVarInfo.second.snps.size() << std::endl;
//						std::vector<uint32_t> variablePositions;
//						//filter snps and indels by occurrence cut off
//						for(const auto & snp : geneVarInfo.second.snps){
//							for(const auto & b : snp.second){
//								//std::cout << "\t" << b.first << " " << b.second << std::endl;
//								if(b.second < occurrenceCutOff){
//									continue;
//								}
//								geneVarInfo.second.snpsFinal[snp.first][b.first] = b.second;
//								if(b.second/static_cast<double>(totalInput) > lowVariantCutOff){
//									variablePositions.emplace_back(snp.first);
//								}
//							}
//						}
//						for(const auto & del : geneVarInfo.second.deletions){
//							for(const auto & d : del.second){
//								if(d.second < occurrenceCutOff){
//									continue;
//								}
//								geneVarInfo.second.deletionsFinal[del.first][d.first] = d.second;
//								if(d.second/static_cast<double>(totalInput) > lowVariantCutOff){
//									variablePositions.emplace_back(del.first);
//								}
//							}
//						}
//						for(const auto & ins : geneVarInfo.second.insertions){
//							for(const auto & i : ins.second){
//								if(i.second < occurrenceCutOff){
//									continue;
//								}
//								geneVarInfo.second.insertionsFinal[ins.first][i.first] = i.second;
//								if(i.second/static_cast<double>(totalInput) > lowVariantCutOff){
//									variablePositions.emplace_back(ins.first);
//								}
//							}
//						}
//
//						auto snpsPositions = getVectorOfMapKeys(geneVarInfo.second.snpsFinal);
//						njh::sort(snpsPositions);
//						OutputStream snpTabOut(OutOptions(njh::files::make_path(setUp.pars_.directoryName_, njh::pasteAsStr(geneVarInfo.first +  "-protein_aminoAcidChanges.tab.txt"))));
//						snpTabOut << "chrom\tposition\tref\tvariant\tcount\tfrequency\talleleDepth\tsamples" << std::endl;
//
//						for(const auto & snpPos : snpsPositions){
//							for(const auto & aa : geneVarInfo.second.snpsFinal[snpPos]){
//								snpTabOut << geneVarInfo.first
//										<< "\t" << snpPos + 1
//										<< "\t" << geneInfos[geneVarInfo.first]->protein_.seq_[snpPos]
//										<< "\t" << aa.first
//										<< "\t" << aa.second
//										<< "\t" << aa.second/static_cast<double>(totalInput)
//										<< "\t" << totalInput
//										<< "\t" << samplesCalled << std::endl;
//							}
//						}
//						if(!variablePositions.empty()){
//							uint32_t variableStart = vectorMinimum(variablePositions);
//							uint32_t variableStop = vectorMaximum(variablePositions);
//
//							GenomicRegion variableRegion;
//							variableRegion.chrom_ = geneVarInfo.first;
//							variableRegion.start_ = variableStart;
//							variableRegion.end_ = variableStop;
//
//							OutputStream bedVariableRegionOut(OutOptions(njh::files::make_path(setUp.pars_.directoryName_, njh::pasteAsStr(geneVarInfo.first +  "-protein_variableRegion.bed"))));
//							bedVariableRegionOut << variableRegion.genBedRecordCore().toDelimStrWithExtra() << std::endl;
//							//initiate denovo typer
//							std::unordered_map<std::string, std::vector<uint32_t>> denovoAmino;
//							denovoAmino[geneVarInfo.first]=variablePositions;
//							GenomicAminoAcidPositionTyper denovoAATyper(denovoAmino);
//							GenomicRegion currentChrom(*genes[geneVarInfo.first]->gene_);
//							auto curAligner = alnPool.popAligner();
//							auto curBReader = bamPool.popReader();
//							TwoBit::TwoBitFile tReader(twoBitFnp);
//							BamTools::BamAlignment bAln;
//							//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
//							setBamFileRegionThrow(*curBReader, currentChrom);
//							//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
//							while (curBReader->GetNextAlignment(bAln)) {
//								if (bAln.IsMapped() && bAln.IsPrimaryAlignment()) {
//									//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
//									auto balnGenomicRegion = GenomicRegion(bAln, refData);
//									if (setUp.pars_.verbose_) {
//										std::cout << balnGenomicRegion.genBedRecordCore().toDelimStr() << std::endl;
//									}
//									if(!njh::in(balnGenomicRegion.createUidFromCoords(), alnRegionToGeneIds)){
//										continue;
//									}
//									//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
//									for (const auto & g : alnRegionToGeneIds.at(balnGenomicRegion.createUidFromCoords())) {
//										//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
//										const auto & currentGene = genes.at(g);
//										const auto & currentGeneInfo = geneInfos.at(g);
//										//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
//										auto aminoTyping = denovoAATyper.typeAlignment(bAln, *currentGene, *currentGeneInfo, tReader, *curAligner, refData);
//										if (!aminoTyping.empty() && njh::in(g, denovoAATyper.aminoPositionsForTyping_)) {
//											//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
//											popHapsTypeddenovo[bAln.Name].emplace(g, GenomicAminoAcidPositionTyper::GeneAminoTyperInfo(g, aminoTyping));
//										}
//									}
//									//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
//								}
//							}
//						}
//					}
//				}
//		}
//
//		if(!allMetaKeys.empty()){
//			OutputStream metaOutPut(OutOptions(njh::files::make_path(setUp.pars_.directoryName_, "popHapMeta.tab.txt")));
//			VecStr columnNames = allMetaKeysVec;
//			if("" != gffFnp){
//				columnNames.emplace_back("knownProteinPositionsTyped");
//				columnNames.emplace_back("denovoProteinPositionsTyped");
//			}
//			metaOutPut << "PopName\tHapPopUIDCount" << njh::conToStr(columnNames, "\t") << std::endl;
//			for(const auto & popSeq : clusters){
//				for(const auto & inputSeq : popSeq.reads_){
//					MetaDataInName outMeta;
//					VecStr missingMetaFields;
//					if (MetaDataInName::nameHasMetaData(inputSeq->seqBase_.name_)) {
//						MetaDataInName seqMeta(inputSeq->seqBase_.name_);
//						outMeta = seqMeta;
//						for (const auto & mf : allMetaKeysVec) {
//							if (!njh::in(mf, seqMeta.meta_)) {
//								missingMetaFields.emplace_back(mf);
//							}
//						}
//					} else {
//						missingMetaFields = allMetaKeysVec;
//					}
//					for(const auto & mf : missingMetaFields){
//						outMeta.addMeta(mf, "NA");
//					}
//					metaOutPut << popSeq.seqBase_.name_;
//					MetaDataInName popSeqMeta(popSeq.seqBase_.name_);
//					metaOutPut << "\t" << popSeqMeta.getMeta("HapPopUIDCount");
//
//					for(const auto & mf : allMetaKeysVec){
//						metaOutPut << '\t' << outMeta.meta_[mf];
//					}
//
//					if(!popHapsTyped.empty()){
//						std::string type = "";
//						for(const auto & gene : popHapsTyped[popSeq.seqBase_.name_]){
//							std::string geneType = gene.first + "=";
//							for(const auto & amino : gene.second.aminos_){
//								if(' ' != amino.second){
//									if(gene.first + "=" != geneType){
//										geneType += "-";
//									}
//									//+1 to make it non zero-based amino acid positioning
//									geneType += njh::pasteAsStr(amino.first + 1, ":", amino.second);
//								}
//							}
//							if(gene.first + "=" != geneType){
//								type += geneType + ";";
//							}
//						}
//						metaOutPut << "\t" << type;
//					}else if("" != gffFnp){
//						metaOutPut << "\t";
//					}
//
//					if(!popHapsTypeddenovo.empty()){
//						std::string type = "";
//						for(const auto & gene : popHapsTypeddenovo[popSeq.seqBase_.name_]){
//							std::string geneType = gene.first + "=";
//							for(const auto & amino : gene.second.aminos_){
//								if(' ' != amino.second){
//									if(gene.first + "=" != geneType){
//										geneType += "-";
//									}
//									//+1 to make it non zero-based amino acid positioning
//									geneType += njh::pasteAsStr(amino.first + 1, ":", amino.second);
//								}
//							}
//							if(gene.first + "=" != geneType){
//								type += geneType + ";";
//							}
//						}
//						metaOutPut << "\t" << type;
//					}else if("" != gffFnp){
//						metaOutPut << "\t";
//					}
//
//					metaOutPut<< std::endl;
//				}
//			}
//		}
//		OutputStream metaLabNamesOut(OutOptions(njh::files::make_path(setUp.pars_.directoryName_, "uniqueSeqs_nonFieldSampleNames.tab.txt")));
//		metaLabNamesOut << "name\tsamples" << std::endl;
//		for(const auto & cIter : clusters){
//			auto inputNames = readVec::getNames(cIter.reads_);
//			nameOut << cIter.seqBase_.name_
//					<< "\t" << cIter.seqBase_.cnt_
//					<< "\t" << njh::conToStr(inputNames, ",") << std::endl;
//			readLens[len(cIter.seqBase_)]+= cIter.seqBase_.cnt_;
//			VecStr nonFieldSampleNames{};
//			for(const auto & name : inputNames){
//				if(MetaDataInName::nameHasMetaData(name)){
//					MetaDataInName meta(name);
//					if(meta.containsMeta("BiologicalSample") && meta.containsMeta("IsFieldSample") && !meta.getMeta<bool>("IsFieldSample")){
//						if(!meta.containsMeta("country") || "NA" != meta.getMeta("country")){
//							nonFieldSampleNames.emplace_back(meta.getMeta("BiologicalSample"));
//						}
//					}
//				}
//			}
//			metaLabNamesOut << cIter.seqBase_.name_ << "\t" << njh::conToStr(nonFieldSampleNames, ",") << std::endl;
//		}
//
//		table readLensTab(readLens, VecStr{"length", "count"});
//		OutputStream readLensOut(njh::files::make_path(setUp.pars_.directoryName_, "readLengthCounts.tab.txt"));
//		readLensTab.outPutContents(readLensOut, "\t");
//		uint32_t mostCommonReadLen = 0;
//		uint32_t mostCommonReadLenCount = 0;
//		for(const auto & rlen : readLens){
//			if(rlen.second > mostCommonReadLenCount){
//				mostCommonReadLen = rlen.first;
//				mostCommonReadLenCount = rlen.second;
//			}
//		}
//		uint32_t countOfNotMostCommonReadLen = 0;
//		for(const auto & rlen : readLens){
//			if(rlen.first != mostCommonReadLen && uAbsdiff(mostCommonReadLen , rlen.first) >= lengthDiffForLengthPoly){
//				countOfNotMostCommonReadLen += rlen.second;
//			}
//		}
//		bool lengthPoly = false;
//		if((countOfNotMostCommonReadLen/static_cast<double>(totalInput)) > lowVariantCutOff){
//			lengthPoly = true;
//		}
//		readVec::allSetFractionByTotalCount(clusters);
//		PopGenCalculator::DiversityMeasures divMeasures  = PopGenCalculator::getGeneralMeasuresOfDiversity(clusters);
//		//PairwisePairFactory pFacIfTotal(totalInput);
//		out << identifier
//				<< "\t" << totalInput
//				<< "\t" << clusters.size()
//				<< "\t" << samplesCalled
//				<< "\t" << divMeasures.singlets_
//				<< "\t" << divMeasures.doublets_
//				<< "\t" << divMeasures.expShannonEntropy_
//				<< "\t" << divMeasures.ShannonEntropyE_
//				<< "\t" << divMeasures.effectiveNumOfAlleles_
//				<< "\t" << divMeasures.heterozygostiy_
//				<< "\t" << (lengthPoly ? "true" : "false")
//				<< "\t" << numSNPs
//				<< "\t" << variableRegionID
//				<< "\t" << avgGCContent
//				//<< "\t" << sumOfPairwiseDist/pFac.totalCompares_
//				<< std::endl;
//				//<< "\t" << sumOfPairwiseDistWeighted/pFacIfTotal.totalCompares_<< std::endl;
//	}
//
//	auto variantCallsDir = njh::files::make_path(setUp.pars_.directoryName_, "variantCalls");
//	njh::files::makeDir(njh::files::MkdirPar{variantCallsDir});
//
//	if(calculatingFst){
//		OutputStream out(njh::files::make_path(setUp.pars_.directoryName_, popMeta + "_basicInfo.tab.txt"));
//		out << "id\t" << popMeta << "\ttotalHaplotypes\tuniqueHaplotypes\tnsamples\tsingletons\tdoublets\texpShannonEntropy\tShannonEntropyE\teffectiveNumOfAlleles\the\tlengthPolymorphism\tnSNPsAbove" << lowVariantCutOff * 100 << "%" << std::endl;
//		std::unordered_map<std::string, double> heForFields;
//
//		for(auto & field : uniqueSeqsByMeta){
//			double sumOfSquares = 0;
//			std::map<uint32_t, uint32_t> readLens;
//			njh::sort(*field.second);
//			uint32_t singletons = 0;
//			uint32_t doublets = 0;
//			double seqCount = readVec::getTotalReadCount(*field.second);
//			for(const auto & cIter : *field.second){
//				if(1 == cIter.seqBase_.cnt_){
//					++singletons;
//				}
//				if(2 == cIter.seqBase_.cnt_){
//					++doublets;
//				}
//				sumOfSquares += std::pow(cIter.seqBase_.cnt_/seqCount, 2);
//				readLens[len(cIter.seqBase_)]+= cIter.seqBase_.cnt_;
//			}
//			std::unordered_map<uint32_t, std::unordered_map<char,        uint32_t>> snps;
//			std::unordered_map<uint32_t, std::unordered_map<std::string, uint32_t>> insertions;
//			std::unordered_map<uint32_t, std::unordered_map<std::string, uint32_t>> deletions;
//
//			std::unordered_map<uint32_t, std::unordered_map<char,        uint32_t>> snpsFinal;
//			std::unordered_map<uint32_t, std::unordered_map<std::string, uint32_t>> insertionsFinal;
//			std::unordered_map<uint32_t, std::unordered_map<std::string, uint32_t>> deletionsFinal;
//
//			for(const auto & seq : *field.second){
//				alignerObj.alignCacheGlobal(refSeq, seq);
//				alignerObj.profileAlignment(refSeq, seq, false, false, false);
//				for(const auto & m : alignerObj.comp_.distances_.mismatches_){
//					snps[m.second.refBasePos][m.second.seqBase]+= seq.seqBase_.cnt_;
//				}
//				for(const auto & g : alignerObj.comp_.distances_.alignmentGaps_){
//					if(g.second.ref_){
//						//insertion
//						insertions[g.second.refPos_ - 1][g.second.gapedSequence_]++;
//					}else{
//						//deletion
//						deletions[g.second.refPos_ - 1][g.second.gapedSequence_]++;
//					}
//				}
//			}
//			std::vector<uint32_t> variablePositions;
//			//filter snps and indels by occurrence cut off
//			for(const auto & snp : snps){
//				for(const auto & b : snp.second){
//					if(b.second < occurrenceCutOff){
//						continue;
//					}
//					snpsFinal[snp.first][b.first] = b.second;
//					if(b.second/static_cast<double>(seqCount) > lowVariantCutOff){
//						variablePositions.emplace_back(snp.first);
//					}
//				}
//			}
//			for(const auto & del : deletions){
//				for(const auto & d : del.second){
//					if(d.second < occurrenceCutOff){
//						continue;
//					}
//					deletionsFinal[del.first][d.first] = d.second;
//					if(d.second/static_cast<double>(seqCount) > lowVariantCutOff){
//						variablePositions.emplace_back(del.first);
//					}
//				}
//			}
//			for(const auto & ins : insertions){
//				for(const auto & i : ins.second){
//					if(i.second < occurrenceCutOff){
//						continue;
//					}
//					insertionsFinal[ins.first][i.first] = i.second;
//					if(i.second/static_cast<double>(seqCount) > lowVariantCutOff){
//						variablePositions.emplace_back(ins.first);
//					}
//				}
//			}
//			uint32_t numSNPs = 0;
//			for(const auto & snp : snpsFinal){
//				uint32_t above = 0;
//				for(const auto & b : snp.second){
//					if(b.second/static_cast<double>(seqCount) > lowVariantCutOff){
//						++above;
//					}
//				}
//				if(above >= 1){
//					++numSNPs;
//				}
//			}
//			uint32_t mostCommonReadLen = 0;
//			uint32_t mostCommonReadLenCount = 0;
//			for(const auto & rlen : readLens){
//				if(rlen.second > mostCommonReadLenCount){
//					mostCommonReadLen = rlen.first;
//					mostCommonReadLenCount = rlen.second;
//				}
//			}
//			uint32_t countOfNotMostCommonReadLen = 0;
//			for(const auto & rlen : readLens){
//				if(rlen.first != mostCommonReadLen && uAbsdiff(mostCommonReadLen , rlen.first) >= lengthDiffForLengthPoly){
//					countOfNotMostCommonReadLen += rlen.second;
//				}
//			}
//
//			bool lengthPoly = false;
//			if((countOfNotMostCommonReadLen/static_cast<double>(seqCount)) > lowVariantCutOff){
//				lengthPoly = true;
//			}
//
//			double he = 1 - sumOfSquares;
//			heForFields[field.first] = he;
//			auto divMeasures = PopGenCalculator::getGeneralMeasuresOfDiversity(*field.second);
//			out << identifier
//					<< "\t" << field.first
//					<< "\t" << seqCount
//					<< "\t" << field.second->size()
//					<< "\t" << samplesByMeta.at(field.first).size()
//					<< "\t" << divMeasures.singlets_
//					<< "\t" << divMeasures.doublets_
//					<< "\t" << divMeasures.expShannonEntropy_
//					<< "\t" << divMeasures.ShannonEntropyE_
//					<< "\t" << divMeasures.effectiveNumOfAlleles_
//					<< "\t" << divMeasures.heterozygostiy_
//					<< "\t" << (lengthPoly ? "true" : "false")
//					<< "\t" << numSNPs << std::endl;
//
//			if(!snpsFinal.empty() || ! deletionsFinal.empty() || !insertionsFinal.empty()){
//				OutputStream vcfOut(njh::files::make_path(variantCallsDir, "variants_" + field.first + ".vcf"));
//				std::unordered_set<uint32_t> positionsSet;
//				for(const auto & snps : snpsFinal){
//					positionsSet.emplace(snps.first);
//				}
//				for(const auto & ins : insertionsFinal){
//					positionsSet.emplace(ins.first);
//				}
//				for(const auto & del : deletionsFinal){
//					positionsSet.emplace(del.first);
//				}
//				vcfOut << "##fileformat=VCFv4.0" << std::endl;
//				vcfOut << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Allele Depth\">" << std::endl;
//				vcfOut << "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">" << std::endl;
//				vcfOut << "##INFO=<ID=AF,Number=.,Type=Float,Description=\"Allele Frequency\">" << std::endl;
//				vcfOut << "##INFO=<ID=AC,Number=.,Type=Integer,Description=\"Allele Count\">" << std::endl;
//				vcfOut << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" << std::endl;
//				std::vector<uint32_t> positions(positionsSet.begin(), positionsSet.end());
//				uint32_t samplesCalledForField = seqCount;
//				if(!samplesByMeta[field.first].empty()){
//					samplesCalledForField = samplesByMeta[field.first].size();
//				}
//				njh::sort(positions);
//				if(refRegion.reverseSrand_){
//					njh::reverse(positions);
//				}
//				for(const auto & pos : positions){
//					vcfOut <<  refRegion.chrom_
//							<< "\t" << refRegion.getRelativePositionFromStartStrandAware(pos) + 1
//							<< "\t" << "."
//							<< "\t";
//					std::vector<std::string> alts;
//					std::vector<uint32_t> altsCounts;
//					std::vector<double> altsFreqs;
//
//					if(refRegion.reverseSrand_){
//						vcfOut << seqUtil::reverseComplement(std::string(1, refSeq.seq_[pos]), "DNA") << "\t";
//						if(njh::in(pos, snpsFinal)){
//							for(const auto & b : snpsFinal[pos]){
//								alts.emplace_back(seqUtil::reverseComplement(std::string(1, b.first), "DNA"));
//								altsCounts.emplace_back(b.second);
//								altsFreqs.emplace_back(b.second/static_cast<double>(seqCount));
//							}
//						}
//						if (njh::in(pos, insertionsFinal)) {
//							for (const auto & ins : insertionsFinal[pos]) {
//								alts.emplace_back(seqUtil::reverseComplement(ins.first, "DNA"));
//								altsCounts.emplace_back(ins.second);
//								altsFreqs.emplace_back(ins.second/static_cast<double>(seqCount));
//							}
//						}
//					}else{
//						vcfOut << refSeq.seq_[pos] << "\t";
//						if(njh::in(pos, snpsFinal)){
//							for(const auto & b : snpsFinal[pos]){
//								alts.emplace_back(std::string(1, b.first));
//								altsCounts.emplace_back(b.second);
//								altsFreqs.emplace_back(b.second/static_cast<double>(seqCount));
//							}
//						}
//						if (njh::in(pos, insertionsFinal)) {
//							for (const auto & ins : insertionsFinal[pos]) {
//								alts.emplace_back(ins.first);
//								altsCounts.emplace_back(ins.second);
//								altsFreqs.emplace_back(ins.second/static_cast<double>(seqCount));
//							}
//						}
//					}
//					vcfOut << njh::conToStr(alts, ",")
//					<< "\t40\tPASS\t";
//					vcfOut
//							<< "DP=" << seqCount << ";"
//							<< "NS=" << samplesCalledForField << ";"
//							<< "AC=" << njh::conToStr(altsCounts, ",") << ";"
//							<< "AF=" << njh::conToStr(altsFreqs, ",")
//					<< std::endl;
//					if (njh::in(pos, deletions)) {
//						for (const auto & d : deletions[pos]) {
//							vcfOut <<  refRegion.chrom_
//									<< "\t" << refRegion.getRelativePositionFromStartStrandAware(pos) + 1
//									<< "\t" << "."
//									<< "\t";
//							if(refRegion.reverseSrand_){
//								vcfOut << seqUtil::reverseComplement(std::string(1, refSeq.seq_[pos]), "DNA") << seqUtil::reverseComplement(d.first, "DNA")
//								<< "\t" << seqUtil::reverseComplement(std::string(1, refSeq.seq_[pos]), "DNA") << "\t";
//							}else{
//								vcfOut << refSeq.seq_[pos] << d.first
//								<< "\t" << refSeq.seq_[pos] << "\t";
//							}
//							vcfOut << "40\tPASS\t";
//							vcfOut
//									<< "DP=" << seqCount << ";"
//									<< "NS=" << samplesCalledForField << ";"
//									<< "AC=" << d.second << ";"
//									<< "AF=" << d.second/static_cast<double>(totalInput)
//							<< std::endl;
//						}
//					}
//				}
//			}
//		}
//		OutputStream pairwiseDiffsOut(njh::files::make_path(setUp.pars_.directoryName_, "pairwisePopDiffMeasures.tab.txt"));
//		pairwiseDiffsOut << "identifier"
//				<< "\t" << popMeta << "1"
//				<< "\t" << "popMeta" << "1_totalHaps"
//				<< "\t" << "popMeta" << "1_uniqueHaps"
//				<< "\t" << "popMeta" << "1_samples"
//				<< "\t" << popMeta << "2"
//				<< "\t" << "popMeta" << "2_totalHaps"
//				<< "\t" << "popMeta" << "2_uniqueHaps"
//				<< "\t" << "popMeta" << "2_samples"
//				<< "\t" << "HsSample"
//								<<"\t"<<"HsEst"
//								<<"\t"<<"HtSample"
//								<<"\t"<<"HtEst"
//								<<"\t"<<"Gst"
//								<<"\t"<<"GstEst"
//								<<"\t"<<"JostD"
//								<<"\t"<<"JostDEst"
//								<<"\t"<<"ChaoA"
//								<<"\t"<<"ChaoB"
//								<<"\t"<<"JostDChaoEst" << std::endl;
//
//		auto pairwiseMeasurements = PopGenCalculator::getPairwisePopDiff(uniqueSeqsByMeta);
//		for(const auto & field1 : pairwiseMeasurements){
//			for(const auto & field2 : field1.second){
//				uint32_t total1 = 0;
//				uint32_t total2 = 0;
//				for(const auto & haps1 : *uniqueSeqsByMeta.at(field1.first)){
//					total1 += haps1.seqBase_.cnt_;
//				}
//				for(const auto & haps2 : *uniqueSeqsByMeta.at(field2.first)){
//					total2 += haps2.seqBase_.cnt_;
//				}
//				pairwiseDiffsOut << identifier
//						<< "\t" << field1.first
//						<< "\t" << total1
//						<< "\t" << uniqueSeqsByMeta.at(field1.first)->size()
//						<< "\t" << samplesByMeta.at(field1.first).size()
//						<< "\t" << field2.first
//						<< "\t" << total2
//						<< "\t" << uniqueSeqsByMeta.at(field2.first)->size()
//						<< "\t" << samplesByMeta.at(field2.first).size()
//				    << "\t" << field2.second.hsSample_
//						<< "\t" << field2.second.hsEst_
//						<< "\t" << field2.second.htSample_
//						<< "\t" << field2.second.htEst_
//						<< "\t" << field2.second.gst_
//						<< "\t" << field2.second.gstEst_
//						<< "\t" << field2.second.jostD_
//						<< "\t" << field2.second.jostDEst_
//						<< "\t" << field2.second.chaoA_
//						<< "\t" << field2.second.chaoB_
//						<< "\t" << field2.second.jostDChaoEst_
//						<< std::endl;
//			}
//		}
//		OutputStream popDiffOut(njh::files::make_path(setUp.pars_.directoryName_, "popDiffMeasures.tab.txt"));
//		popDiffOut << "identifier"
//				<<"\t"<<"totalHaps"
//				<<"\t"<<"uniqueHaps"
//				<<"\t"<<"nsamples"
//				<<"\t"<<"HsSample"
//				<<"\t"<<"HsEst"
//				<<"\t"<<"HtSample"
//				<<"\t"<<"HtEst"
//				<<"\t"<<"Gst"
//				<<"\t"<<"GstEst"
//				<<"\t"<<"JostD"
//				<<"\t"<<"JostDEst"
//				<<"\t"<<"ChaoA"
//				<<"\t"<<"ChaoB"
//				<<"\t"<<"JostDChaoEst" << std::endl;
//
//		//PopDifferentiationMeasures temp;
//		auto diffMeasures = PopGenCalculator::getOverallPopDiffForSeqs(uniqueSeqsByMeta);
//
//		popDiffOut << identifier
//
//				<< "\t" << totalInput
//				<< "\t" << clusters.size()
//				<< "\t" << samplesCalled
//		    << "\t" << diffMeasures.hsSample_
//				<< "\t" << diffMeasures.hsEst_
//				<< "\t" << diffMeasures.htSample_
//				<< "\t" << diffMeasures.htEst_
//				<< "\t" << diffMeasures.gst_
//				<< "\t" << diffMeasures.gstEst_
//				<< "\t" << diffMeasures.jostD_
//				<< "\t" << diffMeasures.jostDEst_
//				<< "\t" << diffMeasures.chaoA_
//				<< "\t" << diffMeasures.chaoB_
//				<< "\t" << diffMeasures.jostDChaoEst_
//				<< std::endl;
//	}
//
//	alignerObj.processAlnInfoOutput(setUp.pars_.outAlnInfoDirName_, setUp.pars_.verbose_);
//	if(!snpsFinal.empty() || ! deletionsFinal.empty() || !insertionsFinal.empty()){
//		OutputStream snpTabOut(njh::files::make_path(variantCallsDir, "allSNPs.tab.txt"));
//		OutputStream vcfOut(njh::files::make_path(variantCallsDir, "allVariants.vcf"));
//		std::unordered_set<uint32_t> positionsSet;
//		for(const auto & snps : snpsFinal){
//			positionsSet.emplace(snps.first);
//		}
//		for(const auto & ins : insertionsFinal){
//			positionsSet.emplace(ins.first);
//		}
//		for(const auto & del : deletionsFinal){
//			positionsSet.emplace(del.first);
//		}
//		vcfOut << "##fileformat=VCFv4.0" << std::endl;
//		vcfOut << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Allele Depth\">" << std::endl;
//		vcfOut << "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">" << std::endl;
//		vcfOut << "##INFO=<ID=AF,Number=.,Type=Float,Description=\"Allele Frequency\">" << std::endl;
//		vcfOut << "##INFO=<ID=AC,Number=.,Type=Integer,Description=\"Allele Count\">" << std::endl;
//		vcfOut << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" << std::endl;
//		snpTabOut << "chrom\tposition\tref\tvariant\tcount\tfrequency\talleleDepth\tsamples" << std::endl;
//
//		std::vector<uint32_t> positions(positionsSet.begin(), positionsSet.end());
//
//		njh::sort(positions);
//		if(refRegion.reverseSrand_){
//			njh::reverse(positions);
//		}
//		for(const auto & pos : positions){
//			vcfOut <<  refRegion.chrom_
//					<< "\t" << refRegion.getRelativePositionFromStartStrandAware(pos) + 1
//					<< "\t" << "."
//					<< "\t";
//			std::vector<std::string> alts;
//			std::vector<uint32_t> altsCounts;
//			std::vector<double> altsFreqs;
//
//			if(refRegion.reverseSrand_){
//				vcfOut << seqUtil::reverseComplement(std::string(1, refSeq.seq_[pos]), "DNA") << "\t";
//				if(njh::in(pos, snpsFinal)){
//					uint32_t snpCount = 0;
//					for(const auto & b : snpsFinal[pos]){
//						snpTabOut << refRegion.chrom_
//								<< "\t" << refRegion.getRelativePositionFromStartStrandAware(pos)
//								<< "\t" << seqUtil::reverseComplement(std::string(1, refSeq.seq_[pos]), "DNA")
//								<< "\t" << seqUtil::reverseComplement(std::string(1, b.first), "DNA")
//								<< "\t" << b.second
//								<< "\t" << b.second/static_cast<double>(totalInput)
//								<< "\t" << totalInput
//								<< "\t" << samplesCalled << std::endl;
//						snpCount+= b.second;
//						alts.emplace_back(seqUtil::reverseComplement(std::string(1, b.first), "DNA"));
//						altsCounts.emplace_back(b.second);
//						altsFreqs.emplace_back(b.second/static_cast<double>(totalInput));
//					}
//					snpTabOut << refRegion.chrom_
//							<< "\t" << refRegion.getRelativePositionFromStartStrandAware(pos)
//							<< "\t" << seqUtil::reverseComplement(std::string(1, refSeq.seq_[pos]), "DNA")
//							<< "\t" << seqUtil::reverseComplement(std::string(1, refSeq.seq_[pos]), "DNA")
//							<< "\t" << totalInput - snpCount
//							<< "\t" << (totalInput - snpCount)/static_cast<double>(totalInput)
//							<< "\t" << totalInput
//							<< "\t" << samplesCalled << std::endl;
//				}
//				if (njh::in(pos, insertionsFinal)) {
//					for (const auto & ins : insertionsFinal[pos]) {
//						alts.emplace_back(seqUtil::reverseComplement(ins.first, "DNA"));
//						altsCounts.emplace_back(ins.second);
//						altsFreqs.emplace_back(ins.second/static_cast<double>(totalInput));
//					}
//				}
//			}else{
//				vcfOut << refSeq.seq_[pos] << "\t";
//				if(njh::in(pos, snpsFinal)){
//					uint32_t snpCount = 0;
//					for(const auto & b : snpsFinal[pos]){
//						snpTabOut << refRegion.chrom_
//								<< "\t" << refRegion.getRelativePositionFromStartStrandAware(pos)
//								<< "\t" << std::string(1, refSeq.seq_[pos])
//								<< "\t" << std::string(1, b.first)
//								<< "\t" << b.second
//								<< "\t" << b.second/static_cast<double>(totalInput)
//								<< "\t" << totalInput
//								<< "\t" << samplesCalled << std::endl;
//						snpCount+= b.second;
//						alts.emplace_back(std::string(1, b.first));
//						altsCounts.emplace_back(b.second);
//						altsFreqs.emplace_back(b.second/static_cast<double>(totalInput));
//					}
//					snpTabOut << refRegion.chrom_
//							<< "\t" << refRegion.getRelativePositionFromStartStrandAware(pos)
//							<< "\t" << std::string(1, refSeq.seq_[pos])
//							<< "\t" << std::string(1, refSeq.seq_[pos])
//							<< "\t" << totalInput - snpCount
//							<< "\t" << (totalInput - snpCount)/static_cast<double>(totalInput)
//							<< "\t" << totalInput
//							<< "\t" << samplesCalled << std::endl;
//				}
//				if (njh::in(pos, insertionsFinal)) {
//					for (const auto & ins : insertionsFinal[pos]) {
//						alts.emplace_back(ins.first);
//						altsCounts.emplace_back(ins.second);
//						altsFreqs.emplace_back(ins.second/static_cast<double>(totalInput));
//					}
//				}
//			}
//			vcfOut << njh::conToStr(alts, ",")
//			<< "\t40\tPASS\t";
//			vcfOut
//					<< "DP=" << totalInput << ";"
//					<< "NS=" << samplesCalled << ";"
//					<< "AC=" << njh::conToStr(altsCounts, ",") << ";"
//					<< "AF=" << njh::conToStr(altsFreqs, ",")
//			<< std::endl;
//
//			if (njh::in(pos, deletions)) {
//				for (const auto & d : deletions[pos]) {
//					vcfOut <<  refRegion.chrom_
//							<< "\t" << refRegion.getRelativePositionFromStartStrandAware(pos) + 1
//							<< "\t" << "."
//							<< "\t";
//					if(refRegion.reverseSrand_){
//						vcfOut << seqUtil::reverseComplement(std::string(1, refSeq.seq_[pos]), "DNA") << seqUtil::reverseComplement(d.first, "DNA")
//						<< "\t" << seqUtil::reverseComplement(std::string(1, refSeq.seq_[pos]), "DNA") << "\t";
//					}else{
//						vcfOut << refSeq.seq_[pos] << d.first
//						<< "\t" << refSeq.seq_[pos] << "\t";
//					}
//					vcfOut << "40\tPASS\t";
//					vcfOut
//							<< "DP=" << totalInput << ";"
//							<< "NS=" << samplesCalled << ";"
//							<< "AC=" << d.second << ";"
//							<< "AF=" << d.second/static_cast<double>(totalInput)
//					<< std::endl;
//				}
//			}
//		}
//	}
//	return 0;
	return 0;
}



}  //namespace njhseq



