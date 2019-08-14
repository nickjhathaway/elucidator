/*
 * genExp_translateSeqsBasedOnGFF.cpp
 *
 *  Created on: Aug 13, 2019
 *      Author: nicholashathaway
 */




#include "genExp.hpp"
#include "elucidator/objects/BioDataObject.h"


#include "elucidator/seqToolsUtils/seqToolsUtils.hpp"

#include <njhseq/GenomeUtils.h>
#include <njhseq/objects/Gene/GeneFromGffs.hpp>

#include "elucidator/BioRecordsUtils/BedUtility.hpp"
#include "elucidator/PopulationGenetics.h"
#include <njhseq/objects/Gene/GenomicAminoAcidPositionTyper.hpp>
#include <njhseq/GenomeUtils/GenomeExtraction/ParsingAlignmentInfo.h>
#include <njhseq/BamToolsUtils.h>


namespace njhseq {





class TranslatorByAlignment{
public:

	struct TranslatorByAlignmentPars{
		TranslatorByAlignmentPars(){
			lzPars_.coverage = 100;
			lzPars_.identity = 70;
		}
		bfs::path gffFnp_ = "";
		BioCmdsUtils::LastZPars lzPars_;
		bool useLastz_ = false;

		bfs::path workingDirtory_;

		bool keepTemporaryFiles_ = false;
	};

	struct TranslateSeqRes {
		seqInfo translation_;
		seqInfo queryAlnTranslation_;
		seqInfo refAlnTranslation_;

		comparison comp_;
	};

	struct TranslatorByAlignmentResult{
		std::set<std::string> geneIds_;
		std::unordered_map<std::string, std::unordered_map<std::string, TranslateSeqRes>> translations_;
	};






	static std::unordered_map<std::string, TranslateSeqRes> translateBasedOnAlignment(
			const BamTools::BamAlignment & bAln,
			const GeneFromGffs & currentGene,
			const std::unordered_map<std::string, std::shared_ptr<GeneSeqInfo>> & transcriptInfosForGene,
			TwoBit::TwoBitFile & tReader,
			aligner & alignerObj,
			const BamTools::RefVector & refData){

		std::unordered_map<std::string, TranslateSeqRes> ret;

		auto results = std::make_shared<AlignmentResults>(bAln, refData, true);
		results->setRefSeq(tReader);
		results->setComparison(true);


		for(const auto & transcript : currentGene.mRNAs_){
			auto currentTranscriptInfo = transcriptInfosForGene.at(transcript->getIDAttr());
			auto genePosInfoByGDna = currentTranscriptInfo->getInfosByGDNAPos();
			bool endsAtStopCodon = false;
			uint32_t transStart = 0;
			seqInfo balnSeq(bAln.Name);
			std::vector<uint32_t> codons;
			std::vector<GFFCore> cDNAIntersectedWith;
			for (const auto & cDna : currentGene.CDS_.at(transcript->getIDAttr())) {
				if (results->gRegion_.overlaps(*cDna)) {
					cDNAIntersectedWith.emplace_back(*cDna);
				}
			}
			if(cDNAIntersectedWith.size() == 0){

			} else {
				if (cDNAIntersectedWith.size() == 1
						&& results->gRegion_.start_ >= cDNAIntersectedWith.front().start_ - 1
						&& results->gRegion_.end_ <= cDNAIntersectedWith.front().end_) {
					balnSeq = *(results->alnSeq_);
					if (currentGene.gene_->isReverseStrand()) {
						if (genePosInfoByGDna.at(results->gRegion_.start_).cDNAPos_
								== currentTranscriptInfo->cDna_.seq_.size() - 1) {
							endsAtStopCodon = true;
						}
						uint32_t gPos = results->gRegion_.end_ - 1;
						auto codon = genePosInfoByGDna.at(gPos).codonPos_;
						while (0 != codon) {
							--gPos;
							codon = genePosInfoByGDna.at(gPos).codonPos_;
							++transStart;
						}
					} else {
						if (genePosInfoByGDna.at(results->gRegion_.end_ - 1).cDNAPos_
								== currentTranscriptInfo->cDna_.seq_.size() - 1) {
							endsAtStopCodon = true;
						}
						uint32_t gPos = results->gRegion_.start_;
						uint32_t codon = genePosInfoByGDna.at(gPos).codonPos_;
						while (0 != codon) {
							++gPos;
							codon = genePosInfoByGDna.at(gPos).codonPos_;
							++transStart;
						}
					}
				} else {
					njh::sort(cDNAIntersectedWith,
							[](const GenomicRegion & reg1, const GenomicRegion & reg2) {
								if(reg1.start_ < reg2.start_) {
									return true;
								}
								return false;
							});

					if (currentGene.gene_->isReverseStrand()) {
						auto cDnaStop = cDNAIntersectedWith.back().end_;
						uint32_t gPos = std::min(cDnaStop, results->gRegion_.end_) - 1;
						auto codon = genePosInfoByGDna.at(gPos).codonPos_;
						while (0 != codon) {
							--gPos;
							codon = genePosInfoByGDna.at(gPos).codonPos_;
							++transStart;
						}
					} else {
						auto cDnaStart = cDNAIntersectedWith.front().start_ - 1;
						uint32_t gPos = std::max(cDnaStart, results->gRegion_.start_);
						uint32_t codon = genePosInfoByGDna.at(gPos).codonPos_;
						while (0 != codon) {
							++gPos;
							codon = genePosInfoByGDna.at(gPos).codonPos_;
							++transStart;
						}
					}
					std::vector<uint32_t> starts;
					std::vector<uint32_t> ends;
					for (const auto & cDna : cDNAIntersectedWith) {
						auto cDnaStart = cDna.start_ - 1;
						auto detStart = std::max(cDnaStart, results->gRegion_.start_);
						auto detStop = std::min(cDna.end_, results->gRegion_.end_);
						ends.emplace_back(detStop);
						starts.emplace_back(detStart);
						detStart -= results->gRegion_.start_;
						detStop -= results->gRegion_.start_;
						auto alnStart = getAlnPosForRealPos(results->refSeqAligned_->seq_,
								detStart);
						auto alnStop = getAlnPosForRealPos(results->refSeqAligned_->seq_,
								detStop - 1);
						balnSeq.append(
								results->alnSeqAligned_->getSubRead(alnStart,
										alnStop - alnStart + 1));
					}
					uint32_t cDnaStart = *std::min_element(starts.begin(), starts.end());
					uint32_t cDnaStop = *std::max_element(ends.begin(), ends.end());
					if (currentGene.gene_->isReverseStrand()) {
						if (genePosInfoByGDna.at(cDnaStart).cDNAPos_
								== currentTranscriptInfo->cDna_.seq_.size() - 1) {
							endsAtStopCodon = true;
						}
					} else {
						if (genePosInfoByGDna.at(cDnaStop - 1).cDNAPos_
								== currentTranscriptInfo->cDna_.seq_.size() - 1) {
							endsAtStopCodon = true;
						}
					}
					balnSeq.removeGaps();
				}
				if (currentGene.gene_->isReverseStrand()) {
					balnSeq.reverseComplementRead(false, true);
				}

				auto balnSeqTrans = balnSeq.translateRet(false, false, transStart);
				MetaDataInName transMeta;
				transMeta.addMeta("transcript", transcript->getIDAttr());
				balnSeqTrans.name_ += transMeta.createMetaName();

				alignerObj.alignCacheGlobal(currentTranscriptInfo->protein_, balnSeqTrans);
				alignerObj.profilePrimerAlignment(currentTranscriptInfo->protein_, balnSeqTrans);
				TranslateSeqRes tRes;

				tRes.translation_ = balnSeqTrans;
				tRes.refAlnTranslation_ = alignerObj.alignObjectA_.seqBase_;
				tRes.queryAlnTranslation_ = alignerObj.alignObjectB_.seqBase_;
				tRes.comp_ = alignerObj.comp_;

				ret[transcript->getIDAttr()] = tRes;
			}
		}
		return ret;
	}


	TranslatorByAlignmentPars pars_;

	TranslatorByAlignment(const TranslatorByAlignmentPars & pars): pars_(pars){
		njh::sys::requireExternalProgramThrow("samtools");
		if(!pars_.useLastz_){
			njh::sys::requireExternalProgramThrow("bowtie2");
		}else{
			njh::sys::requireExternalProgramThrow("lastz");
		}
	}

	TranslatorByAlignmentResult run(const SeqIOOptions & seqOpts){
		TranslatorByAlignmentResult ret;
		std::vector<bfs::path> fnpsToRemove;
		auto seqInputFnp = njh::files::make_path(pars_.workingDirtory_, "inputSeqs.fasta");
		fnpsToRemove.emplace_back(seqInputFnp);
		{
			//write out fasta file of input
			seqInfo seq;
			SeqInput reader(seqOpts);
			reader.openIn();
			SeqOutput writer(SeqIOOptions::genFastaOut(seqInputFnp));
			writer.openOut();
			while(reader.readNextRead(seq)){
				writer.write(seq);
			}
		}

		BioCmdsUtils bRunner(false);
		bRunner.RunFaToTwoBit(pars_.lzPars_.genomeFnp);
		if(!pars_.useLastz_){
			bRunner.RunBowtie2Index(pars_.lzPars_.genomeFnp);
		}

		auto gprefix = bfs::path(pars_.lzPars_.genomeFnp).replace_extension("");
		auto twoBitFnp = gprefix.string() + ".2bit";

		TwoBit::TwoBitFile tReader(twoBitFnp);

		auto uniqueSeqInOpts = SeqIOOptions::genFastaIn(seqInputFnp);
		uniqueSeqInOpts.out_.outFilename_ = njh::files::make_path(pars_.workingDirtory_, "aligned_inputSeqs.sorted.bam");
		uniqueSeqInOpts.out_.outExtention_ = ".sorted.bam";
		fnpsToRemove.emplace_back(uniqueSeqInOpts.out_.outFilename_);
		fnpsToRemove.emplace_back(uniqueSeqInOpts.out_.outFilename_.string() + ".bai");

		uniqueSeqInOpts.out_.transferOverwriteOpts(seqOpts.out_);
		if(!pars_.useLastz_){
			auto bowtieRunOut = bRunner.bowtie2Align(uniqueSeqInOpts, pars_.lzPars_.genomeFnp, "-D 20 -R 3 -N 1 -L 15 -i S,1,0.5 --end-to-end");
			BioCmdsUtils::checkRunOutThrow(bowtieRunOut, __PRETTY_FUNCTION__);
		}else{
			auto lastzRunOut = bRunner.lastzAlign(uniqueSeqInOpts, pars_.lzPars_);
			BioCmdsUtils::checkRunOutThrow(lastzRunOut, __PRETTY_FUNCTION__);
		}

		auto regionsCounter = GenomicRegionCounter::countRegionsInBam(uniqueSeqInOpts.out_.outName());
		auto ids = regionsCounter.getIntersectingGffIds(pars_.gffFnp_);
		ret.geneIds_ = ids;

		// get gene information
		auto geneInfoDir = njh::files::make_path(pars_.workingDirtory_, "geneInfos");
		if(pars_.keepTemporaryFiles_){
			njh::files::makeDir(njh::files::MkdirPar{geneInfoDir});
		}
		OutOptions outOpts(njh::files::make_path(geneInfoDir, "gene"));

		std::unordered_map<std::string, VecStr> idToTranscriptName;
		std::unordered_map<std::string, std::shared_ptr<GeneFromGffs>> genes = GeneFromGffs::getGenesFromGffForIds(pars_.gffFnp_, ids);
		std::unordered_map<std::string, std::unordered_map<std::string, std::shared_ptr<GeneSeqInfo>>> geneTranscriptInfos;

		uint64_t proteinMaxLen = 0;
		std::unordered_map<std::string, std::unordered_map<std::string, std::shared_ptr<GeneFromGffs>>> genesByChrom;

		for(const auto & gene : genes){
			genesByChrom[gene.second->gene_->seqid_].emplace(gene.first, gene.second);
			for(const auto & transcript : gene.second->mRNAs_){
				idToTranscriptName[gene.second->gene_->getIDAttr()].emplace_back(transcript->getIDAttr());
			}
			geneTranscriptInfos[gene.first] = gene.second->generateGeneSeqInfo(tReader, false);
			for(const auto & transcriptInfo : geneTranscriptInfos[gene.first]){
				readVec::getMaxLength(transcriptInfo.second->protein_, proteinMaxLen);
			}
			if(pars_.keepTemporaryFiles_){
				gene.second->writeOutGeneInfo(tReader, outOpts);
			}
		}

		aligner alignObj(proteinMaxLen, gapScoringParameters(5,1,0,0,0,0), substituteMatrix(1,-1));

		std::unordered_map<std::string, std::unordered_map<std::string, std::set<std::string>>> regionsToGeneIds;
		//targetName, GeneID, AA Position
		std::unordered_map<std::string, std::vector<std::string>> alnRegionToGeneIds;
		for(const auto & gCount : regionsCounter.counts_){
			for (const auto & g : genesByChrom[gCount.second.region_.chrom_]) {
				for(const auto & t : g.second->mRNAs_){
					regionsToGeneIds[gCount.second.region_.createUidFromCoords()][g.first].emplace(t->getIDAttr());
					alnRegionToGeneIds[gCount.second.region_.createUidFromCoords()].emplace_back(g.first);
				}
			}
		}
		BamTools::BamReader bReader;
		bReader.Open(uniqueSeqInOpts.out_.outName().string());
		checkBamOpenThrow(bReader, uniqueSeqInOpts.out_.outName());
		auto refData = bReader.GetReferenceData();
		BamTools::BamAlignment bAln;
		while (bReader.GetNextAlignment(bAln)) {
			if (bAln.IsMapped() && bAln.IsPrimaryAlignment()) {
				auto balnGenomicRegion = GenomicRegion(bAln, refData);
				if (!njh::in(balnGenomicRegion.createUidFromCoords(), alnRegionToGeneIds)) {
					continue;
				}
				for (const auto & g : alnRegionToGeneIds.at(balnGenomicRegion.createUidFromCoords())) {
					const auto & currentGene = genes.at(g);
					const auto & currentGeneInfo = geneTranscriptInfos.at(g);
					auto translations = translateBasedOnAlignment(bAln, *currentGene, currentGeneInfo, tReader, alignObj, refData);
					for(const auto & trans : translations){
						ret.translations_[bAln.Name].emplace(trans);
					}
				}
			}
		}
		if(!pars_.keepTemporaryFiles_){
			for(const auto & fnp : fnpsToRemove){
				if(bfs::exists(fnp)){
					bfs::remove(fnp);
				}
			}
		}
		return ret;
	}

};



int genExpRunner::translateSeqsBasedOnGFF(const njh::progutils::CmdArgs & inputCommands){


	TranslatorByAlignment::TranslatorByAlignmentPars transPars;

	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(transPars.keepTemporaryFiles_, "--keepTemporaryFiles", "Keep Temporary Files");
	setUp.setOption(transPars.useLastz_, "--useLastz", "Use lastz instead of bowtie2 to align sequences to the genome");
	setUp.setOption(transPars.gffFnp_, "--gffFnp", "GFF", true);
	setUp.setOption(transPars.lzPars_.genomeFnp, "--genomeFnp", "Genome fasta", true);
	setUp.processReadInNames();
	setUp.processDirectoryOutputName(true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);

	transPars.workingDirtory_ = setUp.pars_.directoryName_;
	TranslatorByAlignment translator(transPars);
	auto results = translator.run(setUp.pars_.ioOptions_);
	SeqOutput writer(SeqIOOptions::genFastaOut(njh::files::make_path(setUp.pars_.directoryName_, "translatedInput.fasta")));
	SeqOutput writerAln(SeqIOOptions::genFastaOut(njh::files::make_path(setUp.pars_.directoryName_, "aln_translatedInput.fasta")));

	for(const auto & seqName : results.translations_){
		for(const auto & transcript : seqName.second){
			writer.openWrite(transcript.second.translation_);
			if(setUp.pars_.debug_){
				writerAln.openWrite(transcript.second.refAlnTranslation_);
				writerAln.openWrite(transcript.second.queryAlnTranslation_);
			}
		}
	}
	return 0;
}

} /* namespace njhseq */


