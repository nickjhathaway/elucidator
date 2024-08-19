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
#include <njhseq/PopulationGenetics.h>
#include <njhseq/objects/Gene/GenomicAminoAcidPositionTyper.hpp>
#include <njhseq/GenomeUtils/GenomeExtraction/ParsingAlignmentInfo.h>
#include <njhseq/BamToolsUtils.h>
#include <njhseq/objects/Gene/TranslatorByAlignment.hpp>


#include <cppitertools/sorted.hpp>

namespace njhseq {









int genExpRunner::translateSeqsBasedOnGFF(const njh::progutils::CmdArgs & inputCommands){

  TranslatorByAlignment::RunPars variantCallerRunPars;
	variantCallerRunPars.occurrenceCutOff = 0;
	TranslatorByAlignment::TranslatorByAlignmentPars transPars;

	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(transPars.keepTemporaryFiles_, "--keepTemporaryFiles", "Keep Temporary Files");
	setUp.setOption(transPars.useLastz_, "--useLastz", "Use lastz instead of bowtie2 to align sequences to the genome");
	setUp.setOption(transPars.gffFnp_, "--gffFnp", "GFF", true);
	setUp.setOption(transPars.lzPars_.genomeFnp, "--genomeFnp", "Genome fasta", true);
  variantCallerRunPars.lowVariantCutOff = 0;
  variantCallerRunPars.occurrenceCutOff = 0;
	setUp.setOption(variantCallerRunPars.realnPars.extendAmount, "--extendAmount", "Amount to extend for re-alignment");
	setUp.processReadInNames();
	setUp.processDirectoryOutputName(true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);

	transPars.workingDirtory_ = setUp.pars_.directoryName_;
	TranslatorByAlignment translator(transPars);
	//key1 = hap, key2 = sample, value = read counts
  std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>> counts;

  SeqInput seqReader(setUp.pars_.ioOptions_);
  auto input = seqReader.readAllReads<seqInfo>(); {
	  uint32_t ongoingCount = 0;
	  for (const auto& seq: input) {
		  auto seqCount = static_cast<uint32_t>(std::ceil(seq.cnt_));
		  if (MetaDataInName::nameHasMetaData(seq.name_)) {
			  MetaDataInName seqMeta(seq.name_);
			  if (seqMeta.containsMeta("sample")) {
				  counts[seq.name_][seqMeta.getMeta("sample")] = seqCount;
			  }
		  } else {
			  counts[seq.name_][njh::pasteAsStr("samp.", ongoingCount)] = seqCount;
		  }
		  ++ongoingCount;
	  }
  }
	auto translatedRes = translator.run(setUp.pars_.ioOptions_, counts, variantCallerRunPars);

	SeqOutput writer(SeqIOOptions::genFastaOut(njh::files::make_path(setUp.pars_.directoryName_, "translatedInput.fasta")));
	SeqOutput writerAln(SeqIOOptions::genFastaOut(njh::files::make_path(setUp.pars_.directoryName_, "aln_translatedInput.fasta")));
	SeqOutput writerCDNA(SeqIOOptions::genFastaOut(njh::files::make_path(setUp.pars_.directoryName_, "cdna_input.fasta")));
	OutputStream outInfo(OutOptions(njh::files::make_path(setUp.pars_.directoryName_, "aa_positionCoveredInfo.tab.txt")));
	OutputStream outChangesInfo(OutOptions(njh::files::make_path(setUp.pars_.directoryName_, "aa_changesInfo.tab.txt")));
	outInfo << "seqName\ttranscript\taminoAcidStart\taminoAcidStop" << std::endl;
	outChangesInfo << "seqName\ttranscript\taaPosition(1-based)\ttype\trefAA\tseqAA\tUID" << std::endl;

	for (const auto & seqName : iter::sorted(getVectorOfMapKeys(translatedRes.translations_))) {
		for (const auto & transcript : translatedRes.translations_.at(seqName)) {

			outInfo << seqName
					<< "\t" << transcript.first
					<< "\t" << std::get<0>(transcript.second.firstAminoInfo_).aaPos_
					<< "\t" << std::get<0>(transcript.second.lastAminoInfo_).aaPos_  << std::endl;
			writer.openWrite(transcript.second.translation_);
			writerCDNA.openWrite(transcript.second.cDna_);
			if(setUp.pars_.debug_){
				writerAln.openWrite(transcript.second.refAlnTranslation_);
				writerAln.openWrite(transcript.second.queryAlnTranslation_);
			}
			for(const auto & aaMismatch : transcript.second.comp_.distances_.mismatches_){
				outChangesInfo<< seqName
						<< "\t" << transcript.first
						<< "\t" << aaMismatch.second.refBasePos + 1
						<< "\t" << "mismatch"
						<< "\t" << aaMismatch.second.seqBase
						<< "\t" << aaMismatch.second.refBase
						<< "\t" << aaMismatch.second.refBase << aaMismatch.second.refBasePos +1<< aaMismatch.second.seqBase
						<< std::endl;
			}

			for(const auto & aaIndel : transcript.second.comp_.distances_.alignmentGaps_){
				outChangesInfo<< seqName
						<< "\t" << transcript.first
						<< "\t" << aaIndel.second.refPos_ + 1
						<< "\t" << (aaIndel.second.ref_ ? "deletion": "insertion")
						<< "\t" << (aaIndel.second.ref_ ? std::string(""): aaIndel.second.gapedSequence_)
						<< "\t" << (aaIndel.second.ref_ ? aaIndel.second.gapedSequence_: std::string(""))
						<< "\t" << ""
						<< std::endl;
			}
		}
	}

  

  OutputStream popBedLocs(njh::files::make_path(setUp.pars_.directoryName_, "inputSeqs.bed"));
  translatedRes.writeOutSeqAlnIndvVars(njh::files::make_path(setUp.pars_.directoryName_, "variantsPerSeqAln.tab.txt.gz"));
  translatedRes.writeSeqLocations(popBedLocs);

  std::unordered_map<std::string, std::set<uint32_t>> knownAAMutsChromPositions;

  OutputStream seqsUnableToBeMappedOut(njh::files::make_path(setUp.pars_.directoryName_, "seqsUnableToBeMapped.txt"));
  seqsUnableToBeMappedOut << njh::conToStr(translatedRes.seqsUnableToBeMapped_, "\n") << std::endl;
  OutputStream seqsTranslationFilteredOut(njh::files::make_path(setUp.pars_.directoryName_, "seqsTranslationFiltered.txt"));
  seqsTranslationFilteredOut << njh::conToStr(translatedRes.seqsTranslationFiltered_, "\n") << std::endl;

//  std::cout << "translatedRes.translations_.size(): " << translatedRes.translations_.size() << std::endl;
  if(!translatedRes.translations_.empty()) {
    SeqOutput transwriter(
        SeqIOOptions::genFastaOutGz(njh::files::make_path(setUp.pars_.directoryName_, "translatedInput.fasta.gz")));
    auto seqNames = njh::getVecOfMapKeys(translatedRes.translations_);
    njh::sort(seqNames);
    for(const auto & seqName : seqNames) {
      for (const auto &transcript: translatedRes.translations_.at(seqName)) {
        transwriter.openWrite(transcript.second.translation_);
      }
    }
  }

	OutputStream variantCalls(njh::files::make_path(setUp.pars_.directoryName_, "varCalls.vcf"));

	for(const auto & seqVar : translatedRes.seqVariants_){
		seqVar.second.writeVCF(variantCalls);
	}

	return 0;
}

} /* namespace njhseq */


