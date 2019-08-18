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
#include <njhseq/objects/Gene/TranslatorByAlignment.hpp>


namespace njhseq {









int genExpRunner::translateSeqsBasedOnGFF(const njh::progutils::CmdArgs & inputCommands){

  TranslatorByAlignment::RunPars variantCallerRunPars;

	TranslatorByAlignment::TranslatorByAlignmentPars transPars;

	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(transPars.keepTemporaryFiles_, "--keepTemporaryFiles", "Keep Temporary Files");
	setUp.setOption(transPars.useLastz_, "--useLastz", "Use lastz instead of bowtie2 to align sequences to the genome");
	setUp.setOption(transPars.gffFnp_, "--gffFnp", "GFF", true);
	setUp.setOption(transPars.lzPars_.genomeFnp, "--genomeFnp", "Genome fasta", true);
	setUp.setOption(variantCallerRunPars.realnPars.extendAmount, "--extendAmount", "Amount to extend for re-alignment");
	setUp.processReadInNames();
	setUp.processDirectoryOutputName(true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);

	transPars.workingDirtory_ = setUp.pars_.directoryName_;
	TranslatorByAlignment translator(transPars);
	SeqInput seqReader(setUp.pars_.ioOptions_);
	auto input = seqReader.readAllReads<seqInfo>();
	std::unordered_map<std::string, uint32_t> counts;
	for(const auto & seq : input ){
		counts[seq.name_] = std::ceil(seq.cnt_);
	}
	auto results = translator.run(setUp.pars_.ioOptions_, counts, variantCallerRunPars);
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


