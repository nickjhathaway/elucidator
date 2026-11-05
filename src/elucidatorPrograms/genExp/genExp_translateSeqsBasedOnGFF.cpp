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









namespace {

/**
 * \brief the read counts for the input seqs, samples come from the sample meta in the name if it has it
 * \return key1 = hap, key2 = sample, value = read counts
 */
std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>> genHapSampleCounts(const std::vector<seqInfo> & seqs) {
	std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>> counts;
	uint32_t ongoingCount = 0;
	for (const auto & seq : seqs) {
		const auto seqCount = static_cast<uint32_t>(std::ceil(seq.cnt_));
		std::string sample = njh::pasteAsStr("samp.", ongoingCount);
		if (MetaDataInName::nameHasMetaData(seq.name_)) {
			MetaDataInName seqMeta(seq.name_);
			if (seqMeta.containsMeta("sample")) {
				sample = seqMeta.getMeta("sample");
			}
		}
		counts[seq.name_][sample] = seqCount;
		++ongoingCount;
	}
	return counts;
}

/**
 * \brief writes the translations, cds, covered amino acids and amino acid changes for translation results
 */
class TranslatedSeqsWriter {
public:
	TranslatedSeqsWriter(const bfs::path & directory, const bool writeAlns, const bool addGenomeColumn) :
			writer_(SeqIOOptions::genFastaOut(njh::files::make_path(directory, "translatedInput.fasta"))),
			writerAln_(SeqIOOptions::genFastaOut(njh::files::make_path(directory, "aln_translatedInput.fasta"))),
			writerCDNA_(SeqIOOptions::genFastaOut(njh::files::make_path(directory, "cdna_input.fasta"))),
			outInfo_(OutOptions(njh::files::make_path(directory, "aa_positionCoveredInfo.tab.txt"))),
			outChangesInfo_(OutOptions(njh::files::make_path(directory, "aa_changesInfo.tab.txt"))),
			writeAlns_(writeAlns),
			addGenomeColumn_(addGenomeColumn) {
		const std::string genomeHeader = addGenomeColumn_ ? "\tgenome" : "";
		outInfo_ << "seqName\ttranscript\taminoAcidStart\taminoAcidStop" << genomeHeader << std::endl;
		outChangesInfo_ << "seqName\ttranscript\taaPosition(1-based)\ttype\trefAA\tseqAA\tUID" << genomeHeader << std::endl;
	}

	void write(const TranslatorByAlignment::TranslatorByAlignmentResult & translatedRes, const std::string & genome = "") {
		const std::string genomeColumn = addGenomeColumn_ ? "\t" + genome : "";
		for (const auto & seqName : iter::sorted(getVectorOfMapKeys(translatedRes.translations_))) {
			for (const auto & transcript : translatedRes.translations_.at(seqName)) {
				outInfo_ << seqName
						<< "\t" << transcript.first
						<< "\t" << std::get<0>(transcript.second.firstAminoInfo_).aaPos_
						<< "\t" << std::get<0>(transcript.second.lastAminoInfo_).aaPos_
						<< genomeColumn << std::endl;
				writer_.openWrite(transcript.second.translation_);
				writerCDNA_.openWrite(transcript.second.cDna_);
				if (writeAlns_) {
					writerAln_.openWrite(transcript.second.refAlnTranslation_);
					writerAln_.openWrite(transcript.second.queryAlnTranslation_);
				}
				for (const auto & aaMismatch : transcript.second.comp_.distances_.mismatches_) {
					outChangesInfo_ << seqName
							<< "\t" << transcript.first
							<< "\t" << aaMismatch.second.refBasePos + 1
							<< "\t" << "mismatch"
							<< "\t" << aaMismatch.second.refBase
							<< "\t" << aaMismatch.second.seqBase
							<< "\t" << aaMismatch.second.refBase << aaMismatch.second.refBasePos + 1 << aaMismatch.second.seqBase
							<< genomeColumn << std::endl;
				}
				for (const auto & aaIndel : transcript.second.comp_.distances_.alignmentGaps_) {
					outChangesInfo_ << seqName
							<< "\t" << transcript.first
							<< "\t" << aaIndel.second.refPos_ + 1
							<< "\t" << (aaIndel.second.ref_ ? "insertion" : "deletion")
							<< "\t" << (aaIndel.second.ref_ ? std::string("") : aaIndel.second.gapedSequence_)
							<< "\t" << (aaIndel.second.ref_ ? aaIndel.second.gapedSequence_ : std::string(""))
							<< "\t" << ""
							<< genomeColumn << std::endl;
				}
			}
		}
	}

private:
	SeqOutput writer_;
	SeqOutput writerAln_;
	SeqOutput writerCDNA_;
	OutputStream outInfo_;
	OutputStream outChangesInfo_;
	bool writeAlns_;
	bool addGenomeColumn_;
};

/**
 * \brief writes the locations, per alignment variants, unmapped and untranslatable seqs and variant calls for one translation run
 */
void writeTranslationRunInfo(const bfs::path & directory, TranslatorByAlignment::TranslatorByAlignmentResult & translatedRes) {
	OutputStream popBedLocs(njh::files::make_path(directory, "inputSeqs.bed"));
	translatedRes.writeOutSeqAlnIndvVars(njh::files::make_path(directory, "variantsPerSeqAln.tab.txt.gz"));
	translatedRes.writeSeqLocations(popBedLocs);

	OutputStream seqsUnableToBeMappedOut(njh::files::make_path(directory, "seqsUnableToBeMapped.txt"));
	seqsUnableToBeMappedOut << njh::conToStr(translatedRes.seqsUnableToBeMapped_, "\n") << std::endl;
	OutputStream seqsTranslationFilteredOut(njh::files::make_path(directory, "seqsTranslationFiltered.txt"));
	seqsTranslationFilteredOut << njh::conToStr(translatedRes.seqsTranslationFiltered_, "\n") << std::endl;

	if (!translatedRes.translations_.empty()) {
		SeqOutput transwriter(SeqIOOptions::genFastaOutGz(njh::files::make_path(directory, "translatedInput.fasta.gz")));
		auto seqNames = njh::getVecOfMapKeys(translatedRes.translations_);
		njh::sort(seqNames);
		for (const auto & seqName : seqNames) {
			for (const auto & transcript : translatedRes.translations_.at(seqName)) {
				transwriter.openWrite(transcript.second.translation_);
			}
		}
	}

	OutputStream variantCalls(njh::files::make_path(directory, "varCalls.vcf"));
	for (const auto & seqVar : translatedRes.seqVariants_) {
		seqVar.second.writeVCF(variantCalls);
	}
}

}  // namespace


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

  SeqInput seqReader(setUp.pars_.ioOptions_);
  const auto input = seqReader.readAllReads<seqInfo>();
  const auto counts = genHapSampleCounts(input);
	auto translatedRes = translator.run(setUp.pars_.ioOptions_, counts, variantCallerRunPars);

	TranslatedSeqsWriter writer(setUp.pars_.directoryName_, setUp.pars_.debug_, false);
	writer.write(translatedRes);
	writeTranslationRunInfo(setUp.pars_.directoryName_, translatedRes);

	return 0;
}


int genExpRunner::translateSeqsAgainstBestGenome(const njh::progutils::CmdArgs & inputCommands){

	TranslatorByAlignment::RunPars variantCallerRunPars;
	variantCallerRunPars.lowVariantCutOff = 0;
	variantCallerRunPars.occurrenceCutOff = 0;
	BioCmdsUtils::LastZPars lzPars;
	lzPars.coverage = 95;
	lzPars.identity = 90;
	uint32_t lencutOffForBowtie = 200;
	uint32_t numThreads = 1;
	std::string genomesStr;
	std::string excludeGenomesStr;
	bool keepTemporaryFiles = false;
	MultiGenomeMapper::inputParameters mapperPars;

	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.description_ = "Align seqs to every genome in a directory, then translate each seq with the gff of the genome it aligns to best (highest alignment score), seqs that fail translation filtering are tried against their next best genome";
	setUp.setOption(mapperPars.genomeDir_, "--genomeDir", "A directory of genomes, only genomes with a gff of the same name in --gffDir are used", true);
	setUp.setOption(mapperPars.gffDir_, "--gffDir", "A directory of gffs named after the genomes, e.g. Pf3D7.gff for Pf3D7.fasta", true);
	setUp.setOption(genomesStr, "--genomes", "Only use these genomes, comma separated, without extensions");
	setUp.setOption(excludeGenomesStr, "--excludeGenomes", "Do not use these genomes, comma separated, without extensions");
	setUp.setOption(mapperPars.primaryGenome_, "--primaryGenome", "Use this genome when it ties for the best score, otherwise ties go to the first genome by name");
	setUp.setOption(lzPars.coverage, "--coverage", "coverage for lastz");
	setUp.setOption(lzPars.identity, "--identity", "identity for lastz");
	setUp.setOption(lencutOffForBowtie, "--lencutOffForBowtie", "Seqs this length or shorter are aligned with bowtie2 instead of lastz");
	setUp.setOption(numThreads, "--numThreads", "number of cpus to use");
	setUp.setOption(keepTemporaryFiles, "--keepTemporaryFiles", "Keep the alignments to every genome");
	setUp.setOption(variantCallerRunPars.realnPars.extendAmount, "--extendAmount", "Amount to extend for re-alignment");
	setUp.processReadInNames();
	setUp.processDirectoryOutputName(true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);

	njh::sys::requireExternalProgramThrow("samtools");

	MultiGenomeMapper gMapper(mapperPars);
	gMapper.pars_.numThreads_ = numThreads;
	if (!genomesStr.empty()) {
		gMapper.setSelectedGenomes(genomesStr);
	}
	gMapper.loadInGenomes();
	gMapper.loadGffFnps();
	for (const auto & genome : njh::tokenizeString(excludeGenomesStr, ",")) {
		gMapper.checkForGenomeThrow(genome, __PRETTY_FUNCTION__);
		gMapper.genomes_.erase(genome);
	}
	VecStr genomesWithoutGffs;
	for (auto genome = gMapper.genomes_.begin(); genome != gMapper.genomes_.end();) {
		if (genome->second->gffFnp_.empty()) {
			genomesWithoutGffs.emplace_back(genome->first);
			genome = gMapper.genomes_.erase(genome);
		} else {
			genome->second->createTwoBit();
			++genome;
		}
	}
	if (gMapper.genomes_.empty()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "no genomes in " << mapperPars.genomeDir_ << " with a gff in " << mapperPars.gffDir_ << "\n";
		throw std::runtime_error{ss.str()};
	}
	if (!genomesWithoutGffs.empty()) {
		njh::sort(genomesWithoutGffs);
		std::cout << "Skipping genomes without a gff: " << njh::conToStr(genomesWithoutGffs, ",") << std::endl;
	}
	if (!mapperPars.primaryGenome_.empty()) {
		gMapper.checkForGenomeThrow(mapperPars.primaryGenome_, __PRETTY_FUNCTION__);
	}

	SeqInput seqReader(setUp.pars_.ioOptions_);
	const auto input = seqReader.readAllReads<seqInfo>();
	const auto counts = genHapSampleCounts(input);

	//write the seqs for aligning with their position as the name so the names are safe for the aligners
	const auto alignmentsDir = njh::files::make_path(setUp.pars_.directoryName_, "alignments");
	njh::files::makeDir(njh::files::MkdirPar{alignmentsDir});
	const auto forLastzFnp = njh::files::make_path(alignmentsDir, "forLastz.fasta");
	const auto forBowtie2Fnp = njh::files::make_path(alignmentsDir, "forBowtie2.fasta");
	std::unordered_map<std::string, uint32_t> inputPositions;
	bool hasLastzSeqs = false;
	bool hasBowtie2Seqs = false;
	{
		SeqOutput lastzWriter(SeqIOOptions::genFastaOut(forLastzFnp));
		SeqOutput bowtie2Writer(SeqIOOptions::genFastaOut(forBowtie2Fnp));
		for (const auto & seq : iter::enumerate(input)) {
			if (!inputPositions.emplace(seq.element.name_, static_cast<uint32_t>(seq.index)).second) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "input names must be unique, found " << seq.element.name_ << " more than once" << "\n";
				throw std::runtime_error{ss.str()};
			}
			seqInfo renamed(estd::to_string(seq.index), seq.element.seq_);
			if (len(seq.element) <= lencutOffForBowtie) {
				bowtie2Writer.openWrite(renamed);
				hasBowtie2Seqs = true;
			} else {
				lastzWriter.openWrite(renamed);
				hasLastzSeqs = true;
			}
		}
	}

	//the bam for each genome per aligner
	struct GenomeBam {
		std::string genome_;
		bfs::path bamFnp_;
	};
	std::vector<GenomeBam> genomeBams;
	if (hasLastzSeqs) {
		njh::sys::requireExternalProgramThrow("lastz");
		const auto lastzOutputs = gMapper.alignToGenomesLastz(SeqIOOptions::genFastaIn(forLastzFnp),
				njh::files::make_path(alignmentsDir, "lastz_"), lzPars);
		for (const auto & output : lastzOutputs) {
			genomeBams.emplace_back(GenomeBam{output.first, output.second.alignedFnp_});
		}
	}
	if (hasBowtie2Seqs) {
		njh::sys::requireExternalProgramThrow("bowtie2");
		for (const auto & genome : gMapper.genomes_) {
			genome.second->buildBowtie2Index();
		}
		const auto bowtie2Outputs = gMapper.alignToGenomes(SeqIOOptions::genFastaIn(forBowtie2Fnp),
				njh::files::make_path(alignmentsDir, "bowtie2_"));
		for (const auto & output : bowtie2Outputs) {
			genomeBams.emplace_back(GenomeBam{output.first, output.second.alignedFnp_});
		}
	}

	//score the alignments, keeping the best alignment of each seq in each genome
	struct GenomeHit {
		double score_{std::numeric_limits<double>::lowest()};
		double identity_{0};
		std::string region_;
	};
	//key1 = seq position, key2 = genome
	std::unordered_map<uint32_t, std::map<std::string, GenomeHit>> hits;
	{
		std::mutex hitsMut;
		njh::concurrent::LockableQueue<GenomeBam> bamQueue(genomeBams);
		std::function<void()> scoreBams = [&bamQueue, &hits, &hitsMut, &gMapper]() {
			GenomeBam genomeBam;
			while (bamQueue.getVal(genomeBam)) {
				TwoBit::TwoBitFile twobitReader(gMapper.genomes_.at(genomeBam.genome_)->fnpTwoBit_);
				BamTools::BamReader bReader;
				bReader.Open(genomeBam.bamFnp_.string());
				checkBamOpenThrow(bReader, genomeBam.bamFnp_.string());
				const auto refData = bReader.GetReferenceData();
				std::unordered_map<uint32_t, GenomeHit> genomeHits;
				BamTools::BamAlignment bAln;
				while (bReader.GetNextAlignment(bAln)) {
					if (bAln.IsMapped()) {
						AlignmentResults results(bAln, refData);
						results.setRefSeq(twobitReader);
						results.setComparison(false);
						auto & hit = genomeHits[njh::StrToNumConverter::stoToNum<uint32_t>(bAln.Name)];
						if (results.comp_.alnScore_ > hit.score_) {
							hit.score_ = results.comp_.alnScore_;
							hit.identity_ = results.comp_.distances_.eventBasedIdentity_;
							hit.region_ = results.gRegion_.createUidFromCoordsStrand();
						}
					}
				}
				std::lock_guard<std::mutex> lock(hitsMut);
				for (const auto & hit : genomeHits) {
					hits[hit.first][genomeBam.genome_] = hit.second;
				}
			}
		};
		njh::concurrent::runVoidFunctionThreaded(scoreBams, numThreads);
	}

	//rank the genomes for each seq by score, ties go to the primary genome then by name
	struct SeqPlacement {
		VecStr rankedGenomes_;
		VecStr tiedGenomes_;
		std::string nextBestGenome_{"NA"};
		double nextBestScore_{std::numeric_limits<double>::lowest()};
		VecStr genomesFiltered_;
		std::string translatedGenome_;
	};
	std::vector<SeqPlacement> placements(input.size());
	std::vector<uint32_t> seqsToTranslate;
	VecStr unmappedToAnyGenome;
	{
		OutputStream allHitsOut(njh::files::make_path(setUp.pars_.directoryName_, "alignmentsPerGenome.tab.txt"));
		allHitsOut << "seqName\tgenome\tregion\talnScore\tidentity" << std::endl;
		for (const auto & seq : iter::enumerate(input)) {
			const auto seqPosition = static_cast<uint32_t>(seq.index);
			const auto seqHits = hits.find(seqPosition);
			if (hits.end() == seqHits) {
				unmappedToAnyGenome.emplace_back(seq.element.name_);
				continue;
			}
			auto & placement = placements[seqPosition];
			double bestScore = std::numeric_limits<double>::lowest();
			for (const auto & hit : seqHits->second) {
				allHitsOut << seq.element.name_ << "\t" << hit.first << "\t" << hit.second.region_ << "\t" << hit.second.score_ << "\t" << hit.second.identity_ << std::endl;
				bestScore = std::max(bestScore, hit.second.score_);
				placement.rankedGenomes_.emplace_back(hit.first);
			}
			for (const auto & hit : seqHits->second) {
				if (hit.second.score_ == bestScore) {
					placement.tiedGenomes_.emplace_back(hit.first);
				} else if (hit.second.score_ > placement.nextBestScore_) {
					placement.nextBestScore_ = hit.second.score_;
					placement.nextBestGenome_ = hit.first;
				}
			}
			std::sort(placement.rankedGenomes_.begin(), placement.rankedGenomes_.end(),
					[&seqHits, &mapperPars](const std::string & genome1, const std::string & genome2) {
						const auto score1 = seqHits->second.at(genome1).score_;
						const auto score2 = seqHits->second.at(genome2).score_;
						if (score1 != score2) {
							return score1 > score2;
						}
						const bool genome1Primary = genome1 == mapperPars.primaryGenome_;
						const bool genome2Primary = genome2 == mapperPars.primaryGenome_;
						if (genome1Primary != genome2Primary) {
							return genome1Primary;
						}
						return genome1 < genome2;
					});
			seqsToTranslate.emplace_back(seqPosition);
		}
	}

	//translate each seq against its best genome, seqs filtered off go on to their next best genome
	struct TranslationGroup {
		std::string genome_;
		bool useLastz_{false};
		std::vector<seqInfo> seqs_;
		TranslatorByAlignment::TranslatorByAlignmentResult result_;
	};
	std::map<std::string, TranslationGroup> groups;
	const auto byGenomeDir = njh::files::make_path(setUp.pars_.directoryName_, "byGenome");
	njh::files::makeDir(njh::files::MkdirPar{byGenomeDir});
	uint32_t retryRound = 0;
	while (!seqsToTranslate.empty()) {
		std::set<std::string> roundGroups;
		for (const auto seqPosition : seqsToTranslate) {
			const auto & placement = placements[seqPosition];
			const auto & genome = placement.rankedGenomes_[placement.genomesFiltered_.size()];
			const bool useLastz = len(input[seqPosition]) > lencutOffForBowtie;
			const auto groupName = njh::pasteAsStr(genome, useLastz ? "-lastz" : "-bowtie2", retryRound > 0 ? njh::pasteAsStr("-retry", retryRound) : std::string(""));
			auto & group = groups[groupName];
			group.genome_ = genome;
			group.useLastz_ = useLastz;
			group.seqs_.emplace_back(input[seqPosition]);
			roundGroups.emplace(groupName);
		}
		njh::concurrent::LockableQueue<std::string> groupQueue(VecStr(roundGroups.begin(), roundGroups.end()));
		std::function<void()> translateGroups = [&groupQueue, &groups, &gMapper, &byGenomeDir, &lzPars, &counts, &variantCallerRunPars, &setUp, &keepTemporaryFiles]() {
			std::string groupName;
			while (groupQueue.getVal(groupName)) {
				auto & group = groups.at(groupName);
				const auto groupDir = njh::files::make_path(byGenomeDir, groupName);
				njh::files::makeDir(njh::files::MkdirPar{groupDir});
				const auto groupInputFnp = njh::files::make_path(groupDir, "input.fasta");
				std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>> groupCounts;
				{
					SeqOutput groupWriter(SeqIOOptions::genFastaOut(groupInputFnp));
					for (const auto & seq : group.seqs_) {
						groupWriter.openWrite(seq);
						groupCounts[seq.name_] = counts.at(seq.name_);
					}
				}
				TranslatorByAlignment::TranslatorByAlignmentPars transPars;
				transPars.keepTemporaryFiles_ = keepTemporaryFiles;
				transPars.useLastz_ = group.useLastz_;
				transPars.gffFnp_ = gMapper.genomes_.at(group.genome_)->gffFnp_;
				transPars.lzPars_.genomeFnp = gMapper.genomes_.at(group.genome_)->fnp_;
				transPars.lzPars_.coverage = lzPars.coverage;
				transPars.lzPars_.identity = lzPars.identity;
				transPars.workingDirtory_ = groupDir;
				TranslatorByAlignment translator(transPars);
				group.result_ = translator.run(SeqIOOptions::genFastaIn(groupInputFnp), groupCounts, variantCallerRunPars);

				TranslatedSeqsWriter groupSeqsWriter(groupDir, setUp.pars_.debug_, false);
				groupSeqsWriter.write(group.result_);
				writeTranslationRunInfo(groupDir, group.result_);
			}
		};
		njh::concurrent::runVoidFunctionThreaded(translateGroups, numThreads);

		seqsToTranslate.clear();
		for (const auto & groupName : roundGroups) {
			const auto & group = groups.at(groupName);
			for (const auto & translation : group.result_.translations_) {
				placements[inputPositions.at(translation.first)].translatedGenome_ = group.genome_;
			}
			for (const auto & seqName : group.result_.seqsTranslationFiltered_) {
				const auto seqPosition = inputPositions.at(seqName);
				auto & placement = placements[seqPosition];
				placement.genomesFiltered_.emplace_back(group.genome_);
				if (placement.genomesFiltered_.size() < placement.rankedGenomes_.size()) {
					seqsToTranslate.emplace_back(seqPosition);
				}
			}
		}
		njh::sort(seqsToTranslate);
		++retryRound;
	}

	//the combined output
	{
		OutputStream bestGenomeOut(njh::files::make_path(setUp.pars_.directoryName_, "bestGenomes.tab.txt"));
		bestGenomeOut << "seqName\tbestGenome\tregion\talnScore\tidentity\ttiedGenomes\tnextBestGenome\tnextBestScore\tgenomesAlignedTo\ttranslatedGenome\ttranslatedAlnScore\tgenomesFiltered" << std::endl;
		for (const auto & seq : iter::enumerate(input)) {
			const auto & placement = placements[seq.index];
			if (placement.rankedGenomes_.empty()) {
				continue;
			}
			const auto & seqHits = hits.at(static_cast<uint32_t>(seq.index));
			const auto & bestHit = seqHits.at(placement.rankedGenomes_.front());
			bestGenomeOut << seq.element.name_
					<< "\t" << placement.rankedGenomes_.front()
					<< "\t" << bestHit.region_
					<< "\t" << bestHit.score_
					<< "\t" << bestHit.identity_
					<< "\t" << njh::conToStr(placement.tiedGenomes_, ",")
					<< "\t" << placement.nextBestGenome_
					<< "\t" << ("NA" == placement.nextBestGenome_ ? std::string("NA") : estd::to_string(placement.nextBestScore_))
					<< "\t" << seqHits.size()
					<< "\t" << (placement.translatedGenome_.empty() ? std::string("NA") : placement.translatedGenome_)
					<< "\t" << (placement.translatedGenome_.empty() ? std::string("NA") : estd::to_string(seqHits.at(placement.translatedGenome_).score_))
					<< "\t" << (placement.genomesFiltered_.empty() ? std::string("NA") : njh::conToStr(placement.genomesFiltered_, ",")) << std::endl;
		}
	}
	TranslatedSeqsWriter writer(setUp.pars_.directoryName_, setUp.pars_.debug_, true);
	VecStr seqsUnableToBeMapped = unmappedToAnyGenome;
	VecStr seqsTranslationFiltered;
	for (const auto & group : groups) {
		writer.write(group.second.result_, group.second.genome_);
		for (const auto & seqName : group.second.result_.seqsUnableToBeMapped_) {
			if (placements[inputPositions.at(seqName)].genomesFiltered_.empty()) {
				seqsUnableToBeMapped.emplace_back(seqName);
			}
		}
	}
	for (const auto & seq : iter::enumerate(input)) {
		const auto & placement = placements[seq.index];
		if (placement.translatedGenome_.empty() && !placement.genomesFiltered_.empty()) {
			seqsTranslationFiltered.emplace_back(seq.element.name_);
		}
	}
	njh::sort(seqsUnableToBeMapped);
	njh::sort(seqsTranslationFiltered);
	OutputStream seqsUnableToBeMappedOut(njh::files::make_path(setUp.pars_.directoryName_, "seqsUnableToBeMapped.txt"));
	seqsUnableToBeMappedOut << njh::conToStr(seqsUnableToBeMapped, "\n") << std::endl;
	OutputStream seqsTranslationFilteredOut(njh::files::make_path(setUp.pars_.directoryName_, "seqsTranslationFiltered.txt"));
	seqsTranslationFilteredOut << njh::conToStr(seqsTranslationFiltered, "\n") << std::endl;

	if (!keepTemporaryFiles) {
		njh::files::rmDirForce(alignmentsDir);
	}
	return 0;
}

} /* namespace njhseq */


