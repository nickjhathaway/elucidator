//
// build_pmo_from_allele_table, a PMO from an allele table and the panel's primers,
// optionally with the targets' genome locations, sample meta, and variant calls already made per target
//

#include "pmo_utils_runner.hpp"

#include <njhseq/objects/BioDataObject/pmo/pmoUtils.hpp>
#include <njhseq/objects/dataContainers/tables/TableReader.hpp>
#include <njhseq/objects/Meta/MetaDataInName.hpp>
#include <njhseq/IO/SeqIO/SeqInput.hpp>
#include <njhseq/helpers/seqUtil.hpp>

namespace njhseq {

namespace {

/**@brief The columns of the allele table, one row per haplotype per sample per target
 *
 */
struct AlleleTableCols {
	std::string sample_{"SampleID"};
	std::string target_{"Locus"};
	std::string seq_{"ASV"};
	std::string reads_{"Reads"};
	std::string mhapName_; /**< optional, a name for each sequence of a target*/
	std::string specimen_; /**< optional, the specimen of each library sample, defaults to the library sample*/
};

/**@brief Primers per target from a SeekDeep id table, the first three columns are target, forward and reverse primer
 *
 * @param fnp the table, with a header
 * @return the targets in the order of the table with their forward and reverse primers
 */
std::vector<std::tuple<std::string, std::string, std::string>> readPrimers(const bfs::path & fnp) {
	std::vector<std::tuple<std::string, std::string, std::string>> ret;
	std::set<std::string> targets;
	TableReader reader(TableIOOpts::genTabFileIn(fnp, true));
	if (reader.header_.nCol() < 3) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error, " << fnp << " should have at least 3 columns, target, forward primer, reverse primer\n";
		throw std::runtime_error{ss.str()};
	}
	VecStr row;
	while (reader.getNextRow(row)) {
		if (njh::in(row[0], targets)) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error, target " << row[0] << " is in " << fnp << " more than once\n";
			throw std::runtime_error{ss.str()};
		}
		targets.emplace(row[0]);
		ret.emplace_back(row[0], njh::strToUpperRet(row[1]), njh::strToUpperRet(row[2]));
	}
	return ret;
}

/**@brief Generate a PMO with the required fields from an allele table, every target in the primers is a target of the one panel
 *
 * Every sample in the table is a library sample, and a specimen of the same name unless cols.specimen_ is set.
 * Rows for the same sample, target and sequence are summed and reads are rounded.
 *
 * @param alleleTableFnp the allele table
 * @param cols the allele table's columns
 * @param primers the targets and their primers, from readPrimers
 * @param panelName the panel's name
 * @param generationMethod the program generating the PMO
 * @return the validated PMO
 */
pmo::PortableMicrohaplotypeObject genPmoFromAlleleTable(const bfs::path & alleleTableFnp,
		const AlleleTableCols & cols,
		const std::vector<std::tuple<std::string, std::string, std::string>> & primers,
		const std::string & panelName,
		const pmo::PmoGenerationMethod & generationMethod) {
	pmo::PortableMicrohaplotypeObject ret;
	auto today = njh::getCurrentDate();
	ret.pmo_header_.pmo_version_ = pmo::PMO_SCHEMA_VERSION;
	ret.pmo_header_.creation_date_ = today.substr(0, today.find('_'));
	ret.pmo_header_.generation_method_ = generationMethod;

	std::unordered_map<std::string, uint32_t> targetIds;
	pmo::ReactionInfo reaction;
	reaction.reaction_name_ = "full";
	for (const auto & primer : primers) {
		pmo::TargetInfo target;
		target.target_name_ = std::get<0>(primer);
		target.forward_primer_.seq_ = std::get<1>(primer);
		target.reverse_primer_.seq_ = std::get<2>(primer);
		targetIds[target.target_name_] = ret.target_info_.size();
		reaction.panel_targets_.emplace_back(ret.target_info_.size());
		ret.target_info_.emplace_back(target);
	}
	pmo::PanelInfo panel;
	panel.panel_name_ = panelName;
	panel.reactions_.emplace_back(reaction);
	ret.panel_info_.emplace_back(panel);

	TableReader reader(TableIOOpts::genTabFileIn(alleleTableFnp, true));
	VecStr requiredCols{cols.sample_, cols.target_, cols.seq_, cols.reads_};
	if (!cols.mhapName_.empty()) {
		requiredCols.emplace_back(cols.mhapName_);
	}
	if (!cols.specimen_.empty()) {
		requiredCols.emplace_back(cols.specimen_);
	}
	reader.header_.checkForColumnsThrow(requiredCols, __PRETTY_FUNCTION__);
	const auto samplePos = reader.header_.getColPos(cols.sample_);
	const auto targetPos = reader.header_.getColPos(cols.target_);
	const auto seqPos = reader.header_.getColPos(cols.seq_);
	const auto readsPos = reader.header_.getColPos(cols.reads_);

	//sequences per target in the order they're first seen
	std::map<uint32_t, VecStr> seqsByTargetId;
	std::unordered_map<uint32_t, std::unordered_map<std::string, uint32_t>> mhapIdBySeqByTargetId;
	std::unordered_map<uint32_t, std::unordered_map<std::string, std::string>> mhapNameBySeqByTargetId;
	//library samples in the order they're first seen
	VecStr samples;
	std::unordered_map<std::string, std::string> specimenBySample;
	//key1 = sample, key2 = target id, key3 = mhap id, value = reads
	std::unordered_map<std::string, std::map<uint32_t, std::map<uint32_t, double>>> readsBySample;
	std::set<std::string> targetsNotInPrimers;
	VecStr row;
	while (reader.getNextRow(row)) {
		const auto & sample = row[samplePos];
		const auto & targetName = row[targetPos];
		const auto seq = njh::strToUpperRet(row[seqPos]);
		auto targetSearch = targetIds.find(targetName);
		if (targetIds.end() == targetSearch) {
			targetsNotInPrimers.emplace(targetName);
			continue;
		}
		const uint32_t targetId = targetSearch->second;
		auto & mhapIdBySeq = mhapIdBySeqByTargetId[targetId];
		auto seqSearch = mhapIdBySeq.find(seq);
		if (mhapIdBySeq.end() == seqSearch) {
			seqSearch = mhapIdBySeq.emplace(seq, seqsByTargetId[targetId].size()).first;
			seqsByTargetId[targetId].emplace_back(seq);
		}
		if (!cols.mhapName_.empty()) {
			const auto & name = row[reader.header_.getColPos(cols.mhapName_)];
			auto & nameBySeq = mhapNameBySeqByTargetId[targetId];
			auto nameSearch = nameBySeq.find(seq);
			if (nameBySeq.end() == nameSearch) {
				nameBySeq[seq] = name;
			} else if (nameSearch->second != name) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error, the same sequence of target " << targetName << " has more than one name in " << cols.mhapName_ << ": " << nameSearch->second << " and " << name << "\n";
				throw std::runtime_error{ss.str()};
			}
		}
		const std::string specimen = cols.specimen_.empty() ? sample : row[reader.header_.getColPos(cols.specimen_)];
		auto specimenSearch = specimenBySample.find(sample);
		if (specimenBySample.end() == specimenSearch) {
			specimenBySample[sample] = specimen;
			samples.emplace_back(sample);
		} else if (specimenSearch->second != specimen) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error, library sample " << sample << " has more than one specimen in " << cols.specimen_ << ": " << specimenSearch->second << " and " << specimen << "\n";
			throw std::runtime_error{ss.str()};
		}
		readsBySample[sample][targetId][seqSearch->second] += njh::StrToNumConverter::stoToNum<double>(row[readsPos]);
	}
	if (!targetsNotInPrimers.empty()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error, targets in " << alleleTableFnp << " don't have primers: " << njh::conToStr(targetsNotInPrimers, ",") << "\n";
		throw std::runtime_error{ss.str()};
	}

	for (const auto & targetSeqs : seqsByTargetId) {
		pmo::RepresentativeMicrohaplotypesForTarget representatives;
		representatives.target_id_ = targetSeqs.first;
		for (const auto & seq : targetSeqs.second) {
			pmo::RepresentativeMicrohaplotype mhap;
			mhap.seq_ = seq;
			if (!cols.mhapName_.empty()) {
				mhap.microhaplotype_name_ = mhapNameBySeqByTargetId.at(targetSeqs.first).at(seq);
			}
			representatives.microhaplotypes_.emplace_back(mhap);
		}
		ret.representative_microhaplotypes_.targets_.emplace_back(representatives);
	}
	//detected haplotypes point to the position of the target in representative_microhaplotypes
	std::unordered_map<uint32_t, uint32_t> repPosByTargetId;
	for (const auto & rep : iter::enumerate(ret.representative_microhaplotypes_.targets_)) {
		repPosByTargetId[rep.element.target_id_] = rep.index;
	}

	std::unordered_map<std::string, uint32_t> specimenIds;
	pmo::DetectedMicrohaplotypes detected;
	for (const auto & sample : samples) {
		const auto & specimenName = specimenBySample.at(sample);
		auto specimenSearch = specimenIds.find(specimenName);
		if (specimenIds.end() == specimenSearch) {
			specimenSearch = specimenIds.emplace(specimenName, ret.specimen_info_.size()).first;
			pmo::SpecimenInfo specimen;
			specimen.specimen_name_ = specimenName;
			ret.specimen_info_.emplace_back(specimen);
		}
		const uint32_t libraryId = ret.library_sample_info_.size();
		pmo::LibrarySampleInfo library;
		library.library_sample_name_ = sample;
		library.specimen_id_ = specimenSearch->second;
		library.panel_id_ = 0;
		ret.library_sample_info_.emplace_back(library);

		pmo::DetectedMicrohaplotypesForSample sampleResult;
		sampleResult.library_sample_id_ = libraryId;
		for (const auto & targetReads : readsBySample.at(sample)) {
			pmo::DetectedMicrohaplotypesForTarget targetResult;
			targetResult.mhaps_target_id_ = repPosByTargetId.at(targetReads.first);
			for (const auto & mhapReads : targetReads.second) {
				pmo::MicrohaplotypeForTarget mhap;
				mhap.mhap_id_ = mhapReads.first;
				mhap.reads_ = static_cast<uint32_t>(std::round(mhapReads.second));
				targetResult.mhaps_.emplace_back(mhap);
			}
			sampleResult.target_results_.emplace_back(targetResult);
		}
		detected.library_samples_.emplace_back(sampleResult);
	}
	ret.detected_microhaplotypes_.emplace_back(detected);
	ret.validate();
	return ret;
}

/**@brief What addVariantCallsFromDir added
 *
 */
struct AddVariantCallsFromDirResults {
	uint32_t targetsWithCalls_{0};
	VecStr targetsWithoutCalls_; /**< targets with haplotypes but no variant call directory*/
	uint32_t mhapsWithCalls_{0};
	uint32_t mhapsWithoutCalls_{0};
	uint32_t proteinVariants_{0};
	uint32_t seqVariants_{0};
	uint32_t codonsWithoutAlt_{0}; /**< protein variants whose haplotype codon couldn't be worked out from the tables*/
};

/**@brief One base of a codon from geneInfos/gene_<transcript>_basePositions.tab.txt
 *
 */
struct CodonBase {
	uint32_t codonPos_;
	uint32_t gPos_;
	char base_; /**< on the transcript's strand*/
};

/**@brief The positions of a transcript's codons, its gene and the gene's Name from the gff
 *
 */
struct TranscriptCodons {
	std::string chrom_;
	char strand_{'+'};
	std::string geneId_;
	std::string altGeneName_; /**< empty when the gene has no Name in the gff or it's the gene id*/
	std::unordered_map<uint32_t, std::vector<CodonBase>> basesByZeroBasedAAPos_;
};

/**@brief Read the codon positions of a transcript from the geneInfos directory collapseAndCallVariants writes
 *
 * @param geneInfosDir the directory
 * @param transcript the transcript
 * @return the codons
 */
TranscriptCodons readTranscriptCodons(const bfs::path & geneInfosDir, const std::string & transcript) {
	TranscriptCodons ret;
	const auto basePositionsFnp = njh::files::make_path(geneInfosDir, "gene_" + transcript + "_basePositions.tab.txt");
	TableReader reader(TableIOOpts::genTabFileIn(basePositionsFnp, true));
	reader.header_.checkForColumnsThrow(VecStr{"chrom", "gDnaPos", "aaPos", "codonPos", "base", "GeneID", "strand"}, __PRETTY_FUNCTION__);
	VecStr row;
	while (reader.getNextRow(row)) {
		ret.chrom_ = row[reader.header_.getColPos("chrom")];
		ret.geneId_ = row[reader.header_.getColPos("GeneID")];
		ret.strand_ = row[reader.header_.getColPos("strand")].front();
		ret.basesByZeroBasedAAPos_[njh::StrToNumConverter::stoToNum<uint32_t>(row[reader.header_.getColPos("aaPos")])].emplace_back(CodonBase{
			njh::StrToNumConverter::stoToNum<uint32_t>(row[reader.header_.getColPos("codonPos")]),
			njh::StrToNumConverter::stoToNum<uint32_t>(row[reader.header_.getColPos("gDnaPos")]),
			row[reader.header_.getColPos("base")].front()});
	}
	for (auto & codon : ret.basesByZeroBasedAAPos_) {
		std::sort(codon.second.begin(), codon.second.end(), [](const CodonBase & a, const CodonBase & b) {
			return a.codonPos_ < b.codonPos_;
		});
	}
	//the gene's Name from the gff is in the extra field of gene_<transcript>.bed, characters the schema doesn't allow become _
	const auto bedFnp = njh::files::make_path(geneInfosDir, "gene_" + transcript + ".bed");
	if (bfs::exists(bedFnp)) {
		InputStream bedIn(bedFnp);
		std::string line;
		while (njh::files::crossPlatGetline(bedIn, line)) {
			auto toks = tokenizeString(line, "\t");
			if (toks.size() < 7 || !MetaDataInName::nameHasMetaData(toks[6])) {
				continue;
			}
			MetaDataInName meta(toks[6]);
			if (!meta.containsMeta("Name")) {
				continue;
			}
			auto name = meta.getMeta("Name");
			if (name.empty() || njh::in(name, pmo::PMO_NA_STRINGS()) || name == ret.geneId_) {
				continue;
			}
			for (auto & letter : name) {
				if (!std::isalnum(static_cast<unsigned char>(letter)) && '-' != letter && '.' != letter && '_' != letter) {
					letter = '_';
				}
			}
			ret.altGeneName_ = name;
			break;
		}
	}
	return ret;
}

/**@brief Add variant calls already made per target, the output of collapseAndCallVariants, to the representative microhaplotypes
 *
 * This gives the same variants as SeekDeep variantCallOnSeqAndProtein --pmo, without calling them again.
 * A target's calls are in <variantCallsDir>/<target>/<subDir>, haplotypes are matched to its uniqueSeqs.fasta.gz by sequence,
 * protein variants come from variantCalls/seqsAATyped.tab.txt.gz (every typed amino acid) with the codons from variantCalls/geneInfos,
 * and sequence variants from variantCalls/variantsPerSeqAln.tab.txt.gz. Any variants the haplotypes already have are replaced.
 * The haplotype's codon is the reference codon with the haplotype's SNPs, or the codon in variantCalls/variantsPerTranslatedSeq.tab.txt.gz
 * when an indel is in the codon or the SNP codon doesn't translate to the haplotype's amino acid, and the codon's alt_seq is left out when neither works.
 *
 * @param pmoObj the PMO to update
 * @param variantCallsDir the directory with a directory per target
 * @param subDir the directory with the calls within each target's directory
 * @param genomeName the genome the variants were called against, added to targeted_genomes if it isn't there
 * @return what was added
 */
AddVariantCallsFromDirResults addVariantCallsFromDir(pmo::PortableMicrohaplotypeObject & pmoObj,
		const bfs::path & variantCallsDir,
		const std::string & subDir,
		const std::string & genomeName) {
	AddVariantCallsFromDirResults ret;
	if (!pmoObj.targeted_genomes_.has_value()) {
		pmoObj.targeted_genomes_ = std::vector<pmo::GenomeInfo>{};
	}
	auto & genomes = pmoObj.targeted_genomes_.value();
	auto genomeSearch = std::find_if(genomes.begin(), genomes.end(), [&genomeName](const pmo::GenomeInfo & genome) {
		return genome.name_ == genomeName;
	});
	uint32_t genomeId = std::distance(genomes.begin(), genomeSearch);
	bool genomeNeeded = genomes.end() == genomeSearch;

	auto genLocation = [&genomeId](const std::string & chrom, const uint32_t start, const uint32_t end,
			const std::string & refSeq, const std::optional<std::string> & altSeq, const std::optional<char> strand) {
		pmo::GenomicLocation loc;
		loc.genome_id_ = genomeId;
		loc.chrom_ = chrom;
		loc.start_ = start;
		loc.end_ = end;
		loc.ref_seq_ = refSeq;
		loc.alt_seq_ = altSeq;
		if (strand.has_value()) {
			loc.strand_ = std::string(1, strand.value());
		}
		return loc;
	};
	auto translate = [](const std::string & codon) {
		if (3 != codon.size() || std::string::npos != codon.find_first_not_of("ACGT")) {
			return '?';
		}
		return seqInfo("codon", codon).translateRet(false, false).seq_.front();
	};

	for (auto & reps : pmoObj.representative_microhaplotypes_.targets_) {
		const auto & targetName = pmoObj.target_info_.at(reps.target_id_).target_name_;
		const auto targetDir = njh::files::make_path(variantCallsDir, targetName, subDir);
		const auto uniqueSeqsFnp = njh::files::make_path(targetDir, "uniqueSeqs.fasta.gz");
		if (!bfs::exists(uniqueSeqsFnp)) {
			ret.targetsWithoutCalls_.emplace_back(targetName);
			continue;
		}
		++ret.targetsWithCalls_;
		if (genomeNeeded) {
			pmo::GenomeInfo genome;
			genome.name_ = genomeName;
			genome.genome_version_ = genomeName;
			genome.url_ = "NA";
			genomes.emplace_back(genome);
			genomeNeeded = false;
		}
		std::unordered_map<std::string, std::string> callNameBySeq;
		{
			SeqInput reader(SeqIOOptions::genFastaInGz(uniqueSeqsFnp));
			reader.openIn();
			seqInfo seq;
			while (reader.readNextRead(seq)) {
				callNameBySeq[seq.seq_] = seq.name_;
			}
		}
		const auto callsDir = njh::files::make_path(targetDir, "variantCalls");
		//the typed amino acids per haplotype, transcripts separated by ; and amino acids by :, e.g. PF3D7_0709000.1--K76T:A220
		std::unordered_map<std::string, std::string> typedByName;
		const auto typedFnp = njh::files::make_path(callsDir, "seqsAATyped.tab.txt.gz");
		if (bfs::exists(typedFnp)) {
			TableReader reader(TableIOOpts::genTabFileIn(typedFnp, true));
			reader.header_.checkForColumnsThrow(VecStr{"name", "fullTyped"}, __PRETTY_FUNCTION__);
			VecStr row;
			while (reader.getNextRow(row)) {
				const auto & typed = row[reader.header_.getColPos("fullTyped")];
				if (!typed.empty() && "Unmappable" != typed && "Untranslatable" != typed) {
					typedByName[row[reader.header_.getColPos("name")]] = typed;
				}
			}
		}
		//the sequence variants per haplotype, in the order they're written so the sort below keeps ties like SeekDeep does
		std::unordered_map<std::string, std::vector<pmo::GenomicLocation>> seqVariantsByName;
		std::unordered_set<std::string> mappedNames;
		const auto seqVarsFnp = njh::files::make_path(callsDir, "variantsPerSeqAln.tab.txt.gz");
		if (bfs::exists(seqVarsFnp)) {
			TableReader reader(TableIOOpts::genTabFileIn(seqVarsFnp, true));
			reader.header_.checkForColumnsThrow(VecStr{"chrom", "start", "end", "queryName", "type", "refSeq", "querySeq"}, __PRETTY_FUNCTION__);
			VecStr row;
			while (reader.getNextRow(row)) {
				seqVariantsByName[row[reader.header_.getColPos("queryName")]].emplace_back(genLocation(
						row[reader.header_.getColPos("chrom")],
						njh::StrToNumConverter::stoToNum<uint32_t>(row[reader.header_.getColPos("start")]),
						njh::StrToNumConverter::stoToNum<uint32_t>(row[reader.header_.getColPos("end")]),
						row[reader.header_.getColPos("refSeq")],
						row[reader.header_.getColPos("querySeq")],
						'+'));
			}
		}
		//haplotypes that mapped to the genome, ones matching the reference have no rows in variantsPerSeqAln
		const auto inputBedFnp = njh::files::make_path(callsDir, "inputSeqs.bed");
		if (bfs::exists(inputBedFnp)) {
			InputStream bedIn(inputBedFnp);
			std::string line;
			while (njh::files::crossPlatGetline(bedIn, line)) {
				auto toks = tokenizeString(line, "\t");
				if (toks.size() >= 4 && !njh::beginsWith(line, "#")) {
					mappedNames.emplace(toks[3]);
				}
			}
		}
		//key1 = name, key2 = transcript, key3 = 1-based amino acid position, value = the haplotype's codon
		std::unordered_map<std::string, std::unordered_map<std::string, std::unordered_map<uint32_t, std::string>>> translatedCodons;
		const auto translatedFnp = njh::files::make_path(callsDir, "variantsPerTranslatedSeq.tab.txt.gz");
		if (bfs::exists(translatedFnp)) {
			TableReader reader(TableIOOpts::genTabFileIn(translatedFnp, true));
			reader.header_.checkForColumnsThrow(VecStr{"chrom", "start", "queryName", "queryCodonSeq"}, __PRETTY_FUNCTION__);
			VecStr row;
			while (reader.getNextRow(row)) {
				translatedCodons[row[reader.header_.getColPos("queryName")]][row[reader.header_.getColPos("chrom")]]
						[njh::StrToNumConverter::stoToNum<uint32_t>(row[reader.header_.getColPos("start")])] = row[reader.header_.getColPos("queryCodonSeq")];
			}
		}
		std::map<std::string, TranscriptCodons> codonsByTranscript;
		auto getCodons = [&codonsByTranscript, &callsDir](const std::string & transcript) -> const TranscriptCodons & {
			auto search = codonsByTranscript.find(transcript);
			if (codonsByTranscript.end() == search) {
				search = codonsByTranscript.emplace(transcript, readTranscriptCodons(njh::files::make_path(callsDir, "geneInfos"), transcript)).first;
			}
			return search->second;
		};

		for (auto & mhap : reps.microhaplotypes_) {
			mhap.associated_protein_variants_.reset();
			mhap.associated_seq_variants_.reset();
			auto nameSearch = callNameBySeq.find(mhap.seq_);
			if (callNameBySeq.end() == nameSearch) {
				++ret.mhapsWithoutCalls_;
				continue;
			}
			const auto & callName = nameSearch->second;
			auto typedSearch = typedByName.find(callName);
			const bool hasProtein = typedByName.end() != typedSearch;
			const bool hasSeq = njh::in(callName, mappedNames) || njh::in(callName, seqVariantsByName);
			if (!hasProtein && !hasSeq) {
				++ret.mhapsWithoutCalls_;
				continue;
			}
			++ret.mhapsWithCalls_;
			//SNPs on the + strand by position and the spans of indels, to work out the haplotype's codons
			std::unordered_map<uint32_t, char> snpsByPos;
			std::vector<std::pair<uint32_t, uint32_t>> indelSpans;
			auto seqVarsSearch = seqVariantsByName.find(callName);
			if (seqVariantsByName.end() != seqVarsSearch) {
				for (const auto & var : seqVarsSearch->second) {
					if (1 == var.ref_seq_->size() && 1 == var.alt_seq_->size() && '-' != var.ref_seq_->front() && '-' != var.alt_seq_->front()) {
						snpsByPos[var.start_] = var.alt_seq_->front();
					} else {
						indelSpans.emplace_back(var.start_, var.end_);
					}
				}
			}
			if (hasProtein) {
				std::vector<pmo::ProteinVariant> proteinVariants;
				auto transcriptsTyped = tokenizeString(typedSearch->second, ";");
				njh::sort(transcriptsTyped);
				for (const auto & transcriptTyped : transcriptsTyped) {
					const auto sepPos = transcriptTyped.find("--");
					if (std::string::npos == sepPos) {
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ", error, couldn't find the transcript in " << transcriptTyped << " of " << callName << " in " << typedFnp << "\n";
						throw std::runtime_error{ss.str()};
					}
					const auto transcript = transcriptTyped.substr(0, sepPos);
					const auto & codons = getCodons(transcript);
					for (const auto & aaTyped : tokenizeString(transcriptTyped.substr(sepPos + 2), ":")) {
						//the reference amino acid, its 1-based position, and the haplotype's amino acid if it's different
						const auto numEnd = aaTyped.find_first_not_of("0123456789", 1);
						if (aaTyped.size() < 2 || 1 == numEnd) {
							std::stringstream ss;
							ss << __PRETTY_FUNCTION__ << ", error, couldn't read amino acid " << aaTyped << " of " << callName << " in " << typedFnp << "\n";
							throw std::runtime_error{ss.str()};
						}
						const char refAA = aaTyped.front();
						const uint32_t pos = njh::StrToNumConverter::stoToNum<uint32_t>(aaTyped.substr(1, numEnd - 1));
						const char altAA = std::string::npos == numEnd ? refAA : aaTyped[numEnd];

						pmo::ProteinVariant variant;
						variant.protein_location_ = genLocation(transcript, pos, pos, std::string(1, refAA), std::string(1, altAA), std::nullopt);
						variant.gene_name_ = codons.geneId_;
						if (!codons.altGeneName_.empty()) {
							variant.alternative_gene_name_ = codons.altGeneName_;
						}
						auto codonSearch = codons.basesByZeroBasedAAPos_.find(pos - 1);
						if (codons.basesByZeroBasedAAPos_.end() == codonSearch || 3 != codonSearch->second.size()) {
							std::stringstream ss;
							ss << __PRETTY_FUNCTION__ << ", error, couldn't find the codon for amino acid " << pos << " of " << transcript << " in " << njh::files::make_path(callsDir, "geneInfos") << "\n";
							throw std::runtime_error{ss.str()};
						}
						std::string refCodon;
						std::string snpCodon;
						uint32_t codonStart = std::numeric_limits<uint32_t>::max();
						uint32_t codonEnd = 0;
						for (const auto & base : codonSearch->second) {
							refCodon.push_back(base.base_);
							auto snpSearch = snpsByPos.find(base.gPos_);
							if (snpsByPos.end() == snpSearch) {
								snpCodon.push_back(base.base_);
							} else {
								snpCodon.push_back('-' == codons.strand_ ? seqUtil::complement(std::string(1, snpSearch->second), "DNA").front() : snpSearch->second);
							}
							codonStart = std::min(codonStart, base.gPos_);
							codonEnd = std::max(codonEnd, base.gPos_ + 1);
						}
						std::optional<std::string> altCodon;
						if ('-' == altAA) {
							altCodon = "---";
						} else {
							const bool indelInCodon = std::any_of(indelSpans.begin(), indelSpans.end(), [&codonStart, &codonEnd](const std::pair<uint32_t, uint32_t> & span) {
								return span.first < codonEnd && span.second > codonStart;
							});
							if (!indelInCodon && translate(snpCodon) == altAA) {
								altCodon = snpCodon;
							} else if (njh::in(callName, translatedCodons) && njh::in(transcript, translatedCodons.at(callName)) &&
									njh::in(pos, translatedCodons.at(callName).at(transcript))) {
								altCodon = translatedCodons.at(callName).at(transcript).at(pos);
							} else {
								++ret.codonsWithoutAlt_;
							}
						}
						variant.codon_genomic_location_ = genLocation(codons.chrom_, codonStart, codonEnd, refCodon, altCodon, codons.strand_);
						proteinVariants.emplace_back(variant);
					}
				}
				ret.proteinVariants_ += proteinVariants.size();
				if (!proteinVariants.empty()) {
					mhap.associated_protein_variants_ = proteinVariants;
				}
			}
			if (seqVariantsByName.end() != seqVarsSearch) {
				auto seqVariants = seqVarsSearch->second;
				std::stable_sort(seqVariants.begin(), seqVariants.end(), [](const pmo::GenomicLocation & a, const pmo::GenomicLocation & b) {
					return std::tie(a.chrom_, a.start_, a.end_) < std::tie(b.chrom_, b.start_, b.end_);
				});
				ret.seqVariants_ += seqVariants.size();
				if (!seqVariants.empty()) {
					mhap.associated_seq_variants_ = seqVariants;
				}
			}
		}
		reps.validate();
	}
	return ret;
}

}  // namespace

int PMOUtilsRunner::build_pmo_from_allele_table(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path alleleTableFnp;
	AlleleTableCols cols;
	bfs::path primersFnp;
	std::string panelName;
	pmo::AddTargetGenomeInfoPars genomePars;
	bfs::path metaFnp;
	std::string metaSampleCol = "sample";
	bfs::path variantCallsDir;
	std::string variantCallsSubDir = "variantCalling";
	std::string variantCallsGenome;
	std::string programVersion = "1";
	OutOptions outOpts(bfs::path(""), ".json.gz");
	seqSetUp setUp(inputCommands);
	setUp.description_ = "Build a PMO from an allele table (one row per haplotype per sample per target) and the panel's primers, "
			"optionally adding each target's genome location and genes (from SeekDeep genTargetInfoFromGenomes), sample meta, "
			"and variant calls already made per target (from SeekDeep variantCallOnSeqAndProtein or processClusters)";
	setUp.examples_.emplace_back("MASTERPROGRAM SUBPROGRAM --alleleTable allele_data.tsv --primers primers.tsv --panelName PfPHAST --out PMO");
	setUp.examples_.emplace_back("MASTERPROGRAM SUBPROGRAM --alleleTable allele_data.tsv --targetCol target_name --primers PfPHAST_Pf3D7_primers.tsv --panelName PfPHAST "
			"--locations PfGenomes_primerLocs_2Errors/locationsByGenome --genomeDir genomes --gffDir info/gff --preferredGenomes Pf3D7 "
			"--meta meta.tsv --variantCallsDir allele_data_variantCalling --out PMO");
	setUp.setOption(alleleTableFnp, "--alleleTable", "Tab delimited table with a header, one row per haplotype per sample per target, e.g. the allele_data.txt of MAD4HatTeR", true, "Allele Table");
	setUp.setOption(cols.sample_, "--sampleCol", "The column with the library sample name", false, "Allele Table");
	setUp.setOption(cols.target_, "--targetCol", "The column with the target name, these need to match the targets of --primers", false, "Allele Table");
	setUp.setOption(cols.seq_, "--seqCol", "The column with the haplotype sequence", false, "Allele Table");
	setUp.setOption(cols.reads_, "--readsCol", "The column with the haplotype's reads in the sample", false, "Allele Table");
	setUp.setOption(cols.mhapName_, "--mhapNameCol", "Optional column with a name for each haplotype, needs to be the same for the same sequence of a target", false, "Allele Table");
	setUp.setOption(cols.specimen_, "--specimenCol", "Optional column with the specimen of each library sample, defaults to the library sample name", false, "Allele Table");
	setUp.setOption(primersFnp, "--primers", "Table with a header of the panel's targets, the first three columns are target, forward primer and reverse primer, every target in it becomes a target of the panel", true, "Panel");
	setUp.setOption(panelName, "--panelName", "Name of the panel", true, "Panel");
	setUp.setOption(genomePars.locationsDir_, "--locations", "Optional, the locationsByGenome directory of SeekDeep genTargetInfoFromGenomes, to add each target's genome, insert and primer locations and genes", false, "Target Locations");
	setUp.setOption(genomePars.gffDir_, "--gffDir", "Directory with <genome>.gff files for the genes the targets overlap, used with --locations", false, "Target Locations");
	setUp.setOption(genomePars.genomeDir_, "--genomeDir", "Directory with <genome>.2bit files for the reference sequence of each target, used with --locations", false, "Target Locations");
	setUp.setOption(genomePars.preferredGenomes_, "--preferredGenomes", "Genomes to use first, in order, for targets found in several genomes, used with --locations", false, "Target Locations");
	setUp.setOption(metaFnp, "--meta", "Optional tab delimited table with a header of sample meta, added as specimen fields, columns that aren't in the PMO schema become custom fields", false, "Meta");
	setUp.setOption(metaSampleCol, "--metaSampleCol", "The column of --meta with the library sample or specimen name", false, "Meta");
	setUp.setOption(variantCallsDir, "--variantCallsDir", "Optional directory with the variant calls per target, <this>/<target>/<--variantCallsSubDir>, e.g. the --dout of SeekDeep variantCallOnSeqAndProtein", false, "Variant Calls");
	setUp.setOption(variantCallsSubDir, "--variantCallsSubDir", "The directory within each target's directory of --variantCallsDir with its calls", false, "Variant Calls");
	setUp.setOption(variantCallsGenome, "--variantCallsGenome", "Name of the genome the variants were called against, defaults to the only genome from --locations", false, "Variant Calls");
	setUp.setOption(programVersion, "--programVersion", "Version of this program to record in the PMO");
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	njh::files::checkExistenceThrow(std::vector<bfs::path>{alleleTableFnp, primersFnp}, __PRETTY_FUNCTION__);

	pmo::PmoGenerationMethod generationMethod;
	generationMethod.program_name_ = njh::pasteAsStr(setUp.commands_.masterProgram_, " build_pmo_from_allele_table");
	generationMethod.program_version_ = programVersion;
	auto pmoObj = genPmoFromAlleleTable(alleleTableFnp, cols, readPrimers(primersFnp), panelName, generationMethod);
	std::cout << "targets: " << pmoObj.target_info_.size() << ", with haplotypes: " << pmoObj.representative_microhaplotypes_.targets_.size() << std::endl;
	std::cout << "library samples: " << pmoObj.library_sample_info_.size() << ", specimens: " << pmoObj.specimen_info_.size() << std::endl;

	if (!genomePars.locationsDir_.empty()) {
		auto results = pmo::addTargetGenomeInfo(pmoObj, genomePars);
		for (const auto & genome : results.targetsPerGenome_) {
			std::cout << "targets in " << genome.first << ": " << genome.second << std::endl;
		}
		std::cout << "targets with genes: " << results.targetsWithGenes_ << std::endl;
		if (!results.targetsNotFound_.empty()) {
			std::cerr << results.targetsNotFound_.size() << " targets aren't in any genome in " << genomePars.locationsDir_ << ": " << njh::conToStr(results.targetsNotFound_, ",") << std::endl;
		}
	}
	if (!metaFnp.empty()) {
		table meta(metaFnp.string(), "\t", true);
		auto results = pmo::addSpecimenMeta(pmoObj, meta, metaSampleCol);
		std::cout << "specimens with meta: " << results.specimensUpdated_ << std::endl;
		if (!results.metaSamplesNotInPmo_.empty()) {
			std::cerr << results.metaSamplesNotInPmo_.size() << " samples in " << metaFnp << " aren't in the PMO: " << njh::conToStr(results.metaSamplesNotInPmo_, ",") << std::endl;
		}
	}
	if (!variantCallsDir.empty()) {
		if (variantCallsGenome.empty()) {
			if (!pmoObj.targeted_genomes_.has_value() || 1 != pmoObj.targeted_genomes_->size()) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error, --variantCallsGenome is needed when --locations doesn't give exactly one genome\n";
				throw std::runtime_error{ss.str()};
			}
			variantCallsGenome = pmoObj.targeted_genomes_->front().name_;
		}
		auto results = addVariantCallsFromDir(pmoObj, variantCallsDir, variantCallsSubDir, variantCallsGenome);
		std::cout << "targets with variant calls: " << results.targetsWithCalls_ << std::endl;
		std::cout << "haplotypes with variant calls: " << results.mhapsWithCalls_ << ", without: " << results.mhapsWithoutCalls_ << std::endl;
		std::cout << "protein variants: " << results.proteinVariants_ << ", sequence variants: " << results.seqVariants_ << std::endl;
		if (!results.targetsWithoutCalls_.empty()) {
			std::cerr << results.targetsWithoutCalls_.size() << " targets have no variant calls in " << variantCallsDir << ": " << njh::conToStr(results.targetsWithoutCalls_, ",") << std::endl;
		}
		if (results.codonsWithoutAlt_ > 0) {
			std::cerr << results.codonsWithoutAlt_ << " protein variants have no haplotype codon, it couldn't be worked out from the tables" << std::endl;
		}
	}
	pmoObj.validate();
	pmo::writePmo(pmoObj, outOpts);
	return 0;
}

}  // namespace njhseq
