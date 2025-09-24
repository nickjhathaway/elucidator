//
// Created by Nicholas Hathaway on 8/20/25.
//

#include "pmo_utils_runner.hpp"


#include <njhseq/IO/SeqIO/SeqIO.hpp>
#include <njhseq/objects/BioDataObject/pmo.h>
#include <njhseq/objects/dataContainers/tables/TableReader.hpp>
#include <nlohmann/json.hpp>

namespace njhseq {
PMOUtilsRunner::PMOUtilsRunner()
        : njh::progutils::ProgramRunner(
        {
        	addFunc("read_pmo", read_pmo, false),
        	addFunc("get_overlap_between_panels_in_pmos", get_overlap_between_panels_in_pmos, false),
        	addFunc("add_protein_variant_info_to_pmo", add_protein_variant_info_to_pmo, false),
        	addFunc("count_protein_variant_info_to_pmo", count_protein_variant_info_to_pmo, false),

        },//
        "PMOUtils") {}


int PMOUtilsRunner::read_pmo(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path pmo_fnp;
	OutOptions out_options;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(pmo_fnp, "--pmo", "pmo file to read", true);
	setUp.processWritingOptions(out_options);
	setUp.finishSetUp(std::cout);

	InputStream in(pmo_fnp);
	auto pmo_obj= pmo::PortableMicrohaplotypeObject::from_json(nlohmann::json::parse(in));
	pmo_obj.validate();
	OutputStream out(out_options);
	out << pmo_obj.to_json();

	return 0;
}

template <class T>
T& ensure_opt_vec(std::optional<T>& opt) {
	if (!opt) opt.emplace();
	return *opt;
}

int PMOUtilsRunner::add_protein_variant_info_to_pmo(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path pmo_fnp;

	//required
	//target_name, seq
	//transcript, amino acid position, ref_aa, seq_aa
	//optional
	//gene_name, alt_gene_name
	//codon info, chrom, start, end, strand, ref_codon, seq_codon

	std::string target_name_col = "target_name";
	std::string seq_col = "seq";
	std::string transcript_col = "transcript";
	std::string aa_position_col = "aa_position";
	std::string ref_aa_col = "ref_aa";
	std::string seq_aa_col = "aa";

	std::string gene_alt_name;
	std::string gene_name;

	uint32_t genome_id = 0;
	OutOptions out_options;
	bfs::path protein_variant_info_fnp;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(pmo_fnp, "--pmo", "pmo file to read", true);
	setUp.setOption(protein_variant_info_fnp, "--protein_variant_info_fnp", "protein variant info fnp", true);

	setUp.setOption(gene_alt_name, "--gene_alt_name", "optional gene_alt_name");
	setUp.setOption(gene_name, "--gene_name", "optional gene_name");

	setUp.setOption(genome_id, "--genome_id", "genome_id that the variant refer to");

	setUp.setOption(target_name_col, "--target_name_col", "target_name_col");
	setUp.setOption(seq_col, "--seq_col", "seq_col");
	setUp.setOption(transcript_col, "--transcript_col", "transcript_col");
	setUp.setOption(aa_position_col, "--aa_position_col", "aa_position_col");
	setUp.setOption(ref_aa_col, "--ref_aa_col", "ref_aa_col");
	setUp.setOption(seq_aa_col, "--seq_aa_col", "seq_aa_col");


	setUp.processWritingOptions(out_options);
	setUp.finishSetUp(std::cout);

	VecStr required_columns{target_name_col, seq_col, transcript_col, aa_position_col, ref_aa_col, seq_aa_col};

	InputStream in(pmo_fnp);
	auto pmo_obj= pmo::PortableMicrohaplotypeObject::from_json(nlohmann::json::parse(in));
	pmo_obj.validate();
	if (genome_id >= pmo_obj.targeted_genomes_.size()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << ", error " << " genome_id: " << genome_id << " out of range of targeted genomes, size: " <<  pmo_obj.targeted_genomes_.size()<< "\n";
		throw std::runtime_error{ss.str()};
	}

	OutputStream out(out_options);

	//build a table of seqs that are present in PMO representative haps
	std::unordered_map<std::string, std::unordered_set<std::string>> seqs_for_targets;
	std::unordered_map<std::string, uint32_t> target_key;
	std::unordered_map<std::string, std::unordered_map<std::string,uint32_t> > target_seq_key;

	for (const auto & target_seqs : iter::enumerate(pmo_obj.representative_microhaplotypes_.targets_)) {
		auto target_name = pmo_obj.target_info_[target_seqs.element.target_id_].target_name_;
		target_key[target_name] = target_seqs.index;
		for (const auto & haps : iter::enumerate(target_seqs.element.microhaplotypes_)) {
			target_seq_key[target_name][haps.element.seq_] = haps.index;
			seqs_for_targets[target_name].insert(haps.element.seq_);
		}
	}

	TableReader reader(TableIOOpts::genTabFileIn(protein_variant_info_fnp));
	reader.header_.checkForColumnsThrow(required_columns, __PRETTY_FUNCTION__);
	uint32_t target_name_col_pos = reader.header_.getColPos(target_name_col);
	uint32_t seq_col_pos = reader.header_.getColPos(seq_col);
	uint32_t transcript_col_pos = reader.header_.getColPos(transcript_col);
	uint32_t aa_position_col_pos = reader.header_.getColPos(aa_position_col);
	uint32_t ref_aa_col_pos = reader.header_.getColPos(ref_aa_col);
	uint32_t seq_aa_col_pos = reader.header_.getColPos(seq_aa_col);
	VecStr row;
	//required
			//target_name, seq
			//transcript, amino acid position, ref_aa, seq_aa
	//optional
			//gene_name, alt_gene_name
			//codon info, chrom, start, end, strand, ref_codon, seq_codon

	while (reader.getNextRow(row)) {
		const auto & target_name = row[target_name_col_pos];
		const auto & seq = row[seq_col_pos];
		if (njh::notIn(target_name, seqs_for_targets)) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << ", error " << " target: "<< target_name << " not found " << "\n";
			ss << "options are: " << njh::conToStr(njh::getSetOfMapKeys(seqs_for_targets), ",") << "\n";
			throw std::runtime_error{ss.str()};
		}
		if (njh::notIn(seq, seqs_for_targets[target_name])) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << ", error " << " seq: "<< seq << " not found " << "\n";
			ss << "options are: " << njh::conToStr(seqs_for_targets[target_name], ",") << "\n";
			throw std::runtime_error{ss.str()};
		}
		pmo::GenomicLocation variant_location;
		variant_location.chrom_ = row[transcript_col_pos];
		variant_location.start_ = njh::StrToNumConverter::stoToNum<uint32_t>(row[aa_position_col_pos]);
		variant_location.end_ = njh::StrToNumConverter::stoToNum<uint32_t>(row[aa_position_col_pos]);//keep it one based sigh
		variant_location.genome_id_ = genome_id;
		variant_location.alt_seq_ = row[seq_aa_col_pos];
		variant_location.ref_seq_ = row[ref_aa_col_pos];
		pmo::ProteinVariant variant;
		variant.protein_location_  = variant_location;
		if (!gene_name.empty() && njh::in(gene_name, reader.header_.columnNames_)) {
			variant.gene_name_ = row[reader.header_.getColPos(gene_name)];
		}
		if (!gene_alt_name.empty() && njh::in(gene_alt_name, reader.header_.columnNames_)) {
			variant.alternative_gene_name_ = row[reader.header_.getColPos(gene_alt_name)];
		}
		ensure_opt_vec(pmo_obj.representative_microhaplotypes_.targets_[target_key[target_name]].microhaplotypes_[target_seq_key[target_name][seq]].associated_protein_variants_).emplace_back(variant);
	}

	pmo_obj.validate();
	out << pmo_obj.to_json();

	return 0;
}



int PMOUtilsRunner::count_protein_variant_info_to_pmo(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path pmo_fnp;

	VecStr specimen_meta;
	std::string transcript_id;
	std::vector<uint32_t> aa_positions;
	bfs::path transcript_table_fnp;
	OutOptions out_options;

	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(pmo_fnp, "--pmo", "pmo file to read", true);
	setUp.setOption(specimen_meta, "--specimen_meta", "specimen_meta");

	setUp.setOption(transcript_id, "--transcript_id", "transcript_id");
	setUp.setOption(aa_positions, "--aa_positions", "aa_positions", !transcript_id.empty());

	setUp.setOption(transcript_table_fnp, "--transcript_table_fnp",
	                "a table of transcripts and positions to count, columns transcript_id and aa_positions need to exist, tab delimited",
	                transcript_id.empty());

	setUp.processWritingOptions(out_options);
	setUp.finishSetUp(std::cout);

	std::map<std::string, std::vector<uint32_t> > aa_positions_for_transcript_ids;
	if (transcript_table_fnp.empty()) {
		aa_positions_for_transcript_ids[transcript_id] = aa_positions;
	} else {
		table transcript_table(transcript_table_fnp, "\t", true);
		transcript_table.checkForColumnsThrow({"transcript_id", "aa_positions"}, __PRETTY_FUNCTION__);
		for (const auto &row: transcript_table) {
			addOtherVec(aa_positions_for_transcript_ids[row[transcript_table.getColPos("transcript_id")]],
			            vecStrToVecNum<uint32_t>(njh::tokenizeString(row[transcript_table.getColPos("aa_positions")], ","))
			);
		}
	}


	InputStream in(pmo_fnp);
	auto pmo_obj= pmo::PortableMicrohaplotypeObject::from_json(nlohmann::json::parse(in));
	pmo_obj.validate();
	OutputStream out(out_options);
	out << "meta_field"
		<< "\t" << "meta_value"
		<< "\t" << "transcript_id"
		<< "\t" << "aa_position"
		<< "\t" << "seq_aa"
		<< "\t" << "allele_count"
		<< "\t" << "allele_freq"
		<< "\t" << "allele_total"
		<< "\t" << "specimen_count"
		<< "\t" << "specimen_prev"
		<< "\t" << "specimen_total" << std::endl;
	std::unordered_set<std::string> targets_with_transcripts;

	for (const auto & target : pmo_obj.representative_microhaplotypes_.targets_) {
		for (const auto & microhap : target.microhaplotypes_) {
			if (microhap.associated_protein_variants_.has_value()) {
				for (const auto & protein_variant : *microhap.associated_protein_variants_) {
					if (njh::in(protein_variant.protein_location_.chrom_, aa_positions_for_transcript_ids) &&
						njh::in(protein_variant.protein_location_.start_, aa_positions_for_transcript_ids[protein_variant.protein_location_.chrom_])) {
						targets_with_transcripts.emplace(pmo_obj.target_info_[target.target_id_].target_name_);
					}
				}
			}
		}
	}

	//check meta
	std::unordered_map<uint32_t, std::unordered_map<std::string, std::string>> spec_metas;
	if (!specimen_meta.empty()) {
		for (const auto & meta : specimen_meta) {
			if (meta == "collection_country") {
				for (const auto & specimen : iter::enumerate(pmo_obj.specimen_info_)) {
					spec_metas[specimen.index][meta] = specimen.element.collection_country_;
				}
			} else if (meta == "collection_year") {
				for (const auto & specimen : iter::enumerate(pmo_obj.specimen_info_)) {
					auto pos = specimen.element.collection_date_.find('-');
					std::string year = pos != std::string::npos ? specimen.element.collection_date_.substr(0, pos) : specimen.element.collection_date_;
					spec_metas[specimen.index][meta] = year;
				}
			} else if (meta == "collection_country::collection_year") {
				for (const auto & specimen : iter::enumerate(pmo_obj.specimen_info_)) {
					auto pos = specimen.element.collection_date_.find('-');
					std::string year = pos != std::string::npos ? specimen.element.collection_date_.substr(0, pos) : specimen.element.collection_date_;
					spec_metas[specimen.index][meta] = njh::pasteAsStr(specimen.element.collection_country_, "-", year);
				}
			}
		}
	}



	if (!targets_with_transcripts.empty()) {
		//key1 = meta_field, key2 = meta_value, key3 = transcript, key4 = aa_position, key5 = seq_aa, value = allele_count
		//key1 = meta_field, key2 = meta_value, key3 = transcript, key4 = aa_position, key5 = seq_aa, value = specimen_ids
		std::map<std::string, std::map<std::string, std::map<std::string, std::map<uint32_t, std::map<std::string, uint32_t>>>>> allele_counts;
		std::map<std::string, std::map<std::string, std::map<std::string, std::map<uint32_t, std::map<std::string, std::unordered_set<uint32_t>>>>>> specimens;

		std::vector<uint32_t> mhaps_ids_for_targets;
		for (const auto & mhaps : iter::enumerate(pmo_obj.representative_microhaplotypes_.targets_)) {
			if (njh::in(pmo_obj.target_info_[mhaps.element.target_id_].target_name_, targets_with_transcripts)) {
				mhaps_ids_for_targets.emplace_back(mhaps.index);
			}
		}
		for (const auto & detected : pmo_obj.detected_microhaplotypes_) {
			for (const auto & library_sample : detected.library_samples_) {
				for (const auto & target_haps : library_sample.target_results_) {
					if (njh::in(target_haps.mhaps_target_id_, mhaps_ids_for_targets)) {
						for (const auto & hap : target_haps.mhaps_) {
							if (pmo_obj.representative_microhaplotypes_.targets_[target_haps.mhaps_target_id_].microhaplotypes_[hap.mhap_id_].associated_protein_variants_.has_value()) {
								for (const auto &protein_variant: *pmo_obj.representative_microhaplotypes_.targets_[target_haps.
									     mhaps_target_id_].microhaplotypes_[hap.mhap_id_].associated_protein_variants_) {
									if (njh::in(protein_variant.protein_location_.chrom_, aa_positions_for_transcript_ids) &&
									    njh::in(protein_variant.protein_location_.start_,
									            aa_positions_for_transcript_ids[protein_variant.protein_location_.chrom_])) {
										auto specimen_id = pmo_obj.library_sample_info_[library_sample.library_sample_id_].specimen_id_;
										++allele_counts["all"]["all"][protein_variant.protein_location_.chrom_][protein_variant.protein_location_.start_][protein_variant.protein_location_.alt_seq_.value()];
										specimens["all"]["all"][protein_variant.protein_location_.chrom_][protein_variant.protein_location_.start_][protein_variant.protein_location_.alt_seq_.value()].emplace(specimen_id);
										if (!specimen_meta.empty()) {
											for (const auto & meta_value : spec_metas[specimen_id]) {
												++allele_counts[meta_value.first][meta_value.second][protein_variant.protein_location_.chrom_][protein_variant.protein_location_.start_][protein_variant.protein_location_.alt_seq_.value()];
												specimens[meta_value.first][meta_value.second][protein_variant.protein_location_.chrom_][protein_variant.protein_location_.start_][protein_variant.protein_location_.alt_seq_.value()].emplace(specimen_id);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		for (const auto & meta_field_counts : allele_counts) {
			for (const auto & meta_value : meta_field_counts.second) {
				for (const auto & transcript : meta_value.second) {
					for (const auto & pos : transcript.second) {
						double allele_total = 0;
						std::unordered_set<uint32_t> total_specimens;
						for (const auto & var : pos.second) {
							allele_total += var.second;
							total_specimens.insert(specimens[meta_field_counts.first][meta_value.first][transcript.first][pos.first][var.first].begin(),
							specimens[meta_field_counts.first][meta_value.first][transcript.first][pos.first][var.first].end()
								);
						}
						for (const auto & var : pos.second) {
							out << meta_field_counts.first
									<< "\t" << meta_value.first
									<< "\t" << transcript.first
									<< "\t" << pos.first
									<< "\t" << var.first
									<< "\t" << var.second
									<< "\t" << var.second/allele_total
									<< "\t" << allele_total
									<< "\t" << specimens[meta_field_counts.first][meta_value.first][transcript.first][pos.first][var.first].size()
									<< "\t" << specimens[meta_field_counts.first][meta_value.first][transcript.first][pos.first][var.first].size()/static_cast<double>(total_specimens.size())
									<< "\t" << total_specimens.size() << std::endl;
						}
					}
				}
			}
		}
	}

	//notes,
	//what take into account untyped microhaps,
	//only will work if one library_sample per specimen for the allele counts
	//what take into account overlapping targets
	return 0;
}

int PMOUtilsRunner::get_overlap_between_panels_in_pmos(const njh::progutils::CmdArgs & inputCommands) {
	std::vector<bfs::path> pmo_fnps;
	bool strands_must_match = false;
	OutOptions out_options;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(pmo_fnps, "--pmos", "A list of file name paths to read in", true);
	setUp.setOption(strands_must_match, "--strands_must_match", "strands must match");

	setUp.processWritingOptions(out_options);
	setUp.finishSetUp(std::cout);
	std::vector<pmo::PortableMicrohaplotypeObject> pmo_objs;
	OutputStream out(out_options);

	for (const auto & pmo_fnp : pmo_fnps) {
		InputStream in(pmo_fnp);
		pmo_objs.emplace_back(pmo::PortableMicrohaplotypeObject::from_json(nlohmann::json::parse(in)));
		pmo_objs.back().validate();
	}
	out << "panel1";
	out << "\t" << "panel1_chrom\tpanel1_start\tpanel1_end\tpanel1_name\tpanel1_len\tpanel1_strand";

	out << "\t" << "panel2";
	out << "\t" << "panel2_chrom\tpanel2_start\tpanel2_end\tpanel2_name\tpanel2_len\tpanel2_strand";

	out << "\t" << "overlap_chrom\toverlap_start\toverlap_end\toverlap_name\toverlap_len";

	out << "\t" << "panel1_coverage" << "\t" << "panel2_coverage" << std::endl;
	std::unordered_map<std::string, std::vector<std::shared_ptr<Bed6RecordCore>>> bedsByPanel;
	for (const auto & pmo_obj : pmo_objs) {
		for (const auto & pmo_panel : pmo_obj.panel_info_) {
			if (njh::notIn(pmo_panel.panel_name_, bedsByPanel)) {
				std::unordered_map<uint32_t, VecStr> targets_to_reactions;
				for (const auto & reaction : pmo_panel.reactions_) {
					for (const auto & target : reaction.panel_targets_) {
						targets_to_reactions[target].emplace_back(reaction.reaction_name_);
					}
				}
				std::vector<std::shared_ptr<Bed6RecordCore>> current_beds;
				for (const auto & target : targets_to_reactions) {
					if (!pmo_obj.target_info_[target.first].insert_location_.has_value()) {
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << ", error " << "target: " << pmo_obj.target_info_[target.first].target_name_ << " doesn't have isnert location loaded" << "\n";
						throw std::runtime_error{ss.str()};
					}
					current_beds.emplace_back(std::make_shared<Bed6RecordCore>(pmo_obj.target_info_[target.first].insert_location_.value().chrom_,
						pmo_obj.target_info_[target.first].insert_location_.value().start_,
						pmo_obj.target_info_[target.first].insert_location_.value().end_,
						pmo_obj.target_info_[target.first].target_name_,
						pmo_obj.target_info_[target.first].insert_location_.value().end_ - pmo_obj.target_info_[target.first].insert_location_.value().start_,
						pmo_obj.target_info_[target.first].insert_location_.value().strand_.has_value() ? pmo_obj.target_info_[target.first].insert_location_.value().strand_.value().front() : '+'
						));
					auto reactions_str = njh::pasteAsStr("[reactions=",  njh::conToStr(target.second, ","), "]");
					current_beds.back()->extraFields_.emplace_back(reactions_str);
				}
				bedsByPanel[pmo_panel.panel_name_] = current_beds;
			}
		}
	}

	std::vector<std::string> panelNames = getVectorOfMapKeys(bedsByPanel);

	njh::sort(panelNames);
	for(const auto & currentPanelPos : iter::range(panelNames.size())){
		const auto & current_panel_name = panelNames[currentPanelPos];
		const auto & current_panel_regions = bedsByPanel.at(current_panel_name);
		for(const auto & otherPanelPos : iter::range(currentPanelPos + 1, panelNames.size())){
			const auto & other_panel_regions = bedsByPanel[panelNames[otherPanelPos]];
			for(const auto current_bed_pos : iter::range(current_panel_regions.size())){
				const auto & current_bed = current_panel_regions[current_bed_pos];
				std::vector<std::shared_ptr<Bed6RecordCore>> overlapping_regions_other_panel;
				for(const auto & other_panel_bed : other_panel_regions){
					if(current_bed->overlaps(*other_panel_bed, 1) && (!strands_must_match || current_bed->strand_ == other_panel_bed->strand_)){
						overlapping_regions_other_panel.emplace_back(other_panel_bed);
					}
				}
				for(const auto & overlapping_region_other_panel : overlapping_regions_other_panel){
					out << current_panel_name;
					out << "\t" << current_bed->toDelimStr();
					out << "\t" << panelNames[otherPanelPos];
					out << "\t" << overlapping_region_other_panel->toDelimStr();
					uint32_t overlapStart = std::max(overlapping_region_other_panel->chromStart_, current_bed->chromStart_);
					uint32_t overlapStop = std::min(overlapping_region_other_panel->chromEnd_, current_bed->chromEnd_);
					uint32_t overlapLen = overlapStop - overlapStart;
					out << "\t" << overlapping_region_other_panel->chrom_ << "\t" << overlapStart << "\t" << overlapStop << "\t"
							<< njh::pasteAsStr(overlapping_region_other_panel->chrom_, "-", overlapStart, "-", overlapStop) << "\t" << overlapLen;
					out << "\t" << static_cast<double>(overlapLen)/current_bed->length() << "\t" << static_cast<double>(overlapLen)/overlapping_region_other_panel->length() << std::endl;
				}
			}
		}
	}

	return 0;
}


} //namespace njhseq


