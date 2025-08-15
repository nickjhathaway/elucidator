//
// Created by Nicholas Hathaway on 8/7/25.
//



#include "seqUtilsModRunner.hpp"
#include <njhseq/IO/SeqIO/SeqIO.hpp>
#include <njhseq/objects/Gene/CodonSampler.hpp>


namespace njhseq {

int seqUtilsModRunner::generateNucleotidePossibleFromProteins(const njh::progutils::CmdArgs & inputCommands) {

	VecStr restriction_sites_to_remove {"GAATTC", "CTCGAG", "AAGCCT"};
	// ecor1 = "GAATTC"
	// xho1 = "CTCGAG"
	// hindIII = "AAGCCT"
	bool use_e_coli_optimized_codons = false;
	bool do_not_set_optimal_stop_codon = false;

	uint64_t seed = std::numeric_limits<uint64_t>::max();
	seqSetUp setUp(inputCommands);
	setUp.description_ = "Generate nucleotide sequences from input protein(s)";
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processDefaultReader(seqSetUp::singleInFormatsAvailable_, true);
	setUp.setOption(use_e_coli_optimized_codons, "--use_e_coli_optimized_codons", "use the e coli codon usage table");
	setUp.setOption(do_not_set_optimal_stop_codon, "--do_not_set_optimal_stop_codon", "use a random selection of stop codons rather than only the AMBER optimal stop codon: TAG");
	setUp.setOption(seed, "--seed", "seed for random generators");
	setUp.setOption(restriction_sites_to_remove, "--restriction_sites_to_remove", "restriction_sites_to_remove");

	setUp.finishSetUp(std::cout);

	VecStr restriction_sites_to_remove_rev_comp;
	for (const auto & site : restriction_sites_to_remove) {
		auto rev_comp_site = seqUtil::reverseComplement(site, "DNA");
		if (rev_comp_site != site) {
			restriction_sites_to_remove_rev_comp.emplace_back(rev_comp_site);
		}
	}
	std::unique_ptr<aminoAcidInfo::CodonSampler> sampler;
	//@todo allow the supplying of a custom usage codon map
	if (use_e_coli_optimized_codons) {
		sampler = std::make_unique<aminoAcidInfo::CodonSampler>(aminoAcidInfo::infos::e_coli_dna_codon_usage);
	} else {
		//uniform codon usage
		sampler = std::make_unique<aminoAcidInfo::CodonSampler>();
	}

	if (!do_not_set_optimal_stop_codon) {
		sampler->set_stop_codon_usage_to_amber_only(true);
	}
	if (std::numeric_limits<uint64_t>::max() != seed) {
		sampler->set_seed(seed);
	}
	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	reader.openOut();
	seqInfo seq;
	auto has_restriction_sites = [&restriction_sites_to_remove](const std::string & out_nuc) {
		bool has_sites = false;
		for (const auto & site : restriction_sites_to_remove) {
			if (std::string::npos != out_nuc.find(site)) {
				has_sites = true;
				break;
			}
		}
		return has_sites;
	};

	auto has_restriction_sites_rev_comp = [&restriction_sites_to_remove_rev_comp](const std::string & out_nuc) {
		bool has_sites = false;
		for (const auto & site : restriction_sites_to_remove_rev_comp) {
			if (std::string::npos != out_nuc.find(site)) {
				has_sites = true;
				break;
			}
		}
		return has_sites;
	};



	auto has_restriction_sites_both_directions_checked = [
		&has_restriction_sites,
		&has_restriction_sites_rev_comp](const std::string & out_nuc) {
		bool has_sites = has_restriction_sites(out_nuc);
		if (!has_sites) {
			//not all sites are palindromes, so search rev comp too for the non-pallinodrome (if any)
			has_sites =  has_restriction_sites_rev_comp(out_nuc);
		}
		return has_sites;
	};

	auto recode_to_remove_sites = [
		// &seed,
		&sampler](seqInfo & out_nuc, const VecStr & restriction_sites_to_remove ) {
		//assumes the input sequence is in frame and is a full reading frame(divisble by 3)
		//not all restriction sites are palidromes so have to check both forward and reverse directions
		for (const auto & site : restriction_sites_to_remove) {
			auto site_positions = findOccurences(out_nuc.seq_, site);
			njh::sort(site_positions);
			uint32_t site_size = site.size();
			for (const auto & pos : site_positions) {
				if (out_nuc.seq_.substr(pos, site.size()) != site) {
					//it's possible that site was changed if multiple occurrences were found
					continue;
				}
				auto end = pos + site.size();
				auto new_start = pos - (pos % 3);
				auto new_end =  (end % 3 == 0) ? end : end + (3 - (end % 3));
				std::string sub_str_to_modify = out_nuc.seq_.substr(new_start, new_end - new_start);
				std::vector<uint32_t> positons_to_attempt_to_recode;
				for (uint32_t i = 0; i < sub_str_to_modify.size(); i += 3) {
					positons_to_attempt_to_recode.emplace_back(i);
				}
				njh::sort(positons_to_attempt_to_recode, [&pos, &new_start,&site_size](uint32_t i1, uint32_t i2) {
					auto i1_from_start = uAbsdiff(i1, pos - new_start);
					auto i1_from_end = uAbsdiff(i1 + 3, pos  - new_start + site_size);
					auto i2_from_start = uAbsdiff(i2, pos - new_start);
					auto I2_from_end = uAbsdiff(i2 + 3, pos  - new_start + site_size);
					if (i1_from_start + i1_from_end == i2_from_start + I2_from_end) {
						return uAbsdiff(i1_from_start, i1_from_end) < uAbsdiff(i2_from_start, I2_from_end);
					}
					return i1_from_start + i1_from_end <= i2_from_start + I2_from_end;
				});
				bool all_sites_no_alts = true;
				std::vector<uint32_t> filt_positons_to_attempt_to_recode;
				for (auto change_pos : positons_to_attempt_to_recode) {
					auto codon = sub_str_to_modify.substr(change_pos, 3);
					auto alts = sampler->get_alt_codons(codon);
					if (!alts.empty()) {
						all_sites_no_alts = false;
						filt_positons_to_attempt_to_recode.emplace_back(change_pos);
					}
				}
				if (all_sites_no_alts) {
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << ", while recoding to remove restriction site, no alt codons for all positions within " << sub_str_to_modify << "\n";
					throw std::runtime_error{ss.str()};
				}

				for (const auto & change_pos : filt_positons_to_attempt_to_recode) {
					auto codon = sub_str_to_modify.substr(change_pos, 3);
					auto alts = sampler->get_alt_codons(codon);
					bool good_change = false;
					for (const auto & alt : alts) {
						sub_str_to_modify.replace(change_pos, alt.size(), alt);
						if (std::string::npos == sub_str_to_modify.find(alt)) {
							good_change = true;
							break;
						}
					}
					if (good_change) {
						break;
					}
				}
				//
				// auto change_pos = filt_positons_to_attempt_to_recode.front();
				// if (filt_positons_to_attempt_to_recode.size() > 1) {
				// 	std::vector<uint32_t> filt_positons_to_attempt_to_recode_weights;
				// 	for (const auto & filt_pos : iter::range(filt_positons_to_attempt_to_recode.size())) {
				// 		filt_positons_to_attempt_to_recode_weights.emplace_back(filt_positons_to_attempt_to_recode.size() + 1 - filt_pos);
				// 	}
				// 	njh::randObjectGen pos_gen(filt_positons_to_attempt_to_recode, filt_positons_to_attempt_to_recode_weights);
				// 	if (std::numeric_limits<uint64_t>::max() != seed) {
				// 		pos_gen.set_seed(seed);
				// 	}
				// 	change_pos = pos_gen.genObj();
				// }
				// auto codon = sub_str_to_modify.substr(change_pos, 3);
				// auto alts = sampler->get_alt_codons(codon);
				// if (alts.size() == 1) {
				// 	sub_str_to_modify.replace(change_pos, alts.front().size(), alts.front());
				// } else {
				// 	std::vector<double> alts_weights;
				// 	for (const auto & alt : alts) {
				// 		alts_weights.emplace_back(sampler->codon_usage_[alt]);
				// 	}
				// 	njh::randObjectGen alt_gen(alts, alts_weights);
				// 	if (std::numeric_limits<uint64_t>::max() != seed) {
				// 		alt_gen.set_seed(seed);
				// 	}
				// 	auto randon_alt = alt_gen.genObj();
				// 	sub_str_to_modify.replace(change_pos, randon_alt.size(), randon_alt);
				// }
				out_nuc.seq_.replace(new_start, new_end - new_start, sub_str_to_modify);
			}
		}
	};

	while (reader.readNextRead(seq)) {
		seqInfo out_seq(seq.name_, sampler->gen(seq.seq_));
		uint32_t attempts = 0;
		while (has_restriction_sites_both_directions_checked(out_seq.seq_)) {
			++attempts;
			if (attempts >= 5) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << ", error " << " attempted 5 times to recode to remove restriction sites and failed for\nnuc:"
				<< out_seq.seq_ << "\n"
				<< "protein: " << seq.seq_ << "\n";
				throw std::runtime_error{ss.str()};
			}
			recode_to_remove_sites(out_seq, restriction_sites_to_remove);
			if (has_restriction_sites_rev_comp(out_seq.seq_) ) {
				recode_to_remove_sites(out_seq, restriction_sites_to_remove_rev_comp);
			}
		}
		if (out_seq.translateRet(false, false).seq_ != seq.seq_) {
			std::stringstream ss;
			auto new_trans = out_seq.translateRet(false, false);
			ss << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << ", error " << " output translated seq:\n" << new_trans.seq_ << "\n";
			ss << "doesn't match input seq:\n" << seq.seq_ << "\n";
			for (const auto pos : iter::range(std::min(new_trans.seq_.size(), seq.seq_.size()))) {
				if (new_trans.seq_[pos] != seq.seq_[pos]) {
					ss << "at pos: " << pos << " input seq: " << seq.seq_[pos] << " new seq: " << new_trans.seq_[pos] << "\n";
				}
			}
			throw std::runtime_error{ss.str()};
		}
		reader.write(out_seq);
	}

	return 0;
}

int seqUtilsModRunner::generateAllPossibleNucleotidePossibleFromProtein(const njh::progutils::CmdArgs & inputCommands) {
	VecStr input_proteins;
	OutOptions outOpts;
	seqSetUp setUp(inputCommands);
	setUp.description_ = "Generate all possible nucleotide sequences possible for input protein(s) using standard codon libary";
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(input_proteins, "--input_proteins", "the input proteins to convert");
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	OutputStream out(outOpts);

	// convert to upper and check for nonstandard amino acids
	std::set<char> non_standard_aas;
	for (auto &protein: input_proteins) {
		njh::strToUpper(protein);
		for (const auto aa : protein) {
			if (njh::notIn(aa, aminoAcidInfo::infos::allInfo)) {
				non_standard_aas.emplace(aa);
			}
		}
	}
	if (!non_standard_aas.empty()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << ", error " << " found the following non-standard amino acids " << "\n";
		ss << njh::conToStr(non_standard_aas, ",") << "\n";
		throw std::runtime_error{ss.str()};
	}

	for (const auto &protein: input_proteins) {
		VecStr output;
		for (const auto aa : protein) {
			if (output.empty()) {
				for (const auto & dna_codon : aminoAcidInfo::infos::allInfo.at(aa).dnaCodons_) {
					output.emplace_back(dna_codon);
				}
			} else {
				VecStr new_output;
				for (auto current_out : output) {
					for (const auto & dna_codon : aminoAcidInfo::infos::allInfo.at(aa).dnaCodons_) {
						new_output.emplace_back(current_out + dna_codon);
					}
				}
				output.swap(new_output);
			}
		}
		out << njh::conToStr(output, "\n") << std::endl;
	}

	return 0;
}


} //namespace njhseq



