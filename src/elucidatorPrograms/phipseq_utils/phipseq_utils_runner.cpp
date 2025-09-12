//
// Created by Nicholas Hathaway on 8/20/25.
//

#include "phipseq_utils_runner.hpp"

#include <njhseq/IO/SeqIO/SeqIO.hpp>
#include <njhseq/objects/counters/hrCounter.hpp>
#include <njhseq/objects/Gene/CodonSampler.hpp>

namespace njhseq {
PhipSeqUtilsRunner::PhipSeqUtilsRunner()
        : njh::progutils::ProgramRunner(
        {
                addFunc("fragmentSequencesForPhipSeq", fragmentSequencesForPhipSeq, false),
                addFunc("generateAllPossibleNucleotidePossibleFromProtein", generateAllPossibleNucleotidePossibleFromProtein, false),
                addFunc("generateNucleotidePossibleFromProteins", generateNucleotidePossibleFromProteins, false),
        	addFunc("appendRandomBarcode", appendRandomBarcode, false),
        	addFunc("countPossiblePhipSeqRandomBarcodes", countPossiblePhipSeqRandomBarcodes, false),
        	addFunc("markGroupsByHammingDistanceCutOff", markGroupsByHammingDistanceCutOff, false),

        },//
        "PhipSeqUtils") {}


int PhipSeqUtilsRunner::fragmentSequencesForPhipSeq(const njh::progutils::CmdArgs & inputCommands) {
	uint32_t step = 4;
	uint32_t window_size = 16;
	double back_seq_overlap_ratio = 0.75;
	bool doNotModifyName = false;
	seqSetUp setUp(inputCommands);
	setUp.description_ = "fragment input with special considerations for preparation of creation of a phipseq library";
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processDefaultReader(true);
	setUp.setOption(step, "--step", "step size");
	setUp.setOption(window_size, "--window_size", "window size");
	setUp.setOption(doNotModifyName, "--doNotModifyName", "do Not Modify output Name");
	setUp.setOption(back_seq_overlap_ratio, "--back_seq_overlap_ratio", "back seq overlap ratio to allow ");

	setUp.finishSetUp(std::cout);

	SeqIO seq_io(setUp.pars_.ioOptions_);
	seq_io.openIn();
	seq_io.openOut();

	seqInfo seq;
	int64_t expected_overlap = window_size - step;
	while (seq_io.readNextRead(seq)) {
		if (seq.seq_.size() > window_size + step) {
			uint32_t seq_count = 0;
			auto back_seq_pos = seq.seq_.size() - window_size;
			for (uint32_t pos = 0; pos + window_size < seq.seq_.size(); pos+=step) {
				auto end = pos + window_size;
				if (end >= back_seq_pos && setUp.pars_.debug_) {
					std::cout << "pos: " << pos << std::endl;
					std::cout << "back_seq_pos: " << back_seq_pos << std::endl;
					std::cout << "end: " << end  << std::endl;
					std::cout << "expected_overlap: " << expected_overlap << std::endl;
					std::cout << "end - back_seq_pos: " << end - back_seq_pos << std::endl;
					std::cout << "static_cast<long double>(end - back_seq_pos)/expected_overlap: " << expected_overlap/static_cast<long double>(end - back_seq_pos) << std::endl<< std::endl;
				}
				if (end < back_seq_pos || (end > back_seq_pos && expected_overlap/static_cast<long double>(end - back_seq_pos)  > back_seq_overlap_ratio)) {
					auto out_seq = seq.getSubRead(pos, window_size);
					out_seq.name_.append(njh::pasteAsStr("_seq", seq_count));
					seq_io.write(out_seq);
					++seq_count;
				}
			}
			//add back
			{
				auto out_seq = seq.getSubRead(seq.seq_.size() - window_size, window_size);
				out_seq.name_.append(njh::pasteAsStr("_seq", seq_count));
				seq_io.write(out_seq);
				++seq_count;
			}
		}
	}
	return 0;
}

int PhipSeqUtilsRunner::generateAllPossibleNucleotidePossibleFromProtein(const njh::progutils::CmdArgs & inputCommands) {
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
				for (const auto & current_out : output) {
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


int PhipSeqUtilsRunner::generateNucleotidePossibleFromProteins(const njh::progutils::CmdArgs & inputCommands) {

	VecStr restriction_sites_to_remove {"GAATTC", "CTCGAG", "AAGCCT"};
	// ecor1 = "GAATTC"
	// xho1 = "CTCGAG"
	// hindIII = "AAGCCT"

	uint32_t max_homopolymer = 8;

	double max_gc_content_window = 0.70;
	double min_gc_content_window = 0.20;
	uint32_t gc_content_window_size = 25;
	uint32_t gc_content_window_step = 10;

	uint32_t max_seq_attempts = 50;

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

	setUp.setOption(max_seq_attempts, "--max_seq_attempts", "max_seq_attempts");

	setUp.setOption(max_homopolymer, "--max_homopolymer", "max_homopolymer");
	setUp.setOption(max_gc_content_window, "--max_gc_content_window", "max_gc_content_window");
	setUp.setOption(min_gc_content_window, "--min_gc_content_window", "min_gc_content_window");

	setUp.setOption(gc_content_window_size, "--gc_content_window_size", "gc_content_window_size");
	setUp.setOption(gc_content_window_step, "--gc_content_window_step", "gc_content_window_step");

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


	auto check_seq = [&max_homopolymer,
		&max_gc_content_window,
		&min_gc_content_window,
		&gc_content_window_size, &gc_content_window_step](const seqInfo & seq){

		readObject seq_info(seq);
		seq_info.setLetterCount();
		seq_info.counter_.calcGcContent();
		bool pass = seq_info.counter_.gcContent_ < max_gc_content_window && seq_info.counter_.gcContent_ > min_gc_content_window;
		if (pass) {
			seq_info.createCondensedSeq();
			for (const auto & count : seq_info.condensedSeqCount) {
				if (count > max_homopolymer) {
					pass = false;
					break;
				}
			}
		}
		if (pass) {
			if (len(seq) > gc_content_window_size) {
				for (const auto pos : iter::range<uint32_t>(0, len(seq) - gc_content_window_size + 1, gc_content_window_step)) {
					charCounter window_count(seq.seq_.substr(pos, gc_content_window_size));
					window_count.calcGcContent();
					if (window_count.gcContent_ > max_gc_content_window || window_count.gcContent_ < min_gc_content_window) {
						pass = false;
						break;
					}
				}
			}
		}
		return pass;
	};


	while (reader.readNextRead(seq)) {
		seqInfo out_seq(seq.name_, sampler->gen(seq.seq_));

		{
			uint32_t seq_attempts = 0;
			while (seq_attempts < max_seq_attempts && !check_seq(out_seq)) {
				out_seq = seqInfo(seq.name_, sampler->gen(seq.seq_));
				{
					uint32_t attempts = 0;
					while (has_restriction_sites_both_directions_checked(out_seq.seq_) ) {
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
				}
				++seq_attempts;
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

} //namespace njhseq


