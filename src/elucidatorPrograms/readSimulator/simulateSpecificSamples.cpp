//
// Created by Nicholas H	athaway on 3/31/25.
//

#include "readSimulatorRunner.hpp"
#include "elucidator/simulation.h"
#include <SeekDeep/objects/TarAmpSetupUtils/PrimersAndMids.hpp>
#include <njhseq/IO/SeqIO/SeqIO.hpp>
#include <njhseq/readVectorManipulation/readVectorHelpers/readVecTrimmer.hpp>

namespace njhseq {



int readSimulatorRunner::simulateSpecificSamples(const njh::progutils::CmdArgs & inputCommands) {
	auto check_for_level_match_up = [](const std::set<std::string> &levels1,
																										const bfs::path& fnp1,
																										const std::set<std::string> &levels2,
																										const bfs::path& fnp2,
																										const std::string& pretty_func_name) {
		VecStr only_in_1;
		VecStr only_in_2;
		VecStr in_both;
		njh::decompose_sets(levels1.begin(), levels1.end(),
												levels2.begin(), levels2.end(),
												std::back_inserter(only_in_1),
												std::back_inserter(only_in_2),
												std::back_inserter(in_both));
		if (!only_in_1.empty() || !only_in_2.empty()) {
			std::string message;
			if (!only_in_1.empty()) {
				message += njh::pasteAsStr("found the following levles only in ", fnp1, " and not in ", fnp2, ": ",
																	 njh::conToStr(only_in_1, ","), ".");
			}
			if (!only_in_2.empty()) {
				message += njh::pasteAsStr("found the following targets only in ", fnp2, " and not in ", fnp1, ": ",
																	 njh::conToStr(only_in_2, ","), ".");
			}
			std::stringstream ss;
			ss << pretty_func_name << ", error " << message << "\n";
			throw std::runtime_error{ss.str()};
		}
	};
	readSimulatorSetUp setUp(inputCommands);
	uint32_t numThreads = 2;
	uint32_t pcrNumThreads = 2;
	bool singleEnd = false;

	uint32_t defaultPcrRounds = 30;
	uint32_t initialPcrRounds = 10;
	std::map<uint32_t, uint32_t> initialPcrRoundsMap;
	bfs::path initialPcrRoundsMapFnp;
	long double errorRate = 3.5e-06;
	double pcrEfficiency = PCRSimulator::defaultPcrEfficiency;
	double templateSwitchRate = PCRSimulator::defaultTemplateSwitchRate;
	std::string templateSwitchScalingStr = PCRSimulator::templateSwitchScalingToStr(PCRSimulator::defaultTemplateSwitchScaling);
	double switchPlateau = PCRSimulator::defaultSwitchPlateau;
	double targetProductMax = PCRSimulator::defaultTargetProductMax;
	double switchExponent = PCRSimulator::defaultSwitchExponent;
	bool keepPCRSeqs = false;
	uint32_t chimeraBasesIn = 5;
	uint32_t templateCap = 500000000;
	bool noChimeras = false;
	double finalReadAmountSDFrac = 0.1;
	LibrarySetup::SimLibrarySetupPars sim_lib_pars;

	bfs::path population_haps_fnp;
	VecStr population_haps_required_cols{"target","seq","hap_id"};
	bfs::path sample_set_up_fnp;
	VecStr sample_set_up_required_cols{"sample_name","starting_template","sample_total_depth","target","hap_id","within_sample_fullhap_id","relative_abundance"};
	// optional, only used when idFile has MIDs, if absent MIDs are assigned automatically and samples are pooled into Subset-N indexes
	VecStr sample_set_up_mid_cols{"index", "mid"};
	// optional, if absent --pcrRounds is used for every sample
	std::string sample_set_up_pcr_rounds_col = "pcr_rounds";
	bfs::path target_read_count_factor_fnp;
	VecStr target_read_count_factor_required_cols{"target", "factor", "factor_stddev"};
	bfs::path idFile;
	bfs::path illuminaProfileDir = njh::pasteAsStr(elucidator_INSTALLDIR, "/etc/illumina_profiles/nextseq_150");

	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(sim_lib_pars.noAddPrimers_, "--noAddPrimers", "don't add primers");
	setUp.setOption(sim_lib_pars.addBluntEndingArtifact_, "--addBluntEndingArtifact", "add Blunt Ending Artifact");
	setUp.setOption(sim_lib_pars.addReverseComplement_, "--addReverseComplement", "add Reverse Complement");
	setUp.setOption(sim_lib_pars.bluntEndingArtifactChance_, "--bluntEndingArtifactChance", "Blunt Ending Artifact Chance");
	setUp.setOption(sim_lib_pars.barcodeRandomPrecedingBases_, "--barcodeRandomPrecedingBases", "Barcode Random Preceding Bases, up to this many random bases are added before the MIDs");
	setUp.setOption(sim_lib_pars.primerRandomPrecedingBases_, "--primerRandomPrecedingBases", "Primer Random Preceding Bases, up to this many random bases are added before the primers");

	setUp.setOption(illuminaProfileDir, "--illuminaProfileDir", "Illumina Profile Dir");
	setUp.setOption(idFile, "--idFile", "Primer MID Fnp, if it contains MIDs, samples are MID-tagged and pooled into index fastqs", true);
	setUp.setOption(population_haps_fnp, "--population_haps_fnp", njh::pasteAsStr("population haplotype sequences file, required columns: ", njh::conToStr(population_haps_required_cols, ",")), true);
	setUp.setOption(sample_set_up_fnp, "--sample_set_up_fnp", njh::pasteAsStr("sample setup file, required columns: ", njh::conToStr(sample_set_up_required_cols, ","), "; optional columns: ", sample_set_up_pcr_rounds_col, " (defaults to --pcrRounds), ", njh::conToStr(sample_set_up_mid_cols, ","), " (when idFile has MIDs)"), true);
	setUp.setOption(target_read_count_factor_fnp, "--target_read_count_factor_fnp", njh::pasteAsStr("target read count factor, used to adjust read depths based on known biases, required columns: ", njh::conToStr(target_read_count_factor_required_cols, ",")), true);


	setUp.setOption(finalReadAmountSDFrac, "--finalReadAmountSDFrac", "final Read Amount SD Frac, the standard deviation of a sample's total read depth as a fraction of sample_total_depth", njh::progutils::ProgramSetUp::CheckCase::GREATERZERO);
	setUp.setOption(noChimeras, "--noChimeras", "Don't simulate chimeras");
	uint64_t chimeraSummaryMinReads = 10;
	double chimeraSummaryMinFrac = 0.005;
	setUp.setOption(chimeraSummaryMinReads, "--chimeraSummaryMinReads", "Reporting only, simulation_summary.tsv also gives parent/chimera ratios for chimeras with more than this many reads in their sample x target (and at least --chimeraSummaryMinFrac)");
	setUp.setOption(chimeraSummaryMinFrac, "--chimeraSummaryMinFrac", "Reporting only, see --chimeraSummaryMinReads");
	setUp.setOption(templateCap, "--templateCap", "Template Cap");
	setUp.setOption(chimeraBasesIn, "--chimeraBasesIn", "The number of bases needed for a template to lay down");
	setUp.setOption(keepPCRSeqs, "--keepPCRSeqs", "Keep PCR Seqs");
	setUp.setOption(errorRate, "--errorRate", "Polymerase Error Rate");
	setUp.setOption(pcrEfficiency, "--pcrEfficiency", "PCR amplification efficiency, between 0-1: chance a template is copied in a round (yield, drift of low-template haplotypes, PCR-error copies); does not set the chimera rate, see --templateSwitchRate");
	setUp.setOption(templateSwitchRate, "--templateSwitchRate", "Template switch rate: with power or product scaling the rate when the target's product reaches --switchPlateau, with constant scaling the chance every round that a newly copied strand stops early and is finished on a different template; set 0 for no chimeras");
	setUp.setOption(templateSwitchScalingStr, "--templateSwitchScaling", "How the template switch rate changes over the PCR, product being all strands made so far for the target: constant (--templateSwitchRate every round), product (--templateSwitchRate * min(1, product / --switchPlateau), so switching happens mostly in late rounds) or power (min(1, --templateSwitchRate * (product / --switchPlateau)^--switchExponent), not capped, so the chimera share keeps rising with starting template)");
	setUp.setOption(switchPlateau, "--switchPlateau", "With --templateSwitchScaling product or power, the amount of product (strands of the target) at which the switch rate reaches --templateSwitchRate");
	setUp.setOption(switchExponent, "--switchExponent", "With --templateSwitchScaling power, the exponent of product / --switchPlateau in the switch rate; the chimera share rises with starting template about as template^switchExponent");
	setUp.setOption(targetProductMax, "--targetProductMax", "Shared per-target plateau for the rounds after --initialPcrRounds: each template is copied with probability --pcrEfficiency * max(0, 1 - product / --targetProductMax), product being all strands of the target; 0 uses --templateCap (a per-haplotype cap) instead");
	const auto templateSwitchScaling = PCRSimulator::templateSwitchScalingFromStr(templateSwitchScalingStr);
	setUp.setOption(defaultPcrRounds, "--pcrRounds", "Number of PCR rounds, used when the sample setup file has no pcr_rounds column");
	setUp.setOption(initialPcrRounds, "--initialPcrRounds", "Number of Initial rounds of PCR before sampling");
	setUp.setOption(initialPcrRoundsMapFnp, "--initialPcrRoundsTable", "Number of Initial rounds of PCR before sampling per starting template amount, columns 1)template, 2) rounds");
	setUp.setOption(numThreads, "--numThreads", "Number of Threads to Use");
	setUp.setOption(pcrNumThreads, "--pcrNumThreads", "Number of Threads to Use for PCR sim");
	setUp.setOption(singleEnd, "--singleEnd", "Single End, the whole read is simulated with the R1 Illumina profile, extended past its end by its last position");
	setUp.setOption(sim_lib_pars.pairedEndLength_, "--pairedEndLength", "Paired End Length", !singleEnd);
	bool extendIlluminaProfile = false;
	setUp.setOption(extendIlluminaProfile, "--extendIlluminaProfile", "Allow a paired end length longer than the Illumina profile, positions past its end use its last position");
	setUp.processDirectoryOutputName("simulateSpecificSamples_TODAY", true);
	setUp.finishSetUp(std::cout);

	setUp.startARunLog(setUp.pars_.directoryName_);
	setUp.writeParametersFile(setUp.pars_.directoryName_ + "parameters.tab.txt", false, true);

	initialPcrRoundsMap[1] = initialPcrRounds;
	if(bfs::exists(initialPcrRoundsMapFnp)){
		initialPcrRoundsMap.clear();
		table initialPcrRoundsMapTab(initialPcrRoundsMapFnp, "\t", true);
		initialPcrRoundsMapTab.checkForColumnsThrow(VecStr{"template", "rounds"}, __PRETTY_FUNCTION__);
		for(const auto & row : initialPcrRoundsMapTab){
			initialPcrRoundsMap[njh::StrToNumConverter::stoToNum<uint32_t>(row[initialPcrRoundsMapTab.getColPos("template")])] =njh::StrToNumConverter::stoToNum<uint32_t>(row[initialPcrRoundsMapTab.getColPos("rounds")]);
		}
	}
	uint64_t intErrorRate = errorRate * std::numeric_limits<uint64_t>::max();

	// get primer names
	PrimersAndMids primers(idFile);
	primers.initPrimerDeterminator();
	auto target_names = njh::getSetOfMapKeys(primers.targets_);
	// get read biases
	std::unordered_map<std::string, std::normal_distribution<double>> read_depth_biases;
	if (target_read_count_factor_fnp.empty()) {
		for (const auto & p : primers.targets_) {
			auto factor = 1.0/static_cast<double>(primers.targets_.size());
			auto factor_stddev = factor * .1;
			read_depth_biases[p.first] = std::normal_distribution<double>(factor, factor_stddev);
		}
	} else {
		table target_read_count_factor(target_read_count_factor_fnp, "\t", true);
		target_read_count_factor.checkForColumnsThrow(target_read_count_factor_required_cols, __PRETTY_FUNCTION__);

		for (const auto & row : target_read_count_factor) {
			auto factor = njh::StrToNumConverter::stoToNum<double>(row[target_read_count_factor.getColPos("factor")]);
			auto factor_stddev = njh::StrToNumConverter::stoToNum<double>(row[target_read_count_factor.getColPos("factor_stddev")]);
			if (factor <= 0 || factor_stddev < 0) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "for target " << row[target_read_count_factor.getColPos("target")]
				<< " factor must be greater than 0 and factor_stddev can't be negative, found factor: " << factor << ", factor_stddev: " << factor_stddev << "\n";
				throw std::runtime_error{ss.str()};
			}
			read_depth_biases[row[target_read_count_factor.getColPos("target")]] = std::normal_distribution<double>(factor, factor_stddev);;
		}
		auto target_names_in_biases = njh::getSetOfMapKeys(read_depth_biases);
		check_for_level_match_up(target_names_in_biases, target_read_count_factor_fnp,
		                         target_names, idFile,
		                         __PRETTY_FUNCTION__);
	}
	//std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
	// population haplotypes
	table population_haps(population_haps_fnp, "\t", true);
	population_haps.checkForColumnsThrow(population_haps_required_cols, __PRETTY_FUNCTION__);
	auto population_haps_target_names = njh::vecToSet(population_haps.getColumnLevels(population_haps.getColPos("target")));
	check_for_level_match_up(population_haps_target_names, population_haps_fnp,
												 target_names, idFile,
												 __PRETTY_FUNCTION__);
	//std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
	std::unordered_map<std::string, std::unordered_map<uint32_t, seqInfo>> population_haps_by_target;
	for (const auto & row : population_haps) {
		auto target_name = row[population_haps.getColPos("target")];
		auto hap_id = njh::StrToNumConverter::stoToNum<uint32_t>(row[population_haps.getColPos("hap_id")]);
		if (njh::in(hap_id, population_haps_by_target[target_name])) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "already have hap_id " << hap_id << " for target " << target_name << "\n";
			throw std::runtime_error{ss.str()};
		}
		population_haps_by_target[target_name][hap_id] = seqInfo(row[population_haps.getColPos("hap_id")], row[population_haps.getColPos("seq")]);
	}
	//std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
	// sample set up
	table sample_setup(sample_set_up_fnp, "\t", true);
	sample_setup.checkForColumnsThrow(sample_set_up_required_cols, __PRETTY_FUNCTION__);
	auto sample_target_names = njh::vecToSet(sample_setup.getColumnLevels(sample_setup.getColPos("target")));
	check_for_level_match_up(sample_target_names, sample_set_up_fnp,
												 target_names, idFile,
												 __PRETTY_FUNCTION__);
	//std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
	VecStr warnings;
	struct SimSampleInfo {
		std::string sample_name_;
		uint32_t starting_template_ = std::numeric_limits<uint32_t>::max();
		uint32_t sample_total_depth_ = std::numeric_limits<uint32_t>::max();
		uint32_t pcr_rounds_ = std::numeric_limits<uint32_t>::max();
		std::string index_name_; // the fastq file this sample is written to, the sample name itself when no MIDs
		std::string mid_name_;   // empty when no MIDs

		struct FullHap {
			FullHap(uint32_t within_sample_fullhap_id, double relative_abundance): within_sample_fullhap_id_(
				                                                                       within_sample_fullhap_id),
			                                                                       relative_abundance_(relative_abundance) {
			}
			FullHap() = default;
			uint32_t within_sample_fullhap_id_ = std::numeric_limits<uint32_t>::max();;
			double relative_abundance_ = std::numeric_limits<double>::max();
			std::unordered_map<std::string, std::vector<uint32_t>> hap_ids_;
		};
		std::map<uint32_t, FullHap> full_haps_;
		// draw from a normal distribution, re-drawing until positive so the tail below zero can't produce negative/zero factors or depths
		static double draw_positive(std::normal_distribution<double> & dist, njh::randomGenerator& rGen) {
			if (dist.mean() <= 0) {
				return 0;
			}
			if (dist.stddev() <= 0) {
				return dist.mean();
			}
			double ret = 0;
			while (ret <= 0) {
				ret = dist(rGen.mtGen_);
			}
			return ret;
		}
		std::map<std::string, uint64_t> gen_read_depths_per_target(std::unordered_map<std::string, std::normal_distribution<double>> & read_depth_biases, njh::randomGenerator& rGen, double total_depth_sd_frac) const {
			std::map<std::string, uint64_t> ret;
			std::map<std::string, double> relative_depths;
			std::map<std::string, uint32_t> hap_counts_per_targets;
			for (const auto & full_hap : full_haps_) {
				for (const auto & target : full_hap.second.hap_ids_) {
					hap_counts_per_targets[target.first] += target.second.size();
				}
			}
			double total_bias_factor = 0;
			for (const auto & hap_counts_per_target : hap_counts_per_targets) {
				double bias_factor = 0;
				for (uint32_t hap = 0; hap < hap_counts_per_target.second; ++hap) {
					bias_factor += draw_positive(read_depth_biases.at(hap_counts_per_target.first), rGen);
				}
				relative_depths[hap_counts_per_target.first] = bias_factor;
				total_bias_factor += bias_factor;
			}
			//renormalize the bias factor
			for (auto & relative_depth : relative_depths) {
				relative_depth.second /= total_bias_factor;
			}
			std::normal_distribution<double> total_depth_dist(sample_total_depth_, sample_total_depth_ * total_depth_sd_frac);
			for (const auto & target_bias : relative_depths) {
				auto raw_read_depth = std::round(target_bias.second * draw_positive(total_depth_dist, rGen));
				uint64_t read_depth = 0;
				if (raw_read_depth >= 1) {
					read_depth = static_cast<uint64_t>(raw_read_depth);
				}
				ret[target_bias.first] = read_depth;
			}
			return ret;
		}

		[[nodiscard]] std::map<std::string, std::unordered_map<uint32_t, double>> gen_relative_depths_per_hap_id() const {
			std::map<std::string, std::unordered_map<uint32_t, double>> ret;
			for (const auto & full_hap : full_haps_) {
				for (const auto & haps : full_hap.second.hap_ids_) {
					for (const auto & hap : haps.second) {
						ret[haps.first][hap] += full_hap.second.relative_abundance_;
					}
				}
			}
			//re-normalize so abundances add up to 1
			for (auto & target : ret) {
				double total = 0;
				for (const auto & hap : target.second) {
					total += hap.second;
				}
				for (auto & hap : target.second) {
					hap.second = hap.second/total;
				}
			}
			return ret;
		}

		[[nodiscard]] std::map<std::string, std::unordered_map<uint32_t, uint32_t>> gen_number_per_hap_id() const {
			std::map<std::string, std::unordered_map<uint32_t, uint32_t>> ret;
			for (const auto & full_hap : full_haps_) {
				for (const auto & haps : full_hap.second.hap_ids_) {
					for (const auto & hap : haps.second) {
						++ret[haps.first][hap];
					}
				}
			}
			return ret;
		}

		[[nodiscard]] std::map<std::string, std::unordered_map<uint32_t, uint64_t>> gen_starting_templates_per_hap_id(njh::randomGenerator & rGen, double starting_template_frac_dev = 0.1) const {
			std::unordered_map<std::string, std::unordered_map<uint32_t, double>> starting_templates;
			for (const auto & full_hap : full_haps_) {
				for (const auto & haps : full_hap.second.hap_ids_) {
					for (const auto & hap : haps.second) {
						std::normal_distribution<double> template_ndist(full_hap.second.relative_abundance_ * starting_template_, full_hap.second.relative_abundance_ * starting_template_ * starting_template_frac_dev);
						starting_templates[haps.first][hap] += template_ndist(rGen.mtGen_);
					}
				}
			}
			std::map<std::string, std::unordered_map<uint32_t, uint64_t>> ret;
			for (auto & target : starting_templates) {
				for (auto & hap : target.second) {
					double hap_template_amount = std::round(hap.second);
					if (hap_template_amount >= 1) {
						ret[target.first][hap.first] = static_cast<uint32_t>(hap_template_amount);
					} else {
						ret[target.first][hap.first] = 0;
					}
				}
			}
			return ret;
		}
	};
	//std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
	// populate samples
	std::map<std::string, SimSampleInfo> sim_samples;
	const bool useMids = primers.containsMids();
	const bool hasMidCols = sample_setup.containsColumns(sample_set_up_mid_cols);
	if (!useMids && (sample_setup.containsColumn("mid") || sample_setup.containsColumn("index"))) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << sample_set_up_fnp << " has mid/index columns but " << idFile << " doesn't contain any MIDs" << "\n";
		throw std::runtime_error{ss.str()};
	}
	if (useMids && !hasMidCols && (sample_setup.containsColumn("mid") || sample_setup.containsColumn("index"))) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << sample_set_up_fnp << " if supplying MIDs need both columns: " << njh::conToStr(sample_set_up_mid_cols, ",") << "\n";
		throw std::runtime_error{ss.str()};
	}
	VecStr sample_only_cols{"sample_name","starting_template","sample_total_depth"};
	const bool hasPcrRoundsCol = sample_setup.containsColumn(sample_set_up_pcr_rounds_col);
	if (hasPcrRoundsCol) {
		sample_only_cols.emplace_back(sample_set_up_pcr_rounds_col);
	}
	if (useMids && hasMidCols) {
		addOtherVec(sample_only_cols, sample_set_up_mid_cols);
	}
	auto sample_only_info = sample_setup.getColumns(sample_only_cols).getUniqueRows();
	for (const auto & row : sample_only_info) {
		auto sample_name = row[sample_only_info.getColPos("sample_name")];
		if (njh::in(sample_name, sim_samples)) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "already have sample info for :" << sample_name << "\n";
			throw std::runtime_error{ss.str()};
		}
		SimSampleInfo sim_sample;
		sim_sample.sample_name_ = sample_name;
		// std::cout << "row[sample_only_info.getColPos(starting_template)]" << row[sample_only_info.getColPos("starting_template")] << std::endl;
		// std::cout << "row[sample_only_info.getColPos(sample_total_depth)]" << row[sample_only_info.getColPos("sample_total_depth")] << std::endl;
		// std::cout << "row[sample_only_info.getColPos(pcr_rounds)]" << row[sample_only_info.getColPos("pcr_rounds")] << std::endl;

		sim_sample.starting_template_ = njh::StrToNumConverter::stoToNum<uint32_t>(row[sample_only_info.getColPos("starting_template")]);
		sim_sample.sample_total_depth_ = njh::StrToNumConverter::stoToNum<uint32_t>(row[sample_only_info.getColPos("sample_total_depth")]);
		sim_sample.pcr_rounds_ = hasPcrRoundsCol ? njh::StrToNumConverter::stoToNum<uint32_t>(row[sample_only_info.getColPos(sample_set_up_pcr_rounds_col)]) : defaultPcrRounds;
		sim_sample.index_name_ = sample_name;
		if (useMids && hasMidCols) {
			sim_sample.index_name_ = row[sample_only_info.getColPos("index")];
			sim_sample.mid_name_ = row[sample_only_info.getColPos("mid")];
		}
		sim_samples.emplace(sample_name, sim_sample);
	}
	// set up MIDs, either checking the given ones or assigning them in sorted order, starting a new index when MIDs run out
	if (useMids) {
		auto midNames = njh::getVecOfMapKeys(primers.mids_);
		njh::sort(midNames);
		if (hasMidCols) {
			VecStr mid_warnings;
			std::map<std::string, std::map<std::string, std::string>> samples_by_index_by_mid;
			for (const auto & sim_sample : sim_samples) {
				if (njh::notIn(sim_sample.second.mid_name_, primers.mids_)) {
					mid_warnings.emplace_back(njh::pasteAsStr("sample: ", sim_sample.first, " has mid: ", sim_sample.second.mid_name_, " which isn't in ", idFile));
				}
				auto & samples_by_mid = samples_by_index_by_mid[sim_sample.second.index_name_];
				if (njh::in(sim_sample.second.mid_name_, samples_by_mid)) {
					mid_warnings.emplace_back(njh::pasteAsStr("samples: ", samples_by_mid[sim_sample.second.mid_name_], " and ", sim_sample.first, " share mid: ", sim_sample.second.mid_name_, " in index: ", sim_sample.second.index_name_));
				}
				samples_by_mid[sim_sample.second.mid_name_] = sim_sample.first;
			}
			if (!mid_warnings.empty()) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << njh::conToStr(mid_warnings, "\n") << "\n";
				throw std::runtime_error{ss.str()};
			}
		} else {
			uint32_t indexTotal = (sim_samples.size() + midNames.size() - 1) / midNames.size();
			uint32_t sampleCount = 0;
			for (auto & sim_sample : sim_samples) {
				sim_sample.second.index_name_ = njh::pasteAsStr("Subset-", njh::leftPadNumStr<uint32_t>(sampleCount / midNames.size(), indexTotal));
				sim_sample.second.mid_name_ = midNames[sampleCount % midNames.size()];
				++sampleCount;
			}
		}
	}
	//std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
	// populate haplotypes
	auto within_sample_hap_only_info = sample_setup.getColumns(VecStr{"sample_name","within_sample_fullhap_id","relative_abundance"}).getUniqueRows();
	for (const auto & row : within_sample_hap_only_info) {
		auto sample_name = row[within_sample_hap_only_info.getColPos("sample_name")];
		auto within_sample_fullhap_id = njh::StrToNumConverter::stoToNum<uint32_t>(row[within_sample_hap_only_info.getColPos("within_sample_fullhap_id")]);
		if (njh::in(within_sample_fullhap_id, sim_samples.at(sample_name).full_haps_)) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__
			<< ", error " << "already have within_sample_fullhap_id: "
			<< within_sample_fullhap_id << " with relative_abundance of: "
			<< sim_samples.at(sample_name).full_haps_[within_sample_fullhap_id].relative_abundance_ << "\n";
			throw std::runtime_error{ss.str()};
		}
		auto relative_abundance = njh::StrToNumConverter::stoToNum<double>(row[within_sample_hap_only_info.getColPos("relative_abundance")]);
		sim_samples.at(sample_name).full_haps_.emplace(within_sample_fullhap_id, SimSampleInfo::FullHap(within_sample_fullhap_id, relative_abundance));
	}
	//std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
	// re -normalize relative_abundances
	for (auto & sim_sample : sim_samples) {
		double total = 0;
		for (const auto & hap : sim_sample.second.full_haps_) {
			total += hap.second.relative_abundance_;
		}
		for (auto & hap : sim_sample.second.full_haps_) {
			hap.second.relative_abundance_ /= total;
		}
	}
	//std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
	VecStr missing_hap_warnings;
	for (const auto & row : sample_setup) {
		auto sample_name = row[sample_setup.getColPos("sample_name")];
		auto target_name = row[sample_setup.getColPos("target")];
		auto hap_id = njh::StrToNumConverter::stoToNum<uint32_t>(row[sample_setup.getColPos("hap_id")]);
		auto within_sample_fullhap_id = njh::StrToNumConverter::stoToNum<uint32_t>(row[sample_setup.getColPos("within_sample_fullhap_id")]);
		if (njh::notIn(hap_id, population_haps_by_target.at(target_name))) {
			missing_hap_warnings.emplace_back(njh::pasteAsStr("missing hap_id: ", hap_id, " for target: ", target_name));
		}
		sim_samples.at(sample_name).full_haps_[within_sample_fullhap_id].hap_ids_[target_name].emplace_back(hap_id);
	}
	if (!missing_hap_warnings.empty()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << njh::conToStr(missing_hap_warnings, "\n") << "\n";
		throw std::runtime_error{ss.str()};
	}
	bfs::path fastqDirectory = njh::files::make_path(setUp.pars_.directoryName_, "fastq");
	bfs::path chimeraDirectory = njh::files::make_path(setUp.pars_.directoryName_, "chimeras");
	njh::files::makeDir(njh::files::MkdirPar{fastqDirectory});
	if(!noChimeras){
		njh::files::makeDir(njh::files::MkdirPar{chimeraDirectory});
	}
	// group samples by the fastq (index) they are written to, without MIDs every sample is its own index
	std::map<std::string, VecStr> samples_by_index;
	for (const auto & sim_sample : sim_samples) {
		samples_by_index[sim_sample.second.index_name_].emplace_back(sim_sample.first);
	}
	auto index_names = njh::getVecOfMapKeys(samples_by_index);
	std::unordered_map<std::string, uint32_t> index_sample_numbers;
	for (const auto & index_name : index_names) {
		index_sample_numbers[index_name] = index_sample_numbers.size() + 1;
	}
	if (useMids) {
		table sample_names_table(VecStr{"#file", "sample", "run1"});
		for (const auto & sim_sample : sim_samples) {
			sample_names_table.addRow(sim_sample.second.index_name_, sim_sample.first, sim_sample.second.mid_name_);
		}
		sample_names_table.sortTable("#file", "run1", false);
		OutputStream sample_names_out(njh::files::make_path(setUp.pars_.directoryName_, "sampleNames.tab.txt"));
		sample_names_table.outPutContents(sample_names_out, "\t");
	}
	std::string extension = ".fastq.gz";
	table output_numbers(VecStr{
		"index",
		"mid",
		"sample_name",
		"COI",
		"expected_total_sample_reads",
		"target",
		"unique_haps_for_target",
		"total_reads_for_target",
		"hap_id",
		"expected_starting_template", "sample_target_starting_template",
		"hap_expected_relative_abundance", "hap_sampled_starting_template",
		"hap_pcred_sample_counts",
		"total_pcr_count",
		"non_chimera_total_pcr_count"
	});
	std::mutex output_numbers_mut;

	// error and chimera reporting
	// per source (true haplotype or chimera) per sample x target read counts, PCR errors, sequencing errors and chimera parents
	// chimera_first_pcr_round is the round the chimera's sequence first formed; the same sequence formed again in later rounds is counted under that first name
	// pcr errors are substitutions vs the source sequence, sequencing errors are substitutions vs the fragment over the bases each mate covers (seq_bases)
	table source_read_counts(VecStr{
		"index", "mid", "sample_name", "target",
		"source", "is_chimera", "chimera_first_pcr_round",
		"reads", "frac_of_target_reads",
		"source_length", "reads_with_pcr_errors", "pcr_error_bases",
		"reads_with_seq_errors", "seq_error_bases", "seq_bases",
		"parent_5prime", "parent_3prime", "parent_5prime_reads", "parent_3prime_reads",
		"weaker_parent_to_chimera_ratio"
	});
	struct SampleTargetSummary {
		std::string index_;
		std::string mid_;
		std::string sample_name_;
		std::string target_;
		uint32_t haps_amplified_{0};
		uint64_t starting_genomes_{0};
		uint64_t finished_copies_{0};
		uint64_t partial_copies_{0};
		uint32_t switch_plateau_round_{0};
		double final_product_{0};
		uint64_t total_reads_{0};
		uint64_t chimeric_reads_{0};
		uint32_t chimera_seqs_{0};
		uint32_t chimera_seqs_two_parent_{0};
		uint64_t pcr_reads_with_errors_{0};
		uint64_t pcr_error_bases_{0};
		uint64_t pcr_bases_{0};
		uint64_t seq_reads_with_errors_{0};
		uint64_t seq_error_bases_{0};
		uint64_t seq_bases_{0};
	};
	struct ChimeraRecord {
		std::string target_;
		uint64_t reads_{0};
		double frac_{0};
		bool two_parent_{false};
		double weaker_parent_ratio_{std::numeric_limits<double>::quiet_NaN()};
	};
	std::vector<SampleTargetSummary> sample_target_summaries;
	std::vector<ChimeraRecord> chimera_records;
	njh::concurrent::LockableQueue<std::string> index_queue(index_names);

	// wall seconds summed over samples and sample x targets, for profiling (phase_timings.tsv): the phases of
	// PCRSimulator::simLibFast plus simLibFast_total, target_total (simulation and reporting per sample x target) and sample_total
	std::map<std::string, double> phase_seconds;
	std::mutex phase_seconds_mut;
	struct PhaseTimer {
		std::map<std::string, double> & seconds_;
		std::mutex & mut_;
		std::string phase_;
		njh::stopWatch watch_;
		PhaseTimer(std::map<std::string, double> & seconds, std::mutex & mut, std::string phase) : seconds_(seconds), mut_(mut), phase_(std::move(phase)) {}
		~PhaseTimer() {
			std::lock_guard<std::mutex> lock(mut_);
			seconds_[phase_] += watch_.totalTime();
		}
	};
	std::function<void()> sim_sample_func = [&index_queue, &sim_samples, &phase_seconds, &phase_seconds_mut,
		&samples_by_index, &index_sample_numbers, &useMids,
		&illuminaProfileDir, &sim_lib_pars,
		&keepPCRSeqs, &singleEnd, &extendIlluminaProfile,
		&fastqDirectory, &extension,
		&read_depth_biases,
		&primers, &setUp,
		&pcrNumThreads,
		&intErrorRate, &pcrEfficiency, &templateSwitchRate, &templateSwitchScaling, &switchPlateau, &switchExponent, &targetProductMax, &templateCap,&noChimeras, &chimeraBasesIn,
		&initialPcrRounds, &initialPcrRoundsMap, &finalReadAmountSDFrac,
		&population_haps_by_target,
		&output_numbers, &output_numbers_mut,
		&source_read_counts, &sample_target_summaries, &chimera_records]() {
		// number of positions that differ, plus any length difference
		auto count_mismatches = [](const std::string & seq, const std::string & ref) -> uint32_t {
			uint32_t ret = 0;
			const auto shared_len = std::min(seq.size(), ref.size());
			for (uint32_t pos = 0; pos < shared_len; ++pos) {
				if (seq[pos] != ref[pos]) {
					++ret;
				}
			}
			ret += std::max(seq.size(), ref.size()) - shared_len;
			return ret;
		};
		// sequencing errors in a simulated read, only over the bases that come from the fragment (not the overhang past its end)
		auto count_seq_errors = [](const std::string & fragment, const std::string & read) -> uint32_t {
			uint32_t ret = 0;
			const auto shared_len = std::min(fragment.size(), read.size());
			for (uint32_t pos = 0; pos < shared_len; ++pos) {
				if (fragment[pos] != read[pos]) {
					++ret;
				}
			}
			return ret;
		};
		auto common_prefix_len = [](const std::string & a, const std::string & b) -> uint32_t {
			uint32_t ret = 0;
			while (ret < a.size() && ret < b.size() && a[ret] == b[ret]) {
				++ret;
			}
			return ret;
		};
		auto common_suffix_len = [](const std::string & a, const std::string & b) -> uint32_t {
			uint32_t ret = 0;
			while (ret < a.size() && ret < b.size() && a[a.size() - 1 - ret] == b[b.size() - 1 - ret]) {
				++ret;
			}
			return ret;
		};
		std::string indexKey;
		//std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
		njh::randObjectGen<char,uint32_t> baseRGen({'A', 'C', 'G', 'T'}, {1,1,1,1});
		RoughIlluminaSimulator simulator(illuminaProfileDir);
		// single end reads are the whole amplicon so can be longer than the profile
		simulator.extendPastProfile_ = singleEnd || extendIlluminaProfile;
		njh::randomGenerator rGen;
		// per thread copy, drawing from a std::normal_distribution modifies it so it can't be shared between threads
		auto thread_read_depth_biases = read_depth_biases;
		auto add_random_bases = [&rGen, &baseRGen](uint32_t maxBases, bool randomOneLength) -> std::string {
			if (0 == maxBases) {
				return "";
			}
			uint32_t minBases = randomOneLength ? maxBases : 0;
			uint32_t numBases = rGen.unifRand<uint32_t>(minBases, maxBases + 1);
			return njh::pasteAsStr(baseRGen.genObjs(numBases));
		};
		while (index_queue.getVal(indexKey)) {
			PhaseTimer sampleTimer(phase_seconds, phase_seconds_mut, "sample_total");
			uint32_t sampleNumber = index_sample_numbers.at(indexKey);
			OutOptions targetOutOpts(njh::files::make_path(fastqDirectory, njh::pasteAsStr(indexKey, "_S", sampleNumber, extension)));
			OutOptions r1OutOpts(njh::files::make_path(fastqDirectory, njh::pasteAsStr(indexKey, "_S", sampleNumber, "_R1_001" + extension)));
			OutOptions r2OutOpts(njh::files::make_path(fastqDirectory, njh::pasteAsStr(indexKey, "_S", sampleNumber, "_R2_001" + extension)));
			std::shared_ptr<OutputStream> targetOut;
			std::shared_ptr<OutputStream> r1Out;
			std::shared_ptr<OutputStream> r2Out;
			if(singleEnd) {
				targetOut = std::make_shared<OutputStream>(targetOutOpts);
			} else {
				r1Out = std::make_shared<OutputStream>(r1OutOpts);
				r2Out = std::make_shared<OutputStream>(r2OutOpts);
			}

			/**@todo add back in the ability to have different Illumina barcodes in overhangs */
	//		std::string sampleAdapter1 = defaultAdapter1;
	//		std::string sampleAdapter2 = defaultAdapter2;
	//		for(const auto pos : iter::range(sampleAdapter1.size())){
	//			if('N' == sampleAdapter1[pos]){
	//				sampleAdapter1[pos] = baseRGen.genObj();
	//			}
	//		}
	//		for(const auto pos : iter::range(sampleAdapter2.size())){
	//			if('N' == sampleAdapter2[pos]){
	//				sampleAdapter2[pos] = baseRGen.genObj();
	//			}
	//		}

			for (const auto & sampleKey : samples_by_index.at(indexKey)) {
				const auto & sim_sample = sim_samples.at(sampleKey);
				auto read_depths_per_target = sim_sample.gen_read_depths_per_target(thread_read_depth_biases, rGen, finalReadAmountSDFrac);
				auto expected_relative_depths_per_hap_id = sim_sample.gen_relative_depths_per_hap_id();
				auto starting_templates_per_hap_id = sim_sample.gen_starting_templates_per_hap_id(rGen);
				// auto copy_number_per_hap_id = sim_sample.gen_number_per_hap_id();



				for (const auto & target : starting_templates_per_hap_id) {
					PhaseTimer targetTimer(phase_seconds, phase_seconds_mut, "target_total");
					std::vector<PCRSimulator::SeqGenomeCnt> seqGCounts;
					// every expected hap including ones that ended up with 0 starting templates, used for reporting only
					std::vector<PCRSimulator::SeqGenomeCnt> allSeqGCounts;
					uint32_t totalGenomes = 0;
					for (const auto & hap_id : target.second) {
						std::string seqName = njh::pasteAsStr(hap_id.first);
						std::string seq = population_haps_by_target.at(target.first).at(hap_id.first).seq_;
						if (!sim_lib_pars.noAddPrimers_) {
							seq.append(seqUtil::reverseComplement(primers.targets_.at(target.first).info_.reversePrimerRaw_, "DNA"));
							seq.insert(0, primers.targets_.at(target.first).info_.forwardPrimerRaw_);
						}
						/**@todo add back in randomizing this per seq, currently random bases are set once per haplotype template like simMultipleMixtureSimPCR*/
						seq.insert(0, add_random_bases(sim_lib_pars.primerRandomPrecedingBases_, false));
						seq.append(add_random_bases(sim_lib_pars.primerRandomPrecedingBases_, false));
						if (useMids) {
							const auto & mid = primers.mids_.at(sim_sample.mid_name_);
							if (nullptr != mid.forwardBar_) {
								seq.insert(0, mid.forwardBar_->bar_->motifOriginal_);
								seq.insert(0, add_random_bases(sim_lib_pars.barcodeRandomPrecedingBases_, false));
							}
							if (nullptr != mid.reverseBar_) {
								seq.append(seqUtil::reverseComplement(mid.reverseBar_->bar_->motifOriginal_, "DNA"));
								seq.append(add_random_bases(sim_lib_pars.barcodeRandomPrecedingBases_, false));
							}
						}
						seqInfo currentSeq(seqName, seq);

						currentSeq.frac_ = expected_relative_depths_per_hap_id.at(target.first).at(hap_id.first);
						currentSeq.cnt_ = expected_relative_depths_per_hap_id.at(target.first).at(hap_id.first);

						allSeqGCounts.emplace_back(currentSeq, hap_id.second);
						// haps with no starting templates can't be amplified, so leave them out of the PCR like simMultipleMixtureSimPCR does
						if (hap_id.second > 0) {
							seqGCounts.emplace_back(currentSeq, hap_id.second);
							totalGenomes += hap_id.second;
						}
					}

					PCRSimulator pcrSim(intErrorRate);
					pcrSim.verbose_ = setUp.pars_.verbose_;
					pcrSim.pcrEfficiency_ = pcrEfficiency;
					pcrSim.templateSwitchRate_ = templateSwitchRate;
					pcrSim.templateSwitchScaling_ = templateSwitchScaling;
					pcrSim.switchPlateau_ = switchPlateau;
					pcrSim.targetProductMax_ = targetProductMax;
					pcrSim.switchExponent_ = switchExponent;
					pcrSim.templateCap_ = templateCap;
					pcrSim.noChimeras_ = noChimeras;
					pcrSim.chimeraPad_ = chimeraBasesIn;


					auto pcrReadsFnp = njh::files::make_path(fastqDirectory, sampleKey  + "_" + target.first + "_pcr_reads.fasta");
					OutOptions pcrReadsOpts(pcrReadsFnp);

					uint32_t currentInitial = initialPcrRounds;
					for(const auto & initialRound : initialPcrRoundsMap){
						if(initialRound.first >= totalGenomes){
							currentInitial = initialRound.second;
							break;
						}
					}
					std::unordered_map<std::string, uint32_t> pcred_read_counts;
					// the sequence each read was amplified from, the input haplotypes plus any chimeras the PCR made
					std::unordered_map<std::string, std::string> source_seqs;
					for (const auto & seqGCount : seqGCounts) {
						source_seqs[seqGCount.seqBase_.name_] = seqGCount.seqBase_.seq_;
					}
					struct SourceStats {
						uint64_t reads_{0};
						uint64_t reads_with_pcr_errors_{0};
						uint64_t pcr_error_bases_{0};
						uint64_t reads_with_seq_errors_{0};
						uint64_t seq_error_bases_{0};
						uint64_t seq_bases_{0};
					};
					std::map<std::string, SourceStats> source_stats;
					// copies made while chimeras were modelled, for the copy-weighted switch rate in sample_target_summary.tsv
					uint64_t simFinishedCopies = 0;
					uint64_t simPartialCopies = 0;
					uint32_t simSwitchPlateauRound = 0;
					double simFinalProduct = 0;
					if (read_depths_per_target.at(target.first) > 0 && totalGenomes > 0) {
						njh::stopWatch simLibFastWatch;
						auto pcrSimAmounts = pcrSim.simLibFast(seqGCounts, pcrReadsOpts, read_depths_per_target.at(target.first), sim_sample.pcr_rounds_, currentInitial, pcrNumThreads);
						{
							std::lock_guard<std::mutex> lock(phase_seconds_mut);
							phase_seconds["simLibFast_total"] += simLibFastWatch.totalTime();
							for (const auto & phase : pcrSimAmounts.phaseSeconds_) {
								phase_seconds[phase.first] += phase.second;
							}
						}
						simFinishedCopies = pcrSimAmounts.finishedCopies_;
						simPartialCopies = pcrSimAmounts.partialCopies_;
						simSwitchPlateauRound = pcrSimAmounts.switchPlateauRound_;
						simFinalProduct = pcrSimAmounts.finalProduct_;
						for (const auto & chi : pcrSimAmounts.chimeraSeqs_) {
							source_seqs[chi.name_] = chi.seq_;
						}
						auto pcrSeqsInOpts = SeqIOOptions::genFastaIn(pcrReadsFnp);
						SeqInput pcrReader(pcrSeqsInOpts);
						pcrReader.openIn();
						seqInfo targetSeq;

						while(pcrReader.readNextRead(targetSeq)){
							MetaDataInName targetSeqMeta(targetSeq.name_);
							const auto source_name = targetSeqMeta.getMeta("hap");
							++pcred_read_counts[source_name];
							auto source_seq_search = source_seqs.find(source_name);
							if (source_seqs.end() == source_seq_search) {
								std::stringstream ss;
								ss << __PRETTY_FUNCTION__ << ", error " << "no source sequence for read " << targetSeq.name_ << " in sample " << sampleKey << " target " << target.first << "\n";
								throw std::runtime_error{ss.str()};
							}
							auto & current_source_stats = source_stats[source_name];
							++current_source_stats.reads_;
							const uint32_t pcr_errors = count_mismatches(targetSeq.seq_, source_seq_search->second);
							if (pcr_errors > 0) {
								++current_source_stats.reads_with_pcr_errors_;
								current_source_stats.pcr_error_bases_ += pcr_errors;
							}
							// hap ids and idNums repeat across targets and samples, so add where the read came from to keep names unique within the output fastq
							targetSeqMeta.addMeta("sample", sampleKey);
							targetSeqMeta.addMeta("target", target.first);
							if (useMids) {
								targetSeqMeta.addMeta("mid", sim_sample.mid_name_);
							}
							targetSeqMeta.addMeta("pcrErrors", pcr_errors);
							targetSeqMeta.resetMetaInName(targetSeq.name_);
							//blunt ending
							if(sim_lib_pars.addBluntEndingArtifact_){
								if('A' == targetSeq.seq_[0] && rGen() < sim_lib_pars.bluntEndingArtifactChance_){
									readVecTrimmer::trimOffForwardBases(targetSeq, 1);
								}
								if ('T' == targetSeq.seq_.back()
										&& rGen() < sim_lib_pars.bluntEndingArtifactChance_) {
									readVecTrimmer::trimOffEndBases(targetSeq, 1);
										}
							}
							//complement
							if(sim_lib_pars.addReverseComplement_){
								bool complement = rGen.unifRand(0,2) == 0;
								if(complement){
									targetSeq.reverseComplementRead(false, true);
								}
							}
							//length
							uint32_t seq_errors = 0;
							if(singleEnd){
								//target, no long read profile yet so errors and qualities come from the R1 profile
								auto simRead = simulator.simR1(targetSeq, len(targetSeq));
								seq_errors = count_seq_errors(targetSeq.seq_, simRead.seq_);
								current_source_stats.seq_bases_ += len(targetSeq);
								targetSeqMeta.addMeta("seqErrors", seq_errors);
								targetSeqMeta.resetMetaInName(simRead.name_);
								simRead.outPutFastq(*targetOut);
							}else{
								//r1
								auto r1SubSeq = len(targetSeq) > sim_lib_pars.pairedEndLength_? targetSeq.getSubRead(0, sim_lib_pars.pairedEndLength_) : targetSeq;
								auto r1Read = simulator.simR1(r1SubSeq, sim_lib_pars.pairedEndLength_);
								seq_errors += count_seq_errors(r1SubSeq.seq_, r1Read.seq_);
								current_source_stats.seq_bases_ += len(r1SubSeq);
								//r2
								targetSeq.reverseComplementRead(false, true);
								auto r2SubSeq = len(targetSeq) > sim_lib_pars.pairedEndLength_? targetSeq.getSubRead(0, sim_lib_pars.pairedEndLength_) : targetSeq;
								auto r2Read = simulator.simR2(r2SubSeq, sim_lib_pars.pairedEndLength_);
								seq_errors += count_seq_errors(r2SubSeq.seq_, r2Read.seq_);
								current_source_stats.seq_bases_ += len(r2SubSeq);
								// seqErrors is the total over both mates
								targetSeqMeta.addMeta("seqErrors", seq_errors);
								targetSeqMeta.resetMetaInName(r1Read.name_);
								targetSeqMeta.resetMetaInName(r2Read.name_);
								r1Read.outPutFastq(*r1Out);
								r2Read.outPutFastq(*r2Out);
							}
							if (seq_errors > 0) {
								++current_source_stats.reads_with_seq_errors_;
								current_source_stats.seq_error_bases_ += seq_errors;
							}
						}
						pcrReader.closeIn();
						if(!keepPCRSeqs) {
							bfs::remove(pcrReadsFnp);
						}
					}
					uint32_t chimera_pcr_count = 0;
					uint32_t total_pcr_count = 0;
					for (const auto & pcred_read_count : pcred_read_counts) {
						total_pcr_count += pcred_read_count.second;
						if (std::string::npos !=  pcred_read_count.first.find("Chi")){
							chimera_pcr_count+= pcred_read_count.second;
						}
					}
					{
						std::lock_guard<std::mutex> lock(output_numbers_mut);
						for (const auto &seqGCount : allSeqGCounts) {
							output_numbers.addRow(
								indexKey,
								sim_sample.mid_name_,
								sampleKey,
								sim_sample.full_haps_.size(),
								sim_sample.sample_total_depth_,
								target.first,
								allSeqGCounts.size(),
								read_depths_per_target.at(target.first),
								seqGCount.seqBase_.name_,
								sim_sample.starting_template_,
								totalGenomes,
								seqGCount.seqBase_.frac_,
								seqGCount.genomeCnt_,
								pcred_read_counts[seqGCount.seqBase_.name_],
								total_pcr_count,
								total_pcr_count - chimera_pcr_count
								);
						}
					}
					// per source counts, errors, and chimera parents
					{
						SampleTargetSummary summary;
						summary.index_ = indexKey;
						summary.mid_ = sim_sample.mid_name_;
						summary.sample_name_ = sampleKey;
						summary.target_ = target.first;
						summary.haps_amplified_ = seqGCounts.size();
						summary.starting_genomes_ = totalGenomes;
						summary.finished_copies_ = simFinishedCopies;
						summary.partial_copies_ = simPartialCopies;
						summary.switch_plateau_round_ = simSwitchPlateauRound;
						summary.final_product_ = simFinalProduct;
						summary.total_reads_ = total_pcr_count;
						std::vector<ChimeraRecord> current_chimera_records;
						auto add_source_row = [&](const std::string & source_name, const SourceStats & stats, std::vector<std::string> && chimera_cols) {
							summary.pcr_reads_with_errors_ += stats.reads_with_pcr_errors_;
							summary.pcr_error_bases_ += stats.pcr_error_bases_;
							summary.pcr_bases_ += stats.reads_ * source_seqs.at(source_name).size();
							summary.seq_reads_with_errors_ += stats.reads_with_seq_errors_;
							summary.seq_error_bases_ += stats.seq_error_bases_;
							summary.seq_bases_ += stats.seq_bases_;
							const bool is_chimera = !chimera_cols.empty();
							std::string chimera_round = "NA";
							if (is_chimera) {
								auto round_start = source_name.find("[PCRRound=");
								if (std::string::npos != round_start) {
									round_start += std::string("[PCRRound=").size();
									chimera_round = source_name.substr(round_start, source_name.find(']', round_start) - round_start);
								}
							} else {
								chimera_cols = VecStr{"NA", "NA", "NA", "NA", "NA"};
							}
							std::lock_guard<std::mutex> lock(output_numbers_mut);
							source_read_counts.addRow(
								indexKey, sim_sample.mid_name_, sampleKey, target.first,
								source_name, njh::boolToStr(is_chimera), chimera_round,
								stats.reads_, 0 == total_pcr_count ? 0.0 : static_cast<double>(stats.reads_) / total_pcr_count,
								source_seqs.at(source_name).size(), stats.reads_with_pcr_errors_, stats.pcr_error_bases_,
								stats.reads_with_seq_errors_, stats.seq_error_bases_, stats.seq_bases_,
								chimera_cols[0], chimera_cols[1], chimera_cols[2], chimera_cols[3], chimera_cols[4]);
						};
						for (const auto & seqGCount : seqGCounts) {
							add_source_row(seqGCount.seqBase_.name_, source_stats[seqGCount.seqBase_.name_], {});
						}
						for (const auto & source : source_stats) {
							if (!njh::beginsWith(source.first, "Chi.")) {
								continue;
							}
							++summary.chimera_seqs_;
							summary.chimeric_reads_ += source.second.reads_;
							// best two-parent explanation among the input haplotypes, a single crossover with the 5' end from one and 3' end from the other,
							// picking the pair with the most reads in its less abundant parent
							const auto & chimera_seq = source_seqs.at(source.first);
							std::string best_5prime;
							std::string best_3prime;
							uint64_t best_weaker_reads = 0;
							bool found_pair = false;
							for (const auto & parent_5prime : seqGCounts) {
								const auto prefix_len = common_prefix_len(chimera_seq, parent_5prime.seqBase_.seq_);
								for (const auto & parent_3prime : seqGCounts) {
									if (parent_5prime.seqBase_.name_ == parent_3prime.seqBase_.name_) {
										continue;
									}
									if (prefix_len + common_suffix_len(chimera_seq, parent_3prime.seqBase_.seq_) >= chimera_seq.size()) {
										const auto weaker_reads = std::min(source_stats[parent_5prime.seqBase_.name_].reads_, source_stats[parent_3prime.seqBase_.name_].reads_);
										if (!found_pair || weaker_reads > best_weaker_reads) {
											found_pair = true;
											best_weaker_reads = weaker_reads;
											best_5prime = parent_5prime.seqBase_.name_;
											best_3prime = parent_3prime.seqBase_.name_;
										}
									}
								}
							}
							ChimeraRecord record;
							record.target_ = target.first;
							record.reads_ = source.second.reads_;
							record.frac_ = 0 == total_pcr_count ? 0.0 : static_cast<double>(source.second.reads_) / total_pcr_count;
							record.two_parent_ = found_pair;
							VecStr chimera_cols;
							if (found_pair) {
								++summary.chimera_seqs_two_parent_;
								record.weaker_parent_ratio_ = static_cast<double>(best_weaker_reads) / source.second.reads_;
								chimera_cols = VecStr{best_5prime, best_3prime,
									njh::pasteAsStr(source_stats[best_5prime].reads_),
									njh::pasteAsStr(source_stats[best_3prime].reads_),
									njh::pasteAsStr(record.weaker_parent_ratio_)};
							} else {
								// made from more than two input haplotypes (a chimera of a chimera)
								chimera_cols = VecStr{"multiple", "multiple", "NA", "NA", "NA"};
							}
							current_chimera_records.emplace_back(record);
							add_source_row(source.first, source.second, std::move(chimera_cols));
						}
						std::lock_guard<std::mutex> lock(output_numbers_mut);
						sample_target_summaries.emplace_back(summary);
						addOtherVec(chimera_records, current_chimera_records);
					}
				}
			}
		}
	};

	njh::stopWatch simulationWatch;
	njh::concurrent::runVoidFunctionThreaded(sim_sample_func, numThreads);
	{
		table phase_timings_tab(VecStr{"phase", "seconds", "frac_of_simLibFast"});
		const double simLibFastTotal = phase_seconds.count("simLibFast_total") ? phase_seconds.at("simLibFast_total") : 0;
		phase_timings_tab.addRow("wall_simulation", simulationWatch.totalTime(), std::numeric_limits<double>::quiet_NaN());
		for (const auto & phase : phase_seconds) {
			phase_timings_tab.addRow(phase.first, phase.second, simLibFastTotal > 0 ? phase.second / simLibFastTotal : std::numeric_limits<double>::quiet_NaN());
		}
		OutputStream phase_timings_out(njh::files::make_path(setUp.pars_.directoryName_, "phase_timings.tsv"));
		phase_timings_tab.outPutContents(phase_timings_out, "\t");
	}
	//std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;

	output_numbers.sortTable("sample_name", false);
	OutputStream out_sample_counts(njh::files::make_path(setUp.pars_.directoryName_, "output_numbers.tsv"));
	output_numbers.outPutContents(out_sample_counts, "\t");

	// error and chimera reports
	{
		source_read_counts.sortTable("sample_name", "target", false);
		OutputStream source_read_counts_out(njh::files::make_path(setUp.pars_.directoryName_, "source_read_counts.tsv"));
		source_read_counts.outPutContents(source_read_counts_out, "\t");
	}
	auto safe_div = [](double num, double denom) -> double {
		return 0 == denom ? 0.0 : num / denom;
	};
	// quantile, type 7 (R's default)
	auto quantile = [](std::vector<double> vals, double prob) -> double {
		if (vals.empty()) {
			return std::numeric_limits<double>::quiet_NaN();
		}
		njh::sort(vals);
		const double h = (vals.size() - 1) * prob;
		const auto lower = static_cast<size_t>(std::floor(h));
		const auto upper = std::min(lower + 1, vals.size() - 1);
		return vals[lower] + (h - lower) * (vals[upper] - vals[lower]);
	};
	auto frac_above = [](const std::vector<double> & vals, double cutoff) -> double {
		if (vals.empty()) {
			return std::numeric_limits<double>::quiet_NaN();
		}
		return static_cast<double>(std::count_if(vals.begin(), vals.end(), [&cutoff](double val) { return val > cutoff; })) / vals.size();
	};
	{
		table sample_target_summary_tab(VecStr{
			"index", "mid", "sample_name", "target", "haps_amplified", "starting_genomes",
			"effective_switch_rate", "switch_plateau_round", "final_product",
			"total_reads", "chimeric_reads", "frac_chimeric", "chimera_seqs", "chimera_seqs_two_parent",
			"reads_with_pcr_errors", "frac_reads_with_pcr_errors", "pcr_error_bases", "pcr_error_rate",
			"reads_with_seq_errors", "frac_reads_with_seq_errors", "seq_error_bases", "seq_error_rate"
		});
		for (const auto & summary : sample_target_summaries) {
			sample_target_summary_tab.addRow(
				summary.index_, summary.mid_, summary.sample_name_, summary.target_, summary.haps_amplified_, summary.starting_genomes_,
				// copy-weighted switch rate, NaN when chimeras weren't modelled (a single haplotype, or no shared stretches)
				0 == summary.finished_copies_ + summary.partial_copies_ ? std::numeric_limits<double>::quiet_NaN() : safe_div(summary.partial_copies_, summary.finished_copies_ + summary.partial_copies_),
				summary.switch_plateau_round_,
				// all strands of the target after the last round, 0 when chimeras weren't modelled
				summary.final_product_,
				summary.total_reads_, summary.chimeric_reads_, safe_div(summary.chimeric_reads_, summary.total_reads_),
				summary.chimera_seqs_, summary.chimera_seqs_two_parent_,
				summary.pcr_reads_with_errors_, safe_div(summary.pcr_reads_with_errors_, summary.total_reads_),
				summary.pcr_error_bases_, safe_div(summary.pcr_error_bases_, summary.pcr_bases_),
				summary.seq_reads_with_errors_, safe_div(summary.seq_reads_with_errors_, summary.total_reads_),
				summary.seq_error_bases_, safe_div(summary.seq_error_bases_, summary.seq_bases_));
		}
		sample_target_summary_tab.sortTable("sample_name", "target", false);
		OutputStream sample_target_summary_out(njh::files::make_path(setUp.pars_.directoryName_, "sample_target_summary.tsv"));
		sample_target_summary_tab.outPutContents(sample_target_summary_out, "\t");
	}
	{
		struct TargetTotals {
			uint32_t sample_targets_{0};
			uint32_t multi_hap_sample_targets_{0};
			uint32_t multi_hap_sample_targets_with_chimeras_{0};
			uint64_t total_reads_{0};
			uint64_t multi_hap_reads_{0};
			uint64_t chimeric_reads_{0};
			uint64_t chimera_seqs_{0};
			uint64_t pcr_error_bases_{0};
			uint64_t pcr_bases_{0};
			uint64_t seq_error_bases_{0};
			uint64_t seq_bases_{0};
			std::vector<double> multi_hap_frac_chimeric_;
			std::vector<double> weaker_parent_ratios_;
		};
		std::map<std::string, TargetTotals> target_totals;
		TargetTotals all_totals;
		for (const auto & summary : sample_target_summaries) {
			for (auto * totals : std::vector<TargetTotals *>{&target_totals[summary.target_], &all_totals}) {
				++totals->sample_targets_;
				totals->total_reads_ += summary.total_reads_;
				totals->chimeric_reads_ += summary.chimeric_reads_;
				totals->chimera_seqs_ += summary.chimera_seqs_;
				totals->pcr_error_bases_ += summary.pcr_error_bases_;
				totals->pcr_bases_ += summary.pcr_bases_;
				totals->seq_error_bases_ += summary.seq_error_bases_;
				totals->seq_bases_ += summary.seq_bases_;
				if (summary.haps_amplified_ > 1 && summary.total_reads_ > 0) {
					++totals->multi_hap_sample_targets_;
					totals->multi_hap_reads_ += summary.total_reads_;
					totals->multi_hap_frac_chimeric_.emplace_back(safe_div(summary.chimeric_reads_, summary.total_reads_));
					if (summary.chimeric_reads_ > 0) {
						++totals->multi_hap_sample_targets_with_chimeras_;
					}
				}
			}
		}
		std::vector<double> all_ratios;
		std::vector<double> passing_ratios;
		std::vector<double> failing_ratios;
		uint64_t two_parent_chimera_seqs = 0;
		for (const auto & record : chimera_records) {
			if (!record.two_parent_) {
				continue;
			}
			++two_parent_chimera_seqs;
			target_totals[record.target_].weaker_parent_ratios_.emplace_back(record.weaker_parent_ratio_);
			all_ratios.emplace_back(record.weaker_parent_ratio_);
			if (record.reads_ > chimeraSummaryMinReads && record.frac_ >= chimeraSummaryMinFrac) {
				passing_ratios.emplace_back(record.weaker_parent_ratio_);
			} else {
				failing_ratios.emplace_back(record.weaker_parent_ratio_);
			}
		}
		table target_summary_tab(VecStr{
			"target", "sample_targets", "multi_hap_sample_targets", "multi_hap_sample_targets_with_chimeras",
			"total_reads", "multi_hap_reads", "chimeric_reads", "frac_chimeric_multi_hap", "median_sample_frac_chimeric_multi_hap",
			"chimera_seqs", "two_parent_chimera_seqs", "median_weaker_parent_to_chimera_ratio",
			"pcr_error_rate", "seq_error_rate"
		});
		for (const auto & totals : target_totals) {
			target_summary_tab.addRow(
				totals.first, totals.second.sample_targets_, totals.second.multi_hap_sample_targets_, totals.second.multi_hap_sample_targets_with_chimeras_,
				totals.second.total_reads_, totals.second.multi_hap_reads_, totals.second.chimeric_reads_,
				safe_div(totals.second.chimeric_reads_, totals.second.multi_hap_reads_),
				quantile(totals.second.multi_hap_frac_chimeric_, 0.5),
				totals.second.chimera_seqs_, totals.second.weaker_parent_ratios_.size(),
				quantile(totals.second.weaker_parent_ratios_, 0.5),
				safe_div(totals.second.pcr_error_bases_, totals.second.pcr_bases_),
				safe_div(totals.second.seq_error_bases_, totals.second.seq_bases_));
		}
		OutputStream target_summary_out(njh::files::make_path(setUp.pars_.directoryName_, "target_summary.tsv"));
		target_summary_tab.outPutContents(target_summary_out, "\t");

		uint64_t pcr_reads_with_errors = 0;
		uint64_t seq_reads_with_errors = 0;
		for (const auto & summary : sample_target_summaries) {
			pcr_reads_with_errors += summary.pcr_reads_with_errors_;
			seq_reads_with_errors += summary.seq_reads_with_errors_;
		}
		table simulation_summary_tab(VecStr{"metric", "value"});
		simulation_summary_tab.addRow("pcr_efficiency", pcrEfficiency);
		simulation_summary_tab.addRow("template_switch_rate", templateSwitchRate);
		simulation_summary_tab.addRow("template_switch_scaling", PCRSimulator::templateSwitchScalingToStr(templateSwitchScaling));
		simulation_summary_tab.addRow("switch_plateau", switchPlateau);
		simulation_summary_tab.addRow("target_product_max", targetProductMax);
		simulation_summary_tab.addRow("switch_exponent", switchExponent);
		simulation_summary_tab.addRow("sample_targets", all_totals.sample_targets_);
		simulation_summary_tab.addRow("multi_hap_sample_targets", all_totals.multi_hap_sample_targets_);
		simulation_summary_tab.addRow("total_reads", all_totals.total_reads_);
		simulation_summary_tab.addRow("frac_reads_with_pcr_errors", safe_div(pcr_reads_with_errors, all_totals.total_reads_));
		simulation_summary_tab.addRow("pcr_error_rate_per_base", safe_div(all_totals.pcr_error_bases_, all_totals.pcr_bases_));
		simulation_summary_tab.addRow("frac_reads_with_seq_errors", safe_div(seq_reads_with_errors, all_totals.total_reads_));
		simulation_summary_tab.addRow("seq_error_rate_per_base", safe_div(all_totals.seq_error_bases_, all_totals.seq_bases_));
		simulation_summary_tab.addRow("chimeric_reads", all_totals.chimeric_reads_);
		simulation_summary_tab.addRow("frac_chimeric_all_reads", safe_div(all_totals.chimeric_reads_, all_totals.total_reads_));
		simulation_summary_tab.addRow("frac_chimeric_multi_hap_reads", safe_div(all_totals.chimeric_reads_, all_totals.multi_hap_reads_));
		simulation_summary_tab.addRow("median_sample_target_frac_chimeric_multi_hap", quantile(all_totals.multi_hap_frac_chimeric_, 0.5));
		simulation_summary_tab.addRow("frac_multi_hap_sample_targets_with_chimeras", safe_div(all_totals.multi_hap_sample_targets_with_chimeras_, all_totals.multi_hap_sample_targets_));
		simulation_summary_tab.addRow("chimera_seqs", all_totals.chimera_seqs_);
		simulation_summary_tab.addRow("two_parent_chimera_seqs", two_parent_chimera_seqs);
		for (const auto & ratios : std::vector<std::pair<std::string, const std::vector<double> *>>{
			     {"all", &all_ratios},
			     {njh::pasteAsStr("passing_gt", chimeraSummaryMinReads, "reads_ge", chimeraSummaryMinFrac, "frac"), &passing_ratios},
			     {"failing", &failing_ratios}}) {
			const std::string prefix = "weaker_parent_to_chimera_ratio_" + ratios.first + "_";
			simulation_summary_tab.addRow(prefix + "n", ratios.second->size());
			for (const auto prob : std::vector<double>{0.05, 0.25, 0.5, 0.75, 0.95}) {
				simulation_summary_tab.addRow(prefix + "q" + njh::leftPadNumStr<uint32_t>(std::round(prob * 100), 99), quantile(*ratios.second, prob));
			}
			for (const auto cutoff : std::vector<double>{2, 8, 16}) {
				simulation_summary_tab.addRow(prefix + "frac_gt" + njh::pasteAsStr(cutoff), frac_above(*ratios.second, cutoff));
			}
		}
		OutputStream simulation_summary_out(njh::files::make_path(setUp.pars_.directoryName_, "simulation_summary.tsv"));
		simulation_summary_tab.outPutContents(simulation_summary_out, "\t");
	}

	return 0;
}

} // namespace njhseq

