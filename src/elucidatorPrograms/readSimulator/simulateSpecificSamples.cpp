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
	double pcrEfficiency = 0.85;
	bool keepPCRSeqs = false;
	uint32_t chimeraBasesIn = 5;
	uint32_t templateCap = 500000000;
	bool noChimeras = false;
	double finalReadAmountSDFrac = 0.1;
	LibrarySetup::SimLibrarySetupPars sim_lib_pars;

	bfs::path population_haps_fnp;
	VecStr population_haps_required_cols{"target","seq","hap_id"};
	bfs::path sample_set_up_fnp;
	VecStr sample_set_up_required_cols{"sample_name","starting_template","sample_total_depth","target","hap_id","within_sample_fullhap_id","relative_abundance", "pcr_rounds"};
	bfs::path target_read_count_factor_fnp;
	VecStr target_read_count_factor_required_cols{"target", "factor", "factor_stddev"};
	bfs::path idFile;
	bfs::path illuminaProfileDir;

	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(sim_lib_pars.noAddPrimers_, "--noAddPrimers", "don't add primers");
	setUp.setOption(sim_lib_pars.addBluntEndingArtifact_, "--addBluntEndingArtifact", "add Blunt Ending Artifact");
	setUp.setOption(sim_lib_pars.addReverseComplement_, "--addReverseComplement", "add Reverse Complement");

	setUp.setOption(illuminaProfileDir, "--illuminaProfileDir", "Illumina Profile Dir", true);
	setUp.setOption(idFile, "--idFile", "Primer MID Fnp", true);
	setUp.setOption(population_haps_fnp, "--population_haps_fnp", njh::pasteAsStr("population haplotype sequences file, required columns: ", njh::conToStr(sample_set_up_required_cols, ",")), true);
	setUp.setOption(sample_set_up_fnp, "--sample_set_up_fnp", njh::pasteAsStr("sample setup file, required columns: ", njh::conToStr(sample_set_up_required_cols, ",")), true);
	setUp.setOption(target_read_count_factor_fnp, "--target_read_count_factor_fnp", njh::pasteAsStr("target read count factor, used to adjust read depths based on known biases, required columns: ", njh::conToStr(target_read_count_factor_required_cols, ",")), true);


	setUp.setOption(sim_lib_pars.pairedEndLength_, "--pairedEndLength", "Paired End Length", !singleEnd);
	setUp.setOption(finalReadAmountSDFrac, "--finalReadAmountSDFrac", "final Read Amount SD Frac", njh::progutils::ProgramSetUp::CheckCase::GREATERZERO);
	setUp.setOption(noChimeras, "--noChimeras", "Don't simulate chimeras");
	setUp.setOption(templateCap, "--templateCap", "Template Cap");
	setUp.setOption(chimeraBasesIn, "--chimeraBasesIn", "The number of bases needed for a template to lay down");
	setUp.setOption(keepPCRSeqs, "--keepPCRSeqs", "Keep PCR Seqs");
	setUp.setOption(errorRate, "--errorRate", "Polymerase Error Rate");
	setUp.setOption(pcrEfficiency, "--pcrEfficiency", "PCR Efficiency, between 0-1, chance a product gets amplified");
	setUp.setOption(defaultPcrRounds, "--pcrRounds", "Number of PCR rounds");
	setUp.setOption(initialPcrRounds, "--initialPcrRounds", "Number of Initial rounds of PCR before sampling");
	setUp.setOption(initialPcrRoundsMapFnp, "--initialPcrRoundsTable", "Number of Initial rounds of PCR before sampling per starting template amount, columns 1)template, 2) rounds");
	setUp.setOption(numThreads, "--numThreads", "Number of Threads to Use");
	setUp.setOption(pcrNumThreads, "--pcrNumThreads", "Number of Threads to Use for PCR sim");
	setUp.setOption(singleEnd, "--singleEnd", "Single End");
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
		std::map<std::string, uint64_t> gen_read_depths_per_target(std::unordered_map<std::string, std::normal_distribution<double>> & read_depth_biases, njh::randomGenerator& rGen) const {
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
					bias_factor += read_depth_biases[hap_counts_per_target.first](rGen.mtGen_);
				}
				relative_depths[hap_counts_per_target.first] = bias_factor;
				total_bias_factor += bias_factor;
			}
			//renormalize the bias factor
			for (auto & relative_depth : relative_depths) {
				relative_depth.second /= total_bias_factor;
			}
			std::normal_distribution<double> total_depth_dist(sample_total_depth_, sample_total_depth_ * 0.10);
			for (const auto & target_bias : relative_depths) {
				auto raw_read_depth = std::round(target_bias.second * total_depth_dist(rGen.mtGen_));
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
	auto sample_only_info = sample_setup.getColumns(VecStr{"sample_name","starting_template","sample_total_depth", "pcr_rounds"}).getUniqueRows();
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
		sim_sample.pcr_rounds_ = njh::StrToNumConverter::stoToNum<uint32_t>(row[sample_only_info.getColPos("pcr_rounds")]);
		sim_samples.emplace(sample_name, sim_sample);
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
		if (njh::notIn(hap_id, population_haps_by_target[target_name])) {
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
	std::atomic_uint sampleCountAtom{0};
	std::string extension = ".fastq.gz";
	table output_numbers(VecStr{
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
	njh::concurrent::LockableQueue<std::string> sample_queue(njh::getVecOfMapKeys(sim_samples));

	std::function<void()> sim_sample_func = [&sample_queue, &sim_samples,
		&illuminaProfileDir, &sampleCountAtom, &sim_lib_pars,
		&keepPCRSeqs, &singleEnd,
		&fastqDirectory, &extension,
		&read_depth_biases,
		&primers, &setUp,
		&pcrNumThreads,
		&intErrorRate, &pcrEfficiency, &templateCap,&noChimeras, &chimeraBasesIn,
		&initialPcrRounds, &initialPcrRoundsMap,
		&population_haps_by_target,
		&output_numbers, &output_numbers_mut]() {
		std::string sampleKey;
		//std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
		njh::randObjectGen<char,uint32_t> baseRGen({'A', 'C', 'G', 'T'}, {1,1,1,1});
		RoughIlluminaSimulator simulator(illuminaProfileDir);
		njh::randomGenerator rGen;
		while (sample_queue.getVal(sampleKey)) {
			const auto & sim_sample = sim_samples.at(sampleKey);

			uint32_t sampleCount = sampleCountAtom++;
			OutOptions targetOutOpts(njh::files::make_path(fastqDirectory, njh::pasteAsStr(sampleKey, "_S", sampleCount + 1, extension)));
			OutOptions r1OutOpts(njh::files::make_path(fastqDirectory, njh::pasteAsStr(sampleKey, "_S", sampleCount + 1, "_R1_001" + extension)));
			OutOptions r2OutOpts(njh::files::make_path(fastqDirectory, njh::pasteAsStr(sampleKey, "_S", sampleCount + 1, "_R2_001" + extension)));
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

			auto read_depths_per_target = sim_sample.gen_read_depths_per_target(read_depth_biases, rGen);
			auto expected_relative_depths_per_hap_id = sim_sample.gen_relative_depths_per_hap_id();
			auto starting_templates_per_hap_id = sim_sample.gen_starting_templates_per_hap_id(rGen);
			// auto copy_number_per_hap_id = sim_sample.gen_number_per_hap_id();



			for (const auto & target : starting_templates_per_hap_id) {
				std::vector<PCRSimulator::SeqGenomeCnt> seqGCounts;
				uint32_t totalGenomes = 0;
				for (const auto & hap_id : target.second) {
					std::string seqName = njh::pasteAsStr(hap_id.first);
					std::string seq = population_haps_by_target[target.first][hap_id.first].seq_;
					if (!sim_lib_pars.noAddPrimers_) {
						seq.append(seqUtil::reverseComplement(primers.targets_.at(target.first).info_.reversePrimerRaw_, "DNA"));
						seq.insert(0, primers.targets_.at(target.first).info_.forwardPrimerRaw_);
					}
					/*
					 * @todo add back in doing MIDs and doing random bases etc
					 */
					seqInfo currentSeq(seqName, seq);

					currentSeq.frac_ = expected_relative_depths_per_hap_id[target.first][hap_id.first];
					currentSeq.cnt_ = expected_relative_depths_per_hap_id[target.first][hap_id.first];

					seqGCounts.emplace_back(currentSeq, hap_id.second);
					totalGenomes += hap_id.second;
				}

				PCRSimulator pcrSim(intErrorRate);
				pcrSim.verbose_ = setUp.pars_.verbose_;
				pcrSim.pcrEfficiency_ = pcrEfficiency;
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
				if (read_depths_per_target[target.first] > 0) {
					auto pcrSimAmounts = pcrSim.simLibFast(seqGCounts, pcrReadsOpts, read_depths_per_target[target.first], sim_sample.pcr_rounds_, currentInitial, pcrNumThreads);
					auto pcrSeqsInOpts = SeqIOOptions::genFastaIn(pcrReadsFnp);
					SeqInput pcrReader(pcrSeqsInOpts);
					pcrReader.openIn();
					seqInfo targetSeq;

					while(pcrReader.readNextRead(targetSeq)){
						MetaDataInName targetSeqMeta(targetSeq.name_);
						++pcred_read_counts[targetSeqMeta.getMeta("hap")];
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
						if(singleEnd){
							//target
							targetSeq.outPutSeq(*targetOut);
						}else{
							{
								//r1
								auto subSeq = len(targetSeq) > sim_lib_pars.pairedEndLength_? targetSeq.getSubRead(0, sim_lib_pars.pairedEndLength_) : targetSeq;
								simulator.simR1(subSeq, sim_lib_pars.pairedEndLength_).outPutFastq(*r1Out);
							}
							{
								//r2
								targetSeq.reverseComplementRead(false, true);
								auto subSeq = len(targetSeq) > sim_lib_pars.pairedEndLength_? targetSeq.getSubRead(0, sim_lib_pars.pairedEndLength_) : targetSeq;
								simulator.simR2(subSeq, sim_lib_pars.pairedEndLength_).outPutFastq(*r2Out);
							}
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
					for (const auto &seqGCount : seqGCounts) {
						output_numbers.addRow(
							sampleKey,
							sim_sample.full_haps_.size(),
							sim_sample.sample_total_depth_,
							target.first,
							seqGCounts.size(),
							read_depths_per_target[target.first],
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
			}
		}
	};

	njh::concurrent::runVoidFunctionThreaded(sim_sample_func, numThreads);
	//std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;

	output_numbers.sortTable("sample_name", false);
	OutputStream out_sample_counts(njh::files::make_path(setUp.pars_.directoryName_, "output_numbers.tsv"));
	output_numbers.outPutContents(out_sample_counts, "\t");

	return 0;
}

} // namespace njhseq

