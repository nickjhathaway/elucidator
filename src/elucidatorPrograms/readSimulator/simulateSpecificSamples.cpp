//
// Created by Nicholas H	athaway on 3/31/25.
//

#include "readSimulatorRunner.hpp"
#include "elucidator/simulation.h"
#include <SeekDeep/objects/TarAmpSetupUtils/PrimersAndMids.hpp>

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
			ss << __PRETTY_FUNCTION__ << ", error " << message << "\n";
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
	uint32_t pairedEndLength = std::numeric_limits<uint32_t>::max();


	bfs::path population_haps_fnp;
	VecStr population_haps_required_cols{"target","seq","hap_id"};
	bfs::path sample_set_up_fnp;
	VecStr sample_set_up_required_cols{"sample_name","starting_template","sample_total_depth","target","hap_id","within_sample_fullhap_id","relative_abundance"};
	bfs::path target_read_count_factor_fnp;
	VecStr target_read_count_factor_required_cols{"target", "factor"};

	bfs::path idFile;
	bfs::path illuminaProfileDir;

	setUp.processVerbose();
	setUp.processDebug();

	setUp.setOption(illuminaProfileDir, "--illuminaProfileDir", "Illumina Profile Dir", true);
	setUp.setOption(idFile, "--idFile", "Primer MID Fnp", true);
	setUp.setOption(population_haps_fnp, "--population_haps_fnp", njh::pasteAsStr("population haplotype sequences file, required columns: ", njh::conToStr(sample_set_up_required_cols, ",")), true);
	setUp.setOption(sample_set_up_fnp, "--sample_set_up_fnp", njh::pasteAsStr("sample setup file, required columns: ", njh::conToStr(sample_set_up_required_cols, ",")), true);
	setUp.setOption(target_read_count_factor_fnp, "--target_read_count_factor_fnp", njh::pasteAsStr("target read count factor, used to adjust read depths based on known biases, required columns: ", njh::conToStr(target_read_count_factor_required_cols, ",")), true);


	setUp.setOption(pairedEndLength, "--pairedEndLength", "Paired End Length");
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
	std::unordered_map<std::string, double> read_depth_biases;
	if (target_read_count_factor_fnp.empty()) {
		for (const auto & p : primers.targets_) {
			read_depth_biases[p.first] = 1.0/static_cast<double>(primers.targets_.size());
		}
	} else {
		table target_read_count_factor(target_read_count_factor_fnp, "\t", true);
		target_read_count_factor.checkForColumnsThrow(target_read_count_factor_required_cols, __PRETTY_FUNCTION__);
		double total_factor = 0.0;
		for (const auto & row : target_read_count_factor) {
			auto factor = njh::StrToNumConverter::stoToNum<double>(row[target_read_count_factor.getColPos("factor")]);
			read_depth_biases[row[target_read_count_factor.getColPos("target")]] = factor;
			total_factor += factor;
		}
		for (auto & tar_factor : read_depth_biases) {
			tar_factor.second = tar_factor.second / total_factor;
		}
		auto target_names_in_biases = njh::getSetOfMapKeys(read_depth_biases);
		check_for_level_match_up(target_names_in_biases, target_read_count_factor_fnp,
		                         target_names, idFile,
		                         __PRETTY_FUNCTION__);
	}
	// population haplotypes
	table population_haps(population_haps_fnp, "\t", true);
	population_haps.checkForColumnsThrow(population_haps_required_cols, __PRETTY_FUNCTION__);
	auto population_haps_target_names = njh::vecToSet(population_haps.getColumnLevels(population_haps.getColPos("target")));
	check_for_level_match_up(population_haps_target_names, population_haps_fnp,
												 target_names, idFile,
												 __PRETTY_FUNCTION__);
	std::unordered_map<std::string, std::unordered_map<uint32_t, seqInfo>> population_haps_by_target;
	for (const auto & row : population_haps) {
		auto target_name = row[population_haps.getColPos("target")];
		auto hap_id = njh::StrToNumConverter::stoToNum<uint32_t>(row[population_haps.getColPos("hap_id")]);
		if (njh::in(hap_id, population_haps_by_target[target_name])) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "already have hap_id " << hap_id << " for target " << target_name << "\n";
			throw std::runtime_error{ss.str()};
		}
	}
	// sample set up
	table sample_setup(sample_set_up_fnp, "\t", true);
	sample_setup.checkForColumnsThrow(sample_set_up_required_cols, __PRETTY_FUNCTION__);
	auto sample_target_names = njh::vecToSet(sample_setup.getColumnLevels(sample_setup.getColPos("target")));
	check_for_level_match_up(sample_target_names, sample_set_up_fnp,
												 target_names, idFile,
												 __PRETTY_FUNCTION__);

	return 0;
}

} // namespace njhseq

