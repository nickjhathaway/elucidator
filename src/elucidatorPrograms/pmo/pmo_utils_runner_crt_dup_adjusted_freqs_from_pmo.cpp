//
// crt_dup_adjusted_freqs_from_pmo, population frequencies of pfcrt variants adjusted for the pfcrt duplication,
// using a phased haplotype (IET 74-76 + 97L by default) as the marker of the duplicated copy
//

#include "pmo_utils_runner.hpp"

#include <njhseq/objects/BioDataObject/pmo.h>
#include <nlohmann/json.hpp>

#include <random>

namespace njhseq {

namespace {

/**@brief The variants of one microhaplotype at the positions being reported
 *
 */
struct MhapCall {
	bool typed_{false};     /**< true if every position was typed*/
	bool isDup_{false};     /**< true if it has every amino acid of the duplication marker haplotype*/
	std::vector<std::string> aas_; /**< the amino acid at each position, in the order of the positions*/
};

/**@brief The within sample frequencies and the copies per parasite of each haplotype for a specimen
 *
 */
struct SpecimenCall {
	uint32_t specimenId_{std::numeric_limits<uint32_t>::max()};
	uint32_t reads_{0};
	double dupWsaf_{0};           /**< within sample frequency of the duplication marker haplotype(s)*/
	bool adjusted_{false};        /**< true if the duplication adjustment was applied*/
	bool singleCopyMarker_{false}; /**< true if the marker frequency is too high for a one to one duplication*/
	std::map<std::string, double> wsaf_;
	std::map<std::string, double> copies_; /**< copies per parasite, non-marker haplotypes sum to 1 if adjusted*/
	double totalCopies_{1};
};

/**@brief Adjust a specimen's within sample frequencies for the duplication
 *
 * Every parasite, duplicated or not, has one non-marker copy, so the non-marker reads count the parasites.
 * Dividing by them turns read frequencies into copies per parasite, marker copies per parasite is the fraction of parasites duplicated (at most 1).
 *
 * @param call the specimen, with wsaf_ set
 * @param isDup whether each haplotype is the duplication marker
 * @param maxDupWsaf above this the marker can't come from a one to one duplication and it's left as single copy
 */
void adjustForDup(SpecimenCall & call, const std::map<std::string, bool> & isDup, double maxDupWsaf) {
	call.dupWsaf_ = 0;
	for (const auto & hap : call.wsaf_) {
		if (isDup.at(hap.first)) {
			call.dupWsaf_ += hap.second;
		}
	}
	call.singleCopyMarker_ = call.dupWsaf_ > maxDupWsaf;
	call.adjusted_ = call.dupWsaf_ > 0 && !call.singleCopyMarker_;
	call.copies_ = call.wsaf_;
	call.totalCopies_ = 1;
	if (call.adjusted_) {
		const double dupFraction = std::min(1.0, call.dupWsaf_ / (1 - call.dupWsaf_));
		call.totalCopies_ = 1 + dupFraction;
		for (auto & hap : call.copies_) {
			hap.second = isDup.at(hap.first) ? hap.second / call.dupWsaf_ * dupFraction : hap.second / (1 - call.dupWsaf_);
		}
	}
}

std::set<std::string> readNames(const bfs::path & fnp) {
	std::set<std::string> ret;
	InputStream in(fnp);
	std::string line;
	while (njh::files::crossPlatGetline(in, line)) {
		if (!line.empty()) {
			ret.emplace(line);
		}
	}
	return ret;
}

std::string specimenMetaValue(const nlohmann::json & specimen, const std::string & field) {
	//:: joins fields
	std::string ret;
	for (const auto & tok : iter::enumerate(njh::tokenizeString(field, "::"))) {
		if (tok.index > 0) {
			ret += "::";
		}
		if (!specimen.contains(tok.element) || specimen.at(tok.element).is_null()) {
			ret += "NA";
		} else if (specimen.at(tok.element).is_string()) {
			ret += specimen.at(tok.element).get<std::string>();
		} else {
			ret += specimen.at(tok.element).dump();
		}
	}
	return ret;
}

std::pair<double, double> percentileInterval(std::vector<double> vals) {
	if (vals.empty()) {
		return {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
	}
	std::sort(vals.begin(), vals.end());
	auto quantile = [&vals](double p) {
		const double pos = p * static_cast<double>(vals.size() - 1);
		const auto low = static_cast<size_t>(std::floor(pos));
		const auto high = std::min(low + 1, vals.size() - 1);
		return vals[low] + (pos - static_cast<double>(low)) * (vals[high] - vals[low]);
	};
	return {quantile(0.025), quantile(0.975)};
}

} // namespace

int PMOUtilsRunner::crt_dup_adjusted_freqs_from_pmo(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path pmoFnp;
	std::string transcriptId = "PF3D7_0709000.1";
	VecStr dupHaplotypeStrs{"74:I", "75:E", "76:T", "97:L"};
	std::vector<uint32_t> aaPositions{72, 73, 74, 75, 76, 97};
	std::string targetName;
	VecStr specimenMeta;
	bfs::path specimensFnp;
	bfs::path excludeSpecimensFnp;
	uint32_t minReads = 100;
	double minWsaf = 0.01;
	double maxDupWsaf = 0.75;
	uint32_t bootstraps = 1000;
	uint64_t seed = 1;
	bfs::path perSpecimenFnp;
	OutOptions outOpts(bfs::path(""), ".tsv");
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.description_ = "Population frequencies of pfcrt variants from a PMO adjusted for the pfcrt duplication. "
			"A duplicated parasite has a marker copy (IET 74-76 + 97L) and a non-marker copy, so its marker reads are diluted by its own non-marker reads. "
			"Every parasite has one non-marker copy, so reads are divided by the non-marker reads to get copies per parasite. "
			"Option 1 is the parasites carrying a variant over all parasites (variants of a position can sum to more than 1), "
			"option 2 is the copies with a variant over all pfcrt copies (variants of a position sum to 1). "
			"Assumes the duplication is one to one and the marker is only on duplicated copies; "
			"for option 1 of variants on the marker (e.g. 76T) also that the duplication's other copy doesn't have the variant when the reads allow it. Rough estimates.";
	setUp.examples_.emplace_back("MASTERPROGRAM SUBPROGRAM --pmo pmo_with_variants.json.gz --out crt_freqs");
	setUp.examples_.emplace_back("MASTERPROGRAM SUBPROGRAM --pmo pmo_with_variants.json.gz --specimen_meta collection_country,site,site::collection_date --exclude_specimens controls.txt --per_specimen_out per_specimen.tsv --out crt_freqs");
	setUp.setOption(pmoFnp, "--pmo", "PMO file to read, needs protein variants on its representative microhaplotypes", true, "Input");
	setUp.setOption(specimensFnp, "--specimens", "Optional, file with the specimen names to use, one per line, default is all", false, "Input");
	setUp.setOption(excludeSpecimensFnp, "--exclude_specimens", "Optional, file with specimen names to leave out (e.g. controls), one per line", false, "Input");
	setUp.setOption(specimenMeta, "--specimen_meta", "Specimen fields to group by, frequencies are reported for all specimens and for each group of each field, join fields with ::", false, "Input");
	setUp.setOption(transcriptId, "--transcript_id", "Transcript id of the protein variants", false, "Variants");
	setUp.setOption(dupHaplotypeStrs, "--dup_haplotype", "Amino acids that mark the duplicated copy, position:amino acid, all have to be on the same microhaplotype", false, "Variants");
	setUp.setOption(aaPositions, "--aa_positions", "Amino acid positions to report, 1-based", false, "Variants");
	setUp.setOption(targetName, "--target", "Target to use, default is the only target with every position typed", false, "Variants");
	setUp.setOption(minReads, "--min_reads", "Minimum reads a specimen needs on the target", false, "Filtering");
	setUp.setOption(minWsaf, "--min_wsaf", "Minimum within sample frequency of a haplotype", false, "Filtering");
	setUp.setOption(maxDupWsaf, "--max_dup_wsaf", "A marker frequency above this can't come from a one to one duplication, so the specimen is left as single copy", false, "Filtering");
	setUp.setOption(bootstraps, "--bootstraps", "Number of bootstraps over specimens for the 95% intervals, 0 for none", false, "Intervals");
	setUp.setOption(seed, "--seed", "Seed for the bootstraps", false, "Intervals");
	setUp.setOption(perSpecimenFnp, "--per_specimen_out", "Optional, file to write the frequencies and copies per parasite for each specimen", false, "Writing Output");
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	//the marker haplotype
	std::map<uint32_t, std::string> dupHaplotype;
	for (const auto & tok : dupHaplotypeStrs) {
		auto toks = njh::tokenizeString(tok, ":");
		if (2 != toks.size() || toks[0].empty() || toks[1].empty()) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error, --dup_haplotype should be position:amino acid, not " << tok << "\n";
			throw std::runtime_error{ss.str()};
		}
		dupHaplotype[njh::StrToNumConverter::stoToNum<uint32_t>(toks[0])] = toks[1];
	}
	std::set<uint32_t> positionsSet(aaPositions.begin(), aaPositions.end());
	for (const auto & aa : dupHaplotype) {
		positionsSet.emplace(aa.first);
	}
	const std::vector<uint32_t> positions(positionsSet.begin(), positionsSet.end());

	InputStream in(pmoFnp);
	const auto pmoObj = pmo::PortableMicrohaplotypeObject::from_json(nlohmann::json::parse(in));
	pmoObj.validate();

	//call the haplotypes of each target
	auto callMhap = [&](const pmo::RepresentativeMicrohaplotype & mhap) {
		MhapCall ret;
		std::map<uint32_t, std::string> aas;
		if (mhap.associated_protein_variants_.has_value()) {
			for (const auto & variant : *mhap.associated_protein_variants_) {
				if (transcriptId == variant.protein_location_.chrom_ && variant.protein_location_.alt_seq_.has_value()) {
					aas[variant.protein_location_.start_] = *variant.protein_location_.alt_seq_;
				}
			}
		}
		ret.typed_ = true;
		for (const auto & pos : positions) {
			if (njh::notIn(pos, aas)) {
				ret.typed_ = false;
				return ret;
			}
			ret.aas_.emplace_back(aas.at(pos));
		}
		ret.isDup_ = std::all_of(dupHaplotype.begin(), dupHaplotype.end(), [&aas](const auto & aa) { return aas.at(aa.first) == aa.second; });
		return ret;
	};
	std::vector<uint32_t> mhapsTargetIds;
	std::vector<MhapCall> mhapCalls;
	for (const auto & target : iter::enumerate(pmoObj.representative_microhaplotypes_.targets_)) {
		const auto & name = pmoObj.target_info_[target.element.target_id_].target_name_;
		if (!targetName.empty() && targetName != name) {
			continue;
		}
		std::vector<MhapCall> calls;
		for (const auto & mhap : target.element.microhaplotypes_) {
			calls.emplace_back(callMhap(mhap));
		}
		if (std::any_of(calls.begin(), calls.end(), [](const MhapCall & call) { return call.typed_; })) {
			mhapsTargetIds.emplace_back(target.index);
			mhapCalls = calls;
		}
	}
	if (1 != mhapsTargetIds.size()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error, " << mhapsTargetIds.size() << " targets have every position typed (" << njh::conToStr(positions, ",") << ") for "
			 << transcriptId << " in " << pmoFnp << (mhapsTargetIds.empty() ? "" : ", choose one with --target") << "\n";
		for (const auto & id : mhapsTargetIds) {
			ss << pmoObj.target_info_[pmoObj.representative_microhaplotypes_.targets_[id].target_id_].target_name_ << "\n";
		}
		throw std::runtime_error{ss.str()};
	}
	const auto mhapsTargetId = mhapsTargetIds.front();
	const auto usedTarget = pmoObj.target_info_[pmoObj.representative_microhaplotypes_.targets_[mhapsTargetId].target_id_].target_name_;

	//reads per haplotype per specimen, library samples of the same specimen are summed
	std::set<std::string> includeNames;
	std::set<std::string> excludeNames;
	if (!specimensFnp.empty()) {
		includeNames = readNames(specimensFnp);
	}
	if (!excludeSpecimensFnp.empty()) {
		excludeNames = readNames(excludeSpecimensFnp);
	}
	std::set<std::string> pmoNames;
	for (const auto & specimen : pmoObj.specimen_info_) {
		pmoNames.emplace(specimen.specimen_name_);
	}
	for (const auto & names : {includeNames, excludeNames}) {
		for (const auto & name : names) {
			if (njh::notIn(name, pmoNames)) {
				std::cerr << "warning, specimen " << name << " isn't in " << pmoFnp << std::endl;
			}
		}
	}
	std::map<uint32_t, std::map<std::string, uint32_t>> readsPerSpecimen;
	std::map<std::string, bool> isDup;
	std::map<std::string, std::vector<std::string>> aasForHap;
	uint64_t untypedReads = 0;
	for (const auto & detected : pmoObj.detected_microhaplotypes_) {
		for (const auto & librarySample : detected.library_samples_) {
			const auto specimenId = pmoObj.library_sample_info_[librarySample.library_sample_id_].specimen_id_;
			const auto & specimenName = pmoObj.specimen_info_[specimenId].specimen_name_;
			if ((!specimensFnp.empty() && njh::notIn(specimenName, includeNames)) || njh::in(specimenName, excludeNames)) {
				continue;
			}
			for (const auto & target : librarySample.target_results_) {
				if (target.mhaps_target_id_ != mhapsTargetId) {
					continue;
				}
				for (const auto & mhap : target.mhaps_) {
					const auto & call = mhapCalls[mhap.mhap_id_];
					if (!call.typed_) {
						untypedReads += mhap.reads_;
						continue;
					}
					const auto hap = njh::conToStr(call.aas_, "");
					isDup[hap] = call.isDup_;
					aasForHap[hap] = call.aas_;
					readsPerSpecimen[specimenId][hap] += mhap.reads_;
				}
			}
		}
	}

	std::vector<SpecimenCall> specimens;
	uint32_t lowReads = 0;
	for (const auto & reads : readsPerSpecimen) {
		SpecimenCall call;
		call.specimenId_ = reads.first;
		for (const auto & hap : reads.second) {
			call.reads_ += hap.second;
		}
		if (call.reads_ < minReads || 0 == call.reads_) {
			++lowReads;
			continue;
		}
		double total = 0;
		for (const auto & hap : reads.second) {
			const double wsaf = hap.second / static_cast<double>(call.reads_);
			if (wsaf >= minWsaf) {
				call.wsaf_[hap.first] = wsaf;
				total += wsaf;
			}
		}
		for (auto & hap : call.wsaf_) {
			hap.second /= total;
		}
		adjustForDup(call, isDup, maxDupWsaf);
		specimens.emplace_back(call);
	}
	std::set<std::string> haps;
	for (const auto & specimen : specimens) {
		for (const auto & hap : specimen.wsaf_) {
			haps.emplace(hap.first);
		}
	}

	//variants to report, each is a set of haplotypes, reported for every group even if it's at 0 in the group
	struct Variant {
		std::string type_;
		std::string name_;
		std::string position_;
		std::string aa_;
		std::set<std::string> haps_;
	};
	std::vector<Variant> variants;
	for (const auto & hap : haps) {
		variants.emplace_back(Variant{"haplotype", hap, "NA", "NA", {hap}});
	}
	for (const auto & pos : iter::enumerate(positions)) {
		std::map<std::string, std::set<std::string>> hapsForAA;
		for (const auto & hap : haps) {
			hapsForAA[aasForHap.at(hap)[pos.index]].emplace(hap);
		}
		for (const auto & aa : hapsForAA) {
			variants.emplace_back(Variant{"aa", njh::pasteAsStr(pos.element, aa.first), estd::to_string(pos.element), aa.first, aa.second});
		}
	}
	auto sumFor = [](const std::map<std::string, double> & vals, const std::set<std::string> & forHaps) {
		double ret = 0;
		for (const auto & hap : forHaps) {
			auto search = vals.find(hap);
			if (vals.end() != search) {
				ret += search->second;
			}
		}
		return ret;
	};

	OutputStream out(outOpts);
	out << "meta_field\tmeta_value\ttranscript_id\ttarget_name\tvariant_type\tvariant\taa_position\taa"
			<< "\tspecimens_total\tspecimens_with_variant\tspecimen_prev\tspecimens_with_dup_marker"
			<< "\traw_freq"
			<< "\tadj_freq_opt1_of_parasites\tadj_freq_opt1_low\tadj_freq_opt1_high"
			<< "\tadj_freq_opt2_of_copies\tadj_freq_opt2_low\tadj_freq_opt2_high"
			<< "\tmean_crt_copies_per_parasite" << std::endl;

	auto reportGroup = [&](const std::string & metaField, const std::string & metaValue, const std::vector<uint32_t> & group) {
		const auto n = group.size();
		std::vector<double> totals;
		double totalCopies = 0;
		uint32_t withDup = 0;
		for (const auto & pos : group) {
			totals.emplace_back(specimens[pos].totalCopies_);
			totalCopies += specimens[pos].totalCopies_;
			withDup += specimens[pos].dupWsaf_ > 0 ? 1 : 0;
		}
		//same resamples for every variant
		std::mt19937_64 gen(seed);
		std::uniform_int_distribution<size_t> pick(0, n - 1);
		std::vector<std::vector<uint32_t>> resamples(bootstraps, std::vector<uint32_t>(n));
		for (auto & resample : resamples) {
			for (auto & pos : resample) {
				pos = pick(gen);
			}
		}
		for (const auto & variant : variants) {
			std::vector<double> carriers(n);
			std::vector<double> copies(n);
			double rawSum = 0;
			double carrierSum = 0;
			double copySum = 0;
			uint32_t withVariant = 0;
			for (const auto & pos : iter::enumerate(group)) {
				const auto & specimen = specimens[pos.element];
				const auto raw = sumFor(specimen.wsaf_, variant.haps_);
				copies[pos.index] = sumFor(specimen.copies_, variant.haps_);
				carriers[pos.index] = std::min(1.0, copies[pos.index]);
				rawSum += raw;
				copySum += copies[pos.index];
				carrierSum += carriers[pos.index];
				withVariant += raw > 0 ? 1 : 0;
			}
			std::vector<double> boot1;
			std::vector<double> boot2;
			for (const auto & resample : resamples) {
				double c1 = 0;
				double c2 = 0;
				double tot = 0;
				for (const auto & pos : resample) {
					c1 += carriers[pos];
					c2 += copies[pos];
					tot += totals[pos];
				}
				boot1.emplace_back(c1 / static_cast<double>(n));
				boot2.emplace_back(c2 / tot);
			}
			const auto interval1 = percentileInterval(boot1);
			const auto interval2 = percentileInterval(boot2);
			auto intervalStr = [](double val) { return std::isnan(val) ? std::string("NA") : estd::to_string(val); };
			out << metaField << "\t" << metaValue << "\t" << transcriptId << "\t" << usedTarget
					<< "\t" << variant.type_ << "\t" << variant.name_ << "\t" << variant.position_ << "\t" << variant.aa_
					<< "\t" << n << "\t" << withVariant << "\t" << withVariant / static_cast<double>(n) << "\t" << withDup
					<< "\t" << rawSum / static_cast<double>(n)
					<< "\t" << carrierSum / static_cast<double>(n) << "\t" << intervalStr(interval1.first) << "\t" << intervalStr(interval1.second)
					<< "\t" << copySum / totalCopies << "\t" << intervalStr(interval2.first) << "\t" << intervalStr(interval2.second)
					<< "\t" << totalCopies / static_cast<double>(n) << "\n";
		}
	};

	if (!specimens.empty()) {
		std::vector<uint32_t> all(specimens.size());
		njh::iota<uint32_t>(all, 0);
		reportGroup("all", "all", all);
		std::vector<nlohmann::json> specimenJsons;
		if (!specimenMeta.empty()) {
			for (const auto & specimen : specimens) {
				specimenJsons.emplace_back(pmoObj.specimen_info_[specimen.specimenId_].to_json());
			}
		}
		for (const auto & field : specimenMeta) {
			std::map<std::string, std::vector<uint32_t>> groups;
			for (const auto & pos : all) {
				groups[specimenMetaValue(specimenJsons[pos], field)].emplace_back(pos);
			}
			for (const auto & group : groups) {
				reportGroup(field, group.first, group.second);
			}
		}
	}

	if (!perSpecimenFnp.empty()) {
		OutOptions perSpecimenOpts(perSpecimenFnp);
		perSpecimenOpts.overWriteFile_ = outOpts.overWriteFile_;
		OutputStream perSpecimenOut(perSpecimenOpts);
		perSpecimenOut << "specimen_name\treads\tdup_marker_wsaf\tdup_fraction_of_parasites\tcrt_copies_per_parasite\tstatus\thaplotype\tis_dup_marker\twsaf\tcopies_per_parasite" << std::endl;
		for (const auto & specimen : specimens) {
			for (const auto & hap : specimen.wsaf_) {
				perSpecimenOut << pmoObj.specimen_info_[specimen.specimenId_].specimen_name_ << "\t" << specimen.reads_ << "\t" << specimen.dupWsaf_
						<< "\t" << (specimen.adjusted_ ? specimen.totalCopies_ - 1 : 0) << "\t" << specimen.totalCopies_
						<< "\t" << (specimen.adjusted_ ? "adjusted" : (specimen.singleCopyMarker_ ? "marker_single_copy" : "no_marker"))
						<< "\t" << hap.first << "\t" << njh::boolToStr(isDup.at(hap.first)) << "\t" << hap.second << "\t" << specimen.copies_.at(hap.first) << "\n";
			}
		}
	}
	const auto countIf = [&specimens](const std::function<bool(const SpecimenCall &)> & pred) { return std::count_if(specimens.begin(), specimens.end(), pred); };
	std::cerr << "target: " << usedTarget << ", specimens: " << specimens.size() << ", below --min_reads: " << lowReads
			<< ", adjusted: " << countIf([](const SpecimenCall & call) { return call.adjusted_; })
			<< ", marker above --max_dup_wsaf: " << countIf([](const SpecimenCall & call) { return call.singleCopyMarker_; })
			<< ", reads on haplotypes without every position typed: " << untypedReads << std::endl;
	return 0;
}

} // namespace njhseq
