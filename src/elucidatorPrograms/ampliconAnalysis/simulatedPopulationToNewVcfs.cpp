//
// Created by Nicholas Hathaway on 4/6/26.
//


#include <njhseq/IO/OutputStream.hpp>
#include <njhseq/objects/dataContainers/tables/TableReader.hpp>
#include <njhseq/objects/BioDataObject/GenomicRegion.hpp>
#include <njhseq/objects/BioDataObject/reading.hpp>
#include <njhseq/objects/Gene/VCFOutput.hpp>

#include "ampliconAnalysisRunner.hpp"


namespace njhseq {
int ampliconAnalysisRunner::simulatedPopulationToNewVcfs(const njh::progutils::CmdArgs &inputCommands) {
	bfs::path input_vcf_fnp;
	bfs::path simulated_population_fnp;
	bool reset_gt_fields = false;
	uint32_t reset_gt_fields_ploidy = 2;
	ampliconAnalysisSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(reset_gt_fields, "--reset_gt_fields", "reset_gt_fields");
	setUp.setOption(reset_gt_fields_ploidy, "--reset_gt_fields_ploidy", "reset_gt_fields_ploidy");

	setUp.setOption(input_vcf_fnp, "--input_vcf_fnp", "input_vcf_fnp", true);
	setUp.setOption(simulated_population_fnp, "--simulated_population_fnp", "simulated_population_fnp", true);
	setUp.processDirectoryOutputName(njh::pasteAsStr(bfs::basename(simulated_population_fnp), "_new_vcfs"), true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_ );
	VCFOutput input_vcf = VCFOutput::readInHeader(input_vcf_fnp);
	{
		InputStream in(input_vcf_fnp);
		input_vcf.addInRecordsFromFile(in);
	}
	if (reset_gt_fields) {
		input_vcf.allAddGTFields(reset_gt_fields_ploidy);
	}
	auto simulated_population = njh::json::parseFile(simulated_population_fnp.string());
	uint32_t pop_index = 0;
	for (const auto & pop : simulated_population) {
		std::string pop_name = njh::leftPadNumStr(pop_index, simulated_population.size());
		++pop_index;
		VCFOutput strains_separately = input_vcf;
		//clean up vcf for output
		//clear samples
		strains_separately.samples_.clear();
		for (auto & rec : strains_separately.records_) {
			rec.sampleFormatInfos_.clear();
		}
		//INFO/AC,INFO/AF,INFO/AN,^FORMAT/GT,FORMAT/AD
		strains_separately.otherHeaderMetaFields_.clear();
		//format entries to remove
		VecStr format_to_remove;
		for (const auto & format : strains_separately.formatEntries_) {
			if (!njh::in(format.first, {"GT", "AD"})) {
				format_to_remove.emplace_back(format.first);
			}
		}
		for (const auto & format : format_to_remove) {
			strains_separately.formatEntries_.erase(format);
		}
		//info entries to remove
		VecStr info_to_remove;
		for (const auto & info : strains_separately.infoEntries_) {
			if (!njh::in(info.first, {"AC", "AF", "AN"})) {
				info_to_remove.emplace_back(info.first);
			}
		}
		for (const auto & info : info_to_remove) {
			strains_separately.infoEntries_.erase(info);
		}


		std::unordered_map<uint32_t, std::string> ancestral_population_index;
		for (const auto & ancestral_pop_index : pop["ancestral_indexes"]) {
			ancestral_population_index[ancestral_pop_index["index"].asUInt() ] = ancestral_pop_index["ancestral_genotype"].asString();
		}
		uint32_t sample_index = 0;
		for (const auto & sim_sample : pop["simulated_samples"]) {
			auto sample_name = njh::leftPadNumStr(sample_index, pop["simulated_samples"].size());
			++sample_index;
			uint32_t genotype_index = 0;
			for (const auto & genotype : sim_sample["genotypes"]) {
				auto genotype_name = njh::leftPadNumStr(genotype_index, sim_sample["genotypes"].size());
				++genotype_index;
				std::string output_sample_name = njh::pasteAsStr(sample_name, "_", genotype_name);
				strains_separately.samples_.emplace_back(output_sample_name);
				for (const auto & segments : genotype["segments"]) {
					Bed3RecordCore current_segment(segments["chrom"].asString(), segments["start"].asUInt(), segments["end"].asUInt());
					auto ancestral_genotype = ancestral_population_index[segments["index"].asUInt()];
					//assumes input is sorted
					for (const auto & recPos : iter::enumerate(strains_separately.records_)) {
						auto zero_based_pos = recPos.element.pos_ - 1;
						if (recPos.element.chrom_ > current_segment.chrom_ || (recPos.element.chrom_ == current_segment.chrom_ && zero_based_pos >= current_segment.chromEnd_)) {
							break;
						}
						if (recPos.element.chrom_ == current_segment.chrom_ && zero_based_pos >= current_segment.chromStart_ && zero_based_pos < current_segment.chromEnd_) {
							MetaDataInName output_sample_info;
							// std::cout << __FILE__ << " " << __LINE__ << std::endl;
							// std::cout << "ancestral_genotype: " << ancestral_genotype << std::endl;
							// std::cout << "input_vcf.records_[recPos.index].sampleFormatInfos_[ancestral_genotype]" << njh::json::toJson(input_vcf.records_[recPos.index].sampleFormatInfos_[ancestral_genotype].meta_) << std::endl;
							output_sample_info.addMeta("GT", input_vcf.records_[recPos.index].sampleFormatInfos_[ancestral_genotype].getMeta("GT"));
							output_sample_info.addMeta("AD", input_vcf.records_[recPos.index].sampleFormatInfos_[ancestral_genotype].getMeta("AD"));
							// std::cout << __FILE__ << " " << __LINE__ << std::endl;
							strains_separately.records_[recPos.index].sampleFormatInfos_.emplace(output_sample_name, output_sample_info);
						}
					}
					// std::cout << segments["chrom"].asString() << " " << segments["start"].asUInt() << " " << segments["end"].asUInt() << " " << segments["index"].asUInt() << std::endl;
				}
			}
		}
		OutputStream pop_out_vcf(njh::files::make_path(setUp.pars_.directoryName_, njh::pasteAsStr(pop_name, "_strains_separately.vcf.gz")));
		strains_separately.allAutoAdd_AN_AC_AF_InfoFields();
		strains_separately.writeOutFixedAndSampleMeta(pop_out_vcf);
		//std::cout << njh::json::toJson(ancestral_population_index) << std::endl;
	}

	if (setUp.pars_.debug_){
		OutputStream out(njh::files::make_path(setUp.pars_.directoryName_, "out.vcf.gz"));
		input_vcf.writeOutFixedAndSampleMeta(out);
	}

	return 0;
}
} //namespace njhseq

