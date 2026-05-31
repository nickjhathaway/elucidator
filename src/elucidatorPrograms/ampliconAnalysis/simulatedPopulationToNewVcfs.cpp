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

	uint32_t reset_gt_fields_ploidy = 2;
  uint32_t index_start = 1;
  uint32_t fake_total_sample_depth = 1000;
  uint32_t num_threads = 1;
	ampliconAnalysisSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
  bool do_not_reset_gt_fields = false;
	setUp.setOption(do_not_reset_gt_fields, "--do_not_reset_gt_fields", "do_not_reset_gt_fields");
  bool reset_gt_fields = ! do_not_reset_gt_fields;
	setUp.setOption(reset_gt_fields_ploidy, "--reset_gt_fields_ploidy", "reset_gt_fields_ploidy");
  setUp.setOption(index_start, "--index_start", "what number to start the indexes from");
  setUp.setOption(fake_total_sample_depth, "--fake_total_sample_depth", "to take into account the proportions, this is used a total detph per loci");
  setUp.setOption(num_threads, "--num_threads", "number of threads");

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
	auto simulated_populations = njh::json::parseFile(simulated_population_fnp.string());
  VecStr population_names = simulated_populations.getMemberNames();
  njh::concurrent::LockableQueue<std::string> pop_queue(population_names);

  std::function<void()> process_population = [&pop_queue,&input_vcf,
    index_start, fake_total_sample_depth, reset_gt_fields_ploidy,
    &setUp,
    &simulated_populations]() {
    std::string pop_name;
    while (pop_queue.getVal(pop_name)) {
      const auto & pop = simulated_populations[pop_name];
		  {
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
		    uint32_t sample_index = index_start;
		    for (const auto & sim_sample : pop["simulated_samples"]) {
			    auto sample_name = njh::leftPadNumStr(sample_index, pop["simulated_samples"].size());
			    ++sample_index;
			    uint32_t genotype_index = index_start;
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

		    OutputStream pop_out_vcf(njh::files::make_path(setUp.pars_.directoryName_, njh::pasteAsStr(pop_name, "_by_strain.vcf.gz")));
		    strains_separately.allAutoAdd_AN_AC_AF_InfoFields();
		    strains_separately.writeOutFixedAndSampleMeta(pop_out_vcf);
		  }
		  {
		    VCFOutput by_sample = input_vcf;
		    //clean up vcf for output
		    //clear samples
		    by_sample.samples_.clear();
		    for (auto & rec : by_sample.records_) {
			    rec.sampleFormatInfos_.clear();
		    }
		    //INFO/AC,INFO/AF,INFO/AN,^FORMAT/GT,FORMAT/AD
		    by_sample.otherHeaderMetaFields_.clear();
		    //format entries to remove
		    VecStr format_to_remove;
		    for (const auto & format : by_sample.formatEntries_) {
			    if (!njh::in(format.first, {"GT", "AD"})) {
				    format_to_remove.emplace_back(format.first);
			    }
		    }
		    for (const auto & format : format_to_remove) {
			    by_sample.formatEntries_.erase(format);
		    }
		    //info entries to remove
		    VecStr info_to_remove;
		    for (const auto & info : by_sample.infoEntries_) {
			    if (!njh::in(info.first, {"AC", "AF", "AN"})) {
				    info_to_remove.emplace_back(info.first);
			    }
		    }
		    for (const auto & info : info_to_remove) {
			    by_sample.infoEntries_.erase(info);
		    }

		    std::unordered_map<uint32_t, std::string> ancestral_population_index;
		    for (const auto & ancestral_pop_index : pop["ancestral_indexes"]) {
			    ancestral_population_index[ancestral_pop_index["index"].asUInt() ] = ancestral_pop_index["ancestral_genotype"].asString();
		    }
		    uint32_t sample_index = index_start;
		    for (const auto & sim_sample : pop["simulated_samples"]) {
			    auto sample_name = njh::leftPadNumStr(sample_index, pop["simulated_samples"].size());
			    ++sample_index;
			    uint32_t genotype_index = 0;
		      by_sample.samples_.emplace_back(sample_name);
		      std::vector<double> proportions( sim_sample["genotypes"].size(), 1.0/sim_sample["genotypes"].size());
		      if (sim_sample.isMember("prop")) {
		        proportions = njh::json::jsonArrayToVec<double>(sim_sample["prop"], [](const Json::Value & val){ return val.asDouble();});
		      }
		      std::vector<uint32_t> depth_per_genotypes;
          depth_per_genotypes.reserve(proportions.size());
          for (const auto &prop: proportions) {
            depth_per_genotypes.emplace_back(std::round(prop * fake_total_sample_depth));
          }
			    for (const auto & genotype : sim_sample["genotypes"]) {

				    for (const auto & segments : genotype["segments"]) {
					    Bed3RecordCore current_segment(segments["chrom"].asString(), segments["start"].asUInt(), segments["end"].asUInt());
					    auto ancestral_genotype = ancestral_population_index[segments["index"].asUInt()];
					    //assumes input is sorted
					    for (const auto & recPos : iter::enumerate(by_sample.records_)) {
						    auto zero_based_pos = recPos.element.pos_ - 1;
						    if (recPos.element.chrom_ > current_segment.chrom_ || (recPos.element.chrom_ == current_segment.chrom_ && zero_based_pos >= current_segment.chromEnd_)) {
							    break;
						    }
						    if (recPos.element.chrom_ == current_segment.chrom_ && zero_based_pos >= current_segment.chromStart_ && zero_based_pos < current_segment.chromEnd_) {

						      auto input_AD_toks = vecStrToVecNum<uint32_t>(tokenizeString(input_vcf.records_[recPos.index].sampleFormatInfos_[ancestral_genotype].getMeta("AD"), ","));
						      auto AD_sum = vectorSum(input_AD_toks);
						      std::vector<uint32_t> output_AD_toks;
						      if (0 == AD_sum) {
						        output_AD_toks = input_AD_toks;
						      } else {
						        for (const auto  AD : input_AD_toks) {
						          output_AD_toks.emplace_back(std::round(AD/AD_sum * depth_per_genotypes[genotype_index]));
						        }
						      }
                  if (0 == genotype_index) {
                    //first genotype so add to the record since none should exist yet for this sample
                    MetaDataInName output_sample_info;
                    output_sample_info.addMeta("AD", njh::conToStr(output_AD_toks, ","));
                    by_sample.records_[recPos.index].sampleFormatInfos_.emplace(sample_name, output_sample_info);
                  } else {
                    //not first genotype so add to the already present genotype
                    auto current_AD_toks = vecStrToVecNum<uint32_t>(tokenizeString(by_sample.records_[recPos.index].sampleFormatInfos_[sample_name].getMeta("AD"), ","));
                    for (const auto & AD_enum : iter::enumerate(output_AD_toks)) {
                      if (AD_enum.index >= current_AD_toks.size()) {
                        std::stringstream ss;
                        ss << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << ", error for sample: "
                        << sample_name << "AD_enum.index " << AD_enum.index
                        << " out of range of current ADs, " << current_AD_toks.size() << "\n";
                        throw std::runtime_error{ss.str()};
                      }
                      current_AD_toks[AD_enum.index]+= AD_enum.element;
                    }
                    by_sample.records_[recPos.index].sampleFormatInfos_[sample_name].addMeta("AD", njh::conToStr(current_AD_toks, ","), true);
                  }
						    }
					    }
				    }
			      ++genotype_index;
			    }
		    }
		    by_sample.allAddGTFields(reset_gt_fields_ploidy);
		    OutputStream pop_out_vcf(njh::files::make_path(setUp.pars_.directoryName_, njh::pasteAsStr(pop_name, "_by_sample.vcf.gz")));
		    by_sample.allAutoAdd_AN_AC_AF_InfoFields();
		    by_sample.writeOutFixedAndSampleMeta(pop_out_vcf);
		  }
    }
  };

  njh::concurrent::runVoidFunctionThreaded(process_population, num_threads);

	if (setUp.pars_.debug_){
		OutputStream out(njh::files::make_path(setUp.pars_.directoryName_, "out.vcf.gz"));
		input_vcf.writeOutFixedAndSampleMeta(out);
	}

	return 0;
}
} //namespace njhseq

