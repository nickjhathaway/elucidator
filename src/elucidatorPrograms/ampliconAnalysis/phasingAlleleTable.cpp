//
// Created by Nicholas Hathaway on 5/6/25.
//


#include <njhseq/objects/BioDataObject/reading.hpp>
#include <njhseq/objects/BioDataObject/BioRecordsUtils/BedUtility.hpp>

#include "ampliconAnalysisRunner.hpp"

namespace njhseq {

template<typename T>
static void genRegCoordSort(std::vector<T> &genomic_regions, bool decending = false) {
  auto bedCoordSorterFunc =
      [](const T &reg1In, const T &reg2In) {
    const auto &reg1 = getRef(reg1In);
    const auto &reg2 = getRef(reg2In);
    if (reg1.chrom_ == reg2.chrom_) {
      if (reg1.start_ == reg2.start_) {
        return reg1.end_ < reg2.end_;
      } else {
        return reg1.start_ < reg2.start_;
      }
    } else {
      return reg1.chrom_ < reg2.chrom_;
    }
  };
  if (decending) {
    std::sort(genomic_regions.rbegin(), genomic_regions.rend(), bedCoordSorterFunc);
  } else {
    njh::sort(genomic_regions, bedCoordSorterFunc);
  }
}


class AllelePhaser {
public:



  class Microhaplotype {
  public:
    Microhaplotype(std::string name, double reads): name_(std::move(name)), reads_(reads) {

    }
    Microhaplotype() = default;
    std::string name_;
    double reads_{0};
    double relative_abundance_{0};

    bool operator==(const Microhaplotype& other) const {
      return name_ == other.name_;
    }

    bool operator==(const std::string & name) const {
      return name_ == name;
    }

    Json::Value toJson() const {
      Json::Value ret;
      ret["class"] = njh::getTypeName(*this);
      ret["name_"] = njh::json::toJson(name_);
      ret["reads_"] =  njh::json::toJson(reads_);
      ret["relative_abundance_"] =  njh::json::toJson(relative_abundance_);
      return ret;
    }
  };

  class Target {

    public:
    Target(std::string name): name_(std::move(name)) {

    }

    Target(std::string name, const std::vector<Microhaplotype> & haps): name_(std::move(name)), haps_(haps) {

    }
    std::string name_;
    std::vector<Microhaplotype> haps_;

    void addHap(Microhaplotype hap) {
      haps_.emplace_back(std::move(hap));
    }

    void set_relative_abundance() {
      double total_reads = 0;
      for (const auto& hap: haps_) {
        total_reads += hap.reads_;
      }
      for (auto & hap : haps_) {
        hap.relative_abundance_ = hap.reads_ / total_reads;
      }
    }

    double get_max_relative_abundance() const {
      double max = 0;
      for (const auto & hap: haps_) {
        if (max < hap.relative_abundance_) {
          max = hap.relative_abundance_;
        }
      }
      return max;
    }

    Microhaplotype get_major_haplotype() const {
      double max = 0;
      Microhaplotype ret;
      for (const auto & hap: haps_) {
        if (max < hap.relative_abundance_) {
          max = hap.relative_abundance_;
          ret = hap;
        }
      }
      return ret;
    }

    void prune_haplotypes(double min_abundance = 0.10) {
      std::vector<uint32_t> to_be_removed;
      for (const auto & hap : iter::enumerate(haps_)) {
        if (hap.element.relative_abundance_ < min_abundance) {
          to_be_removed.emplace_back(hap.index);
        }
      }
      //remove positions backwards to preserve location
      for (const auto & pos : iter::reversed(to_be_removed)) {
        haps_.erase(haps_.begin() + pos);
      }
    }
  };

  class PhasedHaplotype {
    public:
    PhasedHaplotype() = default;

    std::string hap_id_;
    VecStr target_names_;
    std::vector<Microhaplotype> haps_;


    double get_min_relative_abundance() const {
      double min = std::numeric_limits<double>::max();
      for (const auto & hap: haps_) {
        if (hap.relative_abundance_ < min) {
          min = hap.relative_abundance_;
        }
      }
      return min;
    }

    std::string gen_hap_id() const {
      std::string id;
      for (const auto & hap: haps_) {
        id += hap.name_;
      }
      return id;
    }

    Json::Value toJson() const {
      Json::Value ret;
      ret["class"] = njh::getTypeName(*this);
      ret["hap_id_"] = njh::json::toJson(hap_id_);
      ret["target_names_"] = njh::json::toJson(target_names_);
      ret["haps_"] = njh::json::toJson(haps_);
      return ret;
    }
  };

  class Sample {
    public:
    Sample(std::string name): name_(std::move(name)) {

    }
    std::string name_;
    std::unordered_map<std::string, Target> targets_;

    void set_relative_abundance() {
      for (auto & target: targets_) {
        target.second.set_relative_abundance();
      }
    }

    std::vector<PhasedHaplotype> get_possible_major_hap(const VecStr & target_names, double freq_cut_off = 0.70) const {
      std::vector<PhasedHaplotype> ret;
      bool major_possible = true;
      for (const auto & target : target_names) {
        if (njh::mapAt(targets_, target).get_max_relative_abundance() < freq_cut_off) {
          major_possible = false;
          break;
        }
      }
      if (major_possible) {
        PhasedHaplotype hap;
        for (const auto & target : target_names) {
          hap.target_names_.emplace_back(target);
          hap.haps_.emplace_back(njh::mapAt(targets_, target).get_major_haplotype());
        }
        ret.emplace_back(hap);
      }
      return ret;
    }

    void prune_haplotypes(double min_abundance = 0.10) {
      for ( auto & target : targets_) {
        target.second.prune_haplotypes(min_abundance);
      }
    }
  };

  AllelePhaser(const VecStr &sample_order,
               const std::vector<GenomicRegion> &target_regions): sample_order_(sample_order),
                                                                  target_regions_(target_regions) {
    genRegCoordSort(target_regions_);
    std::unordered_map<std::string, uint32_t> region_name_counts;
    for (const auto & region : target_regions_) {
      ++region_name_counts[region.uid_];
    }
    VecStr non_unique_region_names;
    for (const auto & name_count : region_name_counts) {
      if (name_count.second > 1) {
        non_unique_region_names.emplace_back(name_count.first);
      }
    }
    if (!non_unique_region_names.empty()) {
      std::stringstream ss;
      ss << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << ", error " << " region names must be unique, found the following duplicate region name " << "\n";
      ss << njh::conToStr(non_unique_region_names, ",") << "\n";
      throw std::runtime_error{ss.str()};
    }
  }

  void fill_with_target_placeholders() {
    for ( auto & sample: samples_) {
      for (const auto & tar : target_regions_) {
        if (njh::notIn(tar.uid_, sample.second.targets_)) {
          sample.second.targets_.emplace(tar.uid_, Target(tar.uid_, std::vector<Microhaplotype>{{"placeholder", 100}}));
        }
      }
    }
  }

  VecStr sample_order_;
  std::vector<GenomicRegion> target_regions_;
  std::unordered_map<std::string, Sample> samples_;

  std::unordered_map<std::string, uint32_t> gen_region_name_to_index() {
    std::unordered_map<std::string, uint32_t> ret;
    for (const auto & target : iter::enumerate(target_regions_)) {
      ret[target.element.uid_] = target.index;
    }
    return ret;
  }

  void set_relative_abundance() {
    for (auto & sample : samples_) {
      sample.second.set_relative_abundance();
    }
  }

  void prune_haplotypes(double min_abundance = 0.10) {
    for (auto & sample : samples_) {
      sample.second.prune_haplotypes(min_abundance);
    }
  }
  std::vector<PhasedHaplotype> get_all_possible_major_hap(double freq_cut_off = 0.70) {
    VecStr target_names;
    for (const auto & target : target_regions_) {
      target_names.emplace_back(target.uid_);
    }
    std::vector<PhasedHaplotype> ret;
    for (const auto & sample : sample_order_) {
      addOtherVec(ret, njh::mapAt(samples_,sample ).get_possible_major_hap(target_names, freq_cut_off));
    }
    return ret;
  }

};


int ampliconAnalysisRunner::phasingAlleleTable(
        const njh::progutils::CmdArgs & inputCommands) {
  bool skip_missing_allele_table_regions = false;
  bfs::path allele_table_fnp;
  std::string sample_id_col = "library_sample_name";
  std::string target_id_col = "target_name";
  std::string relative_abundance_col = "reads";
  std::string identifier_col = "seq";
  VecStr sample_order;
  double freq_cut_off = 0.70;
  double second_pass_prune_min_abundance = 0.10;
  double second_pass_freq_cut_off = 0.70;
  bfs::path regions_fnp;
  bool no_second_pass = false;
  ampliconAnalysisSetUp setUp(inputCommands);
  setUp.processVerbose();
  setUp.processDebug();
  setUp.setOption(second_pass_prune_min_abundance, "--second_pass_prune_min_abundance", "second_pass_prune_min_abundance");
  setUp.setOption(second_pass_freq_cut_off, "--second_pass_freq_cut_off", "second_pass_freq_cut_off");
  setUp.setOption(skip_missing_allele_table_regions, "--skip_missing_allele_table_regions", "skip missing allele table regions");

  setUp.setOption(no_second_pass, "--no_second_pass", "no_second_pass");


  setUp.setOption(freq_cut_off, "--freq_cut_off", "freq_cut_off");

  setUp.setOption(sample_id_col, "--sample_id_col", "sample_id_col");
  setUp.setOption(target_id_col, "--target_id_col", "target_id_col");
  setUp.setOption(relative_abundance_col, "--relative_abundance_col", "relative_abundance_col");
  setUp.setOption(identifier_col, "--identifier_col", "identifier_col");
  setUp.setOption(allele_table_fnp, "--allele_table_fnp", "allele_table_fnp", true);
  setUp.setOption(sample_order, "--sample_order", "sample_order", true);
  setUp.setOption(regions_fnp, "--regions_fnp", "regions_fnp", true);
  setUp.processDirectoryOutputName(njh::files::removeExtension(basename(allele_table_fnp)) + "_TODAY", true);
  setUp.finishSetUp(std::cout);
  setUp.startARunLog(setUp.pars_.directoryName_ );


  auto input_regions = getBeds(regions_fnp);
  if (input_regions.empty()) {
    std::stringstream ss;
    ss << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << ", error " << " no regions read in from: " <<  regions_fnp << "\n";
    throw std::runtime_error{ss.str()};
  }
  BedUtility::coordSort(input_regions);

  table allele_table(allele_table_fnp, "\t", true);
  allele_table.checkForColumnsThrow({sample_id_col, target_id_col, relative_abundance_col, identifier_col}, __PRETTY_FUNCTION__);

  //check for target and sample names
  VecStr targets_in_regions{};
  for (const auto & reg : input_regions) {
    targets_in_regions.emplace_back(reg->name_);
  }

  auto targets_in_allele_table = allele_table.getColumnLevels(target_id_col);
  njh::sort(targets_in_allele_table);
  njh::sort(targets_in_regions);
  auto targets_decomp = njh::decompose_sets_container(targets_in_regions, targets_in_allele_table);

  std::vector<std::shared_ptr<Bed6RecordCore>> regions;
  if (skip_missing_allele_table_regions) {
    for (const auto & reg : input_regions) {
      if (njh::in(reg->name_, targets_decomp.shared)) {
        regions.emplace_back(reg);
      }
    }
    targets_decomp.only_in_first.clear();
  } else {
    regions = input_regions;
  }

  auto samples_in_allele_table = allele_table.getColumnLevels(sample_id_col);
  njh::sort(sample_order);
  njh::sort(samples_in_allele_table);
  auto samples_decomp = njh::decompose_sets_container(sample_order, samples_in_allele_table);


  if (!samples_decomp.only_in_first.empty() || !samples_decomp.only_in_second.empty() || !targets_decomp.only_in_first.empty() || !targets_decomp.only_in_second.empty()) {
    std::stringstream ss;
    ss << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << ", error " << " " << "\n";
    if (!samples_decomp.only_in_first.empty()) {
      ss << "missing the following sample from allele_table: " << njh::conToStr(samples_decomp.only_in_first, ",") << "\n";
    }
    if (!samples_decomp.only_in_second.empty()) {
      ss << "missing the following sample from sample_order: " << njh::conToStr(samples_decomp.only_in_second, ",") << "\n";
    }
    if (!targets_decomp.only_in_first.empty()) {
      ss << "missing the following targets from allele_table: " << njh::conToStr(targets_decomp.only_in_first, ",") << "\n";
    }
    if (!targets_decomp.only_in_second.empty()) {
      ss << "missing the following targets from regions_fnp: " << njh::conToStr(targets_decomp.only_in_second, ",") << "\n";
    }
    throw std::runtime_error{ss.str()};
  }

  OutputStream phased_haplotypes_in_samples_out(njh::files::make_path(setUp.pars_.directoryName_, "phased_haplotypes_in_samples.tsv.gz"));
  OutputStream phased_haplotypes_out(njh::files::make_path(setUp.pars_.directoryName_, "phased_haplotypes.tsv.gz"));

  auto target_regions = bedPtrsToGenomicRegs(regions);

  AllelePhaser phaser(sample_order, target_regions);

  for (const auto & row : allele_table) {
    auto sample_name = row[allele_table.getColPos(sample_id_col)];
    if (njh::notIn(sample_name, phaser.samples_)) {
      phaser.samples_.emplace(sample_name, AllelePhaser::Sample(sample_name));
    }
    auto target_name = row[allele_table.getColPos(target_id_col)];
    if (njh::notIn(target_name, njh::mapAt(phaser.samples_,sample_name).targets_)) {
      njh::mapAt(phaser.samples_, sample_name).targets_.emplace(target_name, AllelePhaser::Target(target_name));
    }
    auto id = row[allele_table.getColPos(identifier_col)];
    auto rel_abund = njh::StrToNumConverter::stoToNum<double>(row[allele_table.getColPos(relative_abundance_col)]);
    bool found = false;
    for (auto & hap : njh::mapAt(njh::mapAt(phaser.samples_,sample_name).targets_, target_name).haps_) {
      if (id == hap.name_) {
        found = true;
        hap.reads_ += rel_abund;
        break;
      }
    }
    if (!found) {
      njh::mapAt(njh::mapAt(phaser.samples_,sample_name).targets_, target_name).haps_.emplace_back(id, rel_abund);
    }
  }

  phaser.fill_with_target_placeholders();
  phaser.set_relative_abundance();

  auto major_haplotypes = phaser.get_all_possible_major_hap(freq_cut_off);
  std::vector<AllelePhaser::PhasedHaplotype> unique_major_haplotypes;
  std::unordered_map<std::string, uint32_t> major_haplotypes_counts;
  for (const auto & hap : major_haplotypes) {
    auto hap_id = hap.gen_hap_id();
    if (njh::notIn(hap_id, major_haplotypes_counts)) {
      unique_major_haplotypes.emplace_back(hap);
    }
    ++major_haplotypes_counts[hap_id];
  }
  if (setUp.pars_.debug_) {
    std::cout << "unique_major_haplotypes.size(): " << unique_major_haplotypes.size() << std::endl;
    std::cout << "major_haplotypes.size(): " << major_haplotypes.size() << std::endl;
  }
  // second pass
  if (!no_second_pass){
    auto second_past_phaser = phaser;
    for (const auto & hap : iter::enumerate(unique_major_haplotypes)) {
      for (const auto & sample : second_past_phaser.sample_order_) {
        bool detected = true;
        AllelePhaser::PhasedHaplotype detected_hap;
        for (const auto & target : iter::enumerate(hap.element.target_names_)) {
          bool found = false;
          for (const auto & within_sample_hap : second_past_phaser.samples_.at(sample).targets_.at(target.element).haps_) {
            if (within_sample_hap.name_ == hap.element.haps_[target.index].name_) {
              detected_hap.haps_.emplace_back(within_sample_hap);
              detected_hap.target_names_.emplace_back(target.element);
              found = true;
              break;
            }
          }
          if (!found) {
            detected = false;
            break;
          }
        }
        if (detected) {
          auto min_rel_abundance = detected_hap.get_min_relative_abundance();
          for (const auto & detected_hap_target : iter::enumerate(detected_hap.target_names_)) {
            for (auto & within_sample_hap : second_past_phaser.samples_.at(sample).targets_.at(detected_hap_target.element).haps_) {
              if (within_sample_hap.name_ == hap.element.haps_[detected_hap_target.index].name_) {
                within_sample_hap.relative_abundance_ -= min_rel_abundance;
                break;
              }
            }
          }
        }
      }
    }
    second_past_phaser.prune_haplotypes(second_pass_prune_min_abundance);
    second_past_phaser.set_relative_abundance();
    if (setUp.pars_.debug_) {
      for (const auto & sample : njh::naturalSortNameRet(getVectorOfMapKeys(phaser.samples_))) {
        std::cout << "sample: " << sample << std::endl;
        for (const auto & tar : njh::naturalSortNameRet(getVectorOfMapKeys(phaser.samples_.at(sample).targets_)) ) {
          std::cout << "\t" << tar << " " << phaser.samples_.at(sample).targets_.at(tar).haps_.size() << std::endl;
          for (const auto & hap : phaser.samples_.at(sample).targets_.at(tar).haps_) {
            std::cout << "\t\t" << hap.name_ << " " << hap.relative_abundance_ << std::endl;
          }
        }
      }
      std::cout << std::endl;
      std::cout << std::endl;
      for (const auto & sample : njh::naturalSortNameRet(getVectorOfMapKeys(second_past_phaser.samples_))) {
        std::cout << "sample: " << sample << std::endl;
        for (const auto & tar : njh::naturalSortNameRet(getVectorOfMapKeys(second_past_phaser.samples_.at(sample).targets_)) ) {
          std::cout << "\t" << tar << " " << second_past_phaser.samples_.at(sample).targets_.at(tar).haps_.size() << std::endl;
          for (const auto & hap : second_past_phaser.samples_.at(sample).targets_.at(tar).haps_) {
            std::cout << "\t\t" << hap.name_ << " " << hap.relative_abundance_ << std::endl;
          }
        }
      }
    }
    auto second_pass_major_haplotypes = second_past_phaser.get_all_possible_major_hap(second_pass_freq_cut_off);
    for (const auto & hap : second_pass_major_haplotypes) {
      auto hap_id = hap.gen_hap_id();
      if (njh::notIn(hap_id, major_haplotypes_counts)) {
        unique_major_haplotypes.emplace_back(hap);
      }
      ++major_haplotypes_counts[hap_id];
    }
  }

  phased_haplotypes_out << "hap_id"
      << "\t" << target_id_col
      << "\t" << identifier_col << std::endl;
  for (const auto &hap: iter::enumerate(unique_major_haplotypes)) {
    for (const auto &target: iter::enumerate(hap.element.target_names_)) {
      phased_haplotypes_out << hap.index
      << "\t" << target.element
      << "\t" << hap.element.haps_[target.index].name_
      << std::endl;
    }
  }

  if (setUp.pars_.debug_) {
    std::cout << "haplotype\tcount" << std::endl;
    for (const auto & count : iter::enumerate(major_haplotypes_counts)) {
      std::cout << count.index << "\t" << count.element.second << std::endl;
    }
  }
  phased_haplotypes_in_samples_out << "hap_id"
  << "\t" << "hap_within_sample_freq"
  << "\t" << sample_id_col
  << "\t" << target_id_col
  << "\t" << identifier_col
  << "\t" << "reads"
  << "\t" << "within_sample_freq_for_target" << std::endl;

  for (const auto & hap : iter::enumerate(unique_major_haplotypes)) {
    for (const auto & sample : phaser.sample_order_) {
      bool detected = true;
      AllelePhaser::PhasedHaplotype detected_hap;
      for (const auto & target : iter::enumerate(hap.element.target_names_)) {
        bool found = false;
        for (const auto & within_sample_hap : phaser.samples_.at(sample).targets_.at(target.element).haps_) {
          if (within_sample_hap.name_ == hap.element.haps_[target.index].name_) {
            detected_hap.haps_.emplace_back(within_sample_hap);
            detected_hap.target_names_.emplace_back(target.element);
            found = true;
            break;
          }
        }

        if (!found) {
          if (setUp.pars_.debug_) {
            std::cout << "not found: " << std::endl;
            std::cout << "hap_id: " << hap.index << std::endl;
            std::cout << "target: " << target.element << std::endl;
            std::cout << "hap: " << hap.element.haps_[target.index].name_ << std::endl;
          }
          detected = false;
          break;
        }
      }
      if (detected) {
        auto min_rel_abundance = detected_hap.get_min_relative_abundance();
        for (const auto & detected_hap_target : iter::enumerate(detected_hap.target_names_)) {
          phased_haplotypes_in_samples_out << hap.index
              << "\t" << min_rel_abundance
              << "\t" << sample
              << "\t" << detected_hap_target.element
              << "\t" << detected_hap.haps_[detected_hap_target.index].name_
              << "\t" << detected_hap.haps_[detected_hap_target.index].reads_
              << "\t" << detected_hap.haps_[detected_hap_target.index].relative_abundance_
              << std::endl;
        }
      }
    }
  }


  return 0;
}

int ampliconAnalysisRunner::detectPhasedPartialAlleles(const njh::progutils::CmdArgs & inputCommands) {
  bfs::path allele_table_fnp;
  std::string sample_id_col = "library_sample_name";
  std::string target_id_col = "target_name";
  std::string relative_abundance_col = "reads";
  std::string identifier_col = "seq";
  bfs::path regions_fnp;
  uint32_t min_block_size = 2;
  uint64_t min_genomic_block_size = 0;
  bfs::path phased_haplotype_fnp;
  ampliconAnalysisSetUp setUp(inputCommands);
  setUp.processVerbose();
  setUp.processDebug();
  setUp.setOption(min_block_size, "--min_block_size", "min_block_size");
  setUp.setOption(min_genomic_block_size, "--min_genomic_block_size", "min_genomic_block_size");

  setUp.setOption(sample_id_col, "--sample_id_col", "sample_id_col");
  setUp.setOption(target_id_col, "--target_id_col", "target_id_col");
  setUp.setOption(relative_abundance_col, "--relative_abundance_col", "relative_abundance_col");
  setUp.setOption(identifier_col, "--identifier_col", "identifier_col");
  setUp.setOption(allele_table_fnp, "--allele_table_fnp", "allele_table_fnp", true);
  setUp.setOption(phased_haplotype_fnp, "--phased_haplotype_fnp", "previously phased haplotype fnp", true);
  setUp.setOption(regions_fnp, "--regions_fnp", "regions_fnp", true);


  setUp.processDirectoryOutputName(njh::files::removeExtension(basename(allele_table_fnp)) + "_TODAY", true);
  setUp.finishSetUp(std::cout);
  setUp.startARunLog(setUp.pars_.directoryName_ );

  auto regions = getBeds(regions_fnp);
  if (regions.empty()) {
    std::stringstream ss;
    ss << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << ", error " << " no regions read in from: " <<  regions_fnp << "\n";
    throw std::runtime_error{ss.str()};
  }
  BedUtility::coordSort(regions);

  table phased_haps_tab(phased_haplotype_fnp, "\t", true);
  phased_haps_tab.checkForColumnsThrow(VecStr{"hap_id", target_id_col, identifier_col}, __PRETTY_FUNCTION__);
  {
    //check for target in allele table
    VecStr missing_targets;
    auto targets_in_allele_table = phased_haps_tab.getColumnLevels(target_id_col);
    for (const auto & tar : regions) {
      if (njh::notIn(tar->name_, targets_in_allele_table)) {
        missing_targets.emplace_back(tar->name_);
      }
    }
    if (!missing_targets.empty()) {
      std::stringstream ss;
      ss << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << ", error " << " " << "\n";
      if (!missing_targets.empty()) {
        ss << "missing the following targets from phased_haps_tab: " << njh::conToStr(missing_targets, ",") << "\n";
      }
      throw std::runtime_error{ss.str()};
    }
  }
  table allele_table(allele_table_fnp, "\t", true);
  allele_table.checkForColumnsThrow({sample_id_col, target_id_col, relative_abundance_col, identifier_col}, __PRETTY_FUNCTION__);
  {
    //check for target in allele table
    VecStr missing_targets;
    auto targets_in_allele_table = allele_table.getColumnLevels(target_id_col);
    for (const auto & tar : regions) {
      if (njh::notIn(tar->name_, targets_in_allele_table)) {
        missing_targets.emplace_back(tar->name_);
      }
    }
    if (!missing_targets.empty()) {
      std::stringstream ss;
      ss << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << ", error " << " " << "\n";
      if (!missing_targets.empty()) {
        ss << "missing the following targets from allele_table: " << njh::conToStr(missing_targets, ",") << "\n";
      }
      throw std::runtime_error{ss.str()};
    }
  }
  auto target_regions = bedPtrsToGenomicRegs(regions);
  auto samples_in_allele_table = allele_table.getColumnLevels(sample_id_col);

  AllelePhaser phaser(samples_in_allele_table, target_regions);
  //std::cout << __FILE__ << " : " << __LINE__ << std::endl;
  for (const auto & row : allele_table) {
    auto sample_name = row[allele_table.getColPos(sample_id_col)];
    if (njh::notIn(sample_name, phaser.samples_)) {
      phaser.samples_.emplace(sample_name, AllelePhaser::Sample(sample_name));
    }
    auto target_name = row[allele_table.getColPos(target_id_col)];
    if (njh::notIn(target_name, njh::mapAt(phaser.samples_,sample_name).targets_)) {
      njh::mapAt(phaser.samples_, sample_name).targets_.emplace(target_name, AllelePhaser::Target(target_name));
    }
    auto id = row[allele_table.getColPos(identifier_col)];
    auto rel_abund = njh::StrToNumConverter::stoToNum<double>(row[allele_table.getColPos(relative_abundance_col)]);
    bool found = false;
    for (auto & hap : njh::mapAt(njh::mapAt(phaser.samples_,sample_name).targets_, target_name).haps_) {
      if (id == hap.name_) {
        found = true;
        hap.reads_ += rel_abund;
        break;
      }
    }
    if (!found) {
      njh::mapAt(njh::mapAt(phaser.samples_,sample_name).targets_, target_name).haps_.emplace_back(id, rel_abund);
    }
  }
  //std::cout << __FILE__ << " : " << __LINE__ << std::endl;
  phaser.fill_with_target_placeholders();
  phaser.set_relative_abundance();
  //std::cout << __FILE__ << " : " << __LINE__ << std::endl;
  // read in previous major haplotypes
  std::vector<AllelePhaser::PhasedHaplotype> unique_major_haplotypes;
  //std::cout << __FILE__ << " : " << __LINE__ << std::endl;
  {
    auto split_on_hap_id = phased_haps_tab.splitTableOnColumn("hap_id");
    for (const auto & hap_id_tab : split_on_hap_id) {
      //check for target in allele table
      VecStr missing_targets;
      auto targets_in_allele_table = hap_id_tab.second.getColumnLevels(target_id_col);
      for (const auto & tar : regions) {
        if (njh::notIn(tar->name_, targets_in_allele_table)) {
          missing_targets.emplace_back(tar->name_);
        }
      }
      if (!missing_targets.empty()) {
        std::stringstream ss;
        ss << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << ", error " << " " << "\n";
        if (!missing_targets.empty()) {
          ss << "missing the following targets for hap_id " << hap_id_tab.first << " : " << njh::conToStr(missing_targets, ",") << "\n";
        }
        throw std::runtime_error{ss.str()};
      }
      AllelePhaser::PhasedHaplotype phased_hap;
      phased_hap.hap_id_ = hap_id_tab.first;
      std::unordered_map<std::string, uint32_t> hap_target_index;
      for (const auto & target_name : iter::enumerate(hap_id_tab.second.getColumn(target_id_col))) {
        if (njh::in(target_name.element, hap_target_index)) {
          std::stringstream ss;
          ss << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << ", error " << " already have a haplotype for target: " << target_name.first << " for hap_id: " << hap_id_tab.first << "\n";
          throw std::runtime_error{ss.str()};
        }
        hap_target_index[target_name.element] = target_name.index;
      }
      auto hap_names = hap_id_tab.second.getColumn(identifier_col);

      for (const auto & target : phaser.target_regions_) {
        phased_hap.target_names_.emplace_back(target.uid_);
        phased_hap.haps_.emplace_back(hap_names[hap_target_index[target.uid_]], 1);
      }
      unique_major_haplotypes.emplace_back(phased_hap);
    }
  }
  //std::cout << __FILE__ << " : " << __LINE__ << std::endl;
  OutputStream detected_haps_out(njh::files::make_path(setUp.pars_.directoryName_, "detected_haplotypes_within_sample.tsv.gz"));
  OutputStream all_partial_haps(njh::files::make_path(setUp.pars_.directoryName_, "all_partial_haps.tsv.gz"));
  OutputStream detected_haps_target_names_out(njh::files::make_path(setUp.pars_.directoryName_, "detected_haplotypes_within_sample_target_names.tsv.gz"));
  OutputStream all_partial_target_names_haps(njh::files::make_path(setUp.pars_.directoryName_, "all_partial_haps_target_names.tsv.gz"));


  //std::cout << __FILE__ << " : " << __LINE__ << std::endl;
  std::vector<AllelePhaser::PhasedHaplotype> all_unique_partial_haplotypes;
  std::vector<std::vector<uint32_t>> all_unique_partial_haplotypes_target_positions;
  std::vector<std::set<std::string>> samples_per_partial_haps;
  //std::cout << __FILE__ << " : " << __LINE__ << std::endl;
  all_partial_haps << "partial_hap_id"
      << "\t" << "first_region"
      << "\t" << "last_region"
      << "\t" << "block_size"
      << "\t" << "chrom"
      << "\t" << "chrom_start"
      << "\t" << "chrom_end"
      << "\t" << "genomic_size"
      << "\t" << "sample_count"
      << "\t" << "covering_haps_count"
      << std::endl;

  all_partial_target_names_haps << "partial_hap_id"
      << "\t" << target_id_col
      << "\t" << identifier_col
      << std::endl;

  detected_haps_out << sample_id_col
      << "\t" << "partial_hap_id"
      << "\t" << "first_region"
      << "\t" << "last_region"
      << "\t" << "min_within_sample_freq_for_target_for_partial"
      << "\t" << "block_size"
      << "\t" << "chrom"
      << "\t" << "chrom_start"
      << "\t" << "chrom_end"
      << "\t" << "genomic_size"
      << "\t" << "covering_haps_count"
      << std::endl;
  detected_haps_target_names_out << sample_id_col
      << "\t" << "partial_hap_id"
      << "\t" << target_id_col
      << "\t" << identifier_col
      << "\t" << "within_sample_freq_for_target"
      << "\t" << "min_within_sample_freq_for_target_for_partial"
      << std::endl;
  //std::cout << __FILE__ << " : " << __LINE__ << std::endl;
  //now find blocks within samples
  for (const auto & sample : phaser.samples_) {
    //phased haplotypes have been placed in the same order has the target_region in phaser
    std::vector<AllelePhaser::PhasedHaplotype> partial_haplotypes;
    std::vector<std::vector<uint32_t>> partial_haplotypes_target_positions;
    //std::cout << __FILE__ << " : " << __LINE__ << std::endl;
    for (const auto & hap : unique_major_haplotypes) {
      std::vector<uint32_t> growing_haplotype_pos;
      AllelePhaser::PhasedHaplotype partial_hap;
      partial_hap.hap_id_ = hap.hap_id_;
      //std::cout << __FILE__ << " : " << __LINE__ << std::endl;
      auto check_and_reset = [&growing_haplotype_pos,&partial_hap,&phaser,&hap,
        min_block_size, min_genomic_block_size,
        &partial_haplotypes_target_positions, &partial_haplotypes]() {
        //check if current matching block passes
        auto genomic_distance_span = phaser.target_regions_[growing_haplotype_pos.back()].end_ -  phaser.target_regions_[growing_haplotype_pos.front()].start_;
        if (growing_haplotype_pos.size() >= min_block_size && genomic_distance_span >= min_genomic_block_size) {
          partial_haplotypes_target_positions.emplace_back(growing_haplotype_pos);
          partial_haplotypes.emplace_back(partial_hap);
        }
        //reset
        growing_haplotype_pos.clear();
        partial_hap = AllelePhaser::PhasedHaplotype();
        partial_hap.hap_id_ = hap.hap_id_;
      };
      //std::cout << __FILE__ << " : " << __LINE__ << std::endl;
      for (const auto & target_pos : iter::range(phaser.target_regions_.size())) {
        //std::cout << __FILE__ << " : " << __LINE__ << std::endl;
        const auto &target_hap = hap.haps_.at(target_pos);
        //std::cout << __FILE__ << " : " << __LINE__ << std::endl;
        //linking only on the same chromosome
        if (!growing_haplotype_pos.empty() && phaser.target_regions_[growing_haplotype_pos.back()].chrom_ != phaser.target_regions_[target_pos].chrom_) {
          check_and_reset();
        }
        //std::cout << __FILE__ << " : " << __LINE__ << std::endl;
        auto search_hap = std::find(sample.second.targets_.at(phaser.target_regions_[target_pos].uid_).haps_.begin(),
                                    sample.second.targets_.at(phaser.target_regions_[target_pos].uid_).haps_.end(),
                                    target_hap);
        //std::cout << __FILE__ << " : " << __LINE__ << std::endl;
        if (search_hap != sample.second.targets_.at(phaser.target_regions_[target_pos].uid_).haps_.end()) {
          //std::cout << __FILE__ << " : " << __LINE__ << std::endl;
          growing_haplotype_pos.emplace_back(target_pos);
          partial_hap.haps_.emplace_back(*search_hap);
        } else {
          //std::cout << __FILE__ << " : " << __LINE__ << std::endl;
          check_and_reset();
        }
        //std::cout << __FILE__ << " : " << __LINE__ << std::endl;
      }
      if (!growing_haplotype_pos.empty()) {
        check_and_reset();
      }
    }
    //std::cout << __FILE__ << " : " << __LINE__ << std::endl;
    std::vector<AllelePhaser::PhasedHaplotype> unique_partial_haplotypes;
    std::vector<std::vector<uint32_t>> unique_partial_haplotypes_target_positions;
    for (const auto & partial_hap_e : iter::enumerate(partial_haplotypes)) {
      bool found_match = false;
      for (const auto & unique_partial_haplotypes_target_positions_e : iter::enumerate(unique_partial_haplotypes_target_positions)) {
        if (partial_hap_e.element.gen_hap_id() == unique_partial_haplotypes[unique_partial_haplotypes_target_positions_e
              .index].gen_hap_id()
            && unique_partial_haplotypes_target_positions_e.element == partial_haplotypes_target_positions[partial_hap_e
              .index]) {
          found_match = true;
          //found append hap_id to the hap_id
          unique_partial_haplotypes[unique_partial_haplotypes_target_positions_e.index].hap_id_.append("," + partial_hap_e.element.hap_id_);
          break;
        }
      }
      if (!found_match) {
        unique_partial_haplotypes.emplace_back(partial_hap_e.element);
        unique_partial_haplotypes_target_positions.emplace_back(partial_haplotypes_target_positions[partial_hap_e.index]);
      }
    }
    //std::cout << __FILE__ << " : " << __LINE__ << std::endl;
    if (setUp.pars_.verbose_) {
      std::cout << "for sample " << sample.first << " found unique_partial_haplotypes.size(): " << unique_partial_haplotypes.size() << std::endl;
    }
    for (const auto & unique_partial_haplotypes_e : iter::enumerate(unique_partial_haplotypes)) {
      //std::cout << __FILE__ << " : " << __LINE__ << std::endl;
      double min_relative_abundance = std::numeric_limits<double>::max();
      //std::cout << __FILE__ << " : " << __LINE__ << std::endl;
      for (const auto & hap_e : iter::enumerate(unique_partial_haplotypes_e.element.haps_) ) {
        if (hap_e.element.relative_abundance_ < min_relative_abundance) {
          min_relative_abundance = hap_e.element.relative_abundance_;
        }
      }
      //std::cout << __FILE__ << " : " << __LINE__ << std::endl;
      auto genomic_size = phaser.target_regions_[unique_partial_haplotypes_target_positions[unique_partial_haplotypes_e.index].back()].end_ - phaser.target_regions_[unique_partial_haplotypes_target_positions[unique_partial_haplotypes_e.index].front()].start_;
      detected_haps_out << sample.first
              <<  "\t" << unique_partial_haplotypes_e.element.hap_id_
              << "\t" << phaser.target_regions_[unique_partial_haplotypes_target_positions[unique_partial_haplotypes_e.index].front()].uid_
              << "\t" << phaser.target_regions_[unique_partial_haplotypes_target_positions[unique_partial_haplotypes_e.index].back()].uid_
              << "\t" << min_relative_abundance
              << "\t" << unique_partial_haplotypes_target_positions[unique_partial_haplotypes_e.index].size()
      << "\t" << phaser.target_regions_[unique_partial_haplotypes_target_positions[unique_partial_haplotypes_e.index].front()].chrom_
      << "\t" << phaser.target_regions_[unique_partial_haplotypes_target_positions[unique_partial_haplotypes_e.index].front()].start_
      << "\t" << phaser.target_regions_[unique_partial_haplotypes_target_positions[unique_partial_haplotypes_e.index].back()].end_
              << "\t" << genomic_size
      << "\t" << countOccurences(unique_partial_haplotypes_e.element.hap_id_, ",")

      << std::endl;
      for (const auto &hap_e: iter::enumerate(unique_partial_haplotypes_e.element.haps_)) {
        detected_haps_target_names_out << sample.first
            << "\t" << unique_partial_haplotypes_e.element.hap_id_
            << "\t" << phaser.target_regions_[unique_partial_haplotypes_target_positions[unique_partial_haplotypes_e.
              index][hap_e.index]].uid_
            << "\t" << hap_e.element.name_
            << "\t" << hap_e.element.relative_abundance_
            << "\t" << min_relative_abundance
            << std::endl;
      }
      //std::cout << __FILE__ << " : " << __LINE__ << std::endl;
      //add to all unique_haplotypes
      bool found_match = false;
      auto testing_hap_id = unique_partial_haplotypes_e.element.gen_hap_id();
      for (const auto & hap_e : iter::enumerate(all_unique_partial_haplotypes)) {
        if (hap_e.element.gen_hap_id() == testing_hap_id && all_unique_partial_haplotypes_target_positions[hap_e.index] == unique_partial_haplotypes_target_positions[unique_partial_haplotypes_e.index]) {
          found_match = true;
          samples_per_partial_haps[hap_e.index].emplace(sample.first);
          break;
        }
      }
      if (!found_match) {
        all_unique_partial_haplotypes.emplace_back(unique_partial_haplotypes_e.element);
        all_unique_partial_haplotypes_target_positions.emplace_back(unique_partial_haplotypes_target_positions[unique_partial_haplotypes_e.index]);
        samples_per_partial_haps.emplace_back(std::set{sample.first});
      }
    }
    //std::cout << __FILE__ << " : " << __LINE__ << std::endl;
  }
  //std::cout << __FILE__ << " : " << __LINE__ << std::endl;
  for (const auto &all_unique_partial_haplotypes_e: iter::enumerate(all_unique_partial_haplotypes)) {
    auto genomic_size = phaser.target_regions_[all_unique_partial_haplotypes_target_positions[
                          all_unique_partial_haplotypes_e.index].back()].end_ - phaser.target_regions_[
                          all_unique_partial_haplotypes_target_positions[
                            all_unique_partial_haplotypes_e.index].front()].start_;
    all_partial_haps
        << all_unique_partial_haplotypes_e.element.hap_id_
        << "\t" << phaser.target_regions_[
          all_unique_partial_haplotypes_target_positions[
            all_unique_partial_haplotypes_e.index].front()].uid_
        << "\t" << phaser.target_regions_[all_unique_partial_haplotypes_target_positions[
          all_unique_partial_haplotypes_e.index].back()].uid_
        << "\t" << all_unique_partial_haplotypes_target_positions[all_unique_partial_haplotypes_e.index].size()
        << "\t" << phaser.target_regions_[
          all_unique_partial_haplotypes_target_positions[
            all_unique_partial_haplotypes_e.index].front()].chrom_
        << "\t" << phaser.target_regions_[
          all_unique_partial_haplotypes_target_positions[
            all_unique_partial_haplotypes_e.index].front()].start_
        << "\t" << phaser.target_regions_[all_unique_partial_haplotypes_target_positions[
          all_unique_partial_haplotypes_e.index].back()].end_
        << "\t" << genomic_size
        << "\t" << samples_per_partial_haps[all_unique_partial_haplotypes_e.index].size()
    << "\t" << countOccurences(all_unique_partial_haplotypes_e.element.hap_id_, ",")
        << std::endl;
    for (const auto &hap_e: iter::enumerate(all_unique_partial_haplotypes_e.element.haps_)) {
      all_partial_target_names_haps
          << all_unique_partial_haplotypes_e.element.hap_id_
          << "\t" << phaser.target_regions_[all_unique_partial_haplotypes_target_positions[
            all_unique_partial_haplotypes_e.index][hap_e.index]].uid_
          << "\t" << hap_e.element.name_
          << std::endl;
    }
  }
  //std::cout << __FILE__ << " : " << __LINE__ << std::endl;
  return 0;
}



} // namespace njhseq



