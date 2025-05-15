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
  bfs::path allele_table_fnp;
  std::string sample_id_col = "experiment_sample_name";
  std::string target_id_col = "target_name";
  std::string relative_abundance_col = "reads";
  std::string identifier_col = "asv";
  VecStr sample_order;
  double freq_cut_off = 0.70;
  double second_pass_prune_min_abundance = 0.10;
  double second_pass_freq_cut_off = 0.70;
  bfs::path regions_fnp;

  ampliconAnalysisSetUp setUp(inputCommands);
  setUp.processVerbose();
  setUp.processDebug();
  setUp.setOption(second_pass_prune_min_abundance, "--second_pass_prune_min_abundance", "second_pass_prune_min_abundance");
  setUp.setOption(second_pass_freq_cut_off, "--second_pass_freq_cut_off", "second_pass_freq_cut_off");


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
  auto regions = getBeds(regions_fnp);
  BedUtility::coordSort(regions);

  table allele_table(allele_table_fnp, "\t", true);
  allele_table.checkForColumnsThrow({sample_id_col, target_id_col, relative_abundance_col, identifier_col}, __PRETTY_FUNCTION__);

  //check for target and sample names
  VecStr missing_targets;
  auto targets_in_allele_table = allele_table.getColumnLevels(target_id_col);
  for (const auto & tar : regions) {
    if (njh::notIn(tar->name_, targets_in_allele_table)) {
      missing_targets.emplace_back(tar->name_);
    }
  }
  VecStr missing_samples;
  auto samples_in_allele_table = allele_table.getColumnLevels(sample_id_col);
  for (const auto & sample : sample_order) {
    if (njh::notIn(sample, samples_in_allele_table)) {
      missing_samples.emplace_back(sample);
    }
  }

  if (!missing_samples.empty() || !missing_targets.empty()) {
    std::stringstream ss;
    ss << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << ", error " << " " << "\n";
    if (!missing_samples.empty()) {
      ss << "missing the following sample from allele_table: " << njh::conToStr(missing_samples, ",") << "\n";
    }
    if (!missing_targets.empty()) {
      ss << "missing the following targets from allele_table: " << njh::conToStr(missing_targets, ",") << "\n";
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

  // second pass
  {
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

} // namespace njhseq



