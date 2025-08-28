//
// Created by Nicholas Hathaway on 8/20/25.
//

#include <njhseq/concurrency/PairwisePairFactory.hpp>
#include <njhseq/concurrency/pools/AlignerPool.hpp>

#include "phipseq_utils_runner.hpp"

#include <njhseq/IO/SeqIO/SeqIO.hpp>
#include <njhseq/objects/seqObjects/seqKmers/seqWithKmerInfo.hpp>


namespace njhseq {

template <typename F>
void generate_combinations(const std::vector<char>& alphabet,
                           std::size_t length,
                           F&& consume) {
  if (alphabet.empty()) return;

  // Edge-case: length 0 => single empty string
  if (length == 0) { consume(std::string()); return; }

  const std::size_t k = alphabet.size();
  std::string s(length, alphabet[0]);         // reused output buffer
  std::vector<std::size_t> idx(length, 0);    // digits in base k

  bool done = false;
  while (!done) {
    consume(s); // use the current combination

    // increment "odometer" (rightmost digit first)
    for (std::size_t pos = length; pos-- > 0; ) {
      if (++idx[pos] < k) {
        s[pos] = alphabet[idx[pos]];
        break; // no carry; continue generating
      } else {
        idx[pos] = 0;
        s[pos] = alphabet[0]; // carry to the next more-significant digit
        if (pos == 0) done = true; // overflowed the most-significant digit
      }
    }
  }
}

class PhipSeqNucLibraryGraph{
  public:

  struct DBSCAN_PARS {
    uint32_t min_neighbors_inclusive_{3};
    uint32_t epsilon_inclusive_{1};

  };

  class edge {
  public:
    edge(std::unordered_map<std::string, uint32_t> key, const uint32_t edit_dist) : name_to_other_node_pos_(std::move(key)),
                                                                              hamming_dist_(edit_dist) {
    }

    std::unordered_map<std::string, uint32_t> name_to_other_node_pos_;//<! index is current node, value is the other node position
    uint32_t hamming_dist_;
    bool on_{true};
  };
  class node {
    public:
    node(const std::shared_ptr<seqInfo> & seq): seqBase_(seq) {

    }
    std::shared_ptr<seqInfo> seqBase_;
    std::vector<std::shared_ptr<edge>> edges_;//<! other sequences below a specific edit distance away

    std::string barcode_;//<! a random molecular barcode
    uint64_t group_{std::numeric_limits<uint64_t>::max()};

    bool corePoint(const DBSCAN_PARS pars) const {
      uint32_t count = 0;
      for (const auto & e : edges_) {
        if (e->hamming_dist_ <= pars.epsilon_inclusive_) {
          ++count;
        }
      }
      return count <= pars.min_neighbors_inclusive_;
    }

    std::vector<std::shared_ptr<edge>> getEdgesLessEqualEps(const uint32_t & eps) const{
      std::vector<std::shared_ptr<edge>> ret;
      for (const auto & e : edges_) {
        if (e->hamming_dist_ <= eps) {
          ret.emplace_back(e);
        }
      }
      return ret;
    }


    std::vector<uint64_t> getNeighborNodesLessEqualEps(const uint32_t & eps) const {
      std::vector<uint64_t> ret;
      for (const auto & e : edges_) {
        if (e->hamming_dist_ <= eps) {
          ret.emplace_back(e->name_to_other_node_pos_[seqBase_->name_]);
        }
      }
      return ret;
    }
  };

  void reset_groups() {
    njh::for_each(nodes_, [](node & node) {node.group_ = std::numeric_limits<uint64_t>::max();});
  }

  std::vector<node> nodes_;
  std::vector<std::shared_ptr<edge>> all_edges_;

  void determine_groups_dbscan(const DBSCAN_PARS & pars) {
    // reset all groups to std::numeric_limits<uint64_t>::max()
    reset_groups();
    // now link groups
    uint64_t group = 0;
    constexpr uint64_t NOISE = std::numeric_limits<uint64_t>::max() - 1;

    for (const auto node_pos : iter::range(nodes_.size())) {
      //unassigned
      if (std::numeric_limits<uint64_t>::max() == nodes_[node_pos].group_) {
        if (nodes_[node_pos].corePoint(pars)) {
          nodes_[node_pos].group_ = group;
          std::deque<uint64_t> next_core_neighbors;
          {
            //only need to spread to unassigned points or noise points
            std::vector<std::shared_ptr<edge>> filtered_edges;
            for (const auto & e : nodes_[node_pos].getEdgesLessEqualEps(pars.epsilon_inclusive_)) {
              const auto other_node = e->name_to_other_node_pos_[nodes_[node_pos].seqBase_->name_];
              if (std::numeric_limits<uint64_t>::max() == nodes_[other_node].group_ || NOISE == nodes_[other_node].group_) {
                filtered_edges.emplace_back(e);
              }
            }

            for (const auto & e : filtered_edges) {
              nodes_[e->name_to_other_node_pos_[nodes_[node_pos].seqBase_->name_]].group_ = group;
              if (nodes_[e->name_to_other_node_pos_[nodes_[node_pos].seqBase_->name_]].corePoint(pars)) {
                next_core_neighbors.emplace_back(e->name_to_other_node_pos_[nodes_[node_pos].seqBase_->name_]);
              }
            }
          }
          while (!next_core_neighbors.empty()) {
            auto next_node = next_core_neighbors.front();
            next_core_neighbors.pop_front();
            //only need to spread to unassigned points or noise points
            std::vector<std::shared_ptr<edge>> next_filtered_edges;
            for (const auto & e : nodes_[next_node].getEdgesLessEqualEps(pars.epsilon_inclusive_)) {
              const auto other_node = e->name_to_other_node_pos_[nodes_[next_node].seqBase_->name_];
              if (std::numeric_limits<uint64_t>::max() == nodes_[other_node].group_ || NOISE == nodes_[other_node].group_) {
                next_filtered_edges.emplace_back(e);
              }
            }
            for (const auto & e : next_filtered_edges) {
              nodes_[e->name_to_other_node_pos_[nodes_[next_node].seqBase_->name_]].group_ = group;
              if (nodes_[e->name_to_other_node_pos_[nodes_[next_node].seqBase_->name_]].corePoint(pars)) {
                next_core_neighbors.emplace_back(e->name_to_other_node_pos_[nodes_[next_node].seqBase_->name_]);
              }
            }
          }
          ++group;
        } else {
          nodes_[node_pos].group_ = NOISE;
        }
      }
    }
  }

  void determine_groups_single_linkage() {
    // set all groups to std::numeric_limits<uint64_t>::max()
    reset_groups();
    // now link gorups
    uint64_t group = 0;
    for (const auto node_pos : iter::range(nodes_.size())) {
      if (std::numeric_limits<uint64_t>::max() == nodes_[node_pos].group_) {
        nodes_[node_pos].group_ = group;
        if (!nodes_[node_pos].edges_.empty()) {
          std::deque<uint64_t> other_nodes;
          for (const auto & e : nodes_[node_pos].edges_) {
            if (std::numeric_limits<uint64_t>::max() == nodes_[e->name_to_other_node_pos_[nodes_[node_pos].seqBase_->name_]].group_) {
              other_nodes.emplace_back(e->name_to_other_node_pos_[nodes_[node_pos].seqBase_->name_]);
            }
          }
          while (!other_nodes.empty()) {
            auto next_node = other_nodes.front();
            other_nodes.pop_front();
            nodes_[next_node].group_ = group;
            for (const auto & e : nodes_[next_node].edges_) {
              if (std::numeric_limits<uint64_t>::max() == nodes_[e->name_to_other_node_pos_[nodes_[next_node].seqBase_->name_]].group_) {
                other_nodes.emplace_back(e->name_to_other_node_pos_[nodes_[next_node].seqBase_->name_]);
              }
            }
          }
        }
        ++group;
      }
    }
  }

  struct generate_random_barcodes_pars {
    VecStr restriction_sites_to_remove{"GAATTC", "CTCGAG", "AAGCCT"};
    // ecor1 = "GAATTC"
    // xho1 = "CTCGAG"
    // hindIII = "AAGCCT"
    std::string prepend_to_barcode = "TGATAA";
    std::vector<char> barcode_alphabet{'A', 'G', 'C', 'T'};
    double barcode_entropy_filter = 0.52;
    uint32_t barcode_entropy_klen = 3;
    uint32_t barcode_size = 0;
    uint32_t barcode_hamming_distance_cut_off_exclusive = 3;
    uint32_t linker_overlapping_bases_cut_off_inclusive = 9;//!< will reject barcode if linkers overlap by at least this much and pass the error threshold
    uint32_t linker_mismatch_cut_off_inclusive = 2;//!< will reject barcode if linkers match barcodes with this many errors or less
    std::string reverse_linker_seq_3_5;
    std::string forward_linker_seq_5_3;

    VecStr all_restriction_sites_to_check;
    
    void set_all_restriction_sites_to_check() {
      VecStr restriction_sites_to_remove_rev_comp;
      for (const auto & site : restriction_sites_to_remove) {
        auto rev_comp_site = seqUtil::reverseComplement(site, "DNA");
        if (rev_comp_site != site) {
          restriction_sites_to_remove_rev_comp.emplace_back(rev_comp_site);
        }
      }
      all_restriction_sites_to_check = concatVecs(restriction_sites_to_remove, restriction_sites_to_remove_rev_comp);
    }

    void set_options(seqSetUp & setUp, bool barcode_size_required) {
      setUp.setOption(linker_overlapping_bases_cut_off_inclusive, "--linker_overlapping_bases_cut_off_inclusive", "linker_overlapping_bases_cut_off_inclusive");
      setUp.setOption(linker_mismatch_cut_off_inclusive, "--linker_mismatch_cut_off_inclusive", "linker_mismatch_cut_off_inclusive");
      setUp.setOption(restriction_sites_to_remove, "--restriction_sites_to_remove", "restriction_sites_to_remove");
      setUp.setOption(prepend_to_barcode, "--prepend_to_barcode", "sequence to prepend to barcode");
      setUp.setOption(barcode_entropy_klen, "--barcode_entropy_klen", "barcode entropy klen for calculating entropy");
      setUp.setOption(barcode_entropy_filter, "--barcode_entropy_filter", "barcode entropy filter");
      setUp.setOption(barcode_alphabet, "--barcode_alphabet", "alphabet for barcodes");
      setUp.setOption(reverse_linker_seq_3_5, "--reverse_linker_seq_3_5", "reverse linker seq to compare to ensure to have a random barcode that is not similar to this sequence");
      setUp.setOption(forward_linker_seq_5_3, "--forward_linker_seq_5_3", "forward linker seq to compare to ensure to have a random barcode that is not similar to this sequence");
      if (barcode_size_required) {
        setUp.setOption(barcode_size, "--barcode_size", "barcode size", true);
        setUp.setOption(barcode_hamming_distance_cut_off_exclusive, "--barcode_hamming_distance_cut_off_exclusive", "barcode hamming distance difference cut off for similar seqs", true);
        if (barcode_hamming_distance_cut_off_exclusive > barcode_size) {
          setUp.failed_ = true;
          setUp.addWarning(njh::pasteAsStr("--barcode_hamming_distance_cut_off_exclusive: ", barcode_hamming_distance_cut_off_exclusive,
            "can't be more than --barcode_size: ", barcode_size));
        }
      }
      set_all_restriction_sites_to_check();
    }
  };
  static VecStr generate_random_barcodes(const generate_random_barcodes_pars & pars) {
    //generating barcodes with entropy cut off to create non-low complexity barcodes
    std::vector<std::string> barcodes;
    aligner alignObj(2  * (pars.prepend_to_barcode.size() + std::max(pars.reverse_linker_seq_3_5.size(), pars.forward_linker_seq_5_3.size()) + pars.barcode_size),
      gapScoringParameters(5,1,0,0,0,0),
      substituteMatrix::createScoreMatrix(2, -2, false, true, true));


    generate_combinations(pars.barcode_alphabet, pars.barcode_size, [
      &pars,
      &alignObj,
      &barcodes
                          ](const std::string & barcode) {
      auto kinfo = kmerInfo(barcode, pars.barcode_entropy_klen, false);
      if (kinfo.computeKmerEntropy() >= pars.barcode_entropy_filter) {
        std::string pre_barcode_reverse_linker = pars.prepend_to_barcode + barcode + pars.reverse_linker_seq_3_5;
        bool fails_linker_check = false;
        if (!pars.forward_linker_seq_5_3.empty() || !pars.reverse_linker_seq_3_5.empty()) {
          seqInfo pre_barcode_reverse_linker_info("pre_barcode_reverse_linker", pre_barcode_reverse_linker);
          if (!pars.reverse_linker_seq_3_5.empty()) {
            seqInfo surrounding_sequence_info("prepend_to_barcode_barcode", pars.prepend_to_barcode + barcode);
            seqInfo linker_info("linker_seq", pars.reverse_linker_seq_3_5);
            alignObj.alignRegGlobalNoInternalGaps(surrounding_sequence_info, linker_info);
            alignObj.profilePrimerAlignment(surrounding_sequence_info, linker_info);
            if (alignObj.comp_.distances_.basesInAln_ >= pars.linker_overlapping_bases_cut_off_inclusive && alignObj.comp_.hqMismatches_ <= pars.linker_mismatch_cut_off_inclusive) {
              fails_linker_check = true;
            }
            if (!fails_linker_check) {
              //check the reverse complement as well
              linker_info.reverseComplementRead(true, true);
              alignObj.alignRegGlobalNoInternalGaps(pre_barcode_reverse_linker_info, linker_info);
              alignObj.profilePrimerAlignment(pre_barcode_reverse_linker_info, linker_info);
              if (alignObj.comp_.distances_.basesInAln_ >= pars.linker_overlapping_bases_cut_off_inclusive && alignObj.comp_.
                  hqMismatches_ <= pars.linker_mismatch_cut_off_inclusive) {
                fails_linker_check = true;
              }
            }
          }
          if (!fails_linker_check && !pars.forward_linker_seq_5_3.empty()) {
            //check the reverse complement as well
            seqInfo linker_info("forward_linker_seq_5_3", pars.forward_linker_seq_5_3);
            alignObj.alignRegGlobalNoInternalGaps(pre_barcode_reverse_linker_info, linker_info);
            alignObj.profilePrimerAlignment(pre_barcode_reverse_linker_info, linker_info);
            if (alignObj.comp_.distances_.basesInAln_ >= pars.linker_overlapping_bases_cut_off_inclusive && alignObj.comp_.
                hqMismatches_ <= pars.linker_mismatch_cut_off_inclusive) {
              fails_linker_check = true;
            }
            if (!fails_linker_check) {
              linker_info.reverseComplementRead(true, true);
              alignObj.alignRegGlobalNoInternalGaps(pre_barcode_reverse_linker_info, linker_info);
              alignObj.profilePrimerAlignment(pre_barcode_reverse_linker_info, linker_info);
              if (alignObj.comp_.distances_.basesInAln_ >= pars.linker_overlapping_bases_cut_off_inclusive && alignObj.comp_.
                  hqMismatches_ <= pars.linker_mismatch_cut_off_inclusive) {
                fails_linker_check = true;
              }
            }
          }
        }
        bool contains_restriction_sites = false;
        for (const auto & site : pars.all_restriction_sites_to_check) {
          if (pre_barcode_reverse_linker.find(site) != std::string::npos) {
            contains_restriction_sites = true;
            break;
          }
        }
        if (!contains_restriction_sites && !fails_linker_check) {
          barcodes.emplace_back(barcode);
        }
      }
    });
    return barcodes;
  }
  
};

/**
 * @brief Compute the Hamming distance between two equal-length strings (no length check).
 *
 * Counts the number of positions at which the corresponding characters are different.
 *
 * @param s1 First string view. Must be the same length as @p s2.
 * @param s2 Second string view. Must be the same length as @p s1.
 * @return Number of mismatching positions.
 *
 * @pre s1.size() == s2.size()
 * @note No bounds or size checks are performed.
 * @warning Passing strings of different lengths results in undefined behavior.
 * @remark Uses C++17 std::transform_reduce for a single-pass O(n) implementation.
 */
inline uint32_t hamming_distance_no_check(std::string_view s1,
                                          std::string_view s2) noexcept {
    return std::transform_reduce(
        s1.begin(), s1.end(), s2.begin(),
        uint32_t{0},
        std::plus<uint32_t>{},
        [](char a, char b) -> uint32_t { return static_cast<uint32_t>(a != b); }
    );
}

/**
 * @brief Compute the Hamming distance between two strings with a length check.
 *
 * Validates that the inputs have the same length, then computes the
 * number of positions at which the characters differ.
 *
 * @param s1 First string view.
 * @param s2 Second string view.
 * @return Number of mismatching positions.
 *
 * @throws std::invalid_argument if @p s1 and @p s2 have different lengths.
 * @remark Delegates to hamming_distance_no_check() after validation.
 */
inline uint32_t hamming_distance(std::string_view s1, std::string_view s2) {
    if (s1.size() != s2.size()) {
        std::ostringstream oss;
        oss << "hamming_distance: strings must be the same length. "
            << "s1 size: " << s1.size() << ", s2 size: " << s2.size();
        throw std::invalid_argument(oss.str());
    }
    return hamming_distance_no_check(s1, s2);
}


int PhipSeqUtilsRunner::markGroupsByHammingDistanceCutOff(const njh::progutils::CmdArgs & libraryCommands) {
  uint32_t hamming_distance_cut_off_exclusive = 4;
  uint32_t numThreads = 1;
  PhipSeqNucLibraryGraph::DBSCAN_PARS dbscanPars;
  bool do_dbscan = false;

  OutOptions outGroupCounts("", ".tsv");
  seqSetUp setUp(libraryCommands);
  setUp.description_ = "append a random barcode sequence made of AGTC to the end of a phipeq library";
  setUp.processVerbose();
  setUp.processDebug();
  // setUp.processDefaultReader(true);
  setUp.processDefaultReader(true);
  setUp.setOption(outGroupCounts.outFilename_, "--group_counts_out_fnp", "Optional file to write out group counts to");
  outGroupCounts.transferOverwriteOpts(setUp.pars_.ioOptions_.out_);
  setUp.setOption(hamming_distance_cut_off_exclusive, "--hamming_distance_cut_off_exclusive", "hamming_distance_cut_off_exclusive", true);
  setUp.setOption(do_dbscan, "--do_dbscan", "do dbscan groupings rather than single linkage clustering");
  setUp.setOption(dbscanPars.epsilon_inclusive_, "--dbscan_epsilon_inclusive", "dbscan epsilon inclusive");
  setUp.setOption(dbscanPars.min_neighbors_inclusive_, "--dbscan_min_neighbors_inclusive", "dbscan min neighbors inclusive");
  setUp.setOption(numThreads, "--numThreads", "number of threads to use");
  setUp.finishSetUp(std::cout);

  SeqIO seq_io(setUp.pars_.ioOptions_);
  seq_io.openIn();
  seq_io.openOut();
  auto library = seq_io.in_.readAllReadsPtrs<seqInfo>();

  std::unordered_map<std::string, uint32_t> name_counts;
  for (const auto & seq : library) {
    ++name_counts[seq->name_];
    if (len(*seq) != len(*library.front())) {
      std::stringstream ss;
      ss << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << ", error "
      << ", all library has to be the same length, seq: " << seq->name_
      << " has len: " << len(*seq) << " which is different from "
      << library.front()->name_ << " which has len " << len(*library.front()) << "\n";
      throw std::runtime_error{ss.str()};
    }
  }
  {
    VecStr multiple_names;
    for (const auto & name_count : name_counts) {
      if (name_count.second > 1) {
        multiple_names.emplace_back(name_count.first);
      }
    }
    if (!multiple_names.empty()) {
      std::stringstream ss;
      ss << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << ", error " << " the following names were found multiple times, all names must be unique" << "\n";
      ss << njh::conToStr(multiple_names, ",") << "\n";
      throw std::runtime_error{ss.str()};
    }
  }
  std::unique_ptr<OutputStream> group_counts_out;
  if (!outGroupCounts.outFilename_.empty()) {
    group_counts_out = std::make_unique<OutputStream>(outGroupCounts);
  }

  PairwisePairFactory pairFactory(library.size());
  //set up progress bar
  njh::ProgressBar pBar(pairFactory.totalCompares_);

  njhseq::concurrent::AlignerPool aligner_pool(
    len(*library.front()) *2, gapScoringParameters(10, 1), substituteMatrix::createScoreMatrix(1, 0, false, false, true), numThreads
  );
  aligner_pool.initAligners();

  PhipSeqNucLibraryGraph libraryGraph;
  for (const auto & seq : library) {
    libraryGraph.nodes_.emplace_back(seq);
  }
  std::unordered_map<uint32_t, uint32_t> counts;
  std::mutex counts_mut;
  std::function<void()> getEditDistances = [&counts,&counts_mut,&pairFactory,&library,&aligner_pool,&setUp, &pBar,
    &hamming_distance_cut_off_exclusive, &libraryGraph]() {
    PairwisePairFactory::PairwisePairVec pairs;
    auto current_aligner = aligner_pool.popAligner();
    std::vector<std::shared_ptr<PhipSeqNucLibraryGraph::edge>> current_edges;
    std::unordered_map<uint32_t, uint32_t> current_counts;
    while (pairFactory.setNextPairs(pairs, 1000)) {
      if (setUp.pars_.verbose_) {
        pBar.outputProgAdd(std::cout, pairs.pairs_.size(), true);
      }
      for (const auto & pair : pairs.pairs_) {
        current_aligner->noAlignSetAndScore(library[pair.col_], library[pair.row_]);
        auto edit_dist = library[pair.col_]->seq_.size() - current_aligner->parts_.score_;
        ++current_counts[edit_dist];
        if (edit_dist < hamming_distance_cut_off_exclusive) {
          //add to edges to be added later
          current_edges.emplace_back(std::make_shared<PhipSeqNucLibraryGraph::edge>(
            std::unordered_map<std::string, uint32_t>{
              {library[pair.col_]->name_, pair.row_},
              {library[pair.row_]->name_, pair.col_}
            }, edit_dist));
        }
      }
    }
    {
      std::lock_guard lock(counts_mut);
      for (const auto & count : current_counts) {
        counts[count.first] += count.second;
      }
      for ( auto & e : current_edges) {
        //add to nodes
        auto node_positions = getVectorOfMapValues(e->name_to_other_node_pos_);
        for (const auto node_pos : node_positions) {
          libraryGraph.nodes_[node_pos].edges_.emplace_back(e);
        }
        //add to edges
        libraryGraph.all_edges_.emplace_back(e);
      }
    }
  };

  njh::concurrent::runVoidFunctionThreaded(getEditDistances, numThreads);
  std::unordered_map<uint32_t, uint32_t> neighbors_counts;
  std::vector<uint64_t> nodes_with_neighbors;
  for (const auto n_pos : iter::range(libraryGraph.nodes_.size())) {
    const auto & n = libraryGraph.nodes_[n_pos];
    ++neighbors_counts[n.edges_.size()];
    if (n.edges_.size() > 0) {
      nodes_with_neighbors.emplace_back(n_pos);
    }
  }
  std::map<uint32_t, uint32_t> group_counts;
  if (do_dbscan) {
    libraryGraph.determine_groups_dbscan(dbscanPars);
  } else {
    libraryGraph.determine_groups_single_linkage();
  }

  for (const auto & node : libraryGraph.nodes_) {
    MetaDataInName meta;
    if (MetaDataInName::nameHasMetaData(node.seqBase_->name_)) {
      meta = MetaDataInName(node.seqBase_->name_);
    }
    ++group_counts[node.group_];
    meta.addMeta("hamming_distanace_group", node.group_, true);
    meta.resetMetaInName(node.seqBase_->name_);
    seq_io.write(node.seqBase_);
  }
  if (!outGroupCounts.outFilename_.empty()) {
    *group_counts_out << "group\tcount" << std::endl;
    for (const auto & group_count : group_counts) {
      *group_counts_out << group_count.first << "\t" << group_count.second << std::endl;
    }
  }
  if (setUp.pars_.debug_) {
    std::cerr << "neighbor counts:" << std::endl;
    std::cerr << "neighbors\tcount" << std::endl;
    auto counts_key = njh::getSetOfMapKeys(neighbors_counts);
    for (const auto & dist : counts_key) {
      std::cerr << dist << "\t" << neighbors_counts[dist] << std::endl;
    }
  }
  return 0;
}

int PhipSeqUtilsRunner::appendRandomBarcode(const njh::progutils::CmdArgs & libraryCommands) {
  PhipSeqNucLibraryGraph::generate_random_barcodes_pars random_barcodes_pars;
  uint32_t hamming_distance_cut_off_exclusive = 5;
  uint64_t shuffle_seed = std::numeric_limits<uint64_t>::max();
  uint32_t numThreads = 1;
  // OutOptions outOpts("", ".tsv");
  bool mark_with_groups = false;
  PhipSeqNucLibraryGraph::DBSCAN_PARS dbscanPars;
  bool do_dbscan = false;
  seqSetUp setUp(libraryCommands);
  setUp.description_ = "append a random barcode sequence made of AGTC to the end of a phipeq library";
  setUp.processVerbose();
  setUp.processDebug();
  random_barcodes_pars.set_options(setUp, true);
  setUp.setOption(hamming_distance_cut_off_exclusive, "--hamming_distance_cut_off_exclusive", "below this hamming distance cut off exclusive, ensure barcods are signifcnatly different");
  setUp.setOption(mark_with_groups, "--mark_with_groups", "mark with groups");
  setUp.setOption(do_dbscan, "--do_dbscan", "do dbscan groupings rather than single linkage clustering");
  setUp.setOption(dbscanPars.epsilon_inclusive_, "--dbscan_epsilon_inclusive", "dbscan epsilon inclusive");
  setUp.setOption(dbscanPars.min_neighbors_inclusive_, "--dbscan_min_neighbors_inclusive", "dbscan min neighbors inclusive");
  setUp.setOption(shuffle_seed, "--shuffle_seed", "seed for random shuffle of barcodes");
  setUp.processDefaultReader(true);
  setUp.setOption(numThreads, "--numThreads", "number of threads to use");
  setUp.finishSetUp(std::cout);

  SeqIO seq_io(setUp.pars_.ioOptions_);
  seq_io.openIn();
  seq_io.openOut();
  auto library = seq_io.in_.readAllReadsPtrs<seqInfo>();


  //generating barcodes with entropy cut off to create non-low complexity barcodes
  std::vector<std::string> barcodes = PhipSeqNucLibraryGraph::generate_random_barcodes(random_barcodes_pars);
  {
    VecStr barcode_sanity_checks;
    if (barcodes.empty()) {
      barcode_sanity_checks.emplace_back(njh::pasteAsStr("could not create any barcodes with current parameters for size: ", random_barcodes_pars.barcode_size));
    }
    for (const auto & barcode : barcodes) {
      if (barcode.size() != random_barcodes_pars.barcode_size) {
        barcode_sanity_checks.emplace_back(njh::pasteAsStr("barcode: ", barcode, " is sized: ", barcode.size(), "not barcode_size: ", random_barcodes_pars.barcode_size));
      }
    }
    if (!barcode_sanity_checks.empty()) {
      std::stringstream ss;
      ss << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << ", error " << " " << "\n";
      ss << njh::conToStr(barcode_sanity_checks, "\n") << "\n";
      throw std::runtime_error{ss.str()};
    }
  }

  std::unordered_map<std::string, uint32_t> name_counts;
  for (const auto & seq : library) {
    ++name_counts[seq->name_];
    if (len(*seq) != len(*library.front())) {
      std::stringstream ss;
      ss << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << ", error "
      << ", all library has to be the same length, seq: " << seq->name_
      << " has len: " << len(*seq) << " which is different from "
      << library.front()->name_ << " which has len " << len(*library.front()) << "\n";
      throw std::runtime_error{ss.str()};
    }
  }
  {
    VecStr multiple_names;
    for (const auto & name_count : name_counts) {
      if (name_count.second > 1) {
        multiple_names.emplace_back(name_count.first);
      }
    }
    if (!multiple_names.empty()) {
      std::stringstream ss;
      ss << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << ", error " << " the following names were found multiple times, all names must be unique" << "\n";
      ss << njh::conToStr(multiple_names, ",") << "\n";
      throw std::runtime_error{ss.str()};
    }
  }
  PairwisePairFactory pairFactory(library.size());
  //set up progress bar
  njh::ProgressBar pBar(pairFactory.totalCompares_);

  njhseq::concurrent::AlignerPool aligner_pool(
    len(*library.front()) *2, gapScoringParameters(10, 1), substituteMatrix::createScoreMatrix(1, 0, false, false, true), numThreads
  );
  aligner_pool.initAligners();

  PhipSeqNucLibraryGraph libraryGraph;
  for (const auto & seq : library) {
    libraryGraph.nodes_.emplace_back(seq);
  }
  std::unordered_map<uint32_t, uint32_t> counts;
  std::mutex counts_mut;
  std::function<void()> getEditDistances = [&counts,&counts_mut,&pairFactory,&library,&aligner_pool,&setUp, &pBar,
    &hamming_distance_cut_off_exclusive, &libraryGraph]() {
    PairwisePairFactory::PairwisePairVec pairs;
    auto current_aligner = aligner_pool.popAligner();
    std::vector<std::shared_ptr<PhipSeqNucLibraryGraph::edge>> current_edges;
    std::unordered_map<uint32_t, uint32_t> current_counts;
    while (pairFactory.setNextPairs(pairs, 1000)) {
      if (setUp.pars_.verbose_) {
        pBar.outputProgAdd(std::cout, pairs.pairs_.size(), true);
      }
      for (const auto & pair : pairs.pairs_) {
        current_aligner->noAlignSetAndScore(library[pair.col_], library[pair.row_]);
        auto edit_dist = library[pair.col_]->seq_.size() - current_aligner->parts_.score_;
        ++current_counts[edit_dist];
        if (edit_dist < hamming_distance_cut_off_exclusive) {
          //add to edges to be added later
          current_edges.emplace_back(std::make_shared<PhipSeqNucLibraryGraph::edge>(
            std::unordered_map<std::string, uint32_t>{
              {library[pair.col_]->name_, pair.row_},
              {library[pair.row_]->name_, pair.col_}
            }, edit_dist));
        }
      }
    }
    {
      std::lock_guard lock(counts_mut);
      for (const auto & count : current_counts) {
        counts[count.first] += count.second;
      }
      for ( auto & e : current_edges) {
        //add to nodes
        auto node_positions = getVectorOfMapValues(e->name_to_other_node_pos_);
        for (const auto node_pos : node_positions) {
          libraryGraph.nodes_[node_pos].edges_.emplace_back(e);
        }
        //add to edges
        libraryGraph.all_edges_.emplace_back(e);
      }
    }
  };

  njh::concurrent::runVoidFunctionThreaded(getEditDistances, numThreads);
  std::unordered_map<uint32_t, uint32_t> neighbors_counts;
  std::vector<uint64_t> nodes_with_neighbors;
  for (const auto n_pos : iter::range(libraryGraph.nodes_.size())) {
    const auto & n = libraryGraph.nodes_[n_pos];
    ++neighbors_counts[n.edges_.size()];
    if (!n.edges_.empty()) {
      nodes_with_neighbors.emplace_back(n_pos);
    }
  }
  //
  njh::randomGenerator rgen;
  if (std::numeric_limits<uint64_t>::max() != shuffle_seed) {
    rgen.seedNum(shuffle_seed);
  }
  njh::shuffle(barcodes, rgen.mtGen_);
  //first ensure that all nodes with neighbor have distinct barcodes
  std::deque<std::string> barcodes_deque(barcodes.begin(), barcodes.end());
  for (const auto nodes_pos : nodes_with_neighbors) {
    //check to see if we have used up all barcodes, if so replenish
    if (barcodes_deque.empty()) {
      njh::shuffle(barcodes, rgen.mtGen_);
      barcodes_deque = std::deque<std::string>(barcodes.begin(), barcodes.end());
    }
    if (libraryGraph.nodes_[nodes_pos].barcode_.empty()) {
      VecStr neighbor_barcodes;
      for (const auto & e : libraryGraph.nodes_[nodes_pos].edges_) {
        if (!libraryGraph.nodes_[e->name_to_other_node_pos_[libraryGraph.nodes_[nodes_pos].seqBase_->name_]].barcode_.empty()) {
          neighbor_barcodes.emplace_back(libraryGraph.nodes_[e->name_to_other_node_pos_[libraryGraph.nodes_[nodes_pos].seqBase_->name_]].barcode_);
        }
      }
      auto potential_barcode = barcodes_deque.front();
      barcodes_deque.pop_front();
      VecStr failed_barcodes;
      bool matched = false;
      for (const auto & neighbor_barcode : neighbor_barcodes) {
        if (hamming_distance_no_check(potential_barcode, neighbor_barcode) < random_barcodes_pars.barcode_hamming_distance_cut_off_exclusive) {
          matched = true;
          break;
        }
      }
      uint32_t replenish_count = 0;

      while (matched && !barcodes_deque.empty()) {
        failed_barcodes.emplace_back(potential_barcode);
        potential_barcode = barcodes_deque.front();
        barcodes_deque.pop_front();
        matched = false;
        for (const auto & neighbor_barcode : neighbor_barcodes) {
          if (hamming_distance_no_check(potential_barcode, neighbor_barcode) < random_barcodes_pars.barcode_hamming_distance_cut_off_exclusive) {
            matched = true;
            break;
          }
        }
        //if still matching, deque is empty but we haven't replenished yet, first try replenishing
        if (matched && barcodes_deque.empty() && 0 == replenish_count ) {
          ++replenish_count;
          njh::shuffle(barcodes, rgen.mtGen_);
          barcodes_deque = std::deque<std::string>(barcodes.begin(), barcodes.end());
          //since we are adding all barcode, we can reset the failed barcodes
          failed_barcodes.clear();
        }
      }

      //if still matching means the deque was run through
      if (matched) {
        std::stringstream ss;
        ss << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << ", error " << " failed to find a properly distant barcode " << "\n";
        ss << "barcode_pos: " << nodes_pos << "\n";
        ss << "neighbor_barcodes: " << njh::conToStr(neighbor_barcodes, ",") << '\n';
        ss << "failed_barcodes.size(): " << failed_barcodes.size() << "\n";
        throw std::runtime_error{ss.str()};
      }
      // set barcode, add the failed barcodes to the end of the deque
      libraryGraph.nodes_[nodes_pos].barcode_ = potential_barcode;
      if (!failed_barcodes.empty()) {
        barcodes_deque.insert(barcodes_deque.end(), failed_barcodes.begin(), failed_barcodes.end());
      }
    }
  }
  //now set the reset of neighborless nodes
  for (const auto node_pos : iter::range(libraryGraph.nodes_.size())) {
    //check to see if we have used up all barcodes, if so replenish
    if (barcodes_deque.empty()) {
      njh::shuffle(barcodes, rgen.mtGen_);
      barcodes_deque = std::deque<std::string>(barcodes.begin(), barcodes.end());
    }
    if (libraryGraph.nodes_[node_pos].barcode_.empty()) {
      if (!libraryGraph.nodes_[node_pos].edges_.empty()) {
        std::stringstream ss;
        ss << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << ", error " << " nodes with neighbors should have been set already" << "\n";
        throw std::runtime_error{ss.str()};
      }
      libraryGraph.nodes_[node_pos].barcode_ = barcodes_deque.front();
      barcodes_deque.pop_front();
    }
  }
  if (mark_with_groups) {
    if (do_dbscan) {
      libraryGraph.determine_groups_dbscan(dbscanPars);
    }else {
      libraryGraph.determine_groups_single_linkage();
    }
  }
  for (const auto & node : libraryGraph.nodes_) {

    //add the random barcode and write out
    node.seqBase_->append(random_barcodes_pars.prepend_to_barcode + node.barcode_);
    if (mark_with_groups) {
      MetaDataInName meta;
      if (MetaDataInName::nameHasMetaData(node.seqBase_->name_)) {
        meta = MetaDataInName(node.seqBase_->name_);
      }
      meta.addMeta("hamming_distanace_group", node.group_, true);
      meta.resetMetaInName(node.seqBase_->name_);
    }
    seq_io.write(node.seqBase_);
  }

  if (setUp.pars_.debug_) {
    std::cerr << "neighbor counts:" << std::endl;
    std::cerr << "neighbors\tcount" << std::endl;
    auto counts_key = njh::getSetOfMapKeys(neighbors_counts);
    for (const auto & dist : counts_key) {
      std::cerr << dist << "\t" << neighbors_counts[dist] << std::endl;
    }
  }
  if (setUp.pars_.debug_) {
    std::map<uint32_t, uint32_t> hamming_dists_counts_between_neighbors;
    std::cout << "edges:" << std::endl;
    for (const auto & e : libraryGraph.all_edges_) {
      auto node_names = getVectorOfMapKeys(e->name_to_other_node_pos_);
      njh::sort(node_names);
      // std::cerr << njh::conToStr(node_names, " -> ") << ": " << e->hamming_dist_ << std::endl;
      auto nodes_positions = njh::getVecOfMapValues(e->name_to_other_node_pos_);
      ++hamming_dists_counts_between_neighbors[hamming_distance_no_check(
        libraryGraph.nodes_[nodes_positions.front()].barcode_,
        libraryGraph.nodes_[nodes_positions.back()].barcode_)];
    }
    std::cout << "hamming_distance_counts: " << std::endl;
    for (const auto & ham_d : hamming_dists_counts_between_neighbors) {
      std::cout << ham_d.first << '\t' << ham_d.second << std::endl;
    }
  }
  return 0;
}



int PhipSeqUtilsRunner::countPossiblePhipSeqRandomBarcodes(const njh::progutils::CmdArgs & libraryCommands) {
  PhipSeqNucLibraryGraph::generate_random_barcodes_pars random_barcodes_pars;

  OutOptions outOpts("", ".tsv");
  uint32_t barcode_size_start = 6;
  uint32_t barcode_size_end = 10;
  seqSetUp setUp(libraryCommands);
  setUp.description_ = "get the number of random barcodes sequence made of AGTC with various filters";
  setUp.processVerbose();
  setUp.processDebug();
  random_barcodes_pars.set_options(setUp, false);
  setUp.processWritingOptions(outOpts);

  setUp.setOption(barcode_size_start, "--barcode_size_start", "barcode_size_start", true);
  setUp.setOption(barcode_size_end, "--barcode_size_end", "barcode_size_end", true);
  if (barcode_size_end < barcode_size_start) {
    setUp.failed_= true;
    setUp.addWarning(njh::pasteAsStr("--barcode_size_end: ", barcode_size_end, " cannot not be less than --barcode_size_start: ", barcode_size_start));
  }
  setUp.finishSetUp(std::cout);
  OutputStream out(outOpts);
  out << "barcode_size\tbarcode_counts" << std::endl;
  for (const auto & barcode_size : iter::range(barcode_size_start, barcode_size_end + 1)) {
    //generating barcodes with entropy cut off to create non-low complexity barcodes
    random_barcodes_pars.barcode_size = barcode_size;
    std::vector<std::string> barcodes = PhipSeqNucLibraryGraph::generate_random_barcodes(random_barcodes_pars);
    out << barcode_size << "\t" << barcodes.size() << std::endl;
  }
  return 0;
}



} //namespace njhseq
