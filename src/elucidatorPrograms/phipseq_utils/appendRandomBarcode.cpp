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
  class edge {
  public:
    edge(std::unordered_map<std::string, uint32_t> key, const uint32_t edit_dist) : name_to_other_node_pos_(std::move(key)),
                                                                              edit_dist_(edit_dist) {
    }

    std::unordered_map<std::string, uint32_t> name_to_other_node_pos_;//<! index is current node, value is the other node position
    uint32_t edit_dist_;
  };
  class node {
    public:
    node(const std::shared_ptr<seqInfo> & seq): seqBase_(seq) {

    }
    std::shared_ptr<seqInfo> seqBase_;
    std::vector<std::shared_ptr<edge>> edges_;//<! other sequences below a specific edit distance away
  };

  std::vector<node> nodes_;
  std::vector<std::shared_ptr<edge>> all_edges_;
};




int PhipSeqUtilsRunner::countPossiblePhipSeqRandomBarcodes(const njh::progutils::CmdArgs & libraryCommands) {
  std::vector<char> barcode_alphabet{'A', 'G', 'C', 'T'};
  double barcode_entropy_filter = 0.52;
  uint32_t barcode_entropy_klen = 3;
  uint32_t barcode_size_start = 3;
  uint32_t barcode_size_end = 15;

  uint32_t edit_distance_cut_off_exclusive = 3;
  std::string linker_seq;

  OutOptions outOpts("", ".tsv");

  seqSetUp setUp(libraryCommands);
  setUp.description_ = "get the number of random barcodes sequence made of AGTC with various filters";
  setUp.processVerbose();
  setUp.processDebug();

  setUp.processWritingOptions(outOpts);
  setUp.setOption(linker_seq, "--linker_seq", "linker seq to compare to ensure to have a random barcode that is similar");

  setUp.setOption(barcode_size_start, "--barcode_size_start", "barcode_size_start", true);
  setUp.setOption(barcode_size_end, "--barcode_size_end", "barcode_size_end", true);
  if (barcode_size_end < barcode_size_start) {
    setUp.failed_= true;
    setUp.addWarning(njh::pasteAsStr("--barcode_size_end: ", barcode_size_end, " cannot not be less than --barcode_size_start: ", barcode_size_start));
  }

  setUp.setOption(edit_distance_cut_off_exclusive, "--edit_distance_cut_off_exclusive", "edit_distance_cut_off_exclusive", true);
  setUp.setOption(barcode_entropy_klen, "--barcode_entropy_klen", "barcode entropy klen for calculating entropy");
  setUp.setOption(barcode_entropy_filter, "--barcode_entropy_filter", "barcode entropy filter");
  setUp.setOption(barcode_alphabet, "--barcode_alphabet", "alphabet for barcodes");
  setUp.finishSetUp(std::cout);


  OutputStream out(outOpts);
  out << "barcode_size\tbarcode_counts" << std::endl;
  std::unordered_map<uint32_t, uint32_t> barcode_counts;
  for (const auto & barcode_size : iter::range(barcode_size_start, barcode_size_end + 1)) {
    //generating barcodes with entropy cut off to create non-low complexity barcodes
    // std::vector<std::string> barcodes;
    uint32_t barcode_count = 0;
    generate_combinations(barcode_alphabet, barcode_size, [&barcode_count,&barcode_entropy_klen,&barcode_entropy_filter,
      &linker_seq, &edit_distance_cut_off_exclusive](const std::string & barcode) {
      auto kinfo = kmerInfo(barcode, barcode_entropy_klen, false);
      if (kinfo.computeKmerEntropy() >= barcode_entropy_filter) {
        if (!linker_seq.empty()) {
          uint32_t edit_dist = 0;
          //just checking the beginning
          for (const auto pos : iter::range(std::min(linker_seq.size(), barcode.size()))) {
            if (linker_seq[pos] != barcode[pos]) {
              ++edit_dist;
            }
          }
          if (edit_dist >= edit_distance_cut_off_exclusive) {
            // barcodes.emplace_back(barcode);
            ++barcode_count;
          }
        } else {
          // barcodes.emplace_back(barcode);
          ++barcode_count;
        }
      }
    });

    out << barcode_size << "\t" << barcode_count << std::endl;
  }


  return 0;
}


int PhipSeqUtilsRunner::appendRandomBarcode(const njh::progutils::CmdArgs & libraryCommands) {
  std::vector<char> barcode_alphabet{'A', 'G', 'C', 'T'};
  double barcode_entropy_filter = 0.52;
  uint32_t barcode_entropy_klen = 3;
  uint32_t barcode_size = 6;
  uint32_t number_of_stop_codons = 2;
  uint32_t edit_distance_cut_off_exclusive = 4;

  uint32_t numThreads = 1;
  OutOptions outOpts("", ".tsv");

  seqSetUp setUp(libraryCommands);
  setUp.description_ = "append a random barcode sequence made of AGTC to the end of a phipeq library";
  setUp.processVerbose();
  setUp.processDebug();
  // setUp.processDefaultReader(true);
  setUp.processWritingOptions(outOpts);
  setUp.processReadInNames(true);
  setUp.setOption(barcode_size, "--barcode_size", "barcode size", true);
  setUp.setOption(edit_distance_cut_off_exclusive, "--edit_distance_cut_off_exclusive", "edit_distance_cut_off_exclusive", true);
  setUp.setOption(barcode_entropy_klen, "--barcode_entropy_klen", "barcode entropy klen for calculating entropy");
  setUp.setOption(barcode_entropy_filter, "--barcode_entropy_filter", "barcode entropy filter");

  setUp.setOption(number_of_stop_codons, "--number_of_stop_codons", "number of stop codons");
  setUp.setOption(numThreads, "--numThreads", "number of threads to use");
  setUp.setOption(barcode_alphabet, "--barcode_alphabet", "alphabet for barcodes");

  setUp.finishSetUp(std::cout);
  SeqIO seq_io(setUp.pars_.ioOptions_);
  seq_io.openIn();
  auto library = seq_io.in_.readAllReadsPtrs<seqInfo>();
  OutputStream out(outOpts);

  //generating barcodes with entropy cut off to create non-low complexity barcodes
  std::vector<std::string> barcodes;
  generate_combinations(barcode_alphabet, barcode_size, [&barcodes,&barcode_entropy_klen,&barcode_entropy_filter](const std::string & barcode) {
    auto kinfo = kmerInfo(barcode, barcode_entropy_klen, false);
    if (kinfo.computeKmerEntropy() >= barcode_entropy_filter) {
      barcodes.emplace_back(barcode);
    }
  });

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
    &edit_distance_cut_off_exclusive, &libraryGraph]() {
    PairwisePairFactory::PairwisePair pair;
    auto current_aligner = aligner_pool.popAligner();
    std::vector<std::shared_ptr<PhipSeqNucLibraryGraph::edge>> current_edges;
    std::unordered_map<uint32_t, uint32_t> current_counts;
    while (pairFactory.setNextPair(pair)) {
      if (setUp.pars_.verbose_) {
        pBar.outputProgAdd(std::cout, 1, true);
      }
      current_aligner->noAlignSetAndScore(library[pair.col_], library[pair.row_]);
      auto edit_dist = library[pair.col_]->seq_.size() - current_aligner->parts_.score_;
      ++current_counts[edit_dist];
      if (edit_dist < edit_distance_cut_off_exclusive) {
        //add to edges to be added later
        current_edges.emplace_back(std::make_shared<PhipSeqNucLibraryGraph::edge>(
          std::unordered_map<std::string, uint32_t>{
            {library[pair.col_]->name_, pair.row_},
            {library[pair.row_]->name_, pair.col_}
          }, edit_dist));
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
  for (const auto & n : libraryGraph.nodes_) {
    ++neighbors_counts[n.edges_.size()];
  }

  out << "neighbors\tcount" << std::endl;
  auto counts_key = njh::getSetOfMapKeys(neighbors_counts);
  for (const auto & dist : counts_key) {
    out << dist << "\t" << neighbors_counts[dist] << std::endl;
  }
  if (setUp.pars_.debug_) {
    for (const auto & e : libraryGraph.all_edges_) {
      auto node_names = getVectorOfMapKeys(e->name_to_other_node_pos_);
      njh::sort(node_names);
      std::cout << njh::conToStr(node_names, " -> ") << ": " << e->edit_dist_ << std::endl;
    }
  }
  return 0;
}

} //namespace njhseq
