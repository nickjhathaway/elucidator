//
// Created by Nicholas Hathaway on 1/6/26.
//



#include "popGenExp.hpp"
#include "elucidator/objects/BioDataObject.h"
#include <njhseq/objects/counters/DNABaseCounter.hpp>


#include "elucidator/seqToolsUtils/seqToolsUtils.hpp"
#include <njhseq/GenomeUtils.h>
#include <njhseq/objects/Gene/GeneFromGffs.hpp>

#include "elucidator/BioRecordsUtils/BedUtility.hpp"
#include <njhseq/objects/Gene/GenomicAminoAcidPositionTyper.hpp>
#include <njhseq/objects/Gene/TranslatorByAlignment.hpp>

#include <njhseq/objects/seqObjects/Clusters/identicalCluster.hpp>
#include <njhseq/concurrency/pools/BamReaderPool.hpp>
#include <njhseq/concurrency/pools/AlignerPool.hpp>
#include <njhseq/PopulationGenetics.h>
#include <njhseq/objects/seqContainers/CollapsedHaps.hpp>



namespace njhseq {


int popGenExpRunner::calc_lins_concordance_correlation(const njh::progutils::CmdArgs &inputCommands) {
  bfs::path table_fnp;

  std::string col1;
  std::string col2;

  uint32_t bootstraps = 2000;
  double conf_level = 0.95;
  bool bootstrap_ci = false;
  std::string ci_method = "z-transform";
  std::uint64_t bootstrap_seed = 0;
  seqSetUp setUp(inputCommands);
  setUp.processVerbose();
  setUp.processDebug();

  setUp.setOption(table_fnp, "--table_fnp", "The table to read in", true);
  setUp.setOption(col1, "--col1", "column 1 to calculate from", true);
  setUp.setOption(col2, "--col2", "column 2 to calculate from", true);
  setUp.setOption(bootstrap_ci, "--bootstrap_ci", "bootstrap_ci");
  setUp.setOption(bootstrap_seed, "--bootstrap_seed", "bootstrap seed");

  if (bootstrap_ci) {
    ci_method = "bca";
  }
  setUp.setOption(bootstraps, "--bootstraps", "bootstraps");

  setUp.setOption(conf_level, "--conf_level", "conf_level");

  setUp.finishSetUp(std::cout);

  table input(table_fnp, "\t", true);
  input.checkForColumnsThrow(VecStr{col1, col2}, __PRETTY_FUNCTION__);
  std::vector<double> col1_values = vecStrToVecNum<double>(input.getColumn(col1));
  std::vector<double> col2_values = vecStrToVecNum<double>(input.getColumn(col2));
  if (setUp.pars_.debug_) {
    auto ccc = lins_concordance_correlation(col1_values, col2_values);
    std::cout << ccc << std::endl;
  }

  auto ccc_confidence = ConcordanceCalculator::lins_ccc_with_ci(col1_values, col2_values, conf_level, ci_method,
                                                                bootstraps, bootstrap_seed);
  std::cout << ccc_confidence.ccc << " (" << ccc_confidence.conf_level * 100 << "% CI " << ccc_confidence.lower << "-"
      << ccc_confidence.upper << ")" << std::endl;

  return 0;
}

}//namespace njhseq
