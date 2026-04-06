//
// Created by Nicholas Hathaway on 4/3/26.
//
#include "kmerExp.hpp"

#include <njhseq/alignment/aligner/aligner.hpp>
#include <njhseq/concurrency/PairwisePairFactory.hpp>
#include <njhseq/concurrency/pools/AlignerPool.hpp>


#include "elucidator/objects/BioDataObject.h"
#include <njhseq/objects/dataContainers/graphs/UndirWeightedGraph.hpp>
#include <njhseq/objects/seqObjects/seqKmers/KmerVecUtils.hpp>
#include <njhseq/PopulationGenetics/PopGenCalcs.hpp>
#include <njhseq/objects/helperObjects/PeptideLibraryReducer.hpp>


namespace njhseq {

int kmerExpRunner::encodingSeqsByCommonKmers(const njh::progutils::CmdArgs &inputCommands) {
  njh::OutOptions keyOutOpts;
  uint32_t klen = 9;
  seqSetUp setUp(inputCommands);
  setUp.processVerbose();
  setUp.processDebug();
  setUp.processDefaultReader(true);
  setUp.setOption(klen, "--klen", "klen");
  setUp.setOption(keyOutOpts.outFilename_, "--outKey", "output a file with the reduction key");
  keyOutOpts.transferOverwriteOpts(setUp.pars_.ioOptions_.out_);
  setUp.finishSetUp(std::cout);

  // SeqInput
  // @todo, implement, find first the highest counted kmer and then next highest kmer (after masking), and so on and so forth until no more multiples are found

  return 0;
}
} // namespace njhseq

