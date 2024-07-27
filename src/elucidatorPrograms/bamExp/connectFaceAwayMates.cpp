//
// Created by Nicholas Hathaway on 7/26/24.
//


#include "bamExp.hpp"
#include <TwoBit.h>

#include "elucidator/BamToolsUtils.h"
#include <njhseq/objects/BioDataObject.h>


namespace njhseq {


class GenomicRegionGraph {


  public:

  class edge {
  public:
    //key is genomic UID and value is the position of the other node
    std::unordered_map<std::string, uint32_t> nodes_;
    std::unordered_set<std::string> samples_;
    std::unordered_set<std::string> reads_;
  };

  class node {
    public:

    node(GenomicRegion region): region_(std::move(region)) {

    }

    GenomicRegion region_;

    uint32_t group_{std::numeric_limits<uint32_t>::max()};
    std::vector<std::shared_ptr<edge>> edges_;
  };


  void addNode(const GenomicRegion & region) {
    for(const auto & n : nodes_) {
      if(n.region_.uid_ == region.uid_) {
        std::stringstream ss;
        ss << __PRETTY_FUNCTION__ << ", error " << "already have node with UID: " << region.uid_ << "\n";
        throw std::runtime_error{ss.str()};
      }
    }
    nodes_.emplace_back(region);
  }

  std::vector<node> nodes_;
  std::vector<std::shared_ptr<edge>> edges_;
};



int bamExpRunner::connectFaceAwayMatesRegions(
  const njh::progutils::CmdArgs &inputCommands) {

  OutOptions outOpts(bfs::path(""), ".bed");
  bfs::path bedFnp;
  uint32_t minReadAmount = 10;
  uint32_t minInsertSize = 1000;
  uint32_t step = 50;
  uint32_t windowSize = 100;
  uint32_t minLen = 0;
  std::string sample = "sample";
  seqSetUp setUp(inputCommands);

  setUp.setOption(sample, "--sample", "sample");

  setUp.setOption(step, "--step", "step");
  setUp.setOption(windowSize, "--windowSize", "windowSize");
  minLen = windowSize;
  setUp.setOption(minLen, "--minLen", "minimum length");
  setUp.processVerbose();
  setUp.setOption(bedFnp, "--bedFnp", "Bed file to gather the face away reads on", true);
  setUp.setOption(minReadAmount, "--minReadAmount", "min Read Amount");
  setUp.setOption(minInsertSize, "--minInsertSize", "minimum Insert Size to include in merge");


  setUp.processReadInNames( { "--bam" }, true);
  setUp.processWritingOptions(outOpts);
  setUp.finishSetUp(std::cout);

  GenomicRegionGraph graph;

  auto inRegions = getBed3s(bedFnp);
  auto mergedRegions = BedUtility::mergeAndSort(inRegions);

  std::vector<Bed3RecordCore> faceawayRegions;
  {
    BamTools::BamReader bReader;
    bReader.Open(setUp.pars_.ioOptions_.firstName_.string());
    checkBamOpenThrow(bReader, setUp.pars_.ioOptions_.firstName_.string());
    loadBamIndexThrow(bReader);
    auto refData = bReader.GetReferenceData();
    BamTools::BamAlignment bAln;
    for (const auto& reg: mergedRegions) {
      setBamFileRegionThrow(bReader, reg);
      std::vector<Bed6RecordCore> otherRegions;

      while (bReader.GetNextAlignment(bAln)) {
        if (bAln.IsMapped() && bAln.IsMateMapped() &&
            bAln.IsReverseStrand() == bAln.IsMateReverseStrand() &&
            std::abs(bAln.InsertSize) >= minInsertSize) {
          auto bAlnRegion = GenomicRegion(bAln, refData);
          auto bAlnMateRegion = GenomicRegion(bAln.Name,
                                              refData[bAln.MateRefID].RefName, bAln.MatePosition,
                                              bAln.MatePosition + bAln.AlignedBases.size(), bAln.IsMateReverseStrand());
          faceawayRegions.emplace_back(bAlnRegion.genBed3RecordCore());
          faceawayRegions.emplace_back(bAlnMateRegion.genBed3RecordCore());

            }
      }
    }
  }
  auto faceawayRegions_merged = BedUtility::mergeAndSort(faceawayRegions);

  for(const auto & reg : mergedRegions){
    auto windows = BedUtility::createWindowsWithinRegion(reg, windowSize, step, true);
    for(const auto & win : windows) {
      bool overlapWithFace = false;
      for(const auto & face : faceawayRegions_merged) {
        //windows are sorted so if beyond this current window then it won't be in the vector of windows
        if(face.chrom_ > win.chrom_  || face.start_ > face.end_) {
          break;
        }
        if(win.overlaps(face, 1)) {
          overlapWithFace = true;
          break;
        }
      }
      if(overlapWithFace) {
        graph.addNode(win);
      }
    }
  }

  {
    BamTools::BamReader bReader;
    bReader.Open(setUp.pars_.ioOptions_.firstName_.string());
    checkBamOpenThrow(bReader, setUp.pars_.ioOptions_.firstName_.string());
    loadBamIndexThrow(bReader);
    auto refData = bReader.GetReferenceData();
    BamTools::BamAlignment bAln;
    for (const auto& reg: mergedRegions) {
      setBamFileRegionThrow(bReader, reg);
      std::vector<Bed6RecordCore> otherRegions;

      while (bReader.GetNextAlignment(bAln)) {
        if (bAln.IsMapped() && bAln.IsMateMapped() &&
            bAln.IsReverseStrand() == bAln.IsMateReverseStrand() &&
            std::abs(bAln.InsertSize) >= minInsertSize) {
          auto bAlnRegion = GenomicRegion(bAln, refData);
          auto bAlnMateRegion = GenomicRegion(bAln.Name,
                                              refData[bAln.MateRefID].RefName, bAln.MatePosition,
                                              bAln.MatePosition + bAln.AlignedBases.size(), bAln.IsMateReverseStrand());
          std::unordered_set<uint32_t> baln_nodes;
          std::unordered_set<uint32_t> mate_nodes;
          for(const auto & nEnum : iter::enumerate(graph.nodes_)) {
            if(bAlnRegion.overlaps(nEnum.element.region_)) {
              baln_nodes.emplace(nEnum.index);
            }
            if(bAlnMateRegion.overlaps(nEnum.element.region_)) {
              mate_nodes.emplace(nEnum.index);
            }
          }
          if(!baln_nodes.empty()) {
            for(const auto & balnNode: baln_nodes) {
              for(const auto & mateNode : mate_nodes) {
                bool foundEdge = false;
                for(const auto & e : graph.nodes_[balnNode].edges_) {
                  if(njh::in(graph.nodes_[mateNode].region_.uid_, e->nodes_)) {
                    foundEdge = true;
                    e->reads_.emplace(bAln.Name);
                    e->samples_.emplace(sample);
                    break;
                  }
                }
                if(!foundEdge) {
                  auto newEdge = std::make_shared<GenomicRegionGraph::edge>();
                  newEdge->nodes_[graph.nodes_[mateNode].region_.uid_] = balnNode;
                  newEdge->nodes_[graph.nodes_[balnNode].region_.uid_] = mateNode;
                  newEdge->reads_.emplace(bAln.Name);
                  newEdge->samples_.emplace(sample);
                  graph.nodes_[mateNode].edges_.emplace_back(newEdge);
                  graph.nodes_[balnNode].edges_.emplace_back(newEdge);
                  graph.edges_.emplace_back(newEdge);
                }
              }
            }
          }


            }
      }
    }
  }

  std::cout << "regions\tread_count\tsample_count" << std::endl;
  for(const auto & e : graph.edges_) {
    std::cout << njh::conToStr(getVectorOfMapKeys(e->nodes_), ",")
      << "\t" << e->reads_.size()
    << "\t" << e->samples_.size() << std::endl;
  }

  return 0;
}
} //  namespace njhseq {

