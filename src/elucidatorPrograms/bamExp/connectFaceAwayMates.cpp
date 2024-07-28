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
    bool on_{true};
  };

  class node {
    public:

    explicit node(GenomicRegion region): region_(std::move(region)) {

    }

    GenomicRegion region_;

    uint32_t group_{std::numeric_limits<uint32_t>::max()};
    std::vector<std::shared_ptr<edge>> edges_;

    [[nodiscard]] uint32_t getOnEdgeCount() const {
      uint32_t count = 0;
      for(const auto& edge : edges_) {
        if(edge->on_) {
          ++count;
        }
      }
      return count;
    }
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


  void resetNodes() {
    for(auto & n : nodes_) {
      n.group_ = std::numeric_limits<uint32_t>::max();
    }
  }

  void turnOffEdgesReadCountCutOff(const uint32_t readCountCutOff) {
    for(auto & e : edges_) {
      if(e->reads_.size() < readCountCutOff) {
        e->on_ = false;
      }
    }
  }
private:
  void setNodeGroupAndSpread(const uint32_t nodePos, const uint32_t currentGroup) {
    std::stack<uint32_t> addingNodes;
    addingNodes.push(nodePos);
    nodes_[nodePos].group_ = currentGroup;
    for (const auto& e : nodes_[nodePos].edges_) {
      if(e->on_ && nodes_[e->nodes_[nodes_[nodePos].region_.uid_]].group_ == std::numeric_limits<uint32_t>::max()) {
        addingNodes.emplace(e->nodes_[nodes_[nodePos].region_.uid_]);
      }
    }
    while (!addingNodes.empty()) {
      auto currentNodePos = addingNodes.top();
      addingNodes.pop();
      nodes_[currentNodePos].group_ = currentGroup;
      for (const auto& e : nodes_[currentNodePos].edges_) {
        auto nextNodePos = e->nodes_[nodes_[currentNodePos].region_.uid_];
        if (e->on_  && nodes_[nextNodePos].group_ == std::numeric_limits<uint32_t>::max()) {
          addingNodes.emplace(nextNodePos);
        }
      }
    }
  }

public:
  void determineGroups() {
    resetNodes();
    uint32_t currentGroup = 0;
    for(const auto & n : iter::enumerate(nodes_)) {
      if(n.element.group_ == std::numeric_limits<uint32_t>::max()) {
        setNodeGroupAndSpread(n.index, currentGroup);
        ++currentGroup;
      }
    }
  }

  [[nodiscard]] std::map<uint32_t, uint32_t> getGroupCounts() const {
    std::map<uint32_t, uint32_t> groups;
    for(const auto & n : iter::enumerate(nodes_)) {
      ++groups[n.element.group_];
    }
    return groups;
  }

  void writeGraphViz(const OutOptions & outOpts) {
    OutputStream out(outOpts);

	  VecStr colors = {"#006E82","#8214A0","#005AC8","#00A0FA","#FA78FA","#14D2DC","#AA0A3C","#FA7850","#0AB45A","#F0F032","#A0FA82","#FAE6BE",
	    "#ff358f", "#07c652", "#a2009b", "#467f00", "#cb73fd",
			  "#f4be45", "#0157d8", "#ff8d3b", "#0056a1", "#dd0d3a", "#01d5e4",
			  "#b10049", "#7cda97", "#ff76e0", "#018a5a", "#ff87b8", "#4a5b00",
			  "#664092", "#8f7400", "#02aee5", "#9e3500", "#8bd5b8", "#8a306b",
			  "#e4b7f0", "#ff97a2" };

	  out << "digraph graphname {" << std::endl;
	  out << "\t" << "node [fixedsize=true regular=true shape=ellipse]" << std::endl;
    auto groupCounts = getGroupCounts();

	  if(groupCounts.size() >= colors.size()){
		  auto hColors = njh::heatmapColors(groupCounts.size() + 1);
		  njh::reverse(hColors);
		  colors.clear();
		  for(const auto & hColor : hColors){
			  colors.emplace_back(hColor.getHexStr());
		  }
		  out << "\t"
				  << "graph [ bgcolor=black, resolution=128, fontname=Arial, fontcolor=white,  fontsize=12 ]; "
				  << std::endl;
		  out << "\t" << "node [ fontname=Arial, fontcolor=white, fontsize=11];"
				  << std::endl;
		  out << "\t" << "edge [ fontname=Helvetica, fontcolor=white, fontsize=10 ];"
				  << std::endl;
	  }
	  for(const auto & node : iter::enumerate(nodes_)){
		  double diameter = 2 * std::sqrt(node.element.getOnEdgeCount()/M_PI);
		  const std::string& nodeColor = colors[node.element.group_];
		  //out << "\t" << node->uid_  <<"[fixedsize=true,shape=circle,width=" << diameter << ",style=filled,fillcolor=\"" << nodeColor << "\", label=\"" << node->uid_ << "-" << node->inReadNamesIdx_.size()<< "\"]"<< std::endl;
		  out << "\t" << node.index  <<"[fixedsize=true,shape=circle,width=" << diameter << ",style=filled,fillcolor=\"" << nodeColor << "\", label=\"" << node.element.region_.uid_ << "-" << node.element.edges_.size() << "\"]"<< std::endl;
	  }
    for(const auto & e : edges_) {
      if(e->on_ ) {
        out << "\t" << e->nodes_.begin()->second << " -> " << e->nodes_[nodes_[e->nodes_.begin()->second].region_.uid_] << "[penwidth=" << std::log10(e->reads_.size() + 1) << ", label=\"" << e->reads_.size() << "\""<< "]"   << std::endl;
      }
    }
	  out << "}" << std::endl;

  }

};



int bamExpRunner::connectFaceAwayMatesRegions(
  const njh::progutils::CmdArgs &inputCommands) {

  OutOptions outOpts(bfs::path(""), ".bed");
  bfs::path bedFnp;
  uint32_t initialSoftClipCutOff = 20;
  // uint32_t finalSoftClipCutOff = std::numeric_limits<uint32_t>::max();
  uint32_t finalSoftClipCutOff = 10;
  uint32_t softClipCutOffBothSides = 5;
  uint32_t minReadAmount = 2;
  uint32_t minInsertSize = 1000;
  uint32_t step = 50;
  uint32_t windowSize = 100;
  std::string sample = "sample";
  seqSetUp setUp(inputCommands);

  setUp.setOption(sample, "--sample", "sample name", true);
  setUp.setOption(initialSoftClipCutOff, "--initialSoftClipCutOff", "On initial pass if one side has this amount of soft clip don't count");
  setUp.setOption(finalSoftClipCutOff, "--finalSoftClipCutOff", "On final pass if one side has this amount of soft clip don't count");

  setUp.setOption(softClipCutOffBothSides, "--softClipCutOffBothSides", "If both sides are soft clipped more than this soft Clip Cut Off, then don't count it");
  setUp.setOption(step, "--step", "step");
  setUp.setOption(windowSize, "--windowSize", "windowSize");
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
  OutputStream finalOut(outOpts);
  finalOut << "#chrom\tstart\tend\tname\tsize\tstrand\tgroup\tpairCount\tsample\tminStart\tmaxEnd\tspanLength" << std::endl;

  std::vector<Bed3RecordCore> faceawayRegions;
  VecStr expectedRegular_status{
    "true::true::false::true",
    "true::false::true::false",
    "false::true::false::true",
    "false::false::true::false"};
  VecStr faceaway_status{
    "true::true::false::false",
    "true::false::true::true",
    "false::true::false::false",
    "false::false::true::true"};

  {
    BamTools::BamReader bReader;
    bReader.Open(setUp.pars_.ioOptions_.firstName_.string());
    checkBamOpenThrow(bReader, setUp.pars_.ioOptions_.firstName_.string());
    loadBamIndexThrow(bReader);
    auto refData = bReader.GetReferenceData();
    BamTools::BamAlignment bAln;
    // std::map<std::string, uint32_t> counts;
    for (const auto& reg: mergedRegions) {
      setBamFileRegionThrow(bReader, reg);
      std::vector<Bed6RecordCore> otherRegions;

      while (bReader.GetNextAlignment(bAln)) {
        //skip if either mate is not mapped or if the mates are mapped to different chromosomes
        if(!bAln.IsPrimaryAlignment() || !bAln.IsMapped() || !bAln.IsMateMapped() || bAln.RefID != bAln.MateRefID) {
          continue;
        }
        //soft clipping cut off
        if (bAln.CigarData.front().Type == 'S' && bAln.CigarData.front().Length > softClipCutOffBothSides &&
            bAln.CigarData.back().Type == 'S' && bAln.CigarData.back().Length > softClipCutOffBothSides) {
          continue;
        }
        //initialSoftClipCutOff clipping cut off
        if ((bAln.CigarData.front().Type == 'S' && bAln.CigarData.front().Length > initialSoftClipCutOff) ||
            (bAln.CigarData.back().Type == 'S' && bAln.CigarData.back().Length > initialSoftClipCutOff)) {
          continue;
        }


        auto mappingDirectionStatus = njh::pasteAsStr(njh::boolToStr(bAln.IsFirstMate()),
          "::", njh::boolToStr(bAln.IsReverseStrand()),
          "::", njh::boolToStr(bAln.IsMateReverseStrand()),
          "::", njh::boolToStr(bAln.Position > bAln.MatePosition));


        // std::cout << bAln.Name << "\t"
        //         << njh::colorBool(bAln.IsFirstMate())
        // << "\t" << njh::colorBool(bAln.IsReverseStrand())
        // << "\t" << njh::colorBool(bAln.IsMateReverseStrand())
        // << "\t" << njh::colorBool(bAln.Position > bAln.MatePosition) << std::endl;
        // ++counts[mappingDirectionStatus];

        if (njh::in(mappingDirectionStatus, faceaway_status) &&
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


    // table statuses(VecStr{"status", "count", "type"});
    // for(const auto & count : counts) {
    //   std::string type = "other";
    //   if(njh::in(count.first, expectedRegular_status)) {
    //     type = "expected";
    //   } else if (njh::in(count.first, faceaway_status)) {
    //     type = "faceaway";
    //   }
    //   auto addingRow = toVecStr(count.first, count.second, type);
    //   statuses.content_.emplace_back(addingRow);
    // }
    // statuses.outPutContentOrganized(std::cout);
  }

  // std::cout << "faceawayRegions: " << faceawayRegions.size() << std::endl;
  if(!faceawayRegions.empty()) {
    auto faceawayRegions_merged = BedUtility::mergeAndSort(faceawayRegions);
    // add the faceway regions, this will find windows that span outside of the original scanned window
    auto combinedMergedRegions = concatVecs(faceawayRegions_merged, mergedRegions);
    std::vector<Bed6RecordCore> combinedMergedRegions_beds;
    for(const auto & reg : combinedMergedRegions) {
      combinedMergedRegions_beds.emplace_back(reg.genBedRecordCore());
    }
    auto merged_combinedMergedRegions = BedUtility::mergeAndSort(combinedMergedRegions_beds);
    // for(const auto & reg : mergedRegions){
    for(const auto & reg : merged_combinedMergedRegions){
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
    // std::cout << "graph.nodes_.size(): " << graph.nodes_.size() << std::endl;
    // for(const auto & n : graph.nodes_ ) {
    //   std::cout << n.region_.genBedRecordCore().toDelimStr() << std::endl;
    // }
    // std::cout << std::endl;
    {
      BamTools::BamReader bReader;
      bReader.Open(setUp.pars_.ioOptions_.firstName_.string());
      checkBamOpenThrow(bReader, setUp.pars_.ioOptions_.firstName_.string());
      loadBamIndexThrow(bReader);
      auto refData = bReader.GetReferenceData();
      BamTools::BamAlignment bAln;
      // for (const auto& reg: mergedRegions) {
      for (const auto& reg: merged_combinedMergedRegions) {
        setBamFileRegionThrow(bReader, reg);
        std::vector<Bed6RecordCore> otherRegions;

        while (bReader.GetNextAlignment(bAln)) {
          //skip if either mate is not mapped or if the mates are mapped to different chromosomes
          if(!bAln.IsPrimaryAlignment() || !bAln.IsMapped() || !bAln.IsMateMapped() || bAln.RefID != bAln.MateRefID) {
            continue;
          }
          //soft clipping cut off
          if(bAln.CigarData.front().Type == 'S' && bAln.CigarData.front().Length > softClipCutOffBothSides &&
            bAln.CigarData.back().Type == 'S' && bAln.CigarData.back().Length > softClipCutOffBothSides) {
            continue;
          }
          //finalSoftClipCutOff clipping cut off
          if ((bAln.CigarData.front().Type == 'S' && bAln.CigarData.front().Length > finalSoftClipCutOff) ||
              (bAln.CigarData.back().Type == 'S' && bAln.CigarData.back().Length > finalSoftClipCutOff)) {
            continue;
          }
          auto mappingDirectionStatus = njh::pasteAsStr(njh::boolToStr(bAln.IsFirstMate()),
                                                        "::", njh::boolToStr(bAln.IsReverseStrand()),
                                                        "::", njh::boolToStr(bAln.IsMateReverseStrand()),
                                                        "::", njh::boolToStr(bAln.Position > bAln.MatePosition));

          if (njh::in(mappingDirectionStatus, faceaway_status) &&
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
            // std::cout << "baln_nodes: " << njh::conToStr(baln_nodes, ",") << std::endl;
            // std::cout << "mate_nodes: " << njh::conToStr(mate_nodes, ",") << std::endl;
            // std::cout << "genCigarStr: " << genCigarStr(bAln) << std::endl;
            //
            // std::cout << bAlnRegion.genBedRecordCore().toDelimStr() << std::endl;
            // std::cout << bAlnMateRegion.genBedRecordCore().toDelimStr() << std::endl;
            //
            //
            // std::cout << std::endl;
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
    graph.turnOffEdgesReadCountCutOff(minReadAmount);
    graph.determineGroups();
    // std::cout << "regions\tnode1Group\tnode2Group\tread_count\tsample_count" << std::endl;
    // for(const auto & e : graph.edges_) {
    //   std::cout << njh::conToStr(getVectorOfMapKeys(e->nodes_), ",")
    //     << "\t" << graph.nodes_[e->nodes_.begin()->second].group_
    //     << "\t" << graph.nodes_[e->nodes_[e->nodes_.begin()->first]].group_
    //     << "\t" << e->reads_.size()
    //     << "\t" << e->samples_.size() << std::endl;
    // }

    // OutOptions graphVizOut(njh::files::make_path("graph.dot"));
    // graphVizOut.overWriteFile_ = true;
    // graph.writeGraphViz(graphVizOut);

    std::map<uint32_t, std::vector<Bed3RecordCore>> regionsPerGroup;
    std::unordered_map<uint32_t, std::unordered_set<std::string>> pairsPerGroup;
    auto groupCounts = graph.getGroupCounts();
    for(const auto & n : graph.nodes_) {
      if(groupCounts[n.group_] > 1) {
        regionsPerGroup[n.group_].emplace_back(n.region_.genBed3RecordCore());
        for(const auto & e : n.edges_) {
          if(e->on_) {
            pairsPerGroup[n.group_].insert(e->reads_.begin(), e->reads_.end());
          }
        }
      }
    }

    for (auto &regionsForGroup: regionsPerGroup) {
      auto mergedRegionsForGroup = BedUtility::mergeAndSort(regionsForGroup.second);
      size_t minPos = std::numeric_limits<size_t>::max();
      size_t maxPos = 0;

      for(const auto & region : mergedRegionsForGroup) {
        if(region.start_ < minPos) {
          minPos = region.start_;
        }
        if(region.end_ > maxPos) {
          maxPos = region.end_;
        }
      }
      for (const auto &mergedRegion: mergedRegionsForGroup) {
        finalOut << mergedRegion.genBedRecordCore().toDelimStr()
        << "\t" << regionsForGroup.first
        << "\t" << pairsPerGroup[regionsForGroup.first].size()
        << "\t" << sample
        << "\t" << minPos
        << "\t" << maxPos
        << "\t" << maxPos - minPos
        << std::endl;
      }
    }
  }
  return 0;
}
} //  namespace njhseq {

