//
// Created by Nicholas Hathaway on 4/7/23.
//
#include "bamExp.hpp"
#include "elucidator/BamToolsUtils.h"
#include "elucidator/objects/BioDataObject.h"

namespace njhseq {

int bamExpRunner::getMateMapLocationForRegion(const njh::progutils::CmdArgs & inputCommands){
  OutOptions outOpts;

  bfs::path bedFnp = "";
  seqSetUp setUp(inputCommands);
  setUp.processVerbose();
  setUp.processDebug();
  setUp.processReadInNames({"--bam"}, true);

  setUp.setOption(bedFnp, "--bed", "regions to investigate", true);
  setUp.processWritingOptions(outOpts);
  setUp.finishSetUp(std::cout);
  std::vector<GenomicRegion> regions  = bedPtrsToGenomicRegs(getBeds(bedFnp));

  OutputStream out(outOpts);

  BamTools::BamAlignment bamAlignment;
  checkBamFilesForIndexesAndAbilityToOpen({setUp.pars_.ioOptions_.firstName_});
  BamTools::BamReader bReader;
  bReader.Open(setUp.pars_.ioOptions_.firstName_.string());
  checkBamOpenThrow(bReader, setUp.pars_.ioOptions_.firstName_);
  bReader.LocateIndex();
  auto refData = bReader.GetReferenceData();
  out << "chrom\tstart\tend\tname\tlen\tstrand\treadName\treadRef\treadPosition\treadEndPosition\treadRegionName\treadAlignLen\treadStrand\tmateMapped\tproperPair\tmateRef\tmatePosition\tmatePossibleEndPosition\tmateRegionName\tmateRegionLen\tmateStrand" << std::endl;
  for(const auto & region : regions){
    setBamFileRegionThrow(bReader, region);
    BamTools::BamAlignment bAln;
    while (bReader.GetNextAlignment(bAln)) {
      if (bAln.IsMapped() && bAln.IsPaired()) {
        out << region.genBedRecordCore().toDelimStr()
            << "\t" << bAln.Name
            << "\t" << refData[bAln.RefID].RefName
            << "\t" << bAln.Position
            << "\t" << getEndPosition(bAln)
            << "\t" << njh::pasteAsStr(refData[bAln.RefID].RefName,"-", bAln.Position, "-", getEndPosition(bAln))
            << "\t" << getEndPosition(bAln) - bAln.Position
            << "\t" << (bAln.IsReverseStrand() ? "-" : "+")
            << "\t" << njh::boolToStr(bAln.IsMateMapped())
            << "\t" << njh::boolToStr(bAln.IsProperPair());
        if (bAln.IsMateMapped()) {
          out << "\t" << refData[bAln.MateRefID].RefName
              << "\t" << bAln.MatePosition
              << "\t" << bAln.MatePosition + bAln.AlignedBases.size()
              << "\t" << njh::pasteAsStr(refData[bAln.MateRefID].RefName,"-", bAln.MatePosition, "-", bAln.MatePosition + bAln.AlignedBases.size())
              << "\t" << bAln.AlignedBases.size()
              << "\t" << (bAln.IsMateReverseStrand() ? "-" : "+") 
              << std::endl;
        } else {
          out << "\t" << "*"
              << "\t" << "*"
              << "\t" << "*"
              << "\t" << "*"
              << "\t" << "*"
              << "\t" << "*"
              << std::endl;
        }
      }
    }
  }

  return 0;
}

} //namespace njhseq

