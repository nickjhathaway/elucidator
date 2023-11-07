//
// Created by Nicholas Hathaway on 10/26/23.
//
#include <njhseq/IO/OutputStream.hpp>
#include <njhseq/objects/dataContainers/tables/TableReader.hpp>
#include <njhseq/objects/BioDataObject/GenomicRegion.hpp>
#include <njhseq/readVectorManipulation/readVectorHelpers/readVecSorter.hpp>
#include "ampliconAnalysisRunner.hpp"



namespace njhseq {



int ampliconAnalysisRunner::finalClustersFileToJson(const njh::progutils::CmdArgs &inputCommands) {
  OutOptions outOpts("", ".json");
  std::string sequencing_id;
  std::string bioinformatics_id;
  bfs::path finalClustersFnp;
  ampliconAnalysisSetUp setUp(inputCommands);
  setUp.setOption(sequencing_id, "--sequencing_id", "sequencing id", true);
  setUp.setOption(bioinformatics_id, "--bioinformatics_id", "bioinformatics id", true);

  setUp.setOption(finalClustersFnp, "--finalClustersFnp", "Name extracted Info Fnp", true);

  setUp.processWritingOptions(outOpts);
  setUp.finishSetUp(std::cout);

  OutputStream out(outOpts);

  Json::Value outJson;
  Json::Value & haplotypes_detected = outJson["haplotypes_detected"];
  Json::Value & representative_haplotype_sequences = outJson["representative_haplotype_sequences"];

  haplotypes_detected["sequencing_id"] = sequencing_id;
  haplotypes_detected["bioinformatics_id"] = bioinformatics_id;

  Json::Value & samplesJson = haplotypes_detected["samples"];


  std::unordered_map<std::string, std::vector<std::shared_ptr<seqInfo>>> popSeqsByTarget;

  std::unordered_map<std::string, std::unordered_map<std::string,std::vector<std::shared_ptr<seqInfo>>>> popSeqsBySampleByTarget;

  {//
    TableReader reader(TableIOOpts::genTabFileIn(finalClustersFnp));

    VecStr row;
    //read in
    while (reader.getNextRow(row)) {
      bool found = false;
      std::string tarName = row[reader.header_.getColPos("p_name")];
      std::string sampleName = row[reader.header_.getColPos("s_Name")];
      std::shared_ptr<seqInfo> currentSeq = std::make_shared<seqInfo>(row[reader.header_.getColPos("h_popUID")],
                         row[reader.header_.getColPos("h_Consensus")]);
      currentSeq->cnt_ = njh::StrToNumConverter::stoToNum<double>(row[reader.header_.getColPos("c_ReadCnt")]);
      for(const auto & previousSeq : popSeqsByTarget[tarName]){
        if(previousSeq->name_ == currentSeq->name_){
          found = true;
          break;
        }
      }
      popSeqsBySampleByTarget[sampleName][tarName].emplace_back(currentSeq);

      if (!found) {
        popSeqsByTarget[tarName].emplace_back(
            currentSeq
        );
      }
    }
  }
  {
    //haplotypes_detected
    for(const auto & samp : popSeqsBySampleByTarget){
      Json::Value & sampJson = samplesJson[samp.first];
      sampJson["sample_id"] = samp.first;
      Json::Value & target_results = sampJson["target_results"];
      for(const auto & tar : samp.second){
        Json::Value & tarJson = target_results[tar.first];
        tarJson["target_id"] = tar.first;
        VecStr ids;
        std::vector<double> cnts;
        for(const auto & s : tar.second){
          ids.emplace_back(s->name_);
          cnts.emplace_back(s->cnt_);
        }
        tarJson["haplotype_ids"] = njh::json::toJson(ids);
        tarJson["read_counts"] = njh::json::toJson(cnts);
      }
    }
  }
  {
    //representative_haplotype_sequences
    for(auto & seqsForTar : popSeqsByTarget){
      //sort
      readVecSorter::sortByName(seqsForTar.second, true);
      Json::Value & seqsForTarJson = representative_haplotype_sequences[seqsForTar.first];
      seqsForTarJson["target_id"] = seqsForTar.first;
      Json::Value & seqs = seqsForTarJson["seqs"];
      for(const auto & seqForTar : seqsForTar.second){
        Json::Value seqForTarJson;
        seqForTarJson["haplotype_id"] = seqForTar->name_;
        seqForTarJson["seq"] = seqForTar->seq_;
        seqs.append(seqForTarJson);
      }
    }
  }

  out << outJson << std::endl;

  return 0;
}

int ampliconAnalysisRunner::extractedTarAmpInfoFileToJson(const njh::progutils::CmdArgs &inputCommands) {
  bfs::path genomeInfoJsonFnp;
  OutOptions outOpts("", ".json");
  std::string panelName;
  bfs::path extractedInfoFnp;
  ampliconAnalysisSetUp setUp(inputCommands);
  setUp.setOption(panelName, "--panelName", "Name of the panel", true);
  setUp.setOption(extractedInfoFnp, "--extractedInfoFnp", "Name extracted Info Fnp", true);
  setUp.setOption(genomeInfoJsonFnp, "--genomeInfoJsonFnp", "genome Info Json Fnp", true);

  setUp.processWritingOptions(outOpts);
  setUp.finishSetUp(std::cout);

  auto genomeInfo = njh::json::parseFile(genomeInfoJsonFnp.string());

  OutputStream out(outOpts);

  Json::Value outJson;
  outJson["panel_id"] = panelName;
  outJson["target_genome"] = genomeInfo;


  {
    std::unordered_map<std::string, std::vector<VecStr>> extractedRowsPerID;
    TableReader reader(TableIOOpts::genTabFileIn(extractedInfoFnp));

    {
      VecStr row;
      //read in
      while(reader.getNextRow(row)){
        extractedRowsPerID[row[reader.header_.getColPos("target")]].emplace_back(row);
      }
    }

    for(const auto & extractedTargets : extractedRowsPerID){
      uint32_t targetCount = 0;
      std::unordered_set<std::string> insertLocs; //hack to get rid of multi-transcript genes which generate multiple rows per transcript

      for(const auto & row : extractedTargets.second){
        std::string targetName = extractedTargets.first;
        if(targetCount > 0){
          targetName += "." + njh::pasteAsStr(targetCount);
        }
        ++targetCount;
        GenomicRegion insert(row[reader.header_.getColPos("target")],
                             row[reader.header_.getColPos("#chrom")],
                             njh::StrToNumConverter::stoToNum<uint32_t>(row[reader.header_.getColPos("insertStart")]),
                             njh::StrToNumConverter::stoToNum<uint32_t>(row[reader.header_.getColPos("insertStop")]),
                             "-" == row[reader.header_.getColPos("strand")]);
        GenomicRegion fprimer(row[reader.header_.getColPos("target")] + "-forwardPrimer",
                              row[reader.header_.getColPos("#chrom")],
                              njh::StrToNumConverter::stoToNum<uint32_t>(row[reader.header_.getColPos("fPrimerStart")]),
                              njh::StrToNumConverter::stoToNum<uint32_t>(row[reader.header_.getColPos("fPrimerStop")]),
                              "-" == row[reader.header_.getColPos("strand")]);
        GenomicRegion rprimer(row[reader.header_.getColPos("target")] + "-reversePrimer",
                              row[reader.header_.getColPos("#chrom")],
                              njh::StrToNumConverter::stoToNum<uint32_t>(row[reader.header_.getColPos("rPrimerStart")]),
                              njh::StrToNumConverter::stoToNum<uint32_t>(row[reader.header_.getColPos("rPrimerStop")]),
                              "-" != row[reader.header_.getColPos("strand")]);
        if(!njh::in(njh::json::writeAsOneLine(insert.toJsonLocationOnly()), insertLocs)){
          //hack to get rid of multi-transcript genes which generate multiple rows per transcript
          insertLocs.emplace(njh::json::writeAsOneLine(insert.toJsonLocationOnly()));
          Json::Value tarInfo;
          tarInfo["target_id"] = targetName;
          tarInfo["insert_location"] = insert.toJsonLocationOnly();
          if (!row[reader.header_.getColPos("insertGeneDescription")].empty()) {
            tarInfo["gene_id"] = row[reader.header_.getColPos("insertGeneID")];
          }
          Json::Value forwardPrimers;
          Json::Value forwardPrimer;
          forwardPrimer["seq"] = row[reader.header_.getColPos("Fwd_primer")];
          forwardPrimer["location"] = fprimer.toJsonLocationOnly();
          forwardPrimers.append(forwardPrimer);
          tarInfo["forward_primers"] = forwardPrimers;

          Json::Value reversePrimers;
          Json::Value reversePrimer;
          reversePrimer["seq"] = row[reader.header_.getColPos("Rev_primer")];
          reversePrimer["location"] = fprimer.toJsonLocationOnly();
          reversePrimers.append(reversePrimer);
          tarInfo["reverse_primers"] = reversePrimers;
          outJson["targets"][targetName] = tarInfo;
        }
      }
    }


  }

  out << outJson << std::endl;

  return 0;

}


} // namespace njhseq

