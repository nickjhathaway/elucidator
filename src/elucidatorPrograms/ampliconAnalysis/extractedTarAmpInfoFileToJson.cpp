//
// Created by Nicholas Hathaway on 10/26/23.
//
#include <njhseq/IO/OutputStream.hpp>
#include <njhseq/objects/dataContainers/tables/TableReader.hpp>
#include <njhseq/objects/BioDataObject/GenomicRegion.hpp>
#include <njhseq/readVectorManipulation/readVectorHelpers/readVecSorter.hpp>
#include "ampliconAnalysisRunner.hpp"



namespace njhseq {



int ampliconAnalysisRunner::combingAllIntoPMOJson(const njh::progutils::CmdArgs &inputCommands) {
  OutOptions outOpts("", ".json");
  bfs::path bioinformatics_info_input_json_fnp;
  bfs::path specimen_infos_input_json_fnp;
  bfs::path experiment_infos_input_json_fnp;
  bfs::path sequencing_info_input_json_fnp;
  bfs::path panel_info_input_json_fnp;
  bfs::path target_demultiplexed_experiment_samples_json_fnp;
  bfs::path microhaplotypes_detected_and_representative_microhaplotype_sequences_input_json_fnp;

  std::string analysis_name;
  ampliconAnalysisSetUp setUp(inputCommands);
  setUp.setOption(analysis_name, "--analysis_name", "analysis name", true);

  setUp.setOption(bioinformatics_info_input_json_fnp, "--bioinformatics_info_input_json_fnp", "bioinformatics_info_input_json_fnp", true);
  setUp.setOption(specimen_infos_input_json_fnp, "--specimen_infos_input_json_fnp", "specimen_infos_input_json_fnp", true);
  setUp.setOption(experiment_infos_input_json_fnp, "--experiment_infos_input_json_fnp", "experiment_infos_input_json_fnp", true);

  setUp.setOption(sequencing_info_input_json_fnp, "--sequencing_info_input_json_fnp", "sequencing_info_input_json_fnp", true);
  setUp.setOption(panel_info_input_json_fnp, "--panel_info_input_json_fnp", "panel_info_input_json_fnp", true);
  setUp.setOption(target_demultiplexed_experiment_samples_json_fnp, "--target_demultiplexed_experiment_samples_json_fnp", "target_demultiplexed_experiment_samples_json_fnp");

  setUp.setOption(microhaplotypes_detected_and_representative_microhaplotype_sequences_input_json_fnp, "--microhaplotypes_detected_and_representative_microhaplotype_sequences_input_json_fnp", "microhaplotypes_detected_and_representative_microhaplotype_sequences_input_json_fnp", true);

  setUp.processWritingOptions(outOpts);
  setUp.finishSetUp(std::cout);

  OutputStream out(outOpts);

  Json::Value outJson;
  Json::Value bioinformatics_info_input_json = njh::json::parseFile(bioinformatics_info_input_json_fnp.string());
  Json::Value specimen_infos_input_json = njh::json::parseFile(specimen_infos_input_json_fnp.string());
  Json::Value experiment_infos_input_json = njh::json::parseFile(experiment_infos_input_json_fnp.string());

  Json::Value sequencing_info_input_json = njh::json::parseFile(sequencing_info_input_json_fnp.string());
  Json::Value panel_info_input_json = njh::json::parseFile(panel_info_input_json_fnp.string());
  Json::Value microhaplotypes_detected_and_representative_microhaplotype_sequences_input_json = njh::json::parseFile(microhaplotypes_detected_and_representative_microhaplotype_sequences_input_json_fnp.string());
  outJson["analysis_name"] = analysis_name;
  outJson["taramp_bioinformatics_infos"][bioinformatics_info_input_json.get(std::string("tar_amp_bioinformatics_info_id"), "").asString()] = bioinformatics_info_input_json;
  outJson["specimen_infos"] = specimen_infos_input_json;
  outJson["sequencing_infos"] = sequencing_info_input_json;
  outJson["experiment_infos"] = experiment_infos_input_json;
  outJson["panel_info"][panel_info_input_json.get(std::string("panel_id"), "").asString()] = panel_info_input_json;
  outJson["microhaplotypes_detected"][microhaplotypes_detected_and_representative_microhaplotype_sequences_input_json["microhaplotypes_detected"].get(std::string("tar_amp_bioinformatics_info_id"), "").asString()] = microhaplotypes_detected_and_representative_microhaplotype_sequences_input_json["microhaplotypes_detected"];
  outJson["representative_microhaplotype_sequences"][microhaplotypes_detected_and_representative_microhaplotype_sequences_input_json["representative_microhaplotype_sequences"].get(std::string("representative_microhaplotype_id"), "").asString()] = microhaplotypes_detected_and_representative_microhaplotype_sequences_input_json["representative_microhaplotype_sequences"];
  if(!target_demultiplexed_experiment_samples_json_fnp.empty()) {
    Json::Value target_demultiplexed_experiment_samples_json = njh::json::parseFile(target_demultiplexed_experiment_samples_json_fnp.string());
    outJson["target_demultiplexed_experiment_samples"][target_demultiplexed_experiment_samples_json.get(std::string("tar_amp_bioinformatics_info_id"), "").asString()] = target_demultiplexed_experiment_samples_json;
  }

  out << outJson << std::endl;
  return 0;
}

int ampliconAnalysisRunner::specimenInfoFileToJson(const njh::progutils::CmdArgs &inputCommands) {
  OutOptions outOpts("", ".json");
  bfs::path specimenInfoFnp;
  ampliconAnalysisSetUp setUp(inputCommands);
  setUp.setOption(specimenInfoFnp, "--specimenInfoFnp", "Name specimen Info Fnp", true);
  setUp.processWritingOptions(outOpts);
  setUp.finishSetUp(std::cout);

  OutputStream out(outOpts);

  Json::Value outJson;
  VecStr plateInfoCols{"plate_name", "plate_row", "plate_col"};
  {//
    TableReader reader(TableIOOpts::genTabFileIn(specimenInfoFnp));
    reader.header_.checkForColumnsThrow(toVecStr(VecStr{"specimen_id"},plateInfoCols ), __PRETTY_FUNCTION__ );

    VecStr row;
    //read in
    while (reader.getNextRow(row)) {
      std::string specimen_id = row[reader.header_.getColPos("specimen_id")];
      Json::Value & sampleJson = outJson[specimen_id];
      // Json::Value sampleJson;
      for(const auto & colName : reader.header_.columnNames_){
        if(colName == "specimen_id"){
          sampleJson[colName] = row[reader.header_.getColPos(colName)];
        } else {
          const auto & currentColValue = row[reader.header_.getColPos(colName)];
          if(isDoubleStr(currentColValue)){
            if(njh::strAllDigits(currentColValue)) {
              sampleJson[colName] = njh::json::toJson(njh::StrToNumConverter::stoToNum<uint32_t>(currentColValue));
            } else {
              sampleJson[colName] = njh::json::toJson(njh::StrToNumConverter::stoToNum<double>(currentColValue));
            }
          }else{
            sampleJson[colName] = currentColValue;
          }
        }
      }
      if(row[reader.header_.getColPos("plate_name")] == "NA") {
        sampleJson.removeMember("plate_name");
        sampleJson.removeMember("plate_row");
        sampleJson.removeMember("plate_col");
      }else {
        sampleJson["plate_name"] = row[reader.header_.getColPos("plate_name")];
        sampleJson["plate_row"] = row[reader.header_.getColPos("plate_row")];
        sampleJson["plate_col"] = njh::json::toJson(njh::StrToNumConverter::stoToNum<uint32_t>(row[reader.header_.getColPos("plate_col")]));
      }
      // Json::Value & plateJson = sampleJson["plate_info"];
      // plateJson["plate_name"] = row[reader.header_.getColPos("plate_name")];
      // plateJson["row"] = row[reader.header_.getColPos("row")];
      // plateJson["col"] = row[reader.header_.getColPos("col")];
      // outJson.append(sampleJson);
    }
  }

  out << outJson << std::endl;
  return 0;
}



int ampliconAnalysisRunner::demultiplexedExperimentSampleFileToJson(const njh::progutils::CmdArgs &inputCommands) {
  OutOptions outOpts("", ".json");
  bfs::path specimenInfoFnp;
  std::string experiment_sample_id_colName = "experiment_sample_id";
  std::string target_id_colName = "target_id";
  std::string raw_read_count_colName = "raw_read_count";
  std::string tar_amp_bioinformatics_info_id;

  ampliconAnalysisSetUp setUp(inputCommands);
  setUp.setOption(specimenInfoFnp, "--demultiplexedExperimentSampleInfoFnp", "Name demultiplexed Experiment Sample Fnp", true);
  setUp.setOption(tar_amp_bioinformatics_info_id, "--tar_amp_bioinformatics_info_id", "tar_amp_bioinformatics_info_id that is associated with this demultiplexed sample", true);

  setUp.setOption(experiment_sample_id_colName, "--experiment_sample_id_colName", "experiment_sample_id column name");
  setUp.setOption(target_id_colName, "--target_id_colName", "target_id column name");
  setUp.setOption(raw_read_count_colName, "--raw_read_count_colName", "raw_read_count column name");


  setUp.processWritingOptions(outOpts);
  setUp.finishSetUp(std::cout);

  OutputStream out(outOpts);

  Json::Value outJson;

  outJson["tar_amp_bioinformatics_info_id"] = tar_amp_bioinformatics_info_id;
  Json::Value & demultiplexed_experiment_samples_json = outJson["demultiplexed_experiment_samples"];

  {//
    VecStr requiredCols{experiment_sample_id_colName, target_id_colName, raw_read_count_colName};
    TableReader reader(TableIOOpts::genTabFileIn(specimenInfoFnp));
    reader.header_.checkForColumnsThrow(requiredCols, __PRETTY_FUNCTION__ );

    VecStr row;
    //read in
    while (reader.getNextRow(row)) {
      const std::string & experiment_sample_id = row[reader.header_.getColPos(experiment_sample_id_colName)];
      Json::Value & sampleJson = demultiplexed_experiment_samples_json[experiment_sample_id];
      if(njh::notIn(std::string("experiment_sample_id"), sampleJson.getMemberNames())) {
        sampleJson["experiment_sample_id"] = experiment_sample_id;
      }

      Json::Value & demultiplexed_targets_json = sampleJson["demultiplexed_targets"];

      const std::string & target_id = row[reader.header_.getColPos(target_id_colName)];
      const std::string & raw_read_count_str = row[reader.header_.getColPos(raw_read_count_colName)];
      if(njh::in(target_id, demultiplexed_targets_json.getMemberNames())) {
        std::stringstream ss;
        ss << __PRETTY_FUNCTION__ << ", error " << "already have " << target_id << " for " << experiment_sample_id << "\n";
        throw std::runtime_error{ss.str()};
      }
      Json::Value & target_id_json = demultiplexed_targets_json[target_id];

      if(!isDoubleStr(raw_read_count_str)) {
        std::stringstream ss;
        ss << __PRETTY_FUNCTION__ << ", error " << raw_read_count_colName << " should be a number, this does not look like a number: " <<  raw_read_count_str << "\n";
        throw std::runtime_error{ss.str()};
      }
      double raw_read_count = njh::StrToNumConverter::stoToNum<double>(raw_read_count_str);
      target_id_json["target_id"] = target_id;
      target_id_json["raw_read_count"] = raw_read_count;
    }
  }

  out << outJson << std::endl;
  return 0;
}







int ampliconAnalysisRunner::experimentInfoFileToJson(const njh::progutils::CmdArgs &inputCommands) {
  OutOptions outOpts("", ".json");
  bfs::path specimenInfoFnp;
  ampliconAnalysisSetUp setUp(inputCommands);
  setUp.setOption(specimenInfoFnp, "--experimentInfoFnp", "Name specimen Info Fnp", true);
  setUp.processWritingOptions(outOpts);
  setUp.finishSetUp(std::cout);

  OutputStream out(outOpts);

  Json::Value outJson;
  VecStr plateInfoCols{"plate_name", "plate_row", "plate_col"};
  {//
    TableReader reader(TableIOOpts::genTabFileIn(specimenInfoFnp));
    reader.header_.checkForColumnsThrow(toVecStr(VecStr{"experiment_sample_id"},plateInfoCols ), __PRETTY_FUNCTION__ );

    VecStr row;
    //read in
    while (reader.getNextRow(row)) {
      std::string experiment_sample_id = row[reader.header_.getColPos("experiment_sample_id")];
      Json::Value & sampleJson = outJson[experiment_sample_id];
      // Json::Value sampleJson;
      for(const auto & colName : reader.header_.columnNames_){
        if(colName == "experiment_sample_id"){
          sampleJson[colName] = row[reader.header_.getColPos(colName)];
        } else{
          const auto & currentColValue = row[reader.header_.getColPos(colName)];
          if(isDoubleStr(currentColValue) && colName != "specimen_id"){
            if(njh::strAllDigits(currentColValue)) {
              sampleJson[colName] = njh::json::toJson(njh::StrToNumConverter::stoToNum<uint32_t>(currentColValue));
            }else {
              sampleJson[colName] = njh::json::toJson(njh::StrToNumConverter::stoToNum<double>(currentColValue));
            }
          }else{
            sampleJson[colName] = currentColValue;
          }
        }
      }
      if(row[reader.header_.getColPos("plate_name")] == "NA") {
        sampleJson.removeMember("plate_name");
        sampleJson.removeMember("plate_row");
        sampleJson.removeMember("plate_col");
      }else {
        sampleJson["plate_name"] = row[reader.header_.getColPos("plate_name")];
        sampleJson["plate_row"] = row[reader.header_.getColPos("plate_row")];
        sampleJson["plate_col"] = njh::json::toJson(njh::StrToNumConverter::stoToNum<uint32_t>(row[reader.header_.getColPos("plate_col")]));
      }



      // Json::Value & plateJson = sampleJson["plate_info"];
      // plateJson["plate_name"] = row[reader.header_.getColPos("plate_name")];
      // plateJson["row"] = row[reader.header_.getColPos("row")];
      // plateJson["col"] = row[reader.header_.getColPos("col")];
      // outJson.append(sampleJson);
    }
  }

  out << outJson << std::endl;
  return 0;
}


int ampliconAnalysisRunner::finalClustersFileToJson(const njh::progutils::CmdArgs &inputCommands) {
  OutOptions outOpts("", ".json");
  std::string sequencing_id;
  std::string tar_amp_bioinformatics_info_id;
  std::string representative_microhaplotype_id;

  bfs::path finalClustersFnp;
  std::string sampleIDCol = "s_Sample";
  std::string targetIDCol = "p_name";
  std::string microhaplotypeIDCol = "h_popUID";
  std::string readCountCol = "c_ReadCnt";
  std::string umiCountCol = "c_barcodeCnt";

  std::string seqCol = "h_Consensus";

  ampliconAnalysisSetUp setUp(inputCommands);
  setUp.setOption(sequencing_id, "--sequencing_id", "sequencing id", true);
  setUp.setOption(tar_amp_bioinformatics_info_id, "--tar_amp_bioinformatics_info_id", "bioinformatics id", true);
  setUp.setOption(representative_microhaplotype_id, "--representative_microhaplotype_id", "representative microhaplotype id", true);

  setUp.setOption(finalClustersFnp, "--finalClustersFnp", "Name extracted Info Fnp", true);

  setUp.setOption(sampleIDCol, "--sampleIDCol", "sampleIDCol");
  setUp.setOption(targetIDCol, "--targetIDCol", "targetIDCol");
  setUp.setOption(microhaplotypeIDCol, "--microhaplotypeIDCol", "microhaplotypeIDCol");
  setUp.setOption(readCountCol, "--readCountCol", "readCountCol");
  setUp.setOption(seqCol, "--seqCol", "seqCol");

  setUp.setOption(umiCountCol, "--umiCountCol", "umiCountCol");

  setUp.processWritingOptions(outOpts);
  setUp.finishSetUp(std::cout);

  OutputStream out(outOpts);

  Json::Value outJson;
  Json::Value & microhaplotypes_detected = outJson["microhaplotypes_detected"];
  Json::Value & representative_microhaplotype_sequences = outJson["representative_microhaplotype_sequences"];

  microhaplotypes_detected["tar_amp_bioinformatics_info_id"] = tar_amp_bioinformatics_info_id;
  microhaplotypes_detected["representative_microhaplotype_id"] = representative_microhaplotype_id;

  representative_microhaplotype_sequences["representative_microhaplotype_id"] = representative_microhaplotype_id;

  Json::Value & samplesJson = microhaplotypes_detected["experiment_samples"];
  Json::Value & repTargetsJson = representative_microhaplotype_sequences["targets"];


  std::unordered_map<std::string, std::vector<std::shared_ptr<seqInfo>>> popSeqsByTarget;

  std::unordered_map<std::string, std::unordered_map<std::string,std::vector<std::shared_ptr<seqInfo>>>> popSeqsBySampleByTarget;

  bool hasUmi = false;
  {//
    TableReader reader(TableIOOpts::genTabFileIn(finalClustersFnp));
    hasUmi = njh::in(umiCountCol, reader.header_.columnNames_);
    VecStr row;
    //read in
    while (reader.getNextRow(row)) {
      bool found = false;
      std::string tarName = row[reader.header_.getColPos(targetIDCol)];
      std::string sampleName = row[reader.header_.getColPos(sampleIDCol)];
      std::shared_ptr<seqInfo> currentSeq = std::make_shared<seqInfo>(row[reader.header_.getColPos(microhaplotypeIDCol)],
                         row[reader.header_.getColPos(seqCol)]);
      currentSeq->cnt_ = njh::StrToNumConverter::stoToNum<double>(row[reader.header_.getColPos(readCountCol)]);
      if(hasUmi){
        currentSeq->frac_ = njh::StrToNumConverter::stoToNum<double>(row[reader.header_.getColPos(umiCountCol)]);;
      }
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
    //microhaplotypes_detected
    for(const auto & samp : popSeqsBySampleByTarget){
      Json::Value & sampJson = samplesJson[samp.first];
      sampJson["experiment_sample_id"] = samp.first;
      Json::Value & target_results = sampJson["target_results"];
      for(const auto & tar : samp.second){
        Json::Value & tarJson = target_results[tar.first];
        tarJson["target_id"] = tar.first;

        Json::Value & microhaplotypes = tarJson["microhaplotypes"];

        // VecStr ids;
        // std::vector<double> cnts;
        // std::vector<double> umiCnts;
        // for(const auto & s : tar.second){
        //   ids.emplace_back(s->name_);
        //   cnts.emplace_back(s->cnt_);
        //   if(hasUmi){
        //     umiCnts.emplace_back(s->frac_);
        //   }
        // }
        // tarJson["microhaplotype_ids"] = njh::json::toJson(ids);
        // tarJson["read_counts"] = njh::json::toJson(cnts);
        // if(hasUmi){
        //   tarJson["umi_counts"] = njh::json::toJson(umiCnts);
        // }

        for(const auto & s : tar.second){
          Json::Value microhaplotype;
          // microhaplotype["target_id"] = tar.first;
          microhaplotype["microhaplotype_id"] = s->name_;
          microhaplotype["read_count"] = s->cnt_;
          if(hasUmi){
            microhaplotype["umi_count"] = s->frac_;
          }
          microhaplotypes.append(microhaplotype);
        }
      }
    }
  }
  {
    //representative_microhaplotype_sequences
    for(auto & seqsForTar : popSeqsByTarget){
      //sort
      readVecSorter::sortByName(seqsForTar.second, true);
      Json::Value & seqsForTarJson = repTargetsJson[seqsForTar.first];
      seqsForTarJson["target_id"] = seqsForTar.first;
      Json::Value & seqs = seqsForTarJson["seqs"];

      for(const auto & seqForTar : seqsForTar.second){
        Json::Value seqForTarJson;
        seqForTarJson["microhaplotype_id"] = seqForTar->name_;
        seqForTarJson["seq"] = seqForTar->seq_;
        seqs[seqForTar->name_] = seqForTarJson;
        // seqs.append(seqForTarJson);
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

