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
  bfs::path specimen_experiment_infos_input_json_fnp;
  bfs::path sequencing_info_input_json_fnp;
  bfs::path panel_info_input_json_fnp;
  bfs::path reads_by_stage_json_fnp;
  bfs::path microhaplotypes_detected_and_representative_microhaplotype_sequences_input_json_fnp;

  std::string pmo_version = "v1.0.0";
  ampliconAnalysisSetUp setUp(inputCommands);
  setUp.setOption(pmo_version, "--pmo_version", "PMO version for this file");

  setUp.setOption(bioinformatics_info_input_json_fnp, "--bioinformatics_info_input_json_fnp", "bioinformatics_info_input_json_fnp", true);
  setUp.setOption(specimen_experiment_infos_input_json_fnp, "--specimen_experiment_infos_input_json_fnp", "specimen_experiment_infos_input_json_fnp", true);

  setUp.setOption(sequencing_info_input_json_fnp, "--sequencing_info_input_json_fnp", "sequencing_info_input_json_fnp", true);
  setUp.setOption(panel_info_input_json_fnp, "--panel_info_input_json_fnp", "panel_info_input_json_fnp", true);
  setUp.setOption(reads_by_stage_json_fnp, "--reads_by_stage_json_fnp", "reads_by_stage_json_fnp");

  setUp.setOption(microhaplotypes_detected_and_representative_microhaplotype_sequences_input_json_fnp, "--microhaplotypes_detected_and_representative_microhaplotype_sequences_input_json_fnp", "microhaplotypes_detected_and_representative_microhaplotype_sequences_input_json_fnp", true);

  setUp.processWritingOptions(outOpts);
  setUp.finishSetUp(std::cout);

  OutputStream out(outOpts);

  Json::Value outJson;
  outJson["pmo_header"]["pmo_version"] = pmo_version;
  auto today = njh::getCurrentDate();
  outJson["pmo_header"]["creation_date"] = today.substr(0, today.find('_'));
  outJson["pmo_header"]["generation_method"]["program_name"] = njh::pasteAsStr(setUp.commands_.masterProgram_, " ", setUp.commands_.subProgram_);
  outJson["pmo_header"]["generation_method"]["program_version"] = "1.1.1";

  Json::Value bioinformatics_info_input_json = njh::json::parseFile(bioinformatics_info_input_json_fnp.string());
  Json::Value specimen_experiment_infos_input_json = njh::json::parseFile(specimen_experiment_infos_input_json_fnp.string());
  Json::Value sequencing_info_input_json = njh::json::parseFile(sequencing_info_input_json_fnp.string());
  Json::Value panel_info_input_json = njh::json::parseFile(panel_info_input_json_fnp.string());
  Json::Value microhaplotypes_detected_and_representative_microhaplotype_sequences_input_json = njh::json::parseFile(microhaplotypes_detected_and_representative_microhaplotype_sequences_input_json_fnp.string());

  for (const auto & member : bioinformatics_info_input_json.getMemberNames()) {
    outJson[member] = bioinformatics_info_input_json[member];
  }
  for (const auto & member : specimen_experiment_infos_input_json.getMemberNames()) {
    outJson[member] = specimen_experiment_infos_input_json[member];
  }
  for (const auto & member : sequencing_info_input_json.getMemberNames()) {
    outJson[member] = sequencing_info_input_json[member];
  }
  for (const auto & member : panel_info_input_json.getMemberNames()) {
    outJson[member] = panel_info_input_json[member];
  }
  for (const auto & member : microhaplotypes_detected_and_representative_microhaplotype_sequences_input_json.getMemberNames()) {
    outJson[member] = microhaplotypes_detected_and_representative_microhaplotype_sequences_input_json[member];
  }

  if(!reads_by_stage_json_fnp.empty()) {
    Json::Value reads_by_stage_json = njh::json::parseFile(reads_by_stage_json_fnp.string());
    for (const auto & member : reads_by_stage_json.getMemberNames()) {
      outJson[member] = reads_by_stage_json[member];
    }
  }

  Json::StreamWriterBuilder builder;
  builder["indentation"] = "\t";  // or whatever you like
  std::unique_ptr<Json::StreamWriter> writer(
     builder.newStreamWriter());
  writer->write(outJson, &out);
  out << std::endl;
  return 0;
}


int ampliconAnalysisRunner::readsByStageToJson(const njh::progutils::CmdArgs &inputCommands) {
  OutOptions outOpts("", ".json");
  bfs::path readsByStageFnp;
  bfs::path rawCountsFnp;
  std::string experiment_sample_name_colName = "experiment_sample_name";
  std::string target_name_colName = "target_name";
  std::string read_count_colName = "read_count";
  std::string stage_colName = "stage";
  uint32_t bioinformatics_run_id;
  bfs::path panel_target_info_fnp;
  bfs::path experimental_info_fnp;
  ampliconAnalysisSetUp setUp(inputCommands);
  setUp.setOption(readsByStageFnp, "--readsByStageFnp", "reads By Stage Fnp", true);
  setUp.setOption(rawCountsFnp, "--rawCountsFnp", "raw Counts Fnp", true);
  setUp.setOption(bioinformatics_run_id, "--bioinformatics_run_id", "bioinformatics_run_id", true);
  setUp.setOption(panel_target_info_fnp, "--panel_target_info_fnp", "json file containing the information about the panel and target", true);
  setUp.setOption(experimental_info_fnp, "--experimental_info_fnp", "json file containing the information experiment_samples", true);
  setUp.setOption(experiment_sample_name_colName, "--experiment_sample_name_colName", "experiment_sample_name column name");
  setUp.setOption(target_name_colName, "--target_name_colName", "target_name column name");
  setUp.setOption(read_count_colName, "--read_count_colName", "read_count column name");
  setUp.setOption(stage_colName, "--stage_colName", "stage_ colum name");


  setUp.processWritingOptions(outOpts);
  setUp.finishSetUp(std::cout);


  auto panel_target_info = njh::json::parseFile(panel_target_info_fnp.string());
  if (!panel_target_info.isMember("target_info")) {
    std::stringstream ss;
    ss << __FILE__ << " " << __LINE__ << " " << __PRETTY_FUNCTION__ << ", error " << panel_target_info_fnp << " needs to have target_info, only has " << njh::conToStr(panel_target_info.getMemberNames(), ",") << "\n";
    throw std::runtime_error{ss.str()};
  }

  std::unordered_map<std::string, uint32_t> target_indexes;
  VecStr multiple_target_names;
  for (const auto & target_enum : iter::enumerate(panel_target_info["target_info"])) {
    auto target_name = target_enum.second["target_name"].asString();
    if (njh::in(target_name, target_indexes)) {
      multiple_target_names.emplace_back(target_name);
    }
    target_indexes[target_name] = target_enum.index;
  }
  if (!multiple_target_names.empty()) {
    std::stringstream ss;
    ss << __FILE__ << " " << __LINE__ << " " << __PRETTY_FUNCTION__ << ", error " << "found multiple of the same target in " << panel_target_info_fnp << ": " << njh::conToStr(multiple_target_names, ",") << "\n";
    throw std::runtime_error{ss.str()};
  }

  auto experimental_info = njh::json::parseFile(experimental_info_fnp.string());
  if (!experimental_info.isMember("experiment_info")) {
    std::stringstream ss;
    ss << __FILE__ << " " << __LINE__ << " " << __PRETTY_FUNCTION__ << ", error " << experimental_info_fnp << " needs to have experiment_info, only has " << njh::conToStr(experimental_info.getMemberNames(), ",") << "\n";
    throw std::runtime_error{ss.str()};
  }

  std::unordered_map<std::string, uint32_t> experiment_sample_indexes;
  VecStr multiple_experiment_sample_names;
  for (const auto & exp_samp_enum : iter::enumerate(experimental_info["experiment_info"])) {
    auto experiment_sample_name = exp_samp_enum.second["experiment_sample_name"].asString();
    if (njh::in(experiment_sample_name, experiment_sample_indexes)) {
      multiple_experiment_sample_names.emplace_back(experiment_sample_name);
    }
    experiment_sample_indexes[experiment_sample_name] = exp_samp_enum.index;
  }
  if (!multiple_experiment_sample_names.empty()) {
    std::stringstream ss;
    ss << __FILE__ << " " << __LINE__ << " " << __PRETTY_FUNCTION__ << ", error " << "found multiple of the same target in " << experimental_info_fnp << ": " << njh::conToStr(multiple_experiment_sample_names, ",") << "\n";
    throw std::runtime_error{ss.str()};
  }

  OutputStream out(outOpts);

  Json::Value outJson;
  Json::Value current_read_counts_by_stage;

  current_read_counts_by_stage["bioinformatics_run_id"] = bioinformatics_run_id;
  Json::Value & read_counts_by_experimental_sample_by_stage = current_read_counts_by_stage["read_counts_by_experimental_sample_by_stage"];

  std::unordered_map<std::string, uint32_t> experimental_sample_raw_read_counts;
  {
    VecStr requiredCols{experiment_sample_name_colName, read_count_colName};
    TableReader reader(TableIOOpts::genTabFileIn(rawCountsFnp));
    reader.header_.checkForColumnsThrow(requiredCols, __PRETTY_FUNCTION__ );
    VecStr row;
    //read in
    VecStr multiple_experimental_sample_names;
    VecStr missingSamples;

    while (reader.getNextRow(row)) {
      const auto & exp_samp = row[reader.header_.getColPos(experiment_sample_name_colName)];
      auto raw_read_count_str = row[reader.header_.getColPos(read_count_colName)];
      if(!isDoubleStr(raw_read_count_str)) {
        std::stringstream ss;
        ss << __FILE__ << " " << __LINE__ << " " << __PRETTY_FUNCTION__ << ", error " << read_count_colName << " should be a number, this does not look like a number: " <<  raw_read_count_str << "\n";
        throw std::runtime_error{ss.str()};
      }
      uint32_t read_counts = njh::StrToNumConverter::stoToNum<uint32_t>(raw_read_count_str);
      if (njh::in(exp_samp, experimental_sample_raw_read_counts)) {
        multiple_experimental_sample_names.emplace_back(exp_samp);
      }
      if (njh::notIn(exp_samp, experiment_sample_indexes)) {
        missingSamples.emplace_back(exp_samp);
      }
      experimental_sample_raw_read_counts[exp_samp] = read_counts;
    }
    if (!multiple_experimental_sample_names.empty()) {
      std::stringstream ss;
      ss << __FILE__ << " " << __LINE__ << " " << __PRETTY_FUNCTION__ << ", error " << "found multiple of " << njh::conToStr(multiple_experimental_sample_names, ",") << " in " << rawCountsFnp << "\n";
      throw std::runtime_error{ss.str()};
    }
    if (!missingSamples.empty()) {
      std::stringstream ss;
      ss << __FILE__ << " " << __LINE__ << " " << __PRETTY_FUNCTION__ << ", error " << "found " << njh::conToStr(missingSamples, ",") << " in " << rawCountsFnp << " but not in " << experimental_info_fnp << "\n";
      throw std::runtime_error{ss.str()};
    }
  }

  {//
    VecStr requiredCols{experiment_sample_name_colName, target_name_colName, read_count_colName, stage_colName};
    TableReader reader(TableIOOpts::genTabFileIn(readsByStageFnp));
    reader.header_.checkForColumnsThrow(requiredCols, __PRETTY_FUNCTION__ );
    VecStr row;
    //read in
    std::set<std::string> missingSamples;
    std::unordered_map<std::string, std::unordered_map<std::string, VecStr>> multiple_stage_for_experimental_sample_for_target;
    VecStr missing_targets;
    std::unordered_map<std::string, std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>>> exp_samp_stage_counts;
    while (reader.getNextRow(row)) {
      const std::string & experiment_sample_name = row[reader.header_.getColPos(experiment_sample_name_colName)];
      if (njh::notIn(experiment_sample_name, experimental_sample_raw_read_counts)) {
        missingSamples.emplace(experiment_sample_name);
      }

      const std::string & target_name = row[reader.header_.getColPos(target_name_colName)];
      const std::string & read_count_str = row[reader.header_.getColPos(read_count_colName)];
      const std::string & stage = row[reader.header_.getColPos(stage_colName)];
      if (njh::in(stage, exp_samp_stage_counts[experiment_sample_name][target_name])) {
        multiple_stage_for_experimental_sample_for_target[experiment_sample_name][target_name].emplace_back(stage);
      }
      if (njh::notIn(target_name, target_indexes)) {
        missing_targets.emplace_back(target_name);
      }
      if(!isDoubleStr(read_count_str)) {
        std::stringstream ss;
        ss << __FILE__ << " " << __LINE__ << " " << __PRETTY_FUNCTION__ << ", error " << read_count_colName << " should be a number, this does not look like a number: " <<  read_count_str << "\n";
        throw std::runtime_error{ss.str()};
      }
      uint32_t read_counts = njh::StrToNumConverter::stoToNum<uint32_t>(read_count_str);
      exp_samp_stage_counts[experiment_sample_name][target_name][stage] = read_counts;
    }

    if (!missingSamples.empty()) {
      std::stringstream ss;
      ss << __FILE__ << " " << __LINE__ << " " << __PRETTY_FUNCTION__ << ", error " << "found " << njh::conToStr(missingSamples, ",") << " in " << readsByStageFnp << " but not in " << rawCountsFnp << "\n";
      throw std::runtime_error{ss.str()};
    }
    if (!missing_targets.empty()) {
      std::stringstream ss;
      ss << __FILE__ << " " << __LINE__ << " " << __PRETTY_FUNCTION__ << ", error " << "found targets " << njh::conToStr(missing_targets, ",") << " in " << readsByStageFnp << " but not in " << panel_target_info_fnp << "\n";
      throw std::runtime_error{ss.str()};
    }
    if (!multiple_stage_for_experimental_sample_for_target.empty()) {
      std::stringstream ss;
      VecStr warnings;
      for (const auto & exp_name : multiple_stage_for_experimental_sample_for_target) {
        for (const auto & tar_name : exp_name.second) {
          warnings.emplace_back(njh::pasteAsStr(exp_name.first, " ", tar_name.first, ":" , njh::conToStr(tar_name.second, ",")));
        }
      }
      ss << __FILE__ << " " << __LINE__ << " " << __PRETTY_FUNCTION__ << ", error " << "found multiple stage counts for one or more experimental_sample_names  " << njh::conToStr(njh::getVecOfMapKeys(multiple_stage_for_experimental_sample_for_target), ",") << " in " << readsByStageFnp << "\n";
      ss << njh::conToStr(warnings, "\n") << "\n";
      throw std::runtime_error{ss.str()};
    }

    for (const auto & experimental_sample_name : experimental_sample_raw_read_counts) {
      Json::Value sampleJson;
      sampleJson["experiment_sample_id"] = experiment_sample_indexes[experimental_sample_name.first];
      sampleJson["total_raw_count"] = experimental_sample_name.second;
      if (njh::in(experimental_sample_name.first, exp_samp_stage_counts)) {
        auto & read_counts_for_targets = sampleJson["read_counts_for_targets"];
        for (const auto & tar_name : exp_samp_stage_counts[experimental_sample_name.first]) {
          Json::Value tar_json;
          tar_json["target_id"] = target_indexes[tar_name.first];
          auto & stages = tar_json["stages"];
          for (const auto & stage : tar_name.second) {
            Json::Value stageJson;
            stageJson["stage"] = stage.first;
            stageJson["read_count"] = stage.second;
            stages.append(stageJson);
          }
          read_counts_for_targets.append(tar_json);
        }
      }
      read_counts_by_experimental_sample_by_stage.append(sampleJson);
    }
  }
  outJson["read_counts_by_stage"].append(current_read_counts_by_stage);
  Json::StreamWriterBuilder builder;
  builder["indentation"] = "\t";  // or whatever you like
  std::unique_ptr<Json::StreamWriter> writer(
     builder.newStreamWriter());
  writer->write(outJson, &out);
  out << std::endl;
  return 0;
}

int ampliconAnalysisRunner::specimenExperimentInfoFileToJson(const njh::progutils::CmdArgs &inputCommands) {
  OutOptions outOpts("", ".json");
  bfs::path experimentInfoFnp;
  bfs::path specimenInfoFnp;
  uint32_t sequencing_info_id = 0;
  uint32_t panel_id = 0;
  ampliconAnalysisSetUp setUp(inputCommands);
  setUp.setOption(specimenInfoFnp, "--specimenInfoFnp", "Name specimen Info Fnp", true);
  setUp.setOption(experimentInfoFnp, "--experimentInfoFnp", "Name specimen Info Fnp", true);
  setUp.setOption(sequencing_info_id, "--sequencing_info_id", "sequencing_info_id", true);
  setUp.setOption(panel_id, "--panel_id", "panel_id", true);


  setUp.processWritingOptions(outOpts);
  setUp.finishSetUp(std::cout);

  OutputStream out(outOpts);

  Json::Value outJson;
  auto & specimen_info = outJson["specimen_info"];
  auto & experiment_info = outJson["experiment_info"];
  std::unordered_map<std::string, uint32_t> specimen_name_indexes;
  {
    VecStr plateInfoCols{"plate_name", "plate_row", "plate_col"};
    table reader(TableIOOpts::genTabFileIn(specimenInfoFnp));
    auto numeric_cols = reader.getNumericColumnPositions();
    VecStr specimenRequiredCols{"specimen_name", "samp_taxon_id", "collection_date", "collection_country", "collector", "samp_store_loc", "samp_collect_device", "project_name"};
    reader.checkForColumnsThrow(specimenRequiredCols, __PRETTY_FUNCTION__ );


    //read in
    uint32_t specimen_name_index = 0;
    VecStr multiple_specimen_names;
    for (const auto & row : reader) {
      std::string specimen_name = row[reader.getColPos("specimen_name")];
      if (njh::in(specimen_name, specimen_name_indexes)) {
        multiple_specimen_names.emplace_back(specimen_name);
      }
      specimen_name_indexes[specimen_name] = specimen_name_index;
      ++specimen_name_index;
      Json::Value sampleJson;
      for(const auto & colName : reader.columnNames_){
        if(colName == "specimen_name"){
          sampleJson[colName] = row[reader.getColPos(colName)];
        } else if (colName == "parasite_density_method" || colName == "parasite_density") {
          //do nothing
        } else {
          auto colPos = reader.getColPos(colName);
          const auto & currentColValue = row[colPos];
          if(njh::in(colPos, numeric_cols)){
            if(njh::strAllDigits(currentColValue)) {
              sampleJson[colName] = njh::json::toJson(njh::StrToNumConverter::stoToNum<uint32_t>(currentColValue));
            } else {
              sampleJson[colName] = njh::json::toJson(njh::StrToNumConverter::stoToNum<double>(currentColValue));
            }
          } else {
            sampleJson[colName] = currentColValue;
          }
        }
      }
      if (njh::in(std::string("parasite_density_method"), reader.columnNames_) && njh::in(std::string("parasite_density"), reader.columnNames_)) {
        if ("NA" != row[reader.getColPos("parasite_density")]) {
          Json::Value parasite_densityJson;
          parasite_densityJson["method"] = row[reader.getColPos("parasite_density_method")];
          parasite_densityJson["density"] = njh::json::toJson(njh::StrToNumConverter::stoToNum<double>(row[reader.getColPos("parasite_density")]));
          sampleJson["parasite_density_info"].append(parasite_densityJson);
        }
      }

      if(row[reader.getColPos("plate_name")] == "NA") {
        sampleJson.removeMember("plate_name");
        sampleJson.removeMember("plate_row");
        sampleJson.removeMember("plate_col");
      } else {
        sampleJson["plate_name"] = row[reader.getColPos("plate_name")];
        sampleJson["plate_row"] = row[reader.getColPos("plate_row")];
        sampleJson["plate_col"] = njh::json::toJson(njh::StrToNumConverter::stoToNum<uint32_t>(row[reader.getColPos("plate_col")]));
      }
      specimen_info.append(sampleJson);
    }
    if (!multiple_specimen_names.empty()) {
      std::stringstream ss;
      ss << __FILE__ << " " << __LINE__ << " " << __PRETTY_FUNCTION__ << ", error " << "found the following specimens multiple times: " << njh::conToStr(multiple_specimen_names, ",") << "\n";
      throw std::runtime_error{ss.str()};
    }
  }

  {
    table reader(TableIOOpts::genTabFileIn(experimentInfoFnp));
    reader.checkForColumnsThrow(toVecStr(VecStr{"experiment_sample_name", "specimen_name"}), __PRETTY_FUNCTION__ );
    auto numeric_cols = reader.getNumericColumnPositions();


    //read in
    VecStr experimental_sample_names;
    VecStr multiple_experimental_sample_names;
    for (const auto & row : reader) {
      std::string experiment_sample_name = row[reader.getColPos("experiment_sample_name")];
      if (njh::in(experiment_sample_name, experimental_sample_names)) {
        multiple_experimental_sample_names.emplace_back(experiment_sample_name);
      }
      experimental_sample_names.emplace_back(experiment_sample_name);
      Json::Value experimentSampleJson;
      experimentSampleJson["sequencing_info_id"] = sequencing_info_id;
      experimentSampleJson["panel_id"] = panel_id;
      for (const auto& colName: reader.columnNames_) {
        if (colName == "experiment_sample_name") {
          experimentSampleJson[colName] = row[reader.getColPos(colName)];
        } else if (colName == "specimen_name") {
          experimentSampleJson["specimen_id"] = specimen_name_indexes[row[reader.getColPos(colName)]];
        } else {
          auto colPos = reader.getColPos(colName);
          const auto& currentColValue = row[colPos];
          if(njh::in(colPos, numeric_cols)){
            if (njh::strAllDigits(currentColValue)) {
              experimentSampleJson[colName] = njh::json::toJson(
                njh::StrToNumConverter::stoToNum<uint32_t>(currentColValue));
            } else {
              experimentSampleJson[colName] = njh::json::toJson(
                njh::StrToNumConverter::stoToNum<double>(currentColValue));
            }
          } else {
            experimentSampleJson[colName] = currentColValue;
          }
        }
      }
      if(row[reader.getColPos("plate_name")] == "NA") {
        experimentSampleJson.removeMember("plate_name");
        experimentSampleJson.removeMember("plate_row");
        experimentSampleJson.removeMember("plate_col");
      } else {
        experimentSampleJson["plate_name"] = row[reader.getColPos("plate_name")];
        experimentSampleJson["plate_row"] = row[reader.getColPos("plate_row")];
        experimentSampleJson["plate_col"] = njh::json::toJson(njh::StrToNumConverter::stoToNum<uint32_t>(row[reader.getColPos("plate_col")]));
      }
      experiment_info.append(experimentSampleJson);
    }
    if (!multiple_experimental_sample_names.empty()) {
      std::stringstream ss;
      ss << __FILE__ << " " << __LINE__ << " " << __PRETTY_FUNCTION__ << ", error " << "found the following specimens multiple times: " << njh::conToStr(multiple_experimental_sample_names, ",") << "\n";
      throw std::runtime_error{ss.str()};
    }
  }

  Json::StreamWriterBuilder builder;
  builder["indentation"] = "\t";  // or whatever you like
  std::unique_ptr<Json::StreamWriter> writer(
     builder.newStreamWriter());
  writer->write(outJson, &out);
  out << std::endl;
  return 0;
}


int ampliconAnalysisRunner::finalClustersFileToJson(const njh::progutils::CmdArgs &inputCommands) {
  OutOptions outOpts("", ".json");
  // uint32_t sequencing_id;
  uint32_t bioinformatics_run_id;
  bfs::path panel_target_info_fnp;
  bfs::path experimental_info_fnp;

  bfs::path finalClustersFnp;
  std::string sampleIDCol = "s_Sample";
  std::string targetIDCol = "p_name";
  // std::string microhaplotypeIDCol = "h_popUID";
  std::string readCountCol = "c_ReadCnt";
  std::string umiCountCol = "c_barcodeCnt";

  std::string seqCol = "h_Consensus";

  ampliconAnalysisSetUp setUp(inputCommands);
  // setUp.setOption(sequencing_id, "--sequencing_id", "sequencing id", true);
  setUp.setOption(bioinformatics_run_id, "--bioinformatics_run_id", "bioinformatics_run_id", true);
  setUp.setOption(panel_target_info_fnp, "--panel_target_info_fnp", "json file containing the information about the panel and target", true);
  setUp.setOption(experimental_info_fnp, "--experimental_info_fnp", "json file containing the information experiment_samples", true);

  setUp.setOption(finalClustersFnp, "--finalClustersFnp", "Name extracted Info Fnp", true);

  setUp.setOption(sampleIDCol, "--sampleIDCol", "sampleIDCol");
  setUp.setOption(targetIDCol, "--targetIDCol", "targetIDCol");
  // setUp.setOption(microhaplotypeIDCol, "--microhaplotypeIDCol", "microhaplotypeIDCol");
  setUp.setOption(readCountCol, "--readCountCol", "readCountCol");
  setUp.setOption(seqCol, "--seqCol", "seqCol");

  setUp.setOption(umiCountCol, "--umiCountCol", "umiCountCol");

  setUp.processWritingOptions(outOpts);
  setUp.finishSetUp(std::cout);



  auto panel_target_info = njh::json::parseFile(panel_target_info_fnp.string());
  if (!panel_target_info.isMember("target_info")) {
    std::stringstream ss;
    ss << __FILE__ << " " << __LINE__ << " " << __PRETTY_FUNCTION__ << ", error " << panel_target_info_fnp << " needs to have target_info, only has " << njh::conToStr(panel_target_info.getMemberNames(), ",") << "\n";
    throw std::runtime_error{ss.str()};
  }

  std::unordered_map<std::string, uint32_t> target_indexes;
  VecStr multiple_target_names;
  for (const auto & target_enum : iter::enumerate(panel_target_info["target_info"])) {
    auto target_name = target_enum.second["target_name"].asString();
    if (njh::in(target_name, target_indexes)) {
      multiple_target_names.emplace_back(target_name);
    }
    target_indexes[target_name] = target_enum.index;
  }
  if (!multiple_target_names.empty()) {
    std::stringstream ss;
    ss << __FILE__ << " " << __LINE__ << " " << __PRETTY_FUNCTION__ << ", error " << "found multiple of the same target in " << panel_target_info_fnp << ": " << njh::conToStr(multiple_target_names, ",") << "\n";
    throw std::runtime_error{ss.str()};
  }

  auto experimental_info = njh::json::parseFile(experimental_info_fnp.string());
  if (!experimental_info.isMember("experiment_info")) {
    std::stringstream ss;
    ss << __FILE__ << " " << __LINE__ << " " << __PRETTY_FUNCTION__ << ", error " << experimental_info_fnp << " needs to have experiment_info, only has " << njh::conToStr(experimental_info.getMemberNames(), ",") << "\n";
    throw std::runtime_error{ss.str()};
  }

  std::unordered_map<std::string, uint32_t> experiment_sample_indexes;
  VecStr multiple_experiment_sample_names;
  for (const auto & exp_samp_enum : iter::enumerate(experimental_info["experiment_info"])) {
    auto experiment_sample_name = exp_samp_enum.second["experiment_sample_name"].asString();
    if (njh::in(experiment_sample_name, experiment_sample_indexes)) {
      multiple_experiment_sample_names.emplace_back(experiment_sample_name);
    }
    experiment_sample_indexes[experiment_sample_name] = exp_samp_enum.index;
  }
  if (!multiple_experiment_sample_names.empty()) {
    std::stringstream ss;
    ss << __FILE__ << " " << __LINE__ << " " << __PRETTY_FUNCTION__ << ", error " << "found multiple of the same target in " << experimental_info_fnp << ": " << njh::conToStr(multiple_experiment_sample_names, ",") << "\n";
    throw std::runtime_error{ss.str()};
  }
  OutputStream out(outOpts);

  Json::Value outJson;
  Json::Value & microhaplotypes_detected_full = outJson["microhaplotypes_detected"];
  Json::Value & microhaplotypes_info = outJson["microhaplotypes_info"];

  Json::Value microhaplotypes_detected;
  microhaplotypes_detected["bioinformatics_run_id"] = bioinformatics_run_id;

  Json::Value & samplesJson = microhaplotypes_detected["experiment_samples"];
  Json::Value & microhaplotypes_info_targets = microhaplotypes_info["targets"];


  std::map<std::string, std::vector<std::shared_ptr<seqInfo>>> popSeqsByTarget;

  //microhaplotypes_info
  {
    TableReader reader(TableIOOpts::genTabFileIn(finalClustersFnp));
    VecStr row;
    //read in
    while (reader.getNextRow(row)) {
      bool found = false;
      const std::string & tarName = row[reader.header_.getColPos(targetIDCol)];
      std::shared_ptr<seqInfo> currentSeq = std::make_shared<seqInfo>("", row[reader.header_.getColPos(seqCol)]);
      currentSeq->cnt_ = 1;
      for(const auto & previousSeq : popSeqsByTarget[tarName]){
        if(previousSeq->seq_ == currentSeq->seq_){
          found = true;
          ++previousSeq->cnt_;
          break;
        }
      }
      if (!found) {
        popSeqsByTarget[tarName].emplace_back(currentSeq);
      }
    }
  }
  VecStr missingTargetNames;
  for (auto & popSeqsForTarget : popSeqsByTarget) {
    if (njh::notIn(popSeqsForTarget.first, target_indexes)) {
      missingTargetNames.emplace_back(popSeqsForTarget.first);
    }
  }
  if (!missingTargetNames.empty()) {
    std::stringstream ss;
    ss << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << ", error " << "missing the following targets found in " << finalClustersFnp << " from " << panel_target_info_fnp << "\n";
    ss << njh::conToStr(missingTargetNames, "\n") << "\n";
    throw std::runtime_error{ss.str()};
  }

  std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>> popSeqsByTargetMhapIndex;
  std::unordered_map<std::string, uint32_t> popSeqsByTargetIndex;

  uint32_t microhaplotypes_target_index = 0;
  for (auto & popSeqsForTarget : popSeqsByTarget) {
    Json::Value microhaps_for_target_info;
    microhaps_for_target_info["target_id"] = target_indexes[popSeqsForTarget.first];
    readVecSorter::sortReadVectorSimple(popSeqsForTarget.second, "totalCount");
    uint32_t microhaplotypes_in_target_index = 0;
    for (const auto & seq : popSeqsForTarget.second) {
      Json::Value seqForTarJson;
      seqForTarJson["seq"] = seq->seq_;
      popSeqsByTargetMhapIndex[popSeqsForTarget.first][seq->seq_] = microhaplotypes_in_target_index;
      ++microhaplotypes_in_target_index;
      microhaps_for_target_info["microhaplotypes"].append(seqForTarJson);
    }
    popSeqsByTargetIndex[popSeqsForTarget.first] = microhaplotypes_target_index;
    ++microhaplotypes_target_index;
    microhaplotypes_info_targets.append(microhaps_for_target_info);
  }

  {
    //microhaplotypes_detected
    TableReader reader(TableIOOpts::genTabFileIn(finalClustersFnp));
    VecStr row;
    //read in
    bool has_umi = njh::in(umiCountCol, reader.header_.columnNames_);
    std::unordered_map<std::string, std::unordered_map<std::string, Json::Value>> detectedMhaps;
    while (reader.getNextRow(row)) {
      const std::string & tarName = row[reader.header_.getColPos(targetIDCol)];
      const std::string & sampleName = row[reader.header_.getColPos(sampleIDCol)];
      const std::string & seq = row[reader.header_.getColPos(seqCol)];
      Json::Value microhaplotype;
      microhaplotype["mhap_id"] = popSeqsByTargetMhapIndex[tarName][seq];
      microhaplotype["reads"] = njh::json::toJson(static_cast<uint32_t>(std::round(njh::StrToNumConverter::stoToNum<double>(row[reader.header_.getColPos(readCountCol)]))));
      if(has_umi){
        microhaplotype["umis"] = njh::json::toJson(static_cast<uint32_t>(std::round(njh::StrToNumConverter::stoToNum<double>(row[reader.header_.getColPos(umiCountCol)]))));
      }
      detectedMhaps[sampleName][tarName].append(microhaplotype);
    }
    for (const auto & samp : detectedMhaps) {
      Json::Value samp_info;
      samp_info["experiment_sample_id"] = experiment_sample_indexes[samp.first];
      for (const auto & tar : samp.second) {
        Json::Value tar_info;
        tar_info["mhaps_target_id"] = target_indexes[tar.first];
        tar_info["haps"] = tar.second;
        samp_info["target_results"].append(tar_info);
      }
      samplesJson.append(samp_info);
    }
  }

  microhaplotypes_detected_full.append(microhaplotypes_detected);
  Json::StreamWriterBuilder builder;
  builder["indentation"] = "\t";  // or whatever you like
  std::unique_ptr<Json::StreamWriter> writer(
     builder.newStreamWriter());
  writer->write(outJson, &out);
  out << std::endl;
  // out << outJson << std::endl;

  return 0;
}

int ampliconAnalysisRunner::extractedTarAmpInfoFileToJson(const njh::progutils::CmdArgs &inputCommands) {
  bfs::path genomeInfoJsonFnp;
  OutOptions outOpts("", ".json");
  std::string panelName;
  bfs::path reactionNameFnp;
  bfs::path extractedInfoFnp;
  bfs::path genomeTwoBit;
  bfs::path additionalTargetAttributes;
  std::string targetColName = "target";

  ampliconAnalysisSetUp setUp(inputCommands);
  setUp.setOption(reactionNameFnp, "--reactionNameFnp", "table with target column and reaction column", true);
  setUp.setOption(panelName, "--panelName", "Name of the panel", true);
  setUp.setOption(extractedInfoFnp, "--extractedInfoFnp", "Name extracted Info Fnp", true);
  setUp.setOption(genomeInfoJsonFnp, "--genomeInfoJsonFnp", "genome Info Json Fnp", true);
  setUp.setOption(genomeTwoBit, "--2bit", "genome 2bit file, if supplied will add ref_seq to panel info");
  setUp.setOption(targetColName, "--targetColName", "target Column Name");
  setUp.setOption(additionalTargetAttributes, "--additionalTargetAttributes", "additional Target Attributes", false);

  setUp.processWritingOptions(outOpts);
  setUp.finishSetUp(std::cout);

  auto genomeInfo = njh::json::parseFile(genomeInfoJsonFnp.string());

  OutputStream out(outOpts);

  Json::Value outJson;
  Json::Value & panel_info  = outJson["panel_info"];
  Json::Value & target_info = outJson["target_info"];
  Json::Value & targeted_genomes  = outJson["targeted_genomes"];
  targeted_genomes.append(genomeInfo);

  std::shared_ptr<TwoBit::TwoBitFile> treader;
  if (!genomeTwoBit.empty()) {
    treader = std::make_shared<TwoBit::TwoBitFile>(genomeTwoBit);
  }
  std::unordered_map<std::string, VecStr> additionalTargetAttributesMap;
  if (!additionalTargetAttributes.empty()) {
    table additionalTargetAttributesTab;
    additionalTargetAttributesTab = table(additionalTargetAttributes, "\t", true);
    additionalTargetAttributesTab.checkForColumnsThrow({targetColName}, __PRETTY_FUNCTION__);
    for (const auto & row : additionalTargetAttributesTab) {
      for (const auto & col : additionalTargetAttributesTab.columnNames_) {
        if (col != targetColName) {
          additionalTargetAttributesMap[row[additionalTargetAttributesTab.getColPos(targetColName)]].emplace_back(row[additionalTargetAttributesTab.getColPos(col)]);
        }
      }
    }
  }

  std::unordered_map<std::string, VecStr> reactionNameMap;
  table reactionNameTab;
  reactionNameTab = table(reactionNameFnp, "\t", true);
  reactionNameTab.checkForColumnsThrow({targetColName, "reaction"}, __PRETTY_FUNCTION__);
  for (const auto & row : reactionNameTab) {
    auto target_name = row[reactionNameTab.getColPos(targetColName)];
    auto reactions = tokenizeString(row[reactionNameTab.getColPos("reaction")], ",");
    if (njh::in(target_name, reactionNameMap)) {
      std::stringstream ss;
      ss << __FILE__ << " " << __LINE__ << " " << __PRETTY_FUNCTION__ << ", error " << "already have target " << target_name << " in " << reactionNameFnp << "\n";
      throw std::runtime_error{ss.str()};
    }
    reactionNameMap[target_name] = reactions;
  }

  std::unordered_map<std::string, uint32_t> targetIndex;

  {
    std::unordered_map<std::string, std::vector<VecStr>> extractedRowsPerID;
    TableReader reader(TableIOOpts::genTabFileIn(extractedInfoFnp));

    {
      VecStr row;
      //read in
      std::set<std::string> multipleTargets;
      while(reader.getNextRow(row)){
        auto target_name = row[reader.header_.getColPos(targetColName)];
        if (njh::in(target_name, extractedRowsPerID)) {
          multipleTargets.emplace(target_name);
        }
        extractedRowsPerID[target_name].emplace_back(row);
      }
      if (!multipleTargets.empty()) {
        std::stringstream ss;
        ss << __FILE__ << " " << __LINE__ << " " << __PRETTY_FUNCTION__ << ", error " << "already have " << njh::conToStr(multipleTargets, ",") << " in " << extractedInfoFnp << "\n";
        throw std::runtime_error{ss.str()};
      }
    }

    VecStr missingFromReactionTab;
    VecStr missingFromTargetTab;
    for(const auto & extractedTargets : extractedRowsPerID) {
      if (njh::notIn(extractedTargets.first, reactionNameMap)) {
        missingFromReactionTab.emplace_back(extractedTargets.first);
      }
    }
    for(const auto & reactionName : reactionNameMap) {
      if (njh::notIn(reactionName.first, extractedRowsPerID)) {
        missingFromTargetTab.emplace_back(reactionName.first);
      }
    }
    if (!missingFromTargetTab.empty() || !missingFromReactionTab.empty()) {
      VecStr warnings;
      if (!missingFromTargetTab.empty()) {
        warnings.emplace_back(njh::pasteAsStr("found in ", reactionNameFnp ," missing ", njh::conToStr(missingFromTargetTab, ","), " from ", extractedInfoFnp ));
      }
      if (!missingFromReactionTab.empty()) {
        warnings.emplace_back(njh::pasteAsStr("found in ", extractedInfoFnp ," missing ", njh::conToStr(missingFromReactionTab, ","), " from ", reactionNameFnp ));
      }
      std::stringstream ss;
      ss << __FILE__ << " " << __LINE__ << " " << __PRETTY_FUNCTION__ << ", error " << "\n";
      ss << njh::conToStr(warnings, "\n") << "\n";
      throw std::runtime_error{ss.str()};
    }
    uint32_t index = 0;
    for(const auto & extractedTargets : extractedRowsPerID) {
      uint32_t targetCount = 0;
      std::unordered_set<std::string> insertLocs; //hack to get rid of multi-transcript genes which generate multiple rows per transcript

      for(const auto & row : extractedTargets.second){
        std::string targetName = extractedTargets.first;
        if(targetCount > 0){
          targetName += "." + njh::pasteAsStr(targetCount);
        }
        ++targetCount;
        GenomicRegion insert(row[reader.header_.getColPos(targetColName)],
                             row[reader.header_.getColPos("#chrom")],
                             njh::StrToNumConverter::stoToNum<uint32_t>(row[reader.header_.getColPos("insertStart")]),
                             njh::StrToNumConverter::stoToNum<uint32_t>(row[reader.header_.getColPos("insertStop")]),
                             "-" == row[reader.header_.getColPos("strand")]);
        GenomicRegion fprimer(row[reader.header_.getColPos(targetColName)] + "-forwardPrimer",
                              row[reader.header_.getColPos("#chrom")],
                              njh::StrToNumConverter::stoToNum<uint32_t>(row[reader.header_.getColPos("fPrimerStart")]),
                              njh::StrToNumConverter::stoToNum<uint32_t>(row[reader.header_.getColPos("fPrimerStop")]),
                              "-" == row[reader.header_.getColPos("strand")]);
        GenomicRegion rprimer(row[reader.header_.getColPos(targetColName)] + "-reversePrimer",
                              row[reader.header_.getColPos("#chrom")],
                              njh::StrToNumConverter::stoToNum<uint32_t>(row[reader.header_.getColPos("rPrimerStart")]),
                              njh::StrToNumConverter::stoToNum<uint32_t>(row[reader.header_.getColPos("rPrimerStop")]),
                              "-" != row[reader.header_.getColPos("strand")]);
        if(!njh::in(njh::json::writeAsOneLine(insert.toJsonLocationOnly()), insertLocs)){
          //hack to get rid of multi-transcript genes which generate multiple rows per transcript
          insertLocs.emplace(njh::json::writeAsOneLine(insert.toJsonLocationOnly()));
          Json::Value tarInfo;
          tarInfo["target_name"] = targetName;
          if (!genomeTwoBit.empty()) {
            tarInfo["insert_location"] = insert.toJsonLocationOnly(*treader);
          } else {
            tarInfo["insert_location"] = insert.toJsonLocationOnly();
          }
          tarInfo["insert_location"]["genome_id"] = 0;
          if (!row[reader.header_.getColPos("insertGeneDescription")].empty()) {
            tarInfo["gene_name"] = row[reader.header_.getColPos("insertGeneID")];
          }
          Json::Value forwardPrimers;
          Json::Value forwardPrimer;
          forwardPrimer["seq"] = row[reader.header_.getColPos("Fwd_primer")];
          forwardPrimer["location"] = fprimer.toJsonLocationOnly();
          forwardPrimer["location"]["genome_id"] = 0;
          forwardPrimers.append(forwardPrimer);
          tarInfo["forward_primers"] = forwardPrimers;

          Json::Value reversePrimers;
          Json::Value reversePrimer;
          reversePrimer["seq"] = row[reader.header_.getColPos("Rev_primer")];
          reversePrimer["location"] = fprimer.toJsonLocationOnly();
          reversePrimer["location"]["genome_id"] = 0;
          reversePrimers.append(reversePrimer);
          tarInfo["reverse_primers"] = reversePrimers;
          targetIndex[targetName] = index;
          if (njh::in(targetName,  additionalTargetAttributesMap)) {
            tarInfo["target_attributes"] = njh::json::toJson(additionalTargetAttributesMap[targetName]);
          }
          ++index;
          target_info.append(tarInfo);
        }
      }
    }
  }


  Json::Value panel_info_current;
  panel_info_current["panel_name"] = panelName;
  auto & reactions = panel_info_current["reactions"];
  std::unordered_map<std::string, std::set<uint32_t>> targetIndexesForReactions;

  for (const auto & reaction : reactionNameMap) {
    for (const auto & reaction_name : reaction.second) {
      targetIndexesForReactions[reaction_name].emplace(targetIndex[reaction.first]);
    }
  }

  for (const auto & reaction : targetIndexesForReactions) {
    Json::Value reaction_info;
    reaction_info["reaction_name"] = reaction.first;
    reaction_info["panel_targets"] = njh::json::toJson(reaction.second);
    reactions.append(reaction_info);
  }

  panel_info.append(panel_info_current);

  Json::StreamWriterBuilder builder;
  builder["indentation"] = "\t";  // or whatever you like
  std::unique_ptr<Json::StreamWriter> writer(
     builder.newStreamWriter());
  writer->write(outJson, &out);
  out << std::endl;
  // out << outJson << std::endl;

  return 0;

}


} // namespace njhseq

