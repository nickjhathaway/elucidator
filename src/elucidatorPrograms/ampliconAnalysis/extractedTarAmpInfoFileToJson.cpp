//
// Created by Nicholas Hathaway on 10/26/23.
//

#include <njhseq/IO/OutputStream.hpp>
#include <njhseq/objects/dataContainers/tables/TableReader.hpp>
#include <njhseq/objects/BioDataObject/GenomicRegion.hpp>
#include <njhseq/objects/BioDataObject/reading.hpp>
#include <njhseq/readVectorManipulation/readVectorHelpers/readVecSorter.hpp>
#include "ampliconAnalysisRunner.hpp"



namespace njhseq {


int ampliconAnalysisRunner::sraMetaToJson(const njh::progutils::CmdArgs &inputCommands) {
  OutOptions outOpts("", ".json");
  uint32_t panel_id = 0;
  bfs::path sra_meta_fnp;
  std::string specimen_name_col = "sample_alias";
  bool make_run_accession_library_sample_name = false;
  ampliconAnalysisSetUp setUp(inputCommands);
  setUp.setOption(sra_meta_fnp, "--sra_meta_fnp", "sra_meta_fnp", true);
  setUp.setOption(panel_id, "--panel_id", "panel_id", true);
  setUp.setOption(specimen_name_col, "--specimen_name_col", "specimen_name_col");

  setUp.setOption(make_run_accession_library_sample_name, "--make_run_accession_library_sample_name", "make run accession library sample name");

  setUp.processWritingOptions(outOpts);
  outOpts.outFilename_ = njh::files::prependFileBasename(outOpts.outName(), "project_specimen_library_info_");

  setUp.finishSetUp(std::cout);

  OutputStream out_project_specimen_library(outOpts);
  outOpts.outFilename_ = njh::files::prependFileBasename(outOpts.outName(), "sequencing_info_");
  OutputStream out_seq_info(outOpts);


  Json::Value outJson_project_specimen_library;
  Json::Value outJson_seq_info;

  table sra_meta_table(sra_meta_fnp, "\t", true);
  sra_meta_table.checkForColumnsThrow(VecStr{
    "library_name","experiment_accession","fastq_ftp","run_accession",

    specimen_name_col,"tax_id","host_tax_id","collection_date","country","environment_biome",
    "environment_feature","environmental_medium","host_gravidity","age","host_sex","lat","lon",
    "sample_accession","sample_storage",

    "instrument_platform","instrument_model","library_layout","library_strategy","library_source",
    "library_selection",

    "study_title","study_alias","collected_by"}, __PRETTY_FUNCTION__);

  // project_info, ProjectInfo
  auto project_info_tab = sra_meta_table.getColumns(VecStr{"study_accession", "study_title","study_alias","collected_by"}).getUniqueRows();
  Json::Value project_info_jsons;



  auto split_project_info_tab = project_info_tab.splitTableOnColumn("study_accession");
  std::unordered_map<std::string, uint32_t> project_indexes;
  uint32_t project_index = 0;
  for (const auto & project : split_project_info_tab) {
    Json::Value project_info_json;
    project_info_json["project_name"] = project.first;
    project_info_json["project_description"] = njh::conToStr(getUniqueStrings(project.second.getColumn("study_title")), ";");
    project_info_json["BioProject_accession"] = project.first;
    auto project_contributors = getUniqueStrings(project.second.getColumn("collected_by"));
    if (!std::all_of(project_contributors.begin(), project_contributors.end(), [](const std::string & project_contributor) {
      return project_contributor.empty() || njh::allWhiteSpaceStr(project_contributor);
    })) {
      project_info_json["project_contributors"] = njh::json::toJson(project_contributors);
    }
    project_info_jsons.append(project_info_json);
    project_indexes[project.first] = project_index;
    ++project_index;
  }

  outJson_project_specimen_library["project_info"] = project_info_jsons;
  std::unordered_map<std::string, uint32_t> seq_info_key_to_id;

  {
    //sequencing_info, SequencingInfo
    auto sequencing_info_tab = sra_meta_table.getColumns(VecStr{
      "instrument_platform","instrument_model",
      "library_layout","library_strategy","library_source",
      "library_selection"}).getUniqueRows();
    Json::Value sequencing_info_jsons;
    uint32_t seq_info_count = 0;
    for (const auto & row : sequencing_info_tab) {
      Json::Value sequencing_info_json;
      std::string sequencing_info_name = "seq_info";
      if (seq_info_count > 0) {
        sequencing_info_name += std::string("_") + estd::to_string(seq_info_count);
      }
      auto seq_info_key = njh::conToStr(VecStr{
        row[sequencing_info_tab.getColPos("instrument_platform")],
        row[sequencing_info_tab.getColPos("instrument_model")],
        row[sequencing_info_tab.getColPos("library_layout")],
        row[sequencing_info_tab.getColPos("library_strategy")],
        row[sequencing_info_tab.getColPos("library_source")],
        row[sequencing_info_tab.getColPos("library_selection")]
      },
        "-");
      sequencing_info_json["sequencing_info_name"] = sequencing_info_name;
      sequencing_info_json["seq_platform"] = row[sequencing_info_tab.getColPos("instrument_platform")];
      sequencing_info_json["seq_instrument_model"] = row[sequencing_info_tab.getColPos("instrument_model")];
      sequencing_info_json["library_layout"] = row[sequencing_info_tab.getColPos("library_layout")];
      sequencing_info_json["library_strategy"] = row[sequencing_info_tab.getColPos("library_strategy")];
      sequencing_info_json["library_source"] = row[sequencing_info_tab.getColPos("library_source")];
      sequencing_info_json["library_selection"] = row[sequencing_info_tab.getColPos("library_selection")];
      seq_info_key_to_id[seq_info_key] = seq_info_count;
      ++seq_info_count;
      sequencing_info_jsons.append(sequencing_info_json);
    }
    outJson_seq_info["sequencing_info"] = sequencing_info_jsons;
  }
  //SpecimenInfo, specimen_info
  std::unordered_map<std::string, uint32_t> specimen_info_name_to_id;

  {

    auto specimen_info_tab = sra_meta_table.getColumns(VecStr{
      specimen_name_col,"tax_id","host_tax_id","collection_date","country",
      "environment_biome", "environment_feature","environmental_medium",
      "host_gravidity","age",
      "host_sex",
      "lat","lon",
    "sample_accession","sample_storage",
      "study_accession"
    }).getUniqueRows();

    std::unordered_map<std::string, bool> column_testing;
    VecStr na_values_to_avoid{"NA", "missing", "na", "n/a", ""};
    column_testing["environment_biome"] = specimen_info_tab.column_all_na("environment_biome", na_values_to_avoid);
    column_testing["environment_feature"] = specimen_info_tab.column_all_na("environment_feature", na_values_to_avoid);
    column_testing["environmental_medium"] = specimen_info_tab.column_all_na("environmental_medium", na_values_to_avoid);
    column_testing["host_gravidity"] = specimen_info_tab.column_all_na("host_gravidity", na_values_to_avoid);
    column_testing["age"] = specimen_info_tab.column_all_na("age", na_values_to_avoid);
    column_testing["host_sex"] = specimen_info_tab.column_all_na("host_sex", na_values_to_avoid);
    column_testing["lat"] = specimen_info_tab.column_all_na("lat", na_values_to_avoid);
    column_testing["sample_accession"] = specimen_info_tab.column_all_na("sample_accession", na_values_to_avoid);
    column_testing["sample_storage"] = specimen_info_tab.column_all_na("sample_storage", na_values_to_avoid);


    Json::Value specimen_info_jsons;
    uint32_t specimen_info_count = 0;

    for (const auto & row : specimen_info_tab) {
      Json::Value specimen_info_json;

      specimen_info_json["specimen_name"] = row[specimen_info_tab.getColPos(specimen_name_col)];
      specimen_info_json["specimen_taxon_id"].append(njh::json::toJson(njh::StrToNumConverter::stoToNum<uint32_t>(row[specimen_info_tab.getColPos("tax_id")])));
      if ("NA" == row[specimen_info_tab.getColPos("host_tax_id")] || njh::allWhiteSpaceStr(row[specimen_info_tab.getColPos("host_tax_id")])) {
        specimen_info_json["host_taxon_id"] = njh::json::toJson(std::numeric_limits<uint32_t>::max());
      } else {
        specimen_info_json["host_taxon_id"] = njh::json::toJson(njh::StrToNumConverter::stoToNum<uint32_t>(row[specimen_info_tab.getColPos("host_tax_id")]) );
      }
      auto collection_date = row[specimen_info_tab.getColPos("collection_date")];
      VecStr na_values{"n/a", "NA", "missing", "", " "};
      if (njh::in(collection_date, na_values) || njh::allWhiteSpaceStr(collection_date)) {
        collection_date = "NA";
      }

      if (std::string::npos != collection_date.find('/')) {
        collection_date = collection_date.substr(0, collection_date.find('/'));
      }
      specimen_info_json["collection_date"] = collection_date;
      specimen_info_json["collection_country"] = row[specimen_info_tab.getColPos("country")];

      if (column_testing["environment_biome"]) specimen_info_json["env_broad_scale"] = row[specimen_info_tab.getColPos("environment_biome")];
      if (column_testing["environment_feature"]) specimen_info_json["env_local_scale"] = row[specimen_info_tab.getColPos("environment_feature")];
      if (column_testing["environmental_medium"]) specimen_info_json["env_medium"] = row[specimen_info_tab.getColPos("environmental_medium")];
      if (row[specimen_info_tab.getColPos("country")].find(':') != std::string::npos &&
        row[specimen_info_tab.getColPos("country")].find(':') + 1 != row[specimen_info_tab.getColPos("country")].size()) {
        auto geo_admin1 = row[specimen_info_tab.getColPos("country")].substr(row[specimen_info_tab.getColPos("country")].find(':') + 1);
        njh::lstrip(geo_admin1, ' ');
        njh::rstrip(geo_admin1, ' ');
        specimen_info_json["geo_admin1"] = geo_admin1;
        specimen_info_json["collection_country"] = row[specimen_info_tab.getColPos("country")].substr(0, row[specimen_info_tab.getColPos("country")].find(':'));
      }

      if (column_testing["host_gravidity"]) specimen_info_json["gravidity"] = row[specimen_info_tab.getColPos("host_gravidity")];
      if (column_testing["host_age"]) specimen_info_json["host_age"] = row[specimen_info_tab.getColPos("host_age")];
      if (column_testing["host_sex"]) specimen_info_json["host_sex"] = row[specimen_info_tab.getColPos("host_sex")];
      if (column_testing["sample_accession"]) specimen_info_json["specimen_accession"] = row[specimen_info_tab.getColPos("sample_accession")];
      if (column_testing["sample_storage"]) specimen_info_json["specimen_store_loc"] = row[specimen_info_tab.getColPos("sample_storage")];
      if (column_testing["lat"])
        specimen_info_json["lat_lon"] = njh::pasteAsStr(row[specimen_info_tab.getColPos("lat")],
                                                        ",",
                                                        row[specimen_info_tab.getColPos("lon")]);
      specimen_info_json["project_id"] = njh::json::toJson(project_indexes[row[specimen_info_tab.getColPos("study_accession")]]);
      specimen_info_name_to_id[row[specimen_info_tab.getColPos(specimen_name_col)]] = specimen_info_count;
      ++specimen_info_count;
      specimen_info_jsons.append(specimen_info_json);
    }
    outJson_project_specimen_library["specimen_info"] = specimen_info_jsons;
  }
  // LibrarySampleInfo, library_sample_info
  {
    auto library_sample_info_tab = sra_meta_table.getColumns(VecStr{
      "library_name",
      "experiment_accession","fastq_ftp","run_accession",

        specimen_name_col ,
      "instrument_platform","instrument_model","library_layout","library_strategy","library_source",
      "library_selection",
    }).getUniqueRows();
    Json::Value library_sample_info_jsons;

    for (const auto &row: library_sample_info_tab) {
      Json::Value library_sample_info_json;
      auto seq_info_key = njh::conToStr(VecStr{
                                          row[library_sample_info_tab.getColPos("instrument_platform")],
                                          row[library_sample_info_tab.getColPos("instrument_model")],
                                          row[library_sample_info_tab.getColPos("library_layout")],
                                          row[library_sample_info_tab.getColPos("library_strategy")],
                                          row[library_sample_info_tab.getColPos("library_source")],
                                          row[library_sample_info_tab.getColPos("library_selection")]
                                        },
                                        "-");
      auto seq_info_index = seq_info_key_to_id[seq_info_key];
      auto specimen_id = specimen_info_name_to_id[row[library_sample_info_tab.getColPos(specimen_name_col)]];
      library_sample_info_json["specimen_id"] = njh::json::toJson(specimen_id);
      library_sample_info_json["sequencing_info_id"] = njh::json::toJson(seq_info_index);
      library_sample_info_json["panel_id"] = njh::json::toJson(panel_id);
      library_sample_info_json["experiment_accession"] = row[library_sample_info_tab.getColPos("experiment_accession")];
      library_sample_info_json["fastqs_loc"] = row[library_sample_info_tab.getColPos("fastq_ftp")];
      library_sample_info_json["run_accession"] = row[library_sample_info_tab.getColPos("run_accession")];
      if (make_run_accession_library_sample_name) {
        library_sample_info_json["library_sample_name"] = row[library_sample_info_tab.getColPos("run_accession")];
        library_sample_info_json["alternate_identifiers"].append( row[library_sample_info_tab.getColPos("library_name")]);
      } else {
        library_sample_info_json["library_sample_name"] = row[library_sample_info_tab.getColPos("library_name")];
      }
      library_sample_info_jsons.append(library_sample_info_json);
    }
    outJson_project_specimen_library["library_sample_info"] = library_sample_info_jsons;
  }
  {
    Json::StreamWriterBuilder builder;
    builder["indentation"] = "\t";  // or whatever you like
    std::unique_ptr<Json::StreamWriter> writer(
       builder.newStreamWriter());
    writer->write(outJson_project_specimen_library, &out_project_specimen_library);
    out_project_specimen_library << std::endl;
  }
  {
    Json::StreamWriterBuilder builder;
    builder["indentation"] = "\t";  // or whatever you like
    std::unique_ptr<Json::StreamWriter> writer(
       builder.newStreamWriter());
    writer->write(outJson_seq_info, &out_seq_info);
    out_seq_info << std::endl;
  }
  return 0;
}

int ampliconAnalysisRunner::combingAllIntoPMOJson(const njh::progutils::CmdArgs &inputCommands) {
  OutOptions outOpts("", ".json");
  bfs::path bioinformatics_info_input_json_fnp;
  bfs::path specimen_library_sample_infos_input_json_fnp;
  bfs::path sequencing_info_input_json_fnp;
  bfs::path panel_info_input_json_fnp;
  bfs::path reads_by_stage_json_fnp;
  bfs::path detected_microhaplotypes_and_representative_microhaplotype_sequences_input_json_fnp;

  std::string pmo_version = "v1.0.0";
  ampliconAnalysisSetUp setUp(inputCommands);
  setUp.setOption(pmo_version, "--pmo_version", "PMO version for this file");

  setUp.setOption(bioinformatics_info_input_json_fnp, "--bioinformatics_info_input_json_fnp", "bioinformatics_info_input_json_fnp", true);
  setUp.setOption(specimen_library_sample_infos_input_json_fnp, "--specimen_library_sample_infos_input_json_fnp", "specimen_library_sample_infos_input_json_fnp", true);

  setUp.setOption(sequencing_info_input_json_fnp, "--sequencing_info_input_json_fnp", "sequencing_info_input_json_fnp", true);
  setUp.setOption(panel_info_input_json_fnp, "--panel_info_input_json_fnp", "panel_info_input_json_fnp", true);
  setUp.setOption(reads_by_stage_json_fnp, "--reads_by_stage_json_fnp", "reads_by_stage_json_fnp");

  setUp.setOption(detected_microhaplotypes_and_representative_microhaplotype_sequences_input_json_fnp, "--detected_microhaplotypes_and_representative_microhaplotype_sequences_input_json_fnp", "detected_microhaplotypes_and_representative_microhaplotype_sequences_input_json_fnp", true);

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
  Json::Value specimen_library_sample_infos_input_json = njh::json::parseFile(specimen_library_sample_infos_input_json_fnp.string());
  Json::Value sequencing_info_input_json = njh::json::parseFile(sequencing_info_input_json_fnp.string());
  Json::Value panel_info_input_json = njh::json::parseFile(panel_info_input_json_fnp.string());
  Json::Value detected_microhaplotypes_and_representative_microhaplotype_sequences_input_json = njh::json::parseFile(detected_microhaplotypes_and_representative_microhaplotype_sequences_input_json_fnp.string());

  for (const auto & member : bioinformatics_info_input_json.getMemberNames()) {
    outJson[member] = bioinformatics_info_input_json[member];
  }
  for (const auto & member : specimen_library_sample_infos_input_json.getMemberNames()) {
    outJson[member] = specimen_library_sample_infos_input_json[member];
  }
  for (const auto & member : sequencing_info_input_json.getMemberNames()) {
    outJson[member] = sequencing_info_input_json[member];
  }
  for (const auto & member : panel_info_input_json.getMemberNames()) {
    outJson[member] = panel_info_input_json[member];
  }
  for (const auto & member : detected_microhaplotypes_and_representative_microhaplotype_sequences_input_json.getMemberNames()) {
    outJson[member] = detected_microhaplotypes_and_representative_microhaplotype_sequences_input_json[member];
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
  std::string library_sample_name_colName = "library_sample_name";
  std::string target_name_colName = "target_name";
  std::string read_count_colName = "read_count";
  std::string stage_colName = "stage";
  uint32_t bioinformatics_run_id;
  bfs::path panel_target_info_fnp;
  bfs::path library_sample_info_fnp;
  ampliconAnalysisSetUp setUp(inputCommands);
  setUp.setOption(readsByStageFnp, "--readsByStageFnp", "reads By Stage Fnp", true);
  setUp.setOption(rawCountsFnp, "--rawCountsFnp", "raw Counts Fnp", true);
  setUp.setOption(bioinformatics_run_id, "--bioinformatics_run_id", "bioinformatics_run_id", true);
  setUp.setOption(panel_target_info_fnp, "--panel_target_info_fnp", "json file containing the information about the panel and target", true);
  setUp.setOption(library_sample_info_fnp, "--library_sample_info_fnp", "json file containing the information library_samples", true);
  setUp.setOption(library_sample_name_colName, "--library_sample_name_colName", "library_sample_name column name");
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

  auto library_sample_info = njh::json::parseFile(library_sample_info_fnp.string());
  if (!library_sample_info.isMember("library_sample_info")) {
    std::stringstream ss;
    ss << __FILE__ << " " << __LINE__ << " " << __PRETTY_FUNCTION__ << ", error " << library_sample_info_fnp << " needs to have library_sample_info, only has " << njh::conToStr(library_sample_info.getMemberNames(), ",") << "\n";
    throw std::runtime_error{ss.str()};
  }

  std::unordered_map<std::string, uint32_t> library_sample_indexes;
  VecStr multiple_library_sample_names;
  for (const auto & exp_samp_enum : iter::enumerate(library_sample_info["library_sample_info"])) {
    auto library_sample_name = exp_samp_enum.second["library_sample_name"].asString();
    if (njh::in(library_sample_name, library_sample_indexes)) {
      multiple_library_sample_names.emplace_back(library_sample_name);
    }
    library_sample_indexes[library_sample_name] = exp_samp_enum.index;
  }
  if (!multiple_library_sample_names.empty()) {
    std::stringstream ss;
    ss << __FILE__ << " " << __LINE__ << " " << __PRETTY_FUNCTION__ << ", error " << "found multiple of the same target in " << library_sample_info_fnp << ": " << njh::conToStr(multiple_library_sample_names, ",") << "\n";
    throw std::runtime_error{ss.str()};
  }

  OutputStream out(outOpts);

  Json::Value outJson;
  Json::Value current_read_counts_by_stage;

  current_read_counts_by_stage["bioinformatics_run_id"] = bioinformatics_run_id;
  Json::Value & read_counts_by_library_sample_by_stage = current_read_counts_by_stage["read_counts_by_library_sample_by_stage"];

  std::unordered_map<std::string, uint32_t> library_sample_raw_read_counts;
  {
    VecStr requiredCols{library_sample_name_colName, read_count_colName};
    TableReader reader(TableIOOpts::genTabFileIn(rawCountsFnp));
    reader.header_.checkForColumnsThrow(requiredCols, __PRETTY_FUNCTION__ );
    VecStr row;
    //read in
    VecStr multiple_library_sample_names;
    VecStr missingSamples;

    while (reader.getNextRow(row)) {
      const auto & exp_samp = row[reader.header_.getColPos(library_sample_name_colName)];
      auto raw_read_count_str = row[reader.header_.getColPos(read_count_colName)];
      if(!isDoubleStr(raw_read_count_str)) {
        std::stringstream ss;
        ss << __FILE__ << " " << __LINE__ << " " << __PRETTY_FUNCTION__ << ", error " << read_count_colName << " should be a number, this does not look like a number: " <<  raw_read_count_str << "\n";
        throw std::runtime_error{ss.str()};
      }
      uint32_t read_counts = njh::StrToNumConverter::stoToNum<uint32_t>(raw_read_count_str);
      if (njh::in(exp_samp, library_sample_raw_read_counts)) {
        multiple_library_sample_names.emplace_back(exp_samp);
      }
      if (njh::notIn(exp_samp, library_sample_indexes)) {
        missingSamples.emplace_back(exp_samp);
      }
      library_sample_raw_read_counts[exp_samp] = read_counts;
    }
    if (!multiple_library_sample_names.empty()) {
      std::stringstream ss;
      ss << __FILE__ << " " << __LINE__ << " " << __PRETTY_FUNCTION__ << ", error " << "found multiple of " << njh::conToStr(multiple_library_sample_names, ",") << " in " << rawCountsFnp << "\n";
      throw std::runtime_error{ss.str()};
    }
    if (!missingSamples.empty()) {
      std::stringstream ss;
      ss << __FILE__ << " " << __LINE__ << " " << __PRETTY_FUNCTION__ << ", error " << "found " << njh::conToStr(missingSamples, ",") << " in " << rawCountsFnp << " but not in " << library_sample_info_fnp << "\n";
      throw std::runtime_error{ss.str()};
    }
  }

  {//
    VecStr requiredCols{library_sample_name_colName, target_name_colName, read_count_colName, stage_colName};
    TableReader reader(TableIOOpts::genTabFileIn(readsByStageFnp));
    reader.header_.checkForColumnsThrow(requiredCols, __PRETTY_FUNCTION__ );
    VecStr row;
    //read in
    std::set<std::string> missingSamples;
    std::unordered_map<std::string, std::unordered_map<std::string, VecStr>> multiple_stage_for_library_sample_for_target;
    VecStr missing_targets;
    std::unordered_map<std::string, std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>>> exp_samp_stage_counts;
    while (reader.getNextRow(row)) {
      const std::string & library_sample_name = row[reader.header_.getColPos(library_sample_name_colName)];
      if (njh::notIn(library_sample_name, library_sample_raw_read_counts)) {
        missingSamples.emplace(library_sample_name);
      }

      const std::string & target_name = row[reader.header_.getColPos(target_name_colName)];
      const std::string & read_count_str = row[reader.header_.getColPos(read_count_colName)];
      const std::string & stage = row[reader.header_.getColPos(stage_colName)];
      if (njh::in(stage, exp_samp_stage_counts[library_sample_name][target_name])) {
        multiple_stage_for_library_sample_for_target[library_sample_name][target_name].emplace_back(stage);
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
      exp_samp_stage_counts[library_sample_name][target_name][stage] = read_counts;
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
    if (!multiple_stage_for_library_sample_for_target.empty()) {
      std::stringstream ss;
      VecStr warnings;
      for (const auto & exp_name : multiple_stage_for_library_sample_for_target) {
        for (const auto & tar_name : exp_name.second) {
          warnings.emplace_back(njh::pasteAsStr(exp_name.first, " ", tar_name.first, ":" , njh::conToStr(tar_name.second, ",")));
        }
      }
      ss << __FILE__ << " " << __LINE__ << " " << __PRETTY_FUNCTION__ << ", error " << "found multiple stage counts for one or more library_sample_names  " << njh::conToStr(njh::getVecOfMapKeys(multiple_stage_for_library_sample_for_target), ",") << " in " << readsByStageFnp << "\n";
      ss << njh::conToStr(warnings, "\n") << "\n";
      throw std::runtime_error{ss.str()};
    }

    for (const auto & library_sample_name : library_sample_raw_read_counts) {
      Json::Value sampleJson;
      sampleJson["library_sample_id"] = library_sample_indexes[library_sample_name.first];
      sampleJson["total_raw_count"] = library_sample_name.second;
      if (njh::in(library_sample_name.first, exp_samp_stage_counts)) {
        auto & read_counts_for_targets = sampleJson["read_counts_for_targets"];
        for (const auto & tar_name : exp_samp_stage_counts[library_sample_name.first]) {
          Json::Value tar_json;
          tar_json["target_id"] = target_indexes[tar_name.first];
          auto & stages = tar_json["stages"];
          for (const auto & stage : tar_name.second) {
            Json::Value stageJson;
            stageJson["stage"] = stage.first;
            stageJson["reads"] = stage.second;
            stages.append(stageJson);
          }
          read_counts_for_targets.append(tar_json);
        }
      }
      read_counts_by_library_sample_by_stage.append(sampleJson);
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

int ampliconAnalysisRunner::specimenLibraryInfoFileToJson(const njh::progutils::CmdArgs &inputCommands) {
  OutOptions outOpts("", ".json");
  bfs::path libraryInfoFnp;
  bfs::path specimenInfoFnp;
  bfs::path projectInfoFnp;
  uint32_t sequencing_info_id = 0;
  uint32_t panel_id = 0;
  ampliconAnalysisSetUp setUp(inputCommands);
  setUp.setOption(specimenInfoFnp, "--specimenInfoFnp", "Name specimen Info Fnp", true);
  setUp.setOption(libraryInfoFnp, "--libraryInfoFnp", "Name specimen Info Fnp", true);
  setUp.setOption(projectInfoFnp, "--projectInfoFnp", "Name project Info Fnp", true);

  setUp.setOption(sequencing_info_id, "--sequencing_info_id", "sequencing_info_id", true);
  setUp.setOption(panel_id, "--panel_id", "panel_id", true);


  setUp.processWritingOptions(outOpts);
  setUp.finishSetUp(std::cout);

  OutputStream out(outOpts);
  VecStr plateInfoCols{"plate_name", "plate_row", "plate_col"};

  Json::Value outJson;
  auto & project_info = outJson["project_info"];
  auto & specimen_info = outJson["specimen_info"];
  auto & library_sample_info = outJson["library_sample_info"];
  std::unordered_map<std::string, uint32_t> project_name_indexes;

  {
    table reader(TableIOOpts::genTabFileIn(projectInfoFnp));
    VecStr projectRequiredCols{"project_name",
                 "project_description"};
    reader.checkForColumnsThrow(projectRequiredCols, __PRETTY_FUNCTION__ );

    //read in
    uint32_t project_name_index = 0;
    VecStr multiple_project_names;
    for (const auto & row : reader) {
      std::string project_name = row[reader.getColPos("project_name")];
      if (njh::in(project_name, project_name_indexes)) {
        multiple_project_names.emplace_back(project_name);
      }
      project_name_indexes[project_name] = project_name_index;
      ++project_name_index;
      Json::Value projectJson;
      for(const auto & colName : reader.columnNames_){
        if(colName == "project_name"){
          projectJson[colName] = row[reader.getColPos(colName)];
        } else {
          auto colPos = reader.getColPos(colName);
          const auto & currentColValue = row[colPos];
          projectJson[colName] = currentColValue;
        }
      }
      project_info.append(projectJson);
    }
    if (!multiple_project_names.empty()) {
      std::stringstream ss;
      ss << __FILE__ << " " << __LINE__ << " " << __PRETTY_FUNCTION__ << ", error " << "found the following projects multiple times: " << njh::conToStr(multiple_project_names, ",") << "\n";
      throw std::runtime_error{ss.str()};
    }
  }

  std::unordered_map<std::string, uint32_t> specimen_name_indexes;
  {
    table reader(TableIOOpts::genTabFileIn(specimenInfoFnp));
    auto numeric_cols = reader.getNumericColumnPositions();
    // VecStr specimenRequiredCols{"specimen_name", "samp_taxon_id", "collection_date", "collection_country", "collector", "samp_store_loc", "samp_collect_device", "project_name"};
   VecStr specimenRequiredCols{"specimen_name",
                    "specimen_taxon_id",
                    "host_taxon_id",
                    "collection_date",
                    "collection_country",
                    "project_name"};
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
        } else if (colName == "project_name") {
          auto colPos = reader.getColPos(colName);
          const auto & currentColValue = row[colPos];
          if (njh::notIn(currentColValue, project_name_indexes)) {
            std::stringstream ss;
            ss << __PRETTY_FUNCTION__ << ", error " << projectInfoFnp << " does not contain project name: " << currentColValue << "\n";
            throw std::runtime_error{ss.str()};
          }
          Json::Value project_id = njh::json::toJson(project_name_indexes[currentColValue]);
          sampleJson["project_id"] = project_id;
        }  else if (colName == "specimen_taxon_id") {
          auto colPos = reader.getColPos(colName);
          const auto & currentColValue = row[colPos];
          Json::Value specimen_taxon_id = njh::json::toJson(njh::StrToNumConverter::stoToNum<uint32_t>(currentColValue));
          sampleJson["specimen_taxon_id"].append(specimen_taxon_id);
        } else if (njh::in(colName, VecStr{"parasite_density_method", "parasite_density"}) || njh::in(colName, plateInfoCols)) {
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
          parasite_densityJson["parasite_density_method"] = row[reader.getColPos("parasite_density_method")];
          parasite_densityJson["parasite_density"] = njh::json::toJson(njh::StrToNumConverter::stoToNum<double>(row[reader.getColPos("parasite_density")]));
          sampleJson["parasite_density_info"].append(parasite_densityJson);
        }
      }


      if(row[reader.getColPos("plate_name")] != "NA") {
        Json::Value plateJson;
        plateJson["plate_name"] = row[reader.getColPos("plate_name")];
        plateJson["plate_row"] = row[reader.getColPos("plate_row")];
        plateJson["plate_col"] = njh::json::toJson(njh::StrToNumConverter::stoToNum<uint32_t>(row[reader.getColPos("plate_col")]));
        sampleJson["storage_plate_info"] = plateJson;
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
    table reader(TableIOOpts::genTabFileIn(libraryInfoFnp));
    reader.checkForColumnsThrow(toVecStr(VecStr{"library_sample_name", "specimen_name"}), __PRETTY_FUNCTION__ );
    auto numeric_cols = reader.getNumericColumnPositions();


    //read in
    VecStr library_sample_names;
    VecStr multiple_library_sample_names;
    for (const auto & row : reader) {
      std::string library_sample_name = row[reader.getColPos("library_sample_name")];
      if (njh::in(library_sample_name, library_sample_names)) {
        multiple_library_sample_names.emplace_back(library_sample_name);
      }
      library_sample_names.emplace_back(library_sample_name);
      Json::Value librarySampleJson;
      librarySampleJson["sequencing_info_id"] = sequencing_info_id;
      librarySampleJson["panel_id"] = panel_id;
      for (const auto& colName: reader.columnNames_) {
        if (colName == "library_sample_name") {
          librarySampleJson[colName] = row[reader.getColPos(colName)];
        } else if (colName == "specimen_name") {
          librarySampleJson["specimen_id"] = specimen_name_indexes[row[reader.getColPos(colName)]];
        } else if ( njh::in(colName, plateInfoCols)) {
          //do nothing
        } else {
          auto colPos = reader.getColPos(colName);
          const auto& currentColValue = row[colPos];
          if(njh::in(colPos, numeric_cols)){
            if (njh::strAllDigits(currentColValue)) {
              librarySampleJson[colName] = njh::json::toJson(
                njh::StrToNumConverter::stoToNum<uint32_t>(currentColValue));
            } else {
              librarySampleJson[colName] = njh::json::toJson(
                njh::StrToNumConverter::stoToNum<double>(currentColValue));
            }
          } else {
            librarySampleJson[colName] = currentColValue;
          }
        }
      }
      if(row[reader.getColPos("plate_name")] != "NA") {
        Json::Value plateJson;

        plateJson["plate_name"] = row[reader.getColPos("plate_name")];
        plateJson["plate_row"] = row[reader.getColPos("plate_row")];
        plateJson["plate_col"] = njh::json::toJson(njh::StrToNumConverter::stoToNum<uint32_t>(row[reader.getColPos("plate_col")]));
        librarySampleJson["library_prep_plate_info"] = plateJson;
      }
      library_sample_info.append(librarySampleJson);
    }
    if (!multiple_library_sample_names.empty()) {
      std::stringstream ss;
      ss << __FILE__ << " " << __LINE__ << " " << __PRETTY_FUNCTION__ << ", error " << "found the following specimens multiple times: " << njh::conToStr(multiple_library_sample_names, ",") << "\n";
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
  uint32_t bioinformatics_run_id = std::numeric_limits<uint32_t>::max();
  bfs::path sample_to_bio_run_id_fnp;
  bfs::path panel_target_info_fnp;
  bfs::path library_sample_info_fnp;
  bfs::path regions_analyzed_bed_fnp;

  bfs::path finalClustersFnp;
  std::string sampleIDCol = "s_Sample";
  std::string targetIDCol = "p_name";
  // std::string microhaplotypeIDCol = "h_popUID";
  std::string readCountCol = "c_ReadCnt";
  std::string umiCountCol = "c_barcodeCnt";

  std::string seqCol = "h_Consensus";
  uint32_t regions_analyzed_genome_id = 0;
  ampliconAnalysisSetUp setUp(inputCommands);
  // setUp.setOption(sequencing_id, "--sequencing_id", "sequencing id", true);
  setUp.setOption(bioinformatics_run_id, "--bioinformatics_run_id", "bioinformatics_run_id", false);
  setUp.setOption(sample_to_bio_run_id_fnp, "--sample_to_bio_run_id_fnp", "a table with sample to bio run info id, one column is the same sample column as the input data, other is bioinformatics_run_id", false);
  setUp.setOption(panel_target_info_fnp, "--panel_target_info_fnp", "json file containing the information about the panel and target", true);
  setUp.setOption(library_sample_info_fnp, "--library_sample_info_fnp", "json file containing the information library_samples", true);

  setUp.setOption(finalClustersFnp, "--finalClustersFnp", "Name extracted Info Fnp", true);
  setUp.setOption(regions_analyzed_bed_fnp, "--regions_analyzed_bed_fnp", "If the region analyzed was smaller than the full insert region, provide regions here, have to provide all regions for all target if providing any", false);
  setUp.setOption(regions_analyzed_genome_id, "--regions_analyzed_genome_id", "If the region analyzed was smaller than the full insert region, the genome index for the regions in --regions_analyzed_bed_fnp", false);

  setUp.setOption(sampleIDCol, "--sampleIDCol", "sampleIDCol");
  setUp.setOption(targetIDCol, "--targetIDCol", "targetIDCol");
  // setUp.setOption(microhaplotypeIDCol, "--microhaplotypeIDCol", "microhaplotypeIDCol");
  setUp.setOption(readCountCol, "--readCountCol", "readCountCol");
  setUp.setOption(seqCol, "--seqCol", "seqCol");

  setUp.setOption(umiCountCol, "--umiCountCol", "umiCountCol");

  setUp.processWritingOptions(outOpts);
  setUp.finishSetUp(std::cout);

  // get regions analyzed if loaded 
  std::vector<std::shared_ptr<Bed6RecordCore>> regions_analyzed;
  std::unordered_map<std::string, std::shared_ptr<Bed6RecordCore>>regions_analyzed_by_name;
  if (!regions_analyzed_bed_fnp.empty()) {
    regions_analyzed = getBeds(regions_analyzed_bed_fnp);
  }
  VecStr duplicate_targets;
  for (const auto & region : regions_analyzed) {
    if (njh::in(region->name_, regions_analyzed_by_name)) {
      duplicate_targets.emplace_back(region->name_);
    }
    regions_analyzed_by_name[region->name_] = region;
  }
  if (!duplicate_targets.empty()) {
    std::stringstream ss;
    ss << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << ", error " << " loaded in duplicate regions from " << regions_analyzed_bed_fnp << " for targets: " << njh::conToStr(duplicate_targets, ",") << "\n";
    throw std::runtime_error{ss.str()};
  }
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

  auto library_sample_info = njh::json::parseFile(library_sample_info_fnp.string());
  if (!library_sample_info.isMember("library_sample_info")) {
    std::stringstream ss;
    ss << __FILE__ << " " << __LINE__ << " " << __PRETTY_FUNCTION__ << ", error " << library_sample_info_fnp << " needs to have library_sample_info, only has " << njh::conToStr(library_sample_info.getMemberNames(), ",") << "\n";
    throw std::runtime_error{ss.str()};
  }

  std::unordered_map<std::string, uint32_t> library_sample_indexes;
  VecStr multiple_library_sample_names;
  for (const auto & exp_samp_enum : iter::enumerate(library_sample_info["library_sample_info"])) {
    auto library_sample_name = exp_samp_enum.second["library_sample_name"].asString();
    if (njh::in(library_sample_name, library_sample_indexes)) {
      multiple_library_sample_names.emplace_back(library_sample_name);
    }
    library_sample_indexes[library_sample_name] = exp_samp_enum.index;
  }

  std::unordered_map<std::string, uint32_t> sample_to_bio_run_id_map;
  if (std::numeric_limits<uint32_t>::max() != bioinformatics_run_id) {
    for (const auto & library_sample_indx : library_sample_indexes) {
      sample_to_bio_run_id_map[library_sample_indx.first] = bioinformatics_run_id;
    }
  } else {
    table sample_to_bio_run_id_tab(sample_to_bio_run_id_fnp, "\t", true);
    sample_to_bio_run_id_tab.checkForColumnsThrow(VecStr{sampleIDCol, "bioinformatics_run_id"}, __PRETTY_FUNCTION__);
    for (const auto & row : sample_to_bio_run_id_tab) {
      auto samp_name = row[sample_to_bio_run_id_tab.getColPos(sampleIDCol)];
      if (njh::in(samp_name, sample_to_bio_run_id_map)) {
        std::stringstream ss;
        ss << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << ", error " << " already have a bio id for " << samp_name << "\n";
        ss << "have " << sample_to_bio_run_id_map[samp_name] << " and now attempting to add " << row[sample_to_bio_run_id_tab.getColPos("bioinformatics_run_id")] << "\n";
        throw std::runtime_error{ss.str()};
      }
      sample_to_bio_run_id_map[samp_name] = njh::StrToNumConverter::stoToNum<uint32_t>(row[sample_to_bio_run_id_tab.getColPos("bioinformatics_run_id")]);
    }

    VecStr only_in_library_input;
    VecStr only_in_samp_to_bio_id_tab;
    VecStr inboth;
    auto library_input = getVectorOfMapKeys(library_sample_indexes);
    auto samp_to_bio_id = getVectorOfMapKeys(sample_to_bio_run_id_map);
    njh::sort(library_input);
    njh::sort(samp_to_bio_id);

    njh::decompose_sets(library_input.begin(), library_input.end(),
      samp_to_bio_id.begin(), samp_to_bio_id.end(),
      std::back_insert_iterator(only_in_library_input),
      std::back_insert_iterator(only_in_samp_to_bio_id_tab),
      std::back_insert_iterator(inboth)
      );
    if (!only_in_library_input.empty() ||  !only_in_samp_to_bio_id_tab.empty()) {
      std::stringstream ss;
      ss << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << ", error " << " " << "\n";
      if (!only_in_library_input.empty()) {
        ss << "have " << njh::conToStr(only_in_library_input, ",") << " in " << library_sample_info_fnp << " but missing from "  << sample_to_bio_run_id_fnp << "\n";
      }
      if (!only_in_samp_to_bio_id_tab.empty()) {
        ss << "have " << njh::conToStr(only_in_samp_to_bio_id_tab, ",") << " in " << sample_to_bio_run_id_fnp << " but missing from "  << library_sample_info_fnp << "\n";
      }
      throw std::runtime_error{ss.str()};
    }
  }

  if (!multiple_library_sample_names.empty()) {
    std::stringstream ss;
    ss << __FILE__ << " " << __LINE__ << " " << __PRETTY_FUNCTION__ << ", error " << "found multiple of the same target in " << library_sample_info_fnp << ": " << njh::conToStr(multiple_library_sample_names, ",") << "\n";
    throw std::runtime_error{ss.str()};
  }
  OutputStream out(outOpts);

  Json::Value outJson;
  Json::Value & detected_microhaplotypes_full = outJson["detected_microhaplotypes"];
  Json::Value & representative_microhaplotypes = outJson["representative_microhaplotypes"];


  std::unordered_map<uint32_t, Json::Value> map_detected_microhaplotypes;
  auto bioinformatics_run_ids = getVectorOfMapValues(sample_to_bio_run_id_map);
  for (const auto & current_bioinformatics_run_id : bioinformatics_run_ids) {
    map_detected_microhaplotypes[current_bioinformatics_run_id]["bioinformatics_run_id"] = current_bioinformatics_run_id;
  }

  // Json::Value & samplesJson = detected_microhaplotypes["library_samples"];
  Json::Value & representative_microhaplotypes_targets = representative_microhaplotypes["targets"];


  std::map<std::string, std::vector<std::shared_ptr<seqInfo>>> popSeqsByTarget;

  //representative_microhaplotypes
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
  {
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
  }
  if (!regions_analyzed.empty()) {
    VecStr missingTargetNames;
    for (auto & popSeqsForTarget : popSeqsByTarget) {
      if (njh::notIn(popSeqsForTarget.first, regions_analyzed_by_name)) {
        missingTargetNames.emplace_back(popSeqsForTarget.first);
      }
    }
    if (!missingTargetNames.empty()) {
      std::stringstream ss;
      ss << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << ", error " << "missing the following targets found in " << finalClustersFnp << " from " << regions_analyzed_bed_fnp << "\n";
      ss << njh::conToStr(missingTargetNames, "\n") << "\n";
      throw std::runtime_error{ss.str()};
    }
  }
  std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>> popSeqsByTargetMhapIndex;
  std::unordered_map<std::string, uint32_t> popSeqsByTargetIndex;

  uint32_t microhaplotypes_target_index = 0;
  for (auto & popSeqsForTarget : popSeqsByTarget) {
    Json::Value microhaps_for_target_info;
    microhaps_for_target_info["target_id"] = target_indexes[popSeqsForTarget.first];
    if (!regions_analyzed.empty()) {
      Json::Value mhap_location_for_target = GenomicRegion(*regions_analyzed_by_name[popSeqsForTarget.first]).toJsonLocationOnly();
      mhap_location_for_target["genome_id"] = regions_analyzed_genome_id;
      microhaps_for_target_info["mhap_location"] = mhap_location_for_target;
    }
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
    representative_microhaplotypes_targets.append(microhaps_for_target_info);
  }

  {
    //detected_microhaplotypes
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
      samp_info["library_sample_id"] = library_sample_indexes[samp.first];
      for (const auto & tar : samp.second) {
        Json::Value tar_info;
        tar_info["mhaps_target_id"] = popSeqsByTargetIndex[tar.first];
        tar_info["mhaps"] = tar.second;
        samp_info["target_results"].append(tar_info);
      }
      map_detected_microhaplotypes[sample_to_bio_run_id_map[samp.first]]["library_samples"].append(samp_info);
    }
  }
  for (const auto & detected_microhaplotypes : map_detected_microhaplotypes) {
    detected_microhaplotypes_full.append(detected_microhaplotypes.second);
  }
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
  std::string default_reaction;
  bfs::path reactionNameFnp;
  bfs::path extractedInfoFnp;
  bfs::path genomeTwoBit;
  bfs::path additionalTargetAttributes;
  std::string targetColName = "target";

  ampliconAnalysisSetUp setUp(inputCommands);

  auto reaction_set = setUp.setOption(reactionNameFnp, "--reactionNameFnp", "table with target column and reaction column", false);
  setUp.setOption(default_reaction, "--default_reaction", "use this as the default reaction for all targets", !reaction_set);
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
  if (reaction_set) {
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
    if (reaction_set) {
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
    } else {
      for(const auto & extractedTargets : extractedRowsPerID) {
        reactionNameMap[extractedTargets.first].emplace_back(default_reaction);
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
          if (reader.header_.containsColumn("insertGeneDescription") && !row[reader.header_.getColPos("insertGeneDescription")].empty()) {
            tarInfo["gene_name"] = row[reader.header_.getColPos("insertGeneID")];
          }
          // Json::Value forwardPrimers;
          Json::Value forwardPrimer;
          forwardPrimer["seq"] = row[reader.header_.getColPos("Fwd_primer")];
          forwardPrimer["location"] = fprimer.toJsonLocationOnly();
          forwardPrimer["location"]["genome_id"] = 0;
          // forwardPrimers.append(forwardPrimer);
          tarInfo["forward_primer"] = forwardPrimer;

          // Json::Value reversePrimers;
          Json::Value reversePrimer;
          reversePrimer["seq"] = row[reader.header_.getColPos("Rev_primer")];
          reversePrimer["location"] = fprimer.toJsonLocationOnly();
          reversePrimer["location"]["genome_id"] = 0;
          // reversePrimers.append(reversePrimer);
          tarInfo["reverse_primer"] = reversePrimer;
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

