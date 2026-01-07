//
// Created by Nicholas Hathaway on 1/2/26.
//


#include "programWrappers.hpp"
#include <njhseq/objects/BioDataObject/BioDataFileIO.hpp>
#include <njhseq/objects/BioDataObject/BioRecordsUtils/HmmerUtility.hpp>
#include <njhseq/objects/BioDataObject/NucmerRecord.hpp>


namespace njhseq {
//nucmer  ../../refSeqs/PfDd2_vars_allcDNA.fasta vars.fasta --prefix bc07_PfDd2_vars_nucmer 2> nucmer.log &&
//show-coords -T -l  -c -H bc07_PfDd2_vars_nucme.delta | sort | uniq |
//elucidator parseNucmerResultsToBed  --coordsOutput STDIN --overWrite --out bc07_PfDd2_vars_nucmer.delta.bed
int programWrapperRunner::runNucmerToBed(const njh::progutils::CmdArgs & inputCommands) {
  std::string extraParameters;//"--mum -b 100 -l 31";
  bfs::path query_fnp;
  bfs::path ref_fnp;
  std::string prefix;
  bool do_not_use_temp_location = false;
  bfs::path temp_dir_location = "/tmp/";
  bool keep_temp_files = false;
  OutOptions outOpt("", ".bed");

  seqSetUp setUp(inputCommands);
  setUp.processVerbose();
  setUp.processDebug();
  setUp.setOption(query_fnp, "--query_fnp", "query fasta file", true);
  setUp.setOption(ref_fnp, "--ref_fnp", "reference fasta file", true);
  setUp.setOption(prefix, "--prefix", "prefix for output of nucmer results", true);
  setUp.setOption(keep_temp_files, "--keep_temp_files", "keep temp files");
  setUp.setOption(do_not_use_temp_location, "--do_not_use_temp_location", "whether or not to use a temp directory, helpful for when files are in a location with spaces in the filename path cause nucmer can't handle this");
  bool use_temp_location = !do_not_use_temp_location;
  setUp.setOption(temp_dir_location, "--temp_dir_location", "the temp directory to use");
  setUp.processWritingOptions(outOpt);
  if (outOpt.outFilename_.string().empty()) {
    outOpt.outFilename_ = prefix + ".bed";
  }
  setUp.finishSetUp(std::cout);

  bfs::path run_dir = "./";
  bfs::path new_temp_dir;
  bfs::path new_query_fnp;
  bfs::path new_ref_fnp;
  bfs::path new_prefix;
  bfs::path original_prefix = prefix;
  OutputStream bed_out(outOpt);

  std::vector<bfs::path> files_to_remove;
  if (use_temp_location) {
    new_temp_dir = njh::files::make_path(temp_dir_location, "nucmer_run_" + njh::getCurrentDate());
    new_temp_dir = njh::files::findNonexitantFile(new_temp_dir);
    njh::files::makeDir(njh::files::MkdirPar{new_temp_dir});
    new_ref_fnp = njh::files::make_path(new_temp_dir, ref_fnp.filename());
    new_query_fnp = njh::files::make_path(new_temp_dir, query_fnp.filename());
    new_prefix = njh::files::make_path(new_temp_dir, prefix);
    bfs::copy(ref_fnp, new_ref_fnp);
    bfs::copy(query_fnp, new_query_fnp);
    ref_fnp = new_ref_fnp;
    query_fnp = new_query_fnp;
    run_dir = new_temp_dir;
    prefix = new_prefix.string();
    files_to_remove.emplace_back(new_temp_dir);
  }
  // nucmer only works on unzipped fasta files, unzip if necessary
  if (njh::endsWith(ref_fnp.string(), ".gz")) {
    bfs::path unziped_ref_fnp = ref_fnp.string().substr(0, ref_fnp.string().size() - 3);
    if (!bfs::exists(unziped_ref_fnp)) {
      InputStream in(ref_fnp);
      OutputStream out(unziped_ref_fnp);
      out << in.rdbuf();
      files_to_remove.emplace_back(unziped_ref_fnp);
    }
    ref_fnp = unziped_ref_fnp;
  }

  if (njh::endsWith(query_fnp.string(), ".gz")) {
    bfs::path unziped_query_fnp = query_fnp.string().substr(0, query_fnp.string().size() - 3);
    if (!bfs::exists(unziped_query_fnp)) {
      InputStream in(query_fnp);
      OutputStream out(unziped_query_fnp);
      out << in.rdbuf();
      files_to_remove.emplace_back(unziped_query_fnp);
    }
    query_fnp = unziped_query_fnp;
  }

  std::string nucmer_cmd = njh::pasteAsStr("cd ", run_dir.string(), " && ",
    "nucmer ", extraParameters, " ", ref_fnp.string(), " ", query_fnp.string(), " --prefix ", prefix, " 2> /dev/null && ",
    "show-coords -T -l  -c -H ", prefix, ".delta > ", prefix, ".raw.coords && cat ", prefix, ".raw.coords | sort | uniq > ", prefix, ".coords"
    );

  auto run_res = njh::sys::run({nucmer_cmd});
  if (!run_res.success_) {
    std::stringstream ss;
    ss << __PRETTY_FUNCTION__ << ", error " << " running " << nucmer_cmd << "\n";
    ss << "stdout: " << run_res.stdOut_ << "\n";
    ss << "stderr: " << run_res.stdErr_ << "\n";
    ss << "If getting parsing of delta file error, could be the full path name to the delta file has spaces in it, even when given with a relative path, this will cause show-coords to fail, either run in another directory or use --use_temp_location" << "\n";
    throw std::runtime_error{ss.str()};
  }
  bfs::path delta_fnp = njh::files::make_path(prefix + ".delta");
  bfs::path coords_fnp = njh::files::make_path(prefix + ".coords");
  bfs::path output_bed_fnp = original_prefix.string() + ".bed";

  InputStream coords_in(coords_fnp);
  std::string line;
  while(njh::files::crossPlatGetline(coords_in, line)){
    NucmerShowCoordsRecord record(line);
    bed_out << record.genBed6().toDelimStrWithExtra() << std::endl;
  }
  if (use_temp_location) {
    if (exists(delta_fnp.filename())) {
      bfs::remove(delta_fnp.filename());
    }
    if (exists(coords_fnp.filename())) {
      bfs::remove(coords_fnp.filename());
    }
    bfs::copy(delta_fnp, delta_fnp.filename());
    bfs::copy(coords_fnp, coords_fnp.filename());
  }
  if (!keep_temp_files) {
    for (const auto & fnp : files_to_remove) {
      //will work on both files and directories (including non-empty directories)
      njh::files::rmDirForce(fnp);
    }
  }

  return 0;
}


} // namespace njhseq
