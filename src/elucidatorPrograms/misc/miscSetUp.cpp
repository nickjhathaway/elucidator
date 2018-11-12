#include "miscSetUp.hpp"

namespace njhseq {

void miscSetUp::setUpListAllFiles(std::string& directory, VecStr& contains,
                                  bool& specific, bool& recursive,
                                  std::string& filesOrDirectories) {
  if (needsHelp()){

    std::cout << "listAllFiles" << std::endl;
    std::cout << "Commands, order not necessary" << std::endl;
    std::cout << "-dir [option], directory name, will default to current "
                 "directory (.)" << std::endl;
    std::cout << "-contains [options], list only the files that contain this "
                 "keyword will default to all files" << std::endl;
    std::cout << "-r, make the search recurisve, be carefull could take a "
                 "while if lots of subdirectories" << std::endl;
    std::cout << "examples, sequenceTools listAllFiles -dir data -contains "
                 "output -r, sequenceTools listAllFiles -dir ." << std::endl;
    exit(1);
  }
  setOption(directory, "-dir,-directory", "SearchDirectory");
  std::string containsStr = "";
  specific = setOption(containsStr, "-contains", "SpecificFiles");
  if (specific) {
    contains = tokenizeString(containsStr, ",");
  }
  setOption(recursive, "-r", "Recursive");
  bool files = false;
  bool directories = false;
  if (setOption(files, "-onlyfiles", "Onlyfiles")) {
    filesOrDirectories = "file";
  } else if (setOption(directories, "-onlydirectories", "OnlyDirectories")) {
    filesOrDirectories = "directory";
  }

  finishSetUp(std::cout);
}

void miscSetUp::setUpProcessKrecFasta(std::string& filename,
                                      std::string& refFilename,
                                      bool& postProcessed) {
	setOption(filename, "-file", "Filename", true);
  setOption(refFilename, "-ref", "ReferenceFilename", true);
  setOption(postProcessed, "-post", "PostProcessed");
  finishSetUp(std::cout);
}
}  // namespace njh
