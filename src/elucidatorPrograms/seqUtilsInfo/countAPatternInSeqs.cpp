//
// Created by Nicholas Hathaway on 12/24/24.
//



#include "seqUtilsInfoRunner.hpp"
#include <njhseq/IO/SeqIO/SeqIO.hpp>
#include <njhseq/seqToolsUtils/tandemRepeatUtils.hpp>


namespace njhseq {



int seqUtilsInfoRunner::countDiNucleotidePatternsInSeqs(const njh::progutils::CmdArgs & inputCommands) {
  SimpleTandemRepeatFinder::SimpleTRFinderLocsPars pars;

  OutOptions outOpts(bfs::path(""), ".tsv");
  uint32_t proceedingBases = 5;
  uint32_t trailingBases = 5;
  uint32_t minRepeatingAmount = 6;
  pars.minRepeatUnitSize = 2;
  pars.maxRepeatUnitSize = pars.minRepeatUnitSize;
  std::set<std::string> repeats;
  seqSetUp setUp(inputCommands);
  setUp.description_ = "count the pattern surrounding dinucleotide runs";

  setUp.processVerbose();
  setUp.processDebug();
  setUp.processWritingOptions(outOpts);
  setUp.setOption(proceedingBases, "--proceedingBases", "proceeding Bases");
  setUp.setOption(trailingBases, "--trailingBases", "trailing Bases");
  setUp.setOption(minRepeatingAmount, "--minRepeatingAmount", "min repeating amount (length would be 2 x this values, e.g. a min repeat amount of 4 would be 4 x 2 length of 8 bases");

  setUp.setOption(repeats, "--repeats", "repeats");
  setUp.setOption(pars.alphabet, "--bases", "bases to combine to make the dinucleotide repeats");

  setUp.processReadInNames(true);
  setUp.finishSetUp(std::cout);

  SeqIO reader(setUp.pars_.ioOptions_);
  reader.openIn();

  if (!repeats.empty()) {
    VecStr Warnings;
    for (const auto & repeat : repeats) {
      if (repeat.size() != 2) {
        Warnings.push_back("Warning: repeat " + repeat + " is not a valid repeat, must be of length 2");
      }
    }
    if (!Warnings.empty()) {
      std::stringstream ss;
      ss << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << ", error:" << "\n";
      ss << njh::pasteAsStr(Warnings, "\n") << "\n";
      throw std::runtime_error{ss.str()};
    }
  } else {

    SimpleTandemRepeatFinder finder(pars);
    auto minRepeats = finder.genMinimalUnitsNeededForSearch();
    for (const auto & repeat : *minRepeats.allUnits) {
      if (repeat.front() != repeat.back()) {
        repeats.emplace(repeat);
      }
    }
  }
  if (setUp.pars_.verbose_) {
    std::cout << njh::conToStr(repeats, "\n") << std::endl;
  }

  OutputStream out(outOpts);
  std::unordered_map<std::string, std::regex> repeatPatterns;

  for (const auto & repeat : repeats) {

    std::string patStr = njh::pasteAsStr("(.{", proceedingBases, ",", proceedingBases,"})(", njh::pasteAsStr(VecStr(minRepeatingAmount, repeat)), "(?:", repeat, ")+)(.{", trailingBases, ",", trailingBases, "})");
    if (setUp.pars_.verbose_) {
      std::cout << "patStr: " << patStr << std::endl;
    }
    std::regex pattern(patStr);
    repeatPatterns.emplace(repeat, pattern);
  }
  std::map<std::string, std::unordered_map<std::string,std::unordered_map<std::string, uint32_t>>> patternCounts;
  seqInfo seq;
  while(reader.readNextRead(seq)) {
    for (const auto & repeat : repeats) {
      std::smatch match;
      std::string::const_iterator searchStart(seq.seq_.cbegin());
      while (std::regex_search(searchStart, seq.seq_.cend(), match, repeatPatterns.at(repeat))) {
        patternCounts[repeat][njh::pasteAsStr(match[1],"-",match[3])][match[2]]++;
        searchStart = match.suffix().first;
      }
    }
  }
  out << "dinucleotide\tproceeding_trailing_bases\tfull_repeat\tcount\tfull_pattern" << std::endl;
  for (const auto & repeat : repeats) {
    if (njh::in(repeat, patternCounts)) {
      for (const auto & [pattern, patternCountsMap]  : patternCounts.at(repeat)) {
        if (patternCountsMap.size() > 1) {
          double total = 0;
          for (const auto & [pattern2, count] : patternCountsMap) {
            total += count;
          }
          for (const auto & [pattern2, count] : patternCountsMap) {
            out << repeat
                << "\t" << pattern
                << "\t" << pattern2
                << "\t" << count
                << "\t" << count / total
                << "\t" << total
                << "\t" << njh::replaceString(pattern, "-", pattern2)
                << std::endl;
          }
        }
      }
    }
  }
  return 0;
}

int seqUtilsInfoRunner::countHPPatternsInSeqs(const njh::progutils::CmdArgs & inputCommands) {

  OutOptions outOpts(bfs::path(""), ".tsv");
  uint32_t proceedingBases = 5;
  uint32_t trailingBases = 5;
  uint32_t minHomopolymerLength = 5;
  std::vector<char> homopolymerBases = {'A', 'C', 'G', 'T'};
  std::vector<char> allBases = {'A', 'C', 'G', 'T'};

  seqSetUp setUp(inputCommands);
  setUp.description_ = "count the pattern surrounding homopolymer runs";

  setUp.processVerbose();
  setUp.processDebug();
  setUp.processWritingOptions(outOpts);
  setUp.setOption(proceedingBases, "--proceedingBases", "proceeding Bases", njh::progutils::ProgramSetUp::CheckCase::GT1);
  setUp.setOption(trailingBases, "--trailingBases", "trailing Bases", njh::progutils::ProgramSetUp::CheckCase::GT1);
  setUp.setOption(minHomopolymerLength, "--minHomopolymerLength", "min Homopolymer Length", njh::progutils::ProgramSetUp::CheckCase::GT1);
  setUp.setOption(homopolymerBases, "--homopolymerBases", "homopolymer Bases");
  setUp.setOption(allBases, "--allBases", "all Bases");


  setUp.processReadInNames(true);
  setUp.finishSetUp(std::cout);

  SeqIO reader(setUp.pars_.ioOptions_);
  reader.openIn();
  OutputStream out(outOpts);
  std::unordered_map<char, std::regex> basePatterns;
  std::string allBasesStr = njh::pasteAsStr(allBases);

  for (const auto base : homopolymerBases) {
    auto allBasesButHpBase = allBases;
    removeElement(allBasesButHpBase, base);
    std::string allBasesButHpBaseStr = njh::pasteAsStr(allBasesButHpBase);

    std::string patStr = njh::pasteAsStr("([", allBasesStr, "]", "{", proceedingBases - 1, ",", proceedingBases -1,"}","[", allBasesButHpBaseStr,"]",")(", base, "{", minHomopolymerLength, ",})(","[", allBasesButHpBaseStr,"]","[", allBasesStr, "]", "{", trailingBases - 1, ",", trailingBases - 1, "})");
    // std::cout << "patStr: " << patStr << std::endl;
    std::regex pattern(patStr);
    basePatterns.emplace(base, pattern);
  }
  std::map<char, std::unordered_map<std::string,std::unordered_map<std::string, uint32_t>>> patternCounts;
  seqInfo seq;
  while(reader.readNextRead(seq)) {
    if (len(seq) > proceedingBases + trailingBases + minHomopolymerLength) {
      for (const auto base : homopolymerBases) {
        std::smatch match;
        std::string::const_iterator searchStart(seq.seq_.cbegin());
        while (std::regex_search(searchStart, seq.seq_.cend(), match, basePatterns.at(base))) {
          patternCounts[base][njh::pasteAsStr(match[1],"-",match[3])][match[2]]++;
          //searchStart = match.suffix().first;
          if (match[1].str().back() == base || match[3].str().front() == base) {
            std::cout << seq.name_ << std::endl;
            std::cout << seq.seq_ << std::endl;
            std::cout << "searchStart: " << searchStart - seq.seq_.cbegin() << std::endl;
            std::cout << seq.seq_.substr(0, searchStart - seq.seq_.cbegin()) << std::endl;
            std::cout << "\t" << njh::pasteAsStr(match[1],"-",match[3]) << std::endl;
            std::cout << "\t" << njh::pasteAsStr(match[2]) << std::endl;
            exit(1);
          }
          searchStart = seq.seq_.cbegin() + (match.prefix().second - seq.seq_.cbegin()) + match.length(2);

        }
      }
    }
  }

  out << "base\tproceeding_trailing_bases\thomopolymer\tcount\tfreq\ttotal\tfull_pattern" << std::endl;
  for (const auto base : homopolymerBases) {
    if (njh::in(base, patternCounts)) {
      for (const auto & [pattern, patternCountsMap]  : patternCounts.at(base)) {
        if (patternCountsMap.size() > 1) {
          double total = 0;
          for (const auto & [pattern2, count] : patternCountsMap) {
            total += count;
          }
          for (const auto & [pattern2, count] : patternCountsMap) {
            out << base
                << "\t" << pattern
                << "\t" << pattern2
                << "\t" << count
                << "\t" << count / total
                << "\t" << total
                << "\t" << njh::replaceString(pattern, "-", pattern2)
                << std::endl;
          }
        }
      }
    }
  }
  return 0;
}


int seqUtilsInfoRunner::countAPatternInSeqs(const njh::progutils::CmdArgs & inputCommands) {

  OutOptions outOpts(bfs::path(""), ".tsv");
  std::string patternStr = "G{10,}$";
  seqSetUp setUp(inputCommands);
  setUp.description_ = "count the presence of a pattern within input sequences";

  setUp.processVerbose();
  setUp.processDebug();
  setUp.processWritingOptions(outOpts);
  setUp.setOption(patternStr, "--pattern", "pattern to search for", true);
  setUp.processReadInNames(true);
  setUp.finishSetUp(std::cout);

  SeqIO reader(setUp.pars_.ioOptions_);
  reader.openIn();
  OutputStream out(outOpts);

  std::regex pattern{patternStr};
  if (setUp.pars_.ioOptions_.isPairedIn()) {
    out << "inputFnp\tpattern\tr1_foundCount\tr2_coundCount\tinBoth_foundCount\ttotalInputPairedCount" << std::endl;

    PairedRead pseq;
    uint32_t r1_foundCount = 0;
    uint32_t r2_coundCount = 0;
    uint32_t inBoth_foundCount = 0;
    uint32_t totalInputPairedCount = 0;
    while(reader.readNextRead(pseq)){
      ++totalInputPairedCount;
      std::smatch firstMate_match;
      if (std::regex_search(pseq.seqBase_.seq_, firstMate_match, pattern)) {
        ++r1_foundCount;
      }
      std::smatch secondMate_match;
      if (std::regex_search(pseq.mateSeqBase_.seq_, secondMate_match, pattern)) {
        ++r2_coundCount;
      }
      if (firstMate_match.size() == 1 && secondMate_match.size() == 1) {
        ++inBoth_foundCount;
      }
    }
    out << setUp.pars_.ioOptions_.firstName_
        << "\t" << patternStr
        << "\t" << r1_foundCount
        << "\t" << r2_coundCount
        << "\t" << inBoth_foundCount
        << "\t" << totalInputPairedCount << std::endl;
  } else {
    seqInfo seq;
    out << "inputFnp\tpattern\tfoundCount\ttotalInputCount" << std::endl;
    uint32_t totalInputCount = 0;
    uint32_t foundCount = 0;
    while(reader.readNextRead(seq)){
      ++totalInputCount;
      std::smatch match;
      if (std::regex_search(seq.seq_, match, pattern)) {
        ++foundCount;
      }
    }
    out << setUp.pars_.ioOptions_.firstName_
        << "\t" << patternStr
        << "\t" << foundCount
        << "\t" << totalInputCount << std::endl;
  }
  if(setUp.pars_.verbose_){
    setUp.logRunTime(std::cout);
  }
  return 0;
}


} //namespace njhseq

