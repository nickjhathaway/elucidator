//
// Created by Nicholas Hathaway on 2/2/25.
//
#include <njhseq/alignment/aligner/aligner.hpp>

#include "seqUtilsModRunner.hpp"
#include <njhseq/objects/counters/DNABaseCounter.hpp>
#include <njhseq/IO/SeqIO/SeqIO.hpp>
#include <njhseq/objects/seqObjects/seqKmers.h>
#include <njhseq/readVectorManipulation/readVectorHelpers/readVecSorter.hpp>


namespace njhseq {


int seqUtilsModRunner::correctHPRunsBasedOnReference(const njh::progutils::CmdArgs & inputCommands) {

  uint32_t min_hp_run_size_to_correct = 8;
  uint32_t max_hp_run_gap_size = 3;


  seqInfo ref_seq;
  seqSetUp setUp(inputCommands);
  setUp.description_ = "correct possible errors in homopolymer runs by compring to a reference";
  setUp.processVerbose();
  setUp.processDebug();
  setUp.processGap();
  setUp.processScoringPars();
  setUp.processAlnInfoInput();
  setUp.processSeq(ref_seq, "--ref", "Reference sequence to correct with", true, "reference");
  setUp.setOption(max_hp_run_gap_size, "--max_hp_run_gap_size", "the minimum size (inclusive) of the homopolymer run in which homopolymer gaps are corrected in");
  setUp.setOption(min_hp_run_size_to_correct, "--min_hp_run_size_to_correct", "the inserted or deleted homopolymer has to be this size(inclusive) or less to correct");
  setUp.processDefaultReader(true);

  setUp.finishSetUp(std::cout);

  uint64_t max_size = len(ref_seq);
  {

    SeqIO reader(setUp.pars_.ioOptions_);
    reader.openIn();
    seqInfo seq;
    while (reader.in_.readNextRead(seq)) {
      readVec::getMaxLength(seq, max_size);
    }
  }
  if (setUp.pars_.verbose_) {
    std::cout << "gap scoring:" << std::endl;
    std::cout << setUp.pars_.gapInfo_.toJson() << std::endl;
    std::cout << "count_end_gaps: " << setUp.pars_.colOpts_.alignOpts_.countEndGaps_ << std::endl;
  }
  aligner alignerObj(max_size, setUp.pars_.gapInfo_, setUp.pars_.scoring_, setUp.pars_.colOpts_.alignOpts_.countEndGaps_);
  alignerObj.processAlnInfoInput(setUp.pars_.alnInfoDirName_);

  SeqIO reader(setUp.pars_.ioOptions_);
  reader.openIn();
  reader.openOut();

  struct HpCorrection {
    HpCorrection(uint32_t start, uint32_t end, char base, bool delete_hp) : start_(start), end_(end),
      base_(base), delete_hp_(delete_hp) {
    }
    uint32_t start_{std::numeric_limits<uint32_t>::max()};
    uint32_t end_{std::numeric_limits<uint32_t>::max()};
    char base_{' '};
    bool delete_hp_{false}; //! whether or not to delete this postion
    uint32_t size() const {
      return end_ - start_;
    }
  };
  seqInfo seq;
  while (reader.in_.readNextRead(seq)) {
    // alignerObj.parts_.setMaxSize(ref_seq.seq_.size());
    alignerObj.alignCacheGlobal(ref_seq, seq);
    if (setUp.pars_.debug_) {
      alignerObj.alignObjectA_.seqBase_.outPutSeqAnsi(std::cout);
      alignerObj.alignObjectB_.seqBase_.outPutSeqAnsi(std::cout);
    }

    alignerObj.profileAlignment(ref_seq, seq, false, false, false);
    //first correct for homopolymers insertions
    std::vector<HpCorrection> corrections;
    for (const auto &gap: alignerObj.comp_.distances_.alignmentGaps_) {
      // std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
      // std::cout << "gap.second.gapedSequence_.size() <= 2: " <<  njh::colorBool(gap.second.gapedSequence_.size() <= 2)<< std::endl;

      if (gap.second.gapedSequence_.size() <= max_hp_run_gap_size) {
        //check if is homopolymer
        // std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
        // std::cout << "seqUtil::isHomopolyer(gap.second.gapedSequence_): " <<  njh::colorBool(seqUtil::isHomopolyer(gap.second.gapedSequence_))<< std::endl;

        if (seqUtil::isHomopolyer(gap.second.gapedSequence_)) {
          uint32_t size_of_ref_homopolymer = 0;
          uint32_t size_of_query_homopolymer = gap.second.gapedSequence_.size();
          //search backwards
          if (gap.first != 0) {
            //for now let's require that both homopolymers be the same location in both ref and query
            uint32_t cursor = gap.first;
            while (cursor > 0) {
              --cursor;
              if (alignerObj.alignObjectA_.seqBase_.seq_[cursor] == gap.second.gapedSequence_.front() &&
                  alignerObj.alignObjectB_.seqBase_.seq_[cursor] == gap.second.gapedSequence_.front()) {
                ++size_of_ref_homopolymer;
                ++size_of_query_homopolymer;
              } else {
                break;
              }
            }
          }
          if (gap.first + 1 != alignerObj.alignObjectB_.seqBase_.seq_.size()) {
            //for now let's require that both homopolymers be the same location in both ref and query
            uint32_t cursor = gap.first;
            while (cursor + 1 < alignerObj.alignObjectB_.seqBase_.seq_.size()) {
              ++cursor;
              if (alignerObj.alignObjectA_.seqBase_.seq_[cursor] == gap.second.gapedSequence_.front() &&
                  alignerObj.alignObjectB_.seqBase_.seq_[cursor] == gap.second.gapedSequence_.front()) {
                ++size_of_ref_homopolymer;
                ++size_of_query_homopolymer;
              } else {
                break;
              }
            }
          }
          if (size_of_query_homopolymer >= max_hp_run_gap_size && size_of_ref_homopolymer >= max_hp_run_gap_size) {
            corrections.emplace_back(gap.second.seqPos_, gap.second.seqPos_ + gap.second.gapedSequence_.size(),
                                     gap.second.gapedSequence_.front(), gap.second.ref_);
          }
        }
      }
    }
    for (const auto &cor: iter::reversed(corrections)) {
      if (setUp.pars_.verbose_) {
        std::cout << "cor.start_:" << cor.start_ << std::endl;
        std::cout << "cor.end_:" << cor.end_ << std::endl;
        std::cout << "cor.size():" << cor.size() << std::endl;
        std::cout << "cor.base_:" << cor.base_ << std::endl;
        std::cout << "cor.delete_hp_:" << njh::colorBool(cor.delete_hp_) << std::endl;
      }
      if (cor.delete_hp_) {
        seq.removeBases(cor.start_, cor.size());
      } else {
        // std::cout << "std::string(cor.base_, cor.size()): " << std::string(cor.size(), cor.base_) << std::endl;
        seq.insert(cor.start_, std::string(cor.size(), cor.base_));
      }
    }
    reader.write(seq);
    if (setUp.pars_.debug_) {
      alignerObj.alignCacheGlobal(ref_seq, seq);
      alignerObj.alignObjectA_.seqBase_.outPutSeqAnsi(std::cout);
      alignerObj.alignObjectB_.seqBase_.outPutSeqAnsi(std::cout);
      std::cout << std::endl;
    }
  }
  alignerObj.processAlnInfoOutput(setUp.pars_.outAlnInfoDirName_, setUp.pars_.verbose_);

	return 0;
}

int seqUtilsModRunner::correctHPRunsBasedOnSurroundingBaseCounts(const njh::progutils::CmdArgs & inputCommands) {

  // {
    // std::regex pat("(.{5})(T{4,})(.{5})"); // Matches 5 bases, 4+ T's, then 5 bases
    // std::smatch pat_match;
    // std::string seq = "GAATTTTTTAAATTTACTTTTTTAAATGAAGGAAAGTATTGTAAAG";
    //
    // std::string::const_iterator searchStart(seq.cbegin());
    //
    // while (std::regex_search(searchStart, seq.cend(), pat_match, pat)) {
    //   std::cout << "Match found: " << pat_match[0] << std::endl;
    //   std::cout << "Group 1: " << pat_match[1] << std::endl;
    //   std::cout << "Group 2: " << pat_match[2] << std::endl;
    //   std::cout << "Group 3: " << pat_match[3] << std::endl;
    //
    //   // Move iterator to continue searching beyond the current match
    //   std::cout << "(pat_match.prefix().second - seq.cbegin()) + pat_match.length(2): " << (pat_match.prefix().second - seq.cbegin()) + pat_match.length(2) << std::endl;
    //   searchStart = seq.cbegin() + (pat_match.prefix().second - seq.cbegin()) + pat_match.length(2);
    // }
    // std::regex pat("([ACGT]{4}[ACG])(T{4,})([ACG][AGCT]{4})");
    // std::smatch pat_match;
    // std::string seq = "GAATTTTTTAAATTTACTTTTTTAAATGAAGGAAAGTATTGTAAAG";
    //
    // std::string::const_iterator searchStart(seq.cbegin());
    //
    // while (std::regex_search(searchStart, seq.cend(), pat_match, pat)) {
    //   std::cout << "Match found: " << pat_match[0] << std::endl;
    //   std::cout << "Group 1: " << (pat_match[1].matched ? pat_match[1].str() : "N/A") << std::endl;
    //   std::cout << "Group 2: " << pat_match[2] << std::endl;
    //   std::cout << "Group 3: " << (pat_match[3].matched ? pat_match[3].str() : "N/A") << std::endl;
    //   std::cout << "--------------------------" << std::endl;
    //
    //   // Move iterator forward
    //   searchStart = seq.cbegin() + (pat_match.prefix().second - seq.cbegin()) + pat_match.length(2);
    // }
  //   return 0;
  // }


  uint32_t minReadCoverage = 12;
  double minPatternFreq = 0.3;
  OutOptions outOpts(bfs::path(""), ".tsv");
  OutOptions correcting_outOpts(bfs::path(""), ".tsv");

  uint32_t proceedingBases = 5;
  uint32_t trailingBases = 5;
  uint32_t minHomopolymerLength = 5;
  std::vector<char> homopolymerBases = {'A', 'C', 'G', 'T'};
  std::vector<char> allBases = {'A', 'C', 'G', 'T'};

  seqSetUp setUp(inputCommands);
  setUp.description_ = "correct homopolymers based on the pattern surrounding homopolymer runs";

  setUp.processVerbose();
  setUp.processDebug();
  setUp.setOption(correcting_outOpts.outFilename_, "--outCorrecting", "output file for what will be corrected");
  setUp.setOption(outOpts.outFilename_, "--outCounts", "output file for counts of the patterns");
  setUp.setOption(minReadCoverage, "--minReadCoverage", "min Read Coverage");
  setUp.setOption(minPatternFreq, "--minPatternFreq", "min Pattern Freq");

  setUp.setOption(proceedingBases, "--proceedingBases", "proceeding Bases");
  setUp.setOption(trailingBases, "--trailingBases", "trailing Bases");
  setUp.setOption(minHomopolymerLength, "--minHomopolymerLength", "min Homopolymer Length");
  setUp.setOption(homopolymerBases, "--homopolymerBases", "homopolymer Bases");
  setUp.setOption(allBases, "--allBases", "all Bases");

  setUp.processDefaultReader(true);
  outOpts.transferOverwriteOpts(setUp.pars_.ioOptions_.out_);
  correcting_outOpts.transferOverwriteOpts(setUp.pars_.ioOptions_.out_);
  setUp.finishSetUp(std::cout);

  SeqIO reader(setUp.pars_.ioOptions_);
  reader.openIn();
  reader.openOut();

  std::unique_ptr<OutputStream> countsOut;
  if (!outOpts.outFilename_.empty()) {
    countsOut = std::make_unique<OutputStream>(outOpts);
  }
  std::unordered_map<char, std::regex> basePatterns;
  std::string allBasesStr = njh::pasteAsStr(allBases);

  for (const auto base : homopolymerBases) {
    // std::string patStr = njh::pasteAsStr("(.{", proceedingBases, ",", proceedingBases,"})(", base, "{", minHomopolymerLength, ",})(.{", trailingBases, ",", trailingBases, "})");
    auto allBasesButHpBase = allBases;
    removeElement(allBasesButHpBase, base);
    std::string allBasesButHpBaseStr = njh::pasteAsStr(allBasesButHpBase);
    std::string patStr = njh::pasteAsStr("([", allBasesStr, "]", "{", proceedingBases - 1, ",", proceedingBases -1,"}","[", allBasesButHpBaseStr,"]",")(", base, "{", minHomopolymerLength, ",})(","[", allBasesButHpBaseStr,"]","[", allBasesStr, "]", "{", trailingBases - 1, ",", trailingBases - 1, "})");

    // std::cout << "patStr: " << patStr << std::endl;
    std::regex pattern(patStr);
    basePatterns.emplace(base, pattern);
  }

  std::unordered_map<char, std::regex> basePatternsForCorrecting;
  for (const auto base : homopolymerBases) {
    auto allBasesButHpBase = allBases;
    removeElement(allBasesButHpBase, base);
    std::string allBasesButHpBaseStr = njh::pasteAsStr(allBasesButHpBase);
    std::string patStr = njh::pasteAsStr("([", allBasesStr, "]", "{", proceedingBases - 1, ",", proceedingBases -1,"}","[", allBasesButHpBaseStr,"]",")(", base, "{", minHomopolymerLength - 1, ",})(","[", allBasesButHpBaseStr,"]","[", allBasesStr, "]", "{", trailingBases - 1, ",", trailingBases - 1, "})");

    // std::string patStr = njh::pasteAsStr("(.{", proceedingBases, ",", proceedingBases,"})(", base, "{", minHomopolymerLength - 1, ",})(.{", trailingBases, ",", trailingBases, "})");
    // std::cout << "patStr: " << patStr << std::endl;
    std::regex pattern(patStr);
    basePatternsForCorrecting.emplace(base, pattern);
  }
  std::map<char, std::unordered_map<std::string,std::unordered_map<std::string, uint32_t>>> patternCounts;
  {

    seqInfo seq;
    while(reader.readNextRead(seq)) {
      for (const auto base : homopolymerBases) {
        std::smatch match;
        std::string::const_iterator searchStart(seq.seq_.cbegin());
        while (std::regex_search(searchStart, seq.seq_.cend(), match, basePatterns.at(base))) {
          patternCounts[base][njh::pasteAsStr(match[1],"-",match[3])][match[2]]++;
          searchStart = seq.seq_.cbegin() + (match.prefix().second - seq.seq_.cbegin()) + match.length(2);
        }
      }
    }
  }



  if (!outOpts.outFilename_.empty()) {
    *countsOut << "base\tproceeding_trailing_bases\thomopolymer\tcount\tfreq\ttotal\tfull_pattern" << std::endl;
    for (const auto base : homopolymerBases) {
      if (njh::in(base, patternCounts)) {
        for (const auto & [pattern, patternCountsMap]  : patternCounts.at(base)) {
          if (patternCountsMap.size() > 1) {
            double total = 0;
            for (const auto & [pattern2, count] : patternCountsMap) {
              total += count;
            }
            for (const auto & [pattern2, count] : patternCountsMap) {
              *countsOut << base
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
  }

  //key1 = base, key2 = proceeding-trailing bases, value = best pattern to correct to
  std::map<char, std::unordered_map<std::string,std::string>> patternCorrections;

  for (const auto base : homopolymerBases) {
    if (njh::in(base, patternCounts)) {
      for (const auto & [pattern, patternCountsMap]  : patternCounts.at(base)) {
        if (patternCountsMap.size() > 1) {
          double total = 0;
          for (const auto & [pattern2, count] : patternCountsMap) {
            total += count;
          }
          double bestFreq = 0;
          std::string bestPattern;
          for (const auto & [pattern2, count] : patternCountsMap) {
            if (count / total > bestFreq) {
              bestFreq = count / total;
              bestPattern = njh::replaceString(pattern, "-", pattern2);
            }
          }
          if (bestFreq >= minPatternFreq && total >= minReadCoverage) {
            patternCorrections[base][pattern] = bestPattern;
          }
        }
      }
    }
  }
  if (!correcting_outOpts.outFilename_.empty()) {
    OutputStream correctingOut(correcting_outOpts);
    correctingOut << "base\tproceeding_trailing_bases\tcorrection" << std::endl;
    for (const auto & patternCorrection : patternCorrections) {
      for (const auto & [pattern, correction] : patternCorrection.second) {
        correctingOut << patternCorrection.first
            << "\t" << pattern
            << "\t" << correction
            << std::endl;
      }
    }
  }
  {
    //close and re-open
    reader.closeIn();
    reader.openIn();
    seqInfo seq;
    struct HPPatternReplacement {

      uint32_t full_pattern_start = std::numeric_limits<uint32_t>::max();
      uint32_t hp_start= std::numeric_limits<uint32_t>::max();
      uint32_t hp_len= std::numeric_limits<uint32_t>::max();
      std::string previous_full_pattern;
      std::string new_full_pattern;
      std::string surroundingPattern;
      uint8_t qual = 40;
      char base = ' ';
      [[nodiscard]] Json::Value toJson() const {
        Json::Value json;
        json["class"] = njh::json::toJson(njh::getTypeName(*this));
        json["full_pattern_start"] = njh::json::toJson(full_pattern_start);
        json["hp_start"] = njh::json::toJson(hp_start);
        json["hp_len"] = njh::json::toJson(hp_len);
        json["previous_full_pattern"] = njh::json::toJson(previous_full_pattern);
        json["new_full_pattern"] = njh::json::toJson(new_full_pattern);
        json["surroundingPattern"] = njh::json::toJson(surroundingPattern);
        json["qual"] = njh::json::toJson(qual);
        json["base"] = njh::json::toJson(base);

        return json;
      }
    };
    while(reader.readNextRead(seq)) {
      // bool print = "m84127_240426_214046_s2/242290354/ccs/11700_12792" == seq.name_;
      bool print = false;
      if (print) {
        std::cout << seq.name_ << std::endl;
      }
      std::vector<HPPatternReplacement> replacements;
      for (const auto base : homopolymerBases) {
        if (print) {
          std::cout <<  "\t" << base << std::endl;
        }

        std::smatch pat_match;
        std::string::const_iterator searchStart(seq.seq_.cbegin());
        while (std::regex_search(searchStart, seq.seq_.cend(), pat_match, basePatternsForCorrecting.at(base))) {
          auto surroundingPattern= njh::pasteAsStr(pat_match[1],"-",pat_match[3]);
          if (print) {
            std::cout << "\t\tcurrent_search: " << searchStart - seq.seq_.cbegin()  << std::endl;
          }
          searchStart = searchStart = seq.seq_.cbegin() + (pat_match.prefix().second - seq.seq_.cbegin()) + pat_match.length(2);
          if (print) {
            std::cout << "\t\tnext_search: " << searchStart - seq.seq_.cbegin()  << std::endl;
          }
          if (print) {
            std::cout << "\t\tsurroundingPattern: " << surroundingPattern << std::endl;
            std::cout << "\t\tnjh::in(surroundingPattern, patternCorrections[base]): " << njh::colorBool(njh::in(surroundingPattern, patternCorrections[base])) << std::endl;
            if (njh::in(surroundingPattern, patternCorrections[base])) {
              std::cout << "\t\tpat_match[0].str() != patternCorrections[base][surroundingPattern]: " << njh::colorBool(pat_match[0].str() != patternCorrections[base][surroundingPattern]) << std::endl;
            }

            // std::cout << "\t\tpat_match[0].str() != patternCorrections[base][surroundingPattern]: " << njh::colorBool(pat_match[0].str() != patternCorrections[base][surroundingPattern]) << std::endl;
          }
          if (njh::in(surroundingPattern, patternCorrections[base]) && pat_match[0].str() != patternCorrections[base][surroundingPattern]) {
            HPPatternReplacement replacement;
            replacement.surroundingPattern = surroundingPattern;
            replacement.previous_full_pattern = pat_match[0];
            replacement.new_full_pattern = patternCorrections[base][surroundingPattern];
            replacement.full_pattern_start = pat_match.prefix().second - seq.seq_.cbegin();
            replacement.hp_start = pat_match.prefix().second - seq.seq_.cbegin() + proceedingBases;
            replacement.hp_len = pat_match.length(2);
            replacement.qual = seq.qual_[replacement.hp_start + replacement.hp_len - 1];
            replacement.base = base;
            replacements.emplace_back(replacement);
            // if (pat_match[0] == "TAGCGTTTTTTTCCCCA") {
            //   std::cout << "\t\treplacements.size(): " << replacements.size() << std::endl;
            // }
            // if (replacements.empty()) {
            //   replacements.emplace_back(replacement);
            // } else {
            //
            // }
            // //since we are dealing with homopolymer and we are checking the patterns aren't the same, only two scenarios are the correction is longer or shorter, can't be equal
            // if (match[0].str().size() < patternCorrections[base][surroundingPattern].size()) {
            //   // std::cout << __FILE__ << " " << __LINE__ << std::endl;
            //   //adding bases
            //   auto diff = patternCorrections[base][surroundingPattern].size() - match[0].str().size();
            //   // auto hp_start = match.position(2);
            //   auto hp_start = match.prefix().second - seq.seq_.cbegin() + proceedingBases;
            //   auto hp_len = match.length(2);
            //   auto qual = seq.qual_[hp_start + hp_len - 1];
            //   // std::cout << "diff: " << diff << std::endl << " hp_start: " << hp_start << std::endl << " hp_len: " << hp_len << std::endl << " qual: " << static_cast<uint32_t>(qual) << std::endl;
            //   // std::cout << "seq.seq_.find(match[0]): " << seq.seq_.find(match[0].str()) << std::endl;
            //   // std::cout << "match.prefix().first - seq.seq_.cbegin(): " << match.prefix().first - seq.seq_.cbegin() << std::endl;
            //   // std::cout << "match.prefix().second - seq.seq_.cbegin(): " << match.prefix().second - seq.seq_.cbegin() << std::endl;
            //   // std::cout << "match.prefix().second - seq.seq_.cbegin(): " << match.prefix().second - seq.seq_.cbegin() + proceedingBases << std::endl;
            //   // std::cout << "seq.seq_.substr(match.prefix().second - seq.seq_.cbegin(), proceedingBases): " << seq.seq_.substr(match.prefix().second - seq.seq_.cbegin(), proceedingBases) << std::endl;
            //   // std::cout << "seq.seq_.substr(match.prefix().second - seq.seq_.cbegin() + proceedingBases, match.length(2)): " << seq.seq_.substr(match.prefix().second - seq.seq_.cbegin() + proceedingBases, match.length(2)) << std::endl;
            //   // std::cout << "seq.seq_.substr(match.prefix().second - seq.seq_.cbegin() + proceedingBases + match.length(2), trailingBases): " << seq.seq_.substr(match.prefix().second - seq.seq_.cbegin() + proceedingBases + match.length(2), trailingBases) << std::endl;
            //   // std::cout << "match.position(0): " << match.position(0) << std::endl;
            //   // std::cout << "match.position(1): " << match.position(1) << std::endl;
            //   // std::cout << "match.position(2): " << match.position(2) << std::endl;
            //   // std::cout << "match.position(3): " << match.position(3) << std::endl;
            //   // std::cout << "surroundingPattern: " << surroundingPattern << std::endl;
            //   // std::cout << "previous_full_pattern: " << match[0] << std::endl;
            //   // std::cout << "patternCorrections[base][surroundingPattern]: " << patternCorrections[base][surroundingPattern] << std::endl;
            //   // std::cout << "seq: " << seq.seq_ << std::endl;
            //   seq.insert(hp_start + hp_len, seqInfo("", std::string(diff,base), std::vector<uint8_t>(diff,qual)));
            //   // std::cout << "seq: " << seq.seq_ << std::endl;
            //   nextSearchStart += diff;
            //   // exit(1);
            // } else {
            //   //removing bases
            //   auto diff = match[0].str().size() - patternCorrections[base][surroundingPattern].size();
            //   // auto hp_start = match.position(2);
            //   auto hp_start = match.prefix().second - seq.seq_.cbegin() + proceedingBases;
            //   auto hp_len = match.length(2);
            //   std::cout << "diff: " << diff << std::endl << " hp_start: " << hp_start << std::endl << " hp_len: " << hp_len << std::endl;
            //   std::cout << "seq.seq_.find(match[0]): " << seq.seq_.find(match[0].str()) << std::endl;
            //   std::cout << "match.prefix().first - seq.seq_.cbegin(): " << match.prefix().first - seq.seq_.cbegin() << std::endl;
            //   std::cout << "match.prefix().second - seq.seq_.cbegin(): " << match.prefix().second - seq.seq_.cbegin() << std::endl;
            //   std::cout << "match.prefix().second - seq.seq_.cbegin(): " << match.prefix().second - seq.seq_.cbegin() + proceedingBases << std::endl;
            //   std::cout << "surroundingPattern: " << surroundingPattern << std::endl;
            //   std::cout << "previous_full_pattern: " << match[0] << std::endl;
            //
            //   std::cout << "patternCorrections[base][surroundingPattern]: " << patternCorrections[base][surroundingPattern] << std::endl;
            //   std::cout << "seq: " << seq.seq_ << std::endl;
            //   seq.removeBases(hp_start + hp_len - diff, diff);
            //   std::cout << "seq: " << seq.seq_ << std::endl;
            //   nextSearchStart -= diff;
            //   // exit(1);
            // }
          }
        }
        //iterative over backwards so positions don't get messed up
        if (print) {
          std::cout << "\treplacements.size(): " << replacements.size() << std::endl;
          for (const auto & replacement : iter::reversed(replacements)) {
            std::cout << "\t\t" << njh::json::writeAsOneLine(replacement.toJson()) << std::endl;
          }
        }
      }
      for (const auto& replacement: iter::reversed(replacements)) {
        //since we are dealing with homopolymer and we are checking the patterns aren't the same, only two scenarios are the correction is longer or shorter, can't be equal
        if (replacement.previous_full_pattern.size() < replacement.new_full_pattern.size()) {
          //adding bases
          auto diff = replacement.new_full_pattern.size() - replacement.previous_full_pattern.size();
          seq.insert(replacement.hp_start + replacement.hp_len,
                     seqInfo("", std::string(diff, replacement.base), std::vector<uint8_t>(diff, replacement.qual)));
        } else {
          //removing bases
          auto diff = replacement.previous_full_pattern.size() - replacement.new_full_pattern.size();
          seq.removeBases(replacement.hp_start + replacement.hp_len - diff, diff);
        }
      }

      if (print) {
        std::cout << std::endl;
      }
      reader.write(seq);
    }

  }
	return 0;
}


} //namespace njhseq

