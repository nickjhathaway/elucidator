//
// Created by Nicholas Hathaway on 2/2/25.
//
#include "seqUtilsModRunner.hpp"
#include <njhseq/objects/counters/DNABaseCounter.hpp>
#include <njhseq/IO/SeqIO/SeqIO.hpp>
#include <njhseq/objects/seqObjects/seqKmers.h>
#include <njhseq/readVectorManipulation/readVectorHelpers/readVecSorter.hpp>


namespace njhseq {

int seqUtilsModRunner::correctHPRunsBasedOnSurroundingBaseCounts(const njh::progutils::CmdArgs & inputCommands) {
  uint32_t minReadCoverage = 50;
  double minPatternFreq = 0.7;
  OutOptions outOpts(bfs::path(""), ".tsv");
  uint32_t proceedingBases = 5;
  uint32_t trailingBases = 5;
  uint32_t minHomopolymerLength = 5;
  std::vector<char> homopolymerBases = {'A', 'C', 'G', 'T'};
  seqSetUp setUp(inputCommands);
  setUp.description_ = "count the pattern surrounding homopolymer runs";

  setUp.processVerbose();
  setUp.processDebug();
  setUp.setOption(outOpts.outFilename_, "--outCounts", "output file for counts of the patterns");
  setUp.setOption(minReadCoverage, "--minReadCoverage", "min Read Coverage");
  setUp.setOption(minPatternFreq, "--minPatternFreq", "min Pattern Freq");

  setUp.setOption(proceedingBases, "--proceedingBases", "proceeding Bases");
  setUp.setOption(trailingBases, "--trailingBases", "trailing Bases");
  setUp.setOption(minHomopolymerLength, "--minHomopolymerLength", "min Homopolymer Length");
  setUp.setOption(homopolymerBases, "--homopolymerBases", "homopolymer Bases");

  setUp.processDefaultReader(true);
  outOpts.transferOverwriteOpts(setUp.pars_.ioOptions_.out_);
  setUp.finishSetUp(std::cout);

  SeqIO reader(setUp.pars_.ioOptions_);
  reader.openIn();
  reader.openOut();

  std::unique_ptr<OutputStream> countsOut;
  if (!outOpts.outFilename_.empty()) {
    countsOut = std::make_unique<OutputStream>(outOpts);
  }
  std::unordered_map<char, std::regex> basePatterns;
  for (const auto base : homopolymerBases) {
    std::string patStr = njh::pasteAsStr("(.{", proceedingBases, ",", proceedingBases,"})(", base, "{", minHomopolymerLength, ",})(.{", trailingBases, ",", trailingBases, "})");
    // std::cout << "patStr: " << patStr << std::endl;
    std::regex pattern(patStr);
    basePatterns.emplace(base, pattern);
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
          searchStart = match.suffix().first;
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

    };
    while(reader.readNextRead(seq)) {
      for (const auto base : homopolymerBases) {
        std::vector<HPPatternReplacement> replacements;
        std::smatch pat_match;
        std::string::const_iterator searchStart(seq.seq_.cbegin());
        while (std::regex_search(searchStart, seq.seq_.cend(), pat_match, basePatterns.at(base))) {
          auto surroundingPattern= njh::pasteAsStr(pat_match[1],"-",pat_match[3]);
          auto nextSearchStart = pat_match.suffix().first;
          if (njh::in(surroundingPattern, patternCorrections[base]) && pat_match[0].str() != patternCorrections[base][surroundingPattern]) {
            HPPatternReplacement replacement;
            replacement.surroundingPattern = surroundingPattern;
            replacement.previous_full_pattern = pat_match[0];
            replacement.new_full_pattern = patternCorrections[base][surroundingPattern];
            replacement.full_pattern_start = pat_match.prefix().second - seq.seq_.cbegin();
            replacement.hp_start = pat_match.prefix().second - seq.seq_.cbegin() + proceedingBases;
            replacement.hp_len = pat_match.length(2);
            replacement.qual = seq.qual_[replacement.hp_start + replacement.hp_len - 1];
            replacements.emplace_back(replacement);

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
          searchStart = nextSearchStart;
        }
        //iterative over backwards so positions don't get messed up
        for (const auto & replacement : iter::reversed(replacements)) {
            //since we are dealing with homopolymer and we are checking the patterns aren't the same, only two scenarios are the correction is longer or shorter, can't be equal
            if (replacement.previous_full_pattern.size() < replacement.new_full_pattern.size()) {
              // std::cout << __FILE__ << " " << __LINE__ << std::endl;
              //adding bases
              auto diff = replacement.new_full_pattern.size() - replacement.previous_full_pattern.size();
              // auto hp_start = match.position(2);
              // auto hp_start = match.prefix().second - seq.seq_.cbegin() + proceedingBases;
              // auto hp_len = match.length(2);
              // auto qual = seq.qual_[hp_start + hp_len - 1];
              // std::cout << "diff: " << diff << std::endl << " hp_start: " << hp_start << std::endl << " hp_len: " << hp_len << std::endl << " qual: " << static_cast<uint32_t>(qual) << std::endl;
              // std::cout << "seq.seq_.find(match[0]): " << seq.seq_.find(replacement.previous_full_pattern) << std::endl;
              // std::cout << "match.prefix().first - seq.seq_.cbegin(): " << match.prefix().first - seq.seq_.cbegin() << std::endl;
              // std::cout << "match.prefix().second - seq.seq_.cbegin(): " << match.prefix().second - seq.seq_.cbegin() << std::endl;
              // std::cout << "match.prefix().second - seq.seq_.cbegin(): " << match.prefix().second - seq.seq_.cbegin() + proceedingBases << std::endl;
              // std::cout << "seq.seq_.substr(match.prefix().second - seq.seq_.cbegin(), proceedingBases): " << seq.seq_.substr(match.prefix().second - seq.seq_.cbegin(), proceedingBases) << std::endl;
              // std::cout << "seq.seq_.substr(match.prefix().second - seq.seq_.cbegin() + proceedingBases, match.length(2)): " << seq.seq_.substr(match.prefix().second - seq.seq_.cbegin() + proceedingBases, match.length(2)) << std::endl;
              // std::cout << "seq.seq_.substr(match.prefix().second - seq.seq_.cbegin() + proceedingBases + match.length(2), trailingBases): " << seq.seq_.substr(match.prefix().second - seq.seq_.cbegin() + proceedingBases + match.length(2), trailingBases) << std::endl;
              // std::cout << "match.position(0): " << match.position(0) << std::endl;
              // std::cout << "match.position(1): " << match.position(1) << std::endl;
              // std::cout << "match.position(2): " << match.position(2) << std::endl;
              // std::cout << "match.position(3): " << match.position(3) << std::endl;
              // std::cout << "surroundingPattern: " << surroundingPattern << std::endl;
              // std::cout << "previous_full_pattern: " << match[0] << std::endl;
              // std::cout << "replacement.new_full_pattern: " << replacement.new_full_pattern << std::endl;
              // std::cout << "seq: " << seq.seq_ << std::endl;
              seq.insert(replacement.hp_start + replacement.hp_len, seqInfo("", std::string(diff,base), std::vector<uint8_t>(diff,replacement.qual)));
              // std::cout << "seq: " << seq.seq_ << std::endl;
              // exit(1);
            } else {
              //removing bases
              auto diff = replacement.previous_full_pattern.size() - replacement.new_full_pattern.size();
              // auto hp_start = match.position(2);
              // auto hp_start = match.prefix().second - seq.seq_.cbegin() + proceedingBases;
              // auto hp_len = match.length(2);
              // std::cout << "diff: " << diff << std::endl << " hp_start: " << replacement.hp_start << std::endl << " hp_len: " << replacement.hp_len << std::endl;
              // std::cout << "seq.seq_.find(match[0]): " << seq.seq_.find(replacement.previous_full_pattern) << std::endl;
              // std::cout << "match.prefix().first - seq.seq_.cbegin(): " << match.prefix().first - seq.seq_.cbegin() << std::endl;
              // std::cout << "match.prefix().second - seq.seq_.cbegin(): " << match.prefix().second - seq.seq_.cbegin() << std::endl;
              // std::cout << "match.prefix().second - seq.seq_.cbegin(): " << match.prefix().second - seq.seq_.cbegin() + proceedingBases << std::endl;
              // std::cout << "surroundingPattern: " << surroundingPattern << std::endl;
              // std::cout << "previous_full_pattern: " << replacement.previous_full_pattern << std::endl;
              // std::cout << "replacement.new_full_pattern: " << replacement.new_full_pattern << std::endl;
              // std::cout << "seq: " << seq.seq_ << std::endl;
              seq.removeBases(replacement.hp_start + replacement.hp_len - diff, diff);
              // std::cout << "seq: " << seq.seq_ << std::endl;
            }
        }

      }
      reader.write(seq);
    }

  }
	return 0;
}


} //namespace njhseq

