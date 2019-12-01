//
// elucidator - A library for analyzing sequence data
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
//
// This file is part of elucidator.
//
// elucidator is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// elucidator is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with elucidator.  If not, see <http://www.gnu.org/licenses/>.
//

#include "seqToolsUtils.hpp"


namespace njhseq {


VecStr getPossibleDNASubstringsForProtein(const std::string& seq,
                                          const std::string& protein,
                                          const std::string& seqType) {

  VecStr ans;
  std::vector<std::unordered_map<size_t, size_t>> positions;
  std::map<char, VecStr> codonMap;
  if (seqType == "DNA") {
  	for(const auto & aa : aminoAcidInfo::infos::allInfo){
  		codonMap[aa.first] = aa.second.dnaCodons_;
  	}
  } else if (seqType == "RNA") {
  	for(const auto & aa : aminoAcidInfo::infos::allInfo){
  		codonMap[aa.first] = aa.second.rnaCodons_;
  	}
  } else {
  	std::stringstream ss;
    ss << "Unrecognized seqType : " << seqType << std::endl;
    throw std::runtime_error{ss.str()};
  }
  std::vector<size_t> firstCodonPositions;
  for (const auto& codon : codonMap.at(protein.front())) {
    addOtherVec(firstCodonPositions, findOccurences(seq, codon));
  }
  for (const auto& c : protein) {
    if (c == protein.front()) {
      continue;
    }
    std::vector<size_t> currentAminoAcid;
    for (const auto& codon : codonMap.at(c)) {
      addOtherVec(currentAminoAcid, findOccurences(seq, codon));
    }
    std::unordered_map<size_t, size_t> currentPositons;
    for (const auto& pos : currentAminoAcid) {
      currentPositons[pos] = pos;
    }
    positions.push_back(currentPositons);
  }
  std::vector<std::vector<size_t>> possiblePositions;
  for (const auto& firstPos : firstCodonPositions) {
    if (positions.back().find(firstPos + 3 * positions.size()) ==
        positions.back().end()) {
      continue;
    }
    std::vector<size_t> currentPositions;
    currentPositions.push_back(firstPos);
    bool madeItToTheEnd = true;
    size_t count = 1;
    for (const auto& codonPositions : positions) {
      if (codonPositions.find(firstPos + count * 3) != codonPositions.end()) {
        currentPositions.push_back(firstPos + count * 3);
      } else {
        madeItToTheEnd = false;
        break;
      }
      ++count;
    }

    if (madeItToTheEnd) {
      possiblePositions.push_back(currentPositions);
    }
  }
  for (const auto& possible : possiblePositions) {
    ans.emplace_back(getStringFromSubstrings(seq, possible, 3));
  }
  return ans;
}
VecStr findPossibleDNA(const std::string& seq, const std::string& protein,
                       const std::string& seqType, bool checkComplement) {
  VecStr ans = getPossibleDNASubstringsForProtein(seq, protein, seqType);
  if (checkComplement) {
    std::string complementSeq = seqUtil::reverseComplement(seq, "DNA");
    VecStr complementAns =
        getPossibleDNASubstringsForProtein(complementSeq, protein, seqType);
    for (auto& str : complementAns) {
      str = seqUtil::reverseComplement(str, "DNA");
    }
    addOtherVec(ans, complementAns);
  }
  return ans;
}
VecStr getAllCycloProteinFragments(const std::string& protein) {
  VecStr ans;
  for (auto i : iter::range(1, (int)protein.length())) {
    for (auto j : iter::range(protein.length())) {
      std::string currentFragment = "";
      if (j + i > protein.size()) {
        currentFragment = protein.substr(j, protein.size() - j) +
                          protein.substr(0, i - (protein.size() - j));
      } else {
        currentFragment = protein.substr(j, i);
      }
      ans.push_back(currentFragment);
    }
  }
  ans.push_back(protein);
  return ans;
}

VecStr getAllLinearProteinFragments(const std::string& protein) {
  VecStr ans;
  ans.push_back("");
  for (auto i : iter::range(1, (int)protein.length())) {
    for (auto j : iter::range(protein.length())) {
      std::string currentFragment = "";
      if (j + i > protein.size()) {

      } else {
        currentFragment = protein.substr(j, i);
      }
      ans.push_back(currentFragment);
    }
  }
  ans.push_back(protein);
  return ans;
}
std::multimap<int, std::string> getProteinFragmentSpectrum(
    const std::string& protein) {
  std::multimap<int, std::string, std::less<int>> ans;
  ans.insert({0, ""});
  VecStr fragments = getAllCycloProteinFragments(protein);
  for (const auto frag : fragments) {
    ans.insert({seqUtil::calculateWeightOfProteinInt(frag), frag});
  }
  return ans;
}

std::vector<int> getRealPossibleWeights(const std::vector<int>& spectrum) {
  std::vector<int> possibleWeights;
  for (auto i : spectrum) {
    if (i < 187) {
      possibleWeights.push_back(i);
    }
  }
  removeDuplicates(possibleWeights);
  // possibleWeights = getUnique(possibleWeights);
  std::vector<int> realPossibleWeights;
  for (const auto& pw : possibleWeights) {
    if (aminoAcidInfo::infos::weightIntToAminoAcid.find(pw) !=
        aminoAcidInfo::infos::weightIntToAminoAcid.end()) {
      realPossibleWeights.push_back(pw);
    }
  }
  return realPossibleWeights;
}


VecStr organizeLexicallyKmers(const VecStr& input, int colNum) {
  VecStr secondInputAsNumbers;

  MapStrStr stone;
  char count = 'A';
  for (const auto& sIter : input) {
    stone[std::string(1, count)] = sIter;
    secondInputAsNumbers.push_back(std::string(1, count));
    ++count;
  }
  uint64_t reserveSize = 0;
  for (auto i = 1; i <= colNum; ++i) {
    reserveSize += Factorial(i);
  }
  std::vector<VecStr> secondAns;
  secondAns.reserve(reserveSize);

  for (auto i = 1; i <= colNum; ++i) {
    auto temp = permuteVector(secondInputAsNumbers, i);
    secondAns.insert(secondAns.end(), temp.begin(), temp.end());
  }
  VecStr thirdAns;
  for (const auto& iter : secondAns) {
    thirdAns.push_back(vectorToString(iter, ""));
  }
  std::sort(thirdAns.begin(), thirdAns.end());
  for (auto& iter : thirdAns) {
    translateStringWithKey(iter, stone);
  }
  return thirdAns;
}
// doesn't actually work don't use
VecStr organizeLexicallyKmers(const std::string& input, size_t colNum) {
  std::string transformedString = "";

  MapStrStr stone;
  char count = 'A';
  for (const auto& c : input) {
    stone[std::string(1, count)] = c;
    transformedString.push_back(count);
    ++count;
  }
  uint64_t reserveSize = 0;
  for (size_t i = 1; i <= colNum; ++i) {
    reserveSize += Factorial(i);
  }
  VecStr secondAns;
  secondAns.reserve(reserveSize);
  VecStr fragments = getAllCycloProteinFragments(transformedString);
  VecStr trimmedFragments;
  for (const auto& frag : fragments) {
    if (frag.size() <= colNum) {
      trimmedFragments.push_back(frag);
    }
  }
  std::sort(trimmedFragments.begin(), trimmedFragments.end());
  printVector(trimmedFragments);
  for (const auto& frag : trimmedFragments) {
    auto temp = fastPermuteVectorOneLength(frag);
    printVector(temp);
    addOtherVec(secondAns, temp);
  }
  std::sort(secondAns.begin(), secondAns.end());
  for (auto& iter : secondAns) {
    translateStringWithKey(iter, stone);
  }
  return secondAns;
  /*
   exit(1);
   for (auto i = 1; i <= colNum; ++i) {
   auto temp = permuteVector(secondInputAsNumbers, i);
   secondAns.insert(secondAns.end(), temp.begin(), temp.end());
   }
   VecStr thirdAns;
   for (const auto& iter : secondAns) {
   thirdAns.push_back(vectorToString(iter, ""));
   }
   std::sort(thirdAns.begin(), thirdAns.end());

   return thirdAns;*/
}
uint64_t smallestSizePossible(uint64_t weight) {
  return std::floor(weight / 186.00);
}
int64_t getPossibleNumberOfProteins(
    int64_t weight, std::unordered_map<int64_t, int64_t>& cache) {
  if (0 == weight) {
    return 1;
  }
  if (0 > weight) {
    return 0;
  }
  int64_t count = 0;
  for (auto w : aminoAcidInfo::infos::weightIntToAminoAcid) {
    auto q = weight - w.first;
    // std::cout << q << std::endl;
    if (cache.find(q) == cache.end()) {
      auto t = getPossibleNumberOfProteins(q, cache);
      cache[q] = t;
    }
    count += cache[q];
  }
  return count;
}

probabilityProfile randomlyFindBestProfile(const VecStr& dnaStrings,
                                           const std::vector<VecStr>& allKmers,
                                           int numberOfKmers,
                                           njh::randomGenerator& gen) {
  VecStr randomMers;
  // bool needToSeed=true;

  std::vector<int> randomKmersNums =
      gen.unifRandVector(0, numberOfKmers, dnaStrings.size());
  uint32_t pos = 0;
  for (const auto& i : randomKmersNums) {
    randomMers.push_back(allKmers[pos][i]);
    ++pos;
  }
  probabilityProfile bestProfile(randomMers);
  while (true) {
    VecStr currentMers;
    for (const auto& dString : dnaStrings) {
      currentMers.push_back(bestProfile.mostProbableKmers(dString)[0].k_);
    }
    probabilityProfile currentProfile(currentMers);
    if (currentProfile.score_ < bestProfile.score_) {
      bestProfile = currentProfile;
    } else {
      return bestProfile;
    }
  }
  return bestProfile;
}

probabilityProfile randomMotifSearch(const VecStr& dnaStrings, int kLength,
                                     int numberOfRuns, bool gibs, int gibsRuns,
																		 njh::randomGenerator& gen) {
  std::vector<VecStr> allKmers;
  int numOfKmers = 0;
  for (const auto& dString : dnaStrings) {
    allKmers.push_back(kmerCalculator::getAllKmers(dString, kLength));
    numOfKmers = allKmers.back().size();
  }
  probabilityProfile bestProfile =
      randomlyFindBestProfile(dnaStrings, allKmers, numOfKmers, gen);
  for (int i = 0; i < numberOfRuns; ++i) {
    if (gibs) {
    	probabilityProfile currentProfile = randomlyFindBestProfileGibs(dnaStrings, allKmers,
                                                   numOfKmers, gibsRuns, gen);
      if (currentProfile.score_ < bestProfile.score_) {
        bestProfile = currentProfile;
      }
    } else {
    	probabilityProfile currentProfile =
          randomlyFindBestProfile(dnaStrings, allKmers, numOfKmers, gen);
      if (currentProfile.score_ < bestProfile.score_) {
        bestProfile = currentProfile;
      }
    }

  }
  return bestProfile;
}
probabilityProfile randomlyFindBestProfileGibs(
    const VecStr& dnaStrings, const std::vector<VecStr>& allKmers,
    int numberOfKmers, int runTimes, njh::randomGenerator& gen) {
  VecStr randomMers;
  std::vector<int> randomKmersNums =
      gen.unifRandVector(0, numberOfKmers, dnaStrings.size());
  uint32_t pos = 0;
  for (const auto& i : randomKmersNums) {
    randomMers.push_back(allKmers[pos][i]);
    ++pos;
  }
  probabilityProfile bestProfile(randomMers);
  // run random selection j times
  for (int j = 0; j < runTimes; ++j) {
    // pick random kmer to replace and get new profile without that kmer
    size_t randomStringNum =
        gen.unifRand(0, (int)bestProfile.dnaStrings_.size());
    VecStr currentMotifs;
    for (auto i : iter::range(bestProfile.dnaStrings_.size())) {
      if (i != randomStringNum) {
        currentMotifs.push_back(bestProfile.dnaStrings_[i]);
      }
    }
    probabilityProfile currentProfile(currentMotifs);
    // now select a new kmer by preferentially picking a more probable kmer but
    // still randomly
    std::multimap<double, std::string> kmersByProb;
    double cumProb = 0.0;
    for (const auto& ks : allKmers[randomStringNum]) {
      double currentProb = currentProfile.getProbabilityOfKmer(ks);
      kmersByProb.insert({currentProb, ks});
      cumProb += currentProb;
    }
    double randomProb = gen.unifRand(0.0, cumProb);
    double sumOfProbs = 0;
    std::string randomKmer = "";
    for (const auto& kByProb : kmersByProb) {
      sumOfProbs += kByProb.first;
      if (sumOfProbs > randomProb) {
        currentProfile.add(kByProb.second);
        currentProfile.updateScore();
        randomKmer = kByProb.second;
        break;
      }
    }
    // if new randomly replaced kmer creates a better profile repalce
    // best profile but keep order of strings so that they are correctly
    // replaced latter
    if (currentProfile.score_ < bestProfile.score_) {
      bestProfile.dnaStrings_[randomStringNum] = randomKmer;
      bestProfile = probabilityProfile(bestProfile.dnaStrings_);
    }
  }
  return bestProfile;
}
std::vector<std::vector<int>> growNextCycle(
    const std::vector<std::vector<int>>& previousCycle,
    const std::vector<int>& possibleWeights) {
  std::vector<std::vector<int>> nextCycle;
  nextCycle.reserve(previousCycle.size() * possibleWeights.size());
  for (const auto& pw : possibleWeights) {
    for (auto cycle : previousCycle) {
      cycle.push_back(pw);
      nextCycle.push_back(cycle);
    }
  }
  return nextCycle;
}

bool trimCycle(std::vector<std::vector<int>>& nextCycle,
               std::vector<std::vector<int>>& matchesSpectrum,
               const std::map<int, int>& spectrumToWeights,
               const std::vector<int> specVec) {
  int pos = (int)nextCycle.size();
  std::vector<int> positions;
  for (const auto& cycle : iter::reversed(nextCycle)) {
    --pos;
    auto currentFragments = getAllVecLinearFragments(cycle);
    std::map<int, int> fragmentWeightCounts;
    for (const auto& fragment : currentFragments) {
      auto currentSum = vectorSum(fragment);
      ++fragmentWeightCounts[currentSum];
    }
    bool faulted = false;
    for (const auto& weight : fragmentWeightCounts) {
      if (spectrumToWeights.find(weight.first) == spectrumToWeights.end() ||
          spectrumToWeights.at(weight.first) < weight.second) {
        nextCycle.erase(nextCycle.begin() + pos);
        positions.push_back(pos);
        faulted = true;
        break;
      }
    }
    if (faulted) {
      continue;
    }
    std::vector<int> theoSpectrum;
    theoSpectrum.reserve(currentFragments.size());
    auto currentCycloFragments = getAllVecCycloFragments(cycle);
    for (const auto& fragment : currentCycloFragments) {
      auto currentSum = vectorSum(fragment);
      theoSpectrum.push_back(currentSum);
    }
    std::sort(theoSpectrum.begin(), theoSpectrum.end());
    if (theoSpectrum == specVec) {
      matchesSpectrum.push_back(cycle);
      nextCycle.erase(nextCycle.begin() + pos);
    }
  }
  return !nextCycle.empty();
}

std::vector<std::vector<int>> getPossibleProteinsForSpectrum(
    const std::string& spectrum, bool useAllWeights) {
  std::vector<int> specVec = stringToVector<int>(spectrum);
  std::map<int, int> spectrumToWeights;
  for (const auto& i : specVec) {
    ++spectrumToWeights[i];
  }
  std::vector<int> possibleWeights;
  if (useAllWeights) {
    possibleWeights = getVectorOfMapKeys(aminoAcidInfo::infos::weightIntToAminoAcid);
    // for(const auto & weight : )
  } else {
    possibleWeights = getRealPossibleWeights(specVec);
  }

  std::vector<std::vector<int>> initialCycle;
  for (auto pw : possibleWeights) {
    initialCycle.emplace_back(std::vector<int>(1, pw));
  }
  bool keepGrowing = true;
  std::vector<std::vector<int>> matchesSpectrum;
  int count = 0;
  while (keepGrowing) {
    initialCycle = growNextCycle(initialCycle, possibleWeights);
    keepGrowing =
        trimCycle(initialCycle, matchesSpectrum, spectrumToWeights, specVec);
    std::sort(initialCycle.begin(), initialCycle.end());
    ++count;
  }
  return matchesSpectrum;
}

int scoreSpectrumAgreement(const std::map<int, int>& spectrum,
                           const std::map<int, int>& currentSpectrum) {
  int ans = 0;
  for (const auto& weight : currentSpectrum) {
    // std::cout<<"Weight: "<<weight.first<<" count:
    // "<<weight.second<<std::endl;
    if (spectrum.find(weight.first) == spectrum.end()) {

    } else if (spectrum.at(weight.first) <= weight.second) {
      ans += spectrum.at(weight.first);
      // std::cout<<"ansFirstOption: "<<ans<<std::endl;
    } else if (spectrum.at(weight.first) > weight.second) {
      ans += weight.second;
      // std::cout<<"ansSecondOption: "<<ans<<std::endl;
    } else {
      // this shouldn't happen
      std::cout << "ERROR!!" << std::endl;
    }
  }
  return ans;
}

std::multimap<int, std::vector<int>, std::greater<int>> growNextCycleScore(
    const std::multimap<int, std::vector<int>, std::greater<int>>&
        previousCycle,
    const std::vector<int>& possibleWeights,
    const std::map<int, int>& spectrumCounts, int parentMass, bool linear) {
  std::multimap<int, std::vector<int>, std::greater<int>> nextCycle;
  // nextCycle.reserve(previousCycle.size()*possibleWeights.size());
  for (const auto& pw : possibleWeights) {
    for (auto cycle : previousCycle) {
      cycle.second.push_back(pw);
      std::vector<std::vector<int>> currentFragments;
      if (linear) {
        currentFragments = getAllVecLinearFragments(cycle.second);
      } else {
        currentFragments = getAllVecCycloFragments(cycle.second);
      }
      std::map<int, int> fragmentWeightCounts;
      bool tooBig = false;
      // std::cout<<std::endl;
      for (const auto& fragment : currentFragments) {
        // std::cout<<"fragment: ";printVector(fragment);
        auto currentSum = vectorSum(fragment);
        // std::cout<<"Current sum: "<<currentSum<<std::endl;
        if (currentSum > parentMass) {
          tooBig = true;
        }
        ++fragmentWeightCounts[currentSum];
      }
      if (!tooBig) {
        int currentScore =
            scoreSpectrumAgreement(spectrumCounts, fragmentWeightCounts);
        // std::cout<<"currentScore: "<<currentScore<<std::endl;
        nextCycle.insert({currentScore, cycle.second});
      }
      // exit(1);
    }
  }
  return nextCycle;
}

std::multimap<int, std::vector<int>, std::greater<int>> trimCycleScore(
    std::multimap<int, std::vector<int>, std::greater<int>>& nextCycle,
    std::multimap<int, std::vector<int>, std::greater<int>>& matchesSpectrum,
    int parentMass, int leaderBoardNumber, int& currentLeader) {
  std::multimap<int, std::vector<int>, std::greater<int>> ans;
  int pos = 0;
  int lastScore;
  for (const auto& cycle : nextCycle) {
    if (pos == 0) {
      lastScore = cycle.first;
    }
    ++pos;
    if (pos > leaderBoardNumber && cycle.first != lastScore) {

    } else {
      int mass = vectorSum(cycle.second);
      if (mass == parentMass && cycle.first == currentLeader) {
        matchesSpectrum.insert(cycle);
      } else if (mass == parentMass && cycle.first >= currentLeader) {
        matchesSpectrum.clear();
        matchesSpectrum.insert(cycle);
        currentLeader = cycle.first;
      }
      lastScore = cycle.first;
      ans.insert(cycle);
    }
  }
  return ans;
}
std::multimap<int, std::vector<int>, std::greater<int>>
getPossibleProteinsForSpectrum(const std::string& spectrum,
                               int leaderBoardNumber, bool verbose,
                               bool useAllWeights, bool convolution,
                               int convolutionCutOff, bool linear) {
  std::vector<int> specVec = stringToVector<int>(spectrum);
  std::sort(specVec.begin(), specVec.end());
  int parentMass = specVec.back();
  std::map<int, int> spectrumToWeights;
  for (const auto& i : specVec) {
    ++spectrumToWeights[i];
  }
  /*
   for(const auto & specW : spectrumToWeights){
   std::cout<<specW.first<< ":"<<specW.second<<std::endl;
   }*/
  std::vector<int> possibleWeights;
  if (convolution) {
    possibleWeights = topConvolutionWeights(specVec, convolutionCutOff);
    std::cout << "convolutionCutOff: " << convolutionCutOff << std::endl;
    std::cout << "convolutionPossibleWeights: ";
    printVector(possibleWeights);
  } else if (useAllWeights) {
    possibleWeights = std::vector<int>(144);
    std::iota(possibleWeights.begin(), possibleWeights.end(), 57);
  } else {
    possibleWeights = getVectorOfMapKeys(aminoAcidInfo::infos::weightIntToAminoAcid);
  }

  std::multimap<int, std::vector<int>, std::greater<int>> initialCycle;
  for (auto pw : possibleWeights) {
    initialCycle.insert({0, std::vector<int>(1, pw)});
  }
  bool keepGrowing = true;
  std::multimap<int, std::vector<int>, std::greater<int>> matchesSpectrum;
  int currentLeader = 0;
  int count = 0;
  while (keepGrowing) {
    std::cout << "on cycle " << count << std::endl;
    initialCycle = growNextCycleScore(initialCycle, possibleWeights,
                                      spectrumToWeights, parentMass, linear);
    if (verbose) {
      for (const auto& cycle : initialCycle) {
        std::cout << cycle.first << " " << vectorToString(cycle.second)
                  << std::endl;
      }
    }
    initialCycle = trimCycleScore(initialCycle, matchesSpectrum, parentMass,
                                  leaderBoardNumber, currentLeader);
    if (initialCycle.empty()) {
      keepGrowing = false;
    }
    /*if (initialCycle.empty() || !matchesSpectrum.empty()) {
      keepGrowing=false;
    }*/
    ++count;
  }
  return matchesSpectrum;
}
std::map<int, int> getConvolutionWeights(std::vector<int> experimentalSpectrum,
                                         int multiplicityCutOff, int lowerBound,
                                         int upperBound) {
  std::sort(experimentalSpectrum.begin(), experimentalSpectrum.end());
  std::map<int, int, std::greater<int>> countsOfDifferences;
  for (const auto& i : iter::range(experimentalSpectrum.size())) {
    for (const auto& j : iter::range(i + 1, experimentalSpectrum.size())) {
      int currentDifference = experimentalSpectrum[j] - experimentalSpectrum[i];
      if (currentDifference >= lowerBound && currentDifference <= upperBound) {
        ++countsOfDifferences[currentDifference];
      }
    }
  }
  std::map<int, int> ans;
  for (const auto& count : countsOfDifferences) {
    if (count.second >= multiplicityCutOff) {
      ans.insert(count);
    }
  }
  return ans;
}
std::vector<int> convolutionWeights(std::vector<int> experimentalSpectrum,
                                    int multiplicityCutOff, int lowerBound,
                                    int upperBound) {
  std::map<int, int> counts = getConvolutionWeights(
      experimentalSpectrum, multiplicityCutOff, lowerBound, upperBound);
  return getVectorOfMapKeys(counts);
}
std::vector<int> topConvolutionWeights(std::vector<int> experimentalSpectrum,
                                       int mItems, int lowerBound,
                                       int upperBound) {
  std::map<int, int> allCounts =
      getConvolutionWeights(experimentalSpectrum, 1, lowerBound, upperBound);
  std::multimap<int, int, std::greater<int>> byCounts;
  for (const auto& count : allCounts) {
    byCounts.insert({count.second, count.first});
  }
  int counting = 0;
  int lastCount = byCounts.begin()->first;
  // int testCount=0;
  /*for(const auto & count : byCounts){
    testCount++;
    std::cout<<"testCount: "<<testCount<<" currentCount: "<<count.first<<"
  weight: "<<count.second<<std::endl;
  }*/
  for (const auto& count : byCounts) {
    ++counting;
    if (counting > mItems && count.first != lastCount) {
      break;
    }
    lastCount = count.first;
  }
  // std::cout<<"lastcount: "<<lastCount<<std::endl;
  return convolutionWeights(experimentalSpectrum, lastCount, lowerBound,
                            upperBound);
}
int64_t getMinCoins(int64_t change, const std::vector<int64_t>& coins,
                    std::unordered_map<int64_t, int64_t>& cache) {
  if (0 == change) {
    return 0;
  }
  int64_t count = change;
  for (auto coin : coins) {
    if (change - coin < 0) {

    } else {
      int64_t q = change - coin;
      if (cache.find(q) == cache.end()) {
        auto t = getMinCoins(q, coins, cache) + 1;
        cache[q] = t;
      }
      if (cache[q] < count) {
        count = cache[q];
      }
    }
  }
  return count;
}






std::vector<uint32_t> getWindowQuals(const std::vector<uint32_t>& qual,
                                     uint32_t qWindowSize, uint32_t pos) {
  uint32_t lowerBound = 0;
  uint32_t higherBound = qual.size();
  if (pos > qWindowSize) {
    lowerBound = pos - qWindowSize;
  }
  if (pos + qWindowSize + 1 < higherBound) {
    higherBound = pos + qWindowSize + 1;
  }
  return getSubVector(qual, lowerBound, higherBound - lowerBound);
}
std::vector<double> likelihoodForBaseQ(
    const std::vector<uint32_t>& qual,
    std::unordered_map<double, double>& likelihoods) {
  std::vector<double> ans;
  for (const auto& i : iter::range(qual.size())) {
    ans.emplace_back(likelihoods.at(qual[i]));
  }
  return ans;
}
std::vector<double> likelihoodForMeanQ(
    const std::vector<uint32_t>& qual, uint32_t qWindowSize,
    std::unordered_map<double, double>& likelihoods) {
  std::vector<double> ans;
  for (const auto& i : iter::range(qual.size())) {
    double currentMean = vectorMean(getWindowQuals(qual, qWindowSize, i));
    ans.emplace_back(likelihoods.at(roundDecPlaces(currentMean, 2)));
  }
  return ans;
}
std::vector<double> likelihoodForMedianQ(
    const std::vector<uint32_t>& qual, uint32_t qWindowSize,
    std::unordered_map<double, double>& likelihoods) {
  std::vector<double> ans;
  for (const auto& i : iter::range(qual.size())) {
    double currentMedian = vectorMedianCopy(getWindowQuals(qual, qWindowSize, i));
    ans.emplace_back(likelihoods.at(roundDecPlaces(currentMedian, 2)));
  }
  return ans;
}
std::vector<double> likelihoodForMinQ(
    const std::vector<uint32_t>& qual, uint32_t qWindowSize,
    std::unordered_map<double, double>& likelihoods) {
  std::vector<double> ans;
  for (const auto& i : iter::range(qual.size())) {
    double currentMin = vectorMinimum(getWindowQuals(qual, qWindowSize, i));
    ans.emplace_back(likelihoods.at(roundDecPlaces(currentMin, 2)));
  }
  return ans;
}
double getChangeInHydro(const char& firstAA, const char& secondAA) {
  return std::abs(aminoAcidInfo::infos::allInfo.at(firstAA).acidHydrophobicity_ -
                  aminoAcidInfo::infos::allInfo.at(secondAA).acidHydrophobicity_);
}

/*
std::vector<double> getHydroChanges(const std::string& originalCodon,
                                    const VecStr& mutantCodons,
                                    const std::map<std::string, char>& code) {
  std::vector<double> ans;
  if (originalCodon.size() != 3) {
    std::cout << "codon needs to be 3 bases" << std::endl;
    std::cout << originalCodon << std::endl;
    std::cout << originalCodon.size() << std::endl;
    return ans;
  }
  char originalAA = code.at(originalCodon);
  if (originalAA == '*') {
    return ans;
  }
  for (const auto& codon : mutantCodons) {
    char mutantAA = code.at(codon);
    if (mutantAA == '*') {
      continue;
    }
    ans.emplace_back(getChangeInHydro(originalAA, mutantAA));
  }
  return ans;
}*/







void simpleCollapseQueryCov(std::vector<cluster> & consensusReads, aligner & alignerObj,
		const comparison & passableErrors) {
	if(consensusReads.size() <=1){
		return;
	}
	readVecSorter::sort(consensusReads);
	//compare final consensus reads to see if any created very similar consensus and should be collapsed into one
	std::set<uint64_t> needToRemove;
	comparison noErrors;
	uint32_t klenComp = 9;
	double kdistCutOff = .95;
	std::vector<kmerInfo> kinfos;
	for(const auto & seq : consensusReads){
		kinfos.emplace_back(kmerInfo(seq.seqBase_.seq_, klenComp, false));
	}
	//first do a comparison of no errors to collapse almost identical clusters
	for (const auto & firstPos : iter::range(consensusReads.size())) {
		if (consensusReads[firstPos].remove) {
			continue;
		}
		for (const auto & secondPos : iter::range(firstPos + 1,
				consensusReads.size())) {
			if (consensusReads[secondPos].remove) {
				continue;
			}
			//first check for an exact sequence match
			if (consensusReads[secondPos].seqBase_.seq_
					== consensusReads[firstPos].seqBase_.seq_) {
				consensusReads[firstPos].addRead(consensusReads[secondPos]);
				needToRemove.emplace(secondPos);
				consensusReads[secondPos].remove = true;
			} else {
				if(kinfos[firstPos].compareKmers(kinfos[secondPos]).second < kdistCutOff){
					continue;
				}
				//check to see if there are no errors
				//(basically if there is a little bit of overlap, could be dangerous for very different lengths)
				alignerObj.alignCache(consensusReads[secondPos],
						consensusReads[firstPos], false);
				alignerObj.profileAlignment(consensusReads[secondPos],
						consensusReads[firstPos], false, true, false);
				double queryCoverage = passableErrors.distances_.query_.coverage_ != 0 ?  passableErrors.distances_.query_.coverage_: 1;
				if (noErrors.passErrorProfile(alignerObj.comp_)
						&& alignerObj.comp_.distances_.query_.coverage_ > queryCoverage) {
					consensusReads[firstPos].addRead(consensusReads[secondPos]);
					needToRemove.emplace(secondPos);
					consensusReads[secondPos].remove = true;
				}
			}
		}
	}
	//now do a second compare with some errors allowed
	std::vector<uint32_t> positions(consensusReads.size(), 0);
	njh::iota<uint32_t>(positions, 0);
	for (const auto & firstPos : iter::reversed(positions)) {
		if (consensusReads[firstPos].remove) {
			continue;
		}
		for (const auto & secondPos : iter::range(firstPos)) {
			if (consensusReads[secondPos].remove) {
				continue;
			}
			if(kinfos[firstPos].compareKmers(kinfos[secondPos]).second < kdistCutOff){
				continue;
			}
			alignerObj.alignCache(consensusReads[secondPos], consensusReads[firstPos],
					false);
			alignerObj.profileAlignment(consensusReads[secondPos],
					consensusReads[firstPos], false, true, false);
			double queryCoverage = passableErrors.distances_.query_.coverage_ != 0 ?  passableErrors.distances_.query_.coverage_: 1;

			if (passableErrors.passErrorProfile(alignerObj.comp_) &&
					alignerObj.comp_.distances_.query_.coverage_ >= queryCoverage) {
				consensusReads[firstPos].addRead(consensusReads[secondPos]);
				needToRemove.emplace(secondPos);
				consensusReads[secondPos].remove = true;
			}
		}
	}
	//remove the reads that need to be removed
	for (const auto & pos : iter::reversed(needToRemove)) {
		consensusReads.erase(consensusReads.begin() + pos);
	}
	//calculate consensus again
	bool changed = true;
	uint32_t calcCount = 0;
	while(changed){
		++calcCount;
		VecStr oldSeqs = readVec::getSeqs(consensusReads);
		clusterVec::allCalculateConsensus(consensusReads, alignerObj, true);
		readVec::allUpdateName(consensusReads);
		changed = false;
		for(const auto & seqPos : iter::range(consensusReads.size())){
			if(oldSeqs[seqPos] != consensusReads[seqPos].seqBase_.seq_){
				changed = true;
				break;
			}
		}
	}
	//if consensus changed collapse again
	if(calcCount > 1){
		simpleCollapse(consensusReads, alignerObj, passableErrors);
	}
}

void simpleCollapse(std::vector<cluster> & consensusReads, aligner & alignerObj,
		const comparison & passableErrors) {
	if(consensusReads.size() <=1){
		return;
	}
	readVecSorter::sort(consensusReads);
	//compare final consensus reads to see if any created very similar consensus and should be collapsed into one
	std::set<uint64_t> needToRemove;
	comparison noErrors;
	uint32_t klenComp = 9;
	double kdistCutOff = .95;
	std::vector<kmerInfo> kinfos;
	for(const auto & seq : consensusReads){
		kinfos.emplace_back(kmerInfo(seq.seqBase_.seq_, klenComp, false));
	}
	//first do a comparison of no errors to collapse almost identical clusters
	for (const auto & firstPos : iter::range(consensusReads.size())) {
		if (consensusReads[firstPos].remove) {
			continue;
		}
		for (const auto & secondPos : iter::range(firstPos + 1,
				consensusReads.size())) {
			if (consensusReads[secondPos].remove) {
				continue;
			}
			//first check for an exact sequence match
			if (consensusReads[secondPos].seqBase_.seq_
					== consensusReads[firstPos].seqBase_.seq_) {
				consensusReads[firstPos].addRead(consensusReads[secondPos]);
				needToRemove.emplace(secondPos);
				consensusReads[secondPos].remove = true;
			} else {
				if(kinfos[firstPos].compareKmers(kinfos[secondPos]).second < kdistCutOff){
					continue;
				}
				//check to see if there are no errors
				//(basically if there is a little bit of overlap, could be dangerous for very different lengths)
				alignerObj.alignCache(consensusReads[secondPos],
						consensusReads[firstPos], false);
				alignerObj.profileAlignment(consensusReads[secondPos],
						consensusReads[firstPos], false, true, false);
				if (noErrors.passErrorProfile(alignerObj.comp_)) {
					consensusReads[firstPos].addRead(consensusReads[secondPos]);
					needToRemove.emplace(secondPos);
					consensusReads[secondPos].remove = true;
				}
			}
		}
	}
	//now do a second compare with some errors allowed
	//decreacse kDist cut off
	kdistCutOff = .80;
	std::vector<uint32_t> positions(consensusReads.size(), 0);
	njh::iota<uint32_t>(positions, 0);
	for (const auto & firstPos : iter::reversed(positions)) {
		if (consensusReads[firstPos].remove) {
			continue;
		}
		for (const auto & secondPos : iter::range(firstPos)) {
			if (consensusReads[secondPos].remove) {
				continue;
			}
			if(kinfos[firstPos].compareKmers(kinfos[secondPos]).second < kdistCutOff){
				continue;
			}
			alignerObj.alignCache(consensusReads[secondPos], consensusReads[firstPos],
					false);
			alignerObj.profileAlignment(consensusReads[secondPos],
					consensusReads[firstPos], false, true, false);
			if (passableErrors.passErrorProfile(alignerObj.comp_)) {
				consensusReads[firstPos].addRead(consensusReads[secondPos]);
				needToRemove.emplace(secondPos);
				consensusReads[secondPos].remove = true;
			}
		}
	}
	//remove the reads that need to be removed
	for (const auto & pos : iter::reversed(needToRemove)) {
		consensusReads.erase(consensusReads.begin() + pos);
	}
	//calculate consensus again
	bool changed = true;
	uint32_t calcCount = 0;
	while(changed){
		++calcCount;
		VecStr oldSeqs = readVec::getSeqs(consensusReads);
		clusterVec::allCalculateConsensus(consensusReads, alignerObj, true);
		readVec::allUpdateName(consensusReads);
		changed = false;
		for(const auto & seqPos : iter::range(consensusReads.size())){
			if(oldSeqs[seqPos] != consensusReads[seqPos].seqBase_.seq_){
				changed = true;
				break;
			}
		}
	}
	//if consensus changed collapse again
	if(calcCount > 1){
		simpleCollapse(consensusReads, alignerObj, passableErrors);
	}
}


}  // namespace njh
