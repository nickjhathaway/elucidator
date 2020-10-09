#pragma once
//
//  seqToolsUtils.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 4/27/13.
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



#include "elucidator/common.h"
#include <njhseq/objects/BioDataObject/GenomicRegion.hpp>

#include <TwoBit.h>

namespace njhseq {


template<typename SEQTYPE, typename REFTYPE>
std::vector<GenomicRegion> getOverlapingRegions(
		const std::vector<SEQTYPE> & inSeqs, const REFTYPE & refSeq,
		aligner & alignerObj, uint32_t minlen, bool onlyOverlapping) {
	std::vector<GenomicRegion> regions;
	for(const auto & seq : inSeqs){
		alignerObj.alignCacheGlobal(refSeq, seq);
		alignerObj.rearrangeObjsGlobal(refSeq, seq);
		auto refAlnStart =   alignerObj.alignObjectA_.seqBase_.seq_.find_first_not_of("-");
		auto refAlnStop =    alignerObj.alignObjectA_.seqBase_.seq_.find_last_not_of("-");
		auto queryAlnStart = alignerObj.alignObjectB_.seqBase_.seq_.find_first_not_of("-");
		auto queryAlnEnd =   alignerObj.alignObjectB_.seqBase_.seq_.find_last_not_of("-");
		uint32_t start = 0;
		uint32_t end = getSeqBase(refSeq).seq_.size();
		if(queryAlnStart > refAlnStart){
			start =  alignerObj.getSeqPosForAlnAPos(queryAlnStart);
		}
		if(queryAlnEnd < refAlnStop){
			end =  alignerObj.getSeqPosForAlnAPos(queryAlnEnd) + 1;
		}
		regions.emplace_back(getSeqBase(seq).name_, getSeqBase(refSeq).name_, start, end, false);
	}

	sortGRegionsByStart(regions);
	std::vector<GenomicRegion> uniqueRegions = regions;
//	std::vector<GenomicRegion> uniqueRegions;
//	for (const auto & reg : regions) {
//		bool pass = true;
//		for (const auto & uniReg : uniqueRegions) {
//			if (uniReg.sameRegion(reg)) {
//				pass = false;
//				break;
//			}
//		}
//		if (pass) {
//			uniqueRegions.emplace_back(reg);
//		}
//	}

	std::vector<GenomicRegion> trimRegions;
	std::vector<GenomicRegion> overlappingRegions;
	//add in any regions that don't overlap other regions
	for(const auto & uniReg : uniqueRegions){
		bool found = false;
		for(const auto & otherReg : uniqueRegions){
			if(otherReg.uid_ != uniReg.uid_ &&  otherReg.overlaps(uniReg, minlen)){
				found = true;
				break;
			}
		}
		if(!found && !onlyOverlapping && uniReg.getLen() >=minlen){
			auto regionCopy = uniReg;
			regionCopy.uid_ = njh::pasteAsStr(regionCopy.start_, "-", regionCopy.end_);
			trimRegions.emplace_back(regionCopy);
		}else if(found){
			overlappingRegions.emplace_back(uniReg);
		}
	}
	//for debugging;
	if(false){
		std::cout << "overlapping regions" << std::endl;
		for(const auto & reg : overlappingRegions){
			std::cout << reg.genBedRecordCore().toDelimStr() << std::endl;
		}
		std::cout << "non overlapping regions" << std::endl;
		for(const auto & reg : trimRegions){
			std::cout << reg.genBedRecordCore().toDelimStr() << std::endl;
		}
	}
	std::unordered_map<std::string, uint32_t> contigsCoveredCounts;
	uint32_t pos = 0;
	while(pos < len(getSeqBase(refSeq))){
		std::vector<GenomicRegion> currentRegions;
		for(const auto & reg : overlappingRegions){
			if(pos >= reg.start_ && pos < reg.end_){
				currentRegions.emplace_back(reg);
			}
			//can break once region's start is greater than pos as they are sorted and therefore no more regions can overlap
			if(reg.start_ > pos){
				break;
			}
		}
		if (currentRegions.size() > 1) {
			uint32_t start = std::max_element(currentRegions.begin(), currentRegions.end(), [](const GenomicRegion & reg1,const GenomicRegion & reg2 ){
				return reg1.start_ < reg2.start_;
			})->start_;
			uint32_t end = std::min_element(currentRegions.begin(), currentRegions.end(), [](const GenomicRegion & reg1,const GenomicRegion & reg2 ){
				return reg1.end_ < reg2.end_;
			})->end_;
			if(end < start){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error, end, " << end << ", shouldn't be able to be less than start, " << start << "" << "\n";
				throw std::runtime_error{ss.str()};
			}
			uint32_t length = end - start;
			if(length >= minlen){
				trimRegions.emplace_back(njh::pasteAsStr(start, "-", end), getSeqBase(refSeq).name_, start, end , false);
				if(njh::in(trimRegions.back().uid_, contigsCoveredCounts)){
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error, duplicate region for start: " << start  << ", end: " << end << "\n";
					throw std::runtime_error{ss.str()};
				}else{
					contigsCoveredCounts[trimRegions.back().uid_] = currentRegions.size();
				}
				if(start >= overlappingRegions.back().start_){
					//at the last region, end
					pos = len(getSeqBase(refSeq));
				} else {
					for(const auto & nextReg : overlappingRegions){
						if(nextReg.start_ > start){
							pos = nextReg.start_;
							break;
						}
					}
				}
			}else{
				++pos;
			}
		} else {
			++pos;
		}
	}
	sortGRegionsByStart(trimRegions);

	//slim down regions being found completely each other, give preference to which ever region covers more
	std::vector<size_t> toBeRemoved;
	for(const auto pos : iter::range(trimRegions.size())){
		for(const auto pos2 : iter::range(pos + 1,trimRegions.size())){
			if(!njh::in(pos, toBeRemoved) && !njh::in(pos2, toBeRemoved)){
				const auto olen = trimRegions[pos].getOverlapLen(trimRegions[pos2]);
				//if one of the length equals the overlap length, then one region is completely within the other, pick the one that covers more bits;
				if(olen == trimRegions[pos].getLen() || olen == trimRegions[pos2].getLen()){
					if (std::max(trimRegions[pos].getLen(), trimRegions[pos2].getLen())
							/ static_cast<double>(std::min(trimRegions[pos].getLen(),
									trimRegions[pos2].getLen())) > 4.0) {
						//avoid removing regions that are much larger
						continue;
					}
					//optimize # of bits
					if(contigsCoveredCounts[trimRegions[pos].uid_] > contigsCoveredCounts[trimRegions[pos2].uid_]){
						toBeRemoved.emplace_back(pos2);
					} else if(contigsCoveredCounts[trimRegions[pos].uid_] == contigsCoveredCounts[trimRegions[pos2].uid_]){
						//optimize length if same # of bits
						if(trimRegions[pos].getLen() > trimRegions[pos2].getLen()){
							toBeRemoved.emplace_back(pos2);
						} else {
							toBeRemoved.emplace_back(pos);
						}
					} else {
						toBeRemoved.emplace_back(pos);
					}
				}
			}
		}
	}
	//sort and remove latter positions so positions don't become invalidated
	njh::sort(toBeRemoved);
	for(const auto & rem : iter::reversed(toBeRemoved)){
//		//get regions that almost encompass the whole or more of the reference sequence to
//		//avoid removing larger segments that only have a small
//		if(trimRegions[rem].getLen()/static_cast<double>(len(getSeqBase(refSeq))) < .95 ){
//			trimRegions.erase(trimRegions.begin() + rem);
//		}
		trimRegions.erase(trimRegions.begin() + rem);
	}
	sortGRegionsByStart(trimRegions);
	return trimRegions;
}


///
VecStr getPossibleDNASubstringsForProtein(const std::string& seq,
                                          const std::string& protein,
                                          const std::string& seqType = "DNA");
VecStr findPossibleDNA(const std::string& seq, const std::string& protein,
                       const std::string& seqType = "DNA",
                       bool checkComplement = true);

VecStr getAllCycloProteinFragments(const std::string& protein);

std::multimap<int, std::string> getProteinFragmentSpectrum(
    const std::string& protein);

std::vector<int> getRealPossibleWeights(const std::vector<int>& spectrum);



VecStr organizeLexicallyKmers(const std::string& input, size_t colNum);
VecStr organizeLexicallyKmers(const VecStr& input, int colNum);
uint64_t smallestSizePossible(uint64_t weight);


int64_t getPossibleNumberOfProteins(
    int64_t weight, std::unordered_map<int64_t, int64_t>& cache);
probabilityProfile randomMotifSearch(const VecStr& dnaStrings, int kLength,
                                     int numberOfRuns, bool gibs, int gibsRuns,
																		 njh::randomGenerator& gen);
probabilityProfile randomlyFindBestProfile(const VecStr& dnaStrings,
                                           const std::vector<VecStr>& allKmers,
                                           int numberOfKmers,
																					 njh::randomGenerator& gen);
probabilityProfile randomlyFindBestProfileGibs(
    const VecStr& dnaStrings, const std::vector<VecStr>& allKmers,
    int numberOfKmers, int runTimes, njh::randomGenerator& gen);

template <typename T>
std::vector<std::vector<T>> getAllVecCycloFragments(const std::vector<T>& vec) {
  std::vector<std::vector<T>> ans;
  ans.emplace_back(std::vector<int>(0));
  for (auto i : iter::range<uint32_t>(1, vec.size())) {
    for (auto j : iter::range(vec.size())) {
      std::vector<T> currentFragment;
      if (j + i > vec.size()) {
        currentFragment =
            concatVecs(getSubVector(vec, j, vec.size() - j),
                            getSubVector(vec, 0, i - (vec.size() - j)));
      } else {
        currentFragment = getSubVector(vec, j, i);
      }
      ans.push_back(currentFragment);
    }
  }
  ans.push_back(vec);
  std::sort(ans.begin(), ans.end());
  return ans;
}
template <typename T>
std::vector<std::vector<T>> getAllVecLinearFragments(
    const std::vector<T>& vec) {
  std::vector<std::vector<T>> ans;
  ans.emplace_back(std::vector<int>(0));
  for (auto i : iter::range<uint32_t>(1, vec.size())) {
    for (auto j : iter::range(vec.size())) {
      std::vector<T> currentFragment;
      if (j + i > vec.size()) {

      } else {
        currentFragment = getSubVector(vec, j, i);
        ans.push_back(currentFragment);
      }
    }
  }
  ans.push_back(vec);
  std::sort(ans.begin(), ans.end());
  return ans;
}
VecStr getAllLinearProteinFragments(const std::string& protein);
std::vector<std::vector<int>> getPossibleProteinsForSpectrum(
    const std::string& spectrum, bool useAllWeights = false);
bool trimCycle(std::vector<std::vector<int>>& nextCycle,
               std::vector<std::vector<int>>& matchesSpectrum,
               const std::map<int, int>& spectrumToWeights,
               const std::vector<int> specVec);
std::vector<std::vector<int>> growNextCycle(
    const std::vector<std::vector<int>>& previousCycle,
    const std::vector<int>& possibleWeights);
int scoreSpectrumAgreement(const std::map<int, int>& spectrum,
                           const std::map<int, int>& currentSpectrum);
std::multimap<int, std::vector<int>, std::greater<int>> growNextCycleScore(
    const std::multimap<int, std::vector<int>, std::greater<int>>&
        previousCycle,
    const std::vector<int>& possibleWeights,
    const std::map<int, int>& spectrumCounts, int parentMass, bool linear);
std::multimap<int, std::vector<int>, std::greater<int>> trimCycleScore(
    std::multimap<int, std::vector<int>, std::greater<int>>& nextCycle,
    std::multimap<int, std::vector<int>, std::greater<int>>& matchesSpectrum,
    int parentMass, int leaderBoardNumber, int& currentLeader);
std::multimap<int, std::vector<int>, std::greater<int>>
    getPossibleProteinsForSpectrum(const std::string& spectrum,
                                   int leaderBoardNumber, bool verbose,
                                   bool useAllWeights, bool convolution,
                                   int convolutionCutOff, bool linear);
std::map<int, int> getConvolutionWeights(std::vector<int> experimentalSpectrum,
                                         int multiplicityCutOff,
                                         int lowerBound = 57,
                                         int upperBound = 200);
std::vector<int> convolutionWeights(std::vector<int> experimentalSpectrum,
                                    int multiplicityCutOff, int lowerBound = 57,
                                    int upperBound = 200);
std::vector<int> topConvolutionWeights(std::vector<int> experimentalSpectrum,
                                       int mItems, int lowerBound = 57,
                                       int upperBound = 200);
int64_t getMinCoins(int64_t change, const std::vector<int64_t>& coins,
                    std::unordered_map<int64_t, int64_t>& cache);









std::vector<uint32_t> getWindowQuals(const std::vector<uint32_t>& qual,
                                     uint32_t qWindowSize, uint32_t pos);
std::vector<double> likelihoodForBaseQ(
    const std::vector<uint32_t>& qual,
    std::unordered_map<double, double>& likelihoods);
std::vector<double> likelihoodForMeanQ(
    const std::vector<uint32_t>& qual, uint32_t qWindowSize,
    std::unordered_map<double, double>& likelihoods);
std::vector<double> likelihoodForMedianQ(
    const std::vector<uint32_t>& qual, uint32_t qWindowSize,
    std::unordered_map<double, double>& likelihoods);
std::vector<double> likelihoodForMinQ(
    const std::vector<uint32_t>& qual, uint32_t qWindowSize,
    std::unordered_map<double, double>& likelihoods);

double getChangeInHydro(const char& firstAA, const char& secondAA);
template<typename T>
std::vector<double> getHydroChanges(const std::string& originalCodon,
                                    const VecStr& mutantCodons,
                                    const T& code) {
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
}



//seqInfo extractSeq(TwoBit::TwoBitFile & twobitReader, const GenomicRegion & region);


void simpleCollapse(std::vector<cluster> & consensusReads, aligner & alignerObj,
		const comparison & passableErrors);
void simpleCollapseQueryCov(std::vector<cluster> & consensusReads, aligner & alignerObj,
		const comparison & passableErrors);



}  // namespace njh


