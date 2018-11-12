/*
 * MinorVariantsMixture.cpp
 *
 *  Created on: Jan 9, 2016
 *      Author: nick
 */
//
// njhseq - A library for analyzing sequence data
// Copyright (C) 2012-2016 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
// Jeffrey Bailey <Jeffrey.Bailey@umassmed.edu>
//
// This file is part of njhseq.
//
// njhseq is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// njhseq is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with njhseq.  If not, see <http://www.gnu.org/licenses/>.
//

#include "MinorVariantsMixture.hpp"
#include "elucidator/simulation/randomStrGen.hpp"


namespace njhseq {




size_t getUniquePos(njh::randomGenerator & gen,
		std::set<size_t> & allPositions, size_t start, size_t end) {
	if (end <= start) {
		std::stringstream ss;
		ss << "Error in " << __PRETTY_FUNCTION__
				<< " end should be greater than start, start: " << start << ", end: "
				<< end << std::endl;
		throw std::runtime_error { ss.str() };
	}
	if (allPositions.size() >= end - start) {
		std::stringstream ss;
		ss << "Error in " << __PRETTY_FUNCTION__ << " all positions already chosen"
				<< std::endl;
		throw std::runtime_error { ss.str() };
	}
	size_t pos = gen.unifRand<size_t>(start, end);
	while (njh::in(pos, allPositions)) {
		pos = gen.unifRand<size_t>(start, end);
	}
	allPositions.insert(pos);
	return pos;
}

uint32_t MinorVariantsMixture::numOfVars() const {
	return variants_.size();
}

void MinorVariantsMixture::addVariants(const std::vector<std::set<size_t>> & mutPositions){
	std::unordered_map<uint32_t, char> alreadyMutated;
	uint32_t varNum = 0;
	std::string varNameStub = "Minor";
	for (const auto & positions : mutPositions) {
		std::string minorVar = major_->seq_;
		std::unordered_map<size_t, char> varMutations;
		for (const auto pos : positions) {
			if(pos > len(*major_)){
				std::stringstream ss;
				ss << njh::bashCT::red << "Error in " << __PRETTY_FUNCTION__
						<< ": trying to mutate position greater than length of major"
						<< ", len of major: " << len(*major_) << ", pos: " << pos
						<< njh::bashCT::reset << std::endl;
				throw std::runtime_error{ss.str()};
			}
			char mutChar = gen_.unifRandSelection(dnaAlphabet_);
			if (njh::in(pos, alreadyMutated)) {
				mutChar = alreadyMutated[pos];
			} else {
				while (mutChar == major_->seq_[pos]) {
					mutChar = gen_.unifRandSelection(dnaAlphabet_);
				}
				alreadyMutated[pos] = mutChar;
			}
			varMutations.emplace(pos, mutChar);
			minorVar[pos] = mutChar;
		}
		std::string varname = varNameStub + "." + leftPadNumStr<uint32_t>(varNum, mutPositions.size() - 1);
		variants_.emplace_back(std::make_shared<seqInfo>(varname, minorVar));
		variantsMutations_[varname] = varMutations;
		++varNum;
	}
}



void MinorVariantsMixture::createDisparatePairsMixture(uint32_t numberOfMinorStrains){
	std::vector<std::set<size_t>> mutPositions;
	std::vector<size_t> allPossilblePositions(len(*major_));
	njh::iota<size_t>(allPossilblePositions,0);
	std::vector<size_t> positionsPool = gen_.unifRandSelectionVec(allPossilblePositions, 10, false);
	std::set<size_t> allPositions(positionsPool.begin(), positionsPool.end());
	/**@todo default is minor variants are 15 mutations away from major, consider making it so it can be set*/
	for(uint32_t i = 0; i < numberOfMinorStrains; ++i){
		std::vector<size_t> fromPool = gen_.unifRandSelectionVec(positionsPool, 8, false);
		std::vector<size_t> currentPos;
		//get uniuqe postions for the rest
		for(uint32_t j = 0; j < 7; ++j){
			size_t pos = getUniquePos(gen_, allPositions, 0, len(*major_));
			currentPos.emplace_back(pos);
		}
		addOtherVec(currentPos, fromPool);
		mutPositions.emplace_back(currentPos.begin(), currentPos.end());
		//gen_ another minor that's one off of the previous
		{
			size_t pos = getUniquePos(gen_, allPositions, 0, len(*major_));
			currentPos.emplace_back(pos);
			mutPositions.emplace_back(currentPos.begin(), currentPos.end());
		}
	}
	addVariants(mutPositions);
}

void MinorVariantsMixture::createGradientDisparatePairsMixture(uint32_t numberOfMinorStrains){
	std::vector<std::set<size_t>> mutPositions;
	std::vector<size_t> allPossilblePositions(len(*major_));
	njh::iota<size_t>(allPossilblePositions, 0);
	std::vector<size_t> positionsPool = gen_.unifRandSelectionVec(
			allPossilblePositions, 10, false);
	std::set<size_t> allPositions(positionsPool.begin(), positionsPool.end());
	uint32_t mutStep = 1;
	auto varyingMuts = initialNumMutations_;
	if (numberOfMinorStrains > varyingMuts.size()) {
		while (varyingMuts.size() < numberOfMinorStrains) {
			varyingMuts.push_back(varyingMuts.back() + mutStep);
		}
	}
	/**@todo make sure these while loops don't become infinite loops if all positions chosen*/
	for (uint32_t i = 0; i < numberOfMinorStrains; ++i) {
		std::vector<size_t> fromPool = gen_.unifRandSelectionVec(positionsPool,
				8, false);
		std::vector<size_t> currentPos;
		//get uniuqe postions for the rest
		for (uint32_t j = 0; j < 7; ++j) {
			size_t pos = getUniquePos(gen_, allPositions,0, len(*major_));
			currentPos.emplace_back(pos);
		}
		addOtherVec(currentPos, fromPool);
		mutPositions.emplace_back(currentPos.begin(), currentPos.end());
		//gen_ another minor that's one off of the previous
		for(uint32_t j = 0; j <varyingMuts[i]; ++j){
			size_t pos = getUniquePos(gen_, allPositions,0, len(*major_));

			currentPos.emplace_back(pos);
		}
		mutPositions.emplace_back(currentPos.begin(), currentPos.end());
	}
	addVariants(mutPositions);
}



void MinorVariantsMixture::createVaryingMinorVariantsMixture(uint32_t numberOfMinorVariants, double oneVariantPerIdDiff){
	std::vector<std::set<size_t>> mutPositions;
	uint32_t maxDiffs = ::ceil(major_->seq_.size() * oneVariantPerIdDiff);
	std::vector<uint32_t> mutations{maxDiffs};
	auto initialMuts = initialNumMutations_;
	if(njh::in(maxDiffs, initialMuts)){
		removeElement(initialMuts, maxDiffs);
	}
	uint32_t mutStep = 1;
	if(numberOfMinorVariants > 1){
		if(numberOfMinorVariants - 1 < initialMuts.size()){
			initialMuts = getSubVector(initialMuts, 0, numberOfMinorVariants - 1);
		}else{
			while(initialMuts.size() < numberOfMinorVariants - 1){
				initialMuts.push_back(initialMuts.back() + mutStep);
			}
		}
		addOtherVec(mutations, initialMuts);
	}
	njh::sort(mutations);
	std::set<size_t> allPositions;
	for(auto mutNum : mutations){
		std::set<size_t> currentPos;
		for(uint32_t i = 0; i < mutNum; ++i){
			auto pos = getUniquePos(gen_, allPositions, 0, len(*major_));
			currentPos.insert(pos);
		}
		mutPositions.emplace_back(currentPos.begin(), currentPos.end());
	}
	addVariants(mutPositions);
}

void MinorVariantsMixture::setNewMajor(size_t stringLength){
	std::string major = simulation::evenRandStr(stringLength, dnaAlphabet_ , gen_);
	major_ = std::make_shared<seqInfo>("Major", major);
}

MinorVariantsMixture::MinorVariantsMixture(size_t stringLength){
  //sim the major strain
  setNewMajor(stringLength);
}

void MinorVariantsMixture::addAnEqualAbundanceMixture(){
	double abund = 100.0 / (variants_.size() + 1);
	std::string mixName = "mix." + estd::to_string(mixturesAbundances_.size());
	mixturesAbundances_[mixName][major_->name_] = abund;
	for(const auto & var : variants_){
		mixturesAbundances_[mixName][var->name_] = abund;
	}
}

void MinorVariantsMixture::prependPaddingSeq(const std::string & seq){
	major_->prepend(seq);
	for(auto & var : variants_){
		var->prepend(seq);
	}
}

void MinorVariantsMixture::appendPaddingSeq(const std::string & seq){
	major_->append(seq);
	for(auto & var : variants_){
		var->append(seq);
	}
}

MinorVariantsMixture::MinorVariantsMixture(const seqInfo & major) :
		major_(std::make_shared<seqInfo>(major)) {
}

void MinorVariantsMixture::addVariants(const std::string & mutationsString, size_t mutationOffSet){
  //process mutations string
	if(mutationsString.size() == 0){
		std::stringstream ss;
		ss << njh::bashCT::red << "Error in: " << __PRETTY_FUNCTION__ << ": "
				<< "Error with processing mutations: "
				<< "mutation string is empty"
				<< njh::bashCT::reset << std::endl;
		throw std::runtime_error { ss.str() };
	}
	VecStr variants;
	for(const auto & var : variants_){
		variants.emplace_back(var->seq_);
	}
	//generate mutations
	auto mutToks = tokenizeString(mutationsString, ",");
	for (const auto & tok : mutToks) {
		std::string currentVariantName = major_->name_ + "-var" + estd::to_string(variants.size());
		std::unordered_map<size_t, char> currentVariantMutations;
		if (tok.find(":") != std::string::npos) {
			std::string variant = major_->seq_;
			std::vector<size_t> alreadyMutated;
			auto specifcMutToks = tokenizeString(tok, ";");
			for(const auto & smTok : specifcMutToks){
				auto subMutToks = tokenizeString(smTok, ":");
				if (subMutToks.size() != 2) {
					std::stringstream ss;
					ss << njh::bashCT::red  << "Error in: " << __PRETTY_FUNCTION__ << ": "
							<< "Error with processing mutation: " << tok
							<< ", should be a number and character separated by a :"
							<< njh::bashCT::reset << std::endl;
					throw std::runtime_error { ss.str() };
				}
				if (isIntStr(subMutToks[0]) && subMutToks[1].size() == 1
						&& isalpha(subMutToks[1][0])) {
					auto position = njh::lexical_cast<size_t>(subMutToks[0]);
					char base = subMutToks[1][0];
					if(position < len(*major_)){
						if(major_->seq_[position] != base){
							if(!njh::in(position, alreadyMutated )){
								alreadyMutated.push_back(position);
								variant[position] = base;
								currentVariantMutations[position] = base;
								alreadyMutatedPosBase_[position] = base;
							}else{
								std::stringstream ss;
								ss << njh::bashCT::red << "Error in: " << __PRETTY_FUNCTION__ << ": "
										 << "Error with processing mutation: " << tok
										<< ", requested a mutation at position " << position
										<< " but that position was already mutated"
										<< njh::bashCT::reset << std::endl;
								throw std::runtime_error { ss.str() };
							}
						}else{
							std::stringstream ss;
							ss << njh::bashCT::red  << "Error in: " << __PRETTY_FUNCTION__ << ": "
									<< "Error with processing mutation: " << tok
									<< ", requested a mutation at position " << position << " to "
									<< base << " but the reference sequence is already this base"
									<< njh::bashCT::reset << std::endl;
							throw std::runtime_error { ss.str() };
						}
					}else{
						std::stringstream ss;
						ss << njh::bashCT::red  << "Error in: " << __PRETTY_FUNCTION__ << ": "
								<< "Error with processing mutation: " << tok
								<< "position was out of range, seq size: " << len(*major_)
								<< ", position requested, " << position
								<< njh::bashCT::reset << std::endl;
						throw std::out_of_range{ ss.str() };
					}
				} else {
					std::stringstream ss;
					ss << njh::bashCT::red  << "Error in: " << __PRETTY_FUNCTION__ << ": "
							<< "Error with processing mutation: " << tok
							<< ", should be a number and character separated by a :"
							<< njh::bashCT::reset << std::endl;
					throw std::runtime_error { ss.str() };
				}
			}
			if (!njh::in(variant, variants) && variant != major_->seq_) {
				variantsMutations_[currentVariantName] = currentVariantMutations;
				variants.emplace_back(variant);
				variants_.emplace_back(std::make_shared<seqInfo>(currentVariantName, variant));
			} else {
				std::stringstream ss;
				ss << njh::bashCT::red << "Error with creating variant: " << tok
						<< "variant already created with: " << tok
						<< njh::bashCT::reset << std::endl;
				throw std::runtime_error { ss.str() };
			}
		} else {
			if(mutationOffSet * 2 < len(*major_) ){
				currentVariantMutations.clear();
				size_t start = mutationOffSet;
				size_t stop = len(*major_) - mutationOffSet;
				std::vector<size_t> allPositions(stop - start);
				njh::iota<size_t>(allPositions, start);
				std::string variant = major_->seq_;
				while(njh::in(variant, variants) || variant == major_->seq_){
					variant = major_->seq_;
					std::vector<size_t> positions;
					bool alreadyContainsPositions = true;
					while(alreadyContainsPositions){
						positions = gen_.unifRandSelectionVec(allPositions, estd::stou(tok),false);
						bool pass = true;
						for(auto & pos : positions){
							if(njh::in(pos, alreadyMutatedPosBase_)){
								pass = false;
								break;
							}
						}
						if(pass){
							alreadyContainsPositions = false;
						}
					}
					for(const auto pos : positions){
						char refBase = major_->seq_[pos];
						char mutBase = refBase;
						while(mutBase == refBase){
							mutBase = gen_.unifRandSelection(dnaAlphabet_);
						}
						variant[pos] = mutBase;
						currentVariantMutations[pos] = mutBase;
						alreadyMutatedPosBase_[pos] = mutBase;
					}
				}
				variants.emplace_back(variant);
				variantsMutations_[currentVariantName] = currentVariantMutations;
				variants_.emplace_back(std::make_shared<seqInfo>(currentVariantName, variant));
			}else{
				std::stringstream ss;
				ss << njh::bashCT::red << "Error with creating variant: " << tok
						<< "mutation offset for mutation distance from edge is too large for current seq size,"
						<< "distance should be less than half the distance of reference seq, offSet "
						<< mutationOffSet << ", seqSize: " << len(*major_)
						<< njh::bashCT::reset << std::endl;
				throw std::runtime_error { ss.str() };
			}
		}
	}
}

MinorVariantsMixture::MinorVariantsMixture(const seqInfo & major,
		const std::string & mutationsString, size_t mutationOffSet) :
		major_(std::make_shared<seqInfo>(major)) {
	addVariants(mutationsString, mutationOffSet);
}

void MinorVariantsMixture::createAbundancesForMixtures(
		const std::string & abundances) {
	if (abundances.size() == 0) {
		std::stringstream ss;
		ss << njh::bashCT::red << "Error with processing mutations: "
				<< "mutation string is empty" << njh::bashCT::reset << std::endl;
		throw std::runtime_error { ss.str() };
	}
	auto abundToks = tokenizeString(abundances, ",");
	uint32_t mixtureNumber = mixturesAbundances_.size();
	for (const auto & abTok : abundToks) {
		std::string mixName = "mix." + estd::to_string(mixtureNumber);;
		if (abTok.find(":") != std::string::npos) {
			auto spAbToks = tokenizeString(abTok, ":");
			if (spAbToks.size() == numOfVars()) {
				double totalAbundance = 0;
				for (const auto & spAbTokPos : iter::range(spAbToks.size())) {
					const auto & spAbTok = spAbToks[spAbTokPos];
					if (isDoubleStr(spAbTok)) {
						auto abund = njh::lexical_cast<double>(spAbTok);
						totalAbundance += abund;
						mixturesAbundances_[mixName][variants_[spAbTokPos]->name_] = abund;
					} else {
						std::stringstream ss;
						ss << njh::bashCT::red << "Error with processing abundance: "
								<< abTok << ", found non numeric character"
								<< njh::bashCT::reset << std::endl;
						throw std::runtime_error { ss.str() };
					}
				}
				if (totalAbundance >= 100) {
					std::stringstream ss;
					ss << njh::bashCT::red << "Error with processing abundance: " << abTok
							<< ", requesting too much abundance, should be less than 100 total, current total: "
							<< totalAbundance << njh::bashCT::reset << std::endl;
					throw std::runtime_error { ss.str() };
				}
				mixturesAbundances_[mixName][major_->name_] = 100.0 - totalAbundance;
			} else {
				std::stringstream ss;
				ss << njh::bashCT::red << "Error with processing abundance: " << abTok
						<< ", when setting specific abundances for each variant, need to supply the number of variants created, number of variatns: "
						<< numOfVars() << ", number of abundances, " << spAbToks.size()
						<< njh::bashCT::reset << std::endl;
				throw std::runtime_error { ss.str() };
			}
		} else {
			if (isDoubleStr(abTok)) {
				auto abund = njh::lexical_cast<double>(abTok);
				if (abund * numOfVars() < 100) {
					double totalAbundance = 0;
					for (const auto & var : variants_) {
						mixturesAbundances_[mixName][var->name_] = abund;
						totalAbundance+= abund;
					}
					mixturesAbundances_[mixName][major_->name_] = 100.0 - totalAbundance;
				} else {
					std::stringstream ss;
					ss << njh::bashCT::red << "Error with processing abundance: " << abTok
							<< ", requesting too much abundance, minimum amount possible is "
							<< 99.0 / numOfVars() << njh::bashCT::reset << std::endl;
					throw std::runtime_error { ss.str() };
				}
			} else {
				std::stringstream ss;
				ss << njh::bashCT::red << "Error with processing abundance: " << abTok
						<< ", found non numeric character" << njh::bashCT::reset
						<< std::endl;
				throw std::runtime_error { ss.str() };
			}
		}
		++mixtureNumber;
	}
}

void MinorVariantsMixture::outputAbundanceFile(
		const std::string & abundanceFilename, bool dualReplicates) const {
	std::ofstream abundanceFile;
	openTextFile(abundanceFile, OutOptions(abundanceFilename, ".tab.txt"));
	std::map<uint32_t, std::string> mixSorter;
	for(const auto & mix : mixturesAbundances_){
		mixSorter[estd::stou(mix.first.substr(mix.first.find(".") + 1))] = mix.first;
	}
	auto mixKeys = getVectorOfMapValues(mixSorter);
	//log major abndances;
	abundanceFile << major_->name_;
	for(const auto & mixKey : mixKeys ){
		const auto & mix = mixturesAbundances_.at(mixKey);
		abundanceFile << "\t" << mix.at(major_->name_);
		if(dualReplicates){
			abundanceFile << "\t" << mix.at(major_->name_);
		}
	}
	abundanceFile << std::endl;
	auto varNames = readVec::getNames(variants_);
	njh::sort(varNames);
	for(const auto & var : varNames){
		abundanceFile << var;
		for(const auto & mixKey : mixKeys ){
			const auto & mix = mixturesAbundances_.at(mixKey);
			abundanceFile << "\t" << mix.at(var);
			if(dualReplicates){
				abundanceFile << "\t" << mix.at(var);
			}
		}
		abundanceFile << std::endl;
	}
}

void MinorVariantsMixture::outputVariantMutationInfofile(
		const std::string & filename) const {
  //create README for variants information
	std::ofstream variantsInfoFile;
	openTextFile(variantsInfoFile, OutOptions(filename, ".tab.txt"));
	variantsInfoFile << "variantName\trefPos\trefBase\tvariantBase" << std::endl;
	auto varKeys = getVectorOfMapKeys(variantsMutations_);
	njh::sort(varKeys);
	for(const auto & varInfoKey : varKeys){
		auto keys = getVectorOfMapKeys(variantsMutations_.at(varInfoKey));
		njh::sort(keys);
		for(const auto & key : keys){
			variantsInfoFile << varInfoKey
					<< "\t" << key
					<< "\t" << major_->seq_[key]
					<< "\t" << variantsMutations_.at(varInfoKey).at(key) << std::endl;
		}
	}
}

void MinorVariantsMixture::writeFasta(const std::string & dirname)const{
	writeFasta(dirname,  major_->name_ + "_vars.fasta");
}

void MinorVariantsMixture::writeFasta(const std::string & dirname, const std::string & filename)const{
	SeqIOOptions opts = SeqIOOptions::genFastaOut(njh::files::join(dirname, filename));
	SeqIO writer(opts);
	writer.openOut();
	writer.write(major_);
	writer.write(variants_);
}

}  // namespace njhseq
