#include "refVariants.hpp"
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
namespace njhseq {
variant::variant(const seqInfo & seqBase) :
		seqBase_(seqBase) {
}

void variant::outputInfo(std::ofstream & out, const seqInfo & ref) const {
	//mismatches
	for (const auto & m : mismatches_) {
		out << ref.name_ << "\t" << m.first << "\t" << seqBase_.name_ << "\t"
				<< m.second.refBase << "\t" << m.second.seqBase << "\t"
				<< seqInfo::getFastqString( { m.second.seqQual }, SangerQualOffset)
				<< std::endl;
	}
	//insertions
	for (const auto & i : insertions_) {
		out << ref.name_ << "\t" << i.first << "\t" << seqBase_.name_ << "\t" << "-"
				<< "\t" << i.second.gapedSequence_ << "\t"
				<< seqInfo::getFastqString(i.second.qualities_, SangerQualOffset)
				<< std::endl;
	}
	//deletions
	for (const auto & d : deletions_) {
		out << ref.name_ << "\t" << d.first << "\t" << seqBase_.name_ << "\t"
				<< d.second.gapedSequence_ << "\t" << "_" << "\t"
				<< seqInfo::getFastqString(d.second.qualities_, SangerQualOffset)
				<< std::endl;
	}
}


refVariants::refVariants(const seqInfo & seqBase): seqBase_(seqBase){}

void refVariants::addVariant(const seqInfo & var, aligner & alignerObj,
		bool weighHomopolymer) {
	variant varObj(var);
	alignerObj.alignCache(seqBase_, var, false);
	alignerObj.profilePrimerAlignment(seqBase_, var);
	for (const auto & m : alignerObj.comp_.distances_.mismatches_) {
		varObj.mismatches_.emplace(m.second.refBasePos, m.second);
	}
	for (const auto & g : alignerObj.comp_.distances_.alignmentGaps_) {
		if (g.second.ref_) {
			//insertion
			varObj.insertions_.emplace(g.second.refPos_,
					g.second);
		} else {
			//deletion
			varObj.deletions_.emplace(g.second.refPos_,
					g.second);
		}
	}
	variants_.emplace_back(varObj);
}
std::vector<uint32_t> refVariants::getVariantSnpLoci()const{
	std::set<uint32_t> loci;
	for(const auto & v : variants_){
		for(const auto & m : v.mismatches_){
			loci.emplace(m.second.refBasePos);
		}
	}
	return std::vector<uint32_t> {loci.begin(), loci.end()};
}

std::map<uint32_t, std::vector<char>> refVariants::getVariantSnpLociMap()const{
	std::map<uint32_t, std::set<char>> lociChars;
	for(const auto & v : variants_){
		for(const auto & m : v.mismatches_){
			lociChars[m.second.refBasePos].emplace(m.second.seqBase);
		}
	}
	std::map<uint32_t, std::vector<char>> ret;
	for(const auto & l : lociChars){
		ret[l.first] =std::vector<char> {l.second.begin(), l.second.end()};
	}
	return ret;
}


std::vector<uint32_t> refVariants::getUniqueToRefPositions() const {
	std::set<uint32_t> positions;
	for(const auto & v : variants_){
		for(const auto & m : v.mismatches_){
			positions.emplace(m.first);
		}
		for(const auto & d : v.deletions_){
			for(uint32_t pos = 0; pos < d.second.gapedSequence_.size(); ++pos){
				positions.emplace(d.first + pos);
			}
		}
	}
	return std::vector<uint32_t> {positions.begin(), positions.end()};
}

std::vector<uint32_t> refVariants::getVariantSnpLoci(VecStr names, uint32_t expand )const{
	std::set<uint32_t> loci;
	for(const auto & v : variants_){
		if(njh::in(v.seqBase_.name_, names)){
			for(const auto & m : v.mismatches_){
				loci.emplace(m.second.refBasePos);
				if(expand > 0){
					for(const auto & e : iter::range<uint32_t>(1,expand + 1)){
						if(e <=m.second.refBasePos){
							loci.emplace(m.second.refBasePos - e);
						}
						if(e + m.second.refBasePos < seqBase_.seq_.size()){
							loci.emplace(m.second.refBasePos + e);
						}
					}
				}
			}
		}
	}
	return std::vector<uint32_t> {loci.begin(), loci.end()};
}

std::map<uint32_t, std::vector<char>> refVariants::getVariantSnpLociMap(VecStr names, uint32_t expand )const{
	std::map<uint32_t, std::set<char>> lociChars;
	for(const auto & v : variants_){
		if(njh::in(v.seqBase_.name_, names)){
			for(const auto & m : v.mismatches_){
				lociChars[m.second.refBasePos].emplace(m.second.seqBase);
				if(expand > 0){
					for(const auto & e : iter::range<uint32_t>(1,expand + 1)){
						if(e <=m.second.refBasePos){
							lociChars[m.second.refBasePos - e].emplace(m.second.seqBase);
						}
						if(e + m.second.refBasePos < seqBase_.seq_.size()){
							lociChars[m.second.refBasePos + e].emplace(m.second.seqBase);
						}
					}
				}
			}
		}
	}
	std::map<uint32_t, std::vector<char>> ret;
	for(const auto & l : lociChars){
		ret[l.first] =std::vector<char> {l.second.begin(), l.second.end()};
	}
	return ret;
}

void refVariants::outPut(std::ofstream & out)const {
	for(const auto & v : variants_){
		v.outputInfo(out, seqBase_);
	}
}

}  // namespace njhseq
