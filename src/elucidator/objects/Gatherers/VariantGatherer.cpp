/*
 * SnpGatherer.cpp
 *
 *  Created on: Jun 20, 2016
 *      Author: nick
 */
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


#include "VariantGatherer.hpp"
#include "elucidator/filesystem/wrappers/BGZFCPP.hpp"
#include "elucidator/objects/SlimCounter/SlimCounterPos.hpp"



namespace njhseq {

VariantGatherer::VariantGatherer(const std::vector<BamTools::RefData> & refInfos,
		double fracCutOff, uint32_t depthCutOff, double strandBiasCutOff) :
		BaseGatherer(refInfos), fracCutOff_(fracCutOff), depthCutOff_(depthCutOff), strandBiasCutOff_(
				strandBiasCutOff) {
}



VariantGatherer::VarPosRes VariantGatherer::determineIfPosVariant(const SlimCounterPos & pos,
		const char refBase,
		const GathVarPars pars){
	std::vector<char> bases { 'A', 'C', 'G', 'T' };
	//check if there are more than 1 base above the fraccutoff which
	uint32_t numberAboveCutOff = 0;
	double bestBaseFrac = std::numeric_limits<double>::lowest();
	char bestBase = '-';
	for (auto base : bases) {
		//check strand bias
		if (pos.getHqFwBaseFrac(base) > strandBiasCutOff_
				&& pos.getHqFwBaseFrac(base) < (1 - strandBiasCutOff_)) {
			if (pos.getHqBaseFrac(base) > fracCutOff_) {
				++numberAboveCutOff;
				if(bestBaseFrac < pos.getHqBaseFrac(base)){
					bestBaseFrac = pos.getHqBaseFrac(base);
					bestBase = base;
				}
			}
		}
	}
	if (pars.addIndels_) {
		//check strand bias
		if (pos.getInsertionsFrac() > strandBiasCutOff_
				&& pos.getInsertionsFrac() < (1 - strandBiasCutOff_)) {
			if (pos.getInsertionsFrac() > fracCutOff_) {
				++numberAboveCutOff;
				if(bestBaseFrac < pos.getInsertionsFrac()){
					bestBaseFrac = pos.getInsertionsFrac();
					bestBase = '-';
				}
			}
		}
		if (pos.getDeletionsFrac() > strandBiasCutOff_
				&& pos.getDeletionsFrac() < (1 - strandBiasCutOff_)) {
			if (pos.getDeletionsFrac() > fracCutOff_) {
				++numberAboveCutOff;
				if(bestBaseFrac < pos.getDeletionsFrac()){
					bestBaseFrac = pos.getDeletionsFrac();
					bestBase = '-';
				}
			}
		}
	}

	if (numberAboveCutOff == 0 || (numberAboveCutOff == 1 && !pars.addNonRefSnps_)) {
		//check if adding non ref snps which would not be variants > 1
		//this check is done below
		return VarPosRes{false, numberAboveCutOff};
	}
	if(pars.onlyBiallelicVariants_ && numberAboveCutOff > 2){
		return VarPosRes{false, numberAboveCutOff};
	}
	if(pars.addNonRefSnps_ && numberAboveCutOff == 1 && refBase == bestBase){
		//if the number above cut off is 1 and it equals the best base (is reference return false)
		return VarPosRes{false, numberAboveCutOff};
	}
	return VarPosRes{true, numberAboveCutOff};
}

table VariantGatherer::gatherVariantsFromChromFile(
		const std::string & filename,
		const std::string & sampleName,
		const GathVarPars pars,
		TwoBit::TwoBitFile & twoBitfile) {
	BGZFCPP bgzInFile(filename, "r");
	std::string genomicSeq = "";
	uint32_t currentId = std::numeric_limits<uint32_t>::max();
	table outTab(
			concatVecs(VecStr { "sample", "alleleCount" }, SlimCounterPos::hqDetailHeader()));
	std::vector<char> bases { 'A', 'C', 'G', 'T' };
	std::vector<uint32_t> d(SlimCounterPos::NumOfElements);

	while (bgzInFile.readVec(d)) {
		SlimCounterPos pos(d);

		//check depth
		//if(pos.getHqBaseTotal() <= depthCutOff_){
		// continue;
		//}
		//check depth

		if (pos.getHqEventsTotal() <= depthCutOff_) {
			continue;
		}
		uint32_t refPos = d[1];
		auto rName = refInfos_[d[0]].RefName;
		if (currentId != d[0]) {
			twoBitfile[rName]->getSequence(genomicSeq);
			currentId = d[0];
		}
		const char refBase = genomicSeq[refPos];

		auto res = determineIfPosVariant(pos, refBase, pars);
		if(res.variant_){
			for (const auto & row : pos.hqDetailOutVec(rName, refPos,
					refBase)) {
				outTab.addRow(
						concatVecs(toVecStr(sampleName, res.numberOfAlleles_), row));
			}
		}
	}
	return outTab;
}



}  // namespace njhseq

