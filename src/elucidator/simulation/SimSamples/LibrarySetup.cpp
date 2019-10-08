/*
 * LibrarySetup.cpp
 *
 *  Created on: May 2, 2019
 *      Author: nicholashathaway
 */



#include "LibrarySetup.hpp"

namespace njhseq {


void ControlPopulation::Sample::ExperimentRun::setGenomesSampled(const std::map<std::string, Strain> & strainsExpected){
	hapAbundGenomesSampled_.clear();
	hapRegionAmplified_.clear();
	std::vector<seqInfo> seqs;
	for(const auto & strain : strainsExpected){
		seqInfo seq(strain.second.name_);
		seq.cnt_ = strain.second.relativeAbundance_;
		seq.frac_ = strain.second.relativeAbundance_;
		seqs.emplace_back(seq);
	}
	auto genomesSampled = PCRSimulator::randomlySampleGenomes(seqs, expAmounts_.totalGenomesSampled_);
	std::map<std::string, std::map<std::string, double>> subRegionsSampled;
	for (const auto & strainSampled : genomesSampled) {
//		std::cout << strainSampled.seqBase_.name_ << "\t" << strainSampled.genomeCnt_ << std::endl;
		Strain strain(strainSampled.seqBase_.name_,
				strainsExpected.at(strainSampled.seqBase_.name_).genomicRegionToHapNames_,
				strainSampled.genomeCnt_);
		strain.meta_ = strainsExpected.at(strainSampled.seqBase_.name_).meta_;
		hapAbundGenomesSampled_.emplace(strainSampled.seqBase_.name_, strain);

		for(const auto & hap : strain.genomicRegionToHapNames_){
			subRegionsSampled[hap.first][hap.second] += strain.relativeAbundance_;
		}
	}

	for(const auto & subRegionSampled : subRegionsSampled){
		hapRegionAmplified_.emplace(subRegionSampled.first, GenomicRegionAmp(subRegionSampled.first, subRegionSampled.second));
	}

//	for(const auto & regionInfo : hapRegionAmplified_){
//		std::cout << regionInfo.first << std::endl;
//		uint64_t total = 0;
//		for(const auto & subHap: regionInfo.second.hapAbundSampled_){
//			std::cout << '\t' << subHap.first << "\t" << subHap.second << std::endl;
//			total += subHap.second;
//		}
//		std::cout << "\t" << total << std::endl;
//	}
}



table LibrarySetup::createAbundanceTable() const{
	table ret(VecStr{"library", "sample", "mix", "hap",
		"abundanceRaw", "totalRawAbundance", "Percentage",
		"genomesSampled", "totalGenomes", "PercentageGenomes"});
	for(const auto & sample : samples_){
		for(const auto & mix : sample.second->mixtures_){
			double totalRawAbundance = 0.0;
			for(const auto & hap : mix.second->expectedAbundances_){
				totalRawAbundance += hap.second;
			}
			double totalGenomes = 0;
			for(const auto & hap : mix.second->genomeCounts_){
				totalGenomes += hap.second;
			}
			for(const auto & hap : mix.second->expectedAbundances_){
				ret.addRow(name_, sample.second->name_, mix.second->name_, hap.first,
						hap.second, totalRawAbundance, 100 * (hap.second/totalRawAbundance),
						mix.second->genomeCounts_[hap.first], totalGenomes, 100 * (mix.second->genomeCounts_[hap.first]/totalGenomes)
						);
			}
		}
	}
	return ret;
}



}  // namespace njhseq

