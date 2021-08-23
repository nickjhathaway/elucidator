/*
 * programWrappers_parsePrimer3Output.cpp
 *
 *  Created on: Aug 9, 2017
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

#include "programWrappers.hpp"
#include "elucidator/objects/BioDataObject.h"
#include "elucidator/BamToolsUtils.h"

#include <SeekDeep/objects/TarAmpSetupUtils/PrimersAndMids.hpp>

#include <njhseq/concurrency/PairwisePairFactory.hpp>

namespace njhseq {



int programWrapperRunner::genPossiblePrimersWithPrimer3(const njh::progutils::CmdArgs & inputCommands){

	uint32_t PRIMER_MAX_SIZE = 35;
	uint32_t PRIMER_MIN_SIZE = 18;
	uint32_t PRIMER_OPT_SIZE = 33;
	uint32_t minSize = 175;
	uint32_t maxSize = 290;
	uint32_t PRIMER_NUM_RETURN = 50;

	std::string task = "generic";

	bfs::path excludeRegionsFnp;
	bfs::path regionsOfInterestFnp;
	bfs::path bedFnp;
	bfs::path twoBitFnp;
	bfs::path primer3ConfigPath = "/usr/local/Cellar/primer3/2.4.0/share/primer3/primer3_config/";
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(excludeRegionsFnp, "--excludeRegions", "exclude Regions");
	setUp.setOption(regionsOfInterestFnp, "--regionsOfInterest", "regions Of Interest");

	setUp.setOption(PRIMER_MAX_SIZE, "--PRIMER_MAX_SIZE", "PRIMER MAX SIZE");
	setUp.setOption(PRIMER_MIN_SIZE, "--PRIMER_MIN_SIZE", "PRIMER MIN SIZE");
	setUp.setOption(PRIMER_OPT_SIZE, "--PRIMER_OPT_SIZE", "PRIMER OPT SIZE");
	setUp.setOption(minSize, "--minSize", "min Size");
	setUp.setOption(maxSize, "--maxSize", "max Size");
	setUp.setOption(PRIMER_NUM_RETURN, "--PRIMER_NUM_RETURN", "PRIMER NUM RETURN");
	setUp.setOption(task, "--task", "primer picking task, examples include:generic(default), pick_sequencing_primers, pick_primer_list");


	uint32_t extendRegion = std::max<uint32_t>(PRIMER_MAX_SIZE + 5, std::round(maxSize * 1.5));



	setUp.setOption(bedFnp, "--bedFnp", "genomic locations", true);
	setUp.setOption(twoBitFnp, "--twoBit", "two Bit file", true);
	setUp.setOption(primer3ConfigPath, "--primer3ConfigPath", "primer3 Config Path", true);

	setUp.processDirectoryOutputName("genPrimers_TODAY", true);

	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	primer3ConfigPath = njh::appendAsNeededRet(primer3ConfigPath.string(), "/");

	std::vector<std::shared_ptr<Bed6RecordCore>> excludeRegions;
	std::vector<std::shared_ptr<Bed6RecordCore>> regionsOfInterest;
	if(bfs::exists(regionsOfInterestFnp)){
		regionsOfInterest = getBeds(regionsOfInterestFnp);
	}
	if(bfs::exists(excludeRegionsFnp)){
		excludeRegions = getBeds(excludeRegionsFnp);
	}



	auto regions = getBeds(bedFnp);

	std::unordered_map<std::string, std::shared_ptr<Bed6RecordCore>> nameToRegion;
	std::unordered_map<std::string, std::shared_ptr<Bed6RecordCore>> templateToRegion;

	TwoBit::TwoBitFile tReader(twoBitFnp);
	auto chromLens = tReader.getSeqLens();
	std::set<std::string> regionNames;
	for(const auto & reg : regions){
		regionNames.emplace(reg->name_);
	}

	BedUtility::coordSort(excludeRegions, false);
	BedUtility::coordSort(regions, false);
	BedUtility::coordSort(regionsOfInterest, false);

	OutputStream inputRegionsExpanded(njh::files::make_path(setUp.pars_.directoryName_, "inputRegionsExpanded.bed"));


	//SEQUENCE_TARGET, SEQUENCE_EXCLUDED_REGION
	{
		OutputStream primer3Input(njh::files::make_path(setUp.pars_.directoryName_, "primer3File.txt"));
		primer3Input << "PRIMER_THERMODYNAMIC_PARAMETERS_PATH=" << primer3ConfigPath.string() << std::endl;

		for(const auto & originalReg : regions){
			std::shared_ptr<Bed6RecordCore> reg = std::make_shared<Bed6RecordCore>(*originalReg);
			reg->strand_ = '+';
			std::vector<Bed6RecordCore> regionsOfInterest_withinRegion;
			std::vector<Bed6RecordCore> excludeRegions_withinRegion;
			std::vector<Bed6RecordCore> regionsOfInterest_relToRegion;
			std::vector<Bed6RecordCore> excludeRegions_relToRegion;
			{
				for(const auto & inter : regionsOfInterest){
					if(inter->chrom_ < reg->chrom_){
						//bother regions are sorted if we haven't reached this region's chromosome yet
						continue;
					}
					if(inter->chrom_ > reg->chrom_){
						//both regions are sorted so if we run into this we can break
						break;
					}
					if(inter->chrom_ == reg->chrom_ && inter->chromStart_ >= reg->chromEnd_){
						//both regions are sorted so if we run into this we can break
						break;
					}
					if(inter->chrom_ == reg->chrom_ && inter->chromEnd_ < reg->chromStart_){
						continue;
					}
					if(reg->overlaps(*inter, 1) && inter->chromStart_ >= reg->chromStart_ && inter->chromEnd_ <=reg->chromEnd_){
						//take only completely within this region
						regionsOfInterest_withinRegion.emplace_back(*inter);
					}
				}
			}
			BedUtility::extendLeftRight(*reg, extendRegion, extendRegion, njh::mapAt(chromLens, reg->chrom_) );

			inputRegionsExpanded << reg->toDelimStrWithExtra() << std::endl;
			nameToRegion[reg->name_] = reg;

			{
				for(const auto & inter : excludeRegions){
					if(inter->chrom_ < reg->chrom_){
						//bother regions are sorted if we haven't reached this region's chromosome yet
						continue;
					}
					if(inter->chrom_ > reg->chrom_){
						//both regions are sorted so if we run into this we can break
						break;
					}
					if(inter->chrom_ == reg->chrom_ && inter->chromStart_ >= reg->chromEnd_){
						//both regions are sorted so if we run into this we can break
						break;
					}
					if(inter->chrom_ == reg->chrom_ && inter->chromEnd_ < reg->chromStart_){
						continue;
					}
					if(reg->overlaps(*inter, 1)){
						//check first if we already have this exact region for exclusion, this can happen if if input generated from variant calls and there is more than biallelic SNPs etc
						bool alreadyHave = false;
						for(const auto & already : excludeRegions_withinRegion){
							if(already.genUIDFromCoords() == inter->genUIDFromCoords()){
								alreadyHave = true;
								break;
							}
						}
						if(!alreadyHave){
							excludeRegions_withinRegion.emplace_back(*inter);
						}
					}
				}
				//book end the exclusion regions
				for(auto & exclude : excludeRegions_withinRegion){
					exclude.chromStart_ = std::max(exclude.chromStart_, reg->chromStart_);
					exclude.chromEnd_ = std::min(exclude.chromEnd_, reg->chromEnd_);
				}
			}

			//convert to regions relative to the template region
			for(const auto & roi : regionsOfInterest_withinRegion){
				auto modRegion = roi;
				modRegion.chromStart_ = modRegion.chromStart_ - reg->chromStart_;
				modRegion.chromEnd_ =   modRegion.chromEnd_   - reg->chromStart_;
				regionsOfInterest_relToRegion.emplace_back(modRegion);
			}
			for(const auto & exclude : excludeRegions_withinRegion){
				auto modRegion = exclude;
				modRegion.chromStart_ = modRegion.chromStart_ - reg->chromStart_;
				modRegion.chromEnd_ =   modRegion.chromEnd_   - reg->chromStart_;
				excludeRegions_relToRegion.emplace_back(modRegion);
			}

			auto seqTemplate = GenomicRegion(*reg).extractSeq(tReader);
			std::string seqID = reg->name_;
			std::stringstream defaultPars;
			defaultPars << "SEQUENCE_TEMPLATE=" << seqTemplate.seq_<< std::endl;
			defaultPars << "SEQUENCE_ID=" << seqID<< std::endl;
			defaultPars << "PRIMER_FIRST_BASE_INDEX=0" << std::endl;
			defaultPars << "PRIMER_TASK=" << task<< std::endl;
			defaultPars << "PRIMER_SEQUENCING_LEAD=5" << std::endl;
			defaultPars << "PRIMER_EXPLAIN_FLAG=1" << std::endl;
			defaultPars << "PRIMER_PRODUCT_SIZE_RANGE=" << minSize << "-" << maxSize << std::endl;
			defaultPars << "PRIMER_NUM_RETURN=" << PRIMER_NUM_RETURN << std::endl;
			defaultPars << "PRIMER_MAX_SIZE=" << PRIMER_MAX_SIZE << std::endl;
			defaultPars << "PRIMER_MIN_SIZE=" << PRIMER_MIN_SIZE << std::endl;
			defaultPars << "PRIMER_OPT_SIZE=" << PRIMER_OPT_SIZE << std::endl;

			if(!excludeRegions_relToRegion.empty()){
				std::string SEQUENCE_EXCLUDED_REGION = "";
				for(const auto & exclude : excludeRegions_relToRegion){
					if("" != SEQUENCE_EXCLUDED_REGION){
						SEQUENCE_EXCLUDED_REGION += " ";
					}
					SEQUENCE_EXCLUDED_REGION += njh::pasteAsStr(exclude.chromStart_, ",", exclude.chromEnd_ - exclude.chromStart_);
				}
				defaultPars << "SEQUENCE_EXCLUDED_REGION=" << SEQUENCE_EXCLUDED_REGION << std::endl;
			}
			if(!regionsOfInterest_relToRegion.empty()){
				for(const auto & roi : regionsOfInterest_relToRegion){
					std::string SEQUENCE_TARGET = njh::pasteAsStr(roi.chromStart_, ",", roi.chromEnd_ - roi.chromStart_);
					primer3Input << defaultPars.str();
					primer3Input << "SEQUENCE_TARGET=" << SEQUENCE_TARGET << std::endl;
					primer3Input << "="<< std::endl;
				}
			}else{
				primer3Input << defaultPars.str();
				primer3Input << "="<< std::endl;
			}
		}
	}
	std::string primer3Cmd = njh::pasteAsStr("cd ", setUp.pars_.directoryName_, " && ", "primer3_core ", "primer3File.txt > primer3_output.txt 2> primer3_error.txt");
	auto runOutput = njh::sys::run({primer3Cmd});
	OutputStream primer3RunLog(njh::files::make_path(setUp.pars_.directoryName_, "primer3RunLog.json"));
	primer3RunLog << runOutput.toJson() << std::endl;

	auto primer3ResultsFnp = njh::files::make_path(setUp.pars_.directoryName_, "primer3_output.txt");
	auto results = Primer3Results::parsePrimer3OutputResults(primer3ResultsFnp, true);

	OutputStream primer3ResultsOut(njh::files::make_path(setUp.pars_.directoryName_, "primer3_results.tab.txt"));

	OutputStream primer3ResultsPrimerLocs(njh::files::make_path(setUp.pars_.directoryName_, "primer3_results_primerLocs.bed"));
	OutputStream primer3ResultsInsertLocs(njh::files::make_path(setUp.pars_.directoryName_, "primer3_results_insertLocs.bed"));
	OutputStream primer3ResultsAmpLocs   (njh::files::make_path(setUp.pars_.directoryName_, "primer3_results_ampLocs.bed"));
	auto bestDir = njh::files::make_path(setUp.pars_.directoryName_, "bestPrimerPerTarget");
	njh::files::makeDir(njh::files::MkdirPar{bestDir});
	OutputStream best_primer3ResultsPrimerLocs(njh::files::make_path(bestDir, "primer3_results_primerLocs.bed"));
	OutputStream best_primer3ResultsInsertLocs(njh::files::make_path(bestDir, "primer3_results_insertLocs.bed"));
	OutputStream best_primer3ResultsAmpLocs   (njh::files::make_path(bestDir, "primer3_results_ampLocs.bed"));
	OutputStream best_primersTable   (njh::files::make_path(bestDir, "primer3_results_primers.tab.txt"));
	best_primersTable << "target\tforward\treverse" << std::endl;

	primer3ResultsOut << njh::conToStr(VecStr{"seqID", "chrom", "ampStart", "ampEnd", "primerPairName", "compl_any_th", "compl_end_th", "pair_penalty", "product_size", "product_gc_percent",
	"left_name", "left_seq", "left_start", "left_end", "left_size", "left_end_stability", "left_gc_percent", "left_hairpin_th", "left_penalty", "left_self_any_th", "left_self_end_th", "left_tm", "left_tm_hairpin_diff", "left_problems",
	"right_name", "right_seq", "right_start", "right_end", "right_size", "right_end_stability", "right_gc_percent", "right_hairpin_th", "right_penalty", "right_self_any_th", "right_self_end_th", "right_tm", "right_tm_hairpin_diff", "right_problems"
	}, "\t") << std::endl;
	VecStr noTargetsFor;
	std::stringstream noTargetsBed;

	uint32_t maxChromLen = 0;
	for(const auto & chrom : chromLens){
		if(chrom.second > maxChromLen){
			maxChromLen = chrom.second;
		}
	}

	std::vector<std::shared_ptr<njhseq::Primer3Results::PrimerPair>> allBestPrimers;
	std::unordered_set<std::string> alreadyHaveUID;
	for(const auto & res : results){
		const auto & region = nameToRegion[res->sequence_id_];
		bool targeted = !res->sequence_target_.empty();
		if(0 == res->primer_pair_num_returned_){
			//should log which attempts had 0 returned
			std::string targetedName = "";
			if (targeted) {
				targetedName = njh::pasteAsStr(region->chrom_, "-",
						region->chromStart_ + res->sequence_target_.front().start_, "-",
						region->chromStart_ + res->sequence_target_.front().start_ + res->sequence_target_.front().size_);
				noTargetsFor.emplace_back(targetedName);
				noTargetsBed << region->chrom_
						<< "\t" << region->chromStart_ + res->sequence_target_.front().start_
						<< "\t" << region->chromStart_ + res->sequence_target_.front().start_ + res->sequence_target_.front().size_
						<< "\t" << targetedName
						<< "\t" << res->sequence_target_.front().size_
						<< "\t" << "+" << std::endl;
			}
			continue;
		}

		if(targeted){
			auto bestPrimer = res->getLowestPenaltyPair();
			bool passBestPrimer = true;
			for(const auto & bestPrimerAlready :allBestPrimers){
				if(bestPrimer->left_.seq_ == bestPrimerAlready->left_.seq_ &&
					bestPrimer->right_.seq_ == bestPrimerAlready->right_.seq_){
					passBestPrimer = false;
					break;
				}
			}
			if(passBestPrimer){
				allBestPrimers.emplace_back(std::make_shared<njhseq::Primer3Results::PrimerPair>(*bestPrimer));
			}else{
				continue;
			}
			GenomicRegion insertRegion(
					res->sequence_id_,
					region->chrom_,
					region->chromStart_ + bestPrimer->left_.forwardOrientationPos_.start_ + bestPrimer->left_.forwardOrientationPos_.size_,
					region->chromStart_ + bestPrimer->right_.forwardOrientationPos_.start_,
					false);
			insertRegion.uid_ = njh::pasteAsStr(insertRegion.chrom_, "-",
					njh::leftPadNumStr<size_t>(insertRegion.start_, maxChromLen), "-",
					njh::leftPadNumStr<size_t>(insertRegion.end_, maxChromLen));
			uint32_t uidCount = 0;
			while(njh::in(insertRegion.uid_, alreadyHaveUID)){
				++uidCount;
				insertRegion.uid_ = njh::pasteAsStr(insertRegion.chrom_, "-",
						njh::leftPadNumStr<size_t>(insertRegion.start_, maxChromLen), "-",
						njh::leftPadNumStr<size_t>(insertRegion.end_, maxChromLen), "--", uidCount);
			}
			best_primer3ResultsInsertLocs << insertRegion.genBedRecordCore().toDelimStrWithExtra() << std::endl;
			GenomicRegion ampRegion(insertRegion.uid_,
					region->chrom_,
					region->chromStart_ + bestPrimer->getStart(),
					region->chromStart_ + bestPrimer->getEnd(),
					false);
			best_primer3ResultsAmpLocs << ampRegion.genBedRecordCore().toDelimStrWithExtra() << std::endl;

			GenomicRegion leftPrimerRegion(
					insertRegion.uid_ + "-forwardPrimer",
					region->chrom_,
					region->chromStart_ + bestPrimer->left_.forwardOrientationPos_.start_,
					region->chromStart_ + bestPrimer->left_.forwardOrientationPos_.start_ + bestPrimer->left_.forwardOrientationPos_.size_,
					false);
			GenomicRegion rightPrimerRegion(
					insertRegion.uid_ + "-reversePrimer",
					region->chrom_,
					region->chromStart_ + bestPrimer->right_.forwardOrientationPos_.start_,
					region->chromStart_ + bestPrimer->right_.forwardOrientationPos_.start_ + bestPrimer->right_.forwardOrientationPos_.size_,
					true);
			best_primer3ResultsPrimerLocs << leftPrimerRegion.genBedRecordCore().toDelimStrWithExtra() << std::endl;
			best_primer3ResultsPrimerLocs << rightPrimerRegion.genBedRecordCore().toDelimStrWithExtra() << std::endl;

			best_primersTable << insertRegion.uid_
					<< "\t" << bestPrimer->left_.seq_
					<< "\t" << bestPrimer->right_.seq_ << std::endl;
		}else{
			auto bestPrimers = res->getLowestPenaltyNonOverlappingPairs();
			for(const auto & bestPrimer : bestPrimers){
				bool passBestPrimer = true;
				for(const auto & bestPrimerAlready :allBestPrimers){
					if(bestPrimer->left_.seq_ == bestPrimerAlready->left_.seq_ &&
						bestPrimer->right_.seq_ == bestPrimerAlready->right_.seq_){
						passBestPrimer = false;
						break;
					}
				}
				if(passBestPrimer){
					allBestPrimers.emplace_back(std::make_shared<njhseq::Primer3Results::PrimerPair>(*bestPrimer));
				}else{
					continue;
				}
				GenomicRegion insertRegion(
						res->sequence_id_,
						region->chrom_,
						region->chromStart_ + bestPrimer->left_.forwardOrientationPos_.start_ + bestPrimer->left_.forwardOrientationPos_.size_,
						region->chromStart_ + bestPrimer->right_.forwardOrientationPos_.start_,
						false);
				insertRegion.uid_ = njh::pasteAsStr(insertRegion.chrom_, "-",
						njh::leftPadNumStr<size_t>(insertRegion.start_, maxChromLen), "-",
						njh::leftPadNumStr<size_t>(insertRegion.end_, maxChromLen));
				uint32_t uidCount = 0;
				while(njh::in(insertRegion.uid_, alreadyHaveUID)){
					++uidCount;
					insertRegion.uid_ = njh::pasteAsStr(insertRegion.chrom_, "-",
							njh::leftPadNumStr<size_t>(insertRegion.start_, maxChromLen), "-",
							njh::leftPadNumStr<size_t>(insertRegion.end_, maxChromLen), "--", uidCount);
				}
				best_primer3ResultsInsertLocs << insertRegion.genBedRecordCore().toDelimStrWithExtra() << std::endl;
				GenomicRegion ampRegion(insertRegion.uid_,
						region->chrom_,
						region->chromStart_ + bestPrimer->getStart(),
						region->chromStart_ + bestPrimer->getEnd(),
						false);
				best_primer3ResultsAmpLocs << ampRegion.genBedRecordCore().toDelimStrWithExtra() << std::endl;

				GenomicRegion leftPrimerRegion(
						insertRegion.uid_ + "-forwardPrimer",
						region->chrom_,
						region->chromStart_ + bestPrimer->left_.forwardOrientationPos_.start_,
						region->chromStart_ + bestPrimer->left_.forwardOrientationPos_.start_ + bestPrimer->left_.forwardOrientationPos_.size_,
						false);
				GenomicRegion rightPrimerRegion(
						insertRegion.uid_ + "-reversePrimer",
						region->chrom_,
						region->chromStart_ + bestPrimer->right_.forwardOrientationPos_.start_,
						region->chromStart_ + bestPrimer->right_.forwardOrientationPos_.start_ + bestPrimer->right_.forwardOrientationPos_.size_,
						true);
				best_primer3ResultsPrimerLocs << leftPrimerRegion.genBedRecordCore().toDelimStrWithExtra() << std::endl;
				best_primer3ResultsPrimerLocs << rightPrimerRegion.genBedRecordCore().toDelimStrWithExtra() << std::endl;

				best_primersTable << insertRegion.uid_
						<< "\t" << bestPrimer->left_.seq_
						<< "\t" << bestPrimer->right_.seq_ << std::endl;
			}
		}

		for(const auto & primerPair : res->primerPairs_){
			uint32_t targetStart = primerPair->left_.forwardOrientationPos_.start_;
			uint32_t targetEnd = primerPair->right_.forwardOrientationPos_.start_ + primerPair->right_.forwardOrientationPos_.size_;
			std::string target = res->sequence_template_.substr(targetStart, targetEnd - targetStart);
			uint32_t chromStart = region->chromStart_ + targetStart;
			uint32_t chromEnd = region->chromStart_ + targetEnd;

			DNABaseCounter counter;
			counter.increase(target);
			counter.setFractions();
			bool targeted = !res->sequence_target_.empty();
			//since above we split out targets to be just a single at a time will just take the front
			if(res->sequence_target_.size() > 1){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << ", target size shouldn't be more than 1, res->sequence_target_.size(): " << res->sequence_target_.size()<< "\n";
				throw std::runtime_error{ss.str()};
			}
			std::string targetedName = "";
			if (targeted) {
				targetedName = njh::pasteAsStr(region->chrom_, "-",
						region->chromStart_ + res->sequence_target_.front().start_, "-",
						region->chromStart_ + res->sequence_target_.front().start_
								+ res->sequence_target_.front().size_);
			}
			{
				MetaDataInName ampMeta;
				ampMeta.addMeta("SeqID", res->sequence_id_);
				ampMeta.addMeta("PrimerID", primerPair->name_);
				if(targeted){
					ampMeta.addMeta("target", targetedName);
				}
				primer3ResultsAmpLocs << region->chrom_
						<< "\t" << region->chromStart_ + targetStart
						<< "\t" << region->chromStart_ + targetEnd
						<< "\t" << ampMeta.createMetaName()
						<< "\t" << targetEnd - targetStart
						<< "\t" << "+" << std::endl;
			}

			{
				MetaDataInName insertMeta;
				insertMeta.addMeta("SeqID", res->sequence_id_);
				insertMeta.addMeta("PrimerID", primerPair->name_);
				if(targeted){
					insertMeta.addMeta("target", targetedName);
				}
				primer3ResultsInsertLocs << region->chrom_
						<< "\t" << region->chromStart_ + primerPair->left_.forwardOrientationPos_.start_ + primerPair->left_.forwardOrientationPos_.size_
						<< "\t" << region->chromStart_ + primerPair->right_.forwardOrientationPos_.start_
						<< "\t" << insertMeta.createMetaName()
						<< "\t" << primerPair->right_.forwardOrientationPos_.start_ - (primerPair->left_.forwardOrientationPos_.start_ + primerPair->left_.forwardOrientationPos_.size_)
						<< "\t" << "+" << std::endl;
			}

			{
				MetaDataInName leftMeta;
				leftMeta.addMeta("SeqID", res->sequence_id_);
				leftMeta.addMeta("PrimerID", primerPair->name_);
				leftMeta.addMeta("PrimerName", primerPair->left_.name_);
				if(targeted){
					leftMeta.addMeta("target", targetedName);
				}
				primer3ResultsPrimerLocs << region->chrom_
						<< "\t" << region->chromStart_ + primerPair->left_.forwardOrientationPos_.start_
						<< "\t" << region->chromStart_ + primerPair->left_.forwardOrientationPos_.start_ + primerPair->left_.forwardOrientationPos_.size_
						<< "\t" << leftMeta.createMetaName()
						<< "\t" << primerPair->left_.forwardOrientationPos_.size_
						<< "\t" << "+" << std::endl;

				MetaDataInName rightMeta;
				rightMeta.addMeta("SeqID", res->sequence_id_);
				rightMeta.addMeta("PrimerID", primerPair->name_);
				rightMeta.addMeta("PrimerName", primerPair->right_.name_);
				if(targeted){
					rightMeta.addMeta("target", targetedName);
				}
				primer3ResultsPrimerLocs << region->chrom_
						<< "\t" << region->chromStart_ + primerPair->right_.forwardOrientationPos_.start_
						<< "\t" << region->chromStart_ + primerPair->right_.forwardOrientationPos_.start_ + primerPair->right_.forwardOrientationPos_.size_
						<< "\t" << rightMeta.createMetaName()
						<< "\t" << primerPair->right_.forwardOrientationPos_.size_
						<< "\t" << "-" << std::endl;
			}

			primer3ResultsOut << njh::conToStr(toVecStr(
				res->sequence_id_ + (targeted ? "--" + targetedName : ""),
				region->chrom_,



				chromStart,
				chromEnd,
				primerPair->name_,
				primerPair->compl_any_th_,
				primerPair->compl_end_th_,
				primerPair->penalty_,
				primerPair->product_size_,
				roundDecPlaces(counter.calcGcContent() * 100, 2),
				primerPair->left_.name_,
				primerPair->left_.seq_,
				region->chromStart_ + primerPair->left_.forwardOrientationPos_.start_,
				region->chromStart_ + primerPair->left_.forwardOrientationPos_.start_ + primerPair->left_.forwardOrientationPos_.size_,
				primerPair->left_.forwardOrientationPos_.size_,
				primerPair->left_.end_stability_,
				primerPair->left_.gc_percent_,
				primerPair->left_.hairpin_th_,
				primerPair->left_.penalty_,
				primerPair->left_.self_any_th_,
				primerPair->left_.self_end_th_,
				primerPair->left_.tm_,
				primerPair->left_.tm_ - primerPair->left_.hairpin_th_,
				njh::conToStr(primerPair->left_.problems_, ";"),

				primerPair->right_.name_,
				primerPair->right_.seq_,
				region->chromStart_ + primerPair->right_.forwardOrientationPos_.start_,
				region->chromStart_ + primerPair->right_.forwardOrientationPos_.start_ + primerPair->right_.forwardOrientationPos_.size_,
				primerPair->right_.forwardOrientationPos_.size_,
				primerPair->right_.end_stability_,
				primerPair->right_.gc_percent_,
				primerPair->right_.hairpin_th_,
				primerPair->right_.penalty_,
				primerPair->right_.self_any_th_,
				primerPair->right_.self_end_th_,
				primerPair->right_.tm_,
				primerPair->right_.tm_ - primerPair->right_.hairpin_th_,
				njh::conToStr(primerPair->right_.problems_, ";")
			), "\t") << std::endl;
		}
	}
	if(!noTargetsFor.empty()){
		OutputStream noResultsForRegions(njh::files::make_path(setUp.pars_.directoryName_, "noTargetsForRegionsOfInterest.bed"));
		noResultsForRegions << noTargetsBed.str();
	}

	OutputStream primer3ResultsSummaryOut(
			njh::files::make_path(setUp.pars_.directoryName_,
					"primer3_results_summary.tab.txt"));
	primer3ResultsSummaryOut << njh::conToStr(VecStr { "seqID", "chrom", "start",
			"end", "len", "targeted", "targetedID", "targetedChrom", "targetedStart",
			"targetedEnd", "targetedSize", "excludedRegions", "excludedRegionsIDs", "primerTask",
			"pairNumReturned", "leftNumReturned", "rightNumReurned", "pairExplained",
			"leftExplained", "rightExplained", "warnings" }, "\t") << std::endl;

	for (const auto &res : results) {
		const auto &region = nameToRegion[res->sequence_id_];
		bool targeted = !res->sequence_target_.empty();
		std::string targetedName = "";
		if (targeted) {
			targetedName = njh::pasteAsStr(region->chrom_, "-",
					region->chromStart_ + res->sequence_target_.front().start_, "-",
					region->chromStart_ + res->sequence_target_.front().start_
							+ res->sequence_target_.front().size_);
		}
		std::string seqExcludedRegions;
		std::string seqExcludedRegionsIDs;
		for (const auto &exclude : res->sequence_excluded_region_) {
			std::string excludeRel = njh::pasteAsStr(exclude.start_, exclude.size_,";");
			std::string excludeID = njh::pasteAsStr(region->chrom_, "-",
					region->chromStart_ + exclude.start_, "-",
					region->chromStart_ + exclude.start_ + exclude.size_, ";");
			seqExcludedRegions += excludeRel;
			seqExcludedRegionsIDs += excludeID;
		}
		primer3ResultsSummaryOut << njh::conToStr(toVecStr(
			res->sequence_id_,
			region->chrom_,
			region->chromStart_,
			region->chromEnd_,
			region->length(),
			njh::boolToStr(targeted),
			(targeted ? targetedName: std::string("NA")),
			region->chrom_,
			region->chromStart_ + res->sequence_target_.front().start_,
			region->chromStart_ + res->sequence_target_.front().start_ + res->sequence_target_.front().size_,
			res->sequence_target_.front().size_,
			seqExcludedRegions,
			seqExcludedRegionsIDs,
			res->primer_task_,
			res->primer_pair_num_returned_,
			res->primer_left_num_returned_,
			res->primer_right_num_returned_,
			res->primer_pair_explain_,
			res->primer_left_explain_,
			res->primer_right_explain_,
			njh::conToStr(res->warnings_, ";")
		), "\t") << std::endl;
	}
	return 0;
}



int programWrapperRunner::testPrimersWithPrimer3(const njh::progutils::CmdArgs & inputCommands){
	comparison allowablePairwisePrimerComps;
	allowablePairwisePrimerComps.hqMismatches_ = 1;
	bool getPairwisePrimerComps = false;
	uint32_t PRIMER_MAX_SIZE = 35;
	uint32_t PRIMER_MIN_SIZE = 18;
	uint32_t PRIMER_OPT_SIZE = 33;
	uint32_t PRIMER_OPT_GC_PERCENT = 50;

	bfs::path bedFnp;
	bfs::path twoBitFnp;
	bfs::path primersFnp;
	bfs::path primer3ConfigPath = "/usr/local/Cellar/primer3/2.4.0/share/primer3/primer3_config/";
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(PRIMER_MAX_SIZE, "--PRIMER_MAX_SIZE", "PRIMER MAX SIZE");
	setUp.setOption(PRIMER_MIN_SIZE, "--PRIMER_MIN_SIZE", "PRIMER MIN SIZE");
	setUp.setOption(PRIMER_OPT_SIZE, "--PRIMER_OPT_SIZE", "PRIMER OPT SIZE");
	setUp.setOption(PRIMER_OPT_GC_PERCENT, "--PRIMER_OPT_GC_PERCENT", "PRIMER OPT GC PERCENT");

	setUp.setOption(bedFnp, "--bedFnp", "genomic locations", true);
	setUp.setOption(twoBitFnp, "--twoBit", "two Bit file", true);
	setUp.setOption(primersFnp, "--primers", "primers to test", true);
	setUp.setOption(primer3ConfigPath, "--primer3ConfigPath", "primer3 Config Path", true);
	setUp.setOption(getPairwisePrimerComps, "--getPairwisePrimerComps", "get Pairwise Primer Comps");
	setUp.processComparison(allowablePairwisePrimerComps, "--pairwiseComp");


	setUp.processDirectoryOutputName("testPrimers_TODAY", true);

	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	primer3ConfigPath = njh::appendAsNeededRet(primer3ConfigPath.string(), "/");

	auto regions = getBeds(bedFnp);
	PrimersAndMids primers(primersFnp);
	primers.initPrimerDeterminator();
	TwoBit::TwoBitFile tReader(twoBitFnp);

	std::set<std::string> targetNames = njh::getSetOfMapKeys(primers.targets_);
	std::set<std::string> regionNames;
	for(const auto & reg : regions){
		regionNames.emplace(reg->name_);
	}
	VecStr missingFromTargetNames;
	std::set_difference(
			targetNames.begin(), targetNames.end(),
			regionNames.begin(), regionNames.end(),
			std::back_inserter(missingFromTargetNames));

	VecStr missingFromRegionNames;
	std::set_difference(
			regionNames.begin(), regionNames.end(),
			targetNames.begin(), targetNames.end(),
			std::back_inserter(missingFromRegionNames));
	if(!missingFromRegionNames.empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "missing the following targets from " << bedFnp << "\n";
		ss << njh::conToStr(missingFromRegionNames, ",") << std::endl;
		throw std::runtime_error{ss.str()};
	}
	if(!missingFromTargetNames.empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "missing the following regions from " << primersFnp << "\n";
		ss << njh::conToStr(missingFromTargetNames, ",") << std::endl;

		throw std::runtime_error{ss.str()};
	}



	std::unordered_map<std::string, double> gcContentForRegion;
	{
		OutputStream primer3Input(njh::files::make_path(setUp.pars_.directoryName_, "primer3File.txt"));
		primer3Input << "PRIMER_THERMODYNAMIC_PARAMETERS_PATH=" << primer3ConfigPath.string() << std::endl;

		for(const auto & reg : regions){
			auto seqTemplate = GenomicRegion(*reg).extractSeq(tReader);
			DNABaseCounter counter;
			counter.increase(seqTemplate.seq_);
			counter.setFractions();
			gcContentForRegion[reg->name_] = roundDecPlaces(counter.calcGcContent() * 100, 2);
			std::string seqID = reg->name_;
			uint32_t seqCount = 0;
			for(const auto & fwd : njh::mapAt(primers.targets_, reg->name_).info_.fwds_){
				for(const auto & rev : njh::mapAt(primers.targets_, reg->name_).info_.revs_){
					if(seqCount > 0){
						seqID = njh::pasteAsStr(seqID, ".", seqCount);
					}
					std::string leftPrimer = fwd.primer_;
					std::string rightPrimer = rev.primer_;
					primer3Input << "SEQUENCE_TEMPLATE=" << seqTemplate.seq_<< std::endl;
					primer3Input << "SEQUENCE_ID=" << seqID<< std::endl;
					primer3Input << "SEQUENCE_PRIMER=" << leftPrimer<< std::endl;
					primer3Input << "SEQUENCE_PRIMER_REVCOMP=" << rightPrimer<< std::endl;
					primer3Input << "PRIMER_PICK_ANYWAY=1" << std::endl;
					primer3Input << "PRIMER_TASK=check_primers" << std::endl;
					primer3Input << "PRIMER_FIRST_BASE_INDEX=0" << std::endl;
					primer3Input << "PRIMER_EXPLAIN_FLAG=1" << std::endl;

					primer3Input << "PRIMER_PRODUCT_SIZE_RANGE=" << len(seqTemplate) << "-" << len(seqTemplate ) + 1 << std::endl;
					primer3Input << "PRIMER_OPT_GC_PERCENT=" << PRIMER_OPT_GC_PERCENT << std::endl;
					primer3Input << "PRIMER_MAX_SIZE=" << PRIMER_MAX_SIZE << std::endl;
					primer3Input << "PRIMER_MIN_SIZE=" << PRIMER_MIN_SIZE << std::endl;
					primer3Input << "PRIMER_OPT_SIZE=" << PRIMER_OPT_SIZE << std::endl;
					primer3Input << "="<< std::endl;
					++seqCount;
				}
			}
		}
	}
	std::string primer3Cmd = njh::pasteAsStr("cd ", setUp.pars_.directoryName_, " && ", "primer3_core ", "primer3File.txt > primer3_output.txt 2> primer3_error.txt");
	auto runOutput = njh::sys::run({primer3Cmd});
	OutputStream primer3RunLog(njh::files::make_path(setUp.pars_.directoryName_, "primer3RunLog.json"));
	primer3RunLog << runOutput.toJson() << std::endl;

	auto primer3ResultsFnp = njh::files::make_path(setUp.pars_.directoryName_, "primer3_output.txt");
	auto results = Primer3Results::parsePrimer3OutputResults(primer3ResultsFnp, true);

	OutputStream primer3ResultsOut(njh::files::make_path(setUp.pars_.directoryName_, "primer3_results.tab.txt"));

	primer3ResultsOut << njh::conToStr(VecStr{"seqID", "primerPairName", "compl_any_th", "compl_end_th", "pair_penalty", "pair_penalty_noSize", "product_size", "product_gc_percent",
	"left_name", "left_seq", "left_start", "left_end", "left_size", "left_end_stability", "left_gc_percent", "left_hairpin_th", "left_penalty", "left_penalty_noSize", "left_self_any_th", "left_self_end_th", "left_tm", "left_tm_hairpin_diff", "left_problems",
	"right_name", "right_seq", "right_start", "right_end", "right_size", "right_end_stability", "right_gc_percent", "right_hairpin_th", "right_penalty","right_penalty_noSize", "right_self_any_th", "right_self_end_th", "right_tm", "right_tm_hairpin_diff", "right_problems"
	}, "\t") << std::endl;
	for(const auto & res : results){
		for(const auto & primerPair : res->primerPairs_){

			double leftPenaltyWithOutSize = primerPair->left_.penalty_ - uAbsdiff(primerPair->left_.forwardOrientationPos_.size_, PRIMER_OPT_SIZE);
			double rightPenaltyWithOutSize = primerPair->right_.penalty_ - uAbsdiff(primerPair->right_.forwardOrientationPos_.size_, PRIMER_OPT_SIZE);
			double pairPenaltyWithOutSize = leftPenaltyWithOutSize + rightPenaltyWithOutSize;

			primer3ResultsOut << njh::conToStr(toVecStr(
				res->sequence_id_,
				primerPair->name_,

				primerPair->compl_any_th_,
				primerPair->compl_end_th_,
				primerPair->penalty_,
				pairPenaltyWithOutSize,
				primerPair->product_size_,
				gcContentForRegion[res->sequence_id_],

				primerPair->left_.name_,
				primerPair->left_.seq_,
				primerPair->left_.forwardOrientationPos_.start_,
				primerPair->left_.forwardOrientationPos_.start_ + primerPair->left_.forwardOrientationPos_.size_,
				primerPair->left_.forwardOrientationPos_.size_,
				primerPair->left_.end_stability_,
				primerPair->left_.gc_percent_,
				primerPair->left_.hairpin_th_,
				primerPair->left_.penalty_,
				leftPenaltyWithOutSize,
				primerPair->left_.self_any_th_,
				primerPair->left_.self_end_th_,
				primerPair->left_.tm_,
				primerPair->left_.tm_ - primerPair->left_.hairpin_th_,
				njh::conToStr(primerPair->left_.problems_, ";"),

				primerPair->right_.name_,
				primerPair->right_.seq_,
				primerPair->right_.forwardOrientationPos_.start_,
				primerPair->right_.forwardOrientationPos_.start_ + primerPair->right_.forwardOrientationPos_.size_,
				primerPair->right_.forwardOrientationPos_.size_,
				primerPair->right_.end_stability_,
				primerPair->right_.gc_percent_,
				primerPair->right_.hairpin_th_,
				primerPair->right_.penalty_,
				rightPenaltyWithOutSize,
				primerPair->right_.self_any_th_,
				primerPair->right_.self_end_th_,
				primerPair->right_.tm_,
				primerPair->right_.tm_ - primerPair->right_.hairpin_th_,
				njh::conToStr(primerPair->right_.problems_, ";")
			), "\t") << std::endl;
		}
	}
	if(getPairwisePrimerComps){

	//  std::cout << __FILE__ << " " << __LINE__ << std::endl;
		auto maxPrimerLen = primers.pDeterminator_->getMaxPrimerSize();
	//	std::cout << "maxPrimerLen: " << maxPrimerLen << std::endl;

		aligner alnObject(maxPrimerLen, gapScoringParameters(5,1), substituteMatrix::createDegenScoreMatrixLessN(2, -2));
		//alnObject.countEndGaps_ = true;
		PairwisePairFactory pFactor(primers.pDeterminator_->primers_.size());

		auto primerNames = getVectorOfMapKeys(primers.pDeterminator_->primers_);

		PairwisePairFactory::PairwisePair pair;

		OutputStream primerPairWiseCompsOut(njh::files::make_path(setUp.pars_.directoryName_, "primerPairwiseComps.txt"));

		primerPairWiseCompsOut << "target1\ttarget1Primer\ttarget2\ttarget2Primer\tscore\tidentities\tpass" << std::endl;

		while(pFactor.setNextPair(pair)){
			alnObject.alignCacheGlobalDiag(
					primers.pDeterminator_->primers_.at(primerNames[pair.row_]).forwardPrimerRaw_,
					primers.pDeterminator_->primers_.at(primerNames[pair.col_]).forwardPrimerRaw_);
			alnObject.profilePrimerAlignment(
					primers.pDeterminator_->primers_.at(primerNames[pair.row_]).forwardPrimerRaw_,
					primers.pDeterminator_->primers_.at(primerNames[pair.col_]).forwardPrimerRaw_);
			primerPairWiseCompsOut << primerNames[pair.row_] << "\t" << "forward-5-3" << "\t" << primerNames[pair.col_] << "\t" << "forward-5-3";
			primerPairWiseCompsOut << "\t" << alnObject.comp_.distances_.eventBasedIdentity_ << "\t" <<  alnObject.comp_.distances_.query_.identities_ << "\t" << njh::boolToStr(allowablePairwisePrimerComps.passErrorProfile(alnObject.comp_))<< std::endl;

			alnObject.alignCacheGlobalDiag(
					primers.pDeterminator_->primers_.at(primerNames[pair.row_]).forwardPrimerRaw_,
					seqUtil::reverseComplement(primers.pDeterminator_->primers_.at(primerNames[pair.col_]).forwardPrimerRaw_, "DNA"));
			alnObject.profilePrimerAlignment(
					primers.pDeterminator_->primers_.at(primerNames[pair.row_]).forwardPrimerRaw_,
					seqUtil::reverseComplement(primers.pDeterminator_->primers_.at(primerNames[pair.col_]).forwardPrimerRaw_, "DNA"));
			primerPairWiseCompsOut << primerNames[pair.row_] << "\t" << "forward-5-3" << "\t" << primerNames[pair.col_] << "\t" << "forward-3-5";
			primerPairWiseCompsOut << "\t" << alnObject.comp_.distances_.eventBasedIdentity_ << "\t" <<  alnObject.comp_.distances_.query_.identities_ << "\t" << njh::boolToStr(allowablePairwisePrimerComps.passErrorProfile(alnObject.comp_)) << std::endl;

			alnObject.alignCacheGlobalDiag(
					primers.pDeterminator_->primers_.at(primerNames[pair.row_]).forwardPrimerRaw_,
					primers.pDeterminator_->primers_.at(primerNames[pair.col_]).reversePrimerRaw_);
			alnObject.profilePrimerAlignment(
					primers.pDeterminator_->primers_.at(primerNames[pair.row_]).forwardPrimerRaw_,
					primers.pDeterminator_->primers_.at(primerNames[pair.col_]).reversePrimerRaw_);
			primerPairWiseCompsOut << primerNames[pair.row_] << "\t" << "forward-5-3" << "\t" << primerNames[pair.col_] << "\t" << "reverse-5-3";
			primerPairWiseCompsOut << "\t" << alnObject.comp_.distances_.eventBasedIdentity_ << "\t" <<  alnObject.comp_.distances_.query_.identities_ << "\t" << njh::boolToStr(allowablePairwisePrimerComps.passErrorProfile(alnObject.comp_)) << std::endl;

			alnObject.alignCacheGlobalDiag(
					primers.pDeterminator_->primers_.at(primerNames[pair.row_]).forwardPrimerRaw_,
					seqUtil::reverseComplement(primers.pDeterminator_->primers_.at(primerNames[pair.col_]).reversePrimerRaw_, "DNA"));
			alnObject.profilePrimerAlignment(
					primers.pDeterminator_->primers_.at(primerNames[pair.row_]).forwardPrimerRaw_,
					seqUtil::reverseComplement(primers.pDeterminator_->primers_.at(primerNames[pair.col_]).reversePrimerRaw_, "DNA"));
			primerPairWiseCompsOut << primerNames[pair.row_] << "\t" << "forward-5-3" << "\t" << primerNames[pair.col_] << "\t" << "reverse-3-5";
			primerPairWiseCompsOut << "\t" << alnObject.comp_.distances_.eventBasedIdentity_ << "\t" <<  alnObject.comp_.distances_.query_.identities_ << "\t" << njh::boolToStr(allowablePairwisePrimerComps.passErrorProfile(alnObject.comp_)) << std::endl;

			alnObject.alignCacheGlobalDiag(
					primers.pDeterminator_->primers_.at(primerNames[pair.row_]).reversePrimerRaw_,
					primers.pDeterminator_->primers_.at(primerNames[pair.col_]).forwardPrimerRaw_);
			alnObject.profilePrimerAlignment(
					primers.pDeterminator_->primers_.at(primerNames[pair.row_]).reversePrimerRaw_,
					primers.pDeterminator_->primers_.at(primerNames[pair.col_]).forwardPrimerRaw_);
			primerPairWiseCompsOut << primerNames[pair.row_] << "\t" << "reverse-5-3" << "\t" << primerNames[pair.col_] << "\t" << "forward-5-3";
			primerPairWiseCompsOut << "\t" << alnObject.comp_.distances_.eventBasedIdentity_ << "\t" <<  alnObject.comp_.distances_.query_.identities_ << "\t" << njh::boolToStr(allowablePairwisePrimerComps.passErrorProfile(alnObject.comp_)) << std::endl;

			alnObject.alignCacheGlobalDiag(
					primers.pDeterminator_->primers_.at(primerNames[pair.row_]).reversePrimerRaw_,
					seqUtil::reverseComplement(primers.pDeterminator_->primers_.at(primerNames[pair.col_]).forwardPrimerRaw_, "DNA"));
			alnObject.profilePrimerAlignment(
					primers.pDeterminator_->primers_.at(primerNames[pair.row_]).reversePrimerRaw_,
					seqUtil::reverseComplement(primers.pDeterminator_->primers_.at(primerNames[pair.col_]).forwardPrimerRaw_, "DNA"));
			primerPairWiseCompsOut << primerNames[pair.row_] << "\t" << "reverse-5-3" << "\t" << primerNames[pair.col_] << "\t" << "forward-3-5";
			primerPairWiseCompsOut << "\t" << alnObject.comp_.distances_.eventBasedIdentity_ << "\t" <<  alnObject.comp_.distances_.query_.identities_ << "\t" << njh::boolToStr(allowablePairwisePrimerComps.passErrorProfile(alnObject.comp_)) << std::endl;

			alnObject.alignCacheGlobalDiag(
					primers.pDeterminator_->primers_.at(primerNames[pair.row_]).reversePrimerRaw_,
					primers.pDeterminator_->primers_.at(primerNames[pair.col_]).reversePrimerRaw_);
			alnObject.profilePrimerAlignment(
					primers.pDeterminator_->primers_.at(primerNames[pair.row_]).reversePrimerRaw_,
					primers.pDeterminator_->primers_.at(primerNames[pair.col_]).reversePrimerRaw_);
			primerPairWiseCompsOut << primerNames[pair.row_] << "\t" << "reverse-5-3" << "\t" << primerNames[pair.col_] << "\t" << "reverse-5-3";
			primerPairWiseCompsOut << "\t" << alnObject.comp_.distances_.eventBasedIdentity_ << "\t" <<  alnObject.comp_.distances_.query_.identities_ << "\t" << njh::boolToStr(allowablePairwisePrimerComps.passErrorProfile(alnObject.comp_)) << std::endl;

			alnObject.alignCacheGlobalDiag(
					primers.pDeterminator_->primers_.at(primerNames[pair.row_]).reversePrimerRaw_,
					seqUtil::reverseComplement(primers.pDeterminator_->primers_.at(primerNames[pair.col_]).reversePrimerRaw_, "DNA"));
			alnObject.profilePrimerAlignment(
					primers.pDeterminator_->primers_.at(primerNames[pair.row_]).reversePrimerRaw_,
					seqUtil::reverseComplement(primers.pDeterminator_->primers_.at(primerNames[pair.col_]).reversePrimerRaw_, "DNA"));
			primerPairWiseCompsOut << primerNames[pair.row_] << "\t" << "reverse-5-3" << "\t" << primerNames[pair.col_] << "\t" << "reverse-3-5";
			primerPairWiseCompsOut << "\t" << alnObject.comp_.distances_.eventBasedIdentity_ << "\t" <<  alnObject.comp_.distances_.query_.identities_ << "\t" << njh::boolToStr(allowablePairwisePrimerComps.passErrorProfile(alnObject.comp_)) << std::endl;
		}
	}


	//5-3 and 3-5 first base counts

	DNABaseCounter right53_counter;
	DNABaseCounter left53_counter;
	DNABaseCounter combined53_counter;
	DNABaseCounter right35_counter;
	DNABaseCounter left35_counter;
	DNABaseCounter combined35_counter;

	for(const auto & primer : primers.pDeterminator_->primers_){
		left53_counter.increase(primer.second.forwardPrimerRaw_.front());
		left35_counter.increase(primer.second.forwardPrimerRaw_.back());

		right53_counter.increase(primer.second.reversePrimerRaw_.front());
		right35_counter.increase(primer.second.reversePrimerRaw_.back());

		combined53_counter.increase(primer.second.forwardPrimerRaw_.front());
		combined53_counter.increase(primer.second.reversePrimerRaw_.front());
		combined35_counter.increase(primer.second.forwardPrimerRaw_.back());
		combined35_counter.increase(primer.second.reversePrimerRaw_.back());
	}


	OutputStream firstLastPrimerBaseCounts(njh::files::make_path(setUp.pars_.directoryName_, "firstAndLastPrimerBaseCounts.txt"));
	firstLastPrimerBaseCounts << "label\tbase\tcount\tfreq" << std::endl;
	auto writeCounts = [&firstLastPrimerBaseCounts](DNABaseCounter& counter, const std::string & label){
		counter.setFractions();
		std::vector<char> bases{'A', 'C', 'G','T'};

		for(char base : bases){
			firstLastPrimerBaseCounts << label
					<< "\t" << base
					<< "\t" << counter.bases_[base]
					<< "\t" << counter.baseFracs_[base]
					<< std::endl;
		}
	};
	writeCounts(right53_counter, "5-3_right");
	writeCounts(left53_counter, "5-3_left");
	writeCounts(combined53_counter, "5-3_combined");

	writeCounts(right35_counter, "3-5_right");
	writeCounts(left35_counter, "3-5_left");
	writeCounts(combined35_counter, "3-5_combined");


	return 0;
}


int programWrapperRunner::parsePrimer3OutputToPossibleMipArms(const njh::progutils::CmdArgs & inputCommands){
//	bfs::path extInput = "";
//	bfs::path ligInput = "";
//	uint32_t minLen = 100;
//	uint32_t maxLen = 400;
//	OutOptions outOpts(bfs::path(""));
//	seqSetUp setUp(inputCommands);
//	setUp.processVerbose();
//	setUp.processDebug();
//	setUp.processWritingOptions(outOpts);
//	setUp.setOption(minLen, "--minLen", "Minimum length");
//	setUp.setOption(maxLen, "--maxLen", "Maximum length");
//	setUp.setOption(ligInput, "--ligInput", "input file for ligation arm", true);
//	setUp.setOption(extInput, "--extInput", "input file for extension arm", true);
//	setUp.finishSetUp(std::cout);
//
//
//	OutputStream out(outOpts);
//
//	auto ligResults = Primer3Results::parsePrimer3OutputResults(ligInput);
//	auto ligRegions = ligResults.front()->genPrimersRegions();
//
//	auto extResults = Primer3Results::parsePrimer3OutputResults(extInput);
//	auto extRegions = extResults.front()->genPrimersRegions();
//
//	struct PrimerPair {
//		PrimerPair(const std::shared_ptr<Primer3Results::Primer> & extPrimer,
//				       const std::shared_ptr<Primer3Results::Primer> & ligPrimer) :
//				extPrimer_(extPrimer), ligPrimer_(ligPrimer) {
//		}
//
//		std::shared_ptr<Primer3Results::Primer> extPrimer_;
//		std::shared_ptr<Primer3Results::Primer> ligPrimer_;
//
//		GenomicRegion genRegion(const std::string & seqId) const{
//			MetaDataInName meta;
//			meta.addMeta("ext", extPrimer_->name_);
//			meta.addMeta("lig", ligPrimer_->name_);
//			meta.addMeta("seqId", seqId);
//			auto extRegion = extPrimer_->genRegion(seqId);
//			auto ligRegion = ligPrimer_->genRegion(seqId);
//			return GenomicRegion(meta.createMetaName(), seqId,
//					std::min(extRegion.start_, ligRegion.start_),
//					std::max(extRegion.end_, ligRegion.end_), extRegion.reverseSrand_);
//		}
//	};
//	std::vector<PrimerPair> pairs;
//	for(const auto & ext : extRegions){
//		for(const auto & lig : ligRegions){
//			if(!ext.second.reverseSrand_){
//				if(lig.second.reverseSrand_){
//					if(lig.second.start_ > ext.second.end_ &&
//							lig.second.end_ - ext.second.start_ > minLen &&
//							lig.second.end_ - ext.second.start_ < maxLen){
//						pairs.emplace_back(PrimerPair{extResults.front()->primers_.at(ext.first),
//																			   ligResults.front()->primers_.at(lig.first)});
//					}
//				}
//			}else{
//				if(!lig.second.reverseSrand_){
//					if(ext.second.start_ > lig.second.end_ &&
//							ext.second.end_ - lig.second.start_ > minLen &&
//							ext.second.end_ - lig.second.start_ < maxLen){
//						pairs.emplace_back(PrimerPair{extResults.front()->primers_.at(ext.first),
//																			   ligResults.front()->primers_.at(lig.first)});
//					}
//				}
//			}
//		}
//	}
//	for(const auto & pair : pairs){
//		auto region = pair.genRegion(extResults.front()->sequence_id_);
//		out << region.genBedRecordCore().toDelimStr()
//				<< "\t" << (region.reverseSrand_ ? pair.extPrimer_->originalSeq_ : pair.extPrimer_->fowardOrientationSeq_ )
//				<< "\t" << (region.reverseSrand_ ? pair.ligPrimer_->fowardOrientationSeq_ : pair.ligPrimer_->originalSeq_)
//				<< std::endl;
//	}
	return 0;
}

int programWrapperRunner::parsePrimer3OutputToJson(const njh::progutils::CmdArgs & inputCommands){
	bfs::path input = "STDIN";
	bool ignoreUnexpected = false;
	OutOptions outOpts(bfs::path(""));
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processWritingOptions(outOpts);
	setUp.setOption(input, "--input", "input file or STDIN for reading for standard in");
	setUp.setOption(ignoreUnexpected, "--ignoreUnexpected", "ignore Unexpected");

	setUp.finishSetUp(std::cout);


	OutputStream out(outOpts);

	auto results = Primer3Results::parsePrimer3OutputResults(input, ignoreUnexpected);
	out << njh::json::toJson(results) << std::endl;

	return 0;
}

int programWrapperRunner::parsePrimer3OutputToTable(const njh::progutils::CmdArgs & inputCommands){
	bfs::path input = "STDIN";
	bool ignoreUnexpected = false;
	OutOptions outOpts(bfs::path(""));
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processWritingOptions(outOpts);
	setUp.setOption(input, "--input", "input file or STDIN for reading for standard in");
	setUp.setOption(ignoreUnexpected, "--ignoreUnexpected", "ignore Unexpected");

	setUp.finishSetUp(std::cout);


	OutputStream out(outOpts);

	auto results = Primer3Results::parsePrimer3OutputResults(input, ignoreUnexpected);

	out << njh::conToStr(VecStr{"seqID", "primerPairName", "compl_any_th", "compl_end_th", "penalty", "product_size",
	"left_name", "left_seq", "left_start", "left_end", "left_size", "left_end_stability", "left_gc_percent", "left_hairpin_th", "left_penalty", "left_self_any_th", "left_self_end_th", "left_tm", "left_tm_hairpin_diff", "left_problems",
	"right_name", "right_seq", "right_start", "right_end", "right_size", "right_end_stability", "right_gc_percent", "right_hairpin_th", "right_penalty", "right_self_any_th", "right_self_end_th", "right_tm", "right_tm_hairpin_diff", "right_problems"
	}, "\t") << std::endl;
	for(const auto & res : results){
		for(const auto & primerPair : res->primerPairs_){
			out << njh::conToStr(toVecStr(
				res->sequence_id_,
				primerPair->name_,
				primerPair->compl_any_th_,
				primerPair->compl_end_th_,
				primerPair->penalty_,
				primerPair->product_size_,

				primerPair->left_.name_,
				primerPair->left_.seq_,
				primerPair->left_.forwardOrientationPos_.start_,
				primerPair->left_.forwardOrientationPos_.start_ + primerPair->left_.forwardOrientationPos_.size_,
				primerPair->left_.forwardOrientationPos_.size_,
				primerPair->left_.end_stability_,
				primerPair->left_.gc_percent_,
				primerPair->left_.hairpin_th_,
				primerPair->left_.penalty_,
				primerPair->left_.self_any_th_,
				primerPair->left_.self_end_th_,
				primerPair->left_.tm_,
				primerPair->left_.tm_ - primerPair->left_.hairpin_th_,
				njh::conToStr(primerPair->left_.problems_, ";"),

				primerPair->right_.name_,
				primerPair->right_.seq_,
				primerPair->right_.forwardOrientationPos_.start_,
				primerPair->right_.forwardOrientationPos_.start_ + primerPair->right_.forwardOrientationPos_.size_,
				primerPair->right_.forwardOrientationPos_.size_,
				primerPair->right_.end_stability_,
				primerPair->right_.gc_percent_,
				primerPair->right_.hairpin_th_,
				primerPair->right_.penalty_,
				primerPair->right_.self_any_th_,
				primerPair->right_.self_end_th_,
				primerPair->right_.tm_,
				primerPair->right_.tm_ - primerPair->right_.hairpin_th_,
				njh::conToStr(primerPair->right_.problems_, ";")
			), "\t") << std::endl;
		}
	}
	return 0;
}




int programWrapperRunner::parsePrimer3OutputToBed(const njh::progutils::CmdArgs & inputCommands){
//	bfs::path input = "STDIN";
//	OutOptions outOpts(bfs::path(""));
//	outOpts.outExtention_ = ".bed";
//	seqSetUp setUp(inputCommands);
//	setUp.processVerbose();
//	setUp.processDebug();
//	setUp.processWritingOptions(outOpts);
//	setUp.setOption(input, "--input", "input file or STDIN for reading for standard in");
//	setUp.finishSetUp(std::cout);
//
//
//	OutputStream out(outOpts);
//	//addMetaBySampleName
//	auto results = Primer3Results::parsePrimer3OutputResults(input);
//	for(const auto & result : results){
//		for(const auto & primer : result->primers_){
//			out << result->sequence_id_
//					<< "\t" << primer.second->forwardOrientationPos_.start_
//					<< "\t" << primer.second->forwardOrientationPos_.start_ + primer.second->forwardOrientationPos_.size_
//					<< "\t" << primer.second->name_
//					<< "\t" << primer.second->forwardOrientationPos_.size_
//					<< "\t" << (primer.second->right_ ? '-' : '+')
//					<< std::endl;
//		}
//	}


	return 0;
}

int programWrapperRunner::findNonUniquePrimerArms(const njh::progutils::CmdArgs & inputCommands){
	bfs::path inDir = "";
	bfs::path outDir = "";
	std::string pattern = "";
	bool overWriteDir = false;
	uint32_t strainCutOff = 2;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(strainCutOff,   "--strainCutOff",   "Strain Cut Off");
	setUp.setOption(overWriteDir,   "--overWriteDir",   "Over Write Dir");
	setUp.setOption(inDir,   "--inDir",   "Input directory",  true);
	setUp.setOption(outDir,  "--outDir",  "Output directory", true);
	setUp.setOption(pattern, "--pattern", "Pattern",          true);
	setUp.finishSetUp(std::cout);

	njh::files::makeDir(njh::files::MkdirPar(outDir, overWriteDir));

	std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>> counts;
	auto allFiles = njh::files::listAllFiles(inDir, false, std::vector<std::regex>{std::regex{pattern}});

	for (const auto & f : allFiles) {
		if (f.second) {
			continue;
		}
		if (setUp.pars_.verbose_) {
			std::cout << f.first << std::endl;
		}
		BioDataFileIO<Bed6RecordCore> bedReader {
				IoOptions { InOptions { f.first } } };
		bedReader.openIn();
		Bed6RecordCore bRecord;
		while (bedReader.readNextRecord(bRecord)) {
			if(bRecord.extraFields_.size() < 2){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__<< ", error, extra fields should have at least two fields" << "\n";
				throw std::runtime_error{ss.str()};
			}
			std::string strainName = bRecord.chrom_.substr(0, bRecord.chrom_.find("."));
			++counts[bRecord.extraFields_[0] + "-" + bRecord.extraFields_[1]][strainName];
		}
	}
	if (setUp.pars_.debug_) {
		for(const auto & count : counts){
			if(count.second.size() > 1){
				std::cout << count.first << "\t" << count.second.size() << std::endl;
			}
		}
	}

//generatingPrime3TemplatesBasedOnMALN
	for (const auto & f : allFiles) {
		if (f.second) {
			continue;
		}
		if (setUp.pars_.verbose_) {
			std::cout << f.first << std::endl;
		}
		BioDataFileIO<Bed6RecordCore> bedReader {
				IoOptions { InOptions { f.first }, OutOptions(njh::files::make_path(outDir, bfs::basename(f.first))) } };
		bedReader.openIn();
		bedReader.openOut();
		Bed6RecordCore bRecord;
		while (bedReader.readNextRecord(bRecord)) {
			if(bRecord.extraFields_.size() < 2){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__<< ", error, extra fields should have at least two fields" << "\n";
				throw std::runtime_error{ss.str()};
			}
			std::string strainName = bRecord.chrom_.substr(0, bRecord.chrom_.find("."));
			if(counts[bRecord.extraFields_[0] + "-" + bRecord.extraFields_[1]][strainName] >= strainCutOff){
				bRecord.extraFields_.emplace_back(estd::to_string(counts[bRecord.extraFields_[0] + "-" + bRecord.extraFields_[1]].size()));
				bRecord.extraFields_.emplace_back(estd::to_string(counts[bRecord.extraFields_[0] + "-" + bRecord.extraFields_[1]][strainName]));
				bedReader.write(bRecord, [](const Bed6RecordCore & bRecord, std::ostream & out){
					out << bRecord.toDelimStrWithExtra() << std::endl;
				});
			}
		}
	}

	return 0;
}




} // namespace njhseq
