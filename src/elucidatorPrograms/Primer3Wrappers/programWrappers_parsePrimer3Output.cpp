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

#include "elucidatorPrograms/programWrappers/programWrappers.hpp"
#include "elucidator/objects/BioDataObject.h"
#include "elucidator/BamToolsUtils.h"

#include <SeekDeep/objects/TarAmpSetupUtils/PrimersAndMids.hpp>

#include <njhseq/concurrency/PairwisePairFactory.hpp>

namespace njhseq {




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
	std::string task = "generic";
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processWritingOptions(outOpts);
	setUp.setOption(input, "--input", "input file or STDIN for reading for standard in");
	setUp.setOption(ignoreUnexpected, "--ignoreUnexpected", "ignore Unexpected");
	setUp.setOption(task, "--task", "task");

	setUp.finishSetUp(std::cout);


	OutputStream out(outOpts);
	if(njh::in(task, VecStr{"generic", "prick_sequencing_primers"})){
		auto results = Primer3Runner::Primer3ResultsGeneric::parsePrimer3OutputResults(input, ignoreUnexpected);
		out << njh::json::toJson(results) << std::endl;
	}else if("pick_primer_list" == task){
		auto results = Primer3Runner::Primer3ResultsList::parsePrimer3OutputResults(input, ignoreUnexpected);
		out << njh::json::toJson(results) << std::endl;
	}else{
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "unhandled task, " << task<< "\n";
		throw std::runtime_error{ss.str()};
	}
	return 0;
}


int programWrapperRunner::parsePrimer3OutputToTable(const njh::progutils::CmdArgs & inputCommands){
	bfs::path input = "STDIN";
	bool doNotIgnoreUnexpected = false;

	OutOptions outOpts(bfs::path(""));
	std::string task = "generic";
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processWritingOptions(outOpts);
	setUp.setOption(input, "--input", "input file or STDIN for reading for standard in");
	setUp.setOption(doNotIgnoreUnexpected, "--doNotIgnoreUnexpected", "do not ignore Unexpected");
	bool ignoreUnexpected = !doNotIgnoreUnexpected;

	setUp.setOption(task, "--task", "task");

	setUp.finishSetUp(std::cout);


	OutputStream out(outOpts);
	if(njh::in(task, VecStr{"generic", "prick_sequencing_primers"})){
		auto results = Primer3Runner::Primer3ResultsGeneric::parsePrimer3OutputResults(input, ignoreUnexpected);

		out << njh::conToStr(VecStr{"seqID", "primerPairName", "compl_any_th", "compl_end_th", "penalty", "product_size",
		"left_name", "left_seq", "left_start", "left_end", "left_size", "left_end_stability", "left_gc_percent", "left_hairpin_th", "left_penalty",  "left_self_any_th", "left_self_end_th", "left_tm", "left_tm_hairpin_diff", "left_problems",
		"right_name", "right_seq", "right_start", "right_end", "right_size", "right_end_stability", "right_gc_percent", "right_hairpin_th", "right_penalty",  "right_self_any_th", "right_self_end_th", "right_tm", "right_tm_hairpin_diff", "right_problems"
		}, "\t") << std::endl;
		for(const auto & res : results){
			for(const auto & primerPair : res->primerPairs_){
				out << njh::conToStr(toVecStr(
					res->sequence_id_,
					primerPair->name_,
					primerPair->compl_any_th_,
					primerPair->compl_end_th_,
					primerPair->penalty_,
					primerPair->penalty_ ,
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
	}else if("pick_primer_list" == task){
		auto results = Primer3Runner::Primer3ResultsList::parsePrimer3OutputResults(input, ignoreUnexpected);

		out << njh::conToStr(VecStr{"seqID",
		"primer_name", "primer_side", "primer_seq", "primer_start", "primer_end", "primer_size", "primer_end_stability", "primer_gc_percent", "primer_hairpin_th", "primer_penalty",  "primer_self_any_th", "primer_self_end_th", "primer_tm", "primer_tm_hairpin_diff", "primer_problems"
		}, "\t") << std::endl;


		for(const auto & res : results){

			for(const auto & leftPrimer : res->leftPrimers_){
				out << njh::conToStr(toVecStr(
					res->sequence_id_,
					leftPrimer->name_,
					"left",
					leftPrimer->seq_,
					leftPrimer->forwardOrientationPos_.start_,
					leftPrimer->forwardOrientationPos_.start_ + leftPrimer->forwardOrientationPos_.size_,
					leftPrimer->forwardOrientationPos_.size_,
					leftPrimer->end_stability_,
					leftPrimer->gc_percent_,
					leftPrimer->hairpin_th_,
					leftPrimer->penalty_,
					leftPrimer->self_any_th_,
					leftPrimer->self_end_th_,
					leftPrimer->tm_,
					leftPrimer->tm_ - leftPrimer->hairpin_th_,
					njh::conToStr(leftPrimer->problems_, ";")
				), "\t") << std::endl;
			}

			for(const auto & rightPair : res->rightPrimers_){
				out << njh::conToStr(toVecStr(
					res->sequence_id_,
					rightPair->name_,
					"right",
					rightPair->seq_,
					rightPair->forwardOrientationPos_.start_,
					rightPair->forwardOrientationPos_.start_ + rightPair->forwardOrientationPos_.size_,
					rightPair->forwardOrientationPos_.size_,
					rightPair->end_stability_,
					rightPair->gc_percent_,
					rightPair->hairpin_th_,
					rightPair->penalty_,
					rightPair->self_any_th_,
					rightPair->self_end_th_,
					rightPair->tm_,
					rightPair->tm_ - rightPair->hairpin_th_,
					njh::conToStr(rightPair->problems_, ";")
				), "\t") << std::endl;
			}
		}
	}else {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "unhandled task, " << task<< "\n";
		throw std::runtime_error{ss.str()};
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
