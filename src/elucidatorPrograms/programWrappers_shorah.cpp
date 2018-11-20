/*
 * programWrappers_shorah.cpp
 *
 *  Created on: Oct 26, 2017
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
#include <njhseq/GenomeUtils.h>


namespace njhseq {

void processReadObjectForShorahOutput(readObject & seq) {
	//process names
	auto nameToks = tokenizeString(seq.seqBase_.name_, "|");
	if (2 != nameToks.size()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__
				<< ", error, nameToks should be two tokens separated by a |, found "
				<< nameToks.size() << " for " << seq.seqBase_.name_ << "\n";
		throw std::runtime_error { ss.str() };
	}
	auto secondToks = tokenizeString(nameToks[1], "whitespace");
	if (2 != secondToks.size()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__
				<< ", error, nameToks should be two tokens separated by a whitespace, found "
				<< secondToks.size() << " for " << nameToks[1] << "\n";
		throw std::runtime_error { ss.str() };
	}
	std::unordered_map<std::string, double> pars;
	for (const auto & tok : secondToks) {
		auto tokOfToks = tokenizeString(tok, "=");
		if (2 != secondToks.size()) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__
					<< ", error, nameToks should be two tokens separated by a =, found "
					<< tokOfToks.size() << " for " << tok << "\n";
			throw std::runtime_error { ss.str() };
		}
		pars[tokOfToks[0]] = njh::StrToNumConverter::stoToNum<double>(tokOfToks[1]);
	}

	if (!njh::in(std::string("posterior"), pars)
			|| !njh::in(std::string("ave_reads"), pars)) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__
				<< ", need to have both posterior and ave_reads in name for "
				<< seq.seqBase_.name_ << "\n";
		throw std::runtime_error { ss.str() };
	}
	seq.seqBase_.name_ = nameToks[0];
	seq.seqBase_.cnt_ = std::round(pars.at("ave_reads"));
	seq.meta_ = MetaDataInName();
	for (const auto & par : pars) {
		seq.meta_.addMeta(par.first, par.second);
	}
	seq.resetMetaInName();
	seq.updateName();
	//apparently they add gaps for some reason
	seq.seqBase_.removeGaps();
}

readObject processReadObjectForShorahOutputRet(readObject seq) {
	processReadObjectForShorahOutput(seq);
	return seq;
}


std::vector<readObject> processReadObjectVecForShorahOutput(const std::vector<readObject> & seqs, double posteriorCutOff) {
	std::vector<readObject> ret;
	for(const auto & inseq : seqs){
		auto seq = processReadObjectForShorahOutputRet(inseq);
		if (seq.seqBase_.cnt_ > 1
				&& seq.meta_.getMeta<double>("posterior") > posteriorCutOff) {
			ret.emplace_back(seq);
		}
	}
	readVec::allSetFractionByTotalCount(ret);
	readVecSorter::sort(ret);
	return ret;
}



int programWrapperRunner::convertShorahSupportSeqs(
		const njh::progutils::CmdArgs & inputCommands) {
	double posteriorCutOff = 0;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processDefaultReader(true);
	setUp.setOption(posteriorCutOff, "--posteriorCutOff",
			"Get rid of reads with a posterior less than this cut off");
	setUp.finishSetUp(std::cout);

	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	reader.openOut();
	readObject seq;
	std::vector<readObject> input = reader.in_.readAllReads<readObject>();
	auto output = processReadObjectVecForShorahOutput(input, posteriorCutOff);
	reader.write(output);


	return 0;
}


int programWrapperRunner::runShorahAmplian(
		const njh::progutils::CmdArgs & inputCommands) {

	bfs::path additionalOutLocationFnp = "";
	double posteriorCutOff = 0;
	bfs::path genomeFnp = "";
	bfs::path bedFnp = "";
	seqSetUp setUp(inputCommands);
	setUp.pars_.ioOptions_.lowerCaseBases_ = "remove";
	setUp.pars_.ioOptions_.out_.outFilename_ = "output.fasta";
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(genomeFnp, "--genome", "Genome to align to", true);
	setUp.setOption(additionalOutLocationFnp, "--additionalOut", "File to tell where else to write seqs");
	setUp.setOption(bedFnp, "--bed",
			"bed file to analyze a given region, otherwise will be automatically determined from region of highest coverage");
	setUp.processDefaultReader(true);
	setUp.processDirectoryOutputName(true);
	setUp.setOption(posteriorCutOff, "--posteriorCutOff", "Get rid of reads with a posterior less than this cut off");
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	if("" != bedFnp){
		if(!bfs::exists(bedFnp)){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error: " << bedFnp << " doesn't exist" << "\n";
			throw std::runtime_error{ss.str()};
		}
	}
	std::stringstream cmdsUsed;

	auto inputSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(setUp.pars_.directoryName_, "inputSeq.fasta"));
	SeqInput inputReader(setUp.pars_.ioOptions_);
	SeqOutput inputWriter(inputSeqOpts);
	inputReader.openIn();
	inputWriter.openOut();
	seqInfo inputSeq;
	while(inputReader.readNextRead(inputSeq)){
		readVec::handelLowerCaseBases(inputSeq, setUp.pars_.ioOptions_.lowerCaseBases_);
		inputWriter.write(inputSeq);
	}
	inputWriter.closeOut();
	auto inputSeqProcessedInputOpts = SeqIOOptions::genFastaIn(inputSeqOpts.out_.outName());

	//align seqs
	auto initialAlignOpts = inputSeqProcessedInputOpts;
	initialAlignOpts.out_.outFilename_ = njh::files::make_path(setUp.pars_.directoryName_, "initial_aligned.sorted.bam");
	initialAlignOpts.out_.outExtention_ = ".sorted.bam";
	BioCmdsUtils bRunner(setUp.pars_.verbose_);
	bRunner.RunFaToTwoBit(genomeFnp);
	auto outputRunOpts = bRunner.bowtie2Align(initialAlignOpts, genomeFnp);
	cmdsUsed << outputRunOpts.cmd_ << std::endl;
	GenomicRegion reg;
	//determine region
	if("" != bedFnp){
		auto bed = getBed3s(bedFnp);
		reg = GenomicRegion(*(bed.front()));
	}else{
		BamTools::BamReader bReader;
		bReader.Open(initialAlignOpts.out_.outName().string());
		checkBamOpenThrow(bReader, initialAlignOpts.out_.outName().string());
		BamTools::BamAlignment bAln;
		std::vector<BamTools::BamAlignment> bAlns;
		auto refIds = bReader.GetReferenceData();
		while (bReader.GetNextAlignment(bAln)) {
			if (bAln.IsMapped()) {
				bAlns.emplace_back(bAln);
			}
		}
		GenomicRegionCounter gCounter;
		gCounter.increaseCounts(bAlns, refIds);
		reg = gCounter.getRegionsLargestOnTop().front();
	}

	OutOptions regionFileOpts(njh::files::make_path(setUp.pars_.directoryName_, "regionUsed.bed"));
	OutputStream regionFileOut(regionFileOpts);
	regionFileOut << reg.genBedRecordCore().toDelimStr() << std::endl;

	//re-align seqs because shorah is the worst
	auto finalAlignOpts = inputSeqProcessedInputOpts;
	finalAlignOpts.out_.outFilename_ = njh::files::make_path(setUp.pars_.directoryName_,
			"aligned.sorted.bam");
	finalAlignOpts.out_.outExtention_ = ".sorted.bam";

	//extract out just the region covered
	reg.setUidWtihCoords();
	auto twobitFnp = genomeFnp;
	twobitFnp.replace_extension(".2bit");
	TwoBit::TwoBitFile tReader(twobitFnp);
	auto regionSeq = reg.extractSeq(tReader);
	auto regionSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(setUp.pars_.directoryName_, reg.uid_));
	SeqOutput::write(std::vector<seqInfo>{regionSeq},regionSeqOpts);

	auto finalAlignOutputRunOpts = bRunner.bowtie2Align(finalAlignOpts, regionSeqOpts.out_.outName());
	cmdsUsed << finalAlignOutputRunOpts.cmd_ << std::endl;


	std::stringstream shorahCmd;
//	shorahCmd << "cd " << setUp.pars_.directoryName_ << " && amplian.py -b "
//			<< bfs::absolute(finalAlignOpts.out_.outName()) << " -f "
//			<< bfs::absolute(genomeFnp) << " -r " << reg.chrom_ << ":"
//			<< reg.start_ + 1 << "-" << reg.end_;
	shorahCmd << "cd " << setUp.pars_.directoryName_ << " && amplian.py -b "
			<< bfs::absolute(finalAlignOpts.out_.outName()) << " -f "
			<< bfs::absolute(regionSeqOpts.out_.outName());
	auto shorahOutput = njh::sys::run({shorahCmd.str()});
	bRunner.checkRunOutThrow(shorahOutput, __PRETTY_FUNCTION__);
	cmdsUsed << shorahOutput.cmd_;
	std::stringstream shorahSupportReadsFnp;
	//w-Pf3d7-Pf3D7_03_v3-221419-221693-1-274.reads-support.fas
//	shorahSupportReadsFnp <<  "w-" << reg.chrom_ << "-"
//			<< reg.start_ + 1 << "-" << reg.end_  << ".reads-support.fas";
	shorahSupportReadsFnp << "w-" << regionSeq.name_ << "-1" << "-" << reg.getLen() << ".reads-support.fas";

	SeqIOOptions shorahOutputOpts = SeqIOOptions::genFastaIn(njh::files::make_path(setUp.pars_.directoryName_, shorahSupportReadsFnp.str()));
	uint64_t maxLen = 0;
	auto shorahReads = SeqInput::getReferenceSeq(shorahOutputOpts, maxLen);
	auto processedShorahReads = processReadObjectVecForShorahOutput(shorahReads, posteriorCutOff);


	OutOptions infoOpts(
			njh::files::make_path(setUp.pars_.directoryName_, "outputInfo.tab.txt"));
	OutputStream infoOut(infoOpts);
	auto seqOutOpts = SeqIOOptions::genFastaOut(
			njh::files::make_path(setUp.pars_.directoryName_, "output.fasta"));
	SeqOutput seqOutWriter(seqOutOpts);
	seqOutWriter.openOut();
	infoOut << "Name\tcount\tfrac" << std::endl;
	for(const auto & seq : processedShorahReads){
		seqOutWriter.write(seq);
		infoOut << seq.seqBase_.name_
				<< "\t" << seq.seqBase_.cnt_
				<< "\t" << seq.seqBase_.frac_ << std::endl;
	}
	OutOptions cmdsFileOpts(njh::files::make_path(setUp.pars_.directoryName_, "cmdsUsed.txt"));
	OutputStream cmdsFileOut(cmdsFileOpts);
	cmdsFileOut << cmdsUsed.str() << std::endl;

	if ("" != additionalOutLocationFnp) {
		Json::Value metaData;
		auto analysisDirPath = njh::files::bfs::canonical(setUp.pars_.directoryName_);
		metaData["analysisDirPath"] = njh::json::toJson(analysisDirPath.string());
		auto fullPathToInput = njh::files::bfs::canonical(setUp.pars_.ioOptions_.firstName_);
		metaData["inputFile"] = njh::json::toJson(fullPathToInput);
		metaData["extractionDir"] = njh::json::toJson(fullPathToInput.parent_path());
		std::string additionalOutDir = findAdditonalOutLocation(
				additionalOutLocationFnp.string(), setUp.pars_.ioOptions_.firstName_.string());
		if (additionalOutDir == "") {
			std::cerr << njh::bashCT::red << njh::bashCT::bold;
			std::cerr << "No additional out directory found for: "
					<< setUp.pars_.ioOptions_.firstName_ << std::endl;
			std::cerr << njh::bashCT::reset;
		} else {
			setUp.pars_.ioOptions_.out_.outExtention_ = ".fasta";
			SeqOutput::write(processedShorahReads, SeqIOOptions(additionalOutDir + setUp.pars_.ioOptions_.out_.outFilename_.string(),
					SeqIOOptions::outFormats::FASTA,setUp.pars_.ioOptions_.out_));
			std::ofstream metaDataFile;
			openTextFile(metaDataFile, additionalOutDir + "/" + "metaData", ".json",
					setUp.pars_.ioOptions_.out_);
			metaDataFile << metaData;
		}
	}
	return 0;
}



} //namespace njhseq

