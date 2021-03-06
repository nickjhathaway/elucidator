/*
 * seqSearching.cpp
 *
 *  Created on: Nov 12, 2018
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

#include "seqSearching.hpp"

#include "elucidator/objects/BioDataObject.h"
#include "elucidator/BamToolsUtils/BamUtilities.hpp"
#include "elucidator/BioRecordsUtils/BedUtility.hpp"

namespace njhseq {


seqSearchingRunner::seqSearchingRunner()
    : njh::progutils::ProgramRunner(
          {

					 addFunc("chopAndMap", chopAndMap, false),
					 addFunc("chopAndMapAndRefine", chopAndMapAndRefine, false),
					 addFunc("findMotifLocations", findMotifLocations, false),
					 addFunc("findTandemMotifLocations", findTandemMotifLocations, false),
					 addFunc("chopAndMapAndRefineInvidual", chopAndMapAndRefineInvidual, false),
           },//
          "seqSearching") {}



int seqSearchingRunner::findTandemMotifLocations(const njh::progutils::CmdArgs & inputCommands){
	std::string motifstr = "";
	//bfs::path genomeFnp = "";

	uint32_t allowableErrors = 0;
	uint32_t maxAllowableErrors = 0;
	OutOptions outOpts(bfs::path(""));
	outOpts.outExtention_ = ".bed";
	uint32_t repeatCutOff = 3;
	bool noReverse = false;
	seqSetUp setUp(inputCommands);
	setUp.setOption(noReverse, "--noReverse", "Don't look in reverse complement");

	setUp.pars_.ioOptions_.includeWhiteSpaceInName_ = false;
	//setUp.setOption(genomeFnp, "--genomeFnp", "The genome file to look for motifs in", true);
	setUp.processReadInNames({"--fasta", "--fastagz"}, true);
	setUp.setOption(motifstr, "--motif", "The motif to look for", true);
	setUp.setOption(repeatCutOff, "--repeatCutOff", "The minimum number of times the motif repeats to report it");
	setUp.setOption(allowableErrors, "--allowableErrors", "allowable errors in a motif element");
	maxAllowableErrors = allowableErrors;
	setUp.setOption(maxAllowableErrors, "--maxAllowableErrors", "max allowable errors in the whole tandem motif sequence");

	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	seqInfo seq;
	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	motif mot(motifstr);
	OutputStream out(outOpts);
	while(reader.readNextRead(seq)){
		auto locs = mot.findPositionsFull(seq.seq_, allowableErrors);
		njh::sort(locs);
		if(!locs.empty()){
			uint32_t length = 1;
			size_t start = locs.front();
			for(const auto pos : iter::range<uint32_t>(1, locs.size())){
				if(locs[pos] == locs[pos - 1] + mot.size() ){
					++length;
				} else {
					if(length >= repeatCutOff){
						uint32_t numOfErrors = 0;
						for(const auto & seqPos : iter::range(start, start + mot.size() * length, mot.size())){
							numOfErrors += mot.size() - mot.scoreMotif(seq.seq_.begin() + seqPos, seq.seq_.begin() + seqPos + mot.size());
						}
						if(numOfErrors <= maxAllowableErrors){
							out << seq.name_
									<< "\t" << start
									<< "\t" << start + mot.size() * length
									<< "\t" << motifstr << "_x" << length
									<< "\t" << length
									<< "\t" << "+" << '\n';;
						}
					}
					length = 1;
					start = locs[pos];
				}
			}
			if(length >= repeatCutOff){
				uint32_t numOfErrors = 0;
				for(const auto & seqPos : iter::range(start, start + mot.size() * length, mot.size())){
					numOfErrors += mot.size() - mot.scoreMotif(seq.seq_.begin() + seqPos, seq.seq_.begin() + seqPos + mot.size());
				}
				if(numOfErrors <= maxAllowableErrors){
					out << seq.name_
							<< "\t" << start
							<< "\t" << start + mot.size() * length
							<< "\t" << motifstr << "_x" << length
							<< "\t" << length
							<< "\t" << "+" << '\n';;
				}
			}
		}
		if(!noReverse){
			seq.reverseComplementRead(false, true);
			auto revLocs = mot.findPositionsFull(seq.seq_, allowableErrors);
			if(!revLocs.empty()){
				njh::sort(revLocs);
				uint32_t length = 1;
				size_t start = revLocs.front();
				for(const auto  pos : iter::range<uint32_t>(1, revLocs.size())){
					if(revLocs[pos] == revLocs[pos - 1] + mot.size()){
						++length;
					}else{
						if(length >= repeatCutOff){
							uint32_t numOfErrors = 0;
							for(const auto & seqPos : iter::range(start, start + mot.size() * length, mot.size())){
								numOfErrors += mot.size() - mot.scoreMotif(seq.seq_.begin() + seqPos, seq.seq_.begin() + seqPos + mot.size());
							}
							if(numOfErrors <= maxAllowableErrors){
								out << seq.name_
										<< "\t" << len(seq) - (start + mot.size() * length )
										<< "\t" << len(seq) - start
										<< "\t" << motifstr << "_x" << length
										<< "\t" << length
										<< "\t" << "-" << '\n';;
							}
						}
						length = 1;
						start = revLocs[pos];
					}
				}
				if(length >= repeatCutOff){
					uint32_t numOfErrors = 0;
					for(const auto & seqPos : iter::range(start, start + mot.size() * length, mot.size())){
						numOfErrors += mot.size() - mot.scoreMotif(seq.seq_.begin() + seqPos, seq.seq_.begin() + seqPos + mot.size());
					}
					if(numOfErrors <= maxAllowableErrors){
						out << seq.name_
								<< "\t" << len(seq) - (start + mot.size() * length)
								<< "\t" << len(seq) - start
								<< "\t" << motifstr << "_x" << length
								<< "\t" << length
								<< "\t" << "-" << '\n';;
					}
				}
			}
		}
	}

	return 0;
}

int seqSearchingRunner::findMotifLocations(const njh::progutils::CmdArgs & inputCommands){
	std::string motifstr = "";
	seqInfo motifObj;
	//bfs::path genomeFnp = "";
	uint32_t allowableErrors = 0;
	OutOptions outOpts(bfs::path(""));
	outOpts.outExtention_ = ".bed";
	bool noReverse = false;
	seqSetUp setUp(inputCommands);
	setUp.pars_.ioOptions_.includeWhiteSpaceInName_ = false;
	//setUp.setOption(genomeFnp, "--genomeFnp", "The genome file to look for motifs in", true);
	setUp.processReadInNames({"--fasta", "--fastagz"}, true);
	//setUp.setOption(motifstr, "--motif", "The motif to look for", true);
	setUp.processSeq(motifObj, "--motif", "The motif to look for", true);
	setUp.setOption(allowableErrors, "--allowableErrors", "allowable errors in motif");
	setUp.setOption(noReverse, "--noReverse", "Don't look in reverse complement");
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	motifstr = motifObj.seq_;
	seqInfo seq;
	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	motif mot(motifstr);
	OutputStream out(outOpts);
	while(reader.readNextRead(seq)){
		auto locs = mot.findPositionsFull(seq.seq_, allowableErrors);
		for(const auto & loc : locs){
			out << seq.name_
					<< "\t" << loc
					<< "\t" << loc + mot.size()
					<< "\t" << seq.seq_.substr(loc, mot.size())
					<< "\t" << mot.scoreMotif(seq.seq_.begin() + loc, seq.seq_.begin() + loc + mot.size())
					<< "\t" << "+" << '\n';;
		}
		if(!noReverse){
			seq.reverseComplementRead(false, true);
			auto revLocs = mot.findPositionsFull(seq.seq_, allowableErrors);
			for(const auto & loc : revLocs){
				out << seq.name_
						<< "\t" << len(seq) - (loc + mot.size())
						<< "\t" << len(seq) - loc
						<< "\t" << seq.seq_.substr(loc, mot.size())
						<< "\t" << mot.scoreMotif(seq.seq_.begin() + loc, seq.seq_.begin() + loc + mot.size())
						<< "\t" << "-" << '\n';
			}
		}
	}
	return 0;
}


struct ChopAndMapPars{
	uint32_t windowLength = 100;
	uint32_t windowStep = 10;
	uint32_t numThreads = 1;
	uint32_t perFragmentCount = 10;
	bfs::path genomeFnp = "";
	SeqIOOptions inOpts;
	bfs::path outputDirectory;
	bool debug = false;
};

void runChopAndMap(const ChopAndMapPars & pars){
	njh::sys::requireExternalProgramThrow("bwa");
	njh::files::checkExistenceThrow(pars.genomeFnp, __PRETTY_FUNCTION__);
	auto fragOutOpts = SeqIOOptions::genFastaOut(njh::files::make_path(pars.outputDirectory, "fragments.fasta"));

	{
		SeqOutput fragWriter(fragOutOpts);
		fragWriter.openOut();
		seqInfo seq;
		SeqInput reader(pars.inOpts);
		reader.openIn();
		while (reader.readNextRead(seq)) {
			if (len(seq) < pars.windowLength ) {
				std::cerr << "Seq: " << seq.name_ << " is too short" << std::endl;
			} else {
				for(auto const pos : iter::range<uint32_t>(0, len(seq) + 1 - pars.windowLength, pars.windowStep)){
					auto fragment = seq.getSubRead(pos, pars.windowLength);
					MetaDataInName fragMeta;
					fragMeta.addMeta("start", pos);
					fragMeta.addMeta("end", pos + pars.windowLength);
					fragment.name_.append(fragMeta.createMetaName());
					for(uint32_t fragCount = 0; fragCount <= pars.perFragmentCount; ++ fragCount){
						fragWriter.write(fragment);
					}
				}
			}
		}
	}

	bfs::path fragmentBamFnp =  njh::files::make_path(pars.outputDirectory, "fragments.sorted.bam");
	std::stringstream bwaMappCmd;
	bwaMappCmd << "bwa mem -M -t " << pars.numThreads << " "
			<< pars.genomeFnp  << " "
			<< fragOutOpts.out_.outName() << " "
			<< " 2> " << njh::files::make_path(pars.outputDirectory,"bwa.log.txt")
			<< " | samtools sort - -o " << fragmentBamFnp
			<< " && samtools index " << fragmentBamFnp ;
	auto bwaRunOutput = njh::sys::run({bwaMappCmd.str()});
	BioCmdsUtils::checkRunOutThrow(bwaRunOutput, __PRETTY_FUNCTION__);

}

int seqSearchingRunner::chopAndMap(const njh::progutils::CmdArgs & inputCommands) {
	ChopAndMapPars chopPars;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	chopPars.debug = setUp.pars_.debug_;
	setUp.setOption(chopPars.numThreads, "--numThreads", "Number of threads to use");
	setUp.setOption(chopPars.genomeFnp, "--genomeFnp", "Genome to map to", true);
	setUp.setOption(chopPars.perFragmentCount, "--perFragmentCount", "perFragmentCount");
	setUp.setOption(chopPars.windowLength, "--windowLength", "windowLength");
	setUp.setOption(chopPars.windowStep, "--windowStep", "windowStep");
	setUp.processReadInNames(VecStr{"--fasta", "--fastq"});
	setUp.processDirectoryOutputName(true);
	chopPars.inOpts = setUp.pars_.ioOptions_;
	chopPars.outputDirectory = setUp.pars_.directoryName_;
	setUp.finishSetUp(std::cout);

	runChopAndMap(chopPars);

	return 0;
}


int seqSearchingRunner::chopAndMapAndRefineInvidual(const njh::progutils::CmdArgs & inputCommands) {
	ChopAndMapPars globalChopPars;

	CoverageFinderPars covPars;
	RegionRefinementPars refinePars;
	uint32_t expandLeft = 0;
	uint32_t expandRight = 0;
	uint32_t minLength = 0;
	uint32_t inputMinLength = 0;
	bfs::path gff = "";

	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	globalChopPars.debug = setUp.pars_.debug_;
	setUp.setOption(gff, "--genomegff", "Genome gff file");

	setUp.setOption(globalChopPars.numThreads, "--numThreads", "Number of threads to use");
	setUp.setOption(globalChopPars.genomeFnp, "--genomeFnp", "Genome to map to", true);
	setUp.setOption(globalChopPars.perFragmentCount, "--perFragmentCount", "perFragmentCount");
	setUp.setOption(globalChopPars.windowLength, "--windowLength", "windowLength");
	inputMinLength = globalChopPars.windowLength;
	setUp.setOption(globalChopPars.windowStep, "--windowStep", "windowStep");
	setUp.setOption(refinePars.reOrient, "--reOrient", "reOrient");
	setUp.setOption(inputMinLength, "--inputMinLength", "inputMinLength");

	setUp.setOption(expandLeft, "--expandLeft", "expandLeft");
	setUp.setOption(expandRight, "--expandRight", "expandRight");

	setUp.setOption(minLength, "--minLength", "minLength");

	setUp.processReadInNames(VecStr{"--fasta", "--fastq"});
	setUp.processDirectoryOutputName(true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);

	auto seqs = SeqInput::getSeqVec<seqInfo>(setUp.pars_.ioOptions_);
	for(const auto & seq : seqs){
		if(len(seq) < inputMinLength){
			continue;
		}
		auto seqDirFnp = njh::files::make_path(setUp.pars_.directoryName_, seq.name_);
		njh::files::makeDir(njh::files::MkdirPar{seqDirFnp});
		auto seqOut = SeqIOOptions::genFastaOut(njh::files::make_path(seqDirFnp, seq.name_));
		SeqOutput::write(std::vector<seqInfo>{seq}, seqOut);
		auto chopParsCurrent = globalChopPars;

		chopParsCurrent.inOpts = SeqIOOptions::genFastaIn(seqOut.out_.outName());
		chopParsCurrent.outputDirectory = seqDirFnp;
		//chop
		runChopAndMap(chopParsCurrent);
		//determine coverage
		covPars.numThreads = chopParsCurrent.numThreads;
		covPars.coverageCutOff = 5;
		covPars.window = chopParsCurrent.windowLength;
		covPars.step = chopParsCurrent.windowStep;
		covPars.bams = njh::files::make_path(chopParsCurrent.outputDirectory, "fragments.sorted.bam").string();
		covPars.outOpts = OutOptions(njh::files::make_path(chopParsCurrent.outputDirectory, "coverage.bed"));
		RunCoverageFinder(covPars);
		//merge regions
		auto beds = getBed3s(covPars.outOpts.outName());
		njh::for_each(beds,
				[&expandLeft,&expandRight](const std::shared_ptr<Bed3RecordCore> & region) {
					if(0 != expandLeft) {
						BedUtility::extendLeft(*region, expandLeft);
					}
					if(0 != expandRight) {
						BedUtility::extendRight(*region, expandRight);
					}
				});
		njh::sort(beds,
				[](const std::shared_ptr<Bed3RecordCore> & region1, const std::shared_ptr<Bed3RecordCore> & region2) {
					return region1->chrom_ == region2->chrom_ ? (region1->chromStart_ == region2->chromStart_ ? region1->chromEnd_ < region2->chromEnd_ : region1->chromStart_ < region2->chromStart_): region1->chrom_ < region2->chrom_;
				});


		OutOptions mergedCoverageOpts(njh::files::make_path(chopParsCurrent.outputDirectory, "merged_coverage.bed"));
		{
			OutputStream mergedCoverageOut(mergedCoverageOpts);
			Bed3RecordCore currentRegion = *beds.front();
			for(const auto & regPos : iter::range<uint32_t>(1, beds.size())){
				if(currentRegion.chrom_ == beds[regPos]->chrom_ && currentRegion.overlaps(*beds[regPos], 1)){
					currentRegion.chromEnd_ = std::max(currentRegion.chromEnd_, beds[regPos]->chromEnd_);
				}else{
					Bed6RecordCore regOut(currentRegion.chrom_, currentRegion.chromStart_, currentRegion.chromEnd_, "", currentRegion.length(), '+');
					regOut.name_ = GenomicRegion(regOut).createUidFromCoords();

					mergedCoverageOut << regOut.toDelimStr() << std::endl;
					currentRegion = *beds[regPos];
				}
			}
			Bed6RecordCore regOut(currentRegion.chrom_, currentRegion.chromStart_, currentRegion.chromEnd_, "", currentRegion.length(), '+');
			regOut.name_ = GenomicRegion(regOut).createUidFromCoords();
			mergedCoverageOut << regOut.toDelimStr() << std::endl;
		}

		//refine regions
		refinePars.bedFnp = mergedCoverageOpts.outName();
		refinePars.bamFnp = njh::files::make_path(chopParsCurrent.outputDirectory, "fragments.sorted.bam");
		refinePars.outOpts = OutOptions(njh::files::make_path(chopParsCurrent.outputDirectory, "refined_merged.bed"));
		auto refinedRegions = RunRegionRefinement(refinePars);

		//extra fasta file of refined regions if 2bit file exists
		auto genomeTwoBitFnp = chopParsCurrent.genomeFnp;
		genomeTwoBitFnp.replace_extension(".2bit");
		if(bfs::exists(genomeTwoBitFnp)){
			auto refinedGRegions = bedPtrsToGenomicRegs(refinedRegions);
			if("" != gff){
				intersectBedLocsWtihGffRecordsPars interPars(gff,VecStr{"description"}, VecStr{"pseudogene", "gene"});
				refinePars.outOpts.overWriteFile_ = true;
				intersectBedLocsWtihGffRecords(refinedRegions, interPars);
				OutputStream bedOut(refinePars.outOpts);
				for(const auto & b : refinedRegions){
					bedOut << b->toDelimStrWithExtra() << std::endl;
				}
			}
			TwoBit::TwoBitFile tReader(genomeTwoBitFnp);
			auto refindedSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(chopParsCurrent.outputDirectory, "refined_merged.fasta"));
			SeqOutput refinedWriter(refindedSeqOpts);
			refinedWriter.openOut();
			for(const auto & reg : refinedGRegions){
				refinedWriter.write(reg.extractSeq(tReader));
			}
		}

		if(0 != minLength){
			std::vector<GenomicRegion> filteredRefinedRegions;
			{
				auto bedAgain = getBed3s(refinePars.outOpts.outName());
				OutOptions filteredOpts(njh::files::make_path(chopParsCurrent.outputDirectory, "filtered_refined_merged.bed"));
				OutputStream filteredOut(filteredOpts);
				for(const auto & b : bedAgain){
					if(b->length() >= minLength){
						filteredRefinedRegions.emplace_back(*b);
						filteredOut << b->toDelimStrWithExtra() << std::endl;
					}
				}
			}
			if(bfs::exists(genomeTwoBitFnp)){
				TwoBit::TwoBitFile tReader(genomeTwoBitFnp);
				auto refindedSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(chopParsCurrent.outputDirectory, "filtered_refined_merged.fasta"));
				SeqOutput refinedWriter(refindedSeqOpts);
				refinedWriter.openOut();
				for(const auto & reg : filteredRefinedRegions){
					refinedWriter.write(reg.extractSeq(tReader));
				}
			}
		}
	}


	return 0;
}

int seqSearchingRunner::chopAndMapAndRefine(const njh::progutils::CmdArgs & inputCommands) {
	ChopAndMapPars chopPars;
	CoverageFinderPars covPars;
	RegionRefinementPars refinePars;
	uint32_t expandLeft = 0;
	uint32_t expandRight = 0;
	uint32_t minLength = 0;
	covPars.coverageCutOff = 5;

	bfs::path gff = "";

	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	chopPars.debug = setUp.pars_.debug_;
	setUp.setOption(chopPars.numThreads, "--numThreads", "Number of threads to use");
	setUp.setOption(chopPars.genomeFnp, "--genomeFnp", "Genome to map to", true);
	setUp.setOption(gff, "--genomegff", "Genome gff file");

	setUp.setOption(chopPars.perFragmentCount, "--perFragmentCount", "perFragmentCount");
	setUp.setOption(chopPars.windowLength, "--windowLength", "windowLength");
	setUp.setOption(chopPars.windowStep, "--windowStep", "windowStep");
	setUp.setOption(refinePars.reOrient, "--reOrient", "reOrient");

	setUp.setOption(expandLeft, "--expandLeft", "expandLeft");
	setUp.setOption(expandRight, "--expandRight", "expandRight");

	setUp.setOption(covPars.coverageCutOff, "--coverageCutOff", "Coverage Cut Off");


	setUp.setOption(minLength, "--minLength", "minLength");

	setUp.processReadInNames(VecStr{"--fasta", "--fastq"});
	setUp.processDirectoryOutputName(true);
	chopPars.inOpts = setUp.pars_.ioOptions_;
	chopPars.outputDirectory = setUp.pars_.directoryName_;
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	//chop
	runChopAndMap(chopPars);
	//determine coverage
	covPars.numThreads = chopPars.numThreads;
	covPars.window = chopPars.windowLength;
	covPars.step = chopPars.windowStep;
	covPars.bams = njh::files::make_path(chopPars.outputDirectory, "fragments.sorted.bam").string();
	covPars.outOpts = OutOptions(njh::files::make_path(chopPars.outputDirectory, "coverage.bed"));
	RunCoverageFinder(covPars);
	//merge regions
	auto beds = getBed3s(covPars.outOpts.outName());
	njh::for_each(beds,
			[&expandLeft,&expandRight](const std::shared_ptr<Bed3RecordCore> & region) {
				if(0 != expandLeft) {
					BedUtility::extendLeft(*region, expandLeft);
				}
				if(0 != expandRight) {
					BedUtility::extendRight(*region, expandRight);
				}
			});
	njh::sort(beds,
			[](const std::shared_ptr<Bed3RecordCore> & region1, const std::shared_ptr<Bed3RecordCore> & region2) {
				return region1->chrom_ == region2->chrom_ ? (region1->chromStart_ == region2->chromStart_ ? region1->chromEnd_ < region2->chromEnd_ : region1->chromStart_ < region2->chromStart_): region1->chrom_ < region2->chrom_;
			});


	OutOptions mergedCoverageOpts(njh::files::make_path(chopPars.outputDirectory, "merged_coverage.bed"));
	{
		OutputStream mergedCoverageOut(mergedCoverageOpts);
		Bed3RecordCore currentRegion = *beds.front();
		for(const auto & regPos : iter::range<uint32_t>(1, beds.size())){
			if(currentRegion.chrom_ == beds[regPos]->chrom_ && currentRegion.overlaps(*beds[regPos], 1)){
				currentRegion.chromEnd_ = std::max(currentRegion.chromEnd_, beds[regPos]->chromEnd_);
			}else{
				Bed6RecordCore regOut(currentRegion.chrom_, currentRegion.chromStart_, currentRegion.chromEnd_, "", currentRegion.length(), '+');
				regOut.name_ = GenomicRegion(regOut).createUidFromCoords();

				mergedCoverageOut << regOut.toDelimStr() << std::endl;
				currentRegion = *beds[regPos];
			}
		}
		Bed6RecordCore regOut(currentRegion.chrom_, currentRegion.chromStart_, currentRegion.chromEnd_, "", currentRegion.length(), '+');
		regOut.name_ = GenomicRegion(regOut).createUidFromCoords();
		mergedCoverageOut << regOut.toDelimStr() << std::endl;
	}

	//refine regions
	refinePars.bedFnp = mergedCoverageOpts.outName();
	refinePars.bamFnp = njh::files::make_path(chopPars.outputDirectory, "fragments.sorted.bam");
	refinePars.outOpts = OutOptions(njh::files::make_path(chopPars.outputDirectory, "refined_merged.bed"));
	auto refinedRegions = RunRegionRefinement(refinePars);

	/**@todo report number of reads unmapped and what fragments these are*/

	//extra fasta file of refined regions if 2bit file exists
	auto genomeTwoBitFnp = chopPars.genomeFnp;
	genomeTwoBitFnp.replace_extension(".2bit");
	if(bfs::exists(genomeTwoBitFnp)){
		auto refinedGRegions = bedPtrsToGenomicRegs(refinedRegions);
		if("" != gff){
			intersectBedLocsWtihGffRecordsPars interPars(gff,VecStr{"description"}, VecStr{"pseudogene", "gene"});
			refinePars.outOpts.overWriteFile_ = true;
			intersectBedLocsWtihGffRecords(refinedRegions, interPars);
			OutputStream bedOut(refinePars.outOpts);
			for(const auto & b : refinedRegions){
				bedOut << b->toDelimStrWithExtra() << std::endl;
			}
		}
		TwoBit::TwoBitFile tReader(genomeTwoBitFnp);
		auto refindedSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(chopPars.outputDirectory, "refined_merged.fasta"));
		SeqOutput refinedWriter(refindedSeqOpts);
		refinedWriter.openOut();
		for(const auto & reg : refinedGRegions){
			refinedWriter.write(reg.extractSeq(tReader));
		}
	}

	if(0 != minLength){
		std::vector<GenomicRegion> filteredRefinedRegions;
		{
			auto bedAgain = getBed3s(refinePars.outOpts.outName());
			OutOptions filteredOpts(njh::files::make_path(chopPars.outputDirectory, "filtered_refined_merged.bed"));
			OutputStream filteredOut(filteredOpts);
			for(const auto & b : bedAgain){
				if(b->length() >= minLength){
					filteredRefinedRegions.emplace_back(*b);
					filteredOut << b->toDelimStrWithExtra() << std::endl;
				}
			}
		}
		if(bfs::exists(genomeTwoBitFnp)){
			TwoBit::TwoBitFile tReader(genomeTwoBitFnp);
			auto refindedSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(chopPars.outputDirectory, "filtered_refined_merged.fasta"));
			SeqOutput refinedWriter(refindedSeqOpts);
			refinedWriter.openOut();
			for(const auto & reg : filteredRefinedRegions){
				refinedWriter.write(reg.extractSeq(tReader));
			}
		}
	}

	return 0;
}


}  // namespace njhseq

