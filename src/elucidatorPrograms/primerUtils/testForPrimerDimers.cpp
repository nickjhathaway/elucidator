//
// Created by Nicholas Hathaway on 6/15/23.
//

#include "primerUtilsRunner.hpp"
#include <njhseq/IO/SeqIO/SeqInput.hpp>
#include <njhseq/IO/OutputStream.hpp>
#include <SeekDeep/objects/TarAmpSetupUtils/PrimersAndMids.hpp>
#include <SeekDeep/objects/IlluminaUtils/PairedReadProcessor.hpp>
#include <SeekDeep/objects/IlluminaUtils/PrimerDimerDetector.hpp>


namespace njhseq {

int primerUtilsRunner::testForPrimerDimers(const njh::progutils::CmdArgs &inputCommands) {

	OutOptions outOpts;
	bfs::path primersFnp;
	PrimerDimerDetector::DetectorPars detectorPars;
	seqSetUp setUp(inputCommands);
	PairedReadProcessor::ProcessParams stitchPars;
	stitchPars.minOverlap_ = 8;
	//primer dimers are short, the usual 0.02 allows no high quality mismatches in an overlap under 50 bases
	stitchPars.errorAllowed_ = 0.05;
	bool noTrimLowQualWindows = false;
	setUp.processVerbose();
	setUp.processDebug();
	setUp.pars_.ioOptions_.revComplMate_ = true;
	setUp.processReadInNames(concatVecs(njhseq::seqSetUp::singleInFormatsAvailable_, njhseq::seqSetUp::pairedReadInFormatsAvailable_), true);
	setUp.processWritingOptions(outOpts);
	setUp.processComparison(detectorPars.allowableErrors_);
	setUp.setOption(primersFnp, "--primers", "Primers file", true);
	setUp.setOption(detectorPars.maxGapBetweenPrimers_, "--maxGapBetweenPrimers", "Also call front and end primers that don't touch but are separated by up to this many bases, e.g. dimers that keep a piece of adapter between the primers; 0 only calls primers that overlap or abut");
	setUp.setOption(detectorPars.minPrimerMatchBases_, "--minPrimerMatchBases", "The fewest bases of a primer that have to match at a read's end for it to be trusted, 0 to work it out from the number of primers since the more primers there are the more likely a short match is by chance");
	setUp.setOption(detectorPars.minCandidateMatchBases_, "--minCandidateMatchBases", "The fewest bases of a primer that have to match at a read's end for it to be kept as a possible partner of a trusted primer at the other end, which makes the read an unresolved dimer");
	const std::string stitchGroup = "Stitching, only used with paired end input, where only the stitched R1 begins in R2 reads are tested";
	setUp.setOption(stitchPars.minOverlap_, "--minOverlap", "The minimal amount of overlap when stitching", false, stitchGroup);
	setUp.setOption(stitchPars.errorAllowed_, "--overLapErrorAllowed", "The fraction of high quality mismatches to allow in the overlap when stitching, rounds down so an overlap under 1/this bases is allowed none", false, stitchGroup);
	setUp.setOption(stitchPars.hardMismatchCutOff_, "--hardMismatchCutOff", "A hard cut off for number of mismatches in the overlap when stitching", false, stitchGroup);
	setUp.setOption(stitchPars.lqMismatchCutOff, "--lqMismatchCutOff", "A hard cut off for number of low quality mismatches in the overlap when stitching", false, stitchGroup);
	setUp.setOption(stitchPars.hqMismatchCutOff, "--hqMismatchCutOff", "A hard cut off for number of high quality mismatches in the overlap when stitching", false, stitchGroup);
	setUp.setOption(stitchPars.matchScore_, "--stitchMatchScore", "The score for a match when aligning the mates to stitch them", false, stitchGroup);
	setUp.setOption(stitchPars.mismatchPenalty_, "--stitchMismatchPenalty", "The penalty for a mismatch when aligning the mates to stitch them, harsher than the match score so a long low identity placement of the mates doesn't out score a short true overlap", false, stitchGroup);
	setUp.setOption(stitchPars.r1Trim_, "--r1Trim", "Remove this many bases off of the end of r1 reads", false, stitchGroup);
	setUp.setOption(stitchPars.r2Trim_, "--r2Trim", "Remove this many bases off of the end of r2 reads", false, stitchGroup);
	setUp.setOption(stitchPars.qualWindowPar_.avgQualCutOff_, "--qWindowTrimAvgQualCutOff", "Quality Window Trim Avg Qual Cut Off", false, stitchGroup);
	setUp.setOption(stitchPars.qualWindowPar_.windowSize_, "--qWindowSize", "Quality Window Trim Size", false, stitchGroup);
	setUp.setOption(stitchPars.qualWindowPar_.windowStep_, "--qWindowStep", "Quality Window Trim Step", false, stitchGroup);
	setUp.setOption(noTrimLowQualWindows, "--noTrimLowQualWindows", "Don't Trim Low Qual Windows", false, stitchGroup);
	stitchPars.trimLowQaulWindows_ = !noTrimLowQualWindows;
	stitchPars.verbose_ = setUp.pars_.verbose_;
	setUp.finishSetUp(std::cout);

	PrimersAndMids ids(primersFnp);

	if(ids.getTargets().empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error in reading in target primers file " << primersFnp << "\n";
		ss << "Make sure there is a line that starts with target in file" << "\n";
		throw std::runtime_error{ss.str()};
	}
	ids.initPrimerDeterminator();
	PrimerDimerDetector detector(*ids.pDeterminator_, detectorPars);
	const bool pairedInput = setUp.pars_.ioOptions_.isPairedIn();
	uint64_t maxSize = 0;
	uint64_t maxMateSize = 0;
	if(pairedInput){
		PairedRead pair;
		SeqInput reader(setUp.pars_.ioOptions_);
		reader.openIn();
		while(reader.readNextRead(pair)){
			readVec::getMaxLength(pair.seqBase_, maxMateSize);
			readVec::getMaxLength(pair.mateSeqBase_, maxMateSize);
		}
		//a stitched read can't be longer than both of its mates laid end to end
		maxSize = 2 * maxMateSize;
	} else {
		seqInfo seq;
		SeqInput reader(setUp.pars_.ioOptions_);
		reader.openIn();
		while(reader.readNextRead(seq)){
			readVec::getMaxLength(seq, maxSize);
		}
	}

	aligner alignerObj = detector.genAligner(maxSize);

	OutputStream out(outOpts);
	out
					<< "filename"
					<< "\t" << "name"
					<< "\t" << "seq"
					<< "\t" << "possibleDimers"
					<< "\t" << "numberOfFrontMatches"
					<< "\t" << "frontMatchesNames"
					<< "\t" << "numberOfEndMatches"
					<< "\t" << "endMatchesNames"
					<< "\t" << "overlapIDNumber"
					<< "\t" << "frontDimerPrimer-endDimerPrimer"
					<< "\t" << "degreeOfOverlap"
					<< "\t" << "frontDimerPrimer"
					<< "\t" << "frontDimerPrimer_isRevComp"
					<< "\t" << "frontDimerPrimer_seqStart"
					<< "\t" << "frontDimerPrimer_seqEnd"
					<< "\t" << "frontDimerPrimer_primerStart"
					<< "\t" << "frontDimerPrimer_primerEnd"

					<< "\t" << "endDimerPrimer"
					<< "\t" << "endDimerPrimer_isRevComp"
					<< "\t" << "endDimerPrimer_seqStart"
					<< "\t" << "endDimerPrimer_seqEnd"
					<< "\t" << "endDimerPrimer_primerStart"
					<< "\t" << "endDimerPrimer_primerEnd"
					<< "\t" << "gapBetweenPrimers"
					<< "\t" << "readCategory"
					<< "\t" << "candidatePartners"

					<< std::endl;
	auto filenameOut = bfs::basename(bfs::path(setUp.pars_.ioOptions_.firstName_).replace_extension(""));
	std::unordered_map<PrimerDimerDetector::ReadCategory, uint32_t> categoryCounts;
	const auto testForDimers = [&](const seqInfo & seq){
		const auto res = detector.detect(seq, alignerObj);
		++categoryCounts[res.category_];
		VecStr frontNames;
		VecStr endNames;
		for(const auto & front : res.frontMatches_){
			frontNames.emplace_back(front.primerName_);
		}
		for(const auto & end : res.endMatches_){
			endNames.emplace_back(end.primerName_);
		}
		const auto categoryStr = PrimerDimerDetector::getReadCategoryStr(res.category_);
		VecStr candidatePartners;
		for(const auto & candidate : res.candidatePartners_){
			candidatePartners.emplace_back(njh::pasteAsStr(candidate.primerName_, ":", candidate.primerBasesMatched()));
		}
		if(res.dimers_.empty()){
			out << filenameOut << "\t" << seq.name_
					<< "\t" << seq.seq_
					<< "\t" << res.dimers_.size()
					<< "\t" << res.frontMatches_.size()
					<< "\t" << njh::conToStr(frontNames, ";")
					<< "\t" << res.endMatches_.size()
					<< "\t" << njh::conToStr(endNames, ";")
					<< "\t" << ""
					<< "\t" << ""
					<< "\t" << ""
					<< "\t" << ""
					<< "\t" << ""
					<< "\t" << ""
					<< "\t" << ""
					<< "\t" << ""
					<< "\t" << ""

					<< "\t" << ""
					<< "\t" << ""
					<< "\t" << ""
					<< "\t" << ""
					<< "\t" << ""
					<< "\t" << ""
					<< "\t" << ""
					<< "\t" << categoryStr
					<< "\t" << njh::conToStr(candidatePartners, ";")

					<< std::endl;
		} else {
			uint32_t overlapCount = 0;
			for (const auto &dimer: res.dimers_) {
				out << filenameOut << "\t" << seq.name_
						<< "\t" << seq.seq_
						<< "\t" << res.dimers_.size()
						<< "\t" << res.frontMatches_.size()
						<< "\t" << njh::conToStr(frontNames, ";")
						<< "\t" << res.endMatches_.size()
						<< "\t" << njh::conToStr(endNames, ";")
						<< "\t" << overlapCount
						<< "\t" << dimer.front_.primerName_ << "-" << dimer.end_.primerName_
						<< "\t" << dimer.getOverlap()
						<< "\t" << dimer.front_.primerName_
						<< "\t" << njh::boolToStr(dimer.front_.revComp_)
						<< "\t" << dimer.front_.seqStart_
						<< "\t" << dimer.front_.seqEnd_
						<< "\t" << dimer.front_.primerStart_
						<< "\t" << dimer.front_.primerEnd_

						<< "\t" << dimer.end_.primerName_
						<< "\t" << njh::boolToStr(dimer.end_.revComp_)
						<< "\t" << dimer.end_.seqStart_
						<< "\t" << dimer.end_.seqEnd_
						<< "\t" << dimer.end_.primerStart_
						<< "\t" << dimer.end_.primerEnd_
						<< "\t" << dimer.getGap()
						<< "\t" << categoryStr
						<< "\t" << ""

						<< std::endl;
				++overlapCount;
			}
		}
		if(setUp.pars_.debug_){
			seq.outPutSeqAnsi(std::cout);
			std::cout << "category: " << categoryStr << std::endl;
			std::cout << "dimers: " << res.dimers_.size() << std::endl;
			std::cout << "frontMatches: " << res.frontMatches_.size() << std::endl;
			std::cout << "endMatches: " << res.endMatches_.size() << std::endl;
			for(const auto & dimer : res.dimers_){
				seq.outPutSeqAnsi(std::cout);
				std::cout << dimer.front_.alignedPrimer_.getSeqAnsi() << std::endl;
				std::cout << std::string(dimer.end_.seqStart_, ' ') << dimer.end_.alignedPrimer_.getSeqAnsi() << std::endl;
				std::cout << dimer.front_.primerName_ << " " << dimer.front_.seqStart_ << ":" << dimer.front_.seqEnd_ << std::endl;
				std::cout << dimer.end_.primerName_ << " " << dimer.end_.seqStart_ << ":" << dimer.end_.seqEnd_ << std::endl;
				std::cout << "overlap: " << dimer.getOverlap() << std::endl;
				std::cout << "gap: " << dimer.getGap() << std::endl;
				std::cout << "alnScoreSum: " << dimer.getAlnScoreSum() << std::endl;
				std::cout << "identitySum: " << dimer.getIdentitySum() << std::endl;
				std::cout << std::endl;
			}
			std::cout << std::endl;
		}
	};

	if(pairedInput){
		std::unique_ptr<OutputStream> stitchCountsFile;
		if(!outOpts.outFilename_.empty()){
			OutOptions stitchCountsOpts(bfs::path(outOpts.outFilename_.string() + "_stitchCounts.tsv"));
			stitchCountsOpts.transferOverwriteOpts(outOpts);
			stitchCountsFile = std::make_unique<OutputStream>(stitchCountsOpts);
		}
		PairedReadProcessor pairProcessor(stitchPars);
		aligner stitchingAligner = pairProcessor.genStitchingAligner(maxMateSize);
		stitchingAligner.qScorePars_.qualThresWindow_ = 0;
		PairedReadProcessor::ProcessedResultsCounts stitchCounts;
		PairedRead pair;
		SeqInput reader(setUp.pars_.ioOptions_);
		reader.openIn();
		while(reader.readNextRead(pair)){
			++stitchCounts.total;
			auto stitched = pairProcessor.processPairedEnd(pair, stitchCounts, stitchingAligner);
			//primer dimers are short enough for the mates to read through each other
			if(PairedReadProcessor::ReadPairOverLapStatus::R1BEGINSINR2 == stitched.status_){
				testForDimers(*stitched.combinedSeq_);
			}
		}
		const auto dimersCalled = categoryCounts[PrimerDimerDetector::ReadCategory::DIMER];
		std::ostream & stitchCountsOut = stitchCountsFile ? static_cast<std::ostream &>(*stitchCountsFile) : std::cerr;
		stitchCountsOut
						<< "filename"
						<< "\t" << "total"
						<< "\t" << "notCombined"
						<< "\t" << "perfectOverlapCombined"
						<< "\t" << "r1EndsInR2Combined"
						<< "\t" << "r1BeginsInR2Combined"
						<< "\t" << "r1AllInR2Combined"
						<< "\t" << "r2AllInR1Combined"
						<< "\t" << "possibleDimersCalled"
						<< "\t" << "possibleDimersCalledPerc"
						<< "\t" << "primersFarApart"
						<< "\t" << "unresolvedDimers"
						<< "\t" << "unresolvedDimersPerc"
						<< "\t" << "onePrimer"
						<< "\t" << "noPrimer"
						<< "\t" << "minPrimerMatchBases"
						<< std::endl;
		stitchCountsOut
						<< filenameOut
						<< "\t" << stitchCounts.total
						<< "\t" << stitchCounts.overlapFail
						<< "\t" << stitchCounts.perfectOverlapCombined
						<< "\t" << stitchCounts.r1EndsInR2Combined
						<< "\t" << stitchCounts.r1BeginsInR2Combined
						<< "\t" << stitchCounts.r1AllInR2Combined
						<< "\t" << stitchCounts.r2AllInR1Combined
						<< "\t" << dimersCalled
						<< "\t" << (0 == stitchCounts.total ? 0.0 : 100.0 * dimersCalled / stitchCounts.total)
						<< "\t" << categoryCounts[PrimerDimerDetector::ReadCategory::PRIMERS_FAR_APART]
						<< "\t" << categoryCounts[PrimerDimerDetector::ReadCategory::UNRESOLVED_DIMER]
						<< "\t" << (0 == stitchCounts.total ? 0.0 : 100.0 * categoryCounts[PrimerDimerDetector::ReadCategory::UNRESOLVED_DIMER] / stitchCounts.total)
						<< "\t" << categoryCounts[PrimerDimerDetector::ReadCategory::ONE_PRIMER]
						<< "\t" << categoryCounts[PrimerDimerDetector::ReadCategory::NO_PRIMER]
						<< "\t" << detector.pars_.minPrimerMatchBases_
						<< std::endl;
	} else {
		seqInfo seq;
		SeqInput reader(setUp.pars_.ioOptions_);
		reader.openIn();
		while(reader.readNextRead(seq)){
			testForDimers(seq);
		}
	}

	return 0;
}

} // namespace njhseq
