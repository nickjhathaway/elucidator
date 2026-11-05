/*
 * genExp_extractRefSeqsFromGenomesWithFlanks.cpp
 *
 *  Created on: Nov 3, 2026
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

#include <TwoBit.h>
#include "genExp.hpp"
#include "elucidator/objects/BioDataObject.h"
#include "elucidator/objects/dataContainers/graphs.h"
#include <njhseq/GenomeUtils.h>
#include <njhseq/readVectorManipulation/readVectorHelpers/readVecTrimmer.hpp>

namespace njhseq {

//defined in genExp_extractRefSeqsFromGenomes.cpp
std::vector<GenomicRegion> mergeRegionsStrandAware(const std::vector<GenomicRegion> & regions);

namespace {

/**@brief the outcome of trying to pair up the two flank hits for one genome
 */
struct FlankPairResult {
	std::string status_{"ok"};
	bool usable_{false};
	uint32_t leftHits_{0};
	uint32_t rightHits_{0};
	GenomicRegion span_;
	uint32_t fivePrimeTrim_{0};
	uint32_t threePrimeTrim_{0};
};

/**@brief carve a flank out of a region, in genomic (left/right) coordinates
 *
 * outset shifts the flank outwards from the region, 0 puts the flank on the
 * first/last flankLength bases of the region itself
 */
GenomicRegion genFlankRegion(const GenomicRegion & region, uint32_t flankLength,
		uint32_t outset, bool leftSide, uint32_t chromLength) {
	GenomicRegion flank = region;
	if (leftSide) {
		size_t flankStart = region.start_ > outset ? region.start_ - outset : 0;
		flank.start_ = flankStart;
		flank.end_ = flankStart + flankLength;
	} else {
		size_t flankEnd = region.end_ + outset;
		if (flankEnd > chromLength) {
			flankEnd = chromLength;
		}
		flank.end_ = flankEnd;
		flank.start_ = flankEnd > flankLength ? flankEnd - flankLength : 0;
	}
	if (flank.end_ > chromLength) {
		flank.end_ = chromLength;
	}
	flank.uid_ = njh::pasteAsStr(leftSide ? "leftFlank" : "rightFlank");
	return flank;
}

/**@brief strongest tandem periodicity in a sequence
 *
 * the largest fraction, over periods 1..maxPeriod, of positions where the base
 * repeats one period later. A perfect tandem repeat approaches 1, while ordinary
 * AT rich Plasmodium coding sequence sits near 0.4. Only used to skip candidate
 * anchors cheaply before paying for a lastz pass, so the threshold is kept
 * conservative; it is deliberately not the thing the search optimises, because
 * imperfect or long period repeats score low here even when they are unanchorable
 */
double bestTandemPeriodicity(const std::string & seq, uint32_t maxPeriod = 60) {
	double best = 0;
	uint32_t topPeriod = std::min<uint32_t>(maxPeriod, seq.size() / 2);
	for (uint32_t period = 1; period <= topPeriod; ++period) {
		uint32_t matches = 0;
		for (uint32_t pos = 0; pos + period < seq.size(); ++pos) {
			if (seq[pos] == seq[pos + period]) {
				++matches;
			}
		}
		double frac = static_cast<double>(matches) / (seq.size() - period);
		if (frac > best) {
			best = frac;
		}
	}
	return best;
}

/**@brief how many bases to shave off one end of an extracted span to land back on
 * the edge of the original bed record
 *
 * The reference is a window of the primary genome that straddles the edge: it
 * reaches from the outer end of the anchor to some way inside the record, so the
 * edge always sits at a known index within it. Aligning that window against the
 * matching end of the span, and asking where that index landed, puts the record
 * edge back without ever aligning the whole locus. Doing it on the whole locus is the obvious reading of
 * trimToRefWithGlobalAlignmentToRefPositions, but the aligner carries a
 * maxSize_ x maxSize_ score matrix of 16 byte cells and setMaxSize only ever
 * grows it, so a 23 kb record would claim ~8.8 GB and hold it for the rest of the
 * run. Only the edges are moving, so only the edges need aligning.
 *
 * @param refSeq the primary genome window straddling this edge, in record sense
 * @param spanSeq the extracted span, reading in the same direction as the record
 * @param edgeIdx index within refSeq of the record's edge base
 * @param fivePrimeSide whether this is the start of the record rather than its end
 * @param ok set false when the alignment could not place the edge
 * @return bases to remove from that end of spanSeq
 */
uint32_t findEdgeTrim(const seqInfo & refSeq, const seqInfo & spanSeq,
		uint32_t edgeIdx, bool fivePrimeSide, aligner & alignerObj, bool & ok,
		bool verbose = false) {
	ok = false;
	uint32_t anchorLen = len(refSeq);
	if (0 == anchorLen || edgeIdx >= anchorLen) {
		return 0;
	}
	//take only the end of the span that holds this edge, with room for the anchor
	//to have picked up an indel on the way
	uint32_t chunkLen = std::min<uint32_t>(len(spanSeq), 2 * anchorLen);
	if (chunkLen < anchorLen) {
		return 0;
	}
	seqInfo chunk = spanSeq;
	if (fivePrimeSide) {
		chunk.setClip(0, chunkLen - 1);
	} else {
		chunk.setClip(len(spanSeq) - chunkLen, len(spanSeq) - 1);
	}
	readVecTrimmer::GlobalAlnTrimPars trimPars{};
	if (fivePrimeSide) {
		//keep from the record's first base to the inner end of the reference
		trimPars.startInclusive_ = edgeIdx;
		trimPars.endInclusive_ = anchorLen - 1;
	} else {
		//keep from the inner end of the reference up to the record's last base
		trimPars.startInclusive_ = 0;
		trimPars.endInclusive_ = edgeIdx;
	}
	uint64_t needed = std::max<uint64_t>(len(chunk), anchorLen) + 10;
	if (needed >= alignerObj.parts_.maxSize_) {
		alignerObj.parts_.setMaxSize(needed);
	}
	auto trimmed = chunk;
	readVecTrimmer::trimSeqToRefByGlobalAln(trimmed, refSeq, trimPars, alignerObj);
	if (!trimmed.on_ || trimmed.seq_.empty()) {
		if (verbose) {
			std::cerr << "  trim " << (fivePrimeSide ? "5p" : "3p")
					<< ": alignment could not place the edge (on_="
					<< njh::boolToStr(trimmed.on_) << ", trimmedLen=" << trimmed.seq_.size()
					<< ", refLen=" << anchorLen << ", chunkLen=" << chunkLen
					<< ", refRange=" << trimPars.startInclusive_ << "-"
					<< trimPars.endInclusive_ << ")" << std::endl;
		}
		return 0;
	}
	auto pos = chunk.seq_.find(trimmed.seq_);
	if (std::string::npos == pos) {
		if (verbose) {
			std::cerr << "  trim " << (fivePrimeSide ? "5p" : "3p")
					<< ": trimmed piece is not a substring of the chunk" << std::endl;
		}
		return 0;
	}
	ok = true;
	if (fivePrimeSide) {
		return pos;
	}
	return chunkLen - (pos + len(trimmed));
}

} // namespace

int genExpRunner::extractRefSeqsFromGenomesWithFlanks(
		const njh::progutils::CmdArgs & inputCommands) {

	MultiGenomeMapper::inputParameters genomeMappingPars;
	bfs::path bedFile = "";
	bfs::path outputDir = "";
	bool keepRefAlignments = false;
	bool overWriteDirs = false;
	bool combineByGenome = false;
	bool writeOutAllSeqsFile = false;

	uint32_t requestedFlankLength = 200;
	uint32_t minFlankLength = 50;
	uint32_t flankOutset = 0;
	bool autoFlankSearch = false;
	uint32_t autoFlankStep = 25;
	uint32_t autoFlankMaxOutset = 600;
	uint32_t autoFlankMaxAmbiguousGenomes = 0;
	double autoFlankMaxPeriodicity = 0.70;
	double minSpanMultiplier = 0.25;
	double maxSpanMultiplier = 3.0;
	uint32_t maxSpanPad = 2000;
	bool allowAmbiguousFlanks = false;
	bool noTrimToBedEdges = false;

	MultiGenomeMapper::getRefSeqsWithPrimaryGenomePars extractPars;
	extractPars.lzPars.identity = 80;
	extractPars.lzPars.coverage = 90;

	bfs::path gffDir = "";

	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.description_ =
			"Extract a region from several genomes by anchoring on the two flanks of a bed "
			"record rather than aligning the whole record. Intended for regions whose middle "
			"is a tandem repeat that varies in length between genomes, where a whole-region "
			"alignment either fails the lastz coverage cut off or slides within the repeat.";

	setUp.setOption(requestedFlankLength, "--flankLength",
			"Length of each flank anchor to search with; shrunk automatically for a record too short to carry two of them");
	setUp.setOption(minFlankLength, "--minFlankLength",
			"Never shrink an anchor below this; a record too short for two of these is skipped and listed in skippedRegions.tab.txt");
	setUp.setOption(flankOutset, "--flankOutset",
			"Shift each flank this many bases outwards from the bed record; 0 means the anchors are the first and last --flankLength bases of the record itself. Raise it when the repeat runs all the way to the edge of the record");
	setUp.setOption(autoFlankSearch, "--autoFlankSearch",
			"Choose the outset for each side automatically instead of using one fixed --flankOutset: step the anchor outwards until it hits a single locus in every genome. Needed when the repeat runs to the edge of the record, where the two sides usually need different shifts");
	setUp.setOption(autoFlankStep, "--autoFlankStep",
			"Step size when searching outwards with --autoFlankSearch");
	setUp.setOption(autoFlankMaxOutset, "--autoFlankMaxOutset",
			"Stop searching outwards with --autoFlankSearch once the outset reaches this");
	setUp.setOption(autoFlankMaxAmbiguousGenomes, "--autoFlankMaxAmbiguousGenomes",
			"With --autoFlankSearch, accept an anchor once no more than this many genomes still have an ambiguous or missing hit");
	setUp.setOption(autoFlankMaxPeriodicity, "--autoFlankMaxPeriodicity",
			"With --autoFlankSearch, skip a candidate anchor without searching when its strongest tandem periodicity exceeds this; purely a speed up, the acceptance decision is made on the lastz hit counts");
	setUp.setOption(minSpanMultiplier, "--minSpanMultiplier",
			"Reject a genome whose paired span is shorter than this multiple of the primary genome region length");
	setUp.setOption(maxSpanMultiplier, "--maxSpanMultiplier",
			"Reject a genome whose paired span is longer than this multiple of the primary genome region length (plus --maxSpanPad)");
	setUp.setOption(maxSpanPad, "--maxSpanPad",
			"Added to the --maxSpanMultiplier bound before rejecting a span");
	setUp.setOption(allowAmbiguousFlanks, "--allowAmbiguousFlanks",
			"Keep a genome even when a flank hit more than one locus, using the best hit; off by default because paralogous flanks can otherwise pair across two different loci");
	setUp.setOption(noTrimToBedEdges, "--noTrimToBedEdges",
			"Report the whole anchor to anchor span instead of trimming it back to the edges of the bed record. On by default: when an anchor is pushed outwards, by --flankOutset or by --autoFlankSearch, the span starts and ends outside the record, so each extraction is reported with that outset still on it. Trimming aligns each anchor to the matching end of the span to put the record edge back, which also keeps the reported length ratio honest. Turn this off when you want the outset sequence kept");
	const bool trimToBedEdges = !noTrimToBedEdges;

	setUp.setOption(extractPars.shortNames, "--shortNames", "Short Names");
	setUp.setOption(writeOutAllSeqsFile, "--writeOutAllSeqsFile", "Write Out All Seqs File");
	setUp.setOption(combineByGenome, "--combineByGenome",
			"Combine and merge regions extracted and create a bed file for each genome");
	setUp.setOption(extractPars.extendAndTrim, "--extendAndTrim",
			"Refine each individual FLANK HIT by extending it and trimming it back to the flank sequence, which can help when variation falls at the very ends of an anchor. This acts on the anchors, not on the final locus: trimming the locus back to the bed record is --noTrimToBedEdges");
	setUp.setOption(extractPars.extendAndTrimLen, "--extendAndTrimLen",
			"Length to extend by when refining flank hits with --extendAndTrim");
	setUp.setOption(extractPars.lzPars.identity, "--identity",
			"Identity to use for when searching with lastz");
	setUp.setOption(extractPars.lzPars.coverage, "--coverage",
			"Coverage to use for when searching with lastz; this now applies to each flank rather than the whole region, so the default is appropriate");
	setUp.setOption(extractPars.byScore, "--keepBestOnlyByScore",
			"Determine best flank hit by alignment score instead of percent identity");
	setUp.setOption(genomeMappingPars.numThreads_, "--numThreads", "Number of CPUs to utilize");
	setUp.setOption(gffDir, "--gffDir",
			"A directory with a gff for the genomes in --genomeDir, should be named GENOME.gff (for GENOME.fasta)");
	setUp.setOption(keepRefAlignments, "--keepRefAlignments", "Keep Ref Alignments");

	setUp.setOption(bedFile, "--bed", "Bed file of the regions to extract", true);
	setUp.setOption(genomeMappingPars.genomeDir_, "--genomeDir", "Name of the genome file fnp", true);
	setUp.setOption(genomeMappingPars.primaryGenome_, "--primaryGenome",
			"Name of the primary genome, should be in --genomeDir", true);
	std::string selectedGenomesStr;
	setUp.setOption(selectedGenomesStr, "--selectedGenomes",
			"Name of the other genomes in --genomeDir to be read in, leave blank to just do all fastas");
	//comma separated gff attributes to carry onto the hits found in a genome that
	//has a gff in --gffDir. MultiGenomeMapper already defaults this to
	//"description", so an empty value here means keep that, not take nothing.
	//extractRefSeqsFromGenomes also exposes a --gffExtraAttributes that is read
	//into a local string and never used again, which is why it is absent here.
	std::string extraGffAttributes;
	setUp.setOption(extraGffAttributes, "--extraGffAttributes",
			"Comma separated extra attributes to add for a genome with an accompanying gff in --gffDir; leave blank to keep the default of description");
	if (!extraGffAttributes.empty()) {
		genomeMappingPars.gffIntersectPars_.extraAttributes_ = tokenizeString(extraGffAttributes, ",");
	}
	setUp.setOption(overWriteDirs, "--overWriteDirs", "Over write region dirs");
	setUp.setOption(outputDir, "--outputDir",
			"Output Directory to write files, will be created if it doesn't already exist", true);
	setUp.finishSetUp(std::cout);

	if (requestedFlankLength < 20) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error --flankLength of " << requestedFlankLength
				<< " is too short to anchor reliably, use at least 20" << "\n";
		throw std::runtime_error { ss.str() };
	}

	njh::files::makeDirP(njh::files::MkdirPar(outputDir));
	setUp.startARunLog(outputDir.string());

	auto regions = gatherRegions(bedFile.string(), "", setUp.pars_.verbose_);
	if (regions.empty()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error no regions read from " << bedFile << "\n";
		throw std::runtime_error { ss.str() };
	}

	std::unique_ptr<MultiGenomeMapper> gMapper = std::make_unique<MultiGenomeMapper>(genomeMappingPars);
	gMapper->pars_.numThreads_ = genomeMappingPars.numThreads_;
	gMapper->pars_.verbose_ = setUp.pars_.verbose_;
	gMapper->pars_.gffDir_ = gffDir;
	gMapper->setSelectedGenomes(selectedGenomesStr);
	gMapper->init();

	const std::string primaryGenome = gMapper->pars_.primaryGenome_;
	if (!gMapper->hasGenome(primaryGenome)) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error primary genome " << primaryGenome
				<< " was not loaded" << "\n";
		throw std::runtime_error { ss.str() };
	}

	//the flank search always keeps every hit so that ambiguity can be reported,
	//the best hit is still the first one returned
	MultiGenomeMapper::getRefSeqsWithPrimaryGenomePars flankPars = extractPars;
	flankPars.keepBestOnly = false;
	flankPars.shortNames = true;

	auto gapPars = gapScoringParameters(5, 1, 0, 0, 0, 0);
	auto scoring = substituteMatrix(2, -2);
	uint64_t maxAlignLen = requestedFlankLength + extractPars.extendAndTrimLen * 2 + 20;
	aligner alignerObj(maxAlignLen, gapPars, scoring);

	auto genomeNames = getVectorOfMapKeys(gMapper->genomes_);
	njh::sort(genomeNames);

	//regions that could not be processed, so one bad record cannot sink a whole batch
	VecStr skippedRegions;

	for (auto & region : regions) {
		std::string rUid = region.createUidFromCoordsStrand();
		//setUidWtihCoordsStrand throws the bed record's name away, so keep it for the
		//bed output: a regions.bed holding every genome is only readable if each row
		//still says which locus it is
		std::string regionName = region.uid_.empty() ? rUid : region.uid_;
		region.setUidWtihCoordsStrand();
		auto regionDir = njh::files::make_path(outputDir, rUid);
		if (bfs::exists(regionDir) && !overWriteDirs) {
			if (setUp.pars_.verbose_) {
				std::cout << "skipping " << rUid << ", directory already exists" << std::endl;
			}
			continue;
		}
		//one unextractable record should not take the rest of the batch down with it
		try {
		auto regDirPars = njh::files::MkdirPar(regionDir);
		regDirPars.overWriteDir_ = true;
		njh::files::makeDir(regDirPars);
		auto alignsDir = njh::files::makeDir(regionDir, njh::files::MkdirPar("refAlignments"));
		auto bedsDir = njh::files::makeDir(regionDir, njh::files::MkdirPar("beds"));

		uint32_t chromLength = gMapper->genomes_.at(primaryGenome)->chromosomeLengths_.at(region.chrom_);

		//a short region cannot carry two non-overlapping anchors of the requested
		//length, so shrink them for that region rather than failing the whole run
		uint32_t flankLength = requestedFlankLength;
		if (2 * flankLength + 20 > region.getLen()) {
			flankLength = region.getLen() > 20 ? (region.getLen() - 20) / 2 : 0;
			if (flankLength < minFlankLength) {
				std::cerr << "skipping " << rUid << ", only " << region.getLen()
						<< " bp, too short for two anchors of at least --minFlankLength "
						<< minFlankLength << std::endl;
				skippedRegions.emplace_back(njh::pasteAsStr(rUid, "\ttoo_short_for_flanks"));
				continue;
			}
			if (setUp.pars_.verbose_) {
				std::cout << rUid << " is " << region.getLen()
						<< " bp, shrinking the anchors to " << flankLength << " bp" << std::endl;
			}
		}

		//search one flank at a given outset. Each flank has to be searched on its own
		//because determineRegionsLastzAgainstGenomes pools hits from every input
		//sequence into one list per genome, so a single fasta holding both sides would
		//lose which side a hit came from
		auto searchFlank =
				[&gMapper,&region,&flankLength,&chromLength,&primaryGenome,&alignsDir,&flankPars,&alignerObj](
						bool leftSide, uint32_t outset) {
					const std::string side = leftSide ? "left" : "right";
					auto flank = genFlankRegion(region, flankLength, outset, leftSide, chromLength);
					auto sideDir = njh::files::makeDir(alignsDir,
							njh::files::MkdirPar(njh::pasteAsStr(side, "_outset", outset)));
					auto flankSeq = gMapper->extractGenomeRegion(primaryGenome, flank);
					flankSeq.name_ = flank.uid_;
					auto flankFnp = njh::files::make_path(sideDir, side + "Flank.fasta");
					SeqOutput::write(std::vector<seqInfo> { flankSeq },
							SeqIOOptions::genFastaOut(flankFnp));
					return gMapper->determineRegionsLastzAgainstGenomes(
							SeqIOOptions::genFastaIn(flankFnp), sideDir, flankPars, alignerObj);
				};
		//how many genomes this anchor leaves unusable, which is what the search minimises
		auto countBadGenomes =
				[&genomeNames](
						const std::unordered_map<std::string, MultiGenomeMapper::determineRegionsLastzRes> & res) {
					uint32_t bad = 0;
					for (const auto & genome : genomeNames) {
						uint32_t nHits =
								njh::in(genome, res) ? res.at(genome).all_individual_regions.size() : 0;
						if (1 != nHits) {
							++bad;
						}
					}
					return bad;
				};

		std::unordered_map<std::string, MultiGenomeMapper::determineRegionsLastzRes> leftRes;
		std::unordered_map<std::string, MultiGenomeMapper::determineRegionsLastzRes> rightRes;
		uint32_t leftOutsetUsed = flankOutset;
		uint32_t rightOutsetUsed = flankOutset;

		for (const bool leftSide : std::vector<bool> { true, false }) {
			//the sides are searched independently because the repeat often runs to only
			//one edge of the record, so they need different shifts
			uint32_t chosenOutset = flankOutset;
			auto bestRes = searchFlank(leftSide, flankOutset);
			uint32_t bestBad = countBadGenomes(bestRes);
			if (autoFlankSearch && bestBad > autoFlankMaxAmbiguousGenomes) {
				for (uint32_t outset = flankOutset + autoFlankStep;
						outset <= autoFlankMaxOutset; outset += autoFlankStep) {
					auto cand = genFlankRegion(region, flankLength, outset, leftSide, chromLength);
					if (cand.getLen() < flankLength) {
						break; //ran off the end of the chromosome
					}
					auto candSeq = gMapper->extractGenomeRegion(primaryGenome, cand);
					if (bestTandemPeriodicity(candSeq.seq_) > autoFlankMaxPeriodicity) {
						continue; //still inside the repeat, not worth a lastz pass
					}
					auto res = searchFlank(leftSide, outset);
					uint32_t bad = countBadGenomes(res);
					if (bad < bestBad) {
						bestBad = bad;
						bestRes = res;
						chosenOutset = outset;
					}
					if (bestBad <= autoFlankMaxAmbiguousGenomes) {
						break;
					}
				}
				if (setUp.pars_.verbose_) {
					std::cout << rUid << " " << (leftSide ? "left" : "right")
							<< " flank: outset " << chosenOutset << ", " << bestBad
							<< " genome(s) still unusable" << std::endl;
				}
			}
			if (leftSide) {
				leftRes = bestRes;
				leftOutsetUsed = chosenOutset;
			} else {
				rightRes = bestRes;
				rightOutsetUsed = chosenOutset;
			}
		}

		auto leftFlank = genFlankRegion(region, flankLength, leftOutsetUsed, true, chromLength);
		auto rightFlank = genFlankRegion(region, flankLength, rightOutsetUsed, false, chromLength);
		if (leftFlank.end_ > rightFlank.start_) {
			std::cerr << "skipping " << rUid
					<< ", the chosen anchors overlap after the outward search" << std::endl;
			skippedRegions.emplace_back(njh::pasteAsStr(rUid, "\toverlapping_flanks"));
			continue;
		}

		//pair the two sides up per genome
		std::unordered_map<std::string, FlankPairResult> pairedByGenome;
		for (const auto & genome : genomeNames) {
			FlankPairResult pr;
			std::vector<GenomicRegion> lHits;
			std::vector<GenomicRegion> rHits;
			if (njh::in(genome, leftRes)) {
				lHits = leftRes.at(genome).all_individual_regions;
			}
			if (njh::in(genome, rightRes)) {
				rHits = rightRes.at(genome).all_individual_regions;
			}
			pr.leftHits_ = lHits.size();
			pr.rightHits_ = rHits.size();
			if (lHits.empty() && rHits.empty()) {
				pr.status_ = "no_flank_found";
			} else if (lHits.empty()) {
				pr.status_ = "left_flank_missing";
			} else if (rHits.empty()) {
				pr.status_ = "right_flank_missing";
			} else if ((lHits.size() > 1 || rHits.size() > 1) && !allowAmbiguousFlanks) {
				pr.status_ = njh::pasteAsStr("ambiguous_flank_hits");
			} else if (lHits.front().chrom_ != rHits.front().chrom_) {
				pr.status_ = "flanks_on_different_chromosomes";
			} else if (lHits.front().reverseSrand_ != rHits.front().reverseSrand_) {
				//one locus can only align in one orientation, so two anchors landing on
				//opposite strands means they are on two different copies. This is the only
				//thing that catches an inverted duplication, where the span can still look
				//plausible: the two elongation factor 1-alpha copies pair this way in 9 of
				//the 13 Pf references and give a span of 2.4x the record
				pr.status_ = "flanks_on_opposite_strands";
			} else {
				const auto & lh = lHits.front();
				const auto & rh = rHits.front();
				size_t spanStart = std::min(lh.start_, rh.start_);
				size_t spanEnd = std::max(lh.end_, rh.end_);
				uint32_t spanLen = spanEnd - spanStart;
				//the anchors may sit outside the record, so the span legitimately grows
				//by the two outsets before any real length difference is counted
				double lowBound = minSpanMultiplier * region.getLen();
				double highBound = maxSpanMultiplier * region.getLen() + maxSpanPad
						+ leftOutsetUsed + rightOutsetUsed;
				if (spanLen < lowBound || spanLen > highBound) {
					pr.status_ = njh::pasteAsStr("span_out_of_bounds(", spanLen, ")");
				} else {
					GenomicRegion span;
					span.chrom_ = lh.chrom_;
					span.start_ = spanStart;
					span.end_ = spanEnd;
					//The anchors go to lastz already in the record's own sense, because
					//genFlankRegion inherits the record's strand. So a reversed flank hit means
					//the genomic span runs opposite to the sense sequence and has to be reverse
					//complemented to read as the gene does, and a forward hit means it already
					//reads that way. The record's strand is therefore already accounted for and
					//must not be applied a second time: doing that wrote every minus strand
					//record in genomic forward orientation, antisense to the gene, and labelled
					//it '+'.
					span.reverseSrand_ = lh.reverseSrand_;
					span.uid_ = genome;
					pr.span_ = span;
					pr.usable_ = true;
					pr.status_ = "ok";
				}
			}
			pairedByGenome[genome] = pr;
		}

		//An outset puts the anchors outside the record, so every span, including the
		//primary genome's own, starts and ends that far outside the bed record it came
		//from. Put the edges back by aligning each anchor to the matching end of the
		//span. Nothing to do when neither side moved: the anchors are then the first
		//and last bases of the record and the span already lands on its edges.
		if (trimToBedEdges && (leftOutsetUsed > 0 || rightOutsetUsed > 0)) {
			//The trim reference is not the search anchor: an outset larger than the anchor
			//leaves the record edge outside it entirely, with nothing to align the edge to.
			//So build a window per side that reaches from the outer end of the anchor to
			//trimMargin bases inside the record, which always straddles the edge whatever
			//the outset. genFlankRegion already clamps to the chromosome, so read the edge
			//index back off the window's own coordinates rather than assuming the outset.
			uint32_t trimMargin = std::min<uint32_t>(flankLength, region.getLen());
			auto leftTrimReg = genFlankRegion(region, leftOutsetUsed + trimMargin,
					leftOutsetUsed, true, chromLength);
			auto rightTrimReg = genFlankRegion(region, rightOutsetUsed + trimMargin,
					rightOutsetUsed, false, chromLength);
			//genFlankRegion works in genomic left/right while the span reads in the
			//record's direction, so on a minus strand record the genomically right window
			//is the one holding the 5' edge
			const auto & fivePrimeReg = region.reverseSrand_ ? rightTrimReg : leftTrimReg;
			const auto & threePrimeReg = region.reverseSrand_ ? leftTrimReg : rightTrimReg;
			//index of the edge base within each window, read in record sense
			uint32_t fivePrimeEdge = region.reverseSrand_
					? static_cast<uint32_t>(fivePrimeReg.end_ - region.end_)
					: static_cast<uint32_t>(region.start_ - fivePrimeReg.start_);
			uint32_t threePrimeEdge = region.reverseSrand_
					? static_cast<uint32_t>(threePrimeReg.end_ - 1 - region.start_)
					: static_cast<uint32_t>(region.end_ - 1 - threePrimeReg.start_);
			//only a side that was actually pushed outwards has anything to give back
			const bool trimFive = (region.reverseSrand_ ? rightOutsetUsed : leftOutsetUsed) > 0;
			const bool trimThree = (region.reverseSrand_ ? leftOutsetUsed : rightOutsetUsed) > 0;
			auto fivePrimeRef = gMapper->extractGenomeRegion(primaryGenome, fivePrimeReg);
			auto threePrimeRef = gMapper->extractGenomeRegion(primaryGenome, threePrimeReg);
			for (auto & genomePair : pairedByGenome) {
				auto & pr = genomePair.second;
				if (!pr.usable_) {
					continue;
				}
				auto spanSeq = gMapper->extractGenomeRegion(genomePair.first, pr.span_);
				bool fiveOk = true;
				bool threeOk = true;
				uint32_t fiveTrim = 0;
				uint32_t threeTrim = 0;
				if (trimFive) {
					fiveTrim = findEdgeTrim(fivePrimeRef, spanSeq, fivePrimeEdge, true,
							alignerObj, fiveOk, setUp.pars_.verbose_);
				}
				if (trimThree) {
					threeTrim = findEdgeTrim(threePrimeRef, spanSeq, threePrimeEdge, false,
							alignerObj, threeOk, setUp.pars_.verbose_);
				}
				if (setUp.pars_.verbose_) {
					std::cerr << "  " << rUid << " " << genomePair.first << " trim 5p="
							<< fiveTrim << " 3p=" << threeTrim << std::endl;
				}
				if (fiveTrim + threeTrim >= pr.span_.getLen()) {
					//trimming everything away means the alignment is not telling the truth
					pr.status_ = "ok_trim_overran";
					continue;
				}
				if (pr.span_.reverseSrand_) {
					pr.span_.start_ += threeTrim;
					pr.span_.end_ -= fiveTrim;
				} else {
					pr.span_.start_ += fiveTrim;
					pr.span_.end_ -= threeTrim;
				}
				pr.fivePrimeTrim_ = fiveTrim;
				pr.threePrimeTrim_ = threeTrim;
				if (!fiveOk || !threeOk) {
					//the span is still usable, it just keeps the outset on the side that
					//could not be placed
					pr.status_ = njh::pasteAsStr("ok_untrimmed_",
							(!fiveOk && !threeOk) ? "both" : (!fiveOk ? "5p" : "3p"));
				}
			}
		}

		//One bed per region rather than one per genome per region: with 13 genomes the
		//old layout left 26 bed files in every region directory. The genome goes in the
		//first extra field as [genome=X] so rows can be split back out, and different
		//assemblies mostly use different chrom names anyway. Any gff annotation lands in
		//the field after it, since intersectBedLocsWtihGffRecords appends one field of
		//its own.
		auto annotateWithGff = [&gMapper](std::vector<Bed6RecordCore> & beds,
				const std::string & genome) {
			if (beds.empty()) {
				return;
			}
			const auto & gffFnp = gMapper->genomes_.at(genome)->gffFnp_;
			if (gffFnp.empty()) {
				return;
			}
			intersectBedLocsWtihGffRecordsPars intersectPars(gffFnp);
			intersectPars.selectFeatures_ = gMapper->pars_.gffIntersectPars_.selectFeatures_;
			intersectPars.extraAttributes_ = gMapper->pars_.gffIntersectPars_.extraAttributes_;
			intersectBedLocsWtihGffRecords(beds, intersectPars);
		};

		std::vector<seqInfo> refSeqs;
		std::vector<Bed6RecordCore> allRegionBeds;
		std::vector<Bed6RecordCore> allFlankBeds;
		std::unordered_map<std::string, std::string> gffAnnotByGenome;
		for (const auto & genome : genomeNames) {
			const auto & pr = pairedByGenome.at(genome);
			{
				//each genome is annotated against its own gff, so the beds are gathered and
				//intersected a genome at a time before being pooled
				std::vector<Bed6RecordCore> flankBedsForGenome;
				auto addFlankHits = [&flankBedsForGenome,&genome,&regionName](
						const std::vector<GenomicRegion> & hits, const std::string & side) {
					for (const auto & h : hits) {
						auto b = h.genBedRecordCore();
						b.name_ = regionName + "_" + side + "Flank";
						MetaDataInName bedMeta;
						bedMeta.addMeta("genome", genome);
						bedMeta.addMeta("side", side);
						b.extraFields_.emplace_back(bedMeta.createMetaName());
						flankBedsForGenome.emplace_back(b);
					}
				};
				if (njh::in(genome, leftRes)) {
					addFlankHits(leftRes.at(genome).all_individual_regions, "left");
				}
				if (njh::in(genome, rightRes)) {
					addFlankHits(rightRes.at(genome).all_individual_regions, "right");
				}
				annotateWithGff(flankBedsForGenome, genome);
				addOtherVec(allFlankBeds, flankBedsForGenome);
			}
			if (!pr.usable_) {
				continue;
			}
			{
				auto b = pr.span_.genBedRecordCore();
				b.name_ = regionName;
				MetaDataInName bedMeta;
				bedMeta.addMeta("genome", genome);
				b.extraFields_.emplace_back(bedMeta.createMetaName());
				std::vector<Bed6RecordCore> regionBedForGenome { b };
				annotateWithGff(regionBedForGenome, genome);
				//the gff put its hits in the field after the genome meta, carry it to the
				//info table so the annotation is readable without opening the bed
				if (regionBedForGenome.front().extraFields_.size() > 1) {
					gffAnnotByGenome[genome] = regionBedForGenome.front().extraFields_[1];
				}
				addOtherVec(allRegionBeds, regionBedForGenome);
			}
			auto extracted = gMapper->extractGenomeRegion(genome, pr.span_);
			extracted.name_ = genome;
			if (!extractPars.shortNames) {
				MetaDataInName refMeta;
				refMeta.addMeta("genome", genome);
				refMeta.addMeta("chrom", pr.span_.chrom_);
				refMeta.addMeta("start", pr.span_.start_);
				refMeta.addMeta("end", pr.span_.end_);
				refMeta.addMeta("strand", (pr.span_.reverseSrand_ ? '-' : '+'));
				refMeta.addMeta("length", pr.span_.getLen());
				extracted.name_ += " " + refMeta.createMetaName();
			}
			refSeqs.emplace_back(extracted);
		}

		{
			OutputStream regionBedOut(njh::files::make_path(bedsDir, "regions.bed"));
			for (const auto & b : allRegionBeds) {
				regionBedOut << b.toDelimStrWithExtra() << std::endl;
			}
		}
		{
			OutputStream flankBedOut(njh::files::make_path(bedsDir, "flanks.bed"));
			for (const auto & b : allFlankBeds) {
				flankBedOut << b.toDelimStrWithExtra() << std::endl;
			}
		}

		if (writeOutAllSeqsFile) {
			SeqOutput::write(refSeqs,
					SeqIOOptions::genFastaOut(njh::files::make_path(regionDir, "allRefsSeparated.fasta")));
		}
		//collapse identical extractions the way the other extractors do
		std::vector<seqInfo> collapsed;
		for (const auto & seq : refSeqs) {
			bool found = false;
			for (auto & already : collapsed) {
				if (already.seq_ == seq.seq_) {
					already.name_ += "::" + seq.name_;
					found = true;
					break;
				}
			}
			if (!found) {
				collapsed.emplace_back(seq);
			}
		}
		SeqOutput::write(collapsed,
				SeqIOOptions::genFastaOut(njh::files::make_path(regionDir, "allRefs.fasta")));

		{
			OutputStream infoOut(njh::files::make_path(bedsDir, "flankExtractionInfo.tab.txt"));
			//the trim columns go after status so that status stays column 15 for anything
			//already reading these tables
			infoOut << "genome\tregion\tflankLength\tleftOutset\trightOutset\tleftFlankHits\trightFlankHits\tchrom\tstart\tend\tstrand\tspanLength\tprimaryLength\tspanRatio\tstatus\tfivePrimeTrim\tthreePrimeTrim\tgffAnnotation"
					<< std::endl;
			for (const auto & genome : genomeNames) {
				const auto & pr = pairedByGenome.at(genome);
				infoOut << genome << "\t" << rUid << "\t" << flankLength << "\t" << leftOutsetUsed << "\t"
						<< rightOutsetUsed << "\t" << pr.leftHits_ << "\t" << pr.rightHits_;
				if (pr.usable_) {
					infoOut << "\t" << pr.span_.chrom_ << "\t" << pr.span_.start_ << "\t"
							<< pr.span_.end_ << "\t" << (pr.span_.reverseSrand_ ? '-' : '+') << "\t"
							<< pr.span_.getLen() << "\t" << region.getLen() << "\t"
							<< static_cast<double>(pr.span_.getLen()) / region.getLen();
				} else {
					infoOut << "\t*\t*\t*\t*\t*\t" << region.getLen() << "\t*";
				}
				infoOut << "\t" << pr.status_;
				if (pr.usable_) {
					infoOut << "\t" << pr.fivePrimeTrim_ << "\t" << pr.threePrimeTrim_;
				} else {
					infoOut << "\t*\t*";
				}
				infoOut << "\t"
						<< (njh::in(genome, gffAnnotByGenome) ? gffAnnotByGenome.at(genome) : "*");
				infoOut << std::endl;
			}
		}

		if (!keepRefAlignments) {
			njh::files::rmDirForce(alignsDir);
		}
		} catch (const std::exception & e) {
			std::cerr << "failed on " << rUid << ", continuing with the rest: " << e.what()
					<< std::endl;
			skippedRegions.emplace_back(njh::pasteAsStr(rUid, "\tfailed"));
		}
	}

	if (!skippedRegions.empty()) {
		OutputStream skippedOut(njh::files::make_path(outputDir, "skippedRegions.tab.txt"));
		skippedOut << "region\treason" << std::endl;
		for (const auto & skipped : skippedRegions) {
			skippedOut << skipped << std::endl;
		}
		std::cerr << skippedRegions.size() << " of " << regions.size()
				<< " regions were skipped, see "
				<< njh::files::make_path(outputDir, "skippedRegions.tab.txt") << std::endl;
	}

	if (combineByGenome) {
		njh::files::MkdirPar combinedByGenomeMkPars(njh::files::make_path(outputDir, "combinedByGenome"));
		combinedByGenomeMkPars.overWriteDir_ = overWriteDirs;
		njh::files::makeDir(combinedByGenomeMkPars);
		for (const auto & genome : genomeNames) {
			std::vector<GenomicRegion> rawRegionsForGenome;
			for (const auto & region : regions) {
				//regions.bed now holds every genome, so pick this one out by its genome meta
				auto bedFnp = njh::files::make_path(outputDir,
						region.createUidFromCoordsStrand(), "beds", "regions.bed");
				if (!bfs::exists(bedFnp)) {
					continue;
				}
				std::vector<std::shared_ptr<Bed6RecordCore>> bedsForGenome;
				for (const auto & b : getBeds(bedFnp)) {
					if (b->extraFields_.empty()
							|| !MetaDataInName::nameHasMetaData(b->extraFields_.front())) {
						continue;
					}
					MetaDataInName bedMeta(b->extraFields_.front());
					if (bedMeta.containsMeta("genome") && genome == bedMeta.getMeta("genome")) {
						bedsForGenome.emplace_back(b);
					}
				}
				addOtherVec(rawRegionsForGenome, bedPtrsToGenomicRegs(bedsForGenome));
			}
			auto regionsForGenome = mergeRegionsStrandAware(rawRegionsForGenome);
			OutputStream bedOut(njh::files::make_path(combinedByGenomeMkPars.dirName_, genome + ".bed"));
			for (const auto & reg : regionsForGenome) {
				bedOut << reg.genBedRecordCore().toDelimStrWithExtra() << std::endl;
			}
		}
	}

	return 0;
}

} // namespace njhseq
