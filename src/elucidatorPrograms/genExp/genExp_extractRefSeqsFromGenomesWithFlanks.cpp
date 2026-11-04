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

	MultiGenomeMapper::getRefSeqsWithPrimaryGenomePars extractPars;
	extractPars.lzPars.identity = 80;
	extractPars.lzPars.coverage = 90;

	bfs::path gffDir = "";
	std::string gffExtraAttributesStr = "description";

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

	setUp.setOption(extractPars.shortNames, "--shortNames", "Short Names");
	setUp.setOption(writeOutAllSeqsFile, "--writeOutAllSeqsFile", "Write Out All Seqs File");
	setUp.setOption(combineByGenome, "--combineByGenome",
			"Combine and merge regions extracted and create a bed file for each genome");
	setUp.setOption(extractPars.extendAndTrim, "--extendAndTrim",
			"Extend the determine region and then trim back, can be helpful for when variation falls at the very ends of the sequence");
	setUp.setOption(extractPars.extendAndTrimLen, "--extendAndTrimLen",
			"When extending and trimming, use this length");
	setUp.setOption(extractPars.lzPars.identity, "--identity",
			"Identity to use for when searching with lastz");
	setUp.setOption(extractPars.lzPars.coverage, "--coverage",
			"Coverage to use for when searching with lastz; this now applies to each flank rather than the whole region, so the default is appropriate");
	setUp.setOption(extractPars.byScore, "--keepBestOnlyByScore",
			"Determine best flank hit by alignment score instead of percent identity");
	setUp.setOption(genomeMappingPars.numThreads_, "--numThreads", "Number of CPUs to utilize");
	setUp.setOption(gffDir, "--gffDir",
			"A directory with a gff for the genomes in --genomeDir, should be named GENOME.gff (for GENOME.fasta)");
	setUp.setOption(gffExtraAttributesStr, "--gffExtraAttributes",
			"Extra attributes to add to genome that has an accompany gff");
	setUp.setOption(keepRefAlignments, "--keepRefAlignments", "Keep Ref Alignments");

	setUp.setOption(bedFile, "--bed", "Bed file of the regions to extract", true);
	setUp.setOption(genomeMappingPars.genomeDir_, "--genomeDir", "Name of the genome file fnp", true);
	setUp.setOption(genomeMappingPars.primaryGenome_, "--primaryGenome",
			"Name of the primary genome, should be in --genomeDir", true);
	std::string selectedGenomesStr;
	setUp.setOption(selectedGenomesStr, "--selectedGenomes",
			"Name of the other genomes in --genomeDir to be read in, leave blank to just do all fastas");
	std::string extraGffAttributes;
	setUp.setOption(extraGffAttributes, "--extraGffAttributes", "Extra Gff Attributes");
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
					//a flank that maps reversed means the region is inverted relative to the bed
					span.reverseSrand_ =
							lh.reverseSrand_ ? !region.reverseSrand_ : region.reverseSrand_;
					span.uid_ = genome;
					pr.span_ = span;
					pr.usable_ = true;
					pr.status_ = "ok";
				}
			}
			pairedByGenome[genome] = pr;
		}

		//write out the flank beds, the paired region beds and the sequences
		std::vector<seqInfo> refSeqs;
		for (const auto & genome : genomeNames) {
			const auto & pr = pairedByGenome.at(genome);
			{
				OutputStream flankBedOut(njh::files::make_path(bedsDir, genome + "_flanks.bed"));
				if (njh::in(genome, leftRes)) {
					for (const auto & h : leftRes.at(genome).all_individual_regions) {
						auto b = h.genBedRecordCore();
						b.name_ = genome + "_leftFlank";
						flankBedOut << b.toDelimStrWithExtra() << std::endl;
					}
				}
				if (njh::in(genome, rightRes)) {
					for (const auto & h : rightRes.at(genome).all_individual_regions) {
						auto b = h.genBedRecordCore();
						b.name_ = genome + "_rightFlank";
						flankBedOut << b.toDelimStrWithExtra() << std::endl;
					}
				}
			}
			if (!pr.usable_) {
				continue;
			}
			{
				OutputStream regionBedOut(njh::files::make_path(bedsDir, genome + "_region.bed"));
				auto b = pr.span_.genBedRecordCore();
				b.name_ = genome;
				regionBedOut << b.toDelimStrWithExtra() << std::endl;
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
			infoOut << "genome\tregion\tflankLength\tleftOutset\trightOutset\tleftFlankHits\trightFlankHits\tchrom\tstart\tend\tstrand\tspanLength\tprimaryLength\tspanRatio\tstatus"
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
				infoOut << "\t" << pr.status_ << std::endl;
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
				auto bedFnp = njh::files::make_path(outputDir,
						region.createUidFromCoordsStrand(), "beds", genome + "_region.bed");
				if (bfs::exists(bedFnp)) {
					addOtherVec(rawRegionsForGenome, bedPtrsToGenomicRegs(getBeds(bedFnp)));
				}
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
