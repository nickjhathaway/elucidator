//
// Created by Nicholas Hathaway on 5/6/23.
//

#include "kmerExp.hpp"
#include <njhseq/objects/kmer/KmerUtils.hpp>
#include <njhseq/objects/kmer.h>

#include "elucidator/simulation.h"

#include "elucidator/objects/MiscUtility/GenomeSeqSearch.hpp"

#include <njhseq/objects/dataContainers/tables/TableReader.hpp>
#include <njhseq/IO/SeqIO/MultiSeqOutCache.hpp>
#include <njhseq/objects/BioDataObject/reading.hpp>
#include <njhseq/objects/BioDataObject/BioRecordsUtils/BedUtility.hpp>

namespace njhseq {




int kmerExpRunner::findUniqKmersFromGenomeSubRegions(const njh::progutils::CmdArgs & inputCommands){
	KmerGatherer::KmerGathererPars countPars;
	bfs::path genomeFnp;
	bfs::path bedFnp;
	std::string regionName;
	bool getReverseCompOfInputRegions = false;
	bool doNotReverseCompOfGenomeRegions = false;
	OutOptions outOpts;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.description_ = "Get the unique kmers that appear within the bed region compared to the rest of the genome";

	setUp.setOption(genomeFnp, "--genomeFnp", "genome file to extract from, a .2bit file needs to exist for the supplied genome", true);
	setUp.setOption(bedFnp, "--bedFnp", "sub regions to extract to compare to the rest of the genome", true);
	setUp.setOption(regionName, "--regionName", "region name to give to the unique kmer set", true);

	setUp.setOption(getReverseCompOfInputRegions, "--getReverseCompOfInputRegions", "get Reverse Comp Of Input Regions");
	setUp.setOption(doNotReverseCompOfGenomeRegions, "--doNotReverseCompOfGenomeRegions", "do Not Reverse Comp Of Genome Regions");

	setUp.processWritingOptions(outOpts);
	countPars.setOptions(setUp);
	setUp.finishSetUp(std::cout);

	KmerGatherer kGather(countPars);

	//set up 2bit file
	bfs::path genome2bitFnp = bfs::path(genomeFnp).replace_extension("").string() + ".2bit";
	njh::files::checkExistenceThrow(genome2bitFnp, __PRETTY_FUNCTION__ );
	TwoBit::TwoBitFile tReader(genome2bitFnp);
	auto chromLengths = tReader.getSeqLens();

	//get input regions
	auto bedRegions = getBeds(bedFnp);
	BedUtility::coordSort(bedRegions);

	////check in put
	////// check for overlaping
	VecStr overlappingRegionWarnings;
	if(bedRegions.size() > 1){
		for(const auto pos : iter::range(1UL, bedRegions.size())){
			if(bedRegions[pos]->overlaps(*bedRegions[pos - 1], 1)){
				overlappingRegionWarnings.emplace_back(njh::pasteAsStr(bedRegions[pos]->name_, " overlaps ", bedRegions[pos - 1]->name_));
			}
		}
	}
	if(!overlappingRegionWarnings.empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "found the following overlaps" << "\n";
		ss << njh::conToStr(overlappingRegionWarnings, "\n") << "\n";
		throw std::runtime_error{ss.str()};
	}

	////// chroms
	VecStr chromsDoNotExist;
	for(const auto & b : bedRegions){
		if(!njh::in(b->chrom_, chromLengths)){
			chromsDoNotExist.emplace_back(b->chrom_);
		}
	}
	if(!chromsDoNotExist.empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "the following chroms could not be found within " << genome2bitFnp << "\n";
		ss << njh::conToStr(chromsDoNotExist, ", ") << "\n";
		throw std::runtime_error{ss.str()};
	}
	std::unordered_map<std::string, std::vector<std::shared_ptr<Bed6RecordCore>>> bedsByChrom;
	for(const auto & bedRegion : bedRegions){
		bedsByChrom[bedRegion->chrom_].emplace_back(bedRegion);
	}

	//// subtract region to create other regions
	std::vector<Bed6RecordCore> inbetweenRegions;
	for (const auto &chrom: chromLengths) {
		if (njh::in(chrom.first, bedsByChrom)) {
			//add front
			if (0 != bedsByChrom[chrom.first].front()->chromStart_) {
				inbetweenRegions.emplace_back(chrom.first, 0, bedsByChrom[chrom.first].front()->chromStart_,
																			njh::pasteAsStr(chrom.first, "-", 0, "-",
																											bedsByChrom[chrom.first].front()->chromStart_),
																			bedsByChrom[chrom.first].front()->chromStart_, '+');
			}
			if(bedsByChrom[chrom.first].size() > 1){
				for(const auto pos : iter::range(bedsByChrom[chrom.first].size() - 1)){
					auto start = bedsByChrom[chrom.first][pos]->chromEnd_;
					auto end = bedsByChrom[chrom.first][pos + 1]->chromStart_;
					if(start != end){
						inbetweenRegions.emplace_back(chrom.first, start, end,
																					njh::pasteAsStr(chrom.first, "-", start, "-",
																													end),
																					end - start, '+');
					}
				}
			}
			//add back
			if (chrom.second != bedsByChrom[chrom.first].back()->chromEnd_) {
				inbetweenRegions.emplace_back(chrom.first, bedsByChrom[chrom.first].back()->chromEnd_, chrom.second,
																			njh::pasteAsStr(chrom.first, "-", bedsByChrom[chrom.first].back()->chromEnd_, "-",
																											chrom.second),
																			chrom.second - bedsByChrom[chrom.first].back()->chromEnd_, '+');
			}
		} else {
			//no bed regions for this chrom
			inbetweenRegions.emplace_back(chrom.first, 0, chrom.second, chrom.first, chrom.second, '+');
		}
	}

	//// output
	OutputStream out(outOpts);
	//// hasher
	SimpleKmerHash hasher;

	//// get kmers for input bed regions
	std::set<uint64_t> rawKmersPerInput;
	for(const auto & region : bedRegions){
		auto currentSeq = GenomicRegion(*region).extractSeq(tReader);

		for (uint32_t pos = 0; pos < len(currentSeq.seq_) - countPars.kmerLength_ + 1; ++pos) {
			rawKmersPerInput.emplace(
							hasher.hash(currentSeq.seq_.substr(pos, countPars.kmerLength_)));
		}
		if (getReverseCompOfInputRegions) {
			currentSeq.seq_ = seqUtil::reverseComplement(currentSeq.seq_, "DNA");
			for (uint32_t pos = 0; pos < len(currentSeq.seq_) - countPars.kmerLength_ + 1;
					 ++pos) {
				rawKmersPerInput.emplace(
								hasher.hash(currentSeq.seq_.substr(pos, countPars.kmerLength_)));
			}
		}
	}

	std::set<uint64_t> finalKmersPerInput = rawKmersPerInput;

	//// get kmers for the rest of the genome

	for(const auto & region : inbetweenRegions){
		if(setUp.pars_.verbose_){
			std::cout << region.name_ << std::endl;
			std::cout << "\t" << finalKmersPerInput.size() << std::endl;
		}
		std::set<uint64_t> kmersPerInbetween;
		auto currentSeq = GenomicRegion(region).extractSeq(tReader);

		for (uint32_t pos = 0; pos < len(currentSeq.seq_) - countPars.kmerLength_ + 1; ++pos) {
			kmersPerInbetween.emplace(
							hasher.hash(currentSeq.seq_.substr(pos, countPars.kmerLength_)));
		}
		if (!doNotReverseCompOfGenomeRegions) {
			currentSeq.seq_ = seqUtil::reverseComplement(currentSeq.seq_, "DNA");
			for (uint32_t pos = 0; pos < len(currentSeq.seq_) - countPars.kmerLength_ + 1;
					 ++pos) {
				kmersPerInbetween.emplace(
								hasher.hash(currentSeq.seq_.substr(pos, countPars.kmerLength_)));
			}
		}
		//filter kmers
		std::set<uint64_t> filterKmers;
		for(const auto & finalKmer : finalKmersPerInput){
			if(!njh::in(finalKmer, kmersPerInbetween)){
				filterKmers.emplace(finalKmer);
			}
		}
		finalKmersPerInput = filterKmers;
		if(setUp.pars_.verbose_){
			std::cout << "\t" << finalKmersPerInput.size() << std::endl << std::endl;
		}
	}

	for(const auto & finalKmer : finalKmersPerInput){
		out << regionName << "\t" << hasher.reverseHash(finalKmer) << std::endl;
	}
	return 0;
}

}  //namespace njhseq

