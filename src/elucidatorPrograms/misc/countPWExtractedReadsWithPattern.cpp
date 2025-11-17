//
// Created by Nicholas Hathaway on 7/31/23.
//
#include <njhseq/GenomeUtils.h>

#include "miscRunner.hpp"
#include "elucidator/simulation.h"
#include "elucidator/BamToolsUtils.h"
#include "elucidator/objects/dataContainers/graphs/ContigsCompareGraph.hpp"
#include "elucidator/BioRecordsUtils/BedUtility.hpp"


#include <njhseq/PopulationGenetics.h>
#include <njhseq/objects/seqContainers/CollapsedHaps.hpp>
#include <njhseq/objects/Gene/TranslatorByAlignment.hpp>
#include <boost/filesystem.hpp>
#include <PathWeaver/objects/bam/RegionInvestigatorInBam.hpp>


namespace njhseq {

int miscRunner::countPWExtractedReadsWithPattern(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path inputDirectory = "./";
	bfs::path bedFnp = "";
	std::set<std::string> samples;
	std::string pat;
	std::set<std::string> seqPats;
	uint32_t numThreads = 1;
	uint32_t minPatPerRead = 1;
	uint32_t minReadCounts = 2;
	OutOptions outOpts;

	VecStr filesToInvestigation = {"extracted.fastq",
																 "extracted_R1.fastq",
																 "extracted_R2.fastq",
																 "filteredPairs_extracted_R1.fastq",
																 "filteredPairs_extracted_R2.fastq",
																 "filteredSingles_extracted.fastq",
																 "thrownAwayMate_extracted.fastq"};

	seqSetUp setUp(inputCommands);

	setUp.processDebug();
	setUp.processVerbose();

	setUp.setOption(bedFnp, "--bedFnp", "regions", true);
	setUp.setOption(pat, "--pat", "file pattern", true);
	setUp.setOption(seqPats, "--seqPat", "sequence patterns to count", true);
	setUp.setOption(inputDirectory, "--inputDirectory", "Input Directory to search");
	setUp.setOption(samples, "--samples", "Process input from only these samples");
	setUp.setOption(minReadCounts, "--minReadCounts", "min Read Counts to count a file");
	setUp.setOption(minPatPerRead, "--minPatPerRead", "min count of patterns per Read Counts to count a read");
	setUp.setOption(filesToInvestigation, "--filesToInvestigation", "files To Investigation");
	setUp.setOption(numThreads, "--numThreads", "number of threads");

	setUp.processWritingOptions(outOpts);

	setUp.finishSetUp(std::cout);
//	setUp.startARunLog(setUp.pars_.directoryName_);


	std::vector<bfs::path> directories;
	if (samples.empty()) {
		auto allFiles = njh::files::listAllFiles(inputDirectory, false, {std::regex { ".*" + pat + "$" } });
		for (const auto &f : allFiles) {
			if (f.second) {
				directories.emplace_back(f.first);
			}
		}
	} else {
		for (const auto &samp : samples) {
			directories.emplace_back(njh::files::make_path(inputDirectory, samp + pat));
		}
	}

	if(directories.empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error no directories found in " << inputDirectory << " ending with " << pat << "\n";
		throw std::runtime_error{ss.str()};
	}

	auto beds = getBeds(bedFnp);

	OutputStream out(outOpts);
	std::mutex outMut;
	out << "sample\tregion\tfile\tseqPat\tcount" << std::endl;

	njh::concurrent::LockableVec<bfs::path> dirQueue(directories);
	std::function<void()> countSample = [&dirQueue,&beds,&filesToInvestigation,&out,&outMut, &minReadCounts,
																			 &seqPats, &setUp, &pat, &minPatPerRead](){
		bfs::path dir;
		while(dirQueue.getVal(dir)){
			std::string sample = std::regex_replace(dir.filename().string(), std::regex(pat), "");
			for(const auto & bed : beds){
				for(const auto & f : filesToInvestigation){
					bfs::path inputFnp = njh::files::make_path(dir, bed->name_, sample + "_extraction", f);
					if(setUp.pars_.debug_){
						std::cout << "dir: " << dir << ", " << "region: " << bed->name_ << ", " << "f: " << f << std::endl;
						std::cout << "\t" << inputFnp << std::endl;
					}
					if(bfs::exists(inputFnp)){
						for(const auto & seqPat : seqPats){
							std::regex seqPatReg{seqPat};
							seqInfo seq;
							auto inputFnpOpts = SeqIOOptions::genFastqIn(inputFnp);
							SeqInput reader(inputFnpOpts);
							reader.openIn();
							std::stringstream ss;
							uint32_t readCounts = 0;
							while(reader.readNextRead(seq)){
								std::ptrdiff_t const match_count(std::distance(
												std::sregex_iterator(seq.seq_.begin(), seq.seq_.end(), seqPatReg),
												std::sregex_iterator()));
								if(match_count >= minPatPerRead){
									++readCounts;
								}
							}
							if(readCounts >= minReadCounts){
								std::lock_guard<std::mutex> lock(outMut);
								out << sample
										<< "\t" << bed->name_
										<< "\t" << f
										<< "\t" << seqPat
										<< "\t" << readCounts << std::endl;
							}
						}
					}
				}
			}
		}
	};

	njh::concurrent::runVoidFunctionThreaded(countSample, numThreads);

	return 0;

}



std::string getPossibleSampleNameFromFnp(const bfs::path & fnp){
	std::string bName = bfs::basename(fnp);
	if(std::string::npos != bName.find('.')){
		bName = bName.substr(0, bName.find('.'));
	}
	return bName;
}


int miscRunner::scanningRegionsForReadsWithPattern(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bedFnp = "";
	std::set<std::string> seqPats;
	bool write_out_seqs = false;
	uint32_t numThreads = 1;
	uint32_t minPatPerRead = 1;
	uint32_t minReadCounts = 2;
	BamRegionInvestigator::BamRegionInvestigatorPars brInvestPars;
	bfs::path genomeFnp;
	double percInRegion = 0.5;
	bool check_reverse_strand = false;
	BamRegionInvestigator::BamCountSpecficRegionsPars spanningReadsPar;

	seqSetUp setUp(inputCommands);

	setUp.processDebug();
	setUp.processVerbose();

	setUp.setOption(brInvestPars.mapQualityCutOffForMultiMap_, "--mapQualityCutOffForMultiMapForCovInfo", "map Quality Cut Off For Coverage");
	setUp.setOption(brInvestPars.mapQualityCutOff_, "--mapQualityCutOffForCov", "map Quality Cut Off For Coverage");
	setUp.setOption(brInvestPars.countDups_, "--countDups", "count Dups");
	setUp.setOption(percInRegion, "--percInRegion", "read perc In Region in order to be extracted");
	setUp.setOption(check_reverse_strand, "--check_reverse_strand", "check_reverse_strand");
	setUp.setOption(write_out_seqs, "--write_out_seqs", "write_out_seqs");


	setUp.setOption(bedFnp, "--bedFnp", "regions", true);
	setUp.setOption(seqPats, "--seqPat", "sequence patterns to count", true);
	setUp.setOption(genomeFnp, "--genomeFnp", "genome fnp", true);


	setUp.setOption(minReadCounts, "--minReadCounts", "min Read Counts to count a file");
	setUp.setOption(minPatPerRead, "--minPatPerRead", "min count of patterns per Read Counts to count a read");
	setUp.setOption(numThreads, "--numThreads", "number of threads");
	brInvestPars.numThreads_ = numThreads;
	spanningReadsPar.numThreads = numThreads;
	spanningReadsPar.countDuplicates = brInvestPars.countDups_;

	setUp.processReadInNames({"--bam"}, true);
	setUp.processDirectoryOutputName(true);
	setUp.finishSetUp(std::cout);

	setUp.startARunLog(setUp.pars_.directoryName_);

	auto sampName = getPossibleSampleNameFromFnp(setUp.pars_.ioOptions_.firstName_);
	auto beds = getBeds(bedFnp);
	auto inputRegions = bedPtrsToGenomicRegs(beds);
	sortGRegionsByStart(inputRegions);
	// get regex patterns
	std::vector<std::pair<std::string, std::regex>> seqPatRegVec;
	for(const auto & seqPat : seqPats) {
		if (setUp.pars_.debug_) {
			std::cout << "seqPat: " << seqPat << std::endl;
		}
		std::regex seqPatReg{seqPat};
		seqPatRegVec.emplace_back(std::make_pair(seqPat, seqPatReg));
	}

	//get coverage and spanning read info
	auto genomeFnp2bit = genomeFnp;
	genomeFnp2bit.replace_extension(".2bit");
	BamRegionInvestigator brInvestor(brInvestPars);
	auto regInfos = brInvestor.getCoverageAndFullSpanningReads(setUp.pars_.ioOptions_.firstName_, inputRegions, spanningReadsPar, genomeFnp2bit);

	brInvestor.writeBasicInfo(regInfos, sampName, OutOptions(njh::files::make_path(setUp.pars_.directoryName_, "perBaseCoveragePerRegion.bed")));

	//scan across regions
	njh::concurrent::LockableQueue<GenomicRegion> inputRegionsQueue(inputRegions);
	concurrent::BamReaderPool bamPool(setUp.pars_.ioOptions_.firstName_, numThreads);
	bamPool.openBamFile();
	std::mutex counts_mut;
	OutputStream counts_out(njh::files::make_path(setUp.pars_.directoryName_, "region_counts.tsv.gz"));
	counts_out << "#chrom\tstart\tend\tname\tlength\tstrand\tsample\tseq_pat\tseqs_with_pat\tseqs_with_pat_freq";
	if (check_reverse_strand) {
		counts_out << "\tseqs_with_pat_neg_strand\tseqs_with_pat_freq_neg_strand";
	}
	counts_out << "\ttotal_seqs" << std::endl;
	std::function<void()> countPatternPerRegion = [&inputRegionsQueue,
		&minReadCounts,
		&write_out_seqs,
		&counts_mut,&counts_out,
		&setUp, &minPatPerRead,
		&seqPatRegVec,
		&genomeFnp2bit,
		&bamPool, percInRegion, check_reverse_strand, sampName]() {

		GenomicRegion region;
		std::stringstream current_out;
		while (inputRegionsQueue.getVal(region)) {
			BamExtractor bExtractor(setUp.pars_.verbose_);
			bExtractor.debug_ = setUp.pars_.debug_;
			const bfs::path bamFnp = setUp.pars_.ioOptions_.firstName_;
			TwoBit::TwoBitFile tReader(genomeFnp2bit);
			auto bamReader = bamPool.popReader();
			uint32_t total_reads = 0;
			std::unordered_map<std::string, uint32_t> total_reads_with_pats;
			std::unordered_map<std::string, uint32_t> total_reads_with_pats_rev_comp;

			//extract reads
			auto extracted_reads = bExtractor.extractReadsFromBamRegionAlns(*bamReader, region, percInRegion);

			std::unique_ptr<SeqOutput> writer_r1;
			std::unique_ptr<SeqOutput> writer_r2;
			std::unique_ptr<SeqOutput> writer_singles;
			if (write_out_seqs) {
				writer_r1 = std::make_unique<SeqOutput>(SeqIOOptions::genFastqOutGz(njh::files::make_path(setUp.pars_.directoryName_, region.uid_ + "_pairs_R1")));
				writer_r2 = std::make_unique<SeqOutput>(SeqIOOptions::genFastqOutGz(njh::files::make_path(setUp.pars_.directoryName_, region.uid_ + "_pairs_R2")));
				writer_singles = std::make_unique<SeqOutput>(SeqIOOptions::genFastqOutGz(njh::files::make_path(setUp.pars_.directoryName_, region.uid_ + "_singles")));
			}



			//regular pairs
			for (const auto &aln: extracted_reads.pairs_) {
				++total_reads;
				auto r1_seq = bamAlnToSeqInfo(aln.first, true);
				auto r2_seq = bamAlnToSeqInfo(aln.second, true);
				for (const auto &seqPat: seqPatRegVec) {
					std::ptrdiff_t const match_count_r1(std::distance(
						std::sregex_iterator(r1_seq.seq_.begin(), r1_seq.seq_.end(), seqPat.second),
						std::sregex_iterator()));
					std::ptrdiff_t const match_count_r2(std::distance(
						std::sregex_iterator(r2_seq.seq_.begin(), r2_seq.seq_.end(), seqPat.second),
						std::sregex_iterator()));
					if (match_count_r1 >= minPatPerRead || match_count_r2 >= minPatPerRead) {
						++total_reads_with_pats[seqPat.first];
						if (write_out_seqs) {
							auto r1_seq_out = r1_seq;
							auto r2_seq_out = r2_seq;
							MetaDataInName r1_meta;
							r1_meta.addMeta("seq_pat", seqPat.first);
							r1_meta.addMeta("seq_pat_count", match_count_r1);
							r1_seq_out.name_ += r1_meta.createMetaName();

							MetaDataInName r2_meta;
							r2_meta.addMeta("seq_pat", seqPat.first);
							r2_meta.addMeta("seq_pat_count", match_count_r2);
							r2_seq_out.name_ += r2_meta.createMetaName();

							writer_r1->openWrite(r1_seq_out);
							writer_r2->openWrite(r2_seq_out);
						}
					}
				}
				if (check_reverse_strand) {
					r1_seq.reverseComplementRead(true);
					r2_seq.reverseComplementRead(true);
					for (const auto &seqPat: seqPatRegVec) {
						std::ptrdiff_t const match_count_r1(std::distance(
							std::sregex_iterator(r1_seq.seq_.begin(), r1_seq.seq_.end(), seqPat.second),
							std::sregex_iterator()));
						std::ptrdiff_t const match_count_r2(std::distance(
							std::sregex_iterator(r2_seq.seq_.begin(), r2_seq.seq_.end(), seqPat.second),
							std::sregex_iterator()));
						if (match_count_r1 >= minPatPerRead || match_count_r2 >= minPatPerRead) {
							++total_reads_with_pats_rev_comp[seqPat.first];
							if (write_out_seqs) {
								auto r1_seq_out = r1_seq;
								auto r2_seq_out = r2_seq;
								MetaDataInName r1_meta;
								r1_meta.addMeta("seq_pat", seqPat.first);
								r1_meta.addMeta("seq_pat_count", match_count_r1);
								r1_seq_out.name_ += r1_meta.createMetaName();

								MetaDataInName r2_meta;
								r2_meta.addMeta("seq_pat", seqPat.first);
								r2_meta.addMeta("seq_pat_count", match_count_r2);
								r2_seq_out.name_ += r2_meta.createMetaName();

								writer_r1->openWrite(r1_seq_out);
								writer_r2->openWrite(r2_seq_out);
							}
						}
					}
				}
			}
			//mate unmapped
			for (const auto & aln : extracted_reads.pairsMateUnmapped_) {
				++total_reads;
				seqInfo seq;
				if (aln.first.IsMapped()) {
					seq = bamAlnToSeqInfo(aln.first, true);
				} else {
					seq = bamAlnToSeqInfo(aln.second, true);
				}
				for (const auto &seqPat: seqPatRegVec) {
					std::ptrdiff_t const match_count(std::distance(
						std::sregex_iterator(seq.seq_.begin(), seq.seq_.end(), seqPat.second),
						std::sregex_iterator()));
					if (match_count >= minPatPerRead) {
						++total_reads_with_pats[seqPat.first];
						if (write_out_seqs) {
							auto seq_out = seq;
							MetaDataInName meta;
							meta.addMeta("seq_pat", seqPat.first);
							meta.addMeta("seq_pat_count", match_count);
							seq_out.name_ += meta.createMetaName();

							writer_singles->openWrite(seq_out);
						}
					}
				}
				if (check_reverse_strand) {
					seq.reverseComplementRead();
					for (const auto &seqPat: seqPatRegVec) {
						std::ptrdiff_t const match_count(std::distance(
							std::sregex_iterator(seq.seq_.begin(), seq.seq_.end(), seqPat.second),
							std::sregex_iterator()));
						if (match_count >= minPatPerRead) {
							++total_reads_with_pats_rev_comp[seqPat.first];
							if (write_out_seqs) {
								auto seq_out = seq;
								MetaDataInName meta;
								meta.addMeta("seq_pat", seqPat.first);
								meta.addMeta("seq_pat_count", match_count);
								seq_out.name_ += meta.createMetaName();

								writer_singles->openWrite(seq_out);
							}
						}
					}
				}
			}

			// singlets
			for (const auto & aln : extracted_reads.singlets_) {
				++total_reads;
				seqInfo seq = bamAlnToSeqInfo(aln, true);
				for (const auto &seqPat: seqPatRegVec) {
					std::ptrdiff_t const match_count(std::distance(
						std::sregex_iterator(seq.seq_.begin(), seq.seq_.end(), seqPat.second),
						std::sregex_iterator()));
					if (match_count >= minPatPerRead) {
						++total_reads_with_pats[seqPat.first];
						if (write_out_seqs) {
							auto seq_out = seq;
							MetaDataInName meta;
							meta.addMeta("seq_pat", seqPat.first);
							meta.addMeta("seq_pat_count", match_count);
							seq_out.name_ += meta.createMetaName();

							writer_singles->openWrite(seq_out);
						}
					}
				}
				if (check_reverse_strand) {
					seq.reverseComplementRead();
					for (const auto &seqPat: seqPatRegVec) {
						std::ptrdiff_t const match_count(std::distance(
							std::sregex_iterator(seq.seq_.begin(), seq.seq_.end(), seqPat.second),
							std::sregex_iterator()));
						if (match_count >= minPatPerRead) {
							++total_reads_with_pats_rev_comp[seqPat.first];
							if (write_out_seqs) {
								auto seq_out = seq;
								MetaDataInName meta;
								meta.addMeta("seq_pat", seqPat.first);
								meta.addMeta("seq_pat_count", match_count);
								seq_out.name_ += meta.createMetaName();

								writer_singles->openWrite(seq_out);
							}
						}
					}
				}
			}

			for (const auto & pat_count : total_reads_with_pats) {
				current_out << region.genBedRecordCore().toDelimStr()
				<< "\t" << sampName << "\t" << pat_count.first << "\t" << pat_count.second << "\t" << static_cast<double>(pat_count.second) / total_reads ;
				if (check_reverse_strand) {
					counts_out << "\t" << total_reads_with_pats_rev_comp[pat_count.first] << "\t" << static_cast<double>(total_reads_with_pats_rev_comp[pat_count.first]) / total_reads ;
				}
				current_out << "\t" << total_reads << std::endl;
			}
		}
		{
			std::lock_guard lock(counts_mut);
			counts_out << current_out.str();
		}
	};

	njh::concurrent::runVoidFunctionThreaded(countPatternPerRegion, numThreads);

	return 0;

}

}  //namespace njhseq

