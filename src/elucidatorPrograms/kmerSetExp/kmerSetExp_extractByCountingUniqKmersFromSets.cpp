//
// Created by Nicholas Hathaway on 5/26/23.
//
#include "kmerSetExp.hpp"
#include <njhseq/IO/SeqIO/MultiSeqIO.hpp>
#include <njhseq/BamToolsUtils/bamExtractUtils.hpp>
#include "elucidator/helpers/UniqueKmerSetHelper.hpp"

namespace njhseq {

int kmerSetExpRunner::countingUniqKmersFromSetsPerRead(const njh::progutils::CmdArgs & inputCommands){
	uint32_t numThreads = 1;
	bfs::path countTable;
	bfs::path nonUniqueKmerTable;

	UniqueKmerSetHelper::CompareReadToSetPars compPars;
	compPars.kmerLengthForEntropyCalc_ = 2;
	compPars.entropyFilter_ = 1.20;
	OutOptions outOpts(bfs::path(""), ".tab.txt.gz");
	seqSetUp setUp(inputCommands);
	//default options
	setUp.processVerbose();
	setUp.processDebug();
	//read in options
	setUp.pars_.ioOptions_.revComplMate_ = true;
	bool pairedEndSet = setUp.processReadInNames(njhseq::seqSetUp::pairedReadInFormatsAvailable_, false);
	auto singlesOption = SeqIOOptions();
	setUp.processJustReadInNames(singlesOption, njhseq::seqSetUp::singleInFormatsAvailable_, !pairedEndSet);
	singlesOption.lowerCaseBases_ = setUp.pars_.ioOptions_.lowerCaseBases_;
	singlesOption.removeGaps_ = setUp.pars_.ioOptions_.removeGaps_;
	singlesOption.includeWhiteSpaceInName_ = setUp.pars_.ioOptions_.includeWhiteSpaceInName_;

	//comp results
	setUp.setOption(countTable, "--kmerTable", "kmerTable, no header, columns 1)set,2)kmer", true);
	{
		auto possibleNonUniqueTable = njh::files::prependFileBasename(countTable, "nonUniqueKmers_");
		if(exists(possibleNonUniqueTable)){
			nonUniqueKmerTable = possibleNonUniqueTable;
		}
	}
	setUp.setOption(nonUniqueKmerTable, "--nonUniqueKmerTable", "non-unique Kmer Table, 1)set,2)kmer");
	setUp.setOption(compPars.sampleName, "--sampleName", "Name to add to output file", true);
	setUp.setOption(compPars.includeRevComp, "--includeRevComp", "include Rev Comp of the input seqs");
	setUp.setOption(compPars.pairsSeparate, "--pairsSeparate", "count the pairs separately");
	setUp.setOption(compPars.kmerLengthForEntropyCalc_, "--kmerLengthForEntropyCalc", "kmer Length For Entropy Calc");
	setUp.setOption(compPars.entropyFilter_, "--entropyFilter", "entropy Filter");
	bool noInvertedKmerIndexCounting = false;
	setUp.setOption(noInvertedKmerIndexCounting, "--noInvertedKmerIndex",
			"Probe every k-mer set for every read k-mer instead of using one k-mer -> set(s) index. Counts are identical either way; the index is faster, increasingly so with more sets");
	compPars.useInvertedKmerIndex_ = !noInvertedKmerIndexCounting;
	//run options
	setUp.setOption(numThreads, "--numThreads", "numThreads");
	//write options
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);
	njh::stopWatch watch;
	watch.setLapName("initial");

	watch.startNewLap("reading in unique kmer table");
	compPars.klen = UniqueKmerSetHelper::getKmerLenFromUniqueKmerTable(countTable);
	std::unordered_map<std::string, std::unordered_set<uint64_t>> uniqueKmersPerSet = UniqueKmerSetHelper::readInUniqueKmerTablePerSet(countTable);
	{
		//the index handles shared k-mers, but whether the sets are actually exclusive changes how the
		//winner-takes-most classification behaves, so report it rather than assume it
		auto exclusivityIndex = UniqueKmerSetHelper::getInvertedIndex(uniqueKmersPerSet);
		if (nullptr != exclusivityIndex) {
			std::cout << "kmerSetExclusivity\tsets\t" << exclusivityIndex->nSets_
					<< "\ttotalKmers\t" << exclusivityIndex->totalKmers_
					<< "\tsharedKmers\t" << exclusivityIndex->nSharedKmers_
					<< "\texclusive\t" << njh::boolToStr(0 == exclusivityIndex->nSharedKmers_) << std::endl;
		}
	}
	std::string nonUniqueRegionName = "NON_UNIQUE";
	std::unordered_set<uint64_t> nonUniqueKmersPerSet;
	if (!nonUniqueKmerTable.empty()) {
		nonUniqueKmersPerSet = UniqueKmerSetHelper::readInUniqueKmerTableSetsCollapsed(nonUniqueKmerTable);
	}
	if(setUp.pars_.verbose_){
		std::cout << watch.getLapName() << "\t" << watch.timeLapFormatted() <<std::endl;
	}

	OutputStream out(outOpts);
	UniqueKmerSetHelper::CompareReadToSetRes::writeOutputHeader(out, compPars);

	std::mutex outMut;
	std::function<void()> readInComp;
	if(!setUp.pars_.ioOptions_.firstName_.empty()){
		SeqInput reader(setUp.pars_.ioOptions_);
		reader.openIn();
		readInComp = [&reader, &uniqueKmersPerSet,&compPars,&out,&outMut]() {
			SimpleKmerHash hasher;
			PairedRead pseq;
			while(reader.readNextReadLock(pseq)){
				if(compPars.pairsSeparate){
					auto compResFirstMate = UniqueKmerSetHelper::compareReadToSets(pseq.seqBase_, uniqueKmersPerSet, compPars,
																																				 hasher);
					pseq.seqBase_.name_.append("_firstMate");
					auto compResSecondMate = UniqueKmerSetHelper::compareReadToSets(pseq.mateSeqBase_, uniqueKmersPerSet,
																																					compPars,
																																					hasher);
					pseq.mateSeqBase_.name_.append("_secondMate");
					{
						std::lock_guard<std::mutex> lockGuard(outMut);
						compResFirstMate.writeOutput(out,pseq.seqBase_,uniqueKmersPerSet, compPars);
						compResSecondMate.writeOutput(out,pseq.mateSeqBase_,uniqueKmersPerSet, compPars);
					}
				}else{
					auto compRes = UniqueKmerSetHelper::compareReadToSets(pseq, uniqueKmersPerSet, compPars, hasher);
					{
						std::lock_guard<std::mutex> lockGuard(outMut);
						compRes.writeOutput(out,pseq.seqBase_,uniqueKmersPerSet, compPars);
					}
				}
			}
		};
		njh::concurrent::runVoidFunctionThreaded(readInComp, numThreads);
	}
	if (!singlesOption.firstName_.empty()) {
		SeqInput reader(singlesOption);
		reader.openIn();
		readInComp = [&reader, &uniqueKmersPerSet, &compPars, &out, &outMut]() {
			SimpleKmerHash hasher;
			seqInfo seq;
			while (reader.readNextReadLock(seq)) {
				auto compRes = UniqueKmerSetHelper::compareReadToSets(seq, uniqueKmersPerSet, compPars, hasher);
				{
					std::lock_guard<std::mutex> lockGuard(outMut);
					compRes.writeOutput(out, seq, uniqueKmersPerSet, compPars);
				}
			}
		};
		njh::concurrent::runVoidFunctionThreaded(readInComp, numThreads);
	}
	return 0;
}


/*
*						if(extractingPars.markReadsPerIteration){
							MetaDataInName meta;
							if(MetaDataInName::nameHasMetaData(seq.name_)){
								meta = MetaDataInName(seq.name_);
							}
							meta.addMeta("iteration", "initial", true);
							meta.resetMetaInName(seq.name_);
						}
 **/

namespace {

/**@brief Build one strand-bias k-mer corrector per set from the reads the final extraction will write.
 *
 * Read-only pass over the same inputs the final extraction uses, after the k-mer sets are frozen.
 * The winning orientation from compareReadToSets IS the sequencing strand relative to the set,
 * which is the axis the correction needs and which is lost once the read is written. Note the mate
 * convention: the inputs are read with revComplMate_ on, so both mates sit in one frame and share a
 * winnerRevComp; physically they were sequenced in opposite directions, so the second mate's strand
 * is the complement of the first's.
 *
 * In a multigene family strand restriction alone also admits genuine rare paralogs whose coverage
 * happens to be one sided, so the corrector is built with requireQualityDrop_ set: every correction
 * must also show the reference base degraded on the strand carrying the candidate.
 */
std::shared_ptr<const std::unordered_map<std::string, std::shared_ptr<njhseq::StrandBiasKmerCorrector>>>
buildStrandBiasCorrectorsPerSet(const SeqIOOptions & pairedOpts, const SeqIOOptions & singlesOpts,
		const std::unordered_map<std::string, std::unordered_set<uint64_t>> & uniqueKmersPerSet,
		UniqueKmerSetHelper::ProcessReadForExtractingPars extractingPars,
		const SimpleKmerHash & hasher,
		const njhseq::StrandBiasKmerCorrectionPars & sbPars,
		uint32_t maxReadsPerSet,
		std::ostream & logOut) {
	auto ret = std::make_shared<std::unordered_map<std::string, std::shared_ptr<njhseq::StrandBiasKmerCorrector>>>();
	if (!sbPars.on_) { return ret; }
	extractingPars.strandBiasCorrectors = nullptr; //never correct while collecting
	std::unordered_map<std::string, std::vector<njhseq::StrandBiasKmerCorrector::BufferedRead>> perSet;
	std::unordered_set<std::string> skip { "undetermined", "multihit" };
	for (const auto & n : extractingPars.compPars.excludeSetNames) { skip.emplace(n); }

	//fragment: shared id for the two mates of a pair (duplicate detection in frequency mode)
	auto keep = [&perSet, &skip, &maxReadsPerSet](const std::string & setName, const seqInfo & seq, bool reverseStrand,
			uint64_t fragment) {
		if (njh::in(setName, skip)) { return; }
		auto & v = perSet[setName];
		if (v.size() >= maxReadsPerSet) { return; }
		njhseq::StrandBiasKmerCorrector::BufferedRead rd;
		rd.seq_ = seq.seq_;
		rd.reverseStrand_ = reverseStrand;
		rd.fragment_ = fragment;
		rd.quals_.reserve(seq.qual_.size());
		for (const auto q : seq.qual_) { rd.quals_.emplace_back(static_cast<uint8_t>(q)); }
		v.emplace_back(std::move(rd));
	};
	if (!pairedOpts.firstName_.empty()) {
		SeqInput reader(pairedOpts);
		reader.openIn();
		PairedRead pseq;
		uint64_t fragment = 0;
		while (reader.readNextRead(pseq)) {
			auto first = UniqueKmerSetHelper::compareReadToSets(pseq.seqBase_, uniqueKmersPerSet, extractingPars.compPars, hasher);
			auto second = UniqueKmerSetHelper::compareReadToSets(pseq.mateSeqBase_, uniqueKmersPerSet, extractingPars.compPars, hasher);
			//compareReadToSets has already flipped each mate into its winning set's orientation
			keep(first.winnerSet, pseq.seqBase_, first.winnerRevComp, fragment);
			keep(second.winnerSet, pseq.mateSeqBase_, !second.winnerRevComp, fragment);
			++fragment;
		}
	}
	if (!singlesOpts.firstName_.empty()) {
		SeqInput reader(singlesOpts);
		reader.openIn();
		seqInfo seq;
		while (reader.readNextRead(seq)) {
			auto res = UniqueKmerSetHelper::compareReadToSets(seq, uniqueKmersPerSet, extractingPars.compPars, hasher);
			const bool secondMate = njh::endsWith(seq.name_, "_secondMate");
			keep(res.winnerSet, seq, secondMate ? !res.winnerRevComp : res.winnerRevComp, std::numeric_limits<uint64_t>::max());
		}
	}
	for (auto & setReads : perSet) {
		if (setReads.second.size() < 50) { continue; }
		auto corrector = std::make_shared<njhseq::StrandBiasKmerCorrector>(sbPars);
		uint32_t vetoed = 0;
		for (uint32_t pass = 0; pass < sbPars.maxPasses_; ++pass) {
			njhseq::StrandBiasKmerCorrector::KmerTable table;
			corrector->tallyReads(table, setReads.second);
			auto candidates = corrector->findCandidates(table);
			if (candidates.empty()) { break; }
			auto quals = corrector->collectOffsetQuality(setReads.second, candidates);
			uint32_t rescued = 0;
			auto found = corrector->buildPass(candidates, table, quals, rescued);
			vetoed += corrector->vetoSupportedAlleles(found, candidates, table, quals, setReads.second);
			if (found.corrections_.empty()) { break; }
			corrector->addQualityRescued(rescued);
			corrector->addPass(std::move(found));
			for (auto & rd : setReads.second) { corrector->applyToSeq(rd.seq_); }
		}
		if (corrector->hasCorrections()) {
			uint32_t nk = 0;
			for (const auto & pass : corrector->passes_) { nk += pass.corrections_.size(); }
			logOut << "strandBiasKmerCorrection\t" << setReads.first << "\t" << setReads.second.size()
					<< "\t" << corrector->numPasses() << "\t" << nk << "\t" << corrector->numQualityRescued() << "\t" << vetoed << std::endl;
			(*ret)[setReads.first] = corrector;
		}
		setReads.second.clear();
		setReads.second.shrink_to_fit();
	}
	return ret;
}

}  // namespace

int kmerSetExpRunner::extractByCountingUniqKmersFromSets(const njh::progutils::CmdArgs & inputCommands){

	uint32_t numThreads = 1;
	bfs::path countTable;
	bfs::path excludesCountTable;
	bfs::path nonUniqueKmerTable;
	uint32_t maxIterations = std::numeric_limits<uint32_t>::max();

	bool keepTemporaryFiles = false;
	bool doNotDoFinalExtract = false;

	uint32_t iterationsToFilterDissimilar = std::numeric_limits<uint32_t>::max();
	UniqueKmerSetHelper::CompareReadToSetPars dissimilarFilterPars;
	dissimilarFilterPars.hardCountOff = 0;
	dissimilarFilterPars.fracCutOff = 0;

	UniqueKmerSetHelper::ProcessReadForExtractingPars extractingPars;
	extractingPars.compPars.hardCountOff = 20;
	extractingPars.compPars.fracCutOff = 0.12;
	extractingPars.compPars.finalHardCountOff = 20;
	extractingPars.compPars.finalFracCutOff = extractingPars.compPars.fracCutOff;
	// extractingPars.compPars.finalFracCutOff = 0.30;
	extractingPars.compPars.kmerLengthForEntropyCalc_ = 2;
	extractingPars.compPars.entropyFilter_ = 1.20;

	extractingPars.compPars.initialExcludeHardCountOff = 60;
	extractingPars.compPars.initialExcludeFracCutOff = 0.25;

	uint32_t iterationToLowerExcludeCutOff = 6;

	bool finalExtractionPairsSeparate = false;



	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();

	//read in options
	bfs::path bamFnp;
	bfs::path bamExtractBedFnp;
	double bamExtractBedPercInRegion = 0.50;
	setUp.setOption(bamExtractBedPercInRegion, "--bamExtractBedPercInRegion", "bam Extract Bed Perc In Region", false);
	bool bamFnpSet = setUp.setOption(bamFnp, "--bamFnp", "bam Fnp", false);
	setUp.setOption(bamExtractBedFnp, "--bamExtractBedFnp", "bam Extract Bed Fnp regions", false);

	setUp.pars_.ioOptions_.revComplMate_ = true;
	bool pairedEndSet = setUp.processReadInNames(njhseq::seqSetUp::pairedReadInFormatsAvailable_, false);
	auto singlesOption = SeqIOOptions();
	bool singlesEndSet = setUp.processJustReadInNames(singlesOption, njhseq::seqSetUp::singleInFormatsAvailable_, !pairedEndSet && !bamFnpSet);
	singlesOption.lowerCaseBases_ = setUp.pars_.ioOptions_.lowerCaseBases_;
	singlesOption.removeGaps_ = setUp.pars_.ioOptions_.removeGaps_;
	singlesOption.includeWhiteSpaceInName_ = setUp.pars_.ioOptions_.includeWhiteSpaceInName_;

	if((singlesEndSet || pairedEndSet) && bamFnpSet){
		setUp.failed_ = true;
		setUp.addWarning(njh::pasteAsStr("can't set both bam file and fastq inputs"));
	}

	setUp.setOption(countTable, "--kmerTable", "countTable, 1)set,2)kmer", true);
	{
		auto possibleNonUniqueTable = njh::files::prependFileBasename(countTable, "nonUniqueKmers_");
		if(exists(possibleNonUniqueTable)){
			nonUniqueKmerTable = possibleNonUniqueTable;
		}
	}
	setUp.setOption(nonUniqueKmerTable, "--nonUniqueKmerTable", "non-unique Kmer Table, 1)set,2)kmer");
	setUp.setOption(extractingPars.compPars.sampleName, "--sampleName", "Name to add to output file", true);
	bool doNotIncludeRevComp = false;
	setUp.setOption(doNotIncludeRevComp, "--doNotIncludeRevComp", "do not include Rev Comp of the input seqs");
	extractingPars.compPars.includeRevComp = !doNotIncludeRevComp;
	//setUp.setOption(extractingPars.compPars.includeRevComp, "--includeRevComp", "include Rev Comp of the input seqs");

//	bool pairsTogether = false;
//	setUp.setOption(pairsTogether, "--pairsTogether", "extract the pairs together");
//	extractingPars.compPars.pairsSeparate = !pairsTogether;
	setUp.setOption(extractingPars.compPars.pairsSeparate, "--pairsSeparate", "extract the pairs separately");

	bool finalExtractionPairsTogether = false;
	setUp.setOption(finalExtractionPairsTogether, "--finalExtractionPairsTogether", "final extraction the pairs together");
	//strand-bias k-mer correction, applied to the reads written by the final extraction
	njhseq::StrandBiasKmerCorrectionPars sbPars;
	sbPars.kLen_ = 25;              //independent of the set table's k, which is for classification
	sbPars.requireQualityDrop_ = true;  //a multigene family needs quality to corroborate every correction
	bool noStrandBiasKmerCorrection = false;
	setUp.setOption(noStrandBiasKmerCorrection, "--noStrandBiasKmerCorrection",
			"Turn off correcting strand-restricted artifact k-mers in the reads written by the final extraction");
	sbPars.on_ = !noStrandBiasKmerCorrection;
	setUp.setOption(sbPars.kLen_, "--strandBiasKmerCorrectionKLen",
			"K-mer length for the strand-bias correction (max 32). Longer values make a window paralog specific and relieve the shared-paralogy veto, but lower the depth per k-mer; 25 suits most sets, raise it for deep, highly paralogous sets");
	setUp.setOption(sbPars.qualDropMin_, "--strandBiasKmerCorrectionQualityDropMin",
			"Required Phred drop of the reference base at the candidate's position, other strand minus artifact strand; in frequency mode, the partner's base minus the candidate's");
	setUp.setOption(sbPars.minPartnerMissingStrandDepth_, "--strandBiasKmerCorrectionMinPartnerDepth",
			"Balanced partner k-mer must have at least this many reads on the strand the candidate is missing from");
	setUp.setOption(sbPars.fisherPhred_, "--strandBiasKmerCorrectionFisherPhred",
			"Two-sided Fisher exact (Phred scaled) of partner vs candidate strand counts must exceed this");
	setUp.setOption(sbPars.maxPasses_, "--strandBiasKmerCorrectionMaxPasses", "Maximum tally-and-correct passes");
	//alternative criterion for data where the strand axis is unusable (stranded RNA-seq): correct a
	//rare k-mer to a dominant one-off neighbour on a depth ratio alone. Not selective -- it also
	//removes real alleles below roughly 1/(multiplier+1) of their locus -- so it is strictly opt-in.
	setUp.setOption(sbPars.maxMinorStrandFrac_, "--strandBiasKmerCorrectionMaxMinorStrandFrac",
			"Scale the minor-strand veto with each k-mer's depth rather than holding it absolute; 0 keeps the absolute threshold");
	bool noStrandBiasKmerDoubletCorrection = !sbPars.doublet_;
	setUp.setOption(noStrandBiasKmerDoubletCorrection, "--noStrandBiasKmerDoubletCorrection",
			"Turn off correction of adjacent two-base artifacts. The single-base criterion cannot reach these, since a k-mer covering the pair differs at two positions and so has no one-off partner. A doublet is held to stricter evidence than a single: no quality relaxation, and a quality collapse required at both positions");
	sbPars.doublet_ = !noStrandBiasKmerDoubletCorrection;
	setUp.setOption(sbPars.doubletRequireSameSubstitution_, "--strandBiasKmerDoubletRequireSameSubstitution",
			"Correct an adjacent doublet only when both bases change the same way; off by default, since requiring it blocks a real artifact class");
	bool noStrandBiasKmerResultValidation = !sbPars.requireCorrectedKmerObserved_;
	setUp.setOption(noStrandBiasKmerResultValidation, "--noStrandBiasKmerResultValidation",
			"Stop checking that a correction's result is sequence the reads actually carry; on by default, it keeps a correction from importing sequence from another copy of a repeat");
	sbPars.requireCorrectedKmerObserved_ = !noStrandBiasKmerResultValidation;
	setUp.setOption(sbPars.validationFlank_, "--strandBiasKmerValidationFlank",
			"Bases of flanking context on each side of a corrected k-mer that the result must also match in the reads; 0 checks the k-mer alone");
	bool kmerFrequencyCorrection = false;
	setUp.setOption(kmerFrequencyCorrection, "--kmerFrequencyCorrection",
			"Use a k-mer frequency criterion instead of the strand/quality one: correct a rare k-mer to a dominant one-off neighbour. For stranded RNA-seq, where strand is confounded with mate and template errors carry normal quality. WARNING: also erases real alleles below about 1/(multiplier+1) of their locus, so minor strains and rare paralogs below that fraction are lost");
	setUp.setOption(sbPars.minDepthRatio_, "--kmerFrequencyCorrectionMultiplier",
			"In frequency mode, the dominant neighbour must be at least this many fold the k-mer being corrected; 10 corrects below about 9% of a locus, 5 below about 17%");
	setUp.setOption(sbPars.freqMinPartnerDepth_, "--kmerFrequencyCorrectionMinPartnerDepth",
			"In frequency mode, the dominant neighbour needs at least this much support before the ratio is trusted");
	setUp.setOption(sbPars.freqMinCandidateTotal_, "--kmerFrequencyCorrectionMinCandidateDepth",
			"In frequency mode, the smallest k-mer depth eligible for correction (1 corrects single-read errors)");
	setUp.setOption(sbPars.freqMaxCandidateTotal_, "--kmerFrequencyCorrectionMaxCandidateDepth",
			"In frequency mode, the largest k-mer depth eligible for correction, regardless of the ratio; 0 removes the cap. Keeps the correction on the singletons and doubletons that carry the error mass, and so leaves alone tandem-repeat genes, where a genuine repeat-copy variant sits at a fraction of the gene's coverage and satisfies the ratio just as an error does");
	bool kmerFrequencyCorrectionKeepDuplicates = !sbPars.freqCollapseDuplicates_;
	setUp.setOption(kmerFrequencyCorrectionKeepDuplicates, "--kmerFrequencyCorrectionKeepDuplicates",
			"In frequency mode, count duplicate fragments in k-mer depth; by default a fragment whose mates are identical to an earlier fragment's counts once");
	sbPars.freqCollapseDuplicates_ = !kmerFrequencyCorrectionKeepDuplicates;
	bool kmerFrequencyCorrectionNoQualityOverride = !sbPars.freqQualityOverride_;
	setUp.setOption(kmerFrequencyCorrectionNoQualityOverride, "--kmerFrequencyCorrectionNoQualityOverride",
			"In frequency mode, turn off the quality override, which corrects candidates over the depth cap or short of the multiplier when their differing base is low quality (mean --strandBiasKmerCorrectionQualityDropMin below the partner's, with few Q30+ calls), including adjacent two-base candidates");
	sbPars.freqQualityOverride_ = !kmerFrequencyCorrectionNoQualityOverride;
	setUp.setOption(sbPars.freqQualityMinDepthRatio_, "--kmerFrequencyCorrectionQualityMultiplier",
			"In frequency mode, the dominant neighbour fold required instead of --kmerFrequencyCorrectionMultiplier when the quality drop is seen");
	setUp.setOption(sbPars.freqQualityMaxHighQualShare_, "--kmerFrequencyCorrectionQualityMaxHighQualShare",
			"In frequency mode, the quality override requires the candidate's share of Q30+ calls at its differing base to be at most this fraction of the partner's share");
	setUp.setOption(sbPars.freqQualityMaxConfidentCalls_, "--kmerFrequencyCorrectionQualityMaxConfidentCalls",
			"In frequency mode, the quality override requires the candidate to have at most this many Q30+ calls at its differing base");
	if (kmerFrequencyCorrection) {
		sbPars.frequencyOnly_ = true;
		sbPars.requireQualityDrop_ = false;
		sbPars.on_ = !noStrandBiasKmerCorrection;
	}
	uint32_t strandBiasMaxReadsPerSet = 1000000;
	setUp.setOption(strandBiasMaxReadsPerSet, "--strandBiasKmerCorrectionMaxReadsPerSet",
			"Cap on reads buffered per set while building the corrections, which bounds memory");
	finalExtractionPairsSeparate = !finalExtractionPairsTogether;
	//setUp.setOption(finalExtractionPairsSeparate, "--finalExtractionPairsSeparate", "extract the pairs separately for final extraction step");


	setUp.setOption(extractingPars.compPars.allowableCharacters_, "--allowableCharacters", "allowable characters");
	bool noInvertedKmerIndex = false;
	setUp.setOption(noInvertedKmerIndex, "--noInvertedKmerIndex",
			"Probe every k-mer set for every read k-mer instead of using one k-mer -> set(s) index. The index gives identical counts and handles k-mers shared between sets; this exists to fall back from if that is ever in doubt");
	extractingPars.compPars.useInvertedKmerIndex_ = !noInvertedKmerIndex;

	setUp.setOption(extractingPars.compPars.entropyFilter_, "--entropyFilter", "entropy Filter cut off");
	setUp.setOption(extractingPars.compPars.kmerLengthForEntropyCalc_, "--kmerLengthForEntropyCalc", "kmerLength For Entropy Calculation");


	setUp.setOption(dissimilarFilterPars.hardCountOff, "--hardCountOffForDissimilarFiltering", "hard Count Off For Dissimilar Filtering");
	setUp.setOption(dissimilarFilterPars.fracCutOff, "--fracCutOffForDissimilarFiltering", "frac Cut Off For Dissimilar Filtering");


	setUp.setOption(iterationsToFilterDissimilar, "--iterationsToFilterDissimilar", "iterations To Filter Dissimilar");

	setUp.setOption(extractingPars.compPars.initialExcludeFracCutOff, "--initialExcludeFracCutOff", "initial for exclusion sets per read Kmer Frac Cut Off");
	setUp.setOption(extractingPars.compPars.initialExcludeHardCountOff, "--initialExcludeHardCountOff", "initial for exclusion sets hard Count Off, do not count sets unless greater thant his number");
	setUp.setOption(iterationToLowerExcludeCutOff, "--iterationToLowerExcludeCutOff", "iterations To Lower Exclusion set Cut Off to the nomral cut offs");

	if(bfs::exists(countTable)){
		auto row = tokenizeString(njh::files::getFirstLine(countTable), "\t");
		extractingPars.compPars.klen = row[1].size();
	}
	setUp.setOption(extractingPars.compPars.fracCutOff, "--readKmerFracCutOff", "per read Kmer Frac Cut Off");
	setUp.setOption(extractingPars.compPars.hardCountOff, "--hardCountOff", "hard Count Off, do not count sets unless greater thant his number");
	if (extractingPars.compPars.hardCountOff > extractingPars.compPars.finalHardCountOff) {
		extractingPars.compPars.finalHardCountOff = extractingPars.compPars.hardCountOff;
	}
	if (extractingPars.compPars.fracCutOff > extractingPars.compPars.finalFracCutOff) {
		extractingPars.compPars.finalFracCutOff = extractingPars.compPars.fracCutOff;
	}
	setUp.setOption(extractingPars.compPars.finalFracCutOff, "--finalReadKmerFracCutOff", "final round extraction per read Kmer Frac Cut Off");
	setUp.setOption(extractingPars.compPars.finalHardCountOff, "--finalHardCountOff", "final round extraction hard Count Off, do not count sets unless greater thant his number");



	extractingPars.smallLenCutOff = extractingPars.compPars.hardCountOff + extractingPars.compPars.klen + 11;
	setUp.setOption(extractingPars.smallLenCutOff, "--smallLenCutOff", "small Len Cut Off of input sequences, if less than this size will skip over");
	setUp.setOption(extractingPars.qPars.qualCheck_, "--qualCheck", "quality score for filtering");
	extractingPars.qPars.qualCheckCutOff_ = 0.50;
	setUp.setOption(extractingPars.qPars.qualCheckCutOff_, "--qualCheckCutOff", "the fraction of bases that have to be above the qualCheck check");

	// setUp.setOption(extractingPars.qPars.checkingQFrac_, "--checkingQFrac", "filtering on quality");
	bool notCheckingQFrac = false;
	setUp.setOption(notCheckingQFrac, "--noCheckingQFrac", "do not filtering on quality");
	extractingPars.qPars.checkingQFrac_ = ! notCheckingQFrac;
	extractingPars.qual_checker = std::make_shared<ReadCheckerQualCheck>(extractingPars.qPars.qualCheck_, extractingPars.qPars.qualCheckCutOff_, false);
	bool doNotFilterOnNs = false;
	setUp.setOption(doNotFilterOnNs, "--noFilterOnNs", "do not filter reads containing Ns");
	extractingPars.filterOnNs = ! doNotFilterOnNs;
	//setUp.setOption(extractingPars.filterOnNs, "--filterOnNs", "filter reads containing Ns");

	setUp.setOption(extractingPars.filterMultiHitReads, "--filterMultiHitReads", "filter Multi Hit Reads");


	setUp.setOption(extractingPars.markReadsPerIteration, "--markReadsPerIteration", "mark Reads Per Iteration");
	setUp.setOption(extractingPars.writeOutFinalKmerSets, "--writeOutFinalKmerSets", "write Out Final Kmer Sets");

	setUp.setOption(excludesCountTable, "--excludesCountTable", "excludes Count Table, 1)set,2)kmer");
	setUp.setOption(extractingPars.compPars.excludeSetNames, "--excludeSetNames", "names of sets of unique kmers to ignore from the --kmerTable");
	//
	setUp.setOption(keepTemporaryFiles, "--keepTemporaryFiles", "Keep Temporary Files");
	setUp.setOption(extractingPars.addingInKmersCountCutOff, "--addingInKmersCountCutOff", "adding In Kmers Count Cut Off");


	setUp.setOption(extractingPars.doReCheckExcludeSets, "--doReCheckExcludeSets", "do Re Check Exclude Sets on iterations, could lead to some reads being recruit away but will drastically speed up run time");

	// bool doNotFinalRecheckOfExcludeSeqs = false;
	// setUp.setOption(doNotFinalRecheckOfExcludeSeqs, "--doNotFinalRecheckOfExcludeSeqs", "redo a unique kmer filter on the exclusion sets during the final re-extraction process");
	//
	// bool doFinalRecheckOfExcludeSeqs = !doNotFinalRecheckOfExcludeSeqs;
	bool doFinalRecheckOfExcludeSeqs = false;
	setUp.setOption(doFinalRecheckOfExcludeSeqs, "--doFinalRecheckOfExcludeSeqs", "Normal the exclusion kmer set is left alone, if this flag is set, it will be filter during the final re-extraction process");

	setUp.setOption(doNotDoFinalExtract, "--doNotDoFinalExtract", "do Not Do Final Extract");
	//setUp.setOption(extractAfterIterating, "--extractAfterIterating", "Extract After Iterating");
	setUp.setOption(maxIterations, "--maxIterations", "max Iterations to perform");

	bool finalWriteOutMultiHitSeparate = false;
	setUp.setOption(finalWriteOutMultiHitSeparate, "--finalWriteOutMultiHitSeparate", "final Write Out Multi Hit Separate");
	setUp.setOption(extractingPars.writeMultiHitSeparate, "--writeMultiHitSeparate", "Write Out Multi Hit Separate");



	setUp.setOption(extractingPars.doNotWriteUndetermined, "--doNotWriteUndetermined", "do Not Write Undetermined");
	setUp.setOption(extractingPars.writeOutExclude, "--writeOutExclude", "write Out Excluded reads that match the excluded kmer sets");

	setUp.setOption(numThreads, "--numThreads", "numThreads");
	setUp.processDirectoryOutputName(true);

	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	njh::stopWatch watch;
	watch.setLapName("initial");


	UniqueKmerSetHelper::ProcessReadForExtractingCounts initialCounts;

	std::mutex mut;
	OutOptions bamExtractOut(njh::files::make_path(setUp.pars_.directoryName_, "bamExtract"));
	bamExtractOut.append_ = true;
	if (bamFnpSet) {
		watch.startNewLap(njh::pasteAsStr("bam extract files"));
		njh::files::checkExistenceThrow(bamFnp, __PRETTY_FUNCTION__);
		BamExtractor bExtractor(setUp.pars_.verbose_);
		if (bamExtractBedFnp.empty()) {
			bExtractor.writeExtractReadsFromBam(bamFnp, bamExtractOut);
		} else {
			auto regions = gatherRegions(bamExtractBedFnp.string(), "", setUp.pars_.verbose_);
			for (const auto &region: regions) {
				bExtractor.writeExtractReadsFromBamRegion(bamFnp,
																									region, bamExtractBedPercInRegion, bamExtractOut);
			}
		}

		bfs::path bamExtractR1Fnp = bamExtractOut.outFilename_.string() + "_R1.fastq.gz";
		bfs::path bamExtractR2Fnp = bamExtractOut.outFilename_.string() + "_R2.fastq.gz";
		bfs::path bamExtractFnp = bamExtractOut.outFilename_.string() + ".fastq.gz";
		if (exists(bamExtractR1Fnp)) {
			setUp.pars_.ioOptions_ = SeqIOOptions::genPairedInGz(bamExtractR1Fnp, bamExtractR2Fnp);
//			setUp.pars_.ioOptions_.inFormat_ = SeqIOOptions::inFormats::FASTQPAIREDGZ;
//			setUp.pars_.ioOptions_.outFormat_ = SeqIOOptions::outFormats::FASTQPAIREDGZ;
//			setUp.pars_.ioOptions_.firstName_ = bamExtractR1Fnp;
//			setUp.pars_.ioOptions_.firstName_ = bamExtractR2Fnp;
		}
		if (exists(bamExtractFnp)) {
			singlesOption = SeqIOOptions::genFastqInGz(bamExtractFnp);
//			singlesOption.inFormat_ = SeqIOOptions::inFormats::FASTQGZ;
//			singlesOption.outFormat_ = SeqIOOptions::outFormats::FASTQGZ;
//			singlesOption.firstName_ = bamExtractFnp;
		}
		if(setUp.pars_.verbose_){
			std::cout << watch.getLapName() << "\t" << watch.timeLapFormatted() <<std::endl;
		}
	}

	watch.startNewLap("reading in unique kmer table");


	std::unordered_map<std::string, std::unordered_set<uint64_t>> uniqueKmersPerSet = UniqueKmerSetHelper::readInUniqueKmerTablePerSet(
					countTable);
	{
		//the index handles shared k-mers, but whether the sets are actually exclusive changes how the
		//winner-takes-most classification behaves, so report it rather than assume it
		auto exclusivityIndex = UniqueKmerSetHelper::getInvertedIndex(uniqueKmersPerSet);
		if (nullptr != exclusivityIndex) {
			std::cout << "kmerSetExclusivity\tsets\t" << exclusivityIndex->nSets_
					<< "\ttotalKmers\t" << exclusivityIndex->totalKmers_
					<< "\tsharedKmers\t" << exclusivityIndex->nSharedKmers_
					<< "\texclusive\t" << njh::boolToStr(0 == exclusivityIndex->nSharedKmers_) << std::endl;
		}
	}
	std::string nonUniqueRegionName = "NON_UNIQUE";
	std::unordered_set<uint64_t> nonUniqueKmersPerSet;
	if (!nonUniqueKmerTable.empty()) {
		nonUniqueKmersPerSet = UniqueKmerSetHelper::readInUniqueKmerTableSetsCollapsed(nonUniqueKmerTable);
	}

	if(setUp.pars_.verbose_){
		std::cout << watch.getLapName() << "\t" << watch.timeLapFormatted() <<std::endl;
		std::cout << "klen: " << extractingPars.compPars.klen << std::endl;
		auto setNames = njh::getSetOfMapKeys(uniqueKmersPerSet);
		for(const auto & setName : setNames){
			std::cout << setName << "\t" << uniqueKmersPerSet[setName].size() << std::endl;
		}
		std::cout << "NON_UNIQUE" << "\t" << nonUniqueKmersPerSet.size() << std::endl;
	}
	MultiSeqIO initialSeqOut;
	watch.startNewLap("initial scan");

	OutOptions outCountsOpts(njh::files::make_path(setUp.pars_.directoryName_, "extractionCounts.tab.txt"));
	OutOptions inversePairOutCountsOpts(njh::files::make_path(setUp.pars_.directoryName_, "inversePairsCounts.tab.txt"));
	OutOptions multiHitOutCountsOpts(njh::files::make_path(setUp.pars_.directoryName_, "multihitCounts.tab.txt"));

	OutputStream outCounts(outCountsOpts);
	OutputStream inversePairOutCountsOut(inversePairOutCountsOpts);
	OutputStream multiHitOutCountsOut(multiHitOutCountsOpts);

	UniqueKmerSetHelper::ProcessReadForExtractingCounts::writeOutCountsHeader(outCounts, extractingPars);
	inversePairOutCountsOut << njh::conToStr(toVecStr("iteration", "sample", "totalReads", "target", "count", "frac"), "\t") << std::endl;
	multiHitOutCountsOut << njh::conToStr(toVecStr("iteration", "sample", "totalReads", "hits", "count", "frac"), "\t") << std::endl;

	std::unordered_map<std::string, UniqueKmerSetHelper::FilePositons> positionsAfterLastIteration;
	{//initial

		//add outputs
		VecStr names = getVectorOfMapKeys(uniqueKmersPerSet);
		for(const auto & name : names){
			auto seqOutOpts = SeqIOOptions::genPairedOut(njh::files::make_path(setUp.pars_.directoryName_, name));
			initialSeqOut.addReader(njh::pasteAsStr(name, "-paired"), seqOutOpts);
		}
		initialSeqOut.addReader("undetermined-paired", SeqIOOptions::genPairedOut(njh::files::make_path(setUp.pars_.directoryName_, "undetermined")));
		initialSeqOut.addReader("multihit-paired", SeqIOOptions::genPairedOut(njh::files::make_path(setUp.pars_.directoryName_, "multihit")));

		for(const auto & name : names){
			auto seqOutOpts = SeqIOOptions::genFastqOut(njh::files::make_path(setUp.pars_.directoryName_, name));
			initialSeqOut.addReader(njh::pasteAsStr(name, "-single"), seqOutOpts);
		}
		initialSeqOut.addReader("undetermined-single", SeqIOOptions::genFastqOut(njh::files::make_path(setUp.pars_.directoryName_, "undetermined")));
		initialSeqOut.addReader("multihit-single", SeqIOOptions::genFastqOut(njh::files::make_path(setUp.pars_.directoryName_, "multihit")));

		//to make below count correctly, make sure mate is automatically reversed complemented
		setUp.pars_.ioOptions_.revComplMate_ = true;

		if (!setUp.pars_.ioOptions_.firstName_.empty()) {
			SeqInput reader(setUp.pars_.ioOptions_);
			reader.openIn();
			std::function<void()> readInComp;

			if(extractingPars.compPars.pairsSeparate){
				readInComp = [&reader, &uniqueKmersPerSet,
								&mut, &extractingPars, &initialSeqOut, &initialCounts

				]() {
					SimpleKmerHash hasher;
					PairedRead pseq;
					UniqueKmerSetHelper::ProcessReadForExtractingCounts currentCounts;
					while (reader.readNextReadLock(pseq)) {
						bool firstMatePassQcCheck = UniqueKmerSetHelper::readPassesQCChecks(pseq.seqBase_, extractingPars,currentCounts);
						bool secondMatePassQcCheck = UniqueKmerSetHelper::readPassesQCChecks(pseq.mateSeqBase_, extractingPars,currentCounts);
						if(firstMatePassQcCheck && secondMatePassQcCheck) {
							UniqueKmerSetHelper::processReadForExtractingPairsSeparate(pseq, uniqueKmersPerSet,
																							extractingPars, hasher, initialSeqOut,
																							currentCounts
																							);
						} else if(firstMatePassQcCheck) {
							pseq.seqBase_.name_.append("_firstMate");
							UniqueKmerSetHelper::processReadForExtracting(pseq.seqBase_, uniqueKmersPerSet,
																		extractingPars, hasher, initialSeqOut,
																		currentCounts);
						} else if(secondMatePassQcCheck) {
							pseq.mateSeqBase_.name_.append("_secondMate");
							UniqueKmerSetHelper::processReadForExtracting(pseq.mateSeqBase_, uniqueKmersPerSet,
																		extractingPars, hasher, initialSeqOut,
																		currentCounts);
						}
					}
					{
						std::lock_guard<std::mutex> lockGuard(mut);
						initialCounts.addOtherCounts(currentCounts);
					}
				};
			}else{
				readInComp = [&reader, &uniqueKmersPerSet,
								&mut, &extractingPars, &initialSeqOut, &initialCounts

				]() {
					SimpleKmerHash hasher;
					PairedRead pseq;
					UniqueKmerSetHelper::ProcessReadForExtractingCounts currentCounts;
					while (reader.readNextReadLock(pseq)) {
						bool firstMatePassQcCheck = UniqueKmerSetHelper::readPassesQCChecks(pseq.seqBase_, extractingPars,currentCounts);
						bool secondMatePassQcCheck = UniqueKmerSetHelper::readPassesQCChecks(pseq.mateSeqBase_, extractingPars,currentCounts);
			      if(firstMatePassQcCheck && secondMatePassQcCheck) {
				      UniqueKmerSetHelper::processReadForExtractingPairsTogether(pseq, uniqueKmersPerSet,
																							extractingPars, hasher, initialSeqOut,
																							currentCounts
																							);
			      } else if(firstMatePassQcCheck) {
			      	pseq.seqBase_.name_.append("_firstMate");
			      	UniqueKmerSetHelper::processReadForExtracting(pseq.seqBase_, uniqueKmersPerSet,
																		extractingPars, hasher, initialSeqOut,
																		currentCounts);
			      } else if(secondMatePassQcCheck) {
			      	pseq.mateSeqBase_.name_.append("_secondMate");
			      	UniqueKmerSetHelper::processReadForExtracting(pseq.mateSeqBase_, uniqueKmersPerSet,
																		extractingPars, hasher, initialSeqOut,
																		currentCounts);
			      }
					}
					{
						std::lock_guard<std::mutex> lockGuard(mut);
						initialCounts.addOtherCounts(currentCounts);
					}
				};
			}
			njh::concurrent::runVoidFunctionThreaded(readInComp, numThreads);
		}
		if (!singlesOption.firstName_.empty()) {
			SeqInput reader(singlesOption);
			reader.openIn();
			std::function<void()> readInComp = [&reader, &uniqueKmersPerSet,
							&mut, &extractingPars, &initialSeqOut, &initialCounts

//							, &testReads
			]() {
				SimpleKmerHash hasher;
				seqInfo seq;
				UniqueKmerSetHelper::ProcessReadForExtractingCounts currentCounts;
				while (reader.readNextReadLock(seq)) {
					if(UniqueKmerSetHelper::readPassesQCChecks(seq,extractingPars,currentCounts)) {
						UniqueKmerSetHelper::processReadForExtracting(seq, uniqueKmersPerSet,
																								extractingPars, hasher, initialSeqOut,
																								currentCounts);
					}
				}
				{
					std::lock_guard<std::mutex> lockGuard(mut);
					initialCounts.addOtherCounts(currentCounts);
				}
			};
			njh::concurrent::runVoidFunctionThreaded(readInComp, numThreads);
		}
		initialCounts.writeOutCounts(outCounts,extractingPars,uniqueKmersPerSet,"initial");
		{
			auto totalReadCnt = initialCounts.getTotalCounts();
			for(const auto & set : initialCounts.inversePairReadCountsPerSet) {
				inversePairOutCountsOut << "initial"
				<< "\t" << extractingPars.compPars.sampleName
				<< "\t" << totalReadCnt
				<< "\t" << set.first
				<< "\t" << set.second
				<< "\t" << static_cast<double>(set.second)/static_cast<double>(totalReadCnt) << std::endl;
			}
		}
		{
			auto totalReadCnt = initialCounts.getTotalCounts();
			for(const auto & set : initialCounts.multiHitReadCounts) {
				multiHitOutCountsOut << "initial"
				<< "\t" << extractingPars.compPars.sampleName
				<< "\t" << totalReadCnt
				<< "\t" << set.first
				<< "\t" << set.second
				<< "\t" << static_cast<double>(set.second)/static_cast<double>(totalReadCnt) << std::endl;
			}
		}
	}
	initialSeqOut.closeOutForReopeningAll();

	{
		auto kmerSetNames = njh::getSetOfMapKeys(uniqueKmersPerSet);
		for(const auto & kmerSetName : kmerSetNames){
			if(njh::in(kmerSetName, extractingPars.compPars.excludeSetNames)){
				continue;
			}
			UniqueKmerSetHelper::FilePositons positions;
			auto r1Fnp = njh::files::make_path(setUp.pars_.directoryName_, kmerSetName + "_R1.fastq");
			auto r2Fnp =  njh::files::make_path(setUp.pars_.directoryName_, kmerSetName + "_R2.fastq");
			auto singleFnp = njh::files::make_path(setUp.pars_.directoryName_, kmerSetName + ".fastq");
			positions.r1FnpEnd = njh::files::getFilePositionsEndOfFile(r1Fnp);
			positions.r2FnpEnd = njh::files::getFilePositionsEndOfFile(r2Fnp);
			positions.singleFnpEnd = njh::files::getFilePositionsEndOfFile(singleFnp);
			positionsAfterLastIteration[kmerSetName] = positions;
		}
	}
	//
	if(setUp.pars_.verbose_) {
		std::cout << "iterNumber: " << "initial" << std::endl;
		std::cout << "Total Initial Reads: " << initialCounts.getTotalCounts() << std::endl;
		std::cout << "Total Initial Undetermined Reads: " << initialCounts.genTotalUndeterminedCount() << " ("
							<< initialCounts.genTotalUndeterminedCount() * 100 / static_cast<double>(initialCounts.getTotalCounts())
							<< "%)" << std::endl;
		std::cout << "Total Initial Determined Reads: " << initialCounts.genTotalDeterminedCount() << " ("
							<< initialCounts.genTotalDeterminedCount() * 100 / static_cast<double>(initialCounts.getTotalCounts())
							<< "%)" << std::endl;
		std::cout << "Total Initial Less Than Small Len Cut Off Count: " << initialCounts.smallLenCutOffCount << " ("
							<< initialCounts.smallLenCutOffCount * 100 / static_cast<double>(initialCounts.getTotalCounts()) << "%)"
							<< std::endl;
		std::cout << "Total Initial Contains Ns Count: " << initialCounts.containsNs << " ("
							<< initialCounts.containsNs * 100 / static_cast<double>(initialCounts.getTotalCounts()) << "%)"
							<< std::endl;
		std::cout << "Total Initial Poor Quality Count: " << initialCounts.poorQualityCount << " ("
							<< initialCounts.poorQualityCount * 100 / static_cast<double>(initialCounts.getTotalCounts()) << "%)"
							<< std::endl;
		std::cout << "Total Initial Inverse Pair Count: " << initialCounts.genTotalInversePairCount() << " ("
					<< initialCounts.genTotalInversePairCount() * 100 / static_cast<double>(initialCounts.getTotalCounts()) << "%)"
					<< std::endl;
		std::cout << "Total Initial Multihit Count: " << initialCounts.genTotalMultihitCount() << " ("
					<< initialCounts.genTotalMultihitCount() * 100 / static_cast<double>(initialCounts.getTotalCounts()) << "%)"
					<< std::endl;
		std::cout << watch.getLapName() << "\t" << watch.timeLapFormatted() << std::endl;
	}
	//read in the extract reads and add to the unique kmer sets
	bool iterate = true;
	uint64_t totalUndeterminedExtracted = initialCounts.genTotalUndeterminedCount();
	uint32_t iterNumber = 0;
	if(0 == totalUndeterminedExtracted){
		iterate = false;
	}




	auto original_initialCounts = initialCounts;
	initialCounts.readCountsPerSet[false]["undetermined"] = 0;
	initialCounts.readCountsPerSet[true]["undetermined"] = 0;
	uint64_t totalDeterminedReads = initialCounts.genTotalDeterminedCount();
//	std::cout << __FILE__ << " " << __LINE__ << std::endl;
	while(iterate && iterNumber < maxIterations){
//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		watch.startNewLap(njh::pasteAsStr(iterNumber, "- read in new kmers"));
//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		initialSeqOut.closeOutForReopeningAll();
//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		auto rawUniqKmerSets = njh::getVecOfMapKeys(uniqueKmersPerSet);
		removeElement(rawUniqKmerSets, std::string("undetermined"));
		removeElements(rawUniqKmerSets, VecStr(extractingPars.compPars.excludeSetNames.begin(), extractingPars.compPars.excludeSetNames.end()));
		std::unordered_map<std::string, std::unordered_map<uint64_t, uint32_t>> rawKmersPerInput;
		if (iterNumber > 0) {
			rawKmersPerInput = UniqueKmerSetHelper::readInNewKmersFromExtractedReads(
							setUp.pars_.directoryName_, rawUniqKmerSets, extractingPars, positionsAfterLastIteration);
		} else {
//			std::cout << __FILE__ << " " << __LINE__ << std::endl;
			rawKmersPerInput = UniqueKmerSetHelper::readInNewKmersFromExtractedReads(
							setUp.pars_.directoryName_, rawUniqKmerSets, extractingPars);
//			std::cout << __FILE__ << " " << __LINE__ << std::endl;
		}
//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		if (setUp.pars_.verbose_ && setUp.pars_.debug_) {
			std::cout << "new kmer counts for: " << rawKmersPerInput.size() << " sets" << std::endl;
			auto setNames = njh::getSetOfMapKeys(rawKmersPerInput);
			for (const auto &setName: setNames) {
				std::cout << setName << "\t" << rawKmersPerInput[setName].size() << std::endl;
			}
			std::cout << "NON_UNIQUE" << "\t" << nonUniqueKmersPerSet.size() << std::endl << std::endl;
		}

//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		watch.startNewLap(njh::pasteAsStr(iterNumber, "- check new kmers against current set"));
		std::unordered_map<std::string, std::unordered_set<uint64_t>> outputUniqueKmersPerSet = UniqueKmerSetHelper::filterReExtractedKmersForNonUnique(rawKmersPerInput,
																																																																										extractingPars,
																																																																										uniqueKmersPerSet,
																																																																										nonUniqueKmersPerSet);
//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		if(!extractingPars.doReCheckExcludeSets){
			for(const auto & excludeName : extractingPars.compPars.excludeSetNames){
				outputUniqueKmersPerSet[excludeName] = std::move(uniqueKmersPerSet[excludeName]);
			}
		}

//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		watch.startNewLap(njh::pasteAsStr(iterNumber, "- copy new set"));
		uniqueKmersPerSet = std::move(outputUniqueKmersPerSet);
		watch.startNewLap(njh::pasteAsStr(iterNumber, "- scan again"));
//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		if (setUp.pars_.verbose_ && setUp.pars_.debug_) {
			auto setNames = njh::getSetOfMapKeys(uniqueKmersPerSet);
			for (const auto &setName: setNames) {
				std::cout << setName << "\t" << uniqueKmersPerSet[setName].size() << std::endl;
			}
			std::cout << "NON_UNIQUE" << "\t" << nonUniqueKmersPerSet.size() << std::endl;
		}

		if(iterNumber >= iterationToLowerExcludeCutOff){
			extractingPars.compPars.initialExcludeHardCountOff = extractingPars.compPars.hardCountOff;
			extractingPars.compPars.initialExcludeFracCutOff = extractingPars.compPars.fracCutOff;
		}

		bfs::path orig_undeterminedInput = njh::files::make_path(setUp.pars_.directoryName_, "undetermined.fastq");
		bfs::path orig_undeterminedInputR1 = njh::files::make_path(setUp.pars_.directoryName_, "undetermined_R1.fastq");
		bfs::path orig_undeterminedInputR2 = njh::files::make_path(setUp.pars_.directoryName_, "undetermined_R2.fastq");

		bfs::path undeterminedInput = njh::files::make_path(setUp.pars_.directoryName_, njh::pasteAsStr(iterNumber, "_undetermined.fastq"));
		bfs::path undeterminedInputR1 = njh::files::make_path(setUp.pars_.directoryName_, njh::pasteAsStr(iterNumber, "_undetermined_R1.fastq"));
		bfs::path undeterminedInputR2 = njh::files::make_path(setUp.pars_.directoryName_, njh::pasteAsStr(iterNumber, "_undetermined_R2.fastq"));

		SeqIOOptions undeterminedPairedInOpts;
		SeqIOOptions undeterminedSingleInOpts;
		bool pairedInUndeterminedExists = bfs::exists(orig_undeterminedInputR1);
		bool singleInUndeterminedExists = bfs::exists(orig_undeterminedInput);

		if (iterNumber == iterationsToFilterDissimilar) {
			//add outputs
			auto extractingParsForFiltering = extractingPars;
			extractingParsForFiltering.compPars.hardCountOff = dissimilarFilterPars.hardCountOff;
			extractingParsForFiltering.compPars.fracCutOff = dissimilarFilterPars.fracCutOff;
			extractingParsForFiltering.compPars.initialExcludeHardCountOff = dissimilarFilterPars.hardCountOff;
			extractingParsForFiltering.compPars.initialExcludeFracCutOff = dissimilarFilterPars.fracCutOff;

			extractingParsForFiltering.compPars.excludeSetNames = extractingPars.compPars.excludeSetNames;
			MultiSeqIO filteringSeqOut;
			std::string keep = "undetermined";
			std::string filter = "dissimilar";
			filteringSeqOut.addReader("undetermined-paired", SeqIOOptions::genPairedOut(njh::files::make_path(setUp.pars_.directoryName_, njh::pasteAsStr(iterNumber, "_undetermined"))));
			filteringSeqOut.addReader("undetermined-single", SeqIOOptions::genFastqOut(njh::files::make_path(setUp.pars_.directoryName_, njh::pasteAsStr(iterNumber, "_undetermined"))));
			filteringSeqOut.addReader("multihit-paired", SeqIOOptions::genPairedOut(njh::files::make_path(setUp.pars_.directoryName_, njh::pasteAsStr(iterNumber, "_multihit"))));
			filteringSeqOut.addReader("multihit-single", SeqIOOptions::genFastqOut(njh::files::make_path(setUp.pars_.directoryName_, njh::pasteAsStr(iterNumber, "_multihit"))));
			filteringSeqOut.addReader("dissimilar-paired", SeqIOOptions::genPairedOut(njh::files::make_path(setUp.pars_.directoryName_, njh::pasteAsStr(iterNumber, "_dissimilar"))));
			filteringSeqOut.addReader("dissimilar-single", SeqIOOptions::genFastqOut(njh::files::make_path(setUp.pars_.directoryName_, njh::pasteAsStr(iterNumber, "_dissimilar"))));
			if (pairedInUndeterminedExists) {
				SeqIOOptions undeterminedPairedInOptsForFiltering;
				undeterminedPairedInOptsForFiltering.firstName_ = orig_undeterminedInputR1;
				undeterminedPairedInOptsForFiltering.secondName_ = orig_undeterminedInputR2;
				undeterminedPairedInOptsForFiltering.inFormat_ = SeqIOOptions::inFormats::FASTQPAIRED;
				undeterminedPairedInOptsForFiltering.revComplMate_ = true;


				//re-extract paired in
				SeqInput reader(undeterminedPairedInOptsForFiltering);
				reader.openIn();
				std::function<void()> readInComp;

				if(extractingPars.compPars.pairsSeparate){
					readInComp = [&reader, &uniqueKmersPerSet,
									&mut, &extractingParsForFiltering, &filteringSeqOut, &initialCounts,
									&iterNumber,&keep,&filter
					]() {
						SimpleKmerHash hasher;
						PairedRead pseq;
						UniqueKmerSetHelper::ProcessReadForExtractingCounts currentCounts;
						while (reader.readNextReadLock(pseq)) {
							UniqueKmerSetHelper::processReadForFilteringPairsSeparate(pseq, uniqueKmersPerSet,
																																				extractingParsForFiltering, hasher, filteringSeqOut,
																																				 currentCounts,
																																				 njh::pasteAsStr(iterNumber),
																																				 filter,
																																				 keep);
						}
						{
							std::lock_guard<std::mutex> lockGuard(mut);
							initialCounts.addOtherCounts(currentCounts);
						}
					};
				} else {
					readInComp = [&reader, &uniqueKmersPerSet,
									&mut, &extractingParsForFiltering, &filteringSeqOut, &initialCounts,
									&iterNumber,&keep,&filter
					]() {
						SimpleKmerHash hasher;
						PairedRead pseq;
						UniqueKmerSetHelper::ProcessReadForExtractingCounts currentCounts;
						while (reader.readNextReadLock(pseq)) {
							UniqueKmerSetHelper::processReadForFilteringPairsTogether(pseq, uniqueKmersPerSet,
																																				extractingParsForFiltering, hasher, filteringSeqOut,
																																				 currentCounts,
																																				 njh::pasteAsStr(iterNumber),
																																				filter,
																																				keep);
						}
						{
							std::lock_guard<std::mutex> lockGuard(mut);
							initialCounts.addOtherCounts(currentCounts);
						}
					};
				}
				njh::concurrent::runVoidFunctionThreaded(readInComp, numThreads);
			}
			if (singleInUndeterminedExists) {
				//re-extract single in
				SeqIOOptions undeterminedSingleInOptsForFiltering;
				undeterminedSingleInOptsForFiltering.firstName_ = orig_undeterminedInput;
				undeterminedSingleInOptsForFiltering.inFormat_ = SeqIOOptions::inFormats::FASTQ;
				SeqInput reader(undeterminedSingleInOptsForFiltering);
				reader.openIn();
				std::function<void()> readInComp = [&reader, &uniqueKmersPerSet,
								&mut, &extractingParsForFiltering, &filteringSeqOut, &initialCounts
								, &iterNumber,&keep,&filter]() {
					SimpleKmerHash hasher;
					seqInfo seq;
					UniqueKmerSetHelper::ProcessReadForExtractingCounts currentCounts;
					while (reader.readNextReadLock(seq)) {
						UniqueKmerSetHelper::processReadForFiltering(seq, uniqueKmersPerSet,
																												 extractingParsForFiltering, hasher, filteringSeqOut,
																												 currentCounts,
																												 njh::pasteAsStr(iterNumber),
																												 filter,
																												 keep);
					}
					{
						std::lock_guard<std::mutex> lockGuard(mut);
						initialCounts.addOtherCounts(currentCounts);
					}
				};
				njh::concurrent::runVoidFunctionThreaded(readInComp, numThreads);
			}
			if (pairedInUndeterminedExists) {
				bfs::remove(orig_undeterminedInputR1);
				bfs::remove(orig_undeterminedInputR2);
				undeterminedPairedInOpts.firstName_ = undeterminedInputR1;
				undeterminedPairedInOpts.secondName_ = undeterminedInputR2;
				undeterminedPairedInOpts.inFormat_ = SeqIOOptions::inFormats::FASTQPAIRED;
				undeterminedPairedInOpts.revComplMate_ = true;
			}
			if (singleInUndeterminedExists) {
				bfs::remove(orig_undeterminedInput);
				undeterminedSingleInOpts.firstName_ = undeterminedInput;
				undeterminedSingleInOpts.inFormat_ = SeqIOOptions::inFormats::FASTQ;
			}
		} else {
			if (pairedInUndeterminedExists) {
				bfs::rename(orig_undeterminedInputR1, undeterminedInputR1);
				bfs::rename(orig_undeterminedInputR2, undeterminedInputR2);
				undeterminedPairedInOpts.firstName_ = undeterminedInputR1;
				undeterminedPairedInOpts.secondName_ = undeterminedInputR2;
				undeterminedPairedInOpts.inFormat_ = SeqIOOptions::inFormats::FASTQPAIRED;
				undeterminedPairedInOpts.revComplMate_ = true;
				//to make below count correctly, make sure mate is automatically reversed complemented
			}
			//std::cout << __FILE__ << " " << __LINE__ << std::endl;
			if (singleInUndeterminedExists) {
				bfs::rename(orig_undeterminedInput, undeterminedInput);
				undeterminedSingleInOpts.firstName_ = undeterminedInput;
				undeterminedSingleInOpts.inFormat_ = SeqIOOptions::inFormats::FASTQ;
			}
		}

		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		{
			//determine positions before extracting more
			auto kmerSetNames = njh::getSetOfMapKeys(uniqueKmersPerSet);
			for(const auto & kmerSetName : kmerSetNames){
				if(njh::in(kmerSetName, extractingPars.compPars.excludeSetNames)){
					continue;
				}
				UniqueKmerSetHelper::FilePositons positions;
				auto r1Fnp = njh::files::make_path(setUp.pars_.directoryName_, kmerSetName + "_R1.fastq");
				auto r2Fnp =  njh::files::make_path(setUp.pars_.directoryName_, kmerSetName + "_R2.fastq");
				auto singleFnp = njh::files::make_path(setUp.pars_.directoryName_, kmerSetName + ".fastq");
				positions.r1FnpEnd = njh::files::getFilePositionsEndOfFile(r1Fnp);
				positions.r2FnpEnd = njh::files::getFilePositionsEndOfFile(r2Fnp);
				positions.singleFnpEnd = njh::files::getFilePositionsEndOfFile(singleFnp);
				positionsAfterLastIteration[kmerSetName] = positions;
			}
		}
		if (pairedInUndeterminedExists) {
			//std::cout << __FILE__ << " " << __LINE__ << std::endl;
			//re-extract paired in
			SeqInput reader(undeterminedPairedInOpts);
			reader.openIn();
			std::function<void()> readInComp;

			if(extractingPars.compPars.pairsSeparate){
				readInComp = [&reader, &uniqueKmersPerSet,
								&mut, &extractingPars, &initialSeqOut, &initialCounts
				]() {
					SimpleKmerHash hasher;
					PairedRead pseq;
					UniqueKmerSetHelper::ProcessReadForExtractingCounts currentCounts;
					while (reader.readNextReadLock(pseq)) {
						UniqueKmerSetHelper::processReadForExtractingPairsSeparate(pseq, uniqueKmersPerSet,
																																			 extractingPars, hasher, initialSeqOut,
																																			 currentCounts);
					}
					{
						std::lock_guard<std::mutex> lockGuard(mut);
						initialCounts.addOtherCounts(currentCounts);
					}
				};
			}else{
				readInComp = [&reader, &uniqueKmersPerSet,
								&mut, &extractingPars, &initialSeqOut, &initialCounts
				]() {
					SimpleKmerHash hasher;
					PairedRead pseq;
					UniqueKmerSetHelper::ProcessReadForExtractingCounts currentCounts;
					while (reader.readNextReadLock(pseq)) {
						UniqueKmerSetHelper::processReadForExtractingPairsTogether(pseq, uniqueKmersPerSet,
																																			 extractingPars, hasher, initialSeqOut,
																																			 currentCounts);
					}
					{
						std::lock_guard<std::mutex> lockGuard(mut);
						initialCounts.addOtherCounts(currentCounts);
					}
				};
			}
			njh::concurrent::runVoidFunctionThreaded(readInComp, numThreads);
		}
		if(singleInUndeterminedExists) {
			//re-extract single in
			SeqInput reader(undeterminedSingleInOpts);
			reader.openIn();
			std::function<void()> readInComp = [&reader, &uniqueKmersPerSet,
							&mut, &extractingPars, &initialSeqOut, &initialCounts]() {
				SimpleKmerHash hasher;
				seqInfo seq;
				UniqueKmerSetHelper::ProcessReadForExtractingCounts currentCounts;
				while (reader.readNextReadLock(seq)) {
					UniqueKmerSetHelper::processReadForExtracting(seq, uniqueKmersPerSet,
																												extractingPars, hasher, initialSeqOut,
																												currentCounts);
				}

				{
					std::lock_guard<std::mutex> lockGuard(mut);
					initialCounts.addOtherCounts(currentCounts);
				}
			};

			njh::concurrent::runVoidFunctionThreaded(readInComp, numThreads);
		}
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		initialCounts.writeOutCounts(outCounts,extractingPars,uniqueKmersPerSet,njh::pasteAsStr(iterNumber));
		{
			auto totalReadCnt = initialCounts.getTotalCounts();
			for(const auto & set : initialCounts.inversePairReadCountsPerSet) {
				inversePairOutCountsOut << njh::pasteAsStr(iterNumber)
				<< "\t" << extractingPars.compPars.sampleName
				<< "\t" << totalReadCnt
				<< "\t" << set.first
				<< "\t" << set.second
				<< "\t" << static_cast<double>(set.second)/static_cast<double>(totalReadCnt) << std::endl;
			}
		}
		{
			auto totalReadCnt = initialCounts.getTotalCounts();
			for(const auto & set : initialCounts.multiHitReadCounts) {
				multiHitOutCountsOut << njh::pasteAsStr(iterNumber)
				<< "\t" << extractingPars.compPars.sampleName
				<< "\t" << totalReadCnt
				<< "\t" << set.first
				<< "\t" << set.second
				<< "\t" << static_cast<double>(set.second)/static_cast<double>(totalReadCnt) << std::endl;
			}
		}
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		uint64_t currentTotalUndeterminedReads = initialCounts.genTotalUndeterminedCount();
		uint64_t currentTotalReads = initialCounts.getTotalCounts();
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		initialCounts.readCountsPerSet[false]["undetermined"] = 0;
		initialCounts.readCountsPerSet[true]["undetermined"] = 0;
		uint64_t currentTotalDeterminedReads = initialCounts.genTotalDeterminedCount();
		if(setUp.pars_.verbose_){
			std::cout << "iterNumber: " << iterNumber << std::endl;
			std::cout << "Determined: " << currentTotalDeterminedReads - totalDeterminedReads << " more reads" <<std::endl;
			std::cout << "Current counts: " << std::endl;
			std::cout << "\tTotal Reads: " << currentTotalReads << std::endl;
			std::cout << "\tTotal Undetermined Reads: " << currentTotalUndeterminedReads << " (" << currentTotalUndeterminedReads * 100 /static_cast<double>(currentTotalReads)<< "%)" << std::endl;
			std::cout << "\tTotal Determined Reads: " << currentTotalDeterminedReads << " (" << currentTotalDeterminedReads * 100 /static_cast<double>(currentTotalReads)<< "%)"<< std::endl;
			std::cout << "\tTotal Less Than Small Len Cut Off Count: " << initialCounts.smallLenCutOffCount << " (" << initialCounts.smallLenCutOffCount * 100 /static_cast<double>(currentTotalReads)<< "%)"<< std::endl;
			std::cout << "\tTotal Contains Ns Count: " << initialCounts.containsNs << " (" << initialCounts.containsNs * 100 /static_cast<double>(currentTotalReads)<< "%)"<< std::endl;
			std::cout << "\tTotal Poor Quality Count: " << initialCounts.poorQualityCount << " (" << initialCounts.poorQualityCount * 100 /static_cast<double>(currentTotalReads)<< "%)"<< std::endl;
			std::cout << "\tTotal Inverse pair Count: " << initialCounts.genTotalInversePairCount() << " (" << initialCounts.genTotalInversePairCount() * 100 /static_cast<double>(currentTotalReads)<< "%)"<< std::endl;
			std::cout << "\tTotal Multihit Count: " << initialCounts.genTotalMultihitCount() << " (" << initialCounts.genTotalMultihitCount() * 100 /static_cast<double>(currentTotalReads)<< "%)"<< std::endl;
			std::cout << "\tDissimilar Filter: " << initialCounts.filteredDissimilarCount << " (" << initialCounts.filteredDissimilarCount * 100 /static_cast<double>(currentTotalReads)<< "%)"<< std::endl;
			std::cout << std::endl;
		}
		if(currentTotalUndeterminedReads == 0 || totalDeterminedReads == currentTotalDeterminedReads){
			//if there's no more undetermined reads or if the number determined this round was zero then stop iterating
			iterate = false;
		}
		totalDeterminedReads = currentTotalDeterminedReads;
		if (!keepTemporaryFiles) {
			if (bfs::exists(undeterminedPairedInOpts.firstName_)) {
				bfs::remove(undeterminedPairedInOpts.firstName_);
				bfs::remove(undeterminedPairedInOpts.secondName_);
			}
			if (bfs::exists(undeterminedSingleInOpts.firstName_)) {
				bfs::remove(undeterminedSingleInOpts.firstName_);
			}
		}
		++iterNumber;
	}

	if(!doNotDoFinalExtract){
		initialSeqOut.closeOutForReopeningAll();

		if(!extractingPars.doReCheckExcludeSets && doFinalRecheckOfExcludeSeqs){
			watch.startNewLap(njh::pasteAsStr(iterNumber, "- final read in new kmers"));
			auto rawUniqKmerSets = njh::getVecOfMapKeys(uniqueKmersPerSet);
			removeElement(rawUniqKmerSets, std::string("undetermined"));
			std::unordered_map<std::string, std::unordered_map<uint64_t, uint32_t>> rawKmersPerInput = UniqueKmerSetHelper::readInNewKmersFromExtractedReads(
							setUp.pars_.directoryName_, rawUniqKmerSets, extractingPars);


			watch.startNewLap(njh::pasteAsStr(iterNumber, "- final check new kmers against current set"));

			std::unordered_map<std::string, std::unordered_set<uint64_t>> outputUniqueKmersPerSet = UniqueKmerSetHelper::filterReExtractedKmersForNonUniqueIncludeExcludedSets(
							rawKmersPerInput,
							extractingPars,
							uniqueKmersPerSet,
							nonUniqueKmersPerSet);

			watch.startNewLap(njh::pasteAsStr(iterNumber, "- final copy new set"));
			uniqueKmersPerSet = std::move(outputUniqueKmersPerSet);
		}

		UniqueKmerSetHelper::ProcessReadForExtractingCounts finalCounts;


		watch.startNewLap(njh::pasteAsStr(iterNumber, "- scan whole set again"));


		auto finalExtractionDir = njh::files::make_path(setUp.pars_.directoryName_, "finalExtraction");
		njh::files::makeDir(njh::files::MkdirPar{finalExtractionDir});

		//add outputs
		MultiSeqIO finalSeqOut;

		VecStr names = getVectorOfMapKeys(uniqueKmersPerSet);
		for(const auto & name : names){
			auto seqOutOpts = SeqIOOptions::genPairedOutGz(njh::files::make_path(finalExtractionDir, name));
			finalSeqOut.addReader(njh::pasteAsStr(name, "-paired"), seqOutOpts);
		}
		finalSeqOut.addReader("undetermined-paired", SeqIOOptions::genPairedOutGz(njh::files::make_path(finalExtractionDir, "undetermined")));
		finalSeqOut.addReader("multihit-paired", SeqIOOptions::genPairedOutGz(njh::files::make_path(finalExtractionDir, "multihit")));

		for(const auto & name : names){
			auto seqOutOpts = SeqIOOptions::genFastqOutGz(njh::files::make_path(finalExtractionDir, name));
			finalSeqOut.addReader(njh::pasteAsStr(name, "-single"), seqOutOpts);
		}
		finalSeqOut.addReader("undetermined-single", SeqIOOptions::genFastqOutGz(njh::files::make_path(finalExtractionDir, "undetermined")));
		finalSeqOut.addReader("multihit-single", SeqIOOptions::genFastqOutGz(njh::files::make_path(finalExtractionDir, "multihit")));

		//setting final extraction cut offs
		extractingPars.compPars.initialExcludeHardCountOff = extractingPars.compPars.finalHardCountOff;
		extractingPars.compPars.initialExcludeFracCutOff = extractingPars.compPars.finalFracCutOff;

		extractingPars.compPars.hardCountOff = extractingPars.compPars.finalHardCountOff;
		extractingPars.compPars.fracCutOff = extractingPars.compPars.finalFracCutOff;

		//Build the strand-bias corrections before writing anything, using the same final cut offs the
		//write pass will use so the two passes classify identically. This re-reads the inputs, which is
		//affordable because it happens once, after the k-mer sets are frozen, and it is the only point
		//where each read's sequencing strand relative to its set is still knowable.
		if (sbPars.on_) {
			watch.startNewLap("build strand bias corrections");
			SimpleKmerHash sbHasher;
			extractingPars.strandBiasCorrectors = buildStrandBiasCorrectorsPerSet(
					setUp.pars_.ioOptions_, singlesOption, uniqueKmersPerSet, extractingPars, sbHasher,
					sbPars, strandBiasMaxReadsPerSet, std::cout);
			watch.startNewLap("final extraction");
		}

		if (finalWriteOutMultiHitSeparate) {
			extractingPars.writeMultiHitSeparate = true;
		}
		extractingPars.compPars.pairsSeparate = finalExtractionPairsSeparate;
		if (!setUp.pars_.ioOptions_.firstName_.empty()) {
			SeqInput reader(setUp.pars_.ioOptions_);
			reader.openIn();

			std::function<void()> readInComp;

			if(extractingPars.compPars.pairsSeparate){
				readInComp = [&reader, &uniqueKmersPerSet,
								&mut, &extractingPars, &finalSeqOut, &finalCounts
				]() {
					SimpleKmerHash hasher;
					PairedRead pseq;
					UniqueKmerSetHelper::ProcessReadForExtractingCounts currentCounts;
					while (reader.readNextReadLock(pseq)) {
						bool firstMatePassQcCheck = UniqueKmerSetHelper::readPassesQCChecks(pseq.seqBase_, extractingPars,currentCounts);
						bool secondMatePassQcCheck = UniqueKmerSetHelper::readPassesQCChecks(pseq.mateSeqBase_, extractingPars,currentCounts);
						if(firstMatePassQcCheck && secondMatePassQcCheck) {
							UniqueKmerSetHelper::processReadForExtractingPairsSeparate(pseq, uniqueKmersPerSet,
																																			 extractingPars, hasher, finalSeqOut,
																																			 currentCounts
																							);
						} else if(firstMatePassQcCheck) {
							pseq.seqBase_.name_.append("_firstMate");
							UniqueKmerSetHelper::processReadForExtracting(pseq.seqBase_, uniqueKmersPerSet,
																		extractingPars, hasher, finalSeqOut,
																		currentCounts);
						} else if(secondMatePassQcCheck) {
							pseq.mateSeqBase_.name_.append("_secondMate");
							UniqueKmerSetHelper::processReadForExtracting(pseq.mateSeqBase_, uniqueKmersPerSet,
																		extractingPars, hasher, finalSeqOut,
																		currentCounts);
						}
					}
					{
						std::lock_guard<std::mutex> lockGuard(mut);
						finalCounts.addOtherCounts(currentCounts);
					}
				};
			}else{
				readInComp = [&reader, &uniqueKmersPerSet,
								&mut, &extractingPars, &finalSeqOut, &finalCounts
				]() {
					SimpleKmerHash hasher;
					PairedRead pseq;
					UniqueKmerSetHelper::ProcessReadForExtractingCounts currentCounts;
					while (reader.readNextReadLock(pseq)) {
						bool firstMatePassQcCheck = UniqueKmerSetHelper::readPassesQCChecks(pseq.seqBase_, extractingPars,currentCounts);
						bool secondMatePassQcCheck = UniqueKmerSetHelper::readPassesQCChecks(pseq.mateSeqBase_, extractingPars,currentCounts);
						if(firstMatePassQcCheck && secondMatePassQcCheck) {
							UniqueKmerSetHelper::processReadForExtractingPairsTogether(pseq, uniqueKmersPerSet,
																																			 extractingPars, hasher, finalSeqOut,
																																			 currentCounts
																							);
						} else if(firstMatePassQcCheck) {
							pseq.seqBase_.name_.append("_firstMate");
							UniqueKmerSetHelper::processReadForExtracting(pseq.seqBase_, uniqueKmersPerSet,
																		extractingPars, hasher, finalSeqOut,
																		currentCounts);
						} else if(secondMatePassQcCheck) {
							pseq.mateSeqBase_.name_.append("_secondMate");
							UniqueKmerSetHelper::processReadForExtracting(pseq.mateSeqBase_, uniqueKmersPerSet,
																		extractingPars, hasher, finalSeqOut,
																		currentCounts);
						}
					}
					{
						std::lock_guard<std::mutex> lockGuard(mut);
						finalCounts.addOtherCounts(currentCounts);
					}
				};
			}
			njh::concurrent::runVoidFunctionThreaded(readInComp, numThreads);
		}
		if (!singlesOption.firstName_.empty()) {
			SeqInput reader(singlesOption);
			reader.openIn();
			std::function<void()> readInComp = [&reader, &uniqueKmersPerSet,
							&mut, &extractingPars, &finalSeqOut, &finalCounts
//							, &testReads
			]() {
				SimpleKmerHash hasher;
				seqInfo seq;
				UniqueKmerSetHelper::ProcessReadForExtractingCounts currentCounts;
				while (reader.readNextReadLock(seq)) {
					if(UniqueKmerSetHelper::readPassesQCChecks(seq,extractingPars,currentCounts)) {
						UniqueKmerSetHelper::processReadForExtracting(seq, uniqueKmersPerSet,
						                                              extractingPars, hasher, finalSeqOut,
						                                              currentCounts
						);
					}
				}
				{
					std::lock_guard<std::mutex> lockGuard(mut);
					finalCounts.addOtherCounts(currentCounts);
				}
			};
			njh::concurrent::runVoidFunctionThreaded(readInComp, numThreads);
		}
		finalCounts.writeOutCounts(outCounts,extractingPars,uniqueKmersPerSet,"final");
		{
			auto totalReadCnt = finalCounts.getTotalCounts();
			for(const auto & set : finalCounts.inversePairReadCountsPerSet) {
				inversePairOutCountsOut << njh::pasteAsStr("final")
				<< "\t" << extractingPars.compPars.sampleName
				<< "\t" << totalReadCnt
				<< "\t" << set.first
				<< "\t" << set.second
				<< "\t" << static_cast<double>(set.second)/static_cast<double>(totalReadCnt) << std::endl;
			}
		}

		{
			auto totalReadCnt = finalCounts.getTotalCounts();
			for(const auto & set : finalCounts.multiHitReadCounts) {
				multiHitOutCountsOut << njh::pasteAsStr("final")
				<< "\t" << extractingPars.compPars.sampleName
				<< "\t" << totalReadCnt
				<< "\t" << set.first
				<< "\t" << set.second
				<< "\t" << static_cast<double>(set.second)/static_cast<double>(totalReadCnt) << std::endl;
			}
		}
		if(setUp.pars_.verbose_){
			std::cout << "iterNumber: " << "final" << std::endl;
			std::cout << "Total Reads: " << finalCounts.getTotalCounts() << std::endl;
			std::cout << "Total Undetermined Reads: " << finalCounts.genTotalUndeterminedCount() << " (" << finalCounts.genTotalUndeterminedCount() * 100 /static_cast<double>(finalCounts.getTotalCounts())<< "%)" << std::endl;
			std::cout << "Total Determined Reads: " << finalCounts.genTotalDeterminedCount() << " (" << finalCounts.genTotalDeterminedCount() * 100 /static_cast<double>(finalCounts.getTotalCounts())<< "%)"<< std::endl;
			std::cout << "Total Less Than Small Len Cut Off Count: " << finalCounts.smallLenCutOffCount << " (" << finalCounts.smallLenCutOffCount * 100 /static_cast<double>(finalCounts.getTotalCounts())<< "%)"<< std::endl;
			std::cout << "Total Low Quality Count: " << finalCounts.poorQualityCount << " (" << finalCounts.poorQualityCount * 100 /static_cast<double>(finalCounts.getTotalCounts())<< "%)"<< std::endl;
			std::cout << "Total Inverse Pair Count: " << finalCounts.genTotalInversePairCount() << " (" << finalCounts.genTotalInversePairCount() * 100 /static_cast<double>(finalCounts.getTotalCounts())<< "%)"<< std::endl;
			std::cout << "Total Multihit Count: " << finalCounts.genTotalMultihitCount() << " (" << finalCounts.genTotalMultihitCount() * 100 /static_cast<double>(finalCounts.getTotalCounts())<< "%)"<< std::endl;

			if (setUp.pars_.debug_) {
				std::cout << watch.getLapName() << "\t" << watch.timeLapFormatted() << std::endl;
				std::cout << "klen: " << extractingPars.compPars.klen << std::endl;
				auto setNames = njh::getSetOfMapKeys(uniqueKmersPerSet);
				for (const auto &setName: setNames) {
					std::cout << setName << "\t" << uniqueKmersPerSet[setName].size() << std::endl;
				}
				std::cout << "NON_UNIQUE" << "\t" << nonUniqueKmersPerSet.size() << std::endl;
			}
		}
		if(!keepTemporaryFiles){
			//remove all fast[aq] files and fast[aq].gz files within directory
			auto readFiles = njh::files::listAllFiles(setUp.pars_.directoryName_,false,{std::regex{R"((.*)(\.(fast[aq])(\.gz)?)$)"}});
			for(const auto & f : readFiles){
				bfs::remove(f.first);
			}
		} else {
			//zip up any temp files if keeping
			auto readFiles = njh::files::listAllFiles(setUp.pars_.directoryName_,false,{std::regex{R"((.*)(\.(fast[aq]))$)"}});
			for(const auto & f : readFiles){
				njh::gzZipFile(njh::IoOptions{njh::InOptions(f.first), njh::OutOptions(bfs::path(f.first.string() + ".gz"))});
				bfs::remove(f.first);
			}
		}
	} else {
		initialSeqOut.closeOutForReopeningAll();
		auto readFiles = njh::files::listAllFiles(setUp.pars_.directoryName_,false,{std::regex{R"((.*)(\.(fast[aq]))$)"}});
		for(const auto & f : readFiles){
			njh::gzZipFile(njh::IoOptions{njh::InOptions(f.first), njh::OutOptions(bfs::path(f.first.string() + ".gz"))});
			bfs::remove(f.first);
		}
	}

	if(extractingPars.writeOutFinalKmerSets){
		SimpleKmerHash hasher;
		watch.startNewLap(njh::pasteAsStr("write Out Final Kmer Sets"));
		OutputStream finalUniqueKmerSetsOut(njh::files::make_path(setUp.pars_.directoryName_, "finalUniqueKmerSets.tsv.gz"));
		for(const auto & uniqKSet : uniqueKmersPerSet){
			for(const auto & uniqK : uniqKSet.second){
				finalUniqueKmerSetsOut << uniqKSet.first << "\t" << hasher.reverseHash(uniqK)  << "\n";
			}
		}
		OutputStream nonUniqueKmerSets(njh::files::make_path(setUp.pars_.directoryName_, "nonUniqueKmers_finalUniqueKmerSets.tsv.gz"));
		for(const auto & nonUniqK : nonUniqueKmersPerSet){
			nonUniqueKmerSets << nonUniqueRegionName << "\t" << hasher.reverseHash(nonUniqK) << "\n";
		}
	}
	if(!keepTemporaryFiles){
		if (bamFnpSet) {
			watch.startNewLap(njh::pasteAsStr("remove bam extracted files"));

			bfs::path bamExtractR1Fnp = bamExtractOut.outFilename_.string() + "_R1.fastq.gz";
			bfs::path bamExtractR2Fnp = bamExtractOut.outFilename_.string() + "_R2.fastq.gz";
			bfs::path bamExtractFnp = bamExtractOut.outFilename_.string() + ".fastq.gz";
			if (exists(bamExtractR1Fnp)) {
				bfs::remove(bamExtractR1Fnp);
				bfs::remove(bamExtractR2Fnp);
			}
			if (exists(bamExtractFnp)) {
				bfs::remove(bamExtractFnp);
			}
		}
	}
	watch.logLapTimes(setUp.rLog_.runLogFile_, true, 6, true);
	return 0;
}


}  // namespace njhseq

