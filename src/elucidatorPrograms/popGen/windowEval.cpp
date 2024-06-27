//
// Created by Nicholas Hathaway on 6/9/24.
//

/*
 * popGenExpRunner_windowEval.cpp
 *
 *  Created on: Jan 10, 2019
 *      Author: nicholashathaway
 */

#include "popGenExp.hpp"
#include <njhseq/GenomeUtils/GenomeMapping/MultiGenomeMapper.hpp>
#include <njhseq/BamToolsUtils.h>


namespace njhseq {

class SeqWindowEvaluator {
public:

	SeqWindowEvaluator(const MultiGenomeMapper::inputParameters & gmPars) :
			genomeMapper_(gmPars), bioCmdRunner_(gmPars.verbose_) {
		genomeMapper_.init();
	}

	MultiGenomeMapper genomeMapper_;

	bfs::path tempAlignmentsDir_;
	bfs::path fastaDir_;

	BioCmdsUtils bioCmdRunner_;

	bool extendAndTrim_ = false;
	uint32_t extendAndTrimLen_ = 10;

	void setUpWorkingDirectory(const std::string & fallBackDirectoryName) {
		if ("" == genomeMapper_.pars_.workingDirectory_.dirName_) {
			genomeMapper_.pars_.workingDirectory_.dirName_ = fallBackDirectoryName;
		}
		if (bfs::exists(genomeMapper_.pars_.workingDirectory_.dirName_)
				&& !genomeMapper_.pars_.workingDirectory_.overWriteDir_) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error outDir "
					<< genomeMapper_.pars_.workingDirectory_.dirName_
					<< " already exists, use --overWriteDir to over write it" << "\n";
			throw std::runtime_error { ss.str() };
		}
		genomeMapper_.pars_.workingDirectory_.overWriteDir_ = true;
		njh::files::makeDir(genomeMapper_.pars_.workingDirectory_);

		tempAlignmentsDir_ = njh::files::makeDir(
				genomeMapper_.pars_.workingDirectory_.dirName_,
				njh::files::MkdirPar("tempAlignments"));

		fastaDir_ = njh::files::makeDir(
				genomeMapper_.pars_.workingDirectory_.dirName_,
				njh::files::MkdirPar("fastas"));

	}

	class WinowResults {
	public:
		WinowResults() {
		}
		WinowResults(const GenomicRegion & region) :
				region_(region) {
		}
		GenomicRegion region_;

		double avgGcContent_ = 0;
		bool lengthVariation_ = false;
		uint32_t minLen_ = 0;
		uint32_t maxLen_ = 0;
		uint32_t uniqueVarNum_ = 0;
		uint32_t varTotalNum_ = 0;
		bool multiMapping_ = false;
		bool failedSeqEval_ = false;
		bool failedToHitAllGenomes_ = false;


		std::unordered_map<std::string, std::vector<GenomicRegion>> allRegions_;
		std::unordered_map<std::string, std::vector<std::shared_ptr<readObject>>> allSeqs_;
		std::shared_ptr<readObject> primaryExtract_;

		bool addRegions(const std::string & genome,
				const bfs::path & bamFnp) {
			bool addedRegion = false;
			BamTools::BamReader bReader;
			bReader.Open(bamFnp.string());
			checkBamOpenThrow(bReader, bamFnp);
			BamTools::BamAlignment bAln;
			auto refData = bReader.GetReferenceData();
			while (bReader.GetNextAlignmentCore(bAln)) {
				if (bAln.IsMapped()) {
					addRegion(genome,GenomicRegion(region_.uid_,
									refData[bAln.RefID].RefName, bAln.Position, bAln.GetEndPosition(),
									bAln.IsReverseStrand()));
					addedRegion = true;
				}
			}
			return addedRegion;
		}

		void addRegion(const std::string & genome, const GenomicRegion & region){
			allRegions_[genome].emplace_back(region);
		}

		void setPrimaryExtract(const MultiGenomeMapper & gMapper){
			//TwoBit::TwoBitFile tReader(gMapper.genomes_.at(gMapper.pars_.primaryGenome_)->fnpTwoBit_);
			TwoBit::TwoBitFile tReader(njh::mapAt(gMapper.genomes_, gMapper.pars_.primaryGenome_)->fnpTwoBit_);
			primaryExtract_ = std::make_shared<readObject>(region_.extractSeq(tReader));
			primaryExtract_->createCondensedSeq();
			primaryExtract_->setLetterCount();
			primaryExtract_->counter_.resetAlphabet(false);
		}

		void resetCounts(){
			avgGcContent_ = 0;
			lengthVariation_ = false;
			uniqueVarNum_ = 0;
			varTotalNum_ = 0;
		}

		void populateSeqs(const MultiGenomeMapper & gMapper){
			allSeqs_.clear();
			resetCounts();
			for(const auto & gRegions : allRegions_){
				gMapper.checkForGenomeThrow(gRegions.first, __PRETTY_FUNCTION__);
				TwoBit::TwoBitFile tReader(njh::mapAt(gMapper.genomes_,gRegions.first)->fnpTwoBit_);
				uint32_t extractionCount = 0;
				for(const auto & region : gRegions.second){
					MetaDataInName meta;
					meta.addMeta("refRegion", region.uid_);
					meta.addMeta("genome", gRegions.first);
					meta.addMeta("pos", region.start_);
					meta.addMeta("chrom", region.chrom_);
					meta.addMeta("extractionCount", extractionCount);
					std::shared_ptr<readObject> currentExtract = std::make_shared<readObject>(region.extractSeq(tReader));
					currentExtract->createCondensedSeq();
					currentExtract->setLetterCount();
					currentExtract->counter_.resetAlphabet(false);
					currentExtract->seqBase_.name_ = meta.createMetaName();
					allSeqs_[gRegions.first].emplace_back(currentExtract);
					++extractionCount;
				}
			}
			std::vector<seqInfo> uniqueSeqs;
			std::vector<size_t> readLengths;
			std::vector<double> gcContents;
			uint32_t totalSeqCount = 0;
			auto getName =[](const seqInfo & seq){
				MetaDataInName meta(seq.name_);

				std::string outName = meta.getMeta("genome");
				if("0" != meta.getMeta("extractionCount")){
					outName += "." + meta.getMeta("extractionCount");
				}
				return outName;
			};
			for (auto & ref : allSeqs_) {
				for(auto & seq : ref.second){
					++totalSeqCount;
					readLengths.emplace_back(len(*seq));
					seq->counter_.calcGcContent();
					gcContents.emplace_back(seq->counter_.gcContent_);
					bool found = false;
					for (auto & otherRef : uniqueSeqs) {
						if (otherRef.seq_ == seq->seqBase_.seq_) {
							otherRef.name_ += "-" + getName(seq->seqBase_);
							found = true;
							break;
						}
					}
					if (!found) {
						auto seqBaseSeed = seq->seqBase_;
						seqBaseSeed.name_ = getName(seq->seqBase_);
						uniqueSeqs.emplace_back(seq->seqBase_);
					}
				}
			}

			avgGcContent_ = vectorMean(gcContents);
			uniqueVarNum_ = uniqueSeqs.size();
			varTotalNum_ = totalSeqCount;
			lengthVariation_ = vectorMaximum(readLengths) - vectorMinimum(readLengths) >= 2;
			minLen_ = vectorMinimum(readLengths);
			maxLen_ = vectorMaximum(readLengths);
		}

		void write(std::ostream & out, uint32_t totalInputGenomes) {
			std::vector<uint32_t> numberOfHits;
			for(const auto & genome : allRegions_){
				numberOfHits.emplace_back(genome.second.size());
			}
			auto gKeys = getVectorOfMapKeys(allRegions_);
			std::set<std::string> genomesIn(gKeys.begin(), gKeys.end());

			out << region_.genBedRecordCore().toDelimStr()
					<< "\t" << avgGcContent_
					<< "\t" << njh::boolToStr(lengthVariation_)
					<< "\t" << minLen_
					<< "\t" << maxLen_
					<< "\t" << uniqueVarNum_
					<< "\t" << allRegions_.size() << "\t" << totalInputGenomes
					<< "\t" << vectorMean(numberOfHits)
					<< "\t" << njh::conToStr(genomesIn, ",")
					<< std::endl;
		}
	};

	void writeOutResultSeqs(WinowResults & result) const{
		if(result.allSeqs_.empty()){
			result.populateSeqs(genomeMapper_);
		}
		auto allSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(fastaDir_, result.region_.uid_ + "_allSeqs"));
		SeqOutput writer(allSeqOpts);
		writer.openOut();
		for(const auto & gSeqs : result.allSeqs_){
			writer.writeNoCheck(gSeqs.second);
		}
		writer.closeOut();
	}

	std::shared_ptr<WinowResults> evaluateForUniqueWindow(const GenomicRegion & window,
			const BioCmdsUtils::LastZPars & lzPars,
			const std::function<bool(const std::shared_ptr<readObject> &)> & evalFunc,
			aligner & alignerObj) const {


		std::shared_ptr<WinowResults> results = std::make_shared<WinowResults>(window);


		results->setPrimaryExtract(genomeMapper_);
		if(evalFunc(results->primaryExtract_)){
			results->failedSeqEval_ = true;
		}


		if(!results->failedSeqEval_){
			auto currentRegionsTempDir = njh::files::makeDir(tempAlignmentsDir_,
					njh::files::MkdirPar(window.createUidFromCoords()));

			auto primSeqOpts = SeqIOOptions::genFastaOut(
					njh::files::make_path(currentRegionsTempDir,
							genomeMapper_.pars_.primaryGenome_));
			SeqOutput::write(std::vector<readObject> { *(results->primaryExtract_) }, primSeqOpts);


			for (const auto & genome : genomeMapper_.genomes_) {
				auto lastzCurrentPars = lzPars;
				lastzCurrentPars.genomeFnp = njh::mapAt(genomeMapper_.genomes_,genome.first)->fnpTwoBit_;
				auto lastzSeqPars = SeqIOOptions::genFastaIn(primSeqOpts.getPriamryOutName());
				lastzSeqPars.out_ = OutOptions(njh::files::make_path(currentRegionsTempDir, genome.first + ".bam"));
				auto res = bioCmdRunner_.lastzAlignNoSort(lastzSeqPars, lastzCurrentPars);
				if(bioCmdRunner_.verbose_){
					std::cout << res.toJson() << std::endl;
				}
				bool added = results->addRegions(genome.first,lastzSeqPars.out_.outFilename_);
				if(added && extendAndTrim_){
					auto extenedRegions = results->allRegions_.at(genome.first);
					for(auto & reg : results->allRegions_.at(genome.first)){
						auto extenedRegion = reg;
						extenedRegion.start_ = reg.start_ <= extendAndTrimLen_? 0 : reg.start_ - extendAndTrimLen_;
						extenedRegion.end_ = reg.end_ + extendAndTrimLen_ < genome.second->chromosomeLengths_.at(reg.chrom_) ? reg.end_ + extendAndTrimLen_ : genome.second->chromosomeLengths_.at(reg.chrom_);
						TwoBit::TwoBitFile tReader(genome.second->fnpTwoBit_);
						auto extractedSeq = extenedRegion.extractSeq(tReader);
						auto trimmedExtractedSeq = extractedSeq;
						readVecTrimmer::GlobalAlnTrimPars trimPars{};
						trimPars.startInclusive_ = 0;
						trimPars.endInclusive_ = len(*results->primaryExtract_) -1 ;
						if(len(trimmedExtractedSeq) >= alignerObj.parts_.maxSize_){
							alignerObj.parts_.setMaxSize(len(trimmedExtractedSeq));
						}
						readVecTrimmer::trimSeqToRefByGlobalAln(trimmedExtractedSeq, *results->primaryExtract_, trimPars, alignerObj);
						if(trimmedExtractedSeq.on_){
							uint32_t startPos = extractedSeq.seq_.find(trimmedExtractedSeq.seq_);
							uint32_t stopPos = startPos + len(trimmedExtractedSeq);
							uint32_t trimmedOffBack = len(extractedSeq.seq_) - stopPos;
							if(reg.reverseSrand_){
								reg.start_ = extenedRegion.start_ + trimmedOffBack;
								reg.end_ = extenedRegion.end_ - startPos;
							}else{
								reg.start_ = extenedRegion.start_ + startPos;
								reg.end_ = extenedRegion.end_ - trimmedOffBack;
							}
						}
					}
				}
			}
			//check for uniqueness
			for(const auto & regions : results->allRegions_){
				if(regions.second.size() > 1){
					results->multiMapping_ = true;
					break;
				}
			}
			//if each maps uniquely check each extracted seq
			if(!results->multiMapping_){
				results->populateSeqs(genomeMapper_);
				for(const auto & genome : results->allSeqs_){
					if(std::any_of(genome.second.begin(), genome.second.end(), evalFunc)){
						results->failedSeqEval_ = true;
						break;
					}
				}
			}
			if (!genomeMapper_.pars_.keepTempFiles_) {
				njh::files::rmDirForce(currentRegionsTempDir);
			}
		}
		return results;
	}

	std::shared_ptr<WinowResults> collectInfoOnWindow(const GenomicRegion & window,
			const BioCmdsUtils::LastZPars & lzPars) const {
		std::shared_ptr<WinowResults> results = std::make_shared<WinowResults>(window);
		auto currentRegionsTempDir = njh::files::makeDir(tempAlignmentsDir_,
				njh::files::MkdirPar(window.createUidFromCoords()));
		results->setPrimaryExtract(genomeMapper_);
		auto primSeqOpts = SeqIOOptions::genFastaOut(
				njh::files::make_path(currentRegionsTempDir,
						genomeMapper_.pars_.primaryGenome_));
		SeqOutput::write(std::vector<readObject> { *(results->primaryExtract_) }, primSeqOpts);

		for (const auto & genome : genomeMapper_.genomes_) {
			auto lastzCurrentPars = lzPars;
			lastzCurrentPars.genomeFnp = njh::mapAt(genomeMapper_.genomes_,genome.first)->fnpTwoBit_;
			auto lastzSeqPars = SeqIOOptions::genFastaIn(primSeqOpts.getPriamryOutName());
			lastzSeqPars.out_ = OutOptions(njh::files::make_path(currentRegionsTempDir, genome.first + ".sorted.bam"));
			auto res = bioCmdRunner_.lastzAlign(lastzSeqPars, lastzCurrentPars);
			if(bioCmdRunner_.verbose_){
				std::cout << res.toJson() << std::endl;
			}
			results->addRegions(genome.first,lastzSeqPars.out_.outFilename_);
		}
		results->populateSeqs(genomeMapper_);
		if (!genomeMapper_.pars_.keepTempFiles_) {
			njh::files::rmDirForce(currentRegionsTempDir);
		}
		return results;
	}



};



int popGenExpRunner::variationWindowShoppingWithLastz(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bedfile = "";
	bool noSeqBasedFilter = false;
	BioCmdsUtils::LastZPars lzPars;
	lzPars.coverage = 100;
	uint32_t homopolymerCutOff = 15;
	bool allGenomes = false;
	bool extendAndTrim = false;
	bool writeFasta = false;
	uint32_t extendAndTrimLen = 10;
	MultiGenomeMapper::inputParameters gmPars;
	std::string selectedGenomesStr;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(writeFasta, "--writeFasta", "write out a fasta for the final results");
	setUp.setOption(extendAndTrim, "--extendAndTrim", "Extend the determine region and then trim back, can be helpful for when variation falls at the very ends of the sequence");
	setUp.setOption(extendAndTrimLen, "--extendAndTrimLen", "When extending and trimming, use this length");
	setUp.setOption(noSeqBasedFilter, "--noSeqBasedFilter", "No Seq Based Filter");
	setUp.setOption(allGenomes, "--allGenomes", "Must appear in all genomes to be reported on");
	setUp.setOption(gmPars.workingDirectory_.dirName_, "--outDir", "Output directory");
	setUp.setOption(gmPars.keepTempFiles_, "--keepTempFiles", "keep Temp Files");
	setUp.setOption(gmPars.workingDirectory_.overWriteDir_ , "--overWriteDir", "Whether or not to overwrite outDir");
	setUp.setOption(homopolymerCutOff, "--homopolymerCutOff", "homopolymer length Cut Off (inclusive)");
	setUp.setOption(gmPars.numThreads_, "--numThreads", "Number of threads to use");
	setUp.setOption(lzPars.coverage, "--coverage", "Coverage for lastz");
	setUp.setOption(lzPars.identity, "--identity", "Identity for lastz");
	setUp.setOption(lzPars.extraLastzArgs, "--extraLastzArgs", "Extra Lastz Arguments");
	setUp.setOption(bedfile, "--bed", "Input bed to search with", true);
	setUp.setOption(gmPars.genomeDir_, "--genomeDir", "Name of the genome file fnp", true);
	setUp.setOption(gmPars.primaryGenome_, "--primaryGenome", "Name of the primary genome, should be in --genomeDir", true);
	setUp.setOption(selectedGenomesStr, "--selectedGenomes", "Name of the other genomes in --genomeDir to be read in, leave blank to just do all fastas");
	gmPars.selectedGenomes_ = njh::strToSet<std::string>(selectedGenomesStr, ",");
	setUp.finishSetUp(std::cout);
	gmPars.verbose_ = setUp.pars_.verbose_;

	std::vector<char> allowableLetters{'A', 'C', 'G', 'T', 'a', 'c', 'g', 't'};

	auto regions = bedPtrsToGenomicRegs(getBeds(bedfile));
	uint32_t maxlen = 0;
	for(const auto & region : regions){
		if(region.getLen() > maxlen){
			maxlen = region.getLen();
		}
	}

	if(setUp.pars_.debug_){
		std::cout << "Read in " << regions.size() << " regions" << std::endl;
	}

	SeqWindowEvaluator windowEval(gmPars);
	windowEval.extendAndTrimLen_ = extendAndTrimLen;
	windowEval.extendAndTrim_ = extendAndTrim;
	windowEval.setUpWorkingDirectory( bfs::basename(bedfile) + "_variationHuntWithLastz_" + njh::getCurrentDate());
	bfs::copy_file(bfs::absolute(bedfile), njh::files::make_path(windowEval.genomeMapper_.pars_.workingDirectory_.dirName_, "inputRegions.bed"));


	setUp.startARunLog(gmPars.workingDirectory_.dirName_.string());

	auto bedsDir = njh::files::makeDir(gmPars.workingDirectory_.dirName_, njh::files::MkdirPar("beds"));
	auto resultsDir = njh::files::makeDir(gmPars.workingDirectory_.dirName_, njh::files::MkdirPar("results"));

	std::function<bool(const std::shared_ptr<readObject> &)> seqEvalFunc = [&allowableLetters,&homopolymerCutOff](const std::shared_ptr<readObject> & seq){
		bool firstCheckFail = false;
		if (std::any_of(seq->condensedSeqCount.begin(),
				seq->condensedSeqCount.end(), [&homopolymerCutOff](int hpRunLen) {
					return hpRunLen >= homopolymerCutOff;
				})) {
			firstCheckFail = true;
		}
		if(!firstCheckFail){
			for(const auto & base : seq->counter_.alphabet_){
				if(!njh::in(base, allowableLetters)){
					firstCheckFail = true;
					break;
				}
			}
		}
		return firstCheckFail;
	};

	if(noSeqBasedFilter){
		seqEvalFunc = [](const std::shared_ptr<readObject> & seq){
			return false;
		};
	}
	OutOptions resultsOpts(njh::files::make_path(resultsDir, "results.tab.txt"));
	std::ofstream resultsFile;
	resultsOpts.openFile(resultsFile);
	resultsFile << "#chrom\tstart\tend\tname\tlength\tstrand\tavgGcContent\tlengthVariation\tminLen\tmaxLen\tuniqueVariantNumber\tgenomeCount\ttotalPossibleGenomeCount\tavgHits\tgenomesIn" << std::endl;
	OutOptions failedFileOpts(njh::files::make_path(resultsDir, "failed.tab.txt"));
	OutputStream failedFileOut(failedFileOpts);
	failedFileOut << "#chrom\tstart\tend\tname\tlength\tstrand\tmultiMapping\tfailedtoHitAllGenomes\tfailedSeqEval\tgenomeCount\ttotalPossibleGenomeCount\thitRange\tgenomesIn" << std::endl;
	auto genomes = getVectorOfMapKeys(windowEval.genomeMapper_.genomes_);

	OutOptions genomesOpts(njh::files::make_path(gmPars.workingDirectory_.dirName_, "genomesUsed.txt"));
	std::ofstream genomesFile;
	genomesOpts.openFile(genomesFile);
	printVector(genomes, "\n", genomesFile);
	genomesFile.flush();
	njh::concurrent::LockableVec<GenomicRegion> regionsQueue(regions);
	OutOptions allBedOpts(njh::files::make_path(bedsDir, "all.bed"));
	OutputStream allBedfile(allBedOpts);

	std::mutex resultsMut;
	std::vector<SeqWindowEvaluator::WinowResults> results;

	auto gapPars = gapScoringParameters(5, 1, 0, 0, 0, 0);
	auto scoring = substituteMatrix(2, -2);
	uint64_t maxAlignLen = (maxlen + extendAndTrimLen) * 2;

	std::vector<aligner> aligners;
	for(uint32_t t = 0; t < gmPars.numThreads_; ++t){
		aligners.emplace_back(aligner(maxAlignLen, gapPars, scoring));
	}


	auto evaluateRegion = [&resultsMut,&regionsQueue,
												 &seqEvalFunc,&lzPars,&resultsFile,
												 &allBedfile,&allGenomes,&aligners,
												 &failedFileOut,&writeFasta](const SeqWindowEvaluator & windowEval, uint32_t threadNumber){
		GenomicRegion currentRegion;

		while(regionsQueue.getVal(currentRegion)){

			auto result = windowEval.evaluateForUniqueWindow(currentRegion, lzPars, seqEvalFunc, aligners[threadNumber]);
			//will only keep variatns that happend only once so can check to see if var number matches genomes count to see if all where hit
			if(allGenomes && result->varTotalNum_ != windowEval.genomeMapper_.genomes_.size()){
				result->failedToHitAllGenomes_ = true;
			}

			if(result->multiMapping_ || result->failedToHitAllGenomes_ || result->failedSeqEval_){
				std::lock_guard<std::mutex> lock(resultsMut);
				std::vector<uint32_t> genomeHits;

				for(const auto & reg : result->allRegions_){
					genomeHits.emplace_back(reg.second.size());
				}
				uint32_t minHits = 0;
				uint32_t maxHits = 0;
				if(!genomeHits.empty()){
					minHits = *std::min_element(genomeHits.begin(), genomeHits.end());
					maxHits = *std::max_element(genomeHits.begin(), genomeHits.end());
				}
				std::vector<std::string> genomesIn = getVectorOfMapKeys(result->allRegions_);
				njh::sort(genomesIn);
				failedFileOut << currentRegion.genBedRecordCore().toDelimStr()
						<< "\t" << njh::boolToStr(result->multiMapping_)
						<< "\t" << njh::boolToStr(genomesIn.size() != windowEval.genomeMapper_.genomes_.size())
						<< "\t" << njh::boolToStr(result->failedSeqEval_)
						<< "\t" << genomesIn.size()
						<< "\t" << windowEval.genomeMapper_.genomes_.size()
						<< "\t" << minHits << "-" << maxHits
						<< "\t" << njh::conToStr(genomesIn, ",")
				    << std::endl;
				continue;
			}
			{
				std::lock_guard<std::mutex> lock(resultsMut);
				for(auto & genome : result->allRegions_) {
					for(const auto & otherRegion : genome.second){
						allBedfile << otherRegion.genBedRecordCore().toDelimStr() << std::endl;
					}
				}
				if(writeFasta){
					windowEval.writeOutResultSeqs(*result);
				}
				result->write(resultsFile, windowEval.genomeMapper_.genomes_.size());
			}
		}
	};
	std::vector<std::thread> threads;

	for(uint32_t t = 0; t < gmPars.numThreads_; ++t){
		threads.emplace_back(std::thread(evaluateRegion, std::cref(windowEval), t));
	}

	for(auto & t : threads){
		t.join();
	}

	return 0;
}

int popGenExpRunner::collectWindowInfoWithLastz(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bedfile = "";

	BioCmdsUtils::LastZPars lzPars;
	lzPars.coverage = 100;
	MultiGenomeMapper::inputParameters gmPars;
	std::string selectedGenomesStr = "";
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(gmPars.workingDirectory_.dirName_, "--outDir", "Output directory");
	setUp.setOption(gmPars.keepTempFiles_, "--keepTempFiles", "keep Temp Files");
	setUp.setOption(gmPars.workingDirectory_.overWriteDir_ , "--overWriteDir", "Whether or not to overwrite outDir");
	setUp.setOption(gmPars.numThreads_, "--numThreads", "Number of threads to use");
	setUp.setOption(lzPars.coverage, "--coverage", "Coverage for lastz");
	setUp.setOption(lzPars.identity, "--identity", "Identity for lastz");
	setUp.setOption(lzPars.extraLastzArgs, "--extraLastzArgs", "Extra Lastz Arguments");
	setUp.setOption(bedfile, "--bed", "Input bed to search with", true);
	setUp.setOption(gmPars.genomeDir_, "--genomeDir", "Name of the genome file fnp", true);
	setUp.setOption(gmPars.primaryGenome_, "--primaryGenome", "Name of the primary genome, should be in --genomeDir", true);
	setUp.setOption(selectedGenomesStr, "--selectedGenomes", "Name of the other genomes in --genomeDir to be read in, leave blank to just do all fastas");
	gmPars.selectedGenomes_ = njh::strToSet<std::string>(selectedGenomesStr, ",");
	setUp.finishSetUp(std::cout);
	gmPars.verbose_ = setUp.pars_.verbose_;

	auto regions = bedPtrsToGenomicRegs(getBeds(bedfile));

	SeqWindowEvaluator windowEval(gmPars);
	windowEval.setUpWorkingDirectory( bfs::basename(bedfile) + "_variationHuntWithLastz_" + njh::getCurrentDate());


	setUp.startARunLog(gmPars.workingDirectory_.dirName_.string());

	auto bedsDir = njh::files::makeDir(gmPars.workingDirectory_.dirName_, njh::files::MkdirPar("beds"));
	auto resultsDir = njh::files::makeDir(gmPars.workingDirectory_.dirName_, njh::files::MkdirPar("results"));

	OutOptions resultsOpts(njh::files::make_path(resultsDir, "results.tab.txt"));
	std::ofstream resultsFile;
	resultsOpts.openFile(resultsFile);
	resultsFile << "chrom\tstart\tend\tname\tlength\tstrand\tavgGcContent\tlengthVariation\tuniqueVariantNumber\tgenomeCount\ttotalPossibleGenomeCount\tavgHits\tgenomesIn" << std::endl;
	auto genomes = getVectorOfMapKeys(windowEval.genomeMapper_.genomes_);

	OutOptions genomesOpts(njh::files::make_path(gmPars.workingDirectory_.dirName_, "genomesUsed.txt"));
	std::ofstream genomesFile;
	genomesOpts.openFile(genomesFile);
	printVector(genomes, "\n", genomesFile);
	genomesFile.flush();
	njh::concurrent::LockableVec<GenomicRegion> regionsQueue(regions);
	OutOptions allBedOpts(njh::files::make_path(bedsDir, "all.bed"));
	auto allBedfile = allBedOpts.openFile();

	std::mutex resultsMut;
	auto evaluateRegion = [&resultsMut,&regionsQueue,&lzPars,&resultsFile,&allBedfile](const SeqWindowEvaluator & windowEval){
		GenomicRegion currentRegion;
		while(regionsQueue.getVal(currentRegion)){
			auto result = windowEval.collectInfoOnWindow(currentRegion,lzPars);
			{
				std::lock_guard<std::mutex> lock(resultsMut);

				for(auto & genome : result->allRegions_) {
					for(const auto & otherRegion : genome.second){
						(*allBedfile) << otherRegion.genBedRecordCore().toDelimStr() << std::endl;
					}
				}
				windowEval.writeOutResultSeqs(*result);
				result->write(resultsFile, windowEval.genomeMapper_.genomes_.size());
			}
		}
	};
	std::vector<std::thread> threads;
	for(uint32_t t = 0; t < gmPars.numThreads_; ++t){
		threads.emplace_back(std::thread(evaluateRegion, std::cref(windowEval)));
	}
	for(auto & t : threads){
		t.join();
	}

	return 0;
}




}  // namespace njhseq
