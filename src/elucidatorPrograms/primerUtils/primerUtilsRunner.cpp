//
// Created by Nicholas Hathaway on 5/22/23.
//

#include "primerUtilsRunner.hpp"
#include <njhseq/IO/SeqIO/SeqInput.hpp>
#include <njhseq/IO/OutputStream.hpp>
#include <njhseq/PrimerIDUtils/PrimerDimerUtils.hpp>
#include <SeekDeep/objects/TarAmpSetupUtils/PrimersAndMids.hpp>
#include <TwoBit/IO/TwoBitFile.hpp>
#include <njhseq/objects/BioDataObject/BLASTHitTabular.hpp>
#include <njhseq/objects/BioDataObject/BioDataFileIO.hpp>
#include <njhseq/BamToolsUtils/ReAlignedSeq.hpp>
#include <njhseq/GenomeUtils/GenomeExtraction/ParsingAlignmentInfo/GenomeExtractResult.hpp>

namespace njhseq {



primerUtilsRunner::primerUtilsRunner()
				: njh::progutils::ProgramRunner(
				{
								addFunc("computeDimerizationScore", computeDimerizationScore, false),
								addFunc("testWithBlastForUnspecificAmplification", testWithBlastForUnspecificAmplification, false),
								addFunc("creatingMultiplexAmpliconPools", creatingMultiplexAmpliconPools, false),
								addFunc("testForPrimerDimers", testForPrimerDimers, false),
								//
				},
				"primerUtils") {}

int primerUtilsRunner::testWithBlastForUnspecificAmplification(
				const njh::progutils::CmdArgs &inputCommands) {
	uint32_t numThreads = 1;
	uint32_t minLen = 12;
	uint32_t errorAllowed = 0;
	uint32_t maxTargetSize = 6000;
	uint32_t minTargetSize = 100;
	bfs::path primersFnp;
	bfs::path genomeFnp;
	seqSetUp setUp(inputCommands);
	setUp.setOption(numThreads, "--numThreads", "number of Threads");
	setUp.setOption(minLen, "--minLen", "min length of the 3` end of the alignment, (inclusive)");
	setUp.setOption(errorAllowed, "--errorAllowed", "errors allowed");
	setUp.setOption(maxTargetSize, "--maxTargetSize", "max target size");
	setUp.setOption(minTargetSize, "--minTargetSize", "min target size");

	setUp.setOption(primersFnp, "--primersFnp", "Primers", true);
	setUp.setOption(genomeFnp, "--genomeFnp", "genome", true);
	setUp.processDirectoryOutputName(bfs::basename(primersFnp) + std::string("_testUnspecific_TODAY"), true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	PrimersAndMids ids(primersFnp);

	if(ids.getTargets().empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error in reading in target primers file " << primersFnp << "\n";
		ss << "Make sure there is a line that starts with target in file" << "\n";
		throw std::runtime_error{ss.str()};
	}
	ids.initPrimerDeterminator();

	bfs::path genome2bitFnp = bfs::path(genomeFnp).replace_extension(".2bit");

	TwoBit::TwoBitFile genomeReader(genome2bitFnp.string());

	std::vector<seqInfo> primerSeqs;
	std::unordered_map<std::string, std::string> forPrimerNameToTargetName;
	std::unordered_map<std::string, std::string> revPrimerNameToTargetName;
	std::unordered_map<std::string, std::string> primerNameToTargetName;

	for(const auto & primerInfo : ids.pDeterminator_->primers_) {
		uint32_t forwardCount = 0;
		for (const auto &fwd: iter::enumerate(primerInfo.second.fwds_)) {
			auto degens = createDegenStrs(fwd.element.info_.seq_);
			for (const auto &degen: degens) {
				auto outName = njh::pasteAsStr(primerInfo.first, "_F", ".", forwardCount);
				primerSeqs.emplace_back(outName, degen);
				forPrimerNameToTargetName[outName] = primerInfo.first;
				primerNameToTargetName[outName] = primerInfo.first;
			}
			++forwardCount;
		}
		uint32_t reverseCount = 0;
		for (const auto &rev: primerInfo.second.revs_) {
			auto degens = createDegenStrs(rev.info_.seq_);
			for (const auto &degen: degens) {
				auto outName = njh::pasteAsStr(primerInfo.first, "_R", ".", reverseCount);
				primerSeqs.emplace_back(outName, degen);
				revPrimerNameToTargetName[outName] = primerInfo.first;
				primerNameToTargetName[outName] = primerInfo.first;
			}
			++reverseCount;
		}
	}
	std::unordered_map<std::string, uint32_t> primerLengths;
	std::unordered_map<std::string, std::string> primerByName;

	std::vector<seqInfo> primerUniqueSeqs;
	readVecSorter::sortBySeq(primerSeqs, false);
	std::unordered_map<std::string, std::set<std::string>> primerToAllNames;
	for(const auto & primer : primerSeqs){
		primerToAllNames[primer.seq_].emplace(primer.name_);
		if(primerUniqueSeqs.empty()) {
			primerUniqueSeqs.emplace_back(primer);
		} else if(primerUniqueSeqs.back().seq_ != primer.seq_){
			primerUniqueSeqs.emplace_back(primer);
		}
	}

	auto allPrimerFnp = njh::files::make_path(setUp.pars_.directoryName_, "allPrimers.fasta");
	SeqOutput::write(primerUniqueSeqs, SeqIOOptions::genFastaOut(allPrimerFnp));
	for(const auto & seq : primerSeqs){
		primerLengths[seq.name_] = len(seq);
		primerByName[seq.name_] = seq.seq_;
	}
	auto allPrimerBlastHitsArchiveFnp = njh::files::make_path(setUp.pars_.directoryName_, "allPrimersBlast.archive");
	auto allPrimerBlastHitsTableFnp = njh::files::make_path(setUp.pars_.directoryName_, "allPrimersBlast.tsv");

	//_blastdb
	//.nsq
	auto blastDatabaseFnp = bfs::path(njh::rstripRet(bfs::path(genomeFnp).replace_extension("").string(), '.') + "_blastdb");
	std::string blastCmd = "blastn -query " + allPrimerFnp.string() + " -num_threads " + njh::pasteAsStr(numThreads) + " -task blastn-short -word_size 6 -db " + blastDatabaseFnp.string() + " -outfmt 6 -max_target_seqs 100000000 ";
	if(0 == errorAllowed){
		blastCmd += " -perc_identity 100 ";
	}
	blastCmd += " > " + allPrimerBlastHitsTableFnp.string();
	auto blastsCmdsOutput = njh::sys::run(VecStr{blastCmd});
	if(!blastsCmdsOutput.success_){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "error in running: " << blastsCmdsOutput.cmd_ << "\n";
		ss << "stderr: " << blastsCmdsOutput.stdErr_ << "\n";
		throw std::runtime_error{ss.str()};
	}
	std::vector<BLASTHitTab> blastHits;
	std::vector<std::shared_ptr<ReAlignedSeq>> results;
	uint64_t maxLen = 0;
	readVec::getMaxLength(primerSeqs, maxLen);
	ReAlignedSeq::genRealignmentPars reAlnPars;
	reAlnPars.extendAmount = 0;
//	reAlnPars.adjustForSoftClipping = true;
	reAlnPars.adjustForSoftClipping = false;


	aligner alignerObj(maxLen, gapScoringParameters(5, 1, 0, 0, 0, 0, 0, 0, 0, 0),
										 substituteMatrix::createDegenScoreMatrixCaseInsensitive(2, -2));
	auto chromLens = genomeReader.getSeqLens();
	comparison allowableErrors;
	allowableErrors.hqMismatches_ = errorAllowed;
	OutputStream passingBlastHitsInfos(njh::files::make_path(setUp.pars_.directoryName_, "passingBlastHitsInfos.bed"));
	{
		BioDataFileIO<BLASTHitTab> reader{IoOptions(InOptions(allPrimerBlastHitsTableFnp))};
		reader.openIn();
		BLASTHitTab hit;
		while(reader.readNextRecord(hit)){
			if(hit.alignLen_ < minLen){
				continue;
			}
			//check if the 3` end matches, note if even allowing for mismatches, can't just check if the end equals, have to check if adding the amount of error is equal to or more than the expected end length
//			std::cout << hit.queryName_ << std::endl;
//			std::cout << (hit.qEnd_ + errorAllowed) << std::endl;
//			std::cout << njh::mapAt(primerLengths, hit.queryName_) << std::endl;
//			std::cout << "hit.qEnd_ + errorAllowed) < njh::mapAt(primerLengths, hit.queryName_): " << njh::colorBool((hit.qEnd_ + errorAllowed) < njh::mapAt(primerLengths, hit.queryName_)) << std::endl << std::endl;
			if((hit.qEnd_ + errorAllowed) < njh::mapAt(primerLengths, hit.queryName_)){
				continue;
			}
			//check if the 3` end, if in plus strand, check the end, if in the reverseStrand first base
//			if(!hit.reverseStrand()){
//				if(hit.qEnd_ != njh::mapAt(primerLengths, hit.queryName_)){
//					continue;
//				}
//			} else {
//				if(hit.qStart_ > 1){
//					continue;
//				}
//			}

			blastHits.emplace_back(hit);
			//adjust for re-alignment

			if(allowableErrors.hqMismatches_ > 0){
				if(hit.qEnd_ < njh::mapAt(primerLengths, hit.queryName_)){
					auto endDiff = njh::mapAt(primerLengths, hit.queryName_) - hit.qEnd_;
					if (hit.reverseStrand() && hit.sStart_ > endDiff) {
						hit.sEnd_ -= endDiff;
					} else if (hit.sEnd_ + endDiff <= chromLens.at(hit.subjectName_)) {
						hit.sEnd_ += endDiff;
					}
				}
			}


			auto realignedSeq = ReAlignedSeq::genRealignment(hit, njh::mapAt(primerByName, hit.queryName_), alignerObj,
																											 chromLens, genomeReader, reAlnPars);

//			std::cout << "passed: " << njh::colorBool(allowableErrors.passErrorProfile(realignedSeq.comp_)) << std::endl;
			if(allowableErrors.passErrorProfile(realignedSeq.comp_)){
//				std::cout << __FILE__ << " " << __LINE__ << std::endl;
//				realignedSeq.alnRefSeq_.outPutSeqAnsi(std::cout);
//				realignedSeq.alnQuerySeq_.outPutSeqAnsi(std::cout);
//				std::cout << __FILE__ << " " << __LINE__ << std::endl;
				//std::cout << hit.toJson() << std::endl;
				for(const auto & primerName : primerToAllNames[primerByName.at(hit.queryName_)]){
					auto hitCopy = hit;
					auto realignedSeqCopy = realignedSeq;
					hitCopy.queryName_ = primerName;
					passingBlastHitsInfos << hitCopy.genSubjectBed6().toDelimStrWithExtra() << std::endl;
					realignedSeqCopy.gRegion_.uid_ = primerName;
					realignedSeqCopy.querySeq_.name_ = primerName;
					realignedSeqCopy.alnQuerySeq_.name_ = primerName;
					results.emplace_back(std::make_shared<ReAlignedSeq>(realignedSeqCopy));
				}
			}
		}
	}


//	exit(1);
	struct GenomeExtractResultWithHits {
		GenomeExtractResultWithHits(const std::shared_ptr<ReAlignedSeq> &p1, const std::shared_ptr<ReAlignedSeq> &p2) : p1_(
						p1), p2_(p2), extraction_(p1, p2) {
			extraction_.setRegion();
		}
		std::shared_ptr<ReAlignedSeq> p1_;
		std::shared_ptr<ReAlignedSeq> p2_;
		GenomeExtractResult extraction_;

		[[nodiscard]] seqInfo getP1Portion() const {
			auto alnStart = p1_->alnRefSeq_.seq_.find_first_not_of('-');
			auto alnStop = p1_->alnRefSeq_.seq_.find_last_not_of('-') + 1;

			auto portion = p1_->alnQuerySeq_.getSubRead(alnStart, alnStop - alnStart);
			portion.removeGaps();
			return portion;
		}

		[[nodiscard]] seqInfo getP2Portion() const {
			auto alnStart = p2_->alnRefSeq_.seq_.find_first_not_of('-');
			auto alnStop = p2_->alnRefSeq_.seq_.find_last_not_of('-') + 1;

			auto portion = p2_->alnQuerySeq_.getSubRead(alnStart, alnStop - alnStart);
			portion.removeGaps();
			return portion;
		}

		[[nodiscard]] seqInfo getP1PortionOriginal53() const {
			auto alnStart = p1_->alnRefSeq_.seq_.find_first_not_of('-');
			auto alnStop = p1_->alnRefSeq_.seq_.find_last_not_of('-') + 1;

			auto portion = p1_->alnQuerySeq_.getSubRead(alnStart, alnStop - alnStart);
			portion.removeGaps();
			auto start = p1_->querySeq_.seq_.find(portion.seq_);
			auto stop = start + portion.seq_.length();
			if(p1_->gRegion_.reverseSrand_){
				portion.reverseComplementRead(false, true);
				auto revStop = p1_->querySeq_.seq_.size() - start;
				auto revStart = p1_->querySeq_.seq_.size() - stop;
				start = revStart;
				stop = revStop;
			}
			MetaDataInName meta;
			meta.addMeta("start", start);
			meta.addMeta("stop", stop);
			portion.resetMetaInName(meta);
			return portion;
		}

		[[nodiscard]] seqInfo getP2PortionOriginal53() const {
			auto alnStart = p2_->alnRefSeq_.seq_.find_first_not_of('-');
			auto alnStop = p2_->alnRefSeq_.seq_.find_last_not_of('-') + 1;

			auto portion = p2_->alnQuerySeq_.getSubRead(alnStart, alnStop - alnStart);
			portion.removeGaps();
			auto start = p2_->querySeq_.seq_.find(portion.seq_);
			auto stop = start + portion.seq_.length();
			if(p2_->gRegion_.reverseSrand_){
				portion.reverseComplementRead(false, true);
				auto revStop = p2_->querySeq_.seq_.size() - start;
				auto revStart = p2_->querySeq_.seq_.size() - stop;
				start = revStart;
				stop = revStop;
			}
			MetaDataInName meta;
			meta.addMeta("start", start);
			meta.addMeta("stop", stop);
			portion.resetMetaInName(meta);
			return portion;
		}
	};
	std::vector<GenomeExtractResultWithHits> extractions;
	if(results.size() > 1){
		for(const auto pos1 : iter::range(results.size() - 1)){
			for(const auto & pos2 : iter::range(pos1, results.size())){
				//need to be on the same chromosome
				//need to be on opposite strands (should both should be in 5'->3' direction
				//and they shouldn't overlap

				std::shared_ptr<ReAlignedSeq> p1 = results[pos1];
				std::shared_ptr<ReAlignedSeq> p2 = results[pos2];

				if (getRef(p1).gRegion_.chrom_ == getRef(p2).gRegion_.chrom_
						&& getRef(p1).gRegion_.reverseSrand_ != getRef(p2).gRegion_.reverseSrand_
						&& !getRef(p1).gRegion_.overlaps(getRef(p2).gRegion_)
						&& getRef(p1).gRegion_.start_ != getRef(p2).gRegion_.end_
						&& getRef(p1).gRegion_.end_ != getRef(p2).gRegion_.start_ ) {

					if(getRef(p1).gRegion_.reverseSrand_){
						if(getRef(p1).gRegion_.start_ > getRef(p2).gRegion_.start_){
							GenomeExtractResultWithHits extraction(p1, p2);
							if (extraction.extraction_.gRegion_->getLen() <= maxTargetSize && extraction.extraction_.gRegion_->getLen() >=minTargetSize) {
//								std::cout << __FILE__ << " " << __LINE__ << std::endl;
//								p1->refSeq_.outPutSeqAnsi(std::cout);
//								p1->querySeq_.outPutSeqAnsi(std::cout);
//								p1->alnRefSeq_.outPutSeqAnsi(std::cout);
//								p1->alnQuerySeq_.outPutSeqAnsi(std::cout);
//								std::cout << extraction.p1_->gRegion_.genBedRecordCore().toDelimStrWithExtra() << std::endl;
//								std::cout << extraction.getP1PortionOriginal53().seq_.size() << std::endl;
//								extraction.getP1PortionOriginal53().outPutSeqAnsi(std::cout);

//								p2->refSeq_.outPutSeqAnsi(std::cout);
//								p2->querySeq_.outPutSeqAnsi(std::cout);
//								p2->alnRefSeq_.outPutSeqAnsi(std::cout);
//								p2->alnQuerySeq_.outPutSeqAnsi(std::cout);
//								std::cout << extraction.p2_->gRegion_.genBedRecordCore().toDelimStrWithExtra() << std::endl;
//								std::cout << extraction.getP2PortionOriginal53().seq_.size() << std::endl;
//								extraction.getP2PortionOriginal53().outPutSeqAnsi(std::cout);
								extractions.emplace_back(extraction);

							}
						}
					} else {
						if(getRef(p1).gRegion_.start_ < getRef(p2).gRegion_.start_){
							GenomeExtractResultWithHits extraction(p1, p2);
							if (extraction.extraction_.gRegion_->getLen() <= maxTargetSize && extraction.extraction_.gRegion_->getLen() >=minTargetSize) {
//								std::cout << __FILE__ << " " << __LINE__ << std::endl;
//								p1->refSeq_.outPutSeqAnsi(std::cout);
//								p1->querySeq_.outPutSeqAnsi(std::cout);
//								p1->alnRefSeq_.outPutSeqAnsi(std::cout);
//								p1->alnQuerySeq_.outPutSeqAnsi(std::cout);
//								std::cout << extraction.p1_->gRegion_.genBedRecordCore().toDelimStrWithExtra() << std::endl;
//								std::cout << extraction.getP1PortionOriginal53().seq_.size() << std::endl;
//								extraction.getP1PortionOriginal53().outPutSeqAnsi(std::cout);

//								p2->refSeq_.outPutSeqAnsi(std::cout);
//								p2->querySeq_.outPutSeqAnsi(std::cout);
//								p2->alnRefSeq_.outPutSeqAnsi(std::cout);
//								p2->alnQuerySeq_.outPutSeqAnsi(std::cout);
//								std::cout << extraction.p2_->gRegion_.genBedRecordCore().toDelimStrWithExtra() << std::endl;
//								std::cout << extraction.getP2PortionOriginal53().seq_.size() << std::endl;
//								extraction.getP2PortionOriginal53().outPutSeqAnsi(std::cout);
								extractions.emplace_back(extraction);
							}
						}
					}
				}
			}
		}
	}

	njh::sort(extractions, [](const GenomeExtractResultWithHits &reg1In, const GenomeExtractResultWithHits &reg2In) {
		const auto &reg1 = getRef(reg1In);
		const auto &reg2 = getRef(reg2In);
		if (reg1.extraction_.gRegion_->chrom_ == reg2.extraction_.gRegion_->chrom_) {
			if (reg1.extraction_.gRegion_->start_ == reg2.extraction_.gRegion_->start_) {
				return reg1.extraction_.gRegion_->end_ < reg2.extraction_.gRegion_->end_;
			} else {
				return reg1.extraction_.gRegion_->start_ < reg2.extraction_.gRegion_->start_;
			}
		} else {
			return reg1.extraction_.gRegion_->chrom_ < reg2.extraction_.gRegion_->chrom_;
		}
	});


	OutputStream extractionInfos(njh::files::make_path(setUp.pars_.directoryName_, "extractionInfos.tsv"));
	extractionInfos << "chrom\tstart\tend\tname\tlen\tstrand\tgenomicID\texpected\tp1Name\tp1TarName\tp1_5-3_orig\tp1_5-3_portion\tp1ChromStart\tp1Len\tp1LenFracCov\tp1Errors\tp2Name\tp2TarName\tp2_5-3_orig\tp2_5-3_portion\tp2ChromStart\tp2Len\tp2LenFracCov\tp2Errors" << std::endl;

	OutputStream allExtractionBed(njh::files::make_path(setUp.pars_.directoryName_, "allGenomicRegions.bed"));
	OutputStream expectedExtractionBed(njh::files::make_path(setUp.pars_.directoryName_, "expectedExtractionGenomicRegion.bed"));
	OutputStream unexpectedExtractionBed(njh::files::make_path(setUp.pars_.directoryName_, "unexpectedExtractionGenomicRegion.bed"));

	table outputPrimersForTesting(VecStr{"target", "forward", "reverse"});

	for(const auto & extraction : extractions){
		std::string extPrimerTarName;
		bool extForwardPrimer;
		std::string ligPrimerTarName;
		bool ligForwardPrimer;
		if (njh::in(extraction.extraction_.extRegion_.uid_, forPrimerNameToTargetName)) {
			extPrimerTarName = forPrimerNameToTargetName[extraction.extraction_.extRegion_.uid_];
			extForwardPrimer = true;
		} else {
			extPrimerTarName = revPrimerNameToTargetName[extraction.extraction_.extRegion_.uid_];
			extForwardPrimer = false;
		}
		if (njh::in(extraction.extraction_.ligRegion_.uid_, forPrimerNameToTargetName)) {
			ligPrimerTarName = forPrimerNameToTargetName[extraction.extraction_.ligRegion_.uid_];
			ligForwardPrimer = true;
		} else {
			ligPrimerTarName = revPrimerNameToTargetName[extraction.extraction_.ligRegion_.uid_];
			ligForwardPrimer = false;
		}
		//it's an expected amplification if the targets are the same and they are opposite primers
		bool expected = ligPrimerTarName == extPrimerTarName && ligForwardPrimer != extForwardPrimer &&
						primerByName[extraction.extraction_.extRegion_.uid_] == extraction.getP1PortionOriginal53().seq_ &&
						primerByName[extraction.extraction_.ligRegion_.uid_] == extraction.getP2PortionOriginal53().seq_ ;
		auto outRegion = extraction.extraction_.gRegion_->genBedRecordCore();
		outRegion.extraFields_.emplace_back(outRegion.name_);
		outRegion.name_ = outRegion.genUIDFromCoordsWithStrand();
		if (expected) {
			expectedExtractionBed << outRegion.toDelimStrWithExtra() << std::endl;
		} else {
			unexpectedExtractionBed << outRegion.toDelimStrWithExtra() << std::endl;
		}
		allExtractionBed << outRegion.toDelimStrWithExtra() << std::endl;
		outputPrimersForTesting.addRow(
						outRegion.name_,
						extraction.getP1PortionOriginal53().seq_,
						extraction.getP2PortionOriginal53().seq_
						);
		extractionInfos << extraction.extraction_.gRegion_->genBedRecordCore().toDelimStr()
										<< "\t" << extraction.extraction_.gRegion_->genBedRecordCore().genUIDFromCoordsWithStrand()
										<< "\t" << njh::boolToStr(expected)
										<< "\t" << extraction.extraction_.extRegion_.uid_
										<< "\t" << primerNameToTargetName[extraction.extraction_.extRegion_.uid_]
										<< "\t" << primerByName[extraction.extraction_.extRegion_.uid_]
										<< "\t" << extraction.getP1PortionOriginal53().seq_
										<< "\t" << extraction.extraction_.extRegion_.start_
										<< "\t" << extraction.extraction_.extRegion_.getLen()
										<< "\t" << extraction.extraction_.extRegion_.getLen() /
															 static_cast<double>(primerByName[extraction.extraction_.extRegion_.uid_].size())
										<< "\t" << extraction.extraction_.extComp_.distances_.getNumOfEvents(true)
										<< "\t" << extraction.extraction_.ligRegion_.uid_
										<< "\t" << primerNameToTargetName[extraction.extraction_.ligRegion_.uid_]
										<< "\t" << primerByName[extraction.extraction_.ligRegion_.uid_]
										<< "\t" << extraction.getP2PortionOriginal53().seq_
										<< "\t" << extraction.extraction_.ligRegion_.start_
										<< "\t" << extraction.extraction_.ligRegion_.getLen()
										<< "\t" << extraction.extraction_.ligRegion_.getLen() /
															 static_cast<double>(primerByName[extraction.extraction_.ligRegion_.uid_].size())
										<< "\t" << extraction.extraction_.ligComp_.distances_.getNumOfEvents(true)
										<< std::endl;
	}
	outputPrimersForTesting.outPutContents(TableIOOpts::genTabFileOut(njh::files::make_path(setUp.pars_.directoryName_, "primersForTest.tab.txt")));

	return 0;
}



int primerUtilsRunner::computeDimerizationScore(
				const njh::progutils::CmdArgs &inputCommands) {
	OutOptions outOpts("", ".tsv");
	bool outputRawMatrix = false;
	uint32_t numThreads = 1;
	bfs::path primersFnp;
	seqSetUp setUp(inputCommands);
	setUp.setOption(primersFnp, "--primersFnp", "primer table to read in instead of a fasta file");
	setUp.setOption(numThreads, "--numThreads", "number of CPUs to use");

	setUp.processReadInNames(VecStr{"--fasta", "--fastagz"}, primersFnp.empty());
	setUp.processWritingOptions(outOpts);
	setUp.setOption(outputRawMatrix, "--outputRawMatrix", "Output raw matrix without column and row names");
	setUp.finishSetUp(std::cout);
	OutputStream out(outOpts);

	std::vector<seqInfo> primers;
	if (primersFnp.empty()) {
		primers = SeqInput::getSeqVec<seqInfo>(setUp.pars_.ioOptions_);
	} else {
		PrimersAndMids ids(primersFnp);

		if(ids.getTargets().empty()){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error in reading in target primers file " << primersFnp << "\n";
			ss << "Make sure there is a line that starts with target in file" << "\n";
			throw std::runtime_error{ss.str()};
		}
		ids.initPrimerDeterminator();
		//currently only using raw primers, no multi primers
		for(const auto & primerInfo : ids.pDeterminator_->primers_){
			primers.emplace_back(njh::pasteAsStr(primerInfo.first, "_F"), primerInfo.second.forwardPrimerRaw_);
			primers.emplace_back(njh::pasteAsStr(primerInfo.first, "_R"), primerInfo.second.reversePrimerRaw_);
		}
	}
	PrimerDimerUtils dimerScorer;
	std::vector<std::vector<double>> scores = dimerScorer.computeFullScoreMatrix(primers, numThreads);
	PrimerDimerUtils::writeMatrix(scores, out, primers, outputRawMatrix);


//	dimerScorer.debug_  = true;
//	{
//		auto bestScore = dimerScorer.computeDimerScoreTop("TTTTATTTATTTGTTTAACTACAACATTTT", "AATATATTGCACCCATGTCAAAAATGTTGTAGAT");
//		std::cout << bestScore << std::endl;
//	}//ATCTACAACATTTTTTTTATTTATTTGTTT AATATATTGCACCCATGTCAAAAATGTTGTAGAT
//	{// TTTATTTATTTGTTTAATATATTGCACCCATGT
//		auto bestScore = dimerScorer.computeDimerScoreTop( "AATATATTGCACCCATGTCAAAAATGTTGTAGAT", "TTTTATTTATTTGTTTAACTACAACATTTT");
//		std::cout << bestScore << std::endl;
//	}
//	exit(1);
//	{
//		//C867 score = 0
//		auto bestScoreFirst = dimerScorer.computeDimerScoreTop("GTCATGCGCCCCATAACACCCTAAAATCCCTATCA", "CTCGGACGCACCCATTATAATTTTTGGAGATTTTAAAGATG");
//		std::cout << "bestScoreFirst: " << bestScoreFirst << std::endl;
//	}
//	{
//		//C829 score = -2.79
//		auto bestScoreFirst = dimerScorer.computeDimerScoreTop("GTCATGCGCCCCATATATGGGTATAGGGTTTTTAAAAG", "CTCGGACGCACCCACTACCCTCCTCCCTCTAA");
//		std::cout << "bestScoreFirst: " << bestScoreFirst << std::endl;
//	}
//	{
//		//C847 score = -2.83
//		auto bestScoreFirst = dimerScorer.computeDimerScoreTop("GTCATGCGCCCCATATAATATTTTTAACAAAAAACCACATTT", "CTCGGACGCACCCATAAGTTTTAATTGTTGTTTTTGTTTTAG");
//		std::cout << "bestScoreFirst: " << bestScoreFirst << std::endl;
//	}
//	{
//		//C943 score = -7.08
//		auto bestScoreFirst = dimerScorer.computeDimerScoreTop("GTCATGCGCCCCATTTTAAGTTTTTTTTTTTTTGTTTTTTGAAT", "CTCGGACGCACCCATTATATTCAAATAAATAATACATAAAATCA");
//		std::cout << "bestScoreFirst: " << bestScoreFirst << std::endl;
//	}
//	exit(1);
//	auto bestScoreFirst = dimerScorer.computeDimerScore("TGTGCATGCATGTGTGTGTG", "GCTCGTCGTTGATCCACAGA");
//	auto bestScoreSecond = dimerScorer.computeDimerScore("TCTGTGGATCAACGACGAGC", "CACACACACATGCATGCACA");
//	std::cout << "bestScoreFirst: " << bestScoreFirst << std::endl;
//	std::cout << "bestScoreSecond: " << bestScoreSecond << std::endl;
//	exit(1);
	return 0;
}


}  // namespace njhseq

