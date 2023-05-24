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
#include <njhseq/objects/BioDataObject/BioRecordsUtils/BedUtility.hpp>

namespace njhseq {



primerUtilsRunner::primerUtilsRunner()
				: njh::progutils::ProgramRunner(
				{
								addFunc("computeDimerizationScore", computeDimerizationScore, false),
								addFunc("testWithBlastForUnspecificAmplification", testWithBlastForUnspecificAmplification, false),
								//
				},
				"primerUtils") {}

int primerUtilsRunner::testWithBlastForUnspecificAmplification(
				const njh::progutils::CmdArgs &inputCommands) {

	uint32_t minLen = 12;
	uint32_t errorAllowed = 0;
	uint32_t sizeLimit = 6000;
	bfs::path primersFnp;
	bfs::path genomeFnp;
	seqSetUp setUp(inputCommands);
	setUp.setOption(minLen, "--minLen", "min length of the 3` end of the alignment");
	setUp.setOption(errorAllowed, "--errorAllowed", "errors allowed");
	setUp.setOption(sizeLimit, "--sizeLimit", "extraction size limits");

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

	auto allPrimerFnp = njh::files::make_path(setUp.pars_.directoryName_, "allPrimers.fasta");
	SeqOutput::write(primerSeqs, SeqIOOptions::genFastaOut(allPrimerFnp));
	for(const auto & seq : primerSeqs){
		primerLengths[seq.name_] = len(seq);
		primerByName[seq.name_] = seq.seq_;
	}
	auto allPrimerBlastHitsArchiveFnp = njh::files::make_path(setUp.pars_.directoryName_, "allPrimersBlast.archive");
	auto allPrimerBlastHitsTableFnp = njh::files::make_path(setUp.pars_.directoryName_, "allPrimersBlast.tsv");

	//_blastdb
	//.nsq
	auto blastDatabaseFnp = bfs::path(njh::rstripRet(bfs::path(genomeFnp).replace_extension("").string(), '.') + "_blastdb");
	std::string blastCmd = "blastn -query " + allPrimerFnp.string()  + " -task blastn-short -word_size 6 -db " + blastDatabaseFnp.string() + " -outfmt 6 -max_target_seqs 100000000 > " + allPrimerBlastHitsTableFnp.string();
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
			//check if the 3` end, if in plus strand, check the end, if in the reverseStrand first base
			if(!hit.reverseStrand()){
				if(hit.qEnd_ != njh::mapAt(primerLengths, hit.queryName_)){
					continue;
				}
			} else {
				if(hit.qStart_ > 1){
					continue;
				}
			}

			blastHits.emplace_back(hit);
			auto realignedSeq = ReAlignedSeq::genRealignment(hit, njh::mapAt(primerByName, hit.queryName_), alignerObj,
																											 chromLens, genomeReader, reAlnPars);
//			std::cout << "passed: " << njh::colorBool(allowableErrors.passErrorProfile(realignedSeq.comp_)) << std::endl;
			if(allowableErrors.passErrorProfile(realignedSeq.comp_)){
//				std::cout << __FILE__ << " " << __LINE__ << std::endl;
//				realignedSeq.alnRefSeq_.outPutSeqAnsi(std::cout);
//				realignedSeq.alnQuerySeq_.outPutSeqAnsi(std::cout);
//				std::cout << __FILE__ << " " << __LINE__ << std::endl;
				//std::cout << hit.toJson() << std::endl;
				passingBlastHitsInfos << hit.genSubjectBed6().toDelimStrWithExtra() << std::endl;
				results.emplace_back(std::make_shared<ReAlignedSeq>(realignedSeq));
			}
		}
	}
	std::vector<GenomeExtractResult> extractions;
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
							GenomeExtractResult extraction(p1, p2);
							extraction.setRegion();
							if (extraction.gRegion_->getLen() <= sizeLimit) {
								extractions.emplace_back(extraction);
							}
						}
					}else{
						if(getRef(p1).gRegion_.start_ < getRef(p2).gRegion_.start_){
							GenomeExtractResult extraction(p1, p2);
							extraction.setRegion();
							if (extraction.gRegion_->getLen() <= sizeLimit) {
								extractions.emplace_back(extraction);
							}
						}
					}
				}
			}
		}
	}

	njh::sort(extractions, [](const GenomeExtractResult &reg1In, const GenomeExtractResult &reg2In) {
		const auto &reg1 = getRef(reg1In);
		const auto &reg2 = getRef(reg2In);
		if (reg1.gRegion_->chrom_ == reg2.gRegion_->chrom_) {
			if (reg1.gRegion_->start_ == reg2.gRegion_->start_) {
				return reg1.gRegion_->end_ < reg2.gRegion_->end_;
			} else {
				return reg1.gRegion_->start_ < reg2.gRegion_->start_;
			}
		} else {
			return reg1.gRegion_->chrom_ < reg2.gRegion_->chrom_;
		}
	});


	OutputStream extractionInfos(njh::files::make_path(setUp.pars_.directoryName_, "extractionInfos.tsv"));
	extractionInfos << "chrom\tstart\tend\tname\tlen\tstrand\texpected\tp1Name\tp1TarName\tp1_5-3\tp1ChromStart\tp1Len\tp1LenFracCov\tp1Errors\tp2Name\tp2TarName\tp2_5-3\tp2ChromStart\tp2Len\tp2LenFracCov\tp2Errors" << std::endl;

	OutputStream expectedExtractionBed(njh::files::make_path(setUp.pars_.directoryName_, "expectedExtractionGenomicRegion.bed"));
	OutputStream unexpectedExtractionBed(njh::files::make_path(setUp.pars_.directoryName_, "unexpectedExtractionGenomicRegion.bed"));
	for(const auto & extraction : extractions){
		std::string extPrimerTarName;
		bool extForwardPrimer;
		std::string ligPrimerTarName;
		bool ligForwardPrimer;
		if (njh::in(extraction.extRegion_.uid_, forPrimerNameToTargetName)) {
			extPrimerTarName = forPrimerNameToTargetName[extraction.extRegion_.uid_];
			extForwardPrimer = true;
		} else {
			extPrimerTarName = revPrimerNameToTargetName[extraction.extRegion_.uid_];
			extForwardPrimer = false;
		}
		if (njh::in(extraction.ligRegion_.uid_, forPrimerNameToTargetName)) {
			ligPrimerTarName = forPrimerNameToTargetName[extraction.ligRegion_.uid_];
			ligForwardPrimer = true;
		} else {
			ligPrimerTarName = revPrimerNameToTargetName[extraction.ligRegion_.uid_];
			ligForwardPrimer = false;
		}
		//it's an expected amplification if the targets are the same and they are opposite primers
		bool expected = ligPrimerTarName == extPrimerTarName && ligForwardPrimer != extForwardPrimer;
		if(expected){
			expectedExtractionBed << extraction.gRegion_->genBedRecordCore().toDelimStrWithExtra() << std::endl;
		}else{
			unexpectedExtractionBed << extraction.gRegion_->genBedRecordCore().toDelimStrWithExtra() << std::endl;
		}
		extractionInfos << extraction.gRegion_->genBedRecordCore().toDelimStr()
										<< "\t" << njh::boolToStr(expected)
										<< "\t" << extraction.extRegion_.uid_
										<< "\t" << primerNameToTargetName[extraction.extRegion_.uid_]
										<< "\t" << primerByName[extraction.extRegion_.uid_]
										<< "\t" << extraction.extRegion_.start_
										<< "\t" << extraction.extRegion_.getLen()
										<< "\t" << extraction.extRegion_.getLen() /
															 static_cast<double>(primerByName[extraction.extRegion_.uid_].size())
										<< "\t" << extraction.extComp_.distances_.getNumOfEvents(true)
										<< "\t" << extraction.ligRegion_.uid_
										<< "\t" << primerNameToTargetName[extraction.ligRegion_.uid_]
										<< "\t" << primerByName[extraction.ligRegion_.uid_]
										<< "\t" << extraction.ligRegion_.start_
										<< "\t" << extraction.ligRegion_.getLen()
										<< "\t" << extraction.ligRegion_.getLen() /
															 static_cast<double>(primerByName[extraction.ligRegion_.uid_].size())
										<< "\t" << extraction.ligComp_.distances_.getNumOfEvents(true)
										<< std::endl;
	}

	return 0;
}



int primerUtilsRunner::computeDimerizationScore(
				const njh::progutils::CmdArgs &inputCommands) {
	OutOptions outOpts("", ".tsv");
	bool outputRawMatrix = false;
	bfs::path primersFnp;
	seqSetUp setUp(inputCommands);
	setUp.setOption(primersFnp, "--primersFnp", "primer table to read in instead of a fasta file");

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
	std::vector<std::vector<double>> scores = dimerScorer.computeFullScoreMatrix(primers);
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

