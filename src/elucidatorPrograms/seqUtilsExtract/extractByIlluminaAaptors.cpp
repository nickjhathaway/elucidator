/*
 * extractByIlluminaAaptors.cpp
 *
 *  Created on: Nov 28, 2019
 *      Author: nicholashathaway
 */




#include "seqUtilsExtractRunner.hpp"

namespace njhseq {


struct PrimerAlnInfo{
	PrimerAlnInfo(const std::string & primerName, const comparison & comp,
			const seqInfo & seqAln, const seqInfo & primerAln) :
			primerName_(primerName),
			comp_(comp),
			seqAln_(seqAln),
			primerAln_(primerAln) {
		auto firstN = primerAln.seq_.find("N");
		auto lastN = primerAln.seq_.rfind("N") + 1;
		auto alnBarcode = seqAln.seq_.substr(firstN, lastN - firstN);

		barcode_ = njh::replaceString(alnBarcode, "-", "");
	}
	std::string primerName_;
	comparison comp_;
	seqInfo seqAln_;
	seqInfo primerAln_;

	std::string barcode_;
};


/*
 * 		changeSubStrToLower(seq.seq_, ,  - );
		changeSubStrToLower(seq.seq_, fBarPosEnd,  -fBarPosEnd);
		changeSubStrToLower(seq.seq_, , rBarPosStart - );
		changeSubStrToLower(seq.seq_, rBarPosEnd, rPrimerEnd -rBarPosEnd);
 */
struct ProcessedForIlluminaAdaptorRes{
	bool fwdPass_{false};
	bool revPass_{false};

	uint32_t fprimerStart_ = 0;
	uint32_t fprimerEnd_ = 0;
	uint32_t fBarPosStart_ = 0;
	uint32_t fBarPosEnd_ = 0;

	uint32_t rprimerStart_ = 0;
	uint32_t rprimerEnd_ = 0;
	uint32_t rBarPosStart_ = 0;
	uint32_t rBarPosEnd_ = 0;





};

struct ProcessedForIlluminaAdaptorPars{

	bool keepPrimer_ = false;
	bool primerUpper_ = false;

};

ProcessedForIlluminaAdaptorRes processSeqForAdaptors(seqInfo & seq,
		const PrimerAlnInfo & frontSeqAdaptorInfo,
		const PrimerAlnInfo & backSeqAdaptorInfo,
		const uint32_t startOfBackSeq,
		const ProcessedForIlluminaAdaptorPars & pars){
	ProcessedForIlluminaAdaptorRes ret;
	//forward primer start

	if('-' == frontSeqAdaptorInfo.primerAln_.seq_.front()){
		ret.fprimerStart_ = getRealPosForAlnPos(frontSeqAdaptorInfo.seqAln_.seq_, frontSeqAdaptorInfo.primerAln_.seq_.find_first_not_of("-"));
	}
	//forward primer barcode start
	uint32_t firstSeqBase = frontSeqAdaptorInfo.seqAln_.seq_.find_first_not_of("-");
	uint32_t firstNAlnPos = frontSeqAdaptorInfo.primerAln_.seq_.find("N");
	if(firstNAlnPos > firstSeqBase){
		ret.fBarPosStart_ = getRealPosForAlnPos(frontSeqAdaptorInfo.seqAln_.seq_, firstNAlnPos);
		//forward primer barcode end
		ret.fBarPosEnd_= getRealPosForAlnPos(frontSeqAdaptorInfo.seqAln_.seq_, frontSeqAdaptorInfo.primerAln_.seq_.rfind("N")) + 1;
		//check the forward primer
		if('-' != frontSeqAdaptorInfo.seqAln_.seq_.back() && ret.fBarPosStart_ > 4 ){
			ret.fprimerEnd_= getRealPosForAlnPos(frontSeqAdaptorInfo.seqAln_.seq_, frontSeqAdaptorInfo.primerAln_.seq_.find_last_not_of("-")) + 1;
			ret.fwdPass_ = true;
		}
	}

	//the back primer should be found within the back seq
	if(backSeqAdaptorInfo.seqAln_.seq_.front() != '-'){
		uint32_t firstN_AlnPos = backSeqAdaptorInfo.primerAln_.seq_.find("N");
		uint32_t lastN_AlnPos = backSeqAdaptorInfo.primerAln_.seq_.rfind("N");
		uint32_t backSeqAlnEnd = backSeqAdaptorInfo.seqAln_.seq_.find_last_not_of("-");
		if(lastN_AlnPos < backSeqAlnEnd){
			ret.rprimerEnd_= seq.seq_.size();
			if(backSeqAdaptorInfo.seqAln_.seq_.back() != '-'){
				ret.rprimerEnd_= getRealPosForAlnPos(backSeqAdaptorInfo.seqAln_.seq_, backSeqAdaptorInfo.primerAln_.seq_.find_last_not_of("-")) + 1;
			}
			ret.rprimerStart_ = getRealPosForAlnPos(backSeqAdaptorInfo.seqAln_.seq_, backSeqAdaptorInfo.primerAln_.seq_.find_first_not_of("-"));
			ret.rBarPosStart_ = getRealPosForAlnPos(backSeqAdaptorInfo.seqAln_.seq_, firstN_AlnPos);
			ret.rBarPosEnd_ = getRealPosForAlnPos(backSeqAdaptorInfo.seqAln_.seq_, lastN_AlnPos) + 1;
			if( ((seq.seq_.size() - startOfBackSeq) - ret.rBarPosEnd_) > 4){
				ret.revPass_= true;
			}
			ret.rprimerEnd_+=startOfBackSeq;
			ret.rprimerStart_+=startOfBackSeq;
			ret.rBarPosStart_+=startOfBackSeq;
			ret.rBarPosEnd_+=startOfBackSeq;
		}
	}

	if(ret.fwdPass_ && ret.revPass_){
		if(pars.keepPrimer_){
			if(!pars.primerUpper_){
				changeSubStrToLower(seq.seq_, ret.fprimerStart_, ret.fprimerEnd_ - ret.fprimerStart_);
				changeSubStrToLower(seq.seq_, ret.rprimerStart_, ret.rprimerEnd_ - ret.rprimerStart_);
			}
		}
	}

	return ret;
}

int seqUtilsExtractRunner::extractByIlluminaAaptors(const njh::progutils::CmdArgs & inputCommands) {
	std::string illuminaForward = "AATGATACGGCGACCACCGAGATCTACACNNNNNNNNACACTCTTTCCCTACACGACTCGCCAAGCTGA";
	std::string illuminaReverse = "CAAGCAGAAGACGGCATACGAGATNNNNNNNNGTGACTGGAGTTCAGACGTGTGCTCTTCCGAT";
	double primerIdentityCutOff = .70;
	double barIdentityCutOff = .60;
	ProcessedForIlluminaAdaptorPars extractPars;
	bfs::path illuminaBarcodeSampleSheet = "";

	seqUtilsExtractSetUp setUp(inputCommands);
	setUp.setOption(illuminaForward, "--illuminaForward", "Illumina Forward adaptor/primer");
	setUp.setOption(illuminaReverse, "--illuminaReverse", "Illumina Reverse adaptor/primer");

	setUp.setOption(extractPars.primerUpper_, "--primerUpper", "When keeping primer, keep upper case");
	setUp.setOption(extractPars.keepPrimer_, "--keepPrimer", "Keep Primer sequence, otherwise remove it");

	setUp.setOption(barIdentityCutOff, "--barIdentityCutOff", "Bar Identity Cut Off");
	setUp.setOption(primerIdentityCutOff, "--primerIdentityCutOff", "Percent Identity Cut Off for matching the primers");
	setUp.setOption(illuminaBarcodeSampleSheet, "--illuminaBarcodeSampleSheet", "Illumina Barcode Sample Sheet", true);

	setUp.processVerbose();
	setUp.processDebug();
	setUp.processReadInNames(true);
	setUp.processDirectoryOutputName(true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);


	//read in barcodes
	//check for same barcodes
	table barcodesTab(illuminaBarcodeSampleSheet, "\t", true);
	barcodesTab.checkForColumnsThrow(VecStr{"sample_name", "fw", "rev"}, __PRETTY_FUNCTION__);
	struct IlluminaDualBarCode{
		IlluminaDualBarCode(const std::string & samp, const std::string & fwBar5_3,
				const std::string & revBar5_3) :
				sampleName_(samp), fwBar_("fw", fwBar5_3), revBar_("rev", revBar5_3) {

		}
		std::string sampleName_;
		seqInfo fwBar_;
		seqInfo revBar_;
	};

	std::vector<IlluminaDualBarCode> barcodes;

	uint64_t maxBarLen = 0;

	for(const auto & row : barcodesTab){
		auto bar = IlluminaDualBarCode{
			row[barcodesTab.getColPos("sample_name")],
			row[barcodesTab.getColPos("fw")],
			row[barcodesTab.getColPos("rev")]};
		//check other barcodes for same pairing and for same sample name
		for(const auto & prevBar : barcodes){
			if(bar.sampleName_ == prevBar.sampleName_){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "already have sample name: " << bar.sampleName_<< "\n";
				throw std::runtime_error{ss.str()};
			}
			if (bar.fwBar_.seq_ == prevBar.fwBar_.seq_
					&& bar.revBar_.seq_ == prevBar.revBar_.seq_) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error "
						<< "already have barcode pairing: fw: " << bar.fwBar_.seq_
						<< " rev: " << bar.revBar_.seq_ << " for sample: "
						<< prevBar.sampleName_ << "\n";
				throw std::runtime_error { ss.str() };
			}
		}
		readVec::getMaxLength(bar.fwBar_, maxBarLen);
		readVec::getMaxLength(bar.revBar_, maxBarLen);

		barcodes.emplace_back(bar);
	}
	//make outputs
	auto fastqDir = njh::files::make_path(setUp.pars_.directoryName_, "fastq");
	njh::files::makeDir(njh::files::MkdirPar{fastqDir});
	MultiSeqOutCache<seqInfo> outputs;
	for(const auto & bar : barcodes){
		outputs.addReader(bar.sampleName_, SeqIOOptions::genFastqOutGz(njh::files::make_path(fastqDir, bar.sampleName_)));
	}
	outputs.addReader("Undetermined", SeqIOOptions::genFastqOutGz(njh::files::make_path(fastqDir, "Undetermined")));
	//set up for alignment of barcodes for scoring
	maxBarLen = maxBarLen * 2;
	gapScoringParameters gapParsBar(5,1,5,1,5,1);
	auto scoringBar = substituteMatrix::createDegenScoreMatrix(2,-2);
	aligner alignerBarObj(maxBarLen, gapParsBar, scoringBar, true);

	gapScoringParameters gapPars(5,1,0,0,0,0);
	gapPars.gapLeftRefOpen_ = 5;
	gapPars.gapLeftRefExtend_ = 1;
	gapPars.gapRightRefOpen_ = 5;
	gapPars.gapRightRefExtend_ = 1;

	auto scoring = substituteMatrix::createDegenScoreMatrix(2,-2);
	uint64_t maxLen = 500;
	aligner alignerObj(maxLen, gapPars, scoring, false);

	seqInfo fwd5_3("fwd", illuminaForward);
	seqInfo fwd3_5("fwd", seqUtil::reverseComplement(illuminaForward, "DNA"));

	seqInfo rev5_3("rev", illuminaReverse);
	seqInfo rev3_5("rev", seqUtil::reverseComplement(illuminaReverse, "DNA"));


	SeqInput seqReader(setUp.pars_.ioOptions_);
	seqReader.openIn();
	seqInfo seq;


	OutputStream out(njh::files::make_path(setUp.pars_.directoryName_, "info.tab.txt"));
	out << "seqName\tfwd5_3_score\tfwd5_3_bar\trev3_5_score\trev3_5_bar\trev5_3_score\trev5_3_bar\tfwd3_5_score\tfwd3_5_bar" << std::endl;
	OutputStream outBest(njh::files::make_path(setUp.pars_.directoryName_, "infoBest.tab.txt"));
	outBest << "seqName\tdirection\tfbarcode\trbarcode" << std::endl;

	OutputStream outPass(njh::files::make_path(setUp.pars_.directoryName_, "infoPass.tab.txt"));

	outPass << "seqName\tdirection\tfwd5_3_score\tfwd5_3_bar\trev3_5_score\trev3_5_bar\trev5_3_score\trev5_3_bar\tfwd3_5_score\tfwd3_5_bar" << std::endl;

	SeqOutput writer(SeqIOOptions::genFastqOut(njh::files::make_path(setUp.pars_.directoryName_, "processed.fastq")));
	writer.openOut();

	std::map<std::string, std::map<std::string, uint32_t>> bestCounts;
	std::map<std::string, uint32_t> sampleCounts;

	uint32_t smallSeq = 0;


	uint32_t failed_fwd5_3_Identity = 0;
	uint32_t failed_fwd3_5_Identity = 0;
	uint32_t failed_rev5_3_Identity = 0;
	uint32_t failed_rev3_5_Identity = 0;

	uint32_t failed_fwd5_3_PrimerAln = 0;
	uint32_t failed_fwd3_5_PrimerAln = 0;
	uint32_t failed_rev5_3_PrimerAln = 0;
	uint32_t failed_rev3_5_PrimerAln = 0;

	uint32_t mismatchDirections = 0;
	double total = 0;
	uint32_t passFwd = 0;
	uint32_t passRev = 0;
	while(seqReader.readNextRead(seq)){
		++total;
		if(seq.seq_.size() < 251){
			++smallSeq;
			continue;
		}
		seqInfo frontSeq(seq.getSubRead(0, 125));
		seqInfo backSeq(seq.getSubRead(len(seq) -125));


		alignerObj.alignCacheGlobal(frontSeq, fwd5_3);
		alignerObj.profilePrimerAlignment(frontSeq, fwd5_3);
		PrimerAlnInfo fwd5_3_info("fwd5_3", alignerObj.comp_, alignerObj.alignObjectA_.seqBase_, alignerObj.alignObjectB_.seqBase_);

		alignerObj.alignCacheGlobal(backSeq, rev3_5);
		alignerObj.profilePrimerAlignment(backSeq, rev3_5);
		PrimerAlnInfo rev3_5_info("rev3_5", alignerObj.comp_, alignerObj.alignObjectA_.seqBase_, alignerObj.alignObjectB_.seqBase_);

		alignerObj.alignCacheGlobal(frontSeq, rev5_3);
		alignerObj.profilePrimerAlignment(frontSeq, rev5_3);
		PrimerAlnInfo rev5_3_info("rev5_3", alignerObj.comp_, alignerObj.alignObjectA_.seqBase_, alignerObj.alignObjectB_.seqBase_);

		alignerObj.alignCacheGlobal(backSeq, fwd3_5);
		alignerObj.profilePrimerAlignment(backSeq, fwd3_5);
		PrimerAlnInfo fwd3_5_info("fwd3_5", alignerObj.comp_, alignerObj.alignObjectA_.seqBase_, alignerObj.alignObjectB_.seqBase_);

		out << seq.name_
			<< "\t" << fwd5_3_info.comp_.distances_.eventBasedIdentity_
			<< "\t" << fwd5_3_info.barcode_
			<< "\t" << rev3_5_info.comp_.distances_.eventBasedIdentity_
			<< "\t" << rev3_5_info.barcode_
			<< "\t" << rev5_3_info.comp_.distances_.eventBasedIdentity_
			<< "\t" << rev5_3_info.barcode_
			<< "\t" << fwd3_5_info.comp_.distances_.eventBasedIdentity_
			<< "\t" << fwd3_5_info.barcode_
			<< std::endl;

		std::string bestFront = rev5_3_info.comp_.distances_.eventBasedIdentity_ < fwd5_3_info.comp_.distances_.eventBasedIdentity_ ? "rev" : "fwd";
		std::string bestBack = fwd3_5_info.comp_.distances_.eventBasedIdentity_ < rev3_5_info.comp_.distances_.eventBasedIdentity_ ? "rev" : "fwd";
		++bestCounts[bestFront][bestBack];
		auto forwardScore = fwd5_3_info.comp_.distances_.eventBasedIdentity_ + rev3_5_info.comp_.distances_.eventBasedIdentity_;
		auto reverseScore = rev5_3_info.comp_.distances_.eventBasedIdentity_ + fwd3_5_info.comp_.distances_.eventBasedIdentity_;
		std::string direction = "fwd";
		std::string forwardBar = fwd5_3_info.barcode_;
		std::string reverseBar = seqUtil::reverseComplement(rev3_5_info.barcode_, "DNA");
		// check if best score for both forward and reverse
		if(bestFront != bestBack){
			++mismatchDirections;
			continue;
		}
		ProcessedForIlluminaAdaptorRes processRes;
		if(reverseScore > forwardScore){
			direction = "rev";
			forwardBar = seqUtil::reverseComplement(fwd3_5_info.barcode_, "DNA");;
			reverseBar = rev5_3_info.barcode_;
			bool failedIdentity = false;
			if(rev5_3_info.comp_.distances_.eventBasedIdentity_ < primerIdentityCutOff){
				//
				++failed_rev5_3_Identity;
				failedIdentity = true;
			}
			if(fwd3_5_info.comp_.distances_.eventBasedIdentity_ < primerIdentityCutOff){
				//
				++failed_fwd3_5_Identity;
				failedIdentity = true;
			}
			if(failedIdentity){
				continue;
			}
			uint32_t startOfBackSeq = seq.seq_.size() - 125;
			processRes = processSeqForAdaptors(seq, rev5_3_info, fwd3_5_info, startOfBackSeq, extractPars);
			if(!processRes.fwdPass_){
				++failed_rev5_3_PrimerAln;
			}
			if(!processRes.revPass_){
				++failed_fwd3_5_PrimerAln;
			}
			if(processRes.fwdPass_ && processRes.revPass_){
				if(extractPars.keepPrimer_){
					seq = seq.getSubRead(processRes.fprimerStart_, processRes.rprimerEnd_- processRes.fprimerStart_);
				}else{
					seq = seq.getSubRead(processRes.fprimerEnd_, processRes.rprimerStart_- processRes.fprimerEnd_);
				}
				seq.reverseComplementRead(false, true);
				++passRev;
			}
		} else {
			//process seq
			bool failedIdentity = false;
			if(fwd5_3_info.comp_.distances_.eventBasedIdentity_ < primerIdentityCutOff){
				//
				++failed_fwd5_3_Identity;
				failedIdentity = true;
			}
			if(rev3_5_info.comp_.distances_.eventBasedIdentity_ < primerIdentityCutOff){
				//
				++failed_rev3_5_Identity;
				failedIdentity = true;
			}
			if(failedIdentity){
				continue;
			}
			uint32_t startOfBackSeq = seq.seq_.size() - 125;
			processRes = processSeqForAdaptors(seq, fwd5_3_info, rev3_5_info, startOfBackSeq, extractPars);
			if(!processRes.fwdPass_){
				++failed_fwd5_3_PrimerAln;
			}
			if(!processRes.revPass_){
				++failed_rev3_5_PrimerAln;
			}
			if(processRes.fwdPass_ && processRes.revPass_){
				if(extractPars.keepPrimer_){
					seq = seq.getSubRead(processRes.fprimerStart_, processRes.rprimerEnd_- processRes.fprimerStart_);
				}else{
					seq = seq.getSubRead(processRes.fprimerEnd_, processRes.rprimerStart_- processRes.fprimerEnd_);
				}
				++passFwd;
			}
		}

		if(processRes.fwdPass_ && processRes.revPass_){
			std::string determinedSample = "Undetermined";

			{
				VecStr fwdBestSample{"Undetermined"};
				VecStr revBestSample{"Undetermined"};
				uint32_t bestFwdIds = 0;
				uint32_t bestFwdAlnScore = 0;
				uint32_t bestRevIds = 0;
				uint32_t bestRevAlnScore = 0;
				seqInfo fwdBarSeq("determined-fwd", forwardBar);
				seqInfo revBarSeq("determined-rev", reverseBar);

				if(forwardBar.size() < maxBarLen &&
						reverseBar.size() < maxBarLen &&
						forwardBar.size() > 5 &&
						reverseBar.size() > 5){
					for(const auto & bar : barcodes){
						alignerBarObj.alignCacheGlobal(bar.fwBar_, fwdBarSeq);
						alignerBarObj.profilePrimerAlignment(bar.fwBar_, fwdBarSeq);
						if(alignerBarObj.comp_.distances_.eventBasedIdentity_ >=barIdentityCutOff){
							if(alignerBarObj.comp_.distances_.ref_.identities_ > bestFwdIds){
								fwdBestSample.clear();
								bestFwdAlnScore = alignerBarObj.parts_.score_;
								bestFwdIds = alignerBarObj.comp_.distances_.ref_.identities_;
								fwdBestSample.emplace_back(bar.sampleName_);
							}else if(alignerBarObj.comp_.distances_.ref_.identities_ == bestFwdIds){
								if(alignerBarObj.parts_.score_ > bestFwdAlnScore){
									fwdBestSample.clear();
									bestFwdAlnScore = alignerBarObj.parts_.score_;
									bestFwdIds = alignerBarObj.comp_.distances_.ref_.identities_;
									fwdBestSample.emplace_back(bar.sampleName_);
								}else if(alignerBarObj.parts_.score_ == bestFwdAlnScore){
									fwdBestSample.emplace_back(bar.sampleName_);
								}
							}
						}
						alignerBarObj.alignCacheGlobal(bar.revBar_, revBarSeq);
						alignerBarObj.profilePrimerAlignment(bar.revBar_, revBarSeq);
						if(alignerBarObj.comp_.distances_.eventBasedIdentity_ >=barIdentityCutOff){
							if(alignerBarObj.comp_.distances_.ref_.identities_ > bestRevIds){
								revBestSample.clear();
								bestRevAlnScore = alignerBarObj.parts_.score_;
								bestRevIds = alignerBarObj.comp_.distances_.ref_.identities_;
								revBestSample.emplace_back(bar.sampleName_);
							}else if(alignerBarObj.comp_.distances_.ref_.identities_ == bestRevIds){
								if(alignerBarObj.parts_.score_ > bestRevAlnScore){
									revBestSample.clear();
									bestRevAlnScore = alignerBarObj.parts_.score_;
									bestRevIds = alignerBarObj.comp_.distances_.ref_.identities_;
									revBestSample.emplace_back(bar.sampleName_);
								}else if(alignerBarObj.parts_.score_ == bestRevAlnScore){
									revBestSample.emplace_back(bar.sampleName_);
								}
							}
						}
					}
				}
				//need one for each and need to match
				//if multiple in one or both, see if there is only one pair since barcodes can be reused several times
				if (   njh::in(std::string("Undetermined"),fwdBestSample)
						|| njh::in(std::string("Undetermined"),revBestSample)
						|| fwdBestSample.size() == 0
						|| revBestSample.size() == 0) {
					determinedSample = "Undetermined";
				}else{
					if(fwdBestSample.size() == 1 && revBestSample.size() == 1 && fwdBestSample.front() == revBestSample.front()){
						determinedSample = fwdBestSample.front();
					} else {
						VecStr pairs;
						for(const auto & fwdBest : fwdBestSample){
							for(const auto & revBest : revBestSample){
								if(revBest  == fwdBest){
									pairs.emplace_back(fwdBest);
									break;
								}
							}
						}
						if(pairs.size() == 1){
							determinedSample = pairs.front();
						}else{
							determinedSample = "Undetermined";
						}
					}
				}
			}

			MetaDataInName seqMeta;
			seqMeta.addMeta("direction", direction);
			seqMeta.addMeta("fprimerStart", processRes.fprimerStart_);
			seqMeta.addMeta("forwardBar", forwardBar);
			seqMeta.addMeta("reverseBar", reverseBar);
			seqMeta.addMeta("fullBar", njh::pasteAsStr(forwardBar, "+", reverseBar));
			seqMeta.addMeta("rprimerStart", processRes.rprimerStart_);
			seqMeta.addMeta("sample", determinedSample);
			seq.name_ += seqMeta.createMetaName();
			outputs.add(determinedSample, seq);
			++sampleCounts[determinedSample];
			outPass << seq.name_
				<< "\t" << direction
				<< "\t" << fwd5_3_info.comp_.distances_.eventBasedIdentity_
				<< "\t" << fwd5_3_info.barcode_
				<< "\t" << rev3_5_info.comp_.distances_.eventBasedIdentity_
				<< "\t" << rev3_5_info.barcode_
				<< "\t" << rev5_3_info.comp_.distances_.eventBasedIdentity_
				<< "\t" << rev5_3_info.barcode_
				<< "\t" << fwd3_5_info.comp_.distances_.eventBasedIdentity_
				<< "\t" << fwd3_5_info.barcode_
				<< std::endl;
			outBest << seq.name_
					<< "\t" << direction
					<< "\t" << forwardBar
					<< "\t" << reverseBar << std::endl;
		}
	}
	OutputStream outBestCounts(njh::files::make_path(setUp.pars_.directoryName_, "bestDirectionCounts.tab.txt"));
	outBestCounts << "BestFrontPrimer\tBestBackPrimer\tcount\t" << std::endl;
	for(const auto & front : bestCounts){
		for(const auto & back : front.second){
			outBestCounts << front.first << "\t" << back.first << "\t" << back.second << std::endl;
		}
	}

	OutputStream outFailedCounts(njh::files::make_path(setUp.pars_.directoryName_, "failedCounts.tab.txt"));
	outFailedCounts << "Failed\tcount\t" << std::endl;
	{
		outFailedCounts << "failed_fwd5_3_Identity" << "\t" << failed_fwd5_3_Identity << std::endl;
		outFailedCounts << "failed_fwd3_5_Identity" << "\t" << failed_fwd3_5_Identity << std::endl;
		outFailedCounts << "failed_rev5_3_Identity" << "\t" << failed_rev5_3_Identity << std::endl;
		outFailedCounts << "failed_rev3_5_Identity" << "\t" << failed_rev3_5_Identity << std::endl;

		outFailedCounts << "failed_fwd5_3_PrimerAln" << "\t" << failed_fwd5_3_PrimerAln << std::endl;
		outFailedCounts << "failed_fwd3_5_PrimerAln" << "\t" << failed_fwd3_5_PrimerAln << std::endl;
		outFailedCounts << "failed_rev5_3_PrimerAln" << "\t" << failed_rev5_3_PrimerAln << std::endl;
		outFailedCounts << "failed_rev3_5_PrimerAln" << "\t" << failed_rev3_5_PrimerAln << std::endl;

	}

	{
		uint32_t totalFailedIdentity = failed_fwd5_3_Identity + failed_fwd3_5_Identity + failed_rev5_3_Identity + failed_rev3_5_Identity;
		uint32_t totalFailedPrimerAln = failed_fwd5_3_PrimerAln + failed_fwd3_5_PrimerAln + failed_rev5_3_PrimerAln + failed_rev3_5_PrimerAln;

		OutputStream conditionCount(njh::files::make_path(setUp.pars_.directoryName_, "conditionCount.tab.txt"));
		conditionCount << "condition\tcount\tfrac" << std::endl;
		conditionCount << "passFwd" << "\t" << passFwd << "\t" << passFwd/total << std::endl;
		conditionCount << "passRev" << "\t" << passRev << "\t" << passRev/total << std::endl;
		conditionCount << "failedIdentity" << "\t" << totalFailedIdentity << "\t" << totalFailedIdentity/total << std::endl;
		conditionCount << "totalFailedPrimerAln" << "\t" << totalFailedPrimerAln << "\t" << totalFailedPrimerAln/total << std::endl;
		conditionCount << "smallSeq" << "\t" << smallSeq << "\t" << smallSeq/total << std::endl;
		conditionCount << "mismatchDirections" << "\t" << mismatchDirections << "\t" << mismatchDirections/total << std::endl;

	}
	{
		OutputStream sampleCountsOut(njh::files::make_path(setUp.pars_.directoryName_, "sampleCounts.tab.txt"));
		sampleCountsOut << "sample\tcount\tfrac" << std::endl;
		for(const auto & sampleCount : sampleCounts){
			sampleCountsOut << sampleCount.first
					<< "\t" << sampleCount.second
					<< "\t" << sampleCount.second/static_cast<double>(passFwd + passRev) << std::endl;
		}

	}



	return 0;
}


}  // namespace njhseq


