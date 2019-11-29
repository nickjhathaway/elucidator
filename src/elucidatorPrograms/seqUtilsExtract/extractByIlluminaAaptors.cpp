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
			if(ret.fprimerStart_ > 0){
				ret.fwdPass_ = true;
			}
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
			ret.rprimerStart_ = getRealPosForAlnPos(backSeqAdaptorInfo.seqAln_.seq_, backSeqAdaptorInfo.primerAln_.seq_.find_first_not_of("-")) + 1;
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
				changeSubStrToLower(seq.seq_, ret.fprimerStart_, ret.fprimerEnd_ - ret.fprimerEnd_);
				changeSubStrToLower(seq.seq_, ret.rprimerStart_, ret.rprimerEnd_ - ret.rprimerStart_);
			}
		}
	}

	return ret;
}

int seqUtilsExtractRunner::extractByIlluminaAaptors(const njh::progutils::CmdArgs & inputCommands) {
	std::string illuminaForward = "AATGATACGGCGACCACCGAGATCTACACNNNNNNNNACACTCTTTCCCTACACGACTCGCCAAGCTGA";
	std::string illuminaReverse = "CAAGCAGAAGACGGCATACGAGATNNNNNNNNGTGACTGGAGTTCAGACGTGTGCTCTTCCGAT";
	double identityCutOff = .70;
	ProcessedForIlluminaAdaptorPars extractPars;
	seqUtilsExtractSetUp setUp(inputCommands);
	setUp.setOption(extractPars.primerUpper_, "--primerUpper", "When keeping primer, keep upper case");
	setUp.setOption(extractPars.keepPrimer_, "--keepPrimer", "Keep Primer sequence, otherwise remove it");

	setUp.setOption(identityCutOff, "--identityCutOff", "Percent Identity Cut Off for matching the primers");
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processReadInNames(true);
	setUp.processDirectoryOutputName(true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
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
		// add score cut off, possibly 0.60;

		if(reverseScore > forwardScore){
			direction = "rev";
			forwardBar = seqUtil::reverseComplement(fwd3_5_info.barcode_, "DNA");;
			reverseBar = rev5_3_info.barcode_;

			bool failedIdentity = false;
			if(rev5_3_info.comp_.distances_.eventBasedIdentity_ < identityCutOff){
				//
				++failed_rev5_3_Identity;
				failedIdentity = true;
			}
			if(fwd3_5_info.comp_.distances_.eventBasedIdentity_ < identityCutOff){
				//
				++failed_fwd3_5_Identity;
				failedIdentity = true;
			}
			if(failedIdentity){
				continue;
			}
			uint32_t startOfBackSeq = seq.seq_.size() - 125;
			auto processRes = processSeqForAdaptors(seq, rev5_3_info, fwd3_5_info, startOfBackSeq, extractPars);
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
				MetaDataInName seqMeta;
				seqMeta.addMeta("direction", direction);
				seqMeta.addMeta("fprimerStart", processRes.fprimerStart_);
				seqMeta.addMeta("forwardBar", forwardBar);
				seqMeta.addMeta("reverseBar", reverseBar);
				seqMeta.addMeta("fullBar", njh::pasteAsStr(forwardBar, "+", reverseBar));
				seqMeta.addMeta("rprimerStart", processRes.rprimerStart_);
				seq.name_ += seqMeta.createMetaName();
				writer.write(seq);
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
				++passRev;
			}
		} else {
			//process seq

			bool failedIdentity = false;
			if(fwd5_3_info.comp_.distances_.eventBasedIdentity_ < identityCutOff){
				//
				++failed_fwd5_3_Identity;
				failedIdentity = true;
			}
			if(rev3_5_info.comp_.distances_.eventBasedIdentity_ < identityCutOff){
				//
				++failed_rev3_5_Identity;
				failedIdentity = true;
			}
			if(failedIdentity){
				continue;
			}
			uint32_t startOfBackSeq = seq.seq_.size() - 125;
			auto processRes = processSeqForAdaptors(seq, fwd5_3_info, rev3_5_info, startOfBackSeq, extractPars);
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
				MetaDataInName seqMeta;
				seqMeta.addMeta("direction", direction);
				seqMeta.addMeta("fprimerStart", processRes.fprimerStart_);
				seqMeta.addMeta("forwardBar", forwardBar);
				seqMeta.addMeta("reverseBar", reverseBar);
				seqMeta.addMeta("rprimerStart", processRes.rprimerStart_);
				seqMeta.addMeta("fullBar", njh::pasteAsStr(forwardBar, "+", reverseBar));
				seq.name_ += seqMeta.createMetaName();
				writer.write(seq);
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
				++passFwd;
			}
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



	return 0;
}


}  // namespace njhseq


