/*
 * extractByIlluminaAaptors.cpp
 *
 *  Created on: Nov 28, 2019
 *      Author: nicholashathaway
 */




#include "seqUtilsExtractRunner.hpp"

#include <njhseq/IO/SeqIO.h>
#include <njhseq/alignment/aligner/aligner.hpp>
#include <optional>
#include <utility>


namespace njhseq {


struct PrimerAlnInfo{
	PrimerAlnInfo(std::string  primerName, comparison  comp,
			const seqInfo & seqAln, const seqInfo & primerAln) :
			primerName_(std::move(primerName)),
			comp_(std::move(comp)),
			seqAln_(seqAln),
			primerAln_(primerAln) {
		auto firstN = primerAln.seq_.find('N');
		auto lastN = primerAln.seq_.rfind('N') + 1;
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
	//0 requires the barcode in the back of the read to be whole and followed by more than 4 bases, otherwise the barcode may run off the end of the read as long as this many of its bases are left (nanopore reads lose the end of the strand)
	uint32_t minTruncatedBackBarcodeBases_ = 0;

};

//the exclusive end in the un-gapped sequence for an alignment column, i.e. the number of sequence bases at or before that column
//getRealPosForAlnPos(...) + 1 overshoots by one when the column is a gap in the sequence
uint32_t getRealEndForAlnPos(const std::string & seqAln, size_t alnPos){
	return std::count_if(seqAln.begin(), seqAln.begin() + alnPos + 1, [](char base){ return '-' != base; });
}

void checkAdaptorHasBarcode(const std::string & adaptor, const std::string & optionName){
	if(std::string::npos == adaptor.find('N')){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << optionName << " needs to contain the barcode as N's, adaptor: " << adaptor << "\n";
		throw std::runtime_error{ss.str()};
	}
}

ProcessedForIlluminaAdaptorRes processSeqForAdaptors(seqInfo & seq,
		const PrimerAlnInfo & frontSeqAdaptorInfo,
		const PrimerAlnInfo & backSeqAdaptorInfo,
		const uint32_t startOfBackSeq,
		const ProcessedForIlluminaAdaptorPars & pars){
	ProcessedForIlluminaAdaptorRes ret;
	//forward primer start
	if('-' == frontSeqAdaptorInfo.primerAln_.seq_.front()){
		ret.fprimerStart_ = getRealPosForAlnPos(frontSeqAdaptorInfo.seqAln_.seq_, frontSeqAdaptorInfo.primerAln_.seq_.find_first_not_of('-'));
	}
	//forward primer barcode start
	uint32_t firstSeqBase = frontSeqAdaptorInfo.seqAln_.seq_.find_first_not_of('-');
	uint32_t firstNAlnPos = frontSeqAdaptorInfo.primerAln_.seq_.find('N');
	if(firstNAlnPos > firstSeqBase){
		ret.fBarPosStart_ = getRealPosForAlnPos(frontSeqAdaptorInfo.seqAln_.seq_, firstNAlnPos);
		//forward primer barcode end
		ret.fBarPosEnd_= getRealEndForAlnPos(frontSeqAdaptorInfo.seqAln_.seq_, frontSeqAdaptorInfo.primerAln_.seq_.rfind('N'));
		//check the forward primer
		if('-' != frontSeqAdaptorInfo.seqAln_.seq_.back() && ret.fBarPosStart_ > 4 ){
			ret.fprimerEnd_= getRealEndForAlnPos(frontSeqAdaptorInfo.seqAln_.seq_, frontSeqAdaptorInfo.primerAln_.seq_.find_last_not_of('-'));
			ret.fwdPass_ = true;
		}
	}

	//the back primer should be found within the back seq
	if(backSeqAdaptorInfo.seqAln_.seq_.front() != '-'){
		uint32_t firstN_AlnPos = backSeqAdaptorInfo.primerAln_.seq_.find('N');
		uint32_t lastN_AlnPos = backSeqAdaptorInfo.primerAln_.seq_.rfind('N');
		uint32_t backSeqAlnEnd = backSeqAdaptorInfo.seqAln_.seq_.find_last_not_of('-');
		//a barcode cut short by the end of the read still counts when enough of it is left
		bool truncatedBarcodeOk = pars.minTruncatedBackBarcodeBases_ > 0
				&& firstN_AlnPos < backSeqAlnEnd
				&& backSeqAdaptorInfo.barcode_.size() >= pars.minTruncatedBackBarcodeBases_;
		if(lastN_AlnPos < backSeqAlnEnd || truncatedBarcodeOk){
			//positions here are relative to the back seq, shifted by startOfBackSeq below; if the primer runs off the end of the read this is the back seq length
			ret.rprimerEnd_= getRealEndForAlnPos(backSeqAdaptorInfo.seqAln_.seq_, backSeqAdaptorInfo.primerAln_.seq_.find_last_not_of('-'));
			ret.rprimerStart_ = getRealPosForAlnPos(backSeqAdaptorInfo.seqAln_.seq_, backSeqAdaptorInfo.primerAln_.seq_.find_first_not_of('-'));
			ret.rBarPosStart_ = getRealPosForAlnPos(backSeqAdaptorInfo.seqAln_.seq_, firstN_AlnPos);
			ret.rBarPosEnd_ = getRealEndForAlnPos(backSeqAdaptorInfo.seqAln_.seq_, lastN_AlnPos);
			if( (seq.seq_.size() - startOfBackSeq) - ret.rBarPosEnd_ > 4 || truncatedBarcodeOk){
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

//settings shared by the programs that look for a barcoded adaptor at each end of a read
struct AdaptorExtractPars{
	std::string forwardAdaptor_ = "AATGATACGGCGACCACCGAGATCTACACNNNNNNNNACACTCTTTCCCTACACGACTCGCCAAGCTGA";
	std::string reverseAdaptor_ = "CAAGCAGAAGACGGCATACGAGATNNNNNNNNGTGACTGGAGTTCAGACGTGTGCTCTTCCGAT";
	double primerIdentityCutOff_ = .70;
	double barIdentityCutOff_ = .60;
	uint32_t searchWithin_ = 125;
	uint32_t minOutputLen_ = 251;
	uint32_t maxOutputLen_ = std::numeric_limits<uint32_t>::max();
	ProcessedForIlluminaAdaptorPars extractPars_;
};

void setAdaptorExtractOptions(seqUtilsExtractSetUp & setUp, AdaptorExtractPars & pars){
	setUp.setOption(pars.searchWithin_, "--searchWithin", "search Within");
	setUp.setOption(pars.minOutputLen_, "--minOutputLen", "minimum Output sequence Length, otherwise throw these out");
	setUp.setOption(pars.maxOutputLen_, "--maxOutputLen", "maximum Output sequence Length, otherwise throw these out");

	setUp.setOption(pars.forwardAdaptor_, "--illuminaForward,--forwardAdaptor", "Forward adaptor/primer, with the barcode as N's");
	setUp.setOption(pars.reverseAdaptor_, "--illuminaReverse,--reverseAdaptor", "Reverse adaptor/primer, with the barcode as N's");

	setUp.setOption(pars.extractPars_.primerUpper_, "--primerUpper", "When keeping primer, keep upper case");
	setUp.setOption(pars.extractPars_.keepPrimer_, "--keepPrimer", "Keep Primer sequence, otherwise remove it");
	setUp.setOption(pars.extractPars_.minTruncatedBackBarcodeBases_, "--minTruncatedBackBarcodeBases",
			"Let the barcode at the end of the read be cut off by the end of the read as long as this many of its bases are left, 0 requires the whole barcode");

	setUp.setOption(pars.barIdentityCutOff_, "--barIdentityCutOff", "Bar Identity Cut Off");
	setUp.setOption(pars.primerIdentityCutOff_, "--primerIdentityCutOff", "Percent Identity Cut Off for matching the primers");
}

//outcome of looking for the forward and reverse adaptors at the two ends of a read
enum class AdaptorReadOutcome{
	smallSeq,
	largeSeq,
	mismatchDirections,
	failedIdentity,
	failedPrimerAln,
	failedOverlapCheck,
	pass
};

struct AdaptorReadRes{
	AdaptorReadOutcome outcome_ = AdaptorReadOutcome::smallSeq;
	std::string direction_ = "fwd";
	std::string forwardBar_;
	std::string reverseBar_;
	ProcessedForIlluminaAdaptorRes processRes_;
	std::optional<PrimerAlnInfo> fwd5_3_;
	std::optional<PrimerAlnInfo> rev3_5_;
	std::optional<PrimerAlnInfo> rev5_3_;
	std::optional<PrimerAlnInfo> fwd3_5_;
};

/**@brief Find the forward and reverse adaptors at the two ends of reads, pull out the barcodes in them and trim the reads down to what is between them,
 * keeping the counts countIlluminaAaptors and extractByIlluminaAaptors/extractByNanoporeBarcodes report
 */
class AdaptorPairFinder{
public:
	AdaptorPairFinder(const AdaptorExtractPars & pars, bool debug) :
			pars_(pars), debug_(debug),
			alignerObj_(std::max<uint64_t>({500, pars.searchWithin_, pars.forwardAdaptor_.size(), pars.reverseAdaptor_.size()}),
					genGapPars(), substituteMatrix::createDegenScoreMatrixLessN(2,-2), false),
			fwd5_3_("fwd", pars.forwardAdaptor_),
			fwd3_5_("fwd", seqUtil::reverseComplement(pars.forwardAdaptor_, "DNA")),
			rev5_3_("rev", pars.reverseAdaptor_),
			rev3_5_("rev", seqUtil::reverseComplement(pars.reverseAdaptor_, "DNA")) {
		checkAdaptorHasBarcode(pars_.forwardAdaptor_, "--illuminaForward");
		checkAdaptorHasBarcode(pars_.reverseAdaptor_, "--illuminaReverse");
	}

	/**@brief look for the adaptors in seq, on a pass seq is trimmed to between the adaptors (reverse complemented when in the reverse direction)
	 * @param seq the read
	 * @param infoOut where to log the barcode alignment info for every read long enough to be checked
	 */
	AdaptorReadRes process(seqInfo & seq, std::ostream & infoOut){
		AdaptorReadRes ret;
		++total_;
		if(seq.seq_.size() < pars_.minOutputLen_){
			++smallSeq_;
			ret.outcome_ = AdaptorReadOutcome::smallSeq;
			return ret;
		}
		if(seq.seq_.size() > pars_.maxOutputLen_){
			++largeSeq_;
			ret.outcome_ = AdaptorReadOutcome::largeSeq;
			return ret;
		}
		//get a subsequence to align to from the back and front
		seqInfo frontSeq(seq.getSubRead(0, std::min<uint32_t>(pars_.searchWithin_, len(seq))));
		seqInfo backSeq(seq.getSubRead(len(seq) -std::min<uint32_t>(pars_.searchWithin_, len(seq))));

		//create alignment of all primers/adapters
		ret.fwd5_3_ = alignAdaptor("fwd5_3", frontSeq, fwd5_3_);
		ret.rev3_5_ = alignAdaptor("rev3_5", backSeq, rev3_5_);
		ret.rev5_3_ = alignAdaptor("rev5_3", frontSeq, rev5_3_);
		ret.fwd3_5_ = alignAdaptor("fwd3_5", backSeq, fwd3_5_);
		if(debug_) {
			std::cout << std::endl;
		}
		const auto & fwd5_3_info = *ret.fwd5_3_;
		const auto & rev3_5_info = *ret.rev3_5_;
		const auto & rev5_3_info = *ret.rev5_3_;
		const auto & fwd3_5_info = *ret.fwd3_5_;

		//log the info
		infoOut << seq.name_
			<< "\t" << len(seq)
			<< "\t" << fwd5_3_info.comp_.distances_.eventBasedIdentity_
			<< "\t" << fwd5_3_info.barcode_
			<< "\t" << rev3_5_info.comp_.distances_.eventBasedIdentity_
			<< "\t" << rev3_5_info.barcode_
			<< "\t" << rev5_3_info.comp_.distances_.eventBasedIdentity_
			<< "\t" << rev5_3_info.barcode_
			<< "\t" << fwd3_5_info.comp_.distances_.eventBasedIdentity_
			<< "\t" << fwd3_5_info.barcode_
			<< "\n";

		//determine the best alignment for the front of the sequence of either the forward primer or the reverse primer
		std::string bestFront = fwd5_3_info.comp_.distances_.eventBasedIdentity_ < rev5_3_info.comp_.distances_.eventBasedIdentity_ ? "rev" : "fwd";
		std::string bestBack = rev3_5_info.comp_.distances_.eventBasedIdentity_ < fwd3_5_info.comp_.distances_.eventBasedIdentity_ ? "rev" : "fwd";
		++bestCounts_[bestFront][bestBack];

		//get the scores to see which direction should be considered
		auto forwardScore = fwd5_3_info.comp_.distances_.eventBasedIdentity_ + rev3_5_info.comp_.distances_.eventBasedIdentity_;
		auto reverseScore = rev5_3_info.comp_.distances_.eventBasedIdentity_ + fwd3_5_info.comp_.distances_.eventBasedIdentity_;

		// check if best score for both forward and reverse
		if(bestFront != bestBack){
			++mismatchDirections_;
			ret.outcome_ = AdaptorReadOutcome::mismatchDirections;
			return ret;
		}
		//assume the forward position unless the reverse score is better
		ret.forwardBar_ = fwd5_3_info.barcode_;
		ret.reverseBar_ = seqUtil::reverseComplement(rev3_5_info.barcode_, "DNA");
		bool reverse = reverseScore > forwardScore;
		if(reverse){
			ret.direction_ = "rev";
			ret.forwardBar_ = seqUtil::reverseComplement(fwd3_5_info.barcode_, "DNA");
			ret.reverseBar_ = rev5_3_info.barcode_;
		}
		const auto & frontInfo = reverse ? rev5_3_info : fwd5_3_info;
		const auto & backInfo = reverse ? fwd3_5_info : rev3_5_info;
		auto & failedFrontIdentity = reverse ? failed_rev5_3_Identity_ : failed_fwd5_3_Identity_;
		auto & failedBackIdentity = reverse ? failed_fwd3_5_Identity_ : failed_rev3_5_Identity_;
		auto & failedFrontPrimerAln = reverse ? failed_rev5_3_PrimerAln_ : failed_fwd5_3_PrimerAln_;
		auto & failedBackPrimerAln = reverse ? failed_fwd3_5_PrimerAln_ : failed_rev3_5_PrimerAln_;

		bool failedIdentity = false;
		if(frontInfo.comp_.distances_.eventBasedIdentity_ < pars_.primerIdentityCutOff_){
			++failedFrontIdentity;
			failedIdentity = true;
		}
		if(backInfo.comp_.distances_.eventBasedIdentity_ < pars_.primerIdentityCutOff_){
			++failedBackIdentity;
			failedIdentity = true;
		}
		if(failedIdentity){
			++failedIdentityReads_;
			ret.outcome_ = AdaptorReadOutcome::failedIdentity;
			return ret;
		}
		uint32_t startOfBackSeq = seq.seq_.size() - std::min<uint32_t>(pars_.searchWithin_, len(seq));
		ret.processRes_ = processSeqForAdaptors(seq, frontInfo, backInfo, startOfBackSeq, pars_.extractPars_);
		const auto & processRes = ret.processRes_;
		if(!processRes.fwdPass_){
			++failedFrontPrimerAln;
		}
		if(!processRes.revPass_){
			++failedBackPrimerAln;
		}
		if(!processRes.fwdPass_ || !processRes.revPass_){
			++failedPrimerAlnReads_;
			ret.outcome_ = AdaptorReadOutcome::failedPrimerAln;
			return ret;
		}
		if(processRes.fprimerEnd_ >= processRes.rprimerStart_) {
			++failedOverlapCheck_;
			ret.outcome_ = AdaptorReadOutcome::failedOverlapCheck;
			return ret;
		}
		if (pars_.extractPars_.keepPrimer_) {
			seq = seq.getSubRead(processRes.fprimerStart_, processRes.rprimerEnd_ - processRes.fprimerStart_);
		} else {
			seq = seq.getSubRead(processRes.fprimerEnd_, processRes.rprimerStart_ - processRes.fprimerEnd_);
		}
		if(reverse){
			seq.reverseComplementRead(false, true);
			++passRev_;
		}else{
			++passFwd_;
		}
		ret.outcome_ = AdaptorReadOutcome::pass;
		return ret;
	}

	static void writePassLine(std::ostream & out, const std::string & seqName, const AdaptorReadRes & res){
		out << seqName
			<< "\t" << res.direction_
			<< "\t" << res.fwd5_3_->comp_.distances_.eventBasedIdentity_
			<< "\t" << res.fwd5_3_->barcode_
			<< "\t" << res.rev3_5_->comp_.distances_.eventBasedIdentity_
			<< "\t" << res.rev3_5_->barcode_
			<< "\t" << res.rev5_3_->comp_.distances_.eventBasedIdentity_
			<< "\t" << res.rev5_3_->barcode_
			<< "\t" << res.fwd3_5_->comp_.distances_.eventBasedIdentity_
			<< "\t" << res.fwd3_5_->barcode_
			<< "\n";
	}

	void writeBestDirectionCounts(const bfs::path & dir) const {
		OutputStream outBestCounts(njh::files::make_path(dir, "bestDirectionCounts.tab.txt"));
		outBestCounts << "BestFrontPrimer\tBestBackPrimer\tcount\t" << std::endl;
		for(const auto & front : bestCounts_){
			for(const auto & back : front.second){
				outBestCounts << front.first << "\t" << back.first << "\t" << back.second << std::endl;
			}
		}
	}

	void writeFailedCounts(const bfs::path & dir) const {
		OutputStream outFailedCounts(njh::files::make_path(dir, "failedCounts.tab.txt"));
		outFailedCounts << "Failed\tcount\t" << std::endl;
		outFailedCounts << "failed_fwd5_3_Identity" << "\t" << failed_fwd5_3_Identity_ << std::endl;
		outFailedCounts << "failed_fwd3_5_Identity" << "\t" << failed_fwd3_5_Identity_ << std::endl;
		outFailedCounts << "failed_rev5_3_Identity" << "\t" << failed_rev5_3_Identity_ << std::endl;
		outFailedCounts << "failed_rev3_5_Identity" << "\t" << failed_rev3_5_Identity_ << std::endl;

		outFailedCounts << "failed_fwd5_3_PrimerAln" << "\t" << failed_fwd5_3_PrimerAln_ << std::endl;
		outFailedCounts << "failed_fwd3_5_PrimerAln" << "\t" << failed_fwd3_5_PrimerAln_ << std::endl;
		outFailedCounts << "failed_rev5_3_PrimerAln" << "\t" << failed_rev5_3_PrimerAln_ << std::endl;
		outFailedCounts << "failed_rev3_5_PrimerAln" << "\t" << failed_rev3_5_PrimerAln_ << std::endl;
	}

	/**@brief write conditionCount.tab.txt
	 * @param dir output directory
	 * @param afterPasses extra conditions to report right after passFwd and passRev
	 */
	void writeConditionCounts(const bfs::path & dir, const std::vector<std::pair<std::string, uint32_t>> & afterPasses) const {
		OutputStream conditionCount(njh::files::make_path(dir, "conditionCount.tab.txt"));
		conditionCount << "condition\tcount\tfrac" << std::endl;
		conditionCount << "passFwd" << "\t" << passFwd_ << "\t" << passFwd_/total_ << std::endl;
		conditionCount << "passRev" << "\t" << passRev_ << "\t" << passRev_/total_ << std::endl;
		for(const auto & condition : afterPasses){
			conditionCount << condition.first << "\t" << condition.second << "\t" << condition.second/total_ << std::endl;
		}
		conditionCount << "failedIdentity" << "\t" << failedIdentityReads_ << "\t" << failedIdentityReads_/total_ << std::endl;
		conditionCount << "totalFailedPrimerAln" << "\t" << failedPrimerAlnReads_ << "\t" << failedPrimerAlnReads_/total_ << std::endl;
		conditionCount << "smallSeq" << "(len <" << pars_.minOutputLen_ << ")" << "\t" << smallSeq_ << "\t" << smallSeq_/total_ << std::endl;
		conditionCount << "largeSeq" << "(len >" << pars_.maxOutputLen_ << ")" << "\t" << largeSeq_ << "\t" << largeSeq_/total_ << std::endl;
		conditionCount << "mismatchDirections" << "\t" << mismatchDirections_ << "\t" << mismatchDirections_/total_ << std::endl;
		conditionCount << "failedOverlapCheck" << "\t" << failedOverlapCheck_ << "\t" << failedOverlapCheck_/total_ << std::endl;
	}

	uint32_t passes() const {
		return passFwd_ + passRev_;
	}

private:
	static gapScoringParameters genGapPars(){
		gapScoringParameters gapPars(5,1,0,0,0,0);
		gapPars.gapLeftRefOpen_ = 5;
		gapPars.gapLeftRefExtend_ = 1;
		gapPars.gapRightRefOpen_ = 5;
		gapPars.gapRightRefExtend_ = 1;
		return gapPars;
	}

	PrimerAlnInfo alignAdaptor(const std::string & name, const seqInfo & window, const seqInfo & adaptor){
		alignerObj_.alignRegGlobal(window, adaptor);
		alignerObj_.profilePrimerAlignment(window, adaptor);
		if(debug_) {
			alignerObj_.alignObjectA_.seqBase_.outPutSeqAnsi(std::cout);
			alignerObj_.alignObjectB_.seqBase_.outPutSeqAnsi(std::cout);
		}
		return PrimerAlnInfo(name, alignerObj_.comp_, alignerObj_.alignObjectA_.seqBase_, alignerObj_.alignObjectB_.seqBase_);
	}

	AdaptorExtractPars pars_;
	bool debug_;
	//the aligner does not grow its matrix, so it has to fit the search window and the adaptors
	aligner alignerObj_;
	seqInfo fwd5_3_;
	seqInfo fwd3_5_;
	seqInfo rev5_3_;
	seqInfo rev3_5_;

	std::map<std::string, std::map<std::string, uint32_t>> bestCounts_;

	double total_ = 0;
	uint32_t smallSeq_ = 0;
	uint32_t largeSeq_ = 0;

	uint32_t failed_fwd5_3_Identity_ = 0;
	uint32_t failed_fwd3_5_Identity_ = 0;
	uint32_t failed_rev5_3_Identity_ = 0;
	uint32_t failed_rev3_5_Identity_ = 0;

	uint32_t failed_fwd5_3_PrimerAln_ = 0;
	uint32_t failed_fwd3_5_PrimerAln_ = 0;
	uint32_t failed_rev5_3_PrimerAln_ = 0;
	uint32_t failed_rev3_5_PrimerAln_ = 0;

	uint32_t mismatchDirections_ = 0;
	//per read counts, the per primer counts above can count a read twice
	uint32_t failedIdentityReads_ = 0;
	uint32_t failedPrimerAlnReads_ = 0;
	uint32_t failedOverlapCheck_ = 0;

	uint32_t passFwd_ = 0;
	uint32_t passRev_ = 0;
};


int seqUtilsExtractRunner::countIlluminaAaptors(const njh::progutils::CmdArgs & inputCommands) {
	AdaptorExtractPars pars;

	seqUtilsExtractSetUp setUp(inputCommands);
	setAdaptorExtractOptions(setUp, pars);

	setUp.processVerbose();
	setUp.processDebug();
	setUp.processReadInNames(true);
	setUp.processDirectoryOutputName(true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);

	AdaptorPairFinder finder(pars, setUp.pars_.debug_);

	SeqInput seqReader(setUp.pars_.ioOptions_);
	seqReader.openIn();
	seqInfo seq;

	OutputStream out(njh::files::make_path(setUp.pars_.directoryName_, "info.tab.txt"));
	out << "seqName\tseqLen\tfwd5_3_score\tfwd5_3_bar\trev3_5_score\trev3_5_bar\trev5_3_score\trev5_3_bar\tfwd3_5_score\tfwd3_5_bar" << std::endl;
	OutputStream outBest(njh::files::make_path(setUp.pars_.directoryName_, "infoBest.tab.txt"));
	outBest << "seqName\tdirection\tfbarcode\trbarcode" << std::endl;

	OutputStream outPass(njh::files::make_path(setUp.pars_.directoryName_, "infoPass.tab.txt"));
	outPass << "seqName\tdirection\tfwd5_3_score\tfwd5_3_bar\trev3_5_score\trev3_5_bar\trev5_3_score\trev5_3_bar\tfwd3_5_score\tfwd3_5_bar" << std::endl;

	std::map<std::string, std::map<std::string, uint32_t>> bestBarcodePairsCounts;

	while(seqReader.readNextRead(seq)){
		auto res = finder.process(seq, out);
		if(AdaptorReadOutcome::pass != res.outcome_){
			continue;
		}
		AdaptorPairFinder::writePassLine(outPass, seq.name_, res);
		++bestBarcodePairsCounts[res.forwardBar_][res.reverseBar_];
		outBest << seq.name_
						<< "\t" << res.direction_
						<< "\t" << res.forwardBar_
						<< "\t" << res.reverseBar_ << "\n";
	}
	finder.writeBestDirectionCounts(setUp.pars_.directoryName_);
	{
		OutputStream outBestBarcodePairsCounts(njh::files::make_path(setUp.pars_.directoryName_, "bestBarcodePairsCounts.tab.txt"));
		outBestBarcodePairsCounts << "BestFrontPrimer\tBestBackPrimer\tcount\t" << std::endl;
		for(const auto & front : bestBarcodePairsCounts){
			for(const auto & back : front.second){
				outBestBarcodePairsCounts << front.first << "\t" << back.first << "\t" << back.second << std::endl;
			}
		}
	}
	finder.writeFailedCounts(setUp.pars_.directoryName_);
	finder.writeConditionCounts(setUp.pars_.directoryName_, {});
	return 0;
}


//barcode sequences of the ONT native (NB) and PCR (BC) barcoding kits, as listed in dorado's barcode_kits.cpp (v2.1.2); NB01-12 are the reverse complements of BC01-12, the rest are the same
const std::unordered_map<std::string, VecStr> & getOntBarcodeSeqs(){
	static const std::unordered_map<std::string, VecStr> seqs{
			{"NB", {
				"CACAAAGACACCGACAACTTTCTT", "ACAGACGACTACAAACGGAATCGA", "CCTGGTAACTGGGACACAAGACTC", "TAGGGAAACACGATAGAATCCGAA",
				"AAGGTTACACAAACCCTGGACAAG", "GACTACTTTCTGCCTTTGCGAGAA", "AAGGATTCATTCCCACGGTAACAC", "ACGTAACTTGGTTTGTTCCCTGAA",
				"AACCAAGACTCGCTGTGCCTAGTT", "GAGAGGACAAAGGTTTCAACGCTT", "TCCATTCCCTCCGATAGATGAAAC", "TCCGATTCTGCTTCTTTCTACCTG",
				"AGAACGACTTCCATACTCGTGTGA", "AACGAGTCTCTTGGGACCCATAGA", "AGGTCTACCTCGCTAACACCACTG", "CGTCAACTGACAGTGGTTCGTACT",
				"ACCCTCCAGGAAAGTACCTCTGAT", "CCAAACCCAACAACCTAGATAGGC", "GTTCCTCGTGCAGTGTCAAGAGAT", "TTGCGTCCTGTTACGAGAACTCAT",
				"GAGCCTCTCATTGTCCGTTCTCTA", "ACCACTGCCATGTATCAAAGTACG", "CTTACTACCCAGTGAACCTCCTCG", "GCATAGTTCTGCATGATGGGTTAG",
				"GTAAGTTGGGTATGCAACGCAATG", "CATACAGCGACTACGCATTCTCAT", "CGACGGTTAGATTCACCTCTTACA", "TGAAACCTAAGAAGGCACCGTATC",
				"CTAGACACCTTGGGTTGACAGACC", "TCAGTGAGGATCTACTTCGACCCA", "TGCGTACAGCAATCAGTTACATTG", "CCAGTAGAAGTCCGACAACGTCAT",
				"CAGACTTGGTACGGTTGGGTAACT", "GGACGAAGAACTCAAGTCAAAGGC", "CTACTTACGAAGCTGAGGGACTGC", "ATGTCCCAGTTAGAGGAGGAAACA",
				"GCTTGCGATTGATGCTTAGTATCA", "ACCACAGGAGGACGATACAGAGAA", "CCACAGTGTCAACTAGAGCCTCTC", "TAGTTTGGATGACCAAGGATAGCC",
				"GGAGTTCGTCCAGAGAAGTACACG", "CTACGTGTAAGGCATACCTGCCAG", "CTTTCGTTGTTGACTCGACGGTAG", "AGTAGAAAGGGTTCCTTCCCACTC",
				"GATCCAACAGAGATGCCTTCAGTG", "GCTGTGTTCCACTTCATTCTCCTG", "GTGCAACTTTCCCACAGGTAGTTC", "CATCTGGAACGTGGTACACCTGTA",
				"ACTGGTGCAGCTTTGAACATCTAG", "ATGGACTTTGGTAACTTCCTGCGT", "GTTGAATGAGCCTACTGGGTCCTC", "TGAGAGACAAGATTGTTCGTGGAC",
				"AGATTCAGACCGTCTCATGCAAAG", "CAAGAGCTTTGACTAAGGAGCATG", "TGGAAGATGAGACCCTGATCTACG", "TCACTACTCAACAGGTGGCATGAA",
				"GCTAGGTCAATCTCCTTCGGAAGT", "CAGGTTACTCCTCCGTGAGTCTGA", "TCAATCAAGAAGGGAAAGCAAGGT", "CATGTTCAACCAAGGCTTCTATGG",
				"AGAGGGTACTATGTGCCTCAGCAC", "CACCCACACTTACTTCAGGACGTA", "TTCTGAAGTTCCTGGGTCTTGAAC", "GACAGACACCGTTCATCGACTTTC",
				"TTCTCAGTCTTCCTCCAGACAAGG", "CCGATCCTTGTGGCTTCTAACTTC", "GTTTGTCATACTCGTGTGCTCACC", "GAATCTAAGCAAACACGAAGGTGG",
				"TACAGTCCGAGCCTCATGTGATCT", "ACCGAGATCCTACGAATGGAGTGT", "CCTGGGAGCATCAGGTAGTAACAG", "TAGCTGACTGTCTTCCATACCGAC",
				"AAGAAACAGGATGACAGAACCCTC", "TACAAGCATCCCAACACTTCCACT", "GACCATTGTGATGAACCCTGTTGT", "ATGCTTGTTACATCAACCCTGGAC",
				"CGACCTGTTTCTCAGGGATACAAC", "AACAACCGAACCTTTGAATCAGAA", "TCTCGGAGATAGTTCTCACTGCTG", "CGGATGAACATAGGATAGCGATTC",
				"CCTCATCTTGTGAAGTTGTTTCGG", "ACGGTATGTCGAGTTCCAGGACTA", "TGGCTTGATCTAGGTAAGGTCGAA", "GTAGTGGACCTAGAACCTGTGCCA",
				"AACGGAGGAGTTAGTTGGATGATC", "AGGTGATCCCAACAAGCGTAAGTA", "TACATGCTCCTGTTGTTAGGGAGG", "TCTTCTACTACCGATCCGAAGCAG",
				"ACAGCATCAATGTTTGGCTAGTTG", "GATGTAGAGGGTACGGTTTGAGGC", "GGCTCCATAGGAACTCACGCTACT", "TTGTGAGTGGAAAGATACAGGACC",
				"AGTTTCCATCACTTCAGACTTGGG", "GATTGTCCTCAAACTGCCACCTAC", "CCTGTCTGGAAGAAGAATGGACTT", "CTGAACGGTCATAGAGTCCACCAT",
			}},
			{"BC", {
				"AAGAAAGTTGTCGGTGTCTTTGTG", "TCGATTCCGTTTGTAGTCGTCTGT", "GAGTCTTGTGTCCCAGTTACCAGG", "TTCGGATTCTATCGTGTTTCCCTA",
				"CTTGTCCAGGGTTTGTGTAACCTT", "TTCTCGCAAAGGCAGAAAGTAGTC", "GTGTTACCGTGGGAATGAATCCTT", "TTCAGGGAACAAACCAAGTTACGT",
				"AACTAGGCACAGCGAGTCTTGGTT", "AAGCGTTGAAACCTTTGTCCTCTC", "GTTTCATCTATCGGAGGGAATGGA", "CAGGTAGAAAGAAGCAGAATCGGA",
				"AGAACGACTTCCATACTCGTGTGA", "AACGAGTCTCTTGGGACCCATAGA", "AGGTCTACCTCGCTAACACCACTG", "CGTCAACTGACAGTGGTTCGTACT",
				"ACCCTCCAGGAAAGTACCTCTGAT", "CCAAACCCAACAACCTAGATAGGC", "GTTCCTCGTGCAGTGTCAAGAGAT", "TTGCGTCCTGTTACGAGAACTCAT",
				"GAGCCTCTCATTGTCCGTTCTCTA", "ACCACTGCCATGTATCAAAGTACG", "CTTACTACCCAGTGAACCTCCTCG", "GCATAGTTCTGCATGATGGGTTAG",
				"GTAAGTTGGGTATGCAACGCAATG", "CATACAGCGACTACGCATTCTCAT", "CGACGGTTAGATTCACCTCTTACA", "TGAAACCTAAGAAGGCACCGTATC",
				"CTAGACACCTTGGGTTGACAGACC", "TCAGTGAGGATCTACTTCGACCCA", "TGCGTACAGCAATCAGTTACATTG", "CCAGTAGAAGTCCGACAACGTCAT",
				"CAGACTTGGTACGGTTGGGTAACT", "GGACGAAGAACTCAAGTCAAAGGC", "CTACTTACGAAGCTGAGGGACTGC", "ATGTCCCAGTTAGAGGAGGAAACA",
				"GCTTGCGATTGATGCTTAGTATCA", "ACCACAGGAGGACGATACAGAGAA", "CCACAGTGTCAACTAGAGCCTCTC", "TAGTTTGGATGACCAAGGATAGCC",
				"GGAGTTCGTCCAGAGAAGTACACG", "CTACGTGTAAGGCATACCTGCCAG", "CTTTCGTTGTTGACTCGACGGTAG", "AGTAGAAAGGGTTCCTTCCCACTC",
				"GATCCAACAGAGATGCCTTCAGTG", "GCTGTGTTCCACTTCATTCTCCTG", "GTGCAACTTTCCCACAGGTAGTTC", "CATCTGGAACGTGGTACACCTGTA",
				"ACTGGTGCAGCTTTGAACATCTAG", "ATGGACTTTGGTAACTTCCTGCGT", "GTTGAATGAGCCTACTGGGTCCTC", "TGAGAGACAAGATTGTTCGTGGAC",
				"AGATTCAGACCGTCTCATGCAAAG", "CAAGAGCTTTGACTAAGGAGCATG", "TGGAAGATGAGACCCTGATCTACG", "TCACTACTCAACAGGTGGCATGAA",
				"GCTAGGTCAATCTCCTTCGGAAGT", "CAGGTTACTCCTCCGTGAGTCTGA", "TCAATCAAGAAGGGAAAGCAAGGT", "CATGTTCAACCAAGGCTTCTATGG",
				"AGAGGGTACTATGTGCCTCAGCAC", "CACCCACACTTACTTCAGGACGTA", "TTCTGAAGTTCCTGGGTCTTGAAC", "GACAGACACCGTTCATCGACTTTC",
				"TTCTCAGTCTTCCTCCAGACAAGG", "CCGATCCTTGTGGCTTCTAACTTC", "GTTTGTCATACTCGTGTGCTCACC", "GAATCTAAGCAAACACGAAGGTGG",
				"TACAGTCCGAGCCTCATGTGATCT", "ACCGAGATCCTACGAATGGAGTGT", "CCTGGGAGCATCAGGTAGTAACAG", "TAGCTGACTGTCTTCCATACCGAC",
				"AAGAAACAGGATGACAGAACCCTC", "TACAAGCATCCCAACACTTCCACT", "GACCATTGTGATGAACCCTGTTGT", "ATGCTTGTTACATCAACCCTGGAC",
				"CGACCTGTTTCTCAGGGATACAAC", "AACAACCGAACCTTTGAATCAGAA", "TCTCGGAGATAGTTCTCACTGCTG", "CGGATGAACATAGGATAGCGATTC",
				"CCTCATCTTGTGAAGTTGTTTCGG", "ACGGTATGTCGAGTTCCAGGACTA", "TGGCTTGATCTAGGTAAGGTCGAA", "GTAGTGGACCTAGAACCTGTGCCA",
				"AACGGAGGAGTTAGTTGGATGATC", "AGGTGATCCCAACAAGCGTAAGTA", "TACATGCTCCTGTTGTTAGGGAGG", "TCTTCTACTACCGATCCGAAGCAG",
				"ACAGCATCAATGTTTGGCTAGTTG", "GATGTAGAGGGTACGGTTTGAGGC", "GGCTCCATAGGAACTCACGCTACT", "TTGTGAGTGGAAAGATACAGGACC",
				"AGTTTCCATCACTTCAGACTTGGG", "GATTGTCCTCAAACTGCCACCTAC", "CCTGTCTGGAAGAAGAATGGACTT", "CTGAACGGTCATAGAGTCCACCAT",
			}},
	};
	return seqs;
}

//the barcode flanks of an ONT barcoding kit as adaptors with the barcode as N's (dorado's kit arrangements), in the orientation they read at the start of a read
struct OntBarcodeKit{
	std::string barcodePrefix_;
	std::string forwardAdaptor_;
	std::string reverseAdaptor_;
};

OntBarcodeKit getOntBarcodeKit(const std::string & kitName){
	const std::string barcodeNs(24, 'N');
	const VecStr nativeKits{"SQK-NBD114-96", "SQK-NBD114-24", "EXP-NBD114", "EXP-NBD114-24", "EXP-NBD196", "EXP-NBD104", "SQK-NBD111-96", "SQK-NBD111-24"};
	const VecStr pcrKits{"EXP-PBC096", "EXP-PBC001"};
	if(njh::in(kitName, nativeKits)){
		//the native barcode is ligated on as is, so both ends carry the same flanks
		std::string adaptor = "ATTGCTAAGGTTAA" + barcodeNs + "CAGCACCT";
		return OntBarcodeKit{"NB", adaptor, adaptor};
	}
	if(njh::in(kitName, pcrKits)){
		//the barcode PCR primers bind the two different tails put on by the first round primers
		return OntBarcodeKit{"BC",
			"GGTGCTG" + barcodeNs + "TTAACCTTTCTGTTGGTGCTGATATTGC",
			"GGTGCTG" + barcodeNs + "TTAACCTACTTGCCTGTCGCTCTATCTTC"};
	}
	std::stringstream ss;
	ss << __PRETTY_FUNCTION__ << ", error " << "unknown kit: " << kitName << ", options are: "
			<< njh::conToStr(nativeKits, ", ") << ", " << njh::conToStr(pcrKits, ", ") << "\n";
	throw std::runtime_error{ss.str()};
}

//the kit barcode for a name like barcode01, NB01, BC01 or 1
std::string getOntBarcodeSeq(const OntBarcodeKit & kit, const std::string & barcodeName){
	const auto & seqs = getOntBarcodeSeqs().at(kit.barcodePrefix_);
	auto numStart = barcodeName.find_last_not_of("0123456789") + 1;
	if(numStart >= barcodeName.size()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "couldn't get a barcode number from: " << barcodeName << "\n";
		throw std::runtime_error{ss.str()};
	}
	auto barcodeNum = njh::StrToNumConverter::stoToNum<uint32_t>(barcodeName.substr(numStart));
	if(barcodeNum < 1 || barcodeNum > seqs.size()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "barcode number for " << barcodeName << " should be between 1 and " << seqs.size() << "\n";
		throw std::runtime_error{ss.str()};
	}
	return seqs[barcodeNum - 1];
}

struct IlluminaDualBarCode {
	IlluminaDualBarCode(std::string  samp, const std::string& fwBar5_3,
	                    const std::string& revBar5_3) : sampleName_(std::move(samp)), fwBar_("fw", fwBar5_3),
	                                                    revBar_("rev", revBar5_3) {
	}

	std::string sampleName_;
	seqInfo fwBar_;
	seqInfo revBar_;
};

/**@brief read the barcodes for each sample
 *
 * A sample sheet with the columns sample_name, fw and rev (barcode sequences), or for an ONT kit either sample_name and barcode (kit barcode names like barcode01/NB01/BC01)
 * or a two column file without a header of barcode name then sample name (dorado name key)
 */
std::vector<IlluminaDualBarCode> readBarcodeSampleSheet(const bfs::path & sampleSheet, const OntBarcodeKit * kit){
	std::vector<std::pair<std::string, std::pair<std::string, std::string>>> rows;
	table barcodesTab(sampleSheet, "\t", true);
	if(barcodesTab.containsColumns(VecStr{"sample_name", "fw", "rev"})){
		for(const auto & row : barcodesTab){
			rows.emplace_back(row[barcodesTab.getColPos("sample_name")],
					std::make_pair(row[barcodesTab.getColPos("fw")], row[barcodesTab.getColPos("rev")]));
		}
	}else if(nullptr != kit && barcodesTab.containsColumns(VecStr{"sample_name", "barcode"})){
		for(const auto & row : barcodesTab){
			auto bar = getOntBarcodeSeq(*kit, row[barcodesTab.getColPos("barcode")]);
			rows.emplace_back(row[barcodesTab.getColPos("sample_name")], std::make_pair(bar, bar));
		}
	}else if(nullptr != kit && 2 == barcodesTab.columnNames_.size()){
		table nameKey(sampleSheet, "\t", false);
		for(const auto & row : nameKey){
			auto bar = getOntBarcodeSeq(*kit, row[0]);
			rows.emplace_back(row[1], std::make_pair(bar, bar));
		}
	}else{
		barcodesTab.checkForColumnsThrow(VecStr{"sample_name", "fw", "rev"}, __PRETTY_FUNCTION__);
	}
	std::vector<IlluminaDualBarCode> barcodes;
	for(const auto & row : rows){
		auto bar = IlluminaDualBarCode{row.first, row.second.first, row.second.second};
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
		barcodes.emplace_back(bar);
	}
	return barcodes;
}

/**@brief demultiplex reads on a barcoded adaptor at each end, shared by extractByIlluminaAaptors and extractByNanoporeBarcodes
 *
 * @param inputCommands the commands
 * @param nanopore set up for ONT barcoding kits (--kit), reads whose 3' end is cut short, and kit barcodes not in the sample sheet
 */
int runExtractByAdaptors(const njh::progutils::CmdArgs & inputCommands, bool nanopore) {
	AdaptorExtractPars pars;
	bfs::path barcodeSampleSheet = "";
	//when not 0, a barcode within this many matches of the best match is just as likely, see getBestBarSamples
	uint32_t barcodeMatchMargin = 0;
	std::string kitName;
	std::optional<OntBarcodeKit> kit;

	seqUtilsExtractSetUp setUp(inputCommands);
	if(nanopore){
		setUp.description_ = "Demultiplex nanopore reads on the ONT barcodes at both ends of the read, allowing the barcode at the end of the read to be cut short";
		pars.searchWithin_ = 175;
		pars.minOutputLen_ = 100;
		pars.extractPars_.minTruncatedBackBarcodeBases_ = 14;
		//with the margin a lower cut off lets primer/insert sequence match several barcodes and leaves the end undecided
		pars.barIdentityCutOff_ = .70;
		barcodeMatchMargin = 3;
		setUp.setOption(kitName, "--kit", "ONT barcoding kit, sets the adaptors and the barcode names for the sample sheet, e.g. SQK-NBD114-96 or EXP-PBC096", true);
		if(!kitName.empty()){
			kit = getOntBarcodeKit(kitName);
			pars.forwardAdaptor_ = kit->forwardAdaptor_;
			pars.reverseAdaptor_ = kit->reverseAdaptor_;
		}
	}
	setAdaptorExtractOptions(setUp, pars);
	setUp.setOption(barcodeMatchMargin, "--barcodeMatchMargin",
			"A barcode needs this many more matching bases than any other to be called by itself, 0 means just take the best (ties broken by alignment score)");
	setUp.setOption(barcodeSampleSheet, "--illuminaBarcodeSampleSheet,--barcodeSampleSheet",
			"Barcode Sample Sheet, columns sample_name, fw, rev (barcode sequences); with --kit can instead be sample_name and barcode (e.g. barcode01) or a dorado name key (barcode name then sample, no header)", true);

	setUp.processVerbose();
	setUp.processDebug();
	setUp.processReadInNames(true);
	setUp.processDirectoryOutputName(true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	checkAdaptorHasBarcode(pars.forwardAdaptor_, "--illuminaForward");
	checkAdaptorHasBarcode(pars.reverseAdaptor_, "--illuminaReverse");

	//read in barcodes
	//check for same barcodes
	auto barcodes = readBarcodeSampleSheet(barcodeSampleSheet, kit ? &(*kit) : nullptr);
	//kit barcodes no sample used still compete for reads so that a read from one of them isn't given to the closest sample
	std::unordered_set<std::string> unusedBarcodeNames;
	if(kit){
		const auto & kitSeqs = getOntBarcodeSeqs().at(kit->barcodePrefix_);
		for(uint32_t pos = 0; pos < kitSeqs.size(); ++pos){
			bool used = std::any_of(barcodes.begin(), barcodes.end(),
					[&kitSeqs, &pos](const IlluminaDualBarCode & bar){ return bar.fwBar_.seq_ == kitSeqs[pos];});
			if(!used){
				auto name = njh::pasteAsStr("unused_barcode", njh::leftPadNumStr<uint32_t>(pos + 1, 10));
				barcodes.emplace_back(IlluminaDualBarCode{name, kitSeqs[pos], kitSeqs[pos]});
				unusedBarcodeNames.emplace(name);
			}
		}
	}
	uint64_t maxBarLen = 0;
	for(const auto & bar : barcodes){
		readVec::getMaxLength(bar.fwBar_, maxBarLen);
		readVec::getMaxLength(bar.revBar_, maxBarLen);
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
	auto scoringBar = substituteMatrix::createDegenScoreMatrixLessN(2,-2);
	aligner alignerBarObj(maxBarLen, gapParsBar, scoringBar, true);

	AdaptorPairFinder finder(pars, setUp.pars_.debug_);

	SeqInput seqReader(setUp.pars_.ioOptions_);
	seqReader.openIn();
	seqInfo seq;


	OutputStream out(njh::files::make_path(setUp.pars_.directoryName_, "info.tab.txt"));
	out << "seqName\tseqLen\tfwd5_3_score\tfwd5_3_bar\trev3_5_score\trev3_5_bar\trev5_3_score\trev5_3_bar\tfwd3_5_score\tfwd3_5_bar" << std::endl;
	OutputStream outBest(njh::files::make_path(setUp.pars_.directoryName_, "infoBest.tab.txt"));
	outBest << "seqName\tdirection\tfbarcode\trbarcode" << std::endl;

	OutputStream outPass(njh::files::make_path(setUp.pars_.directoryName_, "infoPass.tab.txt"));

	outPass << "seqName\tdirection\tfwd5_3_score\tfwd5_3_bar\trev3_5_score\trev3_5_bar\trev5_3_score\trev5_3_bar\tfwd3_5_score\tfwd3_5_bar" << std::endl;

	std::map<std::string, uint32_t> sampleCounts;

	uint32_t fullPass = 0;
	uint32_t sampleDetermined = 0;

	//best matching samples only depend on the observed barcode, so cache them per observed barcode for each side
	std::unordered_map<std::string, VecStr> fwdBarSampleCache;
	std::unordered_map<std::string, VecStr> revBarSampleCache;
	auto getBestBarSamples = [&barcodes,&alignerBarObj,&pars,&barcodeMatchMargin,&fwdBarSampleCache,&revBarSampleCache](
			const std::string & observedBar, bool fwd) -> const VecStr & {
		auto & cache = fwd ? fwdBarSampleCache : revBarSampleCache;
		auto cached = cache.find(observedBar);
		if(cache.end() != cached){
			return cached->second;
		}
		VecStr bestSample{"Undetermined"};
		uint32_t bestIds = 0;
		int32_t bestAlnScore = 0;
		//identities, sample and barcode of every barcode that passes the identity cut off, for the margin
		std::vector<std::tuple<uint32_t, std::string, std::string>> passingIds;
		seqInfo obsBarSeq(fwd ? "determined-fwd" : "determined-rev", observedBar);
		for(const auto & bar : barcodes){
			const auto & expectedBar = fwd ? bar.fwBar_ : bar.revBar_;
			alignerBarObj.alignCacheGlobal(expectedBar, obsBarSeq);
			alignerBarObj.profilePrimerAlignment(expectedBar, obsBarSeq);
			if(alignerBarObj.comp_.distances_.eventBasedIdentity_ >=pars.barIdentityCutOff_){
				passingIds.emplace_back(alignerBarObj.comp_.distances_.ref_.identities_, bar.sampleName_, expectedBar.seq_);
				if(alignerBarObj.comp_.distances_.ref_.identities_ > bestIds){
					bestSample.clear();
					bestAlnScore = alignerBarObj.parts_.score_;
					bestIds = alignerBarObj.comp_.distances_.ref_.identities_;
					bestSample.emplace_back(bar.sampleName_);
				}else if(alignerBarObj.comp_.distances_.ref_.identities_ == bestIds){
					if(alignerBarObj.parts_.score_ > bestAlnScore){
						bestSample.clear();
						bestAlnScore = alignerBarObj.parts_.score_;
						bestIds = alignerBarObj.comp_.distances_.ref_.identities_;
						bestSample.emplace_back(bar.sampleName_);
					}else if(alignerBarObj.parts_.score_ == bestAlnScore){
						bestSample.emplace_back(bar.sampleName_);
					}
				}
			}
		}
		if(barcodeMatchMargin > 0 && !passingIds.empty()){
			//this end only counts when a single barcode sequence beats every other by the margin (it can still be shared by several samples),
			//an end with no real barcode (e.g. cut off, or primer sequence) matches many barcodes about equally and would otherwise let the other end decide alone
			std::set<std::string> closeBarcodes;
			VecStr closeSamples;
			for(const auto & passing : passingIds){
				if(std::get<0>(passing) + barcodeMatchMargin > bestIds){
					closeBarcodes.emplace(std::get<2>(passing));
					closeSamples.emplace_back(std::get<1>(passing));
				}
			}
			bestSample = 1 == closeBarcodes.size() ? closeSamples : VecStr{"Undetermined"};
		}
		return cache.emplace(observedBar, std::move(bestSample)).first->second;
	};

	while(seqReader.readNextRead(seq)){
		auto res = finder.process(seq, out);
		if(AdaptorReadOutcome::pass != res.outcome_){
			continue;
		}
		const auto & processRes = res.processRes_;
		const auto & forwardBar = res.forwardBar_;
		const auto & reverseBar = res.reverseBar_;
		++fullPass;
		std::string determinedSample = "Undetermined";

		{
			VecStr fwdBestSample{"Undetermined"};
			VecStr revBestSample{"Undetermined"};

			if(forwardBar.size() < maxBarLen &&
					reverseBar.size() < maxBarLen &&
					forwardBar.size() > 5 &&
					reverseBar.size() > 5){
				fwdBestSample = getBestBarSamples(forwardBar, true);
				revBestSample = getBestBarSamples(reverseBar, false);
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
		seqMeta.addMeta("direction", res.direction_);
		seqMeta.addMeta("fprimerStart", processRes.fprimerStart_);
		seqMeta.addMeta("forwardBar", forwardBar);
		seqMeta.addMeta("reverseBar", reverseBar);
		seqMeta.addMeta("fullBar", njh::pasteAsStr(forwardBar, "+", reverseBar));
		seqMeta.addMeta("rprimerStart", processRes.rprimerStart_);
		seqMeta.addMeta("sample", determinedSample);
		seq.name_ += seqMeta.createMetaName();
		outputs.add(determinedSample, seq);
		++sampleCounts[determinedSample];
		if("Undetermined" != determinedSample && 0 == unusedBarcodeNames.count(determinedSample)) {
			++sampleDetermined;
		}
		AdaptorPairFinder::writePassLine(outPass, seq.name_, res);
		outBest << seq.name_
				<< "\t" << res.direction_
				<< "\t" << forwardBar
				<< "\t" << reverseBar << "\n";
	}
	finder.writeBestDirectionCounts(setUp.pars_.directoryName_);
	finder.writeFailedCounts(setUp.pars_.directoryName_);
	{
		std::vector<std::pair<std::string, uint32_t>> afterPasses{
			{"fullPass", fullPass},
			{"Undetermined", sampleCounts["Undetermined"]},
			{"sampleDetermined", sampleDetermined}};
		if(kit){
			uint32_t unusedBarcodeReads = 0;
			for(const auto & sampleCount : sampleCounts){
				if(unusedBarcodeNames.count(sampleCount.first)){
					unusedBarcodeReads += sampleCount.second;
				}
			}
			afterPasses.emplace_back("unusedBarcode", unusedBarcodeReads);
		}
		finder.writeConditionCounts(setUp.pars_.directoryName_, afterPasses);
	}
	{
		OutputStream sampleCountsOut(njh::files::make_path(setUp.pars_.directoryName_, "sampleCounts.tab.txt"));
		sampleCountsOut << "sample\tcount\tfrac" << std::endl;
		for(const auto & sampleCount : sampleCounts){
			sampleCountsOut << sampleCount.first
					<< "\t" << sampleCount.second
					<< "\t" << sampleCount.second/static_cast<double>(finder.passes()) << std::endl;
		}
	}

	return 0;
}

int seqUtilsExtractRunner::extractByIlluminaAaptors(const njh::progutils::CmdArgs & inputCommands) {
	return runExtractByAdaptors(inputCommands, false);
}

int seqUtilsExtractRunner::extractByNanoporeBarcodes(const njh::progutils::CmdArgs & inputCommands) {
	return runExtractByAdaptors(inputCommands, true);
}


}  // namespace njhseq


