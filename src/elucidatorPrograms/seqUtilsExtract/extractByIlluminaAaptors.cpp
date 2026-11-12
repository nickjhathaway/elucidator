/*
 * extractByIlluminaAaptors.cpp
 *
 *  Created on: Nov 28, 2019
 *      Author: nicholashathaway
 */




#include "seqUtilsExtractRunner.hpp"

#include <njhseq/IO/SeqIO.h>
#include <njhseq/alignment/aligner/aligner.hpp>
#include <njhcpp/concurrency/concurrencyUtils.hpp>
#include <optional>
#include <utility>
#include <mutex>
#include <sstream>


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

//gap scoring for aligning an adaptor inside a longer window of the read, the adaptor's own end gaps are free
gapScoringParameters genAdaptorGapPars(){
	gapScoringParameters gapPars(5,1,0,0,0,0);
	gapPars.gapLeftRefOpen_ = 5;
	gapPars.gapLeftRefExtend_ = 1;
	gapPars.gapRightRefOpen_ = 5;
	gapPars.gapRightRefExtend_ = 1;
	return gapPars;
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
	//for mismatchDirections, whether the one adaptor found at both ends is the forward adaptor; forwardBar_ is then the barcode at the front of the read and reverseBar_ the one at the back, both in the orientation of that adaptor
	bool sameAdaptorIsForward_ = true;
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
			//the best front and best back are the same adaptor, e.g. a read that folded back on itself carries one adaptor at both ends
			ret.sameAdaptorIsForward_ = "fwd" == bestFront;
			ret.direction_ = ret.sameAdaptorIsForward_ ? "fwdAdaptorBothEnds" : "revAdaptorBothEnds";
			ret.forwardBar_ = ret.sameAdaptorIsForward_ ? fwd5_3_info.barcode_ : rev5_3_info.barcode_;
			ret.reverseBar_ = seqUtil::reverseComplement(ret.sameAdaptorIsForward_ ? fwd3_5_info.barcode_ : rev3_5_info.barcode_, "DNA");
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
		trimToBetweenAdaptors(seq, processRes, reverse);
		if(reverse){
			++passRev_;
		}else{
			++passFwd_;
		}
		ret.outcome_ = AdaptorReadOutcome::pass;
		return ret;
	}

	/**@brief trim a mismatchDirections read to what is between its two adaptors, holding it to the same checks as a read with its direction called
	 *
	 * A read with the reverse adaptor at the front is reverse complemented, as any other read with that adaptor at the front is
	 *
	 * @param seq the read, only changed if the adaptors pass the checks
	 * @param res the result from process() for this read
	 * @return whether the adaptors passed
	 */
	bool trimMismatchDirectionRead(seqInfo & seq, AdaptorReadRes & res) const {
		const auto & frontInfo = res.sameAdaptorIsForward_ ? *res.fwd5_3_ : *res.rev5_3_;
		const auto & backInfo = res.sameAdaptorIsForward_ ? *res.fwd3_5_ : *res.rev3_5_;
		if(frontInfo.comp_.distances_.eventBasedIdentity_ < pars_.primerIdentityCutOff_
				|| backInfo.comp_.distances_.eventBasedIdentity_ < pars_.primerIdentityCutOff_){
			return false;
		}
		//processSeqForAdaptors can change the case of the primers, so leave seq alone until this is known to pass
		seqInfo trimmed = seq;
		uint32_t startOfBackSeq = seq.seq_.size() - std::min<uint32_t>(pars_.searchWithin_, len(seq));
		res.processRes_ = processSeqForAdaptors(trimmed, frontInfo, backInfo, startOfBackSeq, pars_.extractPars_);
		if(!res.processRes_.fwdPass_ || !res.processRes_.revPass_ || res.processRes_.fprimerEnd_ >= res.processRes_.rprimerStart_){
			return false;
		}
		trimToBetweenAdaptors(trimmed, res.processRes_, !res.sameAdaptorIsForward_);
		seq = trimmed;
		return true;
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
	 * @param afterMismatchDirections extra conditions to report right after mismatchDirections
	 */
	void writeConditionCounts(const bfs::path & dir, const std::vector<std::pair<std::string, uint32_t>> & afterPasses,
			const std::vector<std::pair<std::string, uint32_t>> & afterMismatchDirections = {}) const {
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
		for(const auto & condition : afterMismatchDirections){
			conditionCount << condition.first << "\t" << condition.second << "\t" << condition.second/total_ << std::endl;
		}
		conditionCount << "failedOverlapCheck" << "\t" << failedOverlapCheck_ << "\t" << failedOverlapCheck_/total_ << std::endl;
	}

	uint32_t passes() const {
		return passFwd_ + passRev_;
	}

	/**@brief add another finder's counts to this one's, for combining the finders of several threads
	 */
	void addCounts(const AdaptorPairFinder & other){
		for(const auto & front : other.bestCounts_){
			for(const auto & back : front.second){
				bestCounts_[front.first][back.first] += back.second;
			}
		}
		total_ += other.total_;
		smallSeq_ += other.smallSeq_;
		largeSeq_ += other.largeSeq_;
		failed_fwd5_3_Identity_ += other.failed_fwd5_3_Identity_;
		failed_fwd3_5_Identity_ += other.failed_fwd3_5_Identity_;
		failed_rev5_3_Identity_ += other.failed_rev5_3_Identity_;
		failed_rev3_5_Identity_ += other.failed_rev3_5_Identity_;
		failed_fwd5_3_PrimerAln_ += other.failed_fwd5_3_PrimerAln_;
		failed_fwd3_5_PrimerAln_ += other.failed_fwd3_5_PrimerAln_;
		failed_rev5_3_PrimerAln_ += other.failed_rev5_3_PrimerAln_;
		failed_rev3_5_PrimerAln_ += other.failed_rev3_5_PrimerAln_;
		mismatchDirections_ += other.mismatchDirections_;
		failedIdentityReads_ += other.failedIdentityReads_;
		failedPrimerAlnReads_ += other.failedPrimerAlnReads_;
		failedOverlapCheck_ += other.failedOverlapCheck_;
		passFwd_ += other.passFwd_;
		passRev_ += other.passRev_;
	}

private:
	void trimToBetweenAdaptors(seqInfo & seq, const ProcessedForIlluminaAdaptorRes & processRes, bool reverse) const {
		if (pars_.extractPars_.keepPrimer_) {
			seq = seq.getSubRead(processRes.fprimerStart_, processRes.rprimerEnd_ - processRes.fprimerStart_);
		} else {
			seq = seq.getSubRead(processRes.fprimerEnd_, processRes.rprimerStart_ - processRes.fprimerEnd_);
		}
		if(reverse){
			seq.reverseComplementRead(false, true);
		}
	}

	static gapScoringParameters genGapPars(){
		return genAdaptorGapPars();
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

//what the barcode observed at one end of a read was called as
struct BarEndCall {
	VecStr samples_{"Undetermined"};
	//identity of the called barcode, 0 when the end wasn't called
	double identity_ = 0;
	//matching bases of the called barcode, of the next best different barcode sequence, and the length of the called barcode
	uint32_t matches_ = 0;
	uint32_t nextBestMatches_ = 0;
	uint32_t barcodeLen_ = 0;

	bool undetermined() const {
		return samples_.empty() || njh::in(std::string("Undetermined"), samples_);
	}
};

/**@brief call the sample(s) an observed barcode belongs to by aligning it against every sample's barcode
 *
 * best matching samples only depend on the observed barcode, so they are cached per observed barcode for each side
 */
class BarcodeEndCaller {
public:
	/**
	 * @param barcodes the sample barcodes, has to outlive the caller
	 * @param barIdentityCutOff the identity a barcode needs to be considered
	 * @param barcodeMatchMargin when not 0, a barcode within this many matches of the best match is just as likely
	 * @param cacheAlignments whether to use the aligner's alignment cache
	 */
	BarcodeEndCaller(const std::vector<IlluminaDualBarCode> & barcodes, double barIdentityCutOff, uint32_t barcodeMatchMargin, bool cacheAlignments) :
			barcodes_(barcodes), barIdentityCutOff_(barIdentityCutOff), barcodeMatchMargin_(barcodeMatchMargin), cacheAlignments_(cacheAlignments),
			maxBarLen_(genMaxBarLen(barcodes)),
			alignerBarObj_(maxBarLen_, gapScoringParameters(5,1,5,1,5,1), substituteMatrix::createDegenScoreMatrixLessN(2,-2), true) {
	}

	//whether an observed barcode of this size can be aligned, too short means there is nothing left of the barcode
	bool callable(const std::string & observedBar) const {
		return observedBar.size() < maxBarLen_ && observedBar.size() > 5;
	}

	const BarEndCall & call(const std::string & observedBar, bool fwd){
		auto & cache = fwd ? fwdBarSampleCache_ : revBarSampleCache_;
		auto cached = cache.find(observedBar);
		if(cache.end() != cached){
			return cached->second;
		}
		VecStr bestSample{"Undetermined"};
		uint32_t bestIds = 0;
		int32_t bestAlnScore = 0;
		//identities, sample and barcode of every barcode that passes the identity cut off, for the margin
		std::vector<std::tuple<uint32_t, std::string, std::string, double>> passingIds;
		//identities of all barcodes, for how close the next best barcode was
		std::vector<std::pair<uint32_t, const std::string *>> allIds;
		seqInfo obsBarSeq(fwd ? "determined-fwd" : "determined-rev", observedBar);
		for(const auto & bar : barcodes_){
			const auto & expectedBar = fwd ? bar.fwBar_ : bar.revBar_;
			//nanopore barcodes almost never repeat exactly, so caching their alignments only grows memory (11GB on 4.4M reads)
			if(cacheAlignments_){
				alignerBarObj_.alignCacheGlobal(expectedBar, obsBarSeq);
			}else{
				alignerBarObj_.alignRegGlobal(expectedBar, obsBarSeq);
			}
			alignerBarObj_.profilePrimerAlignment(expectedBar, obsBarSeq);
			allIds.emplace_back(alignerBarObj_.comp_.distances_.ref_.identities_, &expectedBar.seq_);
			if(alignerBarObj_.comp_.distances_.eventBasedIdentity_ >= barIdentityCutOff_){
				passingIds.emplace_back(alignerBarObj_.comp_.distances_.ref_.identities_, bar.sampleName_, expectedBar.seq_, alignerBarObj_.comp_.distances_.eventBasedIdentity_);
				if(alignerBarObj_.comp_.distances_.ref_.identities_ > bestIds){
					bestSample.clear();
					bestAlnScore = alignerBarObj_.parts_.score_;
					bestIds = alignerBarObj_.comp_.distances_.ref_.identities_;
					bestSample.emplace_back(bar.sampleName_);
				}else if(alignerBarObj_.comp_.distances_.ref_.identities_ == bestIds){
					if(alignerBarObj_.parts_.score_ > bestAlnScore){
						bestSample.clear();
						bestAlnScore = alignerBarObj_.parts_.score_;
						bestIds = alignerBarObj_.comp_.distances_.ref_.identities_;
						bestSample.emplace_back(bar.sampleName_);
					}else if(alignerBarObj_.parts_.score_ == bestAlnScore){
						bestSample.emplace_back(bar.sampleName_);
					}
				}
			}
		}
		if(barcodeMatchMargin_ > 0 && !passingIds.empty()){
			//this end only counts when a single barcode sequence beats every other by the margin (it can still be shared by several samples),
			//an end with no real barcode (e.g. cut off, or primer sequence) matches many barcodes about equally and would otherwise let the other end decide alone
			std::set<std::string> closeBarcodes;
			VecStr closeSamples;
			for(const auto & passing : passingIds){
				if(std::get<0>(passing) + barcodeMatchMargin_ > bestIds){
					closeBarcodes.emplace(std::get<2>(passing));
					closeSamples.emplace_back(std::get<1>(passing));
				}
			}
			bestSample = 1 == closeBarcodes.size() ? closeSamples : VecStr{"Undetermined"};
		}
		BarEndCall call;
		std::unordered_set<std::string> calledBarcodes;
		for(const auto & passing : passingIds){
			if(njh::in(std::get<1>(passing), bestSample)){
				call.identity_ = std::max(call.identity_, std::get<3>(passing));
				call.matches_ = std::max(call.matches_, std::get<0>(passing));
				call.barcodeLen_ = std::max<uint32_t>(call.barcodeLen_, std::get<2>(passing).size());
				calledBarcodes.emplace(std::get<2>(passing));
			}
		}
		for(const auto & ids : allIds){
			if(0 == calledBarcodes.count(*ids.second)){
				call.nextBestMatches_ = std::max(call.nextBestMatches_, ids.first);
			}
		}
		call.samples_ = std::move(bestSample);
		//keep the cache bounded, it's only a speed up
		if(cache.size() >= 1000000){
			cache.clear();
		}
		return cache.emplace(observedBar, std::move(call)).first->second;
	}

private:
	static uint64_t genMaxBarLen(const std::vector<IlluminaDualBarCode> & barcodes){
		uint64_t maxBarLen = 0;
		for(const auto & bar : barcodes){
			readVec::getMaxLength(bar.fwBar_, maxBarLen);
			readVec::getMaxLength(bar.revBar_, maxBarLen);
		}
		return maxBarLen * 2;
	}

	const std::vector<IlluminaDualBarCode> & barcodes_;
	double barIdentityCutOff_;
	uint32_t barcodeMatchMargin_;
	bool cacheAlignments_;
	uint64_t maxBarLen_;
	aligner alignerBarObj_;
	std::unordered_map<std::string, BarEndCall> fwdBarSampleCache_;
	std::unordered_map<std::string, BarEndCall> revBarSampleCache_;
};

//the sample(s) each end was called as and how the two ends relate
struct BarcodePairCall {
	std::string determinedSample_ = "Undetermined";
	std::string fwdCall_ = "Undetermined";
	std::string revCall_ = "Undetermined";
	std::string pairStatus_ = "bothEndsUndetermined";
};

const std::string mismatchedSamplesName = "mismatchedSamples";

/**@brief call the samples of the barcodes at the two ends of a read
 * @param caller the barcode caller
 * @param forwardBar the observed barcode of the first end
 * @param reverseBar the observed barcode of the second end
 * @param forwardBarIsFwd whether the first end's barcode is compared against the samples' fw barcodes rather than rev barcodes
 * @param reverseBarIsFwd same for the second end
 * @param mismatchConfirmIdentityCutOff both ends need this identity to count as mismatched, 0 means any two ends that call different samples count as mismatched
 */
BarcodePairCall callBarcodePair(BarcodeEndCaller & caller, const std::string & forwardBar, const std::string & reverseBar,
		bool forwardBarIsFwd, bool reverseBarIsFwd, double mismatchConfirmIdentityCutOff){
	BarcodePairCall ret;
	VecStr fwdBestSample{"Undetermined"};
	VecStr revBestSample{"Undetermined"};
	double fwdIdentity = 0;
	double revIdentity = 0;

	if(caller.callable(forwardBar) && caller.callable(reverseBar)){
		const auto & fwdEndCall = caller.call(forwardBar, forwardBarIsFwd);
		fwdBestSample = fwdEndCall.samples_;
		fwdIdentity = fwdEndCall.identity_;
		const auto & revEndCall = caller.call(reverseBar, reverseBarIsFwd);
		revBestSample = revEndCall.samples_;
		revIdentity = revEndCall.identity_;
	}
	bool fwdUndetermined = njh::in(std::string("Undetermined"),fwdBestSample) || fwdBestSample.empty();
	bool revUndetermined = njh::in(std::string("Undetermined"),revBestSample) || revBestSample.empty();
	if(!fwdUndetermined){
		ret.fwdCall_ = njh::conToStr(fwdBestSample, ",");
	}
	if(!revUndetermined){
		ret.revCall_ = njh::conToStr(revBestSample, ",");
	}
	//need one for each and need to match
	//if multiple in one or both, see if there is only one pair since barcodes can be reused several times
	if (fwdUndetermined || revUndetermined) {
		ret.determinedSample_ = "Undetermined";
		ret.pairStatus_ = fwdUndetermined && revUndetermined ? "bothEndsUndetermined" : "oneEndUndetermined";
	}else{
		if(fwdBestSample.size() == 1 && revBestSample.size() == 1 && fwdBestSample.front() == revBestSample.front()){
			ret.determinedSample_ = fwdBestSample.front();
			ret.pairStatus_ = "matched";
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
				ret.determinedSample_ = pairs.front();
				ret.pairStatus_ = "matched";
			}else if(pairs.empty()){
				//a weak call at either end is more often a misread than a read joining two samples, keep those out of the mismatch counts
				if(fwdIdentity >= mismatchConfirmIdentityCutOff && revIdentity >= mismatchConfirmIdentityCutOff){
					ret.determinedSample_ = mismatchedSamplesName;
					ret.pairStatus_ = "mismatched";
				}else{
					ret.determinedSample_ = "Undetermined";
					ret.pairStatus_ = "mismatchedUnconfirmed";
				}
			}else{
				ret.determinedSample_ = "Undetermined";
				ret.pairStatus_ = "ambiguousPair";
			}
		}
	}
	return ret;
}

//every read that passed both adaptors by the sample(s) each end was called as, most common first; the mismatched rows over
//mismatched plus matched give the rate of reads joining barcodes of two samples
void writeBarcodePairCounts(const bfs::path & dir, const std::map<std::tuple<std::string, std::string, std::string>, uint32_t> & barcodePairCounts, uint32_t fullPass){
	std::vector<std::pair<std::tuple<std::string, std::string, std::string>, uint32_t>> pairCounts(
			barcodePairCounts.begin(), barcodePairCounts.end());
	std::stable_sort(pairCounts.begin(), pairCounts.end(), [](const auto & left, const auto & right){
		return left.second > right.second;
	});
	OutputStream pairCountsOut(njh::files::make_path(dir, "barcodePairCounts.tab.txt"));
	pairCountsOut << "fwdSample\trevSample\tstatus\tcount\tfrac" << std::endl;
	for(const auto & pairCount : pairCounts){
		pairCountsOut << std::get<0>(pairCount.first)
				<< "\t" << std::get<1>(pairCount.first)
				<< "\t" << std::get<2>(pairCount.first)
				<< "\t" << pairCount.second
				<< "\t" << pairCount.second/static_cast<double>(fullPass) << std::endl;
	}
}

//per read info files are a line per read, so only written with --development
struct AdaptorInfoFiles {
	explicit AdaptorInfoFiles(const bfs::path & dir, bool development) : noInfoOut_(nullptr){
		if(development){
			out_ = std::make_unique<OutputStream>(njh::files::make_path(dir, "info.tab.txt"));
			*out_ << "seqName\tseqLen\tfwd5_3_score\tfwd5_3_bar\trev3_5_score\trev3_5_bar\trev5_3_score\trev5_3_bar\tfwd3_5_score\tfwd3_5_bar" << std::endl;
			outBest_ = std::make_unique<OutputStream>(njh::files::make_path(dir, "infoBest.tab.txt"));
			*outBest_ << "seqName\tdirection\tfbarcode\trbarcode" << std::endl;
			outPass_ = std::make_unique<OutputStream>(njh::files::make_path(dir, "infoPass.tab.txt"));
			*outPass_ << "seqName\tdirection\tfwd5_3_score\tfwd5_3_bar\trev3_5_score\trev3_5_bar\trev5_3_score\trev5_3_bar\tfwd3_5_score\tfwd3_5_bar" << std::endl;
		}
	}

	bool on() const {
		return nullptr != out_;
	}

	//lines a thread wrote for info.tab.txt
	void writeInfo(const std::string & lines){
		if(out_ && !lines.empty()){
			std::lock_guard<std::mutex> lock(mut_);
			*out_ << lines;
		}
	}

	void writePass(const std::string & seqName, const AdaptorReadRes & res){
		if(out_){
			std::lock_guard<std::mutex> lock(mut_);
			AdaptorPairFinder::writePassLine(*outPass_, seqName, res);
			*outBest_ << seqName
					<< "\t" << res.direction_
					<< "\t" << res.forwardBar_
					<< "\t" << res.reverseBar_ << "\n";
		}
	}

	std::unique_ptr<OutputStream> out_;
	std::unique_ptr<OutputStream> outBest_;
	std::unique_ptr<OutputStream> outPass_;
	std::ostream noInfoOut_;
	std::mutex mut_;
};

template<typename MAP>
void addMapCounts(MAP & into, const MAP & from){
	for(const auto & count : from){
		into[count.first] += count.second;
	}
}

//read counts the demultiplexers keep, each thread keeps its own and they're added together at the end
struct DemuxCounts {
	std::map<std::string, uint32_t> sampleCounts;
	std::map<std::string, uint32_t> oneEndSampleCounts;
	//fwd end call, rev end call, status of the pair -> read count
	std::map<std::tuple<std::string, std::string, std::string>, uint32_t> barcodePairCounts;
	std::map<std::string, uint32_t> pairStatusCounts;
	//adaptor at both ends, front call, back call, status of the pair -> read count
	std::map<std::tuple<std::string, std::string, std::string, std::string>, uint32_t> mismatchDirectionPairCounts;
	std::map<std::string, uint32_t> mismatchDirectionStatusCounts;
	//what the first pass made of the read -> reads assigned on one end's barcode
	std::map<std::string, uint32_t> oneEndFromCounts;
	std::map<std::string, uint32_t> oneEndUsedCounts;

	uint32_t mismatchDirectionsToSamples = 0;
	uint32_t fullPass = 0;
	uint32_t sampleDetermined = 0;
	uint32_t sampleDeterminedBothEnds = 0;
	uint32_t sampleDeterminedOneEnd = 0;
	uint32_t oneEndToUnusedBarcode = 0;
	uint32_t oneEndConflicting = 0;
	uint32_t oneEndCutAtInnerAdaptor = 0;
	uint32_t oneEndNoInsertBeforeInnerAdaptor = 0;

	void add(const DemuxCounts & other){
		addMapCounts(sampleCounts, other.sampleCounts);
		addMapCounts(oneEndSampleCounts, other.oneEndSampleCounts);
		addMapCounts(barcodePairCounts, other.barcodePairCounts);
		addMapCounts(pairStatusCounts, other.pairStatusCounts);
		addMapCounts(mismatchDirectionPairCounts, other.mismatchDirectionPairCounts);
		addMapCounts(mismatchDirectionStatusCounts, other.mismatchDirectionStatusCounts);
		addMapCounts(oneEndFromCounts, other.oneEndFromCounts);
		addMapCounts(oneEndUsedCounts, other.oneEndUsedCounts);
		mismatchDirectionsToSamples += other.mismatchDirectionsToSamples;
		fullPass += other.fullPass;
		sampleDetermined += other.sampleDetermined;
		sampleDeterminedBothEnds += other.sampleDeterminedBothEnds;
		sampleDeterminedOneEnd += other.sampleDeterminedOneEnd;
		oneEndToUnusedBarcode += other.oneEndToUnusedBarcode;
		oneEndConflicting += other.oneEndConflicting;
		oneEndCutAtInnerAdaptor += other.oneEndCutAtInnerAdaptor;
		oneEndNoInsertBeforeInnerAdaptor += other.oneEndNoInsertBeforeInnerAdaptor;
	}
};

//how the fixed bases of an adaptor aligned, the barcode N's match anything so the whole adaptor's identity passes on nearly any sequence
struct AdaptorFlankStats {
	uint32_t matches_ = 0;
	uint32_t columns_ = 0;

	double identity() const {
		return 0 == columns_ ? 0 : matches_/static_cast<double>(columns_);
	}
};

//only the columns where the read and the adaptor overlap count, so an adaptor cut off by the end of the read isn't held against it
AdaptorFlankStats getAdaptorFlankStats(const PrimerAlnInfo & info){
	AdaptorFlankStats ret;
	const auto & seqAln = info.seqAln_.seq_;
	const auto & adaptorAln = info.primerAln_.seq_;
	auto firstN = adaptorAln.find('N');
	auto lastN = adaptorAln.rfind('N');
	auto start = std::max(seqAln.find_first_not_of('-'), adaptorAln.find_first_not_of('-'));
	auto stop = std::min(seqAln.find_last_not_of('-'), adaptorAln.find_last_not_of('-'));
	for(size_t col = start; col <= stop && col < adaptorAln.size(); ++col){
		if(col >= firstN && col <= lastN){
			continue;
		}
		++ret.columns_;
		if(seqAln[col] == adaptorAln[col]){
			++ret.matches_;
		}
	}
	return ret;
}

//settings for splitting reads at adaptors found inside them, i.e. two library molecules that were read as one
struct InnerAdaptorSplitPars {
	bool split_ = true;
	//k-mers of the adaptors (with each barcode filled in) seed the search, only where they cluster is the adaptor aligned
	uint32_t seedKmerLen_ = 10;
	uint32_t minSeedHits_ = 2;
	//the adaptor's fixed bases (not the barcode) have to match this well for it to count
	uint32_t flankMinMatches_ = 12;
	double flankIdentityCutOff_ = .75;
	//pieces on either side of a junction that share this fraction of k-mers with each other's reverse complement are one molecule folded back on itself, those are left whole
	double foldbackKmerFrac_ = .2;
	//exact k-mers of the adaptor's inner flank (the fixed bases between the barcode and the insert, e.g. the first round primer tail), these find partial or bare tails
	//that don't pass the adaptor check, used to cut a read assigned on one end; 0 turns it off. None of the panel targets has one of these 12-mers (2.4 per Mb of Pf3D7)
	uint32_t innerFlankKmerLen_ = 12;
};

void setInnerAdaptorSplitOptions(seqUtilsExtractSetUp & setUp, InnerAdaptorSplitPars & pars){
	bool noSplit = false;
	setUp.setOption(noSplit, "--noInnerAdaptorSplit",
			"Don't split reads at adaptors found inside them; by default a read holding more than one library molecule is cut at the inner adaptors so each piece keeps the adaptor and barcode of its own molecule, and each piece is demultiplexed by itself (named <read>_split<start>-<end>)");
	pars.split_ = !noSplit;
	setUp.setOption(pars.seedKmerLen_, "--innerAdaptorSeedKmerLen", "Inner adaptor split: k-mer length of the seeds taken from the adaptors and barcodes");
	setUp.setOption(pars.minSeedHits_, "--innerAdaptorMinSeedHits", "Inner adaptor split: seed hits on the same diagonal needed before aligning the adaptor there");
	setUp.setOption(pars.flankMinMatches_, "--innerAdaptorFlankMinMatches", "Inner adaptor split: matching fixed adaptor bases (not barcode) needed to call the adaptor");
	setUp.setOption(pars.flankIdentityCutOff_, "--innerAdaptorFlankIdentityCutOff", "Inner adaptor split: identity of the fixed adaptor bases needed to call the adaptor");
	setUp.setOption(pars.foldbackKmerFrac_, "--innerAdaptorFoldbackKmerFrac",
			"Inner adaptor split: pieces sharing at least this fraction of k-mers with the reverse complement of the piece next to them are a fold back of one molecule and are left joined");
	setUp.setOption(pars.innerFlankKmerLen_, "--innerAdaptorTailKmerLen",
			"Length of the exact k-mers of the adaptor's inner flank (the fixed bases between the barcode and the insert) used to find partial adaptor tails in a read assigned on one end's barcode, the read is cut back to before them; 0 turns this off");
	if(0 == pars.seedKmerLen_ || pars.seedKmerLen_ > 31){
		setUp.failed_ = true;
		setUp.addWarning(njh::pasteAsStr("--innerAdaptorSeedKmerLen should be between 1 and 31, not ", pars.seedKmerLen_));
	}
}

/**@brief split reads that hold more than one library molecule at the adaptors found inside them
 *
 * A molecule reads as [adaptor] insert [reverse complement of an adaptor], so an adaptor in its start orientation inside a read is where a new molecule starts
 * and one in its end orientation is where a molecule stops. The read is cut there, each piece keeping the adaptor (and barcode) of its own molecule,
 * so assigning each piece on its own barcodes keeps every insert with the barcode closest to it.
 * Adaptors closer than the minimum output length to either end of the read are left to the regular search at the ends.
 */
class InnerAdaptorSplitter {
public:
	/**@brief add another splitter's counts to this one's, for combining the splitters of several threads
	 */
	void addCounts(const InnerAdaptorSplitter & other){
		readsLongEnough_ += other.readsLongEnough_;
		readsWithSeeds_ += other.readsWithSeeds_;
		windowsAligned_ += other.windowsAligned_;
		adaptorsFound_[0] += other.adaptorsFound_[0];
		adaptorsFound_[1] += other.adaptorsFound_[1];
		adaptorsNearEnd_ += other.adaptorsNearEnd_;
		foldbackJunctions_ += other.foldbackJunctions_;
		shortPiecesDropped_ += other.shortPiecesDropped_;
		readsSplit_ += other.readsSplit_;
		piecesFromSplits_ += other.piecesFromSplits_;
	}

	InnerAdaptorSplitter(const AdaptorExtractPars & extractPars, const std::vector<IlluminaDualBarCode> & barcodes, const InnerAdaptorSplitPars & pars) :
			pars_(pars), minPieceLen_(std::max<uint32_t>(1, extractPars.minOutputLen_)),
			maxAdaptorLen_(std::max(extractPars.forwardAdaptor_.size(), extractPars.reverseAdaptor_.size())),
			maxWindow_(2 * maxAdaptorLen_ + 2 * windowPad_),
			alignerObj_(std::max<uint64_t>(500, maxWindow_ + 10), genAdaptorGapPars(), substituteMatrix::createDegenScoreMatrixLessN(2,-2), false) {
		bool sameAdaptor = extractPars.forwardAdaptor_ == extractPars.reverseAdaptor_;
		std::vector<std::string> adaptors{extractPars.forwardAdaptor_};
		if(!sameAdaptor){
			adaptors.emplace_back(extractPars.reverseAdaptor_);
		}
		std::unordered_map<uint64_t, std::set<std::tuple<uint32_t, uint32_t, uint32_t>>> seeds;
		for(uint32_t adaptorPos = 0; adaptorPos < adaptors.size(); ++adaptorPos){
			const auto & adaptor = adaptors[adaptorPos];
			adaptorSeqs_.emplace_back(std::array<seqInfo, 2>{
				seqInfo(njh::pasteAsStr("adaptor", adaptorPos, "-start"), adaptor),
				seqInfo(njh::pasteAsStr("adaptor", adaptorPos, "-end"), seqUtil::reverseComplement(adaptor, "DNA"))});
			//the adaptor with each barcode that can be in it filled in, the N's themselves never seed
			std::set<std::string> filled{adaptor};
			auto firstN = adaptor.find('N');
			auto lastN = adaptor.rfind('N');
			if(std::string::npos != firstN){
				for(const auto & bar : barcodes){
					VecStr bars{0 == adaptorPos ? bar.fwBar_.seq_ : bar.revBar_.seq_};
					if(sameAdaptor){
						bars.emplace_back(bar.revBar_.seq_);
					}
					for(const auto & barSeq : bars){
						if(barSeq.size() == lastN - firstN + 1){
							auto adaptorWithBar = adaptor;
							adaptorWithBar.replace(firstN, barSeq.size(), barSeq);
							filled.emplace(adaptorWithBar);
						}
					}
				}
			}
			//the fixed bases between the barcode and the insert
			if(std::string::npos != lastN && pars_.innerFlankKmerLen_ > 0 && adaptor.size() - lastN - 1 >= pars_.innerFlankKmerLen_){
				auto innerFlank = adaptor.substr(lastN + 1);
				for(uint32_t orientation = 0; orientation < 2; ++orientation){
					auto oriented = 0 == orientation ? innerFlank : seqUtil::reverseComplement(innerFlank, "DNA");
					for(uint32_t pos = 0; pos + pars_.innerFlankKmerLen_ <= oriented.size(); ++pos){
						innerFlankKmers_[oriented.substr(pos, pars_.innerFlankKmerLen_)].emplace_back(pos, oriented.size());
					}
				}
			}
			for(const auto & fill : filled){
				for(uint32_t orientation = 0; orientation < 2; ++orientation){
					auto oriented = 0 == orientation ? fill : seqUtil::reverseComplement(fill, "DNA");
					for(uint32_t pos = 0; pos + pars_.seedKmerLen_ <= oriented.size(); ++pos){
						uint64_t code = 0;
						if(encode(oriented, pos, code)){
							seeds[code].emplace(adaptorPos, orientation, pos);
						}
					}
				}
			}
		}
		for(const auto & seed : seeds){
			auto & hits = seeds_[seed.first];
			for(const auto & hit : seed.second){
				hits.emplace_back(SeedHit{std::get<0>(hit), std::get<1>(hit), std::get<2>(hit)});
			}
		}
	}

	/**@brief cut the read at the adaptors inside it
	 * @param read the read
	 * @param pieces filled with the pieces, or just the read when it isn't split
	 * @return whether the read was split
	 */
	bool split(const seqInfo & read, std::vector<seqInfo> & pieces){
		pieces.clear();
		if(!pars_.split_ || len(read) < 2 * minPieceLen_){
			pieces.emplace_back(read);
			return false;
		}
		++readsLongEnough_;
		auto boundaries = findBoundaries(read);
		if(boundaries.empty()){
			pieces.emplace_back(read);
			return false;
		}
		//pieces between the cuts, a piece too short to be output (e.g. what is between the two adaptors of a junction) is dropped
		std::vector<std::pair<uint32_t, uint32_t>> intervals;
		uint32_t start = 0;
		for(const auto & boundary : boundaries){
			if(boundary > start){
				intervals.emplace_back(start, boundary);
			}
			start = boundary;
		}
		intervals.emplace_back(start, len(read));
		//a piece has to hold an insert, not just adaptors, a molecule with its adaptor doubled up (common at the ends of PCR barcoded reads) is not two molecules
		std::vector<std::pair<uint32_t, uint32_t>> kept;
		for(const auto & interval : intervals){
			uint32_t adaptorBases = 0;
			for(const auto & span : adaptorSpans_){
				auto overlapStart = std::max(interval.first, span.first);
				auto overlapEnd = std::min(interval.second, span.second);
				if(overlapEnd > overlapStart){
					adaptorBases += overlapEnd - overlapStart;
				}
			}
			auto pieceLen = interval.second - interval.first;
			if(pieceLen >= adaptorBases + minPieceLen_){
				kept.emplace_back(interval);
			}else{
				++shortPiecesDropped_;
			}
		}
		//a molecule folded back on itself has the adaptors of its end then its start in the middle too, splitting it would count it twice
		for(uint32_t pos = 0; pos + 1 < kept.size();){
			if(isFoldback(read, kept[pos], kept[pos + 1])){
				++foldbackJunctions_;
				kept[pos].second = kept[pos + 1].second;
				kept.erase(kept.begin() + pos + 1);
			}else{
				++pos;
			}
		}
		if(kept.size() < 2){
			pieces.emplace_back(read);
			return false;
		}
		++readsSplit_;
		auto nameEnd = read.name_.find_first_of(" \t");
		if(std::string::npos == nameEnd){
			nameEnd = read.name_.size();
		}
		for(const auto & interval : kept){
			auto piece = read.getSubRead(interval.first, interval.second - interval.first);
			piece.name_ = read.name_;
			piece.name_.insert(nameEnd, njh::pasteAsStr("_split", interval.first, "-", interval.second));
			pieces.emplace_back(piece);
			++piecesFromSplits_;
		}
		return true;
	}

	/**@brief every adaptor (or bare adaptor tail, the fixed bases are what is checked) confirmed in the read, without cutting it
	 *
	 * For a read assigned on the barcode of one end, everything past the first of these belongs to another molecule whose own barcode wasn't read
	 *
	 * @param read the read
	 * @return the start and end of each adaptor in the read, sorted
	 */
	std::vector<std::pair<uint32_t, uint32_t>> findAdaptorSpans(const seqInfo & read){
		std::vector<std::pair<uint32_t, uint32_t>> spans;
		for(const auto & candidate : findCandidates(read)){
			uint32_t boundary = 0;
			std::pair<uint32_t, uint32_t> span;
			if(alignAt(read, candidate.windowStart_, candidate.windowEnd_ - candidate.windowStart_, candidate.orientation_, boundary, span)){
				spans.emplace_back(span);
			}
		}
		std::sort(spans.begin(), spans.end());
		return spans;
	}

	/**@brief where exact k-mers of an adaptor's inner flank (the fixed bases between the barcode and the insert) are in the read, a partial or bare adaptor tail
	 * @param read the read
	 * @return the span of the whole inner flank each hit implies, sorted
	 */
	std::vector<std::pair<uint32_t, uint32_t>> findInnerFlankSpans(const seqInfo & read) const {
		std::vector<std::pair<uint32_t, uint32_t>> spans;
		const auto kLen = pars_.innerFlankKmerLen_;
		if(0 == kLen || len(read) < kLen){
			return spans;
		}
		for(uint32_t pos = 0; pos + kLen <= len(read); ++pos){
			auto hits = innerFlankKmers_.find(read.seq_.substr(pos, kLen));
			if(innerFlankKmers_.end() == hits){
				continue;
			}
			for(const auto & hit : hits->second){
				int64_t start = static_cast<int64_t>(pos) - hit.first;
				spans.emplace_back(static_cast<uint32_t>(std::max<int64_t>(0, start)),
						static_cast<uint32_t>(std::min<int64_t>(len(read), start + hit.second)));
			}
		}
		std::sort(spans.begin(), spans.end());
		spans.erase(std::unique(spans.begin(), spans.end()), spans.end());
		return spans;
	}

	uint32_t readsSplit() const {
		return readsSplit_;
	}

	uint32_t piecesFromSplits() const {
		return piecesFromSplits_;
	}

	void writeCounts(const bfs::path & dir) const {
		OutputStream out(njh::files::make_path(dir, "innerAdaptorSplitCounts.tab.txt"));
		out << "condition\tcount" << std::endl;
		out << "readsLongEnoughToSplit" << "\t" << readsLongEnough_ << std::endl;
		out << "readsWithAdaptorSeeds" << "\t" << readsWithSeeds_ << std::endl;
		out << "windowsAligned" << "\t" << windowsAligned_ << std::endl;
		out << "innerAdaptorsStartingMolecule" << "\t" << adaptorsFound_[0] << std::endl;
		out << "innerAdaptorsEndingMolecule" << "\t" << adaptorsFound_[1] << std::endl;
		//only looked for in reads with an inner adaptor
		out << "adaptorsAtReadEndsOfReadsWithInnerAdaptors" << "\t" << adaptorsNearEnd_ << std::endl;
		out << "foldbackJunctionsLeftJoined" << "\t" << foldbackJunctions_ << std::endl;
		out << "piecesWithoutInsertDropped" << "\t" << shortPiecesDropped_ << std::endl;
		out << "readsSplit" << "\t" << readsSplit_ << std::endl;
		out << "piecesFromSplitReads" << "\t" << piecesFromSplits_ << std::endl;
	}

private:
	struct SeedHit {
		uint32_t adaptor_;
		//0 the adaptor as it reads at the start of a molecule, 1 its reverse complement as it reads at the end of one
		uint32_t orientation_;
		uint32_t offset_;
	};

	bool encode(const std::string & seq, uint32_t pos, uint64_t & code) const {
		code = 0;
		for(uint32_t kPos = pos; kPos < pos + pars_.seedKmerLen_; ++kPos){
			auto base = baseCode(seq[kPos]);
			if(base > 3){
				return false;
			}
			code = (code << 2) | base;
		}
		return true;
	}

	static uint64_t baseCode(char base){
		switch (base) {
			case 'A': case 'a': return 0;
			case 'C': case 'c': return 1;
			case 'G': case 'g': return 2;
			case 'T': case 't': return 3;
			default: return 4;
		}
	}

	//a place in the read the seeds say an adaptor may be
	struct Candidate {
		uint32_t orientation_;
		int64_t windowStart_;
		int64_t windowEnd_;
		//whether the seeds put this adaptor too close to an end of the read to be cut at
		bool nearReadEnd_;
	};

	/**@brief where the seeds say adaptors may be in the read
	 * @param read the read
	 * @return the candidates, seeds on about the same diagonal (allowing for indels) of the same orientation are one candidate adaptor
	 */
	std::vector<Candidate> findCandidates(const seqInfo & read){
		//orientation, diagonal (read position the adaptor would start at), read position of the seed
		std::vector<std::tuple<uint32_t, int64_t, uint32_t>> hits;
		const auto & seq = read.seq_;
		const uint64_t mask = (pars_.seedKmerLen_ >= 32) ? std::numeric_limits<uint64_t>::max() : ((uint64_t(1) << (2 * pars_.seedKmerLen_)) - 1);
		uint64_t code = 0;
		uint32_t validBases = 0;
		for(uint32_t pos = 0; pos < seq.size(); ++pos){
			auto base = baseCode(seq[pos]);
			if(base > 3){
				validBases = 0;
				code = 0;
				continue;
			}
			code = ((code << 2) | base) & mask;
			++validBases;
			if(validBases < pars_.seedKmerLen_){
				continue;
			}
			auto seedHits = seeds_.find(code);
			if(seeds_.end() == seedHits){
				continue;
			}
			uint32_t kmerStart = pos + 1 - pars_.seedKmerLen_;
			for(const auto & hit : seedHits->second){
				hits.emplace_back(hit.orientation_, static_cast<int64_t>(kmerStart) - hit.offset_, kmerStart);
			}
		}
		std::vector<Candidate> candidates;
		if(hits.size() < pars_.minSeedHits_){
			return candidates;
		}
		std::sort(hits.begin(), hits.end());
		//the forward and reverse adaptor share their barcodes so both are aligned and the better one kept
		const int64_t readLen = seq.size();
		for(size_t clusterStart = 0; clusterStart < hits.size();){
			size_t clusterEnd = clusterStart + 1;
			while(clusterEnd < hits.size() && std::get<0>(hits[clusterEnd]) == std::get<0>(hits[clusterStart])
					&& std::get<1>(hits[clusterEnd]) - std::get<1>(hits[clusterEnd - 1]) <= diagonalSlop_){
				++clusterEnd;
			}
			std::set<uint32_t> seedPositions;
			for(size_t hitPos = clusterStart; hitPos < clusterEnd; ++hitPos){
				seedPositions.emplace(std::get<2>(hits[hitPos]));
			}
			if(seedPositions.size() >= pars_.minSeedHits_){
				uint32_t orientation = std::get<0>(hits[clusterStart]);
				int64_t diagStart = std::get<1>(hits[clusterStart]);
				int64_t diagEnd = std::get<1>(hits[clusterEnd - 1]);
				int64_t windowStart = std::max<int64_t>(0, diagStart - windowPad_);
				int64_t windowEnd = std::min<int64_t>({readLen,
					diagEnd + static_cast<int64_t>(maxAdaptorLen_) + windowPad_,
					windowStart + static_cast<int64_t>(maxWindow_)});
				if(windowEnd > windowStart){
					int64_t approxCut = 0 == orientation ? diagStart : diagEnd + static_cast<int64_t>(maxAdaptorLen_);
					bool nearReadEnd = approxCut + nearEndSlop_ < static_cast<int64_t>(minPieceLen_)
							|| approxCut - nearEndSlop_ > readLen - static_cast<int64_t>(minPieceLen_);
					candidates.emplace_back(Candidate{orientation, windowStart, windowEnd, nearReadEnd});
				}
			}
			clusterStart = clusterEnd;
		}
		return candidates;
	}

	//positions in the read to cut at, an adaptor starting a molecule is cut before it and one ending a molecule after it
	std::vector<uint32_t> findBoundaries(const seqInfo & read){
		std::vector<uint32_t> boundaries;
		adaptorSpans_.clear();
		auto candidates = findCandidates(read);
		if(candidates.empty()){
			return boundaries;
		}
		++readsWithSeeds_;
		auto alignCandidate = [&](const Candidate & candidate){
			uint32_t boundary = 0;
			std::pair<uint32_t, uint32_t> span;
			if(alignAt(read, candidate.windowStart_, candidate.windowEnd_ - candidate.windowStart_, candidate.orientation_, boundary, span)){
				//every adaptor found counts against the piece it is in, including the ones at the ends of the read that aren't cut at
				adaptorSpans_.emplace_back(span);
				if(boundary >= minPieceLen_ && len(read) - boundary >= minPieceLen_){
					++adaptorsFound_[candidate.orientation_];
					boundaries.emplace_back(boundary);
				}else{
					++adaptorsNearEnd_;
				}
			}
		};
		//nearly every read has its own adaptors at its ends, those are only aligned (for how much of a piece is adaptor) once an inner adaptor is confirmed
		for(const auto & candidate : candidates){
			if(!candidate.nearReadEnd_){
				alignCandidate(candidate);
			}
		}
		if(!boundaries.empty()){
			for(const auto & candidate : candidates){
				if(candidate.nearReadEnd_){
					alignCandidate(candidate);
				}
			}
		}
		std::sort(boundaries.begin(), boundaries.end());
		boundaries.erase(std::unique(boundaries.begin(), boundaries.end()), boundaries.end());
		return boundaries;
	}

	bool alignAt(const seqInfo & read, uint32_t windowStart, uint32_t windowLen, uint32_t orientation, uint32_t & boundary, std::pair<uint32_t, uint32_t> & span){
		seqInfo window = read.getSubRead(windowStart, windowLen);
		bool found = false;
		AdaptorFlankStats best;
		for(const auto & adaptor : adaptorSeqs_){
			++windowsAligned_;
			const auto & oriented = adaptor[orientation];
			alignerObj_.alignRegGlobal(window, oriented);
			alignerObj_.profilePrimerAlignment(window, oriented);
			PrimerAlnInfo info(oriented.name_, alignerObj_.comp_, alignerObj_.alignObjectA_.seqBase_, alignerObj_.alignObjectB_.seqBase_);
			auto flank = getAdaptorFlankStats(info);
			if(flank.matches_ < pars_.flankMinMatches_ || flank.identity() < pars_.flankIdentityCutOff_ || flank.matches_ <= best.matches_){
				continue;
			}
			//the whole adaptor has to be in the window, otherwise it was only put there by running off the edge
			if('-' == info.seqAln_.seq_.front() || '-' == info.seqAln_.seq_.back()){
				continue;
			}
			best = flank;
			found = true;
			span.first = windowStart + getRealPosForAlnPos(info.seqAln_.seq_, info.primerAln_.seq_.find_first_not_of('-'));
			span.second = windowStart + getRealEndForAlnPos(info.seqAln_.seq_, info.primerAln_.seq_.find_last_not_of('-'));
			boundary = 0 == orientation ? span.first : span.second;
		}
		return found;
	}

	bool isFoldback(const seqInfo & read, const std::pair<uint32_t, uint32_t> & first, const std::pair<uint32_t, uint32_t> & second) const {
		const uint32_t kLen = 15;
		if(first.second - first.first < kLen || second.second - second.first < kLen){
			return false;
		}
		std::unordered_set<std::string> firstKmers;
		for(uint32_t pos = first.first; pos + kLen <= first.second; ++pos){
			firstKmers.emplace(read.seq_.substr(pos, kLen));
		}
		auto secondRevComp = seqUtil::reverseComplement(read.seq_.substr(second.first, second.second - second.first), "DNA");
		uint32_t shared = 0;
		uint32_t secondKmers = 0;
		for(uint32_t pos = 0; pos + kLen <= secondRevComp.size(); ++pos){
			++secondKmers;
			if(firstKmers.count(secondRevComp.substr(pos, kLen))){
				++shared;
			}
		}
		return shared >= pars_.foldbackKmerFrac_ * std::min<uint32_t>(firstKmers.size(), secondKmers);
	}

	InnerAdaptorSplitPars pars_;
	uint32_t minPieceLen_;
	static constexpr int64_t windowPad_ = 20;
	static constexpr int64_t diagonalSlop_ = 10;
	//how far off the cut the seeds point at can be from where aligning puts it
	static constexpr int64_t nearEndSlop_ = 30;
	uint64_t maxAdaptorLen_;
	uint64_t maxWindow_;
	aligner alignerObj_;
	//per adaptor, the start then end orientation
	std::vector<std::array<seqInfo, 2>> adaptorSeqs_;
	std::unordered_map<uint64_t, std::vector<SeedHit>> seeds_;
	//where the adaptors confirmed in the current read are
	std::vector<std::pair<uint32_t, uint32_t>> adaptorSpans_;
	//inner flank k-mer -> its position in the oriented flank and the flank's length
	std::unordered_map<std::string, std::vector<std::pair<uint32_t, uint32_t>>> innerFlankKmers_;

	uint32_t readsLongEnough_ = 0;
	uint32_t readsWithSeeds_ = 0;
	uint32_t windowsAligned_ = 0;
	std::array<uint32_t, 2> adaptorsFound_{0, 0};
	uint32_t adaptorsNearEnd_ = 0;
	uint32_t foldbackJunctions_ = 0;
	uint32_t shortPiecesDropped_ = 0;
	uint32_t readsSplit_ = 0;
	uint32_t piecesFromSplits_ = 0;
};

int seqUtilsExtractRunner::extractByIlluminaAaptors(const njh::progutils::CmdArgs & inputCommands) {
	AdaptorExtractPars pars;
	bfs::path barcodeSampleSheet = "";
	//when not 0, a barcode within this many matches of the best match is just as likely, see BarcodeEndCaller
	uint32_t barcodeMatchMargin = 0;
	bool development = false;
	InnerAdaptorSplitPars splitPars;

	seqUtilsExtractSetUp setUp(inputCommands);
	setAdaptorExtractOptions(setUp, pars);
	setUp.setOption(barcodeMatchMargin, "--barcodeMatchMargin",
			"A barcode needs this many more matching bases than any other to be called by itself, 0 means just take the best (ties broken by alignment score)");
	setUp.setOption(barcodeSampleSheet, "--illuminaBarcodeSampleSheet,--barcodeSampleSheet",
			"Barcode Sample Sheet, columns sample_name, fw, rev (barcode sequences)", true);
	setInnerAdaptorSplitOptions(setUp, splitPars);
	setUp.setOption(development, "--development",
			"Also write the per read barcode info files (info.tab.txt, infoBest.tab.txt, infoPass.tab.txt)");
	uint32_t numThreads = 1;
	setUp.setOption(numThreads, "--numThreads", "Number of threads, reads are written out in a different order with more than one");

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
	auto barcodes = readBarcodeSampleSheet(barcodeSampleSheet, nullptr);
	//make outputs
	auto fastqDir = njh::files::make_path(setUp.pars_.directoryName_, "fastq");
	njh::files::makeDir(njh::files::MkdirPar{fastqDir});
	MultiSeqOutCache<seqInfo> outputs;
	for(const auto & bar : barcodes){
		outputs.addReader(bar.sampleName_, SeqIOOptions::genFastqOutGz(njh::files::make_path(fastqDir, bar.sampleName_)));
	}
	outputs.addReader("Undetermined", SeqIOOptions::genFastqOutGz(njh::files::make_path(fastqDir, "Undetermined")));
	//reads whose two ends were each confidently called as different samples, kept apart from Undetermined (an end that couldn't be called)
	//and out of fastq/, since they measure between sample joining or barcode misreads rather than failed barcode calls
	auto mismatchedSamplesDir = njh::files::make_path(setUp.pars_.directoryName_, "mismatchedSamplesFastq");
	njh::files::makeDir(njh::files::MkdirPar{mismatchedSamplesDir});
	outputs.addReader(mismatchedSamplesName, SeqIOOptions::genFastqOutGz(njh::files::make_path(mismatchedSamplesDir, mismatchedSamplesName)));

	//the per thread finders and splitters are added into these for the reports
	AdaptorPairFinder finderTotals(pars, setUp.pars_.debug_);
	InnerAdaptorSplitter splitterTotals(pars, barcodes, splitPars);

	SeqInput seqReader(setUp.pars_.ioOptions_);
	seqReader.openIn();

	AdaptorInfoFiles infoFiles(setUp.pars_.directoryName_, development);
	DemuxCounts counts;
	std::mutex outMut;
	std::mutex countsMut;

	std::function<void()> demuxReads = [&](){
		//aligners and caches aren't shared between threads
		BarcodeEndCaller caller(barcodes, pars.barIdentityCutOff_, barcodeMatchMargin, true);
		AdaptorPairFinder finder(pars, setUp.pars_.debug_);
		InnerAdaptorSplitter splitter(pars, barcodes, splitPars);
		DemuxCounts local;
		std::ostringstream infoBuf;
		std::ostream noInfo(nullptr);
		std::ostream & infoOut = infoFiles.on() ? static_cast<std::ostream &>(infoBuf) : noInfo;
		auto writeOut = [&outputs, &outMut](const std::string & sample, const seqInfo & read){
			std::lock_guard<std::mutex> lock(outMut);
			outputs.add(sample, read);
		};
		seqInfo inRead;
		std::vector<seqInfo> pieces;
		while(seqReader.readNextReadLock(inRead)){
			//a read holding several library molecules is demultiplexed a molecule at a time
			splitter.split(inRead, pieces);
			for(auto & seq : pieces){
				auto res = finder.process(seq, infoOut);
				if(AdaptorReadOutcome::pass != res.outcome_){
					continue;
				}
				const auto & processRes = res.processRes_;
				const auto & forwardBar = res.forwardBar_;
				const auto & reverseBar = res.reverseBar_;
				++local.fullPass;
				auto call = callBarcodePair(caller, forwardBar, reverseBar, true, false, 0);
				++local.barcodePairCounts[std::make_tuple(call.fwdCall_, call.revCall_, call.pairStatus_)];
				++local.pairStatusCounts[call.pairStatus_];
				MetaDataInName seqMeta;
				seqMeta.addMeta("direction", res.direction_);
				seqMeta.addMeta("fprimerStart", processRes.fprimerStart_);
				seqMeta.addMeta("forwardBar", forwardBar);
				seqMeta.addMeta("reverseBar", reverseBar);
				seqMeta.addMeta("fullBar", njh::pasteAsStr(forwardBar, "+", reverseBar));
				seqMeta.addMeta("rprimerStart", processRes.rprimerStart_);
				seqMeta.addMeta("sample", call.determinedSample_);
				if("matched" != call.pairStatus_){
					seqMeta.addMeta("fwdSample", call.fwdCall_);
					seqMeta.addMeta("revSample", call.revCall_);
					seqMeta.addMeta("pairStatus", call.pairStatus_);
				}
				seq.name_ += seqMeta.createMetaName();
				writeOut(call.determinedSample_, seq);
				++local.sampleCounts[call.determinedSample_];
				if("matched" == call.pairStatus_) {
					++local.sampleDetermined;
				}
				infoFiles.writePass(seq.name_, res);
			}
			if(infoFiles.on()){
				infoFiles.writeInfo(infoBuf.str());
				infoBuf.str("");
			}
		}
		std::lock_guard<std::mutex> lock(countsMut);
		counts.add(local);
		finderTotals.addCounts(finder);
		splitterTotals.addCounts(splitter);
	};
	njh::concurrent::runVoidFunctionThreaded(demuxReads, numThreads);
	auto & finder = finderTotals;
	auto & splitter = splitterTotals;
	auto & sampleCounts = counts.sampleCounts;
	auto & barcodePairCounts = counts.barcodePairCounts;
	auto & pairStatusCounts = counts.pairStatusCounts;
	const uint32_t fullPass = counts.fullPass;
	const uint32_t sampleDetermined = counts.sampleDetermined;
	finder.writeBestDirectionCounts(setUp.pars_.directoryName_);
	finder.writeFailedCounts(setUp.pars_.directoryName_);
	{
		std::vector<std::pair<std::string, uint32_t>> afterPasses{
			{"fullPass", fullPass},
			{"Undetermined", sampleCounts["Undetermined"]},
			{"sampleDetermined", sampleDetermined},
			{mismatchedSamplesName, pairStatusCounts["mismatched"]},
			{"oneEndUndetermined", pairStatusCounts["oneEndUndetermined"]},
			{"bothEndsUndetermined", pairStatusCounts["bothEndsUndetermined"]},
			{"ambiguousSamplePair", pairStatusCounts["ambiguousPair"]}};
		//the fractions in conditionCount are of the reads plus the extra pieces split reads make
		if(splitPars.split_){
			afterPasses.emplace_back("readsSplitAtInnerAdaptors", splitter.readsSplit());
			afterPasses.emplace_back("piecesFromSplitReads", splitter.piecesFromSplits());
			splitter.writeCounts(setUp.pars_.directoryName_);
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
	writeBarcodePairCounts(setUp.pars_.directoryName_, barcodePairCounts, fullPass);

	return 0;
}

//the adaptor found at one end of a nanopore read looked at by itself, without needing the other end
struct NanoporeEndHit {
	//whether this is the forward adaptor rather than the reverse adaptor
	bool forwardAdaptor_ = true;
	//whether the whole inner side of the adaptor is within the search window, otherwise it was only placed there by running off of the window
	bool inWindow_ = false;
	AdaptorFlankStats flank_;
	//the flanks say the adaptor is really here
	bool adaptorDetected_ = false;
	std::string barcode_;
	BarEndCall call_;
	//the call is good enough to assign the read by this end alone
	bool strictCall_ = false;
	//position of the adaptor in the read
	uint32_t adaptorStart_ = 0;
	uint32_t adaptorEnd_ = 0;
};

//settings for assigning a read on the barcode of one end
struct OneEndRescuePars {
	bool bothEndsOnly_ = false;
	//ends confidently calling different samples are more often a real swap than a misread, so a stricter check doesn't say which end is right
	bool fromMismatched_ = false;
	//the other end calling a different sample at the looser cut off is often a real swap too, keeping these doubled the reads put in the wrong panel in PCR barcoding run 3
	bool keepConflicting_ = false;
	//by default a read assigned on one end is cut back to the first adaptor inside it
	bool keepPastInnerAdaptor_ = false;
	//below this the reads called as unused kit barcodes start climbing (21 of 24 bases)
	double barcodeMatchFrac_ = .875;
	uint32_t barcodeMargin_ = 3;
	double flankIdentityCutOff_ = .75;
	uint32_t flankMinMatches_ = 12;
};

/**@brief get the best adaptor for one end of a read
 * @param preferred the adaptor that a read in the forward direction has at this end, wins ties (native kits have the same adaptor on both)
 * @param other the other adaptor at this end
 * @param atFront whether this is the front of the read
 * @param startOfWindow the position of the search window in the read
 * @param caller the barcode caller
 * @param pars rescue settings
 */
NanoporeEndHit getNanoporeEndHit(const PrimerAlnInfo & preferred, const PrimerAlnInfo & other, bool atFront, uint32_t startOfWindow,
		BarcodeEndCaller & caller, const OneEndRescuePars & pars){
	auto preferredFlank = getAdaptorFlankStats(preferred);
	auto otherFlank = getAdaptorFlankStats(other);
	bool usePreferred = preferredFlank.matches_ > otherFlank.matches_ ||
			(preferredFlank.matches_ == otherFlank.matches_ && preferred.comp_.distances_.eventBasedIdentity_ >= other.comp_.distances_.eventBasedIdentity_);
	const auto & best = usePreferred ? preferred : other;
	NanoporeEndHit ret;
	//the preferred adaptor is the forward one at the front and the reverse one at the back
	ret.forwardAdaptor_ = usePreferred == atFront;
	ret.flank_ = usePreferred ? preferredFlank : otherFlank;
	ret.inWindow_ = atFront ? '-' != best.seqAln_.seq_.back() : '-' != best.seqAln_.seq_.front();
	if(!ret.inWindow_){
		return ret;
	}
	ret.adaptorDetected_ = ret.flank_.matches_ >= pars.flankMinMatches_ && ret.flank_.identity() >= pars.flankIdentityCutOff_;
	ret.barcode_ = atFront ? best.barcode_ : seqUtil::reverseComplement(best.barcode_, "DNA");
	if(caller.callable(ret.barcode_)){
		ret.call_ = caller.call(ret.barcode_, ret.forwardAdaptor_);
	}
	ret.strictCall_ = !ret.call_.undetermined() && 1 == ret.call_.samples_.size()
			&& ret.call_.matches_ >= pars.barcodeMatchFrac_ * ret.call_.barcodeLen_
			&& ret.call_.matches_ >= ret.call_.nextBestMatches_ + pars.barcodeMargin_;
	ret.adaptorStart_ = startOfWindow + getRealPosForAlnPos(best.seqAln_.seq_, best.primerAln_.seq_.find_first_not_of('-'));
	ret.adaptorEnd_ = startOfWindow + getRealEndForAlnPos(best.seqAln_.seq_, best.primerAln_.seq_.find_last_not_of('-'));
	return ret;
}

std::string getOutcomeName(AdaptorReadOutcome outcome){
	switch (outcome) {
		case AdaptorReadOutcome::smallSeq:
			return "smallSeq";
		case AdaptorReadOutcome::largeSeq:
			return "largeSeq";
		case AdaptorReadOutcome::mismatchDirections:
			return "mismatchDirections";
		case AdaptorReadOutcome::failedIdentity:
			return "failedIdentity";
		case AdaptorReadOutcome::failedPrimerAln:
			return "failedPrimerAln";
		case AdaptorReadOutcome::failedOverlapCheck:
			return "failedOverlapCheck";
		case AdaptorReadOutcome::pass:
			return "pass";
	}
	return "unknown";
}

int seqUtilsExtractRunner::extractByNanoporeBarcodes(const njh::progutils::CmdArgs & inputCommands) {
	AdaptorExtractPars pars;
	bfs::path barcodeSampleSheet = "";
	//a barcode within this many matches of the best match is just as likely, see BarcodeEndCaller
	uint32_t barcodeMatchMargin = 3;
	std::string kitName;
	bool development = false;
	bool keepMismatchDirectionsWithMatchingBarcodes = false;
	double mismatchConfirmIdentityCutOff = .80;
	OneEndRescuePars rescuePars;

	seqUtilsExtractSetUp setUp(inputCommands);
	setUp.description_ = "Demultiplex nanopore reads on the ONT barcodes at the ends of the read. Reads with the same barcode at both ends are assigned first, then reads that only have one end's barcode are assigned if it passes a stricter check, since the ends of nanopore reads are often cut short";
	pars.searchWithin_ = 175;
	pars.minOutputLen_ = 100;
	pars.extractPars_.minTruncatedBackBarcodeBases_ = 14;
	//with the margin a lower cut off lets primer/insert sequence match several barcodes and leaves the end undecided
	pars.barIdentityCutOff_ = .75;
	setUp.setOption(keepMismatchDirectionsWithMatchingBarcodes, "--keepMismatchDirectionsWithMatchingBarcodes",
			"Reads with the same adaptor at both ends (mismatchDirections, only possible with a kit like PCR barcoding that has different adaptors, native reads can't be called this way) go to the regular sample file, trimmed like the other reads, when both barcodes call that sample and the adaptors pass the usual checks, otherwise all of them are kept untrimmed in mismatchDirectionsFastq/");
	setUp.setOption(mismatchConfirmIdentityCutOff, "--mismatchConfirmIdentityCutOff",
			"Reads whose ends call different samples only count as mismatched (mismatchedSamples) when both called barcodes reach this identity, otherwise they go to Undetermined as unconfirmed; assignment to samples only uses --barIdentityCutOff");
	setUp.setOption(kitName, "--kit", "ONT barcoding kit, sets the adaptors and the barcode names for the sample sheet, e.g. SQK-NBD114-96 or EXP-PBC096", true);
	//empty when just asking for the help message
	OntBarcodeKit kit;
	if(!kitName.empty()){
		kit = getOntBarcodeKit(kitName);
		pars.forwardAdaptor_ = kit.forwardAdaptor_;
		pars.reverseAdaptor_ = kit.reverseAdaptor_;
	}
	setAdaptorExtractOptions(setUp, pars);
	setUp.setOption(barcodeMatchMargin, "--barcodeMatchMargin",
			"A barcode needs this many more matching bases than any other to be called by itself, 0 means just take the best (ties broken by alignment score)");
	setUp.setOption(barcodeSampleSheet, "--illuminaBarcodeSampleSheet,--barcodeSampleSheet",
			"Barcode Sample Sheet, columns sample_name, fw, rev (barcode sequences); can also be sample_name and barcode (e.g. barcode01) or a dorado name key (barcode name then sample, no header)", true);
	setUp.setOption(rescuePars.bothEndsOnly_, "--bothEndsOnly", "Only assign reads with the same barcode called at both ends, no assigning on one end's barcode");
	setUp.setOption(rescuePars.fromMismatched_, "--oneEndFromMismatched", "Also assign on one end's barcode when the two ends confidently called different samples, otherwise those stay in mismatchedSamplesFastq/");
	setUp.setOption(rescuePars.keepConflicting_, "--oneEndKeepConflicting",
			"Also assign on one end's barcode when the other end called a different sample below the one end check, by default these are left out since many of them are real swaps");
	setUp.setOption(rescuePars.keepPastInnerAdaptor_, "--oneEndKeepPastInnerAdaptor",
			"When assigning on one end's barcode, keep the read past an adaptor (or bare adaptor tail) found inside it; by default the read is cut back to the first one, since what is past it is another molecule whose barcode wasn't read");
	setUp.setOption(rescuePars.barcodeMatchFrac_, "--oneEndBarcodeMatchFrac",
			"To assign a read on one end's barcode, the fraction of the barcode's bases that have to match, a barcode cut short by the end of the read counts against this");
	setUp.setOption(rescuePars.barcodeMargin_, "--oneEndBarcodeMargin",
			"To assign a read on one end's barcode, how many more matching bases it needs than any other barcode");
	setUp.setOption(rescuePars.flankIdentityCutOff_, "--adaptorFlankIdentityCutOff",
			"When assigning on one end's barcode, the other end is only trimmed when the adaptor's fixed bases (not the barcode) match at this identity");
	setUp.setOption(rescuePars.flankMinMatches_, "--adaptorFlankMinMatches",
			"When assigning on one end's barcode, the other end is only trimmed when at least this many of the adaptor's fixed bases match");
	setUp.setOption(development, "--development",
			"Also write the per read barcode info files (info.tab.txt, infoBest.tab.txt, infoPass.tab.txt, oneEndInfo.tab.txt)");
	InnerAdaptorSplitPars splitPars;
	setInnerAdaptorSplitOptions(setUp, splitPars);
	uint32_t numThreads = 1;
	setUp.setOption(numThreads, "--numThreads", "Number of threads, reads are written out in a different order with more than one");

	setUp.processVerbose();
	setUp.processDebug();
	setUp.processReadInNames(true);
	setUp.processDirectoryOutputName(true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);

	//read in barcodes
	//check for same barcodes
	auto barcodes = readBarcodeSampleSheet(barcodeSampleSheet, &kit);
	//kit barcodes no sample used still compete for reads so that a read from one of them isn't given to the closest sample
	std::unordered_set<std::string> unusedBarcodeNames;
	{
		const auto & kitSeqs = getOntBarcodeSeqs().at(kit.barcodePrefix_);
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
	//make outputs
	auto fastqDir = njh::files::make_path(setUp.pars_.directoryName_, "fastq");
	njh::files::makeDir(njh::files::MkdirPar{fastqDir});
	//kept out of fastq/ so a pipeline reading every file there doesn't take the unused kit barcodes as samples
	auto unusedBarcodesDir = njh::files::make_path(setUp.pars_.directoryName_, "unusedBarcodes");
	if(!unusedBarcodeNames.empty()){
		njh::files::makeDir(njh::files::MkdirPar{unusedBarcodesDir});
	}
	MultiSeqOutCache<seqInfo> outputs;
	for(const auto & bar : barcodes){
		const auto & outDir = unusedBarcodeNames.count(bar.sampleName_) ? unusedBarcodesDir : fastqDir;
		outputs.addReader(bar.sampleName_, SeqIOOptions::genFastqOutGz(njh::files::make_path(outDir, bar.sampleName_)));
	}
	outputs.addReader("Undetermined", SeqIOOptions::genFastqOutGz(njh::files::make_path(fastqDir, "Undetermined")));
	//reads whose two ends were each confidently called as different samples, kept apart from Undetermined (an end that couldn't be called)
	//and out of fastq/, since they measure between sample joining or barcode misreads rather than failed barcode calls
	auto mismatchedSamplesDir = njh::files::make_path(setUp.pars_.directoryName_, "mismatchedSamplesFastq");
	njh::files::makeDir(njh::files::MkdirPar{mismatchedSamplesDir});
	outputs.addReader(mismatchedSamplesName, SeqIOOptions::genFastqOutGz(njh::files::make_path(mismatchedSamplesDir, mismatchedSamplesName)));
	//reads with the same adaptor at both ends, left untrimmed since this is where the foldback reads of a kit with two different adaptors end up
	const std::string mismatchDirectionsName = "mismatchDirections";
	auto mismatchDirectionsDir = njh::files::make_path(setUp.pars_.directoryName_, "mismatchDirectionsFastq");
	njh::files::makeDir(njh::files::MkdirPar{mismatchDirectionsDir});
	outputs.addReader(mismatchDirectionsName, SeqIOOptions::genFastqOutGz(njh::files::make_path(mismatchDirectionsDir, mismatchDirectionsName)));

	//the per thread finders and splitters are added into these for the reports
	AdaptorPairFinder finderTotals(pars, setUp.pars_.debug_);
	InnerAdaptorSplitter splitterTotals(pars, barcodes, splitPars);

	SeqInput seqReader(setUp.pars_.ioOptions_);
	seqReader.openIn();

	AdaptorInfoFiles infoFiles(setUp.pars_.directoryName_, development);
	std::unique_ptr<OutputStream> oneEndInfoOut;
	if(development){
		oneEndInfoOut = std::make_unique<OutputStream>(njh::files::make_path(setUp.pars_.directoryName_, "oneEndInfo.tab.txt"));
		*oneEndInfoOut << "seqName\tseqLen\tfirstPass"
				<< "\tfrontInWindow\tfrontAdaptor\tfrontFlankMatches\tfrontFlankColumns\tfrontBarLen\tfrontCall\tfrontMatches\tfrontNextBest\tfrontBarcodeLen\tfrontIdentity"
				<< "\tbackInWindow\tbackAdaptor\tbackFlankMatches\tbackFlankColumns\tbackBarLen\tbackCall\tbackMatches\tbackNextBest\tbackBarcodeLen\tbackIdentity"
				<< "\tassignedOn\tsample" << std::endl;
	}

	DemuxCounts counts;
	//an inner adaptor counts as past the assigning one when it starts no more than this before the assigning adaptor's end
	const uint32_t innerAdaptorSlop = 10;
	std::mutex outMut;
	std::mutex infoMut;
	std::mutex countsMut;

	std::function<void()> demuxReads = [&](){
		//aligners and caches aren't shared between threads
		BarcodeEndCaller caller(barcodes, pars.barIdentityCutOff_, barcodeMatchMargin, false);
		AdaptorPairFinder finder(pars, setUp.pars_.debug_);
		//seeded with the unused kit barcodes too, since joined molecules can come from any of them
		InnerAdaptorSplitter splitter(pars, barcodes, splitPars);

		DemuxCounts local;
		auto & sampleCounts = local.sampleCounts;
		auto & oneEndSampleCounts = local.oneEndSampleCounts;
		auto & barcodePairCounts = local.barcodePairCounts;
		auto & pairStatusCounts = local.pairStatusCounts;
		auto & mismatchDirectionPairCounts = local.mismatchDirectionPairCounts;
		auto & mismatchDirectionStatusCounts = local.mismatchDirectionStatusCounts;
		auto & oneEndFromCounts = local.oneEndFromCounts;
		auto & oneEndUsedCounts = local.oneEndUsedCounts;
		auto & mismatchDirectionsToSamples = local.mismatchDirectionsToSamples;
		auto & fullPass = local.fullPass;
		auto & sampleDeterminedBothEnds = local.sampleDeterminedBothEnds;
		auto & sampleDeterminedOneEnd = local.sampleDeterminedOneEnd;
		auto & oneEndToUnusedBarcode = local.oneEndToUnusedBarcode;
		auto & oneEndConflicting = local.oneEndConflicting;
		auto & oneEndCutAtInnerAdaptor = local.oneEndCutAtInnerAdaptor;
		auto & oneEndNoInsertBeforeInnerAdaptor = local.oneEndNoInsertBeforeInnerAdaptor;
		std::ostringstream infoBuf;
		std::ostream noInfo(nullptr);
		std::ostream & infoOut = infoFiles.on() ? static_cast<std::ostream &>(infoBuf) : noInfo;
		auto writeOut = [&outputs, &outMut](const std::string & sample, const seqInfo & read){
			std::lock_guard<std::mutex> lock(outMut);
			outputs.add(sample, read);
		};
		seqInfo inRead;
		std::vector<seqInfo> pieces;

		/**@brief try to assign a read that didn't get the same barcode at both ends on just one end
		 * @param wholeRead the read before any trimming
		 * @param res the result from process() for this read
		 * @param firstPass what the first pass made of this read
		 * @return whether the read was written out to a sample
		 */
		auto assignOnOneEnd = [&](const seqInfo & wholeRead, const AdaptorReadRes & res, const std::string & firstPass){
			uint32_t startOfBackSeq = wholeRead.seq_.size() - std::min<uint32_t>(pars.searchWithin_, len(wholeRead));
			auto front = getNanoporeEndHit(*res.fwd5_3_, *res.rev5_3_, true, 0, caller, rescuePars);
			auto back = getNanoporeEndHit(*res.rev3_5_, *res.fwd3_5_, false, startOfBackSeq, caller, rescuePars);
			std::string assignedOn = "none";
			std::string sample;
			//same adaptor really at both ends is a foldback, those are kept out of the samples
			bool foldback = front.adaptorDetected_ && back.adaptorDetected_ && front.forwardAdaptor_ == back.forwardAdaptor_;
			if(!foldback && !rescuePars.bothEndsOnly_){
				if(front.strictCall_ && back.strictCall_){
					if(front.call_.samples_.front() == back.call_.samples_.front()){
						assignedOn = "both";
						sample = front.call_.samples_.front();
					}
				}else if(front.strictCall_){
					assignedOn = "front";
					sample = front.call_.samples_.front();
				}else if(back.strictCall_){
					assignedOn = "back";
					sample = back.call_.samples_.front();
				}
			}
			bool trimFront = false;
			bool trimBack = false;
			if("none" != assignedOn){
				trimFront = "back" != assignedOn || front.adaptorDetected_;
				trimBack = "front" != assignedOn || back.adaptorDetected_;
				//the adaptors can only overlap when one of them isn't real, the end the read was assigned on is real
				if(trimFront && trimBack && front.adaptorEnd_ >= back.adaptorStart_){
					("back" == assignedOn ? trimFront : trimBack) = false;
				}
			}
			uint32_t keepStart = 0;
			uint32_t keepEnd = wholeRead.seq_.size();
			if(trimFront){
				keepStart = pars.extractPars_.keepPrimer_ ? front.adaptorStart_ : front.adaptorEnd_;
			}
			if(trimBack){
				keepEnd = pars.extractPars_.keepPrimer_ ? back.adaptorEnd_ : back.adaptorStart_;
			}
			if("none" != assignedOn && (keepEnd <= keepStart || keepEnd - keepStart < pars.minOutputLen_)){
				assignedOn = "none";
			}
			//checked after the length so the count is only reads that would have been written out
			if(("front" == assignedOn || "back" == assignedOn) && !rescuePars.keepConflicting_){
				const auto & otherCall = "front" == assignedOn ? back.call_ : front.call_;
				if(!otherCall.undetermined() && !njh::in(sample, otherCall.samples_)){
					assignedOn = "none";
					++oneEndConflicting;
				}
			}
			//an adaptor (or a bare adaptor tail) inside a read assigned on one end means it is joined to another molecule whose own barcode wasn't read,
			//so only what is between the assigning barcode and the first of these is kept (run 3: 122 wrong panel reads per 1000 past a bare tail, 5 with nothing in between,
			//while reads assigned on both ends stay clean with one)
			if(("front" == assignedOn || "back" == assignedOn) && !rescuePars.keepPastInnerAdaptor_){
				auto spans = splitter.findAdaptorSpans(wholeRead);
				//partial or bare tails the adaptor check misses
				auto tailSpans = splitter.findInnerFlankSpans(wholeRead);
				spans.insert(spans.end(), tailSpans.begin(), tailSpans.end());
				std::sort(spans.begin(), spans.end());
				bool cut = false;
				if("front" == assignedOn){
					for(const auto & span : spans){
						//mostly past the assigning adaptor, the assigning adaptor itself is found again here too
						if(span.first + innerAdaptorSlop >= front.adaptorEnd_ && span.second > front.adaptorEnd_ + innerAdaptorSlop && span.first + innerAdaptorSlop < keepEnd){
							keepEnd = span.first;
							trimBack = false;
							cut = true;
							break;
						}
					}
				}else{
					for(auto span = spans.rbegin(); span != spans.rend(); ++span){
						if(span->second <= back.adaptorStart_ + innerAdaptorSlop && span->first + innerAdaptorSlop < back.adaptorStart_ && span->second > keepStart + innerAdaptorSlop){
							keepStart = span->second;
							trimFront = false;
							cut = true;
							break;
						}
					}
				}
				if(cut){
					++oneEndCutAtInnerAdaptor;
					if(keepEnd <= keepStart || keepEnd - keepStart < pars.minOutputLen_){
						assignedOn = "none";
						++oneEndNoInsertBeforeInnerAdaptor;
					}
				}
			}
			if(oneEndInfoOut){
				std::ostringstream line;
				auto writeEnd = [&line](const NanoporeEndHit & hit){
					line << "\t" << njh::boolToStr(hit.inWindow_)
							<< "\t" << (hit.forwardAdaptor_ ? "fwd" : "rev")
							<< "\t" << hit.flank_.matches_
							<< "\t" << hit.flank_.columns_
							<< "\t" << hit.barcode_.size()
							<< "\t" << njh::conToStr(hit.call_.samples_, ",")
							<< "\t" << hit.call_.matches_
							<< "\t" << hit.call_.nextBestMatches_
							<< "\t" << hit.call_.barcodeLen_
							<< "\t" << hit.call_.identity_;
				};
				line << wholeRead.name_ << "\t" << len(wholeRead) << "\t" << firstPass;
				writeEnd(front);
				writeEnd(back);
				line << "\t" << assignedOn << "\t" << ("none" == assignedOn ? "NA" : sample) << "\n";
				std::lock_guard<std::mutex> lock(infoMut);
				*oneEndInfoOut << line.str();
			}
			if("none" == assignedOn){
				return false;
			}
			seqInfo trimmed = wholeRead.getSubRead(keepStart, keepEnd - keepStart);
			if(pars.extractPars_.keepPrimer_ && !pars.extractPars_.primerUpper_){
				if(trimFront){
					changeSubStrToLower(trimmed.seq_, 0, front.adaptorEnd_ - front.adaptorStart_);
				}
				if(trimBack){
					changeSubStrToLower(trimmed.seq_, back.adaptorStart_ - keepStart, back.adaptorEnd_ - back.adaptorStart_);
				}
			}
			const auto & usedEnd = "back" == assignedOn ? back : front;
			bool reverse = "back" == assignedOn ? usedEnd.forwardAdaptor_ : !usedEnd.forwardAdaptor_;
			if(reverse){
				trimmed.reverseComplementRead(false, true);
			}
			MetaDataInName seqMeta;
			seqMeta.addMeta("direction", reverse ? "rev" : "fwd");
			seqMeta.addMeta("assignedOn", assignedOn);
			seqMeta.addMeta("firstPass", firstPass);
			seqMeta.addMeta("frontBar", front.barcode_);
			seqMeta.addMeta("backBar", back.barcode_);
			seqMeta.addMeta("trimmedFront", njh::boolToStr(trimFront));
			seqMeta.addMeta("trimmedBack", njh::boolToStr(trimBack));
			seqMeta.addMeta("sample", sample);
			//how many of the read's ends had the sample's barcode, for clustering to count (SeekDeep processClusters --minBothEndsBarcodeReads)
			seqMeta.addMeta("barcodeEnds", "both" == assignedOn ? 2 : 1);
			trimmed.name_ += seqMeta.createMetaName();
			writeOut(sample, trimmed);
			++sampleCounts[sample];
			++oneEndSampleCounts[sample];
			++oneEndFromCounts[firstPass];
			++oneEndUsedCounts[assignedOn];
			if(unusedBarcodeNames.count(sample)){
				++oneEndToUnusedBarcode;
			}else{
				++sampleDeterminedOneEnd;
			}
			return true;
		};

		while(seqReader.readNextReadLock(inRead)){
			//a read holding several library molecules is demultiplexed a molecule at a time, so each insert goes with the barcodes closest to it
			splitter.split(inRead, pieces);
			for(auto & seq : pieces){
				//process() trims the reads that pass, assigning on one end decides itself what to trim
				seqInfo wholeRead = seq;
				auto res = finder.process(seq, infoOut);
				if(AdaptorReadOutcome::smallSeq == res.outcome_ || AdaptorReadOutcome::largeSeq == res.outcome_){
					continue;
				}
				if(AdaptorReadOutcome::mismatchDirections == res.outcome_){
					//both barcodes come from the one adaptor, so they're both compared against the barcodes for that adaptor
					auto call = callBarcodePair(caller, res.forwardBar_, res.reverseBar_, res.sameAdaptorIsForward_, res.sameAdaptorIsForward_, mismatchConfirmIdentityCutOff);
					++mismatchDirectionPairCounts[std::make_tuple(res.direction_, call.fwdCall_, call.revCall_, call.pairStatus_)];
					++mismatchDirectionStatusCounts[call.pairStatus_];
					if("mismatched" == call.pairStatus_ && (unusedBarcodeNames.count(call.fwdCall_) || unusedBarcodeNames.count(call.revCall_))){
						++mismatchDirectionStatusCounts["mismatchedWithUnusedBarcode"];
					}
					if("matched" == call.pairStatus_ && unusedBarcodeNames.count(call.determinedSample_)){
						++mismatchDirectionStatusCounts["matchedUnusedBarcode"];
					}
					//most of these only look like the same adaptor at both ends because one end has no adaptor and either adaptor fits it about as well
					if("matched" != call.pairStatus_ && ("mismatched" != call.pairStatus_ || rescuePars.fromMismatched_) && assignOnOneEnd(wholeRead, res, "mismatchDirections")){
						continue;
					}
					MetaDataInName seqMeta;
					seqMeta.addMeta("direction", res.direction_);
					seqMeta.addMeta("forwardBar", res.forwardBar_);
					seqMeta.addMeta("reverseBar", res.reverseBar_);
					seqMeta.addMeta("fullBar", njh::pasteAsStr(res.forwardBar_, "+", res.reverseBar_));
					if(keepMismatchDirectionsWithMatchingBarcodes && "matched" == call.pairStatus_ && finder.trimMismatchDirectionRead(seq, res)){
						seqMeta.addMeta("fprimerStart", res.processRes_.fprimerStart_);
						seqMeta.addMeta("rprimerStart", res.processRes_.rprimerStart_);
						seqMeta.addMeta("sample", call.determinedSample_);
						seqMeta.addMeta("barcodeEnds", 2);
						seq.name_ += seqMeta.createMetaName();
						writeOut(call.determinedSample_, seq);
						++sampleCounts[call.determinedSample_];
						++mismatchDirectionsToSamples;
						if(0 == unusedBarcodeNames.count(call.determinedSample_)){
							++sampleDeterminedBothEnds;
						}
					}else{
						seqMeta.addMeta("frontSample", call.fwdCall_);
						seqMeta.addMeta("backSample", call.revCall_);
						seqMeta.addMeta("pairStatus", call.pairStatus_);
						seq.name_ += seqMeta.createMetaName();
						writeOut(mismatchDirectionsName, seq);
					}
					continue;
				}
				if(AdaptorReadOutcome::pass != res.outcome_){
					//these were never written out, but the end that made it through can still have a good barcode
					assignOnOneEnd(wholeRead, res, getOutcomeName(res.outcome_));
					continue;
				}
				const auto & processRes = res.processRes_;
				const auto & forwardBar = res.forwardBar_;
				const auto & reverseBar = res.reverseBar_;
				++fullPass;
				auto call = callBarcodePair(caller, forwardBar, reverseBar, true, false, mismatchConfirmIdentityCutOff);
				const auto & determinedSample = call.determinedSample_;
				const auto & fwdCall = call.fwdCall_;
				const auto & revCall = call.revCall_;
				const auto & pairStatus = call.pairStatus_;
				++barcodePairCounts[std::make_tuple(fwdCall, revCall, pairStatus)];
				++pairStatusCounts[pairStatus];
				//no sample carries an unused kit barcode, so these can only be misread barcodes, the background for the other mismatches
				if("mismatched" == pairStatus && (unusedBarcodeNames.count(fwdCall) || unusedBarcodeNames.count(revCall))){
					++pairStatusCounts["mismatchedWithUnusedBarcode"];
				}
				if("matched" != pairStatus && "ambiguousPair" != pairStatus && ("mismatched" != pairStatus || rescuePars.fromMismatched_) && assignOnOneEnd(wholeRead, res, pairStatus)){
					continue;
				}
				MetaDataInName seqMeta;
				seqMeta.addMeta("direction", res.direction_);
				seqMeta.addMeta("fprimerStart", processRes.fprimerStart_);
				seqMeta.addMeta("forwardBar", forwardBar);
				seqMeta.addMeta("reverseBar", reverseBar);
				seqMeta.addMeta("fullBar", njh::pasteAsStr(forwardBar, "+", reverseBar));
				seqMeta.addMeta("rprimerStart", processRes.rprimerStart_);
				seqMeta.addMeta("sample", determinedSample);
				//how many of the read's ends had the sample's barcode, for clustering to count (SeekDeep processClusters --minBothEndsBarcodeReads)
				if("matched" == pairStatus){
					seqMeta.addMeta("barcodeEnds", 2);
				}
				if("matched" != pairStatus){
					seqMeta.addMeta("fwdSample", fwdCall);
					seqMeta.addMeta("revSample", revCall);
					seqMeta.addMeta("pairStatus", pairStatus);
				}
				seq.name_ += seqMeta.createMetaName();
				writeOut(determinedSample, seq);
				++sampleCounts[determinedSample];
				if("matched" == pairStatus && 0 == unusedBarcodeNames.count(determinedSample)) {
					++sampleDeterminedBothEnds;
				}
				infoFiles.writePass(seq.name_, res);
			}
			if(infoFiles.on()){
				infoFiles.writeInfo(infoBuf.str());
				infoBuf.str("");
			}
		}
		std::lock_guard<std::mutex> lock(countsMut);
		counts.add(local);
		finderTotals.addCounts(finder);
		splitterTotals.addCounts(splitter);
	};
	njh::concurrent::runVoidFunctionThreaded(demuxReads, numThreads);
	auto & finder = finderTotals;
	auto & splitter = splitterTotals;
	auto & sampleCounts = counts.sampleCounts;
	auto & oneEndSampleCounts = counts.oneEndSampleCounts;
	auto & barcodePairCounts = counts.barcodePairCounts;
	auto & pairStatusCounts = counts.pairStatusCounts;
	auto & mismatchDirectionPairCounts = counts.mismatchDirectionPairCounts;
	auto & mismatchDirectionStatusCounts = counts.mismatchDirectionStatusCounts;
	auto & oneEndFromCounts = counts.oneEndFromCounts;
	auto & oneEndUsedCounts = counts.oneEndUsedCounts;
	auto & mismatchDirectionsToSamples = counts.mismatchDirectionsToSamples;
	auto & fullPass = counts.fullPass;
	auto & sampleDeterminedBothEnds = counts.sampleDeterminedBothEnds;
	auto & sampleDeterminedOneEnd = counts.sampleDeterminedOneEnd;
	auto & oneEndToUnusedBarcode = counts.oneEndToUnusedBarcode;
	auto & oneEndConflicting = counts.oneEndConflicting;
	auto & oneEndCutAtInnerAdaptor = counts.oneEndCutAtInnerAdaptor;
	auto & oneEndNoInsertBeforeInnerAdaptor = counts.oneEndNoInsertBeforeInnerAdaptor;
	if(splitPars.split_){
		splitter.writeCounts(setUp.pars_.directoryName_);
	}
	finder.writeBestDirectionCounts(setUp.pars_.directoryName_);
	finder.writeFailedCounts(setUp.pars_.directoryName_);
	{
		uint32_t unusedBarcodeReads = 0;
		for(const auto & sampleCount : sampleCounts){
			if(unusedBarcodeNames.count(sampleCount.first)){
				unusedBarcodeReads += sampleCount.second;
			}
		}
		std::vector<std::pair<std::string, uint32_t>> afterPasses{
			{"fullPass", fullPass},
			{"Undetermined", sampleCounts["Undetermined"]},
			{"sampleDetermined", sampleDeterminedBothEnds + sampleDeterminedOneEnd},
			{"sampleDeterminedBothEnds", sampleDeterminedBothEnds},
			{"sampleDeterminedOneEnd", sampleDeterminedOneEnd},
			{"oneEndToUnusedBarcode", oneEndToUnusedBarcode},
			{"oneEndSkippedOtherEndConflicts", oneEndConflicting},
			{"oneEndCutAtInnerAdaptor", oneEndCutAtInnerAdaptor},
			{"oneEndNoInsertBeforeInnerAdaptor", oneEndNoInsertBeforeInnerAdaptor},
			{"readsSplitAtInnerAdaptors", splitter.readsSplit()},
			{"piecesFromSplitReads", splitter.piecesFromSplits()},
			{"oneEndAssignedOnFront", oneEndUsedCounts["front"]},
			{"oneEndAssignedOnBack", oneEndUsedCounts["back"]},
			{"oneEndAssignedOnBoth", oneEndUsedCounts["both"]}};
		for(const auto & from : oneEndFromCounts){
			afterPasses.emplace_back(njh::pasteAsStr("oneEndFrom_", from.first), from.second);
		}
		afterPasses.emplace_back(mismatchedSamplesName, sampleCounts[mismatchedSamplesName]);
		afterPasses.emplace_back("oneEndUndetermined", pairStatusCounts["oneEndUndetermined"]);
		afterPasses.emplace_back("bothEndsUndetermined", pairStatusCounts["bothEndsUndetermined"]);
		afterPasses.emplace_back("ambiguousSamplePair", pairStatusCounts["ambiguousPair"]);
		afterPasses.emplace_back("unusedBarcode", unusedBarcodeReads);
		afterPasses.emplace_back("mismatchedSamplesWithUnusedBarcode", pairStatusCounts["mismatchedWithUnusedBarcode"]);
		afterPasses.emplace_back("mismatchedSamplesUnconfirmed", pairStatusCounts["mismatchedUnconfirmed"]);
		std::vector<std::pair<std::string, uint32_t>> afterMismatchDirections{
			{"mismatchDirectionsMatchingBarcodes", mismatchDirectionStatusCounts["matched"]},
			{"mismatchDirectionsMatchingUnusedBarcode", mismatchDirectionStatusCounts["matchedUnusedBarcode"]},
			{"mismatchDirectionsDifferentBarcodes", mismatchDirectionStatusCounts["mismatched"]},
			{"mismatchDirectionsDifferentBarcodesWithUnusedBarcode", mismatchDirectionStatusCounts["mismatchedWithUnusedBarcode"]},
			{"mismatchDirectionsDifferentBarcodesUnconfirmed", mismatchDirectionStatusCounts["mismatchedUnconfirmed"]},
			{"mismatchDirectionsOneEndUndetermined", mismatchDirectionStatusCounts["oneEndUndetermined"]},
			{"mismatchDirectionsBothEndsUndetermined", mismatchDirectionStatusCounts["bothEndsUndetermined"]},
			{"mismatchDirectionsAmbiguousSamplePair", mismatchDirectionStatusCounts["ambiguousPair"]},
			{"mismatchDirectionsWrittenToSamples", mismatchDirectionsToSamples}};
		finder.writeConditionCounts(setUp.pars_.directoryName_, afterPasses, afterMismatchDirections);
	}
	{
		double totalWritten = 0;
		for(const auto & sampleCount : sampleCounts){
			totalWritten += sampleCount.second;
		}
		OutputStream sampleCountsOut(njh::files::make_path(setUp.pars_.directoryName_, "sampleCounts.tab.txt"));
		sampleCountsOut << "sample\tcount\tfrac\toneEndCount" << std::endl;
		for(const auto & sampleCount : sampleCounts){
			sampleCountsOut << sampleCount.first
					<< "\t" << sampleCount.second
					<< "\t" << sampleCount.second/totalWritten
					<< "\t" << oneEndSampleCounts[sampleCount.first] << std::endl;
		}
	}
	{
		//every read with the same adaptor at both ends by the sample(s) each end was called as, most common first
		std::vector<std::pair<std::tuple<std::string, std::string, std::string, std::string>, uint32_t>> pairCounts(
				mismatchDirectionPairCounts.begin(), mismatchDirectionPairCounts.end());
		std::stable_sort(pairCounts.begin(), pairCounts.end(), [](const auto & left, const auto & right){
			return left.second > right.second;
		});
		uint32_t totalMismatchDirections = 0;
		for(const auto & pairCount : pairCounts){
			totalMismatchDirections += pairCount.second;
		}
		OutputStream pairCountsOut(njh::files::make_path(setUp.pars_.directoryName_, "mismatchDirectionsBarcodePairCounts.tab.txt"));
		pairCountsOut << "adaptor\tfrontSample\tbackSample\tstatus\tcount\tfrac" << std::endl;
		for(const auto & pairCount : pairCounts){
			pairCountsOut << std::get<0>(pairCount.first)
					<< "\t" << std::get<1>(pairCount.first)
					<< "\t" << std::get<2>(pairCount.first)
					<< "\t" << std::get<3>(pairCount.first)
					<< "\t" << pairCount.second
					<< "\t" << pairCount.second/static_cast<double>(totalMismatchDirections) << std::endl;
		}
	}
	writeBarcodePairCounts(setUp.pars_.directoryName_, barcodePairCounts, fullPass);

	return 0;
}


}  // namespace njhseq


