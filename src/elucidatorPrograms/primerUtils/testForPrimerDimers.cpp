//
// Created by Nicholas Hathaway on 6/15/23.
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
#include <njhseq/concurrency/PairwisePairFactory.hpp>
#include <njhseq/objects/BioDataObject/reading.hpp>
#include <njhseq/objects/dataContainers/tables/TableReader.hpp>
#include <utility>


namespace njhseq {

int primerUtilsRunner::testForPrimerDimers(const njh::progutils::CmdArgs &inputCommands) {

	OutOptions outOpts;
	bfs::path primersFnp;
	seqSetUp setUp(inputCommands);
	comparison allowableErrors;
	allowableErrors.hqMismatches_ = 2;
	allowableErrors.lqMismatches_ = 3;
	allowableErrors.oneBaseIndel_ = 2;
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processReadInNames(njhseq::seqSetUp::singleInFormatsAvailable_, true);
	setUp.processWritingOptions(outOpts);
	setUp.processComparison(allowableErrors);
	setUp.setOption(primersFnp, "--primers", "Primers file", true);
	setUp.finishSetUp(std::cout);

	PrimersAndMids ids(primersFnp);

	if(ids.getTargets().empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error in reading in target primers file " << primersFnp << "\n";
		ss << "Make sure there is a line that starts with target in file" << "\n";
		throw std::runtime_error{ss.str()};
	}
	ids.initPrimerDeterminator();
	uint64_t maxSize = ids.pDeterminator_->getMaxPrimerSize();
	uint32_t totalReads = 0;
	{
		seqInfo seq;
		SeqInput reader(setUp.pars_.ioOptions_);
		reader.openIn();
		while(reader.readNextRead(seq)){
			readVec::getMaxLength(seq, maxSize);
			++totalReads;
		}
	}

	aligner alignerObj(maxSize, gapScoringParameters::genSemiGlobalQueryOnly(5,1), substituteMatrix::createDegenScoreMatrix(2,-4));
	seqInfo seq;
	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();

	struct PrimerLocalLoc{
		PrimerLocalLoc() = default;
		seqInfo alignA_;
		seqInfo alignB_;
		comparison comp_;
		alnInfoLocal infoLocal_;
		std::string primerName_;
		bool revComp_{false};
		[[nodiscard]] uint32_t get_a_end()const {
			return infoLocal_.localAStart_ + infoLocal_.localASize_;
		}

		[[nodiscard]] uint32_t get_b_end()const {
			return infoLocal_.localBStart_ + infoLocal_.localBSize_;
		}

		[[nodiscard]] bool aOverlapsOrAbuds(const PrimerLocalLoc &other) const {
			return ((get_a_end() > other.infoLocal_.localAStart_ && get_a_end() <= other.get_a_end()) ||
							get_a_end() == other.infoLocal_.localAStart_ ||
							(other.get_a_end() > infoLocal_.localAStart_ && other.get_a_end() <=get_a_end()) ||
							other.get_a_end() == infoLocal_.localAStart_ );
		}

		[[nodiscard]] uint32_t get_a_overlapDegree(const PrimerLocalLoc & other) const{
			if (((get_a_end() > other.infoLocal_.localAStart_ && get_a_end() <= other.get_a_end()) ||
					 (other.get_a_end() > infoLocal_.localAStart_ && other.get_a_end() <= get_a_end()))) {
				auto startOverLap = std::max(infoLocal_.localAStart_, other.infoLocal_.localAStart_);
				auto endOverLap = std::min(get_a_end(), other.get_a_end());
				return endOverLap - startOverLap;
			}
			return 0;
		}
	};

	struct PrimerLocalLocPair {
		PrimerLocalLocPair(PrimerLocalLoc front, PrimerLocalLoc end) : front_(std::move(front)), end_(std::move(end)) {

		}

		PrimerLocalLoc front_;
		PrimerLocalLoc end_;

		[[nodiscard]] double getAlnScoreSum()const{
			return front_.comp_.alnScore_ + end_.comp_.alnScore_;
		}

		[[nodiscard]] double getIdentitySum()const {
			return front_.comp_.distances_.eventBasedIdentity_ + end_.comp_.distances_.eventBasedIdentity_;
		}

		[[nodiscard]] double getIdentityHqSum()const {
			return front_.comp_.distances_.eventBasedIdentityHq_ + end_.comp_.distances_.eventBasedIdentityHq_;
		}
	};

	OutputStream out(outOpts);
	out
					<< "filename"
					<< "\t" << "name"
					<< "\t" << "seq"
					<< "\t" << "possibleDimers"
					<< "\t" << "numberOfFrontMatches"
					<< "\t" << "frontMatchesNames"
					<< "\t" << "numberOfEndMatches"
					<< "\t" << "endMatchesNames"
					<< "\t" << "overlapIDNumber"
					<< "\t" << "frontDimerPrimer-endDimerPrimer"
					<< "\t" << "degreeOfOverlap"
					<< "\t" << "frontDimerPrimer"
					<< "\t" << "frontDimerPrimer_isRevComp"
					<< "\t" << "frontDimerPrimer_seqStart"
					<< "\t" << "frontDimerPrimer_seqEnd"
					<< "\t" << "frontDimerPrimer_primerStart"
					<< "\t" << "frontDimerPrimer_primerEnd"

					<< "\t" << "endDimerPrimer"
					<< "\t" << "endDimerPrimer_isRevComp"
					<< "\t" << "endDimerPrimer_seqStart"
					<< "\t" << "endDimerPrimer_seqEnd"
					<< "\t" << "endDimerPrimer_primerStart"
					<< "\t" << "endDimerPrimer_primerEnd"

					<< std::endl;
	auto filenameOut = bfs::basename(bfs::path(setUp.pars_.ioOptions_.firstName_).replace_extension(""));
	while(reader.readNextRead(seq)){
		std::vector<PrimerLocalLoc> frontDeterminedPrimers;
		std::vector<PrimerLocalLoc> endDeterminedPrimers;

		for (const auto& currentPrimer : ids.pDeterminator_->primers_) {
			for(const auto & fwd : currentPrimer.second.fwds_){
				{
					alignerObj.alignCacheLocal(seq, fwd.info_);
					alignerObj.rearrangeObjsLocal(seq, fwd.info_);
					alignerObj.profileAlignment(seq, fwd.info_, false, true, false);
//					if("MN01689:44:000H57G73:1:11102:13521:7834 1:N:0:TGCTCTGTTGTC+AGTGGTTCTGCA"  == seq.name_ && fwd.info_.name_ == "Pf03-0221484-0221679"){
//						seq.outPutSeqAnsi(std::cout);
//						alignerObj.alignObjectA_.seqBase_.outPutSeqAnsi(std::cout);
//						alignerObj.alignObjectB_.seqBase_.outPutSeqAnsi(std::cout);
//						std::cout << alignerObj.comp_.toJson() << std::endl;
//						std::cout << "alignerObj.parts_.lHolder_.localAStart_: " << alignerObj.parts_.lHolder_.localAStart_ << std::endl;
//						std::cout << "alignerObj.parts_.lHolder_.localAStart_: " << alignerObj.parts_.lHolder_.localBStart_ << std::endl;
//						std::cout << "alignerObj.parts_.lHolder_.localAStart_ + alignerObj.parts_.lHolder_.localASize_: " << alignerObj.parts_.lHolder_.localAStart_ + alignerObj.parts_.lHolder_.localASize_ << std::endl;
//						std::cout << "alignerObj.parts_.lHolder_.localBStart_ + alignerObj.parts_.lHolder_.localBSize_: " << alignerObj.parts_.lHolder_.localBStart_ + alignerObj.parts_.lHolder_.localBSize_ << std::endl;
//						std::cout << "seq.seq_.size(): " << seq.seq_.size() << std::endl;
//						std::cout << "rev.infoRC_.seq_.size(): " << fwd.info_.seq_.size() << std::endl;
//						std::cout << ": " << std::endl;
//						exit(1);
//					}
					if(allowableErrors.passErrorProfile(alignerObj.comp_) &&  0 == alignerObj.parts_.lHolder_.localAStart_ && 0 == alignerObj.parts_.lHolder_.localBStart_){
						PrimerLocalLoc loc;
						loc.alignA_ = alignerObj.alignObjectA_.seqBase_;
						loc.alignB_ = alignerObj.alignObjectB_.seqBase_;
						loc.comp_ = alignerObj.comp_;
						loc.infoLocal_ = alignerObj.parts_.lHolder_;
						loc.primerName_ = fwd.info_.name_ + "_F";
						loc.revComp_ = false;
						frontDeterminedPrimers.emplace_back(loc);
					}
					if(allowableErrors.passErrorProfile(alignerObj.comp_) && seq.seq_.size() == (alignerObj.parts_.lHolder_.localAStart_ + alignerObj.parts_.lHolder_.localASize_) && fwd.info_.seq_.size() == (alignerObj.parts_.lHolder_.localBStart_ + alignerObj.parts_.lHolder_.localBSize_)){
						PrimerLocalLoc loc;
						loc.alignA_ = alignerObj.alignObjectA_.seqBase_;
						loc.alignB_ = alignerObj.alignObjectB_.seqBase_;
						loc.comp_ = alignerObj.comp_;
						loc.infoLocal_ = alignerObj.parts_.lHolder_;
						loc.primerName_ = fwd.info_.name_ + "_F";
						loc.revComp_ = false;
						endDeterminedPrimers.emplace_back(loc);
					}
				}
				{
					alignerObj.alignCacheLocal(seq, fwd.infoRC_);
					alignerObj.rearrangeObjsLocal(seq, fwd.infoRC_);

					alignerObj.profileAlignment(seq, fwd.infoRC_, false, true, false);
					if(allowableErrors.passErrorProfile(alignerObj.comp_) && 0 == alignerObj.parts_.lHolder_.localAStart_ && 0 == alignerObj.parts_.lHolder_.localBStart_){
						PrimerLocalLoc loc;
						loc.alignA_ = alignerObj.alignObjectA_.seqBase_;
						loc.alignB_ = alignerObj.alignObjectB_.seqBase_;
						loc.comp_ = alignerObj.comp_;
						loc.infoLocal_ = alignerObj.parts_.lHolder_;
						loc.primerName_ = fwd.infoRC_.name_ + "_F";
						loc.revComp_ = true;
						frontDeterminedPrimers.emplace_back(loc);
					}
					if(allowableErrors.passErrorProfile(alignerObj.comp_) && seq.seq_.size() == (alignerObj.parts_.lHolder_.localAStart_ + alignerObj.parts_.lHolder_.localASize_) && fwd.infoRC_.seq_.size() == (alignerObj.parts_.lHolder_.localBStart_ + alignerObj.parts_.lHolder_.localBSize_)){
						PrimerLocalLoc loc;
						loc.alignA_ = alignerObj.alignObjectA_.seqBase_;
						loc.alignB_ = alignerObj.alignObjectB_.seqBase_;
						loc.comp_ = alignerObj.comp_;
						loc.infoLocal_ = alignerObj.parts_.lHolder_;
						loc.primerName_ = fwd.infoRC_.name_ + "_F";
						loc.revComp_ = true;
						endDeterminedPrimers.emplace_back(loc);
					}
				}
			}
			for(const auto & rev : currentPrimer.second.revs_){
				{
					alignerObj.alignCacheLocal(seq, rev.info_);
					alignerObj.rearrangeObjsLocal(seq, rev.info_);
					alignerObj.profileAlignment(seq, rev.info_, false, true, false);
					if(allowableErrors.passErrorProfile(alignerObj.comp_) && 0 == alignerObj.parts_.lHolder_.localAStart_ && 0 == alignerObj.parts_.lHolder_.localBStart_){
						PrimerLocalLoc loc;
						loc.alignA_ = alignerObj.alignObjectA_.seqBase_;
						loc.alignB_ = alignerObj.alignObjectB_.seqBase_;
						loc.comp_ = alignerObj.comp_;
						loc.infoLocal_ = alignerObj.parts_.lHolder_;
						loc.primerName_ = rev.info_.name_ + "_R";
						loc.revComp_ = false;
						frontDeterminedPrimers.emplace_back(loc);
					}
					if(allowableErrors.passErrorProfile(alignerObj.comp_) && seq.seq_.size() == (alignerObj.parts_.lHolder_.localAStart_ + alignerObj.parts_.lHolder_.localASize_) && rev.info_.seq_.size() == (alignerObj.parts_.lHolder_.localBStart_ + alignerObj.parts_.lHolder_.localBSize_)){
						PrimerLocalLoc loc;
						loc.alignA_ = alignerObj.alignObjectA_.seqBase_;
						loc.alignB_ = alignerObj.alignObjectB_.seqBase_;
						loc.comp_ = alignerObj.comp_;
						loc.infoLocal_ = alignerObj.parts_.lHolder_;
						loc.primerName_ = rev.info_.name_ + "_R";
						loc.revComp_ = false;
						endDeterminedPrimers.emplace_back(loc);
					}
				}
				{
					alignerObj.alignCacheLocal(seq, rev.infoRC_);
					alignerObj.rearrangeObjsLocal(seq, rev.infoRC_);

					alignerObj.profileAlignment(seq, rev.infoRC_, false, true, false);
//					if(seq.name_ == "MN01689:44:000H57G73:1:11102:11431:17330 1:N:0:TGCTCTGTTGTC+AGTGGTTCTGCA" && rev.infoRC_.name_ == "Pf14-0294537-0294730"){
//						seq.outPutSeqAnsi(std::cout);
//						alignerObj.alignObjectA_.seqBase_.outPutSeqAnsi(std::cout);
//						alignerObj.alignObjectB_.seqBase_.outPutSeqAnsi(std::cout);
//						std::cout << alignerObj.comp_.toJson() << std::endl;
//						std::cout << "alignerObj.parts_.lHolder_.localAStart_: " << alignerObj.parts_.lHolder_.localAStart_ << std::endl;
//						std::cout << "alignerObj.parts_.lHolder_.localAStart_: " << alignerObj.parts_.lHolder_.localBStart_ << std::endl;
//						std::cout << "alignerObj.parts_.lHolder_.localAStart_ + alignerObj.parts_.lHolder_.localASize_: " << alignerObj.parts_.lHolder_.localAStart_ + alignerObj.parts_.lHolder_.localASize_ << std::endl;
//						std::cout << "alignerObj.parts_.lHolder_.localBStart_ + alignerObj.parts_.lHolder_.localBSize_: " << alignerObj.parts_.lHolder_.localBStart_ + alignerObj.parts_.lHolder_.localBSize_ << std::endl;
//						std::cout << "seq.seq_.size(): " << seq.seq_.size() << std::endl;
//						std::cout << "rev.infoRC_.seq_.size(): " << rev.infoRC_.seq_.size() << std::endl;
//						std::cout << ": " << std::endl;
//						exit(1);
//					}
					if(allowableErrors.passErrorProfile(alignerObj.comp_) && 0 == alignerObj.parts_.lHolder_.localAStart_ && 0 == alignerObj.parts_.lHolder_.localBStart_){
						PrimerLocalLoc loc;
						loc.alignA_ = alignerObj.alignObjectA_.seqBase_;
						loc.alignB_ = alignerObj.alignObjectB_.seqBase_;
						loc.comp_ = alignerObj.comp_;
						loc.infoLocal_ = alignerObj.parts_.lHolder_;
						loc.primerName_ = rev.infoRC_.name_ + "_R";
						loc.revComp_ = true;
						frontDeterminedPrimers.emplace_back(loc);
					}
					if(allowableErrors.passErrorProfile(alignerObj.comp_) && seq.seq_.size() == (alignerObj.parts_.lHolder_.localAStart_ + alignerObj.parts_.lHolder_.localASize_) && rev.infoRC_.seq_.size() == (alignerObj.parts_.lHolder_.localBStart_ + alignerObj.parts_.lHolder_.localBSize_)){
						PrimerLocalLoc loc;
						loc.alignA_ = alignerObj.alignObjectA_.seqBase_;
						loc.alignB_ = alignerObj.alignObjectB_.seqBase_;
						loc.comp_ = alignerObj.comp_;
						loc.infoLocal_ = alignerObj.parts_.lHolder_;
						loc.primerName_ = rev.infoRC_.name_ + "_R";
						loc.revComp_ = true;
						endDeterminedPrimers.emplace_back(loc);
					}
				}
			}
		}
		std::vector<PrimerLocalLocPair> overlaps;
		for(const auto & front : frontDeterminedPrimers){
			for(const auto & end : endDeterminedPrimers){
				if(front.aOverlapsOrAbuds(end)){
					overlaps.emplace_back(front, end);
				}
			}
		}
		std::vector<PrimerLocalLocPair> bestOverlaps;

		if (overlaps.size() > 1) {
			double bestAlnScore = 0;
			double bestIdentity = 0;
			for (const auto &overlap: overlaps) {
				if (bestAlnScore == overlap.getAlnScoreSum()) {
					if (overlap.getIdentitySum() == bestIdentity) {
						bestOverlaps.emplace_back(overlap);
					} else if (overlap.getIdentitySum() >= bestIdentity) {
						bestAlnScore = overlap.getAlnScoreSum();
						bestIdentity = overlap.getIdentitySum();
						bestOverlaps.clear();
						bestOverlaps.emplace_back(overlap);
					}
				} else if (overlap.getAlnScoreSum() > bestAlnScore) {
					bestAlnScore = overlap.getAlnScoreSum();
					bestIdentity = overlap.getIdentitySum();
					bestOverlaps.clear();
					bestOverlaps.emplace_back(overlap);
				}
			}
		} else if (1 == overlaps.size()) {
			bestOverlaps = overlaps;
		} else {
/*			if(!frontDeterminedPrimers.empty() || !endDeterminedPrimers.empty()){
				seq.outPutSeqAnsi(std::cout);
				std::cout << "bestOverlaps.size(): " << bestOverlaps.size() << std::endl;
				std::cout << "overlaps.size(): " << overlaps.size() << std::endl;
				std::cout << "frontDeterminedPrimers.size(): " << frontDeterminedPrimers.size() << std::endl;
				std::cout << "endDeterminedPrimers.size(): " << endDeterminedPrimers.size() << std::endl;
				for(const auto & front : frontDeterminedPrimers){
					std::cout << seq.getSeqAnsi() << std::endl;
					std::cout << front.alignB_.getSeqAnsi() << std::endl;
				}
				for(const auto & end : endDeterminedPrimers){
					std::cout << seq.getSeqAnsi() << std::endl;
					std::cout << std::string(end.infoLocal_.localAStart_, ' ') << end.alignB_.getSeqAnsi()
										<< std::endl;
				}
			}*/
		}
		VecStr frontNames;
		VecStr ednNames;
		for(const auto & front : frontDeterminedPrimers){
			frontNames.emplace_back(front.primerName_);
		}
		for(const auto & end : endDeterminedPrimers){
			ednNames.emplace_back(end.primerName_);
		}
		if(bestOverlaps.empty()){
			out << filenameOut << "\t" << seq.name_
					<< "\t" << seq.seq_
					<< "\t" << bestOverlaps.size()
					<< "\t" << frontDeterminedPrimers.size()
					<< "\t" << njh::conToStr(frontNames, ";")
					<< "\t" << endDeterminedPrimers.size()
					<< "\t" << njh::conToStr(ednNames, ";")
					<< "\t" << ""
					<< "\t" << ""
					<< "\t" << ""
					<< "\t" << ""
					<< "\t" << ""
					<< "\t" << ""
					<< "\t" << ""
					<< "\t" << ""
					<< "\t" << ""

					<< "\t" << ""
					<< "\t" << ""
					<< "\t" << ""
					<< "\t" << ""
					<< "\t" << ""
					<< "\t" << ""

					<< std::endl;
		} else {
			uint32_t overlapCount = 0;
			for (const auto &overlap: bestOverlaps) {
				out << filenameOut << "\t" << seq.name_
						<< "\t" << seq.seq_
						<< "\t" << bestOverlaps.size()
						<< "\t" << frontDeterminedPrimers.size()
						<< "\t" << njh::conToStr(frontNames, ";")
						<< "\t" << endDeterminedPrimers.size()
						<< "\t" << njh::conToStr(ednNames, ";")
						<< "\t" << overlapCount
						<< "\t" << overlap.front_.primerName_ << "-" << overlap.end_.primerName_
						<< "\t" << overlap.front_.get_a_overlapDegree(overlap.end_)
						<< "\t" << overlap.front_.primerName_
						<< "\t" << njh::boolToStr(overlap.front_.revComp_)
						<< "\t" << overlap.front_.infoLocal_.localAStart_
						<< "\t" << overlap.front_.get_a_end()
						<< "\t" << overlap.front_.infoLocal_.localBStart_
						<< "\t" << overlap.end_.get_b_end()

						<< "\t" << overlap.end_.primerName_
						<< "\t" << njh::boolToStr(overlap.end_.revComp_)
						<< "\t" << overlap.end_.infoLocal_.localAStart_
						<< "\t" << overlap.end_.get_a_end()
						<< "\t" << overlap.end_.infoLocal_.localBStart_
						<< "\t" << overlap.end_.get_b_end()

						<< std::endl;
				++overlapCount;
			}
		}
		if(setUp.pars_.debug_){
			seq.outPutSeqAnsi(std::cout);
			std::cout << "bestOverlaps.size(): " << bestOverlaps.size() << std::endl;
			std::cout << "overlaps.size(): " << overlaps.size() << std::endl;
			std::cout << "frontDeterminedPrimers.size(): " << frontDeterminedPrimers.size() << std::endl;
			std::cout << "endDeterminedPrimers.size(): " << endDeterminedPrimers.size() << std::endl;
			for(const auto & overlap : bestOverlaps){
				seq.outPutSeqAnsi(std::cout);
				std::cout << overlap.front_.alignB_.getSeqAnsi() << std::endl;
				std::cout << std::string(overlap.end_.infoLocal_.localAStart_, ' ') << overlap.end_.alignB_.getSeqAnsi() << std::endl;
				std::cout << overlap.front_.primerName_ << " " << njh::colorBool(overlap.front_.revComp_) << " " << overlap.front_.infoLocal_.localAStart_ << ":" << overlap.front_.get_a_end() << std::endl;
				std::cout << overlap.end_.primerName_ << " " << njh::colorBool(overlap.end_.revComp_)<< " " << overlap.end_.infoLocal_.localAStart_ << ":" << overlap.end_.get_a_end() << std::endl;
				std::cout << "overlap.front_.get_a_overlapDegree(overlap.end_): " << overlap.front_.get_a_overlapDegree(overlap.end_)  << std::endl;
				std::cout << "overlap.end_.comp_.alnScore_ + overlap.front_.comp_.alnScore_: " << overlap.end_.comp_.alnScore_ + overlap.front_.comp_.alnScore_  << std::endl;
				std::cout << "overlap.end_.comp_.distances_.eventBasedIdentityHq_ + overlap.front_.comp_.distances_.eventBasedIdentityHq_: " << overlap.end_.comp_.distances_.eventBasedIdentityHq_ + overlap.front_.comp_.distances_.eventBasedIdentityHq_  << std::endl;
				std::cout << "overlap.end_.comp_.distances_.eventBasedIdentity_ + overlap.front_.comp_.distances_.eventBasedIdentity_: " << overlap.end_.comp_.distances_.eventBasedIdentity_ + overlap.front_.comp_.distances_.eventBasedIdentity_  << std::endl;
				std::cout << std::endl;
			}
			std::cout << std::endl;
		}
	}



	return 0;
}

} // namespace njhseq

