//
// Created by Nicholas Hathaway on 8/1/24.
//

#include <njhseq/IO/OutputStream.hpp>
#include <njhseq/objects/BioDataObject/reading.hpp>

#include "ampliconAnalysisRunner.hpp"

namespace njhseq {


struct RefDeterminedMaskingInfo {
	uint32_t ref_start{std::numeric_limits<uint32_t>::max()};
	uint32_t ref_segment_size{std::numeric_limits<uint32_t>::max()};
	uint32_t replacement_size{std::numeric_limits<uint32_t>::max()};
};

struct MaskingInfo {
	uint32_t seq_start{std::numeric_limits<uint32_t>::max()};
	uint32_t seq_segment_size{std::numeric_limits<uint32_t>::max()};
	uint32_t replacement_size{std::numeric_limits<uint32_t>::max()};

	char replacement = 'N';
};

int ampliconAnalysisRunner::determinePossibleMaskFromSeqs(
				const njh::progutils::CmdArgs & inputCommands) {

	char masking = 'N';
	OutOptions outOpts;
	bfs::path bedFnp;
	bfs::path genome2bitFnp;
	std::string regionMetaField = "regionUID";
	uint32_t numThreads = 2;
	ampliconAnalysisSetUp setUp(inputCommands);

	setUp.processReadInNames(seqSetUp::singleInFormatsAvailable_, true);
	setUp.setOption(regionMetaField, "--regionMetaField", "region Meta Field");
	setUp.setOption(masking, "--masking", "what is being used to mask");
	setUp.setOption(numThreads, "--threads", "Number of Threads");
	setUp.processWritingOptions(outOpts);
	setUp.setOption(bedFnp, "--bed", "The regions for extracting out the reference sequence", true);
	setUp.setOption(genome2bitFnp, "--genome2bitFnp", "genome 2bit Fnp", true);

	setUp.pars_.gapInfo_ = gapScoringParameters::genSemiGlobalQueryOnly(5,1);
	setUp.pars_.scoring_ = substituteMatrix::createDegenScoreMatrix(2, -2);

	setUp.finishSetUp(std::cout);
	if(setUp.pars_.verbose_) {
		std::cout << "Gap info: " << std::endl;
		std::cout << setUp.pars_.gapInfo_.toJson() << std::endl;
	}

	auto regions = getBeds(bedFnp);
	VecStr warnings;

	uint64_t maxLen = 0;

	std::unordered_map<std::string, std::shared_ptr<seqInfo>> ref_seqs;

	TwoBit::TwoBitFile tReader(genome2bitFnp);
	for(const auto & region : regions) {
		auto refSeq = std::make_shared<seqInfo>(GenomicRegion(*region).extractSeq(tReader));
		if(njh::in(region->name_, ref_seqs)) {
			warnings.emplace_back(njh::pasteAsStr("already have region with name: " + region->name_));
		}
		readVec::getMaxLength(refSeq, maxLen);
		ref_seqs[region->name_] = refSeq;
	}

	if(!warnings.empty()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "" << "\n";
		ss << njh::conToStr(warnings, "\n") << "\n";
		throw std::runtime_error{ss.str()};
	}

	aligner alignerObj(maxLen, setUp.pars_.gapInfo_, setUp.pars_.scoring_);

	std::regex pattern("(" + std::string(1, masking) + "+(-+" + std::string(1, masking) + "+)*)");
	{
		seqInfo seq;

		SeqInput reader(setUp.pars_.ioOptions_);
		reader.openIn();
		while(reader.readNextRead(seq)) {
			std::shared_ptr<seqInfo> refSeq;
			if(1 == ref_seqs.size()) {
					refSeq = ref_seqs.begin()->second;
			} else {
				if(!MetaDataInName::nameHasMetaData(seq.name_)) {
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error " << seq.name_ << " doesn't have meta, need to have meta data to match to which seq to determine masking" << "" << "\n";
					throw std::runtime_error{ss.str()};
				}
				MetaDataInName meta(seq.name_);
				if(!meta.containsMeta(regionMetaField)) {
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error " << "need to have " << regionMetaField <<", only found: " << njh::conToStr(njh::getVecOfMapKeys(meta.meta_), ",") << "" << "\n";
					throw std::runtime_error{ss.str()};
				}
				auto regionUID = meta.getMeta(regionMetaField	);
				refSeq = ref_seqs[regionUID];
			}

			alignerObj.alignCacheGlobal(refSeq, seq);
			alignerObj.rearrangeObjsGlobal(*refSeq, seq);
			std::vector<RefDeterminedMaskingInfo> maskingInfos;
			njh::PatPosFinder finder(pattern);
			auto masked_positions = finder.getPatPositions(alignerObj.alignObjectB_.seqBase_.seq_);
			alignerObj.alignObjectA_.seqBase_.outPutSeqAnsi(std::cout);
			alignerObj.alignObjectB_.seqBase_.outPutSeqAnsi(std::cout);
			std::cout << "masked_positions.size(): " << masked_positions.size() << std::endl;
			for(const auto & p : masked_positions) {
				std::cout << "\t" << "p.pos_: " << p.pos_ << std::endl;
				std::cout << "\t" << "p.pat_: " << p.pat_ << std::endl;
				std::cout << "\t" << "p.end(): " << p.end() << std::endl;
			}
		}
	}



	return 0;
}


int ampliconAnalysisRunner::maskRegionBasedOnRefSubRegions(
				const njh::progutils::CmdArgs & inputCommands) {

	ampliconAnalysisSetUp setUp(inputCommands);

	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);


	return 0;
}




}  // namespace njhseq
