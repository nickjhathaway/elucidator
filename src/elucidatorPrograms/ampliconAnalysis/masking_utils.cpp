//
// Created by Nicholas Hathaway on 8/1/24.
//


#include <njhseq/objects/BioDataObject/reading.hpp>
#include <njhseq/objects/BioDataObject/BioRecordsUtils/BedUtility.hpp>

#include "ampliconAnalysisRunner.hpp"

namespace njhseq {


struct RefDeterminedMaskingInfo {
	RefDeterminedMaskingInfo() = default;
	uint32_t ref_start{std::numeric_limits<uint32_t>::max()};
	uint32_t ref_segment_size{std::numeric_limits<uint32_t>::max()};
	uint32_t replacement_size{std::numeric_limits<uint32_t>::max()};

	[[nodiscard]] Json::Value toJson() const {
		Json::Value ret;
		ret["class"] = njh::getTypeName(*this);
		ret["ref_start"] = ref_start;
		ret["ref_segment_size"] = ref_segment_size;
		ret["replacement_size"] = replacement_size;
		return ret;
	}
};

struct MaskingInfo {
	MaskingInfo() = default;
	uint32_t seq_start{std::numeric_limits<uint32_t>::max()};
	uint32_t seq_segment_size{std::numeric_limits<uint32_t>::max()};
	uint32_t replacement_size{std::numeric_limits<uint32_t>::max()};

	char replacement = 'N';

	[[nodiscard]] Json::Value toJson() const {
		Json::Value ret;
		ret["class"] = njh::getTypeName(*this);
		ret["seq_start"] = seq_start;
		ret["seq_segment_size"] = seq_segment_size;
		ret["replacement_size"] = replacement_size;
		ret["replacement"] = replacement;
		return ret;
	}

	static void applyMasks(seqInfo & seq, std::vector<MaskingInfo> maskingInfos) {
		njh::sort(maskingInfos, [](const MaskingInfo & m1, const MaskingInfo & m2) {
			return m1.seq_start < m2.seq_start;
		});
		for (auto & mask : iter::reversed(maskingInfos)) {
			if(mask.seq_start > len(seq)) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "masking start great than length, mask.ref_start: " << mask.seq_start << ", reference length: " << len(seq) << "\n";
				throw std::runtime_error{ss.str()};
			}
			if(mask.seq_start + mask.seq_segment_size > len(seq)) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "masking start + length great than length, mask.ref_start + mask.ref_segment_size: " << mask.seq_start + mask.seq_segment_size << ", reference length: " << len(seq) << "\n";
				throw std::runtime_error{ss.str()};
			}
			seq.replace(mask.seq_start, mask.seq_segment_size, std::string(mask.replacement_size, mask.replacement));
		}
	}
};

int ampliconAnalysisRunner::determinePossibleMaskFromSeqs(
				const njh::progutils::CmdArgs & inputCommands) {
	double freqInclusiveCutOff = 1;
	char masking = 'N';
	OutOptions outOpts("", ".tsv");
	bfs::path bedFnp;
	bfs::path genome2bitFnp;
	std::string regionMetaField = "regionUID";
	uint32_t numThreads = 2;
	ampliconAnalysisSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processReadInNames(seqSetUp::singleInFormatsAvailable_, true);
	setUp.setOption(freqInclusiveCutOff, "--freqInclusiveCutOff", "freq Inclusive Cut Off");
	setUp.setOption(regionMetaField, "--regionMetaField", "region Meta Field");
	setUp.setOption(masking, "--masking", "what is being used to mask");
	setUp.setOption(numThreads, "--threads", "Number of Threads");
	setUp.processWritingOptions(outOpts);
	setUp.setOption(bedFnp, "--bed", "The regions for extracting out the reference sequence", true);
	setUp.setOption(genome2bitFnp, "--genome2bitFnp", "genome 2bit Fnp", true);

	setUp.pars_.gapInfo_ = gapScoringParameters::genSemiGlobalQueryOnly(5,1);
	setUp.pars_.scoring_ = substituteMatrix::createDegenScoreMatrixLessN(2, -2);

	setUp.finishSetUp(std::cout);
	if(setUp.pars_.verbose_) {
		std::cout << "Gap info: " << std::endl;
		std::cout << setUp.pars_.gapInfo_.toJson() << std::endl;
	}

	auto regions = getBeds(bedFnp);
	VecStr warnings;

	uint64_t maxLen = 0;

	std::unordered_map<std::string, std::shared_ptr<seqInfo>> ref_seqs;
	std::unordered_map<std::string, std::shared_ptr<Bed6RecordCore>> ref_regions;

	TwoBit::TwoBitFile tReader(genome2bitFnp);
	for(const auto & region : regions) {
		auto refSeq = std::make_shared<seqInfo>(GenomicRegion(*region).extractSeq(tReader));
		if(njh::in(region->name_, ref_seqs)) {
			warnings.emplace_back(njh::pasteAsStr("already have region with name: " + region->name_));
		}
		readVec::getMaxLength(refSeq, maxLen);
		ref_seqs[region->name_] = refSeq;
		ref_regions[region->name_] = region;
	}

	if(!warnings.empty()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "" << "\n";
		ss << njh::conToStr(warnings, "\n") << "\n";
		throw std::runtime_error{ss.str()};
	}
	OutputStream out(outOpts	);
	aligner alignerObj(maxLen, setUp.pars_.gapInfo_, setUp.pars_.scoring_);


	std::regex pattern("(" + std::string(1, masking) + "+(-+" + std::string(1, masking) + "+)*)");

	std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>> determinedMaskingCount;

	std::unordered_map<std::string, RefDeterminedMaskingInfo> maskings;

	std::unordered_map<std::string, uint32_t> regionTotals;

	{
		seqInfo seq;

		SeqInput reader(setUp.pars_.ioOptions_);
		reader.openIn();
		while(reader.readNextRead(seq)) {
			std::shared_ptr<seqInfo> refSeq;
			std::string regionUID = "region";
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
				regionUID = meta.getMeta(regionMetaField);
				refSeq = ref_seqs[regionUID];
			}

			alignerObj.alignCacheGlobal(refSeq, seq);
			alignerObj.rearrangeObjsGlobal(*refSeq, seq);
			std::vector<RefDeterminedMaskingInfo> maskingInfos;
			njh::PatPosFinder finder(pattern);
			auto masked_positions = finder.getPatPositions(alignerObj.alignObjectB_.seqBase_.seq_);
			if(setUp.pars_.debug_) {
				alignerObj.alignObjectA_.seqBase_.outPutSeqAnsi(std::cout);
				alignerObj.alignObjectB_.seqBase_.outPutSeqAnsi(std::cout);
				std::cout << "masked_positions.size(): " << masked_positions.size() << std::endl;
			}

			for(const auto & p : masked_positions) {
				auto refSegment = alignerObj.alignObjectA_.seqBase_.seq_.substr(p.pos_, p.pat_.size());
				refSegment = njh::replaceString(refSegment,"-", "");
				RefDeterminedMaskingInfo mask;
				mask.ref_segment_size = refSegment.size();
				mask.replacement_size = p.pat_.size();
				mask.ref_start = alignerObj.getSeqPosForAlnAPos(p.pos_);
				if(setUp.pars_.debug_) {
					std::cout << "\t" << "p.pos_: " << p.pos_ << std::endl;
					std::cout << "\t" << "p.pat_: " << p.pat_ << std::endl;
					std::cout << "\t" << "p.end(): " << p.end() << std::endl;
					std::cout << "\t" << njh::json::writeAsOneLine(mask.toJson()) << std::endl;
				}
				auto maskingID = njh::pasteAsStr(mask.ref_start, "-", mask.ref_segment_size, "-", mask.replacement_size);

				maskings[maskingID] = mask;
				determinedMaskingCount[regionUID][maskingID] += static_cast<uint32_t>(std::round(seq.cnt_));
			}
			regionTotals[regionUID] += static_cast<uint32_t>(std::round(seq.cnt_));
		}
	}

	table determinedMaskingTable(determinedMaskingCount, VecStr{regionMetaField, "masking", "count"});
	determinedMaskingTable.sortTable(regionMetaField, false);

	VecStr regionTotalsCol;
	VecStr regionFreqCol;
	for(const auto & row : determinedMaskingTable) {
		auto total = static_cast<double>(regionTotals[row[determinedMaskingTable.getColPos(regionMetaField)]]);
		auto freq = njh::StrToNumConverter::stoToNum<uint32_t>(row[determinedMaskingTable.getColPos("count")])/total;
		regionTotalsCol.emplace_back(estd::to_string(total));
		regionFreqCol.emplace_back(estd::to_string(freq));
	}
	determinedMaskingTable.addColumn(regionTotalsCol, "regionTotals");
	determinedMaskingTable.addColumn(regionFreqCol, "freq");
	if(setUp.pars_.debug_) {
		determinedMaskingTable.outPutContents(std::cout, "\t");
	}
	out <<
		"#chrom\tstart\tend\t" << regionMetaField << "\tlength\tstrand";
	out << "\t" << njh::conToStr(VecStr{"ref_start", "ref_segment_size", "replacement_size", "count", "freq", "regionTotal"}, "\t") << std::endl;

	for(const auto & row : determinedMaskingTable) {
		out << ref_regions[row[determinedMaskingTable.getColPos(regionMetaField)]]->toDelimStr();
		auto maskingToks = tokenizeString(row[determinedMaskingTable.getColPos("masking")], "-");
		out << "\t" << maskingToks[0]
        << "\t" << maskingToks[1]
        << "\t" << maskingToks[2]
        << "\t" << row[determinedMaskingTable.getColPos("count")]
        << "\t" << row[determinedMaskingTable.getColPos("freq")]
        << "\t" << row[determinedMaskingTable.getColPos("regionTotals")]
        << std::endl;
	}
	return 0;
}


int ampliconAnalysisRunner::maskRegionBasedOnRefSubRegions(
				const njh::progutils::CmdArgs & inputCommands) {

	bfs::path maskingInfoRegionFnp;
	char masking = 'N';
	bfs::path genome2bitFnp;
	std::string regionMetaField = "regionUID";
	uint32_t numThreads = 2;

	bfs::path outputMaskingFnp;
	ampliconAnalysisSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processDefaultReader(seqSetUp::singleInFormatsAvailable_, true);
	setUp.setOption(regionMetaField, "--regionMetaField", "region Meta Field");
	setUp.setOption(maskingInfoRegionFnp, "--maskingInfoRegionFnp", "The regions along with masking relative to that region info", true);
	setUp.setOption(genome2bitFnp, "--genome2bitFnp", "genome 2bit Fnp", true);
	setUp.setOption(masking, "--masking", "what is being used to mask");
	setUp.setOption(numThreads, "--threads", "Number of Threads");
	setUp.setOption(outputMaskingFnp, "--outputMaskingFnp", "output Masking Fnp, if left blank, info will not be written");

	setUp.pars_.gapInfo_ = gapScoringParameters::genSemiGlobal(5,1);
	setUp.pars_.scoring_ = substituteMatrix::createDegenScoreMatrixLessN(2, -2);

	setUp.finishSetUp(std::cout);

	OutOptions outOptions(outputMaskingFnp, ".tsv");
	outOptions.transferOverwriteOpts(setUp.pars_.ioOptions_.out_);

	if(setUp.pars_.verbose_) {
		std::cout << "Gap info: " << std::endl;
		std::cout << setUp.pars_.gapInfo_.toJson() << std::endl;
	}

	auto all_regions = getBeds(maskingInfoRegionFnp);
	BedUtility::coordSort(all_regions);
	std::vector<std::shared_ptr<Bed6RecordCore>> uniqueRegions;
	for(const auto & region : all_regions) {
		if(uniqueRegions.empty() || !uniqueRegions.back()->sameRegion(*region)) {
			uniqueRegions.emplace_back(region);
		}
	}
	table maskingInfoRegionTab(maskingInfoRegionFnp, "\t", true);
	maskingInfoRegionTab.checkForColumnsThrow(VecStr{regionMetaField, "ref_start", "ref_segment_size", "replacement_size"}, __PRETTY_FUNCTION__);
	std::unordered_map<std::string, std::vector<RefDeterminedMaskingInfo>> referenceMasking;
	for(const auto & row : maskingInfoRegionTab) {
		RefDeterminedMaskingInfo mask;
		mask.ref_start = njh::StrToNumConverter::stoToNum<uint32_t>(row[maskingInfoRegionTab.getColPos("ref_start")]);
		mask.ref_segment_size = njh::StrToNumConverter::stoToNum<uint32_t>(row[maskingInfoRegionTab.getColPos("ref_segment_size")]);
		mask.replacement_size = njh::StrToNumConverter::stoToNum<uint32_t>(row[maskingInfoRegionTab.getColPos("replacement_size")]);
		auto regionUID = row[maskingInfoRegionTab.getColPos(regionMetaField)];
		referenceMasking[regionUID].emplace_back(mask);
	}

	for(auto & region : referenceMasking) {
		njh::sort(region.second, [](const RefDeterminedMaskingInfo & mask1, const RefDeterminedMaskingInfo & mask2) {
			return mask1.ref_start < mask2.ref_start;
		});
	}

	VecStr warnings;

	uint64_t maxLen = 0;

	std::unordered_map<std::string, std::shared_ptr<seqInfo>> ref_seqs;
	std::unordered_map<std::string, std::shared_ptr<Bed6RecordCore>> ref_regions;

	TwoBit::TwoBitFile tReader(genome2bitFnp);
	for(const auto & region : uniqueRegions) {
		auto refSeq = std::make_shared<seqInfo>(GenomicRegion(*region).extractSeq(tReader));
		if(njh::in(region->name_, ref_seqs)) {
			warnings.emplace_back(njh::pasteAsStr("already have region with name: " + region->name_));
		}
		readVec::getMaxLength(refSeq, maxLen);


		if(njh::in(region->name_, referenceMasking)) {
			for(const auto & mask : iter::reversed(referenceMasking[region->name_])) {
				if(mask.ref_start > len(*refSeq)) {
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error " << "masking start great than length, mask.ref_start: " << mask.ref_start << ", reference length: " << len(*refSeq) << "\n";
					throw std::runtime_error{ss.str()};
				}
				if(mask.ref_start + mask.ref_segment_size > len(*refSeq)) {
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error " << "masking start + length great than length, mask.ref_start + mask.ref_segment_size: " << mask.ref_start + mask.ref_segment_size << ", reference length: " << len(*refSeq) << "\n";
					throw std::runtime_error{ss.str()};
				}
				refSeq->seq_.replace(mask.ref_start, mask.ref_segment_size, std::string(mask.replacement_size, masking));
			}
		}
		if(setUp.pars_.debug_) {
			refSeq->outPutSeqAnsi(std::cout);
		}

		ref_seqs[region->name_] = refSeq;
		ref_regions[region->name_] = region;
	}

	if(!warnings.empty()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "" << "\n";
		ss << njh::conToStr(warnings, "\n") << "\n";
		throw std::runtime_error{ss.str()};
	}
	maxLen *=2;
	aligner alignerObj(maxLen, setUp.pars_.gapInfo_, setUp.pars_.scoring_);
	std::shared_ptr<aligner> debug_alignerObj;
	if(setUp.pars_.debug_) {
		debug_alignerObj = std::make_shared<aligner>(maxLen, setUp.pars_.gapInfo_, substituteMatrix(2,-2));
	}

	// std::regex pattern("(" + std::string(1, masking) + "+(-+" + std::string(1, masking) + "+)*)");
	// std::regex pattern("((" + std::string(1, masking) + "+(-+" + std::string(1, masking) + "+)*|-+" + std::string(1, masking) + "+(-+" + std::string(1, masking) + "+)*))");
	//std::cout << "((" + std::string(1, masking) + "+(-+" + std::string(1, masking) + "+)*|-+" + std::string(1, masking) + "+(-+" + std::string(1, masking) + "+)*))" << std::endl;
	// std::regex pattern(R"((N+(-+N+)*|-+N+(-+N+)*))");
	// std::regex pattern(R"(-*N+(-+N+)*-*)");
	std::regex pattern("(-*" + std::string(1, masking) + "+(-+" + std::string(1, masking) + "+)*-*)");
	table outMaskInfo(VecStr{regionMetaField, "name", "seq_start", "seq_segment_size", "replacement_size", "replacement"});

	{
		seqInfo seq;

		SeqInput reader(setUp.pars_.ioOptions_);
		reader.openIn();

		SeqOutput writer(setUp.pars_.ioOptions_);
		writer.openOut();

		while(reader.readNextRead(seq)) {
			std::shared_ptr<seqInfo> refSeq;
			std::string regionUID = "region";
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
			regionUID = meta.getMeta(regionMetaField);
			//no masking info for this seq, will write out without any masking applied;
			if(njh::in(regionUID, ref_seqs)) {
				refSeq = ref_seqs[regionUID];
				alignerObj.alignCacheGlobal(refSeq, seq);
				alignerObj.rearrangeObjsGlobal(*refSeq, seq);
				std::vector<MaskingInfo> maskingInfos;
				njh::PatPosFinder finder(pattern);
				auto masked_positions = finder.getPatPositions(alignerObj.alignObjectA_.seqBase_.seq_);
				for(const auto & position : masked_positions) {
					MaskingInfo mask;
					mask.replacement = masking;
					mask.seq_start = alignerObj.getSeqPosForAlnBPos(position.pos_);
					auto seqSubSeq = alignerObj.alignObjectB_.seqBase_.getSubRead(position.pos_, position.pat_.size());
					auto refSubSeq = alignerObj.alignObjectA_.seqBase_.getSubRead(position.pos_, position.pat_.size());
					if(setUp.pars_.debug_) {
						std::cout << "position.pat_: " << position.pat_ << std::endl;
						seqSubSeq.outPutSeqAnsi(std::cout);
						refSubSeq.outPutSeqAnsi(std::cout);
					}


					seqSubSeq.removeGaps();
					refSubSeq.removeGaps();

					mask.seq_segment_size = seqSubSeq.seq_.size();
					mask.replacement_size = refSubSeq.seq_.size();
					maskingInfos.emplace_back(mask);
					outMaskInfo.addRow(
						regionUID,
						seq.name_,
						mask.seq_start,
						mask.seq_segment_size,
						mask.replacement_size,
						mask.replacement
					);
				}

				MaskingInfo::applyMasks(seq, maskingInfos);
				if(setUp.pars_.debug_) {
					alignerObj.alignObjectA_.seqBase_.outPutSeqAnsi(std::cout);
					alignerObj.alignObjectB_.seqBase_.outPutSeqAnsi(std::cout);
					debug_alignerObj->alignCacheGlobal(refSeq, seq);
					debug_alignerObj->rearrangeObjsGlobal(*refSeq, seq);

					debug_alignerObj->alignObjectA_.seqBase_.outPutSeqAnsi(std::cout);
					debug_alignerObj->alignObjectB_.seqBase_.outPutSeqAnsi(std::cout);
					std::cout << std::endl;
				}
			}
			writer.write(seq);
		}
	}

	if(!outputMaskingFnp.empty()) {
		OutputStream out(outOptions);
		outMaskInfo.outPutContents(out, "\t");
	}

	return 0;
}




}  // namespace njhseq
