//
// Created by Nicholas Hathaway on 6/18/23.
//
#include "kmerExp.hpp"
#include "elucidator/objects/MiscUtility/GenomeSeqSearch.hpp"
#include "elucidator/objects/BioDataObject.h"
#include <njhseq/objects/Meta.h>
#include <njhseq/objects/kmer/SimpleKmerHash.hpp>
#include <njhseq/PopulationGenetics/PopGenCalcs.hpp>

namespace njhseq {


int kmerExpRunner::countPerKmerPerSeq(const njh::progutils::CmdArgs &inputCommands) {
	bfs::path seqsFnp;
	bool revComp = false;
	uint32_t kmerLength = 19;
	uint32_t kmerStep = 1;
	OutOptions outOpts(bfs::path(""));
	std::string elementStr;
	std::string columnName;

	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.pars_.ioOptions_.revComplMate_ = true;
	bool pairedEndSet = setUp.processReadInNames(njhseq::seqSetUp::pairedReadInFormatsAvailable_);
	auto singlesOption = SeqIOOptions();
	bool singlesEndSet = setUp.processJustReadInNames(singlesOption, njhseq::seqSetUp::singleInFormatsAvailable_,
																										!pairedEndSet);
	singlesOption.lowerCaseBases_ = setUp.pars_.ioOptions_.lowerCaseBases_;
	singlesOption.removeGaps_ = setUp.pars_.ioOptions_.removeGaps_;
	singlesOption.includeWhiteSpaceInName_ = setUp.pars_.ioOptions_.includeWhiteSpaceInName_;
	setUp.setOption(seqsFnp, "--seqsFnp", "seqs Fnp", true);

	setUp.setOption(revComp, "--revComp", "reverse complement");
	setUp.processWritingOptions(outOpts);
	setUp.setOption(columnName, "--newColumnName",
									"Name of a new Column to add to table, can be comma delimited to add multiple");
	setUp.setOption(elementStr, "--newColumnElement",
									"What to Add to the Table Under new Column, can be comma delimited when adding multiple new columns");

	setUp.finishSetUp(std::cout);

	SeqIOOptions seqsInOpts(seqsFnp, SeqIOOptions::getInFormatFromFnp(seqsFnp));
	njh::files::checkExistenceThrow(seqsFnp, __PRETTY_FUNCTION__ );

	std::unordered_map<std::string, uint64_t> kmerCounts;
	std::unordered_map<std::string, uint64_t> kmerCountsRevComp;
	VecStr newColumnEleToks;
	VecStr newColumnToks;
	if (!columnName.empty()) {
		newColumnEleToks = tokenizeString(elementStr, ",");
		newColumnToks = tokenizeString(columnName, ",");
		if (newColumnEleToks.size() != newColumnToks.size()) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__
				 << ", error, when adding multiple columns the number of elements need to match number of new columns" << "\n"
				 << "New Columns #: " << newColumnToks.size() << ", New Elements #: " << newColumnEleToks.size() << "\n";
			throw std::runtime_error{ss.str()};
		}
	}
	OutputStream out(outOpts);
	out << "name\tpos\tkmer\tisUnique\tcount";
	if(revComp){
		out << "\trevCompCount";
	}
	if (!newColumnToks.empty()) {
		for (const auto &col: newColumnToks) {
			out << "\t" << col;
		}
	}

	out << std::endl;

	//count kmers
	if(setUp.pars_.verbose_ && pairedEndSet && !setUp.pars_.ioOptions_.inExists()){
		std::cout << setUp.pars_.ioOptions_.firstName_ << " doesn't exist, skipping " << std::endl;
	}
	if(pairedEndSet && setUp.pars_.ioOptions_.inExists()){
		//paired
		SeqInput reader(setUp.pars_.ioOptions_);
		reader.openIn();
		PairedRead pseq;
		while(reader.readNextRead(pseq)){
			if(len(pseq.seqBase_) >= kmerLength){
				for(const auto pos : iter::range(pseq.seqBase_.seq_.size() + 1-kmerLength)){
					++kmerCounts[pseq.seqBase_.seq_.substr(pos, kmerLength)];
				}
				if(revComp){
					pseq.seqBase_.reverseComplementRead(false, true);
					for(const auto pos : iter::range(pseq.seqBase_.seq_.size() + 1-kmerLength)){
						++kmerCountsRevComp[pseq.seqBase_.seq_.substr(pos, kmerLength)];
					}
				}
			}
			if(len(pseq.mateSeqBase_) >= kmerLength){
				for(const auto pos : iter::range(pseq.mateSeqBase_.seq_.size() + 1-kmerLength)){
					++kmerCounts[pseq.mateSeqBase_.seq_.substr(pos, kmerLength)];
				}
				if(revComp){
					pseq.mateSeqBase_.reverseComplementRead(false, true);
					for(const auto pos : iter::range(pseq.mateSeqBase_.seq_.size() + 1-kmerLength)){
						++kmerCountsRevComp[pseq.mateSeqBase_.seq_.substr(pos, kmerLength)];
					}
				}
			}
		}
	}
	if(setUp.pars_.verbose_ && singlesEndSet && singlesOption.inExists()){
		std::cout << singlesOption.firstName_ << " doesn't exist, skipping " << std::endl;
	}
	if(singlesEndSet && singlesOption.inExists()){
		//single
		SeqInput reader(singlesOption);
		reader.openIn();
		seqInfo seq;
		while(reader.readNextRead(seq)){
			if(len(seq) >= kmerLength){
				for(const auto pos : iter::range(seq.seq_.size() + 1-kmerLength)){
					++kmerCounts[seq.seq_.substr(pos, kmerLength)];
				}
				if(revComp){
					seq.reverseComplementRead(false, true);
					for(const auto pos : iter::range(seq.seq_.size() + 1-kmerLength)){
						++kmerCountsRevComp[seq.seq_.substr(pos, kmerLength)];
					}
				}
			}
		}
	}
	std::unordered_map<std::string, uint32_t> kmerCountsForSeqs;

	{
		seqInfo seq;
		SeqInput reader(seqsInOpts);
		reader.openIn();
		while(reader.readNextRead(seq)){
			if(len(seq)){
				for(const auto pos : iter::range<uint32_t>(0, seq.seq_.size() + 1 - kmerLength, kmerStep)){
					auto k = seq.seq_.substr(pos, kmerLength);
					++kmerCountsForSeqs[k];
				}
				if(revComp){
					seq.reverseComplementRead(false, true);
					for(const auto pos : iter::range<uint32_t>(0, seq.seq_.size() + 1 - kmerLength, kmerStep)){
						auto k = seq.seq_.substr(pos, kmerLength);
						++kmerCountsForSeqs[k];
					}
				}
			}
		}
	}
	seqInfo seq;
	SeqInput reader(seqsInOpts);
	reader.openIn();
	while(reader.readNextRead(seq)){
		if(len(seq)){
			for(const auto pos : iter::range<uint32_t>(0, seq.seq_.size() + 1 - kmerLength, kmerStep)){
				auto k = seq.seq_.substr(pos, kmerLength);
				out << seq.name_
						<< "\t" << pos
						<< "\t" << k
						<< "\t" << njh::boolToStr(1 == kmerCountsForSeqs[k])
						<< "\t" << kmerCounts[k];
				if (revComp) {
					out << "\t" << kmerCountsRevComp[k];
				}
				if (!newColumnToks.empty()) {
					for (const auto &ele: newColumnEleToks) {
						out << "\t" << ele;
					}
				}
				out << std::endl;
			}
		}
	}

	return 0;
}

} // namespace njhseq
