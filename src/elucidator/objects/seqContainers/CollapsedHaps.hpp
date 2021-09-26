#pragma once

/*
 * CollapsedHaps.hpp
 *
 *  Created on: Jun 12, 2021
 *      Author: nick
 */


#include "elucidator/common.h"
#include "elucidator/objects/BioDataObject.h"
#include "elucidator/PopulationGenetics.h"

#include <njhseq/concurrency/pools/AlignerPool.hpp>

#include <njhseq/objects/Meta/MultipleGroupMetaData.hpp>

namespace njhseq {



class CollapsedHaps{
public:

	std::vector<std::shared_ptr<seqInfo>> seqs_;
	std::vector<std::unordered_set<std::string>> names_;
	std::vector<std::string> possibleSampleMetaFields_{"sample", "BiologicalSample"};

	std::unordered_map<std::string, uint32_t> subNamesToMainSeqPos_; /**< The position of the sub seqs in the collapsed unique seq vector */

	bool verbose_{false};

	//post processing
	void setSubNamesToMainSeqPos();
	void revCompSeqs();
	void setFrequencies(uint32_t total);
	void setFrequencies();

	std::unordered_map<std::string, uint32_t> renameBaseOnFreq(
			const std::string &identifier);

	////getting info
	uint32_t getTotalHapCount() const; /**< The total number of input haplotypes */
	uint32_t getTotalUniqueHapCount() const; /**< the total number of unique haplotypes */
	size_t size() const;

	std::unordered_map<std::string, uint32_t> genSeqNameKey() const;

	struct AvgPairwiseMeasures{
		double avgPercentId {0};
		double avgNumOfDiffs {0};
	};

	//population genetics
	struct GenPopMeasuresPar{
		bool getPairwiseComps {false};
		bool diagAlnPairwiseComps {true};

		uint32_t numSegSites_{std::numeric_limits<uint32_t>::max()};
		uint32_t numThreads = 1;
		double lowVarFreq = 0;
		VecStr genHeader() const;
	};

	struct GenPopMeasuresRes {
		PopGenCalculator::DiversityMeasures divMeasures_;
		PopGenCalculator::TajimaTestRes tajimaRes_;
		AvgPairwiseMeasures avgPMeasures_;
		std::vector<std::vector<comparison>> allComps_;

		VecStr getOut(const CollapsedHaps &inputSeqs, const std::string &identifier,
				const GenPopMeasuresPar &pars) const;

		void writeDivMeasures(const OutOptions &outOpts,
				const CollapsedHaps &inputSeqs, const std::string &identifier,
				const GenPopMeasuresPar &pars) const;
	};
	GenPopMeasuresRes getGeneralMeasuresOfDiversity(const GenPopMeasuresPar &pars,
			const std::shared_ptr<aligner> &alignerObj = nullptr) ;


	// getting reads lengths
	std::vector<uint32_t> getReadLenVec() const;
	std::unordered_map<uint32_t, uint32_t> getReadLenMap() const;
	bool hasLengthVariation(const double freqCutOff = 0) const;
	uint32_t getLongestLenDiff(const double freqCutOff = 0) const;


	std::vector<uint32_t> getOrder(const std::function<bool(const seqInfo &,const seqInfo&)> & comparator) const;
	std::vector<uint32_t> getOrderByTopCnt() const;

	//sample names
	static std::string getSampleNameFromSeqName(const std::string & name, const std::vector<std::string> & possibleSampleMetaFields=VecStr{"sample", "BiologicalSample"});
	static std::set<std::string> getPossibleLabIsolateNames(const std::unordered_set<std::string> & names);

	std::set<std::string> getAllSampleNames();
	std::vector<std::unordered_set<std::string>> getSampleNamesPerSeqs();

	//writing out info
	void writeOutSeqsOrdCnt(const SeqIOOptions &seqOpts) const;
	void writeNames(const OutOptions &outOpts) const;
	void writeLabIsolateNames(const OutOptions &outOpts) const;

	void writeOutMetaFields(const OutOptions &outOpts) const;


	//comparisons
	struct CompWithAlnSeqs {
		CompWithAlnSeqs();
		CompWithAlnSeqs(const comparison &comp, const std::string &refAln,
				const std::string &queryAln);
		comparison comp_;
		std::string refAlnSeq_;
		std::string queryAlnSeq_;
	};
	std::vector<CompWithAlnSeqs> getCompsAgainstRef(const seqInfo & refSeq, aligner & alignerObj, uint32_t numThreads = 1) const;
	std::vector<std::vector<comparison>> getPairwiseComps(aligner & alignerObj, uint32_t numThreads = 1) const;
	std::vector<std::vector<comparison>> getPairwiseCompsDiagAln(aligner & alignerObj, uint32_t numThreads = 1) const;


	AvgPairwiseMeasures getAvgPairwiseMeasures(const std::vector<std::vector<comparison>> & allComps) const;


	//factories
	static CollapsedHaps readInReads(const SeqIOOptions & inOpts,
			const std::unique_ptr<MultipleGroupMetaData> & meta = nullptr,
			const std::unordered_map<std::string, std::set<std::string>> & metaValuesToAvoid = std::unordered_map<std::string, std::set<std::string>>{});

	template<typename SEQTYPE>
	static CollapsedHaps collapseReads(const std::vector<SEQTYPE> & seqs,
			const std::unique_ptr<MultipleGroupMetaData> & meta = nullptr,
			const std::unordered_map<std::string, std::set<std::string>> & metaValuesToAvoid = std::unordered_map<std::string, std::set<std::string>>{}){


		CollapsedHaps ret;
		uint32_t seqCount = 0;
		std::unordered_set<std::string> allNames;

		for(const auto & seqObj : seqs) {
			seqInfo seq = getSeqBase(seqObj);
			if(nullptr != meta) {
				meta->attemptToAddSeqMeta(seq);
			}
			//get meta keys if available
			if(MetaDataInName::nameHasMetaData(getSeqBase(seq).name_)){
				MetaDataInName metaData(getSeqBase(seq).name_);
				bool skip = false;
				for(const auto & ignoreField : metaValuesToAvoid){
					if(metaData.containsMeta(ignoreField.first) && njh::in(metaData.getMeta(ignoreField.first), ignoreField.second)){
						skip = true;
						break;
						//skip this seq
					}
				}
				if(skip){
					continue;
				}
			}
			if(njh::in(seq.name_, allNames)){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "can't have seqs with the same name in input"<< "\n";
				ss << seq.name_ << " found more than once" << "\n";
				throw std::runtime_error{ss.str()};
			}
			allNames.emplace(seq.name_);

			seqCount+= std::round(seq.cnt_);
			bool found = false;
			for (const auto  pos : iter::range(ret.seqs_.size())) {
				const auto & otherSeq = ret.seqs_[pos];
				if (otherSeq->seq_ == seq.seq_) {
					otherSeq->cnt_ += 1;
					ret.names_[pos].emplace(seq.name_);
					found = true;
					break;
				}
			}
			if (!found) {
				//since we are considering counts here as haplotypes, will make sure cnt_ is set to 1
				seq.cnt_ = 1;
				ret.seqs_.emplace_back(std::make_shared<seqInfo>(seq));
				ret.names_.emplace_back(std::unordered_set<std::string>{seq.name_});
			}
		}
		ret.setFrequencies();
		ret.setSubNamesToMainSeqPos();
		return ret;

	}

	template<typename SEQTYPE>
	static CollapsedHaps collapseReads(const std::vector<SEQTYPE> & seqs,
			const std::vector<std::unordered_set<std::string>> & names
	){
		if(seqs.size() != names.size()){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "seqs and names must be same length; seqs.size(): " << seqs.size() << ", names.size():" << names.size()<< "\n";
			throw std::runtime_error{ss.str()};
		}

		CollapsedHaps ret;
		uint32_t seqCount = 0;
		std::unordered_set<std::string> allNames;
		uint32_t seqIndex = 0;
		for(const auto & seqObj : seqs) {
			seqInfo seq = getSeqBase(seqObj);
			//get meta keys if available
			if(njh::in(seq.name_, allNames)){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "can't have seqs with the same name in input"<< "\n";
				ss << seq.name_ << " found more than once" << "\n";
				throw std::runtime_error{ss.str()};
			}
			allNames.emplace(seq.name_);
			seqCount+= std::round(seq.cnt_);
			bool found = false;
			for (const auto  pos : iter::range(ret.seqs_.size())) {
				const auto & otherSeq = ret.seqs_[pos];
				if (otherSeq->seq_ == seq.seq_) {
					otherSeq->cnt_ += std::round(seq.cnt_);
					ret.names_[pos].insert(names[seqIndex].begin(), names[seqIndex].end());
					found = true;
					break;
				}
			}
			if (!found) {
				//since we are collapsing further already collapsed seqs, will add counts to total counts
				seq.cnt_ = std::round(seq.cnt_);
				ret.seqs_.emplace_back(std::make_shared<seqInfo>(seq));
				ret.names_.emplace_back(names[seqIndex]);
			}
			++seqIndex;
		}
		ret.setFrequencies();
		ret.setSubNamesToMainSeqPos();
		return ret;

	}


};




}  // namespace njhseq




