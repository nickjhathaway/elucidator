/*
 * GenomeSeqSearch.hpp
 *
 *  Created on: Aug 26, 2017
 *      Author: nick
 */



#include <njhseq/objects/BioDataObject/GenomicRegion.hpp>
#include <njhseq/objects/kmer/kmerInfo.hpp>



namespace njhseq {

class GenomeSeqSearch {
public:
	class SeqStreak {
	public:
		SeqStreak(size_t seqPos, const GenomicRegion & region);

		size_t seqPos_;
		GenomicRegion region_;
		bool growing_ = true;

		bool isNextPos(const std::string & chrom, const size_t gPos,
				const bool isRevCompl, const size_t seqPos, const uint32_t kLen);

		void increaseEndPos(uint32_t size = 1);


	};

	GenomeSeqSearch(const OutOptions & outOpts);

	std::vector<SeqStreak> growingStreaks_;

	std::mutex mut_;
	std::shared_ptr<std::ofstream> out_;

	uint32_t countCapForGrowing_ = 10000;

	void growStreaks(const std::string & kmer, size_t seqPos, bool isRevComp,
			const std::unordered_map<std::string, std::unique_ptr<kmerInfo>> & infos);

	void pruneStreaks(uint32_t minSize);

	void purgeAllStreaks(uint32_t minSize);


};



}  // namespace njhseq
