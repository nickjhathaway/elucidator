/*
 * kmerExp_kmerPositionQualCounts.cpp
 *
 *  Created on: May 12, 2019
 *      Author: nicholashathaway
 */

#include "kmerExp.hpp"

#include <njhseq/IO/SeqIO/SeqIO.hpp>

namespace njhseq {






int kmerExpRunner::kmerPositionQualCounts(const njh::progutils::CmdArgs & inputCommands){
	OutOptions outOpts(bfs::path("out"), ".tab.txt");
	uint32_t klen = 35;
	uint32_t kCut = 1;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processReadInNames(true);
	setUp.processWritingOptions(outOpts);
	setUp.setOption(klen, "--klen", "kmer length");
	setUp.setOption(kCut, "--kCut", "k occurence cut off");

	setUp.finishSetUp(std::cout);

	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();

	seqInfo seq;
	OutputStream out(outOpts);
	std::unordered_map<std::string, uint32_t> kcounts;
	while (reader.readNextRead(seq)) {
		if (seq.seq_.size() > klen) {
			for (const auto pos : iter::range(seq.seq_.size() - klen + 1)) {
				++kcounts[seq.seq_.substr(pos, klen)];
			}
		}
	}
	reader.reOpenIn();
	std::unordered_map<std::string, std::vector<uint8_t>> kCountsQuals;

	while (reader.readNextRead(seq)) {
		if (seq.seq_.size() > klen) {
			for (const auto pos : iter::range(seq.seq_.size() - klen + 1)) {
				auto kmer = seq.seq_.substr(pos, klen) ;
				if( kcounts[kmer] > kCut){
					auto kmerQuals = getSubVector(seq.qual_, pos, klen);
					if(njh::in(kmer, kCountsQuals)){
						for(const auto qpos : iter::range(kmerQuals.size())){
							kCountsQuals[kmer][qpos]+= kmerQuals[qpos];
						}
					}else{
						kCountsQuals[kmer] = kmerQuals;
					}
				}
			}
		}
	}

	out << "kmer\tcount\tqpos\tavgQual" << std::endl;
	auto allKmers = getVectorOfMapKeys(kCountsQuals);
	njh::sort(allKmers);
	for (const auto & k : allKmers) {
		for (const auto pos : iter::range(k.size())) {
			out << k
					<< "\t" << kcounts[k]
					<< "\t"	<< pos
					<< "\t" << kCountsQuals[k][pos] / static_cast<double>(kcounts[k])
					<< std::endl;
		}
	}

	return 0;
}

} //namespace njhseq

