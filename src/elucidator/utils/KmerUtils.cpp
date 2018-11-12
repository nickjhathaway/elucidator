/*
 * KmerUtils.cpp
 *
 *  Created on: Apr 18, 2016
 *      Author: nick
 */


#include "KmerUtils.hpp"

namespace njhseq {

std::vector<std::unique_ptr<seqWithKmerInfo>> createKmerReadVec(const SeqIOOptions & opts,
		uint32_t kLength, bool setReverse){
	std::vector<std::unique_ptr<seqWithKmerInfo>> ret;
	seqInfo info;
	SeqInput reader(opts);
	reader.openIn();
	while(reader.readNextRead(info)){
		ret.emplace_back(std::make_unique<seqWithKmerInfo>(info, kLength, setReverse));
	}
	return ret;
}

void writeDistanceMatrix(std::ostream & out,
		const std::vector<std::vector<double>> & distances){
  for(const auto & row : distances){
  	out << njh::conToStr(row, "\t") << std::endl;;
  }
}

void writeDistanceMatrix(std::ostream & out,
		const std::vector<std::vector<double>> & distances,
		const VecStr & names) {
	printVector(names, "\t", out);
	for (const auto & rowPos : iter::range(distances.size())) {
		out << names[rowPos];
		for (const auto & colPos : iter::range(distances.size())){
			if(rowPos > colPos){
				out << "\t" << distances[rowPos][colPos];
			}else if (rowPos == colPos){
				out << "\t" << 0.00;
			}else{
				out << "\t" << distances[colPos][rowPos];
			}
		}
		out << std::endl;
	}
}

std::vector<std::vector<double>> readDistanceMatrix(std::istream & in){
	std::vector<std::vector<double>> inDist;
  for(std::string line; std::getline(in, line);){
  	auto toks = tokenizeString(line, "\t");
  	if(toks.empty() || toks.front() == ""){
  		inDist.emplace_back(std::vector<double>{});
  	}else{
  		inDist.emplace_back(njh::lexical_cast_con<std::vector<std::string>, std::vector<double>>(toks));
  	}
  }
  return inDist;
}


table getKmerStatsOnFile(const SeqIOOptions & opts,
		const seqWithKmerInfo & compare) {
	table outStats(VecStr{ "file", "compare", "kmerLen", "stats", "percentOrCount", "value" });
	auto kLength = compare.kInfo_.kLen_;
	auto reads = createKmerReadVec(opts, kLength, false);
	std::vector<double> percentShared;
	std::vector<uint32_t> numberShared;
	for (const auto & read : reads) {
		auto diff = compare.compareKmers(*read);
		percentShared.emplace_back(diff.second);
		numberShared.emplace_back(diff.first);
	}
	auto sharedStats = getStatsOnVec(numberShared);
	auto percentStats = getStatsOnVec(percentShared);
	for (const auto & shared : sharedStats) {
		outStats.content_.emplace_back(
				toVecStr(opts.firstName_, compare.seqBase_.name_, kLength, shared.first,
						"count", shared.second));
	}
	for (const auto & dist : percentStats) {
		outStats.content_.emplace_back(
				toVecStr(opts.firstName_, compare.seqBase_.name_, kLength, dist.first,
						"percent", dist.second));
	}
	return outStats;
}


}  // namespace njhseq
