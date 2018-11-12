#pragma once
/*
 * KmerUtils.hpp
 *
 *  Created on: Apr 18, 2016
 *      Author: nick
 */

#include "elucidator/common.h"

namespace njhseq {

std::vector<std::unique_ptr<seqWithKmerInfo>> createKmerReadVec(
		const SeqIOOptions & opts, uint32_t kLength, bool setReverse);

std::vector<std::vector<double>> readDistanceMatrix(std::istream & in);

void writeDistanceMatrix(std::ostream & out,
		const std::vector<std::vector<double>> & distances);

void writeDistanceMatrix(std::ostream & out,
		const std::vector<std::vector<double>> & distances, const VecStr & names);

table getKmerStatsOnFile(const SeqIOOptions & seqFile,
		const seqWithKmerInfo & compare);

}  // namespace njhseq



