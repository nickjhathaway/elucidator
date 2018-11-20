/*
 * kmerExp_writeKmerSimDistanceMatrix.cpp
 *
 *  Created on: Feb 6, 2018
 *      Author: nick
 */
// elucidator - A library for analyzing sequence data
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
//
// This file is part of elucidator.
//
// elucidator is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// elucidator is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with elucidator.  If not, see <http://www.gnu.org/licenses/>.
//



#include "kmerExp.hpp"
#include "elucidator/utils/KmerUtils.hpp"
#include "elucidator/objects/seqObjects/seqKmers.h"

namespace njhseq {

int kmerExpRunner::writeKmerSimDistanceMatrix(const njh::progutils::CmdArgs & inputCommands){
	uint32_t kmerLength = 7;
	bool distance = false;
	bool writeNames = false;
	uint32_t numThreads = 1;
	OutOptions outOpts(bfs::path("out.tab.txt"));
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processReadInNames(true);
	setUp.processWritingOptions(outOpts);
	setUp.setOption(numThreads, "--numThreads", "numThreads");
	setUp.setOption(kmerLength, "--kmerLength", "kmer Length");
	setUp.setOption(distance, "--distance", "By default the the output is similarity score, with --distance it will be 1 - simScore to make the matrix a distance instead");
	setUp.setOption(writeNames, "--writeNames", "write out names");
	setUp.finishSetUp(std::cout);

  setUp.startARunLog(setUp.pars_.directoryName_);
  njh::stopWatch watch;
  watch.setLapName("Index Kmers");

  auto reads = createKmerReadVec(setUp.pars_.ioOptions_, kmerLength, false);
  OutputStream out(outOpts);

  std::function<
  				double(const std::unique_ptr<seqWithKmerInfo> &,
  						const std::unique_ptr<seqWithKmerInfo> &)> disFun =
  				[](const std::unique_ptr<seqWithKmerInfo> & read1,
  						const std::unique_ptr<seqWithKmerInfo> & read2) {
  					auto dist = read1->compareKmers(*read2);
  					return dist.second;
  				};
  watch.startNewLap("Distance");

  auto dist = getDistance(reads, numThreads, disFun);

  uint32_t rowCount = 0;
  if(writeNames){
  		printVector(readVec::getNames(reads), "\t", out);
  }
  if(distance){
		for (auto & row : dist) {
			njh::for_each(row, [](double & res){ res = 1 - res;});
			out << njh::conToStr(row, "\t");
			if(rowCount > 0){
				out << "\t";
			}
			out << 0 << std::endl;
			++rowCount;
		}
  }else{
		for (const auto & row : dist) {
			out << njh::conToStr(row, "\t");
			if(rowCount > 0){
				out << "\t";
			}
			out << 1 << std::endl;
			++rowCount;
		}
  }


  if(setUp.pars_.verbose_	){
  		watch.logLapTimes(std::cout, true, 6, true);
  }

	return 0;
}

} // namespace njhseq

