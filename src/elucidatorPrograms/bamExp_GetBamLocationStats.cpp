/*
 * bamExp_GetBamLocationStats.cpp
 *
 *  Created on: Nov 25, 2017
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


#include "bamExp.hpp"
#include <TwoBit.h>

#include "elucidator/BamToolsUtils.h"
#include "elucidator/objects/BioDataObject.h"
#include "elucidator/objects/dataContainers/graphs.h"
#include <njhseq/GenomeUtils.h>



namespace njhseq {
int bamExpRunner::GetBamLocationStats(const njh::progutils::CmdArgs & inputCommands){
	bfs::path bedFnp = "";
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processDefaultReader({"--bam"}, true);
	setUp.setOption(bedFnp, "--bed", "Bed file for location", true);
	setUp.finishSetUp(std::cout);

	checkBamFilesForIndexesAndAbilityToOpen({setUp.pars_.ioOptions_.firstName_});
	BamTools::BamReader bReader;
	bReader.Open(setUp.pars_.ioOptions_.firstName_.string());
	checkBamOpenThrow(bReader, setUp.pars_.ioOptions_.firstName_);
	bReader.LocateIndex();
	BioDataFileIO <Bed3RecordCore> bedReader{IoOptions(InOptions(bedFnp))};
	bedReader.openIn();
	Bed3RecordCore record;
	if(!bedReader.readNextRecord(record)){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error should have at least 1 record in " << bedFnp << "\n";
		throw std::runtime_error{ss.str()};
	}
	GenomicRegion region(record);
	setBamFileRegionThrow(bReader, region);
	BamTools::BamAlignment bAln;
	uint32_t mapCount = 0;
	uint32_t seqsWithSoftClipping = 0;
	SeqOutput writer(setUp.pars_.ioOptions_);
	writer.openOut();
	while (bReader.GetNextAlignment(bAln)) {
		if (bAln.IsMapped()) {
			++mapCount;
			if (bAln.CigarData.front().Type == 'S'
					|| bAln.CigarData.back().Type == 'S') {
				++seqsWithSoftClipping;
				//CigarData
				auto jInfo = njh::json::toJson(bAln);
				std::cout << jInfo["CigarData"] << std::endl;
				auto seq = bamAlnToSeqInfo(bAln, true);
				writer.write(seq);
			}
		}
	}
	std::cout << seqsWithSoftClipping << "/" << mapCount << " "
			<< getPercentageString(seqsWithSoftClipping, mapCount) << std::endl;

	return 0;
}


} // namespace njhseq

