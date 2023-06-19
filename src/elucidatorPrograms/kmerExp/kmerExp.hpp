#pragma once
/*
 * kmerExpRunner.h
 *
 *  Created on: Feb 4, 2015
 *      Author: nickhathaway
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
#include <njhcpp/progutils.h>
#include <njhseq/common.h>
#include <njhseq/programUtils/seqSetUp.hpp>

namespace njhseq {

class kmerExpRunner : public njh::progutils::ProgramRunner {
 public:
	kmerExpRunner();

	static int getAvgKmerPositionPerSeq(const njh::progutils::CmdArgs & inputCommands);
	static int filterPerSeqAvgKmerPosition(const njh::progutils::CmdArgs & inputCommands);
	static int getAvgKmerPosition(const njh::progutils::CmdArgs & inputCommands);

	static int convertKmerSearchToBinaryMatrix(const njh::progutils::CmdArgs & inputCommands);
	static int kmerSearch(const njh::progutils::CmdArgs & inputCommands);
	static int microsatsKmerSearch(const njh::progutils::CmdArgs & inputCommands);



	static int scaningKmerDist(const njh::progutils::CmdArgs & inputCommands);
	static int profileScanningKmerDist(const njh::progutils::CmdArgs & inputCommands);
	static int getNewScanKmerDist(const njh::progutils::CmdArgs & inputCommands);
	static int profileKmerAccerlation(const njh::progutils::CmdArgs & inputCommands);
	static int readingDistanceCheck(const njh::progutils::CmdArgs & inputCommands);
	static int writingDistanceCheck(const njh::progutils::CmdArgs & inputCommands);
	static int writeKmerAccerlation(const njh::progutils::CmdArgs & inputCommands);
	static int writeKmerSimDistanceMatrix(const njh::progutils::CmdArgs & inputCommands);


	static int kDistVsNucDist(const njh::progutils::CmdArgs & inputCommands);
	static int scoveViaKmers(const njh::progutils::CmdArgs & inputCommands);
	static int pidVsKmers(const njh::progutils::CmdArgs & inputCommands);
	static int kmerRevVsForDist(const njh::progutils::CmdArgs & inputCommands);


	static int randomSampleKmerCompare(const njh::progutils::CmdArgs & inputCommands);
	static int profileLargeKmerIndex(const njh::progutils::CmdArgs & inputCommands);

	static int clostestKmerDist(const njh::progutils::CmdArgs & inputCommands);


	static int getKmerSetDistBetween(const njh::progutils::CmdArgs & inputCommands);


	static int getKmerDistTwoSeqs(const njh::progutils::CmdArgs & inputCommands);
	static int getKmerDistAgainstRef(const njh::progutils::CmdArgs & inputCommands);
	static int getKmerDistStatsMultiple(const njh::progutils::CmdArgs & inputCommands);
	static int getKmerDistStats(const njh::progutils::CmdArgs & inputCommands);
	static int getBestKmerDist(const njh::progutils::CmdArgs & inputCommands);

	static int getKmerDetailedKmerDistAgainstRef(const njh::progutils::CmdArgs & inputCommands);

	static int findingMinimumKLenForNoRedundantKmers(const njh::progutils::CmdArgs & inputCommands);





	static int kmerPositionQualCounts(const njh::progutils::CmdArgs & inputCommands);



	static int pairwiseWithinComparisonOfUniqueKmers(const njh::progutils::CmdArgs & inputCommands);
	static int allByAllComparisonOfUniqueKmers(const njh::progutils::CmdArgs & inputCommands);


	static int getKmerSharedBlocksBetweenGenomes(const njh::progutils::CmdArgs & inputCommands);
	static int getWithinGenomeUniqueKmers(const njh::progutils::CmdArgs & inputCommands);
	static int getUniqKmerBlocksOnGenomeAgainstRef(const njh::progutils::CmdArgs & inputCommands);
	static int findPositionsOfUniqueKmersInEachOther(const njh::progutils::CmdArgs & inputCommands);
	static int chromVsChromUniqueComp(const njh::progutils::CmdArgs & inputCommands);


	static int genomeKmerCompare(const njh::progutils::CmdArgs & inputCommands);
	static int profileSharedKmerBlocks(const njh::progutils::CmdArgs & inputCommands);
	static int kmerCompareTwoSetsOfContigs(const njh::progutils::CmdArgs & inputCommands);
	static int findUniqKmersBetweenSeqs(const njh::progutils::CmdArgs & inputCommands);

	//static int findUniqKmersBetweenSeqSets(const njh::progutils::CmdArgs & inputCommands);
	static int findUniqKmersBetweenSeqSetsMulti(const njh::progutils::CmdArgs & inputCommands);
	static int filterUniqueKmerSetForEntropy(const njh::progutils::CmdArgs & inputCommands);
	static int testingSimpleKmerHasher(const njh::progutils::CmdArgs & inputCommands);

	//static int findUniqKmersFromGenomeSubRegions(const njh::progutils::CmdArgs & inputCommands);
	static int findUniqKmersFromGenomeSubRegionsMultiple(const njh::progutils::CmdArgs & inputCommands);
	static int addToUniqKmersSet(const njh::progutils::CmdArgs & inputCommands);
	static int reportOnUniqKmersSet(const njh::progutils::CmdArgs & inputCommands);
	static int uniqKmersSetToFasta(const njh::progutils::CmdArgs & inputCommands);

	static int findKmersInSets(const njh::progutils::CmdArgs & inputCommands);
	static int countPerKmerPerSeq(const njh::progutils::CmdArgs & inputCommands );


	static int extractByCountingUniqKmersFromSets(const njh::progutils::CmdArgs & inputCommands);
	static int countingUniqKmersFromSetsPerRead(const njh::progutils::CmdArgs & inputCommands);
	static int countingUniqKmersFromSets(const njh::progutils::CmdArgs & inputCommands);
	static int countingUniqKmersFromSetsInUnmappedAlns(const njh::progutils::CmdArgs & inputCommands);

	static int kmerTestingGround(const njh::progutils::CmdArgs & inputCommands);


	static int kmerConnectionGraph(const njh::progutils::CmdArgs & inputCommands);


	static int countKmers(const njh::progutils::CmdArgs & inputCommands);
	static int generateCountsTable(const njh::progutils::CmdArgs & inputCommands);



	static int findingKmerEnrichment(const njh::progutils::CmdArgs & inputCommands);
	static int simpleHashKmer(const njh::progutils::CmdArgs & inputCommands);

};
} /* namespace njhseq */


