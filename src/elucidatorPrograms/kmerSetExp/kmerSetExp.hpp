#pragma once

//
// Created by Nicholas Hathaway on 5/8/24.
//

#include <njhcpp/progutils.h>


namespace njhseq {

class kmerSetExpRunner : public njh::progutils::ProgramRunner {
public:
	kmerSetExpRunner();

	static int getUniqueSequenceRegions(const njh::progutils::CmdArgs & inputCommands);
	static int findKmersInSets(const njh::progutils::CmdArgs & inputCommands);
	static int getKmerSetDistBetween(const njh::progutils::CmdArgs & inputCommands);

	static int findUniqKmersFromGenomeSubRegionsMultiple(const njh::progutils::CmdArgs & inputCommands);
	static int findUniqKmersBetweenSeqSetsMulti(const njh::progutils::CmdArgs & inputCommands);

	static int filterUniqueKmerSetForEntropy(const njh::progutils::CmdArgs & inputCommands);
	static int addToUniqKmersSet(const njh::progutils::CmdArgs & inputCommands);
	static int reportOnUniqKmersSet(const njh::progutils::CmdArgs & inputCommands);
	static int uniqKmersSetToFasta(const njh::progutils::CmdArgs & inputCommands);

	static int extractByKmerUniqueSets(const njh::progutils::CmdArgs & inputCommands);
	static int extractByCountingUniqKmersFromSets(const njh::progutils::CmdArgs & inputCommands);

	static int countingUniqKmersFromSetsPerRead(const njh::progutils::CmdArgs & inputCommands);
	static int countingUniqKmersFromSetsBestSet(const njh::progutils::CmdArgs & inputCommands);
	static int countingUniqKmersFromSetsInUnmappedAlns(const njh::progutils::CmdArgs & inputCommands);
	static int countingUniqKmersFromSetsInUnmappedAlnsBestSet(const njh::progutils::CmdArgs & inputCommands);
	static int countingUniqKmersFromSetsInRegionsAlnsBestSet(const njh::progutils::CmdArgs & inputCommands);
	static int countingUniqKmersFromSets(const njh::progutils::CmdArgs & inputCommands);

	static int getBasicKmerSetCountInfo(const njh::progutils::CmdArgs & inputCommands);
	static int createMinimallyNonRedundantDownSampledSet(const njh::progutils::CmdArgs & inputCommands);


};



} // namespace njhseq

