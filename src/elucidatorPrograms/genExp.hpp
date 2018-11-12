#pragma once
/*
 * genExpRunner.h
 *
 *  Created on: Feb 4, 2015
 *      Author: nickhathaway
 */

#include <njhcpp/progutils.h>
#include <njhseq.h>

namespace njhseq {

class genExpRunner : public njh::progutils::ProgramRunner {
 public:

	genExpRunner();



	static int muscle(const njh::progutils::CmdArgs & inputCommands);
	static int extractByName(const njh::progutils::CmdArgs & inputCommands);

	static int createReadConsensus(const njh::progutils::CmdArgs & inputCommands);
	static int createReadConsensusRandomPicking(const njh::progutils::CmdArgs & inputCommands);

	static int findVariableSites(const njh::progutils::CmdArgs & inputCommands);


	static int tableCountToFasta(const njh::progutils::CmdArgs & inputCommands);
	static int tableCountToDistGraph(const njh::progutils::CmdArgs & inputCommands);


	static int getSnpInfo(const njh::progutils::CmdArgs & inputCommands);

	static int sortOnNucComp(const njh::progutils::CmdArgs & inputCommands);

	static int extractOnLen(const njh::progutils::CmdArgs & inputCommands);

	static int splitUpFile(const njh::progutils::CmdArgs & inputCommands);

	static int makeGraphFromAdjList(const njh::progutils::CmdArgs & inputCommands);


	static int evenOutCompReads(const njh::progutils::CmdArgs & inputCommands);

	static int extractRefGenesFastaFiles(const njh::progutils::CmdArgs & inputCommands);

	static int getMismatchDistances(const njh::progutils::CmdArgs & inputCommands);
	static int getIdDistances(const njh::progutils::CmdArgs & inputCommands);


	static int trimToMostProbableKmer(const njh::progutils::CmdArgs & inputCommands);
	static int trimFromMostProbableKmer(const njh::progutils::CmdArgs & inputCommands);

	static int trimReadNameAtFirstWhiteSpace(const njh::progutils::CmdArgs & inputCommands);

	static int mapReads(const njh::progutils::CmdArgs & inputCommands);


	static int splitFile(const njh::progutils::CmdArgs & inputCommands);
	static int genomePrimerExtract(const njh::progutils::CmdArgs & inputCommands);
	//

	static int findOutliersWithMuscleToRefs(const njh::progutils::CmdArgs & inputCommands);
	static int printMlnScores(const njh::progutils::CmdArgs & inputCommands);
	static int printMlnStreaks(const njh::progutils::CmdArgs & inputCommands);
	static int createAgreementMatrixWithMuscleRef(const njh::progutils::CmdArgs & inputCommands);
	static int createAgreementSegmentsWithMuscleRef(const njh::progutils::CmdArgs & inputCommands);





	static int concatenateDifferentLanes(const njh::progutils::CmdArgs & inputCommands);

	static int extractRefSeqsFromGenomes(const njh::progutils::CmdArgs & inputCommands);
	static int extractRefSeqsFromGenomesWithPrimers(const njh::progutils::CmdArgs & inputCommands);

	static int bioIndexGenomes(const njh::progutils::CmdArgs & inputCommands);
	static int bioIndexGenome(const njh::progutils::CmdArgs & inputCommands);

	static int getSeqFromFile(const njh::progutils::CmdArgs & inputCommands);


	static int bowtie2ExtractAndCompare(const njh::progutils::CmdArgs & inputCommands);
	static int bowtie2ExtractAndCompareMultiple(const njh::progutils::CmdArgs & inputCommands);

	static int lastzExtractAndCompare(const njh::progutils::CmdArgs & inputCommands);




	static int splitLibrary(const njh::progutils::CmdArgs & inputCommands);

};
} /* namespace njhseq */


