//
// Created by Nicholas Hathaway on 5/8/24.
//

#include "kmerSetExp.hpp"
#include "elucidator/objects/dataContainers.h"
#include "elucidator/simulation.h"
#include "elucidator/objects/seqObjects/seqKmers.h"

#include <njhseq/objects/seqObjects/seqKmers.h>
#include <njhseq/objects/kmer/KmerUtils.hpp>
#include <njhseq/concurrency/AllByAllPairFactory.hpp>


namespace njhseq {

kmerSetExpRunner::kmerSetExpRunner()
    : njh::progutils::ProgramRunner(
          {
	          addFunc("getKmerSetDistBetween", getKmerSetDistBetween, false),
	          addFunc("findUniqKmersBetweenSeqSetsMulti", findUniqKmersBetweenSeqSetsMulti, false),
	          addFunc("countingUniqKmersFromSets", countingUniqKmersFromSets, false),
	          addFunc("extractByCountingUniqKmersFromSets", extractByCountingUniqKmersFromSets, false),
	          addFunc("extractByCountingUniqKmersFromSetsIterative", extractByCountingUniqKmersFromSets, false),
	          addFunc("countingUniqKmersFromSetsPerRead", countingUniqKmersFromSetsPerRead, false),
	          addFunc("countingUniqKmersFromSetsInUnmappedAlns", countingUniqKmersFromSetsInUnmappedAlns, false),
	          addFunc("filterUniqueKmerSetForEntropy", filterUniqueKmerSetForEntropy, false),
	          addFunc("addToUniqKmersSet", addToUniqKmersSet, false),
	          addFunc("findUniqKmersFromGenomeSubRegionsMultiple", findUniqKmersFromGenomeSubRegionsMultiple, false),
	          addFunc("reportOnUniqKmersSet", reportOnUniqKmersSet, false),
	          addFunc("findKmersInSets", findKmersInSets, false),
	          addFunc("uniqKmersSetToFasta", uniqKmersSetToFasta, false),
	          addFunc("countingUniqKmersFromSetsBestSet", countingUniqKmersFromSetsBestSet, false),
	          addFunc("countingUniqKmersFromSetsInUnmappedAlnsBestSet", countingUniqKmersFromSetsInUnmappedAlnsBestSet,false),
	          addFunc("getUniqueSequenceRegions", getUniqueSequenceRegions, false),
	          addFunc("countingUniqKmersFromSetsInRegionsAlnsBestSet", countingUniqKmersFromSetsInRegionsAlnsBestSet,false),
          	addFunc("extractByKmerUniqueSets", extractByKmerUniqueSets, false),
          	addFunc("getBasicKmerSetCountInfo", getBasicKmerSetCountInfo, false),
          	addFunc("createMinimallyNonRedundantDownSampledSet", createMinimallyNonRedundantDownSampledSet, false),
          	addFunc("createKmerPresenceMatrixFromSets", createKmerPresenceMatrixFromSets, false),

//,
//
           },
          "kmerSetExp") {}




int kmerSetExpRunner::getKmerSetDistBetween(const njh::progutils::CmdArgs & inputCommands){
	bfs::path file1 = "";
	bfs::path file2 = "";
	std::string name1;
	std::string name2;
	uint32_t kLength = 5;
	bool getRevComp = false;
	OutOptions outOpts(bfs::path(""), ".tab.txt");
	seqSetUp setUp(inputCommands);
	setUp.setOption(file1, "--file1", "file1", true);
	setUp.setOption(file2, "--file2", "file2", true);
	setUp.setOption(name1, "--name1", "name1");
	setUp.setOption(name2, "--name2", "name2");
	setUp.processWritingOptions(outOpts);
	setUp.setOption(kLength, "--kLength", "Kmer Length");
	setUp.setOption(getRevComp, "--getRevComp", "getRevComp");

	setUp.finishSetUp(std::cout);

	SeqIOOptions file1Opts(file1, SeqIOOptions::getInFormatFromFnp(file1));
	SeqIOOptions file2Opts(file2, SeqIOOptions::getInFormatFromFnp(file2));
	if(name1.empty()){
		name1 = bfs::basename(file1);
	}
	if(name2.empty()){
		name2 = bfs::basename(file2);
	}
	OutputStream out(outOpts);
	seqInfo seq;

	kmerInfo file1Info;
	file1Info.kLen_ = kLength;

	SeqInput file1Reader(file1Opts);
	file1Reader.openIn();
	while(file1Reader.readNextRead(seq)){
		file1Info.updateKmers(seq.seq_, false);
	}

	kmerInfo file2Info;
	file2Info.kLen_ = kLength;

	SeqInput file2Reader(file2Opts);
	file2Reader.openIn();
	while(file2Reader.readNextRead(seq)){
		file2Info.updateKmers(seq.seq_, getRevComp);
	}

	std::pair<uint32_t, double> forDist = file1Info.compareKmers(file2Info);
	std::pair<uint32_t, double> revDist;
	if(getRevComp){
		revDist = file1Info.compareKmersRevComp(file2Info);
	}

	out << "set1\tset2\tkLength\tshared\tindex";
	if(getRevComp){
		out << "\tsharedRev\tindexRev";
	}
	out << std::endl;
	out << name1
			<< "\t" << name2
			<< "\t" << kLength
			<< "\t" << forDist.first
			<< "\t" << forDist.second;
	if(getRevComp){
		out << "\t" << revDist.first
				<< "\t" << revDist.second;
	}
	out << std::endl;

	return 0;
}







} // namespace njhseq



