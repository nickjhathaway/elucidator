/*
 * genExp_randomSampConsensusCompareToExpected.cpp
 *
 *  Created on: Dec 1, 2019
 *      Author: nicholashathaway
 */



#include "genExp.hpp"
#include <njhseq/IO/SeqIO/SeqIO.hpp>

#include <njhseq/objects/seqObjects/clusters/cluster.hpp>
#include <njhseq/readVectorManipulation/readVectorHelpers/readVecSorter.hpp>


namespace njhseq {


template<typename T>
cluster getConsensusForSeqs(std::vector<T> & reads, aligner & alignerObj,
		const std::string & name){
	if(reads.empty()){
		return cluster();
	}
	readVecSorter::sortReadVectorFunc<T>(reads, [](const T & first, const T & second) -> bool{
		if (roundDecPlaces(getSeqBase(first).cnt_, 2) == roundDecPlaces(getSeqBase(second).cnt_, 2) ) {
			if (roundDecPlaces(getRef(first).getAverageErrorRate(), 4)  < roundDecPlaces(getRef(second).getAverageErrorRate(), 4) ) {
				return true;
			} else {
				return false;
			}
		} else {
			return roundDecPlaces(getSeqBase(first).cnt_, 4)  > roundDecPlaces(getSeqBase(second).cnt_, 4) ;
		}
	}, true);


	cluster mainCluster(getSeqBase(*reads.begin()));
	if(reads.size() > 1){
		std::vector<cluster> inClusters;
		for(const auto readPos : iter::range<uint64_t>(1,reads.size())){
			inClusters.emplace_back(cluster(getSeqBase(reads[readPos])));
		}
		for(const auto & clus : inClusters){
			mainCluster.addRead(clus);
		}
	}
	mainCluster.calculateConsensus(alignerObj, true);
	//once a consensus has been built, build again against the built consensus
	mainCluster.needToCalculateConsensus_ = true;
	mainCluster.calculateConsensus(alignerObj, true);
	mainCluster.setName(name);
	return mainCluster;
}



int genExpRunner::randomSampConsensusCompareToExpected(const njh::progutils::CmdArgs & inputCommands){
	uint32_t sampleStart = 5;
	uint32_t sampleStep = 5;
	uint32_t sampleStop = 100;

	uint32_t samplings = 5;
	auto gapScoring = gapScoringParameters(5,1);
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processSeq(true);
	setUp.pars_.ioOptions_.lowerCaseBases_ = "remove";
  setUp.processReadInNames();
  setUp.processDirectoryOutputName(true);
  setUp.setOption(samplings,  "--samplings" ,"Number of times to sample");
  setUp.setOption(sampleStart,"--sampleStart" ,"sampleStart");
  setUp.setOption(sampleStep, "--sampleStep" , "sampleStep");
  setUp.setOption(sampleStop, "--sampleStop" , "sampleStop");

	setUp.finishSetUp(std::cout);

	setUp.startARunLog(setUp.pars_.directoryName_);

	uint64_t maxLen = 0;
	SeqInput reader(setUp.pars_.ioOptions_);
	auto allSeqs = reader.readAllReads<seqInfo>();
	readVec::getMaxLength(allSeqs, maxLen);
	readVec::getMaxLength(setUp.pars_.seqObj_, maxLen);
	maxLen = maxLen + 100;
	aligner alignerObj(maxLen, gapScoring, substituteMatrix::createDegenScoreMatrix(2, -2));
	std::vector<uint32_t> readPositions(allSeqs.size());
	njh::iota<uint32_t>(readPositions, 0);
	njh::randomGenerator rGen;
	OutputStream out(njh::files::make_path(setUp.pars_.directoryName_, "identities.tab.txt"));
	out << "samplingAmount\trun\tpercentIdentity" << std::endl;
	for(const auto sampNum : iter::range(sampleStart, sampleStop + sampleStep, sampleStep)){
		if(sampNum >= allSeqs.size()){
			break;
		}

		for(uint32_t s = 0; s < samplings; ++s){
			auto sampledReadPositions = rGen.unifRandSelectionVec(readPositions, sampNum, false);
			std::vector<seqInfo> currentSelection;
			for(const auto  pos : sampledReadPositions){
				currentSelection.emplace_back(allSeqs[pos]);
			}
			auto consensus = getConsensusForSeqs(currentSelection, alignerObj, njh::pasteAsStr(sampNum, "-", s));
			alignerObj.alignCacheGlobal(setUp.pars_.seqObj_, consensus);
			alignerObj.profilePrimerAlignment(setUp.pars_.seqObj_, consensus);
			out << sampNum << "\t" << s << "\t" << alignerObj.comp_.distances_.eventBasedIdentity_ << std::endl;
		}
	}

	return 0;

}

} // namespace njhseq


