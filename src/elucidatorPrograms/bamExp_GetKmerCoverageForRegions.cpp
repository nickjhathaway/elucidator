/*
 * bamExp_GetKmerCoverageForRegions.cpp
 *
 *  Created on: Jan 24, 2018
 *      Author: nick
 */


#include "bamExp.hpp"

#include "elucidator/BamToolsUtils.h"
#include "elucidator/objects/BioDataObject.h"



namespace njhseq {

int bamExpRunner::GetKmerCoverageForRegions(const njh::progutils::CmdArgs & inputCommands){
	bfs::path bedFnp = "";
	bfs::path twoBitFnp = "";
	uint32_t kmerLen = 6;
	OutOptions outOpts(bfs::path(""));
	seqSetUp setUp(inputCommands);

	setUp.processVerbose();
	setUp.processDebug();
	setUp.processWritingOptions(outOpts);
	setUp.setOption(kmerLen, "--kmerLen", "The size of the motif to count");
	setUp.setOption(bedFnp,  "--bed",     "Regions to collection coverage on");
	setUp.setOption(twoBitFnp, "--2bit", "2bit file of the genome mapped to", true);
	setUp.processReadInNames({"--bam"}, true);
	setUp.finishSetUp(std::cout);

	/**
	 * things to consider
	 * - should only plus strand direction be counted or both plus and negative strand
	 *   - going along with the above question should we just count of the actual direction of the query sequence be counted
	 * - should it be only mapped sequences or the full query sequence (here soft clipping could help avoid adapter or other artificial sequences)
	 * -
	 */

	std::unordered_map<std::string, uint64_t> kmerCoverageCount;

	TwoBit::TwoBitFile tReader(twoBitFnp);
	BamTools::BamReader bReader;
	bReader.Open(setUp.pars_.ioOptions_.firstName_.string());
	checkBamOpenThrow(bReader, setUp.pars_.ioOptions_.firstName_);
	bReader.LocateIndex();

	auto regions = gatherRegions(bedFnp.string(), "", setUp.pars_.verbose_);
	OutputStream out(outOpts);
	BamTools::BamAlignment bAln;
	for (const auto & reg : regions) {
		setBamFileRegionThrow(bReader, reg);
		while(bReader.GetNextAlignment(bAln)){
			if(bAln.IsMapped() && bAln.IsPrimaryAlignment()){
				std::string alignedNoGap = njh::replaceString(bAln.AlignedBases, "-", "");
				for(uint32_t pos : iter::range(alignedNoGap.size() + 1 - kmerLen)){
					++kmerCoverageCount[alignedNoGap.substr(pos, kmerLen)];
				}
			}
		}
	}

	std::unordered_map<std::string, uint64_t> genomeKmerCount;
	for(const auto & reg : regions){
		auto extractedSeq = reg.extractSeq(tReader);
		if(reg.reverseSrand_){
			extractedSeq.reverseComplementRead(false, true);
		}
		for(uint32_t pos : iter::range(extractedSeq.seq_.size() + 1 - kmerLen)){
			++genomeKmerCount[extractedSeq.seq_.substr(pos, kmerLen)];
		}
	}

	for(const auto & gK : genomeKmerCount){
		out << gK.first << "\t" << kmerCoverageCount[gK.first]/static_cast<double>(gK.second) << std::endl;
	}

	return 0;
}



}  // namespace njhseq

