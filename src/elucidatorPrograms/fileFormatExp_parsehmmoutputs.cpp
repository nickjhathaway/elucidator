/*
 * fileFormatExp_parsehmmoutputs.cpp
 *
 *  Created on: Feb 19, 2022
 *      Author: nick
 */



#include "fileFormatExp.hpp"
#include "elucidator/objects/BioDataObject.h"
#include <njhseq/objects/BioDataObject/HmmerObjs/nhmmscanOutput.hpp>


namespace njhseq {





int fileFormatExpRunner::parsenhmmscanRaw(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path input = "";
	OutOptions outOpts(bfs::path(""), ".tab.txt");
	seqSetUp setUp(inputCommands);

	setUp.processDebug()	;
	setUp.processVerbose()	;
	setUp.setOption(input, "--file", "Input File Name", true);
	setUp.processWritingOptions(outOpts);

	setUp.finishSetUp(std::cout);

	OutputStream out(outOpts);

	auto ret = nhmmscanOutput::parseRawOutput(input);
	if(setUp.pars_.debug_){
		for(const auto & res : ret.qResults_){
			if(res.hits_.size() > 0){
				std::cout << res.toJson() << std::endl;
			}
		}
	}

	out << "query"
			<< "\t" << "queryLen"
			<< "\t" << "model"
			<< "\t" << "hmm_From"
			<< "\t" << "hmm_to"
			<< "\t" << "hmm_edges"
			<< "\t" << "model_len"
			<< "\t" << "aln_from"
			<< "\t" << "aln_to"
			<< "\t" << "aln_edge"
			<< "\t" << "env_from"
			<< "\t" << "env_to"
			<< "\t" << "env_edge"
			<< "\t" << "strand"
			<< "\t" << "evalue"
			<< "\t" << "score"
			<< "\t" << "bias"
			<< "\t" << "acc"
			<< "\t" << "scoreOverLen"
			<< "\t" << "percentGappedHit"
			<< "\t" << "percentPerfectHit"
			<< "\t" << "averagePP"
			<< std::endl;

	for(const auto & res : ret.qResults_){
		for(const auto & hit : res.hits_){
			out << res.queryName_
					<< "\t" << res.queryLen_
					<< "\t" << hit.targetName_
					<< "\t" << hit.hmmFrom_
					<< "\t" << hit.hmmTo_
					<< "\t" << hit.hmmEdgeInfo_
					<< "\t" << hit.modelLen_
					<< "\t" << hit.alignFrom_
					<< "\t" << hit.alignTo_
					<< "\t" << hit.aliEdgeInfo_
					<< "\t" << hit.envFrom_
					<< "\t" << hit.envTo_
					<< "\t" << hit.envEdgeInfo_

					<< "\t" << hit.strand_
					<< "\t" << hit.modelEvalue_
					<< "\t" << hit.modelScore_
					<< "\t" << hit.modelBias_
					<< "\t" << hit.acc_
					<< "\t" << hit.modelScore_/hit.envLen()
					<< "\t" << hit.percentGappedHit()
					<< "\t" << hit.percentPerfectHit()
					<< "\t" << hit.averagePP()
					<< std::endl;
		}
	}


	return 0;
}


}  // namespace njhseq
