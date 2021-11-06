/*
 * kmerExp_getKmerDetailedKmerDistAgainstRef.cpp
 *
 *  Created on: Oct 30, 2021
 *      Author: nick
 */




#include "kmerExp.hpp"
#include "elucidator/utils/KmerUtils.hpp"
#include "elucidator/objects/dataContainers.h"
#include "elucidator/simulation.h"
#include "elucidator/objects/seqObjects/seqKmers.h"

#include "elucidator/objects/MiscUtility/GenomeSeqSearch.hpp"
#include "elucidator/objects/BioDataObject.h"

#include <njhseq/helpers.h>

namespace njhseq {



int kmerExpRunner::getKmerDetailedKmerDistAgainstRef(const njh::progutils::CmdArgs & inputCommands){
	bfs::path inputFnp = "";
	OutOptions outOpts(bfs::path(""), ".tab.txt");
	uint32_t numThreads = 1;
	uint32_t kmerLength = 7;

	bool getRevComp = false;


	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processReadInNames(true);
	setUp.processRefFilename(true);
	setUp.processWritingOptions(outOpts);

	setUp.setOption(numThreads, "--numThreads", "number of threads to use");
	setUp.setOption(kmerLength, "--kmerLength", "kmer Length");
	setUp.setOption(getRevComp, "--getRevComp", "Compare Reverse Complement as well");

	setUp.finishSetUp(std::cout);

	setUp.pars_.refIoOptions_.includeWhiteSpaceInName_ = setUp.pars_.ioOptions_.includeWhiteSpaceInName_;
	auto refSeqs = createKmerReadVec(setUp.pars_.refIoOptions_, kmerLength, getRevComp);

	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();

	OutputStream out(outOpts);
	out << "name\tref\trevComp\ttotalShared\tdist\tdistLenAdjust\ttotalKmersIn1\ttotalKmersIn2\ttotalUniqueShared\tuniqDist\tuniqDistLenAdjust\ttotalUniqKmersIn1\ttotalUniqKmersIn2\ttotalUniq" << "\n";

	std::mutex outMut;
	std::function<void()> getBestDist = [&reader,&out,&outMut,
																			 &getRevComp,&kmerLength,
																			 &refSeqs](){
		seqInfo seq;
		while(reader.readNextReadLock(seq)){

			kmerInfo kInfo(seq.seq_, kmerLength, false);
			for(const auto pos : iter::range(refSeqs.size())){
				const auto & refSeq = refSeqs[pos];
				auto dist = refSeq->kInfo_.compareKmersDetailed(kInfo);
				{
					std::lock_guard<std::mutex> lock(outMut);
					out << seq.name_
							<< "\t" << refSeq->seqBase_.name_
							<< "\t" << njh::boolToStr(false)
							<< "\t" << dist.totalShared_
							<< "\t" << dist.getDistTotalShared()
							<< "\t" << dist.getDistTotalSharedLenAdjusted()
							<< "\t" << dist.totalKmersIn1_
							<< "\t" << dist.totalKmersIn2_
							<< "\t" << dist.totalUniqShared_
							<< "\t" << dist.getDistUniqueShared()
							<< "\t" << dist.getDistUniqueSharedLenAdjusted()
							<< "\t" << dist.totalUniqKmersIn1_
							<< "\t" << dist.totalUniqKmersIn2_
							<< "\t" << dist.totalUniqBetween_
							<< "\n";
				}
				if(getRevComp){
					auto revDist = kInfo.compareKmersRevCompDetailed(refSeq->kInfo_);
					{
						std::lock_guard<std::mutex> lock(outMut);
						out << seq.name_
								<< "\t" << refSeq->seqBase_.name_
								<< "\t" << njh::boolToStr(true)
								<< "\t" << revDist.totalShared_
								<< "\t" << revDist.getDistTotalShared()
								<< "\t" << revDist.getDistTotalSharedLenAdjusted()
								<< "\t" << revDist.totalKmersIn1_
								<< "\t" << revDist.totalKmersIn2_
								<< "\t" << revDist.totalUniqShared_
								<< "\t" << revDist.getDistUniqueShared()
								<< "\t" << revDist.getDistUniqueSharedLenAdjusted()
								<< "\t" << revDist.totalUniqKmersIn1_
								<< "\t" << revDist.totalUniqKmersIn2_
								<< "\t" << revDist.totalUniqBetween_
								<< "\n";
					}
				}
			}
		}
	};

	njh::concurrent::runVoidFunctionThreaded(getBestDist, numThreads);



	return 0;
}


} // namespace njhseq



