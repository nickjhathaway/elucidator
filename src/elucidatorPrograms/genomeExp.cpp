/*
 * genomeExp.cpp
 *
 *  Created on: Dec 17, 2020
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

#include "genomeExp.hpp"
#include "elucidator/objects/BioDataObject.h"

#include <njhseq/objects/seqObjects/seqKmers.h>
#include <njhseq/concurrency/AllByAllPairFactory.hpp>

namespace njhseq {
genomeExpRunner::genomeExpRunner()
    : njh::progutils::ProgramRunner(
          {
					 addFunc("reorientToRefGenome", reorientToRefGenome, false),
           },
          "genomeExp") {}





int genomeExpRunner::reorientToRefGenome(const njh::progutils::CmdArgs & inputCommands) {
	//currently just turns phase into NA since it's hard to determine what rev complementing will do to the phase without knowing all genes
	bfs::path genomeGff;
	bfs::path genomeTwobit;
	bfs::path refGenomeTwobit;
  uint32_t numThreads = 1;
  uint32_t kLen = std::numeric_limits<uint32_t>::max();
  double cutOff = .75;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.setOption(genomeGff, "--gff", "Input gff file", true);
	setUp.setOption(genomeTwobit, "--genome2bit", "2bit genome file", true);
	setUp.setOption(refGenomeTwobit, "--refGenome2bit", "2bit ref genome file", true);
	setUp.setOption(kLen, "--kLen", "kmer Length Stop", true);
	setUp.setOption(cutOff, "--cutOff", "cutOff");
	setUp.setOption(numThreads, "--numThreads", "number of threads to use");
	setUp.processDirectoryOutputName(njh::pasteAsStr(bfs::path(bfs::basename(genomeTwobit)).replace_extension("").string(),
		"_reorientToRefGenome_TODAY"), true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	TwoBit::TwoBitFile queryReader(genomeTwobit);
	auto queryLens = queryReader.getSeqLens();
	TwoBit::TwoBitFile refReader(refGenomeTwobit);
	auto refLens = refReader.getSeqLens();

  njh::stopWatch watch;
  watch.setLapName("index kmers");
  std::vector<std::unique_ptr<seqWithKmerInfo>> querySeqs;
  std::string buffer;
  for(const auto & queryLen : queryLens){
  	queryReader[queryLen.first]->getSequence(buffer);
  	querySeqs.emplace_back(std::make_unique<seqWithKmerInfo>(seqInfo(queryLen.first, buffer)));
  }
  std::vector<std::unique_ptr<seqWithKmerInfo>> refSeqs;
  for(const auto & refLen : refLens){
  	refReader[refLen.first]->getSequence(buffer);
  	refSeqs.emplace_back(std::make_unique<seqWithKmerInfo>(seqInfo(refLen.first, buffer)));
  }

	auto queryNames = getVectorOfMapKeys(queryLens);
	njh::naturalSortNameSet(queryNames);

	if(setUp.pars_.verbose_){
		std::cout << njh::conToStr(queryNames, "\n") << std::endl;
	}
	allSetKmers(refSeqs, kLen, false);
	allSetKmers(querySeqs, kLen, true);

	struct KmerMatch{

		KmerMatch(const std::pair<uint32_t, double>  & dist, bool revComp):dist_(dist), revComp_(revComp){

		};
		std::pair<uint32_t, double> dist_;
		bool revComp_;
	};
	std::unordered_map<std::string, std::vector<KmerMatch>> matches;
	{
	  OutputStream outFile(njh::files::make_path(setUp.pars_.directoryName_, "kdistShared.tab.txt"));
	  outFile << "seq\tref\tkLen\tkDist\tkShared\n";
		std::mutex outFileMut;
		AllByAllPairFactory allFactory(querySeqs.size(), refSeqs.size());
		std::function<void()> getDist = [&allFactory,&refSeqs,&querySeqs,
										&outFileMut,&outFile,&kLen, &matches,&cutOff](){
			AllByAllPairFactory::AllByAllPair pair;
			while(allFactory.setNextPair(pair)){
				const auto & query = querySeqs[pair.row_];
				const auto & ref = refSeqs[pair.col_];
				auto dist = ref->compareKmers(*query);
				std::pair<uint32_t,double> revDist = ref->compareKmersRevComp(*query);
				{

					std::lock_guard<std::mutex> lock(outFileMut);
					outFile << query->seqBase_.name_
							<< "\t" << ref->seqBase_.name_
							<< "\t" << kLen
							<< "\t" << dist.second
							<< "\t" << dist.first << "\n";
	  			outFile << query->seqBase_.name_ << "_revComp"
	  					<< "\t" << ref->seqBase_.name_
							<< "\t" << kLen
							<< "\t" << revDist.second
	  					<< "\t" << revDist.first << "\n";
	  			if(revDist.second > cutOff){
	  				matches[query->seqBase_.name_].emplace_back(KmerMatch{revDist, true});
	  			}
	  			if(dist.second > cutOff){
	  				matches[query->seqBase_.name_].emplace_back(KmerMatch{dist, false});
	  			}
				}
			}
		};
		njh::concurrent::runVoidFunctionThreaded(getDist, numThreads);
	}
	VecStr seqsToComp;

	for(auto & match : matches){
		njh::sort(match.second, [](const KmerMatch & match1, const KmerMatch & match2){
			return match1.dist_.second > match2.dist_.second;
		});
		if(match.second.front().revComp_){
			seqsToComp.emplace_back(match.first);
		}
	}




	{
	  OutputStream seqsToComplementOutFile(njh::files::make_path(setUp.pars_.directoryName_, "seqsToComplement.txt"));
	  seqsToComplementOutFile <<  njh::conToStr(seqsToComp, "\n") << std::endl;
	}



	bfs::path outGffFnp (njh::files::make_path(setUp.pars_.directoryName_, bfs::basename(genomeGff) + ".gff" ));
	bfs::path outGenomeFnp(njh::files::make_path(setUp.pars_.directoryName_, bfs::basename(genomeGff) + ".fasta"));

	if(seqsToComp.empty()){

	}else{
		OutputStream outGenomeOut(outGenomeFnp);
		std::unordered_map<std::string, uint32_t> nameIndex;
		for(const auto  pos : iter::range(querySeqs.size())){
			nameIndex[querySeqs[pos]->seqBase_.name_] = pos;
		}
		//reverse complement the seqs
		for(const auto & queryName : queryNames){
			if(njh::in(queryName, seqsToComp)){
				querySeqs[nameIndex[queryName]]->seqBase_.reverseComplementRead(false, true);
				querySeqs[nameIndex[queryName]]->seqBase_.outPutSeq(outGenomeOut);
			}else{
				querySeqs[nameIndex[queryName]]->seqBase_.outPutSeq(outGenomeOut);
			}
		}


		//reverse complement the gff
		OutputStream outGffOut(outGffFnp);
		BioDataFileIO<GFFCore> reader{(IoOptions(InOptions(genomeGff)))};
		reader.openIn();
		uint32_t count = 0;
		std::string line = "";
		TwoBit::TwoBitFile tReader(genomeTwobit);
		auto chromLens = tReader.getSeqLens();

		{
			//write header
			std::ifstream infile(genomeGff.string());
			while('#' == infile.peek()){
				njh::files::crossPlatGetline(infile, line);
				outGffOut << line << std::endl;
			}
		}
		std::vector<std::shared_ptr<GFFCore>> allRecords;
		std::shared_ptr<GFFCore> gRecord = reader.readNextRecord();
		while(nullptr != gRecord) {
			if(njh::in(gRecord->seqid_, seqsToComp)){
				gRecord->revCompRecord(chromLens[gRecord->seqid_]);
			}
			gRecord->writeGffRecord(outGffOut);
			bool end = false;
			while ('#' == reader.inFile_->peek()) {
				if (njh::files::nextLineBeginsWith(*reader.inFile_, "##FASTA")) {
					//write out the fasta if there
					while('#' == reader.inFile_->peek()){
						njh::files::crossPlatGetline(*reader.inFile_, line);
						outGffOut << line << "\n";
					}
					end = true;
					break;
				}
				njh::files::crossPlatGetline(*reader.inFile_, line);
			}
			if (end) {
				break;
			}
			gRecord = reader.readNextRecord();
			++count;
		}
	}


	return 0;
}



} /* namespace njhseq */
