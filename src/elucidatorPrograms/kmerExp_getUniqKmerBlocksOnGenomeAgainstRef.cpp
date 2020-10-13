/*
 * kmerExp_getUniqKmerBlocksOnGenomeAgainstRef.cpp
 *
 *  Created on: May 30, 2020
 *      Author: nick
 */




#include "kmerExp.hpp"
#include "elucidator/utils/KmerUtils.hpp"
#include "elucidator/objects/seqObjects/seqKmers.h"
#include <njhseq/IO/SeqIO/SeqIO.hpp>


namespace njhseq {




int kmerExpRunner::getUniqKmerBlocksOnGenomeAgainstRef(const njh::progutils::CmdArgs & inputCommands){
	uint32_t kmerLength = 35;
	bfs::path refGenome = "";
	bfs::path compGenome = "";
	uint32_t numThreads = 1;
	bool trimNameAtWhiteSpace = false;
	OutOptions outOpts(bfs::path("out.bed"));
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(numThreads, "--numThreads", "Number of threads to use");
	setUp.processWritingOptions(outOpts);
	setUp.setOption(kmerLength, "--kmerLength", "kmer Length");
	setUp.setOption(refGenome, "--refGenome", "refGenome", true);
	setUp.setOption(compGenome, "--compGenome", "compGenome", true);

	setUp.setOption(trimNameAtWhiteSpace, "--trimNameAtWhiteSpace", "Trim Name At White Space");



	setUp.finishSetUp(std::cout);

	OutputStream out(outOpts);

	auto compGenomeOpts = SeqIOOptions(compGenome, SeqIOOptions::getInFormatFromFnp(compGenome));
	auto refGenomeOpts = SeqIOOptions(refGenome, SeqIOOptions::getInFormatFromFnp(refGenome));
	compGenomeOpts.includeWhiteSpaceInName_ = !trimNameAtWhiteSpace;


	std::unordered_set<std::string> refKmers;
	std::mutex refKmersMut;
	{

		SeqInput reader(refGenomeOpts);
		reader.openIn();

		std::function<void()> addKmers = [&reader,&kmerLength,&refKmers,&refKmersMut](){
			std::unordered_set<std::string> refKmersCurrent;
			seqInfo seq;
			while(reader.readNextReadLock(seq)){
				if(len(seq) > kmerLength){
					for(uint32_t pos = 0; pos < len(seq) -kmerLength + 1; ++pos){
						refKmersCurrent.emplace(seq.seq_.substr(pos, kmerLength));
					}
				}
			}
			{
				std::lock_guard<std::mutex> lock(refKmersMut);
				refKmers.insert(refKmersCurrent.begin(), refKmersCurrent.end());
			}
		};

		njh::concurrent::runVoidFunctionThreaded(addKmers, numThreads);
	}

	std::mutex compMut;
	SeqInput reader(compGenomeOpts);
	reader.openIn();

	std::function<void()> addKmers = [&reader,&kmerLength,&refKmers,&compMut,&out](){
		std::unordered_set<std::string> refKmersCurrent;
		seqInfo seq;
		while(reader.readNextReadLock(seq)){
			if(len(seq) > kmerLength){
				uint32_t start = 0;
				uint32_t run = 0;
				auto k = seq.seq_.substr(0, kmerLength);
				if(!njh::in(k, refKmers)){
					++run;
					start = 0;
				}
				for(uint32_t pos = 1; pos < len(seq) -kmerLength + 1; ++pos){
					auto k = seq.seq_.substr(pos, kmerLength);
					if(!njh::in(k, refKmers)){
						if(0 == run){
							start = pos;
						}
						++run;
					} else {
						if(run > 1){
							std::lock_guard<std::mutex> lock(compMut);
							auto end = start + run * kmerLength;
							out << seq.name_ << "\t"
									<< start << "\t"
									<< end << "\t"
									<< seq.name_ << "-" << start << "-" << end << "\t"
									<< end - start << "\t"
									<< '+' << std::endl;
						}
						run = 0;
					}
				}
				if(run > 1){
					std::lock_guard<std::mutex> lock(compMut);
					auto end = start + run * kmerLength;
					out << seq.name_ << "\t"
							<< start << "\t"
							<< end << "\t"
							<< seq.name_ << "-" << start << "-" << end << "\t"
							<< end - start << "\t"
							<< '+' << std::endl;
				}
			}
		}
	};
	njh::concurrent::runVoidFunctionThreaded(addKmers, numThreads);




	return 0;


}

} //namespace njhseq


