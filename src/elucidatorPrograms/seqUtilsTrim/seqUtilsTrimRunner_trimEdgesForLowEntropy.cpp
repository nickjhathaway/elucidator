/*
 * seqUtilsTrimRunner_trimEdgesForLowEntropy.cpp
 *
 *  Created on: Mar 2, 2020
 *      Author: nicholashathaway
 */




#include "seqUtilsTrimRunner.hpp"


namespace njhseq {





int seqUtilsTrimRunner::trimEdgesForLowEntropy(const njh::progutils::CmdArgs & inputCommands){
	uint32_t windowStep = 5;
	uint32_t windowSize = 50;
	uint32_t kLen = 1;
	double entropyCutOff = 1.5;
	bool mark = false;
	seqUtilsTrimSetUp setUp(inputCommands);
	setUp.description_ = "Trim front and back of sequences with a sliding window for low entropy, seqs will be removed if all is low entropy";

	FullTrimReadsPars pars;
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(mark, "--mark", "append the seq name with start and end position of the trim");

	setUp.setOption(windowStep, "--windowStep", "Window Step");
	setUp.setOption(windowSize, "--windowSize", "Window Size");
	setUp.setOption(kLen, "--kLen", "Kmer length for entropy calculation");
	setUp.setOption(entropyCutOff, "--entropyCutOff", "Entropy Cut Off");

	setUp.processIoOptions(setUp.readInFormatsAvailable_);
	setUp.finishSetUp(std::cout);

	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	seqInfo seq;

	while(reader.readNextRead(seq)){

		uint32_t start = 0;
		uint32_t end = len(seq);

		if(len(seq) >= windowSize){
			if(setUp.pars_.debug_){
				std::cout << seq.seq_ << std::endl;
			}
			for (auto pos : iter::range<uint32_t>(0, len(seq) - windowSize + 1, windowStep)) {
				if(setUp.pars_.debug_){
					std::cout << "\tpos:" << pos << " " << windowSize << std::endl;
					std::cout << "\tseq.seq_.size(): " << seq.seq_.size() << std::endl;
					std::cout << "\tpos + windowSize: " << pos + windowSize << std::endl;
					std::cout << "\tsubseq: " << seq.seq_.substr(pos, windowSize) << std::endl;
				}

				kmerInfo kInfo(seq.seq_.substr(pos, windowSize), kLen, false);
				if(kInfo.computeKmerEntropy() < entropyCutOff){
					start = pos + windowSize;
					if(setUp.pars_.debug_){
						std::cout << "\tkInfo.computeKmerEntropy(): " << kInfo.computeKmerEntropy() << std::endl;
					}
				}else{
					break;
				}
			}
			if(setUp.pars_.debug_){
				std::cout << seq.seq_ << std::endl;
			}

			for (auto pos : iter::range<uint32_t>(0, len(seq) - windowSize + 1, windowStep)) {
				if(setUp.pars_.debug_){
					std::cout << "\tpos:" << pos << " " << windowSize << std::endl;
					std::cout << "\tseq.seq_.size(): " << seq.seq_.size() << std::endl;
					std::cout << "\t" << "seq.seq_.size() - pos - windowSize: " << seq.seq_.size() - 1 - pos - windowSize  << std::endl;
					std::cout << "\tsubseq: " << seq.seq_.substr(seq.seq_.size() - pos - windowSize, windowSize) << std::endl;
				}

				kmerInfo kInfo(seq.seq_.substr(seq.seq_.size() - pos - windowSize, windowSize), kLen, false);
				if(kInfo.computeKmerEntropy() < entropyCutOff){
					end = seq.seq_.size() -pos - windowSize;
					if(setUp.pars_.debug_){
						std::cout << "\tkInfo.computeKmerEntropy(): " << kInfo.computeKmerEntropy() << std::endl;
					}
				}else{
					break;
				}
			}
		}else{
			kmerInfo kInfo(seq.seq_, kLen, false);
			if(kInfo.computeKmerEntropy() < entropyCutOff){
				end = 0;
			}
		}

		if(end > start){
			seq = seq.getSubRead(start, end - start);
			if(mark){
				MetaDataInName seqMeta;
				if(MetaDataInName::nameHasMetaData(seq.name_)){
					seqMeta= MetaDataInName(seq.name_);
				}
				seqMeta.addMeta("trimStart", start, true);
				seqMeta.addMeta("trimEnd", end, true);
				seqMeta.resetMetaInName(seq.name_);
			}
			reader.openWrite(seq);
		}
	}

	if(setUp.pars_.verbose_){
		setUp.logRunTime(std::cout);
	}
	return 0;
}

} // namespace njhseq


