/*
 * seqUtilsModRunner_sorting.cpp
 *
 *  Created on: Oct 12, 2020
 *      Author: nick
 */




#include "seqUtilsModRunner.hpp"
#include "elucidator/objects/counters/DNABaseCounter.hpp"
#include <njhseq/IO/SeqIO/SeqIO.hpp>
#include <njhseq/objects/seqObjects/seqKmers.h>
#include <njhseq/readVectorManipulation/readVectorHelpers/readVecSorter.hpp>


namespace njhseq {
int seqUtilsModRunner::sortReadsByKmerEntropy(const njh::progutils::CmdArgs & inputCommands) {
	bool mark = false;
	uint32_t kLen = 2;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processDefaultReader(true);
  if (setUp.pars_.ioOptions_.out_.outFilename_ == "out") {
  	setUp.pars_.ioOptions_.out_.outFilename_ = njh::files::prependFileBasename(njh::files::removeExtension(setUp.pars_.ioOptions_.firstName_), "sorted_");
  }
	setUp.setOption(kLen, "--klen", "kmer length");
	setUp.setOption(mark, "--mark", "Add entropy to name");
	setUp.finishSetUp(std::cout);

	SeqInput reader(setUp.pars_.ioOptions_);
	auto inputSeqs = reader.readAllReads<seqInfo>();
	std::vector<seqWithKmerInfo> seqs;
	for(const auto & seq : inputSeqs){
		seqs.emplace_back(seqWithKmerInfo(seq, kLen, false));
	}
	if(mark){
		for(auto & seq : seqs){
			seq.seqBase_.name_.append(njh::pasteAsStr("[entropy=", seq.kInfo_.computeKmerEntropy(), "]"));
		}
	}
	njh::sort(seqs, []( const seqWithKmerInfo & seq1,  const seqWithKmerInfo & seq2){
		return seq1.kInfo_.computeKmerEntropy() < seq2.kInfo_.computeKmerEntropy();
	});

	SeqOutput::write(seqs, setUp.pars_.ioOptions_);

	return 0;
}

int seqUtilsModRunner::sortReadsByEntropy(const njh::progutils::CmdArgs & inputCommands) {
	bool mark = false;
	uint8_t minBase = 4;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processDefaultReader(true);
  if (setUp.pars_.ioOptions_.out_.outFilename_ == "out") {
  	setUp.pars_.ioOptions_.out_.outFilename_ = njh::files::prependFileBasename(njh::files::removeExtension(setUp.pars_.ioOptions_.firstName_), "sorted_");
  }

	setUp.setOption(minBase, "--minBase", "Min Base");
	setUp.setOption(mark, "--mark", "Add entropy to name");
	setUp.finishSetUp(std::cout);

	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	struct SeqInfoWithCountEntropy {
		SeqInfoWithCountEntropy(const seqInfo & seq, uint8_t minBase):seqBase_(seq){
			DNABaseCounter counter;
			counter.increase(seqBase_.seq_);
			entropy_ = counter.computeEntrophyBasedOffAlph(minBase);
		}
		seqInfo seqBase_;
		double entropy_;
	};

	std::vector<SeqInfoWithCountEntropy> seqs;
	seqInfo seq;
	while(reader.readNextRead(seq)){
		seqs.emplace_back(seq, minBase);
	}
	if(mark){
		for(auto & seq : seqs){
			seq.seqBase_.name_.append(njh::pasteAsStr("[entropy=", seq.entropy_, "]"));
		}
	}
	njh::sort(seqs, []( SeqInfoWithCountEntropy & seq1,  SeqInfoWithCountEntropy & seq2){
		return seq1.entropy_ < seq2.entropy_;
	});

	SeqOutput::write(seqs, setUp.pars_.ioOptions_);

	return 0;
}

int seqUtilsModRunner::sortReads(
		const njh::progutils::CmdArgs & inputCommands) {
	// parameters
	std::string sortBy = "averageError";
	bool decending = true;
	seqUtilsModSetUp setUp(inputCommands);
	setUp.setUpSortReads(sortBy, decending);
	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.in_.readAllReads<readObject>();
	if (stringToLowerReturn(sortBy) == "gc") {
		auto gcCompare = [](const readObject & read1, const readObject & read2) {
			if(read1.counter_.gcContent_ < read2.counter_.gcContent_) {
				return true;
			} else if (read1.counter_.gcContent_ == read2.counter_.gcContent_) {
				return read1.seqBase_.seq_ < read2.seqBase_.seq_;
			} else {
				return true;
			}
		};

		njh::for_each(inReads, [](readObject & read) {
			read.setLetterCount();
			read.counter_.calcGcContent();
		});
		njh::sort(inReads, gcCompare);
	} else {
		readVecSorter::sortReadVector(inReads, sortBy, decending);
	}
	reader.openWrite(inReads);
	setUp.logRunTime(std::cout);
	return 0;
}

int seqUtilsModRunner::sortReadsByKCompToTop(
		const njh::progutils::CmdArgs & inputCommands) {
	// parameters

	uint32_t klen = 4;
	bool decending = true;
	seqUtilsModSetUp setUp(inputCommands);
  // input file info
	setUp.processDefaultReader(true);
  if (setUp.pars_.ioOptions_.out_.outFilename_ == "out") {
  	setUp.pars_.ioOptions_.out_.outFilename_ = njh::files::prependFileBasename(njh::files::removeExtension(setUp.pars_.ioOptions_.firstName_), "sorted_");
  }
  setUp.setOption(klen, "--klen", "klen");
  bool ascending = false;
  setUp.setOption(ascending, "--ascending", "Ascending Sort");
  decending = !ascending;
  setUp.finishSetUp(std::cout);

	std::vector<std::shared_ptr<seqWithKmerInfo>> allInputReads;
	seqInfo seq;
	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	while (reader.readNextRead(seq)) {
		auto kSeq = std::make_shared<seqWithKmerInfo>(seq);
		readVec::handelLowerCaseBases(kSeq, setUp.pars_.ioOptions_.lowerCaseBases_);
		allInputReads.emplace_back(kSeq);
	}
	readVecSorter::sortByTotalCountAE(allInputReads, true);
	allSetKmers(allInputReads, klen, false);
	seqWithKmerInfo top = *allInputReads.front();
	readVecSorter::sortReadVectorFunc(allInputReads,[&top](const std::shared_ptr<seqWithKmerInfo>& p1, const std::shared_ptr<seqWithKmerInfo> & p2){
		return top.compareKmers(*p1).second > top.compareKmers(*p2).second;
	});


	reader.openWrite(allInputReads);
	setUp.logRunTime(std::cout);
	return 0;
}




} //namespace njhseq

