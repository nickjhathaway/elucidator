//
// Created by Nicholas Hathaway on 4/1/25.
//


#include "seqUtilsModRunner.hpp"
#include <njhseq/objects/counters/DNABaseCounter.hpp>
#include <njhseq/IO/SeqIO/SeqIO.hpp>
#include <njhseq/objects/seqObjects/seqKmers.h>
#include <njhseq/readVectorManipulation/readVectorHelpers/readVecSorter.hpp>


namespace njhseq {

int seqUtilsModRunner::fragmentSequencesForPhipSeq(const njh::progutils::CmdArgs & inputCommands) {
	uint32_t step = 4;
	uint32_t window_size = 16;
	double back_seq_overlap_ratio = 0.75;
	bool doNotModifyName = false;
	seqSetUp setUp(inputCommands);
	setUp.description_ = "fragment input with special considerations for preparation of creation of a phipseq library";
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processDefaultReader(true);
	setUp.setOption(step, "--step", "step size");
	setUp.setOption(window_size, "--window_size", "window size");
	setUp.setOption(doNotModifyName, "--doNotModifyName", "do Not Modify output Name");
	setUp.setOption(back_seq_overlap_ratio, "--back_seq_overlap_ratio", "back seq overlap ratio to allow ");

	setUp.finishSetUp(std::cout);

	SeqIO seq_io(setUp.pars_.ioOptions_);
	seq_io.openIn();
	seq_io.openOut();

	seqInfo seq;
	int64_t expected_overlap = window_size - step;
	while (seq_io.readNextRead(seq)) {
		if (seq.seq_.size() > window_size + step) {
			uint32_t seq_count = 0;
			auto back_seq_pos = seq.seq_.size() - window_size;
			for (uint32_t pos = 0; pos + window_size < seq.seq_.size(); pos+=step) {
				auto end = pos + window_size;
				if (end >= back_seq_pos && setUp.pars_.debug_) {
					std::cout << "pos: " << pos << std::endl;
					std::cout << "back_seq_pos: " << back_seq_pos << std::endl;
					std::cout << "end: " << end  << std::endl;
					std::cout << "expected_overlap: " << expected_overlap << std::endl;
					std::cout << "end - back_seq_pos: " << end - back_seq_pos << std::endl;
					std::cout << "static_cast<long double>(end - back_seq_pos)/expected_overlap: " << expected_overlap/static_cast<long double>(end - back_seq_pos) << std::endl<< std::endl;
				}
				if (end < back_seq_pos || (end > back_seq_pos && expected_overlap/static_cast<long double>(end - back_seq_pos)  > back_seq_overlap_ratio)) {
					auto out_seq = seq.getSubRead(pos, window_size);
					out_seq.name_.append(njh::pasteAsStr("_seq", seq_count));
					seq_io.write(out_seq);
					++seq_count;
				}
			}
			//add back
			{
				auto out_seq = seq.getSubRead(seq.seq_.size() - window_size, window_size);
				out_seq.name_.append(njh::pasteAsStr("_seq", seq_count));
				seq_io.write(out_seq);
				++seq_count;
			}
		}
	}
	return 0;
}


} //namespace njhseq

