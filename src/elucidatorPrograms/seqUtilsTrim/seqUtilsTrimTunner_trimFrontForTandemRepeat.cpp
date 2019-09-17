/*
 * seqUtilsTrimTunner_trimFrontForTandemRepeat.cpp
 *
 *  Created on: Sep 16, 2019
 *      Author: nicholashathaway
 */




#include "seqUtilsTrimRunner.hpp"


namespace njhseq {


int seqUtilsTrimRunner::trimFrontForTandemRepeat(const njh::progutils::CmdArgs & inputCommands){
	std::string tandem = "";
	bool mark = false;
	seqUtilsTrimSetUp setUp(inputCommands);
	setUp.processDefaultReader();
	if ("" == setUp.pars_.ioOptions_.out_.outFilename_) {
		setUp.pars_.ioOptions_.out_.outFilename_ = njh::files::prependFileBasename(
				setUp.pars_.ioOptions_.firstName_, "trimmed_");
	} else if ("STDIN" == setUp.pars_.ioOptions_.firstName_) {
		setUp.pars_.ioOptions_.out_.outFilename_ = "STDOUT";
	}
	setUp.setOption(tandem, "--tandem", "tandem to trim front for");
	setUp.setOption(mark, "--mark", "mark seq with positions");
	setUp.processVerbose();
	setUp.processDebug();
	setUp.finishSetUp(std::cout);

	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	reader.openOut();

	seqInfo seq;

	uint32_t count = 0;
	while(reader.readNextRead(seq)){
		++count;
		std::cout << "count: " << count << std::endl;
		if(len(seq) > tandem.size()){
			const std::string front = seq.seq_.substr(0, tandem.size());

			if(!checkTwoRotatingStrings(tandem, front, 0).empty()){
				std::cout << __FILE__ << " " << __LINE__ << std::endl;
				uint32_t pos = tandem.size();
				while(pos + tandem.size() < seq.seq_.size() && front == seq.seq_.substr(pos, tandem.size())){
					pos += tandem.size();
				}
				std::cout << "pos: " << pos << std::endl;
				seq = seq.getSubRead(pos);
			}
		}
		reader.write(seq);
	}
	return 0;
}

} //namespace njhseq

