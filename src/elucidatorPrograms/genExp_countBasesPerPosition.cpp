/*
 * genExp_countBasesPerPosition.cpp
 *
 *  Created on: Aug 2, 2019
 *      Author: nicholashathaway
 */



#include "genExp.hpp"

namespace njhseq {


int genExpRunner::countBasesPerPosition(const njh::progutils::CmdArgs & inputCommands){
	OutOptions outOpts(bfs::path(""), ".tab.txt");
	std::string name = "";
	//std::set<char> masterAlphabet{'A', 'G', 'T', 'C'};
	std::set<char> masterAlphabet{};
	seqSetUp setUp(inputCommands);
	setUp.processReadInNames(true);
	setUp.processWritingOptions(outOpts);
	setUp.setOption(name, "--name", "Add an optional name column");
	setUp.setOption(masterAlphabet, "--masterAlphabet", "Report counts even if zero of these bases");
	setUp.finishSetUp(std::cout);

	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	seqInfo seq;
	std::vector<charCounter> countsPerPosition;
	OutputStream out(outOpts);
	std::vector<uint32_t> lens;
	while(reader.readNextRead(seq)){
		lens.emplace_back(len(seq));
		for(const auto basePos : iter::range(seq.seq_.size())){
			if(countsPerPosition.size() <= basePos){
				countsPerPosition.emplace_back(charCounter{});
			}
			countsPerPosition[basePos].increaseCountOfBase(seq.seq_[basePos]);
		}
	}

	for( auto & count : countsPerPosition){
		count.resetAlphabet(false);
		count.setFractions();
		njh::addVecToSet(count.alphabet_, masterAlphabet);
	}
	if("" != name){
		out << "name\t";
	}
	out << "position\tbase\tcount\tfraction\tmeanReadLengthRounded" << std::endl;
	uint32_t meanLen = std::round(vectorMean(lens));
	for(const auto countPos : iter::range(countsPerPosition.size())){
		const auto & count = countsPerPosition[countPos];

		for(const auto base : masterAlphabet){
			if("" != name){
				out << name << '\t';
			}
			out << countPos << "\t";
			out << base
					<< "\t" << count.chars_[base]
					<< "\t" << count.fractions_[base]
					<< "\t" << meanLen << std::endl;
		}
	}

	return 0;
}

}  // namespace njhseq
