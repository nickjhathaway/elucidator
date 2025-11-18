//
// Created by Nicholas Hathaway on 11/17/25.
//

#include "seqUtilsModRunner.hpp"
#include <njhseq/IO/SeqIO/SeqIO.hpp>



namespace njhseq {

static constexpr std::array<char, 256> rc_table = []{
	std::array<char, 256> table{};
	table['A'] = 'T'; table['C'] = 'G'; table['G'] = 'C'; table['T'] = 'A'; table['N'] = 'N';
	table['a'] = 't'; table['c'] = 'g'; table['g'] = 'c'; table['t'] = 'a'; table['n'] = 'n';
	return table;
}();

inline void revcomp_inplace(std::string &s) {
	size_t i = 0, j = s.size() - 1;
	while (i < j) {
		const char ci = rc_table[static_cast<unsigned char>(s[i])];
		const char cj = rc_table[static_cast<unsigned char>(s[j])];
		s[i++] = cj;
		s[j--] = ci;
	}
	if (i == j) {
		s[i] = rc_table[static_cast<unsigned char>(s[i])];
	}
}



int seqUtilsModRunner::revCompFastqFast(const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	setUp.description_ = "reverse complement sequence from a fastq sequence, will only handle fastq and bases A,G,C,T,N";
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processDefaultReader(VecStr{"--fastq", "--fastqgz"}, true);
	setUp.finishSetUp(std::cout);

	OutputStream out(setUp.pars_.ioOptions_.out_);
	InputStream in(setUp.pars_.ioOptions_.firstName_);
	uint64_t line_count = 0;
	std::string line;
	while (njh::files::crossPlatGetline(in, line)) {
		switch (line_count % 4) {
			case 1:
				//reverse complement seq here
				revcomp_inplace(line);
				std::reverse(line.begin(), line.end());
			case 3:
				//just reverse the quality scores
				std::reverse(line.begin(), line.end());
				out << line << "\n";
			default:
				out << line << "\n";
		}
	}
	return 0;
}

} //namespace njhseq



