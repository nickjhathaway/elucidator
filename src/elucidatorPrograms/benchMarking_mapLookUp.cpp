/*
 * benchMarking_mapLookUp.cpp
 *
 *  Created on: Nov 23, 2020
 *      Author: nick
 */

#include "benchMarking.hpp"
#include <njhseq/IO.h>


namespace njhseq {


int benchMarkingRunner::mapLookUp(const njh::progutils::CmdArgs & inputCommands){
	uint32_t kmerLength = 31;
	OutOptions outOpts(bfs::path(""), ".tab.txt");
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processReadInNames({"--fastq", "--fastqgz"});
	setUp.setOption(kmerLength, "--kmerLen", "Kmer length to test");
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);



	if(setUp.pars_.debug_){
		std::vector<uint32_t> numbers(256);
		njh::iota<uint32_t>(numbers, 0);
		std::cout << 'A' << ": " << numbers['A'] << std::endl;
		std::cout << 'C' << ": " << numbers['C'] << std::endl;
		std::cout << 'G' << ": " << numbers['G'] << std::endl;
		std::cout << 'T' << ": " << numbers['T'] << std::endl;
		std::cout << 'a' << ": " << numbers['a'] << std::endl;
		std::cout << 'c' << ": " << numbers['c'] << std::endl;
		std::cout << 'g' << ": " << numbers['g'] << std::endl;
		std::cout << 't' << ": " << numbers['t'] << std::endl;
	}
	OutputStream out(outOpts);

	auto seqs = SeqInput::getSeqVec<seqInfo>(setUp.pars_.ioOptions_	);
	{
		std::unordered_map<std::string, uint32_t> kmerCounts;
		{
			njh::scopedStopWatch watch("unordered_map");

			for(const auto & seq : seqs){
				if(len(seq) > kmerLength){
					for(uint32_t pos = 0; pos < len(seq) - kmerLength + 1; ++pos){
						std::string k = seq.seq_.substr(pos, kmerLength);
						++kmerCounts[k];
					}
				}
			}
			out << "unordered_map"
					 << "\t" << kmerLength
					 << "\t" << 1
					 << "\t" << watch.totalTime() << std::endl;;
		}
	}


	{
		std::unordered_map<std::string, uint32_t> kmerCounts;
		{
			njh::scopedStopWatch watch("unordered_map");

			for(const auto & seq : seqs){
				if(len(seq) > kmerLength){
					for(uint32_t pos = 0; pos < len(seq) - kmerLength + 1; ++pos){
						std::string k = seq.seq_.substr(pos, kmerLength);
						++kmerCounts[k];
					}
				}
			}
			out << "unordered_map"
					 << "\t" << kmerLength
					 << "\t" << 2
					 << "\t" << watch.totalTime() << std::endl;
		}
	}



	//1 depth vector

	{
		std::vector<std::unordered_map<std::string, uint32_t>> kmerContsVector(117);
		{
			njh::scopedStopWatch watch("unordered_map_vector1");

			for(const auto & seq : seqs){
				if(len(seq) > kmerLength){
					for(uint32_t pos = 0; pos < len(seq) - kmerLength + 1; ++pos){
						std::string k = seq.seq_.substr(pos, kmerLength);
						++kmerContsVector[k[0]][k];
					}
				}
			}
			out << "unordered_map_vector1"
					 << "\t" << kmerLength
					 << "\t" << 1
					 << "\t" << watch.totalTime() << std::endl;;
		}
	}

	{
		std::vector<std::unordered_map<std::string, uint32_t>> kmerContsVector(117);
		{
			njh::scopedStopWatch watch("unordered_map_vector1");

			for(const auto & seq : seqs){
				if(len(seq) > kmerLength){
					for(uint32_t pos = 0; pos < len(seq) - kmerLength + 1; ++pos){
						std::string k = seq.seq_.substr(pos, kmerLength);
						++kmerContsVector[k[0]][k];
					}
				}
			}
			out << "unordered_map_vector1"
					 << "\t" << kmerLength
					 << "\t" << 2
					 << "\t" << watch.totalTime() << std::endl;;
		}
	}



	//2 depth vector
	{
		std::vector<std::vector<std::unordered_map<std::string, uint32_t>>> kmerContsVector(117, std::vector<std::unordered_map<std::string, uint32_t>>(117));
		{
			njh::scopedStopWatch watch("unordered_map_vector2");

			for(const auto & seq : seqs){
				if(len(seq) > kmerLength){
					for(uint32_t pos = 0; pos < len(seq) - kmerLength + 1; ++pos){
						std::string k = seq.seq_.substr(pos, kmerLength);
						++kmerContsVector[k[0]][k[1]][k];
					}
				}
			}
			out << "unordered_map_vector2"
					 << "\t" << kmerLength
					 << "\t" << 1
					 << "\t" << watch.totalTime() << std::endl;;
		}
	}

	{
		std::vector<std::vector<std::unordered_map<std::string, uint32_t>>> kmerContsVector(117, std::vector<std::unordered_map<std::string, uint32_t>>(117));
		{
			njh::scopedStopWatch watch("unordered_map_vector2");

			for(const auto & seq : seqs){
				if(len(seq) > kmerLength){
					for(uint32_t pos = 0; pos < len(seq) - kmerLength + 1; ++pos){
						std::string k = seq.seq_.substr(pos, kmerLength);
						++kmerContsVector[k[0]][k[1]][k];
					}
				}
			}
			out << "unordered_map_vector2"
					 << "\t" << kmerLength
					 << "\t" << 2
					 << "\t" << watch.totalTime() << std::endl;;
		}
	}


	//3 depth vector
	{
		std::vector<std::vector<std::vector<std::unordered_map<std::string, uint32_t>>>> kmerContsVector(117, std::vector<std::vector<std::unordered_map<std::string, uint32_t>>>(117, std::vector<std::unordered_map<std::string, uint32_t>>(117)));
		{
			njh::scopedStopWatch watch("unordered_map_vector3");

			for(const auto & seq : seqs){
				if(len(seq) > kmerLength){
					for(uint32_t pos = 0; pos < len(seq) - kmerLength + 1; ++pos){
						std::string k = seq.seq_.substr(pos, kmerLength);
						++kmerContsVector[k[0]][k[1]][k[2]][k];
					}
				}
			}
			out << "unordered_map_vector3"
					 << "\t" << kmerLength
					 << "\t" << 1
					 << "\t" << watch.totalTime() << std::endl;;
		}
	}

	{
		std::vector<std::vector<std::vector<std::unordered_map<std::string, uint32_t>>>> kmerContsVector(117, std::vector<std::vector<std::unordered_map<std::string, uint32_t>>>(117, std::vector<std::unordered_map<std::string, uint32_t>>(117)));
		{
			njh::scopedStopWatch watch("unordered_map_vector3");

			for(const auto & seq : seqs){
				if(len(seq) > kmerLength){
					for(uint32_t pos = 0; pos < len(seq) - kmerLength + 1; ++pos){
						std::string k = seq.seq_.substr(pos, kmerLength);
						++kmerContsVector[k[0]][k[1]][k[2]][k];
					}
				}
			}
			out << "unordered_map_vector3"
					 << "\t" << kmerLength
					 << "\t" << 2
					 << "\t" << watch.totalTime() << std::endl;;
		}
	}

	//3 depth vector
	{
		std::vector<std::vector<std::vector<std::unordered_map<std::string, uint32_t>>>> kmerContsVector(117, std::vector<std::vector<std::unordered_map<std::string, uint32_t>>>(117, std::vector<std::unordered_map<std::string, uint32_t>>(117)));
		{
			njh::scopedStopWatch watch("unordered_map_vector3");

			for(const auto & seq : seqs){
				if(len(seq) > kmerLength){
					for(uint32_t pos = 0; pos < len(seq) - kmerLength + 1; ++pos){
						std::string k = seq.seq_.substr(pos, kmerLength);
						++kmerContsVector[k[0]][k[1]][k[2]][k];
					}
				}
			}
			out << "unordered_map_vector3"
					 << "\t" << kmerLength
					 << "\t" << 3
					 << "\t" << watch.totalTime() << std::endl;;
		}
	}

	{
		std::vector<std::vector<std::vector<std::unordered_map<std::string, uint32_t>>>> kmerContsVector(117, std::vector<std::vector<std::unordered_map<std::string, uint32_t>>>(117, std::vector<std::unordered_map<std::string, uint32_t>>(117)));
		{
			njh::scopedStopWatch watch("unordered_map_vector3");

			for(const auto & seq : seqs){
				if(len(seq) > kmerLength){
					for(uint32_t pos = 0; pos < len(seq) - kmerLength + 1; ++pos){
						std::string k = seq.seq_.substr(pos, kmerLength);
						++kmerContsVector[k[0]][k[1]][k[2]][k];
					}
				}
			}
			out << "unordered_map_vector3"
					 << "\t" << kmerLength
					 << "\t" << 4
					 << "\t" << watch.totalTime() << std::endl;;
		}
	}

	//2 depth vector
	{
		std::vector<std::vector<std::unordered_map<std::string, uint32_t>>> kmerContsVector(117, std::vector<std::unordered_map<std::string, uint32_t>>(117));
		{
			njh::scopedStopWatch watch("unordered_map_vector2");

			for(const auto & seq : seqs){
				if(len(seq) > kmerLength){
					for(uint32_t pos = 0; pos < len(seq) - kmerLength + 1; ++pos){
						std::string k = seq.seq_.substr(pos, kmerLength);
						++kmerContsVector[k[0]][k[1]][k];
					}
				}
			}
			out << "unordered_map_vector2"
					 << "\t" << kmerLength
					 << "\t" << 3
					 << "\t" << watch.totalTime() << std::endl;;
		}
	}

	{
		std::vector<std::vector<std::unordered_map<std::string, uint32_t>>> kmerContsVector(117, std::vector<std::unordered_map<std::string, uint32_t>>(117));
		{
			njh::scopedStopWatch watch("unordered_map_vector2");

			for(const auto & seq : seqs){
				if(len(seq) > kmerLength){
					for(uint32_t pos = 0; pos < len(seq) - kmerLength + 1; ++pos){
						std::string k = seq.seq_.substr(pos, kmerLength);
						++kmerContsVector[k[0]][k[1]][k];
					}
				}
			}
			out << "unordered_map_vector2"
					 << "\t" << kmerLength
					 << "\t" << 4
					 << "\t" << watch.totalTime() << std::endl;;
		}
	}


	//1 depth vector

	{
		std::vector<std::unordered_map<std::string, uint32_t>> kmerContsVector(117);
		{
			njh::scopedStopWatch watch("unordered_map_vector1");

			for(const auto & seq : seqs){
				if(len(seq) > kmerLength){
					for(uint32_t pos = 0; pos < len(seq) - kmerLength + 1; ++pos){
						std::string k = seq.seq_.substr(pos, kmerLength);
						++kmerContsVector[k[0]][k];
					}
				}
			}
			out << "unordered_map_vector1"
					 << "\t" << kmerLength
					 << "\t" << 3
					 << "\t" << watch.totalTime() << std::endl;;
		}
	}

	{
		std::vector<std::unordered_map<std::string, uint32_t>> kmerContsVector(117);
		{
			njh::scopedStopWatch watch("unordered_map_vector1");

			for(const auto & seq : seqs){
				if(len(seq) > kmerLength){
					for(uint32_t pos = 0; pos < len(seq) - kmerLength + 1; ++pos){
						std::string k = seq.seq_.substr(pos, kmerLength);
						++kmerContsVector[k[0]][k];
					}
				}
			}
			out << "unordered_map_vector1"
					 << "\t" << kmerLength
					 << "\t" << 4
					 << "\t" << watch.totalTime() << std::endl;;
		}
	}

	{
		std::unordered_map<std::string, uint32_t> kmerCounts;
		{
			njh::scopedStopWatch watch("unordered_map");

			for(const auto & seq : seqs){
				if(len(seq) > kmerLength){
					for(uint32_t pos = 0; pos < len(seq) - kmerLength + 1; ++pos){
						std::string k = seq.seq_.substr(pos, kmerLength);
						++kmerCounts[k];
					}
				}
			}
			out << "unordered_map"
					 << "\t" << kmerLength
					 << "\t" << 3
					 << "\t" << watch.totalTime() << std::endl;;
		}
	}


	{
		std::unordered_map<std::string, uint32_t> kmerCounts;
		{
			njh::scopedStopWatch watch("unordered_map");

			for(const auto & seq : seqs){
				if(len(seq) > kmerLength){
					for(uint32_t pos = 0; pos < len(seq) - kmerLength + 1; ++pos){
						std::string k = seq.seq_.substr(pos, kmerLength);
						++kmerCounts[k];
					}
				}
			}
			out << "unordered_map"
					 << "\t" << kmerLength
					 << "\t" << 4
					 << "\t" << watch.totalTime() << std::endl;
		}
	}


	return 0;
}




}  // namespace njhseq
