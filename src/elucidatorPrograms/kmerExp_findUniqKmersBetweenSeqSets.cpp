/*
 * kmerExp_findUniqKmersBetweenSeqSets.cpp
 *
 *  Created on: Jul 10, 2021
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

#include "kmerExp.hpp"
#include "elucidator/utils/KmerUtils.hpp"
#include "elucidator/objects/dataContainers.h"
#include "elucidator/simulation.h"
#include "elucidator/objects/seqObjects/seqKmers.h"

#include "elucidator/objects/MiscUtility/GenomeSeqSearch.hpp"


namespace njhseq {

int kmerExpRunner::findUniqKmersBetweenSeqSets(const njh::progutils::CmdArgs & inputCommands){
	bool noRevComp = false;
	uint32_t kmerLength = 21;
	bfs::path genomeFnp1 = "";
	bfs::path genomeFnp2 = "";
	bool writeOutAllCounts = false;
	std::vector<char> allowableCharacters{'A', 'C', 'G', 'T'};
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(genomeFnp1, "--genomeFnp1", "genomeFnp1", true);
	setUp.setOption(genomeFnp2, "--genomeFnp2", "genomeFnp2", true);
	setUp.setOption(kmerLength, "--kmerLength", "kmer Length");
	setUp.setOption(allowableCharacters, "--allowableCharacters", "Only count kmers with these allowable Characters");
	setUp.setOption(writeOutAllCounts, "--writeOutAllCounts", "writeOutAllCounts");
	setUp.setOption(noRevComp, "--noRevComp", "noRevComp");


	std::string genomeFnp1Base = bfs::basename(njh::files::removeExtension(genomeFnp1));
	std::string genomeFnp2Base = bfs::basename(njh::files::removeExtension(genomeFnp2));

	setUp.processDirectoryOutputName(njh::pasteAsStr(genomeFnp1Base, "_vs_", genomeFnp2Base, "_TODAY"), true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);



	std::unordered_map<std::string, uint32_t> genome1Counts;
	std::unordered_map<std::string, uint32_t> genome2Counts;
	setUp.rLog_.setCurrentLapName("initial");
	setUp.rLog_.logCurrentTime("genome1_counting-" + genomeFnp1Base);
	{
		seqInfo seq;
		SeqInput reader(SeqIOOptions::genFastaIn(genomeFnp1));
		reader.openIn();
		while(reader.readNextRead(seq)){
			if(len(seq) > kmerLength){

				for(uint32_t pos = 0; pos < len(seq) - kmerLength + 1; ++pos){
					++genome1Counts[seq.seq_.substr(pos, kmerLength)];
				}
				if(!noRevComp){
					seq.reverseComplementRead(false, true);
					for(uint32_t pos = 0; pos < len(seq) - kmerLength + 1; ++pos){
						++genome1Counts[seq.seq_.substr(pos, kmerLength)];
					}
				}
			}
		}
	}
	setUp.rLog_.logCurrentTime("genome2_counting-" + genomeFnp2Base);

	{
		seqInfo seq;
		SeqInput reader(SeqIOOptions::genFastaIn(genomeFnp2));
		reader.openIn();
		while(reader.readNextRead(seq)){
			if(len(seq) > kmerLength){
				for(uint32_t pos = 0; pos < len(seq) - kmerLength + 1; ++pos){
					++genome2Counts[seq.seq_.substr(pos, kmerLength)];
				}
				if(!noRevComp){
					seq.reverseComplementRead(false, true);
					for(uint32_t pos = 0; pos < len(seq) - kmerLength + 1; ++pos){
						++genome2Counts[seq.seq_.substr(pos, kmerLength)];
					}
				}
			}
		}
	}

	if(writeOutAllCounts){
		setUp.rLog_.logCurrentTime("genome1_writing_all-" + genomeFnp1Base);

		OutputStream genome1CountsOut(njh::files::make_path(setUp.pars_.directoryName_, genomeFnp1Base + "_counts.tab.txt.gz"));
		for(const auto & count : genome1Counts){
			genome1CountsOut << count.first << "\t" << count.second << "\n";
		}
		setUp.rLog_.logCurrentTime("genome2_writing_all-" + genomeFnp2Base);

		OutputStream genome2CountsOut(njh::files::make_path(setUp.pars_.directoryName_, genomeFnp2Base + "_counts.tab.txt.gz"));
		for(const auto & count : genome2Counts){
			genome2CountsOut << count.first << "\t" << count.second << "\n";
		}
	}


	std::function<bool(char)> charCheck = [&allowableCharacters](char base){
		return njh::in(base, allowableCharacters);
	};
	setUp.rLog_.logCurrentTime("genome1_unique_deterination-" + genomeFnp1Base);

	OutputStream genome1UniqueCountsOut(njh::files::make_path(setUp.pars_.directoryName_, genomeFnp1Base + "_uniqueKmers_counts.tab.txt.gz"));
	for(const auto & count : genome1Counts){
		if(std::all_of(count.first.begin(), count.first.end(), charCheck) && !njh::in(count.first, genome2Counts) ){
			genome1UniqueCountsOut << count.first << "\t" << count.second << "\n";
		}
	}
	setUp.rLog_.logCurrentTime("genome2_unique_deterination-" + genomeFnp2Base);
	OutputStream genome2UniqueCountsOut(njh::files::make_path(setUp.pars_.directoryName_, genomeFnp2Base + "_uniqueKmers_counts.tab.txt.gz"));
	for(const auto & count : genome2Counts){
		if(std::all_of(count.first.begin(), count.first.end(), charCheck) && !njh::in(count.first, genome1Counts) ){
			genome2UniqueCountsOut << count.first << "\t" << count.second << "\n";
		}
	}
	setUp.rLog_.logCurrentTime("end");

	return 0;
}

}  //namespace njhseq

