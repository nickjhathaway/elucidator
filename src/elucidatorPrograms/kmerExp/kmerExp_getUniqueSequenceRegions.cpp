//
// Created by Nicholas Hathaway on 7/6/23.
//


#include "kmerExp.hpp"
#include <njhseq/objects/kmer/KmerUtils.hpp>
#include "elucidator/objects/seqObjects/seqKmers.h"
#include <njhseq/IO/SeqIO/SeqIO.hpp>
#include <njhseq/objects/kmer/KmerGatherer.hpp>
#include <njhseq/objects/kmer/SimpleKmerHash.hpp>


namespace njhseq {

int kmerExpRunner::getUniqueSequenceRegions(const njh::progutils::CmdArgs &inputCommands) {
	uint32_t minLen = 50;
	std::vector<bfs::path> inputSequenceFiles;
	KmerGatherer::KmerGathererPars countPars;

	seqSetUp setUp(inputCommands);
	setUp.description_ = "Takes an input of fasta files (e.g. genomes) and then walks across the records within each file to report stretches of sequence of min kmer size completely unique to that fasta file";
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(inputSequenceFiles, "--inputSequenceFiles", "Input list of fasta files", true);
	setUp.setOption(minLen, "--minLen", "minimum segment size");
	countPars.setOptions(setUp);
	setUp.processDirectoryOutputName("uniqueSegments_TODAY", true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);

	std::function<bool(const std::string&)> seqCheck = [&countPars](const std::string & k){
		return std::all_of(k.begin(), k.end(), [&countPars](char base){return njh::in(base, countPars.allowableCharacters_);});
	};
	KmerGatherer kGather(countPars);

	setUp.rLog_.logCurrentTime("count_all");
	setUp.rLog_.runLogFile_.flush();
	std::unordered_map<std::string, std::set<uint64_t>> allKmers = kGather.getUniqueKmersSetHashWithFiltersFromFastas(inputSequenceFiles);
	SimpleKmerHash hasher;

	struct GrowingSegment{
		GrowingSegment() = default;
		uint32_t start_{std::numeric_limits<uint32_t>::max()};
		uint32_t size_{std::numeric_limits<uint32_t>::max()};
	};
	for(const auto & input : inputSequenceFiles){
		setUp.rLog_.logCurrentTime("checking: " + std::string(bfs::basename(bfs::path(input).replace_extension(""))));

		OutputStream out(njh::files::make_path(setUp.pars_.directoryName_, bfs::basename(bfs::path(input).replace_extension("")) + ".bed"));

		auto inputOpts = SeqIOOptions::genFastaIn(input);
		if(countPars.allUpper_){
			inputOpts.lowerCaseBases_ = "upper";
		}
		inputOpts.includeWhiteSpaceInName_ = false;
		seqInfo seq;
		SeqInput reader(inputOpts);
		reader.openIn();


		while(reader.readNextRead(seq)){
			if(len(seq) < countPars.kmerLength_ + 1){
				continue;
			}
			GrowingSegment currentSegment;

			for (uint32_t pos = 0; pos < len(seq) - countPars.kmerLength_ + 1; ++pos) {
				auto k = seq.seq_.substr(pos, countPars.kmerLength_);
				kmerInfo kinfo(k, countPars.kmerLengthForEntropyCalc_, false);
				if (seqCheck(k) && kinfo.computeKmerEntropy() > countPars.entropyFilter_) {
					bool unique = true;
					auto hashk = hasher.hash(k);
					for(const auto & other : allKmers){
						if(other.first != input){
							if(njh::in(hashk, other.second)){
								unique = false;
								break;
							}
						}
					}
					if(unique){
						if(std::numeric_limits<uint32_t>::max() == currentSegment.start_){
							//start
							currentSegment.start_ = pos;
							currentSegment.size_ = countPars.kmerLength_;
						} else {
							if(currentSegment.start_ + currentSegment.size_ - countPars.kmerLength_ + 1 == pos){
								//grow
								++currentSegment.size_;
							} else {
								//broke
								//write out if pass filters
								if (currentSegment.size_ >= minLen) {
									out << seq.name_
											<< "\t" << currentSegment.start_
											<< "\t" << currentSegment.start_ + currentSegment.size_
											<< std::endl;
								}
								//start
								currentSegment.start_ = pos;
								currentSegment.size_ = countPars.kmerLength_;
							}
						}
					}
				}
			}
			if(std::numeric_limits<uint32_t>::max() != currentSegment.start_ && currentSegment.size_ >= minLen){
				out << seq.name_
						<< "\t" << currentSegment.start_
						<< "\t" << currentSegment.start_ + currentSegment.size_
						<< std::endl;
			}
		}
	}
	return 0;
}

} // namespace njhseq