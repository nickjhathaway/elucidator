/*
 * seqSearching.cpp
 *
 *  Created on: Nov 12, 2018
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

#include "seqSearching.hpp"

#include "elucidator/objects/BioDataObject.h"
#include "elucidator/BamToolsUtils/BamUtilities.hpp"
#include "elucidator/BioRecordsUtils/BedUtility.hpp"
#include <njhseq/seqToolsUtils/tandemRepeatUtils.hpp>
#include <njhseq/IO/SeqIO/SeqIO.hpp>
#include <njhseq/system.h>
#include <SeekDeep/objects/TarAmpSetupUtils/PrimersAndMids.hpp>

namespace njhseq {


seqSearchingRunner::seqSearchingRunner()
    : njh::progutils::ProgramRunner(
          {
					 addFunc("findHomopolymerLocations", findHomopolymerLocations, false),
					 addFunc("chopAndMap", chopAndMap, false),
					 addFunc("chopAndMapAndRefine", chopAndMapAndRefine, false),
					 addFunc("findMotifLocations", findMotifLocations, false),
					 addFunc("findTandemMotifLocations", findTandemMotifLocations, false),
					 addFunc("chopAndMapAndRefineInvidual", chopAndMapAndRefineInvidual, false),
					 addFunc("findSimpleTandemRepeatLocations", findSimpleTandemRepeatLocations, false),
          	addFunc("extractBetweenTwoMotifLocations", extractBetweenTwoMotifLocations, false),
           },//
          "seqSearching") {}


int seqSearchingRunner::findHomopolymerLocations(const njh::progutils::CmdArgs & inputCommands){
	OutOptions outOpts(bfs::path(""), ".bed");
	char base = 'N';
	uint32_t minLen = 1;
	uint32_t maxLen = std::numeric_limits<uint32_t>::max();
	seqSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processVerbose();
	setUp.processReadInNames(VecStr{"--fasta", "--fastagz", "--fastqgz", "--fastq"}, true);
	setUp.setOption(base, "--base", "Base to search for");
	setUp.setOption(minLen, "--minLen", "minimum len");
	setUp.setOption(maxLen, "--maxLen", "maximum len");
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

  std::regex e(R"(\d+)");
	std::regex pat(njh::pasteAsStr("(", base, "+)"));
	seqInfo seq;

	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn()	;

	OutputStream out(outOpts);

	while(reader.readNextRead(seq)){
		std::sregex_iterator iter(seq.seq_.begin(), seq.seq_.end(), pat);
		std::sregex_iterator end;
		while (iter != end) {
			auto len = iter->str().size();
			if(len >= minLen && len <= maxLen){
				out << seq.name_
						<< "\t" << iter->position()
						<< "\t" << iter->position() + len
						<< "\t" << base << "x" << len
						<< "\t" << len
						<< "\t" << '+'
						<< std::endl;
			}
			++iter;
		}
	}
	return 0;
}


int seqSearchingRunner::findSimpleTandemRepeatLocations(const njh::progutils::CmdArgs & inputCommands){


	SimpleTandemRepeatFinder::SimpleTRFinderLocsPars pars;

	seqSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processVerbose();

	setUp.processReadInNames(VecStr{"--fasta", "--fastagz", "--fastqgz", "--fastq"}, true);
	pars.verbose = setUp.pars_.verbose_;
	pars.debug = setUp.pars_.debug_;
	pars.setDefaultOpts(setUp);

	setUp.finishSetUp(std::cout);
	SimpleTandemRepeatFinder trfinder(pars);

	trfinder.runSimpleTRFinderLocs(setUp.pars_.ioOptions_);


	return 0;
}


int seqSearchingRunner::findTandemMotifLocations(const njh::progutils::CmdArgs & inputCommands){
	std::string motifstr;
	//bfs::path genomeFnp = "";

	uint32_t allowableErrors = 0;
	uint32_t maxAllowableErrors = 0;
	OutOptions outOpts(bfs::path(""));
	outOpts.outExtention_ = ".bed";
	uint32_t repeatCutOff = 3;
	bool noReverse = false;
	seqSetUp setUp(inputCommands);
	setUp.setOption(noReverse, "--noReverse", "Don't look in reverse complement");

	setUp.pars_.ioOptions_.includeWhiteSpaceInName_ = false;
	//setUp.setOption(genomeFnp, "--genomeFnp", "The genome file to look for motifs in", true);
	setUp.processReadInNames({"--fasta", "--fastagz"}, true);
	setUp.setOption(motifstr, "--motif", "The motif to look for", true);
	setUp.setOption(repeatCutOff, "--repeatCutOff", "The minimum number of times the motif repeats to report it");
	setUp.setOption(allowableErrors, "--allowableErrors", "allowable errors in a motif element");
	maxAllowableErrors = allowableErrors;
	setUp.setOption(maxAllowableErrors, "--maxAllowableErrors", "max allowable errors in the whole tandem motif sequence");

	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	seqInfo seq;
	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	motif mot(motifstr);
	OutputStream out(outOpts);
	while(reader.readNextRead(seq)){
		auto locs = mot.findPositionsFull(seq.seq_, allowableErrors);
		njh::sort(locs);
		if(!locs.empty()){
			uint32_t length = 1;
			size_t start = locs.front();
			for(const auto pos : iter::range<uint32_t>(1, locs.size())){
				if(locs[pos] == locs[pos - 1] + mot.size() ){
					++length;
				} else {
					if(length >= repeatCutOff){
						uint32_t numOfErrors = 0;
						for(const auto seqPos : iter::range(start, start + mot.size() * length, mot.size())){
							numOfErrors += mot.size() - mot.scoreMotif(seq.seq_.begin() + seqPos, seq.seq_.begin() + seqPos + mot.size());
						}
						if(numOfErrors <= maxAllowableErrors){
							out << seq.name_
									<< "\t" << start
									<< "\t" << start + mot.size() * length
									<< "\t" << motifstr << "_x" << length
									<< "\t" << length * mot.size()
									<< "\t" << "+" << '\n';;
						}
					}
					length = 1;
					start = locs[pos];
				}
			}
			if(length >= repeatCutOff){
				uint32_t numOfErrors = 0;
				for(const auto seqPos : iter::range(start, start + mot.size() * length, mot.size())){
					numOfErrors += mot.size() - mot.scoreMotif(seq.seq_.begin() + seqPos, seq.seq_.begin() + seqPos + mot.size());
				}
				if(numOfErrors <= maxAllowableErrors){
					out << seq.name_
							<< "\t" << start
							<< "\t" << start + mot.size() * length
							<< "\t" << motifstr << "_x" << length
							<< "\t" << length * mot.size()
							<< "\t" << "+" << '\n';;
				}
			}
		}
		if(!noReverse){
			seq.reverseComplementRead(false, true);
			auto revLocs = mot.findPositionsFull(seq.seq_, allowableErrors);
			if(!revLocs.empty()){
				njh::sort(revLocs);
				uint32_t length = 1;
				size_t start = revLocs.front();
				for(const auto  pos : iter::range<uint32_t>(1, revLocs.size())){
					if(revLocs[pos] == revLocs[pos - 1] + mot.size()){
						++length;
					}else{
						if(length >= repeatCutOff){
							uint32_t numOfErrors = 0;
							for(const auto seqPos : iter::range(start, start + mot.size() * length, mot.size())){
								numOfErrors += mot.size() - mot.scoreMotif(seq.seq_.begin() + seqPos, seq.seq_.begin() + seqPos + mot.size());
							}
							if(numOfErrors <= maxAllowableErrors){
								out << seq.name_
										<< "\t" << len(seq) - (start + mot.size() * length )
										<< "\t" << len(seq) - start
										<< "\t" << motifstr << "_x" << length
										<< "\t" << length * mot.size()
										<< "\t" << "-" << '\n';;
							}
						}
						length = 1;
						start = revLocs[pos];
					}
				}
				if(length >= repeatCutOff){
					uint32_t numOfErrors = 0;
					for(const auto seqPos : iter::range(start, start + mot.size() * length, mot.size())){
						numOfErrors += mot.size() - mot.scoreMotif(seq.seq_.begin() + seqPos, seq.seq_.begin() + seqPos + mot.size());
					}
					if(numOfErrors <= maxAllowableErrors){
						out << seq.name_
								<< "\t" << len(seq) - (start + mot.size() * length)
								<< "\t" << len(seq) - start
								<< "\t" << motifstr << "_x" << length
								<< "\t" << length * mot.size()
								<< "\t" << "-" << '\n';;
					}
				}
			}
		}
	}

	return 0;
}




int seqSearchingRunner::extractBetweenTwoMotifLocations(const njh::progutils::CmdArgs & inputCommands){
	uint32_t numThreads = 1;
	std::vector<bfs::path> fasta_list;
	bfs::path motif_pair_table;
	std::string motif_pair_name = "motif";
	seqInfo motif1Obj;
	seqInfo motif2Obj;
	//bfs::path genomeFnp = "";
	uint32_t allowableErrors = 0;
	size_t insertSizeCutOff = std::numeric_limits<size_t>::max();
	OutOptions outOpts(bfs::path(""));
	outOpts.outExtention_ = ".bed";
	seqSetUp setUp(inputCommands);
	setUp.pars_.ioOptions_.includeWhiteSpaceInName_ = false;
	bool fasta_list_set = setUp.setOption(fasta_list, "--fasta_list", "a list of fasta files to read from");
	bool motif_pair_table_set = setUp.setOption(motif_pair_table, "--motif_pair_table", "a table with 3 column, 1)target,2)motif1(5`-3` direction), 3)motif2 (5`-3` direction)");

	setUp.processReadInNames({"--fasta", "--fastagz"}, !fasta_list_set);
	setUp.processSeq(motif1Obj, "--motif1", "The first motif to look for", !motif_pair_table_set);
	setUp.processSeq(motif2Obj, "--motif2", "The second motif to look for, should be in reverse complement to motif1", !motif_pair_table_set);
	setUp.setOption(motif_pair_name, "--motif_pair_name", "motif pair name", false);
	setUp.setOption(allowableErrors, "--allowableErrors", "allowable errors in motif");
	setUp.setOption(insertSizeCutOff, "--insertSizeCutOff", "max insert Size Cut Off");
	setUp.setOption(numThreads, "--numThreads", "number of threads to use when processing a list of fasta files");

	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	std::shared_ptr<PrimersAndMids> motifs;
	if (motif_pair_table_set) {
		motifs = std::make_shared<PrimersAndMids>(motif_pair_table);
	} else {
		motifs = std::make_shared<PrimersAndMids>(std::map{
			std::make_pair(motif_pair_name, PrimersAndMids::Target(motif_pair_name, motif1Obj.seq_, motif2Obj.seq_))
		});
	}

	motifs->initPrimerDeterminator();

	OutputStream out(outOpts);


	struct MotifPairSearchResults {
        std::vector<GenomicRegion> fPrimerPositions_;
        std::vector<GenomicRegion> rPrimerPositions_;
        class GenomeExtractResultByMotifPair {
        public:
          GenomeExtractResultByMotifPair(GenomicRegion  fPrimerReg,
                                          GenomicRegion  rPrimerReg): fPrimerReg_(std::move(fPrimerReg)),rPrimerReg_(std::move(rPrimerReg)){
            setRegion();
          }
          GenomicRegion fPrimerReg_;
          GenomicRegion rPrimerReg_;

          std::shared_ptr<GenomicRegion> gRegion_;
          std::shared_ptr<GenomicRegion> gRegionInner_ ;

          void setRegion(){
            if (fPrimerReg_.chrom_ != rPrimerReg_.chrom_) {
              std::stringstream ss;
              ss << __PRETTY_FUNCTION__ << ", error extention chrom, "
                 << fPrimerReg_.chrom_ << "doesn't equal ligation chrom "
                 << rPrimerReg_.chrom_ << "\n";
              throw std::runtime_error { ss.str() };
            }
            if (fPrimerReg_.reverseSrand_ == rPrimerReg_.reverseSrand_) {
              std::stringstream ss;
              ss << __PRETTY_FUNCTION__
                 << ", error extention and ligation are on the same strand, should be mapping to opposite strands"
                 << "\n";
              throw std::runtime_error { ss.str() };
            }
            if (fPrimerReg_.reverseSrand_) {
              if (fPrimerReg_.start_ < rPrimerReg_.start_) {
                std::stringstream ss;
                ss << __PRETTY_FUNCTION__
                   << ", error if extention is mapping to the reverse strand, it's start, "
                   << fPrimerReg_.start_
                   << ", should be greater than ligation start, "
                   << rPrimerReg_.start_ << "\n";
                throw std::runtime_error { ss.str() };
              }
            }else{
              if (fPrimerReg_.start_ > rPrimerReg_.start_) {
                std::stringstream ss;
                ss << __PRETTY_FUNCTION__
                   << ", error if extention is mapping to the plus strand, it's start, "
                   << fPrimerReg_.start_
                   << ", should be less than than ligation start, "
                   << rPrimerReg_.start_ << "\n";
                throw std::runtime_error { ss.str() };
              }
            }
            size_t start = fPrimerReg_.start_;
            size_t end = rPrimerReg_.end_;
            size_t innerStart = fPrimerReg_.end_;
            size_t innereEnd = rPrimerReg_.start_;
            if(fPrimerReg_.reverseSrand_){
              start = rPrimerReg_.start_;
              end = fPrimerReg_.end_;
              innerStart = rPrimerReg_.end_;
              innereEnd = fPrimerReg_.start_;
            }

            gRegion_ = std::make_shared<GenomicRegion>(fPrimerReg_.uid_ + "-" + rPrimerReg_.uid_, fPrimerReg_.chrom_, start, end, fPrimerReg_.reverseSrand_);
            gRegionInner_ = std::make_shared<GenomicRegion>(fPrimerReg_.uid_ + "-" + rPrimerReg_.uid_, fPrimerReg_.chrom_, innerStart, innereEnd, fPrimerReg_.reverseSrand_);
          }

        };


        static std::vector<GenomeExtractResultByMotifPair> getPossibleGenomeExtracts(const std::vector<GenomicRegion> & fPrimerPositions,
                                                                                      const std::vector<GenomicRegion> & rPrimerPositions,
                                                                                      const size_t insertSizeCutOff = std::numeric_limits<size_t>::max()){
          std::vector<GenomeExtractResultByMotifPair> ret;
          //same chrom, opposite strands, less than the insert size
          for (const auto & fwd : fPrimerPositions) {
            for (const auto & rev : rPrimerPositions) {
              //need to be on the same chromosome
              //need to be on opposite strands (should both should be in 5'->3' direction
              //and they shouldn't overlap
              if (fwd.chrom_ == rev.chrom_
                  && fwd.reverseSrand_ != rev.reverseSrand_
                  && !fwd.overlaps(rev)
                  && fwd.start_ != rev.end_
                  && fwd.end_ != rev.start_ ) {

                if(fwd.reverseSrand_){
                  if(fwd.start_ > rev.start_){
                    GenomeExtractResultByMotifPair extraction(fwd, rev);
                    if (extraction.gRegion_->getLen() <= insertSizeCutOff) {
                      ret.emplace_back(extraction);
                    }
                  }
                }else{
                  if(fwd.start_ < rev.start_){
                    GenomeExtractResultByMotifPair extraction(fwd, rev);
                    if (extraction.gRegion_->getLen() <= insertSizeCutOff) {
                      ret.emplace_back(extraction);
                    }
                  }
                }
              }
            }
          }
          return ret;
        }

        std::vector<GenomeExtractResultByMotifPair> regions_;


      };
	if (fasta_list_set) {
		njh::concurrent::LockableVec<bfs::path> fasta_inputs(fasta_list);
		std::mutex out_mut;
		std::function search_for_motifs = [&out_mut, &fasta_inputs,&motifs, insertSizeCutOff, &out, allowableErrors]() {
			bfs::path current_fasta;
			while (fasta_inputs.getVal(current_fasta)) {
				seqInfo fwd_seq;
				SeqInput reader(SeqIOOptions::genFastaIn(current_fasta));
				reader.openIn();
				while(reader.readNextRead(fwd_seq)){
					auto revComp = seqUtil::reverseComplement(fwd_seq.seq_,"DNA");
					for (const auto & motif_pair : motifs->pDeterminator_->primers_) {
						auto locs_1 = motif_pair.second.fwds_.front().mot_.findPositionsFull(fwd_seq.seq_, allowableErrors);
						auto locs_2 = motif_pair.second.revs_.front().mot_.findPositionsFull(fwd_seq.seq_, allowableErrors);
						auto revLocs_1 = motif_pair.second.fwds_.front().mot_.findPositionsFull(revComp, allowableErrors);
						auto revLocs_2 = motif_pair.second.revs_.front().mot_.findPositionsFull(revComp, allowableErrors);
						if ((!locs_1.empty() && !revLocs_2 .empty()) ||
								(!revLocs_1.empty() && !locs_2.empty())) {
							std::vector<GenomicRegion> motif1Positions;
							std::vector<GenomicRegion> motif2Positions;
							for(const auto & loc : locs_1){
								motif1Positions.emplace_back(motif_pair.first, fwd_seq.name_,loc, loc + motif_pair.second.fwds_.front().mot_.size(), false);
							}
							for(const auto & loc : revLocs_1){
								motif1Positions.emplace_back(motif_pair.first, fwd_seq.name_, len(fwd_seq) - (loc + motif_pair.second.fwds_.front().mot_.size()), len(fwd_seq) - loc, true);
							}
							for(const auto & loc : locs_2){
								motif2Positions.emplace_back(motif_pair.first, fwd_seq.name_,loc, loc + motif_pair.second.revs_.front().mot_.size(), false);
							}
							for(const auto & loc : revLocs_2){
								motif2Positions.emplace_back(motif_pair.first, fwd_seq.name_,len(fwd_seq) - (loc + motif_pair.second.revs_.front().mot_.size()), len(fwd_seq) - loc, true);
							}
							auto possible_extractions = MotifPairSearchResults::getPossibleGenomeExtracts(motif1Positions, motif2Positions, insertSizeCutOff);
							{
								std::lock_guard<std::mutex> lock(out_mut);
								for (const auto & extract : possible_extractions) {
									auto out_bed = extract.gRegion_->genBedRecordCore();
									out_bed.extraFields_.emplace_back(current_fasta.string());
									out << out_bed.toDelimStrWithExtra() << std::endl;
								}
							}
						}
					}
				}
			}
		};
		njh::concurrent::runVoidFunctionThreaded(search_for_motifs, numThreads);
	} else {
		seqInfo fwd_seq;
		SeqInput reader(setUp.pars_.ioOptions_);
		reader.openIn();

		while(reader.readNextRead(fwd_seq)){
			auto revComp = seqUtil::reverseComplement(fwd_seq.seq_,"DNA");
			for (const auto & motif_pair : motifs->pDeterminator_->primers_) {
				auto locs_1 = motif_pair.second.fwds_.front().mot_.findPositionsFull(fwd_seq.seq_, allowableErrors);
				auto locs_2 = motif_pair.second.revs_.front().mot_.findPositionsFull(fwd_seq.seq_, allowableErrors);
				auto revLocs_1 = motif_pair.second.fwds_.front().mot_.findPositionsFull(revComp, allowableErrors);
				auto revLocs_2 = motif_pair.second.revs_.front().mot_.findPositionsFull(revComp, allowableErrors);
				if ((!locs_1.empty() && !revLocs_2 .empty()) ||
						(!revLocs_1.empty() && !locs_2.empty())) {
					std::vector<GenomicRegion> motif1Positions;
					std::vector<GenomicRegion> motif2Positions;
					for(const auto & loc : locs_1){
						motif1Positions.emplace_back(motif_pair.first, fwd_seq.name_,loc, loc + motif_pair.second.fwds_.front().mot_.size(), false);
					}
					for(const auto & loc : revLocs_1){
						motif1Positions.emplace_back(motif_pair.first, fwd_seq.name_, len(fwd_seq) - (loc + motif_pair.second.fwds_.front().mot_.size()), len(fwd_seq) - loc, true);
					}
					for(const auto & loc : locs_2){
						motif2Positions.emplace_back(motif_pair.first, fwd_seq.name_,loc, loc + motif_pair.second.revs_.front().mot_.size(), false);
					}
					for(const auto & loc : revLocs_2){
						motif2Positions.emplace_back(motif_pair.first, fwd_seq.name_,len(fwd_seq) - (loc + motif_pair.second.revs_.front().mot_.size()), len(fwd_seq) - loc, true);
					}
					auto possible_extractions = MotifPairSearchResults::getPossibleGenomeExtracts(motif1Positions, motif2Positions, insertSizeCutOff);
					for (const auto & extract : possible_extractions) {
						out << extract.gRegion_->genBedRecordCore().toDelimStrWithExtra() << std::endl;
					}
				}
			}
		}
	}

	return 0;
}


int seqSearchingRunner::findMotifLocations(const njh::progutils::CmdArgs & inputCommands){
	std::string motifstr = "";
	seqInfo motifObj;
	//bfs::path genomeFnp = "";
	uint32_t allowableErrors = 0;
	OutOptions outOpts(bfs::path(""));
	outOpts.outExtention_ = ".bed";
	bool noReverse = false;
	seqSetUp setUp(inputCommands);
	setUp.pars_.ioOptions_.includeWhiteSpaceInName_ = false;
	//setUp.setOption(genomeFnp, "--genomeFnp", "The genome file to look for motifs in", true);
	setUp.processReadInNames({"--fasta", "--fastagz"}, true);
	//setUp.setOption(motifstr, "--motif", "The motif to look for", true);
	setUp.processSeq(motifObj, "--motif", "The motif to look for", true);
	setUp.setOption(allowableErrors, "--allowableErrors", "allowable errors in motif");
	setUp.setOption(noReverse, "--noReverse", "Don't look in reverse complement");
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	motifstr = motifObj.seq_;
	seqInfo seq;
	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	motif mot(motifstr);
	OutputStream out(outOpts);
	while(reader.readNextRead(seq)){
		auto locs = mot.findPositionsFull(seq.seq_, allowableErrors);
		for(const auto & loc : locs){
			out << seq.name_
					<< "\t" << loc
					<< "\t" << loc + mot.size()
					<< "\t" << seq.seq_.substr(loc, mot.size())
					<< "\t" << mot.scoreMotif(seq.seq_.begin() + loc, seq.seq_.begin() + loc + mot.size())
					<< "\t" << "+" << '\n';;
		}
		if(!noReverse){
			seq.reverseComplementRead(false, true);
			auto revLocs = mot.findPositionsFull(seq.seq_, allowableErrors);
			for(const auto & loc : revLocs){
				out << seq.name_
						<< "\t" << len(seq) - (loc + mot.size())
						<< "\t" << len(seq) - loc
						<< "\t" << seq.seq_.substr(loc, mot.size())
						<< "\t" << mot.scoreMotif(seq.seq_.begin() + loc, seq.seq_.begin() + loc + mot.size())
						<< "\t" << "-" << '\n';
			}
		}
	}
	return 0;
}


struct ChopAndMapPars{
	uint32_t minSize = 30;
	uint32_t windowLength = 100;
	uint32_t windowStep = 10;
	uint32_t numThreads = 1;
	uint32_t perFragmentCount = 10;
	bfs::path genomeFnp = "";
	SeqIOOptions inOpts;
	bfs::path outputDirectory;
	bool debug = false;
};

void runChopAndMap(const ChopAndMapPars & pars){
	njh::sys::requireExternalProgramThrow("bwa");
	njh::files::checkExistenceThrow(pars.genomeFnp, __PRETTY_FUNCTION__);
	auto fragOutOpts = SeqIOOptions::genFastaOut(njh::files::make_path(pars.outputDirectory, "fragments.fasta"));

	{
		SeqOutput fragWriter(fragOutOpts);
		fragWriter.openOut();
		seqInfo seq;
		SeqInput reader(pars.inOpts);
		reader.openIn();
		while (reader.readNextRead(seq)) {
			if (len(seq) < pars.windowLength ) {
				if(len(seq) <pars.minSize){
					std::cerr << "Seq: " << seq.name_ << " is too short" << std::endl;
				}else{
					auto fragment = seq;
					MetaDataInName fragMeta;
					fragMeta.addMeta("start", 0);
					fragMeta.addMeta("end", 0 + len(seq));
					fragment.name_.append(fragMeta.createMetaName());
					for(uint32_t fragCount = 0; fragCount <= pars.perFragmentCount; ++ fragCount){
						fragWriter.write(fragment);
					}
				}
			} else {
				for(auto const pos : iter::range<uint32_t>(0, len(seq), pars.windowStep)){
					if(pos + pars.windowLength > len(seq)){
						auto newLength = len(seq) - pos;
						if(newLength > pars.minSize){
							auto fragment = seq.getSubRead(pos);
							MetaDataInName fragMeta;
							fragMeta.addMeta("start", pos);
							fragMeta.addMeta("end", len(seq));
							fragment.name_.append(fragMeta.createMetaName());
							for(uint32_t fragCount = 0; fragCount <= pars.perFragmentCount; ++ fragCount){
								fragWriter.write(fragment);
							}
						}
						break;
					}else{
						auto fragment = seq.getSubRead(pos, pars.windowLength);
						MetaDataInName fragMeta;
						fragMeta.addMeta("start", pos);
						fragMeta.addMeta("end", pos + pars.windowLength);
						fragment.name_.append(fragMeta.createMetaName());
						for(uint32_t fragCount = 0; fragCount <= pars.perFragmentCount; ++ fragCount){
							fragWriter.write(fragment);
						}
					}
				}
			}
		}
	}

	bfs::path fragmentBamFnp =  njh::files::make_path(pars.outputDirectory, "fragments.sorted.bam");
	std::stringstream bwaMappCmd;
	bwaMappCmd << "bwa mem -M -t " << pars.numThreads << " "
			<< pars.genomeFnp  << " "
			<< fragOutOpts.out_.outName() << " "
			<< " 2> " << njh::files::make_path(pars.outputDirectory,"bwa.log.txt")
			<< " | samtools sort - -o " << fragmentBamFnp
			<< " && samtools index " << fragmentBamFnp ;
	auto bwaRunOutput = njh::sys::run({bwaMappCmd.str()});
	BioCmdsUtils::checkRunOutThrow(bwaRunOutput, __PRETTY_FUNCTION__);

}

int seqSearchingRunner::chopAndMap(const njh::progutils::CmdArgs & inputCommands) {
	ChopAndMapPars chopPars;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	chopPars.debug = setUp.pars_.debug_;
	setUp.setOption(chopPars.numThreads, "--numThreads", "Number of threads to use");
	setUp.setOption(chopPars.genomeFnp, "--genomeFnp", "Genome to map to", true);
	setUp.setOption(chopPars.perFragmentCount, "--perFragmentCount", "perFragmentCount");
	setUp.setOption(chopPars.windowLength, "--windowLength", "windowLength");
	setUp.setOption(chopPars.windowStep, "--windowStep", "windowStep");
	setUp.processReadInNames(VecStr{"--fasta", "--fastq"});
	setUp.processDirectoryOutputName(true);
	chopPars.inOpts = setUp.pars_.ioOptions_;
	chopPars.outputDirectory = setUp.pars_.directoryName_;
	setUp.finishSetUp(std::cout);

	runChopAndMap(chopPars);

	return 0;
}


int seqSearchingRunner::chopAndMapAndRefineInvidual(const njh::progutils::CmdArgs & inputCommands) {
	ChopAndMapPars globalChopPars;

	RunCoverageFinderSinglePars covPars;
	RegionRefinementPars refinePars;
	uint32_t expandLeft = 0;
	uint32_t expandRight = 0;
	uint32_t minLength = 0;
	uint32_t inputMinLength = 0;
	double entropyCutOff = 1.5;

	bfs::path gff = "";

	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	globalChopPars.debug = setUp.pars_.debug_;
	setUp.setOption(gff, "--genomegff", "Genome gff file");

	setUp.setOption(globalChopPars.numThreads, "--numThreads", "Number of threads to use");
	setUp.setOption(globalChopPars.genomeFnp, "--genomeFnp", "Genome to map to", true);
	setUp.setOption(globalChopPars.perFragmentCount, "--perFragmentCount", "perFragmentCount");
	setUp.setOption(globalChopPars.windowLength, "--windowLength", "windowLength");
	inputMinLength = globalChopPars.windowLength;
	setUp.setOption(globalChopPars.windowStep, "--windowStep", "windowStep");
	setUp.setOption(refinePars.reOrient, "--reOrient", "reOrient");
	setUp.setOption(inputMinLength, "--inputMinLength", "inputMinLength");

	setUp.setOption(expandLeft, "--expandLeft", "expandLeft");
	setUp.setOption(expandRight, "--expandRight", "expandRight");

	setUp.setOption(minLength, "--minLength", "minLength");
	setUp.setOption(entropyCutOff, "--entropyCutOff", "entropy Cut Off for final refined regions ");

	setUp.processReadInNames(VecStr{"--fasta", "--fastq"});
	setUp.processDirectoryOutputName(true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);

	auto seqs = SeqInput::getSeqVec<seqInfo>(setUp.pars_.ioOptions_);
	for(const auto & seq : seqs){
		if(len(seq) < inputMinLength){
			continue;
		}
		auto seqDirFnp = njh::files::make_path(setUp.pars_.directoryName_, seq.name_);
		njh::files::makeDir(njh::files::MkdirPar{seqDirFnp});
		auto seqOut = SeqIOOptions::genFastaOut(njh::files::make_path(seqDirFnp, seq.name_));
		SeqOutput::write(std::vector<seqInfo>{seq}, seqOut);
		auto chopParsCurrent = globalChopPars;

		chopParsCurrent.inOpts = SeqIOOptions::genFastaIn(seqOut.out_.outName());
		chopParsCurrent.outputDirectory = seqDirFnp;
		//chop
		runChopAndMap(chopParsCurrent);
		//determine coverage
		covPars.numThreads = chopParsCurrent.numThreads;
		covPars.coverageCutOff = 5;
		covPars.window = chopParsCurrent.windowLength;
		covPars.step = chopParsCurrent.windowStep;
		covPars.bamFnp = njh::files::make_path(chopParsCurrent.outputDirectory, "fragments.sorted.bam").string();
		covPars.outOpts = OutOptions(njh::files::make_path(chopParsCurrent.outputDirectory, "coverage.bed"));
		RunCoverageFinderSingle(covPars);
		//merge regions
		auto beds = getBed3s(covPars.outOpts.outName());
		njh::for_each(beds,
				[&expandLeft,&expandRight](const std::shared_ptr<Bed3RecordCore> & region) {
					if(0 != expandLeft) {
						BedUtility::extendLeft(*region, expandLeft);
					}
					if(0 != expandRight) {
						BedUtility::extendRight(*region, expandRight);
					}
				});
		njh::sort(beds,
				[](const std::shared_ptr<Bed3RecordCore> & region1, const std::shared_ptr<Bed3RecordCore> & region2) {
					return region1->chrom_ == region2->chrom_ ? (region1->chromStart_ == region2->chromStart_ ? region1->chromEnd_ < region2->chromEnd_ : region1->chromStart_ < region2->chromStart_): region1->chrom_ < region2->chrom_;
				});


		OutOptions mergedCoverageOpts(njh::files::make_path(chopParsCurrent.outputDirectory, "merged_coverage.bed"));
		{
			OutputStream mergedCoverageOut(mergedCoverageOpts);
			Bed3RecordCore currentRegion = *beds.front();
			for(const auto regPos : iter::range<uint32_t>(1, beds.size())){
				if(currentRegion.chrom_ == beds[regPos]->chrom_ && currentRegion.overlaps(*beds[regPos], 1)){
					currentRegion.chromEnd_ = std::max(currentRegion.chromEnd_, beds[regPos]->chromEnd_);
				}else{
					Bed6RecordCore regOut(currentRegion.chrom_, currentRegion.chromStart_, currentRegion.chromEnd_, "", currentRegion.length(), '+');
					regOut.name_ = GenomicRegion(regOut).createUidFromCoords();

					mergedCoverageOut << regOut.toDelimStr() << std::endl;
					currentRegion = *beds[regPos];
				}
			}
			Bed6RecordCore regOut(currentRegion.chrom_, currentRegion.chromStart_, currentRegion.chromEnd_, "", currentRegion.length(), '+');
			regOut.name_ = GenomicRegion(regOut).createUidFromCoords();
			mergedCoverageOut << regOut.toDelimStr() << std::endl;
		}

		//refine regions
		refinePars.bedFnp = mergedCoverageOpts.outName();
		refinePars.bamFnp = njh::files::make_path(chopParsCurrent.outputDirectory, "fragments.sorted.bam");
		refinePars.outOpts = OutOptions(njh::files::make_path(chopParsCurrent.outputDirectory, "refined_merged.bed"));
		auto refinedRegions = RunRegionRefinement(refinePars);

		//extra fasta file of refined regions if 2bit file exists
		auto genomeTwoBitFnp = chopParsCurrent.genomeFnp;
		genomeTwoBitFnp.replace_extension(".2bit");
		if(bfs::exists(genomeTwoBitFnp)){
			auto refinedGRegions = bedPtrsToGenomicRegs(refinedRegions);
			if("" != gff){
				intersectBedLocsWtihGffRecordsPars interPars(gff, VecStr{"description"}, VecStr{"pseudogene", "gene"});
				refinePars.outOpts.overWriteFile_ = true;
				intersectBedLocsWtihGffRecords(refinedRegions, interPars);
				OutputStream bedOut(refinePars.outOpts);
				for(const auto & b : refinedRegions){
					bedOut << b->toDelimStrWithExtra() << std::endl;
				}
			}
			TwoBit::TwoBitFile tReader(genomeTwoBitFnp);
			auto refindedSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(chopParsCurrent.outputDirectory, "refined_merged.fasta"));
			SeqOutput refinedWriter(refindedSeqOpts);
			refinedWriter.openOut();
			for(const auto & reg : refinedGRegions){
				refinedWriter.write(reg.extractSeq(tReader));
			}
		}

		if(0 != minLength || (entropyCutOff > 0 && bfs::exists(genomeTwoBitFnp))){
			std::vector<GenomicRegion> filteredRefinedRegions;
			std::shared_ptr<TwoBit::TwoBitFile> tReaderPrt;
			if(bfs::exists(genomeTwoBitFnp)){
				tReaderPrt = std::make_shared<TwoBit::TwoBitFile>(genomeTwoBitFnp);
			}
			{
				auto bedAgain = getBed3s(refinePars.outOpts.outName());
				OutOptions filteredOpts(njh::files::make_path(chopParsCurrent.outputDirectory, "filtered_refined_merged.bed"));
				OutputStream filteredOut(filteredOpts);
				for(const auto & b : bedAgain){
					if(b->length() >= minLength){
						bool furtherFilt = false;
						if(nullptr != tReaderPrt &&  entropyCutOff > 0){
							auto eSeq = GenomicRegion(*b).extractSeq(*tReaderPrt);
							charCounter counter(eSeq.seq_);
							if(counter.computeEntrophy() < entropyCutOff){
								furtherFilt = true;
							}
						}
						if(!furtherFilt){
							filteredRefinedRegions.emplace_back(*b);
							filteredOut << b->toDelimStrWithExtra() << std::endl;
						}
					}
				}
			}
			if(bfs::exists(genomeTwoBitFnp)){
				TwoBit::TwoBitFile tReader(genomeTwoBitFnp);
				auto refindedSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(chopParsCurrent.outputDirectory, "filtered_refined_merged.fasta"));
				SeqOutput refinedWriter(refindedSeqOpts);
				refinedWriter.openOut();
				for(const auto & reg : filteredRefinedRegions){
					refinedWriter.write(reg.extractSeq(tReader));
				}
			}
		}
	}


	return 0;
}

int seqSearchingRunner::chopAndMapAndRefine(const njh::progutils::CmdArgs & inputCommands) {
	ChopAndMapPars chopPars;
	RunCoverageFinderSinglePars covPars;
	RegionRefinementPars refinePars;
	uint32_t expandLeft = 0;
	uint32_t expandRight = 0;
	uint32_t minLength = 0;
	covPars.coverageCutOff = 5;
	double entropyCutOff = 1.5;

	bfs::path gff = "";

	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	chopPars.debug = setUp.pars_.debug_;
	setUp.setOption(chopPars.numThreads, "--numThreads", "Number of threads to use");
	setUp.setOption(chopPars.genomeFnp, "--genomeFnp", "Genome to map to", true);
	setUp.setOption(gff, "--genomegff", "Genome gff file");
	setUp.setOption(chopPars.minSize, "--minFragmentSize", "Min Fragment Size");

	setUp.setOption(chopPars.perFragmentCount, "--perFragmentCount", "perFragmentCount");
	setUp.setOption(chopPars.windowLength, "--windowLength", "windowLength");
	setUp.setOption(chopPars.windowStep, "--windowStep", "windowStep");
	setUp.setOption(refinePars.reOrient, "--reOrient", "reOrient");

	setUp.setOption(expandLeft, "--expandLeft", "expandLeft");
	setUp.setOption(expandRight, "--expandRight", "expandRight");

	setUp.setOption(covPars.coverageCutOff, "--coverageCutOff", "Coverage Cut Off");

	setUp.setOption(entropyCutOff, "--entropyCutOff", "entropy Cut Off for final refined regions ");

	setUp.setOption(minLength, "--minLength", "minLength");

	setUp.processReadInNames(VecStr{"--fasta", "--fastq"});
	setUp.processDirectoryOutputName(true);
	chopPars.inOpts = setUp.pars_.ioOptions_;
	chopPars.outputDirectory = setUp.pars_.directoryName_;
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	//chop
	runChopAndMap(chopPars);
	//determine coverage
	covPars.numThreads = chopPars.numThreads;
	covPars.window = chopPars.windowLength;
	covPars.step = chopPars.windowStep;
	covPars.bamFnp = njh::files::make_path(chopPars.outputDirectory, "fragments.sorted.bam").string();
	covPars.outOpts = OutOptions(njh::files::make_path(chopPars.outputDirectory, "coverage.bed"));
	RunCoverageFinderSingle(covPars);
	//merge regions
	auto beds = getBed3s(covPars.outOpts.outName());
	njh::for_each(beds,
			[&expandLeft,&expandRight](const std::shared_ptr<Bed3RecordCore> & region) {
				if(0 != expandLeft) {
					BedUtility::extendLeft(*region, expandLeft);
				}
				if(0 != expandRight) {
					BedUtility::extendRight(*region, expandRight);
				}
			});
	njh::sort(beds,
			[](const std::shared_ptr<Bed3RecordCore> & region1, const std::shared_ptr<Bed3RecordCore> & region2) {
				return region1->chrom_ == region2->chrom_ ? (region1->chromStart_ == region2->chromStart_ ? region1->chromEnd_ < region2->chromEnd_ : region1->chromStart_ < region2->chromStart_): region1->chrom_ < region2->chrom_;
			});


	OutOptions mergedCoverageOpts(njh::files::make_path(chopPars.outputDirectory, "merged_coverage.bed"));
	{
		OutputStream mergedCoverageOut(mergedCoverageOpts);
		Bed3RecordCore currentRegion = *beds.front();
		for(const auto regPos : iter::range<uint32_t>(1, beds.size())){
			if(currentRegion.chrom_ == beds[regPos]->chrom_ && currentRegion.overlaps(*beds[regPos], 1)){
				currentRegion.chromEnd_ = std::max(currentRegion.chromEnd_, beds[regPos]->chromEnd_);
			}else{
				Bed6RecordCore regOut(currentRegion.chrom_, currentRegion.chromStart_, currentRegion.chromEnd_, "", currentRegion.length(), '+');
				regOut.name_ = GenomicRegion(regOut).createUidFromCoords();

				mergedCoverageOut << regOut.toDelimStr() << std::endl;
				currentRegion = *beds[regPos];
			}
		}
		Bed6RecordCore regOut(currentRegion.chrom_, currentRegion.chromStart_, currentRegion.chromEnd_, "", currentRegion.length(), '+');
		regOut.name_ = GenomicRegion(regOut).createUidFromCoords();
		mergedCoverageOut << regOut.toDelimStr() << std::endl;
	}

	//refine regions
	refinePars.bedFnp = mergedCoverageOpts.outName();
	refinePars.bamFnp = njh::files::make_path(chopPars.outputDirectory, "fragments.sorted.bam");
	refinePars.outOpts = OutOptions(njh::files::make_path(chopPars.outputDirectory, "refined_merged.bed"));
	auto refinedRegions = RunRegionRefinement(refinePars);

	/**@todo report number of reads unmapped and what fragments these are*/

	//extra fasta file of refined regions if 2bit file exists
	auto genomeTwoBitFnp = chopPars.genomeFnp;
	genomeTwoBitFnp.replace_extension(".2bit");
	if(bfs::exists(genomeTwoBitFnp)){
		auto refinedGRegions = bedPtrsToGenomicRegs(refinedRegions);
		if("" != gff){
			intersectBedLocsWtihGffRecordsPars interPars(gff,VecStr{"description"}, VecStr{"pseudogene", "gene"});
			refinePars.outOpts.overWriteFile_ = true;
			intersectBedLocsWtihGffRecords(refinedRegions, interPars);
			OutputStream bedOut(refinePars.outOpts);
			for(const auto & b : refinedRegions){
				bedOut << b->toDelimStrWithExtra() << std::endl;
			}
		}
		TwoBit::TwoBitFile tReader(genomeTwoBitFnp);
		auto refindedSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(chopPars.outputDirectory, "refined_merged.fasta"));
		SeqOutput refinedWriter(refindedSeqOpts);
		refinedWriter.openOut();
		for(const auto & reg : refinedGRegions){
			refinedWriter.write(reg.extractSeq(tReader));
		}
	}

//	if(0 != minLength){
//		std::vector<GenomicRegion> filteredRefinedRegions;
//		{
//			auto bedAgain = getBed3s(refinePars.outOpts.outName());
//			OutOptions filteredOpts(njh::files::make_path(chopPars.outputDirectory, "filtered_refined_merged.bed"));
//			OutputStream filteredOut(filteredOpts);
//			for(const auto & b : bedAgain){
//				if(b->length() >= minLength){
//					filteredRefinedRegions.emplace_back(*b);
//					filteredOut << b->toDelimStrWithExtra() << std::endl;
//				}
//			}
//		}
//		if(bfs::exists(genomeTwoBitFnp)){
//			TwoBit::TwoBitFile tReader(genomeTwoBitFnp);
//			auto refindedSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(chopPars.outputDirectory, "filtered_refined_merged.fasta"));
//			SeqOutput refinedWriter(refindedSeqOpts);
//			refinedWriter.openOut();
//			for(const auto & reg : filteredRefinedRegions){
//				refinedWriter.write(reg.extractSeq(tReader));
//			}
//		}
//	}
	if(0 != minLength || (entropyCutOff > 0 && bfs::exists(genomeTwoBitFnp))){
		std::vector<GenomicRegion> filteredRefinedRegions;
		std::shared_ptr<TwoBit::TwoBitFile> tReaderPrt;
		if(bfs::exists(genomeTwoBitFnp)){
			tReaderPrt = std::make_shared<TwoBit::TwoBitFile>(genomeTwoBitFnp);
		}
		{
			auto bedAgain = getBed3s(refinePars.outOpts.outName());
			OutOptions filteredOpts(njh::files::make_path(chopPars.outputDirectory, "filtered_refined_merged.bed"));
			OutputStream filteredOut(filteredOpts);
			for(const auto & b : bedAgain){
				if(b->length() >= minLength){
					bool furtherFilt = false;
					if(nullptr != tReaderPrt &&  entropyCutOff > 0){
						auto eSeq = GenomicRegion(*b).extractSeq(*tReaderPrt);
						charCounter counter(eSeq.seq_);
						if(counter.computeEntrophy() < entropyCutOff){
							furtherFilt = true;
						}
					}
					if(!furtherFilt){
						filteredRefinedRegions.emplace_back(*b);
						filteredOut << b->toDelimStrWithExtra() << std::endl;
					}
				}
			}
		}
		if(bfs::exists(genomeTwoBitFnp)){
			TwoBit::TwoBitFile tReader(genomeTwoBitFnp);
			auto refindedSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(chopPars.outputDirectory, "filtered_refined_merged.fasta"));
			SeqOutput refinedWriter(refindedSeqOpts);
			refinedWriter.openOut();
			for(const auto & reg : filteredRefinedRegions){
				refinedWriter.write(reg.extractSeq(tReader));
			}
		}
	}

	return 0;
}


}  // namespace njhseq

