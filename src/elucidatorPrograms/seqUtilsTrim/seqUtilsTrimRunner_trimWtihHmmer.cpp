/*
 * seqUtilsTrimRunner_trimWtihHmmer.cpp
 *
 *  Created on: May 26, 2021
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

#include "seqUtilsTrimRunner.hpp"
#include <njhseq/IO/SeqIO/SeqIO.hpp>
#include <njhseq/readVectorManipulation/readVectorHelpers/readVecTrimmer.hpp>
#include <njhseq/system/Muscler.hpp>
#include <njhseq/objects/BioDataObject/BioDataFileIO.hpp>

#include <njhseq/objects/BioDataObject/HmmerDomainHitTab.hpp>

namespace njhseq {


seqInfo trimSeqToHmmDomainHit(const seqInfo &seq, const HmmerDomainHitTab &dom,
		bool mark) {
	auto subSeq = seq.getSubRead(dom.envFrom_ - 1,
			(dom.envTo_ - dom.envFrom_ + 1));
	if (mark) {
		MetaDataInName meta;
		if (MetaDataInName::nameHasMetaData(subSeq.name_)) {
			meta = MetaDataInName(subSeq.name_);
		}
		meta.addMeta("hmmFrom", dom.hmmFrom_ - 1, true);
		meta.addMeta("hmmTo", dom.hmmTo_, true);
		meta.addMeta("trimStart", dom.envFrom_ - 1, true);
		meta.addMeta("trimEnd", dom.envTo_, true);
		meta.addMeta("queryName", dom.queryName_, true);
		meta.resetMetaInName(subSeq.name_);
	}

	return subSeq;
}



int seqUtilsTrimRunner::trimWithHmmer(const njh::progutils::CmdArgs & inputCommands) {

	bool useAccuracyForBest = false;
	bool onlyBestDomainHit = false;
	std::string defaultParameters = "--nonull2 --incT 10 --incdomT 10 -T 10 --domT 10";
	bfs::path hmmModel;
	bool mark = false;
	seqUtilsTrimSetUp setUp(inputCommands);
	FullTrimReadsPars pars;
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processIoOptions(VecStr{"fasta"});
	setUp.setOption(pars.keepOnlyOn, "--keepOnlyOn", "Keep Only the Reads that are still on");
	setUp.setOption(onlyBestDomainHit, "--onlyBestDomainHit", "Only trim to best domain hit, if ties will do the first domain found");
	setUp.setOption(useAccuracyForBest, "--useAccuracyForBest", "When only trimming to best domain hit, use best accuracy score rather than best bit score");
	setUp.setOption(defaultParameters, "--defaultParameters", "The default parameters given to hmmsearch");
	setUp.setOption(hmmModel, "--hmmModel", "hmm model database, created by hmmbuild", true);
	setUp.setOption(mark, "--mark", "mark the sequence name with the trim details");


	setUp.finishSetUp(std::cout);

	njh::sys::requireExternalProgramThrow("hmmsearch");


	bfs::path dmtbloutOutfnp = setUp.pars_.ioOptions_.out_.outFilename_.string() + "_domain_hits.txt";
	bfs::path hmmsearchOutfnp = setUp.pars_.ioOptions_.out_.outFilename_.string() + "_hmmsearch_output.txt";

	if(bfs::exists(dmtbloutOutfnp) && !setUp.pars_.ioOptions_.out_.overWriteFile_){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << dmtbloutOutfnp << " exists, use --overWrite to over it"<< "\n";
		throw std::runtime_error{ss.str()};
	}
	if(bfs::exists(hmmsearchOutfnp) && !setUp.pars_.ioOptions_.out_.overWriteFile_){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << hmmsearchOutfnp << " exists, use --overWrite to over it"<< "\n";
		throw std::runtime_error{ss.str()};
	}

	std::stringstream cmdSs;
	cmdSs << "hmmsearch " << defaultParameters
			<< " " << "--domtblout " << dmtbloutOutfnp
			<< " " << hmmModel
			<< " " << setUp.pars_.ioOptions_.firstName_
			<< " " << " > " << hmmsearchOutfnp;
	auto cmdOutput = njh::sys::run({cmdSs.str()});
	if(!cmdOutput.success_){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "failed to run hmmsearch " << "\n";
		ss << cmdOutput.toJson() << "\n";
		throw std::runtime_error{ss.str()};
	}

	std::unordered_map<std::string, std::vector<HmmerDomainHitTab>> domainHits;

	{
		BioDataFileIO<HmmerDomainHitTab> reader{IoOptions(InOptions(dmtbloutOutfnp))};
		HmmerDomainHitTab domain;
		reader.openIn();
		while(reader.readNextRecord(domain)){
			domainHits[domain.targetName_].emplace_back(domain);
		}
	}

	SeqIO seqReader(setUp.pars_.ioOptions_);
	seqReader.openIn();
	seqReader.openOut();
	seqInfo seq;
	VecStr noDomainHits;

	while(seqReader.readNextRead(seq)){
		if(!njh::in(seq.name_, domainHits)){
			seq.on_ = false;
			noDomainHits.emplace_back(seq.name_);
			if(!pars.keepOnlyOn || seq.on_){
				seqReader.write(seq);
			}
		} else {
			std::vector<HmmerDomainHitTab> hits;
			if(!onlyBestDomainHit){
				hits = domainHits[seq.name_];
			}else{
				if(useAccuracyForBest){
					double bestScore = 0;
					for(const auto & dom : domainHits.at(seq.name_)){
						if(dom.acc_ > bestScore){
							bestScore = dom.acc_;
							hits.clear();
							hits.emplace_back(dom);
						}
					}
				}else{
					double bestScore = std::numeric_limits<uint32_t>::max();
					for(const auto & dom : domainHits.at(seq.name_)){
						if(dom.domain_c_evalue_ < bestScore){
							bestScore = dom.domain_c_evalue_;
							hits.clear();
							hits.emplace_back(dom);
						}
					}
				}
			}
			//remember domain hits 1 based sequence positioning
			for(const auto & dom : hits){
				auto subSeq = trimSeqToHmmDomainHit(seq, dom, mark);
				seqReader.write(subSeq);
			}
		}
	}
	if(!noDomainHits.empty()){
		bfs::path noDomainHitsFnp = setUp.pars_.ioOptions_.out_.outFilename_.string() + "_noDomainHits.txt";
		OutOptions noDomainHitsOpts(noDomainHitsFnp);
		noDomainHitsOpts.transferOverwriteOpts(setUp.pars_.ioOptions_.out_);
		OutputStream noDomainHitsOut(noDomainHitsOpts);
		noDomainHitsOut << njh::conToStr(noDomainHits, "\n") << std::endl;
	}
	return 0;
}


} // namespace njhseq

