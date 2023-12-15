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
#include <njhseq/objects/BioDataObject/BioDataFileIO.hpp>
#include <njhseq/objects/BioDataObject/BioRecordsUtils/HmmerUtility.hpp>
#include <njhseq/objects/BioDataObject/BioRecordsUtils/BedUtility.hpp>
#include <njhseq/objects/BioDataObject/reading.hpp>

namespace njhseq {

int seqUtilsTrimRunner::trimWithnhmmscan(const njh::progutils::CmdArgs & inputCommands) {

	bool onlyBestDomainHit = false;
	std::string defaultParameters = "--nonull2 --incT 10 -T 10 --notextw";
	bfs::path hmmModel;
	bool mark = false;
	seqUtilsTrimSetUp setUp(inputCommands);
	FullTrimReadsPars pars;
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processIoOptions(VecStr{"fasta"});
	setUp.setOption(pars.keepOnlyOn, "--keepOnlyOn", "Keep Only the Reads that are still on");
	setUp.setOption(onlyBestDomainHit, "--onlyBestDomainHit", "Only trim to best domain hit, if ties will do the first domain found");
	setUp.setOption(defaultParameters, "--defaultParameters", "The default parameters given to hmmsearch");
	setUp.setOption(hmmModel, "--hmmModel", "hmm model database, created by hmmbuild", true);
	setUp.setOption(mark, "--mark", "mark the sequence name with the trim details");


	setUp.finishSetUp(std::cout);

	njh::sys::requireExternalProgramThrow("nhmmscan");


	bfs::path hitTabOutfnp = setUp.pars_.ioOptions_.out_.outFilename_.string() + "_hits.txt";
	bfs::path hmmsearchOutfnp = setUp.pars_.ioOptions_.out_.outFilename_.string() + "_hmmsearch_output.txt";

	if(bfs::exists(hitTabOutfnp) && !setUp.pars_.ioOptions_.out_.overWriteFile_){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << hitTabOutfnp << " exists, use --overWrite to over it"<< "\n";
		throw std::runtime_error{ss.str()};
	}
	if(bfs::exists(hmmsearchOutfnp) && !setUp.pars_.ioOptions_.out_.overWriteFile_){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << hmmsearchOutfnp << " exists, use --overWrite to over it"<< "\n";
		throw std::runtime_error{ss.str()};
	}

	std::stringstream cmdSs;
	cmdSs << "nhmmscan " << defaultParameters
			<< " " << "--tblout " << hitTabOutfnp
			<< " " << hmmModel
			<< " " << setUp.pars_.ioOptions_.firstName_
			<< " " << " > " << hmmsearchOutfnp;
	auto cmdOutput = njh::sys::run({cmdSs.str()});
	if(!cmdOutput.success_){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "failed to run nhmmscan " << "\n";
		ss << cmdOutput.toJson() << "\n";
		throw std::runtime_error{ss.str()};
	}

	std::unordered_map<std::string, std::vector<HmmerSeqHitTab>> hits;

	{
		BioDataFileIO<HmmerSeqHitTab> reader{IoOptions(InOptions(hitTabOutfnp))};
		HmmerSeqHitTab domain;
		reader.openIn();
		while(reader.readNextRecord(domain)){
			hits[domain.queryName_].emplace_back(domain);
		}
	}

	SeqIO seqReader(setUp.pars_.ioOptions_);
	seqReader.openIn();
	seqReader.openOut();
	seqInfo seq;
	VecStr noDomainHits;

	while(seqReader.readNextRead(seq)){
		if(!njh::in(seq.name_, hits)){
			seq.on_ = false;
			noDomainHits.emplace_back(seq.name_);
			if(!pars.keepOnlyOn || seq.on_){
				seqReader.write(seq);
			}
		} else {
			std::vector<HmmerSeqHitTab> trimHits;
			if(!onlyBestDomainHit){
				trimHits = hits[seq.name_];
			} else {
				//get best nonoverlapping hits
				HmmerUtility::sortHmmSeqHitsByEvalue(hits.at(seq.name_));
				trimHits = HmmerUtility::getNonoverlapingSeqHitsPostSortEnv(hits.at(seq.name_));
			}
			for(const auto & dom : trimHits){
				auto subSeq = HmmerUtility::trimSeqToHmmSeqHitEnv(seq, dom, mark);
				seqReader.write(subSeq);
			}
		}
	}
	if(!noDomainHits.empty()){
		bfs::path noDomainHitsFnp = setUp.pars_.ioOptions_.out_.outFilename_.string() + "_noHits.txt";
		OutOptions noDomainHitsOpts(noDomainHitsFnp);
		noDomainHitsOpts.transferOverwriteOpts(setUp.pars_.ioOptions_.out_);
		OutputStream noDomainHitsOut(noDomainHitsOpts);
		noDomainHitsOut << njh::conToStr(noDomainHits, "\n") << std::endl;
	}
	return 0;
}


int seqUtilsTrimRunner::trimWithhmmsearch(const njh::progutils::CmdArgs & inputCommands) {

	bool useAccuracyForBest = false;
	bool onlyBestDomainHit = false;
	std::string defaultParameters = "--nonull2 --incT 10 --incdomT 10 -T 10 --domT 10 --notextw";
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
					HmmerUtility::sortHmmDomainHitsByAccuracy(domainHits.at(seq.name_));
				}else{
					HmmerUtility::sortHmmDomainHitsByEvalue(domainHits.at(seq.name_));
				}
				hits = HmmerUtility::getNonoverlapingDomainHitsPostSortEnv(domainHits.at(seq.name_));
			}
			for(const auto & dom : hits){
				auto subSeq = HmmerUtility::trimSeqToHmmDomainHitEnv(seq, dom, mark);
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

int seqUtilsTrimRunner::trimBetweenHmmViaHmmsearch(const njh::progutils::CmdArgs & inputCommands) {

	std::string defaultParameters = "--nonull2 --incT 10 --incdomT 10 -T 10 --domT 10";
	bfs::path hmmModelStart;
	bfs::path hmmModelEnd;
	uint32_t hmmStartFilter = 5;
	bool excludeHmmHitRegion = false;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processDefaultReader(VecStr{"fasta", "fastagz", "fastq", "fastqgz"}, true);
	setUp.setOption(defaultParameters, "--defaultParameters", "The default parameters given to hmmsearch");
	setUp.setOption(hmmModelEnd, "--hmmModelEnd", "hmm model to trim to, created by hmmbuild", true);
	setUp.setOption(hmmModelStart, "--hmmModelStart", "hmm model to trim from, created by hmmbuild", true);
	setUp.setOption(hmmStartFilter, "--hmmStartFilter", "Filter partial hmms domain hits if they start or end this far into the model");
	setUp.setOption(excludeHmmHitRegion, "--excludeHmmHitRegion", "exclude Hmm Hit Region");

	setUp.processDirectoryOutputName(true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	auto startDir = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar{"startHmm"});

	std::unordered_map<std::string, seqInfo> allSeqs;
	{
		seqInfo seq;
		SeqInput reader(setUp.pars_.ioOptions_);
		reader.openIn();

		while(reader.readNextRead(seq)){
			if(njh::in(seq.name_, allSeqs)){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "already have: " << seq.name_ << "\n";
				throw std::runtime_error{ss.str()};
			}
			allSeqs[seq.name_] = seq;
		}
	}

	{
		//start

		//convert to fasta and non-gz hmm model
		std::unordered_map<uint32_t, std::string> seqKey;
		auto inputSeqFnp = njh::files::make_path(startDir, "inputSeqs.fasta");
		auto inputSeqSubDomainsFnp = njh::files::make_path(startDir, "inputSeqsSubDomains.fasta");

		auto hmmModelFnp = njh::files::make_path(startDir, "hmmModel.txt");

		{
			seqInfo seq;
			SeqInput reader(setUp.pars_.ioOptions_);
			SeqOutput writer(SeqIOOptions::genFastaOut(inputSeqFnp));
			OutputStream seqNameKeyOut(njh::files::make_path(startDir, "inputSeqNameKey.tab.txt"));
			seqNameKeyOut << "oldName\tnewName" << std::endl;
			reader.openIn();
			writer.openOut();

			uint32_t count = 0;
			while(reader.readNextRead(seq)){
				seqNameKeyOut << seq.name_ << "\t" << count << std::endl;
				seqKey[count] = seq.name_;
				seq.name_ = njh::pasteAsStr(count);
				writer.write(seq);
				++count;
			}
		}
		{
			//write out the model, this way you can supply it in gzip format but this will unzip it
			InputStream inModel(hmmModelStart);
			OutputStream outModel(hmmModelFnp);
			outModel << inModel.rdbuf();
			outModel.flush();
		}

		std::stringstream cmdSs;
		cmdSs << "hmmsearch " << defaultParameters
					<< " " << "--domtblout " << "raw_all_domain_hits_table.txt"
					<< " " << "hmmModel.txt"
					<< " " << "inputSeqs.fasta"
					<< " " << " > " << "hmmsearch_raw_output.txt";
		std::string cdCmd = "cd " + startDir.string() + " && ";
		auto cmdOutput = njh::sys::run({cdCmd, cmdSs.str()});
		if(!cmdOutput.success_){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "failed to run hmmsearch " << "\n";
			ss << cmdOutput.toJson() << "\n";
			throw std::runtime_error{ss.str()};
		}else{
			OutputStream hmmsearchOutput(njh::files::make_path(startDir, "hmmsearchCmdRunDetails.json"));
			hmmsearchOutput << cmdOutput.toJson() << std::endl;
		}

		//convert domain hits table into a real table and
		//filter table
		auto rawDomainHitsFnp = njh::files::make_path(startDir, "raw_all_domain_hits_table.txt");
		auto allDomainHitsFnp = njh::files::make_path(startDir, "all_domain_hits_table.tab.txt");
		auto filtDomainHitsFnp = njh::files::make_path(startDir, "filt_domain_hits_table.tab.txt");
		auto nonOverlappingDomainHitsFnp = njh::files::make_path(startDir, "nonOverlappingBestdomain_hits_table.tab.txt");
		auto bedRegionsFnp = njh::files::make_path(startDir, "nonOverlappingBestdomain_hits.bed");
		auto seqsWithNoDomainHitsFnp = njh::files::make_path(startDir, "seqsWithNoDomainHits.tab.txt");
		auto domainCountsFnp = njh::files::make_path(startDir, "seqDomainCounts.tab.txt");

		{
			OutputStream out(allDomainHitsFnp);
			OutputStream outFilt(filtDomainHitsFnp);

			BioDataFileIO<HmmerDomainHitTab> reader{IoOptions(InOptions(rawDomainHitsFnp))};
			HmmerDomainHitTab domain;
			reader.openIn();
			// uint32_t count = 0;
			out << "#" << njh::conToStr(HmmerDomainHitTab::toDelimStrHeader(), "\t") << std::endl;
			outFilt << "#" << njh::conToStr(HmmerDomainHitTab::toDelimStrHeader(), "\t") << std::endl;

			while(reader.readNextRecord(domain)){
				out << domain.toDelimStr() << std::endl;
				bool passHmmStart = domain.zeroBasedHmmFrom() <= hmmStartFilter;
				bool passHmmEnd = (domain.queryLen_ - domain.hmmTo_) <= hmmStartFilter;
				if(passHmmStart && passHmmEnd) {
					outFilt << domain.toDelimStr() << std::endl;
				}
				// ++count;
			}
		}
		//get best non-overlapping positions
		std::vector<Bed6RecordCore> filteredRegions;

		std::unordered_map<std::string, std::vector<HmmerDomainHitTab>> domainsPerSeq;
		std::map<uint32_t, std::unordered_map<std::string, uint32_t>> domainCountsPerSeq;

		{
			std::vector<HmmerDomainHitTab> domains;
			std::vector<Bed6RecordCore> locations;

			{
				BioDataFileIO<HmmerDomainHitTab> reader{IoOptions(InOptions(filtDomainHitsFnp))};

				HmmerDomainHitTab domain;
				reader.openIn();
				uint32_t count = 0;
				while(reader.readNextRecord(domain)){
					domains.emplace_back(domain);
					locations.emplace_back(Bed6RecordCore(domain.targetName_, domain.envFrom_ -1, domain.envTo_, njh::pasteAsStr(count), domain.domain_c_evalue_, '+'));
					++count;
				}
			}
			njh::sort(locations,[&domains](const Bed6RecordCore & rec1, const Bed6RecordCore & rec2){
				if(rec1.score_ == rec2.score_){
					if(domains[njh::StrToNumConverter::stoToNum<uint32_t>(rec1.name_)].domain_i_evalue_ == domains[njh::StrToNumConverter::stoToNum<uint32_t>(rec2.name_)].domain_i_evalue_){
						return rec1.length() > rec2.length();
					}else{
						return domains[njh::StrToNumConverter::stoToNum<uint32_t>(rec1.name_)].domain_i_evalue_ < domains[njh::StrToNumConverter::stoToNum<uint32_t>(rec2.name_)].domain_i_evalue_;
					}
				} else {
					return rec1.score_ < rec2.score_;
				}
			});
			for(const auto & region : locations){
				bool overlap = false;
				for(const auto & outRegion : filteredRegions){
					if(region.overlaps(outRegion,1)){
						overlap = true;
						break;
					}
				}
				if(!overlap){
					filteredRegions.emplace_back(region);
				}
			}
			OutputStream outFiltNonOverlap(nonOverlappingDomainHitsFnp);
			outFiltNonOverlap << "#" << njh::conToStr(HmmerDomainHitTab::toDelimStrHeader(), "\t") << "\tbedName" << std::endl;
			BedUtility::coordSort(filteredRegions);
			OutputStream bedRegionsOut(bedRegionsFnp);
			for(auto & region : filteredRegions){
				const auto & domForReg  = domains[njh::StrToNumConverter::stoToNum<uint32_t>(region.name_)];
				auto bedNewName = njh::pasteAsStr(domForReg.queryName_, ".", domForReg.domainId_);
				outFiltNonOverlap << domForReg.toDelimStr() << "\t" << bedNewName<< std::endl;
				region.name_ = bedNewName;
				region.chrom_ = seqKey[njh::StrToNumConverter::stoToNum<uint32_t>(region.chrom_)];
				bedRegionsOut << region.toDelimStrWithExtra() << std::endl;
				domainsPerSeq[domForReg.targetName_].emplace_back(domForReg);
				++domainCountsPerSeq[njh::StrToNumConverter::stoToNum<uint32_t>(domForReg.targetName_)][domForReg.queryName_];
			}
		}

		//trim seqs to best overlapping positions sub seqs
		{
			VecStr seqsWithNoDomains;
			seqInfo seq;
			SeqInput reader(SeqIOOptions::genFastaIn(inputSeqFnp));
			SeqOutput writer(SeqIOOptions::genFastaOut(inputSeqSubDomainsFnp));
			reader.openIn();
			writer.openOut();

			while(reader.readNextRead(seq)){
				if(njh::in(seq.name_, domainsPerSeq)){
					for(const auto & domain : domainsPerSeq[seq.name_]){
						Bed6RecordCore region(domain.targetName_, domain.env0BasedPlusStrandStart(), domain.envTo_,njh::pasteAsStr(domain.queryName_, ".", domain.domainId_) , domain.domain_c_evalue_, '+');
						auto subSeq = seq.getSubRead(region.chromStart_, region.length());
						auto oldName = seqKey[njh::StrToNumConverter::stoToNum<uint32_t>(subSeq.name_)];
						MetaDataInName meta;
						if(MetaDataInName::nameHasMetaData(oldName)){
							meta = MetaDataInName(oldName);
						}
						subSeq.name_ = oldName;
						meta.addMeta("hmmFrom", domain.zeroBasedHmmFrom(), true);
						meta.addMeta("hmmTo", domain.hmmTo_, true);

						meta.addMeta("trimStart", domain.env0BasedPlusStrandStart(), true);
						meta.addMeta("trimEnd", domain.envTo_, true);

						meta.addMeta("queryName", domain.queryName_, true);
						meta.addMeta("queryID", region.name_, true);

						meta.resetMetaInName(subSeq.name_);
						writer.write(subSeq);
					}
				}else{
					seqsWithNoDomains.emplace_back(seq.name_);
				}
			}
			OutputStream noDomainHits(seqsWithNoDomainHitsFnp);
			noDomainHits << "newName\toldName" << std::endl;
			for(const auto & name : seqsWithNoDomains){
				noDomainHits << name
										 << "\t" << seqKey[njh::StrToNumConverter::stoToNum<uint32_t>(name)] << std::endl;
			}
			OutputStream domainHitsCouintsOut(domainCountsFnp);
			domainHitsCouintsOut << "newName\toldName\tdomain\tcount" << std::endl;
			for(const auto & name : domainCountsPerSeq){
				for(const auto & count : name.second){
					domainHitsCouintsOut
									<< name.first
									<< "\t" << seqKey[name.first]
									<< "\t" << count.first
									<< "\t" << count.second << std::endl;
				}
			}
		}
	}
	auto endDir = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar{"endHmm"});

	{
		//end
		//convert to fasta and non-gz hmm model
		std::unordered_map<uint32_t, std::string> seqKey;
		auto inputSeqFnp = njh::files::make_path(endDir, "inputSeqs.fasta");
		auto inputSeqSubDomainsFnp = njh::files::make_path(endDir, "inputSeqsSubDomains.fasta");

		auto hmmModelFnp = njh::files::make_path(endDir, "hmmModel.txt");

		{
			seqInfo seq;
			SeqInput reader(setUp.pars_.ioOptions_);
			SeqOutput writer(SeqIOOptions::genFastaOut(inputSeqFnp));
			OutputStream seqNameKeyOut(njh::files::make_path(endDir, "inputSeqNameKey.tab.txt"));
			seqNameKeyOut << "oldName\tnewName" << std::endl;
			reader.openIn();
			writer.openOut();

			uint32_t count = 0;
			while(reader.readNextRead(seq)){
				seqNameKeyOut << seq.name_ << "\t" << count << std::endl;
				seqKey[count] = seq.name_;
				seq.name_ = njh::pasteAsStr(count);
				writer.write(seq);
				++count;
			}
		}
		{
			//write out the model, this way you can supply it in gzip format but this will unzip it
			InputStream inModel(hmmModelEnd);
			OutputStream outModel(hmmModelFnp);
			outModel << inModel.rdbuf();
			outModel.flush();
		}

		std::stringstream cmdSs;
		cmdSs << "hmmsearch " << defaultParameters
					<< " " << "--domtblout " << "raw_all_domain_hits_table.txt"
					<< " " << "hmmModel.txt"
					<< " " << "inputSeqs.fasta"
					<< " " << " > " << "hmmsearch_raw_output.txt";
		std::string cdCmd = "cd " + endDir.string() + " && ";
		auto cmdOutput = njh::sys::run({cdCmd, cmdSs.str()});
		if(!cmdOutput.success_){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "failed to run hmmsearch " << "\n";
			ss << cmdOutput.toJson() << "\n";
			throw std::runtime_error{ss.str()};
		}else{
			OutputStream hmmsearchOutput(njh::files::make_path(endDir, "hmmsearchCmdRunDetails.json"));
			hmmsearchOutput << cmdOutput.toJson() << std::endl;
		}

		//convert domain hits table into a real table and
		//filter table
		auto rawDomainHitsFnp = njh::files::make_path(endDir, "raw_all_domain_hits_table.txt");
		auto allDomainHitsFnp = njh::files::make_path(endDir, "all_domain_hits_table.tab.txt");
		auto filtDomainHitsFnp = njh::files::make_path(endDir, "filt_domain_hits_table.tab.txt");
		auto nonOverlappingDomainHitsFnp = njh::files::make_path(endDir, "nonOverlappingBestdomain_hits_table.tab.txt");
		auto bedRegionsFnp = njh::files::make_path(endDir, "nonOverlappingBestdomain_hits.bed");
		auto seqsWithNoDomainHitsFnp = njh::files::make_path(endDir, "seqsWithNoDomainHits.tab.txt");
		auto domainCountsFnp = njh::files::make_path(endDir, "seqDomainCounts.tab.txt");

		{
			OutputStream out(allDomainHitsFnp);
			OutputStream outFilt(filtDomainHitsFnp);

			BioDataFileIO<HmmerDomainHitTab> reader{IoOptions(InOptions(rawDomainHitsFnp))};
			HmmerDomainHitTab domain;
			reader.openIn();
			// uint32_t count = 0;
			out << "#" << njh::conToStr(HmmerDomainHitTab::toDelimStrHeader(), "\t") << std::endl;
			outFilt << "#" << njh::conToStr(HmmerDomainHitTab::toDelimStrHeader(), "\t") << std::endl;

			while(reader.readNextRecord(domain)){
				out << domain.toDelimStr() << std::endl;
				bool passHmmStart = domain.zeroBasedHmmFrom() <= hmmStartFilter;
				bool passHmmEnd = (domain.queryLen_ - domain.hmmTo_) <= hmmStartFilter;
				if(passHmmStart && passHmmEnd) {
					outFilt << domain.toDelimStr() << std::endl;
				}
				// ++count;
			}
		}
		//get best non-overlapping positions
		std::vector<Bed6RecordCore> filteredRegions;

		std::unordered_map<std::string, std::vector<HmmerDomainHitTab>> domainsPerSeq;
		std::map<uint32_t, std::unordered_map<std::string, uint32_t>> domainCountsPerSeq;

		{
			std::vector<HmmerDomainHitTab> domains;
			std::vector<Bed6RecordCore> locations;

			{
				BioDataFileIO<HmmerDomainHitTab> reader{IoOptions(InOptions(filtDomainHitsFnp))};

				HmmerDomainHitTab domain;
				reader.openIn();
				uint32_t count = 0;
				while(reader.readNextRecord(domain)){
					domains.emplace_back(domain);
					locations.emplace_back(Bed6RecordCore(domain.targetName_, domain.envFrom_ -1, domain.envTo_, njh::pasteAsStr(count), domain.domain_c_evalue_, '+'));
					++count;
				}
			}
			njh::sort(locations,[&domains](const Bed6RecordCore & rec1, const Bed6RecordCore & rec2){
				if(rec1.score_ == rec2.score_){
					if(domains[njh::StrToNumConverter::stoToNum<uint32_t>(rec1.name_)].domain_i_evalue_ == domains[njh::StrToNumConverter::stoToNum<uint32_t>(rec2.name_)].domain_i_evalue_){
						return rec1.length() > rec2.length();
					}else{
						return domains[njh::StrToNumConverter::stoToNum<uint32_t>(rec1.name_)].domain_i_evalue_ < domains[njh::StrToNumConverter::stoToNum<uint32_t>(rec2.name_)].domain_i_evalue_;
					}
				} else {
					return rec1.score_ < rec2.score_;
				}
			});
			for(const auto & region : locations){
				bool overlap = false;
				for(const auto & outRegion : filteredRegions){
					if(region.overlaps(outRegion,1)){
						overlap = true;
						break;
					}
				}
				if(!overlap){
					filteredRegions.emplace_back(region);
				}
			}
			OutputStream outFiltNonOverlap(nonOverlappingDomainHitsFnp);
			outFiltNonOverlap << "#" << njh::conToStr(HmmerDomainHitTab::toDelimStrHeader(), "\t") << "\tbedName" << std::endl;
			BedUtility::coordSort(filteredRegions);
			OutputStream bedRegionsOut(bedRegionsFnp);
			for(auto & region : filteredRegions){
				const auto & domForReg  = domains[njh::StrToNumConverter::stoToNum<uint32_t>(region.name_)];
				auto bedNewName = njh::pasteAsStr(domForReg.queryName_, ".", domForReg.domainId_);
				outFiltNonOverlap << domForReg.toDelimStr() << "\t" << bedNewName<< std::endl;
				region.name_ = bedNewName;
				region.chrom_ = seqKey[njh::StrToNumConverter::stoToNum<uint32_t>(region.chrom_)];
				bedRegionsOut << region.toDelimStrWithExtra() << std::endl;
				domainsPerSeq[domForReg.targetName_].emplace_back(domForReg);
				++domainCountsPerSeq[njh::StrToNumConverter::stoToNum<uint32_t>(domForReg.targetName_)][domForReg.queryName_];
			}
		}

		//trim seqs to best overlapping positions sub seqs
		{
			VecStr seqsWithNoDomains;
			seqInfo seq;
			SeqInput reader(SeqIOOptions::genFastaIn(inputSeqFnp));
			SeqOutput writer(SeqIOOptions::genFastaOut(inputSeqSubDomainsFnp));
			reader.openIn();
			writer.openOut();

			while(reader.readNextRead(seq)){
				if(njh::in(seq.name_, domainsPerSeq)){
					for(const auto & domain : domainsPerSeq[seq.name_]){
						Bed6RecordCore region(domain.targetName_, domain.env0BasedPlusStrandStart(), domain.envTo_,njh::pasteAsStr(domain.queryName_, ".", domain.domainId_) , domain.domain_c_evalue_, '+');
						auto subSeq = seq.getSubRead(region.chromStart_, region.length());
						auto oldName = seqKey[njh::StrToNumConverter::stoToNum<uint32_t>(subSeq.name_)];
						MetaDataInName meta;
						if(MetaDataInName::nameHasMetaData(oldName)){
							meta = MetaDataInName(oldName);
						}
						subSeq.name_ = oldName;
						meta.addMeta("hmmFrom", domain.zeroBasedHmmFrom(), true);
						meta.addMeta("hmmTo", domain.hmmTo_, true);

						meta.addMeta("trimStart", domain.env0BasedPlusStrandStart(), true);
						meta.addMeta("trimEnd", domain.envTo_, true);

						meta.addMeta("queryName", domain.queryName_, true);
						meta.addMeta("queryID", region.name_, true);

						meta.resetMetaInName(subSeq.name_);
						writer.write(subSeq);
					}
				}else{
					seqsWithNoDomains.emplace_back(seq.name_);
				}
			}
			OutputStream noDomainHits(seqsWithNoDomainHitsFnp);
			noDomainHits << "newName\toldName" << std::endl;
			for(const auto & name : seqsWithNoDomains){
				noDomainHits << name
										 << "\t" << seqKey[njh::StrToNumConverter::stoToNum<uint32_t>(name)] << std::endl;
			}
			OutputStream domainHitsCouintsOut(domainCountsFnp);
			domainHitsCouintsOut << "newName\toldName\tdomain\tcount" << std::endl;
			for(const auto & name : domainCountsPerSeq){
				for(const auto & count : name.second){
					domainHitsCouintsOut
									<< name.first
									<< "\t" << seqKey[name.first]
									<< "\t" << count.first
									<< "\t" << count.second << std::endl;
				}
			}
		}
	}

	auto startBedRegionsFnp = njh::files::make_path(startDir, "nonOverlappingBestdomain_hits.bed");
	auto endBedRegionsFnp = njh::files::make_path(endDir, "nonOverlappingBestdomain_hits.bed");

	std::unordered_map<std::string, std::vector<GenomicRegion>> startsBySeq;
	auto startRegions  = bedPtrsToGenomicRegs(getBeds(startBedRegionsFnp));
	for(const auto & start : startRegions){
		startsBySeq[start.chrom_].emplace_back(start);
	}
	std::unordered_map<std::string, std::vector<GenomicRegion>> endsBySeq;
	auto endRegions  = bedPtrsToGenomicRegs(getBeds(endBedRegionsFnp));
	for(const auto & end : endRegions){
		endsBySeq[end.chrom_].emplace_back(end);
	}
	SeqIOOptions trimmedOutOpts(njh::files::make_path(setUp.pars_.directoryName_, "trimmed"), setUp.pars_.ioOptions_.outFormat_);
	SeqOutput seqWriter(trimmedOutOpts);
	seqWriter.openOut();
	for(const auto & seq : allSeqs){
		if(njh::in(seq.second.name_, startsBySeq) && njh::in(seq.second.name_, endsBySeq)){
			for(const auto & start : startsBySeq[seq.second.name_]){
				for(const auto & end : endsBySeq[seq.second.name_]){
					uint32_t startPos = start.start_;
					uint32_t endPos = end.end_;
					if(excludeHmmHitRegion){
						startPos = start.end_;
						endPos = end.start_;
					}
					if(startPos < endPos){
						auto subSeq = seq.second.getSubRead(startPos, endPos - startPos);
						seqWriter.write(subSeq);
					} else {
						//no trim
					}
				}
			}
		} else {
			//no trim
		}
	}

	return 0;
}

int seqUtilsTrimRunner::trimBetweenHmmViaNhmmscan(const njh::progutils::CmdArgs & inputCommands) {

	bfs::path subRegions;
	uint32_t extendSubRegions = 0;
	std::string defaultParameters = "--nonull2 --incT 50 --incdomT 50 -T 50 --notextw";
	bfs::path hmmModelStart;
	bfs::path hmmModelEnd;
	bool excludeHmmHitRegion = false;
	nhmmscanOutput::PostProcessHitsPars postProcessPars;
	postProcessPars.accCutOff = 0.80;
	postProcessPars.scoreCutOff = 200;
	postProcessPars.evalueCutOff = 1e-100;
	postProcessPars.scoreNormCutOff = 0.50;
	nhmmscanOutput::QueryResults::mergeOverlapingHitsPars mergePars;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processDefaultReader(VecStr{"fasta", "fastagz", "fastq", "fastqgz"}, true);
	setUp.setOption(defaultParameters, "--defaultParameters", "The default parameters given to hmmsearch");
	setUp.setOption(extendSubRegions, "--extendSubRegions", "Extend the sub regions by this much in both directions");
	setUp.setOption(excludeHmmHitRegion, "--excludeHmmHitRegion", "exclude Hmm Hit Region");

	setUp.setOption(hmmModelStart, "--hmmModelStart", "hmm model to trim from, created by hmmbuild", true);
	setUp.setOption(hmmModelEnd, "--hmmModelEnd", "hmm model to trim from, created by hmmbuild", true);
	setUp.setOption(postProcessPars.hmmStartFilter, "--hmmStartFilter", "Filter partial hmms domain hits if they start or end this far into the model");
	setUp.setOption(postProcessPars.minLength, "--minLength", "Minimum output domain hit length");
	setUp.setOption(subRegions, "--subRegions", "Run on the subregions as defined by this bed file, needs to have a 2bit file named with same prefix as the input file");

	setUp.setOption(postProcessPars.accCutOff, "--accCutOff", "soft accuracy cut off");
	setUp.setOption(postProcessPars.scoreCutOff, "--scoreCutOff", "soft score cut off");
	setUp.setOption(postProcessPars.evalueCutOff, "--evalueCutOff", "soft evalue cut off");
	setUp.setOption(postProcessPars.scoreNormCutOff, "--scoreNormCutOff", "soft scoreNorm cut off");

	setUp.setOption(postProcessPars.hardAccCutOff, "--hardAccCutOff", "hard accuracy cut off");
	if(postProcessPars.hardAccCutOff > postProcessPars.accCutOff){
		postProcessPars.accCutOff = postProcessPars.hardAccCutOff;
	}
	setUp.setOption(postProcessPars.hardScoreCutOff, "--hardScoreCutOff", "hard score cut off");
	if(postProcessPars.hardScoreCutOff > postProcessPars.scoreCutOff){
		postProcessPars.scoreCutOff = postProcessPars.hardScoreCutOff;
	}
	setUp.setOption(postProcessPars.hardEvalueCutOff, "--hardEvalueCutOff", "hard evalue cut off");
	if(postProcessPars.hardEvalueCutOff > postProcessPars.evalueCutOff){
		postProcessPars.evalueCutOff = postProcessPars.hardEvalueCutOff;
	}
	setUp.setOption(postProcessPars.hardScoreNormCutOff, "--hardScoreNormCutOff", "hard scoreNorm cut off");
	if(postProcessPars.hardScoreNormCutOff > postProcessPars.scoreNormCutOff){
		postProcessPars.scoreNormCutOff = postProcessPars.hardScoreNormCutOff;
	}
	setUp.processDirectoryOutputName(true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);


	auto startDir = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar{"startHmm"});
	auto endDir = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar{"endHmm"});

	std::unordered_map<std::string, seqInfo> allSeqs;
	{
		seqInfo seq;
		SeqInput reader(setUp.pars_.ioOptions_);
		reader.openIn();

		while(reader.readNextRead(seq)){
			if(njh::in(seq.name_, allSeqs)){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "already have: " << seq.name_ << "\n";
				throw std::runtime_error{ss.str()};
			}
			allSeqs[seq.name_] = seq;
		}
	}

	{
		//start
		//convert to fasta and non-gz hmm model

		auto inputSeqFnp = njh::files::make_path(startDir, "inputSeqs.fasta");
		auto noOverlapFiltFnp = njh::files::make_path(startDir, "noOverlapFiltHits.fasta");
		auto noOverlapMergedFiltFnp = njh::files::make_path(startDir, "noOverlapMergedFiltHits.fasta");
		auto hmmModelFnp = njh::files::make_path(startDir, "hmmModel.txt");

		std::unordered_map<std::string, GenomicRegion> regionsByName;
		std::unordered_map<std::string, uint32_t> realQueryLens;
		std::vector<GenomicRegion> regions;
		auto inputOpts = setUp.pars_.ioOptions_;
		if(!subRegions.empty()){
			auto twoBitFnp = bfs::path(setUp.pars_.ioOptions_.firstName_).replace_extension(".2bit");
			TwoBit::TwoBitFile tReader(twoBitFnp);
			realQueryLens = tReader.getSeqLens();
			auto regionsOutOpts = SeqIOOptions::genFastaOut(njh::files::make_path(startDir, "subRegions.fasta"));
			regions = bed3PtrsToGenomicRegs(getBed3s(subRegions));
			if(0 != extendSubRegions){
				njh::for_each(regions, [&realQueryLens,&extendSubRegions](GenomicRegion & region){
					BedUtility::extendLeftRight(region, extendSubRegions, extendSubRegions, njh::mapAt(realQueryLens, region.chrom_));
				});
				sortGRegionsByStart(regions);
				std::vector<GenomicRegion> mergedRegions;
				for(const auto & reg : regions){
					bool overlaps = false;
					if(!mergedRegions.empty()){
						if(mergedRegions.back().overlaps(reg, 1)) {
							overlaps = true;
							mergedRegions.back().end_ = std::max(mergedRegions.back().end_, reg.end_);
						}
					}
					if(!overlaps){
						mergedRegions.emplace_back(reg);
					}
				}
				OutputStream mergedBedOut(njh::files::make_path(startDir, "extendedSubRegions.bed"));
				njh::for_each(mergedRegions,[&mergedBedOut](GenomicRegion & region){
					region.uid_ = region.createUidFromCoords();
					mergedBedOut << region.genBedRecordCore().toDelimStr() << std::endl;
				});
				regions = mergedRegions;
			}


			SeqOutput writer(regionsOutOpts);
			writer.openOut();
			// uint32_t regionCount = 0;
			for(const auto & region : regions){
				regionsByName[region.uid_] = region;
				// ++regionCount;
				auto subRegion = region.extractSeq(tReader);
				subRegion.name_ = region.uid_;
				writer.write(subRegion);
			}
			inputOpts = SeqIOOptions::genFastaIn(regionsOutOpts.out_.outName());
		}

		auto seqKey = SeqIO::rewriteSeqsWithIndexAsName(
						inputOpts,
						SeqIOOptions::genFastaOut(inputSeqFnp),
						njh::files::make_path(startDir, "inputSeqNameKey.tab.txt"));

		//write out the model, this way you can supply it in gzip format but this will unzip it
		njh::files::reWriteFile(hmmModelStart, hmmModelFnp);
		nhmmscanOutput::run_hmmpress_ifNeed(hmmModelFnp);

		std::stringstream cmdSs;
		cmdSs << "nhmmscan " << defaultParameters
					<< " " << "--tblout " << "raw_all_hits_table.txt"
					<< " " << "hmmModel.txt"
					<< " " << "inputSeqs.fasta"
					<< " " << " > " << "nhmmscan_raw_output.txt";
		std::string cdCmd = "cd " + startDir.string() + " && ";
		auto cmdOutput = njh::sys::run( { cdCmd, cmdSs.str() });
		if (!cmdOutput.success_) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "failed to run hmmsearch "
				 << "\n";
			ss << cmdOutput.toJson() << "\n";
			throw std::runtime_error { ss.str() };
		} else {
			OutputStream nhmmscanOutput(
							njh::files::make_path(startDir,
																		"nhmmscanCmdRunDetails.json"));
			nhmmscanOutput << cmdOutput.toJson() << std::endl;
		}



		auto nhmmscan_raw_outputFnp = njh::files::make_path(startDir, "nhmmscan_raw_output.txt");
		auto seqsWithNoDomainHitsFnp = njh::files::make_path(startDir, "seqsWithNoDomainHits.tab.txt");
		nhmmscanOutput outputParsed;
		if(regions.empty()){
			outputParsed = nhmmscanOutput::parseRawOutput(nhmmscan_raw_outputFnp, seqKey);
		}else{
			outputParsed = nhmmscanOutput::parseRawOutput(nhmmscan_raw_outputFnp, seqKey, regionsByName, realQueryLens);
		}
		auto postProcessResults = outputParsed.postProcessHits(postProcessPars);
		//convert hits table into a real table and
		outputParsed.writeInfoFiles(postProcessResults, startDir);

		if(setUp.pars_.debug_){
			//merging debugging
			for(const auto & filteredHits : postProcessResults.filteredHitsByQuery_){
				auto mergedResults = nhmmscanOutput::QueryResults::mergeOverlapingHits(filteredHits.second, mergePars);
				for(const auto & merged : mergedResults){
					if(merged.hits_.size() > 1){
						std::cout << merged.region_.genBedRecordCore().toDelimStrWithExtra() << std::endl;
						for(const auto & hit : merged.hits_){
							std::cout << "\t" << hit.genBed6_env().toDelimStrWithExtra() << std::endl;
						}
						std::cout << std::endl;
					}
				}
			}
		}

		//trim seqs to best overlapping positions sub seqs
		{
			VecStr seqsWithNoDomains;
			seqInfo seq;
			SeqInput reader(setUp.pars_.ioOptions_);
			SeqOutput noOverlapFiltWriter(SeqIOOptions::genFastaOut(noOverlapFiltFnp));
			SeqOutput noOverlapMergedFiltWriter(SeqIOOptions::genFastaOut(noOverlapMergedFiltFnp));

			reader.openIn();
			noOverlapFiltWriter.openOut();
			noOverlapMergedFiltWriter.openOut();

			while(reader.readNextRead(seq)){
				if(njh::in(seq.name_, postProcessResults.filteredHitsMergedNonOverlapByQuery_)){
					for(const auto & hitGroup : postProcessResults.filteredHitsMergedNonOverlapByQuery_[seq.name_]){
						Bed6RecordCore region = hitGroup.region_.genBedRecordCore();
						auto subSeq = seq.getSubRead(region.chromStart_, region.length());
						if(region.reverseStrand()){
							subSeq.reverseComplementRead(false, true);
						}
						MetaDataInName meta = hitGroup.genOutRegion().meta_;
						meta.addMeta("trimStart", region.chromStart_, true);
						meta.addMeta("trimEnd", region.chromEnd_, true);
						meta.addMeta("trimLen", region.length(), true);
						meta.addMeta("trimCov", region.length()/static_cast<double>(len(seq)), true);
						meta.addMeta("revStrand", region.reverseStrand(), true);
						meta.resetMetaInName(subSeq.name_);
						noOverlapMergedFiltWriter.write(subSeq);
					}
				}
				if(njh::in(seq.name_, postProcessResults.filteredNonOverlapHitsByQuery_)){
					for(const auto & hit : postProcessResults.filteredNonOverlapHitsByQuery_[seq.name_]){
						Bed6RecordCore region = hit.genBed6_env();
						auto subSeq = seq.getSubRead(region.chromStart_, region.length());
						if(region.reverseStrand()){
							subSeq.reverseComplementRead(false, true);
						}
						MetaDataInName meta;
						meta.addMeta("acc", hit.acc_, true);
						meta.addMeta("hmmFrom", hit.hmmFrom_, true);
						meta.addMeta("hmmTo", hit.hmmTo_, true);
						meta.addMeta("hmmCovered", hit.modelCoverage());
						meta.addMeta("trimStart", region.chromStart_, true);
						meta.addMeta("trimEnd", region.chromEnd_, true);
						meta.addMeta("trimLen", region.length(), true);
						meta.addMeta("trimCov", region.length()/static_cast<double>(len(seq)), true);
						meta.addMeta("revStrand", region.reverseStrand(), true);
						meta.addMeta("score", hit.modelScore_, true);
						meta.addMeta("scoreNorm", hit.modelScore_/region.length(), true);
						meta.addMeta("evalue", hit.modelEvalue_, true);
						meta.addMeta("model", hit.targetName_);
						meta.addMeta("ID", hit.targetDesc_);
						meta.resetMetaInName(subSeq.name_);
						noOverlapFiltWriter.write(subSeq);
					}
				} else {
					seqsWithNoDomains.emplace_back(seq.name_);
				}
			}
			OutputStream noDomainHits(seqsWithNoDomainHitsFnp);
			for(const auto & name : seqsWithNoDomains){
				noDomainHits << name << std::endl;
			}
		}
	}

	{
		//end
		//convert to fasta and non-gz hmm model

		auto inputSeqFnp = njh::files::make_path(endDir, "inputSeqs.fasta");
		auto noOverlapFiltFnp = njh::files::make_path(endDir, "noOverlapFiltHits.fasta");
		auto noOverlapMergedFiltFnp = njh::files::make_path(endDir, "noOverlapMergedFiltHits.fasta");
		auto hmmModelFnp = njh::files::make_path(endDir, "hmmModel.txt");

		std::unordered_map<std::string, GenomicRegion> regionsByName;
		std::unordered_map<std::string, uint32_t> realQueryLens;
		std::vector<GenomicRegion> regions;
		auto inputOpts = setUp.pars_.ioOptions_;
		if(!subRegions.empty()){
			auto twoBitFnp = bfs::path(setUp.pars_.ioOptions_.firstName_).replace_extension(".2bit");
			TwoBit::TwoBitFile tReader(twoBitFnp);
			realQueryLens = tReader.getSeqLens();
			auto regionsOutOpts = SeqIOOptions::genFastaOut(njh::files::make_path(endDir, "subRegions.fasta"));
			regions = bed3PtrsToGenomicRegs(getBed3s(subRegions));
			if(0 != extendSubRegions){
				njh::for_each(regions, [&realQueryLens,&extendSubRegions](GenomicRegion & region){
					BedUtility::extendLeftRight(region, extendSubRegions, extendSubRegions, njh::mapAt(realQueryLens, region.chrom_));
				});
				sortGRegionsByStart(regions);
				std::vector<GenomicRegion> mergedRegions;
				for(const auto & reg : regions){
					bool overlaps = false;
					if(!mergedRegions.empty()){
						if(mergedRegions.back().overlaps(reg, 1)) {
							overlaps = true;
							mergedRegions.back().end_ = std::max(mergedRegions.back().end_, reg.end_);
						}
					}
					if(!overlaps){
						mergedRegions.emplace_back(reg);
					}
				}
				OutputStream mergedBedOut(njh::files::make_path(endDir, "extendedSubRegions.bed"));
				njh::for_each(mergedRegions,[&mergedBedOut](GenomicRegion & region){
					region.uid_ = region.createUidFromCoords();
					mergedBedOut << region.genBedRecordCore().toDelimStr() << std::endl;
				});
				regions = mergedRegions;
			}


			SeqOutput writer(regionsOutOpts);
			writer.openOut();
			// uint32_t regionCount = 0;
			for(const auto & region : regions){
				regionsByName[region.uid_] = region;
				// ++regionCount;
				auto subRegion = region.extractSeq(tReader);
				subRegion.name_ = region.uid_;
				writer.write(subRegion);
			}
			inputOpts = SeqIOOptions::genFastaIn(regionsOutOpts.out_.outName());
		}

		auto seqKey = SeqIO::rewriteSeqsWithIndexAsName(
						inputOpts,
						SeqIOOptions::genFastaOut(inputSeqFnp),
						njh::files::make_path(endDir, "inputSeqNameKey.tab.txt"));

		//write out the model, this way you can supply it in gzip format but this will unzip it
		njh::files::reWriteFile(hmmModelEnd, hmmModelFnp);
		nhmmscanOutput::run_hmmpress_ifNeed(hmmModelFnp);

		std::stringstream cmdSs;
		cmdSs << "nhmmscan " << defaultParameters
					<< " " << "--tblout " << "raw_all_hits_table.txt"
					<< " " << "hmmModel.txt"
					<< " " << "inputSeqs.fasta"
					<< " " << " > " << "nhmmscan_raw_output.txt";
		std::string cdCmd = "cd " + endDir.string() + " && ";
		auto cmdOutput = njh::sys::run( { cdCmd, cmdSs.str() });
		if (!cmdOutput.success_) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "failed to run hmmsearch "
				 << "\n";
			ss << cmdOutput.toJson() << "\n";
			throw std::runtime_error { ss.str() };
		} else {
			OutputStream nhmmscanOutput(
							njh::files::make_path(endDir,
																		"nhmmscanCmdRunDetails.json"));
			nhmmscanOutput << cmdOutput.toJson() << std::endl;
		}



		auto nhmmscan_raw_outputFnp = njh::files::make_path(endDir, "nhmmscan_raw_output.txt");
		auto seqsWithNoDomainHitsFnp = njh::files::make_path(endDir, "seqsWithNoDomainHits.tab.txt");
		nhmmscanOutput outputParsed;
		if(regions.empty()){
			outputParsed = nhmmscanOutput::parseRawOutput(nhmmscan_raw_outputFnp, seqKey);
		}else{
			outputParsed = nhmmscanOutput::parseRawOutput(nhmmscan_raw_outputFnp, seqKey, regionsByName, realQueryLens);
		}
		auto postProcessResults = outputParsed.postProcessHits(postProcessPars);
		//convert hits table into a real table and
		outputParsed.writeInfoFiles(postProcessResults, endDir);

		if(setUp.pars_.debug_){
			//merging debugging
			for(const auto & filteredHits : postProcessResults.filteredHitsByQuery_){
				auto mergedResults = nhmmscanOutput::QueryResults::mergeOverlapingHits(filteredHits.second, mergePars);
				for(const auto & merged : mergedResults){
					if(merged.hits_.size() > 1){
						std::cout << merged.region_.genBedRecordCore().toDelimStrWithExtra() << std::endl;
						for(const auto & hit : merged.hits_){
							std::cout << "\t" << hit.genBed6_env().toDelimStrWithExtra() << std::endl;
						}
						std::cout << std::endl;
					}
				}
			}
		}

		//trim seqs to best overlapping positions sub seqs
		{
			VecStr seqsWithNoDomains;
			seqInfo seq;
			SeqInput reader(setUp.pars_.ioOptions_);
			SeqOutput noOverlapFiltWriter(SeqIOOptions::genFastaOut(noOverlapFiltFnp));
			SeqOutput noOverlapMergedFiltWriter(SeqIOOptions::genFastaOut(noOverlapMergedFiltFnp));

			reader.openIn();
			noOverlapFiltWriter.openOut();
			noOverlapMergedFiltWriter.openOut();

			while(reader.readNextRead(seq)){
				if(njh::in(seq.name_, postProcessResults.filteredHitsMergedNonOverlapByQuery_)){
					for(const auto & hitGroup : postProcessResults.filteredHitsMergedNonOverlapByQuery_[seq.name_]){
						Bed6RecordCore region = hitGroup.region_.genBedRecordCore();
						auto subSeq = seq.getSubRead(region.chromStart_, region.length());
						if(region.reverseStrand()){
							subSeq.reverseComplementRead(false, true);
						}
						MetaDataInName meta = hitGroup.genOutRegion().meta_;
						meta.addMeta("trimStart", region.chromStart_, true);
						meta.addMeta("trimEnd", region.chromEnd_, true);
						meta.addMeta("trimLen", region.length(), true);
						meta.addMeta("trimCov", region.length()/static_cast<double>(len(seq)), true);
						meta.addMeta("revStrand", region.reverseStrand(), true);
						meta.resetMetaInName(subSeq.name_);
						noOverlapMergedFiltWriter.write(subSeq);
					}
				}
				if(njh::in(seq.name_, postProcessResults.filteredNonOverlapHitsByQuery_)){
					for(const auto & hit : postProcessResults.filteredNonOverlapHitsByQuery_[seq.name_]){
						Bed6RecordCore region = hit.genBed6_env();
						auto subSeq = seq.getSubRead(region.chromStart_, region.length());
						if(region.reverseStrand()){
							subSeq.reverseComplementRead(false, true);
						}
						MetaDataInName meta;
						meta.addMeta("acc", hit.acc_, true);
						meta.addMeta("hmmFrom", hit.hmmFrom_, true);
						meta.addMeta("hmmTo", hit.hmmTo_, true);
						meta.addMeta("hmmCovered", hit.modelCoverage());
						meta.addMeta("trimStart", region.chromStart_, true);
						meta.addMeta("trimEnd", region.chromEnd_, true);
						meta.addMeta("trimLen", region.length(), true);
						meta.addMeta("trimCov", region.length()/static_cast<double>(len(seq)), true);
						meta.addMeta("revStrand", region.reverseStrand(), true);
						meta.addMeta("score", hit.modelScore_, true);
						meta.addMeta("scoreNorm", hit.modelScore_/region.length(), true);
						meta.addMeta("evalue", hit.modelEvalue_, true);
						meta.addMeta("model", hit.targetName_);
						meta.addMeta("ID", hit.targetDesc_);
						meta.resetMetaInName(subSeq.name_);
						noOverlapFiltWriter.write(subSeq);
					}
				} else {
					seqsWithNoDomains.emplace_back(seq.name_);
				}
			}
			OutputStream noDomainHits(seqsWithNoDomainHitsFnp);
			for(const auto & name : seqsWithNoDomains){
				noDomainHits << name << std::endl;
			}
		}
	}

	auto startBedRegionsFnp = njh::files::make_path(startDir, "nhmmscan_hits_nonOverlap_filtered.bed");
	auto endBedRegionsFnp = njh::files::make_path(endDir, "nhmmscan_hits_nonOverlap_filtered.bed");

	std::unordered_map<std::string, std::vector<GenomicRegion>> startsBySeq;
	auto startRegions  = bedPtrsToGenomicRegs(getBeds(startBedRegionsFnp));
	for(const auto & start : startRegions){
		startsBySeq[start.chrom_].emplace_back(start);
	}
	std::unordered_map<std::string, std::vector<GenomicRegion>> endsBySeq;
	auto endRegions  = bedPtrsToGenomicRegs(getBeds(endBedRegionsFnp));
	for(const auto & end : endRegions){
		endsBySeq[end.chrom_].emplace_back(end);
	}
	SeqIOOptions trimmedOutOpts(njh::files::make_path(setUp.pars_.directoryName_, "trimmed"), setUp.pars_.ioOptions_.outFormat_);
	SeqOutput seqWriter(trimmedOutOpts);
	seqWriter.openOut();
	for(const auto & seq : allSeqs){
		if(njh::in(seq.second.name_, startsBySeq) && njh::in(seq.second.name_, endsBySeq)){
			for(const auto & start : startsBySeq[seq.second.name_]){
				for(const auto & end : endsBySeq[seq.second.name_]){
					uint32_t startPos = start.start_;
					uint32_t endPos = end.end_;
					if(excludeHmmHitRegion){
						startPos = start.end_;
						endPos = end.start_;
					}
					if(startPos < endPos){
						auto subSeq = seq.second.getSubRead(startPos, endPos - startPos);
						seqWriter.write(subSeq);
					} else {
						//no trim
					}
				}
			}
		} else {
			//no trim
		}
	}

	return 0;
}

} // namespace njhseq

