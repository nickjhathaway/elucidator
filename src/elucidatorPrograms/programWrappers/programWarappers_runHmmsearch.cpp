/*
 * programWarappers_runHmmsearch.cpp
 *
 *  Created on: May 27, 2021
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



#include "programWrappers.hpp"
#include <njhseq/IO/SeqIO/SeqIO.hpp>
#include <njhseq/objects/BioDataObject/BioDataFileIO.hpp>
#include <njhseq/objects/BioDataObject/BedRecordCore.hpp>
#include <njhseq/objects/BioDataObject/BioRecordsUtils/BedUtility.hpp>
#include <njhseq/objects/BioDataObject/BioRecordsUtils/HmmerUtility.hpp>


namespace njhseq {

namespace bfs = boost::filesystem;


int programWrapperRunner::runHmmsearch(const njh::progutils::CmdArgs & inputCommands){
	std::string defaultParameters = "--nonull2 --incT 10 --incdomT 10 -T 10 --domT 10";
	bfs::path hmmModel;
	uint32_t hmmStartFilter = 25;
	uint32_t targetStartFilter = std::numeric_limits<uint32_t>::max();

	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processDefaultReader(VecStr{"fasta", "fastagz", "fastq", "fastqgz"}, true);
	setUp.setOption(defaultParameters, "--defaultParameters", "The default parameters given to hmmsearch");
	setUp.setOption(hmmModel, "--hmmModel", "hmm model database, created by hmmbuild", true);
	setUp.setOption(hmmStartFilter, "--hmmStartFilter", "Filter partial hmms domain hits if they start or end this far into the model");
	setUp.setOption(targetStartFilter, "--targetStartFilter", "Filter partial hmms domain hits if they start or end this far into the target");

	setUp.processDirectoryOutputName(true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);

	//convert to fasta and non-gz hmm model
	std::unordered_map<uint32_t, std::string> seqKey;
	auto inputSeqFnp = njh::files::make_path(setUp.pars_.directoryName_, "inputSeqs.fasta");
	auto inputSeqSubDomainsFnp = njh::files::make_path(setUp.pars_.directoryName_, "inputSeqsSubDomains.fasta");

	auto hmmModelFnp = njh::files::make_path(setUp.pars_.directoryName_, "hmmModel.txt");

	{
		seqInfo seq;
		SeqInput reader(setUp.pars_.ioOptions_);
		SeqOutput writer(SeqIOOptions::genFastaOut(inputSeqFnp));
		OutputStream seqNameKeyOut(njh::files::make_path(setUp.pars_.directoryName_, "inputSeqNameKey.tab.txt"));
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
		if(njh::endsWith(hmmModel.string(), ".gz")) {
			//write out the model, this way you can supply it in gzip format but this will unzip it
			InputStream inModel(hmmModel);
			OutputStream outModel(hmmModelFnp);
			outModel << inModel.rdbuf();
			outModel.flush();
		} else {
			bfs::create_symlink(hmmModel, hmmModelFnp);
		}
	}

	std::stringstream cmdSs;
	cmdSs << "hmmsearch " << defaultParameters
			<< " " << "--domtblout " << "raw_all_domain_hits_table.txt"
			<< " " << "hmmModel.txt"
			<< " " << "inputSeqs.fasta"
			<< " " << " > " << "hmmsearch_raw_output.txt";
	std::string cdCmd = "cd " + setUp.pars_.directoryName_ + " && ";
	auto cmdOutput = njh::sys::run({cdCmd, cmdSs.str()});
	if(!cmdOutput.success_){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "failed to run hmmsearch " << "\n";
		ss << cmdOutput.toJson() << "\n";
		throw std::runtime_error{ss.str()};
	}else{
		OutputStream hmmsearchOutput(njh::files::make_path(setUp.pars_.directoryName_, "hmmsearchCmdRunDetails.json"));
		hmmsearchOutput << cmdOutput.toJson() << std::endl;
	}

	//convert domain hits table into a real table and
	//filter table
	auto rawDomainHitsFnp = njh::files::make_path(setUp.pars_.directoryName_, "raw_all_domain_hits_table.txt");
	auto allDomainHitsFnp = njh::files::make_path(setUp.pars_.directoryName_, "all_domain_hits_table.tab.txt");
	auto filtDomainHitsFnp = njh::files::make_path(setUp.pars_.directoryName_, "filt_domain_hits_table.tab.txt");
	auto nonOverlappingDomainHitsFnp = njh::files::make_path(setUp.pars_.directoryName_, "nonOverlappingBestdomain_hits_table.tab.txt");
	auto bedRegionsFnp = njh::files::make_path(setUp.pars_.directoryName_, "nonOverlappingBestdomain_hits.bed");
	auto seqsWithNoDomainHitsFnp = njh::files::make_path(setUp.pars_.directoryName_, "seqsWithNoDomainHits.tab.txt");
	auto domainCountsFnp = njh::files::make_path(setUp.pars_.directoryName_, "seqDomainCounts.tab.txt");

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
			bool passTargetStart = domain.env0BasedPlusStrandStart() <= targetStartFilter;
			bool passTargetEnd = domain.targetLen_ - domain.envTo_ <= targetStartFilter;
			if(passHmmStart && passHmmEnd && passTargetStart && passTargetEnd) {
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


	return 0;
}

} // namespace njhseq



