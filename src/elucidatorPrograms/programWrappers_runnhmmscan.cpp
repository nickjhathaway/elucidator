/*
 * programWrappers_runnhmmscan.cpp
 *
 *  Created on: Jan 15, 2022
 *      Author: nick
 */






#include "programWrappers.hpp"
#include <njhseq/IO/SeqIO/SeqIO.hpp>
#include <njhseq/objects/BioDataObject/BioDataFileIO.hpp>
#include <njhseq/objects/BioDataObject/BedRecordCore.hpp>
#include <njhseq/objects/BioDataObject/BioRecordsUtils/HmmerUtility.hpp>

namespace njhseq {



int programWrapperRunner::runnhmmscan(const njh::progutils::CmdArgs & inputCommands){
	std::string defaultParameters = "--nonull2 --incT 20 --incdomT 20 -T 20 --notextw";
	bfs::path hmmModel;
	nhmmscanOutput::PostProcessHitsPars postProcessPars;

	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processDefaultReader(VecStr{"fasta", "fastagz", "fastq", "fastqgz"}, true);
	setUp.setOption(defaultParameters, "--defaultParameters", "The default parameters given to hmmsearch");
	setUp.setOption(hmmModel, "--hmmModel", "hmm model database, created by hmmbuild", true);
	setUp.setOption(postProcessPars.hmmStartFilter, "--hmmStartFilter", "Filter partial hmms domain hits if they start or end this far into the model");
	setUp.setOption(postProcessPars.minLength, "--minLength", "Minimum output domain hit length");
	setUp.setOption(postProcessPars.accCutOff, "--accCutOff", "accuracy cut off");
	setUp.setOption(postProcessPars.scoreCutOff, "--scoreCutOff", "score cut off");

	setUp.processDirectoryOutputName(true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);

	//convert to fasta and non-gz hmm model

	auto inputSeqFnp = njh::files::make_path(setUp.pars_.directoryName_, "inputSeqs.fasta");
	auto inputSeqSubDomainsFnp = njh::files::make_path(setUp.pars_.directoryName_, "nonOverlapFiltHits.fasta");
	auto hmmModelFnp = njh::files::make_path(setUp.pars_.directoryName_, "hmmModel.txt");
	auto seqKey = SeqIO::rewriteSeqsWithIndexAsName(
					setUp.pars_.ioOptions_,
					SeqIOOptions::genFastaOut(inputSeqFnp),
					njh::files::make_path(setUp.pars_.directoryName_, "inputSeqNameKey.tab.txt"));

	//write out the model, this way you can supply it in gzip format but this will unzip it
	njh::files::reWriteFile(hmmModel, hmmModelFnp);
	nhmmscanOutput::run_hmmpress_ifNeed(hmmModelFnp);

	std::stringstream cmdSs;
	cmdSs << "nhmmscan " << defaultParameters
			<< " " << "--tblout " << "raw_all_domain_hits_table.txt"
			<< " " << "hmmModel.txt"
			<< " " << "inputSeqs.fasta"
			<< " " << " > " << "nhmmscan_raw_output.txt";
	std::string cdCmd = "cd " + setUp.pars_.directoryName_ + " && ";
	auto cmdOutput = njh::sys::run( { cdCmd, cmdSs.str() });
	if (!cmdOutput.success_) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "failed to run hmmsearch "
				<< "\n";
		ss << cmdOutput.toJson() << "\n";
		throw std::runtime_error { ss.str() };
	} else {
		OutputStream nhmmscanOutput(
				njh::files::make_path(setUp.pars_.directoryName_,
						"nhmmscanCmdRunDetails.json"));
		nhmmscanOutput << cmdOutput.toJson() << std::endl;
	}

	//convert hits table into a real table and
	//filter table
//	auto rawDomainHitsFnp = njh::files::make_path(setUp.pars_.directoryName_, "raw_all_domain_hits_table.txt");
//
//	auto allDomainHitsFnp = njh::files::make_path(setUp.pars_.directoryName_, "all_domain_hits_table.tab.txt");
//	auto filtDomainHitsFnp = njh::files::make_path(setUp.pars_.directoryName_, "filt_domain_hits_table.tab.txt");
//	auto nonOverlappingDomainHitsFnp = njh::files::make_path(setUp.pars_.directoryName_, "nonOverlappingBestdomain_hits_table.tab.txt");
//	auto bedRegionsFnp = njh::files::make_path(setUp.pars_.directoryName_, "nonOverlappingBestdomain_hits.bed");
//
//	auto domainCountsFnp = njh::files::make_path(setUp.pars_.directoryName_, "seqDomainCounts.tab.txt");

	auto nhmmscan_raw_outputFnp = njh::files::make_path(setUp.pars_.directoryName_, "nhmmscan_raw_output.txt");
	auto seqsWithNoDomainHitsFnp = njh::files::make_path(setUp.pars_.directoryName_, "seqsWithNoDomainHits.tab.txt");



	nhmmscanOutput outputParsed = nhmmscanOutput::parseRawOutput(nhmmscan_raw_outputFnp, seqKey);
	auto postProcessResults = outputParsed.postProcessHits(postProcessPars);
	outputParsed.writeInfoFiles(postProcessResults, setUp.pars_.directoryName_);

	//trim seqs to best overlapping positions sub seqs
	{
		VecStr seqsWithNoDomains;
		seqInfo seq;
		SeqInput reader(setUp.pars_.ioOptions_);
		SeqOutput writer(SeqIOOptions::genFastaOut(inputSeqSubDomainsFnp));
		reader.openIn();
		writer.openOut();

		while(reader.readNextRead(seq)){
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
					writer.write(subSeq);
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

//	{
//		OutputStream out(allDomainHitsFnp);
//		OutputStream outFilt(filtDomainHitsFnp);
//
//		BioDataFileIO<HmmerSeqHitTab> reader{IoOptions(InOptions(rawDomainHitsFnp))};
//		HmmerSeqHitTab domain;
//		reader.openIn();
//		uint32_t count = 0;
//		out << "#" << njh::conToStr(HmmerSeqHitTab::toDelimStrHeader(), "\t") << std::endl;
//		outFilt << "#" << njh::conToStr(HmmerSeqHitTab::toDelimStrHeader(), "\t") << std::endl;
//
//		while(reader.readNextRecord(domain)){
//			out << domain.toDelimStr() << std::endl;
//			if(!(domain.envFrom_ > hmmStartFilter && domain.hmmFrom_ > hmmStartFilter )){
//				outFilt << domain.toDelimStr() << std::endl;
//			}
//			++count;
//		}
//	}
//	//get best non-overlapping positions
//	std::vector<Bed6RecordCore> filteredRegions;
//
//	std::unordered_map<std::string, std::vector<HmmerSeqHitTab>> domainsPerSeq;
//
//	{
//		std::vector<HmmerSeqHitTab> domains;
//		std::vector<Bed6RecordCore> locations;
////		std::cout << __FILE__ << " " << __LINE__ << std::endl;
//		{
//			BioDataFileIO<HmmerSeqHitTab> reader{IoOptions(InOptions(filtDomainHitsFnp))};
//
//			HmmerSeqHitTab domain;
//			reader.openIn();
//			uint32_t count = 0;
//			while(reader.readNextRecord(domain)){
//				if(domain.genBed6_env().length() < minLength){
//					continue;
//				}
//
////				std::cout << __FILE__ << " " << __LINE__ << std::endl;
//				domain.queryName_ = seqKey[njh::StrToNumConverter::stoToNum<uint32_t>(domain.queryName_)];
//				domains.emplace_back(domain);
////				std::cout << __FILE__ << " " << __LINE__ << std::endl;
//				auto genLoc = domain.genBed6_env();
//				genLoc.name_ = njh::pasteAsStr(count);
//				locations.emplace_back(genLoc);
////				std::cout << __FILE__ << " " << __LINE__ << std::endl;
//
//				++count;
//			}
//		}
//
////		std::cout << __FILE__ << " " << __LINE__ << std::endl;
//		//sort by scores, with lowest e-score on top
//		njh::sort(locations,[&domains](const Bed6RecordCore & rec1, const Bed6RecordCore & rec2){
//			if(domains[njh::StrToNumConverter::stoToNum<uint32_t>(rec1.name_)].modelEvalue_ == domains[njh::StrToNumConverter::stoToNum<uint32_t>(rec2.name_)].modelEvalue_){
//				if(domains[njh::StrToNumConverter::stoToNum<uint32_t>(rec1.name_)].modelScore_ == domains[njh::StrToNumConverter::stoToNum<uint32_t>(rec2.name_)].modelScore_){
//					return rec1.length() > rec2.length();
//				}else{
//					return domains[njh::StrToNumConverter::stoToNum<uint32_t>(rec1.name_)].modelScore_ > domains[njh::StrToNumConverter::stoToNum<uint32_t>(rec2.name_)].modelScore_;
//				}
//			}else{
//				return domains[njh::StrToNumConverter::stoToNum<uint32_t>(rec1.name_)].modelEvalue_ < domains[njh::StrToNumConverter::stoToNum<uint32_t>(rec2.name_)].modelEvalue_;
//			}
//		});

////		std::cout << __FILE__ << " " << __LINE__ << std::endl;
//		for(const auto & region : locations){
//			bool overlap = false;
//			for(const auto & outRegion : filteredRegions){
//				if(region.overlaps(outRegion,1)){
//					overlap = true;
//					break;
//				}
//			}
//			if(!overlap){
//				filteredRegions.emplace_back(region);
//			}
//		}
////		std::cout << __FILE__ << " " << __LINE__ << std::endl;
//		OutputStream outFiltNonOverlap(nonOverlappingDomainHitsFnp);
//		outFiltNonOverlap << "#" << njh::conToStr(HmmerSeqHitTab::toDelimStrHeader(), "\t") << "\tbedName" << std::endl;
//		BedUtility::coordSort(filteredRegions);
//		OutputStream bedRegionsOut(bedRegionsFnp);
////		std::cout << __FILE__ << " " << __LINE__ << std::endl;
//		std::unordered_map<std::string, uint32_t> modelCounts;
//		for(auto & region : filteredRegions){
//			const auto & domForReg  = domains[njh::StrToNumConverter::stoToNum<uint32_t>(region.name_)];
//
//			auto bedNewName = njh::pasteAsStr(domForReg.targetName_, ".", modelCounts[domForReg.targetName_]);
//			++modelCounts[domForReg.targetName_];
//			outFiltNonOverlap << domForReg.toDelimStr() << "\t" << bedNewName<< std::endl;
//			region.name_ = bedNewName;
//			bedRegionsOut << region.toDelimStrWithExtra() << std::endl;
//			domainsPerSeq[domForReg.queryName_].emplace_back(domForReg);
//			domainsPerSeq[domForReg.queryName_].back().targetDesc_ = bedNewName;
//		}
////		std::cout << __FILE__ << " " << __LINE__ << std::endl;
//
//
//		BedUtility::coordSort(locations, false);
//		OutputStream domainLocationBedsOut(domainLocationBedsFnp);
//		for( auto & reg : locations){
//			reg.score_ = reg.length();
//			const auto & domForReg  = domains[njh::StrToNumConverter::stoToNum<uint32_t>(reg.name_)];
//			reg.name_ = domForReg.queryName_;
//			domainLocationBedsOut << reg.toDelimStrWithExtra() << std::endl;
// 		}
//	}
//

//

	return 0;
}


} // namespace njhseq

