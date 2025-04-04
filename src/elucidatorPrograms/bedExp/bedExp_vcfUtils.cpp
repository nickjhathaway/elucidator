//
// Created by Nicholas Hathaway on 1/9/25.
//


#include "bedExp.hpp"
#include <njhseq/objects/BioDataObject.h>

#include "elucidator/objects/BioDataObject.h"
#include "elucidator/BioRecordsUtils/BedUtility.hpp"
#include <njhseq/objects/Gene/TranslatorByAlignment.hpp>



namespace njhseq {

int bedExpRunner::vcfToBed(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path vcfFile;
	bfs::path intersectWithBed;
	OutOptions outOpts;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.setOption(vcfFile, "--vcfFile", "vcfFile", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	OutputStream out(outOpts);

	VCFOutput vcf = VCFOutput::readInHeader(vcfFile);
	InputStream in(vcfFile);
	std::string line;
	// uint32_t count = 0;
	while(njh::files::crossPlatGetline(in, line)) {
		if(line.front() != '#') {
			// std::cout << count++ << std::endl;
			out << vcf.processRecordLineForFixedData(line).genRegion().genBedRecordCore().toDelimStrWithExtra() << std::endl;
		}
	}

	return 0;
}

int bedExpRunner::printVcfSamples(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path vcfFile;
	bfs::path intersectWithBed;
	OutOptions outOpts;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.setOption(vcfFile, "--vcfFile", "vcfFile", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	OutputStream out(outOpts);

	VCFOutput vcf = VCFOutput::readInHeader(vcfFile);
	out << njh::conToStr(vcf.samples_, "\n") << std::endl;
	return 0;
}


int bedExpRunner::vcfRenameChroms(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path vcfFile;
	bfs::path nameKeyFnp;
	std::string old_name_column_name;
	std::string new_name_column_name;

	OutOptions outOpts;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.setOption(vcfFile, "--vcfFile", "vcfFile", true);
	setUp.setOption(nameKeyFnp, "--nameKeyFnp", "name Key Fnp, tab-delimited file either no column file with col1 being old name and col2 being new name or can supply which column names are old and new names with flags --oldNameColumnName and --newNameColumnName", true);
	setUp.setOption(old_name_column_name, "--oldNameColumnName", "the name of a column to be the old name");
	setUp.setOption(new_name_column_name, "--newNameColumnName", "the name of the column to be the new/replacement name", "" != old_name_column_name);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	std::unordered_map<std::string, std::string> nameKeyMap;
	std::unordered_map<std::string, std::string> replacementNameToOriginalName;
	if (new_name_column_name.empty()) {
		table keyTab(nameKeyFnp, "\t", false);
		if (keyTab.nCol() != 2) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << ", error " << nameKeyFnp <<
					" should have two columns, not: " << keyTab.nCol() << "\n";
			throw std::runtime_error{ss.str()};
		}
		auto old_name_col_pos = 0;
		auto new_name_col_pos = 1;
		for (const auto & row : keyTab) {
			if (njh::in(row[old_name_col_pos], nameKeyMap)) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "already have " << row[old_name_col_pos] << " in replacement map" << "\n";
				throw std::runtime_error{ss.str()};
			}
			if (njh::in(row[new_name_col_pos], replacementNameToOriginalName)) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "already have replacement name: " << row[new_name_col_pos] << " for " << replacementNameToOriginalName[row[new_name_col_pos]] << "\n";
				throw std::runtime_error{ss.str()};
			}
			nameKeyMap[row[old_name_col_pos]] = row[new_name_col_pos];
			replacementNameToOriginalName[row[new_name_col_pos]] = row[old_name_col_pos];
		}
	} else {
		table keyTab(nameKeyFnp, "\t", true);
		keyTab.checkForColumnsThrow(VecStr{old_name_column_name, new_name_column_name}, __PRETTY_FUNCTION__);
		auto old_name_col_pos = keyTab.getColPos(old_name_column_name);
		auto new_name_col_pos = keyTab.getColPos(new_name_column_name);
		for (const auto & row : keyTab) {
			if (njh::in(row[old_name_col_pos], nameKeyMap)) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "already have " << row[old_name_col_pos] << " in replacement map" << "\n";
				throw std::runtime_error{ss.str()};
			}
			if (njh::in(row[new_name_col_pos], replacementNameToOriginalName)) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "already have replacement name: " << row[new_name_col_pos] << " for " << replacementNameToOriginalName[row[new_name_col_pos]] << "\n";
				throw std::runtime_error{ss.str()};
			}
			nameKeyMap[row[old_name_col_pos]] = row[new_name_col_pos];
			replacementNameToOriginalName[row[new_name_col_pos]] = row[old_name_col_pos];
		}
	}

	auto original_vcf_header = VCFOutput::readInHeader(vcfFile);
	auto output_vcf_header = original_vcf_header;
	output_vcf_header.changeContigNames(nameKeyMap);
	OutputStream out(outOpts);
	output_vcf_header.writeOutFixedAndSampleMeta(out);

	InputStream in(vcfFile);
	std::string line;
	std::string formatOut;
	VecStr formatOutputOrder;
	//force GT to be first field
	if(njh::in(std::string("GT"), output_vcf_header.formatEntries_)) {
		formatOutputOrder.emplace_back("GT");
	}
	for (const auto & infoKey: output_vcf_header.formatEntries_) {
		if(infoKey.first != "GT") {
			formatOutputOrder.emplace_back(infoKey.first);
		}
	}
	for (const auto & infoKey: formatOutputOrder) {
		const auto & format = output_vcf_header.formatEntries_.at(infoKey);
		if(!formatOut.empty()) {
			formatOut +=":";
		}
		formatOut += format.id_;
	}
	while (njh::files::crossPlatGetline(in, line)) {
		if(line.front() != '#') {
			auto record = original_vcf_header.processRecordLineForFixedDataAndSampleMetaData(line);
			if (njh::notIn(record.chrom_, nameKeyMap)) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "don't have chrom " << record.chrom_ << " in name key, options are: " << njh::conToStr(njh::getVecOfMapKeys(nameKeyMap), ",") << "\n";
				throw std::runtime_error{ss.str()};
			}
			record.chrom_ = nameKeyMap[record.chrom_];
			out << record.chrom_
					<< "\t" << record.pos_
					<< "\t" << record.id_
					<< "\t" << record.ref_
					<< "\t" << njh::conToStr(record.alts_, ",")
					<< "\t" << (record.qual_ == std::numeric_limits<uint32_t>::max() ? "." : estd::to_string(record.qual_))
					<< "\t" << record.filter_;
			std::string infoOut;
			for (const auto& infoKey: output_vcf_header.infoEntries_) {
				const auto& info = infoKey.second;
				if (infoKey.second.type_ == "Flag") {
					if (record.info_.containsMeta(info.id_)) {
						if (!infoOut.empty()) {
							infoOut += ";";
						}
						infoOut += info.id_;
					}
				} else {
					if (!infoOut.empty()) {
						infoOut += ";";
					}
					infoOut += info.id_ + "=" + record.info_.getMeta(info.id_);
				}
			}
			out << "\t" << infoOut;
			out << "\t" << formatOut;
			for (const auto& sampleName: output_vcf_header.samples_) {
				const auto& sample = record.sampleFormatInfos_.at(sampleName);
				std::string formatOutForSample;
				for (const auto& infoKey: formatOutputOrder) {
					const auto& format = output_vcf_header.formatEntries_.at(infoKey);
					if (!formatOutForSample.empty()) {
						formatOutForSample += ":";
					}
					formatOutForSample += sample.getMeta(format.id_);
				}
				out << "\t" << formatOutForSample;
			}
			out << std::endl;
		}
	}
	return 0;
}



int bedExpRunner::combineVcfs(const njh::progutils::CmdArgs & inputCommands) {
	std::vector<bfs::path> vcfFnps;
	bfs::path intersectWithBed;
	OutOptions outOpts;
	VCFOutput::comnbineVCFsPars combiningVcfPars;

	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.setOption(vcfFnps, "--vcfFnps", "vcfs", true);
	setUp.setOption(combiningVcfPars.ploidy, "--ploidy", "Ploidy to force for the sample for the vcf files");
	setUp.setOption(combiningVcfPars.doNotRescueVariantCallsAcrossTargets, "--doNotRescueVariantCallsAcrossTargets", "do Not Rescue Variant Calls Across Targets");
	setUp.setOption(combiningVcfPars.combinedOverlappingCallsAcrossTargets, "--combineOverlappingCallsAcrossTargets", "Rather than taking the best variant call for overlapping targets, sum them instead");
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	OutputStream out(outOpts);

	auto combined = VCFOutput::comnbineVCFs(vcfFnps, combiningVcfPars);
	combined.writeOutFixedAndSampleMeta(out);
	return 0;
}

int bedExpRunner::simpleVCFDetermineMonoclonals(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path vcfFnp;
	bfs::path intersectWithBed;
	OutOptions outOpts;
	double minoraf = 0.05;
	uint32_t mindepth = 5;
	uint32_t numThreads = 1;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.setOption(mindepth, "--mindepth", "minimum depth to call a region");
	setUp.setOption(minoraf, "--minoraf", "minor allele frequency cut off");
	setUp.setOption(numThreads, "--numThreads", "number of Threads");

	setUp.setOption(vcfFnp,  "--vcfFnp", "vcf", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	OutputStream out(outOpts);
	auto firstVcf = VCFOutput::readInHeader(vcfFnp);
	if (njh::notIn(std::string("DP"), firstVcf.formatEntries_) || njh::notIn(std::string("AD"), firstVcf.formatEntries_)) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "need to have DP and AD in format entries, only have: " << njh::conToStr(njh::getVecOfMapKeys(firstVcf.formatEntries_), ",") << "\n";
		throw std::runtime_error{ss.str()};
	}



	struct LockedInputStream {
		LockedInputStream(const bfs::path & fnp):in_(fnp) {

		}
		InputStream in_;

		bool getNextLinesLocked(VecStr & lines, const uint32_t batch_count = 5) {
			lines.clear();
			std::lock_guard<std::mutex> lock(in_.mut_);
			std::string line;
			while (njh::files::crossPlatGetline(in_, line) && lines.size() < batch_count) {
				lines.emplace_back(line);
			}
			return !lines.empty();
		}
	};
	LockedInputStream in(vcfFnp);

	std::unordered_map<std::string, std::vector<bool>> samplePositionMixCounts;
	uint32_t variant_count = 0;
	uint32_t biallelic_SNP_count = 0;
	std::mutex samplePositionMixCountsMutex;
	std::function<void()> processLines = [&in, &samplePositionMixCountsMutex, &samplePositionMixCounts,&variant_count,&biallelic_SNP_count,&firstVcf, &mindepth, &minoraf,&setUp]() {
		VecStr lines;
		std::unordered_map<std::string, std::vector<bool>> current_samplePositionMixCounts;
		uint32_t current_variant_count = 0;
		uint32_t current_biallelic_SNP_count = 0;
		while (in.getNextLinesLocked(lines)) {
			for (const auto & line : lines) {
				if(line.front() != '#') {
					++current_variant_count;
					if (setUp.pars_.verbose_) {
						std::cout << current_variant_count << std::endl;
					}
					auto record = firstVcf.processRecordLineForFixedDataAndSampleMetaData(line);
					//check if biallelic SNPs
					if (record.alts_.size() == 1 && record.ref_.size() == 1 && record.alts_.front().size() == 1) {
						++current_biallelic_SNP_count;
						for (const auto & sample : record.sampleFormatInfos_) {
							if ("."  != sample.second.getMeta("DP")) {
								auto DP = sample.second.getMeta<uint32_t>("DP");
								if (DP >= mindepth) {
									uint32_t alleleCounts = 0;
									for (const auto & allele_AD : njh::tokenizeString(sample.second.getMeta("AD"), ",")) {
										auto AD = njh::StrToNumConverter::stoToNum<uint32_t>(allele_AD);
										if (AD/static_cast<double>(DP) > minoraf) {
											++alleleCounts;
										}
									}
									// if (0 == alleleCounts) {
									// 	std::cout << sample.first << std::endl;
									// 	std::cout << "sample.second.getMeta(\"DP\"): " << sample.second.getMeta("DP") << std::endl;
									// 	std::cout << "sample.second.getMeta(\"AD\"): " << sample.second.getMeta("AD") << std::endl;
									// 	exit(1);
									// }
									if (1 == alleleCounts) {
										current_samplePositionMixCounts[sample.first].emplace_back(true);
									} else {
										current_samplePositionMixCounts[sample.first].emplace_back(false);
									}
								}
							}
						}
					}
				}
			}
		}
		{
			//add in counts;
			std::lock_guard<std::mutex> lock(samplePositionMixCountsMutex);
			for (const auto & sample : current_samplePositionMixCounts) {
				addOtherVec(samplePositionMixCounts[sample.first], sample.second);
			}
			biallelic_SNP_count += current_biallelic_SNP_count;
			variant_count += current_variant_count;
		}
	};
	njh::concurrent::runVoidFunctionThreaded(processLines, numThreads);
	// while(njh::files::crossPlatGetline(in, line)) {
	// 	if(line.front() != '#') {
	// 		++variant_count;
	// 		if (setUp.pars_.verbose_) {
	// 			std::cout << variant_count << std::endl;
	// 		}
	// 		auto record = firstVcf.processRecordLineForFixedDataAndSampleMetaData(line);
	// 		//check if biallelic SNPs
	// 		if (record.alts_.size() == 1 && record.ref_.size() == 1 && record.alts_.front().size() == 1) {
	// 			++biallelic_SNP_count;
	// 			for (const auto & sample : record.sampleFormatInfos_) {
	// 				if ("."  != sample.second.getMeta("DP")) {
	// 					auto DP = sample.second.getMeta<uint32_t>("DP");
	// 					if (DP >= mindepth) {
	// 						uint32_t alleleCounts = 0;
	// 						for (const auto & allele_AD : njh::tokenizeString(sample.second.getMeta("AD"), ",")) {
	// 							auto AD = njh::StrToNumConverter::stoToNum<uint32_t>(allele_AD);
	// 							if (AD/static_cast<double>(DP) > minoraf) {
	// 								++alleleCounts;
	// 							}
	// 						}
	// 						// if (0 == alleleCounts) {
	// 						// 	std::cout << sample.first << std::endl;
	// 						// 	std::cout << "sample.second.getMeta(\"DP\"): " << sample.second.getMeta("DP") << std::endl;
	// 						// 	std::cout << "sample.second.getMeta(\"AD\"): " << sample.second.getMeta("AD") << std::endl;
	// 						// 	exit(1);
	// 						// }
	// 						if (1 == alleleCounts) {
	// 							samplePositionMixCounts[sample.first].emplace_back(true);
	// 						} else {
	// 							samplePositionMixCounts[sample.first].emplace_back(false);
	// 						}
	// 					}
	// 				}
	// 			}
	// 		}
	// 	}
	// }
	out << "sample\tfreq_of_monocalls\tcount_of_monoclonals\ttotal_call_for_sample\tfraction_of_sample_covered\ttotal_file_biallelic_SNP_count" << std::endl;
	auto sampleNames = njh::getVecOfMapKeys(samplePositionMixCounts);
	njh::naturalSortNameSet(sampleNames);
	for (const auto & sample_name : sampleNames) {
		const auto & sample = samplePositionMixCounts[sample_name];
		auto sum = vectorSum(sample);
		out << sample_name << "\t" << sum/sample.size() << "\t" << sum << "\t" << sample.size() << "\t" <<  sample.size() /static_cast<double>(biallelic_SNP_count) << "\t" <<biallelic_SNP_count << std::endl;
	}
	return 0;
}

}  //namespace njhseq


