/*
 * programWrappers_parsePrimer3Output.cpp
 *
 *  Created on: Aug 9, 2017
 *      Author: nick
 */



#include "programWrappers.hpp"
#include "elucidator/objects/BioDataObject.h"
#include "elucidator/BamToolsUtils.h"



namespace njhseq {


int programWrapperRunner::parsePrimer3OutputToPossibleMipArms(const njh::progutils::CmdArgs & inputCommands){
	bfs::path extInput = "";
	bfs::path ligInput = "";
	uint32_t minLen = 100;
	uint32_t maxLen = 400;
	OutOptions outOpts(bfs::path(""));
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processWritingOptions(outOpts);
	setUp.setOption(minLen, "--minLen", "Minimum length");
	setUp.setOption(maxLen, "--maxLen", "Maximum length");
	setUp.setOption(ligInput, "--ligInput", "input file for ligation arm", true);
	setUp.setOption(extInput, "--extInput", "input file for extension arm", true);
	setUp.finishSetUp(std::cout);


	OutputStream out(outOpts);

	auto ligResults = Primer3Results::parsePrimer3OutputResults(ligInput);
	auto ligRegions = ligResults.front()->genPrimersRegions();

	auto extResults = Primer3Results::parsePrimer3OutputResults(extInput);
	auto extRegions = extResults.front()->genPrimersRegions();

	struct PrimerPair {
		PrimerPair(const std::shared_ptr<Primer3Results::Primer> & extPrimer,
				       const std::shared_ptr<Primer3Results::Primer> & ligPrimer) :
				extPrimer_(extPrimer), ligPrimer_(ligPrimer) {
		}

		std::shared_ptr<Primer3Results::Primer> extPrimer_;
		std::shared_ptr<Primer3Results::Primer> ligPrimer_;

		GenomicRegion genRegion(const std::string & seqId) const{
			MetaDataInName meta;
			meta.addMeta("ext", extPrimer_->name_);
			meta.addMeta("lig", ligPrimer_->name_);
			meta.addMeta("seqId", seqId);
			auto extRegion = extPrimer_->genRegion(seqId);
			auto ligRegion = ligPrimer_->genRegion(seqId);
			return GenomicRegion(meta.createMetaName(), seqId,
					std::min(extRegion.start_, ligRegion.start_),
					std::max(extRegion.end_, ligRegion.end_), extRegion.reverseSrand_);
		}
	};
	std::vector<PrimerPair> pairs;
	for(const auto & ext : extRegions){
		for(const auto & lig : ligRegions){
			if(!ext.second.reverseSrand_){
				if(lig.second.reverseSrand_){
					if(lig.second.start_ > ext.second.end_ &&
							lig.second.end_ - ext.second.start_ > minLen &&
							lig.second.end_ - ext.second.start_ < maxLen){
						pairs.emplace_back(PrimerPair{extResults.front()->primers_.at(ext.first),
																			   ligResults.front()->primers_.at(lig.first)});
					}
				}
			}else{
				if(!lig.second.reverseSrand_){
					if(ext.second.start_ > lig.second.end_ &&
							ext.second.end_ - lig.second.start_ > minLen &&
							ext.second.end_ - lig.second.start_ < maxLen){
						pairs.emplace_back(PrimerPair{extResults.front()->primers_.at(ext.first),
																			   ligResults.front()->primers_.at(lig.first)});
					}
				}
			}
		}
	}
	for(const auto & pair : pairs){
		auto region = pair.genRegion(extResults.front()->sequence_id_);
		out << region.genBedRecordCore().toDelimStr()
				<< "\t" << (region.reverseSrand_ ? pair.extPrimer_->originalSeq_ : pair.extPrimer_->fowardOrientationSeq_ )
				<< "\t" << (region.reverseSrand_ ? pair.ligPrimer_->fowardOrientationSeq_ : pair.ligPrimer_->originalSeq_)
				<< std::endl;
	}
	return 0;
}

int programWrapperRunner::parsePrimer3OutputToJson(const njh::progutils::CmdArgs & inputCommands){
	bfs::path input = "STDIN";
	OutOptions outOpts(bfs::path(""));
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processWritingOptions(outOpts);
	setUp.setOption(input, "--input", "input file or STDIN for reading for standard in");
	setUp.finishSetUp(std::cout);


	OutputStream out(outOpts);

	auto results = Primer3Results::parsePrimer3OutputResults(input);
	out << njh::json::toJson(results) << std::endl;

	return 0;
}

int programWrapperRunner::parsePrimer3OutputToBed(const njh::progutils::CmdArgs & inputCommands){
	bfs::path input = "STDIN";
	OutOptions outOpts(bfs::path(""));
	outOpts.outExtention_ = ".bed";
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processWritingOptions(outOpts);
	setUp.setOption(input, "--input", "input file or STDIN for reading for standard in");
	setUp.finishSetUp(std::cout);


	OutputStream out(outOpts);
	//addMetaBySampleName
	auto results = Primer3Results::parsePrimer3OutputResults(input);
	for(const auto & result : results){
		for(const auto & primer : result->primers_){
			out << result->sequence_id_
					<< "\t" << primer.second->forwardOrientationPos_.start_
					<< "\t" << primer.second->forwardOrientationPos_.start_ + primer.second->forwardOrientationPos_.size_
					<< "\t" << primer.second->name_
					<< "\t" << primer.second->forwardOrientationPos_.size_
					<< "\t" << (primer.second->right_ ? '-' : '+')
					<< std::endl;
		}
	}


	return 0;
}

int programWrapperRunner::findNonUniquePrimerArms(const njh::progutils::CmdArgs & inputCommands){
	bfs::path inDir = "";
	bfs::path outDir = "";
	std::string pattern = "";
	bool overWriteDir = false;
	uint32_t strainCutOff = 2;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(strainCutOff,   "--strainCutOff",   "Strain Cut Off");
	setUp.setOption(overWriteDir,   "--overWriteDir",   "Over Write Dir");
	setUp.setOption(inDir,   "--inDir",   "Input directory",  true);
	setUp.setOption(outDir,  "--outDir",  "Output directory", true);
	setUp.setOption(pattern, "--pattern", "Pattern",          true);
	setUp.finishSetUp(std::cout);

	njh::files::makeDir(njh::files::MkdirPar(outDir, overWriteDir));

	std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>> counts;
	auto allFiles = njh::files::listAllFiles(inDir, false, std::vector<std::regex>{std::regex{pattern}});

	for (const auto & f : allFiles) {
		if (f.second) {
			continue;
		}
		if (setUp.pars_.verbose_) {
			std::cout << f.first << std::endl;
		}
		BioDataFileIO<Bed6RecordCore> bedReader {
				IoOptions { InOptions { f.first } } };
		bedReader.openIn();
		Bed6RecordCore bRecord;
		while (bedReader.readNextRecord(bRecord)) {
			if(bRecord.extraFields_.size() < 2){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__<< ", error, extra fields should have at least two fields" << "\n";
				throw std::runtime_error{ss.str()};
			}
			std::string strainName = bRecord.chrom_.substr(0, bRecord.chrom_.find("."));
			++counts[bRecord.extraFields_[0] + "-" + bRecord.extraFields_[1]][strainName];
		}
	}
	if (setUp.pars_.debug_) {
		for(const auto & count : counts){
			if(count.second.size() > 1){
				std::cout << count.first << "\t" << count.second.size() << std::endl;
			}
		}
	}

//generatingPrime3TemplatesBasedOnMALN
	for (const auto & f : allFiles) {
		if (f.second) {
			continue;
		}
		if (setUp.pars_.verbose_) {
			std::cout << f.first << std::endl;
		}
		BioDataFileIO<Bed6RecordCore> bedReader {
				IoOptions { InOptions { f.first }, OutOptions(njh::files::make_path(outDir, bfs::basename(f.first))) } };
		bedReader.openIn();
		bedReader.openOut();
		Bed6RecordCore bRecord;
		while (bedReader.readNextRecord(bRecord)) {
			if(bRecord.extraFields_.size() < 2){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__<< ", error, extra fields should have at least two fields" << "\n";
				throw std::runtime_error{ss.str()};
			}
			std::string strainName = bRecord.chrom_.substr(0, bRecord.chrom_.find("."));
			if(counts[bRecord.extraFields_[0] + "-" + bRecord.extraFields_[1]][strainName] >= strainCutOff){
				bRecord.extraFields_.emplace_back(estd::to_string(counts[bRecord.extraFields_[0] + "-" + bRecord.extraFields_[1]].size()));
				bRecord.extraFields_.emplace_back(estd::to_string(counts[bRecord.extraFields_[0] + "-" + bRecord.extraFields_[1]][strainName]));
				bedReader.write(bRecord, [](const Bed6RecordCore & bRecord, std::ostream & out){
					out << bRecord.toDelimStrWithExtra() << std::endl;
				});
			}
		}
	}

	return 0;
}




} // namespace njhseq
