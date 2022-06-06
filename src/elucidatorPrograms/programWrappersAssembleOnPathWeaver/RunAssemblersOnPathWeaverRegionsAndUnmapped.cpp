//
// Created by Nicholas Hathaway on 3/19/22.
//

#include "elucidator/objects/seqContainers.h"
#include "elucidator/objects/BioDataObject.h"
#include "elucidator/seqToolsUtils/seqToolsUtils.hpp"


#include <njhseq/objects/seqObjects/seqKmers.h>
#include <njhseq/system.h>
#include <PathWeaver/objects/bam/RegionInvestigatorInBam.hpp>
#include <utility>

#include "programWrappersAssembleOnPathWeaverRunner.hpp"
#include "otherAssemblersUtils.hpp"


namespace njhseq {

class OtherAssemblersUtility{
public:
	struct InputPars{
		std::string regionUid_;
		bfs::path pwOutputDir_;
		std::string sample_;
		bfs::path bedFile_;
		std::string extraProgramOptions_;
		uint32_t reOrientingKmerLength_ = 9;
		uint32_t minFinalLength_ = 40;
		uint32_t numThreads_ = 1;
		std::string programName_;
		bfs::path outputDir_;

		void setPars(seqSetUp & setUp){
			setUp.setOption(bedFile_, "--bed", "The Regions to analyze", true);
			setUp.setOption(pwOutputDir_, "--pwOutputDir", "The PathWeaver directory", true);
			setUp.setOption(sample_, "--sample", "sample name", true);
			setUp.setOption(regionUid_, "--regionUid", "region ID", true);
			setUp.setOption(numThreads_, "--numThreads", "num Threads");
			setUp.setOption(extraProgramOptions_, "--extraProgramOptions", "extra program Options");
			setUp.setOption(minFinalLength_, "--minFinalLength", "min Final Length");
			setUp.setOption(reOrientingKmerLength_, "--reOrientingKmerLength", "reOrientingKmerLength Length");
			setUp.processDirectoryOutputName(njh::pasteAsStr(bfs::basename(pwOutputDir_), "_" + programName_ + "_TODAY"), true);
			outputDir_ = setUp.pars_.directoryName_;
		}
	};

	explicit OtherAssemblersUtility(InputPars pars): inputPars_(std::move(pars)){
		regInfo_ = std::make_shared<BamRegionInvestigator::RegionInfo>(GenomicRegion(Bed3RecordCore(inputPars_.regionUid_, 0, 1) ) );
		refFnp_ = njh::files::make_path(inputPars_.pwOutputDir_, "inputRegions.fasta");
		extractionFilesDir_ = njh::files::make_path(inputPars_.pwOutputDir_, "originalExtractionFiles");
		//first extract the reads
		pairedR1Fnp_ = njh::files::make_path(extractionFilesDir_, "allRaw_R1.fastq.gz");
		pairedR2Fnp_ = njh::files::make_path(extractionFilesDir_, "allRaw_R2.fastq.gz");
		singlesFnp_ =  njh::files::make_path(extractionFilesDir_, "allRaw.fastq.gz");
		{
			std::vector<uint32_t> readLens;
			if(bfs::exists(pairedR1Fnp_)){
				seqInfo seq;
				SeqInput reader(SeqIOOptions::genFastqInGz(pairedR1Fnp_));
				reader.openIn();
				while(reader.readNextRead(seq)){
					readLens.emplace_back(len(seq));
					++pairedReads_;
				}
			}
			if(bfs::exists(singlesFnp_)){
				seqInfo seq;
				SeqInput reader(SeqIOOptions::genFastqInGz(singlesFnp_));
				reader.openIn();
				while (reader.readNextRead(seq)) {
					readLens.emplace_back(len(seq));
					++singleReads_;
				}
			}
			medianReadLen_ = vectorMedianRef(readLens);
		}

		regInfo_->totalPairedReads_ = pairedReads_;
		regInfo_->totalReads_ = pairedReads_ + singleReads_;
		regInfo_->totalFinalReads_ = regInfo_->totalReads_;
		inputRegions_ = gatherRegions(inputPars_.bedFile_.string(), "", false);
		sortGRegionsByStart(inputRegions_);
		pw_finalPassDir_ = njh::files::make_path(inputPars_.pwOutputDir_, njh::pasteAsStr(inputPars_.sample_, "-finalPass"));
		finalPassDir_ = njh::files::make_path(inputPars_.outputDir_, inputPars_.sample_ + "-finalPass");
		filtStiched_pairedR1Fnp_ = njh::files::make_path(pw_finalPassDir_, "filteredExtractedPairs_R1.fastq");
		filtStiched_pairedR2Fnp_ = njh::files::make_path(pw_finalPassDir_, "filteredExtractedPairs_R2.fastq");
		filtStiched_singlesFnp_ =  njh::files::make_path(pw_finalPassDir_, "filteredSingles.fastq");
		outputFnp_ = njh::files::make_path(finalPassDir_, "output.fasta");
		outputAboveCutOffFnp_ = njh::files::make_path(finalPassDir_, "output_aboveCutOff.fasta");
		njh::files::makeDir(njh::files::MkdirPar(finalPassDir_));
	}

	InputPars inputPars_;
	std::shared_ptr<BamRegionInvestigator::RegionInfo> regInfo_;
	std::vector<GenomicRegion> inputRegions_;
	bfs::path refFnp_;
	bfs::path extractionFilesDir_;
	//first extract the reads
	bfs::path pairedR1Fnp_;
	bfs::path pairedR2Fnp_;
	bfs::path singlesFnp_;

	//the filtered stitched reads
	bfs::path filtStiched_pairedR1Fnp_;
	bfs::path filtStiched_pairedR2Fnp_;
	bfs::path filtStiched_singlesFnp_;

	uint32_t pairedReads_{0};
	uint32_t singleReads_{0};
	double medianReadLen_{1};
	bfs::path pw_finalPassDir_;

	bfs::path finalPassDir_;
	bfs::path outputFnp_;
	bfs::path outputAboveCutOffFnp_;

	uint32_t totalCount() const{
		return pairedReads_ + singleReads_;
	}

};


int programWrappersAssembleOnPathWeaverRunner::runMIRAOnPathWeaverRegionsAndUnmapped(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path MIRAOutDir = "MIRAOut";
	uint32_t miraAttempts = 3 ;
	seqSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processVerbose();
	OtherAssemblersUtility::InputPars inPars;
	inPars.programName_ = "MIRA";
	inPars.setPars(setUp);
	setUp.setOption(miraAttempts, "--miraAttempts", "mira Attempts");
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	njh::sys::requireExternalProgramThrow("mira");
	njh::sys::requireExternalProgramThrow("miraconvert");


	OtherAssemblersUtility utility(inPars);

	auto outputAboveCutOffSeqOpts = SeqIOOptions::genFastaOut(utility.outputAboveCutOffFnp_);
	SeqOutput outputAboveCutOffWriter(outputAboveCutOffSeqOpts);
	outputAboveCutOffWriter.openOut();

	std::string exceptionMess;

	auto regionOutputDir = njh::files::make_path(setUp.pars_.directoryName_, utility.inputPars_.regionUid_, utility.inputPars_.sample_);
	njh::files::makeDirP(njh::files::MkdirPar{regionOutputDir});

	try {
		if(!exists(utility.pairedR1Fnp_) && !exists(utility.singlesFnp_)){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", couldn't find " << utility.pairedR1Fnp_ << " or " << utility.singlesFnp_ << ", need to have at least one of them" << "\n";
			throw std::runtime_error{ss.str()};
		}
		std::stringstream MIRACmdStream;
		MIRACmdStream << "cd " << regionOutputDir;
		MIRACmdStream << " && mira ";

		{
			OutputStream miramanifestOutput(njh::files::make_path(regionOutputDir, "mira_manifest.txt"));
			miramanifestOutput << "project = " << MIRAOutDir.filename().string() << std::endl;
			miramanifestOutput << "job = genome,denovo,accurate"<< std::endl;


			miramanifestOutput << "parameters = -CO:force_nonIUPACconsensus_perseqtype=yes -GENERAL:number_of_threads=" << utility.inputPars_.numThreads_ << " COMMON_SETTINGS -NW:cmrnl=no -NW:cac=warn -NW:csrn=no -NW:cdrn=no"<< std::endl;
			//-EDIT:edit_homopolymer_overcalls=yes
			if(exists(utility.pairedR1Fnp_)){
				if(!exists(utility.pairedR2Fnp_)){
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", found: " << utility.pairedR1Fnp_ << " but couldn't find it's mate file: " << utility.pairedR2Fnp_ << "\n";
					throw std::runtime_error{ss.str()};
				}
				bfs::path pairedR1_app = njh::files::make_path(regionOutputDir, "appended_extracted_R1.fastq");
				bfs::path pairedR2_app = njh::files::make_path(regionOutputDir, "appended_extracted_R2.fastq");
				SeqOutput pairedR1_app_writer(SeqIOOptions::genFastqOut(pairedR1_app));
				pairedR1_app_writer.openOut();
				SeqOutput pairedR2_app_writer(SeqIOOptions::genFastqOut(pairedR2_app));
				pairedR2_app_writer.openOut();
				seqInfo seq;
				SeqInput pairedR1_reading(SeqIOOptions::genFastqIn(utility.pairedR1Fnp_));
				pairedR1_reading.openIn();
				while(pairedR1_reading.readNextRead(seq)){
					seq.name_.append("/1");
					pairedR1_app_writer.write(seq);
				}
				SeqInput pairedR2_reading(SeqIOOptions::genFastqIn(utility.pairedR2Fnp_));
				pairedR2_reading.openIn();
				while(pairedR2_reading.readNextRead(seq)){
					seq.name_.append("/2");
					pairedR2_app_writer.write(seq);
				}

				miramanifestOutput << "readgroup = " << utility.inputPars_.sample_ << "--" << utility.inputPars_.regionUid_ << std::endl;
				miramanifestOutput << "autopairing"<< std::endl;

				miramanifestOutput << "data = appended_extracted_R1.fastq appended_extracted_R2.fastq"<< std::endl;
				miramanifestOutput << "technology = solexa"<< std::endl;
				miramanifestOutput << "template_size = 50 1000 autorefine"<< std::endl;
				miramanifestOutput << "segment_placement = ---> <---"<< std::endl;
			}
			if(exists(utility.singlesFnp_)){
				miramanifestOutput << "readgroup = " << utility.inputPars_.sample_ << "--" << utility.inputPars_.regionUid_ << "-single" << std::endl;
				miramanifestOutput << "data = " << njh::files::normalize(utility.singlesFnp_).string() << std::endl;
				miramanifestOutput << "technology = solexa"<< std::endl;
			}
		}

		MIRACmdStream  << " mira_manifest.txt "
									 << " >> MIRARunLog_" << njh::getCurrentDate() << ".txt 2>&1";
		auto MIRAFullOutputDir = njh::files::make_path(regionOutputDir, MIRAOutDir.filename().string() + "_assembly");

		auto MIRARunOutput = njh::sys::run({MIRACmdStream.str()});
		uint32_t currentMiraAttempt = 1;
		while(!MIRARunOutput.success_ && currentMiraAttempt <= miraAttempts){
			MIRARunOutput = njh::sys::run({MIRACmdStream.str()});

			++currentMiraAttempt;
		}
		BioCmdsUtils::checkRunOutThrow(MIRARunOutput, __PRETTY_FUNCTION__);

		OutOptions MIRARunOutputLogOpts(njh::files::make_path(MIRAFullOutputDir, "MIRARunOutput.json"));
		OutputStream MIRARunOutputLogOut(MIRARunOutputLogOpts);
		MIRARunOutputLogOut << njh::json::toJson(MIRARunOutput) << std::endl;

		auto contigsFnp = njh::files::make_path(regionOutputDir, MIRAOutDir.filename().string() + "_assembly/", MIRAOutDir.filename().string() + "_d_results", MIRAOutDir.filename().string() + "_out.unpadded.fasta");


		auto contigsSeqIoOpts = SeqIOOptions::genFastaIn(contigsFnp);
//				contigsSeqIoOpts.includeWhiteSpaceInName_ = false;
		contigsSeqIoOpts.lowerCaseBases_ = "upper";
		SeqInput contigsReader(contigsSeqIoOpts);
		auto contigsSeqs = contigsReader.readAllReads<seqInfo>();
		std::vector<std::shared_ptr<seqWithKmerInfo>> contigsKmerReads;
		contigsKmerReads.reserve(contigsSeqs.size());
		for (const auto & seq : contigsSeqs) {
			contigsKmerReads.emplace_back(std::make_shared<seqWithKmerInfo>(seq));
		}
		allSetKmers(contigsKmerReads, utility.inputPars_.reOrientingKmerLength_, true);

		RefSeqsWithKmers refSeqs(utility.refFnp_, utility.inputPars_.reOrientingKmerLength_);
		readVec::reorientSeqs(contigsKmerReads, refSeqs.refKmerReads_);
		//sort by sequence length;
		njh::sort(contigsKmerReads, [](const std::shared_ptr<seqWithKmerInfo> & seq1, const std::shared_ptr<seqWithKmerInfo> & seq2){
			return len(seq1->seqBase_) > len(seq2->seqBase_);
		});

		OutOptions contigInfoOpts(njh::files::make_path(MIRAFullOutputDir, "contigs_outputInfo.tab.txt"));
		OutputStream contigInfoOut(contigInfoOpts);
		contigInfoOut << "name\tlength\tcoverage" << std::endl;

		for(const auto & contigsKmerRead : contigsKmerReads){
			auto assembleInfo = MIRAAssembleNameInfo(contigsKmerRead->seqBase_.name_);
			contigInfoOut << contigsKmerRead->seqBase_.name_
										<< "\t" << len(contigsKmerRead->seqBase_)
										<< "\t" << assembleInfo.coverage_ << std::endl;
		}
		auto reOrientedContigsFnp = njh::files::make_path(MIRAFullOutputDir, "reOriented_contigs.fasta");
		for(auto & seq : contigsKmerReads){
			auto assembleInfo = MIRAAssembleNameInfo(seq->seqBase_.name_);
			MetaDataInName seqMeta;
			seqMeta.addMeta("length", len(seq->seqBase_));
			seqMeta.addMeta("estimatedPerBaseCoverage", assembleInfo.coverage_);
			seqMeta.addMeta("regionUID", utility.inputPars_.regionUid_);
			seqMeta.addMeta("sample", utility.inputPars_.sample_);
			seqMeta.resetMetaInName(seq->seqBase_.name_);
			seq->seqBase_.cnt_ = assembleInfo.seqNumber_;
			seq->seqBase_.name_ += njh::pasteAsStr("_t", assembleInfo.seqNumber_);
		}
		SeqOutput::write(contigsKmerReads, SeqIOOptions::genFastaOut(reOrientedContigsFnp));
		SeqOutput::write(contigsKmerReads, SeqIOOptions::genFastaOut(utility.outputFnp_));

		std::vector<std::shared_ptr<seqWithKmerInfo>> finalSeqs = trimToFinalSeqs(contigsKmerReads, refSeqs);
		std::unordered_map<std::string, uint32_t> finalSeqCounts;
		for(const auto & seq : finalSeqs){
			++finalSeqCounts[seq->seqBase_.name_];
		}
		double totalCoverage = 0;
		for(auto & seq : finalSeqs){
			auto assembleInfo = MIRAAssembleNameInfo(seq->seqBase_.name_);
			totalCoverage += assembleInfo.coverage_;
		}
		std::unordered_map<std::string, uint32_t> finalSeqCountsWritten;
		for(auto & seq : finalSeqs){
			auto assembleInfo = MIRAAssembleNameInfo(seq->seqBase_.name_);
			MetaDataInName seqMeta;
			seqMeta.addMeta("trimmedLength", len(seq->seqBase_));
			seqMeta.addMeta("estimatedPerBaseCoverage", assembleInfo.coverage_);
			seqMeta.addMeta("trimStatus", seq->seqBase_.on_);
			seqMeta.addMeta("regionUID", utility.inputPars_.regionUid_);
			seqMeta.addMeta("sample", utility.inputPars_.sample_);
			if(finalSeqCounts[seq->seqBase_.name_] > 1){
				seqMeta.addMeta("seqTrimmedCount", finalSeqCountsWritten[seq->seqBase_.name_]);
				++finalSeqCountsWritten[seq->seqBase_.name_];
			}
			seqMeta.resetMetaInName(seq->seqBase_.name_);
			seq->seqBase_.frac_ = assembleInfo.coverage_/totalCoverage;
			seq->seqBase_.cnt_ = assembleInfo.seqNumber_;
			seq->seqBase_.name_ += njh::pasteAsStr("_t", assembleInfo.seqNumber_);
		}

		OutOptions trimmedContigInfoOpts(njh::files::make_path(MIRAFullOutputDir, "trimmed_reOriented_contigs_outputInfo.tab.txt"));
		OutputStream trimmedContigInfoOut(trimmedContigInfoOpts);
		trimmedContigInfoOut << "name\tlength\tcoverage" << std::endl;
		auto trimmedReOrientedContigsFnp = njh::files::make_path(MIRAFullOutputDir, "trimmed_reOriented_contigs.fasta");
		SeqOutput outputTrimmedWriter(SeqIOOptions::genFastaOut(trimmedReOrientedContigsFnp));
		auto trimmedReOrientedContigsFnp_belowCutOff = njh::files::make_path(MIRAFullOutputDir, "trimmed_reOriented_contigs_belowCutOff.fasta");
		SeqOutput belowCutOffOutputWriter(SeqIOOptions::genFastaOut(trimmedReOrientedContigsFnp_belowCutOff));

		uint32_t belowCutOff = 0;
		uint32_t aboveCutOff = 0;
		bool allPassTrim = true;
		for (const auto & contigsKmerRead : finalSeqs) {
			if (len(contigsKmerRead->seqBase_) < utility.inputPars_.minFinalLength_) {
				++belowCutOff;
				belowCutOffOutputWriter.openWrite(contigsKmerRead);
				contigsKmerRead->seqBase_.on_ = false;
			} else {
				MetaDataInName seqMeta(contigsKmerRead->seqBase_.name_);
				trimmedContigInfoOut << contigsKmerRead->seqBase_.name_
														 << "\t" << len(contigsKmerRead->seqBase_)
														 << "\t" << seqMeta.getMeta("estimatedPerBaseCoverage")
														 << std::endl;
				if(!contigsKmerRead->seqBase_.on_){
					allPassTrim = false;
				}else{
					++aboveCutOff;
				}
				outputAboveCutOffWriter.openWrite(contigsKmerRead);
				outputTrimmedWriter.openWrite(contigsKmerRead);
			}
		}
		if(allPassTrim){
			utility.regInfo_->infoCalled_ = true;
			utility.regInfo_->uniqHaps_ = aboveCutOff;
		}else{
			utility.regInfo_->infoCalled_ = false;
			utility.regInfo_->uniqHaps_ = 0;
		}
	} catch (std::exception & e) {
		exceptionMess = e.what();
		utility.regInfo_->infoCalled_ = false;
		utility.regInfo_->uniqHaps_ = 0;
	}

	outputAboveCutOffWriter.closeOut();
	OutputStream basicInfo(njh::files::make_path(utility.finalPassDir_, "basicInfoPerRegion.tab.txt"));

	basicInfo << "name\tsuccess\tuniqHaps\treadTotal\treadTotalUsed\ttotalPairedReads";
	basicInfo << "\tsample";
	basicInfo << "\n";
	basicInfo << utility.inputPars_.regionUid_;
	basicInfo << "\t" << njh::boolToStr(utility.regInfo_->infoCalled_)
						<< "\t" << utility.regInfo_->uniqHaps_
						<< "\t" << utility.regInfo_->totalReads_
						<< "\t" << utility.regInfo_->totalFinalReads_
						<< "\t" << utility.regInfo_->totalPairedReads_
						<< "\t" << utility.inputPars_.sample_;

	OutputStream exceptionsOut(njh::files::make_path(utility.finalPassDir_, "exceptionsMessages.tab.txt"));
	exceptionsOut << "regionUID\tmessage" << std::endl;
	exceptionsOut << utility.inputPars_.regionUid_ << "\t" << exceptionMess << std::endl;
	return 0;
}
//
int programWrappersAssembleOnPathWeaverRunner::runFermiLiteOnPathWeaverRegionsAndUnmapped(const njh::progutils::CmdArgs & inputCommands) {

	seqSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processVerbose();
	OtherAssemblersUtility::InputPars inPars;
	inPars.programName_ = "fermi-lite";
	inPars.extraProgramOptions_ = "";
	inPars.setPars(setUp);

	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	njh::sys::requireExternalProgramThrow("fml-asm");
	OtherAssemblersUtility utility(inPars);
	auto outputAboveCutOffSeqOpts = SeqIOOptions::genFastaOut(utility.outputAboveCutOffFnp_);
	SeqOutput outputAboveCutOffWriter(outputAboveCutOffSeqOpts);
	outputAboveCutOffWriter.openOut();

	std::string exceptionMess;

	auto regionOutputDir = njh::files::make_path(setUp.pars_.directoryName_, utility.inputPars_.regionUid_, utility.inputPars_.sample_);
	njh::files::makeDirP(njh::files::MkdirPar{regionOutputDir});

	try {
		if(!exists(utility.filtStiched_pairedR1Fnp_) && !exists(utility.filtStiched_singlesFnp_)){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", couldn't find " << utility.filtStiched_pairedR1Fnp_ << " or " << utility.filtStiched_singlesFnp_ << ", need to have at least one of them" << "\n";
			throw std::runtime_error{ss.str()};
		}

		//concatenate into 1 file
		std::vector<bfs::path> filesToCollapse;
		if(bfs::exists(utility.filtStiched_pairedR1Fnp_)){
			filesToCollapse.emplace_back(utility.filtStiched_pairedR1Fnp_);
			filesToCollapse.emplace_back(utility.filtStiched_pairedR2Fnp_);
		}
		if(bfs::exists(utility.filtStiched_singlesFnp_)){
			filesToCollapse.emplace_back(utility.filtStiched_singlesFnp_);
		}
		auto inputFnp = njh::files::make_path(regionOutputDir, "input.fastq.gz");
		auto outputFnp = njh::files::make_path(regionOutputDir, "raw_output.fastq");
		concatenateFiles(filesToCollapse, OutOptions(inputFnp));

		std::stringstream raw_fermiLiteCmdStream;
		raw_fermiLiteCmdStream << "cd " << regionOutputDir;
		raw_fermiLiteCmdStream << " && fml-asm ";

		raw_fermiLiteCmdStream  << " -t " << utility.inputPars_.numThreads_
														<< " " << utility.inputPars_.extraProgramOptions_
														<< " input.fastq.gz "
														<< " > " << "raw_output.fastq";
		std::string raw_fermiLiteCmd = raw_fermiLiteCmdStream.str();
		std::stringstream fermiLiteCmdStream;
		fermiLiteCmdStream << raw_fermiLiteCmd << " 2> fermiLiteRunLog_" << njh::getCurrentDate() << ".txt";
		const auto& fermiLiteFullOutputDir = regionOutputDir;

		auto fermiLiteRunOutput = njh::sys::run({fermiLiteCmdStream.str()});

		OutOptions fermiLiteRunOutputLogOpts(njh::files::make_path(fermiLiteFullOutputDir, "fermiLiteRunOutput.json"));
		OutputStream fermiLiteRunOutputLogOut(fermiLiteRunOutputLogOpts);
		fermiLiteRunOutputLogOut << njh::json::toJson(fermiLiteRunOutput) << std::endl;

		auto contigsFnp = njh::files::make_path(fermiLiteFullOutputDir, "raw_output.fastq");

		auto contigsSeqIoOpts = SeqIOOptions::genFastqIn(contigsFnp);
//				contigsSeqIoOpts.includeWhiteSpaceInName_ = false;
		contigsSeqIoOpts.lowerCaseBases_ = "upper";
		SeqInput contigsReader(contigsSeqIoOpts);
		auto contigsSeqs = contigsReader.readAllReads<seqInfo>();
		std::vector<std::shared_ptr<seqWithKmerInfo>> contigsKmerReads;
		contigsKmerReads.reserve(contigsSeqs.size());
		for (const auto & seq : contigsSeqs) {
			contigsKmerReads.emplace_back(std::make_shared<seqWithKmerInfo>(seq));
		}
		allSetKmers(contigsKmerReads, utility.inputPars_.reOrientingKmerLength_, true);

		RefSeqsWithKmers refSeqs(utility.refFnp_, utility.inputPars_.reOrientingKmerLength_);
		readVec::reorientSeqs(contigsKmerReads, refSeqs.refKmerReads_);
		//sort by sequence length;
		njh::sort(contigsKmerReads, [](const std::shared_ptr<seqWithKmerInfo> & seq1, const std::shared_ptr<seqWithKmerInfo> & seq2){
			return len(seq1->seqBase_) > len(seq2->seqBase_);
		});

		OutOptions contigInfoOpts(njh::files::make_path(fermiLiteFullOutputDir, "contigs_outputInfo.tab.txt"));
		OutputStream contigInfoOut(contigInfoOpts);
		contigInfoOut << "name\tlength\tcoverage" << std::endl;


		for( auto & contigsKmerRead : contigsKmerReads){
			auto assembleInfo = FermiLiteNameParse(contigsKmerRead->seqBase_.name_);
			contigsKmerRead->seqBase_.name_ = assembleInfo.modFullname_; //get rid of the \t characters in the name
			//assembleInfo.coverage_ = (assembleInfo.coverage_ * utility.medianReadLen_)/assembleInfo.len_;
			contigInfoOut << contigsKmerRead->seqBase_.name_
										<< "\t" << len(contigsKmerRead->seqBase_)
										<< "\t" << assembleInfo.coverage_ << std::endl;
		}
		auto reOrientedContigsFnp = njh::files::make_path(fermiLiteFullOutputDir, "reOriented_contigs.fasta");

		std::vector<std::shared_ptr<seqWithKmerInfo>> finalSeqs = trimToFinalSeqs(contigsKmerReads, refSeqs);
		std::unordered_map<std::string, uint32_t> finalSeqCounts;
		for(const auto & seq : finalSeqs){
			++finalSeqCounts[seq->seqBase_.name_];
		}

		double totalCoverage = 0;
		for(auto & seq : finalSeqs){
			auto assembleInfo = FermiLiteNameParse(seq->seqBase_.name_);
			assembleInfo.coverage_ = (assembleInfo.coverage_ * utility.medianReadLen_)/len(seq->seqBase_);
			totalCoverage += assembleInfo.coverage_;
		}

		for(auto & seq : contigsKmerReads){
			auto assembleInfo = FermiLiteNameParse(seq->seqBase_.name_);
			assembleInfo.coverage_ = (assembleInfo.coverage_ * utility.medianReadLen_)/len(seq->seqBase_);
			MetaDataInName seqMeta;
			seqMeta.addMeta("length", len(seq->seqBase_));
			seqMeta.addMeta("estimatedPerBaseCoverage", assembleInfo.coverage_);
			seqMeta.addMeta("regionUID", utility.inputPars_.regionUid_);
			seqMeta.addMeta("sample", utility.inputPars_.sample_);
			seqMeta.resetMetaInName(seq->seqBase_.name_);
			seq->seqBase_.cnt_ = (assembleInfo.coverage_/totalCoverage) * (utility.totalCount());
			seq->seqBase_.name_ += njh::pasteAsStr("_t", seq->seqBase_.cnt_);
		}
		SeqOutput::write(contigsKmerReads, SeqIOOptions::genFastaOut(reOrientedContigsFnp));
		SeqOutput::write(contigsKmerReads, SeqIOOptions::genFastaOut(utility.outputFnp_));

		std::unordered_map<std::string, uint32_t> finalSeqCountsWritten;
		for(auto & seq : finalSeqs){
			auto assembleInfo = FermiLiteNameParse(seq->seqBase_.name_);
			assembleInfo.coverage_ = (assembleInfo.coverage_ * utility.medianReadLen_)/len(seq->seqBase_);
			MetaDataInName seqMeta;
			seqMeta.addMeta("trimmedLength", len(seq->seqBase_));
			seqMeta.addMeta("estimatedPerBaseCoverage", assembleInfo.coverage_);
			seqMeta.addMeta("trimStatus", seq->seqBase_.on_);
			seqMeta.addMeta("regionUID", utility.inputPars_.regionUid_);
			seqMeta.addMeta("sample", utility.inputPars_.sample_);
			if(finalSeqCounts[seq->seqBase_.name_] > 1){
				seqMeta.addMeta("seqTrimmedCount", finalSeqCountsWritten[seq->seqBase_.name_]);
				++finalSeqCountsWritten[seq->seqBase_.name_];
			}
			seqMeta.resetMetaInName(seq->seqBase_.name_);
			seq->seqBase_.cnt_ = (assembleInfo.coverage_/totalCoverage) * (utility.totalCount());
			seq->seqBase_.name_ += njh::pasteAsStr("_t", seq->seqBase_.cnt_);
		}


		OutOptions trimmedContigInfoOpts(njh::files::make_path(fermiLiteFullOutputDir, "trimmed_reOriented_contigs_outputInfo.tab.txt"));
		OutputStream trimmedContigInfoOut(trimmedContigInfoOpts);
		trimmedContigInfoOut << "name\tlength\tcoverage" << std::endl;
		auto trimmedReOrientedContigsFnp = njh::files::make_path(fermiLiteFullOutputDir, "trimmed_reOriented_contigs.fasta");
		SeqOutput outputTrimmedWriter(SeqIOOptions::genFastaOut(trimmedReOrientedContigsFnp));
		auto trimmedReOrientedContigsFnp_belowCutOff = njh::files::make_path(fermiLiteFullOutputDir, "trimmed_reOriented_contigs_belowCutOff.fasta");
		SeqOutput belowCutOffOutputWriter(SeqIOOptions::genFastaOut(trimmedReOrientedContigsFnp_belowCutOff));

		uint32_t belowCutOff = 0;
		uint32_t aboveCutOff = 0;
		bool allPassTrim = true;
		for (const auto & contigsKmerRead : finalSeqs) {
			if (len(contigsKmerRead->seqBase_) < utility.inputPars_.minFinalLength_) {
				++belowCutOff;
				belowCutOffOutputWriter.openWrite(contigsKmerRead);
				contigsKmerRead->seqBase_.on_ = false;
			} else {
				MetaDataInName seqMeta(contigsKmerRead->seqBase_.name_);
				trimmedContigInfoOut << contigsKmerRead->seqBase_.name_
														 << "\t" << len(contigsKmerRead->seqBase_)
														 << "\t" << seqMeta.getMeta("estimatedPerBaseCoverage")
														 << std::endl;
				if(!contigsKmerRead->seqBase_.on_){
					allPassTrim = false;
				}else{
					++aboveCutOff;
				}
				outputAboveCutOffWriter.openWrite(contigsKmerRead);
				outputTrimmedWriter.openWrite(contigsKmerRead);
			}
		}
		if(allPassTrim){
			utility.regInfo_->infoCalled_ = true;
			utility.regInfo_->uniqHaps_ = aboveCutOff;
		}else{
			utility.regInfo_->infoCalled_ = false;
			utility.regInfo_->uniqHaps_ = 0;
		}
	} catch (std::exception & e) {
		exceptionMess = e.what();
		utility.regInfo_->infoCalled_ = false;
		utility.regInfo_->uniqHaps_ = 0;
	}


	outputAboveCutOffWriter.closeOut();
	OutputStream basicInfo(njh::files::make_path(utility.finalPassDir_, "basicInfoPerRegion.tab.txt"));

	basicInfo << "name\tsuccess\tuniqHaps\treadTotal\treadTotalUsed\ttotalPairedReads";
	basicInfo << "\tsample";
	basicInfo << "\n";
	basicInfo << utility.inputPars_.regionUid_;
	basicInfo << "\t" << njh::boolToStr(utility.regInfo_->infoCalled_)
						<< "\t" << utility.regInfo_->uniqHaps_
						<< "\t" << utility.regInfo_->totalReads_
						<< "\t" << utility.regInfo_->totalFinalReads_
						<< "\t" << utility.regInfo_->totalPairedReads_
						<< "\t" << utility.inputPars_.sample_;

	OutputStream exceptionsOut(njh::files::make_path(utility.finalPassDir_, "exceptionsMessages.tab.txt"));
	exceptionsOut << "regionUID\tmessage" << std::endl;
	exceptionsOut << utility.inputPars_.regionUid_ << "\t" << exceptionMess << std::endl;
	return 0;
}


int programWrappersAssembleOnPathWeaverRunner::runUnicyclerOnPathWeaverRegionsAndUnmapped(const njh::progutils::CmdArgs & inputCommands) {

	bfs::path unicyclerOutDir = "unicyclerOut";
	seqSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processVerbose();
	OtherAssemblersUtility::InputPars inPars;
	inPars.programName_ = "unicycler";
	inPars.extraProgramOptions_ = "--no_pilon";
	inPars.setPars(setUp);

	setUp.setOption(unicyclerOutDir,     "--unicyclerOutDir",     "unicycler Out Directory name, will be relative to final pass directory");
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	njh::sys::requireExternalProgramThrow("unicycler");
	OtherAssemblersUtility utility(inPars);
	auto outputAboveCutOffSeqOpts = SeqIOOptions::genFastaOut(utility.outputAboveCutOffFnp_);
	SeqOutput outputAboveCutOffWriter(outputAboveCutOffSeqOpts);
	outputAboveCutOffWriter.openOut();

	std::string exceptionMess;

	auto regionOutputDir = njh::files::make_path(setUp.pars_.directoryName_, utility.inputPars_.regionUid_, utility.inputPars_.sample_);
	njh::files::makeDirP(njh::files::MkdirPar{regionOutputDir});

	try {
		if(!exists(utility.pairedR1Fnp_) && !exists(utility.singlesFnp_)){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", couldn't find " << utility.pairedR1Fnp_ << " or " << utility.singlesFnp_ << ", need to have at least one of them" << "\n";
			throw std::runtime_error{ss.str()};
		}
		std::stringstream raw_unicyclerCmdStream;
		raw_unicyclerCmdStream << "cd " << regionOutputDir;
		raw_unicyclerCmdStream << " && unicycler ";

		if(exists(utility.pairedR1Fnp_)){
			if(!exists(utility.pairedR2Fnp_)){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", found: " << utility.pairedR1Fnp_ << " but couldn't find it's mate file: " << utility.pairedR2Fnp_ << "\n";
				throw std::runtime_error{ss.str()};
			}else{
				raw_unicyclerCmdStream << " -1 " << njh::files::normalize(utility.pairedR1Fnp_) << " -2 " << njh::files::normalize(utility.pairedR2Fnp_) << " ";
			}
		}
		if(exists(utility.singlesFnp_)){
			raw_unicyclerCmdStream << " -s  " << njh::files::normalize(utility.singlesFnp_);
		}
		uint32_t minRegionSize = std::numeric_limits<uint32_t>::max();

		for(const auto & region : utility.inputRegions_){
			if(region.getLen() < minRegionSize){
				minRegionSize = region.getLen();
			}
		}

		raw_unicyclerCmdStream  << " -t " << utility.inputPars_.numThreads_
														<< " " << utility.inputPars_.extraProgramOptions_
														<< " -o " << unicyclerOutDir << " ";
		if(minRegionSize < 1000){
			uint32_t minComponentSize = static_cast<uint32_t>(std::max(minRegionSize * .10, 1.0));
			raw_unicyclerCmdStream << " --min_fasta_length " << minComponentSize << " --min_component_size " << minComponentSize << " --min_dead_end_size " << minComponentSize << " ";
		}
		std::string raw_unicyclerCmd = raw_unicyclerCmdStream.str();
		std::stringstream unicyclerCmdStream;
		unicyclerCmdStream << raw_unicyclerCmd << " > unicyclerRunLog_" << njh::getCurrentDate() << ".txt 2>&1";
		auto unicyclerFullOutputDir = njh::files::make_path(regionOutputDir, unicyclerOutDir);

		auto unicyclerRunOutput = njh::sys::run({unicyclerCmdStream.str()});
		if(!unicyclerRunOutput.success_){
			std::stringstream unicyclerCmdStream_reAttempt;
			//sometimes unicycler fails at the spades correction step, so try re-runing without correcting
			unicyclerCmdStream_reAttempt << raw_unicyclerCmd << " --no_correct  > unicyclerReRunLog_" << njh::getCurrentDate() << ".txt 2>&1";
			if(bfs::exists(unicyclerFullOutputDir)){
				bfs::rename(unicyclerFullOutputDir, njh::files::make_path(regionOutputDir, "failed_" + unicyclerOutDir.string()));
			}
			auto unicyclerReRunOutput = njh::sys::run({unicyclerCmdStream_reAttempt.str()});
			if(!unicyclerReRunOutput.success_ && !unicyclerRunOutput.success_){
				//throw if they both failed
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error:\n" << unicyclerRunOutput.stdErr_ << "\n" << unicyclerReRunOutput.stdErr_ << "\n";
				throw std::runtime_error{ss.str()};
			}
		}


		OutOptions unicyclerRunOutputLogOpts(njh::files::make_path(unicyclerFullOutputDir, "unicyclerRunOutput.json"));
		OutputStream unicyclerRunOutputLogOut(unicyclerRunOutputLogOpts);
		unicyclerRunOutputLogOut << njh::json::toJson(unicyclerRunOutput) << std::endl;

		auto contigsFnp = njh::files::make_path(unicyclerFullOutputDir, "assembly.fasta");
		std::regex unicyclerNamePat{R"(^(\d+) length=(\d+) depth=([0-9.]+)x.*)"};

		auto contigsSeqIoOpts = SeqIOOptions::genFastaIn(contigsFnp);
//				contigsSeqIoOpts.includeWhiteSpaceInName_ = false;
		contigsSeqIoOpts.lowerCaseBases_ = "upper";
		SeqInput contigsReader(contigsSeqIoOpts);
		auto contigsSeqs = contigsReader.readAllReads<seqInfo>();
		std::vector<std::shared_ptr<seqWithKmerInfo>> contigsKmerReads;
		contigsKmerReads.reserve(contigsSeqs.size());
		for (const auto & seq : contigsSeqs) {
			contigsKmerReads.emplace_back(std::make_shared<seqWithKmerInfo>(seq));
		}
		allSetKmers(contigsKmerReads, utility.inputPars_.reOrientingKmerLength_, true);

		RefSeqsWithKmers refSeqs(utility.refFnp_, utility.inputPars_.reOrientingKmerLength_);
		readVec::reorientSeqs(contigsKmerReads, refSeqs.refKmerReads_);
		//sort by sequence length;
		njh::sort(contigsKmerReads, [](const std::shared_ptr<seqWithKmerInfo> & seq1, const std::shared_ptr<seqWithKmerInfo> & seq2){
			return len(seq1->seqBase_) > len(seq2->seqBase_);
		});

		OutOptions contigInfoOpts(njh::files::make_path(unicyclerFullOutputDir, "contigs_outputInfo.tab.txt"));
		OutputStream contigInfoOut(contigInfoOpts);
		contigInfoOut << "name\tlength\tcoverage" << std::endl;

		for(const auto & contigsKmerRead : contigsKmerReads){
			auto assembleInfo = DefaultAssembleNameInfo(contigsKmerRead->seqBase_.name_, unicyclerNamePat);
			assembleInfo.coverage_ = (assembleInfo.coverage_ * utility.totalCount() * utility.medianReadLen_)/assembleInfo.len_;
			contigInfoOut << contigsKmerRead->seqBase_.name_
										<< "\t" << len(contigsKmerRead->seqBase_)
										<< "\t" << assembleInfo.coverage_ << std::endl;
		}
		auto reOrientedContigsFnp = njh::files::make_path(unicyclerFullOutputDir, "reOriented_contigs.fasta");

		std::vector<std::shared_ptr<seqWithKmerInfo>> finalSeqs = trimToFinalSeqs(contigsKmerReads, refSeqs);
		std::unordered_map<std::string, uint32_t> finalSeqCounts;
		for(const auto & seq : finalSeqs){
			++finalSeqCounts[seq->seqBase_.name_];
		}

		double totalCoverage = 0;
		for(auto & seq : finalSeqs){
			auto assembleInfo = DefaultAssembleNameInfo(seq->seqBase_.name_, unicyclerNamePat);
			assembleInfo.coverage_ = (assembleInfo.coverage_ * utility.totalCount() * utility.medianReadLen_)/assembleInfo.len_;
			totalCoverage += assembleInfo.coverage_;
		}

		for(auto & seq : contigsKmerReads){
			auto assembleInfo = DefaultAssembleNameInfo(seq->seqBase_.name_, unicyclerNamePat);
			assembleInfo.coverage_ = (assembleInfo.coverage_ * utility.totalCount() * utility.medianReadLen_)/assembleInfo.len_;
			MetaDataInName seqMeta;
			seqMeta.addMeta("length", len(seq->seqBase_));
			seqMeta.addMeta("estimatedPerBaseCoverage", assembleInfo.coverage_);
			seqMeta.addMeta("regionUID", utility.inputPars_.regionUid_);
			seqMeta.addMeta("sample", utility.inputPars_.sample_);
			seqMeta.resetMetaInName(seq->seqBase_.name_);
			seq->seqBase_.cnt_ = (assembleInfo.coverage_/totalCoverage) * (utility.totalCount());
			seq->seqBase_.name_ += njh::pasteAsStr("_t", seq->seqBase_.cnt_);
		}
		SeqOutput::write(contigsKmerReads, SeqIOOptions::genFastaOut(reOrientedContigsFnp));
		SeqOutput::write(contigsKmerReads, SeqIOOptions::genFastaOut(utility.outputFnp_));

		std::unordered_map<std::string, uint32_t> finalSeqCountsWritten;
		for(auto & seq : finalSeqs){
			auto assembleInfo = DefaultAssembleNameInfo(seq->seqBase_.name_, unicyclerNamePat);
			assembleInfo.coverage_ = (assembleInfo.coverage_ * utility.totalCount() * utility.medianReadLen_)/assembleInfo.len_;
			MetaDataInName seqMeta;
			seqMeta.addMeta("trimmedLength", len(seq->seqBase_));
			seqMeta.addMeta("estimatedPerBaseCoverage", assembleInfo.coverage_);
			seqMeta.addMeta("trimStatus", seq->seqBase_.on_);
			seqMeta.addMeta("regionUID", utility.inputPars_.regionUid_);
			seqMeta.addMeta("sample", utility.inputPars_.sample_);
			if(finalSeqCounts[seq->seqBase_.name_] > 1){
				seqMeta.addMeta("seqTrimmedCount", finalSeqCountsWritten[seq->seqBase_.name_]);
				++finalSeqCountsWritten[seq->seqBase_.name_];
			}
			seqMeta.resetMetaInName(seq->seqBase_.name_);
			seq->seqBase_.cnt_ = (assembleInfo.coverage_/totalCoverage) * (utility.totalCount());
			seq->seqBase_.name_ += njh::pasteAsStr("_t", seq->seqBase_.cnt_);
		}


		OutOptions trimmedContigInfoOpts(njh::files::make_path(unicyclerFullOutputDir, "trimmed_reOriented_contigs_outputInfo.tab.txt"));
		OutputStream trimmedContigInfoOut(trimmedContigInfoOpts);
		trimmedContigInfoOut << "name\tlength\tcoverage" << std::endl;
		auto trimmedReOrientedContigsFnp = njh::files::make_path(unicyclerFullOutputDir, "trimmed_reOriented_contigs.fasta");
		SeqOutput outputTrimmedWriter(SeqIOOptions::genFastaOut(trimmedReOrientedContigsFnp));
		auto trimmedReOrientedContigsFnp_belowCutOff = njh::files::make_path(unicyclerFullOutputDir, "trimmed_reOriented_contigs_belowCutOff.fasta");
		SeqOutput belowCutOffOutputWriter(SeqIOOptions::genFastaOut(trimmedReOrientedContigsFnp_belowCutOff));

		uint32_t belowCutOff = 0;
		uint32_t aboveCutOff = 0;
		bool allPassTrim = true;
		for (const auto & contigsKmerRead : finalSeqs) {
			if (len(contigsKmerRead->seqBase_) < utility.inputPars_.minFinalLength_) {
				++belowCutOff;
				belowCutOffOutputWriter.openWrite(contigsKmerRead);
				contigsKmerRead->seqBase_.on_ = false;
			} else {
				MetaDataInName seqMeta(contigsKmerRead->seqBase_.name_);
				trimmedContigInfoOut << contigsKmerRead->seqBase_.name_
														 << "\t" << len(contigsKmerRead->seqBase_)
														 << "\t" << seqMeta.getMeta("estimatedPerBaseCoverage")
														 << std::endl;
				if(!contigsKmerRead->seqBase_.on_){
					allPassTrim = false;
				}else{
					++aboveCutOff;
				}
				outputAboveCutOffWriter.openWrite(contigsKmerRead);
				outputTrimmedWriter.openWrite(contigsKmerRead);
			}
		}
		if(allPassTrim){
			utility.regInfo_->infoCalled_ = true;
			utility.regInfo_->uniqHaps_ = aboveCutOff;
		}else{
			utility.regInfo_->infoCalled_ = false;
			utility.regInfo_->uniqHaps_ = 0;
		}
	} catch (std::exception & e) {
		exceptionMess = e.what();
		utility.regInfo_->infoCalled_ = false;
		utility.regInfo_->uniqHaps_ = 0;
	}


	outputAboveCutOffWriter.closeOut();
	OutputStream basicInfo(njh::files::make_path(utility.finalPassDir_, "basicInfoPerRegion.tab.txt"));

	basicInfo << "name\tsuccess\tuniqHaps\treadTotal\treadTotalUsed\ttotalPairedReads";
	basicInfo << "\tsample";
	basicInfo << "\n";
	basicInfo << utility.inputPars_.regionUid_;
	basicInfo << "\t" << njh::boolToStr(utility.regInfo_->infoCalled_)
						<< "\t" << utility.regInfo_->uniqHaps_
						<< "\t" << utility.regInfo_->totalReads_
						<< "\t" << utility.regInfo_->totalFinalReads_
						<< "\t" << utility.regInfo_->totalPairedReads_
						<< "\t" << utility.inputPars_.sample_ << std::endl;

	OutputStream exceptionsOut(njh::files::make_path(utility.finalPassDir_, "exceptionsMessages.tab.txt"));
	exceptionsOut << "regionUID\tmessage" << std::endl;
	exceptionsOut << utility.inputPars_.regionUid_ << "\t" << exceptionMess << std::endl;
	return 0;
}

int programWrappersAssembleOnPathWeaverRunner::runSpadesOnPathWeaverRegionsAndUnmapped(const njh::progutils::CmdArgs & inputCommands) {

	bool runMeta = false;
	bfs::path spadesOutDir = "spadesOut";
	seqSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processVerbose();
	OtherAssemblersUtility::InputPars inPars;
	inPars.programName_ = "spades";
	inPars.setPars(setUp);
	setUp.setOption(runMeta, "--runMeta", "Run Meta");
	if(runMeta){
		spadesOutDir = "metaspadesOut";
	}
	setUp.setOption(spadesOutDir,     "--spadesOutDir",     "spades Out Directory name, will be relative to final pass directory");


	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	njh::sys::requireExternalProgramThrow("spades.py");

	OtherAssemblersUtility utility(inPars);
	auto outputAboveCutOffSeqOpts = SeqIOOptions::genFastaOut(utility.outputAboveCutOffFnp_);
	SeqOutput outputAboveCutOffWriter(outputAboveCutOffSeqOpts);
	outputAboveCutOffWriter.openOut();

	std::string exceptionMess;

	auto regionOutputDir = njh::files::make_path(setUp.pars_.directoryName_, utility.inputPars_.regionUid_, utility.inputPars_.sample_);
	njh::files::makeDirP(njh::files::MkdirPar{regionOutputDir});

	try {
		if(!exists(utility.pairedR1Fnp_) && !exists(utility.singlesFnp_)){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", couldn't find " << utility.pairedR1Fnp_ << " or " << utility.singlesFnp_ << ", need to have at least one of them" << "\n";
			throw std::runtime_error{ss.str()};
		}
		std::stringstream spadesCmdStream;
		spadesCmdStream << "cd " << regionOutputDir;
		if(runMeta){
			spadesCmdStream << " && metaspades.py ";
		}else{
			spadesCmdStream << " && spades.py ";
		}

		if(exists(utility.pairedR1Fnp_)){
			if(!exists(utility.pairedR2Fnp_)){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", found: " << utility.pairedR1Fnp_ << " but couldn't find it's mate file: " << utility.pairedR2Fnp_ << "\n";
				throw std::runtime_error{ss.str()};
			}else{
				spadesCmdStream << " -1 " << njh::files::normalize(utility.pairedR1Fnp_) << " -2 " << njh::files::normalize(utility.pairedR2Fnp_) << " ";
			}
		}
		if(exists(utility.singlesFnp_)){
			spadesCmdStream << " -s  " << njh::files::normalize(utility.singlesFnp_);
		}
		spadesCmdStream  << " -t " << utility.inputPars_.numThreads_
										 << " " << utility.inputPars_.extraProgramOptions_
										 << " -o " << spadesOutDir
										 << " > spadesRunLog_" << njh::getCurrentDate() << ".txt 2>&1";
		auto spadesFullOutputDir = njh::files::make_path(regionOutputDir, spadesOutDir);

		auto spadesRunOutput = njh::sys::run({spadesCmdStream.str()});

		BioCmdsUtils::checkRunOutThrow(spadesRunOutput, __PRETTY_FUNCTION__);

		OutOptions spadesRunOutputLogOpts(njh::files::make_path(spadesFullOutputDir, "spadesRunOutput.json"));
		OutputStream spadesRunOutputLogOut(spadesRunOutputLogOpts);
		spadesRunOutputLogOut << njh::json::toJson(spadesRunOutput) << std::endl;

		auto contigsFnp = njh::files::make_path(spadesFullOutputDir, "contigs.fasta");


		auto contigsSeqIoOpts = SeqIOOptions::genFastaIn(contigsFnp);
//				contigsSeqIoOpts.includeWhiteSpaceInName_ = false;
		contigsSeqIoOpts.lowerCaseBases_ = "upper";
		SeqInput contigsReader(contigsSeqIoOpts);
		auto contigsSeqs = contigsReader.readAllReads<seqInfo>();
		std::vector<std::shared_ptr<seqWithKmerInfo>> contigsKmerReads;
		contigsKmerReads.reserve(contigsSeqs.size());
		for (const auto & seq : contigsSeqs) {
			contigsKmerReads.emplace_back(std::make_shared<seqWithKmerInfo>(seq));
		}
		allSetKmers(contigsKmerReads, utility.inputPars_.reOrientingKmerLength_, true);

		RefSeqsWithKmers refSeqs(utility.refFnp_, utility.inputPars_.reOrientingKmerLength_);
		readVec::reorientSeqs(contigsKmerReads, refSeqs.refKmerReads_);
		//sort by sequence length;
		njh::sort(contigsKmerReads, [](const std::shared_ptr<seqWithKmerInfo> & seq1, const std::shared_ptr<seqWithKmerInfo> & seq2){
			return len(seq1->seqBase_) > len(seq2->seqBase_);
		});

		OutOptions contigInfoOpts(njh::files::make_path(spadesFullOutputDir, "contigs_outputInfo.tab.txt"));
		OutputStream contigInfoOut(contigInfoOpts);
		contigInfoOut << "name\tlength\tcoverage" << std::endl;

		for(const auto & contigsKmerRead : contigsKmerReads){
			auto assembleInfo = DefaultAssembleNameInfo(contigsKmerRead->seqBase_.name_);
			contigInfoOut << contigsKmerRead->seqBase_.name_
										<< "\t" << len(contigsKmerRead->seqBase_)
										<< "\t" << assembleInfo.coverage_ << std::endl;
		}
		auto reOrientedContigsFnp = njh::files::make_path(spadesFullOutputDir, "reOriented_contigs.fasta");


		//readVecTrimmer::trimSeqToRefByGlobalAln(contigsKmerReads, refSeqs, refSeqsKmerInfos, alignerObj);
		std::vector<std::shared_ptr<seqWithKmerInfo>> finalSeqs = trimToFinalSeqs(contigsKmerReads, refSeqs);
		std::unordered_map<std::string, uint32_t> finalSeqCounts;
		for(const auto & seq : finalSeqs){
			++finalSeqCounts[seq->seqBase_.name_];
		}
		double totalCoverage = 0;
		for(auto & seq : finalSeqs){
			auto assembleInfo = DefaultAssembleNameInfo(seq->seqBase_.name_);
			totalCoverage += assembleInfo.coverage_;
		}

		for(auto & seq : contigsKmerReads){
			auto assembleInfo = DefaultAssembleNameInfo(seq->seqBase_.name_);
			MetaDataInName seqMeta;
			seqMeta.addMeta("length", len(seq->seqBase_));
			seqMeta.addMeta("estimatedPerBaseCoverage", assembleInfo.coverage_);
			seqMeta.addMeta("regionUID", utility.inputPars_.regionUid_);
			seqMeta.addMeta("sample", utility.inputPars_.sample_);
			seqMeta.resetMetaInName(seq->seqBase_.name_);
			seq->seqBase_.cnt_ = (assembleInfo.coverage_/totalCoverage) * (utility.totalCount());
			seq->seqBase_.name_ += njh::pasteAsStr("_t", seq->seqBase_.cnt_);
		}
		SeqOutput::write(contigsKmerReads, SeqIOOptions::genFastaOut(reOrientedContigsFnp));
		SeqOutput::write(contigsKmerReads, SeqIOOptions::genFastaOut(utility.outputFnp_));

		std::unordered_map<std::string, uint32_t> finalSeqCountsWritten;
		for(auto & seq : finalSeqs){
			auto assembleInfo = DefaultAssembleNameInfo(seq->seqBase_.name_);
			MetaDataInName seqMeta;
			seqMeta.addMeta("trimmedLength", len(seq->seqBase_));
			seqMeta.addMeta("estimatedPerBaseCoverage", assembleInfo.coverage_);
			seqMeta.addMeta("trimStatus", seq->seqBase_.on_);
			seqMeta.addMeta("regionUID", utility.inputPars_.regionUid_);
			seqMeta.addMeta("sample", utility.inputPars_.sample_);
			if(finalSeqCounts[seq->seqBase_.name_] > 1){
				seqMeta.addMeta("seqTrimmedCount", finalSeqCountsWritten[seq->seqBase_.name_]);
				++finalSeqCountsWritten[seq->seqBase_.name_];
			}
			seqMeta.resetMetaInName(seq->seqBase_.name_);
			seq->seqBase_.cnt_ = (assembleInfo.coverage_/totalCoverage) * (utility.totalCount());
			seq->seqBase_.name_ += njh::pasteAsStr("_t", seq->seqBase_.cnt_);
		}

		OutOptions trimmedContigInfoOpts(njh::files::make_path(spadesFullOutputDir, "trimmed_reOriented_contigs_outputInfo.tab.txt"));
		OutputStream trimmedContigInfoOut(trimmedContigInfoOpts);
		trimmedContigInfoOut << "name\tlength\tcoverage" << std::endl;
		auto trimmedReOrientedContigsFnp = njh::files::make_path(spadesFullOutputDir, "trimmed_reOriented_contigs.fasta");
		SeqOutput outputTrimmedWriter(SeqIOOptions::genFastaOut(trimmedReOrientedContigsFnp));
		auto trimmedReOrientedContigsFnp_belowCutOff = njh::files::make_path(spadesFullOutputDir, "trimmed_reOriented_contigs_belowCutOff.fasta");
		SeqOutput belowCutOffOutputWriter(SeqIOOptions::genFastaOut(trimmedReOrientedContigsFnp_belowCutOff));

		uint32_t belowCutOff = 0;
		uint32_t aboveCutOff = 0;
		bool allPassTrim = true;
		for (const auto & contigsKmerRead : finalSeqs) {
			if (len(contigsKmerRead->seqBase_) < utility.inputPars_.minFinalLength_) {
				++belowCutOff;
				belowCutOffOutputWriter.openWrite(contigsKmerRead);
				contigsKmerRead->seqBase_.on_ = false;
			} else {
				MetaDataInName seqMeta(contigsKmerRead->seqBase_.name_);
				trimmedContigInfoOut << contigsKmerRead->seqBase_.name_
														 << "\t" << len(contigsKmerRead->seqBase_)
														 << "\t" << seqMeta.getMeta("estimatedPerBaseCoverage")
														 << std::endl;
				if(!contigsKmerRead->seqBase_.on_){
					allPassTrim = false;
				}else{
					++aboveCutOff;
				}
				outputAboveCutOffWriter.openWrite(contigsKmerRead);
				outputTrimmedWriter.openWrite(contigsKmerRead);
			}
		}
		if(allPassTrim){
			utility.regInfo_->infoCalled_ = true;
			utility.regInfo_->uniqHaps_ = aboveCutOff;
		}else{
			utility.regInfo_->infoCalled_ = false;
			utility.regInfo_->uniqHaps_ = 0;
		}
	} catch (std::exception & e) {
		exceptionMess = e.what();
		utility.regInfo_->infoCalled_ = false;
		utility.regInfo_->uniqHaps_ = 0;
	}

	outputAboveCutOffWriter.closeOut();
	OutputStream basicInfo(njh::files::make_path(utility.finalPassDir_, "basicInfoPerRegion.tab.txt"));

	basicInfo << "name\tsuccess\tuniqHaps\treadTotal\treadTotalUsed\ttotalPairedReads";
	basicInfo << "\tsample";
	basicInfo << "\n";
	basicInfo << utility.inputPars_.regionUid_;
	basicInfo << "\t" << njh::boolToStr(utility.regInfo_->infoCalled_)
						<< "\t" << utility.regInfo_->uniqHaps_
						<< "\t" << utility.regInfo_->totalReads_
						<< "\t" << utility.regInfo_->totalFinalReads_
						<< "\t" << utility.regInfo_->totalPairedReads_
						<< "\t" << utility.inputPars_.sample_ << std::endl;

	OutputStream exceptionsOut(njh::files::make_path(utility.finalPassDir_, "exceptionsMessages.tab.txt"));
	exceptionsOut << "regionUID\tmessage" << std::endl;
	exceptionsOut << utility.inputPars_.regionUid_ << "\t" << exceptionMess << std::endl;
	return 0;
}

//


//


int programWrappersAssembleOnPathWeaverRunner::runRayOnPathWeaverRegionsAndUnmapped(const njh::progutils::CmdArgs & inputCommands) {

	uint32_t RayKmerLength = std::numeric_limits<uint32_t>::max();
	bfs::path RayOutDir = "RayOut";
	seqSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processVerbose();
	OtherAssemblersUtility::InputPars inPars;
	inPars.programName_ = "Ray";
	inPars.setPars(setUp);
	setUp.setOption(RayKmerLength, "--RayKmerLength", "Ray Kmer Length");



	//setUp.setOption(RayNumThreads, "--RayNumThreads", "Ray Num Threads");


	setUp.setOption(RayOutDir,     "--RayOutDir",     "Ray Out Directory name, will be relative to final pass directory");


	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	njh::sys::requireExternalProgramThrow("Ray");

	OtherAssemblersUtility utility(inPars);

	auto outputAboveCutOffSeqOpts = SeqIOOptions::genFastaOut(utility.outputAboveCutOffFnp_);
	SeqOutput outputAboveCutOffWriter(outputAboveCutOffSeqOpts);
	outputAboveCutOffWriter.openOut();

	std::string exceptionMess;

	auto regionOutputDir = njh::files::make_path(setUp.pars_.directoryName_, utility.inputPars_.regionUid_, utility.inputPars_.sample_);
	njh::files::makeDirP(njh::files::MkdirPar{regionOutputDir});


	try {
		bfs::path optimJsonFnp = njh::files::make_path(utility.inputPars_.pwOutputDir_, utility.inputPars_.sample_ + "-finalPass", "optimizationInfoBest.json");
		Json::Value optimJson = njh::json::parseFile(optimJsonFnp.string());
		if(std::numeric_limits<uint32_t>::max() == RayKmerLength){
			RayKmerLength = optimJson["runParams_"]["klen_"].asUInt64() ;
		}
//				std::cout << njh::conToStr(optimJson.getMemberNames(), "\n") << std::endl;
//
//				std::cout << "optimJson[\"runParams_\"][\"klen_\"].asUInt64(): " << optimJson["runParams_"]["klen_"].asUInt64() << std::endl;
//
//				exit(1);
		if(!exists(utility.pairedR1Fnp_) && !exists(utility.singlesFnp_)){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", couldn't find " << utility.pairedR1Fnp_ << " or " << utility.singlesFnp_ << ", need to have at least one of them" << "\n";
			throw std::runtime_error{ss.str()};
		}
		std::stringstream RayCmdStream;
		RayCmdStream << "cd " << regionOutputDir << " && Ray ";
		if(!exists(utility.pairedR1Fnp_)){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "Ray requires paired reads"<< "\n";
			throw std::runtime_error{ss.str()};
		}
		//unzip files since Ray is the worst
		{
			concatenateFiles({utility.pairedR1Fnp_}, njh::files::make_path(regionOutputDir, utility.pairedR1Fnp_.filename().replace_extension("")));
			concatenateFiles({utility.pairedR2Fnp_}, njh::files::make_path(regionOutputDir, utility.pairedR2Fnp_.filename().replace_extension("")));
		}
//
//
		RayCmdStream    << " -k " << RayKmerLength
										<< " -p " << utility.pairedR1Fnp_.filename().replace_extension("") << " " <<  utility.pairedR2Fnp_.filename().replace_extension("")
										<< " " << utility.inputPars_.extraProgramOptions_
										<< " -o " << RayOutDir
										<< " > RayRunLog_" << njh::getCurrentDate() << ".txt 2>&1";
		auto RayFullOutputDir = njh::files::make_path(regionOutputDir, RayOutDir);

		auto RayRunOutput = njh::sys::run({RayCmdStream.str()});

		BioCmdsUtils::checkRunOutThrow(RayRunOutput, __PRETTY_FUNCTION__);

		OutOptions RayRunOutputLogOpts(njh::files::make_path(RayFullOutputDir, "RayRunOutput.json"));
		OutputStream RayRunOutputLogOut(RayRunOutputLogOpts);
		RayRunOutputLogOut << njh::json::toJson(RayRunOutput) << std::endl;

		auto contigsFnp = njh::files::make_path(RayFullOutputDir, "Contigs.fasta");


		std::unordered_map<std::string, double> kmerCoverage;
		auto coverageInfo = njh::files::make_path(RayFullOutputDir, "BiologicalAbundances/_DeNovoAssembly/Contigs.tsv");
		table covTab(coverageInfo, "\t", true);
		for(const auto & row : covTab){
			kmerCoverage[row[covTab.getColPos("#Contig name")]] = njh::StrToNumConverter::stoToNum<double>(row[covTab.getColPos("Mode k-mer coverage depth")]);
			kmerCoverage[row[covTab.getColPos("#Contig name")] + "_Comp"] = njh::StrToNumConverter::stoToNum<double>(row[covTab.getColPos("Mode k-mer coverage depth")]);
		}
		auto contigsSeqIoOpts = SeqIOOptions::genFastaIn(contigsFnp);
		contigsSeqIoOpts.includeWhiteSpaceInName_ = false;
		contigsSeqIoOpts.lowerCaseBases_ = "upper";
		SeqInput contigsReader(contigsSeqIoOpts);
		auto contigsSeqs = contigsReader.readAllReads<seqInfo>();
		std::vector<std::shared_ptr<seqWithKmerInfo>> contigsKmerReads;
		contigsKmerReads.reserve(contigsSeqs.size());
		for (const auto & seq : contigsSeqs) {
			contigsKmerReads.emplace_back(std::make_shared<seqWithKmerInfo>(seq));
		}
		allSetKmers(contigsKmerReads, utility.inputPars_.reOrientingKmerLength_, true);

		RefSeqsWithKmers refSeqs(utility.refFnp_, utility.inputPars_.reOrientingKmerLength_);
		readVec::reorientSeqs(contigsKmerReads, refSeqs.refKmerReads_);
		//sort by sequence length;
		njh::sort(contigsKmerReads, [](const std::shared_ptr<seqWithKmerInfo> & seq1, const std::shared_ptr<seqWithKmerInfo> & seq2){
			return len(seq1->seqBase_) > len(seq2->seqBase_);
		});

		OutOptions contigInfoOpts(njh::files::make_path(RayFullOutputDir, "contigs_outputInfo.tab.txt"));
		OutputStream contigInfoOut(contigInfoOpts);
		contigInfoOut << "name\tlength\tcoverage" << std::endl;

		for(const auto & contigsKmerRead : contigsKmerReads){
			//auto assembleInfo = DefaultAssembleNameInfo(contigsKmerRead->seqBase_.name_, true);
			contigInfoOut << contigsKmerRead->seqBase_.name_
										<< "\t" << len(contigsKmerRead->seqBase_)
										<< "\t" << kmerCoverage[contigsKmerRead->seqBase_.name_] << std::endl;
		}
		auto reOrientedContigsFnp = njh::files::make_path(RayFullOutputDir, "reOriented_contigs.fasta");


		std::vector<std::shared_ptr<seqWithKmerInfo>> finalSeqs = trimToFinalSeqs(contigsKmerReads, refSeqs);
		std::unordered_map<std::string, uint32_t> finalSeqCounts;
		for(const auto & seq : finalSeqs){
			++finalSeqCounts[seq->seqBase_.name_];
		}

		double totalCoverage = 0;
		for(auto & seq : finalSeqs){
			//auto assembleInfo = DefaultAssembleNameInfo(seq->seqBase_.name_, true);
			totalCoverage += kmerCoverage[seq->seqBase_.name_];
		}
		for(auto & seq : contigsKmerReads){
			//auto assembleInfo = DefaultAssembleNameInfo(seq->seqBase_.name_);
			MetaDataInName seqMeta;
			seqMeta.addMeta("length", len(seq->seqBase_));
			seqMeta.addMeta("estimatedPerBaseCoverage", kmerCoverage[seq->seqBase_.name_]);
			seqMeta.addMeta("regionUID", utility.inputPars_.regionUid_);
			seqMeta.addMeta("sample", utility.inputPars_.sample_);


//			std::cout << "utility.totalCount(): " << utility.totalCount() << std::endl;
//			std::cout << "totalCoverage: " << totalCoverage << std::endl;
//			std::cout << "kmerCoverage[seq->seqBase_.name_]: " << kmerCoverage[seq->seqBase_.name_] << std::endl;
//			std::cout << "kmerCoverage[seq->seqBase_.name_]/totalCoverage: " << kmerCoverage[seq->seqBase_.name_]/totalCoverage << std::endl;
			seq->seqBase_.cnt_ = (kmerCoverage[seq->seqBase_.name_]/totalCoverage) * (utility.totalCount());
			std::string oldName = seq->seqBase_.name_;
			seqMeta.resetMetaInName(seq->seqBase_.name_);
			kmerCoverage[seq->seqBase_.name_] = kmerCoverage[oldName];
			seq->seqBase_.name_ += njh::pasteAsStr("_t", seq->seqBase_.cnt_);

		}
		SeqOutput::write(contigsKmerReads, SeqIOOptions::genFastaOut(reOrientedContigsFnp));
		SeqOutput::write(contigsKmerReads, SeqIOOptions::genFastaOut(utility.outputFnp_));

		std::unordered_map<std::string, uint32_t> finalSeqCountsWritten;
		for(auto & seq : finalSeqs){
			//auto assembleInfo = DefaultAssembleNameInfo(seq->seqBase_.name_, true);
			MetaDataInName seqMeta;
			seqMeta.addMeta("trimmedLength", len(seq->seqBase_));
			seqMeta.addMeta("estimatedPerBaseCoverage", kmerCoverage[seq->seqBase_.name_]);
			seqMeta.addMeta("trimStatus", seq->seqBase_.on_);
			seqMeta.addMeta("regionUID", utility.inputPars_.regionUid_);
			seqMeta.addMeta("sample", utility.inputPars_.sample_);
			if(finalSeqCounts[seq->seqBase_.name_] > 1){
				seqMeta.addMeta("seqTrimmedCount", finalSeqCountsWritten[seq->seqBase_.name_]);
				++finalSeqCountsWritten[seq->seqBase_.name_];
			}
			//seq->seqBase_.cnt_ = (defaultCoverage/totalCoverage) * (readCounts.pairedReads_ + readCounts.unpaiedReads_ + readCounts.orphans_);
//					std::cout << "kmerCoverage[seq->seqBase_.name_]: " << kmerCoverage[seq->seqBase_.name_] << std::endl;
//					std::cout << "totalCoverage: " << totalCoverage << std::endl;
//					std::cout << "reads: " << readCounts.pairedReads_ + readCounts.unpaiedReads_ + readCounts.orphans_ << std::endl;

			seq->seqBase_.cnt_ = (kmerCoverage[seq->seqBase_.name_]/totalCoverage) * (utility.totalCount());
			std::string oldName = seq->seqBase_.name_;
			seqMeta.resetMetaInName(seq->seqBase_.name_);
			kmerCoverage[seq->seqBase_.name_] = kmerCoverage[oldName];
			seq->seqBase_.name_ += njh::pasteAsStr("_t", seq->seqBase_.cnt_);
		}

		OutOptions trimmedContigInfoOpts(njh::files::make_path(RayFullOutputDir, "trimmed_reOriented_contigs_outputInfo.tab.txt"));
		OutputStream trimmedContigInfoOut(trimmedContigInfoOpts);
		trimmedContigInfoOut << "name\tlength\tcoverage" << std::endl;
		auto trimmedReOrientedContigsFnp = njh::files::make_path(RayFullOutputDir, "trimmed_reOriented_contigs.fasta");
		SeqOutput outputTrimmedWriter(SeqIOOptions::genFastaOut(trimmedReOrientedContigsFnp));
		auto trimmedReOrientedContigsFnp_belowCutOff = njh::files::make_path(RayFullOutputDir, "trimmed_reOriented_contigs_belowCutOff.fasta");
		SeqOutput belowCutOffOutputWriter(SeqIOOptions::genFastaOut(trimmedReOrientedContigsFnp_belowCutOff));

		uint32_t belowCutOff = 0;
		uint32_t aboveCutOff = 0;
		bool allPassTrim = true;
		for (const auto & contigsKmerRead : finalSeqs) {
			if (len(contigsKmerRead->seqBase_) < utility.inputPars_.minFinalLength_) {
				++belowCutOff;
				belowCutOffOutputWriter.openWrite(contigsKmerRead);
				contigsKmerRead->seqBase_.on_ = false;
			} else {
				MetaDataInName seqMeta(contigsKmerRead->seqBase_.name_);
				trimmedContigInfoOut << contigsKmerRead->seqBase_.name_
														 << "\t" << len(contigsKmerRead->seqBase_)
														 << "\t" << seqMeta.getMeta("estimatedPerBaseCoverage")
														 << std::endl;
				if(!contigsKmerRead->seqBase_.on_){
					allPassTrim = false;
				}else{
					++aboveCutOff;
				}
				outputAboveCutOffWriter.openWrite(contigsKmerRead);
				outputTrimmedWriter.openWrite(contigsKmerRead);
			}
		}
		if(allPassTrim){
			utility.regInfo_->infoCalled_ = true;
			utility.regInfo_->uniqHaps_ = aboveCutOff;
		}else{
			utility.regInfo_->infoCalled_ = false;
			utility.regInfo_->uniqHaps_ = 0;
		}
	} catch (std::exception & e) {
		exceptionMess = e.what();
		utility.regInfo_->infoCalled_ = false;
		utility.regInfo_->uniqHaps_ = 0;
	}

	outputAboveCutOffWriter.closeOut();
	OutputStream basicInfo(njh::files::make_path(utility.finalPassDir_, "basicInfoPerRegion.tab.txt"));

	basicInfo << "name\tsuccess\tuniqHaps\treadTotal\treadTotalUsed\ttotalPairedReads";
	basicInfo << "\tsample";
	basicInfo << "\n";
	basicInfo << utility.inputPars_.regionUid_;
	basicInfo << "\t" << njh::boolToStr(utility.regInfo_->infoCalled_)
						<< "\t" << utility.regInfo_->uniqHaps_
						<< "\t" << utility.regInfo_->totalReads_
						<< "\t" << utility.regInfo_->totalFinalReads_
						<< "\t" << utility.regInfo_->totalPairedReads_
						<< "\t" << utility.inputPars_.sample_ << std::endl;

	OutputStream exceptionsOut(njh::files::make_path(utility.finalPassDir_, "exceptionsMessages.tab.txt"));
	exceptionsOut << "regionUID\tmessage" << std::endl;
	exceptionsOut << utility.inputPars_.regionUid_ << "\t" << exceptionMess << std::endl;
	return 0;
}


int programWrappersAssembleOnPathWeaverRunner::runIDBAUDOnPathWeaverRegionsAndUnmapped(const njh::progutils::CmdArgs & inputCommands) {



	uint32_t IDBAUDKmerLength = std::numeric_limits<uint32_t>::max();
	bfs::path IDBAUDOutDir = "IDBAUDOut";
	seqSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processVerbose();
	OtherAssemblersUtility::InputPars inPars;
	inPars.programName_ = "IDBAUD";
	inPars.setPars(setUp);
	setUp.setOption(IDBAUDKmerLength, "--IDBAUDKmerLength", "IDBAUD Kmer Length");

	setUp.setOption(IDBAUDOutDir,     "--IDBAUDOutDir",     "IDBAUD Out Directory name, will be relative to final pass directory");


	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	njh::sys::requireExternalProgramThrow("fq2fa");
	njh::sys::requireExternalProgramThrow("idba_ud");
	OtherAssemblersUtility utility(inPars);

	auto outputAboveCutOffSeqOpts = SeqIOOptions::genFastaOut(utility.outputAboveCutOffFnp_);
	SeqOutput outputAboveCutOffWriter(outputAboveCutOffSeqOpts);
	outputAboveCutOffWriter.openOut();

	std::string exceptionMess;

	auto regionOutputDir = njh::files::make_path(setUp.pars_.directoryName_, utility.inputPars_.regionUid_, utility.inputPars_.sample_);
	njh::files::makeDirP(njh::files::MkdirPar{regionOutputDir});


	try {
		bfs::path optimJsonFnp = njh::files::make_path(utility.inputPars_.pwOutputDir_, utility.inputPars_.sample_ + "-finalPass", "optimizationInfoBest.json");
		Json::Value optimJson = njh::json::parseFile(optimJsonFnp.string());
		if(std::numeric_limits<uint32_t>::max() == IDBAUDKmerLength){
			IDBAUDKmerLength = optimJson["runParams_"]["klen_"].asUInt64() ;
		}
//				std::cout << njh::conToStr(optimJson.getMemberNames(), "\n") << std::endl;
//
//				std::cout << "optimJson[\"runParams_\"][\"klen_\"].asUInt64(): " << optimJson["runParams_"]["klen_"].asUInt64() << std::endl;
//
//				exit(1);
		if(!exists(utility.pairedR1Fnp_) && !exists(utility.singlesFnp_)){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", couldn't find " << utility.pairedR1Fnp_ << " or " << utility.singlesFnp_ << ", need to have at least one of them" << "\n";
			throw std::runtime_error{ss.str()};
		}
		if(!exists(utility.pairedR1Fnp_)){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "IDBAUD requires paired reads"<< "\n";
			throw std::runtime_error{ss.str()};
		}
		{
			concatenateFiles({utility.pairedR1Fnp_}, njh::files::make_path(regionOutputDir, "extracted_R1.fastq"));
			concatenateFiles({utility.pairedR2Fnp_}, njh::files::make_path(regionOutputDir, "extracted_R2.fastq"));
		}
		std::stringstream IDBAUDCmdStream;

		IDBAUDCmdStream << "cd " << regionOutputDir ;

		IDBAUDCmdStream << " && fq2fa --merge extracted_R1.fastq extracted_R2.fastq extracted.fasta";

		IDBAUDCmdStream << " && idba_ud ";

//
//
		//
		//idba_ud  -o ibda_testing

		IDBAUDCmdStream
						<< " -r extracted.fasta "
						<< " --num_threads " << utility.inputPars_.numThreads_
						<< " " << utility.inputPars_.extraProgramOptions_
						<< " -o " << IDBAUDOutDir
						<< " > IDBAUDRunLog_" << njh::getCurrentDate() << ".txt 2>&1";
		auto IDBAUDFullOutputDir = njh::files::make_path(regionOutputDir, IDBAUDOutDir);

		auto IDBAUDRunOutput = njh::sys::run({IDBAUDCmdStream.str()});

		BioCmdsUtils::checkRunOutThrow(IDBAUDRunOutput, __PRETTY_FUNCTION__);

		OutOptions IDBAUDRunOutputLogOpts(njh::files::make_path(IDBAUDFullOutputDir, "IDBAUDRunOutput.json"));
		OutputStream IDBAUDRunOutputLogOut(IDBAUDRunOutputLogOpts);
		IDBAUDRunOutputLogOut << njh::json::toJson(IDBAUDRunOutput) << std::endl;

		auto contigsFnp = njh::files::make_path(IDBAUDFullOutputDir, "contig.fa");


		auto contigsSeqIoOpts = SeqIOOptions::genFastaIn(contigsFnp);
		contigsSeqIoOpts.includeWhiteSpaceInName_ = false;
		contigsSeqIoOpts.lowerCaseBases_ = "upper";
		SeqInput contigsReader(contigsSeqIoOpts);
		auto contigsSeqs = contigsReader.readAllReads<seqInfo>();
		std::vector<std::shared_ptr<seqWithKmerInfo>> contigsKmerReads;
		contigsKmerReads.reserve(contigsSeqs.size());
		for (const auto & seq : contigsSeqs) {
			contigsKmerReads.emplace_back(std::make_shared<seqWithKmerInfo>(seq));
		}
		allSetKmers(contigsKmerReads, utility.inputPars_.reOrientingKmerLength_, true);

		RefSeqsWithKmers refSeqs(utility.refFnp_, utility.inputPars_.reOrientingKmerLength_);
		readVec::reorientSeqs(contigsKmerReads, refSeqs.refKmerReads_);
		//sort by sequence length;
		njh::sort(contigsKmerReads, [](const std::shared_ptr<seqWithKmerInfo> & seq1, const std::shared_ptr<seqWithKmerInfo> & seq2){
			return len(seq1->seqBase_) > len(seq2->seqBase_);
		});

		uint32_t defaultCoverage = 10;

		//sort by sequence length;
		njh::sort(contigsKmerReads, [](const std::shared_ptr<seqWithKmerInfo> & seq1, const std::shared_ptr<seqWithKmerInfo> & seq2){
			return len(seq1->seqBase_) > len(seq2->seqBase_);
		});

		OutOptions contigInfoOpts(njh::files::make_path(IDBAUDFullOutputDir, "contigs_outputInfo.tab.txt"));
		OutputStream contigInfoOut(contigInfoOpts);
		contigInfoOut << "name\tlength\tcoverage" << std::endl;

		for(const auto & contigsKmerRead : contigsKmerReads){
			//auto assembleInfo = DefaultAssembleNameInfo(contigsKmerRead->seqBase_.name_, true);
			contigInfoOut << contigsKmerRead->seqBase_.name_
										<< "\t" << len(contigsKmerRead->seqBase_)
										<< "\t" << defaultCoverage << std::endl;
		}
		auto reOrientedContigsFnp = njh::files::make_path(IDBAUDFullOutputDir, "reOriented_contigs.fasta");

		std::vector<std::shared_ptr<seqWithKmerInfo>> finalSeqs = trimToFinalSeqs(contigsKmerReads, refSeqs);
		std::unordered_map<std::string, uint32_t> finalSeqCounts;
		for(const auto & seq : finalSeqs){
			++finalSeqCounts[seq->seqBase_.name_];
		}

		double totalCoverage = defaultCoverage * finalSeqs.size();

		for(auto & seq : contigsKmerReads){
			//auto assembleInfo = DefaultAssembleNameInfo(seq->seqBase_.name_);
			MetaDataInName seqMeta;
			seqMeta.addMeta("length", len(seq->seqBase_));
			seqMeta.addMeta("estimatedPerBaseCoverage", defaultCoverage);
			seqMeta.addMeta("regionUID", utility.inputPars_.regionUid_);
			seqMeta.addMeta("sample", utility.inputPars_.sample_);
			seqMeta.resetMetaInName(seq->seqBase_.name_);
			seq->seqBase_.cnt_ = (defaultCoverage/totalCoverage) * (utility.totalCount());
			seq->seqBase_.name_ += njh::pasteAsStr("_t", seq->seqBase_.cnt_);
		}
		SeqOutput::write(contigsKmerReads, SeqIOOptions::genFastaOut(reOrientedContigsFnp));
		SeqOutput::write(contigsKmerReads, SeqIOOptions::genFastaOut(utility.outputFnp_));

		std::unordered_map<std::string, uint32_t> finalSeqCountsWritten;

		for(auto & seq : finalSeqs){
			//auto assembleInfo = DefaultAssembleNameInfo(seq->seqBase_.name_, true);
			MetaDataInName seqMeta;
			seqMeta.addMeta("trimmedLength", len(seq->seqBase_));
//					seqMeta.addMeta("estimatedPerBaseCoverage", kmerCoverage[seq->seqBase_.name_]);
			seqMeta.addMeta("estimatedPerBaseCoverage", defaultCoverage);

			seqMeta.addMeta("trimStatus", seq->seqBase_.on_);
			seqMeta.addMeta("regionUID", utility.inputPars_.regionUid_);
			seqMeta.addMeta("sample", utility.inputPars_.sample_);
			if(finalSeqCounts[seq->seqBase_.name_] > 1){
				seqMeta.addMeta("seqTrimmedCount", finalSeqCountsWritten[seq->seqBase_.name_]);
				++finalSeqCountsWritten[seq->seqBase_.name_];
			}
			//seq->seqBase_.cnt_ = (defaultCoverage/totalCoverage) * (readCounts.pairedReads_ + readCounts.unpaiedReads_ + readCounts.orphans_);
//					std::cout << "kmerCoverage[seq->seqBase_.name_]: " << kmerCoverage[seq->seqBase_.name_] << std::endl;
//					std::cout << "totalCoverage: " << totalCoverage << std::endl;
//					std::cout << "reads: " << readCounts.pairedReads_ + readCounts.unpaiedReads_ + readCounts.orphans_ << std::endl;
			seq->seqBase_.cnt_ = (defaultCoverage/totalCoverage) * (utility.totalCount());
			seqMeta.resetMetaInName(seq->seqBase_.name_);
			seq->seqBase_.name_ += njh::pasteAsStr("_t", seq->seqBase_.cnt_);
		}

		OutOptions trimmedContigInfoOpts(njh::files::make_path(IDBAUDFullOutputDir, "trimmed_reOriented_contigs_outputInfo.tab.txt"));
		OutputStream trimmedContigInfoOut(trimmedContigInfoOpts);
		trimmedContigInfoOut << "name\tlength\tcoverage" << std::endl;
		auto trimmedReOrientedContigsFnp = njh::files::make_path(IDBAUDFullOutputDir, "trimmed_reOriented_contigs.fasta");
		SeqOutput outputTrimmedWriter(SeqIOOptions::genFastaOut(trimmedReOrientedContigsFnp));
		auto trimmedReOrientedContigsFnp_belowCutOff = njh::files::make_path(IDBAUDFullOutputDir, "trimmed_reOriented_contigs_belowCutOff.fasta");
		SeqOutput belowCutOffOutputWriter(SeqIOOptions::genFastaOut(trimmedReOrientedContigsFnp_belowCutOff));

		uint32_t belowCutOff = 0;
		uint32_t aboveCutOff = 0;
		bool allPassTrim = true;
		for (const auto & contigsKmerRead : finalSeqs) {
			if (len(contigsKmerRead->seqBase_) < utility.inputPars_.minFinalLength_) {
				++belowCutOff;
				belowCutOffOutputWriter.openWrite(contigsKmerRead);
				contigsKmerRead->seqBase_.on_ = false;
			} else {
				MetaDataInName seqMeta(contigsKmerRead->seqBase_.name_);
				trimmedContigInfoOut << contigsKmerRead->seqBase_.name_
														 << "\t" << len(contigsKmerRead->seqBase_)
														 << "\t" << seqMeta.getMeta("estimatedPerBaseCoverage")
														 << std::endl;
				if(!contigsKmerRead->seqBase_.on_){
					allPassTrim = false;
				}else{
					++aboveCutOff;
				}
				outputAboveCutOffWriter.openWrite(contigsKmerRead);
				outputTrimmedWriter.openWrite(contigsKmerRead);
			}
		}
		if(allPassTrim){
			utility.regInfo_->infoCalled_ = true;
			utility.regInfo_->uniqHaps_ = aboveCutOff;
		}else{
			utility.regInfo_->infoCalled_ = false;
			utility.regInfo_->uniqHaps_ = 0;
		}
	} catch (std::exception & e) {
		exceptionMess = e.what();
		utility.regInfo_->infoCalled_ = false;
		utility.regInfo_->uniqHaps_ = 0;
	}

	outputAboveCutOffWriter.closeOut();
	OutputStream basicInfo(njh::files::make_path(utility.finalPassDir_, "basicInfoPerRegion.tab.txt"));

	basicInfo << "name\tsuccess\tuniqHaps\treadTotal\treadTotalUsed\ttotalPairedReads";
	basicInfo << "\tsample";
	basicInfo << "\n";
	basicInfo << utility.inputPars_.regionUid_;
	basicInfo << "\t" << njh::boolToStr(utility.regInfo_->infoCalled_)
						<< "\t" << utility.regInfo_->uniqHaps_
						<< "\t" << utility.regInfo_->totalReads_
						<< "\t" << utility.regInfo_->totalFinalReads_
						<< "\t" << utility.regInfo_->totalPairedReads_
						<< "\t" << utility.inputPars_.sample_ << std::endl;

	OutputStream exceptionsOut(njh::files::make_path(utility.finalPassDir_, "exceptionsMessages.tab.txt"));
	exceptionsOut << "regionUID\tmessage" << std::endl;
	exceptionsOut << utility.inputPars_.regionUid_ << "\t" << exceptionMess << std::endl;
	return 0;
}



int programWrappersAssembleOnPathWeaverRunner::runTrinityOnPathWeaverRegionsAndUnmapped(const njh::progutils::CmdArgs & inputCommands) {
	uint32_t TrinityMaxMemory = 10;

	bfs::path TrinityOutDir = "TrinityOut";
	seqSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processVerbose();
	OtherAssemblersUtility::InputPars inPars;
	inPars.programName_ = "Trinity";
	inPars.setPars(setUp);
	setUp.setOption(TrinityMaxMemory, "--TrinityMaxMemory", "Trinity Max Memory (in gigabytes)");
	setUp.setOption(TrinityOutDir,     "--TrinityOutDir",     "Trinity Out Directory name, will be relative to final pass directory");
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	njh::sys::requireExternalProgramThrow("Trinity");
	OtherAssemblersUtility utility(inPars);


	auto outputAboveCutOffSeqOpts = SeqIOOptions::genFastaOut(utility.outputAboveCutOffFnp_);
	SeqOutput outputAboveCutOffWriter(outputAboveCutOffSeqOpts);
	outputAboveCutOffWriter.openOut();

	std::string exceptionMess;

	auto regionOutputDir = njh::files::make_path(setUp.pars_.directoryName_, utility.inputPars_.regionUid_, utility.inputPars_.sample_);
	njh::files::makeDirP(njh::files::MkdirPar{regionOutputDir});

	try {
		if(!exists(utility.pairedR1Fnp_) && !exists(utility.singlesFnp_)){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", couldn't find " << utility.pairedR1Fnp_ << " or " << utility.singlesFnp_ << ", need to have at least one of them" << "\n";
			throw std::runtime_error{ss.str()};
		}
		std::stringstream TrinityCmdStream;
		TrinityCmdStream << "cd " << regionOutputDir << " && Trinity ";

		if(exists(utility.pairedR1Fnp_)){
			if(!exists(utility.pairedR2Fnp_)){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", found: " << utility.pairedR1Fnp_ << " but couldn't find it's mate file: " << utility.pairedR2Fnp_ << "\n";
				throw std::runtime_error{ss.str()};
			}else{
				TrinityCmdStream << " --left " << njh::files::normalize(utility.pairedR1Fnp_) << " --right " << njh::files::normalize(utility.pairedR2Fnp_) << " ";
			}
		}else if(exists(utility.singlesFnp_)){
			TrinityCmdStream << " --single  " << njh::files::normalize(utility.singlesFnp_);
		}

		TrinityCmdStream  << " --seqType fq  --CPU " << utility.inputPars_.numThreads_
											<< " --max_memory " << TrinityMaxMemory << "G"
											<< " --min_contig_length " << utility.inputPars_.minFinalLength_
											<< " " << utility.inputPars_.extraProgramOptions_
											<< " --output " << TrinityOutDir
											<< " > TrinityRunLog_" << njh::getCurrentDate() << ".txt 2>&1";
		auto TrinityFullOutputDir = njh::files::make_path(regionOutputDir, TrinityOutDir);

		auto TrinityRunOutput = njh::sys::run({TrinityCmdStream.str()});

		BioCmdsUtils::checkRunOutThrow(TrinityRunOutput, __PRETTY_FUNCTION__);

		OutOptions TrinityRunOutputLogOpts(njh::files::make_path(TrinityFullOutputDir, "TrinityRunOutput.json"));
		OutputStream TrinityRunOutputLogOut(TrinityRunOutputLogOpts);
		TrinityRunOutputLogOut << njh::json::toJson(TrinityRunOutput) << std::endl;

		auto contigsFnp = njh::files::make_path(TrinityFullOutputDir, "Trinity.fasta");

		auto contigsSeqIoOpts = SeqIOOptions::genFastaIn(contigsFnp);
		contigsSeqIoOpts.includeWhiteSpaceInName_ = false;
		contigsSeqIoOpts.lowerCaseBases_ = "upper";
		SeqInput contigsReader(contigsSeqIoOpts);
		auto contigsSeqs = contigsReader.readAllReads<seqInfo>();
		std::vector<std::shared_ptr<seqWithKmerInfo>> contigsKmerReads;
		contigsKmerReads.reserve(contigsSeqs.size());
		for (const auto & seq : contigsSeqs) {
			contigsKmerReads.emplace_back(std::make_shared<seqWithKmerInfo>(seq));
		}
		allSetKmers(contigsKmerReads, utility.inputPars_.reOrientingKmerLength_, true);

		RefSeqsWithKmers refSeqs(utility.refFnp_, utility.inputPars_.reOrientingKmerLength_);
		readVec::reorientSeqs(contigsKmerReads, refSeqs.refKmerReads_);
		//sort by sequence length;
		njh::sort(contigsKmerReads, [](const std::shared_ptr<seqWithKmerInfo> & seq1, const std::shared_ptr<seqWithKmerInfo> & seq2){
			return len(seq1->seqBase_) > len(seq2->seqBase_);
		});

		OutOptions contigInfoOpts(njh::files::make_path(TrinityFullOutputDir, "contigs_outputInfo.tab.txt"));
		OutputStream contigInfoOut(contigInfoOpts);
		contigInfoOut << "name\tlength" << std::endl;
		uint32_t defaultCoverage = 10;
		for(const auto & contigsKmerRead : contigsKmerReads){
			//auto assembleInfo = DefaultAssembleNameInfo(contigsKmerRead->seqBase_.name_, true);
			contigInfoOut << contigsKmerRead->seqBase_.name_
										<< "\t" << len(contigsKmerRead->seqBase_)
										<< std::endl;
			//<< "\t" << assembleInfo.coverage_ << std::endl;
		}
		auto reOrientedContigsFnp = njh::files::make_path(TrinityFullOutputDir, "reOriented_contigs.fasta");

		std::vector<std::shared_ptr<seqWithKmerInfo>> finalSeqs = trimToFinalSeqs(contigsKmerReads, refSeqs);
		std::unordered_map<std::string, uint32_t> finalSeqCounts;
		for(const auto & seq : finalSeqs){
			++finalSeqCounts[seq->seqBase_.name_];
		}

		double totalCoverage = defaultCoverage * finalSeqs.size();

		for(auto & seq : contigsKmerReads){
			//auto assembleInfo = DefaultAssembleNameInfo(seq->seqBase_.name_);
			MetaDataInName seqMeta;
			seqMeta.addMeta("length", len(seq->seqBase_));
			seqMeta.addMeta("estimatedPerBaseCoverage", defaultCoverage);
			seqMeta.addMeta("regionUID", utility.inputPars_.regionUid_);
			seqMeta.addMeta("sample", utility.inputPars_.sample_);
			seqMeta.resetMetaInName(seq->seqBase_.name_);
			seq->seqBase_.cnt_ = (defaultCoverage/totalCoverage) * (utility.totalCount());
			seq->seqBase_.name_ += njh::pasteAsStr("_t", seq->seqBase_.cnt_);
		}
		SeqOutput::write(contigsKmerReads, SeqIOOptions::genFastaOut(reOrientedContigsFnp));
		SeqOutput::write(contigsKmerReads, SeqIOOptions::genFastaOut(utility.outputFnp_));

		std::unordered_map<std::string, uint32_t> finalSeqCountsWritten;


		for(auto & seq : finalSeqs){
			//auto assembleInfo = DefaultAssembleNameInfo(seq->seqBase_.name_, true);
			MetaDataInName seqMeta;
			seqMeta.addMeta("trimmedLength", len(seq->seqBase_));
			//seqMeta.addMeta("estimatedPerBaseCoverage", assembleInfo.coverage_);
			seqMeta.addMeta("estimatedPerBaseCoverage", defaultCoverage);
			seqMeta.addMeta("trimStatus", seq->seqBase_.on_);
			seqMeta.addMeta("regionUID", utility.inputPars_.regionUid_);
			seqMeta.addMeta("sample", utility.inputPars_.sample_);
			if(finalSeqCounts[seq->seqBase_.name_] > 1){
				seqMeta.addMeta("seqTrimmedCount", finalSeqCountsWritten[seq->seqBase_.name_]);
				++finalSeqCountsWritten[seq->seqBase_.name_];
			}
			seqMeta.resetMetaInName(seq->seqBase_.name_);
			seq->seqBase_.cnt_ = (defaultCoverage/totalCoverage) * (utility.totalCount());
			//seq->seqBase_.cnt_ = (assembleInfo.coverage_/totalCoverage) * (readCounts.pairedReads_ + readCounts.unpaiedReads_ + readCounts.orphans_);
			seq->seqBase_.name_ += njh::pasteAsStr("_t", seq->seqBase_.cnt_);
		}

		OutOptions trimmedContigInfoOpts(njh::files::make_path(TrinityFullOutputDir, "trimmed_reOriented_contigs_outputInfo.tab.txt"));
		OutputStream trimmedContigInfoOut(trimmedContigInfoOpts);
		trimmedContigInfoOut << "name\tlength\tcoverage" << std::endl;
		auto trimmedReOrientedContigsFnp = njh::files::make_path(TrinityFullOutputDir, "trimmed_reOriented_contigs.fasta");
		SeqOutput outputTrimmedWriter(SeqIOOptions::genFastaOut(trimmedReOrientedContigsFnp));
		auto trimmedReOrientedContigsFnp_belowCutOff = njh::files::make_path(TrinityFullOutputDir, "trimmed_reOriented_contigs_belowCutOff.fasta");
		SeqOutput belowCutOffOutputWriter(SeqIOOptions::genFastaOut(trimmedReOrientedContigsFnp_belowCutOff));

		uint32_t belowCutOff = 0;
		uint32_t aboveCutOff = 0;
		bool allPassTrim = true;
		for (const auto & contigsKmerRead : finalSeqs) {
			if (len(contigsKmerRead->seqBase_) < utility.inputPars_.minFinalLength_) {
				++belowCutOff;
				belowCutOffOutputWriter.openWrite(contigsKmerRead);
				contigsKmerRead->seqBase_.on_ = false;
			} else {
				MetaDataInName seqMeta(contigsKmerRead->seqBase_.name_);
				trimmedContigInfoOut << contigsKmerRead->seqBase_.name_
														 << "\t" << len(contigsKmerRead->seqBase_)
														 << "\t" << seqMeta.getMeta("estimatedPerBaseCoverage")
														 << std::endl;
				if(!contigsKmerRead->seqBase_.on_){
					allPassTrim = false;
				}else{
					++aboveCutOff;
				}
				outputAboveCutOffWriter.openWrite(contigsKmerRead);
				outputTrimmedWriter.openWrite(contigsKmerRead);
			}
		}
		if(allPassTrim){
			utility.regInfo_->infoCalled_ = true;
			utility.regInfo_->uniqHaps_ = aboveCutOff;
		}else{
			utility.regInfo_->infoCalled_ = false;
			utility.regInfo_->uniqHaps_ = 0;
		}
	} catch (std::exception & e) {
		exceptionMess = e.what();
		utility.regInfo_->infoCalled_ = false;
		utility.regInfo_->uniqHaps_ = 0;
	}

	outputAboveCutOffWriter.closeOut();
	OutputStream basicInfo(njh::files::make_path(utility.finalPassDir_, "basicInfoPerRegion.tab.txt"));

	basicInfo << "name\tsuccess\tuniqHaps\treadTotal\treadTotalUsed\ttotalPairedReads";
	basicInfo << "\tsample";
	basicInfo << "\n";
	basicInfo << utility.inputPars_.regionUid_;
	basicInfo << "\t" << njh::boolToStr(utility.regInfo_->infoCalled_)
						<< "\t" << utility.regInfo_->uniqHaps_
						<< "\t" << utility.regInfo_->totalReads_
						<< "\t" << utility.regInfo_->totalFinalReads_
						<< "\t" << utility.regInfo_->totalPairedReads_
						<< "\t" << utility.inputPars_.sample_ << std::endl;

	OutputStream exceptionsOut(njh::files::make_path(utility.finalPassDir_, "exceptionsMessages.tab.txt"));
	exceptionsOut << "regionUID\tmessage" << std::endl;
	exceptionsOut << utility.inputPars_.regionUid_ << "\t" << exceptionMess << std::endl;

	return 0;
}




int programWrappersAssembleOnPathWeaverRunner::runMegahitOnPathWeaverRegionsAndUnmapped(const njh::progutils::CmdArgs & inputCommands) {

	bfs::path megahitOutDir = "megahitOut";
	seqSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processVerbose();
	OtherAssemblersUtility::InputPars inPars;
	inPars.programName_ = "megahit";
	inPars.setPars(setUp);
	setUp.setOption(megahitOutDir,     "--megahitOutDir",     "megahit Out Directory name, will be relative to final pass directory");


	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	njh::sys::requireExternalProgramThrow("megahit");
	OtherAssemblersUtility utility(inPars);

	auto outputAboveCutOffSeqOpts = SeqIOOptions::genFastaOut(utility.outputAboveCutOffFnp_);
	SeqOutput outputAboveCutOffWriter(outputAboveCutOffSeqOpts);
	outputAboveCutOffWriter.openOut();

	std::string exceptionMess;

	auto regionOutputDir = njh::files::make_path(setUp.pars_.directoryName_, utility.inputPars_.regionUid_, utility.inputPars_.sample_);
	njh::files::makeDirP(njh::files::MkdirPar{regionOutputDir});

	try {
		if(!exists(utility.pairedR1Fnp_) && !exists(utility.singlesFnp_)){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", couldn't find " << utility.pairedR1Fnp_ << " or " << utility.singlesFnp_ << ", need to have at least one of them" << "\n";
			throw std::runtime_error{ss.str()};
		}
		std::stringstream megahitCmdStream;
		megahitCmdStream << "cd " << regionOutputDir << " && megahit ";

		if(exists(utility.pairedR1Fnp_)){
			if(!exists(utility.pairedR2Fnp_)){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", found: " << utility.pairedR1Fnp_ << " but couldn't find it's mate file: " << utility.pairedR2Fnp_ << "\n";
				throw std::runtime_error{ss.str()};
			}else{
				megahitCmdStream << " -1 " << njh::files::normalize(utility.pairedR1Fnp_) << " -2 " << njh::files::normalize(utility.pairedR2Fnp_) << " ";
			}
		}else if(exists(utility.singlesFnp_)){
			megahitCmdStream << " -r  " << njh::files::normalize(utility.singlesFnp_);
		}
		megahitCmdStream  << " -t " << utility.inputPars_.numThreads_
											<< " " << utility.inputPars_.extraProgramOptions_
											<< " -o " << megahitOutDir
											<< " > megahitRunLog_" << njh::getCurrentDate() << ".txt 2>&1";
		auto megahitFullOutputDir = njh::files::make_path(regionOutputDir, megahitOutDir);

		auto megahitRunOutput = njh::sys::run({megahitCmdStream.str()});

		BioCmdsUtils::checkRunOutThrow(megahitRunOutput, __PRETTY_FUNCTION__);

		OutOptions megahitRunOutputLogOpts(njh::files::make_path(megahitFullOutputDir, "megahitRunOutput.json"));
		OutputStream megahitRunOutputLogOut(megahitRunOutputLogOpts);
		megahitRunOutputLogOut << njh::json::toJson(megahitRunOutput) << std::endl;

		auto contigsFnp = njh::files::make_path(megahitFullOutputDir, "final.contigs.fa");

		auto contigsSeqIoOpts = SeqIOOptions::genFastaIn(contigsFnp);
//				contigsSeqIoOpts.includeWhiteSpaceInName_ = false;
		contigsSeqIoOpts.lowerCaseBases_ = "upper";
		SeqInput contigsReader(contigsSeqIoOpts);
		auto contigsSeqs = contigsReader.readAllReads<seqInfo>();
		std::vector<std::shared_ptr<seqWithKmerInfo>> contigsKmerReads;
		contigsKmerReads.reserve(contigsSeqs.size());
		for (const auto & seq : contigsSeqs) {
			contigsKmerReads.emplace_back(std::make_shared<seqWithKmerInfo>(seq));
		}
		allSetKmers(contigsKmerReads, utility.inputPars_.reOrientingKmerLength_, true);
		RefSeqsWithKmers refSeqs(utility.refFnp_, utility.inputPars_.reOrientingKmerLength_);
		readVec::reorientSeqs(contigsKmerReads, refSeqs.refKmerReads_);
		//sort by sequence length;
		njh::sort(contigsKmerReads, [](const std::shared_ptr<seqWithKmerInfo> & seq1, const std::shared_ptr<seqWithKmerInfo> & seq2){
			return len(seq1->seqBase_) > len(seq2->seqBase_);
		});

		OutOptions contigInfoOpts(njh::files::make_path(megahitFullOutputDir, "contigs_outputInfo.tab.txt"));
		OutputStream contigInfoOut(contigInfoOpts);
		contigInfoOut << "name\tlength\tcoverage" << std::endl;

		for(const auto & contigsKmerRead : contigsKmerReads){
			auto assembleInfo = DefaultAssembleNameInfo(contigsKmerRead->seqBase_.name_, true);
			contigInfoOut << contigsKmerRead->seqBase_.name_
										<< "\t" << len(contigsKmerRead->seqBase_)
										<< "\t" << assembleInfo.coverage_ << std::endl;
		}
		auto reOrientedContigsFnp = njh::files::make_path(megahitFullOutputDir, "reOriented_contigs.fasta");
		std::vector<std::shared_ptr<seqWithKmerInfo>> finalSeqs = trimToFinalSeqs(contigsKmerReads, refSeqs);
		std::unordered_map<std::string, uint32_t> finalSeqCounts;
		for(const auto & seq : finalSeqs){
			++finalSeqCounts[seq->seqBase_.name_];
		}

		double totalCoverage = 0;
		for(auto & seq : finalSeqs){
			auto assembleInfo = DefaultAssembleNameInfo(seq->seqBase_.name_, true);
			totalCoverage += assembleInfo.coverage_;
		}
		for(auto & seq : contigsKmerReads){
			auto assembleInfo = DefaultAssembleNameInfo(seq->seqBase_.name_, true);
			MetaDataInName seqMeta;
			seqMeta.addMeta("length", len(seq->seqBase_));
			seqMeta.addMeta("estimatedPerBaseCoverage", assembleInfo.coverage_);
			seqMeta.addMeta("regionUID", utility.inputPars_.regionUid_);
			seqMeta.addMeta("sample", utility.inputPars_.sample_);
			seqMeta.resetMetaInName(seq->seqBase_.name_);
			seq->seqBase_.cnt_ = (assembleInfo.coverage_/totalCoverage) * (utility.totalCount());
			seq->seqBase_.name_ += njh::pasteAsStr("_t", seq->seqBase_.cnt_);
		}
		SeqOutput::write(contigsKmerReads, SeqIOOptions::genFastaOut(reOrientedContigsFnp));
		SeqOutput::write(contigsKmerReads, SeqIOOptions::genFastaOut(utility.outputFnp_));

		std::unordered_map<std::string, uint32_t> finalSeqCountsWritten;

		for(auto & seq : finalSeqs){
			auto assembleInfo = DefaultAssembleNameInfo(seq->seqBase_.name_, true);
			MetaDataInName seqMeta;
			seqMeta.addMeta("trimmedLength", len(seq->seqBase_));
			seqMeta.addMeta("estimatedPerBaseCoverage", assembleInfo.coverage_);
			seqMeta.addMeta("trimStatus", seq->seqBase_.on_);
			seqMeta.addMeta("regionUID", utility.inputPars_.regionUid_);
			seqMeta.addMeta("sample", utility.inputPars_.sample_);
			if(finalSeqCounts[seq->seqBase_.name_] > 1){
				seqMeta.addMeta("seqTrimmedCount", finalSeqCountsWritten[seq->seqBase_.name_]);
				++finalSeqCountsWritten[seq->seqBase_.name_];
			}
			seqMeta.resetMetaInName(seq->seqBase_.name_);
			seq->seqBase_.cnt_ = (assembleInfo.coverage_/totalCoverage) * (utility.totalCount());
			seq->seqBase_.name_ += njh::pasteAsStr("_t", seq->seqBase_.cnt_);
		}

		OutOptions trimmedContigInfoOpts(njh::files::make_path(megahitFullOutputDir, "trimmed_reOriented_contigs_outputInfo.tab.txt"));
		OutputStream trimmedContigInfoOut(trimmedContigInfoOpts);
		trimmedContigInfoOut << "name\tlength\tcoverage" << std::endl;
		auto trimmedReOrientedContigsFnp = njh::files::make_path(megahitFullOutputDir, "trimmed_reOriented_contigs.fasta");
		SeqOutput outputTrimmedWriter(SeqIOOptions::genFastaOut(trimmedReOrientedContigsFnp));
		auto trimmedReOrientedContigsFnp_belowCutOff = njh::files::make_path(megahitFullOutputDir, "trimmed_reOriented_contigs_belowCutOff.fasta");
		SeqOutput belowCutOffOutputWriter(SeqIOOptions::genFastaOut(trimmedReOrientedContigsFnp_belowCutOff));

		uint32_t belowCutOff = 0;
		uint32_t aboveCutOff = 0;
		bool allPassTrim = true;
		for (const auto & contigsKmerRead : finalSeqs) {
			if (len(contigsKmerRead->seqBase_) < utility.inputPars_.minFinalLength_) {
				++belowCutOff;
				belowCutOffOutputWriter.openWrite(contigsKmerRead);
				contigsKmerRead->seqBase_.on_ = false;
			} else {
				MetaDataInName seqMeta(contigsKmerRead->seqBase_.name_);
				trimmedContigInfoOut << contigsKmerRead->seqBase_.name_
														 << "\t" << len(contigsKmerRead->seqBase_)
														 << "\t" << seqMeta.getMeta("estimatedPerBaseCoverage")
														 << std::endl;
				if(!contigsKmerRead->seqBase_.on_){
					allPassTrim = false;
				}else{
					++aboveCutOff;
				}
				outputAboveCutOffWriter.openWrite(contigsKmerRead);
				outputTrimmedWriter.openWrite(contigsKmerRead);
			}
		}
		if(allPassTrim){
			utility.regInfo_->infoCalled_ = true;
			utility.regInfo_->uniqHaps_ = aboveCutOff;
		}else{
			utility.regInfo_->infoCalled_ = false;
			utility.regInfo_->uniqHaps_ = 0;
		}
	} catch (std::exception & e) {
		exceptionMess = e.what();
		utility.regInfo_->infoCalled_ = false;
		utility.regInfo_->uniqHaps_ = 0;
	}

	outputAboveCutOffWriter.closeOut();
	OutputStream basicInfo(njh::files::make_path(utility.finalPassDir_, "basicInfoPerRegion.tab.txt"));

	basicInfo << "name\tsuccess\tuniqHaps\treadTotal\treadTotalUsed\ttotalPairedReads";
	basicInfo << "\tsample";
	basicInfo << "\n";
	basicInfo << utility.inputPars_.regionUid_;
	basicInfo << "\t" << njh::boolToStr(utility.regInfo_->infoCalled_)
						<< "\t" << utility.regInfo_->uniqHaps_
						<< "\t" << utility.regInfo_->totalReads_
						<< "\t" << utility.regInfo_->totalFinalReads_
						<< "\t" << utility.regInfo_->totalPairedReads_
						<< "\t" << utility.inputPars_.sample_ << std::endl;

	OutputStream exceptionsOut(njh::files::make_path(utility.finalPassDir_, "exceptionsMessages.tab.txt"));
	exceptionsOut << "regionUID\tmessage" << std::endl;
	exceptionsOut << utility.inputPars_.regionUid_ << "\t" << exceptionMess << std::endl;
	return 0;
}



int programWrappersAssembleOnPathWeaverRunner::runSavageOnPathWeaverRegionsAndUnmapped(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path haploconductPath = "/home/hathawan/sourceCodes/savage/HaploConduct/haploconduct";
	bfs::path savageOutDir = "savageOut";
	seqSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processVerbose();
	OtherAssemblersUtility::InputPars inPars;
	inPars.programName_ = "savage";
	inPars.setPars(setUp);

	setUp.setOption(haploconductPath, "--haploconductPath", "haploconductPath", true);
	setUp.setOption(savageOutDir,     "--savageOutDir",     "savage Out Directory name, will be relative to final pass directory");


	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	//njh::sys::requireExternalProgramThrow("savage");
	OtherAssemblersUtility utility(inPars);
	auto outputAboveCutOffSeqOpts = SeqIOOptions::genFastaOut(utility.outputAboveCutOffFnp_);
	SeqOutput outputAboveCutOffWriter(outputAboveCutOffSeqOpts);
	outputAboveCutOffWriter.openOut();

	std::string exceptionMess;

	auto regionOutputDir = njh::files::make_path(setUp.pars_.directoryName_, utility.inputPars_.regionUid_, utility.inputPars_.sample_);
	njh::files::makeDirP(njh::files::MkdirPar{regionOutputDir});

	try {
		if(!exists(utility.pairedR1Fnp_) && !exists(utility.singlesFnp_)){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", couldn't find " << utility.pairedR1Fnp_ << " or " << utility.singlesFnp_ << ", need to have at least one of them" << "\n";
			throw std::runtime_error{ss.str()};
		}
		std::stringstream savageCmdStream;
		savageCmdStream << "cd " << regionOutputDir << " && " << haploconductPath << " savage ";
		if(bfs::exists(utility.pairedR1Fnp_)){
			concatenateFiles({utility.pairedR1Fnp_}, njh::files::make_path(regionOutputDir, utility.pairedR1Fnp_.filename().replace_extension("")));
			concatenateFiles({utility.pairedR2Fnp_}, njh::files::make_path(regionOutputDir, utility.pairedR2Fnp_.filename().replace_extension("")));
		}
		if(bfs::exists(utility.singlesFnp_)){
			concatenateFiles({utility.singlesFnp_}, njh::files::make_path(regionOutputDir, utility.singlesFnp_.filename().replace_extension("")));
		}
		if(exists(utility.pairedR1Fnp_)){
			if(!exists(utility.pairedR2Fnp_)){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", found: " << utility.pairedR1Fnp_ << " but couldn't find it's mate file: " << utility.pairedR2Fnp_ << "\n";
				throw std::runtime_error{ss.str()};
			}else{
				savageCmdStream << " -p1 " << utility.pairedR1Fnp_.filename().replace_extension("") << " -p2 " << utility.pairedR2Fnp_.filename().replace_extension("") << " ";
			}
		}
		if(exists(utility.singlesFnp_)){
			savageCmdStream << " -s  " << utility.singlesFnp_.filename().replace_extension("");
		}
		savageCmdStream  << " -t " << utility.inputPars_.numThreads_
										 << " --split 1 "
										 << " " << utility.inputPars_.extraProgramOptions_
										 << " -o " << savageOutDir
										 << " > savageRunLog_" << njh::getCurrentDate() << ".txt 2>&1";

		// ~/sourceCodes/savage/HaploConduct/haploconduct  -p1 extracted_R1.fastq -p2 extracted_R2.fastq -s extracted.fastq -t 10 --split 1


		auto savageFullOutputDir = njh::files::make_path(regionOutputDir, savageOutDir);

		auto savageRunOutput = njh::sys::run({savageCmdStream.str()});

		BioCmdsUtils::checkRunOutThrow(savageRunOutput, __PRETTY_FUNCTION__);

		OutOptions savageRunOutputLogOpts(njh::files::make_path(savageFullOutputDir, "savageRunOutput.json"));
		OutputStream savageRunOutputLogOut(savageRunOutputLogOpts);
		savageRunOutputLogOut << njh::json::toJson(savageRunOutput) << std::endl;

		auto contigsFnp = njh::files::make_path(savageFullOutputDir, "contigs_stage_c.fasta");

		auto contigsSeqIoOpts = SeqIOOptions::genFastaIn(contigsFnp);
//				contigsSeqIoOpts.includeWhiteSpaceInName_ = false;
		contigsSeqIoOpts.lowerCaseBases_ = "upper";
		SeqInput contigsReader(contigsSeqIoOpts);
		auto contigsSeqs = contigsReader.readAllReads<seqInfo>();
		std::vector<std::shared_ptr<seqWithKmerInfo>> contigsKmerReads;
		contigsKmerReads.reserve(contigsSeqs.size());
		for (const auto & seq : contigsSeqs) {
			contigsKmerReads.emplace_back(std::make_shared<seqWithKmerInfo>(seq));
		}
		allSetKmers(contigsKmerReads, utility.inputPars_.reOrientingKmerLength_, true);

		RefSeqsWithKmers refSeqs(utility.refFnp_, utility.inputPars_.reOrientingKmerLength_);
		readVec::reorientSeqs(contigsKmerReads, refSeqs.refKmerReads_);
		//sort by sequence length;
		njh::sort(contigsKmerReads, [](const std::shared_ptr<seqWithKmerInfo> & seq1, const std::shared_ptr<seqWithKmerInfo> & seq2){
			return len(seq1->seqBase_) > len(seq2->seqBase_);
		});

		OutOptions contigInfoOpts(njh::files::make_path(savageFullOutputDir, "contigs_outputInfo.tab.txt"));
		OutputStream contigInfoOut(contigInfoOpts);
		contigInfoOut << "name\tlength" << std::endl;

		for(const auto & contigsKmerRead : contigsKmerReads){
			//auto assembleInfo = DefaultAssembleNameInfo(contigsKmerRead->seqBase_.name_, true);
			contigInfoOut << contigsKmerRead->seqBase_.name_
										<< "\t" << len(contigsKmerRead->seqBase_) << std::endl;
		}
		auto reOrientedContigsFnp = njh::files::make_path(savageFullOutputDir, "reOriented_contigs.fasta");

		std::vector<std::shared_ptr<seqWithKmerInfo>> finalSeqs = trimToFinalSeqs(contigsKmerReads, refSeqs);
		std::unordered_map<std::string, uint32_t> finalSeqCounts;
		for(const auto & seq : finalSeqs){
			++finalSeqCounts[seq->seqBase_.name_];
		}

		auto totalCoverage = static_cast<double>(finalSeqs.size());

		for(auto & seq : contigsKmerReads){
			//auto assembleInfo = DefaultAssembleNameInfo(seq->seqBase_.name_);
			MetaDataInName seqMeta;
			seqMeta.addMeta("length", len(seq->seqBase_));
			seqMeta.addMeta("estimatedPerBaseCoverage", 10);
			seqMeta.addMeta("regionUID", utility.inputPars_.regionUid_);
			seqMeta.addMeta("sample", utility.inputPars_.sample_);
			seqMeta.resetMetaInName(seq->seqBase_.name_);
			seq->seqBase_.cnt_ = (1/totalCoverage) * (utility.totalCount());
			seq->seqBase_.name_ += njh::pasteAsStr("_t", seq->seqBase_.cnt_);
		}
		SeqOutput::write(contigsKmerReads, SeqIOOptions::genFastaOut(reOrientedContigsFnp));
		SeqOutput::write(contigsKmerReads, SeqIOOptions::genFastaOut(utility.outputFnp_));

		std::unordered_map<std::string, uint32_t> finalSeqCountsWritten;



		for(auto & seq : finalSeqs){
			//auto assembleInfo = DefaultAssembleNameInfo(seq->seqBase_.name_, true);
			MetaDataInName seqMeta;
			seqMeta.addMeta("trimmedLength", len(seq->seqBase_));
			seqMeta.addMeta("estimatedPerBaseCoverage", 10);
			seqMeta.addMeta("trimStatus", seq->seqBase_.on_);
			seqMeta.addMeta("regionUID", utility.inputPars_.regionUid_);
			seqMeta.addMeta("sample", utility.inputPars_.sample_);
			if(finalSeqCounts[seq->seqBase_.name_] > 1){
				seqMeta.addMeta("seqTrimmedCount", finalSeqCountsWritten[seq->seqBase_.name_]);
				++finalSeqCountsWritten[seq->seqBase_.name_];
			}
			seqMeta.resetMetaInName(seq->seqBase_.name_);
			seq->seqBase_.cnt_ = (1/totalCoverage) * (utility.totalCount());
			seq->seqBase_.name_ += njh::pasteAsStr("_t", seq->seqBase_.cnt_);
		}


		OutOptions trimmedContigInfoOpts(njh::files::make_path(savageFullOutputDir, "trimmed_reOriented_contigs_outputInfo.tab.txt"));
		OutputStream trimmedContigInfoOut(trimmedContigInfoOpts);
		trimmedContigInfoOut << "name\tlength\tcoverage" << std::endl;
		auto trimmedReOrientedContigsFnp = njh::files::make_path(savageFullOutputDir, "trimmed_reOriented_contigs.fasta");
		SeqOutput outputTrimmedWriter(SeqIOOptions::genFastaOut(trimmedReOrientedContigsFnp));
		auto trimmedReOrientedContigsFnp_belowCutOff = njh::files::make_path(savageFullOutputDir, "trimmed_reOriented_contigs_belowCutOff.fasta");
		SeqOutput belowCutOffOutputWriter(SeqIOOptions::genFastaOut(trimmedReOrientedContigsFnp_belowCutOff));

		uint32_t belowCutOff = 0;
		uint32_t aboveCutOff = 0;
		bool allPassTrim = true;
		for (const auto & contigsKmerRead : finalSeqs) {
			if (len(contigsKmerRead->seqBase_) < utility.inputPars_.minFinalLength_) {
				++belowCutOff;
				belowCutOffOutputWriter.openWrite(contigsKmerRead);
				contigsKmerRead->seqBase_.on_ = false;
			} else {
				MetaDataInName seqMeta(contigsKmerRead->seqBase_.name_);
				trimmedContigInfoOut << contigsKmerRead->seqBase_.name_
														 << "\t" << len(contigsKmerRead->seqBase_)
														 << "\t" << seqMeta.getMeta("estimatedPerBaseCoverage")
														 << std::endl;
				if(!contigsKmerRead->seqBase_.on_){
					allPassTrim = false;
				}else{
					++aboveCutOff;
				}
				outputAboveCutOffWriter.openWrite(contigsKmerRead);
				outputTrimmedWriter.openWrite(contigsKmerRead);
			}
		}
		if(allPassTrim){
			utility.regInfo_->infoCalled_ = true;
			utility.regInfo_->uniqHaps_ = aboveCutOff;
		}else{
			utility.regInfo_->infoCalled_ = false;
			utility.regInfo_->uniqHaps_ = 0;
		}
	} catch (std::exception & e) {
		exceptionMess = e.what();
		utility.regInfo_->infoCalled_ = false;
		utility.regInfo_->uniqHaps_ = 0;
	}

	outputAboveCutOffWriter.closeOut();
	OutputStream basicInfo(njh::files::make_path(utility.finalPassDir_, "basicInfoPerRegion.tab.txt"));

	basicInfo << "name\tsuccess\tuniqHaps\treadTotal\treadTotalUsed\ttotalPairedReads";
	basicInfo << "\tsample";
	basicInfo << "\n";
	basicInfo << utility.inputPars_.regionUid_;
	basicInfo << "\t" << njh::boolToStr(utility.regInfo_->infoCalled_)
						<< "\t" << utility.regInfo_->uniqHaps_
						<< "\t" << utility.regInfo_->totalReads_
						<< "\t" << utility.regInfo_->totalFinalReads_
						<< "\t" << utility.regInfo_->totalPairedReads_
						<< "\t" << utility.inputPars_.sample_ << std::endl;

	OutputStream exceptionsOut(njh::files::make_path(utility.finalPassDir_, "exceptionsMessages.tab.txt"));
	exceptionsOut << "regionUID\tmessage" << std::endl;
	exceptionsOut << utility.inputPars_.regionUid_ << "\t" << exceptionMess << std::endl;

	return 0;
}
int programWrappersAssembleOnPathWeaverRunner::runPolyteOnPathWeaverRegionsAndUnmapped(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path haploconductPath = "/home/hathawan/sourceCodes/savage/HaploConduct/haploconduct";
	uint32_t hardInsertSizeCutOff = 10000;
	uint32_t mapQualityCutOff = 20;

	seqSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processVerbose();
	OtherAssemblersUtility::InputPars inPars;
	inPars.programName_ = "polyte";
	inPars.setPars(setUp);

	setUp.setOption(haploconductPath, "--haploconductPath", "haploconductPath", true);


	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	//njh::sys::requireExternalProgramThrow("polyte");
	OtherAssemblersUtility utility(inPars);
	auto outputAboveCutOffSeqOpts = SeqIOOptions::genFastaOut(utility.outputAboveCutOffFnp_);
	SeqOutput outputAboveCutOffWriter(outputAboveCutOffSeqOpts);
	outputAboveCutOffWriter.openOut();

	bfs::path extractBam = njh::files::make_path(utility.inputPars_.pwOutputDir_, "originalExtractionFiles", "extracted.bam");

	uint32_t insertSize = 0;
	double insertSizeSD = 0;
	{
		BamTools::BamReader bReader;
		bReader.Open(extractBam.string());
		//checkBamOpenThrow(bReader, extractionBam);
		std::vector<uint32_t> currentInsertSizes;
		BamTools::BamAlignment aln;

		while (bReader.GetNextAlignmentCore(aln)) {
			//skip alignments that don't start in this region
			//this way if regions are close to each other it will avoid counting twice
			if (aln.IsMapped() &&
					aln.IsMateMapped() &&
					aln.IsPrimaryAlignment() &&
					aln.IsReverseStrand() != aln.IsMateReverseStrand()) {
				if (aln.RefID == aln.MateRefID) {
					if(std::abs(aln.InsertSize) > hardInsertSizeCutOff || aln.MapQuality < mapQualityCutOff){
						continue;
					}
					currentInsertSizes.emplace_back(std::abs(aln.InsertSize));
				}
			}
		}
		insertSize = vectorMedianRef(currentInsertSizes);
		insertSizeSD = vectorStandardDeviationPop(currentInsertSizes);
	}
	double possibleAvgCoverage = 0;
	{
		auto perBaseCovFnp = njh::files::make_path(utility.inputPars_.pwOutputDir_, "extractionStats/perBaseCoveragePerRegion.bed");
		table perBaseCovTab(perBaseCovFnp, "\t", true);
		std::vector<double> covs =vecStrToVecNum<double>(perBaseCovTab.getColumn("perBaseCoverage"));
		possibleAvgCoverage = vectorMean(covs);
	}
	std::string exceptionMess;

	auto regionOutputDir = njh::files::make_path(setUp.pars_.directoryName_, utility.inputPars_.regionUid_, utility.inputPars_.sample_);
	njh::files::makeDirP(njh::files::MkdirPar{regionOutputDir});

	try {
		if(!exists(utility.pairedR1Fnp_) && !exists(utility.singlesFnp_)){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", couldn't find " << utility.pairedR1Fnp_ << " or " << utility.singlesFnp_ << ", need to have at least one of them" << "\n";
			throw std::runtime_error{ss.str()};
		}
		std::stringstream polyteCmdStream;
		polyteCmdStream << "cd " << regionOutputDir << " && " << haploconductPath << " polyte ";
		if(bfs::exists(utility.pairedR1Fnp_)){
			concatenateFiles({utility.pairedR1Fnp_}, njh::files::make_path(regionOutputDir, utility.pairedR1Fnp_.filename().replace_extension("")));
			concatenateFiles({utility.pairedR2Fnp_}, njh::files::make_path(regionOutputDir, utility.pairedR2Fnp_.filename().replace_extension("")));
		}
		if(bfs::exists(utility.singlesFnp_)){
			concatenateFiles({utility.singlesFnp_}, njh::files::make_path(regionOutputDir, utility.singlesFnp_.filename().replace_extension("")));
		}
		if(exists(utility.pairedR1Fnp_)){
			if(!exists(utility.pairedR2Fnp_)){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", found: " << utility.pairedR1Fnp_ << " but couldn't find it's mate file: " << utility.pairedR2Fnp_ << "\n";
				throw std::runtime_error{ss.str()};
			}else{
				polyteCmdStream << " -p1 " << utility.pairedR1Fnp_.filename().replace_extension("") << " -p2 " << utility.pairedR2Fnp_.filename().replace_extension("") << " ";
			}
		}
		if(exists(utility.singlesFnp_)){
			polyteCmdStream << " -s  " << utility.singlesFnp_.filename().replace_extension("");
		}
		polyteCmdStream  << " -t " << utility.inputPars_.numThreads_
										 << " " << utility.inputPars_.extraProgramOptions_
						         << " --hap_cov " << possibleAvgCoverage << " --insert_size " << insertSize << "  --stddev " << insertSizeSD
										 << " > polyteRunLog_" << njh::getCurrentDate() << ".txt 2>&1";

		// ~/sourceCodes/polyte/HaploConduct/haploconduct  -p1 extracted_R1.fastq -p2 extracted_R2.fastq -s extracted.fastq -t 10 --split 1


		auto polyteFullOutputDir = njh::files::make_path(regionOutputDir);

		auto polyteRunOutput = njh::sys::run({polyteCmdStream.str()});

		BioCmdsUtils::checkRunOutThrow(polyteRunOutput, __PRETTY_FUNCTION__);

		OutOptions polyteRunOutputLogOpts(njh::files::make_path(polyteFullOutputDir, "polyteRunOutput.json"));
		OutputStream polyteRunOutputLogOut(polyteRunOutputLogOpts);
		polyteRunOutputLogOut << njh::json::toJson(polyteRunOutput) << std::endl;

		auto contigsFnp = njh::files::make_path(polyteFullOutputDir, "contigs.fasta");

		auto contigsSeqIoOpts = SeqIOOptions::genFastaIn(contigsFnp);
//				contigsSeqIoOpts.includeWhiteSpaceInName_ = false;
		contigsSeqIoOpts.lowerCaseBases_ = "upper";
		SeqInput contigsReader(contigsSeqIoOpts);
		auto contigsSeqs = contigsReader.readAllReads<seqInfo>();
		std::vector<std::shared_ptr<seqWithKmerInfo>> contigsKmerReads;
		contigsKmerReads.reserve(contigsSeqs.size());
		for (const auto & seq : contigsSeqs) {
			contigsKmerReads.emplace_back(std::make_shared<seqWithKmerInfo>(seq));
		}
		allSetKmers(contigsKmerReads, utility.inputPars_.reOrientingKmerLength_, true);

		RefSeqsWithKmers refSeqs(utility.refFnp_, utility.inputPars_.reOrientingKmerLength_);
		readVec::reorientSeqs(contigsKmerReads, refSeqs.refKmerReads_);
		//sort by sequence length;
		njh::sort(contigsKmerReads, [](const std::shared_ptr<seqWithKmerInfo> & seq1, const std::shared_ptr<seqWithKmerInfo> & seq2){
			return len(seq1->seqBase_) > len(seq2->seqBase_);
		});

		OutOptions contigInfoOpts(njh::files::make_path(polyteFullOutputDir, "contigs_outputInfo.tab.txt"));
		OutputStream contigInfoOut(contigInfoOpts);
		contigInfoOut << "name\tlength" << std::endl;

		for(const auto & contigsKmerRead : contigsKmerReads){
			//auto assembleInfo = DefaultAssembleNameInfo(contigsKmerRead->seqBase_.name_, true);
			contigInfoOut << contigsKmerRead->seqBase_.name_
										<< "\t" << len(contigsKmerRead->seqBase_) << std::endl;
		}
		auto reOrientedContigsFnp = njh::files::make_path(polyteFullOutputDir, "reOriented_contigs.fasta");

		std::vector<std::shared_ptr<seqWithKmerInfo>> finalSeqs = trimToFinalSeqs(contigsKmerReads, refSeqs);
		std::unordered_map<std::string, uint32_t> finalSeqCounts;
		for(const auto & seq : finalSeqs){
			++finalSeqCounts[seq->seqBase_.name_];
		}

		auto totalCoverage = static_cast<double>(finalSeqs.size());

		for(auto & seq : contigsKmerReads){
			//auto assembleInfo = DefaultAssembleNameInfo(seq->seqBase_.name_);
			MetaDataInName seqMeta;
			seqMeta.addMeta("length", len(seq->seqBase_));
			seqMeta.addMeta("estimatedPerBaseCoverage", 10);
			seqMeta.addMeta("regionUID", utility.inputPars_.regionUid_);
			seqMeta.addMeta("sample", utility.inputPars_.sample_);
			seqMeta.resetMetaInName(seq->seqBase_.name_);
			seq->seqBase_.cnt_ = (1/totalCoverage) * (utility.totalCount());
			seq->seqBase_.name_ += njh::pasteAsStr("_t", seq->seqBase_.cnt_);
		}
		SeqOutput::write(contigsKmerReads, SeqIOOptions::genFastaOut(reOrientedContigsFnp));
		SeqOutput::write(contigsKmerReads, SeqIOOptions::genFastaOut(utility.outputFnp_));

		std::unordered_map<std::string, uint32_t> finalSeqCountsWritten;



		for(auto & seq : finalSeqs){
			//auto assembleInfo = DefaultAssembleNameInfo(seq->seqBase_.name_, true);
			MetaDataInName seqMeta;
			seqMeta.addMeta("trimmedLength", len(seq->seqBase_));
			seqMeta.addMeta("estimatedPerBaseCoverage", 10);
			seqMeta.addMeta("trimStatus", seq->seqBase_.on_);
			seqMeta.addMeta("regionUID", utility.inputPars_.regionUid_);
			seqMeta.addMeta("sample", utility.inputPars_.sample_);
			if(finalSeqCounts[seq->seqBase_.name_] > 1){
				seqMeta.addMeta("seqTrimmedCount", finalSeqCountsWritten[seq->seqBase_.name_]);
				++finalSeqCountsWritten[seq->seqBase_.name_];
			}
			seqMeta.resetMetaInName(seq->seqBase_.name_);
			seq->seqBase_.cnt_ = (1/totalCoverage) * (utility.totalCount());
			seq->seqBase_.name_ += njh::pasteAsStr("_t", seq->seqBase_.cnt_);
		}


		OutOptions trimmedContigInfoOpts(njh::files::make_path(polyteFullOutputDir, "trimmed_reOriented_contigs_outputInfo.tab.txt"));
		OutputStream trimmedContigInfoOut(trimmedContigInfoOpts);
		trimmedContigInfoOut << "name\tlength\tcoverage" << std::endl;
		auto trimmedReOrientedContigsFnp = njh::files::make_path(polyteFullOutputDir, "trimmed_reOriented_contigs.fasta");
		SeqOutput outputTrimmedWriter(SeqIOOptions::genFastaOut(trimmedReOrientedContigsFnp));
		auto trimmedReOrientedContigsFnp_belowCutOff = njh::files::make_path(polyteFullOutputDir, "trimmed_reOriented_contigs_belowCutOff.fasta");
		SeqOutput belowCutOffOutputWriter(SeqIOOptions::genFastaOut(trimmedReOrientedContigsFnp_belowCutOff));

		uint32_t belowCutOff = 0;
		uint32_t aboveCutOff = 0;
		bool allPassTrim = true;
		for (const auto & contigsKmerRead : finalSeqs) {
			if (len(contigsKmerRead->seqBase_) < utility.inputPars_.minFinalLength_) {
				++belowCutOff;
				belowCutOffOutputWriter.openWrite(contigsKmerRead);
				contigsKmerRead->seqBase_.on_ = false;
			} else {
				MetaDataInName seqMeta(contigsKmerRead->seqBase_.name_);
				trimmedContigInfoOut << contigsKmerRead->seqBase_.name_
														 << "\t" << len(contigsKmerRead->seqBase_)
														 << "\t" << seqMeta.getMeta("estimatedPerBaseCoverage")
														 << std::endl;
				if(!contigsKmerRead->seqBase_.on_){
					allPassTrim = false;
				}else{
					++aboveCutOff;
				}
				outputAboveCutOffWriter.openWrite(contigsKmerRead);
				outputTrimmedWriter.openWrite(contigsKmerRead);
			}
		}
		if(allPassTrim){
			utility.regInfo_->infoCalled_ = true;
			utility.regInfo_->uniqHaps_ = aboveCutOff;
		}else{
			utility.regInfo_->infoCalled_ = false;
			utility.regInfo_->uniqHaps_ = 0;
		}
	} catch (std::exception & e) {
		exceptionMess = e.what();
		utility.regInfo_->infoCalled_ = false;
		utility.regInfo_->uniqHaps_ = 0;
	}

	outputAboveCutOffWriter.closeOut();
	OutputStream basicInfo(njh::files::make_path(utility.finalPassDir_, "basicInfoPerRegion.tab.txt"));

	basicInfo << "name\tsuccess\tuniqHaps\treadTotal\treadTotalUsed\ttotalPairedReads";
	basicInfo << "\tsample";
	basicInfo << "\n";
	basicInfo << utility.inputPars_.regionUid_;
	basicInfo << "\t" << njh::boolToStr(utility.regInfo_->infoCalled_)
						<< "\t" << utility.regInfo_->uniqHaps_
						<< "\t" << utility.regInfo_->totalReads_
						<< "\t" << utility.regInfo_->totalFinalReads_
						<< "\t" << utility.regInfo_->totalPairedReads_
						<< "\t" << utility.inputPars_.sample_ << std::endl;

	OutputStream exceptionsOut(njh::files::make_path(utility.finalPassDir_, "exceptionsMessages.tab.txt"));
	exceptionsOut << "regionUID\tmessage" << std::endl;
	exceptionsOut << utility.inputPars_.regionUid_ << "\t" << exceptionMess << std::endl;

	return 0;
}


int programWrappersAssembleOnPathWeaverRunner::runPRICEOnPathWeaverRegionsAndUnmapped(const njh::progutils::CmdArgs & inputCommands) {

	/// add options
	uint32_t hardInsertSizeCutOff = 10000;
	uint32_t mapQualityCutOff = 20;
	uint32_t numberOfCycles = 10;
	bfs::path priceOutputDir = "price_output";
	///
	bfs::path PRICEOutDir = "PRICEOut";
	seqSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processVerbose();
	OtherAssemblersUtility::InputPars inPars;
	inPars.programName_ = "PRICE";
	inPars.setPars(setUp);

	setUp.setOption(PRICEOutDir,     "--PRICEOutDir",     "PRICE Out Directory name, will be relative to final pass directory");
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	njh::sys::requireExternalProgramThrow("PriceTI");
	OtherAssemblersUtility utility(inPars);

	auto outputAboveCutOffSeqOpts = SeqIOOptions::genFastaOut(utility.outputAboveCutOffFnp_);
	SeqOutput outputAboveCutOffWriter(outputAboveCutOffSeqOpts);
	outputAboveCutOffWriter.openOut();

	std::string exceptionMess;

	auto regionOutputDir = njh::files::make_path(setUp.pars_.directoryName_, utility.inputPars_.regionUid_, utility.inputPars_.sample_);
	njh::files::makeDirP(njh::files::MkdirPar{regionOutputDir});


	// get insert size
	bfs::path extractBam = njh::files::make_path(utility.inputPars_.pwOutputDir_, "originalExtractionFiles", "extracted.bam");
	uint32_t insertSize = 0;
	{
		BamTools::BamReader bReader;
		bReader.Open(extractBam.string());
		//checkBamOpenThrow(bReader, extractionBam);
		std::vector<uint32_t> currentInsertSizes;
		BamTools::BamAlignment aln;

		while (bReader.GetNextAlignmentCore(aln)) {
			//skip alignments that don't start in this region
			//this way if regions are close to each other it will avoid counting twice
			if (aln.IsMapped() &&
					aln.IsMateMapped() &&
					aln.IsPrimaryAlignment() &&
					aln.IsReverseStrand() != aln.IsMateReverseStrand()) {
				if (aln.RefID == aln.MateRefID) {
					if(std::abs(aln.InsertSize) > hardInsertSizeCutOff || aln.MapQuality < mapQualityCutOff){
						continue;
					}
					currentInsertSizes.emplace_back(std::abs(aln.InsertSize));
				}
			}
		}
		insertSize = vectorMedianRef(currentInsertSizes);
	}
	try {
		if(!exists(utility.pairedR1Fnp_) && !exists(utility.singlesFnp_)){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", couldn't find " << utility.pairedR1Fnp_ << " or " << utility.singlesFnp_ << ", need to have at least one of them" << "\n";
			throw std::runtime_error{ss.str()};
		}
		std::stringstream PriceTICmdStream;
		PriceTICmdStream << "cd " << regionOutputDir << " && PriceTI ";
		if (!exists(utility.pairedR1Fnp_)){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "PriceTI only works on paired end reads"<< "\n";
			throw std::runtime_error{ss.str()};
		}
		if(utility.inputPars_.numThreads_ > 1){
			PriceTICmdStream << " -a " << utility.inputPars_.numThreads_ << " ";
		}
		if (!exists(utility.pairedR2Fnp_)) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", found: " << utility.pairedR1Fnp_
				 << " but couldn't find it's mate file: " << utility.pairedR2Fnp_ << "\n";
			throw std::runtime_error { ss.str() };
		} else {
			PriceTICmdStream
							<< " -fpp " << "extracted_R1.fastq" << " " << "extracted_R2.fastq" << " " << insertSize << " 95 ";
			PriceTICmdStream << " -nc " << numberOfCycles;
		}
		if(exists(utility.singlesFnp_)){
			PriceTICmdStream << " -icf " << "extracted.fastq" << " 2 1 1 ";
		}else{
			PriceTICmdStream << " -picf 400 " << "extracted_R1.fastq" << " 1 1 1 -picf 400 " << "extracted_R2.fastq" << " 1 1 1 ";
		}
		PriceTICmdStream << " -o " << priceOutputDir << "/price_out.fasta";

		if(bfs::exists(utility.singlesFnp_)){
			concatenateFiles({utility.singlesFnp_}, njh::files::make_path(regionOutputDir, "extracted.fastq"));
		}

		if(bfs::exists(utility.pairedR1Fnp_)){
			concatenateFiles({utility.pairedR1Fnp_}, njh::files::make_path(regionOutputDir, "extracted_R1.fastq"));
			concatenateFiles({utility.pairedR2Fnp_}, njh::files::make_path(regionOutputDir, "extracted_R2.fastq"));
		}


		auto PRICEFullOutputDir = njh::files::make_path(regionOutputDir, priceOutputDir);

		PriceTICmdStream
						<< " > " << priceOutputDir << "/priceOutputDirRunLog_" << njh::getCurrentDate() << ".txt 2>&1";
		njh::files::makeDir(njh::files::MkdirPar{PRICEFullOutputDir});

		auto PRICERunOutput = njh::sys::run({PriceTICmdStream.str()});

		BioCmdsUtils::checkRunOutThrow(PRICERunOutput, __PRETTY_FUNCTION__);

		OutOptions PRICERunOutputLogOpts(njh::files::make_path(PRICEFullOutputDir, "PRICERunOutput.json"));
		OutputStream PRICERunOutputLogOut(PRICERunOutputLogOpts);
		PRICERunOutputLogOut << njh::json::toJson(PRICERunOutput) << std::endl;

		auto contigsFnp = njh::files::make_path(PRICEFullOutputDir, njh::pasteAsStr("price_out.cycle", numberOfCycles, ".fasta"));

		auto contigsSeqIoOpts = SeqIOOptions::genFastaIn(contigsFnp);
//				contigsSeqIoOpts.includeWhiteSpaceInName_ = false;
		contigsSeqIoOpts.lowerCaseBases_ = "upper";
		SeqInput contigsReader(contigsSeqIoOpts);
		auto contigsSeqs = contigsReader.readAllReads<seqInfo>();
		std::vector<std::shared_ptr<seqWithKmerInfo>> contigsKmerReads;
		contigsKmerReads.reserve(contigsSeqs.size());
		for (const auto & seq : contigsSeqs) {
			contigsKmerReads.emplace_back(std::make_shared<seqWithKmerInfo>(seq));
		}
		allSetKmers(contigsKmerReads, utility.inputPars_.reOrientingKmerLength_, true);

		RefSeqsWithKmers refSeqs(utility.refFnp_, utility.inputPars_.reOrientingKmerLength_);
		readVec::reorientSeqs(contigsKmerReads, refSeqs.refKmerReads_);
		//sort by sequence length;
		njh::sort(contigsKmerReads, [](const std::shared_ptr<seqWithKmerInfo> & seq1, const std::shared_ptr<seqWithKmerInfo> & seq2){
			return len(seq1->seqBase_) > len(seq2->seqBase_);
		});

		OutOptions contigInfoOpts(njh::files::make_path(PRICEFullOutputDir, "contigs_outputInfo.tab.txt"));
		OutputStream contigInfoOut(contigInfoOpts);
		contigInfoOut << "name\tlength" << std::endl;
		uint32_t defaultCoverage = 10;

		for(const auto & contigsKmerRead : contigsKmerReads){
			//auto assembleInfo = DefaultAssembleNameInfo(contigsKmerRead->seqBase_.name_, true);
			contigInfoOut << contigsKmerRead->seqBase_.name_
										<< "\t" << len(contigsKmerRead->seqBase_)
										<< std::endl;
		}
		auto reOrientedContigsFnp = njh::files::make_path(PRICEFullOutputDir, "reOriented_contigs.fasta");


		std::vector<std::shared_ptr<seqWithKmerInfo>> finalSeqs = trimToFinalSeqs(contigsKmerReads, refSeqs);
		std::unordered_map<std::string, uint32_t> finalSeqCounts;
		for(const auto & seq : finalSeqs){
			++finalSeqCounts[seq->seqBase_.name_];
		}

		double totalCoverage = finalSeqs.size() * defaultCoverage;

		for(auto & seq : contigsKmerReads){
			//auto assembleInfo = DefaultAssembleNameInfo(seq->seqBase_.name_);
			MetaDataInName seqMeta;
			seqMeta.addMeta("length", len(seq->seqBase_));
			seqMeta.addMeta("estimatedPerBaseCoverage", 10);
			seqMeta.addMeta("regionUID", utility.inputPars_.regionUid_);
			seqMeta.addMeta("sample", utility.inputPars_.sample_);
			seqMeta.resetMetaInName(seq->seqBase_.name_);
			seq->seqBase_.cnt_ = (defaultCoverage/totalCoverage) * (utility.totalCount());
			seq->seqBase_.name_ += njh::pasteAsStr("_t", seq->seqBase_.cnt_);
		}
		SeqOutput::write(contigsKmerReads, SeqIOOptions::genFastaOut(reOrientedContigsFnp));
		SeqOutput::write(contigsKmerReads, SeqIOOptions::genFastaOut(utility.outputFnp_));

		std::unordered_map<std::string, uint32_t> finalSeqCountsWritten;


		for(auto & seq : finalSeqs){
			//auto assembleInfo = DefaultAssembleNameInfo(seq->seqBase_.name_, true);
			MetaDataInName seqMeta;
			seqMeta.addMeta("trimmedLength", len(seq->seqBase_));
			seqMeta.addMeta("estimatedPerBaseCoverage", defaultCoverage);
			seqMeta.addMeta("trimStatus", seq->seqBase_.on_);
			seqMeta.addMeta("regionUID", utility.inputPars_.regionUid_);
			seqMeta.addMeta("sample", utility.inputPars_.sample_);
			if(finalSeqCounts[seq->seqBase_.name_] > 1){
				seqMeta.addMeta("seqTrimmedCount", finalSeqCountsWritten[seq->seqBase_.name_]);
				++finalSeqCountsWritten[seq->seqBase_.name_];
			}
			seqMeta.resetMetaInName(seq->seqBase_.name_);
			seq->seqBase_.cnt_ = (defaultCoverage/totalCoverage) * (utility.totalCount());
			seq->seqBase_.name_ += njh::pasteAsStr("_t", seq->seqBase_.cnt_);
		}


		OutOptions trimmedContigInfoOpts(njh::files::make_path(PRICEFullOutputDir, "trimmed_reOriented_contigs_outputInfo.tab.txt"));
		OutputStream trimmedContigInfoOut(trimmedContigInfoOpts);
		trimmedContigInfoOut << "name\tlength\tcoverage" << std::endl;
		auto trimmedReOrientedContigsFnp = njh::files::make_path(PRICEFullOutputDir, "trimmed_reOriented_contigs.fasta");
		SeqOutput outputTrimmedWriter(SeqIOOptions::genFastaOut(trimmedReOrientedContigsFnp));
		auto trimmedReOrientedContigsFnp_belowCutOff = njh::files::make_path(PRICEFullOutputDir, "trimmed_reOriented_contigs_belowCutOff.fasta");
		SeqOutput belowCutOffOutputWriter(SeqIOOptions::genFastaOut(trimmedReOrientedContigsFnp_belowCutOff));

		uint32_t belowCutOff = 0;
		uint32_t aboveCutOff = 0;
		bool allPassTrim = true;
		for (const auto & contigsKmerRead : finalSeqs) {
			if (len(contigsKmerRead->seqBase_) < utility.inputPars_.minFinalLength_) {
				++belowCutOff;
				belowCutOffOutputWriter.openWrite(contigsKmerRead);
				contigsKmerRead->seqBase_.on_ = false;
			} else {
				MetaDataInName seqMeta(contigsKmerRead->seqBase_.name_);
				trimmedContigInfoOut << contigsKmerRead->seqBase_.name_
														 << "\t" << len(contigsKmerRead->seqBase_)
														 << "\t" << seqMeta.getMeta("estimatedPerBaseCoverage")
														 << std::endl;
				if(!contigsKmerRead->seqBase_.on_){
					allPassTrim = false;
				}else{
					++aboveCutOff;
				}
				outputAboveCutOffWriter.openWrite(contigsKmerRead);
				outputTrimmedWriter.openWrite(contigsKmerRead);
			}
		}
		if(allPassTrim){
			utility.regInfo_->infoCalled_ = true;
			utility.regInfo_->uniqHaps_ = aboveCutOff;
		}else{
			utility.regInfo_->infoCalled_ = false;
			utility.regInfo_->uniqHaps_ = 0;
		}
	} catch (std::exception & e) {
		exceptionMess = e.what();
		utility.regInfo_->infoCalled_ = false;
		utility.regInfo_->uniqHaps_ = 0;
	}

	outputAboveCutOffWriter.closeOut();
	OutputStream basicInfo(njh::files::make_path(utility.finalPassDir_, "basicInfoPerRegion.tab.txt"));

	basicInfo << "name\tsuccess\tuniqHaps\treadTotal\treadTotalUsed\ttotalPairedReads";
	basicInfo << "\tsample";
	basicInfo << "\n";
	basicInfo << utility.inputPars_.regionUid_;
	basicInfo << "\t" << njh::boolToStr(utility.regInfo_->infoCalled_)
						<< "\t" << utility.regInfo_->uniqHaps_
						<< "\t" << utility.regInfo_->totalReads_
						<< "\t" << utility.regInfo_->totalFinalReads_
						<< "\t" << utility.regInfo_->totalPairedReads_
						<< "\t" << utility.inputPars_.sample_ << std::endl;

	OutputStream exceptionsOut(njh::files::make_path(utility.finalPassDir_, "exceptionsMessages.tab.txt"));
	exceptionsOut << "regionUID\tmessage" << std::endl;
	exceptionsOut << utility.inputPars_.regionUid_ << "\t" << exceptionMess << std::endl;
	return 0;
}


int programWrappersAssembleOnPathWeaverRunner::runVelvetOptimizerAndMetaVelvetOnPathWeaverRegionsAndUnmapped(const njh::progutils::CmdArgs & inputCommands) {


	uint32_t velvetStartKmer = 31;
	uint32_t velvetEndKmer = 71;
	uint32_t velvetKmerStep = 10;
//	std::string optFuncKmer = "n50*Lcon/tbp+log(Lbp)";
//	std::string optFuncCov = "n50*Lcon/tbp+log(Lbp)";

	std::string optFuncKmer = "n50";
	std::string optFuncCov = "n50";


	VecStr optimizerFuncsAvail{"LNbp", "Lbp", "Lcon", "max", "n50", "ncon", "tbp"};

	bool overWriteDir = false;
	bool breakUpAmbigousContigs = true;
	double coverageCutOff = 2.0;

	std::string extraVelvetOptimiserOptions;
	bfs::path VelvetOptimiserOutDir = "VelvetOptimiserOutDir";

	seqSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processVerbose();
	OtherAssemblersUtility::InputPars inPars;
	inPars.programName_ = "velvet";
	inPars.setPars(setUp);
	OtherAssemblersUtility::InputPars inParsMetaVelvet = inPars;

	setUp.setOption(velvetStartKmer, "--velvetStartKmer", "velvet Start Kmer Size");
	setUp.setOption(velvetEndKmer, "--velvetEndKmer", "velvet End Kmer");
	setUp.setOption(velvetKmerStep, "--velvetKmerStep", "velvet Kmer Step");
	setUp.setOption(optFuncKmer, "--optFuncKmer", "optimizer Func Kmer");
	setUp.setOption(optFuncCov, "--optFuncCov", "optimizer Func Coverage");


	bool doNotBreakUpAmbigousContigs = false;
	setUp.setOption(doNotBreakUpAmbigousContigs, "--doBreakUpAmbigousContigs", "Do not Break Up Ambigous Contigs at Ns");
	breakUpAmbigousContigs = !doNotBreakUpAmbigousContigs;

	setUp.setOption(coverageCutOff, "--coverageCutOff", "Don't include these sequences in final output");
//	if(!njh::in(optFuncKmer, optimizerFuncsAvail)){
//		setUp.failed_ = true;
//		setUp.addWarning("Error for --optFuncKmer, value was set as " + optFuncKmer + " but doesn't match available options " + njh::conToStr(optimizerFuncsAvail, ", "));
//	}
//	if(!njh::in(optFuncCov, optimizerFuncsAvail)){
//		setUp.failed_ = true;
//		setUp.addWarning("Error for --optFuncCov, value was set as " + optFuncCov + " but doesn't match available options " + njh::conToStr(optimizerFuncsAvail, ", "));
//	}

	setUp.setOption(extraVelvetOptimiserOptions, "--extraVelvetOptimiserOptions", "Extra options to give to spades");
	setUp.setOption(VelvetOptimiserOutDir, "--VelvetOptimiserOutDir",
									"VelvetOptimiser.pl Out Directory name, will be relative to final pass directory");
	if ("VelvetOptimiserOutDir" == VelvetOptimiserOutDir) {
		if (!njh::in(optFuncKmer, optimizerFuncsAvail) || !njh::in(optFuncCov, optimizerFuncsAvail)) {
			VelvetOptimiserOutDir = VelvetOptimiserOutDir.string() + "_complex";
		} else {
			VelvetOptimiserOutDir = VelvetOptimiserOutDir.string() + "_kfunc" + optFuncKmer + "_covfunc" + optFuncCov;
		}
	}
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	OtherAssemblersUtility utility(inPars);

	bfs::path metaVelvetFinalDir = njh::replaceString(setUp.pars_.directoryName_, "Velvet", "MetaVelvet");
	njh::files::MkdirPar metaVelvetFinalDirPar(metaVelvetFinalDir, setUp.pars_.overWriteDir_);
	if (bfs::exists(metaVelvetFinalDir) && setUp.pars_.overWriteDir_) {
		njh::files::rmDirForce(metaVelvetFinalDir);
	}
	njh::files::makeDirP(metaVelvetFinalDirPar);
	inParsMetaVelvet.outputDir_ = metaVelvetFinalDir;
	OtherAssemblersUtility utilityMetaVelvet(inParsMetaVelvet);


	njh::sys::requireExternalProgramThrow("meta-velvetg");
	njh::sys::requireExternalProgramThrow("velveth");
	njh::sys::requireExternalProgramThrow("velvetg");
	njh::sys::requireExternalProgramThrow("VelvetOptimiser.pl");

	auto outputAboveCutOffSeqOptsvOpt = SeqIOOptions::genFastaOut(utility.outputAboveCutOffFnp_);
	SeqOutput outputAboveCutOffvOptWriter(outputAboveCutOffSeqOptsvOpt);
	outputAboveCutOffvOptWriter.openOut();

	auto outputAboveCutOffSeqOptsMetaVel = SeqIOOptions::genFastaOut(utilityMetaVelvet.outputAboveCutOffFnp_);
	SeqOutput outputAboveCutOffMetaVelWriter(outputAboveCutOffSeqOptsMetaVel);
	outputAboveCutOffMetaVelWriter.openOut();

	std::string exceptionMess;

	auto regionOutputDir = njh::files::make_path(setUp.pars_.directoryName_, utility.inputPars_.regionUid_,
																							 utility.inputPars_.sample_);
	njh::files::makeDirP(njh::files::MkdirPar{regionOutputDir});


	auto vOptFullOutputDir = njh::files::make_path(regionOutputDir, VelvetOptimiserOutDir);

	try {

		if (!exists(utility.pairedR1Fnp_) && !exists(utility.singlesFnp_)) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", couldn't find " << utility.pairedR1Fnp_ << " or " << utility.singlesFnp_
				 << ", need to have at least one of them" << "\n";
			throw std::runtime_error{ss.str()};
		}

		std::stringstream vOptCmdStream;
		vOptCmdStream << "cd " << regionOutputDir << " && VelvetOptimiser.pl -f '";
		if (exists(utility.pairedR1Fnp_)) {
			if (!exists(utility.pairedR2Fnp_)) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", found: " << utility.pairedR1Fnp_ << " but couldn't find it's mate file: "
					 << utility.pairedR2Fnp_ << "\n";
				throw std::runtime_error{ss.str()};
			} else {
				vOptCmdStream << " -fastq -shortPaired -separate " << njh::files::normalize(utility.pairedR1Fnp_) << " "
											<< njh::files::normalize(utility.pairedR2Fnp_) << " ";
			}
		}
		if (exists(utility.singlesFnp_)) {
			vOptCmdStream << " -fastq -short " << njh::files::normalize(utility.singlesFnp_) << " ";
		}
		vOptCmdStream << "'";
		//VelvetOptimiser.pl  -x 2 -f ' '  --d withRevComp_optimize_n50

		vOptCmdStream << " -t " << utility.inputPars_.numThreads_
									<< " -s " << velvetStartKmer
									<< " -e " << velvetEndKmer
									<< " -x " << velvetKmerStep
									<< " -optFuncKmer '" << optFuncKmer << "'"
									<< " -optFuncCov '" << optFuncCov << "'"
									<< " " << extraVelvetOptimiserOptions
									<< " --d " << VelvetOptimiserOutDir
									<< " -o '-exp_cov auto -cov_cutoff 300' ";
		std::string vOptCmdPreCovCutOff = vOptCmdStream.str();
		vOptCmdStream
						<< " -m " << coverageCutOff
						<< " > VelvetOptimzerRunLog_" << njh::getCurrentDate() << ".txt 2>&1";
		if (exists(vOptFullOutputDir)) {
			if (overWriteDir) {
				njh::files::rmDirForce(vOptFullOutputDir);
			} else {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << vOptFullOutputDir
					 << " already exists, use --overWriteDir to overwrite" << "\n";
				throw std::runtime_error{ss.str()};
			}
		}

		auto vOptRunOutput = njh::sys::run({vOptCmdStream.str()});
		auto firstRunOutputJson = njh::json::toJson(vOptRunOutput);
		bool failedFirstRun = false;
		if (!vOptRunOutput.success_) {
			failedFirstRun = true;
			if (exists(vOptFullOutputDir)) {
				bfs::rename(vOptFullOutputDir,
										njh::files::prependFileBasename(vOptFullOutputDir, "failed_" + njh::getCurrentDate() + "_"));
			}
			std::stringstream vOptCmdStreamAgain;
			vOptCmdStreamAgain << vOptCmdPreCovCutOff
												 << " > VelvetOptimzerRunLog_" << njh::getCurrentDate() << ".txt 2>&1";
			vOptRunOutput = njh::sys::run({vOptCmdStreamAgain.str()});
		}

		//failed a second time
		auto secondRunOutputJson = njh::json::toJson(vOptRunOutput);
		bool failedSecondRun = false;
		if (!vOptRunOutput.success_) {
			failedSecondRun = true;
			if (exists(vOptFullOutputDir)) {
				bfs::rename(vOptFullOutputDir,
										njh::files::prependFileBasename(vOptFullOutputDir,
																										"failed2_" + njh::getCurrentDate() + "_"));
			}
			std::stringstream vOptCmdStreamAgain;
			vOptCmdStreamAgain << vOptCmdPreCovCutOff
												 << " -m " << "-1" << " > VelvetOptimzerRunLog_"
												 << njh::getCurrentDate() << ".txt 2>&1";
			vOptRunOutput = njh::sys::run({vOptCmdStreamAgain.str()});
		}


		BioCmdsUtils::checkRunOutThrow(vOptRunOutput, __PRETTY_FUNCTION__);
		//setUp.startARunLog(vOptFullOutputDir.string());

		OutOptions spadesRunOutputLogOpts(njh::files::make_path(vOptFullOutputDir, "VelvetOptimiserRunOutput.json"));
		OutputStream spadesRunOutputLogOut(spadesRunOutputLogOpts);
		spadesRunOutputLogOut << njh::json::toJson(vOptRunOutput) << std::endl;

		if (failedFirstRun) {
			OutOptions failedVelvetRunOutputLogOpts(
							njh::files::make_path(vOptFullOutputDir, "failed_VelvetOptimiserRunOutput.json"));
			OutputStream failedVelvetRunOutputLogOut(failedVelvetRunOutputLogOpts);
			failedVelvetRunOutputLogOut << firstRunOutputJson << std::endl;
		}

		if (failedSecondRun) {
			OutOptions failedVelvetRunOutputLogOpts(
							njh::files::make_path(vOptFullOutputDir, "failed_VelvetOptimiserSecondTimeRunOutput.json"));
			OutputStream failedVelvetRunOutputLogOut(failedVelvetRunOutputLogOpts);
			failedVelvetRunOutputLogOut << secondRunOutputJson << std::endl;
		}

		{ //velvet optimizer run
			auto contigsFnp = njh::files::make_path(vOptFullOutputDir, "contigs.fasta");
			{
				//make contigs upper case
				auto originalContigsFnp = njh::files::make_path(vOptFullOutputDir, "contigs.fa");
				auto originalContigsOpts = SeqIOOptions::genFastaInOut(originalContigsFnp, contigsFnp);
				originalContigsOpts.lowerCaseBases_ = "upper";
				SeqInput originalContigsReader(originalContigsOpts);
				auto originalContigsSeqs = originalContigsReader.readAllReads<seqInfo>();
				SeqOutput::write(originalContigsSeqs, originalContigsOpts);
			}

			if (breakUpAmbigousContigs) {
				std::regex pat{"N+"};
				bool mark = true;
				//break up
				auto beforeBreakupContigsFnp = njh::files::make_path(vOptFullOutputDir, "beforeBreakUp_contigs.fasta");
				bfs::copy_file(contigsFnp, beforeBreakupContigsFnp);
				auto contigsOpts = SeqIOOptions::genFastaInOut(beforeBreakupContigsFnp, contigsFnp);
				contigsOpts.out_.overWriteFile_ = true;
				SeqIO contigsReader(contigsOpts);
				auto contigsSeqs = contigsReader.in_.readAllReads<seqInfo>();
				contigsReader.out_.openOut();
				for (const auto &seq: contigsSeqs) {
					auto trimmedSeqs = readVecTrimmer::breakUpSeqOnPat(seq, pat);
					for (auto &trimmedSeq: trimmedSeqs) {
						if (mark && 0 != trimmedSeq.start_ && len(seq) != trimmedSeq.end_) {
							trimmedSeq.seqBase_.name_.append(njh::pasteAsStr("-s", trimmedSeq.start_, "-e", trimmedSeq.end_));
						}
						contigsReader.write(trimmedSeq.seqBase_);
					}
				}
			}


//					auto contigsOpts = SeqIOOptions::genFastaIn(contigsFnp);
			auto contigsSeqIoOpts = SeqIOOptions::genFastaIn(contigsFnp);
//					contigsSeqIoOpts.includeWhiteSpaceInName_ = false;
			contigsSeqIoOpts.lowerCaseBases_ = "upper";
			SeqInput contigsReader(contigsSeqIoOpts);
			auto contigsSeqs = contigsReader.readAllReads<seqInfo>();

			std::vector<std::shared_ptr<seqWithKmerInfo>> contigsKmerReads;
			contigsKmerReads.reserve(contigsSeqs.size());
			for (const auto &seq: contigsSeqs) {
				contigsKmerReads.emplace_back(std::make_shared<seqWithKmerInfo>(seq));
			}

			allSetKmers(contigsKmerReads, utility.inputPars_.reOrientingKmerLength_, true);

			RefSeqsWithKmers refSeqs(utility.refFnp_, utility.inputPars_.reOrientingKmerLength_);
			readVec::reorientSeqs(contigsKmerReads, refSeqs.refKmerReads_);
			//sort by sequence length;
			njh::sort(contigsKmerReads,
								[](const std::shared_ptr<seqWithKmerInfo> &seq1, const std::shared_ptr<seqWithKmerInfo> &seq2) {
									return len(seq1->seqBase_) > len(seq2->seqBase_);
								});


			OutOptions contigInfoOpts(njh::files::make_path(vOptFullOutputDir, "contigs_outputInfo.tab.txt"));
			OutputStream contigInfoOut(contigInfoOpts);
			contigInfoOut << "name\tlength\tcoverage" << std::endl;

			for (const auto &contigsKmerRead: contigsKmerReads) {
				auto assembleInfo = DefaultAssembleNameInfo(contigsKmerRead->seqBase_.name_);
				contigInfoOut << contigsKmerRead->seqBase_.name_
											<< "\t" << len(contigsKmerRead->seqBase_)
											<< "\t" << assembleInfo.coverage_ << std::endl;
			}

			auto reOrientedContigsFnp = njh::files::make_path(vOptFullOutputDir, "reOriented_contigs.fasta");


			std::vector<std::shared_ptr<seqWithKmerInfo>> finalSeqs = trimToFinalSeqs(contigsKmerReads, refSeqs);
			std::unordered_map<std::string, uint32_t> finalSeqCounts;
			for (const auto &seq: finalSeqs) {
				++finalSeqCounts[seq->seqBase_.name_];
			}

			double totalCoverage = 0;
			for (auto &seq: finalSeqs) {
				auto assembleInfo = DefaultAssembleNameInfo(seq->seqBase_.name_);
				totalCoverage += assembleInfo.coverage_;
			}
			for (auto &seq: contigsKmerReads) {
				auto assembleInfo = DefaultAssembleNameInfo(seq->seqBase_.name_);
				MetaDataInName seqMeta;
				seqMeta.addMeta("length", len(seq->seqBase_));
				seqMeta.addMeta("estimatedPerBaseCoverage", 10);
				seqMeta.addMeta("regionUID", utility.inputPars_.regionUid_);
				seqMeta.addMeta("sample", utility.inputPars_.sample_);
				seqMeta.resetMetaInName(seq->seqBase_.name_);
				seq->seqBase_.cnt_ = (assembleInfo.coverage_ / totalCoverage) * (utility.totalCount());
				seq->seqBase_.name_ += njh::pasteAsStr("_t", seq->seqBase_.cnt_);
			}
			SeqOutput::write(contigsKmerReads, SeqIOOptions::genFastaOut(reOrientedContigsFnp));
			SeqOutput::write(contigsKmerReads, SeqIOOptions::genFastaOut(utility.outputFnp_));


			std::unordered_map<std::string, uint32_t> finalSeqCountsWritten;

			for (auto &seq: finalSeqs) {
				auto assembleInfo = DefaultAssembleNameInfo(seq->seqBase_.name_);
				MetaDataInName seqMeta;
				seqMeta.addMeta("trimmedLength", len(seq->seqBase_));
				seqMeta.addMeta("estimatedPerBaseCoverage", assembleInfo.coverage_);
				seqMeta.addMeta("trimStatus", seq->seqBase_.on_);
				seqMeta.addMeta("regionUID", utility.inputPars_.regionUid_);
				seqMeta.addMeta("sample", utility.inputPars_.sample_);
				if (finalSeqCounts[seq->seqBase_.name_] > 1) {
					seqMeta.addMeta("seqTrimmedCount", finalSeqCountsWritten[seq->seqBase_.name_]);
					++finalSeqCountsWritten[seq->seqBase_.name_];
				}
				seqMeta.resetMetaInName(seq->seqBase_.name_);
				seq->seqBase_.cnt_ = (assembleInfo.coverage_ / totalCoverage) * (utility.totalCount());
				seq->seqBase_.name_ += njh::pasteAsStr("_t", seq->seqBase_.cnt_);
			}

			OutOptions trimmedContigInfoOpts(
							njh::files::make_path(vOptFullOutputDir, "trimmed_reOriented_contigs_outputInfo.tab.txt"));
			OutputStream trimmedContigInfoOut(trimmedContigInfoOpts);
			trimmedContigInfoOut << "name\tlength\tcoverage" << std::endl;
			auto trimmedReOrientedContigsFnp = njh::files::make_path(vOptFullOutputDir, "trimmed_reOriented_contigs.fasta");
			SeqOutput outputTrimmedWriter(SeqIOOptions::genFastaOut(trimmedReOrientedContigsFnp));
			auto trimmedReOrientedContigsFnp_belowCutOff = njh::files::make_path(vOptFullOutputDir,
																																					 "trimmed_reOriented_contigs_belowCutOff.fasta");
			SeqOutput belowCutOffOutputWriter(SeqIOOptions::genFastaOut(trimmedReOrientedContigsFnp_belowCutOff));

			uint32_t belowCutOff = 0;
			uint32_t aboveCutOff = 0;
			bool allPassTrim = true;
			for (const auto &contigsKmerRead: finalSeqs) {
				if (len(contigsKmerRead->seqBase_) < utility.inputPars_.minFinalLength_) {
					++belowCutOff;
					belowCutOffOutputWriter.openWrite(contigsKmerRead);
					contigsKmerRead->seqBase_.on_ = false;
				} else {
					MetaDataInName seqMeta(contigsKmerRead->seqBase_.name_);
					trimmedContigInfoOut << contigsKmerRead->seqBase_.name_
															 << "\t" << len(contigsKmerRead->seqBase_)
															 << "\t" << seqMeta.getMeta("estimatedPerBaseCoverage")
															 << std::endl;
					if (!contigsKmerRead->seqBase_.on_) {
						allPassTrim = false;
					} else {
						++aboveCutOff;
					}
					outputAboveCutOffvOptWriter.openWrite(contigsKmerRead);
					outputTrimmedWriter.openWrite(contigsKmerRead);
				}
			}
			if (allPassTrim) {
				utility.regInfo_->infoCalled_ = true;
				utility.regInfo_->uniqHaps_ = aboveCutOff;
			} else {
				utility.regInfo_->infoCalled_ = false;
				utility.regInfo_->uniqHaps_ = 0;
			}
		}
	} catch (std::exception &e) {
		exceptionMess = e.what();
		utility.regInfo_->infoCalled_ = false;
		utility.regInfo_->uniqHaps_ = 0;
	}

	{
		outputAboveCutOffvOptWriter.closeOut();
		OutputStream basicInfo(njh::files::make_path(utility.finalPassDir_, "basicInfoPerRegion.tab.txt"));

		basicInfo << "name\tsuccess\tuniqHaps\treadTotal\treadTotalUsed\ttotalPairedReads";
		basicInfo << "\tsample";
		basicInfo << "\n";
		basicInfo << utility.inputPars_.regionUid_;
		basicInfo << "\t" << njh::boolToStr(utility.regInfo_->infoCalled_)
							<< "\t" << utility.regInfo_->uniqHaps_
							<< "\t" << utility.regInfo_->totalReads_
							<< "\t" << utility.regInfo_->totalFinalReads_
							<< "\t" << utility.regInfo_->totalPairedReads_
							<< "\t" << utility.inputPars_.sample_ << std::endl;

		OutputStream exceptionsOut(njh::files::make_path(utility.finalPassDir_, "exceptionsMessages.tab.txt"));
		exceptionsOut << "regionUID\tmessage" << std::endl;
		exceptionsOut << utility.inputPars_.regionUid_ << "\t" << exceptionMess << std::endl;
	}

	try {



		// now run meta-velvetg
		bfs::path MetaVelvetOutDir = njh::files::make_path(vOptFullOutputDir, "MetaVelvetOutDir");
		njh::files::makeDir(njh::files::MkdirPar{MetaVelvetOutDir});
		auto logfiles = njh::files::gatherFiles(vOptFullOutputDir, "_Logfile.txt");
		if(logfiles.empty() || 1 != logfiles.size()){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "couldn't determine log file in " << vOptFullOutputDir << "\n";
			if(logfiles.size() > 1){
				ss << "Multiple log files found: " << njh::conToStr(logfiles, ", ") << "\n";
			}
			throw std::runtime_error{ss.str()};
		}
		InputStream logFileIn(InOptions(logfiles.front()));
		std::string line;
		bool encounteredFinalInfo = false;
		std::string velvethArgs;
		std::string velvetgArgs;
		while(njh::files::crossPlatGetline(logFileIn, line)){
			if(njh::beginsWith(line, "Final optimised assembly details")){
				encounteredFinalInfo = true;
			}
			if(encounteredFinalInfo && std::string::npos != line.find(':')){
				auto toks = tokenizeString(line, ":");
				if(toks.size() == 2 && "Velveth parameter string" == toks[0]){
					njh::trim(toks[1]);
					velvethArgs = toks[1];
				}
				if(toks.size() == 2 && "Velvetg parameter string" == toks[0]){
					njh::trim(toks[1]);
					velvetgArgs = njh::replaceString(toks[1], "-clean yes", "");
				}
			}
		}

		if(velvethArgs.empty()){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "couldn't determine velveth args from " << logfiles.front()<< "\n";
			throw std::runtime_error{ss.str()};
		}

		if(velvetgArgs.empty()){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "couldn't determine velvetg args from " << logfiles.front()<< "\n";
			throw std::runtime_error{ss.str()};
		}
		//	std::string velvetOutputDirStr = velvetgArgs;
		//	trimAtFirstWhitespace(velvetOutputDirStr);
		std::string velvetOutputDirStr = "optimized_velvet_run";
		bfs::path velvetOutputDir = njh::files::make_path(VelvetOptimiserOutDir, "MetaVelvetOutDir", velvetOutputDirStr);

		std::stringstream velvethCmd;
		velvethCmd << "cd " << vOptFullOutputDir << "/../ && velveth " << velvetOutputDir << " " << velvethArgs.substr(velvethArgs.find(' '))
						<< " > VelvetOptimzer_velveth_RunLog_" << njh::getCurrentDate() << ".txt 2>&1 ";

		std::stringstream velvetgCmd;
		velvetgCmd << "cd " << vOptFullOutputDir << "/../ && velvetg " << velvetOutputDir << " " << velvetgArgs.substr(velvetgArgs.find(' '))
						<< " > VelvetOptimzer_velvetg_RunLog_" << njh::getCurrentDate() << ".txt 2>&1 ";

		std::stringstream metavelvetgCmd;
		metavelvetgCmd << "cd " << vOptFullOutputDir << "/../ && meta-velvetg " << velvetOutputDir
						<< " > VelvetOptimzer_meta-velvetg_RunLog_" << njh::getCurrentDate() << ".txt 2>&1 ";

		auto velvethCmdOutut = njh::sys::run(VecStr{velvethCmd.str()});
		BioCmdsUtils::checkRunOutThrow(velvethCmdOutut, __PRETTY_FUNCTION__);

		auto velvetgCmdOutut = njh::sys::run(VecStr{velvetgCmd.str()});
		BioCmdsUtils::checkRunOutThrow(velvetgCmdOutut, __PRETTY_FUNCTION__);

		auto metavelvetgCmdOutput = njh::sys::run(VecStr{metavelvetgCmd.str()});
		BioCmdsUtils::checkRunOutThrow(metavelvetgCmdOutput, __PRETTY_FUNCTION__);

		Json::Value runOutputs;
		runOutputs["velvethCmd"] = velvethCmdOutut.toJson();
		runOutputs["velvetgCmd"] = velvetgCmdOutut.toJson();
		runOutputs["metavelvetgCmdOutput"] = metavelvetgCmdOutput.toJson();
		{
			OutputStream runOutputsLogOut(njh::files::make_path(vOptFullOutputDir, "MetaVelvetOutDir", velvetOutputDirStr, "runOutputsLog.json"));
			runOutputsLogOut << runOutputs << std::endl;
		}



		{
			auto metaVelvetDir = njh::files::make_path(vOptFullOutputDir, "MetaVelvetOutDir", velvetOutputDirStr);
			bfs::path metaVelvetContigsfile = njh::files::make_path(metaVelvetDir, "meta-velvetg.contigs.fasta");
			{
				//make contigs upper case
				auto originalContigsFnp = njh::files::make_path(metaVelvetDir, "meta-velvetg.contigs.fa");
				auto originalContigsOpts = SeqIOOptions::genFastaInOut(originalContigsFnp, metaVelvetContigsfile);
				originalContigsOpts.lowerCaseBases_ = "upper";
				SeqInput originalContigsReader(originalContigsOpts);
				auto originalContigsSeqs = originalContigsReader.readAllReads<seqInfo>();
				SeqOutput::write(originalContigsSeqs, originalContigsOpts);
			}

			if (breakUpAmbigousContigs) {
				std::regex pat{"N+"};
				bool mark = true;
				//break up
				auto beforeBreakupContigsFnp = njh::files::make_path(metaVelvetDir, "beforeBreakUp_meta-velvetg.contigs.fasta");
				bfs::copy_file(metaVelvetContigsfile, beforeBreakupContigsFnp);
				auto contigsOpts = SeqIOOptions::genFastaInOut(beforeBreakupContigsFnp, metaVelvetContigsfile);
				contigsOpts.out_.overWriteFile_ = true;
				SeqIO contigsReader(contigsOpts);
				auto contigsSeqs = contigsReader.in_.readAllReads<seqInfo>();
				contigsReader.out_.openOut();
				for (const auto & seq : contigsSeqs) {
					auto trimmedSeqs = readVecTrimmer::breakUpSeqOnPat(seq, pat);
					for(auto & trimmedSeq : trimmedSeqs){
						if(mark && 0 != trimmedSeq.start_ && len(seq) != trimmedSeq.end_){
							trimmedSeq.seqBase_.name_.append(njh::pasteAsStr("-s", trimmedSeq.start_, "-e", trimmedSeq.end_));
						}
						contigsReader.write(trimmedSeq.seqBase_);
					}
				}
			}


//					auto contigsOpts = SeqIOOptions::genFastaIn(metaVelvetContigsfile);
			auto contigsSeqIoOpts = SeqIOOptions::genFastaIn(metaVelvetContigsfile);
//					contigsSeqIoOpts.includeWhiteSpaceInName_ = false;
			contigsSeqIoOpts.lowerCaseBases_ = "upper";
			SeqInput contigsReader(contigsSeqIoOpts);
			auto contigsSeqs = contigsReader.readAllReads<seqInfo>();

			std::vector<std::shared_ptr<seqWithKmerInfo>> contigsKmerReads;
			contigsKmerReads.reserve(contigsSeqs.size());
			for (const auto & seq : contigsSeqs) {
				contigsKmerReads.emplace_back(std::make_shared<seqWithKmerInfo>(seq));
			}

			allSetKmers(contigsKmerReads, utilityMetaVelvet.inputPars_.reOrientingKmerLength_, true);

			RefSeqsWithKmers refSeqs(utilityMetaVelvet.refFnp_, utilityMetaVelvet.inputPars_.reOrientingKmerLength_);
			readVec::reorientSeqs(contigsKmerReads, refSeqs.refKmerReads_);
			//sort by sequence length;
			njh::sort(contigsKmerReads, [](const std::shared_ptr<seqWithKmerInfo> & seq1, const std::shared_ptr<seqWithKmerInfo> & seq2){
				return len(seq1->seqBase_) > len(seq2->seqBase_);
			});


			OutOptions contigInfoOpts(njh::files::make_path(metaVelvetDir, "meta-velvetg.contigs_outputInfo.tab.txt"));
			OutputStream contigInfoOut(contigInfoOpts);
			contigInfoOut << "name\tlength\tcoverage" << std::endl;

			for(const auto & contigsKmerRead : contigsKmerReads){
				auto assembleInfo = DefaultAssembleNameInfo(contigsKmerRead->seqBase_.name_);
				contigInfoOut << contigsKmerRead->seqBase_.name_
											<< "\t" << len(contigsKmerRead->seqBase_)
											<< "\t" << assembleInfo.coverage_ << std::endl;
			}

			auto reOrientedContigsFnp = njh::files::make_path(metaVelvetDir, "reOriented_meta-velvetg.contigs.fasta");

			std::vector<std::shared_ptr<seqWithKmerInfo>> finalSeqs = trimToFinalSeqs(contigsKmerReads, refSeqs);
			std::unordered_map<std::string, uint32_t> finalSeqCounts;
			for (const auto &seq: finalSeqs) {
				++finalSeqCounts[seq->seqBase_.name_];
			}

			double totalCoverage = 0;
			for(auto & seq : finalSeqs){
				auto assembleInfo = DefaultAssembleNameInfo(seq->seqBase_.name_);
				totalCoverage += assembleInfo.coverage_;
			}
			for (auto &seq: contigsKmerReads) {
				auto assembleInfo = DefaultAssembleNameInfo(seq->seqBase_.name_);
				MetaDataInName seqMeta;
				seqMeta.addMeta("length", len(seq->seqBase_));
				seqMeta.addMeta("estimatedPerBaseCoverage", 10);
				seqMeta.addMeta("regionUID", utilityMetaVelvet.inputPars_.regionUid_);
				seqMeta.addMeta("sample", utilityMetaVelvet.inputPars_.sample_);
				seqMeta.resetMetaInName(seq->seqBase_.name_);
				seq->seqBase_.cnt_ = (assembleInfo.coverage_ / totalCoverage) * (utilityMetaVelvet.totalCount());
				seq->seqBase_.name_ += njh::pasteAsStr("_t", seq->seqBase_.cnt_);
			}
			SeqOutput::write(contigsKmerReads, SeqIOOptions::genFastaOut(reOrientedContigsFnp));
			SeqOutput::write(contigsKmerReads, SeqIOOptions::genFastaOut(utilityMetaVelvet.outputFnp_));

			std::unordered_map<std::string, uint32_t> finalSeqCountsWritten;

			for(auto & seq : finalSeqs){
				auto assembleInfo = DefaultAssembleNameInfo(seq->seqBase_.name_);
				MetaDataInName seqMeta;
				seqMeta.addMeta("trimmedLength", len(seq->seqBase_));
				seqMeta.addMeta("estimatedPerBaseCoverage", assembleInfo.coverage_);
				seqMeta.addMeta("trimStatus", seq->seqBase_.on_);
				seqMeta.addMeta("regionUID", utilityMetaVelvet.inputPars_.regionUid_);
				seqMeta.addMeta("sample", utilityMetaVelvet.inputPars_.sample_);
				if(finalSeqCounts[seq->seqBase_.name_] > 1){
					seqMeta.addMeta("seqTrimmedCount", finalSeqCountsWritten[seq->seqBase_.name_]);
					++finalSeqCountsWritten[seq->seqBase_.name_];
				}
				seqMeta.resetMetaInName(seq->seqBase_.name_);
				seq->seqBase_.cnt_ = (assembleInfo.coverage_/totalCoverage) * (utilityMetaVelvet.totalCount());
				seq->seqBase_.name_ += njh::pasteAsStr("_t", seq->seqBase_.cnt_);
			}

			OutOptions trimmedContigInfoOpts(
							njh::files::make_path(metaVelvetDir, "trimmed_reOriented_contigs_outputInfo.tab.txt"));
			OutputStream trimmedContigInfoOut(trimmedContigInfoOpts);
			trimmedContigInfoOut << "name\tlength\tcoverage" << std::endl;
			auto trimmedReOrientedContigsFnp = njh::files::make_path(metaVelvetDir, "trimmed_reOriented_contigs.fasta");
			SeqOutput outputTrimmedWriter(SeqIOOptions::genFastaOut(trimmedReOrientedContigsFnp));
			auto trimmedReOrientedContigsFnp_belowCutOff = njh::files::make_path(metaVelvetDir,
																																					 "trimmed_reOriented_contigs_belowCutOff.fasta");
			SeqOutput belowCutOffOutputWriter(SeqIOOptions::genFastaOut(trimmedReOrientedContigsFnp_belowCutOff));

			uint32_t belowCutOff = 0;
			uint32_t aboveCutOff = 0;
			bool allPassTrim = true;
			for (const auto &contigsKmerRead: finalSeqs) {
				if (len(contigsKmerRead->seqBase_) < utilityMetaVelvet.inputPars_.minFinalLength_) {
					++belowCutOff;
					belowCutOffOutputWriter.openWrite(contigsKmerRead);
					contigsKmerRead->seqBase_.on_ = false;
				} else {
					MetaDataInName seqMeta(contigsKmerRead->seqBase_.name_);
					trimmedContigInfoOut << contigsKmerRead->seqBase_.name_
															 << "\t" << len(contigsKmerRead->seqBase_)
															 << "\t" << seqMeta.getMeta("estimatedPerBaseCoverage")
															 << std::endl;
					if (!contigsKmerRead->seqBase_.on_) {
						allPassTrim = false;
					} else {
						++aboveCutOff;
					}
					outputAboveCutOffMetaVelWriter.openWrite(contigsKmerRead);
					outputTrimmedWriter.openWrite(contigsKmerRead);
				}
			}
			if (allPassTrim) {
				utilityMetaVelvet.regInfo_->infoCalled_ = true;
				utilityMetaVelvet.regInfo_->uniqHaps_ = aboveCutOff;
			} else {
				utilityMetaVelvet.regInfo_->infoCalled_ = false;
				utilityMetaVelvet.regInfo_->uniqHaps_ = 0;
			}
		}
	} catch (std::exception &e) {
		exceptionMess = e.what();
		utilityMetaVelvet.regInfo_->infoCalled_ = false;
		utilityMetaVelvet.regInfo_->uniqHaps_ = 0;
	}

	{
		outputAboveCutOffMetaVelWriter.closeOut();
		OutputStream basicInfo(njh::files::make_path(utilityMetaVelvet.finalPassDir_, "basicInfoPerRegion.tab.txt"));

		basicInfo << "name\tsuccess\tuniqHaps\treadTotal\treadTotalUsed\ttotalPairedReads";
		basicInfo << "\tsample";
		basicInfo << "\n";
		basicInfo << utilityMetaVelvet.inputPars_.regionUid_;
		basicInfo << "\t" << njh::boolToStr(utilityMetaVelvet.regInfo_->infoCalled_)
							<< "\t" << utilityMetaVelvet.regInfo_->uniqHaps_
							<< "\t" << utilityMetaVelvet.regInfo_->totalReads_
							<< "\t" << utilityMetaVelvet.regInfo_->totalFinalReads_
							<< "\t" << utilityMetaVelvet.regInfo_->totalPairedReads_
							<< "\t" << utilityMetaVelvet.inputPars_.sample_;

		OutputStream exceptionsOut(njh::files::make_path(utilityMetaVelvet.finalPassDir_, "exceptionsMessages.tab.txt"));
		exceptionsOut << "regionUID\tmessage" << std::endl;
		exceptionsOut << utilityMetaVelvet.inputPars_.regionUid_ << "\t" << exceptionMess << std::endl;
	}

	return 0;
}


} // namespace njhseq
