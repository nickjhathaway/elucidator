//
// Created by Nicholas Hathaway on 6/12/23.
//

#include "primerUtilsRunner.hpp"
#include "elucidator/primerUtils/PossibleAmpliconPanelGraph.hpp"
#include <njhseq/IO/SeqIO/SeqInput.hpp>
#include <njhseq/IO/OutputStream.hpp>
#include <njhseq/PrimerIDUtils/PrimerDimerUtils.hpp>
#include <SeekDeep/objects/TarAmpSetupUtils/PrimersAndMids.hpp>
#include <TwoBit/IO/TwoBitFile.hpp>
#include <njhseq/objects/BioDataObject/BLASTHitTabular.hpp>
#include <njhseq/objects/BioDataObject/BioDataFileIO.hpp>

#include <njhseq/concurrency/PairwisePairFactory.hpp>
#include <njhseq/objects/BioDataObject/reading.hpp>
#include <njhseq/objects/dataContainers/tables/TableReader.hpp>


namespace njhseq {





int primerUtilsRunner::creatingMultiplexAmpliconPools(
				const njh::progutils::CmdArgs &inputCommands) {
	uint32_t numThreads = 1;
	uint32_t sliceTop = 10;
	uint32_t unspecificPrimersErrors = 1;
	uint32_t unspecificInsertMinLen = 500;
	uint32_t unspecificPrimerPortionSize = 15;

	uint32_t testIterMax = 0;

	double defaultExpP5 = 0.80;
	std::set <bfs::path> additionalNonDiversityDirs;

	std::set <bfs::path> subsegmentDirs;
	std::set <bfs::path> unspecificAmplificationGenomes;
	std::string primer3ResultDirName;

	double diversityFactor = 2;
	double dimerFactor = 0.5;
	double pairPenaltyFactor = 1;
	double unspecificAmpFactor = 0.50;

	uint32_t hardUnspecAmpCutOff = 20;

	bfs::path genomeGff;
	bfs::path genomeFnp;

	seqSetUp setUp(inputCommands);

	setUp.processVerbose();
	setUp.processDebug();

	setUp.setOption(additionalNonDiversityDirs, "--additionalNonDiversityDirs", "additional Non Diversity Dirs, could be vaccine, drug resistance, diagnostic etc");
	setUp.setOption(defaultExpP5, "--defaultExpP5", "defaultExpP5");

	bool gffSet = setUp.setOption(genomeGff, "--genomeGff", "A path to a genome annotation gff file to compare final panel against");
	setUp.setOption(genomeFnp, "--genomeFnp", "A path to a genome annotation gff file to compare final panel against", gffSet);


	setUp.setOption(diversityFactor, "--diversityFactor", "diversity Factor");
	setUp.setOption(dimerFactor, "--dimerFactor", "dimer Factor");
	setUp.setOption(pairPenaltyFactor, "--pairPenaltyFactor", "pair Penalty Factor");
	setUp.setOption(unspecificAmpFactor, "--unspecificAmpFactor", "unspecific Amp Factor");
	setUp.setOption(hardUnspecAmpCutOff, "--hardUnspecAmpCutOff", "if a primer pair or primer pair match up have more than this amount of unspecific amplification turn completely off");

	setUp.setOption(sliceTop, "--sliceTop", "take this many of the top pairs from each segment");
	setUp.setOption(numThreads, "--numThreads", "number of Threads");
	setUp.setOption(unspecificPrimersErrors, "--unspecificPrimersErrors", "unspecific Primers Errors allowed when searching");
	setUp.setOption(unspecificInsertMinLen, "--unspecificInsertMinLen", "unspecific Insert Min target Length");
	setUp.setOption(unspecificPrimerPortionSize, "--unspecificPrimerPortionSize", "the amount of the 3` end to search for");
	setUp.setOption(testIterMax, "--testIterMax", "how many randomly selected pools to check");

	setUp.setOption(primer3ResultDirName, "--primer3ResultDirName", "primer3 Result Directory Name", true);
	setUp.setOption(subsegmentDirs, "--subsegmentDirs", "subsegment Dirs", true);
	setUp.setOption(unspecificAmplificationGenomes, "--unspecificAmplificationGenomes", "genomes fasta fnps to test unspecific amplification against");

	setUp.processDirectoryOutputName(bfs::basename(primer3ResultDirName) + std::string("_testUnspecific_TODAY"), true);

	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);

	njh::files::checkExistenceThrow(std::vector<bfs::path>(unspecificAmplificationGenomes.begin(), unspecificAmplificationGenomes.end()), __PRETTY_FUNCTION__ );
	njh::files::checkExistenceThrow(std::vector<bfs::path>(subsegmentDirs.begin(), subsegmentDirs.end()), __PRETTY_FUNCTION__ );
	{
		VecStr warnings;
		for(const auto & d : subsegmentDirs){
			auto p3ResDir = njh::files::make_path(d, primer3ResultDirName);
			if(!exists(p3ResDir)){
				warnings.emplace_back(njh::pasteAsStr("directory ", p3ResDir, " doesn't exist"));
			}
		}
		if(!warnings.empty()){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "the following directories don't exist" << "\n";
			ss << njh::conToStr(warnings, "\n") << "\n";
			throw std::runtime_error{ss.str()};
		}
	}
	std::set<bfs::path> additionalNonDiversityDirsWithResults;
	std::set<bfs::path> subsegmentDirsWithResults;
	std::vector<seqInfo> allPrimers;
	std::vector<seqInfo> allUniquePrimers;

	//genomic locations
	std::vector<std::shared_ptr<Bed6RecordCore>> ampLocs;
	std::vector<std::shared_ptr<Bed6RecordCore>> insertLocs;
	std::vector<std::shared_ptr<Bed6RecordCore>> primerLocs;


	{
		//gather all possible primers
		for (const auto &d: subsegmentDirs) {
			auto allPrimersFnp = njh::files::make_path(d, primer3ResultDirName, "allPrimers.fasta");
			njh::files::checkExistenceThrow(allPrimersFnp, __PRETTY_FUNCTION__ );
			if(!njh::files::isFileEmpty(allPrimersFnp)){
				subsegmentDirsWithResults.emplace(d);
				//if empty that means no primers were generated for the input to primer3
				SeqInput reader(SeqIOOptions::genFastaIn(allPrimersFnp));
				reader.openIn();
				seqInfo seq;
				while(reader.readNextRead(seq)){
					seq.name_ = seq.seq_;
					allPrimers.emplace_back(seq);
				}
				addOtherVec(ampLocs, getBeds(njh::files::make_path(d, primer3ResultDirName, "primer3_results_ampLocs.bed")));
				addOtherVec(insertLocs, getBeds(njh::files::make_path(d, primer3ResultDirName, "primer3_results_insertLocs.bed")));
				addOtherVec(primerLocs, getBeds(njh::files::make_path(d, primer3ResultDirName, "primer3_results_primerLocs.bed")));
			}
		}

		//gather any potential non-div targets
		for (const auto &d: additionalNonDiversityDirs) {
			auto allPrimersFnp = njh::files::make_path(d, primer3ResultDirName, "allPrimers.fasta");
			njh::files::checkExistenceThrow(allPrimersFnp, __PRETTY_FUNCTION__ );
			if(!njh::files::isFileEmpty(allPrimersFnp)){
				additionalNonDiversityDirsWithResults.emplace(d);
				//if empty that means no primers were generated for the input to primer3
				SeqInput reader(SeqIOOptions::genFastaIn(allPrimersFnp));
				reader.openIn();
				seqInfo seq;
				while(reader.readNextRead(seq)){
					seq.name_ = seq.seq_;
					allPrimers.emplace_back(seq);
				}
				addOtherVec(ampLocs, getBeds(njh::files::make_path(d, primer3ResultDirName, "primer3_results_ampLocs.bed")));
				addOtherVec(insertLocs, getBeds(njh::files::make_path(d, primer3ResultDirName, "primer3_results_insertLocs.bed")));
				addOtherVec(primerLocs, getBeds(njh::files::make_path(d, primer3ResultDirName, "primer3_results_primerLocs.bed")));
			}
		}
		readVecSorter::sortBySeq(allPrimers, false);
		for (const auto &primer: allPrimers) {
			if (allUniquePrimers.empty() || (allUniquePrimers.back().seq_ != primer.seq_)) {
				allUniquePrimers.emplace_back(primer);
			}
		}
	}

	std::unordered_map<std::string, std::shared_ptr<Bed6RecordCore>> primerPairNameToAmpLoc;
	for(const auto & amp : ampLocs){
		primerPairNameToAmpLoc[amp->name_] = amp;
	}

	std::unordered_map<std::string, uint32_t> primerToIdx;
	for(const auto & en : iter::enumerate(allUniquePrimers)){
		primerToIdx[en.element.seq_] = en.index;
	}
	{
		SeqOutput::write(allUniquePrimers, SeqIOOptions::genFastaOutGz(njh::files::make_path(setUp.pars_.directoryName_, "allUniquePrimers.fasta.gz")));
	}
	//compute dimerization dimerScores
	PrimerDimerUtils dimerScorer;
	std::vector<std::vector<double>> dimerScores = dimerScorer.computeFullScoreMatrix(allUniquePrimers, numThreads);
	{
		OutputStream dimerScoresOut(njh::files::make_path(setUp.pars_.directoryName_, "dimerScores.txt"));
		PrimerDimerUtils::writeMatrix(dimerScores, dimerScoresOut, allUniquePrimers, false);
	}
	{
		OutputStream subsegmentDirsWithResultsOut(njh::files::make_path(setUp.pars_.directoryName_, "subsegmentDirsWithResults.txt"));
		subsegmentDirsWithResultsOut << njh::conToStr(subsegmentDirsWithResults, "\n") << std::endl;
		std::set<bfs::path> dirsWithNoResults;
		for(const auto & d : subsegmentDirs){
			if(!njh::in(d, subsegmentDirsWithResults)){
				dirsWithNoResults.emplace(d);
			}
		}
		OutputStream subsegmentDirsWithNoResultsOut(njh::files::make_path(setUp.pars_.directoryName_, "subsegmentDirsWithNoResults.txt"));
		subsegmentDirsWithNoResultsOut << njh::conToStr(dirsWithNoResults, "\n") << std::endl;
	}

	if(!additionalNonDiversityDirs.empty()){
		OutputStream additionalNonDiversityDirsWithResultsOut(njh::files::make_path(setUp.pars_.directoryName_, "additionalNonDiversityDirsWithResults.txt"));
		additionalNonDiversityDirsWithResultsOut << njh::conToStr(additionalNonDiversityDirsWithResults, "\n") << std::endl;
		std::set<bfs::path> dirsWithNoResults;
		for(const auto & d : additionalNonDiversityDirs){
			if(!njh::in(d, additionalNonDiversityDirsWithResults)){
				dirsWithNoResults.emplace(d);
			}
		}
		OutputStream additionalNonDiversityDirsWithNoResultsOut(njh::files::make_path(setUp.pars_.directoryName_, "additionalNonDiversityDirsWithNoResults.txt"));
		additionalNonDiversityDirsWithNoResultsOut << njh::conToStr(dirsWithNoResults, "\n") << std::endl;
	}

	//summarizing finalResults
	std::map<std::string, uint32_t> primerPairsPerTopRegion;
	std::map<std::string, std::unordered_set<std::string>> primerPairsNamesPerTopRegion;
	std::map<std::string, std::string> primerPairNameToTopRegionKey;

	//testing


	table top;

	for(const auto & d : subsegmentDirsWithResults){
		//diversity
		auto divResFnp = njh::files::make_path(d, primer3ResultDirName, "diversityPerPair/diversityPerPrimerPair.tab.txt");
		table divResTab(divResFnp, "\t", true);

		primerPairsPerTopRegion[d.string()] = divResTab.nRow();
		//primer3 finalResults
		auto primer3ResFnp = njh::files::make_path(d, primer3ResultDirName, "primer3_results.tab.txt");
		table primer3ResTab(primer3ResFnp, "\t", true);

		for(const auto & primerPairName : primer3ResTab.getColumn("SeqIDPrimerPairName")){
			primerPairsNamesPerTopRegion[d.string()].emplace(primerPairName);
			primerPairNameToTopRegionKey[primerPairName] = d.string();
		}
		//add in diversity
		//primer3ResTab.cbind(divResTab, false);
		primer3ResTab = table::cbind(std::vector<table>{primer3ResTab, divResTab}, "SeqIDPrimerPairName", false);
		//1-Expected ploidy 5
		std::vector<double> OneMinusExpP5(primer3ResTab.nRow(), 0);
		for(const auto & row : iter::enumerate(primer3ResTab.content_)){
			OneMinusExpP5[row.index] = 1.0 - njh::StrToNumConverter::stoToNum<double>(row.element[primer3ResTab.getColPos("ExpP5")]);
		}
		primer3ResTab.addColumn(OneMinusExpP5, "1-ExpP5");
		//norm primer penality
		std::vector<double> normPairPenalty(primer3ResTab.nRow(), 0);
		for(const auto & row : iter::enumerate(primer3ResTab.content_)){
			normPairPenalty[row.index] = njh::StrToNumConverter::stoToNum<double>(row.element[primer3ResTab.getColPos("pair_penalty_noSize")])/4.0;
		}
		primer3ResTab.addColumn(normPairPenalty, "norm_pair_penalty_noSize");
		njh::for_each(normPairPenalty,[](double & score){score *=-1;});
		primer3ResTab.addColumn(normPairPenalty, "negative_norm_pair_penalty_noSize");

		std::vector<double> OneMinusExpP5PlusPairPenalty(primer3ResTab.nRow(), 0);
		for(const auto & row : iter::enumerate(primer3ResTab.content_)){
			OneMinusExpP5PlusPairPenalty[row.index] = njh::StrToNumConverter::stoToNum<double>(row.element[primer3ResTab.getColPos("1-ExpP5")]) + njh::StrToNumConverter::stoToNum<double>(row.element[primer3ResTab.getColPos("pair_penalty_noSize")]);
		}
		primer3ResTab.addColumn(OneMinusExpP5PlusPairPenalty, "1-ExpP5+pair_penalty_noSize");
		//add in dimer scores for this pair
		table dimerScoresTab = table(VecStr{"SeqIDPrimerPairName", "forwardPrimerSelfDimer", "reversePrimerSelfDimer", "forwardReverseDimer", "sumDimerScores", "forwardReverseDimerNormScore"});
		for(const auto & row : iter::enumerate(primer3ResTab.content_)){
			auto forwardPrimer = row.element[primer3ResTab.getColPos("left_seq")];
			auto reversePrimer = row.element[primer3ResTab.getColPos("right_seq")];
			auto forwardPrimerSelfDimer = dimerScores[primerToIdx[forwardPrimer]][primerToIdx[forwardPrimer]];
			auto reversePrimerSelfDimer = dimerScores[primerToIdx[reversePrimer]][primerToIdx[reversePrimer]];
			auto forwardReverseDimer = dimerScores[primerToIdx[forwardPrimer]][primerToIdx[reversePrimer]];
			auto forwardReverseDimerNormScore = -1.0 *(forwardReverseDimer + reversePrimerSelfDimer + forwardReverseDimer)/20.0;
			dimerScoresTab.addRow(
							row.element[primer3ResTab.getColPos("SeqIDPrimerPairName")],
							forwardPrimerSelfDimer,
							reversePrimerSelfDimer,
							forwardReverseDimer,
							forwardReverseDimer + reversePrimerSelfDimer + forwardReverseDimer,
							forwardReverseDimerNormScore
							);
			//OneMinusExpP5PlusPairPenalty[row.index] = njh::StrToNumConverter::stoToNum<double>(row.element[primer3ResTab.getColPos("1-ExpP5")]) + njh::StrToNumConverter::stoToNum<double>(row.element[primer3ResTab.getColPos("pair_penalty_noSize")]);
		}
		primer3ResTab = table::cbind(std::vector<table>{primer3ResTab, dimerScoresTab}, "SeqIDPrimerPairName", false);
		//primer3ResTab.cbind(dimerScoresTab, false);
		//add in a score based on 1-ExpP5+pair_penalty_noSize+normDimer
		std::vector<double> OneMinusExpP5PlusPairPenaltyPlusDimer(primer3ResTab.nRow(), 0);
		for (const auto &row: iter::enumerate(primer3ResTab.content_)) {
			OneMinusExpP5PlusPairPenaltyPlusDimer[row.index] =
							njh::StrToNumConverter::stoToNum<double>(row.element[primer3ResTab.getColPos("1-ExpP5")]) +
							njh::StrToNumConverter::stoToNum<double>(row.element[primer3ResTab.getColPos("pair_penalty_noSize")]) -
							njh::StrToNumConverter::stoToNum<double>(row.element[primer3ResTab.getColPos("forwardReverseDimerNormScore")]);

		}
		primer3ResTab.addColumn(OneMinusExpP5PlusPairPenaltyPlusDimer, "1-ExpP5+pair_penalty_noSize+normDimer");

		//add in a score based on ExpP5+norm_pair_penalty_noSize+-normDimer
		std::vector<double> ExpP5MinusNormPairPenaltyMinusDimer(primer3ResTab.nRow(), 0);
		for (const auto &row: iter::enumerate(primer3ResTab.content_)) {
			ExpP5MinusNormPairPenaltyMinusDimer[row.index] =
							diversityFactor * njh::StrToNumConverter::stoToNum<double>(row.element[primer3ResTab.getColPos("ExpP5")]) -
							pairPenaltyFactor * njh::StrToNumConverter::stoToNum<double>(row.element[primer3ResTab.getColPos("norm_pair_penalty_noSize")]) -
							dimerFactor * njh::StrToNumConverter::stoToNum<double>(row.element[primer3ResTab.getColPos("forwardReverseDimerNormScore")]);
		}
		primer3ResTab.addColumn(ExpP5MinusNormPairPenaltyMinusDimer, "ExpP5-norm_pair_penalty_noSize-normDimer");

		primer3ResTab.sortTable("ExpP5-norm_pair_penalty_noSize-normDimer","ExpP5", "negative_norm_pair_penalty_noSize", "forwardReverseDimerNormScore", true);
		primer3ResTab.hasHeader_ = true;
		primer3ResTab.addColumn(VecStr{"diversity"}, "class");
		{
			primer3ResTab.outPutContents(TableIOOpts::genTabFileOut(
							njh::files::make_path(setUp.pars_.directoryName_, bfs::basename(d) + "_primer3_results.tab.txt.gz")));

		}
		std::vector<uint32_t> rowsToSelect(std::min(sliceTop, primer3ResTab.nRow()));
		njh::iota(rowsToSelect, 0U);
		if (top.empty()) {
			top = primer3ResTab.getRows(rowsToSelect);
		} else {
			top.rbind(primer3ResTab.getRows(rowsToSelect), false);
		}
	}


	for(const auto & d : additionalNonDiversityDirsWithResults){
		//diversity


		//primer3 finalResults
		auto primer3ResFnp = njh::files::make_path(d, primer3ResultDirName, "primer3_results.tab.txt");
		table primer3ResTab(primer3ResFnp, "\t", true);
		primerPairsPerTopRegion[d.string()] = primer3ResTab.nRow();

		for(const auto & primerPairName : primer3ResTab.getColumn("SeqIDPrimerPairName")){
			primerPairsNamesPerTopRegion[d.string()].emplace(primerPairName);
			primerPairNameToTopRegionKey[primerPairName] = d.string();
		}
		//add in default ExpP5
		primer3ResTab.addColumn(VecStr{estd::to_string(defaultExpP5)}, "ExpP5");
		primer3ResTab.addColumn(VecStr{estd::to_string(defaultExpP5)}, "he");

		//primer3ResTab.cbind(divResTab, false);
		//1-Expected ploidy 5
		std::vector<double> OneMinusExpP5(primer3ResTab.nRow(), 0);
		for(const auto & row : iter::enumerate(primer3ResTab.content_)){
			OneMinusExpP5[row.index] = 1.0 - njh::StrToNumConverter::stoToNum<double>(row.element[primer3ResTab.getColPos("ExpP5")]);
		}
		primer3ResTab.addColumn(OneMinusExpP5, "1-ExpP5");
		//norm primer penality
		std::vector<double> normPairPenalty(primer3ResTab.nRow(), 0);
		for(const auto & row : iter::enumerate(primer3ResTab.content_)){
			normPairPenalty[row.index] = njh::StrToNumConverter::stoToNum<double>(row.element[primer3ResTab.getColPos("pair_penalty_noSize")])/4.0;
		}
		primer3ResTab.addColumn(normPairPenalty, "norm_pair_penalty_noSize");
		njh::for_each(normPairPenalty,[](double & score){score *=-1;});
		primer3ResTab.addColumn(normPairPenalty, "negative_norm_pair_penalty_noSize");

		std::vector<double> OneMinusExpP5PlusPairPenalty(primer3ResTab.nRow(), 0);
		for(const auto & row : iter::enumerate(primer3ResTab.content_)){
			OneMinusExpP5PlusPairPenalty[row.index] = njh::StrToNumConverter::stoToNum<double>(row.element[primer3ResTab.getColPos("1-ExpP5")]) + njh::StrToNumConverter::stoToNum<double>(row.element[primer3ResTab.getColPos("pair_penalty_noSize")]);
		}
		primer3ResTab.addColumn(OneMinusExpP5PlusPairPenalty, "1-ExpP5+pair_penalty_noSize");
		//add in dimer scores for this pair
		table dimerScoresTab = table(VecStr{"SeqIDPrimerPairName", "forwardPrimerSelfDimer", "reversePrimerSelfDimer", "forwardReverseDimer", "sumDimerScores", "forwardReverseDimerNormScore"});
		for(const auto & row : iter::enumerate(primer3ResTab.content_)){
			auto forwardPrimer = row.element[primer3ResTab.getColPos("left_seq")];
			auto reversePrimer = row.element[primer3ResTab.getColPos("right_seq")];
			auto forwardPrimerSelfDimer = dimerScores[primerToIdx[forwardPrimer]][primerToIdx[forwardPrimer]];
			auto reversePrimerSelfDimer = dimerScores[primerToIdx[reversePrimer]][primerToIdx[reversePrimer]];
			auto forwardReverseDimer = dimerScores[primerToIdx[forwardPrimer]][primerToIdx[reversePrimer]];
			auto forwardReverseDimerNormScore = -1.0 *(forwardReverseDimer + reversePrimerSelfDimer + forwardReverseDimer)/20.0;
			dimerScoresTab.addRow(
							row.element[primer3ResTab.getColPos("SeqIDPrimerPairName")],
							forwardPrimerSelfDimer,
							reversePrimerSelfDimer,
							forwardReverseDimer,
							forwardReverseDimer + reversePrimerSelfDimer + forwardReverseDimer,
							forwardReverseDimerNormScore
			);
			//OneMinusExpP5PlusPairPenalty[row.index] = njh::StrToNumConverter::stoToNum<double>(row.element[primer3ResTab.getColPos("1-ExpP5")]) + njh::StrToNumConverter::stoToNum<double>(row.element[primer3ResTab.getColPos("pair_penalty_noSize")]);
		}
		primer3ResTab = table::cbind(std::vector<table>{primer3ResTab, dimerScoresTab}, "SeqIDPrimerPairName", false);
		//primer3ResTab.cbind(dimerScoresTab, false);
		//add in a score based on 1-ExpP5+pair_penalty_noSize+normDimer
		std::vector<double> OneMinusExpP5PlusPairPenaltyPlusDimer(primer3ResTab.nRow(), 0);
		for (const auto &row: iter::enumerate(primer3ResTab.content_)) {
			OneMinusExpP5PlusPairPenaltyPlusDimer[row.index] =
							njh::StrToNumConverter::stoToNum<double>(row.element[primer3ResTab.getColPos("1-ExpP5")]) +
							njh::StrToNumConverter::stoToNum<double>(row.element[primer3ResTab.getColPos("pair_penalty_noSize")]) -
							njh::StrToNumConverter::stoToNum<double>(row.element[primer3ResTab.getColPos("forwardReverseDimerNormScore")]);

		}
		primer3ResTab.addColumn(OneMinusExpP5PlusPairPenaltyPlusDimer, "1-ExpP5+pair_penalty_noSize+normDimer");

		//add in a score based on ExpP5+norm_pair_penalty_noSize+-normDimer
		std::vector<double> ExpP5MinusNormPairPenaltyMinusDimer(primer3ResTab.nRow(), 0);
		for (const auto &row: iter::enumerate(primer3ResTab.content_)) {
			ExpP5MinusNormPairPenaltyMinusDimer[row.index] =
							diversityFactor * njh::StrToNumConverter::stoToNum<double>(row.element[primer3ResTab.getColPos("ExpP5")]) -
							pairPenaltyFactor * njh::StrToNumConverter::stoToNum<double>(row.element[primer3ResTab.getColPos("norm_pair_penalty_noSize")]) -
							dimerFactor * njh::StrToNumConverter::stoToNum<double>(row.element[primer3ResTab.getColPos("forwardReverseDimerNormScore")]);
		}
		primer3ResTab.addColumn(ExpP5MinusNormPairPenaltyMinusDimer, "ExpP5-norm_pair_penalty_noSize-normDimer");

		primer3ResTab.sortTable("ExpP5-norm_pair_penalty_noSize-normDimer","ExpP5", "negative_norm_pair_penalty_noSize", "forwardReverseDimerNormScore", true);
		primer3ResTab.hasHeader_ = true;
		primer3ResTab.addColumn(VecStr{"diversity"}, "add");
		primer3ResTab.outPutContents(TableIOOpts::genTabFileOut(
						njh::files::make_path(setUp.pars_.directoryName_, bfs::basename(d) + "_primer3_results.tab.txt.gz")));

		std::vector<uint32_t> rowsToSelect(std::min(sliceTop, primer3ResTab.nRow()));
		njh::iota(rowsToSelect, 0U);
		if (top.empty()) {
			top = primer3ResTab.getRows(rowsToSelect);
		} else {
			top.rbind(primer3ResTab.getRows(rowsToSelect), true);
		}
	}




	top.outPutContents(TableIOOpts::genTabFileOut(
					njh::files::make_path(setUp.pars_.directoryName_, "top_primer3_results.tab.txt.gz")));

	auto leftPrimers = top.getColumn("left_seq");
	auto rightPrimers = top.getColumn("right_seq");
	auto topPrimers = getUniqueStrings(njh::catVecs(leftPrimers,rightPrimers));
	{
		auto topPrimersOut = SeqIOOptions::genFastaOut(njh::files::make_path(setUp.pars_.directoryName_, "topPrimers.fasta"));
		SeqOutput writer(topPrimersOut);
		writer.openOut();
		for(const auto & p : topPrimers){
			seqInfo primerSeq(p, p);
			writer.write(primerSeq);
		}
	}

	std::unordered_map<std::string,std::string> primerPairToForward;
	std::unordered_map<std::string,std::string> primerPairToReverse;

	{
		auto primerTab = top.getColumns(VecStr{"SeqIDPrimerPairName", "left_seq", "right_seq"});
		primerTab.columnNames_ = VecStr {"target", "forward", "reverse"};
		primerTab.setColNamePositions();
		for(const auto & row : primerTab){
			primerPairToForward[row[primerTab.getColPos("target")]] = row[primerTab.getColPos("forward")];
			primerPairToReverse[row[primerTab.getColPos("target")]] = row[primerTab.getColPos("reverse")];
		}
		primerTab.outPutContents(TableIOOpts::genTabFileOut(njh::files::make_path(setUp.pars_.directoryName_, "topPrimers.tab.txt")));
	}

	std::unordered_map<std::string, uint32_t> numberOfIndvUnspecificAmps;
	std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>> numberOfBetweenPairsUnspecificAmps;



	table unspecificAmps;
	{
		//testing for unspecific amplification

		std::vector<bfs::path> outputDirs;

		for (const auto &genome: unspecificAmplificationGenomes) {
			auto outputDir = njh::files::make_path(setUp.pars_.directoryName_, "testingPrimerMapping_CheckAgainst_" +
																																				 bfs::path(bfs::basename(genome)).replace_extension("").string());
			std::string cmd = njh::pasteAsStr("elucidator testWithBlastForUnspecificAmplification --errorAllowed ",
																				unspecificPrimersErrors, " --primersFnp ",
																				njh::files::make_path(setUp.pars_.directoryName_, "topPrimers.tab.txt"),
																				" --genomeFnp ", genome,
																				" --overWriteDir --dout ", outputDir,
																				" --minTargetSize ", unspecificInsertMinLen, " --numThreads ",
																				numThreads, " --minLen ", unspecificPrimerPortionSize);
			auto cmdOutput = njh::sys::run({cmd});
			OutputStream cmdOutFnp(njh::files::make_path(setUp.pars_.directoryName_, bfs::path(bfs::basename(genome)).replace_extension("").string()+ "_unspecificRun.json"));
			cmdOutFnp << cmdOutput.toJson() << std::endl;
			if(!cmdOutput.success_){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "in running: " << cmd  << "\n";
				ss << "stdcer: " << "\n";
				ss << cmdOutput.stdErr_ << "\n";
				throw std::runtime_error{ss.str()};
			}
			outputDirs.emplace_back(outputDir);
		}
		for(const auto & output : outputDirs){
			auto extractionInfoFnp = njh::files::make_path(output, "extractionInfos.tsv");
			TableReader treader(TableIOOpts::genTabFileIn(extractionInfoFnp));
			VecStr row;
			while(treader.getNextRow(row)){
				auto p1Name = row[treader.header_.getColPos("p1TarName")];
				auto p2Name = row[treader.header_.getColPos("p2TarName")];
				auto genomicID = njh::pasteAsStr(row[treader.header_.getColPos("chrom")], "-",row[treader.header_.getColPos("start")] , "-",row[treader.header_.getColPos("end")]);
				bool expected = p1Name == p2Name && njh::mapAt(primerPairNameToAmpLoc, p1Name)->genUIDFromCoords() == genomicID;
				if(expected && setUp.pars_.debug_){
					std::cout << p1Name << std::endl;
					std::cout << "\t" << "genomicID: " << genomicID << std::endl;
				}
				if (!expected) {
					if(unspecificAmps.content_.empty()){
						unspecificAmps = treader.header_;
					}
					unspecificAmps.addRow(row);
					if (p1Name == p2Name) {
						++numberOfIndvUnspecificAmps[p1Name];
					} else {
						++numberOfBetweenPairsUnspecificAmps[p1Name][p2Name];
						++numberOfBetweenPairsUnspecificAmps[p2Name][p1Name];
					}
				}
			}
		}
		for(auto & row : unspecificAmps){
			row[unspecificAmps.getColPos("expected")] = "false";
		}
	}
	{
		VecStr primersID;
		for(auto & row : unspecificAmps){
			VecStr primers{row[unspecificAmps.getColPos("p1TarName")], row[unspecificAmps.getColPos("p2TarName")]};
			njh::sort(primers);
			primersID.emplace_back(njh::pasteAsStr(primers.front(), "__", primers.back()));
		}
		unspecificAmps.addColumn(primersID, "primersID");
	}

	if(setUp.pars_.debug_){
		std::cout << "numberOfBetweenPairsUnspecificAmps: " << std::endl;
		for(const auto & numberOfBetweenPairsUnspecificAmp : numberOfBetweenPairsUnspecificAmps){
			for(const auto & other : numberOfBetweenPairsUnspecificAmp.second){
				std::cout << numberOfBetweenPairsUnspecificAmp.first << " " << other.first << " " << other.second << std::endl;
			}
		}
	}
	VecStr unspecificAmpsCnt;
	for(const auto & row : top){
		unspecificAmpsCnt.emplace_back(njh::pasteAsStr(numberOfIndvUnspecificAmps[row[top.getColPos("SeqIDPrimerPairName")]]));
	}
	top.addColumn(unspecificAmpsCnt, "unspecificAmpsCnt");
	std::unordered_map<std::string, double> primerPairsToInvCost;

	{
		//add in a score based on ExpP5+norm_pair_penalty_noSize+-normDimer
		std::vector<double> ExpP5MinusNormPairPenaltyMinusDimerMinusUnspecAmp(top.nRow(), 0);
		for (const auto &row: iter::enumerate(top.content_)) {
			ExpP5MinusNormPairPenaltyMinusDimerMinusUnspecAmp[row.index] =
							diversityFactor * njh::StrToNumConverter::stoToNum<double>(row.element[top.getColPos("ExpP5")]) -
							pairPenaltyFactor * njh::StrToNumConverter::stoToNum<double>(row.element[top.getColPos("norm_pair_penalty_noSize")]) -
							dimerFactor * njh::StrToNumConverter::stoToNum<double>(row.element[top.getColPos("forwardReverseDimerNormScore")]) -
							unspecificAmpFactor * njh::StrToNumConverter::stoToNum<double>(row.element[top.getColPos("unspecificAmpsCnt")]) ;
		}
		top.addColumn(ExpP5MinusNormPairPenaltyMinusDimerMinusUnspecAmp, "ExpP5-norm_pair_penalty_noSize-normDimer-unspecAmp");
	}

	for (const auto &row: top) {
//		primerPairsToInvCost[row[top.getColPos("SeqIDPrimerPairName")]] = njh::StrToNumConverter::stoToNum<double>(
//						row[top.getColPos("ExpP5-norm_pair_penalty_noSize-normDimer")]);
		primerPairsToInvCost[row[top.getColPos("SeqIDPrimerPairName")]] = njh::StrToNumConverter::stoToNum<double>(
						row[top.getColPos("ExpP5-norm_pair_penalty_noSize-normDimer-unspecAmp")]);
	}

	OutputStream primerPairsPerTopRegionOut(njh::files::make_path(setUp.pars_.directoryName_, "primerPairsPerTopRegion.tsv"));
	primerPairsPerTopRegionOut << "regionDir\tpossiblePrimerPairs" << std::endl;
	for(const auto & primerCount : primerPairsPerTopRegion){
		primerPairsPerTopRegionOut << primerCount.first << '\t' << primerCount.second << std::endl;
	}

	VecStr subsegmentDirsWithResultsNames;
	std::transform(subsegmentDirsWithResults.begin(), subsegmentDirsWithResults.end(),
								 std::back_inserter(subsegmentDirsWithResultsNames), [](const bfs::path &fnp) { return fnp.string(); });
	if(!additionalNonDiversityDirsWithResults.empty()){
		std::transform(additionalNonDiversityDirsWithResults.begin(), additionalNonDiversityDirsWithResults.end(),
									 std::back_inserter(subsegmentDirsWithResultsNames), [](const bfs::path &fnp) { return fnp.string(); });
	}

	PossibleAmpliconPanelGraph pool(subsegmentDirsWithResultsNames);
	//add in edges

	auto allPrimerPairs = top.getColumn("SeqIDPrimerPairName");
	{
		PairwisePairFactory pfac(allPrimerPairs.size());
		PairwisePairFactory::PairwisePair pair;
		while(pfac.setNextPair(pair)){
			auto regionRow = primerPairNameToTopRegionKey[allPrimerPairs[pair.row_]];
			auto regionCol = primerPairNameToTopRegionKey[allPrimerPairs[pair.col_]];
			if(regionRow != regionCol){
				auto e = std::make_shared<PossibleAmpliconPanelGraph::edge>();
				e->regionToPrimerName_[regionRow] = allPrimerPairs[pair.row_];
				e->regionToPrimerName_[regionCol] = allPrimerPairs[pair.col_];
				e->primerPairIndvWeights_[allPrimerPairs[pair.row_]] = primerPairsToInvCost[allPrimerPairs[pair.row_]];
				e->primerPairIndvWeights_[allPrimerPairs[pair.col_]] = primerPairsToInvCost[allPrimerPairs[pair.col_]];

				//dimer cost, dimerization between self and between the two primers within the pair are already taken into account with the cost of the primer so can just do the between pairs interactions
				//so the costs would be p1-F -> p2-F, p1-F -< p2-R, p1-R -> p2-F, p1-R -> p2-R
				//so 4 cost in total
				auto p1_forwardPrimer = primerPairToForward[regionRow];
				auto p1_reversePrimer = primerPairToReverse[regionRow];
				auto p2_forwardPrimer = primerPairToForward[regionCol];
				auto p2_reversePrimer = primerPairToReverse[regionCol];
				auto p1_forwardPrimer_vs_p2_forwardPrimer = dimerScores[primerToIdx[p1_forwardPrimer]][primerToIdx[p2_forwardPrimer]];
				auto p1_forwardPrimer_vs_p2_reversePrimer = dimerScores[primerToIdx[p1_forwardPrimer]][primerToIdx[p2_reversePrimer]];
				auto p1_reversePrimer_vs_p2_forwardPrimer = dimerScores[primerToIdx[p1_reversePrimer]][primerToIdx[p2_forwardPrimer]];
				auto p1_reversePrimer_vs_p2_reversePrimer = dimerScores[primerToIdx[p1_reversePrimer]][primerToIdx[p2_reversePrimer]];
				auto dimerNormScore = dimerFactor * ((p1_forwardPrimer_vs_p2_forwardPrimer + p1_forwardPrimer_vs_p2_reversePrimer + p1_reversePrimer_vs_p2_forwardPrimer + p1_reversePrimer_vs_p2_reversePrimer)/26.666667);
				e->primerPairVsPrimerPairWeight_ = dimerNormScore -unspecificAmpFactor * numberOfBetweenPairsUnspecificAmps[allPrimerPairs[pair.row_]][allPrimerPairs[pair.col_]];
				e->regionToOtherRegion_[regionRow] = regionCol;
				e->regionToOtherRegion_[regionCol] = regionRow;
				e->primerPairVsPrimerPairUnSpecAmpCnt_ = numberOfBetweenPairsUnspecificAmps[allPrimerPairs[pair.row_]][allPrimerPairs[pair.col_]];
				if(e->primerPairVsPrimerPairUnSpecAmpCnt_ > hardUnspecAmpCutOff){
					e->on_ = false;
				}
				pool.nodes_[pool.nodeIdx_[regionRow]].edges_.emplace_back(e);
				pool.nodes_[pool.nodeIdx_[regionCol]].edges_.emplace_back(e);
				pool.edges_.emplace_back(e);
			}
		}
	}

	//turn off unspecified hard cut-offs;
	for(const auto & primer : numberOfIndvUnspecificAmps){
		if(primer.second > hardUnspecAmpCutOff){
			pool.turnOffPrimer(primer.first);
//			std::cout << "primer.first: " << primer.first << ", " << primer.second << std::endl;
		}
	}


	pool.greedyDetermineHeaviestPool();
	auto determinedPoolDir = njh::files::make_path(setUp.pars_.directoryName_, "determinedPool");
	njh::files::makeDir(njh::files::MkdirPar{determinedPoolDir});
	Json::Value finalResults;
	{
		auto finalPrimerPairs = pool.getCurrentOnPrimerPairs();
		finalResults["primerPairs"] = njh::json::toJson(finalPrimerPairs);
		finalResults["totalWeight"] = njh::json::toJson(pool.getTotalOnWeight());
		finalResults["subsegmentDirsWithResults"] = njh::json::toJson(subsegmentDirsWithResults);
		finalResults["additionalNonDiversityDirsWithResults"] = njh::json::toJson(additionalNonDiversityDirsWithResults);
		finalResults["subsegmentDirs"] = njh::json::toJson(subsegmentDirs);
		if(setUp.pars_.debug_){
			for(const auto & finalPrimerPair : finalPrimerPairs){
				std::cout << "finalPrimerPair: " << finalPrimerPair << ", " << numberOfIndvUnspecificAmps[finalPrimerPair] << std::endl;
			}

			for(const auto & n : pool.nodes_){
				for(const auto & e : n.edges_){
					std::cout << e->getUid() << ", on: " << njh::colorBool(e->on_) << ", unspecs: " << e->primerPairVsPrimerPairUnSpecAmpCnt_ << std::endl;
				}
			}
			top.outPutContents(TableIOOpts::genTabFileOut(njh::files::make_path(determinedPoolDir, "top_primer3_and_diversity_measures.tab.txt")));
		}
		//big primer3 table
		auto finalTable = top.extractByComp("SeqIDPrimerPairName", [&finalPrimerPairs](const std::string & str){
			return njh::in(str, finalPrimerPairs);
		});
		finalTable.outPutContents(TableIOOpts::genTabFileOut(njh::files::make_path(determinedPoolDir, "primer3_and_diversity_measures.tab.txt")));
		auto pair_penalty_noSize_vec = vecStrToVecNum<double>(finalTable.getColumn("pair_penalty_noSize"));
		auto ExpP5_vec = vecStrToVecNum<double>(finalTable.getColumn("ExpP5"));
		auto he_vec = vecStrToVecNum<double>(finalTable.getColumn("he"));

		//unspecific anps
		VecStr possiblePrimersIds;

		//unspecific amps for self unspecific amps
		auto unspecificAmpsCnt_vec = vecStrToVecNum<double>(finalTable.getColumn("unspecificAmpsCnt"));
		//add in the pairwise
		for(const auto & e : pool.edges_){
			if(e->on_){
				VecStr primers = njh::getVecOfMapValues(e->regionToPrimerName_);
				njh::sort(primers);
				possiblePrimersIds.emplace_back(njh::pasteAsStr(primers.front(), "__", primers.back()));
				unspecificAmpsCnt_vec.emplace_back(e->primerPairVsPrimerPairUnSpecAmpCnt_);
			}
		}
		auto finalUnspecificAmpsTable = unspecificAmps.extractByComp("primersID", [&possiblePrimersIds](const std::string & str){
			return njh::in(str, possiblePrimersIds);
		});

		finalUnspecificAmpsTable.outPutContents(TableIOOpts::genTabFileOut(njh::files::make_path(determinedPoolDir, "unspecificAmplificationInfos.tab.txt")));

		finalResults["pair_penalty_noSize_stats"] = njh::json::toJson(getStatsOnVec(pair_penalty_noSize_vec));
		finalResults["ExpP5_stats"] = njh::json::toJson(getStatsOnVec(ExpP5_vec));
		finalResults["he_stats"] = njh::json::toJson(getStatsOnVec(he_vec));
		finalResults["unspecificAmpsCnt_stats"] = njh::json::toJson(getStatsOnVec(unspecificAmpsCnt_vec));


		//genomic locations for amp, insert and primers
		OutputStream ampLocsOut(njh::files::make_path(determinedPoolDir, "ampLocs.bed"));
		for(const auto & reg : ampLocs){
			if(njh::in(reg->name_, finalPrimerPairs)){
				ampLocsOut << reg->toDelimStrWithExtra() << std::endl;
			}
		}
		OutputStream insertLocsOut(njh::files::make_path(determinedPoolDir, "insertLocs.bed"));
		for(const auto & reg : insertLocs){
			if(njh::in(reg->name_, finalPrimerPairs)){
				insertLocsOut << reg->toDelimStrWithExtra() << std::endl;
			}
		}
		OutputStream primerLocsOut(njh::files::make_path(determinedPoolDir, "primerLocs.bed"));
		for(const auto & reg : primerLocs){
			if(njh::in(reg->name_, finalPrimerPairs)){
				auto outReg = *reg;
				if(outReg.reverseStrand()){
					outReg.name_ += "_R";
				}else{
					outReg.name_ += "_F";
				}
				primerLocsOut << outReg.toDelimStrWithExtra() << std::endl;
			}
		}



		//primer table
		{
			auto primerTab = finalTable.getColumns(VecStr{"SeqIDPrimerPairName", "left_seq", "right_seq"});
			primerTab.columnNames_ = VecStr {"target", "forward", "reverse"};
			primerTab.setColNamePositions();
			primerTab.outPutContents(TableIOOpts::genTabFileOut(njh::files::make_path(determinedPoolDir, "primers.tab.txt")));
			//matrix of primer dimer scores
			std::vector<seqInfo> allFinalPrimerSeqs;
			for(const auto & row : primerTab) {
				allFinalPrimerSeqs.emplace_back(njh::pasteAsStr(row[primerTab.getColPos("target")], "--F"),
																				row[primerTab.getColPos("forward")]);
				allFinalPrimerSeqs.emplace_back(njh::pasteAsStr(row[primerTab.getColPos("target")], "--R"),
																				row[primerTab.getColPos("reverse")]);
			}

			//compute dimerization dimerScores
			std::vector<std::vector<double>> finalDimerScores = dimerScorer.computeFullScoreMatrix(allFinalPrimerSeqs,
																																														 numThreads);

			std::vector<double> all_finalDimerScores;
			for(const auto & scoreVec : finalDimerScores){
				addOtherVec(all_finalDimerScores, scoreVec);
			}
			finalResults["dimerizationScores_stats"] = njh::json::toJson(getStatsOnVec(all_finalDimerScores));

			{
				OutputStream dimerScoresOut(njh::files::make_path(determinedPoolDir, "dimerScores.txt"));
				PrimerDimerUtils::writeMatrix(finalDimerScores, dimerScoresOut, allFinalPrimerSeqs, false);
			}
		}

		OutputStream resultsOut(njh::files::make_path(determinedPoolDir, "determinedPoolResults.json"));
		resultsOut << finalResults << std::endl;
	}

	if(testIterMax > 0){
		auto randomDeterminedPoolDir = njh::files::make_path(determinedPoolDir, "randomTestPools");
		njh::files::makeDir(njh::files::MkdirPar{randomDeterminedPoolDir});
		Json::Value randomResults;
		OutputStream testingStatsOut(njh::files::make_path(randomDeterminedPoolDir, "poolStats.tab.txt"));
		testingStatsOut << "pool\tweight\tExpP5_min\tExpP5_max\tExpP5_mean\tExpP5_median\tExpP5_sum\the_min\the_max\the_mean\the_median\the_sum\tpair_penalty_noSize_min\tpair_penalty_noSize_max\tpair_penalty_noSize_mean\tpair_penalty_noSize_median\tpair_penalty_noSize_sum\tdimerization_min\tdimerization_max\tdimerization_mean\tdimerization_median\tdimerization_sum\tunspecAmpCnts_min\tunspecAmpCnts_max\tunspecAmpCnts_mean\tunspecAmpCnts_median\tunspecAmpCnts_sum" << std::endl;
		randomResults["final"] = finalResults;
		testingStatsOut << "final"

						<< "\t" << finalResults["totalWeight"].asDouble()

						<< "\t" << finalResults["ExpP5_stats"]["min"].asDouble()
						<< "\t" << finalResults["ExpP5_stats"]["max"].asDouble()
						<< "\t" << finalResults["ExpP5_stats"]["mean"].asDouble()
						<< "\t" << finalResults["ExpP5_stats"]["median"].asDouble()
						<< "\t" << finalResults["ExpP5_stats"]["sum"].asDouble()

						<< "\t" << finalResults["he_stats"]["min"].asDouble()
						<< "\t" << finalResults["he_stats"]["max"].asDouble()
						<< "\t" << finalResults["he_stats"]["mean"].asDouble()
						<< "\t" << finalResults["he_stats"]["median"].asDouble()
						<< "\t" << finalResults["he_stats"]["sum"].asDouble()

						<< "\t" << finalResults["pair_penalty_noSize_stats"]["min"].asDouble()
						<< "\t" << finalResults["pair_penalty_noSize_stats"]["max"].asDouble()
						<< "\t" << finalResults["pair_penalty_noSize_stats"]["mean"].asDouble()
						<< "\t" << finalResults["pair_penalty_noSize_stats"]["median"].asDouble()
						<< "\t" << finalResults["pair_penalty_noSize_stats"]["sum"].asDouble()

						<< "\t" << finalResults["dimerizationScores_stats"]["min"].asDouble()
						<< "\t" << finalResults["dimerizationScores_stats"]["max"].asDouble()
						<< "\t" << finalResults["dimerizationScores_stats"]["mean"].asDouble()
						<< "\t" << finalResults["dimerizationScores_stats"]["median"].asDouble()
						<< "\t" << finalResults["dimerizationScores_stats"]["sum"].asDouble()

						<< "\t" << finalResults["unspecificAmpsCnt_stats"]["min"].asDouble()
						<< "\t" << finalResults["unspecificAmpsCnt_stats"]["max"].asDouble()
						<< "\t" << finalResults["unspecificAmpsCnt_stats"]["mean"].asDouble()
						<< "\t" << finalResults["unspecificAmpsCnt_stats"]["median"].asDouble()
						<< "\t" << finalResults["unspecificAmpsCnt_stats"]["sum"].asDouble()

						<< std::endl;

		for(uint32_t testIter = 0; testIter < testIterMax; ++testIter){
			if(setUp.pars_.verbose_){
				std::cout << "testIter: " << testIter << std::endl;
			}
			pool.resetAndGenerateRandomPool();
			auto finalPrimerPairs = pool.getCurrentOnPrimerPairs();
			Json::Value results;

			results["primerPairs"] = njh::json::toJson(finalPrimerPairs);
			results["totalWeight"] = njh::json::toJson(pool.getTotalOnWeight());
			//big primer3 table
			auto finalTable = top.extractByComp("SeqIDPrimerPairName", [&finalPrimerPairs](const std::string & str){
				return njh::in(str, finalPrimerPairs);
			});
			auto norm_pair_penalty_noSize_vec = vecStrToVecNum<double>(finalTable.getColumn("pair_penalty_noSize"));
			auto ExpP5_vec = vecStrToVecNum<double>(finalTable.getColumn("ExpP5"));
			auto he_vec = vecStrToVecNum<double>(finalTable.getColumn("he"));
			//unspecific amps for self unspecific amps
			auto unspecificAmpsCnt_vec = vecStrToVecNum<double>(finalTable.getColumn("unspecificAmpsCnt"));
			//add in the pairwise
			for(const auto & e : pool.edges_){
				if(e->on_){
					unspecificAmpsCnt_vec.emplace_back(e->primerPairVsPrimerPairUnSpecAmpCnt_);
				}
			}
			results["pair_penalty_noSize_stats"] = njh::json::toJson(getStatsOnVec(norm_pair_penalty_noSize_vec));
			results["ExpP5_stats"] = njh::json::toJson(getStatsOnVec(ExpP5_vec));
			results["he_stats"] = njh::json::toJson(getStatsOnVec(he_vec));
			results["unspecificAmpsCnt_stats"] = njh::json::toJson(getStatsOnVec(unspecificAmpsCnt_vec));

			//primer table
			{
				auto primerTab = finalTable.getColumns(VecStr{"SeqIDPrimerPairName", "left_seq", "right_seq"});
				primerTab.columnNames_ = VecStr{"target", "forward", "reverse"};
				primerTab.setColNamePositions();
				//matrix of primer dimer scores
				std::vector<seqInfo> allFinalPrimerSeqs;
				for (const auto &row: primerTab) {
					allFinalPrimerSeqs.emplace_back(njh::pasteAsStr(row[primerTab.getColPos("target")], "--F"),
																					row[primerTab.getColPos("forward")]);
					allFinalPrimerSeqs.emplace_back(njh::pasteAsStr(row[primerTab.getColPos("target")], "--R"),
																					row[primerTab.getColPos("reverse")]);
				}

				//compute dimerization dimerScores
				std::vector<std::vector<double>> finalDimerScores = dimerScorer.computeFullScoreMatrix(allFinalPrimerSeqs,
																																															 numThreads);

				std::vector<double> all_finalDimerScores;
				for (const auto &scoreVec: finalDimerScores) {
					addOtherVec(all_finalDimerScores, scoreVec);
				}
				results["dimerizationScores_stats"] = njh::json::toJson(getStatsOnVec(all_finalDimerScores));

			}
			randomResults[njh::pasteAsStr(testIter)] = results;
			testingStatsOut << testIter
											<< "\t" << results["totalWeight"].asDouble()

											<< "\t" << results["ExpP5_stats"]["min"].asDouble()
											<< "\t" << results["ExpP5_stats"]["max"].asDouble()
											<< "\t" << results["ExpP5_stats"]["mean"].asDouble()
											<< "\t" << results["ExpP5_stats"]["median"].asDouble()
											<< "\t" << results["ExpP5_stats"]["sum"].asDouble()

											<< "\t" << results["he_stats"]["min"].asDouble()
											<< "\t" << results["he_stats"]["max"].asDouble()
											<< "\t" << results["he_stats"]["mean"].asDouble()
											<< "\t" << results["he_stats"]["median"].asDouble()
											<< "\t" << results["he_stats"]["sum"].asDouble()

											<< "\t" << results["pair_penalty_noSize_stats"]["min"].asDouble()
											<< "\t" << results["pair_penalty_noSize_stats"]["max"].asDouble()
											<< "\t" << results["pair_penalty_noSize_stats"]["mean"].asDouble()
											<< "\t" << results["pair_penalty_noSize_stats"]["median"].asDouble()
											<< "\t" << results["pair_penalty_noSize_stats"]["sum"].asDouble()

											<< "\t" << results["dimerizationScores_stats"]["min"].asDouble()
											<< "\t" << results["dimerizationScores_stats"]["max"].asDouble()
											<< "\t" << results["dimerizationScores_stats"]["mean"].asDouble()
											<< "\t" << results["dimerizationScores_stats"]["median"].asDouble()
											<< "\t" << results["dimerizationScores_stats"]["sum"].asDouble()

											<< "\t" << results["unspecificAmpsCnt_stats"]["min"].asDouble()
											<< "\t" << results["unspecificAmpsCnt_stats"]["max"].asDouble()
											<< "\t" << results["unspecificAmpsCnt_stats"]["mean"].asDouble()
											<< "\t" << results["unspecificAmpsCnt_stats"]["median"].asDouble()
							        << "\t" << results["unspecificAmpsCnt_stats"]["sum"].asDouble()

							<< std::endl;
		}
		OutputStream randomResOut(njh::files::make_path(randomDeterminedPoolDir, "poolResults.json"));
		randomResOut << randomResults << std::endl;
	}

	if(gffSet){
		std::stringstream cmd;
		cmd << "cd " << determinedPoolDir << " && elucidator genTargetInfoFromGenomes --longRangeAmplicon --primers primers.tab.txt --genomeDir " << genomeFnp << " --gffDir " << genomeGff << " --dout targetExtraction --useBlast --numThreads " << numThreads << " --errors 1";
		njh::sys::run(VecStr{cmd.str()});
	}

	return 0;
}

} // namespace njhseq
