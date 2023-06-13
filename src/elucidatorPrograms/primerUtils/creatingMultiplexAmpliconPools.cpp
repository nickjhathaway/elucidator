//
// Created by Nicholas Hathaway on 6/12/23.
//

#include "primerUtilsRunner.hpp"
#include <njhseq/IO/SeqIO/SeqInput.hpp>
#include <njhseq/IO/OutputStream.hpp>
#include <njhseq/PrimerIDUtils/PrimerDimerUtils.hpp>
#include <SeekDeep/objects/TarAmpSetupUtils/PrimersAndMids.hpp>
#include <TwoBit/IO/TwoBitFile.hpp>
#include <njhseq/objects/BioDataObject/BLASTHitTabular.hpp>
#include <njhseq/objects/BioDataObject/BioDataFileIO.hpp>
#include <njhseq/BamToolsUtils/ReAlignedSeq.hpp>
#include <njhseq/GenomeUtils/GenomeExtraction/ParsingAlignmentInfo/GenomeExtractResult.hpp>
#include <njhseq/objects/BioDataObject/BioRecordsUtils/BedUtility.hpp>


namespace njhseq {
int primerUtilsRunner::creatingMultiplexAmpliconPools(
				const njh::progutils::CmdArgs &inputCommands) {
	uint32_t numThreads = 1;
	uint32_t sliceTop = 10;
	std::set <bfs::path> subsegmentDirs;
	std::string primer3ResultDirName;
	seqSetUp setUp(inputCommands);
	setUp.setOption(sliceTop, "--sliceTop", "take this many of the top pairs from each segment");
	setUp.setOption(numThreads, "--numThreads", "number of Threads");
	setUp.setOption(primer3ResultDirName, "--primer3ResultDirName", "primer3 Result Directory Name", true);
	setUp.setOption(subsegmentDirs, "--subsegmentDirs", "subsegment Dirs", true);

	setUp.processDirectoryOutputName(bfs::basename(primer3ResultDirName) + std::string("_testUnspecific_TODAY"), true);

	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);

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
	std::set<bfs::path> subsegmentDirsWithResults;
	std::vector<seqInfo> allPrimers;
	std::vector<seqInfo> allUniquePrimers;
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
			}
		}
		readVecSorter::sortBySeq(allPrimers, false);
		for (const auto &primer: allPrimers) {
			if (allUniquePrimers.empty()) {
				allUniquePrimers.emplace_back(primer);
			} else if (allUniquePrimers.back().seq_ != primer.seq_) {
				allUniquePrimers.emplace_back(primer);
			}
		}
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

	//summarizing results
	std::map<std::string, uint32_t> primerPairsPerTopRegion;
	//testing
	double diversityFactor = 2;
	double dimerFactor = 0.5;
	double pairPenaltyFactor = 1;

	table top;

	for(const auto & d : subsegmentDirsWithResults){
		//diversity
		auto divResFnp = njh::files::make_path(d, primer3ResultDirName, "diversityPerPair/diversityPerPrimerPair.tab.txt");
		table divResTab(divResFnp, "\t", true);

		primerPairsPerTopRegion[d.string()] = divResTab.nRow();
		//primer3 results
		auto primer3ResFnp = njh::files::make_path(d, primer3ResultDirName, "primer3_results.tab.txt");
		table primer3ResTab(primer3ResFnp, "\t", true);
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
			auto forwardReverseDimerNormScore = (forwardReverseDimer + reversePrimerSelfDimer + forwardReverseDimer)/20.0;
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
		primer3ResTab.outPutContents(TableIOOpts::genTabFileOut(
						njh::files::make_path(setUp.pars_.directoryName_, bfs::basename(d) + "_primer3_results.tab.txt.gz")));
		std::vector<uint32_t> rowsToSelect(std::min(sliceTop, primer3ResTab.nRow()));
		njh::iota(rowsToSelect, 0U);
		if (top.empty()) {
			top = primer3ResTab.getRows(rowsToSelect);
		} else {
			top.rbind(primer3ResTab.getRows(rowsToSelect), false);
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

	{
		auto primerTab = top.getColumns(VecStr{"SeqIDPrimerPairName", "left_seq", "right_seq"});
		primerTab.columnNames_ = VecStr {"target", "forward", "reverse"};
		primerTab.setColNamePositions();
		primerTab.outPutContents(TableIOOpts::genTabFileOut(njh::files::make_path(setUp.pars_.directoryName_, "topPrimers.tsv")));
	}

	OutputStream primerPairsPerTopRegionOut(njh::files::make_path(setUp.pars_.directoryName_, "primerPairsPerTopRegion.tsv"));
	primerPairsPerTopRegionOut << "regionDir\tpossiblePrimerPairs" << std::endl;
	for(const auto & primerCount : primerPairsPerTopRegion){
		primerPairsPerTopRegionOut << primerCount.first << '\t' << primerCount.second << std::endl;
	}
	return 0;
}

} // namespace njhseq
