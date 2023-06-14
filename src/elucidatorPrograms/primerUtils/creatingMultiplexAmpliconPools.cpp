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
#include <njhseq/concurrency/PairwisePairFactory.hpp>


namespace njhseq {


/**
 * @brief Undirected weighted graph, each node is a different region, each edge is a different set of primer pairs combo with a cost of having those primers together
 */
class PossibleAmpliconPanelGraph{

public:

	class edge;

	class node{
	public:
		explicit node(std::string region): region_(std::move(region)){

		}
		std::string region_;
		bool visited_{false};

		std::vector<std::shared_ptr<edge>> edges_;

		void turnOnAllEdges() {
			for (auto &e: edges_) {
				e->on_ = true;
			}
		}

		void reset() {
			visited_ = false;
			turnOnAllEdges();
		}

		void sortEdgesMaxTotalWeightTop(){
			njh::sort(edges_, [](const std::shared_ptr<edge> & e1, const std::shared_ptr<edge> & e2){
				return e1->getTotalWeights() > e2->getTotalWeights();
			});
		}

		void sortEdgesMaxRegionPrimerWeightTop(){
			njh::sort(edges_, [this](const std::shared_ptr<edge> & e1, const std::shared_ptr<edge> & e2){
				if(e1->primerPairIndvWeights_[e1->regionToPrimerName_[this->region_]] == e2->primerPairIndvWeights_[e2->regionToPrimerName_[this->region_]]){
					return e1->getTotalWeights() > e2->getTotalWeights();
				}else{
					return e1->primerPairIndvWeights_[e1->regionToPrimerName_[this->region_]] > e2->primerPairIndvWeights_[e2->regionToPrimerName_[this->region_]];
				}
			});
		}

		std::string getTopOnPrimerMaxTotalWeightOnTop(){
			sortEdgesMaxTotalWeightTop();
			if(edges_.empty()){
				return {""};
			}
			return edges_.front()->regionToPrimerName_[region_];
		}

		void toggleOnEdgesForPrimer(const std::string & name){
			for(const auto & e : edges_){
				e->on_ = e->regionToPrimerName_[region_] == name;
			}
		}

		void turnOffAllOtherPrimers(const std::string & name){
			for(const auto & e : edges_){
				if(e->regionToPrimerName_[region_] != name){
					e->on_ = false;
				}
			}
		}

		/**
		 * @brief turn off all edges to the other region in the keepEdge expect for the keepEdge
		 * @param keepEdge the edge to keep, turn off all other edges going to this region
		 */
		void turnOfOtherEdgesForOtherRegion(const edge &keepEdge) {
			for (auto &e: edges_) {
				if (njh::in(keepEdge.regionToOtherRegion_.at(region_), e->regionToPrimerName_) &&
						e->regionToPrimerName_[keepEdge.regionToOtherRegion_.at(region_)] !=
						keepEdge.regionToPrimerName_.at(keepEdge.regionToOtherRegion_.at(region_))) {
					e->on_ = false;
				}
			}
		}

		[[nodiscard]] VecStr getUniquePrimersOn() const {
			auto onCounts = getTotalOnEdgesForPrimers();
			return njh::getVecOfMapKeys(onCounts);
		}

		[[nodiscard]] uint32_t getTotalUniquePrimersOn() const {
			auto onCounts = getTotalOnEdgesForPrimers();
			return onCounts.size();
		}

		[[nodiscard]] std::unordered_map<std::string, uint32_t> getTotalOnEdgesForPrimers() const {
			std::unordered_map<std::string, uint32_t> ret;
			for (const auto &e: edges_) {
				if (e->on_) {
					++ret[e->regionToPrimerName_[region_]];
				}
			}
			return ret;
		}

		[[nodiscard]] std::unordered_map<std::string, uint32_t> getTotalOnEdgesForOtherRegions() const {
			std::unordered_map<std::string, uint32_t> ret;
			for (const auto &e: edges_) {
				if (e->on_) {
					++ret[e->regionToOtherRegion_[region_]];
				}
			}
			return ret;
		}
	};

	class edge {
	public:
		std::unordered_map<std::string, std::string> regionToPrimerName_;/*< key is region name, key is primer pair name*/
		std::unordered_map<std::string, std::string> regionToOtherRegion_;/*< key is region name, key is the other region in this edge*/

		bool on_{true};
		double primerPairVsPrimerPairWeight_{0};
		std::unordered_map<std::string, double> primerPairIndvWeights_;

		[[nodiscard]] std::string getUid() const{
			auto pairs = njh::getVecOfMapValues(regionToPrimerName_);
			njh::sort(pairs);
			return njh::pasteAsStr(pairs);
		}

		[[nodiscard]] double sumInvWeights() const {
			double ret = 0;
			for (const auto &pweight: primerPairIndvWeights_) {
				ret += pweight.second;
			}
			return ret;
		}

		[[nodiscard]] double getTotalWeights() const {
			return primerPairVsPrimerPairWeight_ + sumInvWeights();
		}
	};

	explicit PossibleAmpliconPanelGraph(const VecStr & regions){
		for(const auto & r : regions){
			nodes_.emplace_back(r);
		}
		setNodeIdx();
	}
	void setNodeIdx(){
		nodeIdx_.clear();
		for(const auto & n : iter::enumerate(nodes_)){
			nodeIdx_[n.element.region_] = n.index;
		}
	}

	void resetAll(){
		njh::for_each(nodes_,[](auto & n){
			n.reset();
		});
	}

	void resetAndGenerateRandomPool(){
		resetAll();
		njh::randomGenerator rgen;
		for(auto & n : nodes_){
			auto randomPrimer = rgen.unifRandSelection(n.getUniquePrimersOn());
			n.turnOffAllOtherPrimers(randomPrimer);
		}
	}

	std::vector<std::shared_ptr<PossibleAmpliconPanelGraph::edge>> greedyDetermineHeaviestPool(){
		auto keepGoingFunc = [this](){
			bool keepGoing = false;
			for(const auto & n : nodes_){
				if(n.getTotalUniquePrimersOn() != 1){
					keepGoing = true;
					break;
				}
			}
			return keepGoing;
		};
		std::vector<std::shared_ptr<PossibleAmpliconPanelGraph::edge>> topEdges;
		std::unordered_set<std::string> topEdgeUids;
		std::shared_ptr<PossibleAmpliconPanelGraph::edge> topEdge = nullptr;
		while(keepGoingFunc()){
			topEdge = nullptr;
			for(auto & n : nodes_){
				for(const auto & e : n.edges_){
					if(e->on_ && !njh::in(e->getUid(), topEdgeUids)){
						if(nullptr == topEdge || e->getTotalWeights() > topEdge->getTotalWeights()){
							topEdge = e;
						}
					}
				}
				//turn off all other primers
				nodes_[nodeIdx_[topEdge->regionToOtherRegion_.begin()->first]].turnOffAllOtherPrimers(topEdge->regionToPrimerName_[topEdge->regionToOtherRegion_.begin()->first]);
				nodes_[nodeIdx_[topEdge->regionToOtherRegion_[topEdge->regionToOtherRegion_.begin()->first]]].turnOffAllOtherPrimers(topEdge->regionToPrimerName_[topEdge->regionToOtherRegion_[topEdge->regionToOtherRegion_.begin()->first]]);
				//turn off all other connections between these two regions
				nodes_[nodeIdx_[topEdge->regionToOtherRegion_.begin()->first]].turnOfOtherEdgesForOtherRegion(*topEdge);
				nodes_[nodeIdx_[topEdge->regionToOtherRegion_[topEdge->regionToOtherRegion_.begin()->first]]].turnOfOtherEdgesForOtherRegion(*topEdge);
			}
/*			for (const auto &n: nodes_) {
				auto onCounts = n.getTotalOnEdgesForOtherRegions();
				auto primersOn = n.getTotalUniquePrimersOn();
				std::cout << n.region_ << ", primers on: " << primersOn << std::endl;
				for (const auto &on: onCounts) {
					std::cout << "\t" << on.first << " " << on.second << std::endl;
				}
				std::cout << std::endl;
			}*/
			topEdges.emplace_back(topEdge);
			topEdgeUids.emplace(topEdge->getUid());
		}
		return topEdges;
	}

	double getTotalOnWeight(){
		double ret = 0;
		for(const auto & n : nodes_){
			for(const auto & e : n.edges_){
				if(e->on_){
					ret += e->getTotalWeights();
				}
			}
		}
		return ret;
	}

	std::vector<node> nodes_;
	std::unordered_map<std::string, uint32_t> nodeIdx_;
};


int primerUtilsRunner::creatingMultiplexAmpliconPools(
				const njh::progutils::CmdArgs &inputCommands) {
	uint32_t numThreads = 1;
	uint32_t sliceTop = 10;
	uint32_t unspecificPrimersErrors = 1;
	uint32_t unspecificInsertMinLen = 500;
	uint32_t unspecificPrimerPortionSize = 14;
	std::set <bfs::path> subsegmentDirs;
	std::set <bfs::path> unspecificAmplificationGenomes;
	std::string primer3ResultDirName;
	seqSetUp setUp(inputCommands);
	setUp.setOption(sliceTop, "--sliceTop", "take this many of the top pairs from each segment");
	setUp.setOption(numThreads, "--numThreads", "number of Threads");
	setUp.setOption(unspecificPrimersErrors, "--unspecificPrimersErrors", "unspecific Primers Errors allowed when searching");
	setUp.setOption(unspecificInsertMinLen, "--unspecificInsertMinLen", "unspecific Insert Min target Length");
	setUp.setOption(unspecificPrimerPortionSize, "--unspecificPrimerPortionSize", "the amount of the 3` end to search for");

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
	std::map<std::string, std::unordered_set<std::string>> primerPairsNamesPerTopRegion;
	std::map<std::string, std::string> primerPairNameToTopRegionKey;

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

	{
		//testing for unspecific amplification
		for (const auto &genome: unspecificAmplificationGenomes) {
			std::string cmd = njh::pasteAsStr("elucidator testWithBlastForUnspecificAmplification --errorAllowed ",
																				unspecificPrimersErrors, " --primersFnp ",
																				njh::files::make_path(setUp.pars_.directoryName_, "topPrimers.tab.txt"),
																				" --genomeFnp ", genome,
																				" --overWriteDir --dout ", njh::files::make_path(setUp.pars_.directoryName_, "testingPrimerMapping_CheckAgainst_" +
																				bfs::path(bfs::basename(genome)).replace_extension("").string()),
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
		}
	}

	std::unordered_map<std::string, double> primerPairsToInvCost;

	for (const auto &row: top) {
		primerPairsToInvCost[row[top.getColPos("SeqIDPrimerPairName")]] = njh::StrToNumConverter::stoToNum<double>(
						row[top.getColPos("ExpP5-norm_pair_penalty_noSize-normDimer")]);
	}

	OutputStream primerPairsPerTopRegionOut(njh::files::make_path(setUp.pars_.directoryName_, "primerPairsPerTopRegion.tsv"));
	primerPairsPerTopRegionOut << "regionDir\tpossiblePrimerPairs" << std::endl;
	for(const auto & primerCount : primerPairsPerTopRegion){
		primerPairsPerTopRegionOut << primerCount.first << '\t' << primerCount.second << std::endl;
	}

	VecStr subsegmentDirsWithResultsNames;
	std::transform(subsegmentDirsWithResults.begin(), subsegmentDirsWithResults.end(),
								 std::back_inserter(subsegmentDirsWithResultsNames), [](const bfs::path &fnp) { return fnp.string(); });

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
				e->primerPairVsPrimerPairWeight_ = -dimerNormScore;
				e->regionToOtherRegion_[regionRow] = regionCol;
				e->regionToOtherRegion_[regionCol] = regionRow;
				pool.nodes_[pool.nodeIdx_[regionRow]].edges_.emplace_back(e);
				pool.nodes_[pool.nodeIdx_[regionCol]].edges_.emplace_back(e);
			}
		}
	}
/*
	std::shared_ptr<PossibleAmpliconPanelGraph::edge> topEdge = nullptr;
	for(auto & n : pool.nodes_){
		for(const auto & e : n.edges_){
			if(e->on_){
				if(nullptr == topEdge || e->getTotalWeights() > topEdge->getTotalWeights()){
					topEdge = e;
				}
			}
		}
	}
	std::cout << std::endl;
	std::cout << "top edge: " << std::endl;
	{
		for(const auto & name : topEdge->regionToPrimerName_){
			std::cout << "\tregion to pair: " << name.first << ": " << name.second << std::endl;
		}
		for(const auto & cost : topEdge->primerPairIndvWeights_){
			std::cout << "\tpair: " << cost.first << ": " << cost.second << std::endl;
		}
		std::cout << "\ttotalWeight: " << topEdge->getTotalWeights() << std::endl;
		std::cout << "\ton: " << njh::colorBool(topEdge->on_) << std::endl;
	}

	//turn off all other primers
	pool.nodes_[pool.nodeIdx_[topEdge->regionToOtherRegion_.begin()->first]].toggleOnEdgesForPrimer(topEdge->regionToPrimerName_[topEdge->regionToOtherRegion_.begin()->first]);
	pool.nodes_[pool.nodeIdx_[topEdge->regionToOtherRegion_[topEdge->regionToOtherRegion_.begin()->first]]].toggleOnEdgesForPrimer(topEdge->regionToPrimerName_[topEdge->regionToOtherRegion_[topEdge->regionToOtherRegion_.begin()->first]]);
	//turn off all other connections between these two regions
	pool.nodes_[pool.nodeIdx_[topEdge->regionToOtherRegion_.begin()->first]].turnOfOtherEdgesForOtherRegion(*topEdge);
	pool.nodes_[pool.nodeIdx_[topEdge->regionToOtherRegion_[topEdge->regionToOtherRegion_.begin()->first]]].turnOfOtherEdgesForOtherRegion(*topEdge);

	for (const auto &n: pool.nodes_) {
		auto onCounts = n.getTotalOnEdgesForOtherRegions();
		auto primersOn = n.getTotalUniquePrimersOn();
		std::cout << n.region_ << ", primers on: " << primersOn << std::endl;
		for (const auto &on: onCounts) {
			std::cout << "\t" << on.first << " " << on.second << std::endl;
		}
		std::cout << std::endl;
	}

	std::shared_ptr<PossibleAmpliconPanelGraph::edge> nextTopEdge = nullptr;
	for(auto & n : pool.nodes_){
		for(const auto & e : n.edges_){
			if(e->on_ && e->getUid() != topEdge->getUid()){
				if(nullptr == nextTopEdge || e->getTotalWeights() > nextTopEdge->getTotalWeights()){
					nextTopEdge = e;
				}
			}
		}
	}
	std::cout << std::endl;
	std::cout << "next top edge: " << std::endl;
	{
		for(const auto & name : nextTopEdge->regionToPrimerName_){
			std::cout << "\tregion to pair: " << name.first << ": " << name.second << std::endl;
		}
		for(const auto & cost : nextTopEdge->primerPairIndvWeights_){
			std::cout << "\tpair: " << cost.first << ": " << cost.second << std::endl;
		}
		std::cout << "\ttotalWeight: " << nextTopEdge->getTotalWeights() << std::endl;
		std::cout << "\ton: " << njh::colorBool(nextTopEdge->on_) << std::endl;
	}

	//turn off all other primers
	pool.nodes_[pool.nodeIdx_[nextTopEdge->regionToOtherRegion_.begin()->first]].toggleOnEdgesForPrimer(nextTopEdge->regionToPrimerName_[nextTopEdge->regionToOtherRegion_.begin()->first]);
	pool.nodes_[pool.nodeIdx_[nextTopEdge->regionToOtherRegion_[nextTopEdge->regionToOtherRegion_.begin()->first]]].toggleOnEdgesForPrimer(nextTopEdge->regionToPrimerName_[nextTopEdge->regionToOtherRegion_[nextTopEdge->regionToOtherRegion_.begin()->first]]);
	//turn off all other connections between these two regions
	pool.nodes_[pool.nodeIdx_[nextTopEdge->regionToOtherRegion_.begin()->first]].turnOfOtherEdgesForOtherRegion(*nextTopEdge);
	pool.nodes_[pool.nodeIdx_[nextTopEdge->regionToOtherRegion_[nextTopEdge->regionToOtherRegion_.begin()->first]]].turnOfOtherEdgesForOtherRegion(*nextTopEdge);
*/

	auto topEdges = pool.greedyDetermineHeaviestPool();
	for (const auto &n: pool.nodes_) {
		auto onCounts = n.getTotalOnEdgesForOtherRegions();
		auto primersOn = n.getTotalUniquePrimersOn();
		std::cout << n.region_ << ", primers on: " << primersOn << std::endl;
		for (const auto &on: onCounts) {
			std::cout << "\t" << on.first << " " << on.second << std::endl;
		}
		std::cout << std::endl;
	}
	uint32_t topEdgeCount = 0;
	std::set<std::string> finalPrimerPairsFromTopEdges;
	for(const auto & topEdge : topEdges){
		njh::addVecToSet(njh::getVecOfMapValues(topEdge->regionToPrimerName_), finalPrimerPairsFromTopEdges);
		std::cout << topEdgeCount++ << std::endl;
		for(const auto & name : topEdge->regionToPrimerName_){
			std::cout << "\tregion to pair: " << name.first << ": " << name.second << std::endl;
		}
		for(const auto & cost : topEdge->primerPairIndvWeights_){
			std::cout << "\tpair: " << cost.first << ": " << cost.second << std::endl;
		}
		std::cout << "\ttotalWeight: " << topEdge->getTotalWeights() << std::endl;
		std::cout << "\ton: " << njh::colorBool(topEdge->on_) << std::endl;
	}
	std::set<std::string> finalPrimerPairsFrom;
	topEdgeCount = 0;
	for(const auto & n : pool.nodes_){
		for(const auto & e : n.edges_){
			if(e->on_){
				njh::addVecToSet(njh::getVecOfMapValues(e->regionToPrimerName_), finalPrimerPairsFrom);
				std::cout << topEdgeCount++ << std::endl;
				for(const auto & name : e->regionToPrimerName_){
					std::cout << "\tregion to pair: " << name.first << ": " << name.second << std::endl;
				}
				for(const auto & cost : e->primerPairIndvWeights_){
					std::cout << "\tpair: " << cost.first << ": " << cost.second << std::endl;
				}
				std::cout << "\ttotalWeight: " << e->getTotalWeights() << std::endl;
				std::cout << "\ton: " << njh::colorBool(e->on_) << std::endl;
			}
		}
	}
	std::cout << "pool.nodes_.size() : " << pool.nodes_.size()  << std::endl;
	std::cout << "subsegmentDirsWithResults.size() : " << subsegmentDirsWithResults.size()  << std::endl;
	std::cout << "finalPrimerPairsFromTopEdges.size() : " << finalPrimerPairsFromTopEdges.size()  << std::endl;
	std::cout << njh::conToStr(finalPrimerPairsFromTopEdges, "\n") << std::endl;
	std::cout << "finalPrimerPairsFrom.size() : " << finalPrimerPairsFrom.size()  << std::endl;
	std::cout << njh::conToStr(finalPrimerPairsFrom, "\n") << std::endl;
	std::cout << "pool.getTotalOnWeight(): " << pool.getTotalOnWeight() << std::endl;


	for(uint32_t testIter = 0; testIter < 10; ++testIter){
		std::cout << "test random pool selection: " << testIter << std::endl;
		pool.resetAndGenerateRandomPool();
		std::cout << "pool.getTotalOnWeight(): " << pool.getTotalOnWeight() << std::endl;
	}

	return 0;
}

} // namespace njhseq
