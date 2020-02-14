/*
 * genExp_evaluateContigsWithBowtie2.cpp
 *
 *  Created on: Feb 5, 2020
 *      Author: nicholashathaway
 */



#include <TwoBit.h>
#include "genExp.hpp"
#include "elucidator/BamToolsUtils.h"
#include "elucidator/objects/BioDataObject.h"
#include "elucidator/objects/dataContainers/graphs.h"
#include <njhseq/GenomeUtils.h>
#include "elucidator/concurrency/LockableJsonLog.hpp"




namespace njhseq {




int genExpRunner::evaluateContigsWithBowtie2(const njh::progutils::CmdArgs & inputCommands){
	std::string genomesStr = "";
	uint32_t numThreads = 1;
	bfs::path bedFnp = "";
	bfs::path idsRequired = "";
	std::string program = "";
	std::string sample  = "";
	comparison amountOfErrorForCoverageCalc;
	MultiGenomeMapper::inputParameters mapperPars;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processComparison(amountOfErrorForCoverageCalc);
  setUp.setOption(program, "--program", "Name of the program to output with the report", true);
  setUp.setOption(sample,  "--sample",  "Name of the sample to output with the report",  true);

	setUp.setOption(numThreads, "--numThreads", "number of cpus to use");
  setUp.setOption(genomesStr, "--genomes", "Names of the genomes to extract from (should not have extension, e.g. Pf3d7 for Pf3d7.fasta");
  setUp.setOption(mapperPars.genomeDir_, "--genomeDir", "Names of the genome directory where the genomes are stored", true);
  setUp.setOption(bedFnp, "--bed", "A bed file to report coverage of");
	setUp.setOption(idsRequired,      "--idsRequired",       "IDs Required, 2 columns, 1)genome,2)ids, ids should be comma separated", "" == bedFnp);

  setUp.setOption(mapperPars.gffDir_, "--gffDir", "A directory of gffs that go with the genomes", "" == bedFnp);
  setUp.setOption(mapperPars.gffIntersectPars_.extraAttributes_, "--gffExtraAttributes", "gff Extra Attributes");

  setUp.processReadInNames();
  setUp.processDirectoryOutputName(true);
	setUp.finishSetUp(std::cout);

	setUp.startARunLog(setUp.pars_.directoryName_);
	BioCmdsUtils bRunner(setUp.pars_.verbose_);

	MultiGenomeMapper gMapper(mapperPars);
	gMapper.pars_.numThreads_ = numThreads;
	if("" != genomesStr){
		gMapper.setSelectedGenomes(njh::strToSet<std::string>(genomesStr, ","));
	}
	gMapper.loadInGenomes();
	gMapper.setUpGenomes();

	std::unordered_map<std::string, std::string> chromosomeToGenome;
	{
		for (const auto & genome : gMapper.genomes_) {
			TwoBit::TwoBitFile treader(genome.second->fnpTwoBit_);
			auto seqNames = treader.sequenceNames();
			for(const auto & seqName : seqNames){
				chromosomeToGenome[seqName] = genome.first;
			}
		}
	}


	std::vector<GenomicRegion> requiredRegions;
	{
		//gathering the expected regions;
		if("" != bedFnp){
			auto beds = getBeds(bedFnp);
			requiredRegions = bedPtrsToGenomicRegs(beds);
			//check to make sure they are regions that can be found in the input genomes;
			std::set<std::string> chromosomes;
			for(const auto & region : requiredRegions){
				chromosomes.emplace(region.chrom_);
			}
			std::set<std::string> missingChromosomes;
			for(const auto & chromosome : chromosomes){
				bool found = false;
				for(const auto & genome : gMapper.genomes_){
					TwoBit::TwoBitFile treader(genome.second->fnpTwoBit_);
					auto seqNames = treader.sequenceNames();
					if(njh::in(chromosome, seqNames)){
						found = true;
						break;
					}
				}
				if(!found){
					missingChromosomes.emplace(chromosome);
				}
			}
			if(!missingChromosomes.empty()){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "missing the following chromosomes: " << njh::conToStrEndSpecial(missingChromosomes, ", ", " and ")<< "\n";
				throw std::runtime_error{ss.str()};
			}
		} else {
			// read in the required
			table requiredTab(idsRequired, "\t", true);
			requiredTab.checkForColumnsThrow(VecStr{"genome", "ids"}, __PRETTY_FUNCTION__);
			std::set<std::string> genomesSet;
			njh::addVecToSet(requiredTab.getColumnLevels("genome"), genomesSet);
			VecStr ids;
			VecStr foundIds;
			for(const auto & row : requiredTab){
				addOtherVec(ids, tokenizeString(row[requiredTab.getColPos("ids")], ","));
			}

			VecStr missingGenomes;
			for(const auto & genome : genomesSet){
				if(!njh::in(genome, gMapper.genomes_)){
					missingGenomes.emplace_back(genome);
				}
			}

			if(!missingGenomes.empty()){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "missing the following genomes: " << njh::conToStrEndSpecial(missingGenomes, ", ", " and ")<< "\n";
				throw std::runtime_error{ss.str()};
			}

			//get the genes;
			for(const auto & genoome : genomesSet){
				auto gffFnp = gMapper.genomes_.at(genoome)->gffFnp_;
				BioDataFileIO<GFFCore> gffIo{IoOptions(InOptions(gffFnp))};
				gffIo.openIn();
				GFFCore gff;
				std::string line = "";
				while(gffIo.readNextRecord(gff)){
					if("gene" == gff.type_&& gff.hasAttr("ID") && njh::in(gff.getAttr("ID"), ids)){
						requiredRegions.emplace_back(GenomicRegion(gff));
						foundIds.emplace_back(gff.getAttr("ID"));
					}
					bool end = false;
					while('#' == gffIo.inFile_->peek()){
						if (njh::files::nextLineBeginsWith(*gffIo.inFile_, "##FASTA")) {
							end = true;
							break;
						}
						njh::files::crossPlatGetline(*gffIo.inFile_, line);
					}
					if(end){
						break;
					}
				}
			}
			VecStr missingIDs;
			for(const auto & id : ids){
				if(!njh::in(id, foundIds)){
					missingIDs.emplace_back(id);
				}
			}
			if(!missingIDs.empty()){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "missing the following IDs: " << njh::conToStrEndSpecial(missingIDs, ", ", " and ")<< "\n";
				throw std::runtime_error{ss.str()};
			}
		}
	}




	std::map<std::string, std::unordered_map<std::string, std::vector<std::shared_ptr<AlignmentResults>>>> allAlnResults;
	std::map<std::string, std::unordered_map<std::string, std::vector<std::shared_ptr<AlignmentResults>>>> bestAlnResults;

	std::unordered_map<std::string, uint32_t> unmappedCounts;
	std::unordered_map<std::string, uint32_t> readLengths;
	std::unordered_map<std::string, std::string> nameKey;
	std::unordered_map<std::string, uint32_t> nameToPositionKey;
	auto tempSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(setUp.pars_.directoryName_, "temp_" + bfs::basename(setUp.pars_.ioOptions_.firstName_) ));

	{
		seqInfo seq;
		SeqInput reader(setUp.pars_.ioOptions_);
		SeqOutput tempWriter(tempSeqOpts);
		tempWriter.openOut();
		reader.openIn();
		uint32_t count = 0;

		while(reader.readNextRead(seq)){
			nameToPositionKey[seq.name_] = count;
			readLengths[seq.name_] = len(seq);
			nameKey[estd::to_string(count)] = seq.name_;
			seq.name_ = estd::to_string(count);
			tempWriter.write(seq);
			++count;
		}
	}
	auto tempInOpts = SeqIOOptions::genFastaIn(tempSeqOpts.out_.outName());
	tempInOpts.out_.transferOverwriteOpts(setUp.pars_.ioOptions_.out_);

	{
		table nameKeyTab(nameKey, VecStr{"key", "name"});
		OutOptions nameKeyOpts(njh::files::make_path(setUp.pars_.directoryName_, "nameKey.tab.txt"));
		OutputStream nameKeyOut(nameKeyOpts);
		nameKeyTab.outPutContents(nameKeyOut, "\t");
	}
	//align to genomes in parallel
	njh::concurrent::LockableQueue<std::string> genomesQueue(getVectorOfMapKeys(gMapper.genomes_));
	std::function<void()> alignToGenome = [&setUp,&genomesQueue,&gMapper,&bRunner,&tempInOpts](){
		std::string genome = "";
		while(genomesQueue.getVal(genome)){
			bfs::path genomeDir = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar{genome});
//			//align seqs
//			auto seqOpts = tempInOpts;
//			seqOpts.out_.outFilename_ = njh::files::make_path(genomeDir, "alignedSeqs.sorted.bam");
//
//			BioCmdsUtils::LastZPars copyLzPars = lzPars;
//			copyLzPars.genomeFnp = gMapper.genomes_.at(genome)->fnpTwoBit_;
//			auto runOut = bRunner.lastzAlign(seqOpts, copyLzPars);
//			OutOptions lastzAlignLogOpts(njh::files::make_path(genomeDir, "lastzLog.json"));
//			OutputStream lastzAlignLogOut(lastzAlignLogOpts);
//			lastzAlignLogOut << njh::json::toJson(runOut) << std::endl;


			//align seqs
			auto seqOpts = tempInOpts;
			seqOpts.out_.outFilename_ = njh::files::make_path(genomeDir, "alignedSeqs.sorted.bam");

			auto runOut = bRunner.bowtie2Align(seqOpts, gMapper.genomes_.at(genome)->fnp_, "-a");
			OutOptions bowtie2AlignLogOpts(njh::files::make_path(genomeDir, "bowtie2Log.json"));
			OutputStream bowtie2AlignLogOut(bowtie2AlignLogOpts);
			bowtie2AlignLogOut << njh::json::toJson(runOut) << std::endl;
		}
	};


	njh::concurrent::runVoidFunctionThreaded(alignToGenome, numThreads);

	for(const auto & genome : gMapper.genomes_){
		bfs::path genomeDir = njh::files::make_path(setUp.pars_.directoryName_, genome.first);
		//align seqs
		auto seqOpts = setUp.pars_.ioOptions_;
		seqOpts.out_.outFilename_ = njh::files::make_path(genomeDir, "alignedSeqs.sorted.bam");


		//extract locations and mapping stats
		BamTools::BamAlignment bAln;
		std::unordered_map<std::string, std::vector<BamTools::BamAlignment>> bamAligns;
		std::unordered_set<std::string> mappedReads;
		BamTools::BamReader bReader;
		bReader.Open(seqOpts.out_.outFilename_.string());
		checkBamOpenThrow(bReader, seqOpts.out_.outFilename_.string());
		auto refData = bReader.GetReferenceData();

		while (bReader.GetNextAlignment(bAln)) {
			if (bAln.IsMapped()) {
				bAln.Name = nameKey[bAln.Name];
				bamAligns[bAln.Name].emplace_back(bAln);
				mappedReads.emplace(bAln.Name);
			}
		}

		std::set<std::string> unmappedReads;
		for(const auto & readName : readLengths){
			if(!njh::in(readName.first, mappedReads)){
				unmappedReads.emplace(readName.first);
				++unmappedCounts[readName.first];
			}
		}

		TwoBit::TwoBitFile twobitReader(genome.second->fnpTwoBit_);
		std::vector<Bed6RecordCore> alignedRegions;
		std::map<uint32_t, uint32_t> mapCounts;
		mapCounts[0] = mappedReads.size();

		auto regionsExtractedOpts = SeqIOOptions::genFastaOut(njh::files::make_path(genomeDir, "regions"));
		SeqOutput extractedWriter(regionsExtractedOpts);
		extractedWriter.openOut();

		//alignment comparisons
		OutOptions comparisonOpts(njh::files::make_path(genomeDir, "refComparisonInfo.tab.txt"));
		OutputStream comparisonOut(comparisonOpts);
		comparisonOut << "ReadNumber\tReadId\tBestRef\tscore\talnScore\thqScore\tkDist-5\t1bIndel\t2bIndel\t>2bIndel\tlqMismatch\thqMismatch" << std::endl;
		//alignments
		auto alnOut = SeqIOOptions::genFastaOut(njh::files::make_path(genomeDir, "refAlignments.fasta"));
		SeqOutput writer(alnOut);
		writer.openOut();
		uint32_t readNumber = 0;
		std::vector<std::shared_ptr<seqInfo>> refSeqs;
		std::unordered_map<std::string, VecStr> readNamesToRefSeqs;
		auto bamAlignKeys = getVectorOfMapKeys(bamAligns);
		njh::sort(bamAlignKeys, [&nameToPositionKey](const std::string & name1, const std::string & name2){
			return nameToPositionKey[name1] < nameToPositionKey[name2];
		});
		for (const auto & bamAlignKey : bamAlignKeys) {
			const auto & alnForRead = bamAligns[bamAlignKey];
			++mapCounts[alnForRead.size()];
			uint32_t extractionCount = 0;
			for (const auto & aln : alnForRead) {

				auto results = std::make_shared<AlignmentResults>(aln, refData);
				results->setRefSeq(twobitReader);
				MetaDataInName refMeta;
				refMeta.addMeta("genome", genome.first);
				refMeta.addMeta("chrom", results->gRegion_.chrom_);
				refMeta.addMeta("start", results->gRegion_.start_);
				refMeta.addMeta("end", results->gRegion_.end_);
				results->refSeq_->name_.append(refMeta.createMetaName());
				if(!njh::in(results->refSeq_->name_, readNamesToRefSeqs)){
					refSeqs.emplace_back(std::make_shared<seqInfo>(*(results->refSeq_)));
				}
				readNamesToRefSeqs[results->refSeq_->name_].emplace_back(aln.Name);
				kmerInfo refInfo(results->refSeq_->seq_, 5, false);
				kmerInfo seqKInfo(results->alnSeq_->seq_, 5, false);

				//results->setComparison(false);
				results->setComparison(true);
				writer.write(results->refSeqAligned_);
				writer.write(results->alnSeqAligned_);

				alignedRegions.emplace_back(results->gRegion_.genBedRecordCore());
				std::string appName = "";
				MetaDataInName rangeMeta;
				if('H' == results->bAln_.CigarData.front().Type ){
					rangeMeta.addMeta("start",results->bAln_.CigarData.front().Length);
				}
				if('H' == results->bAln_.CigarData.back().Type ){
					rangeMeta.addMeta("end", readLengths[results->bAln_.Name] - results->bAln_.CigarData.back().Length);
				}
				if(!rangeMeta.meta_.empty()){
					appName = rangeMeta.createMetaName();
				}
				comparisonOut << readNumber
						<< '\t' << aln.Name << appName
						<< '\t' << results->gRegion_.createUidFromCoordsStrand()
						<< '\t' << results->comp_.distances_.eventBasedIdentityHq_
						<< '\t' << results->comp_.alnScore_
						<< '\t' << refInfo.compareKmers(seqKInfo).second
						<< '\t' << results->comp_.oneBaseIndel_
						<< '\t' << results->comp_.twoBaseIndel_
						<< '\t' << results->comp_.largeBaseIndel_
						<< '\t' << results->comp_.lqMismatches_
						<< '\t' << results->comp_.hqMismatches_ << std::endl;
				allAlnResults[aln.Name][genome.first].emplace_back(results);
				++extractionCount;
			}
			++readNumber;
		}
		//write out ref seqs;
		//read names for ref seqs
		OutOptions readNamesOpts(njh::files::make_path(genomeDir, "readNamesForRefSeqs.tab.txt"));
		OutputStream readNamesOut(readNamesOpts);
		readNamesOut << "refName\treadNames" << std::endl;
		VecStr refNames = njh::getVecOfMapKeys(readNamesToRefSeqs);
		njh::sort(refNames, [&readNamesToRefSeqs](const std::string & ref1, const std::string & ref2){
			return readNamesToRefSeqs[ref1].size() == readNamesToRefSeqs[ref2].size() ? ref1 < ref2 :  readNamesToRefSeqs[ref1].size() > readNamesToRefSeqs[ref2].size();
		});
		for(const auto & refName : refNames){
			readNamesOut << refName
					<< "\t" << njh::conToStr(readNamesToRefSeqs[refName], ",") << std::endl;
		}
		//collapse refseqs
		for(auto & refSeq : refSeqs){
			refSeq->cnt_ = readNamesToRefSeqs[refSeq->name_].size();
			MetaDataInName refMeta(refSeq->name_);
			refMeta.addMeta("extractCount", refSeq->cnt_);
			refMeta.resetMetaInName(refSeq->name_);
		}
		njh::sort(refSeqs, [](const std::shared_ptr<seqInfo> & ref1, const std::shared_ptr<seqInfo> & ref2){
			return ref1->cnt_ == ref2->cnt_ ? ref1->name_ < ref2->name_ : ref1->cnt_ > ref2->cnt_;
		});
		extractedWriter.write(refSeqs);


		//genomic regions hit
		OutOptions regionsOpts(njh::files::make_path(genomeDir, "regions.bed"));
		OutputStream regionsOut(regionsOpts);

		BedUtility::coordSort(alignedRegions);
		for(const auto & reg : alignedRegions){
			regionsOut << reg.toDelimStrWithExtra() << std::endl;
		}

		//map counts
		OutOptions mapCountsOpts(njh::files::make_path(genomeDir, "mapCounts.tab.txt"));
		OutputStream mapCountsOut(mapCountsOpts);
		table mapCountsTab(mapCounts, VecStr{"hits", "total"});
		mapCountsTab.outPutContents(mapCountsOut, "\t");

		//names of unmapped sequences
		OutOptions unmmapedOpts(njh::files::make_path(genomeDir, "unmappedReads.txt"));
		OutputStream unmappedOut(unmmapedOpts);
		for(const auto & unmappedAln : unmappedReads){
			unmappedOut << unmappedAln << std::endl;
		}
	}


	//get best hits only



	auto regionsExtractedOpts = SeqIOOptions::genFastaOut(njh::files::make_path(setUp.pars_.directoryName_, "bestRegions"));
	SeqOutput extractedWriter(regionsExtractedOpts);
	extractedWriter.openOut();
	std::vector<std::shared_ptr<seqInfo>> refSeqs;
	std::unordered_map<std::string, VecStr> readNamesToRefSeqs;

	//alignment comparisons
	OutOptions comparisonOpts(njh::files::make_path(setUp.pars_.directoryName_, "refComparisonInfo.tab.txt"));
	OutputStream comparisonOut(comparisonOpts);

	comparisonOut << "ReadNumber\tReadId\tBestRef\tscore\talnScore\tkDist-5\t1bIndel\t2bIndel\t>2bIndel\tlqMismatch\thqMismatch" << std::endl;
	uint32_t readNumber = 0;
	std::unordered_map<std::string, std::vector<Bed6RecordCore>> bestRegionsByGenome;

	auto allAlnResultsKeys = getVectorOfMapKeys(allAlnResults);
	njh::sort(allAlnResultsKeys, [&nameToPositionKey](const std::string & name1, const std::string & name2){
		return nameToPositionKey[name1] < nameToPositionKey[name2];
	});

	for(const auto & allAlnResultsKey : allAlnResultsKeys){
		const auto & alnResults = allAlnResults[allAlnResultsKey];
		double bestScore = std::numeric_limits<double>::lowest();
		std::vector<std::shared_ptr<AlignmentResults>> bestResults;
		std::unordered_map<std::string, std::string> regionNameToGenome;

		for(const auto & genomeRes : alnResults){
			for(const auto & res : genomeRes.second){
				if(res->comp_.alnScore_ > bestScore){
					bestScore = res->comp_.alnScore_;
					bestResults.clear();
					bestResults.emplace_back(res);
					regionNameToGenome[res->gRegion_.genBedRecordCore().toDelimStrWithExtra()] = genomeRes.first;
				}else if(res->comp_.alnScore_  == bestScore){
					bestResults.emplace_back(res);
					regionNameToGenome[res->gRegion_.genBedRecordCore().toDelimStrWithExtra()] = genomeRes.first;
				}
			}
		}

		for(const auto & results : bestResults){
			bestRegionsByGenome[regionNameToGenome[results->gRegion_.genBedRecordCore().toDelimStrWithExtra()]].emplace_back(results->gRegion_.genBedRecordCore());

			bestAlnResults[allAlnResultsKey][regionNameToGenome[results->gRegion_.genBedRecordCore().toDelimStrWithExtra()]].emplace_back(results);
		}

		for(const auto & results : bestResults){
			if(!njh::in(results->refSeq_->name_, readNamesToRefSeqs)){
				refSeqs.emplace_back(results->refSeq_);
			}

			readNamesToRefSeqs[results->refSeq_->name_].emplace_back(results->bAln_.Name);

			kmerInfo refInfo(results->refSeq_->seq_, 5, false);
			kmerInfo seqKInfo(results->alnSeq_->seq_, 5, false);
			std::string appName = "";
			MetaDataInName rangeMeta;
			if('H' == results->bAln_.CigarData.front().Type ){
				rangeMeta.addMeta("start",results->bAln_.CigarData.front().Length);
			}
			if('H' == results->bAln_.CigarData.back().Type ){
				rangeMeta.addMeta("end", readLengths[results->bAln_.Name] - results->bAln_.CigarData.back().Length);
			}
			if(!rangeMeta.meta_.empty()){
				appName = rangeMeta.createMetaName();
			}
			comparisonOut << readNumber
					<< '\t' << results->bAln_.Name << appName
					<< '\t' << results->gRegion_.createUidFromCoordsStrand()
					<< '\t' << results->comp_.distances_.eventBasedIdentityHq_
					<< '\t' << results->comp_.alnScore_
					<< '\t' << refInfo.compareKmers(seqKInfo).second
					<< '\t' << results->comp_.oneBaseIndel_
					<< '\t' << results->comp_.twoBaseIndel_
					<< '\t' << results->comp_.largeBaseIndel_
					<< '\t' << results->comp_.lqMismatches_
					<< '\t' << results->comp_.hqMismatches_ << std::endl;
		}
		++readNumber;
	}



	for(auto & best : bestRegionsByGenome){
		OutOptions bestRegionsBedOpts(njh::files::make_path(setUp.pars_.directoryName_, njh::pasteAsStr("bestRegions_", best.first, ".bed")));
		OutputStream bestRegionsBedOut(bestRegionsBedOpts);
		if(bfs::exists(gMapper.genomes_.at(best.first)->gffFnp_)){
			intersectBedLocsWtihGffRecordsPars gffPars = gMapper.pars_.gffIntersectPars_;
			gffPars.gffFnp_ = gMapper.genomes_.at(best.first)->gffFnp_;
			intersectBedLocsWtihGffRecords(best.second, gffPars);
		}
		BedUtility::coordSort(best.second);
		for(const auto & reg : best.second){
			bestRegionsBedOut << reg.toDelimStrWithExtra() << std::endl;
		}
	}





	//write out ref seqs;
	//read names for ref seqs
	OutOptions readNamesOpts(
			njh::files::make_path(setUp.pars_.directoryName_, "readNamesForRefSeqs.tab.txt"));
	OutputStream readNamesOut(readNamesOpts);
	readNamesOut << "refName\treadNames" << std::endl;
	VecStr refNames = njh::getVecOfMapKeys(readNamesToRefSeqs);
	njh::sort(refNames,
			[&readNamesToRefSeqs](const std::string & ref1, const std::string & ref2) {
				return readNamesToRefSeqs[ref1].size() == readNamesToRefSeqs[ref2].size() ? ref1 < ref2 : readNamesToRefSeqs[ref1].size() > readNamesToRefSeqs[ref2].size();
			});
	for (const auto & refName : refNames) {
		readNamesOut << refName << "\t"
				<< njh::conToStr(readNamesToRefSeqs[refName], ",") << std::endl;
	}
	//collapse refseqs
	for (auto & refSeq : refSeqs) {
		refSeq->cnt_ = readNamesToRefSeqs[refSeq->name_].size();
		MetaDataInName refMeta(refSeq->name_);
		refMeta.addMeta("extractCount", refSeq->cnt_, true);
		refMeta.resetMetaInName(refSeq->name_);
	}
	njh::sort(refSeqs,
			[](const std::shared_ptr<seqInfo> & ref1, const std::shared_ptr<seqInfo> & ref2) {
				return ref1->cnt_ == ref2->cnt_ ? ref1->name_ < ref2->name_ : ref1->cnt_ > ref2->cnt_;
			});
	extractedWriter.write(refSeqs);

	//names of unmapped sequences
	VecStr unmappedToAllGenomes;
	for(const auto & unmappedCount : unmappedCounts){
		if(gMapper.genomes_.size() == unmappedCount.second ){
			unmappedToAllGenomes.emplace_back(unmappedCount.first);
			comparisonOut << readNumber
					<< '\t' << unmappedCount.first
					<< '\t' << "*"
					<< '\t' << "*"
					<< '\t' << "*"
					<< '\t' << "*"
					<< '\t' << "*"
					<< '\t' << "*"
					<< '\t' << "*"
					<< '\t' << "*"
					<< '\t' << "*" << std::endl;
			++readNumber;
		}
	}
	if(!unmappedToAllGenomes.empty()){
		OutOptions unmmapedOpts(njh::files::make_path(setUp.pars_.directoryName_, "unmappedReads.txt"));
		OutputStream unmappedOut(unmmapedOpts);
		for(const auto & unmappedAln : unmappedToAllGenomes){
			unmappedOut << unmappedAln << std::endl;
		}
	}
	bfs::remove(tempSeqOpts.out_.outName());



	{
		// determine regions covered
		std::unordered_map<std::string, std::map<uint32_t, uint32_t>> simpleCoverageCounts;
		for(const auto & seqAlnResults : bestAlnResults){
			for(const auto & genome : seqAlnResults.second){
				for(const auto & res : genome.second){
					if(amountOfErrorForCoverageCalc.passErrorProfile(res->comp_)){
						for(uint32_t pos = res->gRegion_.start_; pos <res->gRegion_.end_; ++pos){
							simpleCoverageCounts[res->gRegion_.chrom_][pos] += 1;
						}
					}
				}
			}
		}

		//getting regions that were covered but not expected
		{
			std::unordered_map<std::string, std::map<uint32_t, uint32_t>> requiredRegionsPositions;
			for(const auto & reg : requiredRegions){
				for(uint32_t pos = reg.start_; pos < reg.end_; ++pos){
					requiredRegionsPositions[reg.chrom_][pos] +=1;
				}
			}
			std::vector<Bed3RecordCore> allRegionsNotExpected;

			for(const auto & chrom : simpleCoverageCounts){
				std::vector<Bed3RecordCore> regionsNotExpectedRaw;
				for(const auto & pos : chrom.second){

					if(0 == requiredRegionsPositions[chrom.first][pos.first] ){
						regionsNotExpectedRaw.emplace_back(Bed3RecordCore(chrom.first, pos.first, pos.first + 1));
					}
				}
				BedUtility::coordSort(regionsNotExpectedRaw, false);
				std::vector<Bed3RecordCore> regionsNotExpected;
				for(const auto & region : regionsNotExpectedRaw){
					if(regionsNotExpected.empty()){
						regionsNotExpected.emplace_back(region);
					}else{
						if(regionsNotExpected.back().chromEnd_ == region.chromStart_){
							regionsNotExpected.back().chromEnd_ = region.chromEnd_;
						}else{
							regionsNotExpected.emplace_back(region);
						}
					}
				}
				addOtherVec(allRegionsNotExpected, regionsNotExpected);
			}
			std::vector<std::shared_ptr<Bed6RecordCore>> allRegionsNotExpectedB6;

			for(const auto & notExp : allRegionsNotExpected){
				allRegionsNotExpectedB6.emplace_back(std::make_shared<Bed6RecordCore>(GenomicRegion(notExp).genBedRecordCore()));
			}
			std::unordered_map<std::string, std::vector<std::shared_ptr<Bed6RecordCore>>> allRegionsByGenome;

			for(const auto & region : allRegionsNotExpectedB6){
				allRegionsByGenome[chromosomeToGenome[region->chrom_]].push_back(region);
			}
			for( auto & genome : allRegionsByGenome){
				if("" != gMapper.genomes_.at(genome.first)->gffFnp_){
					intersectBedLocsWtihGffRecordsPars interPars;
					interPars.extraAttributes_ = mapperPars.gffIntersectPars_.extraAttributes_;
					interPars.gffFnp_ = gMapper.genomes_.at(genome.first)->gffFnp_;
					intersectBedLocsWtihGffRecords(genome.second, interPars	);
				}
			}
			BedUtility::coordSort(allRegionsNotExpectedB6);
			if(!allRegionsNotExpectedB6.empty()){
				OutputStream out(njh::files::make_path(setUp.pars_.directoryName_, "regionsNotExpected.bed"));
				for(const auto & region : allRegionsNotExpectedB6){
					out << region->toDelimStrWithExtra() << std::endl;
				}
			}
		}
		//

		table coveredCountsTab(VecStr{"Region", "basesCovered", "totalBases", "fractionCovered"});
		for(const auto & gene: requiredRegions){
			uint32_t covered = 0;
			std::vector<uint32_t> positionsNotCovered;
			for(const auto & pos : iter::range(gene.start_, gene.end_)){
				if(simpleCoverageCounts[gene.chrom_][pos] > 0){
					++covered;
				}else{
					positionsNotCovered.emplace_back(pos);
				}
			}
			coveredCountsTab.addRow(gene.uid_, covered, gene.getLen(), static_cast<double>(covered)/gene.getLen());

			if(positionsNotCovered.size() > 0){
				std::vector<Bed6RecordCore> regionsNotCoveredRaw;
				for(const auto & pos : positionsNotCovered){
					regionsNotCoveredRaw.emplace_back(Bed6RecordCore(gene.chrom_, pos, pos + 1, gene.uid_, 1, gene.reverseSrand_? '-':'+'));
				}
				BedUtility::coordSort(regionsNotCoveredRaw, false);
				std::vector<Bed6RecordCore> regionsNotCovered;
				for(const auto & region : regionsNotCoveredRaw){
					if(regionsNotCovered.empty()){
						regionsNotCovered.emplace_back(region);
					}else{
						if(regionsNotCovered.back().chromEnd_ == region.chromStart_){
							regionsNotCovered.back().chromEnd_ = region.chromEnd_;
						}else{
							regionsNotCovered.emplace_back(region);
						}
					}
				}
				for( auto & region : regionsNotCovered){
					region.score_ = region.length();
					region.name_ = njh::pasteAsStr(gene.uid_, ":", GenomicRegion(region).createUidFromCoords());
				}
				OutputStream out(njh::files::make_path(setUp.pars_.directoryName_, gene.uid_ + "_notCoveredRegions.bed"));
				for(const auto & region : regionsNotCovered){
					out << region.toDelimStr() << std::endl;
				}
			}
		}



		coveredCountsTab.outPutContents(
				TableIOOpts::genTabFileOut(
						njh::files::make_path(setUp.pars_.directoryName_,
								"coveragedInfo.tab.txt"), true));


	}
	return 0;
}



}  // namespace njhseq

