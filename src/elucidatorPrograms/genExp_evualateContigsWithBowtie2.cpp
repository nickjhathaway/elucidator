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
#include <njhseq/objects/seqContainers/refVariants.hpp>
#include <utility>


namespace njhseq {





int genExpRunner::extractFromGenomesAndCompare(const njh::progutils::CmdArgs & inputCommands){

	std::string genomesStr = "";
	uint32_t numThreads = 1;
	std::string program = "";
	std::string sample  = "";
	std::string target = "";
	uint32_t lencutOffForBowtie = 200;
	BioCmdsUtils::LastZPars lzPars;
	lzPars.coverage = 95;
	lzPars.identity = 90;

	uint32_t nBowtie2Alns = 100;

	comparison amountOfErrorForCoverageCalc;
	MultiGenomeMapper::inputParameters mapperPars;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processComparison(amountOfErrorForCoverageCalc);
	setUp.setOption(lzPars.coverage, "--coverage", "coverage for lastz");
	setUp.setOption(lzPars.identity, "--identity", "identity for lastz");

  setUp.setOption(target, "--target", "Name of the target region to output with the report", true);
  setUp.setOption(program, "--program", "Name of the program to output with the report", true);
  setUp.setOption(sample,  "--sample",  "Name of the sample to output with the report",  true);
	setUp.setOption(lencutOffForBowtie, "--lencutOffForBowtie", "Contigs below this length will be evaulated by bowtie instead of lastz");
	setUp.setOption(nBowtie2Alns, "--nBowtie2Alns", "The max number of bowtie alignments to check");

	setUp.setOption(numThreads, "--numThreads", "number of cpus to use");
  setUp.setOption(genomesStr, "--genomes", "Names of the genomes to extract from (should not have extension, e.g. Pf3d7 for Pf3d7.fasta");
  setUp.setOption(mapperPars.genomeDir_, "--genomeDir", "Names of the genome directory where the genomes are stored", true);
  setUp.setOption(mapperPars.gffDir_, "--gffDir", "A directory of gffs that go with the genomes");
  setUp.setOption(mapperPars.gffIntersectPars_.extraAttributes_, "--gffExtraAttributes", "gff Extra Attributes");

  setUp.processReadInNames();
  setUp.processDirectoryOutputName(true);
	setUp.finishSetUp(std::cout);

	setUp.startARunLog(setUp.pars_.directoryName_);

	njh::sys::requireExternalProgramsThrow(VecStr{"bowtie2", "lastz", "samtools"});

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



	std::map<std::string, std::unordered_map<std::string, std::vector<std::shared_ptr<AlignmentResults>>>> allAlnResults;
	std::map<std::string, std::unordered_map<std::string, std::vector<std::shared_ptr<AlignmentResults>>>> bestAlnResults;

	std::unordered_map<std::string, uint32_t> unmappedCounts;
	std::unordered_map<std::string, uint32_t> readLengths;
	std::unordered_map<std::string, std::string> nameKey;
	std::unordered_map<std::string, uint32_t> nameToPositionKey;
	auto tempSeqBowtie2Opts = SeqIOOptions::genFastaOut(njh::files::make_path(setUp.pars_.directoryName_, "tempForBowtie2_" + bfs::basename(setUp.pars_.ioOptions_.firstName_) ));
	auto tempSeqLastzOpts = SeqIOOptions::genFastaOut(njh::files::make_path(setUp.pars_.directoryName_, "tempForLastz_" + bfs::basename(setUp.pars_.ioOptions_.firstName_) ));

	{
		seqInfo seq;
		SeqInput reader(setUp.pars_.ioOptions_);
		SeqOutput tempForBowtie2Writer(tempSeqBowtie2Opts);

		SeqOutput tempForLastzWriter(tempSeqLastzOpts);

		reader.openIn();
		uint32_t count = 0;

		while(reader.readNextRead(seq)){
			nameToPositionKey[seq.name_] = count;
			readLengths[seq.name_] = len(seq);
			nameKey[estd::to_string(count)] = seq.name_;
			seq.name_ = estd::to_string(count);
			if (len(seq) <= lencutOffForBowtie) {
				tempForBowtie2Writer.openWrite(seq);
			} else {
				tempForLastzWriter.openWrite(seq);
			}
			++count;
		}
	}
	auto tempInBowtie2Opts = SeqIOOptions::genFastaIn(tempSeqBowtie2Opts.out_.outName());
	tempInBowtie2Opts.out_.transferOverwriteOpts(setUp.pars_.ioOptions_.out_);

	auto tempLastzInOpts = SeqIOOptions::genFastaIn(tempSeqLastzOpts.out_.outName());
	tempLastzInOpts.out_.transferOverwriteOpts(setUp.pars_.ioOptions_.out_);

	{
		table nameKeyTab(nameKey, VecStr{"key", "name"});
		OutOptions nameKeyOpts(njh::files::make_path(setUp.pars_.directoryName_, "nameKey.tab.txt"));
		OutputStream nameKeyOut(nameKeyOpts);
		nameKeyTab.outPutContents(nameKeyOut, "\t");
	}
	//align to genomes in parallel

	struct GenomeWithProgram {
		GenomeWithProgram() = default;
		GenomeWithProgram(bool lastz, std::string  genome) :
				lastz_(lastz), genome_(std::move(genome)) {
		}
		bool lastz_ { false };
		std::string genome_;
	};
	std::vector<GenomeWithProgram> genomesWithPrograms;
	auto inputGenomes = getVectorOfMapKeys(gMapper.genomes_);
	for(const auto & genome : inputGenomes){
		genomesWithPrograms.emplace_back(false, genome);
		genomesWithPrograms.emplace_back(true, genome);
		njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar{genome});
	}
	njh::concurrent::LockableQueue<GenomeWithProgram> genomesQueue(genomesWithPrograms);
	std::function<void()> alignToGenome = [&setUp,&genomesQueue,&gMapper,&bRunner,&tempLastzInOpts,&tempInBowtie2Opts,&lzPars,&nBowtie2Alns](){
		GenomeWithProgram gAndProgram;
		while(genomesQueue.getVal(gAndProgram)){
			bfs::path genomeDir = njh::files::make_path(setUp.pars_.directoryName_, gAndProgram.genome_);
			if(gAndProgram.lastz_){
				if(tempLastzInOpts.inExists()){
					//align seqs
					auto seqOpts = tempLastzInOpts;
					seqOpts.out_.outFilename_ = njh::files::make_path(genomeDir, "alignedSeqsLastz.sorted.bam");

					BioCmdsUtils::LastZPars copyLzPars = lzPars;
					copyLzPars.genomeFnp = gMapper.genomes_.at(gAndProgram.genome_)->fnpTwoBit_;

					auto runOut = bRunner.lastzAlign(seqOpts, copyLzPars);
					OutOptions lastzAlignLogOpts(njh::files::make_path(genomeDir, "lastzLog.json"));
					OutputStream lastzAlignLogOut(lastzAlignLogOpts);
					lastzAlignLogOut << njh::json::toJson(runOut) << std::endl;
				}
			} else {
				//align seqs
				if(tempInBowtie2Opts.inExists()){
					auto seqOpts = tempInBowtie2Opts;
					seqOpts.out_.outFilename_ = njh::files::make_path(genomeDir, "alignedSeqsBowtie2.sorted.bam");

					auto runOut = bRunner.bowtie2Align(seqOpts, gMapper.genomes_.at(gAndProgram.genome_)->fnp_, njh::pasteAsStr("-k ", nBowtie2Alns));
					OutOptions bowtie2AlignLogOpts(njh::files::make_path(genomeDir, "bowtie2Log.json"));
					OutputStream bowtie2AlignLogOut(bowtie2AlignLogOpts);
					bowtie2AlignLogOut << njh::json::toJson(runOut) << std::endl;
				}
			}
		}
	};


	njh::concurrent::runVoidFunctionThreaded(alignToGenome, numThreads);

	for(const auto & genome : gMapper.genomes_){
		bfs::path genomeDir = njh::files::make_path(setUp.pars_.directoryName_, genome.first);
		auto bowtie2AlignedFnp = njh::files::make_path(genomeDir, "alignedSeqsBowtie2.sorted.bam");
		auto lastzAlignedFnp = njh::files::make_path(genomeDir, "alignedSeqsLastz.sorted.bam");
		std::unordered_map<std::string, std::vector<BamTools::BamAlignment>> bamAligns;
		std::unordered_set<std::string> mappedReads;
		BamTools::RefVector refData;
		if(bfs::exists(bowtie2AlignedFnp)){
			auto seqOpts = setUp.pars_.ioOptions_;
			seqOpts.out_.outFilename_ = bowtie2AlignedFnp;
			//extract locations and mapping stats
			BamTools::BamAlignment bAln;
			BamTools::BamReader bReader;
			bReader.Open(seqOpts.out_.outFilename_.string());
			checkBamOpenThrow(bReader, seqOpts.out_.outFilename_.string());
			refData = bReader.GetReferenceData();

			while (bReader.GetNextAlignment(bAln)) {
				if (bAln.IsMapped()) {
					bAln.Name = nameKey[bAln.Name];
					bamAligns[bAln.Name].emplace_back(bAln);
					mappedReads.emplace(bAln.Name);
				}
			}
		}
		if(bfs::exists(lastzAlignedFnp)){
			auto seqOpts = setUp.pars_.ioOptions_;
			seqOpts.out_.outFilename_ = lastzAlignedFnp;
			//extract locations and mapping stats
			BamTools::BamAlignment bAln;
			BamTools::BamReader bReader;
			bReader.Open(seqOpts.out_.outFilename_.string());
			checkBamOpenThrow(bReader, seqOpts.out_.outFilename_.string());
			refData = bReader.GetReferenceData();

			while (bReader.GetNextAlignment(bAln)) {
				if (bAln.IsMapped()) {
					bAln.Name = nameKey[bAln.Name];
					bamAligns[bAln.Name].emplace_back(bAln);
					mappedReads.emplace(bAln.Name);
				}
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
						<< '\t' << results->comp_.distances_.eventBasedIdentity_
						<< '\t' << results->comp_.alnScore_
						<< '\t' << results->comp_.distances_.eventBasedIdentityHq_
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

	comparisonOut << "ReadNumber\tReadId\tBestRef\tscore\talnScore\tkDist-5\t1bIndel\t2bIndel\t>2bIndel\tlqMismatch\thqMismatch"
			<< '\t' << "program"
			<< '\t' << "sample"
			<< '\t' << "target"
			<< std::endl;




	uint32_t readNumber = 0;
	std::unordered_map<std::string, std::vector<Bed6RecordCore>> bestRegionsByGenome;

	auto allAlnResultsKeys = getVectorOfMapKeys(allAlnResults);
	njh::sort(allAlnResultsKeys, [&nameToPositionKey](const std::string & name1, const std::string & name2){
		return nameToPositionKey[name1] < nameToPositionKey[name2];
	});


	auto refAlnOutErrorOpts = SeqIOOptions::genFastaOut(njh::files::make_path(setUp.pars_.directoryName_, "refAlignmentsWithErrors.fasta"));
	SeqOutput refAlnOutErrorOut(refAlnOutErrorOpts);
	refAlnOutErrorOut.openOut();
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
			if (1 != results->comp_.distances_.eventBasedIdentityHq_) {
				refAlnOutErrorOut.write(results->refSeqAligned_);
				refAlnOutErrorOut.write(results->alnSeqAligned_);
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
					<< '\t' << results->comp_.hqMismatches_
					<< '\t' << program
					<< '\t' << sample
					<< '\t' << target
					<< std::endl;
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
					<< '\t' << "*"
					<< '\t' << program
					<< '\t' << sample
					<< '\t' << target
					<< std::endl;
			/*
			 * 		coveredCountsTab.addColumn(VecStr{program}, "program");
		coveredCountsTab.addColumn(VecStr{sample}, "sample");
			 */
			++readNumber;
		}
	}

	std::set<std::string> matchingContigs;
	std::set<std::string> notMatchingContigs(unmappedToAllGenomes.begin(), unmappedToAllGenomes.end());

	if(!unmappedToAllGenomes.empty()){
		OutOptions unmmapedOpts(njh::files::make_path(setUp.pars_.directoryName_, "unmappedReads.txt"));
		OutputStream unmappedOut(unmmapedOpts);
		for(const auto & unmappedAln : unmappedToAllGenomes){
			unmappedOut << unmappedAln << std::endl;
		}
	}
	bfs::remove(tempSeqBowtie2Opts.out_.outName());
	bfs::remove(tempSeqLastzOpts.out_.outName());



	{

		// determine regions covered
		std::unordered_map<std::string, std::map<uint32_t, uint32_t>> simpleCoverageCounts;
		for(const auto & seqAlnResults : bestAlnResults){
			for(const auto & genome : seqAlnResults.second){
				for(const auto & res : genome.second){
					if(amountOfErrorForCoverageCalc.passErrorProfile(res->comp_)){
						matchingContigs.emplace(res->bAln_.Name);
						for(uint32_t pos = res->gRegion_.start_; pos <res->gRegion_.end_; ++pos){
							simpleCoverageCounts[res->gRegion_.chrom_][pos] += 1;
						}
					}else{
						notMatchingContigs.emplace(res->bAln_.Name);
					}
				}
			}
		}

		{
			//write out matching contigs info
			OutputStream matchInfo(njh::files::make_path(setUp.pars_.directoryName_, "seqsMatchingExpectedInfo.tab.txt"));
			matchInfo << "program\tsample\ttarget";
			matchInfo
					<< "\t" << "SeqsMatchingExpectedCnt"
					<< "\t" << "SeqsMatchingExpectedFrac"
					<< "\t" << "SeqsMatchingNotExpectedCnt"
					<< "\t" << "SeqsMatchingNotExpectedFrac"
					<< "\t" << "TotalSeqs" << std::endl;
			matchInfo << program
					<< "\t" << sample
					<< "\t" << target
					<< "\t" << matchingContigs.size()
					<< "\t" << matchingContigs.size()/static_cast<double>(readNumber)
					<< "\t" << notMatchingContigs.size()
					<< "\t" << notMatchingContigs.size()/static_cast<double>(readNumber)
					<< "\t" << readNumber
					<< std::endl;

		}

		//getting regions that were covered but not expected
		{
			std::vector<Bed3RecordCore> allRegionsNotExpected;

			for(const auto & chrom : simpleCoverageCounts){
				std::vector<Bed3RecordCore> regionsNotExpectedRaw;
				for(const auto & pos : chrom.second){
					regionsNotExpectedRaw.emplace_back(Bed3RecordCore(chrom.first, pos.first, pos.first + 1));
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
					intersectBedLocsWtihGffRecordsPars interPars = gMapper.pars_.gffIntersectPars_;
					interPars.gffFnp_ = gMapper.genomes_.at(genome.first)->gffFnp_;
					intersectBedLocsWtihGffRecords(genome.second, interPars	);
				}
			}
			BedUtility::coordSort(allRegionsNotExpectedB6);
			if(!allRegionsNotExpectedB6.empty()){
				OutputStream out(njh::files::make_path(setUp.pars_.directoryName_, "regionsCovered.bed"));
				for(const auto & region : allRegionsNotExpectedB6){
					out << region->toDelimStrWithExtra() << std::endl;
				}
			}
		}
	}

	{
		//read lengths
		auto allRLens = getVectorOfMapValues(readLengths);
		auto stats = getStatsOnVec(allRLens);
		table statsTable;
		statsTable.columnNames_ = getVectorOfMapKeys(stats);
		statsTable.hasHeader_ = true;
		auto nums = getVectorOfMapValues(stats);
		statsTable.content_.emplace_back(numVecToVecStr(nums));
		addOtherVec(statsTable.columnNames_, VecStr{"n","program", "sample", "target"});
		addOtherVec(statsTable.content_.front(), toVecStr(allRLens.size(), program, sample, target));
		statsTable.outPutContents(
				TableIOOpts::genTabFileOut(
						njh::files::make_path(setUp.pars_.directoryName_,
								"contigsLengthsInfo.tab.txt"), true));
	}
	return 0;
}




int genExpRunner::evaluateContigsAgainstExpected(const njh::progutils::CmdArgs & inputCommands){
	std::string genomesStr;
	uint32_t numThreads = 1;
	bfs::path bedFnp = "";
	bfs::path idsRequired = "";
	std::string program;
	std::string sample;
	std::string target;
	uint32_t lencutOffForBowtie = 200;
	BioCmdsUtils::LastZPars lzPars;
	lzPars.coverage = 95;
	lzPars.identity = 90;
	//VecStr features{"genes","pseudogene"};
	uint32_t nBowtie2Alns = 100;
	aligner alignerObj(10, gapScoringParameters(5,1), substituteMatrix::createDegenScoreMatrix(1,-1));
	alignerObj.weighHomopolymers_ = false;
	alignerObj.countEndGaps_ = false;

	bool calcSpecificCoverage = false;
	comparison amountOfErrorForCoverageCalc;
	MultiGenomeMapper::inputParameters mapperPars;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processComparison(amountOfErrorForCoverageCalc);
	setUp.setOption(lzPars.coverage, "--coverage", "coverage for lastz");
	setUp.setOption(lzPars.identity, "--identity", "identity for lastz");
	setUp.setOption(calcSpecificCoverage, "--calcSpecificCoverage", "calculate coverage of regions specific to each input reference");
	//setUp.setOption(features, "--features", "features to look for");
	setUp.setOption(alignerObj.weighHomopolymers_, "--weighHomopolymers", "weigh Homopolymers");

  setUp.setOption(program, "--program", "Name of the program to output with the report", true);
  setUp.setOption(sample,  "--sample",  "Name of the sample to output with the report",  true);
  setUp.setOption(target,  "--target", "Name of the target region to output with the report", true);

	setUp.setOption(lencutOffForBowtie, "--lencutOffForBowtie", "Contigs below this length will be evaulated by bowtie instead of lastz");
	setUp.setOption(nBowtie2Alns, "--nBowtie2Alns", "The max number of bowtie alignments to check");

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

	njh::sys::requireExternalProgramsThrow(VecStr{"bowtie2", "lastz", "samtools"});

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
					//if("gene" == gff.type_&& gff.hasAttr("ID") && njh::in(gff.getAttr("ID"), ids)){
					if(gff.hasAttr("ID") && njh::in(gff.getAttr("ID"), ids)){
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


	double sumTotalRequired = 0;
	for(const auto & region : requiredRegions){
		sumTotalRequired += region.getLen();
	}

	std::vector<seqInfo> requiredRegionsSeqs;
	for(const auto  & region : requiredRegions){
		TwoBit::TwoBitFile treader(gMapper.genomes_.at(chromosomeToGenome[region.chrom_])->fnpTwoBit_);
		auto refSeq = region.extractSeq(treader);
		if(region.reverseSrand_){
			refSeq.reverseComplementRead(false, true);
		}
		requiredRegionsSeqs.emplace_back(refSeq);
	}

	{
		auto expectedRegionsDir = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar("expectedRegionsInputInfo"));
		SeqOutput::write(requiredRegionsSeqs, SeqIOOptions::genFastaOutGz(njh::files::make_path(expectedRegionsDir, "expectedRegionsSeqs.fasta.gz")));
		OutputStream expectedRegionsOut(njh::files::make_path(expectedRegionsDir, "expectedRegions.bed"));
		for(const auto & region : requiredRegions){
			expectedRegionsOut << region.genBedRecordCore().toDelimStrWithExtra() << std::endl;
		}
	}

	std::unordered_map<std::string, std::vector<GenomicRegion>> specificRegions;
	if(calcSpecificCoverage){
		if(1 == requiredRegions.size()){
			specificRegions[requiredRegions[0].uid_] = requiredRegions;
		}else{
			uint64_t maxLen = 0;
			readVec::getMaxLength(requiredRegionsSeqs, maxLen);
			std::vector<refVariants> refVariationInfo;
			aligner alignerObjForComps(maxLen, gapScoringParameters(5, 1, 0, 0, 0, 0), substituteMatrix(2, -2), false);
			std::unordered_map<std::string, GenomicRegion> regionByName;


			for (const auto refPos : iter::range(requiredRegionsSeqs.size())) {
				regionByName[requiredRegions[refPos].uid_] = requiredRegions[refPos];
				refVariationInfo.emplace_back(requiredRegionsSeqs[refPos]);
			}
			for (const auto refPos : iter::range(requiredRegionsSeqs.size())) {
				for (const auto refSubPos : iter::range(requiredRegionsSeqs.size())) {
					if (refPos == refSubPos) {
						continue;
					}
					refVariationInfo[refPos].addVariant(requiredRegionsSeqs[refSubPos],
																							alignerObjForComps, false);
				}
			}
			for (const auto refPos : iter::range(requiredRegionsSeqs.size())) {
				auto specificPositions = refVariationInfo[refPos].getUniqueToRefPositions();
				for(const auto & pos : specificPositions){
					auto specReg = regionByName[requiredRegions[refPos].uid_];
					specReg.start_ = specReg.start_ + pos;
					specReg.end_ = specReg.start_ + 1;
					specificRegions[requiredRegions[refPos].uid_].emplace_back(specReg);
				}
			}
		}
	}


	std::map<std::string, std::unordered_map<std::string, std::vector<std::shared_ptr<AlignmentResults>>>> allAlnResults;
	std::map<std::string, std::unordered_map<std::string, std::vector<std::shared_ptr<AlignmentResults>>>> bestAlnResults;

	std::unordered_map<std::string, uint32_t> unmappedCounts;
	std::unordered_map<std::string, uint32_t> readLengths;
	std::unordered_map<std::string, std::string> nameKey;
	std::unordered_map<std::string, uint32_t> nameToPositionKey;
	auto tempSeqBowtie2Opts = SeqIOOptions::genFastaOut(njh::files::make_path(setUp.pars_.directoryName_, "tempForBowtie2_" + bfs::basename(setUp.pars_.ioOptions_.firstName_) ));
	auto tempSeqLastzOpts = SeqIOOptions::genFastaOut(njh::files::make_path(setUp.pars_.directoryName_, "tempForLastz_" + bfs::basename(setUp.pars_.ioOptions_.firstName_) ));

	std::vector<uint32_t> allContigsReadLengths;
	{
		seqInfo seq;
		SeqInput reader(setUp.pars_.ioOptions_);
		SeqOutput tempForBowtie2Writer(tempSeqBowtie2Opts);

		SeqOutput tempForLastzWriter(tempSeqLastzOpts);

		reader.openIn();
		uint32_t count = 0;

		while(reader.readNextRead(seq)){
			allContigsReadLengths.emplace_back(len(seq));
			nameToPositionKey[seq.name_] = count;
			readLengths[seq.name_] = len(seq);
			nameKey[estd::to_string(count)] = seq.name_;
			seq.name_ = estd::to_string(count);
			if (len(seq) <= lencutOffForBowtie) {
				tempForBowtie2Writer.openWrite(seq);
			} else {
				tempForLastzWriter.openWrite(seq);
			}
			++count;
		}
	}
	njh::sort(allContigsReadLengths, [](uint32_t len1, uint32_t len2){
		return len1 > len2;
	});
	//calculating assembly stats
	uint32_t l50 = 0;
	uint32_t n50 = 0;
	uint32_t lg50 = 0;
	uint32_t ng50 = 0;
	{
		auto allContigsSum = vectorSum(allContigsReadLengths);
		{
			uint32_t contigSum = 0;
			uint32_t contigCount = 0;
			for(const auto & leng : allContigsReadLengths){
				contigSum += leng;
				++contigCount;
				if(contigSum > 0.5* allContigsSum	){
					n50 = contigSum;
					l50 = contigCount;
					break;
				}
			}
		}
		{
			uint32_t contigGenomeSum = 0;
			uint32_t contigGenomeCount = 0;
			for(const auto & leng : allContigsReadLengths){
				contigGenomeSum += leng;
				++contigGenomeCount;
				if(contigGenomeSum > 0.5* sumTotalRequired	){
					ng50 = contigGenomeSum;
					lg50 = contigGenomeCount;
					break;
				}
			}
		}
	}
	auto tempInBowtie2Opts = SeqIOOptions::genFastaIn(tempSeqBowtie2Opts.out_.outName());
	tempInBowtie2Opts.out_.transferOverwriteOpts(setUp.pars_.ioOptions_.out_);

	auto tempLastzInOpts = SeqIOOptions::genFastaIn(tempSeqLastzOpts.out_.outName());
	tempLastzInOpts.out_.transferOverwriteOpts(setUp.pars_.ioOptions_.out_);

	{
		table nameKeyTab(nameKey, VecStr{"key", "name"});
		OutOptions nameKeyOpts(njh::files::make_path(setUp.pars_.directoryName_, "nameKey.tab.txt"));
		OutputStream nameKeyOut(nameKeyOpts);
		nameKeyTab.outPutContents(nameKeyOut, "\t");
	}
	//align to genomes in parallel

	struct GenomeWithProgram {
		GenomeWithProgram() = default;
		GenomeWithProgram(bool lastz, std::string genome) :
				lastz_(lastz), genome_(std::move(genome)) {
		}
		bool lastz_ { false };
		std::string genome_;
	};
	std::vector<GenomeWithProgram> genomesWithPrograms;
	auto inputGenomes = getVectorOfMapKeys(gMapper.genomes_);
	for(const auto & genome : inputGenomes){
		genomesWithPrograms.emplace_back(false, genome);
		genomesWithPrograms.emplace_back(true, genome);
		njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar{genome});
	}
	njh::concurrent::LockableQueue<GenomeWithProgram> genomesQueue(genomesWithPrograms);
	std::function<void()> alignToGenome = [&setUp,&genomesQueue,&gMapper,&bRunner,&tempLastzInOpts,&tempInBowtie2Opts,&lzPars,&nBowtie2Alns](){
		GenomeWithProgram gAndProgram;
		while(genomesQueue.getVal(gAndProgram)){
			bfs::path genomeDir = njh::files::make_path(setUp.pars_.directoryName_, gAndProgram.genome_);
			if(gAndProgram.lastz_){
				if(tempLastzInOpts.inExists()){
					//align seqs
					auto seqOpts = tempLastzInOpts;
					seqOpts.out_.outFilename_ = njh::files::make_path(genomeDir, "alignedSeqsLastz.sorted.bam");

					BioCmdsUtils::LastZPars copyLzPars = lzPars;
					copyLzPars.genomeFnp = gMapper.genomes_.at(gAndProgram.genome_)->fnpTwoBit_;

					auto runOut = bRunner.lastzAlign(seqOpts, copyLzPars);
					OutOptions lastzAlignLogOpts(njh::files::make_path(genomeDir, "lastzLog.json"));
					OutputStream lastzAlignLogOut(lastzAlignLogOpts);
					lastzAlignLogOut << njh::json::toJson(runOut) << std::endl;
				}
			} else {
				if(tempInBowtie2Opts.inExists()){
					//align seqs
					auto seqOpts = tempInBowtie2Opts;
					seqOpts.out_.outFilename_ = njh::files::make_path(genomeDir, "alignedSeqsBowtie2.sorted.bam");

					auto runOut = bRunner.bowtie2Align(seqOpts, gMapper.genomes_.at(gAndProgram.genome_)->fnp_, njh::pasteAsStr("-k ", nBowtie2Alns));
					OutOptions bowtie2AlignLogOpts(njh::files::make_path(genomeDir, "bowtie2Log.json"));
					OutputStream bowtie2AlignLogOut(bowtie2AlignLogOpts);
					bowtie2AlignLogOut << njh::json::toJson(runOut) << std::endl;
				}
			}
		}
	};

	njh::concurrent::runVoidFunctionThreaded(alignToGenome, numThreads);
	for(const auto & genome : gMapper.genomes_){
		bfs::path genomeDir = njh::files::make_path(setUp.pars_.directoryName_, genome.first);
		auto bowtie2AlignedFnp = njh::files::make_path(genomeDir, "alignedSeqsBowtie2.sorted.bam");
		auto lastzAlignedFnp = njh::files::make_path(genomeDir, "alignedSeqsLastz.sorted.bam");
		std::unordered_map<std::string, std::vector<BamTools::BamAlignment>> bamAligns;
		std::unordered_set<std::string> mappedReads;
		BamTools::RefVector refData;
		if(bfs::exists(bowtie2AlignedFnp)){
			auto seqOpts = setUp.pars_.ioOptions_;
			seqOpts.out_.outFilename_ = bowtie2AlignedFnp;
			//extract locations and mapping stats
			BamTools::BamAlignment bAln;
			BamTools::BamReader bReader;
			bReader.Open(seqOpts.out_.outFilename_.string());
			checkBamOpenThrow(bReader, seqOpts.out_.outFilename_.string());
			refData = bReader.GetReferenceData();
			while (bReader.GetNextAlignment(bAln)) {
				if (bAln.IsMapped()) {
					bAln.Name = nameKey[bAln.Name];
					bamAligns[bAln.Name].emplace_back(bAln);
					mappedReads.emplace(bAln.Name);
				}
			}
		}
		if(bfs::exists(lastzAlignedFnp)){
			auto seqOpts = setUp.pars_.ioOptions_;
			seqOpts.out_.outFilename_ = lastzAlignedFnp;
			//extract locations and mapping stats
			BamTools::BamAlignment bAln;
			BamTools::BamReader bReader;
			bReader.Open(seqOpts.out_.outFilename_.string());
			checkBamOpenThrow(bReader, seqOpts.out_.outFilename_.string());

			uint32_t lastzAlnsReadIn = 0;
			while (bReader.GetNextAlignment(bAln)) {
				++lastzAlnsReadIn;

				if (bAln.IsMapped()) {
					bAln.Name = nameKey[bAln.Name];
					bamAligns[bAln.Name].emplace_back(bAln);
					mappedReads.emplace(bAln.Name);
				}
			}
			if(lastzAlnsReadIn > 0){
				//if lastz doesn't align anyting then no samview is created so there's no header
				refData = bReader.GetReferenceData();
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
				results->setComparison(true, alignerObj);
				writer.write(results->refSeqAligned_);
				writer.write(results->alnSeqAligned_);
				alignedRegions.emplace_back(results->gRegion_.genBedRecordCore());
				std::string appName;
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

	comparisonOut << "ReadNumber\tReadId\tBestRef\tscore\talnScore\tkDist-5\t1bIndel\t2bIndel\t>2bIndel\tlqMismatch\thqMismatch"
			<< '\t' << "totalErrors"
			<< '\t' << "length"
			<< '\t' << "program"
			<< '\t' << "sample"
			<< '\t' << "target"
			<< std::endl;




	uint32_t readNumber = 0;
	std::unordered_map<std::string, std::vector<Bed6RecordCore>> bestRegionsByGenome;
	std::unordered_map<std::string, std::string> regionNameToInputName;

	auto allAlnResultsKeys = getVectorOfMapKeys(allAlnResults);
	njh::sort(allAlnResultsKeys, [&nameToPositionKey](const std::string & name1, const std::string & name2){
		return nameToPositionKey[name1] < nameToPositionKey[name2];
	});

	auto refAlnOutErrorOpts = SeqIOOptions::genFastaOut(njh::files::make_path(setUp.pars_.directoryName_, "refAlignmentsWithErrors.fasta"));
	SeqOutput refAlnOutErrorOut(refAlnOutErrorOpts);
	refAlnOutErrorOut.openOut();

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
					regionNameToInputName[res->gRegion_.genBedRecordCore().genUIDFromCoordsWithStrand()] = res->bAln_.Name;
				}else if(res->comp_.alnScore_  == bestScore){
					bestResults.emplace_back(res);
					regionNameToGenome[res->gRegion_.genBedRecordCore().toDelimStrWithExtra()] = genomeRes.first;
					regionNameToInputName[res->gRegion_.genBedRecordCore().genUIDFromCoordsWithStrand()] = res->bAln_.Name;
				}
			}
		}
		for(const auto & results : bestResults){
			if (1 != results->comp_.distances_.eventBasedIdentityHq_) {
				refAlnOutErrorOut.write(results->refSeqAligned_);
				refAlnOutErrorOut.write(results->alnSeqAligned_);
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
			std::string appName;
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
					<< '\t' << results->comp_.hqMismatches_
					<< '\t' << results->comp_.distances_.getNumOfEvents(true)
					<< '\t' << readLengths[results->bAln_.Name]
					<< '\t' << program
					<< '\t' << sample
					<< '\t' << target
					<< std::endl;
		}
		++readNumber;
	}


	{
		std::unordered_map<std::string, std::set<std::string>> interceptedIDs;
		OutputStream allBestRegionsBedOut(njh::files::make_path(setUp.pars_.directoryName_, njh::pasteAsStr("bestRegions_", "all", ".bed")));
		for(auto & best : bestRegionsByGenome){
			OutOptions bestRegionsBedOpts(njh::files::make_path(setUp.pars_.directoryName_, njh::pasteAsStr("bestRegions_", best.first, ".bed")));
			OutputStream bestRegionsBedOut(bestRegionsBedOpts);
			if(bfs::exists(gMapper.genomes_.at(best.first)->gffFnp_)){
				intersectBedLocsWtihGffRecordsPars gffPars = gMapper.pars_.gffIntersectPars_;
				gffPars.gffFnp_ = gMapper.genomes_.at(best.first)->gffFnp_;
				intersectBedLocsWtihGffRecords(best.second, gffPars);
				for(const auto & region : best.second){
					if(!region.extraFields_.empty() && MetaDataInName::nameHasMetaData(region.extraFields_[0])){
						auto regionMeta = MetaDataInName(region.extraFields_[0]);
						if(regionMeta.containsMeta("ID")){
							interceptedIDs[best.first].emplace(regionMeta.getMeta("ID"));
						}
					}
				}
			}
			BedUtility::coordSort(best.second);
			for(const auto & reg : best.second){
				bestRegionsBedOut << reg.toDelimStrWithExtra() << std::endl;
				allBestRegionsBedOut << reg.toDelimStrWithExtra() << std::endl;
			}
		}
		for(const auto & genomeIDs : interceptedIDs){
			OutputStream interceptedIDsOut(njh::files::make_path(setUp.pars_.directoryName_, njh::pasteAsStr("bestRegions_", genomeIDs.first, "_IDs.txt")));
			interceptedIDsOut << njh::conToStr(genomeIDs.second, "\n") << std::endl;
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
					<< '\t' << "NA"
					<< '\t' << "NA"
					<< '\t' << "NA"
					<< '\t' << "NA"
					<< '\t' << "NA"
					<< '\t' << "NA"
					<< '\t' << "NA"
					<< '\t' << "NA"
					<< '\t' << "NA"
					<< '\t' << "NA"
					<< '\t' << readLengths[unmappedCount.first]
					<< '\t' << program
					<< '\t' << sample
					<< '\t' << target
					<< std::endl;
			/*
			 * 		coveredCountsTab.addColumn(VecStr{program}, "program");
		coveredCountsTab.addColumn(VecStr{sample}, "sample");
			 */
			++readNumber;
		}
	}

	std::set<std::string> matchingContigs;
	std::set<std::string> notMatchingContigs(unmappedToAllGenomes.begin(), unmappedToAllGenomes.end());

	if(!unmappedToAllGenomes.empty()){
		OutOptions unmmapedOpts(njh::files::make_path(setUp.pars_.directoryName_, "unmappedReads.txt"));
		OutputStream unmappedOut(unmmapedOpts);
		for(const auto & unmappedAln : unmappedToAllGenomes){
			unmappedOut << unmappedAln << std::endl;
		}
	}
	bfs::remove(tempSeqBowtie2Opts.out_.outName());
	bfs::remove(tempSeqLastzOpts.out_.outName());

	{
		// determine regions covered
		std::unordered_map<std::string, std::map<uint32_t, uint32_t>> simpleCoverageCounts;
		for(const auto & seqAlnResults : bestAlnResults){
			for(const auto & genome : seqAlnResults.second){
				for(const auto & res : genome.second){
					if(amountOfErrorForCoverageCalc.passErrorProfile(res->comp_)){
						matchingContigs.emplace(res->bAln_.Name);
						for(uint32_t pos = res->gRegion_.start_; pos <res->gRegion_.end_; ++pos){
							simpleCoverageCounts[res->gRegion_.chrom_][pos] += 1;
						}
					}else{
						notMatchingContigs.emplace(res->bAln_.Name);
					}
				}
			}
		}

		{
			//write out matching contigs info
			OutputStream matchInfo(njh::files::make_path(setUp.pars_.directoryName_, "contigsMatchingExpectedInfo.tab.txt"));
			matchInfo << "program\tsample\ttarget";
			matchInfo
					<< "\t" << "ContigsMatchingExpectedCnt"
					<< "\t" << "ContigsMatchingExpectedFrac"
					<< "\t" << "ContigsMatchingNotExpectedCnt"
					<< "\t" << "ContigsMatchingNotExpectedFrac"
					<< "\t" << "TotalContigs"
					<< "\t" << "TotalBases"
					<< "\t" << "TotalRequiredContigs"
					<< "\t" << "TotalRequiredBases"
					<< "\t" << "n50"
					<< "\t" << "l50"
					<< "\t" << "ng50"
					<< "\t" << "lg50"
					<< std::endl;
			matchInfo << program
					<< "\t" << sample
					<< "\t" << target
					<< "\t" << matchingContigs.size()
					<< "\t" << matchingContigs.size()/static_cast<double>(readNumber)
					<< "\t" << notMatchingContigs.size()
					<< "\t" << notMatchingContigs.size()/static_cast<double>(readNumber)
					<< "\t" << readNumber
					<< "\t" << vectorSum(allContigsReadLengths)
					<< "\t" << requiredRegions.size()
					<< "\t" << sumTotalRequired
					<< "\t" << n50
					<< "\t" << l50
					<< "\t" << ng50
					<< "\t" << lg50
					<< std::endl;

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

			allRegionsNotExpectedB6.reserve(allRegionsNotExpected.size());
			for(const auto & notExp : allRegionsNotExpected){
				allRegionsNotExpectedB6.emplace_back(std::make_shared<Bed6RecordCore>(GenomicRegion(notExp).genBedRecordCore()));
			}
			std::unordered_map<std::string, std::vector<std::shared_ptr<Bed6RecordCore>>> allRegionsByGenome;

			for(const auto & region : allRegionsNotExpectedB6){
				allRegionsByGenome[chromosomeToGenome[region->chrom_]].push_back(region);
			}
			for( auto & genome : allRegionsByGenome){
				if("" != gMapper.genomes_.at(genome.first)->gffFnp_){
					intersectBedLocsWtihGffRecordsPars interPars = gMapper.pars_.gffIntersectPars_;
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
		{
			std::vector<Bed6RecordCore> allRegionsNotCovered;
			table coveredCountsTab(VecStr{"Region", "basesCovered", "totalBases", "fractionCovered"});
			for(const auto & gene: requiredRegions){
				uint32_t covered = 0;
				std::vector<uint32_t> positionsNotCovered;
				for(const auto pos : iter::range(gene.start_, gene.end_)){
					if(simpleCoverageCounts[gene.chrom_][pos] > 0){
						++covered;
					}else{
						positionsNotCovered.emplace_back(pos);
					}
				}
				coveredCountsTab.addRow(gene.uid_, covered, gene.getLen(), static_cast<double>(covered)/gene.getLen());

				if(!positionsNotCovered.empty()){

					std::vector<Bed6RecordCore> regionsNotCoveredRaw;
					regionsNotCoveredRaw.reserve(positionsNotCovered.size());
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
						allRegionsNotCovered.emplace_back(region);
					}
				}
			}
			OutputStream allRegionsNotCoveredOut(njh::files::make_path(setUp.pars_.directoryName_, "allNotCoveredRegions.bed"));
			for(const auto & region : allRegionsNotCovered) {
				allRegionsNotCoveredOut << region.toDelimStr() << std::endl;
			}

			coveredCountsTab.addColumn(VecStr{program}, "program");
			coveredCountsTab.addColumn(VecStr{sample}, "sample");
			coveredCountsTab.addColumn(VecStr{target}, "target");

			coveredCountsTab.outPutContents(
					TableIOOpts::genTabFileOut(
							njh::files::make_path(setUp.pars_.directoryName_,
									"coveragedInfo.tab.txt"), true));
		}




		if(calcSpecificCoverage){
			table specifcCoveredCountsTab(VecStr{"Region", "basesCovered", "totalBases", "fractionCovered"});
			std::vector<uint32_t> positionsNotCovered;
			for(const auto & regions : specificRegions){
				uint32_t totalRegionBases = 0;
				uint32_t totalCovered = 0;
				for(const auto & reg : regions.second){
					for(const auto pos : iter::range(reg.start_, reg.end_)){
						++totalRegionBases;
						if(simpleCoverageCounts[reg.chrom_][pos] > 0){
							++totalCovered;
						}else{
							positionsNotCovered.emplace_back(pos);
						}
					}
				}
				specifcCoveredCountsTab.addRow(regions.first, totalCovered, totalRegionBases, static_cast<double>(totalCovered)/totalRegionBases);
				if(!positionsNotCovered.empty()){
					std::vector<Bed6RecordCore> regionsNotCoveredRaw;
					regionsNotCoveredRaw.reserve(positionsNotCovered.size());
					for(const auto & pos : positionsNotCovered){
						regionsNotCoveredRaw.emplace_back(Bed6RecordCore(regions.second.front().chrom_, pos, pos + 1, regions.second.front().uid_, 1, regions.second.front().reverseSrand_? '-':'+'));
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
						region.name_ = njh::pasteAsStr(regions.first, ":", GenomicRegion(region).createUidFromCoords());
					}
					OutputStream out(njh::files::make_path(setUp.pars_.directoryName_, regions.first + "_specificNotCoveredRegions.bed"));
					for(const auto & region : regionsNotCovered){
						out << region.toDelimStr() << std::endl;
					}
				}
			} //
			specifcCoveredCountsTab.addColumn(VecStr{program}, "program");
			specifcCoveredCountsTab.addColumn(VecStr{sample}, "sample");
			specifcCoveredCountsTab.addColumn(VecStr{target}, "target");
			specifcCoveredCountsTab.outPutContents(
					TableIOOpts::genTabFileOut(
							njh::files::make_path(setUp.pars_.directoryName_,
									"specifcRegionsCoveragedInfo.tab.txt"), true));
		}
	}

	{
		//read lengths
		auto allRLens = getVectorOfMapValues(readLengths);
		auto stats = getStatsOnVec(allRLens);
		table statsTable;
		statsTable.columnNames_ = getVectorOfMapKeys(stats);
		statsTable.hasHeader_ = true;
		auto nums = getVectorOfMapValues(stats);
		statsTable.content_.emplace_back(numVecToVecStr(nums));
		addOtherVec(statsTable.columnNames_, VecStr{"n","program", "sample", "target"});
		addOtherVec(statsTable.content_.front(), toVecStr(allRLens.size(), program, sample, target));
		statsTable.outPutContents(
				TableIOOpts::genTabFileOut(
						njh::files::make_path(setUp.pars_.directoryName_,
								"contigsLengthsInfo.tab.txt"), true));
	}

	//grouping
	OutputStream groupedRegionsOut(njh::files::make_path(setUp.pars_.directoryName_, "groupedRegions.bed"));

	uint32_t distanceWithin = 100;
	for(const auto & expectedRegion : requiredRegions){
		std::vector<Bed6RecordCore> associatedRegions;
		for(const auto & best : bestRegionsByGenome){
			for(const auto & region : best.second){
				if(region.getDistanceBetween(expectedRegion.genBed3RecordCore()) < distanceWithin){
					associatedRegions.emplace_back(region);
				}
			}
		}
		uint32_t count = 0;
		for(const auto & region : associatedRegions){
			groupedRegionsOut << region.toDelimStr() << "\t" << expectedRegion.uid_ << "\t" << regionNameToInputName[region.genUIDFromCoordsWithStrand()] << "\t" << count << std::endl;
			++count;
		}
	}

	return 0;
}



}  // namespace njhseq

