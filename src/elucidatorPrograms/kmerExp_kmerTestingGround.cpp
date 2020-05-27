/*
 * kmerExp_kmerTestingGround.cpp
 *
 *  Created on: May 24, 2020
 *      Author: nick
 */




#include "kmerExp.hpp"
#include "elucidator/utils/KmerUtils.hpp"
#include "elucidator/objects/seqObjects/seqKmers.h"

#include <njhseq/concurrency/AllByAllPairFactory.hpp>

namespace njhseq {


int kmerExpRunner::kmerTestingGround(const njh::progutils::CmdArgs & inputCommands){
	uint32_t klen = 35;
	MultiGenomeMapper::inputParameters pars;




	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(pars.numThreads_, "--numThreads", "Number of CPUs to utilize");
	setUp.setOption(pars.genomeDir_, "--genomeDir", "Name of the genome file fnp", true);
	setUp.setOption(pars.primaryGenome_, "--primaryGenome", "The primary reference genome");
  setUp.setOption(pars.selectedGenomes_, "--selectedGenomes", "Name of the other genomes in --genomeDir to be read in, leave blank to just do all fastas");


	setUp.setOption(pars.gffDir_, "--gffDir", "A directory with a gff for the genomes in --genomeDir, should be named GENOME.gff (for GENOME.fasta)");
  setUp.setOption(pars.gffIntersectPars_.extraAttributes_, "--gffExtraAttributes", "Extra attributes to add to genome that has an accompany gff");


	setUp.setOption(klen, "--klen", "kmer Length", true);
	setUp.processDirectoryOutputName(pars.primaryGenome_ + "_kmerTestingGround" + "_TODAY", true);
	setUp.finishSetUp(std::cout);

	setUp.startARunLog(setUp.pars_.directoryName_);
  njh::stopWatch watch;

  watch.setLapName("Set up genomes");

	//set up genome mapper;
	auto gMapper = std::make_unique<MultiGenomeMapper>(pars);
	//set up selected genomes
	if(!pars.selectedGenomes_.empty()){
		gMapper->setSelectedGenomes(pars.selectedGenomes_);
	}
	//init is threaded
	gMapper->init();

	auto genomeNames = getVectorOfMapKeys(gMapper->genomes_);

	std::unordered_map<std::string, std::unordered_map<std::string, std::shared_ptr<KmersSharedBlocks>>> uniqueKmersPerGenomePerRecord;
	std::mutex uniqKsMut;

	njh::concurrent::LockableQueue<std::string> genomeNamesQueue(genomeNames);
	std::function<void()> indexGenomeForUniqueKmers = [&uniqueKmersPerGenomePerRecord,&uniqKsMut,&genomeNamesQueue,&gMapper,&klen](){
		std::string genome = "";
		while(genomeNamesQueue.getVal(genome)){
		  std::unordered_map<std::string, std::shared_ptr<KmersSharedBlocks>> uniqueKmersPerRecord;
		  {
		  	std::unordered_map<std::string, uint32_t> genomeWideKmers;
		  	auto genomeFnpInOpts = SeqIOOptions::genFastaIn(gMapper->genomes_.at(genome)->fnp_);
		  	genomeFnpInOpts.includeWhiteSpaceInName_ = false;
		  	{
		  		//get genome wide kmers
		    	SeqInput reader(genomeFnpInOpts);
		    	reader.openIn();
		    	seqInfo seq;
		    	while (reader.readNextRead(seq)) {
		    		if (len(seq) > klen + 1) {
		    			for (uint32_t pos = 0; pos < len(seq) - klen + 1; ++pos) {
		    				++genomeWideKmers[seq.seq_.substr(pos, klen)];
		    			}
		    		}
		    	}
		  	}
		  	{
		  		//get genome wide kmers
		  		SeqInput reader(genomeFnpInOpts);
		    	reader.openIn();
		    	seqInfo seq;
		    	while (reader.readNextRead(seq)) {
		    		if (len(seq) > klen + 1) {
		    			kmerInfo uniqueKmers;
		    			uniqueKmers.kLen_ = klen;
		    			for (uint32_t pos = 0; pos < len(seq) - klen + 1; ++pos) {
		    				auto k = seq.seq_.substr(pos, klen);
		    				if(1 == genomeWideKmers[k]){
		    					uniqueKmers.kmers_[k] = kmer(k,pos);
		    				}
		    			}
		    			uniqueKmersPerRecord[seq.name_] = std::make_shared<KmersSharedBlocks>(seq, uniqueKmers);
		    		}
		    	}
		  	}
		  }
		  {
		  	std::lock_guard<std::mutex> lock(uniqKsMut);
		  	uniqueKmersPerGenomePerRecord[genome] = uniqueKmersPerRecord;
		  }
		}
	};
	watch.startNewLap("Get unique kmers for all genomes");
	njh::concurrent::runVoidFunctionThreaded(indexGenomeForUniqueKmers, gMapper->pars_.numThreads_);


	auto perGenomeUniCountsDir = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar{"uniqueKmersCounts"});
	watch.startNewLap("Outputing info per genome");

	for(const auto & genome : uniqueKmersPerGenomePerRecord){
		OutputStream out(njh::files::make_path(perGenomeUniCountsDir, genome.first + "_uniqueKmerNumbers.tab.txt"));
		out << "record1\tnumberOfUniqueKmers" << std::endl;
		auto recordNames = getVectorOfMapKeys(genome.second);
		njh::sort(recordNames);
		for(const auto & name : recordNames){
			out << name
					<< "\t" << genome.second.at(name)->kInfo_->kmers_.size() << std::endl;
		}

		OutputStream outBed(njh::files::make_path(perGenomeUniCountsDir, genome.first + "_uniqueKmerNumbers.bed"));
		for(const auto & name : recordNames){
			std::vector<uint32_t> positions;
			for(const auto & k : genome.second.at(name)->kInfo_->kmers_){
				positions.emplace_back(k.second.positions_.front());
			}
			njh::sort(positions);
			if(positions.size() > 0){
				std::vector<Bed6RecordCore> regions;
				regions.emplace_back(Bed6RecordCore(name, positions[0], positions[0] +1, "", 1, '+'));
				if(positions.size() > 1){
					for(uint32_t pos = 1; pos < positions.size(); ++pos){
						auto gPos = positions[pos];
						if(regions.back().chromEnd_ == gPos){
							regions.back().chromEnd_ += 1;
						}else{
							regions.emplace_back(Bed6RecordCore(name, positions[pos], positions[pos] +1, "", 1, '+'));
						}
					}
				}
				njh::for_each(regions,[&klen](Bed6RecordCore & record){
					record.chromEnd_ += klen - 1;
					record.score_ = record.length();
					record.name_ = record.genUIDFromCoords();
				});
				for(const auto & rec : regions){
					outBed << rec.toDelimStr() << std::endl;
				}
			}
		}
	}
	watch.startNewLap("Compare Against Ref");

	auto refRecNames = getVectorOfMapKeys(uniqueKmersPerGenomePerRecord[pars.primaryGenome_]);
	struct UniKmersCompRes {
		UniKmersCompRes(const std::string &gRec, const std::string &rRec,
				uint32_t genomeUniqKs, uint32_t refUniKs) :
				genomeRec_(gRec), refRec_(rRec), genomeUniqKs_(genomeUniqKs), refUniKs_(
						refUniKs) {

		}
		std::string genomeRec_;
		std::string refRec_;
		uint32_t genomeUniqKs_ { 0 };
		uint32_t refUniKs_ { 0 };
		uint32_t kmersShared_ { 0 };

		double proportionOfGenomeUniqShared() const {
			return static_cast<double>(kmersShared_) / genomeUniqKs_;
		}
		double indexScore() const {
			return static_cast<double>(kmersShared_)
					/ (genomeUniqKs_ + refUniKs_ - kmersShared_);
		}

	};
	std::unordered_map<std::string, std::vector<UniKmersCompRes>> compsAgainstRef;
	for(const auto & genome : uniqueKmersPerGenomePerRecord){
		if(genome.first != pars.primaryGenome_){
			OutputStream genomeCompFile(njh::files::make_path(setUp.pars_.directoryName_, genome.first + "_againstRefGenome.tab.txt"));
			genomeCompFile << "record\trefRecord\ttotalUniqueKs\ttotalUniqueKsRef\tkShared\tindexScore\tproportionOfKsShared" << std::endl;
			std::mutex genomeCompMut;
			auto genomeRecNames = getVectorOfMapKeys(genome.second);
			AllByAllPairFactory pairFactory(genome.second.size(), uniqueKmersPerGenomePerRecord.at(pars.primaryGenome_).size());
			std::function<void()> compToRef = [&refRecNames,&genomeRecNames,&pairFactory,&genome,&genomeCompMut,&compsAgainstRef,&pars,&uniqueKmersPerGenomePerRecord](){
				AllByAllPairFactory::AllByAllPair pair;
				while(pairFactory.setNextPair(pair)){
					std::string refRec = refRecNames[pair.col_];
					std::string genomeRec = genomeRecNames[pair.row_];
					UniKmersCompRes compRes(genomeRec, refRec, genome.second.at(genomeRec)->kInfo_->kmers_.size(), uniqueKmersPerGenomePerRecord[pars.primaryGenome_][refRec]->kInfo_->kmers_.size());
					for(const auto & k : genome.second.at(genomeRec)->kInfo_->kmers_){
						if(njh::in(k.first, uniqueKmersPerGenomePerRecord[pars.primaryGenome_][refRec]->kInfo_->kmers_)){
							++compRes.kmersShared_;
						}
					}
					{
						std::lock_guard<std::mutex> lock(genomeCompMut);
						compsAgainstRef[genome.first].emplace_back(compRes);
					}
				}
			};
			njh::concurrent::runVoidFunctionThreaded(compToRef, pars.numThreads_);
			njh::sort(compsAgainstRef[genome.first], [](const UniKmersCompRes & comp1, const UniKmersCompRes & comp2){
				if(comp1.genomeRec_ == comp2.genomeRec_){
					return comp1.refRec_ < comp2.refRec_;
				}
				return comp1.genomeRec_ < comp2.genomeRec_;
			});
			for(const auto & comp : compsAgainstRef[genome.first]){
				genomeCompFile << comp.genomeRec_
						<< "\t" << comp.refRec_
						<< "\t" << comp.genomeUniqKs_
						<< "\t" << comp.refUniKs_
						<< "\t" << comp.kmersShared_
						<< "\t" << comp.indexScore()
						<< "\t" << comp.proportionOfGenomeUniqShared() << std::endl;
			}
		}
	}

//
//  OutputStream out(njh::files::make_path(setUp.pars_.directoryName_, "uniqueKmerNumbers.tab.txt"));
//  std::unordered_map<std::string, std::shared_ptr<KmersSharedBlocks>> uniqueKmersPerRecord;
//  {
//  	std::unordered_map<std::string, uint32_t> genomeWideKmers;
//
//  	{
//  		//get genome wide kmers
//    	SeqInput reader(setUp.pars_.ioOptions_);
//    	reader.openIn();
//    	seqInfo seq;
//    	while (reader.readNextRead(seq)) {
//    		if (len(seq) > klen + 1) {
//    			for (uint32_t pos = 0; pos < len(seq) - klen + 1; ++pos) {
//    				++genomeWideKmers[seq.seq_.substr(pos, klen)];
//    			}
//    		}
//    	}
//  	}
//  	watch.startNewLap("Get unique kmers");
//  	{
//  		//get genome wide kmers
//    	SeqInput reader(setUp.pars_.ioOptions_);
//    	reader.openIn();
//    	seqInfo seq;
//    	while (reader.readNextRead(seq)) {
//    		if (len(seq) > klen + 1) {
//    			kmerInfo uniqueKmers;
//    			uniqueKmers.kLen_ = klen;
//    			for (uint32_t pos = 0; pos < len(seq) - klen + 1; ++pos) {
//    				auto k = seq.seq_.substr(pos, klen);
//    				if(1 == genomeWideKmers[k]){
//    					uniqueKmers.kmers_[k] = kmer(k,pos);
//    				}
//    			}
//    			uniqueKmersPerRecord[seq.name_] = std::make_shared<KmersSharedBlocks>(seq, uniqueKmers);
//    		}
//    	}
//  	}
//  }
//
//
//	out << "record1\tnumberOfUniqueKmers" << std::endl;
//	auto recordNames = getVectorOfMapKeys(uniqueKmersPerRecord);
//	njh::sort(recordNames);
//	for(const auto & name : recordNames){
//		out << name
//				<< "\t" << uniqueKmersPerRecord[name]->kInfo_->kmers_.size() << std::endl;
//	}
//
//  OutputStream outBed(njh::files::make_path(setUp.pars_.directoryName_, "uniqueKmerNumbers.bed"));
//
//
//	for(const auto & name : recordNames){
//		std::vector<uint32_t> positions;
//		for(const auto & k : uniqueKmersPerRecord[name]->kInfo_->kmers_){
//			positions.emplace_back(k.second.positions_.front());
//		}
//		njh::sort(positions);
//		if(positions.size() > 0){
//			std::vector<Bed6RecordCore> regions;
//			regions.emplace_back(Bed6RecordCore(name, positions[0], positions[0] +1, "", 1, '+'));
//			if(positions.size() > 1){
//				for(uint32_t pos = 1; pos < positions.size(); ++pos){
//					auto gPos = positions[pos];
//					if(regions.back().chromEnd_ == gPos){
//						regions.back().chromEnd_ += 1;
//					}else{
//						regions.emplace_back(Bed6RecordCore(name, positions[pos], positions[pos] +1, "", 1, '+'));
//					}
//				}
//			}
//			njh::for_each(regions,[&klen](Bed6RecordCore & record){
//				record.chromEnd_ += klen - 1;
//				record.score_ = record.length();
//				record.name_ = record.genUIDFromCoords();
//			});
//			for(const auto & rec : regions){
//				outBed << rec.toDelimStr() << std::endl;
//			}
//		}
//	}
//

  if(setUp.pars_.verbose_){
  	watch.logLapTimes(std::cout, true, 6, true);
  	std::cout << "totalTime: " << watch.totalTimeFormatted(6) << std::endl;
  }else{
  	OutputStream watchLogTime(njh::files::make_path(setUp.pars_.directoryName_, "time.txt"));
  	watch.logLapTimes(watchLogTime, true, 6, true);
  	watchLogTime << "totalTime: " << watch.totalTimeFormatted(6) << std::endl;
  }

	return 0;
}

int kmerExpRunner::getWithinGenomeUniqueKmers(const njh::progutils::CmdArgs & inputCommands){
	uint32_t klen = 35;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processReadInNames(true);
	setUp.setOption(klen, "--klen", "kmer Length", true);
	setUp.processDirectoryOutputName(true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
  njh::stopWatch watch;
  OutputStream out(njh::files::make_path(setUp.pars_.directoryName_, "uniqueKmerNumbers.tab.txt"));
  std::unordered_map<std::string, std::shared_ptr<KmersSharedBlocks>> uniqueKmersPerRecord;
  {
  	std::unordered_map<std::string, uint32_t> genomeWideKmers;
  	watch.setLapName("Get genome wide kmers");
  	{
  		//get genome wide kmers
    	SeqInput reader(setUp.pars_.ioOptions_);
    	reader.openIn();
    	seqInfo seq;
    	while (reader.readNextRead(seq)) {
    		if (len(seq) > klen + 1) {
    			for (uint32_t pos = 0; pos < len(seq) - klen + 1; ++pos) {
    				++genomeWideKmers[seq.seq_.substr(pos, klen)];
    			}
    		}
    	}
  	}
  	watch.startNewLap("Get unique kmers");
  	{
  		//get genome wide kmers
    	SeqInput reader(setUp.pars_.ioOptions_);
    	reader.openIn();
    	seqInfo seq;
    	while (reader.readNextRead(seq)) {
    		if (len(seq) > klen + 1) {
    			kmerInfo uniqueKmers;
    			uniqueKmers.kLen_ = klen;
    			for (uint32_t pos = 0; pos < len(seq) - klen + 1; ++pos) {
    				auto k = seq.seq_.substr(pos, klen);
    				if(1 == genomeWideKmers[k]){
    					uniqueKmers.kmers_[k] = kmer(k,pos);
    				}
    			}
    			uniqueKmersPerRecord[seq.name_] = std::make_shared<KmersSharedBlocks>(seq, uniqueKmers);
    		}
    	}
  	}
  }


	out << "record1\tnumberOfUniqueKmers" << std::endl;
	auto recordNames = getVectorOfMapKeys(uniqueKmersPerRecord);
	njh::sort(recordNames);
	for(const auto & name : recordNames){
		out << name
				<< "\t" << uniqueKmersPerRecord[name]->kInfo_->kmers_.size() << std::endl;
	}

  OutputStream outBed(njh::files::make_path(setUp.pars_.directoryName_, "uniqueKmerNumbers.bed"));


	for(const auto & name : recordNames){
		std::vector<uint32_t> positions;
		for(const auto & k : uniqueKmersPerRecord[name]->kInfo_->kmers_){
			positions.emplace_back(k.second.positions_.front());
		}
		njh::sort(positions);
		if(positions.size() > 0){
			std::vector<Bed6RecordCore> regions;
			regions.emplace_back(Bed6RecordCore(name, positions[0], positions[0] +1, "", 1, '+'));
			if(positions.size() > 1){
				for(uint32_t pos = 1; pos < positions.size(); ++pos){
					auto gPos = positions[pos];
					if(regions.back().chromEnd_ == gPos){
						regions.back().chromEnd_ += 1;
					}else{
						regions.emplace_back(Bed6RecordCore(name, positions[pos], positions[pos] +1, "", 1, '+'));
					}
				}
			}
			njh::for_each(regions,[&klen](Bed6RecordCore & record){
				record.chromEnd_ += klen - 1;
				record.score_ = record.length();
				record.name_ = record.genUIDFromCoords();
			});
			for(const auto & rec : regions){
				outBed << rec.toDelimStr() << std::endl;
			}
		}
	}


  if(setUp.pars_.verbose_){
  	watch.logLapTimes(std::cout, true, 6, true);
  	std::cout << "totalTime: " << watch.totalTimeFormatted(6) << std::endl;
  }
	return 0;
}

int kmerExpRunner::pairwiseComparisonOfUniqueKmers(const njh::progutils::CmdArgs & inputCommands){
	uint32_t klen = 35;
	uint32_t numThreads = 1;
	OutOptions outOpts(bfs::path(""), ".tab.txt");
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processReadInNames(true);
	setUp.processWritingOptions(outOpts);
	setUp.setOption(numThreads, "--numThreads", "Number of threads to use", njh::progutils::ProgramSetUp::CheckCase::GTEQ1);
	setUp.setOption(klen, "--klen", "kmer Length", true);
	setUp.finishSetUp(std::cout);

  njh::stopWatch watch;
  watch.setLapName("Get Unique K-mers");
  OutputStream out(outOpts);

	SeqInput reader(setUp.pars_.ioOptions_);
	std::unordered_map<std::string, std::unordered_map<std::string, kmer>> uniqueKmersPerRecord;
	reader.openIn();
	seqInfo seq;
	while (reader.readNextRead(seq)) {
		std::unordered_map<std::string, kmer> allKmers;
		if (len(seq) > klen + 1) {
			for (uint32_t pos = 0; pos < len(seq) - klen + 1; ++pos) {
				auto k = seq.seq_.substr(pos, klen);
				if(njh::in(k, allKmers)){
					allKmers[k].addPosition(pos);
				}else{
					allKmers[k] = kmer(k, pos);
				}
			}
		}
		for (const auto &k : allKmers) {
			if (1 == k.second.count_) {
				uniqueKmersPerRecord[seq.name_].emplace(k);
			}
		}
	}
	watch.startNewLap("comparing");
	PairwisePairFactory pFac(uniqueKmersPerRecord.size());


	out << "record1\trecord2\tr1TotalUniqueKmers\tr2TotalUniqueKmers\ttotalUniqueKmersShared\tshareScore\tshareScore2" << std::endl;
	auto recordNames = getVectorOfMapKeys(uniqueKmersPerRecord);
	njh::sort(recordNames);

	std::mutex outMut;

	std::function<void()> runComp = [&out,&outMut,&pFac,&uniqueKmersPerRecord,recordNames](){
		PairwisePairFactory::PairwisePair pair;
		while(pFac.setNextPair(pair)){
			uint32_t shared = 0;
			for(const auto & k :uniqueKmersPerRecord[recordNames[pair.col_]]){
				if(njh::in(k.first, uniqueKmersPerRecord[recordNames[pair.row_]])){
					++shared;
				}
			}
			{
				std::lock_guard<std::mutex> lock(outMut);
				out << recordNames[pair.col_]
						<< "\t" << recordNames[pair.row_]
						<< "\t" << uniqueKmersPerRecord[recordNames[pair.col_]].size()
						<< "\t" << uniqueKmersPerRecord[recordNames[pair.row_]].size()
						<< "\t" << shared
						<< "\t" << static_cast<double>(shared)/std::min(uniqueKmersPerRecord[recordNames[pair.col_]].size(), uniqueKmersPerRecord[recordNames[pair.row_]].size())
				<< "\t" << (static_cast<double>(shared))/(uniqueKmersPerRecord[recordNames[pair.col_]].size() +  uniqueKmersPerRecord[recordNames[pair.row_]].size() - shared)
				<< std::endl;
			}
		}
	};

	njh::concurrent::runVoidFunctionThreaded(runComp, numThreads);



  if(setUp.pars_.verbose_){
  	watch.logLapTimes(std::cout, true, 6, true);
  	std::cout << "totalTime: " << watch.totalTimeFormatted(6) << std::endl;
  }
	return 0;
}



}  // namespace njhseq
