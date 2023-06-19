/*
 * bamExp_BamFilterAlnsOnLenAndQual.cpp
 *
 *  Created on: Feb 17, 2020
 *      Author: nick
 */





#include "bamExp.hpp"
#include <TwoBit.h>

#include "elucidator/BamToolsUtils.h"
#include <njhseq/objects/BioDataObject.h>


namespace njhseq {



int bamExpRunner::BamFilterAlnsOnLenAndQual(const njh::progutils::CmdArgs & inputCommands){
	/**@todo not finished yet */
	OutOptions outOpts(bfs::path("out"));
	uint32_t lenCutOff = 50;
	double qualCheckCutOff = 0.50;
	uint32_t qualCheck = 30;

	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processReadInNames({"--bam"}, true);
	setUp.setOption(qualCheckCutOff, "--qualCheckCutOff", "Fraction of a read has to be above this of the given quality by --qualCheck");
	setUp.setOption(qualCheck, "--qualCheck", "Per base quality check to for filtering based off of quality");
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	ReadCheckerQualCheck qualChecker(qualCheck, qualCheckCutOff, false);

	BamTools::BamReader bReader;
	BamTools::BamAlignment bAln;
	bReader.Open(setUp.pars_.ioOptions_.firstName_.string());
	checkBamOpenThrow(bReader, setUp.pars_.ioOptions_.firstName_);
	loadBamIndexThrow(bReader);
	auto refs = bReader.GetReferenceData();
	BamAlnsCache alnCache;
	auto refData = bReader.GetReferenceData();
	std::unordered_map<std::string, uint32_t> refNameToId;
	for (auto pos : iter::range(refData.size())) {
		refNameToId[refData[pos].RefName] = pos;
	}
	//pair writer
	SeqIOOptions outOptsPaired(outOpts.outFilename_, SeqIOOptions::outFormats::FASTQPAIREDGZ, outOpts);
	SeqOutput pairWriter(outOptsPaired);

	//non paired writer
	SeqIOOptions outOptsSingle(outOpts.outFilename_, SeqIOOptions::outFormats::FASTQGZ, outOpts);
	SeqOutput writer(outOptsSingle);



	while (bReader.GetNextAlignment(bAln)) {
		//skip secondary alignments
		if (!bAln.IsPrimaryAlignment()) {
			continue;
		}

		if (bAln.IsPaired()) {
			if ((bAln.MateRefID != bAln.RefID)
					|| (!bAln.IsMapped() || !bAln.IsMateMapped())) {
				// do non-concordant chromosome mapping operation or non-mapping mates
				if (!alnCache.has(bAln.Name)) {
					//pair hasn't been added to cache yet so add to cache
					//this only works if mate and first mate have the same name
					alnCache.add(bAln);
					continue;
				} else {
					auto search = alnCache.get(bAln.Name);
					auto balnSeq = bamAlnToSeqInfo(bAln);
					auto searchSeq = bamAlnToSeqInfo(*search);


					bool balnFailed = false;
					if(len(balnSeq) < lenCutOff){
						balnFailed = true;
					}else{
						qualChecker.checkRead(balnSeq);
						balnFailed = balnSeq.on_;
					}

					bool searchFailed = false;
					if(len(searchSeq) < lenCutOff){
						searchFailed = true;
					}else{
						qualChecker.checkRead(searchSeq);
						searchFailed = searchSeq.on_;
					}

					if(!balnFailed && !searchFailed){
						if (bAln.IsFirstMate()) {
							pairWriter.openWrite(PairedRead(balnSeq, searchSeq, false));
						} else {
							pairWriter.openWrite(PairedRead(searchSeq, balnSeq, false));
						}
					} else if (!balnFailed) {
						writer.openWrite(balnSeq);
					} else if (!searchFailed) {
						writer.openWrite(searchSeq);
					}
					// now that operations have been computed, remove first mate found from cache
					alnCache.remove(search->Name);
					continue;
				}
			}
			if (bAln.MatePosition == bAln.Position) {
				if (!alnCache.has(bAln.Name)) {
					//if mapped to the same place and the mate is yet to be encountered
					//enter into cache for until mate is encountered
					alnCache.add(bAln);
					continue;
				}
			}
			if (bAln.MatePosition <= bAln.Position) {
				if (!alnCache.has(bAln.Name)) {
					//since input should be sorted if matePosition is less than this position
					//it should be in the cache therefore program was unable to find mate
					//do orphaned operation
					auto balnSeq = bamAlnToSeqInfo(bAln);
					bool balnFailed = false;
					if(len(balnSeq) < lenCutOff){
						balnFailed = true;
					}else{
						qualChecker.checkRead(balnSeq);
						balnFailed = balnSeq.on_;
					}
					if(!balnFailed){
						writer.openWrite(balnSeq);
					}
					continue;
				} else {
					auto search = alnCache.get(bAln.Name);
					auto balnSeq = bamAlnToSeqInfo(bAln);
					auto searchSeq = bamAlnToSeqInfo(*search);


					bool balnFailed = false;
					if(len(balnSeq) < lenCutOff){
						balnFailed = true;
					}else{
						qualChecker.checkRead(balnSeq);
						balnFailed = balnSeq.on_;
					}

					bool searchFailed = false;
					if(len(searchSeq) < lenCutOff){
						searchFailed = true;
					}else{
						qualChecker.checkRead(searchSeq);
						searchFailed = searchSeq.on_;
					}

					if(!balnFailed && !searchFailed){
						if (bAln.IsFirstMate()) {
							pairWriter.openWrite(PairedRead(balnSeq, searchSeq, false));
						} else {
							pairWriter.openWrite(PairedRead(searchSeq, balnSeq, false));
						}
					} else if (!balnFailed) {
						writer.openWrite(balnSeq);
					} else if (!searchFailed) {
						writer.openWrite(searchSeq);
					}
					// now that operations have been computed, remove first mate found from cache
					alnCache.remove(search->Name);
				}
			} else {
				//enter into cache for until mate is encountered
				alnCache.add(bAln);
			}
		} else {
			//unpaired read
			if (!bAln.IsMapped()) {
				// do non-mapping operation
				auto balnSeq = bamAlnToSeqInfo(bAln);
				bool balnFailed = false;
				if(len(balnSeq) < lenCutOff){
					balnFailed = true;
				}else{
					qualChecker.checkRead(balnSeq);
					balnFailed = balnSeq.on_;
				}
				if(!balnFailed){
					writer.openWrite(balnSeq);
				}
			} else {
				//do unpaired read operation
				auto balnSeq = bamAlnToSeqInfo(bAln);
				bool balnFailed = false;
				if(len(balnSeq) < lenCutOff){
					balnFailed = true;
				}else{
					qualChecker.checkRead(balnSeq);
					balnFailed = balnSeq.on_;
				}
				if(!balnFailed){
					writer.openWrite(balnSeq);
				}
			}
		}
	}
	//save the orphans;
	if (len(alnCache) > 0) {
		auto names = alnCache.getNames();
		for (const auto & name : names) {
			auto search = alnCache.get(name);
			auto searchSeq = bamAlnToSeqInfo(*search);

			bool searchFailed = false;
			if(len(searchSeq) < lenCutOff){
				searchFailed = true;
			}else{
				qualChecker.checkRead(searchSeq);
				searchFailed = searchSeq.on_;
			}
			if(!searchFailed){
				writer.openWrite(searchSeq);
			}
			alnCache.remove(name);
		}
	}

	return 0;
}


}  //namespace njhseq
