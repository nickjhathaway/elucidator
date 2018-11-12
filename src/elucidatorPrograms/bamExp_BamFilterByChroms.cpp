/*
 * bamExp_BamFilterByChroms.cpp
 *
 *  Created on: Jul 21, 2018
 *      Author: nick
 */




#include "bamExp.hpp"
#include <TwoBit.h>

#include "elucidator/BamToolsUtils.h"
#include <njhseq/objects/BioDataObject.h>


namespace njhseq {


int bamExpRunner::BamFilterByChroms(const njh::progutils::CmdArgs & inputCommands){
	std::string chromFnp = "";
	OutOptions outOpts(bfs::path("out"));
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.setOption(chromFnp, "--chroms", "chromosomes to filter", true);
	setUp.processReadInNames({"--bam"}, true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	auto chroms = getInputValues(chromFnp, ",");


	BamTools::BamReader bReader;
	bReader.Open(setUp.pars_.ioOptions_.firstName_.string());
	checkBamOpenThrow(bReader, setUp.pars_.ioOptions_.firstName_.string());

	BamTools::BamAlignment bAln;

	auto singlesOpts = SeqIOOptions::genFastqOutGz(outOpts.outFilename_.string() + "_singles");
	auto pairedOpts = SeqIOOptions::genPairedOutGz(outOpts.outFilename_.string() + "_pairs");
	auto filteredSinglesOpts = SeqIOOptions::genFastqOutGz(outOpts.outFilename_.string() + "_filteredOffSingles");
	auto filteredPairedOpts = SeqIOOptions::genPairedOutGz(outOpts.outFilename_.string() + "_filteredOffPairs");

	singlesOpts.out_.transferOverwriteOpts(outOpts);
	pairedOpts.out_.transferOverwriteOpts(outOpts);
	filteredSinglesOpts.out_.transferOverwriteOpts(outOpts);
	filteredPairedOpts.out_.transferOverwriteOpts(outOpts);

	SeqOutput singlesWriter(singlesOpts);
	SeqOutput pairedWriter(pairedOpts);
	SeqOutput filteredSinglesWriter(filteredSinglesOpts);
	SeqOutput filteredPairedWriter(filteredPairedOpts);

	BamAlnsCache alnCache;
	BamAlnsCache filterAlnCache;

	auto refData = bReader.GetReferenceData();
	while (bReader.GetNextAlignment(bAln)) {
		//skip secondary alignments
		if (!bAln.IsPrimaryAlignment()) {
			continue;
		}
		if (!bAln.IsPaired()) {
			if (!bAln.IsMapped()) {
				singlesWriter.openWrite(bamAlnToSeqInfo(bAln));
			} else {
				if (njh::in(refData[bAln.RefID].RefName, chroms)) {
					filteredSinglesWriter.openWrite(bamAlnToSeqInfo(bAln));
				} else {
					singlesWriter.openWrite(bamAlnToSeqInfo(bAln));
				}
			}
		} else {
			if (bAln.IsMapped() &&
					bAln.IsMateMapped() &&
					njh::in(refData[bAln.RefID].RefName, chroms) &&
					njh::in(refData[bAln.MateRefID].RefName, chroms)){
				if (!filterAlnCache.has(bAln.Name)) {
					//pair hasn't been added to cache yet so add to cache
					//this only works if mate and first mate have the same name
					filterAlnCache.add(bAln);
					continue;
				} else {
					auto search = filterAlnCache.get(bAln.Name);
					if (bAln.IsFirstMate()) {
						filteredPairedWriter.openWrite(PairedRead(bamAlnToSeqInfo(bAln), bamAlnToSeqInfo(*search),false));
					} else {
						filteredPairedWriter.openWrite(PairedRead(bamAlnToSeqInfo(*search), bamAlnToSeqInfo(bAln),false));
					}
					// now that operations have been computed, remove ther other mate found from cache
					filterAlnCache.remove(search->Name);
					continue;
				}
			}else{
				if (!alnCache.has(bAln.Name)) {
					//pair hasn't been added to cache yet so add to cache
					//this only works if mate and first mate have the same name
					alnCache.add(bAln);
					continue;
				} else {
					auto search = alnCache.get(bAln.Name);
					if (bAln.IsFirstMate()) {
						pairedWriter.openWrite(PairedRead(bamAlnToSeqInfo(bAln), bamAlnToSeqInfo(*search),false));
					} else {
						pairedWriter.openWrite(PairedRead(bamAlnToSeqInfo(*search), bamAlnToSeqInfo(bAln),false));
					}
					// now that operations have been computed, remove ther other mate found from cache
					alnCache.remove(search->Name);
					continue;
				}
			}
		}
	}

	//save the orphans;
	if (len(alnCache) > 0) {
		auto names = alnCache.getNames();
		for (const auto & name : names) {
			auto search = alnCache.get(name);
			singlesWriter.openWrite(bamAlnToSeqInfo(*search));
			alnCache.remove(name);
		}
	}
	if (len(filterAlnCache) > 0) {
		auto names = filterAlnCache.getNames();
		for (const auto & name : names) {
			auto search = filterAlnCache.get(name);
			filteredSinglesWriter.openWrite(bamAlnToSeqInfo(*search));
			filterAlnCache.remove(name);
		}
	}
	return 0;
}


} // namespace njhseq

