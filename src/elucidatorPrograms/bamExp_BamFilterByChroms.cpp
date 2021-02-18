/*
 * bamExp_BamFilterByChroms.cpp
 *
 *  Created on: Jul 21, 2018
 *      Author: nick
 */
// elucidator - A library for analyzing sequence data
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
//
// This file is part of elucidator.
//
// elucidator is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// elucidator is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with elucidator.  If not, see <http://www.gnu.org/licenses/>.
//



#include "bamExp.hpp"
#include <TwoBit.h>

#include "elucidator/BamToolsUtils.h"
#include <njhseq/objects/BioDataObject.h>


namespace njhseq {


uint32_t getTotalSoftClippedBases(const BamTools::BamAlignment & baln){
	uint32_t ret = 0;
	if('S' == baln.CigarData.front().Type){
		ret += baln.CigarData.front().Length;
	}
	if(baln.CigarData.size() > 1 && 'S' == baln.CigarData.back().Type){
		ret += baln.CigarData.front().Length;
	}
	return ret;
}


int bamExpRunner::BamFilterByChromsToBam(const njh::progutils::CmdArgs & inputCommands){
	std::string chromFnp = "";
	bool writeFilteredBam = false;
	OutOptions outOpts(bfs::path("out"), ".bam");
	uint32_t allowableSoftClip = std::numeric_limits<uint32_t>::max();
	bool requireProperPair = false;
	bool skipWritingCounts = false;
	bool writeOnlyFilteredBam = false;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(allowableSoftClip, "--allowableSoftClip", "Number of bases that can be soft clipped in order to be included in the filtered off sequences, keep this zero to be more conservative in what gets filtered");
	setUp.setOption(chromFnp, "--chroms", "chromosomes to filter off", true);
	setUp.setOption(requireProperPair, "--requireProperPair", "Require Proper Pair to be filtered off");
	setUp.setOption(writeFilteredBam, "--writeFilteredBam", "Write Filtered Bam");
	setUp.setOption(writeOnlyFilteredBam, "--writeOnlyFilteredBam", "Write Only Filtered Bam");

	setUp.setOption(skipWritingCounts, "--skipWritingCounts", "Skip Writing Counts");
	setUp.processReadInNames({"--bam"}, true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	auto chroms = getInputValues(chromFnp, ",");

	BamTools::BamReader bReader;
	bReader.Open(setUp.pars_.ioOptions_.firstName_.string());
	checkBamOpenThrow(bReader, setUp.pars_.ioOptions_.firstName_.string());
	auto refData = bReader.GetReferenceData();
	BamTools::BamAlignment bAln;

	BamTools::BamWriter bWriter;
	BamTools::BamWriter bWriterFiltered;
	outOpts.throwIfOutExistsNoOverWrite(__PRETTY_FUNCTION__);

	bfs::path bamOut = outOpts.outName();
	bfs::path bamFilterOut = njh::files::prependFileBasename(outOpts.outName(), "filtered_");

	if(!writeOnlyFilteredBam){
		bWriter.Open(bamOut.string(), bReader.GetConstSamHeader(), refData);
	}
	if(writeFilteredBam || writeOnlyFilteredBam){
		bWriterFiltered.Open(bamFilterOut.string(), bReader.GetConstSamHeader(), refData);
	}


	auto totalsCountsOpts = OutOptions(njh::files::make_path(outOpts.outFilename_.string() + "_totalReadCounts"), ".tab.txt");
	auto filteredChromCountsOpts = OutOptions(njh::files::make_path(outOpts.outFilename_.string() + "_filteredByChrom"), ".tab.txt");

	totalsCountsOpts.transferOverwriteOpts(outOpts);
	filteredChromCountsOpts.transferOverwriteOpts(outOpts);


	std::shared_ptr<OutputStream> totalsCountsOut;
	std::shared_ptr<OutputStream>  filteredCountsOut;
	if(!skipWritingCounts){
		totalsCountsOut = std::make_shared<OutputStream>(totalsCountsOpts);
		filteredCountsOut = std::make_shared<OutputStream>(filteredChromCountsOpts);
	}



	//check to make sure chroms contains chromosome from the input bam file
	VecStr missing;
	for(const auto & chrom : chroms){
		bool found = false;
		for(const auto & ref : refData){
			if(ref.RefName == chrom){
				found = true;
				break;
			}
		}
		if(!found){
			missing.emplace_back(chrom);
		}
	}
	if(!missing.empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "the following chromosomes were not found in the input bam file " << njh::conToStr(missing)<< "\n";
		throw std::runtime_error{ss.str()};
	}

	struct ReadCounts{
		uint64_t singles_ = 0;
		uint64_t pairs_ = 0;
	};


	uint64_t filteredOrphans_ = 0;
	uint64_t keptOrphans_ = 0;


	ReadCounts input;
	ReadCounts kept;
	BamAlnsCache alnCache;
	BamAlnsCache filterAlnCache;

	std::unordered_map<std::string, ReadCounts> filteredCountsByChrom;

	while (bReader.GetNextAlignment(bAln)) {
		//skip secondary alignments
		if (!bAln.IsPrimaryAlignment()) {
			continue;
		}
		if (!bAln.IsPaired()) {
			++input.singles_;
			if (!bAln.IsMapped()) {
				++kept.singles_;
				if(!writeOnlyFilteredBam){
					bWriter.SaveAlignment(bAln);
				}
			} else {
				if (njh::in(refData[bAln.RefID].RefName, chroms)) {
					if(getTotalSoftClippedBases(bAln)   <= allowableSoftClip){
						if(writeFilteredBam || writeOnlyFilteredBam){
							bWriterFiltered.SaveAlignment(bAln);
						}
						++filteredCountsByChrom[refData[bAln.RefID].RefName].singles_;
					}else{
						++kept.singles_;
						if(!writeOnlyFilteredBam){
							bWriter.SaveAlignment(bAln);
						}
					}
				} else {
					++kept.singles_;
					if(!writeOnlyFilteredBam){
						bWriter.SaveAlignment(bAln);
					}
				}
			}
		} else {
			++input.pairs_;
			if (bAln.IsMapped() &&
					bAln.IsMateMapped() &&
					njh::in(refData[bAln.RefID].RefName, chroms) &&
					njh::in(refData[bAln.MateRefID].RefName, chroms) &&
					(!requireProperPair || bAln.IsProperPair())){
				if (!filterAlnCache.has(bAln.Name)) {
					//pair hasn't been added to cache yet so add to cache
					//this only works if mate and first mate have the same name
					filterAlnCache.add(bAln);
					continue;
				} else {
					auto search = filterAlnCache.get(bAln.Name);
					if(getTotalSoftClippedBases(*search) <= allowableSoftClip &&
							getTotalSoftClippedBases(bAln)   <= allowableSoftClip){
						++filteredCountsByChrom[njh::pasteAsStr(refData[search->RefID].RefName, "--", refData[bAln.RefID].RefName)].pairs_;
						++filteredCountsByChrom[njh::pasteAsStr(refData[search->RefID].RefName, "--", refData[bAln.RefID].RefName)].pairs_;
						if(writeFilteredBam || writeOnlyFilteredBam){
							bWriterFiltered.SaveAlignment(bAln);
							bWriterFiltered.SaveAlignment(*search);
						}
					}else{
						++kept.pairs_;++kept.pairs_;
						if(!writeOnlyFilteredBam){
							bWriter.SaveAlignment(bAln);
							bWriter.SaveAlignment(*search);
						}
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
					++kept.pairs_;++kept.pairs_;
					if(!writeOnlyFilteredBam){
						bWriter.SaveAlignment(bAln);
						bWriter.SaveAlignment(*search);
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
			++keptOrphans_;
			auto search = alnCache.get(name);
			if(!writeOnlyFilteredBam){
				bWriter.SaveAlignment(*search);
			}
			alnCache.remove(name);
		}
	}

	if (len(filterAlnCache) > 0) {
		auto names = filterAlnCache.getNames();
		for (const auto &name : names) {
			++filteredOrphans_;
			auto search = filterAlnCache.get(name);
			if (writeFilteredBam || writeOnlyFilteredBam) {
				bWriterFiltered.SaveAlignment(*search);
			}
			filterAlnCache.remove(name);
		}
	}


	ReadCounts filtered;
	for(const auto & filt : filteredCountsByChrom){
		filtered.pairs_ += filt.second.pairs_;
		filtered.singles_ += filt.second.singles_;
	}


	if(!skipWritingCounts){
		*totalsCountsOut << "bam\tcondition\tcount\tfrac\ttotal" << std::endl;;
		*totalsCountsOut << bfs::basename(setUp.pars_.ioOptions_.firstName_.filename())
				<< "\t" << "keptPairs"
				<< "\t" << kept.pairs_
				<< "\t" << kept.pairs_/static_cast<double>(input.pairs_)
				<< "\t" << input.pairs_ << std::endl;
		*totalsCountsOut << bfs::basename(setUp.pars_.ioOptions_.firstName_.filename())
				<< "\t" << "keptSingles"
				<< "\t" << kept.singles_
				<< "\t" << kept.singles_/static_cast<double>(input.singles_)
				<< "\t" << input.singles_ << std::endl;

		*totalsCountsOut << bfs::basename(setUp.pars_.ioOptions_.firstName_.filename())
				<< "\t" << "filteredPairs"
				<< "\t" << filtered.pairs_
				<< "\t" << filtered.pairs_/static_cast<double>(input.pairs_)
				<< "\t" << input.pairs_ << std::endl;
		*totalsCountsOut << bfs::basename(setUp.pars_.ioOptions_.firstName_.filename())
				<< "\t" << "filteredSingles"
				<< "\t" << filtered.singles_
				<< "\t" << filtered.singles_/static_cast<double>(input.singles_)
				<< "\t" << input.singles_ << std::endl;

		*totalsCountsOut << bfs::basename(setUp.pars_.ioOptions_.firstName_.filename())
				<< "\t" << "keptOrphans"
				<< "\t" << keptOrphans_
				<< "\t"
				<< "\t" << std::endl;

		*totalsCountsOut << bfs::basename(setUp.pars_.ioOptions_.firstName_.filename())
				<< "\t" << "filteredOrphans"
				<< "\t" << filteredOrphans_
				<< "\t"
				<< "\t" << std::endl;

		auto names = getVectorOfMapKeys(filteredCountsByChrom);
		njh::sort(names);
		*filteredCountsOut << "bam\tchrom\tpairs\tpairsFrac\tsingles\tsinglesFrac" << std::endl;
		for(const auto & name : names){
			*filteredCountsOut << bfs::basename(setUp.pars_.ioOptions_.firstName_.filename())
					<< "\t" << name
					<< "\t" << filteredCountsByChrom[name].pairs_
					<< "\t" << filteredCountsByChrom[name].pairs_/static_cast<double>(filtered.pairs_)
					<< "\t" << filteredCountsByChrom[name].singles_
					<< "\t" << filteredCountsByChrom[name].singles_/static_cast<double>(filtered.singles_) << std::endl;

		}
	}

	return 0;
}



int bamExpRunner::BamGetImproperPairsOnChroms(const njh::progutils::CmdArgs & inputCommands){
	std::string chromFnp = "";
	OutOptions outOpts(bfs::path("out"));
	uint32_t allowableSoftClip = std::numeric_limits<uint32_t>::max();
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(allowableSoftClip, "--allowableSoftClip", "Number of bases that can be soft clipped in order to be included in the filtered off sequences, keep this zero to be more conservative in what gets filtered");
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


	singlesOpts.out_.transferOverwriteOpts(outOpts);
	pairedOpts.out_.transferOverwriteOpts(outOpts);


	SeqOutput singlesWriter(singlesOpts);
	SeqOutput pairedWriter(pairedOpts);


	BamAlnsCache alnCache;
	BamAlnsCache filterAlnCache;

	auto refData = bReader.GetReferenceData();
	//check to make sure chroms contains chromosome from the input bam file
	VecStr missing;
	for(const auto & chrom : chroms){
		bool found = false;
		for(const auto & ref : refData){
			if(ref.RefName == chrom){
				found = true;
				break;
			}
		}
		if(!found){
			missing.emplace_back(chrom);
		}
	}
	if(!missing.empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "the following chromosomes were not found in the input bam file " << njh::conToStr(missing)<< "\n";
		throw std::runtime_error{ss.str()};
	}
	while (bReader.GetNextAlignment(bAln)) {
		//skip secondary alignments
		if (!bAln.IsPrimaryAlignment()) {
			continue;
		}
		if (bAln.IsPaired()) {
			if (bAln.IsMapped() &&
					bAln.IsMateMapped() &&
					njh::in(refData[bAln.RefID].RefName, chroms) &&
					njh::in(refData[bAln.MateRefID].RefName, chroms) &&
					!bAln.IsProperPair()){
				if (!filterAlnCache.has(bAln.Name)) {
					//pair hasn't been added to cache yet so add to cache
					//this only works if mate and first mate have the same name
					filterAlnCache.add(bAln);
					continue;
				} else {
					auto search = filterAlnCache.get(bAln.Name);
					if(getTotalSoftClippedBases(*search) <= allowableSoftClip &&
							getTotalSoftClippedBases(bAln)   <= allowableSoftClip){
						if (bAln.IsFirstMate()) {
							pairedWriter.openWrite(PairedRead(bamAlnToSeqInfo(bAln), bamAlnToSeqInfo(*search),false));
						} else {
							pairedWriter.openWrite(PairedRead(bamAlnToSeqInfo(*search), bamAlnToSeqInfo(bAln),false));
						}
					}
					// now that operations have been computed, remove ther other mate found from cache
					filterAlnCache.remove(search->Name);
					continue;
				}
			}
		}
	}
	return 0;
}


int bamExpRunner::BamFilterByChroms(const njh::progutils::CmdArgs & inputCommands){
	std::string chromFnp = "";
	OutOptions outOpts(bfs::path("out"));
	uint32_t allowableSoftClip = 10;
	bool requireProperPair = false;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(allowableSoftClip, "--allowableSoftClip", "Number of bases that can be soft clipped in order to be included in the filtered off sequences, keep this zero to be more conservative in what gets filtered");
	setUp.setOption(chromFnp, "--chroms", "chromosomes to filter off", true);
	setUp.setOption(requireProperPair, "--requireProperPair", "Require Proper Pair to be filtered off");
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

	auto totalsCountsOpts = OutOptions(njh::files::make_path(outOpts.outFilename_.string() + "_totalReadCounts"), ".tab.txt");
	auto filteredChromCountsOpts = OutOptions(njh::files::make_path(outOpts.outFilename_.string() + "_filteredByChrom"), ".tab.txt");

	totalsCountsOpts.transferOverwriteOpts(outOpts);
	filteredChromCountsOpts.transferOverwriteOpts(outOpts);
	singlesOpts.out_.transferOverwriteOpts(outOpts);
	pairedOpts.out_.transferOverwriteOpts(outOpts);
	filteredSinglesOpts.out_.transferOverwriteOpts(outOpts);
	filteredPairedOpts.out_.transferOverwriteOpts(outOpts);

	OutputStream totalsCountsOut(totalsCountsOpts);
	OutputStream filteredCountsOut(filteredChromCountsOpts);


	SeqOutput singlesWriter(singlesOpts);
	SeqOutput pairedWriter(pairedOpts);
	SeqOutput filteredSinglesWriter(filteredSinglesOpts);
	SeqOutput filteredPairedWriter(filteredPairedOpts);

	BamAlnsCache alnCache;
	BamAlnsCache filterAlnCache;

	auto refData = bReader.GetReferenceData();
	//check to make sure chroms contains chromosome from the input bam file
	VecStr missing;
	for(const auto & chrom : chroms){
		bool found = false;
		for(const auto & ref : refData){
			if(ref.RefName == chrom){
				found = true;
				break;
			}
		}
		if(!found){
			missing.emplace_back(chrom);
		}
	}
	if(!missing.empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "the following chromosomes were not found in the input bam file " << njh::conToStr(missing)<< "\n";
		throw std::runtime_error{ss.str()};
	}

	struct ReadCounts{
		uint64_t singles_ = 0;
		uint64_t pairs_ = 0;
	};


	uint64_t filteredOrphans_ = 0;
	uint64_t keptOrphans_ = 0;


	ReadCounts input;
	ReadCounts kept;

	std::unordered_map<std::string, ReadCounts> filteredCountsByChrom;

	while (bReader.GetNextAlignment(bAln)) {
		//skip secondary alignments
		if (!bAln.IsPrimaryAlignment()) {
			continue;
		}
		if (!bAln.IsPaired()) {
			++input.singles_;
			if (!bAln.IsMapped()) {
				++kept.singles_;
				singlesWriter.openWrite(bamAlnToSeqInfo(bAln));
			} else {
				if (njh::in(refData[bAln.RefID].RefName, chroms)) {
					if(getTotalSoftClippedBases(bAln)   <= allowableSoftClip){
						filteredSinglesWriter.openWrite(bamAlnToSeqInfo(bAln));
						++filteredCountsByChrom[refData[bAln.RefID].RefName].singles_;
					}else{
						++kept.singles_;
						singlesWriter.openWrite(bamAlnToSeqInfo(bAln));
					}
				} else {
					++kept.singles_;
					singlesWriter.openWrite(bamAlnToSeqInfo(bAln));
				}
			}
		} else {
			++input.pairs_;
			if (bAln.IsMapped() &&
					bAln.IsMateMapped() &&
					njh::in(refData[bAln.RefID].RefName, chroms) &&
					njh::in(refData[bAln.MateRefID].RefName, chroms) &&
					(!requireProperPair || bAln.IsProperPair())){
				if (!filterAlnCache.has(bAln.Name)) {
					//pair hasn't been added to cache yet so add to cache
					//this only works if mate and first mate have the same name
					filterAlnCache.add(bAln);
					continue;
				} else {
					auto search = filterAlnCache.get(bAln.Name);
					if(getTotalSoftClippedBases(*search) <= allowableSoftClip &&
							getTotalSoftClippedBases(bAln)   <= allowableSoftClip){
						++filteredCountsByChrom[njh::pasteAsStr(refData[search->RefID].RefName, "--", refData[bAln.RefID].RefName)].pairs_;
						++filteredCountsByChrom[njh::pasteAsStr(refData[search->RefID].RefName, "--", refData[bAln.RefID].RefName)].pairs_;
						if (bAln.IsFirstMate()) {
							filteredPairedWriter.openWrite(PairedRead(bamAlnToSeqInfo(bAln), bamAlnToSeqInfo(*search),false));
						} else {
							filteredPairedWriter.openWrite(PairedRead(bamAlnToSeqInfo(*search), bamAlnToSeqInfo(bAln),false));
						}
					}else{
						++kept.pairs_;++kept.pairs_;
						if (bAln.IsFirstMate()) {
							pairedWriter.openWrite(PairedRead(bamAlnToSeqInfo(bAln), bamAlnToSeqInfo(*search),false));
						} else {
							pairedWriter.openWrite(PairedRead(bamAlnToSeqInfo(*search), bamAlnToSeqInfo(bAln),false));
						}
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
					++kept.pairs_;++kept.pairs_;
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
			++keptOrphans_;
			auto search = alnCache.get(name);
			singlesWriter.openWrite(bamAlnToSeqInfo(*search));
			alnCache.remove(name);
		}
	}
	if (len(filterAlnCache) > 0) {
		auto names = filterAlnCache.getNames();
		for (const auto & name : names) {
			++filteredOrphans_;
			auto search = filterAlnCache.get(name);
			filteredSinglesWriter.openWrite(bamAlnToSeqInfo(*search));
			filterAlnCache.remove(name);
		}
	}

	ReadCounts filtered;
	for(const auto & filt : filteredCountsByChrom){
		filtered.pairs_ += filt.second.pairs_;
		filtered.singles_ += filt.second.singles_;
	}



	totalsCountsOut << "bam\tcondition\tcount\tfrac\ttotal" << std::endl;;
	totalsCountsOut << bfs::basename(setUp.pars_.ioOptions_.firstName_.filename())
			<< "\t" << "keptPairs"
			<< "\t" << kept.pairs_
			<< "\t" << kept.pairs_/static_cast<double>(input.pairs_)
			<< "\t" << input.pairs_ << std::endl;
	totalsCountsOut << bfs::basename(setUp.pars_.ioOptions_.firstName_.filename())
			<< "\t" << "keptSingles"
			<< "\t" << kept.singles_
			<< "\t" << kept.singles_/static_cast<double>(input.singles_)
			<< "\t" << input.singles_ << std::endl;

	totalsCountsOut << bfs::basename(setUp.pars_.ioOptions_.firstName_.filename())
			<< "\t" << "filteredPairs"
			<< "\t" << filtered.pairs_
			<< "\t" << filtered.pairs_/static_cast<double>(input.pairs_)
			<< "\t" << input.pairs_ << std::endl;
	totalsCountsOut << bfs::basename(setUp.pars_.ioOptions_.firstName_.filename())
			<< "\t" << "filteredSingles"
			<< "\t" << filtered.singles_
			<< "\t" << filtered.singles_/static_cast<double>(input.singles_)
			<< "\t" << input.singles_ << std::endl;

	totalsCountsOut << bfs::basename(setUp.pars_.ioOptions_.firstName_.filename())
			<< "\t" << "keptOrphans"
			<< "\t" << keptOrphans_
			<< "\t"
			<< "\t" << std::endl;

	totalsCountsOut << bfs::basename(setUp.pars_.ioOptions_.firstName_.filename())
			<< "\t" << "filteredOrphans"
			<< "\t" << filteredOrphans_
			<< "\t"
			<< "\t" << std::endl;

	auto names = getVectorOfMapKeys(filteredCountsByChrom);
	njh::sort(names);
	filteredCountsOut << "bam\tchrom\tpairs\tpairsFrac\tsingles\tsinglesFrac" << std::endl;
	for(const auto & name : names){
		filteredCountsOut << bfs::basename(setUp.pars_.ioOptions_.firstName_.filename())
				<< "\t" << name
				<< "\t" << filteredCountsByChrom[name].pairs_
				<< "\t" << filteredCountsByChrom[name].pairs_/static_cast<double>(filtered.pairs_)
				<< "\t" << filteredCountsByChrom[name].singles_
				<< "\t" << filteredCountsByChrom[name].singles_/static_cast<double>(filtered.singles_) << std::endl;

	}

	return 0;
}


} // namespace njhseq

