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





int bamExpRunner::BamFilterByChromsToBam(const njh::progutils::CmdArgs & inputCommands){
	std::string chromFnp;
	bool writeFilteredBam = false;
	OutOptions outOpts(bfs::path("out"), ".bam");
	uint32_t allowableSoftClipInAln = std::numeric_limits<uint32_t>::max();
	bool requireProperPair = false;
	bool skipWritingCounts = false;
	bool writeOnlyFilteredBam = false;
	bool any_mate = false;
	bool writeOutUnmappedSeparately = false;
	bool doNotWriteFilterOff = false;
	uint32_t minMappingQuality = 0;
	bool filterWithUnmappedMate = false;


	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(allowableSoftClipInAln, "--allowableSoftClipInAln", "Number of bases that can be soft clipped in order to be included in the filtered off sequences, keep this zero to be more conservative in what gets filtered");
	setUp.setOption(chromFnp, "--chroms", "chromosomes to filter off", true);
	setUp.setOption(requireProperPair, "--requireProperPair", "Require Proper Pair to be filtered off");
	setUp.setOption(any_mate, "--any_mate", "filter off if even one mate maps to the filter chroms");
	setUp.setOption(doNotWriteFilterOff, "--doNotWriteFilterOff", "do Not Write Filter Off");
	setUp.setOption(minMappingQuality, "--minMappingQuality", "min Mapping Quality");
	setUp.setOption(filterWithUnmappedMate, "--filterWithUnmappedMate", "by default requires both mates to map to a filter chromosome, this will filter if one mate maps and the other is unmapped");

	setUp.setOption(requireProperPair, "--requireProperPair", "Require Proper Pair to be filtered off");
	setUp.setOption(writeFilteredBam, "--writeFilteredBam", "Write Filtered Bam");
	setUp.setOption(writeOnlyFilteredBam, "--writeOnlyFilteredBam", "Write Only Filtered Bam");
	setUp.setOption(writeOutUnmappedSeparately, "--writeOutUnmappedSeparately", "write Out Unmapped Separately");


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
	BamTools::BamWriter bWriterUnmapped;

	outOpts.throwIfOutExistsNoOverWrite(__PRETTY_FUNCTION__);

	bfs::path bamOut = outOpts.outName();
	bfs::path bamFilterOut = njh::files::prependFileBasename(outOpts.outName(), "filtered_");
	bfs::path bamUnmappedOut = njh::files::prependFileBasename(outOpts.outName(), "unmapped_");


	if(!writeOnlyFilteredBam){
		bWriter.Open(bamOut.string(), bReader.GetConstSamHeader(), refData);
	}
	if (writeOutUnmappedSeparately) {
		bWriterUnmapped.Open(bamUnmappedOut.string(), bReader.GetConstSamHeader(), refData);
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
	uint64_t unmappedOrphans_ = 0;

	auto doesAlnPassSoftClipFilt = [&allowableSoftClipInAln,&refData](const BamTools::BamAlignment & bamAln){
		uint32_t softClipSum = 0;
		if(!bamAln.CigarData.empty() && 'S' == bamAln.CigarData.front().Type && 0 != bamAln.Position){
			softClipSum += bamAln.CigarData.front().Length;
		}
		if(bamAln.CigarData.size() > 1 && 'S' == bamAln.CigarData.back().Type && bamAln.GetEndPosition() != refData[bamAln.RefID].RefLength){
			softClipSum += bamAln.CigarData.back().Length;
		}
		return softClipSum <= allowableSoftClipInAln;
	};

	BamAlnsCache alnCache;
	BamAlnsCache filterAlnCache;
	ReadCounts input;
	ReadCounts kept;
	ReadCounts unmapped;

	std::unordered_map<std::string, ReadCounts> filteredCountsByChrom;

	while (bReader.GetNextAlignment(bAln)) {
		//skip secondary alignments
		if (!bAln.IsPrimaryAlignment()) {
			continue;
		}
		if (!bAln.IsPaired()) {
			++input.singles_;
			if (!bAln.IsMapped()) {
				if (writeOutUnmappedSeparately) {
					++unmapped.singles_;
					bWriterUnmapped.SaveAlignment(bAln);
				} else {
					++kept.singles_;
					bWriter.SaveAlignment(bAln);
				}
			} else {
				if (njh::in(refData[bAln.RefID].RefName, chroms)) {
					if(doesAlnPassSoftClipFilt(bAln) && bAln.MapQuality >= minMappingQuality){
						if(!doNotWriteFilterOff){
							bWriterFiltered.SaveAlignment(bAln);
						}
						++filteredCountsByChrom[refData[bAln.RefID].RefName].singles_;
					} else {
						if (writeOutUnmappedSeparately) {
							++unmapped.singles_;
							bWriterUnmapped.SaveAlignment(bAln);
						} else {
							++kept.singles_;
							bWriter.SaveAlignment(bAln);
						}
					}
				} else {
					++kept.singles_;
					bWriter.SaveAlignment(bAln);
				}
			}
		} else {
			++input.pairs_;
			if (bAln.IsMapped() &&
					bAln.IsMateMapped() &&
					((any_mate && (njh::in(refData[bAln.RefID].RefName, chroms) ||
					njh::in(refData[bAln.MateRefID].RefName, chroms))) || (njh::in(refData[bAln.RefID].RefName, chroms) &&
					njh::in(refData[bAln.MateRefID].RefName, chroms))) &&
					(!requireProperPair || bAln.IsProperPair())){
				if (!filterAlnCache.has(bAln.Name)) {
					//pair hasn't been added to cache yet so add to cache
					//this only works if mate and first mate have the same name
					filterAlnCache.add(bAln);
				} else {
					auto search = filterAlnCache.get(bAln.Name);
					bool baln_pass = true;
					bool search_pass = true;
					if (njh::in(refData[bAln.RefID].RefName, chroms)) {
						baln_pass = !(doesAlnPassSoftClipFilt(bAln) && bAln.MapQuality >= minMappingQuality);
					}
					if (njh::in(refData[search->RefID].RefName, chroms)) {
						search_pass = !(doesAlnPassSoftClipFilt(*search) && search->MapQuality >= minMappingQuality);
					}
					if(baln_pass && search_pass){
						++filteredCountsByChrom[njh::pasteAsStr(refData[search->RefID].RefName, "--", refData[bAln.RefID].RefName)].pairs_;
						++filteredCountsByChrom[njh::pasteAsStr(refData[search->RefID].RefName, "--", refData[bAln.RefID].RefName)].pairs_;
						if(!doNotWriteFilterOff){
							bWriterFiltered.SaveAlignment(*search);
							bWriterFiltered.SaveAlignment(bAln);
						}
					} else {
						if (writeOutUnmappedSeparately) {
							++unmapped.pairs_;
							++unmapped.pairs_;
							bWriterUnmapped.SaveAlignment(*search);
							bWriterUnmapped.SaveAlignment(bAln);
						} else {
							++kept.pairs_;
							++kept.pairs_;
							bWriter.SaveAlignment(*search);
							bWriter.SaveAlignment(bAln);
						}
					}
					// now that operations have been computed, remove their other mate found from cache
					filterAlnCache.remove(search->Name);
				}
			}else{
				if (!alnCache.has(bAln.Name)) {
					//pair hasn't been added to cache yet so add to cache
					//this only works if mate and first mate have the same name
					alnCache.add(bAln);
				} else {
					auto search = alnCache.get(bAln.Name);
					if (filterWithUnmappedMate &&
					    (
						    (bAln.IsMapped() && !bAln.IsMateMapped() && njh::in(refData[bAln.RefID].RefName, chroms) &&
						     doesAlnPassSoftClipFilt(bAln) && bAln.MapQuality >= minMappingQuality) ||
						    (!bAln.IsMapped() && bAln.IsMateMapped() && njh::in(refData[bAln.MateRefID].RefName, chroms) &&
						     doesAlnPassSoftClipFilt(*search) && search->MapQuality >= minMappingQuality)
					    )
					) {
						std::string filterChromName;
						if (bAln.IsMapped() && !bAln.IsMateMapped()) {
							filterChromName = njh::pasteAsStr("unmapped", "--", refData[bAln.RefID].RefName);
						} else {
							filterChromName = njh::pasteAsStr(refData[search->RefID].RefName, "--", "unmapped");
						}
						++filteredCountsByChrom[filterChromName].pairs_;
						++filteredCountsByChrom[filterChromName].pairs_;
						if(!doNotWriteFilterOff){
							bWriterFiltered.SaveAlignment(*search);
							bWriterFiltered.SaveAlignment(bAln);
						}
					} else {
						if (writeOutUnmappedSeparately &&
							(
								(!bAln.IsMapped() && !bAln.IsMateMapped()) ||
							(
								( bAln.IsMapped() &&!bAln.IsMateMapped()  && njh::in(refData[bAln.RefID].RefName, chroms)) ||
							  (!bAln.IsMapped() && bAln.IsMateMapped()  && njh::in(refData[bAln.MateRefID].RefName, chroms))
							 )
							 )
							 ) {
							++unmapped.pairs_;
							++unmapped.pairs_;
							bWriterUnmapped.SaveAlignment(*search);
							bWriterUnmapped.SaveAlignment(bAln);
						} else {
							++kept.pairs_;
							++kept.pairs_;
							bWriter.SaveAlignment(*search);
							bWriter.SaveAlignment(bAln);
						}
					}
					// now that operations have been computed, remove ther other mate found from cache
					alnCache.remove(search->Name);
				}
			}
		}
	}

	//save the orphans;
	if (len(alnCache) > 0) {
		auto names = alnCache.getNames();
		for (const auto & name : names) {

			auto search = alnCache.get(name);
			if (writeOutUnmappedSeparately && !search->IsMapped()) {
				++unmappedOrphans_;
				bWriterUnmapped.SaveAlignment(*search);
			} else {
				++keptOrphans_;
				bWriter.SaveAlignment(*search);
			}
			alnCache.remove(name);
		}
	}
	if (len(filterAlnCache) > 0) {
		auto names = filterAlnCache.getNames();
		for (const auto & name : names) {
			auto search = filterAlnCache.get(name);
			if(doesAlnPassSoftClipFilt(*search) && search->MapQuality >= minMappingQuality) {
				++filteredOrphans_;
				if(!doNotWriteFilterOff){
					bWriterFiltered.SaveAlignment(*search);
				}
			} else {
				if (writeOutUnmappedSeparately) {
					++unmappedOrphans_;
					bWriterUnmapped.SaveAlignment(*search);
				} else {
					++keptOrphans_;
					bWriter.SaveAlignment(*search);
				}
			}
			filterAlnCache.remove(name);
		}
	}

	if (!skipWritingCounts) {
		ReadCounts filtered;
		for(const auto & filt : filteredCountsByChrom){
			filtered.pairs_ += filt.second.pairs_;
			filtered.singles_ += filt.second.singles_;
		}

		auto bname = bfs::basename(setUp.pars_.ioOptions_.firstName_.filename());
		*totalsCountsOut << "bam\tcondition\tcount\tfrac\ttotal" << std::endl;
		*totalsCountsOut << bname
				<< "\t" << "keptPairs"
				<< "\t" << kept.pairs_
				<< "\t" << kept.pairs_/static_cast<long double>(input.pairs_)
				<< "\t" << input.pairs_ << std::endl;
		*totalsCountsOut << bname
				<< "\t" << "keptSingles"
				<< "\t" << kept.singles_
				<< "\t" << kept.singles_/static_cast<long double>(input.singles_)
				<< "\t" << input.singles_ << std::endl;

		*totalsCountsOut << bname
				<< "\t" << "filteredPairs"
				<< "\t" << filtered.pairs_
				<< "\t" << filtered.pairs_/static_cast<long double>(input.pairs_)
				<< "\t" << input.pairs_ << std::endl;
		*totalsCountsOut << bname
				<< "\t" << "filteredSingles"
				<< "\t" << filtered.singles_
				<< "\t" << filtered.singles_/static_cast<long double>(input.singles_)
				<< "\t" << input.singles_ << std::endl;

		if (writeOutUnmappedSeparately) {
			*totalsCountsOut << bname
					<< "\t" << "unmappedPairs"
					<< "\t" << unmapped.pairs_
					<< "\t" << unmapped.pairs_ / static_cast<long double>(input.pairs_)
					<< "\t" << input.pairs_ << std::endl;
			*totalsCountsOut << bname
					<< "\t" << "unmappedSingles"
					<< "\t" << unmapped.singles_
					<< "\t" << unmapped.singles_ / static_cast<long double>(input.singles_)
					<< "\t" << input.singles_ << std::endl;
		}

		*totalsCountsOut << bname
				<< "\t" << "keptOrphans"
				<< "\t" << keptOrphans_
				<< "\t"
				<< "\t" << std::endl;

		*totalsCountsOut << bname
				<< "\t" << "filteredOrphans"
				<< "\t" << filteredOrphans_
				<< "\t"
				<< "\t" << std::endl;
		if (writeOutUnmappedSeparately) {
			*totalsCountsOut << bname
					<< "\t" << "unmappedOrphans"
					<< "\t" << unmappedOrphans_
					<< "\t"
					<< "\t" << std::endl;
		}

		auto names = getVectorOfMapKeys(filteredCountsByChrom);
		njh::sort(names);
		*filteredCountsOut << "bam\tchrom\tpairs\tpairsFrac\tsingles\tsinglesFrac" << std::endl;
		for(const auto & name : names){
			*filteredCountsOut << bname
					<< "\t" << name
					<< "\t" << filteredCountsByChrom[name].pairs_
					<< "\t" << filteredCountsByChrom[name].pairs_/static_cast<long double>(filtered.pairs_)
					<< "\t" << filteredCountsByChrom[name].singles_
					<< "\t" << filteredCountsByChrom[name].singles_/static_cast<long double>(filtered.singles_) << std::endl;
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
					if(getSoftClipAmount(*search) <= allowableSoftClip &&
							getSoftClipAmount(bAln)   <= allowableSoftClip){
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
	std::string chromFnp;
	OutOptions outOpts(bfs::path("out"));
	uint32_t minMappingQuality = 0;
	uint32_t allowableSoftClipInAln = std::numeric_limits<uint32_t>::max();
	bool any_mate = false;
	bool requireProperPair = false;
	bool doNotWriteFilterOff = false;
	bool filterWithUnmappedMate = false;
	bool writeOutUnmappedSeparately = false;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(allowableSoftClipInAln, "--allowableSoftClip", "Number of bases that can be soft clipped in order to be included in the filtered off sequences, keep this zero to be more conservative in what gets filtered");
	setUp.setOption(chromFnp, "--chroms", "chromosomes to filter off", true);
	setUp.setOption(any_mate, "--any_mate", "filter off if even one mate maps to the filter chroms");
	setUp.setOption(requireProperPair, "--requireProperPair", "Require Proper Pair to be filtered off");
	setUp.setOption(doNotWriteFilterOff, "--doNotWriteFilterOff", "do Not Write Filter Off");
	setUp.setOption(minMappingQuality, "--minMappingQuality", "min Mapping Quality");
	setUp.setOption(filterWithUnmappedMate, "--filterWithUnmappedMate", "by default requires both mates to map to a filter chromosome, this will filter if one mate maps and the other is unmapped");
	setUp.setOption(writeOutUnmappedSeparately, "--writeOutUnmappedSeparately", "write out completely unmapped sequences (both mates unmapped) to a different file from the kept ones");

	setUp.processReadInNames({"--bam"}, true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	auto filter_chroms = getInputValues(chromFnp, ",");

	BamTools::BamReader bReader;
	bReader.Open(setUp.pars_.ioOptions_.firstName_.string());
	checkBamOpenThrow(bReader, setUp.pars_.ioOptions_.firstName_.string());

	BamTools::BamAlignment bAln;

	auto singlesOpts = SeqIOOptions::genFastqOutGz(outOpts.outFilename_.string() + "_singles");
	auto pairedOpts = SeqIOOptions::genPairedOutGz(outOpts.outFilename_.string() + "_pairs");
	auto filteredSinglesOpts = SeqIOOptions::genFastqOutGz(outOpts.outFilename_.string() + "_filteredOffSingles");
	auto filteredPairedOpts = SeqIOOptions::genPairedOutGz(outOpts.outFilename_.string() + "_filteredOffPairs");
	auto unmappedSinglesOpts = SeqIOOptions::genFastqOutGz(outOpts.outFilename_.string() + "_unmappedSingles");
	auto unmappedPairedOpts = SeqIOOptions::genPairedOutGz(outOpts.outFilename_.string() + "_unmappedPairs");

	auto totalsCountsOpts = OutOptions(njh::files::make_path(outOpts.outFilename_.string() + "_totalReadCounts"), ".tab.txt");
	auto filteredChromCountsOpts = OutOptions(njh::files::make_path(outOpts.outFilename_.string() + "_filteredByChrom"), ".tab.txt");

	totalsCountsOpts.transferOverwriteOpts(outOpts);
	filteredChromCountsOpts.transferOverwriteOpts(outOpts);
	singlesOpts.out_.transferOverwriteOpts(outOpts);
	pairedOpts.out_.transferOverwriteOpts(outOpts);
	filteredSinglesOpts.out_.transferOverwriteOpts(outOpts);
	filteredPairedOpts.out_.transferOverwriteOpts(outOpts);
	unmappedSinglesOpts.out_.transferOverwriteOpts(outOpts);
	unmappedPairedOpts.out_.transferOverwriteOpts(outOpts);

	OutputStream totalsCountsOut(totalsCountsOpts);
	OutputStream filteredCountsOut(filteredChromCountsOpts);


	SeqOutput singlesWriter(singlesOpts);
	SeqOutput pairedWriter(pairedOpts);

	SeqOutput filteredSinglesWriter(filteredSinglesOpts);
	SeqOutput filteredPairedWriter(filteredPairedOpts);
	std::unique_ptr<SeqOutput> unmappedSinglesWriter;
	std::unique_ptr<SeqOutput> unmappedPairedWriter;

	if (writeOutUnmappedSeparately) {
		unmappedSinglesWriter = std::make_unique<SeqOutput>(unmappedSinglesOpts);
		unmappedPairedWriter = std::make_unique<SeqOutput>(unmappedPairedOpts);
	}

	BamAlnsCache alnCache;
	BamAlnsCache filterAlnCache;

	auto refData = bReader.GetReferenceData();
	std::vector<int32_t> filter_chroms_ref_ids;
	//check to make sure chroms contains chromosome from the input bam file
	VecStr missing;
	for(const auto & chrom : filter_chroms){
		bool found = false;
		for(const auto & ref : iter::enumerate(refData)){
			if(ref.element.RefName == chrom){
				found = true;
				filter_chroms_ref_ids.emplace_back(ref.index);
				break;
			}
		}
		if(!found){
			missing.emplace_back(chrom);
		}
	}

	auto onFilterChrom = [&filter_chroms_ref_ids](const int32_t refId){return njh::in(refId, filter_chroms_ref_ids);};

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
	uint64_t unmappedOrphans_ = 0;

	auto doesAlnPassSoftClipFilt = [&allowableSoftClipInAln,&refData](const BamTools::BamAlignment & bamAln){
		uint32_t softClipSum = 0;
		if(!bamAln.CigarData.empty() && 'S' == bamAln.CigarData.front().Type && 0 != bamAln.Position){
			softClipSum += bamAln.CigarData.front().Length;
		}
		if(bamAln.CigarData.size() > 1 && 'S' == bamAln.CigarData.back().Type && bamAln.GetEndPosition() != refData[bamAln.RefID].RefLength){
			softClipSum += bamAln.CigarData.back().Length;
		}
		return softClipSum <= allowableSoftClipInAln;
	};


	ReadCounts input;
	ReadCounts kept;
	ReadCounts unmapped;

	std::unordered_map<std::string, ReadCounts> filteredCountsByChrom;

	while (bReader.GetNextAlignment(bAln)) {
		//skip secondary alignments
		if (!bAln.IsPrimaryAlignment()) {
			continue;
		}

		const bool bAln_isPaired = bAln.IsPaired();
		const bool bAln_isMapped = bAln.IsMapped();
		const bool bAln_mateMapped = bAln.IsMateMapped();
		const bool bAln_isFirst = bAln.IsFirstMate();
		const bool bAln_isProper = bAln.IsProperPair();
		const bool bAln_onChrom      = bAln_isMapped && onFilterChrom(bAln.RefID);
		const bool bAln_mateOnChrom  = bAln_mateMapped && onFilterChrom(bAln.MateRefID);


		if (!bAln_isPaired) {
			++input.singles_;
			if (!bAln_isMapped) {
				if (writeOutUnmappedSeparately) {
					++unmapped.singles_;
					unmappedSinglesWriter->openWrite(bamAlnToSeqInfo(bAln));
				} else {
					++kept.singles_;
					singlesWriter.openWrite(bamAlnToSeqInfo(bAln));
				}
			} else {
				if (bAln_onChrom) {
					if(doesAlnPassSoftClipFilt(bAln) && bAln.MapQuality >= minMappingQuality){
						if(!doNotWriteFilterOff){
							filteredSinglesWriter.openWrite(bamAlnToSeqInfo(bAln));
						}
						++filteredCountsByChrom[refData[bAln.RefID].RefName].singles_;
					} else {
						if (writeOutUnmappedSeparately) {
							++unmapped.singles_;
							unmappedSinglesWriter->openWrite(bamAlnToSeqInfo(bAln));
						} else {
							++kept.singles_;
							singlesWriter.openWrite(bamAlnToSeqInfo(bAln));
						}
					}
				} else {
					++kept.singles_;
					singlesWriter.openWrite(bamAlnToSeqInfo(bAln));
				}
			}
		} else {
			++input.pairs_;
			if (bAln_isMapped &&
					bAln_mateMapped &&
					((any_mate && (bAln_onChrom ||
					bAln_mateOnChrom)) || (bAln_onChrom &&
					bAln_mateOnChrom)) &&
					(!requireProperPair || bAln_isProper)){
				if (!filterAlnCache.has(bAln.Name)) {
					//pair hasn't been added to cache yet so add to cache
					//this only works if mate and first mate have the same name
					filterAlnCache.add(bAln);
				} else {
					auto search = filterAlnCache.get(bAln.Name);
					//the below checks are in case any_mate is being used then only check for the soft clip filter and min mapp quality on the filter mapping alignments
					//if any_mate is not being use then both will be checked 
					bool baln_pass = !bAln_onChrom || (doesAlnPassSoftClipFilt(bAln) && bAln.MapQuality >= minMappingQuality);
					bool search_pass = !onFilterChrom(search->RefID) || (doesAlnPassSoftClipFilt(*search) && search->MapQuality >= minMappingQuality);
					if(baln_pass && search_pass){
						filteredCountsByChrom[njh::pasteAsStr(refData[search->RefID].RefName, "--", refData[bAln.RefID].RefName)].pairs_ += 2;
						if(!doNotWriteFilterOff){
							if (bAln_isFirst) {
								filteredPairedWriter.openWrite(PairedRead(bamAlnToSeqInfo(bAln), bamAlnToSeqInfo(*search),false));
							} else {
								filteredPairedWriter.openWrite(PairedRead(bamAlnToSeqInfo(*search), bamAlnToSeqInfo(bAln),false));
							}
						}
					} else {
						if (writeOutUnmappedSeparately) {
							unmapped.pairs_ += 2;
							if (bAln_isFirst) {
								unmappedPairedWriter->openWrite(PairedRead(bamAlnToSeqInfo(bAln), bamAlnToSeqInfo(*search),false));
							} else {
								unmappedPairedWriter->openWrite(PairedRead(bamAlnToSeqInfo(*search), bamAlnToSeqInfo(bAln),false));
							}
						} else {
							kept.pairs_ += 2;
							if (bAln_isFirst) {
								pairedWriter.openWrite(PairedRead(bamAlnToSeqInfo(bAln), bamAlnToSeqInfo(*search),false));
							} else {
								pairedWriter.openWrite(PairedRead(bamAlnToSeqInfo(*search), bamAlnToSeqInfo(bAln),false));
							}
						}
					}
					// now that operations have been computed, remove their other mate found from cache
					filterAlnCache.remove(search->Name);
				}
			}else{
				if (!alnCache.has(bAln.Name)) {
					//pair hasn't been added to cache yet so add to cache
					//this only works if mate and first mate have the same name
					alnCache.add(bAln);
				} else {
					auto search = alnCache.get(bAln.Name);
					if (filterWithUnmappedMate &&
					    (
						    (bAln_isMapped && !bAln_mateMapped && bAln_onChrom &&
						     doesAlnPassSoftClipFilt(bAln) && bAln.MapQuality >= minMappingQuality) ||
						    (!bAln_isMapped && bAln_mateMapped && bAln_mateOnChrom &&
						     doesAlnPassSoftClipFilt(*search) && search->MapQuality >= minMappingQuality)
					    )
					) {
						std::string filterChromName;
						if (bAln_isMapped && !bAln_mateMapped) {
							filterChromName = njh::pasteAsStr("unmapped", "--", refData[bAln.RefID].RefName);
						} else {
							filterChromName = njh::pasteAsStr(refData[search->RefID].RefName, "--", "unmapped");
						}
						filteredCountsByChrom[filterChromName].pairs_ += 2;
						if(!doNotWriteFilterOff){
							if (bAln_isFirst) {
								filteredPairedWriter.openWrite(PairedRead(bamAlnToSeqInfo(bAln), bamAlnToSeqInfo(*search),false));
							} else {
								filteredPairedWriter.openWrite(PairedRead(bamAlnToSeqInfo(*search), bamAlnToSeqInfo(bAln),false));
							}
						}
					} else {
						if (writeOutUnmappedSeparately &&
							(
								(!bAln_isMapped && !bAln_mateMapped) ||
							(
								( bAln_isMapped &&!bAln_mateMapped  && bAln_onChrom) ||
							  (!bAln_isMapped && bAln_mateMapped  && bAln_mateOnChrom)
							 )
							 )
							 ) {
							unmapped.pairs_ += 2;
							if (bAln_isFirst) {
								unmappedPairedWriter->openWrite(PairedRead(bamAlnToSeqInfo(bAln), bamAlnToSeqInfo(*search),false));
							} else {
								unmappedPairedWriter->openWrite(PairedRead(bamAlnToSeqInfo(*search), bamAlnToSeqInfo(bAln),false));
							}
						} else {
							kept.pairs_ += 2;
							if (bAln_isFirst) {
								pairedWriter.openWrite(PairedRead(bamAlnToSeqInfo(bAln), bamAlnToSeqInfo(*search),false));
							} else {
								pairedWriter.openWrite(PairedRead(bamAlnToSeqInfo(*search), bamAlnToSeqInfo(bAln),false));
							}
						}
					}
					// now that operations have been computed, remove ther other mate found from cache
					alnCache.remove(search->Name);
				}
			}
		}
	}

	//save the orphans;
	if (len(alnCache) > 0) {
		auto names = alnCache.getNames();
		for (const auto & name : names) {

			auto search = alnCache.get(name);
			if (writeOutUnmappedSeparately && !search->IsMapped()) {
				++unmappedOrphans_;
				unmappedSinglesWriter->openWrite(bamAlnToSeqInfo(*search));
			} else {
				++keptOrphans_;
				singlesWriter.openWrite(bamAlnToSeqInfo(*search));
			}
			alnCache.remove(name);
		}
	}
	if (len(filterAlnCache) > 0) {
		auto names = filterAlnCache.getNames();
		for (const auto & name : names) {
			auto search = filterAlnCache.get(name);
			if(doesAlnPassSoftClipFilt(*search) && search->MapQuality >= minMappingQuality) {
				++filteredOrphans_;
				if(!doNotWriteFilterOff){
					filteredSinglesWriter.openWrite(bamAlnToSeqInfo(*search));
				}
			} else {
				if (writeOutUnmappedSeparately) {
					++unmappedOrphans_;
					unmappedSinglesWriter->openWrite(bamAlnToSeqInfo(*search));
				} else {
					++keptOrphans_;
					singlesWriter.openWrite(bamAlnToSeqInfo(*search));
				}
			}
			filterAlnCache.remove(name);
		}
	}

	ReadCounts filtered;
	for(const auto & filt : filteredCountsByChrom){
		filtered.pairs_ += filt.second.pairs_;
		filtered.singles_ += filt.second.singles_;
	}

	auto bname = bfs::basename(setUp.pars_.ioOptions_.firstName_.filename());
	totalsCountsOut << "bam\tcondition\tcount\tfrac\ttotal" << std::endl;
	totalsCountsOut << bname
			<< "\t" << "keptPairs"
			<< "\t" << kept.pairs_
			<< "\t" << (input.pairs_ >0 ? kept.pairs_/static_cast<long double>(input.pairs_) : 0.0)
			<< "\t" << input.pairs_ << std::endl;
	totalsCountsOut << bname
			<< "\t" << "keptSingles"
			<< "\t" << kept.singles_
			<< "\t" << (input.singles_ > 0 ? kept.singles_/static_cast<long double>(input.singles_): 0.0)
			<< "\t" << input.singles_ << std::endl;

	totalsCountsOut << bname
			<< "\t" << "filteredPairs"
			<< "\t" << filtered.pairs_
			<< "\t" << (input.pairs_ >0 ? filtered.pairs_/static_cast<long double>(input.pairs_): 0.0)
			<< "\t" << input.pairs_ << std::endl;
	totalsCountsOut << bname
			<< "\t" << "filteredSingles"
			<< "\t" << filtered.singles_
			<< "\t" << (input.singles_ > 0 ? filtered.singles_/static_cast<long double>(input.singles_) : 0.0)
			<< "\t" << input.singles_ << std::endl;

	if (writeOutUnmappedSeparately) {
		totalsCountsOut << bname
				<< "\t" << "unmappedPairs"
				<< "\t" << unmapped.pairs_
				<< "\t" << (input.pairs_ >0 ? unmapped.pairs_ / static_cast<long double>(input.pairs_) : 0.0)
				<< "\t" << input.pairs_ << std::endl;
		totalsCountsOut << bname
				<< "\t" << "unmappedSingles"
				<< "\t" << unmapped.singles_
				<< "\t" << (input.singles_ > 0 ? unmapped.singles_ / static_cast<long double>(input.singles_) : 0.0)
				<< "\t" << input.singles_ << std::endl;
	}

	totalsCountsOut << bname
			<< "\t" << "keptOrphans"
			<< "\t" << keptOrphans_
			<< "\t" << (input.pairs_ >0 ? keptOrphans_ / static_cast<long double>(input.pairs_) : 0.0)
			<< "\t" << input.pairs_ << std::endl;

	totalsCountsOut << bname
			<< "\t" << "filteredOrphans"
			<< "\t" << filteredOrphans_
			<< "\t" << (input.pairs_ >0 ? filteredOrphans_ / static_cast<long double>(input.pairs_) : 0.0)
			<< "\t" << input.pairs_ << std::endl;
	if (writeOutUnmappedSeparately) {
		totalsCountsOut << bname
				<< "\t" << "unmappedOrphans"
				<< "\t" << unmappedOrphans_
				<< "\t" << (input.pairs_ >0 ? unmappedOrphans_ / static_cast<long double>(input.pairs_) : 0.0)
				<< "\t" << input.pairs_ << std::endl;
	}

	auto names = getVectorOfMapKeys(filteredCountsByChrom);
	njh::sort(names);
	filteredCountsOut << "bam\tchrom\tpairs\tpairsFrac\tsingles\tsinglesFrac" << std::endl;
	for(const auto & name : names){
		filteredCountsOut << bname
				<< "\t" << name
				<< "\t" << filteredCountsByChrom[name].pairs_
				<< "\t" << filteredCountsByChrom[name].pairs_/static_cast<long double>(filtered.pairs_)
				<< "\t" << filteredCountsByChrom[name].singles_
				<< "\t" << filteredCountsByChrom[name].singles_/static_cast<long double>(filtered.singles_) << std::endl;
	}

	return 0;
}


} // namespace njhseq

