/*
 * bamExp_BamGetReadInfoForRegions.cpp
 *
 *  Created on: Jun 9, 2020
 *      Author: nick
 */
#include "bamExp.hpp"
#include <TwoBit.h>

#include "elucidator/BamToolsUtils.h"
#include <njhseq/objects/BioDataObject.h>




namespace njhseq {
int bamExpRunner::BamGetPairedReadInfoForRegions(const njh::progutils::CmdArgs & inputCommands){
	bfs::path bedFnp = "";
	OutOptions outOpts(bfs::path(""), ".tab.txt");


	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(bedFnp, "--bedFnp", "bedFnp", true);
	setUp.processReadInNames({"--bam"}, true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);


	BamTools::BamAlignment bAln;
	BamTools::BamReader bReader;
	bReader.Open(setUp.pars_.ioOptions_.firstName_.string());
	checkBamOpenThrow(bReader, setUp.pars_.ioOptions_.firstName_.string());
	loadBamIndexThrow(bReader);
	auto refs = bReader.GetReferenceData();
	BamAlnsCacheWithRegion alnCache;


	auto regions = bedPtrsToGenomicRegs(getBeds(bedFnp));
	OutputStream out(outOpts);
	out << "name\tfirstMate\tchrom\tstart\tend\tstrand\tinsertSize\tintersectedRegion" << std::endl;
	std::vector<GenomicRegion> regionsSearched;
	for(const auto & reg : regions){

		setBamFileRegionThrow(bReader, reg);
		while (bReader.GetNextAlignment(bAln)) {
			if (!bAln.IsPrimaryAlignment()) {
				continue;
			}
			bool overlapsWithPreviousRegion = false;
			for(const auto & previous : regionsSearched){
				if(previous.overlaps(GenomicRegion(bAln, refs))){
					overlapsWithPreviousRegion = true;
					break;
				}
			}
			if(overlapsWithPreviousRegion){
				continue;
				//skip, should have already have been grabbed in the previous region search
				//this will prevent getting a read twice if multiple regions given are close to each other.
			}
			if (!bAln.IsPaired() || !bAln.IsMapped() || !bAln.IsMateMapped()) {
				//skip single, or if one of the mates is unmmaped
				continue;
			} else {
				if (!alnCache.has(bAln.Name)) {
					//pair hasn't been added to cache yet so add to cache
					//this only works if mate and first mate have the same name
					alnCache.addWithRegion(bAln, reg);
					continue;
				} else {
					auto search = alnCache.get(bAln.Name);
					out << bAln.Name
							<< "\t" << njh::boolToStr(bAln.IsFirstMate())
							<< "\t" << refs[bAln.RefID].RefName
							<< "\t" << bAln.Position
							<< "\t" << bAln.GetEndPosition()
							<< "\t" << (bAln.IsReverseStrand() ? '-': '+')
							<< "\t" << bAln.InsertSize
							<< "\t" << reg.uid_ << std::endl;
					out << search->Name
							<< "\t" << njh::boolToStr(search->IsFirstMate())
							<< "\t" << refs[search->RefID].RefName
							<< "\t" << search->Position
							<< "\t" << search->GetEndPosition()
							<< "\t" << (search->IsReverseStrand() ? '-': '+')
							<< "\t" << search->InsertSize
							<< "\t" << reg.uid_ << std::endl;
					// now that operations have been computed, remove ther other mate found from cache
					alnCache.remove(search->Name);
					continue;
				}
			}
		}
		regionsSearched.emplace_back(reg);
	}


	//find the orphan mates which have mates that fall outside of the regions given
	if (len(alnCache) > 0) {
		auto names = alnCache.getNames();

		//find orphans' mates if possible;
		BamTools::BamReader bReaderMateFinder;
		bReaderMateFinder.Open(setUp.pars_.ioOptions_.firstName_.string());
		checkBamOpenThrow(bReaderMateFinder, setUp.pars_.ioOptions_.firstName_.string());
		loadBamIndexThrow(bReaderMateFinder);
		auto refData = bReaderMateFinder.GetReferenceData();
		//gather all the orphans regions
		std::unordered_map<std::string, std::set<uint32_t>> orphanPositions;
		for (const auto & name : names) {
			auto search = alnCache.get(name);
			if(search->IsPaired()){
				if (search->IsMateMapped()) {
					orphanPositions[refData[search->MateRefID].RefName].emplace(search->MatePosition);
				}
			}
		}
		std::vector<GenomicRegion> orphanMateRegions;
		for(const auto & orPos : orphanPositions){
			for(const auto & pos  : orPos.second){
				orphanMateRegions.emplace_back(GenomicRegion("", orPos.first, pos, pos + 1, false));
			}
		}
		for(const auto & reg : orphanMateRegions){
			setBamFileRegionThrow(bReaderMateFinder, reg);
			while (bReaderMateFinder.GetNextAlignment(bAln)) {
				//skip secondary alignments
				if (!bAln.IsPrimaryAlignment()) {
					continue;
				}
				if(static_cast<uint32_t>(bAln.Position) < reg.start_){
					continue;
				}
				if (bAln.IsPaired()) {
					if (alnCache.has(bAln.Name)) {
						auto search = alnCache.get(bAln.Name);
						auto searchRegion = alnCache.getRegion(bAln.Name);
						out << bAln.Name
								<< "\t" << njh::boolToStr(bAln.IsFirstMate())
								<< "\t" << refs[bAln.RefID].RefName
								<< "\t" << bAln.Position
								<< "\t" << bAln.GetEndPosition()
								<< "\t" << (bAln.IsReverseStrand() ? '-': '+')
								<< "\t" << bAln.InsertSize
								<< "\t" << std::endl;
						out << search->Name
								<< "\t" << njh::boolToStr(search->IsFirstMate())
								<< "\t" << refs[search->RefID].RefName
								<< "\t" << search->Position
								<< "\t" << search->GetEndPosition()
								<< "\t" << (search->IsReverseStrand() ? '-': '+')
								<< "\t" << search->InsertSize
								<< "\t" << searchRegion->uid_ << std::endl;
						// now that operations have been computed, remove ther other mate found from cache
						alnCache.remove(search->Name);
					}
				}
			}
		}
		if(setUp.pars_.debug_){
			auto names = alnCache.getNames();
			std::cout << "names.size(): " << names.size() << std::endl;
		}
	}


	return 0;

}


}  // namespace njhseq
