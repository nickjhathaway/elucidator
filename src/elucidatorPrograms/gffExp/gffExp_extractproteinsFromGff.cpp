//
// Created by Nicholas Hathaway on 9/14/23.
//
#include "gffExp.hpp"
#include "elucidator/objects/BioDataObject.h"
#include "elucidator/seqToolsUtils/seqToolsUtils.hpp"

#include <njhseq/objects/Gene.h>
#include <TwoBit.h>

namespace njhseq {

int gffExpRunner::extractProteinsFromGff(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path gffFnp = "";
	bfs::path twoBitFnp = "";
	uint32_t numThreads = 1;
	bool forceUntranslatableProteins = false;
	bool addDescriptionToName = false;
	VecStr derivedFromRecordFeatures {"polypeptide"};
	VecStr acceptableMrnaLikeRecords = {"mRNA", "CDS", "transcript"};
	//VecStr acceptableMrnaLikeRecords = {"mRNA"};
	VecStr allowableFeatureType {"gene", "protein_coding_gene"};
	std::set<char> allowableStarts = {'M'};
	std::set<std::string> selectedGeneIDs;
	auto seqOutOpts = SeqIOOptions::genFastaOutGz("out");

	OutOptions infoKeyOutOpts(bfs::path(""), ".tsv.gz");
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.setOption(gffFnp, "--gff", "GFF (gene feature format) file", true);
	setUp.setOption(twoBitFnp, "--2bit", "2bit file of genome", true);
	setUp.setOption(selectedGeneIDs, "--selectedGeneIDs", "only extract these selected gene IDs");
	setUp.setOption(allowableFeatureType, "--allowableFeatureType", "allowable Feature Type for extracting");
	setUp.setOption(acceptableMrnaLikeRecords, "--acceptableMrnaLikeRecords", "allowable mRNA like Type for extracting");
	setUp.setOption(forceUntranslatableProteins, "--forceUntranslatableProteins", "force Untranslatable Proteins, even when feature is gene and not pseudogene the record is still a pseudogene and has premature stop codons");
	setUp.setOption(infoKeyOutOpts.outFilename_, "--infoOutname", "name of output file to contain description info on proteins");
	setUp.setOption(numThreads, "--numThreads", "number of threads to use");
	setUp.setOption(allowableStarts, "--allowableStarts", "allowable Start codons");
	setUp.setOption(addDescriptionToName, "--addDescriptionToName", "add Description To Name");

	setUp.processWritingOptions(seqOutOpts.out_);
	infoKeyOutOpts.transferOverwriteOpts(seqOutOpts.out_);
	if(infoKeyOutOpts.outFilename_.string().empty()){
		infoKeyOutOpts.outFilename_ = njh::files::prependFileBasename(seqOutOpts.out_.outFilename_, "info");
	}
	setUp.finishSetUp(std::cout);

	TwoBit::TwoBitFile globaltReader(twoBitFnp);

	OutputStream infoOut(infoKeyOutOpts);
	infoOut << "GeneID\ttranscriptID\tavailable_description" << std::endl;
	SeqOutput out(seqOutOpts);
	out.openOut();

	//scan for genes

	std::unordered_set<std::string> geneIds;
	std::unordered_set<std::string> extratedGeneIds;

	//first gather all GeneIDs
	{
		BioDataFileIO<GFFCore> reader{IoOptions(InOptions(gffFnp))};
		std::vector<std::shared_ptr<GFFCore>> cache;
		reader.openIn();
		std::string line = "";
		std::shared_ptr<GFFCore> gRecord = reader.readNextRecord();
		while (nullptr != gRecord) {
			if(njh::in(gRecord->type_, allowableFeatureType)) {
				//possibly geneID, add, if doing only specific gene IDs then add only if it's in the selected IDs
				if(selectedGeneIDs.empty() || njh::in(gRecord->getAttr("ID"), selectedGeneIDs)){
					geneIds.emplace(gRecord->getAttr("ID"));
				}
			}
			bool end = false;
			while ('#' == reader.inFile_->peek()) {
				if (njh::files::nextLineBeginsWith(*reader.inFile_, "##FASTA")) {
					end = true;
					break;
				}
				njh::files::crossPlatGetline(*reader.inFile_, line);
			}
			if (end) {
				break;
			}
			gRecord = reader.readNextRecord();
		}
	}
	if(!selectedGeneIDs.empty()){
		VecStr missingIDs;
		for(const auto & id : selectedGeneIDs){
			if(!njh::in(id, geneIds)){
				missingIDs.emplace_back(id);
			}
		}
		if(!missingIDs.empty()){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "missing the following ids: " << njh::conToStr(missingIDs, ",") << "\n";
			throw std::runtime_error{ss.str()};
		}
	}

	std::unordered_map<std::string, std::set<std::string>> parents;
	std::unordered_map<std::string, std::vector<std::shared_ptr<GFFCore>>> gffRecs;
	std::unordered_map<std::string, std::shared_ptr<GeneFromGffs>> genes;


	BioDataFileIO<GFFCore> reader{IoOptions(InOptions(gffFnp))};
	std::vector<std::shared_ptr<GFFCore>> cache;
	reader.openIn();
	// uint32_t count = 0;
	std::string line;
	std::shared_ptr<GFFCore> gRecord = reader.readNextRecord();
	while (nullptr != gRecord) {
		if(gRecord->hasAttr("ID") && njh::in(gRecord->getAttr("ID"), geneIds) ){
			if(!njh::in(gRecord->type_, allowableFeatureType)) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error feature type needs to be gene, not " << gRecord->type_ << " for: " << gRecord->getAttr("ID") << "\n";
				throw std::runtime_error{ss.str()};
			}
			auto currentId = gRecord->getAttr("ID");
			parents[currentId].insert(gRecord->getAttr("ID"));
			gffRecs[currentId].emplace_back(std::make_shared<GFFCore>(*gRecord));
		} else if(gRecord->hasAttr("Parent")){
			auto currentParents = tokenizeString(gRecord->getAttr("Parent"), ",");
			for(const auto & currentParent : currentParents){
				for(auto & p : parents){
					if(njh::in(currentParent, p.second)){
						p.second.insert(gRecord->getAttr("ID"));
						gffRecs[p.first].emplace_back(std::make_shared<GFFCore>(*gRecord));
					}
				}
			}
		} else if(njh::in(gRecord->type_, derivedFromRecordFeatures) ){
			cache.emplace_back(std::make_shared<GFFCore>(*gRecord));
		}
		if(njh::in(gRecord->type_, allowableFeatureType)){
			for(const auto & fromCache : cache){
				//grab any possible polypeptides or misc records,
				//I have found that sometimes these come before the gene record so you have to keep a cache around since the last cache
				//there might be a better way of doing this
				if(fromCache->hasAttr("Derives_from")){
					auto currentParent = fromCache->getAttr("Derives_from");
					for(auto & p : parents){
						if(njh::in(currentParent, p.second)){
							gffRecs[p.first].emplace_back(std::make_shared<GFFCore>(*fromCache));
						}
					}
				}
			}
			cache.clear();
		}
		bool end = false;
		while ('#' == reader.inFile_->peek()) {
			if (njh::files::nextLineBeginsWith(*reader.inFile_, "##FASTA")) {
				end = true;
				break;
			}
			njh::files::crossPlatGetline(*reader.inFile_, line);
		}
		if (end) {
			break;
		}
		gRecord = reader.readNextRecord();
		// ++count;
	}

	auto allGeneIds = njh::getVecOfMapKeys(gffRecs);
	njh::sort(allGeneIds);
	njh::concurrent::LockableVec<std::string> geneRecordName(allGeneIds);
	std::mutex outMut;

	std::function<void()> extractProteinSeq = [&geneRecordName,&outMut,
														&gffRecs, &out, &extratedGeneIds, &forceUntranslatableProteins, &twoBitFnp, &acceptableMrnaLikeRecords, &infoOut, &allowableStarts, &addDescriptionToName](){
		TwoBit::TwoBitFile tReader(twoBitFnp);

		std::string geneID;
		while(geneRecordName.getVal(geneID)){
			//check to make sure these "gene" typed features actually contain mRNA as some annotation tools will
			//annotate rRNA, ncRNA as "gene" cause why not
			bool contains_mRNA = false;

			const auto & geneGffs = gffRecs.at(geneID);

			for (const auto &rec: geneGffs) {
				if (njh::in(rec->type_, acceptableMrnaLikeRecords)) {
					contains_mRNA = true;
				}
			}

			bool contains_skip_region = false;
			VecStr filterSubRegionFeatures_{"rrna", "trna", "snorna","snrna","ncrna"};
			for (const auto &rec: geneGffs) {
				if (njh::in(njh::strToLowerRet(rec->type_), filterSubRegionFeatures_)) {
					contains_skip_region = true;
				}
			}
			if (contains_mRNA && !contains_skip_region) {
				auto geneInfo = std::make_shared<GeneFromGffs>(geneGffs);
				auto gsInfos = geneInfo->generateGeneSeqInfo(tReader, false);
				for (const auto &genGeneInfo: gsInfos) {
					auto modProtein = genGeneInfo.second->protein_;
					bool endsWithStop = modProtein.seq_.back() == '*';
					bool startsWithStartCodon = njh::in(modProtein.seq_.front(), allowableStarts);
					modProtein.rstrip('*');
					//check for proteins with premature stop codons
					if (forceUntranslatableProteins || (startsWithStartCodon && endsWithStop && std::string::npos == modProtein.seq_.find('*'))) {
						modProtein.name_ = genGeneInfo.first; //rename to transcript name
						if(addDescriptionToName){
							modProtein.name_.append(" ");
							modProtein.name_.append(geneInfo->getOneGeneDetailedName());
						}
						std::lock_guard<std::mutex> lock(outMut);
						out.write(modProtein);
						extratedGeneIds.emplace(geneID);

						infoOut << geneID << "\t" << genGeneInfo.first << "\t" << geneInfo->getOneGeneDetailedName() << std::endl;
					}
				}
			}
		}
	};

	njh::concurrent::runVoidFunctionThreaded(extractProteinSeq, numThreads);

	if(!selectedGeneIDs.empty()){
		VecStr missingIDs;
		for(const auto & id : selectedGeneIDs){
			if(!njh::in(id, extratedGeneIds)){
				missingIDs.emplace_back(id);
			}
		}
		if(!missingIDs.empty()){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "missing the following ids: " << njh::conToStr(missingIDs, ",") << "\n";
			throw std::runtime_error{ss.str()};
		}
	}
	return 0;
}



} //namespace njhseq

