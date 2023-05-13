/*
 * parsingFileExpRunner.cpp
 *
 *  Created on: May 18, 2015
 *      Author: nickhathaway
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


#include "parsingFileExp.hpp"
#include <njhseq/objects/BioDataObject/BioDataFileIO.hpp>
#include <njhseq/objects/BioDataObject/BLASTHitTabular.hpp>
#include <SeekDeep/objects/TarAmpSetupUtils/PrimersAndMids.hpp>

namespace njhseq {
parsingFileExpRunner::parsingFileExpRunner()
    : njh::progutils::ProgramRunner(
          {
					 addFunc("parsingBioSampleSetNCBIJson", parsingBioSampleSetNCBIJson, false),
					 addFunc("getAttributeLevelsBioSampleSetNCBIJson", getAttributeLevelsBioSampleSetNCBIJson, false),
					 addFunc("getStrainInfoTableBioSampleSetNCBIJson", getStrainInfoTableBioSampleSetNCBIJson, false),
					 addFunc("parseBlastpHitsTab", parseBlastpHitsTab, false),
					 addFunc("parseSTOCKHOLM", parseSTOCKHOLM, false),
					 addFunc("parseSTOCKHOLMToFasta", parseSTOCKHOLMToFasta, false),
					 addFunc("parsehmmerDomainHitTab", parsehmmerDomainHitTab, false),
					 addFunc("quickCountFastq", quickCountFastq, false),
					 addFunc("quickCountFasta", quickCountFasta, false),
					 addFunc("quickCountDirectory", quickCountDirectory, false),
					 addFunc("BlastpHitsTabToBed", BlastpHitsTabToBed, false),
					 addFunc("parsePrimerFastaToPrimerTxt", parsePrimerFastaToPrimerTxt, false),
					 addFunc("parseNucmerResultsToBed", parseNucmerResultsToBed, false),
					 addFunc("parseMashTriangleResults", parseMashTriangleResults, false),
					 addFunc("parseMummberResultsToBed", parseMummberResultsToBed, false),
           },
          "parsingFileExp") {}
//,

class BioSampleSetNCBIJsonParse {
public:
	BioSampleSetNCBIJsonParse(const bfs::path & fnp) :
			fnp_(fnp), info_(njh::json::parseFile(fnp_.string())) {
	}

	const bfs::path fnp_;

	const Json::Value info_;

	std::set<std::string> attributeLevels_;

	void setAttributeLevelsAll(){
		attributeLevels_.clear();
		const auto & bioSampleList = info_["BioSampleSet"]["BioSample"];
		for(const auto & sample : bioSampleList){
			auto levels = getAttrLevels(sample);
			attributeLevels_.insert(levels.begin(), levels.end());
		}
	}

	static VecStr getAttrLevels(const Json::Value & sample){
		VecStr ret;
		if(sample.isMember("Attributes") && sample["Attributes"].isMember("Attribute") && sample["Attributes"]["Attribute"].isArray()){
			for(const auto & attr : sample["Attributes"]["Attribute"]){
				if(attr.isMember("attribute_name")){
					ret.emplace_back(attr["attribute_name"].asString());
				}
			}
		}
		return ret;
	}

	static std::string getLevelForAttrName(const Json::Value & sample, const std::string & attrName){
		std::string ret = "NA";
		if(sample.isMember("Attributes") && sample["Attributes"].isMember("Attribute") && sample["Attributes"]["Attribute"].isArray()){
			for(const auto & attr : sample["Attributes"]["Attribute"]){
				if(attr.isMember("attribute_name")){
					if(attrName == attr["attribute_name"].asString()){
						ret = attr["attribute_name"].asString();
					}
				}
			}
		}
		return ret;
	};

	static std::string getMemeberValue(const Json::Value & sample,
			const std::string & member) {
		std::string ret = "NA";
		if(sample.isMember(member)){
			ret = sample[member].asString();
		}
		return ret;
	}

	static std::string getMemeberValueNested1(const Json::Value & sample,
			const std::string & member, const std::string & nestedMember) {
		std::string ret = "NA";
		if(sample.isMember(member) && sample[member].isMember(nestedMember)){
			ret = sample[member][nestedMember].asString();
		}
		return ret;
	}

	static std::string getAttributeValues(const Json::Value & sample,
			std::function<std::string(const Json::Value &)> extractorFunc) {

		std::string ret = "NA";

		if (sample.isMember("Attributes")
				&& sample["Attributes"].isMember("Attribute")
				&& sample["Attributes"]["Attribute"].isArray()) {
			ret = extractorFunc(sample["Attributes"]["Attribute"]);
		}
		return ret;
	}

};

int parsingFileExpRunner::getStrainInfoTableBioSampleSetNCBIJson(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path jsonFile = "";
	auto outOpts = TableIOOpts::genTabFileOut("out.tab.txt", true);
	seqSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processVerbose();
	setUp.processWritingOptions(outOpts.out_);
	setUp.setOption(jsonFile, "--jsonFile",
				"A JSON file convert by xml2json from xml downloaded from NCBI search", true);
	setUp.finishSetUp(std::cout);

	BioSampleSetNCBIJsonParse parser(jsonFile);
	const auto & bioSampleList = parser.info_["BioSampleSet"]["BioSample"];
	table outTable(VecStr{"id", "accession", "Title", "strain"});

	for(const auto & sample : bioSampleList){
		outTable.addRow(
									 BioSampleSetNCBIJsonParse::getMemeberValue(sample, "id"),
									 BioSampleSetNCBIJsonParse::getMemeberValue(sample, "accession"),
									 BioSampleSetNCBIJsonParse::getMemeberValueNested1(sample, "Description", "Title"),
									 BioSampleSetNCBIJsonParse::getAttributeValues(sample, [](const Json::Value & val){
			std::string ret = "NA";
			for(const auto & attr : val){
				if(attr.isMember("attribute_name") && "strain" == njh::strToLowerRet(attr["attribute_name"].asString())){
					if(attr.isMember("$t")){
						if("NA" == ret){
							ret = attr["$t"].asString();
						}else{
							ret.append(",");
							ret.append(attr["$t"].asString());
						}
					}
				}
			}
			return ret;
		}));
	}
	outTable.outPutContents(outOpts);
	return 0;
}

int parsingFileExpRunner::parsingBioSampleSetNCBIJson(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path jsonFile = "";
	auto outOpts = TableIOOpts::genTabFileOut("out.tab.txt", true);
	seqSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processVerbose();
	setUp.processWritingOptions(outOpts.out_);
	setUp.setOption(jsonFile, "--jsonFile",
				"A JSON file convert by xml2json from xml downloaded from NCBI search", true);
	setUp.finishSetUp(std::cout);

	BioSampleSetNCBIJsonParse parser(jsonFile);
	const auto & bioSampleList = parser.info_["BioSampleSet"]["BioSample"];
	uint32_t noAccession = 0;
	uint32_t noId = 0;
	uint32_t noDescription = 0;
	for(const auto & sample : bioSampleList){
		if(sample.isMember("accession")){
			std::cout << "accession" << std::endl;
			std::cout << sample["accession"] << std::endl;
		}else{
			++	noAccession;
		}
		if(sample.isMember("id")){
			std::cout << "id" << std::endl;
			std::cout << sample["id"] << std::endl;
		}else{
			++	noId;
		}
		if(sample.isMember("Description")){
			std::cout << "Description" << std::endl;
			std::cout << sample["Description"] << std::endl;
		}else{
			++	noDescription;
		}
	}
	std::cout << "noId: " << noId << std::endl;
	std::cout << "noAccession: " << noAccession << std::endl;
	std::cout << "noDescription: " << noDescription << std::endl;
	return 0;
}



int parsingFileExpRunner::getAttributeLevelsBioSampleSetNCBIJson(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path jsonFile = "";
	auto outOpts = TableIOOpts::genTabFileOut("out.tab.txt", true);
	seqSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processVerbose();
	setUp.processWritingOptions(outOpts.out_);
	setUp.setOption(jsonFile, "--jsonFile",
				"A JSON file convert by xml2json from xml downloaded from NCBI search", true);
	setUp.finishSetUp(std::cout);

	BioSampleSetNCBIJsonParse parser(jsonFile);
	parser.setAttributeLevelsAll();
	std::cout << njh::conToStr(parser.attributeLevels_, "\n") << std::endl;
	return 0;
}

int parsingFileExpRunner::parseBlastpHitsTab(const njh::progutils::CmdArgs & inputCommands) {
	IoOptions ioOpts;
	ioOpts.out_.outExtention_ = ".json";
	seqSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processVerbose();
	setUp.processWritingOptions(ioOpts.out_);
	setUp.setOption(ioOpts.in_.inFilename_, "--hitFnp", "output from blastp filename", true);
	setUp.finishSetUp(std::cout);

	BioDataFileIO<BLASTHitTab> reader(ioOpts);
	reader.openIn();
	reader.openOut();
	BLASTHitTab hit;
	while(reader.readNextRecord(hit)){
		(*reader.out_) << hit.toJson() << std::endl;
	}

	return 0;
}

int parsingFileExpRunner::BlastpHitsTabToBed(const njh::progutils::CmdArgs & inputCommands) {
	bool renameWithHitCount = false;

	IoOptions ioOpts;
	ioOpts.out_.outExtention_ = ".bed";
	seqSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processVerbose();
	setUp.processWritingOptions(ioOpts.out_);
	setUp.setOption(ioOpts.in_.inFilename_, "--hitFnp", "output from blastp filename", true);
	setUp.setOption(renameWithHitCount, "--renameWithHitCount", "rename With Hit Count");

	setUp.finishSetUp(std::cout);

	BioDataFileIO<BLASTHitTab> reader(ioOpts);
	reader.openIn();
	reader.openOut();
	BLASTHitTab hit;
	std::unordered_map<std::string, uint32_t> hitCounts;

	while(reader.readNextRecord(hit)){
		auto bedLoc = hit.genSubjectBed6();
		if(renameWithHitCount){
			++hitCounts[bedLoc.name_];
			if(hitCounts[bedLoc.name_] > 1){
				bedLoc.name_.append(njh::pasteAsStr(".",hitCounts[bedLoc.name_] ));
			}
		}
		(*reader.out_) << bedLoc.toDelimStrWithExtra() << std::endl;
	}

	return 0;
}


class QuickSeqFileCounter {
public:
	QuickSeqFileCounter() = default;


	static uint32_t quickCountFasta(const InOptions & inputOpts){
		if(inputOpts.inExists()){
			InputStream  in(inputOpts);
			return quickCountFasta(in);
		}
		return 0;
	}

	static uint32_t quickCountFasta(InputStream & in){
		uint32_t count = 0;
		std::string line;
		while(njh::files::crossPlatGetline(in, line)){
			if('>' == line.front()){
				++count;
			}
		}
		return count;
	}

	static uint32_t quickCountFastq(const InOptions & inputOpts){
		if(inputOpts.inExists()){
			InputStream  in(inputOpts);
			return quickCountFastq(in);
		}
		return 0;
	}

	static uint32_t quickCountFastq(InputStream & in){
		std::string line;
		uint32_t lineCount = 0;
		while(njh::files::crossPlatGetline(in, line)){
			++lineCount;
		}
		return lineCount/4;
	}

};


int parsingFileExpRunner::quickCountFastq(const njh::progutils::CmdArgs & inputCommands) {
	IoOptions ioOpts;
	seqSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processVerbose();
	setUp.processWritingOptions(ioOpts.out_);
	setUp.setOption(ioOpts.in_.inFilename_, "--fastq", "fastq file", true);
	setUp.finishSetUp(std::cout);
	OutputStream out(ioOpts.out_);
	out << QuickSeqFileCounter::quickCountFastq(ioOpts.in_) << std::endl;
	return 0;
}


int parsingFileExpRunner::quickCountFasta(const njh::progutils::CmdArgs & inputCommands) {
	IoOptions ioOpts;
	seqSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processVerbose();
	setUp.processWritingOptions(ioOpts.out_);
	setUp.setOption(ioOpts.in_.inFilename_, "--fasta", "fasta file", true);
	setUp.finishSetUp(std::cout);
	OutputStream out(ioOpts.out_);
	out << QuickSeqFileCounter::quickCountFasta(ioOpts.in_) << std::endl;
	return 0;
}

int parsingFileExpRunner::quickCountDirectory(const njh::progutils::CmdArgs & inputCommands) {
	OutOptions outOpts;
	bfs::path dirName;
	uint32_t numThreads = 1;
	seqSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processVerbose();
	setUp.processWritingOptions(outOpts);
	setUp.setOption(numThreads, "--numThreads", "number of threads to use");
	setUp.setOption(dirName, "--dirName", "Directory to count all fastas and/or fastqs", true);
	setUp.finishSetUp(std::cout);

	std::unordered_map<std::string, std::set<std::string>> fileExtensionLookUps;
	fileExtensionLookUps["fastq"] = std::set<std::string>{".fq", ".fq.gz", ".fastq", ".fastq.gz", ".fnq", ".fnq.gz"};
	fileExtensionLookUps["fasta"] = std::set<std::string>{".fa", ".fa.gz", ".fasta", ".fasta.gz", ".fna", ".fna.gz"};

	std::unordered_map<std::string, std::function<uint32_t(const bfs::path & fnp)>> fileCountingFuncs;
	fileCountingFuncs.emplace("fasta", [](const bfs::path & fnp){
		return QuickSeqFileCounter::quickCountFasta(fnp);
	});
	fileCountingFuncs.emplace("fastq", [](const bfs::path & fnp){
		return QuickSeqFileCounter::quickCountFastq(fnp);
	});

	OutputStream out(outOpts);
	out << "file\tguessedFormat\tcount" << std::endl;

	auto allFiles = njh::files::filesInFolder(dirName);
	njh::sort(allFiles);
	njh::concurrent::LockableQueue<bfs::path> pathsQueue(allFiles);
	std::mutex outMut;
	std::function<void()> countFiles = [&fileExtensionLookUps,&pathsQueue,&outMut,&out,&setUp,&fileCountingFuncs](){
		bfs::path f;
		while(pathsQueue.getVal(f)){
			std::string guessFormat = "unknown";
			for(const auto & fileExtens : fileExtensionLookUps){
				for(const auto & ext : fileExtens.second){
					if(njh::endsWith(f.string(), ext))	{
						guessFormat = fileExtens.first;
						break;
					}
				}
			}

			if("unknown" != guessFormat){
				auto count = fileCountingFuncs.at(guessFormat)(f);
				{
					std::lock_guard<std::mutex> lock(outMut);
					out << f.filename().string() << "\t" << guessFormat << "\t" << count << std::endl;
				}
			} else if(setUp.pars_.verbose_){
				std::cerr << "Unknown file type for " << f << std::endl;
			}
		}
	};
	njh::concurrent::runVoidFunctionThreaded(countFiles, numThreads);
	return 0;
}

int parsingFileExpRunner::parsePrimerFastaToPrimerTxt(const njh::progutils::CmdArgs &inputCommands){
	IoOptions ioOpts;
	std::string forwardRegexStr = "(.+)_[Ff]$";
	std::string reverseRegexStr = "(.+)_[Rr]$";

	seqSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processVerbose();
	setUp.processWritingOptions(ioOpts.out_);
	setUp.setOption(ioOpts.in_.inFilename_, "--fasta", "fasta file", true);
	setUp.setOption(forwardRegexStr, "--forwardRegex", "Regex pattern for forward primer");
	setUp.setOption(reverseRegexStr, "--reverseRegex", "Regex pattern for reverse primer");

	setUp.finishSetUp(std::cout);
	std::regex forwardRegex{forwardRegexStr};
	std::regex reverseRegex{reverseRegexStr};

	seqInfo seq;
	SeqInput reader(SeqIOOptions::genFastaIn(ioOpts.in_.inFilename_));
	reader.openIn();
	OutputStream out(ioOpts.out_);

	struct SimplePair{
		SimplePair() = default;
		std::string for_;
		std::string rev_;
	};
	std::unordered_map<std::string, SimplePair> primerPairs;
	while(reader.readNextRead(seq)){
		std::smatch forwardMatch;
		std::smatch reverseMatch;
		if(std::regex_match(seq.name_,forwardMatch,  forwardRegex)){
			primerPairs[forwardMatch[1]].for_ = seq.seq_;
		}else if(std::regex_match(seq.name_,reverseMatch,  reverseRegex)){
			primerPairs[reverseMatch[1]].rev_ = seq.seq_;
		}else{
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "couldn't match " << seq.name_ << "to either forward pattern: " << forwardRegexStr << " or reversePattern: " << reverseRegexStr << "\n";
			throw std::runtime_error{ss.str()};
		}
	}

	std::map<std::string, PrimersAndMids::Target>	targets;

	//check to see if a forward and reverse was found for all primer pairs;
	for(const auto & primerPair : primerPairs){
		if(primerPair.second.for_.empty()){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "no forward primer for " << primerPair.first << "\n";
			throw std::runtime_error{ss.str()};
		}
		if(primerPair.second.rev_.empty()){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "no reverse primer for " << primerPair.first << "\n";
			throw std::runtime_error{ss.str()};
		}
		targets.emplace(primerPair.first, PrimersAndMids::Target(primerPair.first, primerPair.second.for_, primerPair.second.rev_));
	}

	out << "target\tforward\treverse" << std::endl;
	for(const auto & target : targets){
		out << target.first
		<< "\t" << target.second.info_.forwardPrimerRaw_
		<< "\t" << target.second.info_.reversePrimerRaw_
		<< std::endl;
	}
	return 0;
}


class MummerOutput {
public:
	//excepts mummer to be run like with the -F option, so expects 4 columns each time regardless of the number of seqs in reference or query
	//e.g. mummer -b -c -F ref.fasta query.fasta

	// mummer -b -c -F -l 200 -maxmatch, maxmatch does all matches regarless of how unique they are in each sequence
	class MummerTargetHits {
	public:
		class MummerRecord {
		public:
			MummerRecord(std::string query, bool reverse, std::string target, uint32_t targetStart, uint32_t queryStart, uint32_t size)
							:  query_(std::move(query)), reverse_(reverse), target_(std::move(target)),
								 targetStart_(targetStart), queryStart_(queryStart), size_(size) {

			}
			MummerRecord(std::string query, bool reverse, std::string line): query_(std::move(query)), reverse_(reverse){

				auto toks = tokenizeString(njh::trim(line), "whitespace");
				if(toks.size()!= 4){
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error " << "error, expected 4 items but got: " << toks.size() << "\n";
					throw std::runtime_error{ss.str()};
				}
				target_ = toks[0];
				targetStart_ = njh::StrToNumConverter::stoToNum<uint32_t>(toks[1]);
				queryStart_ =  njh::StrToNumConverter::stoToNum<uint32_t>(toks[2]);
				size_ =        njh::StrToNumConverter::stoToNum<uint32_t>(toks[3]);
			}
			std::string query_;
			bool reverse_{false};
			std::string target_;
			uint32_t targetStart_{
							std::numeric_limits<uint32_t>::max()}; /**< 1 based positing since this is from the file specs */
			uint32_t queryStart_{
							std::numeric_limits<uint32_t>::max()};/**<1 based positing since this is from the file specs */
			uint32_t size_{std::numeric_limits<uint32_t>::max()};

			[[nodiscard]] Bed3RecordCore genBed3()const{
				Bed3RecordCore ret{target_, targetStart_ -1, targetStart_ -1 + size_};
				MetaDataInName meta;
				meta.addMeta("query", query_);
				if(reverse_){
					meta.addMeta("queryStart", queryStart_ - size_);
					meta.addMeta("queryEnd", queryStart_);
				}else{
					meta.addMeta("queryStart", queryStart_ - 1);
					meta.addMeta("queryEnd", queryStart_ - 1 + size_);
				}
				ret.extraFields_.emplace_back(meta.createMetaName());
				return ret;
			}
			[[nodiscard]] Bed6RecordCore genBed6() const {
				Bed6RecordCore ret{target_, targetStart_ -1, targetStart_ -1 + size_, query_, static_cast<double>(size_), reverse_ ? '-' : '+'};
				MetaDataInName meta;
				meta.addMeta("query", query_);
				if(reverse_){
					meta.addMeta("queryStart", queryStart_ - size_);
					meta.addMeta("queryEnd", queryStart_);
				}else{
					meta.addMeta("queryStart", queryStart_ - 1);
					meta.addMeta("queryEnd", queryStart_ - 1 + size_);
				}
				ret.extraFields_.emplace_back(meta.createMetaName());
				return ret;
			}
		};
		explicit MummerTargetHits(std::string target): target_(std::move(target)){

		};
		std::string target_;
		std::vector<MummerRecord> hits_;
	};

	std::vector<MummerTargetHits> targets_;
};

class NucmerShowCoordsRecord{
public:

	explicit NucmerShowCoordsRecord(const std::string & line){
		auto toks = tokenizeString(line, "\t");

		if(toks.size() < 13){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "should be at least 13 toks, " << toks.size() << "\n";
			throw std::runtime_error{ss.str()};
		}
		refStart_ = njh::StrToNumConverter::stoToNum<uint32_t>(toks[0]);
		refEnd_ = njh::StrToNumConverter::stoToNum<uint32_t>(toks[1]);
		queryStart_ = njh::StrToNumConverter::stoToNum<uint32_t>(toks[2]);
		queryEnd_ = njh::StrToNumConverter::stoToNum<uint32_t>(toks[3]);
		refHitLen_ = njh::StrToNumConverter::stoToNum<uint32_t>(toks[4]);
		queryHitLen_ = njh::StrToNumConverter::stoToNum<uint32_t>(toks[5]);
		perId_ = njh::StrToNumConverter::stoToNum<double>(toks[6]);
		refFullLen_ = njh::StrToNumConverter::stoToNum<uint32_t>(toks[7]);
		queryFullLen_ = njh::StrToNumConverter::stoToNum<uint32_t>(toks[8]);
		refLenCov_ = njh::StrToNumConverter::stoToNum<double>(toks[9]);
		queryLenCov_ = njh::StrToNumConverter::stoToNum<double>(toks[10]);
		refName_ = toks[11];
		queryName_ = toks[12];

	}
	//13 toks
	//show-coords -T -l  -c -H out.delta
	//[S1]	[E1]	[S2]	[E2]	[LEN 1]	[LEN 2]	[% IDY]	[LEN R]	[LEN Q]	[COV R]	[COV Q]	[TAGS]
	// positions are 1 based
	//0 ref start
	//1 ref end
	//2 query start
	//3 query end
	//4 ref hit len
	//5 query hit len
	//6 identity
	//7 full ref length
	//8 full query length
	//9 ref coverage percentage
	//10 query coverage percentage
	//11 reference name
	//12 query name

	uint32_t refStart_{std::numeric_limits<uint32_t>::max()};
	uint32_t refEnd_{std::numeric_limits<uint32_t>::max()};
	uint32_t queryStart_{std::numeric_limits<uint32_t>::max()};
	uint32_t queryEnd_{std::numeric_limits<uint32_t>::max()};
	uint32_t refHitLen_{std::numeric_limits<uint32_t>::max()};
	uint32_t queryHitLen_{std::numeric_limits<uint32_t>::max()};
	double perId_{std::numeric_limits<double>::max()};
	uint32_t refFullLen_{std::numeric_limits<uint32_t>::max()};
	uint32_t queryFullLen_{std::numeric_limits<uint32_t>::max()};
	double refLenCov_{std::numeric_limits<double>::max()};
	double queryLenCov_{std::numeric_limits<double>::max()};
	std::string refName_;
	std::string queryName_;

	[[nodiscard]] bool reverseStrand() const{
		return queryEnd_ < queryStart_;
	}
	[[nodiscard]] Bed6RecordCore genBed6() const{
		MetaDataInName meta;
		meta.addMeta("perID", perId_);
		meta.addMeta("refLenCov", refLenCov_);
		meta.addMeta("queryLenCov", queryLenCov_);
		uint32_t actualStart = reverseStrand() ? queryEnd_ : queryStart_;
		--actualStart;
		uint32_t actualEnd = reverseStrand() ? queryStart_ : queryEnd_;
		meta.addMeta("actualStart", actualStart);
		meta.addMeta("actualEnd", actualEnd);
    meta.addMeta("queryName", queryName_);
		meta.addMeta("queryStart", queryStart_);
		meta.addMeta("queryEnd", queryEnd_);

//		meta.addMeta("queryStart", queryStart_);
//		meta.addMeta("queryEnd", queryEnd_);

		Bed6RecordCore ret(refName_,
											 refStart_ -1 ,
											 refEnd_ ,
//											 reverseStrand() ? refEnd_ -1 : refStart_ -1,
//											 reverseStrand() ? refStart_ : refEnd_,
											 njh::pasteAsStr(queryName_, "-",actualStart , "-", actualEnd),
											 uAbsdiff(refStart_, refEnd_) + 1,
											 reverseStrand() ? '-' : '+');
		ret.extraFields_.emplace_back(meta.createMetaName());
		return ret;
	}

};

int parsingFileExpRunner::parseMummberResultsToBed(const njh::progutils::CmdArgs & inputCommands) {

	IoOptions ioOpts;
	ioOpts.out_.outExtention_ = ".bed";
	seqSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processVerbose();
	setUp.processWritingOptions(ioOpts.out_);
	setUp.setOption(ioOpts.in_.inFilename_, "--mummerOut", "output from mummer", true);

	setUp.finishSetUp(std::cout);

	InputStream  in(ioOpts.in_);
	OutputStream out(ioOpts.out_);
	std::string line;
	bool currentlyRevStrand = false;
	std::string currentQuery;
	out << "#target\ttargetStart\ttargetEnd\tqueryName\tlength\tstrand\tmeta" << std::endl;
	while(njh::files::crossPlatGetline(in, line)){
		if(!njh::beginsWith(line, "#")){
			if(njh::beginsWith(line, ">")){
				//new record
				currentlyRevStrand = njh::endsWith(line, "Reverse");
				auto toks = tokenizeString(line, "whitespace");
				if(toks.size() < 2){
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error " << "expected record line to be at least 2 items, not: " << toks.size() << "\n";
					throw std::runtime_error{ss.str()};
				}
				currentQuery = toks[1];
			} else {
				MummerOutput::MummerTargetHits::MummerRecord record(currentQuery, currentlyRevStrand, line);
				out << record.genBed6().toDelimStrWithExtra() << std::endl;
			}
		}
	}
	return 0;
}

int parsingFileExpRunner::parseNucmerResultsToBed(const njh::progutils::CmdArgs & inputCommands) {
	bool renameWithHitCount = false;

	IoOptions ioOpts;
	ioOpts.out_.outExtention_ = ".bed";
	seqSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processVerbose();
	setUp.processWritingOptions(ioOpts.out_);
	setUp.setOption(ioOpts.in_.inFilename_, "--coordsOutput", "output from show-coords filename", true);
	setUp.setOption(renameWithHitCount, "--renameWithHitCount", "rename With Hit Count");

	setUp.finishSetUp(std::cout);

	InputStream  in(ioOpts.in_);
	OutputStream out(ioOpts.out_);
	std::string line;
	while(njh::files::crossPlatGetline(in, line)){
		NucmerShowCoordsRecord record(line);
		out << record.genBed6().toDelimStrWithExtra() << std::endl;
	}
	return 0;
}

} /* namespace njhseq */
