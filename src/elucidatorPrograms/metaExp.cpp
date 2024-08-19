/*
 * metaExp.cpp
 *
 *  Created on: Mar 14, 2018
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

#include "metaExp.hpp"
#include <njhseq/IO/SeqIO.h>
#include <njhseq/objects/Meta.h>


namespace njhseq {


metaExpRunner::metaExpRunner()
    : njh::progutils::ProgramRunner(
          {
					 addFunc("addMetaBySampleName", addMetaBySampleName, false),
					 addFunc("addMetaByMetaField", addMetaByMetaField, false),
					 addFunc("excludeSeqsFileWithNumericMetaCutOff", excludeSeqsFileWithNumericMetaCutOff, false),
					 addFunc("excludeSeqsFileWithMatchingMeta", excludeSeqsFileWithMatchingMeta, false),
					 addFunc("selectMetaFieldsToKeep", selectMetaFieldsToKeep, false),
					 addFunc("splitSeqFileWithMeta", splitSeqFileWithMeta, false),
					 addFunc("splitSeqFileWithExternalMeta", splitSeqFileWithExternalMeta, false),
					 addFunc("createTableFromSeqs", createTableFromSeqs, false),
					 addFunc("printMetaFieldsFromSeqs", printMetaFieldsFromSeqs, false),
					 addFunc("addMetaFieldToAll", addMetaFieldToAll, false),
					 addFunc("renameSeqsWithMetaField", renameSeqsWithMetaField, false),
					 addFunc("addSeqNameAsSampleMeta", addSeqNameAsSampleMeta, false),

           },//
          "metaExp") {}


int metaExpRunner::addSeqNameAsSampleMeta(const njh::progutils::CmdArgs & inputCommands){
	bool overWriteMeta = false;
	bool addLengthMeta = false;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(addLengthMeta, "--addLengthMeta", "Add length meta field, the length of the sequence");
	setUp.setOption(overWriteMeta, "--overWriteMeta", "Over Write Meta Fields in names");
	setUp.processDefaultReader(true);
	setUp.description_ = "Set sample meta field as the seq name";

	setUp.finishSetUp(std::cout);


	seqInfo seq;
	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	reader.openOut();
	while (reader.readNextRead(seq)) {
		MetaDataInName seqMeta;
		std::string sample = seq.name_;
		if (MetaDataInName::nameHasMetaData(seq.name_)){
			seqMeta = MetaDataInName(seq.name_);
			MetaDataInName::removeMetaDataInName(sample);
		}
		seqMeta.addMeta("sample", sample, overWriteMeta);
		if(addLengthMeta){
			seqMeta.addMeta("length", len(seq), overWriteMeta);
		}
		seqMeta.resetMetaInName(seq.name_);
		reader.write(seq);
	}
	//
	return 0;
}



int metaExpRunner::renameSeqsWithMetaField(const njh::progutils::CmdArgs & inputCommands){
	VecStr metaFields;
	bool makeNamesUnique = false;
	bool keepMeta = false;
	std::string sep = "-";
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(sep, "--sep", "Separator between meta fields if using multiple");
	setUp.setOption(metaFields, "--metaFields", "Meta fields within name to use", true);
	setUp.setOption(makeNamesUnique, "--makeNamesUnique", "Make Names unique by appending to the name a number if name already take");
	setUp.setOption(keepMeta, "--keepMeta", "Keep Meta fields used to rename");

	setUp.processDefaultReader(true);
  if (setUp.pars_.ioOptions_.out_.outFilename_ == "out") {
  	setUp.pars_.ioOptions_.out_.outFilename_ = njh::files::prependFileBasename(njh::files::removeExtension(setUp.pars_.ioOptions_.firstName_), "renamed_");
  }
	setUp.finishSetUp(std::cout);
	seqInfo seq;
	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	reader.openOut();
	std::unordered_map<std::string, uint32_t> nameCounts;
	while (reader.readNextRead(seq)) {
		if (MetaDataInName::nameHasMetaData(seq.name_)) {
			MetaDataInName meta(seq.name_);
			MetaDataInName keptMeta;
			VecStr metaData;
			for (const auto & field : metaFields) {
				metaData.emplace_back(meta.getMeta(field));
				keptMeta.addMeta(field, meta.getMeta(field));
			}
			auto newName = njh::conToStr(metaData, sep);
			if(makeNamesUnique){
				++nameCounts[newName];
				if(nameCounts[newName] > 1){
					newName = njh::pasteAsStr(newName, ".", nameCounts[newName]);
				}
			}
			if(keepMeta){
				newName+=keptMeta.createMetaName();
			}
			seq.name_ = newName;
		} else {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << seq.name_
					<< " doesn't contain meta data" << "\n";
			throw std::runtime_error { ss.str() };
		}
		reader.write(seq);
	}
	return 0;
}


int metaExpRunner::addMetaBySampleName(const njh::progutils::CmdArgs & inputCommands){
	bfs::path metaDataFnp = "";
	bool overWriteMeta = false;
	bool useName = false;
	bool addLengthMeta = false;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(addLengthMeta, "--addLengthMeta", "Add length meta field, the length of the sequence");
	setUp.setOption(useName, "--useName", "Just use the name of the seq rather than trying to guess sample name");
	setUp.setOption(metaDataFnp, "--metaData", "Name of the meta data field, must have a column named sample", true);
	setUp.setOption(overWriteMeta, "--overWriteMeta", "Over Write Meta Fields in names");
	setUp.processDefaultReader(true);
	setUp.description_ = "Set meta data in sequence names by matching either the meta filed sample or just by seq name";
	setUp.examples_.emplace_back("MASTERPROGRAM SUBPROGRAM --metaDataFnp meta.tab.txt --fasta infile.fasta");
	setUp.examples_.emplace_back("MASTERPROGRAM SUBPROGRAM --metaDataFnp meta.tab.txt --fasta infile.fasta --overWriteMeta #over write meta fields already present");
	setUp.examples_.emplace_back("MASTERPROGRAM SUBPROGRAM --metaDataFnp meta.tab.txt --fasta infile.fasta --out outfile.fasta");

	setUp.finishSetUp(std::cout);

	MultipleGroupMetaData groupMetaData(metaDataFnp);

	seqInfo seq;
	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	reader.openOut();
	while (reader.readNextRead(seq)) {
		if (MetaDataInName::nameHasMetaData(seq.name_) && !useName) {
			MetaDataInName meta(seq.name_);
			std::string possibleSampName = seq.name_;
			if (meta.containsMeta("sample")) {
				possibleSampName = meta.getMeta("sample");
			} else {
				MetaDataInName::removeMetaDataInName(possibleSampName);
			}
			if (groupMetaData.hasSample(possibleSampName)) {
				auto inputMeta = groupMetaData.getMetaForSample(possibleSampName,
						getVectorOfMapKeys(groupMetaData.groupData_));
				meta.addMeta(inputMeta, overWriteMeta);
				if(addLengthMeta){
					meta.addMeta("length",len(seq), overWriteMeta);
				}
				meta.resetMetaInName(seq.name_);
			} else {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error, couldn't find sample name, " << possibleSampName <<  " in "
						<< metaDataFnp
						<< " by looking for meta field sample or use name with meta removed"
						<< "\n";
				ss << "Options are: " << njh::conToStr(groupMetaData.samples_, ",") << "\n";
				throw std::runtime_error { ss.str() };
			}
		} else {
			if (groupMetaData.hasSample(seq.name_)) {
				auto inputMeta = groupMetaData.getMetaForSample(seq.name_,
						getVectorOfMapKeys(groupMetaData.groupData_));
				inputMeta.resetMetaInName(seq.name_);
				if(addLengthMeta){
					inputMeta.addMeta("length",len(seq), overWriteMeta);
				}
			} else {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error, couldn't find sample name " << seq.name_ << " in "
						<< metaDataFnp
						<< " by looking for meta field sample or use name with meta removed"
						<< "\n";
				ss << "Options are: " << njh::conToStr(groupMetaData.samples_, ",") << "\n";

				throw std::runtime_error { ss.str() };
			}
		}
		reader.write(seq);
	}
	//
	return 0;
}



int metaExpRunner::addMetaFieldToAll(const njh::progutils::CmdArgs & inputCommands){
	bool overWriteMeta = false;
	std::string metaTokenizer = ",";
	std::string metaField;
	std::string metaData;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(metaData, "--metaData", "The meta value to add for the field", true);
	setUp.setOption(metaField, "--metaField", "Name of the meta data field", true);
	setUp.setOption(overWriteMeta, "--overWriteMeta", "Over Write Meta Fields in names");
	setUp.setOption(metaTokenizer, "--metaTokenizer", "meta Tokenizer for adding multiple at once");

	setUp.processDefaultReader(true);
	setUp.description_ = "Set meta field to a set value for all sequences";
	setUp.examples_.emplace_back("MASTERPROGRAM SUBPROGRAM --metaData condition1 --metaField experiment --fasta infile.fasta");

	setUp.finishSetUp(std::cout);


	auto metaDataToks = tokenizeString(metaData, metaTokenizer);
	auto metaFieldToks = tokenizeString(metaField, metaTokenizer);
	if(metaDataToks.size() != metaFieldToks.size()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "metaDataToks.size(): " << metaDataToks.size() << " must equal metaFieldToks.size(): " << metaFieldToks.size() << "" << "\n";
		throw std::runtime_error{ss.str()};
	}

	seqInfo seq;
	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	reader.openOut();
	while (reader.readNextRead(seq)) {
		MetaDataInName seqMeta;
		if (MetaDataInName::nameHasMetaData(seq.name_)){
			seqMeta = MetaDataInName(seq.name_);
		}
		for(const auto & tokEnum : iter::enumerate(metaFieldToks)) {
			seqMeta.addMeta(tokEnum.element, metaDataToks[tokEnum.index], overWriteMeta);
		}
		if (MetaDataInName::nameHasMetaData(seq.name_)){
			seqMeta.resetMetaInName(seq.name_);
		}else{
			seq.name_.append(seqMeta.createMetaName());
		}
		reader.write(seq);
	}

	return 0;
}

int metaExpRunner::addMetaByMetaField(const njh::progutils::CmdArgs & inputCommands){
	bfs::path metaDataFnp = "";
	bool overWriteMeta = false;
	bool useName = false;
	bool addLengthMeta = false;
	std::string metaField;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(addLengthMeta, "--addLengthMeta", "Add length meta field, the length of the sequence");
	setUp.setOption(useName, "--useName", "Just use the name of the seq rather than trying to guess sample name");
	setUp.setOption(metaDataFnp, "--metaData", "Name of the meta data file, must have a column named for the given meta field", true);
	setUp.setOption(metaField, "--metaField", "Name of the meta data field", true);

	setUp.setOption(overWriteMeta, "--overWriteMeta", "Over Write Meta Fields in names");
	setUp.processDefaultReader(true);
	setUp.description_ = "Set meta data in sequence names by matching either the meta filed sample or just by seq name";
	setUp.examples_.emplace_back("MASTERPROGRAM SUBPROGRAM --metaDataFnp meta.tab.txt --fasta infile.fasta");
	setUp.examples_.emplace_back("MASTERPROGRAM SUBPROGRAM --metaDataFnp meta.tab.txt --fasta infile.fasta --overWriteMeta #over write meta fields already present");
	setUp.examples_.emplace_back("MASTERPROGRAM SUBPROGRAM --metaDataFnp meta.tab.txt --fasta infile.fasta --out outfile.fasta");

	setUp.finishSetUp(std::cout);

	table inTab(metaDataFnp, "\t", true);
	inTab.checkForColumnsThrow(VecStr{metaField},__PRETTY_FUNCTION__);
	auto values = inTab.getColumnLevels(metaField);
	std::set<std::string> valuesSet(values.begin(), values.end());
	inTab.columnNames_[inTab.getColPos(metaField)] = "sample";
	inTab.setColNamePositions();
	MultipleGroupMetaData groupMetaData(inTab, valuesSet);

	seqInfo seq;
	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	reader.openOut();
	while (reader.readNextRead(seq)) {
		if (MetaDataInName::nameHasMetaData(seq.name_) && !useName) {
			MetaDataInName meta(seq.name_);
			std::string possibleSampName = seq.name_;
			if (meta.containsMeta(metaField)) {
				possibleSampName = meta.getMeta(metaField);
			} else {
				MetaDataInName::removeMetaDataInName(possibleSampName);
			}
			if (groupMetaData.hasSample(possibleSampName)) {
				auto inputMeta = groupMetaData.getMetaForSample(possibleSampName,
						getVectorOfMapKeys(groupMetaData.groupData_));
				meta.addMeta(inputMeta, overWriteMeta);
				if(addLengthMeta){
					meta.addMeta("length",len(seq), overWriteMeta);
				}
				meta.resetMetaInName(seq.name_);
			} else {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error, couldn't find sample name, " << possibleSampName <<  " in "
						<< metaDataFnp
						<< " by looking for meta field sample or use name with meta removed"
						<< "\n";
				ss << "Options are: " << njh::conToStr(groupMetaData.samples_, ",") << "\n";
				throw std::runtime_error { ss.str() };
			}
		} else {
			if (groupMetaData.hasSample(seq.name_)) {
				auto inputMeta = groupMetaData.getMetaForSample(seq.name_,
						getVectorOfMapKeys(groupMetaData.groupData_));
				inputMeta.resetMetaInName(seq.name_);
				if(addLengthMeta){
					inputMeta.addMeta("length",len(seq), overWriteMeta);
				}
			} else {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error, couldn't find sample name " << seq.name_ << " in "
						<< metaDataFnp
						<< " by looking for meta field sample or use name with meta removed"
						<< "\n";
				ss << "Options are: " << njh::conToStr(groupMetaData.samples_, ",") << "\n";

				throw std::runtime_error { ss.str() };
			}
		}
		reader.write(seq);
	}

	return 0;
}

int metaExpRunner::selectMetaFieldsToKeep(const njh::progutils::CmdArgs & inputCommands) {
	std::string metaFields ;
	seqSetUp setUp(inputCommands);
	setUp.description_ = "Take a sequence file with meta data in names and keep just the given fields";
	setUp.processDebug();
	setUp.processVerbose();
	setUp.pars_.ioOptions_.out_.outFilename_ = "";
	setUp.processDefaultReader(VecStr { "-fastq", "-fasta", "--fastq1" }, true);
	setUp.setOption(metaFields, "--metaFields",
			"The meta fields to keep, separated by a comma", true);

	setUp.finishSetUp(std::cout);

	auto metaToks = getInputValues(metaFields, ",");

	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	std::set<std::string> inputMetaFields;
	if (reader.ioOptions_.isPairedIn()) {
		PairedRead seq;
		while (reader.readNextRead(seq)) {
			if (MetaDataInName::nameHasMetaData(seq.seqBase_.name_)) {
				MetaDataInName seqMeta(seq.seqBase_.name_);
				njh::addVecToSet(getVectorOfMapKeys(seqMeta.meta_), inputMetaFields);
			}
		}
	} else {
		seqInfo seq;
		while (reader.readNextRead(seq)) {
			if (MetaDataInName::nameHasMetaData(seq.name_)) {
				MetaDataInName seqMeta(seq.name_);
				njh::addVecToSet(getVectorOfMapKeys(seqMeta.meta_), inputMetaFields);
			}
		}
	}

	VecStr missingMeta;
	for(const auto & meta : metaToks ){
		if(!njh::in(meta, inputMetaFields)){
			missingMeta.emplace_back(meta)	;
		}
	}
	if(!missingMeta.empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", missing " << njh::conToStr(missingMeta, ",") << "\n";
		throw std::runtime_error{ss.str()};
	}

	reader.in_.reOpenIn();
	if (reader.ioOptions_.isPairedIn()) {
		auto seqOpts = SeqIOOptions(setUp.pars_.ioOptions_.out_.outFilename_.string(), setUp.pars_.ioOptions_.outFormat_);
		seqOpts.out_.transferOverwriteOpts(setUp.pars_.ioOptions_.out_);
		SeqOutput pairWriter(seqOpts);
		pairWriter.openOut();
		PairedRead seq;
		while(reader.readNextRead(seq)){
			MetaDataInName outputMeta;
			VecStr missingFields;
			if(MetaDataInName::nameHasMetaData(seq.seqBase_.name_)){
				MetaDataInName seqMeta(seq.seqBase_.name_);
				for(const auto & m : metaToks){
					if(njh::in(m, seqMeta.meta_)){
						outputMeta.addMeta(m, seqMeta.getMeta(m));
					}else{
						missingFields.emplace_back(m);
					}
				}
			}else{
				missingFields = metaToks;
			}
			for(const auto & m : missingFields){
				outputMeta.addMeta(m, "NA");
			}
			outputMeta.resetMetaInName(seq.seqBase_.name_);
			outputMeta.resetMetaInName(seq.mateSeqBase_.name_);
			pairWriter.write(seq);
		}
	} else {
		auto seqOpts = SeqIOOptions(setUp.pars_.ioOptions_.out_.outFilename_.string(), setUp.pars_.ioOptions_.outFormat_);
		seqOpts.out_.transferOverwriteOpts(setUp.pars_.ioOptions_.out_);
		SeqOutput writer(seqOpts);
		writer.openOut();
		seqInfo seq;
		while(reader.readNextRead(seq)){
			MetaDataInName outputMeta;
			VecStr missingFields;
			if(MetaDataInName::nameHasMetaData(seq.name_)){
				MetaDataInName seqMeta(seq.name_);
				for(const auto & m : metaToks){
					if(njh::in(m, seqMeta.meta_)){
						outputMeta.addMeta(m, seqMeta.getMeta(m));
					}else{
						missingFields.emplace_back(m);
					}
				}
			}else{
				missingFields = metaToks;
			}
			for(const auto & m : missingFields){
				outputMeta.addMeta(m, "NA");
			}
			outputMeta.resetMetaInName(seq.name_);
			writer.write(seq);
		}
	}
	return 0;
}

int metaExpRunner::splitSeqFileWithMeta(const njh::progutils::CmdArgs & inputCommands) {
	std::string metaField;
	seqSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processVerbose();
	setUp.pars_.ioOptions_.out_.outFilename_ = "";
	setUp.processDefaultReader(VecStr { "--fastq", "--fasta", "--fastq1", "--fastqgz", "--fastagz", "--fastq1gz" }, true);
	setUp.setOption(metaField, "--metaField",
			"Meta Field to split on, could be multiple separated by a comma", true);

	setUp.finishSetUp(std::cout);

	auto metaToks = getInputValues(metaField, ",");

	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	std::set<std::string> inputMetaFields;
	if (reader.ioOptions_.isPairedIn()) {
		PairedRead seq;
		while (reader.readNextRead(seq)) {
			if (MetaDataInName::nameHasMetaData(seq.seqBase_.name_)) {
				MetaDataInName seqMeta(seq.seqBase_.name_);
				njh::addVecToSet(getVectorOfMapKeys(seqMeta.meta_), inputMetaFields);
			}
		}
	} else {
		seqInfo seq;
		while (reader.readNextRead(seq)) {
			if (MetaDataInName::nameHasMetaData(seq.name_)) {
				MetaDataInName seqMeta(seq.name_);
				njh::addVecToSet(getVectorOfMapKeys(seqMeta.meta_), inputMetaFields);
			}
		}
	}

	VecStr missingMeta;
	for(const auto & meta : metaToks ){
		if(!njh::in(meta, inputMetaFields)){
			missingMeta.emplace_back(meta)	;
		}
	}
	if(!missingMeta.empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", missing " << njh::conToStr(missingMeta, ",") << "\n";
		throw std::runtime_error{ss.str()};
	}

	reader.in_.reOpenIn();

	if (reader.ioOptions_.isPairedIn()) {
		MultiSeqOutCache<PairedRead> writers;
		PairedRead seq;
		while(reader.readNextRead(seq)){
			std::string uid;
			if(MetaDataInName::nameHasMetaData(seq.seqBase_.name_)){
				MetaDataInName seqMeta(seq.seqBase_.name_);
				for(const auto & m : metaToks){
					if(!uid.empty()){
						uid += "_";
					}
					if(njh::in(m, seqMeta.meta_)){
						uid += seqMeta.meta_[m];
					}else{
						uid += "NA";
					}
				}
			}else{
				uid = njh::conToStr(VecStr(metaToks.size(), "NA"), "_");
			}
			if(!writers.containsReader(uid)){
				auto seqOpts = SeqIOOptions(setUp.pars_.ioOptions_.out_.outFilename_.string() + "_" + uid, setUp.pars_.ioOptions_.outFormat_) ;
				if(setUp.pars_.ioOptions_.out_.outFilename_.empty()){
					seqOpts = SeqIOOptions(uid, setUp.pars_.ioOptions_.outFormat_);
				}
				seqOpts.out_.transferOverwriteOpts(setUp.pars_.ioOptions_.out_);
				writers.addReader(uid, seqOpts);
			}
			writers.add(uid, seq);
		}
	} else {
		MultiSeqOutCache<seqInfo> writers;
		seqInfo seq;
		while(reader.readNextRead(seq)){
			std::string uid;
			if(MetaDataInName::nameHasMetaData(seq.name_)){
				MetaDataInName seqMeta(seq.name_);
				for(const auto & m : metaToks){
					if(!uid.empty()){
						uid += "_";
					}
					if(njh::in(m, seqMeta.meta_)){
						uid += seqMeta.meta_[m];
					}else{
						uid += "NA";
					}
				}
			}else{
				uid = njh::conToStr(VecStr(metaToks.size(), "NA"), "_");
			}
			if(!writers.containsReader(uid)){
				auto opts = SeqIOOptions(setUp.pars_.ioOptions_.out_.outFilename_.string() + "_" + uid, setUp.pars_.ioOptions_.outFormat_);
				if(setUp.pars_.ioOptions_.out_.outFilename_.empty()){
					opts = SeqIOOptions(uid, setUp.pars_.ioOptions_.outFormat_);
				}
				opts.out_.transferOverwriteOpts(setUp.pars_.ioOptions_.out_);
				writers.addReader(uid, opts );
			}
			writers.add(uid, seq);
		}
	}


	return 0;
}

int metaExpRunner::excludeSeqsFileWithNumericMetaCutOff(const njh::progutils::CmdArgs & inputCommands) {
	std::string metaField;
	double metaValue = std::numeric_limits<double>::max();
	bool aboveCutOff = false;
	seqSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processVerbose();
	setUp.pars_.ioOptions_.out_.outFilename_ = "out";
	setUp.processDefaultReader(VecStr { "-fastq", "-fasta", "-fastqgz", "-fastagz", "--fastq1", "--fastq1gz" }, true);
	setUp.setOption(metaField, "--metaField",
			"Meta Field compare", true);
	setUp.setOption(metaValue, "--metaValue",
			"Meta value to compare to, exclusion be equal to this value and below (default) or above (set with --aboveCutOff)", true);
	setUp.setOption(aboveCutOff, "--aboveCutOff",
			"Exclude above cut off instead of the default excluding below cut off");
	setUp.finishSetUp(std::cout);


	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();

	OutOptions includeOpts(bfs::path(setUp.pars_.ioOptions_.out_.outFilename_.string() + "_include" + SeqIOOptions::getOutExtension(setUp.pars_.ioOptions_.outFormat_)));
	OutOptions excludeOpts(bfs::path(setUp.pars_.ioOptions_.out_.outFilename_.string() + "_exclude" + SeqIOOptions::getOutExtension(setUp.pars_.ioOptions_.outFormat_)));

	includeOpts.transferOverwriteOpts(setUp.pars_.ioOptions_.out_);
	excludeOpts.transferOverwriteOpts(setUp.pars_.ioOptions_.out_);

	SeqIOOptions includeSeqOpts(includeOpts, setUp.pars_.ioOptions_.outFormat_);
	SeqIOOptions excludeSeqOpts(excludeOpts, setUp.pars_.ioOptions_.outFormat_);

	SeqOutput includeWriter(includeSeqOpts);
	SeqOutput excludeWriter(excludeSeqOpts);
	includeWriter.openOut();
	excludeWriter.openOut();

	std::function<bool(const MetaDataInName &)> exclusionCheck =
			[&metaValue,&metaField](const MetaDataInName & mvalues) {
				return "NA" == mvalues.getMeta(metaField) || mvalues.getMeta<double>(metaField) <=metaValue;
			};

	if (aboveCutOff) {
		exclusionCheck = [&metaValue,&metaField](const MetaDataInName & mvalues) {
			return "NA" == mvalues.getMeta(metaField) || mvalues.getMeta<double>(metaField) >=metaValue;
		};
	}


	if (reader.ioOptions_.isPairedIn()) {
		PairedRead seq;
		while (reader.readNextRead(seq)) {
			if (MetaDataInName::nameHasMetaData(seq.seqBase_.name_)) {
				MetaDataInName seqMeta(seq.seqBase_.name_);
				if(seqMeta.containsMeta(metaField) && exclusionCheck(seqMeta) ){
					excludeWriter.write(seq);
				}else{
					includeWriter.write(seq);
				}
			}else{
				includeWriter.write(seq);
			}
		}
	} else {
		seqInfo seq;
		while (reader.readNextRead(seq)) {
			if (MetaDataInName::nameHasMetaData(seq.name_)) {
				MetaDataInName seqMeta(seq.name_);
				if(seqMeta.containsMeta(metaField) && exclusionCheck(seqMeta) ){
					excludeWriter.write(seq);
				}else{
					includeWriter.write(seq);
				}
			}else{
				includeWriter.write(seq);
			}
		}
	}



	return 0;
}

int metaExpRunner::excludeSeqsFileWithMatchingMeta(const njh::progutils::CmdArgs & inputCommands) {
	std::string metaField;
	std::string metaValue;
	seqSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processVerbose();
	setUp.pars_.ioOptions_.out_.outFilename_ = "out";
	setUp.processDefaultReader(VecStr { "-fastq", "-fasta", "--fastq1" }, true);
	setUp.setOption(metaField, "--metaField",
			"Meta Field compare", true);
	setUp.setOption(metaValue, "--metaValue",
			"Meta value to match, could be multiple separated by a comma", true);
	setUp.finishSetUp(std::cout);

	auto metaToks = getInputValues(metaValue, ",");

	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();

	OutOptions includeOpts(bfs::path(setUp.pars_.ioOptions_.out_.outFilename_.string() + "_include" + SeqIOOptions::getOutExtension(setUp.pars_.ioOptions_.outFormat_)));
	OutOptions excludeOpts(bfs::path(setUp.pars_.ioOptions_.out_.outFilename_.string() + "_exclude" + SeqIOOptions::getOutExtension(setUp.pars_.ioOptions_.outFormat_)));

	includeOpts.transferOverwriteOpts(setUp.pars_.ioOptions_.out_);
	excludeOpts.transferOverwriteOpts(setUp.pars_.ioOptions_.out_);

	SeqIOOptions includeSeqOpts(includeOpts, setUp.pars_.ioOptions_.outFormat_);
	SeqIOOptions excludeSeqOpts(excludeOpts, setUp.pars_.ioOptions_.outFormat_);

	SeqOutput includeWriter(includeSeqOpts);
	SeqOutput excludeWriter(excludeSeqOpts);
	includeWriter.openOut();
	excludeWriter.openOut();


	if (reader.ioOptions_.isPairedIn()) {
		PairedRead seq;
		while (reader.readNextRead(seq)) {
			if (MetaDataInName::nameHasMetaData(seq.seqBase_.name_)) {
				MetaDataInName seqMeta(seq.seqBase_.name_);
				if(seqMeta.containsMeta(metaField) && njh::in(seqMeta.getMeta(metaField),metaToks)){
					excludeWriter.write(seq);
				}else{
					includeWriter.write(seq);
				}
			}else{
				includeWriter.write(seq);
			}
		}
	} else {
		seqInfo seq;
		while (reader.readNextRead(seq)) {
			if (MetaDataInName::nameHasMetaData(seq.name_)) {
				MetaDataInName seqMeta(seq.name_);
				if(seqMeta.containsMeta(metaField) && njh::in(seqMeta.getMeta(metaField),metaToks)){
					excludeWriter.write(seq);
				}else{
					includeWriter.write(seq);
				}
			}else{
				includeWriter.write(seq);
			}
		}
	}



	return 0;
}

int metaExpRunner::splitSeqFileWithExternalMeta(const njh::progutils::CmdArgs & inputCommands) {
	std::string metaField;
	bfs::path metaFnp = "";
	seqSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processVerbose();
	setUp.pars_.ioOptions_.out_.outFilename_ = "";
	setUp.processDefaultReader(VecStr { "-fastq", "-fasta", "--fastq1" }, true);
	setUp.setOption(metaField, "--metaField",
			"Meta Field to split on, could be multiple separated by a comma", true);
	setUp.setOption(metaFnp, "--metaFnp",
				"Meta file, read names should match to the sample column", true);
	setUp.finishSetUp(std::cout);


	auto metaToks = getInputValues(metaField, ",");
	MultipleGroupMetaData externalMeta(metaFnp);
	externalMeta.checkForFieldsThrow(metaToks);
	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	std::set<std::string> inputMetaFields;
	if (reader.ioOptions_.isPairedIn()) {
		PairedRead seq;
		while (reader.readNextRead(seq)) {
			MetaDataInName seqMeta = externalMeta.getMetaForSample(seq.seqBase_.name_, metaToks);
			njh::addVecToSet(getVectorOfMapKeys(seqMeta.meta_), inputMetaFields);
		}
	} else {
		seqInfo seq;
		while (reader.readNextRead(seq)) {
			MetaDataInName seqMeta = externalMeta.getMetaForSample(seq.name_, metaToks);
			njh::addVecToSet(getVectorOfMapKeys(seqMeta.meta_), inputMetaFields);
		}
	}

	VecStr missingMeta;
	for(const auto & meta : metaToks ){
		if(!njh::in(meta, inputMetaFields)){
			missingMeta.emplace_back(meta)	;
		}
	}
	if(!missingMeta.empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", missing " << njh::conToStr(missingMeta, ",") << "\n";
		throw std::runtime_error{ss.str()};
	}

	reader.in_.reOpenIn();

	if (reader.ioOptions_.isPairedIn()) {
		MultiSeqOutCache<PairedRead> writers;
		PairedRead seq;
		while(reader.readNextRead(seq)){
			std::string uid;
			MetaDataInName seqMeta = externalMeta.getMetaForSample(seq.seqBase_.name_, metaToks);
			for(const auto & m : metaToks){
				if(!uid.empty()){
					uid += "_";
				}
				if(njh::in(m, seqMeta.meta_)){
					uid += seqMeta.meta_[m];
				}else{
					uid += "NA";
				}
			}
			if(!writers.containsReader(uid)){
				writers.addReader(uid, SeqIOOptions(setUp.pars_.ioOptions_.out_.outFilename_.string() + "_" + uid, setUp.pars_.ioOptions_.outFormat_) );
			}
			writers.add(uid, seq);
		}
	} else {
		MultiSeqOutCache<seqInfo> writers;
		seqInfo seq;
		while(reader.readNextRead(seq)){
			std::string uid;
			MetaDataInName seqMeta = externalMeta.getMetaForSample(seq.name_, metaToks);
			for(const auto & m : metaToks){
				if(!uid.empty()){
					uid += "_";
				}
				if(njh::in(m, seqMeta.meta_)){
					uid += seqMeta.meta_[m];
				}else{
					uid += "NA";
				}
			}
			if(!writers.containsReader(uid)){
				auto opts = SeqIOOptions(setUp.pars_.ioOptions_.out_.outFilename_.string() + "_" + uid, setUp.pars_.ioOptions_.outFormat_);
				if("" == setUp.pars_.ioOptions_.out_.outFilename_){
					opts = SeqIOOptions(uid, setUp.pars_.ioOptions_.outFormat_);
				}
				opts.out_.transferOverwriteOpts(setUp.pars_.ioOptions_.out_);
				writers.addReader(uid, opts );
			}
			writers.add(uid, seq);
		}
	}


	return 0;
}

int metaExpRunner::printMetaFieldsFromSeqs(const njh::progutils::CmdArgs & inputCommands) {
	auto tabOpts = TableIOOpts::genTabFileOut("", true);
	std::string fields;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processReadInNames(true);
	setUp.processWritingOptions(tabOpts.out_);
	setUp.finishSetUp(std::cout);

	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	OutputStream out(tabOpts.out_);
	seqInfo seq;
	std::set<std::string> allMetaKeys;
	while(reader.readNextRead(seq)){
		if(MetaDataInName::nameHasMetaData(getSeqBase(seq).name_)){
			MetaDataInName metaData(getSeqBase(seq).name_);
			for(const auto & meta : metaData.meta_){
				allMetaKeys.emplace(meta.first);
			}
		}
	}

	out << njh::conToStr(allMetaKeys, "\t") << std::endl;

	return 0;
}



int metaExpRunner::createTableFromSeqs(const njh::progutils::CmdArgs & inputCommands) {
	auto tabOpts = TableIOOpts::genTabFileOut("", true);
	std::string fields;
	std::string excludeFields;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processReadInNames(true);
	setUp.processWritingOptions(tabOpts.out_);
	setUp.setOption(fields, "--fields", "Only export these fields");
	setUp.setOption(excludeFields, "--excludeFields", "Exclude these fields");
	setUp.finishSetUp(std::cout);

	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<readObject>();

	std::set<std::string> allMetaKeys;
	for(const auto & seq : inReads){
		if(MetaDataInName::nameHasMetaData(getSeqBase(seq).name_)){
			MetaDataInName metaData(getSeqBase(seq).name_);
			for(const auto & meta : metaData.meta_){
				allMetaKeys.emplace(meta.first);
			}
		}
	}

	allMetaKeys.emplace("name");
	allMetaKeys.emplace("seq");

	OutputStream out(tabOpts.out_);
	table outTab;
	auto fieldToks = tokenizeString(fields, ",");
	auto excludeFieldToks = tokenizeString(excludeFields, ",");
	if(!fields.empty() || !excludeFields.empty()){
		std::set<std::string> exportFields;
		for(const auto & col : allMetaKeys){
			if( (fieldToks.empty()        ||  njh::in(col, fieldToks) ) &&
					(excludeFieldToks.empty() || !njh::in(col, excludeFieldToks) ) ){
				exportFields.emplace(col);
				outTab = seqsToMetaTable(inReads, exportFields);
			}
		}
	} else {
		outTab = seqsToMetaTable(inReads);
	}
	outTab.outPutContents(out, "\t");
	return 0;
}



} // namespace njhseq

