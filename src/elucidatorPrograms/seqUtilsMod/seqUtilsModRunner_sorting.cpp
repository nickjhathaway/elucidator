/*
 * seqUtilsModRunner_sorting.cpp
 *
 *  Created on: Oct 12, 2020
 *      Author: nick
 */




#include "seqUtilsModRunner.hpp"
#include <njhseq//objects/counters/DNABaseCounter.hpp>
#include <njhseq/IO/SeqIO/SeqIO.hpp>
#include <njhseq/objects/seqObjects/seqKmers.h>
#include <njhseq/readVectorManipulation/readVectorHelpers/readVecSorter.hpp>


namespace njhseq {




int seqUtilsModRunner::sortReadsByNameNaturalSort(const njh::progutils::CmdArgs & inputCommands) {

	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processDefaultReader(true);
  if (setUp.pars_.ioOptions_.out_.outFilename_ == "out") {
  	setUp.pars_.ioOptions_.out_.outFilename_ = njh::files::prependFileBasename(njh::files::removeExtension(setUp.pars_.ioOptions_.firstName_), "sorted_");
  }
	setUp.finishSetUp(std::cout);
	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	std::regex regPat{"([A-Za-z0-9\\.]+)"};
	std::regex subPat{"([A-Za-z]*)([0-9\\.]*)"};
	std::regex subPatDoubleFirst{"([0-9\\.]*)([A-Za-z]*)"};

	struct SeqWithNameSplit {
		SeqWithNameSplit(const seqInfo &seq, const std::regex &regPat, const std::regex &subPat, const std::regex &subPatDoubleFirst) :
						seqBase_(seq),
						regPat_(regPat),
						subPat_(subPat),
						subPatDoubleFirst_(subPatDoubleFirst),
						nameToks_(njh::tokStrOnMatchRegex(seqBase_.name_, regPat)) {
			for (const auto &name: nameToks_) {
				std::smatch nameMatch;
				std::smatch secondNameMatch;

				if (std::regex_match(name.begin(), name.end(), secondNameMatch, subPatDoubleFirst_)) {
//					std::cout << __FILE__ << " " << __LINE__ << std::endl;
//					std::cout << njh::bashCT::red;
//					std::cout << "name: " << name << std::endl;
//					std::cout << "secondNameMatch[1]: " << secondNameMatch[1] << std::endl;
//					std::cout << "secondNameMatch[2]: " << secondNameMatch[2] << std::endl;
//					std::cout << njh::bashCT::reset << std::endl;
					subNameToks_.emplace_back(std::make_pair(secondNameMatch[2],
																									 ("" == secondNameMatch[1] ? std::numeric_limits<double>::min() : std::stod(
																													 secondNameMatch[1]))));
				} else if (!std::regex_match(name.begin(), name.end(), nameMatch, subPat_)) {
//					std::stringstream ss;
//					ss << __PRETTY_FUNCTION__ << ", error " << "name didn't match pattern"<< "\n";
//					throw std::runtime_error{ss.str()};
//					std::cout << __FILE__ << " " << __LINE__ << std::endl;
//					std::cout << njh::bashCT::purple;
//					std::cout << "name: " << name << std::endl;
//					std::cout << njh::bashCT::reset << std::endl;
					subNameToks_.emplace_back(std::make_pair(name, std::numeric_limits<double>::min()));
				} else {
//					std::cout << __FILE__ << " " << __LINE__ << std::endl;
//					std::cout << njh::bashCT::blue;
//					std::cout << "name: " << name << std::endl;
//					std::cout << "nameMatch[1]: " << nameMatch[1] << std::endl;
//					std::cout << "nameMatch[2]: " << nameMatch[2] << std::endl;
//					std::cout << njh::bashCT::reset << std::endl;
					subNameToks_.emplace_back(std::make_pair(nameMatch[1],
																									 ("" == nameMatch[2] ? std::numeric_limits<double>::min() : std::stod(
																													 nameMatch[2]))));
				}
//				std::smatch nameMatch;
//				if (!std::regex_match(name.begin(), name.end(), nameMatch, subPat_)) {
////					std::stringstream ss;
////					ss << __PRETTY_FUNCTION__ << ", error " << "name didn't match pattern"<< "\n";
////					throw std::runtime_error{ss.str()};
//					subNameToks_.emplace_back(std::make_pair(name, std::numeric_limits<double>::min()));
//				} else if (!std::regex_match(name.begin(), name.end(), nameMatch, subPatDoubleFirst_)) {
//					std::cout << __FILE__ << " " << __LINE__ << std::endl;
//					std::cout << njh::bashCT::red;
//					std::cout << nameMatch[1] << std::endl;
//					std::cout << nameMatch[2] << std::endl;
//					std::cout << njh::bashCT::reset << std::endl;
//					subNameToks_.emplace_back(std::make_pair(nameMatch[2],
//																									 ("" == nameMatch[1] ? std::numeric_limits<double>::min() : std::stod(
//																													 nameMatch[1]))));
//				} else {
//					std::cout << __FILE__ << " " << __LINE__ << std::endl;
//					std::cout << njh::bashCT::blue;
//					std::cout << nameMatch[1] << std::endl;
//					std::cout << nameMatch[2] << std::endl;
//					std::cout << njh::bashCT::reset << std::endl;
//					subNameToks_.emplace_back(std::make_pair(nameMatch[1],
//																									 ("" == nameMatch[2] ? std::numeric_limits<double>::min() : std::stod(
//																													 nameMatch[2]))));
//				}

			}

//			std::cout << seqBase_.name_ << std::endl;
//			std::cout << "nameToks_.size():    " << nameToks_.size() << std::endl;
//			for(const auto & tok : nameToks_){
//				std::cout << "\t" << tok << std::endl;
//			}
//			std::cout << "subNameToks_.size(): " << subNameToks_.size() << std::endl;
//			for(const auto & subtok : subNameToks_){
//				std::cout << "\t" << subtok.first << ", " << subtok.second << std::endl;
//			}
		}

		seqInfo seqBase_;
		std::regex regPat_;
		std::regex subPat_;
		std::regex subPatDoubleFirst_;
		VecStr nameToks_;
		std::vector<std::pair<std::string, double>> subNameToks_;
	};
	std::vector<SeqWithNameSplit> seqs;
	seqInfo seq;
	while(reader.readNextRead(seq)){
		seqs.emplace_back(seq, regPat, subPat, subPatDoubleFirst);
	}

	njh::sort(seqs, []( const SeqWithNameSplit & seq1,  const SeqWithNameSplit & seq2){
		auto smallest = std::min(seq1.nameToks_.size(), seq2.nameToks_.size());

		for(uint32_t pos = 0; pos < smallest; ++pos){
			if(seq1.subNameToks_[pos].first == seq2.subNameToks_[pos].first){
				if(seq1.subNameToks_[pos].second != seq2.subNameToks_[pos].second){
					return seq1.subNameToks_[pos].second < seq2.subNameToks_[pos].second;
				}
			} else {
				return seq1.subNameToks_[pos].first < seq2.subNameToks_[pos].first;
			}
		}
		return seq1.subNameToks_.size() < seq2.subNameToks_.size();
	});

	SeqOutput::write(seqs, setUp.pars_.ioOptions_);

	return 0;
}

int seqUtilsModRunner::sortReadsByKmerEntropy(const njh::progutils::CmdArgs & inputCommands) {
	bool mark = false;
	uint32_t kLen = 2;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processDefaultReader(true);
  if (setUp.pars_.ioOptions_.out_.outFilename_ == "out") {
  	setUp.pars_.ioOptions_.out_.outFilename_ = njh::files::prependFileBasename(njh::files::removeExtension(setUp.pars_.ioOptions_.firstName_), "sorted_");
  }
	setUp.setOption(kLen, "--klen", "kmer length");
	setUp.setOption(mark, "--mark", "Add entropy to name");
	setUp.finishSetUp(std::cout);

	SeqInput reader(setUp.pars_.ioOptions_);
	auto inputSeqs = reader.readAllReads<seqInfo>();
	std::vector<seqWithKmerInfo> seqs;
	for(const auto & seq : inputSeqs){
		seqs.emplace_back(seqWithKmerInfo(seq, kLen, false));
	}
	if(mark){
		for(auto & seq : seqs){
			seq.seqBase_.name_.append(njh::pasteAsStr("[entropy=", seq.kInfo_.computeKmerEntropy(), "]"));
		}
	}
	njh::sort(seqs, []( const seqWithKmerInfo & seq1,  const seqWithKmerInfo & seq2){
		return seq1.kInfo_.computeKmerEntropy() < seq2.kInfo_.computeKmerEntropy();
	});

	SeqOutput::write(seqs, setUp.pars_.ioOptions_);

	return 0;
}

int seqUtilsModRunner::sortReadsByEntropy(const njh::progutils::CmdArgs & inputCommands) {
	bool mark = false;
	uint8_t minBase = 4;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processDefaultReader(true);
  if (setUp.pars_.ioOptions_.out_.outFilename_ == "out") {
  	setUp.pars_.ioOptions_.out_.outFilename_ = njh::files::prependFileBasename(njh::files::removeExtension(setUp.pars_.ioOptions_.firstName_), "sorted_");
  }

	setUp.setOption(minBase, "--minBase", "Min Base");
	setUp.setOption(mark, "--mark", "Add entropy to name");
	setUp.finishSetUp(std::cout);

	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	struct SeqInfoWithCountEntropy {
		SeqInfoWithCountEntropy(const seqInfo & seq, uint8_t minBase):seqBase_(seq){
			DNABaseCounter counter;
			counter.increase(seqBase_.seq_);
			entropy_ = counter.computeEntrophyBasedOffAlph(minBase);
		}
		seqInfo seqBase_;
		double entropy_;
	};

	std::vector<SeqInfoWithCountEntropy> seqs;
	seqInfo seq;
	while(reader.readNextRead(seq)){
	  readVec::handelLowerCaseBases(seq, setUp.pars_.ioOptions_.lowerCaseBases_);
	  if (setUp.pars_.ioOptions_.removeGaps_) {
	    readVec::removeGapsFromReads(seq);
	  }
		seqs.emplace_back(seq, minBase);
	}
	if(mark){
		for(auto & seq : seqs){
			seq.seqBase_.name_.append(njh::pasteAsStr("[entropy=", seq.entropy_, "]"));
		}
	}
	njh::sort(seqs, []( SeqInfoWithCountEntropy & seq1,  SeqInfoWithCountEntropy & seq2){
		return seq1.entropy_ < seq2.entropy_;
	});

	SeqOutput::write(seqs, setUp.pars_.ioOptions_);

	return 0;
}

int seqUtilsModRunner::sortReads(
		const njh::progutils::CmdArgs & inputCommands) {
	// parameters
	std::string sortBy = "averageError";
	bool decending = true;
	seqUtilsModSetUp setUp(inputCommands);
	setUp.setUpSortReads(sortBy, decending);
	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.in_.readAllReads<readObject>();
	if (stringToLowerReturn(sortBy) == "gc") {
		auto gcCompare = [](const readObject & read1, const readObject & read2) {
			if(read1.counter_.gcContent_ < read2.counter_.gcContent_) {
				return true;
			} else if (read1.counter_.gcContent_ == read2.counter_.gcContent_) {
				return read1.seqBase_.seq_ < read2.seqBase_.seq_;
			} else {
				return true;
			}
		};

		njh::for_each(inReads, [](readObject & read) {
			read.setLetterCount();
			read.counter_.calcGcContent();
		});
		njh::sort(inReads, gcCompare);
	} else {
		readVecSorter::sortReadVector(inReads, sortBy, decending);
	}
	reader.openWrite(inReads);
	setUp.logRunTime(std::cout);
	return 0;
}

int seqUtilsModRunner::sortReadsByKCompToTop(
		const njh::progutils::CmdArgs & inputCommands) {
	// parameters

	uint32_t klen = 4;
	// bool decending = true;
	seqUtilsModSetUp setUp(inputCommands);
  // input file info
	setUp.processDefaultReader(true);
  if (setUp.pars_.ioOptions_.out_.outFilename_ == "out") {
  	setUp.pars_.ioOptions_.out_.outFilename_ = njh::files::prependFileBasename(njh::files::removeExtension(setUp.pars_.ioOptions_.firstName_), "sorted_");
  }
  setUp.setOption(klen, "--klen", "klen");
  // bool ascending = false;
  // setUp.setOption(ascending, "--ascending", "Ascending Sort");
  // decending = !ascending;
  setUp.finishSetUp(std::cout);

	std::vector<std::shared_ptr<seqWithKmerInfo>> allInputReads;
	seqInfo seq;
	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	while (reader.readNextRead(seq)) {
		auto kSeq = std::make_shared<seqWithKmerInfo>(seq);
		readVec::handelLowerCaseBases(kSeq, setUp.pars_.ioOptions_.lowerCaseBases_);
		allInputReads.emplace_back(kSeq);
	}
	readVecSorter::sortByTotalCountAE(allInputReads, true);
	allSetKmers(allInputReads, klen, false);
	seqWithKmerInfo top = *allInputReads.front();
	readVecSorter::sortReadVectorFunc(allInputReads,[&top](const std::shared_ptr<seqWithKmerInfo>& p1, const std::shared_ptr<seqWithKmerInfo> & p2){
		return top.compareKmers(*p1).second > top.compareKmers(*p2).second;
	});
	reader.openWrite(allInputReads);
	setUp.logRunTime(std::cout);
	return 0;
}


int seqUtilsModRunner::sortReadsPairedEnd(
				const njh::progutils::CmdArgs & inputCommands) {
	// parameters
	seqUtilsModSetUp setUp(inputCommands);
	// input file info
	setUp.processDefaultReader(seqSetUp::pairedReadInFormatsAvailable_, true);
	if (setUp.pars_.ioOptions_.out_.outFilename_ == "out") {
		setUp.pars_.ioOptions_.out_.outFilename_ = njh::files::prependFileBasename(
						njh::files::removeExtension(setUp.pars_.ioOptions_.firstName_), "sorted_");
	}
	bool ascending = false;
	setUp.setOption(ascending, "--ascending", "Ascending Sort");

	setUp.finishSetUp(std::cout);
	std::vector<PairedRead> allInputReads = SeqInput::getSeqVec<PairedRead>(setUp.pars_.ioOptions_);
	readVecSorter::sortBySeq(allInputReads, !ascending);
	SeqOutput::write(allInputReads, setUp.pars_.ioOptions_);
	if (setUp.pars_.verbose_) {
		setUp.logRunTime(std::cout);
	}
	return 0;
}




} //namespace njhseq

