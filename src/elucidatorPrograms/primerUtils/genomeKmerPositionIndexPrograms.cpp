//
// Created by Nicholas Hathaway on 11/2/26.
//
// Build and query the on-disk k-mer position index: a primer only primes where its 3'-terminal bases match
// exactly, so the places it could prime are a lookup rather than an alignment search.
//

#include "primerUtilsRunner.hpp"

#include "elucidator/objects/MiscUtility/GenomeKmerPositionIndex.hpp"

#include <njhseq/IO/SeqIO.h>


namespace njhseq {

int primerUtilsRunner::buildGenomeKmerPositionIndex(const njh::progutils::CmdArgs &inputCommands) {
	bfs::path genomeFnp;
	bfs::path indexFnp;
	uint32_t kmerLength = 12;
	double buildBufferGb = 4;
	bool overWrite = false;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(genomeFnp, "--genomeFnp", "genome fasta to index", true);
	setUp.setOption(kmerLength, "--kmerLength", "k to build at. An index built at k answers every query from 1 to k, and since the position array holds one entry per genome position whatever k is, a larger k costs only a larger offset table");
	setUp.setOption(indexFnp, "--indexFnp", "where to write the index; defaults to sitting beside the genome, the same way the blast database does");
	setUp.setOption(buildBufferGb, "--buildBufferGb", "how much memory the build may use to assemble the position array, which it writes one sequential block at a time. Each block costs one sweep of the genome, so a larger buffer means fewer sweeps; peak memory is about this plus the offset table");
	setUp.setOption(overWrite, "--overWrite", "overwrite an index that is already there");
	setUp.finishSetUp(std::cout);

	if (indexFnp.empty()) {
		indexFnp = GenomeKmerPositionIndex::genIndexFnpForGenome(genomeFnp, kmerLength);
	}
	if (bfs::exists(indexFnp) && !overWrite) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << indexFnp << " already exists, use --overWrite to replace it"
				<< "\n";
		throw std::runtime_error{ss.str()};
	}

	njh::stopWatch watch;
	watch.setLapName("building");
	GenomeKmerPositionIndex::build(genomeFnp, indexFnp, kmerLength,
																		 static_cast<uint64_t>(buildBufferGb * 1073741824.0), setUp.pars_.verbose_);

	//open it again and report what it holds, which is also a check that it reads back
	GenomeKmerPositionIndex index(indexFnp);
	std::cout << "wrote " << indexFnp << std::endl;
	std::cout << "\t" << "k: " << index.k() << std::endl;
	std::cout << "\t" << "chromosomes: " << index.chroms().size() << std::endl;
	std::cout << "\t" << "indexed positions: " << index.numPositions() << std::endl;
	std::cout << "\t" << "size on disk: " << (bfs::file_size(indexFnp) / 1073741824.0) << " GB" << std::endl;
	if (setUp.pars_.verbose_) {
		watch.logLapTimes(std::cout, true, 6, true);
	}
	return 0;
}

int primerUtilsRunner::queryGenomeKmerPositionIndex(const njh::progutils::CmdArgs &inputCommands) {
	bfs::path indexFnp;
	bfs::path genomeFnp;
	uint32_t kmerLength = 12;
	std::string query;
	bfs::path queryFastaFnp;
	uint32_t threePrimeLen = 0;
	bool countOnly = false;
	bool alsoRevComp = false;
	bool verify = false;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(genomeFnp, "--genomeFnp", "genome the index was built from; used to find the index and, with --verify, to check the hits");
	setUp.setOption(kmerLength, "--kmerLength", "k the index was built at, used to find it beside the genome");
	setUp.setOption(indexFnp, "--indexFnp", "the index to query, if not the one beside --genomeFnp");
	setUp.setOption(query, "--query", "a single sequence to look up");
	setUp.setOption(queryFastaFnp, "--queryFastaFnp", "a fasta of sequences to look up");
	setUp.setOption(threePrimeLen, "--threePrimeLen", "look up only this many 3'-terminal bases of each query, which is the part that has to match for a primer to prime; 0 uses the whole query");
	setUp.setOption(countOnly, "--countOnly", "report how many hits each query has without reading the positions");
	setUp.setOption(alsoRevComp, "--alsoRevComp", "also look up the reverse complement, since only the forward strand is indexed");
	setUp.setOption(verify, "--verify", "read the genome back and confirm every reported hit really carries the query, and that no occurrence was missed");
	setUp.finishSetUp(std::cout);

	if (indexFnp.empty()) {
		if (genomeFnp.empty()) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error supply --indexFnp or --genomeFnp" << "\n";
			throw std::runtime_error{ss.str()};
		}
		indexFnp = GenomeKmerPositionIndex::genIndexFnpForGenome(genomeFnp, kmerLength);
	}
	GenomeKmerPositionIndex index(indexFnp);

	std::vector<std::pair<std::string, std::string>> queries;  //name, sequence
	if (!query.empty()) {
		queries.emplace_back(query, query);
	}
	if (!queryFastaFnp.empty()) {
		SeqInput reader(SeqIOOptions::genFastaIn(queryFastaFnp));
		reader.openIn();
		seqInfo seq;
		while (reader.readNextRead(seq)) {
			queries.emplace_back(seq.name_, seq.seq_);
		}
	}
	if (queries.empty()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error supply --query or --queryFastaFnp" << "\n";
		throw std::runtime_error{ss.str()};
	}

	//trim to the 3' end and, when asked, add the reverse complement as its own query
	std::vector<std::pair<std::string, std::string>> toLookUp;
	for (const auto &q: queries) {
		auto seq = q.second;
		stringToUpper(seq);  //the index is case insensitive, so normalise before trimming and looking up
		if (threePrimeLen > 0 && seq.size() > threePrimeLen) {
			seq = seq.substr(seq.size() - threePrimeLen);
		}
		toLookUp.emplace_back(q.first, seq);
		if (alsoRevComp) {
			toLookUp.emplace_back(q.first + "_revComp", seqUtil::reverseComplement(seq, "DNA"));
		}
	}

	if (countOnly) {
		std::cout << "name\tquery\thits" << std::endl;
		for (const auto &q: toLookUp) {
			std::cout << q.first << "\t" << q.second << "\t" << index.count(q.second) << std::endl;
		}
		return 0;
	}

	std::unordered_map<std::string, std::vector<GenomeKmerPositionIndex::Hit>> hitsByQuery;
	std::cout << "name\tquery\tchrom\tpos" << std::endl;
	for (const auto &q: toLookUp) {
		auto hits = index.find(q.second);
		for (const auto &hit: hits) {
			std::cout << q.first << "\t" << q.second << "\t" << hit.chrom_ << "\t" << hit.pos_ << std::endl;
		}
		hitsByQuery[q.second] = std::move(hits);
	}

	if (verify) {
		if (genomeFnp.empty()) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error --verify needs --genomeFnp" << "\n";
			throw std::runtime_error{ss.str()};
		}
		//independently scan the genome for each query and require the two answers to agree exactly, which
		//catches both a wrong position and a missed occurrence. The scan is case insensitive because the index
		//is: a soft-masked genome marks repeat content in lower case, and a primer primes there like anywhere
		//else, so those positions belong in the index and must be matched here too.
		std::unordered_map<std::string, std::vector<std::pair<std::string, uint64_t>>> truth;
		SeqInput reader(SeqIOOptions::genFastaIn(genomeFnp));
		reader.openIn();
		seqInfo seq;
		while (reader.readNextRead(seq)) {
			auto chromName = seq.name_.substr(0, seq.name_.find_first_of(" \t"));
			stringToUpper(seq.seq_);
			for (const auto &q: hitsByQuery) {
				const auto &needle = q.first;
				if (seq.seq_.size() < needle.size()) {
					continue;
				}
				auto found = seq.seq_.find(needle);
				while (std::string::npos != found) {
					truth[needle].emplace_back(chromName, found);
					found = seq.seq_.find(needle, found + 1);
				}
			}
		}
		uint32_t disagreements = 0;
		for (const auto &q: hitsByQuery) {
			std::set<std::pair<std::string, uint64_t>> fromIndex;
			for (const auto &hit: q.second) {
				fromIndex.emplace(hit.chrom_, hit.pos_);
			}
			std::set<std::pair<std::string, uint64_t>> fromGenome(truth[q.first].begin(), truth[q.first].end());
			if (fromIndex != fromGenome) {
				++disagreements;
				std::cerr << "MISMATCH for " << q.first << ": index found " << fromIndex.size()
						<< ", scanning the genome found " << fromGenome.size() << std::endl;
				for (const auto &g: fromGenome) {
					if (!njh::in(g, fromIndex)) {
						std::cerr << "\t" << "missed by the index: " << g.first << ":" << g.second << std::endl;
					}
				}
				for (const auto &i: fromIndex) {
					if (!njh::in(i, fromGenome)) {
						std::cerr << "\t" << "not actually there: " << i.first << ":" << i.second << std::endl;
					}
				}
			} else {
				std::cout << "#verified " << q.first << ": " << fromIndex.size()
						<< " hits agree with a direct scan of the genome" << std::endl;
			}
		}
		if (disagreements > 0) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << disagreements << " queries disagreed with the genome" << "\n";
			throw std::runtime_error{ss.str()};
		}
	}
	return 0;
}

}  // namespace njhseq
