/*
 * GenomeKmerPositionIndex.cpp
 *
 *  Created on: Nov 2, 2026
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

#include "GenomeKmerPositionIndex.hpp"

#include <njhseq/IO/SeqIO.h>

#include <functional>

#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>


namespace njhseq {

namespace {

const char INDEX_MAGIC[8] = {'N', 'J', 'H', 'K', 'I', 'D', 'X', '1'};

//-1 for anything that is not A/C/G/T, which is what makes N windows skippable
constexpr int8_t baseCode(char c) {
	switch (c) {
		case 'A': case 'a': return 0;
		case 'C': case 'c': return 1;
		case 'G': case 'g': return 2;
		case 'T': case 't': return 3;
		default: return -1;
	}
}

template<typename T>
void writeRaw(std::ostream &out, const T &val) {
	out.write(reinterpret_cast<const char *>(&val), sizeof(T));
}

template<typename T>
void readRaw(std::istream &in, T &val) {
	in.read(reinterpret_cast<char *>(&val), sizeof(T));
}

}  // namespace


bool GenomeKmerPositionIndex::encode(const char *seq, uint32_t len, uint64_t &code) {
	code = 0;
	for (uint32_t pos = 0; pos < len; ++pos) {
		const auto base = baseCode(seq[pos]);
		if (base < 0) {
			return false;
		}
		//first base ends up in the most significant bits, which is what makes a prefix a contiguous range
		code = (code << 2) | static_cast<uint64_t>(base);
	}
	return true;
}

bfs::path GenomeKmerPositionIndex::genIndexFnpForGenome(const bfs::path &genomeFnp, uint32_t k) {
	//matches how the blast database is derived from the genome path
	return bfs::path(njh::rstripRet(bfs::path(genomeFnp).replace_extension("").string(), '.')
									 + "_" + estd::to_string(k) + "merPosIndex.bin");
}

void GenomeKmerPositionIndex::build(const bfs::path &genomeFnp, const bfs::path &indexFnp, uint32_t k,
																		uint64_t buildBufferBytes, bool verbose) {
	if (0 == k || k > MAX_K) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error k must be between 1 and " << MAX_K << ", got " << k << "\n";
		throw std::runtime_error{ss.str()};
	}
	if (0 == buildBufferBytes) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error the build buffer must be larger than zero" << "\n";
		throw std::runtime_error{ss.str()};
	}
	const uint64_t numKmers = static_cast<uint64_t>(1) << (2 * k);
	const uint64_t codeMask = (64 == 2 * k) ? std::numeric_limits<uint64_t>::max() : (numKmers - 1);

	//Walk a sequence's k-mers with a rolling code rather than re-reading k bases at every position. A base
	//that is not A/C/G/T breaks the run, so the next k-1 positions are skipped and N stretches cost nothing.
	auto forEachKmer = [&k, &codeMask](const std::string &seq, const std::function<void(size_t, uint64_t)> &handle) {
		uint64_t code = 0;
		uint32_t runLen = 0;
		for (size_t pos = 0; pos < seq.size(); ++pos) {
			const auto base = baseCode(seq[pos]);
			if (base < 0) {
				runLen = 0;
				continue;
			}
			code = ((code << 2) | static_cast<uint64_t>(base)) & codeMask;
			++runLen;
			if (runLen >= k) {
				handle(pos - k + 1, code);
			}
		}
	};

	//PASS 1 - count each k-mer and record where each chromosome sits in the flat coordinate space. The
	//sequence is read a record at a time, so the peak here is one chromosome plus the count table.
	std::vector<ChromInfo> chroms;
	std::vector<uint64_t> offsets(numKmers + 1, 0);
	uint64_t flatLen = 0;
	uint64_t numPositions = 0;
	{
		SeqInput reader(SeqIOOptions::genFastaIn(genomeFnp));
		reader.openIn();
		seqInfo seq;
		while (reader.readNextRead(seq)) {
			ChromInfo info;
			//the fasta header can carry a description; the identifier is the first token, as elsewhere
			info.name_ = seq.name_.substr(0, seq.name_.find_first_of(" \t"));
			info.flatStart_ = flatLen;
			info.length_ = seq.seq_.size();
			chroms.emplace_back(info);
			flatLen += seq.seq_.size();
			forEachKmer(seq.seq_, [&offsets, &numPositions](size_t, uint64_t code) {
				++offsets[code];
				++numPositions;
			});
			if (verbose) {
				std::cout << "counted " << info.name_ << " (" << info.length_ << " bp), "
						<< numPositions << " indexable positions so far" << std::endl;
			}
		}
	}
	if (flatLen > std::numeric_limits<uint32_t>::max()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error genome is " << flatLen
				<< " bp, which does not fit a 32 bit position; this index holds genomes up to "
				<< std::numeric_limits<uint32_t>::max() << " bp" << "\n";
		throw std::runtime_error{ss.str()};
	}

	//cumulative offsets, in place so a second table of 4^k is not needed
	uint64_t running = 0;
	for (uint64_t kmer = 0; kmer < numKmers; ++kmer) {
		const auto thisCount = offsets[kmer];
		offsets[kmer] = running;
		running += thisCount;
	}
	offsets[numKmers] = running;

	//Split the position array into chunks that each fit in the build buffer. Buckets are laid out in code
	//order, so a chunk is a contiguous range of codes AND a contiguous range of the output file. Each chunk
	//is filled in memory and written in one sequential block, which keeps the build's cost independent of how
	//much of the index would fit in the page cache. A single bucket larger than the buffer gets a chunk of
	//its own and is the one case that can exceed it.
	const uint64_t perChunk = std::max<uint64_t>(1, buildBufferBytes / sizeof(uint32_t));
	std::vector<std::pair<uint64_t, uint64_t>> chunks;  //half open code ranges
	{
		uint64_t chunkStartCode = 0;
		for (uint64_t kmer = 0; kmer < numKmers; ++kmer) {
			if (offsets[kmer + 1] - offsets[chunkStartCode] >= perChunk) {
				chunks.emplace_back(chunkStartCode, kmer + 1);
				chunkStartCode = kmer + 1;
			}
		}
		if (chunkStartCode < numKmers) {
			chunks.emplace_back(chunkStartCode, numKmers);
		}
	}
	if (verbose) {
		std::cout << numPositions << " positions over " << chunks.size()
				<< " sequential write chunks of at most " << (buildBufferBytes / 1048576.0) << " MB" << std::endl;
	}

	std::ofstream out(indexFnp.string(), std::ios::binary | std::ios::trunc);
	if (!out) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error could not open " << indexFnp << " for writing" << "\n";
		throw std::runtime_error{ss.str()};
	}
	out.write(INDEX_MAGIC, sizeof(INDEX_MAGIC));
	writeRaw(out, VERSION);
	writeRaw(out, k);
	writeRaw(out, numKmers);
	writeRaw(out, numPositions);
	const auto numChroms = static_cast<uint32_t>(chroms.size());
	writeRaw(out, numChroms);
	for (const auto &chrom: chroms) {
		const auto nameLen = static_cast<uint32_t>(chrom.name_.size());
		writeRaw(out, nameLen);
		out.write(chrom.name_.data(), nameLen);
		writeRaw(out, chrom.flatStart_);
		writeRaw(out, chrom.length_);
	}
	out.write(reinterpret_cast<const char *>(offsets.data()),
						static_cast<std::streamsize>((numKmers + 1) * sizeof(uint64_t)));

	//PASS 2 - one sweep of the sequence per chunk, placing only the positions whose k-mer falls in that
	//chunk. Appending in sequence order leaves every bucket's run ascending, which is what lets a query hand
	//back a sorted range.
	std::vector<uint32_t> buffer;
	std::vector<uint64_t> cursors;
	for (const auto &chunk: iter::enumerate(chunks)) {
		const auto codeLo = chunk.element.first;
		const auto codeHi = chunk.element.second;
		const auto posLo = offsets[codeLo];
		const auto posHi = offsets[codeHi];
		buffer.assign(posHi - posLo, 0);
		cursors.assign(offsets.begin() + static_cast<int64_t>(codeLo),
									 offsets.begin() + static_cast<int64_t>(codeHi));
		SeqInput reader(SeqIOOptions::genFastaIn(genomeFnp));
		reader.openIn();
		seqInfo seq;
		uint32_t chromPos = 0;
		while (reader.readNextRead(seq)) {
			const auto flatStart = chroms[chromPos].flatStart_;
			forEachKmer(seq.seq_, [&](size_t pos, uint64_t code) {
				if (code >= codeLo && code < codeHi) {
					buffer[cursors[code - codeLo] - posLo] = static_cast<uint32_t>(flatStart + pos);
					++cursors[code - codeLo];
				}
			});
			++chromPos;
		}
		out.write(reinterpret_cast<const char *>(buffer.data()),
							static_cast<std::streamsize>(buffer.size() * sizeof(uint32_t)));
		if (verbose) {
			std::cout << "wrote chunk " << (chunk.index + 1) << " of " << chunks.size() << ", "
					<< buffer.size() << " positions" << std::endl;
		}
	}
	out.close();
	if (!out) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error writing " << indexFnp << "\n";
		throw std::runtime_error{ss.str()};
	}
	if (verbose) {
		std::cout << "wrote " << indexFnp << ": k=" << k << ", " << numPositions << " positions" << std::endl;
	}
}

GenomeKmerPositionIndex::GenomeKmerPositionIndex(const bfs::path &indexFnp) : indexFnp_(indexFnp) {
	uint64_t positionsStart = 0;
	{
		std::ifstream in(indexFnp.string(), std::ios::binary);
		if (!in) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error could not open " << indexFnp << "\n";
			throw std::runtime_error{ss.str()};
		}
		char magic[sizeof(INDEX_MAGIC)];
		in.read(magic, sizeof(magic));
		if (0 != std::memcmp(magic, INDEX_MAGIC, sizeof(magic))) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << indexFnp << " is not a k-mer position index" << "\n";
			throw std::runtime_error{ss.str()};
		}
		uint32_t version = 0;
		readRaw(in, version);
		if (VERSION != version) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << indexFnp << " is version " << version
					<< ", this build reads version " << VERSION << "\n";
			throw std::runtime_error{ss.str()};
		}
		readRaw(in, k_);
		readRaw(in, numKmers_);
		readRaw(in, numPositions_);
		uint32_t numChroms = 0;
		readRaw(in, numChroms);
		for (uint32_t chromPos = 0; chromPos < numChroms; ++chromPos) {
			ChromInfo info;
			uint32_t nameLen = 0;
			readRaw(in, nameLen);
			info.name_.resize(nameLen);
			in.read(info.name_.data(), nameLen);
			readRaw(in, info.flatStart_);
			readRaw(in, info.length_);
			chroms_.emplace_back(info);
		}
		offsets_.resize(numKmers_ + 1);
		in.read(reinterpret_cast<char *>(offsets_.data()),
						static_cast<std::streamsize>((numKmers_ + 1) * sizeof(uint64_t)));
		if (!in) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << indexFnp << " ended before its offset table did" << "\n";
			throw std::runtime_error{ss.str()};
		}
		positionsStart = static_cast<uint64_t>(in.tellg());
	}

	fd_ = ::open(indexFnp.string().c_str(), O_RDONLY);
	if (fd_ < 0) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error could not open " << indexFnp << " for mapping" << "\n";
		throw std::runtime_error{ss.str()};
	}
	const auto pageSize = static_cast<uint64_t>(::sysconf(_SC_PAGESIZE));
	const uint64_t mapStart = (positionsStart / pageSize) * pageSize;
	mapLen_ = positionsStart + numPositions_ * sizeof(uint32_t) - mapStart;
	map_ = ::mmap(nullptr, mapLen_, PROT_READ, MAP_SHARED, fd_, static_cast<off_t>(mapStart));
	if (MAP_FAILED == map_) {
		::close(fd_);
		fd_ = -1;
		map_ = nullptr;
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error could not map " << indexFnp << "\n";
		throw std::runtime_error{ss.str()};
	}
	positions_ = reinterpret_cast<const uint32_t *>(static_cast<char *>(map_) + (positionsStart - mapStart));
}

GenomeKmerPositionIndex::~GenomeKmerPositionIndex() {
	if (nullptr != map_) {
		::munmap(map_, mapLen_);
	}
	if (fd_ >= 0) {
		::close(fd_);
	}
}

void GenomeKmerPositionIndex::queryRange(const std::string &query, uint64_t &lowKmer, uint64_t &highKmer) const {
	if (query.empty() || query.size() > k_) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error query must be between 1 and " << k_ << " long, got "
				<< query.size() << " (\"" << query << "\")" << "\n";
		throw std::runtime_error{ss.str()};
	}
	uint64_t prefix = 0;
	if (!encode(query.data(), static_cast<uint32_t>(query.size()), prefix)) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error query \"" << query << "\" holds a base other than A/C/G/T" << "\n";
		throw std::runtime_error{ss.str()};
	}
	//a q-base prefix covers every k-mer that starts with it, which is 4^(k-q) codes laid out consecutively
	const auto shift = 2 * (k_ - static_cast<uint32_t>(query.size()));
	lowKmer = prefix << shift;
	highKmer = lowKmer + (static_cast<uint64_t>(1) << shift);  //half open
}

uint64_t GenomeKmerPositionIndex::count(const std::string &query) const {
	uint64_t lowKmer = 0;
	uint64_t highKmer = 0;
	queryRange(query, lowKmer, highKmer);
	return offsets_[highKmer] - offsets_[lowKmer];
}

std::vector<uint32_t> GenomeKmerPositionIndex::findFlatPositions(const std::string &query) const {
	uint64_t lowKmer = 0;
	uint64_t highKmer = 0;
	queryRange(query, lowKmer, highKmer);
	const auto start = offsets_[lowKmer];
	const auto stop = offsets_[highKmer];
	std::vector<uint32_t> ret(positions_ + start, positions_ + stop);
	//within one bucket the run is already ascending, but a query shorter than k spans several buckets
	if (highKmer - lowKmer > 1) {
		njh::sort(ret);
	}
	return ret;
}

GenomeKmerPositionIndex::Hit GenomeKmerPositionIndex::resolve(uint32_t flatPos) const {
	//chroms_ is in flat order, so the owning chromosome is the last one starting at or before flatPos
	auto search = std::upper_bound(chroms_.begin(), chroms_.end(), static_cast<uint64_t>(flatPos),
																 [](uint64_t pos, const ChromInfo &info) {
																	 return pos < info.flatStart_;
																 });
	if (chroms_.begin() == search) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error flat position " << flatPos << " is before the first chromosome"
				<< "\n";
		throw std::runtime_error{ss.str()};
	}
	--search;
	Hit hit;
	hit.chrom_ = search->name_;
	hit.flatPos_ = flatPos;
	hit.pos_ = flatPos - search->flatStart_;
	return hit;
}

std::vector<GenomeKmerPositionIndex::Hit> GenomeKmerPositionIndex::find(const std::string &query) const {
	std::vector<Hit> ret;
	for (const auto flatPos: findFlatPositions(query)) {
		ret.emplace_back(resolve(flatPos));
	}
	return ret;
}

}  // namespace njhseq
