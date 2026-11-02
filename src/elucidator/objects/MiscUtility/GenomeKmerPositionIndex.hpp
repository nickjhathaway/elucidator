#pragma once

/*
 * GenomeKmerPositionIndex.hpp
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

#include <njhseq/objects/BioDataObject/GenomicRegion.hpp>


namespace njhseq {

/**@brief every position of every k-mer in a genome, on disk, so an exact-match query is a lookup
 *
 * A primer only primes where its 3'-terminal bases match exactly, so finding the places it could prime is a
 * k-mer lookup rather than an alignment search. Answering it with a short-read aligner means reporting
 * matches anywhere in the primer at a permissive e-value and then discarding most of them for not reaching
 * the 3' end; this answers it directly, in time proportional to the number of hits.
 *
 * Only the forward strand is indexed; query with the reverse complement to find minus-strand sites.
 * Any window containing a character other than A/C/G/T is skipped, so N runs cost nothing.
 *
 * The index is built at one k but answers any query from 1 to k, because a k-mer code places the FIRST base
 * in the most significant bits: every k-mer sharing a q-base prefix therefore occupies a contiguous block of
 * 4^(k-q) buckets, and since buckets are laid out in code order that is a single contiguous run of positions.
 * The position array holds one entry per genome position whatever k is, so building at a larger k costs only
 * a larger offset table and keeps the shorter query lengths available.
 *
 * File layout, all little-endian:
 *   header    magic[8] "NJHKIDX1", version, k, numKmers (4^k), numPositions, numChroms
 *   chroms    numChroms x { nameLen, name bytes, flat start offset, length }
 *   offsets   (numKmers + 1) x uint64, cumulative start of each k-mer's run in the position array
 *   positions numPositions x uint32, flat genome coordinates, ascending within each run
 */
class GenomeKmerPositionIndex {
public:
	static constexpr uint32_t VERSION = 1;
	static constexpr uint32_t MAX_K = 15;   /**< beyond this the offset table itself becomes unreasonable */

	struct ChromInfo {
		std::string name_;
		uint64_t flatStart_{0};  /**< where this chromosome begins in the flat coordinate space */
		uint64_t length_{0};
	};

	/**@brief a hit: which chromosome, where in it, and the flat coordinate it came from */
	struct Hit {
		std::string chrom_;
		uint64_t pos_{0};
		uint64_t flatPos_{0};
	};

	/**@brief open an existing index for querying; the position array is memory mapped, not read in
	 *
	 * Only the buckets actually queried are touched, so a panel's worth of lookups keeps a few hundred pages
	 * resident rather than the whole file.
	 */
	explicit GenomeKmerPositionIndex(const bfs::path &indexFnp);

	~GenomeKmerPositionIndex();

	GenomeKmerPositionIndex(const GenomeKmerPositionIndex &) = delete;
	GenomeKmerPositionIndex & operator=(const GenomeKmerPositionIndex &) = delete;

	/**@brief build an index for a genome fasta and write it to indexFnp
	 *
	 * One pass counts every k-mer, then the position array is filled a chunk at a time, each chunk being a
	 * contiguous range of k-mer codes that fits in buildBufferBytes. Every chunk is assembled in memory and
	 * written as one sequential block, so peak memory is the buffer plus the offset table and the cost does
	 * not depend on how much of the index the page cache can hold. Filling a chunk costs one sweep of the
	 * sequence, so a larger buffer means fewer sweeps.
	 */
	static void build(const bfs::path &genomeFnp, const bfs::path &indexFnp, uint32_t k,
										uint64_t buildBufferBytes = 4294967296, bool verbose = false);

	/**@brief where this genome's index lives, by the same convention the blast database uses */
	static bfs::path genIndexFnpForGenome(const bfs::path &genomeFnp, uint32_t k);

	/**@brief positions where query occurs exactly, as flat coordinates
	 *
	 * query must be at most k long and contain only A/C/G/T; anything else throws. The returned run is
	 * ascending.
	 */
	[[nodiscard]] std::vector<uint32_t> findFlatPositions(const std::string &query) const;

	/**@brief the same, resolved to chromosome and offset */
	[[nodiscard]] std::vector<Hit> find(const std::string &query) const;

	/**@brief how many times query occurs, without reading the positions */
	[[nodiscard]] uint64_t count(const std::string &query) const;

	[[nodiscard]] uint32_t k() const { return k_; }
	[[nodiscard]] uint64_t numPositions() const { return numPositions_; }
	[[nodiscard]] const std::vector<ChromInfo> & chroms() const { return chroms_; }

	/**@brief turn a flat coordinate back into a chromosome and an offset within it */
	[[nodiscard]] Hit resolve(uint32_t flatPos) const;

	/**@brief 2-bit code of a sequence, first base in the most significant bits
	 *
	 * Returns false when the sequence holds anything other than A/C/G/T, which is how N windows are skipped.
	 */
	static bool encode(const char *seq, uint32_t len, uint64_t &code);

private:
	bfs::path indexFnp_;
	uint32_t k_{0};
	uint64_t numKmers_{0};
	uint64_t numPositions_{0};
	std::vector<ChromInfo> chroms_;

	std::vector<uint64_t> offsets_;   /**< 4^k + 1 entries, read in full; the positions stay mapped */
	int fd_{-1};                      /**< the mapping of the position array */
	void *map_{nullptr};
	size_t mapLen_{0};
	const uint32_t *positions_{nullptr};

	/**@brief the half-open bucket range a query of length <= k covers */
	void queryRange(const std::string &query, uint64_t &lowKmer, uint64_t &highKmer) const;
};

}  // namespace njhseq
