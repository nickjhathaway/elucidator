//
// Created by Nicholas Hathaway on 5/25/23.
//

#include "UniqueKmerSetHelper.hpp"


namespace njhseq {
uint32_t UniqueKmerSetHelper::getKmerLenFromUniqueKmerTable(const bfs::path &uniqueKmerTable) {
	auto toks = tokenizeString(njh::files::getFirstLine(uniqueKmerTable), "\t");
	return toks.at(1).size();
}


std::unordered_map<std::string, std::unordered_set<uint64_t>> UniqueKmerSetHelper::readInUniqueKmerTablePerSet(const bfs::path &uniqueKmerTableFnp) {
	std::unordered_map<std::string, std::unordered_set<uint64_t>> ret;
	SimpleKmerHash hasher;
	TableReader uniqKmers(TableIOOpts::genTabFileIn(uniqueKmerTableFnp, false));
	if (uniqKmers.header_.nCol() < 2) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "need to have 2 columns" << "\n";
		throw std::runtime_error{ss.str()};
	}
	VecStr row;
	while (uniqKmers.getNextRow(row)) {
		ret[row[0]].emplace(hasher.hash(row[1]));
	}
	return ret;
}

std::unordered_set<uint64_t> UniqueKmerSetHelper::readInUniqueKmerTableSetsCollapsed(const bfs::path &uniqueKmerTableFnp) {
	std::unordered_set<uint64_t> ret;
	SimpleKmerHash hasher;
	TableReader uniqKmers(TableIOOpts::genTabFileIn(uniqueKmerTableFnp, false));
	if (uniqKmers.header_.nCol() < 2) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "need to have 2 columns" << "\n";
		throw std::runtime_error{ss.str()};
	}
	VecStr row;
	while (uniqKmers.getNextRow(row)) {
		ret.emplace(hasher.hash(row[1]));
	}
	return ret;
}

}  // namespace njhseq

