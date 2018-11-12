#pragma once
/*
 * SlimCounterMaster.hpp
 *
 *  Created on: Jun 16, 2016
 *      Author: nick
 */



#include "elucidator/common.h"
#include "elucidator/objects/SlimCounter/SlimCounterRef.hpp"

namespace njhseq {
class SlimCounterMaster {
public:
	SlimCounterMaster(bfs::path outputDir,
			const std::vector<BamTools::RefData> & refInfos);

	const bfs::path mainOutputDir_;
	const bfs::path countsOutputDir_;

	std::vector<BamTools::RefData> refInfos_;

	std::unordered_map<uint32_t, std::shared_ptr<SlimCounterRef>> counts_;

	std::shared_ptr<SlimCounterRef> getSlimCounterRef(uint32_t refId);

	std::shared_ptr<SlimCounterRef> getSlimCounterRef(
			const std::string & refName);

	void increaseCount(uint32_t refId, uint32_t pos, char base, bool reverseStand,
			bool highQuality, uint32_t count);

	std::mutex mut_;


	static void checkDirectoryStructure(const bfs::path & dir);
	static std::unordered_map<uint64_t, uint64_t> readGzChromIndex(const bfs::path & filename);

};

}  // namespace njhseq





