#pragma once
//

//  seqUtilsModSetUp.hpp
//
//  Created by Nicholas Hathaway on 2015/05/28.
//  Copyright (c) 2015 Nicholas Hathaway. All rights reserved.
//

#include <njhseq.h>

namespace njhseq {

class seqUtilsModSetUp: public seqSetUp {

public:
	using seqSetUp::seqSetUp;

	void setUpSortReads(std::string& sortBy, bool& decending);
	void setUpRenameIDs(std::string& stub, std::string& sortBy,
			bool& keepChimeraFlag);

	void setUpComplementSeq(std::string &seqType);

	void setUpTranslate(uint64_t &start, bool &complement, bool &reverse);
};
} // namespace njhseq
