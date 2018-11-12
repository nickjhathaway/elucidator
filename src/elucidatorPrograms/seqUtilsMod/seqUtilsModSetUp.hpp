#pragma once
//

//  seqUtilsModSetUp.hpp
//
//  Created by Nicholas Hathaway on 2015/05/28.
//  Copyright (c) 2015 Nicholas Hathaway. All rights reserved.
//

#include <njhseq.h>

namespace njhseq {



class  seqUtilsModSetUp: public seqSetUp {

 public:
    using seqSetUp::seqSetUp;

    void setUpSortReads(std::string& sortBy, bool& decending);
    void setUpRenameIDs(std::string& stub, std::string& sortBy,
                                       bool& keepChimeraFlag);

    void setUpComplementSeq(std::string &seqType);

    void setUpSplit(std::string& splitOption, uint32_t& minLen,
    		uint32_t& within, std::string& runCutoffString,
                                   std::string& nameContains,
                                   std::string& seqContains, uint32_t& occurences,
																	 uint32_t& maxLength, uint32_t& qualWindowSize,
																	 uint32_t& qualWindowStep, uint32_t& qualWindowThres);
	void setUpTranslate(uint64_t &start, bool &complement, bool &reverse);
};
} // namespace njhseq
