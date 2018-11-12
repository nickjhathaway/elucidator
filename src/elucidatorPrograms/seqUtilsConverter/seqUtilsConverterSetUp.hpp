#pragma once
//

//  seqUtilsConverterSetUp.hpp
//
//  Created by Nicholas Hathaway on 2015/05/28.
//  Copyright (c) 2015 Nicholas Hathaway. All rights reserved.
//

#include <njhseq.h>

namespace njhseq {

class seqUtilsConverterSetUp : public seqSetUp {

 public:
    using seqSetUp::seqSetUp;

    void setUpConvertFiles();
};
} // namespace njhseq
