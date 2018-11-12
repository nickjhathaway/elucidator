#pragma once
//
//  miscSetUp.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 11/3/13.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
//

#include <njhseq.h>

namespace njhseq {

class miscSetUp : public seqSetUp {

 public:
	using njhseq::seqSetUp::seqSetUp;
  void setUpListAllFiles(std::string& directory, VecStr& contains,
                         bool& specific, bool& recursive,
                         std::string& filesOrDirectories);
  void setUpProcessKrecFasta(std::string& filename, std::string& refFilename,
                             bool& postProcessed);
};
}  // namespace njhseq

