#pragma once
//
//  seqUtilsSetUp.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 10/17/13.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
//

#include <njhseq.h>

namespace njhseq {

class seqUtilsSetUp : public seqSetUp {

 public:
	using seqSetUp::seqSetUp;


  void setUpMapToReferenceCount(bool& extra);
  void setUpCompareToRef();

};
}  // namespace njhseq
