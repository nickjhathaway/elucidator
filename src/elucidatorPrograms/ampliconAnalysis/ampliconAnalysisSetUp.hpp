#pragma once
//
//  ampliconAnalysisSetUp.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 10/24/13.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
//

#include <njhseq.h>


namespace njhseq {

class ampliconAnalysisSetUp : public seqSetUp {

 public:
  // constructors
	using njhseq::seqSetUp::seqSetUp;


  void setUpCollapseTandems(double& freqCutoff, bool& extra,
                            bool& additionalOut,
                            std::string& additionalOutLocationFile);
  void setUpMarkChimeras();

};
}  // namespace njhseq

