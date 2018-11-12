#pragma once
//

//  seqUtilsTrimSetUp.hpp
//
//
//  Created by Nicholas Hathaway on 2016/10/05.
//  Copyright (c) 2016 Nicholas Hathaway. All rights reserved.
//
//

#include <njhseq.h>

namespace njhseq {



class seqUtilsTrimSetUp: public seqSetUp {

public:
	using seqSetUp::seqSetUp;

	void setUpTrimFront(FullTrimReadsPars & pars);
	void setUpTrimAtFirstQual(FullTrimReadsPars & pars);
	void setUpTrimAtFirstBase(FullTrimReadsPars & pars);
	void setUpTrimAtLastBase(FullTrimReadsPars & pars);
	void setUpTrimEnd(FullTrimReadsPars & pars);
	void setUpTrimEnds(FullTrimReadsPars & pars);
	void setUpTrimBeforeSeq(FullTrimReadsPars & pars);
	void setUpTrimFromSeq(FullTrimReadsPars & pars);
	void setUpTrimBetweenSeqs(FullTrimReadsPars & pars);
	void setUpTrimToLen(FullTrimReadsPars & pars);
	void setUpTrimToSimilarSeq(FullTrimReadsPars & pars);


	void processIoOptions(const VecStr & formats);
	void processIoOptions();

	void procesingTrimmingWithSeqsOpts(FullTrimReadsPars & pars);


};




} // namespace njhseq
