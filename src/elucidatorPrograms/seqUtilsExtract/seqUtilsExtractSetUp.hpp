#pragma once
//

//  seqUtilsExtractSetUp.hpp
//
//  Created by Nicholas Hathaway on 2015/05/28.
//  Copyright (c) 2015 Nicholas Hathaway. All rights reserved.
//

#include <njhseq.h>

namespace njhseq {


//struct ExtractByMIDPars {
//  bfs::path idFilename = "";
//  std::string idFileDelim = "tab";
//  uint32_t barcodeErrors = 0;
//  MidDeterminator::MidDeterminePars mDetPars;
//	uint32_t smallFragmentCutoff = 10;
//};

class seqUtilsExtractSetUp : public seqSetUp {

 public:
    using seqSetUp::seqSetUp;


//    void setUpExtractByMID(ExtractByMIDPars & pars);

    void setUpGetSimilarSequences(bool& useNucComp, double& maxNucCompDiff, double& idCutOff,
    		double& gapCutoff, double& queryCutoff);

    void setUpExtractSameSeqs(std::string& compareFileName);

};
} // namespace njhseq
