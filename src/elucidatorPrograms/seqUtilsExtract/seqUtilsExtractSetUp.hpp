#pragma once
//

//  seqUtilsExtractSetUp.hpp
//
//  Created by Nicholas Hathaway on 2015/05/28.
// elucidator - A library for analyzing sequence data
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
//
// This file is part of elucidator.
//
// elucidator is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// elucidator is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with elucidator.  If not, see <http://www.gnu.org/licenses/>.
//

#include <njhseq/common.h>
#include <njhseq/programUtils/seqSetUp.hpp>

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



    void setUpExtractSameSeqs(std::string& compareFileName);

};
} // namespace njhseq
