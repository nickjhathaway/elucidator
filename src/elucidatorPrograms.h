#pragma once

// Created on 2014/09/11
// Including headers in common
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
#include "elucidatorPrograms/ampliconAnalysis.h"
#include "elucidatorPrograms/bamCountingExp.hpp"
#include "elucidatorPrograms/bamExp.hpp"
#include "elucidatorPrograms/bedExp.hpp"
#include "elucidatorPrograms/clusteringExp.hpp"
#include "elucidatorPrograms/fileFormatExp.hpp"
#include "elucidatorPrograms/genExp.hpp"
#include "elucidatorPrograms/geneExp.hpp"
#include "elucidatorPrograms/gffExp.hpp"
#include "elucidatorPrograms/graphicsUtils.h"
#include "elucidatorPrograms/kmerExp.hpp"
#include "elucidatorPrograms/jsonExp.hpp"
#include "elucidatorPrograms/misc.h"
#include "elucidatorPrograms/pacbioExp.hpp"
#include "elucidatorPrograms/printInfo.h"
#include "elucidatorPrograms/programWrappers.hpp"
#include "elucidatorPrograms/readSimulator.h"
#include "elucidatorPrograms/repelin.h"
#include "elucidatorPrograms/seqUtils.h"
#include "elucidatorPrograms/seqUtilsConverter.h"
#include "elucidatorPrograms/seqUtilsExtract.h"
#include "elucidatorPrograms/seqUtilsInfo.h"
#include "elucidatorPrograms/seqUtilsMod.h"
#include "elucidatorPrograms/seqUtilsTrim.h"
#include "elucidatorPrograms/seqUtilsSplit.h"
#include "elucidatorPrograms/simulator.h"
#include "elucidatorPrograms/parsingFileExp.hpp"
#include "elucidatorPrograms/benchMarking.hpp"
#include "elucidatorPrograms/metaExp.hpp"
#include "elucidatorPrograms/pairProcessing.hpp"

#include "elucidatorPrograms/seqSearching.hpp"


#include <njhcpp/progutils/oneRing.hpp>

#include <seqServerPrograms.h>
#include <SeekDeepPrograms.h>
#include <TwoBitPrograms.h>
#include <njhseq/ProgramRunners.h>

namespace njhseq {



class elucidatorRunner: public njh::progutils::OneRing {
public:
	elucidatorRunner();
};
//,createReadConsensusRandomPicking  BamExtractPathawaysFromRegion
elucidatorRunner::elucidatorRunner() :
		njh::progutils::OneRing(
				{ addRing<pacbioExpRunner>(),      addRing<genExpRunner>(),
					addRing<kmerExpRunner>(),        addRing<clusteringExpRunner>(),
					addRing<fileFormatExpRunner>(),
					addRing<bamExpRunner>(),         addRing<bamCountingExpRunner>(),
					addRing<programWrapperRunner>(),
					addRing<seqUtilsTrimRunner>(),
					addRing<seqUtilsModRunner>(),
					addRing<gffExpRunner>(),   			 addRing<bedExpRunner>(),
					addRing<geneExpRunner>(),        addRing<SeqServerRunner>(),
					addRing<SeekDeepRunner>(),       addRing<TwoBit::TwoBitRunner>(),
					addRing<repelinRunner>(),        addRing<ampliconAnalysisRunner>(),
					addRing<readSimulatorRunner>(),  addRing<graphicsUtilsRunner>(),
					addRing<miscRunner>(),           addRing<printInfoRunner>(),
					addRing<seqUtilsRunner>(),       addRing<seqUtilsConverterRunner>(),
					addRing<seqUtilsExtractRunner>(),addRing<seqUtilsInfoRunner>(),
					addRing<simulatorRunner>(),      addRing<parsingFileExpRunner>(),
					addRing<benchMarkingRunner>(),   addRing<jsonExpRunner>(),
					addRing<ManipulateTableRunner>(),addRing<seqUtilsSplitRunner>(),
					addRing<metaExpRunner>(),
					addRing<seqSearchingRunner>(),
					addRing<pairProcessingRunner>(),
				},//
				{}, "elucidator", "1", "1", "0") {
}

}  // namespace njhseq
