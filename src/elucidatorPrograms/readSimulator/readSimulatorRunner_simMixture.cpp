/*
 * readSimulatorRunner_simMixture.cpp
 *
 *  Created on: Sep 16, 2018
 *      Author: nick
 */
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
#include "readSimulatorRunner.hpp"
#include "elucidator/simulation.h"
#include <SeekDeep/objects/TarAmpSetupUtils/PrimersAndMids.hpp>
namespace njhseq {


int readSimulatorRunner::createIlluminaErrorProfile(
		const njh::progutils::CmdArgs & inputCommands) {
	bfs::path bedFnp = "";
	bfs::path twoBitFnp = "";
	uint32_t numThreads = 1;
	uint32_t lengthLimit = 0;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processReadInNames(VecStr{"--bam"});
	setUp.processDirectoryOutputName(true);
	setUp.setOption(lengthLimit, "--lengthLimit", "Length Limit");
	setUp.setOption(numThreads, "--numThreads", "Number Threads");
	setUp.setOption(bedFnp, "--bed", "Bed file");
	setUp.setOption(twoBitFnp, "--twoBitFnp", "2bit file for genome aligned to", true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	BamTools::BamReader bReader;
	bReader.Open(setUp.pars_.ioOptions_.firstName_.string());
	checkBamOpenThrow(bReader, setUp.pars_.ioOptions_.firstName_);
	loadBamIndexThrow(bReader);
	auto refs = bReader.GetReferenceData();
	BamAlnsCache alnCache;
	auto refData = bReader.GetReferenceData();
	std::unordered_map<std::string, uint32_t> refNameToId;
	for (auto pos : iter::range(refData.size())) {
		refNameToId[refData[pos].RefName] = pos;
	}

	TwoBit::TwoBitFile tReader(twoBitFnp);
	RoughIlluminaProfiler profiler;
	if("" != bedFnp){
		auto regions = bedPtrsToGenomicRegs(getBeds(bedFnp));
		njh::concurrent::LockableQueue<GenomicRegion> regionsQueue(regions);
		concurrent::BamReaderPool bamPools(setUp.pars_.ioOptions_.firstName_, numThreads);
		bamPools.openBamFile();
		std::mutex mut;
		auto increaseCount = [&regionsQueue,&bamPools,&mut,&profiler,&twoBitFnp,&refData,&setUp,&lengthLimit](){
			GenomicRegion region;
			TwoBit::TwoBitFile tReader(twoBitFnp);
			auto bReader = bamPools.popReader();
			BamTools::BamAlignment bAln;
			RoughIlluminaProfiler currentProfiler;
			while(regionsQueue.getVal(region)){
				if(setUp.pars_.verbose_){
					std::lock_guard<std::mutex> lock(mut);
					std::cout << region.uid_ << std::endl;
				}
				setBamFileRegionThrow(*bReader, region);
				while (bReader->GetNextAlignment(bAln)) {
					if(bAln.IsMapped() && bAln.IsPrimaryAlignment()){
						if(bAln.GetEndPosition() - bAln.Position > lengthLimit){
							currentProfiler.increaseCounts(bAln, refData, tReader);
						}
					}
				}
			}
			{
				std::lock_guard<std::mutex> lock(mut);
				profiler.addOther(currentProfiler);
			}
		};
		std::vector<std::thread> threads;
		for(uint32_t t = 0; t < numThreads; ++t){
			threads.emplace_back(std::thread(increaseCount));
		}
		njh::concurrent::joinAllJoinableThreads(threads);
	}else{
		BamTools::BamAlignment bAln;
		while (bReader.GetNextAlignment(bAln)) {
			if(bAln.IsMapped() && bAln.IsPrimaryAlignment()){
				if(bAln.GetEndPosition() - bAln.Position > lengthLimit){
					profiler.increaseCounts(bAln, refData, tReader);
				}
			}
		}
	}
	profiler.r1_counts.writeProfiles(njh::files::make_path(setUp.pars_.directoryName_, "r1").string(), false);
	profiler.r2_counts.writeProfiles(njh::files::make_path(setUp.pars_.directoryName_, "r2").string(), false);
	return 0;
}



int createLibrarySimMultipleMixtureExampleTesting(
		const njh::progutils::CmdArgs & inputCommands) {
	LibrarySetup::SimLibrarySetupPars simPars;
	OutOptions outOpts(bfs::path(""), ".json");
	readSimulatorSetUp setUp(inputCommands);
	setUp.processWritingOptions(outOpts);
	setUp.setOption(simPars.barcodeRandomPrecedingBases_, "--barcodeRandomPrecedingBases", "Barcode Random Preceding Bases");
	setUp.setOption(simPars.addReverseComplement_, "--addReverseComplement", "Add Reverse Complement");
	setUp.setOption(simPars.addBluntEndingArtifact_, "--addBluntEndingArtifact", "Add Blunt Ending Artifact");
	setUp.setOption(simPars.bluntEndingArtifactChance_, "--bluntEndingArtifactChance", "Blunt Ending Artifact Chance");
	setUp.finishSetUp(std::cout);
	OutputStream out(outOpts);
	LibrarySetup lSetup("testing_barcode_scheme", simPars);
	{
		auto sampleSet = std::make_shared<SampleSetup>("Sample1-FrontBarcode");
		{
			auto mixture = std::make_shared<MixtureSetUp>("Mixture1");
			mixture->setPrimers("ama1_S0_Sub0_mip2", "CAATATTTAAAAGATGGAGGTTTTGCT",
					"ATCTTCACAATTTCCATCGACCCAT");
			mixture->addAbundance("PfML01", 10);
			mixture->addAbundance("PfDd2", 10);
			mixture->addAbundance("PfGB4-PfTG01", 10);
			mixture->addAbundance("PfKH02", 10);
			mixture->setForwardBarcode("MID07", "CGTGTCTCTA");
			mixture->forwardBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
			sampleSet->addMixture(mixture);
		}
		{
			auto mixture = std::make_shared<MixtureSetUp>("Mixture2");
			mixture->setPrimers("ama1_S0_Sub0_mip2", "CAATATTTAAAAGATGGAGGTTTTGCT",
					"ATCTTCACAATTTCCATCGACCCAT");
			mixture->addAbundance("PfML01", 10);
			mixture->addAbundance("PfDd2", 20);
			mixture->addAbundance("PfGB4-PfTG01", 30);
			mixture->addAbundance("PfKH02", 40);
			mixture->setForwardBarcode("MID08", "CTCGCGTGTC");
			mixture->forwardBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
			sampleSet->addMixture(mixture);
		}

		{
			auto mixture = std::make_shared<MixtureSetUp>("Mixture3");
			mixture->setPrimers("ama1_S0_Sub0_mip2", "CAATATTTAAAAGATGGAGGTTTTGCT",
					"ATCTTCACAATTTCCATCGACCCAT");
			mixture->addAbundance("PfML01", 10);
			mixture->addAbundance("PfDd2", 20);
			mixture->addAbundance("PfGB4-PfTG01", 30);
			mixture->addAbundance("PfKH02", 40);
			mixture->setForwardBarcode("MID10", "TCTCTATGCG");
			mixture->forwardBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
			sampleSet->addMixture(mixture);
		}
		{
			auto mixture = std::make_shared<MixtureSetUp>("Mixture4");
			mixture->setPrimers("ama1_S0_Sub0_mip2", "CAATATTTAAAAGATGGAGGTTTTGCT",
					"ATCTTCACAATTTCCATCGACCCAT");
			mixture->addAbundance("PfML01", 10);
			mixture->addAbundance("PfDd2", 20);
			mixture->addAbundance("PfGB4-PfTG01", 30);
			mixture->addAbundance("PfKH02", 40);
			mixture->setForwardBarcode("MID11", "TGATACGTCT");
			mixture->forwardBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
			sampleSet->addMixture(mixture);
		}
		{
			auto mixture = std::make_shared<MixtureSetUp>("Mixture5");
			mixture->setPrimers("ama1_S0_Sub0_mip2", "CAATATTTAAAAGATGGAGGTTTTGCT",
					"ATCTTCACAATTTCCATCGACCCAT");
			mixture->addAbundance("PfML01", 2);
			mixture->addAbundance("PfDd2", 3);
			mixture->addAbundance("PfGB4-PfTG01", 5);
			mixture->addAbundance("PfKH02", 90);
			mixture->setForwardBarcode("MID01", "ACGAGTGCGT");
			mixture->forwardBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
			sampleSet->addMixture(mixture);
		}
		{
			auto mixture = std::make_shared<MixtureSetUp>("Mixture6");
			mixture->setPrimers("ama1_S0_Sub0_mip2", "CAATATTTAAAAGATGGAGGTTTTGCT",
					"ATCTTCACAATTTCCATCGACCCAT");
			mixture->addAbundance("PfML01", 1);
			mixture->addAbundance("PfDd2", 1);
			mixture->addAbundance("PfGB4-PfTG01", 1);
			mixture->addAbundance("PfKH02", 99);
			mixture->setForwardBarcode("MID29", "ACTGTACAGT");
			mixture->forwardBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
			sampleSet->addMixture(mixture);
		}
		lSetup.addSample(sampleSet);
	}
//	{
//		auto sampleSet = std::make_shared<SampleSetup>("Sample2-ReverseBarcode");
//		{
//			auto mixture = std::make_shared<MixtureSetUp>("Mixture1");
//			mixture->setPrimers("ama1_S0_Sub0_mip2", "CAATATTTAAAAGATGGAGGTTTTGCT",
//					"ATCTTCACAATTTCCATCGACCCAT");
//			mixture->addAbundance("PfML01", 10);
//			mixture->addAbundance("PfDd2", 10);
//			mixture->addAbundance("PfGB4-PfTG01", 10);
//			mixture->addAbundance("PfKH02", 10);
//
//			mixture->setReverseBarcode("MID07", "CGTGTCTCTA");
//			sampleSet->addMixture(mixture);
//		}
//
//		{
//			auto mixture = std::make_shared<MixtureSetUp>("Mixture2");
//			mixture->setPrimers("ama1_S0_Sub0_mip2", "CAATATTTAAAAGATGGAGGTTTTGCT",
//					"ATCTTCACAATTTCCATCGACCCAT");
//			mixture->addAbundance("PfML01", 10);
//			mixture->addAbundance("PfDd2", 20);
//			mixture->addAbundance("PfGB4-PfTG01", 30);
//			mixture->addAbundance("PfKH02", 40);
//			mixture->setReverseBarcode("MID08", "CTCGCGTGTC");
//			sampleSet->addMixture(mixture);
//		}
//		{
//			auto mixture = std::make_shared<MixtureSetUp>("Mixture3");
//			mixture->setPrimers("ama1_S0_Sub0_mip2", "CAATATTTAAAAGATGGAGGTTTTGCT",
//					"ATCTTCACAATTTCCATCGACCCAT");
//			mixture->addAbundance("PfML01", 10);
//			mixture->addAbundance("PfDd2", 20);
//			mixture->addAbundance("PfGB4-PfTG01", 30);
//			mixture->addAbundance("PfKH02", 40);
//			mixture->setReverseBarcode("MID10", "TCTCTATGCG");
//
//			sampleSet->addMixture(mixture);
//		}
//
//		{
//			auto mixture = std::make_shared<MixtureSetUp>("Mixture4");
//			mixture->setPrimers("ama1_S0_Sub0_mip2", "CAATATTTAAAAGATGGAGGTTTTGCT",
//					"ATCTTCACAATTTCCATCGACCCAT");
//			mixture->addAbundance("PfML01", 10);
//			mixture->addAbundance("PfDd2", 20);
//			mixture->addAbundance("PfGB4-PfTG01", 30);
//			mixture->addAbundance("PfKH02", 40);
//			mixture->setReverseBarcode("MID11", "TGATACGTCT");
//			sampleSet->addMixture(mixture);
//		}
//		lSetup.addSample(sampleSet);
//	}
	{
		auto sampleSet = std::make_shared<SampleSetup>("Sample3-BarcodeBothEndsDifferent");
		{
			auto mixture = std::make_shared<MixtureSetUp>("Mixture1");
			mixture->setPrimers("ama1_S0_Sub0_mip2", "CAATATTTAAAAGATGGAGGTTTTGCT",
					"ATCTTCACAATTTCCATCGACCCAT");
			mixture->addAbundance("PfML01", 10);
			mixture->addAbundance("PfDd2", 10);
			mixture->addAbundance("PfGB4-PfTG01", 10);
			mixture->addAbundance("PfKH02", 10);
			mixture->setForwardBarcode("MID07", "CGAGAGATAC");
			mixture->setReverseBarcode("MID07", "CGTGTCTCTA");
			mixture->forwardBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
			mixture->reverseBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
			sampleSet->addMixture(mixture);
		}
		{
			auto mixture = std::make_shared<MixtureSetUp>("Mixture2");
			mixture->setPrimers("ama1_S0_Sub0_mip2", "CAATATTTAAAAGATGGAGGTTTTGCT",
					"ATCTTCACAATTTCCATCGACCCAT");
			mixture->addAbundance("PfML01", 10);
			mixture->addAbundance("PfDd2", 20);
			mixture->addAbundance("PfGB4-PfTG01", 30);
			mixture->addAbundance("PfKH02", 40);
			mixture->setForwardBarcode("MID08", "TCACGTACTA");
			mixture->setReverseBarcode("MID08", "CTCGCGTGTC");
			mixture->forwardBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
			mixture->reverseBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
			sampleSet->addMixture(mixture);
		}
		{
			auto mixture = std::make_shared<MixtureSetUp>("Mixture3");
			mixture->setPrimers("ama1_S0_Sub0_mip2", "CAATATTTAAAAGATGGAGGTTTTGCT",
					"ATCTTCACAATTTCCATCGACCCAT");
			mixture->addAbundance("PfML01", 10);
			mixture->addAbundance("PfDd2", 20);
			mixture->addAbundance("PfGB4-PfTG01", 30);
			mixture->addAbundance("PfKH02", 40);
			mixture->setForwardBarcode("MID10", "CGTCTAGTAC");
			mixture->setReverseBarcode("MID10", "TCTCTATGCG");
			mixture->forwardBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
						mixture->reverseBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
			sampleSet->addMixture(mixture);
		}
		{
			auto mixture = std::make_shared<MixtureSetUp>("Mixture4");
			mixture->setPrimers("ama1_S0_Sub0_mip2", "CAATATTTAAAAGATGGAGGTTTTGCT",
					"ATCTTCACAATTTCCATCGACCCAT");
			mixture->addAbundance("PfML01", 10);
			mixture->addAbundance("PfDd2", 20);
			mixture->addAbundance("PfGB4-PfTG01", 30);
			mixture->addAbundance("PfKH02", 40);
			mixture->setForwardBarcode("MID11", "TCTACGTAGC");
			mixture->setReverseBarcode("MID11", "TGATACGTCT");
			mixture->forwardBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
						mixture->reverseBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
			sampleSet->addMixture(mixture);
		}
		{
			auto mixture = std::make_shared<MixtureSetUp>("Mixture5");
			mixture->setPrimers("ama1_S0_Sub0_mip2", "CAATATTTAAAAGATGGAGGTTTTGCT",
					"ATCTTCACAATTTCCATCGACCCAT");
			mixture->addAbundance("PfML01", 2);
			mixture->addAbundance("PfDd2", 3);
			mixture->addAbundance("PfGB4-PfTG01", 5);
			mixture->addAbundance("PfKH02", 90);
			mixture->setForwardBarcode("MID01", "ACGAGTGCGT");
			mixture->setReverseBarcode("MID01", "GAGTGCGTCT");
			mixture->forwardBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
			sampleSet->addMixture(mixture);
		}

		{
			auto mixture = std::make_shared<MixtureSetUp>("Mixture6");
			mixture->setPrimers("ama1_S0_Sub0_mip2", "CAATATTTAAAAGATGGAGGTTTTGCT",
					"ATCTTCACAATTTCCATCGACCCAT");
			mixture->addAbundance("PfML01", 1);
			mixture->addAbundance("PfDd2", 1);
			mixture->addAbundance("PfGB4-PfTG01", 1);
			mixture->addAbundance("PfKH02", 99);
			mixture->setForwardBarcode("MID29", "ACTGTACAGT");
			mixture->setReverseBarcode("MID29", "CACGCTACGT");
			mixture->forwardBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
			sampleSet->addMixture(mixture);
		}
		lSetup.addSample(sampleSet);
	}
	{
		auto sampleSet = std::make_shared<SampleSetup>("Sample4-BarcodeBothEndsSame");
		{
			auto mixture = std::make_shared<MixtureSetUp>("Mixture1");
			mixture->setPrimers("ama1_S0_Sub0_mip2", "CAATATTTAAAAGATGGAGGTTTTGCT",
					"ATCTTCACAATTTCCATCGACCCAT");
			mixture->addAbundance("PfML01", 10);
			mixture->addAbundance("PfDd2", 10);
			mixture->addAbundance("PfGB4-PfTG01", 10);
			mixture->addAbundance("PfKH02", 10);
			mixture->setForwardBarcode("MID07", "CGTGTCTCTA");
			mixture->setReverseBarcode("MID07", "CGTGTCTCTA");
			mixture->forwardBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
						mixture->reverseBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
			sampleSet->addMixture(mixture);
		}
		{
			auto mixture = std::make_shared<MixtureSetUp>("Mixture2");
			mixture->setPrimers("ama1_S0_Sub0_mip2", "CAATATTTAAAAGATGGAGGTTTTGCT",
					"ATCTTCACAATTTCCATCGACCCAT");
			mixture->addAbundance("PfML01", 10);
			mixture->addAbundance("PfDd2", 20);
			mixture->addAbundance("PfGB4-PfTG01", 30);
			mixture->addAbundance("PfKH02", 40);
			mixture->setForwardBarcode("MID08", "CTCGCGTGTC");
			mixture->setReverseBarcode("MID08", "CTCGCGTGTC");
			mixture->forwardBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
						mixture->reverseBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
			sampleSet->addMixture(mixture);
		}
		{
			auto mixture = std::make_shared<MixtureSetUp>("Mixture3");
			mixture->setPrimers("ama1_S0_Sub0_mip2", "CAATATTTAAAAGATGGAGGTTTTGCT",
					"ATCTTCACAATTTCCATCGACCCAT");
			mixture->addAbundance("PfML01", 10);
			mixture->addAbundance("PfDd2", 20);
			mixture->addAbundance("PfGB4-PfTG01", 30);
			mixture->addAbundance("PfKH02", 40);
			mixture->setForwardBarcode("MID10", "TCTCTATGCG");
			mixture->setReverseBarcode("MID10", "TCTCTATGCG");
			mixture->forwardBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
						mixture->reverseBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
			sampleSet->addMixture(mixture);
		}
		{
			auto mixture = std::make_shared<MixtureSetUp>("Mixture4");
			mixture->setPrimers("ama1_S0_Sub0_mip2", "CAATATTTAAAAGATGGAGGTTTTGCT",
					"ATCTTCACAATTTCCATCGACCCAT");
			mixture->addAbundance("PfML01", 10);
			mixture->addAbundance("PfDd2", 20);
			mixture->addAbundance("PfGB4-PfTG01", 30);
			mixture->addAbundance("PfKH02", 40);
			mixture->setForwardBarcode("MID11", "TGATACGTCT");
			mixture->setReverseBarcode("MID11", "TGATACGTCT");
			mixture->forwardBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
						mixture->reverseBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
			sampleSet->addMixture(mixture);
		}
		{
			auto mixture = std::make_shared<MixtureSetUp>("Mixture5");
			mixture->setPrimers("ama1_S0_Sub0_mip2", "CAATATTTAAAAGATGGAGGTTTTGCT",
					"ATCTTCACAATTTCCATCGACCCAT");
			mixture->addAbundance("PfML01", 2);
			mixture->addAbundance("PfDd2", 3);
			mixture->addAbundance("PfGB4-PfTG01", 5);
			mixture->addAbundance("PfKH02", 90);
			mixture->setForwardBarcode("MID01", "ACGAGTGCGT");
			mixture->setReverseBarcode("MID01", "ACGAGTGCGT");
			mixture->forwardBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
			sampleSet->addMixture(mixture);
		}
		{
			auto mixture = std::make_shared<MixtureSetUp>("Mixture6");
			mixture->setPrimers("ama1_S0_Sub0_mip2", "CAATATTTAAAAGATGGAGGTTTTGCT",
					"ATCTTCACAATTTCCATCGACCCAT");
			mixture->addAbundance("PfML01", 1);
			mixture->addAbundance("PfDd2", 1);
			mixture->addAbundance("PfGB4-PfTG01", 1);
			mixture->addAbundance("PfKH02", 99);
			mixture->setForwardBarcode("MID29", "ACTGTACAGT");
			mixture->setReverseBarcode("MID29", "ACTGTACAGT");
			mixture->forwardBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
			sampleSet->addMixture(mixture);
		}
		lSetup.addSample(sampleSet);
	}

	{
		auto sampleSet = std::make_shared<SampleSetup>("Sample5-BarcodeBothEndsRComp");
		{
			auto mixture = std::make_shared<MixtureSetUp>("Mixture1");
			mixture->setPrimers("ama1_S0_Sub0_mip2", "CAATATTTAAAAGATGGAGGTTTTGCT",
					"ATCTTCACAATTTCCATCGACCCAT");
			mixture->addAbundance("PfML01", 10);
			mixture->addAbundance("PfDd2", 10);
			mixture->addAbundance("PfGB4-PfTG01", 10);
			mixture->addAbundance("PfKH02", 10);
			mixture->setForwardBarcode("MID07", "CGTGTCTCTA");
			mixture->setReverseBarcode("MID07", "TAGAGACACG");
			mixture->forwardBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
						mixture->reverseBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
			sampleSet->addMixture(mixture);
		}
		{
			auto mixture = std::make_shared<MixtureSetUp>("Mixture2");
			mixture->setPrimers("ama1_S0_Sub0_mip2", "CAATATTTAAAAGATGGAGGTTTTGCT",
					"ATCTTCACAATTTCCATCGACCCAT");
			mixture->addAbundance("PfML01", 10);
			mixture->addAbundance("PfDd2", 20);
			mixture->addAbundance("PfGB4-PfTG01", 30);
			mixture->addAbundance("PfKH02", 40);
			mixture->setForwardBarcode("MID08", "CTCGCGTGTC");
			mixture->setReverseBarcode("MID08", "GACACGCGAG");
			mixture->forwardBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
						mixture->reverseBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
			sampleSet->addMixture(mixture);
		}
		{
			auto mixture = std::make_shared<MixtureSetUp>("Mixture3");
			mixture->setPrimers("ama1_S0_Sub0_mip2", "CAATATTTAAAAGATGGAGGTTTTGCT",
					"ATCTTCACAATTTCCATCGACCCAT");
			mixture->addAbundance("PfML01", 10);
			mixture->addAbundance("PfDd2", 20);
			mixture->addAbundance("PfGB4-PfTG01", 30);
			mixture->addAbundance("PfKH02", 40);
			mixture->setForwardBarcode("MID10", "TCTCTATGCG");
			mixture->setReverseBarcode("MID10", "CGCATAGAGA");
			mixture->forwardBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
						mixture->reverseBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
			sampleSet->addMixture(mixture);
		}
		{
			auto mixture = std::make_shared<MixtureSetUp>("Mixture4");
			mixture->setPrimers("ama1_S0_Sub0_mip2", "CAATATTTAAAAGATGGAGGTTTTGCT",
					"ATCTTCACAATTTCCATCGACCCAT");
			mixture->addAbundance("PfML01", 10);
			mixture->addAbundance("PfDd2", 20);
			mixture->addAbundance("PfGB4-PfTG01", 30);
			mixture->addAbundance("PfKH02", 40);
			mixture->setForwardBarcode("MID11", "TGATACGTCT");
			mixture->setReverseBarcode("MID11", "AGACGTATCA");
			mixture->forwardBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
						mixture->reverseBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
			sampleSet->addMixture(mixture);
		}
		{
			auto mixture = std::make_shared<MixtureSetUp>("Mixture5");
			mixture->setPrimers("ama1_S0_Sub0_mip2", "CAATATTTAAAAGATGGAGGTTTTGCT",
					"ATCTTCACAATTTCCATCGACCCAT");
			mixture->addAbundance("PfML01", 2);
			mixture->addAbundance("PfDd2", 3);
			mixture->addAbundance("PfGB4-PfTG01", 5);
			mixture->addAbundance("PfKH02", 90);
			mixture->setForwardBarcode("MID01", "ACGAGTGCGT");
			mixture->setReverseBarcode("MID01", "ACGCACTCGT");
			mixture->forwardBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
			sampleSet->addMixture(mixture);
		}
		{
			auto mixture = std::make_shared<MixtureSetUp>("Mixture6");
			mixture->setPrimers("ama1_S0_Sub0_mip2", "CAATATTTAAAAGATGGAGGTTTTGCT",
					"ATCTTCACAATTTCCATCGACCCAT");
			mixture->addAbundance("PfML01", 1);
			mixture->addAbundance("PfDd2", 1);
			mixture->addAbundance("PfGB4-PfTG01", 1);
			mixture->addAbundance("PfKH02", 99);
			mixture->setForwardBarcode("MID29", "ACTGTACAGT");
			mixture->setReverseBarcode("MID29", "ACTGTACAGT");
			mixture->forwardBarcode_->randomPrecedingBases_ = simPars.barcodeRandomPrecedingBases_;
			sampleSet->addMixture(mixture);
		}

		lSetup.addSample(sampleSet);
	}
	out << lSetup.toJson() << std::endl;

//	PfML01 vs PfDd2 1
//	PfGB4-PfTG01 vs PfKH02 3
//	PfGA01 vs Pf7G8 5
	return 0;
}

//
//
//int readSimulatorRunner::simMultipleMixture(const njh::progutils::CmdArgs & inputCommands) {
//
//  readSimulatorSetUp setUp(inputCommands);
//
//	uint32_t numThreads = 2;
//	bool singleEnd = false;
//	bool nonGz = false;
//	bfs::path idFile = "";
//	bfs::path referenceFile = "";
//	bfs::path librarySetUpFile = "";
//
//	bfs::path illuminaProfileDir = "";
//	setUp.processVerbose();
//	setUp.processDebug();
//	setUp.setOption(librarySetUpFile, "--librarySetUpFile", "Library Set Up File", true);
//	setUp.processReadInNames(VecStr{"--fasta", "--fastagz"},true);
//	setUp.setOption(numThreads, "--numThreads", "Number of Threads to Use");
//	setUp.setOption(singleEnd, "--singleEnd", "Single End");
//	setUp.setOption(nonGz, "--nonGz", "do not compress the output fastqs");
//	setUp.setOption(illuminaProfileDir, "--illuminaProfileDir", "Illumina Profile Dir", true);
//	setUp.processDirectoryOutputName("simMultipleMixture_TODAY", true);
//	setUp.finishSetUp(std::cout);
//	setUp.startARunLog(setUp.pars_.directoryName_);
//	setUp.writeParametersFile(setUp.pars_.directoryName_ + "parameters.tab.txt", false, true);
//
//	std::unordered_map<std::string, std::shared_ptr<seqInfo>> refSeqs;
//
//
//	{
//		SeqInput reader(setUp.pars_.ioOptions_);
//		reader.openIn();
//		seqInfo seq;
//		while(reader.readNextRead(seq)){
//			refSeqs[seq.name_] = std::make_shared<seqInfo>(seq);
//		}
//	}
//	Json::Value librarySetUp = njh::json::parseFile(librarySetUpFile.string());
//
//	njh::json::MemberChecker libraryChecker(librarySetUp);
//	libraryChecker.failMemberCheckThrow(LibrarySetup::jsonMembers(), __PRETTY_FUNCTION__);
//	LibrarySetup::SimLibrarySetupPars simPars(librarySetUp["pars"]);
//	LibrarySetup lSetup(librarySetUp["name"].asString(), simPars, librarySetUp, getVectorOfMapKeys(refSeqs));
//
//
//	OutOptions idOutOpts(njh::files::make_path(setUp.pars_.directoryName_, "ids.tab.txt"));
//
//	lSetup.ids_->writeIdFile(OutOptions(njh::files::make_path(setUp.pars_.directoryName_, "ids.tab.txt")));
//
//	OutOptions outOptsSetup(njh::files::make_path(setUp.pars_.directoryName_, "setup.json"));
//	OutputStream outOptsSetupOut(outOptsSetup);
//	outOptsSetupOut << lSetup.toJson() << std::endl;
//
//	bfs::path fastqDirectory =njh::files::make_path(setUp.pars_.directoryName_, "fastq");
//	njh::files::makeDir(njh::files::MkdirPar{fastqDirectory});
//	auto sampleKeys = njh::getVecOfMapKeys(lSetup.samples_);
//	njh::sort(sampleKeys);
//
//	std::atomic_uint sampleCountAtom{0};
//
//	njh::concurrent::LockableQueue<std::string> samplesQueue(sampleKeys);
//	std::string extension = nonGz ? ".fastq" : ".fastq.gz";
//
//	auto simSample = [&samplesQueue,&extension,&sampleCountAtom,&fastqDirectory,&lSetup,&illuminaProfileDir,&singleEnd,&refSeqs](){
//		std::string sampleKey = "";
//		njh::randObjectGen<char,uint32_t> baseRGen({'A', 'C', 'G', 'T'}, {1,1,1,1});
//		RoughIlluminaSimulator simulator(illuminaProfileDir);
//		njh::randomGenerator rGen;
//
//		while(samplesQueue.getVal(sampleKey)){
//			uint32_t sampleCount = sampleCountAtom++;
//			OutOptions targetOutOpts(njh::files::make_path(fastqDirectory, njh::pasteAsStr(sampleKey, "_S", sampleCount + 1, extension)));
//			OutOptions r1OutOpts(njh::files::make_path(fastqDirectory, njh::pasteAsStr(sampleKey, "_S", sampleCount + 1, "_R1_001" + extension)));
//			OutOptions r2OutOpts(njh::files::make_path(fastqDirectory, njh::pasteAsStr(sampleKey, "_S", sampleCount + 1, "_R2_001" + extension)));
//			std::shared_ptr<OutputStream> targetOut;
//
//			std::shared_ptr<OutputStream> r1Out;
//			std::shared_ptr<OutputStream> r2Out;
//			if(singleEnd) {
//				targetOut = std::make_shared<OutputStream>(targetOutOpts);
//			} else {
//				r1Out = std::make_shared<OutputStream>(r1OutOpts);
//				r2Out = std::make_shared<OutputStream>(r2OutOpts);
//			}
//
//
//			/**@todo add back in the ability to have different Illumina barcodes in overhangs */
//	//		std::string sampleAdapter1 = defaultAdapter1;
//	//		std::string sampleAdapter2 = defaultAdapter2;
//	//		for(const auto pos : iter::range(sampleAdapter1.size())){
//	//			if('N' == sampleAdapter1[pos]){
//	//				sampleAdapter1[pos] = baseRGen.genObj();
//	//			}
//	//		}
//	//		for(const auto pos : iter::range(sampleAdapter2.size())){
//	//			if('N' == sampleAdapter2[pos]){
//	//				sampleAdapter2[pos] = baseRGen.genObj();
//	//			}
//	//		}
//			for(const auto & mixture : lSetup.samples_.at(sampleKey)->mixtures_ ){
//				VecStr refNames;
//				std::vector<double> amounts;
//				for(const auto & abun : mixture.second->expectedAbundances_){
//					refNames.emplace_back(abun.first);
//					amounts.emplace_back(abun.second);
//				}
//				njh::randObjectGen<std::string,double> refNameGen(refNames, amounts);
//				for(uint32_t readNumber = 0; readNumber < mixture.second->finalReadAmount_; ++readNumber){
//				//for(uint32_t readNumber = 0; readNumber < lSetup.pars_.sampleReadAmount_; ++readNumber){
//					MetaDataInName nameMeta;
//					auto refName = refNameGen.genObj();
//
//					nameMeta.addMeta("Mixture", mixture.second->name_);
//					nameMeta.addMeta("Sample", sampleKey);
//					nameMeta.addMeta("readNumber", readNumber);
//					nameMeta.addMeta("refName", refName);
//					seqInfo targetSeq("", refSeqs.at(refName)->seq_);
//
//					uint32_t forwardPrimerPosition = 0;
//					uint32_t reversePrimerPosition = 0;
//
//					uint32_t forwardBarcodePosition = 0;
//					uint32_t reverseBarcodePosition = 0;
//
//					if(nullptr != mixture.second->primers_){
////						if(!njh::in(mixture.second->primers_->name_, primerNames)){
////							primerOut << mixture.second->primers_->name_
////									<< "\t" << mixture.second->primers_->forward_
////									<< "\t" << mixture.second->primers_->reverse_
////									<< std::endl;
////							primerNames.emplace_back(mixture.second->primers_->name_);
////						}
//						if(!lSetup.pars_.noAddPrimers_){
//							/**@todo add allowing for ambigious bases in primers*/
//							targetSeq.prepend(mixture.second->primers_->forward_);
//							reversePrimerPosition = targetSeq.seq_.size();
//							targetSeq.append(mixture.second->primers_->reverse_3_5_);
//						}
//
//
//
//						if(mixture.second->primers_->forward_randomPrecedingBases_ > 0){
//							auto minBaseAmount = mixture.second->primers_->randomOneLength_ ? mixture.second->primers_->forward_randomPrecedingBases_: 0;
//							uint32_t numFBase = rGen.unifRand<uint32_t>(minBaseAmount, mixture.second->primers_->forward_randomPrecedingBases_ + 1);
//							if(numFBase > 0){
//								targetSeq.prepend(njh::pasteAsStr(baseRGen.genObjs(numFBase)));
//								reversePrimerPosition += numFBase;
//								forwardPrimerPosition += numFBase;
//							}
//						}
//						if(mixture.second->primers_->reverse_randomPrecedingBases_ > 0){
//							auto minBaseAmount = mixture.second->primers_->randomOneLength_ ? mixture.second->primers_->reverse_randomPrecedingBases_: 0;
//							uint32_t numRBase = rGen.unifRand<uint32_t>(minBaseAmount, mixture.second->primers_->reverse_randomPrecedingBases_ + 1);
//							if(numRBase > 0){
//								targetSeq.append(njh::pasteAsStr(baseRGen.genObjs(numRBase)));
//							}
//						}
//					}
//					if(nullptr != mixture.second->forwardBarcode_ || nullptr != mixture.second->reverseBarcode_){
//						std::string midName = nullptr != mixture.second->forwardBarcode_ ? mixture.second->forwardBarcode_->name_ : mixture.second->reverseBarcode_->name_;
////						if(!njh::in(midName, midNames)){
////							midNames.emplace_back(midName);
////							midOut << midName;
////							if(nullptr != mixture.second->forwardBarcode_){
////								midOut << "\t" << mixture.second->forwardBarcode_->barcode_;
////							}else{
////								midOut << "\t";
////							}
////							if(nullptr != mixture.second->reverseBarcode_){
////								midOut << "\t" << mixture.second->reverseBarcode_->barcode_;
////							}
////							midOut << std::endl;
////						}
//					}
//					if(nullptr != mixture.second->forwardBarcode_){
//						targetSeq.prepend(mixture.second->forwardBarcode_->barcode_);
//						reversePrimerPosition += mixture.second->forwardBarcode_->barcode_.size();
//						forwardPrimerPosition += mixture.second->forwardBarcode_->barcode_.size();
//						if(mixture.second->forwardBarcode_->randomPrecedingBases_ > 0){
//							auto minBaseAmount = mixture.second->forwardBarcode_->randomOneLength_ ? mixture.second->forwardBarcode_->randomPrecedingBases_: 0;
//							uint32_t numFBase = rGen.unifRand<uint32_t>(minBaseAmount, mixture.second->forwardBarcode_->randomPrecedingBases_ + 1);
//							if(numFBase > 0){
//								targetSeq.prepend(njh::pasteAsStr(baseRGen.genObjs(numFBase)));
//								reversePrimerPosition += numFBase;
//								forwardPrimerPosition += numFBase;
//								forwardBarcodePosition+= numFBase;
//							}
//						}
//					}
//					if(nullptr != mixture.second->reverseBarcode_){
//						reverseBarcodePosition = len(targetSeq);
//						targetSeq.append(mixture.second->reverseBarcode_->barcode_3_5_);
//						if(mixture.second->reverseBarcode_->randomPrecedingBases_ > 0){
//							auto minBaseAmount = mixture.second->reverseBarcode_->randomOneLength_ ? mixture.second->reverseBarcode_->randomPrecedingBases_: 0;
//							uint32_t numFBase = rGen.unifRand<uint32_t>(minBaseAmount, mixture.second->reverseBarcode_->randomPrecedingBases_ + 1);
//							if(numFBase > 0){
//								targetSeq.append(njh::pasteAsStr(baseRGen.genObjs(numFBase)));
//							}
//						}
//					}
//					//blunt ending
//					if(lSetup.pars_.addBluntEndingArtifact_){
//						if('A' == targetSeq.seq_[0] && rGen() < lSetup.pars_.bluntEndingArtifactChance_){
//							readVecTrimmer::trimOffForwardBases(targetSeq, 1);
//							forwardPrimerPosition = forwardPrimerPosition == 0 ? forwardPrimerPosition: forwardPrimerPosition - 1;
//							forwardBarcodePosition = forwardBarcodePosition == 0 ? forwardBarcodePosition: forwardBarcodePosition - 1;
//							reverseBarcodePosition -= 1;
//							reversePrimerPosition -= 1;
//							nameMeta.addMeta("frontEndBluntEndArtifact", true);
//						}else{
//							nameMeta.addMeta("frontEndBluntEndArtifact", false);
//						}
//						if ('T' == targetSeq.seq_.back()
//								&& rGen() < lSetup.pars_.bluntEndingArtifactChance_) {
//							readVecTrimmer::trimOffEndBases(targetSeq, 1);
//							nameMeta.addMeta("backEndBluntEndArtifact", true);
//						} else {
//							nameMeta.addMeta("backEndBluntEndArtifact", false);
//						}
//					}
//					//positions
//					if(nullptr != mixture.second->primers_){
//						nameMeta.addMeta("forwardPrimerPosition", forwardPrimerPosition);
//						nameMeta.addMeta("reversePrimerPosition", reversePrimerPosition);
//					}
//					if(nullptr != mixture.second->forwardBarcode_){
//						nameMeta.addMeta("forwardBarcodePosition", forwardBarcodePosition);
//					}
//					if(nullptr != mixture.second->reverseBarcode_){
//						nameMeta.addMeta("reverseBarcodePosition", reverseBarcodePosition);
//					}
//					//complement
//					if(lSetup.pars_.addReverseComplement_){
//						bool complement = rGen.unifRand(0,2) == 0;
//						nameMeta.addMeta("complement", complement);
//						if(complement){
//							targetSeq.reverseComplementRead(false, true);
//						}
//					}
//					//length
//					nameMeta.addMeta("targetSeqLength", len(targetSeq));
//
//					targetSeq.name_ = nameMeta.createMetaName();
//
//					if(singleEnd){
//						//target
//						targetSeq.outPutSeq(*targetOut);
//					}else{
//						{
//							//r1
//							auto subSeq = len(targetSeq) > lSetup.pars_.pairedEndLength_? targetSeq.getSubRead(0, lSetup.pars_.pairedEndLength_) : targetSeq;
//							simulator.simR1(subSeq, lSetup.pars_.pairedEndLength_).outPutFastq(*r1Out);
//						}
//						{
//							//r2
//							targetSeq.reverseComplementRead(false, true);
//							auto subSeq = len(targetSeq) > lSetup.pars_.pairedEndLength_? targetSeq.getSubRead(0, lSetup.pars_.pairedEndLength_) : targetSeq;
//							simulator.simR2(subSeq, lSetup.pars_.pairedEndLength_).outPutFastq(*r2Out);
//						}
//					}
//				}
//			}
//		}
//	};
//
//	std::vector<std::thread> threads;
//	for(uint32_t t = 0; t < numThreads; ++t){
//		threads.emplace_back(std::thread(simSample));
//	}
//	njh::concurrent::joinAllJoinableThreads(threads);
//
//	return 0;
//}
//




} // namespace njhseq
