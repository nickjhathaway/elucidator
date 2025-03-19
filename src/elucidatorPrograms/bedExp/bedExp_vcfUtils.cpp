//
// Created by Nicholas Hathaway on 1/9/25.
//


#include "bedExp.hpp"
#include <njhseq/objects/BioDataObject.h>

#include "elucidator/objects/BioDataObject.h"
#include "elucidator/BioRecordsUtils/BedUtility.hpp"
#include <njhseq/objects/Gene/TranslatorByAlignment.hpp>



namespace njhseq {

int bedExpRunner::vcfToBed(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path vcfFile;
	bfs::path intersectWithBed;
	OutOptions outOpts;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.setOption(vcfFile, "--vcfFile", "vcfFile", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	OutputStream out(outOpts);

	VCFOutput vcf = VCFOutput::readInHeader(vcfFile);
	InputStream in(vcfFile);
	std::string line;
	// uint32_t count = 0;
	while(njh::files::crossPlatGetline(in, line)) {
		if(line.front() != '#') {
			// std::cout << count++ << std::endl;
			out << vcf.processRecordLineForFixedData(line).genRegion().genBedRecordCore().toDelimStrWithExtra() << std::endl;
		}
	}

	return 0;
}

int bedExpRunner::printVcfSamples(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path vcfFile;
	bfs::path intersectWithBed;
	OutOptions outOpts;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.setOption(vcfFile, "--vcfFile", "vcfFile", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	OutputStream out(outOpts);

	VCFOutput vcf = VCFOutput::readInHeader(vcfFile);
	out << njh::conToStr(vcf.samples_, "\n") << std::endl;
	return 0;
}




int bedExpRunner::combineVcfs(const njh::progutils::CmdArgs & inputCommands) {
	std::vector<bfs::path> vcfFnps;
	bfs::path intersectWithBed;
	OutOptions outOpts;
	VCFOutput::comnbineVCFsPars combiningVcfPars;

	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.setOption(vcfFnps, "--vcfFnps", "vcfs", true);
	setUp.setOption(combiningVcfPars.ploidy, "--ploidy", "Ploidy to force for the sample for the vcf files");
	setUp.setOption(combiningVcfPars.doNotRescueVariantCallsAcrossTargets, "--doNotRescueVariantCallsAcrossTargets", "do Not Rescue Variant Calls Across Targets");
	setUp.setOption(combiningVcfPars.combinedOverlappingCallsAcrossTargets, "--combineOverlappingCallsAcrossTargets", "Rather than taking the best variant call for overlapping targets, sum them instead");
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	OutputStream out(outOpts);

	auto combined = VCFOutput::comnbineVCFs(vcfFnps, combiningVcfPars);
	combined.writeOutFixedAndSampleMeta(out);
	return 0;
}



}  //namespace njhseq


