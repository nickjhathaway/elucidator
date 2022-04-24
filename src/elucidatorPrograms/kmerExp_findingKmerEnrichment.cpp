

#include "kmerExp.hpp"
#include "elucidator/objects/MiscUtility/GenomeSeqSearch.hpp"
#include "elucidator/objects/BioDataObject.h"
#include <njhseq/objects/Meta.h>
#include "elucidator/objects/kmerUtils/SimpleKmerHash.hpp"
#include <njhseq/PopulationGenetics/PopGenCalcs.hpp>

namespace njhseq {

int kmerExpRunner::findingKmerEnrichment(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path sampleMetaFnp;
	std::vector<bfs::path> inputFiles;
	uint32_t kmerLength = 19;
	uint32_t numThreads = 1;
	uint32_t minGroupSize = 1;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(numThreads, "--numThreads", "number Threads");
	setUp.setOption(minGroupSize, "--minGroupSize", "minGroupSize");

	setUp.setOption(kmerLength, "--kmerLength", "kmer Length");
	setUp.setOption(inputFiles, "--inputFiles", "input files", true);
	setUp.setOption(sampleMetaFnp, "--sampleMetaFnp", "Sample Meta Fnp to find enrichment of", true);

	setUp.processDirectoryOutputName("findingKmerEnrichment_TODAY", true);
	setUp.finishSetUp(std::cout);

	setUp.startARunLog(setUp.pars_.directoryName_);

	std::set<std::string> allSamples;
	for(const auto & inputFile : inputFiles){
		std::string sample = inputFile.filename().replace_extension("").string();
		if(njh::endsWith(inputFile.string(), ".gz")){
			sample = bfs::path(sample).replace_extension("").string();
		}
		allSamples.emplace(sample);
	}
	MultipleGroupMetaData meta(sampleMetaFnp);
	std::set<std::string> groupsWith2Levels;
	for(const auto & group : meta.groupData_){
		uint32_t count = 0;
		for(const auto & group : group.second->subGroupsLevels_){
			if("NA" != group){
				++count;
			}
		}
		if(count == 2){
			groupsWith2Levels.emplace(group.first);
		}
	}

	std::cout << "groupsWith2Levels: " << njh::conToStr(groupsWith2Levels, ",") << std::endl;
	std::unordered_map<uint64_t, std::unordered_set<std::string>> kmersSamples;
	std::mutex kmersSampleMut;

	njh::concurrent::LockableQueue<bfs::path> inputFilesQueue(inputFiles);

	std::function<void()> count = [&inputFilesQueue,&kmersSamples,&kmersSampleMut,&kmerLength](){
		std::unordered_map<uint64_t, std::unordered_set<std::string>> currentThread_kmersSamples;

		seqInfo seq;
		SimpleKmerHash khasher;
		bfs::path inputFile;
		while(inputFilesQueue.getVal(inputFile)){
			auto inOpts = SeqIOOptions(inputFile, SeqIOOptions::getInFormatFromFnp(inputFile));
			SeqInput reader(inOpts);
			reader.openIn();
			std::string sample = inputFile.filename().replace_extension("").string();
			if(njh::endsWith(inputFile.string(), ".gz")){
				sample = bfs::path(sample).replace_extension("").string();
			}
			while(reader.readNextRead(seq)){
				for (const auto pos : iter::range(seq.seq_.size() + 1 - kmerLength)) {
					currentThread_kmersSamples[khasher.hash(seq.seq_, pos, kmerLength)].emplace(sample);
					currentThread_kmersSamples[khasher.revCompHash(seq.seq_, pos, kmerLength)].emplace(sample);
				}
			}
		}
		{
			std::lock_guard<std::mutex> lock(kmersSampleMut);
			for(const auto & kmers : currentThread_kmersSamples){
				kmersSamples[kmers.first].insert(kmers.second.begin(), kmers.second.end())	;
			}
		}
	};

	njh::concurrent::runVoidFunctionThreaded(count, numThreads);

	OutputStream outCounts(njh::files::make_path(setUp.pars_.directoryName_, "counts.tab.txt.gz"));
	outCounts << "kmer\tsampleCount" << std::endl;

	OutputStream fisherTestsOut(njh::files::make_path(setUp.pars_.directoryName_, "fisherTests.tab.txt.gz"));
	fisherTestsOut << "kmer\tgroup\tpositiveGroup\toddsRatio\tlowerConf\tupperConf\tpValue\tTP\tFP\tFN\tTN" << std::endl;
	for(const auto & kmers : kmersSamples){
		outCounts << kmers.first << "\t" << kmers.second.size() << std::endl;
		if(kmers.second.size() + minGroupSize <= inputFiles.size()){
			for(const auto & group : groupsWith2Levels){

				std::unordered_map<std::string, uint32_t> groupIndex;
				std::unordered_map<uint32_t, std::string> indexToGroup;
				{
					uint32_t groupLevelCount = 0;
					for(const auto & groupLevel : meta.groupData_.at(group)->subGroupsLevels_){
						if("NA"!= groupLevel){
							groupIndex[groupLevel] = groupLevelCount;
							indexToGroup[groupLevelCount] = groupLevel;
							++groupLevelCount;
						}
					}
				}
				std::vector<std::vector<uint32_t>> counts(2, std::vector<uint32_t>(2, 0));
//				std::cout << "counts: " << std::endl;
//				for(const auto & count : counts){
//					printVector(count, " ", std::cout);
//				}
				for(const auto & samp : kmers.second){
					++counts[0][groupIndex[meta.groupData_.at(group)->getGroupForSample(samp)]];
				}
				for(const auto & samp : allSamples){
					if(!njh::in(samp, kmers.second)){
						++counts[1][groupIndex[meta.groupData_.at(group)->getGroupForSample(samp)]];
					}
				}


				PopGenCalculator::FisherExactFor2x2::FisherExactFor2x2Input fisherInput;
				fisherInput.TP = counts[0][0];
				fisherInput.FP = counts[1][0];
				fisherInput.FN = counts[0][1];
				fisherInput.TN = counts[1][1];

//				std::cout << "kmers: " << kmers.first << std::endl;
//				std::cout << "allSamples: ";
//				std::cout << njh::conToStr(allSamples, ", ") << std::endl;
//				std::cout << "samplesWithKmers: ";
//				std::cout << njh::conToStr(kmers.second, ", ") << std::endl;
//
//				std::cout << "counts: " << std::endl;
//				for(const auto & count : counts){
//					printVector(count, " ", std::cout);
//				}
				if(fisherInput.TP > 2){
					auto res = PopGenCalculator::FisherExactFor2x2::runFisherExactOn2x2(fisherInput);
//				std::cout << "odds ratio: " << res.oddsRatio_ << std::endl;
//				std::cout << "lower " << roundDecPlaces(100*fisherInput.confInterval, 2)<< "% confidence: " << res.lowerConfInterval_ << std::endl;
//				std::cout << "upper " << roundDecPlaces(100*fisherInput.confInterval, 2)<< "% confidence: " << res.upperConfInterval_ << std::endl;
//				std::cout << "p-value: " << res.pValue_ << std::endl; 	TP	FP	FN	TN
					fisherTestsOut << kmers.first
												 << "\t" << group
												 << "\t" << indexToGroup[0]
												 << "\t" << res.oddsRatio_
												 << "\t" << res.lowerConfInterval_
												 << "\t" << res.upperConfInterval_
												 << "\t" << res.pValue_
												 << "\t" << fisherInput.TP
												 << "\t" << fisherInput.FP
												 << "\t" << fisherInput.FN
												 << "\t" << fisherInput.TN
												 << std::endl;
				}
			}
		}
	}

	return 0;
}

} // namespace njhseq