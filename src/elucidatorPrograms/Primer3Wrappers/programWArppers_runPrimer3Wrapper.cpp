/*
 * programWArppers_runPrimer3Wrapper.cpp
 *
 *  Created on: Sep 7, 2021
 *      Author: nick
 */



#include "elucidatorPrograms/programWrappers.hpp"
#include "elucidator/objects/BioDataObject.h"
#include "elucidator/BamToolsUtils.h"

#include <SeekDeep/objects/TarAmpSetupUtils/PrimersAndMids.hpp>

#include <njhseq/concurrency/PairwisePairFactory.hpp>

namespace njhseq {




int programWrapperRunner::genPossiblePrimersWithPrimer3(const njh::progutils::CmdArgs & inputCommands){


	Primer3Runner::Primer3Options p3Opts;
	bfs::path vcfFnp;
	double majorAFVCFFreq = 0.99;


	bfs::path excludeRegionsFnp;
	bfs::path bedFnp;
	bfs::path twoBitFnp;

	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(excludeRegionsFnp, "--excludeRegions", "exclude Regions");
	setUp.setOption(majorAFVCFFreq, "--majorAFVCFFreq", "If supplying a VCF file than change SNPs with this frequency or more before designing primers");
	setUp.setOption(vcfFnp, "--vcfFnp", "VCF Fnp");

	p3Opts.setPrimaryOptions(setUp);
	p3Opts.setPrimerSizeOpts(setUp);
	p3Opts.setReturnOptions(setUp);
	p3Opts.task = "pick_primer_list";
	//setUp.setOption(task, "--task", "primer picking task, examples include:generic(default), pick_sequencing_primers, pick_primer_list");


	uint32_t extendRegion = std::max<uint32_t>(p3Opts.PRIMER_MAX_SIZE + 5, std::round(p3Opts.maxSize * 1.5));



	setUp.setOption(bedFnp, "--bedFnp", "genomic locations", true);
	setUp.setOption(twoBitFnp, "--twoBit", "two Bit file", true);

	setUp.processDirectoryOutputName("genPrimers_TODAY", true);

	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);

	std::vector<std::shared_ptr<Bed6RecordCore>> excludeRegions;

	if(bfs::exists(excludeRegionsFnp)){
		excludeRegions = getBeds(excludeRegionsFnp);
	}

	std::vector<VCFVariant> allVariants;
	if(bfs::exists(vcfFnp)){
		InputStream vcfIn(vcfFnp);
		std::string line;
		while(njh::files::crossPlatGetline(vcfIn, line)){
			if(njh::beginsWith(line, "#")){
				continue;
			}
			auto variants = VCFVariant::readVCFLine(line);
			njh::addConToVec(allVariants, variants);
		}
	}

	auto regions = getBeds(bedFnp);

	std::unordered_map<std::string, std::shared_ptr<Bed6RecordCore>> nameToRegion;
	std::unordered_map<std::string, std::shared_ptr<Bed6RecordCore>> templateToRegion;

	TwoBit::TwoBitFile tReader(twoBitFnp);
	auto chromLens = tReader.getSeqLens();
	std::set<std::string> regionNames;
	for(const auto & reg : regions){
		regionNames.emplace(reg->name_);
	}

	BedUtility::coordSort(excludeRegions, false);
	BedUtility::coordSort(regions, false);
	njh::sort(allVariants, [](const VCFVariant & reg1, const VCFVariant & reg2) {
				if(reg1.region_.chrom_ == reg2.region_.chrom_) {
					if(reg1.region_.start_ == reg2.region_.start_) {
						return reg1.region_.end_ < reg2.region_.end_;
					} else {
						return reg1.region_.start_ < reg2.region_.start_;
					}
				} else {
					return reg1.region_.chrom_ < reg2.region_.chrom_;
				}
			});
	OutputStream inputRegionsExpanded(njh::files::make_path(setUp.pars_.directoryName_, "inputRegionsExpanded.bed"));


	//SEQUENCE_TARGET, SEQUENCE_EXCLUDED_REGION
	{
		OutputStream primer3Input(njh::files::make_path(setUp.pars_.directoryName_, "primer3File.txt"));
		primer3Input << "PRIMER_THERMODYNAMIC_PARAMETERS_PATH=" << p3Opts.primer3ConfigPath.string() << std::endl;

		for(const auto & originalReg : regions){
			std::shared_ptr<Bed6RecordCore> reg = std::make_shared<Bed6RecordCore>(*originalReg);
			reg->strand_ = '+';
			std::vector<Bed6RecordCore> excludeRegions_withinRegion;
			std::vector<Bed6RecordCore> excludeRegions_relToRegion;

			std::vector<VCFVariant> variants_withinRegion;
			std::vector<VCFVariant> variants_relToRegion;
			bool lenIsScore = reg->score_ == reg->length();
			BedUtility::extendLeftRight(*reg, extendRegion, extendRegion, njh::mapAt(chromLens, reg->chrom_) );
			if(lenIsScore){
				reg->score_ = reg->length();
			}
			inputRegionsExpanded << reg->toDelimStrWithExtra() << std::endl;
			nameToRegion[reg->name_] = reg;

			{
				for(const auto & inter : excludeRegions){
					if(inter->chrom_ < reg->chrom_){
						//bother regions are sorted if we haven't reached this region's chromosome yet
						continue;
					}
					if(inter->chrom_ > reg->chrom_){
						//both regions are sorted so if we run into this we can break
						break;
					}
					if(inter->chrom_ == reg->chrom_ && inter->chromStart_ >= reg->chromEnd_){
						//both regions are sorted so if we run into this we can break
						break;
					}
					if(inter->chrom_ == reg->chrom_ && inter->chromEnd_ < reg->chromStart_){
						continue;
					}
					if(reg->overlaps(*inter, 1)){
						//check first if we already have this exact region for exclusion, this can happen if if input generated from variant calls and there is more than biallelic SNPs etc
						bool alreadyHave = false;
						for(const auto & already : excludeRegions_withinRegion){
							if(already.genUIDFromCoords() == inter->genUIDFromCoords()){
								alreadyHave = true;
								break;
							}
						}
						if(!alreadyHave){
							excludeRegions_withinRegion.emplace_back(*inter);
						}
					}
				}
				//book end the exclusion regions
				for(auto & exclude : excludeRegions_withinRegion){
					exclude.chromStart_ = std::max(exclude.chromStart_, reg->chromStart_);
					exclude.chromEnd_ = std::min(exclude.chromEnd_, reg->chromEnd_);
				}
			}

			//variant regions;
			{
				for(const auto & var : allVariants){
					if(var.region_.chrom_ < reg->chrom_){
						//bother regions are sorted if we haven't reached this region's chromosome yet
						continue;
					}
					if(var.region_.chrom_ > reg->chrom_){
						//both regions are sorted so if we run into this we can break
						break;
					}
					if(var.region_.chrom_ == reg->chrom_ && var.region_.start_ >= reg->chromEnd_){
						//both regions are sorted so if we run into this we can break
						break;
					}
					if(var.region_.chrom_ == reg->chrom_ && var.region_.end_ < reg->chromStart_){
						continue;
					}
					if(reg->overlaps(var.region_.genBed3RecordCore(), 1)){
						variants_withinRegion.emplace_back(var);
					}
				}
				//book end the exclusion regions
				for(auto & var : variants_withinRegion){
					var.region_.start_ = std::max<size_t>(var.region_.start_, reg->chromStart_);
					var.region_.end_ = std::min<size_t>(var.region_.end_, reg->chromEnd_);
				}
			}

			//convert to regions relative to the template region

			for(const auto & var : variants_withinRegion){
				auto modRegion = var;
				modRegion.region_.start_ = modRegion.region_.start_ - reg->chromStart_;
				modRegion.region_.end_ =   modRegion.region_.end_   - reg->chromStart_;
				variants_relToRegion.emplace_back(modRegion);
				if(var.freq_ < majorAFVCFFreq){
					//add variants
					excludeRegions_relToRegion.emplace_back(modRegion.region_.genBedRecordCore());
				}
			}



			for(const auto & exclude : excludeRegions_withinRegion){
				bool withinVCF = false;
				for(const auto & var : variants_withinRegion){
					if(var.region_.createUidFromCoords() == exclude.genUIDFromCoords()){
						withinVCF = true;
						break;
					}
				}
				if(!withinVCF){
					auto modRegion = exclude;
					modRegion.chromStart_ = modRegion.chromStart_ - reg->chromStart_;
					modRegion.chromEnd_ =   modRegion.chromEnd_   - reg->chromStart_;
					excludeRegions_relToRegion.emplace_back(modRegion);
				}
			}

			auto seqTemplate = GenomicRegion(*reg).extractSeq(tReader);
			//make any changes to the templates as needed
			for(const auto & var : variants_relToRegion){
				//only doing SNPs right now cause otherwise handling INDELs and the changes to the genomic locations would be kind of a nightmare
				if(var.vtype_ == VCFVariant::VarType::SNP && var.freq_ >= majorAFVCFFreq){
					seqTemplate.seq_[var.region_.start_] = var.variant_[0]; //change base
				}
			}
			std::string seqID = reg->name_;
			std::stringstream defaultPars;
			defaultPars << "SEQUENCE_TEMPLATE=" << seqTemplate.seq_<< std::endl;
			defaultPars << "SEQUENCE_ID=" << seqID<< std::endl;
			defaultPars << "PRIMER_FIRST_BASE_INDEX=0" << std::endl;
			defaultPars << "PRIMER_TASK=" << p3Opts.task<< std::endl;
			defaultPars << "PRIMER_SEQUENCING_LEAD=5" << std::endl;
			defaultPars << "PRIMER_EXPLAIN_FLAG=1" << std::endl;
			defaultPars << "PRIMER_PRODUCT_SIZE_RANGE=" << p3Opts.minSize << "-" << p3Opts.maxSize << std::endl;
			defaultPars << "PRIMER_NUM_RETURN=" << p3Opts.PRIMER_NUM_RETURN << std::endl;
			defaultPars << "PRIMER_MAX_SIZE=" << p3Opts.PRIMER_MAX_SIZE << std::endl;
			defaultPars << "PRIMER_MIN_SIZE=" << p3Opts.PRIMER_MIN_SIZE << std::endl;
			defaultPars << "PRIMER_OPT_SIZE=" << p3Opts.PRIMER_OPT_SIZE << std::endl;

			if(!excludeRegions_relToRegion.empty()){
				std::string SEQUENCE_EXCLUDED_REGION = "";
				for(const auto & exclude : excludeRegions_relToRegion){
					if("" != SEQUENCE_EXCLUDED_REGION){
						SEQUENCE_EXCLUDED_REGION += " ";
					}
					SEQUENCE_EXCLUDED_REGION += njh::pasteAsStr(exclude.chromStart_, ",", exclude.chromEnd_ - exclude.chromStart_);
				}
				defaultPars << "SEQUENCE_EXCLUDED_REGION=" << SEQUENCE_EXCLUDED_REGION << std::endl;
			}
			primer3Input << defaultPars.str();
			primer3Input << "="<< std::endl;
		}
	}
	std::string primer3Cmd = njh::pasteAsStr("cd ", setUp.pars_.directoryName_, " && ", "primer3_core ", "primer3File.txt > primer3_output.txt 2> primer3_error.txt");
	auto runOutput = njh::sys::run({primer3Cmd});
	OutputStream primer3RunLog(njh::files::make_path(setUp.pars_.directoryName_, "primer3RunLog.json"));
	primer3RunLog << runOutput.toJson() << std::endl;

	auto primer3ResultsFnp = njh::files::make_path(setUp.pars_.directoryName_, "primer3_output.txt");




	auto results = Primer3Runner::Primer3ResultsList::parsePrimer3OutputResults(primer3ResultsFnp, true);

	OutputStream primer3ResultsOut(njh::files::make_path(setUp.pars_.directoryName_, "primer3_results.tab.txt"));
	OutputStream primer3ResultsLeftPrimerLocs(njh::files::make_path(setUp.pars_.directoryName_, "primer3_results_leftPrimerLocs.bed"));
	OutputStream primer3ResultsRightPrimerLocs(njh::files::make_path(setUp.pars_.directoryName_, "primer3_results_rightPrimerLocs.bed"));


	primer3ResultsOut << njh::conToStr(VecStr{"seqID",
		"primer_name", "primer_side", "primer_seq", "primer_start", "primer_end", "primer_size", "primer_end_stability", "primer_gc_percent", "primer_hairpin_th", "primer_penalty", "primer_penalty_noSize",  "primer_self_any_th", "primer_self_end_th", "primer_tm", "primer_tm_hairpin_diff", "primer_problems"
		}, "\t") << std::endl;

	uint32_t maxChromLen = 0;
	for(const auto & chrom : chromLens){
		if(chrom.second > maxChromLen){
			maxChromLen = chrom.second;
		}
	}

  for(const auto & res : results){
    const auto & region = nameToRegion[res->sequence_id_];
    for(const auto & primer : res->leftPrimers_){
      double penaltyWithOutSize = primer->penalty_ - uAbsdiff(primer->forwardOrientationPos_.size_, p3Opts.PRIMER_OPT_SIZE);

      MetaDataInName primerMeta;
      primerMeta.addMeta("SeqID", res->sequence_id_);
      primerMeta.addMeta("Name", primer->name_);
      primerMeta.addMeta("RightSide", false);
      primerMeta.addMeta("Penalty", primer->penalty_);
      primerMeta.addMeta("PenaltyWithOutSize", penaltyWithOutSize);

      primer3ResultsLeftPrimerLocs << region->chrom_
                                   << "\t" << region->chromStart_ + primer->forwardOrientationPos_.start_
                                   << "\t" << region->chromStart_ + primer->forwardOrientationPos_.start_ + primer->forwardOrientationPos_.size_
                                   << "\t" << primerMeta.createMetaName()
                                   << "\t" << primer->forwardOrientationPos_.size_
                                   << "\t" << "+" << std::endl;

      primer3ResultsOut << njh::conToStr(toVecStr(
          res->sequence_id_,
          primer->name_,
          "left",
          primer->seq_,
          primer->forwardOrientationPos_.start_,
          primer->forwardOrientationPos_.start_ + primer->forwardOrientationPos_.size_,
          primer->forwardOrientationPos_.size_,
          primer->end_stability_,
          primer->gc_percent_,
          primer->hairpin_th_,
          primer->penalty_,
          penaltyWithOutSize,
          primer->self_any_th_,
          primer->self_end_th_,
          primer->tm_,
          primer->tm_ - primer->hairpin_th_,
          njh::conToStr(primer->problems_, ";")
      ), "\t") << std::endl;
    }

    for(const auto & primer : res->rightPrimers_){
      double penaltyWithOutSize = primer->penalty_ - uAbsdiff(primer->forwardOrientationPos_.size_, p3Opts.PRIMER_OPT_SIZE);

      MetaDataInName primerMeta;
      primerMeta.addMeta("SeqID", res->sequence_id_);
      primerMeta.addMeta("Name", primer->name_);
      primerMeta.addMeta("RightSide", true);
      primerMeta.addMeta("Penalty", primer->penalty_);
      primerMeta.addMeta("PenaltyWithOutSize", penaltyWithOutSize);

      primer3ResultsRightPrimerLocs << region->chrom_
                                    << "\t" << region->chromStart_ + primer->forwardOrientationPos_.start_
                                    << "\t" << region->chromStart_ + primer->forwardOrientationPos_.start_ + primer->forwardOrientationPos_.size_
                                    << "\t" << primerMeta.createMetaName()
                                    << "\t" << primer->forwardOrientationPos_.size_
                                    << "\t" << "+" << std::endl;

      primer3ResultsOut << njh::conToStr(toVecStr(
          res->sequence_id_,
          primer->name_,
          "right",
          primer->seq_,
          primer->forwardOrientationPos_.start_,
          primer->forwardOrientationPos_.start_ + primer->forwardOrientationPos_.size_,
          primer->forwardOrientationPos_.size_,
          primer->end_stability_,
          primer->gc_percent_,
          primer->hairpin_th_,
          primer->penalty_,
          penaltyWithOutSize,
          primer->self_any_th_,
          primer->self_end_th_,
          primer->tm_,
          primer->tm_ - primer->hairpin_th_,
          njh::conToStr(primer->problems_, ";")
      ), "\t") << std::endl;
    }
  }

	return 0;
}






int programWrapperRunner::testPrimersWithPrimer3(const njh::progutils::CmdArgs & inputCommands){
	comparison allowablePairwisePrimerComps;
	allowablePairwisePrimerComps.hqMismatches_ = 1;
	bool getPairwisePrimerComps = false;

	Primer3Runner::Primer3Options p3Opts;

	bfs::path bedFnp;
	bfs::path twoBitFnp;
	bfs::path primersFnp;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	p3Opts.setPrimaryOptions(setUp);
	p3Opts.setPrimerSizeOpts(setUp);

	setUp.setOption(bedFnp, "--bedFnp", "genomic locations", true);
	setUp.setOption(twoBitFnp, "--twoBit", "two Bit file", true);
	setUp.setOption(primersFnp, "--primers", "primers to test", true);
	setUp.setOption(getPairwisePrimerComps, "--getPairwisePrimerComps", "get Pairwise Primer Comps");
	setUp.processComparison(allowablePairwisePrimerComps, "--pairwiseComp");


	setUp.processDirectoryOutputName("testPrimers_TODAY", true);

	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);

	auto regions = getBeds(bedFnp);
	PrimersAndMids primers(primersFnp);
	primers.initPrimerDeterminator();
	TwoBit::TwoBitFile tReader(twoBitFnp);

	std::set<std::string> targetNames = njh::getSetOfMapKeys(primers.targets_);
	std::set<std::string> regionNames;
	for(const auto & reg : regions){
		regionNames.emplace(reg->name_);
	}
	VecStr missingFromTargetNames;
	std::set_difference(
			targetNames.begin(), targetNames.end(),
			regionNames.begin(), regionNames.end(),
			std::back_inserter(missingFromTargetNames));

	VecStr missingFromRegionNames;
	std::set_difference(
			regionNames.begin(), regionNames.end(),
			targetNames.begin(), targetNames.end(),
			std::back_inserter(missingFromRegionNames));
	if(!missingFromRegionNames.empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "missing the following targets from " << bedFnp << "\n";
		ss << njh::conToStr(missingFromRegionNames, ",") << std::endl;
		throw std::runtime_error{ss.str()};
	}
	if(!missingFromTargetNames.empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "missing the following regions from " << primersFnp << "\n";
		ss << njh::conToStr(missingFromTargetNames, ",") << std::endl;

		throw std::runtime_error{ss.str()};
	}



	std::unordered_map<std::string, double> gcContentForRegion;
	{
		OutputStream primer3Input(njh::files::make_path(setUp.pars_.directoryName_, "primer3File.txt"));
		primer3Input << "PRIMER_THERMODYNAMIC_PARAMETERS_PATH=" << p3Opts.primer3ConfigPath.string() << std::endl;

		for(const auto & reg : regions){
			auto seqTemplate = GenomicRegion(*reg).extractSeq(tReader);
			njh::strToUpper(seqTemplate.seq_);
			DNABaseCounter counter;
			counter.increase(seqTemplate.seq_);
			counter.setFractions();
			gcContentForRegion[reg->name_] = roundDecPlaces(counter.calcGcContent() * 100, 2);
			std::string seqID = reg->name_;
			uint32_t seqCount = 0;
			for(const auto & fwd : njh::mapAt(primers.targets_, reg->name_).info_.fwds_){
				for(const auto & rev : njh::mapAt(primers.targets_, reg->name_).info_.revs_){
					if(seqCount > 0){
						seqID = njh::pasteAsStr(seqID, ".", seqCount);
					}
					std::string leftPrimer = fwd.primer_;
					std::string rightPrimer = rev.primer_;
					primer3Input << "SEQUENCE_TEMPLATE=" << seqTemplate.seq_<< std::endl;
					primer3Input << "SEQUENCE_ID=" << seqID<< std::endl;
					primer3Input << "SEQUENCE_PRIMER=" << leftPrimer<< std::endl;
					primer3Input << "SEQUENCE_PRIMER_REVCOMP=" << rightPrimer<< std::endl;
					primer3Input << "PRIMER_PICK_ANYWAY=1" << std::endl;
					primer3Input << "PRIMER_TASK=check_primers" << std::endl;
					primer3Input << "PRIMER_FIRST_BASE_INDEX=0" << std::endl;
					primer3Input << "PRIMER_EXPLAIN_FLAG=1" << std::endl;

					primer3Input << "PRIMER_PRODUCT_SIZE_RANGE=" << len(seqTemplate) << "-" << len(seqTemplate ) + 1 << std::endl;
					primer3Input << "PRIMER_OPT_GC_PERCENT=" << p3Opts.PRIMER_OPT_GC_PERCENT << std::endl;
//					primer3Input << "PRIMER_MAX_SIZE=" << std::min(leftPrimer.size(), rightPrimer.size()) << std::endl;
//					primer3Input << "PRIMER_MIN_SIZE=" << std::max(leftPrimer.size(), rightPrimer.size())  << std::endl;

					primer3Input << "PRIMER_MAX_SIZE=" << p3Opts.PRIMER_MAX_SIZE << std::endl;
					primer3Input << "PRIMER_MIN_SIZE=" << p3Opts.PRIMER_MIN_SIZE << std::endl;
					primer3Input << "PRIMER_OPT_SIZE=" << p3Opts.PRIMER_OPT_SIZE << std::endl;
					primer3Input << "="<< std::endl;
					++seqCount;
				}
			}
		}
	}
	std::string primer3Cmd = njh::pasteAsStr("cd ", setUp.pars_.directoryName_, " && ", "primer3_core ", "primer3File.txt > primer3_output.txt 2> primer3_error.txt");
	auto runOutput = njh::sys::run({primer3Cmd});
	OutputStream primer3RunLog(njh::files::make_path(setUp.pars_.directoryName_, "primer3RunLog.json"));
	primer3RunLog << runOutput.toJson() << std::endl;

	auto primer3ResultsFnp = njh::files::make_path(setUp.pars_.directoryName_, "primer3_output.txt");
	auto results = Primer3Runner::Primer3ResultsGeneric::parsePrimer3OutputResults(primer3ResultsFnp, true);

	OutputStream primer3ResultsOut(njh::files::make_path(setUp.pars_.directoryName_, "primer3_results.tab.txt"));

	primer3ResultsOut << njh::conToStr(VecStr{"seqID", "primerPairName", "compl_any_th", "compl_end_th", "pair_penalty", "pair_penalty_noSize", "product_size", "product_gc_percent",
	"left_name", "left_seq", "left_start", "left_end", "left_size", "left_end_stability", "left_gc_percent", "left_hairpin_th", "left_penalty", "left_penalty_noSize", "left_self_any_th", "left_self_end_th", "left_tm", "left_tm_hairpin_diff", "left_problems",
	"right_name", "right_seq", "right_start", "right_end", "right_size", "right_end_stability", "right_gc_percent", "right_hairpin_th", "right_penalty","right_penalty_noSize", "right_self_any_th", "right_self_end_th", "right_tm", "right_tm_hairpin_diff", "right_problems"
	}, "\t") << std::endl;
	for(const auto & res : results){
		for(const auto & primerPair : res->primerPairs_){

			double leftPenaltyWithOutSize = primerPair->left_.penalty_ - uAbsdiff(primerPair->left_.forwardOrientationPos_.size_, p3Opts.PRIMER_OPT_SIZE);
			double rightPenaltyWithOutSize = primerPair->right_.penalty_ - uAbsdiff(primerPair->right_.forwardOrientationPos_.size_, p3Opts.PRIMER_OPT_SIZE);
			double pairPenaltyWithOutSize = leftPenaltyWithOutSize + rightPenaltyWithOutSize;

			primer3ResultsOut << njh::conToStr(toVecStr(
				res->sequence_id_,
				primerPair->name_,

				primerPair->compl_any_th_,
				primerPair->compl_end_th_,
				primerPair->penalty_,
				pairPenaltyWithOutSize,
				primerPair->product_size_,
				gcContentForRegion[res->sequence_id_],

				primerPair->left_.name_,
				primerPair->left_.seq_,
				primerPair->left_.forwardOrientationPos_.start_,
				primerPair->left_.forwardOrientationPos_.start_ + primerPair->left_.forwardOrientationPos_.size_,
				primerPair->left_.forwardOrientationPos_.size_,
				primerPair->left_.end_stability_,
				primerPair->left_.gc_percent_,
				primerPair->left_.hairpin_th_,
				primerPair->left_.penalty_,
				leftPenaltyWithOutSize,
				primerPair->left_.self_any_th_,
				primerPair->left_.self_end_th_,
				primerPair->left_.tm_,
				primerPair->left_.tm_ - primerPair->left_.hairpin_th_,
				njh::conToStr(primerPair->left_.problems_, ";"),

				primerPair->right_.name_,
				primerPair->right_.seq_,
				primerPair->right_.forwardOrientationPos_.start_,
				primerPair->right_.forwardOrientationPos_.start_ + primerPair->right_.forwardOrientationPos_.size_,
				primerPair->right_.forwardOrientationPos_.size_,
				primerPair->right_.end_stability_,
				primerPair->right_.gc_percent_,
				primerPair->right_.hairpin_th_,
				primerPair->right_.penalty_,
				rightPenaltyWithOutSize,
				primerPair->right_.self_any_th_,
				primerPair->right_.self_end_th_,
				primerPair->right_.tm_,
				primerPair->right_.tm_ - primerPair->right_.hairpin_th_,
				njh::conToStr(primerPair->right_.problems_, ";")
			), "\t") << std::endl;
		}
	}
	if(getPairwisePrimerComps){

	//  std::cout << __FILE__ << " " << __LINE__ << std::endl;
		auto maxPrimerLen = primers.pDeterminator_->getMaxPrimerSize();
	//	std::cout << "maxPrimerLen: " << maxPrimerLen << std::endl;

		aligner alnObject(maxPrimerLen, gapScoringParameters(5,1), substituteMatrix::createDegenScoreMatrixLessN(2, -2));
		//alnObject.countEndGaps_ = true;
		PairwisePairFactory pFactor(primers.pDeterminator_->primers_.size());

		auto primerNames = getVectorOfMapKeys(primers.pDeterminator_->primers_);

		PairwisePairFactory::PairwisePair pair;

		OutputStream primerPairWiseCompsOut(njh::files::make_path(setUp.pars_.directoryName_, "primerPairwiseComps.txt"));

		primerPairWiseCompsOut << "target1\ttarget1Primer\ttarget2\ttarget2Primer\tscore\tidentities\tpass" << std::endl;

		while(pFactor.setNextPair(pair)){
			alnObject.alignCacheGlobalDiag(
					primers.pDeterminator_->primers_.at(primerNames[pair.row_]).forwardPrimerRaw_,
					primers.pDeterminator_->primers_.at(primerNames[pair.col_]).forwardPrimerRaw_);
			alnObject.profilePrimerAlignment(
					primers.pDeterminator_->primers_.at(primerNames[pair.row_]).forwardPrimerRaw_,
					primers.pDeterminator_->primers_.at(primerNames[pair.col_]).forwardPrimerRaw_);
			primerPairWiseCompsOut << primerNames[pair.row_] << "\t" << "forward-5-3" << "\t" << primerNames[pair.col_] << "\t" << "forward-5-3";
			primerPairWiseCompsOut << "\t" << alnObject.comp_.distances_.eventBasedIdentity_ << "\t" <<  alnObject.comp_.distances_.query_.identities_ << "\t" << njh::boolToStr(allowablePairwisePrimerComps.passErrorProfile(alnObject.comp_))<< std::endl;

			alnObject.alignCacheGlobalDiag(
					primers.pDeterminator_->primers_.at(primerNames[pair.row_]).forwardPrimerRaw_,
					seqUtil::reverseComplement(primers.pDeterminator_->primers_.at(primerNames[pair.col_]).forwardPrimerRaw_, "DNA"));
			alnObject.profilePrimerAlignment(
					primers.pDeterminator_->primers_.at(primerNames[pair.row_]).forwardPrimerRaw_,
					seqUtil::reverseComplement(primers.pDeterminator_->primers_.at(primerNames[pair.col_]).forwardPrimerRaw_, "DNA"));
			primerPairWiseCompsOut << primerNames[pair.row_] << "\t" << "forward-5-3" << "\t" << primerNames[pair.col_] << "\t" << "forward-3-5";
			primerPairWiseCompsOut << "\t" << alnObject.comp_.distances_.eventBasedIdentity_ << "\t" <<  alnObject.comp_.distances_.query_.identities_ << "\t" << njh::boolToStr(allowablePairwisePrimerComps.passErrorProfile(alnObject.comp_)) << std::endl;

			alnObject.alignCacheGlobalDiag(
					primers.pDeterminator_->primers_.at(primerNames[pair.row_]).forwardPrimerRaw_,
					primers.pDeterminator_->primers_.at(primerNames[pair.col_]).reversePrimerRaw_);
			alnObject.profilePrimerAlignment(
					primers.pDeterminator_->primers_.at(primerNames[pair.row_]).forwardPrimerRaw_,
					primers.pDeterminator_->primers_.at(primerNames[pair.col_]).reversePrimerRaw_);
			primerPairWiseCompsOut << primerNames[pair.row_] << "\t" << "forward-5-3" << "\t" << primerNames[pair.col_] << "\t" << "reverse-5-3";
			primerPairWiseCompsOut << "\t" << alnObject.comp_.distances_.eventBasedIdentity_ << "\t" <<  alnObject.comp_.distances_.query_.identities_ << "\t" << njh::boolToStr(allowablePairwisePrimerComps.passErrorProfile(alnObject.comp_)) << std::endl;

			alnObject.alignCacheGlobalDiag(
					primers.pDeterminator_->primers_.at(primerNames[pair.row_]).forwardPrimerRaw_,
					seqUtil::reverseComplement(primers.pDeterminator_->primers_.at(primerNames[pair.col_]).reversePrimerRaw_, "DNA"));
			alnObject.profilePrimerAlignment(
					primers.pDeterminator_->primers_.at(primerNames[pair.row_]).forwardPrimerRaw_,
					seqUtil::reverseComplement(primers.pDeterminator_->primers_.at(primerNames[pair.col_]).reversePrimerRaw_, "DNA"));
			primerPairWiseCompsOut << primerNames[pair.row_] << "\t" << "forward-5-3" << "\t" << primerNames[pair.col_] << "\t" << "reverse-3-5";
			primerPairWiseCompsOut << "\t" << alnObject.comp_.distances_.eventBasedIdentity_ << "\t" <<  alnObject.comp_.distances_.query_.identities_ << "\t" << njh::boolToStr(allowablePairwisePrimerComps.passErrorProfile(alnObject.comp_)) << std::endl;

			alnObject.alignCacheGlobalDiag(
					primers.pDeterminator_->primers_.at(primerNames[pair.row_]).reversePrimerRaw_,
					primers.pDeterminator_->primers_.at(primerNames[pair.col_]).forwardPrimerRaw_);
			alnObject.profilePrimerAlignment(
					primers.pDeterminator_->primers_.at(primerNames[pair.row_]).reversePrimerRaw_,
					primers.pDeterminator_->primers_.at(primerNames[pair.col_]).forwardPrimerRaw_);
			primerPairWiseCompsOut << primerNames[pair.row_] << "\t" << "reverse-5-3" << "\t" << primerNames[pair.col_] << "\t" << "forward-5-3";
			primerPairWiseCompsOut << "\t" << alnObject.comp_.distances_.eventBasedIdentity_ << "\t" <<  alnObject.comp_.distances_.query_.identities_ << "\t" << njh::boolToStr(allowablePairwisePrimerComps.passErrorProfile(alnObject.comp_)) << std::endl;

			alnObject.alignCacheGlobalDiag(
					primers.pDeterminator_->primers_.at(primerNames[pair.row_]).reversePrimerRaw_,
					seqUtil::reverseComplement(primers.pDeterminator_->primers_.at(primerNames[pair.col_]).forwardPrimerRaw_, "DNA"));
			alnObject.profilePrimerAlignment(
					primers.pDeterminator_->primers_.at(primerNames[pair.row_]).reversePrimerRaw_,
					seqUtil::reverseComplement(primers.pDeterminator_->primers_.at(primerNames[pair.col_]).forwardPrimerRaw_, "DNA"));
			primerPairWiseCompsOut << primerNames[pair.row_] << "\t" << "reverse-5-3" << "\t" << primerNames[pair.col_] << "\t" << "forward-3-5";
			primerPairWiseCompsOut << "\t" << alnObject.comp_.distances_.eventBasedIdentity_ << "\t" <<  alnObject.comp_.distances_.query_.identities_ << "\t" << njh::boolToStr(allowablePairwisePrimerComps.passErrorProfile(alnObject.comp_)) << std::endl;

			alnObject.alignCacheGlobalDiag(
					primers.pDeterminator_->primers_.at(primerNames[pair.row_]).reversePrimerRaw_,
					primers.pDeterminator_->primers_.at(primerNames[pair.col_]).reversePrimerRaw_);
			alnObject.profilePrimerAlignment(
					primers.pDeterminator_->primers_.at(primerNames[pair.row_]).reversePrimerRaw_,
					primers.pDeterminator_->primers_.at(primerNames[pair.col_]).reversePrimerRaw_);
			primerPairWiseCompsOut << primerNames[pair.row_] << "\t" << "reverse-5-3" << "\t" << primerNames[pair.col_] << "\t" << "reverse-5-3";
			primerPairWiseCompsOut << "\t" << alnObject.comp_.distances_.eventBasedIdentity_ << "\t" <<  alnObject.comp_.distances_.query_.identities_ << "\t" << njh::boolToStr(allowablePairwisePrimerComps.passErrorProfile(alnObject.comp_)) << std::endl;

			alnObject.alignCacheGlobalDiag(
					primers.pDeterminator_->primers_.at(primerNames[pair.row_]).reversePrimerRaw_,
					seqUtil::reverseComplement(primers.pDeterminator_->primers_.at(primerNames[pair.col_]).reversePrimerRaw_, "DNA"));
			alnObject.profilePrimerAlignment(
					primers.pDeterminator_->primers_.at(primerNames[pair.row_]).reversePrimerRaw_,
					seqUtil::reverseComplement(primers.pDeterminator_->primers_.at(primerNames[pair.col_]).reversePrimerRaw_, "DNA"));
			primerPairWiseCompsOut << primerNames[pair.row_] << "\t" << "reverse-5-3" << "\t" << primerNames[pair.col_] << "\t" << "reverse-3-5";
			primerPairWiseCompsOut << "\t" << alnObject.comp_.distances_.eventBasedIdentity_ << "\t" <<  alnObject.comp_.distances_.query_.identities_ << "\t" << njh::boolToStr(allowablePairwisePrimerComps.passErrorProfile(alnObject.comp_)) << std::endl;
		}
	}


	//5-3 and 3-5 first base counts

	DNABaseCounter right53_counter;
	DNABaseCounter left53_counter;
	DNABaseCounter combined53_counter;
	DNABaseCounter right35_counter;
	DNABaseCounter left35_counter;
	DNABaseCounter combined35_counter;

	for(const auto & primer : primers.pDeterminator_->primers_){
		left53_counter.increase(primer.second.forwardPrimerRaw_.front());
		left35_counter.increase(primer.second.forwardPrimerRaw_.back());

		right53_counter.increase(primer.second.reversePrimerRaw_.front());
		right35_counter.increase(primer.second.reversePrimerRaw_.back());

		combined53_counter.increase(primer.second.forwardPrimerRaw_.front());
		combined53_counter.increase(primer.second.reversePrimerRaw_.front());
		combined35_counter.increase(primer.second.forwardPrimerRaw_.back());
		combined35_counter.increase(primer.second.reversePrimerRaw_.back());
	}


	OutputStream firstLastPrimerBaseCounts(njh::files::make_path(setUp.pars_.directoryName_, "firstAndLastPrimerBaseCounts.txt"));
	firstLastPrimerBaseCounts << "label\tbase\tcount\tfreq" << std::endl;
	auto writeCounts = [&firstLastPrimerBaseCounts](DNABaseCounter& counter, const std::string & label){
		counter.setFractions();
		std::vector<char> bases{'A', 'C', 'G','T'};

		for(char base : bases){
			firstLastPrimerBaseCounts << label
					<< "\t" << base
					<< "\t" << counter.bases_[base]
					<< "\t" << counter.baseFracs_[base]
					<< std::endl;
		}
	};
	writeCounts(right53_counter, "5-3_right");
	writeCounts(left53_counter, "5-3_left");
	writeCounts(combined53_counter, "5-3_combined");

	writeCounts(right35_counter, "3-5_right");
	writeCounts(left35_counter, "3-5_left");
	writeCounts(combined35_counter, "3-5_combined");


	return 0;
}


} // namespace njhseq
