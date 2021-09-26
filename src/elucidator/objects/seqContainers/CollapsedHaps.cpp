/*
 * CollapsedHaps.cpp
 *
 *  Created on: Jun 12, 2021
 *      Author: nick
 */


#include "CollapsedHaps.hpp"

#include <njhseq/concurrency/PairwisePairFactory.hpp>

namespace njhseq {


void CollapsedHaps::setSubNamesToMainSeqPos(){
	subNamesToMainSeqPos_.clear();
	for(const auto  names : iter::enumerate(names_)){
		for(const auto & name : names.second){
			subNamesToMainSeqPos_[name] = names.first;
		}
	}
}

void CollapsedHaps::revCompSeqs(){
	njh::for_each(seqs_, [](std::shared_ptr<seqInfo> & seq){ seq->reverseComplementRead(false, true);});
}

void CollapsedHaps::setFrequencies(uint32_t total){
	for(auto & seq : seqs_){
		seq->frac_ = seq->cnt_/total;
	}
}

void CollapsedHaps::setFrequencies(){
	auto total = getTotalHapCount();
	setFrequencies(total);
}

std::unordered_map<std::string, uint32_t> CollapsedHaps::renameBaseOnFreq(
		const std::string &identifier) {
	//rename based on freq
	std::vector<uint32_t> orderByCnt = getOrderByTopCnt();
	uint32_t seqId = 0;
	for (const auto pos : orderByCnt) {
		MetaDataInName popMeta;
		popMeta.addMeta("HapPopUIDCount",
				static_cast<uint32_t>(std::round(seqs_[pos]->cnt_)));
		seqs_[pos]->name_ = njh::pasteAsStr(identifier, ".",
				leftPadNumStr<uint32_t>(seqId, size()), popMeta.createMetaName());
		++seqId;
	}
	return genSeqNameKey();
}

uint32_t CollapsedHaps::getTotalHapCount() const{
	uint32_t ret = 0;
	for(const auto & seq : seqs_){
		ret += seq->cnt_;
	}
	return ret;
}

uint32_t CollapsedHaps::getTotalUniqueHapCount() const{
	return seqs_.size();
}

size_t CollapsedHaps::size() const{
	return seqs_.size();
}

std::unordered_map<std::string, uint32_t> CollapsedHaps::genSeqNameKey() const {
	std::unordered_map<std::string, uint32_t> ret;
	for (const auto pos : iter::range(seqs_.size())) {
		ret[seqs_[pos]->name_] = pos;
	}
	return ret;
}


VecStr CollapsedHaps::GenPopMeasuresPar::genHeader() const {
	VecStr header { "id", "totalHaplotypes", "uniqueHaplotypes", "singlets",
			"doublets", "expShannonEntropy", "ShannonEntropyE",
			"effectiveNumOfAlleles", "SimpsonIndex", "he", "ExpP3", "ExpP4", "ExpP5",
			"lengthPolymorphism" };
	if (getPairwiseComps) {
		header.emplace_back("avgPercentID");
		header.emplace_back("avgNumOfDiffs");
		header.emplace_back("simpleAvalance");
		header.emplace_back("completeAvalance");
		if (numSegSites_ != std::numeric_limits < uint32_t > ::max()) {
			header.emplace_back("nSegratingSites");
			header.emplace_back("TajimaD");
			header.emplace_back("TajimaDPVal");
		}
	}
	return header;
}



VecStr CollapsedHaps::GenPopMeasuresRes::getOut(const CollapsedHaps & inputSeqs, const std::string & identifier, const GenPopMeasuresPar & pars) const{
	VecStr ret;
	ret = toVecStr(
			identifier,
			inputSeqs.getTotalHapCount(),
			inputSeqs.seqs_.size(),
			divMeasures_.singlets_,
			divMeasures_.doublets_,
			divMeasures_.expShannonEntropy_,
			divMeasures_.ShannonEntropyE_,
			divMeasures_.effectiveNumOfAlleles_,
			divMeasures_.simpsonIndex_,
			divMeasures_.heterozygostiy_,
			divMeasures_.ploidy3_.expectedCOIForPloidy_.at(3),
			divMeasures_.ploidy4_.expectedCOIForPloidy_.at(4),
			divMeasures_.ploidy5_.expectedCOIForPloidy_.at(5),

			njh::boolToStr(inputSeqs.hasLengthVariation(pars.lowVarFreq))
			);
	if(pars.getPairwiseComps){
		njh::addConToVec(ret, toVecStr(avgPMeasures_.avgPercentId, avgPMeasures_.avgNumOfDiffs));
		njh::addConToVec(ret, toVecStr(avgPMeasures_.simpleAvalance_, avgPMeasures_.completeAvalance_));
		if(pars.numSegSites_ != std::numeric_limits<uint32_t>::max()){
			njh::addConToVec(ret, toVecStr(pars.numSegSites_, tajimaRes_.d_, tajimaRes_.pval_beta_));
		}
	}
	return ret;
}

void CollapsedHaps::GenPopMeasuresRes::writeDivMeasures(const OutOptions & outOpts, const CollapsedHaps & inputSeqs, const std::string & identifier, const GenPopMeasuresPar & pars) const {
	OutputStream divMeasuresOut(outOpts);
	divMeasuresOut << njh::conToStr(pars.genHeader(), "\t") << std::endl;
	divMeasuresOut << njh::conToStr(getOut(inputSeqs, identifier, pars), "\t") << std::endl;
}


CollapsedHaps::GenPopMeasuresRes CollapsedHaps::getGeneralMeasuresOfDiversity(const GenPopMeasuresPar &pars,
		const std::shared_ptr<aligner> &alignerObj) {
	GenPopMeasuresRes ret;
	if (pars.getPairwiseComps) {
		if(nullptr == alignerObj){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "getting paired wise comparisons but didn't supply an aligner" << "\n";
			throw std::runtime_error{ss.str()};
		}
		if (size() > 1) {
			if (pars.diagAlnPairwiseComps) {
				ret.allComps_ = getPairwiseCompsDiagAln(*alignerObj, pars.numThreads);
				ret.avgPMeasures_ = getAvgPairwiseMeasures(ret.allComps_);
			} else {
				ret.allComps_ = getPairwiseComps(*alignerObj, pars.numThreads);
				ret.avgPMeasures_ = getAvgPairwiseMeasures(ret.allComps_);
			}
		} else {
			ret.avgPMeasures_.avgNumOfDiffs = 0;
			ret.avgPMeasures_.avgPercentId = 1;
		}
	}
	setFrequencies();
	ret.divMeasures_ = PopGenCalculator::getGeneralMeasuresOfDiversity(seqs_);
	if (pars.getPairwiseComps && size() > 1
			&& std::numeric_limits < uint32_t > ::max() != pars.numSegSites_) {
		if(pars.numSegSites_ == 0){
			ret.tajimaRes_.d_ = 0;
			ret.tajimaRes_.pval_beta_ = 1;
			ret.tajimaRes_.pval_normal_ = 1;
		}else{
			try {
				ret.tajimaRes_ = PopGenCalculator::calcTajimaTest(getTotalHapCount(),
						pars.numSegSites_, ret.avgPMeasures_.avgNumOfDiffs);
			} catch (std::exception &e) {
				//currently doing nothing, some times due to frequency filtering etc the calc throw an exception;
			}
		}
	}
	return ret;
}


std::vector<uint32_t> CollapsedHaps::getReadLenVec() const{
	std::vector<uint32_t> ret;
	for(const auto & seq : seqs_){
		for(uint32_t count = 0; count < seq->cnt_; ++count){
			ret.emplace_back(len(*seq));
		}
	}
	return ret;
}

std::unordered_map<uint32_t, uint32_t> CollapsedHaps::getReadLenMap() const{
	std::unordered_map<uint32_t, uint32_t> ret;
	for(const auto & seq : seqs_){
		ret[len(*seq)] += seq->cnt_;
	}
	return ret;
}


bool CollapsedHaps::hasLengthVariation(const double freqCutOff) const {
	auto readLens = getReadLenMap();
	if (0 == freqCutOff) {
		return readLens.size() > 1;
	}

	bool lengthVariation = false;

	uint32_t mode = 0;
	uint32_t lenMode = std::numeric_limits < uint32_t > ::max();

	for (const auto &lenCount : readLens) {
		if (lenCount.second > mode) {
			mode = lenCount.second;
			lenMode = lenCount.first;
		}
	}
	uint32_t lengthNotMode = 0;
	for (const auto &lenCount : readLens) {
		if (lenCount.first != lenMode) {
			lengthNotMode += lenCount.second;
		}
	}
	if (static_cast<double>(lengthNotMode) / getTotalHapCount() >= freqCutOff) {
		lengthVariation = true;
	}
	return lengthVariation;
}

uint32_t CollapsedHaps::getLongestLenDiff(const double freqCutOff) const {
	if (0 == freqCutOff) {
		std::unordered_map<uint32_t, uint32_t> readLens = getReadLenMap();
		uint32_t longestLen = vectorMaximum(njh::getVecOfMapKeys(readLens));
		uint32_t shortestLen = vectorMinimum(njh::getVecOfMapKeys(readLens));
		uint32_t biggestLenDiff = longestLen - shortestLen;
		return biggestLenDiff;
	}
	std::map<uint32_t, uint32_t> lens;
	double total = 0;
	for (const auto &seq : seqs_) {
		lens[len(*seq)] += seq->cnt_;
		total += seq->cnt_;
	}
	if (lens.size() == 1) {
		return 0;
	}
	uint32_t shortestLen = 0;
	uint32_t longestLen = std::numeric_limits<uint32_t>::max();
	{
		//shortest len
		uint32_t cumTotal = 0;
		for (const auto &len : lens) {
			cumTotal += len.second;
			if (cumTotal / total > freqCutOff) {
				shortestLen = len.first;
				break;
			}
		}
	}
	{
		//longest len
		uint32_t cumTotal = 0;
		for (const auto &len : iter::reversed(lens)) {
			cumTotal += len.second;
			if (cumTotal / total > freqCutOff) {
				longestLen = len.first;
				break;
			}
		}
	}
	if (longestLen > shortestLen) {
		return longestLen - shortestLen;
	} else {
		return 0;
	}
}



std::vector<uint32_t> CollapsedHaps::getOrder(const std::function<bool(const seqInfo &,const seqInfo&)> & comparator) const{
	std::vector<uint32_t> order(seqs_.size());
	njh::iota<uint32_t>(order, 0);
	njh::sort(order,[&comparator,this](uint32_t idx1, uint32_t idx2){
		return comparator(*seqs_[idx1], *seqs_[idx2]);
	});
	return order;
}

std::vector<uint32_t> CollapsedHaps::getOrderByTopCnt() const {
	std::function<bool(const seqInfo &,const seqInfo&)> byCnt = [](const seqInfo & seq1, const seqInfo & seq2){
		return seq1.cnt_ > seq2.cnt_;
	};
	return getOrder(byCnt);
}

std::string CollapsedHaps::getSampleNameFromSeqName(const std::string &name,
		const std::vector<std::string> &possibleSampleMetaFields) {
	std::string sample = name;
	if (MetaDataInName::nameHasMetaData(name)) {
		MetaDataInName seqMeta(name);
		for (const auto &metaField : possibleSampleMetaFields) {
			if (seqMeta.containsMeta(metaField)) {
				sample = seqMeta.getMeta(metaField);
				break;
			}
		}
	}
	return sample;
}

std::set<std::string> CollapsedHaps::getPossibleLabIsolateNames(const std::unordered_set<std::string> & names){
	std::set<std::string> nonFieldSampleNames;
	for(const auto & name : names){
		if(MetaDataInName::nameHasMetaData(name)){
			MetaDataInName meta(name);
			std::string sampleField = "sample";
			if(meta.containsMeta("BiologicalSample")){
				sampleField = "BiologicalSample";
			}
			if(meta.containsMeta(sampleField) && meta.containsMeta("IsFieldSample") && !meta.getMeta<bool>("IsFieldSample")){
				if(meta.containsMeta("site") && "LabIsolate" == meta.getMeta("site")){
					nonFieldSampleNames.emplace(meta.getMeta(sampleField));
				}
			}
		}
	}
	return nonFieldSampleNames;
}

std::set<std::string> CollapsedHaps::getAllSampleNames(){
	std::set<std::string> ret;
	for(const auto & subNames : names_){
		for(const auto & name : subNames){
			ret.emplace(getSampleNameFromSeqName(name, possibleSampleMetaFields_));
		}
	}
	return ret;
}

std::vector<std::unordered_set<std::string>> CollapsedHaps::getSampleNamesPerSeqs(){
	std::vector<std::unordered_set<std::string>> ret;
	for(const auto & subNames : names_){
		std::unordered_set<std::string> samps;
		for(const auto & name : subNames){
			samps.emplace(getSampleNameFromSeqName(name, possibleSampleMetaFields_));
		}
		ret.emplace_back(samps);
	}
	return ret;
}


void CollapsedHaps::writeOutSeqsOrdCnt(const SeqIOOptions &seqOpts) const {
	auto orderByCnt = getOrderByTopCnt();
	SeqOutput uniqueWriter(seqOpts);
	uniqueWriter.openOut();
	for (const auto pos : orderByCnt) {
		uniqueWriter.write(seqs_[pos]);
	}
	uniqueWriter.closeOut();
}

void CollapsedHaps::writeNames(const OutOptions & outOpts) const {
	auto orderByCnt = getOrderByTopCnt();
	OutputStream nameOut(outOpts);
	nameOut << "name\tnumber\tinputNames"	<< std::endl;
	for(const auto pos : orderByCnt){
		nameOut << seqs_[pos]->name_
				<< "\t" << seqs_[pos]->cnt_
				<< "\t" << njh::conToStr(names_[pos], ",") << std::endl;
	}
}


void CollapsedHaps::writeOutMetaFields(const OutOptions &outOpts) const{
	OutputStream out(outOpts);
	std::set<std::string> allMetaFields;
	for(const auto & inputNamesPerSeq : names_){
		for(const auto & name : inputNamesPerSeq){
			if(MetaDataInName::nameHasMetaData(name)){
				MetaDataInName meta(name);
				njh::addVecToSet(njh::getVecOfMapKeys(meta.meta_), allMetaFields);
			}
		}
	}
	out << "CollapsedName";
	out << "\t" << njh::conToStr(allMetaFields, "\t") << std::endl;
	for(const auto pos : iter::range(names_.size())){
		for(const auto & name : names_[pos]){
			out << seqs_[pos]->name_;
			MetaDataInName seqMeta;
			if(MetaDataInName::nameHasMetaData(name)){
				seqMeta = MetaDataInName(name);
			}
			for(const auto & field : allMetaFields){
				std::string val = "NA";
				if(seqMeta.containsMeta(field)){
					val = seqMeta.getMeta(field);
				}
				out << "\t" << val;
			}
			out << std::endl;
		}
	}
}



void CollapsedHaps::writeLabIsolateNames(const OutOptions & outOpts) const {
	auto orderByCnt = getOrderByTopCnt();
	OutputStream metaLabNamesOut(outOpts);
	metaLabNamesOut << "name\tsamples" << std::endl;
	for(const auto pos : orderByCnt){
		std::set<std::string> nonFieldSampleNames = CollapsedHaps::getPossibleLabIsolateNames(names_[pos]);
		metaLabNamesOut << seqs_[pos]->name_ << "\t" << njh::conToStr(nonFieldSampleNames, ",") << std::endl;
	}
}

CollapsedHaps CollapsedHaps::readInReads(const SeqIOOptions & inOpts, const std::unique_ptr<MultipleGroupMetaData> & meta,
		const std::unordered_map<std::string, std::set<std::string>> & metaValuesToAvoid){
	CollapsedHaps ret;
	SeqInput reader(inOpts);
	reader.openIn();
	seqInfo seq;
	uint32_t seqCount = 0;
	std::unordered_set<std::string> allNames;

	while(reader.readNextRead(seq)) {
		if(nullptr != meta) {
			meta->attemptToAddSeqMeta(seq);
		}
		//get meta keys if available
		if(MetaDataInName::nameHasMetaData(getSeqBase(seq).name_)){
			MetaDataInName metaData(getSeqBase(seq).name_);
			bool skip = false;
			for(const auto & ignoreField : metaValuesToAvoid){
				if(metaData.containsMeta(ignoreField.first) && njh::in(metaData.getMeta(ignoreField.first), ignoreField.second)){
					skip = true;
					break;
					//skip this seq
				}
			}
			if(skip){
				continue;
			}
		}
		if(njh::in(seq.name_, allNames)){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "can't have seqs with the same name in input"<< "\n";
			ss << seq.name_ << " found more than once" << "\n";
			throw std::runtime_error{ss.str()};
		}
		allNames.emplace(seq.name_);
		readVec::handelLowerCaseBases(seq, inOpts.lowerCaseBases_);
		if(inOpts.removeGaps_){
			seq.removeGaps();
		}
		seqCount+= std::round(seq.cnt_);
		bool found = false;
		for (const auto pos : iter::range(ret.seqs_.size())) {
			const auto & otherSeq = ret.seqs_[pos];
			if (otherSeq->seq_ == seq.seq_) {
				otherSeq->cnt_ += 1;
				ret.names_[pos].emplace(seq.name_);
				found = true;
				break;
			}
		}
		if (!found) {
			//since we are considering counts here as haplotypes, will make sure cnt_ is set to 1
			seq.cnt_ = 1;
			ret.seqs_.emplace_back(std::make_shared<seqInfo>(seq));
			ret.names_.emplace_back(std::unordered_set<std::string>{seq.name_});
		}
	}
	ret.setFrequencies();
	ret.setSubNamesToMainSeqPos();
	return ret;
}
CollapsedHaps::CompWithAlnSeqs::CompWithAlnSeqs(){

}

CollapsedHaps::CompWithAlnSeqs::CompWithAlnSeqs(const comparison & comp, const std::string & refAln, const std::string & queryAln):comp_(comp), refAlnSeq_(refAln), queryAlnSeq_(queryAln){

}

std::vector<CollapsedHaps::CompWithAlnSeqs> CollapsedHaps::getCompsAgainstRef(const seqInfo & refSeq, aligner & alignerObj, uint32_t numThreads) const {
	std::vector<CompWithAlnSeqs> ret(seqs_.size());
	concurrent::AlignerPool alnPool(alignerObj, numThreads);
	alnPool.initAligners();
	std::vector<uint32_t> positions(seqs_.size());
	njh::iota<uint32_t>(positions, 0);
	njh::concurrent::LockableQueue<uint32_t> posQueue(positions);
	std::mutex mut;
	std::function<void()> alignComp = [this,&posQueue,&refSeq,&alnPool,&ret,&mut,&alignerObj](){
		auto currentAligner = alnPool.popAligner();
		uint32_t pos = std::numeric_limits<uint32_t>::max();
		while(posQueue.getVal(pos)){
			currentAligner->alignCacheGlobal(refSeq, seqs_[pos]);
			currentAligner->profileAlignment(refSeq, seqs_[pos], false, false, false);
			ret[pos] = CompWithAlnSeqs(currentAligner->comp_, currentAligner->alignObjectA_.seqBase_.seq_, currentAligner->alignObjectB_.seqBase_.seq_);
		}
		{
			std::lock_guard<std::mutex> lock(mut);
			alignerObj.alnHolder_.mergeOtherHolder(currentAligner->alnHolder_);
		}
	};
	njh::concurrent::runVoidFunctionThreaded(alignComp, numThreads);
	return ret;
}

std::vector<std::vector<comparison>> CollapsedHaps::getPairwiseComps(aligner & alignerObj, uint32_t numThreads) const{
	PairwisePairFactory pFac(seqs_.size());
	std::vector<std::vector<comparison>> ret;
	for(const auto pos : iter::range(seqs_.size())){
		ret.emplace_back(std::vector<comparison>(pos));
	}
	concurrent::AlignerPool alnPool(alignerObj, numThreads);
	alnPool.initAligners();
	njh::ProgressBar pBar(pFac.totalCompares_);
	std::mutex mut;
	std::function<void()> alignComp = [this,&pFac,&alnPool,&ret,&mut,&alignerObj,&pBar](){
		auto currentAligner = alnPool.popAligner();
		PairwisePairFactory::PairwisePair pair;
		while(pFac.setNextPair(pair)){
			currentAligner->alignCacheGlobal(seqs_[pair.row_], seqs_[pair.col_]);
			currentAligner->profileAlignment(seqs_[pair.row_], seqs_[pair.col_], false, false, false);
			ret[pair.row_][pair.col_] = currentAligner->comp_;
			if(verbose_){
				pBar.outputProgAdd(std::cout, 1, true);
			}
		}
		{
			std::lock_guard<std::mutex> lock(mut);
			alignerObj.alnHolder_.mergeOtherHolder(currentAligner->alnHolder_);
		}
	};
	njh::concurrent::runVoidFunctionThreaded(alignComp, numThreads);
	return ret;
}

std::vector<std::vector<comparison>> CollapsedHaps::getPairwiseCompsDiagAln(aligner & alignerObj, uint32_t numThreads) const{
	PairwisePairFactory pFac(seqs_.size());
	std::vector<std::vector<comparison>> ret;
	for(const auto pos : iter::range(seqs_.size())){
		ret.emplace_back(std::vector<comparison>(pos));
	}

	concurrent::AlignerPool alnPool(alignerObj, numThreads);
	alnPool.initAligners();
	njh::ProgressBar pBar(pFac.totalCompares_);
	std::mutex mut;
	std::function<void()> alignComp = [this,&pFac,&alnPool,&ret,&mut,&alignerObj,&pBar](){
		auto currentAligner = alnPool.popAligner();
		PairwisePairFactory::PairwisePair pair;
		while(pFac.setNextPair(pair)){
			uint32_t blockSize = 100 + uAbsdiff(len(*seqs_[pair.row_]), len(*seqs_[pair.col_]));
			currentAligner->inputAlignmentBlockSize_ = std::max<uint32_t>(100, blockSize);
			currentAligner->inputAlignmentBlockWalkbackSize_ = std::max<uint32_t>(100, blockSize);
			if(blockSize *2 > std::min(len(*seqs_[pair.row_]), len(*seqs_[pair.col_]))){
				currentAligner->alignCacheGlobal(seqs_[pair.row_], seqs_[pair.col_]);
			}else{
				currentAligner->alignCacheGlobalDiag(seqs_[pair.row_], seqs_[pair.col_]);
			}
			currentAligner->profileAlignment(seqs_[pair.row_], seqs_[pair.col_], false, false, false);
			ret[pair.row_][pair.col_] = currentAligner->comp_;
			if(verbose_){
				pBar.outputProgAdd(std::cout, 1, true);
			}
		}
		{
			std::lock_guard<std::mutex> lock(mut);
			alignerObj.alnHolder_.mergeOtherHolder(currentAligner->alnHolder_);
		}
	};
	njh::concurrent::runVoidFunctionThreaded(alignComp, numThreads);
	return ret;
}





CollapsedHaps::AvgPairwiseMeasures CollapsedHaps::getAvgPairwiseMeasures(const std::vector<std::vector<comparison>> & allComps) const{
	AvgPairwiseMeasures ret;
	PairwisePairFactory pFac(seqs_.size());
	PairwisePairFactory::PairwisePair pair;



	while(pFac.setNextPair(pair)){
		uint32_t toalSeqInComp = seqs_[pair.row_]->cnt_ + seqs_[pair.col_]->cnt_;

		uint32_t totalBetweenComps =
				  PairwisePairFactory::getTotalPairwiseComps(toalSeqInComp)
				- PairwisePairFactory::getTotalPairwiseComps(seqs_[pair.row_]->cnt_)
				- PairwisePairFactory::getTotalPairwiseComps(seqs_[pair.col_]->cnt_);
		ret.simpleAvalance_   += (1 - allComps[pair.row_][pair.col_].distances_.eventBasedIdentityHq_) * 2;
		ret.completeAvalance_ += (1 - allComps[pair.row_][pair.col_].distances_.eventBasedIdentityHq_) * seqs_[pair.row_]->frac_ * seqs_[pair.col_]->frac_ * 2;
		ret.avgPercentId  += allComps[pair.row_][pair.col_].distances_.eventBasedIdentityHq_ * totalBetweenComps;
		ret.avgNumOfDiffs += allComps[pair.row_][pair.col_].distances_.getNumOfEvents(true) * totalBetweenComps;
	}

	uint64_t total = getTotalHapCount();
	for(const auto pos : iter::range(getTotalUniqueHapCount() )){
		//add in the 100% percent matches between the same exact seqs for the pairwise comps
		ret.avgPercentId += PairwisePairFactory::getTotalPairwiseComps(seqs_[pos]->cnt_);
	}

	ret.simpleAvalance_ /= (getTotalUniqueHapCount() * (getTotalUniqueHapCount() - 1));
	ret.avgPercentId  /= PairwisePairFactory::getTotalPairwiseComps(total);
	ret.avgNumOfDiffs /= PairwisePairFactory::getTotalPairwiseComps(total);
	return ret;
}


}  // namespace njhseq

