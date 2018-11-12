/*
 * SlimCounterPos.cpp
 *
 *  Created on: Jun 16, 2016
 *      Author: nick
 */

#include "SlimCounterPos.hpp"

namespace njhseq {

SlimCounterPos::SlimCounterPos(){
	baseCounts_.fill(0);
	/////
	/// baseCounts_ array setup
	/////
	//0-4   high quality forward strand
	//5-9   high quality reverse strand
	//10-14 low  quality forward strand
	//15-19 low  quality reverse strand
}

SlimCounterPos::SlimCounterPos(const std::array<uint32_t, NumOfElements > & data){
	for(auto pos : iter::range(2,22)){
		baseCounts_[pos - 2] = data[pos];
	}

	fwInsertions_ =  data[22];
	revInsertions_ = data[23];
	fwDeletions_  =  data[24];
	revDeletions_  = data[25];
	/////
	/// baseCounts_ array setup
	/////
	//data:2-6    array:0-4   high quality forward strand
	//data:7-11   array:5-9   high quality reverse strand
	//data:12-16  array:10-14 low  quality forward strand
	//data:17-21  array:15-19 low  quality reverse strand
}

SlimCounterPos::SlimCounterPos(const std::vector<uint32_t> & data){
	if(NumOfElements != data.size()){
		std::stringstream ss;
		ss << "Error in : " << __PRETTY_FUNCTION__ << std::endl;
		ss << "Size of input vector is not 26, it's " << data.size() << std::endl;
		throw std::runtime_error{ss.str()};
	}
	for(auto pos : iter::range(2,22)){
		baseCounts_[pos - 2] = data[pos];
	}

	fwInsertions_ =  data[22];
	revInsertions_ = data[23];
	fwDeletions_  =  data[24];
	revDeletions_  = data[25];
	/////
	/// baseCounts_ array setup
	/////
	//data:2-6    array:0-4   high quality forward strand
	//data:7-11   array:5-9   high quality reverse strand
	//data:12-16  array:10-14 low  quality forward strand
	//data:17-21  array:15-19 low  quality reverse strand
}


uint32_t SlimCounterPos::getFwInsertions() const {
	return fwInsertions_;
}

double SlimCounterPos::getFwInsertionsFrac() const {
	return fwInsertions_ / static_cast<double>(getInsertions());
}

uint32_t SlimCounterPos::getRevInsertions() const {
	return revInsertions_;
}

uint32_t SlimCounterPos::getInsertions() const {
	return getFwInsertions() + getRevInsertions();
}

uint32_t SlimCounterPos::getFwDeletions() const{
	return fwDeletions_;
}

double SlimCounterPos::getFwDeletionsFrac() const{
	return fwInsertions_ / static_cast<double>(getDeletions());
}

uint32_t SlimCounterPos::getRevDeletions() const{
	return revDeletions_;
}

uint32_t SlimCounterPos::getDeletions() const{
	return getFwDeletions() + getRevDeletions();
}


void SlimCounterPos::increaseCount(char base, bool reverseStand, bool highQuality,
		uint32_t count) {
	baseCounts_[getBaseCountIndex(reverseStand, highQuality)
			+ SlimCounter::colIndex[base]] += count;
}

void SlimCounterPos::increaseInsertionCount(bool reverseStrand,uint32_t count) {
	if(reverseStrand){
		revInsertions_ += count;
	}else{
		fwInsertions_ += count;
	}
}

void SlimCounterPos::increaseDeletionCount(bool reverseStrand,uint32_t count) {
	if(reverseStrand){
		revDeletions_ += count;
	}else{
		fwDeletions_ += count;
	}
}

std::vector<uint32_t> SlimCounterPos::genOutVec()const{
	std::vector<uint32_t> ret;
	ret.reserve(26);
	ret.push_back(0); // placeholder for refId_
	ret.push_back(0); // placeholder for pos_
	ret.insert(ret.begin() + 2, baseCounts_.begin(), baseCounts_.end());
	ret.push_back(fwInsertions_);
	ret.push_back(revInsertions_);
	ret.push_back(fwDeletions_);
	ret.push_back(revDeletions_);

	return ret;
}

VecStr SlimCounterPos::fullDetailHeader() {
	return toVecStr("refId", "pos","refBase",
			"base",
			"hqCount", "hqFraction","hqTotal",
			"hqFwCount", "hqFwFraction", "hqRevCount",
			"lqTotal", "lqFraction","lqTotal",
			"lqFwCount", "lqFwFraction", "lqRevCount");
}



void SlimCounterPos::fullDetailOut(const std::string & refName, uint32_t pos,
		char refBase, std::ostream & out){
	for(auto base : {'A', 'C', 'G', 'T', 'N'}){
		out << refName
				<< "\t" << pos
				<< "\t" << refBase
				<< "\t" << base

				<< "\t" << getHqBase(base)
				<< "\t" << getHqBaseFrac(base)
				<< "\t" << getHqBaseTotal()

				<< "\t" << getHqFwBase(base)
				<< "\t" << getHqFwBaseFrac(base)
				<< "\t" << getHqRevBase(base)

				<< "\t" << getLqBase(base)
				<< "\t" << getLqBaseFrac(base)
				<< "\t" << getLqBaseTotal()

				<< "\t" << getLqFwBase(base)
				<< "\t" << getLqFwBaseFrac(base)
				<< "\t" << getLqRevBase(base)

				<< std::endl;
	}
}

std::vector<VecStr> SlimCounterPos::fullDetailOutVec(const std::string & refName, uint32_t pos,
		char refBase){
	std::vector<VecStr> ret;
	for(auto base : {'A', 'C', 'G', 'T', 'N'}){
		ret.emplace_back(toVecStr(refName, pos, refBase, base, getHqBase(base),getHqBaseFrac(base),getHqBaseTotal()
				,getHqFwBase(base),getHqFwBaseFrac(base),getHqRevBase(base)
				,getLqBase(base),getLqBaseFrac(base),getLqBaseTotal()
				,getLqFwBase(base),getLqFwBaseFrac(base),getLqRevBase(base)));
	}
	return ret;
}


VecStr SlimCounterPos::hqDetailHeader(){
	return toVecStr("refId", "pos","refBase",
			"base",
			"count", "fraction","total",
			"fwCount", "fwFraction", "revCount",
			"hqFrac");
}

double SlimCounterPos::getInsertionsFrac() const {
	return static_cast<double>(getInsertions()) / getHqEventsTotal();
}

double SlimCounterPos::getDeletionsFrac() const {
	return static_cast<double>(getDeletions()) / getHqEventsTotal();
}


std::vector<VecStr> SlimCounterPos::hqDetailOutVec(const std::string & refName,
		uint32_t pos, char refBase) {
	std::vector<VecStr> ret;
	for (auto base : { 'A', 'C', 'G', 'T', 'N' }) {
		ret.emplace_back(
				toVecStr(refName, pos, refBase, base, getHqBase(base),
						getHqBaseFrac(base), getHqBaseTotal(), getHqFwBase(base),
						getHqFwBaseFrac(base), getHqRevBase(base), getHqFraction()));
	}
	ret.emplace_back(
					toVecStr(refName, pos, refBase, "INS", getInsertions(),
							getInsertionsFrac(), getHqBaseTotal(), getFwInsertions(),
							getFwInsertionsFrac(), getRevInsertions(), getHqFraction()));
	ret.emplace_back(
					toVecStr(refName, pos, refBase, "DEL", getDeletions(),
							getDeletionsFrac(), getHqBaseTotal(), getFwDeletions(),
							getFwDeletionsFrac(), getRevDeletions(), getHqFraction()));
	return ret;
}

uint32_t SlimCounterPos::getBaseTotal() const {
	return getHqBaseTotal() + getLqBaseTotal();
}

uint32_t SlimCounterPos::getEventsTotal() const {
	return getHqBaseTotal() + getLqBaseTotal() + getDeletions() + getInsertions();
}

uint32_t SlimCounterPos::getHqBaseTotal() const {
	return getHqFwBaseTotal() + getHqRevBaseTotal();
}

uint32_t SlimCounterPos::getHqEventsTotal() const {
	return getHqBaseTotal() + getDeletions() + getInsertions();
}

uint32_t SlimCounterPos::getHqFwBaseTotal() const {
	return baseCounts_[0] + baseCounts_[1] + baseCounts_[2] + baseCounts_[3] + baseCounts_[4];
}

uint32_t SlimCounterPos::getHqRevBaseTotal() const {
	return baseCounts_[5] + baseCounts_[6] + baseCounts_[7] + baseCounts_[8] + baseCounts_[9];
}

uint32_t SlimCounterPos::getHqFwBase(char base) const {
	return baseCounts_[SlimCounter::colIndex[base] + getBaseCountIndex(false,true)];
}

uint32_t SlimCounterPos::getHqRevBase(char base) const {
	return baseCounts_[SlimCounter::colIndex[base] + getBaseCountIndex(true,true)];
}

uint32_t SlimCounterPos::getHqBase(char base) const {
	return getHqFwBase(base) + getHqRevBase(base);
}

double SlimCounterPos::getHqFwBaseFrac(char base) const{
	return getHqFwBase(base)/static_cast<double>(getHqBase(base));
}

double SlimCounterPos::getHqBaseFrac(char base) const{
	return getHqBase(base)/static_cast<double>(getHqBaseTotal());
}

uint32_t SlimCounterPos::getLqBaseTotal() const {
	return getLqFwBaseTotal() + getLqRevBaseTotal();
}

uint32_t SlimCounterPos::getLqFwBaseTotal() const {
	return baseCounts_[10] + baseCounts_[11] + baseCounts_[12] + baseCounts_[13];
}

uint32_t SlimCounterPos::getLqRevBaseTotal() const {
	return baseCounts_[15] + baseCounts_[16] + baseCounts_[17] + baseCounts_[18];
}

uint32_t SlimCounterPos::getLqFwBase(char base) const {
	return baseCounts_[SlimCounter::colIndex[base] + getBaseCountIndex(false,false)];
}

uint32_t SlimCounterPos::getLqRevBase(char base) const {
	return baseCounts_[SlimCounter::colIndex[base] + getBaseCountIndex(true,false)];
}

uint32_t SlimCounterPos::getLqBase(char base) const {
	return getLqFwBase(base) + getLqRevBase(base);
}

double SlimCounterPos::getLqFwBaseFrac(char base) const{
	return getLqFwBase(base)/static_cast<double>(getLqBase(base));
}

double SlimCounterPos::getLqBaseFrac(char base) const{
	return getLqBase(base)/static_cast<double>(getLqBaseTotal());
}

double SlimCounterPos::getHqFraction() const{
	return static_cast<double>(getHqBaseTotal())/getBaseTotal();
}


}  // namespace njhseq

