/*
 * kmerExp_testingSimpleKmerHasher.cpp
 *
 *  Created on: Jul 13, 2021
 *      Author: nick
 */




#include "kmerExp.hpp"
#include <njhseq/objects/kmer.h>
#include "elucidator/objects/dataContainers.h"
#include "elucidator/simulation.h"
#include "elucidator/objects/seqObjects/seqKmers.h"

#include "elucidator/objects/MiscUtility/GenomeSeqSearch.hpp"

#include <njhseq/objects/dataContainers/tables/TableReader.hpp>




namespace njhseq {


namespace StrToNumConverter {

/**@brief Function for converting a string to a number, which is just njh::lexical_cast by default and then several specific int conversions are defined for faster converting
 *
 * @param str the string to convert
 * @return the string convert to a number
 */
	template<typename T>
	T stoToNum(const std::string & str){
		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		return njh::lexical_cast<T>(str);
	}

	template<>
	unsigned short stoToNum(const std::string & str){
		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		return estd::stous(str);
	}

	template<>
	unsigned stoToNum(const std::string & str){
		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		return estd::stou(str);
	}

	template<>
	unsigned long stoToNum(const std::string & str){
		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		return std::stoul(str);
	}

	template<>
	unsigned long long stoToNum(const std::string & str){
		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		return std::stoull(str);
	}

	template<>
	short stoToNum(const std::string & str){
		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		return estd::stos(str);
	}

	template<>
	int stoToNum(const std::string & str){
		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		return std::stoi(str);
	}

	template<>
	long int stoToNum(const std::string & str){
		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		return std::stol(str);
	}

	template<>
	long long int stoToNum(const std::string & str){
		return std::stoll(str);
	}

}  // namespace StrToNumConverter




int testingOfWeirdStod(){

	{
		std::cout << "13313441321414123" << std::endl;
		double testDoub = 13313441321414123;
		std::cout << std::setprecision(20) << testDoub << std::endl;
		std::cout << "13313441321414124" << std::endl;
		double testDoub2 = 13313441321414124;
		std::cout << std::setprecision(20) << testDoub2 << std::endl;
		std::cout << njh::colorBool(testDoub == testDoub2) << std::endl;
		std::cout << __FILE__ << " " << __LINE__ << std::endl;
	}
	{
		std::cout << "13313441321414125" << std::endl;
		double testDoub = 13313441321414125;
		std::cout << std::setprecision(20) << testDoub << std::endl;
		std::cout << "13313441321414126" << std::endl;
		double testDoub2 = 13313441321414126;
		std::cout << std::setprecision(20) << testDoub2 << std::endl;
		std::cout << njh::colorBool(testDoub == testDoub2) << std::endl;
		std::cout << __FILE__ << " " << __LINE__ << std::endl;
	}
	{
		std::cout << "13313441321414127" << std::endl;
		double testDoub = 13313441321414127;
		std::cout << std::setprecision(20) << testDoub << std::endl;
		std::cout << "13313441321414128" << std::endl;
		double testDoub2 = 13313441321414128;
		std::cout << std::setprecision(20) << testDoub2 << std::endl;
		std::cout << njh::colorBool(testDoub == testDoub2) << std::endl;
		std::cout << __FILE__ << " " << __LINE__ << std::endl;
	}


	std::cout << std::endl;

	{
		std::string testStr = "3313441321414121";
		std::cout << "places: " << testStr.size() << std::endl;
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "3313441321414122";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "3313441321414123";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "3313441321414124";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "3313441321414125";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "3313441321414126";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "3313441321414127";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "3313441321414129";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "3313441321414129";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}

	std::cout << __FILE__ << " " << __LINE__ << std::endl;

	{
		std::string testStr = "13313441321414121";
		std::cout << "places: " << testStr.size() << std::endl;
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "13313441321414122";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "13313441321414123";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "13313441321414124";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "13313441321414125";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "13313441321414126";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "13313441321414127";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "13313441321414128";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "13313441321414129";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	std::cout << __FILE__ << " " << __LINE__ << std::endl;

	std::cout << std::endl;

	{
		std::string testStr = "113313441321414121";
		std::cout << "places: " << testStr.size() << std::endl;
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "113313441321414122";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "113313441321414123";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "113313441321414124";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "113313441321414125";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "113313441321414126";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "113313441321414127";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "113313441321414128";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "113313441321414129";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}



	std::cout << __FILE__ << " " << __LINE__ << std::endl;

	{
		std::string testStr = "1213313441321414121";
		std::cout << "places: " << testStr.size() << std::endl;
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "1213313441321414122";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "1213313441321414123";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "1213313441321414124";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "1213313441321414125";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "1213313441321414126";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "1213313441321414127";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "1213313441321414128";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << std::endl;
	{
		std::string testStr = "1213313441321414129";
		std::cout << testStr << std::endl;
		std::cout << "njh::StrToNumConverter::stoToNum<uint64_t>(testStr) :" << njh::StrToNumConverter::stoToNum<uint64_t>(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stod(testStr)            :" << std::setprecision(20)  << stod(testStr) << std::endl;
		std::cout << "std::setprecision(20) std::stold(testStr)           :" << std::setprecision(20)  << std::stold(testStr) << std::endl;
	}
	std::cout << __FILE__ << " " << __LINE__ << std::endl;


	SimpleKmerHash hashifier;

	{
		std::string kmerTest = "AGGAGTTAGCATATACN";
		uint64_t hash = hashifier.hash(kmerTest);
		std::string hashback = hashifier.reverseHash(hash);
		std::cout << "hash    : " << hash << std::endl;
		std::cout << "original: " << kmerTest << std::endl;
		std::cout << "convert : " << hashback << std::endl;

	}

	return 0;

}

int kmerExpRunner::testingSimpleKmerHasher(const njh::progutils::CmdArgs & inputCommands){
	SimpleKmerHash hasher;
	std::string testSeq = "AGGAGTTAGCATATACN";
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(testSeq, "--testSeq", "testSeq");
	setUp.finishSetUp(std::cout);
	{
		auto testSeqHash = hasher.hash(testSeq);
		auto hashBack = hasher.reverseHash(testSeqHash);

		std::cout << "testSeq    : " << testSeq << std::endl;
		std::cout << "testSeqHash: " << testSeqHash << std::endl;
		std::cout << "hashBack   : " << hashBack << std::endl;
		std::cout << "check      : " << njh::colorBool(hashBack == testSeq) << std::endl;
		std::cout << std::endl;
		auto testSeqHashRevComp = hasher.revCompHash(testSeq);
		auto hashBackRevComp = hasher.revCompReverseHash(testSeqHashRevComp);
		auto hashBackNoRevComp = hasher.reverseHash(testSeqHashRevComp);

		std::cout << "testSeq           : " << testSeq << std::endl;
		std::cout << "testSeqHashRevComp: " << testSeqHashRevComp << std::endl;
		std::cout << "hashBackRevComp   : " << hashBackRevComp << std::endl;
		std::cout << "hashBackNoRevComp : " << hashBackNoRevComp << std::endl;
		std::cout << "regRevcomp        : " << seqUtil::reverseComplement(testSeq, "DNA") << std::endl;
		std::cout << "checkHashBack     : " << njh::colorBool(hashBackRevComp == testSeq) << std::endl;
		std::cout << "checkRevComp      : " << njh::colorBool(hashBackNoRevComp == seqUtil::reverseComplement(testSeq, "DNA")) << std::endl;


	}

	return 0;
}


}  //namespace njhseq

