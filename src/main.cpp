
// Created on 2015/01/13
// main.cpp
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


#include "elucidatorPrograms.h"


namespace njhseq {

//void test(){
//	std::map<std::string, double>singleMismatchScores = std::map<std::string, double>{
//					{"AG/TT", 0.71},
//					{"GT/CG", -0.59},
//					{"CA/GC", 0.75},
//					{"TC/AA", 1.33},
//					{"GC/CT", 0.62},
//					{"AG/TA", 0.02},
//					{"TA/AG", 0.42},
//					{"TA/AA", 0.69},
//					{"AG/TG", -0.13},
//					{"CT/GT", -0.12},
//					{"AT/TG", 0.07},
//					{"TG/AT", 0.43},
//					{"CC/GA", 0.79},
//					{"AC/TT", 0.64},
//					{"GT/CC", 0.98},
//					{"CA/GG", 0.03},
//					{"TG/AA", 0.74},
//					{"AC/TC", 1.33},
//					{"CG/GG", -0.11},
//					{"GT/CT", 0.45},
//					{"CG/GT", -0.47},
//					{"TT/AG", 0.34},
//					{"GA/CC", 0.81},
//					{"AT/TC", 0.73},
//					{"TC/AT", 0.97},
//					{"CG/GA", 0.11},
//					{"AA/TA", 0.61},
//					{"CC/GC", 0.7},
//					{"GG/CG", -1.11},
//					{"TT/AT", 0.68},
//					{"CT/GG", -0.32},
//					{"AA/TC", 0.88},
//					{"GC/CA", 0.47},
//					{"CC/GT", 0.62},
//					{"TT/AC", 0.75},
//					{"GA/CG", -0.25},
//					{"CA/GA", 0.43},
//					{"GC/CC", 0.79},
//					{"TG/AG", 0.44},
//					{"GG/CT", 0.08},
//					{"AC/TA", 0.77},
//					{"TA/AC", 0.92},
//					{"CT/GC", 0.4},
//					{"AA/TG", 0.14},
//					{"GG/CA", -0.52},
//					{"GA/CA", 0.17},
//					{"TC/AC", 1.05},
//					{"AT/TT", 0.69}};
//
//	std::cout << "{" << std::endl;
//	for(const auto & p : singleMismatchScores){
//		std::cout << "{" << "\"" << p.first.substr(0,3) << seqUtil::complement(p.first.substr(3,2), "DNA") << "\"" << ", " << p.second << "},"<< std::endl;
//	}
//	std::cout << "}" << std::endl;
//}

}  // namespace njhseq


int main(int argc, char* argv[]) {
//	njhseq::test();
//	return 0;

	try {
		njhseq::elucidatorRunner runner;
		return runner.run(argc, argv);
	} catch (std::exception & e) {
		std::cerr << e.what() << std::endl;
		return 1;
	}
}
