
// Created on 2015/01/13
// main.cpp



#include "elucidatorPrograms.h"




int main(int argc, char* argv[]) {
	try {
		njhseq::elucidatorRunner runner;
		return runner.run(argc, argv);
	} catch (std::exception & e) {
		std::cerr << e.what() << std::endl;
		return 1;
	}
}
