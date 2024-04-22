#!/usr/bin/env bash

if [ $# -ne 1 ]; then
	echo "Need to supply number of CPUs to use, e.g. ./upgradeDevelop.sh 7"
	exit
fi

#make sure on the master branch and then pull for any new commits
git checkout master
git pull
#re-run install
./configure.py
rm -fr external/build/bamtools/ external/build/jsoncpp/ external/build/restbed/ external/build/hts/
rm -fr external/local/bamtools/ external/local/jsoncpp/ external/local/restbed/ external/local/hts/
rm -fr external/local/pathweaver/ external/local/seekdeep/ external/local/njhseq/ external/local/seqServer/ external/local/TwoBit/ external/local/njhcpp/

./setup.py --compfile compfile.mk --outMakefile makefile-common.mk --overWrite

make clean

make -j $1
