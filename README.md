# elucidator
Experimental code for njhseq libraries

To configure developmental branch

```bash
git checkout develop
./configure.py
```

To configure with specific compilers set the environmental CC and CXX 

e.g. for g++-6

```bash
CC=gcc-6 CXX=g++-6 ./configure.py
```

To download and compile dependencies 

```bash
./setup.py --compfile compfile.mk --outMakefile makefile-common.mk 
```

Then finally to compile 

```bash
make
```

