elucidator
=======================
version 1.0.0  

Utilies for manipulating sequencing data.  

To configure developmental branch

```bash
git clone https://github.com/nickjhathaway/elucidator
cd elucidator
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

# Install script
All the aboce commands can be found within ./install.sh and can be ran all at once. This will compile elucidator in elucidator/bin/ which should be added to your path afterwards 

```bash
git clone https://github.com/nickjhathaway/elucidator
cd elucidator
./install.sh

```



# Bash Completion  

elucidator tends to have long flags so that they can be clear what they do but it's somewhat annoying to type them out so bash completion has been added.  Put the content of the file at bashCompletion/SeekDeep into a file ~/.bash_completion and it will be source on your next login or use the bellow command while in the SeekDeep directory  

```bash
./setup.py --addBashCompletion  
```

