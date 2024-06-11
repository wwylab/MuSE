#!/bin/bash

make clean

mkdir -p lib

# compile boost library
cd boost_1_70_0/
./bootstrap.sh
./b2 --clean
./b2

cp stage/lib/libboost_iostreams.a ../lib/
cd ..

# compile libtcmalloc
if [ -d "gperftools" ]
then
    chmod -R 755 gperftools
    rm -rf gperftools
fi

git clone https://github.com/gperftools/gperftools.git
cd gperftools
git checkout gperftools-2.9.1
./autogen.sh
./configure --libdir="$PWD"
make -j4
make install
cp libtcmalloc_minimal.a ../lib/
cd ..

# combine htslib
if [ -d "htslib" ]
then
    chmod -R 755 htslib
    rm -rf htslib
fi

git clone https://github.com/samtools/htslib.git
cd htslib
git checkout 1.9
#git submodule update --init --recursive
autoreconf -i
./configure
make -j4
cp libhts.a ../lib/

cd ..

make

#clean htslib and gperftools and boost libraries
rm -r boost_1_70_0/stage

chmod -R 755 gperftools
rm -r gperftools

chmod -R 755 htslib
rm -r htslib








