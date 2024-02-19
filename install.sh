#! /usr/bin/env bash

#INSTALL lemon
cd lemon-development
mkdir build
cd build
cmake -DLEMON_ENABLE_GLPK=NO -DLEMON_ENABLE_COIN=NO -DLEMON_ENABLE_ILOG=NO ..
make
#make check


cd ../..
#DOWNLOAD boost, netCDF, eigen
wget https://boostorg.jfrog.io/artifactory/main/release/1.84.0/source/boost_1_84_0.tar.gz
tar -xvzf boost_1_84_0.tar.gz
rm boost_1_84_0.tar.gz

git clone https://gitlab.com/libeigen/eigen.git

wget https://downloads.unidata.ucar.edu/netcdf-c/4.9.2/netcdf-c-4.9.2.tar.gz
tar -xvzf netcdf-c-4.9.2.tar.gz
rm netcdf-c-4.9.2.tar.gz

git clone https://github.com/frangio68/Min-Cost-Flow-Class.git
#COPY all files .h and .a to LemonSolver

mkdir -p lib
mkdir -p include
mkdir -p include/lemon
mkdir -p include/lemon/bits
mkdir -p include/lemon/concepts

cp lemon-development/build/lemon/libemon.a lib
cp lemon-development/build/lemon/config.h include/lemon


cp lemon-development/lemon/*.h include/lemon
cp lemon-development/lemon/bits/*.h include/lemon/bits
cp lemon-development/lemon/concepts/*.h include/lemon/concepts
