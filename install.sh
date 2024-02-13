#! /usr/bin/env bash

#INSTALL lemon
cd lemon-development
mkdir build
cd build
cmake -DLEMON_ENABLE_GLPK=NO -DLEMON_ENABLE_COIN=NO -DLEMON_ENABLE_ILOG=NO ..
make
#make check

#COPY all files .h and .a to LemonSolver
cd ../..

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