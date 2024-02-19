
# Thesis's project for my bachelor's degree at University of Pisa

A C++ interface between the SMS++ library and the Lemon library for solving
block-structured MinCostFlow problems.

The lemon-development folder contains a copy of the development branch made on 07/29/2020 by Lemon; since there is no git repository for that, in case of newer versions, you'll have to manually update it.

## INSTALL SECTION 

The install.sh bash script allow users to create or insert into the /include and /lib folders.
The include folder contains a copy of the /lemon-development/lemon sub-folder in which header files are stored.
The lib folder contains a copy of the libemon.a file, obtained from the /lemon-development/build folder.

In order to run this script, you need the following dependencies installed:

- CMake   
- Make   
- Clang    
- g++
  
