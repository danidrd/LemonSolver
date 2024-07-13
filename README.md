# LEMONSolver
This project implements `MCFLEMONSolver`, a SMS++ :Solver
for [MCFBlock](https://gitlab.com/smspp/mcfblock) based
on interfacing solvers from the
[LEMON PROJECT](https://lemon.cs.elte.hu/trac/lemon).
In fact, `MCFLEMONSolver` is not a single solver but no less
than 16 classes obtained instantiating 4 different "base"
solvers

Each solver is instantiated on two different types of graphs (SmartDigraph and MCFListDigraph) and on the types of costs and capacities on the edges, which can be Double or Long. Although int is supported, for very large graphs it risks overflowing the size of an int.

The algorithms implemented by the LEMON project are:

- NetworkSimplex
- CostScaling
- CycleCanceling
- CapacityScaling

For each of these algorithms, a class has been created, derived from a base class MCFLemonSolver that unites them, to manage method calls correctly using polymorphism.

The interface with the solvers provided by the LEMON project works thanks to the modifications made to the include/lemon/array_map.h file. Since it has not been updated for some time, it was not compatible with C++20 versions, causing problems, especially with the compilation of Concepts. This file is provided within the include/lemon directory, modified to ensure correctness, and a patch has been proposed to the LEMON developers, who I assume will modify this file to make it available in their repositories.

The files are provided "working," so the user does not need to download LEMON themselves. If they did, the array_map.h file downloaded from LEMON might not have been updated yet, potentially compromising the compilation of MCFLemonSolver.

If LEMON provides a new release and the user wants to update, it is necessary to check the correctness of array_map.h, particularly the correct use of std::allocator.

This repository depends, of course, on MCFBlock, as its proper functioning depends on it. The reasons why MCFLemonSolver strictly depends on MCFBlock are:

- All LEMON structures contained within MCFLemonSolver depend on the relative structures of MCFBlock (the cost map of LEMON depends on the cost vector of MCFBlock, and so on ...)

- For the Modifications, an instance of MCFBlock is modified to reflect its changes to the LEMON structures contained in MCFLemonSolver.

Currently, LEMON provides solvers that are compatible with integer costs and capacities, so tests with double values for these would not work.


## Getting started

These instructions will let you build MCFBlock and MCFSolver on your system.

### Requirements

- The [SMS++ core library](https://gitlab.com/smspp/smspp) and
  its requirements

- [MCFBlock](https://gitlab.com/smspp/mcfblock) and its
  requirements (but no actual `MCFSolver` are needed, since
  `LEMONSolver` provides its own)

Visit [LEMONSolver](https://gitlab.com/smspp/lemonsolver) for instruction for compiling and using this project
