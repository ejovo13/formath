# Formath

Formath is a Fortran library that implements a framework for computing fundamental linear algebra routines like QR decomposition, vector operations (projections, householder transformations), Iterative and Direct eigenvalue methods. 

Formath contains two base derived types, the vector (rank1 array), and matrix (rank1 array of vectors). These based types allow for the rigourous definition of a vector, including some basic operators defined on vectors.

A Matrix object allows certain Linear Algebra algorithms to be called from and object-oriented interface. Implementing these fundamental routines from the ground up is a necessary exercise to intimately know these algorithms.

# Build

Formath has no dependencies and can be built using cmake. Formath was succesfully built on my linux machine when using gfortran.

During the build process, `bin`, `lib`, and `modules` directories will get built. Upon installation (calling `make install`) the Formath library will install the binaries and libraries in standard GNU locations under /usr/local/ as the `CMAKE_INSTALL_PREFIX`. Installation will also install Cmake package config files so that Formath can be used flawlessly in other Cmake projects. An example gist of using Formath in a different project is [here](https://gist.github.com/ejovo13/f2773b441482a6bcbd9471cbd88b0301).

# Documentation

Complete documentation generated with [FORD](https://github.com/Fortran-FOSS-Programmers/ford) is hosted [here](https://ejovo13.github.io/formath/).

# Tests

After having build Formath with cmake, call make test to assure that all the programs compile and run properly.