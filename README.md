# Formath

Formath is a Fortran library that implements a framework for computing fundamental linear algebra routines like QR decomposition, vector operations (projections, householder transformations), Iterative and Direct eigenvalue methods. 

Formath contains two base derived types, the vector (rank1 array), and matrix (rank1 array of vectors). These based types allow for the rigourous definition of a vector, including some basic operators defined on vectors.

A Matrix object allows certain Linear Algebra algorithms to be called from and object-oriented interface. Implementing these fundamental routines from the ground up is a necessary exercise to intimately know these algorithms.