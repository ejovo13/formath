title: Fortran Mathematics (Formath)
src_dir: ./src
project: Formath
output_dir: ./docs
project_dir: ./
media_dir: ./media
project_github: https://github.com/ejovo13/formath
project_website: https://github.com/ejovo13
summary: A Fortran library implementing basic Linear Algebra routines.
author: Evan Voyles
author_description: Learning Fortran by writing Fortran
github: https://github.com/ejovo13 
email: Evan.Voyles@kzoo.edu
author_pic: ./media/thomas_1.png
display: public
exclude: vec_test.f90
    gram_schmidt_test.f90

Formath is a linear algebra framework that is focused on implementing the fundamental numerical linear algebra procedures. Formath includes a vector and matrix types which model rank1 and rank2 arrays, respectively. These types allow for the modularization of vector and matrix operations.

@Note This Library is a work in process and is the process of me learning how to implement fundamental numerical algorithms from the ground up. This package is currently configured to my system, as it is intended principally as a learning experience.

## Vector 

A vector object is just an allocatale rank 1 Fortran array with double precision elements. I wanted to learn modern Fortran's approach to object-oriented programming so I decided wrap up a variety of functionality within the vector class. Vectors can be added to one another, they can be orthogonally projected onto eachother, and they can be transformed about eachother (Householder Transformation)

The vector class is my attempt at creating a programming object that conforms to the mathematical definition.

## Matrix

A matrix object is an allocatable rank 1 Fortran array of **vectors**. We can perform matrix multiplication, access rows and columns, and execute a variety of matrix decompositions. 