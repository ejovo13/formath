cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

enable_language(Fortran)

project(Mathemathics)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_SOURCE_DIR}/.././modules)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/.././bin)

add_executable(vec_test vector.f90 vec_test.f90 matrix.f90)
add_executable(gs_test vector.f90 gram_schmidt_test.f90 matrix.f90)

target_compile_options(gs_test PRIVATE -Ofast)
