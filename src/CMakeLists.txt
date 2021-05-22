cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

enable_language(Fortran)
enable_testing(ON)

project(Mathemathics)



set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_SOURCE_DIR}/../modules)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/../lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/../lib)
set(CMAKE_INSTALL_PREFIX ~/lib)

add_compile_options(-Wall -Werror)

add_library(formath vector.f90 matrix.f90)

add_subdirectory(test)
