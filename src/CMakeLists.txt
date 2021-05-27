cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

enable_language(Fortran)
enable_testing(ON)

project(Formath VERSION 1.0)

set(FORD_CONFIG ${CMAKE_SOURCE_DIR}/../docs)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_SOURCE_DIR}/../modules)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/../lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/../lib)
set(CMAKE_INSTALL_PREFIX /usr/local/Formath-${PROJECT_VERSION}) # This will later get added to the "Import_prefix" defined in the FormathTargets.cmake file

add_compile_options(-Wall) #-Werror)

add_library(Formath vector.f90 matrix.f90)

add_subdirectory(test)

message("Updated CMAKE")

# execute_process(COMMAND ford formath-project-file.md COMMENT "Running Ford" WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/../")
execute_process(COMMAND ford formath-project-file.md WORKING_DIRECTORY ../)

# add_custom_command(OUTPUT ${FORD_CONFIG} COMMAND ${CMAKE_COMMAND} -E ford formath-project-file.md
# WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/../" COMMENT "Running Ford Documentation")

# ===================================== #
# Start configuring the package exports #
# ===================================== #

# make cache variables for install destinations
include(GNUInstallDirs)

# message("Cmake install libdir = ${CMAKE_INSTALL_LIBDIR}")
# message("Cmake install bindir = ${CMAKE_INSTALL_BINDIR}")
# message("Cmake install include = ${CMAKE_INSTALL_INCLUDEDIR}")

# ===================================== #
# Define export locations for Formath   #
# ===================================== #
install(TARGETS Formath 
    EXPORT FormathTargets # Associate the target Formath with the export FormathTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# ===================================== #
# Explicitly install the Targets        #
# ===================================== #
install(EXPORT FormathTargets
        FILE FormathTargets.cmake
        NAMESPACE Formath::
        DESTINATION /usr/local/lib/cmake/Formath
)

install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ 
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
) # Trailing / makes it so just the files are installed

# ===================================== #
# Start creating a Package Config file  #
# ===================================== #

include(CMakePackageConfigHelpers)

configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
        "${CMAKE_CURRENT_BINARY_DIR}/FormathConfig.cmake"
        INSTALL_DESTINATION /usr/local/lib/cmake/Formath
)

install(FILES

        "${CMAKE_CURRENT_BINARY_DIR}/FormathConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/FormathConfigVersion.cmake"
        DESTINATION /usr/local/lib/cmake/Formath

)

# ===================================== #
#    Create a Package Version File      #
# ===================================== #

set(version 1.0)

set_property(TARGET Formath PROPERTY VERSION ${version})
set_property(TARGET Formath PROPERTY SOVERSION 1)
set_property(TARGET Formath PROPERTY INTERFACE_Formath_MAJOR_VERSION 1)
set_property(TARGET Formath APPEND PROPERTY COMPATIBLE_INTERFACE_STRING Formath_MAJOR_VERSION)

write_basic_package_version_file("${CMAKE_CURRENT_BINARY_DIR}/FormathConfigVersion.cmake"
VERSION "${version}"
COMPATIBILITY AnyNewerVersion)