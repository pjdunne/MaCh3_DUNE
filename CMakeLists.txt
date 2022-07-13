set(CMAKE_VERBOSE_MAKEFILE ON)

# CMake version check
cmake_minimum_required(VERSION 3.17 FATAL_ERROR)  


# Create project (todo: set configurable CUDA)
project(DuneMaCh3 LANGUAGES CXX)

# creates a custom target to drive download, update/patch, configure, build, 
# install and test steps of an external project
include(ExternalProject)

# Changes default install path to be a subdirectory of the build dir.
# Can set build dir at configure time with -DCMAKE_INSTALL_PREFIX=/install/path
if(CMAKE_INSTALL_PREFIX STREQUAL "" OR CMAKE_INSTALL_PREFIX STREQUAL
  "/usr/local")
  set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}")
elseif(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}")
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Use the compilers found in the path
find_program(CMAKE_C_COMPILER NAMES $ENV{CC} gcc PATHS ENV PATH NO_DEFAULT_PATH)
find_program(CMAKE_CXX_COMPILER NAMES $ENV{CXX} g++ PATHS ENV PATH NO_DEFAULT_PATH)

# To make output lines about environment vars
include(${CMAKE_SOURCE_DIR}/cmake/cmessage.cmake)

# Check the where to install 
cmessage(STATUS "CMAKE_INSTALL_PREFIX: \"${CMAKE_INSTALL_PREFIX}\"")

# add cmake script files
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
#So that CORE code can check if this has been done 
set(CMAKE_SCRIPT_SETUP) 

#CUDA
if(NOT ${CPU_ONLY})
  include(${CMAKE_SOURCE_DIR}/cmake/CUDASetup.cmake)
endif()

# ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
include(${CMAKE_SOURCE_DIR}/cmake/ROOTSetup.cmake)

# Custom commands
set( CMAKE_EXPORT_COMPILE_COMMANDS ON )

#Compile and link options 
add_compile_options(-I${CUDA_SAMPLES}/common/inc -L${ROOT_LIBDIR})
add_compile_options(-g -fopenmp -std=c++11 -O3 -Wextra)
add_compile_options(-I$ENV{MACH3CORE})
add_link_options(-fopenmp)

# Build components
add_subdirectory(../MaCh3_refactor ../bMaCh3_re)
add_subdirectory(samplePDFDUNE)
add_subdirectory(src)
