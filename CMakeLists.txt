set(CMAKE_VERBOSE_MAKEFILE ON)

# CMake version check
cmake_minimum_required(VERSION 3.17 FATAL_ERROR)  

SET(MaCh3_DUNE_VERSION 0.0.1)

# Create project (todo: set configurable CUDA)
project(MaCh3_DUNE VERSION ${MaCh3_DUNE_VERSION} LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)

# Changes default install path to be a subdirectory of the build dir.
# Can set build dir at configure time with -DCMAKE_INSTALL_PREFIX=/install/path
if(CMAKE_INSTALL_PREFIX STREQUAL "" OR CMAKE_INSTALL_PREFIX STREQUAL
  "/usr/local")
  set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}")
elseif(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}")
endif()

# Use the compilers found in the path
find_program(CMAKE_C_COMPILER NAMES $ENV{CC} gcc PATHS ENV PATH NO_DEFAULT_PATH)
find_program(CMAKE_CXX_COMPILER NAMES $ENV{CXX} g++ PATHS ENV PATH NO_DEFAULT_PATH)

LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/Modules)

include(CPM)

#Go get the core MAch3, this will be public one day....
#This is over-riden with -DCPM_MaCh3_SOURCE=/path/to/MaCh3/folder
CPMFindPackage(
    NAME MaCh3
    GIT_TAG "feature/removelibconfig_addyamlcpp"
    GITHUB_REPOSITORY mach3-software/MaCh3
)

#If Single thread is chosen MultiThread flags are omitted
if(NOT ${SINGLE_THREAD_ONLY})

  #Compile and link options (Should be in experiment project but copied in case CORE is built alone) 
  add_compile_options(-I${CUDA_SAMPLES}/common/inc)
  add_compile_options(-g -fopenmp -O3 -Wextra -DMULTITHREAD)
  add_link_options(-fopenmp -DMULTITHREAD)

else()

  #Compile and link options (Should be in experiment project but copied in case CORE is built alone) 
  add_compile_options(-I${CUDA_SAMPLES}/common/inc)
  add_compile_options(-g -O3 -Wextra)

endif()

add_subdirectory(splines)
add_subdirectory(samplePDFDUNE)
add_subdirectory(src)
