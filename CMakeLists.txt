set(CMAKE_VERBOSE_MAKEFILE ON)

# CMake version check
cmake_minimum_required(VERSION 3.17 FATAL_ERROR)  

SET(MaCh3_DUNE_VERSION 0.0.1)

# Create project (todo: set configurable CUDA)
if(${CPU_ONLY})
  project(MaCh3_DUNE VERSION ${MaCh3_DUNE_VERSION} LANGUAGES CXX)
else()
  project(MaCh3_DUNE VERSION ${MaCh3_DUNE_VERSION} LANGUAGES CXX CUDA)
endif()

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
CPMFindPackage(
    NAME CMakeModules
    GIT_TAG v0.0.11
    GITHUB_REPOSITORY NuHepMC/CMakeModules
    DOWNLOAD_ONLY
)
include(${CMakeModules_SOURCE_DIR}/NuHepMCModules.cmake)
include(NuHepMCUtils)

include(ROOT)

#Go get the core MAch3, this will be public one day....
#This is over-riden with -DCPM_MaCh3_SOURCE=/path/to/MaCh3/folder
CPMFindPackage(
    NAME MaCh3
	  GIT_TAG "feature_DUNECore2024-Patches_pyMaCh3_compat"
    GITHUB_REPOSITORY mach3-software/MaCh3
)

# HW: Python nonsense!
set(MACH3_PYVMAJMIN "${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}")
set(MACH3_PYTHONPATH "python/${MACH3_PYVMAJMIN}/pyMaCh3")

CPMFindPackage(
    NAME pybind11
    VERSION 2.9.2
    GITHUB_REPOSITORY pybind/pybind11
    GIT_TAG "v2.9.2"
)


#Compile and link options (Should be in experiment project but copied in case CORE is built alone) 
add_compile_options(-g -O3 -Wextra -fopenmp)
add_link_options(-fopenmp)


#If Single thread is chosen MultiThread flags are omitted
if(NOT ${SINGLE_THREAD_ONLY})
  add_compile_options(-DMULTITHREAD)
  add_link_options(-DMULTITHREAD)
endif()

if(${CPU_ONLY})
  add_compile_options(-DCPU_ONLY)
else()
  add_compile_options(-I${CUDA_SAMPLES}/common/inc)
  add_compile_options(-I$ENV{CUDAPATH}/targets/x86_64-linux/include)
endif()

if(${USE_PROB3})
  add_compile_options(-DUSE_PROB3)
endif()

add_subdirectory(splines)
add_subdirectory(samplePDFDUNE)
add_subdirectory(src)
add_subdirectory(python)

configure_file(cmake/Templates/setup.MaCh3DUNE.sh.in
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/setup.MaCh3DUNE.sh" @ONLY)
install(FILES
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/setup.MaCh3DUNE.sh" DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
