# MaCh3_DUNE

##################################
# CMAKE #########
##################################

Dependencies

- CMake (version > 3.8). 

Getting the code

To use this repository, you first need to get the core MaCh3 code which contains the main implementations:

$ git clone ./MaCh3_core -b liban_develop https://github.com/mach3-software/MaCh3.git # !! Make sure you use the name the folder MaCh3_core, cmake will look for this directory !!

Then for this repository:

$ git clone -b dummy_branch https://github.com/DUNE/MaCh3_DUNE.git

Building:

$ cd MaCh3_DUNE
$ source setup.sh # !! Here you need to make sure that ROOTSYS and Cuda libraries are also set !!
$ source setup_dune_env.sh
$ cd ../
$ mkdir build build_core; cd build !! # !! Make sure you use the name build_core, cmake will look for this directory !!

Optional flags are described briefly below, options are shown grouped by square brackets and delimited by vertical lines. Default is on the left.

$ cmake ../MaCh3_DUNE -DCPU_ONLY=[OFF|ON] -DSINGLE_THREAD_ONLY=[OFF|ON] -DCUDA_SAMPLES=<path_to_cuda>/CentOS/samples 

CUDA_SAMPLES not necessary if using CPU_ONLY=ON

$ make
