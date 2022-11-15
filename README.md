# MaCh3_DUNE

##################################
# CMAKE #########
##################################

Dependencies

- CMake (version > 3.8). 
- MaCh3 Core tag: v1.0-alpha
- ROOT (currently tested on 6.18)

Building:

$ source setup.sh # !! Here you need to make sure that ROOTSYS and Cuda libraries are also set !!
$ source setup_dune_env.sh
$ cd ../
$ mkdir build;
$ cd build

Optional flags are described briefly below, options are shown grouped by square brackets and delimited by vertical lines. Default is on the left.

$ cmake .. -DCPU_ONLY=[OFF|ON] -DSINGLE_THREAD_ONLY=[OFF|ON] -DCUDA_SAMPLES=<path_to_cuda>/CentOS/samples 

CUDA_SAMPLES not necessary if using CPU_ONLY=ON

$ make

If you want to simultaneously develop both the MaCh3 core code and the MaCh3 DUNE code then you can build against a local version of MaCh3 by adding
$ -DCPM_MaCh3_SOURCE=/path/to/MaCh3/folder
this is overrule the CPMFindPackage command in the CMakeList.txt and will tell CPM to build that instead.


Current (November 2022) event rates using DUNE TDR inputs are below. These are made using xsec systematics at their prior central value. Oscillation parameter values used here are:
sin2th12 = 0.307
sin2th23 = 0.52
sin2th13 = 0.0218
dm2_32 = 7.53E-5 eV^2
dm2_12 = 2.509E-3 eV^2 
dCP = -1.601 radians


Integrals of nominal hists: 
FHC_numu unosc:      25941.57467
FHC_numu   osc:      7979.64829
~~~~~~~~~~~~~~~~
Integrals of nominal hists: 
FHC_nue unosc:      391.59946
FHC_nue   osc:      1702.14708
~~~~~~~~~~~~~~~~
Integrals of nominal hists: 
RHC_numu unosc:      12492.61743
RHC_numu   osc:      4219.10087
~~~~~~~~~~~~~~~~
Integrals of nominal hists: 
RHC_nue unosc:      208.80159
RHC_nue   osc:      447.97657
~~~~~~~~~~~~~~~~
