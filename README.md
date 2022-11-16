# MaCh3_DUNE

##################################
# CMAKE #########
##################################

Dependencies

- CMake (version > 3.8). 
- MaCh3 Core tag: v1.0.1-alpha
- ROOT (currently tested on 6.18)

Building:

~~~~~~~~~~~~~~
$ mkdir MaCh3_DUNE
$ git clone git@github.com:DUNE/MaCh3_DUNE.git MaCh3_DUNE
$ cd MaCh3_DUNE
~~~~~~~~~~~~~~

Currently need to first checkout a tag of the core code and set it up manually. When the repo is public CPM will do this for you but for now we'll do it this way.
You can build core whereever, you'll just need to make sure that CPM_MACH3_SOURCE is set to it correctly later.

~~~~~~~~~~~~~~
$ mkdir MaCh3_core
$ git clone git@github.com:mach3-software/MaCh3.git MaCh3_core
$ cd MaCh3_core
$ git checkout tags/v1.0.1-alpha
$ cd ../
~~~~~~~~~~~~~~

Now setup some dependencies and then actually build MaCh3_DUNE

~~~~~~~~~~~~~~~
$ source setup.sh # !! Here you need to make sure that ROOTSYS and Cuda libraries are also set !!
$ source setup_dune_env.sh
$ mkdir build;
$ cd build
~~~~~~~~~~~~~~~

Optional flags are described briefly below, options are shown grouped by square brackets and delimited by vertical lines. Default is on the left.

~~~~~~~~~~~~~~
$ cmake .. -DCPU_ONLY=[OFF|ON] -DSINGLE_THREAD_ONLY=[OFF|ON] -DCUDA_SAMPLES=<path_to_cuda>/CentOS/samples -DCPM_MaCh3_SOURCE=/path/to/MaCh3_core
$ make
~~~~~~~~~~~~~~

A few notes:
CUDA_SAMPLES not necessary if using CPU_ONLY=ON

If you want to simultaneously develop both the MaCh3 core code and the MaCh3 DUNE code then you can build against a local version of MaCh3 by adding:

~~~~~~~~~~~~~~
$ -DCPM_MaCh3_SOURCE=/path/to/MaCh3/folder
~~~~~~~~~~~~~~

this is overrule the CPMFindPackage command in the CMakeList.txt and will tell CPM to build that instead.
As already described this is the default way to build for now. Eventually we won't need this and CPM will find the tag we give it in the CMakeList

###################################
# Event Rates ######
###################################

Once you've got setup you'll then need to setup some symlinks to point to your MC and spline files. You can do this by modifying scripts/link_files.sh script. 

Current (November 2022) event rates using DUNE TDR inputs are below. These are made using xsec systematics at their prior central value. Oscillation parameter values used here are:

sin2th12 = 0.307

sin2th23 = 0.52

sin2th13 = 0.0218

dm2_32 = 7.53E-5 eV^2

dm2_12 = 2.509E-3 eV^2 

dCP = -1.601 radians

~~~~~~~~~~~~~~~~
Integrals of nominal hists: 
FHC_numu unosc:      25941.57467
FHC_numu   osc:      7979.64829

FHC_nue unosc:      391.59946
FHC_nue   osc:      1702.14708

RHC_numu unosc:      12492.61743
RHC_numu   osc:      4219.10087

RHC_nue unosc:      208.80159
RHC_nue   osc:      447.97657
~~~~~~~~~~~~~~~~
