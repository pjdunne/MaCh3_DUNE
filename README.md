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
$ cmake .. -DCPU_ONLY=[OFF|ON] -DUSE_PROB3=[OFF|ON] -DSINGLE_THREAD_ONLY=[OFF|ON] -DCUDA_SAMPLES=<path_to_cuda>/CentOS/samples -DCPM_MaCh3_SOURCE=/path/to/MaCh3_core
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
FHC_numu   osc:      7977.36421
 
FHC_nue unosc:      390.85150
FHC_nue   osc:      1698.28079
 
RHC_numu unosc:      12492.61743
RHC_numu   osc:      4217.78765
 
RHC_nue unosc:      208.31873
RHC_nue   osc:      447.09673
~~~~~~~~~~~~~~~~


####################################
# NDGAr Branch #######
####################################

In order to set up the NDGAr branch and read in ND-GAr CAF files you need to follow these instructions:
You will need to check the paths in scripts/link_files.sh to make sure they point to directories where your files are saved. Likewise, you will need to check the paths in the setup.sh script. If running on SL7 instead of Alma9, then you need to change line 18 on setup.sh to "source setup_dune_env.sh" 

~~~~~~~~~~~~~~
$ mkdir MaCh3_DUNE
$ git clone -b NDGAr https://github.com/DUNE/MaCh3_DUNE.git MaCh3_DUNE
$ cd MaCh3_DUNE
$ source setup.sh
~~~~~~~~~~~~~~
Then need to get a duneanaobj library to read the ND-GAr CAFs. If the CAF format changes then this library will have to be updated

~~~~~~~~~~~~~~
$ cd ..
$ mkdir duneanaobj
$ git clone https://github.com/naseemkhan99/duneanaobj.git duneanaobj
~~~~~~~~~~~~~~

If you need to build MaCh3 core locally

~~~~~~~~~~~~~~
$ mkdir MaCh3_core
$ git clone -b NDGAr https://github.com/mach3-software/MaCh3.git MaCh3_core
~~~~~~~~~~~~~~

Now make a build directory and build MaCh3
~~~~~~~~~~~~~~
$ mkdir build
$ cd build
$ cmake ../MaCh3_DUNE/ -DCPU_ONLY=ON -DUSE_PROB3=OFF -DSINGLE_THREAD_ONLY=OFF -DCPM_duneanaobj_SOURCE=../duneanaobj -DSTANDALONE_BUILD=ON -DCPM_DOWNLOAD_LOCATION=../MaCh3_DUNE/cmake/Modules/CPM.cmake -DCPM_MaCh3_SOURCE=path/to/MaCh3core
$ make
$ cd ../MaCh3_DUNE
~~~~~~~~~~~~~~

Now you can run executables from MaCh3_DUNE directory using config files in configs/ folder.
