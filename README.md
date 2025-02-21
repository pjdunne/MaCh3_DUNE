# MaCh3_DUNE

##################################
# Building MaCh3 DUNE    #########
##################################

Dependencies:

- gcc (tested on 12.2.0)
- CMake (tested on 3.27.7) 
- ROOT (tested on 6.28.06)

A setup script which pulls cvmfs dependancies is included here:
$ source setup_dune_env.sh

Cloning:

~~~~~~~~~~~~~~
$ mkdir MaCh3_DUNE
$ git clone git@github.com:DUNE/MaCh3_DUNE.git MaCh3_DUNE
$ cd MaCh3_DUNE
$ mkdir build;
$ cd build
~~~~~~~~~~~~~~~

Then perform the cmake build command:

~~~~~~~~~~~~~~
$ cmake .. -DCUDAProb3_ENABLED=ON -DCUDAProb3Linear_ENABLED=ON
$ make install
~~~~~~~~~~~~~~
Additional cmake options are available in the MaCh3-Core README

Then source the installation of MaCh3:
~~~~~~~~~~~~~~
source build/bin/setup.MaCh3DUNE.sh
~~~~~~~~~~~~~~

This sets everything needed, and needs to be re-sourced on each terminal session when using MaCh3 (Along with any dependancies)

###################################
# Event Rates ######
###################################

Once you've got setup you'll then need to setup some symlinks to point to your MC and spline files. You can do this by modifying scripts/link_files.sh script. You'll need to change the FILESDIR variable to point to the relevant folder on your machine. The places these files currently live are listed here:

Imperial College London lx:
~~~~~~~~~~~~~~
/vols/dune/ljw20/
~~~~~~~~~~~~~~

FNAL cluster:
~~~~~~~~~~~~~~
/dune/data/users/lwarsame
~~~~~~~~~~~~~~

ComputeCanada Cedar:
~~~~~~~~~~~~~~
/scratch/liban
~~~~~~~~~~~~~~

NERSC Perlmutter:
~~~~~~~~~~~~~~
/pscratch/sd/l/lwarsame
~~~~~~~~~~~~~~

CVMFS:
~~~~~~~~~~~~~~
/cvmfs/dune.osgstorage.org/pnfs/fnal.gov/usr/dune/persistent/stash/MaCh3/inputs/TDR/v2
~~~~~~~~~~~~~~

Current (Feburary 2024) FD event rates using DUNE FD TDR inputs are below (ND is still under-development). These are made using xsec systematics at their prior central value. Oscillation parameter values used here are:

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

