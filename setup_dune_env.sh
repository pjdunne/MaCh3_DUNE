source /cvmfs/larsoft.opensciencegrid.org/spack-packages/setup-env.sh
spack load root@6.28.06
spack load cmake
spack load gcc@12.2.0
export CXX=`which g++` # this might be specific for Fermilab?
export CC=`which gcc` # this might be specific for Fermilab?
