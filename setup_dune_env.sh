source /cvmfs/larsoft.opensciencegrid.org/spack-packages/setup-env.sh
spack load root@6.28.06
spack load cmake
spack load gcc@12.2.0

source /cvmfs/larsoft.opensciencegrid.org/spack-packages/opt/spack/linux-almalinux9-x86_64_v2/gcc-12.2.0/root-6.28.06-jhpj2jsdlwoxbvpnwmxvzkntrxcgw5of/bin/thisroot.sh

spack load hdf5
spack load gsl

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cvmfs/larsoft.opensciencegrid.org/spack-packages/opt/spack/linux-almalinux9-x86_64_v2/gcc-11.4.1/gsl-2.7.1-e6a2ujxmdumzrhuug4b5yk6om4mefz2u/lib
