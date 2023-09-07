source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup cmake v3_19_6
setup ifdhc
# setup root v6_18_02a -q e17:prof
setup root v6_18_02a -q e17:prof
export CXX=`which g++` # this might be specific for Fermilab?
export CC=`which gcc` # this might be specific for Fermilab?
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/build/lib:$PWD/build/_deps/yaml-cpp-build/
