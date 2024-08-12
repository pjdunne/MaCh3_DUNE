#!/bin/bash

CURDIR=$PWD

cd build
make
cd src

for VAR in DUNE_atm_chain DUNE_LLHscan_osc DUNE_LLHscan_systs; do
    /cvmfs/larsoft.opensciencegrid.org/spack-packages/opt/spack/linux-almalinux9-x86_64_v2/gcc-11.4.1/gcc-12.2.0-ojrjuib44dbyxvxgvwugxw77gsrkv3yc/bin/g++   -std=c++17 -pipe -fsigned-char -pthread -O3 -DNDEBUG -flto CMakeFiles/test.dir/test.cpp.o -o test   -L/cvmfs/larsoft.opensciencegrid.org/spack-packages/opt/spack/linux-almalinux9-x86_64_v2/gcc-12.2.0/root-6.28.06-jhpj2jsdlwoxbvpnwmxvzkntrxcgw5of/lib/root  -Wl,-rpath,/cvmfs/larsoft.opensciencegrid.org/spack-packages/opt/spack/linux-almalinux9-x86_64_v2/gcc-12.2.0/root-6.28.06-jhpj2jsdlwoxbvpnwmxvzkntrxcgw5of/lib/root:/exp/dune/app/users/barrowd/MaCh3/Atmospherics/build/samplePDFDUNE:/exp/dune/app/users/barrowd/MaCh3/Atmospherics/build/splines:/exp/dune/app/users/barrowd/MaCh3/Atmospherics/build/_deps/mach3-build/mcmc:/exp/dune/app/users/barrowd/MaCh3/Atmospherics/build/_deps/mach3-build/samplePDF:/exp/dune/app/users/barrowd/MaCh3/Atmospherics/build/_deps/mach3-build/OscClass:/exp/dune/app/users/barrowd/MaCh3/Atmospherics/build/_deps/mach3-build/covariance:/exp/dune/app/users/barrowd/MaCh3/Atmospherics/build/_deps/mach3-build/splines:/exp/dune/app/users/barrowd/MaCh3/Atmospherics/build/_deps/mach3-build/manager:/exp/dune/app/users/barrowd/MaCh3/Atmospherics/build/_deps/yaml-cpp-build: -lMinuit ../splines/libsplinesDUNE.so ../_deps/mach3-build/mcmc/libMCMC.so ../_deps/mach3-build/OscClass/libOscClass.so ../_deps/mach3-build/covariance/libCovariance.so ../_deps/mach3-build/splines/libSplines.so ../_deps/mach3-build/manager/libManager.so ../_deps/yaml-cpp-build/libyaml-cpp.so.0.7.0 ../_deps/spdlog-build/libspdlog.a -lCore -lRIO -lXMLIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathMore -lMathCore -lThread -lEG -lGeom -lGenVector -lgomp ../samplePDFDUNE/libSamplePDFDUNE.so ../_deps/mach3-build/samplePDF/libSamplePDF.so
done

cd $CURDIR
